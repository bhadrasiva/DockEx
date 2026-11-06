"""Run target preparation and ligand-by-target docking."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import logging
from pathlib import Path
import subprocess
from typing import Callable

import pandas as pd

from interactions import analyze_pose_interactions
from engine import dock_targets as docking
from engine import prepare_targets as preparation
from version import APP_VERSION
from workbook import LigandInput, TargetInput

logger = logging.getLogger(__name__)
ProgressCallback = Callable[[int, int, str], None]


def _notify(
    callback: ProgressCallback | None, completed: int, total: int, message: str
) -> None:
    if callback is not None:
        callback(completed, total, message)


def _vina_version(executable: Path) -> str:
    try:
        result = subprocess.run(
            [str(executable), "--version"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"unavailable: {exc}"
    output = (result.stdout or result.stderr).strip()
    return output or f"version command exited with code {result.returncode}"


def _package_versions() -> dict[str, str]:
    packages = {
        "biopython": "biopython",
        "google-genai": "google-genai",
        "meeko": "meeko",
        "numpy": "numpy",
        "openbabel": "openbabel",
        "pandas": "pandas",
        "plotly": "plotly",
        "prolif": "prolif",
        "rdkit": "rdkit",
        "streamlit": "streamlit",
    }
    versions = {}
    for label, distribution in packages.items():
        try:
            versions[label] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[label] = "not-installed"
    return versions


def _relative_path(path: str | Path | None, root: Path) -> str | None:
    if path is None:
        return None
    value = Path(path)
    if not value.is_absolute():
        return value.as_posix()
    try:
        return value.resolve().relative_to(root.resolve()).as_posix()
    except (OSError, ValueError):
        return None


def run_docking_job(
    ligands: tuple[LigandInput, ...],
    targets: tuple[TargetInput, ...],
    output_dir: Path,
    vina_executable: Path,
    exhaustiveness: int = 8,
    seed: int | None = None,
    progress: ProgressCallback | None = None,
    *,
    input_workbook: bytes | None = None,
    input_filename: str | None = None,
) -> tuple[Path, Path]:
    """Run all requested pairs and report individual failures without aborting the batch."""
    if not vina_executable.is_file():
        raise FileNotFoundError(f"Vina executable not found: {vina_executable}")
    if exhaustiveness < 1:
        raise ValueError("Exhaustiveness must be a positive integer.")
    if not ligands or not targets:
        raise ValueError("At least one ligand and one target are required.")

    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = output_dir / "input"
    prepared_dir = output_dir / "targets_prepared"
    ligand_dir = output_dir / "ligands"
    config_dir = output_dir / "config"
    docking_log_dir = output_dir / "logs"
    pose_dir = output_dir / "poses"
    for directory in (
        input_dir,
        ligand_dir,
        config_dir,
        docking_log_dir,
        pose_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)
    if input_workbook is not None:
        (input_dir / "workbook.xlsx").write_bytes(input_workbook)

    log_path = output_dir / "application.log"
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
    )
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)

    pair_count = len(targets) * len(ligands)
    total = len(targets) + len(ligands) + 2 * pair_count
    completed = 0
    try:
        target_rows: list[dict] = []
        for target in targets:
            spec = preparation.TargetSpec(
                name=target.name,
                pdb_id=target.pdb_id,
                uniprot_id=target.uniprot_id,
                ligand_override=target.ligand_resname,
                manual_pocket_residues=list(target.pocket_residues),
                forced_chain=target.chain,
            )
            try:
                row = preparation.process_target(
                    spec,
                    prepared_dir,
                    preparation.DEFAULT_POCKET_RADIUS,
                    preparation.DEFAULT_BOX_PADDING,
                    preparation.DEFAULT_PH,
                )
            except Exception as exc:
                logger.exception("Target preparation failed for %s", target.name)
                row = {
                    "target": target.name,
                    "pdb_id": target.pdb_id,
                    "uniprot_id": target.uniprot_id,
                    "status": "FAILED",
                    "notes": str(exc),
                }
            target_rows.append(row)
            completed += 1
            _notify(progress, completed, total, f"Prepared target {target.name}")

        prepared_by_name = {}
        if any(
            row.get("status") in {"OK", "BLIND_WHOLE_PROTEIN"}
            for row in target_rows
        ):
            try:
                prepared_by_name = {
                    target.target: target
                    for target in docking.load_prepared_targets(prepared_dir)
                }
            except (Exception, SystemExit):
                logger.exception("Could not load prepared target records")

        target_rows_by_name = {row["target"]: row for row in target_rows}
        for target in targets:
            if target.name not in prepared_by_name:
                row = target_rows_by_name[target.name]
                if row.get("status") in {"OK", "BLIND_WHOLE_PROTEIN"}:
                    row["status"] = "FAILED"
                    row["notes"] = (
                        str(row.get("notes", "")).strip()
                        or "Prepared receptor or pocket record was not found."
                    )

        ligand_pdbqts: dict[str, Path | None] = {}
        ligand_errors: dict[str, str] = {}
        for ligand in ligands:
            spec = docking.LigandSpec(name=ligand.name, smiles=ligand.smiles)
            try:
                ligand_pdbqts[ligand.name] = docking.prepare_ligand_pdbqt(
                    spec, ligand_dir
                )
                if ligand_pdbqts[ligand.name] is None:
                    ligand_errors[ligand.name] = (
                        "Ligand preparation failed; see application.log."
                    )
            except Exception as exc:
                logger.exception("Ligand preparation failed for %s", ligand.name)
                ligand_pdbqts[ligand.name] = None
                ligand_errors[ligand.name] = str(exc)
            completed += 1
            _notify(progress, completed, total, f"Prepared ligand {ligand.name}")

        score_rows: list[dict] = []
        pair_rows: list[dict] = []
        interaction_jobs: list[tuple[str, str, Path]] = []
        for requested_target in targets:
            target_record = prepared_by_name.get(requested_target.name)
            target_info = target_rows_by_name[requested_target.name]
            row: dict[str, object] = {
                "target": requested_target.name,
                "pdb_id": target_info.get("pdb_id", requested_target.pdb_id),
                "uniprot_id": target_info.get(
                    "uniprot_id", requested_target.uniprot_id
                ),
                "chain": target_info.get("chain"),
                "box_basis": target_info.get("box_basis"),
            }

            for ligand in ligands:
                score = None
                pose_path = None
                error = ""
                if target_record is None:
                    pair_status = "TARGET_PREPARATION_FAILED"
                    error = str(target_info.get("notes", "Target is unavailable."))
                elif ligand_pdbqts[ligand.name] is None:
                    pair_status = "LIGAND_PREPARATION_FAILED"
                    error = ligand_errors.get(ligand.name, "Ligand is unavailable.")
                else:
                    run_id = (
                        f"{target_record.base_id}_"
                        f"{docking.LigandSpec(ligand.name, ligand.smiles).safe_name}"
                    )
                    try:
                        result = docking.dock(
                            vina_executable,
                            target_record.receptor_pdbqt,
                            ligand_pdbqts[ligand.name],
                            target_record.box_center,
                            target_record.box_size,
                            run_id,
                            exhaustiveness,
                            seed,
                            config_dir,
                            docking_log_dir,
                            pose_dir,
                        )
                        score = result.docking_score
                        pose_path = result.pose_path
                        if score is None:
                            pair_status = "DOCKING_FAILED"
                            error = getattr(
                                result, "error", None
                            ) or "Vina did not return a docking score; see its log."
                        else:
                            pair_status = "SUCCEEDED"
                            if pose_path:
                                interaction_jobs.append(
                                    (requested_target.name, ligand.name, Path(pose_path))
                                )
                    except Exception as exc:
                        logger.exception(
                            "Docking failed for %s against %s",
                            ligand.name,
                            requested_target.name,
                        )
                        pair_status = "DOCKING_FAILED"
                        error = str(exc)

                row[ligand.name] = score
                pair_rows.append(
                    {
                        "target": requested_target.name,
                        "ligand": ligand.name,
                        "status": pair_status,
                        "docking_score": score,
                        "pose_path": _relative_path(pose_path, output_dir),
                        "error": error,
                    }
                )
                completed += 1
                _notify(
                    progress,
                    completed,
                    total,
                    f"Docked {ligand.name} against {requested_target.name}"
                    if pair_status == "SUCCEEDED"
                    else f"{pair_status.replace('_', ' ').title()}: "
                    f"{ligand.name} vs {requested_target.name}",
                )
            score_rows.append(row)

        interaction_jobs_by_pair = {
            (target_name, ligand_name): pose
            for target_name, ligand_name, pose in interaction_jobs
        }
        interaction_rows: list[dict[str, str]] = []
        interaction_errors: list[dict[str, str]] = []
        for target in targets:
            for ligand in ligands:
                pose_path = interaction_jobs_by_pair.get((target.name, ligand.name))
                if pose_path is None:
                    completed += 1
                    _notify(
                        progress,
                        completed,
                        total,
                        f"Interaction analysis skipped for {ligand.name} vs {target.name}",
                    )
                    continue
                receptor_path = prepared_by_name[target.name].receptor_pdbqt.with_name(
                    f"{prepared_by_name[target.name].receptor_pdbqt.stem}_clean.pdb"
                )
                try:
                    rows = analyze_pose_interactions(receptor_path, pose_path)
                    interaction_rows.extend(
                        {
                            "target": target.name,
                            "ligand": ligand.name,
                            "residue": str(item["residue"]),
                            "interaction": str(item["interaction"]),
                        }
                        for item in rows
                    )
                except Exception as exc:
                    logger.exception(
                        "Interaction analysis failed for %s vs %s",
                        ligand.name,
                        target.name,
                    )
                    interaction_errors.append(
                        {"target": target.name, "ligand": ligand.name, "error": str(exc)}
                    )
                completed += 1
                _notify(
                    progress,
                    completed,
                    total,
                    f"Analyzed interactions for {ligand.name} vs {target.name}",
                )

        summary_path = output_dir / "targets_summary.xlsx"
        results_path = output_dir / "docking_scores.xlsx"
        results_frame = pd.DataFrame(score_rows)
        results_frame.to_excel(results_path, index=False)
        results_frame.to_csv(output_dir / "docking_scores.csv", index=False)
        pd.DataFrame(pair_rows).to_csv(output_dir / "docking_pairs.csv", index=False)
        pd.DataFrame(
            interaction_rows,
            columns=["target", "ligand", "residue", "interaction"],
        ).to_csv(output_dir / "interactions.csv", index=False)
        pd.DataFrame(
            interaction_errors, columns=["target", "ligand", "error"]
        ).to_csv(output_dir / "interaction_analysis_errors.csv", index=False)
        pd.DataFrame(target_rows).to_excel(summary_path, index=False)

        succeeded = sum(row["status"] == "SUCCEEDED" for row in pair_rows)
        target_failures = sum(
            row.get("status") == "FAILED" for row in target_rows
        )
        pair_failures = len(pair_rows) - succeeded
        interaction_failure_count = len(interaction_errors)
        run_status = (
            "COMPLETED"
            if not pair_failures and not target_failures and not interaction_failure_count
            else "COMPLETED_WITH_ERRORS"
        )
        source_hash = (
            hashlib.sha256(input_workbook).hexdigest()
            if input_workbook is not None
            else None
        )
        metadata = {
            "schema_version": 1,
            "application": "DockEx",
            "application_version": APP_VERSION,
            "run_id": output_dir.name,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "status": run_status,
            "source_workbook": {
                "filename": input_filename,
                "stored_path": "input/workbook.xlsx"
                if input_workbook is not None
                else None,
                "sha256": source_hash,
            },
            "inputs": {
                "ligands": [
                    {"name": ligand.name, "smiles": ligand.smiles}
                    for ligand in ligands
                ],
                "targets": [
                    {
                        "name": target.name,
                        "pdb_id": target.pdb_id,
                        "uniprot_id": target.uniprot_id,
                        "ligand_resname": target.ligand_resname,
                        "pocket_residues": list(target.pocket_residues),
                        "chain": target.chain,
                    }
                    for target in targets
                ],
            },
            "settings": {
                "vina_executable": vina_executable.name,
                "vina_version": _vina_version(vina_executable),
                "exhaustiveness": exhaustiveness,
                "seed": seed,
                "pocket_radius_angstrom": preparation.DEFAULT_POCKET_RADIUS,
                "box_padding_angstrom": preparation.DEFAULT_BOX_PADDING,
                "protonation_ph": preparation.DEFAULT_PH,
            },
            "software": _package_versions(),
            "targets": [
                {
                    field: row[field]
                    for field in (
                        "target",
                        "pdb_id",
                        "uniprot_id",
                        "source",
                        "status",
                        "chain",
                        "ligand_resname",
                        "pocket_residue_count",
                        "box_basis",
                        "box_center",
                        "box_size",
                        "notes",
                    )
                    if field in row
                }
                for row in target_rows
            ],
            "pairs": {
                "requested": pair_count,
                "succeeded": succeeded,
                "failed_or_skipped": pair_failures,
                "interaction_analysis_failures": interaction_failure_count,
            },
            "artifacts": {
                "input_workbook": "input/workbook.xlsx"
                if input_workbook is not None
                else None,
                "docking_scores": "docking_scores.csv",
                "docking_pairs": "docking_pairs.csv",
                "interactions": "interactions.csv",
                "interaction_errors": "interaction_analysis_errors.csv",
                "target_summary": "targets_summary.xlsx",
            },
        }
        (output_dir / "run_metadata.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False),
            encoding="utf-8",
        )
        (output_dir / "run_summary.txt").write_text(
            f"Status: {run_status}\n"
            f"Input ligands: {len(ligands)}\n"
            f"Input targets: {len(targets)}\n"
            f"Requested ligand-target pairs: {pair_count}\n"
            f"Successful docking pairs: {succeeded}\n"
            f"Failed or skipped pairs: {pair_failures}\n"
            f"Interaction analysis failures: {interaction_failure_count}\n",
            encoding="utf-8",
        )
        return results_path, summary_path
    finally:
        root_logger.removeHandler(file_handler)
        file_handler.close()
