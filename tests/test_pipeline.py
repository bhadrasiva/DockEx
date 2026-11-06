from pathlib import Path
import json
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pipeline
import dock_targets
from workbook import LigandInput, TargetInput


def test_run_docking_job_executes_every_ligand_target_pair(tmp_path, monkeypatch):
    ligands = (
        LigandInput("ligand-a", "C"),
        LigandInput("ligand-b", "CC"),
    )
    targets = (
        TargetInput("target-a", pdb_id="1M17"),
        TargetInput("target-b", pdb_id="1P4O"),
    )
    docking_calls = []

    def fake_process_target(spec, out_dir, radius, padding, ph):
        return {
            "target": spec.name,
            "pdb_id": spec.pdb_id,
            "uniprot_id": None,
            "status": "OK",
            "notes": "",
        }

    def fake_load_prepared_targets(targets_dir):
        return [
            dock_targets.TargetRecord(
                target=target.name,
                pdb_id=target.pdb_id,
                uniprot_id=None,
                chain="A",
                receptor_pdbqt=targets_dir / f"{target.name}.pdbqt",
                box_center=(0.0, 0.0, 0.0),
                box_size=(20.0, 20.0, 20.0),
                box_basis="ligand_centroid",
                base_id=target.name,
            )
            for target in targets
        ]

    def fake_prepare_ligand(ligand, ligand_dir):
        path = ligand_dir / f"{ligand.safe_name}.pdbqt"
        path.write_text("ligand")
        return path

    def fake_dock(*args):
        docking_calls.append(args[5])
        return dock_targets.DockingResult(-7.0, "pose.pdbqt")

    monkeypatch.setattr(pipeline.preparation, "process_target", fake_process_target)
    monkeypatch.setattr(pipeline.docking, "load_prepared_targets", fake_load_prepared_targets)
    monkeypatch.setattr(pipeline.docking, "prepare_ligand_pdbqt", fake_prepare_ligand)
    monkeypatch.setattr(pipeline.docking, "dock", fake_dock)
    monkeypatch.setattr(pipeline, "analyze_pose_interactions", lambda receptor, pose: [])

    progress_updates = []
    results_path, summary_path = pipeline.run_docking_job(
        ligands,
        targets,
        tmp_path / "run",
        Path(sys.executable),
        progress=lambda completed, total, message: progress_updates.append((completed, total)),
    )

    assert len(docking_calls) == len(ligands) * len(targets)
    assert results_path.exists()
    assert summary_path.exists()
    assert (tmp_path / "run" / "docking_scores.csv").exists()
    assert (tmp_path / "run" / "interactions.csv").exists()
    assert (tmp_path / "run" / "application.log").exists()
    assert progress_updates[-1][0] == progress_updates[-1][1]
    scores = pd.read_excel(results_path)
    assert scores.shape == (2, 7)
    assert scores["ligand-a"].tolist() == [-7.0, -7.0]


def test_run_docking_job_isolates_pair_failures_and_records_every_pair(
    tmp_path, monkeypatch
):
    ligands = (
        LigandInput("ligand-a", "C"),
        LigandInput("ligand-b", "CC"),
    )
    targets = (
        TargetInput("target-a", pdb_id="1M17"),
        TargetInput("target-b", pdb_id="1P4O"),
    )
    output_dir = tmp_path / "run"
    vina = Path(sys.executable)

    def fake_process_target(spec, out_dir, radius, padding, ph):
        return {
            "target": spec.name,
            "pdb_id": spec.pdb_id,
            "uniprot_id": None,
            "status": "OK" if spec.name == "target-a" else "FAILED",
            "notes": "structure unavailable" if spec.name == "target-b" else "",
        }

    def fake_load_prepared_targets(targets_dir):
        return [
            dock_targets.TargetRecord(
                target="target-a",
                pdb_id="1M17",
                uniprot_id=None,
                chain="A",
                receptor_pdbqt=targets_dir / "target-a.pdbqt",
                box_center=(0.0, 0.0, 0.0),
                box_size=(20.0, 20.0, 20.0),
                box_basis="ligand_centroid",
                base_id="target-a",
            )
        ]

    def fake_prepare_ligand(ligand, ligand_dir):
        path = ligand_dir / f"{ligand.safe_name}.pdbqt"
        path.write_text("ligand")
        return path

    def fake_dock(*args):
        if args[5].endswith("ligand-b"):
            raise RuntimeError("simulated Vina failure")
        return dock_targets.DockingResult(-7.0, "pose.pdbqt")

    monkeypatch.setattr(pipeline.preparation, "process_target", fake_process_target)
    monkeypatch.setattr(
        pipeline.docking, "load_prepared_targets", fake_load_prepared_targets
    )
    monkeypatch.setattr(pipeline.docking, "prepare_ligand_pdbqt", fake_prepare_ligand)
    monkeypatch.setattr(pipeline.docking, "dock", fake_dock)
    def fail_interaction_analysis(receptor, pose):
        raise ValueError("simulated ProLIF failure")

    monkeypatch.setattr(
        pipeline, "analyze_pose_interactions", fail_interaction_analysis
    )

    workbook_bytes = b"workbook-content"
    results_path, _ = pipeline.run_docking_job(
        ligands,
        targets,
        output_dir,
        vina,
        input_workbook=workbook_bytes,
        input_filename="input.xlsx",
    )

    pairs = pd.read_csv(output_dir / "docking_pairs.csv").set_index(
        ["target", "ligand"]
    )
    assert len(pairs) == 4
    assert pairs.loc[("target-a", "ligand-a"), "status"] == "SUCCEEDED"
    assert pairs.loc[("target-a", "ligand-b"), "status"] == "DOCKING_FAILED"
    assert pairs.loc[("target-b", "ligand-a"), "status"] == "TARGET_PREPARATION_FAILED"
    assert pairs.loc[("target-b", "ligand-b"), "status"] == "TARGET_PREPARATION_FAILED"

    scores = pd.read_excel(results_path).set_index("target")
    assert pd.isna(scores.loc["target-b", "ligand-a"])
    metadata = json.loads((output_dir / "run_metadata.json").read_text())
    assert metadata["pairs"] == {
        "requested": 4,
        "succeeded": 1,
        "failed_or_skipped": 3,
        "interaction_analysis_failures": 1,
    }
    interaction_errors = pd.read_csv(output_dir / "interaction_analysis_errors.csv")
    assert interaction_errors.loc[0, "error"] == "simulated ProLIF failure"
    assert metadata["source_workbook"]["filename"] == "input.xlsx"
    assert metadata["source_workbook"]["sha256"]
    assert (output_dir / "input" / "workbook.xlsx").read_bytes() == workbook_bytes
