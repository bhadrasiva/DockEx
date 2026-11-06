"""
dock_targets.py

Docking driver for multi-target molecular docking. Second stage of
the pipeline -- run this after prepare_targets.py.

Takes:
    - the directory produced by prepare_targets.py
      (receptors/*.pdbqt + pockets/*_pocket.json -- one pocket record
      per receptor, each holding its own docking box)
    - an Excel of ligands (name, smiles)

and docks every ligand against every prepared target, using each
target's own box (ligand-centroid, manual-residue, or whole-protein
blind -- whatever prepare_targets.py derived for it), via
AutoDock Vina.

Ligand prep (RDKit embed -> Meeko PDBQT) and the Vina subprocess call
follow the same pattern as the original single-target docking script:
one structured output per run (config/log/pose), so every docking run
is fully reproducible from disk alone. Targets that prepare_targets.py marked FAILED have no pocket JSON and are skipped
automatically; everything else (OK or BLIND_WHOLE_PROTEIN) is docked.

--------------------------------------------------------------------
INPUT LIGAND EXCEL SCHEMA
--------------------------------------------------------------------
    name / ligand / compound      ligand label -- used as the results
                                    column header and in output filenames
    smiles / smi                  ligand SMILES string

--------------------------------------------------------------------
OUTPUT LAYOUT (under --outdir, default "docking_results")
--------------------------------------------------------------------
    ligands/<name>.pdbqt                  prepared ligand (once per ligand)
    config/<target>_<name>.txt            Vina parameters used for that run
    logs/<target>_<name>.log              full Vina console output
    poses/<target>_<name>_docked.pdbqt    docked poses
    results.xlsx                          rows = targets, columns = ligand scores

--------------------------------------------------------------------
USAGE
--------------------------------------------------------------------
    python dock_targets.py data\\drugs.xlsx --targets-dir targets_prepared ^
        --vina ..\\vina.exe --outdir docking_results

--------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------
    pip install rdkit meeko pandas openpyxl numpy
    (Open Babel not needed here -- receptors are already prepared by
    prepare_targets.py. Vina itself is the external vina.exe
    binary, located via --vina.)
"""

from __future__ import annotations

import argparse
import contextlib
import json
import logging
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, Optional

import pandas as pd
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

DEFAULT_EXHAUSTIVENESS = 8
RESULTS_FILENAME = "results.xlsx"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class LigandSpec:
    name: str
    smiles: str

    @property
    def safe_name(self) -> str:
        return re.sub(r"[^A-Za-z0-9_-]", "_", self.name.strip())


@dataclass
class TargetRecord:
    target: str
    pdb_id: Optional[str]
    uniprot_id: Optional[str]
    chain: Optional[str]
    receptor_pdbqt: Path
    box_center: tuple[float, float, float]
    box_size: tuple[float, float, float]
    box_basis: str
    base_id: str  # matches the <target>_<id>_chain<C> naming from stage 1


class DockingResult(NamedTuple):
    docking_score: Optional[float]
    pose_path: Optional[str]
    error: Optional[str] = None


LIGAND_COLUMN_ALIASES = {"name", "ligand", "compound"}
SMILES_COLUMN_ALIASES = {"smiles", "smi"}


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def _find_column(columns: list[str], aliases: set[str]) -> Optional[str]:
    lowered = {c.lower().strip(): c for c in columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    return None


def read_ligands_excel(path: Path) -> list[LigandSpec]:
    df = pd.read_excel(path, dtype=str)
    columns = list(df.columns)
    col_name = _find_column(columns, LIGAND_COLUMN_ALIASES)
    col_smiles = _find_column(columns, SMILES_COLUMN_ALIASES)

    if col_name is None or col_smiles is None:
        raise SystemExit(
            f"Ligand Excel must have a name and smiles column. Found columns: {columns}"
        )

    ligands = []
    for _, row in df.iterrows():
        name = str(row[col_name]).strip()
        smiles = str(row[col_smiles]).strip()
        if not name or name.lower() == "nan" or not smiles or smiles.lower() == "nan":
            continue
        ligands.append(LigandSpec(name=name, smiles=smiles))
    return ligands


def load_prepared_targets(targets_dir: Path) -> list[TargetRecord]:
    """Scan the output of prepare_targets.py for every target
    that has both a pocket record and its receptor PDBQT. Targets that
    stage 1 marked FAILED have no pocket JSON and are skipped here
    automatically."""
    pockets_dir = targets_dir / "pockets"
    receptors_dir = targets_dir / "receptors"

    if not pockets_dir.exists():
        raise SystemExit(f"No pockets/ directory found under {targets_dir} -- run prepare_targets.py first")

    records: list[TargetRecord] = []
    for json_path in sorted(pockets_dir.glob("*_pocket.json")):
        base = json_path.stem[: -len("_pocket")]
        receptor_pdbqt = receptors_dir / f"{base}.pdbqt"
        if not receptor_pdbqt.exists():
            logger.info(f"[{base}] SKIPPED: pocket record found but no matching receptor PDBQT")
            continue

        data = json.loads(json_path.read_text())
        box = data["box"]
        records.append(
            TargetRecord(
                target=data.get("target", base),
                pdb_id=data.get("pdb_id"),
                uniprot_id=data.get("uniprot_id"),
                chain=data.get("chain"),
                receptor_pdbqt=receptor_pdbqt,
                box_center=tuple(box["center"]),
                box_size=tuple(box["size"]),
                box_basis=data.get("box_basis", "unknown"),
                base_id=base,
            )
        )

    return records


# ---------------------------------------------------------------------------
# Ligand prep (RDKit embed + Meeko PDBQT write) -- same approach as the
# original single-target docking script
# ---------------------------------------------------------------------------

def embed_to_sdf(smiles: str) -> Optional[Path]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xF00D) != 0:
        return None
    with contextlib.suppress(Exception):
        AllChem.MMFFOptimizeMolecule(mol)
    fd, name = tempfile.mkstemp(suffix=".sdf")
    os.close(fd)
    path = Path(name)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return path


def prepare_ligand(sdf_path: Path, out_path: Path) -> Path:
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = supplier[0] if supplier else None
    if mol is None:
        raise ValueError(f"RDKit failed to parse ligand: {sdf_path}")

    preparer = MoleculePreparation()
    setups = preparer.prepare(mol)
    pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(setups[0])
    out_path.write_text(pdbqt_string)
    return out_path


def prepare_ligand_pdbqt(ligand: LigandSpec, ligand_dir: Path) -> Optional[Path]:
    out_path = ligand_dir / f"{ligand.safe_name}.pdbqt"
    if out_path.exists():
        return out_path

    sdf = embed_to_sdf(ligand.smiles)
    if sdf is None:
        logger.info(f"[{ligand.name}] FAILED to embed/parse SMILES")
        return None
    try:
        prepare_ligand(sdf, out_path)
    except Exception as exc:
        logger.info(f"[{ligand.name}] FAILED ligand prep: {exc}")
        return None
    finally:
        sdf.unlink(missing_ok=True)

    return out_path


# ---------------------------------------------------------------------------
# Docking (vina.exe subprocess), with structured config/log/pose output --
# same pattern as the original single-target docking script
# ---------------------------------------------------------------------------

def dock(
    vina_exe: Path,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    run_id: str,
    exhaustiveness: int,
    seed: Optional[int],
    config_dir: Path,
    log_dir: Path,
    pose_dir: Path,
) -> DockingResult:
    out_pose = pose_dir / f"{run_id}_docked.pdbqt"
    log_path = log_dir / f"{run_id}.log"
    config_path = config_dir / f"{run_id}.txt"

    config_text = (
        f"receptor = {receptor_pdbqt}\n"
        f"ligand = {ligand_pdbqt}\n"
        f"center_x = {center[0]}\ncenter_y = {center[1]}\ncenter_z = {center[2]}\n"
        f"size_x = {size[0]}\nsize_y = {size[1]}\nsize_z = {size[2]}\n"
        f"exhaustiveness = {exhaustiveness}\n"
        f"seed = {seed if seed is not None else 'random'}\n"
        f"out = {out_pose}\n"
    )
    config_path.write_text(config_text)

    cmd = [
        str(vina_exe),
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--out", str(out_pose),
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    stdout = result.stdout or ""
    stderr = result.stderr or ""
    log_path.write_text(stdout + "\n" + stderr, encoding="utf-8")

    if result.returncode != 0:
        details = (stderr or stdout).strip()
        return DockingResult(
            None,
            None,
            details or f"Vina exited with code {result.returncode}.",
        )

    score = None
    for line in stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.?\d*)", line)
        if m:
            score = float(m.group(1))
            break

    if score is None:
        return DockingResult(
            None, None, "Vina completed without reporting a best-mode score."
        )
    if not out_pose.is_file():
        return DockingResult(
            None, None, "Vina reported a score but did not create the pose file."
        )
    return DockingResult(score, str(out_pose))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ligands_excel", type=Path, help="Excel file with ligand name + smiles columns")
    ap.add_argument("--targets-dir", type=Path, default=Path("targets_prepared"), help="Output dir from prepare_targets.py")
    ap.add_argument("--vina", type=Path, default=Path("vina.exe"), help="Path to the Vina executable")
    ap.add_argument("--outdir", type=Path, default=Path("docking_results"), help="Output directory for this stage")
    ap.add_argument("--exhaustiveness", type=int, default=DEFAULT_EXHAUSTIVENESS)
    ap.add_argument("--seed", type=int, default=None, help="Fix the Vina seed for reproducibility (default: random)")
    args = ap.parse_args()

    if not args.vina.exists():
        raise SystemExit(f"Vina executable not found at {args.vina.resolve()} -- pass --vina if it lives elsewhere")

    ligand_dir = args.outdir / "ligands"
    config_dir = args.outdir / "config"
    log_dir = args.outdir / "logs"
    pose_dir = args.outdir / "poses"
    for d in (ligand_dir, config_dir, log_dir, pose_dir):
        d.mkdir(parents=True, exist_ok=True)

    ligands = read_ligands_excel(args.ligands_excel)
    logger.info(f"Loaded {len(ligands)} ligand(s) from {args.ligands_excel}")

    targets = load_prepared_targets(args.targets_dir)
    logger.info(f"Loaded {len(targets)} prepared target(s) from {args.targets_dir}")

    # Prepare every ligand once, up front, and reuse across all targets.
    ligand_pdbqts: dict[str, Optional[Path]] = {}
    for ligand in ligands:
        ligand_pdbqts[ligand.name] = prepare_ligand_pdbqt(ligand, ligand_dir)

    rows = []
    for target in targets:
        row = {
            "target": target.target,
            "pdb_id": target.pdb_id,
            "uniprot_id": target.uniprot_id,
            "chain": target.chain,
            "box_basis": target.box_basis,
        }

        for ligand in ligands:
            ligand_pdbqt = ligand_pdbqts.get(ligand.name)
            if ligand_pdbqt is None:
                row[ligand.name] = None
                continue

            run_id = f"{target.base_id}_{ligand.safe_name}"
            result = dock(
                args.vina,
                target.receptor_pdbqt,
                ligand_pdbqt,
                target.box_center,
                target.box_size,
                run_id,
                args.exhaustiveness,
                args.seed,
                config_dir,
                log_dir,
                pose_dir,
            )
            row[ligand.name] = result.docking_score
            status = f"{result.docking_score:.2f}" if result.docking_score is not None else "FAILED"
            logger.info(f"[{target.target}] {ligand.name} -> {status} kcal/mol")

        rows.append(row)
        # Write results after every target so a crash midway doesn't
        # lose the docking runs already completed.
        pd.DataFrame(rows).to_excel(args.outdir / RESULTS_FILENAME, index=False)

    logger.info(f"\nDone. Results written to {args.outdir / RESULTS_FILENAME}")


if __name__ == "__main__":
    main()