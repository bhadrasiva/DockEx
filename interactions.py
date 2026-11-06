"""Protein-ligand interaction analysis for docked poses using ProLIF."""

from __future__ import annotations

from pathlib import Path

from openbabel import openbabel
import pandas as pd
from rdkit import Chem
import prolif as plf


def _read_pose_as_rdkit(pose_path: Path) -> Chem.Mol:
    text = pose_path.read_text(encoding="utf-8")
    if "MODEL" in text:
        lines: list[str] = []
        in_first_model = False
        for line in text.splitlines():
            if line.startswith("MODEL"):
                if in_first_model:
                    break
                in_first_model = True
                continue
            if line.startswith("ENDMDL"):
                break
            if in_first_model:
                lines.append(line)
        text = "\n".join(lines)

    conversion = openbabel.OBConversion()
    conversion.SetInAndOutFormats("pdbqt", "mol")
    molecule = openbabel.OBMol()
    if not conversion.ReadString(molecule, text):
        raise ValueError(f"Open Babel could not parse docked pose {pose_path}")
    mol_block = conversion.WriteString(molecule)
    rdkit_molecule = Chem.MolFromMolBlock(mol_block, removeHs=False)
    if rdkit_molecule is None:
        raise ValueError(f"RDKit could not parse docked pose {pose_path}")
    return rdkit_molecule


def analyze_pose_interactions(receptor_path: Path, pose_path: Path) -> list[dict[str, str]]:
    """Return one row per detected interaction in a Vina pose."""
    protein = Chem.MolFromPDBFile(str(receptor_path), removeHs=False, sanitize=True)
    if protein is None:
        raise ValueError(f"RDKit could not parse prepared receptor {receptor_path}")

    ligand = _read_pose_as_rdkit(pose_path)
    protein_molecule = plf.Molecule.from_rdkit(protein)
    ligand_molecule = plf.Molecule.from_rdkit(ligand)

    fingerprint = plf.Fingerprint()
    fingerprint.run_from_iterable([ligand_molecule], protein_molecule, n_jobs=1, progress=False)
    interactions = fingerprint.to_dataframe()
    if interactions.empty:
        return []

    if isinstance(interactions.columns, pd.MultiIndex):
        interactions.columns = interactions.columns.droplevel(0)

    rows: list[dict[str, str]] = []
    for (residue, interaction), present in interactions.iloc[0].items():
        if bool(present):
            rows.append({"residue": str(residue), "interaction": str(interaction)})
    return rows
