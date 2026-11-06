"""Validation and parsing for the app's two-sheet Excel input workbook."""

from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
import re
from typing import BinaryIO

import pandas as pd
from rdkit import Chem


@dataclass(frozen=True)
class LigandInput:
    name: str
    smiles: str


@dataclass(frozen=True)
class TargetInput:
    name: str
    pdb_id: str | None = None
    uniprot_id: str | None = None
    ligand_resname: str | None = None
    pocket_residues: tuple[int, ...] = ()
    chain: str | None = None


@dataclass(frozen=True)
class ParsedWorkbook:
    ligands: tuple[LigandInput, ...]
    targets: tuple[TargetInput, ...]


class WorkbookValidationError(ValueError):
    """Raised when a workbook does not satisfy the input schema."""


def _column(columns: list[str], aliases: tuple[str, ...]) -> str | None:
    normalized = {str(column).strip().lower(): str(column) for column in columns}
    return next((normalized[alias] for alias in aliases if alias in normalized), None)


def _cell_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", value.strip())


def _read_frame(workbook: BytesIO | BinaryIO | Path | str, sheet_name: str) -> pd.DataFrame:
    return pd.read_excel(workbook, sheet_name=sheet_name, dtype=object)


def parse_workbook(workbook: bytes | BytesIO | BinaryIO | Path | str) -> ParsedWorkbook:
    """Parse and validate the required ``Ligands`` and ``Targets`` sheets."""
    source: BytesIO | BinaryIO | Path | str
    source = BytesIO(workbook) if isinstance(workbook, bytes) else workbook

    try:
        excel = pd.ExcelFile(source)
    except Exception as exc:
        raise WorkbookValidationError(f"Could not read the uploaded Excel workbook: {exc}") from exc

    sheet_names = {name.strip().lower(): name for name in excel.sheet_names}
    missing_sheets = [name for name in ("ligands", "targets") if name not in sheet_names]
    if missing_sheets:
        raise WorkbookValidationError(
            "Workbook must contain sheets named 'Ligands' and 'Targets'. "
            f"Missing: {', '.join(missing_sheets)}."
        )

    ligand_frame = _read_frame(source, sheet_names["ligands"])
    target_frame = _read_frame(source, sheet_names["targets"])
    errors: list[str] = []

    ligand_name_col = _column(list(ligand_frame.columns), ("name", "ligand", "compound"))
    smiles_col = _column(list(ligand_frame.columns), ("smiles", "smi"))
    if ligand_name_col is None or smiles_col is None:
        raise WorkbookValidationError(
            "The Ligands sheet must have a name/ligand/compound column and a smiles/smi column."
        )

    ligands: list[LigandInput] = []
    ligand_names: set[str] = set()
    safe_ligand_names: set[str] = set()
    for row_index, row in ligand_frame.iterrows():
        name = _cell_text(row[ligand_name_col])
        smiles = _cell_text(row[smiles_col])
        if not name and not smiles:
            continue
        excel_row = row_index + 2
        if not name or not smiles:
            errors.append(f"Ligands row {excel_row}: both name and SMILES are required.")
            continue
        key = name.casefold()
        if key in ligand_names:
            errors.append(f"Ligands row {excel_row}: duplicate ligand name {name!r}.")
            continue
        ligand_names.add(key)
        safe_name = _safe_name(name)
        if safe_name in safe_ligand_names:
            errors.append(
                f"Ligands row {excel_row}: name {name!r} collides with another output filename."
            )
            continue
        safe_ligand_names.add(safe_name)
        if Chem.MolFromSmiles(smiles) is None:
            errors.append(f"Ligands row {excel_row}: invalid SMILES for {name!r}.")
            continue
        ligands.append(LigandInput(name=name, smiles=smiles))

    target_name_col = _column(list(target_frame.columns), ("target", "name", "protein"))
    pdb_col = _column(list(target_frame.columns), ("pdb_id", "pdb", "pdbid"))
    uniprot_col = _column(list(target_frame.columns), ("uniprot_id", "uniprot", "uniprotid"))
    ligand_col = _column(list(target_frame.columns), ("ligand_resname", "ligand", "ligand_override"))
    pocket_col = _column(list(target_frame.columns), ("pocket_residues", "manual_pocket_residues"))
    chain_col = _column(list(target_frame.columns), ("chain", "chain_id"))

    if target_name_col is None or (pdb_col is None and uniprot_col is None):
        raise WorkbookValidationError(
            "The Targets sheet must have a target/name/protein column and a pdb_id/pdb "
            "or uniprot_id/uniprot column."
        )

    targets: list[TargetInput] = []
    target_names: set[str] = set()
    safe_target_names: set[str] = set()
    for row_index, row in target_frame.iterrows():
        name = _cell_text(row[target_name_col])
        pdb_id = _cell_text(row[pdb_col]).upper() if pdb_col else ""
        uniprot_id = _cell_text(row[uniprot_col]).upper() if uniprot_col else ""
        ligand_resname = _cell_text(row[ligand_col]).upper() if ligand_col else ""
        chain = _cell_text(row[chain_col]) if chain_col else ""
        raw_pocket = _cell_text(row[pocket_col]) if pocket_col else ""

        if not any((name, pdb_id, uniprot_id, ligand_resname, chain, raw_pocket)):
            continue
        excel_row = row_index + 2
        if not name:
            errors.append(f"Targets row {excel_row}: target name is required.")
            continue
        if not pdb_id and not uniprot_id:
            errors.append(f"Targets row {excel_row}: provide a PDB ID or UniProt ID.")
            continue
        if pdb_id and (len(pdb_id) != 4 or not re.fullmatch(r"[A-Z0-9]{4}", pdb_id)):
            errors.append(f"Targets row {excel_row}: PDB ID {pdb_id!r} must be four letters/digits.")
            continue
        if uniprot_id and not re.fullmatch(r"[A-Z0-9]{6,10}", uniprot_id):
            errors.append(f"Targets row {excel_row}: invalid UniProt ID {uniprot_id!r}.")
            continue

        key = name.casefold()
        if key in target_names:
            errors.append(f"Targets row {excel_row}: duplicate target name {name!r}.")
            continue
        target_names.add(key)
        safe_name = _safe_name(name)
        if safe_name in safe_target_names:
            errors.append(
                f"Targets row {excel_row}: name {name!r} collides with another output filename."
            )
            continue
        safe_target_names.add(safe_name)

        pocket_residues: list[int] = []
        if raw_pocket:
            tokens = [token for token in re.split(r"[;,\s]+", raw_pocket) if token]
            invalid = [token for token in tokens if not token.isdigit() or int(token) < 1]
            if invalid:
                errors.append(
                    f"Targets row {excel_row}: pocket_residues must be positive residue numbers; "
                    f"invalid values: {', '.join(invalid)}."
                )
                continue
            pocket_residues = [int(token) for token in tokens]

        targets.append(
            TargetInput(
                name=name,
                pdb_id=pdb_id or None,
                uniprot_id=uniprot_id or None,
                ligand_resname=ligand_resname or None,
                pocket_residues=tuple(pocket_residues),
                chain=chain or None,
            )
        )

    if not ligands:
        errors.append("The Ligands sheet must contain at least one valid ligand.")
    if not targets:
        errors.append("The Targets sheet must contain at least one valid target.")
    if errors:
        raise WorkbookValidationError("\n".join(errors))

    return ParsedWorkbook(ligands=tuple(ligands), targets=tuple(targets))
