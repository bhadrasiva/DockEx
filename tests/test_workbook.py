from io import BytesIO
from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from workbook import TargetInput, WorkbookValidationError, parse_workbook


def _workbook(ligands: pd.DataFrame, targets: pd.DataFrame) -> bytes:
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        ligands.to_excel(writer, sheet_name="Ligands", index=False)
        targets.to_excel(writer, sheet_name="Targets", index=False)
    return buffer.getvalue()


def test_parse_workbook_accepts_required_and_optional_fields():
    content = _workbook(
        pd.DataFrame({"name": ["test-ligand"], "smiles": ["C"]}),
        pd.DataFrame(
            {
                "target": ["EGFR"],
                "pdb_id": ["1m17"],
                "ligand_resname": ["aq4"],
                "pocket_residues": ["790; 791,797"],
                "chain": ["A"],
            }
        ),
    )

    parsed = parse_workbook(content)

    assert parsed.ligands[0].name == "test-ligand"
    assert parsed.targets[0].pdb_id == "1M17"
    assert parsed.targets[0].ligand_resname == "AQ4"
    assert parsed.targets[0].pocket_residues == (790, 791, 797)
    assert parsed.targets[0].chain == "A"


def test_parse_workbook_requires_both_sheets():
    buffer = BytesIO()
    pd.DataFrame({"name": ["test-ligand"], "smiles": ["C"]}).to_excel(
        buffer, sheet_name="Ligands", index=False
    )

    with pytest.raises(WorkbookValidationError, match="Targets"):
        parse_workbook(buffer.getvalue())


def test_parse_workbook_rejects_invalid_smiles_and_duplicate_names():
    content = _workbook(
        pd.DataFrame({"name": ["ligand", "LIGAND", "invalid"], "smiles": ["C", "CC", "not smiles"]}),
        pd.DataFrame({"target": ["target alpha"], "pdb_id": ["1M17"]}),
    )

    with pytest.raises(WorkbookValidationError) as error:
        parse_workbook(content)

    assert "duplicate ligand name" in str(error.value)
    assert "invalid SMILES" in str(error.value)


def test_parse_workbook_requires_target_identifier():
    content = _workbook(
        pd.DataFrame({"name": ["test-ligand"], "smiles": ["C"]}),
        pd.DataFrame({"target": ["target alpha"], "pdb_id": [None]}),
    )

    with pytest.raises(WorkbookValidationError, match="PDB ID or UniProt ID"):
        parse_workbook(content)


def test_parse_workbook_rejects_non_numeric_pocket_residues():
    content = _workbook(
        pd.DataFrame({"name": ["test-ligand"], "smiles": ["C"]}),
        pd.DataFrame(
            {"target": ["target alpha"], "pdb_id": ["1M17"], "pocket_residues": ["10;unknown"]}
        ),
    )

    with pytest.raises(WorkbookValidationError, match="positive residue numbers"):
        parse_workbook(content)


def test_parse_workbook_rejects_names_that_collide_after_filename_sanitizing():
    content = _workbook(
        pd.DataFrame({"name": ["ligand one", "ligand/one"], "smiles": ["C", "CC"]}),
        pd.DataFrame({"target": ["target alpha"], "pdb_id": ["1M17"]}),
    )

    with pytest.raises(WorkbookValidationError, match="collides"):
        parse_workbook(content)


def test_example_workbook_is_valid():
    example_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "dockex_sample_input.xlsx"
    )

    parsed = parse_workbook(example_path)

    assert len(parsed.ligands) == 2
    assert len(parsed.targets) == 1
    assert [(ligand.name, ligand.smiles) for ligand in parsed.ligands] == [
        ("Erlotinib", "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC"),
        ("Gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4"),
    ]
    assert parsed.targets == (
        TargetInput(
            name="EGFR",
            pdb_id="1M17",
            uniprot_id="P00533",
            ligand_resname="AQ4",
            pocket_residues=(),
            chain="A",
        ),
    )
    assert pd.ExcelFile(example_path).sheet_names == ["Ligands", "Targets"]
    assert list(pd.read_excel(example_path, sheet_name="Ligands").columns) == [
        "name",
        "smiles",
    ]
    assert list(pd.read_excel(example_path, sheet_name="Targets").columns) == [
        "target",
        "pdb_id",
        "uniprot_id",
        "ligand_resname",
        "pocket_residues",
        "chain",
    ]
