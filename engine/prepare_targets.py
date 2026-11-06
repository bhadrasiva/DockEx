"""
prepare_targets.py

Batch target-preparation pipeline for structure-based docking.

Given an Excel sheet listing docking targets, this script:

    1. Fetches each structure (RCSB PDB, or AlphaFold DB as a fallback
       for targets with no experimental structure).
    2. Auto-detects the co-crystallized ligand and derives a per-target
       binding pocket (residues + a ligand-centroid docking box) --
       so every protein gets its own box instead of one box for all.
    3. Selects a single, best-resolved protein chain, resolves alternate
       conformations, and protonates the receptor (Open Babel, pH 7.4).
    4. Writes a receptor PDBQT ready for docking, a per-target pocket
       JSON record, and a master summary spreadsheet.

Targets with no bound ligand (true apo structures) fall back in this
order:
    1. a bound small-molecule ligand -> ligand-centroid box (preferred)
    2. a manually supplied pocket_residues column -> CA-centroid box
    3. neither -> a whole-protein box covering the entire selected
       chain, for blind docking. This is flagged BLIND_WHOLE_PROTEIN
       in the summary so it's clear the box wasn't pocket-derived.

--------------------------------------------------------------------
INPUT EXCEL SCHEMA
--------------------------------------------------------------------
One row per target. Column names are matched case-insensitively;
accepted aliases are listed for each.

Required (one of):
    target / name / protein            human-readable label, used in filenames

Required (one of, per row):
    pdb_id / pdb                       4-character RCSB accession -> fetched
                                        from files.rcsb.org
    uniprot_id / uniprot                UniProt accession -> fetched as an
                                        AlphaFold model (used when no
                                        experimental structure exists)

Optional:
    ligand_resname / ligand             force a specific 3-letter HETATM
                                        code as the pocket-defining ligand,
                                        for structures with multiple bound
                                        small molecules
    pocket_residues                     semicolon- or comma-separated
                                        residue numbers, e.g. "790;791;797".
                                        Used to build a manual CA-based box
                                        when a target has no usable bound
                                        ligand (apo structures, AlphaFold
                                        models).
    chain                               force a specific chain ID instead
                                        of letting the script pick one

--------------------------------------------------------------------
OUTPUT LAYOUT (under --outdir, default "targets_prepared")
--------------------------------------------------------------------
    raw/<PDBID or AF-UNIPROT>.cif       as-downloaded structure
    receptors/<target>_<id>_chain<C>.pdb          protein-only, single chain
    receptors/<target>_<id>_chain<C>_clean.pdb    altloc-resolved
    receptors/<target>_<id>_chain<C>.pdbqt        docking-ready receptor
    pockets/<target>_<id>_pocket.json             residues + box + metadata
    targets_summary.xlsx                          one row per target

--------------------------------------------------------------------
USAGE
--------------------------------------------------------------------
    python prepare_targets.py targets.xlsx --outdir targets_prepared

--------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------
    pip install biopython pandas openpyxl numpy requests
    Open Babel with Python bindings (conda install -c conda-forge openbabel,
    or however it is already installed in your environment -- this script
    uses the same `from openbabel import openbabel` API as the docking driver).
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
from Bio.PDB import PDBIO, MMCIFParser, NeighborSearch, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from openbabel import openbabel

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
ALPHAFOLD_CIF_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"

DEFAULT_POCKET_RADIUS = 6.0   # Angstrom, ligand -> pocket residue cutoff
DEFAULT_BOX_PADDING = 8.0     # Angstrom, added around ligand/residue extent
MIN_LIGAND_HEAVY_ATOMS = 6    # below this, treat as a crystallization additive
DEFAULT_PH = 7.4

WATER_RESNAMES = {"HOH", "DOD"}

ION_RESNAMES = {
    "NA", "CL", "MG", "MN", "ZN", "CA", "K", "FE", "FE2", "FE3", "CO",
    "NI", "CU", "CU1", "CD", "CS", "LI", "BR", "IOD", "SR", "BA",
}

CRYO_BUFFER_RESNAMES = {
    "GOL", "EDO", "PEG", "PG4", "1PE", "2PE", "MPD", "DMS", "ACT", "FMT",
    "SO4", "PO4", "TRS", "HEPES", "BME", "EPE", "IMD", "CIT", "BOG",
    "P6G", "PGE", "MES", "BCT", "NO3", "ACY", "TAR", "CO3", "DTT",
    "MRD", "TLA", "OXL", "UNL", "UNX",
}

MODIFIED_RESIDUE_RESNAMES = {
    "SEP", "TPO", "PTR", "MSE", "CSO", "CSD", "KCX", "MLY", "M3L",
    "CME", "OCS", "PCA", "HYP", "SAH", "SAM",
}

EXCLUDED_LIGAND_RESNAMES = (
    WATER_RESNAMES | ION_RESNAMES | CRYO_BUFFER_RESNAMES | MODIFIED_RESIDUE_RESNAMES
)

RAW_DIR_NAME = "raw"
RECEPTOR_DIR_NAME = "receptors"
POCKET_DIR_NAME = "pockets"
SUMMARY_FILENAME = "targets_summary.xlsx"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class TargetSpec:
    name: str
    pdb_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    ligand_override: Optional[str] = None
    manual_pocket_residues: list[int] = field(default_factory=list)
    forced_chain: Optional[str] = None

    @property
    def safe_name(self) -> str:
        return re.sub(r"[^A-Za-z0-9_-]", "_", self.name.strip())

    @property
    def source_id(self) -> str:
        """Identifier used in output filenames."""
        return (self.pdb_id or self.uniprot_id or "unknown").upper()


COLUMN_ALIASES = {
    "name": {"target", "name", "protein"},
    "pdb_id": {"pdb_id", "pdb", "pdbid"},
    "uniprot_id": {"uniprot_id", "uniprot", "uniprotid"},
    "ligand_override": {"ligand_resname", "ligand", "ligand_override"},
    "pocket_residues": {"pocket_residues", "manual_pocket_residues"},
    "chain": {"chain", "chain_id"},
}


def _find_column(columns: list[str], aliases: set[str]) -> Optional[str]:
    lowered = {c.lower().strip(): c for c in columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    return None


def read_targets_excel(path: Path) -> list[TargetSpec]:
    """Read the input Excel sheet into a list of TargetSpec, tolerating
    a handful of common column-name spellings (see module docstring)."""
    df = pd.read_excel(path, dtype=str)
    columns = list(df.columns)

    col_name = _find_column(columns, COLUMN_ALIASES["name"])
    col_pdb = _find_column(columns, COLUMN_ALIASES["pdb_id"])
    col_uniprot = _find_column(columns, COLUMN_ALIASES["uniprot_id"])
    col_ligand = _find_column(columns, COLUMN_ALIASES["ligand_override"])
    col_pocket = _find_column(columns, COLUMN_ALIASES["pocket_residues"])
    col_chain = _find_column(columns, COLUMN_ALIASES["chain"])

    if col_name is None or (col_pdb is None and col_uniprot is None):
        raise SystemExit(
            "Input Excel must have a target/name column and at least one of "
            f"pdb_id or uniprot_id. Found columns: {columns}"
        )

    targets: list[TargetSpec] = []
    for _, row in df.iterrows():
        name = str(row[col_name]).strip()
        if not name or name.lower() == "nan":
            continue

        pdb_id = str(row[col_pdb]).strip() if col_pdb and pd.notna(row[col_pdb]) else None
        uniprot_id = (
            str(row[col_uniprot]).strip() if col_uniprot and pd.notna(row[col_uniprot]) else None
        )
        ligand_override = (
            str(row[col_ligand]).strip().upper()
            if col_ligand and pd.notna(row[col_ligand])
            else None
        )
        forced_chain = (
            str(row[col_chain]).strip() if col_chain and pd.notna(row[col_chain]) else None
        )

        pocket_residues: list[int] = []
        if col_pocket and pd.notna(row[col_pocket]):
            raw = str(row[col_pocket])
            for token in re.split(r"[;,\s]+", raw.strip()):
                if token.isdigit():
                    pocket_residues.append(int(token))

        if not pdb_id and not uniprot_id:
            logger.info(f"[{name}] SKIPPED: no pdb_id or uniprot_id given")
            continue

        targets.append(
            TargetSpec(
                name=name,
                pdb_id=pdb_id or None,
                uniprot_id=uniprot_id or None,
                ligand_override=ligand_override,
                manual_pocket_residues=pocket_residues,
                forced_chain=forced_chain,
            )
        )

    return targets


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

def fetch_structure(spec: TargetSpec, raw_dir: Path) -> tuple[Path, str]:
    """Download the structure for a target. Tries RCSB first if a pdb_id
    is given, otherwise falls back to an AlphaFold model. Returns the
    local mmCIF path and a source tag ("rcsb" or "alphafold")."""
    raw_dir.mkdir(parents=True, exist_ok=True)

    if spec.pdb_id:
        out_path = raw_dir / f"{spec.pdb_id.upper()}.cif"
        if not out_path.exists():
            url = RCSB_CIF_URL.format(pdb_id=spec.pdb_id.upper())
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
            out_path.write_text(resp.text)
        return out_path, "rcsb"

    if spec.uniprot_id:
        out_path = raw_dir / f"AF-{spec.uniprot_id.upper()}.cif"
        if not out_path.exists():
            url = ALPHAFOLD_CIF_URL.format(uniprot_id=spec.uniprot_id.upper())
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
            out_path.write_text(resp.text)
        return out_path, "alphafold"

    raise ValueError(f"[{spec.name}] no pdb_id or uniprot_id to fetch")


# ---------------------------------------------------------------------------
# Structure parsing helpers
# ---------------------------------------------------------------------------

def load_structure(cif_path: Path):
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(cif_path.stem, str(cif_path))


def remap_multichar_chains(structure) -> None:
    """mmCIF entries may use multi-character chain IDs (e.g. author chain
    'AAA'); the legacy PDB format Biopython writes only supports single
    characters, so remap deterministically before doing anything else."""
    pool = iter(string.ascii_uppercase + string.ascii_lowercase + string.digits)
    used = set()
    for model in structure:
        for chain in model:
            if len(chain.id) != 1:
                new_id = next(c for c in pool if c not in used)
                used.add(new_id)
                chain.id = new_id
            else:
                used.add(chain.id)


def _is_ligand_candidate(residue) -> bool:
    hetfield = residue.id[0]
    if hetfield == " ":
        return False  # standard polymer residue
    resname = residue.resname.strip().upper()
    if resname in EXCLUDED_LIGAND_RESNAMES:
        return False
    heavy_atoms = sum(1 for atom in residue if atom.element != "H")
    return heavy_atoms >= MIN_LIGAND_HEAVY_ATOMS


def _chain_protein_length(chain) -> int:
    return sum(1 for res in chain if is_aa(res, standard=True))


def identify_ligand_and_chain(
    structure, ligand_override: Optional[str], forced_chain: Optional[str]
):
    """Pick the ligand instance that best represents the co-crystallized
    small molecule, and the chain it sits in. Returns (residue_or_None,
    chain_id). residue is None for apo structures / AlphaFold models."""
    model = next(iter(structure))

    candidates = []
    for chain in model:
        if forced_chain and chain.id != forced_chain:
            continue
        for residue in chain:
            if _is_ligand_candidate(residue):
                if ligand_override and residue.resname.strip().upper() != ligand_override:
                    continue
                heavy_atoms = sum(1 for atom in residue if atom.element != "H")
                candidates.append((chain.id, residue, heavy_atoms))

    if candidates:
        # Prefer the chain with the most complete protein, then the
        # largest ligand instance, so we land on a well-resolved protomer.
        candidates.sort(
            key=lambda c: (_chain_protein_length(model[c[0]]), c[2]), reverse=True
        )
        chain_id, ligand_residue, _ = candidates[0]
        return ligand_residue, chain_id

    # No ligand found: pick chain by protein completeness alone.
    if forced_chain:
        return None, forced_chain
    best_chain = max(model, key=_chain_protein_length)
    return None, best_chain.id


# ---------------------------------------------------------------------------
# Pocket + box derivation
# ---------------------------------------------------------------------------

def derive_pocket_from_ligand(structure, ligand_residue, radius: float) -> dict:
    model = next(iter(structure))
    ns = NeighborSearch(list(model.get_atoms()))

    seen = set()
    pocket_residues = []
    for atom in ligand_residue:
        for near_atom in ns.search(atom.coord, radius):
            res = near_atom.get_parent()
            if not is_aa(res, standard=True):
                continue
            key = (res.get_parent().id, res.id[1])
            if key not in seen:
                seen.add(key)
                pocket_residues.append(
                    {"chain": res.get_parent().id, "resnum": res.id[1], "resname": res.resname}
                )

    coords = np.array([atom.coord for atom in ligand_residue])
    center = coords.mean(axis=0)
    size = np.ptp(coords, axis=0) + DEFAULT_BOX_PADDING

    return {
        "pocket_residues": pocket_residues,
        "box": {
            "center": [round(float(v), 2) for v in center],
            "size": [round(float(v), 2) for v in size],
        },
        "box_basis": "ligand_centroid",
    }


def derive_box_whole_protein(structure, chain_id: str, padding: float) -> dict:
    """Fallback for targets with no bound ligand and no manually supplied
    pocket residues: box around every atom of the selected chain, for
    blind docking rather than a guessed-at pocket."""
    model = next(iter(structure))
    chain = model[chain_id]

    coords = [atom.coord for res in chain if is_aa(res, standard=True) for atom in res]
    if not coords:
        raise ValueError(f"No protein atoms found in chain {chain_id} to build a whole-protein box")

    coords = np.array(coords)
    center = coords.mean(axis=0)
    size = np.ptp(coords, axis=0) + padding

    return {
        "pocket_residues": [],  # not a pocket -- whole chain is in the box
        "box": {
            "center": [round(float(v), 2) for v in center],
            "size": [round(float(v), 2) for v in size],
        },
        "box_basis": "whole_protein_blind",
    }


def derive_pocket_from_residue_list(
    structure, chain_id: str, residue_numbers: list[int], padding: float
) -> dict:
    model = next(iter(structure))
    chain = model[chain_id]

    coords = []
    pocket_residues = []
    for resnum in residue_numbers:
        for res in chain:
            if res.id[1] == resnum and "CA" in res:
                coords.append(res["CA"].coord)
                pocket_residues.append(
                    {"chain": chain_id, "resnum": resnum, "resname": res.resname}
                )

    if not coords:
        raise ValueError(f"None of the supplied pocket residues {residue_numbers} found in chain {chain_id}")

    coords = np.array(coords)
    center = coords.mean(axis=0)
    size = np.ptp(coords, axis=0) + padding

    return {
        "pocket_residues": pocket_residues,
        "box": {
            "center": [round(float(v), 2) for v in center],
            "size": [round(float(v), 2) for v in size],
        },
        "box_basis": "manual_residue_ca",
    }


# ---------------------------------------------------------------------------
# Receptor extraction + cleaning + protonation
# ---------------------------------------------------------------------------

class _ChainProteinOnly(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue) -> bool:
        return bool(is_aa(residue, standard=True))


def write_single_chain_protein(structure, chain_id: str, out_pdb: Path) -> Path:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), _ChainProteinOnly(chain_id))
    return out_pdb


def clean_altlocs(pdb_in: Path, pdb_out: Path) -> Path:
    """Resolve alternate conformations: keep the highest-occupancy altloc only."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_in))

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in list(residue):
                    if atom.is_disordered():
                        atom.disordered_select(
                            max(
                                atom.disordered_get_id_list(),
                                key=lambda k: atom.disordered_get(k).get_occupancy(),
                            )
                        )

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_out))
    return pdb_out


def protonate_to_pdbqt(clean_pdb: Path, out_pdbqt: Path, ph: float = DEFAULT_PH) -> Path:
    """Strip existing hydrogens, add fresh ones at the given pH, and write
    a rigid receptor PDBQT with Gasteiger charges via Open Babel."""
    in_conv = openbabel.OBConversion()
    in_conv.SetInFormat("pdb")

    mol = openbabel.OBMol()
    if not in_conv.ReadFile(mol, str(clean_pdb)):
        raise RuntimeError(f"Open Babel failed to read {clean_pdb}")

    mol.DeleteHydrogens()
    mol.AddHydrogens(False, True, ph)  # polaronly=False, correctForPH=True

    out_conv = openbabel.OBConversion()
    out_conv.SetOutFormat("pdbqt")
    out_conv.AddOption("r", openbabel.OBConversion.OUTOPTIONS)  # rigid receptor

    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    if not out_conv.WriteFile(mol, str(out_pdbqt)):
        raise RuntimeError(f"Open Babel failed to write {out_pdbqt}")

    return out_pdbqt


# ---------------------------------------------------------------------------
# Per-target pipeline
# ---------------------------------------------------------------------------

def process_target(spec: TargetSpec, out_dir: Path, radius: float, padding: float, ph: float) -> dict:
    raw_dir = out_dir / RAW_DIR_NAME
    receptor_dir = out_dir / RECEPTOR_DIR_NAME
    pocket_dir = out_dir / POCKET_DIR_NAME

    row: dict = {
        "target": spec.name,
        "pdb_id": spec.pdb_id,
        "uniprot_id": spec.uniprot_id,
        "source": None,
        "chain": None,
        "ligand_resname": None,
        "pocket_residue_count": None,
        "box_center": None,
        "box_size": None,
        "box_basis": None,
        "receptor_pdbqt": None,
        "status": "FAILED",
        "notes": "",
    }

    try:
        cif_path, source = fetch_structure(spec, raw_dir)
        row["source"] = source

        structure = load_structure(cif_path)
        remap_multichar_chains(structure)

        ligand_residue, chain_id = identify_ligand_and_chain(
            structure, spec.ligand_override, spec.forced_chain
        )
        row["chain"] = chain_id

        pocket_info = None
        if ligand_residue is not None:
            row["ligand_resname"] = ligand_residue.resname.strip()
            pocket_info = derive_pocket_from_ligand(structure, ligand_residue, radius)
        elif spec.manual_pocket_residues:
            pocket_info = derive_pocket_from_residue_list(
                structure, chain_id, spec.manual_pocket_residues, padding
            )
        else:
            pocket_info = derive_box_whole_protein(structure, chain_id, padding)
            row["status"] = "BLIND_WHOLE_PROTEIN"
            row["notes"] = (
                "No bound ligand detected and no pocket_residues supplied; "
                "box covers the entire chain (blind docking) rather than a derived pocket."
            )

        # Receptor extraction + cleaning + protonation happens regardless.
        base = f"{spec.safe_name}_{row['pdb_id'] or row['uniprot_id']}_chain{chain_id}"
        raw_chain_pdb = receptor_dir / f"{base}.pdb"
        clean_pdb = receptor_dir / f"{base}_clean.pdb"
        pdbqt_path = receptor_dir / f"{base}.pdbqt"

        write_single_chain_protein(structure, chain_id, raw_chain_pdb)
        clean_altlocs(raw_chain_pdb, clean_pdb)
        protonate_to_pdbqt(clean_pdb, pdbqt_path, ph)

        row["receptor_pdbqt"] = str(pdbqt_path)
        if row["status"] == "FAILED":
            row["status"] = "OK"

        if pocket_info is not None:
            row["pocket_residue_count"] = len(pocket_info["pocket_residues"])
            row["box_center"] = pocket_info["box"]["center"]
            row["box_size"] = pocket_info["box"]["size"]
            row["box_basis"] = pocket_info["box_basis"]

            pocket_dir.mkdir(parents=True, exist_ok=True)
            pocket_record = {
                "target": spec.name,
                "pdb_id": spec.pdb_id,
                "uniprot_id": spec.uniprot_id,
                "source": source,
                "chain": chain_id,
                "ligand_resname": row["ligand_resname"],
                **pocket_info,
            }
            (pocket_dir / f"{base}_pocket.json").write_text(json.dumps(pocket_record, indent=2))

        logger.info(f"[{spec.name}] {row['status']} (chain {chain_id}, ligand {row['ligand_resname']})")

    except Exception as exc:  # noqa: BLE001 - keep the batch going on a per-target failure
        row["status"] = "FAILED"
        row["notes"] = str(exc)
        logger.info(f"[{spec.name}] FAILED: {exc}")

    return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("excel", type=Path, help="Input Excel file listing targets (see module docstring)")
    ap.add_argument("--outdir", type=Path, default=Path("targets_prepared"), help="Output directory")
    ap.add_argument("--radius", type=float, default=DEFAULT_POCKET_RADIUS, help="Ligand -> pocket residue cutoff (A)")
    ap.add_argument("--padding", type=float, default=DEFAULT_BOX_PADDING, help="Box padding around ligand/residue extent (A)")
    ap.add_argument("--ph", type=float, default=DEFAULT_PH, help="Protonation pH for receptor prep")
    args = ap.parse_args()

    targets = read_targets_excel(args.excel)
    logger.info(f"Loaded {len(targets)} target(s) from {args.excel}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for spec in targets:
        row = process_target(spec, args.outdir, args.radius, args.padding, args.ph)
        rows.append(row)
        # Write the summary after every target so a crash midway
        # doesn't lose progress on the ones already done.
        pd.DataFrame(rows).to_excel(args.outdir / SUMMARY_FILENAME, index=False)

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_blind = sum(1 for r in rows if r["status"] == "BLIND_WHOLE_PROTEIN")
    n_failed = sum(1 for r in rows if r["status"] == "FAILED")
    logger.info(
        f"\nDone. {n_ok} ready with a derived pocket, {n_blind} ready with a whole-protein "
        f"blind box, {n_failed} failed. Summary: {args.outdir / SUMMARY_FILENAME}"
    )


if __name__ == "__main__":
    main()
