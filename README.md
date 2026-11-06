# DockEx

**DockEx — Molecular Docking & Virtual Screening** is a local Streamlit app for
preparing protein structures, docking every ligand against every target in an
Excel workbook, inspecting ProLIF interaction fingerprints, and retaining the
inputs, reports, poses, and logs for each run.

DockEx uses independent target-specific docking boxes. It is **not** an
ensemble-docking or ensemble-scoring model. Vina scores are approximate
computational scores; compare compounds primarily within a single target, not
as directly comparable affinities across different proteins.

## What it does

1. Validates a workbook containing `Ligands` and `Targets` sheets.
2. Retrieves structures from RCSB PDB or AlphaFold DB, selects/prepares a
   protein chain, and derives a pocket from a bound ligand, user-supplied
   residue numbers, or (when neither is available) a clearly flagged
   whole-protein blind-docking box.
3. Prepares each ligand and runs AutoDock Vina for every usable ligand-target
   pair. A failed target, ligand, or docking pair is recorded without
   discarding results from other pairs.
4. Analyzes the first returned Vina pose with ProLIF and records detected
   residue-level interaction classifications. ProLIF failures are reported
   separately and do not erase scores or poses.
5. Saves a durable local run folder with progress reports, scores, pair-level
   statuses, structures, poses, configs, logs, and provenance metadata.
6. Optionally lets a user discuss a selected run with a read-only Gemini
   assistant after explicit consent to share run evidence with Google.

## Workbook

Download [`examples/dockex_sample_input.xlsx`](./examples/dockex_sample_input.xlsx)
from the app or use it as a template. It contains exactly two sheets:

### `Ligands`

| name | smiles |
|---|---|
| Erlotinib | `COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC` |
| Gefitinib | `COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4` |

Accepted name headers: `name`, `ligand`, or `compound`; SMILES headers:
`smiles` or `smi`. Names must be unique and SMILES must be valid.

### `Targets`

| target | pdb_id | uniprot_id | ligand_resname | pocket_residues | chain |
|---|---|---|---|---|---|
| EGFR | 1M17 | P00533 | AQ4 | | A |

Each row needs a target name and either a four-character PDB ID or a UniProt
ID. Optional columns select a ligand residue, provide comma/semicolon-separated
positive pocket residue numbers, or force a chain. A bound ligand is preferred
for pocket derivation; supplied residues are used when no suitable bound ligand
is available; otherwise the run is flagged for blind whole-protein docking.

The sample is an EGFR inhibitor screening example: 1M17 is an experimentally
determined EGFR structure with co-crystallized erlotinib in residue AQ4 on chain
A. The target entry specifies AQ4 so DockEx can use that ligand to define the
docking pocket. Erlotinib is included as the co-crystallized reference compound,
alongside gefitinib as a second EGFR inhibitor. This is a scientifically
relevant example input, not a validated docking benchmark; predicted scores
remain approximate and require scientific interpretation.

## Installation and launch

Python dependencies are installed once per environment; cloning alone does not
install them. AutoDock Vina is a separate executable and must also be installed
or downloaded for the operating system.

1. Install Python 3.10 and create an environment from this directory:

   ```powershell
   py -3.10 -m venv .venv
   .\.venv\Scripts\python.exe -m pip install --upgrade pip
   .\.venv\Scripts\python.exe -m pip install -r requirements.txt
   ```

2. Install AutoDock Vina from the
   [official AutoDock Vina project](https://github.com/ccsb-scripps/AutoDock-Vina/releases).
   Set `VINA_EXECUTABLE` to the executable path, or put `vina.exe` in the
   DockEx directory. For example, in PowerShell:

   ```powershell
   $env:VINA_EXECUTABLE = "C:\tools\vina.exe"
   ```

3. Launch with the same environment's Python:

   ```powershell
   .\.venv\Scripts\python.exe -m streamlit run app.py
   ```

   Use the environment's Python for both installation and launch. A bare
   `streamlit` command can resolve to a different Python installation and
   produce missing-module errors.

4. Upload the workbook, review the parsed inputs and L×T workload, set Vina
   search exhaustiveness/seed, and start the run. Network access is required
   to retrieve target structures. Runtime grows with the number of ligand-target
   pairs.

## Results and retained data

The app shows preparation status, score heatmaps and within-target rankings,
ProLIF interaction classifications, pair-level failures, and downloadable
reports. Each run is stored beneath the local data directory:

```text
<data directory>/
  history.sqlite3
  runs/<run UUID>/
    input/workbook.xlsx
    targets_prepared/       # retrieved structures, prepared receptors, pockets
    ligands/                # prepared ligands
    config/                 # Vina settings per pair
    logs/                   # Vina output per pair
    poses/                  # docked poses
    application.log
    docking_scores.csv
    docking_scores.xlsx
    docking_pairs.csv       # one status row per requested ligand-target pair
    interactions.csv
    interaction_analysis_errors.csv
    targets_summary.xlsx
    run_metadata.json       # versions, settings, input hash, target/pocket data
```

By default, the per-user data directory is
`%USERPROFILE%\.molecular_docking` on Windows and `~/.molecular_docking` on
macOS/Linux. Set `DOCKEX_DATA_DIR` to use another location; the legacy
`MOLECULAR_DOCKING_DATA_DIR` environment variable is also recognized. Run
history and chat are stored in SQLite; result files are retained in the
corresponding run folder, not embedded in the database. Deleting a run in the
app removes its local run folder and associated chat history.

The local database and artifacts are not encrypted. Do not retain sensitive
input data in an untrusted or shared account. Structure retrieval and Gemini
are external network services; the app does not send data to Gemini unless the
user supplies an API key, accepts the disclosure, and submits a question.

## Optional Gemini discussion

The post-run assistant uses the `gemini-3.5-flash-lite` model. It is read-only:
it can discuss saved run evidence, but it does not launch docking jobs or edit
results. The API key is used only in the active app session and is not stored
in SQLite. After consent, the user's question, recent chat, ligand names and
SMILES, target identifiers, scores, and detected interactions are sent to
Google Gemini. API access, quota, and charges are subject to the user's Google
AI account and Google's terms.

## Architecture

```text
Excel upload
  └─ workbook.py: schema validation and typed ligand/target records
      └─ pipeline.py: per-target preparation, per-ligand preparation,
                      L×T Vina calls, ProLIF, reports, provenance
          ├─ engine/prepare_targets.py: structure retrieval and receptor/pocket prep
          ├─ engine/dock_targets.py: ligand prep, Vina calls, per-pair logs/poses
          └─ interactions.py: first-pose conversion and ProLIF fingerprints
              └─ app.py: progress, history, plots, downloads, optional Q&A
                  ├─ chat_store.py: local SQLite run index/chat history
                  └─ gemini_agent.py: consent-gated, read-only Gemini requests
```

Vina score values and ProLIF contacts are hypotheses for prioritization and
follow-up, not experimental measurements. DockEx does not calculate binding
free energies, establish binding, or make cross-target scores equivalent.
Structure selection, protonation, pocket definition, and pose quality require
scientific review.

## Tests

Run the local tests using the same environment:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-dev.txt
.\.venv\Scripts\python.exe -m pytest tests -q
```

The tests cover workbook validation, all-pairs execution and failure isolation,
Vina result handling, local history/artifact deletion, and Gemini request
construction without making live network or API calls.

## Repository hygiene

Do not commit local databases, run outputs, credentials, virtual environments,
or downloaded Vina binaries. `.gitignore` excludes these local artifacts.
DockEx is distributed under the MIT License; see [`LICENSE`](./LICENSE) for
the terms and copyright notice.
