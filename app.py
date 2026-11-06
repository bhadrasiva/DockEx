from __future__ import annotations

from io import BytesIO
import json
import os
from pathlib import Path
import shutil
import uuid
import zipfile

import pandas as pd
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw
import streamlit as st

from chat_store import (
    add_chat_turn,
    artifact_directory,
    database_path,
    delete_run,
    get_chat_history,
    initialize_database,
    list_runs,
    load_run,
    save_run,
)
from gemini_agent import MODEL_ID, answer_question
from pipeline import run_docking_job
from workbook import WorkbookValidationError, parse_workbook

APP_DIR = Path(__file__).resolve().parent
EXAMPLE_WORKBOOK = APP_DIR / "examples" / "dockex_sample_input.xlsx"


def _find_vina() -> Path | None:
    configured = os.environ.get("VINA_EXECUTABLE")
    candidates = [Path(configured)] if configured else []
    located = shutil.which("vina")
    if located:
        candidates.append(Path(located))
    candidates.extend((APP_DIR / "vina.exe", APP_DIR / "vina"))
    return next((candidate.resolve() for candidate in candidates if candidate.is_file()), None)


def _zip_job(output_dir: Path) -> bytes:
    archive = BytesIO()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zipped:
        for path in output_dir.rglob("*"):
            if path.is_file() and not path.is_symlink():
                zipped.write(path, path.relative_to(output_dir))
    return archive.getvalue()


def _frame_records(frame: pd.DataFrame) -> list[dict]:
    return json.loads(frame.to_json(orient="records"))


def _xlsx_bytes(frame: pd.DataFrame, sheet_name: str = "Results") -> bytes:
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        frame.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    return buffer.getvalue()


def _load_saved_payload(db_path: Path, run_id: str) -> tuple[dict | None, Path | None]:
    payload = load_run(db_path, run_id)
    if payload is None:
        return None, None
    artifacts = artifact_directory(db_path, payload, run_id)
    if artifacts is None or not artifacts.is_dir():
        payload.setdefault("metadata", {"status": "LEGACY"})
        return payload, None

    report_files = {
        "scores": ("docking_scores.csv", pd.read_csv),
        "interactions": ("interactions.csv", pd.read_csv),
        "interaction_errors": ("interaction_analysis_errors.csv", pd.read_csv),
        "pair_results": ("docking_pairs.csv", pd.read_csv),
        "target_preparations": ("targets_summary.xlsx", pd.read_excel),
    }
    for key, (filename, reader) in report_files.items():
        report_path = artifacts / filename
        if report_path.is_file():
            payload[key] = _frame_records(reader(report_path))
        elif key not in payload:
            payload[key] = []
    metadata_path = artifacts / "run_metadata.json"
    if metadata_path.is_file():
        payload["metadata"] = json.loads(metadata_path.read_text(encoding="utf-8"))
    else:
        payload.setdefault("metadata", {"status": "UNKNOWN"})
    return payload, artifacts


def _render_results(
    payload: dict, selected_run: dict[str, str], artifact_dir: Path | None
) -> None:
    scores = pd.DataFrame(payload.get("scores", []))
    interactions = pd.DataFrame(payload.get("interactions", []))
    target_summary = pd.DataFrame(payload.get("target_preparations", []))

    st.caption(
        f"Run {selected_run['run_id'][:8]} · "
        f"{selected_run['created_at'][:19].replace('T', ' ')} UTC"
    )
    run_status = payload.get("metadata", {}).get("status", "UNKNOWN")
    if run_status == "COMPLETED":
        st.success(f"Run status: {run_status}")
    elif run_status == "COMPLETED_WITH_ERRORS":
        st.warning(f"Run status: {run_status}. Review pair-level failures below.")
    elif run_status == "LEGACY":
        st.info("Legacy run: its saved summaries may be available, but its artifact folder is not.")
    else:
        st.error(f"Run status: {run_status}")
    if "status" in target_summary:
        targets_count = int(
            target_summary["status"].isin({"OK", "BLIND_WHOLE_PROTEIN"}).sum()
        )
    else:
        targets_count = len(target_summary)
    ligand_columns = [
        column
        for column in scores.columns
        if column not in {"target", "pdb_id", "uniprot_id", "chain", "box_basis"}
    ]
    docked_count = int(scores[ligand_columns].notna().sum().sum()) if ligand_columns else 0
    c1, c2, c3 = st.columns(3)
    c1.metric("Prepared targets", targets_count)
    c2.metric("Ligands", len(ligand_columns))
    c3.metric("Successful docking scores", docked_count)
    st.info(
        "Vina scores are approximate and target-specific. Compare compounds primarily "
        "within the same target; do not treat scores across different proteins as "
        "directly comparable binding affinities."
    )

    score_tab, interaction_tab, qc_tab, files_tab = st.tabs(
        ["Score landscape", "Interaction analysis", "Target preparation", "Downloads"]
    )
    with score_tab:
        if scores.empty or not ligand_columns:
            st.warning("No successful docking scores are available for this run.")
        else:
            st.subheader("Docking score matrix")
            st.caption("Hover over a tile for its target, ligand, and score.")
            heatmap = px.imshow(
                scores.set_index("target")[ligand_columns],
                labels={"x": "Ligand", "y": "Target", "color": "Vina score (kcal/mol)"},
                color_continuous_scale="RdYlGn_r",
                aspect="auto",
                text_auto=".2f",
            )
            heatmap.update_layout(height=max(380, 40 * len(scores)), margin=dict(l=20, r=20, t=30, b=20))
            st.plotly_chart(heatmap, width="stretch")
            st.subheader("Rank ligands within one target")
            target_names = scores["target"].astype(str).tolist()
            selected_target = st.selectbox("Target", target_names, key="rank_target")
            target_scores = scores.loc[scores["target"].astype(str) == selected_target, ligand_columns]
            ranked = (
                target_scores.iloc[0]
                .rename_axis("Ligand")
                .reset_index(name="Vina score (kcal/mol)")
                .dropna()
                .sort_values("Vina score (kcal/mol)")
            )
            if not ranked.empty:
                chart = px.bar(
                    ranked,
                    x="Vina score (kcal/mol)",
                    y="Ligand",
                    orientation="h",
                    color="Vina score (kcal/mol)",
                    color_continuous_scale="RdYlGn_r",
                    title=f"Within-target score ranking: {selected_target}",
                )
                chart.update_layout(yaxis={"categoryorder": "total ascending"}, height=max(320, 36 * len(ranked)))
                st.plotly_chart(chart, width="stretch")
            st.dataframe(scores, width="stretch", hide_index=True)

            if len(ligand_columns) > 0:
                st.subheader("Ligand structure")
                selected_ligand = st.selectbox("Ligand", ligand_columns, key="structure_ligand")
                smiles = next(
                    (item["smiles"] for item in payload.get("ligands", []) if item["name"] == selected_ligand),
                    None,
                )
                molecule = Chem.MolFromSmiles(smiles) if smiles else None
                if molecule is not None:
                    st.image(Draw.MolToImage(molecule, size=(520, 360)), caption=selected_ligand)
                    st.code(smiles, language="text")

    with interaction_tab:
        st.subheader("ProLIF interaction fingerprint")
        st.caption(
            "These are geometry-based interaction classifications for docked poses; "
            "they are hypotheses for follow-up, not experimental evidence."
        )
        if interactions.empty:
            st.info("No interaction records are available in this saved run.")
        else:
            pair_columns = st.columns(2)
            interaction_targets = ["All targets"] + sorted(interactions["target"].unique().tolist())
            interaction_ligands = ["All ligands"] + sorted(interactions["ligand"].unique().tolist())
            selected_interaction_target = pair_columns[0].selectbox(
                "Filter target", interaction_targets, key="interaction_target"
            )
            selected_interaction_ligand = pair_columns[1].selectbox(
                "Filter ligand", interaction_ligands, key="interaction_ligand"
            )
            filtered = interactions
            if selected_interaction_target != "All targets":
                filtered = filtered[filtered["target"] == selected_interaction_target]
            if selected_interaction_ligand != "All ligands":
                filtered = filtered[filtered["ligand"] == selected_interaction_ligand]
            if filtered.empty:
                st.info("No detected interactions match these filters.")
            else:
                counts = (
                    filtered.groupby("interaction", as_index=False)
                    .size()
                    .sort_values("size", ascending=False)
                )
                st.plotly_chart(
                    px.bar(
                        counts,
                        x="interaction",
                        y="size",
                        color="interaction",
                        labels={"size": "Detected residue contacts"},
                        title="Detected interaction types",
                    ),
                    width="stretch",
                )
                st.dataframe(filtered, width="stretch", hide_index=True)
        errors = pd.DataFrame(payload.get("interaction_errors", []))
        if not errors.empty:
            with st.expander(f"Interaction analysis issues ({len(errors)})"):
                st.dataframe(errors, width="stretch", hide_index=True)
        pair_results = pd.DataFrame(payload.get("pair_results", []))
        pair_failures = (
            pair_results.loc[pair_results["status"] != "SUCCEEDED"]
            if not pair_results.empty and "status" in pair_results.columns
            else pd.DataFrame()
        )
        if not pair_failures.empty:
            with st.expander(
                f"Failed or skipped ligand-target pairs ({len(pair_failures)})"
            ):
                st.dataframe(pair_failures, width="stretch", hide_index=True)

    with qc_tab:
        st.subheader("Target preparation quality control")
        st.write(
            "Preparation status and pocket derivation method help interpret which "
            "docking boxes were used. Whole-protein blind docking is less targeted "
            "than a box derived from a bound ligand or supplied pocket."
        )
        st.dataframe(target_summary, width="stretch", hide_index=True)

    with files_tab:
        st.subheader("Download run outputs")
        score_frame = scores.copy()
        st.download_button(
            "Download docking scores (.xlsx)",
            _xlsx_bytes(score_frame, "Docking scores"),
            file_name="docking_scores.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
        st.download_button(
            "Download interaction report (.csv)",
            interactions.to_csv(index=False).encode("utf-8"),
            file_name="interactions.csv",
            mime="text/csv",
        )
        summary_bytes = _xlsx_bytes(target_summary, "Target preparation")
        st.download_button(
            "Download target preparation summary (.xlsx)",
            summary_bytes,
            file_name="targets_summary.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
        if artifact_dir is not None and artifact_dir.is_dir():
            st.download_button(
                "Download structures, poses, logs and configs (.zip)",
                _zip_job(artifact_dir),
                file_name="dockex_run_results.zip",
                mime="application/zip",
            )
            metadata = payload.get("metadata", {})
            with st.expander("Run provenance and settings"):
                st.json(metadata)
        else:
            st.warning(
                "Run artifacts are unavailable. The local history entry remains, but "
                "the original workbook, structures, poses, and logs cannot be downloaded."
            )


def _render_assistant(payload: dict, selected_run: dict[str, str], db_path: Path) -> None:
    st.subheader("Ask about this docking run")
    st.write(
        "The assistant interprets the selected run’s scores, target preparation "
        "summary, and ProLIF interactions. It does not dock molecules or replace "
        "scientific review."
    )
    api_key = st.text_input(
        "Gemini API key",
        type="password",
        key="gemini_api_key",
        help="Used only in this active app session; it is not saved in the local database.",
    )
    st.caption(f"Reasoning model: `{MODEL_ID}`")
    consent = st.checkbox(
        "I understand that my question, chat history, ligand names/SMILES, target identifiers, "
        "scores, and detected interactions are sent to Google Gemini.",
        key="gemini_data_consent",
    )
    run_id = selected_run["run_id"]
    for message in get_chat_history(db_path, run_id):
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

    question = st.chat_input(
        "Ask about score patterns, a residue contact, or analysis limitations",
        disabled=not (api_key and consent),
    )
    if question:
        if not consent:
            st.error("Confirm the data-sharing notice before sending a question.")
            return
        history = get_chat_history(db_path, run_id)
        try:
            with st.spinner("Reviewing this run with Gemini…"):
                answer = answer_question(api_key, payload, history, question)
            add_chat_turn(db_path, run_id, question, answer)
            st.rerun()
        except Exception as exc:
            st.error(f"Gemini request failed: {exc}")


DB_PATH = database_path()
initialize_database(DB_PATH)

st.set_page_config(page_title="DockEx", page_icon="🧪", layout="wide")
st.markdown(
    """
    <style>
    .block-container { max-width: 1500px; padding-top: 1.8rem; }
    [data-testid="stMetric"] {
      background: linear-gradient(135deg, rgba(20, 85, 100, .10), rgba(38, 166, 154, .08));
      border: 1px solid rgba(38, 166, 154, .25);
      border-radius: 12px;
      padding: 14px 18px;
    }
    div[data-testid="stFileUploader"] { border-radius: 12px; }
    </style>
    """,
    unsafe_allow_html=True,
)
st.title("DockEx · Molecular Docking & Virtual Screening")
st.markdown(
    "Prepare protein structures, dock ligand libraries, inspect interaction fingerprints, "
    "and review the evidence and provenance for each ligand-target pair."
)

with st.expander("How to use this workbench", expanded=True):
    instructions = st.columns(4)
    instructions[0].markdown("**1 · Prepare input**  \nUse the DockEx example workbook or your own two-sheet `.xlsx`.")
    instructions[1].markdown("**2 · Validate**  \nCheck ligand SMILES, target identifiers, and the estimated L×T workload.")
    instructions[2].markdown("**3 · Run & inspect**  \nTrack each preparation and docking step, then inspect scores and ProLIF contacts.")
    instructions[3].markdown("**4 · Ask**  \nOptionally enter a Gemini key to discuss this run; the key is session-only.")
    st.caption(
        "Target structures are fetched from RCSB PDB or AlphaFold over the internet. "
        "Each run and its outputs are retained in the local DockEx data directory. "
        "A larger ligand × target matrix increases runtime."
    )
st.subheader("New Screening")
if EXAMPLE_WORKBOOK.is_file():
    st.download_button(
        "Download EGFR inhibitor screening example (.xlsx)",
        EXAMPLE_WORKBOOK.read_bytes(),
        file_name="dockex_sample_input.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

uploaded_workbook = st.file_uploader("Upload your Excel workbook (.xlsx)", type=["xlsx"])
parsed = None
if uploaded_workbook is not None:
    try:
        parsed = parse_workbook(uploaded_workbook.getvalue())
    except WorkbookValidationError as exc:
        st.error(str(exc))
    except Exception as exc:
        st.error(f"Could not validate workbook: {exc}")

if parsed is not None:
    matrix_size = len(parsed.ligands) * len(parsed.targets)
    col1, col2, col3 = st.columns(3)
    col1.metric("Ligands", len(parsed.ligands))
    col2.metric("Targets", len(parsed.targets))
    col3.metric("Planned ligand × target pairs", matrix_size)
    with st.expander("Review uploaded molecules and targets", expanded=True):
        ligand_tab, target_tab = st.tabs(["Ligands", "Targets"])
        with ligand_tab:
            st.dataframe(
                [{"name": ligand.name, "smiles": ligand.smiles} for ligand in parsed.ligands],
                width="stretch",
                hide_index=True,
            )
        with target_tab:
            st.dataframe(
                [
                    {
                        "target": target.name,
                        "pdb_id": target.pdb_id or "",
                        "uniprot_id": target.uniprot_id or "",
                        "chain": target.chain or "",
                        "manual pocket residues": ", ".join(map(str, target.pocket_residues)),
                    }
                    for target in parsed.targets
                ],
                width="stretch",
                hide_index=True,
            )

    vina = _find_vina()
    if vina is None:
        st.error(
            "AutoDock Vina was not found. Install Vina and set VINA_EXECUTABLE to its path, "
            "or place its executable in this app directory."
        )
    else:
        st.caption(f"Docking engine: AutoDock Vina · {vina.name}")
        with st.form("run_settings"):
            settings1, settings2 = st.columns(2)
            exhaustiveness = settings1.number_input(
                "Search exhaustiveness", min_value=1, max_value=128, value=8, step=1,
                help="Higher values spend more search effort per pose and take longer.",
            )
            seed = settings2.number_input(
                "Random seed (0 uses Vina default)", min_value=0, value=0, step=1
            )
            run_clicked = st.form_submit_button("Prepare structures and run docking", type="primary")

        if run_clicked:
            run_id = str(uuid.uuid4())
            job_dir = DB_PATH.parent / "runs" / run_id
            job_dir.mkdir(parents=True, exist_ok=False)
            input_dir = job_dir / "input"
            input_dir.mkdir(parents=True, exist_ok=True)
            workbook_bytes = uploaded_workbook.getvalue()
            (input_dir / "workbook.xlsx").write_bytes(workbook_bytes)
            label = f"{len(parsed.ligands)} ligands × {len(parsed.targets)} targets"
            progress_bar = st.progress(0.0, text="Starting preparation")
            status = st.empty()

            def update_progress(completed: int, total: int, message: str) -> None:
                progress_bar.progress(completed / total, text=message)
                status.write(f"Step {completed} of {total} · {message}")

            try:
                results_path, summary_path = run_docking_job(
                    parsed.ligands,
                    parsed.targets,
                    job_dir,
                    vina,
                    int(exhaustiveness),
                    int(seed) if seed else None,
                    update_progress,
                    input_workbook=workbook_bytes,
                    input_filename=Path(uploaded_workbook.name).name,
                )
                metadata = json.loads(
                    (job_dir / "run_metadata.json").read_text(encoding="utf-8")
                )
                run_payload = {
                    "metadata": metadata,
                    "ligands": [{"name": item.name, "smiles": item.smiles} for item in parsed.ligands],
                    "targets": [
                        {
                            "name": item.name,
                            "pdb_id": item.pdb_id,
                            "uniprot_id": item.uniprot_id,
                        }
                        for item in parsed.targets
                    ],
                    "artifact_dir": (Path("runs") / run_id).as_posix(),
                }
                save_run(DB_PATH, run_id, label, run_payload)
                st.session_state["current_run_id"] = run_id
                if metadata["status"] == "COMPLETED":
                    st.success("DockEx run completed successfully.")
                else:
                    st.warning("DockEx run completed with pair-level issues. Review its reports.")
                st.rerun()
            except Exception as exc:
                st.error(f"Docking run failed: {exc}")
                failure_payload = {
                    "metadata": {
                        "application": "DockEx",
                        "status": "FAILED",
                        "error": str(exc),
                    },
                    "ligands": [{"name": item.name, "smiles": item.smiles} for item in parsed.ligands],
                    "targets": [
                        {
                            "name": item.name,
                            "pdb_id": item.pdb_id,
                            "uniprot_id": item.uniprot_id,
                        }
                        for item in parsed.targets
                    ],
                    "artifact_dir": (Path("runs") / run_id).as_posix(),
                }
                (job_dir / "run_error.txt").write_text(str(exc), encoding="utf-8")
                try:
                    save_run(DB_PATH, run_id, label, failure_payload)
                    st.session_state["current_run_id"] = run_id
                except Exception as save_exc:
                    st.error(f"Could not save failed-run history: {save_exc}")
                log_file = job_dir / "application.log"
                if log_file.exists():
                    st.download_button(
                        "Download run log",
                        log_file.read_bytes(),
                        file_name="application.log",
                        mime="text/plain",
                    )
                summary_file = job_dir / "targets_summary.xlsx"
                if summary_file.exists():
                    st.download_button(
                        "Download target preparation summary",
                        summary_file.read_bytes(),
                        file_name="targets_summary.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    )
                st.download_button(
                    "Download partial run outputs (ZIP)",
                    _zip_job(job_dir),
                    file_name="dockex_partial_run_results.zip",
                    mime="application/zip",
                )

runs = list_runs(DB_PATH)
if runs:
    current_run_id = st.session_state.get("current_run_id")
    run_ids = [run["run_id"] for run in runs]
    default_index = run_ids.index(current_run_id) if current_run_id in run_ids else 0
    selected_run_id = st.selectbox(
        "DockEx run history",
        run_ids,
        index=default_index,
        format_func=lambda value: next(
            f"{item['label']} · {item['created_at'][:19].replace('T', ' ')} UTC"
            for item in runs
            if item["run_id"] == value
        ),
    )
    selected_run = next(item for item in runs if item["run_id"] == selected_run_id)
    delete_confirmed = st.checkbox(
        "Permanently delete this run, its local artifacts, and its chat history",
        key=f"delete_run_confirm_{selected_run_id}",
    )
    if st.button("Delete selected run", disabled=not delete_confirmed):
        try:
            if delete_run(DB_PATH, selected_run_id):
                if st.session_state.get("current_run_id") == selected_run_id:
                    del st.session_state["current_run_id"]
                st.success("Run and retained artifacts deleted.")
                st.rerun()
            else:
                st.error("The selected run no longer exists in local history.")
        except OSError as exc:
            st.error(f"Could not delete the run artifacts: {exc}")

    try:
        payload, saved_artifacts = _load_saved_payload(DB_PATH, selected_run_id)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        st.error(f"Could not load saved run reports: {exc}")
        payload, saved_artifacts = None, None
    if payload is not None:
        results_tab, assistant_tab = st.tabs(["Results & visualizations", "Scientific Q&A"])
        with results_tab:
            _render_results(payload, selected_run, saved_artifacts)
        with assistant_tab:
            _render_assistant(payload, selected_run, DB_PATH)
