import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chat_store import (
    add_chat_turn,
    artifact_directory,
    delete_run,
    get_chat_history,
    initialize_database,
    list_runs,
    load_run,
    save_run,
)


def test_runs_and_chat_turns_persist_in_sqlite(tmp_path):
    database = tmp_path / "history.sqlite3"
    initialize_database(database)
    payload = {
        "scores": [{"target": "target-a", "ligand-a": -7.2}],
        "interactions": [{"target": "target-a", "ligand": "ligand-a", "interaction": "HBDonor"}],
    }

    save_run(database, "run-1", "1 ligand × 1 target", payload)
    add_chat_turn(database, "run-1", "What does this contact mean?", "It is a predicted donor contact.")

    assert load_run(database, "run-1") == payload
    assert list_runs(database)[0]["run_id"] == "run-1"
    history = get_chat_history(database, "run-1")
    assert [message["role"] for message in history] == ["user", "assistant"]
    assert history[0]["content"] == "What does this contact mean?"


def test_run_payload_rejects_non_json_nan_values(tmp_path):
    database = tmp_path / "history.sqlite3"
    initialize_database(database)

    try:
        save_run(database, "run-1", "bad payload", {"score": float("nan")})
    except ValueError as exc:
        assert "JSON compliant" in str(exc)
    else:
        raise AssertionError("NaN should not be stored in JSON")


def test_run_artifacts_are_resolved_and_deleted_with_run_history(tmp_path):
    database = tmp_path / "history.sqlite3"
    initialize_database(database)
    run_directory = tmp_path / "runs" / "run-2"
    run_directory.mkdir(parents=True)
    (run_directory / "pose.pdbqt").write_text("pose")
    payload = {"artifact_dir": "runs/run-2", "metadata": {"status": "COMPLETED"}}
    save_run(database, "run-2", "test run", payload)
    add_chat_turn(database, "run-2", "question", "answer")

    assert artifact_directory(database, payload) == run_directory.resolve()
    assert load_run(database, "run-2")["metadata"]["status"] == "COMPLETED"
    assert list_runs(database)[0]["run_id"] == "run-2"
    assert delete_run(database, "run-2")
    assert not run_directory.exists()
    assert load_run(database, "run-2") is None
    assert get_chat_history(database, "run-2") == []


def test_artifact_resolution_rejects_paths_outside_run_storage(tmp_path):
    database = tmp_path / "history.sqlite3"
    initialize_database(database)
    outside = tmp_path / "private"
    outside.mkdir()

    assert artifact_directory(database, {"artifact_dir": str(outside)}) is None
    assert artifact_directory(database, {"artifact_dir": "runs"}) is None
    assert (
        artifact_directory(database, {"artifact_dir": "runs/other-run"}, "run-2")
        is None
    )
