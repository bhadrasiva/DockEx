"""Local SQLite history and safe access to persistent run artifacts."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sqlite3
from typing import Any


def database_path() -> Path:
    configured = os.environ.get("DOCKEX_DATA_DIR") or os.environ.get(
        "MOLECULAR_DOCKING_DATA_DIR"
    )
    data_dir = Path(configured).expanduser() if configured else Path.home() / ".molecular_docking"
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir / "history.sqlite3"


def artifact_directory(
    path: Path, payload: dict[str, Any], run_id: str | None = None
) -> Path | None:
    """Resolve a saved run directory only when it is inside this database's runs folder."""
    stored_path = payload.get("artifact_dir")
    if not isinstance(stored_path, str) or not stored_path:
        return None
    candidate = Path(stored_path)
    if not candidate.is_absolute():
        candidate = path.parent / candidate
    if candidate.is_symlink():
        return None
    root = (path.parent / "runs").resolve()
    resolved = candidate.resolve()
    if resolved == root or resolved.parent != root:
        return None
    if run_id is not None and resolved.name != run_id:
        return None
    try:
        resolved.relative_to(root)
    except ValueError:
        return None
    return resolved


def _connect(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path, timeout=10)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def initialize_database(path: Path) -> None:
    with _connect(path) as connection:
        connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS runs (
                run_id TEXT PRIMARY KEY,
                label TEXT NOT NULL,
                created_at TEXT NOT NULL,
                payload_json TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS chat_messages (
                message_id INTEGER PRIMARY KEY AUTOINCREMENT,
                run_id TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
                role TEXT NOT NULL CHECK (role IN ('user', 'assistant')),
                content TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS idx_chat_messages_run
                ON chat_messages(run_id, message_id);
            """
        )


def save_run(path: Path, run_id: str, label: str, payload: dict[str, Any]) -> None:
    created_at = datetime.now(timezone.utc).isoformat()
    payload_json = json.dumps(payload, ensure_ascii=False, allow_nan=False)
    with _connect(path) as connection:
        connection.execute(
            """
            INSERT INTO runs(run_id, label, created_at, payload_json)
            VALUES (?, ?, ?, ?)
            ON CONFLICT(run_id) DO UPDATE SET
                label = excluded.label,
                payload_json = excluded.payload_json
            """,
            (run_id, label, created_at, payload_json),
        )


def list_runs(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    with _connect(path) as connection:
        query = (
            "SELECT run_id, label, created_at "
            "FROM runs ORDER BY created_at DESC"
        )
        if limit is None:
            rows = connection.execute(query).fetchall()
        else:
            rows = connection.execute(query + " LIMIT ?", (limit,)).fetchall()
    return [dict(row) for row in rows]


def load_run(path: Path, run_id: str) -> dict[str, Any] | None:
    with _connect(path) as connection:
        row = connection.execute(
            "SELECT payload_json FROM runs WHERE run_id = ?", (run_id,)
        ).fetchone()
    return json.loads(row["payload_json"]) if row else None


def delete_run(path: Path, run_id: str) -> bool:
    payload = load_run(path, run_id)
    if payload is None:
        return False
    artifacts = artifact_directory(path, payload, run_id)
    if artifacts is not None and artifacts.is_dir():
        shutil.rmtree(artifacts)
    with _connect(path) as connection:
        cursor = connection.execute("DELETE FROM runs WHERE run_id = ?", (run_id,))
    return cursor.rowcount == 1


def add_chat_turn(path: Path, run_id: str, question: str, answer: str) -> None:
    created_at = datetime.now(timezone.utc).isoformat()
    with _connect(path) as connection:
        connection.executemany(
            """
            INSERT INTO chat_messages(run_id, role, content, created_at)
            VALUES (?, ?, ?, ?)
            """,
            (
                (run_id, "user", question, created_at),
                (run_id, "assistant", answer, created_at),
            ),
        )


def get_chat_history(path: Path, run_id: str, limit: int = 40) -> list[dict[str, str]]:
    with _connect(path) as connection:
        rows = connection.execute(
            """
            SELECT role, content, created_at FROM (
                SELECT message_id, role, content, created_at
                FROM chat_messages
                WHERE run_id = ?
                ORDER BY message_id DESC
                LIMIT ?
            )
            ORDER BY message_id
            """,
            (run_id, limit),
        ).fetchall()
    return [dict(row) for row in rows]
