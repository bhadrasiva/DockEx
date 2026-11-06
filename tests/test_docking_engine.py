from pathlib import Path
import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "engine"))

import dock_targets


def _dock(tmp_path, monkeypatch, process_result):
    config_dir = tmp_path / "config"
    log_dir = tmp_path / "logs"
    pose_dir = tmp_path / "poses"
    for directory in (config_dir, log_dir, pose_dir):
        directory.mkdir()

    monkeypatch.setattr(dock_targets.subprocess, "run", lambda *args, **kwargs: process_result)
    return dock_targets.dock(
        Path("vina.exe"),
        Path("receptor.pdbqt"),
        Path("ligand.pdbqt"),
        (0.0, 0.0, 0.0),
        (20.0, 20.0, 20.0),
        "target_ligand",
        8,
        None,
        config_dir,
        log_dir,
        pose_dir,
    )


def test_dock_reports_vina_process_failure(tmp_path, monkeypatch):
    result = _dock(
        tmp_path,
        monkeypatch,
        SimpleNamespace(returncode=1, stdout="", stderr="invalid receptor"),
    )

    assert result.docking_score is None
    assert result.pose_path is None
    assert result.error == "invalid receptor"
    assert "invalid receptor" in (
        tmp_path / "logs" / "target_ligand.log"
    ).read_text()


def test_dock_rejects_score_without_pose_file(tmp_path, monkeypatch):
    result = _dock(
        tmp_path,
        monkeypatch,
        SimpleNamespace(returncode=0, stdout="  1 -7.50 0.0 0.0\n", stderr=""),
    )

    assert result.docking_score is None
    assert result.pose_path is None
    assert "did not create the pose file" in result.error


def test_dock_returns_score_and_pose_when_both_exist(tmp_path, monkeypatch):
    def run_vina(command, **kwargs):
        pose_path = Path(command[command.index("--out") + 1])
        pose_path.write_text("pose")
        return SimpleNamespace(
            returncode=0,
            stdout="  1 -7.50 0.0 0.0\n",
            stderr="",
        )

    monkeypatch.setattr(dock_targets.subprocess, "run", run_vina)
    config_dir = tmp_path / "config"
    log_dir = tmp_path / "logs"
    pose_dir = tmp_path / "poses"
    for directory in (config_dir, log_dir, pose_dir):
        directory.mkdir()

    result = dock_targets.dock(
        Path("vina.exe"),
        Path("receptor.pdbqt"),
        Path("ligand.pdbqt"),
        (0.0, 0.0, 0.0),
        (20.0, 20.0, 20.0),
        "target_ligand",
        8,
        17,
        config_dir,
        log_dir,
        pose_dir,
    )

    assert result.docking_score == -7.5
    assert Path(result.pose_path).is_file()
    assert result.error is None
