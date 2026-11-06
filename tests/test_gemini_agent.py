import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import gemini_agent
from google import genai


def test_prompt_includes_run_evidence_and_recent_chat():
    prompt = gemini_agent.build_agent_prompt(
        {
            "scores": [{"target": "target-a", "compound": -8.1}],
            "artifact_dir": "C:\\private\\temporary\\path",
        },
        [{"role": "user", "content": "What is the score?"}],
        "Which contacts were detected?",
    )

    assert "target-a" in prompt
    assert "What is the score?" in prompt
    assert "Which contacts were detected?" in prompt
    assert "private" not in prompt


def test_answer_question_uses_selected_supported_model_and_thinking(monkeypatch):
    calls = {}

    class FakeModels:
        def generate_content(self, **kwargs):
            calls.update(kwargs)

            class Response:
                text = "The score is an approximate ranking signal."

            return Response()

    class FakeClient:
        def __init__(self, api_key):
            calls["api_key"] = api_key
            self.models = FakeModels()

    monkeypatch.setattr(genai, "Client", FakeClient)
    answer = gemini_agent.answer_question("test-key", {"scores": []}, [], "Explain scores")

    assert answer.startswith("The score")
    assert calls["api_key"] == "test-key"
    assert calls["model"] == "gemini-3.5-flash-lite"
    assert calls["config"].thinking_config.thinking_level == "MEDIUM"
    assert calls["contents"].endswith("Explain scores")


def test_answer_question_requires_api_key():
    with pytest.raises(ValueError, match="API key"):
        gemini_agent.answer_question("", {}, [], "Question")
