"""Context-grounded Gemini assistant for interpreting a completed docking run."""

from __future__ import annotations

import json
from typing import Any

MODEL_ID = "gemini-3.5-flash-lite"
MAX_CONTEXT_CHARS = 60_000
MAX_HISTORY_MESSAGES = 20

SYSTEM_INSTRUCTION = """You are a scientific assistant for interpreting one molecular-docking run.
Use only the supplied run data for claims about its scores and detected interactions.
Explain that Vina scores are approximate scoring-function outputs, not binding
affinities, and that ProLIF contacts are computational classifications, not
experimental evidence. Compare ligand scores only within the same target unless
the user explicitly asks for a cautious cross-target discussion. Distinguish
observations in the uploaded run from general scientific background. Do not
invent missing contacts, experimental results, or confidence. Explain methods
and limitations clearly; do not reveal hidden chain-of-thought. This is not
medical advice or a substitute for experimental validation."""


def build_agent_prompt(
    run_payload: dict[str, Any],
    history: list[dict[str, str]],
    question: str,
) -> str:
    evidence_keys = (
        "ligands",
        "targets",
        "scores",
        "interactions",
        "interaction_errors",
        "target_preparations",
    )
    evidence = {key: run_payload.get(key) for key in evidence_keys if key in run_payload}
    context = json.dumps(evidence, ensure_ascii=False, allow_nan=False)
    if len(context) > MAX_CONTEXT_CHARS:
        context = context[:MAX_CONTEXT_CHARS] + "\n[Run context truncated to fit the assistant request.]"

    conversation = history[-MAX_HISTORY_MESSAGES:]
    transcript = "\n".join(
        f"{message['role'].upper()}: {message['content']}" for message in conversation
    )
    return (
        f"DOCKING RUN DATA (JSON):\n{context}\n\n"
        f"RECENT CONVERSATION:\n{transcript or '[No earlier messages]'}\n\n"
        f"USER QUESTION:\n{question.strip()}"
    )


def answer_question(
    api_key: str,
    run_payload: dict[str, Any],
    history: list[dict[str, str]],
    question: str,
) -> str:
    if not api_key.strip():
        raise ValueError("Enter a Gemini API key to use the scientific assistant.")
    if not question.strip():
        raise ValueError("Enter a question about the docking run.")

    from google import genai
    from google.genai import types

    client = genai.Client(api_key=api_key.strip())
    response = client.models.generate_content(
        model=MODEL_ID,
        contents=build_agent_prompt(run_payload, history, question),
        config=types.GenerateContentConfig(
            system_instruction=SYSTEM_INSTRUCTION,
            temperature=0.2,
            thinking_config=types.ThinkingConfig(thinking_level="medium"),
        ),
    )
    answer = response.text
    if not answer or not answer.strip():
        raise RuntimeError("Gemini returned an empty response.")
    return answer.strip()
