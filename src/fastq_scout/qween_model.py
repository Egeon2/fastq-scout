from __future__ import annotations

import re
from abc import ABC, abstractmethod

from fastq_scout.explain import build_explain_payload, payload_to_prompt_text

SYSTEM_PROMPT = (
    "You are the FastqScout assistant for wet-lab scientists. "
    "Write only in English. Use only facts from the JSON. Do not invent numbers. "
    "Do not change the verdict. Keep plain language. Use exactly four ## headings."
)

USER_PROMPT_TEMPLATE = """Explain this FastqScout pre-flight QC report in simple English.

Use ONLY numbers and facts from the JSON below.

Format (exactly 4 sections):

## Summary
1–2 sentences: verdict + main issue or that the sample looks fine.

## What looks good
2–4 bullet points starting with "- ".

## What to watch
Bullet points from issues and metrics (adapter %, duplicates, quality drop).

## Next steps
Bullet points from recommendations.

JSON:
{payload_json}"""

_BAD_PATTERNS = re.compile(
    r"grammar fragment|DNA fragment member|chromosome member|"
    r"геном|грамматик|фрагмент|член",
    re.IGNORECASE,
)

_REQUIRED_HEADINGS = (
    "## Summary",
    "## Next steps",
)


def is_llm_response_usable(text: str, payload: dict) -> bool:
    if not text or len(text) < 80:
        return False
    if _BAD_PATTERNS.search(text):
        return False

    # Reject Cyrillic — we want English output only.
    if re.search(r"[\u0400-\u04FF]", text):
        return False

    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if len(lines) != len(set(lines)) and len(lines) > 5:
        from collections import Counter

        counts = Counter(lines)
        if any(c >= 3 for c in counts.values()):
            return False

    if not all(h in text for h in _REQUIRED_HEADINGS):
        return False

    verdict = payload.get("verdict", "")
    if verdict and verdict not in text.upper():
        return False

    return True


class BaseModel(ABC):
    @abstractmethod
    def name(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def generate(self, payload: dict) -> str:
        raise NotImplementedError


class QwenModel(BaseModel):
    def __init__(self, model_name: str = "Qwen/Qwen2.5-0.5B-Instruct"):
        self.model_name = model_name
        self._model = None
        self._tokenizer = None

    def _load(self) -> None:
        if self._model is not None:
            return

        from transformers import AutoModelForCausalLM, AutoTokenizer

        print(f"Loading {self.model_name}...")
        self._model = AutoModelForCausalLM.from_pretrained(
            self.model_name,
            torch_dtype="auto",
            device_map="auto",
        )
        self._tokenizer = AutoTokenizer.from_pretrained(self.model_name)

    def name(self) -> str:
        return self.model_name

    def generate(self, payload: dict) -> str:
        self._load()

        if "reads" not in payload and "results" in payload:
            payload = build_explain_payload(
                metrics=payload.get("results") or payload.get("results_r1"),
                scout_report=payload["scout_report"],
                sample_plan=payload.get("sample_plan", {}),
                fastq_path=payload["fastq_path"],
                r2_metrics=payload.get("r2_metrics"),
            )

        payload_json = payload_to_prompt_text(payload)
        user_content = USER_PROMPT_TEMPLATE.format(payload_json=payload_json)

        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_content},
        ]
        text = self._tokenizer.apply_chat_template(
            messages,
            tokenize=False,
            add_generation_prompt=True,
        )
        model_inputs = self._tokenizer([text], return_tensors="pt").to(self._model.device)
        input_len = model_inputs["input_ids"].shape[1]

        generated_ids = self._model.generate(
            **model_inputs,
            max_new_tokens=400,
            do_sample=False,
            repetition_penalty=1.15,
        )
        new_tokens = generated_ids[0][input_len:]
        return self._tokenizer.decode(new_tokens, skip_special_tokens=True).strip()
