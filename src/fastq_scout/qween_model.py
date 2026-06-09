from __future__ import annotations

from abc import ABC, abstractmethod

from fastq_scout.explain import build_explain_payload, payload_to_prompt_text

SYSTEM_PROMPT = """You are the FastqScout assistant at the Bioinfomics laboratory.

Your audience: wet-lab biologists who are not bioinformaticians. They need a clear,
reassuring, practical explanation of a pre-flight FASTQ QC report.

Rules (strict):
1. Do NOT change or question the verdict field in the JSON — treat it as final.
2. Quote numbers ONLY from the JSON. Never invent metrics, percentages, or tool names.
3. Do not claim the full file was analyzed if sampling shows only a fraction was read.
4. Use plain language; briefly explain jargon (PHRED, Q20, adapter, duplicate rate).
5. Structure your answer with these headings:
   ## Краткий вывод
   ## Что выглядит хорошо
   ## На что обратить внимание
   ## Что делать дальше
6. Keep the answer under 350 words. Write in Russian.
7. If data is missing, say so — do not guess.
8. This is QC guidance, not a clinical or diagnostic conclusion."""


USER_PROMPT_TEMPLATE = """Explain the following FastqScout report to a wet-lab scientist.

The verdict and all numbers below are authoritative. Your job is to translate them
into plain language and practical next steps.

```json
{payload_json}
```"""


class BaseModel(ABC):
    @abstractmethod
    def name(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def generate(self, payload: dict) -> str:
        raise NotImplementedError


class QwenModel(BaseModel):
    def __init__(self, model_name: str = "Qwen/Qwen2.5-1.5B-Instruct"):
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
            max_new_tokens=512,
            do_sample=False,
        )
        new_tokens = generated_ids[0][input_len:]
        response = self._tokenizer.decode(new_tokens, skip_special_tokens=True).strip()
        return response
