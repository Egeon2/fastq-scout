from __future__ import annotations

from abc import ABC, abstractmethod

from fastq_scout.explain import build_explain_payload, extract_key_facts

SYSTEM_PROMPT = (
    "You rewrite FASTQ QC summaries for wet-lab scientists. "
    "Write ONLY in English. Use simple words. "
    "Never change the verdict or any number. Never invent facts."
)

# Short few-shot pairs: draft → polished (teaches format + tone for small models).
FEW_SHOT_DRAFT = """## Summary
Verdict: TRIM. Preprocessing needed before alignment.
Analyzed 10,000 reads (15% sample).

## What looks good
- R1: good average read quality (PHRED 34)
- R1: 100% of bases at Q≥20

## What to watch
- Adapter sequence detected on read tails (13.65%)
- Quality drops at read tail (head 36.5 vs tail 31.0 PHRED)

## Next steps
- Run fastp --adapter_sequence GTCTGAACTCCAGTCAC"""

FEW_SHOT_REWRITE = """## Summary
Verdict: TRIM. This sample should be trimmed before you run alignment — adapter sequences were found on many read tails.
FastqScout checked 10,000 reads (15% of the file).

## What looks good
- Overall read quality is strong (mean PHRED 34).
- Essentially all bases pass Q≥20 (100%).

## What to watch
- Illumina adapter signal on 13.65% of read tails — trim these before mapping.
- Quality falls toward the end of reads (36.5 at the start vs 31.0 at the tail).

## Next steps
- Run fastp with: --adapter_sequence GTCTGAACTCCAGTCAC
- Re-run QC after trimming."""

USER_PROMPT_TEMPLATE = """Rewrite the draft summary below in friendly plain English.

Rules:
- English only
- Verdict must stay: {verdict}
- Keep every number exactly (PHRED scores, percentages, read counts)
- Keep exactly 4 sections with these headings:
  ## Summary
  ## What looks good
  ## What to watch
  ## Next steps
- Use "- " bullet points where appropriate
- Do not add facts that are not in the draft

Example draft:
{few_shot_draft}

Example rewrite:
{few_shot_rewrite}

Key facts that MUST appear unchanged:
{key_facts}

Now rewrite this draft:
{template}"""


class BaseModel(ABC):
    @abstractmethod
    def name(self) -> str:
        raise NotImplementedError

    @abstractmethod
    def generate(self, payload: dict, *, template: str) -> str:
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

    def generate(self, payload: dict, *, template: str) -> str:
        self._load()

        if "reads" not in payload and "results" in payload:
            payload = build_explain_payload(
                metrics=payload.get("results") or payload.get("results_r1"),
                scout_report=payload["scout_report"],
                sample_plan=payload.get("sample_plan", {}),
                fastq_path=payload["fastq_path"],
                r2_metrics=payload.get("r2_metrics"),
            )

        verdict = payload.get("verdict", "UNKNOWN")
        key_facts = "\n".join(f"- {fact}" for fact in extract_key_facts(payload)) or "- (see draft)"

        user_content = USER_PROMPT_TEMPLATE.format(
            verdict=verdict,
            few_shot_draft=FEW_SHOT_DRAFT,
            few_shot_rewrite=FEW_SHOT_REWRITE,
            key_facts=key_facts,
            template=template,
        )

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
            repetition_penalty=1.2,
        )
        new_tokens = generated_ids[0][input_len:]
        return self._tokenizer.decode(new_tokens, skip_special_tokens=True).strip()
