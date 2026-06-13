"""Build a compact JSON payload and plain-language report summaries."""

from __future__ import annotations

import json
import re
from pathlib import Path


def _quality_profile(quality: dict) -> dict:
    per_position = quality.get("per_position_mean", [])
    if len(per_position) < 10:
        return {}

    head_mean = round(sum(per_position[:10]) / 10, 2)
    tail_start = int(len(per_position) * 0.8)
    tail_mean = round(sum(per_position[tail_start:]) / len(per_position[tail_start:]), 2)
    return {
        "overall_mean": quality.get("overall_mean"),
        "head_mean_first_10bp": head_mean,
        "tail_mean_last_20pct": tail_mean,
        "tail_drop_phred": round(head_mean - tail_mean, 2),
    }


def _gc_profile(gc: dict) -> dict:
    histogram = gc.get("gc_histogram", [])
    if not histogram:
        return {"mean_gc": gc.get("mean_gc"), "std_gc": gc.get("std_gc")}

    peak_bin = max(range(len(histogram)), key=lambda i: histogram[i])
    significant_peaks = sum(1 for count in histogram if count > sum(histogram) * 0.05)
    return {
        "mean_gc": gc.get("mean_gc"),
        "std_gc": gc.get("std_gc"),
        "main_peak_pct_gc": peak_bin,
        "peak_count_significant": significant_peaks,
        "shape": "bimodal" if significant_peaks >= 2 else "unimodal",
    }


def _summarize_read(metrics: dict, label: str = "R1") -> dict:
    quality = metrics.get("Per position quality", {})
    seq_quality = metrics.get("Per sequence quality scores", {})
    length = metrics.get("Length distribution", {})
    gc = metrics.get("GC content", {})
    base_content = metrics.get("Per base sequence content", {})
    adapter = metrics.get("Adapter discovery", {})
    cpg = metrics.get("CpG O/E ratio", {})

    summary = {
        "read_label": label,
        "mean_quality_phred": quality.get("overall_mean"),
        "q20_pct": seq_quality.get("q20_pct"),
        "q30_pct": seq_quality.get("q30_pct"),
        "mean_read_length": length.get("mean_length"),
        "read_length_range_bp": [
            length.get("min_length"),
            length.get("max_length"),
        ],
        "duplicate_pct": metrics.get("Duplicates rate"),
        "quality_profile": _quality_profile(quality),
        "gc_profile": _gc_profile(gc),
        "cpg_oe_ratio": cpg.get("cpg_oe_ratio"),
        "composition_class": cpg.get("composition_class"),
        "base_content_summary": base_content.get("summary", {}),
    }

    if adapter:
        summary["adapter"] = {
            "detection_method": adapter.get("detection_method"),
            "reference_name": adapter.get("reference_name"),
            "adapter_content_pct": adapter.get("adapter_content_pct"),
            "trim_sequence": adapter.get("trim_sequence") or adapter.get("consensus"),
        }

    return summary


def build_explain_payload(
    metrics: dict,
    scout_report: dict,
    sample_plan: dict,
    fastq_path: Path,
    r2_metrics: dict | None = None,
) -> dict:
    """Structured facts for the narrator (no raw reads or plot images)."""
    reads = [_summarize_read(metrics, "R1")]
    if r2_metrics:
        reads.append(_summarize_read(r2_metrics, "R2"))

    sampling = {
        "reads_analyzed": sample_plan.get("reads_processed"),
        "total_reads_in_file": sample_plan.get("total_reads"),
        "sample_fraction_pct": sample_plan.get("sample_fraction_pct"),
        "mode": sample_plan.get("mode"),
    }
    if "R1" in sample_plan:
        sampling["R1"] = {
            "reads_analyzed": sample_plan["R1"].get("reads_processed"),
            "total_reads_in_file": sample_plan["R1"].get("total_reads"),
        }
        sampling["R2"] = {
            "reads_analyzed": sample_plan["R2"].get("reads_processed"),
            "total_reads_in_file": sample_plan["R2"].get("total_reads"),
        }

    return {
        "tool": "FastqScout",
        "input_fastq": fastq_path.name,
        "layout": scout_report.get("layout", sample_plan.get("layout", "single")),
        "library_type": scout_report.get("library_type", sample_plan.get("library_type", "genome")),
        "expected_species": scout_report.get("expected_species"),
        "verdict": scout_report.get("verdict"),
        "issues": scout_report.get("issues", []),
        "recommendations": scout_report.get("recommendations", []),
        "composition": scout_report.get("composition", {}),
        "sampling": sampling,
        "reads": reads,
    }


def payload_to_prompt_text(payload: dict) -> str:
    return json.dumps(payload, indent=2, ensure_ascii=False)


_VERDICT_TEXT = {
    "PROCEED": (
        "The sample passed pre-flight QC — you can continue with downstream analysis "
        "(alignment, variant calling, etc.) without mandatory preprocessing."
    ),
    "TRIM": (
        "The sample needs preprocessing before main analysis — "
        "typically adapter and/or low-quality end trimming (e.g. fastp)."
    ),
    "REJECT": (
        "Quality is too low to proceed confidently. "
        "Investigate the cause or consider re-sequencing."
    ),
}


def _read_good_points(read: dict) -> list[str]:
    points = []
    label = read.get("read_label", "R1")
    prefix = f"{label}: " if label else ""

    mean_q = read.get("mean_quality_phred")
    if mean_q is not None and mean_q >= 25:
        points.append(f"{prefix}good average read quality (PHRED {mean_q})")

    q20 = read.get("q20_pct")
    if q20 is not None and q20 >= 90:
        points.append(f"{prefix}{q20}% of bases at Q≥20 (error rate ≤1%)")

    dup = read.get("duplicate_pct")
    if dup is not None and dup < 20:
        points.append(f"{prefix}duplicate rate is acceptable ({dup}%)")

    gc = read.get("gc_profile", {})
    mean_gc = gc.get("mean_gc")
    if mean_gc is not None and 30 <= mean_gc <= 70 and gc.get("shape") != "bimodal":
        points.append(f"{prefix}typical GC content ({mean_gc}%)")

    qp = read.get("quality_profile", {})
    tail_drop = qp.get("tail_drop_phred", 0)
    if tail_drop is not None and tail_drop < 3:
        points.append(f"{prefix}quality is stable along read length")

    adapter = read.get("adapter")
    if adapter and adapter.get("adapter_content_pct", 0) < 1:
        points.append(f"{prefix}no strong adapter signal detected")

    return points


def _read_concern_points(read: dict) -> list[str]:
    points = []
    label = read.get("read_label", "R1")
    prefix = f"{label}: " if label else ""

    qp = read.get("quality_profile", {})
    tail_drop = qp.get("tail_drop_phred", 0)
    if tail_drop and tail_drop >= 5:
        points.append(
            f"{prefix}quality drops toward read ends "
            f"(start ~{qp.get('head_mean_first_10bp')}, "
            f"tail ~{qp.get('tail_mean_last_20pct')} PHRED)"
        )

    adapter = read.get("adapter")
    if adapter:
        pct = adapter.get("adapter_content_pct", 0)
        if pct >= 5:
            name = adapter.get("reference_name") or "adapter"
            points.append(
                f"{prefix}{name} detected on read tails ({pct}% of reads) — "
                "trim before alignment"
            )
        elif pct >= 0.5:
            points.append(f"{prefix}low-level adapter signal ({pct}%)")

    dup = read.get("duplicate_pct")
    if dup is not None and dup >= 30:
        points.append(f"{prefix}high duplicate rate ({dup}%)")

    gc = read.get("gc_profile", {})
    if gc.get("shape") == "bimodal":
        points.append(f"{prefix}bimodal GC distribution — possible mix or contamination")

    mean_gc = gc.get("mean_gc")
    if mean_gc is not None and (mean_gc < 30 or mean_gc > 70):
        points.append(f"{prefix}unusual mean GC ({mean_gc}%)")

    base = read.get("base_content_summary", {})
    if base.get("mean_n_pct", 0) > 1:
        points.append(f"{prefix}elevated unknown-base (N) content ({base['mean_n_pct']}%)")

    cpg_oe = read.get("cpg_oe_ratio")
    if cpg_oe is not None and cpg_oe >= 0.85:
        points.append(
            f"{prefix}high CpG O/E ({cpg_oe}) — not typical of a mammalian methylated genome"
        )

    return points


def _metrics_for_plot_title(
    plot_title: str,
    metrics_r1: dict,
    metrics_r2: dict | None,
) -> tuple[str, dict]:
    """Return (read label, metrics dict) for a report plot title."""
    title = plot_title
    label = "R1"
    if title.startswith("R2: "):
        label = "R2"
        title = title[4:]
        source = metrics_r2 or {}
    else:
        source = metrics_r1

    if title == "Duplicate rate":
        title = "Duplicates rate"
    elif title == "Adapter enrichment":
        title = "Adapter discovery"

    return label, source


def _caption_per_position_quality(quality: dict, label: str) -> str:
    overall = quality.get("overall_mean")
    qp = _quality_profile(quality)
    parts = [
        f"This chart shows the average base quality (PHRED score) at each position along the "
        f"{label} reads — higher is better; PHRED 30 means roughly 1 error per 1,000 bases.",
    ]
    if overall is not None:
        parts.append(f"The overall mean quality is {overall}.")
    tail_drop = qp.get("tail_drop_phred", 0)
    if tail_drop and tail_drop >= 5:
        parts.append(
            f"Quality falls toward the read ends (about {qp.get('head_mean_first_10bp')} at the "
            f"start vs {qp.get('tail_mean_last_20pct')} near the tail) — trimming low-quality "
            f"ends often helps before alignment."
        )
    elif tail_drop is not None and tail_drop < 3:
        parts.append("Quality stays fairly stable along the read — no strong end decay.")
    return " ".join(parts)


def _caption_per_sequence_quality(seq_quality: dict, label: str) -> str:
    q20 = seq_quality.get("q20_pct")
    q30 = seq_quality.get("q30_pct")
    parts = [
        f"This histogram summarizes the average quality of whole {label} reads. "
        f"Q20 and Q30 tell you what fraction of bases meet common quality cutoffs "
        f"(≤1% and ≤0.1% error rate)."
    ]
    if q20 is not None and q30 is not None:
        parts.append(f"In this sample, {q20}% of bases are Q≥20 and {q30}% are Q≥30.")
    if q20 is not None and q20 >= 90:
        parts.append("Most bases pass typical Q20 filters — good for standard pipelines.")
    elif q20 is not None and q20 < 80:
        parts.append("A sizable fraction of bases are below Q20 — consider trimming or filtering.")
    return " ".join(parts)


def _caption_length_distribution(length: dict, label: str) -> str:
    mean_len = length.get("mean_length")
    min_len = length.get("min_length")
    max_len = length.get("max_length")
    parts = [
        f"This plot shows how long the {label} reads are. "
        f"Unexpectedly short reads can mean adapter dimers or library problems."
    ]
    if mean_len is not None:
        parts.append(f"Mean read length is {mean_len} bp")
        if min_len is not None and max_len is not None:
            parts.append(f"(range {min_len}–{max_len} bp).")
        else:
            parts.append(".")
    return " ".join(parts)


def _caption_gc_content(gc: dict, label: str) -> str:
    profile = _gc_profile(gc)
    mean_gc = profile.get("mean_gc")
    parts = [
        f"GC content is the percentage of G and C bases in each {label} read. "
        f"Most genomes have a characteristic GC profile — odd shapes can hint at contamination or mixed libraries."
    ]
    if mean_gc is not None:
        parts.append(f"Mean GC here is {mean_gc}%.")
    if profile.get("shape") == "bimodal":
        parts.append(
            "The distribution has two peaks (bimodal) — worth checking whether the sample "
            "or metadata matches what you expect."
        )
    else:
        parts.append("The distribution looks single-peaked (unimodal).")
    return " ".join(parts)


def _caption_base_content(base_content: dict, label: str) -> str:
    summary = base_content.get("summary", {})
    mean_n = summary.get("mean_n_pct", 0)
    parts = [
        f"This chart tracks A, T, G, C (and N / other symbols) across read positions for {label}. "
        f"Flat lines are normal; spikes in N mean base-calling problems at those cycles."
    ]
    if mean_n and mean_n > 1:
        parts.append(f"Unknown bases (N) average {mean_n}% — filtering or trimming may be needed.")
    else:
        parts.append("N content looks low — base composition is stable.")
    return " ".join(parts)


def _caption_duplicates(rate, label: str) -> str:
    parts = [
        f"Duplicate rate estimates how many {label} reads are identical or near-identical "
        f"(often from PCR over-amplification or very high sequencing depth)."
    ]
    if rate is not None:
        parts.append(f"Here it is {rate}%.")
    if rate is not None and rate < 20:
        parts.append("This level is usually acceptable for whole-genome work.")
    elif rate is not None and rate >= 30:
        parts.append("High duplicates can skew variant and expression calls — deduplication may help.")
    else:
        parts.append("Moderate duplication — note it when interpreting downstream results.")
    return " ".join(parts)


def _caption_adapter(adapter: dict, label: str) -> str:
    if not adapter:
        return (
            f"Adapter enrichment was not run for {label}. "
            f"Re-run FastqScout with --mode with_adapter to scan read tails for sequencing adapters."
        )
    pct = adapter.get("adapter_content_pct", 0)
    name = adapter.get("reference_name") or "a sequencing adapter"
    parts = [
        f"This plot shows how often adapter sequence appears on {label} read tails. "
        f"Adapters are lab oligos that should be removed before mapping reads to a genome."
    ]
    parts.append(f"{name} was detected on about {pct}% of tails.")
    if pct >= 5:
        parts.append("Trim these with fastp (or similar) using the suggested sequence in the report.")
    elif pct >= 0.5:
        parts.append("Signal is low but present — trimming is optional but often still worthwhile.")
    else:
        parts.append("Little or no adapter signal — trimming may not be necessary.")
    return " ".join(parts)


def _caption_for_metric(metric_name: str, metrics: dict, label: str) -> str:
    if metric_name == "Per position quality":
        return _caption_per_position_quality(metrics.get("Per position quality", {}), label)
    if metric_name == "Per sequence quality scores":
        return _caption_per_sequence_quality(metrics.get("Per sequence quality scores", {}), label)
    if metric_name == "Length distribution":
        return _caption_length_distribution(metrics.get("Length distribution", {}), label)
    if metric_name == "GC content":
        return _caption_gc_content(metrics.get("GC content", {}), label)
    if metric_name == "Per base sequence content":
        return _caption_base_content(metrics.get("Per base sequence content", {}), label)
    if metric_name == "Duplicates rate":
        return _caption_duplicates(metrics.get("Duplicates rate"), label)
    if metric_name == "Adapter discovery":
        return _caption_adapter(metrics.get("Adapter discovery", {}), label)
    return ""


def build_plot_captions(
    metrics_r1: dict,
    plot_titles: list[str],
    metrics_r2: dict | None = None,
) -> dict[str, str]:
    """Plain-language captions keyed by report plot title (always English)."""
    captions: dict[str, str] = {}
    for plot_title in plot_titles:
        label, source = _metrics_for_plot_title(plot_title, metrics_r1, metrics_r2)
        title = plot_title[4:] if plot_title.startswith("R2: ") else plot_title
        if title == "Duplicate rate":
            metric_name = "Duplicates rate"
        elif title == "Adapter enrichment":
            metric_name = "Adapter discovery"
        else:
            metric_name = title
        caption = _caption_for_metric(metric_name, source, label)
        if caption:
            captions[plot_title] = caption
    return captions


_VERDICT_EXTRA = {
    "PROCEED": (
        "In practice, you can hand this FASTQ to your usual pipeline without mandatory cleanup."
    ),
    "TRIM": (
        "Plan a trimming step (usually fastp) before alignment — the plots below show what to fix."
    ),
    "REJECT": (
        "Do not start a long downstream run until you understand why quality failed."
    ),
}


def build_template_explanation(payload: dict) -> str:
    """Rule-based plain-language summary — reliable without LLM."""
    verdict = payload.get("verdict", "UNKNOWN")
    layout = payload.get("layout", "single")
    library = payload.get("library_type", "genome")
    sampling = payload.get("sampling", {})
    reads_analyzed = sampling.get("reads_analyzed")
    total_reads = sampling.get("total_reads_in_file")
    fraction = sampling.get("sample_fraction_pct")

    lines = ["## Summary", ""]

    verdict_line = _VERDICT_TEXT.get(verdict, f"FastqScout verdict: {verdict}.")
    extra = _VERDICT_EXTRA.get(verdict, "")
    lines.append(f"Verdict: {verdict}. {verdict_line}")
    if extra:
        lines.append(extra)
    lines.append(
        "This is a quick pre-flight check on a sample of reads — use it to catch adapter, "
        "quality, or metadata problems before expensive analysis."
    )

    if reads_analyzed and total_reads and fraction is not None:
        lines.append(
            f"Analyzed {reads_analyzed:,} reads out of {total_reads:,} in the file "
            f"({fraction}% sample) — layout: {layout}, library: {library}."
        )
    else:
        lines.append(f"Layout: {layout}, library type: {library}.")

    composition = payload.get("composition", {})
    if composition.get("summary"):
        lines.extend(["", "## Genome composition (CpG O/E)", ""])
        lines.append(composition["summary"])
        species_check = composition.get("expected_species_check", {})
        if species_check.get("message"):
            lines.append(species_check["message"])
        lines.append(
            "This is a composition hint, not species identification — "
            "use Kraken/sourmash for taxonomy."
        )
        lines.append(
            "CpG O/E compares observed CG dinucleotides to random expectation — "
            "mammalian genomes are often depleted (~0.2); bacteria and some invertebrates are closer to 1.0."
        )

    lines.extend(["", "## What looks good", ""])
    lines.append(
        "These metrics passed basic checks and support continuing — they do not mean the sample is perfect."
    )
    good = []
    for read in payload.get("reads", []):
        good.extend(_read_good_points(read))
    if good:
        lines.extend(f"- {p}" for p in good)
    else:
        lines.append("- Few strong positives in the metrics — see Issues above.")

    lines.extend(["", "## What to watch", ""])
    lines.append(
        "Each point below may affect mapping or variant calls if left unaddressed."
    )
    concerns = list(payload.get("issues", []))
    if concerns:
        lines.extend(f"- {c}" for c in concerns)
    else:
        lines.append("- No major issues detected.")

    lines.extend(["", "## Next steps", ""])
    lines.append(
        "Follow these in order; share this report with bioinformatics if anything is unclear."
    )
    recs = payload.get("recommendations", [])
    if recs:
        lines.extend(f"- {r}" for r in recs)
    elif verdict == "PROCEED":
        lines.append("- Continue with your main pipeline (alignment, downstream QC).")
    elif verdict == "TRIM":
        lines.append("- Run fastp (or similar) as recommended above, then re-check QC.")
    else:
        lines.append("- Discuss the result with bioinformatics or your core facility.")

    return "\n".join(lines)


def extract_key_facts(payload: dict) -> list[str]:
    """Numbers and labels the LLM must preserve (for validation)."""
    facts: list[str] = []
    verdict = payload.get("verdict")
    if verdict:
        facts.append(verdict)

    for read in payload.get("reads", []):
        for key in ("mean_quality_phred", "q20_pct", "q30_pct", "duplicate_pct"):
            val = read.get(key)
            if val is not None:
                facts.append(str(val))
        adapter = read.get("adapter") or {}
        pct = adapter.get("adapter_content_pct")
        if pct is not None:
            facts.append(str(pct))

    composition = payload.get("composition", {})
    cpg_oe = composition.get("cpg_oe_ratio")
    if cpg_oe is not None:
        facts.append(str(cpg_oe))

    sampling = payload.get("sampling", {})
    fraction = sampling.get("sample_fraction_pct")
    if fraction is not None:
        facts.append(str(fraction))

    return facts


_BAD_LLM_PATTERNS = re.compile(
    r"grammar fragment|DNA fragment member|chromosome member|"
    r"as an AI|I cannot|I'm sorry",
    re.IGNORECASE,
)

_REQUIRED_LLM_HEADINGS = (
    "## Summary",
    "## What looks good",
    "## What to watch",
    "## Next steps",
)


def is_llm_response_usable(text: str, payload: dict) -> bool:
    if not text or len(text) < 100:
        return False
    if _BAD_LLM_PATTERNS.search(text):
        return False
    if re.search(r"[\u0400-\u04FF]", text):
        return False

    if not all(h in text for h in _REQUIRED_LLM_HEADINGS):
        return False

    verdict = payload.get("verdict", "")
    if verdict and verdict not in text.upper():
        return False

    facts = extract_key_facts(payload)
    if facts:
        hits = sum(1 for fact in facts if fact in text)
        if hits < max(2, len(facts) // 2):
            return False

    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if len(lines) > 5:
        from collections import Counter

        counts = Counter(lines)
        if any(c >= 3 for c in counts.values()):
            return False

    return True


def generate_explanation(
    payload: dict,
    *,
    use_llm: bool = False,
    model_name: str = "Qwen/Qwen2.5-0.5B-Instruct",
) -> tuple[str, str]:
    """
    Return (explanation_text, source) where source is 'template' or 'llm'.
    Default: template. LLM is optional and falls back to template if output is bad.
    """
    template = build_template_explanation(payload)
    if not use_llm:
        return template, "template"

    from fastq_scout.qween_model import QwenModel

    try:
        llm_text = QwenModel(model_name=model_name).generate(payload, template=template)
        if is_llm_response_usable(llm_text, payload):
            return llm_text, "llm"
        print("LLM output looked unreliable — using rule-based summary instead.")
    except Exception as exc:
        print(f"LLM explanation failed ({exc}) — using rule-based summary.")

    return template, "template"
