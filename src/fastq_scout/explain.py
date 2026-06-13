"""Build a compact JSON payload and plain-language report summaries."""

from __future__ import annotations

import json
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
    lines.append(f"Verdict: {verdict}. {verdict_line}")

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

    lines.extend(["", "## What looks good", ""])
    good = []
    for read in payload.get("reads", []):
        good.extend(_read_good_points(read))
    if good:
        lines.extend(f"- {p}" for p in good)
    else:
        lines.append("- Few strong positives in the metrics — see Issues above.")

    lines.extend(["", "## What to watch", ""])
    concerns = list(payload.get("issues", []))
    for read in payload.get("reads", []):
        for point in _read_concern_points(read):
            if point not in concerns:
                concerns.append(point)
    if concerns:
        lines.extend(f"- {c}" for c in concerns)
    else:
        lines.append("- No major issues detected.")

    lines.extend(["", "## Next steps", ""])
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

    from fastq_scout.qween_model import QwenModel, is_llm_response_usable

    try:
        llm_text = QwenModel(model_name=model_name).generate(payload)
        if is_llm_response_usable(llm_text, payload):
            return llm_text, "llm"
        print("LLM output looked unreliable — using rule-based summary instead.")
    except Exception as exc:
        print(f"LLM explanation failed ({exc}) — using rule-based summary.")

    return template, "template"
