"""Build a compact JSON payload for LLM report explanations."""

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
    """Structured facts for the LLM narrator (no raw reads or plot images)."""
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
