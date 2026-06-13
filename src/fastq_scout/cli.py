import argparse
import json
import shutil
import tempfile
from pathlib import Path
import time

from fastq_scout.metrics import (
    PerPositionQuality,
    LengthDistribution,
    GCContent,
    PerBaseSequenceContent,
    DuplicateRate,
    PerSequenceQuality,
    AdapterDiscovery,
    CpGObservedExpected,
)
from fastq_scout.pipeline import Pipeline
from fastq_scout.profiles import effective_min_reads
from fastq_scout.run_context import RunContext, build_context
from fastq_scout.reader import FastqReader
from fastq_scout.report import HtmlReport, build_plot_paths
from fastq_scout.sampling import count_fastq_reads, resolve_sample_budget
from fastq_scout.scout import FastqScout
from fastq_scout.explain import build_explain_payload, generate_explanation


VERDICT_EXIT_CODES = {
    "PROCEED": 0,
    "TRIM": 1,
    "REJECT": 2,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fastq-scout",
        description="FastqScout: lightweight QC scout for FASTQ files",
    )
    parser.add_argument(
        "fastq",
        type=Path,
        help="Path to R1 FASTQ file (or single-end FASTQ)",
    )
    parser.add_argument(
        "fastq_r2",
        type=Path,
        nargs="?",
        help="R2 FASTQ (required when --layout paired)",
    )
    parser.add_argument(
        "--layout",
        choices=["single", "paired"],
        default="single",
        help="Read layout (default: single)",
    )
    parser.add_argument(
        "--library-type",
        choices=["genome", "transcriptome"],
        default="genome",
        help="Library type (default: genome)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Path for HTML report (default: <input_stem>_scout.html next to input file)",
    )
    parser.add_argument(
        "--json",
        type=Path,
        help="Optional path for JSON metrics export",
    )
    parser.add_argument(
        "-c", "--chunk-size",
        type=int,
        default=10_000,
        help="Number of reads per chunk (default: 10000)",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        help="Keep PNG plots in this directory (plots are also embedded in HTML)",
    )
    parser.add_argument(
        "--full-file",
        action="store_true",
        help="Analyze the entire FASTQ file (disable auto sampling)",
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        help="Hard cap on reads to analyze (overrides auto sample budget)",
    )
    parser.add_argument(
        "--max-fraction",
        type=float,
        default=0.15,
        help="Maximum fraction of reads to analyze (default: 0.15)",
    )
    parser.add_argument(
        "--min-reads",
        type=int,
        default=10_000,
        help="Minimum reads for auto sampling (default: 10000)",
    )
    parser.add_argument(
        "--margin-rate",
        type=float,
        default=0.05,
        help="Target margin for rate metrics, e.g. duplicate rate (default: 0.05)",
    )
    parser.add_argument(
        "--margin-mean",
        type=float,
        default=1.0,
        help="Target margin for mean PHRED score (default: 1.0)",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="base",
        help="Sampling mode: base or with_adapter (default: base)",
    )
    parser.add_argument(
        "--expected-species",
        choices=["human", "mouse", "drosophila", "ecoli"],
        default=None,
        help=(
            "Optional metadata sanity check: warn if CpG O/E and GC "
            "do not match typical genome composition for this species"
        ),
    )
    parser.add_argument(
        "--explain",
        action="store_true",
        help="Add plain-language summary to HTML (rule-based template, English)",
    )
    parser.add_argument(
        "--explain-llm",
        action="store_true",
        help="Try local Qwen for summary; falls back to template if output is bad",
    )
    parser.add_argument(
        "--llm-model",
        default="Qwen/Qwen2.5-0.5B-Instruct",
        help="Hugging Face model for --explain-llm (default: Qwen2.5-0.5B-Instruct)",
    )
    return parser


def _build_metrics(mode: str, read_label: str) -> list:
    metrics = [
        PerPositionQuality(),
        PerSequenceQuality(),
        LengthDistribution(),
        GCContent(),
        PerBaseSequenceContent(),
        DuplicateRate(),
        CpGObservedExpected(),
    ]
    if mode == "with_adapter":
        adapter_set = "universal" if read_label == "R1" else "read2"
        metrics.append(AdapterDiscovery(adapter_set=adapter_set))
    return metrics


def _print_adapter_summary(results: dict, read_label: str = "") -> None:
    adapter = results.get("Adapter discovery", {})
    if not adapter or adapter.get("detection_method", "none") == "none":
        return

    prefix = f"[{read_label}] " if read_label else ""
    pct = adapter.get("adapter_content_pct", 0)
    ref_name = adapter.get("reference_name") or "Unknown"
    trim_seq = adapter.get("trim_sequence") or adapter.get("consensus", "")
    ref_seq = adapter.get("reference_sequence", "")
    method = adapter.get("detection_method", "none")

    print(f"{prefix}Adapter detected ({method}): {ref_name} — {pct}% of read tails")
    if trim_seq:
        print(f"{prefix}  fastp trim sequence: {trim_seq}")
    if ref_seq:
        print(f"{prefix}  full reference:      {ref_seq}")
    elif trim_seq:
        print(f"{prefix}  detected motif:      {trim_seq}")


def _print_composition_summary(results: dict, scout_report: dict, read_label: str = "") -> None:
    cpg = results.get("CpG O/E ratio", {})
    ratio = cpg.get("cpg_oe_ratio")
    if ratio is None:
        return

    prefix = f"[{read_label}] " if read_label else ""
    print(f"{prefix}CpG O/E ratio: {ratio} (observed/expected CG dinucleotides)")
    composition = scout_report.get("composition", {})
    if composition.get("summary"):
        print(f"{prefix}  {composition['summary']}")
    species_check = composition.get("expected_species_check", {})
    if species_check.get("message"):
        print(f"{prefix}  Species check: {species_check['message']}")


def _build_sample_plan(
    args: argparse.Namespace,
    ctx: RunContext,
    fastq: Path,
) -> tuple[dict, int | None]:
    total_reads = count_fastq_reads(fastq)
    min_reads = effective_min_reads(ctx, args.min_reads)

    if args.full_file:
        plan = {
            "total_reads": total_reads,
            "sample_budget": total_reads,
            "sample_budget_base": total_reads,
            "sample_budget_for_adapters": total_reads,
            "sample_fraction_pct": 100.0 if total_reads else 0.0,
            "sample_fraction_pct_for_adapters": 100.0 if total_reads else 0.0,
            "mode": "full",
        }
        return plan, None

    if args.max_reads is not None:
        if args.max_reads <= 0:
            raise SystemExit("Error: --max-reads must be a positive integer")
        sample_budget = min(args.max_reads, total_reads)
        plan = {
            "total_reads": total_reads,
            "sample_budget": sample_budget,
            "sample_budget_base": sample_budget,
            "sample_budget_for_adapters": sample_budget,
            "sample_fraction_pct": round(sample_budget / total_reads * 100, 2) if total_reads else 0.0,
            "sample_fraction_pct_for_adapters": round(sample_budget / total_reads * 100, 2) if total_reads else 0.0,
            "mode": "max_reads",
        }
        return plan, sample_budget

    plan = resolve_sample_budget(
        total_reads,
        margin_rate=args.margin_rate,
        margin_mean=args.margin_mean,
        min_reads=min_reads,
        max_fraction=args.max_fraction,
        mode=args.mode,
    )
    return plan, plan["sample_budget"]


def _annotate_sample_plan(plan: dict, ctx: RunContext, fastq: Path, read_label: str) -> dict:
    annotated = dict(plan)
    annotated["layout"] = ctx.layout
    annotated["library_type"] = ctx.library_type
    annotated["read_label"] = read_label
    annotated["fastq"] = str(fastq)
    return annotated


def _run_fastq(
    fastq: Path,
    read_label: str,
    args: argparse.Namespace,
    ctx: RunContext,
    sample_plan: dict,
    sample_budget: int | None,
) -> tuple[dict, int, dict]:
    reader = FastqReader(fastq, chunk_size=args.chunk_size, sample_budget=sample_budget)
    pipeline = Pipeline(metrics=_build_metrics(args.mode, read_label))
    results = pipeline.run(reader)
    plan = _annotate_sample_plan(sample_plan, ctx, fastq, read_label)
    plan["reads_processed"] = reader.reads_processed
    return results, reader.reads_processed, plan


def main(argv: list[str] | None = None) -> int:
    
    args = build_parser().parse_args(argv)

    if not args.fastq.is_file():
        raise SystemExit(f"Error: FASTQ file not found: {args.fastq}")

    if args.chunk_size <= 0:
        raise SystemExit("Error: --chunk-size must be a positive integer")

    if args.max_fraction <= 0 or args.max_fraction > 1:
        raise SystemExit("Error: --max-fraction must be between 0 and 1")

    if args.mode not in ("base", "with_adapter"):
        raise SystemExit("Error: --mode must be 'base' or 'with_adapter'")

    if args.layout == "paired":
        if args.fastq_r2 is None:
            raise SystemExit("Error: --layout paired requires fastq_r2")
        if not args.fastq_r2.is_file():
            raise SystemExit(f"Error: R2 FASTQ file not found: {args.fastq_r2}")
    elif args.fastq_r2 is not None:
        raise SystemExit("Error: fastq_r2 provided but layout is single")

    ctx = build_context(args)
    scout = FastqScout(ctx)
    html_path = args.output or args.fastq.with_name(f"{args.fastq.stem}_scout.html")

    print(f"Layout: {ctx.layout} | Library: {ctx.library_type}")
    if ctx.expected_species:
        print(f"Expected species (composition check): {ctx.expected_species}")

    if ctx.is_paired:
        plan_r1, sample_budget_r1 = _build_sample_plan(args, ctx, args.fastq)
        plan_r2, sample_budget_r2 = _build_sample_plan(args, ctx, args.fastq_r2)

        print(f"R1: {args.fastq.name}")
        print(f"R2: {args.fastq_r2.name}")
        print("Start processing R1...")
        results_r1, processed_r1, plan_r1 = _run_fastq(
            args.fastq, "R1", args, ctx, plan_r1, sample_budget_r1
        )
        print(f"Processed R1: {processed_r1:,} reads")

        print("Start processing R2...")
        results_r2, processed_r2, plan_r2 = _run_fastq(
            args.fastq_r2, "R2", args, ctx, plan_r2, sample_budget_r2
        )
        print(f"Processed R2: {processed_r2:,} reads")

        verdict, scout_report = scout.result_paired(results_r1, results_r2)
        _print_adapter_summary(results_r1, "R1")
        _print_adapter_summary(results_r2, "R2")
        _print_composition_summary(results_r1, scout_report, "R1")

        combined_plan = {
            "layout": ctx.layout,
            "library_type": ctx.library_type,
            "R1": plan_r1,
            "R2": plan_r2,
        }

        with tempfile.TemporaryDirectory(prefix="fastq_scout_") as temp_dir:
            temp_path = Path(temp_dir)
            plot_paths = build_plot_paths(results_r1, temp_path, name_suffix="_r1")
            r2_plot_paths = build_plot_paths(results_r2, temp_path, name_suffix="_r2")
            plot_paths.update({f"R2: {title}": path for title, path in r2_plot_paths.items()})

            if args.plot_dir:
                args.plot_dir.mkdir(parents=True, exist_ok=True)
                saved_plots = {}
                for title, plot_path in plot_paths.items():
                    target = args.plot_dir / plot_path.name
                    shutil.copy2(plot_path, target)
                    saved_plots[title] = target
                plot_paths = saved_plots

            explanation = None
            explanation_source = None
            if args.explain or args.explain_llm:
                print("Generating plain-language explanation...")
                explain_payload = build_explain_payload(
                    results_r1,
                    scout_report,
                    combined_plan,
                    args.fastq,
                    r2_metrics=results_r2,
                )
                explanation, explanation_source = generate_explanation(
                    explain_payload,
                    use_llm=args.explain_llm,
                    model_name=args.llm_model,
                )
                print(f"Summary source: {explanation_source}")

            HtmlReport(
                args.fastq,
                results_r1,
                scout_report,
                plot_paths,
                sample_plan=combined_plan,
                fastq_r2=args.fastq_r2,
                r2_metrics=results_r2,
                explanation=explanation,
                explanation_source=explanation_source,
            ).save(html_path)

        json_metrics = {"R1": results_r1, "R2": results_r2}
        json_sampling = combined_plan
    else:
        sample_plan, sample_budget = _build_sample_plan(args, ctx, args.fastq)

        print(f"Total reads in file: {sample_plan['total_reads']:,}")
        print(f"Sampling mode: {sample_plan.get('mode', 'auto')}")
        if sample_budget is not None:
            print(
                f"Scout sample budget: {sample_plan['sample_budget']:,} reads "
                f"({sample_plan['sample_fraction_pct']}% of file)"
            )
            if sample_plan.get("mode") == "with_adapter":
                print(
                    f"Adapter-oriented budget: {sample_plan['sample_budget_for_adapters']:,} reads "
                    f"({sample_plan['sample_fraction_pct_for_adapters']}% of file)"
                )
        else:
            print("Scout mode: full file")

        print("Start processing...")
        results, processed, sample_plan = _run_fastq(
            args.fastq, "R1", args, ctx, sample_plan, sample_budget
        )
        print(f"Processed {processed:,} reads")

        verdict, scout_report = scout.result(results)
        if args.mode == "with_adapter":
            _print_adapter_summary(results)
        else:
            print("Adapter detection: off (use --mode with_adapter to enable)")
        _print_composition_summary(results, scout_report)

        with tempfile.TemporaryDirectory(prefix="fastq_scout_") as temp_dir:
            plot_paths = build_plot_paths(results, Path(temp_dir))

            if args.plot_dir:
                args.plot_dir.mkdir(parents=True, exist_ok=True)
                saved_plots = {}
                for title, plot_path in plot_paths.items():
                    target = args.plot_dir / plot_path.name
                    shutil.copy2(plot_path, target)
                    saved_plots[title] = target
                plot_paths = saved_plots

            explanation = None
            explanation_source = None
            if args.explain or args.explain_llm:
                print("Generating plain-language explanation...")
                explain_payload = build_explain_payload(
                    results,
                    scout_report,
                    sample_plan,
                    args.fastq,
                )
                explanation, explanation_source = generate_explanation(
                    explain_payload,
                    use_llm=args.explain_llm,
                    model_name=args.llm_model,
                )
                print(f"Summary source: {explanation_source}")

            HtmlReport(
                args.fastq,
                results,
                scout_report,
                plot_paths,
                sample_plan=sample_plan,
                explanation=explanation,
                explanation_source=explanation_source,
            ).save(html_path)

        json_metrics = results
        json_sampling = sample_plan

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(
                {
                    "run": {
                        "layout": ctx.layout,
                        "library_type": ctx.library_type,
                        "expected_species": ctx.expected_species,
                    },
                    "sampling": json_sampling,
                    "metrics": json_metrics,
                    "scout": scout_report,
                },
                f,
                indent=2,
            )
        print(f"JSON saved to {args.json}")

    print(f"HTML report saved to {html_path}")
    return VERDICT_EXIT_CODES.get(verdict, 0)

if __name__ == "__main__":
    start_time = time.time()
    exit_code = main()
    print(f"Time taken: {time.time() - start_time:.2f} seconds")
    raise SystemExit(exit_code)
