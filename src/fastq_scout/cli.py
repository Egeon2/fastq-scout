import argparse
import json
import shutil
import tempfile
from pathlib import Path

from fastq_scout.metrics import (
    PerPositionQuality,
    LengthDistribution,
    GCContent,
    PerBaseSequenceContent,
    DuplicateRate,
    PerSequenceQuality,
    AdapterDiscovery,
)
from fastq_scout.pipeline import Pipeline
from fastq_scout.reader import FastqReader
from fastq_scout.report import HtmlReport, build_plot_paths
from fastq_scout.sampling import count_fastq_reads, resolve_sample_budget
from fastq_scout.scout import FastqScout


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
        help="Path to input FASTQ file",
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
        help="Minimum reads for auto sampling (default: 100000)",
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
    return parser


def _build_sample_plan(args: argparse.Namespace) -> tuple[dict, int | None]:
    total_reads = count_fastq_reads(args.fastq)

    if args.full_file:
        return {
            "total_reads": total_reads,
            "sample_budget": total_reads,
            "sample_budget_base": total_reads,
            "sample_budget_for_adapters": total_reads,
            "sample_fraction_pct": 100.0 if total_reads else 0.0,
            "sample_fraction_pct_for_adapters": 100.0 if total_reads else 0.0,
            "mode": "full",
        }, None

    if args.max_reads is not None:
        if args.max_reads <= 0:
            raise SystemExit("Error: --max-reads must be a positive integer")
        sample_budget = min(args.max_reads, total_reads)
        return {
            "total_reads": total_reads,
            "sample_budget": sample_budget,
            "sample_budget_base": sample_budget,
            "sample_budget_for_adapters": sample_budget,
            "sample_fraction_pct": round(sample_budget / total_reads * 100, 2) if total_reads else 0.0,
            "sample_fraction_pct_for_adapters": round(sample_budget / total_reads * 100, 2) if total_reads else 0.0,
            "mode": "max_reads",
        }, sample_budget

    sample_plan = resolve_sample_budget(
        total_reads,
        margin_rate=args.margin_rate,
        margin_mean=args.margin_mean,
        min_reads=args.min_reads,
        max_fraction=args.max_fraction,
        mode=args.mode,
    )
    return sample_plan, sample_plan["sample_budget"]


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

    html_path = args.output or args.fastq.with_name(f"{args.fastq.stem}_scout.html")

    sample_plan, sample_budget = _build_sample_plan(args)
    reader = FastqReader(args.fastq, chunk_size=args.chunk_size, sample_budget=sample_budget)

    metrics = [
        PerPositionQuality(),
        PerSequenceQuality(),
        LengthDistribution(),
        GCContent(),
        PerBaseSequenceContent(),
        DuplicateRate(),
    ]
    if args.mode == "with_adapter":
        metrics.append(AdapterDiscovery())

    pipeline = Pipeline(metrics=metrics)

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
    results = pipeline.run(reader)

    sample_plan["reads_processed"] = reader.reads_processed
    print(f"Processed {reader.reads_processed:,} reads")

    scout = FastqScout()
    verdict, scout_report = scout.result(results)

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

        HtmlReport(
            args.fastq,
            results,
            scout_report,
            plot_paths,
            sample_plan=sample_plan,
        ).save(html_path)

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(
                {
                    "sampling": sample_plan,
                    "metrics": results,
                    "scout": scout_report,
                },
                f,
                indent=2,
            )
        print(f"JSON saved to {args.json}")

    print(f"HTML report saved to {html_path}")
    return VERDICT_EXIT_CODES.get(verdict, 0)


if __name__ == "__main__":
    raise SystemExit(main())
