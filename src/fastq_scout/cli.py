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
)
from fastq_scout.pipeline import Pipeline
from fastq_scout.reader import FastqReader
from fastq_scout.report import HtmlReport, build_plot_paths
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
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.fastq.is_file():
        raise SystemExit(f"Error: FASTQ file not found: {args.fastq}")

    if args.chunk_size <= 0:
        raise SystemExit("Error: --chunk-size must be a positive integer")

    html_path = args.output or args.fastq.with_name(f"{args.fastq.stem}_scout.html")

    reader = FastqReader(args.fastq, chunk_size=args.chunk_size)
    pipeline = Pipeline(metrics=[
        PerPositionQuality(),
        PerSequenceQuality(),
        LengthDistribution(),
        GCContent(),
        PerBaseSequenceContent(),
        DuplicateRate(),
    ])

    print("Start processing...")
    results = pipeline.run(reader)

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

        HtmlReport(args.fastq, results, scout_report, plot_paths).save(html_path)

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(
                {
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
