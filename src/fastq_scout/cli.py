import argparse
import json
from pathlib import Path

from fastq_scout.metrics import PerPositionQuality, LengthDistribution, GCContent, DuplicateRate
from fastq_scout.pipeline import Pipeline
from fastq_scout.plot import MetricPlotter
from fastq_scout.reader import FastqReader


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
        help="Path for JSON report (default: <input_stem>_scout.json next to input file)",
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
        help="Directory for PNG plots (skip plotting if omitted)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.fastq.is_file():
        raise SystemExit(f"Error: FASTQ file not found: {args.fastq}")

    if args.chunk_size <= 0:
        raise SystemExit("Error: --chunk-size must be a positive integer")

    output_path = args.output or args.fastq.with_name(f"{args.fastq.stem}_scout.json")

    reader = FastqReader(args.fastq, chunk_size=args.chunk_size)
    pipeline = Pipeline(metrics=[
        PerPositionQuality(),
        LengthDistribution(),
        GCContent(),
        DuplicateRate(),
    ])

    print("Start processing...")
    results = pipeline.run(reader)

    print("\n--- Processing Results ---")
    for metric_name, data in results.items():
        print(f"\n{metric_name}:")
        print(data)

    if args.plot_dir:
        print(f"\n--- Plotting Results ({args.plot_dir}) ---")
        for metric_name, data in results.items():
            plot_path = MetricPlotter(metric_name, data).plot(args.plot_dir)
            print(f"  {plot_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
