# FastqScout

FastqScout is a lightweight Python toolkit for **pre-flight quality control** of FASTQ sequencing data. It scans reads in streaming chunks, calculates essential QC metrics, builds an HTML report with plots, and returns a **verdict** (`PROCEED`, `TRIM`, or `REJECT`) with actionable recommendations.

Unlike full diagnostic tools such as FastQC, FastqScout is designed as a **scout**: a fast first check that tells you whether a sample is ready for downstream analysis and what to do next.

## Features

- **Memory efficient** — stream-based chunking for large FASTQ files (GB-scale)
- **QC metrics** — per-position quality, GC content, read length distribution, duplicate rate
- **Pre-flight verdict** — automatic `PROCEED` / `TRIM` / `REJECT` decision with recommendations
- **HTML report** — single self-contained file with metrics, plots, issues, and next steps
- **CLI** — argparse interface with optional JSON export and PNG output
- **Pipeline-friendly exit codes** — `0` / `1` / `2` based on verdict

## Project structure

```text
fastq-scout/
├── docs/
│   └── images/              # Example plots for README
├── src/
│   └── fastq_scout/
│       ├── cli.py           # CLI entry point
│       ├── reader.py        # FASTQ streaming reader
│       ├── pipeline.py      # Analysis orchestration
│       ├── metrics.py       # Metric calculation logic
│       ├── models.py        # Data structures
│       ├── plot.py          # PNG plot generation
│       ├── scout.py         # Pre-flight verdict engine
│       └── report.py        # HTML report builder
└── README.md
```

## Installation

Requires Python 3.10+ and matplotlib:

```bash
pip install matplotlib
```

Clone the repository and run from the project root:

```bash
git clone <repo-url>
cd fastq-scout
```

## Usage

### Basic run

Generates `<input_stem>_scout.html` next to the input file:

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py your_sample.fastq
```

### Full run with custom output paths

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py your_sample.fastq \
  -o reports/sample_scout.html \
  --json reports/sample_scout.json \
  --plot-dir reports/plots \
  -c 10000
```

### CLI options

| Option | Description |
|--------|-------------|
| `fastq` | Path to input FASTQ file (required) |
| `-o, --output` | Path for HTML report (default: `<input_stem>_scout.html`) |
| `--json` | Optional path for JSON export (metrics + scout verdict) |
| `-c, --chunk-size` | Reads per chunk (default: `10000`) |
| `--plot-dir` | Save PNG plots to this directory (also embedded in HTML) |

### Exit codes

| Code | Verdict | Meaning |
|------|---------|---------|
| `0` | `PROCEED` | Sample looks good for downstream analysis |
| `1` | `TRIM` | Preprocessing recommended before continuing |
| `2` | `REJECT` | Sample quality is too low — investigate or re-sequence |

## HTML report

The main output is a single HTML file containing:

1. **Pre-flight verdict** — color-coded banner (`PROCEED` / `TRIM` / `REJECT`)
2. **Summary cards** — mean quality, read length, GC content, duplicate rate
3. **Issues** — detected problems
4. **What to do next** — concrete recommendations (e.g. fastp commands)
5. **QC plots** — embedded PNG charts
6. **Raw metrics summary** — JSON block for programmatic use

Open the HTML file in any browser — no server required.

## Example output

Example run on a bacterial isolate FASTQ file (~11.3M reads, 151 bp):

| Metric | Value |
|--------|-------|
| Mean quality (PHRED) | 34.56 |
| Read length | 151 bp (uniform) |
| Mean GC | 48.4% |
| Duplicate rate | 63.98% |
| **Verdict** | **TRIM** |

**Issues:**
- GC distribution looks bimodal
- Elevated duplicate rate (63.98%)

**Recommendations:**
- Screen for contamination with Kraken or sourmash
- Check PCR cycles; consider `fastp --dedup` for WGS

## Example plots

Generated from a real FASTQ run:

### Per-position quality

![Per-position quality](docs/images/per_position_quality.png)

### Read length distribution

![Length distribution](docs/images/length_distribution.png)

### GC content

![GC content](docs/images/gc_content.png)

### Duplicate rate

![Duplicate rate](docs/images/duplicate_rate.png)

## Programmatic usage

```python
from pathlib import Path

from fastq_scout.reader import FastqReader
from fastq_scout.pipeline import Pipeline
from fastq_scout.metrics import PerPositionQuality, LengthDistribution, GCContent, DuplicateRate
from fastq_scout.scout import FastqScout
from fastq_scout.report import HtmlReport, build_plot_paths

reader = FastqReader("sample.fastq", chunk_size=10_000)
pipeline = Pipeline(metrics=[
    PerPositionQuality(),
    LengthDistribution(),
    GCContent(),
    DuplicateRate(),
])

results = pipeline.run(reader)
verdict, scout_report = FastqScout().result(results)

plot_paths = build_plot_paths(results, Path("plots"))
HtmlReport(Path("sample.fastq"), results, scout_report, plot_paths).save(Path("sample_scout.html"))

print(verdict, scout_report)
```

## How it works

```text
FASTQ file
    │
    ▼
FastqReader (chunks)
    │
    ▼
Pipeline → PerPositionQuality, LengthDistribution, GCContent, DuplicateRate
    │
    ▼
FastqScout → verdict + issues + recommendations
    │
    ▼
HtmlReport → single HTML file with embedded plots
```

## Roadmap

- [ ] Sampling mode (`--max-reads`) for faster pre-flight on huge files
- [ ] Paired-end and gzip support
- [ ] Adapter content detection
- [ ] Multi-sample cohort comparison

---

*FastqScout — lightweight pre-flight QC for FASTQ data.*
