# FastqScout

**FastqScout** is a lightweight Python tool for **pre-flight quality control** of FASTQ sequencing data. It analyzes a statistical sample of reads, computes essential QC metrics, detects Illumina adapter signal, and returns an actionable verdict — **`PROCEED`**, **`TRIM`**, or **`REJECT`** — before you spend time on full FastQC / fastp runs.

```
FASTQ  →  FastqScout (sample + verdict)  →  FastQC / fastp / pipeline
```

> **Quick demo:** open [`reports/sample_scout.html`](reports/sample_scout.html) in a browser — a full report from a real 11.3M-read bacterial WGS library (10K-read sample).

---

## Why FastqScout?

| Tool | Role |
|------|------|
| **FastQC** | Full diagnostic QC — thorough, but slow on large files |
| **fastp** | Trimming + filtering — runs on the whole file |
| **FastqScout** | **Pre-flight gate** — fast sample-based decision: proceed, trim first, or reject |

FastqScout answers: *"Is this library worth processing? What should I run next?"*

---

## Features

### Core QC metrics
- **Per-position quality** — PHRED by read position (head vs tail drop)
- **Per-sequence quality** — Q20 / Q30 rates, quality histogram
- **Length distribution** — min / mean / max read length
- **GC content** — histogram + main-peak Gaussian fit
- **Per-base sequence content** — A/C/G/T bias along reads
- **Duplicate rate** — fraction of duplicated sequences

### Statistical sampling
- Auto sample budget from **Cochran / mean formulas** (95% confidence)
- Separate budgets for base QC vs adapter detection (`--mode with_adapter`)
- Caps: `--min-reads`, `--max-fraction`, `--max-reads`, `--full-file`
- Sampling details in HTML report and JSON export

### Adapter detection (reference-first)
- Matches read tails against **known Illumina adapters** (TruSeq Universal, Read2, Nextera, Small RNA)
- Outputs a **short trim sequence** ready for fastp (`--adapter_sequence …`)
- **De novo fallback** if no reference matches (short k-mer motif, ≤20 bp)
- R1 → Universal adapter; R2 → Read2 adapter (paired-end mode)

### Library-aware verdict engine
- **`--library-type genome`** — stricter duplicate thresholds; high duplicates → REJECT
- **`--library-type transcriptome`** — relaxed duplicate rules; higher default sample size
- Context-specific recommendations (WGS dedup vs RNA-seq rRNA note)

### Layouts
- **Single-end** — one FASTQ file
- **Paired-end** — R1 + R2; separate metrics, adapters, and R2 summary section

### Outputs
- **Self-contained HTML report** — verdict, metrics, issues, recommendations, embedded plots
- **JSON export** — metrics, sampling plan, scout verdict (pipeline-friendly)
- **PNG plots** — optional `--plot-dir`
- **Exit codes** — `0` PROCEED / `1` TRIM / `2` REJECT

---

## Project structure

```text
fastq-scout/
├── reports/
│   ├── sample_scout.html      # Example HTML report (open in browser)
│   ├── sample_scout.json      # Example JSON export
│   └── plots/                 # Example PNG plots
├── docs/images/               # Plot previews for README
├── src/fastq_scout/
│   ├── cli.py                 # CLI entry point
│   ├── reader.py              # Streaming FASTQ reader + sample budget
│   ├── pipeline.py            # Metric orchestration
│   ├── metrics.py             # QC metrics + AdapterDiscovery
│   ├── sampling.py            # Statistical sample size
│   ├── adapter_references.py  # Known Illumina adapter DB
│   ├── adapter_searching.py   # Reference matching + de novo fallback
│   ├── scout.py               # Verdict engine (FastqScout)
│   ├── profiles.py            # genome / transcriptome thresholds
│   ├── run_context.py         # Layout + library type context
│   ├── plot.py                # PNG plot generation
│   └── report.py              # HTML report builder
└── README.md
```

---

## Installation

**Requirements:** Python 3.10+, matplotlib

```bash
pip install matplotlib
git clone <repo-url>
cd fastq-scout
export PYTHONPATH=src
```

---

## Usage examples

### 1. Minimal run (single-end, genome)

Generates `<input_stem>_scout.html` next to the FASTQ:

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py your_sample.fastq
```

### 2. Full run with adapter detection (recommended for unknown libraries)

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py your_sample.fastq \
  --mode with_adapter \
  --library-type genome \
  -o reports/sample_scout.html \
  --json reports/sample_scout.json \
  --plot-dir reports/plots
```

### 3. Transcriptome (RNA-seq)

Relaxed duplicate thresholds; minimum sample size ≥ 50K reads:

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py rnaseq_sample.fastq \
  --library-type transcriptome \
  --mode with_adapter \
  -o reports/rnaseq_scout.html
```

### 4. Paired-end WGS

R1 uses Illumina Universal adapter; R2 uses Read2 adapter:

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py sample_R1.fastq sample_R2.fastq \
  --layout paired \
  --library-type genome \
  --mode with_adapter \
  -o reports/pe_scout.html \
  --json reports/pe_scout.json
```

### 5. Large file — control sampling

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py huge_sample.fastq \
  --min-reads 100000 \
  --max-fraction 0.05 \
  --mode with_adapter
```

### 6. Full-file analysis (no sampling)

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py your_sample.fastq --full-file
```

---

## CLI reference

| Option | Default | Description |
|--------|---------|-------------|
| `fastq` | — | R1 or single-end FASTQ (required) |
| `fastq_r2` | — | R2 FASTQ (required with `--layout paired`) |
| `--layout` | `single` | `single` or `paired` |
| `--library-type` | `genome` | `genome` or `transcriptome` |
| `-o, --output` | `<stem>_scout.html` | HTML report path |
| `--json` | — | JSON export path |
| `--plot-dir` | — | Save PNG plots (also embedded in HTML) |
| `-c, --chunk-size` | `10000` | Reads per processing chunk |
| `--mode` | `base` | `base` or `with_adapter` |
| `--min-reads` | `10000` | Minimum reads in auto sample |
| `--max-fraction` | `0.15` | Max fraction of file to analyze |
| `--max-reads` | — | Hard cap on reads analyzed |
| `--full-file` | off | Analyze entire file |
| `--margin-rate` | `0.05` | CI margin for rate metrics |
| `--margin-mean` | `1.0` | CI margin for mean PHRED |

### Exit codes

| Code | Verdict | Meaning |
|------|---------|---------|
| `0` | `PROCEED` | Sample looks good for downstream analysis |
| `1` | `TRIM` | Preprocessing recommended (trimming, filtering) |
| `2` | `REJECT` | Quality too low — investigate or re-sequence |

---

## HTML report

Open [`reports/sample_scout.html`](reports/sample_scout.html) — no server needed.

The report includes:

1. **Run context** — layout (single/paired), library type (genome/transcriptome)
2. **Sampling section** — total reads, sample budget, confidence margins, formula-derived *n*
3. **Verdict banner** — `PROCEED` / `TRIM` / `REJECT` with explanation
4. **Summary cards** — quality, Q20/Q30, GC, duplicates, adapter match
5. **Adapter discovery** — matched reference, fastp trim sequence, identity %
6. **Issues & recommendations** — e.g. `fastp --adapter_sequence GTCTGAACTCCAGTCAC`
7. **QC plots** — embedded PNG (see below)
8. **Raw JSON summary** — for programmatic parsing

For paired-end: separate sampling tables for R1/R2, R2 summary section, R2 plots.

---

## Example results

Real run: **SRR37361581** bacterial isolate, **11,387,076 reads**, 151 bp, **10,000-read sample** (0.09%), `--mode with_adapter`:

| Metric | Value |
|--------|-------|
| Mean quality (PHRED) | 34.1 |
| Q20 / Q30 | 95.3% / 87.0% |
| Mean GC | 48.5% |
| Duplicate rate | 5.7% |
| Adapter content | 13.65% |
| Matched reference | Illumina TruSeq Universal |
| fastp trim sequence | `GTCTGAACTCCAGTCAC` |
| **Verdict** | **TRIM** |

**Issues:** quality drop at read tail; bimodal GC; adapter signal on tails

**Recommendations:**
```bash
fastp --adapter_sequence GTCTGAACTCCAGTCAC  # Illumina TruSeq Universal
fastp --cut_tail --cut_window_size 4 --cut_mean_quality 20
```

Full data: [`reports/sample_scout.json`](reports/sample_scout.json)

---

## Plots

Generated plots (also in [`reports/plots/`](reports/plots/)):

### Per-position quality
![Per-position quality](docs/images/per_position_quality.png)

### Per-sequence quality scores
![Per-sequence quality](docs/images/per_sequence_quality_scores.png)

### Read length distribution
![Length distribution](docs/images/length_distribution.png)

### GC content (histogram + main-peak fit)
![GC content](docs/images/gc_content.png)

### Per-base sequence content
![Per-base content](docs/images/per_base_sequence_content.png)

### Duplicate rate
![Duplicate rate](docs/images/duplicate_rate.png)

### Adapter reference match rate
![Adapter enrichment](docs/images/adapter_enrichment.png)

---

## How it works

```text
                    ┌─────────────────────────────────┐
                    │  Statistical sample budget      │
                    │  (Cochran + mean formulas)      │
                    └───────────────┬─────────────────┘
                                    │
FASTQ ──► FastqReader (stream) ──► Pipeline ──► metrics dict
                                    │
                    ┌───────────────┼───────────────┐
                    │               │               │
              PerPositionQ    GCContent      AdapterDiscovery
              LengthDist      Duplicates     (reference-first)
                    │               │               │
                    └───────────────┼───────────────┘
                                    │
                              FastqScout
                         (profiles: genome /
                          transcriptome)
                                    │
                    ┌───────────────┴───────────────┐
                    │  PROCEED / TRIM / REJECT      │
                    │  + fastp recommendations      │
                    └───────────────┬───────────────┘
                                    │
                              HtmlReport + JSON
```

---

## Programmatic usage

```python
from pathlib import Path

from fastq_scout.reader import FastqReader
from fastq_scout.pipeline import Pipeline
from fastq_scout.metrics import (
    PerPositionQuality, PerSequenceQuality, LengthDistribution,
    GCContent, DuplicateRate, AdapterDiscovery,
)
from fastq_scout.run_context import RunContext
from fastq_scout.scout import FastqScout
from fastq_scout.report import HtmlReport, build_plot_paths

reader = FastqReader("sample.fastq", chunk_size=10_000, sample_budget=10_000)
pipeline = Pipeline(metrics=[
    PerPositionQuality(),
    PerSequenceQuality(),
    LengthDistribution(),
    GCContent(),
    DuplicateRate(),
    AdapterDiscovery(adapter_set="universal"),
])

results = pipeline.run(reader)
ctx = RunContext(layout="single", library_type="genome")
verdict, scout_report = FastqScout(ctx).result(results)

plot_paths = build_plot_paths(results, Path("plots"))
HtmlReport(Path("sample.fastq"), results, scout_report, plot_paths).save(
    Path("sample_scout.html")
)
print(verdict, scout_report["recommendations"])
```

---

## Positioning vs existing tools

| | FastqScout | FastQC | fastp |
|---|-----------|--------|-------|
| Speed on 10M+ reads | Fast (sampled) | Slow (full file) | Moderate (full file) |
| Verdict | **Yes** | No | No |
| Adapter → fastp command | **Yes** | Partial (module) | Auto-detect |
| Statistical sample plan | **Yes** | No | No |
| Library-type profiles | **Yes** | No | No |

FastqScout is **complementary** — run it first, then FastQC/fastp on libraries that pass the gate.

---

## Roadmap

- [x] Statistical sampling with confidence intervals
- [x] Adapter detection (reference-first + de novo fallback)
- [x] Single / paired-end layouts
- [x] Genome / transcriptome profiles
- [x] HTML + JSON report
- [ ] Benchmark vs FastQC/fastp on public SRA datasets
- [ ] gzip-compressed FASTQ support
- [ ] Random sampling (currently reads from file start)
- [ ] Multi-sample cohort comparison

---

## Citation & collaboration

This is an active research prototype. If you use FastqScout in a study or are interested in co-developing the statistical pre-flight framework (sampling + verdict policy + benchmark), please open an issue or contact the author.

---

*FastqScout — lightweight pre-flight QC gate for FASTQ sequencing data.*
