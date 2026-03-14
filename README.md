# 🧬 FastqScout

FastqScout is a lightweight, high-performance Python toolkit for Quality Control (QC) analysis of FASTQ sequencing data. It processes large genomics files in chunks to minimize memory footprint while generating essential metrics and visualizations.

## 🚀 Features

- **Memory Efficient**: Processes large FASTQ files (GBs) using stream-based chunking.
- **Comprehensive Metrics**:
  - **Per-position Quality**: Mean PHRED scores across read positions.
  - **GC Content**: Histogram of GC percentage distribution.
  - **Length Distribution**: Analysis of read lengths across the dataset.
  - **Duplicate Rate**: Estimation of sequence duplication levels.
- **Visual Reports**: Automatically generates PNG plots for all analyzed metrics.

## 🛠 Project Structure

```text
fastq-scout/
├── src/
│   └── fastq_scout/
│       ├── cli.py        # Entry point
│       ├── reader.py     # Fastq streaming reader
│       ├── pipeline.py   # Analysis orchestration
│       ├── metrics.py    # Metric calculation logic
│       ├── models.py     # Data structures
│       └── plot.py       # Visualization engine
└── README.md
```

## 📦 Installation

Ensure you have the required dependencies:

```bash
pip install matplotlib
```

## 🖥 Usage

To run the analysis on your FASTQ file, use the following command from the project root:

```bash
export PYTHONPATH=src
python src/fastq_scout/cli.py
```

### Example Code (`cli.py`)

```python
from fastq_scout.reader import FastqReader
from fastq_scout.pipeline import Pipeline
from fastq_scout.metrics import PerPositionQuality, LengthDistribution, GCContent, DuplicateRate
from fastq_scout.plot import MetricPlotter

# Initialize reader and pipeline
reader = FastqReader("data.fastq", chunk_size=10_000)
pipeline = Pipeline(metrics=[
    PerPositionQuality(), 
    LengthDistribution(), 
    GCContent(), 
    DuplicateRate()
])

# Run analysis
results = pipeline.run(reader)

# Generate plots
for metric_name, data in results.items():
    plotter = MetricPlotter(metric_name, data)
    plotter.plot()
```

## 📊 Visualizations

The library generates the following plots:
- `per_position_quality.png`
- `gc_content.png`
- `length_distribution.png`
- `duplicate_rate.png`

---
*Developed for efficient genomics data scouting.*
