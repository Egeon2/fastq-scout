import math
from pathlib import Path

import matplotlib.pyplot as plt


def _smooth(values: list[int | float], window: int = 2) -> list[float]:
    smoothed = []
    for i in range(len(values)):
        start = max(0, i - window)
        end = min(len(values), i + window + 1)
        smoothed.append(sum(values[start:end]) / (end - start))
    return smoothed


def _theoretical_gc_curve(mean_gc: float, std_gc: float, total_reads: int, bins: list[int]) -> list[float]:
    if std_gc <= 0 or total_reads == 0:
        return [0.0] * len(bins)

    curve = []
    for x in bins:
        exponent = -0.5 * ((x - mean_gc) / std_gc) ** 2
        pdf = math.exp(exponent) / (std_gc * math.sqrt(2 * math.pi))
        curve.append(pdf * total_reads)
    return curve


class MetricPlotter:
    def __init__(self, metric_name, data):
        self.metric_name = metric_name
        self.data = data

    def plot(self, output_dir: Path = Path(".")) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(10, 6))

        if self.metric_name == "Per position quality":
            plt.plot(self.data["per_position_mean"], color='#2ecc71', linewidth=2)
            plt.title("Per Position Mean Quality", fontsize=14)
            plt.xlabel("Position in Read")
            plt.ylabel("PHRED Score")
            plt.grid(True, linestyle='--', alpha=0.7)
            output_path = output_dir / "per_position_quality.png"

        elif self.metric_name == "Length distribution":
            dist = self.data["distribution"]
            lengths = sorted(dist.keys())
            counts = [dist[l] for l in lengths]
            plt.bar(lengths, counts, color='#3498db', alpha=0.8)
            plt.title("Read Length Distribution", fontsize=14)
            plt.xlabel("Length (bp)")
            plt.ylabel("Number of Reads")
            output_path = output_dir / "length_distribution.png"

        elif self.metric_name == "GC content":
            bins = list(range(101))
            counts = self.data["gc_histogram"]
            total_reads = sum(counts)
            mean_gc = self.data["mean_gc"]
            std_gc = self.data["std_gc"]
            observed = _smooth(counts)
            theoretical = _theoretical_gc_curve(mean_gc, std_gc, total_reads, bins)

            plt.plot(bins, observed, color='#e74c3c', linewidth=2, label='GC distribution')
            plt.plot(bins, theoretical, color='#3498db', linewidth=2, label='Theoretical distribution')
            plt.title("Per Sequence GC Content", fontsize=14)
            plt.xlabel("% GC")
            plt.ylabel("Count")
            plt.legend(loc='upper right')
            plt.grid(True, linestyle='--', alpha=0.7)
            output_path = output_dir / "gc_content.png"

        elif self.metric_name == "Per base sequence content":
            positions = list(range(1, len(self.data["A"]) + 1))
            plt.plot(positions, self.data["A"], color='#27ae60', linewidth=1.5, label='A')
            plt.plot(positions, self.data["C"], color='#3498db', linewidth=1.5, label='C')
            plt.plot(positions, self.data["T"], color='#e74c3c', linewidth=1.5, label='T')
            plt.plot(positions, self.data["G"], color='#2c3e50', linewidth=1.5, label='G')
            plt.title("Per Base Sequence Content", fontsize=14)
            plt.xlabel("Position in Read")
            plt.ylabel("% Sequence Content")
            plt.legend(loc='upper right')
            plt.grid(True, linestyle='--', alpha=0.7)
            output_path = output_dir / "per_base_sequence_content.png"

        elif self.metric_name == "Duplicates rate":
            labels = ['Duplicates', 'Unique']
            vals = [self.data, 100 - self.data]
            plt.pie(vals, labels=labels, autopct='%1.1f%%', colors=['#e74c3c', '#95a5a6'], startangle=140)
            plt.title("Duplicate Rate Overview", fontsize=14)
            output_path = output_dir / "duplicate_rate.png"
        
        elif self.metric_name == "Per sequence quality scores":
            bins = list(range(41))
            counts = self.data["histogram"]
            plt.bar(bins, counts, color='#27ae60', alpha=0.85, width=1.0)
            plt.title("Per Sequence Quality Scores", fontsize=14)
            plt.xlabel("Mean PHRED Score")
            plt.ylabel("Number of Reads")
            plt.grid(True, linestyle='--', alpha=0.7)
            output_path = output_dir / "per_sequence_quality_scores.png"

        elif self.metric_name == "Adapter discovery":
            enrichment_top = self.data.get("enrichment_top", [])
            method = self.data.get("detection_method", "none")
            if not enrichment_top:
                plt.text(0.5, 0.5, "No adapter candidates", ha="center", va="center")
                plt.axis("off")
            elif method == "reference":
                labels = [item["kmer"] for item in enrichment_top]
                values = [item["enrichment"] for item in enrichment_top]
                plt.barh(labels[::-1], values[::-1], color="#8e44ad", alpha=0.85)
                plt.title("Reference Adapter Match Rate", fontsize=14)
                plt.xlabel("Reads matched (%)")
                plt.ylabel("Reference")
                plt.grid(True, linestyle='--', alpha=0.7)
            else:
                kmers = [item["kmer"] for item in enrichment_top]
                values = [item["enrichment"] for item in enrichment_top]
                plt.barh(kmers[::-1], values[::-1], color="#9b59b6", alpha=0.85)
                plt.title("Tail k-mer Enrichment (de novo fallback)", fontsize=14)
                plt.xlabel("Enrichment (tail / middle)")
                plt.ylabel("k-mer")
                plt.grid(True, linestyle='--', alpha=0.7)
            output_path = output_dir / "adapter_enrichment.png"

        else:
            plt.close()
            raise ValueError(f"Unknown metric: {self.metric_name}")

        plt.savefig(output_path)
        plt.close()
        return output_path
