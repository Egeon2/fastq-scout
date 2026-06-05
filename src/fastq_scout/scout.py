from abc import ABC, abstractmethod


class Scout(ABC):

    @property
    @abstractmethod
    def name(self):
        raise NotImplementedError("Subclass must implement name()")

    @abstractmethod
    def result(self, metrics: dict) -> tuple[str, dict]:
        raise NotImplementedError("Subclass must implement result()")


class FastqScout(Scout):

    @property
    def name(self):
        return "Pre-flight scout"

    def result(self, metrics: dict) -> tuple[str, dict]:
        issues = []
        recommendations = []

        quality = metrics.get("Per position quality", {})
        length = metrics.get("Length distribution", {})
        gc = metrics.get("GC content", {})
        duplicates = metrics.get("Duplicates rate", 0)

        overall_mean = quality.get("overall_mean", 0)
        per_position = quality.get("per_position_mean", [])

        if overall_mean < 20:
            issues.append(f"Overall mean quality is too low ({overall_mean})")
            recommendations.append("Reject sample or repeat sequencing")
        elif overall_mean < 25:
            issues.append(f"Overall mean quality is below recommended threshold ({overall_mean})")
            recommendations.append("Run fastp --cut_mean_quality 20 --cut_tail")

        if len(per_position) >= 10:
            head_mean = sum(per_position[:10]) / 10
            tail_start = int(len(per_position) * 0.8)
            tail_mean = sum(per_position[tail_start:]) / len(per_position[tail_start:])
            if head_mean - tail_mean > 5:
                issues.append(
                    f"Quality drops at read tail (head {round(head_mean, 1)} vs tail {round(tail_mean, 1)})"
                )
                recommendations.append("Run fastp --cut_tail --cut_window_size 4 --cut_mean_quality 20")

        min_length = length.get("min_length")
        mean_length = length.get("mean_length")
        if min_length is not None and min_length < 20:
            issues.append(f"Very short reads detected (min length {min_length})")
            recommendations.append("Filter short reads: fastp --length_required 50")
        if mean_length is not None and mean_length < 50:
            issues.append(f"Mean read length is low ({mean_length})")
            recommendations.append("Check library prep and adapter trimming settings")

        mean_gc = gc.get("mean_gc", 0)
        if mean_gc < 30 or mean_gc > 70:
            issues.append(f"Unusual mean GC content ({mean_gc}%)")
            recommendations.append("Check for contamination or mixed sample")

        gc_histogram = gc.get("gc_histogram", [])
        if gc_histogram:
            peak_count = sum(1 for count in gc_histogram if count > sum(gc_histogram) * 0.05)
            if peak_count >= 2:
                issues.append("GC distribution looks bimodal")
                recommendations.append("Screen for contamination with Kraken or sourmash")

        if duplicates > 70:
            issues.append(f"Very high duplicate rate ({duplicates}%)")
            recommendations.append("Consider deduplication: fastp --dedup")
        elif duplicates > 40:
            issues.append(f"Elevated duplicate rate ({duplicates}%)")
            recommendations.append("Check PCR cycles; consider fastp --dedup for WGS")

        adapter = metrics.get("Adapter discovery", {})
        adapter_pct = adapter.get("adapter_content_pct", 0)
        consensus = adapter.get("consensus", "")
        if adapter_pct > 5:
            issues.append(f"Adapter sequence detected on read tails ({adapter_pct}%)")
            if consensus:
                recommendations.append(
                    f"Run fastp --adapter_sequence {consensus} --cut_tail"
                )
            else:
                recommendations.append("Run fastp with adapter trimming enabled")
        elif adapter_pct > 0.5 and consensus:
            issues.append(f"Low-level adapter signal detected ({adapter_pct}%)")
            recommendations.append(
                f"Consider trimming: fastp --adapter_sequence {consensus}"
            )

        if overall_mean < 20 or duplicates > 70:
            verdict = "REJECT"
        elif issues:
            verdict = "TRIM"
        else:
            verdict = "PROCEED"

        report = {
            "verdict": verdict,
            "issues": issues,
            "recommendations": recommendations,
        }

        return verdict, report
