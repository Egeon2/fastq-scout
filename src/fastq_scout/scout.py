from abc import ABC, abstractmethod

from fastq_scout.profiles import ScoutProfile, get_scout_profile
from fastq_scout.run_context import RunContext


class Scout(ABC):

    @property
    @abstractmethod
    def name(self):
        raise NotImplementedError("Subclass must implement name()")

    @abstractmethod
    def result(self, metrics: dict) -> tuple[str, dict]:
        raise NotImplementedError("Subclass must implement result()")


class FastqScout(Scout):

    def __init__(self, ctx: RunContext | None = None):
        self.ctx = ctx

    @property
    def name(self):
        return "Pre-flight scout"

    def _profile(self) -> ScoutProfile:
        if self.ctx is None:
            return get_scout_profile(RunContext(layout="single", library_type="genome"))
        return get_scout_profile(self.ctx)

    def _evaluate(
        self,
        metrics: dict,
        profile: ScoutProfile,
        read_label: str = "",
    ) -> tuple[list[str], list[str], str]:
        prefix = f"[{read_label}] " if read_label else ""
        issues: list[str] = []
        recommendations: list[str] = []

        quality = metrics.get("Per position quality", {})
        length = metrics.get("Length distribution", {})
        gc = metrics.get("GC content", {})
        duplicates = metrics.get("Duplicates rate", 0)
        base_content = metrics.get("Per base sequence content", {})
        base_summary = base_content.get("summary", {})

        overall_mean = quality.get("overall_mean", 0)
        per_position = quality.get("per_position_mean", [])

        if overall_mean < 20:
            issues.append(f"{prefix}Overall mean quality is too low ({overall_mean})")
            recommendations.append("Reject sample or repeat sequencing")
        elif overall_mean < 25:
            issues.append(
                f"{prefix}Overall mean quality is below recommended threshold ({overall_mean})"
            )
            recommendations.append("Run fastp --cut_mean_quality 20 --cut_tail")

        if len(per_position) >= 10:
            head_mean = sum(per_position[:10]) / 10
            tail_start = int(len(per_position) * 0.8)
            tail_mean = sum(per_position[tail_start:]) / len(per_position[tail_start:])
            if head_mean - tail_mean > 5:
                issues.append(
                    f"{prefix}Quality drops at read tail "
                    f"(head {round(head_mean, 1)} vs tail {round(tail_mean, 1)})"
                )
                recommendations.append(
                    "Run fastp --cut_tail --cut_window_size 4 --cut_mean_quality 20"
                )

        min_length = length.get("min_length")
        mean_length = length.get("mean_length")
        if min_length is not None and min_length < profile.min_length_warn:
            issues.append(f"{prefix}Very short reads detected (min length {min_length})")
            recommendations.append(
                f"Filter short reads: fastp --length_required {profile.length_required}"
            )
        if mean_length is not None and mean_length < profile.mean_length_warn:
            issues.append(f"{prefix}Mean read length is low ({mean_length})")
            recommendations.append("Check library prep and adapter trimming settings")

        mean_gc = gc.get("mean_gc", 0)
        if mean_gc < 30 or mean_gc > 70:
            issues.append(f"{prefix}Unusual mean GC content ({mean_gc}%)")
            recommendations.append("Check for contamination or mixed sample")

        gc_histogram = gc.get("gc_histogram", [])
        if gc_histogram:
            peak_count = sum(
                1 for count in gc_histogram if count > sum(gc_histogram) * 0.05
            )
            if peak_count >= 2:
                issues.append(f"{prefix}GC distribution looks bimodal")
                recommendations.append(profile.gc_bimodal_note)

        if duplicates > profile.duplicate_reject:
            issues.append(f"{prefix}Very high duplicate rate ({duplicates}%)")
            recommendations.append(f"Consider deduplication: fastp --dedup")
        elif duplicates > profile.duplicate_warn:
            issues.append(f"{prefix}Elevated duplicate rate ({duplicates}%)")
            recommendations.append(profile.dedup_recommendation)

        mean_n_pct = base_summary.get("mean_n_pct", 0)
        max_n_pct = base_summary.get("max_n_pct", 0)
        if max_n_pct > 5:
            issues.append(
                f"{prefix}High unknown-base (N) content at some positions "
                f"(up to {max_n_pct}%)"
            )
            recommendations.append(
                "Filter or trim reads with N: fastp --n_base_limit 5"
            )
        elif mean_n_pct > 1:
            issues.append(
                f"{prefix}Elevated unknown-base (N) content ({mean_n_pct}% on average)"
            )

        mean_iupac_pct = base_summary.get("mean_iupac_pct", 0)
        if mean_iupac_pct > 0.5:
            issues.append(
                f"{prefix}IUPAC ambiguous bases detected ({mean_iupac_pct}% on average)"
            )
            recommendations.append(
                "Check base-calling quality or convert/filter ambiguous bases before alignment"
            )

        mean_u_pct = base_summary.get("mean_u_pct", 0)
        if mean_u_pct > 0.5:
            issues.append(
                f"{prefix}RNA uracil (U) bases detected ({mean_u_pct}% on average)"
            )
            recommendations.append(
                "Confirm library type; for DNA data consider U→T normalization before mapping"
            )

        adapter = metrics.get("Adapter discovery", {})
        adapter_pct = adapter.get("adapter_content_pct", 0)
        trim_sequence = adapter.get("trim_sequence") or adapter.get("consensus", "")
        reference_name = adapter.get("reference_name", "")
        if adapter_pct > 5:
            issues.append(f"{prefix}Adapter sequence detected on read tails ({adapter_pct}%)")
            if trim_sequence:
                label = reference_name or "detected adapter"
                fastp_flag = (
                    f"--adapter_sequence {trim_sequence}"
                    if read_label != "R2"
                    else f"--adapter_sequence_r2 {trim_sequence}"
                )
                recommendations.append(f"Run fastp {fastp_flag}  # {label}")
            else:
                recommendations.append("Run fastp with adapter trimming enabled")
        elif adapter_pct > 0.5 and trim_sequence:
            label = reference_name or "detected adapter"
            fastp_flag = (
                f"--adapter_sequence {trim_sequence}"
                if read_label != "R2"
                else f"--adapter_sequence_r2 {trim_sequence}"
            )
            issues.append(f"{prefix}Low-level adapter signal detected ({adapter_pct}%)")
            recommendations.append(f"Consider trimming: fastp {fastp_flag}  # {label}")

        if overall_mean < 20:
            partial = "REJECT"
        elif duplicates > profile.duplicate_reject and profile.reject_on_high_duplicates:
            partial = "REJECT"
        elif issues:
            partial = "TRIM"
        else:
            partial = "PROCEED"

        return issues, recommendations, partial

    def _merge_verdicts(self, *verdicts: str) -> str:
        if "REJECT" in verdicts:
            return "REJECT"
        if "TRIM" in verdicts:
            return "TRIM"
        return "PROCEED"

    def result(self, metrics: dict) -> tuple[str, dict]:
        profile = self._profile()
        issues, recommendations, partial = self._evaluate(metrics, profile)
        verdict = partial

        report = {
            "verdict": verdict,
            "issues": issues,
            "recommendations": _dedupe_preserve_order(recommendations),
            "library_type": self.ctx.library_type if self.ctx else "genome",
            "layout": self.ctx.layout if self.ctx else "single",
        }
        return verdict, report

    def result_paired(
        self,
        metrics_r1: dict,
        metrics_r2: dict,
    ) -> tuple[str, dict]:
        profile = self._profile()
        issues_r1, recs_r1, verdict_r1 = self._evaluate(metrics_r1, profile, "R1")
        issues_r2, recs_r2, verdict_r2 = self._evaluate(metrics_r2, profile, "R2")
        verdict = self._merge_verdicts(verdict_r1, verdict_r2)

        report = {
            "verdict": verdict,
            "issues": issues_r1 + issues_r2,
            "recommendations": _dedupe_preserve_order(recs_r1 + recs_r2),
            "library_type": self.ctx.library_type if self.ctx else "genome",
            "layout": "paired",
            "R1": {"verdict": verdict_r1},
            "R2": {"verdict": verdict_r2},
        }
        return verdict, report


def _dedupe_preserve_order(items: list[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result
