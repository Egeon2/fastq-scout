import math
from pathlib import Path

Z_SCORE_95 = 1.96
VALID_MODES = ("base", "with_adapter")


def count_fastq_reads(file: Path) -> int:
    count = 0
    with open(file, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("@"):
                count += 1
    return count


def sample_size_for_proportion(
    margin: float = 0.05,
    p: float = 0.5,
    z: float = Z_SCORE_95,
) -> int:
    if margin <= 0:
        raise ValueError("margin must be positive")
    return math.ceil((z ** 2) * p * (1 - p) / (margin ** 2))


def sample_size_for_mean(
    margin: float = 1.0,
    std: float = 5.0,
    z: float = Z_SCORE_95,
) -> int:
    if margin <= 0:
        raise ValueError("margin must be positive")
    return math.ceil((z * std / margin) ** 2)


def _budget_from_stat(
    total_reads: int,
    n_stat: int,
    min_reads: int,
    max_fraction: float,
) -> int:
    if total_reads <= min_reads:
        return total_reads

    n_required = max(n_stat, min_reads)
    fraction = min(n_required / total_reads, max_fraction)
    sample_budget = min(math.ceil(total_reads * fraction), total_reads)
    sample_budget = max(sample_budget, min(n_required, total_reads))
    return min(sample_budget, total_reads)


def resolve_sample_budget(
    total_reads: int,
    margin_rate: float = 0.05,
    margin_mean: float = 1.0,
    std_mean: float = 5.0,
    proportion_p: float = 0.5,
    min_reads: int = 100_000,
    max_fraction: float = 0.15,
    max_reads: int | None = None,
    mode: str = "base",
) -> dict:
    if mode not in VALID_MODES:
        raise ValueError(f"mode must be one of {VALID_MODES}")

    if total_reads <= 0:
        return {
            "total_reads": 0,
            "sample_budget": 0,
            "sample_budget_for_adapters": 0,
            "sample_fraction_pct": 0.0,
            "sample_fraction_pct_for_adapters": 0.0,
            "confidence_pct": 95.0,
            "margin_rate": margin_rate,
            "margin_mean": margin_mean,
            "mode": mode,
        }

    n_proportion = sample_size_for_proportion(margin=margin_rate, p=proportion_p)
    n_mean = sample_size_for_mean(margin=margin_mean, std=std_mean)
    n_stat_base = max(n_proportion, n_mean)

    n_proportion_adapter = sample_size_for_proportion(margin=0.01, p=0.01)
    n_stat_adapter = max(n_stat_base, n_proportion_adapter)

    sample_budget_base = _budget_from_stat(total_reads, n_stat_base, min_reads, max_fraction)
    sample_budget_for_adapters = _budget_from_stat(
        total_reads, n_stat_adapter, min_reads, max_fraction
    )

    if mode == "with_adapter":
        sample_budget = max(sample_budget_base, sample_budget_for_adapters)
    else:
        sample_budget = sample_budget_base

    if max_reads is not None:
        sample_budget = min(max_reads, total_reads)
        sample_budget_for_adapters = min(max_reads, total_reads)

    sample_fraction_pct = round(sample_budget / total_reads * 100, 2)
    sample_fraction_pct_for_adapters = round(
        sample_budget_for_adapters / total_reads * 100, 2
    )

    return {
        "total_reads": total_reads,
        "sample_budget": sample_budget,
        "sample_budget_base": sample_budget_base,
        "sample_budget_for_adapters": sample_budget_for_adapters,
        "sample_fraction_pct": sample_fraction_pct,
        "sample_fraction_pct_for_adapters": sample_fraction_pct_for_adapters,
        "n_proportion": n_proportion,
        "n_mean": n_mean,
        "n_proportion_adapter": n_proportion_adapter,
        "n_stat_base": n_stat_base,
        "n_stat_adapter": n_stat_adapter,
        "n_required_stats": max(n_stat_base, min_reads),
        "confidence_pct": 95.0,
        "margin_rate": margin_rate,
        "margin_mean": margin_mean,
        "mode": mode,
    }
