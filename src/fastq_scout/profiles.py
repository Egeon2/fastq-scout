from dataclasses import dataclass

from fastq_scout.run_context import RunContext


@dataclass(frozen=True)
class ScoutProfile:
    duplicate_warn: float
    duplicate_reject: float
    mean_length_warn: float
    min_length_warn: int
    length_required: int
    reject_on_high_duplicates: bool
    dedup_recommendation: str
    gc_bimodal_note: str


PROFILES: dict[tuple[str, str], ScoutProfile] = {
    ("genome", "single"): ScoutProfile(
        duplicate_warn=40,
        duplicate_reject=70,
        mean_length_warn=50,
        min_length_warn=20,
        length_required=50,
        reject_on_high_duplicates=True,
        dedup_recommendation="Check PCR cycles; consider fastp --dedup for WGS",
        gc_bimodal_note="Screen for contamination with Kraken or sourmash",
    ),
    ("genome", "paired"): ScoutProfile(
        duplicate_warn=40,
        duplicate_reject=70,
        mean_length_warn=50,
        min_length_warn=20,
        length_required=50,
        reject_on_high_duplicates=True,
        dedup_recommendation="Check PCR cycles; consider fastp --dedup for WGS",
        gc_bimodal_note="Screen for contamination with Kraken or sourmash",
    ),
    ("transcriptome", "single"): ScoutProfile(
        duplicate_warn=60,
        duplicate_reject=85,
        mean_length_warn=30,
        min_length_warn=15,
        length_required=30,
        reject_on_high_duplicates=False,
        dedup_recommendation="High duplicates are common in RNA-seq; check rRNA depletion",
        gc_bimodal_note="Bimodal GC may reflect isoforms or rRNA; check with FastQC",
    ),
    ("transcriptome", "paired"): ScoutProfile(
        duplicate_warn=60,
        duplicate_reject=85,
        mean_length_warn=30,
        min_length_warn=15,
        length_required=30,
        reject_on_high_duplicates=False,
        dedup_recommendation="High duplicates are common in RNA-seq; check rRNA depletion",
        gc_bimodal_note="Bimodal GC may reflect isoforms or rRNA; check with FastQC",
    ),
}


def get_scout_profile(ctx: RunContext) -> ScoutProfile:
    return PROFILES[(ctx.library_type, ctx.layout)]


def effective_min_reads(ctx: RunContext, min_reads: int) -> int:
    if ctx.is_transcriptome:
        return max(min_reads, 50_000)
    return min_reads
