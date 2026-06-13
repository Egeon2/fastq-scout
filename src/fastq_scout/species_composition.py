"""CpG O/E composition hints and optional expected-species sanity checks."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SpeciesCpGProfile:
    label: str
    cpg_oe_min: float
    cpg_oe_max: float
    gc_min: float | None = None
    gc_max: float | None = None
    description: str = ""


SPECIES_PROFILES: dict[str, SpeciesCpGProfile] = {
    "human": SpeciesCpGProfile(
        "human",
        0.15,
        0.35,
        38.0,
        48.0,
        "Mammalian genome with CpG methylation suppression",
    ),
    "mouse": SpeciesCpGProfile(
        "mouse",
        0.15,
        0.35,
        38.0,
        45.0,
        "Mammalian genome with CpG methylation suppression",
    ),
    "drosophila": SpeciesCpGProfile(
        "drosophila",
        0.65,
        1.15,
        40.0,
        48.0,
        "Insect genome without strong CpG suppression",
    ),
    "ecoli": SpeciesCpGProfile(
        "ecoli",
        0.85,
        1.15,
        48.0,
        52.0,
        "Prokaryotic genome without CpG methylation suppression",
    ),
}

COMPOSITION_DISCLAIMER = (
    "Composition hint only — not taxonomic identification. "
    "Use Kraken, sourmash, or reference alignment for species confirmation."
)


def count_cpg_dinucleotides(sequence: str) -> int:
    sequence = sequence.upper()
    return sum(1 for i in range(len(sequence) - 1) if sequence[i : i + 2] == "CG")


def compute_cpg_oe(c_count: int, g_count: int, cpg_observed: int, total_bases: int) -> float:
    if total_bases <= 0 or c_count == 0 or g_count == 0:
        return 0.0
    expected = (c_count * g_count) / total_bases
    if expected <= 0:
        return 0.0
    return cpg_observed / expected


def classify_composition(cpg_oe: float) -> str:
    if cpg_oe < 0.45:
        return "mammalian_methylated"
    if cpg_oe < 0.75:
        return "intermediate"
    return "unmethylated_high_oe"


def composition_summary(cpg_oe: float, mean_gc: float | None = None) -> str:
    cls = classify_composition(cpg_oe)
    gc_part = f", mean GC {mean_gc}%" if mean_gc is not None else ""
    if cls == "mammalian_methylated":
        return (
            f"CpG O/E {cpg_oe:.2f}{gc_part} — low CpG depletion typical of "
            "mammalian methylated genomes (human/mouse-like)."
        )
    if cls == "intermediate":
        return (
            f"CpG O/E {cpg_oe:.2f}{gc_part} — intermediate profile; "
            "interpret with library type and expected organism."
        )
    return (
        f"CpG O/E {cpg_oe:.2f}{gc_part} — high CpG O/E typical of "
        "unmethylated or prokaryotic / invertebrate-like DNA."
    )


def assess_expected_species(
    cpg_oe: float,
    mean_gc: float | None,
    expected_species: str,
    *,
    library_type: str = "genome",
) -> dict:
    profile = SPECIES_PROFILES.get(expected_species.lower())
    if profile is None:
        return {"status": "unknown_species", "mismatch": False}

    if library_type == "transcriptome":
        return {
            "status": "skipped_transcriptome",
            "mismatch": False,
            "message": (
                "CpG O/E species check skipped for transcriptome libraries "
                "(read composition reflects transcripts, not whole genome)."
            ),
        }

    cpg_ok = profile.cpg_oe_min <= cpg_oe <= profile.cpg_oe_max
    gc_ok = True
    if mean_gc is not None and profile.gc_min is not None and profile.gc_max is not None:
        gc_ok = profile.gc_min <= mean_gc <= profile.gc_max

    if cpg_ok and gc_ok:
        return {
            "status": "consistent",
            "mismatch": False,
            "message": (
                f"CpG O/E {cpg_oe:.2f} is consistent with expected species "
                f"({profile.label}; typical range {profile.cpg_oe_min}–{profile.cpg_oe_max})."
            ),
        }

    parts = []
    if not cpg_ok:
        parts.append(
            f"CpG O/E {cpg_oe:.2f} outside expected range "
            f"{profile.cpg_oe_min}–{profile.cpg_oe_max} for {profile.label}"
        )
    if not gc_ok and mean_gc is not None:
        parts.append(
            f"mean GC {mean_gc}% outside typical range "
            f"{profile.gc_min}–{profile.gc_max}% for {profile.label}"
        )

    return {
        "status": "mismatch",
        "mismatch": True,
        "issue": (
            f"Composition mismatch for --expected-species {profile.label}: "
            + "; ".join(parts)
        ),
        "message": "; ".join(parts),
        "recommendation": (
            "Verify sample metadata and screen for contamination "
            "(Kraken, sourmash, or align to expected reference)."
        ),
    }
