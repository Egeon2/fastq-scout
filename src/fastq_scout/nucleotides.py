"""IUPAC nucleotide classification for sequence QC."""

STANDARD_DNA = frozenset("ACGT")
RNA_BASE = "U"
UNKNOWN_BASES = frozenset("NX.")

IUPAC_AMBIGUOUS = frozenset("MRWSKYVHDB")

# Plot / metric column order: standard DNA, RNA, unknown, then IUPAC codes.
TRACKED_BASES = ("A", "C", "G", "T", "U", "N", *sorted(IUPAC_AMBIGUOUS), "Other")

IUPAC_MEANINGS = {
    "M": "A or C",
    "R": "A or G",
    "W": "A or T",
    "S": "C or G",
    "K": "G or T",
    "Y": "C or T",
    "V": "A or C or G",
    "H": "A or C or T",
    "D": "A or G or T",
    "B": "C or G or T",
}


def classify_base(base: str) -> str:
    """Map a single base symbol to a tracked category."""
    base = base.upper()
    if base in STANDARD_DNA:
        return base
    if base == RNA_BASE:
        return RNA_BASE
    if base in UNKNOWN_BASES:
        return "N"
    if base in IUPAC_AMBIGUOUS:
        return base
    return "Other"


def is_unambiguous_dna(base: str) -> bool:
    return base.upper() in STANDARD_DNA
