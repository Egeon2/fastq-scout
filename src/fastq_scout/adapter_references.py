"""Known adapter sequences and short trim motifs for fastp."""

from dataclasses import dataclass


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@dataclass(frozen=True)
class AdapterReference:
    name: str
    sequence: str
    trim_motifs: tuple[str, ...]
    kit: str = "Illumina"


def _motifs(*parts: str) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for part in parts:
        for motif in (part, revcomp(part)):
            if motif and motif not in seen:
                seen.add(motif)
                ordered.append(motif)
    return tuple(ordered)


ADAPTER_REFERENCES: tuple[AdapterReference, ...] = (
    AdapterReference(
        name="Illumina TruSeq Universal",
        kit="Illumina",
        sequence="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        trim_motifs=_motifs(
            "AGATCGGAAGAGCACACGTCTGAA",
            "AGATCGGAAGAGCACACGTCTGA",
            "GATCGGAAGAGCACACGTCTGAA",
            "CACACGTCTGAACTCCAGTCAC",
            "CACGTCTGAACTCCAGTCAC",
            "GTCTGAACTCCAGTCAC",
        ),
    ),
    AdapterReference(
        name="Illumina TruSeq Read2",
        kit="Illumina",
        sequence="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        trim_motifs=_motifs(
            "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
            "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC",
            "AGACGTGTGCTCTTCCGATCT",
            "GACGTGTGCTCTTCCGATCT",
        ),
    ),
    AdapterReference(
        name="Illumina Nextera Transposase",
        kit="Illumina",
        sequence="CTGTCTCTTATACACATCT",
        trim_motifs=_motifs(
            "CTGTCTCTTATACACATCT",
            "GTCTCTTATACACATCT",
            "TCTCTTATACACATCT",
        ),
    ),
    AdapterReference(
        name="Illumina Small RNA 3'",
        kit="Illumina",
        sequence="TGGAATTCTCGGGTGCCAAGG",
        trim_motifs=_motifs(
            "TGGAATTCTCGGGTGCCAAGG",
            "GGAATTCTCGGGTGCCAAGG",
            "AATTCTCGGGTGCCAAGG",
        ),
    ),
)

ADAPTER_SETS: dict[str, tuple[AdapterReference, ...]] = {
    "universal": (ADAPTER_REFERENCES[0],),
    "read2": (ADAPTER_REFERENCES[1],),
    "all": ADAPTER_REFERENCES,
}

