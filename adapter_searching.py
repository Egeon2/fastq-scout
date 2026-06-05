from pathlib import Path

def raw_sequences(file: Path) -> str:
    with open(file, "r") as f:
        for line in f:
            if line.startswith(">"):
                yield line.strip()[1:]

def kmer_search(seq: str, kmer: int) -> dict[str, int]:
    kmer_dict = {}
    for i in range(len(seq) - kmer + 1):
        kmer = seq[i:i+kmer]
        if kmer in kmer_dict:
            kmer_dict[kmer] += 1
        else:
            kmer_dict[kmer] = 1
    return kmer_dict