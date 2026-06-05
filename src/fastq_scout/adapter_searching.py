VALID_BASES = set("ACGT")


def is_valid_kmer(kmer: str) -> bool:
    return bool(kmer) and all(base in VALID_BASES for base in kmer)


def is_homopolymer_kmer(kmer: str, min_run: int = 4) -> bool:
    run = 1
    for i in range(1, len(kmer)):
        if kmer[i] == kmer[i - 1]:
            run += 1
            if run >= min_run:
                return True
        else:
            run = 1
    return False


def is_low_complexity(kmer: str, threshold: float = 0.7) -> bool:
    if not kmer:
        return True
    max_base_count = max(kmer.count(base) for base in "ACGT")
    return max_base_count / len(kmer) > threshold


def should_filter_kmer(kmer: str) -> bool:
    return (
        not is_valid_kmer(kmer)
        or is_homopolymer_kmer(kmer)
        or is_low_complexity(kmer)
    )


def count_kmers_in_seq(seq: str, k: int, counter: dict[str, int]) -> None:
    if len(seq) < k:
        return
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if should_filter_kmer(kmer):
            continue
        counter[kmer] = counter.get(kmer, 0) + 1


def _hamming_mismatches(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def fuzzy_suffix_match(tail: str, adapter: str, max_mismatches: int = 1) -> bool:
    """Return True if adapter appears in the read tail (fuzzy Hamming match)."""
    if not tail or not adapter:
        return False
    if len(adapter) > len(tail):
        return _hamming_mismatches(tail, adapter[-len(tail):]) <= max_mismatches

    for start in range(len(tail) - len(adapter) + 1):
        segment = tail[start:start + len(adapter)]
        if _hamming_mismatches(segment, adapter) <= max_mismatches:
            return True
    return False


def match_fraction(
    sequence: str,
    tails: list[str],
    max_mismatches: int = 1,
) -> float:
    if not tails:
        return 0.0
    matches = sum(
        1 for tail in tails if fuzzy_suffix_match(tail, sequence, max_mismatches)
    )
    return matches / len(tails)


def extend_seed(
    seed: str,
    tails: list[str],
    max_extension: int = 20,
    max_mismatches: int = 1,
) -> tuple[str, float]:
    consensus = seed
    best_score = match_fraction(consensus, tails, max_mismatches)

    for _ in range(max_extension):
        improved = False
        for base in "ACGT":
            candidate = base + consensus
            score = match_fraction(candidate, tails, max_mismatches)
            if score > best_score + 0.001:
                consensus = candidate
                best_score = score
                improved = True
        for base in "ACGT":
            candidate = consensus + base
            score = match_fraction(candidate, tails, max_mismatches)
            if score > best_score + 0.001:
                consensus = candidate
                best_score = score
                improved = True
        if not improved:
            break

    return consensus, best_score


def top_kmer_candidates(
    tail_kmers: dict[str, int],
    middle_kmers: dict[str, int],
    top_n: int = 5,
) -> list[dict]:
    scored = []

    for kmer, tail_count in tail_kmers.items():
        if should_filter_kmer(kmer):
            continue
        middle_count = middle_kmers.get(kmer, 0)
        enrichment = (tail_count + 1) / (middle_count + 1)
        scored.append({
            "kmer": kmer,
            "tail_count": tail_count,
            "middle_count": middle_count,
            "enrichment": round(enrichment, 2),
            "rank_score": enrichment * tail_count,
        })

    scored.sort(key=lambda item: item["rank_score"], reverse=True)
    return scored[:top_n]


def discover_adapters(
    tail_kmers: dict[str, int],
    middle_kmers: dict[str, int],
    read_tails: list[str],
    top_n: int = 5,
) -> dict:
    empty = {
        "adapter_content_pct": 0.0,
        "consensus": "",
        "candidates": [],
        "enrichment_top": [],
        "reads_analyzed": len(read_tails),
    }

    if not read_tails:
        return empty

    enrichment_ranked = top_kmer_candidates(tail_kmers, middle_kmers, top_n=10)
    seed_candidates = top_kmer_candidates(tail_kmers, middle_kmers, top_n=top_n)

    extended_candidates = []
    for candidate in seed_candidates:
        sequence, match_rate = extend_seed(candidate["kmer"], read_tails)
        reads_matched = round(match_rate * len(read_tails))
        extended_candidates.append({
            "seed": candidate["kmer"],
            "sequence": sequence,
            "reads_matched": reads_matched,
            "reads_pct": round(match_rate * 100, 2),
            "tail_count": candidate["tail_count"],
            "middle_count": candidate["middle_count"],
            "enrichment": candidate["enrichment"],
        })

    extended_candidates.sort(
        key=lambda item: (item["reads_pct"], item["enrichment"]),
        reverse=True,
    )

    best = extended_candidates[0] if extended_candidates else None
    if best is None:
        return empty

    return {
        "adapter_content_pct": best["reads_pct"],
        "consensus": best["sequence"],
        "candidates": extended_candidates[:3],
        "enrichment_top": enrichment_ranked,
        "reads_analyzed": len(read_tails),
    }
