from fastq_scout.adapter_references import ADAPTER_REFERENCES, AdapterReference

VALID_BASES = set("ACGT")

REFERENCE_MIN_PCT = 0.5
DENOVO_MIN_PCT = 0.5
DENOVO_MAX_TRIM_LEN = 20
MIN_TRIM_LEN = 16


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


def _best_local_identity(tail: str, reference: str) -> float:
    best = 0.0
    max_len = min(len(tail), len(reference))
    for length in range(8, max_len + 1):
        ref_segments = {reference[i:i + length] for i in range(len(reference) - length + 1)}
        for start in range(len(tail) - length + 1):
            segment = tail[start:start + length]
            if segment in ref_segments:
                best = max(best, 100.0)
                continue
            for ref_seg in ref_segments:
                mismatches = _hamming_mismatches(segment, ref_seg)
                if mismatches <= 2:
                    identity = (length - mismatches) / length * 100
                    best = max(best, identity)
    return round(best, 1)


def _reference_identity(read_tails: list[str], reference: str, trim_motif: str) -> float:
    matched_tails = [
        tail for tail in read_tails if fuzzy_suffix_match(tail, trim_motif, max_mismatches=2)
    ]
    if not matched_tails:
        return 0.0
    return round(
        sum(_best_local_identity(tail, reference) for tail in matched_tails) / len(matched_tails),
        1,
    )


def _pick_trim_motif(
    read_tails: list[str],
    adapter: AdapterReference,
    max_mismatches: int = 2,
) -> tuple[str, float]:
    scored = [
        (motif, match_fraction(motif, read_tails, max_mismatches))
        for motif in adapter.trim_motifs
        if len(motif) >= MIN_TRIM_LEN
    ]
    if not scored:
        scored = [
            (motif, match_fraction(motif, read_tails, max_mismatches))
            for motif in adapter.trim_motifs
        ]
    if not scored:
        return "", 0.0

    best_rate = max(rate for _, rate in scored)
    qualifying = [
        (motif, rate) for motif, rate in scored if rate >= best_rate - 0.005
    ]
    qualifying.sort(key=lambda item: (-item[1], len(item[0])))
    return qualifying[0]


def match_adapter_references(
    read_tails: list[str],
    max_mismatches: int = 2,
    references: tuple[AdapterReference, ...] | None = None,
) -> list[dict]:
    candidates = []
    adapter_pool = references if references is not None else ADAPTER_REFERENCES

    for adapter in ADAPTER_REFERENCES:
        trim_motif, match_rate = _pick_trim_motif(read_tails, adapter, max_mismatches)
        if match_rate <= 0:
            continue

        reads_matched = round(match_rate * len(read_tails))
        candidates.append({
            "reference_name": adapter.name,
            "reference_kit": adapter.kit,
            "reference_sequence": adapter.sequence,
            "trim_sequence": trim_motif,
            "sequence": trim_motif,
            "reads_matched": reads_matched,
            "reads_pct": round(match_rate * 100, 2),
            "identity_pct": _reference_identity(read_tails, adapter.sequence, trim_motif),
        })

    candidates.sort(
        key=lambda item: (item["reads_pct"], item["identity_pct"]),
        reverse=True,
    )
    return candidates


def top_kmer_candidates(
    tail_kmers: dict[str, int],
    middle_kmers: dict[str, int],
    top_n: int = 5,
    max_middle_ratio: float = 0.08,
) -> list[dict]:
    scored = []

    for kmer, tail_count in tail_kmers.items():
        if should_filter_kmer(kmer):
            continue
        middle_count = middle_kmers.get(kmer, 0)
        if tail_count > 0 and middle_count / tail_count > max_middle_ratio:
            continue
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


def extend_seed_short(
    seed: str,
    tails: list[str],
    max_extension: int = 8,
    max_mismatches: int = 2,
    max_length: int = DENOVO_MAX_TRIM_LEN,
) -> tuple[str, float]:
    consensus = seed[:max_length]
    best_score = match_fraction(consensus, tails, max_mismatches)

    for _ in range(max_extension):
        improved = False
        for base in "ACGT":
            candidate = (base + consensus)[:max_length]
            if len(candidate) <= len(consensus):
                continue
            score = match_fraction(candidate, tails, max_mismatches)
            if score > best_score + 0.001:
                consensus = candidate
                best_score = score
                improved = True
        for base in "ACGT":
            candidate = (consensus + base)[:max_length]
            if len(candidate) <= len(consensus):
                continue
            score = match_fraction(candidate, tails, max_mismatches)
            if score > best_score + 0.001:
                consensus = candidate
                best_score = score
                improved = True
        if not improved:
            break

    return consensus, best_score


def _discover_denovo(
    tail_kmers: dict[str, int],
    middle_kmers: dict[str, int],
    read_tails: list[str],
    top_n: int = 3,
) -> dict:
    empty = {
        "detection_method": "de_novo",
        "adapter_content_pct": 0.0,
        "consensus": "",
        "trim_sequence": "",
        "reference_name": "",
        "reference_sequence": "",
        "identity_pct": 0.0,
        "candidates": [],
        "enrichment_top": [],
        "reads_analyzed": len(read_tails),
    }

    enrichment_ranked = top_kmer_candidates(tail_kmers, middle_kmers, top_n=10)
    if not enrichment_ranked:
        return empty

    candidates = []
    for item in enrichment_ranked[:top_n]:
        sequence, match_rate = extend_seed_short(item["kmer"], read_tails)
        candidates.append({
            "reference_name": "Unknown (de novo)",
            "reference_kit": "—",
            "reference_sequence": "",
            "trim_sequence": sequence,
            "sequence": sequence,
            "seed": item["kmer"],
            "reads_matched": round(match_rate * len(read_tails)),
            "reads_pct": round(match_rate * 100, 2),
            "tail_count": item["tail_count"],
            "middle_count": item["middle_count"],
            "enrichment": item["enrichment"],
            "identity_pct": 0.0,
        })

    candidates.sort(key=lambda item: item["reads_pct"], reverse=True)
    best = candidates[0]
    if best["reads_pct"] < DENOVO_MIN_PCT:
        return empty

    return {
        "detection_method": "de_novo",
        "adapter_content_pct": best["reads_pct"],
        "consensus": best["trim_sequence"],
        "trim_sequence": best["trim_sequence"],
        "reference_name": best["reference_name"],
        "reference_sequence": "",
        "identity_pct": 0.0,
        "candidates": candidates,
        "enrichment_top": enrichment_ranked,
        "reads_analyzed": len(read_tails),
    }


def _discover_reference(
    read_tails: list[str],
    references: tuple[AdapterReference, ...] | None = None,
) -> dict | None:
    candidates = match_adapter_references(read_tails, references=references)
    if not candidates:
        return None

    best = candidates[0]
    if best["reads_pct"] < REFERENCE_MIN_PCT:
        return None

    return {
        "detection_method": "reference",
        "adapter_content_pct": best["reads_pct"],
        "consensus": best["trim_sequence"],
        "trim_sequence": best["trim_sequence"],
        "reference_name": best["reference_name"],
        "reference_kit": best.get("reference_kit", ""),
        "reference_sequence": best["reference_sequence"],
        "identity_pct": best["identity_pct"],
        "candidates": candidates[:3],
        "enrichment_top": [
            {
                "kmer": item["reference_name"],
                "tail_count": item["reads_matched"],
                "middle_count": 0,
                "enrichment": item["reads_pct"],
                "rank_score": item["reads_pct"],
            }
            for item in candidates[:10]
        ],
        "reads_analyzed": len(read_tails),
    }


def discover_adapters(
    read_tails: list[str],
    tail_kmers: dict[str, int] | None = None,
    middle_kmers: dict[str, int] | None = None,
    references: tuple[AdapterReference, ...] | None = None,
) -> dict:
    empty = {
        "detection_method": "none",
        "adapter_content_pct": 0.0,
        "consensus": "",
        "trim_sequence": "",
        "reference_name": "",
        "reference_sequence": "",
        "identity_pct": 0.0,
        "candidates": [],
        "enrichment_top": [],
        "reads_analyzed": len(read_tails),
    }

    if not read_tails:
        return empty

    reference_result = _discover_reference(read_tails, references=references)
    if reference_result is not None:
        return reference_result

    if tail_kmers is None:
        tail_kmers = {}
    if middle_kmers is None:
        middle_kmers = {}

    denovo_result = _discover_denovo(tail_kmers, middle_kmers, read_tails)
    return denovo_result if denovo_result["adapter_content_pct"] > 0 else empty
