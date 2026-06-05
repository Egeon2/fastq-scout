from abc import ABC, abstractmethod

from fastq_scout.adapter_searching import count_kmers_in_seq, discover_adapters


class BaseMetric(ABC):

    @property
    @abstractmethod
    def name(self):
        raise NotImplementedError("Subclass must implement name()")

    @abstractmethod
    def update(self, chunk):
        raise NotImplementedError("Subclass must implement update()")

    @abstractmethod
    def result(self):
        raise NotImplementedError("Subclass must implement result()")


class PerPositionQuality(BaseMetric):
    def __init__(self):
        self._quality_sum = []
        self._count = []
    
    @property
    def name(self):
        return "Per position quality"

    def update(self, chunk):
        for read in chunk:
            scores = read.phred_scores

            if len(scores) > len(self._quality_sum):
                extra = len(scores) - len(self._quality_sum)
                self._quality_sum.extend([0.0] * extra)
                self._count.extend([0] * extra)

            for i, score in enumerate(scores):
                self._quality_sum[i] += score
                self._count[i] += 1

    def result(self):
        mean = []

        for s, c in zip(self._quality_sum, self._count):
            mean.append(round(s / c, 2) if c > 0 else 0.0)

        return {
            "per_position_mean": mean,
            "overall_mean": round(sum(self._quality_sum) / sum(self._count), 2) if sum(self._count) > 0 else 0.0
        }

class PerSequenceQuality(BaseMetric):
    def __init__(self):
        self._quality_bins = [0] * 41
        self._total_reads = 0

    @property
    def name(self):
        return "Per sequence quality scores"

    def update(self, chunk):
        for read in chunk:
            bin_idx = min(max(0, round(read.mean_quality)), 40)
            self._quality_bins[bin_idx] += 1
            self._total_reads += 1

    def result(self):
        if self._total_reads == 0:
            return {"histogram": self._quality_bins, "q20_pct": 0.0, "q30_pct": 0.0}

        reads_q20_plus = sum(self._quality_bins[20:])
        reads_q30_plus = sum(self._quality_bins[30:])

        return {
            "histogram": self._quality_bins,
            "q20_pct": round(reads_q20_plus / self._total_reads * 100, 2),
            "q30_pct": round(reads_q30_plus / self._total_reads * 100, 2),
        }

class LengthDistribution(BaseMetric):
    def __init__(self):
        self._length_counter = {}
        self._total_read = 0

    @property
    def name(self):
        return "Length distribution"

    def update(self, chunk):
        for read in chunk:
            length = read.length
            self._length_counter[length] = self._length_counter.get(length, 0) + 1
            self._total_read += 1

    def result(self):
        if not self._length_counter:
            return {}
        lengths = list(self._length_counter.keys())
        return {
            "min_length": min(lengths),
            "max_length": max(lengths),
            "mean_length": round(sum(length * count for length, count in self._length_counter.items()) / self._total_read, 2),
            "distribution": self._length_counter
        }

class GCContent(BaseMetric):
    def __init__(self):
        self._gc_bins = [0] * 101
        self._total_reads = 0

    @property
    def name(self):
        return "GC content"

    def update(self, chunk):
        for read in chunk:
            bin_idx = round(read.gc_content * 100)
            self._gc_bins[bin_idx] += 1
            self._total_reads += 1

    def result(self):
        if self._total_reads == 0:
            return {
                "mean_gc": 0.0,
                "std_gc": 0.0,
                "gc_histogram": self._gc_bins,
            }

        mean_gc = sum(i * count for i, count in enumerate(self._gc_bins)) / self._total_reads
        variance = sum((i - mean_gc) ** 2 * count for i, count in enumerate(self._gc_bins)) / self._total_reads
        std_gc = variance ** 0.5

        return {
            "mean_gc": round(mean_gc, 2),
            "std_gc": round(std_gc, 2),
            "gc_histogram": self._gc_bins,
        }

class PerBaseSequenceContent(BaseMetric):
    _BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}

    def __init__(self):
        self._base_counts = []
        self._position_count = []

    @property
    def name(self):
        return "Per base sequence content"

    def update(self, chunk):
        for read in chunk:
            seq = read.sequence.upper()

            if len(seq) > len(self._base_counts):
                extra = len(seq) - len(self._base_counts)
                for _ in range(extra):
                    self._base_counts.append([0, 0, 0, 0])
                    self._position_count.append(0)

            for i, base in enumerate(seq):
                base_idx = self._BASE_INDEX.get(base)
                if base_idx is not None:
                    self._base_counts[i][base_idx] += 1
                self._position_count[i] += 1

    def result(self):
        content = {"A": [], "C": [], "G": [], "T": []}

        for counts, total in zip(self._base_counts, self._position_count):
            if total == 0:
                for base in content:
                    content[base].append(0.0)
            else:
                content["A"].append(round(counts[0] / total * 100, 2))
                content["C"].append(round(counts[1] / total * 100, 2))
                content["G"].append(round(counts[2] / total * 100, 2))
                content["T"].append(round(counts[3] / total * 100, 2))

        return content

class DuplicateRate(BaseMetric):
    def __init__(self):
        self._seen = set()
        self._total = 0
        self._duplicates = 0

    @property
    def name(self):
        return "Duplicates rate"
    
    def update(self, chunk):
        for read in chunk:
            self._total += 1
            h = hash(read.sequence[:50])

            if h in self._seen:
                self._duplicates += 1
            else:
                self._seen.add(h)

    def result(self):
        if self._total == 0:
            return 0.0
        
        rate = (self._duplicates / self._total) * 100
        return round(rate, 2)


class AdapterDiscovery(BaseMetric):
    def __init__(
        self,
        k: int = 10,
        tail_len: int = 15,
        middle_start: int = 20,
        middle_end: int = 80,
    ):
        self.k = k
        self.tail_len = tail_len
        self.middle_start = middle_start
        self.middle_end = middle_end
        self._tail_kmers: dict[str, int] = {}
        self._middle_kmers: dict[str, int] = {}
        self._read_tails: list[str] = []

    @property
    def name(self):
        return "Adapter discovery"

    def update(self, chunk):
        for read in chunk:
            seq = read.sequence.upper()
            if len(seq) < self.k:
                continue

            tail = seq[-self.tail_len:] if len(seq) >= self.tail_len else seq
            self._read_tails.append(tail)
            count_kmers_in_seq(tail, self.k, self._tail_kmers)

            if len(seq) >= self.middle_end:
                middle = seq[self.middle_start:self.middle_end]
                count_kmers_in_seq(middle, self.k, self._middle_kmers)

    def result(self):
        return discover_adapters(
            self._tail_kmers,
            self._middle_kmers,
            self._read_tails,
        )