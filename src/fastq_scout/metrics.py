from abc import ABC, abstractmethod

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
        self._gc_sum = 0.0

    @property
    def name(self):
        return "GC content"

    def update(self, chunk):
        for read in chunk:
            bin_idx = round(read.gc_content * 100)
            self._gc_bins[bin_idx] += 1
            self._total_reads += 1
            self._gc_sum += read.gc_content

    def result(self):
        mean_gc = round(self._gc_sum / self._total_reads * 100, 2) if self._total_reads else 0.0
        return {
            "mean_gc": mean_gc,
            "gc_histogram": self._gc_bins
        }

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