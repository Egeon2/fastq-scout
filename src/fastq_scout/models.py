class SequenceRecord:
    __slots__ = ("read_id", "sequence", "quality_scores", "_gc_content")

    def __init__(self, read_id, sequence, quality_scores, _gc_content=None):
        self.read_id = read_id
        self.sequence = sequence
        self.quality_scores = quality_scores
        self._gc_content = _gc_content

    @property
    def length(self):
        return len(self.sequence)

    @property
    def phred_scores(self):
        return [ord(symb) - 33 for symb in self.quality_scores]

    @property
    def gc_content(self):
        if self._gc_content is None:
            gc = sum(1 for a in self.sequence if a in "GCgc")
            self._gc_content = gc / len(self.sequence) if len(self.sequence) > 0 else 0
        return self._gc_content

    @property
    def mean_quality(self):
        scores = self.phred_scores
        return sum(scores) / len(scores) if scores else 0