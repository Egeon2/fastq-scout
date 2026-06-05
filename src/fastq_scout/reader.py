from pathlib import Path
from typing import Iterator

from fastq_scout.models import SequenceRecord


class FastqReader:
    def __init__(self, file: Path, chunk_size: int, sample_budget: int | None = None):
        self.file = file
        self.chunk_size = chunk_size
        self.sample_budget = sample_budget
        self.reads_processed = 0

    def __iter__(self) -> Iterator[list[SequenceRecord]]:
        chunk = []
        with open(self.file, encoding="utf-8", errors="replace") as f:
            while True:
                if self.sample_budget is not None and self.reads_processed >= self.sample_budget:
                    break

                header = f.readline().strip().lstrip("@")
                if not header:
                    break
                sequence = f.readline().strip()
                _ = f.readline()
                quality = f.readline().strip()

                chunk.append(SequenceRecord(header, sequence, quality))
                self.reads_processed += 1

                if len(chunk) >= self.chunk_size:
                    yield chunk
                    chunk = []

                if self.sample_budget is not None and self.reads_processed >= self.sample_budget:
                    if chunk:
                        yield chunk
                    return

        if chunk:
            yield chunk
