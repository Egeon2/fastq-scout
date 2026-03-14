from pathlib import Path
from typing import Iterator
from fastq_scout.models import SequenceRecord

class FastqReader:
    def __init__(self, file: Path, chunk_size: int):
        self.file = file
        self.chunk_size = chunk_size

    def __iter__(self) -> Iterator[list[SequenceRecord]]:
        chunk = []
        with open(self.file, "r") as f:
            while True:
                header = f.readline().strip().lstrip('@')
                if not header:
                    break
                sequence = f.readline().strip()
                _ = f.readline()
                quality = f.readline().strip()

                chunk.append(SequenceRecord(header, sequence, quality))

                if len(chunk) >= self.chunk_size:
                    yield chunk
                    chunk = []

        if chunk:
            yield chunk