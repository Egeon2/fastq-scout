from fastq_scout.reader import FastqReader
from fastq_scout.metrics import BaseMetric
from typing import List

class Pipeline:
    def __init__(self, metrics: List[BaseMetric]) -> None:
        self.metrics = metrics 

    def run(self, reader: FastqReader) -> dict:
        print(f"Reading file in chunks of {reader.chunk_size} reads...")
        count = 0
        for chunk in reader:
            for metric in self.metrics:
                metric.update(chunk)
            count += len(chunk)
            if count % 100_000 == 0:
                print(f"Processed {count:,} reads...")

        print("Calculation finished. Generating results...")

        results = {}
        for metric in self.metrics:
            results[metric.name] = metric.result()
        
        return results
