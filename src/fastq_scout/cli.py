from fastq_scout.reader import FastqReader
from fastq_scout.models import SequenceRecord
from fastq_scout.pipeline import Pipeline
from fastq_scout.metrics import PerPositionQuality, LengthDistribution, GCContent, DuplicateRate
from fastq_scout.plot import MetricPlotter

reader = FastqReader("SRR37361581_PMSS1_single_colony_isolate_1.fastq", chunk_size=10_000)
pipeline = Pipeline(metrics=[PerPositionQuality(), LengthDistribution(), GCContent(), DuplicateRate()])
print("Start processing...")
results = pipeline.run(reader)
print("\n--- Processing Results ---")
for metric_name, data in results.items():
    print(f"\n{metric_name}:")
    print(data)

print("\n--- Plotting Results ---")
for metric_name, data in results.items():
    plotter = MetricPlotter(metric_name, data)
    plotter.plot()

