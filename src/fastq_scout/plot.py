import matplotlib.pyplot as plt

class MetricPlotter:
    def __init__(self, metric_name, data):
        self.metric_name = metric_name
        self.data = data

    def plot(self):
        plt.figure(figsize=(10, 6))
        
        if self.metric_name == "Per position quality":
            plt.plot(self.data["pre_position_mean"], color='#2ecc71', linewidth=2)
            plt.title("Per Position Mean Quality", fontsize=14)
            plt.xlabel("Position in Read")
            plt.ylabel("PHRED Score")
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.savefig("per_position_quality.png")
            
        elif self.metric_name == "Length distribution":
            dist = self.data["distribution"]
            lengths = sorted(dist.keys())
            counts = [dist[l] for l in lengths]
            plt.bar(lengths, counts, color='#3498db', alpha=0.8)
            plt.title("Read Length Distribution", fontsize=14)
            plt.xlabel("Length (bp)")
            plt.ylabel("Number of Reads")
            plt.savefig("length_distribution.png")
            
        elif self.metric_name == "GC content":
            bins = list(range(101))
            counts = self.data["gc_histogram"]
            plt.bar(bins, counts, color='#e67e22', width=1.0)
            plt.title("GC Content Distribution", fontsize=14)
            plt.xlabel("GC %")
            plt.ylabel("Number of Reads")
            plt.savefig("gc_content.png")
            
        elif self.metric_name == "Duplicates rate":
            # For duplicate rate, it's just a single number, maybe a pie chart or just skip
            labels = ['Duplicates', 'Unique']
            vals = [self.data, 100 - self.data]
            plt.pie(vals, labels=labels, autopct='%1.1f%%', colors=['#e74c3c', '#95a5a6'], startangle=140)
            plt.title("Duplicate Rate Overview", fontsize=14)
            plt.savefig("duplicate_rate.png")
            
        plt.close()