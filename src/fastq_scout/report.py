import base64
import html
import json
from datetime import datetime, timezone
from pathlib import Path


VERDICT_STYLES = {
    "PROCEED": ("#1e8449", "#eafaf1", "Sample looks good — you can proceed with downstream analysis."),
    "TRIM": ("#b7950b", "#fef9e7", "Sample needs preprocessing before downstream analysis."),
    "REJECT": ("#c0392b", "#fdedec", "Sample quality is too low — do not proceed without re-sequencing or investigation."),
}


class HtmlReport:

    def __init__(self, fastq_path: Path, metrics: dict, scout_report: dict, plot_paths: dict[str, Path]):
        self.fastq_path = fastq_path
        self.metrics = metrics
        self.scout_report = scout_report
        self.plot_paths = plot_paths

    def save(self, output_path: Path) -> Path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(self.render(), encoding="utf-8")
        return output_path

    def render(self) -> str:
        verdict = self.scout_report.get("verdict", "UNKNOWN")
        accent, background, verdict_hint = VERDICT_STYLES.get(
            verdict, ("#566573", "#f4f6f7", "Review the report before continuing.")
        )

        quality = self.metrics.get("Per position quality", {})
        length = self.metrics.get("Length distribution", {})
        gc = self.metrics.get("GC content", {})
        duplicates = self.metrics.get("Duplicates rate", 0)
        mean_gc = gc.get("mean_gc", 0)

        issues = self.scout_report.get("issues", [])
        recommendations = self.scout_report.get("recommendations", [])

        generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>FastqScout Report — {html.escape(self.fastq_path.name)}</title>
    <style>
        :root {{
            --accent: {accent};
            --accent-bg: {background};
            --text: #1c2833;
            --muted: #5d6d7e;
            --border: #d5dbdb;
            --card: #ffffff;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            margin: 0;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
            background: #f8f9fa;
            color: var(--text);
            line-height: 1.5;
        }}
        .container {{
            max-width: 1100px;
            margin: 0 auto;
            padding: 32px 20px 48px;
        }}
        header {{
            margin-bottom: 24px;
        }}
        header h1 {{
            margin: 0 0 8px;
            font-size: 1.9rem;
        }}
        header p {{
            margin: 0;
            color: var(--muted);
        }}
        .verdict {{
            background: var(--accent-bg);
            border: 1px solid var(--accent);
            border-left: 6px solid var(--accent);
            border-radius: 10px;
            padding: 20px 24px;
            margin-bottom: 24px;
        }}
        .verdict-label {{
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            color: var(--muted);
            margin-bottom: 6px;
        }}
        .verdict-value {{
            font-size: 2rem;
            font-weight: 700;
            color: var(--accent);
            margin-bottom: 8px;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .card {{
            background: var(--card);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 16px;
        }}
        .card-label {{
            font-size: 0.85rem;
            color: var(--muted);
            margin-bottom: 6px;
        }}
        .card-value {{
            font-size: 1.5rem;
            font-weight: 700;
        }}
        section {{
            background: var(--card);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 20px 24px;
            margin-bottom: 24px;
        }}
        section h2 {{
            margin: 0 0 16px;
            font-size: 1.2rem;
        }}
        ul {{
            margin: 0;
            padding-left: 20px;
        }}
        li {{
            margin-bottom: 8px;
        }}
        .plots {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 20px;
        }}
        .plot-card img {{
            width: 100%;
            height: auto;
            border: 1px solid var(--border);
            border-radius: 8px;
            background: #fff;
        }}
        .plot-card h3 {{
            margin: 0 0 12px;
            font-size: 1rem;
        }}
        .empty {{
            color: var(--muted);
            font-style: italic;
        }}
        footer {{
            color: var(--muted);
            font-size: 0.85rem;
            text-align: center;
            margin-top: 16px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>FastqScout Report</h1>
            <p>Input: <strong>{html.escape(str(self.fastq_path))}</strong></p>
            <p>Generated: {generated_at}</p>
        </header>

        <div class="verdict">
            <div class="verdict-label">Pre-flight verdict</div>
            <div class="verdict-value">{html.escape(verdict)}</div>
            <p>{html.escape(verdict_hint)}</p>
        </div>

        <div class="grid">
            <div class="card">
                <div class="card-label">Mean quality (PHRED)</div>
                <div class="card-value">{quality.get("overall_mean", "—")}</div>
            </div>
            <div class="card">
                <div class="card-label">Mean read length</div>
                <div class="card-value">{length.get("mean_length", "—")}</div>
            </div>
            <div class="card">
                <div class="card-label">Read length range</div>
                <div class="card-value">{self._length_range(length)}</div>
            </div>
            <div class="card">
                <div class="card-label">Mean GC content</div>
                <div class="card-value">{self._format_gc(mean_gc)}</div>
            </div>
            <div class="card">
                <div class="card-label">Duplicate rate</div>
                <div class="card-value">{self._format_duplicate(duplicates)}</div>
            </div>
        </div>

        <section>
            <h2>Issues</h2>
            {self._list_block(issues, "No issues detected.")}
        </section>

        <section>
            <h2>What to do next</h2>
            {self._list_block(recommendations, "No extra steps required — proceed with your pipeline.")}
        </section>

        <section>
            <h2>QC plots</h2>
            <div class="plots">
                {self._plots_block()}
            </div>
        </section>

        <section>
            <h2>Raw metrics summary</h2>
            <pre>{html.escape(json.dumps(self._metrics_summary(), indent=2))}</pre>
        </section>

        <footer>FastqScout — lightweight pre-flight QC for FASTQ data</footer>
    </div>
</body>
</html>"""

    def _length_range(self, length: dict) -> str:
        min_length = length.get("min_length")
        max_length = length.get("max_length")
        if min_length is None or max_length is None:
            return "—"
        return f"{min_length}–{max_length} bp"

    def _format_gc(self, mean_gc) -> str:
        if mean_gc == 0 and not self.metrics.get("GC content"):
            return "—"
        return f"{mean_gc}%"

    def _format_duplicate(self, duplicates) -> str:
        if duplicates == 0 and "Duplicates rate" not in self.metrics:
            return "—"
        return f"{duplicates}%"

    def _list_block(self, items: list[str], empty_text: str) -> str:
        if not items:
            return f'<p class="empty">{html.escape(empty_text)}</p>'
        rows = "".join(f"<li>{html.escape(item)}</li>" for item in items)
        return f"<ul>{rows}</ul>"

    def _plots_block(self) -> str:
        if not self.plot_paths:
            return '<p class="empty">No plots available.</p>'

        blocks = []
        for title, path in self.plot_paths.items():
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            blocks.append(
                f"""<div class="plot-card">
                    <h3>{html.escape(title)}</h3>
                    <img src="data:image/png;base64,{encoded}" alt="{html.escape(title)}">
                </div>"""
            )
        return "\n".join(blocks)

    def _metrics_summary(self) -> dict:
        summary = {}
        for name, data in self.metrics.items():
            if name == "Length distribution" and isinstance(data, dict):
                summary[name] = {
                    key: value for key, value in data.items() if key != "distribution"
                }
            else:
                summary[name] = data
        summary["scout"] = self.scout_report
        return summary


def build_plot_paths(metrics: dict, output_dir: Path) -> dict[str, Path]:
    from fastq_scout.plot import MetricPlotter

    plot_titles = {
        "Per position quality": "Per position quality",
        "Length distribution": "Length distribution",
        "GC content": "GC content",
        "Duplicates rate": "Duplicate rate",
    }

    plot_paths = {}
    for metric_name, data in metrics.items():
        plot_path = MetricPlotter(metric_name, data).plot(output_dir)
        plot_paths[plot_titles.get(metric_name, metric_name)] = plot_path

    return plot_paths
