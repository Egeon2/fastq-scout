import base64
import html
from datetime import datetime, timezone
from pathlib import Path


VERDICT_STYLES = {
    "PROCEED": ("#1e8449", "#eafaf1", "Sample looks good — you can proceed with downstream analysis."),
    "TRIM": ("#b7950b", "#fef9e7", "Sample needs preprocessing before downstream analysis."),
    "REJECT": ("#c0392b", "#fdedec", "Sample quality is too low — do not proceed without re-sequencing or investigation."),
}


class HtmlReport:

    def __init__(
        self,
        fastq_path: Path,
        metrics: dict,
        scout_report: dict,
        plot_paths: dict[str, Path],
        sample_plan: dict | None = None,
        fastq_r2: Path | None = None,
        r2_metrics: dict | None = None,
    ):
        self.fastq_path = fastq_path
        self.metrics = metrics
        self.scout_report = scout_report
        self.plot_paths = plot_paths
        self.sample_plan = sample_plan or {}
        self.fastq_r2 = fastq_r2
        self.r2_metrics = r2_metrics or {}

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
        seq_quality = self.metrics.get("Per sequence quality scores", {})
        length = self.metrics.get("Length distribution", {})
        gc = self.metrics.get("GC content", {})
        duplicates = self.metrics.get("Duplicates rate", 0)
        mean_gc = gc.get("mean_gc", 0)
        q20_pct = seq_quality.get("q20_pct", "—")
        q30_pct = seq_quality.get("q30_pct", "—")
        adapter = self.metrics.get("Adapter discovery", {})
        adapter_pct = adapter.get("adapter_content_pct", "—")
        adapter_trim = adapter.get("trim_sequence") or adapter.get("consensus", "")
        adapter_reference = adapter.get("reference_name", "")
        adapter_ref_seq = adapter.get("reference_sequence", "")
        base_content = self.metrics.get("Per base sequence content", {})
        base_summary = base_content.get("summary", {})

        issues = self.scout_report.get("issues", [])
        recommendations = self.scout_report.get("recommendations", [])

        generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        sampling_section = self._sampling_section()
        adapter_section = self._adapter_section(adapter)
        r2_section = self._r2_summary_section()
        run_meta = self._run_meta_lines()

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
        .grid-adapter {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .card {{
            background: var(--card);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 16px;
            min-width: 0;
        }}
        .card--seq {{
            grid-column: 1 / -1;
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
        .card-value--text {{
            font-size: 1rem;
            font-weight: 600;
            line-height: 1.45;
            word-break: break-all;
            overflow-wrap: anywhere;
            white-space: normal;
        }}
        .card-value--mono {{
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
            font-size: 0.92rem;
            font-weight: 500;
            line-height: 1.5;
            word-break: break-all;
            overflow-wrap: anywhere;
            white-space: normal;
            letter-spacing: 0.02em;
        }}
        .adapter-table {{
            width: 100%;
            border-collapse: collapse;
            table-layout: fixed;
        }}
        .adapter-table th,
        .adapter-table td {{
            padding: 10px 12px;
            border-bottom: 1px solid var(--border);
            vertical-align: top;
            text-align: left;
        }}
        .adapter-table th {{
            font-size: 0.85rem;
            color: var(--muted);
        }}
        .adapter-table .seq-cell {{
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
            font-size: 0.88rem;
            word-break: break-all;
            overflow-wrap: anywhere;
            white-space: normal;
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
            grid-template-columns: 1fr;
            gap: 24px;
        }}
        .plot-card {{
            cursor: zoom-in;
        }}
        .plot-card a {{
            display: block;
            text-decoration: none;
        }}
        .plot-card img {{
            width: 100%;
            height: auto;
            min-height: 280px;
            max-height: 520px;
            object-fit: contain;
            border: 1px solid var(--border);
            border-radius: 8px;
            background: #fff;
            transition: box-shadow 0.15s ease, transform 0.15s ease;
        }}
        .plot-card img:hover {{
            box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12);
            transform: translateY(-1px);
        }}
        .plot-card h3 {{
            margin: 0 0 12px;
            font-size: 1.05rem;
        }}
        .plot-hint {{
            margin: 0 0 16px;
            color: var(--muted);
            font-size: 0.9rem;
        }}
        .plot-modal {{
            display: none;
            position: fixed;
            inset: 0;
            z-index: 1000;
            background: rgba(0, 0, 0, 0.88);
            padding: 24px;
            box-sizing: border-box;
        }}
        .plot-modal:target {{
            display: flex;
            align-items: center;
            justify-content: center;
        }}
        .plot-modal img {{
            max-width: 96vw;
            max-height: 92vh;
            object-fit: contain;
            background: #fff;
            border-radius: 8px;
        }}
        .plot-modal-close {{
            position: fixed;
            top: 16px;
            right: 24px;
            color: #fff;
            font-size: 2rem;
            text-decoration: none;
            line-height: 1;
        }}
        .plot-modal-title {{
            position: fixed;
            top: 20px;
            left: 24px;
            color: #fff;
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
            {run_meta}
            <p>Generated: {generated_at}</p>
        </header>

        {sampling_section}

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
                <div class="card-label">Q20 rate</div>
                <div class="card-value">{self._format_pct(q20_pct)}</div>
            </div>
            <div class="card">
                <div class="card-label">Q30 rate</div>
                <div class="card-value">{self._format_pct(q30_pct)}</div>
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
            {self._nonstandard_base_cards(base_summary)}
        </div>

        {self._adapter_cards_block(adapter, adapter_pct, adapter_trim, adapter_reference, adapter_ref_seq)}

        {adapter_section}

        {r2_section}

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
            <p class="plot-hint">Click any plot to open it full size.</p>
            <div class="plots">
                {self._plots_block()}
            </div>
        </section>
        {self._plot_modals()}

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

    def _format_pct(self, value) -> str:
        if value == "—":
            return "—"
        return f"{value}%"

    def _nonstandard_base_cards(self, summary: dict) -> str:
        if not summary:
            return ""

        cards = [
            (
                "Unknown bases (N)",
                summary.get("mean_n_pct", 0),
                f"max {summary.get('max_n_pct', 0)}%",
            ),
        ]

        if summary.get("mean_u_pct", 0) > 0:
            cards.append(("RNA uracil (U)", summary["mean_u_pct"], "avg across positions"))

        if summary.get("mean_iupac_pct", 0) > 0:
            cards.append(
                ("IUPAC ambiguous", summary["mean_iupac_pct"], "M R W S K Y V H D B")
            )

        if summary.get("mean_other_pct", 0) > 0:
            cards.append(("Unrecognized symbols", summary["mean_other_pct"], "not in IUPAC set"))

        blocks = []
        for label, mean_pct, hint in cards:
            blocks.append(
                f"""<div class="card">
                <div class="card-label">{html.escape(label)}</div>
                <div class="card-value">{mean_pct}%</div>
                <div class="card-label">{html.escape(hint)}</div>
            </div>"""
            )
        return "\n".join(blocks)

    def _adapter_cards_block(
        self,
        adapter: dict,
        adapter_pct,
        adapter_trim: str,
        adapter_reference: str,
        adapter_ref_seq: str = "",
    ) -> str:
        cards = self._adapter_cards(
            adapter, adapter_pct, adapter_trim, adapter_reference, adapter_ref_seq
        )
        if not cards:
            return ""
        return f'<div class="grid-adapter">{cards}</div>'

    def _adapter_cards(
        self,
        adapter: dict,
        adapter_pct,
        adapter_trim: str,
        adapter_reference: str,
        adapter_ref_seq: str = "",
    ) -> str:
        if adapter_pct == "—" and not adapter_trim:
            return ""

        trim_value = adapter_trim or "—"
        ref_seq_value = adapter_ref_seq or trim_value

        method = adapter.get("detection_method", "none")
        method_label = {
            "reference": "Reference match",
            "de_novo": "De novo fallback",
            "none": "—",
        }.get(method, method)

        reference_value = adapter_reference or "—"
        identity = adapter.get("identity_pct", 0)
        identity_value = f"{identity}%" if identity else "—"

        return f"""
            <div class="card">
                <div class="card-label">Adapter content</div>
                <div class="card-value">{self._format_pct(adapter_pct)}</div>
            </div>
            <div class="card">
                <div class="card-label">Matched reference</div>
                <div class="card-value card-value--text">{html.escape(reference_value)}</div>
            </div>
            <div class="card">
                <div class="card-label">Detection / identity</div>
                <div class="card-value card-value--text">{html.escape(method_label)} / {identity_value}</div>
            </div>
            <div class="card card--seq">
                <div class="card-label">fastp trim sequence</div>
                <div class="card-value card-value--mono">{html.escape(trim_value)}</div>
            </div>
            <div class="card card--seq">
                <div class="card-label">Full adapter / motif</div>
                <div class="card-value card-value--mono">{html.escape(ref_seq_value)}</div>
            </div>"""

    def _run_meta_lines(self) -> str:
        layout = self.sample_plan.get("layout") or self.scout_report.get("layout", "single")
        library_type = self.sample_plan.get("library_type") or self.scout_report.get(
            "library_type", "genome"
        )
        lines = [
            f"<p>Layout: <strong>{html.escape(layout)}</strong> | "
            f"Library: <strong>{html.escape(library_type)}</strong></p>"
        ]
        if self.fastq_r2 is not None:
            lines.append(f"<p>R2: <strong>{html.escape(str(self.fastq_r2))}</strong></p>")
        return "\n".join(lines)

    def _r2_summary_section(self) -> str:
        if not self.r2_metrics:
            return ""

        quality = self.r2_metrics.get("Per position quality", {})
        adapter = self.r2_metrics.get("Adapter discovery", {})
        duplicates = self.r2_metrics.get("Duplicates rate", 0)
        adapter_pct = adapter.get("adapter_content_pct", "—")
        adapter_trim = adapter.get("trim_sequence") or adapter.get("consensus", "")
        adapter_reference = adapter.get("reference_name", "—")

        return f"""
        <section>
            <h2>R2 summary</h2>
            <div class="grid-adapter">
                <div class="card">
                    <div class="card-label">R2 mean quality</div>
                    <div class="card-value">{quality.get("overall_mean", "—")}</div>
                </div>
                <div class="card">
                    <div class="card-label">R2 duplicate rate</div>
                    <div class="card-value">{duplicates}%</div>
                </div>
                <div class="card">
                    <div class="card-label">R2 adapter content</div>
                    <div class="card-value">{self._format_pct(adapter_pct)}</div>
                </div>
                <div class="card card--seq">
                    <div class="card-label">R2 fastp sequence</div>
                    <div class="card-value card-value--mono">{html.escape(adapter_trim or "—")}</div>
                </div>
                <div class="card card--seq">
                    <div class="card-label">R2 matched reference</div>
                    <div class="card-value card-value--text">{html.escape(adapter_reference)}</div>
                </div>
            </div>
        </section>"""

    def _adapter_section(self, adapter: dict) -> str:
        if not adapter:
            return ""

        candidates = adapter.get("candidates", [])
        if not candidates and not adapter.get("consensus"):
            return ""

        method = adapter.get("detection_method", "none")
        if method == "reference":
            intro = (
                "Matched read tails against known adapter references. "
                "Trim sequence is the shortest motif suitable for fastp."
            )
        elif method == "de_novo":
            intro = (
                "No known reference matched; short de novo motif from enriched tail k-mers."
            )
        else:
            intro = "Analyzed read tails for adapter signal."

        rows = []
        for candidate in candidates:
            ref_name = candidate.get("reference_name") or "Unknown"
            trim_seq = candidate.get("trim_sequence") or candidate.get("sequence", "")
            ref_seq = candidate.get("reference_sequence") or trim_seq
            identity = candidate.get("identity_pct", 0)
            identity_label = f"{identity}%" if identity else "—"
            rows.append(
                "<tr>"
                f"<td>{html.escape(ref_name)}</td>"
                f'<td class="seq-cell">{html.escape(trim_seq)}</td>'
                f'<td class="seq-cell">{html.escape(ref_seq)}</td>'
                f"<td>{candidate.get('reads_pct', 0)}%</td>"
                f"<td>{identity_label}</td>"
                "</tr>"
            )

        table_rows = "".join(rows) if rows else (
            f"<tr><td colspan='5'>{html.escape(adapter.get('consensus', 'No candidates'))}</td></tr>"
        )

        return f"""
        <section>
            <h2>Adapter discovery</h2>
            <p>{intro} Analyzed {adapter.get('reads_analyzed', 0):,} read tails.</p>
            <table class="adapter-table">
                <colgroup>
                    <col style="width:18%">
                    <col style="width:28%">
                    <col style="width:28%">
                    <col style="width:12%">
                    <col style="width:14%">
                </colgroup>
                <thead>
                    <tr>
                        <th align="left">Reference</th>
                        <th align="left">fastp sequence</th>
                        <th align="left">Full adapter</th>
                        <th align="left">Reads matched</th>
                        <th align="left">Identity</th>
                    </tr>
                </thead>
                <tbody>{table_rows}</tbody>
            </table>
        </section>"""

    def _sampling_section(self) -> str:
        if not self.sample_plan:
            return ""

        if "R1" in self.sample_plan and "R2" in self.sample_plan:
            blocks = []
            for label in ("R1", "R2"):
                blocks.append(self._sampling_block(self.sample_plan[label], label))
            return "\n".join(blocks)

        return self._sampling_block(self.sample_plan, self.sample_plan.get("read_label", "R1"))

    def _sampling_block(self, plan: dict, title: str) -> str:
        total_reads = plan.get("total_reads")
        sample_budget = plan.get("sample_budget")
        if total_reads is None or sample_budget is None:
            return ""

        reads_processed = plan.get("reads_processed", sample_budget)
        mode = plan.get("mode", "auto")
        rows = [
            ("File", plan.get("fastq", "—")),
            ("Layout", plan.get("layout", self.sample_plan.get("layout", "—"))),
            ("Library type", plan.get("library_type", self.sample_plan.get("library_type", "—"))),
            ("Mode", mode),
            ("Total reads in file", f"{total_reads:,}"),
            ("Reads analyzed", f"{reads_processed:,}"),
            ("Sample fraction", f"{plan.get('sample_fraction_pct', '—')}%"),
            ("Confidence target", f"{plan.get('confidence_pct', 95.0)}%"),
            ("Margin (rates)", f"±{plan.get('margin_rate', '—')}"),
            ("Margin (mean PHRED)", f"±{plan.get('margin_mean', '—')}"),
            ("n (proportion formula)", plan.get("n_proportion", "—")),
            ("n (mean formula)", plan.get("n_mean", "—")),
            ("n (base stat)", plan.get("n_stat_base", "—")),
        ]

        if mode == "with_adapter":
            rows.extend([
                ("n (adapter formula)", plan.get("n_proportion_adapter", "—")),
                ("n (adapter stat)", plan.get("n_stat_adapter", "—")),
                (
                    "Adapter-oriented budget",
                    f"{plan.get('sample_budget_for_adapters', '—'):,} reads "
                    f"({plan.get('sample_fraction_pct_for_adapters', '—')}%)",
                ),
            ])

        if mode in ("base", "with_adapter") and plan.get("sample_budget_base") is not None:
            rows.append(("Base QC budget", f"{plan['sample_budget_base']:,} reads"))

        table_rows = "".join(
            f"<tr><td>{html.escape(str(label))}</td><td><strong>{html.escape(str(value))}</strong></td></tr>"
            for label, value in rows
        )

        return f"""
        <section>
            <h2>Sampling — {html.escape(title)}</h2>
            <table style="width:100%; border-collapse: collapse;">
                <tbody>{table_rows}</tbody>
            </table>
        </section>"""

    def _list_block(self, items: list[str], empty_text: str) -> str:
        if not items:
            return f'<p class="empty">{html.escape(empty_text)}</p>'
        rows = "".join(f"<li>{html.escape(item)}</li>" for item in items)
        return f"<ul>{rows}</ul>"

    def _plot_id(self, title: str) -> str:
        safe = "".join(ch if ch.isalnum() else "-" for ch in title.lower()).strip("-")
        return f"plot-{safe}"

    def _plots_block(self) -> str:
        if not self.plot_paths:
            return '<p class="empty">No plots available.</p>'

        blocks = []
        for title, path in self.plot_paths.items():
            plot_id = self._plot_id(title)
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            blocks.append(
                f"""<div class="plot-card">
                    <h3>{html.escape(title)}</h3>
                    <a href="#{plot_id}">
                        <img src="data:image/png;base64,{encoded}" alt="{html.escape(title)}">
                    </a>
                </div>"""
            )
        return "\n".join(blocks)

    def _plot_modals(self) -> str:
        if not self.plot_paths:
            return ""

        modals = []
        for title, path in self.plot_paths.items():
            plot_id = self._plot_id(title)
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            modals.append(
                f"""<div class="plot-modal" id="{plot_id}">
                    <a class="plot-modal-close" href="#" aria-label="Close">&times;</a>
                    <div class="plot-modal-title">{html.escape(title)}</div>
                    <img src="data:image/png;base64,{encoded}" alt="{html.escape(title)}">
                </div>"""
            )
        return "\n".join(modals)


def build_plot_paths(
    metrics: dict,
    output_dir: Path,
    name_suffix: str = "",
) -> dict[str, Path]:
    from fastq_scout.plot import MetricPlotter

    plot_titles = {
        "Per position quality": "Per position quality",
        "Per sequence quality scores": "Per sequence quality scores",
        "Length distribution": "Length distribution",
        "GC content": "GC content",
        "Per base sequence content": "Per base sequence content",
        "Duplicates rate": "Duplicate rate",
        "Adapter discovery": "Adapter enrichment",
    }

    plot_paths = {}
    for metric_name, data in metrics.items():
        plot_path = MetricPlotter(metric_name, data).plot(output_dir)
        if name_suffix:
            renamed = plot_path.with_name(f"{plot_path.stem}{name_suffix}{plot_path.suffix}")
            plot_path.rename(renamed)
            plot_path = renamed
        plot_paths[plot_titles.get(metric_name, metric_name)] = plot_path

    return plot_paths
