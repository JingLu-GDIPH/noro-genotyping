#!/usr/bin/env python3
"""Render a self-contained HTML table report for norovirus genotyping results.

The report presents, for every query sequence, its genogroup, RdRp and VP1
genotypes with SH-like support and confidence, the genomic position (ORF and
specific protein such as NTPase / VPg / 3CLpro / 5'UTR), and — for sequences
that do not cover RdRp or VP1 — the closest genome references.

Everything is written into a single ``*.html`` file with no external
dependencies, so it can be emailed or opened offline.
"""

import argparse
import csv
import html
from pathlib import Path


# Threshold below which a genotype call is flagged "unsure". Kept in sync with
# genotype_norovirus.py's default so the report's colouring matches the TSV.
BOOTSTRAP_THRESHOLD = 0.75

REGION_LABELS = {"rdrp": "RdRp", "vp1": "VP1"}


def load_results(tsv_path):
    """Read the genotyping TSV into a list of dicts (preserve column order)."""
    with tsv_path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def confidence_tier(value):
    """Map a numeric bootstrap onto a CSS colour tier label."""
    if value in (None, "", "NA"):
        return "none"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "none"
    if v >= BOOTSTRAP_THRESHOLD:
        return "high"
    if v >= 0.5:
        return "medium"
    return "low"


def status_badge_class(status):
    """Return a CSS class for the status badge."""
    mapping = {
        "typed": "ok",
        "partial_genotype": "warn",
        "genotype_unsure": "warn",
        "genogroup_only": "muted",
        "off_region": "muted",
        "unclassified": "bad",
    }
    return mapping.get(status, "muted")


CSS = """
:root {
  --bg: #0f1115;
  --panel: #171a21;
  --panel-2: #1e222b;
  --text: #e6e9ef;
  --muted: #9aa3ad;
  --border: #2a2f3a;
  --accent: #4c9aff;
  --ok: #2ea043;
  --warn: #d29922;
  --bad: #f85149;
  --high: #2ea043;
  --medium: #d29922;
  --low: #f85149;
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 0;
  background: var(--bg); color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  font-size: 14px; line-height: 1.5;
}
.container { max-width: 1280px; margin: 0 auto; padding: 32px 24px 80px; }
header { border-bottom: 1px solid var(--border); padding-bottom: 20px; margin-bottom: 28px; }
header h1 { margin: 0 0 6px; font-size: 24px; font-weight: 600; }
header .subtitle { color: var(--muted); font-size: 13px; }
.summary { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin: 20px 0 32px; }
.summary .card { background: var(--panel); border: 1px solid var(--border); border-radius: 8px; padding: 14px 16px; }
.summary .card .label { color: var(--muted); font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; }
.summary .card .value { font-size: 22px; font-weight: 600; margin-top: 4px; }
h2 { font-size: 18px; margin: 36px 0 14px; font-weight: 600; border-bottom: 1px solid var(--border); padding-bottom: 8px; }
h2 .count { color: var(--muted); font-weight: 400; font-size: 14px; margin-left: 8px; }
table { border-collapse: collapse; width: 100%; margin: 6px 0 14px; font-size: 13px; }
th, td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--border); vertical-align: top; }
th { color: var(--muted); font-weight: 500; font-size: 11px; text-transform: uppercase; letter-spacing: 0.4px; }
td .mono, .mono { font-family: "SF Mono", Menlo, Consolas, monospace; font-size: 12px; }
tbody tr:hover { background: var(--panel-2); }
.geno { font-weight: 600; }
.support.high { color: var(--high); font-weight: 600; }
.support.medium { color: var(--medium); font-weight: 600; }
.support.low { color: var(--low); font-weight: 600; }
.badge { display: inline-block; padding: 2px 9px; border-radius: 12px; font-size: 11px; font-weight: 600; letter-spacing: 0.3px; text-transform: uppercase; white-space: nowrap; }
.badge.ok { background: rgba(46,160,67,0.18); color: var(--ok); border: 1px solid rgba(46,160,67,0.4); }
.badge.warn { background: rgba(210,153,34,0.18); color: var(--warn); border: 1px solid rgba(210,153,34,0.4); }
.badge.muted { background: rgba(154,163,173,0.15); color: var(--muted); border: 1px solid rgba(154,163,173,0.35); }
.badge.bad { background: rgba(248,81,73,0.18); color: var(--bad); border: 1px solid rgba(248,81,73,0.4); }
.protein-tag { display: inline-block; padding: 1px 7px; border-radius: 4px; font-size: 11px; background: var(--panel-2); border: 1px solid var(--border); color: var(--accent); font-family: "SF Mono", Menlo, monospace; margin: 1px 2px; }
.offref { font-size: 11px; color: var(--muted); }
.offref .best { color: var(--text); font-weight: 600; }
footer { color: var(--muted); font-size: 12px; margin-top: 40px; border-top: 1px solid var(--border); padding-top: 16px; }
"""


def fmt(value, default="—"):
    if value is None or value == "":
        return default
    return html.escape(str(value))


def summary_cards(rows):
    """Return HTML for the summary stat cards."""
    counts = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    order = ["typed", "partial_genotype", "genotype_unsure",
             "genogroup_only", "off_region", "unclassified"]
    cards = [("<strong>Total queries</strong>", str(len(rows)))]
    for status in order:
        if status in counts:
            label = status.replace("_", " ")
            cards.append((label, str(counts[status])))
    return "".join(
        f'<div class="card"><div class="label">{html.escape(label)}</div>'
        f'<div class="value">{value}</div></div>'
        for label, value in cards
    )


def confidence_cell(confidence):
    """Render the typed/unsure confidence as a badge."""
    if confidence == "unsure":
        return '<span class="badge warn">unsure</span>'
    if confidence == "typed":
        return '<span class="badge ok">typed</span>'
    return fmt(confidence)


def support_cell(bootstrap):
    """Render the SH support value coloured by tier."""
    if not bootstrap:
        return "—"
    tier = confidence_tier(bootstrap)
    return f'<span class="support {tier}">{html.escape(bootstrap)}</span>'


def genotype_cell(row, region):
    """Render a genotype with its confidence badge inline."""
    genotype = row.get(f"{region}_genotype", "")
    confidence = row.get(f"{region}_genotype_confidence", "")
    if not genotype:
        return "—"
    badge = ""
    if confidence == "unsure":
        badge = ' <span class="badge warn">unsure</span>'
    elif confidence == "typed":
        badge = ' <span class="badge ok">typed</span>'
    return f'<span class="geno">{html.escape(genotype)}</span>{badge}'


def position_cell(row):
    """Render the genomic position (region / proteins / coords) compactly."""
    region = row.get("query_genome_region", "")
    proteins = row.get("query_genome_proteins", "")
    coords = row.get("query_genome_coords", "")
    parts = []
    if region:
        parts.append(f'<span class="badge muted">{html.escape(region)}</span>')
    if proteins and proteins != "intergenic":
        protein_list = proteins.split("+")
        # Cap the number of tags shown inline; for whole-genome queries that
        # overlap every protein the full list would overwhelm the cell, so we
        # show the first few and summarise the rest.
        max_inline = 5
        shown = protein_list[:max_inline]
        tags = "".join(
            f'<span class="protein-tag">{html.escape(p)}</span>'
            for p in shown
        )
        if len(protein_list) > max_inline:
            tags += (f'<span class="offref"> +{len(protein_list) - max_inline} '
                     "more</span>")
        parts.append(tags)
    if coords:
        parts.append(f'<span class="mono">{html.escape(coords)}</span>')
    return " ".join(parts) if parts else "—"


def offref_cell(row):
    """Render the closest references for off-region queries."""
    refs = row.get("off_region_top_references", "")
    if not refs:
        return "—"
    items = [r.strip() for r in refs.split(",") if r.strip()]
    if not items:
        return "—"
    rendered = []
    for idx, item in enumerate(items[:3]):
        cls = "best" if idx == 0 else ""
        rendered.append(f'<span class="{cls}">{html.escape(item)}</span>')
    return '<div class="offref">' + ", ".join(rendered) + "</div>"


def genotype_summary_table(rows, region):
    """Build the per-region summary table for the report top."""
    label = REGION_LABELS[region]
    header = (
        "<table><thead><tr>"
        "<th>Query</th><th>Genogroup</th>"
        f"<th>{label} genotype</th><th>Confidence</th><th>SH support</th>"
        "<th>Nearest reference</th><th>Tree distance</th>"
        "<th>BLAST identity</th><th>Ref coverage</th>"
        "</tr></thead><tbody>"
    )
    body = []
    for row in rows:
        genotype = row.get(f"{region}_genotype", "")
        nearest = row.get(f"{region}_nearest_reference", "")
        identity = row.get(f"{region}_blast_identity_pct", "")
        coverage = row.get(f"{region}_reference_coverage_pct", "")
        body.append(
            "<tr>"
            f'<td class="mono">{fmt(row["sequence_name"])}</td>'
            f'<td>{fmt(row["genogroup"])}</td>'
            f"<td>{genotype_cell(row, region)}</td>"
            f"<td>{confidence_cell(row.get(f'{region}_genotype_confidence',''))}</td>"
            f"<td>{support_cell(row.get(f'{region}_bootstrap',''))}</td>"
            f'<td class="mono">{fmt(nearest)}</td>'
            f'<td class="mono">{fmt(row.get(f"{region}_tree_distance",""))}</td>'
            f'<td class="mono">{fmt(identity)}{ "%" if identity else "" }</td>'
            f'<td class="mono">{fmt(coverage)}{ "%" if coverage else "" }</td>'
            "</tr>"
        )
    return header + "".join(body) + "</tbody></table>"


def position_table(rows):
    """Build the genomic-position / off-region table."""
    header = (
        "<table><thead><tr>"
        "<th>Query</th><th>Length</th><th>Genogroup</th><th>Status</th>"
        "<th>Genomic position</th><th>Coarse reference (identity)</th>"
        "<th>Closest references (off-region)</th>"
        "</tr></thead><tbody>"
    )
    body = []
    for row in rows:
        coarse_ref = row.get("coarse_reference", "")
        coarse_id = row.get("coarse_identity_pct", "")
        coarse_cov = row.get("coarse_query_coverage_pct", "")
        coarse_html = "—"
        if coarse_ref:
            coarse_html = (
                f'<span class="mono">{fmt(coarse_ref)}</span>'
                f' <span class="offref">({fmt(coarse_id)}% id / {fmt(coarse_cov)}% cov)</span>'
            )
        body.append(
            "<tr>"
            f'<td class="mono">{fmt(row["sequence_name"])}</td>'
            f'<td class="mono">{fmt(row["sequence_length"])} nt</td>'
            f"<td>{fmt(row['genogroup'])}</td>"
            f'<td><span class="badge {status_badge_class(row["status"])}">{html.escape(row["status"])}</span></td>'
            f"<td>{position_cell(row)}</td>"
            f"<td>{coarse_html}</td>"
            f"<td>{offref_cell(row)}</td>"
            "</tr>"
        )
    return header + "".join(body) + "</tbody></table>"


def render_report(rows, tree_dir, reference_names, output_path, title, input_name):
    """Assemble the full HTML document and write it to ``output_path``.

    ``tree_dir`` and ``reference_names`` are accepted for backward compatibility
    with the ``genotype_norovirus.py`` caller but are no longer used: the report
    is table-only by design.
    """
    del tree_dir, reference_names  # table-only report; kept for API stability

    summary_html = summary_cards(rows)

    # RdRp / VP1 sections: only show queries that attempted that region.
    rdrp_rows = [r for r in rows if r.get("rdrp_genotype")]
    vp1_rows = [r for r in rows if r.get("vp1_genotype")]
    rdrp_section = (
        f"<h2>RdRp (P-type) genotypes <span class='count'>"
        f"{len(rdrp_rows)} of {len(rows)} queries</span></h2>"
        + genotype_summary_table(rows, "rdrp")
    )
    vp1_section = (
        f"<h2>VP1 (capsid) genotypes <span class='count'>"
        f"{len(vp1_rows)} of {len(rows)} queries</span></h2>"
        + genotype_summary_table(rows, "vp1")
    )

    position_section = (
        f"<h2>Genomic positions <span class='count'>{len(rows)} queries</span></h2>"
        + position_table(rows)
    )

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="container">
<header>
<h1>{html.escape(title)}</h1>
<div class="subtitle">Input: <code>{html.escape(input_name)}</code> &middot;
{len(rows)} query sequence(s)</div>
</header>
<section class="summary">{summary_html}</section>
{rdrp_section}
{vp1_section}
{position_section}
<footer>
SH-like branch support values are FastTree local supports (&minus;nome
&minus;mllen) at the query's placement clade. Genotypes are flagged
<em>unsure</em> when the placement support is below {BOOTSTRAP_THRESHOLD}.
Genomic positions are projected from the best genome BLAST hit onto the gene
coordinate table: ORF1 cleavage products (p48 / NTPase / p22 / VPg / 3CLpro /
RdRp) use conserved relative boundaries, VP1/VP2 and the UTRs come from the
GenBank annotation, and RdRp coordinates come from the region-library BLAST.
</footer>
</div>
</body>
</html>
"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Render a self-contained HTML table report from the "
                    "norovirus genotyping TSV (genotypes, SH support, "
                    "genomic positions, closest references)."
    )
    parser.add_argument("--tsv", required=True, type=Path,
                        help="genotype_norovirus.py *_genotyping.tsv output")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output HTML path")
    parser.add_argument("--title", default="Norovirus Genotyping Report")
    parser.add_argument("--bootstrap-threshold", type=float,
                        default=BOOTSTRAP_THRESHOLD,
                        help="SH support threshold used to colour tiers")
    parser.add_argument("--trees", type=Path, default=None,
                        help="(ignored) kept for backward compatibility")
    parser.add_argument("--ref-dir", type=Path, default=None,
                        help="(ignored) kept for backward compatibility")
    args = parser.parse_args()

    rows = load_results(args.tsv)
    input_name = rows[0]["input_file"] if rows else args.tsv.name
    render_report(rows, args.trees, {}, args.output, args.title, input_name)
    print(f"Wrote HTML report to {args.output}")
    print(f"  {len(rows)} queries")


if __name__ == "__main__":
    main()
