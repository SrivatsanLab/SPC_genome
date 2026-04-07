#!/usr/bin/env python3
"""
Generate a self-contained HTML dashboard for CapWGS pipeline results.

Reads compiled QC metrics, preprocessing stats, and existing plots from
a results directory and produces a single HTML file with interactive
Plotly.js charts, styled after 10x Genomics Cell Ranger web_summary.html.

Usage:
    python3 generate_dashboard.py <results_dir> [--output dashboard.html] [--title "Sample"]

    # For old runs: first parse SLURM logs, then generate dashboard
    python3 parse_slurm_logs.py SLURM_outs/array_outs results/sample/ --job-id 12345
    python3 generate_dashboard.py results/sample/
"""

import argparse
import base64
import csv
import json
import os
import sys


def load_json(filepath):
    if not os.path.exists(filepath):
        return None
    with open(filepath) as f:
        return json.load(f)


def load_csv(filepath):
    if not os.path.exists(filepath):
        return None
    rows = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def load_csv_raw(filepath):
    """Load CSV and return header + rows."""
    if not os.path.exists(filepath):
        return None, None
    with open(filepath) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    return header, rows


def encode_image_base64(filepath):
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode('utf-8')


def load_flagstat(filepath):
    if not os.path.exists(filepath):
        return None
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
    from preprocessing_parsers import parse_flagstat
    with open(filepath) as f:
        text = f.read()
    result = parse_flagstat(text)
    return result if result else None


def parse_barcode_assignment(filepath):
    if not os.path.exists(filepath):
        return None
    result = {}
    with open(filepath) as f:
        text = f.read()
    import re
    m = re.search(r'Total input reads.*?:\s*([\d,]+)', text)
    if m:
        result['total_input_reads'] = int(m.group(1).replace(',', ''))
    m = re.search(r'Reads assigned.*?:\s*([\d,]+)', text)
    if m:
        result['reads_assigned'] = int(m.group(1).replace(',', ''))
    m = re.search(r'Assignment rate:\s*([\d.]+)%', text)
    if m:
        result['assignment_rate'] = float(m.group(1))
    return result


def fmt_number(n):
    """Format a number with commas."""
    if n is None:
        return 'N/A'
    if isinstance(n, float):
        if n < 0.01:
            return '{:.4f}'.format(n)
        return '{:,.2f}'.format(n)
    return '{:,}'.format(int(n))


def fmt_pct(n):
    if n is None:
        return 'N/A'
    return '{:.1f}%'.format(float(n))


def generate_dashboard(results_dir, title='CapWGS Pipeline', output_path='dashboard.html',
                       conditions_json=None, barcode_map=None):
    # Load all available data
    preproc = load_json(os.path.join(results_dir, 'preprocessing_summary.json'))
    flagstat = load_flagstat(os.path.join(results_dir, 'flagstat.txt'))
    assignment = parse_barcode_assignment(os.path.join(results_dir, 'barcode_assignment_stats.txt'))
    qc_metrics = load_csv(os.path.join(results_dir, 'compiled_qc_metrics.csv'))
    readcounts = load_csv(os.path.join(results_dir, 'readcounts.csv'))
    kneeplot_b64 = encode_image_base64(os.path.join(results_dir, 'kneeplot.png'))
    lw_plot_b64 = encode_image_base64(os.path.join(results_dir, 'lander_waterman_coverage.png'))

    # Load adapter histograms
    hist_r1_header, hist_r1_rows = load_csv_raw(os.path.join(results_dir, 'adapter_histogram_r1.csv'))
    hist_r2_header, hist_r2_rows = load_csv_raw(os.path.join(results_dir, 'adapter_histogram_r2.csv'))

    # Load lorenz curves
    lorenz_header, lorenz_rows = load_csv_raw(os.path.join(results_dir, 'compiled_lorenz_curves.csv'))

    # Load condition assignments if provided
    condition_map = {}  # barcode -> condition
    if conditions_json and barcode_map:
        from assign_conditions import assign_conditions_to_barcodes
        with open(conditions_json) as f:
            cond_config = json.load(f)
        real_cells_path = os.path.join(results_dir, 'real_cells.txt')
        if os.path.exists(real_cells_path):
            with open(real_cells_path) as f:
                cell_barcodes = [line.strip() for line in f if line.strip()]
            condition_map = assign_conditions_to_barcodes(cond_config, barcode_map, cell_barcodes)

    # Count cells
    real_cells_path = os.path.join(results_dir, 'real_cells.txt')
    n_cells = 0
    if os.path.exists(real_cells_path):
        with open(real_cells_path) as f:
            n_cells = sum(1 for line in f if line.strip())

    # Compute summary stats from QC metrics
    qc_summary = {}
    if qc_metrics:
        def safe_median(values):
            vals = sorted([v for v in values if v is not None])
            if not vals:
                return None
            mid = len(vals) // 2
            if len(vals) % 2 == 0:
                return (vals[mid - 1] + vals[mid]) / 2
            return vals[mid]

        def safe_float(v):
            try:
                return float(v)
            except (ValueError, TypeError):
                return None

        pct_aligned = [safe_float(r.get('pct_reads_aligned')) for r in qc_metrics]
        pct_dup = [safe_float(r.get('pct_duplicates')) for r in qc_metrics]
        mean_cov = [safe_float(r.get('mean_coverage')) for r in qc_metrics]
        pct_1x = [safe_float(r.get('pct_1x')) for r in qc_metrics]
        gini = [safe_float(r.get('gini_coefficient')) for r in qc_metrics]
        total_reads_per_cell = [safe_float(r.get('total_reads')) for r in qc_metrics]

        qc_summary['median_pct_aligned'] = safe_median([v for v in pct_aligned if v is not None])
        qc_summary['median_pct_duplicates'] = safe_median([v for v in pct_dup if v is not None])
        qc_summary['median_mean_coverage'] = safe_median([v for v in mean_cov if v is not None])
        qc_summary['median_pct_1x'] = safe_median([v for v in pct_1x if v is not None])
        qc_summary['median_gini'] = safe_median([v for v in gini if v is not None])
        valid_reads = [v for v in total_reads_per_cell if v is not None]
        qc_summary['mean_reads_per_cell'] = sum(valid_reads) / len(valid_reads) if valid_reads else None

    # Build JavaScript data for plots
    js_data = {}

    # Adapter histograms
    if hist_r1_rows:
        js_data['adapter_r1'] = {
            'lengths': [int(r[0]) for r in hist_r1_rows],
            'counts': [int(r[1]) for r in hist_r1_rows],
        }
    if hist_r2_rows:
        js_data['adapter_r2'] = {
            'lengths': [int(r[0]) for r in hist_r2_rows],
            'counts': [int(r[1]) for r in hist_r2_rows],
        }

    # Readcounts for knee plot
    if readcounts:
        # Sort by count descending, preserving barcode for condition lookup
        sorted_rc = sorted(readcounts, key=lambda r: int(r['read_count']), reverse=True)
        counts = [int(r['read_count']) for r in sorted_rc]
        rc_conditions = []
        if condition_map:
            rc_conditions = [condition_map.get(r['barcode'], 'Unassigned') for r in sorted_rc]
        js_data['readcounts'] = {
            'ranks': list(range(1, len(counts) + 1)),
            'counts': counts,
            'barcodes': [r['barcode'] for r in sorted_rc],
            'conditions': rc_conditions,
            'n_cells': n_cells,
        }

    # QC metric distributions
    if qc_metrics:
        def extract_values(key):
            return [safe_float(r.get(key)) for r in qc_metrics if safe_float(r.get(key)) is not None]

        # Add condition labels if available
        conditions_list = []
        if condition_map:
            conditions_list = [condition_map.get(r.get('barcode', ''), 'Unassigned')
                               for r in qc_metrics]

        js_data['qc'] = {
            'barcodes': [r.get('barcode', '') for r in qc_metrics],
            'pct_aligned': extract_values('pct_reads_aligned'),
            'pct_duplicates': extract_values('pct_duplicates'),
            'mean_coverage': extract_values('mean_coverage'),
            'pct_1x': extract_values('pct_1x'),
            'gini': extract_values('gini_coefficient'),
            'conditions': conditions_list,
        }

    # Lorenz curves (sample up to 50 cells for readability)
    if lorenz_header and lorenz_rows:
        x_vals = [float(r[0]) for r in lorenz_rows]
        cell_names = lorenz_header[1:]  # first column is cumulative_fraction_genome
        step = max(1, len(cell_names) // 50)
        sampled_cells = cell_names[::step]
        lorenz_data = {'x': x_vals, 'cells': {}}
        for i, name in enumerate(cell_names):
            if name in sampled_cells:
                short_name = name[:12] + '...' if len(name) > 12 else name
                lorenz_data['cells'][short_name] = [float(r[i + 1]) for r in lorenz_rows]
        js_data['lorenz'] = lorenz_data

    # Read funnel data for alluvial/sankey plot
    if preproc:
        demux_data = preproc.get('demux', {})
        trim_data = preproc.get('trimming', {})
        # Demux counts are read pairs; convert everything to individual reads (*2)
        total_input_pairs = demux_data.get('total_input_reads', 0)
        after_demux_pairs = demux_data.get('total_mapped_reads', 0)
        pairs_removed = trim_data.get('pairs_removed', 0)
        after_trim_pairs = after_demux_pairs - pairs_removed

        # Use individual reads (x2) for consistent units with flagstat
        total_input = total_input_pairs * 2
        after_demux = after_demux_pairs * 2
        after_trim = after_trim_pairs * 2

        funnel = {
            'labels': ['Input Reads', 'After Demux', 'After Trimming', 'Mapped', 'Properly Paired'],
            'values': [total_input, after_demux, after_trim],
        }

        if flagstat:
            funnel['values'].append(flagstat.get('mapped', after_trim))
            funnel['values'].append(flagstat.get('properly_paired', 0))
        else:
            # Without flagstat, stop at after trimming
            funnel['labels'] = funnel['labels'][:3]

        js_data['funnel'] = funnel

    # Reference genome
    ref_genome = ''
    if preproc:
        ref_genome = preproc.get('alignment', {}).get('reference_genome', '')

    # Build HTML
    html = _build_html(
        title=title,
        n_cells=n_cells,
        assignment=assignment,
        preproc=preproc,
        flagstat=flagstat,
        qc_summary=qc_summary,
        ref_genome=ref_genome,
        kneeplot_b64=kneeplot_b64,
        lw_plot_b64=lw_plot_b64,
        js_data=js_data,
    )

    with open(output_path, 'w') as f:
        f.write(html)
    print('Dashboard written to: {}'.format(output_path))


def _build_html(title, n_cells, assignment, preproc, flagstat, qc_summary, ref_genome,
                kneeplot_b64, lw_plot_b64, js_data):

    # Summary card values
    cards = []
    cards.append(('Estimated Cells', fmt_number(n_cells), ''))
    if qc_summary.get('mean_reads_per_cell') is not None:
        cards.append(('Mean Reads per Cell', fmt_number(int(qc_summary['mean_reads_per_cell'])), ''))

    if preproc:
        rate = preproc.get('demux', {}).get('barcode_mapping_rate')
        if rate is not None:
            cards.append(('Barcode Mapping Rate', fmt_pct(rate * 100), ''))
    elif assignment:
        cards.append(('Barcode Assignment Rate', fmt_pct(assignment.get('assignment_rate')), ''))

    if flagstat and 'mapped_pct' in flagstat:
        cards.append(('Reads Mapped to Genome', fmt_pct(flagstat['mapped_pct']), ''))
    if flagstat and 'properly_paired_pct' in flagstat:
        cards.append(('Properly Paired', fmt_pct(flagstat['properly_paired_pct']), ''))

    if preproc:
        r1_pct = preproc.get('trimming', {}).get('r1', {}).get('pct_reads_with_adapters')
        if r1_pct is not None:
            cards.append(('Reads with Adapters (R1)', fmt_pct(r1_pct), ''))

    if qc_summary.get('median_pct_aligned') is not None:
        val = qc_summary['median_pct_aligned']
        cards.append(('Median Alignment Rate', fmt_pct(val * 100) if val <= 1 else fmt_pct(val), ''))
    if qc_summary.get('median_pct_duplicates') is not None:
        cards.append(('Median Duplication Rate', fmt_pct(qc_summary['median_pct_duplicates'] * 100), ''))
    if qc_summary.get('median_mean_coverage') is not None:
        cards.append(('Median Coverage', '{:.4f}x'.format(qc_summary['median_mean_coverage']), ''))
    if qc_summary.get('median_gini') is not None:
        cards.append(('Median Gini Coefficient', '{:.3f}'.format(qc_summary['median_gini']), ''))
    cards_html = ''
    for label, value, css_class in cards:
        cards_html += '''
        <div class="card">
            <div class="card-label">{label}</div>
            <div class="card-value {css_class}">{value}</div>
        </div>'''.format(label=label, value=value, css_class=css_class)

    # Preprocessing table
    preproc_table = ''
    if preproc:
        demux = preproc.get('demux', {})
        trim = preproc.get('trimming', {})
        aln = preproc.get('alignment', {})
        r1 = trim.get('r1', {})
        r2 = trim.get('r2', {})

        # Compute basepair percentages
        r1_bp_pct = ''
        if r1.get('bp_written') and r1.get('total_bp_processed'):
            r1_bp_pct = ' ({})'.format(fmt_pct(100 * r1['bp_written'] / r1['total_bp_processed']))
        r2_bp_pct = ''
        if r2.get('bp_written') and r2.get('total_bp_processed'):
            r2_bp_pct = ' ({})'.format(fmt_pct(100 * r2['bp_written'] / r2['total_bp_processed']))

        # Reads sent to alignment = pairs after trim * 2 (both mates)
        total_reads_to_bwa = aln.get('total_reads_processed', 0)
        pairs_after_trim = r1.get('total_reads', 0) - trim.get('pairs_removed', 0)
        reads_to_bwa_str = fmt_number(total_reads_to_bwa) if total_reads_to_bwa else fmt_number(pairs_after_trim * 2)

        rows = [
            ('Reference Genome', aln.get('reference_genome', 'N/A')),
            ('', ''),
            ('Total Input Reads (paired-end)', fmt_number(demux.get('total_input_reads'))),
            ('Reads with Valid Barcode', fmt_number(demux.get('total_mapped_reads'))),
            ('Barcode Mapping Rate', fmt_pct((demux.get('barcode_mapping_rate', 0)) * 100)),
            ('', ''),
            ('R1 Reads with Adapters', '{} ({})'.format(
                fmt_number(r1.get('reads_with_adapters')),
                fmt_pct(r1.get('pct_reads_with_adapters')))),
            ('R1 Basepairs Written', '{} / {}{}'.format(
                fmt_number(r1.get('bp_written')),
                fmt_number(r1.get('total_bp_processed')),
                r1_bp_pct)),
            ('R2 Reads with Adapters', '{} ({})'.format(
                fmt_number(r2.get('reads_with_adapters')),
                fmt_pct(r2.get('pct_reads_with_adapters')))),
            ('R2 Basepairs Written', '{} / {}{}'.format(
                fmt_number(r2.get('bp_written')),
                fmt_number(r2.get('total_bp_processed')),
                r2_bp_pct)),
            ('Pairs Removed (length cutoff)', '{} ({})'.format(
                fmt_number(trim.get('pairs_removed')),
                fmt_pct(trim.get('pct_pairs_removed')))),
            ('', ''),
            ('Reads Sent to Alignment', reads_to_bwa_str),
            ('BWA Version', aln.get('bwa_version', 'N/A')),
            ('Insert Size (mean / median)', '{} / {}'.format(
                aln.get('insert_size_mean', 'N/A'),
                aln.get('insert_size_median', 'N/A'))),
        ]

        # Add flagstat rows if available
        if flagstat:
            rows.append(('', ''))
            rows.append(('Reads Mapped to Genome', '{} ({})'.format(
                fmt_number(flagstat.get('mapped')),
                fmt_pct(flagstat.get('mapped_pct')))))
            rows.append(('Properly Paired', '{} ({})'.format(
                fmt_number(flagstat.get('properly_paired')),
                fmt_pct(flagstat.get('properly_paired_pct')))))
            rows.append(('Singletons', '{} ({})'.format(
                fmt_number(flagstat.get('singletons')),
                fmt_pct(flagstat.get('singletons_pct')))))
            rows.append(('Secondary Alignments', fmt_number(flagstat.get('secondary'))))
            rows.append(('Supplementary Alignments', fmt_number(flagstat.get('supplementary'))))

        preproc_table = '<table class="stats-table">'
        for label, value in rows:
            if not label:
                preproc_table += '<tr><td colspan="2" class="spacer"></td></tr>'
            else:
                preproc_table += '<tr><td class="stat-label">{}</td><td class="stat-value">{}</td></tr>'.format(label, value)
        preproc_table += '</table>'

    # Embedded images
    kneeplot_img = ''
    if kneeplot_b64:
        kneeplot_img = '<img src="data:image/png;base64,{}" alt="Knee Plot" class="plot-img">'.format(kneeplot_b64)

    lw_img = ''
    if lw_plot_b64:
        lw_img = '<img src="data:image/png;base64,{}" alt="Lander-Waterman Coverage" class="plot-img">'.format(lw_plot_b64)

    js_data_json = json.dumps(js_data)

    html = '''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title} - Pipeline Dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=DM+Sans:ital,wght@0,400;0,500;0,700&family=DM+Mono:wght@300;400;500&display=swap" rel="stylesheet">
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: 'DM Sans', -apple-system, BlinkMacSystemFont, sans-serif; background: #fafafa; color: #111; }}
.header {{ background: #fff; color: #111; padding: 32px 40px 24px; border-bottom: 1px solid #e0e0e0; }}
.header h1 {{ font-size: 26px; font-weight: 700; letter-spacing: -0.5px; }}
.header .subtitle {{ font-family: 'DM Mono', monospace; font-size: 12px; color: #888; margin-top: 6px; letter-spacing: 0.5px; }}
.tabs {{ background: #fff; border-bottom: 1px solid #e0e0e0; padding: 0 40px; display: flex; gap: 0; position: sticky; top: 0; z-index: 100; }}
.tab {{ padding: 14px 24px; cursor: pointer; font-size: 13px; font-weight: 500; color: #999; border-bottom: 2px solid transparent; transition: all 0.2s; text-transform: uppercase; letter-spacing: 0.5px; }}
.tab:hover {{ color: #111; }}
.tab.active {{ color: #111; border-bottom-color: #111; }}
.content {{ max-width: 1200px; margin: 0 auto; padding: 28px 40px; }}
.section {{ display: none; }}
.section.active {{ display: block; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 14px; margin-bottom: 28px; }}
.card {{ background: #fff; border: 1px solid #e8e8e8; border-radius: 6px; padding: 18px 20px; }}
.card-label {{ font-family: 'DM Mono', monospace; font-size: 10px; text-transform: uppercase; letter-spacing: 0.8px; color: #999; margin-bottom: 8px; }}
.card-value {{ font-size: 26px; font-weight: 700; color: #111; }}
.card-value.small {{ font-family: 'DM Mono', monospace; font-size: 12px; font-weight: 400; color: #555; word-break: break-all; }}
.panel {{ background: #fff; border: 1px solid #e8e8e8; border-radius: 6px; padding: 24px; margin-bottom: 20px; }}
.panel h2 {{ font-size: 14px; font-weight: 700; margin-bottom: 16px; color: #111; text-transform: uppercase; letter-spacing: 0.3px; }}
.panel h3 {{ font-size: 13px; font-weight: 500; margin-bottom: 12px; color: #555; }}
.stats-table {{ width: 100%; border-collapse: collapse; table-layout: fixed; }}
.stats-table td {{ padding: 7px 12px; border-bottom: 1px solid #f0f0f0; font-size: 13px; word-wrap: break-word; overflow-wrap: break-word; }}
.stat-label {{ color: #666; width: 50%; }}
.stat-value {{ font-family: 'DM Mono', monospace; color: #111; font-weight: 400; text-align: right; font-size: 12px; }}
.spacer {{ height: 6px; }}
.plot-container {{ width: 100%; min-height: 400px; }}
.plot-container-sm {{ width: 100%; min-height: 280px; }}
.plot-img {{ max-width: 100%; height: auto; border-radius: 4px; }}
.row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
@media (max-width: 768px) {{ .row {{ grid-template-columns: 1fr; }} }}
.no-data {{ color: #999; font-style: italic; padding: 40px; text-align: center; }}
.footer {{ text-align: center; padding: 28px; color: #bbb; font-family: 'DM Mono', monospace; font-size: 11px; }}
</style>
</head>
<body>

<div class="header">
    <h1>{title}</h1>
    <div class="subtitle">SPC Genome &middot; Pipeline Dashboard</div>
</div>

<div class="tabs">
    <div class="tab active" data-tab="summary">Summary</div>
    <div class="tab" data-tab="qc">Per-Cell QC</div>
    <div class="tab" data-tab="coverage">Coverage</div>
</div>

<div class="content">

<!-- SUMMARY TAB (merged preprocessing + summary cards) -->
<div class="section active" id="tab-summary">
    <div class="cards">{cards_html}</div>
    <div class="panel">
        <h2>Read Retention Through Pipeline</h2>
        <div id="read-funnel-summary" class="plot-container"></div>
    </div>
    {preproc_section}
    <div class="panel">
        <h2>Barcode Rank Plot (Knee Plot)</h2>
        <div id="kneeplot-interactive" class="plot-container"></div>
    </div>
</div>

<!-- PER-CELL QC TAB -->
<div class="section" id="tab-qc">
    <div class="panel">
        <h2>Alignment Rate</h2>
        <div id="plot-alignment" class="plot-container"></div>
    </div>
    <div class="panel">
        <h2>Duplication Rate</h2>
        <div id="plot-duplication" class="plot-container"></div>
    </div>
    <div class="panel">
        <h2>Mean Coverage</h2>
        <div id="plot-coverage" class="plot-container"></div>
    </div>
    <div class="panel">
        <h2>Gini Coefficient</h2>
        <div id="plot-gini" class="plot-container"></div>
    </div>
    <div class="panel">
        <h2>Coverage vs Genome Fraction (1x)</h2>
        <div id="plot-pct1x" class="plot-container"></div>
    </div>
</div>

<!-- COVERAGE TAB -->
<div class="section" id="tab-coverage">
    <div class="panel">
        <h2>Lorenz Curves</h2>
        <div id="plot-lorenz" class="plot-container"></div>
    </div>
    <div class="panel">
        <h2>Lander-Waterman Coverage</h2>
        {lw_section}
    </div>
</div>

</div>

<div class="footer">
    SPC Genome Pipeline &middot; Srivatsan Lab
</div>

<script>
var DATA = {js_data};

// Tab switching
document.querySelectorAll('.tab').forEach(function(tab) {{
    tab.addEventListener('click', function() {{
        document.querySelectorAll('.tab').forEach(function(t) {{ t.classList.remove('active'); }});
        document.querySelectorAll('.section').forEach(function(s) {{ s.classList.remove('active'); }});
        tab.classList.add('active');
        document.getElementById('tab-' + tab.dataset.tab).classList.add('active');
        window.dispatchEvent(new Event('resize'));
    }});
}});

var plotLayout = {{
    margin: {{ t: 20, r: 20, b: 50, l: 60 }},
    font: {{ family: 'DM Sans, -apple-system, sans-serif', size: 12, color: '#333' }},
    paper_bgcolor: '#fff',
    plot_bgcolor: '#fff',
    xaxis: {{ gridcolor: '#f0f0f0' }},
    yaxis: {{ gridcolor: '#f0f0f0' }},
}};

var plotConfig = {{ responsive: true, displayModeBar: false }};

// Read funnel (Sankey diagram)
if (DATA.funnel) {{
    var labels = DATA.funnel.labels;
    var values = DATA.funnel.values;
    // Build sankey: each step flows to the next, with a "lost" node
    var sankeyLabels = [];
    var source = [];
    var target = [];
    var sankeyValues = [];
    var nodeColors = [];

    // Create nodes: each step + "lost" nodes between steps
    for (var i = 0; i < labels.length; i++) {{
        sankeyLabels.push(labels[i] + '\\n(' + values[i].toLocaleString() + ')');
        nodeColors.push('#333');
    }}
    for (var i = 0; i < labels.length - 1; i++) {{
        var lost = values[i] - values[i + 1];
        var lostLabel = 'Lost: ' + labels[i] + ' → ' + labels[i + 1] + '\\n(' + lost.toLocaleString() + ')';
        sankeyLabels.push(lostLabel);
        nodeColors.push('#e0e0e0');
    }}

    // Links: step[i] -> step[i+1] (retained) and step[i] -> lost_node (lost)
    var nSteps = labels.length;
    for (var i = 0; i < nSteps - 1; i++) {{
        // Retained
        source.push(i);
        target.push(i + 1);
        sankeyValues.push(values[i + 1]);
        // Lost
        var lost = values[i] - values[i + 1];
        if (lost > 0) {{
            source.push(i);
            target.push(nSteps + i);
            sankeyValues.push(lost);
        }}
    }}

    var sankeyTrace = {{
        type: 'sankey',
        orientation: 'h',
        node: {{
            pad: 20,
            thickness: 25,
            line: {{ color: '#ccc', width: 1 }},
            label: sankeyLabels,
            color: nodeColors
        }},
        link: {{
            source: source,
            target: target,
            value: sankeyValues,
            color: source.map(function(s, i) {{
                return target[i] < nSteps ? 'rgba(50,50,50,0.25)' : 'rgba(200,200,200,0.3)';
            }})
        }}
    }};
    var sankeyLayout = Object.assign({{}}, plotLayout, {{
        margin: {{ t: 20, r: 20, b: 20, l: 20 }},
    }});
    Plotly.newPlot('read-funnel-summary', [sankeyTrace], sankeyLayout, plotConfig);
}}

// Adapter histograms (combined R1 + R2)
if (DATA.adapter_r1 || DATA.adapter_r2) {{
    var adapterTraces = [];
    if (DATA.adapter_r1) {{
        adapterTraces.push({{
            x: DATA.adapter_r1.lengths,
            y: DATA.adapter_r1.counts,
            type: 'bar',
            name: 'R1',
            marker: {{ color: '#333', opacity: 0.7 }}
        }});
    }}
    if (DATA.adapter_r2) {{
        adapterTraces.push({{
            x: DATA.adapter_r2.lengths,
            y: DATA.adapter_r2.counts,
            type: 'bar',
            name: 'R2',
            marker: {{ color: '#e8710a', opacity: 0.7 }}
        }});
    }}
    var adapterLayout = Object.assign({{}}, plotLayout, {{
        height: 300,
        barmode: 'overlay',
        xaxis: {{ title: 'Adapter Position (bp)', gridcolor: '#eee' }},
        yaxis: {{ title: 'Read Count', gridcolor: '#eee' }},
        legend: {{ x: 0.85, y: 0.95 }}
    }});
    Plotly.newPlot('adapter-hist-combined', adapterTraces, adapterLayout, plotConfig);
}}

// Interactive knee plot
var kneeCondColors = ['#111', '#e8710a', '#0d904f', '#9c27b0', '#e91e63', '#00bcd4', '#ff9800', '#795548'];
if (DATA.readcounts) {{
    var hasKneeConditions = DATA.readcounts.conditions && DATA.readcounts.conditions.length > 0;
    var kneeTraces = [];

    if (hasKneeConditions) {{
        // Color points by condition (scatter, not line)
        var kneeUnique = [];
        DATA.readcounts.conditions.forEach(function(c) {{
            if (kneeUnique.indexOf(c) === -1) kneeUnique.push(c);
        }});
        kneeUnique.forEach(function(cond, ci) {{
            var xv = [], yv = [];
            for (var i = 0; i < DATA.readcounts.ranks.length; i++) {{
                if (DATA.readcounts.conditions[i] === cond) {{
                    xv.push(DATA.readcounts.ranks[i]);
                    yv.push(DATA.readcounts.counts[i]);
                }}
            }}
            kneeTraces.push({{
                x: xv, y: yv, type: 'scatter', mode: 'markers',
                marker: {{ color: kneeCondColors[ci % kneeCondColors.length], size: 4, opacity: 0.7 }},
                name: cond
            }});
        }});
    }} else {{
        kneeTraces.push({{
            x: DATA.readcounts.ranks,
            y: DATA.readcounts.counts,
            type: 'scatter',
            mode: 'lines',
            line: {{ color: '#333', width: 2 }},
            name: 'Barcodes'
        }});
    }}

    var kneeShapes = [];
    if (DATA.readcounts.n_cells > 0) {{
        kneeShapes.push({{
            type: 'line',
            x0: DATA.readcounts.n_cells, x1: DATA.readcounts.n_cells,
            y0: 0, y1: 1, yref: 'paper',
            line: {{ color: 'red', width: 2, dash: 'dash' }}
        }});
    }}
    var kneeLayout = Object.assign({{}}, plotLayout, {{
        xaxis: {{ title: 'Barcode Rank', type: 'log', gridcolor: '#eee' }},
        yaxis: {{ title: 'Read Count', type: 'log', gridcolor: '#eee' }},
        shapes: kneeShapes,
        annotations: DATA.readcounts.n_cells > 0 ? [{{
            x: Math.log10(DATA.readcounts.n_cells),
            y: 1, yref: 'paper',
            text: DATA.readcounts.n_cells + ' cells',
            showarrow: true, arrowhead: 2,
            ax: 40, ay: -30,
            font: {{ color: 'red', size: 13 }}
        }}] : []
    }});
    Plotly.newPlot('kneeplot-interactive', kneeTraces, kneeLayout, plotConfig);
}}

// QC violin/histogram plots
var condColors = ['#111', '#e8710a', '#0d904f', '#9c27b0', '#e91e63', '#00bcd4', '#ff9800', '#795548'];
var hasConditions = DATA.qc && DATA.qc.conditions && DATA.qc.conditions.length > 0;
var uniqueConditions = [];
if (hasConditions) {{
    DATA.qc.conditions.forEach(function(c) {{
        if (uniqueConditions.indexOf(c) === -1) uniqueConditions.push(c);
    }});
}}

function makeConditionPlot(divId, values, xLabel, defaultColor) {{
    if (!values || values.length === 0) return;
    var traces = [];
    if (hasConditions && uniqueConditions.length > 1) {{
        uniqueConditions.forEach(function(cond, ci) {{
            var vals = [];
            for (var i = 0; i < values.length; i++) {{
                if (DATA.qc.conditions[i] === cond) vals.push(values[i]);
            }}
            traces.push({{
                x: vals,
                type: 'histogram',
                name: cond,
                marker: {{ color: condColors[ci % condColors.length], opacity: 0.6 }},
                nbinsx: 25
            }});
        }});
    }} else {{
        traces.push({{
            x: values,
            type: 'histogram',
            marker: {{ color: defaultColor, opacity: 0.7 }},
            nbinsx: 30
        }});
    }}
    var layout = Object.assign({{}}, plotLayout, {{
        xaxis: {{ title: xLabel, gridcolor: '#eee' }},
        yaxis: {{ title: 'Count', gridcolor: '#eee' }},
        barmode: 'overlay',
        bargap: 0.05
    }});
    Plotly.newPlot(divId, traces, layout, plotConfig);
}}

if (DATA.qc) {{
    makeConditionPlot('plot-alignment', DATA.qc.pct_aligned, 'Fraction Aligned', '#333');
    makeConditionPlot('plot-duplication', DATA.qc.pct_duplicates, 'Duplication Rate', '#e8710a');
    makeConditionPlot('plot-coverage', DATA.qc.mean_coverage, 'Mean Coverage', '#0d904f');
    makeConditionPlot('plot-gini', DATA.qc.gini, 'Gini Coefficient', '#9c27b0');

    // Coverage vs pct_1x scatter
    if (DATA.qc.mean_coverage && DATA.qc.pct_1x) {{
        var covTraces = [];
        if (hasConditions && uniqueConditions.length > 1) {{
            uniqueConditions.forEach(function(cond, ci) {{
                var xv = [], yv = [];
                for (var i = 0; i < DATA.qc.mean_coverage.length; i++) {{
                    if (DATA.qc.conditions[i] === cond) {{
                        xv.push(DATA.qc.mean_coverage[i]);
                        yv.push(DATA.qc.pct_1x[i]);
                    }}
                }}
                covTraces.push({{
                    x: xv, y: yv, type: 'scatter', mode: 'markers',
                    marker: {{ color: condColors[ci % condColors.length], size: 6, opacity: 0.6 }},
                    name: cond
                }});
            }});
        }} else {{
            covTraces.push({{
                x: DATA.qc.mean_coverage, y: DATA.qc.pct_1x,
                type: 'scatter', mode: 'markers',
                marker: {{ color: '#333', size: 6, opacity: 0.6 }},
                name: 'Cells'
            }});
        }}
        var covLayout = Object.assign({{}}, plotLayout, {{
            xaxis: {{ title: 'Mean Coverage', gridcolor: '#eee' }},
            yaxis: {{ title: 'Fraction Genome at 1x', gridcolor: '#eee' }}
        }});
        Plotly.newPlot('plot-pct1x', covTraces, covLayout, plotConfig);
    }}
}}

// Lorenz curves
if (DATA.lorenz) {{
    var lorenzTraces = [];
    var cellKeys = Object.keys(DATA.lorenz.cells);
    cellKeys.forEach(function(name, i) {{
        lorenzTraces.push({{
            x: DATA.lorenz.x,
            y: DATA.lorenz.cells[name],
            type: 'scatter',
            mode: 'lines',
            line: {{ width: 1, color: 'rgba(50,50,50,0.25)' }},
            name: name,
            showlegend: false
        }});
    }});
    // Perfect equality line
    lorenzTraces.push({{
        x: [0, 1],
        y: [0, 1],
        type: 'scatter',
        mode: 'lines',
        line: {{ color: '#333', width: 2, dash: 'dash' }},
        name: 'Perfect Uniformity'
    }});
    var lorenzLayout = Object.assign({{}}, plotLayout, {{
        xaxis: {{ title: 'Cumulative Fraction of Genome', gridcolor: '#eee', range: [0, 1] }},
        yaxis: {{ title: 'Cumulative Fraction of Reads', gridcolor: '#eee', range: [0, 1] }},
        showlegend: true,
        legend: {{ x: 0.02, y: 0.98 }}
    }});
    Plotly.newPlot('plot-lorenz', lorenzTraces, lorenzLayout, plotConfig);
}}

</script>
</body>
</html>'''.format(
        title=title,
        cards_html=cards_html,
        kneeplot_panel='<div class="panel"><h2>Barcode Rank Plot</h2>{}</div>'.format(kneeplot_img) if kneeplot_img else '',
        preproc_section=_preproc_section(preproc_table, assignment),
        kneeplot_static=kneeplot_img if kneeplot_img else '<div class="no-data">No knee plot image available</div>',
        lw_section=lw_img if lw_img else '<div class="no-data">No Lander-Waterman plot available</div>',
        js_data=js_data_json,
    )

    return html


def _preproc_section(preproc_table, assignment):
    parts = []

    if preproc_table:
        parts.append('''
        <div class="panel">
            <h2>Preprocessing Statistics</h2>
            {table}
        </div>'''.format(table=preproc_table))

    if not preproc_table and assignment:
        parts.append('''
        <div class="panel">
            <h2>Barcode Assignment</h2>
            <table class="stats-table">
                <tr><td class="stat-label">Total Input Reads</td><td class="stat-value">{total}</td></tr>
                <tr><td class="stat-label">Reads Assigned</td><td class="stat-value">{assigned}</td></tr>
                <tr><td class="stat-label">Assignment Rate</td><td class="stat-value">{rate}</td></tr>
            </table>
        </div>'''.format(
            total=fmt_number(assignment.get('total_input_reads')),
            assigned=fmt_number(assignment.get('reads_assigned')),
            rate=fmt_pct(assignment.get('assignment_rate')),
        ))

    parts.append('''
    <div class="panel">
        <h2>Adapter Position Histogram</h2>
        <div id="adapter-hist-combined" class="plot-container-sm"></div>
    </div>''')

    if not preproc_table and not assignment:
        parts.insert(0, '''
        <div class="panel">
            <div class="no-data">No preprocessing statistics available. Run parse_slurm_logs.py to extract stats from old SLURM logs.</div>
        </div>''')

    return '\n'.join(parts)


def main():
    parser = argparse.ArgumentParser(
        description='Generate an HTML dashboard for CapWGS pipeline results.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('results_dir', help='Results directory containing pipeline outputs')
    parser.add_argument('--output', '-o', default=None,
                        help='Output HTML file (default: <results_dir>/dashboard.html)')
    parser.add_argument('--title', '-t', default=None,
                        help='Dashboard title (default: directory name)')
    parser.add_argument('--conditions', '-c', default=None,
                        help='Conditions JSON file mapping condition names to wells')
    parser.add_argument('--barcode-map', default=None,
                        help='Barcode plate map CSV (default: barcodes/map.csv relative to script)')
    args = parser.parse_args()

    if args.output is None:
        args.output = os.path.join(args.results_dir, 'dashboard.html')
    if args.title is None:
        args.title = os.path.basename(os.path.normpath(args.results_dir))
    if args.barcode_map is None and args.conditions is not None:
        # Default barcode map location relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        args.barcode_map = os.path.join(script_dir, '..', '..', 'barcodes', 'map.csv')

    generate_dashboard(args.results_dir, title=args.title, output_path=args.output,
                       conditions_json=args.conditions, barcode_map=args.barcode_map)


if __name__ == '__main__':
    main()
