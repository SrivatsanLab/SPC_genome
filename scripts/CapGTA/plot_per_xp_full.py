#!/usr/bin/env python3
"""Run the per-XP mate-offset analysis + assignment-breakdown on the full
bidirectional pipeline (2402 cells).

Walks:
  - FASTQ R1 of each cell → read-name → XP tag lookup
  - fc-annotated BAMs → per-read XS/XT tags + R1/R2 mate positions

Writes:
  r2a_mate_offsets.csv              — one row per gene-assigned pair
  xp_assignment_breakdown.tsv       — Assigned/NoFeatures/Overlapping/Ambiguity by XP
  offset_from_gene_ends_<XP>.png    — standard QC plot per XP
  offset_from_gene_ends_all_xp.png  — 4-row composite
"""
from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


XP_PAT = re.compile(r'\bXP:Z:(\w+)')
UMI_IN_NAME_PAT = re.compile(r'_UMI_[ACGTN]+$')


def build_gene_lookup(counts_raw: str) -> dict:
    df = pd.read_csv(counts_raw, sep='\t', comment='#').rename(columns={'Geneid': 'gene_id'})
    out = {}
    for _, row in df.iterrows():
        starts = [int(s) for s in str(row['Start']).split(';')]
        ends   = [int(s) for s in str(row['End']).split(';')]
        strand = str(row['Strand']).split(';')[0]
        out[row['gene_id']] = (min(starts), max(ends), strand)
    return out


def _strip_umi(name: str) -> str:
    return UMI_IN_NAME_PAT.sub('', name)


def load_xp_for_cell(fastq_path: Path) -> dict:
    xp = {}
    if not fastq_path.exists():
        return xp
    with gzip.open(fastq_path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0: continue
            m = XP_PAT.search(line)
            if m is None: continue
            name = line.split()[0][1:]
            xp[name] = m.group(1)
    return xp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq-root', required=True)
    ap.add_argument('--cells-list', required=True)
    ap.add_argument('--counts-raw', required=True)
    ap.add_argument('--anno-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--max-offset', type=int, default=3000)
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    print('Loading gene coords...', flush=True)
    gene_lookup = build_gene_lookup(args.counts_raw)

    cells = [l.strip() for l in open(args.cells_list) if l.strip()]
    anno_dir = Path(args.anno_dir)
    fastq_root = Path(args.fastq_root)

    # Breakdown counters: xp -> {Assigned, NoFeatures, Overlapping_Length, Ambiguity, Unmapped, Other}
    xs_cats = ['Assigned', 'Unassigned_NoFeatures',
               'Unassigned_Overlapping_Length', 'Unassigned_Ambiguity',
               'Unassigned_Unmapped', 'other']
    xp_xs_counts = {xp: Counter() for xp in ('R1T', 'R1A', 'R2A', 'R2T')}
    xp_total = {xp: 0 for xp in xp_xs_counts}

    rows = []
    n_cells_done = 0
    for cell in cells:
        udi = '_'.join(cell.split('_')[-2:])
        fq_r1 = fastq_root / udi / f'{cell}_R1.polyAT_trimmed.fastq.gz'
        bam_path = anno_dir / f'{cell}.bam.featureCounts.bam'
        if not bam_path.exists():
            continue
        xp_lookup = load_xp_for_cell(fq_r1)
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception:
            continue

        pair_data = defaultdict(dict)
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            raw_name = _strip_umi(read.query_name)
            xp = xp_lookup.get(raw_name)
            if xp is None:
                continue
            try:
                xs = read.get_tag('XS')
            except KeyError:
                xs = 'other'
            # Pair-level tally — only R1 (so each pair counted once per fc assignment)
            if read.is_read1:
                xp_total[xp] = xp_total.get(xp, 0) + 1
                xp_xs_counts.setdefault(xp, Counter())[xs if xs in xs_cats else 'other'] += 1
            if read.is_unmapped:
                continue
            d = pair_data[raw_name]
            d['xp'] = xp
            d['chrom'] = read.reference_name
            if xs == 'Assigned':
                try:
                    d['gene_id'] = read.get_tag('XT')
                except KeyError:
                    pass
            if read.is_read1:
                d['r1_start'] = read.reference_start + 1
                d['r1_end']   = read.reference_end
            elif read.is_read2:
                d['r2_start'] = read.reference_start + 1
                d['r2_end']   = read.reference_end
        bam.close()

        for name, d in pair_data.items():
            if 'gene_id' not in d or 'r1_start' not in d or 'r2_start' not in d:
                continue
            gl = gene_lookup.get(d['gene_id'])
            if gl is None:
                continue
            gs, ge, strand = gl
            if strand == '+':
                r1_off_5p = d['r1_start'] - gs
                r1_off_3p = d['r1_end']   - ge
                r2_off_5p = d['r2_start'] - gs
                r2_off_3p = d['r2_end']   - ge
            else:
                r1_off_5p = ge - d['r1_end']
                r1_off_3p = gs - d['r1_start']
                r2_off_5p = ge - d['r2_end']
                r2_off_3p = gs - d['r2_start']
            rows.append({
                'xp': d['xp'], 'gene_id': d['gene_id'], 'gene_strand': strand,
                'r1_off_5p': r1_off_5p, 'r1_off_3p': r1_off_3p,
                'r2_off_5p': r2_off_5p, 'r2_off_3p': r2_off_3p,
            })

        n_cells_done += 1
        if n_cells_done % 100 == 0:
            print(f'  [{n_cells_done}] offsets_so_far={len(rows):,}', flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out / 'r2a_mate_offsets.csv', index=False)
    print(f'\nTotal gene-assigned pairs with XP: {len(df):,}')
    print(df.xp.value_counts().to_string())

    # Assignment breakdown
    print('\n=== featureCounts assignment breakdown per XP ===')
    bd_rows = []
    for xp in xp_xs_counts:
        tot = xp_total[xp]
        if tot == 0: continue
        print(f'\n{xp} (n_pairs={tot:,}):')
        for cat in xs_cats:
            n = xp_xs_counts[xp].get(cat, 0)
            frac = n / tot
            print(f'  {cat:<32s}: {n:>10,}  ({frac:.2%})')
            bd_rows.append({'xp': xp, 'category': cat, 'n': n, 'frac': frac})
    pd.DataFrame(bd_rows).to_csv(out / 'xp_assignment_breakdown.tsv', sep='\t', index=False)

    # Standard-QC-style plot per XP (mate at 5' + primer-side at 3' end)
    XP_SIDES = {
        'R1T': {'mate_5p': 'r2_off_5p', 'primer_3p': 'r1_off_3p',
                'mate_lbl': 'R2', 'primer_lbl': 'R1'},
        'R1A': {'mate_5p': 'r2_off_5p', 'primer_3p': 'r1_off_3p',
                'mate_lbl': 'R2', 'primer_lbl': 'R1'},
        'R2A': {'mate_5p': 'r1_off_5p', 'primer_3p': 'r2_off_3p',
                'mate_lbl': 'R1', 'primer_lbl': 'R2'},
        'R2T': {'mate_5p': 'r1_off_5p', 'primer_3p': 'r2_off_3p',
                'mate_lbl': 'R1', 'primer_lbl': 'R2'},
    }
    xp_order = ['R1T', 'R1A', 'R2A', 'R2T']
    xp_colors = {'R1T': '#2a7', 'R1A': '#7a2', 'R2A': '#a27', 'R2T': '#27a'}
    M = args.max_offset

    def plot_xp(ax_tl, ax_tr, ax_bl, ax_br, sub, xp, color):
        s = XP_SIDES[xp]
        n = len(sub)
        bins_full = np.linspace(-M, M, 121)
        # top-left: mate vs 5' end (full)
        ax_tl.hist(sub[s['mate_5p']].clip(-M, M), bins=bins_full, color=color, edgecolor='white')
        ax_tl.axvline(0, color='k', ls='--', alpha=0.5, label="gene 5' end (TSS)")
        ax_tl.set_title(f'XP={xp}  mate ({s["mate_lbl"]}) vs 5\' end  (n={n:,})')
        ax_tl.set_xlabel(f'{s["mate_lbl"]} start offset from gene 5\' end (bp)')
        ax_tl.set_ylabel(f'# {s["mate_lbl"]} reads')
        ax_tl.legend(loc='upper right')
        # top-right: zoom
        ax_tr.hist(sub[s['mate_5p']].clip(-200, 2000),
                   bins=np.linspace(-200, 2000, 111), color=color, edgecolor='white')
        ax_tr.axvline(0, color='k', ls='--', alpha=0.5)
        ax_tr.set_title(f'XP={xp}  mate ({s["mate_lbl"]}) vs 5\' end — Zoom')
        ax_tr.set_xlabel(f'{s["mate_lbl"]} start offset from gene 5\' end (bp)')
        ax_tr.set_ylabel(f'# {s["mate_lbl"]} reads')
        # bottom-left: primer-side vs 3' end
        ax_bl.hist(sub[s['primer_3p']].clip(-M, M), bins=bins_full, color=color, edgecolor='white')
        ax_bl.axvline(0, color='k', ls='--', alpha=0.5, label="gene 3' end (polyA site)")
        ax_bl.set_title(f'XP={xp}  primer-side ({s["primer_lbl"]}) vs 3\' end')
        ax_bl.set_xlabel(f'{s["primer_lbl"]} offset from gene 3\' end (bp)')
        ax_bl.set_ylabel(f'# {s["primer_lbl"]} reads')
        ax_bl.legend(loc='upper right')
        # bottom-right: zoom
        ax_br.hist(sub[s['primer_3p']].clip(-2000, 500),
                   bins=np.linspace(-2000, 500, 111), color=color, edgecolor='white')
        ax_br.axvline(0, color='k', ls='--', alpha=0.5)
        ax_br.set_title(f'XP={xp}  primer-side ({s["primer_lbl"]}) vs 3\' end — Zoom')
        ax_br.set_xlabel(f'{s["primer_lbl"]} offset from gene 3\' end (bp)')
        ax_br.set_ylabel(f'# {s["primer_lbl"]} reads')

    for xp in xp_order:
        sub = df[df.xp == xp]
        if not len(sub): continue
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        plot_xp(axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1],
                sub, xp, xp_colors[xp])
        plt.tight_layout()
        fp = out / f'offset_from_gene_ends_{xp}.png'
        plt.savefig(fp, dpi=140, bbox_inches='tight')
        plt.close(fig)
        print(f'Wrote {fp}')


if __name__ == '__main__':
    raise SystemExit(main())
