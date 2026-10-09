#!/usr/bin/env python3
"""Standard-QC-style plot per XP primer type, with the mate correctly identified.

Original QC plot shows:
  Top: MATE (non-primer read) start offset from gene 5' end
  Bottom: PRIMER-ADJACENT read offset from gene 3' end

For bidirectional-extracted reads:
  R1T, R1A → primer on R1 → mate = R2 → top: R2 offset from 5' end
                                      → bottom: R1 offset from 3' end
  R2A, R2T → primer on R2 → mate = R1 → top: R1 offset from 5' end
                                      → bottom: R2 offset from 3' end
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Map XP -> (mate side at 5', primer side at 3')
# For R1T/R1A: primer=R1, mate=R2
# For R2A/R2T: primer=R2, mate=R1
XP_SIDES = {
    'R1T': {'mate_5p': 'r2_off_5p', 'primer_3p': 'r1_off_3p',
            'mate_label': 'R2', 'primer_label': 'R1'},
    'R1A': {'mate_5p': 'r2_off_5p', 'primer_3p': 'r1_off_3p',
            'mate_label': 'R2', 'primer_label': 'R1'},
    'R2A': {'mate_5p': 'r1_off_5p', 'primer_3p': 'r2_off_3p',
            'mate_label': 'R1', 'primer_label': 'R2'},
    'R2T': {'mate_5p': 'r1_off_5p', 'primer_3p': 'r2_off_3p',
            'mate_label': 'R1', 'primer_label': 'R2'},
}


def plot_one_xp(ax_mate_full, ax_mate_zoom, ax_primer_full, ax_primer_zoom,
                sub: pd.DataFrame, xp: str, max_offset: int, color: str):
    n = len(sub)
    sides = XP_SIDES[xp]
    mate_col   = sides['mate_5p']
    primer_col = sides['primer_3p']
    mate_lbl   = sides['mate_label']
    prim_lbl   = sides['primer_label']

    bins_full = np.linspace(-max_offset, max_offset, 121)
    # Mate vs gene 5' end (full)
    ax_mate_full.hist(sub[mate_col].clip(-max_offset, max_offset),
                      bins=bins_full, color=color, edgecolor='white')
    ax_mate_full.axvline(0, color='k', ls='--', alpha=0.5, label="gene 5' end (TSS)")
    ax_mate_full.set_xlabel(f"{mate_lbl} start offset from gene 5' end (bp)")
    ax_mate_full.set_ylabel(f'# {mate_lbl} reads')
    ax_mate_full.set_title(f'XP={xp}  mate ({mate_lbl}) vs 5\' end  (n={n:,})')
    ax_mate_full.legend(loc='upper right')
    # Mate vs 5' end (zoom)
    ax_mate_zoom.hist(sub[mate_col].clip(-200, 2000),
                      bins=np.linspace(-200, 2000, 111), color=color, edgecolor='white')
    ax_mate_zoom.axvline(0, color='k', ls='--', alpha=0.5)
    ax_mate_zoom.set_xlabel(f"{mate_lbl} start offset from gene 5' end (bp)")
    ax_mate_zoom.set_ylabel(f'# {mate_lbl} reads')
    ax_mate_zoom.set_title(f'XP={xp}  mate ({mate_lbl}) vs 5\' end — Zoom')
    # Primer-adjacent read vs gene 3' end (full)
    ax_primer_full.hist(sub[primer_col].clip(-max_offset, max_offset),
                        bins=bins_full, color=color, edgecolor='white')
    ax_primer_full.axvline(0, color='k', ls='--', alpha=0.5, label="gene 3' end (polyA site)")
    ax_primer_full.set_xlabel(f"{prim_lbl} offset from gene 3' end (bp)")
    ax_primer_full.set_ylabel(f'# {prim_lbl} reads')
    ax_primer_full.set_title(f'XP={xp}  primer ({prim_lbl}) vs 3\' end  (n={n:,})')
    ax_primer_full.legend(loc='upper right')
    # Primer-adjacent read vs 3' end (zoom)
    ax_primer_zoom.hist(sub[primer_col].clip(-2000, 500),
                        bins=np.linspace(-2000, 500, 111), color=color, edgecolor='white')
    ax_primer_zoom.axvline(0, color='k', ls='--', alpha=0.5)
    ax_primer_zoom.set_xlabel(f"{prim_lbl} offset from gene 3' end (bp)")
    ax_primer_zoom.set_ylabel(f'# {prim_lbl} reads')
    ax_primer_zoom.set_title(f'XP={xp}  primer ({prim_lbl}) vs 3\' end — Zoom')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--offsets-csv', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--max-offset', type=int, default=3000)
    args = ap.parse_args()

    df = pd.read_csv(args.offsets_csv)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    xp_order = ['R1T', 'R1A', 'R2A', 'R2T']
    xp_colors = {'R1T': '#2a7', 'R1A': '#7a2', 'R2A': '#a27', 'R2T': '#27a'}

    for xp in xp_order:
        sub = df[df.xp == xp]
        if not len(sub): continue
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        plot_one_xp(
            ax_mate_full=axes[0, 0], ax_mate_zoom=axes[0, 1],
            ax_primer_full=axes[1, 0], ax_primer_zoom=axes[1, 1],
            sub=sub, xp=xp, max_offset=args.max_offset, color=xp_colors[xp],
        )
        plt.tight_layout()
        fp = out / f'offset_from_gene_ends_{xp}.png'
        plt.savefig(fp, dpi=140, bbox_inches='tight')
        print(f'Wrote {fp}')
        plt.close(fig)

    # Composite: 4 rows (XP) × 4 cols
    fig, axes = plt.subplots(len(xp_order), 4, figsize=(22, 4.2 * len(xp_order)))
    for i, xp in enumerate(xp_order):
        sub = df[df.xp == xp]
        if not len(sub):
            for j in range(4):
                axes[i, j].text(0.5, 0.5, f'XP={xp}: no data', ha='center', va='center')
            continue
        plot_one_xp(
            ax_mate_full=axes[i, 0], ax_mate_zoom=axes[i, 1],
            ax_primer_full=axes[i, 2], ax_primer_zoom=axes[i, 3],
            sub=sub, xp=xp, max_offset=args.max_offset, color=xp_colors[xp],
        )
    plt.tight_layout()
    plt.savefig(out / 'offset_from_gene_ends_all_xp.png', dpi=120, bbox_inches='tight')
    print(f'Wrote {out}/offset_from_gene_ends_all_xp.png')


if __name__ == '__main__':
    raise SystemExit(main())
