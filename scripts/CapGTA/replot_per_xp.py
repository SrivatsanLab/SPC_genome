#!/usr/bin/env python3
"""Re-plot the per-XP offset_from_gene_ends.png panels using FILTERING instead
of clipping (matches qc_primer_trimmed.py behavior)."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


XP_SIDES = {
    'R1T': ('r2_off_5p', 'r1_off_3p', 'R2', 'R1'),
    'R1A': ('r2_off_5p', 'r1_off_3p', 'R2', 'R1'),
    'R2A': ('r1_off_5p', 'r2_off_3p', 'R1', 'R2'),
    'R2T': ('r1_off_5p', 'r2_off_3p', 'R1', 'R2'),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--offsets-csv', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--max-offset', type=int, default=3000)
    args = ap.parse_args()

    df = pd.read_csv(args.offsets_csv)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    M = args.max_offset

    xp_order = ['R1T', 'R1A', 'R2A', 'R2T']
    xp_colors = {'R1T': '#2a7', 'R1A': '#7a2', 'R2A': '#a27', 'R2T': '#27a'}

    for xp in xp_order:
        sub = df[df.xp == xp]
        if not len(sub): continue
        mate_col, prim_col, mate_lbl, prim_lbl = XP_SIDES[xp]
        c = xp_colors[xp]

        mate = sub[mate_col].to_numpy()
        prim = sub[prim_col].to_numpy()

        # Filter: match the original QC convention (drop |offset| > M)
        mate_full = mate[np.abs(mate) <= M]
        prim_full = prim[np.abs(prim) <= M]
        # Zoom windows: match the original QC too (R2@5' zoom was -200..+2000; R1@3' zoom was -2000..+500)
        mate_zoom = mate[(mate >= -200) & (mate <= 2000)]
        prim_zoom = prim[(prim >= -2000) & (prim <= 500)]

        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        bins_full = np.linspace(-M, M, 121)
        # top-left: mate vs 5' end — full
        axes[0, 0].hist(mate_full, bins=bins_full, color=c, edgecolor='white')
        axes[0, 0].axvline(0, color='k', ls='--', alpha=0.5, label="gene 5' end (TSS)")
        axes[0, 0].set_title(f"XP={xp}  mate ({mate_lbl}) vs 5' end  "
                             f"(n={len(mate_full):,} of {len(sub):,})")
        axes[0, 0].set_xlabel(f"{mate_lbl} start offset from gene 5' end (bp)")
        axes[0, 0].set_ylabel(f'# {mate_lbl} reads')
        axes[0, 0].legend(loc='upper right')
        # top-right: zoom
        axes[0, 1].hist(mate_zoom, bins=np.linspace(-200, 2000, 111),
                        color=c, edgecolor='white')
        axes[0, 1].axvline(0, color='k', ls='--', alpha=0.5)
        axes[0, 1].set_title(f"XP={xp}  mate ({mate_lbl}) vs 5' end — Zoom")
        axes[0, 1].set_xlabel(f"{mate_lbl} start offset from gene 5' end (bp)")
        axes[0, 1].set_ylabel(f'# {mate_lbl} reads')
        # bottom-left: primer-side vs 3' end — full
        axes[1, 0].hist(prim_full, bins=bins_full, color=c, edgecolor='white')
        axes[1, 0].axvline(0, color='k', ls='--', alpha=0.5, label="gene 3' end (polyA site)")
        axes[1, 0].set_title(f"XP={xp}  primer-side ({prim_lbl}) vs 3' end  "
                             f"(n={len(prim_full):,} of {len(sub):,})")
        axes[1, 0].set_xlabel(f"{prim_lbl} offset from gene 3' end (bp)")
        axes[1, 0].set_ylabel(f'# {prim_lbl} reads')
        axes[1, 0].legend(loc='upper right')
        # bottom-right: zoom
        axes[1, 1].hist(prim_zoom, bins=np.linspace(-2000, 500, 111),
                        color=c, edgecolor='white')
        axes[1, 1].axvline(0, color='k', ls='--', alpha=0.5)
        axes[1, 1].set_title(f"XP={xp}  primer-side ({prim_lbl}) vs 3' end — Zoom")
        axes[1, 1].set_xlabel(f"{prim_lbl} offset from gene 3' end (bp)")
        axes[1, 1].set_ylabel(f'# {prim_lbl} reads')

        plt.tight_layout()
        fp = out / f'offset_from_gene_ends_{xp}.png'
        plt.savefig(fp, dpi=140, bbox_inches='tight')
        plt.close(fig)
        n_dropped_mate = len(sub) - len(mate_full)
        n_dropped_prim = len(sub) - len(prim_full)
        print(f'Wrote {fp}')
        print(f'  dropped |mate 5p| > {M}: {n_dropped_mate:,}/{len(sub):,} '
              f'({n_dropped_mate/len(sub):.1%})')
        print(f'  dropped |primer 3p| > {M}: {n_dropped_prim:,}/{len(sub):,} '
              f'({n_dropped_prim/len(sub):.1%})')


if __name__ == '__main__':
    raise SystemExit(main())
