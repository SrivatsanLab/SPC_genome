#!/usr/bin/env python
"""Knee plot of per-barcode E. coli-mapped read counts.

Input TSV: header line `barcode\tecoli_mapped_reads`, sorted descending by count.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="TSV: barcode<TAB>ecoli_mapped_reads")
    ap.add_argument("--sample", required=True, help="Sample name (used in plot title)")
    ap.add_argument("--output", required=True, help="Output PNG path")
    args = ap.parse_args()

    df = pd.read_csv(args.counts, sep="\t")
    if df.empty:
        sys.stderr.write(f"[plot_ecoli_kneeplot] no rows in {args.counts}; skipping plot\n")
        Path(args.output).touch()
        return

    counts = df["ecoli_mapped_reads"].sort_values(ascending=False).to_numpy()
    rank = np.arange(1, counts.size + 1)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(rank, counts, lw=1.2, color="#1f77b4")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("E. coli-mapped reads")
    ax.set_title(f"{args.sample}  |  E. coli reads per barcode\n"
                 f"{counts.size:,} barcodes, {counts.sum():,} mapped reads")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
