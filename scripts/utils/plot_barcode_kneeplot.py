#!/usr/bin/env python
"""Generic barcode-rank knee plot from a 2-column counts TSV.

Input TSV: header `barcode\t<count_col>`, sorted descending by count.
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
    ap.add_argument("--counts", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--ylabel", default="Reads per barcode")
    ap.add_argument("--title-suffix", default="")
    args = ap.parse_args()

    df = pd.read_csv(args.counts, sep="\t")
    if df.empty or df.shape[1] < 2:
        sys.stderr.write(f"[plot_barcode_kneeplot] no usable rows in {args.counts}; touching empty output\n")
        Path(args.output).touch()
        return

    counts = df.iloc[:, 1].sort_values(ascending=False).to_numpy()
    rank = np.arange(1, counts.size + 1)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(rank, counts, lw=1.2, color="#1f77b4")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel(args.ylabel)
    title = f"{args.sample}"
    if args.title_suffix:
        title += f"  |  {args.title_suffix}"
    title += f"\n{counts.size:,} barcodes, {counts.sum():,} reads"
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
