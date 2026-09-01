"""For pseudobulk variants below a pooled-VAF threshold, count how many worms
carry the variant (GT in {0/1, 1/1}) and plot the distribution.

Purpose: get a first look at how many segregating (non-fixed-vs-N2) variants
are private to one worm vs shared across many. Fixed strain differences from
the N2 reference sit near pooled VAF ≈ 1 and are excluded by the < 0.8 cut.

Usage:
    python variant_sharing_histogram.py <regenotyped_tsv> <output_dir> \\
        [--vaf-cut 0.8]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("regenotyped_tsv")
    p.add_argument("output_dir")
    p.add_argument("--vaf-cut", type=float, default=0.8,
                   help="Drop sites with pooled VAF >= this (default 0.8)")
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.regenotyped_tsv}…")
    tab = pd.read_csv(args.regenotyped_tsv, sep="\t")

    worms = sorted({c.rsplit("_", 1)[0] for c in tab.columns if c.endswith("_gt")})
    print(f"  {len(tab)} sites, {len(worms)} worms")

    dp = tab[[f"{w}_dp" for w in worms]].to_numpy()
    ad = tab[[f"{w}_ad" for w in worms]].to_numpy()
    gt = tab[[f"{w}_gt" for w in worms]].to_numpy()

    total_dp = dp.sum(axis=1)
    total_ad = ad.sum(axis=1)
    pooled_vaf = np.divide(total_ad, total_dp,
                           out=np.full_like(total_ad, np.nan, dtype=float),
                           where=total_dp > 0)

    n_carriers = ((gt == "0/1") | (gt == "1/1")).sum(axis=1)

    keep = (pooled_vaf < args.vaf_cut) & (n_carriers > 0)
    n_dropped_fixed = int((pooled_vaf >= args.vaf_cut).sum())
    n_dropped_none = int(((pooled_vaf < args.vaf_cut) & (n_carriers == 0)).sum())
    print(f"  dropped {n_dropped_fixed} sites with pooled VAF >= {args.vaf_cut}")
    print(f"  dropped {n_dropped_none} sites with 0 carriers among the rest")
    print(f"  keeping  {int(keep.sum())} sites for sharing histogram")

    carriers = n_carriers[keep]

    # Counts at every integer 1..N so the CSV has explicit zero bins if any.
    n = len(worms)
    counts = pd.Series(0, index=range(1, n + 1), dtype=int)
    counts.update(pd.Series(carriers).value_counts())
    counts.index.name = "n_worms"
    counts.name = "n_variants"
    csv_path = out / "variant_sharing_counts.csv"
    counts.to_csv(csv_path)
    print(f"\n  wrote {csv_path}")
    print(counts.to_string())

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(counts.index, counts.values, color="steelblue", edgecolor="black")
    ax.set_xticks(range(1, n + 1))
    ax.set_xlabel("number of worms carrying the variant (GT ∈ {0/1, 1/1})")
    ax.set_ylabel("number of variants")
    ax.set_title(f"Variant sharing across worms (pooled VAF < {args.vaf_cut})\n"
                 f"n={int(keep.sum())} sites, {n} worms")
    for x, y in zip(counts.index, counts.values):
        if y > 0:
            ax.text(x, y, str(int(y)), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig_path = out / "variant_sharing_histogram.png"
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)
    print(f"  wrote {fig_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
