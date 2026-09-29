#!/usr/bin/env python3
"""Coverage-breadth diagnostic to pick top-N cutoff.

For each panel h5ad in a tag directory, plot:
  - breadth spectrum (# cells with DP>0 per variant)
  - VAF (pooled) spectrum of top-N vs full panel, for a sweep of N

Writes: <out_dir>/breadth_summary.tsv, breadth_spectrum.png/svg,
        topN_vaf_preservation.png/svg
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import issparse


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--ns", nargs="+", type=int, default=[50, 100, 200, 500, 1000],
                    help="Top-N sweep values")
    return ap.parse_args()


def compute_stats(h5ad: Path):
    a = ad.read_h5ad(h5ad)
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else a.layers["DP"]
    AD = a.layers["AD"].toarray() if issparse(a.layers["AD"]) else a.layers["AD"]
    n_cov = (DP > 0).sum(axis=0)
    ad_sum, dp_sum = AD.sum(axis=0), DP.sum(axis=0)
    pooled_vaf = np.where(dp_sum > 0, ad_sum / dp_sum, np.nan)
    breadth = n_cov / a.n_obs
    return {
        "worm": a.uns.get("panel_worm", h5ad.stem.replace("worm_", "")),
        "n_cells": a.n_obs, "n_vars": a.n_vars,
        "n_cov": n_cov, "breadth": breadth, "pooled_vaf": pooled_vaf,
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    stats = [compute_stats(p) for p in sorted(args.panels_dir.glob("worm_worm*.h5ad"))]
    if not stats:
        print("no panels"); return

    # Summary table per worm at each N.
    rows = []
    for s in stats:
        for N in args.ns:
            idx = np.argsort(-s["n_cov"], kind="stable")[:min(N, s["n_vars"])]
            n_kept = len(idx)
            med_breadth = float(np.median(s["breadth"][idx]))
            min_breadth = float(np.min(s["breadth"][idx]))
            vaf_med_full = float(np.nanmedian(s["pooled_vaf"]))
            vaf_med_kept = float(np.nanmedian(s["pooled_vaf"][idx]))
            rows.append({
                "worm": s["worm"], "N": N, "n_vars_total": s["n_vars"],
                "n_kept": n_kept, "min_breadth_kept": min_breadth,
                "med_breadth_kept": med_breadth,
                "med_pooled_vaf_full": vaf_med_full,
                "med_pooled_vaf_kept": vaf_med_kept,
            })
    df = pd.DataFrame(rows)
    df.to_csv(args.out_dir / "breadth_summary.tsv", sep="\t", index=False)
    print(f"wrote breadth_summary.tsv")

    # Breadth spectrum, one panel per worm.
    n = len(stats)
    ncols = 4; nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 2.3 * nrows))
    ax_flat = axes.ravel() if n > 1 else [axes]
    for i, s in enumerate(stats):
        ax = ax_flat[i]
        ax.hist(s["breadth"], bins=40, range=(0, 1), color="tab:blue", alpha=0.75)
        ax.set_title(f"{s['worm']}: {s['n_vars']} vars", fontsize=9)
        ax.set_xlabel("coverage breadth", fontsize=8); ax.set_xlim(0, 1)
        ax.set_ylabel("count", fontsize=8); ax.tick_params(labelsize=7)
    for j in range(n, len(ax_flat)):
        ax_flat[j].axis("off")
    fig.suptitle("Coverage breadth per variant", fontsize=11)
    fig.tight_layout()
    for suf in (".png", ".svg"):
        fig.savefig(args.out_dir / f"breadth_spectrum{suf}", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # VAF preservation under top-N cutoffs, per worm.
    Ns_use = args.ns
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 2.3 * nrows))
    ax_flat = axes.ravel() if n > 1 else [axes]
    for i, s in enumerate(stats):
        ax = ax_flat[i]
        v_all = s["pooled_vaf"][~np.isnan(s["pooled_vaf"])]
        bins = np.linspace(0, 0.5, 40)
        ax.hist(v_all, bins=bins, color="0.8", alpha=0.8, label=f"all ({s['n_vars']})", density=True)
        colors = plt.cm.viridis(np.linspace(0, 0.9, len(Ns_use)))
        for N, col in zip(Ns_use, colors):
            idx = np.argsort(-s["n_cov"], kind="stable")[:min(N, s["n_vars"])]
            v = s["pooled_vaf"][idx]
            v = v[~np.isnan(v)]
            ax.hist(v, bins=bins, histtype="step", color=col, lw=1.2,
                    label=f"top {N}", density=True)
        ax.set_title(f"{s['worm']}", fontsize=9)
        ax.set_xlabel("pooled VAF", fontsize=8); ax.set_xlim(0, 0.5)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=6, loc="upper right")
    for j in range(n, len(ax_flat)):
        ax_flat[j].axis("off")
    fig.suptitle("Top-N (by coverage breadth) — VAF spectrum preservation", fontsize=11)
    fig.tight_layout()
    for suf in (".png", ".svg"):
        fig.savefig(args.out_dir / f"topN_vaf_preservation{suf}", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote breadth_spectrum + topN_vaf_preservation png/svg")


if __name__ == "__main__":
    main()
