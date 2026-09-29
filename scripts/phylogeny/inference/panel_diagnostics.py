#!/usr/bin/env python3
"""D1–D3 panel diagnostics per plan §1.

For each per-worm panel h5ad, compute:

  variants.tsv            per-variant stats: pooled_VAF, carrier_frac, mean_DP,
                          implied_f, and the D3 detection-ratio for artifact ID.
  d1_spectrum.png/svg     carrier_frac + implied_f histograms (D1).
  d2_vaf_spectrum.png/svg pooled VAF spectrum with band overlays (D2).
  d3_ratio_vs_dp.png/svg  carrier_frac / pooled_VAF vs mean_DP, overlaid on
                          the §2.1 depth curve (1.00× at DP=1, 1.30× at 1.5×,
                          1.99× at 3×).

The point of D1–D3: does the panel actually mark clades, or is it dominated by
artifacts / detection-driven noise?

Usage
-----
  panel_diagnostics.py \
      --panels-dir results/worm6_final/DNA_analysis/phylogeny/panels/<tag> \
      --out-dir    results/worm6_final/DNA_analysis/phylogeny/panels/<tag>/diagnostics
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


# Alpha-hap × contamination adjustment from HAPLOTYPE_ANALYSIS_SUMMARY.md §2.
POOLED_VAF_AT_F1 = 0.345         # pooled VAF when clade fraction f = 1
D3_DEPTH_CURVE = [(1.0, 1.00), (1.5, 1.30), (3.0, 1.99)]  # (mean DP, ratio)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels-dir", type=Path, required=True,
                    help="phylogeny/panels/<panel_tag>")
    ap.add_argument("--out-dir", type=Path, required=True)
    return ap.parse_args()


def compute_variant_stats(h5ad_path: Path) -> pd.DataFrame:
    a = ad.read_h5ad(h5ad_path)
    AD = a.layers["AD"].toarray() if issparse(a.layers["AD"]) else a.layers["AD"]
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else a.layers["DP"]
    covered = DP > 0
    n_covered = covered.sum(axis=0)
    n_carrier = (AD >= 1).sum(axis=0)
    ad_sum = AD.sum(axis=0)
    dp_sum = DP.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled_vaf = np.where(dp_sum > 0, ad_sum / dp_sum, np.nan)
        carrier_frac = np.where(n_covered > 0, n_carrier / n_covered, np.nan)
        mean_dp = np.where(n_covered > 0, dp_sum / n_covered, np.nan)
        implied_f = pooled_vaf / POOLED_VAF_AT_F1
        ratio = np.where((pooled_vaf > 0), carrier_frac / pooled_vaf, np.nan)
    return pd.DataFrame({
        "variant_id": a.var_names,
        "chrom": a.var["chrom"].to_numpy(),
        "pos": a.var["pos"].to_numpy(),
        "trinuc_type": a.var["trinuc_type"].astype(str).to_numpy() if "trinuc_type" in a.var else "",
        "n_covered": n_covered,
        "n_carriers": n_carrier,
        "pooled_vaf": pooled_vaf,
        "carrier_frac": carrier_frac,
        "mean_dp": mean_dp,
        "implied_f": implied_f,
        "detection_ratio": ratio,
        "worm": a.uns.get("panel_worm", h5ad_path.stem.replace("worm_", "")),
        "worm_n_cells": a.n_obs,
    })


def plot_d1(df: pd.DataFrame, out_stem: Path) -> None:
    """Per-worm carrier_frac + implied_f spectrum."""
    worms = sorted(df["worm"].unique())
    n = len(worms)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.4 * nrows),
                             sharex=False, sharey=False)
    axes_flat = axes.ravel() if n > 1 else [axes]
    for i, w in enumerate(worms):
        ax = axes_flat[i]
        sub = df[df["worm"] == w]
        ax.hist(sub["carrier_frac"].dropna(), bins=40, range=(0, 1),
                color="tab:blue", alpha=0.7, label="carrier_frac")
        ax.hist(sub["implied_f"].dropna().clip(0, 1), bins=40, range=(0, 1),
                color="tab:orange", alpha=0.5, label="implied_f")
        n_cells = int(sub["worm_n_cells"].iloc[0])
        n_vars = len(sub)
        ax.set_title(f"{w}: {n_vars} var × {n_cells} cells", fontsize=9)
        ax.set_xlim(0, 1)
        ax.set_xlabel("frac / clade fraction", fontsize=8)
        ax.set_ylabel("count", fontsize=8)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7, loc="upper right")
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle("D1: per-variant carrier_frac (blue) and implied clade fraction (orange)",
                 fontsize=11)
    fig.tight_layout()
    for suf in (".png", ".svg"):
        fig.savefig(out_stem.with_suffix(suf), dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_d2(df: pd.DataFrame, out_stem: Path) -> None:
    """Per-worm pooled_VAF spectrum with reference lines."""
    worms = sorted(df["worm"].unique())
    n = len(worms)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.4 * nrows))
    axes_flat = axes.ravel() if n > 1 else [axes]
    for i, w in enumerate(worms):
        ax = axes_flat[i]
        sub = df[df["worm"] == w]
        ax.hist(sub["pooled_vaf"].dropna(), bins=50, range=(0, 0.5), color="0.4")
        # Reference lines
        ax.axvline(0.20, color="tab:red", ls="--", lw=1, label="0.20 (f≈0.58)")
        ax.axvline(POOLED_VAF_AT_F1, color="tab:orange", ls=":", lw=1, label="0.345 (f=1)")
        # Fraction of panel above 0.20
        n_above = int((sub["pooled_vaf"] > 0.20).sum())
        pct = 100 * n_above / len(sub) if len(sub) else 0
        ax.set_title(f"{w}: {pct:.0f}% > 0.20", fontsize=9)
        ax.set_xlim(0, 0.5)
        ax.set_xlabel("pooled VAF", fontsize=8)
        ax.set_ylabel("count", fontsize=8)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7, loc="upper right")
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle("D2: pooled VAF spectrum (red = f≈0.58 threshold; orange = f=1)",
                 fontsize=11)
    fig.tight_layout()
    for suf in (".png", ".svg"):
        fig.savefig(out_stem.with_suffix(suf), dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_d3(df: pd.DataFrame, out_stem: Path) -> None:
    """Detection ratio (carrier_frac / pooled_VAF) vs mean_DP with §2.1 curve."""
    worms = sorted(df["worm"].unique())
    n = len(worms)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.4 * nrows))
    axes_flat = axes.ravel() if n > 1 else [axes]
    curve_x = np.array([p[0] for p in D3_DEPTH_CURVE])
    curve_y = np.array([p[1] for p in D3_DEPTH_CURVE])
    for i, w in enumerate(worms):
        ax = axes_flat[i]
        sub = df[df["worm"] == w]
        mask = sub["detection_ratio"].notna() & sub["mean_dp"].notna()
        s = sub[mask]
        ax.scatter(s["mean_dp"], s["detection_ratio"], s=4, alpha=0.35, color="0.35")
        ax.plot(curve_x, curve_y, "o-", color="tab:red", label="theory (§2.1)")
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 5)
        ax.set_title(f"{w}", fontsize=9)
        ax.set_xlabel("mean DP over covered cells", fontsize=8)
        ax.set_ylabel("carrier_frac / pooled_VAF", fontsize=8)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7, loc="upper right")
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle("D3: detection ratio vs mean depth (red = calibrated curve)",
                 fontsize=11)
    fig.tight_layout()
    for suf in (".png", ".svg"):
        fig.savefig(out_stem.with_suffix(suf), dpi=140, bbox_inches="tight")
    plt.close(fig)


def summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """Per-worm one-line summary."""
    rows = []
    for w, sub in df.groupby("worm"):
        rows.append({
            "worm": w,
            "n_cells": int(sub["worm_n_cells"].iloc[0]),
            "n_variants": len(sub),
            "med_carrier_frac": float(sub["carrier_frac"].median()),
            "med_pooled_vaf": float(sub["pooled_vaf"].median()),
            "med_implied_f": float(sub["implied_f"].median()),
            "frac_pooled_vaf_gt_0.20": float((sub["pooled_vaf"] > 0.20).mean()),
            "frac_implied_f_gt_0.5": float((sub["implied_f"] > 0.5).mean()),
            "med_mean_dp": float(sub["mean_dp"].median()),
            "med_detection_ratio": float(sub["detection_ratio"].median()),
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    h5ads = sorted(args.panels_dir.glob("worm_worm*.h5ad"))
    if not h5ads:
        print(f"[panel_diagnostics] no worm_*.h5ad under {args.panels_dir}"); return

    frames = []
    for p in h5ads:
        print(f"  computing stats for {p.name}")
        frames.append(compute_variant_stats(p))
    df = pd.concat(frames, ignore_index=True)

    df.to_csv(args.out_dir / "variants.tsv", sep="\t", index=False)
    summary_table(df).to_csv(args.out_dir / "summary.tsv", sep="\t", index=False)
    print(f"[panel_diagnostics] wrote variants.tsv ({len(df)} rows), summary.tsv")

    plot_d1(df, args.out_dir / "d1_spectrum")
    plot_d2(df, args.out_dir / "d2_vaf_spectrum")
    plot_d3(df, args.out_dir / "d3_ratio_vs_dp")
    print(f"[panel_diagnostics] wrote d1/d2/d3 png+svg")


if __name__ == "__main__":
    main()
