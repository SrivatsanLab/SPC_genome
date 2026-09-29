#!/usr/bin/env python3
"""Audit how build_panels' filter chain drops variants, per worm.

For a given config, load the joint AnnData + apply QC, then step through the
panels filter cascade (matching phylo.pp.panels.build_panels) and count what
each stage drops. Also cross-check against the notebook's per-worm VAF spectrum
(``spc.tl.compute_bulk_vaf`` at target_dp = median).

Usage
-----
  panel_filter_audit.py --config scripts/phylogeny/configs/tct_kept.yaml \
                        --out-dir results/worm6_final/DNA_analysis/phylogeny/panels/tct_kept/audit
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
sys.path.insert(0, str(REPO / "scripts" / "phylogeny"))

from phylo.pp.panels import (  # noqa: E402
    _compute_per_worm_vaf,
    _n_carriers_per_worm,
    _pool_ad_dp_per_worm,
    _sub_types_from_trinuc,
    canonical_sub_type,
)
from phylo.pp.qc import apply_qc  # noqa: E402
from phylo.utils.config import load_config  # noqa: E402


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    return ap.parse_args()


def compute_notebook_vaf(adata: ad.AnnData) -> pd.DataFrame:
    """Recreate the notebook's per-worm bulk VAF (hypergeometric at per-worm median DP).

    Matches ``worm6_somatic_panels.ipynb`` cells 5–8: for each donor, subset
    cells → ``spc.tl.compute_bulk_vaf(sub, target_dp=median(bulk_dp))``.
    """
    from cellspec.utils.context import compute_vaf
    donors = sorted(adata.obs["donor"].astype(str).unique())
    out = pd.DataFrame(index=adata.var_names, columns=donors, dtype=np.float32)
    for d in donors:
        mask = (adata.obs["donor"].astype(str) == d).to_numpy()
        AD = adata.layers["AD"][mask]
        DP = adata.layers["DP"][mask]
        if issparse(AD):
            AD = AD.toarray()
        if issparse(DP):
            DP = DP.toarray()
        ad_sum = AD.sum(axis=0).astype(np.int64)
        dp_sum = DP.sum(axis=0).astype(np.int64)
        nz = dp_sum[dp_sum > 0]
        med = int(np.median(nz)) if nz.size else 1
        out[d] = compute_vaf(ad_sum, dp_sum, target_dp=med)
    return out


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(args.config)

    adata_path = Path(cfg["inputs"]["adata"])
    if not adata_path.is_absolute():
        adata_path = (args.config.parent / adata_path).resolve()
    print(f"[audit] loading {adata_path}")
    adata = ad.read_h5ad(adata_path)
    print(f"[audit] loaded {adata.n_obs:,} cells × {adata.n_vars:,} variants")

    apply_qc(adata, cfg["qc"]["predicates"], inplace=True)
    n_qc = int(adata.obs["qc_pass"].sum())
    print(f"[audit] qc_pass: {n_qc}/{adata.n_obs} cells")

    # --- config knobs --------------------------------------------------------
    p_cfg = cfg["panels"]
    band_raw = p_cfg["somatic_vaf_band"]
    lo = -np.inf if band_raw[0] is None else float(band_raw[0])
    hi = np.inf if band_raw[1] is None else float(band_raw[1])
    other_max = float(p_cfg["other_worm_max_vaf"])
    min_carriers = int(p_cfg["min_carriers"])
    bl_min_worms = int(p_cfg["blacklist_min_worms"])
    bl_vaf_min = float(p_cfg["blacklist_vaf_min"])
    min_dp_per_worm = int(p_cfg.get("min_bulk_dp_per_worm", 5))
    excluded = set(p_cfg.get("exclude_mutation_types") or [])
    excluded = {canonical_sub_type(s.split(">")[0], s.split(">")[1]) for s in excluded}
    retain_ctx = set(str(s).upper() for s in (p_cfg.get("retain_trinuc_contexts") or []))

    # --- filter 1: substitution type ----------------------------------------
    n_before = adata.n_vars
    if excluded:
        sub_types = _sub_types_from_trinuc(adata.var["anc"], adata.var["der"])
        type_pass = ~pd.Series(sub_types).isin(excluded).to_numpy()
        if retain_ctx:
            pairs = (adata.var["anc"].astype(str).str.upper() + ">"
                     + adata.var["der"].astype(str).str.upper()).to_numpy()
            type_pass |= np.isin(pairs, list(retain_ctx))
        adata = adata[:, type_pass].copy()
        print(f"[audit] filter 1 (sub-type + retain): {n_before:,} → {adata.n_vars:,}")

    # --- pool AD/DP per worm; compute VAFs ---------------------------------
    print("[audit] pooling AD/DP per worm ...")
    AD_w, DP_w = _pool_ad_dp_per_worm(adata)
    VAF_w = _compute_per_worm_vaf(adata, AD_w, DP_w, method="hypergeometric",
                                  target_dp=p_cfg.get("target_dp"))
    N_w = _n_carriers_per_worm(adata)

    shared = (VAF_w >= bl_vaf_min).sum(axis=1)
    bl_mask = (shared >= bl_min_worms).to_numpy()
    print(f"[audit] blacklist (VAF≥{bl_vaf_min} in ≥{bl_min_worms} worms): "
          f"{int(bl_mask.sum()):,} / {len(bl_mask):,}")

    # --- per-worm cascade ---------------------------------------------------
    audit_rows = []
    worms = list(VAF_w.columns)
    for w in worms:
        own_vaf = VAF_w[w].to_numpy()
        own_dp = DP_w[w].to_numpy()
        others = [c for c in worms if c != w]
        other_max_vaf = VAF_w[others].to_numpy().max(axis=1)
        own_ncarriers = N_w[w].to_numpy()

        m_band = (own_vaf >= lo) & (own_vaf < hi)
        m_other = other_max_vaf <= other_max
        m_dp = own_dp >= min_dp_per_worm
        m_carriers = own_ncarriers >= min_carriers
        m_blacklist = ~bl_mask

        row = {
            "worm": w,
            "n_variants_after_type_filter": int(adata.n_vars),
            "own_vaf_in_band": int(m_band.sum()),
            "and_other_vaf_le_0.005": int((m_band & m_other).sum()),
            "and_own_dp_ge_5": int((m_band & m_other & m_dp).sum()),
            "and_min_carriers_ge_2": int((m_band & m_other & m_dp & m_carriers).sum()),
            "and_not_blacklisted (FINAL)": int((m_band & m_other & m_dp & m_carriers & m_blacklist).sum()),
        }
        # Also count how many would-be-panel-vars would be admitted if we
        # relaxed each filter individually (from the intersection of the others).
        row["would_add_if_other_max_relaxed"] = int(
            ((m_band & m_dp & m_carriers & m_blacklist) & ~m_other).sum())
        row["would_add_if_min_carriers_relaxed"] = int(
            ((m_band & m_other & m_dp & m_blacklist) & ~m_carriers).sum())
        row["would_add_if_blacklist_relaxed"] = int(
            ((m_band & m_other & m_dp & m_carriers) & ~m_blacklist).sum())
        row["would_add_if_min_dp_relaxed"] = int(
            ((m_band & m_other & m_carriers & m_blacklist) & ~m_dp).sum())
        audit_rows.append(row)

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(args.out_dir / "filter_cascade.tsv", sep="\t", index=False)
    print(f"[audit] wrote filter_cascade.tsv")

    # --- notebook-style VAF spectrum overlaid with panel-membership --------
    print("[audit] computing notebook-style per-worm VAF ...")
    nb_vaf = compute_notebook_vaf(adata)
    nb_vaf.to_csv(args.out_dir / "notebook_vaf.tsv", sep="\t")

    # For each worm, plot histograms of nb_vaf coloured by whether the variant
    # eventually made it into the panel. This lets you visually see which
    # subset of the "somatic band" gets dropped and where.
    print("[audit] plotting per-worm VAF overlays ...")
    n = len(worms)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 2.5 * nrows))
    axes_flat = axes.ravel() if n > 1 else [axes]
    for i, w in enumerate(worms):
        ax = axes_flat[i]
        v_all = nb_vaf[w].to_numpy()
        # Build the "in-panel" mask same as the row above
        own_vaf = VAF_w[w].to_numpy()
        own_dp = DP_w[w].to_numpy()
        others = [c for c in worms if c != w]
        other_max_vaf = VAF_w[others].to_numpy().max(axis=1)
        own_ncarriers = N_w[w].to_numpy()
        in_panel = (
            (own_vaf >= lo) & (own_vaf < hi)
            & (other_max_vaf <= other_max)
            & (own_dp >= min_dp_per_worm)
            & (own_ncarriers >= min_carriers)
            & (~bl_mask)
        )
        # Everything in the type-filtered joint set with per-worm VAF > 0
        m_pos = v_all > 0
        v_in = v_all[in_panel & m_pos]
        v_out = v_all[~in_panel & m_pos]
        bins = np.linspace(0, 0.5, 51)
        ax.hist(v_out, bins=bins, color="0.6", alpha=0.75, label=f"dropped n={len(v_out)}")
        ax.hist(v_in,  bins=bins, color="tab:blue", alpha=0.8, label=f"in panel n={len(v_in)}")
        ax.set_title(w, fontsize=9)
        ax.set_xlim(0, 0.5)
        ax.set_xlabel("notebook bulk_vaf", fontsize=8)
        ax.set_ylabel("count", fontsize=8)
        ax.set_yscale("log")
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7, loc="upper right")
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle("Notebook-style bulk VAF spectrum: dropped (grey) vs in-panel (blue)",
                 fontsize=11)
    fig.tight_layout()
    for suf in (".png", ".svg"):
        fig.savefig(args.out_dir / f"vaf_dropped_vs_kept{suf}", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[audit] wrote vaf_dropped_vs_kept.{{png,svg}}")


if __name__ == "__main__":
    main()
