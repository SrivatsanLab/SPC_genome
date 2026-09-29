#!/usr/bin/env python3
"""Set-based shared/private mutation diagnostics for top-down SVD trees.

For every internal node of every tree.json in the sweep, walks the panel
h5ad and classifies the variants that had ≥1 carrier among this node's
cells into four bins by SET MEMBERSHIP (not distributional inference):

  * shared_across_split   — carriers in BOTH L and R subtrees
  * L_specific            — carriers only in L  (potential L-defining muts)
  * R_specific            — carriers only in R  (potential R-defining muts)
  * (no_carrier)          — dropped: no carrier at parent

For each side-specific bin we further split by carrier count:
  * *_ge2                 — variants with ≥2 same-side carriers  (real signal)
  * *_single              — variants with exactly 1 carrier      (singleton, likely noise)

We also compare against what SVD's K-means partition assigned to each
side (module_variants[L/R]) so we can measure per-split concordance:

  * n_assigned_L, n_assigned_R                          — from module_variants
  * n_assigned_L_that_are_set_L_specific                — TP for L
  * n_assigned_L_that_are_set_shared                    — clonal contamination
  * n_assigned_L_that_are_set_R_specific                — misassigned to L
  * concordance_L = n_L_that_are_set_L / n_assigned_L   — fraction of L pool that is truly L-side

Carrier definition (configurable): AD ≥ min_ad AND DP ≥ min_dp.  Default 1/1.

Usage
-----
  topdown_shared_private_diag.py --panel-tag all_variants_relax03 \
      [--panel-root results/worm6_final/DNA_analysis/phylogeny/panels] \
      [--min-ad 1 --min-dp 1]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
TD_ROOT = REPO / "results/worm6_final/DNA_analysis/phylogeny/topdown"
PANEL_ROOT_DEFAULT = REPO / "results/worm6_final/DNA_analysis/phylogeny/panels"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-tag", required=True)
    ap.add_argument("--panel-root", type=Path, default=PANEL_ROOT_DEFAULT)
    ap.add_argument("--out-dir", type=Path,
                    default=TD_ROOT / "_summaries")
    ap.add_argument("--min-ad", type=int, default=1,
                    help="Minimum AD to call a cell a carrier for a variant.")
    ap.add_argument("--min-dp", type=int, default=1,
                    help="Minimum DP for a cell/variant to be considered.")
    return ap.parse_args()


def _dense(x):
    if sp.issparse(x): x = x.toarray()
    return np.asarray(x)


def _walk_splits(node, out, path=""):
    if node.get("split") is None:
        return
    out.append((path or "root", node))
    if node.get("left"):
        _walk_splits(node["left"], out, path=(f"{path}.L" if path else "L"))
    if node.get("right"):
        _walk_splits(node["right"], out, path=(f"{path}.R" if path else "R"))


def _classify_split(node: dict, AD_full: np.ndarray, DP_full: np.ndarray,
                    cell_ix: dict, var_ix: dict,
                    min_ad: int, min_dp: int) -> dict:
    """Set-based shared/private counts + concordance vs SVD's K-means assignment."""
    sp_info = node["split"]
    parent_cells = node["cells"]
    L_cells = parent_cells  # placeholder, need to derive from split's module_variants? no — from tree structure
    # tree nodes have their cell partition on the child nodes
    L_cells = node["left"]["cells"]
    R_cells = node["right"]["cells"]

    Lc_ix = np.array([cell_ix[c] for c in L_cells if c in cell_ix], dtype=np.int64)
    Rc_ix = np.array([cell_ix[c] for c in R_cells if c in cell_ix], dtype=np.int64)
    Pc_ix = np.concatenate([Lc_ix, Rc_ix])

    # Carrier bool per cell/variant (bounded by min_ad + min_dp).
    ADp = AD_full[Pc_ix, :]
    DPp = DP_full[Pc_ix, :]
    carrier = (ADp >= min_ad) & (DPp >= min_dp)
    # Split carrier by side
    ADL = AD_full[Lc_ix, :]; DPL = DP_full[Lc_ix, :]
    ADR = AD_full[Rc_ix, :]; DPR = DP_full[Rc_ix, :]
    carrierL = (ADL >= min_ad) & (DPL >= min_dp)
    carrierR = (ADR >= min_ad) & (DPR >= min_dp)
    n_carr_L = carrierL.sum(0)     # per-variant carrier count in L
    n_carr_R = carrierR.sum(0)     # per-variant carrier count in R
    has_L = n_carr_L > 0
    has_R = n_carr_R > 0
    any_carrier = has_L | has_R

    # Set-based classification (over variants with any carrier at parent)
    set_shared      = has_L & has_R
    set_L_specific  = has_L & (~has_R)
    set_R_specific  = has_R & (~has_L)

    L_ge2   = set_L_specific & (n_carr_L >= 2)
    L_singl = set_L_specific & (n_carr_L == 1)
    R_ge2   = set_R_specific & (n_carr_R >= 2)
    R_singl = set_R_specific & (n_carr_R == 1)

    # Concordance vs SVD K-means assignment
    assigned_L = np.array([var_ix[v] for v in sp_info["module_variants"]["L"] if v in var_ix], dtype=np.int64)
    assigned_R = np.array([var_ix[v] for v in sp_info["module_variants"]["R"] if v in var_ix], dtype=np.int64)

    def _bins(idx_arr):
        return dict(
            total=int(idx_arr.size),
            set_L_specific=int(set_L_specific[idx_arr].sum()),
            set_R_specific=int(set_R_specific[idx_arr].sum()),
            set_shared=int(set_shared[idx_arr].sum()),
            set_none=int((~any_carrier[idx_arr]).sum()),
        )

    L_bins = _bins(assigned_L)
    R_bins = _bins(assigned_R)

    n_any = int(any_carrier.sum())
    return dict(
        n_variants_any_carrier=n_any,
        n_shared=int(set_shared.sum()),
        n_L_specific=int(set_L_specific.sum()),
        n_R_specific=int(set_R_specific.sum()),
        n_L_ge2=int(L_ge2.sum()),
        n_L_single=int(L_singl.sum()),
        n_R_ge2=int(R_ge2.sum()),
        n_R_single=int(R_singl.sum()),
        # SVD-assignment concordance
        svd_assigned_L_total=L_bins["total"],
        svd_assigned_L_set_L=L_bins["set_L_specific"],
        svd_assigned_L_set_R=L_bins["set_R_specific"],
        svd_assigned_L_set_shared=L_bins["set_shared"],
        svd_assigned_L_set_none=L_bins["set_none"],
        svd_assigned_R_total=R_bins["total"],
        svd_assigned_R_set_L=R_bins["set_L_specific"],
        svd_assigned_R_set_R=R_bins["set_R_specific"],
        svd_assigned_R_set_shared=R_bins["set_shared"],
        svd_assigned_R_set_none=R_bins["set_none"],
    )


def main():
    args = parse_args()
    panel_root = args.panel_root / args.panel_tag
    if not panel_root.exists():
        raise SystemExit(f"panels not found: {panel_root}")
    td_root = TD_ROOT / args.panel_tag
    args.out_dir.mkdir(parents=True, exist_ok=True)

    trees = sorted(td_root.glob("*/worm_*/tree.json"))
    print(f"[diag] found {len(trees)} trees under {td_root}")

    # Preload per-worm panel h5ads once
    worm_cache: dict[str, tuple[np.ndarray, np.ndarray, dict, dict]] = {}
    def _get_worm(worm: str):
        if worm in worm_cache: return worm_cache[worm]
        p = panel_root / f"worm_{worm}.h5ad"
        a = ad.read_h5ad(p)
        AD = _dense(a.layers["AD"]).astype(np.int32)
        DP = _dense(a.layers["DP"]).astype(np.int32)
        cell_ix = {c: i for i, c in enumerate(a.obs_names.astype(str))}
        var_ix  = {v: i for i, v in enumerate(a.var_names.astype(str))}
        worm_cache[worm] = (AD, DP, cell_ix, var_ix)
        return worm_cache[worm]

    rows = []
    for tj in trees:
        config = tj.parent.parent.name
        worm = tj.parent.name.replace("worm_", "")
        try:
            AD, DP, cell_ix, var_ix = _get_worm(worm)
        except FileNotFoundError:
            print(f"  [skip] no panel h5ad for {worm}")
            continue
        d = json.loads(tj.read_text())
        splits = []
        _walk_splits(d["root"], splits)
        for path, node in splits:
            cls = _classify_split(node, AD, DP, cell_ix, var_ix,
                                  args.min_ad, args.min_dp)
            cls.update(dict(
                config=config, worm=worm, path=path,
                depth=node["depth"], n_cells=node["n_cells"],
                n_cells_L=node["left"]["n_cells"],
                n_cells_R=node["right"]["n_cells"],
            ))
            rows.append(cls)

    df = pd.DataFrame(rows)
    col_order = [
        "config","worm","path","depth","n_cells","n_cells_L","n_cells_R",
        "n_variants_any_carrier",
        "n_shared","n_L_specific","n_R_specific",
        "n_L_ge2","n_L_single","n_R_ge2","n_R_single",
        "svd_assigned_L_total","svd_assigned_L_set_L","svd_assigned_L_set_shared",
        "svd_assigned_L_set_R","svd_assigned_L_set_none",
        "svd_assigned_R_total","svd_assigned_R_set_R","svd_assigned_R_set_shared",
        "svd_assigned_R_set_L","svd_assigned_R_set_none",
    ]
    df = df[col_order]
    out = args.out_dir / f"{args.panel_tag}__shared_private.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[diag] wrote {out}  ({len(df)} splits)")

    # Console summary
    print("\n=== set-based mutation counts by depth (median across all splits, all configs) ===")
    print(df.groupby("depth").agg(
        n_splits=("depth","count"),
        med_any=("n_variants_any_carrier","median"),
        med_shared=("n_shared","median"),
        med_L_spec=("n_L_specific","median"),
        med_R_spec=("n_R_specific","median"),
        med_L_ge2=("n_L_ge2","median"),
        med_R_ge2=("n_R_ge2","median"),
        med_L_singleton=("n_L_single","median"),
        med_R_singleton=("n_R_single","median"),
    ).round(1).to_string())

    print("\n=== SVD-assignment concordance by depth (mean fraction, all configs) ===")
    dd = df.copy()
    dd["concord_L"] = dd["svd_assigned_L_set_L"]  / dd["svd_assigned_L_total"].clip(1)
    dd["concord_R"] = dd["svd_assigned_R_set_R"]  / dd["svd_assigned_R_total"].clip(1)
    dd["contam_L"]  = dd["svd_assigned_L_set_R"]  / dd["svd_assigned_L_total"].clip(1)
    dd["contam_R"]  = dd["svd_assigned_R_set_L"]  / dd["svd_assigned_R_total"].clip(1)
    dd["shared_L"]  = dd["svd_assigned_L_set_shared"] / dd["svd_assigned_L_total"].clip(1)
    dd["shared_R"]  = dd["svd_assigned_R_set_shared"] / dd["svd_assigned_R_total"].clip(1)
    print(dd.groupby("depth").agg(
        n_splits=("depth","count"),
        concord_L=("concord_L","mean"),
        concord_R=("concord_R","mean"),
        contam_L =("contam_L","mean"),
        contam_R =("contam_R","mean"),
        shared_L =("shared_L","mean"),
        shared_R =("shared_R","mean"),
    ).round(3).to_string())

    print("\n=== per-config: median L_specific-ge2 + R_specific-ge2 across all splits ===")
    dd["true_side_ge2"] = dd["n_L_ge2"] + dd["n_R_ge2"]
    print(dd.groupby("config").agg(
        n_splits=("depth","count"),
        med_true_side_ge2=("true_side_ge2","median"),
        med_shared=("n_shared","median"),
        med_singleton=("n_L_single","median"),
    ).round(1).to_string())


if __name__ == "__main__":
    main()
