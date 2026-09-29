#!/usr/bin/env python3
"""Concatenate per-worm panel h5ads into one combined AnnData.

Layers preserved: AD, DP (imputed) and AD_raw, DP_raw (pre-imputation).
Variant matrix = UNION across worms. Per-cell entries for variants not in
that worm's panel are filled with NaN (float32), so downstream code can
distinguish "not called in this worm" from "0 reads observed".

Usage
-----
  combine_panel_h5ads.py --panel-tag <tag> --out-path <path.h5ad>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
PANEL_ROOT = REPO / "results/worm6_final/DNA_analysis/phylogeny/panels"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel-tag", required=True)
    p.add_argument("--out-path", type=Path, required=True)
    p.add_argument("--layers", nargs="+", default=["AD","DP","AD_raw","DP_raw"])
    return p.parse_args()


def _dense_layer(a: ad.AnnData, name: str) -> np.ndarray:
    x = a.layers[name]
    if sp.issparse(x): x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def main():
    args = parse_args()
    pdir = PANEL_ROOT / args.panel_tag
    if not pdir.exists():
        raise SystemExit(f"panel dir not found: {pdir}")

    files = sorted(pdir.glob("worm_*.h5ad"))
    if not files:
        raise SystemExit(f"no worm h5ads in {pdir}")
    print(f"[combine] {len(files)} worm files under {pdir}")

    # First pass: build var union (preserve first-seen var metadata)
    var_order: list[str] = []
    var_seen: set[str] = set()
    var_meta_frames: list[pd.DataFrame] = []

    per_worm: dict[str, ad.AnnData] = {}
    for f in files:
        worm = f.stem.replace("worm_", "")
        a = ad.read_h5ad(f)
        per_worm[worm] = a
        new_vars = [v for v in a.var_names if v not in var_seen]
        if new_vars:
            var_seen.update(new_vars)
            var_order.extend(new_vars)
            df = a.var.loc[new_vars].copy()
            var_meta_frames.append(df)
        print(f"  {worm}: cells={a.n_obs}  vars={a.n_vars}  (union so far: {len(var_order)})")

    var_meta = pd.concat(var_meta_frames, axis=0)
    # Keep only columns that appear in the first frame (avoid schema drift)
    ref_cols = var_meta_frames[0].columns
    var_meta = var_meta[[c for c in ref_cols if c in var_meta.columns]]
    var_meta = var_meta.loc[var_order]

    N_vars = len(var_order)
    var_ix = {v: i for i, v in enumerate(var_order)}

    # Build obs frame + layer matrices
    obs_frames = []
    layer_arrays = {name: [] for name in args.layers}

    for worm, a in per_worm.items():
        obs = a.obs.copy()
        obs["worm"] = worm
        obs_frames.append(obs)
        n = a.n_obs
        col_ix = np.array([var_ix[v] for v in a.var_names], dtype=np.int64)
        for name in args.layers:
            if name not in a.layers:
                # Fill NaN for missing layer
                arr = np.full((n, N_vars), np.nan, dtype=np.float32)
            else:
                block = _dense_layer(a, name)
                arr = np.full((n, N_vars), np.nan, dtype=np.float32)
                arr[:, col_ix] = block
            layer_arrays[name].append(arr)

    obs = pd.concat(obs_frames, axis=0)
    layers = {name: np.vstack(chunks) for name, chunks in layer_arrays.items()}
    # X = imputed DP so a bare view is informative; caller can layer-swap.
    X = layers.get("DP", layers[args.layers[0]])

    combined = ad.AnnData(
        X=X,
        obs=obs,
        var=var_meta,
        layers=layers,
        uns={
            "panel_tag": args.panel_tag,
            "n_worms": len(per_worm),
            "combined_from": [str(f) for f in files],
            "variant_axis": "union",
            "notes": "AD/DP are imputed panel values; AD_raw/DP_raw are pre-imputation. NaN = variant not in that worm's panel.",
        },
    )
    print(f"[combine] combined shape: {combined.n_obs} cells x {combined.n_vars} vars, "
          f"layers={list(combined.layers.keys())}")

    args.out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(args.out_path, compression="gzip")
    size_mb = args.out_path.stat().st_size / 1e6
    print(f"[combine] wrote {args.out_path} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
