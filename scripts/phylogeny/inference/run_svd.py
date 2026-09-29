#!/usr/bin/env python3
"""Run production SVDBipartitionBuilder on a per-worm panel.

Writes tree.json in the phylo.tl.trees format.

Usage
-----
  run_svd.py --panel-h5ad phylogeny/panels/<tag>/worm_<W>.h5ad \
             --out-dir  phylogeny/inference/<tag>/svd/worm_<W> \
             [--ncomp 3]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
sys.path.insert(0, str(REPO / "scripts" / "phylogeny"))

from phylo.tl.trees import SVDBipartitionBuilder, build_worm_tree  # noqa: E402


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-h5ad", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--ncomp", type=int, default=3)
    ap.add_argument("--min-variants-to-split", type=int, default=20)
    ap.add_argument("--min-clade-size", type=int, default=10)
    ap.add_argument("--max-depth", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--core-selection", default="none",
                    choices=("none", "centroid", "mquad", "setop", "hybrid"),
                    help="Per-node variant filter before final K-means refit.")
    ap.add_argument("--n-core", type=int, default=None,
                    help="For centroid core selection: keep top-N variants closest to K-means centroid.")
    ap.add_argument("--mquad-delta-bic-threshold", type=float, default=5.0,
                    help="For mquad core selection: minimum delta-BIC to retain a variant.")
    ap.add_argument("--mquad-nproc", type=int, default=1)
    ap.add_argument("--hybrid-switch-cells", type=int, default=20,
                    help="Hybrid mode: use mquad when clade size >= this, else set-op filter.")
    ap.add_argument("--setop-min-carriers", type=int, default=2,
                    help="Set-op filter: keep variants with >= this many carriers on one side and 0 on the other.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(args.panel_h5ad)
    worm = a.uns.get("panel_worm", args.panel_h5ad.stem.replace("worm_", ""))
    print(f"[svd] worm={worm}  n_cells={a.n_obs}  n_vars={a.n_vars}  ncomp={args.ncomp}")

    builder = SVDBipartitionBuilder(
        ncomp=args.ncomp,
        min_variants_to_split=args.min_variants_to_split,
        seed=args.seed,
        core_selection=args.core_selection,
        n_core=args.n_core,
        mquad_delta_bic_threshold=args.mquad_delta_bic_threshold,
        mquad_nproc=args.mquad_nproc,
        hybrid_switch_cells=args.hybrid_switch_cells,
        setop_min_carriers=args.setop_min_carriers,
    )
    tree = build_worm_tree(
        a, builder=builder,
        min_clade_size=args.min_clade_size,
        max_depth=args.max_depth,
        worm=str(worm),
        params={
            "method": "svd_bipart",
            "ncomp": args.ncomp,
            "core_selection": args.core_selection,
            "n_core": args.n_core,
            "mquad_delta_bic_threshold": args.mquad_delta_bic_threshold,
            "hybrid_switch_cells": args.hybrid_switch_cells,
            "setop_min_carriers": args.setop_min_carriers,
        },
    )
    tree.to_json(args.out_dir / "tree.json")
    print(f"[svd] wrote {args.out_dir}/tree.json")


if __name__ == "__main__":
    main()
