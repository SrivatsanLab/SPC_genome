#!/usr/bin/env python3
"""Encode a phangorn-ready MP matrix (naive 0/1/?) from a per-worm panel.

Writes at <out_dir>/<tag>.tsv (default tag = "real"):
  cell_id \t char_string
where char_string is a per-variant string in {'0', '1', '?'}:
  AD >= 1                       -> '1'
  DP >= <dp_threshold>, AD == 0 -> '0'   (default dp_threshold = 1)
  otherwise                     -> '?'

Under Camin-Sokal parsimony (Sankoff with cost 0→1 = 1, 1→0 = ∞), the
observed '0' states are the anchor that forces the tree to explain '1's
as gain events. A presence-only encoding (chars in {1, ?}) would make Fitch
score identically 0 (any tree explains a single-observed-state character
with 0 changes when '?' is ambiguous), producing a star polytomy — hence
the naive encoding here.

The --dp-threshold flag raises the confidence bar for calling '0'. At
DP=1 a covered-ref observation is one Bernoulli draw at ~50% under true
heterozygosity, so dropout is common and '0' calls are noisy. Raising to
DP≥2 or ≥3 pushes low-confidence covered-refs to '?', at the cost of
losing informative characters.

Also writes cells.txt (row order) and variants.txt (column order).

Usage
-----
  encode_mp_matrix.py --panel-h5ad phylogeny/panels/<tag>/worm_<W>.h5ad \
                      --out-dir phylogeny/inference/<tag>/mp/worm_<W> \
                      [--dp-threshold 1] [--tag real]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.sparse import issparse


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-h5ad", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--dp-threshold", type=int, default=1,
                    help="Minimum raw DP required to call '0' (covered-ref). "
                         "Positions with 0 < DP < threshold and AD == 0 become '?'. "
                         "Default 1 (backwards compatible).")
    ap.add_argument("--tag", type=str, default="real",
                    help="Output filename tag (default 'real'). "
                         "Writes <out_dir>/<tag>.tsv.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(args.panel_h5ad)
    AD = a.layers["AD"].toarray() if issparse(a.layers["AD"]) else np.asarray(a.layers["AD"])
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else np.asarray(a.layers["DP"])
    cells = list(a.obs_names.astype(str))
    variants = list(a.var_names.astype(str))
    n_c, n_v = AD.shape
    K = args.dp_threshold
    print(f"[mp] worm={a.uns.get('panel_worm', '?')}  {n_c} cells × {n_v} variants  "
          f"dp_threshold=K={K}  tag={args.tag}")

    out_tsv = args.out_dir / f"{args.tag}.tsv"
    with open(out_tsv, "w") as f:
        f.write("cell_id\tchars\n")
        for i, c in enumerate(cells):
            # '1' if any alt read; '0' only if covered at >= K without alt; else '?'
            chars = np.where(AD[i] >= 1, "1",
                             np.where(DP[i] >= K, "0", "?"))
            f.write(f"{c}\t{''.join(chars)}\n")
    with open(args.out_dir / "cells.txt", "w") as f:
        f.write("\n".join(cells))
    with open(args.out_dir / "variants.txt", "w") as f:
        f.write("\n".join(variants))
    # Quick sanity summary
    total = n_c * n_v
    n_ad = int((AD >= 1).sum())
    n_ref = int(((AD == 0) & (DP >= K)).sum())
    n_q   = int(((AD == 0) & (DP < K)).sum())
    print(f"[mp] chars: '1'={n_ad:,} ({100*n_ad/total:.1f}%)  "
          f"'0'={n_ref:,} ({100*n_ref/total:.1f}%)  "
          f"'?'={n_q:,} ({100*n_q/total:.1f}%)")
    print(f"[mp] wrote {out_tsv.name}, cells.txt, variants.txt")


if __name__ == "__main__":
    main()
