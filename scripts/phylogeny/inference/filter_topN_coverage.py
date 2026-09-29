#!/usr/bin/env python3
"""Rank variants by within-worm cell coverage and keep the top N.

For each per-worm panel in <panel_dir>/, compute for every variant
  n_covered = (# cells with DP > 0)
rank descending, keep the top N. If a worm has fewer than N variants, keep
all of them. Writes derived panels to <panel_dir>__topN<N>/. All var/uns/obs
metadata is preserved; ``uns['panel_config']`` gains ``top_N``.

Usage
-----
  filter_topN_coverage.py \
      --panels-dir results/worm6_final/DNA_analysis/phylogeny/panels/<tag> \
      --out-dir    results/worm6_final/DNA_analysis/phylogeny/panels/<tag>__topN<N> \
      --top-n 500
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--top-n", type=int, required=True)
    return ap.parse_args()


def filter_panel(h5ad_path: Path, top_n: int, out_path: Path) -> tuple[int, int]:
    a = ad.read_h5ad(h5ad_path)
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else a.layers["DP"]
    n_cov = (DP > 0).sum(axis=0)
    order = np.argsort(-n_cov, kind="stable")
    keep_idx = order[:min(top_n, a.n_vars)]
    keep_idx = np.sort(keep_idx)  # preserve original var order
    sub = a[:, keep_idx].copy()
    sub.uns["panel_config"] = {**dict(a.uns.get("panel_config", {})), "top_N": int(top_n)}
    sub.write_h5ad(out_path)
    return a.n_vars, sub.n_vars


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for p in sorted(args.panels_dir.glob("worm_worm*.h5ad")):
        n_before, n_after = filter_panel(p, args.top_n, args.out_dir / p.name)
        summary_rows.append({"worm": p.stem.replace("worm_", ""),
                             "n_before": n_before, "n_after": n_after})
        print(f"  {p.name}: {n_before} → {n_after}")

    # Carry over provenance + config
    for aux in ("summary.csv", "cells.csv", "blacklist.tsv",
                "per_worm_bulk.h5ad", "panel_config.resolved.yaml"):
        src = args.panels_dir / aux
        if src.exists():
            (args.out_dir / aux).write_bytes(src.read_bytes())

    pd.DataFrame(summary_rows).to_csv(args.out_dir / "topN_summary.tsv", sep="\t", index=False)
    print(f"[filter] wrote {len(summary_rows)} panels + topN_summary.tsv")


if __name__ == "__main__":
    main()
