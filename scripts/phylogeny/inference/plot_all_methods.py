#!/usr/bin/env python3
"""Render per-worm trees for cellphy / distance-NJ / distance-UPGMA / MP.

For a given panel tag, writes into
  inference/<tag>/figures/trees/
per-worm PNG+SVG for every available tree source, plus overview grids.

Usage
-----
  plot_all_methods.py \
      --inference-dir results/worm6_final/DNA_analysis/phylogeny/inference/<tag> \
      --out-dir       results/worm6_final/DNA_analysis/phylogeny/inference/<tag>/figures/trees \
      [--methods cellphy distance mp]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
sys.path.insert(0, str(REPO / "scripts" / "phylogeny"))

from phylo.pl.cellphy_trees import (  # noqa: E402
    plot_all_cellphy_trees,
    plot_all_distance_trees,
    plot_all_mp_trees,
)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inference-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--methods", nargs="+",
                    default=["cellphy", "distance", "mp"],
                    choices=["cellphy", "distance", "mp"])
    ap.add_argument("--cellphy-matrices", nargs="+",
                    default=["real", "covmask", "shuffled"])
    ap.add_argument("--distance-metric", default="soft")
    ap.add_argument("--support-threshold", type=float, default=50.0)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if "cellphy" in args.methods:
        for m in args.cellphy_matrices:
            plot_all_cellphy_trees(
                cellphy_dir=args.inference_dir / "cellphy",
                out_dir=args.out_dir,
                matrix=m,
                support_threshold=args.support_threshold,
            )
    if "distance" in args.methods:
        plot_all_distance_trees(
            inference_dir=args.inference_dir,
            out_dir=args.out_dir,
            metric=args.distance_metric,
            support_threshold=args.support_threshold,
        )
    if "mp" in args.methods:
        plot_all_mp_trees(
            inference_dir=args.inference_dir,
            out_dir=args.out_dir,
            support_threshold=args.support_threshold,
        )


if __name__ == "__main__":
    main()
