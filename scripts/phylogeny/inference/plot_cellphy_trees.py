#!/usr/bin/env python3
"""Render per-worm CellPhy trees with nodes colored by FBP bootstrap support.

Usage
-----
  plot_cellphy_trees.py \
      --inference-dir results/worm6_final/DNA_analysis/phylogeny/inference/<tag> \
      --matrix real \
      --out-dir       results/worm6_final/DNA_analysis/phylogeny/inference/<tag>/figures
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
sys.path.insert(0, str(REPO / "scripts" / "phylogeny"))

from phylo.pl.cellphy_trees import plot_all_cellphy_trees  # noqa: E402


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inference-dir", type=Path, required=True,
                    help="phylogeny/inference/<panel_tag>")
    ap.add_argument("--matrix", default="real",
                    choices=("real", "covmask", "shuffled"))
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--cmap", default="viridis")
    ap.add_argument("--support-threshold", type=float, default=50.0,
                    help="Edges below this FBP fade to light gray (default 50).")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    plot_all_cellphy_trees(
        cellphy_dir=args.inference_dir / "cellphy",
        out_dir=args.out_dir,
        matrix=args.matrix,
        cmap_name=args.cmap,
        support_threshold=args.support_threshold,
    )


if __name__ == "__main__":
    main()
