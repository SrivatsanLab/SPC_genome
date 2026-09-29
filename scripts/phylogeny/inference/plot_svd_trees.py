#!/usr/bin/env python3
"""Render SVD bipartition trees with per-split VAF + purity annotations.

Wraps ``phylo.pl.trees.plot_tree``: with a panel h5ad, each split shows
vaf=X / n_mod=Y / pur=±Z. Purity is the L-vs-R VAF contrast at that
side's module variants — scale-invariant to panel size.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
sys.path.insert(0, str(REPO / "scripts" / "phylogeny"))

from phylo.pl.trees import _clade_confidences, load_tree, plot_tree  # noqa: E402


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inference-dir", type=Path, required=True)
    ap.add_argument("--svd-subdir", default="svd", help="name of SVD subdir under inference-dir (svd | svd_centroid500 | svd_mquad5)")
    ap.add_argument("--panels-dir", type=Path, required=True,
                    help="Panel h5ads for the same tag (per-worm h5ads used for purity computation).")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--cmap", default="viridis")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    svd_root = args.inference_dir / args.svd_subdir
    trees = sorted(svd_root.glob("worm_*/tree.json"))
    if not trees:
        print(f"no tree.json under {svd_root}"); return

    run_vmax = 0.0
    for tj in trees:
        _, root, _ = load_tree(tj)
        c = _clade_confidences(root)
        if c: run_vmax = max(run_vmax, max(c.values()))
    if run_vmax <= 0.0: run_vmax = 1.0

    def _panel_for(worm_dirname):
        # worm_dirname like "worm_worm07" → panel at panels_dir/worm_worm07.h5ad
        return args.panels_dir / f"{worm_dirname}.h5ad"

    rendered = []
    for tj in trees:
        worm = tj.parent.name
        panel_h5 = _panel_for(worm)
        panel_arg = panel_h5 if panel_h5.exists() else None
        for ext in (".png", ".svg"):
            plot_tree(tj, args.out_dir / f"{worm}__{args.svd_subdir}{ext}",
                      cmap_name=args.cmap, vmax=run_vmax,
                      panel_h5ad=panel_arg)
        rendered.append((worm, tj, panel_arg))
        print(f"[svd_plot] {worm}__{args.svd_subdir}.{{png,svg}}  panel={'yes' if panel_arg else 'no'}")

    # Overview grid
    n = len(rendered)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 4.0 * nrows))
    axes_flat = axes.ravel() if hasattr(axes, "ravel") else [axes]
    for i, (worm, tj, panel_arg) in enumerate(rendered):
        plot_tree(tj, out_path=None, ax=axes_flat[i], cmap_name=args.cmap,
                  vmax=run_vmax, show_colorbar=False, title=worm,
                  panel_h5ad=panel_arg)
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle(f"SVD trees ({args.svd_subdir}) — {args.inference_dir.name}", fontsize=14)
    fig.tight_layout()
    for suf in (".png", ".svg"):
        fig.savefig(args.out_dir / f"overview__{args.svd_subdir}{suf}", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[svd_plot] overview__svd.{{png,svg}}")


if __name__ == "__main__":
    main()
