#!/usr/bin/env python3
"""Render the MP encoding + parsimony-model sweep trees for one panel.

Layout of the sweep on disk:
    scratch/mp_sweep/<panel>__<worm>__<config>/real.majority.newick

For each panel, emits under <out-dir>:
  - by_worm/<worm>__<panel>.png    : configs-in-a-row grid (method comparison)
  - by_config/<config>__<panel>.png: worms-in-a-grid (per-method worm sweep)
  - overview__<panel>.png          : worms × configs matrix

Usage
-----
  plot_mp_sweep.py \
      --panel-tag <panel> \
      --sweep-root scratch/mp_sweep \
      --out-dir   results/worm6_final/DNA_analysis/phylogeny/figures/mp_sweep
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

from phylo.pl.cellphy_trees import plot_cellphy_tree  # noqa: E402


WORMS = [
    "worm01","worm03","worm05","worm06","worm07","worm08","worm09","worm10",
    "worm11","worm12","worm13","worm15","worm16","worm17","worm18","worm20",
]
# Preferred config order for grid columns (only rendered if present on disk).
CONFIG_ORDER = [
    "baseline_K1_kInf",
    "relax_K1_k20", "relax_K1_k10", "relax_K1_k5", "relax_K1_k3",
    "fitch_K1_k1",
    "stringent_K3_kInf", "combined_K3_k5",
    "dollo_g10_l1", "dollo_g100_l1", "dollo_g1000_l1",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-tag", required=True)
    ap.add_argument("--sweep-root", type=Path, default=REPO / "scratch/mp_sweep")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--support-threshold", type=float, default=50.0)
    ap.add_argument("--cmap", default="viridis")
    return ap.parse_args()


def _tree_path(sweep_root: Path, panel: str, worm: str, config: str) -> Path | None:
    d = sweep_root / f"{panel}__{worm}__{config}"
    for name in ("real.majority.newick", "real.mp_support.newick", "real.mp.newick"):
        p = d / name
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def _discover(sweep_root: Path, panel: str) -> tuple[list[str], list[str]]:
    """Return (present_worms, present_configs) for this panel, ordered."""
    import re
    worm_re = re.compile(r"^worm\d{2}$")
    dirs = list(sweep_root.glob(f"{panel}__*__*"))
    worms, configs = set(), set()
    prefix = f"{panel}__"
    for d in dirs:
        rest = d.name[len(prefix):]
        # rest is <worm>__<config>. worms are wormNN.
        parts = rest.split("__", 1)
        if len(parts) != 2:
            continue
        w, c = parts
        # Guard against panel prefix collisions (e.g. "tct_kept_relax03" is a
        # prefix of "tct_kept_relax03__svdImp"). Only accept wormNN tokens.
        if not worm_re.match(w):
            continue
        worms.add(w); configs.add(c)
    present_worms = [w for w in WORMS if w in worms]
    present_configs = [c for c in CONFIG_ORDER if c in configs]
    # append any unknown configs at end (won't happen normally)
    for c in sorted(configs):
        if c not in present_configs:
            present_configs.append(c)
    return present_worms, present_configs


def _plot_grid(sweep_root: Path, panel: str, worms: list[str], configs: list[str],
               out_path: Path, support_threshold: float, cmap: str,
               title: str, per_cell_w: float = 3.2, per_cell_h: float = 2.6):
    n_rows, n_cols = len(worms), len(configs)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(per_cell_w * n_cols, per_cell_h * n_rows),
                             squeeze=False)
    for i, w in enumerate(worms):
        for j, c in enumerate(configs):
            ax = axes[i, j]
            tp = _tree_path(sweep_root, panel, w, c)
            if tp is None:
                ax.text(0.5, 0.5, "missing", ha="center", va="center",
                        transform=ax.transAxes, color="0.6")
                ax.set_xticks([]); ax.set_yticks([])
                if i == 0: ax.set_title(c, fontsize=9)
                if j == 0: ax.set_ylabel(w, fontsize=9)
                continue
            try:
                plot_cellphy_tree(tp, out_path=None, ax=ax,
                                  title=None,
                                  cmap_name=cmap,
                                  support_threshold=support_threshold,
                                  show_colorbar=False)
            except Exception as e:
                ax.text(0.5, 0.5, f"err\n{type(e).__name__}",
                        ha="center", va="center",
                        transform=ax.transAxes, color="red", fontsize=6)
                ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(c, fontsize=9)
            if j == 0:
                ax.set_ylabel(w, fontsize=9, rotation=0, ha="right", va="center", labelpad=15)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    for suf in (".png",):
        fig.savefig(out_path.with_suffix(suf), dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"[grid] {out_path.with_suffix('.png')}")


def _plot_per_worm(sweep_root: Path, panel: str, worm: str, configs: list[str],
                   out_dir: Path, support_threshold: float, cmap: str):
    """One row of trees per config, for a single worm."""
    n = len(configs)
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.2 * nrows),
                             squeeze=False)
    axes_flat = axes.ravel()
    for i, c in enumerate(configs):
        ax = axes_flat[i]
        tp = _tree_path(sweep_root, panel, worm, c)
        if tp is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes, color="0.6")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(c, fontsize=10)
            continue
        try:
            plot_cellphy_tree(tp, out_path=None, ax=ax,
                              title=c, cmap_name=cmap,
                              support_threshold=support_threshold,
                              show_colorbar=False)
        except Exception as e:
            ax.text(0.5, 0.5, f"err\n{type(e).__name__}",
                    ha="center", va="center",
                    transform=ax.transAxes, color="red")
            ax.set_xticks([]); ax.set_yticks([])
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle(f"{worm} · panel={panel}  — MP-sweep configs", fontsize=12)
    fig.tight_layout()
    p = out_dir / f"{worm}__{panel}.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"[per_worm] {p}")


def _plot_per_config(sweep_root: Path, panel: str, config: str, worms: list[str],
                     out_dir: Path, support_threshold: float, cmap: str):
    """A 4×4 (or similar) grid of worms for a single (panel, config)."""
    n = len(worms)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.2 * nrows),
                             squeeze=False)
    axes_flat = axes.ravel()
    for i, w in enumerate(worms):
        ax = axes_flat[i]
        tp = _tree_path(sweep_root, panel, w, config)
        if tp is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes, color="0.6")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(w, fontsize=10)
            continue
        try:
            plot_cellphy_tree(tp, out_path=None, ax=ax,
                              title=w, cmap_name=cmap,
                              support_threshold=support_threshold,
                              show_colorbar=False)
        except Exception as e:
            ax.text(0.5, 0.5, f"err\n{type(e).__name__}",
                    ha="center", va="center",
                    transform=ax.transAxes, color="red")
            ax.set_xticks([]); ax.set_yticks([])
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle(f"panel={panel} · config={config}", fontsize=12)
    fig.tight_layout()
    p = out_dir / f"{config}__{panel}.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"[per_config] {p}")


def main() -> None:
    args = parse_args()
    panel = args.panel_tag
    sweep_root = Path(args.sweep_root)
    out_root = Path(args.out_dir) / panel
    (out_root / "by_worm").mkdir(parents=True, exist_ok=True)
    (out_root / "by_config").mkdir(parents=True, exist_ok=True)

    worms, configs = _discover(sweep_root, panel)
    if not worms or not configs:
        print(f"[plot_mp_sweep] nothing to plot for {panel}")
        return
    print(f"[plot_mp_sweep] panel={panel}  worms={len(worms)}  configs={len(configs)}")

    for w in worms:
        _plot_per_worm(sweep_root, panel, w, configs, out_root / "by_worm",
                       args.support_threshold, args.cmap)
    for c in configs:
        _plot_per_config(sweep_root, panel, c, worms, out_root / "by_config",
                         args.support_threshold, args.cmap)

    _plot_grid(sweep_root, panel, worms, configs,
               out_root / f"overview__{panel}",
               args.support_threshold, args.cmap,
               title=f"MP-sweep overview · panel={panel}")


if __name__ == "__main__":
    main()
