#!/usr/bin/env python
"""Sweep panel VAF band × tree ncomp; compare per-variant module VAF distributions.

For each combination:
- write a config that overrides `panels.somatic_vaf_band` and `trees.ncomp`
- run `phylo panels --force` and `phylo trees --force`
- iterate the resulting trees; for each split's L/R module, look up per-variant
  own_bulk_vaf from the worm's panel h5ad and record distribution stats.

Outputs go under
``results/worm6_final/DNA_analysis/phylogeny/sweeps/panel_ncomp/``:

- ``summary.csv``           one row per (config, worm, node_path, side)
- ``config_summary.csv``    per-config aggregate: median module VAF at root
- ``vaf_distributions.png`` box-plot grid, one panel per config

The metric of interest: does *any* config produce split modules whose
per-variant own_bulk_vaf concentrates around the plan §6.1 target
(0.20–0.45 for a real somatic het branch at α≈0.3)? If not, the SVD
builder is picking up structural noise regardless of panel choice.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
PHYLO_ROOT = HERE.parent
OUTPUT_ROOT = REPO / "results/worm6_final/DNA_analysis/phylogeny"
SWEEP_DIR = OUTPUT_ROOT / "sweeps/panel_ncomp"
CONFIGS_DIR = SWEEP_DIR / "configs"

# Sweep grid --------------------------------------------------------------
BANDS = [
    ("wide_lo05", [None, 0.05]),
    ("wide_lo10", [None, 0.10]),
    ("mid",       [0.02, 0.30]),
    ("wide_lo30", [None, 0.30]),
]
NCOMPS = [2, 4, 8]


def _write_config(band_name: str, band: list, ncomp: int) -> tuple[Path, str]:
    """Write a config overriding just band + ncomp; return (path, run_name)."""
    run_name = f"sw_band-{band_name}_ncomp-{ncomp:02d}"
    cfg = {
        "name": run_name,
        "inputs": {
            "adata": str(REPO / "results/worm6_final/DNA_analysis/joint_variants_process.h5ad"),
        },
        "output_root": str(OUTPUT_ROOT),
        "panels": {"somatic_vaf_band": band},
        "trees": {"ncomp": ncomp},
    }
    CONFIGS_DIR.mkdir(parents=True, exist_ok=True)
    p = CONFIGS_DIR / f"{run_name}.yaml"
    with open(p, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    return p, run_name


def _run_stage(stage: str, cfg_path: Path) -> None:
    subprocess.run(
        [sys.executable, "-m", "phylo", stage, "--config", str(cfg_path), "--force"],
        check=True,
        capture_output=True,
        cwd=str(PHYLO_ROOT),
    )


def _collect_run(run_name: str, band_name: str, band: list, ncomp: int) -> list[dict]:
    """For every split in every worm tree, extract per-variant own_bulk_vaf stats."""
    run_dir = OUTPUT_ROOT / "runs" / run_name
    panels_dir = run_dir / "panels"
    trees_dir = run_dir / "trees"
    rows: list[dict] = []
    for tree_json in sorted(trees_dir.glob("worm_*/tree.json")):
        worm = tree_json.parent.name.replace("worm_", "")
        panel_h5 = panels_dir / f"worm_{worm}.h5ad"
        if not panel_h5.exists():
            continue
        panel = ad.read_h5ad(panel_h5)
        vaf_lookup = dict(zip(panel.var_names.astype(str), panel.var["own_bulk_vaf"].to_numpy()))

        with open(tree_json) as f:
            tree = json.load(f)

        for node in _walk(tree["root"]):
            if "split" not in node:
                continue
            for side in ("L", "R"):
                mvars = node["split"]["module_variants"][side]
                vafs = np.array([vaf_lookup.get(v, np.nan) for v in mvars])
                vafs = vafs[~np.isnan(vafs)]
                if vafs.size == 0:
                    continue
                rows.append(
                    {
                        "band_name": band_name,
                        "band_lo": band[0],
                        "band_hi": band[1],
                        "ncomp": ncomp,
                        "worm": worm,
                        "node_path": node["path"] or "root",
                        "depth": node["depth"],
                        "side": side,
                        "n_variants": int(vafs.size),
                        "median_vaf": float(np.median(vafs)),
                        "p25_vaf": float(np.percentile(vafs, 25)),
                        "p75_vaf": float(np.percentile(vafs, 75)),
                        "frac_ge_0.10": float((vafs >= 0.10).mean()),
                        "frac_ge_0.20": float((vafs >= 0.20).mean()),
                        "frac_ge_0.35": float((vafs >= 0.35).mean()),
                        "pooled_vaf": float(node["split"]["module_pooled_vaf"][side]),
                    }
                )
    return rows


def _walk(node: dict):
    yield node
    if "left" in node:
        yield from _walk(node["left"])
    if "right" in node:
        yield from _walk(node["right"])


def _plot_distributions(df: pd.DataFrame, out_path: Path) -> None:
    """Box-plot of per-variant own_bulk_vaf per module, faceted by (band, ncomp)."""
    bands = df["band_name"].unique()
    ncomps = sorted(df["ncomp"].unique())
    fig, axes = plt.subplots(len(ncomps), len(bands), figsize=(4.5 * len(bands), 3.2 * len(ncomps)),
                             sharey=True, sharex=True)
    axes = np.atleast_2d(axes)
    ylim_hi = float(df[["p75_vaf", "median_vaf"]].max().max()) * 1.15
    for i, nc in enumerate(ncomps):
        for j, band in enumerate(bands):
            ax = axes[i, j]
            sub = df[(df["ncomp"] == nc) & (df["band_name"] == band)]
            if sub.empty:
                ax.axis("off")
                continue
            # Group by depth; box the median_vaf across modules
            depths = sorted(sub["depth"].unique())
            data = [sub[sub["depth"] == d]["median_vaf"].to_numpy() for d in depths]
            positions = list(range(len(depths)))
            ax.boxplot(data, positions=positions, widths=0.6, showfliers=True)
            ax.set_xticks(positions)
            ax.set_xticklabels([f"d={d}" for d in depths], fontsize=8)
            ax.axhline(0.20, color="green", ls="--", lw=0.8, alpha=0.6, label="lower plausible")
            ax.axhline(0.35, color="crimson", ls="--", lw=0.8, alpha=0.6, label="§6.1 target")
            ax.set_ylim(-0.02, ylim_hi)
            if i == 0:
                ax.set_title(f"band={band}", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"ncomp={nc}\nmed(module VAF)", fontsize=9)
            ax.tick_params(axis="y", labelsize=8)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
    fig.suptitle("Per-module median of per-variant own_bulk_vaf, by depth × config", y=1.005, fontsize=11)
    axes[0, -1].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    for band_name, band in BANDS:
        for ncomp in NCOMPS:
            cfg_path, run_name = _write_config(band_name, band, ncomp)
            print(f"[sweep] === {run_name} ===")
            for stage in ("panels", "trees"):
                _run_stage(stage, cfg_path)
            rows = _collect_run(run_name, band_name, band, ncomp)
            print(f"[sweep]   {len(rows)} modules recorded")
            all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    df.to_csv(SWEEP_DIR / "summary.csv", index=False)

    # Aggregate per config
    agg = (
        df.groupby(["band_name", "ncomp"])
        .agg(
            n_modules=("median_vaf", "size"),
            median_module_vaf=("median_vaf", "median"),
            p75_module_vaf=("median_vaf", lambda s: np.percentile(s, 75)),
            frac_modules_ge_0p20=("frac_ge_0.20", "mean"),
            median_pooled_vaf=("pooled_vaf", "median"),
            root_median_vaf=("median_vaf", lambda s: np.nan),
        )
        .reset_index()
    )
    # root_median_vaf: median per-variant VAF for root (depth 0) modules only
    root = df[df["depth"] == 0].groupby(["band_name", "ncomp"])["median_vaf"].median().rename("root_median_vaf")
    agg = agg.drop(columns="root_median_vaf").merge(root, on=["band_name", "ncomp"], how="left")
    agg = agg.sort_values("frac_modules_ge_0p20", ascending=False)
    agg.to_csv(SWEEP_DIR / "config_summary.csv", index=False)

    _plot_distributions(df, SWEEP_DIR / "vaf_distributions.png")

    print("\n[sweep] config ranking (best = highest fraction of modules with median VAF ≥ 0.20):")
    print(agg.to_string(index=False))
    print(f"\n[sweep] outputs → {SWEEP_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
