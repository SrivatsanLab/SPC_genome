#!/usr/bin/env python3
"""Summarize top-down SVD sweep trees.

For each (panel, config, worm) tree.json, walks every internal node and
records per-clade properties:

  * depth
  * n_cells (parent), n_cells_L, n_cells_R    -- clade sizes at the split
  * n_variants_used, n_core                    -- carrier variants pre/post core selection
  * n_muts_L, n_muts_R                         -- module_variants counts (mutations "assigned" per side)
  * vaf_L, vaf_R                               -- within-clade pooled VAF at that side's module vars
  * purity_L, purity_R                         -- median purity contrast per side

Emits:
  * one long-format TSV: <panel>__nodes.tsv
  * one config x worm summary TSV: <panel>__summary.tsv

Usage
-----
  summarize_topdown_sweep.py --panel-tag all_variants_relax03
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
TD_ROOT = REPO / "results/worm6_final/DNA_analysis/phylogeny/topdown"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-tag", required=True)
    ap.add_argument("--out-dir", type=Path,
                    default=REPO / "results/worm6_final/DNA_analysis/phylogeny/topdown/_summaries")
    return ap.parse_args()


def _walk(node, config, worm, rows, path=""):
    if node.get("split") is None:
        return
    sp = node["split"]
    ncL = sp["module_n_cells"]["L"]
    ncR = sp["module_n_cells"]["R"]
    rows.append(dict(
        config=config, worm=worm, path=path or "root", depth=node["depth"],
        n_cells=node["n_cells"],
        n_cells_L=ncL, n_cells_R=ncR,
        min_child_size=min(ncL, ncR),
        n_variants_used=sp["n_variants_used"],
        n_core=sp["n_core"],
        n_muts_L=len(sp["module_variants"]["L"]),
        n_muts_R=len(sp["module_variants"]["R"]),
        vaf_L=sp["module_pooled_vaf"]["L"],
        vaf_R=sp["module_pooled_vaf"]["R"],
        purity_L=sp["module_purity_median"]["L"],
        purity_R=sp["module_purity_median"]["R"],
    ))
    if node.get("left"):
        _walk(node["left"], config, worm, rows,
              path=(f"{path}.L" if path else "L"))
    if node.get("right"):
        _walk(node["right"], config, worm, rows,
              path=(f"{path}.R" if path else "R"))


def main():
    args = parse_args()
    panel_root = TD_ROOT / args.panel_tag
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    per_tree = []
    trees = sorted(panel_root.glob("*/worm_*/tree.json"))
    print(f"[summarize] {len(trees)} trees under {panel_root}")
    for tj in trees:
        try:
            d = json.loads(tj.read_text())
        except Exception as e:
            print(f"  bad json: {tj}: {e}")
            continue
        config = tj.parent.parent.name
        worm = tj.parent.name.replace("worm_", "")
        _walk(d["root"], config, worm, rows)
        # per-tree summary
        internal = [r for r in rows if r["config"] == config and r["worm"] == worm]
        n_leaves = sum(1 for _ in _leaves(d["root"]))
        max_depth_reached = max((r["depth"] for r in internal), default=0)
        singleton_leaves = sum(1 for L in _leaves(d["root"]) if L["n_cells"] == 1)
        two_cell_leaves = sum(1 for L in _leaves(d["root"]) if L["n_cells"] == 2)
        per_tree.append(dict(
            config=config, worm=worm,
            n_cells_total=d["root"]["n_cells"],
            n_splits=len(internal),
            n_leaves=n_leaves,
            max_depth=max_depth_reached,
            n_singleton_leaves=singleton_leaves,
            n_two_cell_leaves=two_cell_leaves,
            fraction_fully_resolved=(singleton_leaves + two_cell_leaves) / max(n_leaves, 1),
        ))

    long_df = pd.DataFrame(rows)
    summary_df = pd.DataFrame(per_tree)

    out_long = args.out_dir / f"{args.panel_tag}__nodes.tsv"
    out_sum  = args.out_dir / f"{args.panel_tag}__summary.tsv"
    long_df.to_csv(out_long, sep="\t", index=False)
    summary_df.to_csv(out_sum, sep="\t", index=False)
    print(f"[summarize] wrote {out_long}  ({len(long_df)} splits)")
    print(f"[summarize] wrote {out_sum}  ({len(summary_df)} trees)")

    # Console summary: per-config aggregates
    print("\n=== per-config aggregate over all worms ===")
    agg = summary_df.groupby("config").agg(
        n_worms=("worm", "count"),
        mean_splits=("n_splits", "mean"),
        mean_leaves=("n_leaves", "mean"),
        mean_max_depth=("max_depth", "mean"),
        mean_frac_resolved=("fraction_fully_resolved", "mean"),
    ).round(2)
    print(agg.to_string())

    print("\n=== per-config aggregate over all splits ===")
    long_agg = long_df.groupby("config").agg(
        n_splits=("depth", "count"),
        mean_depth=("depth", "mean"),
        med_muts_L=("n_muts_L", "median"),
        med_muts_R=("n_muts_R", "median"),
        med_vaf_L=("vaf_L", "median"),
        med_vaf_R=("vaf_R", "median"),
        med_purity_L=("purity_L", "median"),
        med_purity_R=("purity_R", "median"),
        med_n_vars_used=("n_variants_used", "median"),
    ).round(3)
    print(long_agg.to_string())


def _leaves(node):
    if node.get("left") is None and node.get("right") is None:
        yield node
        return
    if node.get("left") is not None:
        yield from _leaves(node["left"])
    if node.get("right") is not None:
        yield from _leaves(node["right"])


if __name__ == "__main__":
    main()
