#!/usr/bin/env python3
"""Compare distance-based tree inference across the imputation sweep.

Mirrors compare_svd_depth_sweep.py but on the distance/*.support outputs.

For each source panel × distance pipeline × imputation variant:
  1. Per-worm bootstrap support stats (mean, median, %>=50/70/90).
  2. Per-worm tree geometry (total branch length, tree height).
  3. RF distance between paired trees:
     - each variant vs. unimputed baseline (does imputation change topology?)
     - each variant vs. its own CellPhy tree (cross-method concordance)
  4. Wilcoxon paired test: variant > unimputed on per-worm mean support.

No files written; tables printed to stdout.

Usage
-----
  micromamba run -n cellphy python \
      scripts/phylogeny/inference/compare_distance_imputation.py
"""
from __future__ import annotations

from pathlib import Path

import ete3
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
INFER = REPO / "results/worm6_final/DNA_analysis/phylogeny/inference"

PANEL_ROOTS = ("tct_kept_relax03", "all_variants_relax03")
# label -> suffix (empty for unimputed baseline)
VARIANTS = [
    ("unimputed",     ""),
    ("svdImp",        "__svdImp"),
    ("mquad5_cs5",    "__svdImp_mquad5_cs5"),
    ("mquad10_cs5",   "__svdImp_mquad10_cs5"),
    ("mquad10_cs3",   "__svdImp_mquad10_cs3"),
    ("Sym_mquad10_cs3", "__svdImpSym_mquad10_cs3"),
]
PIPELINES = ("soft__nj_nni__mad", "soft__nj_nni__og", "soft__bme_nni__og")
WORMS = "worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20".split()


# ---------------------------------------------------------------------------

def load_tree(path: Path) -> ete3.Tree:
    t = ete3.Tree(str(path), format=0)
    for lf in t.get_leaves():
        lf.name = lf.name.strip("'")
    return t


def support_values(tree) -> np.ndarray:
    out = []
    for n in tree.traverse():
        if n.is_leaf() or n.is_root():
            continue
        s = getattr(n, "support", None)
        try:
            f = float(s)
            if 0.0 <= f <= 100.0 and not np.isnan(f):
                out.append(f)
        except (TypeError, ValueError):
            pass
    return np.array(out)


def total_bl(tree) -> float:
    return float(sum(float(n.dist or 0.0) for n in tree.traverse() if not n.is_root()))


def tree_height(tree) -> float:
    m = 0.0
    for lf in tree.get_leaves():
        d = 0.0
        n = lf
        while n.up is not None:
            d += float(n.dist or 0.0)
            n = n.up
        m = max(m, d)
    return m


def rf_norm(a: ete3.Tree, b: ete3.Tree) -> float:
    try:
        rf, max_rf, *_ = a.robinson_foulds(b, unrooted_trees=True)
        return rf / max_rf if max_rf else np.nan
    except Exception:
        return np.nan


def dist_tree_path(panel_root: str, variant_suffix: str, pipeline: str, worm: str) -> Path:
    return INFER / f"{panel_root}{variant_suffix}" / "distance" / f"worm_{worm}" / f"{pipeline}.support"


def cellphy_tree_path(panel_root: str, variant_suffix: str, worm: str) -> Path:
    return INFER / f"{panel_root}{variant_suffix}" / "cellphy" / f"worm_{worm}" / "real" / "sup.raxml.supportFBP"


# ---------------------------------------------------------------------------

def main():
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.float_format", "{:.2f}".format)

    for panel_root in PANEL_ROOTS:
        print("\n" + "=" * 130)
        print(f"PANEL: {panel_root}")
        print("=" * 130)

        for pipeline in PIPELINES:
            print(f"\n--- Pipeline: {pipeline} ---")

            # ------------------------------------------------------------------
            # (1) Per-worm bootstrap support
            # ------------------------------------------------------------------
            rows = []
            for variant, suffix in VARIANTS:
                for w in WORMS:
                    p = dist_tree_path(panel_root, suffix, pipeline, w)
                    if not p.exists():
                        continue
                    t = load_tree(p)
                    s = support_values(t)
                    rows.append(dict(
                        variant=variant, worm=w,
                        n_leaves=len(t),
                        n_int=len(s),
                        mean_sup=float(s.mean()) if len(s) else np.nan,
                        median_sup=float(np.median(s)) if len(s) else np.nan,
                        pct_ge50=float((s >= 50).mean() * 100) if len(s) else np.nan,
                        pct_ge70=float((s >= 70).mean() * 100) if len(s) else np.nan,
                        pct_ge90=float((s >= 90).mean() * 100) if len(s) else np.nan,
                        tbl=total_bl(t),
                        ht=tree_height(t),
                    ))
            df = pd.DataFrame(rows)
            if df.empty:
                print("  (no trees for this pipeline)")
                continue

            agg = df.groupby("variant").agg(
                n_worms=("worm", "count"),
                mean_support=("mean_sup", "mean"),
                median_support=("median_sup", "mean"),
                pct_ge50=("pct_ge50", "mean"),
                pct_ge70=("pct_ge70", "mean"),
                pct_ge90=("pct_ge90", "mean"),
                mean_tbl=("tbl", "mean"),
                mean_height=("ht", "mean"),
            )
            order = [v[0] for v in VARIANTS if v[0] in agg.index]
            agg = agg.reindex(order)
            print("Bootstrap support & geometry (mean across worms):")
            print(agg.round(3).to_string())

            # (1b) Per-worm delta vs unimputed
            wide = df.pivot(index="worm", columns="variant", values="mean_sup")
            wide = wide[[v for v in order if v in wide.columns]]
            if "unimputed" in wide.columns:
                print(f"\nPer-worm mean-support delta vs unimputed:")
                deltas = pd.DataFrame(index=wide.index)
                for v in [x for x in order if x != "unimputed"]:
                    if v in wide.columns:
                        deltas[f"d_{v}"] = wide[v] - wide["unimputed"]
                cols = ["unimputed"] + [c for c in wide.columns if c != "unimputed"]
                combined = wide[cols].join(deltas)
                print(combined.round(2).to_string())

                # (1c) Wilcoxon paired vs unimputed
                print(f"\nWilcoxon paired (variant > unimputed on per-worm mean support):")
                for v in [x for x in order if x != "unimputed"]:
                    if v not in wide.columns:
                        continue
                    a = wide[v].dropna()
                    b = wide["unimputed"].reindex(a.index).dropna()
                    idx = a.index.intersection(b.index)
                    if len(idx) < 3:
                        print(f"  {v}: n<3 skipped")
                        continue
                    try:
                        st = wilcoxon(a.loc[idx], b.loc[idx], alternative="greater")
                        n_wins = int((a.loc[idx] > b.loc[idx]).sum())
                        print(f"  {v:<20s}: wins {n_wins}/{len(idx)}  W={st.statistic:.1f}  p={st.pvalue:.3g}")
                    except ValueError as e:
                        print(f"  {v}: {e}")

            # ------------------------------------------------------------------
            # (2) Robinson-Foulds distance vs unimputed and vs baseline svdImp
            # ------------------------------------------------------------------
            print(f"\nRF distance (normalized 0-1, unrooted) vs unimputed and svdImp:")
            rf_rows = []
            for w in WORMS:
                unimp_t = None
                base_t = None
                p_un = dist_tree_path(panel_root, "", pipeline, w)
                p_sv = dist_tree_path(panel_root, "__svdImp", pipeline, w)
                if p_un.exists(): unimp_t = load_tree(p_un)
                if p_sv.exists(): base_t = load_tree(p_sv)
                for variant, suffix in VARIANTS:
                    if variant in ("unimputed", "svdImp"):
                        continue
                    p = dist_tree_path(panel_root, suffix, pipeline, w)
                    if not p.exists():
                        continue
                    var_t = load_tree(p)
                    rf_rows.append(dict(
                        worm=w, variant=variant,
                        rf_vs_unimp=rf_norm(var_t, unimp_t) if unimp_t is not None else np.nan,
                        rf_vs_svdImp=rf_norm(var_t, base_t) if base_t is not None else np.nan,
                    ))
            rf_df = pd.DataFrame(rf_rows)
            if not rf_df.empty:
                rf_agg = rf_df.groupby("variant").agg(
                    mean_rf_vs_unimp=("rf_vs_unimp", "mean"),
                    median_rf_vs_unimp=("rf_vs_unimp", "median"),
                    mean_rf_vs_svdImp=("rf_vs_svdImp", "mean"),
                    median_rf_vs_svdImp=("rf_vs_svdImp", "median"),
                ).round(3)
                order2 = [v[0] for v in VARIANTS if v[0] not in ("unimputed","svdImp") and v[0] in rf_agg.index]
                rf_agg = rf_agg.reindex(order2)
                print(rf_agg.to_string())

            # (2b) Reference: baseline svdImp vs unimputed on this pipeline
            base_rows = []
            for w in WORMS:
                p_un = dist_tree_path(panel_root, "", pipeline, w)
                p_sv = dist_tree_path(panel_root, "__svdImp", pipeline, w)
                if p_un.exists() and p_sv.exists():
                    base_rows.append(rf_norm(load_tree(p_un), load_tree(p_sv)))
            if base_rows:
                arr = np.array(base_rows)
                print(f"\nReference (baseline svdImp vs unimputed):"
                      f" mean={arr.mean():.3f} median={np.median(arr):.3f} "
                      f"range=[{arr.min():.3f}, {arr.max():.3f}]")

            # ------------------------------------------------------------------
            # (3) Cross-method: distance tree vs its own CellPhy tree
            # ------------------------------------------------------------------
            print(f"\nCross-method RF (distance {pipeline} vs CellPhy real):")
            cm_rows = []
            for variant, suffix in VARIANTS:
                if variant == "unimputed":
                    # No CellPhy tree on the unimputed panel via this workflow;
                    # skip if it wasn't run.
                    continue
                for w in WORMS:
                    p_d = dist_tree_path(panel_root, suffix, pipeline, w)
                    p_c = cellphy_tree_path(panel_root, suffix, w)
                    if not (p_d.exists() and p_c.exists()):
                        continue
                    d_t = load_tree(p_d); c_t = load_tree(p_c)
                    cm_rows.append(dict(variant=variant, worm=w,
                                        rf_dist_vs_cellphy=rf_norm(d_t, c_t)))
            cm_df = pd.DataFrame(cm_rows)
            if not cm_df.empty:
                cm_agg = cm_df.groupby("variant").agg(
                    n_worms=("worm", "count"),
                    mean_rf_dist_vs_cellphy=("rf_dist_vs_cellphy", "mean"),
                    median_rf_dist_vs_cellphy=("rf_dist_vs_cellphy", "median"),
                ).round(3)
                order3 = [v[0] for v in VARIANTS if v[0] in cm_agg.index]
                cm_agg = cm_agg.reindex(order3)
                print(cm_agg.to_string())


if __name__ == "__main__":
    main()
