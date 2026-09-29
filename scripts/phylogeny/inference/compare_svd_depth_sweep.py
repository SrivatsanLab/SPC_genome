#!/usr/bin/env python3
"""Compare ML inference across the deep-SVD imputation sweep.

For each of two source panels (tct_kept_relax03 and all_variants_relax03),
we have four imputed variants:
  - baseline:   svdImp                 (min_clade=10, BIC=5,  original)
  - variant A:  svdImp_mquad5_cs5     (min_clade=5,  BIC=5,  deeper)
  - variant B:  svdImp_mquad10_cs5    (min_clade=5,  BIC=10, deeper + stricter markers)
  - variant C:  svdImp_mquad10_cs3    (min_clade=3,  BIC=10, deepest + stricter markers)

For each variant we compare against the unimputed panel and against each
other on:
  1. SVD-tree depth achieved (from tree.json).
  2. Number of pseudo-alt reads injected during imputation (uns['imputation']).
  3. Per-worm bootstrap support stats on the CellPhy real tree.
  4. Robinson–Foulds distance between paired trees (variant vs baseline, and
     variant vs unimputed).

Output is a series of tables printed to stdout. No files written.

Usage
-----
  micromamba run -n cellphy python scripts/phylogeny/inference/compare_svd_depth_sweep.py
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import ete3
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import wilcoxon

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
INFER = REPO / "results/worm6_final/DNA_analysis/phylogeny/inference"
PANELS = REPO / "results/worm6_final/DNA_analysis/phylogeny/panels"

PANEL_ROOTS = ("tct_kept_relax03", "all_variants_relax03")
# label -> (svd_subdir, imputed_tag_suffix)
VARIANTS = [
    ("baseline",        "svd_mquad5",        "__svdImp"),
    ("mquad5_cs5",      "svd_mquad5_cs5",    "__svdImp_mquad5_cs5"),
    ("mquad10_cs5",     "svd_mquad10_cs5",   "__svdImp_mquad10_cs5"),
    ("mquad10_cs3",     "svd_mquad10_cs3",   "__svdImp_mquad10_cs3"),
    ("Sym_mquad10_cs3", "svd_mquad10_cs3",   "__svdImpSym_mquad10_cs3"),
]
WORMS = "worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20".split()


# ---------------------------------------------------------------------------
# Tree loading + metrics
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


def rf_norm(a: ete3.Tree, b: ete3.Tree) -> float | float:
    """Unrooted normalized Robinson-Foulds. Returns NaN if incomparable."""
    try:
        rf, max_rf, common, *_ = a.robinson_foulds(b, unrooted_trees=True)
        return rf / max_rf if max_rf else np.nan
    except Exception:
        return np.nan


def svd_tree_stats(tree_json_path: Path) -> dict:
    """Return summary of an SVD tree.json: n_splits, max_depth, per-depth n_splits."""
    if not tree_json_path.exists():
        return dict(n_splits=0, max_depth=0, leaves=0)
    tree = json.loads(tree_json_path.read_text())
    depths = []
    n_leaves = 0

    def walk(n):
        nonlocal n_leaves
        if "split" in n and n["split"] is not None:
            depths.append(int(n.get("depth", 0)))
            walk(n["left"]); walk(n["right"])
        else:
            n_leaves += 1

    walk(tree["root"])
    return dict(n_splits=len(depths),
                max_depth=max(depths) if depths else 0,
                leaves=n_leaves)


def imputation_stats(panel_h5_path: Path) -> dict:
    """Return the uns['imputation'] dict from an imputed panel."""
    if not panel_h5_path.exists():
        return dict()
    a = ad.read_h5ad(panel_h5_path, backed="r")
    d = dict(a.uns.get("imputation", {}))
    d["n_cells"] = a.n_obs
    d["n_vars"] = a.n_vars
    return d


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_ml_tree_for(panel_root: str, variant_suffix: str, worm: str) -> ete3.Tree | None:
    p = INFER / f"{panel_root}{variant_suffix}" / "cellphy" / f"worm_{worm}" / "real" / "sup.raxml.supportFBP"
    if not p.exists():
        return None
    return load_tree(p)


def load_unimp_tree(panel_root: str, worm: str) -> ete3.Tree | None:
    p = INFER / panel_root / "cellphy" / f"worm_{worm}" / "real" / "sup.raxml.supportFBP"
    if not p.exists():
        return None
    return load_tree(p)


def main() -> None:
    pd.set_option("display.width", 240)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.float_format", "{:.2f}".format)

    for panel_root in PANEL_ROOTS:
        print("\n" + "=" * 100)
        print(f"PANEL: {panel_root}")
        print("=" * 100)

        # ---- (1) SVD tree depths achieved --------------------------------
        print(f"\n--- (1) SVD tree depths achieved (n_splits, max_depth per worm) ---")
        rows = []
        for variant, svd_subdir, _ in VARIANTS:
            for w in WORMS:
                p = INFER / panel_root / svd_subdir / f"worm_{w}" / "tree.json"
                s = svd_tree_stats(p)
                rows.append(dict(variant=variant, worm=w, **s))
        svd_df = pd.DataFrame(rows)
        # Compact wide view: mean n_splits and max_depth per variant.
        agg = svd_df.groupby("variant").agg(
            mean_n_splits=("n_splits", "mean"),
            median_n_splits=("n_splits", "median"),
            max_n_splits=("n_splits", "max"),
            mean_max_depth=("max_depth", "mean"),
            median_max_depth=("max_depth", "median"),
            max_max_depth=("max_depth", "max"),
        ).round(2)
        # Preserve variant order
        agg = agg.reindex([v[0] for v in VARIANTS])
        print(agg.to_string())

        # ---- (2) Imputation payload -------------------------------------
        print(f"\n--- (2) Pseudo-alt entries injected per worm ---")
        rows = []
        for variant, _, imp_suffix in VARIANTS:
            for w in WORMS:
                p = PANELS / f"{panel_root}{imp_suffix}" / f"worm_{w}.h5ad"
                s = imputation_stats(p)
                rows.append(dict(variant=variant, worm=w,
                                 n_imputed=s.get("n_imputed_entries", np.nan),
                                 n_cells=s.get("n_cells", np.nan),
                                 n_vars=s.get("n_vars", np.nan)))
        imp_df = pd.DataFrame(rows)
        imp_agg = imp_df.groupby("variant").agg(
            mean_n_imputed=("n_imputed", "mean"),
            median_n_imputed=("n_imputed", "median"),
            total_n_imputed=("n_imputed", "sum"),
        ).round(0)
        imp_agg = imp_agg.reindex([v[0] for v in VARIANTS])
        print(imp_agg.to_string())

        # ---- (3) Bootstrap support on CellPhy real trees -----------------
        print(f"\n--- (3) Bootstrap support on CellPhy real trees ---")
        support_rows = []
        for variant, _, imp_suffix in VARIANTS:
            for w in WORMS:
                t = load_ml_tree_for(panel_root, imp_suffix, w)
                if t is None:
                    continue
                s = support_values(t)
                support_rows.append(dict(
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
        sup_df = pd.DataFrame(support_rows)
        sup_agg = sup_df.groupby("variant").agg(
            n_worms=("worm", "count"),
            mean_support=("mean_sup", "mean"),
            median_support=("median_sup", "mean"),
            pct_ge50=("pct_ge50", "mean"),
            pct_ge70=("pct_ge70", "mean"),
            pct_ge90=("pct_ge90", "mean"),
            mean_tbl=("tbl", "mean"),
            mean_height=("ht", "mean"),
        )
        sup_agg = sup_agg.reindex([v[0] for v in VARIANTS])
        print(sup_agg.round(3).to_string())

        # Per-worm delta vs baseline on mean support
        print(f"\n--- (3b) Per-worm mean-support delta vs baseline ---")
        wide = sup_df.pivot(index="worm", columns="variant", values="mean_sup")
        wide = wide[[v[0] for v in VARIANTS if v[0] in wide.columns]]
        if "baseline" in wide.columns:
            for v in wide.columns:
                if v == "baseline":
                    continue
                wide[f"d_{v}"] = wide[v] - wide["baseline"]
            print(wide.round(2).to_string())

            # Wilcoxon paired vs baseline
            print(f"\n--- (3c) Wilcoxon paired (variant > baseline on per-worm mean support) ---")
            for v in [x[0] for x in VARIANTS if x[0] != "baseline"]:
                a = wide[v].dropna()
                b = wide["baseline"].reindex(a.index).dropna()
                idx = a.index.intersection(b.index)
                if len(idx) < 3:
                    print(f"  {v}: n<3, skipped")
                    continue
                try:
                    st = wilcoxon(a.loc[idx], b.loc[idx], alternative="greater")
                    n_wins = int((a.loc[idx] > b.loc[idx]).sum())
                    print(f"  {v}: wins {n_wins}/{len(idx)}  W={st.statistic:.1f}  p={st.pvalue:.3g}")
                except ValueError as e:
                    print(f"  {v}: {e}")

        # ---- (4) Robinson-Foulds distances ------------------------------
        print(f"\n--- (4) RF distance vs baseline (normalized 0-1; 0=identical, 1=maximally different) ---")
        rf_rows = []
        for w in WORMS:
            base_t = load_ml_tree_for(panel_root, "__svdImp", w)
            unimp_t = load_unimp_tree(panel_root, w)
            for variant, _, imp_suffix in VARIANTS:
                if variant == "baseline":
                    continue
                var_t = load_ml_tree_for(panel_root, imp_suffix, w)
                if var_t is None or base_t is None:
                    continue
                rf_rows.append(dict(
                    worm=w,
                    variant=variant,
                    rf_vs_baseline=rf_norm(var_t, base_t),
                    rf_vs_unimputed=rf_norm(var_t, unimp_t) if unimp_t is not None else np.nan,
                ))
        rf_df = pd.DataFrame(rf_rows)
        if rf_df.empty:
            print("  (no variants to compare)")
        else:
            rf_agg = rf_df.groupby("variant").agg(
                mean_rf_vs_baseline=("rf_vs_baseline", "mean"),
                median_rf_vs_baseline=("rf_vs_baseline", "median"),
                mean_rf_vs_unimputed=("rf_vs_unimputed", "mean"),
                median_rf_vs_unimputed=("rf_vs_unimputed", "median"),
            ).round(3)
            rf_agg = rf_agg.reindex([v[0] for v in VARIANTS if v[0] != "baseline"])
            print(rf_agg.to_string())

        # ---- (5) baseline vs unimputed RF for context -------------------
        print(f"\n--- (5) Baseline (svdImp) RF vs unimputed for reference ---")
        rows = []
        for w in WORMS:
            base_t = load_ml_tree_for(panel_root, "__svdImp", w)
            unimp_t = load_unimp_tree(panel_root, w)
            if base_t is None or unimp_t is None:
                continue
            rows.append(dict(worm=w, rf_norm=rf_norm(base_t, unimp_t)))
        base_unimp_df = pd.DataFrame(rows)
        if not base_unimp_df.empty:
            print(f"  mean={base_unimp_df.rf_norm.mean():.3f}  "
                  f"median={base_unimp_df.rf_norm.median():.3f}  "
                  f"range=[{base_unimp_df.rf_norm.min():.3f}, {base_unimp_df.rf_norm.max():.3f}]")


if __name__ == "__main__":
    main()
