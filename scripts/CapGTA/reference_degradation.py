#!/usr/bin/env python3
"""Reference degradation: fit CapGTA detectability surface, thin CENGEN.

Reproduces the winning-config degraded CENGEN reference used by
`notebooks/worm6_final_GEX.ipynb`. Pipeline (single script):

    1. Load CENGEN L4 (Taylor 2020), downsample neurons to `--neuron-frac`
       (default 0.32; matches whole-worm somatic composition).
    2. Pseudobulk both sides per anchor tissue on gene proportions.
    3. GAM fit of log(cap_prop / ref_prop) ~ s(log_exon_bp) + s(n_junctions)
       + β·polyA_ind + β·SL_ind (pooled across tissues).
    4. Predict p_g per reference gene, rescale by `total_preserve` (pick c
       such that E[binom(X, c·p_g)] = X_total, then clip to [0,1]).
    5. Binomial-thin the (downsampled) reference N times → N h5ads.

Anchor CapGTA clusters use the arm B (leiden 0.6) mapping in DEFAULT_ANCHOR_CAPGTA;
override with `--anchor-json` if calling on a different annotation.

Usage:
    reference_degradation.py \\
        --capgta   <capgta_annotated.h5ad> \\
        --reference <cengen_taylor2020.h5ad> \\
        --calibration <calibration.h5ad> \\
        --kg-source <capgta_scaled_none.h5ad> \\
        --out <out_dir> \\
        [--neuron-frac 0.32] [--n-seeds 5] [--seed 0]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


DEFAULT_ANCHOR_CAPGTA = {
    "3":  "excretory_gland",
    "8":  "excretory_gland",
    "10": "excretory_gland",
    "5":  "peptidergic_neurons",
    "4":  "glia_sheath",
    "6":  "coelomocyte",
}

DEFAULT_ANCHOR_REF = {
    "excretory_gland":     ["Excretory_gland_cell"],
    "peptidergic_neurons": ["_ALL_NEURONS_"],
    "glia_sheath":         ["AMsh", "PHsh", "CEPsh"],
    "coelomocyte":         ["Coelomocyte"],
}

NON_NEURON_TAGS = {
    "Epidermis", "Coelomocyte", "Spermatheca", "Germline", "Body_wall_muscle",
    "Body_wall_muscle_anterior", "Intestine", "Excretory_cell",
    "Excretory_gland_cell", "Rectal_gland", "Anal_muscle", "AMsh", "PHsh",
    "CEPsh", "AMso", "PHso", "Pharyngeal_muscle", "Pharyngeal_gland_cell",
    "Vulval_muscle", "Vulval_cells", "Uterine_cell",
    "Spermathecal-uterine_junction_or_uterine_toroid", "Gonadal_sheath_cell",
    "Distal_tip_cell", "Sperm", "Seam_cell", "Arcade_cell", "Marginal_cell",
    "hmc", "XXX", "Unannotated", "Unknown_non_neuronal",
    "Unknown_non_neuronal_2", "Unknown_non-neuronal_3",
    "Glia_1", "Glia_2", "Glia_3", "Glia_4", "Glia_5",
}


def _downsample_neurons(A: ad.AnnData, target_frac: float, seed: int,
                        label_col: str = "cell_type") -> ad.AnnData:
    ct = A.obs[label_col].astype(str)
    is_neuron = ~ct.isin(NON_NEURON_TAGS)
    n_neuron = int(is_neuron.sum())
    n_non_neuron = int((~is_neuron).sum())
    target = int(round(target_frac * n_non_neuron / (1 - target_frac)))
    if target >= n_neuron:
        print(f"  target {target} >= available {n_neuron}; keeping all")
        return A
    rng = np.random.default_rng(seed)
    neuron_idx = np.where(is_neuron.to_numpy())[0]
    keep_neuron = rng.choice(neuron_idx, size=target, replace=False)
    keep = np.concatenate([np.where(~is_neuron.to_numpy())[0], keep_neuron])
    keep.sort()
    out = A[keep].copy()
    print(f"  neuron-downsample: {n_neuron} → {target} neurons "
          f"(non-neuron kept {n_non_neuron}); output {out.n_obs} cells")
    return out


def _is_short_neuron_code(name: str) -> bool:
    return (2 <= len(name) <= 5) and name.replace("_", "").isalnum() and (
        name.isupper() or (len(name) >= 3 and name[0].isupper()
                           and any(c.isdigit() for c in name))
    )


def _resolve_ref_anchor_types(anchor_ref: dict, ref_types: list[str]) -> dict:
    out: dict[str, list[str]] = {}
    non_neuron = {t for spec in anchor_ref.values()
                  for t in spec if t != "_ALL_NEURONS_"} | NON_NEURON_TAGS
    for label, spec in anchor_ref.items():
        if any(t == "_ALL_NEURONS_" for t in spec):
            out[label] = [t for t in ref_types
                          if _is_short_neuron_code(t) and t not in non_neuron]
        else:
            out[label] = spec
    return out


def _pseudobulk_proportions(A: ad.AnnData, tissue_col: str) -> pd.DataFrame:
    X = A.X
    tissues = A.obs[tissue_col].astype(str).values
    unique = sorted(t for t in np.unique(tissues) if t and t != "nan")
    props = {}
    for t in unique:
        mask = tissues == t
        if mask.sum() == 0:
            continue
        totals = np.asarray(X[mask].sum(axis=0)).ravel()
        s = totals.sum()
        if s > 0:
            props[t] = totals / s
    return pd.DataFrame(props, index=A.var_names)


def _build_gene_features(ref_var_names: pd.Index, kg_source: Path,
                         calibration: Path) -> pd.DataFrame:
    K = ad.read_h5ad(kg_source, backed="r")
    wb = K.var["wormbase_id"].astype(str).values
    df_k = pd.DataFrame({
        "wb": wb,
        "n_junctions": K.var["n_junctions"].astype(float).values,
        "exon_bp": K.var["exonic_length"].astype(float).values,
    }).groupby("wb").max()

    C = ad.read_h5ad(calibration)
    polyat = C.layers["polyat"]; sl = C.layers["sl"]
    polyA_cells = (np.asarray((polyat > 0).sum(axis=0)).ravel()
                   if sp.issparse(polyat) else (polyat > 0).sum(axis=0))
    sl_cells = (np.asarray((sl > 0).sum(axis=0)).ravel()
                if sp.issparse(sl) else (sl > 0).sum(axis=0))
    df_c = pd.DataFrame({"polyA_cells": polyA_cells, "sl_cells": sl_cells},
                        index=C.var_names)

    feats = df_k.reindex(ref_var_names).fillna(
        {"n_junctions": 0, "exon_bp": 0}
    )
    feats = feats.join(df_c, how="left").fillna(
        {"polyA_cells": 0, "sl_cells": 0}
    )
    feats["log_exon_bp"] = np.log10(feats["exon_bp"].clip(lower=1))
    feats["polyA_ind"] = (feats["polyA_cells"] >= 5).astype(int)
    feats["sl_ind"] = (feats["sl_cells"] >= 5).astype(int)
    return feats


def _fit_gam(fit_frame: pd.DataFrame, n_splines: int = 8):
    from pygam import LinearGAM, s, l
    X = fit_frame[["log_exon_bp", "n_junctions", "polyA_ind", "sl_ind"]].values
    y = fit_frame["log_ratio"].values
    gam = LinearGAM(s(0, n_splines=n_splines) + s(1, n_splines=n_splines)
                    + l(2) + l(3))
    gam.fit(X, y)
    return gam


def _predict_p_g(gam, feats: pd.DataFrame) -> np.ndarray:
    X = feats[["log_exon_bp", "n_junctions", "polyA_ind", "sl_ind"]].values
    return np.exp(gam.predict(X))


def _rescale_total_preserve(p_raw: np.ndarray, X_ref,
                            max_iter: int = 20, tol: float = 1e-3) -> np.ndarray:
    """Pick c iteratively so E[sum(binom(X, c·p_g))] = sum(X), clipping to 1."""
    gene_totals = (np.asarray(X_ref.sum(axis=0)).ravel() if sp.issparse(X_ref)
                   else X_ref.sum(axis=0))
    target_total = gene_totals.sum()
    denom = (gene_totals * p_raw).sum()
    c = target_total / max(denom, 1e-12)
    p = np.clip(c * p_raw, 0.0, 1.0)
    for _ in range(max_iter):
        est = (gene_totals * p).sum()
        if abs(est - target_total) / max(target_total, 1e-12) < tol:
            break
        c *= target_total / max(est, 1e-12)
        p = np.clip(c * p_raw, 0.0, 1.0)
    return p


def _binomial_thin(X, p_g: np.ndarray, rng: np.random.Generator):
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    Xc = X.tocoo(copy=True)
    n = Xc.data.astype(np.int64)
    p = np.clip(p_g[Xc.col], 0.0, 1.0)
    thinned = rng.binomial(n, p)
    keep = thinned > 0
    return sp.csr_matrix(
        (thinned[keep].astype(np.int64), (Xc.row[keep], Xc.col[keep])),
        shape=X.shape,
    )


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--capgta", required=True, type=Path,
                   help="CapGTA arm B annotated h5ad "
                        "(has leiden_0.6 and wormbase_id in .var).")
    p.add_argument("--reference", required=True, type=Path,
                   help="CENGEN L4 (Taylor 2020) full h5ad, WBGene-indexed.")
    p.add_argument("--calibration", required=True, type=Path,
                   help="Calibration h5ad with polyat + sl layers, WBGene-indexed.")
    p.add_argument("--kg-source", required=True, type=Path,
                   help="Preprocessed CapGTA h5ad with n_junctions + exonic_length "
                        "in .var + wormbase_id.")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--capgta-cluster-key", default="leiden_0.6")
    p.add_argument("--capgta-counts-layer", default="counts")
    p.add_argument("--ref-label-key", default="cell_type")
    p.add_argument("--anchor-json", type=Path, default=None,
                   help="Override CapGTA cluster→anchor mapping "
                        "(JSON dict {cluster: anchor_label}).")
    p.add_argument("--neuron-frac", type=float, default=0.32,
                   help="Target neuron fraction in downsampled reference "
                        "(0 to skip downsampling).")
    p.add_argument("--n-seeds", type=int, default=5)
    p.add_argument("--r-floor", type=float, default=1e-5)
    p.add_argument("--min-genes-per-tissue", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    # ─────────── Reference load + neuron downsampling ───────────
    print(f"[degrade] Loading reference {args.reference}")
    R = ad.read_h5ad(args.reference)
    if args.neuron_frac > 0:
        print(f"[degrade] Downsampling neurons to {args.neuron_frac:.0%}")
        R = _downsample_neurons(R, args.neuron_frac, args.seed,
                                label_col=args.ref_label_key)

    # ─────────── Anchor matching ───────────
    print("[degrade] Anchor matching")
    A_cap = ad.read_h5ad(args.capgta)
    A_cap.X = A_cap.layers[args.capgta_counts_layer].copy()
    anchor_map = (json.loads(args.anchor_json.read_text())
                  if args.anchor_json else DEFAULT_ANCHOR_CAPGTA)
    A_cap.obs["anchor_tissue"] = (
        A_cap.obs[args.capgta_cluster_key].astype(str).map(anchor_map)
    )
    A_cap = A_cap[A_cap.obs["anchor_tissue"].notna()].copy()
    A_cap.var_names = A_cap.var["wormbase_id"].astype(str).values
    if not A_cap.var_names.is_unique:
        A_cap = A_cap[:, ~A_cap.var_names.duplicated()].copy()
    print(f"  CapGTA anchor cells: {A_cap.n_obs}")

    ref_types = R.obs[args.ref_label_key].astype(str).unique().tolist()
    anchor_ref = _resolve_ref_anchor_types(DEFAULT_ANCHOR_REF, ref_types)
    tissue_of = {t: label for label, ts in anchor_ref.items() for t in ts}
    R.obs["anchor_tissue"] = R.obs[args.ref_label_key].astype(str).map(tissue_of)
    R_anchor = R[R.obs["anchor_tissue"].notna()].copy()
    print(f"  reference anchor cells: {R_anchor.n_obs}")

    anchors = sorted(set(A_cap.obs["anchor_tissue"].unique())
                     & set(R_anchor.obs["anchor_tissue"].unique()))
    print(f"  common anchors: {anchors}")

    # ─────────── Pseudobulk both sides ───────────
    pb_cap = _pseudobulk_proportions(A_cap, "anchor_tissue")
    pb_ref = _pseudobulk_proportions(R_anchor, "anchor_tissue")
    common = pb_cap.index.intersection(pb_ref.index)
    pb_cap = pb_cap.loc[common]; pb_ref = pb_ref.loc[common]
    print(f"[degrade] gene intersection: {len(common)}")

    # ─────────── GAM fit ───────────
    print("[degrade] Building fit frame")
    feats_all = _build_gene_features(common, args.kg_source, args.calibration)
    rows = []
    for t in anchors:
        o = pb_cap[t].values; r = pb_ref[t].values
        keep = (r > args.r_floor) & (o > 0)
        if keep.sum() < args.min_genes_per_tissue:
            print(f"  WARN: tissue '{t}' only {keep.sum()} genes above floor — skipping")
            continue
        sub = feats_all.loc[common[keep]].copy()
        sub["log_ratio"] = np.log(o[keep] / r[keep])
        rows.append(sub)
    fit_frame = pd.concat(rows, ignore_index=True)
    print(f"  fit_frame: {len(fit_frame)} rows")
    gam = _fit_gam(fit_frame)
    r2 = gam.statistics_["pseudo_r2"]["explained_deviance"]
    print(f"[degrade] pooled GAM R² = {r2:.3f}")

    # ─────────── p_g on all reference genes ───────────
    print("[degrade] Computing p_g on all reference genes")
    feats_ref = _build_gene_features(R.var_names, args.kg_source, args.calibration)
    p_g_raw = _predict_p_g(gam, feats_ref)
    p_g = _rescale_total_preserve(p_g_raw, R.X)
    print(f"  p_g mean={p_g.mean():.4f}, median={np.median(p_g):.4f}, "
          f"frac==1={(p_g >= 0.999).mean():.4f}")
    pd.DataFrame({
        "wormbase_id": R.var_names,
        "n_junctions": feats_ref["n_junctions"].values,
        "exon_bp": feats_ref["exon_bp"].values,
        "polyA_ind": feats_ref["polyA_ind"].values,
        "sl_ind": feats_ref["sl_ind"].values,
        "p_g": p_g,
    }).to_csv(out / "p_g_per_gene.csv", index=False)

    # ─────────── Binomial thinning × N seeds ───────────
    print(f"[degrade] Binomial thinning × {args.n_seeds} seeds")
    seed_dir = out / "degraded_seeds"; seed_dir.mkdir(exist_ok=True)
    X_ref = R.X.tocsr() if sp.issparse(R.X) else sp.csr_matrix(R.X)
    for s in range(args.n_seeds):
        rng_s = np.random.default_rng(args.seed + 100 + s)
        Xd = _binomial_thin(X_ref, p_g, rng_s)
        Ad = ad.AnnData(Xd, obs=R.obs.copy(), var=R.var.copy())
        Ad.var["p_g"] = p_g
        Ad.uns["degradation_config"] = json.dumps({
            "seed": args.seed + 100 + s,
            "neuron_frac": args.neuron_frac,
            "anchors": anchors,
            "n_genes_ref": int(R.n_vars),
            "gam_pooled_r2": float(r2),
            "rescale": "total_preserve",
        })
        path = seed_dir / f"taylor2020_degraded_seed{s}.h5ad"
        Ad.write_h5ad(path)
        pre = int(X_ref.sum()); post = int(Xd.sum())
        print(f"  seed {s}: {path.name} — total counts {pre:.2e} → {post:.2e} "
              f"(kept {post/pre:.1%})")

    with (out / "report.md").open("w") as f:
        f.write("# Reference degradation report\n\n")
        f.write(f"Reference: `{args.reference}`  \n")
        f.write(f"CapGTA (anchor cells): `{args.capgta}` — {A_cap.n_obs}  \n")
        f.write(f"Neuron downsampling target: {args.neuron_frac:.0%} "
                f"({R.n_obs} cells)  \n")
        f.write(f"Common anchor tissues: {anchors}  \n")
        f.write(f"Gene intersection: {len(common)}  \n")
        f.write(f"Pooled GAM pseudo-R²: **{r2:.3f}**  \n")
        f.write(f"p_g mean={p_g.mean():.3f}, median={np.median(p_g):.3f}, "
                f"frac ≥ 0.5 = {(p_g >= 0.5).mean():.3f}  \n\n")
        f.write("## Artifacts\n\n")
        f.write("- `p_g_per_gene.csv` — fitted p_g per reference gene\n")
        f.write(f"- `degraded_seeds/taylor2020_degraded_seed{{0..{args.n_seeds-1}}}.h5ad`\n")
    print(f"[degrade] Done. Report → {out / 'report.md'}")


if __name__ == "__main__":
    main()
