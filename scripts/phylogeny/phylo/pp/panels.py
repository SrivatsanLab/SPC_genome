"""Per-worm somatic panel construction.

For each passing worm:

1. Pool AD, DP over cells assigned to that worm to get per-worm bulk VAF.
2. Keep variants where own-worm VAF sits in the somatic band and every other
   worm's VAF is at or below ``other_worm_max_vaf``.
3. Require ``n_carriers >= min_carriers`` in the own worm.
4. Blacklist variants shared across ``>= blacklist_min_worms`` worms at
   ``bulk_vaf >= blacklist_vaf_min`` (recurrent artifact / unfiltered common
   variant).

Panel = self-contained AnnData subset per worm, written to
``run/panels/worm_{XX}.h5ad``. Downstream stages call :func:`load_worm_panel`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse

from ..utils.run import RunDir

# --- substitution-type canonicalization ------------------------------------

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
# Canonical form collapses to a C or T reference base.
_CANONICAL_REF = {"C", "T"}


def canonical_sub_type(ref: str, alt: str) -> str:
    """Return ``'C>A'``, ``'C>G'``, ..., collapsing to the C/T-referenced strand.

    ``G>T`` and ``C>A`` both return ``'C>A'`` (same substitution, opposite strands).
    """
    if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        return "?"
    if ref not in _CANONICAL_REF:
        ref, alt = _COMPLEMENT[ref], _COMPLEMENT[alt]
    return f"{ref}>{alt}"


def _sub_types_from_trinuc(anc: pd.Series, der: pd.Series) -> np.ndarray:
    """Extract canonical substitution type per variant from trinuc context strings.

    Expects three-letter strings like ``ATA`` / ``AGA``; the middle base is the
    substitution site.
    """
    ref = anc.astype(str).str[1]
    alt = der.astype(str).str[1]
    return np.array([canonical_sub_type(r, a) for r, a in zip(ref, alt)])


def _pool_ad_dp_per_worm(adata: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (AD_worm, DP_worm) as (n_variants × n_worms) DataFrames.

    Pooled read-level, per plan §HAPLOTYPE 2.1 (contamination-invariant).
    Only cells passing QC (``obs['qc_pass']``) contribute; cells whose donor
    fails ``worm_size_gte`` will simply produce a small-sum column.
    """
    obs = adata.obs
    donors = obs["donor"].astype(str).to_numpy()
    keep = obs["qc_pass"].to_numpy() if "qc_pass" in obs else np.ones(adata.n_obs, dtype=bool)
    donors_kept = donors[keep]
    unique = sorted(set(donors_kept))

    AD, DP = adata.layers["AD"], adata.layers["DP"]
    if not issparse(AD):
        AD = csr_matrix(AD)
    if not issparse(DP):
        DP = csr_matrix(DP)

    n_var = adata.n_vars
    ad_cols, dp_cols = {}, {}
    for w in unique:
        rows = np.where(keep & (donors == w))[0]
        if rows.size == 0:
            ad_cols[w] = np.zeros(n_var, dtype=np.int64)
            dp_cols[w] = np.zeros(n_var, dtype=np.int64)
            continue
        ad_cols[w] = np.asarray(AD[rows].sum(axis=0)).ravel().astype(np.int64)
        dp_cols[w] = np.asarray(DP[rows].sum(axis=0)).ravel().astype(np.int64)

    idx = adata.var_names
    return pd.DataFrame(ad_cols, index=idx), pd.DataFrame(dp_cols, index=idx)


def _compute_per_worm_vaf(
    adata: ad.AnnData,
    AD_w: pd.DataFrame,
    DP_w: pd.DataFrame,
    *,
    method: str = "pooled",
    target_dp: int | None = None,
) -> pd.DataFrame:
    """Per-worm bulk VAF matrix under a chosen normalization.

    Modes
    -----
    ``pooled``
        ``Σ AD / Σ DP`` — depth-weighted (current behavior). A few
        very-well-covered cells dominate.
    ``hypergeometric``
        Hypergeometric downsampling to ``target_dp`` via
        ``cellspec.utils.context.compute_vaf``. Reduces the influence of
        high-DP outliers; matches how ``.var['bulk_vaf']`` is computed at
        the pool level.
    ``colmean``
        Unweighted mean of per-cell VAFs (``AD/DP``) over cells with
        ``DP > 0``. Each covered cell contributes equally regardless of DP.
    """
    method = method.lower()
    if method == "pooled":
        vaf = AD_w.astype(float) / DP_w.replace(0, np.nan)
        return vaf.fillna(0.0)
    if method == "hypergeometric":
        try:
            from cellspec.utils.context import compute_vaf
        except ImportError as e:
            raise ImportError("vaf_normalize='hypergeometric' requires cellspec.") from e
        # Default target_dp: median positive DP per worm. Leaves most sites unchanged
        # and pulls only outlier-high-DP sites in; using min(positive DP) sends the
        # target to 1, making VAF a Bernoulli(true p) at every site.
        vaf = pd.DataFrame(index=AD_w.index, columns=AD_w.columns, dtype=np.float32)
        for worm in AD_w.columns:
            ad_col = AD_w[worm].to_numpy(dtype=np.int64)
            dp_col = DP_w[worm].to_numpy(dtype=np.int64)
            if target_dp is None:
                nz = dp_col[dp_col > 0]
                td = int(np.median(nz)) if nz.size else 1
            else:
                td = int(target_dp)
            vaf[worm] = compute_vaf(ad_col, dp_col, target_dp=td)
        return vaf.fillna(0.0)
    if method == "colmean":
        obs = adata.obs
        donors = obs["donor"].astype(str).to_numpy()
        keep = obs["qc_pass"].to_numpy() if "qc_pass" in obs else np.ones(adata.n_obs, dtype=bool)
        AD, DP = adata.layers["AD"], adata.layers["DP"]
        if not issparse(AD):
            AD = csr_matrix(AD)
        if not issparse(DP):
            DP = csr_matrix(DP)
        unique = sorted(set(donors[keep]))
        cols = {}
        for w in unique:
            rows = np.where(keep & (donors == w))[0]
            if rows.size == 0:
                cols[w] = np.zeros(adata.n_vars)
                continue
            AD_w_mat = AD[rows].toarray().astype(np.float32)
            DP_w_mat = DP[rows].toarray().astype(np.float32)
            covered = DP_w_mat > 0
            with np.errstate(invalid="ignore", divide="ignore"):
                per_cell = np.where(covered, AD_w_mat / np.maximum(DP_w_mat, 1), np.nan)
            cols[w] = np.nanmean(per_cell, axis=0)
            cols[w] = np.where(np.isnan(cols[w]), 0.0, cols[w])
        return pd.DataFrame(cols, index=adata.var_names)
    raise ValueError(f"unknown vaf_normalize={method!r} (want pooled/hypergeometric/colmean)")


def _n_carriers_per_worm(adata: ad.AnnData) -> pd.DataFrame:
    """Cells-with-≥1-alt-read count per (variant × worm)."""
    obs = adata.obs
    donors = obs["donor"].astype(str).to_numpy()
    keep = obs["qc_pass"].to_numpy() if "qc_pass" in obs else np.ones(adata.n_obs, dtype=bool)
    AD = adata.layers["AD"]
    if not issparse(AD):
        AD = csr_matrix(AD)

    unique = sorted(set(donors[keep]))
    cols = {}
    for w in unique:
        rows = np.where(keep & (donors == w))[0]
        if rows.size == 0:
            cols[w] = np.zeros(adata.n_vars, dtype=np.int32)
        else:
            cols[w] = np.asarray((AD[rows] > 0).sum(axis=0)).ravel().astype(np.int32)
    return pd.DataFrame(cols, index=adata.var_names)


def build_panels(adata: ad.AnnData, cfg: dict[str, Any], run: RunDir) -> dict[str, Path]:
    """Build per-worm somatic panels.

    Requires ``adata.obs['qc_pass']`` (from :func:`phylo.pp.apply_qc`).

    Outputs (under ``run/panels/``):
      - ``worm_{donor}.h5ad`` per passing worm (AnnData subset, self-contained)
      - ``summary.csv`` — one row per worm
      - ``blacklist.tsv`` — cross-worm shared variants
      - ``cells.csv`` — per-cell QC state
      - ``per_worm_bulk_vaf.parquet`` — variant × worm VAF matrix (kept for §6.1 amplitude check)

    Returns a dict of ``{worm_id: h5ad_path}``.
    """
    p_cfg = cfg["panels"]
    band_raw = p_cfg["somatic_vaf_band"]
    # Nullable bounds: None → no bound on that side. Written back as a two-tuple
    # of (lo, hi) with -inf/+inf sentinels for clean comparisons.
    lo = -np.inf if band_raw[0] is None else float(band_raw[0])
    hi = np.inf if band_raw[1] is None else float(band_raw[1])
    band = (lo, hi)
    other_max = float(p_cfg["other_worm_max_vaf"])
    min_carriers = int(p_cfg["min_carriers"])
    bl_min_worms = int(p_cfg["blacklist_min_worms"])
    bl_vaf_min = float(p_cfg["blacklist_vaf_min"])
    min_dp_per_worm = int(p_cfg.get("min_bulk_dp_per_worm", 5))
    vaf_normalize = str(p_cfg.get("vaf_normalize") or "pooled").lower()
    target_dp = p_cfg.get("target_dp")                              # for hypergeometric
    excluded_types = set(p_cfg.get("exclude_mutation_types") or [])
    # Canonicalize user-supplied types (accept G>T as C>A alias etc.).
    excluded_types = {canonical_sub_type(s.split(">")[0], s.split(">")[1]) for s in excluded_types}
    # Trinuc-context whitelist: variants matching one of these `anc>der` 3-mer
    # pairs are retained even if their substitution type is in `excluded_types`.
    # E.g. exclude_mutation_types=[C>A], retain_trinuc_contexts=[TCT>TAT] keeps
    # the pole-1 hotspot while dropping the rest of the 8-oxoG-prone C>A pool.
    retain_contexts = set(str(s).upper() for s in (p_cfg.get("retain_trinuc_contexts") or []))

    out = run.subdir("panels")

    if "qc_pass" not in adata.obs:
        raise ValueError("Run QC first: adata.obs['qc_pass'] is missing.")

    adata.obs.reset_index().rename(columns={"index": "cell_id"})[
        ["cell_id", "donor", "purity", "alpha_hap", "eff_modules", "capsule_type", "qc_pass", "qc_reason"]
    ].to_csv(out / "cells.csv", index=False)

    # Substitution-type filter (§10.9 A of TREE_FORWARD_STATUS.md).
    if excluded_types:
        sub_types = _sub_types_from_trinuc(adata.var["anc"], adata.var["der"])
        type_pass = ~pd.Series(sub_types).isin(excluded_types).to_numpy()
        if retain_contexts:
            trinuc_pairs = (adata.var["anc"].astype(str).str.upper() + ">"
                            + adata.var["der"].astype(str).str.upper()).to_numpy()
            retain_mask = np.isin(trinuc_pairs, list(retain_contexts))
            type_pass = type_pass | retain_mask
            n_retained = int(retain_mask.sum())
        else:
            n_retained = 0
        n_before = adata.n_vars
        adata = adata[:, type_pass].copy()
        msg = (f"[panels] excluded {list(sorted(excluded_types))}: "
               f"{n_before} → {adata.n_vars} variants")
        if retain_contexts:
            msg += f" (retained {n_retained} matching {sorted(retain_contexts)})"
        print(msg)

    print(f"[panels] pooling AD/DP per worm (vaf_normalize={vaf_normalize})…")
    AD_w, DP_w = _pool_ad_dp_per_worm(adata)
    VAF_w = _compute_per_worm_vaf(
        adata, AD_w, DP_w, method=vaf_normalize, target_dp=target_dp,
    )

    # Persist per-worm bulk AD/DP/VAF as a tiny AnnData (obs=variants, var=worms).
    per_worm = ad.AnnData(
        X=VAF_w.to_numpy().astype(np.float32),
        obs=adata.var[["chrom", "pos", "anc", "der", "bulk_vaf"]].copy(),
        var=pd.DataFrame(index=VAF_w.columns),
        layers={"AD": AD_w.to_numpy().astype(np.int64), "DP": DP_w.to_numpy().astype(np.int64)},
    )
    per_worm.write_h5ad(out / "per_worm_bulk.h5ad")

    N_w = _n_carriers_per_worm(adata)

    # Cross-worm blacklist: variants with bulk_vaf ≥ bl_vaf_min in ≥ bl_min_worms worms.
    shared = (VAF_w >= bl_vaf_min).sum(axis=1)
    bl_mask = shared >= bl_min_worms
    blacklist = adata.var.loc[bl_mask.values, ["chrom", "pos", "anc", "der"]].copy()
    blacklist["n_worms_sharing"] = shared[bl_mask].values
    blacklist.to_csv(out / "blacklist.tsv", sep="\t")
    print(f"[panels] blacklist: {int(bl_mask.sum())} variants shared in ≥{bl_min_worms} worms")

    summary_rows = []
    written: dict[str, Path] = {}

    for worm in VAF_w.columns:
        n_cells = int(((adata.obs["donor"].astype(str) == worm) & adata.obs["qc_pass"]).sum())
        if n_cells < int(cfg["qc"].get("min_cells_per_worm_for_panel", 30)):
            print(f"[panels] {worm}: skipped (n_cells={n_cells})")
            continue

        own_vaf = VAF_w[worm].to_numpy()
        own_dp = DP_w[worm].to_numpy()
        others = [c for c in VAF_w.columns if c != worm]
        other_max_vaf = VAF_w[others].to_numpy().max(axis=1)

        panel_mask = (
            (own_vaf >= band[0])
            & (own_vaf < band[1])
            & (other_max_vaf <= other_max)
            & (own_dp >= min_dp_per_worm)
            & (N_w[worm].to_numpy() >= min_carriers)
            & (~bl_mask.to_numpy())
        )
        n_panel = int(panel_mask.sum())

        cell_mask = (adata.obs["donor"].astype(str) == worm).to_numpy() & adata.obs["qc_pass"].to_numpy()
        sub = adata[cell_mask, panel_mask].copy()
        sub.var["own_bulk_vaf"] = own_vaf[panel_mask]
        sub.var["own_bulk_dp"] = own_dp[panel_mask]
        sub.var["other_worm_max_vaf"] = other_max_vaf[panel_mask]
        sub.var["own_n_carriers"] = N_w[worm].to_numpy()[panel_mask]
        sub.uns["panel_worm"] = worm
        sub.uns["panel_config"] = {
            # AnnData .uns can't hold None; write nullable bounds as strings.
            "somatic_vaf_band_lo": ("null" if band_raw[0] is None else float(band_raw[0])),
            "somatic_vaf_band_hi": ("null" if band_raw[1] is None else float(band_raw[1])),
            "vaf_normalize": vaf_normalize,
            "target_dp": ("null" if target_dp is None else int(target_dp)),
            "excluded_mutation_types": sorted(excluded_types) if excluded_types else "none",
            "other_worm_max_vaf": other_max,
            "min_carriers": min_carriers,
            "min_bulk_dp_per_worm": min_dp_per_worm,
            "blacklist_min_worms": bl_min_worms,
            "blacklist_vaf_min": bl_vaf_min,
        }

        path = out / f"worm_{worm}.h5ad"
        sub.write_h5ad(path)
        written[worm] = path

        summary_rows.append(
            {
                "worm": worm,
                "n_cells": n_cells,
                "n_panel_variants": n_panel,
                "med_own_vaf": float(np.median(own_vaf[panel_mask])) if n_panel else np.nan,
                "med_own_dp": float(np.median(own_dp[panel_mask])) if n_panel else np.nan,
            }
        )
        print(f"[panels] {worm}: {n_cells} cells × {n_panel} panel variants → {path.name}")

    if summary_rows:
        summary = pd.DataFrame(summary_rows).sort_values("n_cells", ascending=False)
    else:
        summary = pd.DataFrame(columns=["worm", "n_cells", "n_panel_variants", "med_own_vaf", "med_own_dp"])
        print("[panels] WARNING: no worm produced a panel. Loosen QC or lower min_cells_per_worm_for_panel.")
    summary.to_csv(out / "summary.csv", index=False)
    print(f"[panels] {len(written)} worm panels written → {out}")
    return written


def load_worm_panel(run: RunDir, worm: str) -> ad.AnnData:
    """Load a per-worm panel written by :func:`build_panels`."""
    p = run.path / "panels" / f"worm_{worm}.h5ad"
    if not p.exists():
        raise FileNotFoundError(f"No panel for {worm} at {p}")
    return ad.read_h5ad(p)
