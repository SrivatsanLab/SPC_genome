"""G0 diagnostics: μ_div, within-vs-cross-worm sharing, RNA budget.

These are the feasibility gates in
``docs/T2_tree_validation_plan.md`` §G0.1–G0.4 (and PLAN_AMENDMENTS.md §6.G0.1).
Each function is standalone and can be called from a notebook; :func:`run_g0`
is the driver used by the CLI.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse

from ..utils.run import RunDir


# --- G0.1: burden -----------------------------------------------------------


def compute_burden(
    adata: ad.AnnData,
    *,
    somatic_vaf_band: tuple[float, float] = (0.02, 0.1),
    inplace: bool = True,
) -> pd.DataFrame:
    """Per-cell somatic-detection count → burden → μ_div.

    Filters variants to the somatic VAF band (default (0.02, 0.1) — the mid
    band from ``soup_signature/spectra_by_vaf_and_region.csv`` where C>A
    damage has fallen to ~12%). Counts non-zero AD entries per cell.

    ``burden ≈ detections / (b * 0.5 * (1 - alpha_hap))`` per plan §G0.1;
    ``mu_div ≈ median(burden) / n_divisions``.

    Parameters
    ----------
    adata : AnnData
        Must have ``.layers['AD']``, ``.layers['DP']``, ``.var['bulk_vaf']``,
        ``.obs['alpha_hap']``, ``.obs['sites_covered']``, ``.obs['donor']``.
    somatic_vaf_band : (lo, hi)
        Bulk VAF band treated as somatic.
    inplace : bool
        Write ``.obs['somatic_detections']`` and ``.obs['burden_est']``.

    Returns
    -------
    DataFrame
        Per-cell: ``cell_id, donor, alpha_hap, b_i, detections, burden``.
    """
    lo, hi = somatic_vaf_band
    var_mask = (adata.var["bulk_vaf"].to_numpy() >= lo) & (adata.var["bulk_vaf"].to_numpy() < hi)
    AD = adata.layers["AD"]
    if issparse(AD):
        AD_sub = AD[:, var_mask]
        detections = np.asarray((AD_sub > 0).sum(axis=1)).ravel()
    else:
        detections = (AD[:, var_mask] > 0).sum(axis=1)

    alpha = adata.obs["alpha_hap"].to_numpy()
    # coverage breadth b_i: fraction of variants with DP>=1 for this cell.
    # sites_covered is genome-wide; here we want the panel breadth.
    DP = adata.layers["DP"]
    DP_sub = DP[:, var_mask]
    if issparse(DP_sub):
        b_i = np.asarray((DP_sub > 0).sum(axis=1)).ravel() / max(int(var_mask.sum()), 1)
    else:
        b_i = (DP_sub > 0).sum(axis=1) / max(int(var_mask.sum()), 1)

    denom = np.clip(b_i * 0.5 * (1.0 - alpha), 1e-6, None)
    burden = detections / denom

    out = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "donor": adata.obs["donor"].astype(str).to_numpy(),
            "alpha_hap": alpha,
            "b_i": b_i,
            "detections": detections,
            "burden": burden,
        }
    )
    if inplace:
        adata.obs["somatic_detections"] = detections
        adata.obs["burden_est"] = burden
    return out


# --- G0.2: within vs cross-worm variant sharing -----------------------------


def compute_within_vs_cross_sharing(
    adata: ad.AnnData,
    *,
    somatic_vaf_band: tuple[float, float] = (0.02, 0.1),
    max_pairs: int = 5000,
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:
    """Rate of shared non-germline alt calls, within vs across worms.

    For random pairs of cells, restrict to co-covered somatic-band sites and
    compute ``sum(AD_i * AD_j > 0) / co_covered``. The cross-worm mean sets
    the recurrent-artifact floor; excess within-worm is the lineage budget.

    Parameters
    ----------
    max_pairs : int
        Cap total pair evaluations (per side). 5,000 pairs at ~10k sites is
        sub-second.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    lo, hi = somatic_vaf_band
    var_mask = (adata.var["bulk_vaf"].to_numpy() >= lo) & (adata.var["bulk_vaf"].to_numpy() < hi)
    AD = adata.layers["AD"][:, var_mask]
    DP = adata.layers["DP"][:, var_mask]
    if not issparse(AD):
        AD = csr_matrix(AD)
    if not issparse(DP):
        DP = csr_matrix(DP)
    AD_bool = AD.astype(bool).astype(np.int8)
    DP_bool = (DP > 0).astype(np.int8)

    donors = adata.obs["donor"].astype(str).to_numpy()
    n = adata.n_obs

    def _pair_stat(i: int, j: int) -> tuple[float, int]:
        co = DP_bool[i].multiply(DP_bool[j]).sum()
        if co == 0:
            return (np.nan, 0)
        shared = AD_bool[i].multiply(AD_bool[j]).sum()
        return (shared / co, int(co))

    within, cross = [], []
    within_n = min(max_pairs, n * (n - 1) // 2)
    for _ in range(within_n):
        i, j = rng.integers(0, n, size=2)
        if i == j:
            continue
        if donors[i] == donors[j]:
            s, co = _pair_stat(int(i), int(j))
            if not np.isnan(s):
                within.append((donors[i], s, co))
    for _ in range(max_pairs):
        i, j = rng.integers(0, n, size=2)
        if i == j:
            continue
        if donors[i] != donors[j]:
            s, co = _pair_stat(int(i), int(j))
            if not np.isnan(s):
                cross.append((f"{donors[i]}|{donors[j]}", s, co))

    df_w = pd.DataFrame(within, columns=["worm", "shared_frac", "co_covered"])
    df_c = pd.DataFrame(cross, columns=["worm_pair", "shared_frac", "co_covered"])
    df_w["side"] = "within"
    df_c["side"] = "cross"
    df_c = df_c.rename(columns={"worm_pair": "worm"})
    out = pd.concat([df_w, df_c], ignore_index=True)
    return out


# --- driver -----------------------------------------------------------------


def run_g0(adata: ad.AnnData, cfg: dict[str, Any], run: RunDir) -> None:
    """Run G0 diagnostics and write outputs to ``run/g0/``.

    Runs on the QC-passing subset (``adata.obs['qc_pass']``); background
    capsules and low-purity cells otherwise dominate the detection counts
    and inflate the sharing rate.

    Writes:
      - ``g0/burden_per_cell.csv``
      - ``g0/within_vs_cross_worm_sharing.csv``
      - ``g0/summary.txt`` (headline μ_div, within/cross median)
    """
    out = run.subdir("g0")
    band_raw = cfg["panels"]["somatic_vaf_band"]
    band = (
        -np.inf if band_raw[0] is None else float(band_raw[0]),
        np.inf if band_raw[1] is None else float(band_raw[1]),
    )
    n_div = float(cfg["g0"].get("n_divisions_l4_soma", 10.0))

    if "qc_pass" in adata.obs:
        n_before = adata.n_obs
        adata = adata[adata.obs["qc_pass"].to_numpy()].copy()
        print(f"[g0] restricted to QC-pass cells: {adata.n_obs}/{n_before}")

    print(f"[g0] burden: somatic VAF band = {band}")
    burden = compute_burden(adata, somatic_vaf_band=band, inplace=False)
    burden.to_csv(out / "burden_per_cell.csv", index=False)
    mu_div = float(np.nanmedian(burden["burden"])) / n_div

    print("[g0] within-vs-cross-worm sharing")
    share = compute_within_vs_cross_sharing(
        adata,
        somatic_vaf_band=band,
        max_pairs=int(cfg["g0"].get("sharing_max_pairs", 5000)),
        rng=np.random.default_rng(int(cfg["g0"].get("seed", 0))),
    )
    share.to_csv(out / "within_vs_cross_worm_sharing.csv", index=False)
    within_med = float(share.query("side == 'within'")["shared_frac"].median())
    cross_med = float(share.query("side == 'cross'")["shared_frac"].median())

    with open(out / "summary.txt", "w") as f:
        f.write(f"somatic_vaf_band       {band}\n")
        f.write(f"n_cells                {adata.n_obs}\n")
        f.write(f"n_variants_in_band     {int(((adata.var['bulk_vaf'] >= band[0]) & (adata.var['bulk_vaf'] < band[1])).sum())}\n")
        f.write(f"median_burden          {np.nanmedian(burden['burden']):.1f}\n")
        f.write(f"median_mu_div          {mu_div:.2f}   (assuming n_divisions={n_div})\n")
        f.write(f"within_med_shared_frac {within_med:.4f}\n")
        f.write(f"cross_med_shared_frac  {cross_med:.4f}\n")
        f.write(f"within/cross ratio     {within_med / max(cross_med, 1e-9):.2f}\n")
        gate_pass = mu_div >= float(cfg["g0"].get("min_mu_div", 3.0)) and within_med > cross_med
        f.write(f"g0_gate_pass           {gate_pass}\n")

    print(f"[g0] μ_div ≈ {mu_div:.2f}; within/cross = {within_med:.4f} / {cross_med:.4f}")
    print(f"[g0] outputs → {out}")
