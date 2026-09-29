"""Primitives shared by the demultiplexer and the tree builder.

Extracted verbatim from ``bin/haplotype_sweep.py`` so the K=20 pool sweep and
the K=2 recursive tree run identical code. The 8 bug retractions in
``HAPLOTYPE_ANALYSIS_SUMMARY.md`` §9 all lived in this loop; keeping one copy
guarantees they stay fixed once.

Three functions:

- :func:`build_matrix` — per-variant centering, uncovered entries impute to
  the variant mean (which is zero after centering, so they drop out).
- :func:`fit_modules` — truncated SVD → L2-normalized variant loadings →
  k-means → labels + centroid distances.
- :func:`score_cells` — depth-weighted pooled VAF per (cell × module) →
  argmax assignment + purity/ratio.
"""

from __future__ import annotations

import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD


def build_matrix(
    AD: np.ndarray,
    DP: np.ndarray,
    AF: np.ndarray,
    vaf_lo: float,
    vaf_hi: float,
    min_dp_var: int = 1,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Centered ``cells × variants`` matrix over a VAF window.

    Parameters
    ----------
    AD, DP : ndarray, shape (n_cells, n_variants)
    AF : ndarray, shape (n_variants,)
        Bulk VAF per variant (from :func:`cellspec.tl.compute_bulk_vaf`).
    vaf_lo, vaf_hi : float
        Exclusive bounds. Pass (-inf, inf) or (0, 1.1) to disable the filter.
    min_dp_var : int
        Minimum per-cell depth to treat a variant as observed.

    Returns
    -------
    X : ndarray, shape (n_cells, n_informative), float32
        Centered VAF matrix. Uncovered entries are exactly zero.
    inf_idx : ndarray, shape (n_informative,), int
        Column indices into the original variant axis.
    stats : dict
        ``n_var``, ``sites_per_cell``, ``cells_per_var``.
    """
    inf_idx = np.flatnonzero((AF > vaf_lo) & (AF < vaf_hi))
    ad_i = AD[:, inf_idx].astype(np.float32)
    dp_i = DP[:, inf_idx].astype(np.float32)
    m_i = dp_i >= min_dp_var
    M_i = m_i.astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        f_i = np.where(m_i, ad_i / np.maximum(dp_i, 1), np.nan)
    colmean = np.nansum(np.where(m_i, f_i, 0), 0) / np.maximum(M_i.sum(0), 1)
    X = np.where(m_i, f_i - colmean, 0.0).astype(np.float32)
    stats = dict(
        n_var=int(len(inf_idx)),
        sites_per_cell=float(np.median(M_i.sum(1))),
        cells_per_var=float(np.median(M_i.sum(0))),
    )
    return X, inf_idx, stats


def _l2_normalize(V: np.ndarray) -> np.ndarray:
    return V / np.maximum(np.linalg.norm(V, axis=1, keepdims=True), 1e-9)


def fit_modules(
    X: np.ndarray,
    ncomp: int,
    K: int,
    *,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Truncated SVD + K-means on L2-normalized variant loadings.

    Returns 1-based labels (matching haplotype_sweep convention) and
    per-variant distance to its assigned centroid.

    Parameters
    ----------
    X : ndarray, shape (n_cells, n_variants)
        Centered VAF matrix from :func:`build_matrix`.
    ncomp : int
        SVD components. Scree cliff sits at ``K - 1`` after centering
        (``HAPLOTYPE_ANALYSIS_SUMMARY.md`` §4).
    K : int
        Number of modules.

    Returns
    -------
    labels : ndarray, shape (n_variants,), int
        1-based cluster labels.
    dists : ndarray, shape (n_variants,), float
        L2 distance from each variant to its assigned centroid in the
        normalized loading space.
    """
    V = TruncatedSVD(n_components=ncomp, random_state=seed).fit(X).components_.T
    Vn = _l2_normalize(V)
    km = KMeans(K, n_init=10, random_state=0).fit(Vn)
    d = np.linalg.norm(Vn - km.cluster_centers_[km.labels_], axis=1)
    return km.labels_ + 1, d


def score_cells(
    AD: np.ndarray,
    DP: np.ndarray,
    var_idx: np.ndarray,
    labels: np.ndarray,
    min_dp_var: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Depth-weighted pooled VAF per cell per module.

    For cell *i* and module *j*::

        S[i, j] = sum(AD[i, v] for v in module_j) / sum(DP[i, v] for v in module_j)

    ``S`` rows are normalized to ``W`` (weights sum to 1). Assignment is
    ``argmax(W)``; purity is ``max(W)``.

    Parameters
    ----------
    AD, DP : ndarray, shape (n_cells, n_variants)
        Full matrices (not the informative subset).
    var_idx : ndarray, shape (n_informative,)
        Column indices into AD/DP for the module-carrying variants.
    labels : ndarray, shape (n_informative,)
        Module label per variant. Matches ``var_idx`` in order.
    min_dp_var : int
        Minimum per-cell depth to include a variant.

    Returns
    -------
    W : ndarray, shape (n_cells, n_modules), float
        Normalized weights (rows sum to ~1).
    modules : ndarray, shape (n_modules,)
        Unique module labels in the same column order as W.
    assignment : ndarray, shape (n_cells,)
        Module label with maximum weight per cell.
    """
    modules = np.unique(labels)
    S = np.zeros((AD.shape[0], len(modules)))
    for j, k in enumerate(modules):
        c = var_idx[labels == k]
        a = AD[:, c].astype(np.float64)
        dd = DP[:, c].astype(np.float64)
        ok = dd >= min_dp_var
        S[:, j] = (a * ok).sum(1) / np.maximum((dd * ok).sum(1), 1)
    W = S / np.maximum(S.sum(1, keepdims=True), 1e-9)
    return W, modules, modules[S.argmax(1)]
