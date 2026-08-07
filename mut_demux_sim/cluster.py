"""
Demultiplexing: site selection, pairwise distances, clustering, assignment.

All pairwise quantities are computed *pairwise-complete* - restricted to sites
covered in both cells of a pair. This is the single most important detail: if
uncovered entries are treated as reference, distances track sequencing depth
instead of genotype, which is the classic failure mode for this problem.

Everything is vectorised via matrix products; there are no per-pair loops.
"""
from __future__ import annotations
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

__all__ = [
    "colsum", "carrier_fraction", "select_sites", "binary_matrices",
    "abcd", "dist_hamming", "dist_jaccard", "sim_phi", "corr_masked",
    "upgma", "consensus", "assign_em", "demultiplex", "within_cross",
]


def colsum(X):
    """Column sums as a 1-D array, for dense or scipy-sparse input."""
    return np.asarray(X.sum(0)).ravel()


def carrier_fraction(AD, DP, min_dp=1, min_alt=1):
    """Per-site (carrier fraction, n_carriers, n_covered).

    carrier fraction = cells with an alt read / cells with coverage. For a
    homozygous germline variant carried by j of K worms this is ~ j/K, which
    is what makes the histogram of this quantity diagnostic (see script 06).
    """
    cov = DP >= min_dp
    n_cov = colsum(cov)
    n_alt = colsum((AD >= min_alt).multiply(cov) if hasattr(AD, "multiply")
                   else (AD >= min_alt) & cov)
    return n_alt / np.maximum(n_cov, 1), n_alt, n_cov


def select_sites(AD, DP, frac_lo=0.03, frac_hi=0.97, min_cells=2,
                 max_sites=None, min_dp=1, min_alt=1, rng=None):
    """Indices of informative sites, filtered on carrier fraction.

    frac_lo excludes somatic clones and error; frac_hi excludes founder /
    strain-background sites fixed across the pool. max_sites subsamples the
    survivors (highest carrier count first) to bound memory.
    """
    frac, n_alt, n_cov = carrier_fraction(AD, DP, min_dp, min_alt)
    keep = (n_alt >= min_cells) & (frac > frac_lo) & (frac < frac_hi)
    idx = np.flatnonzero(keep)
    if max_sites is not None and len(idx) > max_sites:
        idx = idx[np.argsort(-n_alt[idx])[:max_sites]]
        idx.sort()
    return idx


def binary_matrices(AD, DP, idx, min_dp=1, min_alt=1, dtype=np.float32):
    """(A, M) presence and coverage matrices over the selected sites.

    A[i, s] = 1 if cell i shows an alt read at site s (and is covered)
    M[i, s] = 1 if cell i is covered at site s
    """
    ADs, DPs = AD[:, idx], DP[:, idx]
    if hasattr(ADs, "toarray"):
        ADs, DPs = ADs.toarray(), DPs.toarray()
    M = (DPs >= min_dp)
    A = (ADs >= min_alt) & M
    return A.astype(dtype), M.astype(dtype)


def abcd(A, M):
    """Per-pair 2x2 contingency counts over jointly covered sites.

    a = both present, b = present in i only, c = present in j only,
    d = both absent, n = jointly covered.
    """
    a = A @ A.T
    n = M @ M.T
    b = A @ M.T - a
    c = M @ A.T - a
    d = n - a - b - c
    return a, b, c, d, n


def _finish_distance(D, n, min_shared=0):
    D = np.asarray(D, dtype=np.float64)
    if min_shared:
        D[n < min_shared] = 1.0
    D[n == 0] = 1.0
    D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    return D


def dist_hamming(A, M, min_shared=0):
    """Fraction of jointly covered sites at which the two cells disagree."""
    a, b, c, d, n = abcd(A, M)
    return _finish_distance((b + c) / np.maximum(n, 1), n, min_shared)


def dist_jaccard(A, M, min_shared=0):
    """1 - |both present| / |present in either|. Ignores co-absence."""
    a, b, c, d, n = abcd(A, M)
    return _finish_distance((b + c) / np.maximum(a + b + c, 1e-9), n, min_shared)


def sim_phi(A, M, min_shared=0):
    """Phi coefficient (Pearson on binary calls), pairwise-complete.

    Uses co-absence but discounts it against the marginal frequencies, so a
    shared rare variant counts for more than a shared common absence.
    Returns a *similarity* in [-1, 1]; convert with D = 1 - S.
    """
    a, b, c, d, n = abcd(A, M)
    den = np.sqrt((a + b) * (c + d) * (a + c) * (b + d))
    S = np.divide(a * d - b * c, den, out=np.zeros_like(den), where=den > 0)
    if min_shared:
        S[n < min_shared] = 0.0
    S[n == 0] = 0.0
    S = np.clip(0.5 * (S + S.T), -1, 1)
    np.fill_diagonal(S, 1.0)
    return S


def corr_masked(AD, DP, idx, min_dp=1, min_shared=0, binary=True):
    """Pairwise-complete correlation of per-variant-centred values.

    binary=True centres presence/absence calls (right for homozygous worms);
    binary=False centres alt fraction AD/DP (right when heterozygotes carry
    information, e.g. outbred donors).

    Missing entries are set to the centred mean (0), so they drop out of the
    numerator automatically; the mask restores the correct per-pair norms.
    """
    if binary:
        A, M = binary_matrices(AD, DP, idx, min_dp=min_dp)
        num = A
    else:
        ADs, DPs = AD[:, idx], DP[:, idx]
        if hasattr(ADs, "toarray"):
            ADs, DPs = ADs.toarray(), DPs.toarray()
        M = (DPs >= min_dp).astype(np.float32)
        with np.errstate(invalid="ignore", divide="ignore"):
            num = np.where(M > 0, ADs / np.maximum(DPs, 1), 0).astype(np.float32)
    colmean = num.sum(0) / np.maximum(M.sum(0), 1)
    X = np.where(M > 0, num - colmean, 0).astype(np.float32)
    n = M @ M.T
    den = np.sqrt(((X * X) @ M.T) * (M @ (X * X).T))
    C = np.divide(X @ X.T, den, out=np.zeros_like(den), where=den > 0)
    if min_shared:
        C[n < min_shared] = 0.0
    C = np.clip(0.5 * (C + C.T), -1, 1).astype(np.float64)
    np.fill_diagonal(C, 1.0)
    return C, n


def upgma(D, K, method="average"):
    """Average-linkage clustering of a square distance matrix, cut into K groups.

    'average' (UPGMA) is the default because 1 - correlation is not a Euclidean
    metric, so Ward is unjustified; and single linkage chains across the
    intermediate cells that bridge clusters.
    """
    Z = linkage(squareform(D, checks=False), method=method)
    return fcluster(Z, K, "maxclust"), Z


def consensus(AD, DP, idx, labels, eps=0.01):
    """Pooled alt fraction per cluster per site (depth-summed, so high depth)."""
    ADs, DPs = AD[:, idx], DP[:, idx]
    if hasattr(ADs, "toarray"):
        ADs, DPs = ADs.toarray(), DPs.toarray()
    ids = np.unique(labels)
    Q = np.stack([
        ADs[labels == c].sum(0) / np.maximum(DPs[labels == c].sum(0), 1)
        for c in ids
    ])
    return ids, np.clip(Q, eps, 1 - eps)


def assign_em(AD, DP, idx, Q, n_iter=25, ambient=True, eps=0.01):
    """Refine consensus genotypes and assign every cell by binomial likelihood.

    Handles low-coverage cells properly: a site with one read contributes one
    read's worth of evidence rather than being thresholded away. Optionally
    includes a fixed 'ambient' component at the pooled allele frequency, which
    absorbs capsules containing DNA from many individuals.
    """
    ADs, DPs = AD[:, idx], DP[:, idx]
    if hasattr(ADs, "toarray"):
        ADs, DPs = ADs.toarray(), DPs.toarray()
    ad = ADs.astype(np.float32)
    dp = DPs.astype(np.float32)
    ref = dp - ad
    pool = np.clip(ad.sum(0) / np.maximum(dp.sum(0), 1), eps, 1 - eps)

    for _ in range(n_iter):
        Qa = np.vstack([Q, pool]) if ambient else Q
        L = ad @ np.log(Qa).T + ref @ np.log(1 - Qa).T
        R = np.exp(L - L.max(1, keepdims=True))
        R /= R.sum(1, keepdims=True)
        W = R[:, :len(Q)]
        Q = np.clip((W.T @ ad) / np.maximum(W.T @ dp, 1e-9), eps, 1 - eps)

    LL = ad @ np.log(Q).T + ref @ np.log(1 - Q).T
    out = dict(Q=Q, LL=LL, labels=LL.argmax(1), n_reads=dp.sum(1))
    if ambient:
        out["ll_ambient"] = ad @ np.log(pool) + ref @ np.log(1 - pool)
    return out


def demultiplex(AD, DP, K, frac_lo=0.03, frac_hi=0.97, min_cells=2,
                max_sites=20000, metric="hamming", min_shared=0,
                refine_em=False, n_iter=25):
    """End-to-end: select sites -> pairwise distance -> UPGMA -> (optional EM).

    metric: 'hamming' | 'jaccard' | 'phi' | 'corr'
    """
    idx = select_sites(AD, DP, frac_lo, frac_hi, min_cells, max_sites)
    if len(idx) < 10:
        raise ValueError(f"only {len(idx)} sites survived filtering; "
                         "widen frac_lo/frac_hi or lower min_cells")
    if metric == "corr":
        C, n = corr_masked(AD, DP, idx, min_shared=min_shared)
        D = 1 - C
        np.fill_diagonal(D, 0.0)
        S = C
    else:
        A, M = binary_matrices(AD, DP, idx)
        n = M @ M.T
        if metric == "hamming":
            D = dist_hamming(A, M, min_shared); S = 1 - D
        elif metric == "jaccard":
            D = dist_jaccard(A, M, min_shared); S = 1 - D
        elif metric == "phi":
            S = sim_phi(A, M, min_shared); D = 1 - S
            np.fill_diagonal(D, 0.0)
        else:
            raise ValueError(f"unknown metric: {metric}")
    labels, Z = upgma(D, K)
    res = dict(labels=labels, linkage=Z, distance=D, similarity=S,
               n_shared=n, sites=idx)
    if refine_em:
        ids, Q = consensus(AD, DP, idx, labels)
        em = assign_em(AD, DP, idx, Q, n_iter=n_iter, ambient=False)
        res["labels_em"] = ids[em["labels"]]
        res["em"] = em
    return res


def within_cross(S, truth):
    """Median within-group and between-group similarity, and their gap.

    The absolute within-group value is far less important than this gap:
    clustering uses the separation, not the level.
    """
    off = ~np.eye(len(S), dtype=bool)
    same = (truth[:, None] == truth[None, :]) & off
    diff = (truth[:, None] != truth[None, :]) & off
    w, b = np.median(S[same]), np.median(S[diff])
    return dict(within=float(w), cross=float(b), gap=float(w - b))
