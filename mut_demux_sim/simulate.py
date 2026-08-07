"""
Simulator for pooled single-cell WGS of selfing C. elegans lines.

Models the design in question: N worms, all selfing descendants of one
engineered founder, pooled and sequenced at low coverage. Because selfing
drives heterozygosity to ~0, worms are treated as homozygous (genotype 0/1
per site) unless `het_frac` is set.

Four classes of site are generated, distinguished by *carrier fraction*
(fraction of covered cells showing an alt read):

  founder   present in every cell                    carrier frac ~ 1.0
  germline  present in all cells of >=1 worm         carrier frac ~ j/K
  somatic   present in a clonal subset of one worm   carrier frac < 1/K
  (errors)  false-positive calls, sprinkled          carrier frac ~ fp_rate

The germline sites are the demultiplexing signal; somatic sites are the
lineage-tracing signal and act as noise for this task.
"""
from __future__ import annotations
import numpy as np

__all__ = ["simulate_pool", "lineage_blocks"]


def lineage_blocks(cell_idx, rng, min_size=2):
    """Recursively split cells into a binary lineage; return clade index arrays.

    Each returned array is the set of cells descending from one node, i.e. the
    cells that would share a somatic mutation arising on that branch.
    """
    blocks = []

    def rec(idx):
        if len(idx) < min_size:
            return
        blocks.append(idx)
        perm = rng.permutation(idx)
        h = len(idx) // 2
        rec(perm[:h])
        rec(perm[h:])

    rec(np.asarray(cell_idx))
    return blocks


def simulate_pool(
    n_worms=20,
    cells_per_worm=150,
    germline_per_worm=300,
    somatic_per_worm=0,
    n_founder=0,
    sibling_shared_frac=0.0,
    breadth=0.25,
    breadth_heterogeneous=False,
    n_high_coverage=0,
    mean_depth=0.8,
    err=0.002,
    fp_rate=0.0,
    het_frac=0.0,
    seed=0,
    dtype=np.int16,
):
    """Simulate a pooled multi-worm single-cell WGS experiment.

    Parameters
    ----------
    n_worms, cells_per_worm
        Pool composition. Total cells = n_worms * cells_per_worm.
    germline_per_worm
        De novo germline mutations accumulated per lineage since the founder.
        This is the parameter that sets feasibility (see script 01).
    somatic_per_worm
        Total somatic mutations per worm, spread over a binary cell lineage so
        that clones of every size carry some. 0 disables somatic mutations.
    n_founder
        Sites fixed in every worm (strain background vs the reference). These
        carry no information; included to test that filtering them is optional.
    sibling_shared_frac
        If > 0, worms are paired and each pair shares this fraction of its
        germline mutations - i.e. recent common ancestry. 1.0 makes the pair
        genetically identical and unresolvable.
    breadth
        Mean fraction of sites covered per cell.
    breadth_heterogeneous
        If True, per-cell breadth is drawn from a right-skewed Beta rather than
        being constant, reproducing the wide coverage range of real data.
    n_high_coverage
        Number of cells forced to 0.80-0.95 breadth (only if heterogeneous).
    mean_depth
        Poisson mean for read depth at covered sites; actual depth is
        Poisson(mean_depth) + 1, so >=1 read wherever covered.
    err
        Per-read base error: P(alt read | non-carrier) and P(ref read | carrier).
    fp_rate
        Per-cell-per-site probability of a spurious alt call (WGA/caller
        artifact). Applied on top of the genotype.
    het_frac
        Fraction of germline mutations not yet fixed by selfing. These are
        heterozygous, so each read is a coin flip - a strong attenuator of
        within-worm similarity at low depth.
    seed
        RNG seed.

    Returns
    -------
    dict with keys:
        AD, DP     : (cells x sites) int arrays, alt depth and total depth
        truth      : (cells,) worm index for each cell
        site_kind  : (sites,) one of 'founder' / 'germline' / 'somatic'
        is_het     : (sites,) bool, unfixed germline sites
        breadth    : (cells,) realised per-cell breadth parameter
        params     : dict of the arguments used
    """
    rng = np.random.default_rng(seed)
    N = n_worms * cells_per_worm

    # ---- assemble site blocks: (cell_mask, n_sites, kind) ----
    blocks = []
    if n_founder:
        blocks.append((np.ones(N, bool), n_founder, "founder"))

    worm_cells = [np.arange(w * cells_per_worm, (w + 1) * cells_per_worm)
                  for w in range(n_worms)]

    # germline, with optional sibling-pair sharing
    if sibling_shared_frac > 0:
        n_shared = int(germline_per_worm * sibling_shared_frac)
        n_priv = germline_per_worm - n_shared
        for p in range(n_worms // 2):
            a, b = 2 * p, 2 * p + 1
            pair = np.zeros(N, bool)
            pair[worm_cells[a]] = True
            pair[worm_cells[b]] = True
            if n_shared:
                blocks.append((pair, n_shared, "germline"))
            for w in (a, b):
                if n_priv:
                    m = np.zeros(N, bool); m[worm_cells[w]] = True
                    blocks.append((m, n_priv, "germline"))
        if n_worms % 2:                      # odd worm out gets all private
            m = np.zeros(N, bool); m[worm_cells[-1]] = True
            blocks.append((m, germline_per_worm, "germline"))
    else:
        for w in range(n_worms):
            m = np.zeros(N, bool); m[worm_cells[w]] = True
            blocks.append((m, germline_per_worm, "germline"))

    # somatic: nested clones along a per-worm binary cell lineage
    if somatic_per_worm:
        for w in range(n_worms):
            clades = lineage_blocks(worm_cells[w], rng)
            clades = [c for c in clades if len(c) < cells_per_worm]  # exclude root
            if not clades:
                continue
            per = max(somatic_per_worm // len(clades), 1)
            for c in clades:
                m = np.zeros(N, bool); m[c] = True
                blocks.append((m, per, "somatic"))

    V = sum(b[1] for b in blocks)
    if N * V > 4e8:
        raise MemoryError(
            f"{N} cells x {V} sites is too large to materialise densely. "
            "Reduce cells_per_worm, germline_per_worm, or somatic_per_worm."
        )

    # ---- genotype matrix ----
    G = np.zeros((N, V), bool)
    site_kind = np.empty(V, dtype=object)
    is_het = np.zeros(V, bool)
    col = 0
    for mask, k, kind in blocks:
        G[mask, col:col + k] = True
        site_kind[col:col + k] = kind
        if kind == "germline" and het_frac:
            is_het[col:col + k] = rng.random(k) < het_frac
        col += k

    # ---- coverage ----
    if breadth_heterogeneous:
        b = np.clip(rng.beta(1.2, 4.0, N) * (breadth / 0.23), 0.01, 0.99)
        if n_high_coverage:
            hi = rng.choice(N, min(n_high_coverage, N), replace=False)
            b[hi] = rng.uniform(0.80, 0.95, len(hi))
    else:
        b = np.full(N, float(breadth))

    # ---- reads ----
    AD = np.zeros((N, V), dtype)
    DP = np.zeros((N, V), dtype)
    for i in range(N):
        covered = rng.random(V) < b[i]
        dp = np.where(covered, rng.poisson(mean_depth, V) + 1, 0)
        q = np.where(G[i], 1.0 - err, err)
        q = np.where(G[i] & is_het, 0.5, q)          # unfixed germline: coin flip
        if fp_rate:
            q = np.maximum(q, (rng.random(V) < fp_rate).astype(float))
        AD[i] = rng.binomial(dp, q)
        DP[i] = dp

    return dict(
        AD=AD, DP=DP,
        truth=np.repeat(np.arange(n_worms), cells_per_worm),
        site_kind=site_kind, is_het=is_het, breadth=b,
        params=dict(
            n_worms=n_worms, cells_per_worm=cells_per_worm,
            germline_per_worm=germline_per_worm, somatic_per_worm=somatic_per_worm,
            n_founder=n_founder, sibling_shared_frac=sibling_shared_frac,
            breadth=breadth, breadth_heterogeneous=breadth_heterogeneous,
            n_high_coverage=n_high_coverage, mean_depth=mean_depth, err=err,
            fp_rate=fp_rate, het_frac=het_frac, seed=seed,
        ),
    )
