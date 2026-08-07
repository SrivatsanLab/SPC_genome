"""
Diagnostics to run on real data before (and after) clustering.

These answer the questions that determine whether demultiplexing is feasible
at all, which no amount of parameter tuning can substitute for.
"""
from __future__ import annotations
import numpy as np

from .cluster import colsum, carrier_fraction, select_sites, binary_matrices

__all__ = ["site_summary", "shared_site_stats", "plot_carrier_fraction",
           "feasibility_report"]


def site_summary(AD, DP, min_dp=1, min_alt=1, frac_bands=None):
    """Counts of sites by carrier-fraction band.

    The bands separate the four classes of variant: error/somatic at the
    bottom, worm-level germline in the middle, founder background at the top.
    """
    frac, n_alt, n_cov = carrier_fraction(AD, DP, min_dp, min_alt)
    if frac_bands is None:
        frac_bands = [
            (0.0, 0.005, "somatic (small clones) / error"),
            (0.005, 0.03, "somatic (large clones) / small clades"),
            (0.03, 0.97, "germline, worm-level"),
            (0.97, 1.01, "founder / strain background"),
        ]
    ok = n_alt >= 2
    rows = [(lbl, int(((frac >= lo) & (frac < hi) & ok).sum()))
            for lo, hi, lbl in frac_bands]
    return dict(
        total_sites=int(AD.shape[1]),
        never_covered=int((n_cov == 0).sum()),
        singleton_carriers=int((n_alt == 1).sum()),
        bands=rows, frac=frac, n_alt=n_alt, n_cov=n_cov,
    )


def shared_site_stats(AD, DP, idx, n_sub=500, seed=0):
    """Distribution of jointly covered informative sites per cell pair.

    This is the feasibility number. Clean separation in simulation needed a
    median of roughly 150+; around 40 accuracy degraded badly; below that the
    information simply is not present.
    """
    rng = np.random.default_rng(seed)
    rows = rng.choice(AD.shape[0], min(n_sub, AD.shape[0]), replace=False)
    _, M = binary_matrices(AD[rows], DP[rows], idx)
    shared = M @ M.T
    off = shared[~np.eye(len(rows), dtype=bool)]
    return dict(
        per_cell_median=float(np.median(M.sum(1))),
        shared_median=float(np.median(off)),
        shared_p10=float(np.percentile(off, 10)),
        frac_pairs_under_40=float(np.mean(off < 40)),
    )


def plot_carrier_fraction(AD, DP, K=None, min_covered=50, ax=None, out=None):
    """Histogram of carrier fraction, with 1/K gridlines.

    Interpretation: homozygous germline variants carried by j of K worms sit at
    j/K, so discrete peaks at multiples of 1/K confirm that the sites track
    worm identity and independently estimate K. A smooth featureless decay
    instead suggests recurrent artefacts rather than germline structure. Peak
    *heights* read out the pedigree: mass at 1/K means well-separated lines,
    mass at 2/K or 3/K means clades sharing recent ancestry.
    """
    import matplotlib.pyplot as plt
    frac, n_alt, n_cov = carrier_fraction(AD, DP)
    ok = (n_alt >= 2) & (n_cov >= min_covered)
    if ax is None:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].hist(frac[ok], bins=200, range=(0, 1), color="0.3")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("carrier fraction"); ax[0].set_ylabel("sites")
    ax[0].set_title("all sites with >=2 carriers (log y)")
    ax[1].hist(frac[ok & (frac > 0.02)], bins=120, range=(0, 1), color="0.3")
    if K:
        for j in range(1, K + 1):
            ax[1].axvline(j / K, color="crimson", lw=0.6, alpha=0.6)
        ax[1].set_title(f"vs multiples of 1/{K} (red)")
    ax[1].set_xlabel("carrier fraction")
    if out:
        ax[0].figure.savefig(out, dpi=130, bbox_inches="tight")
    return ax


def feasibility_report(AD, DP, K=None, frac_lo=0.03, frac_hi=0.97,
                       max_sites=20000, n_sub=500):
    """Print the numbers that decide whether to proceed. Returns them too."""
    s = site_summary(AD, DP)
    idx = select_sites(AD, DP, frac_lo, frac_hi, max_sites=max_sites)
    st = shared_site_stats(AD, DP, idx, n_sub=n_sub) if len(idx) else None

    print(f"total called sites        : {s['total_sites']}")
    print(f"never covered             : {s['never_covered']}")
    print(f"singleton carriers        : {s['singleton_carriers']}"
          "   <- cell-private somatic / error")
    for lbl, cnt in s["bands"]:
        print(f"  {lbl:<42} {cnt:>10}")
    print(f"sites selected for clustering: {len(idx)}")
    if st:
        print(f"informative sites per cell   : median {st['per_cell_median']:.0f}")
        print(f"shared per pair              : median {st['shared_median']:.0f} "
              f"| p10 {st['shared_p10']:.0f} "
              f"| pairs <40: {st['frac_pairs_under_40']:.2%}")
        m = st["shared_median"]
        verdict = ("comfortable" if m >= 150 else
                   "marginal" if m >= 40 else "not feasible")
        print(f"verdict: {verdict} (simulation: >=150 clean, ~40 degraded, <40 fails)")
    return dict(summary=s, sites=idx, shared=st)
