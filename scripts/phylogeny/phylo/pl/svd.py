"""Diagnostic plots for SVD-based decomposition.

Moved from ``bin/haplotype_sweep.py`` verbatim (with the seaborn palette
inlined) so both the pool-level sweep and per-worm tree stages share one
implementation.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import TruncatedSVD

_COLORS = sns.color_palette("colorblind")
_N_PROBE = 40


def plot_scree(X: np.ndarray, out_dir: str | Path, title: str) -> None:
    """Singular values and consecutive gaps for a centered VAF matrix.

    The scree cliff sits at ``K - 1`` after per-variant centering
    (``HAPLOTYPE_ANALYSIS_SUMMARY.md`` §4). It flattens under unequal cell
    counts; read as suggestive.

    Writes ``scree.png`` and ``scree.csv`` under ``out_dir``.
    """
    out_dir = Path(out_dir)
    n = min(_N_PROBE, X.shape[1] - 1, X.shape[0] - 1)
    probe = TruncatedSVD(n_components=n, random_state=0).fit(X)
    sv = probe.singular_values_
    gaps = sv[:-1] - sv[1:]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(np.arange(1, n + 1), sv, "o-", color=_COLORS[2], linewidth=1.5)
    ax[0].set_xlabel("component", size=18)
    ax[0].set_ylabel("singular value", size=18)
    ax[0].tick_params(labelsize=18)
    ax[1].plot(np.arange(1, n), gaps, "o-", color=_COLORS[2], linewidth=2)
    ax[1].set_xlabel("component", size=18)
    ax[1].set_ylabel("gap to next", size=18)
    ax[1].tick_params(labelsize=18)
    sns.despine(fig)
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "scree.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        {
            "component": np.arange(1, n + 1),
            "singular_value": sv,
            "gap_to_next": np.append(gaps, np.nan),
            "cum_var": np.cumsum(probe.explained_variance_ratio_),
        }
    ).to_csv(out_dir / "scree.csv", index=False)


def plot_distances(
    d_all: np.ndarray,
    d_selected: np.ndarray,
    d_refit: np.ndarray,
    n_core: int,
    out_dir: str | Path,
    title: str,
) -> None:
    """Variant-to-centroid distances at three stages.

    - ``d_all``: every informative variant, provisional fit
    - ``d_selected``: same values restricted to the chosen core (truncated)
    - ``d_refit``: after refitting on the core (tighter by construction)

    Distance scale is not comparable across ``ncomp`` — unit-sphere distance
    grows with dimensionality.
    """
    out_dir = Path(out_dir)
    thresh = np.sort(d_all)[min(n_core, len(d_all)) - 1]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4))
    panels = [
        (d_all, "d_all: all informative variants", True),
        (d_selected, "d_selected: core only (truncated)", False),
        (d_refit, "d_refit: after refitting on the core", False),
    ]
    for a, (vals, ttl, mark) in zip(ax, panels):
        sns.histplot(vals, bins=200, alpha=0.75, color=_COLORS[2], linewidth=0, ax=a)
        if mark:
            a.axvline(thresh, color="crimson", ls="--", lw=1.5, label=f"core cut ({n_core} variants)")
            a.legend(fontsize=10)
        a.set_title(f"{ttl}\nmedian {np.median(vals):.3f}", fontsize=11)
        a.set_xlabel("distance to assigned module centroid", fontsize=12)
        a.set_ylabel("count of variants", fontsize=12)
        a.tick_params(labelsize=11)
    sns.despine(fig)
    fig.suptitle(title, y=1.03)
    fig.tight_layout()
    fig.savefig(out_dir / "distances.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
