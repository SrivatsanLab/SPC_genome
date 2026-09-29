"""MQuad per-variant informativeness (ΔBIC of 1- vs 2-component binomial mixture).

Wraps :class:`mquad.mquad.Mquad` for use inside the recursive tree builder.
Returns a per-variant ΔBIC array (positive → 2-component preferred → informative)
so callers can pick their own threshold or top-N.

Notes
-----
- MQuad's default ``minDP=10`` is inappropriate for 1× per-cell coverage; we default
  to ``minDP=1`` and rely on aggregate qualified-cell counts per variant.
- MQuad writes CSVs to disk when ``export_csv=True``; we disable that.
- MQuad's internal multiprocessing spawns a pool on every call. For per-node use in
  the recursion driver, ``nproc=1`` avoids pool-creation overhead; set higher only
  for whole-panel one-shot calls.

Reference: Kwok, Ho, Lin, Yeung, Huang, McCarthy & Sham, "MQuad enables clonal
substructure discovery using single cell mitochondrial variants",
Nat Commun 13:1205 (2022), doi:10.1038/s41467-022-28845-0.
"""

from __future__ import annotations

import contextlib
import io
import os
import tempfile

import numpy as np
from scipy.sparse import csr_matrix, issparse


def compute_delta_bic(
    AD: np.ndarray,
    DP: np.ndarray,
    *,
    min_dp: int = 1,
    min_ad: int = 1,
    nproc: int = 1,
    quiet: bool = True,
    seed: int | None = 0,
) -> np.ndarray:
    """Per-variant ΔBIC = BIC(1-component) − BIC(2-component) binomial mixture.

    Parameters
    ----------
    AD, DP : ndarray, shape (n_cells, n_variants)
        Alt-read counts and total depth. Sparse or dense.
    min_dp : int
        Minimum per-cell depth to include a cell in the fit for a given variant.
    min_ad : int
        Minimum per-cell alt-read count to include a cell.
    nproc : int
        MQuad's internal pool size. Use 1 for per-node calls, higher for one-shot.
        **Only nproc=1 gives fully deterministic results** — under multiprocessing
        BBMix's random init in worker processes is not reproducible via a global
        seed. This wrapper forces ``nproc=1`` when ``seed`` is set.
    quiet : bool
        Suppress MQuad's stdout chatter (one line per variant).
    seed : int or None
        Global numpy seed applied before the fit. Set to None to opt out.
        BBMix's EM init uses ``np.random`` internally without accepting a seed
        argument; seeding globally makes single-process runs reproducible.

    Returns
    -------
    delta_bic : ndarray, shape (n_variants,), float
        Positive → 2-component preferred → variant is informative.
        NaN for variants MQuad could not fit (all cells filtered out).
    """
    from mquad.mquad import Mquad

    if seed is not None:
        np.random.seed(int(seed))
        if nproc != 1:
            # Determinism only holds in-process; force single-threaded.
            nproc = 1

    if not issparse(AD):
        AD = csr_matrix(AD.astype(np.int64))
    if not issparse(DP):
        DP = csr_matrix(DP.astype(np.int64))
    # MQuad expects (n_variants, n_cells)
    AD_T = AD.T.astype(np.int64)
    DP_T = DP.T.astype(np.int64)

    with tempfile.TemporaryDirectory() as tmp:
        stream = io.StringIO() if quiet else None
        ctx = contextlib.redirect_stdout(stream) if quiet else contextlib.nullcontext()
        with ctx:
            m = Mquad(AD=AD_T, DP=DP_T)
            m.fit_deltaBIC(out_dir=tmp, nproc=nproc, minDP=min_dp, minAD=min_ad, export_csv=False)
    # MQuad returns deltaBIC as object dtype; coerce.
    return m.df["deltaBIC"].to_numpy(dtype=float)


def select_informative(
    AD: np.ndarray,
    DP: np.ndarray,
    *,
    delta_bic_threshold: float = 10.0,
    min_dp: int = 1,
    min_ad: int = 1,
    nproc: int = 1,
    seed: int | None = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Boolean mask over variants keeping those with ΔBIC ≥ threshold.

    Returns
    -------
    keep : ndarray, shape (n_variants,), bool
    delta_bic : ndarray, shape (n_variants,), float
        The raw scores (retained for diagnostics; NaN kept as NaN and excluded from ``keep``).
    """
    delta_bic = compute_delta_bic(
        AD, DP, min_dp=min_dp, min_ad=min_ad, nproc=nproc, quiet=True, seed=seed
    )
    keep = np.isfinite(delta_bic) & (delta_bic >= delta_bic_threshold)
    return keep, delta_bic
