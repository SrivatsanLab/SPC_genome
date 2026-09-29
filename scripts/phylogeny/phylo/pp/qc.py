"""Composable cell-level QC predicates.

Each predicate is a callable ``(adata, cfg) -> np.ndarray[bool]`` returning
one bool per cell (True = pass). Register with :func:`register_predicate` and
list them in the config's ``qc.predicates`` array; :func:`apply_qc` ANDs them
and records the first-failing reason per cell.

The registry is the extension point: swap or add filters without touching
downstream code.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

Predicate = Callable[[ad.AnnData, dict[str, Any]], np.ndarray]
PREDICATES: dict[str, Predicate] = {}


def register_predicate(name: str) -> Callable[[Predicate], Predicate]:
    """Decorator: register a QC predicate under ``name``."""

    def _wrap(fn: Predicate) -> Predicate:
        PREDICATES[name] = fn
        return fn

    return _wrap


# --- built-in predicates ----------------------------------------------------


@register_predicate("not_background")
def _not_background(adata: ad.AnnData, cfg: dict[str, Any]) -> np.ndarray:
    """``capsule_type == 'cell'`` (drops 591 DNA-background capsules)."""
    return (adata.obs["capsule_type"].astype(str) == "cell").to_numpy()


@register_predicate("purity_gte")
def _purity_gte(adata: ad.AnnData, cfg: dict[str, Any]) -> np.ndarray:
    """Demultiplexing purity ≥ ``cfg['threshold']`` (default 0.7)."""
    t = float(cfg.get("threshold", 0.7))
    return adata.obs["purity"].to_numpy() >= t


@register_predicate("not_doublet")
def _not_doublet(adata: ad.AnnData, cfg: dict[str, Any]) -> np.ndarray:
    """``eff_modules <= cfg['max_eff_modules']`` (default 1.5).

    A capsule with mixed DNA has ``eff_modules ≈ 2`` and will float between
    clades in the tree.
    """
    t = float(cfg.get("max_eff_modules", 1.5))
    return adata.obs["eff_modules"].to_numpy() <= t


@register_predicate("worm_size_gte")
def _worm_size_gte(adata: ad.AnnData, cfg: dict[str, Any]) -> np.ndarray:
    """Cells belonging to worms with at least ``cfg['min_cells_per_worm']`` cells (after previous filters)."""
    t = int(cfg.get("min_cells_per_worm", 30))
    donors = adata.obs["donor"].astype(str).to_numpy()
    counts = pd.Series(donors).value_counts()
    ok_donors = set(counts[counts >= t].index)
    return np.array([d in ok_donors for d in donors])


@register_predicate("not_abu_module")
def _not_abu_module(adata: ad.AnnData, cfg: dict[str, Any]) -> np.ndarray:
    """Drop cells listed in ``cfg['csv']`` (from GEX; abu-module carriers).

    Optional: if ``cfg['csv']`` is None or missing, the predicate is a no-op
    (all True) with a warning printed by :func:`apply_qc`.
    """
    csv = cfg.get("csv")
    if not csv or not Path(csv).exists():
        return np.ones(adata.n_obs, dtype=bool)
    ids = set(pd.read_csv(csv)["cell_id"].astype(str))
    return np.array([c not in ids for c in adata.obs_names.astype(str)])


# --- application ------------------------------------------------------------


def apply_qc(
    adata: ad.AnnData,
    predicates: list[dict[str, Any]],
    *,
    inplace: bool = True,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Evaluate every predicate; AND the results; log per-cell exclusion reasons.

    Parameters
    ----------
    adata : AnnData
    predicates : list of dict
        Each dict has ``name`` (registry key) plus any predicate-specific args.
    inplace : bool
        Write ``adata.obs['qc_pass']`` and ``adata.obs['qc_reason']`` if True.

    Returns
    -------
    mask : np.ndarray, shape (n_obs,), bool
    reasons : DataFrame
        One row per cell: ``cell_id, qc_pass, qc_reason`` (first failing predicate).
    """
    n = adata.n_obs
    mask = np.ones(n, dtype=bool)
    reason = np.array(["pass"] * n, dtype=object)
    for spec in predicates:
        name = spec["name"]
        if name not in PREDICATES:
            raise KeyError(f"Unknown QC predicate '{name}'. Registered: {sorted(PREDICATES)}")
        sub = {k: v for k, v in spec.items() if k != "name"}
        if name == "not_abu_module" and not sub.get("csv"):
            print(f"[qc] warning: '{name}' predicate has no csv path — skipping.")
        this = PREDICATES[name](adata, sub)
        newly_failing = mask & ~this
        reason[newly_failing] = name
        mask &= this

    reasons = pd.DataFrame(
        {"cell_id": adata.obs_names.astype(str), "qc_pass": mask, "qc_reason": reason}
    )

    if inplace:
        adata.obs["qc_pass"] = mask
        adata.obs["qc_reason"] = pd.Categorical(reason)

    return mask, reasons
