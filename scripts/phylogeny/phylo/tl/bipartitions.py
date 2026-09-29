"""Enumerate candidate clades from per-worm ``tree.json`` files.

Every internal node in a tree defines *two* candidate clades — the L child
and the R child. Each row of the resulting :func:`extract_bipartitions`
DataFrame is one such candidate, annotated with the metadata a downstream
consumer needs for filtering: mutation count, within-clade VAF, sibling
stats, and a composite ``trustworthy`` flag matching the §6.1 amplitude
gate and the "MQuad actually selected a subset" rule.

Downstream uses:

- clade pseudobulking (join ``cell_ids`` against the RNA AnnData)
- gene-set enrichment (one row = one candidate clade to test)
- DNA-internal validation gates (VAF amplitude, split-half, jackknife —
  the ``trustworthy`` column is the amplitude gate applied on both sides)

See ``docs/phylo_pipeline.md`` §7 for the tree.json schema and downstream
recipes.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# --- amplitude gate defaults (§6.1 of T2_TREE_VALIDATION_PLAN.md) ----------
VAF_BAND_LO = 0.20
VAF_BAND_HI = 0.45
MIN_CLADE_SIZE_TRUSTWORTHY = 10

# --- continuous confidence defaults ---------------------------------------
# The five subcomponents each map an input to [0, 1] and are combined as a
# weighted arithmetic mean. Defaults are chosen so the binary
# `trustworthy` criterion corresponds roughly to composite ≥ 0.60.
CONFIDENCE_WEIGHTS = {
    "mquad": 1.0,       # per-variant informativeness (§10.2 item 4)
    "vaf": 1.0,         # own-side pooled VAF matches somatic-heterozygous target
    "sibling_vaf": 1.0, # sibling side also in-band (real splits are symmetric)
    "size": 1.0,        # both sides large enough to be statistically usable
    "purity": 1.0,      # cells clearly picked one module (§10.2 item 2)
}


def q_mquad(n_core: float, n_variants_used: float, saturation: int = 20) -> float:
    """MQuad selection quality.

    0 when the SVD builder fell back to full-panel k-means labels (n_core ==
    n_variants_used, meaning MQuad found too few informative variants to select).
    Otherwise n_core / (n_core + saturation) — saturates smoothly to 1 as the
    number of informative variants grows past ``saturation``.
    """
    if not np.isfinite(n_core) or not np.isfinite(n_variants_used):
        return 0.0
    if n_core >= n_variants_used:
        return 0.0
    return float(n_core / (n_core + saturation))


def q_vaf(vaf: float, target: float = 0.35, plateau: float = 0.10, edge: float = 0.10) -> float:
    """Amplitude quality: trapezoidal window centred on the §6.1 target.

    1.0 over ``[target-plateau, target+plateau]`` (default [0.25, 0.45]),
    linearly ramps to 0 over the ``edge`` band on either side (default width
    0.10 → zero outside [0.15, 0.55]).
    """
    if not np.isfinite(vaf):
        return 0.0
    dist = abs(vaf - target)
    if dist <= plateau:
        return 1.0
    if dist <= plateau + edge:
        return float(1.0 - (dist - plateau) / edge)
    return 0.0


def q_size(min_side_size: float, half_point: int = 10) -> float:
    """Clade-size quality: n / (n + half_point).

    Passes through 0.5 at ``half_point`` (matching the binary ≥ 10 threshold)
    and saturates toward 1 as clades grow.
    """
    if not np.isfinite(min_side_size) or min_side_size < 0:
        return 0.0
    return float(min_side_size / (min_side_size + half_point))


def q_purity(purity_median: float, floor: float = 0.5) -> float:
    """Module-assignment purity, rescaled.

    Purity at K=2 is bounded below by 0.5 (max of two weights that sum to 1).
    Maps [0.5, 1.0] → [0, 1]; NaN → 0.
    """
    if not np.isfinite(purity_median):
        return 0.0
    return float(max(0.0, (purity_median - floor) / (1.0 - floor)))


def compute_confidence(
    n_core: float,
    n_variants_used: float,
    pooled_vaf: float,
    sibling_vaf: float,
    min_side_size: float,
    purity_median: float,
    *,
    weights: dict[str, float] | None = None,
    mode: str = "mean",
) -> tuple[float, dict[str, float]]:
    """Composite clade confidence in [0, 1] and its subcomponents.

    Modes
    -----
    ``mean`` (default)
        Weighted arithmetic mean. Interpretable ("40% of the way to
        maximally confident"). Tolerant of a single weak subscore.
    ``geomean``
        Weighted geometric mean. "Weakest-link" semantics — one 0
        component zeros the whole score.
    ``min``
        Minimum subcomponent. Extreme weakest-link.
    """
    subs = {
        "mquad": q_mquad(n_core, n_variants_used),
        "vaf": q_vaf(pooled_vaf),
        "sibling_vaf": q_vaf(sibling_vaf),
        "size": q_size(min_side_size),
        "purity": q_purity(purity_median),
    }
    w = weights or CONFIDENCE_WEIGHTS
    keys = list(subs.keys())
    scores = np.array([subs[k] for k in keys])
    weights_arr = np.array([w.get(k, 1.0) for k in keys], dtype=float)
    weights_arr = weights_arr / weights_arr.sum()
    if mode == "mean":
        conf = float(np.dot(scores, weights_arr))
    elif mode == "geomean":
        # Clip to avoid log(0); when a score is 0 the geomean is 0.
        if (scores == 0).any():
            conf = 0.0
        else:
            conf = float(np.exp(np.dot(weights_arr, np.log(scores))))
    elif mode == "min":
        conf = float(scores.min())
    else:
        raise ValueError(f"unknown confidence mode {mode!r}")
    return conf, subs


def _walk(node: dict) -> Iterable[dict]:
    yield node
    if node.get("split") is not None:
        yield from _walk(node["left"])
        yield from _walk(node["right"])


def _bipartitions_for_tree(
    tree: dict,
    *,
    vaf_band_lo: float,
    vaf_band_hi: float,
    min_clade_size: int,
) -> list[dict]:
    worm = tree["worm"]
    rows: list[dict] = []
    for node in _walk(tree["root"]):
        s = node.get("split")
        if s is None:
            continue
        mquad_selected = bool(s["n_core"] < s["n_variants_used"])
        vaf = {"L": s["module_pooled_vaf"]["L"], "R": s["module_pooled_vaf"]["R"]}
        n_cells = {"L": s["module_n_cells"]["L"], "R": s["module_n_cells"]["R"]}
        purity_med = {
            "L": s["module_purity_median"]["L"],
            "R": s["module_purity_median"]["R"],
        }
        n_mut = {
            "L": len(s["module_variants"]["L"]),
            "R": len(s["module_variants"]["R"]),
        }
        for side in ("L", "R"):
            other = "R" if side == "L" else "L"
            child = node["left"] if side == "L" else node["right"]
            in_band = vaf_band_lo <= vaf[side] <= vaf_band_hi
            sibling_in_band = vaf_band_lo <= vaf[other] <= vaf_band_hi
            min_size = min(n_cells[side], n_cells[other])
            trustworthy = bool(
                mquad_selected
                and in_band
                and sibling_in_band
                and min_size >= min_clade_size
            )
            confidence, subs = compute_confidence(
                n_core=s["n_core"],
                n_variants_used=s["n_variants_used"],
                pooled_vaf=vaf[side],
                sibling_vaf=vaf[other],
                min_side_size=min_size,
                purity_median=purity_med[side],
            )
            rows.append(
                {
                    "worm": worm,
                    "node_path": node["path"] or "root",
                    "side": side,
                    "clade_path": child["path"],
                    "depth": child["depth"],
                    "n_cells": n_cells[side],
                    "n_mutations": n_mut[side],
                    "pooled_vaf": vaf[side],
                    "module_purity_median": purity_med[side],
                    "sibling_n_cells": n_cells[other],
                    "sibling_vaf": vaf[other],
                    "n_core": s["n_core"],
                    "n_variants_used": s["n_variants_used"],
                    "mquad_selected": mquad_selected,
                    "in_vaf_band": in_band,
                    "sibling_in_vaf_band": sibling_in_band,
                    "min_side_size": min_size,
                    "trustworthy": trustworthy,
                    "confidence": confidence,
                    "conf_mquad": subs["mquad"],
                    "conf_vaf": subs["vaf"],
                    "conf_sibling_vaf": subs["sibling_vaf"],
                    "conf_size": subs["size"],
                    "conf_purity": subs["purity"],
                    "cell_ids": ",".join(child["cells"]),  # comma-joined for CSV
                }
            )
    return rows


def extract_bipartitions(
    run_dir: str | Path,
    *,
    vaf_band_lo: float = VAF_BAND_LO,
    vaf_band_hi: float = VAF_BAND_HI,
    min_clade_size: int = MIN_CLADE_SIZE_TRUSTWORTHY,
    write: bool = True,
) -> pd.DataFrame:
    """Walk every ``tree.json`` under ``run_dir/trees/`` and return one row per clade.

    Parameters
    ----------
    run_dir : path
        Path to a phylo run directory (contains ``trees/worm_*/tree.json``).
    vaf_band_lo, vaf_band_hi : float
        Amplitude gate for the ``in_vaf_band`` / ``sibling_in_vaf_band`` /
        ``trustworthy`` columns. Defaults to §6.1's [0.20, 0.45].
    min_clade_size : int
        Both sides must have at least this many cells for ``trustworthy=True``.
    write : bool
        Write ``trees/bipartitions.csv`` inside ``run_dir``.

    Returns
    -------
    DataFrame
        See ``docs/phylo_pipeline.md`` for the column list. Two rows per
        internal node (one per side).
    """
    run_dir = Path(run_dir)
    trees_dir = run_dir / "trees"
    all_rows: list[dict] = []
    for tj in sorted(trees_dir.glob("worm_*/tree.json")):
        tree = json.load(open(tj))
        all_rows.extend(
            _bipartitions_for_tree(
                tree,
                vaf_band_lo=vaf_band_lo,
                vaf_band_hi=vaf_band_hi,
                min_clade_size=min_clade_size,
            )
        )
    df = pd.DataFrame(all_rows)
    if write:
        out = trees_dir / "bipartitions.csv"
        df.to_csv(out, index=False)
        print(f"[bipartitions] wrote {len(df)} rows → {out}")
        n_ok = int(df["trustworthy"].sum())
        n_worms_ok = int(df.loc[df["trustworthy"], "worm"].nunique())
        print(f"[bipartitions] {n_ok} trustworthy clade-sides across {n_worms_ok} worms")
    return df
