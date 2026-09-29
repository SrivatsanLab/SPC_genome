"""Depth-matched cross-worm clade-expression matching.

Question addressed: *does the per-worm bipartition tree carry RNA-coherent
signal that generalises across worms?* If yes, depth-matched clades across
worms should have higher pseudobulk-expression similarity than expected
under random bipartitions.

Pipeline stages:

1. **Pseudobulks** (:func:`build_pseudobulks`) — sum raw counts across the
   cells of each clade (both sides of every internal node); attach
   ``log_cpm`` layer, one clade per row.
2. **Pairwise similarity** (:func:`pairwise_similarity`) — long-form
   clade × clade table for Pearson / Spearman / cosine on the log_cpm
   matrix.
3. **Orientation matching** (:func:`orientation_match_at_depth`) — at
   each depth ≥ 1, for every worm pair, exhaustively try all
   bipartite matchings between their clades at that depth (Hungarian
   assignment for degree ≥ 2) and take the matching that maximises
   summed similarity. Records the per-pair best-orientation similarity.
4. **Null distribution** (:func:`null_matched_similarity`) — 500 random
   bipartitions per worm matched on cell count. For each draw, stages
   1–3 are rerun to build a null distribution of matched similarities.
5. **Test + report** (:func:`summarize`) — per-depth: observed matched
   similarity vs the null; empirical p-values.
6. **Cross-worm groups** (:func:`emit_cross_worm_groups`) — matched
   clades → group assignments for downstream pooled tissue enrichment.

Consume via ``python -m phylo clade_matching --config <cfg>`` or by calling
:func:`run` from a notebook.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats
from scipy.optimize import linear_sum_assignment


# ---------------------------------------------------------------------------
# 1. Pseudobulk construction
# ---------------------------------------------------------------------------


def _walk_tree(node: dict):
    yield node
    if node.get("split") is not None:
        yield from _walk_tree(node["left"])
        yield from _walk_tree(node["right"])


def _collect_clades_from_trees(trees_dir: Path) -> list[dict]:
    """One row per clade (each side of every internal node)."""
    out: list[dict] = []
    for tj in sorted(trees_dir.glob("worm_*/tree.json")):
        tree = json.load(open(tj))
        worm = tree["worm"]
        for node in _walk_tree(tree["root"]):
            s = node.get("split")
            if s is None:
                continue
            for side in ("L", "R"):
                child = node["left"] if side == "L" else node["right"]
                other = "R" if side == "L" else "L"
                out.append(
                    {
                        "worm": worm,
                        "node_path": node["path"] or "root",
                        "side": side,
                        "clade_path": child["path"],
                        "depth": child["depth"],
                        "n_cells_tree": child["n_cells"],
                        "cells": child["cells"],
                        "n_mutations": len(s["module_variants"][side]),
                        "pooled_vaf": s["module_pooled_vaf"][side],
                        "module_purity_median": s["module_purity_median"][side],
                        "n_core": s["n_core"],
                        "n_variants_used": s["n_variants_used"],
                        "mquad_selected": bool(s["n_core"] < s["n_variants_used"]),
                        "sibling_vaf": s["module_pooled_vaf"][other],
                        "sibling_n_cells": s["module_n_cells"][other],
                    }
                )
    return out


def _clade_id(worm: str, clade_path: str) -> str:
    return f"{worm}.{clade_path}"


def _log_cpm(X: np.ndarray) -> np.ndarray:
    """log((counts / total) * 1e6 + 1) per row."""
    row_sums = X.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1)
    return np.log1p((X / row_sums) * 1e6)


def build_pseudobulks(
    run_dir: str | Path,
    rna_path: str | Path,
    *,
    counts_layer: str = "counts",
    write: bool = True,
) -> ad.AnnData:
    """Per-clade pseudobulk counts + log-CPM.

    Parameters
    ----------
    run_dir : path
        Phylo run directory (contains ``trees/worm_*/tree.json``).
    rna_path : path
        RNA AnnData with raw integer counts in ``layers[counts_layer]``.
    counts_layer : str
        Layer name to sum. Defaults to ``"counts"``.
    write : bool
        Write ``clade_matching/pseudobulks.h5ad`` under ``run_dir``.

    Returns
    -------
    AnnData
        ``obs`` = one row per clade with worm/node_path/side/depth/n_cells_tree/
        n_cells_matched/n_mutations/pooled_vaf/trustworthy. ``X`` = summed raw
        counts. ``layers['log_cpm']`` = log(CPM+1).
    """
    run_dir = Path(run_dir)
    trees_dir = run_dir / "trees"
    clades = _collect_clades_from_trees(trees_dir)
    print(f"[pseudobulks] {len(clades)} clades from {len(list(trees_dir.glob('worm_*/tree.json')))} trees")

    rna = ad.read_h5ad(rna_path)
    if counts_layer not in rna.layers:
        raise KeyError(f"'{counts_layer}' not in rna.layers; have {list(rna.layers)}")
    counts = rna.layers[counts_layer]
    if sp.issparse(counts):
        counts = counts.toarray()
    rna_index: dict[str, int] = {c: i for i, c in enumerate(rna.obs_names.astype(str))}

    n_clades = len(clades)
    n_genes = rna.n_vars
    X = np.zeros((n_clades, n_genes), dtype=np.float64)
    obs_rows = []
    for i, cl in enumerate(clades):
        rows = [rna_index[c] for c in cl["cells"] if c in rna_index]
        X[i, :] = counts[rows, :].sum(axis=0) if rows else 0
        obs_rows.append(
            {
                **{k: v for k, v in cl.items() if k != "cells"},
                "clade_id": _clade_id(cl["worm"], cl["clade_path"]),
                "n_cells_matched": len(rows),
                "total_counts": float(X[i, :].sum()),
            }
        )
    obs = pd.DataFrame(obs_rows).set_index("clade_id")

    # trustworthy composite (matches phylo.tl.bipartitions defaults).
    obs["in_vaf_band"] = (obs["pooled_vaf"] >= 0.20) & (obs["pooled_vaf"] <= 0.45)
    obs["sibling_in_vaf_band"] = (obs["sibling_vaf"] >= 0.20) & (obs["sibling_vaf"] <= 0.45)
    obs["min_side_size"] = np.minimum(obs["n_cells_tree"], obs["sibling_n_cells"])
    obs["trustworthy"] = (
        obs["mquad_selected"]
        & obs["in_vaf_band"]
        & obs["sibling_in_vaf_band"]
        & (obs["min_side_size"] >= 10)
    )
    # Continuous confidence + subcomponents.
    from .bipartitions import compute_confidence
    confs = obs.apply(
        lambda r: compute_confidence(
            n_core=r["n_core"], n_variants_used=r["n_variants_used"],
            pooled_vaf=r["pooled_vaf"], sibling_vaf=r["sibling_vaf"],
            min_side_size=r["min_side_size"], purity_median=r["module_purity_median"],
        ),
        axis=1,
    )
    obs["confidence"] = confs.map(lambda t: t[0])
    for k in ("mquad", "vaf", "sibling_vaf", "size", "purity"):
        obs[f"conf_{k}"] = confs.map(lambda t, k=k: t[1][k])

    log_cpm = _log_cpm(X)
    pseudobulks = ad.AnnData(X=X.astype(np.float32), obs=obs, var=rna.var.copy())
    pseudobulks.layers["log_cpm"] = log_cpm.astype(np.float32)

    if write:
        out_dir = run_dir / "clade_matching"
        out_dir.mkdir(parents=True, exist_ok=True)
        pseudobulks.write_h5ad(out_dir / "pseudobulks.h5ad")
        print(f"[pseudobulks] wrote {n_clades} clades × {n_genes} genes → {out_dir/'pseudobulks.h5ad'}")

    return pseudobulks


# ---------------------------------------------------------------------------
# 2. Pairwise similarity
# ---------------------------------------------------------------------------


def _rowwise_similarity(M: np.ndarray, method: str) -> np.ndarray:
    """M shape (n_clades, n_genes) → (n_clades, n_clades) similarity matrix."""
    method = method.lower()
    if method == "pearson":
        # Row-center and unit-scale; then M @ M.T / n_genes gives Pearson.
        Mc = M - M.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(Mc, axis=1, keepdims=True)
        norm = np.where(norm > 0, norm, 1)
        Mn = Mc / norm
        return Mn @ Mn.T
    if method == "cosine":
        norm = np.linalg.norm(M, axis=1, keepdims=True)
        norm = np.where(norm > 0, norm, 1)
        Mn = M / norm
        return Mn @ Mn.T
    if method == "spearman":
        # Rank each row, then Pearson on ranks.
        R = stats.rankdata(M, axis=1).astype(np.float64)
        Rc = R - R.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(Rc, axis=1, keepdims=True)
        norm = np.where(norm > 0, norm, 1)
        return (Rc / norm) @ (Rc / norm).T
    raise ValueError(f"unknown similarity method {method!r}")


def _select_genes(pseudobulks: ad.AnnData, min_total_counts: int = 20) -> np.ndarray:
    """Boolean mask over genes with enough total pseudobulk expression."""
    total = pseudobulks.X.sum(axis=0)
    return total >= min_total_counts


def _residualize_by_worm(M: np.ndarray, worms: np.ndarray) -> np.ndarray:
    """Subtract per-worm mean expression from each clade's expression vector.

    Removes the shared worm-level baseline (dissection batch, RNA quality,
    ambient composition) that dominates raw log-CPM similarity — otherwise
    any two clades from the same worm look nearly identical regardless of
    biology.
    """
    out = M.copy()
    for w in np.unique(worms):
        idx = worms == w
        out[idx] = out[idx] - out[idx].mean(axis=0, keepdims=True)
    return out


def pairwise_similarity(
    pseudobulks: ad.AnnData,
    method: str = "pearson",
    *,
    min_total_counts: int = 20,
    residualize_by_worm: bool = True,
) -> pd.DataFrame:
    """Long-form clade × clade similarity table.

    Parameters
    ----------
    method : str
        ``pearson``, ``spearman``, or ``cosine`` on ``layers['log_cpm']``.
    min_total_counts : int
        Drop genes with total pseudobulk expression below this threshold.
    residualize_by_worm : bool
        If True (default), subtract each worm's mean per-gene expression from
        its clades before computing similarity. Necessary at 1× coverage
        because shared worm-level signal dominates raw log-CPM; without this
        step observed ≈ null median at every depth.

    Returns columns: ``clade_A``, ``clade_B``, ``sim`` (float).
    """
    gene_mask = _select_genes(pseudobulks, min_total_counts=min_total_counts)
    M = pseudobulks.layers["log_cpm"][:, gene_mask].astype(np.float64)
    if residualize_by_worm:
        M = _residualize_by_worm(M, pseudobulks.obs["worm"].astype(str).to_numpy())
    sim = _rowwise_similarity(M, method=method)
    ids = pseudobulks.obs_names.astype(str).to_numpy()
    n = len(ids)
    iu, ju = np.triu_indices(n, k=1)
    return pd.DataFrame(
        {
            "clade_A": ids[iu],
            "clade_B": ids[ju],
            "sim": sim[iu, ju].astype(np.float32),
        }
    )


# ---------------------------------------------------------------------------
# 3. Orientation matching at a depth
# ---------------------------------------------------------------------------


def _pair_sim(sim_lookup: dict, a: str, b: str) -> float:
    key = (a, b) if a < b else (b, a)
    return sim_lookup.get(key, np.nan)


def orientation_match_at_depth(
    pseudobulks: ad.AnnData,
    sim_long: pd.DataFrame,
    depth: int,
) -> pd.DataFrame:
    """Best-orientation matching per worm pair for all clades at ``depth``.

    For each pair (A, B) of worms with ≥ 1 clade each at this depth, tries
    every bipartite matching (Hungarian for degrees ≥ 2) and keeps the
    matching that maximises summed similarity. Reports:

    - ``worm_A, worm_B``
    - ``matched_pairs`` — list of ``(clade_A_side, clade_B_side)``
    - ``matched_sim_sum`` — summed similarity of the chosen matching
    - ``matched_sim_mean`` — mean (per matched pair)
    - ``max_over_alternative`` — gap vs the best *alternative* orientation
    - ``n_a, n_b`` — number of clades at this depth in each worm
    """
    obs = pseudobulks.obs
    cur = obs[obs["depth"] == depth]
    if cur.empty:
        return pd.DataFrame()

    sim_lookup = {
        (min(a, b), max(a, b)): s
        for a, b, s in zip(sim_long["clade_A"], sim_long["clade_B"], sim_long["sim"])
    }

    worm_to_ids = cur.groupby("worm", observed=True).apply(lambda g: g.index.tolist(), include_groups=False)
    worm_to_ids = worm_to_ids[worm_to_ids.map(len) > 0]
    worms = worm_to_ids.index.tolist()
    rows = []
    for i in range(len(worms)):
        for j in range(i + 1, len(worms)):
            A, B = worms[i], worms[j]
            ids_a = worm_to_ids[A]
            ids_b = worm_to_ids[B]
            na, nb = len(ids_a), len(ids_b)
            if na == 0 or nb == 0:
                continue
            # cost = -sim (linear_sum_assignment minimises); NaN → 0
            cost = np.zeros((na, nb))
            sim_mat = np.zeros((na, nb))
            for ia, ca in enumerate(ids_a):
                for ib, cb in enumerate(ids_b):
                    s = _pair_sim(sim_lookup, ca, cb)
                    sim_mat[ia, ib] = 0.0 if np.isnan(s) else s
                    cost[ia, ib] = -sim_mat[ia, ib]
            row_ind, col_ind = linear_sum_assignment(cost)
            matched_sim = sim_mat[row_ind, col_ind]
            # alternative: consider only *other* orderings for the "gap"
            # (only meaningful at na=nb=2 where there are 2 orderings).
            if na == nb == 2:
                alt_sim = sim_mat[[0, 1], [1, 0]] if list(col_ind) == [0, 1] else sim_mat[[0, 1], [0, 1]]
                alt_sum = float(alt_sim.sum())
            else:
                alt_sum = np.nan
            rows.append(
                {
                    "depth": depth,
                    "worm_A": A,
                    "worm_B": B,
                    "n_a": na,
                    "n_b": nb,
                    "matched_pairs": [(ids_a[a], ids_b[b]) for a, b in zip(row_ind, col_ind)],
                    "matched_sim_sum": float(matched_sim.sum()),
                    "matched_sim_mean": float(matched_sim.mean()),
                    "max_over_alternative": float(matched_sim.sum() - alt_sum) if not np.isnan(alt_sum) else np.nan,
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. Null distribution
# ---------------------------------------------------------------------------


def _random_bipartition_pseudobulks(
    rna_counts: np.ndarray,
    rna_index: dict,
    clades: list[dict],
    rng: np.random.Generator,
) -> ad.AnnData:
    """For each observed clade, replace its cells with a random subset of its
    worm's cells matched on cell count. Return an AnnData of the same shape
    as the observed pseudobulks.
    """
    # Group observed clades by worm and identify the cell pool per worm.
    worm_cells: dict[str, list[str]] = {}
    for cl in clades:
        worm_cells.setdefault(cl["worm"], []).extend(cl["cells"])
    # Deduplicate per worm.
    for w, cs in worm_cells.items():
        worm_cells[w] = list(dict.fromkeys(cs))  # order-preserving unique

    n_clades = len(clades)
    n_genes = rna_counts.shape[1]
    X = np.zeros((n_clades, n_genes), dtype=np.float64)
    obs_rows = []
    for i, cl in enumerate(clades):
        pool = worm_cells[cl["worm"]]
        pool_idx = np.array([rna_index[c] for c in pool if c in rna_index])
        n_want = min(cl["n_cells_tree"], pool_idx.size)
        pick = rng.choice(pool_idx, size=n_want, replace=False) if n_want else np.array([], dtype=int)
        X[i, :] = rna_counts[pick, :].sum(axis=0) if pick.size else 0
        obs_rows.append({
            "worm": cl["worm"],
            "node_path": cl["node_path"],
            "side": cl["side"],
            "depth": cl["depth"],
            "n_cells_tree": cl["n_cells_tree"],
            "n_cells_matched": int(pick.size),
            "clade_id": _clade_id(cl["worm"], cl["clade_path"]),
        })
    obs = pd.DataFrame(obs_rows).set_index("clade_id")
    pb = ad.AnnData(X=X.astype(np.float32), obs=obs)
    pb.layers["log_cpm"] = _log_cpm(X).astype(np.float32)
    return pb


def null_matched_similarity(
    run_dir: str | Path,
    rna_path: str | Path,
    *,
    method: str = "pearson",
    n_null: int = 500,
    counts_layer: str = "counts",
    min_total_counts: int = 20,
    residualize_by_worm: bool = True,
    trustworthy_only: bool = False,
    seed: int = 0,
) -> pd.DataFrame:
    """Null distribution of orientation-matched cross-worm similarity per depth.

    For each of ``n_null`` random draws, replace each observed clade's cells
    with a random size-matched subset of that worm's cells; recompute
    pseudobulks + log-CPM + pairwise similarity + orientation matching;
    aggregate the median matched similarity per depth.

    Returns a long DataFrame with columns
    ``(draw, depth, worm_A, worm_B, matched_sim_mean)``.
    """
    run_dir = Path(run_dir)
    clades = _collect_clades_from_trees(run_dir / "trees")
    if trustworthy_only:
        # Recompute trustworthy per clade (matches build_pseudobulks logic).
        def _trust(cl):
            in_band = 0.20 <= cl["pooled_vaf"] <= 0.45
            sib_in = 0.20 <= cl["sibling_vaf"] <= 0.45
            min_size = min(cl["n_cells_tree"], cl["sibling_n_cells"])
            return cl["mquad_selected"] and in_band and sib_in and min_size >= 10
        clades = [cl for cl in clades if _trust(cl)]
        print(f"[null] trustworthy-only: kept {len(clades)} clades")
    rna = ad.read_h5ad(rna_path)
    counts = rna.layers[counts_layer]
    if sp.issparse(counts):
        counts = counts.toarray()
    rna_index = {c: i for i, c in enumerate(rna.obs_names.astype(str))}

    depths = sorted({cl["depth"] for cl in clades})

    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for draw in range(n_null):
        pb = _random_bipartition_pseudobulks(counts, rna_index, clades, rng)
        gene_mask = _select_genes(pb, min_total_counts=min_total_counts)
        M = pb.layers["log_cpm"][:, gene_mask].astype(np.float64)
        if residualize_by_worm:
            M = _residualize_by_worm(M, pb.obs["worm"].astype(str).to_numpy())
        S = _rowwise_similarity(M, method=method)
        ids = pb.obs_names.astype(str).to_numpy()
        iu, ju = np.triu_indices(len(ids), k=1)
        sim_long = pd.DataFrame({"clade_A": ids[iu], "clade_B": ids[ju], "sim": S[iu, ju]})
        for d in depths:
            m = orientation_match_at_depth(pb, sim_long, d)
            if m.empty:
                continue
            for _, r in m.iterrows():
                rows.append({
                    "draw": draw,
                    "depth": d,
                    "worm_A": r["worm_A"],
                    "worm_B": r["worm_B"],
                    "matched_sim_mean": r["matched_sim_mean"],
                })
        if (draw + 1) % 50 == 0:
            print(f"[null] {method}: {draw + 1}/{n_null} draws")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. Test + report
# ---------------------------------------------------------------------------


def summarize(
    observed: pd.DataFrame,
    null: pd.DataFrame,
    out_dir: Path,
    method: str,
) -> pd.DataFrame:
    """Per-depth test result: observed median matched sim vs null median.

    Empirical p-value: fraction of null draws whose median matched sim (over
    worm pairs at that depth) meets or exceeds the observed median.
    """
    rows = []
    for depth in sorted(observed["depth"].unique()):
        obs = observed[observed["depth"] == depth]["matched_sim_mean"]
        nul = null[null["depth"] == depth]
        obs_med = float(obs.median())
        # median matched sim over worm pairs, per draw
        null_meds = nul.groupby("draw")["matched_sim_mean"].median()
        p = float((null_meds >= obs_med).mean())
        rows.append({
            "method": method,
            "depth": depth,
            "n_worm_pairs": len(obs),
            "n_worms": int(pd.concat([observed[observed.depth == depth]["worm_A"],
                                        observed[observed.depth == depth]["worm_B"]]).nunique()),
            "obs_median_matched_sim": obs_med,
            "null_median": float(null_meds.median()),
            "null_p95": float(null_meds.quantile(0.95)),
            "obs_gt_null_p95": obs_med > float(null_meds.quantile(0.95)),
            "empirical_p": p,
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / f"test_results_{method}.csv", index=False)
    print(f"[summary/{method}]"); print(df.to_string(index=False))
    return df


# ---------------------------------------------------------------------------
# 6. Cross-worm groups
# ---------------------------------------------------------------------------


def emit_cross_worm_groups(
    pseudobulks: ad.AnnData,
    sim_long: pd.DataFrame,
    out_dir: Path,
    method: str,
) -> pd.DataFrame:
    """Cluster clades at each depth into K groups by expression similarity.

    K is chosen per depth as the median number of clades per worm at that
    depth, clamped to [2, 2**depth]. Spectral clustering on the similarity
    matrix; group confidence = mean within-group pairwise similarity.

    At depth 1: K=2 → two groups corresponding (if the tree captures real
    lineage signal) to the two sides of the ancestral split.
    """
    from sklearn.cluster import SpectralClustering

    lookup: dict[tuple[str, str], float] = {
        (min(a, b), max(a, b)): s
        for a, b, s in zip(sim_long["clade_A"], sim_long["clade_B"], sim_long["sim"])
    }

    all_rows: list[dict] = []
    for depth in sorted(pseudobulks.obs["depth"].unique()):
        cur = pseudobulks.obs[pseudobulks.obs["depth"] == depth]
        if cur["worm"].nunique() < 2:
            continue
        clades_per_worm = cur.groupby("worm", observed=True).size()
        K_target = int(clades_per_worm.median())
        K = max(2, min(K_target, 2 ** int(depth)))
        clades_at_d = cur.index.astype(str).tolist()
        n = len(clades_at_d)
        if n < K:
            continue

        # Build similarity matrix, shift to non-negative for spectral affinity.
        S = np.zeros((n, n), dtype=np.float64)
        for i, ci in enumerate(clades_at_d):
            for j, cj in enumerate(clades_at_d):
                if i == j:
                    S[i, j] = 1.0
                    continue
                key = (min(ci, cj), max(ci, cj))
                S[i, j] = lookup.get(key, 0.0)
        # Shift so min becomes 0; add tiny epsilon to avoid degenerate all-zero.
        A = S - S.min() + 1e-6

        try:
            sc = SpectralClustering(
                n_clusters=K, affinity="precomputed", random_state=0, assign_labels="discretize"
            )
            labels = sc.fit_predict(A)
        except Exception as e:
            print(f"[groups/{method}] depth={depth}: spectral clustering failed ({e}); skipping.")
            continue

        for k in range(K):
            mask = labels == k
            group_ids = [clades_at_d[i] for i in range(n) if mask[i]]
            if len(group_ids) < 2:
                within_conf = np.nan
            else:
                sub = S[np.ix_(mask, mask)]
                iu = np.triu_indices(mask.sum(), k=1)
                within_conf = float(sub[iu].mean())
            for clade_id in group_ids:
                row = pseudobulks.obs.loc[clade_id]
                all_rows.append({
                    "depth": int(depth),
                    "K": int(K),
                    "group_id": f"d{int(depth)}_g{k:02d}",
                    "group_size": int(mask.sum()),
                    "group_confidence_mean_sim": within_conf,
                    "clade_id": clade_id,
                    "worm": row["worm"],
                    "node_path": row["node_path"],
                    "side": row["side"],
                    "n_cells": int(row["n_cells_matched"]),
                    "n_mutations": int(row["n_mutations"]),
                    "pooled_vaf": float(row["pooled_vaf"]),
                    "trustworthy": bool(row["trustworthy"]),
                })

    df = pd.DataFrame(all_rows)
    df.to_csv(out_dir / f"cross_worm_groups_{method}.csv", index=False)
    print(f"[groups/{method}] wrote {len(df)} clade rows → cross_worm_groups_{method}.csv")
    return df


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run(
    run_dir: str | Path,
    rna_path: str | Path,
    *,
    similarity_methods: tuple[str, ...] = ("pearson", "spearman", "cosine"),
    n_null: int = 500,
    min_total_counts: int = 20,
    counts_layer: str = "counts",
    min_worms_per_depth: int = 2,
    residualize_by_worm: bool = True,
    trustworthy_only: bool = False,
    output_subdir: str | None = None,
    seed: int = 0,
) -> dict[str, Any]:
    """End-to-end: build pseudobulks, run matching + null per method, emit outputs.

    Writes under ``run_dir/clade_matching/``:

    - ``pseudobulks.h5ad``
    - ``pairwise_similarity_<method>.csv`` for each method
    - ``matches_depth<d>_<method>.csv``, ``null_<method>.csv``,
      ``test_results_<method>.csv``
    - ``cross_worm_groups_<method>.csv``
    """
    run_dir = Path(run_dir)
    subdir = output_subdir or ("clade_matching_trustworthy" if trustworthy_only else "clade_matching")
    out_dir = run_dir / subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    pb = build_pseudobulks(run_dir, rna_path, counts_layer=counts_layer, write=False)
    if trustworthy_only:
        n_before = pb.n_obs
        pb = pb[pb.obs["trustworthy"]].copy()
        print(f"[run] trustworthy-only filter: {pb.n_obs}/{n_before} clades kept")
    pb.write_h5ad(out_dir / "pseudobulks.h5ad")
    print(f"[pseudobulks] wrote → {out_dir/'pseudobulks.h5ad'}")

    # Depths that at least ``min_worms_per_depth`` worms achieved.
    counts_per_depth = pb.obs.groupby("depth")["worm"].nunique()
    usable_depths = sorted(int(d) for d in counts_per_depth.index if counts_per_depth[d] >= min_worms_per_depth)
    print(f"[run] usable depths (worms ≥ {min_worms_per_depth}): {usable_depths}")

    results: dict[str, Any] = {"depths": usable_depths}
    for method in similarity_methods:
        print(f"\n[run] ==== method = {method} ====")
        sim_long = pairwise_similarity(
            pb, method=method,
            min_total_counts=min_total_counts,
            residualize_by_worm=residualize_by_worm,
        )
        sim_long.to_csv(out_dir / f"pairwise_similarity_{method}.csv", index=False)

        matches_per_depth: dict[int, pd.DataFrame] = {}
        observed_rows = []
        for d in usable_depths:
            m = orientation_match_at_depth(pb, sim_long, d)
            if m.empty:
                continue
            m.to_csv(out_dir / f"matches_depth{d}_{method}.csv", index=False)
            matches_per_depth[d] = m
            observed_rows.append(m)
        if not observed_rows:
            print(f"[run/{method}] no matches at any usable depth; skipping null.")
            continue
        observed = pd.concat(observed_rows, ignore_index=True)

        null = null_matched_similarity(
            run_dir, rna_path,
            method=method, n_null=n_null,
            counts_layer=counts_layer,
            min_total_counts=min_total_counts,
            residualize_by_worm=residualize_by_worm,
            trustworthy_only=trustworthy_only,
            seed=seed,
        )
        null.to_csv(out_dir / f"null_{method}.csv", index=False)

        summary = summarize(observed, null, out_dir, method=method)
        results[method] = {"summary": summary, "matches": matches_per_depth}

        groups = emit_cross_worm_groups(pb, sim_long, out_dir, method=method)
        results[method]["groups"] = groups

    print(f"\n[run] outputs → {out_dir}")
    return results
