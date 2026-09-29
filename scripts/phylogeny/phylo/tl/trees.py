"""Tree data model + builders + recursion driver.

Design:

- **Data model** (`Tree`, `TreeNode`, `SplitInfo`) — plain dataclasses,
  serializable to JSON. No topology library dependency; export to Newick if
  needed.
- **Builder ABC** (`TreeBuilder`) — one method: ``bipart(adata) -> SplitInfo``
  or None if no defensible split. Recursion driver knows nothing about the
  builder's internals.
- **`SVDBipartitionBuilder`** — default. Wraps :mod:`phylo.tl.svd_kmeans` at
  ``K=2`` with core-variant selection.
- **`build_worm_tree`** — recursion driver. Stop conditions are size + depth
  only; validation gates run in a separate stage and annotate nodes.
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.sparse import issparse

from .mquad_select import select_informative
from .svd_kmeans import build_matrix, fit_modules, score_cells


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class SplitInfo:
    """Metadata for one bipartition."""

    n_variants_used: int
    n_core: int
    module_variants: dict[str, list[str]]
    module_pooled_vaf: dict[str, float]         # within-clade pooled VAF, for §6.1 amplitude check
    module_n_cells: dict[str, int]
    module_purity_median: dict[str, float]
    method: str = "svd_bipart"
    filter_mode: str | None = None              # which core-variant filter actually ran at this node


@dataclass
class TreeNode:
    """A node in a per-worm bipartition tree.

    ``path`` is a dotted string: root = ``""``; children = ``"L"``, ``"R"``;
    grandchildren = ``"L.L"``, ``"L.R"``, etc. ``depth == len(path.split('.'))``
    for non-root; root depth is 0.
    """

    path: str
    cells: list[str]
    depth: int
    left: TreeNode | None = None
    right: TreeNode | None = None
    split: SplitInfo | None = None                # None ⇒ leaf
    stop_reason: str | None = None                # populated on leaves

    @property
    def is_leaf(self) -> bool:
        return self.split is None

    def leaves(self) -> list[TreeNode]:
        if self.is_leaf:
            return [self]
        return self.left.leaves() + self.right.leaves()

    def walk(self):
        yield self
        if not self.is_leaf:
            yield from self.left.walk()
            yield from self.right.walk()


@dataclass
class Tree:
    """A per-worm bipartition tree."""

    worm: str
    root: TreeNode
    builder: str
    params: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        def _serialize(node: TreeNode) -> dict:
            d = {
                "path": node.path,
                "depth": node.depth,
                "n_cells": len(node.cells),
                "cells": node.cells,
                "stop_reason": node.stop_reason,
            }
            if node.split is not None:
                d["split"] = asdict(node.split)
                d["left"] = _serialize(node.left)
                d["right"] = _serialize(node.right)
            return d

        return {
            "worm": self.worm,
            "builder": self.builder,
            "params": self.params,
            "root": _serialize(self.root),
        }

    def to_json(self, path: str | Path) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, default=str)

    def cell_assignments(self) -> list[tuple[str, str]]:
        """List of ``(cell_id, leaf_path)`` for every cell in the tree."""
        out = []
        for leaf in self.root.leaves():
            for c in leaf.cells:
                out.append((c, leaf.path or "root"))
        return out

    def stats(self) -> dict:
        leaves = self.root.leaves()
        splits = [n for n in self.root.walk() if not n.is_leaf]
        return {
            "worm": self.worm,
            "n_cells": len(self.root.cells),
            "n_leaves": len(leaves),
            "n_splits": len(splits),
            "max_depth": max((n.depth for n in self.root.walk()), default=0),
            "smallest_leaf": min((len(n.cells) for n in leaves), default=0),
        }


# ---------------------------------------------------------------------------
# Builder ABC + SVD bipartition implementation
# ---------------------------------------------------------------------------


class TreeBuilder(ABC):
    """One split. The recursion driver calls this once per node."""

    name: str

    @abstractmethod
    def bipart(self, adata: ad.AnnData, node: TreeNode) -> SplitInfo | None:
        """Compute a bipartition of ``adata`` (already subset to this node's cells).

        Return None if no defensible split exists (e.g. too few variants).
        Otherwise return a :class:`SplitInfo`. The recursion driver reads
        ``module_variants`` to assign cells to sides.
        """
        raise NotImplementedError


class SVDBipartitionBuilder(TreeBuilder):
    """Recursive K=2 wrapper around :mod:`phylo.tl.svd_kmeans`.

    At each node:

    1. Re-center the AD/DP matrix over this node's cells (variant means
       shift as cells subset).
    2. Truncated SVD (`ncomp` components) → K-means at K=2.
    3. Optional core-variant selection: keep the ``n_core`` variants closest
       to their assigned centroid, refit.
    4. Score cells against each module (depth-weighted pooled VAF); assign
       by argmax.
    5. Compute within-clade pooled VAF per module (used later by the §6.1
       VAF-amplitude check).
    """

    name = "svd_bipart"

    def __init__(
        self,
        *,
        ncomp: int = 2,
        n_core: int | None = None,
        min_dp_var: int = 1,
        min_variants_to_split: int = 20,
        seed: int = 0,
        core_selection: str = "none",
        mquad_delta_bic_threshold: float = 10.0,
        min_carriers_at_node: int = 2,
        mquad_nproc: int = 1,
        hybrid_switch_cells: int = 20,
        setop_min_carriers: int = 2,
    ):
        if ncomp < 1:
            raise ValueError("ncomp must be ≥ 1")
        # ncomp=1 is defensible at K=2 per HAPLOTYPE_ANALYSIS_SUMMARY.md §4
        # (per-variant centering drops the signal rank to K-1 = 1); higher
        # ncomp adds nuisance dimensions (coverage / α gradient / ploidy).
        # See TREE_BUILDER_METHOD_REVIEW_RESULTS_1.md §2.4.
        if core_selection not in ("none", "centroid", "mquad", "setop", "hybrid"):
            raise ValueError(f"core_selection must be one of none/centroid/mquad/setop/hybrid, got {core_selection!r}")
        self.ncomp = ncomp
        self.n_core = n_core
        self.min_dp_var = min_dp_var
        self.min_variants_to_split = min_variants_to_split
        self.seed = seed
        self.core_selection = core_selection
        self.mquad_delta_bic_threshold = mquad_delta_bic_threshold
        self.min_carriers_at_node = min_carriers_at_node
        self.mquad_nproc = mquad_nproc
        self.hybrid_switch_cells = hybrid_switch_cells
        self.setop_min_carriers = setop_min_carriers

    def bipart(self, adata: ad.AnnData, node: TreeNode) -> SplitInfo | None:
        AD = adata.layers["AD"]
        DP = adata.layers["DP"]
        if issparse(AD):
            AD = AD.toarray()
        if issparse(DP):
            DP = DP.toarray()

        # own_bulk_vaf is precomputed by panels stage; if absent, recompute
        # on-the-fly (allows use on ad-hoc AnnData subsets).
        if "own_bulk_vaf" in adata.var:
            AF = adata.var["own_bulk_vaf"].to_numpy()
        else:
            AF = np.asarray(DP.sum(0))  # non-zero to pass build_matrix's mask
            AF = np.where(AF > 0, AD.sum(0) / np.maximum(AF, 1), 0.0)

        # Per-node re-filter: drop variants that fail the ≥ min_carriers_at_node
        # rule *among this node's cells*. At depth ≥ 1 most panel variants are
        # markers for the other subtree and contribute pure noise.
        n_carriers_at_node = (AD > 0).sum(0)
        node_carrier_mask = n_carriers_at_node >= self.min_carriers_at_node
        if node_carrier_mask.sum() < self.min_variants_to_split:
            return None

        AD_n = AD[:, node_carrier_mask]
        DP_n = DP[:, node_carrier_mask]
        AF_n = AF[node_carrier_mask] if AF is not None else None
        node_var_offsets = np.flatnonzero(node_carrier_mask)

        # No further VAF filtering — the panel is already filtered.
        X, inf_idx_local, stats = build_matrix(
            AD_n, DP_n, AF_n, vaf_lo=-1.0, vaf_hi=2.0, min_dp_var=self.min_dp_var
        )
        if stats["n_var"] < self.min_variants_to_split:
            return None
        # Translate inf_idx (into AD_n columns) back to panel-variant indices.
        inf_idx = node_var_offsets[inf_idx_local]

        labels, d_all = fit_modules(X, ncomp=self.ncomp, K=2, seed=self.seed)

        # Core selection. `filter_mode` records which branch actually ran (useful for
        # hybrid diagnostics + fallback tracking).
        filter_mode: str
        # `hybrid` dispatches by clade size: mquad has no power below ~20 cells
        # (BBMix can't distinguish 1- vs 2-component under the log(n) BIC penalty),
        # so we swap to a set-op informativeness filter down there.
        effective_mode = self.core_selection
        if effective_mode == "hybrid":
            effective_mode = "mquad" if adata.n_obs >= self.hybrid_switch_cells else "setop"

        if effective_mode == "mquad":
            # Per-node seed: base seed XOR a *deterministic* hash of the node path
            # so sibling nodes don't share BBMix init noise AND repeat runs of the
            # same config give identical trees. Python's built-in ``hash()`` on
            # strings is randomised by PYTHONHASHSEED, which was the source of
            # cross-run non-determinism before this fix.
            import zlib
            node_seed = self.seed ^ (zlib.crc32(node.path.encode()) & 0x7FFFFFFF)
            keep_local, delta_bic = select_informative(
                AD[:, inf_idx],
                DP[:, inf_idx],
                delta_bic_threshold=self.mquad_delta_bic_threshold,
                min_dp=1,
                min_ad=1,
                nproc=self.mquad_nproc,
                seed=int(node_seed),
            )
            if keep_local.sum() < self.min_variants_to_split:
                labels_final = labels
                var_idx_final = inf_idx
                n_core_used = int(len(inf_idx))
                filter_mode = "mquad_fallback"
            else:
                X_core = X[:, keep_local]
                labels_final, _ = fit_modules(X_core, ncomp=self.ncomp, K=2, seed=self.seed)
                var_idx_final = inf_idx[keep_local]
                n_core_used = int(keep_local.sum())
                filter_mode = "mquad"
        elif effective_mode == "setop":
            # Tentative cell partition from initial K-means labels.
            _, mods_t, assign_t = score_cells(AD, DP, inf_idx, labels, self.min_dp_var)
            L_mask = assign_t == mods_t[0]
            R_mask = assign_t == mods_t[1]
            # Per-variant carrier counts on each tentative side.
            AD_inf = AD[:, inf_idx]
            n_L = (AD_inf[L_mask] > 0).sum(0)
            n_R = (AD_inf[R_mask] > 0).sum(0)
            k = self.setop_min_carriers
            keep_local = ((n_L >= k) & (n_R == 0)) | ((n_R >= k) & (n_L == 0))
            if keep_local.sum() < self.min_variants_to_split:
                labels_final = labels
                var_idx_final = inf_idx
                n_core_used = int(len(inf_idx))
                filter_mode = "setop_fallback"
            else:
                X_core = X[:, keep_local]
                labels_final, _ = fit_modules(X_core, ncomp=self.ncomp, K=2, seed=self.seed)
                var_idx_final = inf_idx[keep_local]
                n_core_used = int(keep_local.sum())
                filter_mode = "setop"
        elif effective_mode == "centroid" and self.n_core is not None and self.n_core < len(inf_idx):
            core = np.sort(np.argsort(d_all)[: self.n_core])
            X_core = X[:, core]
            labels_final, _ = fit_modules(X_core, ncomp=self.ncomp, K=2, seed=self.seed)
            var_idx_final = inf_idx[core]
            n_core_used = int(core.size)
            filter_mode = "centroid"
        else:
            labels_final = labels
            var_idx_final = inf_idx
            n_core_used = int(len(inf_idx))
            filter_mode = "none"

        W, modules, assignment = score_cells(
            AD, DP, var_idx_final, labels_final, self.min_dp_var
        )
        # Everyone got some weight? If a module is empty, no split.
        if len(np.unique(assignment)) < 2:
            return None

        var_names = adata.var_names.to_numpy().astype(str)
        module_variants: dict[str, list[str]] = {}
        module_pooled_vaf: dict[str, float] = {}
        module_n_cells: dict[str, int] = {}
        module_purity_median: dict[str, float] = {}
        for side_label, mod in zip(("L", "R"), modules):
            mvars = var_idx_final[labels_final == mod]
            in_side = assignment == mod
            # within-clade pooled VAF over this module's variants
            ad_sum = AD[in_side][:, mvars].sum()
            dp_sum = DP[in_side][:, mvars].sum()
            module_variants[side_label] = var_names[mvars].tolist()
            module_pooled_vaf[side_label] = float(ad_sum / max(dp_sum, 1))
            module_n_cells[side_label] = int(in_side.sum())
            module_purity_median[side_label] = float(np.median(W[in_side, np.where(modules == mod)[0][0]]))

        return SplitInfo(
            n_variants_used=int(len(inf_idx)),
            n_core=n_core_used,
            module_variants=module_variants,
            module_pooled_vaf=module_pooled_vaf,
            module_n_cells=module_n_cells,
            module_purity_median=module_purity_median,
            method=self.name,
            filter_mode=filter_mode,
        )


# ---------------------------------------------------------------------------
# Recursion driver
# ---------------------------------------------------------------------------


def _assign_cells_to_sides(
    adata: ad.AnnData, split: SplitInfo, min_dp_var: int
) -> tuple[list[str], list[str]]:
    """Deterministic cell → side assignment using the split's module variants."""
    AD = adata.layers["AD"]
    DP = adata.layers["DP"]
    if issparse(AD):
        AD = AD.toarray()
    if issparse(DP):
        DP = DP.toarray()
    var_names = adata.var_names.to_numpy().astype(str)
    idx = {v: i for i, v in enumerate(var_names)}
    L_vars = np.array([idx[v] for v in split.module_variants["L"] if v in idx])
    R_vars = np.array([idx[v] for v in split.module_variants["R"] if v in idx])
    labels_lab = np.concatenate(
        [np.full(L_vars.size, 1, dtype=int), np.full(R_vars.size, 2, dtype=int)]
    )
    var_idx = np.concatenate([L_vars, R_vars])
    _, modules, assignment = score_cells(AD, DP, var_idx, labels_lab, min_dp_var)
    cells = adata.obs_names.to_numpy().astype(str)
    L_side = cells[assignment == modules[0]].tolist()
    R_side = cells[assignment == modules[1]].tolist()
    return L_side, R_side


def build_worm_tree(
    adata: ad.AnnData,
    builder: TreeBuilder,
    *,
    min_clade_size: int = 10,
    max_depth: int = 6,
    worm: str | None = None,
    params: dict | None = None,
) -> Tree:
    """Recursively bipartition ``adata`` into a per-worm tree.

    Parameters
    ----------
    adata : AnnData
        Per-worm somatic panel from :func:`phylo.pp.build_panels`. All cells
        assumed already QC-filtered and assigned to this worm.
    builder : TreeBuilder
        Any implementation. :class:`SVDBipartitionBuilder` is the default.
    min_clade_size : int
        Do not split below this cell count; the node becomes a leaf.
    max_depth : int
        Hard depth cap.
    worm : str, optional
        Worm ID for record-keeping. Falls back to ``adata.uns['panel_worm']``.
    """
    worm = worm or adata.uns.get("panel_worm") or "unknown"
    all_cells = adata.obs_names.astype(str).tolist()
    root = TreeNode(path="", cells=all_cells, depth=0)

    def _recurse(node: TreeNode) -> None:
        if len(node.cells) < min_clade_size:
            node.stop_reason = f"clade_size<{min_clade_size}"
            return
        if node.depth >= max_depth:
            node.stop_reason = f"depth>={max_depth}"
            return
        sub = adata[node.cells].copy()
        split = builder.bipart(sub, node)
        if split is None:
            node.stop_reason = "no_split"
            return
        L_cells, R_cells = _assign_cells_to_sides(sub, split, min_dp_var=1)
        if min(len(L_cells), len(R_cells)) < min_clade_size:
            node.stop_reason = f"child_size<{min_clade_size}"
            return
        node.split = split
        node.left = TreeNode(
            path=(f"{node.path}.L" if node.path else "L"),
            cells=L_cells,
            depth=node.depth + 1,
        )
        node.right = TreeNode(
            path=(f"{node.path}.R" if node.path else "R"),
            cells=R_cells,
            depth=node.depth + 1,
        )
        _recurse(node.left)
        _recurse(node.right)

    _recurse(root)
    return Tree(worm=worm, root=root, builder=builder.name, params=params or {})
