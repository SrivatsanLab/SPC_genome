"""Per-worm bipartition tree plots.

Layout: leaves are placed left-to-right in traversal order; internal nodes sit
at the midpoint of their children; y = depth (root at top). Each split edge
is annotated with the within-clade pooled VAF (``module_pooled_vaf``); each
node with its cell count. Node markers are colored by pooled VAF using a
diverging colormap centered near the plan §6.1 target (0.35).
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ..tl.bipartitions import compute_confidence
from ..tl.trees import Tree, TreeNode


def _from_dict(d: dict) -> TreeNode:
    n = TreeNode(path=d["path"], cells=d.get("cells", []), depth=d["depth"])
    n.stop_reason = d.get("stop_reason")
    if "split" in d:
        # Not reconstructing SplitInfo dataclass; keep the dict on ._split_dict.
        n._split_dict = d["split"]                                            # type: ignore[attr-defined]
        n.left = _from_dict(d["left"])
        n.right = _from_dict(d["right"])
    return n


def load_tree(path: str | Path) -> tuple[str, TreeNode, dict]:
    """Reload a Tree from its JSON. Returns ``(worm, root, params)``."""
    with open(path) as f:
        d = json.load(f)
    root = _from_dict(d["root"])
    return d["worm"], root, d.get("params", {})


def _layout(root: TreeNode) -> tuple[dict[str, tuple[float, float]], list[TreeNode]]:
    """Assign (x, y) to every node. Leaves get integer x; internals average children."""
    leaves: list[TreeNode] = []

    def _collect_leaves(n: TreeNode) -> None:
        if getattr(n, "_split_dict", None) is None and n.left is None:
            leaves.append(n)
            return
        _collect_leaves(n.left)
        _collect_leaves(n.right)

    _collect_leaves(root)
    xs = {leaf.path or "root": float(i) for i, leaf in enumerate(leaves)}

    def _place(n: TreeNode) -> float:
        key = n.path or "root"
        if key in xs:
            return xs[key]
        xl = _place(n.left)
        xr = _place(n.right)
        xs[key] = (xl + xr) / 2
        return xs[key]

    _place(root)
    pos = {k: (x, -float(next_n.depth)) for k, x in xs.items() for next_n in [_find(root, k)]}
    return pos, leaves


def _find(root: TreeNode, path_key: str) -> TreeNode:
    """Locate node by its stored path key (root path is ``''``, we key as ``'root'``)."""
    target = "" if path_key == "root" else path_key
    for n in _walk(root):
        if n.path == target:
            return n
    raise KeyError(path_key)


def _walk(n: TreeNode):
    yield n
    if n.left is not None:
        yield from _walk(n.left)
    if n.right is not None:
        yield from _walk(n.right)


def _clade_confidences(root: TreeNode) -> dict[str, float]:
    """Attach a confidence to every non-root node, computed from its parent's split.

    Returns ``{node_path: confidence}``; root gets no entry (it isn't a child
    of anything). Uses the same formula as :func:`phylo.tl.bipartitions.compute_confidence`
    with default equal weights.
    """
    out: dict[str, float] = {}

    def walk(node: TreeNode) -> None:
        s = getattr(node, "_split_dict", None)
        if s is None:
            return
        vaf = s["module_pooled_vaf"]
        n_cells = s["module_n_cells"]
        purity = s["module_purity_median"]
        min_size = min(n_cells["L"], n_cells["R"])
        for side, child in (("L", node.left), ("R", node.right)):
            other = "R" if side == "L" else "L"
            conf, _ = compute_confidence(
                n_core=s["n_core"],
                n_variants_used=s["n_variants_used"],
                pooled_vaf=vaf[side],
                sibling_vaf=vaf[other],
                min_side_size=min_size,
                purity_median=purity[side],
            )
            out[child.path] = conf
            walk(child)

    walk(root)
    return out



def _annotate_purity(root: TreeNode, panel_h5ad: str | Path) -> None:
    """Mutate each internal node's ``_split_dict`` to add cross-side pooled VAF.

    Adds fields per split:
      vaf_L_at_moduleL  = pooled VAF of module-L variants restricted to L-side cells
      vaf_R_at_moduleL  = pooled VAF of module-L variants restricted to R-side cells
      vaf_R_at_moduleR / vaf_L_at_moduleR (symmetric)
      purity_L = vaf_L_at_moduleL - vaf_R_at_moduleL   # signal contrast for L module
      purity_R = vaf_R_at_moduleR - vaf_L_at_moduleR   # signal contrast for R module
    """
    import anndata as ad
    from scipy.sparse import issparse

    a = ad.read_h5ad(panel_h5ad)
    AD = a.layers["AD"].toarray() if issparse(a.layers["AD"]) else a.layers["AD"]
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else a.layers["DP"]
    cell_idx = {c: i for i, c in enumerate(a.obs_names.astype(str))}
    var_idx = {v: i for i, v in enumerate(a.var_names.astype(str))}

    def _pooled(cell_ids, var_ids):
        cs = [cell_idx[c] for c in cell_ids if c in cell_idx]
        vs = [var_idx[v] for v in var_ids if v in var_idx]
        if not cs or not vs:
            return float("nan")
        a_sum = AD[np.ix_(cs, vs)].sum()
        d_sum = DP[np.ix_(cs, vs)].sum()
        return float(a_sum / d_sum) if d_sum else 0.0

    for node in _walk(root):
        s = getattr(node, "_split_dict", None)
        if s is None:
            continue
        L_cells = node.left.cells if node.left else []
        R_cells = node.right.cells if node.right else []
        mL_vars = s["module_variants"]["L"]
        mR_vars = s["module_variants"]["R"]
        s["vaf_L_at_moduleL"] = _pooled(L_cells, mL_vars)
        s["vaf_R_at_moduleL"] = _pooled(R_cells, mL_vars)
        s["vaf_R_at_moduleR"] = _pooled(R_cells, mR_vars)
        s["vaf_L_at_moduleR"] = _pooled(L_cells, mR_vars)
        s["purity_L"] = s["vaf_L_at_moduleL"] - s["vaf_R_at_moduleL"]
        s["purity_R"] = s["vaf_R_at_moduleR"] - s["vaf_L_at_moduleR"]


def plot_tree(
    tree_or_json: Tree | str | Path,
    out_path: str | Path | None = None,
    *,
    title: str | None = None,
    ax: plt.Axes | None = None,
    cmap_name: str = "viridis",
    show_colorbar: bool = True,
    vmax: float | None = None,
    panel_h5ad: str | Path | None = None,
) -> plt.Axes:
    """Plot a per-worm bipartition tree with confidence-colored nodes.

    Every non-root node marker is colored by the *child clade's* confidence
    score (see :func:`phylo.tl.bipartitions.compute_confidence`). Root node is
    white. Branch edges are gray to keep the figure readable.

    Parameters
    ----------
    vmax : float or None
        Colormap upper bound. ``None`` (default) uses the max confidence
        observed in this tree, so the palette spans the full range of
        meaningful values. Pass an explicit value (e.g. the max across a set
        of trees) for cross-tree comparability.
    """
    if isinstance(tree_or_json, (str, Path)):
        worm, root, params = load_tree(tree_or_json)
    else:
        worm = tree_or_json.worm
        # Convert to dict-form nodes to reuse layout code
        d = tree_or_json.to_dict()
        root = _from_dict(d["root"])
        params = d.get("params", {})

    if panel_h5ad is not None:
        _annotate_purity(root, panel_h5ad)

    pos, leaves = _layout(root)
    confidences = _clade_confidences(root)

    cmap = plt.get_cmap(cmap_name)
    from matplotlib.colors import Normalize
    if vmax is None:
        vmax = max(confidences.values()) if confidences else 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(6.0, 0.7 * len(leaves) + 2), 4.5))
    else:
        fig = ax.figure

    # Edges (uniform gray) + branch annotations (within-clade pooled VAF + module size).
    # Horizontal bar of each split sits at the parent's y so the parent node
    # marker overlays the T-junction instead of hovering above it.
    for n in _walk(root):
        if getattr(n, "_split_dict", None) is None:
            continue
        s = n._split_dict                                                     # type: ignore[attr-defined]
        n_mod_L = len(s["module_variants"]["L"])
        n_mod_R = len(s["module_variants"]["R"])
        vaf_L = s["module_pooled_vaf"]["L"]
        vaf_R = s["module_pooled_vaf"]["R"]
        px, py = pos[n.path or "root"]
        lx, ly = pos[n.left.path]
        rx, ry = pos[n.right.path]
        ax.plot([lx, rx], [py, py], color="0.4", lw=1.5)
        ax.plot([lx, lx], [py, ly], color="0.4", lw=1.5)
        ax.plot([rx, rx], [py, ry], color="0.4", lw=1.5)
        pur_L = s.get("purity_L")
        pur_R = s.get("purity_R")
        lab_L = f"vaf={vaf_L:.2f}\nn_mod={n_mod_L}"
        lab_R = f"vaf={vaf_R:.2f}\nn_mod={n_mod_R}"
        if pur_L is not None and not (isinstance(pur_L, float) and pur_L != pur_L):
            lab_L += f"\npur={pur_L:+.2f}"
        if pur_R is not None and not (isinstance(pur_R, float) and pur_R != pur_R):
            lab_R += f"\npur={pur_R:+.2f}"
        ax.text(lx - 0.05, (py + ly) / 2, lab_L,
                ha="right", va="center", fontsize=7, color="0.15")
        ax.text(rx + 0.05, (py + ry) / 2, lab_R,
                ha="left", va="center", fontsize=7, color="0.15")

    # Nodes colored by their own (child-of-parent) confidence.
    # Root has no parent → colored white with black outline.
    for n in _walk(root):
        x, y = pos[n.path or "root"]
        is_leaf = getattr(n, "_split_dict", None) is None
        is_root = (n.path == "")
        if is_root:
            color = "white"
        else:
            color = cmap(norm(confidences.get(n.path, 0.0)))
        size = 90 if not is_leaf else 70
        ax.scatter([x], [y], s=size, color=color, edgecolor="black", linewidth=0.8, zorder=5)
        if is_leaf:
            ax.text(x, y - 0.25, f"{len(n.cells)}", ha="center", va="top",
                    fontsize=9, color="0.15", fontweight="bold")

    # Cosmetics
    ax.set_xticks([])
    ax.set_yticks(list(range(0, -int(max(-y for _, y in pos.values())) - 1, -1)))
    ax.set_yticklabels([f"d={-i}" for i in ax.get_yticks()])
    ax.set_xlim(-0.7, max(x for x, _ in pos.values()) + 0.7)
    ax.set_ylim(min(y for _, y in pos.values()) - 1.2, 0.6)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title(title or f"{worm}: {len(leaves)} leaves, {len(list(_walk(root))) - len(leaves)} splits",
                 fontsize=11)

    if show_colorbar:
        from matplotlib.cm import ScalarMappable
        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02, aspect=25)
        cbar.set_label("clade confidence", fontsize=9)
        cbar.ax.tick_params(labelsize=8)

    if out_path is not None:
        fig.tight_layout()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return ax


def plot_all_trees(run_dir: str | Path, cmap_name: str = "viridis") -> None:
    """Render one PNG per worm plus an overview grid, into ``run_dir/trees/``.

    The colormap upper bound is the maximum clade confidence observed across
    all trees in this run, so the palette spans the full range of meaningful
    values and per-worm PNGs remain visually comparable to each other and to
    the grid view.
    """
    run_dir = Path(run_dir)
    trees_dir = run_dir / "trees"
    tree_jsons = sorted(trees_dir.glob("worm_*/tree.json"))
    if not tree_jsons:
        print(f"[plot_all_trees] no tree.json under {trees_dir}")
        return

    # Compute a run-wide vmax so all trees share the same scale.
    run_vmax = 0.0
    for tj in tree_jsons:
        _, root, _ = load_tree(tj)
        confs = _clade_confidences(root)
        if confs:
            run_vmax = max(run_vmax, max(confs.values()))
    if run_vmax <= 0.0:
        run_vmax = 1.0

    for tj in tree_jsons:
        plot_tree(tj, tj.parent / "tree.png", cmap_name=cmap_name, vmax=run_vmax)

    # Overview grid — suppress per-panel colorbars, add one shared bar.
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    n = len(tree_jsons)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows))
    axes = np.atleast_2d(axes)
    for i, tj in enumerate(tree_jsons):
        ax = axes[i // ncols, i % ncols]
        plot_tree(tj, out_path=None, ax=ax, cmap_name=cmap_name, show_colorbar=False, vmax=run_vmax)
    for j in range(n, nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")
    fig.tight_layout()
    # Shared colorbar for the whole grid.
    norm = Normalize(vmin=0.0, vmax=run_vmax)
    sm = ScalarMappable(norm=norm, cmap=plt.get_cmap(cmap_name))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.6, aspect=40, pad=0.02)
    cbar.set_label("clade confidence", fontsize=11)
    cbar.ax.tick_params(labelsize=10)
    fig.savefig(trees_dir / "all_trees.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_all_trees] wrote {len(tree_jsons)} per-worm PNGs + all_trees.png (vmax={run_vmax:.2f})")
