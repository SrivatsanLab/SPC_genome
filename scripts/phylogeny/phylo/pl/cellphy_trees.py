"""CellPhy phylogram renderer with nodes colored by bootstrap support.

Reads a raxml-ng ``sup.raxml.supportFBP`` (or ``run.raxml.support``) newick —
internal node labels are the FBP bootstrap support values (0–100). Renders a
rectangular phylogram: root at left, leaves at right, branch lengths
proportional to substitutions per site. Internal-node scatter markers are
colored by support; low-support branches (default < 50) fade to gray so the
robust backbone stands out.

Companion to :mod:`phylo.pl.trees` (SVD trees), same styling conventions.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _load_tree(newick_path: str | Path):
    """Return an ete3 tree, rooted for display.

    Preference order:
      1. If a leaf named "ROOT" is present (synthetic outgroup added by the
         distance-tree runner), root by that outgroup.
      2. Otherwise, midpoint-root — appropriate default for unrooted ML/MP.
    """
    import ete3
    t = ete3.Tree(str(newick_path), format=0)
    root_leaves = [lf for lf in t.get_leaves() if lf.name == "ROOT"]
    if root_leaves:
        try:
            t.set_outgroup(root_leaves[0])
            return t
        except Exception:
            pass
    try:
        og = t.get_midpoint_outgroup()
        if og is not None:
            t.set_outgroup(og)
    except Exception:
        pass
    return t


def _layout(root) -> tuple[dict, list, float]:
    """Assign (x, y) to every node.

    ``x`` = cumulative branch length from the root (phylogram).
    ``y`` = leaf order (evenly spaced 0..n_leaves-1); internals get the
    mean of their children.
    """
    leaves = [lf for lf in root.get_leaves()]
    yorder = {id(lf): float(i) for i, lf in enumerate(leaves)}

    pos = {}

    def _walk(node, x_parent):
        x = x_parent + float(node.dist or 0.0)
        if node.is_leaf():
            y = yorder[id(node)]
            pos[id(node)] = (x, y)
            return y, x
        child_ys, child_xs = [], []
        for c in node.children:
            cy, cx = _walk(c, x)
            child_ys.append(cy)
            child_xs.append(cx)
        y = float(np.mean(child_ys))
        pos[id(node)] = (x, y)
        return y, max(child_xs)

    _, xmax = _walk(root, 0.0)
    return pos, leaves, xmax


def _support(node) -> float | None:
    """Return the FBP support attached to this internal node, or None."""
    if node.is_leaf() or node.is_root():
        return None
    s = getattr(node, "support", None)
    if s is None:
        return None
    try:
        f = float(s)
    except (TypeError, ValueError):
        return None
    if np.isnan(f):
        return None
    return f


def plot_cellphy_tree(
    newick_path: str | Path,
    out_path: str | Path | None = None,
    *,
    title: str | None = None,
    ax: plt.Axes | None = None,
    cmap_name: str = "viridis",
    support_threshold: float = 50.0,
    show_colorbar: bool = True,
) -> plt.Axes:
    """Plot a CellPhy tree with internal nodes colored by FBP support.

    Parameters
    ----------
    newick_path
        Path to a raxml-ng support newick (``sup.raxml.supportFBP``).
    support_threshold
        Internal edges with support < threshold are drawn light gray to
        emphasise the robust backbone; nodes are still colored by their
        raw support so the colormap remains informative.
    """
    root = _load_tree(newick_path)
    pos, leaves, xmax = _layout(root)
    n_leaves = len(leaves)

    from matplotlib.colors import Normalize
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=0.0, vmax=100.0)

    if ax is None:
        h = max(4.5, 0.16 * n_leaves + 1.5)
        fig, ax = plt.subplots(figsize=(9.0, h))
    else:
        fig = ax.figure

    def _edge_color(sup: float | None) -> str:
        if sup is None or sup < support_threshold:
            return "0.75"
        return "0.25"

    def _draw(node):
        px, py = pos[id(node)]
        for c in node.children:
            cx, cy = pos[id(c)]
            sup = _support(c)
            col = _edge_color(sup)
            lw = 2.6 if col != "0.75" else 1.8
            ax.plot([px, px], [py, cy], color=col, lw=lw, solid_capstyle="round")
            ax.plot([px, cx], [cy, cy], color=col, lw=lw, solid_capstyle="round")
            _draw(c)

    _draw(root)

    # Internal-node scatter, colored by support.
    sup_x, sup_y, sup_c = [], [], []
    for n in root.traverse():
        if n.is_leaf() or n.is_root():
            continue
        s = _support(n)
        if s is None:
            continue
        x, y = pos[id(n)]
        sup_x.append(x); sup_y.append(y); sup_c.append(s)
    if sup_x:
        ax.scatter(sup_x, sup_y, c=sup_c, cmap=cmap, norm=norm,
                   s=16, edgecolor="black", linewidth=0.35, zorder=5)

    # Leaves: small tick marks; no labels (dense).
    # Highlight a synthetic outgroup leaf (named "ROOT") with a distinct marker
    # so it doesn't disappear into the surrounding tick strip.
    non_root_leaves = [lf for lf in leaves if lf.name != "ROOT"]
    lx = [pos[id(lf)][0] for lf in non_root_leaves]
    ly = [pos[id(lf)][1] for lf in non_root_leaves]
    ax.scatter(lx, ly, color="black", s=5, marker="|", zorder=4)
    root_leaves = [lf for lf in leaves if lf.name == "ROOT"]
    for r in root_leaves:
        rx, ry = pos[id(r)]
        ax.scatter([rx], [ry], color="tab:red", s=60, marker="s",
                   edgecolor="black", linewidth=0.6, zorder=6)

    # Cosmetics.
    ax.set_xlim(-0.02 * xmax, xmax * 1.05)
    ax.set_ylim(-1, n_leaves)
    ax.set_yticks([])
    ax.set_xlabel("substitutions per site", fontsize=9)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    if title is None:
        title = f"{Path(newick_path).parent.parent.name} / {Path(newick_path).parent.name}: {n_leaves} cells"
    ax.set_title(title, fontsize=11)

    if show_colorbar:
        from matplotlib.cm import ScalarMappable
        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, aspect=25)
        cbar.set_label("FBP bootstrap support", fontsize=9)
        cbar.ax.tick_params(labelsize=8)

    if out_path is not None:
        fig.tight_layout()
        out_path = Path(out_path)
        # Save both a raster (PNG) and a vector (SVG) so figures are both
        # inline-viewable and cleanly rescalable for publication.
        for suffix in (".png", ".svg"):
            fig.savefig(out_path.with_suffix(suffix), dpi=150, bbox_inches="tight")
        plt.close(fig)
    return ax


def plot_all_distance_trees(
    inference_dir: str | Path,
    out_dir: str | Path,
    *,
    metric: str = "soft",
    pipelines: tuple[str, ...] = ("nj", "nj_nni", "bme_nni", "upgma"),
    rootings: tuple[str, ...] = ("og", "mad"),
    cmap_name: str = "viridis",
    support_threshold: float = 50.0,
) -> None:
    """Render per-worm distance trees for the v2 filename convention.

    Looks for newicks at
      ``<inference_dir>/distance/worm_*/{metric}__{pipeline}__{rooting}.support``
    (falls back to the legacy v1 layout when v2 outputs are missing).
    """
    inference_dir = Path(inference_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dist_root = inference_dir / "distance"

    for pipeline in pipelines:
        for rooting in rootings:
            rendered: list[tuple[str, Path]] = []
            for wd in sorted(dist_root.glob("worm_*")):
                cand = [wd / f"{metric}__{pipeline}__{rooting}.support"]
                # v1 fallback (only for the vanilla nj / upgma names)
                if pipeline in ("nj", "upgma") and rooting == "og":
                    cand.append(wd / f"{metric}__{pipeline}.support")
                    cand.append(wd / f"{metric}__{pipeline}.newick")
                tree_path = next((p for p in cand if p.exists()), None)
                if tree_path is None:
                    continue
                worm = wd.name
                out_png = out_dir / f"{worm}__dist_{metric}_{pipeline}_{rooting}.png"
                plot_cellphy_tree(
                    tree_path, out_png,
                    title=f"{worm} · dist:{metric}:{pipeline}:{rooting}",
                    cmap_name=cmap_name,
                    support_threshold=support_threshold,
                )
                rendered.append((worm, tree_path))
            if not rendered:
                continue
            n = len(rendered)
            ncols = 4
            nrows = int(np.ceil(n / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.2 * nrows))
            axes_flat = axes.ravel() if hasattr(axes, "ravel") else [axes]
            for i, (worm, tp) in enumerate(rendered):
                plot_cellphy_tree(
                    tp, out_path=None, ax=axes_flat[i],
                    title=worm, cmap_name=cmap_name,
                    support_threshold=support_threshold,
                    show_colorbar=False,
                )
            for j in range(n, len(axes_flat)):
                axes_flat[j].axis("off")
            fig.suptitle(f"Distance {pipeline} · rooting={rooting} · metric={metric}",
                         fontsize=13)
            fig.tight_layout()
            for suf in (".png", ".svg"):
                fig.savefig(out_dir / f"overview__dist_{metric}_{pipeline}_{rooting}{suf}",
                            dpi=140, bbox_inches="tight")
            plt.close(fig)
            print(f"[plot_dist] overview__dist_{metric}_{pipeline}_{rooting}.{{png,svg}}")


def plot_all_mp_trees(
    inference_dir: str | Path,
    out_dir: str | Path,
    *,
    cmap_name: str = "viridis",
    support_threshold: float = 50.0,
) -> None:
    """Render per-worm MP majority-consensus trees, PNG+SVG plus overview.

    Looks for ``<inference_dir>/mp/worm_*/real.majority.newick`` (bootstrap
    support annotated as internal-node labels).
    """
    inference_dir = Path(inference_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mp_root = inference_dir / "mp"
    rendered: list[tuple[str, Path]] = []
    for wd in sorted(mp_root.glob("worm_*")):
        cand = [wd / "real.majority.newick",
                wd / "real.mp_support.newick",
                wd / "real.mp.newick"]
        tree_path = next((p for p in cand if p.exists()), None)
        if tree_path is None:
            continue
        worm = wd.name
        out_png = out_dir / f"{worm}__mp.png"
        plot_cellphy_tree(
            tree_path, out_png,
            title=f"{worm} · MP · {tree_path.name}",
            cmap_name=cmap_name,
            support_threshold=support_threshold,
        )
        rendered.append((worm, tree_path))
        print(f"[plot_mp] {out_png}")
    if not rendered:
        return
    n = len(rendered)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.2 * nrows))
    axes_flat = axes.ravel() if hasattr(axes, "ravel") else [axes]
    for i, (worm, tp) in enumerate(rendered):
        plot_cellphy_tree(
            tp, out_path=None, ax=axes_flat[i],
            title=worm, cmap_name=cmap_name,
            support_threshold=support_threshold,
            show_colorbar=False,
        )
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle("MP-CS majority-rule trees", fontsize=13)
    fig.tight_layout()
    for suffix in (".png", ".svg"):
        fig.savefig(out_dir / f"overview__mp{suffix}", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_mp] overview__mp.{{png,svg}}")


def plot_all_cellphy_trees(
    cellphy_dir: str | Path,
    out_dir: str | Path,
    *,
    matrix: str = "real",
    cmap_name: str = "viridis",
    support_threshold: float = 50.0,
) -> None:
    """Render one PNG per worm under ``cellphy_dir/worm_*/<matrix>/``.

    Also emits an ``overview.png`` grid.
    """
    cellphy_dir = Path(cellphy_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    worm_dirs = sorted(cellphy_dir.glob(f"worm_*/{matrix}"))
    if not worm_dirs:
        print(f"[plot_all_cellphy_trees] no {matrix} runs under {cellphy_dir}")
        return

    rendered = []
    for wd in worm_dirs:
        cand = [wd / "sup.raxml.supportFBP", wd / "run.raxml.support", wd / "run.raxml.bestTree"]
        tree_path = next((p for p in cand if p.exists()), None)
        if tree_path is None:
            print(f"[plot_all_cellphy_trees] no tree in {wd}")
            continue
        worm = wd.parent.name  # worm_worm07
        out_png = out_dir / f"{worm}__{matrix}.png"
        plot_cellphy_tree(
            tree_path, out_png,
            title=f"{worm} · {matrix} · {tree_path.name}",
            cmap_name=cmap_name,
            support_threshold=support_threshold,
        )
        print(f"[plot_all_cellphy_trees] {out_png}")
        rendered.append((worm, tree_path))

    if not rendered:
        return

    # Overview grid.
    n = len(rendered)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.2 * nrows))
    axes_flat = axes.ravel() if hasattr(axes, "ravel") else [axes]
    for i, (worm, tp) in enumerate(rendered):
        plot_cellphy_tree(
            tp, out_path=None, ax=axes_flat[i],
            title=worm, cmap_name=cmap_name,
            support_threshold=support_threshold,
            show_colorbar=False,
        )
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle(f"CellPhy trees · matrix={matrix}", fontsize=13)
    fig.tight_layout()
    for suffix in (".png", ".svg"):
        fig.savefig(out_dir / f"overview__{matrix}{suffix}", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_all_cellphy_trees] overview__{matrix}.{{png,svg}}")
