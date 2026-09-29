#!/usr/bin/env python3
"""Cross-tree diagnostics for the CellPhy sweep.

For one worm, given the three CellPhy trees ({real, covmask, shuffled}), compute:

  bipartitions.tsv        every non-trivial split in each tree (with support)
  overlap.tsv             pairwise exact-match + max-Jaccard against `real`
  permutation_null.tsv    cell-label permutation null for max-Jaccard

The null: within this worm's cell set, permute leaf labels of the tree under
comparison K times and recompute the mean max-Jaccard against the reference.
Reports the observed value, the null mean, and the fraction of null draws
≥ observed (empirical p).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import ete3
import numpy as np
import pandas as pd


def _load_newick(p: Path) -> ete3.Tree | None:
    if not p.exists():
        return None
    nw = p.read_text().strip()
    if not nw or nw == ";":
        return None
    for fmt in (0, 1):
        try:
            return ete3.Tree(nw, format=fmt)
        except Exception:
            continue
    return None


def bipartitions(tree: ete3.Tree, all_taxa: set[str]) -> list[tuple[frozenset, float]]:
    out = []
    for node in tree.traverse():
        if node.is_leaf() or node.is_root():
            continue
        below = frozenset(n.name for n in node.get_leaves()) & all_taxa
        above = all_taxa - below
        if len(below) < 2 or len(above) < 2:
            continue
        small = below if len(below) <= len(above) else above
        support = float(node.support) if hasattr(node, "support") and node.support is not None else np.nan
        out.append((frozenset(small), support))
    return out


def max_jaccard(query: list[tuple[frozenset, float]],
                target: list[tuple[frozenset, float]],
                all_taxa: set[str]) -> tuple[float, int]:
    """Mean max-Jaccard of query bipartitions against target set, and n_exact."""
    if not query or not target:
        return 0.0, 0
    target_sets = [b for b, _ in target]
    target_lookup = set(target_sets)
    best_js = []
    n_exact = 0
    for q, _ in query:
        qc = all_taxa - q
        if q in target_lookup or qc in target_lookup:
            n_exact += 1
        best = 0.0
        for t in target_sets:
            for qq in (q, qc):
                inter = len(qq & t)
                uni = len(qq | t)
                j = inter / uni if uni else 0.0
                if j > best:
                    best = j
        best_js.append(best)
    return float(np.mean(best_js)), n_exact


def permutation_null(query_tree: ete3.Tree, target_biparts, all_taxa: set[str],
                     n_perm: int = 200, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    leaves = list(all_taxa)
    # For each permutation, remap leaf names of a copy of query_tree, recompute.
    null_vals = []
    for k in range(n_perm):
        perm = rng.permutation(leaves)
        name_map = dict(zip(leaves, perm))
        t_copy = query_tree.copy()
        for lf in t_copy.get_leaves():
            if lf.name in name_map:
                lf.name = name_map[lf.name]
        biparts = bipartitions(t_copy, all_taxa)
        j, _ = max_jaccard(biparts, target_biparts, all_taxa)
        null_vals.append(j)
    return {
        "n_perm": n_perm,
        "null_mean": float(np.mean(null_vals)),
        "null_p95": float(np.quantile(null_vals, 0.95)),
    }, null_vals


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--worm", required=True)
    ap.add_argument("--inference-dir", type=Path, required=True,
                    help="phylogeny/inference/<panel_tag>")
    ap.add_argument("--panel-h5ad", type=Path, required=True,
                    help="Used to pin the cell set (all_taxa).")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    a = ad.read_h5ad(args.panel_h5ad)
    all_taxa = set(a.obs_names.astype(str))

    cellphy_dir = args.inference_dir / "cellphy" / f"worm_{args.worm}"
    trees = {}
    for m in ("real", "covmask", "shuffled"):
        # Prefer support tree if bootstrap ran, else best tree.
        cand = [
            cellphy_dir / m / "sup.raxml.supportFBP",
            cellphy_dir / m / "run.raxml.support",
            cellphy_dir / m / "run.raxml.bestTree",
        ]
        for p in cand:
            t = _load_newick(p)
            if t is not None:
                trees[f"cellphy:{m}"] = (t, p)
                break

    print(f"[diagnostics] worm={args.worm}  loaded {len(trees)} trees:")
    for k, (_, p) in trees.items():
        print(f"    {k} <- {p}")

    if "cellphy:real" not in trees:
        print("[diagnostics] cellphy:real missing — nothing to compare against")
        return

    # bipartitions table
    rows = []
    biparts_by_src = {}
    for src, (tree, _) in trees.items():
        b = bipartitions(tree, all_taxa)
        biparts_by_src[src] = b
        for taxa, sup in b:
            rows.append({"worm": args.worm, "source": src, "small_size": len(taxa),
                         "support": sup, "taxa_small": ",".join(sorted(taxa))})
    pd.DataFrame(rows).to_csv(args.out_dir / "bipartitions.tsv", sep="\t", index=False)

    # overlap table: everything vs cellphy:real
    ref = biparts_by_src["cellphy:real"]
    over_rows = []
    for src, b in biparts_by_src.items():
        if src == "cellphy:real":
            continue
        mj, nex = max_jaccard(b, ref, all_taxa)
        over_rows.append({"worm": args.worm, "source": src, "n_biparts": len(b),
                          "n_exact_vs_real": nex, "mean_max_jaccard_vs_real": mj})
    pd.DataFrame(over_rows).to_csv(args.out_dir / "overlap.tsv", sep="\t", index=False)

    # permutation null: shuffle covmask leaves against real, and shuffled leaves against real
    null_rows = []
    for src in ("cellphy:covmask", "cellphy:shuffled"):
        if src not in trees:
            continue
        q_tree = trees[src][0]
        summary, null_vals = permutation_null(q_tree, ref, all_taxa, args.n_perm, args.seed)
        obs, _ = max_jaccard(biparts_by_src[src], ref, all_taxa)
        p_emp = float(np.mean(np.asarray(null_vals) >= obs)) if null_vals else np.nan
        null_rows.append({"worm": args.worm, "source": src, "observed_max_jaccard": obs,
                          "null_mean": summary["null_mean"], "null_p95": summary["null_p95"],
                          "empirical_p": p_emp, "n_perm": summary["n_perm"]})
    pd.DataFrame(null_rows).to_csv(args.out_dir / "permutation_null.tsv", sep="\t", index=False)
    print(f"[diagnostics] wrote {args.out_dir}/{{bipartitions,overlap,permutation_null}}.tsv")


if __name__ == "__main__":
    main()
