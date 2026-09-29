#!/usr/bin/env python3
"""Compute per-cell division mutation counts from CellPhy Mapped-tree output.

For each cell (leaf), walks the root→leaf path on the midpoint-rooted tree
and records:
  * n_muts_from_root: total mutations across all branches on the path
  * n_divisions_from_root: number of internal nodes traversed (= path depth)
  * n_muts_terminal: mutations on the terminal branch only
  * n_muts_internal: mutations on the internal branches only (from_root minus terminal)
  * mean_muts_per_division: n_muts_from_root / max(1, n_divisions_from_root)

Combined across all worms in the panel. Also writes the vectors back as obs
columns on the combined panel h5ad if it exists.

Usage
-----
  cellphy_percell_division_muts.py \
      --panel-tag all_variants_relax03__svdImpSym_mquad10_cs3 \
      --out-csv results/worm6_final/DNA_analysis/cellphy_percell_division_muts__<tag>.tsv \
      [--panel-h5ad results/worm6_final/DNA_analysis/panel__<tag>.h5ad]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
INF_ROOT = REPO / "results/worm6_final/DNA_analysis/phylogeny/inference"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel-tag", required=True)
    p.add_argument("--out-csv", type=Path, required=True)
    p.add_argument("--panel-h5ad", type=Path, default=None,
                   help="If given, also write per-cell columns back as obs.")
    return p.parse_args()


def _parse_mut_list(path: Path) -> dict[int, int]:
    """Return {branch_id: n_muts}."""
    counts = {}
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 2: continue
        try:
            bid = int(parts[0])
            n = int(parts[1])
        except ValueError:
            continue
        counts[bid] = n
    return counts


def _parse_tagged_tree(path: Path):
    """Return an ete3.Tree where each non-root node has ``.branch_id``
    from the ``[N]`` tag after its branch length in the mutationMapTree.
    """
    import ete3
    text = path.read_text()
    # Strip [N] tags before ete3 parses, keeping a map of leaf/label → branch_id.
    # Pattern: after ':<float>' comes '[<int>]'. We rewrite each ':X[N]' to ':X'
    # and remember the sequence order matches a post-order walk in raxml-ng
    # output — but that's fragile. Safer: parse via regex to associate the
    # NAME (or Node<M>) immediately preceding ':<float>[<id>]' with the id.
    branch_id: dict[str, int] = {}
    # matches:  <name>:<len>[<id>]   or  ):<len>[<id>]  or  )<label>:<len>[<id>]
    for m in re.finditer(r"([A-Za-z0-9_.]+)?:[0-9.eE+-]+\[(\d+)\]", text):
        name = m.group(1)
        bid  = int(m.group(2))
        if name is None:
            continue
        branch_id[name] = bid

    # Clean the tree of [N] and reparse
    clean = re.sub(r"\[\d+\]", "", text)
    t = ete3.Tree(clean, format=1)
    for node in t.traverse():
        if node.is_root(): continue
        nm = node.name
        if nm and nm in branch_id:
            node.branch_id = branch_id[nm]
        else:
            node.branch_id = None
    return t


def _midpoint_root(t):
    try:
        og = t.get_midpoint_outgroup()
        if og is not None:
            t.set_outgroup(og)
    except Exception:
        pass
    return t


def per_worm(worm_dir: Path) -> pd.DataFrame:
    real = worm_dir / "real"
    tree_p = real / "run.Mapped.raxml.mutationMapTree"
    list_p = real / "run.Mapped.raxml.mutationMapList"
    if not (tree_p.exists() and list_p.exists()):
        return pd.DataFrame()
    counts = _parse_mut_list(list_p)
    t = _parse_tagged_tree(tree_p)
    _midpoint_root(t)

    rows = []
    for leaf in t.get_leaves():
        # walk root → leaf, summing per-branch mutation counts
        path = leaf.get_ancestors()[::-1] + [leaf]  # root first, leaf last
        # only nodes with a valid branch_id contribute a branch
        n_muts_total = 0
        n_divs = 0
        per_branch = []
        for node in path:
            if node.is_root(): continue
            bid = getattr(node, "branch_id", None)
            n = counts.get(bid, 0) if bid is not None else 0
            per_branch.append(n)
            n_muts_total += n
            n_divs += 1
        terminal = per_branch[-1] if per_branch else 0
        internal = n_muts_total - terminal
        rows.append(dict(
            cell_id=leaf.name,
            n_muts_from_root=n_muts_total,
            n_divisions_from_root=n_divs,
            n_muts_terminal=terminal,
            n_muts_internal=internal,
            mean_muts_per_division=n_muts_total / max(1, n_divs),
        ))
    df = pd.DataFrame(rows)
    df["worm"] = worm_dir.name.replace("worm_", "")
    return df


def main():
    args = parse_args()
    panel_root = INF_ROOT / args.panel_tag / "cellphy"
    if not panel_root.exists():
        raise SystemExit(f"no cellphy dir under {panel_root}")

    dfs = []
    for wd in sorted(panel_root.glob("worm_*")):
        d = per_worm(wd)
        if not d.empty:
            print(f"  {wd.name}: {len(d)} cells")
            dfs.append(d)
    if not dfs:
        raise SystemExit("no per-worm data assembled")
    all_df = pd.concat(dfs, axis=0, ignore_index=True)
    all_df = all_df[["worm","cell_id","n_muts_from_root","n_divisions_from_root",
                     "n_muts_terminal","n_muts_internal","mean_muts_per_division"]]
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    all_df.to_csv(args.out_csv, sep="\t", index=False)
    print(f"[done] wrote {args.out_csv}  ({len(all_df)} cells)")

    # Optional: attach to combined panel h5ad as obs columns
    if args.panel_h5ad and args.panel_h5ad.exists():
        a = ad.read_h5ad(args.panel_h5ad)
        df_idx = all_df.set_index("cell_id")
        for col in ("n_muts_from_root","n_divisions_from_root",
                    "n_muts_terminal","n_muts_internal","mean_muts_per_division"):
            a.obs[col] = df_idx[col].reindex(a.obs_names).values
        a.write_h5ad(args.panel_h5ad, compression="gzip")
        n_match = int(a.obs["n_muts_from_root"].notna().sum())
        print(f"[done] wrote back to {args.panel_h5ad} "
              f"({n_match}/{a.n_obs} cells matched)")


if __name__ == "__main__":
    main()
