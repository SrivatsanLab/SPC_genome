#!/usr/bin/env python3
"""Convert top-down SVD ``tree.json`` bipartition trees to Newick.

The recursive bipartition tree in ``tree.json`` has no branch lengths — depth is
an integer recursion level — and its per-split quality is measured by
``module_purity_median`` (median softmax margin from K-means on the filtered
core; 0-1 scale). To make these trees usable in downstream Newick-based
plotters (the notebook's ``ete3.Tree(..., format=0)`` + ``node.support``),
emit each internal node with:

  * ``support`` ← mean(purity_L, purity_R) × 100  (0-100 scale, matches the
    notebook's bootstrap support colorbar).
  * branch length ← 1.0  (topology only; the notebook's ``max_depth``
    normalizer handles depth scaling).
  * leaf names ← cell IDs from the tree.json.

Usage (single tree):
    topdown_tree_to_newick.py --tree-json path/to/tree.json --out path/to/tree.support.newick

Usage (walk a whole panel dir, emit .support.newick next to each tree.json):
    topdown_tree_to_newick.py --panel-tag all_variants_relax03
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
TD_ROOT = REPO / "results/worm6_final/DNA_analysis/phylogeny/topdown"

_UNSAFE = re.compile(r"[\s(),:;\[\]']")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--tree-json", type=Path)
    grp.add_argument("--panel-tag")
    ap.add_argument("--out", type=Path, default=None,
                    help="Only valid with --tree-json. Default: <tree.json>.support.newick")
    ap.add_argument("--td-root", type=Path, default=TD_ROOT,
                    help="Root of the top-down sweep output tree.")
    ap.add_argument("--branch-length", type=float, default=1.0,
                    help="Uniform branch length for every edge (default 1.0).")
    return ap.parse_args()


def _safe(name: str) -> str:
    """Return a Newick-safe leaf label — strip whitespace / structural chars."""
    return _UNSAFE.sub("_", str(name))


# Sentinel support value for polytomies (multi-cell leaves where recursion
# couldn't split further). Falls below the notebook's viridis vmin=0 so it
# renders as dark = "unresolved" rather than being confused with a low-support
# real split.
POLYTOMY_SUPPORT = -1.0


def _to_newick(node: dict, br_len: float) -> str:
    """Recursively serialize a bipartition tree node to a Newick clade string.

    Returns the clade *without* the trailing branch length; the caller
    (parent) appends ``:<br_len>``.
    """
    if node.get("split") is None:
        cells = node.get("cells", [])
        if len(cells) == 1:
            return _safe(cells[0])
        # Multi-cell leaf: recursion stopped by size/depth/no_split. Render as
        # an unresolved star clade with an explicit sentinel support so ete3
        # doesn't default to 1.0 and pollute the support color scale.
        inner = ",".join(f"{_safe(c)}:{br_len:g}" for c in cells)
        return f"({inner}){POLYTOMY_SUPPORT:g}"
    left  = _to_newick(node["left"],  br_len)
    right = _to_newick(node["right"], br_len)
    sp = node["split"]
    pur_L = float(sp["module_purity_median"]["L"])
    pur_R = float(sp["module_purity_median"]["R"])
    support = 100.0 * (pur_L + pur_R) / 2.0
    # ete3 format 0 reads an integer/float token after the closing paren as node.support
    return f"({left}:{br_len:g},{right}:{br_len:g}){support:.1f}"


def tree_json_to_newick(tree_json: dict, branch_length: float = 1.0) -> str:
    root = tree_json["root"]
    body = _to_newick(root, branch_length)
    # If root itself was a leaf (only 1 cell in panel — pathological), still emit valid newick
    if not body.startswith("("):
        return f"{body}:{branch_length:g};"
    return f"{body};"


def convert_one(tree_json_path: Path, out_path: Path, branch_length: float = 1.0) -> None:
    d = json.loads(tree_json_path.read_text())
    nwk = tree_json_to_newick(d, branch_length=branch_length)
    out_path.write_text(nwk + "\n")


def main() -> None:
    args = parse_args()
    if args.tree_json is not None:
        out = args.out or args.tree_json.with_name(args.tree_json.stem + ".support.newick")
        convert_one(args.tree_json, out, args.branch_length)
        print(f"[to_newick] {args.tree_json}  ->  {out}")
        return

    # panel-tag walk
    panel_dir = args.td_root / args.panel_tag
    if not panel_dir.exists():
        raise SystemExit(f"panel dir not found: {panel_dir}")
    tjs = sorted(panel_dir.glob("*/worm_*/tree.json"))
    print(f"[to_newick] {len(tjs)} tree.json files under {panel_dir}")
    for tj in tjs:
        out = tj.with_name("tree.support.newick")
        try:
            convert_one(tj, out, args.branch_length)
        except Exception as e:
            print(f"  [err] {tj}: {e}")
    print(f"[to_newick] wrote tree.support.newick next to each tree.json")


if __name__ == "__main__":
    main()
