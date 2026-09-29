#!/usr/bin/env python3
"""Tree-guided haplotype imputation from an SVD mquad tree.

Takes a per-worm panel h5ad and its corresponding SVD-mquad tree.json,
imputes AD/DP layers per the plan (docs/tree_guided_imputation_plan.md).

Carrier-side rule (always on when --pseudo-alt > 0):
  Modify (cell, variant) entries where:
    - the variant is in module_variants[m] at some split S (mquad-selected)
    - the cell's ancestor path passes through S on side m
    - the observation doesn't already confirm carrier status (AD_v < 1)
  New AD = old AD + pseudo_alt
  New DP = old DP + pseudo_alt

Non-carrier-side rule (on when --pseudo-ref > 0):
  For the same module variants, cells whose ancestor path took the
  *opposite* side of split S are predicted non-carriers. For those cells,
  where the observation is entirely absent (DP_v == 0), inject a pseudo-ref
  observation:
    New AD = old AD (unchanged, still 0)
    New DP = old DP + pseudo_ref
  Skipped if DP_v >= 1 (cell already has real observation of any state).

Preserves the original AD/DP as ``AD_raw``/``DP_raw`` layers.

Usage
-----
  impute_panel.py \
      --panel-h5ad phylogeny/panels/<tag>/worm_<W>.h5ad \
      --tree-json  phylogeny/inference/<tag>/svd_mquad5/worm_<W>/tree.json \
      --out-h5ad   phylogeny/panels/<tag>__svdImp/worm_<W>.h5ad \
      [--pseudo-alt 1] [--pseudo-ref 0] [--trust-max-depth <int>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-h5ad", type=Path, required=True)
    ap.add_argument("--tree-json", type=Path, required=True)
    ap.add_argument("--out-h5ad", type=Path, required=True)
    ap.add_argument("--pseudo-alt", type=int, default=1,
                    help="Pseudo-alt read(s) to inject at tree-predicted carrier positions "
                         "where AD < 1. Set to 0 to disable carrier-side imputation.")
    ap.add_argument("--pseudo-ref", type=int, default=0,
                    help="Pseudo-ref read(s) to inject at tree-predicted non-carrier positions "
                         "where DP == 0 (i.e., uncovered). Adds to DP only, leaves AD unchanged. "
                         "Default 0 = disabled (asymmetric imputation). Set to 1 for symmetric.")
    ap.add_argument("--trust-max-depth", type=int, default=None,
                    help="Skip splits deeper than this. None (default) trusts every split in the tree.")
    return ap.parse_args()


def build_cell_ancestor_paths(root: dict) -> dict[str, dict[str, str]]:
    """For each leaf cell, return {split_path_key: side_taken_at_that_split}.

    ``split_path_key`` is the path of the parent internal node — the split
    the cell passed through on its way down. ``side_taken`` is 'L' or 'R'.
    """
    cell_paths: dict[str, dict[str, str]] = {}

    def walk(node, current_path: dict[str, str]):
        if "split" in node and node["split"] is not None:
            # This node has a split; recurse into left and right children.
            parent_path = node.get("path", "") or "root"
            for side, child_key in (("L", "left"), ("R", "right")):
                if child_key not in node:
                    continue
                child = node[child_key]
                new_path = dict(current_path)
                new_path[parent_path] = side
                walk(child, new_path)
        else:
            # Leaf — record ancestor path for every cell in this clade.
            for c in node.get("cells", []):
                cell_paths[str(c)] = dict(current_path)

    walk(root, {})
    return cell_paths


def build_marker_pool(root: dict, trust_max_depth: int | None) -> pd.DataFrame:
    """Walk the tree, extract every (variant, split_path, module, depth).

    Returns a DataFrame with columns: variant_id, split_path, module, depth.
    """
    rows = []

    def walk(node):
        if "split" not in node or node["split"] is None:
            return
        depth = node.get("depth", 0)
        if trust_max_depth is not None and depth > trust_max_depth:
            return
        split_path = node.get("path", "") or "root"
        s = node["split"]
        for module in ("L", "R"):
            for v in s["module_variants"][module]:
                rows.append({"variant_id": v, "split_path": split_path,
                             "module": module, "depth": depth})
        if "left" in node and "right" in node:
            walk(node["left"]); walk(node["right"])

    walk(root)
    # Explicit columns so callers can access .variant_id etc. even when the
    # tree produced no splits (empty markers is a valid no-op imputation).
    return pd.DataFrame(rows, columns=["variant_id", "split_path", "module", "depth"])


def main() -> None:
    args = parse_args()
    args.out_h5ad.parent.mkdir(parents=True, exist_ok=True)

    print(f"[impute] loading panel {args.panel_h5ad}")
    a = ad.read_h5ad(args.panel_h5ad)
    AD = a.layers["AD"].toarray() if issparse(a.layers["AD"]) else a.layers["AD"].copy()
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else a.layers["DP"].copy()
    cell_names = list(a.obs_names.astype(str))
    var_names = list(a.var_names.astype(str))
    cell_idx = {c: i for i, c in enumerate(cell_names)}
    var_idx = {v: i for i, v in enumerate(var_names)}

    print(f"[impute] loading tree {args.tree_json}")
    with open(args.tree_json) as f:
        tree = json.load(f)

    ancestor_paths = build_cell_ancestor_paths(tree["root"])
    print(f"[impute] {len(ancestor_paths)}/{len(cell_names)} cells have an ancestor path in the tree")

    markers = build_marker_pool(tree["root"], args.trust_max_depth)
    if markers.empty:
        print(f"[impute] marker pool: 0 rows — tree has no valid splits; writing pass-through h5ad")
    else:
        print(f"[impute] marker pool: {len(markers):,} rows, "
              f"{markers.variant_id.nunique():,} unique variants, "
              f"{markers.split_path.nunique()} splits, "
              f"depths {markers.depth.min()}..{markers.depth.max()}")

    # For each (cell, variant) tuple: check if tree predicts carrier.
    # Use an efficient loop: for each marker row, iterate over cells whose
    # ancestor path passes through this split on the same side.
    manifest_rows = []
    n_imputed_carrier = 0
    n_imputed_noncarrier = 0
    for (split_path, module), sub in markers.groupby(["split_path", "module"]):
        v_ids = sub.variant_id.to_numpy()
        v_indices = np.array([var_idx[v] for v in v_ids if v in var_idx], dtype=np.int64)
        if v_indices.size == 0:
            continue

        # --- Carrier-side: cells whose ancestor path took `module` at this split.
        if args.pseudo_alt > 0:
            carrier_cells = [c for c, path in ancestor_paths.items()
                             if path.get(split_path) == module]
            if carrier_cells:
                c_indices = np.array([cell_idx[c] for c in carrier_cells], dtype=np.int64)
                sub_AD = AD[np.ix_(c_indices, v_indices)]
                need_impute = sub_AD < 1
                if need_impute.any():
                    for i_local, c_i in enumerate(c_indices):
                        for j_local, v_i in enumerate(v_indices):
                            if need_impute[i_local, j_local]:
                                AD[c_i, v_i] += args.pseudo_alt
                                DP[c_i, v_i] += args.pseudo_alt
                                n_imputed_carrier += 1
                                manifest_rows.append({
                                    "cell_id": cell_names[c_i],
                                    "variant_id": var_names[v_i],
                                    "split_path": split_path,
                                    "module": module,
                                    "depth": int(sub.iloc[0].depth),
                                    "imputation_type": "carrier",
                                })

        # --- Non-carrier-side: cells whose ancestor path took the OPPOSITE side.
        # Only fires for uncovered cells (DP == 0); real evidence stays untouched.
        if args.pseudo_ref > 0:
            opposite = "L" if module == "R" else "R"
            noncar_cells = [c for c, path in ancestor_paths.items()
                            if path.get(split_path) == opposite]
            if noncar_cells:
                nc_indices = np.array([cell_idx[c] for c in noncar_cells], dtype=np.int64)
                sub_DP_nc = DP[np.ix_(nc_indices, v_indices)]
                need_ref = sub_DP_nc == 0
                if need_ref.any():
                    for i_local, c_i in enumerate(nc_indices):
                        for j_local, v_i in enumerate(v_indices):
                            if need_ref[i_local, j_local]:
                                # Add ref evidence: DP += pseudo_ref, AD unchanged.
                                DP[c_i, v_i] += args.pseudo_ref
                                n_imputed_noncarrier += 1
                                manifest_rows.append({
                                    "cell_id": cell_names[c_i],
                                    "variant_id": var_names[v_i],
                                    "split_path": split_path,
                                    "module": module,
                                    "depth": int(sub.iloc[0].depth),
                                    "imputation_type": "non_carrier",
                                })

    n_imputed = n_imputed_carrier + n_imputed_noncarrier
    print(f"[impute] imputed {n_imputed:,} (cell, variant) entries "
          f"({n_imputed_carrier:,} carrier + {n_imputed_noncarrier:,} non-carrier)")

    # Build new AnnData with the original AD/DP preserved.
    orig_AD_layer = a.layers["AD"]
    orig_DP_layer = a.layers["DP"]
    new_a = a.copy()
    new_a.layers["AD"] = csr_matrix(AD) if issparse(orig_AD_layer) else AD
    new_a.layers["DP"] = csr_matrix(DP) if issparse(orig_DP_layer) else DP
    new_a.layers["AD_raw"] = orig_AD_layer
    new_a.layers["DP_raw"] = orig_DP_layer
    new_a.uns["imputation"] = {
        "source_tree": str(args.tree_json),
        "source_panel": str(args.panel_h5ad),
        "pseudo_alt": int(args.pseudo_alt),
        "pseudo_ref": int(args.pseudo_ref),
        "trust_max_depth": int(args.trust_max_depth) if args.trust_max_depth is not None else -1,
        "n_imputed_entries": int(n_imputed),
        "n_imputed_carrier": int(n_imputed_carrier),
        "n_imputed_noncarrier": int(n_imputed_noncarrier),
    }
    new_a.write_h5ad(args.out_h5ad)
    print(f"[impute] wrote {args.out_h5ad}")

    # Sidecar TSVs
    outdir = args.out_h5ad.parent
    manifest = pd.DataFrame(manifest_rows)
    if manifest.empty:
        manifest = pd.DataFrame(columns=["cell_id","variant_id","split_path","module","depth","imputation_type"])
    manifest.to_csv(outdir / f"{args.out_h5ad.stem}.imputation_manifest.tsv",
                    sep="\t", index=False)
    markers.to_csv(outdir / f"{args.out_h5ad.stem}.mquad_marker_pool.tsv",
                   sep="\t", index=False)
    print(f"[impute] wrote {args.out_h5ad.stem}.imputation_manifest.tsv "
          f"and .mquad_marker_pool.tsv")


if __name__ == "__main__":
    main()
