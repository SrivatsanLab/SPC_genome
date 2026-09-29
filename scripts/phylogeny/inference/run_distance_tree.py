#!/usr/bin/env python3
"""Distance-based tree inference (v2): NJ + BME × optional NNI × outgroup + MinVar rooting.

Distance metric (choose via --metric):

  hard   :  fraction of shared-covered sites where hard-carrier calls disagree
  soft   :  same shape, per-site disagreement uses P(alt|reads) from PL

Rooting protocols emitted (all off the same base topology per start-method):
  og     :  synthetic ROOT outgroup added to the distance matrix, tree
            rooted by that outgroup (K562-style)
  mad    :  MinVar rooting (Mai & Mirarab 2017) — pick root position to
            minimise variance of root-to-tip distances; requires no outgroup

Bootstrap: 100 variant-resample replicates. Each replicate uses the same
start-method (+ NNI if enabled). Bipartitions are computed over the ingroup
only (excluding ROOT so support values are consistent across rootings).

Outputs at <out_dir>/, per start-method:
  distance_matrix.tsv
  {metric}__nj__{og,mad}.support             — NJ starting tree
  {metric}__nj_nni__{og,mad}.support         — NJ + NNI polish
  {metric}__bme_nni__{og,mad}.support        — BME + NNI polish
  {metric}__upgma__{og,mad}.support          — UPGMA (rooted natively)
  {metric}__{...}.bootstrap.newick           — 100 replicate trees per method
"""
from __future__ import annotations

import argparse
import io
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import issparse
from skbio import DistanceMatrix, TreeNode
from skbio.tree import bme, nj, nni, upgma


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-h5ad", type=Path, required=True)
    ap.add_argument("--vcf", type=Path, default=None,
                    help="Per-worm VCF with PL. Required for --metric soft.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--metric", choices=("hard", "soft"), default="soft")
    ap.add_argument("--n-bootstraps", type=int, default=100)
    ap.add_argument("--min-shared-sites", type=int, default=30)
    ap.add_argument("--outgroup-label", default="ROOT")
    ap.add_argument("--no-outgroup", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    return ap.parse_args()


# ---------- distance construction (unchanged from v1) ------------------------

def _load_ad_dp(a: ad.AnnData):
    AD = a.layers["AD"].toarray() if issparse(a.layers["AD"]) else a.layers["AD"]
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else a.layers["DP"]
    return AD, DP


def _load_p_alt_from_vcf(vcf: Path, cells: list[str], var_ids: list[str]) -> np.ndarray:
    v = pysam.VariantFile(str(vcf))
    idx_cell = {c: i for i, c in enumerate(cells)}
    idx_var = {vid: i for i, vid in enumerate(var_ids)}
    sample_names = list(v.header.samples)
    P = np.full((len(cells), len(var_ids)), np.nan, dtype=np.float32)
    for rec in v:
        vid = f"{rec.chrom}-{rec.pos}-{rec.ref}>{rec.alts[0]}"
        j = idx_var.get(vid)
        if j is None: continue
        for s in sample_names:
            i = idx_cell.get(s)
            if i is None: continue
            pl = rec.samples[s].get("PL")
            if pl is None: continue
            p = np.array([10 ** (-x / 10) for x in pl], dtype=np.float64)
            tot = p.sum()
            if tot > 0:
                P[i, j] = float((p[1] + p[2]) / tot)
    v.close()
    return P


def _hard_distance(carrier, covered, min_shared):
    n = carrier.shape[0]
    D = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            m = covered[i] & covered[j]
            nz = m.sum()
            if nz < min_shared:
                D[i, j] = D[j, i] = np.nan; continue
            d = np.logical_xor(carrier[i, m], carrier[j, m]).sum() / nz
            D[i, j] = D[j, i] = d
    return D


def _soft_distance(p_alt, covered, min_shared):
    n = p_alt.shape[0]
    D = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        pi = p_alt[i]
        for j in range(i + 1, n):
            m = covered[i] & covered[j]
            nz = m.sum()
            if nz < min_shared:
                D[i, j] = D[j, i] = np.nan; continue
            D[i, j] = D[j, i] = float(np.abs(pi[m] - p_alt[j, m]).mean())
    return D


def _fill_nan_with_median(D):
    if not np.any(np.isnan(D)): return D
    off = D[~np.eye(D.shape[0], dtype=bool)]
    med = float(np.nanmedian(off))
    D = D.copy(); D[np.isnan(D)] = med
    return D


def _root_distance_from_signal(signal, covered, min_shared):
    n = signal.shape[0]
    out = np.zeros(n, dtype=np.float32)
    for i in range(n):
        m = covered[i]
        out[i] = float(signal[i, m].mean()) if m.sum() >= min_shared else np.nan
    return out


def _append_root(D, cells, root_dists, label):
    n = D.shape[0]
    D2 = np.zeros((n + 1, n + 1), dtype=np.float32)
    D2[:n, :n] = D
    D2[:n, n] = root_dists; D2[n, :n] = root_dists
    return D2, cells + [label]


# ---------- tree construction ------------------------------------------------

def _start_tree(D, ids, method):
    dm = DistanceMatrix(D, ids=ids)
    if method == "nj": return nj(dm)
    if method == "bme": return bme(dm)
    if method == "upgma": return upgma(dm)
    raise ValueError(method)


def _polish_nni(tree, D, ids, max_iter=10):
    """Apply skbio NNI to improve the tree. Silently skips on failure."""
    try:
        dm = DistanceMatrix(D, ids=ids)
        return nni(tree, dm)
    except Exception as e:
        print(f"  NNI failed ({e}); returning unpolished tree.")
        return tree


def _root_by_outgroup(tree, outgroup):
    try:
        return tree.root_by_outgroup([outgroup])
    except Exception:
        return tree


# ---------- MinVar rooting (Mai & Mirarab 2017) ------------------------------

def _minvar_root(tree_ete, exclude_names: set | None = None):
    """MinVar rooting (snap-to-node): pick the node that, when used as the
    outgroup, minimises variance of root-to-tip distances. Robust to tree
    topology quirks (unlike edge-insertion approaches that can produce
    malformed newicks after ete3 rerooting).

    ``exclude_names`` skips leaves with those names (e.g., the synthetic
    outgroup ROOT — picking it would trivially give outgroup rooting).
    """
    import ete3
    from copy import deepcopy

    exclude_names = exclude_names or set()
    best_var = float("inf")
    best_node = None
    leaves = tree_ete.get_leaves()

    for candidate in list(tree_ete.traverse()):
        if candidate.is_root():
            continue
        if candidate.is_leaf() and candidate.name in exclude_names:
            continue
        # Reroot a *copy* to avoid destabilising the traversal.
        t_copy = deepcopy(tree_ete)
        # Locate the corresponding node in the copy by path from root.
        # Simplest: find by name if leaf, else by leaves-under signature.
        if candidate.is_leaf():
            match = next((lf for lf in t_copy.get_leaves()
                          if lf.name == candidate.name), None)
        else:
            target_set = frozenset(l.name for l in candidate.get_leaves())
            match = None
            for n in t_copy.traverse():
                if frozenset(l.name for l in n.get_leaves()) == target_set:
                    match = n; break
        if match is None:
            continue
        try:
            t_copy.set_outgroup(match)
        except Exception:
            continue
        # Variance of root-to-tip distances.
        dists = np.array([t_copy.get_distance(lf) for lf in t_copy.get_leaves()])
        v = float(np.var(dists))
        if v < best_var:
            best_var = v
            best_node = candidate

    if best_node is None:
        return tree_ete

    # Apply the best rooting to the actual tree.
    try:
        if best_node.is_leaf():
            target = next((lf for lf in tree_ete.get_leaves()
                           if lf.name == best_node.name), None)
        else:
            leaf_sig = frozenset(l.name for l in best_node.get_leaves())
            target = next((n for n in tree_ete.traverse()
                           if frozenset(l.name for l in n.get_leaves()) == leaf_sig), None)
        if target is not None:
            tree_ete.set_outgroup(target)
    except Exception:
        pass
    return tree_ete


def _apply_mad_rooting(nw: str, drop_outgroup: str | None = None) -> str:
    """Reroot a newick string via MinVar (snap-to-node).

    If ``drop_outgroup`` is set, that leaf is excluded from the candidate
    outgroup set (we don't want to pick the synthetic root as the MAD root —
    that would collapse MAD to the outgroup rooting). ROOT is removed from
    the final tree via pruning after rerooting.
    """
    import ete3
    t = ete3.Tree(nw, format=0)
    t = _minvar_root(t, exclude_names=({drop_outgroup} if drop_outgroup else set()))
    if drop_outgroup is not None:
        # Prune ROOT: keep every leaf except ROOT.
        keep = [lf.name for lf in t.get_leaves() if lf.name != drop_outgroup]
        if keep:
            t.prune(keep, preserve_branch_length=True)
    return t.write(format=0)


# ---------- bootstrap + support ---------------------------------------------

def _tree_to_newick(tree: TreeNode) -> str:
    buf = io.StringIO()
    tree.write(buf)
    return buf.getvalue().strip()


def _bipartitions(tree: TreeNode, all_leaves: frozenset[str]) -> set[frozenset[str]]:
    out = set()
    for node in tree.non_tips(include_self=False):
        below = frozenset(l.name for l in node.tips())
        other = all_leaves - below
        small = below if len(below) <= len(other) else other
        if len(small) >= 2:
            out.add(small)
    return out


def _annotate_support(base_tree: TreeNode, boot_trees: list[str]) -> TreeNode:
    all_leaves = frozenset(l.name for l in base_tree.tips())
    boot_bipsets = []
    for nw in boot_trees:
        try:
            t = TreeNode.read(io.StringIO(nw), format="newick")
        except Exception:
            continue
        boot_bipsets.append(_bipartitions(t, all_leaves))
    n = len(boot_bipsets)
    if n == 0:
        return base_tree
    for node in base_tree.non_tips(include_self=False):
        below = frozenset(l.name for l in node.tips())
        other = all_leaves - below
        small = below if len(below) <= len(other) else other
        if len(small) < 2:
            continue
        hits = sum(1 for s in boot_bipsets if small in s)
        node.name = f"{100 * hits / n:.0f}"
    return base_tree


# ---------- driver -----------------------------------------------------------

def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    a = ad.read_h5ad(args.panel_h5ad)
    cells = list(a.obs_names.astype(str))
    var_ids = list(a.var_names.astype(str))
    AD, DP = _load_ad_dp(a)
    covered = DP > 0
    carrier = AD >= 1
    n_c, n_v = AD.shape

    if args.metric == "soft":
        if args.vcf is None:
            raise SystemExit("--metric soft requires --vcf")
        print(f"[dist] loading PL from {args.vcf}")
        P = _load_p_alt_from_vcf(args.vcf, cells, var_ids)
        n_miss = int(np.isnan(P).sum())
        if n_miss:
            print(f"[dist] PL missing at {n_miss}/{P.size}; falling back to carrier")
            P = np.where(np.isnan(P), carrier.astype(np.float32), P)

    def pairwise(idx):
        if args.metric == "hard":
            return _hard_distance(carrier[:, idx], covered[:, idx], args.min_shared_sites)
        return _soft_distance(P[:, idx], covered[:, idx], args.min_shared_sites)

    def root_dists(idx):
        sig = carrier[:, idx].astype(np.float32) if args.metric == "hard" else P[:, idx]
        return _root_distance_from_signal(sig, covered[:, idx], args.min_shared_sites)

    def build_dm(idx):
        D = pairwise(idx)
        ids = list(cells)
        og = None
        if not args.no_outgroup:
            rd = root_dists(idx)
            if np.any(np.isnan(rd)):
                rd = np.where(np.isnan(rd), float(np.nanmedian(rd[~np.isnan(rd)])) if np.any(~np.isnan(rd)) else 0.0, rd)
            D, ids = _append_root(D, cells, rd, args.outgroup_label)
            og = args.outgroup_label
        return _fill_nan_with_median(D), ids, og

    print(f"[dist] worm={a.uns.get('panel_worm', '?')}  {n_c} cells × {n_v} vars  "
          f"metric={args.metric}")
    D_base, ids_base, og_base = build_dm(np.arange(n_v))
    pd.DataFrame(D_base, index=ids_base, columns=ids_base).to_csv(
        args.out_dir / "distance_matrix.tsv", sep="\t")

    # Method pipelines to run: (name, start, apply_nni)
    pipelines = [
        ("nj",       "nj",    False),
        ("nj_nni",   "nj",    True),
        ("bme_nni",  "bme",   True),
        ("upgma",    "upgma", False),
    ]

    # Build base tree per pipeline
    base_by_pipe: dict[str, TreeNode] = {}
    for pname, start, do_nni in pipelines:
        print(f"[dist] {pname}: build ...")
        t = _start_tree(D_base, ids_base, start)
        if do_nni:
            t = _polish_nni(t, D_base, ids_base)
        base_by_pipe[pname] = t

    # Bootstrap for support (bipartitions only; rooting-invariant)
    per_pipe_boots: dict[str, list[str]] = {pname: [] for pname, _, _ in pipelines}
    if args.n_bootstraps > 0:
        print(f"[dist] running {args.n_bootstraps} bootstraps ...")
        for k in range(args.n_bootstraps):
            idx = rng.integers(0, n_v, size=n_v)
            D_k, ids_k, _ = build_dm(idx)
            for pname, start, do_nni in pipelines:
                try:
                    tk = _start_tree(D_k, ids_k, start)
                    if do_nni:
                        tk = _polish_nni(tk, D_k, ids_k)
                    per_pipe_boots[pname].append(_tree_to_newick(tk))
                except Exception as e:
                    print(f"  boot {k} {pname} failed: {e}")

    # For each pipeline: annotate support, then emit both rootings.
    for pname, _, _ in pipelines:
        base = base_by_pipe[pname]
        boots = per_pipe_boots[pname]
        # Save raw bootstrap newicks
        with open(args.out_dir / f"{args.metric}__{pname}.bootstrap.newick", "w") as f:
            f.write("\n".join(boots))
        # Annotate the base tree with FBP support.
        sup_tree = _annotate_support(base, boots)
        # Outgroup-rooted version
        if og_base is not None:
            sup_og = _root_by_outgroup(sup_tree, og_base)
            with open(args.out_dir / f"{args.metric}__{pname}__og.support", "w") as f:
                f.write(_tree_to_newick(sup_og))
        # MinVar (MAD-family) rooting; drop the outgroup first if present so
        # MinVar isn't fooled by ROOT sitting near the "ancestral" side.
        try:
            nw_for_mad = _tree_to_newick(sup_tree)
            nw_mad = _apply_mad_rooting(nw_for_mad, drop_outgroup=og_base)
            with open(args.out_dir / f"{args.metric}__{pname}__mad.support", "w") as f:
                f.write(nw_mad)
        except Exception as e:
            print(f"  MAD rooting failed on {pname}: {e}")

    print("[dist] done")


if __name__ == "__main__":
    main()
