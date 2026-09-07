#!/usr/bin/env python3
"""Scanpy preprocessing for CapGTA spliced count h5ads.

Pipeline:
    load h5ad
    rename var_names WBID → gene symbol (WBID kept in .var['wormbase_id'])
    [if --drop-cuticle: drop col-* + dpy-7/sqt-1/sqt-3/rol-6/bli-1/bli-2]
    calculate_qc_metrics
    save raw integer counts to .layers['counts']   (unconditional; never overwritten)
    populate .var['n_junctions'], .var['exonic_length'], .var['detectability_k']
        (when the corresponding CSVs are provided; diagnostics-only in 'none' mode)
    apply --count-scaling arm to .X
    write .obs[assay_key] = assay_value
    normalize_total → log1p → HVG → PCA → neighbors → UMAP → Leiden
    write h5ad

Count-scaling arms (see docs/C3.3_count_scaling_review.md):
    none         .X = raw ints                                     (arm A, default)
    current      .X @ diag(1 / (n_junctions · exonic_length_kb))   (arm B, legacy)
    median_rint  X * median(k)/k, then np.rint                     (arm C, tutorial)
    median_rrint X * median(k)/k, then randomized rounding         (arm D)
Where k_g = n_junctions · exonic_length / 1000. Genes with k_g == 0 collapse to 0
in the scaling arms; in 'none' they are kept intact.

Legacy compatibility: if --count-scaling is not passed and either --junction-csv
or --gene-lengths-csv IS passed, arm defaults to 'current' so the existing SoupX
pipeline (which passes both CSVs) continues to reproduce arm B without an
orchestrator change.

.layers['counts'] always holds the raw integer counts as-loaded, which is what
SoupX and downstream count-model tools (scVI) expect.

Usage:
    preprocess_scanpy.py <input.h5ad> <output.h5ad> --gtf <annotation.gtf> \
        [--drop-cuticle] \
        [--junction-csv <gene_junctions.csv>] \
        [--gene-lengths-csv <gene_lengths.csv>] \
        [--count-scaling {none,current,median_rint,median_rrint}] \
        [--seed 0] \
        [--assay-key assay --assay-value capgta] \
        [--resolution 2] [--n-hvg 2000] [--batch-key sample] \
        [--min-dist 0.1] [--spread 5]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp


def gtf_gid_to_symbol(gtf_path: Path) -> dict[str, str]:
    """Return gene_id → gene_name map from a GTF (first-seen wins)."""
    m: dict[str, str] = {}
    with gtf_path.open() as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            gid = sym = None
            for kv in fields[8].strip().rstrip(';').split(';'):
                kv = kv.strip()
                if kv.startswith('gene_id '):
                    gid = kv.split('"')[1]
                elif kv.startswith('gene_name '):
                    sym = kv.split('"')[1]
            if gid and sym and gid not in m:
                m[gid] = sym
    return m


CUTICLE_PREFIXES = (
    # cuticle-structural (collagens, cuticulins)
    'col-', 'cut-', 'idpa-',
    # molt-oscillator gene families (Meeuse 2020 Mol Syst Biol; ~24% of expressed
    # genes oscillate with a ~7h period; cuticle collagens + regulators peak at
    # phased offsets around each molt)
    'bli-', 'ptr-', 'mlt-', 'noah-',
    # all dpy-* — drop broadly since the molt oscillator drives them regardless
    # of whether the specific family member is cuticle-structural or regulatory
    'dpy-',
)
CUTICLE_SYMBOLS = {
    # squat / roller cuticle families (bli-* now covered by prefix above)
    'sqt-1', 'sqt-2', 'sqt-3',
    'rol-1', 'rol-3', 'rol-4', 'rol-5', 'rol-6', 'rol-8',
    # additional molt-oscillator markers cited in Meeuse 2020 and WormBook cuticle
    'fbn-1', 'qua-1',
}


def cuticle_exclude_symbols(gid2sym: dict[str, str]) -> set[str]:
    """Genes excluded by the --drop-cuticle flag: col-/cut-/idpa- structural cuticle
    families + molt-oscillator families (bli-, ptr-, mlt-, noah-, dpy-) + curated
    sqt/rol/fbn-1/qua-1. See docs/C3.3_count_scaling_review.md §F.2 for the
    Meeuse 2020 justification — the molt oscillator drives ~24% of the transcriptome
    with a ~7h period and can look like lineage/subtype structure if left in."""
    syms = {s for s in gid2sym.values() if s.lower().startswith(CUTICLE_PREFIXES)}
    syms.update(CUTICLE_SYMBOLS)
    return syms


def _apply_scaling(X, k, arm: str, rng: np.random.Generator):
    """Return a new .X scaled per arm. X is CSR; k is the per-gene detectability
    vector; arm ∈ {none, current, median_rint, median_rrint}."""
    if arm == 'none':
        return X
    X = X if sp.issparse(X) else sp.csr_matrix(X)
    has_k = k > 0
    if arm == 'current':
        scale = np.where(has_k, 1.0 / np.where(has_k, k, 1.0), 0.0)
        return (X @ sp.diags(scale)).tocsr()
    if arm in ('median_rint', 'median_rrint'):
        k_med = float(np.median(k[has_k])) if has_k.any() else 1.0
        scale = np.where(has_k, k_med / np.where(has_k, k, 1.0), 0.0)
        Xs = (X @ sp.diags(scale)).tocsr()
        data = Xs.data
        if arm == 'median_rint':
            data = np.rint(data)
        else:
            floor = np.floor(data)
            frac = data - floor
            data = floor + (rng.random(size=frac.shape) < frac).astype(floor.dtype)
        Xs.data = data.astype(np.int64, copy=False)
        Xs.eliminate_zeros()
        return Xs
    raise ValueError(f'unknown count-scaling arm: {arm}')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input_h5ad', type=Path)
    p.add_argument('output_h5ad', type=Path)
    p.add_argument('--gtf', type=Path, required=True,
                   help='GTF used to map gene_id → gene_name for var rename')
    p.add_argument('--drop-cuticle', action='store_true',
                   help='drop col-* collagens plus dpy-7, sqt-1, sqt-3, rol-6, bli-1, bli-2')
    p.add_argument('--junction-csv', type=Path, default=None,
                   help='gene_junctions.csv (columns: gene_id,n_junctions); always populates .var[n_junctions]')
    p.add_argument('--gene-lengths-csv', type=Path, default=None,
                   help='gene_lengths.csv (columns: gene_id,exonic_length); always populates .var[exonic_length]')
    p.add_argument('--count-scaling', choices=['none', 'current', 'median_rint', 'median_rrint'],
                   default=None,
                   help="how to scale .X: none (raw ints, default when no CSVs), current (legacy 1/(J·L)), "
                        "median_rint (X*median(k)/k → rint), median_rrint (X*median(k)/k → randomized rounding). "
                        "If unset and either CSV is passed, defaults to 'current' for SoupX-pipeline compatibility.")
    p.add_argument('--seed', type=int, default=0,
                   help='RNG seed for median_rrint randomized rounding (default 0)')
    p.add_argument('--assay-key', default='assay',
                   help='name of the .obs column written to identify this assay (default: assay)')
    p.add_argument('--assay-value', default='capgta',
                   help='value written to .obs[assay_key] for every cell (default: capgta)')
    p.add_argument('--resolution', type=float, default=2.0)
    p.add_argument('--n-hvg', type=int, default=2000)
    p.add_argument('--batch-key', default='sample')
    p.add_argument('--min-dist', type=float, default=0.1)
    p.add_argument('--spread', type=float, default=5.0)
    args = p.parse_args()

    if not args.input_h5ad.is_file():
        print(f'Error: input not found: {args.input_h5ad}', file=sys.stderr)
        return 1
    if not args.gtf.is_file():
        print(f'Error: GTF not found: {args.gtf}', file=sys.stderr)
        return 1

    # Resolve count-scaling arm with legacy fallback.
    has_csvs = args.junction_csv is not None or args.gene_lengths_csv is not None
    if args.count_scaling is None:
        arm = 'current' if has_csvs else 'none'
    else:
        arm = args.count_scaling
    if arm != 'none' and not has_csvs:
        print(f"Error: --count-scaling={arm} requires --junction-csv and --gene-lengths-csv", file=sys.stderr)
        return 1

    print(f'Loading {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)
    print(f'  {adata.n_obs} cells x {adata.n_vars} genes')

    # --- Rename var_names WBID → gene symbol ---------------------------------
    # Prefer .var['wormbase_id'] when present (e.g., post-SoupX reimport already
    # holds the original WBIDs; var_names may already be symbols).
    gid2sym = gtf_gid_to_symbol(args.gtf)
    if 'wormbase_id' in adata.var.columns:
        wbids = adata.var['wormbase_id'].astype(str).to_numpy()
    else:
        wbids = adata.var_names.to_numpy()
    new_names = np.array([gid2sym.get(g, g) for g in wbids])
    n_renamed = int(sum(1 for g in wbids if g in gid2sym))
    adata.var['wormbase_id'] = wbids
    adata.var_names = new_names
    adata.var_names_make_unique()
    print(f'Renamed {n_renamed}/{len(wbids)} var_names to gene symbols'
          f' (remaining kept as WBID; wormbase_id stored in .var)')

    if args.drop_cuticle:
        exclude = cuticle_exclude_symbols(gid2sym)
        mask = ~adata.var_names.isin(exclude)
        n_drop = int((~mask).sum())
        adata = adata[:, mask].copy()
        print(f'Dropped {n_drop} cuticle genes  → {adata.n_vars} genes remain')

    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Always save the loaded integer counts before any transformation.
    adata.layers['counts'] = adata.X.copy()
    if adata.layers['counts'].dtype.kind not in 'iu':
        print(f"Warning: input .X dtype is {adata.layers['counts'].dtype}, not integer; "
              "SoupX / scVI assume raw integer counts in .layers['counts']", file=sys.stderr)

    # --- Populate diagnostics into .var (always when CSVs provided) -----------
    import pandas as pd
    wbids = adata.var['wormbase_id']
    have_junc = args.junction_csv is not None
    have_len = args.gene_lengths_csv is not None
    if have_junc:
        if not args.junction_csv.is_file():
            print(f'Error: junction CSV not found: {args.junction_csv}', file=sys.stderr)
            return 1
        junc = pd.read_csv(args.junction_csv).set_index('gene_id')
        adata.var['n_junctions'] = junc['n_junctions'].reindex(wbids).fillna(0).to_numpy()
    if have_len:
        if not args.gene_lengths_csv.is_file():
            print(f'Error: gene lengths CSV not found: {args.gene_lengths_csv}', file=sys.stderr)
            return 1
        lens = pd.read_csv(args.gene_lengths_csv).set_index('gene_id')
        adata.var['exonic_length'] = lens['exonic_length'].reindex(wbids).fillna(0).to_numpy()
    if have_junc and have_len:
        adata.var['detectability_k'] = (
            adata.var['n_junctions'].to_numpy() * adata.var['exonic_length'].to_numpy() / 1000.0
        )

    # --- Apply count-scaling arm to .X ----------------------------------------
    print(f'Count-scaling arm: {arm}')
    if arm != 'none':
        if 'detectability_k' not in adata.var.columns:
            print("Error: --count-scaling != 'none' requires both --junction-csv and --gene-lengths-csv",
                  file=sys.stderr)
            return 1
        rng = np.random.default_rng(args.seed)
        adata.X = _apply_scaling(adata.X, adata.var['detectability_k'].to_numpy(), arm, rng)
        n_has = int((adata.var['detectability_k'] > 0).sum())
        print(f'  scaled with k = n_junctions · exonic_length_kb  '
              f'({n_has}/{adata.n_vars} genes with k > 0)')

    # Defensive: raw ints preserved regardless of arm.
    counts = adata.layers['counts']
    counts_dtype = counts.dtype
    if counts_dtype.kind not in 'iu':
        print(f"Warning: .layers['counts'] dtype is {counts_dtype} after scaling; "
              "downstream scVI expects integer counts", file=sys.stderr)

    # --- Assay label ---------------------------------------------------------
    adata.obs[args.assay_key] = args.assay_value
    adata.obs[args.assay_key] = adata.obs[args.assay_key].astype('category')
    print(f'Set .obs[{args.assay_key!r}] = {args.assay_value!r}')

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    batch_key = args.batch_key if args.batch_key in adata.obs else None
    if batch_key is None and args.batch_key:
        print(f"Warning: batch_key '{args.batch_key}' not in .obs — HVGs computed without batching")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, batch_key=batch_key)

    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, min_dist=args.min_dist, spread=args.spread)
    sc.tl.leiden(adata, flavor='igraph', n_iterations=2, resolution=args.resolution)

    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output_h5ad)
    print(f'Wrote {args.output_h5ad}')
    print(f'  {adata.n_obs} cells x {adata.n_vars} genes  |  {adata.obs["leiden"].nunique()} leiden clusters')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
