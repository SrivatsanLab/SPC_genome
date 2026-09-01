#!/usr/bin/env python3
"""Scanpy preprocessing for CapGTA spliced count h5ads.

Pipeline:
    load h5ad
    rename var_names WBID → gene symbol (WBID kept in .var['wormbase_id'])
    [if --drop-cuticle: drop col-* + dpy-7/sqt-1/sqt-3/rol-6/bli-1/bli-2]
    calculate_qc_metrics
    save raw integer counts to .layers['counts']   (unconditional; overwrites)
    [scale .X columns by 1 / (n_junctions * exonic_length_kb) when the matching
        --junction-csv / --gene-lengths-csv are provided; either factor omitted
        when its CSV is not passed; single-exon (n_junctions = 0) genes → 0]
    normalize_total → log1p
    highly_variable_genes (n_top_genes, batch_key)
    pca
    neighbors
    umap (min_dist, spread)
    leiden (flavor='igraph', n_iterations=2, resolution)
    write h5ad

Junction + length correction is applied to .X only; .layers['counts'] always
holds the raw integer counts as-loaded, which is what SoupX and other
count-model tools expect. Both CSVs are keyed on gene_id (WBID), so reindexing
uses .var['wormbase_id'] not var_names.

Usage:
    preprocess_scanpy.py <input.h5ad> <output.h5ad> --gtf <annotation.gtf> \
        [--drop-cuticle] \
        [--junction-csv <gene_junctions.csv>] \
        [--gene-lengths-csv <gene_lengths.csv>] \
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


CUTICLE_PREFIXES = ('col-', 'cut-', 'idpa-')
CUTICLE_SYMBOLS = {
    # cuticle-collagen dpy genes (dpy-1, dpy-11, dpy-18+ have non-cuticle roles → excluded)
    'dpy-2', 'dpy-3', 'dpy-4', 'dpy-5', 'dpy-6', 'dpy-7',
    'dpy-8', 'dpy-9', 'dpy-10', 'dpy-13', 'dpy-14', 'dpy-17',
    # squat / roller / blister cuticle families
    'sqt-1', 'sqt-2', 'sqt-3',
    'rol-1', 'rol-3', 'rol-4', 'rol-5', 'rol-6', 'rol-8',
    'bli-1', 'bli-2', 'bli-3', 'bli-4', 'bli-5', 'bli-6',
}


def cuticle_exclude_symbols(gid2sym: dict[str, str]) -> set[str]:
    """Cuticle gene symbols: all col-/cut-/idpa- family members + curated cuticle set."""
    syms = {s for s in gid2sym.values() if s.lower().startswith(CUTICLE_PREFIXES)}
    syms.update(CUTICLE_SYMBOLS)
    return syms


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input_h5ad', type=Path)
    p.add_argument('output_h5ad', type=Path)
    p.add_argument('--gtf', type=Path, required=True,
                   help='GTF used to map gene_id → gene_name for var rename')
    p.add_argument('--drop-cuticle', action='store_true',
                   help='drop col-* collagens plus dpy-7, sqt-1, sqt-3, rol-6, bli-1, bli-2')
    p.add_argument('--junction-csv', type=Path, default=None,
                   help='gene_junctions.csv (columns: gene_id,n_junctions) for splice-count correction')
    p.add_argument('--gene-lengths-csv', type=Path, default=None,
                   help='gene_lengths.csv (columns: gene_id,exonic_length) for cDNA-length correction')
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

    if args.junction_csv is not None or args.gene_lengths_csv is not None:
        import pandas as pd
        wbids = adata.var['wormbase_id']
        n_genes = adata.n_vars
        scale = np.ones(n_genes, dtype=float)
        keep = np.ones(n_genes, dtype=bool)
        parts = []

        if args.junction_csv is not None:
            if not args.junction_csv.is_file():
                print(f'Error: junction CSV not found: {args.junction_csv}', file=sys.stderr)
                return 1
            junc = pd.read_csv(args.junction_csv).set_index('gene_id')
            n_junc = junc['n_junctions'].reindex(wbids).fillna(0).to_numpy()
            adata.var['n_junctions'] = n_junc
            has_junc = n_junc > 0
            scale = np.where(has_junc, scale / np.where(has_junc, n_junc, 1), 0.0)
            keep &= has_junc
            parts.append(f'{int(has_junc.sum())}/{n_genes} with >=1 junction')

        if args.gene_lengths_csv is not None:
            if not args.gene_lengths_csv.is_file():
                print(f'Error: gene lengths CSV not found: {args.gene_lengths_csv}', file=sys.stderr)
                return 1
            lens = pd.read_csv(args.gene_lengths_csv).set_index('gene_id')
            exonic_bp = lens['exonic_length'].reindex(wbids).fillna(0).to_numpy()
            adata.var['exonic_length'] = exonic_bp
            has_len = exonic_bp > 0
            length_kb = exonic_bp / 1000.0
            scale = np.where(has_len, scale / np.where(has_len, length_kb, 1), 0.0)
            keep &= has_len
            parts.append(f'{int(has_len.sum())}/{n_genes} with exonic length')

        scale[~keep] = 0.0
        X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
        adata.X = (X @ sp.diags(scale)).tocsr()
        factor_names = []
        if args.junction_csv is not None:
            factor_names.append('n_junctions')
        if args.gene_lengths_csv is not None:
            factor_names.append('exonic_length_kb')
        print(f'Scaled .X by 1 / ({" * ".join(factor_names)})  ({"; ".join(parts)})')

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
