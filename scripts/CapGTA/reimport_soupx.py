#!/usr/bin/env python3
"""Fold SoupX decontaminated counts back into an AnnData.

The R soupx_decontaminate.R script writes counts_decontam.mtx + genes.tsv +
barcodes.tsv + soupx_summary.csv into <soupx_dir>. This script builds a fresh
h5ad from those outputs plus a preprocessed cells h5ad (source of obs/var
metadata and pre-SoupX raw counts).

The output h5ad has:
    .X                      SoupX decontam integer counts (int32 sparse)
    .layers['counts_raw']   pre-SoupX raw integer counts (from cells.layers['counts'])
    .obs                    cells.obs minus stale clustering (leiden, soupx_cluster
                            preserved) + SoupX per-cell metadata
    .var                    cells.var (schema only — QC/HVG columns dropped)
    .uns['soupx']           run metadata

Clustering (obsm['X_pca'], ['X_umap'], obsp neighbours, uns keys) is NOT
carried over; recompute downstream with preprocess_scanpy.py.

Usage:
    reimport_soupx.py <cells_h5ad> <soupx_dir> <output_h5ad>
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cells_h5ad', type=Path)
    p.add_argument('soupx_dir', type=Path)
    p.add_argument('output_h5ad', type=Path)
    args = p.parse_args()

    if not args.cells_h5ad.is_file():
        print(f'Error: cells h5ad not found: {args.cells_h5ad}', file=sys.stderr)
        return 1

    cells = sc.read_h5ad(args.cells_h5ad)
    if 'counts' not in cells.layers:
        print("Error: cells h5ad must have .layers['counts'] with raw integer counts", file=sys.stderr)
        return 1

    mtx_path = args.soupx_dir / 'counts_decontam.mtx'
    for f in (mtx_path, args.soupx_dir / 'genes.tsv',
              args.soupx_dir / 'barcodes.tsv',
              args.soupx_dir / 'soupx_summary.csv'):
        if not f.is_file():
            print(f'Error: SoupX output not found: {f}', file=sys.stderr)
            return 1

    X_dec = sio.mmread(str(mtx_path)).T.tocsr().astype(np.int32)   # cells x genes
    genes = (args.soupx_dir / 'genes.tsv').read_text().splitlines()
    barcodes = (args.soupx_dir / 'barcodes.tsv').read_text().splitlines()
    if genes != list(cells.var_names):
        print('Error: gene order in SoupX output does not match cells h5ad', file=sys.stderr)
        return 1
    if barcodes != list(cells.obs_names):
        print('Error: barcode order in SoupX output does not match cells h5ad', file=sys.stderr)
        return 1

    summary = pd.read_csv(args.soupx_dir / 'soupx_summary.csv') \
                .set_index('barcode').loc[cells.obs_names]

    # obs: keep upstream metadata (QC, sample, etc.), drop stale clustering.
    drop_cols = [c for c in ('leiden',) if c in cells.obs]
    obs_out = cells.obs.drop(columns=drop_cols).copy()
    obs_out['soupx_rho']        = summary['rho'].astype(np.float32).values
    obs_out['soupx_nUMIs_pre']  = summary['nUMIs'].astype(np.int32).values
    obs_out['soupx_nUMIs_post'] = summary['nUMIs_after'].astype(np.int32).values
    obs_out['soupx_cluster']    = summary['cluster'].astype(str).values

    keep_var_cols = [c for c in ('wormbase_id',) if c in cells.var.columns]
    new = sc.AnnData(
        X=X_dec.copy(),
        obs=obs_out,
        var=cells.var[keep_var_cols].copy(),
        layers={'counts_raw': cells.layers['counts'].tocsr().astype(np.int32)},
    )
    new.uns['soupx'] = {
        'method':          'SoupX autoEstCont',
        'global_rho_mean': float(summary['rho'].mean()),
        'global_rho_med':  float(summary['rho'].median()),
        'n_cells':         int(new.n_obs),
    }

    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    new.write_h5ad(args.output_h5ad)

    n_pre  = np.asarray(new.layers['counts_raw'].sum(axis=1)).ravel()
    n_post = np.asarray(new.X.sum(axis=1)).ravel()
    print(f'Wrote {args.output_h5ad}')
    print(f'  {new.n_obs} cells x {new.n_vars} genes')
    print(f'  Median counts/cell: {int(np.median(n_pre)):,} → {int(np.median(n_post)):,}'
          f'  ({100*(1 - np.median(n_post)/max(np.median(n_pre),1)):.1f}% reduction)')
    print(f'  rho: median={new.obs["soupx_rho"].median():.3f}  mean={new.obs["soupx_rho"].mean():.3f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
