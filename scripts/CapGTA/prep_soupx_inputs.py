#!/usr/bin/env python3
"""Prep SoupX inputs from a cells h5ad and an empties h5ad.

Writes into <output_dir>:
    counts.mtx / genes.tsv / barcodes.tsv    — real cells (from adata.layers['counts'].T)
    clusters.csv                              — barcode,cluster (from --cluster-key)
    empties/counts.mtx / genes.tsv / barcodes.tsv  — filtered empties, aligned to cell genes

Empty droplet QC filter (WGA-aware):
    - drop empties with total counts < --min-counts
    - drop empties with total counts > <percentile>th percentile (--max-percentile)
    - drop empties where a single gene contributes > --max-gene-frac of that empty's counts
      (guards against WGA-amplified single transcripts poisoning the soup profile)

Usage:
    prep_soupx_inputs.py <cells_h5ad> <empties_h5ad> <output_dir> \
        [--cluster-key leiden] [--min-counts 5] [--max-percentile 95] [--max-gene-frac 0.30]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp


def _get_counts(adata):
    """Return raw integer counts as CSR (cells x genes)."""
    if 'counts' in adata.layers:
        X = adata.layers['counts']
    else:
        X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    return X.tocsr().astype(np.int32)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cells_h5ad', type=Path)
    p.add_argument('empties_h5ad', type=Path)
    p.add_argument('output_dir', type=Path)
    p.add_argument('--cluster-key', default='leiden')
    p.add_argument('--min-counts', type=int, default=5)
    p.add_argument('--max-percentile', type=float, default=95.0)
    p.add_argument('--max-gene-frac', type=float, default=0.30)
    args = p.parse_args()

    for f in (args.cells_h5ad, args.empties_h5ad):
        if not f.is_file():
            print(f'Error: not found: {f}', file=sys.stderr)
            return 1

    print(f'Loading cells: {args.cells_h5ad}')
    cells = sc.read_h5ad(args.cells_h5ad)
    print(f'  {cells.n_obs} cells x {cells.n_vars} genes')

    if args.cluster_key not in cells.obs:
        print(f'Error: --cluster-key "{args.cluster_key}" not in cells.obs', file=sys.stderr)
        return 1

    print(f'Loading empties: {args.empties_h5ad}')
    empties = sc.read_h5ad(args.empties_h5ad)
    print(f'  {empties.n_obs} droplets x {empties.n_vars} genes')

    # If cells were renamed WBID → symbol upstream, remap empties the same way
    # so downstream alignment on var_names works.
    if 'wormbase_id' in cells.var.columns:
        wbid_to_name = dict(zip(cells.var['wormbase_id'].astype(str), cells.var_names))
        new_names = np.array([wbid_to_name.get(g, g) for g in empties.var_names])
        empties.var['wormbase_id'] = empties.var_names.to_numpy()
        empties.var_names = new_names
        empties.var_names_make_unique()
        n_mapped = int(sum(g in wbid_to_name for g in empties.var['wormbase_id']))
        print(f'  Remapped {n_mapped}/{empties.n_vars} empties var_names via cells wormbase_id')

    # Align empties to cells gene order
    common = cells.var_names
    missing = set(common) - set(empties.var_names)
    if missing:
        print(f'  Note: {len(missing)} cell genes absent in empties — will be zero-padded')
    empties_aligned = empties[:, [g for g in common if g in empties.var_names]].copy()
    # Reindex to full cell gene set with zeros
    empties_full = sc.AnnData(
        X=sp.lil_matrix((empties.n_obs, len(common)), dtype=np.int32),
        obs=empties.obs.copy(),
        var=cells.var[[]].copy(),
    )
    idx = {g: i for i, g in enumerate(common)}
    cols = [idx[g] for g in empties_aligned.var_names]
    Xa = _get_counts(empties_aligned).tolil()
    empties_full.X = sp.lil_matrix(empties_full.shape, dtype=np.int32)
    empties_full.X[:, cols] = Xa
    empties_full.X = empties_full.X.tocsr()

    # QC
    E = empties_full.X
    total = np.asarray(E.sum(axis=1)).ravel()
    max_gene = np.asarray(E.max(axis=1).todense()).ravel()
    max_frac = np.where(total > 0, max_gene / np.maximum(total, 1), 0)

    upper = np.percentile(total, args.max_percentile)
    keep = (total >= args.min_counts) & (total <= upper) & (max_frac <= args.max_gene_frac)

    print()
    print('Empty droplet QC:')
    print(f'  starting empties:               {len(total)}')
    print(f'  total counts: min={total.min():.0f}  median={int(np.median(total))}  max={int(total.max())}  {args.max_percentile:.0f}th pct={int(upper)}')
    print(f'  max-single-gene fraction: median={np.median(max_frac):.3f}  90th={np.percentile(max_frac, 90):.3f}  max={max_frac.max():.3f}')
    print(f'  drop (total < {args.min_counts}):         {int((total < args.min_counts).sum())}')
    print(f'  drop (total > {int(upper)}):           {int((total > upper).sum())}')
    print(f'  drop (max gene frac > {args.max_gene_frac}):  {int((max_frac > args.max_gene_frac).sum())}')
    print(f'  KEPT:                           {int(keep.sum())}')

    # Write cells
    args.output_dir.mkdir(parents=True, exist_ok=True)
    C = _get_counts(cells)                       # cells x genes
    Ct = C.T.tocoo()                             # genes x cells
    sio.mmwrite(str(args.output_dir / 'counts.mtx'), Ct, field='integer')
    (args.output_dir / 'genes.tsv').write_text('\n'.join(cells.var_names) + '\n')
    (args.output_dir / 'barcodes.tsv').write_text('\n'.join(cells.obs_names) + '\n')
    cells.obs[[args.cluster_key]].rename(columns={args.cluster_key: 'cluster'}) \
         .rename_axis('barcode').to_csv(args.output_dir / 'clusters.csv')
    print(f'\nWrote cells: {args.output_dir}/counts.mtx  ({Ct.shape[0]} x {Ct.shape[1]})')

    # Write filtered empties
    emp_dir = args.output_dir / 'empties'
    emp_dir.mkdir(exist_ok=True)
    Ef = E[keep]                                  # empties x genes
    Et = Ef.T.tocoo()                             # genes x empties
    sio.mmwrite(str(emp_dir / 'counts.mtx'), Et, field='integer')
    (emp_dir / 'genes.tsv').write_text('\n'.join(cells.var_names) + '\n')
    kept_barcodes = [b for b, k in zip(empties_full.obs_names, keep) if k]
    (emp_dir / 'barcodes.tsv').write_text('\n'.join(kept_barcodes) + '\n')
    print(f'Wrote empties: {emp_dir}/counts.mtx  ({Et.shape[0]} x {Et.shape[1]})')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
