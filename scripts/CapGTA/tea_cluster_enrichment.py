#!/usr/bin/env python3
"""Run WormBase tissue / phenotype / GO enrichment (TEA) on top-N DEGs per cluster.

Loads a CapGTA h5ad, runs sc.tl.rank_genes_groups on the given cluster key,
takes the top-N genes per cluster, maps them to WBIDs (via .var['wormbase_id']
if present), and submits each cluster's list to WormBase's enrichment tool via
the `tissue_enrichment_analysis` Python package (same backend as
https://www.wormbase.org/tools/enrichment/tea/tea.cgi).

Requires:  pip install tissue_enrichment_analysis

The package's built-in fetch_dictionary() is broken (caltech.wormbase.org URL
404s). Build a local dict CSV with build_wormbase_anatomy_dict.py and pass
--dict-csv, OR let it try the package's fetch first.

Usage:
    tea_cluster_enrichment.py <input.h5ad> <output.csv> \
        [--cluster-key leiden] [--top-n 50] \
        [--analysis tissue|phenotype|go] [--alpha 0.05] \
        [--dict-csv <local anatomy_dict.csv>]

Output CSV columns:
    cluster, Term, Expected, Observed, Enrichment Fold Change, P value, Q value
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input_h5ad', type=Path)
    p.add_argument('output_csv', type=Path)
    p.add_argument('--cluster-key', default='leiden')
    p.add_argument('--top-n', type=int, default=50)
    p.add_argument('--analysis', choices=['tissue', 'phenotype', 'go'], default='tissue')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--dict-csv', type=Path, default=None,
                   help='local anatomy_dict.csv path (skips package fetch_dictionary)')
    args = p.parse_args()

    try:
        import tissue_enrichment_analysis as tea
    except ImportError:
        print('Error: tissue_enrichment_analysis not installed.', file=sys.stderr)
        print('Install with:  pip install tissue_enrichment_analysis', file=sys.stderr)
        return 1

    if not args.input_h5ad.is_file():
        print(f'Error: input not found: {args.input_h5ad}', file=sys.stderr)
        return 1

    print(f'Loading {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)
    print(f'  {adata.n_obs} cells x {adata.n_vars} genes')

    if args.cluster_key not in adata.obs:
        print(f'Error: cluster key {args.cluster_key!r} not in .obs. '
              f'Available: {list(adata.obs.columns)}', file=sys.stderr)
        return 1

    if 'wormbase_id' in adata.var.columns:
        sym_to_wbid = dict(zip(adata.var_names, adata.var['wormbase_id'].astype(str)))
        print('Mapping var_names → WBIDs via .var[wormbase_id]')
    else:
        sym_to_wbid = {n: n for n in adata.var_names}

    xmax = float(adata.X.max())
    if xmax > 50:
        print(f'.X looks like counts (max={xmax:.0f}); normalizing + log1p')
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    print(f'Running rank_genes_groups on {args.cluster_key!r} (wilcoxon)')
    sc.tl.rank_genes_groups(adata, groupby=args.cluster_key, method='wilcoxon')

    names = adata.uns['rank_genes_groups']['names']
    clusters = list(names.dtype.names)
    print(f'  {len(clusters)} clusters, top-{args.top_n} genes each')

    if args.dict_csv is not None:
        if not args.dict_csv.is_file():
            print(f'Error: dict CSV not found: {args.dict_csv}', file=sys.stderr)
            return 1
        print(f'Loading local dictionary from {args.dict_csv}')
        ann_df = pd.read_csv(args.dict_csv)
    else:
        print(f'Fetching WormBase {args.analysis!r} dictionary…')
        ann_df = tea.fetch_dictionary(analysis=args.analysis)
    if ann_df is None or len(ann_df) == 0:
        print('Error: failed to load annotation dictionary (fetch returned None). '
              'Build a local one with build_wormbase_anatomy_dict.py and pass --dict-csv.',
              file=sys.stderr)
        return 1
    print(f'  dictionary shape: {ann_df.shape[0]} genes x {ann_df.shape[1] - 1} terms')

    all_rows = []
    for c in clusters:
        n_available = len(names[c])
        symbols = [names[c][r] for r in range(min(args.top_n, n_available))]
        wbids = [sym_to_wbid.get(s, s) for s in symbols]
        try:
            enrich_df = tea.enrichment_analysis(wbids, ann_df, alpha=args.alpha, show=False)
        except Exception as e:
            print(f'  cluster {c}: enrichment failed: {e}', file=sys.stderr)
            continue
        if enrich_df is None or len(enrich_df) == 0:
            print(f'  cluster {c}: no enriched terms at alpha={args.alpha}')
            continue
        enrich_df = enrich_df.copy()
        enrich_df.insert(0, 'cluster', c)
        all_rows.append(enrich_df)
        print(f'  cluster {c}: {len(enrich_df)} enriched terms')

    if not all_rows:
        print('No enriched terms found for any cluster', file=sys.stderr)
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        args.output_csv.write_text('')
        return 1

    combined = pd.concat(all_rows, ignore_index=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output_csv, index=False)
    print(f'Wrote {args.output_csv}  ({len(combined)} rows)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
