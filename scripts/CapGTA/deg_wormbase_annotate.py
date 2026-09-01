#!/usr/bin/env python3
"""Rank genes per cluster and annotate top DEGs from WormBase.

Loads a CapGTA h5ad, normalizes+log if .X looks like counts, runs
sc.tl.rank_genes_groups on the given cluster key, pulls the top N genes per
cluster, and queries the WormBase REST API for each unique WBGene ID to fetch
gene name, gene class, and concise description.

Usage:
    deg_wormbase_annotate.py <input.h5ad> <output.csv> \
        [--cluster-key leiden] [--resolution 0.6] [--top-n 10] [--method wilcoxon]

If --resolution is given, computes leiden at that resolution (into
.obs['leiden_<res>']) and uses that as the grouping. When var_names are
gene symbols, the WBID for each is pulled from .var['wormbase_id'].
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import requests
import scanpy as sc


WB_BASE = 'https://rest.wormbase.org/rest/field/gene'


def wb_get(wbid: str, field: str, session: requests.Session, retries: int = 3) -> dict:
    url = f'{WB_BASE}/{wbid}/{field}'
    for attempt in range(retries):
        try:
            r = session.get(url, headers={'Accept': 'application/json'}, timeout=20)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if attempt == retries - 1:
                print(f'  WormBase fetch failed for {wbid}/{field}: {e}', file=sys.stderr)
                return {}
            time.sleep(1.0 * (attempt + 1))
    return {}


def annotate_gene(wbid: str, session: requests.Session) -> dict:
    name_js = wb_get(wbid, 'name', session)
    desc_js = wb_get(wbid, 'concise_description', session)
    class_js = wb_get(wbid, 'gene_class', session)

    name = ''
    try:
        name = name_js['name']['data']['label'] or ''
    except (KeyError, TypeError):
        pass

    description = ''
    try:
        description = desc_js['concise_description']['data']['text'] or ''
    except (KeyError, TypeError):
        pass

    gene_class = ''
    try:
        d = class_js['gene_class']['data']
        if d:
            tag = d.get('tag') or {}
            label = tag.get('label') or ''
            full = d.get('description') or ''
            gene_class = f'{label} ({full})' if label and full else (label or full)
    except (KeyError, TypeError):
        pass

    return {'gene_name': name, 'gene_class': gene_class, 'description': description}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input_h5ad', type=Path)
    p.add_argument('output_csv', type=Path)
    p.add_argument('--cluster-key', default='leiden')
    p.add_argument('--resolution', type=float, default=None,
                   help='if set, compute leiden at this resolution into .obs and use it')
    p.add_argument('--top-n', type=int, default=10)
    p.add_argument('--method', default='wilcoxon')
    p.add_argument('--sleep', type=float, default=0.1, help='seconds between WormBase requests')
    args = p.parse_args()

    if not args.input_h5ad.is_file():
        print(f'Error: input not found: {args.input_h5ad}', file=sys.stderr)
        return 1

    print(f'Loading {args.input_h5ad}')
    adata = sc.read_h5ad(args.input_h5ad)
    print(f'  {adata.n_obs} cells x {adata.n_vars} genes')

    if args.resolution is not None:
        cluster_key = f'leiden_{args.resolution}'
        print(f'Computing leiden at resolution {args.resolution} → .obs[{cluster_key!r}]')
        sc.tl.leiden(adata, flavor='igraph', n_iterations=2,
                     resolution=args.resolution, key_added=cluster_key)
    else:
        cluster_key = args.cluster_key

    if cluster_key not in adata.obs:
        print(f'Error: cluster key {cluster_key!r} not in .obs. Available: {list(adata.obs.columns)}', file=sys.stderr)
        return 1

    xmax = float(adata.X.max())
    if xmax > 50:
        print(f'.X looks like counts (max={xmax:.0f}); normalizing + log1p in place')
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    else:
        print(f'.X looks already log-transformed (max={xmax:.2f}); skipping normalization')

    # Map var_names → WBID. If var_names look like symbols and .var['wormbase_id']
    # exists, use it; otherwise assume var_names ARE WBIDs.
    if 'wormbase_id' in adata.var.columns and not str(adata.var_names[0]).startswith('WBGene'):
        name_to_wbid = dict(zip(adata.var_names, adata.var['wormbase_id'].astype(str)))
        print('Var_names are gene symbols; mapping to WBIDs via .var[wormbase_id]')
    else:
        name_to_wbid = {n: n for n in adata.var_names}

    print(f'Running rank_genes_groups on {cluster_key!r} ({args.method})')
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=args.method)

    names = adata.uns['rank_genes_groups']['names']
    clusters = list(names.dtype.names)
    print(f'  {len(clusters)} clusters')

    rows = []
    for c in clusters:
        for rank in range(args.top_n):
            symbol = names[c][rank]
            wbid = name_to_wbid.get(symbol, symbol)
            rows.append({'cluster': c, 'rank': rank + 1, 'wormbase_id': wbid, 'var_name': symbol})

    unique_ids = sorted({r['wormbase_id'] for r in rows})
    print(f'Querying WormBase for {len(unique_ids)} unique genes')

    session = requests.Session()
    cache = {}
    for i, wbid in enumerate(unique_ids, 1):
        cache[wbid] = annotate_gene(wbid, session)
        if i % 25 == 0 or i == len(unique_ids):
            print(f'  {i}/{len(unique_ids)}')
        time.sleep(args.sleep)

    for r in rows:
        r.update(cache[r['wormbase_id']])

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['cluster', 'rank', 'var_name', 'wormbase_id', 'gene_name', 'gene_class', 'description'])
        w.writeheader()
        w.writerows(rows)
    print(f'Wrote {args.output_csv}  ({len(rows)} rows)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
