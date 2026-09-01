#!/usr/bin/env python3
"""Build the wide gene × anatomy-term dictionary CSV that the
`tissue_enrichment_analysis` package expects.

Downloads the current WormBase anatomy association (GAF) and anatomy ontology
(OBO), pivots to a WBID × WBbt-name binary matrix, and writes anatomy_dict.csv.
This is a drop-in replacement for the package's (now-404) fetch_dictionary().

Usage:
    build_wormbase_anatomy_dict.py <output.csv> [--release WS298]
"""

import argparse
import gzip
import io
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd


BASE = 'https://downloads.wormbase.org/releases/current-production-release/ONTOLOGY'


def download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def parse_obo_names(obo_bytes: bytes) -> dict[str, str]:
    """Parse OBO for id → name mappings (Term stanzas only)."""
    text = gzip.decompress(obo_bytes).decode('utf-8', errors='replace')
    names: dict[str, str] = {}
    cur_id = None
    in_term = False
    for line in text.splitlines():
        line = line.strip()
        if line == '[Term]':
            in_term = True
            cur_id = None
            continue
        if line.startswith('['):
            in_term = False
            continue
        if not in_term:
            continue
        if line.startswith('id: '):
            cur_id = line[4:].strip()
        elif line.startswith('name: ') and cur_id is not None:
            names[cur_id] = line[6:].strip()
    return names


def parse_gaf(gaf_bytes: bytes) -> dict[str, set[str]]:
    """WBID → set(WBbt term IDs). Excludes qualifier 'Not' negative annotations."""
    text = gzip.decompress(gaf_bytes).decode('utf-8', errors='replace')
    m: dict[str, set[str]] = defaultdict(set)
    for line in text.splitlines():
        if not line or line.startswith('!'):
            continue
        fields = line.split('\t')
        if len(fields) < 5:
            continue
        wbid = fields[1]
        qualifier = fields[3]
        term = fields[4]
        if qualifier and qualifier.lower().startswith('not'):
            continue
        if not wbid.startswith('WBGene') or not term.startswith('WBbt:'):
            continue
        m[wbid].add(term)
    return m


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('output_csv', type=Path)
    p.add_argument('--release', default='WS298', help='WormBase release (default: WS298)')
    args = p.parse_args()

    gaf_url = f'{BASE}/anatomy_association.{args.release}.wb.gz'
    obo_url = f'{BASE}/anatomy_ontology.{args.release}.obo.gz'

    print(f'Downloading {gaf_url}')
    gaf_bytes = download(gaf_url)
    print(f'Downloading {obo_url}')
    obo_bytes = download(obo_url)

    print('Parsing OBO for term names…')
    id2name = parse_obo_names(obo_bytes)
    print(f'  {len(id2name)} terms')

    print('Parsing GAF associations…')
    wbid2terms = parse_gaf(gaf_bytes)
    print(f'  {len(wbid2terms)} genes with annotations')

    # Only keep terms that appear in the associations
    used_term_ids = sorted({t for ts in wbid2terms.values() for t in ts})
    used_term_names = [id2name.get(t, t) for t in used_term_ids]
    term_to_col = {t: i for i, t in enumerate(used_term_ids)}
    print(f'  {len(used_term_ids)} unique terms after filtering')

    print('Building wide binary matrix…')
    genes = sorted(wbid2terms.keys())
    rows = []
    for wbid in genes:
        row = [0] * len(used_term_ids)
        for t in wbid2terms[wbid]:
            row[term_to_col[t]] = 1
        rows.append(row)

    df = pd.DataFrame(rows, columns=used_term_names)
    df.insert(0, 'wbid', genes)

    # Collapse duplicate column names (multiple WBbt IDs sharing the same name)
    df = df.groupby(level=0, axis=1).max() if df.columns.duplicated().any() else df
    if 'wbid' not in df.columns[:1]:
        # groupby may have reordered; re-insert wbid as first column
        wbids = df['wbid']
        df = df.drop(columns=['wbid'])
        df.insert(0, 'wbid', wbids)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f'Wrote {args.output_csv}  ({df.shape[0]} genes x {df.shape[1] - 1} terms)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
