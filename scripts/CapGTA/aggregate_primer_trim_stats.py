#!/usr/bin/env python3
"""Collapse per-cell <barcode>_stats.json into a single TSV across all cells.

Usage:
    python aggregate_primer_trim_stats.py <output_root> [--out all_cells_stats.tsv]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


SCALARS = [
    'barcode', 'udi', 'total_pairs', 'pairs_no_polyT', 'pairs_bad_umi',
    'pairs_dropped_short_after_trim', 'pairs_kept', 'frac_kept',
    'run_min', 'max_mm', 'qual_threshold', 'min_cdna_len',
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('output_root', help='Root that holds <UDI>/<barcode>_stats.json files')
    ap.add_argument('--out', default=None,
                    help='Output TSV path (default: <output_root>/all_cells_stats.tsv)')
    args = ap.parse_args()

    root = Path(args.output_root)
    out_tsv = Path(args.out) if args.out else root / 'all_cells_stats.tsv'

    rows = []
    for stats_file in sorted(root.rglob('*_stats.json')):
        barcode = stats_file.name[: -len('_stats.json')]
        udi = stats_file.parent.name
        try:
            with open(stats_file) as fh:
                d = json.load(fh)
        except Exception as e:
            print(f'warn: skipping {stats_file}: {e}', file=sys.stderr)
            continue
        d['barcode'] = barcode
        d['udi'] = udi
        rows.append({k: d.get(k, '') for k in SCALARS})

    with open(out_tsv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=SCALARS, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'wrote {len(rows)} cells → {out_tsv}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
