#!/usr/bin/env python3
"""Enumerate empty droplet BAMs for SoupX ambient RNA estimation.

Reads a real_cells.csv (columns: barcode,bam_path), discovers the parent
directories of the real-cell BAMs, lists every *.bam in those directories,
subtracts real-cell barcodes, and random-samples N empties.

Usage:
    enumerate_empties.py <real_cells.csv> <output.csv> [--n-empties 5000] [--seed 0]
"""

import argparse
import csv
import random
import sys
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('real_cells_csv', type=Path)
    p.add_argument('output_csv', type=Path)
    p.add_argument('--n-empties', type=int, default=5000)
    p.add_argument('--seed', type=int, default=0)
    args = p.parse_args()

    if not args.real_cells_csv.is_file():
        print(f'Error: not found: {args.real_cells_csv}', file=sys.stderr)
        return 1

    real_barcodes: set[str] = set()
    parent_dirs: set[str] = set()
    with args.real_cells_csv.open() as fh:
        reader = csv.reader(fh)
        next(reader)  # header
        for row in reader:
            real_barcodes.add(row[0])
            parent_dirs.add(str(Path(row[1]).parent))

    print(f'Real cells: {len(real_barcodes)}   Parent BAM dirs: {len(parent_dirs)}')

    all_empties: list[tuple[str, str]] = []
    for d in sorted(parent_dirs):
        for p_bam in Path(d).glob('*.bam'):
            bc = p_bam.stem
            if bc not in real_barcodes:
                all_empties.append((bc, str(p_bam)))
    print(f'Candidate empty droplets: {len(all_empties)}')

    if not all_empties:
        print('Error: no empty droplet BAMs found', file=sys.stderr)
        return 1

    random.seed(args.seed)
    sample = random.sample(all_empties, min(args.n_empties, len(all_empties)))
    print(f'Sampled: {len(sample)}')

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['barcode', 'bam_path'])
        for bc, p in sample:
            w.writerow([bc, p])
    print(f'Wrote {args.output_csv}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
