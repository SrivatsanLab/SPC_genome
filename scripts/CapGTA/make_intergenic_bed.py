#!/usr/bin/env python3
"""Build an intergenic BED from a GTF + FASTA index.

Intergenic = regions of the genome that are >= flank bp away from any gene
body (min exon start to max exon end per gene, unioned across isoforms) AND
at least min_size bp wide. Used as the DNA-only baseline for per-cell
background rate estimation in the RNA expression pipeline.

Usage:
    make_intergenic_bed.py <annotation.gtf> <reference.fa.fai> <output.bed>
        [--flank 5000] [--min-size 10000]
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_attributes(attr_str: str) -> dict:
    out = {}
    for kv in attr_str.strip().rstrip(';').split(';'):
        kv = kv.strip()
        if not kv:
            continue
        parts = kv.split(' ', 1)
        if len(parts) == 2:
            key, val = parts
            out[key] = val.strip().strip('"')
    return out


def load_gene_bodies(gtf_path: Path) -> dict:
    """Return {chrom: [(start, end), ...]} of gene bodies (1-based inclusive)."""
    gene_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    gene_chrom: dict[str, str] = {}
    with gtf_path.open() as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            chrom, _, _, start, end, _, _, _, attrs = fields
            gid = parse_attributes(attrs).get('gene_id')
            if gid is None:
                continue
            gene_exons[gid].append((int(start), int(end)))
            gene_chrom.setdefault(gid, chrom)

    bodies_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for gid, exons in gene_exons.items():
        starts = [s for s, _ in exons]
        ends = [e for _, e in exons]
        bodies_by_chrom[gene_chrom[gid]].append((min(starts), max(ends)))
    return bodies_by_chrom


def load_chrom_sizes(fai_path: Path) -> dict:
    sizes = {}
    with fai_path.open() as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            sizes[fields[0]] = int(fields[1])
    return sizes


def merge_intervals(intervals: list) -> list:
    """Merge overlapping/adjacent 1-based inclusive intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def complement_within(intervals: list, chrom_size: int) -> list:
    """Return gaps in the 1-based inclusive interval list within [1, chrom_size]."""
    if not intervals:
        return [(1, chrom_size)]
    intervals = sorted(intervals)
    gaps = []
    prev_end = 0
    for s, e in intervals:
        if s > prev_end + 1:
            gaps.append((prev_end + 1, s - 1))
        prev_end = max(prev_end, e)
    if prev_end < chrom_size:
        gaps.append((prev_end + 1, chrom_size))
    return gaps


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('gtf', type=Path)
    p.add_argument('fai', type=Path, help='reference FASTA index (.fai)')
    p.add_argument('output_bed', type=Path)
    p.add_argument('--flank', type=int, default=5000, help='bp to extend around each gene body (default: 5000)')
    p.add_argument('--min-size', type=int, default=10000, help='drop intergenic regions smaller than this (default: 10000)')
    args = p.parse_args()

    for f in (args.gtf, args.fai):
        if not f.is_file():
            print(f'Error: not found: {f}', file=sys.stderr)
            return 1

    chrom_sizes = load_chrom_sizes(args.fai)
    gene_bodies = load_gene_bodies(args.gtf)

    args.output_bed.parent.mkdir(parents=True, exist_ok=True)
    total_bp = 0
    n_regions = 0
    with args.output_bed.open('w') as out:
        for chrom in sorted(chrom_sizes.keys()):
            csize = chrom_sizes[chrom]
            expanded = [(max(1, s - args.flank), min(csize, e + args.flank))
                        for s, e in gene_bodies.get(chrom, [])]
            expanded = merge_intervals(expanded)
            for s, e in complement_within(expanded, csize):
                width = e - s + 1
                if width >= args.min_size:
                    # BED is 0-based half-open.
                    out.write(f'{chrom}\t{s - 1}\t{e}\n')
                    total_bp += width
                    n_regions += 1

    print(f'Wrote {args.output_bed}')
    print(f'  flank:            {args.flank} bp')
    print(f'  min_size:         {args.min_size} bp')
    print(f'  regions:          {n_regions}')
    print(f'  total intergenic: {total_bp:,} bp')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
