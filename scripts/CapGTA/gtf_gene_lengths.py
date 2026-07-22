#!/usr/bin/env python3
"""Compute per-gene lengths from a GTF file.

For each gene, collect exons from every isoform, merge overlapping intervals,
and sum the merged widths. This "union exonic length" matches what
featureCounts reports in the Length column when features are grouped by gene_id.

Output CSV columns:
    gene_id
    chrom
    gene_start        min exon start (1-based inclusive)
    gene_end          max exon end
    gene_span         gene_end - gene_start + 1 (includes introns)
    exonic_length     sum of merged exon widths (excludes introns)
    n_transcripts     number of annotated isoforms

Usage:
    gtf_gene_lengths.py <annotation.gtf> <output.csv>
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_attributes(attr_str: str) -> dict:
    """Parse a GTF attribute column: 'key "value"; key "value"; ...'."""
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


def merged_length(intervals: list[tuple[int, int]]) -> int:
    """Sum of widths after merging overlapping/adjacent 1-based inclusive intervals."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            if end > cur_end:
                cur_end = end
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    total += cur_end - cur_start + 1
    return total


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('gtf', type=Path, help='GTF annotation file')
    parser.add_argument('output_csv', type=Path, help='output CSV path')
    args = parser.parse_args()

    if not args.gtf.is_file():
        print(f"Error: GTF not found: {args.gtf}", file=sys.stderr)
        return 1

    gene_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    gene_transcripts: dict[str, set[str]] = defaultdict(set)
    gene_chrom: dict[str, str] = {}

    with args.gtf.open() as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            chrom, _, _, start, end, _, _, _, attrs = fields
            a = parse_attributes(attrs)
            gid = a.get('gene_id')
            tid = a.get('transcript_id')
            if gid is None:
                continue
            gene_exons[gid].append((int(start), int(end)))
            if tid is not None:
                gene_transcripts[gid].add(tid)
            gene_chrom.setdefault(gid, chrom)

    all_genes = sorted(gene_exons.keys())
    rows = []
    for gid in all_genes:
        exons = gene_exons[gid]
        starts = [s for s, _ in exons]
        ends = [e for _, e in exons]
        rows.append({
            'gene_id':       gid,
            'chrom':         gene_chrom.get(gid, ''),
            'gene_start':    min(starts),
            'gene_end':      max(ends),
            'gene_span':     max(ends) - min(starts) + 1,
            'exonic_length': merged_length(exons),
            'n_transcripts': len(gene_transcripts.get(gid, ())),
        })

    df = pd.DataFrame(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)

    print(f"Wrote {args.output_csv}")
    print(f"  genes:                    {len(df)}")
    print(f"  median exonic_length:     {int(df['exonic_length'].median())}")
    print(f"  median gene_span:         {int(df['gene_span'].median())}")
    print(f"  max exonic_length:        {int(df['exonic_length'].max())}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
