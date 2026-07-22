#!/usr/bin/env python3
"""Enumerate splice junctions per gene from a GTF file.

For each transcript, sort its exons and treat the gap between consecutive
exons as one intron / one splice junction. Junctions shared across isoforms
of the same gene (same donor+acceptor coordinates) are counted once.

Output CSV columns:
    gene_id
    chrom
    n_junctions       unique introns for the gene (0 for single-exon genes)
    n_transcripts     number of annotated isoforms
    n_exons_max       exon count of the largest isoform

Usage:
    gtf_gene_junctions.py <annotation.gtf> <output.csv>
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

    transcript_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    transcript_gene: dict[str, str] = {}
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
            tid = a.get('transcript_id')
            gid = a.get('gene_id')
            if tid is None or gid is None:
                continue
            transcript_exons[tid].append((int(start), int(end)))
            transcript_gene[tid] = gid
            gene_chrom.setdefault(gid, chrom)

    gene_introns: dict[str, set[tuple[int, int]]] = defaultdict(set)
    gene_transcripts: dict[str, set[str]] = defaultdict(set)
    gene_max_exons: dict[str, int] = defaultdict(int)

    for tid, exons in transcript_exons.items():
        gid = transcript_gene[tid]
        gene_transcripts[gid].add(tid)
        exons_sorted = sorted(exons)
        gene_max_exons[gid] = max(gene_max_exons[gid], len(exons_sorted))
        for (_, e1_end), (e2_start, _) in zip(exons_sorted, exons_sorted[1:]):
            # Intron span in 1-based inclusive: (e1_end + 1) .. (e2_start - 1)
            gene_introns[gid].add((e1_end + 1, e2_start - 1))

    all_genes = sorted(gene_transcripts.keys())
    df = pd.DataFrame({
        'gene_id':       all_genes,
        'chrom':         [gene_chrom.get(g, '') for g in all_genes],
        'n_junctions':   [len(gene_introns.get(g, ())) for g in all_genes],
        'n_transcripts': [len(gene_transcripts.get(g, ())) for g in all_genes],
        'n_exons_max':   [gene_max_exons.get(g, 0) for g in all_genes],
    })

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)

    print(f"Wrote {args.output_csv}")
    print(f"  genes:                 {len(df)}")
    print(f"  single-exon (0 junc):  {(df['n_junctions'] == 0).sum()}")
    print(f"  median junctions/gene: {df['n_junctions'].median():.0f}")
    print(f"  max junctions/gene:    {df['n_junctions'].max()}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
