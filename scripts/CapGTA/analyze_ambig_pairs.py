#!/usr/bin/env python3
"""For each Unassigned_Ambiguity R1 in the fc-annotated BAMs, find which gene
pair(s) it falls between. Tally pair frequency and summarize overlap types
(nested, antisense, head-to-head, head-to-tail, tandem)."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import pysam


def build_arrays(counts_raw: str, ext: int) -> dict:
    df = pd.read_csv(counts_raw, sep='\t', comment='#')
    df = df.rename(columns={'Geneid': 'gene_id'})
    rows_per_chrom: dict = defaultdict(list)
    strand_lookup: dict = {}
    for _, row in df.iterrows():
        starts = [int(s) for s in str(row['Start']).split(';')]
        ends   = [int(s) for s in str(row['End']).split(';')]
        chrom  = str(row['Chr']).split(';')[0]
        strand = str(row['Strand']).split(';')[0]
        gs = max(1, min(starts) - ext)
        ge = max(ends) + 1 + ext
        if ge > gs:
            rows_per_chrom[chrom].append((gs, ge, row['gene_id']))
            strand_lookup[row['gene_id']] = (strand,
                                             min(starts),
                                             max(ends))
    out = {}
    for chrom, rows in rows_per_chrom.items():
        rows.sort(key=lambda r: r[0])
        out[chrom] = (np.array([r[0] for r in rows], dtype=np.int64),
                      np.array([r[1] for r in rows], dtype=np.int64),
                      np.array([r[2] for r in rows], dtype=object))
    return out, strand_lookup


def overlapping(arr, lo, hi):
    starts, ends, ids = arr
    hi_idx = np.searchsorted(starts, hi, side='right')
    if hi_idx == 0:
        return []
    cand = np.arange(hi_idx)
    mask = ends[cand] > lo
    return [str(ids[i]) for i in cand[mask]]


def classify_pair(g1: str, g2: str, strand_lookup: dict) -> str:
    s1, start1, end1 = strand_lookup[g1]
    s2, start2, end2 = strand_lookup[g2]
    # Order so g1 starts first
    if start1 > start2:
        g1, g2 = g2, g1
        s1, s2 = s2, s1
        start1, start2 = start2, start1
        end1, end2 = end2, end1
    # Nested?
    if start2 >= start1 and end2 <= end1:
        return 'g2_nested_in_g1' + ('_antisense' if s1 != s2 else '_sense')
    if start1 >= start2 and end1 <= end2:
        return 'g1_nested_in_g2' + ('_antisense' if s1 != s2 else '_sense')
    # Partial overlap. Head/tail depends on strand.
    if s1 == s2:
        return 'tandem_same_strand'   # arranged in-line, same strand
    # Antisense partial overlap
    # if g1 is +, g2 is - and they overlap: g1 end abuts g2 end → tails face each other
    # if g1 is -, g2 is + and overlap: g1 start abuts g2 start → heads face each other
    if s1 == '+' and s2 == '-':
        return 'convergent_tails_antisense'   # 3' ends overlap
    if s1 == '-' and s2 == '+':
        return 'divergent_heads_antisense'    # 5' ends overlap (shared TSS region)
    return 'other'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--anno-dir', required=True)
    ap.add_argument('--counts-raw', required=True)
    ap.add_argument('--n-cells', type=int, default=200)
    ap.add_argument('--gene-extension', type=int, default=200,
                    help='Match the featureCounts --readExtension used')
    args = ap.parse_args()

    print(f'Building gene arrays with extension {args.gene_extension} bp...', flush=True)
    arrays, strand_lookup = build_arrays(args.counts_raw, args.gene_extension)
    print(f'  {sum(len(a[0]) for a in arrays.values())} genes', flush=True)

    bams = sorted(Path(args.anno_dir).glob('*.bam.featureCounts.bam'))[:args.n_cells]
    print(f'Scanning {len(bams)} cells...', flush=True)

    pair_counts: Counter = Counter()
    n_ambig_total = 0
    n_resolvable = 0
    n_ambig_with_mate_hint = 0

    for i, bam_path in enumerate(bams):
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception:
            continue
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped: continue
            if read.is_secondary or read.is_supplementary: continue
            if not read.is_read1: continue
            try:
                xs = read.get_tag('XS')
            except KeyError:
                continue
            if xs != 'Unassigned_Ambiguity':
                continue
            n_ambig_total += 1
            arr = arrays.get(read.reference_name)
            if arr is None: continue
            genes = overlapping(arr, read.reference_start + 1, read.reference_end)
            if len(genes) < 2:
                continue  # ambiguity resolved by extensions alone
            # Normalize pair as sorted tuple
            key = tuple(sorted(genes))
            pair_counts[key] += 1
            n_resolvable += 1
        bam.close()
        if (i + 1) % 50 == 0:
            print(f'  [{i+1}/{len(bams)}] ambig_r1={n_ambig_total} resolvable={n_resolvable}',
                  flush=True)

    print(f'\nTotal Ambiguity R1 reads scanned: {n_ambig_total}')
    print(f'  Resolvable (overlap ≥2 genes with current extension): {n_resolvable}')
    print(f'  Unique gene-pair combinations:                        {len(pair_counts)}')

    # Classify top gene pairs
    print(f'\n=== TOP 30 recurring gene-pair conflicts ===')
    print(f'{"count":>6} | {"geneA":<16} {"geneB":<16} | overlap type')
    for pair, c in pair_counts.most_common(30):
        if len(pair) != 2: continue
        if pair[0] not in strand_lookup or pair[1] not in strand_lookup: continue
        kind = classify_pair(pair[0], pair[1], strand_lookup)
        print(f'{c:>6} | {pair[0]:<16} {pair[1]:<16} | {kind}')

    # Aggregate classification across all pairs
    print(f'\n=== Classification summary (all 2-gene pairs, weighted by read count) ===')
    kind_counts: Counter = Counter()
    kind_pairs: Counter = Counter()
    for pair, c in pair_counts.items():
        if len(pair) != 2: continue
        if pair[0] not in strand_lookup or pair[1] not in strand_lookup: continue
        kind = classify_pair(pair[0], pair[1], strand_lookup)
        kind_counts[kind] += c
        kind_pairs[kind] += 1
    total_c = sum(kind_counts.values())
    total_p = sum(kind_pairs.values())
    for kind in sorted(kind_counts, key=lambda k: -kind_counts[k]):
        c = kind_counts[kind]
        p = kind_pairs[kind]
        print(f'  {kind:<36s}: {c:>6d} reads ({c/total_c:.1%})  across {p:>5d} unique pairs ({p/total_p:.1%})')

    print(f'\n=== >2-gene overlaps (triple+ ambiguity) ===')
    triples = sum(c for pair, c in pair_counts.items() if len(pair) >= 3)
    print(f'  {triples} reads in {sum(1 for p in pair_counts if len(p) >= 3)} gene triplets+')


if __name__ == '__main__':
    main()
