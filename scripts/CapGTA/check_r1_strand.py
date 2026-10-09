#!/usr/bin/env python3
"""For reads already featureCounts-Assigned, check R1 orientation vs gene strand.

If the polyT RT primer gives us first-strand cDNA (= antisense of mRNA), we
expect R1 to be on the OPPOSITE strand of the gene it was assigned to.
  - gene '+', R1 is_reverse=True   → antisense (expected if `-s 2`)
  - gene '-', R1 is_reverse=False  → antisense
  - gene '+', R1 is_reverse=False  → sense
  - gene '-', R1 is_reverse=True   → sense

Prints the breakdown and recommends -s 1 vs -s 2 for featureCounts.
"""
import argparse
from collections import Counter
from pathlib import Path
import pandas as pd
import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--anno-dir', required=True)
    ap.add_argument('--counts-raw', required=True)
    ap.add_argument('--n-cells', type=int, default=200)
    args = ap.parse_args()

    df = pd.read_csv(args.counts_raw, sep='\t', comment='#')
    df = df.rename(columns={'Geneid': 'gene_id'})
    strand_lookup = {row['gene_id']: str(row['Strand']).split(';')[0]
                     for _, row in df.iterrows()}
    print(f'Loaded strand info for {len(strand_lookup)} genes', flush=True)

    bams = sorted(Path(args.anno_dir).glob('*.bam.featureCounts.bam'))[:args.n_cells]
    print(f'Scanning {len(bams)} cells...', flush=True)

    counts = Counter()   # (gene_strand, r1_reverse) -> n
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
                if read.get_tag('XS') != 'Assigned': continue
                gene = read.get_tag('XT')
            except KeyError:
                continue
            s = strand_lookup.get(gene)
            if s is None: continue
            counts[(s, read.is_reverse)] += 1
        bam.close()
        if (i + 1) % 50 == 0:
            print(f'  [{i+1}/{len(bams)}] cumulative={sum(counts.values())}', flush=True)

    total = sum(counts.values())
    print(f'\n=== Total Assigned R1 reads examined: {total} ===')
    for s in ('+', '-'):
        for rev in (False, True):
            n = counts.get((s, rev), 0)
            label = f'gene={s}  R1.is_reverse={rev}'
            orient = ('sense' if (s == '+' and not rev) or (s == '-' and rev)
                      else 'antisense')
            print(f'  {label:<35s}: {n:>8d}  ({n/max(total,1):.1%})  [{orient}]')

    n_antisense = counts[('+', True)] + counts[('-', False)]
    n_sense     = counts[('+', False)] + counts[('-', True)]
    print(f'\nSummary:')
    print(f'  sense  (R1 matches gene strand):      {n_sense:>8d}  ({n_sense/max(total,1):.1%})')
    print(f'  antisense (R1 opposite gene strand):  {n_antisense:>8d}  ({n_antisense/max(total,1):.1%})')

    print(f'\n>>> Recommendation for featureCounts: ', end='')
    if n_antisense > 3 * n_sense:
        print('-s 2  (reverse-stranded; R1 is antisense to the mRNA)')
    elif n_sense > 3 * n_antisense:
        print('-s 1  (forward-stranded; R1 matches the mRNA sense)')
    else:
        print('UNSTRANDED (ratios too close — library may not be directional)')


if __name__ == '__main__':
    main()
