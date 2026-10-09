#!/usr/bin/env python3
"""Count pairs where BOTH mates sit mid-gene (far from both the 5' and 3' end).

For each pair:
    min_dist_any_end = min(|r1_off_5p|, |r1_off_3p|, |r2_off_5p|, |r2_off_3p|)
    using the nearest-transcript-3' end when available

Those with min_dist > threshold are candidates for genomic-hexamer-primed reads
that happened to carry a polyT/polyA motif the extractor mis-identified as a
primer tail. Overlap with genomic A/T tracts is reported to confirm."""
from __future__ import annotations
import re
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

OFFS = '/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline/qc/per_xp_full/r2a_mate_offsets_with_tx_nearest.csv'
COUNTS = '/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline/counts_raw.txt'
FASTA = '/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa'


def find_tracts(fasta: str, min_run: int) -> dict:
    pat_A = re.compile(rb'[Aa]{' + str(min_run).encode() + rb',}')
    pat_T = re.compile(rb'[Tt]{' + str(min_run).encode() + rb',}')
    out = {}
    cur, buf = None, []
    def flush():
        if cur is None: return
        seq = b''.join(buf)
        a = [(m.start()+1, m.end()) for m in pat_A.finditer(seq)]
        t = [(m.start()+1, m.end()) for m in pat_T.finditer(seq)]
        out[cur] = {
            'A': (np.array([x[0] for x in a], dtype=np.int64),
                  np.array([x[1] for x in a], dtype=np.int64)),
            'T': (np.array([x[0] for x in t], dtype=np.int64),
                  np.array([x[1] for x in t], dtype=np.int64)),
            'len': len(seq),
        }
    with open(fasta, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                flush()
                cur = line[1:].strip().split()[0].decode()
                buf = []
            else:
                buf.append(line.strip())
        flush()
    return out


def within(starts, ends, pos, win):
    if starts.size == 0: return False
    hi = pos + win
    lo = pos - win
    idx = np.searchsorted(starts, hi + 1, side='right')
    if idx == 0: return False
    return bool((ends[:idx] >= lo).any())


print('Finding tracts...', flush=True)
tracts = find_tracts(FASTA, 10)

print('Loading gene coords...', flush=True)
cr = pd.read_csv(COUNTS, sep='\t', comment='#', low_memory=False).rename(columns={'Geneid': 'gene_id'})
gene_meta = {}
for _, row in cr.iterrows():
    gs = min(int(s) for s in str(row['Start']).split(';'))
    ge = max(int(s) for s in str(row['End']).split(';'))
    strand = str(row['Strand']).split(';')[0]
    chrom = str(row['Chr']).split(';')[0]
    gene_meta[row['gene_id']] = (chrom, gs, ge, strand)

print('Loading offsets CSV...', flush=True)
df = pd.read_csv(OFFS, low_memory=False)
print(f'  {len(df):,} rows', flush=True)

# Minimum |offset to any gene end| across all 4 mate offsets
#   r1_off_5p, r1_off_3p, r2_off_5p, r2_off_3p
# Also compute 5' end offsets (gene body-based already in CSV).
# For the 3' end, use nearest-tx offsets when they're valid (not the sentinel ±1e7).
r1_3p_near = df.r1_off_nearest_3p.where(df.r1_off_nearest_3p.abs() < 1_000_000, df.r1_off_3p)
r2_3p_near = df.r2_off_nearest_3p.where(df.r2_off_nearest_3p.abs() < 1_000_000, df.r2_off_3p)

min_dist = np.minimum.reduce([
    df.r1_off_5p.abs().to_numpy(),
    df.r2_off_5p.abs().to_numpy(),
    r1_3p_near.abs().to_numpy(),
    r2_3p_near.abs().to_numpy(),
])

df['min_dist_any_end'] = min_dist

# Primer-side position for tract overlap check
XP_PRIMER = {'R1T': 'r1', 'R1A': 'r1', 'R2A': 'r2', 'R2T': 'r2'}
def primer_pos(row):
    xp = row['xp']
    side = XP_PRIMER.get(xp)
    if side is None: return None
    meta = gene_meta.get(row['gene_id'])
    if meta is None: return None
    chrom, gs, ge, strand = meta
    if strand == '+':
        return chrom, int(row[f'{side}_off_3p']) + ge
    else:
        return chrom, gs - int(row[f'{side}_off_3p'])

bins = [(0, 100, '0-100 bp'), (100, 200, '100-200'),
        (200, 500, '200-500'), (500, 1000, '500-1k'),
        (1000, 2000, '1k-2k'), (2000, 10_000_000, '>2k')]

print('\n=== Pairs by closest distance of any mate to any gene end (per XP) ===')
rows = []
for xp in ('R1T', 'R1A', 'R2A', 'R2T'):
    sub = df[df.xp == xp]
    n_total = len(sub)
    if n_total == 0: continue
    print(f'\n{xp} (n_assigned={n_total:,}):')
    for lo, hi, label in bins:
        in_bin = sub[(sub.min_dist_any_end >= lo) & (sub.min_dist_any_end < hi)]
        n_bin = len(in_bin)
        frac = n_bin / n_total
        # tract overlap of primer base
        hits = 0
        for _, row in in_bin.iterrows():
            r = primer_pos(row)
            if r is None: continue
            chrom, pos = r
            t = tracts.get(chrom)
            if t is None: continue
            if within(t['A'][0], t['A'][1], pos, 30) or within(t['T'][0], t['T'][1], pos, 30):
                hits += 1
        tract_frac = hits / max(n_bin, 1)
        print(f'  {label:<12s}: {n_bin:>8,}  ({frac:.2%})  '
              f'tract overlap of primer base: {hits:,} ({tract_frac:.2%})')
        rows.append({'xp': xp, 'bin': label, 'n': n_bin, 'frac': frac,
                     'tract_hits': hits, 'tract_frac': tract_frac})

pd.DataFrame(rows).to_csv(
    '/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline/qc/per_xp_full/midgene_distance_bins.tsv',
    sep='\t', index=False)
