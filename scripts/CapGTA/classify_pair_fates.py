#!/usr/bin/env python3
"""Classify each pair's primer-read and mate-read into one of:
    can_map       - unmapped read
    intergenic    - mapped but no overlap with any gene body
    at_5prime     - mapped, overlaps a gene, within N bp of the gene 5' end
    at_3prime     - mapped, overlaps a gene, within N bp of nearest tx 3' end
    internal      - mapped, overlaps a gene, far from both ends

Prefers the gene-level info in XS/XT tags from featureCounts; falls back to
our own interval tree for intergenic positions.
"""
from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


UMI_IN_NAME_PAT = re.compile(r'_UMI_[ACGTN]+$')
XP_PAT = re.compile(r'\bXP:Z:(\w+)')


def _strip_umi(name): return UMI_IN_NAME_PAT.sub('', name)


def load_xp(fastq):
    xp = {}
    if not fastq.exists(): return xp
    with gzip.open(str(fastq), 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0: continue
            m = XP_PAT.search(line)
            if m is None: continue
            name = line.split()[0][1:]
            xp[name] = m.group(1)
    return xp


def load_gene_meta_and_tx(counts_raw: str, gtf: str) -> tuple[dict, dict]:
    """
    gene_meta[gene_id] = (chrom, min_start, max_end, strand)
    tx_3p[gene_id] = sorted np array of unique transcript 3' end coords
    chrom_genes[chrom] = (starts, ends) sorted arrays for intergenic check
    """
    cr = pd.read_csv(counts_raw, sep='\t', comment='#',
                     low_memory=False).rename(columns={'Geneid': 'gene_id'})
    gene_meta = {}
    chrom_genes = defaultdict(list)
    for _, row in cr.iterrows():
        gs = min(int(s) for s in str(row['Start']).split(';'))
        ge = max(int(s) for s in str(row['End']).split(';'))
        strand = str(row['Strand']).split(';')[0]
        chrom = str(row['Chr']).split(';')[0]
        gene_meta[row['gene_id']] = (chrom, gs, ge, strand)
        chrom_genes[chrom].append((gs, ge))
    chrom_arr = {}
    for c, rows in chrom_genes.items():
        rows.sort()
        chrom_arr[c] = (np.array([r[0] for r in rows], dtype=np.int64),
                        np.array([r[1] for r in rows], dtype=np.int64))

    tx_3p: dict = {}
    attr = re.compile(r'gene_id "([^"]+)"')
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'transcript': continue
            m = attr.search(p[8])
            if not m: continue
            gid = m.group(1)
            start = int(p[3]); end = int(p[4]); strand = p[6]
            tp = end if strand == '+' else start
            tx_3p.setdefault(gid, set()).add(tp)
    tx_3p = {g: np.array(sorted(v), dtype=np.int64) for g, v in tx_3p.items()}
    return gene_meta, tx_3p, chrom_arr


def overlaps_any_gene(chrom_arr: tuple, pos_lo: int, pos_hi: int) -> bool:
    starts, ends = chrom_arr
    hi = np.searchsorted(starts, pos_hi, side='right')
    if hi == 0: return False
    return bool((ends[:hi] >= pos_lo).any())


def classify(read, gene_meta: dict, tx_3p: dict, chrom_arr_lookup: dict,
             near_window: int) -> str:
    """Return one of: 'can_map', 'intergenic', 'at_5prime', 'at_3prime', 'internal'."""
    if read.is_unmapped:
        return 'can_map'
    # Check XS tag for a quick lookup
    try:
        xs = read.get_tag('XS')
    except KeyError:
        xs = None
    gene_id = None
    if xs == 'Assigned':
        try:
            gene_id = read.get_tag('XT')
        except KeyError:
            gene_id = None

    if gene_id is not None:
        meta = gene_meta.get(gene_id)
        if meta is None:
            return 'intergenic'
        chrom, gs, ge, strand = meta
        read_start_1 = read.reference_start + 1
        read_end_1   = read.reference_end
        if strand == '+':
            off_5p = read_start_1 - gs
            off_3p_body = read_end_1 - ge
        else:
            off_5p = ge - read_end_1
            off_3p_body = gs - read_start_1
        # Prefer nearest transcript 3' end
        tx = tx_3p.get(gene_id)
        if tx is not None and tx.size:
            if strand == '+':
                diffs = read_end_1 - tx
            else:
                diffs = tx - read_start_1
            off_3p = int(diffs[np.argmin(np.abs(diffs))])
        else:
            off_3p = off_3p_body
        if abs(off_3p) <= near_window:
            return 'at_3prime'
        if abs(off_5p) <= near_window:
            return 'at_5prime'
        return 'internal'

    # Not assigned - check if the read overlaps ANY gene body
    chrom = read.reference_name
    arr = chrom_arr_lookup.get(chrom)
    if arr is None:
        return 'intergenic'
    read_start_1 = read.reference_start + 1
    read_end_1 = read.reference_end
    if overlaps_any_gene(arr, read_start_1, read_end_1):
        return 'internal'  # overlaps a gene but fc rejected (maybe Ambiguity / Overlapping_Length)
    return 'intergenic'


XP_PRIMER_SIDE = {'R1T': 'r1', 'R1A': 'r1', 'R2A': 'r2', 'R2T': 'r2'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq-root', required=True)
    ap.add_argument('--cells-list', required=True)
    ap.add_argument('--counts-raw', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--anno-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--near-window', type=int, default=200)
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    print('Loading gene meta + transcript 3-prime ends...', flush=True)
    gene_meta, tx_3p, chrom_arr_lookup = load_gene_meta_and_tx(args.counts_raw, args.gtf)
    print(f'  {len(gene_meta):,} genes', flush=True)

    cells = [l.strip() for l in open(args.cells_list) if l.strip()]
    fqroot = Path(args.fastq_root); anno = Path(args.anno_dir)

    # xp -> Counter((primer_class, mate_class))
    counts = defaultdict(Counter)
    totals = defaultdict(int)
    cell_i = 0
    for cell in cells:
        udi = '_'.join(cell.split('_')[-2:])
        fq_r1 = fqroot / udi / f'{cell}_R1.polyAT_trimmed.fastq.gz'
        bam_path = anno / f'{cell}.bam.featureCounts.bam'
        if not bam_path.exists(): continue
        xp_lu = load_xp(fq_r1)
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception:
            continue
        pair_reads = {}
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary: continue
            raw = _strip_umi(read.query_name)
            xp = xp_lu.get(raw)
            if xp is None: continue
            d = pair_reads.setdefault(raw, {'xp': xp})
            if read.is_read1:
                d['r1_class'] = classify(read, gene_meta, tx_3p, chrom_arr_lookup, args.near_window)
            elif read.is_read2:
                d['r2_class'] = classify(read, gene_meta, tx_3p, chrom_arr_lookup, args.near_window)
        bam.close()

        for name, d in pair_reads.items():
            if 'r1_class' not in d or 'r2_class' not in d:
                continue  # one mate missing from BAM — shouldn't happen but skip
            xp = d['xp']
            side = XP_PRIMER_SIDE.get(xp)
            if side is None: continue
            primer = d['r1_class'] if side == 'r1' else d['r2_class']
            mate   = d['r2_class'] if side == 'r1' else d['r1_class']
            counts[xp][(primer, mate)] += 1
            totals[xp] += 1
        cell_i += 1
        if cell_i % 200 == 0:
            print(f'  [{cell_i}] totals: ' +
                  ' '.join(f'{xp}={totals[xp]:,}' for xp in ('R1T','R1A','R2A','R2T')),
                  flush=True)

    # Report
    classes = ['can_map', 'intergenic', 'at_5prime', 'internal', 'at_3prime']
    all_rows = []
    for xp in ('R1T', 'R1A', 'R2A', 'R2T'):
        tot = totals[xp]
        if tot == 0: continue
        print(f'\n=== XP={xp}  (n_pairs={tot:,}) ===')
        print(f'{"primer_class":<14s} {"mate_class":<14s} {"n_pairs":>12s} {"frac":>8s}')
        rows_for_xp = []
        for pr in classes:
            for mt in classes:
                n = counts[xp].get((pr, mt), 0)
                if n == 0 and pr != mt and pr not in ('internal', 'at_3prime'):
                    continue
                frac = n / tot
                print(f'{pr:<14s} {mt:<14s} {n:>12,d} {frac:>8.2%}')
                rows_for_xp.append({'xp': xp, 'primer': pr, 'mate': mt, 'n': n, 'frac': frac})
                all_rows.append({'xp': xp, 'primer': pr, 'mate': mt, 'n': n, 'frac': frac})

    pd.DataFrame(all_rows).to_csv(out / 'pair_fate_classification.tsv', sep='\t', index=False)
    print(f'\nWrote {out}/pair_fate_classification.tsv')


if __name__ == '__main__':
    raise SystemExit(main())
