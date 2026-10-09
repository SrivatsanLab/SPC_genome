#!/usr/bin/env python3
"""Count genomic polyA/polyT tracts and overlap each XP's reads with them.

Finds all polyA or polyT runs >= MIN_RUN on the + strand of the C. elegans
genome, builds a sorted interval table per chromosome, then for each
primary R1 read in each fc-annotated BAM: records whether the alignment
overlaps a polyA tract (+-strand run of A's) or a polyT tract.

Report per XP:
  n_reads, n_overlap_polyA_tract, n_overlap_polyT_tract, n_overlap_either
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


def _strip_umi(name: str) -> str:
    return UMI_IN_NAME_PAT.sub('', name)


def find_tracts(fasta_path: str, min_run: int) -> tuple[dict, dict, dict]:
    """Returns ({chrom: (starts, ends)} for polyA tracts, same for polyT, chrom lengths).
       Coordinates are 1-based inclusive on the + strand."""
    pat_A = re.compile(rb'A{' + str(min_run).encode() + rb',}')
    pat_T = re.compile(rb'T{' + str(min_run).encode() + rb',}')
    polyA = {}; polyT = {}; chrom_lens = {}
    cur_chrom = None; cur_buf = []
    def flush():
        if cur_chrom is None: return
        seq = b''.join(cur_buf).upper()
        chrom_lens[cur_chrom] = len(seq)
        a_runs = [(m.start() + 1, m.end()) for m in pat_A.finditer(seq)]  # 1-based inclusive end
        t_runs = [(m.start() + 1, m.end()) for m in pat_T.finditer(seq)]
        polyA[cur_chrom] = (np.array([r[0] for r in a_runs], dtype=np.int64),
                            np.array([r[1] for r in a_runs], dtype=np.int64))
        polyT[cur_chrom] = (np.array([r[0] for r in t_runs], dtype=np.int64),
                            np.array([r[1] for r in t_runs], dtype=np.int64))
    with open(fasta_path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                flush()
                cur_chrom = line[1:].strip().split()[0].decode()
                cur_buf = []
            else:
                cur_buf.append(line.strip())
        flush()
    return polyA, polyT, chrom_lens


def _overlaps(starts: np.ndarray, ends: np.ndarray, lo: int, hi: int) -> bool:
    """Does the half-open [lo, hi) overlap any closed [starts, ends]?"""
    if starts.size == 0: return False
    # First tract with start < hi
    idx = np.searchsorted(starts, hi, side='right')
    if idx == 0: return False
    # Any of indices 0..idx-1 whose end >= lo?
    return bool((ends[:idx] >= lo).any())


def read_xp_lookup(fastq_path: Path) -> dict:
    xp = {}
    if not fastq_path.exists(): return xp
    with gzip.open(fastq_path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0: continue
            m = XP_PAT.search(line)
            if m is None: continue
            name = line.split()[0][1:]
            xp[name] = m.group(1)
    return xp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--fastq-root', required=True)
    ap.add_argument('--cells-list', required=True)
    ap.add_argument('--anno-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--min-run', type=int, default=15)
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    print(f'Finding polyA/polyT tracts >= {args.min_run} bp in genome...', flush=True)
    polyA, polyT, chrom_lens = find_tracts(args.fasta, args.min_run)
    tot_A = sum(len(v[0]) for v in polyA.values())
    tot_T = sum(len(v[0]) for v in polyT.values())
    tot_bp = sum(chrom_lens.values())
    tot_A_bp = sum(int((v[1] - v[0] + 1).sum()) for v in polyA.values())
    tot_T_bp = sum(int((v[1] - v[0] + 1).sum()) for v in polyT.values())
    print(f'\n=== Genome polyA/polyT tract summary (min run = {args.min_run} bp) ===')
    print(f'  Chromosomes: {len(chrom_lens)}  total bp: {tot_bp:,}')
    print(f'  polyA tracts (+ strand runs of A):  {tot_A:,}  '
          f'total {tot_A_bp:,} bp  ({tot_A_bp/tot_bp:.3%} of genome)')
    print(f'  polyT tracts (+ strand runs of T):  {tot_T:,}  '
          f'total {tot_T_bp:,} bp  ({tot_T_bp/tot_bp:.3%} of genome)')
    print(f'  ratio polyA/polyT: {tot_A/max(tot_T,1):.3f} by count, '
          f'{tot_A_bp/max(tot_T_bp,1):.3f} by bp')

    # Walk BAMs + FASTQs
    cells = [l.strip() for l in open(args.cells_list) if l.strip()]
    anno = Path(args.anno_dir); fqroot = Path(args.fastq_root)
    counts: dict = defaultdict(Counter)  # xp -> Counter with keys 'total','polyA','polyT','either'
    n_done = 0
    for cell in cells:
        udi = '_'.join(cell.split('_')[-2:])
        bam_path = anno / f'{cell}.bam.featureCounts.bam'
        if not bam_path.exists(): continue
        fq_r1 = fqroot / udi / f'{cell}_R1.polyAT_trimmed.fastq.gz'
        xp_lu = read_xp_lookup(fq_r1)
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception:
            continue
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
            if not read.is_read1: continue  # count per pair
            raw_name = _strip_umi(read.query_name)
            xp = xp_lu.get(raw_name)
            if xp is None: continue
            counts[xp]['total'] += 1
            try:
                if read.get_tag('XS') == 'Assigned':
                    counts[xp]['assigned'] += 1
            except KeyError:
                pass
            chrom = read.reference_name
            a_arr = polyA.get(chrom); t_arr = polyT.get(chrom)
            r_lo = read.reference_start + 1
            r_hi = read.reference_end
            has_A = a_arr is not None and _overlaps(a_arr[0], a_arr[1], r_lo, r_hi)
            has_T = t_arr is not None and _overlaps(t_arr[0], t_arr[1], r_lo, r_hi)
            if has_A: counts[xp]['polyA'] += 1
            if has_T: counts[xp]['polyT'] += 1
            if has_A or has_T: counts[xp]['either'] += 1
        bam.close()
        n_done += 1
        if n_done % 100 == 0:
            print(f'  [{n_done}] cumulative:', {k: sum(v.values()) for k, v in counts.items()}, flush=True)

    rows = []
    for xp in sorted(counts):
        d = counts[xp]
        total = d.get('total', 0)
        if total == 0: continue
        rows.append({
            'xp': xp,
            'n_r1_reads': total,
            'n_assigned': d.get('assigned', 0),
            'n_overlap_polyA_tract': d.get('polyA', 0),
            'n_overlap_polyT_tract': d.get('polyT', 0),
            'n_overlap_either': d.get('either', 0),
            'frac_polyA': d.get('polyA', 0) / total,
            'frac_polyT': d.get('polyT', 0) / total,
            'frac_either': d.get('either', 0) / total,
        })
    df = pd.DataFrame(rows)
    df.to_csv(out / 'xp_tract_overlap.tsv', sep='\t', index=False)
    print('\n=== Per-XP overlap with genomic polyA/polyT tracts ===')
    print(df.to_string(index=False))


if __name__ == '__main__':
    raise SystemExit(main())
