#!/usr/bin/env python3
"""Build per-gene local gDNA-baseline intervals and count reads over them.

Two baseline flavors:
    intron    per-gene merged intron span (single-exon genes → empty)
    flanking  ±window bp around gene body, minus overlap with ANY annotated
              exon/intron of any gene (i.e. clean intergenic in the
              gene's neighborhood)

For each flavor we emit a featureCounts-compatible SAF and then run
featureCounts on the (unfiltered) per-cell BAM list. Result: a table of
`gene_id × cell` counts for each baseline flavor, plus the total interval
length per gene (needed to convert counts → rate).

Usage:
    build_local_baselines.py \\
        --gtf <annotation.gtf> \\
        --bam-list <bam_list.txt> \\
        --output-dir <dir> \\
        [--flanking-window 5000] \\
        [--threads 8]

Outputs under <output_dir>:
    intron.saf, flanking.saf                  interval definitions
    intron_counts.tsv, flanking_counts.tsv    featureCounts raw output
    intron_lengths.tsv, flanking_lengths.tsv  gene_id, total_bp
"""

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def parse_gtf_exons(gtf_path: Path):
    """Yield (chrom, start, end, strand, gene_id) for every exon record."""
    with gtf_path.open() as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            attrs = {}
            for kv in f[8].strip().rstrip(';').split(';'):
                kv = kv.strip()
                if ' ' in kv:
                    k, v = kv.split(' ', 1)
                    attrs[k] = v.strip().strip('"')
            gid = attrs.get('gene_id')
            if gid is None:
                continue
            yield f[0], int(f[3]), int(f[4]), f[6], gid


def merge_intervals(ivs):
    """Merge overlapping [start,end] (1-based inclusive) intervals."""
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def subtract_intervals(a, b):
    """Return a - b (both lists of sorted, non-overlapping tuples)."""
    result = []
    bi = 0
    for s, e in a:
        cur_s = s
        while bi < len(b) and b[bi][1] < cur_s:
            bi += 1
        j = bi
        while j < len(b) and b[j][0] <= e:
            bs, be = b[j]
            if bs > cur_s:
                result.append((cur_s, min(bs - 1, e)))
            cur_s = max(cur_s, be + 1)
            if cur_s > e:
                break
            j += 1
        if cur_s <= e:
            result.append((cur_s, e))
    return result


def build_gene_intervals(gtf_path: Path, flanking_window: int):
    """Return per-gene (chrom, strand, exons_merged, introns, flanking)."""
    exons_by_gene = defaultdict(list)
    gene_chrom = {}
    gene_strand = {}
    exons_by_chrom = defaultdict(list)  # for masking flanking regions

    for chrom, s, e, strand, gid in parse_gtf_exons(gtf_path):
        exons_by_gene[gid].append((s, e))
        gene_chrom.setdefault(gid, chrom)
        gene_strand.setdefault(gid, strand)
        exons_by_chrom[chrom].append((s, e))

    # Global exon-mask per chrom (union across all genes)
    exon_mask = {c: merge_intervals(iv) for c, iv in exons_by_chrom.items()}

    out = {}
    for gid, ivs in exons_by_gene.items():
        merged = merge_intervals(ivs)
        gs = merged[0][0]
        ge = merged[-1][1]
        # Introns = gaps between merged exons
        introns = [(merged[i][1] + 1, merged[i + 1][0] - 1)
                   for i in range(len(merged) - 1)
                   if merged[i + 1][0] - merged[i][1] > 1]
        # Flanking = [gs-W, gs-1] ∪ [ge+1, ge+W], minus any exon on the chrom
        flank_raw = []
        if flanking_window > 0:
            flank_raw.append((max(1, gs - flanking_window), gs - 1))
            flank_raw.append((ge + 1, ge + flanking_window))
        flank_raw = [(a, b) for a, b in flank_raw if b >= a]
        flanking = subtract_intervals(flank_raw, exon_mask[gene_chrom[gid]])
        out[gid] = {
            'chrom':    gene_chrom[gid],
            'strand':   gene_strand[gid],
            'introns':  introns,
            'flanking': flanking,
        }
    return out


def write_saf(gene_intervals: dict, kind: str, out_path: Path):
    """Write a featureCounts SAF for one baseline flavor.

    SAF columns: GeneID Chr Start End Strand
    We emit one row per (gene, interval) so featureCounts sums intervals
    that share a GeneID.
    """
    with out_path.open('w') as fh:
        fh.write('GeneID\tChr\tStart\tEnd\tStrand\n')
        for gid, meta in gene_intervals.items():
            for s, e in meta[kind]:
                fh.write(f'{gid}\t{meta["chrom"]}\t{s}\t{e}\t{meta["strand"]}\n')


def write_lengths(gene_intervals: dict, kind: str, out_path: Path):
    with out_path.open('w') as fh:
        fh.write('gene_id\ttotal_bp\n')
        for gid, meta in gene_intervals.items():
            bp = sum(e - s + 1 for s, e in meta[kind])
            fh.write(f'{gid}\t{bp}\n')


def run_featurecounts(saf: Path, bam_list: Path, out_tsv: Path, threads: int):
    """Run featureCounts with the same PE-aware options as the main pipeline
    (see scripts/CapGTA/create_rna_count_matrix.sh) so counts are directly
    comparable across exon / intron / flanking passes.
    """
    with bam_list.open() as fh:
        bams = [ln.strip() for ln in fh if ln.strip()]
    cmd = [
        'featureCounts',
        '-F', 'SAF',
        '-a', str(saf),
        '-o', str(out_tsv),
        '-T', str(threads),
        '--fracOverlap', '0.5',
        '--primary',
        '-p',
        '-B',
        '-C',
        *bams,
    ]
    print('$', ' '.join(cmd[:14]), f'... ({len(bams)} BAMs)')
    subprocess.run(cmd, check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gtf', type=Path, required=True)
    ap.add_argument('--bam-list', type=Path, required=True,
                    help='one BAM path per line')
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--flanking-window', type=int, default=5000)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--skip-featurecounts', action='store_true',
                    help='write SAFs only; skip running featureCounts')
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Parsing GTF: {args.gtf}')
    gi = build_gene_intervals(args.gtf, args.flanking_window)
    n_genes = len(gi)
    n_with_introns = sum(1 for g in gi.values() if g['introns'])
    print(f'  genes:                {n_genes}')
    print(f'  with introns:         {n_with_introns}')
    print(f'  flanking window:      {args.flanking_window} bp')

    for kind in ('introns', 'flanking'):
        saf = args.output_dir / f'{kind.rstrip("s")}.saf'
        lens = args.output_dir / f'{kind.rstrip("s")}_lengths.tsv'
        write_saf(gi, kind, saf)
        write_lengths(gi, kind, lens)
        print(f'  wrote {saf}, {lens}')

    if args.skip_featurecounts:
        print('--skip-featurecounts set; done.')
        return 0

    for kind in ('intron', 'flanking'):
        saf = args.output_dir / f'{kind}.saf'
        out_tsv = args.output_dir / f'{kind}_counts.tsv'
        run_featurecounts(saf, args.bam_list, out_tsv, args.threads)

    return 0


if __name__ == '__main__':
    sys.exit(main())
