#!/usr/bin/env python
"""Per-cell false-positive rate vs matched bulk WGS truth.

For each cell in a donor's joint VCF, identify "cell calls" and query the
matched bulk BAM(s) for read evidence at those sites. Each call is classified:

  TP     : bulk DP > MIN_BULK_DP and bulk has >= 1 read of the same alt.
           (Single-read bulk support is treated as low-VAF somatic truth, not
           coincident sequencing error -- the joint probability of the bulk
           and a cell making the same error at the same base is negligible.
           n_TP_singleton tracks this subset for transparency.)
  FP     : bulk DP > MIN_BULK_DP and bulk has 0 reads of that alt
  Lowcov : bulk DP <= MIN_BULK_DP (call dropped from FPR denominator)

Two definitions of "cell call":
  --call-mode gt   (default) HC's GT contains the alt allele index
                   (i.e. genotype-level variant call from HaplotypeCaller).
  --call-mode ad   AD[alt] >= --min-alt (raw read-evidence floor; ignores
                   HC's quality model and is much more permissive).

SNV-only; indel candidates are counted but skipped.
Output: one TSV row per cell with TP/FP/ambig/lowcov counts and FP fraction.
"""
import argparse
import gzip
import sys
import time
from collections import defaultdict

import pysam


def is_snv(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and ref != alt and alt in 'ACGT'


def good_read(read, min_mq):
    return (read.mapping_quality >= min_mq
            and not read.is_unmapped
            and not read.is_duplicate
            and not read.is_secondary
            and not read.is_supplementary
            and not read.is_qcfail)


def cell_calls_alt(sample_data, alt_idx, mode, min_alt):
    """Return True if this cell calls this alt at this site under `mode`."""
    if mode == 'gt':
        gt = sample_data.get('GT')
        if gt is None:
            return False
        return (alt_idx + 1) in gt
    if mode == 'ad':
        ad = sample_data.get('AD')
        if ad is None or len(ad) <= alt_idx + 1:
            return False
        ac = ad[alt_idx + 1]
        return ac is not None and ac >= min_alt
    raise ValueError(f"unknown call mode: {mode}")


def collect_candidate_sites(joint_vcf, mode, min_alt, log):
    """Pass 1: stream joint VCF, collect unique SNV candidate sites
    (chrom, pos) -> set of alt bases that any cell calls under `mode`.
    """
    candidates = defaultdict(set)
    indel_counts = defaultdict(int)
    vcf = pysam.VariantFile(joint_vcf)
    samples = list(vcf.header.samples)
    n_records = 0
    t0 = time.time()
    for rec in vcf:
        n_records += 1
        if n_records % 1_000_000 == 0:
            log(f"  pass1: {n_records:,} records in {time.time()-t0:.0f}s")
        if rec.alts is None:
            continue
        for alt_idx, alt in enumerate(rec.alts):
            snv = is_snv(rec.ref, alt)
            for sample in samples:
                if not cell_calls_alt(rec.samples[sample], alt_idx, mode, min_alt):
                    continue
                if snv:
                    candidates[(rec.chrom, rec.pos)].add(alt)
                else:
                    indel_counts[sample] += 1
    log(f"  pass1 done: {n_records:,} records, "
        f"{len(candidates):,} unique SNV positions, "
        f"{len(samples)} cells, {time.time()-t0:.0f}s")
    return candidates, samples, indel_counts


def cluster_positions(positions, max_gap=20000):
    """Cluster sorted positions into (start, end) ranges where consecutive
    positions are within max_gap of each other. Avoids wasting BAM iteration
    over long stretches with no candidate sites.
    """
    if not positions:
        return []
    clusters = []
    start = prev = positions[0]
    for p in positions[1:]:
        if p - prev > max_gap:
            clusters.append((start, prev))
            start = p
        prev = p
    clusters.append((start, prev))
    return clusters


def query_bulk(bulk_bams, candidates, min_mq, min_bq, threads, log):
    """Streaming pileup per chromosome, summed across bulk BAMs.

    Groups candidates by chromosome, clusters them so we don't iterate over
    long empty stretches, then runs a single pileup per cluster per BAM.
    Returns dict (chrom, pos, alt) -> (alt_count, dp).
    """
    by_chrom = defaultdict(dict)
    for (chrom, pos), alts in candidates.items():
        by_chrom[chrom][pos] = alts

    handles = [pysam.AlignmentFile(b, 'rb', threads=threads) for b in bulk_bams]
    log(f"  opened {len(handles)} bulk BAM(s) with threads={threads}")

    result = {}
    t0 = time.time()
    n_done = 0
    n_total = sum(len(d) for d in by_chrom.values())

    for chrom in sorted(by_chrom.keys()):
        site_pos = by_chrom[chrom]
        positions = sorted(site_pos.keys())
        positions_set = set(positions)
        per_pos = {pos: {'A': 0, 'C': 0, 'G': 0, 'T': 0} for pos in positions}
        clusters = cluster_positions(positions, max_gap=20000)

        for h in handles:
            for (lo, hi) in clusters:
                for col in h.pileup(chrom, lo - 1, hi,
                                    truncate=True,
                                    min_base_quality=0,
                                    min_mapping_quality=min_mq,
                                    stepper='samtools'):
                    p = col.reference_pos + 1
                    if p not in positions_set:
                        continue
                    cnt = per_pos[p]
                    for read in col.pileups:
                        if read.is_del or read.is_refskip:
                            continue
                        qpos = read.query_position
                        if qpos is None:
                            continue
                        if read.alignment.query_qualities[qpos] < min_bq:
                            continue
                        base = read.alignment.query_sequence[qpos]
                        if base in 'ACGT':
                            cnt[base] += 1

        for pos in positions:
            cnt = per_pos[pos]
            dp = cnt['A'] + cnt['C'] + cnt['G'] + cnt['T']
            for alt in site_pos[pos]:
                result[(chrom, pos, alt)] = (cnt.get(alt, 0), dp)

        n_done += len(positions)
        log(f"    {chrom}: {len(positions):,} sites, "
            f"{len(clusters)} clusters, total {n_done:,}/{n_total:,} "
            f"in {time.time()-t0:.0f}s")

    for h in handles:
        h.close()
    log(f"  pileup done: {n_total:,} sites, {time.time()-t0:.0f}s")
    return result


def classify_per_cell(joint_vcf, bulk_evidence, mode, min_alt, min_bulk_dp, log,
                      donor=None, per_call_path=None):
    """Pass 2: stream joint VCF again, accumulate per-cell counts.

    If per_call_path is set, also stream a gzipped per-call TSV with one row
    per SNV cell call (TP/FP/lowcov) for downstream stratified analysis
    (e.g. plotting bulk_dp distributions).
    """
    counts = defaultdict(lambda: {
        'n_calls_total': 0,
        'n_snv_calls': 0,
        'n_indel_calls': 0,
        'n_evaluable': 0,
        'n_TP': 0,
        'n_FP': 0,
        'n_TP_singleton': 0,
        'n_lowcov': 0,
    })
    pc_handle = None
    if per_call_path:
        pc_handle = gzip.open(per_call_path, 'wt')
        pc_handle.write('\t'.join(['donor', 'sample', 'chrom', 'pos', 'alt',
                                   'bulk_dp', 'bulk_alt', 'class']) + '\n')
        log(f"  per-call TSV: {per_call_path}")
    vcf = pysam.VariantFile(joint_vcf)
    samples = list(vcf.header.samples)
    n_records = 0
    t0 = time.time()
    for rec in vcf:
        n_records += 1
        if n_records % 1_000_000 == 0:
            log(f"  pass2: {n_records:,} records in {time.time()-t0:.0f}s")
        if rec.alts is None:
            continue
        for alt_idx, alt in enumerate(rec.alts):
            snv = is_snv(rec.ref, alt)
            for sample in samples:
                if not cell_calls_alt(rec.samples[sample], alt_idx, mode, min_alt):
                    continue
                c = counts[sample]
                c['n_calls_total'] += 1
                if not snv:
                    c['n_indel_calls'] += 1
                    continue
                c['n_snv_calls'] += 1
                ev = bulk_evidence.get((rec.chrom, rec.pos, alt))
                if ev is None:
                    c['n_lowcov'] += 1
                    if pc_handle:
                        pc_handle.write(f'{donor}\t{sample}\t{rec.chrom}\t{rec.pos}\t{alt}\t-1\t-1\tno_pileup\n')
                    continue
                bulk_alt, bulk_dp = ev
                if bulk_dp <= min_bulk_dp:
                    c['n_lowcov'] += 1
                    if pc_handle:
                        pc_handle.write(f'{donor}\t{sample}\t{rec.chrom}\t{rec.pos}\t{alt}\t{bulk_dp}\t{bulk_alt}\tlowcov\n')
                    continue
                c['n_evaluable'] += 1
                if bulk_alt >= 1:
                    c['n_TP'] += 1
                    if bulk_alt == 1:
                        c['n_TP_singleton'] += 1
                        klass = 'TP_singleton'
                    else:
                        klass = 'TP'
                else:
                    c['n_FP'] += 1
                    klass = 'FP'
                if pc_handle:
                    pc_handle.write(f'{donor}\t{sample}\t{rec.chrom}\t{rec.pos}\t{alt}\t{bulk_dp}\t{bulk_alt}\t{klass}\n')
    if pc_handle:
        pc_handle.close()
    log(f"  pass2 done: {n_records:,} records, {time.time()-t0:.0f}s")
    return counts, samples


def write_tsv(out_path, donor, samples, counts):
    cols = ['donor', 'sample',
            'n_calls_total', 'n_snv_calls', 'n_indel_calls',
            'n_evaluable', 'n_TP', 'n_FP', 'n_TP_singleton', 'n_lowcov',
            'fp_fraction']
    with open(out_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for s in samples:
            c = counts.get(s, {})
            tp = c.get('n_TP', 0)
            fp = c.get('n_FP', 0)
            denom = tp + fp
            fp_frac = (fp / denom) if denom > 0 else float('nan')
            row = [donor, s,
                   c.get('n_calls_total', 0),
                   c.get('n_snv_calls', 0),
                   c.get('n_indel_calls', 0),
                   c.get('n_evaluable', 0),
                   tp, fp,
                   c.get('n_TP_singleton', 0),
                   c.get('n_lowcov', 0),
                   f'{fp_frac:.6f}']
            f.write('\t'.join(str(x) for x in row) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--donor', required=True)
    p.add_argument('--joint-vcf', required=True)
    p.add_argument('--bulk-bam', required=True, nargs='+',
                   help='One or more matched bulk BAMs (indexed). Counts are summed.')
    p.add_argument('--output', required=True)
    p.add_argument('--call-mode', choices=['gt', 'ad'], default='gt',
                   help='How to define a cell call: gt (HC genotype) | ad (raw AD>=min-alt)')
    p.add_argument('--min-alt', type=int, default=2,
                   help='Min alt-supporting reads in cell (only used in ad mode)')
    p.add_argument('--min-bulk-dp', type=int, default=10,
                   help='Min bulk DP to evaluate a candidate (strict: > MIN_BULK_DP)')
    p.add_argument('--min-mq', type=int, default=20)
    p.add_argument('--min-bq', type=int, default=20)
    p.add_argument('--threads', type=int, default=4,
                   help='BAM decompression threads (per BAM handle)')
    p.add_argument('--per-call-tsv', default=None,
                   help='Optional path to write a gzipped per-call TSV: one row per '
                        'SNV cell call with chrom/pos/alt/bulk_dp/bulk_alt/class. '
                        'Useful for DP-stratified analysis (plot bulk_dp distribution, '
                        'singleton-fraction by DP bin, etc).')
    args = p.parse_args()

    def log(msg):
        print(msg, file=sys.stderr, flush=True)

    log("=== per_cell_fpr.py ===")
    log(f"Donor:           {args.donor}")
    log(f"Joint VCF:       {args.joint_vcf}")
    log(f"Bulk BAM(s):     {args.bulk_bam}")
    log(f"Call mode:       {args.call_mode}"
        + (f" (min AD[alt] >= {args.min_alt})" if args.call_mode == 'ad' else ""))
    log(f"Min bulk DP:     >{args.min_bulk_dp}")
    log(f"Min MQ / BQ:     {args.min_mq} / {args.min_bq}")
    log(f"Output:          {args.output}")
    log("")

    log("Step 1: collecting unique SNV candidate sites...")
    candidates, samples, indel_counts = collect_candidate_sites(
        args.joint_vcf, mode=args.call_mode, min_alt=args.min_alt, log=log)
    n_unique_alts = sum(len(v) for v in candidates.values())
    log(f"  {len(candidates):,} positions, {n_unique_alts:,} unique alt alleles")

    log("Step 2: querying bulk BAM(s) at SNV candidates...")
    bulk_evidence = query_bulk(args.bulk_bam, candidates,
                                min_mq=args.min_mq, min_bq=args.min_bq,
                                threads=args.threads, log=log)
    log(f"  {len(bulk_evidence):,} (chrom, pos, alt) entries")

    log("Step 3: classifying per-cell calls...")
    counts, samples = classify_per_cell(
        args.joint_vcf, bulk_evidence,
        mode=args.call_mode, min_alt=args.min_alt,
        min_bulk_dp=args.min_bulk_dp, log=log,
        donor=args.donor, per_call_path=args.per_call_tsv)

    log(f"Step 4: writing TSV to {args.output}")
    write_tsv(args.output, args.donor, samples, counts)

    log("")
    log("=== Per-cell summary ===")
    log(f"{'sample':<50} {'n_call':>10} {'n_eval':>10} {'TP':>10} {'FP':>10} {'FP%':>7}")
    for s in samples:
        c = counts.get(s, {})
        tp, fp = c.get('n_TP', 0), c.get('n_FP', 0)
        denom = tp + fp
        fp_pct = (fp / denom * 100) if denom > 0 else float('nan')
        log(f"{s:<50} {c.get('n_calls_total', 0):>10,} "
            f"{c.get('n_evaluable', 0):>10,} {tp:>10,} {fp:>10,} {fp_pct:>6.2f}%")
    log("")
    log("Done.")


if __name__ == '__main__':
    main()
