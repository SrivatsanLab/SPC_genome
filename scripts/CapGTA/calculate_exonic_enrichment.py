#!/usr/bin/env python3
"""
Calculate exonic enrichment for a BAM file.

Exonic Enrichment = (Exonic Reads / Total Reads) - (Exonic Bases / Total Reference Bases)

A positive value indicates enrichment for exonic regions (expected for RNA).
A value near zero indicates no enrichment (expected for DNA).

Usage:
    python calculate_exonic_enrichment.py <bam_file> <gtf_file> <reference_fasta_or_dir>
"""

import sys
import os
import pysam
import argparse


def parse_gtf_exons(gtf_file):
    """
    Parse GTF file and extract exon coordinates.

    Returns:
        dict: {chrom: [(start, end), ...]} - exon intervals by chromosome
        int: total exonic bases (non-overlapping)
    """
    from collections import defaultdict

    exons_by_chrom = defaultdict(list)

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'exon':
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # GTF is 1-based, convert to 0-based
            end = int(fields[4])  # GTF end is inclusive, keep as is for 0-based half-open

            exons_by_chrom[chrom].append((start, end))

    # Merge overlapping exons per chromosome to get total unique exonic bases
    total_exonic_bases = 0
    for chrom in exons_by_chrom:
        # Sort intervals
        intervals = sorted(exons_by_chrom[chrom])

        # Merge overlapping intervals
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1]:
                # Overlapping - extend the previous interval
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                # Non-overlapping - add new interval
                merged.append((start, end))

        exons_by_chrom[chrom] = merged

        # Sum up merged interval lengths
        for start, end in merged:
            total_exonic_bases += (end - start)

    return exons_by_chrom, total_exonic_bases


def get_genome_size(reference_input):
    """
    Get total genome size from reference FASTA or directory.

    Args:
        reference_input: Path to FASTA file or directory containing FASTA

    Returns:
        int: Total genome size in bases
    """
    # Determine if input is file or directory
    if os.path.isfile(reference_input):
        fasta_path = reference_input
    elif os.path.isdir(reference_input):
        # Search for FASTA file in standard locations
        search_paths = [
            os.path.join(reference_input, 'BWAIndex', 'genome.fa'),
            os.path.join(reference_input, 'genome.fa'),
        ]

        # Also search for any .fa or .fna files
        for subdir in ['BWAIndex', '']:
            search_dir = os.path.join(reference_input, subdir) if subdir else reference_input
            if os.path.isdir(search_dir):
                for fname in os.listdir(search_dir):
                    if fname.endswith('.fa') or fname.endswith('.fna'):
                        search_paths.append(os.path.join(search_dir, fname))

        # Find first existing FASTA
        fasta_path = None
        for path in search_paths:
            if os.path.isfile(path):
                fasta_path = path
                break

        if fasta_path is None:
            raise FileNotFoundError(f"Could not find FASTA file in {reference_input}")
    else:
        raise ValueError(f"Reference input is neither file nor directory: {reference_input}")

    # Check if FASTA index exists
    fai_path = fasta_path + '.fai'
    if not os.path.isfile(fai_path):
        print(f"Warning: FASTA index not found at {fai_path}, creating...", file=sys.stderr)
        pysam.faidx(fasta_path)

    # Read genome size from FASTA index
    total_size = 0
    with open(fai_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom_length = int(fields[1])
            total_size += chrom_length

    return total_size


def count_exonic_reads(bam_file, exons_by_chrom):
    """
    Count reads overlapping exons.

    Args:
        bam_file: Path to BAM file
        exons_by_chrom: dict of {chrom: [(start, end), ...]}

    Returns:
        tuple: (total_reads, exonic_reads)
    """
    total_reads = 0
    exonic_reads = 0

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch():
            # Skip secondary and supplementary alignments
            if read.is_secondary or read.is_supplementary:
                continue

            total_reads += 1

            # Check if read overlaps any exon on its chromosome
            chrom = read.reference_name
            if chrom not in exons_by_chrom:
                continue

            read_start = read.reference_start
            read_end = read.reference_end

            # Check for overlap with any exon
            for exon_start, exon_end in exons_by_chrom[chrom]:
                # Check if intervals overlap
                if read_start < exon_end and read_end > exon_start:
                    exonic_reads += 1
                    break  # Count read only once even if it overlaps multiple exons

    return total_reads, exonic_reads


def calculate_enrichment(bam_file, gtf_file, reference_input):
    """
    Calculate exonic enrichment for a BAM file.

    Returns:
        dict: {
            'total_reads': int,
            'exonic_reads': int,
            'exonic_read_fraction': float,
            'total_genome_size': int,
            'total_exonic_bases': int,
            'exonic_base_fraction': float,
            'exonic_enrichment': float
        }
    """
    # Parse GTF to get exon coordinates and total exonic bases
    print(f"Parsing GTF file: {gtf_file}", file=sys.stderr)
    exons_by_chrom, total_exonic_bases = parse_gtf_exons(gtf_file)

    # Get genome size
    print(f"Getting genome size from: {reference_input}", file=sys.stderr)
    total_genome_size = get_genome_size(reference_input)

    # Count exonic reads
    print(f"Counting exonic reads in: {bam_file}", file=sys.stderr)
    total_reads, exonic_reads = count_exonic_reads(bam_file, exons_by_chrom)

    # Calculate fractions
    exonic_read_fraction = exonic_reads / total_reads if total_reads > 0 else 0
    exonic_base_fraction = total_exonic_bases / total_genome_size

    # Calculate enrichment
    exonic_enrichment = exonic_read_fraction - exonic_base_fraction

    return {
        'total_reads': total_reads,
        'exonic_reads': exonic_reads,
        'exonic_read_fraction': exonic_read_fraction,
        'total_genome_size': total_genome_size,
        'total_exonic_bases': total_exonic_bases,
        'exonic_base_fraction': exonic_base_fraction,
        'exonic_enrichment': exonic_enrichment
    }


def main():
    parser = argparse.ArgumentParser(
        description='Calculate exonic enrichment for a BAM file'
    )
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('gtf_file', help='GTF annotation file')
    parser.add_argument('reference', help='Reference FASTA file or directory')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print detailed statistics to stderr')

    args = parser.parse_args()

    # Calculate enrichment
    result = calculate_enrichment(args.bam_file, args.gtf_file, args.reference)

    # Print results
    if args.verbose:
        print("", file=sys.stderr)
        print("Exonic Enrichment Analysis", file=sys.stderr)
        print("=" * 50, file=sys.stderr)
        print(f"Total reads:              {result['total_reads']:,}", file=sys.stderr)
        print(f"Exonic reads:             {result['exonic_reads']:,}", file=sys.stderr)
        print(f"Exonic read fraction:     {result['exonic_read_fraction']:.4f}", file=sys.stderr)
        print(f"Total genome size:        {result['total_genome_size']:,} bp", file=sys.stderr)
        print(f"Total exonic bases:       {result['total_exonic_bases']:,} bp", file=sys.stderr)
        print(f"Exonic base fraction:     {result['exonic_base_fraction']:.4f}", file=sys.stderr)
        print(f"Exonic enrichment:        {result['exonic_enrichment']:.4f}", file=sys.stderr)
        print("=" * 50, file=sys.stderr)
        print("", file=sys.stderr)

    # Print enrichment value to stdout (for easy capture)
    print(f"{result['exonic_enrichment']:.6f}")


if __name__ == '__main__':
    main()
