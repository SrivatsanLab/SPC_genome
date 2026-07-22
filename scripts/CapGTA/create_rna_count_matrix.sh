#!/bin/bash
#SBATCH -J rna_counts
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# RNA count matrix (genes x cells) from a list of per-cell BAMs.
#
# Runs featureCounts once with all BAMs as inputs (single job, multi-threaded) and reformats the
# output into a clean `gene_id x barcode` CSV. General-purpose: the caller decides which BAMs go
# in (e.g. spliced-only for coassay data — see filter_spliced_reads_array.sh).
#
# Usage:
#   sbatch scripts/CapGTA/create_rna_count_matrix.sh <bam_list.txt> <annotation.gtf> <output_prefix>
#
# Writes:
#   <output_prefix>_raw.txt         featureCounts raw output (with summary in *.summary)
#   <output_prefix>_matrix.csv      gene_id x barcode counts
#   <output_prefix>_summary.csv     per-cell total counts + genes detected
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_list.txt> <annotation.gtf> <output_prefix>" >&2
    exit 1
fi

BAM_LIST="$1"
GTF="$2"
OUTPUT_PREFIX="$3"

for f in "${BAM_LIST}" "${GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "${OUTPUT_PREFIX}")"

N_BAMS=$(wc -l < "${BAM_LIST}")

echo "=========================================="
echo "RNA count matrix"
echo "=========================================="
echo "BAM list:       ${BAM_LIST} (${N_BAMS} BAMs)"
echo "GTF:            ${GTF}"
echo "Output prefix:  ${OUTPUT_PREFIX}"
echo "Started:        $(date -Iseconds)"
echo "=========================================="

module load Subread

# Read BAM paths into an array for featureCounts invocation.
mapfile -t BAM_PATHS < "${BAM_LIST}"

# featureCounts args:
#   -a GTF, -o output
#   -T 16 threads
#   -t exon, -g gene_id: count reads overlapping exons, aggregate by gene
#   --fracOverlap 0.5:  require 50% of read length to overlap exon
#   --primary:          only primary alignments
#   -p:                 paired-end input (BAMs are PE from spc-align)
#   -B:                 count only pairs where both ends mapped
#   -C:                 exclude chimeric pairs (different chromosomes)
# Note: no -M (multi-mapper), because STAR was run with --outFilterMultimapNmax 1 upstream.
featureCounts \
    -a "${GTF}" \
    -o "${OUTPUT_PREFIX}_raw.txt" \
    -T 16 \
    -t exon \
    -g gene_id \
    --fracOverlap 0.5 \
    --primary \
    -p \
    -B \
    -C \
    "${BAM_PATHS[@]}"

echo ""
echo "featureCounts complete — reformatting into matrix..."

python3 - "${OUTPUT_PREFIX}_raw.txt" "${OUTPUT_PREFIX}" <<'PYTHON_SCRIPT'
import os
import sys

import pandas as pd

raw_path = sys.argv[1]
out_prefix = sys.argv[2]

counts = pd.read_csv(raw_path, sep='\t', comment='#')

# featureCounts columns: Geneid, Chr, Start, End, Strand, Length, <bam_path>, <bam_path>, ...
gene_ids = counts['Geneid']
bam_cols = counts.columns[6:]

# Barcode = BAM basename without the .bam suffix.
barcodes = [os.path.basename(c).removesuffix('.bam') for c in bam_cols]

matrix = counts[bam_cols].copy()
matrix.columns = barcodes
matrix.insert(0, 'gene_id', gene_ids)

matrix_csv = f'{out_prefix}_matrix.csv'
matrix.to_csv(matrix_csv, index=False)
print(f'Wrote {matrix_csv}  ({len(gene_ids)} genes x {len(barcodes)} cells)')

per_cell_totals = matrix[barcodes].sum(axis=0)
per_cell_detected = (matrix[barcodes] > 0).sum(axis=0)

summary = pd.DataFrame({
    'barcode': barcodes,
    'total_counts': per_cell_totals.values,
    'genes_detected': per_cell_detected.values,
})
summary_csv = f'{out_prefix}_summary.csv'
summary.to_csv(summary_csv, index=False)
print(f'Wrote {summary_csv}')

print()
print('Per-cell summary:')
print(f'  Mean counts/cell:   {per_cell_totals.mean():.0f}')
print(f'  Median counts/cell: {per_cell_totals.median():.0f}')
print(f'  Mean genes/cell:    {per_cell_detected.mean():.0f}')
print(f'  Median genes/cell:  {per_cell_detected.median():.0f}')
PYTHON_SCRIPT

echo ""
echo "Finished: $(date -Iseconds)"
