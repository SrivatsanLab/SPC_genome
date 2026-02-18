#!/bin/bash
#SBATCH -J extract_sc_bam
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -c 1
#SBATCH --mem=4G

##########################################################################################################################
# Extract reads for a single cell from a bulk BAM file
# This is a unified utility script used by all pipelines (CapWGS, CapGTA, etc.)
#
# Usage: sbatch --array=1-N extract_sc_from_bam_array.sh <bulk_bam> <barcode_file> <output_dir> [suffix]
#
# Arguments:
#   bulk_bam      - Input bulk BAM file (coordinate sorted, indexed)
#   barcode_file  - File with list of cell barcodes (one per line)
#   output_dir    - Output directory for single-cell BAMs
#   suffix        - Optional suffix for output filename (e.g., "_dna", "_rna")
#                   If provided, output is: ${barcode}${suffix}.bam
#                   If not provided, output is: ${barcode}.bam
#
# This script:
#   - Extracts reads with CB:Z tag matching the barcode
#   - Outputs coordinate-sorted BAM
#   - Indexes the output BAM
#
# This is an array job - each task extracts one cell based on SLURM_ARRAY_TASK_ID
##########################################################################################################################

set -euo pipefail

module load SAMtools

BULK_BAM="$1"           # Input bulk BAM file
BARCODE_FILE="$2"       # File with list of cell barcodes
SC_OUTPUTS_DIR="$3"     # Output directory for single-cell BAMs
SUFFIX="${4:-}"         # Optional suffix for output filename

# Get the barcode for this array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

if [ -z "$barcode" ]; then
    echo "ERROR: Could not read barcode for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Construct output BAM filename
if [ -n "$SUFFIX" ]; then
    output_bam="${SC_OUTPUTS_DIR}/${barcode}${SUFFIX}.bam"
else
    output_bam="${SC_OUTPUTS_DIR}/${barcode}.bam"
fi

echo "Extracting cell: ${barcode}"
echo "Input BAM: ${BULK_BAM}"
echo "Output BAM: ${output_bam}"
if [ -n "$SUFFIX" ]; then
    echo "Output suffix: ${SUFFIX}"
fi

# Extract reads with this cell's barcode and sort output
# The CB:Z: tag contains the cell barcode
# Using samtools sort instead of view -b ensures output is coordinate sorted
samtools view -h "${BULK_BAM}" | \
    grep -E "^@|CB:Z:${barcode}" | \
    samtools sort -o "${output_bam}"

# Index the output BAM
samtools index "${output_bam}"

echo "Extraction complete for cell: ${barcode}"
