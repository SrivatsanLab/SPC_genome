#!/bin/bash
#SBATCH -J extract_sc_bam
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH -c 1
#SBATCH --mem=4G

##########################################################################################################################
# Extract reads for a single cell from preprocessed bulk BAM
# This is an array job - each task extracts one cell
##########################################################################################################################

set -euo pipefail

module load SAMtools

BULK_BAM="$1"           # Input bulk BAM file
BARCODE_FILE="$2"       # File with list of cell barcodes
SC_OUTPUTS_DIR="$3"     # Output directory for single-cell BAMs

# Get the barcode for this array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

if [ -z "$barcode" ]; then
    echo "ERROR: Could not read barcode for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Output BAM file
output_bam="${SC_OUTPUTS_DIR}/${barcode}.bam"

echo "Extracting cell: ${barcode}"
echo "Input BAM: ${BULK_BAM}"
echo "Output BAM: ${output_bam}"

# Extract reads with this cell's barcode
# The CB:Z: tag contains the cell barcode
samtools view -h "${BULK_BAM}" | \
    grep -E "^@|CB:Z:${barcode}" | \
    samtools view -b -o "${output_bam}"

# Index the output BAM
samtools index "${output_bam}"

echo "Extraction complete for cell: ${barcode}"
