#!/bin/bash
#SBATCH -J index_bam
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 1:00:00

set -euo pipefail

BARCODE_FILE="$1"
BAM_DIR="$2"

barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")
BAM_FILE="${BAM_DIR}/${barcode}.bam"

module load SAMtools/1.19.2-GCC-13.2.0
samtools index "${BAM_FILE}"

echo "Indexed ${barcode}.bam"
