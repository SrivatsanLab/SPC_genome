#!/bin/bash
#SBATCH -J sex_counts
#SBATCH -o sex_array_outs/%A_%a.out
#SBATCH -p short
#SBATCH -t 02:00:00
#SBATCH -c 8

mkdir -p sex_array_outs/

module load SAMtools

BARCODES=$1
DIR=$2

BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BARCODES}")
BAM_FILE="${DIR}/${BARCODE}.bam"

# Count reads aligned to the X chromosome, Y chromosome, and total reads
X_COUNT=$(samtools view -@ 8 -c "$BAM_FILE" "chrX")
Y_COUNT=$(samtools view -@ 8 -c "$BAM_FILE" "chrY")
TOTAL_COUNT=$(samtools view -@ 8 -c "$BAM_FILE")

# Write counts to a CSV file
echo "${BARCODE},${X_COUNT},${Y_COUNT},${TOTAL_COUNT}" > "$DIR/${BARCODE}_sex_counts.csv"