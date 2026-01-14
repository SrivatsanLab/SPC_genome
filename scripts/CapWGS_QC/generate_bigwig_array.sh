#!/bin/bash
#SBATCH -J bw_gen
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 2:00:00

###########################################################################################################################
# Generate bigwig files from BAM files using deepTools bamCoverage
###########################################################################################################################

set -euo pipefail

# Arguments
SAMPLE_LIST="$1"      # File containing sample names (one per line)
BAM_DIR="$2"          # Directory containing BAM files
OUTPUT_DIR="$3"       # Output directory for bigwig files

# Get sample name from array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# Input and output files
BAM_FILE="${BAM_DIR}/${SAMPLE}.bam"
BIGWIG_FILE="${OUTPUT_DIR}/${SAMPLE}.bw"

echo "Processing sample: ${SAMPLE}"
echo "BAM file: ${BAM_FILE}"
echo "Output bigwig: ${BIGWIG_FILE}"

# Check if BAM file exists
if [ ! -f "${BAM_FILE}" ]; then
    echo "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

# Load deepTools module
module load deepTools

# Generate bigwig file with 1000bp bins
bamCoverage -b "${BAM_FILE}" -o "${BIGWIG_FILE}" --binSize 1000 -of bigwig -p 4

echo "Bigwig generation complete for ${SAMPLE}"
