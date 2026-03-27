#!/bin/bash
#SBATCH -J preseq
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 2:00:00

set -euo pipefail

# Arguments
BAM_FILE="$1"
OUTPUT_DIR="$2"
EXTRAP_READS="${3:-100000000}"  # Default 100M if not specified
SAMPLE_NAME=$(basename "${BAM_FILE}" | sed 's/\.bam$//' | sed 's/_dna$//')

echo "Running preseq on sample: ${SAMPLE_NAME}"
echo "Extrapolating to ${EXTRAP_READS} reads"

# Run preseq lc_extrap (library complexity extrapolation)
# -B: input is BAM
# -e: extrapolate to this many reads
# -s: step size
# -v: verbose
preseq lc_extrap \
    -B \
    -e ${EXTRAP_READS} \
    -s 1000000 \
    -o "${OUTPUT_DIR}/${SAMPLE_NAME}_lc_extrap.txt" \
    -v \
    "${BAM_FILE}"

# Run preseq c_curve (complexity curve - actual observed)
# -B: input is BAM
# -o: output file
# -v: verbose
preseq c_curve \
    -B \
    -o "${OUTPUT_DIR}/${SAMPLE_NAME}_c_curve.txt" \
    -v \
    "${BAM_FILE}"

echo "Preseq complete for ${SAMPLE_NAME}"
