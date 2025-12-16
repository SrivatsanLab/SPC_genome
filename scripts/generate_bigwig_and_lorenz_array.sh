#!/bin/bash
#SBATCH -J bw_lorenz
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 3:00:00

###########################################################################################################################
# Generate bigwig files and Lorenz curves from BAM files
# Outputs: bigwig and Lorenz curve CSV to same directory as BAM
###########################################################################################################################

set -euo pipefail

# Arguments
SAMPLE_LIST="$1"      # File containing sample names (one per line)
BAM_DIR="$2"          # Directory containing BAM files
OUTPUT_DIR="$3"       # Output directory for bigwig and Lorenz curve files
BINSIZE="${4:-1000}"  # Bin size for bigwig (default: 1000)

# Get sample name from array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# Input and output files
BAM_FILE="${BAM_DIR}/${SAMPLE}.bam"
BIGWIG_FILE="${OUTPUT_DIR}/${SAMPLE}.bw"
LORENZ_FILE="${OUTPUT_DIR}/${SAMPLE}_lorenz.csv"
GINI_FILE="${OUTPUT_DIR}/${SAMPLE}_gini.txt"

echo "Processing sample: ${SAMPLE}"
echo "BAM file: ${BAM_FILE}"
echo "Output bigwig: ${BIGWIG_FILE}"
echo "Output Lorenz curve: ${LORENZ_FILE}"

# Check if BAM file exists
if [ ! -f "${BAM_FILE}" ]; then
    echo "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

# Generate bigwig file
echo "Generating bigwig with ${BINSIZE}bp bins..."
module load deepTools
bamCoverage -b "${BAM_FILE}" -o "${BIGWIG_FILE}" --binSize ${BINSIZE} -of bigwig -p 4
module unload deepTools

echo "Bigwig generation complete"

# Compute Lorenz curve
echo "Computing Lorenz curve..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

python scripts/lorenz.py "${BIGWIG_FILE}" \
  -o "${LORENZ_FILE}" \
  -b ${BINSIZE} \
  --gini-output "${GINI_FILE}"

echo "Complete for ${SAMPLE}!"
