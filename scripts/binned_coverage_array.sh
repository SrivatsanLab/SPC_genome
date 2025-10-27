#!/bin/bash
#SBATCH -J binned_cov
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 1:00:00

###########################################################################################################################
# Generate binned coverage CSV files from bigwig files
###########################################################################################################################

set -euo pipefail

# Arguments
SAMPLE_LIST="$1"      # File containing sample names (one per line)
BIGWIG_DIR="$2"       # Directory containing bigwig files
OUTPUT_DIR="$3"       # Output directory for CSV files

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/envs/default_jupyter

# Run binned coverage Python script
# Note: The script expects 0-indexed task IDs internally
python scripts/binned_coverage.py "${SAMPLE_LIST}" "${BIGWIG_DIR}" "${OUTPUT_DIR}"

echo "Binned coverage complete for task ${SLURM_ARRAY_TASK_ID}"
