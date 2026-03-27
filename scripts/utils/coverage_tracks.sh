#!/bin/bash
#SBATCH -J coverage_track
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 2:00:00

###########################################################################################################################
# Generate coverage track plots from bigwig files
# This is an array job wrapper for the coverage_tracks.py script
#
# Usage:
#   sbatch --array=1-N scripts/coverage_tracks.sh <barcode_file> <bigwig_dir> <output_dir> [options]
#
# Arguments:
#   barcode_file: File containing sample names/barcodes (one per line)
#   bigwig_dir: Directory containing bigwig files
#   output_dir: Directory for output plots
#   options: Additional options to pass to coverage_tracks.py
#
# Examples of options:
#   --color '#ff1a5e'                    # Custom color
#   --dark-mode                           # Use dark mode
#   -c chr1                               # Plot only chr1
#   --linewidth 1.5                       # Thicker lines
#   --no-log-scale                        # Linear scale
#   --ylim 0 50                           # Custom y-axis limits
#   --gap-size 5000000                    # Smaller gaps between contigs
#   --no-contig-labels                    # Hide contig labels
#   --width 20 --height 3                 # Larger figure
###########################################################################################################################

set -euo pipefail

BARCODE_FILE="$1"
BIGWIG_DIR="$2"
OUTPUT_DIR="$3"
shift 3  # Remove first 3 arguments, rest are options

# Get barcode for this array task
BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

# Define input and output paths
BIGWIG_FILE="${BIGWIG_DIR}/${BARCODE}.bw"
OUTPUT_FILE="${OUTPUT_DIR}/${BARCODE}_coverage_track.png"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Activate Python environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

echo "Generating coverage track for: ${BARCODE}"
echo "Input: ${BIGWIG_FILE}"
echo "Output: ${OUTPUT_FILE}"

# Run coverage_tracks.py with provided options
python scripts/CapWGS_QC/coverage_tracks.py "${BIGWIG_FILE}" -o "${OUTPUT_FILE}" "$@"

echo "Complete!"
