#!/bin/bash
#SBATCH -J generate_intervals
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 30

##########################################################################################################################
# Generate balanced genomic intervals for parallelized joint calling
# Uses GATK SplitIntervals to create ~100 balanced intervals for faster parallelization
# Intervals are stored in temp directory and cleaned up after pipeline completion
##########################################################################################################################

set -euo pipefail

REFERENCE_DIR="$1"      # Reference directory (e.g., /shared/biodata/reference/GATK/hg38)
OUTPUT_DIR="$2"         # Output directory for interval files (temp directory)
SCATTER_COUNT="${3:-100}"  # Number of intervals to generate (default: 100)

echo "=========================================="
echo "Generating Genomic Intervals"
echo "=========================================="
echo "Reference: ${REFERENCE_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Scatter count: ${SCATTER_COUNT}"
echo ""

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
INTERVALS="${REFERENCE_DIR}/wgs_calling_regions.hg38.interval_list"
INTERVAL_DIR="${OUTPUT_DIR}/intervals"
INTERVAL_LIST="${OUTPUT_DIR}/interval_list.txt"

mkdir -p "${INTERVAL_DIR}"

echo "Starting interval generation..."
echo "Started at: $(date)"
echo ""

module load GATK

# Split intervals into balanced chunks
# This creates files like: 0000-scattered.interval_list, 0001-scattered.interval_list, etc.
gatk SplitIntervals \
    -R "${REFERENCE}" \
    -L "${INTERVALS}" \
    --scatter-count "${SCATTER_COUNT}" \
    -O "${INTERVAL_DIR}" \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION

module unload GATK

echo ""
echo "Interval generation complete at: $(date)"
echo ""

# Create interval list file (one interval file per line)
# Sort by filename to ensure consistent ordering
ls -1 "${INTERVAL_DIR}"/*.interval_list | sort > "${INTERVAL_LIST}"

N_INTERVALS=$(wc -l < "${INTERVAL_LIST}")

echo "=========================================="
echo "Interval generation successful!"
echo "=========================================="
echo ""
echo "Generated ${N_INTERVALS} intervals"
echo "Interval files: ${INTERVAL_DIR}/"
echo "Interval list: ${INTERVAL_LIST}"
echo ""
echo "Each interval is roughly balanced for computational load"
echo "These intervals will be used for parallel GenomicsDB import and joint calling"
echo ""
