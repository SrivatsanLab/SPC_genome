#!/bin/bash
#
# Integration test for preprocessing pipeline
# This script runs the full CapWGS_PP.sh pipeline on test data
#

set -euo pipefail

# Default values (can be overridden by command-line arguments)
OUTPUT_NAME="sc_PolE_test"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
READ1="$PROJECT_ROOT/data/K562_tree/raw/k562_tree_1_S1_R1_001.fastq.gz"
READ2="$PROJECT_ROOT/data/K562_tree/raw/k562_tree_1_S1_R2_001.fastq.gz"
READ_COUNT=9000000000  # Estimate - adjust if you have exact count
REFERENCE_GENOME="/shared/biodata/reference/GATK/hg38/"
SCRIPTS_DIR="$PROJECT_ROOT"
OUTPUT_DIR="$PROJECT_ROOT/data/K562_tree"
N_CHUNKS=500
TMP_DIR="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"  # Default temp directory

# Function to display help
show_help() {
    cat << EOF
Convenience script to run the preprocessing pipeline

Usage: $0 [OPTIONS]

OPTIONS:
  -o    Output name (default: sc_PolE_test)
  -r    Read count (default: 9000000000)
  -n    Number of chunks (default: 500)
  -O    Output directory (default: .)
  -t    Temp directory (default: /hpc/temp/srivatsan_s/SPC_genome_preprocessing)
  -h    Show this help message

EXAMPLES:
  # Run with defaults
  $0

  # Run with custom output name and chunk count
  $0 -o my_test -n 1000

  # Run with specific temp directory (if not using default)
  $0 -t /path/to/alternative/temp/dir

  # Create a small subset test (10M reads)
  $0 -r 10000000 -n 50 -o subset_test

NOTES:
  - This script submits the pipeline as a SLURM job
  - Monitor progress with: squeue -u $USER
  - Check logs in: SLURM_outs/
  - See PREPROCESSING_TEST.md for detailed documentation

EOF
}

# Parse command-line arguments
while getopts "o:r:n:O:t:h" opt; do
  case $opt in
    o) OUTPUT_NAME="$OPTARG" ;;
    r) READ_COUNT="$OPTARG" ;;
    n) N_CHUNKS="$OPTARG" ;;
    O) OUTPUT_DIR="$OPTARG" ;;
    t) TMP_DIR="$OPTARG" ;;
    h) show_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
  esac
done

# Verify input files exist
if [ ! -L "$READ1" ] && [ ! -f "$READ1" ]; then
    echo "ERROR: Read 1 file not found: $READ1"
    echo "Have you run the setup to create data symlinks?"
    exit 1
fi

if [ ! -L "$READ2" ] && [ ! -f "$READ2" ]; then
    echo "ERROR: Read 2 file not found: $READ2"
    echo "Have you run the setup to create data symlinks?"
    exit 1
fi

# Verify reference genome exists
if [ ! -d "$REFERENCE_GENOME" ]; then
    echo "ERROR: Reference genome directory not found: $REFERENCE_GENOME"
    exit 1
fi

# Verify scripts directory exists
if [ ! -d "$SCRIPTS_DIR" ]; then
    echo "ERROR: Scripts directory not found: $SCRIPTS_DIR"
    exit 1
fi

# Create output directory if needed
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/SLURM_outs"

# Display configuration
echo "========================================="
echo "Preprocessing Pipeline Configuration"
echo "========================================="
echo "Output name:      $OUTPUT_NAME"
echo "Read 1:           $READ1"
echo "Read 2:           $READ2"
echo "Read count:       $READ_COUNT"
echo "Reference:        $REFERENCE_GENOME"
echo "Scripts dir:      $SCRIPTS_DIR"
echo "Output dir:       $OUTPUT_DIR"
echo "Number of chunks: $N_CHUNKS"
echo "Temp directory:   $TMP_DIR"

# Ensure temp directory exists
if [ ! -d "$TMP_DIR" ]; then
    echo "Creating temp directory: $TMP_DIR"
    mkdir -p "$TMP_DIR"
fi
echo "========================================="
echo ""

# Confirm before submitting
read -p "Submit preprocessing job? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

# Build command
CMD="sbatch $PROJECT_ROOT/CapWGS_PP.sh -o $OUTPUT_NAME -1 $READ1 -2 $READ2 -r $READ_COUNT -g $REFERENCE_GENOME -s $SCRIPTS_DIR -O $OUTPUT_DIR -n $N_CHUNKS"

if [ -n "$TMP_DIR" ]; then
    CMD="$CMD -t $TMP_DIR"
fi

# Submit job
echo "Submitting job..."
echo "Command: $CMD"
echo ""

JOB_ID=$($CMD)

echo "========================================="
echo "Job submitted successfully!"
echo "Job ID: $JOB_ID"
echo "========================================="
echo ""
echo "Monitor progress with:"
echo "  squeue -u $USER"
echo "  tail -f SLURM_outs/Preprocessing_*.out"
echo ""
echo "See PREPROCESSING_TEST.md for more details on monitoring and troubleshooting."
