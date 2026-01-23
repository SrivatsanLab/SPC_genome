#!/bin/bash

##########################################################################################################################
# Generate sample list for bulk alignment array job
# Extracts unique clone names from fastq files (merging all passages per clone)
##########################################################################################################################

set -euo pipefail

FASTQ_DIR="$1"
OUTPUT_FILE="$2"
PATTERN="${3:-*_Clone_*}"  # Default pattern to match AAVS_Clone and PolE_Clone

echo "Generating sample list from: ${FASTQ_DIR}"
echo "Pattern: ${PATTERN}"
echo "Output: ${OUTPUT_FILE}"
echo ""

# Find all fastq files matching the pattern and extract unique sample names
# Sample names are in the format: AAVS_Clone_4_P1_S2_R1_001.fastq.gz
# We want to extract: AAVS_Clone_4, AAVS_Clone_5, etc. (without passage info)

# Use ls with pattern matching (more reliable with symlinks/permissions than find)
ls "${FASTQ_DIR}"/${PATTERN}_P*_R1_001.fastq.gz 2>/dev/null | \
    sed 's|.*/||' | \
    sed 's/_P[0-9]*_S[0-9]*_R[0-9]*_001.fastq.gz//' | \
    sort -u > "${OUTPUT_FILE}"

# Display the samples found
n_samples=$(wc -l < "${OUTPUT_FILE}")

echo "Found ${n_samples} unique samples:"
cat "${OUTPUT_FILE}"
echo ""
echo "Sample list saved to: ${OUTPUT_FILE}"
