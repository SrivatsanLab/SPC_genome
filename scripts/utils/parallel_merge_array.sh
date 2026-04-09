#!/bin/bash
#SBATCH --job-name=parallel_merge
#SBATCH --output=SLURM_outs/array_outs/parallel_merge_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH -p short

# Parallel SAM to BAM merge array job
# Each task merges a group of SAM files into one intermediate BAM
#
# Arguments:
#   $1: GROUP_DIR - Directory containing group_XX files (lists of SAM files)
#   $2: INTERMEDIATE_DIR - Output directory for intermediate BAMs
#   $3: SAMPLE_NAME - Sample name for output file naming

GROUP_DIR="$1"
INTERMEDIATE_DIR="$2"
SAMPLE_NAME="$3"

# Get group file for this array task
# Note: split creates files starting from 00, but SLURM arrays start from 1
GROUP_ID=$(printf "%02d" $((SLURM_ARRAY_TASK_ID - 1)))
GROUP_FILE="${GROUP_DIR}/group_${GROUP_ID}"

if [ ! -f "${GROUP_FILE}" ]; then
    echo "ERROR: Group file not found: ${GROUP_FILE}"
    exit 1
fi

# Output intermediate BAM for this group
INTERMEDIATE_BAM="${INTERMEDIATE_DIR}/${SAMPLE_NAME}_group_${GROUP_ID}.bam"

echo "=========================================="
echo "Parallel Merge - Group ${GROUP_ID}"
echo "=========================================="
echo "Group file: ${GROUP_FILE}"
echo "Output BAM: ${INTERMEDIATE_BAM}"
echo "=========================================="
echo ""

# Load SAMtools module
module load SAMtools

# Count SAM files in this group
SAM_COUNT=$(wc -l < "${GROUP_FILE}")
echo "Merging ${SAM_COUNT} SAM files..."

# Merge SAM files for this group
# -c: Combine @RG headers with colliding IDs (preserve read groups)
# -f 0x2: Filter for properly paired reads (removes duplicates from chunking)
# -O BAM: Output as BAM format
# -b: Input file list
samtools merge -@ 4 -c -b "${GROUP_FILE}" -O SAM - | \
    samtools view -@ 4 -f 0x2 -b - | \
    samtools sort -@ 4 -o "${INTERMEDIATE_BAM}"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Intermediate BAM created successfully"
    echo "  ${INTERMEDIATE_BAM}"

    # Report file size
    BAM_SIZE=$(du -h "${INTERMEDIATE_BAM}" | cut -f1)
    echo "  Size: ${BAM_SIZE}"
else
    echo ""
    echo "✗ ERROR: Merge failed for group ${GROUP_ID}"
    exit 1
fi

echo ""
echo "=========================================="
echo "Group ${GROUP_ID} merge complete!"
echo "=========================================="
