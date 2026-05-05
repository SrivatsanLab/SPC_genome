#!/bin/bash
#SBATCH -J bcf_jc_interval
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00

##########################################################################################################################
# Per-interval bcftools joint calling for CapWGS pipeline.
# Mirrors the GATK interval-scatter pattern: for each balanced interval produced by generate_intervals.sh,
# run bcftools merge -> norm -> filter (ALT=".") -> index, restricted to that interval.
# Outputs ${INTERVAL_NAME}.vcf.gz which the existing gatk_joint_calling_parallel_merge.sh concats.
##########################################################################################################################

set -euo pipefail

# Arguments
VCF_DIR="$1"            # Directory containing per-cell *.g.vcf.gz files
REFERENCE_DIR="$2"      # Reference directory (e.g., /shared/biodata/reference/GATK/hg38)
PER_INTERVAL_DIR="$3"   # Output directory for per-interval VCFs
INTERVAL_LIST="$4"      # File with one interval_list path per line

# Get the interval file for this array task
INTERVAL_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INTERVAL_LIST")
INTERVAL_NAME=$(basename "$INTERVAL_FILE" .interval_list)

echo "============================================"
echo "BCFtools per-interval joint calling: ${INTERVAL_NAME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "============================================"
echo ""

# Resolve reference
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
if [ ! -f "${REFERENCE}" ]; then
    REFERENCE=$(find "${REFERENCE_DIR}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) | head -n 1)
    if [ -z "${REFERENCE}" ] || [ ! -f "${REFERENCE}" ]; then
        echo "ERROR: Reference genome not found in ${REFERENCE_DIR}"
        exit 1
    fi
fi

if [ ! -f "${INTERVAL_FILE}" ]; then
    echo "ERROR: Interval file not found: ${INTERVAL_FILE}"
    exit 1
fi

mkdir -p "${PER_INTERVAL_DIR}"

OUTPUT_VCF="${PER_INTERVAL_DIR}/${INTERVAL_NAME}.vcf.gz"
REGIONS_FILE="${PER_INTERVAL_DIR}/${INTERVAL_NAME}.regions.txt"
VCF_LIST="${PER_INTERVAL_DIR}/${INTERVAL_NAME}_vcf_list.txt"

echo "Reference: ${REFERENCE}"
echo "Interval file: ${INTERVAL_FILE}"
echo "VCF directory: ${VCF_DIR}"
echo "Output VCF: ${OUTPUT_VCF}"
echo ""

##########################################################################################################################
# Step 1: Convert GATK interval_list -> bcftools regions file (1-based inclusive, three columns)
##########################################################################################################################

echo "Step 1: Converting interval_list to bcftools regions file..."
awk '!/^@/ {print $1"\t"$2"\t"$3}' "${INTERVAL_FILE}" > "${REGIONS_FILE}"

n_regions=$(wc -l < "${REGIONS_FILE}")
if [ "${n_regions}" -eq 0 ]; then
    echo "ERROR: No regions extracted from ${INTERVAL_FILE}"
    exit 1
fi
echo "  ${n_regions} regions"
echo ""

##########################################################################################################################
# Step 2: Build VCF input list
##########################################################################################################################

echo "Step 2: Listing per-cell GVCFs..."
find "${VCF_DIR}" -maxdepth 1 -name "*.g.vcf.gz" | sort > "${VCF_LIST}"
n_cells=$(wc -l < "${VCF_LIST}")
if [ "${n_cells}" -eq 0 ]; then
    echo "ERROR: No GVCFs found in ${VCF_DIR}"
    exit 1
fi
echo "  ${n_cells} GVCFs"
echo ""

##########################################################################################################################
# Step 3: Stream merge -> norm -> filter (ALT=".") -> bgzip
##########################################################################################################################

echo "Step 3: Running bcftools merge | norm | view for ${INTERVAL_NAME}..."
echo "Started at: $(date)"

module load BCFtools

bcftools merge \
    --threads 3 \
    --regions-file "${REGIONS_FILE}" \
    -m none \
    -O u \
    --file-list "${VCF_LIST}" | \
bcftools norm \
    --threads 3 \
    -c x \
    -m -both \
    -f "${REFERENCE}" \
    -O u | \
bcftools view \
    --threads 3 \
    -e 'ALT="."' \
    -O z \
    -o "${OUTPUT_VCF}"

echo "Pipeline complete at: $(date)"
echo ""

##########################################################################################################################
# Step 4: Index
##########################################################################################################################

echo "Step 4: Indexing ${OUTPUT_VCF}..."
bcftools index --threads 3 -f "${OUTPUT_VCF}"

module unload BCFtools

echo ""
echo "============================================"
echo "Done: ${INTERVAL_NAME}"
echo "============================================"
echo "Output: ${OUTPUT_VCF}"
echo ""
