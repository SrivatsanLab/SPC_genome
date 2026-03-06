#!/bin/bash
#SBATCH -J genomicsdb_import
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00

##########################################################################################################################
# Create GenomicsDB from single-cell GVCFs for CapWGS pipeline
# This needs to complete before per-chromosome joint calling can begin
##########################################################################################################################

set -euo pipefail

# Arguments
GVCF_DIR="$1"           # Directory containing single-cell GVCFs (e.g., data/sample/sc_outputs)
REFERENCE_DIR="$2"      # Reference directory (e.g., /shared/biodata/reference/GATK/hg38)
GENOMICS_DB="$3"        # Output GenomicsDB path (e.g., results/sample/genomicsdb_sample)
SAMPLE_MAP="$4"         # Output sample map path (e.g., results/sample/sample_map.txt)

echo "============================================"
echo "GenomicsDB Import for CapWGS"
echo "============================================"
echo ""
echo "Configuration:"
echo "  GVCF directory: ${GVCF_DIR}"
echo "  Sample map: ${SAMPLE_MAP}"
echo "  GenomicsDB: ${GENOMICS_DB}"
echo "  Reference dir: ${REFERENCE_DIR}"
echo ""

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
INTERVALS="${REFERENCE_DIR}/wgs_calling_regions.hg38.interval_list"
TMP_DIR="$(dirname ${GENOMICS_DB})/tmp"

mkdir -p "${TMP_DIR}"

# Create sample map from existing GVCFs
echo "Creating sample map from GVCFs..."
> "${SAMPLE_MAP}"

for gvcf in "${GVCF_DIR}"/*.g.vcf.gz; do
    if [ -f "$gvcf" ]; then
        # Extract sample name from filename (remove .g.vcf.gz extension)
        sample_name=$(basename "$gvcf" .g.vcf.gz)
        echo -e "${sample_name}\t${gvcf}" >> "${SAMPLE_MAP}"
    fi
done

# Check if any samples were found
n_samples=$(wc -l < "${SAMPLE_MAP}")
if [ "$n_samples" -eq 0 ]; then
    echo "ERROR: No GVCF files found in ${GVCF_DIR}"
    echo "Expected files matching: ${GVCF_DIR}/*.g.vcf.gz"
    exit 1
fi

echo "Found ${n_samples} samples for GenomicsDB import"
echo ""

# Remove incomplete GenomicsDB if it exists
if [ -d "${GENOMICS_DB}" ]; then
    if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
        echo "Removing incomplete GenomicsDB..."
        chmod -R +w "${GENOMICS_DB}" 2>/dev/null || true
        find "${GENOMICS_DB}" -type f -delete 2>/dev/null || true
        find "${GENOMICS_DB}" -depth -type d -exec rmdir {} \; 2>/dev/null || true
        rm -rf "${GENOMICS_DB}" 2>/dev/null || true
        # If directory still exists, just rename it
        if [ -d "${GENOMICS_DB}" ]; then
            echo "Could not fully remove GenomicsDB, renaming it..."
            mv "${GENOMICS_DB}" "${GENOMICS_DB}_old_$(date +%s)" || true
        fi
    else
        echo "ERROR: GenomicsDB already exists and appears complete."
        echo "Remove it manually if you want to recreate it:"
        echo "  rm -rf ${GENOMICS_DB}"
        exit 1
    fi
fi

echo "Starting GenomicsDB import..."
echo "Started at: $(date)"
echo ""

module load GATK

# Import all GVCFs into GenomicsDB
gatk --java-options "-Xmx56g -Xms56g" GenomicsDBImport \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${GENOMICS_DB}" \
    -R "${REFERENCE}" \
    -L "${INTERVALS}" \
    --tmp-dir "${TMP_DIR}" \
    --reader-threads 8

module unload GATK

echo ""
echo "GenomicsDB import complete at: $(date)"
echo ""
echo "============================================"
echo "GenomicsDB creation successful!"
echo "============================================"
echo ""
echo "Output:"
echo "  GenomicsDB: ${GENOMICS_DB}"
echo "  Sample map: ${SAMPLE_MAP}"
echo "  Samples: ${n_samples}"
echo ""
