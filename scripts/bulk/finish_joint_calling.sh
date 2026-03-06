#!/bin/bash
#SBATCH -J bulk_joint_call
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 7-00:00:00

##########################################################################################################################
# Finish joint genotyping of bulk K562 samples using existing GenomicsDB
# This script runs GenotypeGVCFs on the already-created GenomicsDB workspace
# Increased time limit to 7 days to ensure completion
##########################################################################################################################

set -euo pipefail

# Configuration
COHORT_NAME="K562_mut_accumulation"
OUTPUT_DIR="data/${COHORT_NAME}"
REFERENCE_DIR="/shared/biodata/reference/GATK/hg38"

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
INTERVALS="${REFERENCE_DIR}/wgs_calling_regions.hg38.interval_list"
GENOMICS_DB="${OUTPUT_DIR}/vcfs/genomicsdb_${COHORT_NAME}"
OUTPUT_VCF="${OUTPUT_DIR}/vcfs/${COHORT_NAME}_joint.vcf.gz"
TMP_DIR="${OUTPUT_DIR}/vcfs/tmp"

echo "============================================"
echo "Resuming joint calling for: ${COHORT_NAME}"
echo "============================================"
echo ""
echo "Configuration:"
echo "  GenomicsDB: ${GENOMICS_DB}"
echo "  Output VCF: ${OUTPUT_VCF}"
echo "  Reference: ${REFERENCE}"
echo "  Time limit: 7 days"
echo ""

mkdir -p "${TMP_DIR}"

# Check if GenomicsDB is valid
GVCF_DIR="${OUTPUT_DIR}/gvcfs"
SAMPLE_MAP="${OUTPUT_DIR}/vcfs/${COHORT_NAME}_sample_map.txt"

if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
    echo "WARNING: GenomicsDB is incomplete or corrupted (missing callset.json)"
    echo "Regenerating GenomicsDB from scratch..."
    echo ""

    # Remove incomplete GenomicsDB
    if [ -d "${GENOMICS_DB}" ]; then
        echo "Removing incomplete GenomicsDB..."
        rm -rf "${GENOMICS_DB}"
    fi

    ##########################################################################################################################
    # Step 1a: Import GVCFs into GenomicsDB
    ##########################################################################################################################

    echo "Step 1a: Importing GVCFs into GenomicsDB..."
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
        --reader-threads 4

    echo ""
    echo "GenomicsDB import complete at: $(date)"
    echo ""
fi

##########################################################################################################################
# Step 1b: Joint genotyping with GenotypeGVCFs
##########################################################################################################################

echo "Step 1b: Running GenotypeGVCFs..."
echo "Started at: $(date)"
echo ""

module load GATK

# Perform joint genotyping
gatk --java-options "-Xmx56g" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V gendb://"${GENOMICS_DB}" \
    -L "${INTERVALS}" \
    -O "${OUTPUT_VCF}" \
    --tmp-dir "${TMP_DIR}"

echo ""
echo "GenotypeGVCFs completed at: $(date)"
echo ""

##########################################################################################################################
# Step 2: Generate statistics
##########################################################################################################################

echo "Step 2: Generating variant statistics..."

module load BCFtools

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "============================================"
echo "Joint calling complete!"
echo "============================================"
echo ""
echo "Output files:"
echo "  Joint VCF: ${OUTPUT_VCF}"
echo "  Joint VCF index: ${OUTPUT_VCF}.tbi"
echo "  Statistics: ${OUTPUT_VCF%.vcf.gz}_stats.txt"
echo ""

# Variant summary
echo "=== Variant Summary ==="
n_samples=$(bcftools query -l "${OUTPUT_VCF}" | wc -l)
echo "Total samples: ${n_samples}"
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

module unload BCFtools GATK

echo ""
echo "Done!"
