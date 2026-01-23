#!/bin/bash
#SBATCH -J bulk_joint_call
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00

##########################################################################################################################
# Joint genotyping of bulk samples using GATK GenotypeGVCFs
# This script combines all GVCFs from individual samples into a single unified VCF
# Following GATK best practices for joint calling
##########################################################################################################################

set -euo pipefail

# Arguments
GVCF_DIR="$1"           # Directory containing individual GVCF files
OUTPUT_VCF="$2"         # Output joint-called VCF file
REFERENCE_DIR="$3"      # Reference directory (e.g., /shared/biodata/reference/GATK/hg38)
COHORT_NAME="$4"        # Name for this cohort (e.g., "AAVS_PolE_clones")

echo "============================================"
echo "Joint calling for cohort: ${COHORT_NAME}"
echo "============================================"
echo ""

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
INTERVALS="${REFERENCE_DIR}/wgs_calling_regions.hg38.interval_list"
OUTPUT_DIR=$(dirname "${OUTPUT_VCF}")
GENOMICS_DB="${OUTPUT_DIR}/genomicsdb_${COHORT_NAME}"
SAMPLE_MAP="${OUTPUT_DIR}/${COHORT_NAME}_sample_map.txt"
TMP_DIR="${OUTPUT_DIR}/tmp"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMP_DIR}"

##########################################################################################################################
# Step 1: Create sample map file for GenomicsDBImport
##########################################################################################################################

echo "Step 1: Creating sample map file..."

# Find all GVCF files and create sample map
# Format: sample_name\t/path/to/sample.g.vcf.gz
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
    exit 1
fi

echo "Found ${n_samples} samples for joint calling:"
cat "${SAMPLE_MAP}"
echo ""

##########################################################################################################################
# Step 2: Import GVCFs into GenomicsDB
##########################################################################################################################

echo "Step 2: Importing GVCFs into GenomicsDB..."

module load GATK

# Remove existing GenomicsDB if it exists
if [ -d "${GENOMICS_DB}" ]; then
    echo "Removing existing GenomicsDB: ${GENOMICS_DB}"
    rm -rf "${GENOMICS_DB}"
fi

# Import all GVCFs into GenomicsDB
# Using WGS calling regions intervals
gatk --java-options "-Xmx32g -Xms32g" GenomicsDBImport \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${GENOMICS_DB}" \
    -R "${REFERENCE}" \
    -L "${INTERVALS}" \
    --tmp-dir "${TMP_DIR}" \
    --reader-threads 4

echo "GenomicsDB import complete."
echo ""

##########################################################################################################################
# Step 3: Joint genotyping with GenotypeGVCFs
##########################################################################################################################

echo "Step 3: Joint genotyping with GenotypeGVCFs..."

# Perform joint genotyping
gatk --java-options "-Xmx32g" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V gendb://"${GENOMICS_DB}" \
    -L "${INTERVALS}" \
    -O "${OUTPUT_VCF}" \
    --tmp-dir "${TMP_DIR}"

echo "Joint genotyping complete."
echo ""

##########################################################################################################################
# Step 4: Variant Quality Score Recalibration (VQSR) - Optional
##########################################################################################################################

# Note: VQSR is typically recommended for large cohorts (>30 samples)
# For smaller cohorts, hard filtering is more appropriate
# Uncommenting the VQSR steps below if you have a large enough cohort

echo "Step 4: Skipping VQSR (recommended for cohorts >30 samples)"
echo "For smaller cohorts, use hard filtering based on GATK recommendations"
echo ""

# VQSR code (commented out):
# DBSNP="${REFERENCE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf"
# HAPMAP="${REFERENCE_DIR}/hapmap_3.3.hg38.vcf.gz"
# OMNI="${REFERENCE_DIR}/1000G_omni2.5.hg38.vcf.gz"
# THOUSAND_G="${REFERENCE_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
# MILLS_INDELS="${REFERENCE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
#
# # SNP recalibration
# gatk --java-options "-Xmx32g" VariantRecalibrator \
#     -R "${REFERENCE}" \
#     -V "${OUTPUT_VCF}" \
#     --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${HAPMAP}" \
#     --resource:omni,known=false,training=true,truth=false,prior=12.0 "${OMNI}" \
#     --resource:1000G,known=false,training=true,truth=false,prior=10.0 "${THOUSAND_G}" \
#     --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${DBSNP}" \
#     -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
#     -mode SNP \
#     -O "${OUTPUT_DIR}/${COHORT_NAME}_snp.recal" \
#     --tranches-file "${OUTPUT_DIR}/${COHORT_NAME}_snp.tranches"

##########################################################################################################################
# Step 5: Generate statistics
##########################################################################################################################

echo "Step 5: Generating variant statistics..."

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
echo "  Sample map: ${SAMPLE_MAP}"
echo "  GenomicsDB: ${GENOMICS_DB}"
echo "  Statistics: ${OUTPUT_VCF%.vcf.gz}_stats.txt"
echo ""

# Variant summary
echo "=== Variant Summary ==="
echo "Total samples: ${n_samples}"
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

module unload BCFtools GATK

echo ""
echo "Done!"
