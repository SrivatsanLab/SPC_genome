#!/bin/bash
#SBATCH -J jc_donor
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00

##########################################################################################################################
# Per-donor joint calling for a benchmarking donor group.
#
# Adapted from scripts/CapWGS/sc_var_array_gatk.sh (single-cell HaplotypeCaller pattern)
# and the existing GATK joint-calling pipeline (gatk_genomicsdb_import_array.sh +
# gatk_joint_calling_parallel_array.sh + gatk_joint_calling_parallel_merge.sh), but
# collapsed into a single-job, whole-genome pipeline appropriate for the small
# per-donor cohorts used here (10-50 cells).
#
# For one donor it:
#   1. Builds a sample map from the donor's per-cell GVCFs
#   2. Runs GATK GenomicsDBImport (whole genome, no scatter)
#   3. Runs GATK GenotypeGVCFs to produce a per-donor multi-sample VCF
#   4. Computes bcftools stats on the joint VCF
#
# Bulks are intentionally NOT joint-called with cells: bulks are used as the truth set
# for benchmarking, so they get called separately (single-sample HC) to keep the
# comparison non-circular.
#
# Usage as an array job, one task per donor:
#   sbatch --array=1-N joint_call_donor.sh <donor_list.txt> <gvcf_root> <output_root> <reference_dir>
#
# donor_list.txt is one donor_id per line. <gvcf_root>/<donor_id>/*.g.vcf.gz are inputs.
##########################################################################################################################

set -euo pipefail

DONOR_LIST="$1"          # File: one donor_id per line
GVCF_ROOT="$2"           # Root dir holding <donor>/<sample>.g.vcf.gz
OUTPUT_ROOT="$3"         # Output root for joint VCFs (per-donor subdirs created)
REFERENCE_DIR="$4"       # GATK reference dir (must contain Homo_sapiens_assembly38.fasta + dict)
INTERVALS="${5:-${REFERENCE_DIR}/wgs_calling_regions.hg38.interval_list}"

DONOR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DONOR_LIST")
if [ -z "$DONOR_ID" ]; then
    echo "ERROR: no donor for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
GVCF_DIR="${GVCF_ROOT}/${DONOR_ID}"
DONOR_OUT="${OUTPUT_ROOT}/${DONOR_ID}"
SAMPLE_MAP="${DONOR_OUT}/sample_map.txt"
GENOMICS_DB="${DONOR_OUT}/genomicsdb"
TMP_DIR="${DONOR_OUT}/tmp"
JOINT_VCF="${DONOR_OUT}/${DONOR_ID}_joint.vcf.gz"

mkdir -p "${DONOR_OUT}" "${TMP_DIR}"

echo "=========================================="
echo "Joint calling donor: ${DONOR_ID}"
echo "GVCF dir:            ${GVCF_DIR}"
echo "Reference:           ${REFERENCE}"
echo "Calling intervals:   ${INTERVALS}"
echo "Joint VCF:           ${JOINT_VCF}"
echo "=========================================="
echo ""

if [ ! -d "$GVCF_DIR" ]; then
    echo "ERROR: GVCF directory not found: $GVCF_DIR"
    exit 1
fi

##########################################################################################################################
# Step 1: Build sample map (one row per cell GVCF)
##########################################################################################################################

> "${SAMPLE_MAP}"
for gvcf in "${GVCF_DIR}"/*.g.vcf.gz; do
    [ -f "$gvcf" ] || continue
    sample_name=$(basename "$gvcf" .g.vcf.gz)
    printf "%s\t%s\n" "$sample_name" "$gvcf" >> "${SAMPLE_MAP}"
done

n_samples=$(wc -l < "${SAMPLE_MAP}")
if [ "$n_samples" -eq 0 ]; then
    echo "ERROR: no GVCFs found in ${GVCF_DIR}"
    exit 1
fi

echo "Found ${n_samples} GVCFs for donor ${DONOR_ID}"
echo ""

##########################################################################################################################
# Step 2: GenomicsDBImport (whole-genome, single shot)
##########################################################################################################################

# Remove an incomplete GenomicsDB if a previous run was interrupted
if [ -d "${GENOMICS_DB}" ] && [ ! -f "${GENOMICS_DB}/callset.json" ]; then
    echo "Removing incomplete GenomicsDB at ${GENOMICS_DB}..."
    rm -rf "${GENOMICS_DB}"
fi

module load GATK

if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
    echo "Step 2: GenomicsDBImport..."
    gatk --java-options "-Xmx28g -Xms28g" GenomicsDBImport \
        --sample-name-map "${SAMPLE_MAP}" \
        --genomicsdb-workspace-path "${GENOMICS_DB}" \
        -R "${REFERENCE}" \
        -L "${INTERVALS}" \
        --tmp-dir "${TMP_DIR}" \
        --reader-threads 4 \
        --batch-size 50
else
    echo "Step 2: GenomicsDB already exists; reusing ${GENOMICS_DB}"
fi

##########################################################################################################################
# Step 3: GenotypeGVCFs
##########################################################################################################################

echo "Step 3: GenotypeGVCFs..."
gatk --java-options "-Xmx28g" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V gendb://"${GENOMICS_DB}" \
    -L "${INTERVALS}" \
    -O "${JOINT_VCF}" \
    --tmp-dir "${TMP_DIR}"

module unload GATK

##########################################################################################################################
# Step 4: bcftools stats
##########################################################################################################################

echo "Step 4: bcftools stats..."
module load BCFtools/1.18-GCC-12.2.0
bcftools stats "${JOINT_VCF}" > "${JOINT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "Variant summary:"
bcftools view -H "${JOINT_VCF}" | wc -l | awk '{print "  Total variant sites: " $1}'
bcftools view -H -v snps   "${JOINT_VCF}" | wc -l | awk '{print "  SNP sites:           " $1}'
bcftools view -H -v indels "${JOINT_VCF}" | wc -l | awk '{print "  Indel sites:         " $1}'

module unload BCFtools

# Clean up tmp; keep GenomicsDB for reproducibility/restart
rm -rf "${TMP_DIR}"

echo ""
echo "=========================================="
echo "Done: ${DONOR_ID} -> ${JOINT_VCF}"
echo "=========================================="
