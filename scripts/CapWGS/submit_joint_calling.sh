#!/bin/bash
#SBATCH -J submit_jc
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 30

###########################################################################################################################
# Wrapper script that submits joint calling
# This script runs AFTER single-cell variant calling completes
###########################################################################################################################

set -euo pipefail

REFERENCE_GENOME="$1"
BIN_DIR="$2"
ALIGNED_DIR="$3"
RESULTS_DIR="$4"
SCRIPTS_DIR="$5"
OUTPUT_NAME="$6"
VARIANT_CALLER="${7:-gatk}"  # Default to gatk for backward compatibility
TMP_DIR="${8:-/hpc/temp/srivatsan_s/joint_calling_${OUTPUT_NAME}}"  # Temp directory for intermediate files

echo "=========================================="
echo "Joint Calling Submission"
echo "=========================================="
echo "Sample: ${OUTPUT_NAME}"
echo "Variant caller: ${VARIANT_CALLER}"
echo "=========================================="
echo ""

# Paths
REFERENCE="${REFERENCE_GENOME}/Homo_sapiens_assembly38.fasta"

mkdir -p "${RESULTS_DIR}"

##########################################################################################################################
# Handle variant caller-specific joint calling
##########################################################################################################################

if [ "$VARIANT_CALLER" = "bcftools" ]; then
    ##################################################################################################################
    # BCFtools mode: Merge and normalize VCFs
    ##################################################################################################################

    echo "BCFtools mode: Merging and normalizing VCFs..."
    echo ""

    OUTPUT_VCF="${RESULTS_DIR}/${OUTPUT_NAME}_joint.vcf.gz"

    # Submit BCFtools merge and normalize job
    bcf_jc_ID=$(sbatch --parsable \
        "${SCRIPTS_DIR}/scripts/CapWGS/bcftools_joint_calling.sh" \
        "${ALIGNED_DIR}" \
        "${REFERENCE_GENOME}" \
        "${OUTPUT_VCF}" \
        "${OUTPUT_NAME}")

    echo "BCFtools joint calling job submitted: ${bcf_jc_ID}"
    echo "Output VCF will be: ${OUTPUT_VCF}"
    echo ""

elif [ "$VARIANT_CALLER" = "gatk" ]; then
    ##################################################################################################################
    # GATK mode: Per-chromosome GenomicsDBImport + GenotypeGVCFs + merge
    ##################################################################################################################

    echo "GATK mode: Running parallelized joint calling with per-chromosome GenomicsDB..."
    echo ""

    # Define paths
    INTERVALS="${REFERENCE_GENOME}/wgs_calling_regions.hg38.interval_list"
    GENOMICSDB_DIR="${RESULTS_DIR}/genomicsdb"  # Base directory for per-chromosome databases
    SAMPLE_MAP="${RESULTS_DIR}/sample_map.txt"
    CHROMOSOME_FILE="${RESULTS_DIR}/chromosome_list.txt"
    PER_CHR_DIR="${TMP_DIR}/per_chromosome_vcfs"  # Intermediate files go in temp directory
    FINAL_VCF="${RESULTS_DIR}/${OUTPUT_NAME}_joint.vcf.gz"

    mkdir -p "${GENOMICSDB_DIR}"
    mkdir -p "${PER_CHR_DIR}"

    ##############################################################################################################
    # Step 1: Extract chromosome list from GATK interval list
    ##############################################################################################################

    echo "Step 1: Extracting chromosome list from GATK interval list..."

    # Extract unique chromosome names from wgs_calling_regions.hg38.interval_list
    grep -v "^@" "${INTERVALS}" | cut -f1 | sort -u > "${CHROMOSOME_FILE}"

    N_CHROMOSOMES=$(wc -l < "${CHROMOSOME_FILE}")
    echo "Extracted ${N_CHROMOSOMES} chromosomes from ${INTERVALS}"
    echo ""

    ##############################################################################################################
    # Step 2: Create sample map from GVCFs
    ##############################################################################################################

    echo "Step 2: Creating sample map from GVCFs..."

    > "${SAMPLE_MAP}"
    for gvcf in "${ALIGNED_DIR}"/*.g.vcf.gz; do
        if [ -f "$gvcf" ]; then
            sample_name=$(basename "$gvcf" .g.vcf.gz)
            echo -e "${sample_name}\t${gvcf}" >> "${SAMPLE_MAP}"
        fi
    done

    n_cells=$(wc -l < "${SAMPLE_MAP}")
    if [ "$n_cells" -eq 0 ]; then
        echo "ERROR: No GVCF files found in ${ALIGNED_DIR}"
        exit 1
    fi

    echo "Found ${n_cells} cells"
    echo ""

    ##############################################################################################################
    # Step 3: Submit per-chromosome GenomicsDB import array job
    ##############################################################################################################

    echo "Step 3: Submitting per-chromosome GenomicsDB import array..."

    genomicsdb_array_ID=$(sbatch --parsable \
        --array=1-${N_CHROMOSOMES} \
        "${SCRIPTS_DIR}/scripts/CapWGS/gatk_genomicsdb_import_array.sh" \
        "${ALIGNED_DIR}" \
        "${GENOMICSDB_DIR}" \
        "${SAMPLE_MAP}" \
        "${REFERENCE_GENOME}" \
        "${CHROMOSOME_FILE}")

    echo "GenomicsDB import array submitted: ${genomicsdb_array_ID}"
    echo "  Importing ${N_CHROMOSOMES} chromosomes in parallel"
    echo ""

    ##############################################################################################################
    # Step 4: Submit per-chromosome joint calling array job (depends on GenomicsDB)
    ##############################################################################################################

    echo "Step 4: Submitting per-chromosome joint calling array..."

    jc_array_ID=$(sbatch --parsable \
        --dependency=afterok:${genomicsdb_array_ID} \
        --array=1-${N_CHROMOSOMES} \
        "${SCRIPTS_DIR}/scripts/CapWGS/gatk_joint_calling_parallel_array.sh" \
        "${GENOMICSDB_DIR}" \
        "${PER_CHR_DIR}" \
        "${REFERENCE_GENOME}" \
        "${CHROMOSOME_FILE}")

    echo "Joint calling array submitted: ${jc_array_ID}"
    echo "  Calling ${N_CHROMOSOMES} chromosomes in parallel"
    echo "  (Will run after GenomicsDB import completes)"
    echo ""

    ##############################################################################################################
    # Step 5: Submit merge job (depends on joint calling completion)
    ##############################################################################################################

    echo "Step 5: Submitting merge job..."

    merge_job_ID=$(sbatch --parsable \
        --dependency=afterok:${jc_array_ID} \
        "${SCRIPTS_DIR}/scripts/CapWGS/gatk_joint_calling_parallel_merge.sh" \
        "${PER_CHR_DIR}" \
        "${FINAL_VCF}" \
        "${CHROMOSOME_FILE}")

    echo "Merge job submitted: ${merge_job_ID}"
    echo "  (Will run after joint calling completes)"
    echo ""

    ##############################################################################################################
    # Summary
    ##############################################################################################################

    echo "GATK joint calling pipeline submitted!"
    echo ""
    echo "Job chain:"
    echo "  1. GenomicsDB import: ${genomicsdb_array_ID} (${N_CHROMOSOMES} chromosomes)"
    echo "  2. Joint calling: ${jc_array_ID} (${N_CHROMOSOMES} chromosomes)"
    echo "  3. Merge VCFs: ${merge_job_ID}"
    echo ""
    echo "Output files:"
    echo "  GenomicsDB (per-chromosome): ${GENOMICSDB_DIR}/"
    echo "  Per-chromosome VCFs (temp): ${PER_CHR_DIR}/"
    echo "  Final merged VCF: ${FINAL_VCF}"
    echo ""
    echo "Note: Per-chromosome VCFs are intermediate files in temp directory"
    echo "      They will be automatically deleted when temp directory is cleaned"
    echo ""

else
    echo "ERROR: Invalid variant caller: ${VARIANT_CALLER}"
    echo "Must be 'gatk' or 'bcftools'"
    exit 1
fi

echo "=========================================="
echo "Joint calling submission complete!"
echo "=========================================="
echo ""
