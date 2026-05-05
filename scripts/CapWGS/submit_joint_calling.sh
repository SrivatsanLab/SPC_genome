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
    # BCFtools mode: Interval-scattered merge + norm + filter, then concat
    ##################################################################################################################

    echo "BCFtools mode: Running parallelized joint calling with balanced intervals..."
    echo ""

    INTERVAL_LIST="${TMP_DIR}/interval_list.txt"
    PER_INTERVAL_DIR="${TMP_DIR}/per_interval_vcfs"
    FINAL_VCF="${RESULTS_DIR}/${OUTPUT_NAME}_joint.vcf.gz"
    SCATTER_COUNT=100

    mkdir -p "${PER_INTERVAL_DIR}"

    ##############################################################################################################
    # Step 1: Generate balanced genomic intervals
    ##############################################################################################################

    echo "Step 1: Generating ${SCATTER_COUNT} balanced genomic intervals..."

    interval_gen_job_ID=$(sbatch --parsable \
        "${SCRIPTS_DIR}/scripts/CapWGS/generate_intervals.sh" \
        "${REFERENCE_GENOME}" \
        "${TMP_DIR}" \
        "${SCATTER_COUNT}")

    echo "Interval generation job submitted: ${interval_gen_job_ID}"
    echo ""

    ##############################################################################################################
    # Step 2: After intervals are generated, dispatch the per-interval array + concat-merge
    # (Sized dynamically to the actual interval count, since SplitIntervals with
    #  BALANCING_WITHOUT_INTERVAL_SUBDIVISION may produce fewer bins than requested.)
    ##############################################################################################################

    echo "Step 2: Submitting dispatch job (sizes array + merge after intervals exist)..."

    dispatch_job_ID=$(sbatch --parsable \
        --dependency=afterok:${interval_gen_job_ID} \
        "${SCRIPTS_DIR}/scripts/CapWGS/bcftools_joint_calling_parallel_dispatch.sh" \
        "${SCRIPTS_DIR}" \
        "${ALIGNED_DIR}" \
        "${REFERENCE_GENOME}" \
        "${PER_INTERVAL_DIR}" \
        "${INTERVAL_LIST}" \
        "${FINAL_VCF}")

    echo "Dispatch job submitted: ${dispatch_job_ID}"
    echo ""

    echo "BCFtools joint calling pipeline submitted!"
    echo ""
    echo "Job chain:"
    echo "  1. Generate intervals: ${interval_gen_job_ID}"
    echo "  2. Dispatch + array + concat: ${dispatch_job_ID} (children submitted after intervals exist)"
    echo ""
    echo "Output files:"
    echo "  Intervals: ${TMP_DIR}/intervals/"
    echo "  Per-interval VCFs (temp): ${PER_INTERVAL_DIR}/"
    echo "  Final merged VCF: ${FINAL_VCF}"
    echo ""

elif [ "$VARIANT_CALLER" = "gatk" ]; then
    ##################################################################################################################
    # GATK mode: Interval-based GenomicsDBImport + GenotypeGVCFs + merge
    ##################################################################################################################

    echo "GATK mode: Running parallelized joint calling with interval-based GenomicsDB..."
    echo ""

    # Define paths
    GENOMICSDB_DIR="${RESULTS_DIR}/genomicsdb"  # Base directory for per-interval databases
    SAMPLE_MAP="${RESULTS_DIR}/sample_map.txt"
    INTERVAL_LIST="${TMP_DIR}/interval_list.txt"  # List of interval files
    PER_INTERVAL_DIR="${TMP_DIR}/per_interval_vcfs"  # Intermediate files go in temp directory
    FINAL_VCF="${RESULTS_DIR}/${OUTPUT_NAME}_joint.vcf.gz"
    SCATTER_COUNT=100  # Number of intervals to generate

    mkdir -p "${GENOMICSDB_DIR}"
    mkdir -p "${PER_INTERVAL_DIR}"

    ##############################################################################################################
    # Step 1: Generate balanced genomic intervals
    ##############################################################################################################

    echo "Step 1: Generating balanced genomic intervals..."

    interval_gen_job_ID=$(sbatch --parsable \
        "${SCRIPTS_DIR}/scripts/CapWGS/generate_intervals.sh" \
        "${REFERENCE_GENOME}" \
        "${TMP_DIR}" \
        "${SCATTER_COUNT}")

    echo "Interval generation job submitted: ${interval_gen_job_ID}"
    echo "  Creating ~${SCATTER_COUNT} balanced intervals"
    echo ""

    ##############################################################################################################
    # Step 2: Create sample map from GVCFs (runs while intervals are being generated)
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
    # Step 3: Submit per-interval GenomicsDB import array job (depends on interval generation)
    ##############################################################################################################

    echo "Step 3: Submitting per-interval GenomicsDB import array..."
    echo "  (Will run after interval generation completes)"

    genomicsdb_array_ID=$(sbatch --parsable \
        --dependency=afterok:${interval_gen_job_ID} \
        --array=1-${SCATTER_COUNT} \
        "${SCRIPTS_DIR}/scripts/CapWGS/gatk_genomicsdb_import_array.sh" \
        "${ALIGNED_DIR}" \
        "${GENOMICSDB_DIR}" \
        "${SAMPLE_MAP}" \
        "${REFERENCE_GENOME}" \
        "${INTERVAL_LIST}")

    echo "GenomicsDB import array submitted: ${genomicsdb_array_ID}"
    echo "  Importing ${SCATTER_COUNT} intervals in parallel"
    echo ""

    ##############################################################################################################
    # Step 4: Submit per-interval joint calling array job (depends on GenomicsDB)
    ##############################################################################################################

    echo "Step 4: Submitting per-interval joint calling array..."

    jc_array_ID=$(sbatch --parsable \
        --dependency=afterok:${genomicsdb_array_ID} \
        --array=1-${SCATTER_COUNT} \
        "${SCRIPTS_DIR}/scripts/CapWGS/gatk_joint_calling_parallel_array.sh" \
        "${GENOMICSDB_DIR}" \
        "${PER_INTERVAL_DIR}" \
        "${REFERENCE_GENOME}" \
        "${INTERVAL_LIST}")

    echo "Joint calling array submitted: ${jc_array_ID}"
    echo "  Calling ${SCATTER_COUNT} intervals in parallel"
    echo "  (Will run after GenomicsDB import completes)"
    echo ""

    ##############################################################################################################
    # Step 5: Submit merge job (depends on joint calling completion)
    ##############################################################################################################

    echo "Step 5: Submitting merge job..."

    merge_job_ID=$(sbatch --parsable \
        --dependency=afterok:${jc_array_ID} \
        "${SCRIPTS_DIR}/scripts/CapWGS/gatk_joint_calling_parallel_merge.sh" \
        "${PER_INTERVAL_DIR}" \
        "${FINAL_VCF}" \
        "${INTERVAL_LIST}")

    echo "Merge job submitted: ${merge_job_ID}"
    echo "  (Will run after joint calling completes)"
    echo ""

    ##############################################################################################################
    # Summary
    ##############################################################################################################

    echo "GATK joint calling pipeline submitted!"
    echo ""
    echo "Job chain:"
    echo "  0. Generate intervals: ${interval_gen_job_ID} (~${SCATTER_COUNT} balanced intervals)"
    echo "  1. GenomicsDB import: ${genomicsdb_array_ID} (${SCATTER_COUNT} intervals)"
    echo "  2. Joint calling: ${jc_array_ID} (${SCATTER_COUNT} intervals)"
    echo "  3. Merge VCFs: ${merge_job_ID}"
    echo ""
    echo "Output files:"
    echo "  Intervals: ${TMP_DIR}/intervals/"
    echo "  GenomicsDB (per-interval): ${GENOMICSDB_DIR}/"
    echo "  Per-interval VCFs (temp): ${PER_INTERVAL_DIR}/"
    echo "  Final merged VCF: ${FINAL_VCF}"
    echo ""
    echo "Note: Intervals and per-interval VCFs are intermediate files in temp directory"
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
