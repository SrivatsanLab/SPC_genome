#!/bin/bash
#SBATCH -J submit_final
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 1:00:00

##########################################################################################################################
# Submit final processing steps for CapGTA pipeline:
# 1. Single-cell variant calling (bcftools) on DNA BAMs
# 2. Merge single-cell VCFs
# 3. Generate RNA count matrix from RNA BAMs
##########################################################################################################################

set -euo pipefail

SC_OUTPUTS_DIR="$1"
BARCODE_FILE="$2"
REFERENCE="$3"
GTF_FILE="$4"
RESULTS_DIR="$5"
SCRIPTS_DIR="$6"
SAMPLE_NAME="$7"

# Function to wait for a SLURM job to complete
wait_for_job() {
    local job_id=$1
    while squeue -j "$job_id" -h -t PENDING,RUNNING 2>/dev/null | grep -q .; do
        sleep 60
    done
    # Check if any tasks failed
    if sacct -j "$job_id" --format=State -n | grep -q FAILED; then
        echo "Job $job_id had failures" >&2
        return 1
    fi
}

# Determine reference FASTA path (handle directory or file input)
if [ -f "${REFERENCE}" ]; then
    REFERENCE_FA="${REFERENCE}"
elif [ -d "${REFERENCE}" ]; then
    # Search for FASTA file in standard locations
    if [ -f "${REFERENCE}/BWAIndex/genome.fa" ]; then
        REFERENCE_FA="${REFERENCE}/BWAIndex/genome.fa"
    elif [ -f "${REFERENCE}/genome.fa" ]; then
        REFERENCE_FA="${REFERENCE}/genome.fa"
    elif compgen -G "${REFERENCE}/BWAIndex/"*.fa > /dev/null 2>&1; then
        REFERENCE_FA=$(ls "${REFERENCE}/BWAIndex/"*.fa 2>/dev/null | head -1)
    elif compgen -G "${REFERENCE}/BWAIndex/"*.fna > /dev/null 2>&1; then
        REFERENCE_FA=$(ls "${REFERENCE}/BWAIndex/"*.fna 2>/dev/null | head -1)
    elif compgen -G "${REFERENCE}/"*.fa > /dev/null 2>&1; then
        REFERENCE_FA=$(ls "${REFERENCE}/"*.fa 2>/dev/null | head -1)
    elif compgen -G "${REFERENCE}/"*.fna > /dev/null 2>&1; then
        REFERENCE_FA=$(ls "${REFERENCE}/"*.fna 2>/dev/null | head -1)
    else
        echo "Error: Could not find FASTA file in ${REFERENCE}"
        exit 1
    fi
else
    echo "Error: Reference is neither file nor directory: ${REFERENCE}"
    exit 1
fi

# Count number of cells
NUM_CELLS=$(wc -l < "${BARCODE_FILE}")

echo "=========================================="
echo "CapGTA Final Processing"
echo "=========================================="
echo "Sample: ${SAMPLE_NAME}"
echo "Number of cells: ${NUM_CELLS}"
echo "SC outputs: ${SC_OUTPUTS_DIR}"
echo "Results: ${RESULTS_DIR}"
echo "Reference FASTA: ${REFERENCE_FA}"
echo "GTF: ${GTF_FILE}"
echo "=========================================="
echo ""

######################################################################################################
#### Step 1: Single-cell variant calling (bcftools)

echo "Submitting single-cell variant calling array (${NUM_CELLS} cells)..."

sc_vcall_job=$(sbatch --parsable \
    --array=1-${NUM_CELLS} \
    "${SCRIPTS_DIR}/scripts/CapGTA/sc_variant_calling_bcftools_array.sh" \
    "${BARCODE_FILE}" \
    "${SC_OUTPUTS_DIR}" \
    "${REFERENCE_FA}" \
    "${SC_OUTPUTS_DIR}")

echo "  Variant calling job ID: ${sc_vcall_job}"
echo "  Waiting for variant calling to complete..."
wait_for_job "${sc_vcall_job}"
echo "  ✓ Variant calling complete!"
echo ""

######################################################################################################
#### Step 2: Merge single-cell VCFs

echo "Submitting VCF merge job..."

vcf_merge_job=$(sbatch --parsable \
    "${SCRIPTS_DIR}/scripts/CapGTA/merge_sc_vcfs.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${RESULTS_DIR}/sc_variants_merged.vcf.gz")

echo "  VCF merge job ID: ${vcf_merge_job}"
echo "  (Not waiting - final output)"
echo ""

######################################################################################################
#### Step 3: Generate RNA count matrix

echo "Submitting RNA count matrix generation..."

rna_counts_job=$(sbatch --parsable \
    "${SCRIPTS_DIR}/scripts/CapGTA/create_rna_count_matrix.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${GTF_FILE}" \
    "${RESULTS_DIR}/rna_counts")

echo "  RNA count matrix job ID: ${rna_counts_job}"
echo "  (Not waiting - final output)"
echo ""

echo "=========================================="
echo "Final processing jobs submitted!"
echo "=========================================="
echo "  1. Variant calling: ${sc_vcall_job} (completed)"
echo "  2. VCF merge: ${vcf_merge_job} (running)"
echo "  3. RNA counts: ${rna_counts_job} (running)"
echo ""
echo "Final outputs:"
echo "  - Merged VCF: ${RESULTS_DIR}/sc_variants_merged.vcf.gz"
echo "  - RNA count matrix: ${RESULTS_DIR}/rna_counts.matrix"
echo ""
