#!/bin/bash
#SBATCH -J sc_bam_gta
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1

##########################################################################################################################
# Extract single-cell DNA and RNA BAMs from concatenated bulk BAMs, then run variant calling and RNA counting
#
# This script:
#   1. Submits array jobs to extract DNA and RNA BAMs for each cell from bulk BAMs
#   2. Submits variant calling jobs on DNA BAMs
#   3. Submits RNA counting jobs on RNA BAMs
##########################################################################################################################

set -euo pipefail

BARCODE_FILE="$1"       # File with list of real cell barcodes
DNA_BAM="$2"            # Bulk DNA BAM file
RNA_BAM="$3"            # Bulk RNA BAM file
SC_OUTPUTS_DIR="$4"     # Output directory for single-cell BAMs and VCFs
RESULTS_DIR="$5"        # Results directory for merged outputs
SCRIPTS_DIR="$6"        # Scripts directory
OUTPUT_NAME="$7"        # Sample name

echo "=========================================="
echo "CapGTA Single Cell Extraction and Analysis"
echo "=========================================="
echo "Barcode file: ${BARCODE_FILE}"
echo "DNA BAM: ${DNA_BAM}"
echo "RNA BAM: ${RNA_BAM}"
echo "SC outputs: ${SC_OUTPUTS_DIR}"
echo "Results: ${RESULTS_DIR}"
echo "Sample: ${OUTPUT_NAME}"
echo "=========================================="
echo ""

# Verify inputs
if [ ! -f "${BARCODE_FILE}" ]; then
    echo "ERROR: Barcode file not found: ${BARCODE_FILE}"
    exit 1
fi

if [ ! -f "${DNA_BAM}" ]; then
    echo "ERROR: DNA BAM file not found: ${DNA_BAM}"
    exit 1
fi

if [ ! -f "${RNA_BAM}" ]; then
    echo "ERROR: RNA BAM file not found: ${RNA_BAM}"
    exit 1
fi

# Create output directory
mkdir -p "${SC_OUTPUTS_DIR}"

# Count cells
cell_count=$(wc -l < "${BARCODE_FILE}")

echo "Found ${cell_count} cells to extract"
echo ""

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells to extract. Barcode file is empty."
    exit 1
fi

######################################################################################################
#### Extract DNA BAMs
echo "Submitting DNA extraction array job (${cell_count} cells)..."

dna_extraction_job_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    "${SCRIPTS_DIR}/scripts/utils/extract_sc_from_bam_array.sh" \
    "${DNA_BAM}" \
    "${BARCODE_FILE}" \
    "${SC_OUTPUTS_DIR}" \
    "_dna")

echo "DNA extraction job ID: ${dna_extraction_job_ID}"

######################################################################################################
#### Extract RNA BAMs
echo "Submitting RNA extraction array job (${cell_count} cells)..."

rna_extraction_job_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    "${SCRIPTS_DIR}/scripts/utils/extract_sc_from_bam_array.sh" \
    "${RNA_BAM}" \
    "${BARCODE_FILE}" \
    "${SC_OUTPUTS_DIR}" \
    "_rna")

echo "RNA extraction job ID: ${rna_extraction_job_ID}"
echo ""

######################################################################################################
#### Submit single-cell variant calling (BCFtools)
echo "Submitting variant calling array job..."

REFERENCE_FA="data/reference/worm_GCA_028201515.1_STAR/GCA_028201515.1_genomic_renamed.fna"

vcall_job_ID=$(sbatch --parsable \
    --dependency=afterok:${dna_extraction_job_ID} \
    --array=1-${cell_count} \
    "${SCRIPTS_DIR}/scripts/CapGTA/sc_variant_calling_bcftools_array.sh" \
    "${BARCODE_FILE}" \
    "${SC_OUTPUTS_DIR}" \
    "${REFERENCE_FA}" \
    "${SC_OUTPUTS_DIR}")

echo "Single-cell variant calling job ID: ${vcall_job_ID}"

######################################################################################################
#### Merge single-cell VCFs
MERGED_VCF="${RESULTS_DIR}/sc_variants_merged.vcf.gz"

merge_vcf_job_ID=$(sbatch --parsable \
    --dependency=afterok:${vcall_job_ID} \
    "${SCRIPTS_DIR}/scripts/CapGTA/merge_sc_vcfs.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${MERGED_VCF}")

echo "VCF merge job ID: ${merge_vcf_job_ID}"

######################################################################################################
#### Submit RNA count matrix generation
echo "Submitting RNA count matrix job..."

GTF="/shared/biodata/ngs/Reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"

count_job_ID=$(sbatch --parsable \
    --dependency=afterok:${rna_extraction_job_ID} \
    "${SCRIPTS_DIR}/scripts/CapGTA/create_rna_count_matrix.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${GTF}" \
    "${RESULTS_DIR}/rna_counts")

echo "Count matrix job ID: ${count_job_ID}"

echo ""
echo "================================"
echo "Pipeline submission complete!"
echo "================================"
echo "Final outputs:"
echo "  - SC DNA BAMs: ${SC_OUTPUTS_DIR}/*_dna.bam"
echo "  - SC RNA BAMs: ${SC_OUTPUTS_DIR}/*_rna.bam"
echo "  - SC VCFs: ${SC_OUTPUTS_DIR}/*.g.vcf.gz"
echo "  - Merged VCF: ${MERGED_VCF}"
echo "  - RNA count matrix: ${RESULTS_DIR}/rna_counts_matrix.csv"
echo ""
