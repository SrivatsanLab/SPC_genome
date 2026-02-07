#!/bin/bash
#SBATCH -J submit_sc_var
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

###########################################################################################################################
# Wrapper script that submits single cell variant calling array job
# This script runs AFTER concatenation completes and reads real_cells.txt to determine array size
###########################################################################################################################

set -euo pipefail

TMP_DIR="$1"
REAL_CELLS_FILE="$2"
REFERENCE_GENOME="$3"
OUTPUT_DIR="$4"
SCRIPTS_DIR="$5"
VARIANT_CALLER="${6:-gatk}"  # Default to gatk for backward compatibility

# Verify real_cells.txt exists
if [ ! -f "$REAL_CELLS_FILE" ]; then
    echo "ERROR: real_cells.txt not found at: $REAL_CELLS_FILE"
    exit 1
fi

# Count cells
cell_count=$(wc -l < "$REAL_CELLS_FILE")

echo "Found ${cell_count} cells in ${REAL_CELLS_FILE}"

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells detected. Cannot submit variant calling array."
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Submit appropriate variant calling array based on caller
if [ "$VARIANT_CALLER" = "gatk" ]; then
    echo "Using GATK HaplotypeCaller for variant calling (GVCF mode)"

    # Submit GATK HaplotypeCaller array (generates GVCFs from preprocessed single-cell BAMs)
    sc_var_array_ID=$(sbatch --parsable \
        --array=1-${cell_count} \
        "${SCRIPTS_DIR}/scripts/CapWGS/sc_var_array_gatk.sh" \
        "${REAL_CELLS_FILE}" \
        "${REFERENCE_GENOME}" \
        "${OUTPUT_DIR}")

    echo "GATK HaplotypeCaller array submitted: ${sc_var_array_ID}"
    echo "Array size: 1-${cell_count}"
    echo "Output: GVCF files for joint calling"
elif [ "$VARIANT_CALLER" = "bcftools" ]; then
    echo "Using BCFtools for variant calling"

    # Submit BCFtools mpileup/call array (generates VCFs)
    sc_var_array_ID=$(sbatch --parsable \
        --array=1-${cell_count} \
        "${SCRIPTS_DIR}/scripts/CapWGS/sc_var_array_bcftools.sh" \
        "${REAL_CELLS_FILE}" \
        "${REFERENCE_GENOME}" \
        "${OUTPUT_DIR}")

    echo "BCFtools variant calling array submitted: ${sc_var_array_ID}"
    echo "Array size: 1-${cell_count}"
    echo "Output: VCF files for merging"
else
    echo "ERROR: Invalid variant caller: ${VARIANT_CALLER}"
    echo "Must be 'gatk' or 'bcftools'"
    exit 1
fi
