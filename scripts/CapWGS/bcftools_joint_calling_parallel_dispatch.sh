#!/bin/bash
#SBATCH -J bcf_jc_dispatch
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 10

##########################################################################################################################
# After interval generation completes, read the actual interval count and submit the
# bcftools per-interval array + concat-merge with the correct sizing.
#
# This exists because GATK SplitIntervals (with BALANCING_WITHOUT_INTERVAL_SUBDIVISION)
# may produce fewer bins than the requested SCATTER_COUNT, and the array size has to
# match the produced count or tasks above the count fail with empty interval lines.
##########################################################################################################################

set -euo pipefail

SCRIPTS_DIR="$1"
ALIGNED_DIR="$2"
REFERENCE_GENOME="$3"
PER_INTERVAL_DIR="$4"
INTERVAL_LIST="$5"
FINAL_VCF="$6"

if [ ! -f "${INTERVAL_LIST}" ]; then
    echo "ERROR: Interval list not found: ${INTERVAL_LIST}"
    exit 1
fi

n_intervals=$(wc -l < "${INTERVAL_LIST}")
if [ "${n_intervals}" -eq 0 ]; then
    echo "ERROR: Interval list is empty: ${INTERVAL_LIST}"
    exit 1
fi

echo "============================================"
echo "Dispatching bcftools array + merge"
echo "============================================"
echo "Intervals produced: ${n_intervals}"
echo ""

mkdir -p "${PER_INTERVAL_DIR}"

bcf_array_ID=$(sbatch --parsable \
    --array=1-${n_intervals} \
    "${SCRIPTS_DIR}/scripts/CapWGS/bcftools_joint_calling_parallel_array.sh" \
    "${ALIGNED_DIR}" \
    "${REFERENCE_GENOME}" \
    "${PER_INTERVAL_DIR}" \
    "${INTERVAL_LIST}")

echo "Per-interval array submitted: ${bcf_array_ID} (1-${n_intervals})"

merge_job_ID=$(sbatch --parsable \
    --dependency=afterok:${bcf_array_ID} \
    "${SCRIPTS_DIR}/scripts/CapWGS/gatk_joint_calling_parallel_merge.sh" \
    "${PER_INTERVAL_DIR}" \
    "${FINAL_VCF}" \
    "${INTERVAL_LIST}")

echo "Concat-merge submitted: ${merge_job_ID}"
echo ""
echo "Final VCF: ${FINAL_VCF}"
