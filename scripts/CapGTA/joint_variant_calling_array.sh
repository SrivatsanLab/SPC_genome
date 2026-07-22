#!/bin/bash
#SBATCH -J joint_vcall
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 12:00:00

###########################################################################################################################
# Joint single-cell variant calling — one genomic region per SLURM array task.
#
# Runs `bcftools mpileup` across ALL cells in a bam-list for the region assigned to this task,
# calls variants jointly with the multiallelic model, splits multiallelic sites and left-aligns
# indels. Output is one bgzipped, indexed VCF per region; downstream `concat_region_vcfs.sh`
# stitches them into a single joint VCF.
#
# Usage (submit as an array — array size == number of regions):
#   sbatch --array=1-$(wc -l < regions.txt) \
#       scripts/CapGTA/joint_variant_calling_array.sh \
#       <bam_list.txt> <regions.txt> <reference.fa> <output_dir>
#
# Inputs:
#   bam_list.txt   one BAM path per line (all coord-sorted, indexed, RG set with SM)
#   regions.txt    one `chrom:start-end` region per line (see make_regions.sh)
#   reference.fa   reference FASTA (must have .fai next to it)
#   output_dir     per-region VCFs are written here as region_<taskid>.vcf.gz
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: sbatch --array=1-N $0 <bam_list.txt> <regions.txt> <reference.fa> <output_dir>" >&2
    exit 1
fi

BAM_LIST="$1"
REGIONS_FILE="$2"
REFERENCE_FA="$3"
OUTPUT_DIR="$4"

for f in "${BAM_LIST}" "${REGIONS_FILE}" "${REFERENCE_FA}" "${REFERENCE_FA}.fai"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — this script must be submitted as an array." >&2
    exit 1
fi

REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${REGIONS_FILE}")
if [ -z "${REGION}" ]; then
    echo "Error: no region at line ${SLURM_ARRAY_TASK_ID} of ${REGIONS_FILE}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Zero-pad output name for stable lexical ordering
REGION_TAG=$(printf "region_%06d" "${SLURM_ARRAY_TASK_ID}")
OUTPUT_VCF="${OUTPUT_DIR}/${REGION_TAG}.vcf.gz"

# 1660+ BAMs opened at once — raise the file-descriptor limit.
ulimit -n 8192 || echo "Warning: could not raise ulimit -n" >&2

module load BCFtools SAMtools

N_BAMS=$(wc -l < "${BAM_LIST}")

echo "=========================================="
echo "Joint variant calling"
echo "=========================================="
echo "Task ID:    ${SLURM_ARRAY_TASK_ID}"
echo "Region:     ${REGION}"
echo "BAMs:       ${N_BAMS} (${BAM_LIST})"
echo "Reference:  ${REFERENCE_FA}"
echo "Output:     ${OUTPUT_VCF}"
echo "Started:    $(date -Iseconds)"
echo "=========================================="

# Joint mpileup → call → normalize.
# -a FORMAT/AD,FORMAT/DP: keep per-sample allele depths for downstream cellspec analysis.
# -q 20 -Q 20: minimum mapping and base quality.
# Duplicate reads are skipped by default (bcftools mpileup default skip-flags include BAM_FDUP).
# bcftools call -mv: multiallelic caller, variant sites only.
# bcftools norm -m -both: split multiallelic sites; -f enables left-align/normalization of indels.
bcftools mpileup \
    --threads 2 \
    -f "${REFERENCE_FA}" \
    -b "${BAM_LIST}" \
    -r "${REGION}" \
    -a FORMAT/AD,FORMAT/DP \
    -q 20 -Q 20 \
    -O u | \
bcftools call \
    --threads 2 \
    -mv \
    -O u | \
bcftools norm \
    --threads 2 \
    -f "${REFERENCE_FA}" \
    -m -both \
    -O z \
    -o "${OUTPUT_VCF}"

bcftools index "${OUTPUT_VCF}"

N_VARIANTS=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
echo "Finished:  $(date -Iseconds)"
echo "Variants:  ${N_VARIANTS}"
