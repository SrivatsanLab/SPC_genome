#!/bin/bash
#SBATCH -J intergenic_ct
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH -t 30:00

###########################################################################################################################
# Count intergenic fragments per cell — one BAM per SLURM array task.
#
# Counts first-in-pair reads (each proper pair once) that overlap the intergenic BED,
# excluding duplicates / secondary / supplementary / QC-fail / unmapped:
#     samtools view -c -f 66 -F 3844 -L <intergenic.bed> <cell.bam>
# (-f 66 = paired + first-in-pair; -F 3844 excludes standard skip flags)
#
# The output is one TSV per cell: `<barcode>\t<intergenic_fragments>`
# The reducer (assemble_expression_matrix.py) turns these into per-cell lambdas.
#
# Usage:
#   sbatch --array=1-$(wc -l < bam_list.txt) \
#       scripts/CapGTA/count_intergenic_reads_array.sh \
#       <bam_list.txt> <intergenic.bed> <output_dir>
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch --array=1-N $0 <bam_list.txt> <intergenic.bed> <output_dir>" >&2
    exit 1
fi

BAM_LIST="$1"
INTERGENIC_BED="$2"
OUTPUT_DIR="$3"

for f in "${BAM_LIST}" "${INTERGENIC_BED}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: not found: ${f}" >&2
        exit 1
    fi
done

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BAM_LIST}")
if [ -z "${BAM}" ] || [ ! -f "${BAM}" ]; then
    echo "Error: BAM at line ${SLURM_ARRAY_TASK_ID} missing: '${BAM}'" >&2
    exit 1
fi

barcode=$(basename "${BAM}" .bam)
mkdir -p "${OUTPUT_DIR}"
OUT_TSV="${OUTPUT_DIR}/${barcode}.tsv"

module load SAMtools

count=$(samtools view --threads 2 -c -f 66 -F 3844 -L "${INTERGENIC_BED}" "${BAM}")

printf '%s\t%s\n' "${barcode}" "${count}" > "${OUT_TSV}"
echo "intergenic fragments for ${barcode}: ${count}"
