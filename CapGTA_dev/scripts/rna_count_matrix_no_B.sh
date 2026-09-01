#!/bin/bash
#SBATCH -J rna_counts_no_B
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# featureCounts on a list of per-cell BAMs, WITHOUT the -B flag.
#
# The main-pipeline version (scripts/CapGTA/create_rna_count_matrix.sh) uses -B (both mates required).
# For BAMs that have been PRE-FILTERED (e.g. spliced-only, polyA/T-only, SL-only), -B discards any
# fragment whose partner mate was filtered out — a systematic undercount. This variant keeps
# fragments whose surviving mate still passes fracOverlap/primary/etc.
#
# Usage:
#   sbatch CapGTA_dev/scripts/rna_count_matrix_no_B.sh <bam_list.txt> <annotation.gtf> <output_prefix>
#
# Writes:
#   <output_prefix>_raw.txt         featureCounts raw output (with .summary)
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_list.txt> <annotation.gtf> <output_prefix>" >&2
    exit 1
fi

BAM_LIST="$1"
GTF="$2"
OUTPUT_PREFIX="$3"

for f in "${BAM_LIST}" "${GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "${OUTPUT_PREFIX}")"

N_BAMS=$(wc -l < "${BAM_LIST}")

echo "=========================================="
echo "featureCounts (no -B)"
echo "=========================================="
echo "BAM list:       ${BAM_LIST} (${N_BAMS} BAMs)"
echo "GTF:            ${GTF}"
echo "Output prefix:  ${OUTPUT_PREFIX}"
echo "Started:        $(date -Iseconds)"
echo "=========================================="

module load Subread

mapfile -t BAM_PATHS < "${BAM_LIST}"

featureCounts \
    -a "${GTF}" \
    -o "${OUTPUT_PREFIX}_raw.txt" \
    -T 16 \
    -t exon \
    -g gene_id \
    --fracOverlap 0.5 \
    --primary \
    -p \
    -C \
    "${BAM_PATHS[@]}"

echo ""
echo "Finished: $(date -Iseconds)"
