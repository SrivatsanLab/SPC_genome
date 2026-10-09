#!/bin/bash
#SBATCH -J fc_pt_relaxed
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# featureCounts v2 for primer-trimmed BAMs with EVEN MORE relaxed gene-assignment knobs:
#
#   --fracOverlap 0.01         (was 0.10) — recover Overlapping_Length reads
#   --readExtension5 500       (new)      — treat 500 bp upstream of 5' gene edge as the gene
#   --readExtension3 500       (new)      — treat 500 bp downstream of 3' gene edge as the gene
#                                            — rescues ~21% of NoFeatures (reads within 500 bp of a gene)
#   --largestOverlap           (new)      — on multi-overlap, pick the gene with the longest overlap
#                                            — rescues Ambiguity (two-gene cases are the majority)
#
# All other settings mirror featurecounts_primer_trimmed.sh (-t transcript -g gene_id -p --countReadPairs
# --primary -R BAM --Rpath).
#
# Usage:
#   sbatch scripts/CapGTA/featurecounts_primer_trimmed_relaxed.sh <bams_dir_or_list> <annotation.gtf> <output_dir>
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bams_dir_or_list> <annotation.gtf> <output_dir>" >&2
    exit 1
fi

BAMS_INPUT="$1"
GTF="$2"
OUTPUT_DIR="$3"

if [ ! -f "${GTF}" ]; then
    echo "Error: GTF not found: ${GTF}" >&2
    exit 1
fi

if [ -d "${BAMS_INPUT}" ]; then
    mapfile -t BAM_PATHS < <(ls "${BAMS_INPUT}"/*.bam 2>/dev/null)
elif [ -f "${BAMS_INPUT}" ]; then
    mapfile -t BAM_PATHS < "${BAMS_INPUT}"
else
    echo "Error: input is neither a BAM directory nor a list file: ${BAMS_INPUT}" >&2
    exit 1
fi

N_BAMS="${#BAM_PATHS[@]}"
if [ "${N_BAMS}" -eq 0 ]; then
    echo "Error: no BAMs found under ${BAMS_INPUT}" >&2
    exit 1
fi

ANNO_DIR="${OUTPUT_DIR}/featurecounts_annotated"
mkdir -p "${ANNO_DIR}"

echo "=========================================="
echo "featureCounts v2 (relaxed knobs)"
echo "=========================================="
echo "BAM input:        ${BAMS_INPUT} (${N_BAMS} BAMs)"
echo "GTF:              ${GTF}"
echo "Output:           ${OUTPUT_DIR}"
echo "Annotated:        ${ANNO_DIR}"
echo "fracOverlap:      0.01 (was 0.10)"
echo "readExtension5:   500 bp"
echo "readExtension3:   500 bp"
echo "largestOverlap:   enabled"
echo "Started:          $(date -Iseconds)"
echo "=========================================="

module load Subread

featureCounts \
    -a "${GTF}" \
    -o "${OUTPUT_DIR}/counts_raw.txt" \
    -T 16 \
    -t transcript \
    -g gene_id \
    --fracOverlap 0.01 \
    --readExtension5 500 \
    --readExtension3 500 \
    --largestOverlap \
    --primary \
    -p --countReadPairs \
    -R BAM \
    --Rpath "${ANNO_DIR}" \
    "${BAM_PATHS[@]}"

echo ""
echo "Finished: $(date -Iseconds)"
echo "Annotated BAMs under: ${ANNO_DIR}"
