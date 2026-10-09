#!/bin/bash
#SBATCH -J fc_assign
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# featureCounts gene assignment for UMI-tagged per-cell BAMs.
#
# Diffs vs the classic create_rna_count_matrix.sh:
#   - -R BAM --Rpath <dir>:    emit per-cell annotated BAMs (reads tagged with XT:Z:<gene_id>)
#                              so umi_count.py can dedup by (cell, gene, UMI).
#   - -p (declare PE) but NO --countReadPairs / -B / -C:
#                              the source BAMs are paired-end so the per-read SAM flag has the PE
#                              bit set. We need -p for featureCounts to accept that, but we drop
#                              --countReadPairs so each kept read counts as 1 observation (the
#                              primer-mate's partner was discarded upstream and isn't in the BAM).
#   - No --ignoreDup:          UMI dedup replaces coordinate-based dup collapse. MarkDuplicates
#                              marks positional dups which conflates distinct UMIs at the same
#                              position; umi_count.py uses the UB tag directly.
#   - Keep --primary --fracOverlap 0.5: same as before.
#
# Also writes the standard gene × cell matrix of read counts (pre-dedup) via -o. That matrix is
# useful as a QC comparison vs the dedup'd UMI matrix, and is persisted in h5ad as layers['read_counts'].
#
# Usage:
#   sbatch scripts/CapGTA/featureCounts_assign.sh <bam_list.txt> <annotation.gtf> <output_dir>
#
# Writes:
#   <output_dir>/counts_raw.txt              # featureCounts summary matrix (gene × cell)
#   <output_dir>/counts_raw.txt.summary      # featureCounts per-cell QC
#   <output_dir>/featurecounts_annotated/<barcode>.bam.featureCounts.bam
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_list.txt> <annotation.gtf> <output_dir>" >&2
    exit 1
fi

BAM_LIST="$1"
GTF="$2"
OUTPUT_DIR="$3"

for f in "${BAM_LIST}" "${GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

ANNO_DIR="${OUTPUT_DIR}/featurecounts_annotated"
mkdir -p "${ANNO_DIR}"

N_BAMS=$(wc -l < "${BAM_LIST}")

echo "=========================================="
echo "featureCounts gene assignment (UMI mode)"
echo "=========================================="
echo "BAM list:       ${BAM_LIST} (${N_BAMS} BAMs)"
echo "GTF:            ${GTF}"
echo "Output dir:     ${OUTPUT_DIR}"
echo "Annotated BAMs: ${ANNO_DIR}"
echo "Started:        $(date -Iseconds)"
echo "=========================================="

module load Subread

mapfile -t BAM_PATHS < "${BAM_LIST}"

featureCounts \
    -a "${GTF}" \
    -o "${OUTPUT_DIR}/counts_raw.txt" \
    -T 16 \
    -t exon \
    -g gene_id \
    --fracOverlap 0.5 \
    --primary \
    -p \
    -R BAM \
    --Rpath "${ANNO_DIR}" \
    "${BAM_PATHS[@]}"

echo ""
echo "Finished featureCounts: $(date -Iseconds)"
echo "Annotated BAMs under: ${ANNO_DIR}"
