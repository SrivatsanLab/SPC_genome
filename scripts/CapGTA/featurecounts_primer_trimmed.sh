#!/bin/bash
#SBATCH -J fc_primer_trim
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# featureCounts for primer-trimmed paired-end BAMs with relaxed gene-overlap settings.
#
# Differences vs featureCounts_assign.sh:
#   -t transcript (not exon):   count reads overlapping the WHOLE gene body (exon + intron). The iGenomes
#                               WBcel235 GTF lacks a `gene` feature row; `transcript` entries span the full
#                               genomic locus for each transcript, and meta-feature aggregation by -g gene_id
#                               ensures each read is counted at most once per gene.
#   --fracOverlap 0.1 (not 0.5): lenient to short cDNA reads that only partially overlap a gene.
#   -p --countReadPairs:        each kept PAIR is one observation (UMI dedup happens downstream).
#   -R BAM --Rpath:             per-read gene assignment preserved so build_umi_matrix.py can read XT/XS.
#
# Usage:
#   sbatch scripts/CapGTA/featurecounts_primer_trimmed.sh <bams_dir_or_list> <annotation.gtf> <output_dir>
#
# <bams_dir_or_list> — directory containing <cell>.bam files, OR a text file with one BAM path per line.
#
# Writes:
#   <output_dir>/counts_raw.txt
#   <output_dir>/counts_raw.txt.summary
#   <output_dir>/featurecounts_annotated/<cell>.bam.featureCounts.bam   (per-read XS/XT tags)
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_list.txt> <annotation.gtf> <output_dir>" >&2
    exit 1
fi

BAMS_INPUT="$1"
GTF="$2"
OUTPUT_DIR="$3"

if [ ! -f "${GTF}" ]; then
    echo "Error: GTF not found: ${GTF}" >&2
    exit 1
fi

# Build BAM path array from either a directory or a list file
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
echo "featureCounts (primer-trimmed, relaxed)"
echo "=========================================="
echo "BAM input:  ${BAMS_INPUT} (${N_BAMS} BAMs)"
echo "GTF:        ${GTF}"
echo "Output:     ${OUTPUT_DIR}"
echo "Annotated:  ${ANNO_DIR}"
echo "Started:    $(date -Iseconds)"
echo "=========================================="

module load Subread

featureCounts \
    -a "${GTF}" \
    -o "${OUTPUT_DIR}/counts_raw.txt" \
    -T 16 \
    -t transcript \
    -g gene_id \
    --fracOverlap 0.1 \
    --primary \
    -p --countReadPairs \
    -R BAM \
    --Rpath "${ANNO_DIR}" \
    "${BAM_PATHS[@]}"

echo ""
echo "Finished: $(date -Iseconds)"
echo "Annotated BAMs under: ${ANNO_DIR}"
