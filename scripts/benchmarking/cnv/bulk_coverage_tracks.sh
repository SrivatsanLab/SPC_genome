#!/bin/bash
#SBATCH -J bulkcov
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 6:00:00

# Build 1Mb-bin coverage bigWigs for each bulk BAM using deepTools bamCoverage.
# Drop-in replacement for AneuFinder bulk runs (which mis-call CN at high depth).
# Output is CPM-normalized, MarkDup'd reads excluded, MAPQ>=20.
#
# Usage:
#   sbatch --array=1-4 bulk_coverage_tracks.sh <task_list.tsv>
#
# TASK_LIST format (TSV, header):
#   donor_id  sample_id  bam_path
#
# Outputs land in data/benchmarking/cnv/<donor>/<sample>.coverage.1Mb.bw

set -euo pipefail

TASK_LIST="${1:?usage: bulk_coverage_tracks.sh <task_list.tsv>}"
BINSIZE="${2:-1000000}"
OUT_ROOT="${3:-data/benchmarking/cnv}"

ROW=$(awk -v idx=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==idx' "$TASK_LIST")
[ -z "$ROW" ] && { echo "ERROR: no row for idx $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
DONOR=$(echo "$ROW"  | awk -F'\t' '{print $1}')
SAMPLE=$(echo "$ROW" | awk -F'\t' '{print $2}')
BAM=$(echo "$ROW"    | awk -F'\t' '{print $3}')
[ -f "$BAM" ] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }

OUT_DIR="${OUT_ROOT}/${DONOR}"
mkdir -p "$OUT_DIR"
BW="${OUT_DIR}/${SAMPLE}.coverage.${BINSIZE}.bw"

if [ -s "$BW" ]; then
    echo "Already exists: $BW"; exit 0
fi

module purge >/dev/null 2>&1
module load deepTools/3.5.4.post1-gfbf-2022b

echo "=== bamCoverage task $SLURM_ARRAY_TASK_ID ==="
echo "Donor:    $DONOR"
echo "Sample:   $SAMPLE"
echo "BAM:      $BAM"
echo "Binsize:  $BINSIZE"
echo "Out:      $BW"
echo

time bamCoverage \
    --bam "$BAM" \
    --outFileName "$BW" \
    --outFileFormat bigwig \
    --binSize "$BINSIZE" \
    --normalizeUsing CPM \
    --minMappingQuality 20 \
    --ignoreDuplicates \
    --numberOfProcessors 16

echo "Done: $BW"
