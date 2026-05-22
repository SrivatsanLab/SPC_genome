#!/bin/bash
#SBATCH -J aneufinder
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00

# Per-sample AneuFinder CNV calling.
# Array index = row in TASK_LIST (1-indexed, header stripped).
#
# TASK_LIST format (TSV, header):
#   donor_id  sample_id  bam_path  is_bulk
#
# Usage:
#   sbatch --array=1-N aneufinder_array.sh <task_list.tsv>
#
# Outputs land in data/benchmarking/cnv/<donor>/<sample>.{rds,segments.tsv,bins.tsv.gz,qc.tsv}

set -euo pipefail

TASK_LIST="${1:?usage: aneufinder_array.sh <task_list.tsv>}"
BINSIZE="${2:-1000000}"
STEPSIZE="${3:-1000000}"
OUT_ROOT="${4:-data/benchmarking/cnv}"

# Pick task row (skip header)
ROW=$(awk -v idx=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==idx' "$TASK_LIST")
if [ -z "$ROW" ]; then
    echo "ERROR: no row for array index $SLURM_ARRAY_TASK_ID in $TASK_LIST" >&2
    exit 1
fi
DONOR=$(echo "$ROW"  | awk -F'\t' '{print $1}')
SAMPLE=$(echo "$ROW" | awk -F'\t' '{print $2}')
BAM=$(echo "$ROW"    | awk -F'\t' '{print $3}')
IS_BULK=$(echo "$ROW"| awk -F'\t' '{print $4}')

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found: $BAM" >&2
    exit 1
fi

OUT_DIR="${OUT_ROOT}/${DONOR}"
mkdir -p "$OUT_DIR"

# Skip if all outputs already exist (so reruns are cheap)
RDS="${OUT_DIR}/${SAMPLE}.aneufinder.rds"
QC="${OUT_DIR}/${SAMPLE}.qc.tsv"
SEG="${OUT_DIR}/${SAMPLE}.segments.tsv"
BINS="${OUT_DIR}/${SAMPLE}.bins.tsv.gz"
if [ -s "$RDS" ] && [ -s "$QC" ] && [ -s "$SEG" ] && [ -s "$BINS" ]; then
    echo "All outputs exist for $SAMPLE, skipping."
    exit 0
fi

# Detect paired-end status from first ~10k reads (flag 0x1 = paired)
module load SAMtools/1.19.2-GCC-13.2.0 2>/dev/null || module load SAMtools
N_PAIRED=$(samtools view -c -f 1 "$BAM" 2>/dev/null | head -1)
if [ "${N_PAIRED:-0}" -gt 0 ]; then
    PAIRED=TRUE
else
    PAIRED=FALSE
fi
module unload SAMtools 2>/dev/null || true

echo "=== AneuFinder task $SLURM_ARRAY_TASK_ID ==="
echo "Donor:    $DONOR"
echo "Sample:   $SAMPLE"
echo "BAM:      $BAM"
echo "Bulk:     $IS_BULK"
echo "Paired:   $PAIRED"
echo "Binsize:  $BINSIZE"
echo "Stepsize: $STEPSIZE"
echo "Out dir:  $OUT_DIR"
echo

module load fhR

time Rscript scripts/benchmarking/cnv/aneufinder_one_sample.R \
    --bam "$BAM" \
    --sample "$SAMPLE" \
    --out-dir "$OUT_DIR" \
    --paired-end "$PAIRED" \
    --binsize "$BINSIZE" \
    --stepsize "$STEPSIZE" \
    --chromosomes autosomes+X

module unload fhR 2>/dev/null || true

echo
echo "Done: $SAMPLE"
