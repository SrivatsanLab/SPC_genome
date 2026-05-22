#!/bin/bash
#SBATCH -J per_cell_fpr
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 12:00:00

# Per-cell FPR vs matched bulk WGS truth.
# Array task = donor index in DONORS array below.

set -euo pipefail

# Activate the user's micromamba env that has pysam/pandas
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

DONORS=(HSC4 CD34_cord_blood BJ_fibroblast)
DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID - 1))]}"

JOINT_VCF="data/benchmarking/joint/${DONOR}/${DONOR}_joint.vcf.gz"

case "$DONOR" in
    HSC4)
        BULK_BAMS=("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/HSC4_bulk/bams/bulkHSC4.bqsr.bam")
        ;;
    CD34_cord_blood)
        BULK_BAMS=("data/benchmarking/bulks/CD34_cord_blood/SRR8438271.bqsr.bam")
        ;;
    BJ_fibroblast)
        BULK_BAMS=("data/benchmarking/bulks/BJ_fibroblast/SRR5365377.bqsr.bam"
                   "data/benchmarking/bulks/BJ_fibroblast/SRR5365378.bqsr.bam")
        ;;
    *)
        echo "ERROR: unknown donor $DONOR" >&2; exit 1 ;;
esac

OUT_DIR="data/benchmarking/fpr"
mkdir -p "$OUT_DIR"
OUT_TSV="${OUT_DIR}/${DONOR}_per_cell_fpr.tsv"
PER_CALL_TSV="${OUT_DIR}/${DONOR}_per_call.tsv.gz"

echo "=== per_cell_fpr array task $SLURM_ARRAY_TASK_ID ==="
echo "Donor:        $DONOR"
echo "Joint VCF:    $JOINT_VCF"
echo "Bulks:        ${BULK_BAMS[*]}"
echo "Summary TSV:  $OUT_TSV"
echo "Per-call TSV: $PER_CALL_TSV"
echo

time python scripts/benchmarking/per_cell_fpr.py \
    --donor "$DONOR" \
    --joint-vcf "$JOINT_VCF" \
    --bulk-bam "${BULK_BAMS[@]}" \
    --output "$OUT_TSV" \
    --per-call-tsv "$PER_CALL_TSV" \
    --min-alt 2 \
    --min-bulk-dp 10 \
    --min-mq 20 \
    --min-bq 20

echo
echo "Done."
