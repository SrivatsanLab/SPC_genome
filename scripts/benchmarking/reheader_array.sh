#!/bin/bash
#SBATCH -J reheader
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 30

# Strip the over-long ##GATKCommandLine line (~8330 chars) from each GVCF
# header so GATK htsjdk can parse it. Uses bcftools reheader which only
# rewrites the BGZF header blocks (data blocks copied unchanged).

set -euo pipefail

TASK_FILE="${TASK_FILE:-data/benchmarking/manifests/reheader_tasks.tsv}"
TASK_LINE=$(tail -n +2 "$TASK_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
IFS=$'\t' read -r DONOR_ID GVCF_PATH <<< "$TASK_LINE"
GVCF_PATH="${GVCF_PATH%$'\r'}"

if [ -z "$GVCF_PATH" ] || [ ! -f "$GVCF_PATH" ]; then
    echo "ERROR: GVCF not found at '$GVCF_PATH'" >&2
    exit 1
fi

module load BCFtools/1.18-GCC-12.2.0

WORK=/hpc/temp/srivatsan_s/reheader_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK"
HDR="$WORK/new_header.txt"

OUT="${GVCF_PATH%.g.vcf.gz}.fixed.g.vcf.gz"

echo "=== Reheader task: $SLURM_ARRAY_TASK_ID ==="
echo "Donor:  $DONOR_ID"
echo "Input:  $GVCF_PATH"
echo "Output: $OUT"

echo "Extracting header (excluding ##GATKCommandLine)..."
bcftools view -h "$GVCF_PATH" | grep -v '^##GATKCommandLine' > "$HDR"
echo "  Lines: $(wc -l < $HDR), longest: $(awk '{print length}' $HDR | sort -n | tail -1) chars"

echo "Reheadering..."
time bcftools reheader -h "$HDR" -o "$OUT" "$GVCF_PATH"

echo "Indexing..."
time tabix -p vcf "$OUT"

# Atomic-ish swap: rename original to .orig, fixed to canonical, then remove .orig
echo "Swapping in fixed file..."
mv "$GVCF_PATH" "${GVCF_PATH}.orig"
mv "${GVCF_PATH}.tbi" "${GVCF_PATH}.tbi.orig"
mv "$OUT" "$GVCF_PATH"
mv "${OUT}.tbi" "${GVCF_PATH}.tbi"
ls -la "${GVCF_PATH}" "${GVCF_PATH}.tbi"

# Quick sanity: bcftools view -h works on swapped file
N=$(bcftools view -h "$GVCF_PATH" | wc -l)
echo "Swapped file header lines: $N"

# Cleanup
rm -rf "$WORK"
echo "Done."
