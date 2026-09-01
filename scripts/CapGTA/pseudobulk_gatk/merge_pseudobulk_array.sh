#!/bin/bash
#SBATCH -J merge_pb
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 8:00:00

###########################################################################################################################
# Merge per-cell BAMs into one pseudobulk per worm, strip RNA-artifact reads, collapse RG, index, measure depth.
#
# Per SLURM array task (one worm):
#   1. samtools merge -b <worm>.bampaths -> raw.bam
#   2. samtools view -e '!(cigar =~ "N") && !(seq =~ "A{20,}") && !(seq =~ "T{20,}")' -> filtered.bam
#      (drops spliced reads and poly-A/T reads that would seed spurious HC reassemblies)
#   3. samtools addreplacerg to rewrite all RGs to SM=<worm> (one VCF column per worm)
#   4. samtools index
#   5. samtools depth summary -> depth_report/<worm>.txt
#
# Duplicates are NOT re-marked at the worm level: cross-cell "duplicates" from different WGA
# capsules are independent observations, not PCR duplicates.
#
# Usage:
#   sbatch --array=1-$(wc -l < <output_dir>/worm_lists/passing_worms.txt) \
#       scripts/CapGTA/pseudobulk_gatk/merge_pseudobulk_array.sh <output_dir>
#
# Writes:
#   <output_dir>/worm_bams/<worm>.bam[.bai]
#   <output_dir>/depth_report/<worm>.txt
#   <output_dir>/filter_report.tsv   (appended: worm  n_in  n_out  frac_kept)
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch --array=1-N $0 <output_dir>" >&2
    exit 1
fi

OUT="$1"
WORM_LIST="${OUT}/worm_lists/passing_worms.txt"

if [ ! -f "${WORM_LIST}" ]; then
    echo "Error: passing_worms.txt not found: ${WORM_LIST}" >&2
    exit 1
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

WORM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORM_LIST}")
if [ -z "${WORM}" ]; then
    echo "Error: no worm at line ${SLURM_ARRAY_TASK_ID} of ${WORM_LIST}" >&2
    exit 1
fi

BAMPATHS="${OUT}/worm_lists/${WORM}.bampaths"
if [ ! -s "${BAMPATHS}" ]; then
    echo "Error: bampaths file missing or empty: ${BAMPATHS}" >&2
    exit 1
fi

WORM_BAM_DIR="${OUT}/worm_bams"
DEPTH_DIR="${OUT}/depth_report"
mkdir -p "${WORM_BAM_DIR}" "${DEPTH_DIR}"

RAW_BAM="${WORM_BAM_DIR}/${WORM}.raw.bam"
FILT_BAM="${WORM_BAM_DIR}/${WORM}.filt.bam"
OUT_BAM="${WORM_BAM_DIR}/${WORM}.bam"
DEPTH_OUT="${DEPTH_DIR}/${WORM}.txt"

N_CELLS=$(wc -l < "${BAMPATHS}")

echo "=========================================="
echo "merge_pseudobulk  worm=${WORM}  cells=${N_CELLS}"
echo "task_id=${SLURM_ARRAY_TASK_ID}  started=$(date -Iseconds)"
echo "=========================================="

module load SAMtools

# 1660+ open file descriptors during merge — raise the limit as in joint_variant_calling_array.sh.
ulimit -n 8192 || echo "Warning: could not raise ulimit -n" >&2

# 1. Merge
echo ""
echo "[1/5] samtools merge…"
samtools merge -@ 8 -f -b "${BAMPATHS}" "${RAW_BAM}"

# 2. Filter: drop spliced reads (CIGAR N), poly-A (>=20 A) and poly-T (>=20 T) reads.
#    Requires samtools >= 1.15 (POSIX-ERE filter expressions).
echo ""
echo "[2/5] samtools view (filter spliced + polyA/T)…"
samtools view -@ 8 -b -h \
    -e '!(cigar =~ "N") && !(seq =~ "A{20,}") && !(seq =~ "T{20,}")' \
    -o "${FILT_BAM}" "${RAW_BAM}"

N_IN=$(samtools view -c "${RAW_BAM}")
N_OUT=$(samtools view -c "${FILT_BAM}")
FRAC=$(awk -v a="${N_OUT}" -v b="${N_IN}" 'BEGIN{printf "%.4f", (b>0)?a/b:0}')
printf "%s\t%d\t%d\t%s\n" "${WORM}" "${N_IN}" "${N_OUT}" "${FRAC}" \
    >> "${OUT}/filter_report.tsv"
echo "  kept ${N_OUT} / ${N_IN} (${FRAC})"

rm "${RAW_BAM}"

# 3. Collapse all per-cell RGs to one worm-level RG so the joint VCF gets ONE column per worm.
echo ""
echo "[3/5] samtools addreplacerg (SM=${WORM})…"
samtools addreplacerg -@ 8 \
    -r "ID:${WORM}\tSM:${WORM}\tLB:${WORM}\tPL:ILLUMINA" \
    -m overwrite_all \
    -o "${OUT_BAM}" "${FILT_BAM}"

rm "${FILT_BAM}"

# 4. Index
echo ""
echo "[4/5] samtools index…"
samtools index -@ 8 "${OUT_BAM}"

# 5. Depth summary (skip flagged duplicates and secondary alignments; MAPQ >= 20).
echo ""
echo "[5/5] samtools depth summary…"
samtools depth -a -Q 20 -G 0x400 "${OUT_BAM}" \
    | awk -v w="${WORM}" '
        {s+=$3; n++; if($3>0){c++; cs+=$3}}
        END {
            printf "worm\t%s\n", w;
            printf "mean_depth_all_sites\t%.4f\n", (n>0)?s/n:0;
            printf "breadth_frac_covered\t%.4f\n", (n>0)?c/n:0;
            printf "mean_depth_covered_sites\t%.4f\n", (c>0)?cs/c:0;
            printf "n_sites\t%d\n", n;
        }' > "${DEPTH_OUT}"

cat "${DEPTH_OUT}"

echo ""
echo "=========================================="
echo "Done: ${WORM}  finished=$(date -Iseconds)"
echo "Output: ${OUT_BAM}"
echo "=========================================="
