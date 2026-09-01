#!/bin/bash
###########################################################################################################################
# CapGTA pseudobulk GATK pipeline — worm6 strain-germline / zygosity calling.
#
# Stages (chained via SLURM --dependency=afterok):
#   0. [inline] build_worm_lists.py       — worm_summary.csv + worm_lists/
#   1. [array]  merge_pseudobulk_array.sh — one task per worm: merge + spliced/polyA/T filter + RG collapse + depth
#   2. [array]  haplotype_caller_array.sh — worms × regions: HC in GVCF mode with per-worm --contamination-fraction-to-filter
#   3. [array]  combine_gvcfs_array.sh    — one task per region: CombineGVCFs + GenotypeGVCFs
#   4. [job]    concat_and_final.sh       — bcftools concat per-region joint VCFs → joint.vcf.gz
#   5. [job]    regenotype_from_ad.py     — contamination-aware GT from AD/DP + per-worm VAF histograms
#   6. [job]    analyze_worms_vs_subclones.py — private-variant zygosity, fused-pair detection, blacklist, REPORT.md
#
# Usage (run from repo root so SLURM_outs/ resolves next to the launch cwd):
#   bash scripts/CapGTA/CapGTA_pseudobulk_gatk_pipeline.sh \
#       <h5ad> <reference_fa> <regions_file> <existing_vcf> <output_dir> [--min-cells 50]
#
# Example (worm6):
#   bash scripts/CapGTA/CapGTA_pseudobulk_gatk_pipeline.sh \
#       results/worm6_final/DNA_analysis/joint_variants_process.h5ad \
#       /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
#       results/worm6_final/joint_variants/regions.txt \
#       results/worm6_final/joint_variants/joint_variants.vcf.gz \
#       results/worm6_final/pseudobulk_gatk
#
# Terminal outputs:
#   <output_dir>/joint.vcf.gz
#   <output_dir>/regenotyped.tsv
#   <output_dir>/analyses/REPORT.md
###########################################################################################################################

set -euo pipefail

MIN_CELLS=50

POS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --min-cells)  MIN_CELLS="$2"; shift 2 ;;
        --)           shift; POS+=("$@"); break ;;
        *)            POS+=("$1"); shift ;;
    esac
done

if [ "${#POS[@]}" -ne 5 ]; then
    echo "Usage: $0 <h5ad> <reference_fa> <regions_file> <existing_vcf> <output_dir> [--min-cells N]" >&2
    exit 1
fi

H5AD="${POS[0]}"
REF="${POS[1]}"
REGIONS_FILE="${POS[2]}"
EXISTING_VCF="${POS[3]}"
OUT="${POS[4]}"

# The stage scripts live in a sibling `pseudobulk_gatk/` directory next to this file.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PB="${SCRIPT_DIR}/pseudobulk_gatk"

# Capture the currently-active micromamba env — SLURM jobs don't inherit it and
# need to re-activate before running any Python step. Same convention as
# scripts/CapWGS/PP_array.sh.
MMENV="${CONDA_DEFAULT_ENV:-${MAMBA_DEFAULT_ENV:-}}"
if [ -z "${MMENV}" ] || [ "${MMENV}" = "base" ]; then
    echo "Error: no micromamba env active (CONDA_DEFAULT_ENV empty or 'base')." >&2
    echo "       Activate the env with cellspec/anndata installed, then rerun:" >&2
    echo "         micromamba activate <env>" >&2
    exit 1
fi
echo "Micromamba env: ${MMENV}"

for f in "${H5AD}" "${REF}" "${REF}.fai" "${REGIONS_FILE}" "${EXISTING_VCF}"; do
    [ -f "${f}" ] || { echo "Error: required input not found: ${f}" >&2; exit 1; }
done

mkdir -p "${OUT}" SLURM_outs SLURM_outs/array_outs

echo "=========================================="
echo "CapGTA pseudobulk GATK pipeline"
echo "=========================================="
echo "H5AD:          ${H5AD}"
echo "Reference:     ${REF}"
echo "Regions:       ${REGIONS_FILE}"
echo "Existing VCF:  ${EXISTING_VCF}"
echo "Output dir:    ${OUT}"
echo "Min cells:     ${MIN_CELLS}"
echo "=========================================="

# --- Stage 0: build worm lists (inline; needed to know array shapes for stages 1-2) ------------
echo ""
echo "[Stage 0] Building worm lists (inline)…"
python3 "${PB}/build_worm_lists.py" "${H5AD}" "${OUT}" --min-cells "${MIN_CELLS}"

N_WORMS=$(wc -l < "${OUT}/worm_lists/passing_worms.txt")
N_REG=$(wc -l < "${REGIONS_FILE}")
N_HC=$((N_WORMS * N_REG))

if [ "${N_WORMS}" -lt 1 ] || [ "${N_REG}" -lt 1 ]; then
    echo "Error: N_WORMS=${N_WORMS} N_REG=${N_REG}" >&2
    exit 1
fi

echo "  N_WORMS=${N_WORMS}  N_REG=${N_REG}  HC tasks=${N_HC}"

# --- Stage 1: merge pseudobulk (one task per worm) --------------------------------------------
echo ""
echo "[Stage 1] Submitting merge_pseudobulk_array (${N_WORMS} tasks)…"
J1=$(sbatch --parsable \
        --array=1-"${N_WORMS}" \
        "${PB}/merge_pseudobulk_array.sh" "${OUT}")
echo "  job ${J1}"

# --- Stage 2: HaplotypeCaller (worms × regions) -----------------------------------------------
echo ""
echo "[Stage 2] Submitting haplotype_caller_array (${N_HC} tasks)…"
J2=$(sbatch --parsable \
        --dependency=afterok:"${J1}" \
        --array=1-"${N_HC}" \
        "${PB}/haplotype_caller_array.sh" "${OUT}" "${REF}" "${REGIONS_FILE}")
echo "  job ${J2} (depends on ${J1})"

# --- Stage 3: CombineGVCFs + GenotypeGVCFs (one task per region) ------------------------------
echo ""
echo "[Stage 3] Submitting combine_gvcfs_array (${N_REG} tasks)…"
J3=$(sbatch --parsable \
        --dependency=afterok:"${J2}" \
        --array=1-"${N_REG}" \
        "${PB}/combine_gvcfs_array.sh" "${OUT}" "${REF}" "${REGIONS_FILE}")
echo "  job ${J3} (depends on ${J2})"

# --- Stage 4: concat per-region joint VCFs → joint.vcf.gz -------------------------------------
echo ""
echo "[Stage 4] Submitting concat_and_final…"
J4=$(sbatch --parsable \
        --dependency=afterok:"${J3}" \
        "${PB}/concat_and_final.sh" "${OUT}" "${REGIONS_FILE}")
echo "  job ${J4} (depends on ${J3})"

# --- Stage 5: contamination-aware regenotyping + VAF hists ------------------------------------
# Use dedicated .sh wrappers (not --wrap): sbatch --wrap writes a /bin/sh script, but we need
# `module load` (a bash function from /etc/profile.d/modules.sh). The wrapper scripts have
# `#!/bin/bash` shebangs so module + micromamba both work.
echo ""
echo "[Stage 5] Submitting regenotype_from_ad…"
J5=$(sbatch --parsable \
        --dependency=afterok:"${J4}" \
        "${PB}/regenotype_job.sh" "${OUT}" "${MMENV}")
echo "  job ${J5} (depends on ${J4}, env=${MMENV})"

# --- Stage 6: downstream analyses ------------------------------------------------------------
echo ""
echo "[Stage 6] Submitting analyze_worms_vs_subclones…"
J6=$(sbatch --parsable \
        --dependency=afterok:"${J5}" \
        "${PB}/analyze_job.sh" "${OUT}" "${EXISTING_VCF}" "${MMENV}")
echo "  job ${J6} (depends on ${J5}, env=${MMENV})"

echo ""
echo "=========================================="
echo "Pipeline submitted. Job chain:"
echo "  ${J1} (merge, ${N_WORMS} tasks)"
echo "  → ${J2} (HC, ${N_HC} tasks)"
echo "  → ${J3} (combine+genotype, ${N_REG} tasks)"
echo "  → ${J4} (concat)"
echo "  → ${J5} (regenotype)"
echo "  → ${J6} (analyze)"
echo ""
echo "Monitor:   squeue -u \$USER"
echo "Terminal:  ${OUT}/joint.vcf.gz"
echo "           ${OUT}/regenotyped.tsv"
echo "           ${OUT}/analyses/REPORT.md"
echo "=========================================="
