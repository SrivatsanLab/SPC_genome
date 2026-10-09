#!/bin/bash
###########################################################################################################################
# Full SoupX pipeline for the polyA/T RNA count matrices.
#
# For a given mode (umi | reads), chain:
#   1. submit_rna_counts_umi.sh or submit_rna_counts_polyAT.sh with empty_cells.csv → empties h5ad
#   2. prep_soupx_inputs.py: processed cells h5ad + empties h5ad → soupx_in/ (mtx + clusters)
#   3. soupx_decontaminate.R → soupx_out/
#   4. reimport_soupx.py → <mode>_soupx.h5ad
#   5. preprocess_scanpy.py → <mode>_soupx_processed.h5ad
#
# All stages are sbatch'd with afterok dependencies.
#
# Usage:
#   scripts/CapGTA/submit_soupx_polyAT.sh <mode>
#     mode = umi | reads
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <mode:umi|reads>" >&2
    exit 1
fi

MODE="$1"
case "${MODE}" in
    umi|reads) ;;
    *) echo "Error: mode must be 'umi' or 'reads' (got '${MODE}')" >&2; exit 1 ;;
esac

REPO="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome"
SCRIPTS_DIR="${REPO}/scripts/CapGTA"
EMPTY_CSV="${REPO}/results/worm6_final/GEX/soupx_empties_polyAT/empty_cells.csv"
GTF="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"

if [ "${MODE}" = "umi" ]; then
    CELLS_PROCESSED="${REPO}/results/worm6_final/GEX/adata_rna_polyAT_umi_processed.h5ad"
    EMPTIES_DIR="${REPO}/results/worm6_final/rna_umi_counts/L6_empties"
    EMPTIES_H5AD="${EMPTIES_DIR}/rna_umi_matrix.h5ad"
    SOUPX_ROOT="${REPO}/results/worm6_final/GEX/soupx_polyAT_umi"
    SOUPX_H5AD="${REPO}/results/worm6_final/GEX/adata_rna_polyAT_umi_soupx.h5ad"
    SOUPX_PROCESSED="${REPO}/results/worm6_final/GEX/adata_rna_polyAT_umi_soupx_processed.h5ad"
    EMPTIES_SUBMIT="${SCRIPTS_DIR}/submit_rna_counts_umi.sh"
else
    CELLS_PROCESSED="${REPO}/results/worm6_final/GEX/adata_rna_polyAT_reads_processed.h5ad"
    EMPTIES_DIR="${REPO}/results/worm6_final/rna_umi_counts/L6_polyAT_empties"
    EMPTIES_H5AD="${EMPTIES_DIR}/rna_polyAT_matrix.h5ad"
    SOUPX_ROOT="${REPO}/results/worm6_final/GEX/soupx_polyAT_reads"
    SOUPX_H5AD="${REPO}/results/worm6_final/GEX/adata_rna_polyAT_reads_soupx.h5ad"
    SOUPX_PROCESSED="${REPO}/results/worm6_final/GEX/adata_rna_polyAT_reads_soupx_processed.h5ad"
    EMPTIES_SUBMIT="${SCRIPTS_DIR}/submit_rna_counts_polyAT.sh"
fi

for f in "${CELLS_PROCESSED}" "${EMPTY_CSV}" "${GTF}" "${EMPTIES_SUBMIT}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

mkdir -p "${EMPTIES_DIR}" "${SOUPX_ROOT}" SLURM_outs SLURM_outs/array_outs

echo "=========================================="
echo "SoupX pipeline on polyAT ${MODE} matrix"
echo "=========================================="
echo "Cells processed:  ${CELLS_PROCESSED}"
echo "Empties CSV:      ${EMPTY_CSV}"
echo "Empties out:      ${EMPTIES_H5AD}"
echo "SoupX root:       ${SOUPX_ROOT}"
echo "Final processed:  ${SOUPX_PROCESSED}"
echo "=========================================="

# --- 1) Build empties counts matrix via the normal counting pipeline ----------
# The nested submitter prints "Submitted ... job: NNN" for each sub-step and chains
# them internally; we capture the LAST (count+h5ad) job id to depend on here.
EMPTIES_LOG=$("${EMPTIES_SUBMIT}" "${EMPTY_CSV}" "${GTF}" "${EMPTIES_DIR}" 2>&1)
echo "${EMPTIES_LOG}"
EMPTIES_LAST_JOB=$(echo "${EMPTIES_LOG}" | grep -oE 'job:[[:space:]]+[0-9]+' | tail -1 | grep -oE '[0-9]+')
if [ -z "${EMPTIES_LAST_JOB}" ]; then
    echo "Error: could not parse final empties job id from submitter output" >&2
    exit 1
fi
echo "[chain] empties count+h5ad job: ${EMPTIES_LAST_JOB}"

# --- 2) prep_soupx_inputs: assemble cells + empties into soupx_in/ ------------
prep_job=$(sbatch --parsable \
    --dependency="afterok:${EMPTIES_LAST_JOB}" \
    -J soupx_prep_${MODE} \
    -o SLURM_outs/%x_%j.out \
    -p campus-new -c 2 --mem=24G -t 30:00 \
    --wrap="eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate default_jupyter && \
            python '${SCRIPTS_DIR}/prep_soupx_inputs.py' \
                '${CELLS_PROCESSED}' '${EMPTIES_H5AD}' '${SOUPX_ROOT}/soupx_in' \
                --cluster-key leiden")
echo "[chain] prep_soupx_inputs job:   ${prep_job} (depends on ${EMPTIES_LAST_JOB})"

# --- 3) soupx_decontaminate.R -------------------------------------------------
soupx_job=$(sbatch --parsable \
    --dependency="afterok:${prep_job}" \
    -J soupx_decontam_${MODE} \
    -o SLURM_outs/%x_%j.out \
    -p campus-new -c 2 --mem=24G -t 30:00 \
    --wrap="module load fhR && \
            Rscript '${SCRIPTS_DIR}/soupx_decontaminate.R' \
                '${SOUPX_ROOT}/soupx_in' '${SOUPX_ROOT}/soupx_out'")
echo "[chain] soupx_decontaminate job: ${soupx_job} (depends on ${prep_job})"

# --- 4) reimport_soupx --------------------------------------------------------
reimport_job=$(sbatch --parsable \
    --dependency="afterok:${soupx_job}" \
    -J soupx_reimport_${MODE} \
    -o SLURM_outs/%x_%j.out \
    -p campus-new -c 2 --mem=16G -t 15:00 \
    --wrap="eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate default_jupyter && \
            python '${SCRIPTS_DIR}/reimport_soupx.py' \
                '${CELLS_PROCESSED}' '${SOUPX_ROOT}/soupx_out' '${SOUPX_H5AD}'")
echo "[chain] reimport_soupx job:      ${reimport_job} (depends on ${soupx_job})"

# --- 5) preprocess_scanpy on the SoupX-decontaminated h5ad --------------------
pp_job=$(sbatch --parsable \
    --dependency="afterok:${reimport_job}" \
    "${SCRIPTS_DIR}/preprocess_polyAT.sh" \
    "${SOUPX_H5AD}" "${SOUPX_PROCESSED}")
echo "[chain] preprocess_scanpy job:   ${pp_job} (depends on ${reimport_job})"

echo ""
echo "Final SoupX-decontaminated + preprocessed h5ad will land at:"
echo "  ${SOUPX_PROCESSED}"
echo ""
echo "Monitor with:  squeue -u \$USER | grep -E 'soupx|preprocess|rna_counts|umi_filter|polyAT'"
