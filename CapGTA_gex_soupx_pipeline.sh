#!/bin/bash
###########################################################################################################################
# Full GEX + SoupX pipeline: raw BAMs → SoupX-decontaminated preprocessed h5ad.
#
# Stages (chained via SLURM afterok):
#   0. [inline] Build per-gene splice-junction table from the GTF.
#   1. [array]  Real-cell spliced counts (submit_rna_counts.sh).
#   2. [job]    Preprocess real cells → adata_gex_spliced.h5ad (Leiden clusters for SoupX).
#   3. [inline] Enumerate empty droplets (random subset of non-real BAMs).
#   4. [array]  Empty-droplet spliced counts (submit_rna_counts.sh).
#   5. [job]    Prep SoupX inputs (WGA-aware QC on empties + gene-order alignment).
#   6. [job]    SoupX autoEstCont + adjustCounts (R via fhR module).
#   7. [job]    Reimport SoupX output → adata_gex_spliced_soupx.h5ad.
#   8. [job]    Preprocess post-SoupX → adata_gex_spliced_soupx_processed.h5ad (final).
#
# Usage:
#   bash CapGTA_gex_soupx_pipeline.sh \
#       <real_cells.csv> \
#       <annotation.gtf> \
#       <output_dir> \
#       [--n-empties 5000] \
#       [--resolution 2]
#
# Terminal output:
#   <output_dir>/adata_gex_spliced_soupx_processed.h5ad
###########################################################################################################################

set -euo pipefail

N_EMPTIES=5000
RESOLUTION=2

POS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --n-empties)   N_EMPTIES="$2";   shift 2 ;;
        --resolution)  RESOLUTION="$2";  shift 2 ;;
        --)            shift; POS+=("$@"); break ;;
        *)             POS+=("$1");      shift ;;
    esac
done

if [ "${#POS[@]}" -ne 3 ]; then
    echo "Usage: $0 <real_cells.csv> <annotation.gtf> <output_dir> [--n-empties N] [--resolution R]" >&2
    exit 1
fi

REAL_CELLS_CSV="${POS[0]}"
GTF="${POS[1]}"
OUTPUT_DIR="${POS[2]}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${REPO_ROOT}/scripts/CapGTA"

for f in "${REAL_CELLS_CSV}" "${GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}" SLURM_outs

CELLS_RNA_DIR="${OUTPUT_DIR}/rna_counts"
EMPTIES_DIR="${OUTPUT_DIR}/soupx_empties"
EMPTIES_RNA_DIR="${EMPTIES_DIR}/rna_counts"
SOUPX_IN="${OUTPUT_DIR}/soupx_in"
SOUPX_OUT="${OUTPUT_DIR}/soupx_out"

CELLS_RAW_H5AD="${CELLS_RNA_DIR}/rna_counts.h5ad"
CELLS_PROCESSED_H5AD="${OUTPUT_DIR}/adata_gex_spliced.h5ad"
EMPTIES_H5AD="${EMPTIES_RNA_DIR}/rna_counts.h5ad"
JUNCTIONS_CSV="${OUTPUT_DIR}/gene_junctions.csv"
SOUPX_H5AD="${OUTPUT_DIR}/adata_gex_spliced_soupx.h5ad"
FINAL_H5AD="${OUTPUT_DIR}/adata_gex_spliced_soupx_processed.h5ad"
EMPTIES_CSV="${EMPTIES_DIR}/empty_cells.csv"

echo "=========================================="
echo "GEX + SoupX pipeline"
echo "=========================================="
echo "Real cells:     ${REAL_CELLS_CSV}"
echo "GTF:            ${GTF}"
echo "Output dir:     ${OUTPUT_DIR}"
echo "N empties:      ${N_EMPTIES}"
echo "Leiden res:     ${RESOLUTION}"
echo "Final h5ad:     ${FINAL_H5AD}"
echo "=========================================="

# --- Stage 0: junction table (inline, fast) --------------------------------
if [ ! -f "${JUNCTIONS_CSV}" ]; then
    echo "[Stage 0] Building junction table…"
    python3 "${SCRIPTS_DIR}/gtf_gene_junctions.py" "${GTF}" "${JUNCTIONS_CSV}"
else
    echo "[Stage 0] Junction table exists, skipping: ${JUNCTIONS_CSV}"
fi

# --- Stage 1: real-cell spliced counts (returns terminal h5ad job id) ------
echo ""
echo "[Stage 1] Submitting real-cell spliced counts…"
mkdir -p "${CELLS_RNA_DIR}"
STAGE1_OUT=$("${SCRIPTS_DIR}/submit_rna_counts.sh" \
    "${REAL_CELLS_CSV}" "${GTF}" "${CELLS_RNA_DIR}")
echo "${STAGE1_OUT}"
J1=$(echo "${STAGE1_OUT}" | grep -oE "Submitted h5ad conversion:\s+[0-9]+" | grep -oE "[0-9]+$")
if [ -z "${J1}" ]; then
    echo "Error: could not parse h5ad job id from Stage 1 output" >&2
    exit 1
fi

# --- Stage 2: preprocess real cells ----------------------------------------
echo ""
echo "[Stage 2] Submitting preprocess (real cells)…"
J2=$(sbatch --parsable --dependency=afterok:${J1} \
    -J gex_preprocess_cells -o SLURM_outs/%x_%j.out \
    -p campus-new -c 4 --mem=32G -t 1:00:00 \
    --wrap="python3 '${SCRIPTS_DIR}/preprocess_scanpy.py' \
        '${CELLS_RAW_H5AD}' '${CELLS_PROCESSED_H5AD}' \
        --gtf '${GTF}' \
        --junction-csv '${JUNCTIONS_CSV}' --resolution ${RESOLUTION}")
echo "  job ${J2} (depends on ${J1})"

# --- Stage 3: enumerate empties (inline, fast) -----------------------------
echo ""
echo "[Stage 3] Enumerating empty droplets…"
mkdir -p "${EMPTIES_DIR}"
python3 "${SCRIPTS_DIR}/enumerate_empties.py" \
    "${REAL_CELLS_CSV}" "${EMPTIES_CSV}" --n-empties "${N_EMPTIES}"

# --- Stage 4: empty-droplet spliced counts (returns terminal h5ad job id) --
echo ""
echo "[Stage 4] Submitting empty-droplet spliced counts…"
mkdir -p "${EMPTIES_RNA_DIR}"
STAGE4_OUT=$("${SCRIPTS_DIR}/submit_rna_counts.sh" \
    "${EMPTIES_CSV}" "${GTF}" "${EMPTIES_RNA_DIR}")
echo "${STAGE4_OUT}"
J4=$(echo "${STAGE4_OUT}" | grep -oE "Submitted h5ad conversion:\s+[0-9]+" | grep -oE "[0-9]+$")
if [ -z "${J4}" ]; then
    echo "Error: could not parse h5ad job id from Stage 4 output" >&2
    exit 1
fi

# --- Stage 5: prep SoupX inputs (depends on J2 + J4) -----------------------
echo ""
echo "[Stage 5] Submitting prep_soupx_inputs…"
J5=$(sbatch --parsable --dependency=afterok:${J2}:${J4} \
    -J gex_prep_soupx -o SLURM_outs/%x_%j.out \
    -p campus-new -c 4 --mem=32G -t 30:00 \
    --wrap="python3 '${SCRIPTS_DIR}/prep_soupx_inputs.py' \
        '${CELLS_PROCESSED_H5AD}' '${EMPTIES_H5AD}' '${SOUPX_IN}'")
echo "  job ${J5} (depends on ${J2} + ${J4})"

# --- Stage 6: SoupX (R) ----------------------------------------------------
echo ""
echo "[Stage 6] Submitting SoupX…"
J6=$(sbatch --parsable --dependency=afterok:${J5} \
    -J gex_soupx -o SLURM_outs/%x_%j.out \
    -p campus-new -c 4 --mem=32G -t 1:00:00 \
    --wrap="source /app/lmod/lmod/init/bash && ml fhR && \
        Rscript '${SCRIPTS_DIR}/soupx_decontaminate.R' '${SOUPX_IN}' '${SOUPX_OUT}'")
echo "  job ${J6} (depends on ${J5})"

# --- Stage 7: reimport SoupX → h5ad ----------------------------------------
echo ""
echo "[Stage 7] Submitting reimport_soupx…"
J7=$(sbatch --parsable --dependency=afterok:${J6} \
    -J gex_reimport_soupx -o SLURM_outs/%x_%j.out \
    -p campus-new -c 2 --mem=16G -t 30:00 \
    --wrap="python3 '${SCRIPTS_DIR}/reimport_soupx.py' \
        '${CELLS_PROCESSED_H5AD}' '${SOUPX_OUT}' '${SOUPX_H5AD}'")
echo "  job ${J7} (depends on ${J6})"

# --- Stage 8: preprocess post-SoupX ----------------------------------------
echo ""
echo "[Stage 8] Submitting preprocess (post-SoupX)…"
J8=$(sbatch --parsable --dependency=afterok:${J7} \
    -J gex_preprocess_soupx -o SLURM_outs/%x_%j.out \
    -p campus-new -c 4 --mem=32G -t 1:00:00 \
    --wrap="python3 '${SCRIPTS_DIR}/preprocess_scanpy.py' \
        '${SOUPX_H5AD}' '${FINAL_H5AD}' \
        --gtf '${GTF}' --drop-cuticle \
        --junction-csv '${JUNCTIONS_CSV}' --resolution ${RESOLUTION}")
echo "  job ${J8} (depends on ${J7})"

echo ""
echo "=========================================="
echo "Pipeline submitted. Job chain:"
echo "  ${J1} → ${J2} ─┐"
echo "                 ├→ ${J5} → ${J6} → ${J7} → ${J8}"
echo "  ${J4}         ─┘"
echo ""
echo "Monitor:   squeue -u \$USER"
echo "Terminal:  ${FINAL_H5AD}"
echo "=========================================="
