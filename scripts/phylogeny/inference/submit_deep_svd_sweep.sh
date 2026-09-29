#!/bin/bash
# Deep SVD → imputation → CellPhy sweep on tct_kept_relax03.
#
# Three variants explore whether pushing the SVD bipartition tree past its
# current min-clade-size=10 floor yields more useful imputation. For each
# variant we chain:
#
#   svd (16-task array)  →  impute (16)  →  prepare_vcf (16)  →  cellphy (48)
#
# Only the "real" matrix runs in CellPhy for this sweep (skipping the
# covmask/shuffled controls to cut CellPhy compute 3x). Add MATRICES="real
# covmask shuffled" back if you want the controls.
#
# Variants:
#   A: mquad5_cs5    — MQUAD_BIC=5,  min_clade_size=5   (current BIC, deeper)
#   B: mquad10_cs5   — MQUAD_BIC=10, min_clade_size=5   (stricter markers)
#   C: mquad10_cs3   — MQUAD_BIC=10, min_clade_size=3   (deepest, stricter markers)
#
# Usage:
#   scripts/phylogeny/inference/submit_deep_svd_sweep.sh [SRC_TAG]
#
# SRC_TAG defaults to tct_kept_relax03. Pass all_variants_relax03 (or any
# other panel tag under panels/) to run the sweep on a different source panel.
#
# Prints one line per (variant, stage): VARIANT/STAGE=JOBID

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome

SRC_TAG="${1:-tct_kept_relax03}"

# variant_name | MQUAD_BIC | MIN_CLADE_SIZE | MAX_DEPTH
VARIANTS=(
    "mquad5_cs5   5.0   5   10"
    "mquad10_cs5  10.0  5   10"
    "mquad10_cs3  10.0  3   12"
)

for row in "${VARIANTS[@]}"; do
    read -r CORE_TAG MQUAD_BIC MIN_CLADE MAX_DEPTH <<< "$row"
    IMPUTED_TAG="${SRC_TAG}__svdImp_${CORE_TAG}"

    echo ""
    echo "=== ${CORE_TAG}: BIC=${MQUAD_BIC}  min_clade=${MIN_CLADE}  max_depth=${MAX_DEPTH} ==="
    echo "    svd_subdir=svd_${CORE_TAG}"
    echo "    imputed_tag=${IMPUTED_TAG}"

    # 1. SVD tree
    SVD=$(sbatch --job-name="svd_${CORE_TAG}" \
        --export=ALL,PANEL_TAG="${SRC_TAG}",CORE_SEL="mquad",CORE_TAG="${CORE_TAG}",MQUAD_BIC="${MQUAD_BIC}",MIN_CLADE_SIZE="${MIN_CLADE}",MAX_DEPTH="${MAX_DEPTH}" \
        --parsable scripts/phylogeny/inference/sbatch_svd.sh)
    echo "${CORE_TAG}/svd=${SVD}"

    # 2. Impute panel using the new tree.
    IMP=$(sbatch --job-name="imp_${CORE_TAG}" \
        --export=ALL,SRC_TAG="${SRC_TAG}",IMPUTED_TAG="${IMPUTED_TAG}",TREE_SUBDIR="svd_${CORE_TAG}" \
        --dependency=afterok:${SVD} \
        --parsable scripts/phylogeny/inference/sbatch_impute.sh)
    echo "${CORE_TAG}/impute=${IMP}"

    # 3. Prepare VCFs from the imputed panel.
    PREP=$(sbatch --job-name="prep_${CORE_TAG}" \
        --export=ALL,PANEL_TAG="${IMPUTED_TAG}" \
        --dependency=afterok:${IMP} \
        --parsable scripts/phylogeny/inference/sbatch_prepare_imputed.sh)
    echo "${CORE_TAG}/prepare=${PREP}"

    # 4. CellPhy ML — real matrix only for the sweep.
    CP=$(sbatch --job-name="cp_${CORE_TAG}" \
        --export=ALL,PANEL_TAG="${IMPUTED_TAG}",MATRICES="real" \
        --dependency=afterok:${PREP} \
        --parsable scripts/phylogeny/inference/sbatch_cellphy_imputed.sh)
    echo "${CORE_TAG}/cellphy=${CP}"
done

echo ""
echo "All three variants submitted. Monitor with: squeue -u \$USER --format='%.10i %.20j %.8T %.10M'"
