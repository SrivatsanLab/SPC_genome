#!/bin/bash
#SBATCH -J compile_h5ad
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 4
#SBATCH --mem=32G

###########################################################################################################################
# Compile interval h5ad files into final AnnData object
###########################################################################################################################

set -euo pipefail

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

TMP_DIR="$1"
OUTPUT_H5AD="$2"

# Check if h5ad files exist
h5ad_count=$(find "${TMP_DIR}" -name "*.h5ad" | wc -l)

if [ "$h5ad_count" -eq 0 ]; then
    echo "ERROR: No h5ad files found in ${TMP_DIR}"
    exit 1
fi

echo "Found ${h5ad_count} h5ad files to compile"

# Use SLURM_SUBMIT_DIR to find the repository root
# (this is where sbatch was called from)
REPO_DIR="${SLURM_SUBMIT_DIR}"

# Compile all h5ad files
python "${REPO_DIR}/scripts/compile_adatas.py" "${TMP_DIR}" "${OUTPUT_H5AD}"

echo "Compiled h5ad saved to: ${OUTPUT_H5AD}"

# Cleanup temp files
echo "Cleaning up intermediate files..."
rm -rf "${TMP_DIR}"

echo "Done!"
