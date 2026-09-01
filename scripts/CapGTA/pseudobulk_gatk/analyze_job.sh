#!/bin/bash
#SBATCH -J analyze_pb
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 2:00:00

###########################################################################################################################
# Wrap for analyze_worms_vs_subclones.py — same reason as regenotype_job.sh (needs bash for `module`).
#
# Usage:
#   sbatch scripts/CapGTA/pseudobulk_gatk/analyze_job.sh <output_dir> <existing_vcf> <mm_env>
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch $0 <output_dir> <existing_vcf> <mm_env>" >&2
    exit 1
fi

OUT="$1"
EXISTING_VCF="$2"
MMENV="$3"
# sbatch copies this file into a per-job tmpdir before running, so BASH_SOURCE[0] is unreliable.
# SLURM_SUBMIT_DIR = directory sbatch was launched from = repo root in this pipeline.
PY="${SLURM_SUBMIT_DIR}/scripts/CapGTA/pseudobulk_gatk/analyze_worms_vs_subclones.py"

module load BCFtools

micromamba run -n "${MMENV}" python3 "${PY}" "${OUT}" "${EXISTING_VCF}"
