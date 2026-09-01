#!/bin/bash
#SBATCH -J regenotype
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# Wrap for regenotype_from_ad.py. `#!/bin/bash` (not the /bin/sh that sbatch --wrap uses) so that
# `module load` — a bash-only function from /etc/profile.d/modules.sh — is defined.
#
# Usage:
#   sbatch scripts/CapGTA/pseudobulk_gatk/regenotype_job.sh <output_dir> <mm_env>
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 <output_dir> <mm_env>" >&2
    exit 1
fi

OUT="$1"
MMENV="$2"
# sbatch copies this file into a per-job tmpdir before running, so BASH_SOURCE[0] is unreliable.
# SLURM_SUBMIT_DIR = directory sbatch was launched from = repo root in this pipeline.
PY="${SLURM_SUBMIT_DIR}/scripts/CapGTA/pseudobulk_gatk/regenotype_from_ad.py"

module load GATK

micromamba run -n "${MMENV}" python3 "${PY}" "${OUT}"
