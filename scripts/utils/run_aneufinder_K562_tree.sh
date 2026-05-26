#!/bin/bash
#SBATCH --job-name=aneufinder_K562_tree
#SBATCH --output=SLURM_outs/aneufinder_K562_tree_%j.out
#SBATCH -c 36
#SBATCH --mem=180G
#SBATCH -t 48:00:00

set -euo pipefail

module load fhR/4.4.1-foss-2023b

cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome

Rscript scripts/utils/run_aneufinder_K562_tree.R
