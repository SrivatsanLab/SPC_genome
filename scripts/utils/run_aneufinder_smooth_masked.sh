#!/bin/bash
#SBATCH --job-name=aneufinder_smooth_masked
#SBATCH --output=SLURM_outs/aneufinder_smooth_masked_%j.out
#SBATCH -c 36
#SBATCH --mem=120G
#SBATCH -t 12:00:00

# Repeat-masked, smoothed (500 kb windows, 250 kb step) re-segmentation of the
# 1000 sc_PolE K562 cells from the GC-corrected 250 kb binned counts; see
# scripts/utils/aneufinder_smooth_masked.R. Run run_aneufinder_K562_sc_PolE_gc.sh 2.5e5 first.

set -euo pipefail
module load fhR/4.4.1-foss-2023b
export R_LIBS_USER=/home/sanjay/R/aneufinder-bioc3.19
cd /home/sanjay/SPC_genome
Rscript scripts/utils/aneufinder_smooth_masked.R
