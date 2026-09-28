#!/bin/bash
#SBATCH --job-name=aneufinder_gc_toggle
#SBATCH --output=SLURM_outs/aneufinder_gc_toggle_%j.out
#SBATCH -c 36
#SBATCH --mem=64G
#SBATCH -t 2:00:00

# Copy-number calls with and without GC correction on the same saved bins;
# see scripts/utils/aneufinder_gc_toggle.R.

set -euo pipefail

module load fhR/4.4.1-foss-2023b
# see run_aneufinder_K562_sc_PolE_gc.sh for why AneuFinder has its own library
export R_LIBS_USER=/home/sanjay/R/aneufinder-bioc3.19

cd /home/sanjay/SPC_genome

export ANEUFINDER_OUTPUT=/home/sanjay/SPC_genome/results/K562_tree/aneufinder_gc/output

Rscript scripts/utils/aneufinder_gc_toggle.R
