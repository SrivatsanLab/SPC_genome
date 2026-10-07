#!/bin/bash
#SBATCH --job-name=fragment_gc_segment
#SBATCH --output=SLURM_outs/fragment_gc_segment_%j.out
#SBATCH -c 36
#SBATCH --mem=120G
#SBATCH -t 6:00:00

# Mask, mappability-correct, smooth and segment the fragment-level GC-corrected
# depth (scripts/utils/fragment_gc_segment.R). Run after run_fragment_gc_counts.sh.

set -euo pipefail
module load fhR/4.4.1-foss-2023b
export R_LIBS_USER=/home/sanjay/R/aneufinder-bioc3.19
cd /home/sanjay/SPC_genome
Rscript scripts/utils/fragment_gc_segment.R
