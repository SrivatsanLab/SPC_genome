#!/bin/bash
#SBATCH --job-name=rt_regress_segment
#SBATCH --output=SLURM_outs/rt_regress_segment_%j.out
#SBATCH -c 36
#SBATCH --mem=120G
#SBATCH -t 6:00:00

# Fit and remove each cell's replication (cell-cycle) component, then re-call
# copy number (scripts/utils/replication_regress_segment.R). Needs
# run_fragment_gc_segment.sh and paper_figures/scripts/K562_replication_timing_coverage.R first.

set -euo pipefail
module load fhR/4.4.1-foss-2023b
export R_LIBS_USER=/home/sanjay/R/aneufinder-bioc3.19
cd /home/sanjay/SPC_genome
Rscript scripts/utils/replication_regress_segment.R
