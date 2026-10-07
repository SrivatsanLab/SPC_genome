#!/bin/bash
#SBATCH --job-name=fragment_gc
#SBATCH --output=SLURM_outs/array_outs/fragment_gc_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH --array=1-100

# Fragment-level GC-corrected 250 kb depth for the 1000 sc_PolE K562 cells,
# 10 cells per task (scripts/utils/fragment_gc_counts.R). Run
# scripts/utils/fragment_gc_setup.R first; cell list from
# ls results/K562_tree/aneufinder_gc/input/*.bam > results/K562_tree/fragment_gc/cell_bams.txt

set -euo pipefail
module load fhR/4.4.1-foss-2023b
export R_LIBS_USER=/home/sanjay/R/aneufinder-bioc3.19
cd /home/sanjay/SPC_genome
Rscript scripts/utils/fragment_gc_counts.R results/K562_tree/fragment_gc/cell_bams.txt "$SLURM_ARRAY_TASK_ID" 10
