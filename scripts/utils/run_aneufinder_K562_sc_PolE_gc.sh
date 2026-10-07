#!/bin/bash
#SBATCH --job-name=aneufinder_K562_sc_PolE_gc
#SBATCH --output=SLURM_outs/aneufinder_K562_sc_PolE_gc_%j.out
#SBATCH -c 36
#SBATCH --mem=180G
#SBATCH -t 48:00:00

# Re-run AneuFinder on the 1000 sc_PolE_novaseq K562 cells behind the per-cell
# SBS tree (sc_og_test.newick), this time with GC correction and the ENCODE
# blacklist; the original run (sc_PolE_novaseq/AneuFinder_output) had neither.
# Input is a folder of symlinks to those cells' BAMs; see
# paper_figures/scripts/K562_cnv_event_tree.R for the tree built from the output.
#
# Optional argument: bin size in bp (default 1e6). Other bin sizes write to
# aneufinder_gc_<binsize>/output and skip the per-cell profile and karyogram
# plots, e.g.
#   sbatch scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh 5e5

set -euo pipefail

module load fhR/4.4.1-foss-2023b
# AneuFinder lives in its own library: the default user library holds older
# Bioconductor 3.18 GenomicRanges/GenomeInfoDb that mask the module's 3.19
# builds and break binning ("object 'normarg_seqnames1' not found")
export R_LIBS_USER=/home/sanjay/R/aneufinder-bioc3.19

cd /home/sanjay/SPC_genome

BINSIZE=${1:-1e6}
run_dir=/home/sanjay/SPC_genome/results/K562_tree/aneufinder_gc
if [[ "$BINSIZE" != "1e6" ]]; then
    run_dir=${run_dir}_${BINSIZE}
    export ANEUFINDER_SKIP_SLOW_PLOTS=1
fi

export ANEUFINDER_BINSIZE=$BINSIZE
export ANEUFINDER_INPUT=/home/sanjay/SPC_genome/results/K562_tree/aneufinder_gc/input
export ANEUFINDER_OUTPUT=${run_dir}/output

Rscript scripts/utils/run_aneufinder_K562_tree.R
