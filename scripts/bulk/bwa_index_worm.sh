#!/bin/bash
#SBATCH --job-name=bwa_index_worm
#SBATCH --output=SLURM_outs/bwa_index_worm_%j.out
#SBATCH -c 4
#SBATCH -t 2:00:00

set -euo pipefail

# Load BWA module
module load BWA/0.7.18-GCCcore-13.3.0

# Reference genome path
REF_GENOME="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/reference/worm_GCA_028201515.1/GCA_028201515.1_genomic.fna"

echo "Creating BWA index for ${REF_GENOME}"
echo "Start time: $(date)"

# Create BWA index
bwa index "${REF_GENOME}"

echo "BWA index creation complete"
echo "End time: $(date)"

# List index files
echo "Index files created:"
ls -lh "${REF_GENOME}"*
