#!/bin/bash
#SBATCH -J concat_sc_bams
#SBATCH --output=SLURM_outs/%x_%j.out
#SBATCH -c 16

module load SAMtools

samtools cat -b bam_list.txt -o sc_PolE_bulked.bam
