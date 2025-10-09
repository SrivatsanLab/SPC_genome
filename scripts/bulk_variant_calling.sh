#!/bin/bash
#SBATCH -J bulk_calling
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4

bam="$1"
intervals="$2"
genome="$3"

temp="temp_bulk_vcf"
mkdir $temp

GATK_output="${temp}/${INTERVAL}.vcf.gz"

norm_temp="${temp}/${INTERVAL}_norm_temp.bcf"

final_output="${temp}/${INTERVAL}.bcf"