#!/bin/bash
#SBATCH -J jc
#SBATCH --output=SLURM_outs/%x_%j.out
#SBATCH -c 36

module load GATK

database="$1"
map="$2"
intervals="$3"
genome="$4"
output="$5"

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport \
  --genomicsdb-workspace-path "${database}" \
  --intervals "${intervals}" \
  --sample-name-map "${map}" \
  --reader-threads 4 \
  --max-num-intervals-to-import-in-parallel 9 \
  --batch-size 100 
  
gatk GenotypeGVCFs \
  -R "${genome}" \
  -V "gendb://${database}" \
  -O "${output}"