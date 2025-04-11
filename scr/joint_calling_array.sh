#!/bin/bash
#SBATCH -J jc
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4

set -euo pipefail

# Arguments
temp="$1"
map="$2"
intervals="$3"
genome="$4"

# Prepare workspace
mkdir -p "${temp}"

# Convert BED-style interval to GATK format (1-based, inclusive)
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${intervals}" | awk '{print $1 ":" $2+1 "-" $3}')
DB_PATH="${temp}/database_${SLURM_ARRAY_TASK_ID}"
GATK_output="${temp}/${INTERVAL}.vcf.gz"
VEP_output="${temp}/${INTERVAL}_anno.vcf.gz"

echo "Running for interval: ${INTERVAL}" 
echo "GenomicsDB path: ${DB_PATH}"
echo "VCF output: ${GATK_output}"

# Load GATK and run GenomicsDBImport + GenotypeGVCFs
module load GATK 

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport \
  --genomicsdb-workspace-path "${DB_PATH}" \
  -L "${INTERVAL}" \
  --sample-name-map "${map}" \
  --reader-threads 4 \
  --batch-size 100

gatk GenotypeGVCFs \
  -R "${genome}" \
  -V "gendb://${DB_PATH}" \
  -O "${GATK_output}"

module unload GATK

# Load VEP and annotate
module load VEP

vep \
  -i "${GATK_output}" \
  -o "${VEP_output}" \
  --vcf \
  --compress_output bgzip \
  --cache \
  --offline \
  --dir_cache ~/.vep \
  --assembly GRCh38 \
  --species homo_sapiens \
  --everything \
  --fork "${SLURM_CPUS_PER_TASK}" \
  --verbose

# Overwrite unannotated VCF with annotated version
mv "${VEP_output}" "${GATK_output}"

# Index the annotated VCF
tabix -p vcf "${GATK_output}"