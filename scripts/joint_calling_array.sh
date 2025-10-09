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
scripts_dir="$5"

# Prepare workspace
mkdir -p "${temp}"

# Convert BED-style interval to GATK format (1-based, inclusive)
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${intervals}" | awk '{print $1 ":" $2+1 "-" $3}')
DB_PATH="${temp}/database_${SLURM_ARRAY_TASK_ID}"

GATK_output="${temp}/${INTERVAL}.vcf.gz"

norm_temp="${temp}/${INTERVAL}_norm_temp.bcf"

final_output="${temp}/${INTERVAL}.bcf"
anndata="${temp}/${INTERVAL}.h5ad"

echo "Running for interval: ${INTERVAL}" 

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

# Run Mutyper
bcftools norm -m - "${GATK_output}" -Ob -o "${norm_temp}" 

mutyper variants "${genome}" "${norm_temp}" | bcftools view --threads 8 -Ob -o "${final_output}"

# Make Anndata Object
python "${scripts_dir}/make_variant_anndata.py" "${final_output}" "${anndata}"

# cleanup
rm "${norm_temp}"
rm "${GATK_output}"