#!/bin/bash
#SBATCH -J ad
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4

set -euo pipefail

module load BCFtools

# Arguments
temp="$1"
map="$2"
intervals="$3"
genome="$4"
scripts_dir="$5"

# Convert BED-style interval to GATK format (1-based, inclusive)
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${intervals}" | awk '{print $1 ":" $2+1 "-" $3}')

DB_PATH="${temp}/database_${SLURM_ARRAY_TASK_ID}"

GATK_output="${temp}/${INTERVAL}.vcf.gz"

norm_temp="${temp}/${INTERVAL}_norm_temp.bcf"

final_output="${temp}/${INTERVAL}.bcf"
anndata="${temp}/${INTERVAL}.h5ad"

# Run Mutyper
#bcftools norm -m - "${GATK_output}" -Ob -o "${norm_temp}" 
#bcftools norm -m - "${GATK_output}" -Ob -o "${final_output}"

#mutyper variants "${genome}" "${norm_temp}" | bcftools view --threads 4 -Ob -o "${final_output}"

# Make Anndata Object
python "${scripts_dir}/make_variant_anndata.py" "${final_output}" "${anndata}"

# cleanup
rm "${norm_temp}"
#rm "${GATK_output}"
