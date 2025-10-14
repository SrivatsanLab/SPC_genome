#!/bin/bash
#SBATCH -J jc
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4

set -euo pipefail

# Activate conda environment for mutyper and Python scripts
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

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
  --batch-size 100 \
  --overwrite-existing-genomicsdb-workspace true

gatk GenotypeGVCFs \
  -R "${genome}" \
  -V "gendb://${DB_PATH}" \
  -O "${GATK_output}"

module unload GATK

# Load BCFtools for normalization and filtering (use version compatible with GATK's GCC)
module load BCFtools/1.18-GCC-12.2.0

# Normalize VCF (split multi-allelic sites)
bcftools norm -m - "${GATK_output}" -Ob -o "${final_output}"

# Unload modules to avoid Python path conflicts, then run Python script
module purge

# Make Anndata Object from normalized VCF
$CONDA_PREFIX/bin/python "${scripts_dir}/scripts/make_variant_anndata.py" "${final_output}" "${anndata}"

# cleanup
rm "${GATK_output}"