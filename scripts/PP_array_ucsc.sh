#!/bin/bash
#SBATCH -J PP_array
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH -p short
#SBATCH -t 8:00:00

set -euo pipefail

##########################################################################################################################
# This job performs custom de-multiplexing and standard bioinformatics (trimming, alignment) on chunks of paired fastq's #
# Modified version for UCSC hg38 reference genome                                                                        #
##########################################################################################################################

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

mkdir -p SLURM_outs/array_outs

chunk_indices="$1"
genome="$2"
scripts_DIR="$3"
TMP_DIR="$4"

barcodes="${scripts_DIR}/barcodes"
demux_scr="${scripts_DIR}/scripts/atrandi_demux.py"

chunk=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$chunk_indices")

READ1="${TMP_DIR}/read1_chunk_${chunk}"
READ2="${TMP_DIR}/read2_chunk_${chunk}"

##########################################################################################################################
# Demultiplexing- extract barcode from read2, add to headers, delete reads lacking a legitimate barcode
##########################################################################################################################

python $demux_scr $READ1 $READ2 $barcodes \
  --R1_output "${TMP_DIR}/corr_read1_chunk_${chunk}" \
  --R2_output "${TMP_DIR}/corr_read2_chunk_${chunk}" \
  --gzip False

#delete uncorrected fastqs
rm $READ1 $READ2

##########################################################################################################################
# Trimming
##########################################################################################################################

# Reset read1 and read2 to corrected version
READ1="${TMP_DIR}/corr_read1_chunk_${chunk}"
READ2="${TMP_DIR}/corr_read2_chunk_${chunk}"

module load cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

trim_galore --paired --cores 4 -o "${TMP_DIR}" --illumina --gzip "${READ1}" "${READ2}"

module unload cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

# Delete untrimmed fastqs
rm $READ1 $READ2

##########################################################################################################################
# Alignment
##########################################################################################################################

# Reset read1 and read2 to trimmed version
READ1="${TMP_DIR}/corr_read1_chunk_${chunk}_val_1.fq.gz"
READ2="${TMP_DIR}/corr_read2_chunk_${chunk}_val_2.fq.gz"

module load BWA

# Construct BWA index path - handle both UCSC and GATK reference structures
if [ -f "${genome}" ]; then
    # If genome is a direct path to fasta file
    BWA_INDEX="${genome}"
elif [ -f "${genome}/BWAIndex/genome.fa" ]; then
    # UCSC iGenomes structure
    BWA_INDEX="${genome}/BWAIndex/genome.fa"
elif [ -f "${genome}/BWAIndex/Homo_sapiens_assembly38.fasta.64" ]; then
    # GATK bundle structure
    BWA_INDEX="${genome}/BWAIndex/Homo_sapiens_assembly38.fasta.64"
else
    echo "Error: Could not find BWA index at expected locations"
    echo "Tried:"
    echo "  ${genome}"
    echo "  ${genome}/BWAIndex/genome.fa"
    echo "  ${genome}/BWAIndex/Homo_sapiens_assembly38.fasta.64"
    exit 1
fi

# Output file name
SAM_FILE="${TMP_DIR}/${chunk}.sam"

echo "Aligning with BWA-MEM using index: ${BWA_INDEX}"
bwa mem -t 4 -c 1 -C -o "${SAM_FILE}" "${BWA_INDEX}" "${READ1}" "${READ2}"

# Delete fastq's
rm "${READ1}" "${READ2}"
