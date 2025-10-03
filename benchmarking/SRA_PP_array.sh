#!/bin/bash
#SBATCH -J PP_array
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH -p short

mkdir -p SLURM_outs/array_outs

SRRs="$1"
path="$2"
genome="$3"
TMP_DIR="$4"
OUT_DIR="$5"

name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRRs")

##########################################################################################################################
# fastq dump 
##########################################################################################################################

fasterq-dump "${path}/${name}" -O "${TMP_DIR}" -t "${TMP_DIR}" -m 16G

READ1="${TMP_DIR}/${name}_1.fastq"
READ2="${TMP_DIR}/${name}_2.fastq"

##########################################################################################################################
# Trimming 
##########################################################################################################################

module load cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

trim_galore --paired --cores 4 -o "${TMP_DIR}" --gzip "${READ1}" "${READ2}" 

module unload cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

# Delete untrimmed fastqs
rm $READ1 $READ2

##########################################################################################################################
# Alignment 
##########################################################################################################################

# Reset read1 and read2 to trimmed version
READ1="${TMP_DIR}/${name}_1_val_1.fq.gz"
READ2="${TMP_DIR}/${name}_2_val_2.fq.gz"

module load BWA SAMtools

# Output file name
SAM_FILE="${TMP_DIR}/${name}.sam"
BAM_FILE="${OUT_DIR}/${name}.bam"

echo "Aligning with BWA-MEM"
bwa mem -t 4 -c 1 -C -o "${SAM_FILE}" "${genome}" "${READ1}" "${READ2}"

samtools view "${SAM_FILE}" -b | samtools sort -o "${BAM_FILE}" && samtools index "${BAM_FILE}"

# Delete temp files:
#rm "${READ1}" "${READ2}" "${SAM_FILE}"

##########################################################################################################################
# Coverage 
##########################################################################################################################

module unload BWA SAMtools
module load deepTools

BIGWIG="${OUT_DIR}/${name}.bw"

bamCoverage -b "${BAM_FILE}" -o "${BIGWIG}" --binSize 1000 -of bigwig -p 4
