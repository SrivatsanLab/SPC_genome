#!/bin/bash
#SBATCH --job-name=bulk_alignment
#SBATCH --output=SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH -t 4:00:00

######################################################################################################
# Bulk sample alignment pipeline (no demultiplexing, no variant calling)
# Performs: Trimming -> Alignment -> BAM sorting/indexing -> QC metrics
######################################################################################################

mkdir -p SLURM_outs/

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

# Function to display help message
show_help() {
    cat << EOF
Usage: $0 -o <OUTPUT_NAME> -1 <read1.fastq.gz> -2 <read2.fastq.gz> -g <reference_genome>

This script aligns bulk samples without demultiplexing or variant calling.

Required arguments:
  -o    <output_name>           Desired sample name (prefix for outputs)
  -1    <read1.fastq.gz>        Read 1 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R1*.fastq.gz")
  -2    <read2.fastq.gz>        Read 2 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R2*.fastq.gz")
  -g    <reference_genome>      Path to BWA index (e.g., /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa)

Optional arguments:
  -O    <output_dir>            Output directory (default: ./data/bulk_samples/aligned)
  -s    <scripts_dir>           Scripts directory (default: .)
  -h                            Show this help message and exit

EOF
}

# Default values
SCRIPTS_DIR="."
OUTPUT_DIR="./data/bulk_samples/aligned"

# Parse command-line options
while getopts ":o:1:2:g:O:s:h" option; do
  case $option in
    o) OUTPUT_NAME=$OPTARG ;;
    1) READ1=$OPTARG ;;
    2) READ2=$OPTARG ;;
    g) REFERENCE_GENOME=$OPTARG ;;
    O) OUTPUT_DIR=$OPTARG ;;
    s) SCRIPTS_DIR=$OPTARG ;;
    h) show_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; show_help; exit 1 ;;
  esac
done

# Check that required arguments are provided
if [ -z "$OUTPUT_NAME" ] || [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$REFERENCE_GENOME" ]; then
    echo "Missing required arguments." >&2
    show_help
    exit 1
fi

# Create output directories
mkdir -p "${OUTPUT_DIR}"
BIN_DIR="${SCRIPTS_DIR}/bin/${OUTPUT_NAME}"
mkdir -p "${BIN_DIR}"

# Create temporary directory for intermediate files
TMP_DIR="/hpc/temp/srivatsan_s/bulk_alignment_${OUTPUT_NAME}_$$"
mkdir -p "${TMP_DIR}"

# Print inputs
echo "=========================================="
echo "Bulk Alignment Pipeline"
echo "=========================================="
echo "Sample Name: ${OUTPUT_NAME}"
echo "Read 1: ${READ1}"
echo "Read 2: ${READ2}"
echo "Reference Genome: ${REFERENCE_GENOME}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Bin Directory: ${BIN_DIR}"
echo "Temp Directory: ${TMP_DIR}"
echo "=========================================="

######################################################################################################
#### Trimming
######################################################################################################

echo ""
echo "Step 1: Trimming adapters with Trim Galore..."

module load cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

# For multi-lane support: concatenate all R1 and R2 files, then trim
# Use process substitution to avoid creating large intermediate concatenated files
CONCAT_R1="${TMP_DIR}/concatenated_R1.fastq.gz"
CONCAT_R2="${TMP_DIR}/concatenated_R2.fastq.gz"

echo "Concatenating input files..."
# Note: READ1 and READ2 are intentionally unquoted to allow shell expansion for multi-lane patterns
cat $READ1 > "${CONCAT_R1}"
cat $READ2 > "${CONCAT_R2}"

trim_galore --paired --cores 4 -o "${TMP_DIR}" --illumina --gzip "${CONCAT_R1}" "${CONCAT_R2}"

# Delete concatenated files to save space
rm "${CONCAT_R1}" "${CONCAT_R2}"

module unload cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

# Get trimmed file names
TRIMMED_R1="${TMP_DIR}/concatenated_R1_val_1.fq.gz"
TRIMMED_R2="${TMP_DIR}/concatenated_R2_val_2.fq.gz"

echo "Trimmed files created:"
echo "  R1: ${TRIMMED_R1}"
echo "  R2: ${TRIMMED_R2}"

######################################################################################################
#### Alignment
######################################################################################################

echo ""
echo "Step 2: Aligning with BWA-MEM..."

module load BWA SAMtools

SAM_FILE="${TMP_DIR}/${OUTPUT_NAME}.sam"

echo "Running BWA-MEM alignment..."
bwa mem -t 8 -o "${SAM_FILE}" "${REFERENCE_GENOME}" "${TRIMMED_R1}" "${TRIMMED_R2}"

# Delete trimmed fastqs to save space
rm "${TRIMMED_R1}" "${TRIMMED_R2}"

######################################################################################################
#### Convert to BAM, sort, and index
######################################################################################################

echo ""
echo "Step 3: Converting to BAM, sorting, and indexing..."

BAM_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.bam"

samtools view -@ 8 -bS "${SAM_FILE}" | samtools sort -@ 8 -o "${BAM_FILE}"
samtools index -@ 8 "${BAM_FILE}"

# Delete SAM file
rm "${SAM_FILE}"

echo "BAM file created and indexed: ${BAM_FILE}"

######################################################################################################
#### Generate QC metrics
######################################################################################################

echo ""
echo "Step 4: Generating QC metrics..."

# Basic alignment statistics
samtools flagstat "${BAM_FILE}" > "${BIN_DIR}/${OUTPUT_NAME}_flagstat.txt"
samtools stats "${BAM_FILE}" > "${BIN_DIR}/${OUTPUT_NAME}_stats.txt"

echo "QC metrics saved to:"
echo "  Flagstat: ${BIN_DIR}/${OUTPUT_NAME}_flagstat.txt"
echo "  Stats: ${BIN_DIR}/${OUTPUT_NAME}_stats.txt"

######################################################################################################
#### Cleanup
######################################################################################################

echo ""
echo "Step 5: Cleaning up temporary files..."

rm -rf "${TMP_DIR}"

echo ""
echo "=========================================="
echo "Bulk alignment pipeline completed!"
echo "=========================================="
echo "Output files:"
echo "  BAM: ${BAM_FILE}"
echo "  BAM index: ${BAM_FILE}.bai"
echo "  Flagstat: ${BIN_DIR}/${OUTPUT_NAME}_flagstat.txt"
echo "  Stats: ${BIN_DIR}/${OUTPUT_NAME}_stats.txt"
echo "=========================================="
