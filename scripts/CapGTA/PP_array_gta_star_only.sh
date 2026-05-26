#!/bin/bash
#SBATCH -J PP_array_gta_star
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 8
#SBATCH -p short
#SBATCH -t 8:00:00
#SBATCH --mem=32G

set -euo pipefail

##########################################################################################################################
# STAR-only alignment for genome-transcriptome coassay data:                                                           #
# 1. Demultiplex and trim reads (FASTQ format)                                                                         #
# 2. Convert trimmed FASTQ to unaligned SAM with CB:Z tags (preserves barcodes as proper SAM tags)                    #
# 3. STAR aligns from SAM input (preserves CB:Z tags automatically)                                                    #
# 4. Split aligned reads by splice junctions:                                                                          #
#    - RNA BAM: Reads with splice junctions (N in CIGAR) = high-confidence transcripts                                #
#    - DNA BAM: All other reads (includes unspliced RNA, but provides conservative RNA calls)                         #
##########################################################################################################################

# Activate conda environment for Python scripts (atrandi_demux.py)
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

mkdir -p SLURM_outs/array_outs

chunk_indices="$1"
genome="$2"
scripts_DIR="$3"
TMP_DIR="$4"
DEBUG="${5:-false}"

barcodes="${scripts_DIR}/barcodes"
demux_scr="${scripts_DIR}/scripts/utils/atrandi_demux.py"

chunk=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$chunk_indices")

READ1="${TMP_DIR}/read1_chunk_${chunk}"
READ2="${TMP_DIR}/read2_chunk_${chunk}"

# In debug mode, keep intermediates instead of deleting them.
debug_rm() { [ "$DEBUG" = "true" ] && return 0; rm "$@"; }

##########################################################################################################################
# Demultiplexing - extract barcode from read2, add to headers, delete reads lacking a legitimate barcode
##########################################################################################################################

python $demux_scr $READ1 $READ2 $barcodes \
  --R1_output "${TMP_DIR}/corr_read1_chunk_${chunk}" \
  --R2_output "${TMP_DIR}/corr_read2_chunk_${chunk}" \
  --gzip False \
  --ignore-overhangs

# Delete uncorrected fastqs
debug_rm $READ1 $READ2

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
debug_rm $READ1 $READ2

##########################################################################################################################
# Convert trimmed FASTQ to unaligned SAM with CB:Z tags
##########################################################################################################################

# Reset read1 and read2 to trimmed version
READ1="${TMP_DIR}/corr_read1_chunk_${chunk}_val_1.fq.gz"
READ2="${TMP_DIR}/corr_read2_chunk_${chunk}_val_2.fq.gz"

UNALIGNED_SAM="${TMP_DIR}/${chunk}_unaligned.sam"

echo "Converting trimmed FASTQ to unaligned SAM with CB tags..."
python "${scripts_DIR}/scripts/utils/fastq_to_unaligned_sam.py" "${READ1}" "${READ2}" "${UNALIGNED_SAM}"

# Delete trimmed fastqs to save space
debug_rm "${READ1}" "${READ2}"

##########################################################################################################################
# STAR Alignment (all reads) - using SAM input to preserve CB tags
##########################################################################################################################

# Load STAR 2.7.9a
module load STAR/2.7.9a-GCC-11.2.0

# Construct STAR index path from genome parameter
if [ -d "${genome}/STARIndex" ]; then
    STAR_INDEX="${genome}/STARIndex"
    echo "Using STAR index: ${STAR_INDEX}"
else
    echo "Error: Could not find STAR index at ${genome}/STARIndex"
    exit 1
fi

echo "Aligning all reads with STAR (SAM input with CB tags)..."

# In debug mode, retain unmapped reads inline so they can be routed to a 3rd SAM stream below.
STAR_UNMAPPED_OPT=""
if [ "$DEBUG" = "true" ]; then
    STAR_UNMAPPED_OPT="--outSAMunmapped Within"
fi

# STAR alignment - using unaligned SAM as input preserves CB tags automatically
# Don't include CB in --outSAMattributes (causes error), STAR preserves all input tags by default
STAR --runThreadN 8 \
     --genomeDir "${STAR_INDEX}" \
     --readFilesIn "${UNALIGNED_SAM}" \
     --readFilesType SAM PE \
     --outFileNamePrefix "${TMP_DIR}/${chunk}_" \
     --outSAMtype BAM Unsorted \
     --outSAMattributes NH HI AS nM NM MD \
     --outFilterMultimapNmax 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 15000 \
     --alignMatesGapMax 15000 \
     ${STAR_UNMAPPED_OPT}

STAR_BAM="${TMP_DIR}/${chunk}_Aligned.out.bam"

# Delete unaligned SAM to save space
debug_rm "${UNALIGNED_SAM}"

##########################################################################################################################
# CB tags are now proper SAM tags - STAR preserves them automatically
##########################################################################################################################

module load SAMtools

##########################################################################################################################
# Split reads by splice junctions
# RNA: Reads with N in CIGAR (splice junctions)
# DNA: All other reads
##########################################################################################################################

DNA_SAM="${TMP_DIR}/${chunk}_dna.sam"
RNA_SAM="${TMP_DIR}/${chunk}_rna.sam"
UNMAPPED_SAM="${TMP_DIR}/${chunk}_unmapped.sam"

echo "Splitting reads by splice junctions..."

# Extract header
samtools view -H "${STAR_BAM}" > "${DNA_SAM}"
samtools view -H "${STAR_BAM}" > "${RNA_SAM}"
if [ "$DEBUG" = "true" ]; then
    samtools view -H "${STAR_BAM}" > "${UNMAPPED_SAM}"
fi

# Split reads:
#   - unmapped (FLAG & 4) -> unmapped (debug mode only; in normal mode STAR drops them)
#   - mapped + N in CIGAR -> RNA
#   - else                -> DNA
# CB tags are already preserved in the BAM by STAR (from unaligned SAM input).
samtools view "${STAR_BAM}" | awk -v debug="${DEBUG}" '
BEGIN {
    dna_sam = "'"${DNA_SAM}"'"
    rna_sam = "'"${RNA_SAM}"'"
    unmapped_sam = "'"${UNMAPPED_SAM}"'"
}
{
    flag = $2
    cigar = $6

    # Bit 4 (0x4) of FLAG = "unmapped". Use int arithmetic (portable across awk implementations).
    is_unmapped = (int(flag / 4) % 2 == 1)

    if (is_unmapped) {
        if (debug == "true") {
            print $0 >> unmapped_sam
        }
        next
    }

    if (cigar ~ /N/) {
        print $0 >> rna_sam
    } else {
        print $0 >> dna_sam
    }
}
'

# Clean up
debug_rm "${STAR_BAM}"

# Clean up STAR temporary files
if [ "$DEBUG" != "true" ]; then
    rm -f "${TMP_DIR}/${chunk}_Log.out"
    rm -f "${TMP_DIR}/${chunk}_Log.progress.out"
    rm -f "${TMP_DIR}/${chunk}_Log.final.out"
    rm -f "${TMP_DIR}/${chunk}_SJ.out.tab"
    rm -rf "${TMP_DIR}/${chunk}__STARtmp"
fi

module unload STAR SAMtools

echo "Alignment complete!"
echo "DNA SAM: ${DNA_SAM}"
echo "RNA SAM: ${RNA_SAM}"
if [ "$DEBUG" = "true" ]; then
    echo "Unmapped SAM: ${UNMAPPED_SAM}"
fi
