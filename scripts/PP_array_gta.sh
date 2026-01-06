#!/bin/bash
#SBATCH -J PP_array_gta
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 8
#SBATCH -p short
#SBATCH -t 8:00:00
#SBATCH --mem=32G

set -euo pipefail

##########################################################################################################################
# This job performs dual-alignment for genome-transcriptome coassay data:                                               #
# 1. BWA-MEM alignment for DNA reads                                                                                    #
# 2. STAR alignment for RNA reads (splice-aware)                                                                        #
# Reads are classified based on BWA alignment quality (MAPQ and soft-clipping)                                          #
##########################################################################################################################

# Activate conda environment for Python scripts (atrandi_demux.py)
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
# Demultiplexing - extract barcode from read2, add to headers, delete reads lacking a legitimate barcode
##########################################################################################################################

python $demux_scr $READ1 $READ2 $barcodes \
  --R1_output "${TMP_DIR}/corr_read1_chunk_${chunk}" \
  --R2_output "${TMP_DIR}/corr_read2_chunk_${chunk}" \
  --gzip False

# Delete uncorrected fastqs
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
# PASS 1: BWA-MEM Alignment
##########################################################################################################################

# Reset read1 and read2 to trimmed version
READ1="${TMP_DIR}/corr_read1_chunk_${chunk}_val_1.fq.gz"
READ2="${TMP_DIR}/corr_read2_chunk_${chunk}_val_2.fq.gz"

module load BWA SAMtools

# Construct BWA index path for C. elegans
if [ -f "${genome}/BWAIndex/genome.fa" ]; then
    BWA_INDEX="${genome}/BWAIndex/genome.fa"
elif [ -f "${genome}/genome.fa" ]; then
    BWA_INDEX="${genome}/genome.fa"
else
    echo "Error: Could not find BWA index at expected locations"
    echo "Tried:"
    echo "  ${genome}/BWAIndex/genome.fa"
    echo "  ${genome}/genome.fa"
    exit 1
fi

# Temporary SAM file for BWA alignment
BWA_SAM="${TMP_DIR}/${chunk}_bwa_all.sam"

echo "Aligning with BWA-MEM using index: ${BWA_INDEX}"
bwa mem -t 8 -c 1 -C -o "${BWA_SAM}" "${BWA_INDEX}" "${READ1}" "${READ2}"

##########################################################################################################################
# Classify reads: DNA vs RNA candidates
# DNA: MAPQ >= 20 AND soft-clipping < 10bp
# RNA candidates: Unmapped OR MAPQ < 20 OR soft-clipping >= 10bp
##########################################################################################################################

DNA_SAM="${TMP_DIR}/${chunk}_dna.sam"
RNA_QNAMES="${TMP_DIR}/${chunk}_rna_qnames.txt"

echo "Classifying reads into DNA and RNA candidates..."

# Extract header
samtools view -H "${BWA_SAM}" > "${DNA_SAM}"

# Process alignments
# For DNA: keep reads with MAPQ >= 20 and soft-clipping < 10bp
# For RNA: extract read names for everything else
samtools view "${BWA_SAM}" | awk '
BEGIN {
    dna_sam = "'"${DNA_SAM}"'"
    rna_qnames = "'"${RNA_QNAMES}"'"
}
{
    qname = $1
    flag = $3
    mapq = $5
    cigar = $6

    # Parse CIGAR string to count soft-clipped bases
    soft_clip = 0
    len = length(cigar)
    num = ""
    for (i = 1; i <= len; i++) {
        c = substr(cigar, i, 1)
        if (c ~ /[0-9]/) {
            num = num c
        } else {
            if (c == "S" && num != "") {
                soft_clip += num
            }
            num = ""
        }
    }

    # Classify read
    if (mapq >= 20 && soft_clip < 10) {
        # DNA read
        print $0 >> dna_sam
    } else {
        # RNA candidate - save read name
        print qname >> rna_qnames
    }
}
'

echo "DNA reads extracted to: ${DNA_SAM}"
echo "RNA candidate read names extracted to: ${RNA_QNAMES}"

# Remove full BWA SAM to save space
rm "${BWA_SAM}"

module unload BWA SAMtools

##########################################################################################################################
# Extract RNA candidate reads as FASTQ for STAR alignment
##########################################################################################################################

echo "Extracting RNA candidate reads to FASTQ..."

# Sort RNA read names for faster lookup
sort -u "${RNA_QNAMES}" > "${RNA_QNAMES}.sorted"
mv "${RNA_QNAMES}.sorted" "${RNA_QNAMES}"

RNA_R1="${TMP_DIR}/${chunk}_rna_R1.fastq"
RNA_R2="${TMP_DIR}/${chunk}_rna_R2.fastq"
CB_MAP="${TMP_DIR}/${chunk}_rna_cb_map.txt"

# Extract reads using seqtk (faster than custom script)
module load seqtk

seqtk subseq "${READ1}" "${RNA_QNAMES}" > "${RNA_R1}"
seqtk subseq "${READ2}" "${RNA_QNAMES}" > "${RNA_R2}"

module unload seqtk

# Create a lookup table of read names to CB tags from the FASTQ file
echo "Creating CB tag lookup table from RNA FASTQ..."
awk 'NR % 4 == 1 {
    # Parse FASTQ header: @READNAME TAGS
    split($0, header, " ")
    read_name = substr(header[1], 2)  # Remove @ prefix

    # Extract CB tag from remaining fields
    for (i = 2; i <= NF; i++) {
        if ($i ~ /^CB:Z:/) {
            cb_tag = $i
            print read_name "\t" cb_tag
            break
        }
    }
}' "${RNA_R1}" > "${CB_MAP}"

echo "RNA FASTQ files created: ${RNA_R1}, ${RNA_R2}"
echo "CB lookup table created: ${CB_MAP}"

# Delete original trimmed fastqs to save space
rm "${READ1}" "${READ2}"

##########################################################################################################################
# PASS 2: STAR Alignment for RNA candidates
##########################################################################################################################

RNA_SAM="${TMP_DIR}/${chunk}_rna.sam"

# Check if there are any RNA candidates
rna_count=$(wc -l < "${RNA_QNAMES}")

if [ "$rna_count" -eq 0 ]; then
    echo "No RNA candidate reads found. Creating empty RNA SAM file."
    # Create empty SAM with just header
    module load STAR
    STAR_INDEX="${genome}/STARIndex"

    # Get genome file for header
    GENOME_FA="${genome}/WholeGenomeFasta/genome.fa"

    module load SAMtools
    samtools view -H -T "${GENOME_FA}" > "${RNA_SAM}"
    module unload SAMtools STAR
else
    echo "Aligning ${rna_count} RNA candidate read pairs with STAR..."

    # Load STAR 2.7.1a - compatible with versionGenome 20201 format
    module load STAR/2.7.1a-foss-2019b

    # Construct STAR index path
    if [ -d "${genome}/STARIndex" ]; then
        STAR_INDEX="${genome}/STARIndex"
    else
        echo "Error: Could not find STAR index at: ${genome}/STARIndex"
        exit 1
    fi

    # STAR alignment
    # Output directly as SAM
    STAR --runThreadN 8 \
         --genomeDir "${STAR_INDEX}" \
         --readFilesIn "${RNA_R1}" "${RNA_R2}" \
         --outFileNamePrefix "${TMP_DIR}/${chunk}_rna_" \
         --outSAMtype SAM \
         --outSAMattributes NH HI AS nM NM MD \
         --outFilterMultimapNmax 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 15000 \
         --alignMatesGapMax 15000

    STAR_SAM="${TMP_DIR}/${chunk}_rna_Aligned.out.sam"

    # Add CB tags to STAR output using the lookup table
    echo "Adding CB tags to STAR alignments..."
    module load SAMtools

    # Extract header
    samtools view -H "${STAR_SAM}" > "${RNA_SAM}"

    # Add CB tags to alignments by matching read names
    samtools view "${STAR_SAM}" | awk '
    BEGIN {
        # Load CB tag lookup table
        cb_map_file = "'"${CB_MAP}"'"
        while ((getline < cb_map_file) > 0) {
            cb_tags[$1] = $2
        }
        close(cb_map_file)
        rna_sam = "'"${RNA_SAM}"'"
    }
    {
        read_name = $1
        if (read_name in cb_tags) {
            # Append CB tag to the line
            print $0 "\t" cb_tags[read_name] >> rna_sam
        } else {
            # No CB tag found (shouldn'"'"'t happen, but handle gracefully)
            print $0 >> rna_sam
        }
    }
    '

    module unload SAMtools

    # Clean up STAR temporary files
    rm -f "${STAR_SAM}"
    rm -f "${TMP_DIR}/${chunk}_rna_Log.out"
    rm -f "${TMP_DIR}/${chunk}_rna_Log.progress.out"
    rm -f "${TMP_DIR}/${chunk}_rna_Log.final.out"
    rm -f "${TMP_DIR}/${chunk}_rna_SJ.out.tab"
    rm -rf "${TMP_DIR}/${chunk}_rna__STARtmp"

    module unload STAR
fi

# Clean up intermediate files
rm "${RNA_R1}" "${RNA_R2}" "${RNA_QNAMES}"

if [ -f "${CB_MAP}" ]; then
    rm "${CB_MAP}"
fi

echo "Alignment complete!"
echo "DNA SAM: ${DNA_SAM}"
echo "RNA SAM: ${RNA_SAM}"
