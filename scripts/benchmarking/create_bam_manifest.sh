#!/bin/bash
##########################################################################################################################
# Manifest Builder for Downsampling Analysis
#
# Purpose: Process user-provided BAM list and extract metadata
# Input:   Text file with one BAM path per line
# Output:  TSV with columns: dataset, cell_id, bam_path, reference, total_mapped_reads
##########################################################################################################################

set -euo pipefail

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_bam_list.txt> <output_manifest.tsv>"
    echo ""
    echo "Example:"
    echo "  $0 data/benchmarking/manifests/input_bams.txt data/benchmarking/manifests/bam_manifest.tsv"
    exit 1
fi

INPUT_LIST="$1"
OUTPUT_TSV="$2"

# Validate input file exists
if [ ! -f "$INPUT_LIST" ]; then
    echo "Error: Input file not found: $INPUT_LIST"
    exit 1
fi

# Load samtools
module load SAMtools

echo "==========================================="
echo "BAM Manifest Builder"
echo "==========================================="
echo "Input list: $INPUT_LIST"
echo "Output TSV: $OUTPUT_TSV"
echo ""

# Create output directory if needed
mkdir -p "$(dirname "$OUTPUT_TSV")"

# Write header
echo -e "dataset\tcell_id\tbam_path\treference\ttotal_mapped_reads" > "$OUTPUT_TSV"

# Counter
TOTAL_BAMS=0
SUCCESS_COUNT=0
FAIL_COUNT=0

# Process each BAM
while IFS= read -r BAM_PATH || [ -n "$BAM_PATH" ]; do
    # Skip empty lines and comments
    [[ -z "$BAM_PATH" || "$BAM_PATH" =~ ^[[:space:]]*# ]] && continue

    TOTAL_BAMS=$((TOTAL_BAMS + 1))

    echo "Processing BAM $TOTAL_BAMS: $BAM_PATH"

    # Validate BAM exists
    if [ ! -f "$BAM_PATH" ]; then
        echo "  WARNING: BAM file not found, skipping"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    # Extract dataset (parent directory name)
    PARENT_DIR=$(dirname "$BAM_PATH")
    DATASET=$(basename "$PARENT_DIR")

    # If parent is 'sc_outputs', go up one more level
    if [ "$DATASET" = "sc_outputs" ]; then
        DATASET=$(basename "$(dirname "$PARENT_DIR")")
    fi

    # Extract cell_id (basename without .bam)
    CELL_ID=$(basename "$BAM_PATH" .bam)

    # Extract reference from BAM header
    # Look for the first @SQ line and try to reconstruct the reference path
    # We'll extract the sequence name and try to find the original FASTA
    echo "  Extracting reference from BAM header..."

    # Get full BAM header
    HEADER=$(samtools view -H "$BAM_PATH" 2>/dev/null || echo "")

    if [ -z "$HEADER" ]; then
        echo "  WARNING: Could not read BAM header, skipping"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    # Try to find @PG lines with reference information
    # Often the aligner command line contains the reference path
    REFERENCE=$(echo "$HEADER" | grep "^@PG" | grep -oP '(?<=-x\s|--reference\s|-R\s)\S+' | head -1 || echo "")

    # If that didn't work, try to extract from @SQ lines and make an educated guess
    if [ -z "$REFERENCE" ]; then
        # Get first sequence name
        FIRST_SEQ=$(echo "$HEADER" | grep "^@SQ" | head -1 | grep -oP 'SN:\K\S+' || echo "")

        # Common reference patterns based on sequence names
        if [[ "$FIRST_SEQ" =~ ^chr ]]; then
            # Likely human, try common paths
            if [ -f "/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta" ]; then
                REFERENCE="/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta"
            elif [ -f "/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa" ]; then
                REFERENCE="/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
            else
                REFERENCE="UNKNOWN_hg38"
            fi
        else
            # Unknown reference
            REFERENCE="UNKNOWN_${FIRST_SEQ}"
        fi
    fi

    echo "  Reference: $REFERENCE"

    # Count primary mapped reads
    echo "  Counting mapped reads..."
    MAPPED_READS=$(samtools view -c -F 260 "$BAM_PATH" 2>/dev/null || echo "0")

    if [ "$MAPPED_READS" -eq 0 ]; then
        echo "  WARNING: No mapped reads found, skipping"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    echo "  Mapped reads: $MAPPED_READS"

    # Write to output
    echo -e "${DATASET}\t${CELL_ID}\t${BAM_PATH}\t${REFERENCE}\t${MAPPED_READS}" >> "$OUTPUT_TSV"

    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    echo "  ✓ Success"
    echo ""

done < "$INPUT_LIST"

echo "==========================================="
echo "Manifest building complete!"
echo "==========================================="
echo "Total BAMs processed: $TOTAL_BAMS"
echo "Successfully added: $SUCCESS_COUNT"
echo "Failed/skipped: $FAIL_COUNT"
echo ""
echo "Output written to: $OUTPUT_TSV"
echo ""

# Show first few lines of output
echo "Preview of manifest:"
head -5 "$OUTPUT_TSV" | column -t -s $'\t'
echo ""
