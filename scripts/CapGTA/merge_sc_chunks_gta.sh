#!/bin/bash
#SBATCH -J merge_sc_gta
#SBATCH -o SLURM_outs/%x_%A.out
#SBATCH -c 4
#SBATCH -t 4:00:00
#SBATCH --mem=16G

##########################################################################################################################
# Merge per-chunk single-cell BAMs into final per-cell DNA and RNA BAMs
##########################################################################################################################

barcode_file="$1"
sc_output_dir="$2"
final_output_dir="$3"

module load SAMtools

mkdir -p "${final_output_dir}"

# Count barcodes
cell_count=$(wc -l < "$barcode_file")

echo "Merging chunked BAMs for ${cell_count} cells..."
echo "Source: ${sc_output_dir}"
echo "Destination: ${final_output_dir}"
echo ""

# Process each barcode
while IFS= read -r barcode; do
    echo "Processing barcode: ${barcode}"

    # Find all DNA chunk BAMs for this barcode
    dna_chunks=("${sc_output_dir}/${barcode}"_*_dna.bam)
    rna_chunks=("${sc_output_dir}/${barcode}"_*_rna.bam)

    # Output files
    final_dna="${final_output_dir}/${barcode}_dna.bam"
    final_rna="${final_output_dir}/${barcode}_rna.bam"

    # Merge DNA chunks
    if [ ${#dna_chunks[@]} -gt 0 ] && [ -f "${dna_chunks[0]}" ]; then
        echo "  Merging ${#dna_chunks[@]} DNA chunks..."
        samtools merge -@ 4 -f "${final_dna}" "${dna_chunks[@]}"
        samtools index "${final_dna}"

        # Clean up chunk files
        rm -f "${dna_chunks[@]}"
    else
        echo "  Warning: No DNA chunks found for ${barcode}"
    fi

    # Merge RNA chunks
    if [ ${#rna_chunks[@]} -gt 0 ] && [ -f "${rna_chunks[0]}" ]; then
        echo "  Merging ${#rna_chunks[@]} RNA chunks..."
        samtools merge -@ 4 -f "${final_rna}" "${rna_chunks[@]}"
        samtools index "${final_rna}"

        # Clean up chunk files
        rm -f "${rna_chunks[@]}"
    else
        echo "  Warning: No RNA chunks found for ${barcode}"
    fi

    echo "  Complete: ${barcode}"
    echo ""
done < "$barcode_file"

echo "====================================="
echo "Single-cell BAM merging complete!"
echo "====================================="
echo "Output directory: ${final_output_dir}"
echo "Files per cell: 2 (DNA + RNA BAMs)"
