#!/bin/bash
#SBATCH -J feature_dist
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 2:00:00

###########################################################################################################################
# Quantify the proportion of reads aligning to exons, introns, and intergenic regions
# Compares DNA vs RNA BAMs to verify RNA enrichment for exonic reads
###########################################################################################################################

set -euo pipefail

BAM_DNA="$1"
BAM_RNA="$2"
GTF="$3"
OUTPUT_DIR="$4"

mkdir -p "${OUTPUT_DIR}"

echo "Analyzing genomic feature distribution..."
echo "DNA BAM: ${BAM_DNA}"
echo "RNA BAM: ${BAM_RNA}"
echo "GTF: ${GTF}"
echo "Output: ${OUTPUT_DIR}"
echo ""

module load BEDTools SAMtools

# Create BED files for genomic features
EXONS_BED="${OUTPUT_DIR}/exons.bed"
GENES_BED="${OUTPUT_DIR}/genes.bed"
INTRONS_BED="${OUTPUT_DIR}/introns.bed"

echo "Extracting genomic features from GTF..."

# Extract exon coordinates
awk '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' "${GTF}" | \
    sed 's/"//g' | sed 's/;//g' | \
    sort -k1,1 -k2,2n > "${EXONS_BED}"

# Extract gene coordinates by grouping transcripts by gene_id
# This GTF doesn't have "gene" features, so we group transcripts
awk '$3 == "transcript" {
    # Extract gene_id
    match($0, /gene_id "([^"]+)"/, arr)
    gene_id = arr[1]
    chr = $1
    start = $4 - 1
    end = $5
    strand = $7

    # Track min/max coords for each gene
    if (gene_id in gene_chr) {
        if (start < gene_start[gene_id]) gene_start[gene_id] = start
        if (end > gene_end[gene_id]) gene_end[gene_id] = end
    } else {
        gene_chr[gene_id] = chr
        gene_start[gene_id] = start
        gene_end[gene_id] = end
        gene_strand[gene_id] = strand
    }
}
END {
    for (g in gene_chr) {
        print gene_chr[g]"\t"gene_start[g]"\t"gene_end[g]"\t"g"\t.\t"gene_strand[g]
    }
}' "${GTF}" | sort -k1,1 -k2,2n > "${GENES_BED}"

# Create intron regions (gene minus exons)
bedtools subtract -a "${GENES_BED}" -b "${EXONS_BED}" > "${INTRONS_BED}"

echo "Genomic features extracted:"
echo "  Exons: $(wc -l < ${EXONS_BED})"
echo "  Genes: $(wc -l < ${GENES_BED})"
echo "  Introns: $(wc -l < ${INTRONS_BED})"
echo ""

# Function to analyze a BAM file
analyze_bam() {
    local bam="$1"
    local prefix="$2"

    echo "Analyzing ${prefix}..."

    # Convert BAM to BED (using proper paired-end handling)
    local reads_bed="${OUTPUT_DIR}/${prefix}_reads.bed"
    samtools view -b -F 4 "${bam}" | \
        bedtools bamtobed -i stdin | \
        sort -k1,1 -k2,2n > "${reads_bed}"

    local total_reads=$(wc -l < "${reads_bed}")

    # Count reads overlapping exons
    local exon_reads=$(bedtools intersect -a "${reads_bed}" -b "${EXONS_BED}" -u | wc -l)

    # Count reads overlapping introns
    local intron_reads=$(bedtools intersect -a "${reads_bed}" -b "${INTRONS_BED}" -u | wc -l)

    # Intergenic = total - (exonic + intronic)
    # Note: some reads may overlap both exon and intron, so we need to be careful
    # Better approach: get reads in genes, then subtract to get intergenic
    local genic_reads=$(bedtools intersect -a "${reads_bed}" -b "${GENES_BED}" -u | wc -l)
    local intergenic_reads=$((total_reads - genic_reads))

    # Calculate percentages
    local exon_pct=$(awk "BEGIN {printf \"%.2f\", (${exon_reads}/${total_reads})*100}")
    local intron_pct=$(awk "BEGIN {printf \"%.2f\", (${intron_reads}/${total_reads})*100}")
    local intergenic_pct=$(awk "BEGIN {printf \"%.2f\", (${intergenic_reads}/${total_reads})*100}")

    echo "Results for ${prefix}:"
    echo "  Total aligned reads: ${total_reads}"
    echo "  Exonic reads: ${exon_reads} (${exon_pct}%)"
    echo "  Intronic reads: ${intron_reads} (${intron_pct}%)"
    echo "  Intergenic reads: ${intergenic_reads} (${intergenic_pct}%)"
    echo ""

    # Save to file
    cat > "${OUTPUT_DIR}/${prefix}_feature_distribution.txt" <<EOF
Feature Distribution for ${prefix}
================================
Total aligned reads: ${total_reads}

Exonic reads:       ${exon_reads} (${exon_pct}%)
Intronic reads:     ${intron_reads} (${intron_pct}%)
Intergenic reads:   ${intergenic_reads} (${intergenic_pct}%)
EOF

    # Clean up temporary BED
    rm "${reads_bed}"
}

# Analyze both BAMs
analyze_bam "${BAM_DNA}" "dna"
analyze_bam "${BAM_RNA}" "rna"

echo ""
echo "=== Summary ==="
echo "DNA BAM:"
cat "${OUTPUT_DIR}/dna_feature_distribution.txt" | grep -E "Exonic|Intronic|Intergenic"
echo ""
echo "RNA BAM:"
cat "${OUTPUT_DIR}/rna_feature_distribution.txt" | grep -E "Exonic|Intronic|Intergenic"

# Create a comparison CSV
cat > "${OUTPUT_DIR}/feature_distribution_comparison.csv" <<EOF
Feature,DNA_Count,DNA_Percent,RNA_Count,RNA_Percent
EOF

# Extract values and add to CSV
dna_exon=$(grep "Exonic reads:" "${OUTPUT_DIR}/dna_feature_distribution.txt" | awk '{print $3}')
dna_exon_pct=$(grep "Exonic reads:" "${OUTPUT_DIR}/dna_feature_distribution.txt" | awk '{print $4}' | tr -d '(%)' )
rna_exon=$(grep "Exonic reads:" "${OUTPUT_DIR}/rna_feature_distribution.txt" | awk '{print $3}')
rna_exon_pct=$(grep "Exonic reads:" "${OUTPUT_DIR}/rna_feature_distribution.txt" | awk '{print $4}' | tr -d '(%)' )

dna_intron=$(grep "Intronic reads:" "${OUTPUT_DIR}/dna_feature_distribution.txt" | awk '{print $3}')
dna_intron_pct=$(grep "Intronic reads:" "${OUTPUT_DIR}/dna_feature_distribution.txt" | awk '{print $4}' | tr -d '(%)' )
rna_intron=$(grep "Intronic reads:" "${OUTPUT_DIR}/rna_feature_distribution.txt" | awk '{print $3}')
rna_intron_pct=$(grep "Intronic reads:" "${OUTPUT_DIR}/rna_feature_distribution.txt" | awk '{print $4}' | tr -d '(%)' )

dna_inter=$(grep "Intergenic reads:" "${OUTPUT_DIR}/dna_feature_distribution.txt" | awk '{print $3}')
dna_inter_pct=$(grep "Intergenic reads:" "${OUTPUT_DIR}/dna_feature_distribution.txt" | awk '{print $4}' | tr -d '(%)' )
rna_inter=$(grep "Intergenic reads:" "${OUTPUT_DIR}/rna_feature_distribution.txt" | awk '{print $3}')
rna_inter_pct=$(grep "Intronic reads:" "${OUTPUT_DIR}/rna_feature_distribution.txt" | awk '{print $4}' | tr -d '(%)' )

cat >> "${OUTPUT_DIR}/feature_distribution_comparison.csv" <<EOF
Exonic,${dna_exon},${dna_exon_pct},${rna_exon},${rna_exon_pct}
Intronic,${dna_intron},${dna_intron_pct},${rna_intron},${rna_intron_pct}
Intergenic,${dna_inter},${dna_inter_pct},${rna_inter},${rna_inter_pct}
EOF

echo ""
echo "Analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}"

module unload BEDTools SAMtools
