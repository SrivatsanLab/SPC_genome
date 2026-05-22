#!/bin/bash
#SBATCH -J test_reheader
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 30

set -uo pipefail

module load GATK
module load BCFtools/1.18-GCC-12.2.0

F=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/benchmarking/gvcfs/B_lymphocyte/SRR8438251.g.vcf.gz
REF=/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta
WORK=/hpc/temp/srivatsan_s/reheader_test_$SLURM_JOB_ID
mkdir -p "$WORK"

NEW_HDR="$WORK/new_header.txt"
FIXED="$WORK/SRR8438251.fixed.g.vcf.gz"

echo "=== Step 1: extract header, strip ##GATKCommandLine ==="
bcftools view -h "$F" | grep -v '^##GATKCommandLine' > "$NEW_HDR"
echo "Lines in new header: $(wc -l < $NEW_HDR)"
echo "Longest line: $(awk '{print length}' $NEW_HDR | sort -n | tail -1)"
echo

echo "=== Step 2: bcftools reheader (does NOT decompress data) ==="
time bcftools reheader -h "$NEW_HDR" -o "$FIXED" "$F"
ls -la "$FIXED"
echo

echo "=== Step 3: tabix the fixed file ==="
time tabix -p vcf "$FIXED"
ls -la "$FIXED"*
echo

echo "=== Step 4: htslib still reads it ==="
bcftools view -h "$FIXED" | wc -l
echo

echo "=== Step 5: GATK ValidateVariants on fixed file ==="
time gatk ValidateVariants -V "$FIXED" -gvcf -R "$REF" 2>&1 | tail -5
echo "Validate exit: $?"
echo

echo "=== Cleanup ==="
ls -la "$WORK"
rm -rf "$WORK"
