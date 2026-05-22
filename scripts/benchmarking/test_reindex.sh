#!/bin/bash
#SBATCH -J test_diag
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 30

set -uo pipefail

module load GATK
module load BCFtools/1.18-GCC-12.2.0

F=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/benchmarking/gvcfs/B_lymphocyte/SRR8438251.g.vcf.gz
REF=/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta
WORK=/hpc/temp/srivatsan_s/htsjdk_diag_$SLURM_JOB_ID
mkdir -p "$WORK"

echo "=== Tool versions ==="
which bgzip bcftools tabix
which gatk
echo

echo "=== Test 4: Bcftools header read on ORIGINAL (control) ==="
bcftools view -h "$F" 2>&1 | head -3
echo "..."
bcftools view -h "$F" 2>&1 | wc -l
echo "Bcftools view -h exit: $?"
echo

# Build small header-only test files in TMP
HEADER_FULL="$WORK/header_full.vcf"
zcat "$F" | head -287 > "$HEADER_FULL"
ls -la "$HEADER_FULL"
echo "Lines in header_full: $(wc -l < $HEADER_FULL)"
echo "Longest line in header: $(awk '{print length}' $HEADER_FULL | sort -n | tail -1)"
echo

echo "=== Test 3a: Header + 1 record, BGZF-compressed; index; ValidateVariants ==="
HDR_BGZ="$WORK/header_full.vcf.gz"
bgzip -c "$HEADER_FULL" > "$HDR_BGZ"
tabix -p vcf "$HDR_BGZ"
ls -la "$HDR_BGZ" "$HDR_BGZ.tbi"
gatk ValidateVariants -V "$HDR_BGZ" -gvcf -R "$REF" 2>&1 | tail -5
echo "Header+1rec validate exit: $?"
echo

echo "=== Test 2: Strip ##GATKCommandLine line, header-only repack ==="
HDR_STRIPPED="$WORK/header_stripped.vcf"
grep -v '^##GATKCommandLine' "$HEADER_FULL" > "$HDR_STRIPPED"
echo "Lines in stripped header: $(wc -l < $HDR_STRIPPED)"
echo "Longest line in stripped header: $(awk '{print length}' $HDR_STRIPPED | sort -n | tail -1)"
HDR_STRIP_BGZ="$WORK/header_stripped.vcf.gz"
bgzip -c "$HDR_STRIPPED" > "$HDR_STRIP_BGZ"
tabix -p vcf "$HDR_STRIP_BGZ"
gatk ValidateVariants -V "$HDR_STRIP_BGZ" -gvcf -R "$REF" 2>&1 | tail -5
echo "Stripped validate exit: $?"
echo

echo "=== Test 5: ValidateVariants on ORIGINAL (control) ==="
gatk ValidateVariants -V "$F" -gvcf -R "$REF" 2>&1 | tail -5
echo "Original validate exit: $?"

echo "=== Cleanup ==="
ls -la "$WORK"
rm -rf "$WORK"
