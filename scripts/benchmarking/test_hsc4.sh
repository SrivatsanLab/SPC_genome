#!/bin/bash
#SBATCH -J test_hsc4
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 30

set -uo pipefail

module load GATK
module load BCFtools/1.18-GCC-12.2.0

REF=/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta

echo "=== Testing 3 HSC4 GVCFs with htsjdk via gatk IndexFeatureFile (header read only) ==="
COUNT=0
for f in $(ls /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/benchmarking/gvcfs/HSC4/*.g.vcf.gz | head -3); do
    echo "--- File: $f"
    echo "Resolved: $(readlink -f $f)"
    echo "Longest header line: $(bcftools view -h $f 2>/dev/null | awk '{print length}' | sort -n | tail -1)"
    # Try IndexFeatureFile (writes nothing if --output isn't passed; reads header)
    echo "Test: gatk ValidateVariants header-only (will fail at coverage but parse header)"
    timeout 60 gatk ValidateVariants -V "$f" -gvcf -R "$REF" 2>&1 | grep -E 'never saw|CHROM|Validation passed|covering it' | head -3
    COUNT=$((COUNT+1))
done

echo
echo "=== Compare: same test for one CD34, one B_lymph, one BJ_fib GVCF ==="
for f in /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/benchmarking/gvcfs/CD34_cord_blood/SRR8438255.g.vcf.gz \
         /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/benchmarking/gvcfs/B_lymphocyte/SRR8438251.g.vcf.gz \
         /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/benchmarking/gvcfs/BJ_fibroblast/SRR5365331.g.vcf.gz; do
    echo "--- File: $f"
    echo "Longest header line: $(bcftools view -h $f 2>/dev/null | awk '{print length}' | sort -n | tail -1)"
    timeout 60 gatk ValidateVariants -V "$f" -gvcf -R "$REF" 2>&1 | grep -E 'never saw|CHROM|Validation passed|covering it' | head -3
done
