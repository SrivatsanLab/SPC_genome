#!/bin/bash
#SBATCH --job-name=K562_joint_calling_test       # Job name
#SBATCH --output=K562_joint_calling_test_%j.out  # Output file
#SBATCH --error=K562_joint_calling_test_%j.err   # Error file
#SBATCH --cpus-per-task=16               # Number of CPUs
#SBATCH --mem=16G                        # Memory allocation
#SBATCH --time=12:00:00                  # Maximum runtime

gatk GenotypeGVCFs \
    -R ../hg38_new/hg38_new.fasta \
    -V gendb://K562_genomicsdb \
    -O K562_joint_genotyped_test.vcf

