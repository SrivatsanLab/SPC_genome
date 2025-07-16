#!/bin/bash
#SBATCH --job-name=K562_createDB_test       # Job name
#SBATCH --output=K562_createDB_test_%j.out  # Output file
#SBATCH --error=K562_createDB_test_%j.err   # Error file
#SBATCH --cpus-per-task=16               # Number of CPUs
#SBATCH --mem=16G                        # Memory allocation
#SBATCH --time=12:00:00                  # Maximum runtime

gatk GenomicsDBImport \
    --genomicsdb-workspace-path K562_genomicsdb \
    --sample-name-map sample_map_K562.txt \
    --reader-threads 16 \
    -L ../BIRC6.bed \
    -L ../FANCC.bed \
    -L ../BRCA1.bed \
    -L ../AKT3_chr1.bed
