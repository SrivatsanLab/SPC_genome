#!/usr/bin/env python3

import pyBigWig
import pandas as pd
import numpy as np
import sys
import os

# Function to calculate genome coverage for a single BigWig file
def calculate_coverage(bigwig_file, chrom_lengths):
    bw = pyBigWig.open(bigwig_file)

    total_covered = 0
    total_bases = 0

    for chrom, length in chrom_lengths.items():
    
    
        if chrom in bw.chroms():
            coverage = bw.values(chrom, 0, length, numpy=True)
            covered_bases = np.count_nonzero(coverage > 0)
            total_covered += covered_bases
            total_bases += length

    proportion_covered = total_covered / total_bases if total_bases > 0 else 0
    return proportion_covered

# Get SLURM task ID to determine which file to process
# task_id = int(os.getenv('SLURM_ARRAY_TASK_ID', '0'))
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID', '0')) - 1
cell_list_file = sys.argv[1]  # File with the list of cells passed as argument
chrom_sizes = sys.argv[2]  # path to UCSC genome browser formatted chromosome lengths file
temp_dir = sys.argv[3]

# Get Chrom lengths
chrom_lengths = {}
with open(f'{chrom_sizes}', 'r') as fin:
    for line in fin:
        line = line.strip()
        line = line.split()
        chrom_lengths[line[0]] = int(line[1])
# print(chrom_lengths['chrX'], chrom_lengths['chrY'])

# Read the file listing BigWig file paths
cells_df = pd.read_csv(cell_list_file, sep='\t', header=None)

# Extract the BigWig file for this job based on task_id
barcode = cells_df.iloc[task_id][0]
bigwig_file = f'{temp_dir}/{barcode}.bw'

print(f"Processing barcode {barcode} at task ID {task_id + 1}")

# Calculate coverage for the current BigWig file
proportion_covered = calculate_coverage(bigwig_file, chrom_lengths)

# Save results to a file, using the BigWig filename for output
output_file = f"{temp_dir}/{os.path.basename(bigwig_file).replace('.bw', '')}_coverage.csv"
with open(output_file, 'w') as f:
    f.write(f"{barcode},{proportion_covered}\n")