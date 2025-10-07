#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyBigWig
import sys
from tqdm import tqdm
import os

def get_chrom_lengths(chrom_sizes_file):
    chrom_lengths = {}
    with open(f'{chrom_sizes_file}', 'r') as fin:
        for line in fin:
            line = line.strip()
            line = line.split()
            chrom_lengths[line[0]] = int(line[1])
    return chrom_lengths

def get_coverage(bigwig_file, chrom_lengths):
    bin_coverage = [0]
    bw = pyBigWig.open(bigwig_file)
    for chrom in chrom_lengths.keys():
        intervals = bw.intervals(chrom)
        for start, end, coverage in intervals:
            bin_coverage.append(coverage)
    bw.close()
    return bin_coverage

def create_uniform_bins(chrom_lengths, binsize):
    uniform_bins = []
    for chrom, length in chrom_lengths.items():
        for start in range(0, length, binsize):
            end = min(start + binsize, length)
            uniform_bins.append((chrom, start, end))
    return uniform_bins

def map_coverage_to_bins(bigwig_file, uniform_bins):
    coverage = []
    bw = pyBigWig.open(bigwig_file)
    for chrom, start, end in tqdm(uniform_bins):
        interval_value = bw.values(chrom, start, end, numpy=False)[0]
        coverage.append(interval_value if interval_value is not None else 0)  # Replace None with 0
    bw.close()
    return coverage

# Get SLURM task ID to determine which file to process
# task_id = int(os.getenv('SLURM_ARRAY_TASK_ID', '0'))
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID', '0')) - 1
cell_list_file = sys.argv[1]  # File with the list of cells
working_dir = sys.argv[2]
out_dir = sys.argv[3]

human_chroms = get_chrom_lengths('hg38.chrom.sizes')

cells_df = pd.read_csv(cell_list_file, sep='\t', header=None)

# Extract the BigWig file for this job based on task_id
barcode = cells_df.iloc[task_id][0]

bigwig_file = f'{working_dir}/{barcode}.bw'

bins = create_uniform_bins(human_chroms, 1000)
coverage = map_coverage_to_bins(f"{bigwig_file}", bins)

output = pd.Series(coverage)
output.to_csv(f"{out_dir}/{barcode}.csv", index=False)
