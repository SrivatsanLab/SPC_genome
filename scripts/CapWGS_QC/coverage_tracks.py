#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pyBigWig
import pysam
from tqdm import tqdm
from scipy import stats
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000
import sys


cell = sys.argv[1]
species = sys.argv[2]
amplification = sys.argv[3]

print("plotting:")
print(f"cell: {cell}")
print(f"species: {species}")
print(f"amplification: {amplification}")

color_dict={
    'mouse' : '#ff1a5e',
    'human' : '#80f15d'
}

color = color_dict[species]

print(f"using color: {color}")

bw = pyBigWig.open(f"{species}_sc_outs/sc_output/{cell}.bw")

# Retrieve chromosome sizes and order them
chrom_sizes = bw.chroms()
if species == 'human':
    chrom_order = [f"chr{i}_hg38" for i in range(1,23)]
    chrom_order += ['chrX_hg38']

elif species == 'mouse':
    chrom_order = [f"chr{i}_mm10" for i in range(1,20)]
    chrom_order += ['chrX_mm10']

bin_size = 1000
gap_size = 10000000  # Define a gap between chromosomes
binned_positions = []
binned_coverage = []
chrom_offsets = {}
curr_offset = 0

for chrom in chrom_order:
    chrom_len = chrom_sizes[chrom]
    print(chrom, ',', chrom_len)
    coverage = bw.values(chrom, 0, chrom_len)
    
    num_bins = len(coverage) // bin_size
    binned_pos = np.arange(curr_offset, curr_offset + chrom_len, bin_size)
    binned_cov = [np.nanmean(coverage[i:i + bin_size]) for i in range(0, len(coverage), bin_size)]
    
    if len(binned_cov) > len(binned_pos):
        binned_pos = np.append(binned_pos, curr_offset + chrom_len)
    
    chrom_offsets[chrom] = (curr_offset, curr_offset + chrom_len)
    curr_offset += chrom_len + gap_size  # Add a gap between chromosomes
    
    binned_positions.extend(binned_pos)
    binned_coverage.extend(binned_cov)

# bw.close()

fig, ax = plt.subplots(1, 1, figsize=(16, 2), dpi=600)

sns.lineplot(x=binned_positions[:len(binned_coverage)], y=binned_coverage, color=color, lw=1, ax=ax)
ax.set_yscale('log')
ax.set_ylim(0,15)
# ax.set_xlabel('Genomic Position', fontsize=10, fontweight='bold')
# fig.text(0.075, 0.5, 'Coverage', va='center', rotation='vertical', fontsize=12, fontweight='bold')

# Add chromosome labels
# for chrom, (start, end) in chrom_offsets.items():
    # mid = (start + end) / 2
    # ax.text(mid, ax.get_ylim()[1], chrom[:-5], fontsize=8, verticalalignment='top', rotation=45,horizontalalignment='center')
for chrom, (start,end) in chrom_offsets.items():
    ax.axvline(x=end+5000000, color='gray', linestyle='--', lw=0.5)
ax.set_xticks(ticks = [(start + end) / 2 for (start, end) in chrom_offsets.values()], labels = [chrom[:-5] for chrom in chrom_order],rotation=45)
#ax.set_yticks([])
# plt.savefig(f'Barn10/{cell}_track.png', dpi=600, bbox_inches='tight')
sns.despine()
#fig.savefig(f'coverage_tracks/WG_{cell}_{species}_{amplification}.png', bbox_inches='tight')
fig.savefig(f'scale_WG_{cell}_{species}_{amplification}.png', bbox_inches='tight')

if species == 'human':
    chrom='chr1_hg38'
elif species == 'mouse':
    chrom='chr1_mm10'

chrom_len = bw.chroms()[chrom]
coverage = bw.values(chrom, 0, chrom_len)
bw.close()

# Define the bin size (e.g., 1000 base pairs per bin)
bin_size = 1000

# Compute the number of bins
num_bins = len(coverage) // bin_size

# Downsample the coverage array by computing the mean for each bin
binned_positions = np.arange(0, chrom_len, bin_size)
binned_coverage = [
    np.nanmean(coverage[i:i + bin_size]) for i in range(0, len(coverage), bin_size)
]

# If there's an extra bin for the remainder, extend binned_positions accordingly
if len(binned_coverage) > len(binned_positions):
    binned_positions = np.append(binned_positions, chrom_length)
    
fig,ax = plt.subplots(1,1,figsize=(10,2))

sns.lineplot(x=binned_positions[:len(binned_coverage)+1], y=binned_coverage, color=color, lw=1, ax=ax)
ax.set_yscale('log')
ax.set_xlabel('Chr1', fontsize=12)
ax.set_ylim(0,15)

#ax.set_yticks([])
ax.set_xticks([])

sns.despine()
#fig.savefig(f'coverage_tracks/chr1_{cell}_{species}_{amplification}.png', bbox_inches='tight')
fig.savefig(f'scale_chr1_{cell}_{species}_{amplification}.png', bbox_inches='tight')
