#!/usr/bin/env python
"""
Generate all possible combinations of 3 cells from 3 adjacent subpopulations
for coalescence testing.
"""
import argparse
import pandas as pd
import numpy as np
from itertools import product
import anndata as ad

parser = argparse.ArgumentParser(description="Generate cell combinations for coalescence testing")
parser.add_argument("-d", "--adata", required=True, help="Path to the AnnData object (.h5ad)")
parser.add_argument("-o", "--output", required=True, help="Output CSV file for combinations")
parser.add_argument("-p", "--populations", required=True, nargs=3, help="Three adjacent populations (e.g., P_0 P_1 P_2)")
args = parser.parse_args()

# Load AnnData object
print(f"Loading AnnData from {args.adata}...")
adata = ad.read_h5ad(args.adata)

# Extract cells from each population
pop0, pop1, pop2 = args.populations
cells_pop0 = adata.obs[adata.obs['pop'] == pop0].index.tolist()
cells_pop1 = adata.obs[adata.obs['pop'] == pop1].index.tolist()
cells_pop2 = adata.obs[adata.obs['pop'] == pop2].index.tolist()

print(f"Population {pop0}: {len(cells_pop0)} cells")
print(f"Population {pop1}: {len(cells_pop1)} cells")
print(f"Population {pop2}: {len(cells_pop2)} cells")

# Generate all combinations
print("Generating all combinations...")
combinations = list(product(cells_pop0, cells_pop1, cells_pop2))
print(f"Total combinations: {len(combinations)}")

# Save to CSV
df = pd.DataFrame(combinations, columns=['C0', 'C1', 'C2'])
df.to_csv(args.output, index=False)
print(f"Saved combinations to {args.output}")
