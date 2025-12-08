#!/usr/bin/env python
import argparse 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from itertools import product
from cyvcf2 import VCF
from tqdm import tqdm
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process cell combinations and classify coalescence types.")
parser.add_argument("-c", "--comb", required=True, help="Path to the file containing combinations of cells from adjacent populations (CSV format)")
parser.add_argument("-d", "--adata", required=True, help="Path to the AnnData object (.h5ad)")
parser.add_argument("-o", "--output", default=".", help="Output directory for results (default: current directory)")
parser.add_argument("-t", "--task_id", required=True, type=int, help="SLURM task ID")
parser.add_argument("-n", "--num_tasks", required=True, type=int, help="Total number of tasks (SLURM job array size)")
args = parser.parse_args()

# Load combinations
combinations = pd.read_csv(args.comb).values.tolist()  # Expecting a CSV with columns C0, C1, C2

# Divide combinations into chunks for parallel processing
chunk_size = len(combinations) // args.num_tasks
start_idx = args.task_id * chunk_size
end_idx = start_idx + chunk_size if args.task_id != args.num_tasks - 1 else len(combinations)

# Slice combinations for this task
combinations_subset = combinations[start_idx:end_idx]

# Load AnnData object
adata = ad.read_h5ad(args.adata)

tree_counts = {
    "CO_210": 0,  # Coalescence
    "NC_120": 0,  # No coalescence
    "NC_012": 0,  # No coalescence
}
coalescent_combos = []

with tqdm(total=len(combinations_subset), desc="Processing Cell Combinations", unit="comb") as pbar:
    for combo in combinations_subset:
        subset = adata[combo, :].copy()
        shared_observed_variants = (subset.layers["DP"] > 0).sum(axis=0) == len(combo)
        subset = subset[:, shared_observed_variants.A1 if hasattr(shared_observed_variants, "A1") else shared_observed_variants]
        
        X = subset.X.toarray() if hasattr(subset.X, "toarray") else subset.X

        # Compute shared genes for each pair
        shared_C0_C1 = ((X[0] > 0) & (X[1] > 0)).sum()
        shared_C0_C2 = ((X[0] > 0) & (X[2] > 0)).sum()
        shared_C1_C2 = ((X[1] > 0) & (X[2] > 0)).sum()

        if shared_C0_C1 > shared_C1_C2 and shared_C0_C1 > shared_C0_C2:  # supports CO_210
            tree_counts["CO_210"] += 1
            coalescent_combos.append(combo)
        elif shared_C0_C2 > shared_C0_C1 and shared_C0_C2 > shared_C1_C2:  # supports NC_120
            tree_counts["NC_120"] += 1
        elif shared_C1_C2 > shared_C0_C1 and shared_C1_C2 > shared_C0_C2:  # supports NC_012
            tree_counts["NC_012"] += 1

        pbar.update(1)

# Print final counts
print("Tree structure counts:", tree_counts)

# Ensure output directory exists
os.makedirs(args.output, exist_ok=True)

# Save tree counts to CSV
counts_df = pd.DataFrame.from_dict(tree_counts, orient='index', columns=['Count'])
counts_df.to_csv(os.path.join(args.output, f"tree_counts_{args.task_id}.csv"))

# Save coalescent combinations to TXT
with open(os.path.join(args.output, f"coalescent_combos_{args.task_id}.txt"), "w") as f:
    for combo in coalescent_combos:
        f.write(",".join(combo) + "\n")
