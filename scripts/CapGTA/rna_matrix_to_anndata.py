#!/usr/bin/env python3
"""
Convert RNA count matrix to AnnData object with matching obs_names to variant AnnData
"""

import sys
import pandas as pd
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix

if len(sys.argv) != 3:
    print("Usage: python rna_matrix_to_anndata.py <count_matrix.csv> <output.h5ad>")
    sys.exit(1)

count_matrix_file = sys.argv[1]
output_file = sys.argv[2]

print(f"Loading RNA count matrix from: {count_matrix_file}")

# Load count matrix (genes x cells)
counts = pd.read_csv(count_matrix_file, index_col=0)

print(f"Loaded matrix: {counts.shape[0]} genes x {counts.shape[1]} cells")

# Transpose to cells x genes (AnnData convention)
counts_T = counts.T

# Create AnnData object
adata = ad.AnnData(csr_matrix(counts_T.values, dtype=np.float32))
adata.obs_names = counts_T.index
adata.var_names = counts_T.columns

# Add basic metadata
adata.obs['n_counts'] = np.array(counts_T.sum(axis=1))
adata.obs['n_genes'] = np.array((counts_T > 0).sum(axis=1))

print(f"\nAnnData summary:")
print(f"  Observations (cells): {adata.n_obs}")
print(f"  Variables (genes): {adata.n_vars}")
print(f"  Mean counts per cell: {adata.obs['n_counts'].mean():.1f}")
print(f"  Mean genes per cell: {adata.obs['n_genes'].mean():.1f}")

# Save to h5ad
print(f"\nSaving to: {output_file}")
adata.write_h5ad(output_file)

print("Done!")
