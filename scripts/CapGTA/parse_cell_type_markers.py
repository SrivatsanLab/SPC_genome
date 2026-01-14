#!/usr/bin/env python3
"""
Parse C. elegans cell type marker genes from Packer et al. 2019 supplementary table
and create a Python dictionary for cell type annotation
"""

import pandas as pd
import sys

# Read the Excel file
excel_file = "bin/PolE_worm_pilot/aax1971_packer_tables_s1_to_s6_s9_s12_s13_s15_s16.xlsx"

print("Reading Excel file...")
print(f"File: {excel_file}")
print()

# Read first sheet
df = pd.read_excel(excel_file, sheet_name=0)

print(f"First sheet shape: {df.shape}")
print(f"Columns: {list(df.columns)}")
print()
print("First few rows:")
print(df.head(10))
print()

# Try to identify the structure
print("\n=== Analyzing table structure ===")
print(f"Number of rows: {len(df)}")
print(f"Number of columns: {len(df.columns)}")
print()

# Check for common patterns in cell type marker tables
print("Checking for cell type/cluster columns...")
for col in df.columns:
    print(f"  {col}: {df[col].dtype}, unique values: {df[col].nunique()}")
    if df[col].nunique() < 200:  # Likely categorical
        print(f"    Sample values: {df[col].dropna().unique()[:5].tolist()}")
print()

# Try to create a dictionary structure
print("\n=== Creating marker dictionary ===")

# Attempt 1: If there's a cell type column and gene column
if 'Cell type' in df.columns or 'Celltype' in df.columns or 'cell_type' in df.columns:
    celltype_col = [c for c in df.columns if 'type' in c.lower()][0]
    gene_col = [c for c in df.columns if 'gene' in c.lower()][0] if any('gene' in c.lower() for c in df.columns) else df.columns[0]

    print(f"Detected cell type column: {celltype_col}")
    print(f"Detected gene column: {gene_col}")

    # Group genes by cell type
    marker_dict = {}
    for cell_type in df[celltype_col].dropna().unique():
        genes = df[df[celltype_col] == cell_type][gene_col].dropna().tolist()
        marker_dict[cell_type] = genes

    print("\nMarker dictionary created!")
    print(f"Number of cell types: {len(marker_dict)}")
    for ct, genes in list(marker_dict.items())[:5]:
        print(f"  {ct}: {len(genes)} genes (e.g., {genes[:3]})")

# Attempt 2: If structure is different, show more info
else:
    print("Could not auto-detect structure. Here's more detailed info:")
    print("\nFirst 20 rows of all columns:")
    print(df.head(20).to_string())

print("\n" + "="*80)
print("To create a dictionary for your notebook, use one of these formats:")
print("="*80)
print("""
# Format 1: Dictionary of lists
marker_genes = {
    'Neuron': ['unc-86', 'egl-3', ...],
    'Muscle': ['myo-3', 'unc-54', ...],
    ...
}

# Format 2: Dictionary with scores
marker_genes = {
    'Neuron': {
        'unc-86': 0.95,
        'egl-3': 0.88,
        ...
    },
    ...
}
""")
