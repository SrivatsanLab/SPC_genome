#!/usr/bin/env python3
"""
Convert cell type markers to CSV format with gene IDs from GTF
"""

import pandas as pd
import re
from collections import defaultdict

# Read the cell type markers file
print("Reading cell type markers...")
with open('bin/PolE_worm_pilot/cell_type_markers.txt', 'r') as f:
    lines = f.readlines()

# Parse the file
marker_data = []
for line in lines[1:]:  # Skip header
    parts = line.strip().split('\t')
    if len(parts) >= 3:
        subtype = parts[0]
        celltype = parts[1]
        genes_str = parts[2]

        # Split genes by comma
        genes = [g.strip() for g in genes_str.split(',')]

        for gene in genes:
            if gene:  # Skip empty strings
                marker_data.append({
                    'gene': gene,
                    'celltype': celltype,
                    'subtype': subtype
                })

print(f"Found {len(marker_data)} gene-cell type associations")

# Create DataFrame
df = pd.DataFrame(marker_data)

# Now read GTF to map gene names to WBGene IDs
print("\nReading GTF file to map gene names to WBGene IDs...")
gtf_file = '/shared/biodata/ngs/Reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf'

# Parse GTF to extract gene_id and gene_name mappings
gene_name_to_id = {}
gene_id_to_name = {}

with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue

        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue

        attributes = parts[8]

        # Extract gene_id and gene_name from any feature type
        gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
        gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)

        if gene_id_match and gene_name_match:
            gene_id = gene_id_match.group(1)
            gene_name = gene_name_match.group(1)

            # Store the mapping
            if gene_name not in gene_name_to_id:
                gene_name_to_id[gene_name] = gene_id
            if gene_id not in gene_id_to_name:
                gene_id_to_name[gene_id] = gene_name

print(f"Found {len(gene_name_to_id)} gene name to ID mappings")

# Map genes to their WBGene IDs
def get_gene_id(gene_name):
    # Try exact match
    if gene_name in gene_name_to_id:
        return gene_name_to_id[gene_name]

    # Try case-insensitive match
    gene_name_lower = gene_name.lower()
    for gn, gid in gene_name_to_id.items():
        if gn.lower() == gene_name_lower:
            return gid

    return None

df['gene_id'] = df['gene'].apply(get_gene_id)

# Report mapping stats
mapped = df['gene_id'].notna().sum()
total = len(df)
print(f"\nMapping results: {mapped}/{total} genes mapped to WBGene IDs ({100*mapped/total:.1f}%)")

# Show unmapped genes
unmapped = df[df['gene_id'].isna()]['gene'].unique()
if len(unmapped) > 0:
    print(f"\nUnmapped genes ({len(unmapped)}):")
    for gene in unmapped[:20]:
        print(f"  {gene}")
    if len(unmapped) > 20:
        print(f"  ... and {len(unmapped)-20} more")

# Set gene as index and save
df_out = df.set_index('gene')
output_file = 'bin/PolE_worm_pilot/cell_type_markers.csv'
df_out.to_csv(output_file)

print(f"\nSaved to: {output_file}")
print(f"Shape: {df_out.shape}")
print("\nFirst few rows:")
print(df_out.head(10))
