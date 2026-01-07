#!/bin/bash

###########################################################################################################################
# Compile Lorenz curves for coverage benchmarking datasets
# Usage: ./scripts/compile_lorenz_for_benchmarking.sh
###########################################################################################################################

set -euo pipefail

# Create results directory
RESULTS_DIR="results/coverage_benchmarking"
mkdir -p "${RESULTS_DIR}"

# Activate Python environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

echo "=========================================="
echo "Compiling Lorenz curves for benchmarking"
echo "=========================================="
echo ""

# Function to compile Lorenz curves
compile_lorenz() {
    local dataset_name=$1
    local input_dir=$2
    local output_file="${RESULTS_DIR}/${dataset_name}_lorenz_curves.csv"

    echo "Dataset: ${dataset_name}"
    echo "Input directory: ${input_dir}"
    echo "Output file: ${output_file}"

    python << EOFPYTHON
import sys
import os
import pandas as pd
import glob

input_dir = "${input_dir}"
output_file = "${output_file}"

# Find all Lorenz curve files
lorenz_files = sorted(glob.glob(f"{input_dir}/*_lorenz.csv"))

if len(lorenz_files) == 0:
    print(f"ERROR: No Lorenz curve files found in {input_dir}")
    sys.exit(1)

print(f"Found {len(lorenz_files)} Lorenz curve files")

# Read all files and combine
dfs = {}
for filepath in lorenz_files:
    barcode = os.path.basename(filepath).replace('_lorenz.csv', '')
    df = pd.read_csv(filepath)

    # Use the correct column names from lorenz.py output
    if 'cumulative_fraction_genome' in df.columns and 'cumulative_fraction_reads' in df.columns:
        dfs[barcode] = df.set_index('cumulative_fraction_genome')['cumulative_fraction_reads']
    else:
        print(f"WARNING: Unexpected column names in {filepath}: {df.columns.tolist()}")

# Combine into single dataframe
if dfs:
    combined_df = pd.DataFrame(dfs)
    combined_df.index.name = 'cumulative_fraction_genome'
    combined_df.to_csv(output_file)
    print(f"Successfully compiled {len(dfs)} Lorenz curves to {output_file}")
    print(f"Output shape: {combined_df.shape}")
else:
    print("ERROR: No valid Lorenz curves found")
    sys.exit(1)

EOFPYTHON

    echo "✓ Completed: ${dataset_name}"
    echo ""
}

# Compile HSC_enzyme_coverage Lorenz curves
compile_lorenz "HSC_enzyme_coverage" "data/HSC_enzyme_coverage/sc_outputs"

# Compile public data Lorenz curves
compile_lorenz "public_data" "data/coverage_benchmarking"

echo "=========================================="
echo "All Lorenz curve compilations complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - ${RESULTS_DIR}/HSC_enzyme_coverage_lorenz_curves.csv"
echo "  - ${RESULTS_DIR}/public_data_lorenz_curves.csv"
