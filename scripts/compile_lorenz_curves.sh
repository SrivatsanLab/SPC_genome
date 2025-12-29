#!/bin/bash

###########################################################################################################################
# Compile individual Lorenz curve CSV files into a single CSV with columns for each cell
# Usage: ./scripts/compile_lorenz_curves.sh <sample_name> <sc_outputs_dir>
###########################################################################################################################

set -euo pipefail

SAMPLE="$1"
SC_OUTPUTS_DIR="$2"

# Create results directory
RESULTS_DIR="results/${SAMPLE}"
mkdir -p "${RESULTS_DIR}"

OUTPUT_FILE="${RESULTS_DIR}/lorenz_curves_compiled.csv"

echo "Compiling Lorenz curves for sample: ${SAMPLE}"
echo "Reading from: ${SC_OUTPUTS_DIR}"
echo "Output file: ${OUTPUT_FILE}"

# Activate Python environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

# Python script to compile Lorenz curves
python << 'EOFPYTHON'
import sys
import os
import pandas as pd
import glob

sc_outputs_dir = sys.argv[1]
output_file = sys.argv[2]

# Find all Lorenz curve files
lorenz_files = sorted(glob.glob(f"{sc_outputs_dir}/*_lorenz.csv"))

if len(lorenz_files) == 0:
    print(f"ERROR: No Lorenz curve files found in {sc_outputs_dir}")
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

EOFPYTHON "${SC_OUTPUTS_DIR}" "${OUTPUT_FILE}"

echo "Compilation complete!"
