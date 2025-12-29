#!/bin/bash

###########################################################################################################################
# Compile individual Gini coefficient files into a single CSV
# Usage: ./scripts/compile_gini_coefficients.sh <sample_name> <sc_outputs_dir>
###########################################################################################################################

set -euo pipefail

SAMPLE="$1"
SC_OUTPUTS_DIR="$2"

# Create results directory
RESULTS_DIR="results/${SAMPLE}"
mkdir -p "${RESULTS_DIR}"

OUTPUT_FILE="${RESULTS_DIR}/gini_coefficients.csv"

echo "Compiling Gini coefficients for sample: ${SAMPLE}"
echo "Reading from: ${SC_OUTPUTS_DIR}"
echo "Output file: ${OUTPUT_FILE}"

# Activate Python environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

# Python script to compile Gini coefficients
python << 'EOFPYTHON'
import sys
import os
import pandas as pd
import glob

sc_outputs_dir = sys.argv[1]
output_file = sys.argv[2]

# Find all Gini coefficient files
gini_files = sorted(glob.glob(f"{sc_outputs_dir}/*_gini.txt"))

if len(gini_files) == 0:
    print(f"ERROR: No Gini coefficient files found in {sc_outputs_dir}")
    sys.exit(1)

print(f"Found {len(gini_files)} Gini coefficient files")

# Read all files and compile
data = []
for filepath in gini_files:
    barcode = os.path.basename(filepath).replace('_gini.txt', '')

    # Read the Gini coefficient from the file
    with open(filepath, 'r') as f:
        gini_value = float(f.read().strip())

    data.append({'barcode': barcode, 'gini_coefficient': gini_value})

# Create dataframe and save
df = pd.DataFrame(data)
df.to_csv(output_file, index=False)

print(f"Successfully compiled {len(data)} Gini coefficients to {output_file}")
print(f"Gini coefficient range: {df['gini_coefficient'].min():.4f} - {df['gini_coefficient'].max():.4f}")

EOFPYTHON "${SC_OUTPUTS_DIR}" "${OUTPUT_FILE}"

echo "Compilation complete!"
