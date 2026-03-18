#!/bin/bash
#SBATCH -J compile_qc
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 30

###########################################################################################################################
# Compile QC results from multiple sources into summary CSV files
#
# This script aggregates:
# 1. Lorenz curves (coverage uniformity)
# 2. Benchmarking QC metrics (alignment, GC bias, duplicates, coverage)
#
# Outputs summary CSV files to results directory for easy analysis
###########################################################################################################################

set -euo pipefail

SC_OUTPUTS_DIR="$1"
QC_METRICS_DIR="$2"
RESULTS_DIR="$3"

# Activate python environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

echo "=========================================="
echo "Compiling QC Results"
echo "=========================================="
echo "Single-cell BAM directory: ${SC_OUTPUTS_DIR}"
echo "QC metrics directory: ${QC_METRICS_DIR}"
echo "Results directory: ${RESULTS_DIR}"
echo "=========================================="
echo ""

##########################################################################################################################
# Compile Lorenz curves
##########################################################################################################################

echo "Compiling Lorenz curves..."

LORENZ_OUTPUT="${RESULTS_DIR}/compiled_lorenz_curves.csv"

if ls "${SC_OUTPUTS_DIR}"/*_lorenz.csv 1> /dev/null 2>&1; then
    python scripts/CapWGS_QC/compile_lorenz.py \
        "${SC_OUTPUTS_DIR}" \
        "${LORENZ_OUTPUT}"

    echo "✓ Lorenz curves compiled: ${LORENZ_OUTPUT}"
else
    echo "⚠ No Lorenz curve files found in ${SC_OUTPUTS_DIR}"
fi

echo ""

##########################################################################################################################
# Compile benchmarking QC metrics
##########################################################################################################################

echo "Compiling benchmarking QC metrics..."

QC_OUTPUT="${RESULTS_DIR}/compiled_qc_metrics.csv"

if ls "${QC_METRICS_DIR}"/*_alignment_metrics.txt 1> /dev/null 2>&1; then
    python scripts/CapWGS_QC/compile_qc_metrics.py \
        "${QC_METRICS_DIR}" \
        "${QC_OUTPUT}"

    echo "✓ QC metrics compiled: ${QC_OUTPUT}"
else
    echo "⚠ No QC metrics files found in ${QC_METRICS_DIR}"
fi

echo ""

##########################################################################################################################
# Summary
##########################################################################################################################

echo "=========================================="
echo "QC Compilation Complete!"
echo "=========================================="
echo "Output files:"
if [ -f "${LORENZ_OUTPUT}" ]; then
    echo "  ✓ ${LORENZ_OUTPUT}"
fi
if [ -f "${QC_OUTPUT}" ]; then
    echo "  ✓ ${QC_OUTPUT}"
fi
echo "=========================================="
