#!/bin/bash
# SPC Genome Project Environment Setup
# This script verifies dependencies and provides setup instructions

set -e  # Exit on error

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="spc_genome"

echo "========================================="
echo "SPC Genome Project Setup"
echo "========================================="
echo ""

# Check if micromamba is installed
if ! command -v micromamba &> /dev/null; then
    echo "ERROR: micromamba is not installed or not in PATH"
    echo ""
    echo "To install micromamba, run:"
    echo "  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba"
    echo "  ./bin/micromamba shell init -s bash -p ~/micromamba"
    echo ""
    exit 1
fi

echo "✓ Found micromamba: $(which micromamba)"
echo ""

# Check if environment already exists
if micromamba env list | grep -q "^${ENV_NAME} "; then
    echo "✓ Environment '${ENV_NAME}' already exists"
    echo ""
    echo "To activate it:"
    echo "  micromamba activate ${ENV_NAME}"
    echo ""
    echo "To update it:"
    echo "  micromamba update -n ${ENV_NAME} --all"
    echo ""
else
    echo "Environment '${ENV_NAME}' not found."
    echo ""
    echo "To create it, run:"
    echo "  micromamba create -n ${ENV_NAME} -f environment.yml"
    echo ""
    echo "Then activate it:"
    echo "  micromamba activate ${ENV_NAME}"
    echo ""
fi

echo "========================================="
echo "Next steps:"
echo "========================================="
echo "  1. Edit config.yaml with your data paths"
echo "  2. Run preprocessing: ./WGS_PP.sh -o sample_name -1 read1.fq.gz -2 read2.fq.gz -r read_count"
echo "  3. Run analysis notebooks in JupyterLab"
echo ""
