#!/bin/bash
# SPC Genome Project Installation Script
# This script sets up the computational environment for the project

set -e  # Exit on error

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="spc_genome"
ENV_DIR="${PROJECT_ROOT}/env"

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
if [ -d "${ENV_DIR}" ]; then
    echo "WARNING: Environment directory already exists at ${ENV_DIR}"
    read -p "Do you want to remove it and reinstall? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        rm -rf "${ENV_DIR}"
    else
        echo "Keeping existing environment. Skipping installation."
        echo ""
        echo "To activate the environment, run:"
        echo "  micromamba activate ${ENV_DIR}"
        exit 0
    fi
fi

# Create environment
echo "Creating micromamba environment from environment.yml..."
echo "This may take several minutes..."
echo ""

micromamba create -y -p "${ENV_DIR}" -f "${PROJECT_ROOT}/environment.yml"

echo ""
echo "========================================="
echo "✓ Installation complete!"
echo "========================================="
echo ""
echo "To activate the environment, run:"
echo "  micromamba activate ${ENV_DIR}"
echo ""
echo "Or add this to your script:"
echo "  eval \"\$(micromamba shell hook --shell bash)\""
echo "  micromamba activate ${ENV_DIR}"
echo ""
echo "Next steps:"
echo "  1. Edit config.yaml with your data paths"
echo "  2. Run preprocessing: ./WGS_PP.sh -o sample_name -1 read1.fq.gz -2 read2.fq.gz -r read_count"
echo "  3. Run analysis notebooks in JupyterLab"
echo ""
