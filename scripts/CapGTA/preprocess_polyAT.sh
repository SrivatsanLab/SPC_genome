#!/bin/bash
#SBATCH -J preprocess_polyAT
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 1:00:00

###########################################################################################################################
# Run preprocess_scanpy.py on a polyA/T RNA count matrix (parity with the arm='current' SoupX-input pipeline).
#
# Usage:
#   sbatch scripts/CapGTA/preprocess_polyAT.sh <input.h5ad> <output.h5ad>
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.h5ad> <output.h5ad>" >&2
    exit 1
fi

IN_H5AD="$1"
OUT_H5AD="$2"

REPO="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome"
GTF="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
JUNC_CSV="${REPO}/results/worm6_final/GEX/gene_junctions.csv"
LEN_CSV="${REPO}/data/reference/gene_lengths.csv"

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

python "${REPO}/scripts/CapGTA/preprocess_scanpy.py" \
    "${IN_H5AD}" "${OUT_H5AD}" \
    --gtf "${GTF}" \
    --drop-cuticle \
    --junction-csv "${JUNC_CSV}" \
    --gene-lengths-csv "${LEN_CSV}"
