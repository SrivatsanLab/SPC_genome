#!/bin/bash
###########################################################################################################################
# Split a reference FASTA index (.fai) into fixed-size regions, one per line.
#
# Usage:
#   scripts/CapGTA/make_regions.sh <reference.fa[.fai]> <region_size_bp> <output_regions.txt>
#
# Output format (one region per line, bcftools -r compatible):
#   I:1-1000000
#   I:1000001-2000000
#   II:1-1000000
#   ...
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference.fa[.fai]> <region_size_bp> <output_regions.txt>" >&2
    exit 1
fi

REFERENCE="$1"
REGION_SIZE="$2"
OUTPUT="$3"

# Accept either the FASTA or the .fai directly
if [ "${REFERENCE##*.}" = "fai" ]; then
    FAI="${REFERENCE}"
else
    FAI="${REFERENCE}.fai"
fi

if [ ! -f "${FAI}" ]; then
    echo "Error: .fai not found at ${FAI}. Run 'samtools faidx' on the FASTA first." >&2
    exit 1
fi

if ! [[ "${REGION_SIZE}" =~ ^[0-9]+$ ]] || [ "${REGION_SIZE}" -lt 1 ]; then
    echo "Error: region_size_bp must be a positive integer (got '${REGION_SIZE}')" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUTPUT}")"

awk -v size="${REGION_SIZE}" 'BEGIN {OFS=""} {
    chrom = $1
    length_bp = $2
    for (start = 1; start <= length_bp; start += size) {
        end = start + size - 1
        if (end > length_bp) end = length_bp
        print chrom, ":", start, "-", end
    }
}' "${FAI}" > "${OUTPUT}"

n_regions=$(wc -l < "${OUTPUT}")
echo "Wrote ${n_regions} regions to ${OUTPUT} (region size: ${REGION_SIZE} bp)"
