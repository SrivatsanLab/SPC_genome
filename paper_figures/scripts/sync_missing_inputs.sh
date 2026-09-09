#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
legacy_root="$(cd "$repo_root/.." && pwd)"

mkdir -p "$repo_root/paper_figures/data/external"

if [[ -f "$legacy_root/data/encode_variant_annotation/somatic_max_peaks.tsv" ]]; then
  cp "$legacy_root/data/encode_variant_annotation/somatic_max_peaks.tsv" \
     "$repo_root/paper_figures/data/encode_variant_annotation/somatic_max_peaks.tsv"
  echo "Copied somatic_max_peaks.tsv"
else
  echo "Missing in legacy data: data/encode_variant_annotation/somatic_max_peaks.tsv"
fi

if [[ -f "$legacy_root/data/Single_cell_bottlenecking_summary_statistics/filtered_bulk_vaf.csv" && ! -f "$repo_root/paper_figures/data/Single_cell_bottlenecking_summary_statistics/filtered_bulk_vaf.csv.gz" ]]; then
  gzip -c "$legacy_root/data/Single_cell_bottlenecking_summary_statistics/filtered_bulk_vaf.csv" \
    > "$repo_root/paper_figures/data/Single_cell_bottlenecking_summary_statistics/filtered_bulk_vaf.csv.gz"
  echo "Copied filtered_bulk_vaf.csv.gz"
fi

if [[ -f "$legacy_root/data/external/COSMIC_v3.4_SBS_GRCh38.txt" ]]; then
  cp "$legacy_root/data/external/COSMIC_v3.4_SBS_GRCh38.txt" \
     "$repo_root/paper_figures/data/external/COSMIC_v3.4_SBS_GRCh38.txt"
  echo "Copied COSMIC_v3.4_SBS_GRCh38.txt"
else
  echo "Provide COSMIC_v3.4_SBS_GRCh38.txt at paper_figures/data/external/"
fi
