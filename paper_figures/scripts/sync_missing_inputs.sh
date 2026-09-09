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

# K562 mutation accumulation: derived tables written by
# notebooks/K562_mut_accumulation.ipynb into <repo>/results/K562_mut_accumulation.
# Point K562_MUT_ACCUM_RESULTS at another checkout to pull them from there.
k562_src="${K562_MUT_ACCUM_RESULTS:-$repo_root/results/K562_mut_accumulation}"
k562_dst="$repo_root/paper_figures/data/K562_mut_accumulation"

if [[ -d "$k562_src" ]]; then
  mkdir -p "$k562_dst"
  for f in mutation_accumulation.csv spectrum.csv spectrum_background.csv \
           top10_EDT.csv full_EDT.csv; do
    if [[ -f "$k562_src/$f" ]]; then
      cp "$k562_src/$f" "$k562_dst/$f"
      echo "Copied K562_mut_accumulation/$f"
    else
      echo "Missing in notebook results: K562_mut_accumulation/$f"
    fi
  done
else
  echo "Run notebooks/K562_mut_accumulation.ipynb first, or set" \
       "K562_MUT_ACCUM_RESULTS to a checkout that has results/K562_mut_accumulation"
fi
