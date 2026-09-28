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
  for f in mutation_accumulation.csv spectrum.csv spectrum_background.csv; do
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

# K562 consensus trees: newick written by notebooks/K562_tree.ipynb into
# <repo>/results/K562_tree/trees. Point K562_TREE_RESULTS at another checkout
# to pull them from there.
tree_src="${K562_TREE_RESULTS:-$repo_root/results/K562_tree/trees}"
tree_dst="$repo_root/paper_figures/data/K562_tree"

if [[ -d "$tree_src" ]]; then
  mkdir -p "$tree_dst"
  for f in grouped_bootstrap_consensus.newick grouped_bootstrap_consensus_upgma.newick; do
    if [[ -f "$tree_src/$f" ]]; then
      cp "$tree_src/$f" "$tree_dst/$f"
      echo "Copied K562_tree/$f"
    else
      echo "Missing in notebook results: K562_tree/trees/$f"
    fi
  done
else
  echo "Run notebooks/K562_tree.ipynb first, or set K562_TREE_RESULTS to a" \
       "checkout that has results/K562_tree/trees"
fi

# worm6 haplotype UMAP: coordinates exported from
# notebooks/worm6_final_haplotype_assignment.ipynb. The notebook does not write
# these by default - see the header of worm6_haplotype_umap.R for the one line
# that does. Point WORM6_RESULTS at another checkout to pull from there.
worm6_src="${WORM6_RESULTS:-$repo_root/results/worm6_final/figures}"
worm6_dst="$repo_root/paper_figures/data/worm6_final"

mkdir -p "$worm6_dst"
for f in clean_ind_umap_coords.csv coassay_readcount_counts.csv; do
  if [[ -f "$worm6_src/$f" ]]; then
    cp "$worm6_src/$f" "$worm6_dst/$f"
    echo "Copied worm6_final/$f"
  else
    echo "Missing worm6 input: $worm6_src/$f" \
         "(generate it with paper_figures/scripts/export_worm6_*.py)"
  fi
done

# K562 single-cell CNV (sc_PolE_novaseq cells): the GC-corrected, blacklisted
# AneuFinder run (scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh) and the
# per-cell ploidy and CNV event tree built from it
# (paper_figures/scripts/K562_cnv_event_tree.R). These replace the original
# sc_PolE_novaseq AneuFinder_output/result.csv and full_meta.csv ploidy, which
# had no GC correction and mislabelled columns. Point K562_CNV_RESULTS at
# another checkout to pull from there.
cnv_src="${K562_CNV_RESULTS:-$repo_root/results/K562_tree}"
cnv_dst="$repo_root/paper_figures/data/Anneufinder"

mkdir -p "$cnv_dst"
for pair in "aneufinder_gc/output/result.csv:result.csv" \
            "sc_trees/cnv_cell_ploidy_gc.csv:cell_ploidy.csv" \
            "sc_trees/cnv_event_tree_gc.newick:cnv_event_tree_gc.newick"; do
  src="$cnv_src/${pair%%:*}"
  dst="$cnv_dst/${pair##*:}"
  if [[ -f "$src" ]]; then
    cp "$src" "$dst"
    echo "Copied Anneufinder/${pair##*:}"
  else
    echo "Missing K562 CNV input: $src" \
         "(run scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh, then K562_cnv_event_tree.R)"
  fi
done

# Per-cell SBS tree used by the CNV section of draw_trees.R
sc_trees_src="${K562_SC_TREES:-/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/results}"
sc_trees_dst="$repo_root/paper_figures/data/SingleCellTrees"
mkdir -p "$sc_trees_dst"
if [[ -f "$sc_trees_src/sc_test.newick" ]]; then
  cp "$sc_trees_src/sc_test.newick" "$sc_trees_dst/sc_test.newick"
  echo "Copied SingleCellTrees/sc_test.newick"
else
  echo "Missing $sc_trees_src/sc_test.newick; set K562_SC_TREES"
fi
