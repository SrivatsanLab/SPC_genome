# Paper Figures Reproducibility Bundle

This directory consolidates figure-related scripts and input data from legacy `../code` and `../data` into a GitHub-friendly workflow.

## Structure
- `paper_figures/scripts/`: plotting and analysis scripts (ported from legacy code with repo-relative paths).
- `paper_figures/scripts/spectrum_utils.R`: shared 96-context spectrum plotting (`annotate_spectrum()`, `plot_spectra()`, COSMIC ordering and palette), extracted from `mutation_spectra.R`. Sourced by other scripts, not run directly.
- `paper_figures/data/`: curated figure input data.
- `paper_figures/output/`: generated plots/tables.
- `paper_figures/run_all_figures.R`: orchestrator to run all scripts and write `output/run_summary.csv`.

## Run
From repository root:

```bash
Rscript paper_figures/run_all_figures.R
```

## Notebook-derived inputs
Some panels are re-renders of analysis done in the Python notebooks under
`notebooks/`. The notebook owns the variant filtering and writes derived tables
to `<repo>/results/<experiment>/`; the R script here consumes those tables and
renders the panels in the house style. Filtering thresholds are documented in
the header of each script.

- `K562_consensus_tree.R` <- `notebooks/K562_tree.ipynb` cell 97
  (`grouped_bootstrap_consensus.newick`). Re-draws the baltic/matplotlib tree
  with ggtree in the `draw_trees.R` style.
- `K562_mut_accumulation.R` <- `notebooks/K562_mut_accumulation.ipynb`
  (`spectrum.csv`, `spectrum_background.csv`). Renders the de novo spectra
  pooled by construct (AAVS/PolE), the per-sample facet grid, and the
  ancestral background spectrum.

- K562 single-cell copy number (`anneufinder_plot.R`, the CNV section of
  `draw_trees.R`, and the `K562_cnv_*` / `K562_tree_*` scripts) uses the
  GC-corrected, blacklisted AneuFinder run of the 1000 sc_PolE_novaseq cells
  (`scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh`, 1 Mb bins, edivisive),
  staged into `data/Anneufinder/` as `result.csv`, `cell_ploidy.csv` and
  `cnv_event_tree_gc.newick`. It replaces the original sc_PolE_novaseq
  `AneuFinder_output/result.csv` and the `ploidy` column of `full_meta.csv`:
  that run had no GC correction, which put about a quarter of cells on the
  wrong overall scale (the apparent genome doublings; see
  `K562_gc_toggle_ploidy.R`), and its cell columns were mislabelled (column i
  held the i-th cell in alphabetical order).

- `K562_offpeak_ploidy.R` classifies the 110 cells off the dominant ploidy
  as artefacts or real copy-number changes, from the run's QC and GC-corrected
  read depth, and recomputes agreement with the K562 reference profile. It
  needs the AneuFinder R library (see the run script) and the run's
  `MODELS/` and `binned-GC/` folders. Results: `docs/K562_aneufinder_gc_review.md`.

- `K562_cnv_reference_correlation.R` reproduces the notebook's per-cell
  correlation with the K562 reference profile (cells 134-148) on the
  GC-corrected run, and writes the unmappable-bin mask
  (`results/K562_tree/masks/unmappable_bins_<binsize>.csv`) that the notebook
  now reads.
- `K562_cnv_heatmap_by_bin.R` draws the copy-number heatmap with bins as rows
  (genome order) and hierarchically clustered cells as columns; legends go to a
  separate `*_legends.svg`.
- `K562_cnv_population_scan.R` looks for copy-number changes specific to one
  population (P_0-P_5).
- `K562_cnv_cell_similarity.R` asks whether highly similar cells (GC-adjusted
  depth) share a population, controlling for barcode-well effects, and lists
  groups of near-identical cells.
- These take `K562_AF_DIR` to run on the 500 kb / 250 kb AneuFinder runs
  (`run_aneufinder_K562_sc_PolE_gc.sh 5e5`).

- Cell cycle and replication-corrected copy number (full write-up:
  `docs/K562_cell_cycle_cnv.md`): `K562_replication_timing_coverage.R`
  (phase per cell from coverage vs K562 Repli-seq, each chromosome a
  replicate), `K562_cell_cycle_profile.R` (FACS-style profile),
  `K562_replication_ideogram.R` (single-cell ideogram, zooms, arm scatter;
  uses ggideogram), `K562_cnv_distance_tree.R` and `K562_cnv_tree_sisters.R`
  (distance tree and same-population sister test on the corrected calls),
  `K562_sbs_tree_cell_cycle.R` (SBS tree vs replication state),
  `K562_aneufinder_clean_heatmap.R` (bare AneuFinder heatmap without artefact
  cells). Upstream inputs are built by `scripts/utils/download_replication_timing.sh`,
  `replication_timing_to_hg38.R`, `fragment_gc_*.R` and
  `replication_regress_segment.R`.
- SBS vs replication timing: `K562_sbs_tree_sites.R` (recovers the 700
  sites behind `sc_og_test.newick` from `sc_PolE_filtered_no_singletons.h5ad`),
  `K562_bulk_sbs_sites.R` (de novo SBS from bulk WGS of the
  mutation-accumulation clones) and `K562_sbs_rt_enrichment.R`
  (trinucleotide-matched enrichment in early/mid/late replicating thirds of
  the callable genome). `K562_sbs_tree_germline.R` flags germline/ancestral,
  chrM and unplaced-contig sites; `K562_sbs_tree_rebuild.R` rebuilds the
  single-cell and population SBS trees without them.

Stage them with:

```bash
bash paper_figures/scripts/sync_missing_inputs.sh
# or, to pull from another checkout that has already run the notebook:
K562_MUT_ACCUM_RESULTS=/path/to/SPC_genome/results/K562_mut_accumulation \
  bash paper_figures/scripts/sync_missing_inputs.sh
```

## External/large inputs
Some files were excluded from tracked data due GitHub file-size constraints.

- `paper_figures/data/encode_variant_annotation/somatic_max_peaks.tsv` (very large)
- `paper_figures/data/external/COSMIC_v3.4_SBS_GRCh38.txt` (required by `mutation_spectra.R`)

Place missing files at those exact paths before running all scripts.

## Notes
- Script behavior is preserved as much as possible; pathing and outputs were standardized.
- `bulk_VAF_bottleneck.R` supports `filtered_bulk_vaf.csv` or `filtered_bulk_vaf.csv.gz`.
- Known current blockers:
  - `plot_spc_sizes.R` expects non-empty `mask_summaries.csv` (source file currently empty).
  - `draw_trees.R` has unresolved legacy assumptions in downstream sections (`depth` field in derived objects).
- `mutation_spectra.R` still carries its own copies of the spectrum helpers;
  it can be switched over to `spectrum_utils.R` once its inputs are staged
  and it can be re-run end to end.
