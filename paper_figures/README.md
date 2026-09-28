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
