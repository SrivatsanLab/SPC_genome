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

- `K562_mut_accumulation.R` <- `notebooks/K562_mut_accumulation.ipynb`
  (`spectrum.csv`, `spectrum_background.csv`). Renders the de novo spectra
  pooled by construct (AAVS/PolE), the per-sample facet grid, and the
  ancestral background spectrum.

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
