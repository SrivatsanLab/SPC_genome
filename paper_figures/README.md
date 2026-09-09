# Paper Figures Reproducibility Bundle

This directory consolidates figure-related scripts and input data from legacy `../code` and `../data` into a GitHub-friendly workflow.

## Structure
- `paper_figures/scripts/`: plotting and analysis scripts (ported from legacy code with repo-relative paths).
- `paper_figures/data/`: curated figure input data.
- `paper_figures/output/`: generated plots/tables.
- `paper_figures/run_all_figures.R`: orchestrator to run all scripts and write `output/run_summary.csv`.

## Run
From repository root:

```bash
Rscript paper_figures/run_all_figures.R
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
