# Worm6 co-assay read count scatter.
#
# Ported from notebooks/worm6_final_GEX.ipynb cell 16, which draws it with
# seaborn. Transcriptomic reads against genes detected, one point per cell,
# coloured by the matched WGS depth.
#
# The inputs are per-cell values already stored upstream, so unlike the UMAP
# nothing here is re-derived and the figure reproduces the published one
# exactly. Stage the CSV with export_worm6_coassay_counts.py, which joins
# GEX total_counts / genes_detected to mean_coverage from the DNA obs table,
# sorts by total_counts descending and drops the single highest cell, as the
# notebook does.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(scales)

data_dir <- file.path(project_root, "paper_figures/data/worm6_final/")
output_dir <- file.path(project_root, "paper_figures/output/worm6_coassay_scatter/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(data_dir, "coassay_readcount_counts.csv")

if (!file.exists(counts_path)) {
  message("Skipping worm6 co-assay scatter: ", counts_path, " not found.\n",
          "  Run export_worm6_coassay_counts.py first.")
} else {

  counts <- read.csv(counts_path, row.names = 1)

  # The notebook colours points on a continuous viridis normalised to [0, 3]
  # with clipping, then hand-builds a three-swatch legend at 0.1 / 1.5 / 3.0.
  # Passing guide_legend() to a continuous scale does the same thing: discrete
  # keys drawn at those break values, coloured by the scale itself.
  ggplot(counts, aes(x = total_counts, y = genes_detected, fill = mean_coverage)) +
    geom_point(shape = 21, size = 1.8, colour = "black", stroke = 0.25) +
    scale_fill_viridis_c(
      limits = c(0, 3),
      oob = scales::squish,
      breaks = c(0.1, 1.5, 3.0),
      labels = c("0–1x", "1–3x", ">3x"),
      name = "WGS Depth",
      guide = guide_legend(override.aes = list(size = 3))
    ) +
    # matplotlib's default break spacing, which ggplot would otherwise put at 3000
    scale_y_continuous(breaks = seq(0, 10000, by = 2000)) +
    scale_x_continuous(breaks = seq(0, 20000, by = 10000)) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    xlab("Transcriptomic Reads") +
    ylab("Genes Detected")
  ggsave(file.path(output_dir, "coassay_readcount_scatter.svg"),
         bg = "transparent",
         height = 4, width = 6)

  cat("Wrote scatter to:", output_dir, "\n")
}
