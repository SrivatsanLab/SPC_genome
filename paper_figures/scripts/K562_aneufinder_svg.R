# Re-draw AneuFinder's own plots as PDF (K562_FORMAT=svg for SVG, which gets
# very large: one element per segment) for one run (K562_AF_DIR, default the
# GC-corrected 1 Mb run), with AneuFinder's plotting functions and sizes:
#   genomeHeatmap      - heatmapGenomewide() of all cells, clustered (as in
#                        PLOTS/method-edivisive/genomeHeatmap_*.pdf)
#   aneuploidyHeatmap  - heatmapAneuploidies(), clustered, 30 in wide
#   qualityHeatmap     - heatmapGenomewide() of the best clusterByQuality()
#                        class (as in heatmapGenomewide.pdf from
#                        scripts/utils/run_aneufinder_K562_tree.R)
# Both genome-wide heatmaps carry a class bar with each cell's population
# (P_0-P_5, colours as in K562_cnv_heatmap.R); it does not change the
# clustering. K562_KARYOGRAM_CELLS (comma-separated cell barcodes) adds one karyogram SVG
# per listed cell.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(AneuFinder)
  library(ggplot2)
})

af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_gc/output"))
model_files <- list.files(file.path(af_dir, "MODELS/method-edivisive"), "\\.RData$", full.names = TRUE)
tag <- Sys.getenv("K562_OUT_TAG", sub(".*_binsize_([^_]+)_.*", "\\1", basename(model_files[1])))   # K562_OUT_TAG names the outputs
fmt <- Sys.getenv("K562_FORMAT", "pdf")
output_dir <- file.path(project_root, "paper_figures/output/K562_aneufinder_svg", tag)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("%d models, bin size %s\n", length(model_files), tag))

pop_colors <- c("P_0" = "#260091", "P_1" = "#1E90FF", "P_2" = "#FFDB58",
                "P_3" = "#FF9D71", "P_4" = "#FF1B5E", "P_5" = "#E3E6E6")
meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
pop <- setNames(meta$pop, meta$cell)[sub("(\\.bam)?_binsize.*", "", basename(model_files))]
stopifnot(!anyNA(pop))

suppressMessages(heatmapGenomewide(model_files, classes = pop,
                                   classes.color = pop_colors[sort(unique(pop))],
                                   cluster = TRUE, plot.breakpoints = FALSE,
                                   file = file.path(output_dir, paste0("genomeHeatmap_", tag, ".", fmt))))

p <- suppressMessages(heatmapAneuploidies(model_files, cluster = TRUE))
ggsave(file.path(output_dir, paste0("aneuploidyHeatmap_", tag, ".", fmt)), p,
       width = 30, height = max(0.3 * length(model_files), 2 / 2.54), limitsize = FALSE)

models <- loadFromFiles(model_files)
cl <- clusterByQuality(models)
best <- cl$classification[[1]]
# clusterByQuality() returns model file paths here, or models when given models
best_ids <- if (is.character(best)) sub("(\\.bam)?_binsize.*", "", basename(best)) else
  sub("\\.bam$", "", sapply(best, function(m) m$ID))
best_pop <- setNames(meta$pop, meta$cell)[best_ids]
suppressMessages(heatmapGenomewide(best, classes = best_pop,
                                   classes.color = pop_colors[sort(unique(best_pop))], cluster = TRUE,
                                   file = file.path(output_dir, paste0("qualityHeatmap_", tag, ".", fmt))))

karyo_cells <- strsplit(Sys.getenv("K562_KARYOGRAM_CELLS", ""), ",")[[1]]
if (length(karyo_cells)) {
  ids <- sub("\\.bam$", "", sapply(models, function(m) m$ID))
  for (cell in karyo_cells) {
    i <- match(cell, ids)
    if (is.na(i)) { warning("no model for ", cell); next }
    ggsave(file.path(output_dir, paste0("karyogram_", cell, ".", fmt)),
           plot(models[[i]], type = "karyogram", plot.breakpoints = TRUE),
           width = 12 * 1.4, height = 2 * 4.6)
  }
}
cat("Wrote", output_dir, "\n")
