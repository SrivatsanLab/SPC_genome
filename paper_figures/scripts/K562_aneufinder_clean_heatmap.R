# AneuFinder genome-wide heatmap of the K562 sc_PolE cells without the cells
# called artefacts, drawn bare for figure assembly: no cell names, axis text,
# legends or panel backgrounds. Dendrogram, population bar and chromosome
# boundary lines are kept; legends go to a separate *_legends.pdf.
#
# Input: an AneuFinder run (K562_AF_DIR, default the GC-corrected 1 Mb run).
# Cells dropped: the off-peak cells classed as artefacts by
# K562_offpeak_ploidy.R (per_cell_metrics.csv; noisy or off-integer depth, or
# depth matching the consensus with only the ploidy mis-scaled), identified on
# the 1 Mb run and dropped from every bin size.
#
# Same steps as AneuFinder::heatmapGenomewide(): clusterHMMs() order and
# dendrogram, segment states in AneuFinder's colours (stateColors()), genome
# coordinates from transCoord(). K562_DROP_NOISY=1 also drops cells flagged
# noisy (within-segment coverage CV above the threshold) by K562_cell_noise.R.
# Each cell is a tile (geom_rect) rather than a
# fixed-width line, so rows do not overlap at a normal page size.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
  library(ggplot2)
})

af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_gc/output"))
model_files <- list.files(file.path(af_dir, "MODELS/method-edivisive"), "\\.RData$", full.names = TRUE)
tag <- Sys.getenv("K562_OUT_TAG", sub(".*_binsize_([^_]+)_.*", "\\1", basename(model_files[1])))   # K562_OUT_TAG names the outputs
output_dir <- file.path(project_root, "paper_figures/output/K562_aneufinder_clean_heatmap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pop_colors <- c("P_0" = "#260091", "P_1" = "#1E90FF", "P_2" = "#FFDB58",
                "P_3" = "#FF9D71", "P_4" = "#FF1B5E", "P_5" = "#E3E6E6")
width_in <- 10
height_in <- 12

# ---- cells ------------------------------------------------------------------

metrics <- read.csv(file.path(project_root, "paper_figures/output/K562_offpeak_ploidy/per_cell_metrics.csv"))
artefact <- metrics$cell[grepl("^artefact", metrics$class)]
cell_of <- sub("(\\.bam)?_binsize.*", "", basename(model_files))
keep <- !(cell_of %in% artefact)
if (Sys.getenv("K562_DROP_NOISY") == "1") {
  nz <- read.csv(file.path(project_root, "paper_figures/output/K562_cell_noise/per_cell_noise.csv"))
  keep <- keep & !(cell_of %in% nz$cell[nz$noisy])
}
model_files <- model_files[keep]
meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
pop <- setNames(meta$pop, meta$cell)[cell_of[keep]]
stopifnot(!anyNA(pop))
cat(sprintf("%s: %d artefact%s cells dropped, %d kept\n", tag, sum(!keep),
            if (Sys.getenv("K562_DROP_NOISY") == "1") " or noisy" else "", sum(keep)))

# ---- order, segments, dendrogram (as in heatmapGenomewide) ------------------

hmms <- loadFromFiles(model_files, check.class = c("aneuHMM", "aneuBiHMM"))
ids <- sapply(hmms, `[[`, "ID")
names(pop) <- ids
cl <- clusterHMMs(hmms, cluster = TRUE)
hmms <- hmms[cl$IDorder]
ids <- sapply(hmms, `[[`, "ID")

segs <- lapply(hmms, function(h) AneuFinder:::transCoord(h$segments))
df <- do.call(rbind, lapply(seq_along(segs), function(i) {
  s <- segs[[i]]
  data.frame(x = i, start = s$start.genome, end = s$end.genome, state = as.character(s$state))
}))
lev <- names(sort(AneuFinder:::initializeStates(unique(df$state))$multiplicity))
df$state <- factor(df$state, levels = lev)
state_cols <- stateColors(levels(df$state))

sl <- seqlengths(segs[[1]])
chrom_lines <- data.frame(y = c(0, cumsum(as.numeric(sl))))
n <- length(hmms)

heat <- ggplot(df) +
  geom_rect(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = start, ymax = end, fill = state)) +
  geom_hline(data = chrom_lines, aes(yintercept = y), linewidth = 0.2, colour = "black") +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  coord_flip() + theme_void() + theme(legend.position = "none")

classbar <- ggplot(data.frame(x = seq_len(n), pop = pop[ids])) +
  geom_rect(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 0, ymax = 1, fill = pop)) +
  scale_fill_manual(values = pop_colors) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  coord_flip() + theme_void() + theme(legend.position = "none")

dd <- ggdendro::dendro_data(stats::as.dendrogram(cl$hclust), type = "rectangle")
dendro <- ggplot(dd$segments) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.15) +
  scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0), limits = c(0.5, n + 0.5)) +
  coord_flip() + theme_void()

# AneuFinder's relative widths: dendrogram 20, class bar 5, heatmap 150 * genome / 3 Gb
widths <- c(20, 5, sum(as.numeric(sl)) / 3e9 * 150)
fig <- cowplot::plot_grid(dendro, classbar, heat, nrow = 1, align = "h", rel_widths = widths)
out <- file.path(output_dir, paste0("genomeHeatmap_no_artefacts_", tag, ".pdf"))
ggsave(out, fig, width = width_in, height = height_in, bg = "transparent")

# ---- legends (separate file) ------------------------------------------------

# copy number 0-10 in AneuFinder's colours; every state above 10 is black there
leg_states <- paste0(0:11, "-somy")
leg_cols <- setNames(stateColors(leg_states), c(0:10, "11+"))
legend_of <- function(values, name) {
  d <- data.frame(k = factor(names(values), levels = names(values)))
  cowplot::get_plot_component(
    ggplot(d) + geom_tile(aes(1, 1, fill = k)) + scale_fill_manual(values = values, name = name) +
      theme(legend.position = "right"),
    "guide-box-right")
}
ggsave(file.path(output_dir, paste0("genomeHeatmap_no_artefacts_", tag, "_legends.pdf")),
       cowplot::plot_grid(legend_of(leg_cols, "Copy number"), legend_of(pop_colors, "Population"),
                          ncol = 1, rel_heights = c(1.3, 1)), width = 1.5, height = 6.5)

write.csv(data.frame(row = seq_len(n), cell = sub("\\.bam$", "", ids), pop = pop[ids]),
          file.path(output_dir, paste0("genomeHeatmap_no_artefacts_", tag, "_row_order.csv")),
          row.names = FALSE)
cat("Wrote", out, "\n")
