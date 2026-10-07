# Genome-wide copy-number heatmap of the 1000 sc_PolE K562 cells with genomic
# bins as rows (chr1 at the top) and cells as columns.
#
# Input: a GC-corrected AneuFinder run (scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh
# [binsize]); K562_AF_DIR picks the run, default the 1 Mb one. chr1-22 and X.
# Rows stay in genome order; columns are ordered by hierarchical clustering of
# all cells (Ward.D2 on Manhattan distance of copy number). Top strips: population, ploidy group (mean CN; [2.5, 3.0) is
# the dominant near-triploid peak) and Pearson r with the K562 reference
# profile. Copy number on a diverging scale centred on 3, capped at 8.
#
# Legends are written to a separate *_legends.svg.
# Drawn with base graphics (image(useRaster = TRUE)) so 250 kb bins stay fast.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

sc_dir <- Sys.getenv("K562_SC_PROJECT",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_gc/output"))

pop_colors <- c("P_0" = "#260091", "P_1" = "#1E90FF", "P_2" = "#FFDB58",
                "P_3" = "#FF9D71", "P_4" = "#FF1B5E", "P_5" = "#E3E6E6")
ploidy_colors <- c("<2.5" = "#5b8fd6", "2.5-3.0" = "#d9d9d9", ">3.0" = "#d6604d")
# same scale as K562_cnv_heatmap.R: these colours at copy number 0, 2, 3, 4, 6, 8
cn_stops <- c(0, 2, 3, 4, 6, 8)
cn_stop_colors <- c("#08306b", "#6baed6", "#f0f0f0", "#fc9272", "#cb181d", "#67000d")
cn_palette <- sapply(0:8, function(v) {
  i <- min(findInterval(v, cn_stops), length(cn_stops) - 1)
  f <- (v - cn_stops[i]) / (cn_stops[i + 1] - cn_stops[i])
  rgb(colorRamp(cn_stop_colors[i:(i + 1)])(f), maxColorValue = 255)
})

# ---- data -------------------------------------------------------------------

cnv <- read.csv(file.path(af_dir, "result.csv"), row.names = 1, check.names = FALSE)
binsize <- cnv$end[1] - cnv$start[1] + 1
tag <- format(binsize, scientific = TRUE, trim = TRUE)
output_dir <- file.path(project_root, "paper_figures/output/K562_cnv_heatmap_by_bin/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

keep <- cnv$seqnames %in% c(paste0("chr", 1:22), "chrX")
bins <- cnv[keep, c("seqnames", "start", "end")]
cells <- colnames(cnv)[-(1:3)]
cn <- as.matrix(cnv[keep, cells])

meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
pop <- setNames(meta$pop, meta$cell)[cells]
ploidy <- colMeans(as.matrix(cnv[, cells]))
ploidy_group <- as.character(cut(ploidy, c(0, 2.5, 3.0, 99), labels = names(ploidy_colors),
                                 right = FALSE))

# Pearson r with the reference: each bin takes the reference segment at its
# midpoint; bins the reference does not cover, and unmappable bins, are dropped
ref <- read.table(file.path(sc_dir, "K562_ref_ploidy_hg38.bed"),
                  col.names = c("chr", "start", "end", "cn"))
mid <- (bins$start + bins$end) %/% 2
ref_cn <- sapply(seq_len(nrow(bins)), function(i) {
  h <- which(ref$chr == bins$seqnames[i] & ref$start < mid[i] & ref$end >= mid[i])
  if (length(h)) ref$cn[h[1]] else NA
})
covered <- !is.na(ref_cn)
# drop unmappable bins if the mask for this bin size exists (K562_cnv_reference_correlation.R)
mask_csv <- file.path(project_root, "results/K562_tree/masks", paste0("unmappable_bins_", tag, ".csv"))
if (file.exists(mask_csv)) {
  mk <- read.csv(mask_csv)
  covered <- covered & !mk$masked[match(paste(bins$seqnames, bins$start), paste(mk$seqnames, mk$start))]
}
ref_r <- apply(cn[covered, ], 2, function(v) if (sd(v) > 0) cor(v, ref_cn[covered]) else NA)

# ---- column order -----------------------------------------------------------

hc <- hclust(dist(t(cn), method = "manhattan"), method = "ward.D2")
ord <- hc$order

# ---- draw -------------------------------------------------------------------

chr_runs <- rle(bins$seqnames)
chr_end <- cumsum(chr_runs$lengths)
chr_start <- chr_end - chr_runs$lengths + 1
n_bin <- nrow(bins)
n_cell <- length(ord)

draw <- function() {
  layout(matrix(c(0, 1,
                  0, 2,
                  0, 3,
                  4, 5), ncol = 2, byrow = TRUE),
         widths = c(0.9, 10), heights = c(0.25, 0.25, 0.25, 10))
  strip <- function(cols, label) {
    par(mar = c(0.1, 0.5, 0.1, 0.5))
    image(seq_len(n_cell), 1, matrix(seq_len(n_cell), ncol = 1), col = cols,
          axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
    mtext(label, side = 2, las = 1, cex = 0.55, line = 0.2, xpd = NA)
  }
  strip(pop_colors[pop[ord]], "Population")
  strip(ploidy_colors[ploidy_group[ord]], "Ploidy")
  grey_ramp <- grey(seq(1, 0, length.out = 101))
  strip(grey_ramp[1 + round(100 * pmax(0, ref_r[ord]))], "Reference r")

  # chromosome labels, alternating shading
  par(mar = c(0.5, 0.1, 0.1, 0))
  plot(NA, xlim = c(0, 1), ylim = c(n_bin + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(0.55, chr_start - 0.5, 1, chr_end + 0.5,
       col = rep(c("grey85", "white"), length.out = length(chr_start)), border = NA)
  text(0.5, (chr_start + chr_end) / 2, sub("chr", "", chr_runs$values), adj = 1, cex = 0.5)

  # heatmap: rows = bins (top to bottom), columns = cells
  par(mar = c(0.5, 0.5, 0.1, 0.5))
  z <- t(pmin(cn[, ord], 8))[, n_bin:1]
  image(seq_len(n_cell), seq_len(n_bin), z, col = cn_palette, breaks = seq(-0.5, 8.5, 1),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  abline(h = n_bin - chr_start[-1] + 1.5, col = "black", lwd = 0.3)
  mtext(sprintf("%d cells, hierarchically clustered; %s bp bins", n_cell, tag), side = 1,
        line = -0.2, cex = 0.55)
}

# legends go in their own file, to be placed by hand
draw_legends <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot.new()
  y <- 1
  for (l in list(list("Copy number", c(0:7, "8+"), cn_palette, NA),
                 list("Population", names(pop_colors), pop_colors, NA),
                 list("Ploidy", names(ploidy_colors), ploidy_colors, NA),
                 list("Reference r", c("0", "0.5", "1"), grey(c(1, 0.5, 0)), "grey60"))) {
    lg <- legend(0, y, legend = l[[2]], fill = l[[3]], border = l[[4]], title = l[[1]],
                 title.adj = 0, bty = "n", cex = 0.6, xpd = NA)
    y <- y - lg$rect$h - 0.02
  }
}

base <- file.path(output_dir, paste0("cnv_heatmap_by_bin_", tag))
svg(paste0(base, "_legends.svg"), width = 1.5, height = 5, bg = "white")
draw_legends(); invisible(dev.off())
png(paste0(base, ".png"), width = 8, height = 11, units = "in", res = 300, bg = "white",
    type = "cairo")
draw(); invisible(dev.off())
svg(paste0(base, ".svg"), width = 8, height = 11, bg = "white")
draw(); invisible(dev.off())
pdf(paste0(base, ".pdf"), width = 8, height = 11)
draw(); invisible(dev.off())

write.csv(data.frame(cell = cells, pop = pop, ploidy = ploidy, ploidy_group = ploidy_group,
                     reference_r = ref_r, column = match(seq_along(cells), ord)),
          paste0(base, "_cells.csv"), row.names = FALSE)
cat("Wrote", paste0(base, ".svg/.png/.pdf"), "and", paste0(base, "_legends.svg"), "\n")
