# Per-cell agreement of K562 sc_PolE copy number with the K562 reference
# profile, reproducing notebooks/K562_tree.ipynb cells 135-148 on the
# GC-corrected AneuFinder run (K562_AF_DIR, default the 1 Mb run).
#
# Unmappable-bin mask: bins >= 50% N (hg38 gaps) or ENCODE blacklist, the
# regions AneuFinder has no reads for after blacklisting (they read CN0 in
# every cell). Written to results/K562_tree/masks/unmappable_bins_<binsize>.csv,
# which the notebook also reads. Defined from the genome, not the calls, so the
# real 9p21 (CDKN2A) homozygous deletion is kept.
#
# Bin-to-reference mappings:
#   notebook        - every bin joined to each reference segment it overlaps (a
#                     bin spanning a breakpoint appears once per segment), bins
#                     with no reference segment filled with 3, as in the notebook
#   notebook_masked - the same after dropping masked bins and chrY (the
#                     notebook's current method; used for the figure)
#   strict          - each bin takes the segment at its midpoint; masked,
#                     uncovered and chrY bins dropped
# Per cell: Pearson r, cosine distance (1 - cosine similarity of the raw
# copy-number vectors) and mean |CN - reference|. Significance: the notebook's
# null, the mean r (and mean cosine distance) across cells after shuffling the
# reference over bins, 1000 permutations (seed 42); empirical p = fraction of
# null means at least as extreme as observed. Also the mean pairwise cosine
# distance between cells (each bin counted once). Reported for all 1000 cells
# and for the 890 at the dominant ploidy [2.5, 3.0).

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_gc/output"))
n_perm <- 1000

cnv <- read.csv(file.path(af_dir, "result.csv"), row.names = 1, check.names = FALSE)
# 1 Mb outputs keep the original folder; other bin sizes get a suffix
bin_tag <- format(cnv$end[1] - cnv$start[1] + 1, scientific = TRUE, trim = TRUE)
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_cnv_reference_correlation",
                                             if (bin_tag == "1e+06") "" else paste0("_", bin_tag), "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cells <- colnames(cnv)[-(1:3)]
cn <- as.matrix(cnv[, cells])
ploidy <- colMeans(cn)
on_peak <- ploidy >= 2.5 & ploidy < 3.0
ref <- read.table(file.path(sc_dir, "K562_ref_ploidy_hg38.bed"),
                  col.names = c("chr", "start", "end", "cn"))

# ---- unmappable-bin mask ---------------------------------------------------

binsize <- cnv$end[1] - cnv$start[1] + 1
mask_csv <- file.path(project_root, "results/K562_tree/masks",
                      paste0("unmappable_bins_", format(binsize, scientific = TRUE, trim = TRUE), ".csv"))
gr <- GRanges(cnv$seqnames, IRanges(cnv$start, cnv$end))
n_frac <- letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, gr), "N", as.prob = TRUE)[, 1]
bl <- read.table(file.path(af_dir, "blacklist-hg38.raw.bed.gz"))
bl <- reduce(GRanges(bl$V1, IRanges(bl$V2 + 1, bl$V3)))
ov <- findOverlaps(gr, bl)
bl_bp <- tapply(width(pintersect(gr[queryHits(ov)], bl[subjectHits(ov)])), queryHits(ov), sum)
bl_frac <- numeric(length(gr))
bl_frac[as.integer(names(bl_bp))] <- bl_bp / width(gr)[as.integer(names(bl_bp))]
masked <- n_frac + bl_frac >= 0.5
write.csv(data.frame(seqnames = cnv$seqnames, start = cnv$start, end = cnv$end,
                     n_frac = n_frac, blacklist_frac = bl_frac, masked = masked),
          mask_csv, row.names = FALSE)
cat(sprintf("%d of %d bins masked (>= 50%% N or blacklist) -> %s\n", sum(masked), length(masked), mask_csv))

# ---- bin-to-reference mappings ---------------------------------------------

# notebook: overlap join, one row per (bin, overlapping segment), fill 3
nb <- do.call(rbind, lapply(seq_len(nrow(cnv)), function(i) {
  h <- which(ref$chr == cnv$seqnames[i] & ref$start < cnv$end[i] & ref$end > cnv$start[i] - 1)
  if (length(h)) data.frame(bin = i, ref = ref$cn[h], filled = FALSE)
  else data.frame(bin = i, ref = 3, filled = TRUE)
}))
# strict: midpoint, uncovered and chrY dropped
mid <- (cnv$start + cnv$end) %/% 2
st_ref <- sapply(seq_len(nrow(cnv)), function(i) {
  h <- which(ref$chr == cnv$seqnames[i] & ref$start < mid[i] & ref$end >= mid[i])
  if (length(h)) ref$cn[h[1]] else NA
})
st <- data.frame(bin = which(!is.na(st_ref) & cnv$seqnames != "chrY" & !masked))
st$ref <- st_ref[st$bin]
nbm <- nb[!masked[nb$bin] & cnv$seqnames[nb$bin] != "chrY", ]

mappings <- list(notebook = nb, notebook_masked = nbm, strict = st)
cat(sprintf("notebook: %d rows (%d filled with 3); notebook_masked: %d rows (%d filled); strict: %d bins\n",
            nrow(nb), sum(nb$filled), nrow(nbm), sum(nbm$filled), nrow(st)))

# ---- per-cell scores and permutation null ----------------------------------

set.seed(42)
per_cell <- list(); summary_rows <- list(); null_rows <- list()
for (m in names(mappings)) {
  x <- cn[mappings[[m]]$bin, , drop = FALSE]
  r <- mappings[[m]]$ref
  rho <- suppressWarnings(as.vector(cor(x, r)))
  cos_dist <- function(xx, rr) 1 - as.vector(crossprod(xx, rr)) / (sqrt(colSums(xx^2)) * sqrt(sum(rr^2)))
  cosd <- cos_dist(x, r)
  mad <- colMeans(abs(x - r))
  per_cell[[m]] <- data.frame(cell = cells, mapping = m, pearson = rho, cosine_distance = cosd,
                              mean_abs_diff = mad, ploidy = ploidy, on_peak = on_peak)
  # cell-cell cosine distance, each bin once
  xu <- cn[unique(mappings[[m]]$bin), , drop = FALSE]
  xn <- sweep(xu, 2, sqrt(colSums(xu^2)), "/")
  cell_cos <- 1 - crossprod(xn)
  for (set in c("all", "on_peak")) {
    keep <- if (set == "all") rep(TRUE, length(cells)) else on_peak
    obs <- mean(rho[keep], na.rm = TRUE)
    obs_cos <- mean(cosd[keep])
    perms <- replicate(n_perm, {
      rs <- sample(r)
      c(mean(suppressWarnings(cor(x[, keep], rs)), na.rm = TRUE), mean(cos_dist(x[, keep], rs)))
    })
    null <- perms[1, ]
    pw <- cell_cos[keep, keep][upper.tri(diag(sum(keep)))]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      mapping = m, cells = set, n_cells = sum(keep), n_bins = length(r),
      mean_r = obs, sd_r = sd(rho[keep], na.rm = TRUE), median_r = median(rho[keep], na.rm = TRUE),
      mean_abs_diff = mean(mad[keep]), null_mean = mean(null), null_max = max(null),
      empirical_p = mean(null >= obs),
      mean_cos_dist_to_ref = obs_cos, sd_cos_dist_to_ref = sd(cosd[keep]),
      null_mean_cos_dist = mean(perms[2, ]), null_min_cos_dist = min(perms[2, ]),
      empirical_p_cos = mean(perms[2, ] <= obs_cos),
      mean_pairwise_cos_dist = mean(pw), sd_pairwise_cos_dist = sd(pw),
      median_pairwise_cos_dist = median(pw))
    null_rows[[length(null_rows) + 1]] <- data.frame(mapping = m, cells = set, null_mean_r = null)
  }
}
per_cell <- do.call(rbind, per_cell)
summ <- do.call(rbind, summary_rows)
nulls <- do.call(rbind, null_rows)
write.csv(per_cell, file.path(output_dir, "per_cell_reference_agreement.csv"), row.names = FALSE)
write.csv(summ, file.path(output_dir, "reference_correlation_summary.csv"), row.names = FALSE)
print(summ, digits = 3)

# ---- figure (notebook_masked, all cells; legend in a separate file) ---------

pc <- per_cell %>% filter(mapping == "notebook_masked")
s <- summ %>% filter(mapping == "notebook_masked", cells == "all")

p_scatter <- ggplot(pc, aes(mean_abs_diff, pearson, color = ploidy)) +
  geom_point(size = 0.9, alpha = 0.8) +
  scale_color_viridis_c(name = "Average\nploidy") +
  labs(x = "Mean absolute copy-number difference", y = "Pearson r with K562 reference") +
  theme_classic(base_size = 8)

nl <- nulls %>% filter(mapping == "notebook_masked", cells == "all")
p_hist <- ggplot(pc, aes(pearson)) +
  geom_histogram(binwidth = 0.02, fill = "grey55", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = s$mean_r, linewidth = 0.4) +
  geom_vline(xintercept = max(nl$null_mean_r), linetype = "dashed", linewidth = 0.4, color = "grey40") +
  annotate("text", x = s$mean_r, y = Inf, vjust = 1.5, hjust = 1.05, size = 2.5,
           label = sprintf("mean r = %.2f (s.d. %.2f)", s$mean_r, s$sd_r)) +
  annotate("text", x = max(nl$null_mean_r), y = Inf, vjust = 4, hjust = -0.05, size = 2.5,
           color = "grey30",
           label = sprintf("shuffled max = %.2f\np < %.3f", max(nl$null_mean_r), 1 / n_perm)) +
  labs(x = "Pearson r with K562 reference", y = "Cells") +
  theme_classic(base_size = 8)

fig <- p_hist + (p_scatter + theme(legend.position = "none")) + plot_layout(widths = c(1, 1))
# the ploidy colour bar on its own (invisible points keep the scale), to be
# placed by hand
legend_plot <- ggplot(pc, aes(mean_abs_diff, pearson, color = ploidy)) +
  geom_point(alpha = 0) + scale_color_viridis_c(name = "Average\nploidy") +
  theme_void(base_size = 8) +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5))
ggsave(file.path(output_dir, "reference_correlation_legend.svg"), legend_plot, width = 1, height = 1.5)
ggsave(file.path(output_dir, "reference_correlation.svg"), fig, width = 7, height = 3)
ggsave(file.path(output_dir, "reference_correlation.png"), fig, width = 7, height = 3, dpi = 300,
       bg = "white")
cat("Outputs written to", output_dir, "\n")
