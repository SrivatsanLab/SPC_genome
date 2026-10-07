# Cell-cycle phase of the K562 sc_PolE cells from how their coverage follows
# replication timing, with each chromosome as an independent reporter.
#
# Replication timing: K562 UW Repli-seq on hg38
# (scripts/utils/download_replication_timing.sh, replication_timing_to_hg38.R),
# RT = early-replication score of 1 kb windows (1 = G1b/earliest, 0 = G2),
# averaged over each kept 250 kb bin.
#
# Coverage: fragment-level GC-weighted unique fragments per kept 250 kb bin
# (fragment_gc_counts.R / fragment_gc_segment.R: repeat, unmappable and
# low-mappability bins dropped; depth / mappable fraction). Fragment-level GC
# weights act on each fragment's own GC, so the regional replication signal is
# not regressed away as it would be by bin-level GC correction. Divided by the
# consensus copy number (median integer CN of on-peak cells, GC-corrected 1 Mb
# run) and by the cell's median, on a log2 scale: a G1 or G2 cell is flat
# (copy number aside), an S-phase cell is raised where replication is done.
#
# Per chromosome (chr1-22, X; >= 40 bins with RT):
#   slope  - least-squares slope of log2 coverage on RT
#   step   - best two-level fit y = a + b * 1[RT > t], b >= 0, over thresholds t
#            at RT percentiles 5-95; replicated fraction = share of the
#            chromosome's bins with RT > t
# Per cell: median slope across chromosomes; number of chromosomes with
# slope > 0 and a two-sided sign test.
#
# Two passes. Against consensus copy number alone, non-replicating cells are
# not flat: they share a negative slope (early, GC-rich, open regions are
# under-covered in amplification), forming a tight peak below 0. So, as in
# scRepli-seq / Kronos, the core of that peak (median slope < ref_max_slope and
# <= 3 chromosomes with positive slope) is taken as the non-replicating
# reference, and every cell is re-measured against the reference cells' median
# profile (which also carries copy number). In this second pass S phase = sign
# test p < 0.01 and median slope above the reference cells' 99th percentile
# (the sign test alone also passes small consistent offsets of G1 cells); G1
# or G2 otherwise. For S-phase cells,
# the replicated fraction is the median across chromosomes with b > 0.
# G1 and G2 cannot be told apart: per-cell normalisation removes total DNA.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})

rt_bed <- Sys.getenv("K562_RT_BED", file.path(project_root,
                     "data/reference/replication_timing/K562_replication_timing_hg38.bed.gz"))
run_dir <- file.path(project_root, "results/K562_tree/aneufinder_fragment_gc/output")
output_dir <- file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
min_bins <- 40
s_phase_p <- 0.01
ref_max_slope <- -0.4

# ---- bins, replication timing, consensus copy number --------------------------

kb <- read.csv(file.path(run_dir, "kept_bins.csv"))
use <- kb$kept & kb$seqnames != "chrY"
bins <- GRanges(kb$seqnames[use], IRanges(kb$start[use], kb$end[use]))

rt <- read.table(rt_bed, col.names = c("chrom", "start", "end", "G1", "S1", "S2", "S3", "S4", "G2", "RT"))
rt_gr <- GRanges(rt$chrom, IRanges(rt$start + 1, rt$end), RT = rt$RT)
ov <- findOverlaps(bins, rt_gr)
w <- width(pintersect(bins[queryHits(ov)], rt_gr[subjectHits(ov)]))
bin_rt <- rep(NA_real_, length(bins))
s <- tapply(rt_gr$RT[subjectHits(ov)] * w, queryHits(ov), sum) / tapply(w, queryHits(ov), sum)
bin_rt[as.integer(names(s))] <- s

cnv1 <- read.csv(file.path(project_root, "results/K562_tree/aneufinder_gc/output/result.csv"),
                 row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
cons_gr <- GRanges(cnv1$seqnames, IRanges(cnv1$start, cnv1$end),
                   cn = apply(as.matrix(cnv1[, meta$cell[meta$kept]]), 1, median))
bin_cn <- cons_gr$cn[findOverlaps(resize(bins, 1, fix = "center"), cons_gr, select = "first")]

ok <- is.finite(bin_rt) & !is.na(bin_cn) & bin_cn > 0
chrom <- as.character(seqnames(bins))[ok]
rt_ok <- bin_rt[ok]
cat(sprintf("%d of %d kept bins with replication timing and consensus CN > 0\n", sum(ok), length(ok)))
write.csv(data.frame(seqnames = chrom, start = start(bins)[ok], end = end(bins)[ok], rt = rt_ok,
                     consensus_cn = bin_cn[ok], gc = kb$GC[use][ok]),
          file.path(output_dir, "bin_replication_timing.csv"), row.names = FALSE)

# ---- coverage -------------------------------------------------------------------

fc <- readRDS(file.path(run_dir, "fragment_counts.rds"))
cells <- colnames(fc$counts_gc)
depth <- fc$counts_gc[use, ][ok, ] / kb$mappable_frac[use][ok] / bin_cn[ok]
y_all <- log2(sweep(depth, 2, apply(depth, 2, median), "/") + 1e-3)

# ---- per chromosome, per cell ----------------------------------------------------

chroms <- intersect(c(paste0("chr", 1:22), "chrX"), names(which(table(chrom) >= min_bins)))
step_fit <- function(y, x) {
  ts <- quantile(x, seq(0.05, 0.95, 0.025))
  best <- list(rss = sum((y - mean(y))^2), b = 0, t = NA, frac = NA)
  for (t in ts) {
    hi <- x > t
    b <- mean(y[hi]) - mean(y[!hi])
    if (b <= 0) next
    rss <- sum((y[hi] - mean(y[hi]))^2) + sum((y[!hi] - mean(y[!hi]))^2)
    if (rss < best$rss) best <- list(rss = rss, b = b, t = t, frac = mean(hi))
  }
  best
}

chromosome_fits <- function(ymat) do.call(rbind, lapply(seq_along(cells), function(j) {
  y <- ymat[, j]
  do.call(rbind, lapply(chroms, function(ch) {
    i <- chrom == ch
    st <- step_fit(y[i], rt_ok[i])
    data.frame(cell = cells[j], chrom = ch, slope = unname(coef(lm(y[i] ~ rt_ok[i]))[2]),
               step_b = st$b, replicated_frac = st$frac)
  }))
}))

# pass 1: against consensus copy number; pick the non-replicating reference
pass1 <- chromosome_fits(y_all)
p1_cell <- do.call(rbind, lapply(split(pass1, pass1$cell), function(d)
  data.frame(cell = d$cell[1], median_slope = median(d$slope), n_positive = sum(d$slope > 0))))
ref_cells <- p1_cell$cell[p1_cell$median_slope < ref_max_slope & p1_cell$n_positive <= 3]
cat(sprintf("pass 1: %d reference (non-replicating) cells, median slope < %.2f and <= 3 chromosomes positive\n",
            length(ref_cells), ref_max_slope))
write.csv(p1_cell, file.path(output_dir, "pass1_against_copy_number.csv"), row.names = FALSE)

# pass 2: against the reference cells' median profile
ref_profile <- apply(y_all[, ref_cells], 1, median)
y_ref <- y_all - ref_profile
per_chr <- chromosome_fits(y_ref)
write.csv(per_chr, file.path(output_dir, "per_chromosome_replication.csv"), row.names = FALSE)

per_cell <- do.call(rbind, lapply(split(per_chr, per_chr$cell), function(d) {
  n_pos <- sum(d$slope > 0); n <- nrow(d)
  frac <- d$replicated_frac[d$step_b > 0]
  data.frame(cell = d$cell[1], n_chrom = n, n_positive = n_pos, median_slope = median(d$slope),
             iqr_slope = IQR(d$slope), sign_p = binom.test(n_pos, n)$p.value,
             replicated_frac = if (length(frac)) median(frac) else NA,
             replicated_frac_iqr = if (length(frac) > 1) IQR(frac) else NA)
}))
slope_cut <- quantile(per_cell$median_slope[per_cell$cell %in% ref_cells], 0.99)
cat(sprintf("S-phase slope threshold (99th percentile of reference cells): %.3f\n", slope_cut))
per_cell$phase <- ifelse(per_cell$sign_p < s_phase_p & per_cell$median_slope > slope_cut, "S", "G1/G2")
per_cell$replicated_frac[per_cell$phase != "S"] <- NA
per_cell$replicated_frac_iqr[per_cell$phase != "S"] <- NA
per_cell$s_stage <- cut(per_cell$replicated_frac, c(0, 1 / 3, 2 / 3, 1), labels = c("early S", "mid S", "late S"),
                        include.lowest = TRUE)

metrics <- read.csv(file.path(project_root, "paper_figures/output/K562_offpeak_ploidy/per_cell_metrics.csv"))
gcb <- read.csv(file.path(project_root, "paper_figures/output/K562_gc_toggle_ploidy/ploidy_with_without_gc.csv"))
per_cell$reference_cell <- per_cell$cell %in% ref_cells
per_cell$pop <- setNames(meta$pop, meta$cell)[per_cell$cell]
per_cell$class <- setNames(metrics$class, metrics$cell)[per_cell$cell]
per_cell$spikiness <- setNames(metrics$spikiness, metrics$cell)[per_cell$cell]
per_cell$gc_bias <- setNames(gcb$gc_bias, gcb$cell)[per_cell$cell]
per_cell$n_fragments <- setNames(fc$n_fragments, cells)[per_cell$cell]
write.csv(per_cell, file.path(output_dir, "per_cell_phase.csv"), row.names = FALSE)

# ---- summaries ---------------------------------------------------------------------

cat(sprintf("\n%d cells, %d chromosomes each\n", nrow(per_cell), length(chroms)))
cat("phase calls:\n"); print(table(per_cell$phase, useNA = "ifany"))
cat("S-phase stage:\n"); print(table(per_cell$s_stage))
cat("\nmedian slope quantiles 5/25/50/75/95%:",
    round(quantile(per_cell$median_slope, c(.05, .25, .5, .75, .95)), 3), "\n")
cat("chromosomes with positive slope (of", length(chroms), "), histogram:\n")
print(table(cut(per_cell$n_positive, c(-1, 3, 7, 11, 15, 19, 23))))
sp <- function(a, b) round(cor(a, b, method = "spearman", use = "complete"), 3)
cat(sprintf("\nSpearman of median slope with: GC-bias score %.3f, spikiness %.3f, fragments %.3f\n",
            sp(per_cell$median_slope, per_cell$gc_bias), sp(per_cell$median_slope, per_cell$spikiness),
            sp(per_cell$median_slope, per_cell$n_fragments)))
cat("\nS-phase fraction by population:\n")
tabp <- table(per_cell$pop, per_cell$phase); print(cbind(tabp, frac_S = round(tabp[, "S"] / rowSums(tabp), 3)))
cat(sprintf("chi-square across populations p = %.3g\n", chisq.test(tabp)$p.value))
cat("\nS-phase fraction by cell class:\n")
tabc <- table(per_cell$class, per_cell$phase); print(cbind(tabc, frac_S = round(tabc[, "S"] / rowSums(tabc), 3)))

# ---- figures (no legends inside; text kept minimal) ----------------------------------

p1 <- ggplot(per_cell, aes(median_slope, fill = phase)) +
  geom_histogram(binwidth = 0.025, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = c("G1/G2" = "grey60", "S" = "#2a6fb0")) +
  labs(x = "Median per-chromosome slope of coverage on replication timing", y = "Cells") +
  theme_classic(base_size = 8) + theme(legend.position = "none")
p2 <- ggplot(subset(per_cell, phase == "S"), aes(replicated_frac)) +
  geom_histogram(binwidth = 0.05, fill = "#2a6fb0", colour = "white", linewidth = 0.15) +
  labs(x = "Fraction of genome replicated (S-phase cells)", y = "Cells") +
  theme_classic(base_size = 8)
ex <- per_cell$cell[order(per_cell$median_slope)][c(50, round(nrow(per_cell) / 2), nrow(per_cell) - 50)]
ex_df <- do.call(rbind, lapply(ex, function(cl) data.frame(cell = cl, rt = rt_ok, y = y_ref[, cl])))
ex_df$cell <- factor(ex_df$cell, levels = ex, labels = c("low", "median", "high"))
p3 <- ggplot(ex_df, aes(rt, y)) + geom_point(size = 0.2, alpha = 0.25) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, colour = "#2a6fb0") +
  facet_wrap(~cell, nrow = 1) + coord_cartesian(ylim = c(-1.5, 1.5)) +
  labs(x = "Replication timing (1 = early)", y = "log2 coverage vs non-replicating reference") +
  theme_classic(base_size = 8)
fig <- (p1 + p2) / p3
ggsave(file.path(output_dir, "cell_cycle_from_coverage.svg"), fig, width = 7, height = 5)
ggsave(file.path(output_dir, "cell_cycle_from_coverage.png"), fig, width = 7, height = 5, dpi = 200, bg = "white")
cat("Outputs written to", output_dir, "\n")
