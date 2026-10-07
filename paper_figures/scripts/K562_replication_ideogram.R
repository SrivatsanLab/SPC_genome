# One S-phase K562 cell on an ideogram: where its coverage is raised by
# replication versus where K562 replication timing says it should be, plus
# zoomed coverage tracks. Drawn with ggideogram (github.com/dxsbiocc/ggideogram).
#
# Cell: K562_CELL, default the best-fitting (replication_fit.csv r2) on-peak
# cell with inferred DNA content closest to K562_DNA_CONTENT (1.25), i.e. about
# a quarter of the genome replicated.
#
# Coverage y (log2, per 250 kb kept bin) as in
# scripts/utils/replication_regress_segment.R: fragment-level GC-weighted depth
# / mappable fraction / consensus copy number / cell median, minus the
# non-replicating reference profile. The cell's fit gives
# y = c + log2(1 + a * p(RT)), p(RT) = 1 / (1 + exp(-k (RT - t))).
#   observed replicated - 3-bin running mean of y above c + log2(1 + a) / 2
#   expected replicated - RT > t (p > 0.5), from K562 Repli-seq
# Ideogram: chr1-22 and X, two tracks per chromosome (left expected, right
# observed), only the expected-early (blue) and replicated (red) segments
# filled; late, unreplicated and masked bins (repeats, centromeres, low
# mappability) are white. Zoom: three 25 Mb windows on different chromosomes with the most
# early/late switches, showing y (grey points, red 3-bin mean) and RT (blue, right
# axis, scaled so the cell's threshold t meets the coverage midpoint on the
# dashed line), the expected-early regions shaded, and the window's place on
# the banded chromosome. Legends in a separate file.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggideogram)
})

rt_out <- file.path(project_root, "paper_figures/output/K562_replication_timing_coverage")
run_dir <- file.path(project_root, "results/K562_tree/aneufinder_fragment_gc")
fit_csv <- file.path(project_root, "results/K562_tree/aneufinder_fragment_gc_rtcorr/output/replication_fit.csv")
output_dir <- file.path(project_root, "paper_figures/output/K562_replication_ideogram/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
target <- as.numeric(Sys.getenv("K562_DNA_CONTENT", "1.25"))
label <- sprintf("dna%.2f", target)   # output file suffix
zoom_mb <- 25

col_state <- c("expected early" = "#2a6fb0", "expected late" = "#d9e3ef",
               "replicated" = "#c0392b", "not replicated" = "#f2d7d5", "no data" = "white")

# ---- cell and coverage ------------------------------------------------------------

fit <- read.csv(fit_csv)
dna <- read.csv(file.path(project_root, "paper_figures/output/K562_cell_cycle_profile/inferred_dna_content.csv"))
phase <- read.csv(file.path(rt_out, "per_cell_phase.csv"))
cell <- Sys.getenv("K562_CELL", "")
if (cell == "") {
  cand <- merge(dna, phase[, c("cell", "class")])
  cand <- cand[cand$class == "on-peak" & abs(cand$dna_content - target) <= 0.015, ]
  cell <- cand$cell[which.max(cand$r2)]
}
f <- fit[fit$cell == cell, ]
cat(sprintf("cell %s: DNA content %.3f, t = %.2f, k = %d, a = %.2f, r2 = %.3f\n", cell,
            dna$dna_content[dna$cell == cell], f$t, f$k, f$a, f$r2))

bt <- read.csv(file.path(rt_out, "bin_replication_timing.csv"))
kb <- read.csv(file.path(run_dir, "output/kept_bins.csv"))
fc <- readRDS(file.path(run_dir, "output/fragment_counts.rds"))
idx <- match(paste(bt$seqnames, bt$start), paste(kb$seqnames, kb$start))
d_all <- fc$counts_gc[idx, ] / kb$mappable_frac[idx] / bt$consensus_cn
y_all <- log2(sweep(d_all, 2, apply(d_all, 2, median), "/") + 1e-3)
ref <- apply(y_all[, phase$cell[phase$reference_cell]], 1, median)
bt$y <- y_all[, cell] - ref
bt$y <- bt$y - median(bt$y[bt$rt < f$t])   # unreplicated level at 0
run_mean <- function(v, ch) ave(v, ch, FUN = function(x) as.numeric(stats::filter(x, rep(1 / 3, 3), sides = 2)))
bt$y_smooth <- run_mean(bt$y, bt$seqnames)
bt$y_smooth[is.na(bt$y_smooth)] <- bt$y[is.na(bt$y_smooth)]
half <- log2(1 + f$a) / 2
bt$observed <- ifelse(bt$y_smooth > half, "replicated", "not replicated")
bt$expected <- ifelse(bt$rt > f$t, "expected early", "expected late")
agree <- mean((bt$observed == "replicated") == (bt$expected == "expected early"))
cat(sprintf("observed vs expected agreement: %.3f of %d bins\n", agree, nrow(bt)))

# ---- ideogram ----------------------------------------------------------------------

data(karyotype_hg38, package = "ggideogram")
data(chromosome.size.hg38, package = "ggideogram")
chroms <- c(paste0("chr", 1:22), "chrX")
sizes <- setNames(as.data.frame(chromosome.size.hg38)$End, as.character(as.data.frame(chromosome.size.hg38)$Chrom))
cen <- sapply(chroms, function(ch) {
  a <- karyotype_hg38[karyotype_hg38$Chrom == ch & karyotype_hg38$Stain == "acen", ]
  mean(c(min(a$Start), max(a$End)))
})

segments_for <- function(state_col) do.call(rbind, lapply(chroms, function(ch) {
  b <- bt[bt$seqnames == ch, ]
  r <- rle(b[[state_col]])
  e <- cumsum(r$lengths); s <- e - r$lengths + 1
  seg <- data.frame(Chrom = ch, Start = b$start[s] - 1, End = b$end[e], state = r$values)
  # gaps between kept bins and at the chromosome ends: no data
  edges <- sort(unique(c(0, seg$Start, seg$End, sizes[[ch]])))
  covered <- IRanges::IRanges(seg$Start + 1, seg$End)
  gaps <- IRanges::gaps(covered, start = 1, end = sizes[[ch]])
  rbind(seg, if (length(gaps)) data.frame(Chrom = ch, Start = IRanges::start(gaps) - 1,
                                          End = IRanges::end(gaps), state = "no data"))
}))
ideo <- rbind(cbind(segments_for("expected"), track = "expected"),
              cbind(segments_for("observed"), track = "observed"))
ideo <- ideo[ideo$state %in% c("expected early", "replicated"), ]   # only the replicated/early segments are filled

# x: discrete levels per chromosome (expected, observed, spacer), which keeps
# ggideogram's outlines aligned; y negative so chromosomes hang from the top
# (a reversed y scale distorts ggideogram's rounded ends)
track_w <- 0.75
lv <- as.vector(t(outer(chroms, c("expected", "observed", "spacer"), paste)))
lv <- head(lv, -1)
ideo$xi <- match(paste(ideo$Chrom, ideo$track), lv)
outline <- data.frame(Chrom = rep(chroms, each = 2), track = rep(c("expected", "observed"), length(chroms)))
outline$key <- factor(paste(outline$Chrom, outline$track), levels = lv)
outline$Start <- 0; outline$End <- sizes[outline$Chrom]

p_ideo <- ggplot() +
  geom_ideogram(data = outline, aes(x = key, ymin = -End, ymax = -Start, chrom = key),
                fill = "white", colour = NA, radius = unit(3, "pt"), width = track_w, show.legend = FALSE) +
  geom_rect(data = ideo, aes(xmin = xi - track_w / 2 * 0.9, xmax = xi + track_w / 2 * 0.9,
                             ymin = -End, ymax = -Start, fill = state)) +
  geom_ideogram(data = outline, aes(x = key, ymin = -End, ymax = -Start, chrom = key),
                fill = NA, colour = "black", radius = unit(3, "pt"), width = track_w,
                linewidth = 0.25, show.legend = FALSE) +
  annotate("text", x = match(paste(chroms, "expected"), lv) + 0.5, y = -max(sizes[chroms]) * 1.03,
           label = sub("chr", "", chroms), size = 2.2) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = col_state, guide = "none") +
  theme_void()
ggsave(file.path(output_dir, paste0("replication_ideogram_", label, ".pdf")), p_ideo, width = 9, height = 4)
ggsave(file.path(output_dir, paste0("replication_ideogram_", label, ".svg")), p_ideo, width = 9, height = 4)
ggsave(file.path(output_dir, paste0("replication_ideogram_", label, ".png")), p_ideo, width = 9, height = 4, dpi = 200, bg = "white")

# ---- zoomed tracks ----------------------------------------------------------------

windows <- do.call(rbind, lapply(chroms, function(ch) {
  b <- bt[bt$seqnames == ch, ]
  starts <- seq(0, max(b$end) - zoom_mb * 1e6, by = 5e6)
  do.call(rbind, lapply(starts, function(s) {
    w <- b[b$start >= s & b$end <= s + zoom_mb * 1e6, ]
    if (nrow(w) < 0.8 * zoom_mb * 4) return(NULL)
    data.frame(chrom = ch, start = s, switches = sum(diff(w$rt > f$t) != 0),
               agree = mean((w$observed == "replicated") == (w$expected == "expected early")))
  }))
}))
windows <- windows[order(-windows$switches, -windows$agree), ]
zoom <- windows[!duplicated(windows$chrom), ][1:3, ]
write.csv(zoom, file.path(output_dir, paste0("zoom_windows_", label, ".csv")), row.names = FALSE)

zoom_panel <- function(ch, s) {
  e <- s + zoom_mb * 1e6
  w <- bt[bt$seqnames == ch & bt$start >= s & bt$end <= e, ]
  mb <- function(x) x / 1e6
  runs <- rle(w$rt > f$t); en <- cumsum(runs$lengths); st <- en - runs$lengths + 1
  shade <- data.frame(xmin = mb(w$start[st]), xmax = mb(w$end[en]))[runs$values, ]
  # replication timing on a second axis, scaled so the cell's threshold t sits
  # on the coverage midpoint (half-way to doubled) and RT 0-1 fits the panel
  ylim <- c(-0.8, 1.6)
  sc <- min((ylim[2] - half) / (1 - f$t), (half - ylim[1]) / f$t)
  to_cov <- function(r) half + sc * (r - f$t)
  to_rt <- function(y) (y - half) / sc + f$t
  p_cov <- ggplot(w) +
    geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = col_state[["expected early"]], alpha = 0.18, inherit.aes = FALSE) +
    geom_hline(yintercept = c(0, log2(1 + f$a)), linetype = "dotted", linewidth = 0.25, colour = "grey50") +
    geom_hline(yintercept = half, linetype = "dashed", linewidth = 0.25, colour = "grey40") +
    geom_line(aes(mb((start + end) / 2), to_cov(rt)), linewidth = 0.45, colour = col_state[["expected early"]]) +
    geom_point(aes(mb((start + end) / 2), y), size = 0.4, colour = "grey30") +
    geom_line(aes(mb((start + end) / 2), y_smooth), linewidth = 0.45, colour = col_state[["replicated"]]) +
    scale_y_continuous(name = "log2 coverage", sec.axis = sec_axis(to_rt, name = "Replication timing",
                                                                   breaks = c(0, 0.25, 0.5, 0.75, 1))) +
    coord_cartesian(xlim = c(mb(s), mb(e)), ylim = ylim) +
    labs(x = sprintf("%s (Mb)", ch)) + theme_classic(base_size = 7)
  kt <- karyotype_hg38[karyotype_hg38$Chrom == ch, ]
  p_chr <- ggplot(kt) +
    geom_ideogram(aes(x = Chrom, ymin = Start, ymax = End, chrom = Chrom, fill = Stain, arm = Arm),
                  radius = unit(3, "pt"), width = 0.6, linewidth = 0.15, colour = "black", show.legend = FALSE) +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = s, ymax = e, fill = NA, colour = "#c0392b", linewidth = 0.5) +
    scale_fill_manual(values = cytoband_colors) + coord_flip() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.background = element_blank())
  p_chr / p_cov + plot_layout(heights = c(0.5, 3.6))
}
p_zoom <- wrap_plots(lapply(seq_len(nrow(zoom)), function(i) zoom_panel(zoom$chrom[i], zoom$start[i])), nrow = 1)
ggsave(file.path(output_dir, paste0("replication_zoom_", label, ".pdf")), p_zoom, width = 9, height = 3.5)
ggsave(file.path(output_dir, paste0("replication_zoom_", label, ".svg")), p_zoom, width = 9, height = 3.5)
ggsave(file.path(output_dir, paste0("replication_zoom_", label, ".png")), p_zoom, width = 9, height = 3.5, dpi = 200, bg = "white")

# ---- arm level: called above baseline vs expected early ---------------------------
# per arm (>= 20 kept bins): fraction of bins called replicated (coverage above
# baseline) against fraction expected early (RT > t)

bt$arm <- paste0(bt$seqnames, ifelse((bt$start + bt$end) / 2 < cen[bt$seqnames], "p", "q"))
arms <- do.call(rbind, lapply(split(bt, bt$arm), function(a) {
  if (nrow(a) < 20) return(NULL)
  data.frame(arm = a$arm[1], chrom = a$seqnames[1], n_bins = nrow(a),
             frac_expected_early = mean(a$expected == "expected early"),
             frac_observed_replicated = mean(a$observed == "replicated"),
             bin_agreement = mean((a$expected == "expected early") == (a$observed == "replicated")),
             mean_rt = mean(a$rt), mean_coverage = mean(a$y))
}))
arms$chrom <- factor(arms$chrom, levels = chroms)
write.csv(arms, file.path(output_dir, paste0("arms_", cell, ".csv")), row.names = FALSE)
r_arm <- cor(arms$frac_expected_early, arms$frac_observed_replicated)
cat(sprintf("%d arms: Pearson r(fraction expected early, fraction called above baseline) = %.3f; per-arm bin agreement median %.3f (range %.2f-%.2f)\n",
            nrow(arms), r_arm, median(arms$bin_agreement), min(arms$bin_agreement), max(arms$bin_agreement)))
chrom_cols <- setNames(viridisLite::viridis(length(chroms)), chroms)
lim <- c(0, max(c(arms$frac_expected_early, arms$frac_observed_replicated)) * 1.05)
p_arm <- ggplot(arms, aes(frac_expected_early, frac_observed_replicated, colour = chrom)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.8) +
  scale_colour_manual(values = chrom_cols, guide = "none") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Fraction of arm expected early (Repli-seq)", y = "Fraction of arm above baseline") +
  theme_classic(base_size = 8)
ggsave(file.path(output_dir, paste0("arm_coverage_vs_rt_", label, ".pdf")), p_arm, width = 3.2, height = 3)
ggsave(file.path(output_dir, paste0("arm_coverage_vs_rt_", label, ".svg")), p_arm, width = 3.2, height = 3)
ggsave(file.path(output_dir, paste0("arm_coverage_vs_rt_", label, ".png")), p_arm, width = 3.2, height = 3, dpi = 200, bg = "white")
leg_chr <- cowplot::get_plot_component(
  ggplot(data.frame(c = factor(chroms, chroms))) + geom_point(aes(1, 1, colour = c), size = 2) +
    scale_colour_manual(values = chrom_cols, labels = sub("chr", "", chroms), name = "Chromosome") +
    guides(colour = guide_legend(ncol = 2)) + theme(legend.position = "right"),
  "guide-box-right")
ggsave(file.path(output_dir, paste0("arm_coverage_vs_rt_legend_", label, ".pdf")), leg_chr, width = 1.4, height = 3)

# ---- legends (separate file) -------------------------------------------------------

leg <- cowplot::get_plot_component(
  ggplot(data.frame(s = factor(names(col_state), names(col_state)))) + geom_tile(aes(1, 1, fill = s)) +
    scale_fill_manual(values = col_state, name = NULL) + theme(legend.position = "right"),
  "guide-box-right")
ggsave(file.path(output_dir, paste0("replication_ideogram_legend_", label, ".pdf")), leg, width = 1.6, height = 1.4)
write.csv(bt[, c("seqnames", "start", "end", "rt", "y", "y_smooth", "expected", "observed")],
          file.path(output_dir, paste0("bins_", cell, ".csv")), row.names = FALSE)
cat("Outputs written to", output_dir, "\n")
