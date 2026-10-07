# Per-cell coverage noise within copy-number segments, to flag poorly
# amplified K562 sc_PolE cells.
#
# For each cell, unique fragments per kept 250 kb bin (fragment_gc_counts.R;
# divided by the bin's mappable fraction) are compared with the mean of the
# cell's own copy-number segment (replication-corrected calls,
# scripts/utils/replication_regress_segment.R), so copy-number structure does
# not count as noise. Reported:
#   cv_within        - sqrt(mean((x - mu)^2 / mu) / mean(x)): coefficient of
#                      variation within segments
#   cv_poisson       - 1 / sqrt(mean count): what counting noise alone gives
#   dispersion       - within-segment variance / Poisson variance
#   autocor          - correlation of neighbouring within-segment residuals
#   dropout          - fraction of bins below 25% of their segment mean
# noisy = cv_within > K562_NOISE_MAX (default 0.3). Writes
# paper_figures/output/K562_cell_noise/per_cell_noise.csv and a histogram.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
})

noise_max <- as.numeric(Sys.getenv("K562_NOISE_MAX", "0.3"))
run_dir <- file.path(project_root, "results/K562_tree/aneufinder_fragment_gc_rtcorr/output")
output_dir <- file.path(project_root, "paper_figures/output/K562_cell_noise/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fc <- readRDS(file.path(project_root, "results/K562_tree/aneufinder_fragment_gc/output/fragment_counts.rds"))
kb <- read.csv(file.path(project_root, "results/K562_tree/aneufinder_fragment_gc/output/kept_bins.csv"))
key_all <- paste(kb$seqnames, kb$start)
seg_files <- list.files(file.path(run_dir, "MODELS/method-edivisive"), full.names = TRUE)
names(seg_files) <- sub("_binsize.*", "", basename(seg_files))
cells <- intersect(names(seg_files), colnames(fc$counts_raw))

noise <- do.call(rbind, lapply(cells, function(cl) {
  e <- new.env(); load(seg_files[[cl]], envir = e); b <- e$model$bins
  i <- match(paste(seqnames(b), start(b)), key_all)
  x <- fc$counts_raw[i, cl] / kb$mappable_frac[i]
  chr <- as.character(seqnames(b))
  seg <- cumsum(c(TRUE, b$copy.number[-1] != b$copy.number[-length(b)] | chr[-1] != chr[-length(b)]))
  mu <- ave(x, seg, FUN = mean); ok <- mu > 0
  disp <- sum((x - mu)[ok]^2 / mu[ok]) / sum(ok)
  r <- (x - mu) / sqrt(pmax(mu, 1)); same <- seg[-1] == seg[-length(seg)]
  data.frame(cell = cl, mean_count = mean(x), cv_poisson = 1 / sqrt(mean(x)),
             cv_within = sqrt(disp / mean(x)), dispersion = disp,
             autocor = cor(r[-1][same], r[-length(r)][same]),
             dropout = mean(x[ok] < 0.25 * mu[ok]), segments = max(seg))
}))
noise$noisy <- noise$cv_within > noise_max

metrics <- read.csv(file.path(project_root, "paper_figures/output/K562_offpeak_ploidy/per_cell_metrics.csv"))
ph <- read.csv(file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/per_cell_phase.csv"))
meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
noise$class <- metrics$class[match(noise$cell, metrics$cell)]
noise$phase <- ph$phase[match(noise$cell, ph$cell)]
noise$pop <- meta$pop[match(noise$cell, meta$cell)]
write.csv(noise, file.path(output_dir, "per_cell_noise.csv"), row.names = FALSE)

cat(sprintf("%d cells; within-segment CV quantiles 5/25/50/75/95/99%%: %s\n", nrow(noise),
            paste(round(quantile(noise$cv_within, c(.05, .25, .5, .75, .95, .99)), 3), collapse = " ")))
cat(sprintf("noisy (CV > %.2f): %d cells\n", noise_max, sum(noise$noisy)))
print(table(class = noise$class, noisy = noise$noisy))
print(table(phase = noise$phase, noisy = noise$noisy))
print(table(pop = noise$pop, noisy = noise$noisy))

p <- ggplot(noise, aes(cv_within)) +
  geom_histogram(binwidth = 0.01, fill = "grey60", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = noise_max, linetype = "dashed", colour = "#c0392b") +
  labs(x = "Within-segment coverage CV", y = "Cells") + theme_classic(base_size = 8)
ggsave(file.path(output_dir, "cv_within_histogram.pdf"), p, width = 3.2, height = 2.4)
ggsave(file.path(output_dir, "cv_within_histogram.png"), p, width = 3.2, height = 2.4, dpi = 200, bg = "white")
