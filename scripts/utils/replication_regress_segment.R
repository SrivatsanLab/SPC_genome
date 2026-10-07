#!/usr/bin/env Rscript

# Remove the cell-cycle (replication) component from each K562 sc_PolE cell's
# coverage, then re-call copy number.
#
# Inputs: fragment-level GC-weighted depth on masked 250 kb bins
# (fragment_gc_counts.R / fragment_gc_segment.R), K562 replication timing per
# bin and the non-replicating reference cells
# (paper_figures/scripts/K562_replication_timing_coverage.R).
#
# Per cell, on bins with replication timing:
#   y = log2(depth / mappable / consensus CN / cell median) - reference profile
#   where the reference profile is the median y of the non-replicating
#   reference cells (shared regional amplification bias).
#   Model: y = c + log2(1 + a * p(RT)),  p(RT) = 1 / (1 + exp(-k (RT - t)))
#   p is the share of a region already replicated at the cell's point in S;
#   a (0-1) how strongly the doubling shows. Fitted by least squares on the
#   medians of y in 50 RT classes (equal bin counts), so copy-number changes
#   confined to a minority of each class do not move it; grid over t (RT
#   0.05-0.95), k (5, 10, 20, 40) and a (0-1, step 0.02), c the weighted mean
#   residual. Replicated fraction = mean p over bins.
# Correction applied to the cell's own depth (its copy-number changes kept):
#   depth_corrected = depth / mappable / 2^reference / (1 + a * p(RT))
# Then the usual smoothing (500 kb windows, 250 kb step) and AneuFinder
# edivisive segmentation (findCNVs, R = 10, sig.lvl = 0.1).
# Writes results/K562_tree/aneufinder_fragment_gc_rtcorr/output/:
# replication_fit.csv, MODELS/method-edivisive/*.RData, result.csv,
# corrected_depth.rds, kept_bins.csv.

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
  library(parallel)
})

repo <- "/home/sanjay/SPC_genome"
src <- file.path(repo, "results/K562_tree/aneufinder_fragment_gc/output")
rt_out <- file.path(repo, "paper_figures/output/K562_replication_timing_coverage")
out_dir <- Sys.getenv("ANEUFINDER_OUTPUT", file.path(repo, "results/K562_tree/aneufinder_fragment_gc_rtcorr/output"))
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
models_dir <- file.path(out_dir, "MODELS", "method-edivisive")
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

# ---- inputs -----------------------------------------------------------------------

kb <- read.csv(file.path(src, "kept_bins.csv"))
fc <- readRDS(file.path(src, "fragment_counts.rds"))
bt <- read.csv(file.path(rt_out, "bin_replication_timing.csv"))   # kept bins with RT and consensus CN
phase <- read.csv(file.path(rt_out, "per_cell_phase.csv"))
cells <- colnames(fc$counts_gc)

idx <- match(paste(bt$seqnames, bt$start), paste(kb$seqnames, kb$start))
depth <- fc$counts_gc[idx, ] / kb$mappable_frac[idx]
rt <- bt$rt
y <- log2(sweep(depth / bt$consensus_cn, 2, apply(depth / bt$consensus_cn, 2, median), "/") + 1e-3)
ref_cells <- phase$cell[phase$reference_cell]
ref_profile <- apply(y[, ref_cells], 1, median)
y <- y - ref_profile
cat(sprintf("%d bins with replication timing; %d reference cells\n", length(rt), length(ref_cells)))

# ---- per-cell replication fit -------------------------------------------------------

rt_class <- cut(rank(rt, ties.method = "first"), 50, labels = FALSE)
rt_mid <- as.numeric(tapply(rt, rt_class, median))
n_class <- tabulate(rt_class, 50)
grid <- expand.grid(t = seq(0.05, 0.95, 0.01), k = c(5, 10, 20, 40))
P <- sapply(seq_len(nrow(grid)), function(g) 1 / (1 + exp(-grid$k[g] * (rt_mid - grid$t[g]))))
a_grid <- seq(0, 1, 0.02)

fit_cell <- function(j) {
  m <- as.numeric(tapply(y[, j], rt_class, median))
  best <- list(rss = Inf)
  for (g in seq_len(ncol(P))) {
    F <- log2(1 + outer(P[, g], a_grid))              # classes x a values
    c0 <- colSums(n_class * (m - F)) / sum(n_class)
    rss <- colSums(n_class * sweep(m - F, 2, c0)^2)
    i <- which.min(rss)
    if (rss[i] < best$rss) best <- list(rss = rss[i], t = grid$t[g], k = grid$k[g], a = a_grid[i], c = c0[i])
  }
  null_rss <- sum(n_class * (m - weighted.mean(m, n_class))^2)
  p_bins <- 1 / (1 + exp(-best$k * (rt - best$t)))
  data.frame(cell = cells[j], t = best$t, k = best$k, a = best$a,
             replicated_frac = if (best$a > 0) mean(p_bins) else NA,
             r2 = 1 - best$rss / null_rss)
}
fits <- do.call(rbind, mclapply(seq_along(cells), fit_cell, mc.cores = ncpu))
fits$phase <- phase$phase[match(fits$cell, phase$cell)]
write.csv(fits, file.path(out_dir, "replication_fit.csv"), row.names = FALSE)
cat("fitted amplitude a by phase call (median, IQR):\n")
print(aggregate(a ~ phase, fits, function(x) round(c(median = median(x), q25 = quantile(x, .25), q75 = quantile(x, .75)), 3)))

# ---- corrected depth --------------------------------------------------------------

rep_factor <- sapply(seq_along(cells), function(j) {
  f <- fits[j, ]
  1 + f$a / (1 + exp(-f$k * (rt - f$t)))
})
corrected <- depth / 2^ref_profile / rep_factor
colnames(corrected) <- cells
bins_out <- bt[, c("seqnames", "start", "end")]
saveRDS(list(bins = bins_out, corrected = corrected, ref_profile = ref_profile, rep_factor = rep_factor),
        file.path(out_dir, "corrected_depth.rds"))
write.csv(bt, file.path(out_dir, "kept_bins.csv"), row.names = FALSE)

# ---- smoothing and segmentation ---------------------------------------------------

gr <- GRanges(bins_out$seqnames, IRanges(bins_out$start, bins_out$end))
seqlevels(gr) <- unique(bins_out$seqnames)
# chromosome lengths, which AneuFinder's heatmaps need to lay out the genome
seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]
chr <- bins_out$seqnames
left_ok <- c(FALSE, chr[-1] == chr[-length(chr)] & bins_out$start[-1] == bins_out$end[-length(chr)] + 1)
right_ok <- c(left_ok[-1], FALSE)
smooth <- function(v) {
  l <- c(v[1], v[-length(v)]); r <- c(v[-1], v[length(v)])
  l[!left_ok] <- v[!left_ok]; r[!right_ok] <- v[!right_ok]
  (0.5 * l + v + 0.5 * r) / 2
}
# segmentation needs count-like magnitudes: rescale each cell to its raw mean depth
scale_to <- colMeans(depth)
run_cell <- function(j) {
  id <- cells[j]
  v <- smooth(corrected[, j]); v <- v / mean(v) * scale_to[j]
  b <- gr; b$counts <- v; b$mcounts <- v / 2; b$pcounts <- v / 2
  attr(b, "ID") <- id
  model <- suppressMessages(suppressWarnings(
    findCNVs(b, ID = id, method = "edivisive", R = 10, sig.lvl = 0.1, verbosity = 0)))
  save(model, file = file.path(models_dir, paste0(id, "_binsize_2.5e+05_fragmentGC_rtcorr.RData")))
  model$bins$copy.number
}
t0 <- Sys.time()
cn <- mclapply(seq_along(cells), function(j) tryCatch(run_cell(j), error = function(e) {
  message(cells[j], ": ", conditionMessage(e)); NULL }), mc.cores = ncpu)
ok <- !sapply(cn, is.null)
cn <- do.call(cbind, cn[ok]); colnames(cn) <- cells[ok]
cat(sprintf("segmented %d cells in %.1f min\n", ncol(cn), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(cbind(bins_out, cn), file.path(out_dir, "result.csv"))
cat(sprintf("Wrote %s (%d bins x %d cells); copy number > 10 in %.4f of cell-bins\n",
            file.path(out_dir, "result.csv"), nrow(cn), ncol(cn), mean(cn > 10)))
