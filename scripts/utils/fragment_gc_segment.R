#!/usr/bin/env Rscript

# Segment the K562 sc_PolE cells on fragment-level GC-corrected depth
# (scripts/utils/fragment_gc_counts.R, per-cell RDS in
# results/K562_tree/fragment_gc/cells/).
#
#   1. bins: the 250 kb grid, dropping bins >= 50% repeat-masked (centromeres,
#      gaps, Satellite repeats, segmental duplications) or >= 50% N gap +
#      ENCODE blacklist, or < MIN_MAPPABLE (default 0.9) uniquely mappable
#      (Umap k100); kept bins' depth is divided by their mappable fraction
#   2. smoothing: each kept bin takes the mean of the two overlapping 500 kb
#      windows (250 kb step) containing it, 0.5 * left + self + 0.5 * right,
#      over kept, adjacent neighbours (a missing neighbour is replaced by the
#      bin itself)
#   3. AneuFinder edivisive segmentation (findCNVs, R = 10, sig.lvl = 0.1)
# Writes to results/K562_tree/aneufinder_fragment_gc/output/: counts matrices
# (raw and GC-weighted, all bins), kept_bins.csv, MODELS/method-edivisive/*.RData
# and result.csv (kept bins x cells, same layout as the other runs).

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
  library(parallel)
})

repo <- "/home/sanjay/SPC_genome"
gc_dir <- file.path(repo, "results/K562_tree/fragment_gc")
out_dir <- Sys.getenv("ANEUFINDER_OUTPUT", file.path(repo, "results/K562_tree/aneufinder_fragment_gc/output"))
min_mappable <- as.numeric(Sys.getenv("MIN_MAPPABLE", "0.9"))
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
models_dir <- file.path(out_dir, "MODELS", "method-edivisive")
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

bins <- readRDS(file.path(gc_dir, "bins.rds"))
cell_files <- list.files(file.path(gc_dir, "cells"), "\\.rds$", full.names = TRUE)
cells <- lapply(cell_files, readRDS)
ids <- sapply(cells, `[[`, "cell")
counts_gc <- sapply(cells, `[[`, "counts_gc"); colnames(counts_gc) <- ids
counts_raw <- sapply(cells, `[[`, "counts_raw"); colnames(counts_raw) <- ids
bin_df <- data.frame(seqnames = as.character(seqnames(bins)), start = start(bins), end = end(bins))
saveRDS(list(bins = bin_df, counts_gc = counts_gc, counts_raw = counts_raw,
             n_fragments = sapply(cells, `[[`, "n_fragments"),
             gc_tables = lapply(cells, `[[`, "gc_table")),
        file.path(out_dir, "fragment_counts.rds"))
cat(sprintf("%d cells; median %.0f unique fragments per cell\n", length(ids),
            median(sapply(cells, `[[`, "n_fragments"))))

# ---- mask and mappability ---------------------------------------------------

kept <- bins$repeat_frac < 0.5 & bins$unmappable_frac < 0.5 & bins$mappable_frac >= min_mappable
write.csv(cbind(bin_df, as.data.frame(mcols(bins)), kept = kept), file.path(out_dir, "kept_bins.csv"),
          row.names = FALSE)
cat(sprintf("%d of %d bins kept (repeat < 0.5, unmappable < 0.5, mappable >= %.2f)\n",
            sum(kept), length(kept), min_mappable))
kept_bins <- bins[kept]
mcols(kept_bins) <- NULL
depth <- counts_gc[kept, , drop = FALSE] / bins$mappable_frac[kept]

chr <- as.character(seqnames(kept_bins))
left_ok <- c(FALSE, chr[-1] == chr[-length(chr)] & start(kept_bins)[-1] == end(kept_bins)[-length(chr)] + 1)
right_ok <- c(left_ok[-1], FALSE)
smooth <- function(v) {
  l <- c(v[1], v[-length(v)]); r <- c(v[-1], v[length(v)])
  l[!left_ok] <- v[!left_ok]; r[!right_ok] <- v[!right_ok]
  (0.5 * l + v + 0.5 * r) / 2
}

# ---- segment ----------------------------------------------------------------

run_cell <- function(j) {
  id <- ids[j]
  b <- kept_bins
  v <- smooth(depth[, j])
  b$counts <- v; b$mcounts <- v / 2; b$pcounts <- v / 2
  attr(b, "ID") <- id
  model <- suppressMessages(suppressWarnings(
    findCNVs(b, ID = id, method = "edivisive", R = 10, sig.lvl = 0.1, verbosity = 0)))
  save(model, file = file.path(models_dir, paste0(id, "_binsize_2.5e+05_fragmentGC_smooth500k_masked.RData")))
  model$bins$copy.number
}
t0 <- Sys.time()
cn <- mclapply(seq_along(ids), function(j) tryCatch(run_cell(j), error = function(e) {
  message(ids[j], ": ", conditionMessage(e)); NULL }), mc.cores = ncpu)
ok <- !sapply(cn, is.null)
cn <- do.call(cbind, cn[ok]); colnames(cn) <- ids[ok]
cat(sprintf("segmented %d cells in %.1f min\n", ncol(cn), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
result <- cbind(as.data.frame(kept_bins)[, c("seqnames", "start", "end")], cn)
write.csv(result, file.path(out_dir, "result.csv"))
cat(sprintf("Wrote %s (%d bins x %d cells); copy number > 10 in %.4f of cell-bins\n",
            file.path(out_dir, "result.csv"), nrow(cn), ncol(cn), mean(cn > 10)))
