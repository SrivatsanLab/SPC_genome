#!/usr/bin/env Rscript

# Re-segment the K562 sc_PolE cells on repeat-masked, smoothed 250 kb bins.
#
# Starts from the GC-corrected 250 kb binned counts of an existing AneuFinder
# run (ANEUFINDER_SOURCE, default results/K562_tree/aneufinder_gc_2.5e5/output;
# MAPQ >= 10, duplicates removed, ENCODE blacklist, GC-corrected). Then, per cell:
#   1. drop bins >= MASK_MIN covered by the hg38 repeat mask (centromeres, gaps,
#      Satellite repeats, segmental duplications; make_hg38_repeat_mask.R) or
#      by N gaps + ENCODE blacklist; masking whole bins before segmentation
#      avoids the false losses that read-level blacklisting gives partly
#      masked bins
#   1b. (GC_RECORRECT=1) second, per-cell GC correction: each cell's relative
#      depth divided by the consensus (median of on-peak cells) is fitted
#      against bin GC with a robust loess (span 0.3), and the cell's counts are
#      divided by the fitted trend; dividing by the consensus first keeps the
#      shared K562 karyotype out of the fit, and the robust fit keeps small
#      real CNVs from moving it. Removes GC dependence left after AneuFinder's
#      correction so it cannot drive the clustering
#   2. smooth: each kept bin takes the mean of the two overlapping 500 kb
#      windows (250 kb step) that contain it, i.e. 0.5 * left + self +
#      0.5 * right, using only kept, adjacent neighbours (a missing neighbour
#      is replaced by the bin itself)
#   3. AneuFinder edivisive segmentation (findCNVs, R = 10, sig.lvl = 0.1, as in
#      the original runs) on the kept bins, chromosome by chromosome
# Writes MODELS/method-edivisive/*.RData (aneuHMM objects), result.csv (kept
# bins x cells, same layout as the other runs) and kept_bins.csv.

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
  library(parallel)
})

repo <- "/home/sanjay/SPC_genome"
source_dir <- Sys.getenv("ANEUFINDER_SOURCE", file.path(repo, "results/K562_tree/aneufinder_gc_2.5e5/output"))
out_dir <- Sys.getenv("ANEUFINDER_OUTPUT", file.path(repo, "results/K562_tree/aneufinder_gc_smooth500k/output"))
mask_min <- as.numeric(Sys.getenv("MASK_MIN", "0.5"))
gc_recorrect <- Sys.getenv("GC_RECORRECT", "0") == "1"
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
models_dir <- file.path(out_dir, "MODELS", "method-edivisive")
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

load_gr <- function(f) {
  e <- new.env(); load(f, envir = e); g <- get(ls(e)[1], envir = e)
  if (!is(g, "GRanges")) g <- g[[1]]
  g
}
files <- list.files(file.path(source_dir, "binned-GC"), "\\.RData$", full.names = TRUE)
n_test <- as.integer(Sys.getenv("N_TEST", "0"))   # > 0: only the first n cells, for testing
if (n_test > 0) files <- head(files, n_test)
bins <- load_gr(files[1])
bins <- bins[seqnames(bins) %in% c(paste0("chr", 1:22), "chrX", "chrY")]

# ---- mask --------------------------------------------------------------------

covered_frac <- function(gr, bed) {
  ov <- findOverlaps(gr, bed)
  bp <- tapply(width(pintersect(gr[queryHits(ov)], bed[subjectHits(ov)])), queryHits(ov), sum)
  f <- numeric(length(gr)); f[as.integer(names(bp))] <- bp / width(gr)[as.integer(names(bp))]
  f
}
rep_bed <- read.table(file.path(repo, "results/K562_tree/masks/hg38_repeat_mask.bed"))
rep_gr <- reduce(GRanges(rep_bed$V1, IRanges(rep_bed$V2 + 1, rep_bed$V3)))
unmap <- read.csv(file.path(repo, "results/K562_tree/masks/unmappable_bins_2.5e+05.csv"))
key <- paste(seqnames(bins), start(bins))
unmap_frac <- (unmap$n_frac + unmap$blacklist_frac)[match(key, paste(unmap$seqnames, unmap$start))]
repeat_frac <- covered_frac(bins, rep_gr)
kept <- repeat_frac < mask_min & unmap_frac < mask_min
write.csv(data.frame(seqnames = as.character(seqnames(bins)), start = start(bins), end = end(bins),
                     repeat_frac = repeat_frac, unmappable_frac = unmap_frac, kept = kept),
          file.path(out_dir, "kept_bins.csv"), row.names = FALSE)
cat(sprintf("%d of %d bins kept (repeat or unmappable fraction < %.2f)\n", sum(kept), length(kept), mask_min))
kept_bins <- bins[kept]
chr <- as.character(seqnames(kept_bins))
# neighbours count only if kept, on the same chromosome and directly adjacent
left_ok <- c(FALSE, chr[-1] == chr[-length(chr)] & start(kept_bins)[-1] == end(kept_bins)[-length(chr)] + 1)
right_ok <- c(left_ok[-1], FALSE)

# ---- optional second GC correction --------------------------------------------

cell_ids <- sub("_binsize.*", "", basename(files))
gc_factor <- NULL
if (gc_recorrect) {
  counts_all <- do.call(cbind, mclapply(files, function(f) {
    g <- load_gr(f); mcols(g)$counts[match(key[kept], paste(seqnames(g), start(g)))]
  }, mc.cores = ncpu))
  colnames(counts_all) <- cell_ids
  rel <- sweep(counts_all, 2, colMeans(counts_all), "/")
  meta <- read.csv(file.path(repo, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
  on_peak <- sub("\\.bam$", "", cell_ids) %in% meta$cell[meta$kept]
  consensus <- apply(rel[, on_peak, drop = FALSE], 1, median)
  gc <- mcols(kept_bins)$GC
  fit_ok <- consensus > 0 & is.finite(gc)
  gc_factor <- do.call(cbind, mclapply(seq_len(ncol(rel)), function(j) {
    y <- log2((rel[, j] + 1e-3) / (consensus + 1e-3))
    fit <- loess(y ~ gc, subset = fit_ok, span = 0.3, family = "symmetric", degree = 1)
    2^predict(fit, newdata = data.frame(gc = pmin(pmax(gc, min(gc[fit_ok])), max(gc[fit_ok]))))
  }, mc.cores = ncpu))
  colnames(gc_factor) <- cell_ids
  cat(sprintf("second GC correction: median per-cell range of fitted GC factor %.3f\n",
              median(apply(gc_factor, 2, function(x) diff(range(x))))))
}

smooth <- function(v) {
  l <- c(v[1], v[-length(v)]); r <- c(v[-1], v[length(v)])
  l[!left_ok] <- v[!left_ok]; r[!right_ok] <- v[!right_ok]
  (0.5 * l + v + 0.5 * r) / 2
}

# ---- per cell ----------------------------------------------------------------

run_cell <- function(f) {
  id <- sub("_binsize.*", "", basename(f))
  g <- load_gr(f)
  g <- g[match(key[kept], paste(seqnames(g), start(g)))]
  b <- kept_bins
  mcols(b) <- NULL
  fac <- if (is.null(gc_factor)) 1 else gc_factor[, id]
  for (col in c("counts", "mcounts", "pcounts")) mcols(b)[[col]] <- smooth(mcols(g)[[col]] / fac)
  attr(b, "ID") <- id
  model <- suppressMessages(suppressWarnings(
    findCNVs(b, ID = id, method = "edivisive", R = 10, sig.lvl = 0.1, verbosity = 0)))
  save(model, file = file.path(models_dir, paste0(id, "_binsize_2.5e+05_smooth500k_masked.RData")))
  data.frame(id = id, cn = I(list(model$bins$copy.number)))
}
t0 <- Sys.time()
res <- mclapply(files, function(f) tryCatch(run_cell(f), error = function(e) {
  message(basename(f), ": ", conditionMessage(e)); NULL }), mc.cores = ncpu)
res <- do.call(rbind, res[!sapply(res, is.null)])
cat(sprintf("segmented %d cells in %.1f min\n", nrow(res), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cn <- do.call(cbind, res$cn)
colnames(cn) <- sub("\\.bam$", "", res$id)
result <- cbind(as.data.frame(kept_bins)[, c("seqnames", "start", "end")], cn)
write.csv(result, file.path(out_dir, "result.csv"))
cat(sprintf("Wrote %s (%d bins x %d cells); copy number > 10 in %.4f of cell-bins\n",
            file.path(out_dir, "result.csv"), nrow(cn), ncol(cn), mean(cn > 10)))
