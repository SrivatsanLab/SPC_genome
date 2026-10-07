#!/usr/bin/env Rscript

# Copy number of the K562 sc_PolE cells on 1 Mb bins from the repeat- and
# mappability-masked, fragment-GC-weighted, replication-corrected 250 kb depth
# (scripts/utils/replication_regress_segment.R, corrected_depth.rds).
#
# 1 Mb bins are the AneuFinder grid (1-1,000,000, ...), each the union of four
# 250 kb bins. A 1 Mb bin's depth is the mean corrected depth of its kept
# 250 kb sub-bins (masked sub-bins - >= 50% repeat or unmappable, < 90%
# uniquely mappable, no replication timing - are left out rather than counted
# as zero); bins with fewer than 2 kept sub-bins are dropped. No smoothing.
# Segmentation as in the original 1 Mb run: AneuFinder edivisive (findCNVs,
# R = 10, sig.lvl = 0.1), with each cell's depth scaled to its raw fragment
# depth per 1 Mb so the counts are on a count-like scale.
# Writes results/K562_tree/aneufinder_fragment_gc_rtcorr_1mb/output/:
# MODELS/method-edivisive/*.RData, result.csv, bins_1mb.csv.

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
  library(parallel)
})

repo <- "/home/sanjay/SPC_genome"
src <- file.path(repo, "results/K562_tree/aneufinder_fragment_gc_rtcorr/output")
out_dir <- Sys.getenv("ANEUFINDER_OUTPUT", file.path(repo, "results/K562_tree/aneufinder_fragment_gc_rtcorr_1mb/output"))
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
min_sub <- 2
models_dir <- file.path(out_dir, "MODELS", "method-edivisive")
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

cd <- readRDS(file.path(src, "corrected_depth.rds"))
b250 <- cd$bins
cells <- colnames(cd$corrected)
mb_start <- floor((b250$start - 1) / 1e6) * 1e6 + 1
key <- paste(b250$seqnames, mb_start)
groups <- split(seq_len(nrow(b250)), factor(key, levels = unique(key)))
n_sub <- lengths(groups)
groups <- groups[n_sub >= min_sub]
bins <- data.frame(seqnames = sub(" .*", "", names(groups)),
                   start = as.numeric(sub(".* ", "", names(groups))))
bins$end <- bins$start + 1e6 - 1
bins$n_sub_kept <- lengths(groups)
cat(sprintf("%d 1 Mb bins with >= %d kept 250 kb sub-bins (from %d kept 250 kb bins)\n",
            nrow(bins), min_sub, nrow(b250)))

depth1 <- t(sapply(groups, function(i) colMeans(cd$corrected[i, , drop = FALSE])))
colnames(depth1) <- cells
fc <- readRDS(file.path(repo, "results/K562_tree/aneufinder_fragment_gc/output/fragment_counts.rds"))
kb <- read.csv(file.path(repo, "results/K562_tree/aneufinder_fragment_gc/output/kept_bins.csv"))
raw_per_mb <- 4 * colMeans(fc$counts_gc[kb$kept, cells, drop = FALSE] / kb$mappable_frac[kb$kept])

gr <- GRanges(bins$seqnames, IRanges(bins$start, bins$end))
seqlevels(gr) <- unique(bins$seqnames)
seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]
end(gr) <- pmin(end(gr), seqlengths(gr)[as.character(seqnames(gr))])

run_cell <- function(j) {
  id <- cells[j]
  v <- depth1[, j] / mean(depth1[, j]) * raw_per_mb[j]
  b <- gr; b$counts <- v; b$mcounts <- v / 2; b$pcounts <- v / 2
  attr(b, "ID") <- id
  model <- suppressMessages(suppressWarnings(
    findCNVs(b, ID = id, method = "edivisive", R = 10, sig.lvl = 0.1, verbosity = 0)))
  save(model, file = file.path(models_dir, paste0(id, "_binsize_1e+06_fragmentGC_rtcorr.RData")))
  model$bins$copy.number
}
t0 <- Sys.time()
cn <- mclapply(seq_along(cells), function(j) tryCatch(run_cell(j), error = function(e) {
  message(cells[j], ": ", conditionMessage(e)); NULL }), mc.cores = ncpu)
ok <- !sapply(cn, is.null)
cn <- do.call(cbind, cn[ok]); colnames(cn) <- cells[ok]
cat(sprintf("segmented %d cells in %.1f min\n", ncol(cn), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(bins, file.path(out_dir, "bins_1mb.csv"), row.names = FALSE)
write.csv(cbind(bins[, c("seqnames", "start", "end")], cn), file.path(out_dir, "result.csv"))
saveRDS(list(bins = bins, depth = depth1), file.path(out_dir, "corrected_depth_1mb.rds"))
cat(sprintf("Wrote %s (%d bins x %d cells); copy number > 10 in %.4f of cell-bins\n",
            file.path(out_dir, "result.csv"), nrow(cn), ncol(cn), mean(cn > 10)))
