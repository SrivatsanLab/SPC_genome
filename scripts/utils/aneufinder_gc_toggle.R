#!/usr/bin/env Rscript

# Controlled test of GC correction in AneuFinder: call copy number on the same
# cells, bins and blacklist, with GC correction as the only difference.
#
# run_aneufinder_K562_sc_PolE_gc.sh saved, for every cell, the blacklisted
# 1 Mb read counts before (binned/) and after (binned-GC/) GC correction. This
# runs AneuFinder's copy-number step (findCNVs, edivisive, R = 10, sig.lvl =
# 0.1, as in that run) on each, with the same random seed per cell in both arms
# (edivisive's significance test permutes). The GC arm should reproduce the
# run's own models; the check is written alongside.
#
# Output (ANEUFINDER_OUTPUT/gc_toggle/): per-cell ploidy in each arm and the
# bins x cells copy-number matrix for each arm.

suppressPackageStartupMessages({
    library(AneuFinder)
    library(parallel)
})

outputfolder <- Sys.getenv("ANEUFINDER_OUTPUT",
    "/home/sanjay/SPC_genome/results/K562_tree/aneufinder_gc/output")
out <- file.path(outputfolder, "gc_toggle")
dir.create(out, showWarnings = FALSE)
ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

files <- list.files(file.path(outputfolder, "binned"), "\\.RData$")
cells <- sub("\\.bam_.*", "", files)
cat(sprintf("%d cells, %d CPUs\n", length(cells), ncpu))

load_bins <- function(dir, cell) {
    f <- list.files(file.path(outputfolder, dir), paste0("^", cell, "\\.bam_"), full.names = TRUE)
    stopifnot(length(f) == 1)
    x <- get(load(f))
    if (is(x, "GRangesList")) x <- x[[1]]
    x
}
call_cn <- function(bins, cell, seed) {
    set.seed(seed)
    m <- suppressMessages(findCNVs(bins, ID = cell, method = "edivisive", R = 10, sig.lvl = 0.1,
                                   verbosity = 0))
    m$bins$copy.number
}

res <- mclapply(seq_along(cells), function(i) {
    cell <- cells[i]
    raw <- load_bins("binned", cell)
    gc <- load_bins("binned-GC", cell)
    list(cell = cell,
         raw = call_cn(raw, cell, seed = i),
         gc = call_cn(gc, cell, seed = i))
}, mc.cores = ncpu)

failed <- sapply(res, inherits, "try-error")
if (any(failed)) stop(sum(failed), " cells failed, e.g. ", as.character(res[[which(failed)[1]]]))

bins <- load_bins("binned", cells[1])
coords <- as.data.frame(bins)[, c("seqnames", "start", "end")]
cn_raw <- sapply(res, `[[`, "raw"); cn_gc <- sapply(res, `[[`, "gc")
colnames(cn_raw) <- colnames(cn_gc) <- sapply(res, `[[`, "cell")
write.csv(cbind(coords, cn_raw), gzfile(file.path(out, "cn_without_gc_correction.csv.gz")), row.names = FALSE)
write.csv(cbind(coords, cn_gc), gzfile(file.path(out, "cn_with_gc_correction.csv.gz")), row.names = FALSE)

ploidy <- data.frame(cell = colnames(cn_raw),
                     ploidy_without_gc = colMeans(cn_raw),
                     ploidy_with_gc = colMeans(cn_gc))
write.csv(ploidy, file.path(out, "ploidy_by_arm.csv"), row.names = FALSE)

# the GC arm should match the full run's own models
orig <- read.csv(file.path(outputfolder, "result.csv"), row.names = 1, check.names = FALSE)
same <- colMeans(as.matrix(orig[, colnames(cn_gc)]) == cn_gc)
writeLines(sprintf("GC arm vs full-run models: median %.4f of bins identical; %d of %d cells identical in every bin",
                   median(same), sum(same == 1), length(same)),
           file.path(out, "reproduction_check.txt"))
cat(readLines(file.path(out, "reproduction_check.txt")), sep = "\n")
cat("Wrote", out, "\n")
