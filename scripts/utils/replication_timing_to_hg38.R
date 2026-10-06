#!/usr/bin/env Rscript

# K562 replication timing on hg38, from the UW Repli-seq fractions downloaded by
# scripts/utils/download_replication_timing.sh.
#
# The six fractions share one hg19 window grid (mostly 1 kb). Per window: each
# fraction's percentage-normalised signal and an early-replication score
#   RT = sum(w_f * s_f) / sum(s_f),  w = 1, 0.8, 0.6, 0.4, 0.2, 0 for G1b, S1-S4, G2
# (1 = earliest). Windows are lifted to hg38 with the UCSC chain; windows that
# split or fail to lift are dropped, as are windows with no signal.
#
# Usage: replication_timing_to_hg38.R [rt_dir]
# Writes <rt_dir>/K562_replication_timing_hg38.bed.gz (chrom, start, end
# 0-based, G1..G2 signal, RT), sorted.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
rt_dir <- if (length(args)) args[1] else "data/reference/replication_timing"
fractions <- c(G1 = "ENCFF001GRX", S1 = "ENCFF001GSF", S2 = "ENCFF001GSJ",
               S3 = "ENCFF001GSN", S4 = "ENCFF001GSP", G2 = "ENCFF001GSB")
weights <- c(G1 = 1, S1 = 0.8, S2 = 0.6, S3 = 0.4, S4 = 0.2, G2 = 0)
chroms <- c(paste0("chr", 1:22), "chrX")

tracks <- lapply(names(fractions), function(fr) {
  bw <- import(file.path(rt_dir, sprintf("K562_%s_%s.bigWig", fr, fractions[[fr]])))
  bw[seqnames(bw) %in% chroms]
})
names(tracks) <- names(fractions)

# common hg19 window grid: windows present, with identical coordinates, in all six
key <- function(t) paste(seqnames(t), start(t), end(t))
grid <- granges(tracks[[1]])
sig <- sapply(tracks, function(t) t$score[match(key(grid), key(t))])
keep <- rowSums(sig, na.rm = TRUE) > 0 & !apply(is.na(sig), 1, any)
grid <- grid[keep]; sig <- sig[keep, , drop = FALSE]
mcols(grid) <- as.data.frame(sig)
grid$RT <- as.vector(sig %*% weights) / rowSums(sig)
cat(sprintf("hg19: %d windows (%.2f Gb) with signal in all six fractions\n",
            length(grid), sum(as.numeric(width(grid))) / 1e9))

chain_gz <- file.path(rt_dir, "hg19ToHg38.over.chain.gz")
chain_file <- file.path(tempdir(), "hg19ToHg38.over.chain")
R.utils::gunzip(chain_gz, chain_file, remove = FALSE, overwrite = TRUE)
lifted <- liftOver(grid, import.chain(chain_file))
one <- lengths(lifted) == 1
hg38 <- unlist(lifted[one])
hg38 <- sort(hg38[seqnames(hg38) %in% chroms])
cat(sprintf("hg38: %d windows lifted one-to-one (%.1f%% of hg19 windows)\n",
            length(hg38), 100 * mean(one)))

out <- data.frame(chrom = as.character(seqnames(hg38)), start = start(hg38) - 1, end = end(hg38),
                  round(as.data.frame(mcols(hg38))[, c(names(fractions), "RT")], 4))
gz <- gzfile(file.path(rt_dir, "K562_replication_timing_hg38.bed.gz"), "w")
write.table(out, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(gz)
cat("Wrote", file.path(rt_dir, "K562_replication_timing_hg38.bed.gz"), "\n")
