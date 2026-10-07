#!/usr/bin/env Rscript

# Fragment-level GC-corrected read depth for K562 sc_PolE cells on the 250 kb
# grid (inputs from scripts/utils/fragment_gc_setup.R).
#
# Usage: fragment_gc_counts.R <cell_list.txt> <task_index> <cells_per_task>
# cell_list.txt has one BAM path per line; this task takes lines
# (task_index - 1) * cells_per_task + 1 ... task_index * cells_per_task.
#
# Per cell:
#   1. fragments: properly paired primary alignments with both mates
#      MAPQ >= 10, fragment = leftmost to rightmost mapped base (includes the
#      unsequenced insert), width <= 1000 bp, chr1-22/X/Y, not overlapping the
#      ENCODE blacklist; duplicates (same chromosome, start and end) removed
#   2. fragment GC from BSgenome.Hsapiens.UCSC.hg38 (the UCSC hg38 sequence the
#      cells were aligned to), as an integer percent 0-100
#   3. expected GC: 1e6 fragments built from the setup pool (uniquely mappable,
#      unmasked, copy-number-weighted start positions) with lengths drawn from
#      this cell's own fragment lengths
#   4. per-GC weight = expected / observed fraction (both smoothed over 5 GC
#      percentage points, pseudocount 0.5); GC values with < 20 fragments in
#      either take the nearest well-supported weight; weights scaled to mean 1
#      over the cell's fragments and clamped to [1/8, 8]
#   5. depth per bin = sum of fragment weights, by fragment midpoint (also the
#      unweighted count)
# Writes results/K562_tree/fragment_gc/cells/<cell>.rds.

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(Rsamtools)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

args <- commandArgs(trailingOnly = TRUE)
cell_list <- readLines(args[1])
task <- as.integer(args[2]); per_task <- as.integer(args[3])
todo <- cell_list[seq((task - 1) * per_task + 1, min(task * per_task, length(cell_list)))]

repo <- "/home/sanjay/SPC_genome"
gc_dir <- file.path(repo, "results/K562_tree/fragment_gc")
out_dir <- file.path(gc_dir, "cells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bins <- readRDS(file.path(gc_dir, "bins.rds"))
pool <- readRDS(file.path(gc_dir, "pool.rds"))
chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
bl <- read.table(file.path(repo, "results/K562_tree/aneufinder_gc/output/blacklist-hg38.raw.bed.gz"))
bl <- reduce(GRanges(bl$V1, IRanges(bl$V2 + 1, bl$V3)))
genome <- BSgenome.Hsapiens.UCSC.hg38

gc_percent <- function(gr) {
  as.integer(round(100 * letterFrequency(getSeq(genome, gr), "GC", as.prob = TRUE)[, 1]))
}

for (bam in todo) {
  cell <- sub("\\.bam$", "", basename(bam))
  out <- file.path(out_dir, paste0(cell, ".rds"))
  if (file.exists(out)) next
  t0 <- Sys.time()

  # 1. fragments, chromosome by chromosome
  flag <- scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE,
                      isSupplementaryAlignment = FALSE, isUnmappedQuery = FALSE)
  n_pairs <- 0
  frags <- do.call(c, lapply(chroms, function(ch) {
    p <- readGAlignmentPairs(bam, param = ScanBamParam(flag = flag, mapqFilter = 10,
                                                       which = GRanges(ch, IRanges(1, seqlengths(genome)[ch]))))
    n_pairs <<- n_pairs + length(p)
    g <- granges(p, on.discordant.seqnames = "drop")
    strand(g) <- "*"
    g <- g[width(g) <= 1000]
    unique(g)
  }))
  frags <- frags[!overlapsAny(frags, bl)]

  # 2-3. observed and expected fragment GC
  frag_gc <- gc_percent(frags)
  exp_idx <- sample(length(pool), 1e6, replace = TRUE)
  exp_frag <- GRanges(seqnames(pool)[exp_idx],
                      IRanges(start(pool)[exp_idx], width = sample(width(frags), 1e6, replace = TRUE)))
  exp_frag <- trim(suppressWarnings({ seqlengths(exp_frag) <- seqlengths(genome)[seqlevels(exp_frag)]; exp_frag }))
  exp_gc <- gc_percent(exp_frag)

  # 4. weights per GC percent
  obs_n <- tabulate(frag_gc + 1, 101); exp_n <- tabulate(exp_gc + 1, 101)
  sm <- function(x) as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2))
  obs_f <- sm(obs_n + 0.5); exp_f <- sm(exp_n + 0.5)
  obs_f[is.na(obs_f)] <- (obs_n + 0.5)[is.na(obs_f)]; exp_f[is.na(exp_f)] <- (exp_n + 0.5)[is.na(exp_f)]
  w <- (exp_f / sum(exp_f)) / (obs_f / sum(obs_f))
  ok <- obs_n >= 20 & exp_n >= 20
  w[!ok] <- w[ok][sapply(which(!ok), function(i) which.min(abs(which(ok) - i)))]
  w <- w / mean(w[frag_gc + 1])
  w <- pmin(pmax(w, 1 / 8), 8)

  # 5. per-bin depth by fragment midpoint
  mid <- resize(frags, 1, fix = "center")
  hit <- findOverlaps(mid, bins, select = "first")
  keep <- !is.na(hit)
  counts_raw <- tabulate(hit[keep], length(bins))
  counts_gc <- as.numeric(tapply(w[frag_gc[keep] + 1], factor(hit[keep], levels = seq_along(bins)), sum))
  counts_gc[is.na(counts_gc)] <- 0

  saveRDS(list(cell = cell, counts_raw = counts_raw, counts_gc = counts_gc,
               gc_table = data.frame(gc = 0:100, observed = obs_n, expected = exp_n, weight = w),
               n_pairs = n_pairs, n_fragments = length(frags),
               median_fragment = median(width(frags))), out)
  cat(sprintf("%s: %d pairs -> %d unique fragments, weight range %.2f-%.2f, %.1f min\n", cell, n_pairs,
              length(frags), min(w[ok]), max(w[ok]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
