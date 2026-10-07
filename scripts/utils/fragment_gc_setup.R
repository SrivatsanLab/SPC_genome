#!/usr/bin/env Rscript

# One-off inputs for fragment-level GC correction of the K562 sc_PolE cells
# (scripts/utils/fragment_gc_counts.R). Writes to results/K562_tree/fragment_gc/:
#   bins.rds  - the 250 kb grid of the AneuFinder runs (chr1-22, X, Y) with
#               each bin's fraction uniquely mappable (Umap k100,
#               results/K562_tree/masks/k100.umap.bed.gz), covered by the hg38
#               repeat mask (centromeres, gaps, Satellite repeats, segmental
#               duplications) and by N gaps + ENCODE blacklist
#   pool.rds  - candidate fragment start positions for the expected fragment GC
#               distribution: drawn uniformly from uniquely mappable sequence
#               outside the repeat mask and blacklist, each position weighted by
#               the K562 consensus copy number at that point (median of the
#               on-peak cells in the GC-corrected 1 Mb run), so the karyotype is
#               not read as GC bias

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

repo <- "/home/sanjay/SPC_genome"
out_dir <- file.path(repo, "results/K562_tree/fragment_gc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
n_pool <- 5e6
set.seed(1)
chroms <- c(paste0("chr", 1:22), "chrX", "chrY")

read_bed <- function(f, skip = 0) {
  b <- read.table(f, skip = skip, sep = "\t", quote = "", comment.char = "")
  b <- b[b$V1 %in% chroms, ]
  reduce(GRanges(b$V1, IRanges(b$V2 + 1, b$V3)))
}
covered_frac <- function(gr, bed) {
  ov <- findOverlaps(gr, bed)
  bp <- tapply(width(pintersect(gr[queryHits(ov)], bed[subjectHits(ov)])), queryHits(ov), sum)
  f <- numeric(length(gr)); f[as.integer(names(bp))] <- bp / width(gr)[as.integer(names(bp))]
  f
}

# ---- 250 kb grid, identical to the AneuFinder runs --------------------------

e <- new.env()
load(list.files(file.path(repo, "results/K562_tree/aneufinder_gc_2.5e5/output/binned-GC"),
                full.names = TRUE)[1], envir = e)
bins <- get(ls(e)[1], envir = e)[[1]]
bins <- bins[seqnames(bins) %in% chroms]
mcols(bins) <- mcols(bins)["GC"]

umap <- read_bed(file.path(repo, "results/K562_tree/masks/k100.umap.bed.gz"), skip = 1)
rep_mask <- read_bed(file.path(repo, "results/K562_tree/masks/hg38_repeat_mask.bed"))
blacklist <- read_bed(file.path(repo, "results/K562_tree/aneufinder_gc/output/blacklist-hg38.raw.bed.gz"))
n_frac <- letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, bins), "N", as.prob = TRUE)[, 1]

bins$mappable_frac <- covered_frac(bins, umap)
bins$repeat_frac <- covered_frac(bins, rep_mask)
bins$unmappable_frac <- pmin(1, n_frac + covered_frac(bins, blacklist))
saveRDS(bins, file.path(out_dir, "bins.rds"))
cat(sprintf("%d bins; median mappable fraction %.3f; %d bins < 0.9 mappable\n",
            length(bins), median(bins$mappable_frac), sum(bins$mappable_frac < 0.9)))

# ---- pool of fragment start positions for the expected GC ------------------

allowed <- setdiff(setdiff(umap, rep_mask), blacklist)
cnv <- read.csv(file.path(repo, "results/K562_tree/aneufinder_gc/output/result.csv"),
                row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(repo, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
cons <- apply(as.matrix(cnv[, meta$cell[meta$kept]]), 1, median)
cons_gr <- GRanges(cnv$seqnames, IRanges(cnv$start, cnv$end), cn = cons)
cons_gr <- cons_gr[cons_gr$cn > 0]
allowed <- unlist(tile(allowed, width = 1e6))   # pieces short enough to take one copy number
hit <- findOverlaps(resize(allowed, 1, fix = "center"), cons_gr, select = "first")
allowed <- allowed[!is.na(hit)]
w <- as.numeric(width(allowed)) * cons_gr$cn[hit[!is.na(hit)]]
pick <- sample(length(allowed), n_pool, replace = TRUE, prob = w)
pos <- start(allowed)[pick] + floor(runif(n_pool) * width(allowed)[pick])
pool <- GRanges(seqnames(allowed)[pick], IRanges(pos, width = 1))
saveRDS(pool, file.path(out_dir, "pool.rds"))
cat(sprintf("pool: %d positions from %.2f Gb allowed sequence\n", length(pool), sum(as.numeric(width(allowed))) / 1e9))
