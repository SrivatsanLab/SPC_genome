#!/usr/bin/env Rscript

# Re-anchor each K562 cell's copy-number scale to the K562 karyotype.
#
# AneuFinder's edivisive step scales each cell's depth by the multiplier
# (1.5-6) that puts its segments closest to whole numbers. A profile and the
# same profile x 4/3, x 3/2 or x 2 fit about equally well, so in some cells the
# dominant K562 state (3 copies) is called 2, 4, 5 or 6 and every state shifts
# with it (on the replication-corrected 1 Mb calls: 13 good cells at 6, 8 at 4,
# 3 at 5, 29 at 2). Coverage cannot detect true whole-genome changes, so the
# scale is fixed from the karyotype instead:
#   baseline = the cell's median depth over bins where the consensus is 3
#   copy number of each AneuFinder segment = round(3 * segment mean depth / baseline)
# Segment boundaries are AneuFinder's; only the integer levels change.
# Consensus = per-bin median copy number of good cells (not artefact, not
# noisy) already at the dominant ploidy [2.5, 3.2).
#
# Input/output: K562_AF_DIR (default the replication-corrected 1 Mb run) with
# corrected_depth_1mb.rds; writes <run>_anchored/output with MODELS,
# result.csv and anchoring.csv (per-cell scale change).

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
})

repo <- "/home/sanjay/SPC_genome"
af_dir <- Sys.getenv("K562_AF_DIR", file.path(repo, "results/K562_tree/aneufinder_fragment_gc_rtcorr_1mb/output"))
out_dir <- file.path(dirname(af_dir) |> paste0("_anchored"), "output")
models_out <- file.path(out_dir, "MODELS", "method-edivisive")
dir.create(models_out, recursive = TRUE, showWarnings = FALSE)

cnv <- read.csv(file.path(af_dir, "result.csv"), row.names = 1, check.names = FALSE)
cells <- setdiff(colnames(cnv), c("seqnames", "start", "end"))
x <- readRDS(file.path(af_dir, "corrected_depth_1mb.rds"))
stopifnot(all(paste(x$bins$seqnames, x$bins$start) == paste(cnv$seqnames, cnv$start)))
metrics <- read.csv(file.path(repo, "paper_figures/output/K562_offpeak_ploidy/per_cell_metrics.csv"))
noise <- read.csv(file.path(repo, "paper_figures/output/K562_cell_noise/per_cell_noise.csv"))
good <- !(cells %in% c(metrics$cell[grepl("^artefact", metrics$class)], noise$cell[noise$noisy]))
auto <- cnv$seqnames != "chrY"
cn <- as.matrix(cnv[, cells])
ploidy <- colMeans(cn[auto, ])
consensus <- apply(cn[, good & ploidy >= 2.5 & ploidy < 3.2], 1, median)
base_bins <- auto & consensus == 3
cat(sprintf("consensus from %d cells; %d bins at consensus 3\n", sum(good & ploidy >= 2.5 & ploidy < 3.2), sum(base_bins)))

model_files <- list.files(file.path(af_dir, "MODELS", "method-edivisive"), full.names = TRUE)
names(model_files) <- sub("_binsize.*", "", basename(model_files))
key <- paste(cnv$seqnames, cnv$start)

anchor <- do.call(rbind, lapply(cells, function(cl) {
  e <- new.env(); load(model_files[[cl]], envir = e); model <- e$model
  b <- model$bins
  dep <- x$depth[match(paste(seqnames(b), start(b)), key), cl]
  baseline <- median(dep[match(paste(seqnames(b), start(b)), key) %in% which(base_bins)])
  old_cn <- b$copy.number
  # segments = runs of AneuFinder's own calls along each chromosome
  chr <- as.character(seqnames(b))
  seg <- cumsum(c(TRUE, old_cn[-1] != old_cn[-length(old_cn)] | chr[-1] != chr[-length(chr)]))
  new_cn <- round(3 * ave(dep, seg, FUN = mean) / baseline)
  new_cn <- pmax(new_cn, 0)
  levs <- paste0(sort(unique(new_cn)), "-somy")
  b$copy.number <- new_cn
  b$state <- factor(paste0(new_cn, "-somy"), levels = levs)
  model$bins <- b
  segs <- unlist(reduce(split(granges(b), paste(chr, seg, sep = "_"))))
  first <- match(start(segs), start(b))
  first <- sapply(seq_along(segs), function(i) which(chr == as.character(seqnames(segs))[i] & start(b) == start(segs)[i])[1])
  segs$copy.number <- new_cn[first]
  segs$state <- factor(paste0(segs$copy.number, "-somy"), levels = levs)
  segs <- sort(segs)
  model$segments <- segs
  save(model, file = file.path(models_out, basename(model_files[[cl]])))
  data.frame(cell = cl, old_cn_at_consensus3 = median(old_cn[match(paste(seqnames(b), start(b)), key) %in% which(base_bins)]),
             new_cn_at_consensus3 = median(new_cn[match(paste(seqnames(b), start(b)), key) %in% which(base_bins)]),
             old_ploidy = mean(old_cn[chr != "chrY"]), new_ploidy = mean(new_cn[chr != "chrY"]),
             changed_bins = mean(new_cn != old_cn), good = cl %in% cells[good])
}))
write.csv(anchor, file.path(out_dir, "anchoring.csv"), row.names = FALSE)

new_cn <- sapply(cells, function(cl) { e <- new.env(); load(file.path(models_out, basename(model_files[[cl]])), envir = e)
  e$model$bins$copy.number[match(key, paste(seqnames(e$model$bins), start(e$model$bins)))] })
write.csv(cbind(cnv[, c("seqnames", "start", "end")], new_cn), file.path(out_dir, "result.csv"))
file.copy(file.path(af_dir, "corrected_depth_1mb.rds"), out_dir, overwrite = TRUE)

cat("called copy number at consensus-3 bins, before -> after (good cells):\n")
print(table(before = anchor$old_cn_at_consensus3[anchor$good], after = anchor$new_cn_at_consensus3[anchor$good]))
cat(sprintf("cells with any bin changed: %d; median fraction of bins changed among those: %.3f\n",
            sum(anchor$changed_bins > 0), median(anchor$changed_bins[anchor$changed_bins > 0])))
cat("Wrote", out_dir, "\n")
