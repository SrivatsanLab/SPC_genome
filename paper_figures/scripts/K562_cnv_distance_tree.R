# Distance-based copy-number tree of the K562 sc_PolE cells from the
# replication-corrected calls (scripts/utils/replication_regress_segment.R:
# fragment-level GC weights, repeat/unmappable/low-mappability bins masked,
# cell-cycle replication component removed, smoothed, AneuFinder edivisive).
#
# Cells: all except the off-peak artefacts (K562_offpeak_ploidy.R). Each cell's
# copy number is rescaled to the median ploidy of the non-artefact cells and
# rounded, so a whole-genome scale offset is not counted as distance (as in
# K562_cnv_event_tree.R). Distance = mean absolute copy-number difference over
# bins. Neighbour joining (ape::nj) with the K562 consensus (per-bin median of
# on-peak cells) added as an outgroup and used as the root.
# K562_AF_DIR picks another run with the same result.csv layout; K562_CELLS=G1G2
# keeps only cells called G1/G2 (K562_replication_timing_coverage.R), which have
# no replication component, and writes to <run>_G1G2; K562_CELLS=G1G2strict also
# drops G1/G2 cells with a fitted replication amplitude a >= 0.1
# (replication_fit.csv; likely at the edges of S), writing to <run>_G1G2strict. K562_MIN_EVENT_MB (default
# 0, off) keeps only copy-number events of at least that size: per cell, runs of
# bins on one chromosome at the same deviation from the consensus that span
# less than the minimum are set back to the consensus before distances are
# computed (output suffix _ev<size>Mb).

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(ape))

af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_fragment_gc_rtcorr/output"))
cell_set <- Sys.getenv("K562_CELLS", "all")
min_event_mb <- as.numeric(Sys.getenv("K562_MIN_EVENT_MB", "0"))
run_tag <- paste0(basename(dirname(af_dir)), if (cell_set != "all") paste0("_", cell_set) else "",
                  if (min_event_mb > 0) sprintf("_ev%gMb", min_event_mb) else "")
output_dir <- file.path(project_root, "paper_figures/output/K562_cnv_distance_tree", run_tag)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cnv <- read.csv(file.path(af_dir, "result.csv"), row.names = 1, check.names = FALSE)
cells <- setdiff(colnames(cnv), c("seqnames", "start", "end"))
metrics <- read.csv(file.path(project_root, "paper_figures/output/K562_offpeak_ploidy/per_cell_metrics.csv"))
cells <- setdiff(cells, metrics$cell[grepl("^artefact", metrics$class)])
if (cell_set %in% c("G1G2", "G1G2strict")) {
  ph <- read.csv(file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/per_cell_phase.csv"))
  cells <- intersect(cells, ph$cell[ph$phase == "G1/G2"])
}
if (cell_set == "G1G2strict") {
  rf <- read.csv(file.path(project_root, "results/K562_tree/aneufinder_fragment_gc_rtcorr/output/replication_fit.csv"))
  cells <- intersect(cells, rf$cell[rf$a < 0.1])
}
auto <- cnv$seqnames %in% c(paste0("chr", 1:22), "chrX")
cn <- as.matrix(cnv[auto, cells])

ploidy <- colMeans(cn)
cn <- round(sweep(cn, 2, median(ploidy) / ploidy, "*"))
on_peak <- cells %in% metrics$cell[metrics$class == "on-peak"]
consensus <- apply(cn[, on_peak], 1, median)

if (min_event_mb > 0) {
  chrom <- cnv$seqnames[auto]; bstart <- cnv$start[auto]; bend <- cnv$end[auto]
  n_dropped <- 0
  for (j in seq_len(ncol(cn))) {
    dev <- cn[, j] - consensus
    for (ch in unique(chrom)) {
      i <- which(chrom == ch)
      r <- rle(dev[i]); e <- cumsum(r$lengths); st <- e - r$lengths + 1
      span <- (bend[i][e] - bstart[i][st] + 1) / 1e6
      short <- which(r$values != 0 & span < min_event_mb)
      for (q in short) cn[i[st[q]:e[q]], j] <- consensus[i[st[q]:e[q]]]
      n_dropped <- n_dropped + length(short)
    }
  }
  cat(sprintf("dropped %d events < %g Mb (%.1f per cell)\n", n_dropped, min_event_mb, n_dropped / ncol(cn)))
}
x <- cbind(cn, K562_consensus = consensus)

d <- dist(t(x), method = "manhattan") / nrow(x)   # mean |CN difference| per bin
tree <- root(nj(d), outgroup = "K562_consensus", resolve.root = TRUE)
tree$edge.length[tree$edge.length < 0] <- 0
write.tree(tree, file.path(output_dir, "cnv_distance_tree.newick"))
saveRDS(as.matrix(d), file.path(output_dir, "cnv_distance_matrix.rds"))
cat(sprintf("%s: %d cells, %d bins; tree written to %s\n", run_tag, length(cells), nrow(cn),
            file.path(output_dir, "cnv_distance_tree.newick")))
cat(sprintf("pairwise distance (mean |CN difference| per bin): median %.3f, IQR %.3f-%.3f\n",
            median(d), quantile(d, 0.25), quantile(d, 0.75)))
