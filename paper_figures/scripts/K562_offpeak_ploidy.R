# Which K562 cells off the dominant ploidy are artefacts, and which carry real
# copy-number variation?
#
# Input: the GC-corrected, blacklisted AneuFinder run of the 1000 sc_PolE cells
# (scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh): result.csv, the per-cell
# models (QC) and binned-GC read counts (depth before copy-number calling).
#
# Off-peak = mean copy number outside [ploidy_min, ploidy_max), as in
# K562_cnv_event_tree.R. Every test works on GC-corrected relative depth over
# informative autosomal bins (not CN0 in most cells, not amplicons at median
# CN >= 6), compared with the median depth profile of the on-peak cells:
#   r_depth     Pearson r with that consensus profile
#   noise       median |bin-to-bin change| of log2(depth / consensus)
#   sdev        fraction of 10-bin windows where log2(depth / consensus)
#               averages beyond +/-0.35 (large structured deviation)
#   lattice_err weighted distance of segment depths (runs of constant CN,
#               >= 3 bins) from whole-number copy states at the best scale for
#               mean ploidy 1.5-6
#   r_gc        r of log2(depth / consensus) with bin GC (residual GC bias)
# Cut-offs are the 95th (5th for r_depth) percentile of the on-peak cells.
#
# Classes: artefact if noisy or not on whole-number levels; "scale only" if
# not structurally deviant (depth matches the consensus, ploidy is mis-scaled);
# otherwise a real-looking copy-number change.
#
# Also: per-cell Pearson r with the K562 reference profile and a bin-shuffle
# permutation p, as in notebooks/K562_tree.ipynb.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
af_dir <- file.path(project_root, "results/K562_tree/aneufinder_gc/output")
output_dir <- file.path(project_root, "paper_figures/output/K562_offpeak_ploidy/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ploidy_min <- 2.5
ploidy_max <- 3.0
set.seed(1)

# ---- copy number, QC and depth ----------------------------------------------

cnv <- read.csv(file.path(af_dir, "result.csv"), row.names = 1, check.names = FALSE)
cells <- colnames(cnv)[-(1:3)]
cn_all <- as.matrix(cnv[, cells])
ploidy <- colMeans(cn_all)
kept <- ploidy >= ploidy_min & ploidy < ploidy_max
cat(sprintf("%d of %d cells at ploidy [%.2f, %.2f)\n", sum(kept), length(cells),
            ploidy_min, ploidy_max))

cell_of <- function(f) sub("\\.bam_binsize.*", "", basename(f))
model_files <- list.files(file.path(af_dir, "MODELS/method-edivisive"), "\\.RData$",
                          full.names = TRUE)
qc <- suppressWarnings(getQC(model_files))
rownames(qc) <- cell_of(model_files)
qc <- qc[cells, ]

load_gr <- function(f) {
  e <- new.env()
  load(f, envir = e)
  g <- get(ls(e)[1], envir = e)
  if (!is(g, "GRanges")) g <- g[[1]]
  g
}
binned_files <- list.files(file.path(af_dir, "binned-GC"), "\\.RData$", full.names = TRUE)
bins <- load_gr(binned_files[1])
stopifnot(all(start(bins) == cnv$start), all(as.character(seqnames(bins)) == cnv$seqnames))
depth_all <- sapply(binned_files, function(f) mcols(load_gr(f))$counts)
colnames(depth_all) <- cell_of(binned_files)
depth_all <- depth_all[, cells]

autosome <- cnv$seqnames %in% paste0("chr", 1:22)
informative <- autosome &
  rowMeans(cn_all == 0) <= 0.5 &
  apply(cn_all, 1, median) < 6 &
  apply(depth_all, 1, median) > 0
cat(sprintf("%d informative autosomal bins\n", sum(informative)))

chrom <- cnv$seqnames[informative]
cn <- cn_all[informative, ]
rel <- sweep(depth_all[informative, ], 2, colMeans(depth_all[informative, ]), "/")
consensus <- apply(rel[, kept], 1, median)
resid <- log2((rel + 1e-3) / (consensus + 1e-3))
gc <- mcols(bins)$GC[informative]

# ---- per-cell metrics -------------------------------------------------------

window_dev <- function(v) {
  w <- unlist(lapply(split(v, chrom), function(x) {
    if (length(x) < 10) return(NULL)
    stats::filter(x, rep(1 / 10, 10), sides = 1)[10:length(x)]
  }))
  mean(abs(w) > 0.35)
}

lattice_err <- function(i) {
  seg_mean <- c(); seg_len <- c()
  for (ch in unique(chrom)) {
    idx <- which(chrom == ch)
    r <- rle(cn[idx, i])
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    long <- r$lengths >= 3
    seg_mean <- c(seg_mean, mapply(function(a, b) mean(rel[idx[a:b], i]),
                                   starts[long], ends[long]))
    seg_len <- c(seg_len, r$lengths[long])
  }
  min(sapply(seq(1.5, 6, by = 0.01), function(p) {
    y <- seg_mean * p
    weighted.mean(abs(y - round(y)), seg_len)
  }))
}

metrics <- data.frame(
  cell = cells,
  ploidy = ploidy,
  kept = kept,
  spikiness = qc$spikiness,
  bhattacharyya = qc$bhattacharyya,
  total_reads = qc$total.read.count,
  r_depth = apply(rel, 2, cor, consensus),
  noise = apply(abs(diff(resid)), 2, median),
  sdev = apply(resid, 2, window_dev),
  lattice_err = sapply(seq_along(cells), lattice_err),
  r_gc = apply(resid, 2, cor, gc)
)

on <- metrics[metrics$kept, ]
cut_noise <- quantile(on$noise, 0.95)
cut_sdev <- quantile(on$sdev, 0.95)
cut_lattice <- quantile(on$lattice_err, 0.95)

metrics$class <- with(metrics, ifelse(
  kept, "on-peak",
  ifelse(noise > cut_noise | lattice_err > cut_lattice, "artefact: noisy or off-integer",
  ifelse(sdev <= cut_sdev, "artefact: scale only", "real")))))
metrics$side <- ifelse(metrics$kept, "on", ifelse(metrics$ploidy < ploidy_min, "low", "high"))

write.csv(metrics, file.path(output_dir, "per_cell_metrics.csv"), row.names = FALSE)
class_table <- as.data.frame.matrix(table(metrics$class, metrics$side))
write.csv(class_table, file.path(output_dir, "class_counts.csv"))
print(class_table)

summ <- aggregate(cbind(spikiness, r_depth, noise, sdev, lattice_err, r_gc) ~ side,
                  metrics, median)
write.csv(summ, file.path(output_dir, "metric_medians_by_side.csv"), row.names = FALSE)
print(summ, digits = 3)

# ---- agreement with the K562 reference profile ------------------------------

ref <- read.table(file.path(sc_dir, "K562_ref_ploidy_hg38.bed"),
                  col.names = c("chr", "start", "end", "cn"))
ref_gr <- GRanges(ref$chr, IRanges(ref$start + 1, ref$end))
mid <- GRanges(cnv$seqnames, IRanges((cnv$start + cnv$end) %/% 2, width = 1))
hit <- findOverlaps(mid, ref_gr, select = "first")
ref_bins <- !is.na(hit) & cnv$seqnames != "chrY"
ref_cn <- ref$cn[hit[ref_bins]]
x <- cn_all[ref_bins, ]
rho <- apply(x, 2, function(v) if (sd(v) > 0) cor(v, ref_cn) else NA)
metrics$ref_pearson <- rho

ref_rows <- lapply(list(all = rep(TRUE, length(cells)), on_peak = kept), function(m) {
  obs <- mean(rho[m], na.rm = TRUE)
  null <- replicate(1000, mean(cor(x[, m], sample(ref_cn)), na.rm = TRUE))
  data.frame(n_cells = sum(m), bins = sum(ref_bins), mean_r = obs,
             sd_r = sd(rho[m], na.rm = TRUE), median_r = median(rho[m], na.rm = TRUE),
             null_max = max(null), empirical_p = mean(null >= obs))
})
ref_summary <- cbind(cells = names(ref_rows), do.call(rbind, ref_rows))
write.csv(ref_summary, file.path(output_dir, "reference_correlation.csv"), row.names = FALSE)
print(ref_summary, digits = 3)
write.csv(metrics, file.path(output_dir, "per_cell_metrics.csv"), row.names = FALSE)

cat("\nOutputs written to", output_dir, "\n")
