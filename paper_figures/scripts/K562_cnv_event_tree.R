# Per-cell K562 lineage tree from called CNV events.
#
# Built to replace cnv_sc_test.newick (NJ on 1 - Pearson r of whole copy-number
# profiles) in the SBS vs CNV tree comparison. Correlation ignores copy-number
# level and weighs every bin of a large event separately; here each cell is
# instead scored for discrete gain/loss events, and the tree is built the same
# way as the SBS tree (Hamming distance, NJ, rooted on an unaltered K562 tip).
#
# Steps:
#   1. keep cells at the dominant (near-triploid) ploidy, mean copy number in
#      [ploidy_min, ploidy_max); off-peak cells are mostly scaling errors or
#      doublets and would read as genome-wide gains or losses. To keep them,
#      widen the window and set K562_NORMALIZE_PLOIDY=1, which rescales every
#      cell to the median ploidy first
#   2. reference = per-bin modal copy number across those cells; bins whose
#      reference is 0 (unmappable) and chrY are dropped, and optionally bins
#      mostly covered by large repeats (K562_MASK_FRACTION, K562_MASK_MIN)
#   3. CNV call = a run of >= min_bins consecutive bins on one chromosome where
#      a cell is above (gain) or below (loss) the reference; shorter runs are
#      treated as noise
#   4. segment the genome at recurrent breakpoints (>= min_bp_cells cells
#      change state within +/- 1 bin), so a shared event counts once rather
#      than once per Mb; a cell is gained/lost in a segment when most of its
#      bins are
#   5. features = segment x {gain, loss}; drop features seen in < 2 cells
#   6. Hamming distance -> NJ, rooted on an all-zero K562 tip
#
# Input: an AneuFinder result.csv (1 Mb bins x cells; K562_CNV_RESULT, default
# the GC-corrected, blacklisted run from run_aneufinder_K562_sc_PolE_gc.sh; the
# original sc_PolE_novaseq run had neither, and its columns are mislabelled) and
# full_meta.csv for the set of cells. Ploidy is each cell's mean copy number in
# that result.csv, as in the notebook. K562_CNV_TAG suffixes the output files.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
})

input_dir <- Sys.getenv(
  "K562_SC_PROJECT",
  "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq"
)
result_csv <- Sys.getenv("K562_CNV_RESULT",
                         file.path(project_root, "results/K562_tree/aneufinder_gc/output/result.csv"))
tag <- Sys.getenv("K562_CNV_TAG", "_gc")
output_dir <- file.path(project_root, "results/K562_tree/sc_trees/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ploidy_min <- as.numeric(Sys.getenv("K562_PLOIDY_MIN", "2.5"))
ploidy_max <- as.numeric(Sys.getenv("K562_PLOIDY_MAX", "3.0"))
min_bins <- 3
min_cells <- 2
min_bp_cells <- 5

# ---- cells and copy-number matrix ------------------------------------------

meta <- read.csv(file.path(input_dir, "results/full_meta.csv"))
cnv <- read.csv(result_csv, row.names = 1, check.names = FALSE)
stopifnot(all(meta$X %in% colnames(cnv)))
bins <- cnv[, c("seqnames", "start", "end")]

ploidy <- colMeans(as.matrix(cnv[, meta$X]))
keep_cells <- meta$X[ploidy >= ploidy_min & ploidy < ploidy_max]
cat(sprintf("%d of %d cells at ploidy [%.2f, %.2f)\n",
            length(keep_cells), nrow(meta), ploidy_min, ploidy_max))
write.csv(data.frame(cell = meta$X, pop = meta$pop, ploidy = ploidy,
                     kept = meta$X %in% keep_cells),
          file.path(output_dir, paste0("cnv_cell_ploidy", tag, ".csv")), row.names = FALSE)
cn <- as.matrix(cnv[, keep_cells])

# Optionally rescale each cell to the median ploidy before calling events, so a
# uniform offset (AneuFinder picking the wrong overall scale, or a whole-genome
# doubling) is not read as a gain or loss in every bin
if (Sys.getenv("K562_NORMALIZE_PLOIDY", "0") == "1") {
  target <- median(ploidy[keep_cells])
  cn <- round(sweep(cn, 2, target / ploidy[keep_cells], "*"))
  storage.mode(cn) <- "integer"
  cat(sprintf("rescaled each cell to ploidy %.2f\n", target))
}

modal <- function(x) as.integer(names(which.max(table(x))))
reference <- apply(cn, 1, modal)

keep_bins <- reference > 0 & bins$seqnames != "chrY"
# Optionally drop bins dominated by large repeats (centromeres, gaps, satellites,
# segmental duplications; see make_hg38_repeat_mask.R)
mask_csv <- Sys.getenv("K562_MASK_FRACTION", "")
if (nzchar(mask_csv)) {
  mask_min <- as.numeric(Sys.getenv("K562_MASK_MIN", "0.25"))
  mf <- read.csv(mask_csv)
  mf <- mf$mask_frac[match(paste(bins$seqnames, bins$start), paste(mf$seqnames, mf$start))]
  stopifnot(!anyNA(mf))
  cat(sprintf("masking %d bins >= %.0f%% repeat (%d of them otherwise kept)\n",
              sum(mf >= mask_min), 100 * mask_min, sum(keep_bins & mf >= mask_min)))
  keep_bins <- keep_bins & mf < mask_min
}
bins <- bins[keep_bins, ]
cn <- cn[keep_bins, ]
reference <- reference[keep_bins]
cat(sprintf("%d bins after dropping chrY and reference-0 bins\n", nrow(bins)))

# ---- per-cell CNV calls ----------------------------------------------------

# sign of the deviation from reference, with runs shorter than min_bins zeroed
dev <- sign(cn - reference)
for (chr in unique(bins$seqnames)) {
  rows <- which(bins$seqnames == chr)
  dev[rows, ] <- apply(dev[rows, , drop = FALSE], 2, function(s) {
    r <- rle(s)
    r$values[r$values != 0 & r$lengths < min_bins] <- 0
    inverse.rle(r)
  })
}

# ---- shared segments -> event features --------------------------------------

# Splitting wherever any cell changes state leaves ~one segment per bin, since
# noise in any one of hundreds of cells moves the boundary. Instead segment at
# recurrent breakpoints: bin boundaries where >= min_bp_cells cells change
# state, counted within +/- 1 bin, keeping the local maximum of each cluster.
chr_start <- c(TRUE, bins$seqnames[-1] != bins$seqnames[-nrow(bins)])
bp <- rbind(FALSE, dev[-1, , drop = FALSE] != dev[-nrow(dev), , drop = FALSE])
bp[chr_start, ] <- FALSE
bp_count <- rowSums(bp)
window <- stats::filter(bp_count, rep(1, 3), sides = 2)
window[is.na(window)] <- bp_count[is.na(window)]
is_peak <- window >= min_bp_cells &
  bp_count == pmax(bp_count, c(0, head(bp_count, -1)), c(tail(bp_count, -1), 0)) &
  bp_count > 0
changed <- chr_start | is_peak
# fold segments shorter than min_bins into the preceding one on the chromosome
repeat {
  seg <- cumsum(changed)
  short <- which(changed & !chr_start & tabulate(seg)[seg] < min_bins)
  if (length(short) == 0) break
  changed[short[1]] <- FALSE
}
segment <- cumsum(changed)
first_bin <- which(changed)
cat(sprintf("%d recurrent breakpoints (>= %d cells) -> %d segments\n",
            sum(changed & !chr_start), min_bp_cells, max(segment)))

# a cell is gained (lost) in a segment when most of its bins are gained (lost)
seg_mean <- rowsum(dev, segment) / tabulate(segment)
seg_dev <- sign(seg_mean) * (abs(seg_mean) > 0.5)
seg_info <- data.frame(
  segment = seq_along(first_bin),
  chr = bins$seqnames[first_bin],
  start = bins$start[first_bin],
  end = tapply(bins$end, segment, max),
  n_bins = as.vector(table(segment)),
  reference_cn = reference[first_bin]
)

features <- rbind(
  (seg_dev > 0) * 1,
  (seg_dev < 0) * 1
)
feature_info <- rbind(
  mutate(seg_info, type = "gain", n_cells = rowSums(seg_dev > 0)),
  mutate(seg_info, type = "loss", n_cells = rowSums(seg_dev < 0))
)
informative <- feature_info$n_cells >= min_cells
features <- features[informative, , drop = FALSE]
feature_info <- feature_info[informative, ]
rownames(features) <- with(feature_info, sprintf("%s:%d-%d_%s", chr, start, end, type))
cat(sprintf("%d segments; %d informative gain/loss features (in >= %d cells)\n",
            nrow(seg_info), nrow(features), min_cells))

events_per_cell <- colSums(features)
cat("events per cell:\n")
print(summary(events_per_cell))
cat(sprintf("%d cells carry no informative event\n", sum(events_per_cell == 0)))

# ---- tree ------------------------------------------------------------------

x <- cbind(features, K562 = 0)
d <- dist(t(x), method = "manhattan") / nrow(x)
tree <- root(nj(d), outgroup = "K562", resolve.root = TRUE)
# NJ can return small negative branch lengths; clamp them, as is conventional
tree$edge.length <- pmax(tree$edge.length, 0)

write.tree(tree, file.path(output_dir, paste0("cnv_event_tree", tag, ".newick")))
write.csv(feature_info, file.path(output_dir, paste0("cnv_event_features", tag, ".csv")), row.names = FALSE)
write.csv(data.frame(cell = colnames(features), events = events_per_cell),
          file.path(output_dir, paste0("cnv_events_per_cell", tag, ".csv")), row.names = FALSE)
saveRDS(list(features = features, reference = reference, bins = bins),
        file.path(output_dir, paste0("cnv_event_matrix", tag, ".rds")))

cat("\nWrote tree to:", output_dir, "\n")
