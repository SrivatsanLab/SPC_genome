# Copy-number changes specific to one K562 population (P_0-P_5, successive
# splits of the culture, P_0 first).
#
# Input: a GC-corrected AneuFinder run of the 1000 sc_PolE cells
# (scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh [binsize]); K562_AF_DIR picks
# the run, default the 1 Mb one. Cells: on-peak (mean CN in [2.5, 3.0)) with
# Bhattacharyya >= the on-peak 10th percentile and spikiness <= the on-peak
# 95th percentile. Percentiles rather than fixed cut-offs, because both metrics
# shift with bin size (fewer reads per bin, more overlapping states): the 10th
# percentile is ~0.5 at 1 Mb, where a fixed 0.5 was used before.
#
# Caveat: hits are not adjusted for per-cell residual GC bias, which can mimic
# gains in GC-rich and losses in AT-rich regions (chr7q11, chr4q28; see
# docs/K562_aneufinder_gc_review.md). Check any hit against GC bias.
#
# Each bin's gain or loss is relative to the median CN of those cells. A hit is
# a run of >= min_bins consecutive bins on one chromosome where the change is
# carried by >= min_in of one population and <= max_out of the others. For each
# hit: carriers (cells changed in >= 80% of the run's bins), their population
# counts, a hypergeometric enrichment p (x6 populations), and where each
# carrier's changed segment around the hit starts and ends (breakpoints).
#
# Population labels come from cell_ploidy.csv; the per-cell SBS tree is not used.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(AneuFinder)
})

af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_gc/output"))
min_bins <- as.integer(Sys.getenv("K562_MIN_BINS", "3"))
min_in <- 0.15
max_out <- 0.05
ploidy_min <- 2.5
ploidy_max <- 3.0

cnv <- read.csv(file.path(af_dir, "result.csv"), row.names = 1, check.names = FALSE)
binsize <- cnv$end[1] - cnv$start[1] + 1
tag <- format(binsize, scientific = TRUE, trim = TRUE)
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_cnv_population_scan/", tag, "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- colnames(cnv)[-(1:3)]
cn <- as.matrix(cnv[, cells])
pop <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
pop <- setNames(pop$pop, pop$cell)[cells]

model_files <- list.files(file.path(af_dir, "MODELS/method-edivisive"), "\\.RData$", full.names = TRUE)
qc <- suppressWarnings(getQC(model_files))
rownames(qc) <- sub("\\.bam_binsize.*", "", basename(model_files))
qc <- qc[cells, ]

ploidy <- colMeans(cn)
on_peak <- ploidy >= ploidy_min & ploidy < ploidy_max
use <- on_peak & qc$bhattacharyya >= quantile(qc$bhattacharyya[on_peak], 0.10) &
  qc$spikiness <= quantile(qc$spikiness[on_peak], 0.95)
cat(sprintf("%s bins: %d cells (%s)\n", tag, sum(use),
            paste(names(table(pop[use])), table(pop[use]), sep = "=", collapse = " ")))

informative <- cnv$seqnames %in% c(paste0("chr", 1:22), "chrX") & rowMeans(cn == 0) <= 0.5
x <- cn[informative, use]
p <- pop[use]
pops <- sort(unique(p))
bins <- cnv[informative, c("seqnames", "start", "end")]
consensus <- apply(x, 1, median)

# ---- population median profiles --------------------------------------------

med_diff <- do.call(rbind, lapply(pops, function(q) {
  d <- apply(x[, p == q, drop = FALSE], 1, median) != consensus
  data.frame(pop = q, n_bins = sum(d), mb = sum(d) * binsize / 1e6)
}))
write.csv(med_diff, file.path(output_dir, "population_median_differences.csv"), row.names = FALSE)

# ---- scan -------------------------------------------------------------------

runs_of <- function(ok) {
  r <- rle(ok)
  ends <- cumsum(r$lengths)
  data.frame(from = ends - r$lengths + 1, to = ends)[r$values & r$lengths >= min_bins, ]
}

hits <- list()
for (kind in c("gain", "loss")) {
  changed <- if (kind == "gain") x > consensus else x < consensus
  for (q in pops) {
    f_in <- rowMeans(changed[, p == q, drop = FALSE])
    f_out <- rowMeans(changed[, p != q, drop = FALSE])
    for (ch in unique(bins$seqnames)) {
      idx <- which(bins$seqnames == ch)
      r <- runs_of(f_in[idx] >= min_in & f_out[idx] <= max_out)
      for (k in seq_len(nrow(r))) {
        b <- idx[r$from[k]:r$to[k]]
        carrier <- colMeans(changed[b, , drop = FALSE]) >= 0.8
        n <- sum(carrier); hit_in <- sum(carrier & p == q); n_q <- sum(p == q)
        pval <- phyper(hit_in - 1, n_q, length(p) - n_q, n, lower.tail = FALSE)
        hits[[length(hits) + 1]] <- data.frame(
          kind = kind, pop = q, chr = ch,
          start_mb = bins$start[min(b)] / 1e6, end_mb = bins$end[max(b)] / 1e6,
          f_in = mean(f_in[b]), f_out = mean(f_out[b]),
          carriers = n, carriers_in_pop = hit_in, frac_pop = hit_in / n_q,
          frac_other = (n - hit_in) / (length(p) - n_q),
          p_adj = min(1, pval * length(pops)),
          pops = paste(sapply(pops, function(z) sum(carrier & p == z)), collapse = "/"),
          carrier_cells = paste(colnames(x)[carrier], collapse = ","))
      }
    }
  }
}
hits <- if (length(hits)) do.call(rbind, hits) else data.frame()
write.csv(hits, file.path(output_dir, "population_hits.csv"), row.names = FALSE)

# ---- carrier breakpoints ----------------------------------------------------

# for each hit, the extent of each carrier's changed segment containing the hit
bp <- list()
for (h in seq_len(nrow(hits))) {
  idx <- which(bins$seqnames == hits$chr[h])
  core <- idx[bins$start[idx] / 1e6 >= hits$start_mb[h] & bins$end[idx] / 1e6 <= hits$end_mb[h]]
  mid <- core[ceiling(length(core) / 2)]
  for (cell in strsplit(hits$carrier_cells[h], ",")[[1]]) {
    v <- x[idx, cell]
    ch <- if (hits$kind[h] == "gain") v > consensus[idx] else v < consensus[idx]
    i <- match(mid, idx)
    if (!ch[i]) next
    lo <- i; while (lo > 1 && ch[lo - 1]) lo <- lo - 1
    hi <- i; while (hi < length(ch) && ch[hi + 1]) hi <- hi + 1
    bp[[length(bp) + 1]] <- data.frame(
      hit = h, chr = hits$chr[h], cell = cell, pop = p[cell],
      seg_start_mb = bins$start[idx[lo]] / 1e6, seg_end_mb = bins$end[idx[hi]] / 1e6,
      cn = median(v[lo:hi]))
  }
}
bp <- if (length(bp)) do.call(rbind, bp) else data.frame()
write.csv(bp, file.path(output_dir, "carrier_breakpoints.csv"), row.names = FALSE)

print(med_diff)
if (nrow(hits)) print(hits[, setdiff(names(hits), "carrier_cells")], digits = 3) else cat("no hits\n")
cat("Outputs written to", output_dir, "\n")
