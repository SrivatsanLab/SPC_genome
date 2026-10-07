# Small clones: are highly similar pairs of K562 cells concordant with the
# populations (P_0-P_5, successive splits of the culture)?
#
# Input: the GC-corrected 1 Mb AneuFinder run (K562_AF_DIR), its binned-GC read
# counts, the unmappable-bin mask (K562_cnv_reference_correlation.R) and the
# per-cell GC bias from K562_gc_toggle_ploidy.R. Cells: on-peak with
# Bhattacharyya >= on-peak 10th percentile and spikiness <= 95th.
#
# Similarity between two cells = Pearson r of their GC-adjusted depth
# profiles: log2 relative depth per informative bin (autosomes + X, unmasked,
# median CN < 6), with each bin's dependence on per-cell GC bias removed by a
# cubic regression across cells. (Similarity on copy-number deviations is also
# reported; it tracks GC-bias similarity, so it is not used for inference.)
#
# Tests: same-population fraction among the top-k pairs and nearest neighbours;
# barcode-well effects (segment 4 of the cell barcode is the round-1 well and
# encodes population; segments 1-3 are later rounds shared across
# populations); a well-preserving permutation (whole round-1 wells shuffled
# between populations) of the same-population, different-well excess.
# Groups: connected components of pairs in the top 0.1%, with population and
# well make-up, QC, and shared regions (3-bin smoothed mean residual beyond
# +/-0.25 over >= 3 bins).

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
})

af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_gc/output"))
output_dir <- file.path(project_root, "paper_figures/output/K562_cnv_cell_similarity/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(1)
n_perm <- 2000

# ---- cells, bins, depth -----------------------------------------------------

cnv <- read.csv(file.path(af_dir, "result.csv"), row.names = 1, check.names = FALSE)
all_cells <- colnames(cnv)[-(1:3)]
tag <- format(cnv$end[1] - cnv$start[1] + 1, scientific = TRUE, trim = TRUE)
cn_all <- as.matrix(cnv[, all_cells])
ploidy <- colMeans(cn_all)
on_peak <- ploidy >= 2.5 & ploidy < 3.0

cell_of <- function(f) sub("\\.bam_binsize.*", "", basename(f))
model_files <- list.files(file.path(af_dir, "MODELS/method-edivisive"), "\\.RData$", full.names = TRUE)
qc <- suppressWarnings(getQC(model_files))
rownames(qc) <- cell_of(model_files)
qc <- qc[all_cells, ]
use <- on_peak & qc$bhattacharyya >= quantile(qc$bhattacharyya[on_peak], 0.10) &
  qc$spikiness <= quantile(qc$spikiness[on_peak], 0.95)
cells <- all_cells[use]
n <- length(cells)

meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
pop <- setNames(meta$pop, meta$cell)[cells]
gcb <- read.csv(file.path(project_root, "paper_figures/output/K562_gc_toggle_ploidy/ploidy_with_without_gc.csv"))
gc_bias <- setNames(gcb$gc_bias, gcb$cell)[cells]
seg <- do.call(rbind, regmatches(cells, regexec("^(\\w{8})AGGA(\\w{8})ACTC(\\w{8})AAGG(\\w+)$", cells)))[, 2:5]
well <- seg[, 4]

mask <- read.csv(file.path(project_root, "results/K562_tree/masks", paste0("unmappable_bins_", tag, ".csv")))
cn <- cn_all[, cells]
bins <- cnv$seqnames %in% c(paste0("chr", 1:22), "chrX") & !mask$masked & apply(cn, 1, median) < 6

load_gr <- function(f) {
  e <- new.env(); load(f, envir = e); g <- get(ls(e)[1], envir = e)
  if (!is(g, "GRanges")) g <- g[[1]]
  g
}
binned_files <- list.files(file.path(af_dir, "binned-GC"), "\\.RData$", full.names = TRUE)
names(binned_files) <- cell_of(binned_files)
depth <- sapply(binned_files[cells], function(f) mcols(load_gr(f))$counts)[bins, ]
rel <- log2(sweep(depth, 2, colMeans(depth), "/") + 1e-3)
X <- cbind(1, gc_bias, gc_bias^2, gc_bias^3)
resid <- t(t(rel) - X %*% qr.solve(X, t(rel)))   # per-bin GC-bias trend removed
dev_cn <- cn[bins, ] - apply(cn[bins, ], 1, median)

sim <- cor(resid)
sim_cn <- cor(dev_cn)
chrom <- cnv$seqnames[bins]
start <- cnv$start[bins]
cat(sprintf("%d cells, %d bins\n", n, sum(bins)))

# ---- population concordance of similar pairs ---------------------------------

ut <- upper.tri(sim)
same_pop <- outer(pop, pop, "==")[ut]
same_well <- outer(well, well, "==")[ut]
expected <- mean(same_pop)
gc_close <- -abs(outer(gc_bias, gc_bias, "-"))[ut]

concord <- do.call(rbind, lapply(c("depth", "cn"), function(which_sim) {
  S <- if (which_sim == "depth") sim else sim_cn
  v <- S[ut]; o <- order(v, decreasing = TRUE)
  S2 <- S; diag(S2) <- -Inf; nn <- apply(S2, 1, which.max)
  nn_obs <- mean(pop[nn] == pop)
  nn_null <- replicate(n_perm, { p <- sample(pop); mean(p[nn] == p) })
  data.frame(
    similarity = which_sim,
    measure = c(paste0("top ", c(25, 50, 100, 250, 1000, 5000), " pairs"), "nearest neighbour"),
    same_population = c(sapply(c(25, 50, 100, 250, 1000, 5000), function(k) mean(same_pop[o[1:k]])), nn_obs),
    expected = expected,
    p = c(rep(NA, 6), mean(nn_null >= nn_obs)),
    rho_with_gc_bias_similarity = cor(v, gc_close, method = "spearman"))
}))
write.csv(concord, file.path(output_dir, "population_concordance.csv"), row.names = FALSE)
print(concord[, c("similarity", "measure", "same_population", "expected", "p", "rho_with_gc_bias_similarity")],
      digits = 3, row.names = FALSE)

# ---- barcode wells ----------------------------------------------------------

v <- sim[ut]
top <- v >= quantile(v, 0.999)
well_rows <- list(
  "same population, same round-1 well" = same_pop & same_well,
  "same population, different round-1 well" = same_pop & !same_well,
  "different population" = !same_pop)
for (i in 1:3) {
  s <- outer(seg[, i], seg[, i], "==")[ut]
  well_rows[[sprintf("different population, share round %d well", 5 - i)]] <- !same_pop & s
  well_rows[[sprintf("different population, different round %d well", 5 - i)]] <- !same_pop & !s
}
wells_tab <- data.frame(pairs = names(well_rows), n = sapply(well_rows, sum),
                        mean_r = sapply(well_rows, function(m) mean(v[m])),
                        frac_top_0.1pct = sapply(well_rows, function(m) mean(top[m])))

stat <- function(p) {
  sp <- outer(p, p, "==")[ut]
  c(mean(v[sp & !same_well]) - mean(v[!sp]), mean(top[sp & !same_well]) / mean(top[!sp]))
}
obs <- stat(pop)
w2p <- tapply(pop, well, `[`, 1)
null <- t(replicate(n_perm, { m <- setNames(sample(w2p), names(w2p)); stat(m[well]) }))
perm <- data.frame(statistic = c("mean r: same pop / different well minus different pop",
                                 "top-0.1% rate ratio: same pop / different well vs different pop"),
                   observed = obs, null_mean = colMeans(null),
                   p = c(mean(null[, 1] >= obs[1]), mean(null[, 2] >= obs[2])))
write.csv(wells_tab, file.path(output_dir, "barcode_well_effects.csv"), row.names = FALSE)
write.csv(perm, file.path(output_dir, "well_preserving_permutation.csv"), row.names = FALSE)
print(wells_tab, digits = 3, row.names = FALSE)
print(perm, digits = 3, row.names = FALSE)

# ---- groups of highly similar cells -----------------------------------------

adj <- sim >= quantile(v, 0.999); diag(adj) <- FALSE
group <- rep(NA_integer_, n); k <- 0
for (i in seq_len(n)) {
  if (!is.na(group[i]) || !any(adj[i, ])) next
  k <- k + 1; stack <- i; group[i] <- k
  while (length(stack)) {
    j <- stack[1]; stack <- stack[-1]
    nb <- which(adj[j, ] & is.na(group))
    group[nb] <- k; stack <- c(stack, nb)
  }
}

shared_regions <- function(idx) {
  m <- rowMeans(resid[, idx, drop = FALSE]); out <- c()
  for (ch in unique(chrom)) {
    s <- which(chrom == ch)
    sm <- as.numeric(stats::filter(m[s], rep(1 / 3, 3), sides = 2)); sm[is.na(sm)] <- m[s][is.na(sm)]
    r <- rle(ifelse(sm > 0.25, 1, ifelse(sm < -0.25, -1, 0)))
    ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    for (q in which(r$values != 0 & r$lengths >= 3)) {
      out <- c(out, sprintf("%s:%.0f-%.0fMb %s", ch, start[s][starts[q]] / 1e6,
                            start[s][ends[q]] / 1e6 + 1, if (r$values[q] > 0) "gain" else "loss"))
    }
  }
  paste(out, collapse = "; ")
}

groups <- do.call(rbind, lapply(sort(unique(na.omit(group))), function(g) {
  idx <- which(group == g)
  data.frame(group = g, n_cells = length(idx),
             mean_r = mean(sim[idx, idx][upper.tri(diag(length(idx)))]),
             pops = paste(names(table(pop[idx])), table(pop[idx]), sep = ":", collapse = " "),
             top_pop_frac = max(table(pop[idx])) / length(idx),
             round1_wells = length(unique(well[idx])),
             median_spikiness = median(qc[cells[idx], "spikiness"]),
             median_gc_bias = median(gc_bias[idx]),
             shared_regions = shared_regions(idx))
}))
groups <- groups[order(-groups$n_cells), ]
write.csv(groups, file.path(output_dir, "similarity_groups.csv"), row.names = FALSE)
write.csv(data.frame(cell = cells, pop = pop, round1_well = well, group = group,
                     spikiness = qc[cells, "spikiness"], gc_bias = gc_bias),
          file.path(output_dir, "cells.csv"), row.names = FALSE)
cat(sprintf("\nthreshold r >= %.3f (top 0.1%% of pairs): %d cells in %d groups; median spikiness %.3f overall\n",
            quantile(v, 0.999), sum(!is.na(group)), k, median(qc[cells, "spikiness"])))
print(groups[groups$n_cells >= 3, setdiff(names(groups), "shared_regions")], digits = 3, row.names = FALSE)
pairs2 <- groups$group[groups$n_cells == 2]
cat(sprintf("size-2 groups: %d, same population %.2f (random %.3f)\n", length(pairs2),
            mean(sapply(pairs2, function(g) length(unique(pop[which(group == g)])) == 1)), expected))
cat("Outputs written to", output_dir, "\n")
