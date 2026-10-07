# Do K562 cells that sit next to each other on the copy-number dendrogram come
# from the same population (P_0-P_5) more often than chance?
#
# Input: an AneuFinder run (K562_AF_DIR; default the repeat-masked, smoothed
# 250 kb run from scripts/utils/run_aneufinder_smooth_masked.sh). Cells: all
# except the off-peak artefacts (K562_offpeak_ploidy.R). Tree: AneuFinder's
# clusterHMMs(), the clustering drawn on its genome-wide heatmap.
#
# Statistics, on the dendrogram's leaf order:
#   adjacent_k - fraction of cells k leaves apart (k = 1-5) in the same population
#   cherries   - fraction of sister pairs (two leaves joined directly) in the
#                same population
# Nulls (n_perm each): population labels shuffled across cells; and whole
# round-1 barcode wells (barcode segment 4, which encodes population) shuffled
# between populations, so any well effect is kept in the null.
# Controls: the same statistics for per-cell GC-bias quintiles
# (K562_gc_toggle_ploidy.R; largely a noise measure, rho -0.77 with spikiness)
# and for spikiness quintiles (AneuFinder QC, 1 Mb run), to show how much the
# tree groups cells by GC bias or noise; cell-cycle stage (G1/G2, early, mid,
# late S; K562_replication_timing_coverage.R) to show whether it groups cells by
# replication state; and the population statistics with labels shuffled within
# spikiness quintiles.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(AneuFinder)
})

af_dir <- Sys.getenv("K562_AF_DIR", file.path(project_root, "results/K562_tree/aneufinder_gc_smooth500k/output"))
run_tag <- basename(dirname(af_dir))
output_dir <- file.path(project_root, "paper_figures/output/K562_cnv_tree_neighbours", run_tag)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
n_perm <- 10000
set.seed(1)

model_files <- list.files(file.path(af_dir, "MODELS/method-edivisive"), "\\.RData$", full.names = TRUE)
cell_of <- sub("(\\.bam)?_binsize.*", "", basename(model_files))
metrics <- read.csv(file.path(project_root, "paper_figures/output/K562_offpeak_ploidy/per_cell_metrics.csv"))
keep <- !(cell_of %in% metrics$cell[grepl("^artefact", metrics$class)])
model_files <- model_files[keep]; cell_of <- cell_of[keep]

hmms <- loadFromFiles(model_files, check.class = c("aneuHMM", "aneuBiHMM"))
cl <- clusterHMMs(hmms, cluster = TRUE)
ids <- sub("\\.bam$", "", sapply(hmms, `[[`, "ID"))
ord <- ids[cl$IDorder]
n <- length(ord)

meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
pop <- setNames(meta$pop, meta$cell)[ord]
# barcode = seg1 AGGA seg2 ACTC seg3 AAGG seg4; seg4 is the round-1 well
well <- sapply(regmatches(ord, regexec("^(\\w{8})AGGA(\\w{8})ACTC(\\w{8})AAGG(\\w+)$", ord)), `[`, 5)
gcb <- read.csv(file.path(project_root, "paper_figures/output/K562_gc_toggle_ploidy/ploidy_with_without_gc.csv"))
gc_q <- cut(setNames(gcb$gc_bias, gcb$cell)[ord], quantile(gcb$gc_bias, 0:5 / 5), include.lowest = TRUE, labels = FALSE)
spik <- setNames(metrics$spikiness, metrics$cell)[ord]
spik_q <- cut(spik, quantile(spik, 0:5 / 5), include.lowest = TRUE, labels = FALSE)
ph <- read.csv(file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/per_cell_phase.csv"))
stage <- ifelse(ph$phase == "S", as.character(ph$s_stage), "G1/G2")
stage <- setNames(stage, ph$cell)[ord]
stopifnot(!anyNA(pop), !anyNA(well), !anyNA(gc_q), !anyNA(spik_q), !anyNA(stage))
cat(sprintf("%s: %d cells on the tree\n", run_tag, n))

# sister pairs: hclust merges of two leaves, as positions in the leaf order
hc <- cl$hclust
pos <- match(seq_len(n), hc$order)
cherry <- hc$merge[hc$merge[, 1] < 0 & hc$merge[, 2] < 0, , drop = FALSE]
cherry_pos <- cbind(pos[-cherry[, 1]], pos[-cherry[, 2]])

stats_of <- function(lab) {
  c(setNames(sapply(1:5, function(k) mean(lab[-(1:k)] == lab[1:(n - k)])), paste0("adjacent_", 1:5)),
    cherries = mean(lab[cherry_pos[, 1]] == lab[cherry_pos[, 2]]))
}

test <- function(lab, label_name) {
  obs <- stats_of(lab)
  null <- t(replicate(n_perm, stats_of(sample(lab))))
  data.frame(labels = label_name, statistic = names(obs), observed = obs,
             null_mean = colMeans(null), enrichment = obs / colMeans(null),
             p = (colSums(sweep(null, 2, obs, ">=")) + 1) / (n_perm + 1))
}

res <- rbind(test(pop, "population"), test(gc_q, "GC-bias quintile (control)"),
             test(spik_q, "spikiness quintile (control)"),
             test(stage, "cell-cycle stage (control)"))

# population, labels shuffled within spikiness quintiles
obs <- stats_of(pop)
snull <- t(replicate(n_perm, { p <- pop; for (s in 1:5) { i <- which(spik_q == s); p[i] <- p[i][sample.int(length(i))] }; stats_of(p) }))
res <- rbind(res, data.frame(labels = "population (shuffled within spikiness quintiles)", statistic = names(obs),
                             observed = obs, null_mean = colMeans(snull), enrichment = obs / colMeans(snull),
                             p = (colSums(sweep(snull, 2, obs, ">=")) + 1) / (n_perm + 1)))

# well-preserving null for population
w2p <- tapply(pop, well, `[`, 1)
obs <- stats_of(pop)
wnull <- t(replicate(n_perm, { m <- setNames(sample(w2p), names(w2p)); stats_of(m[well]) }))
res <- rbind(res, data.frame(labels = "population (round-1 wells shuffled)", statistic = names(obs),
                             observed = obs, null_mean = colMeans(wnull), enrichment = obs / colMeans(wnull),
                             p = (colSums(sweep(wnull, 2, obs, ">=")) + 1) / (n_perm + 1)))

# adjacent same-population pairs split by whether the two cells share a
# round-1 well (populations were each indexed in their own wells, so a
# population-wide difference predicts enrichment in both; a per-well effect
# only in the first)
split_counts <- function(p, w) {
  sp <- p[-1] == p[-n]; sw <- w[-1] == w[-n]
  c(same_well = sum(sp & sw), different_well = sum(sp & !sw))
}
obs_split <- split_counts(pop, well)
null_split <- t(replicate(n_perm, { i <- sample(n); split_counts(pop[i], well[i]) }))
well_split <- data.frame(pairs = paste("adjacent, same population,", c("same well", "different well")),
                         observed = obs_split, expected = colMeans(null_split),
                         enrichment = obs_split / colMeans(null_split),
                         p = (colSums(sweep(null_split, 2, obs_split, ">=")) + 1) / (n_perm + 1))
write.csv(well_split, file.path(output_dir, "adjacent_same_vs_different_well.csv"), row.names = FALSE)
print(well_split, digits = 3, row.names = FALSE)

# per population: adjacent same-population pairs, observed vs expected
per_pop <- do.call(rbind, lapply(sort(unique(pop)), function(q) {
  is_q <- pop == q
  obs_q <- sum(is_q[-1] & is_q[-n])
  null_q <- replicate(n_perm, { s <- sample(is_q); sum(s[-1] & s[-n]) })
  data.frame(pop = q, cells = sum(is_q), adjacent_pairs = obs_q, expected = mean(null_q),
             enrichment = obs_q / mean(null_q), p = (sum(null_q >= obs_q) + 1) / (n_perm + 1))
}))

write.csv(res, file.path(output_dir, "neighbour_concordance.csv"), row.names = FALSE)
write.csv(per_pop, file.path(output_dir, "per_population_adjacency.csv"), row.names = FALSE)
write.csv(data.frame(leaf = seq_len(n), cell = ord, pop = pop, round1_well = well, gc_bias_quintile = gc_q),
          file.path(output_dir, "leaf_order.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
print(per_pop, digits = 3, row.names = FALSE)
cat("Outputs written to", output_dir, "\n")
