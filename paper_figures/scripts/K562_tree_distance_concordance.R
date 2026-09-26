# Do the per-cell SBS and CNV trees agree on which K562 cells are close?
#
# Inputs (Dustin's sc_PolE_novaseq run, 1000 cells + a K562 root tip):
#   sc_og_test.newick  - NJ on Hamming distance over selected SBS sites
#                        (notebooks/K562_tree.ipynb cells 115-117)
#   cnv_sc_test.newick - NJ on 1 - Pearson r of AneuFinder 1 Mb copy-number
#                        profiles (cells 150-152); override with K562_CNV_TREE,
#                        e.g. the tree from K562_cnv_event_tree.R
#   full_meta.csv      - per-cell population labels
#
# Distances are patristic (sum of branch lengths between two tips). Tests:
#   global - Spearman r over all cell pairs (Mantel), against a full label
#            shuffle and a shuffle restricted to within each population, so
#            shared population structure alone cannot pass the second test
#   local  - for each cell, its k nearest SBS-tree neighbours: the fraction
#            that are also among its k nearest CNV-tree neighbours, and their
#            median CNV-tree neighbour rank, against the same two shuffles
#   per population - both tests within each population's cells only

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

input_dir <- Sys.getenv(
  "K562_SC_TREES",
  "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/results"
)
cnv_tree_path <- Sys.getenv("K562_CNV_TREE", file.path(input_dir, "cnv_sc_test.newick"))
# suffix for the output folder, so runs against different CNV trees coexist
output_tag <- Sys.getenv("K562_CONCORDANCE_TAG", "")
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_tree_distance_concordance", output_tag, "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
n_perm <- 999
ks <- c(5, 10, 20, 50)

pop_colors <- c("P_0" = "#260091", "P_1" = "#1E90FF", "P_2" = "#FFDB58",
                "P_3" = "#FF9D71", "P_4" = "#FF1B5E", "P_5" = "#E3E6E6")

# ---- load ------------------------------------------------------------------

patristic <- function(path) {
  d <- cophenetic.phylo(read.tree(path))
  d[rownames(d) != "K562", colnames(d) != "K562"]
}
d_sbs <- patristic(file.path(input_dir, "sc_og_test.newick"))
d_cnv <- patristic(cnv_tree_path)

# the CNV tree may cover a subset of cells (e.g. dominant ploidy only); patristic
# distances within a subset are unchanged by pruning the other tips
cells <- sort(intersect(rownames(d_sbs), rownames(d_cnv)))
d_sbs <- d_sbs[cells, cells]
d_cnv <- d_cnv[cells, cells]

meta <- read.csv(file.path(input_dir, "full_meta.csv"))
pop <- setNames(meta$pop, meta$X)[cells]
stopifnot(!anyNA(pop))
cat(sprintf("%d cells in both trees\n", length(cells)))

# ---- permutation helpers ---------------------------------------------------

# A permutation relabels cells in the CNV tree only; the SBS tree stays fixed.
shuffle_all <- function(idx, groups) sample(idx)
shuffle_within <- function(idx, groups) {
  out <- idx
  for (g in split(seq_along(idx), groups[idx])) {
    out[g] <- if (length(g) > 1) sample(idx[g]) else idx[g]
  }
  out
}

# Spearman over pairs: rank the CNV distances once; relabelling cells only
# moves those ranks around, so each permutation is a re-index plus cor().
mantel_spearman <- function(idx, shuffle) {
  lt <- lower.tri(d_sbs[idx, idx])
  r_sbs <- rank(d_sbs[idx, idx][lt])
  r_cnv_mat <- matrix(0, length(idx), length(idx))
  r_cnv_mat[lt] <- rank(d_cnv[idx, idx][lt])
  r_cnv_mat <- r_cnv_mat + t(r_cnv_mat)
  pos <- seq_along(idx)
  obs <- cor(r_sbs, r_cnv_mat[lt])
  null <- replicate(n_perm, {
    p <- shuffle(pos, pop[idx])
    cor(r_sbs, r_cnv_mat[p, p][lt])
  })
  c(r = obs, null_mean = mean(null), null_sd = sd(null),
    p = (1 + sum(null >= obs)) / (n_perm + 1))
}

# Row-wise neighbour ranks: row i ranks every other cell by distance to i
# (1 = nearest). Random tie-breaking, since NJ leaves some zero-length edges.
neighbour_ranks <- function(d) {
  diag(d) <- Inf
  t(apply(d, 1, rank, ties.method = "random"))
}

knn_concordance <- function(idx, shuffle) {
  n <- length(idx)
  rk_sbs <- neighbour_ranks(d_sbs[idx, idx])
  rk_cnv <- neighbour_ranks(d_cnv[idx, idx])
  pos <- seq_len(n)
  do.call(rbind, lapply(ks[ks < n], function(k) {
    # (i, j) pairs where j is among i's k nearest cells in the SBS tree
    nn <- which(rk_sbs <= k, arr.ind = TRUE)
    stat <- function(p) {
      r <- rk_cnv[cbind(p[nn[, 1]], p[nn[, 2]])]
      c(overlap = mean(r <= k), median_rank = median(r) / (n - 1))
    }
    obs <- stat(pos)
    null <- replicate(n_perm, stat(shuffle(pos, pop[idx])))
    data.frame(
      k = k,
      overlap = obs[["overlap"]],
      overlap_null_mean = mean(null["overlap", ]),
      overlap_null_lo = unname(quantile(null["overlap", ], 0.025)),
      overlap_null_hi = unname(quantile(null["overlap", ], 0.975)),
      overlap_p = (1 + sum(null["overlap", ] >= obs[["overlap"]])) / (n_perm + 1),
      median_rank_pct = obs[["median_rank"]],
      median_rank_null_mean = mean(null["median_rank", ]),
      median_rank_p = (1 + sum(null["median_rank", ] <= obs[["median_rank"]])) / (n_perm + 1)
    )
  }))
}

# ---- run ---------------------------------------------------------------------

all_idx <- seq_along(cells)
subsets <- c(list(all = all_idx), split(all_idx, pop))

mantel <- list()
knn <- list()
for (s in names(subsets)) {
  idx <- subsets[[s]]
  nulls <- if (s == "all") list(shuffle_all = shuffle_all, shuffle_within_pop = shuffle_within)
           else list(shuffle_all = shuffle_all)
  for (nm in names(nulls)) {
    cat(sprintf("%s (%d cells), %s\n", s, length(idx), nm))
    mantel[[length(mantel) + 1]] <- data.frame(
      subset = s, n_cells = length(idx), null = nm,
      t(mantel_spearman(idx, nulls[[nm]]))
    )
    knn[[length(knn) + 1]] <- data.frame(
      subset = s, n_cells = length(idx), null = nm,
      knn_concordance(idx, nulls[[nm]])
    )
  }
}
mantel <- bind_rows(mantel)
knn <- bind_rows(knn)
write.csv(mantel, file.path(output_dir, "mantel_spearman.csv"), row.names = FALSE)
write.csv(knn, file.path(output_dir, "knn_concordance.csv"), row.names = FALSE)

print(mantel, row.names = FALSE, digits = 3)
print(knn %>% select(subset, null, k, overlap, overlap_null_mean, overlap_p,
                     median_rank_pct, median_rank_p),
      row.names = FALSE, digits = 3)

# ---- plots -----------------------------------------------------------------

# 1. every cell pair, SBS vs CNV patristic distance
lt <- lower.tri(d_sbs)
pair_pop <- outer(pop, pop, function(a, b) ifelse(a == b, "same population", "different populations"))
pairs <- data.frame(sbs = d_sbs[lt], cnv = d_cnv[lt], type = pair_pop[lt])
r_all <- mantel$r[mantel$subset == "all"][1]

ggplot(pairs, aes(sbs, cnv)) +
  geom_bin2d(bins = 80) +
  scale_fill_gradient(low = "#dbe4f3", high = "#1b3a6b", trans = "log10", name = "Pairs") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
           label = sprintf("Spearman r = %.3f", r_all)) +
  theme_classic() +
  xlab("SBS tree distance") +
  ylab("CNV tree distance")
ggsave(file.path(output_dir, "pairwise_distance_hexbin.svg"),
       bg = "transparent", height = 3, width = 3.75)

# 2. local neighbour overlap vs k: observed against each null's 95% band
knn_all <- knn %>% filter(subset == "all") %>%
  mutate(null = recode(null, shuffle_all = "Shuffle all cells",
                       shuffle_within_pop = "Shuffle within population"))
ggplot(knn_all, aes(x = factor(k))) +
  geom_linerange(aes(ymin = overlap_null_lo, ymax = overlap_null_hi),
                 linewidth = 3, color = "grey80") +
  geom_point(aes(y = overlap_null_mean), shape = 95, size = 5, color = "grey40") +
  geom_point(aes(y = overlap), size = 2.5) +
  facet_wrap(~null) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  xlab("k nearest SBS-tree neighbours") +
  ylab("Also k-nearest in CNV tree")
ggsave(file.path(output_dir, "knn_overlap_vs_k.svg"),
       bg = "transparent", height = 2.75, width = 4.5)

# 3. per-population Mantel r with its shuffle null
per_pop <- mantel %>% filter(subset != "all")
ggplot(per_pop, aes(x = subset)) +
  geom_linerange(aes(ymin = null_mean - 1.96 * null_sd, ymax = null_mean + 1.96 * null_sd),
                 linewidth = 3, color = "grey80") +
  geom_point(aes(y = r, fill = subset), shape = 21, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = pop_colors) +
  scale_x_discrete(labels = c("P0", "P1", "P2", "P3", "P4", "P5")) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Population") +
  ylab("Spearman r, SBS vs CNV")
ggsave(file.path(output_dir, "per_population_mantel.svg"),
       bg = "transparent", height = 2.75, width = 3)

cat("\nWrote results to:", output_dir, "\n")
