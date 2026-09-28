# Are the per-cell SBS and CNV trees topologically similar?
#
# Same inputs as K562_tree_distance_concordance.R: sc_og_test.newick and the
# GC-corrected CNV event tree (override with K562_CNV_TREE); the K562 root tip
# is dropped. Each statistic is compared with 999 random relabellings of the
# CNV tree's tips, which keep both tree shapes and destroy any correspondence:
#   Robinson-Foulds  - normalised RF distance over all splits (unrooted)
#   shared cherries  - sister-tip pairs present in both trees
#   clade matching   - for each SBS-tree clade, the best Jaccard overlap with
#                      any CNV-tree clade, summarised by clade size, so local
#                      (small clade) and deep agreement are reported separately

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

input_dir <- Sys.getenv(
  "K562_SC_TREES",
  "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/results"
)
cnv_tree_path <- Sys.getenv("K562_CNV_TREE",
                           file.path(project_root, "results/K562_tree/sc_trees/cnv_event_tree_gc.newick"))
# suffix for the output folder, so runs against different CNV trees coexist
output_tag <- Sys.getenv("K562_CONCORDANCE_TAG", "_cnv_events_gc")
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_tree_topology_concordance", output_tag, "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
n_perm <- 999
size_bins <- c(2, 5, 10, 20, 50, 100, 1000)

read_cells <- function(path) drop.tip(read.tree(path), "K562")
t_sbs <- read_cells(file.path(input_dir, "sc_og_test.newick"))
t_cnv <- read_cells(cnv_tree_path)
# the CNV tree may cover a subset of cells (e.g. dominant ploidy only); compare
# against the SBS tree pruned to the same cells
t_sbs <- drop.tip(t_sbs, setdiff(t_sbs$tip.label, t_cnv$tip.label))
stopifnot(setequal(t_sbs$tip.label, t_cnv$tip.label))
cells <- sort(t_sbs$tip.label)
n <- length(cells)
cat(sprintf("%d cells in both trees\n", n))

relabel <- function(tree, perm) {
  tree$tip.label <- cells[perm][match(tree$tip.label, cells)]
  tree
}

# ---- clade membership ------------------------------------------------------

# Clades of a rooted tree as a sparse clades x cells 0/1 matrix, dropping the
# root (all cells) and single tips.
clade_matrix <- function(tree) {
  parts <- prop.part(tree)
  labs <- attr(parts, "labels")
  keep <- lengths(parts) >= 2 & lengths(parts) < n
  parts <- parts[keep]
  sparseMatrix(
    i = rep(seq_along(parts), lengths(parts)),
    j = match(labs[unlist(parts)], cells),
    x = 1, dims = c(length(parts), n)
  )
}
m_sbs <- clade_matrix(t_sbs)
m_cnv <- clade_matrix(t_cnv)
size_sbs <- rowSums(m_sbs)
size_cnv <- rowSums(m_cnv)
size_bin <- cut(size_sbs, size_bins, right = FALSE,
                labels = paste0(head(size_bins, -1), "-", tail(size_bins, -1) - 1))

# Best Jaccard of each SBS clade against every CNV clade, with the CNV tree's
# cells relabelled by perm (a column permutation of its clade matrix).
best_jaccard <- function(perm) {
  m <- m_cnv[, order(perm)]
  inter <- as.matrix(tcrossprod(m_sbs, m))
  union <- outer(size_sbs, size_cnv, "+") - inter
  apply(inter / union, 1, max)
}

cherries <- function(m, sizes) {
  pairs <- m[sizes == 2, , drop = FALSE]
  apply(pairs, 1, function(r) paste(which(r > 0), collapse = "-"))
}
cherry_sbs <- cherries(m_sbs, size_sbs)

# ---- observed and null -----------------------------------------------------

observed <- function(perm) {
  t_cnv_p <- relabel(t_cnv, perm)
  jac <- best_jaccard(perm)
  cherry_cnv <- cherries(m_cnv[, order(perm)], size_cnv)
  c(
    rf_norm = RF.dist(unroot(t_sbs), unroot(t_cnv_p), normalize = TRUE),
    shared_cherries = sum(cherry_sbs %in% cherry_cnv),
    tapply(jac, size_bin, mean)
  )
}

obs <- observed(seq_len(n))
null <- replicate(n_perm, observed(sample(n)))

summary_tbl <- data.frame(
  statistic = names(obs),
  observed = obs,
  null_mean = rowMeans(null),
  null_lo = apply(null, 1, quantile, 0.025),
  null_hi = apply(null, 1, quantile, 0.975),
  # RF is a distance (similar = low); everything else is a similarity
  p = ifelse(names(obs) == "rf_norm",
             (1 + rowSums(null <= obs)) / (n_perm + 1),
             (1 + rowSums(null >= obs)) / (n_perm + 1)),
  row.names = NULL
)
summary_tbl$n_sbs_clades <- c(NA, sum(size_sbs == 2), as.vector(table(size_bin)))
write.csv(summary_tbl, file.path(output_dir, "topology_concordance.csv"), row.names = FALSE)
print(summary_tbl, digits = 3)
cat(sprintf("\nCherries: %d in the SBS tree, %d in the CNV tree\n",
            sum(size_sbs == 2), sum(size_cnv == 2)))

# ---- plot ------------------------------------------------------------------

clade_rows <- summary_tbl %>%
  filter(!statistic %in% c("rf_norm", "shared_cherries")) %>%
  mutate(statistic = factor(statistic, levels = levels(size_bin)))

# Plotted as the excess over the null mean: the null bands are too narrow to
# see on the raw Jaccard scale, which is set mostly by clade size.
ggplot(clade_rows, aes(x = statistic)) +
  geom_linerange(aes(ymin = null_lo - null_mean, ymax = null_hi - null_mean),
                 linewidth = 3, color = "grey80") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(y = observed - null_mean), size = 2.5) +
  theme_classic() +
  xlab("SBS-tree clade size (cells)") +
  ylab("Best-match Jaccard\nminus shuffled mean")
ggsave(file.path(output_dir, "clade_match_by_size.svg"),
       bg = "transparent", height = 2.75, width = 4.25)

cat("\nWrote results to:", output_dir, "\n")
