# Sister pairs (cherries) in the per-cell SBS and CNV trees: where does each
# tree's sister pair sit in the other tree?
#
# Topological, not branch-length based: for a pair of cells, the size of the
# smallest clade containing both (2 = also sisters; n = split at the root).
# Null: random pairs of cells in the same trees. Also reports how often sister
# pairs share a population.
#
# CNV tree defaults to the GC-corrected event tree (K562_cnv_event_tree.R with
# K562_CNV_TAG=_gc); the SBS tree is pruned to the same cells.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(dplyr)
})

input_dir <- Sys.getenv(
  "K562_SC_TREES",
  "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/results"
)
cnv_tree_path <- Sys.getenv(
  "K562_CNV_TREE",
  file.path(project_root, "results/K562_tree/sc_trees/cnv_event_tree_gc.newick")
)
output_tag <- Sys.getenv("K562_CONCORDANCE_TAG", "_cnv_events_gc")
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_tree_sister_pairs", output_tag, "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
n_random <- 20000

read_cells <- function(path) drop.tip(read.tree(path), "K562")
t_cnv <- read_cells(cnv_tree_path)
t_sbs <- read_cells(file.path(input_dir, "sc_og_test.newick"))
t_sbs <- drop.tip(t_sbs, setdiff(t_sbs$tip.label, t_cnv$tip.label))
stopifnot(setequal(t_sbs$tip.label, t_cnv$tip.label))
n <- Ntip(t_sbs)

meta <- read.csv(file.path(input_dir, "full_meta.csv"))
pop <- setNames(meta$pop, meta$X)

# MRCA clade size for every pair of cells in a tree, as a cells x cells matrix
mrca_size <- function(tree) {
  m <- mrca(tree)
  sizes <- lengths(Descendants(tree, seq_len(n + tree$Nnode), type = "tips"))
  s <- matrix(sizes[m], nrow(m), dimnames = dimnames(m))
  s
}
size_sbs <- mrca_size(t_sbs)
size_cnv <- mrca_size(t_cnv)
cells <- rownames(size_sbs)
size_cnv <- size_cnv[cells, cells]

cherries <- function(tree) {
  tips <- seq_len(Ntip(tree))
  parent <- tree$edge[match(tips, tree$edge[, 2]), 1]
  pairs <- split(tree$tip.label[tips], parent)
  do.call(rbind, pairs[lengths(pairs) == 2])
}

describe_pairs <- function(pairs, in_tree, other_size) {
  data.frame(
    cell_a = pairs[, 1], cell_b = pairs[, 2], cherry_in = in_tree,
    clade_size_other_tree = other_size[cbind(pairs[, 1], pairs[, 2])],
    same_pop = pop[pairs[, 1]] == pop[pairs[, 2]]
  )
}
pairs <- bind_rows(
  describe_pairs(cherries(t_cnv), "CNV tree", size_sbs),
  describe_pairs(cherries(t_sbs), "SBS tree", size_cnv)
)

# random pairs, sized in each tree
ri <- matrix(sample(n, 2 * n_random, replace = TRUE), ncol = 2)
ri <- ri[ri[, 1] != ri[, 2], ]
random <- bind_rows(
  data.frame(cherry_in = "CNV tree", clade_size_other_tree = size_sbs[ri]),
  data.frame(cherry_in = "SBS tree", clade_size_other_tree = size_cnv[ri])
)
same_pop_random <- mean(pop[cells[ri[, 1]]] == pop[cells[ri[, 2]]])

summ <- pairs %>%
  group_by(cherry_in) %>%
  summarise(
    n_cherries = n(),
    also_cherry_in_other = sum(clade_size_other_tree == 2),
    within_clade_le_10 = sum(clade_size_other_tree <= 10),
    within_clade_le_50 = sum(clade_size_other_tree <= 50),
    median_clade_size_other = median(clade_size_other_tree),
    same_pop_frac = mean(same_pop),
    .groups = "drop"
  ) %>%
  left_join(
    random %>% group_by(cherry_in) %>%
      summarise(random_le_10_frac = mean(clade_size_other_tree <= 10),
                random_le_50_frac = mean(clade_size_other_tree <= 50),
                random_median_clade_size = median(clade_size_other_tree),
                .groups = "drop"),
    by = "cherry_in"
  ) %>%
  mutate(same_pop_random = same_pop_random)

# one-sided rank test: are sister pairs closer in the other tree than random pairs?
summ$wilcoxon_p <- sapply(summ$cherry_in, function(t) {
  wilcox.test(pairs$clade_size_other_tree[pairs$cherry_in == t],
              random$clade_size_other_tree[random$cherry_in == t],
              alternative = "less")$p.value
})
# binomial test: do sister pairs share a population more often than random pairs?
summ$same_pop_p <- sapply(seq_len(nrow(summ)), function(i) {
  binom.test(round(summ$same_pop_frac[i] * summ$n_cherries[i]), summ$n_cherries[i],
             same_pop_random, alternative = "greater")$p.value
})

write.csv(pairs, file.path(output_dir, "sister_pairs.csv"), row.names = FALSE)
write.csv(summ, file.path(output_dir, "sister_pairs_summary.csv"), row.names = FALSE)
cat(sprintf("%d cells in both trees\n", n))
print(as.data.frame(summ), digits = 3)

# ---- plot ------------------------------------------------------------------

# where each tree's sister pairs sit in the other tree, against random pairs
bind_rows(
  mutate(pairs, set = "Sister pairs"),
  mutate(random, set = "Random pairs")
) %>%
  ggplot(aes(clade_size_other_tree, color = set)) +
  stat_ecdf(linewidth = 0.6) +
  scale_x_log10() +
  scale_color_manual(values = c("Sister pairs" = "black", "Random pairs" = "grey60"), name = NULL) +
  facet_wrap(~cherry_in, labeller = labeller(cherry_in = c(
    "CNV tree" = "CNV-tree sisters, in SBS tree",
    "SBS tree" = "SBS-tree sisters, in CNV tree"))) +
  theme_classic() +
  theme(strip.background = element_blank(), legend.position = "bottom") +
  xlab("Smallest clade containing the pair (cells)") +
  ylab("Cumulative fraction of pairs")
ggsave(file.path(output_dir, "sister_pairs_clade_size.svg"),
       bg = "transparent", height = 3.25, width = 5.5)

cat("\nWrote results to:", output_dir, "\n")
