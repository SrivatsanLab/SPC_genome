# Do small, clearly called CNV clades segregate by population, and to P5 in
# particular?
#
# Clades come from K562_cnv_clade_sbs_test.R (gain/loss events >= 5 Mb carried
# by 3-50 high-quality cells, merged at Jaccard >= 0.8). Per clade:
#   P5 enrichment  - one-sided hypergeometric test of the number of P5 members,
#                    given P5's share of the high-quality cells
#   any population - permutation test of the chi-square statistic of member
#                    population counts against random same-size groups
# Across clades: is the fraction of members from the clade's own most common
# population higher than for random groups of the same sizes?

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

sc_dir <- Sys.getenv("K562_SC_TREES",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/results")
tag <- Sys.getenv("K562_CNV_TAG", "_gc")
clade_dir <- file.path(project_root, paste0("paper_figures/output/K562_cnv_clade_sbs_test", tag, "/"))
ref_csv <- file.path(project_root, "paper_figures/output/K562_cnv_vs_reference/per_cell_agreement.csv")
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_cnv_clade_population", tag, "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

focus_pop <- "P_5"
min_ref_r <- 0.7
n_perm <- 9999
set.seed(1)

tbl <- read.csv(file.path(clade_dir, "cnv_clades_in_sbs_tree.csv"))
members <- strsplit(readLines(file.path(clade_dir, "cnv_clade_members.txt")), " ")
members <- setNames(lapply(members, function(x) strsplit(x[2], ",")[[1]]),
                    sapply(members, `[`, 1))

# the same high-quality cell pool the clades were drawn from
em <- readRDS(file.path(project_root, "results/K562_tree/sc_trees", paste0("cnv_event_matrix", tag, ".rds")))
ref <- read.csv(ref_csv) %>% filter(run == "new", bins == "All bins")
ref_r <- setNames(ref$pearson, ref$cell)
pool <- colnames(em$features)[ref_r[colnames(em$features)] >= min_ref_r]
meta <- read.csv(file.path(sc_dir, "full_meta.csv"))
pop <- setNames(meta$pop, meta$X)[pool]
pops <- sort(unique(pop))
pool_frac <- as.vector(table(factor(pop, levels = pops))) / length(pool)
n_focus <- sum(pop == focus_pop)
cat(sprintf("%d cells in the pool; %s = %d (%.1f%%)\n", length(pool), focus_pop, n_focus,
            100 * n_focus / length(pool)))

chisq_stat <- function(counts) {
  e <- sum(counts) * pool_frac
  sum((counts - e)^2 / e)
}
res <- bind_rows(lapply(names(members), function(k) {
  cells <- members[[k]]
  sz <- length(cells)
  counts <- as.vector(table(factor(pop[cells], levels = pops)))
  k_focus <- counts[pops == focus_pop]
  obs <- chisq_stat(counts)
  null <- replicate(n_perm, chisq_stat(as.vector(table(factor(sample(pop, sz), levels = pops)))))
  data.frame(
    clade = as.integer(k), n_cells = sz,
    focus_members = k_focus,
    focus_expected = sz * n_focus / length(pool),
    p_focus = phyper(k_focus - 1, n_focus, length(pool) - n_focus, sz, lower.tail = FALSE),
    top_pop = pops[which.max(counts)], top_pop_frac = max(counts) / sz,
    p_any_pop = (1 + sum(null >= obs)) / (n_perm + 1),
    pops = paste(counts, collapse = "/")
  )
})) %>%
  left_join(select(tbl, clade, events, event_mb, p_sbs = p), by = "clade") %>%
  mutate(q_focus = p.adjust(p_focus, "BH"), q_any_pop = p.adjust(p_any_pop, "BH")) %>%
  arrange(p_focus)

# across clades: purity (top-population fraction) vs random groups of the same sizes
purity <- function(cells_pop) max(table(cells_pop)) / length(cells_pop)
obs_purity <- mean(res$top_pop_frac)
null_purity <- replicate(999, mean(sapply(res$n_cells, function(sz) purity(sample(pop, sz)))))
obs_focus <- sum(res$focus_members)
null_focus <- replicate(9999, sum(sapply(res$n_cells, function(sz) sum(sample(pop, sz) == focus_pop))))

write.csv(res, file.path(output_dir, "cnv_clade_population.csv"), row.names = FALSE)

cat(sprintf("\n%d clades. %s-enriched at p < 0.05: %d (%.1f expected); BH q < 0.1: %d\n",
            nrow(res), focus_pop, sum(res$p_focus < 0.05), 0.05 * nrow(res), sum(res$q_focus < 0.1)))
cat(sprintf("any-population association at p < 0.05: %d (%.1f expected); BH q < 0.1: %d\n",
            sum(res$p_any_pop < 0.05), 0.05 * nrow(res), sum(res$q_any_pop < 0.1)))
cat(sprintf("mean top-population fraction: observed %.3f, random %.3f (p = %.3f)\n",
            obs_purity, mean(null_purity), (1 + sum(null_purity >= obs_purity)) / 1000))
cat(sprintf("total %s memberships across clades: observed %d, random %.1f (p = %.4f)\n",
            focus_pop, obs_focus, mean(null_focus), (1 + sum(null_focus >= obs_focus)) / 10000))
cat("\npops column: P0/P1/P2/P3/P4/P5 member counts\n")
print(head(res %>% select(events, n_cells, focus_members, focus_expected, p_focus, q_focus,
                          p_any_pop, pops, p_sbs), 15), row.names = FALSE, digits = 3, right = FALSE)

# ---- plot ------------------------------------------------------------------

ggplot(res, aes(n_cells, focus_members / n_cells)) +
  geom_hline(yintercept = n_focus / length(pool), linetype = "dashed", color = "grey50") +
  geom_point(aes(fill = q_focus < 0.1), shape = 21, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"), name = "BH q < 0.1") +
  scale_x_log10() +
  theme_classic() +
  xlab("CNV clade size (cells)") +
  ylab(sprintf("Fraction of members from %s", sub("_", "", focus_pop)))
ggsave(file.path(output_dir, "cnv_clade_P5_fraction.svg"), bg = "transparent", height = 3, width = 4)

cat("\nWrote results to:", output_dir, "\n")
