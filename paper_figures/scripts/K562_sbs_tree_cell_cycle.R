# Does the per-cell K562 SBS tree (sc_og_test.newick; NJ on Hamming distance
# over selected SBS sites, notebooks/K562_tree.ipynb) follow cell-cycle /
# replication state inferred from coverage?
#
# Replication state: K562_replication_timing_coverage.R (phase, S stage) and
# K562_cell_cycle_profile.R (inferred DNA content = 1 + fraction replicated).
# Tests, 10,000 label shuffles unless noted:
#   sister pairs  - fraction of SBS-tree cherries sharing stage (G1/G2, early,
#                   mid, late S) and sharing phase (S vs G1/G2)
#   distance      - Spearman of patristic distance with |difference in inferred
#                   DNA content| over cell pairs, vs 999 label shuffles
#   depth         - root-to-tip distance by phase (Wilcoxon) and its Spearman
#                   with DNA content and with spikiness
#   SBS counts    - per-cell SNP_count_de_novo and tot_observed_sites
#                   (full_meta.csv) by phase
# Population is reported alongside for comparison.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(ape))

sc_dir <- Sys.getenv("K562_SC_TREES", "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/results")
output_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_cell_cycle/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
n_perm <- 10000
set.seed(1)

tree_k <- read.tree(file.path(sc_dir, "sc_og_test.newick"))
root_tip <- "K562"
depth_all <- node.depth.edgelength(tree_k)[seq_len(Ntip(tree_k))]
names(depth_all) <- tree_k$tip.label
tree <- drop.tip(tree_k, root_tip)
tips <- tree$tip.label; n <- length(tips)

ph <- read.csv(file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/per_cell_phase.csv"))
dna <- read.csv(file.path(project_root, "paper_figures/output/K562_cell_cycle_profile/inferred_dna_content.csv"))
meta <- read.csv(file.path(sc_dir, "full_meta.csv"))
stage <- setNames(ifelse(ph$phase == "S", as.character(ph$s_stage), "G1/G2"), ph$cell)[tips]
phase <- setNames(ph$phase, ph$cell)[tips]
spik <- setNames(ph$spikiness, ph$cell)[tips]
dnac <- setNames(dna$dna_content, dna$cell)[tips]
pop <- setNames(meta$pop, meta$X)[tips]
stopifnot(!anyNA(stage), !anyNA(dnac), !anyNA(pop))

# ---- sister pairs ----------------------------------------------------------------
kids <- split(tree$edge[, 2], tree$edge[, 1])
cherry <- do.call(rbind, Filter(Negate(is.null), lapply(kids, function(k) if (length(k) == 2 && all(k <= n)) k else NULL)))
same <- function(lab) mean(lab[cherry[, 1]] == lab[cherry[, 2]])
sis <- do.call(rbind, lapply(list("cell-cycle stage" = stage, "S vs G1/G2" = phase, "population" = pop), function(lab) {
  o <- same(lab); nl <- replicate(n_perm, same(sample(lab)))
  data.frame(observed = o, expected = mean(nl), enrichment = o / mean(nl), p = (sum(nl >= o) + 1) / (n_perm + 1))
}))
sis <- cbind(label = rownames(sis), cherries = nrow(cherry), sis)

# ---- distance vs DNA content difference ---------------------------------------------
pd <- cophenetic.phylo(tree)[tips, tips]
ut <- upper.tri(pd)
dd <- abs(outer(dnac, dnac, "-"))[ut]
rho <- cor(pd[ut], dd, method = "spearman")
rho_null <- replicate(999, { s <- sample(dnac); cor(pd[ut], abs(outer(s, s, "-"))[ut], method = "spearman") })
dist_res <- data.frame(statistic = "Spearman(patristic distance, |delta DNA content|)", observed = rho,
                       null_mean = mean(rho_null), null_sd = sd(rho_null), p = (sum(rho_null >= rho) + 1) / 1000)

# ---- depth and SBS counts ---------------------------------------------------------
depth <- depth_all[tips]
cnt <- meta[match(tips, meta$X), c("SNP_count_de_novo", "tot_observed_sites")]
depth_res <- data.frame(
  statistic = c("root-to-tip: median G1/G2", "root-to-tip: median S", "root-to-tip: Wilcoxon p (S vs G1/G2)",
                "Spearman(root-to-tip, DNA content)", "Spearman(root-to-tip, spikiness)",
                "de novo SBS: median G1/G2", "de novo SBS: median S", "de novo SBS: Wilcoxon p",
                "observed sites: median G1/G2", "observed sites: median S", "observed sites: Wilcoxon p",
                "Spearman(de novo SBS, DNA content)"),
  value = c(median(depth[phase != "S"]), median(depth[phase == "S"]), wilcox.test(depth ~ phase)$p.value,
            cor(depth, dnac, method = "spearman"), cor(depth, spik, method = "spearman", use = "complete"),
            median(cnt$SNP_count_de_novo[phase != "S"]), median(cnt$SNP_count_de_novo[phase == "S"]),
            wilcox.test(cnt$SNP_count_de_novo ~ phase)$p.value,
            median(cnt$tot_observed_sites[phase != "S"]), median(cnt$tot_observed_sites[phase == "S"]),
            wilcox.test(cnt$tot_observed_sites ~ phase)$p.value,
            cor(cnt$SNP_count_de_novo, dnac, method = "spearman")))

write.csv(sis, file.path(output_dir, "sbs_tree_sisters_by_cell_cycle.csv"), row.names = FALSE)
write.csv(dist_res, file.path(output_dir, "sbs_tree_distance_vs_dna_content.csv"), row.names = FALSE)
write.csv(depth_res, file.path(output_dir, "sbs_tree_depth_and_counts_by_phase.csv"), row.names = FALSE)
cat(sprintf("%d cells, %d SBS-tree sister pairs\n", n, nrow(cherry)))
print(sis, digits = 3, row.names = FALSE); print(dist_res, digits = 3, row.names = FALSE)
print(depth_res, digits = 4, row.names = FALSE)
