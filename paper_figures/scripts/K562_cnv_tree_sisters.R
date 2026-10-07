# Do K562 cells from the same population (P_0-P_5) end up as sisters in the
# copy-number distance tree more often than chance?
#
# Tree: K562_cnv_distance_tree.R (neighbour joining on mean |CN difference|,
# rooted on the K562 consensus; K562_TREE_RUN picks the run, default the
# replication-corrected calls; e.g. aneufinder_fragment_gc_rtcorr_G1G2 for the
# tree of G1/G2 cells only, where the cell-cycle null and control are skipped). Sister pairs = cherries: two cells joined
# directly by one internal node. Statistic: fraction of sister pairs from the
# same population. Nulls, 10,000 shuffles each:
#   labels       - population labels shuffled across cells
#   within stage - shuffled within cell-cycle stage (G1/G2, early, mid, late S;
#                  K562_replication_timing_coverage.R), so pairing of cells in
#                  the same stage cannot create population concordance
#   within noise - shuffled within spikiness quintiles
#   whole wells  - whole round-1 barcode wells (barcode segment 4, which marks
#                  population) shuffled between populations
# Also: same-population sisters split by shared vs different round-1 well;
# the same statistic for cell-cycle stage and spikiness quintile as controls;
# per-population counts.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(ape))

run_tag <- Sys.getenv("K562_TREE_RUN", "aneufinder_fragment_gc_rtcorr")
tree_dir <- file.path(project_root, "paper_figures/output/K562_cnv_distance_tree", run_tag)
output_dir <- file.path(project_root, "paper_figures/output/K562_cnv_tree_sisters", run_tag)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
n_perm <- 10000
set.seed(1)

tree <- drop.tip(read.tree(file.path(tree_dir, "cnv_distance_tree.newick")), "K562_consensus")
n <- Ntip(tree)
tips <- tree$tip.label

# cherries: internal nodes whose two children are both tips
kids <- split(tree$edge[, 2], tree$edge[, 1])
cherry <- do.call(rbind, Filter(Negate(is.null), lapply(kids, function(k)
  if (length(k) == 2 && all(k <= n)) k else NULL)))
cat(sprintf("%s: %d cells, %d sister pairs (%d cells in a pair)\n", run_tag, n, nrow(cherry), 2 * nrow(cherry)))

meta <- read.csv(file.path(project_root, "paper_figures/data/Anneufinder/cell_ploidy.csv"))
pop <- setNames(meta$pop, meta$cell)[tips]
well <- sapply(regmatches(tips, regexec("^(\\w{8})AGGA(\\w{8})ACTC(\\w{8})AAGG(\\w+)$", tips)), `[`, 5)
ph <- read.csv(file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/per_cell_phase.csv"))
stage <- setNames(ifelse(ph$phase == "S", as.character(ph$s_stage), "G1/G2"), ph$cell)[tips]
metrics <- read.csv(file.path(project_root, "paper_figures/output/K562_offpeak_ploidy/per_cell_metrics.csv"))
spik <- setNames(metrics$spikiness, metrics$cell)[tips]
spik_q <- cut(spik, quantile(spik, 0:5 / 5), include.lowest = TRUE, labels = FALSE)
stopifnot(!anyNA(pop), !anyNA(well), !anyNA(stage), !anyNA(spik_q))

same <- function(lab) mean(lab[cherry[, 1]] == lab[cherry[, 2]])
shuffle_within <- function(lab, strata) {
  for (s in unique(strata)) { i <- which(strata == s); lab[i] <- lab[i][sample.int(length(i))] }
  lab
}
w2p <- tapply(pop, well, `[`, 1)
obs <- same(pop)
nulls <- list(
  "labels shuffled" = replicate(n_perm, same(sample(pop))))
if (length(unique(stage)) > 1) nulls[["shuffled within cell-cycle stage"]] <- replicate(n_perm, same(shuffle_within(pop, stage)))
nulls <- c(nulls, list(
  "shuffled within spikiness quintile" = replicate(n_perm, same(shuffle_within(pop, spik_q))),
  "whole round-1 wells shuffled" = replicate(n_perm, { m <- setNames(sample(w2p), names(w2p)); same(m[well]) })))
res <- do.call(rbind, lapply(names(nulls), function(k) data.frame(
  null = k, same_population = obs, expected = mean(nulls[[k]]), enrichment = obs / mean(nulls[[k]]),
  p = (sum(nulls[[k]] >= obs) + 1) / (n_perm + 1))))

ctrl_labels <- list("cell-cycle stage" = stage, "spikiness quintile" = spik_q)
ctrl_labels <- ctrl_labels[sapply(ctrl_labels, function(l) length(unique(l)) > 1)]
ctrl <- do.call(rbind, lapply(ctrl_labels, function(lab) {
  o <- same(lab); nl <- replicate(n_perm, same(sample(lab)))
  data.frame(same = o, expected = mean(nl), enrichment = o / mean(nl), p = (sum(nl >= o) + 1) / (n_perm + 1))
}))
ctrl <- cbind(label = rownames(ctrl), ctrl)

sp <- pop[cherry[, 1]] == pop[cherry[, 2]]; sw <- well[cherry[, 1]] == well[cherry[, 2]]
split_null <- t(replicate(n_perm, { i <- sample(n); p <- pop[i]; w <- well[i]
  s1 <- p[cherry[, 1]] == p[cherry[, 2]]; s2 <- w[cherry[, 1]] == w[cherry[, 2]]
  c(sum(s1 & s2), sum(s1 & !s2)) }))
well_split <- data.frame(sisters = c("same population, same well", "same population, different well"),
                         observed = c(sum(sp & sw), sum(sp & !sw)), expected = colMeans(split_null))
well_split$enrichment <- well_split$observed / well_split$expected
well_split$p <- c((sum(split_null[, 1] >= well_split$observed[1]) + 1) / (n_perm + 1),
                  (sum(split_null[, 2] >= well_split$observed[2]) + 1) / (n_perm + 1))

per_pop <- do.call(rbind, lapply(sort(unique(pop)), function(q) {
  both <- sum(pop[cherry[, 1]] == q & pop[cherry[, 2]] == q)
  nl <- replicate(n_perm, { p <- sample(pop); sum(p[cherry[, 1]] == q & p[cherry[, 2]] == q) })
  data.frame(pop = q, cells = sum(pop == q), sister_pairs_within = both, expected = mean(nl),
             enrichment = both / mean(nl), p = (sum(nl >= both) + 1) / (n_perm + 1))
}))

write.csv(res, file.path(output_dir, "sister_population_concordance.csv"), row.names = FALSE)
write.csv(ctrl, file.path(output_dir, "sister_controls.csv"), row.names = FALSE)
write.csv(well_split, file.path(output_dir, "sister_same_vs_different_well.csv"), row.names = FALSE)
write.csv(per_pop, file.path(output_dir, "sister_per_population.csv"), row.names = FALSE)
write.csv(data.frame(cell1 = tips[cherry[, 1]], cell2 = tips[cherry[, 2]],
                     pop1 = pop[cherry[, 1]], pop2 = pop[cherry[, 2]]),
          file.path(output_dir, "sister_pairs.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE); print(well_split, digits = 3, row.names = FALSE)
print(ctrl, digits = 3, row.names = FALSE); print(per_pop, digits = 3, row.names = FALSE)
