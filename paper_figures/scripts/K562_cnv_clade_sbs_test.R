# Do small, clearly called CNV clades hold together in the SBS tree?
#
# The per-cell CNV tree cannot resolve close relatives, but a large event shared
# by a few cells is readily distinguishable: at ~2.4 SD per copy per 1 Mb bin,
# a >= 5 Mb event is ~5 SD. If such a group is a real lineage, its members
# should also sit close together in the per-cell SBS tree (sc_og_test.newick).
#
# CNV clades: from the GC-corrected event matrix (K562_cnv_event_tree.R,
# K562_CNV_TAG=_gc), gain/loss features of >= min_mb, carried by
# min_cells..max_cells high-quality cells (main-peak ploidy and Pearson r >=
# min_ref_r with the K562 reference). Features whose carrier sets overlap with
# Jaccard >= merge_jaccard are merged into one clade (carriers = cells in any).
#
# Test per clade: mean over member pairs of the size of the smallest SBS-tree
# clade containing the pair (topology only, not branch length), against
# n_perm random sets of the same size from the same high-quality cells.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(dplyr)
})

sc_dir <- Sys.getenv("K562_SC_TREES",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/results")
tag <- Sys.getenv("K562_CNV_TAG", "_gc")
tree_dir <- file.path(project_root, "results/K562_tree/sc_trees")
ref_csv <- file.path(project_root, "paper_figures/output/K562_cnv_vs_reference/per_cell_agreement.csv")
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_cnv_clade_sbs_test", tag, "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_mb <- 5
min_cells <- 3
max_cells <- 50
min_ref_r <- 0.7
merge_jaccard <- 0.8
n_perm <- 999
set.seed(1)

# ---- high-quality cells and CNV clades ---------------------------------------

em <- readRDS(file.path(tree_dir, paste0("cnv_event_matrix", tag, ".rds")))
feat <- read.csv(file.path(tree_dir, paste0("cnv_event_features", tag, ".csv")))
X <- em$features
stopifnot(nrow(X) == nrow(feat))

ref <- read.csv(ref_csv) %>% filter(run == "new", bins == "All bins")
ref_r <- setNames(ref$pearson, ref$cell)
good <- colnames(X)[ref_r[colnames(X)] >= min_ref_r]
X <- X[, good, drop = FALSE]
cat(sprintf("%d high-quality cells (main-peak ploidy, reference r >= %.2f)\n", length(good), min_ref_r))

feat$mb <- (feat$end - feat$start + 1) / 1e6
feat$carriers <- rowSums(X)
cand <- which(feat$mb >= min_mb & feat$carriers >= min_cells & feat$carriers <= max_cells)
cat(sprintf("%d gain/loss features >= %g Mb carried by %d-%d cells\n", length(cand), min_mb, min_cells, max_cells))

# greedy merge of features with near-identical carrier sets, largest first
cand <- cand[order(-feat$mb[cand])]
sets <- lapply(cand, function(i) colnames(X)[X[i, ] == 1])
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
clade_of <- rep(NA_integer_, length(cand))
for (i in seq_along(cand)) {
  if (!is.na(clade_of[i])) next
  clade_of[i] <- i
  for (j in seq_along(cand)) {
    if (is.na(clade_of[j]) && jac(sets[[i]], sets[[j]]) >= merge_jaccard) clade_of[j] <- i
  }
}
clades <- lapply(split(seq_along(cand), clade_of), function(ix) {
  list(events = with(feat[cand[ix], ], sprintf("%s:%d-%dMb %s", chr, round(start / 1e6), round(end / 1e6), type)),
       mb = sum(feat$mb[cand[ix]]),
       cells = sort(unique(unlist(sets[ix]))))
})
cat(sprintf("%d CNV clades after merging (Jaccard >= %.1f)\n", length(clades), merge_jaccard))

# ---- SBS tree: pairwise smallest-shared-clade size ---------------------------

sbs <- read.tree(file.path(sc_dir, "sc_og_test.newick"))
sbs <- drop.tip(sbs, setdiff(sbs$tip.label, good))
n <- Ntip(sbs)
m <- mrca(sbs)
sizes <- lengths(Descendants(sbs, seq_len(n + sbs$Nnode), type = "tips"))
S <- matrix(sizes[m], nrow(m), dimnames = dimnames(m))
S <- S[good, good]

meta <- read.csv(file.path(sc_dir, "full_meta.csv"))
pop <- setNames(meta$pop, meta$X)

mean_pair <- function(cells) { s <- S[cells, cells]; mean(s[lower.tri(s)]) }

res <- bind_rows(lapply(seq_along(clades), function(k) {
  cl <- clades[[k]]; sz <- length(cl$cells)
  obs <- mean_pair(cl$cells)
  null <- replicate(n_perm, mean_pair(sample(good, sz)))
  pops <- table(factor(pop[cl$cells], levels = sort(unique(meta$pop))))
  data.frame(
    clade = k, events = paste(cl$events, collapse = "; "), event_mb = cl$mb,
    n_cells = sz,
    sbs_mean_pair_clade = obs, null_mean = mean(null),
    null_lo = unname(quantile(null, 0.025)),
    p = (1 + sum(null <= obs)) / (n_perm + 1),
    pops = paste(pops, collapse = "/"),
    top_pop_frac = max(pops) / sz
  )
})) %>% arrange(p)
res$q_bh <- p.adjust(res$p, "BH")

write.csv(res, file.path(output_dir, "cnv_clades_in_sbs_tree.csv"), row.names = FALSE)
writeLines(unlist(lapply(seq_along(clades), function(k) paste(k, paste(clades[[k]]$cells, collapse = ",")))),
           file.path(output_dir, "cnv_clade_members.txt"))

cat(sprintf("\nclades with p < 0.05: %d of %d (%.1f expected by chance); BH q < 0.1: %d\n",
            sum(res$p < 0.05), nrow(res), 0.05 * nrow(res), sum(res$q_bh < 0.1)))
fisher <- pchisq(-2 * sum(log(res$p)), df = 2 * nrow(res), lower.tail = FALSE)
cat(sprintf("Fisher's combined p over all clades: %.3g\n", fisher))
cat("pops column: P0/P1/P2/P3/P4/P5 member counts\n")
print(res %>% select(clade, events, n_cells, sbs_mean_pair_clade, null_mean, p, q_bh, pops),
      row.names = FALSE, digits = 3, right = FALSE)

# ---- plot ----------------------------------------------------------------------

res %>%
  mutate(ratio = sbs_mean_pair_clade / null_mean, lo = null_lo / null_mean,
         label = reorder(sprintf("%d (%d cells)", clade, n_cells), ratio)) %>%
  ggplot(aes(ratio, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_segment(aes(x = lo, xend = 1, yend = label), color = "grey80", linewidth = 2) +
  geom_point(aes(fill = p < 0.05), shape = 21, size = 2) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"), name = "p < 0.05") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 5)) +
  xlab("Mean SBS-tree shared-clade size,\nrelative to random groups (grey: 2.5% quantile)") +
  ylab("CNV clade")
ggsave(file.path(output_dir, "cnv_clades_in_sbs_tree.svg"),
       bg = "transparent", height = max(3, 0.12 * nrow(res) + 1), width = 4.5)

cat("\nWrote results to:", output_dir, "\n")
