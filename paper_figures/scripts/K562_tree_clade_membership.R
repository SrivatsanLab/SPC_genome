# Do the per-cell SBS and CNV trees group the same cells together?
#
# Each tree is cut into k clades by topology alone: start from the root's
# children and repeatedly split the largest clade into its children until there
# are k groups. The two k-group partitions are compared with the adjusted Rand
# index (0 = chance, 1 = identical membership), against 999 random relabellings
# of the CNV tree's tips, for k = 2..50. Cross-tabulations of clade membership
# are written for small k. Only clades of >= min_clade cells count towards k
# (see partitions()), so each comparison covers the cells assigned in both.
#
# CNV tree defaults to the GC-corrected event tree; the SBS tree is pruned to
# the same cells.

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
output_dir <- file.path(project_root, paste0("paper_figures/output/K562_tree_clade_membership", output_tag, "/"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
n_perm <- 999
ks <- c(2:10, 15, 20, 30, 40, 50)
min_clade <- 10
table_ks <- c(2:6, 15, 20)

read_cells <- function(path) drop.tip(read.tree(path), "K562")
t_cnv <- read_cells(cnv_tree_path)
t_sbs <- read_cells(file.path(input_dir, "sc_og_test.newick"))
t_sbs <- drop.tip(t_sbs, setdiff(t_sbs$tip.label, t_cnv$tip.label))
stopifnot(setequal(t_sbs$tip.label, t_cnv$tip.label))
cells <- sort(t_sbs$tip.label)
n <- length(cells)

# k-clade partitions for every k up to kmax, as a cells x k membership matrix.
# Both trees are ladder-like near the root, so splitting from the root mostly
# peels off single cells; only clades of >= min_clade cells count towards k,
# and cells in smaller fragments are left unassigned (NA).
partitions <- function(tree, kmax) {
  tips_of <- function(node) tree$tip.label[unlist(Descendants(tree, node, type = "tips"))]
  groups <- Children(tree, Ntip(tree) + 1)
  out <- matrix(NA_integer_, n, kmax, dimnames = list(cells, NULL))
  repeat {
    sizes <- sapply(groups, function(g) length(tips_of(g)))
    big <- groups[sizes >= min_clade]
    k <- length(big)
    if (k >= 1 && k <= kmax && all(is.na(out[, k]))) {
      memb <- setNames(rep(NA_integer_, n), cells)
      for (g in seq_along(big)) memb[tips_of(big[g])] <- g
      out[, k] <- memb
    }
    splittable <- sizes >= min_clade & groups > Ntip(tree)
    if (k >= kmax || !any(splittable)) break
    split <- which(splittable)[which.max(sizes[splittable])]
    groups <- c(groups[-split], Children(tree, groups[split]))
  }
  out
}
p_sbs <- partitions(t_sbs, max(ks))
p_cnv <- partitions(t_cnv, max(ks))

ari <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(x) sum(x * (x - 1) / 2)
  s_ij <- comb2(tab); s_a <- comb2(rowSums(tab)); s_b <- comb2(colSums(tab))
  expected <- s_a * s_b / comb2(length(a))
  (s_ij - expected) / ((s_a + s_b) / 2 - expected)
}

res <- bind_rows(lapply(ks, function(k) {
  a <- p_sbs[, k]; b <- p_cnv[, k]
  if (all(is.na(a)) || all(is.na(b))) return(NULL)
  both <- !is.na(a) & !is.na(b)
  obs <- ari(a[both], b[both])
  # relabel the CNV tree's cells, then score the cells assigned in both
  null <- replicate(n_perm, {
    bp <- b[sample(n)]
    ok <- !is.na(a) & !is.na(bp)
    ari(a[ok], bp[ok])
  })
  data.frame(k = k, cells_compared = sum(both), ari = obs, null_mean = mean(null),
             null_hi = unname(quantile(null, 0.975)),
             p = (1 + sum(null >= obs)) / (n_perm + 1),
             sbs_largest = max(table(a)), cnv_largest = max(table(b)))
}))
write.csv(res, file.path(output_dir, "clade_membership_ari.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)

# cross-tabulations: rows = SBS clades, columns = CNV clades, cells = shared members
sink(file.path(output_dir, "clade_membership_tables.txt"))
for (k in table_ks) {
  cat(sprintf("\n== k = %d (rows: SBS-tree clades, columns: CNV-tree clades)\n", k))
  if (all(is.na(p_sbs[, k])) || all(is.na(p_cnv[, k]))) next
  tab <- table(SBS = p_sbs[, k], CNV = p_cnv[, k])
  print(addmargins(tab))
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  cat("expected under independence:\n")
  print(round(expected, 1))
  # clade pairs sharing at least twice the expected number of cells (>= 5)
  enriched <- which(tab >= 5 & tab >= 2 * expected, arr.ind = TRUE)
  if (nrow(enriched) > 0) {
    cat("enriched SBS x CNV clade pairs:\n")
    print(data.frame(sbs_clade = rownames(tab)[enriched[, 1]],
                     sbs_size = rowSums(tab)[enriched[, 1]],
                     cnv_clade = colnames(tab)[enriched[, 2]],
                     cnv_size = colSums(tab)[enriched[, 2]],
                     shared = tab[enriched],
                     expected = round(expected[enriched], 1)), row.names = FALSE)
  }
}
sink()

# ---- plot ------------------------------------------------------------------

ggplot(res, aes(k)) +
  geom_ribbon(aes(ymin = pmin(0, null_mean), ymax = null_hi), fill = "grey85") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(y = ari), linewidth = 0.6) +
  geom_point(aes(y = ari), size = 1.5) +
  theme_classic() +
  xlab("Clades per tree (k)") +
  ylab("Adjusted Rand index, SBS vs CNV")
ggsave(file.path(output_dir, "clade_membership_ari.svg"),
       bg = "transparent", height = 2.75, width = 3.75)

cat("\nWrote results to:", output_dir, "\n")
