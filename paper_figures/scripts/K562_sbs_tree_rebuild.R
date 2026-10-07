# Rebuild the K562 SBS trees without germline, mitochondrial and unplaced-contig
# sites.
#
# Inputs: sc_PolE_filtered_no_singletons.h5ad (1000 cells x 3.1 M sites, 0/1),
# per-site flags from K562_sbs_tree_germline.R (all_sites_germline.csv.gz, in
# the h5ad's site order). A site is flagged if it is germline/ancestral in bulk
# (pooled VAF >= 0.15 in the mutation-accumulation WGS) or in ENCODE K562
# germline calls, or lies on chrM or an unplaced/alt contig.
#
# Site sets:
#   original    - the 100-per-set selection of sc_PolE_novaseq/analysis.ipynb
#                 (100 most frequent + 100 most enriched per population, P_5..P_0;
#                 reproduces sc_og_test.newick, see K562_sbs_tree_sites.R)
#   pruned      - original minus flagged sites
#   reselected  - the same selection run on unflagged sites only
# Trees per set, as in the notebook:
#   single-cell - Hamming distance over the sites (+ all-zero K562 root),
#                 neighbour joining, rooted on K562
#   population  - per population, a site is present if any cell has it; sites
#                 present in >= 2 populations; Hamming, NJ, rooted on P_5
# Comparisons: patristic r and Robinson-Foulds distance with the original
# single-cell tree; same-population sister pairs (cherries) vs 10,000 label
# shuffles, with cell-cycle stage and spikiness quintile as controls; Mantel
# (Spearman) correlation with the CNV distance tree of G1/G2 cells with events
# >= 3 Mb (K562_cnv_distance_tree.R), 999 permutations.
# Writes paper_figures/output/K562_sbs_tree_rebuild/.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(ape)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT", "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
h5 <- file.path(sc_dir, "sc_PolE_filtered_no_singletons.h5ad")
output_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_rebuild/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cnv_tree_path <- file.path(project_root,
  "paper_figures/output/K562_cnv_distance_tree/aneufinder_fragment_gc_rtcorr_G1G2_ev3Mb/cnv_distance_tree.newick")
n_perm <- 10000
set.seed(1)

# ---- matrix and flags ----------------------------------------------------------------
shape <- h5readAttributes(h5, "/X")$shape
X <- sparseMatrix(j = as.integer(h5read(h5, "/X/indices")) + 1L, p = as.integer(h5read(h5, "/X/indptr")),
                  x = as.numeric(h5read(h5, "/X/data")), dims = shape, index1 = TRUE, repr = "R")
X <- as(X, "CsparseMatrix")
cells <- as.character(h5read(h5, "/obs/_index"))
sites <- as.character(h5read(h5, "/var/_index"))
pop <- as.character(h5read(h5, "/obs/pop/categories"))[as.integer(h5read(h5, "/obs/pop/codes")) + 1]
h5closeAll()
fl <- read.csv(gzfile(file.path(project_root, "paper_figures/output/K562_sbs_tree_germline/all_sites_germline.csv.gz")))
stopifnot(identical(fl$site, sites))
flagged <- fl$germline_any | fl$location != "chr1-22/X/Y"
cat(sprintf("%d of %d sites flagged (germline in bulk or ENCODE, chrM, unplaced/alt contigs)\n", sum(flagged), length(flagged)))

# ---- selection as in analysis.ipynb, optionally on a restricted universe ----------------
ordered_pops <- c("P_5", "P_4", "P_3", "P_2", "P_1", "P_0")
col_tot <- colSums(X)
select_features <- function(n, allowed = rep(TRUE, ncol(X))) {
  universe <- which(allowed)
  top_all <- universe[order(col_tot[universe], decreasing = TRUE)[seq_len(n)]]
  rest <- setdiff(universe, top_all)
  features <- data.frame(idx = top_all, selected_for = "top_all")
  for (i in seq_along(ordered_pops)) {
    in_q <- pop == ordered_pops[i]
    count_with <- colSums(X[in_q, rest, drop = FALSE])
    score <- if (i == 1) count_with else {
      prev <- pop %in% ordered_pops[seq_len(i - 1)]
      prev_with <- colSums(X[prev, rest, drop = FALSE])
      ((count_with + 0.5) * (sum(prev) - prev_with + 0.5)) / ((sum(in_q) - count_with + 0.5) * (prev_with + 0.5))
    }
    features <- rbind(features, data.frame(idx = rest[order(-score)[seq_len(n)]], selected_for = ordered_pops[i]))
  }
  features
}

sc_tree <- function(features) {
  w <- tabulate(features$idx, nbins = ncol(X)); used <- which(w > 0)
  B <- X[, used, drop = FALSE]; wu <- w[used]
  n_i <- as.vector(B %*% wu); S <- as.matrix(B %*% Diagonal(x = wu) %*% t(B)); f <- sum(wu)
  D <- (outer(n_i, n_i, "+") - 2 * S) / f
  D <- rbind(cbind(D, n_i / f), c(n_i / f, 0))
  dimnames(D) <- list(c(cells, "K562"), c(cells, "K562"))
  root(nj(as.dist(D)), outgroup = "K562", resolve.root = TRUE)
}
pop_tree <- function(features) {
  G <- t(sapply(ordered_pops, function(q) as.numeric(colSums(X[pop == q, features$idx, drop = FALSE]) > 0)))
  G <- G[, colSums(G) >= 2, drop = FALSE]
  D <- as.matrix(dist(G, method = "manhattan")) / ncol(G)
  list(tree = root(nj(as.dist(D)), outgroup = "P_5", resolve.root = TRUE), n_sites = ncol(G))
}

sets <- list(original = select_features(100))
sets$pruned <- sets$original[!flagged[sets$original$idx], ]
sets$reselected <- select_features(100, allowed = !flagged)

# ---- trees and comparisons ------------------------------------------------------------
ph <- read.csv(file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/per_cell_phase.csv"))
stage <- setNames(ifelse(ph$phase == "S", as.character(ph$s_stage), "G1/G2"), ph$cell)
spik <- setNames(ph$spikiness, ph$cell)
spik_q <- setNames(cut(spik, quantile(spik, 0:5 / 5, na.rm = TRUE), include.lowest = TRUE, labels = FALSE), names(spik))
cnv_tree <- drop.tip(read.tree(cnv_tree_path), "K562_consensus")

cherry_stats <- function(tr, lab) {
  tr <- drop.tip(tr, "K562"); n <- Ntip(tr); tips <- tr$tip.label
  kids <- split(tr$edge[, 2], tr$edge[, 1])
  ch <- do.call(rbind, Filter(Negate(is.null), lapply(kids, function(k) if (length(k) == 2 && all(k <= n)) k else NULL)))
  l <- lab[tips]; ok <- !is.na(l[ch[, 1]]) & !is.na(l[ch[, 2]]); ch <- ch[ok, , drop = FALSE]
  obs <- mean(l[ch[, 1]] == l[ch[, 2]])
  nl <- replicate(n_perm, { s <- sample(l); mean(s[ch[, 1]] == s[ch[, 2]]) })
  c(pairs = nrow(ch), same = obs, expected = mean(nl), enrichment = obs / mean(nl), p = (sum(nl >= obs) + 1) / (n_perm + 1))
}
mantel_cnv <- function(tr) {
  common <- intersect(cnv_tree$tip.label, tr$tip.label)
  a <- cophenetic.phylo(drop.tip(tr, setdiff(tr$tip.label, common)))[common, common]
  b <- cophenetic.phylo(drop.tip(cnv_tree, setdiff(cnv_tree$tip.label, common)))[common, common]
  ut <- upper.tri(a); r <- cor(a[ut], b[ut], method = "spearman")
  nl <- replicate(999, { s <- sample(common); cor(a[ut], b[s, s][ut], method = "spearman") })
  c(cells = length(common), r = r, p = (sum(nl >= r) + 1) / 1000)
}

pop_lab <- setNames(pop, cells)
orig_tree <- NULL
res <- list(); pop_trees <- list()
for (k in names(sets)) {
  ft <- sets[[k]]
  tr <- sc_tree(ft); write.tree(tr, file.path(output_dir, paste0("sc_tree_", k, ".newick")))
  if (k == "original") orig_tree <- tr
  oc <- cophenetic.phylo(orig_tree)[cells, cells]; tc <- cophenetic.phylo(tr)[cells, cells]; ut <- upper.tri(oc)
  pt <- pop_tree(ft); write.tree(pt$tree, file.path(output_dir, paste0("population_tree_", k, ".newick")))
  pop_trees[[k]] <- write.tree(pt$tree)
  sp <- cherry_stats(tr, pop_lab); st <- cherry_stats(tr, stage); sq <- cherry_stats(tr, spik_q); mc <- mantel_cnv(tr)
  res[[k]] <- data.frame(set = k, n_sites = length(unique(ft$idx)),
                         n_flagged = sum(flagged[unique(ft$idx)]),
                         patristic_r_vs_original = cor(oc[ut], tc[ut]),
                         rf_vs_original = dist.topo(unroot(drop.tip(tr, "K562")), unroot(drop.tip(orig_tree, "K562"))),
                         sister_pairs = sp["pairs"], same_pop = sp["same"], same_pop_expected = sp["expected"],
                         same_pop_enrichment = sp["enrichment"], same_pop_p = sp["p"],
                         same_stage_enrichment = st["enrichment"], same_stage_p = st["p"],
                         same_noise_enrichment = sq["enrichment"], same_noise_p = sq["p"],
                         cnv_mantel_cells = mc["cells"], cnv_mantel_r = mc["r"], cnv_mantel_p = mc["p"],
                         pop_tree_sites = pt$n_sites)
  print(res[[k]], digits = 3, row.names = FALSE)
}
res <- do.call(rbind, res)
write.csv(res, file.path(output_dir, "sbs_tree_rebuild_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(sets), function(k) data.frame(set = k, site = sites[sets[[k]]$idx],
                                                                    selected_for = sets[[k]]$selected_for))),
          file.path(output_dir, "sbs_tree_rebuild_sites.csv"), row.names = FALSE)
writeLines(paste(names(pop_trees), unlist(pop_trees)), file.path(output_dir, "population_trees.txt"))
cat("\npopulation trees:\n"); cat(paste(names(pop_trees), unlist(pop_trees)), sep = "\n")
