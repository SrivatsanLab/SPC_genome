# Recover the SBS sites used to build the per-cell K562 SBS tree
# (sc_og_test.newick).
#
# The tree was built in sc_PolE_novaseq/analysis.ipynb (cell 143): Hamming
# distance over a selected set of sites in sc_PolE_filtered_no_singletons.h5ad
# (1000 cells x 3,095,002 sites, 0/1 = variant seen in the cell), plus an
# all-zero K562 root, neighbour joining (skbio), rooted on K562. The site list
# was not saved, and the notebook has three selection cells run in several
# sessions (121: n = 5000, 126: n = 100, 147: n = 1000), so each is rebuilt
# here and the one that reproduces the tree is kept. Selection, as in the
# notebook: the n most frequent sites, then for each population in the order
# P_5, P_4, ..., P_0 (excluding those top sites) the n sites most enriched over
# the earlier populations (P_5: highest count; others: odds ratio with a 0.5
# continuity correction). Sites chosen for more than one population appear
# more than once and are weighted accordingly.
#
# Writes paper_figures/output/K562_sbs_tree_sites/: tree_match.csv (patristic
# distance correlation and topological distance per candidate) and
# sbs_tree_sites.csv for the matching candidate (site, chrom, pos, ref, alt,
# weight, selected_for, n_cells; pos 1-based - the site names in the h5ad are
# 0-based).

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(ape)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT", "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
h5 <- file.path(sc_dir, "sc_PolE_filtered_no_singletons.h5ad")
output_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_sites/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- matrix ----------------------------------------------------------------------

shape <- h5readAttributes(h5, "/X")$shape
X <- sparseMatrix(j = as.integer(h5read(h5, "/X/indices")) + 1L,
                  p = as.integer(h5read(h5, "/X/indptr")),
                  x = as.numeric(h5read(h5, "/X/data")),
                  dims = shape, index1 = TRUE, repr = "R")
X <- as(X, "CsparseMatrix")
cells <- as.character(h5read(h5, "/obs/_index"))
sites <- as.character(h5read(h5, "/var/_index"))
pop <- as.character(h5read(h5, "/obs/pop/categories"))[as.integer(h5read(h5, "/obs/pop/codes")) + 1]
h5closeAll()
cat(sprintf("%d cells x %d sites, %d nonzero\n", nrow(X), ncol(X), length(X@x)))

# ---- selection, as in analysis.ipynb ------------------------------------------------

ordered_pops <- c("P_5", "P_4", "P_3", "P_2", "P_1", "P_0")
col_tot <- colSums(X)
select_features <- function(n) {
  top_all <- order(col_tot, decreasing = TRUE)[seq_len(n)]
  rest <- setdiff(seq_len(ncol(X)), top_all)
  features <- data.frame(idx = top_all, selected_for = "top_all")
  for (i in seq_along(ordered_pops)) {
    q <- ordered_pops[i]
    in_q <- pop == q
    count_with <- colSums(X[in_q, rest, drop = FALSE])
    if (i == 1) {
      score <- count_with
    } else {
      prev <- pop %in% ordered_pops[seq_len(i - 1)]
      prev_with <- colSums(X[prev, rest, drop = FALSE])
      a <- count_with + 0.5; b <- sum(in_q) - count_with + 0.5
      c <- prev_with + 0.5; d <- sum(prev) - prev_with + 0.5
      score <- (a * d) / (b * c)
    }
    # pandas sort_values(ascending = False) is a stable quicksort; ties keep
    # the original column order, as order(-score) does here
    top <- rest[order(-score)[seq_len(n)]]
    features <- rbind(features, data.frame(idx = top, selected_for = q))
  }
  features
}

hamming_tree <- function(features) {
  w <- tabulate(features$idx, nbins = ncol(X))
  used <- which(w > 0)
  B <- X[, used, drop = FALSE]; wu <- w[used]
  n_i <- as.vector(B %*% wu)
  S <- as.matrix(B %*% Diagonal(x = wu) %*% t(B))
  f <- sum(wu)
  D <- (outer(n_i, n_i, "+") - 2 * S) / f
  D <- rbind(cbind(D, n_i / f), c(n_i / f, 0))     # K562 root: all zeros
  dimnames(D) <- list(c(cells, "K562"), c(cells, "K562"))
  root(nj(as.dist(D)), outgroup = "K562", resolve.root = TRUE)
}

ref_tree <- read.tree(file.path(sc_dir, "results/sc_og_test.newick"))
ref_coph <- cophenetic.phylo(ref_tree)[cells, cells]
ut <- upper.tri(ref_coph)

match <- list(); sets <- list()
for (n in c(100, 1000, 5000)) {
  ft <- select_features(n)
  tr <- hamming_tree(ft)
  coph <- cophenetic.phylo(tr)[cells, cells]
  match[[length(match) + 1]] <- data.frame(
    n_per_set = n, n_selected = nrow(ft), n_unique_sites = length(unique(ft$idx)),
    patristic_r = cor(coph[ut], ref_coph[ut]),
    rf_distance = dist.topo(unroot(drop.tip(tr, "K562")), unroot(drop.tip(ref_tree, "K562"))))
  sets[[as.character(n)]] <- ft
  print(tail(match[[length(match)]], 1))
}
match <- do.call(rbind, match)
write.csv(match, file.path(output_dir, "tree_match.csv"), row.names = FALSE)

best <- match$n_per_set[which.max(match$patristic_r)]
ft <- sets[[as.character(best)]]
agg <- aggregate(selected_for ~ idx, ft, function(s) paste(unique(s), collapse = ";"))
agg$weight <- tabulate(ft$idx, nbins = ncol(X))[agg$idx]
parts <- do.call(rbind, strsplit(sites[agg$idx], "[->]"))
# site names carry 0-based positions (the reference base matches at name + 1
# for every site); pos below is 1-based
out <- data.frame(site = sites[agg$idx], chrom = parts[, 1], pos = as.integer(parts[, 2]) + 1L,
                  ref = parts[, 3], alt = parts[, 4], weight = agg$weight,
                  selected_for = agg$selected_for, n_cells = col_tot[agg$idx])
write.csv(out, file.path(output_dir, "sbs_tree_sites.csv"), row.names = FALSE)
# all sites in the filtered no-singleton set, for comparison
parts_all <- do.call(rbind, strsplit(sites, "[->]"))
write.csv(data.frame(site = sites, chrom = parts_all[, 1], pos = as.integer(parts_all[, 2]) + 1L,
                     ref = parts_all[, 3], alt = parts_all[, 4], n_cells = col_tot),
          gzfile(file.path(output_dir, "all_called_sites.csv.gz")), row.names = FALSE)
cat(sprintf("best match: n = %d per set (%d unique sites, patristic r = %.4f, RF = %g)\n",
            best, nrow(out), max(match$patristic_r), match$rf_distance[match$n_per_set == best]))
