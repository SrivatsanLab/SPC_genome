# Bootstrap support for the population-level SBS trees of
# K562_sbs_tree_rebuild.R (original 700 sites; pruned of germline / chrM /
# unplaced sites, 333; reselected on unflagged sites, 700).
#
# Tree, as in the rebuild and sc_PolE_novaseq/analysis.ipynb: per population a
# site is present if any of its cells has it; sites present in >= 2 populations;
# Hamming distance; neighbour joining; rooted on P_5. Two resamplings, 1000
# replicates each:
#   sites - sites resampled with replacement (standard phylogenetic bootstrap,
#           as in the notebook)
#   cells - cells resampled with replacement within each population, then
#           presence recomputed (sensitivity to which cells were sequenced)
#   equal - each population subsampled without replacement to the size of the
#           smallest (P_1, 88 cells), so presence ("any cell has it") is not
#           driven by how many cells a population has
# K562_RESAMPLING picks which of these to run (default "sites,cells,equal").
# Root: K562_TREE_ROOT = P_5 (default, as in the notebook) or WT - an all-zero
# profile (the reference genome, carrying none of the sites) added to the
# population matrix and used as the outgroup, so the data decide which
# population is most basal; WT outputs go to bootstrap_rootWT/.
# Reported: support for each clade of the observed tree (ape::prop.clades), the
# frequency of each rooted topology, and of the P_0/P_1/P_2 groupings; trees
# drawn with support values (PDF, consensus-tree style).
# The per-cell matrix for the sites in the three sets is cached after the first
# read of sc_PolE_filtered_no_singletons.h5ad.
# K562_ALL_CELLS=1 uses all 1441 cells with SNV calls instead of the notebook's
# 1000: the sites are read from the per-region AnnData chunks
# (sc_PolE_novaseq/vcf_chunks/*.h5ad, which reproduce the 1000-cell matrix for
# those cells), populations come from the round-1 barcode (analysis.ipynb), and
# over-calling cells (> 10x the median number of calls in
# results/K562_tree/sc_snv_1441/sc_snv_1441.rds; the cells analysis.ipynb
# dropped with SNP_count >= 1e6) are removed. Outputs go to *_1441 folders.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT", "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
tree_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_rebuild")
tree_root <- Sys.getenv("K562_TREE_ROOT", "P_5")
kinds <- strsplit(Sys.getenv("K562_RESAMPLING", "sites,cells,equal"), ",")[[1]]
all_cells <- Sys.getenv("K562_ALL_CELLS", "0") == "1"
output_dir <- file.path(tree_dir, paste0(if (tree_root == "WT") "bootstrap_rootWT" else "bootstrap", if (all_cells) "_1441" else ""))
cache_dir <- file.path(tree_dir, "bootstrap")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
n_boot <- 1000
set.seed(1)
ordered_pops <- c("P_5", "P_4", "P_3", "P_2", "P_1", "P_0")
pop_colors <- c("P_0" = "#260091", "P_1" = "#1E90FF", "P_2" = "#FFDB58",
                "P_3" = "#FF9D71", "P_4" = "#FF1B5E", "P_5" = "#E3E6E6", "WT" = "black")

# ---- per-cell matrix for the sites in the three sets (cached) --------------------------
set_sites <- read.csv(file.path(tree_dir, "sbs_tree_rebuild_sites.csv"))
cache <- file.path(cache_dir, if (all_cells) "cell_site_matrix_1441.rds" else "cell_site_matrix.rds")
pop_wells <- list(P_0 = c("CTGTGAAC", "ACATCGTC", "TTGGATGC", "GTTACGGT"), P_1 = c("TGGTCTCA", "TCCAGGAT", "CATACCGT", "TTCCACTC"),
                  P_2 = c("GTGACTCT", "GGCATCAA", "AGATACCG", "TTGATGGC"), P_3 = c("CGACCTAT", "TCCTCAAC", "TTGACCAC", "GATGGTCT"),
                  P_4 = c("CGGTAAGT", "TGACAGTG", "CACTGACA", "GTCGGTAA"), P_5 = c("GTAACCGA", "GTCCGATT", "GACAGCAT", "CGGTTGAT"))
pop_from_barcode <- function(cells) {
  out <- rep(NA_character_, length(cells))
  for (q in names(pop_wells)) for (w in pop_wells[[q]]) out[is.na(out) & grepl(w, cells, fixed = TRUE)] <- q
  out
}
if (file.exists(cache)) {
  cm <- readRDS(cache)
} else if (all_cells) {
  keep <- unique(set_sites$site)
  chunks <- sort(list.files(file.path(sc_dir, "vcf_chunks"), "\\.h5ad$", full.names = TRUE))
  cells <- as.character(h5read(chunks[1], "/obs/_index"))
  read_chunk <- function(f) {
    v <- as.character(h5read(f, "/var/_index")); j <- which(v %in% keep)
    if (length(j) == 0) return(NULL)
    i <- match(cells, as.character(h5read(f, "/obs/_index"))); stopifnot(!anyNA(i))
    sh <- h5readAttributes(f, "/X")$shape
    X <- as(sparseMatrix(j = as.integer(h5read(f, "/X/indices")) + 1L, p = as.integer(h5read(f, "/X/indptr")),
                         x = as.numeric(h5read(f, "/X/data")), dims = sh, index1 = TRUE, repr = "R"), "CsparseMatrix")
    h5closeAll()
    m <- X[i, j, drop = FALSE]; colnames(m) <- v[j]; m
  }
  parts <- parallel::mclapply(chunks, read_chunk, mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
  X <- do.call(cbind, Filter(Negate(is.null), parts))
  stopifnot(setequal(colnames(X), keep))
  X <- X[, keep]; X@x[] <- 1; rownames(X) <- cells
  n_calls <- rowSums(readRDS(file.path(project_root, "results/K562_tree/sc_snv_1441/sc_snv_1441.rds"))$X[cells, ])
  over <- n_calls > 10 * median(n_calls)
  cat(sprintf("dropping %d over-calling cells\n", sum(over)))
  cm <- list(X = X[!over, ], sites = keep, cells = cells[!over], pop = pop_from_barcode(cells[!over]))
  stopifnot(!anyNA(cm$pop))
  saveRDS(cm, cache)
  cm <- readRDS(cache)
} else {
  h5 <- file.path(sc_dir, "sc_PolE_filtered_no_singletons.h5ad")
  shape <- h5readAttributes(h5, "/X")$shape
  X <- sparseMatrix(j = as.integer(h5read(h5, "/X/indices")) + 1L, p = as.integer(h5read(h5, "/X/indptr")),
                    x = as.numeric(h5read(h5, "/X/data")), dims = shape, index1 = TRUE, repr = "R")
  X <- as(X, "CsparseMatrix")
  sites <- as.character(h5read(h5, "/var/_index"))
  cells <- as.character(h5read(h5, "/obs/_index"))
  pop <- as.character(h5read(h5, "/obs/pop/categories"))[as.integer(h5read(h5, "/obs/pop/codes")) + 1]
  h5closeAll()
  keep <- unique(set_sites$site)
  cm <- list(X = X[, match(keep, sites)], sites = keep, cells = cells, pop = pop)
  colnames(cm$X) <- keep
  saveRDS(cm, cache)
}
cat(sprintf("%d cells x %d sites in the three sets\n", nrow(cm$X), ncol(cm$X)))

# ---- tree from a presence matrix --------------------------------------------------------
presence <- function(Xs, pop) t(sapply(ordered_pops, function(q) as.numeric(colSums(Xs[pop == q, , drop = FALSE]) > 0)))
nj_tree <- function(G) {
  G <- G[, colSums(G) >= 2, drop = FALSE]
  D <- as.matrix(dist(G, method = "manhattan")) / ncol(G)
  if (tree_root == "WT") {
    G <- rbind(G, WT = 0)
    D <- as.matrix(dist(G, method = "manhattan")) / ncol(G)
  }
  root(nj(as.dist(D)), outgroup = tree_root, resolve.root = TRUE)
}
topo_string <- function(tr) {
  # canonical rooted topology string: sort children alphabetically
  tr$edge.length <- NULL
  canon <- function(node) {
    kids <- tr$edge[tr$edge[, 1] == node, 2]
    labs <- sapply(kids, function(k) if (k <= Ntip(tr)) tr$tip.label[k] else canon(k))
    paste0("(", paste(sort(labs), collapse = ","), ")")
  }
  canon(Ntip(tr) + 1)
}
p012_group <- function(tr) {
  # which pair among P_0, P_1, P_2 is joined first (smallest clade with two of them)
  pairs <- list("P_0+P_1" = c("P_0", "P_1"), "P_0+P_2" = c("P_0", "P_2"), "P_1+P_2" = c("P_1", "P_2"))
  sz <- sapply(pairs, function(pr) length(extract.clade(tr, getMRCA(tr, pr))$tip.label))
  if (sum(sz == min(sz)) > 1) "unresolved" else names(which.min(sz))
}

plot_support <- function(tr, support, title) {
  tr$node.label <- c("", round(100 * support[-1]))
  tr$edge.length <- pmax(tr$edge.length, 0)   # NJ rounding (~1e-17)
  if (tree_root == "WT") {
    # WT (all zeros) is the ancestor: put it at depth 0 and move the whole
    # root-to-WT length onto the ingroup stem
    r <- Ntip(tr) + 1L; o <- which(tr$tip.label == "WT")
    e_out <- which(tr$edge[, 1] == r & tr$edge[, 2] == o); e_in <- which(tr$edge[, 1] == r & tr$edge[, 2] != o)
    tr$edge.length[e_in] <- tr$edge.length[e_in] + tr$edge.length[e_out]; tr$edge.length[e_out] <- 0
  }
  depth <- max(node.depth.edgelength(tr)); breaks <- seq(0, ceiling(depth * 10) / 10, by = 0.1)
  d <- ggtree(tr)$data; tips <- d[d$isTip, ]; nodes <- d[!d$isTip & d$label != "", ]
  p <- ggtree(tr, layout = "roundrect", size = 1) +
    geom_point(data = tips, aes(x = x, y = y, fill = label), shape = 21, size = 6, stroke = 1, color = "black") +
    geom_text(data = nodes, aes(x = x, y = y, label = label), vjust = -0.6, hjust = 1.2, size = 3) +
    scale_fill_manual(values = pop_colors) + scale_x_reverse(breaks = breaks) + coord_flip() + theme_tree2() +
    theme(legend.position = "none", axis.text.y = element_text(size = 10, color = "black"), axis.line.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          plot.title = element_text(size = 9)) + ggtitle(title)
  p$layers <- c(geom_vline(xintercept = breaks, color = "grey", linetype = "dashed", linewidth = 0.5, alpha = 0.6), p$layers)
  p
}

# ---- bootstrap per site set --------------------------------------------------------------
summary_rows <- list(); topo_rows <- list(); plots <- list()
for (k in c("original", "pruned", "reselected")) {
  s <- set_sites$site[set_sites$set == k]                       # with duplicates, as selected
  Xs <- cm$X[, s, drop = FALSE]
  G <- presence(Xs, cm$pop)
  obs <- nj_tree(G)
  # site bootstrap: resample columns of the presence matrix
  site_trees <- cell_trees <- equal_trees <- list()
  if ("sites" %in% kinds)
    site_trees <- lapply(seq_len(n_boot), function(b) nj_tree(G[, sample(ncol(G), replace = TRUE), drop = FALSE]))
  # cell bootstrap: resample cells within each population
  if ("cells" %in% kinds) cell_trees <- lapply(seq_len(n_boot), function(b) {
    idx <- unlist(lapply(ordered_pops, function(q) { w <- which(cm$pop == q); w[sample.int(length(w), replace = TRUE)] }))
    nj_tree(presence(Xs[idx, , drop = FALSE], cm$pop[idx]))
  })
  n_min <- min(table(cm$pop))
  if ("equal" %in% kinds) equal_trees <- lapply(seq_len(n_boot), function(b) {
    idx <- unlist(lapply(ordered_pops, function(q) { w <- which(cm$pop == q); w[sample.int(length(w), n_min)] }))
    nj_tree(presence(Xs[idx, , drop = FALSE], cm$pop[idx]))
  })
  class(site_trees) <- class(cell_trees) <- class(equal_trees) <- "multiPhylo"
  for (kind in kinds) {
    trs <- switch(kind, sites = site_trees, cells = cell_trees, equal = equal_trees)
    supp <- prop.clades(obs, trs, rooted = TRUE) / n_boot
    clades <- sapply(seq_len(obs$Nnode), function(i) paste(sort(extract.clade(obs, Ntip(obs) + i)$tip.label), collapse = ","))
    obs_topo <- topo_string(obs)
    topos <- table(sapply(trs, topo_string))
    grp <- table(factor(sapply(trs, p012_group), levels = c("P_0+P_1", "P_0+P_2", "P_1+P_2", "unresolved")))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(set = k, resampling = kind, clade = clades,
                                                           support = round(supp, 3))
    topo_rows[[length(topo_rows) + 1]] <- data.frame(set = k, resampling = kind, topology = names(topos),
                                                     frequency = as.vector(topos) / n_boot,
                                                     observed = names(topos) == obs_topo)
    topo_rows[[length(topo_rows)]] <- topo_rows[[length(topo_rows)]][order(-topo_rows[[length(topo_rows)]]$frequency), ]
    cat(sprintf("\n%s, %s bootstrap: observed topology %s in %.1f%%; first-joined pair among P_0/P_1/P_2: %s\n",
                k, kind, obs_topo, 100 * mean(sapply(trs, topo_string) == obs_topo),
                paste(sprintf("%s %.1f%%", names(grp), 100 * as.vector(grp) / n_boot), collapse = ", ")))
    print(data.frame(clade = clades, support = round(supp, 3)), row.names = FALSE)
    plots[[paste(k, kind)]] <- plot_support(obs, supp, sprintf("%s - %s bootstrap", k, kind))
  }
  write.tree(obs, file.path(output_dir, paste0("population_tree_", k, ".newick")))
  if ("sites" %in% kinds) write.tree(site_trees, file.path(output_dir, paste0("bootstrap_sites_", k, ".newick")))
  if ("cells" %in% kinds) write.tree(cell_trees, file.path(output_dir, paste0("bootstrap_cells_", k, ".newick")))
  if ("equal" %in% kinds) write.tree(equal_trees, file.path(output_dir, paste0("subsample_equal_", k, ".newick")))
}
write.csv(do.call(rbind, summary_rows), file.path(output_dir, "clade_support.csv"), row.names = FALSE)
write.csv(do.call(rbind, topo_rows), file.path(output_dir, "topology_frequencies.csv"), row.names = FALSE)
for (nm in names(plots)) ggsave(file.path(output_dir, paste0("population_tree_", sub(" ", "_", nm), "_support.pdf")),
                                plots[[nm]], width = 3, height = 3)
combined <- wrap_plots(plots, ncol = length(kinds), byrow = TRUE)
ggsave(file.path(output_dir, "population_trees_support.pdf"), combined, width = 3.2 * length(kinds), height = 9.5)
ggsave(file.path(output_dir, "population_trees_support.png"), combined, width = 3.2 * length(kinds), height = 9.5, dpi = 200, bg = "white")
cat("\nOutputs written to", output_dir, "\n")
