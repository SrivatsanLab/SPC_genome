#!/usr/bin/env Rscript
# MP tree search + bootstrap under Camin-Sokal (irreversible gain) presence-only
# encoding. Runs on one (worm, matrix) at a time so it composes cleanly with a
# Python-side parallel launcher.
#
# Usage:
#   mp_run.R --tsv <matrix.tsv> --out-dir <out>  \
#            --tag <real|covmask|shuffled>       \
#            --bootstrap 100 --k 10 --seed 0     \
#            --threads 4
#
# Inputs:
#   matrix.tsv : cell_id \t char_string, chars in {0,1,?} (we only emit 1,?)
#
# Outputs at <out>/:
#   <tag>.mp.newick        best MP tree
#   <tag>.strict.newick    strict consensus of all MP trees found
#   <tag>.majority.newick  majority-rule consensus of bootstrap replicates,
#                          internal-node labels = support (0-100)
#   <tag>.bootstrap.newick multi-Newick of all bootstrap trees
#   <tag>.stats.json       {n_mpts, mp_score, ci, ri, hi, n_informative_chars,
#                           n_cells, n_variants, elapsed_search_s,
#                           elapsed_bootstrap_s}
#   <tag>.log              stdout+stderr

suppressPackageStartupMessages({
  library(optparse)
  library(phangorn)
  library(ape)
  library(jsonlite)
  library(parallel)
})

opt_list <- list(
  make_option("--tsv", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--tag", type = "character"),
  make_option("--bootstrap", type = "integer", default = 100L),
  make_option("--k", type = "integer", default = 10L,
              help = "pratchet ratchet iterations"),
  make_option("--seed", type = "integer", default = 0L),
  make_option("--threads", type = "integer", default = 1L),
  make_option("--gain-cost", type = "double", default = 1.0,
              help = "cost of 0->1 transition. Default 1 (Fitch/CS). Set high (e.g. 100) for Dollo-approximation."),
  make_option("--cs-cost", type = "double", default = 1e9,
              help = "cost of 1->0 transition. 1e9 = strict CS (irreversible); 1 = Fitch; low values allow reversal.")
)
opt <- parse_args(OptionParser(option_list = opt_list))

set.seed(opt$seed)
options(mc.cores = opt$threads)

stopifnot(!is.null(opt$tsv), !is.null(opt$`out-dir`), !is.null(opt$tag))

out_dir <- opt$`out-dir`
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(out_dir, paste0(opt$tag, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(NULL), add = TRUE)

cat("mp_run.R", format(Sys.time()), "\n")
cat("tsv:", opt$tsv, "\n")
cat("tag:", opt$tag, " bootstrap:", opt$bootstrap,
    " k:", opt$k, " threads:", opt$threads,
    " gain-cost(0->1):", opt$`gain-cost`,
    " cs-cost(1->0):", opt$`cs-cost`, "\n")

# ---- load matrix -------------------------------------------------------------
tsv <- read.table(opt$tsv, sep = "\t", header = TRUE,
                  colClasses = c("character", "character"), quote = "")
cells <- tsv$cell_id
char_mat <- do.call(rbind, strsplit(tsv$chars, split = ""))
rownames(char_mat) <- cells
n_cells <- nrow(char_mat)
n_vars  <- ncol(char_mat)
cat("loaded ", n_cells, " cells x ", n_vars, " chars\n", sep = "")

# Reject constant columns (all-?, all-1, all-0) before phyDat — they are
# uninformative and phangorn drops them anyway, but pratchet-with-sankoff can
# choke on them.
u <- apply(char_mat, 2, function(x) length(unique(x)))
inf_col <- u > 1
if (any(!inf_col)) {
  cat("removing ", sum(!inf_col), " constant columns; ",
      sum(inf_col), " remain\n", sep = "")
  char_mat <- char_mat[, inf_col, drop = FALSE]
}
n_informative <- ncol(char_mat)
stopifnot(n_informative >= 3)  # phangorn needs some signal

# ---- build phyDat with USER type + ? ambiguity -------------------------------
# We restrict to 2 states {0,1}; "?" resolves to {0,1} (both states).
dat <- phyDat(char_mat, type = "USER", levels = c("0", "1"),
              ambiguity = c("?"))
attr(dat, "contrast")           # sanity
cat("phyDat: ", length(dat), " taxa, ", attr(dat, "nr"), " site patterns\n",
    sep = "")

# ---- Parsimony method: Sankoff Camin-Sokal ---------------------------------
# We use the naive 0/1/? encoding (DP=0 → '?', DP>0 & AD=0 → '0', AD≥1 → '1').
# Observed '0' states anchor Fitch/Sankoff to real bipartitions. Sankoff with
# an asymmetric cost matrix (0→1 free, 1→0 forbidden) implements Camin-Sokal:
# gains only, no reversions. If phangorn's Sankoff-ratchet crashes on small
# subtrees, fall back to Fitch (still meaningful on the naive encoding, just
# symmetric on state transitions).
cs_cost <- matrix(c(0, opt$`gain-cost`, opt$`cs-cost`, 0),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(c("0", "1"), c("0", "1")))
cat("Sankoff cost matrix:\n")
print(cs_cost)

# ---- MP search ---------------------------------------------------------------
t0 <- Sys.time()
best_tree <- tryCatch({
    pratchet(dat, k = opt$k, method = "sankoff", cost = cs_cost,
             trace = 0, all = TRUE)
}, error = function(e) {
    cat("Sankoff ratchet failed (", conditionMessage(e),
        "); falling back to Fitch.\n", sep = "")
    pratchet(dat, k = opt$k, method = "fitch", trace = 0, all = TRUE)
})
elapsed_search <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (inherits(best_tree, "phylo")) {
  best_tree <- list(best_tree)
  class(best_tree) <- "multiPhylo"
}
n_mpts <- length(best_tree)
# Score on Sankoff C-S if available, else Fitch.
mp_score <- tryCatch(
    as.numeric(parsimony(best_tree[[1]], dat, method = "sankoff", cost = cs_cost)),
    error = function(e) as.numeric(parsimony(best_tree[[1]], dat, method = "fitch"))
)
cat("search: ", n_mpts, " MPT(s) in ", round(elapsed_search, 1), "s ",
    "(score = ", mp_score, ")\n", sep = "")

# Strict consensus.
strict_tree <- consensus(best_tree, p = 1.0)

# CI/RI/HI — computed via phangorn on the first MP tree with Fitch (standard
# formulation). Camin-Sokal + missing data makes CI/RI conservative but the
# statistic is still comparable across the three matrices.
ci_val <- as.numeric(CI(best_tree[[1]], dat, sitewise = FALSE))
ri_val <- as.numeric(RI(best_tree[[1]], dat, sitewise = FALSE))
hi_val <- 1 - ci_val

# ---- Bootstrap ---------------------------------------------------------------
t0 <- Sys.time()
if (opt$bootstrap > 0) {
  bs_fn <- function(x) {
    tr <- tryCatch({
        pratchet(x, k = max(3L, opt$k %/% 2L),
                 method = "sankoff", cost = cs_cost, trace = 0, all = FALSE)
    }, error = function(e) {
        pratchet(x, k = max(3L, opt$k %/% 2L),
                 method = "fitch", trace = 0, all = FALSE)
    })
    if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
    tr
  }
  bs_trees <- bootstrap.phyDat(dat, bs_fn, bs = opt$bootstrap,
                               multicore = opt$threads > 1, mc.cores = opt$threads)
  elapsed_boot <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  majority_tree <- consensus(bs_trees, p = 0.5)
  # Support values on the best_tree topology (plotBS returns tree with support).
  best_with_support <- plotBS(best_tree[[1]], bs_trees, type = "n",
                              method = "FBP")
  cat("bootstrap: ", length(bs_trees), " reps in ",
      round(elapsed_boot, 1), "s\n", sep = "")
} else {
  bs_trees <- NULL
  majority_tree <- strict_tree
  best_with_support <- best_tree[[1]]
  elapsed_boot <- 0
}

# ---- Write outputs -----------------------------------------------------------
write.tree(best_tree[[1]],  file.path(out_dir, paste0(opt$tag, ".mp.newick")))
write.tree(strict_tree,     file.path(out_dir, paste0(opt$tag, ".strict.newick")))
write.tree(majority_tree,   file.path(out_dir, paste0(opt$tag, ".majority.newick")))
write.tree(best_with_support,
           file.path(out_dir, paste0(opt$tag, ".mp_support.newick")))
if (!is.null(bs_trees)) {
  write.tree(bs_trees,
             file.path(out_dir, paste0(opt$tag, ".bootstrap.newick")))
}

# Resolution: fraction of internal nodes in the consensus (rooted at the
# center) that are resolved (i.e., not part of a polytomy). For a rooted tree
# with n tips, a fully resolved binary tree has n-1 internal nodes.
resolution <- function(tree) {
  if (!inherits(tree, "phylo")) return(NA_real_)
  n_tips <- length(tree$tip.label)
  # count internal edges (not tip edges)
  n_internal <- tree$Nnode
  # for an unrooted binary tree Nnode = n-2; for rooted = n-1
  # phangorn consensus() returns rooted trees.
  max_internal <- n_tips - 1L
  n_internal / max_internal
}

stats <- list(
  worm = NA,  # populated by launcher
  tag = opt$tag,
  n_cells = n_cells,
  n_variants = n_vars,
  n_informative_chars = n_informative,
  n_mpts = n_mpts,
  mp_score = mp_score,
  ci = ci_val,
  ri = ri_val,
  hi = hi_val,
  resolution_strict = resolution(strict_tree),
  resolution_majority = resolution(majority_tree),
  bootstrap_reps = if (is.null(bs_trees)) 0L else length(bs_trees),
  elapsed_search_s = elapsed_search,
  elapsed_bootstrap_s = elapsed_boot,
  cs_reversal_cost = opt$`cs-cost`,
  gain_cost = opt$`gain-cost`,
  ratchet_k = opt$k,
  seed = opt$seed
)
write_json(stats,
           file.path(out_dir, paste0(opt$tag, ".stats.json")),
           auto_unbox = TRUE, pretty = TRUE)

cat("done. wrote ", opt$tag, ".{mp,strict,majority,mp_support,bootstrap}.newick, ",
    opt$tag, ".stats.json\n", sep = "")
