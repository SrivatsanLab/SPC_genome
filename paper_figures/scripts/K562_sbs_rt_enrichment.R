# Are the K562 SBS sites used for the SBS tree (and the called SBS overall)
# enriched in early, mid or late replicating parts of the genome?
#
# Sites: K562_sbs_tree_sites.R (sbs_tree_sites.csv - the sites that reproduce
# sc_og_test.newick - and all_called_sites.csv.gz - every site in
# sc_PolE_filtered_no_singletons.h5ad).
# Replication timing: K562 UW Repli-seq on hg38, 1 kb windows
# (download_replication_timing.sh, replication_timing_to_hg38.R); RT = 1 for
# earliest, 0 for latest.
# Callable genome: uniquely mappable 100-mers (Umap k100), outside the hg38
# repeat mask (centromeres, gaps, Satellite, segmental duplications) and the
# ENCODE blacklist; sites and genome are both restricted to it. Early / mid /
# late = thirds of the callable genome (by bp) ranked by RT.
#
# Expected counts are context-matched: each site contributes the fraction of
# callable occurrences of its trinucleotide (pyrimidine-centred, from
# BSgenome.Hsapiens.UCSC.hg38) that fall in each RT class, since Pol e
# signature contexts are not uniform over replication timing. Enrichment =
# observed / expected per class, 95% CI from 1000 bootstrap resamples of sites,
# chi-square goodness of fit; mean RT of the sites is reported alongside. Site sets: all tree sites; tree sites chosen as the most
# frequent ("top_all") vs population-enriched; all called sites; Pol e
# signature (TCT>TAT, TTT>TGT on the pyrimidine strand) vs other SBS; and, if
# K562_bulk_sbs_sites.R has been run, de novo SBS from bulk WGS of the
# mutation-accumulation clones (P286R private / persistent, by signature; WT).

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
})

site_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_sites")
rt_bed <- file.path(project_root, "data/reference/replication_timing/K562_replication_timing_hg38.bed.gz")
mask_dir <- file.path(project_root, "results/K562_tree/masks")
output_dir <- file.path(project_root, "paper_figures/output/K562_sbs_rt_enrichment/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(1)
genome <- BSgenome.Hsapiens.UCSC.hg38
chroms <- c(paste0("chr", 1:22), "chrX")

read_bed <- function(f, skip = 0) {
  b <- read.table(f, skip = skip, sep = "\t", quote = "", comment.char = "")
  b <- b[b$V1 %in% chroms, ]
  reduce(GRanges(b$V1, IRanges(b$V2 + 1, b$V3)))
}

# ---- callable genome with RT classes -------------------------------------------------

rt <- read.table(rt_bed, col.names = c("chrom", "start", "end", "G1", "S1", "S2", "S3", "S4", "G2", "RT"))
rt <- rt[rt$chrom %in% chroms, ]
rt_gr <- GRanges(rt$chrom, IRanges(rt$start + 1, rt$end), RT = rt$RT)
callable <- read_bed(file.path(mask_dir, "k100.umap.bed.gz"), skip = 1)
callable <- setdiff(callable, read_bed(file.path(mask_dir, "hg38_repeat_mask.bed")))
callable <- setdiff(callable, read_bed(file.path(project_root, "results/K562_tree/aneufinder_gc/output/blacklist-hg38.raw.bed.gz")))
ov <- findOverlaps(rt_gr, callable)
pieces <- pintersect(rt_gr[queryHits(ov)], callable[subjectHits(ov)])
pieces$RT <- rt_gr$RT[queryHits(ov)]
w <- as.numeric(width(pieces))
o <- order(pieces$RT); cw <- cumsum(w[o]) / sum(w)
cuts <- c(pieces$RT[o][which(cw >= 1 / 3)[1]], pieces$RT[o][which(cw >= 2 / 3)[1]])
classes <- c("late", "mid", "early")
pieces$class <- factor(classes[findInterval(pieces$RT, cuts, left.open = TRUE) + 1], levels = classes)
cat(sprintf("callable genome with RT: %.2f Gb in %d pieces; RT cut-offs %.3f / %.3f\n",
            sum(w) / 1e9, length(pieces), cuts[1], cuts[2]))

# ---- genome trinucleotide counts per class (pyrimidine-centred) ------------------------

tri <- names(trinucleotideFrequency(DNAString("ACGT")))
pyr_of <- function(k) {
  mid <- substr(k, 2, 2)
  ifelse(mid %in% c("C", "T"), k, as.character(reverseComplement(DNAStringSet(k))))
}
tri_pyr <- pyr_of(tri)
tri_cache <- file.path(output_dir, "callable_trinucleotide_counts_by_rt.csv")
genome_tri <- if (file.exists(tri_cache)) as.matrix(read.csv(tri_cache, row.names = 1))[, classes] else sapply(classes, function(cl) {
  p <- pieces[pieces$class == cl]
  counts <- numeric(length(tri))
  for (chunk in split(seq_along(p), ceiling(seq_along(p) / 2e5))) {
    counts <- counts + colSums(trinucleotideFrequency(getSeq(genome, p[chunk])))
  }
  tapply(counts, tri_pyr, sum)
})
p_class_given_ctx <- genome_tri / rowSums(genome_tri)
write.csv(genome_tri, tri_cache)   # reused on later runs (same callable genome and RT classes)

# ---- sites ---------------------------------------------------------------------------

annotate_sites <- function(d) {
  d <- d[d$chrom %in% chroms & nchar(d$ref) == 1 & nchar(d$alt) == 1, ]
  g <- GRanges(d$chrom, IRanges(d$pos, width = 1))
  h <- findOverlaps(g, pieces, select = "first")
  d$rt <- pieces$RT[h]; d$class <- pieces$class[h]
  d <- d[!is.na(h), ]
  g <- GRanges(d$chrom, IRanges(d$pos - 1, width = 3))
  ctx <- as.character(getSeq(genome, g))
  ok <- substr(ctx, 2, 2) == d$ref & !grepl("N", ctx)
  d <- d[ok, ]; ctx <- ctx[ok]
  flip <- d$ref %in% c("A", "G")
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  d$context <- ifelse(flip, as.character(reverseComplement(DNAStringSet(ctx))), ctx)
  d$alt_pyr <- ifelse(flip, comp[d$alt], d$alt)
  d$type <- paste0(d$context, ">", substr(d$context, 1, 1), d$alt_pyr, substr(d$context, 3, 3))
  d$pole_signature <- d$type %in% c("TCT>TAT", "TTT>TGT")
  d
}

tree <- annotate_sites(read.csv(file.path(site_dir, "sbs_tree_sites.csv")))
all_sites <- annotate_sites(read.csv(gzfile(file.path(site_dir, "all_called_sites.csv.gz"))))
cat(sprintf("tree sites in callable genome with RT: %d; all called sites: %d\n", nrow(tree), nrow(all_sites)))

enrich <- function(d, label, n_boot = 1000) {
  if (nrow(d) < 30) return(NULL)
  obs <- table(factor(d$class, levels = classes))
  pc <- p_class_given_ctx[d$context, , drop = FALSE]
  expct <- colSums(pc)
  boot <- replicate(n_boot, { i <- sample(nrow(d), replace = TRUE)
    table(factor(d$class[i], levels = classes)) / colSums(pc[i, , drop = FALSE]) })
  data.frame(set = label, n_sites = nrow(d), class = classes,
             observed = as.vector(obs), expected = as.vector(expct),
             enrichment = as.vector(obs / expct),
             ci_lo = apply(boot, 1, quantile, 0.025), ci_hi = apply(boot, 1, quantile, 0.975),
             chisq_p = suppressWarnings(chisq.test(as.vector(obs), p = expct / sum(expct))$p.value),
             mean_rt = mean(d$rt))
}

bulk_csv <- file.path(project_root, "paper_figures/output/K562_bulk_sbs_sites/bulk_denovo_sites.csv.gz")
bulk <- if (file.exists(bulk_csv)) annotate_sites(read.csv(gzfile(bulk_csv))) else NULL
sets <- list(
  "tree sites (all)" = tree,
  "tree sites: most frequent" = tree[grepl("top_all", tree$selected_for), ],
  "tree sites: population-enriched" = tree[!grepl("top_all", tree$selected_for), ],
  "tree sites: Pol e signature" = tree[tree$pole_signature, ],
  "tree sites: other SBS" = tree[!tree$pole_signature, ],
  "all called sites" = all_sites,
  "all called: Pol e signature" = all_sites[all_sites$pole_signature, ],
  "all called: other SBS" = all_sites[!all_sites$pole_signature, ])
if (!is.null(bulk)) sets <- c(sets, list(
  "bulk P286R private" = bulk[bulk$P286R_private, ],
  "bulk P286R persistent" = bulk[bulk$P286R_persistent, ],
  "bulk P286R persistent: Pol e signature" = bulk[bulk$P286R_persistent & bulk$pole_signature, ],
  "bulk P286R persistent: other SBS" = bulk[bulk$P286R_persistent & !bulk$pole_signature, ],
  "bulk WT private" = bulk[bulk$WT_private, ],
  "bulk WT persistent" = bulk[bulk$WT_persistent, ]))
res <- do.call(rbind, lapply(names(sets), function(k) enrich(sets[[k]], k)))
write.csv(res, file.path(output_dir, "rt_enrichment.csv"), row.names = FALSE)
write.csv(tree, file.path(output_dir, "tree_sites_annotated.csv"), row.names = FALSE)

# expected class shares if sites followed the genome without context matching (1/3 each),
# and tree sites vs all called sites directly
direct <- chisq.test(rbind(table(factor(tree$class, levels = classes)),
                           table(factor(all_sites$class, levels = classes))))
cat(sprintf("tree vs all called sites (class counts): chi-square p = %.3g\n", direct$p.value))
print(res[, c("set", "n_sites", "class", "observed", "expected", "enrichment", "ci_lo", "ci_hi", "chisq_p")],
      digits = 3, row.names = FALSE)

# ---- figure ----------------------------------------------------------------------------

res$class <- factor(res$class, levels = c("early", "mid", "late"))
res$set <- factor(res$set, levels = names(sets))
p <- ggplot(res, aes(class, enrichment)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_col(fill = "grey60", width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.3) +
  facet_wrap(~set, ncol = 4) +
  labs(x = "Replication timing (thirds of callable genome)", y = "Observed / expected SBS (context-matched)") +
  theme_classic(base_size = 7)
h <- 1.8 * ceiling(length(unique(res$set)) / 4)
ggsave(file.path(output_dir, "rt_enrichment.pdf"), p, width = 7, height = h)
ggsave(file.path(output_dir, "rt_enrichment.png"), p, width = 7, height = h, dpi = 200, bg = "white")
cat("Outputs written to", output_dir, "\n")
