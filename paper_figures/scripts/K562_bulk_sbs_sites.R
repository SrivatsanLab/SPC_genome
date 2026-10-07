# De novo SBS called in bulk WGS of the K562 mutation-accumulation clones, for
# testing against replication timing (K562_sbs_rt_enrichment.R).
#
# Input: K562_mut_accumulation_raw.h5ad (17 bulk samples: AAVS (WT) and PolE
# P286R clones 4-6 at passages 1-3; layers AD and DP). Definition follows
# notebooks/K562_mut_accumulation.ipynb (cells 27-36):
#   - single-base substitutions only
#   - pooled VAF across samples (sum AD / sum DP) < 0.15, excluding ancestral
#     and germline K562 variants
#   - "high" in a sample: 10 <= DP <= 100, AD >= 2, VAF > 0.15
# The notebook notes its counts were broken after a change and gives millions
# of sites per clone, so the sets are made stricter here:
#   P286R_private    - high in >= 1 PolE sample and in no AAVS sample
#   P286R_persistent - P286R_private and high in >= 2 passages of the same
#                      PolE clone (a mutation carried forward)
#   WT_private       - high in >= 1 AAVS sample and in no PolE sample
#   WT_persistent    - the same persistence rule within AAVS clones
# Positions are checked against the reference (the base-offset that matches is
# used) and written 1-based.
# Writes paper_figures/output/K562_bulk_sbs_sites/bulk_denovo_sites.csv.gz.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

h5 <- Sys.getenv("K562_MUTACC_H5AD",
                 "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/results/K562_mut_accumulation/K562_mut_accumulation_raw.h5ad")
output_dir <- file.path(project_root, "paper_figures/output/K562_bulk_sbs_sites/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_csr <- function(group) {
  shape <- h5readAttributes(h5, "/X")$shape
  m <- sparseMatrix(j = as.integer(h5read(h5, paste0(group, "/indices"))) + 1L,
                    p = as.integer(h5read(h5, paste0(group, "/indptr"))),
                    x = as.numeric(h5read(h5, paste0(group, "/data"))),
                    dims = shape, index1 = TRUE, repr = "R")
  as(m, "CsparseMatrix")
}
samples <- as.character(h5read(h5, "/obs/_index"))
sites <- as.character(h5read(h5, "/var/_index"))
AD <- read_csr("/layers/AD"); DP <- read_csr("/layers/DP")
h5closeAll()
dimnames(AD) <- dimnames(DP) <- list(samples, NULL)
cat(sprintf("%d samples x %d sites\n", nrow(AD), ncol(AD)))

parts <- do.call(rbind, strsplit(sites, "[->]"))
chrom <- parts[, 1]; pos <- as.integer(parts[, 2]); ref <- parts[, 3]; alt <- parts[, 4]
snv <- nchar(ref) == 1 & nchar(alt) == 1 & alt %in% c("A", "C", "G", "T") &
  chrom %in% c(paste0("chr", 1:22), "chrX", "chrY")
cat(sprintf("%d single-base substitutions on chr1-22/X/Y\n", sum(snv)))

# position base: which offset puts the reference base on the genome?
chk <- sample(which(snv), 1e5)
offset <- c(-1, 0, 1)[which.max(sapply(c(-1, 0, 1), function(o)
  mean(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(chrom[chk], IRanges(pos[chk] + o, width = 1)))) == ref[chk])))]
cat(sprintf("position offset to 1-based: %+d\n", offset))

keep <- which(snv)
ad <- AD[, keep]; dp <- DP[, keep]
pooled <- colSums(ad) / pmax(colSums(dp), 1)
ok <- pooled < 0.15
keep <- keep[ok]; ad <- ad[, ok]; dp <- dp[, ok]
vaf <- ad / pmax(dp, 1)
high <- (dp >= 10) & (dp <= 100) & (ad >= 2) & (vaf > 0.15)
high <- as(high, "lMatrix")
cat(sprintf("%d SNVs with pooled VAF < 0.15\n", length(keep)))

pole <- grepl("^PolE", samples); aavs <- grepl("^AAVS", samples)
clone <- sub("_P[0-9]+$", "", samples)
n_pole <- colSums(high[pole, , drop = FALSE]); n_aavs <- colSums(high[aavs, , drop = FALSE])
persist <- function(sel) {
  per_clone <- sapply(unique(clone[sel]), function(cl) colSums(high[clone == cl, , drop = FALSE]))
  apply(per_clone, 1, max) >= 2
}
out <- data.frame(site = sites[keep], chrom = chrom[keep], pos = pos[keep] + offset,
                  ref = ref[keep], alt = alt[keep], pooled_vaf = pooled[ok],
                  n_P286R_samples = n_pole, n_WT_samples = n_aavs)
out$P286R_private <- n_pole >= 1 & n_aavs == 0
out$P286R_persistent <- out$P286R_private & persist(pole)
out$WT_private <- n_aavs >= 1 & n_pole == 0
out$WT_persistent <- out$WT_private & persist(aavs)
out <- out[out$P286R_private | out$WT_private, ]
write.csv(out, gzfile(file.path(output_dir, "bulk_denovo_sites.csv.gz")), row.names = FALSE)
cat(sprintf("P286R private %d (persistent %d); WT private %d (persistent %d)\n",
            sum(out$P286R_private), sum(out$P286R_persistent), sum(out$WT_private), sum(out$WT_persistent)))
