# How many of the single-cell SBS behind the per-cell SBS tree are germline
# (or ancestral K562) variants according to bulk calls?
#
# Sites: K562_sbs_tree_sites.R (the 700 tree sites, and all 3.1 M sites in
# sc_PolE_filtered_no_singletons.h5ad), 1-based.
# Germline references:
#   bulk     - K562_mut_accumulation_raw.h5ad (17 bulk WGS samples of the
#              mutation-accumulation clones): pooled VAF = sum AD / sum DP over
#              samples; germline/ancestral if pooled VAF >= 0.15 (the cut the
#              mutation-accumulation analysis uses to set them aside), low-VAF
#              (de novo-like) if < 0.15, "not in bulk" if the site/allele is
#              absent from the bulk call set
#   ENCODE   - ENCFF538YDL (K562 germline variants, hg38), matched on
#              chromosome, position and alternate allele
# Writes paper_figures/output/K562_sbs_tree_germline/: per-site tables for the
# tree sites and for all called sites, and summary counts for tree sites (by selection group) and for all
# called sites.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
})

site_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_sites")
bulk_h5 <- Sys.getenv("K562_MUTACC_H5AD",
                      "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/results/K562_mut_accumulation/K562_mut_accumulation_raw.h5ad")
encode_vcf <- Sys.getenv("K562_ENCODE_VCF",
                         "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/ENCFF538YDL_hg38.vcf.gz")
output_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_germline/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tree <- read.csv(file.path(site_dir, "sbs_tree_sites.csv"))
all_sites <- read.csv(gzfile(file.path(site_dir, "all_called_sites.csv.gz")))
key <- function(chrom, pos, alt) paste(chrom, pos, alt, sep = ":")
tree$key <- key(tree$chrom, tree$pos, tree$alt)
all_sites$key <- key(all_sites$chrom, all_sites$pos, all_sites$alt)

# ---- bulk pooled VAF (positions in this h5ad are 1-based) -------------------------------
bsites <- as.character(h5read(bulk_h5, "/var/_index"))
bp <- do.call(rbind, strsplit(bsites, "[->]"))
bkey <- key(bp[, 1], as.integer(bp[, 2]), bp[, 4])
read_csr <- function(group) {
  shape <- h5readAttributes(bulk_h5, "/X")$shape
  m <- sparseMatrix(j = as.integer(h5read(bulk_h5, paste0(group, "/indices"))) + 1L,
                    p = as.integer(h5read(bulk_h5, paste0(group, "/indptr"))),
                    x = as.numeric(h5read(bulk_h5, paste0(group, "/data"))),
                    dims = shape, index1 = TRUE, repr = "R")
  as(m, "CsparseMatrix")
}
AD <- read_csr("/layers/AD"); DP <- read_csr("/layers/DP"); h5closeAll()
need <- which(bkey %in% all_sites$key)
pooled <- colSums(AD[, need]) / pmax(colSums(DP[, need]), 1)
depth <- colSums(DP[, need])
bulk_vaf <- setNames(pooled, bkey[need]); bulk_depth <- setNames(depth, bkey[need])
cat(sprintf("bulk: %d sites; %d of the single-cell sites present\n", length(bkey), length(need)))

# ---- ENCODE K562 germline ----------------------------------------------------------------
enc <- read.table(gzfile(encode_vcf), comment.char = "#", sep = "\t", quote = "",
                  colClasses = c("character", "integer", "NULL", "character", "character", rep("NULL", 5)),
                  fill = TRUE)
enc <- enc[nchar(enc$V4) == 1, ]
enc_alts <- strsplit(enc$V5, ",")
enc_key <- key(rep(enc$V1, lengths(enc_alts)), rep(enc$V2, lengths(enc_alts)), unlist(enc_alts))
cat(sprintf("ENCODE K562: %d SNV alleles\n", length(enc_key)))

classify <- function(d) {
  v <- bulk_vaf[d$key]
  d$bulk_vaf <- unname(v); d$bulk_depth <- unname(bulk_depth[d$key])
  d$bulk_class <- ifelse(is.na(v), "not in bulk", ifelse(v >= 0.15, "germline/ancestral (bulk VAF >= 0.15)", "low VAF in bulk (< 0.15)"))
  d$encode_germline <- d$key %in% enc_key
  d$germline_any <- d$bulk_class == "germline/ancestral (bulk VAF >= 0.15)" | d$encode_germline
  d$location <- ifelse(d$chrom == "chrM", "chrM", ifelse(grepl("_", d$chrom), "unplaced/alt contig", "chr1-22/X/Y"))
  d
}
tree <- classify(tree)
all_sites <- classify(all_sites)
tree$group <- ifelse(grepl("top_all", tree$selected_for), "most frequent", "population-enriched")
write.csv(tree, file.path(output_dir, "tree_sites_germline.csv"), row.names = FALSE)
# every called site, for re-selecting tree sites without the flagged ones
write.csv(all_sites[, c("site", "chrom", "pos", "ref", "alt", "n_cells", "bulk_vaf", "bulk_class",
                        "encode_germline", "germline_any", "location")],
          gzfile(file.path(output_dir, "all_sites_germline.csv.gz")), row.names = FALSE)

summ <- function(d, label) {
  data.frame(set = label, n = nrow(d),
             germline_bulk = sum(d$bulk_class == "germline/ancestral (bulk VAF >= 0.15)"),
             germline_encode = sum(d$encode_germline),
             germline_either = sum(d$germline_any),
             low_vaf_bulk = sum(d$bulk_class == "low VAF in bulk (< 0.15)"),
             not_in_bulk = sum(d$bulk_class == "not in bulk"))
}
res <- rbind(summ(tree, "tree sites (700)"),
             summ(tree[tree$group == "most frequent", ], "tree: most frequent"),
             summ(tree[tree$group == "population-enriched", ], "tree: population-enriched"),
             summ(tree[tree$location == "chr1-22/X/Y", ], "tree: nuclear main chromosomes"),
             summ(all_sites, "all called sites"))
res$germline_either_pct <- round(100 * res$germline_either / res$n, 1)
write.csv(res, file.path(output_dir, "germline_summary.csv"), row.names = FALSE)
print(res, row.names = FALSE)
cat("\ntree sites by location and bulk class:\n")
print(table(tree$location, tree$bulk_class))
cat("\ntree sites flagged germline by ENCODE, by bulk class:\n")
print(table(tree$encode_germline, tree$bulk_class))
