#!/usr/bin/env Rscript
###############################################################################
# Estimate ambient RNA contamination with SoupX and write decontaminated counts.
#
# INPUT LAYOUT (under <input_dir>):
#   counts.mtx        genes x cells sparse (real cells, from adata.layers['counts'].T)
#   genes.tsv         one gene ID per line, matching counts row order
#   barcodes.tsv      one cell barcode per line, matching counts column order
#   clusters.csv      columns: barcode,cluster  (Leiden ids as strings)
#   empties/          OPTIONAL subdirectory — if present, contents are used as
#                     the soup source (tod). Same three-file layout inside:
#     empties/counts.mtx     genes x empties
#     empties/genes.tsv      MUST match ../genes.tsv order (aligned in Python)
#     empties/barcodes.tsv   one empty-droplet barcode per line
#
# If empties/ is absent, we fall back to using the cell aggregate as the soup
# profile (SoupChannel with tod == toc), which is an approximation and can
# fail autoEstCont on sparse data.
#
# USAGE:
#   Rscript scripts/CapGTA/soupx_decontaminate.R <input_dir> <output_dir> [manual_rho]
#     manual_rho: optional float in (0,1). If set, skip autoEstCont and use
#                 this contamination fraction for every cell.
#
# OUTPUT LAYOUT (under <output_dir>):
#   counts_decontam.mtx   integer-rounded corrected counts (genes x cells)
#   genes.tsv             gene IDs (copied — same order)
#   barcodes.tsv          cell barcodes (copied — same order)
#   soupx_summary.csv     per-cell: cluster, rho, nUMIs, nUMIs_after
#   soupx_soup_profile.csv per-gene soup fraction (est. ambient composition)
###############################################################################

suppressPackageStartupMessages({
  library(SoupX)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: soupx_decontaminate.R <input_dir> <output_dir> [manual_rho]")
}
input_dir  <- args[1]
output_dir <- args[2]
manual_rho <- if (length(args) >= 3) as.numeric(args[3]) else NA_real_

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading counts…")
counts <- readMM(file.path(input_dir, "counts.mtx"))
rownames(counts) <- readLines(file.path(input_dir, "genes.tsv"))
colnames(counts) <- readLines(file.path(input_dir, "barcodes.tsv"))
counts <- as(counts, "CsparseMatrix")

clusters_df <- read.csv(file.path(input_dir, "clusters.csv"),
                        stringsAsFactors = FALSE)
stopifnot(c("barcode", "cluster") %in% colnames(clusters_df))
rownames(clusters_df) <- clusters_df$barcode
stopifnot(all(colnames(counts) %in% rownames(clusters_df)))
clusters <- as.character(clusters_df[colnames(counts), "cluster"])
names(clusters) <- colnames(counts)

message(sprintf("Genes: %d   Cells: %d   Clusters: %d",
                nrow(counts), ncol(counts), length(unique(clusters))))

empties_dir <- file.path(input_dir, "empties")
if (dir.exists(empties_dir)) {
  message("Loading empty droplet counts (soup source)…")
  tod <- readMM(file.path(empties_dir, "counts.mtx"))
  rownames(tod) <- readLines(file.path(empties_dir, "genes.tsv"))
  colnames(tod) <- readLines(file.path(empties_dir, "barcodes.tsv"))
  tod <- as(tod, "CsparseMatrix")
  stopifnot(identical(rownames(tod), rownames(counts)))
  message(sprintf("  Empties: %d droplets, total %s counts",
                  ncol(tod), format(sum(tod), big.mark=",")))
  sc <- SoupChannel(tod = tod, toc = counts, calcSoupProfile = TRUE)
} else {
  message("No empties/ subdir found — falling back to tod = toc (approximate)")
  sc <- SoupChannel(tod = counts, toc = counts, calcSoupProfile = TRUE)
}
sc <- setClusters(sc, clusters)

if (!is.na(manual_rho)) {
  message(sprintf("Applying manual contamination fraction rho = %.3f", manual_rho))
  sc <- setContaminationFraction(sc, manual_rho, forceAccept = TRUE)
} else {
  message("Estimating contamination fraction via autoEstCont…")
  sc <- autoEstCont(sc, forceAccept = TRUE, doPlot = FALSE)
}
message(sprintf("Median rho: %.3f   Mean rho: %.3f",
                median(sc$metaData$rho), mean(sc$metaData$rho)))

message("Adjusting counts (roundToInt = TRUE)…")
adj <- adjustCounts(sc, roundToInt = TRUE)

message("Writing outputs…")
writeMM(as(adj, "CsparseMatrix"),
        file.path(output_dir, "counts_decontam.mtx"))
writeLines(rownames(adj), file.path(output_dir, "genes.tsv"))
writeLines(colnames(adj), file.path(output_dir, "barcodes.tsv"))

nUMIs_before <- colSums(counts)
nUMIs_after  <- colSums(adj)
summary_df <- data.frame(
  barcode      = rownames(sc$metaData),
  cluster      = sc$metaData$clusters,
  rho          = sc$metaData$rho,
  nUMIs        = nUMIs_before[rownames(sc$metaData)],
  nUMIs_after  = nUMIs_after[rownames(sc$metaData)]
)
write.csv(summary_df,
          file.path(output_dir, "soupx_summary.csv"),
          row.names = FALSE)

soup_df <- data.frame(
  gene_id  = rownames(sc$soupProfile),
  est_frac = sc$soupProfile$est,
  counts   = sc$soupProfile$counts
)
write.csv(soup_df[order(-soup_df$est_frac), ],
          file.path(output_dir, "soupx_soup_profile.csv"),
          row.names = FALSE)

message(sprintf("Done. Wrote counts_decontam.mtx (%d x %d) and summaries to %s",
                nrow(adj), ncol(adj), output_dir))
