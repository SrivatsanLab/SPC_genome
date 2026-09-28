# hg38 mask of large repetitive regions for the per-cell CNV event tree.
#
# Merges four UCSC hg38 tables, downloaded into results/K562_tree/masks/ from
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ :
#   centromeres.txt.gz       centromere models
#   gap.txt.gz               assembly gaps (telomeres, heterochromatin, etc.)
#   genomicSuperDups.txt.gz  segmental duplications
#   rmsk.txt.gz              RepeatMasker, Satellite class only
# Interspersed repeats (LINE/SINE/etc.) are left in: they cover about half of
# every 1 Mb bin, so masking them would remove the whole genome.
#
# Writes hg38_repeat_mask.bed (merged, 0-based) and the fraction of each 1 Mb
# AneuFinder bin it covers, for K562_cnv_event_tree.R (K562_MASK_BED).

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
})

mask_dir <- file.path(project_root, "results/K562_tree/masks")
bins_csv <- file.path(project_root, "results/K562_tree/aneufinder_gc/output/result.csv")

ucsc <- function(file, chrom, start, end, filter = NULL) {
  d <- fread(file.path(mask_dir, file), header = FALSE, sep = "\t", quote = "")
  if (!is.null(filter)) d <- d[filter(d)]
  # UCSC tables are 0-based half-open
  GRanges(d[[chrom]], IRanges(d[[start]] + 1, d[[end]]), source = file)
}
parts <- list(
  centromere = ucsc("centromeres.txt.gz", 2, 3, 4),
  gap = ucsc("gap.txt.gz", 2, 3, 4),
  segdup = ucsc("genomicSuperDups.txt.gz", 2, 3, 4),
  satellite = ucsc("rmsk.txt.gz", 6, 7, 8, filter = function(d) d[[12]] == "Satellite")
)
for (nm in names(parts)) {
  cat(sprintf("%-10s %7d intervals, %6.1f Mb\n", nm, length(parts[[nm]]),
              sum(width(reduce(parts[[nm]]))) / 1e6))
}
# the tables span different unplaced/alt contigs; the seqinfo mismatch is harmless
mask <- reduce(suppressWarnings(do.call(c, unname(lapply(parts, granges)))))
mask <- mask[seqnames(mask) %in% c(paste0("chr", 1:22), "chrX", "chrY")]
cat(sprintf("merged mask: %.1f Mb\n", sum(width(mask)) / 1e6))

write.table(data.frame(seqnames(mask), start(mask) - 1, end(mask)),
            file.path(mask_dir, "hg38_repeat_mask.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# fraction of each AneuFinder bin covered, overall and by source
bins <- fread(bins_csv, select = c("seqnames", "start", "end"))
bins_gr <- GRanges(bins$seqnames, IRanges(bins$start, bins$end))
covered <- function(gr) {
  hits <- findOverlaps(bins_gr, reduce(gr))
  ov <- pintersect(bins_gr[queryHits(hits)], reduce(gr)[subjectHits(hits)])
  out <- numeric(length(bins_gr))
  agg <- tapply(width(ov), queryHits(hits), sum)
  out[as.integer(names(agg))] <- agg
  out / width(bins_gr)
}
frac <- data.frame(bins, mask_frac = covered(mask),
                   sapply(parts, covered))
fwrite(frac, file.path(mask_dir, "bin_mask_fraction.csv"))
for (t in c(0.1, 0.25, 0.5)) {
  cat(sprintf("bins with >= %.0f%% masked: %d of %d\n", 100 * t,
              sum(frac$mask_frac >= t), nrow(frac)))
}
cat("Wrote", file.path(mask_dir, "hg38_repeat_mask.bed"), "\n")
