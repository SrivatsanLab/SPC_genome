#!/usr/bin/env Rscript
# Per-sample AneuFinder CNV calling for the benchmarking pipeline.
#
# Inputs : a MarkDup'd BAM, sample id, output dir, paired-end flag.
# Outputs: <out>/<sample>.aneufinder.rds  (HMM + edivisive model objects)
#          <out>/<sample>.segments.tsv     (flattened segments, both methods)
#          <out>/<sample>.bins.tsv.gz      (per-bin copy state, HMM)
#          <out>/<sample>.qc.tsv           (one-row spikiness / aneuploidy / etc.)
#
# Run via:
#   ml fhR
#   Rscript aneufinder_one_sample.R --bam <bam> --sample <id> --out-dir <dir> \
#       [--paired-end TRUE|FALSE] [--binsize 1000000] [--stepsize 500000]

suppressPackageStartupMessages({
  library(AneuFinder)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(optparse)
  library(GenomicRanges)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bam",         type = "character", help = "Input BAM"),
  make_option("--sample",      type = "character", help = "Sample ID"),
  make_option("--out-dir",     type = "character", help = "Output directory"),
  make_option("--paired-end",  type = "logical",   default = TRUE,
              help = "TRUE for paired-end BAMs, FALSE for single-end"),
  make_option("--binsize",     type = "integer",   default = 1000000L),
  make_option("--stepsize",    type = "integer",   default = 1000000L,
              help = "Default = binsize (non-overlapping bins, what HMM uses)."),
  make_option("--chromosomes", type = "character", default = "autosomes+X",
              help = "'autosomes', 'autosomes+X', or comma-separated chr names")
)))

stopifnot(!is.null(opt$bam), !is.null(opt$sample), !is.null(opt$`out-dir`))
if (!file.exists(opt$bam)) stop("BAM not found: ", opt$bam)
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

chroms <- switch(opt$chromosomes,
  "autosomes"   = paste0("chr", 1:22),
  "autosomes+X" = c(paste0("chr", 1:22), "chrX"),
  strsplit(opt$chromosomes, ",")[[1]]
)

cat("[", format(Sys.time()), "] binReads (binsize=", opt$binsize,
    ", paired=", opt$`paired-end`, ")\n", sep = "")
binned_list <- binReads(
  file           = opt$bam,
  assembly       = "hg38",
  binsizes       = opt$binsize,
  stepsizes      = opt$stepsize,
  chromosomes    = chroms,
  pairedEndReads = opt$`paired-end`
)
# binReads returns a list keyed by paste0("binsize_", bs, "_stepsize_", ss). One element here.
binned <- binned_list[[1]]

cat("[", format(Sys.time()), "] correctGC\n", sep = "")
binned_list <- correctGC(list(binned), GC.BSgenome = BSgenome.Hsapiens.UCSC.hg38)
binned <- binned_list[[1]]

cat("[", format(Sys.time()), "] findCNVs (HMM)\n", sep = "")
hmm <- findCNVs(binned, method = "HMM", ID = opt$sample)

cat("[", format(Sys.time()), "] findCNVs (edivisive)\n", sep = "")
edi <- findCNVs(binned, method = "edivisive", ID = opt$sample)

# --- Flatten segments to TSV --------------------------------------------------
segs_to_df <- function(model, label) {
  segs <- model$segments
  if (is.null(segs) || length(segs) == 0L) return(data.frame())
  df <- as.data.frame(segs)
  df$method <- label
  df$sample <- opt$sample
  df
}
align_cols <- function(dfs) {
  dfs <- dfs[vapply(dfs, function(d) nrow(d) > 0L, logical(1))]
  if (length(dfs) == 0L) return(data.frame())
  cols <- unique(unlist(lapply(dfs, colnames)))
  for (i in seq_along(dfs)) {
    miss <- setdiff(cols, colnames(dfs[[i]]))
    for (m in miss) dfs[[i]][[m]] <- NA
    dfs[[i]] <- dfs[[i]][, cols, drop = FALSE]
  }
  do.call(rbind, dfs)
}
hmm_df <- segs_to_df(hmm, "HMM")
edi_df <- segs_to_df(edi, "edivisive")
seg_df <- align_cols(list(hmm_df, edi_df))
out_seg <- file.path(opt$`out-dir`, paste0(opt$sample, ".segments.tsv"))
write.table(seg_df, out_seg, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote ", nrow(seg_df), " segments to ", out_seg, "\n", sep = "")

# --- Per-bin TSV (counts + GC from `binned`, state + copy.number from `hmm`) --
bins_df <- as.data.frame(binned)
state_df <- as.data.frame(hmm$bins)
if (nrow(bins_df) == nrow(state_df)) {
  bins_df$state <- state_df$state
  bins_df$copy.number <- state_df$copy.number
} else {
  warning(sprintf("Row mismatch binned (%d) vs hmm$bins (%d); skipping state merge.",
                  nrow(bins_df), nrow(state_df)))
}
bins_df$sample <- opt$sample
out_bins <- file.path(opt$`out-dir`, paste0(opt$sample, ".bins.tsv.gz"))
gz <- gzfile(out_bins, "w")
write.table(bins_df, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
cat("Wrote ", nrow(bins_df), " bins to ", out_bins, "\n", sep = "")

# --- One-row QC TSV -----------------------------------------------------------
qi <- if (!is.null(hmm$qualityInfo)) hmm$qualityInfo else list()
get1 <- function(name) {
  v <- qi[[name]]
  if (is.null(v) || length(v) == 0L) return(NA_real_)
  x <- suppressWarnings(as.numeric(v[[1]]))
  if (length(x) == 0L || is.na(x[1])) NA_real_ else x[1]
}
# Force counts to numeric (AneuFinder may store as integer or list)
bin_counts <- suppressWarnings(as.numeric(bins_df$counts))
bin_widths <- suppressWarnings(as.numeric(bins_df$width))

total_bp <- sum(bin_widths, na.rm = TRUE)
bp_nondip <- if (nrow(hmm_df) > 0 && "copy.number" %in% names(hmm_df))
  sum(as.numeric(hmm_df$width[hmm_df$copy.number != 2]), na.rm = TRUE) else 0
n_seg_nondip <- if (nrow(hmm_df) > 0 && "copy.number" %in% names(hmm_df))
  sum(hmm_df$copy.number != 2, na.rm = TRUE) else 0

qc_df <- data.frame(
  sample                   = opt$sample,
  bam                      = opt$bam,
  paired_end               = opt$`paired-end`,
  binsize                  = opt$binsize,
  n_bins                   = nrow(bins_df),
  median_bin_count         = median(bin_counts, na.rm = TRUE),
  mean_bin_count           = mean(bin_counts, na.rm = TRUE),
  spikiness                = get1("spikiness"),
  entropy                  = get1("entropy"),
  bhattacharyya            = get1("bhattacharyya"),
  num.segments             = get1("num.segments"),
  loglik                   = get1("loglik"),
  n_segments_hmm           = nrow(hmm_df),
  n_segments_hmm_nondiploid= n_seg_nondip,
  n_segments_edi           = nrow(edi_df),
  total_bp                 = total_bp,
  bp_nondiploid_hmm        = bp_nondip,
  aneuploidy_index_hmm     = if (total_bp > 0) bp_nondip / total_bp else NA_real_,
  stringsAsFactors         = FALSE
)
out_qc <- file.path(opt$`out-dir`, paste0(opt$sample, ".qc.tsv"))
write.table(qc_df, out_qc, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote QC: ", out_qc, "\n", sep = "")

# --- Save model objects (RDS) -------------------------------------------------
out_rds <- file.path(opt$`out-dir`, paste0(opt$sample, ".aneufinder.rds"))
saveRDS(list(
  sample    = opt$sample,
  binsize   = opt$binsize,
  stepsize  = opt$stepsize,
  paired    = opt$`paired-end`,
  hmm       = hmm,
  edivisive = edi
), out_rds)
cat("Wrote RDS: ", out_rds, "\n", sep = "")
cat("[", format(Sys.time()), "] done: ", opt$sample, "\n", sep = "")
