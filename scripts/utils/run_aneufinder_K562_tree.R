#!/usr/bin/env Rscript

# Run AneuFinder jointly on K562_tree single-cell BAMs.
# Inputs aligned to GATK hg38 (chr-prefixed). Uses 1 Mb bins, edivisive
# segmentation, GC correction via BSgenome.Hsapiens.UCSC.hg38, and the
# ENCODE hg38 blacklist (downloaded on first run).

suppressPackageStartupMessages({
    library(AneuFinder)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

repo_root    <- "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome"
inputfolder  <- file.path(repo_root, "data/K562_tree/sc_outputs")
outputfolder <- file.path(repo_root, "results/K562_tree/aneufinder")
dir.create(outputfolder, showWarnings = FALSE, recursive = TRUE)

ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
cat(sprintf("Using %d CPUs\n", ncpu))

blacklist_url     <- "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
blacklist_raw     <- file.path(outputfolder, "blacklist-hg38.raw.bed.gz")
blacklist_path    <- file.path(outputfolder, "blacklist-hg38.bed.gz")
if (!file.exists(blacklist_raw)) {
    cat("Downloading ENCODE hg38 blacklist...\n")
    download.file(blacklist_url, blacklist_raw, mode = "wb", quiet = TRUE)
}
# AneuFinder's importBed() requires 6 columns (chr,start,end,name,score,strand);
# ENCODE blacklist ships as 3 columns, so pad with placeholder values.
if (!file.exists(blacklist_path)) {
    bl <- read.table(blacklist_raw, sep = "\t", stringsAsFactors = FALSE,
                     colClasses = c("character", "integer", "integer"))
    bl$name   <- "."
    bl$score  <- 0L
    bl$strand <- "."
    gz <- gzfile(blacklist_path, "w")
    write.table(bl, gz, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    close(gz)
}

chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

cat(sprintf("Starting AneuFinder on %d BAMs\n",
            length(list.files(inputfolder, pattern = "\\.bam$"))))
t0 <- Sys.time()

Aneufinder(
    inputfolder       = inputfolder,
    outputfolder      = outputfolder,
    assembly          = "hg38",
    numCPU            = ncpu,
    binsizes          = 1e6,
    chromosomes       = chromosomes,
    blacklist         = blacklist_path,
    correction.method = "GC",
    GC.BSgenome       = BSgenome.Hsapiens.UCSC.hg38,
    method            = "edivisive",
    refine.breakpoints = FALSE
)

cat(sprintf("Aneufinder finished in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

models_dir <- file.path(outputfolder, "MODELS", "method-edivisive")
model_files <- list.files(models_dir, pattern = "\\.RData$", full.names = TRUE)
cat(sprintf("Loading %d model files for heatmap\n", length(model_files)))

models <- loadFromFiles(model_files)
cl     <- clusterByQuality(models)

heatmap_pdf <- file.path(outputfolder, "heatmapGenomewide.pdf")
pdf(heatmap_pdf, width = 16, height = max(8, length(models) * 0.05))
heatmapGenomewide(cl$classification[[1]])
dev.off()
cat(sprintf("Wrote %s\n", heatmap_pdf))
