# How well do per-cell K562 copy-number calls match the K562 reference
# karyotype, before and after GC correction?
#
# old: sc_PolE_novaseq/AneuFinder_output/result.csv - no GC correction or
#      blacklist; rows are segments disjoined across cells (1-3 Mb). Its
#      columns are relabelled below; see the note where it is read.
# new: results/K562_tree/aneufinder_gc/output/result.csv - GC-corrected,
#      ENCODE blacklist, 1 Mb bins (scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh)
# reference: K562_ref_ploidy_hg38.bed, the reference profile used in
#      notebooks/K562_tree.ipynb cells 135-148
#
# Both call sets are put on the new 1 Mb grid. Each bin takes the reference
# segment covering its midpoint; bins the reference does not cover are dropped
# rather than filled with 3 as in the notebook, so they cannot count as
# agreement. chrY is dropped. Per cell: Pearson r with the reference (blind to
# overall scale, so it also works for mis-scaled cells) and the fraction of bins
# at exactly the reference copy number. Each is computed on all bins and again
# with the most GC-rich bins removed (top gc_drop_frac by hg38 GC content).

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
new_csv <- file.path(project_root, "results/K562_tree/aneufinder_gc/output/result.csv")
output_dir <- file.path(project_root, "paper_figures/output/K562_cnv_vs_reference/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gc_drop_frac <- 0.10

# ---- common 1 Mb grid --------------------------------------------------------

new <- read.csv(new_csv, row.names = 1, check.names = FALSE)
old <- read.csv(file.path(sc_dir, "AneuFinder_output/result.csv"),
                row.names = 1, check.names = FALSE)
# The old result.csv is mislabelled: sc_PolE_novaseq/bin/aneufinder.R loaded the
# models in alphabetical file order but named the columns from adata.obs_names,
# which is in a different order. Column i holds the i-th cell alphabetically
# (checked against the per-cell MODELS files: 100 of 100 columns match exactly).
cell_cols <- setdiff(colnames(old), c("seqnames", "start", "end"))
colnames(old)[match(cell_cols, colnames(old))] <- sort(cell_cols)
meta <- read.csv(file.path(sc_dir, "results/full_meta.csv"))
cells <- meta$X
stopifnot(all(cells %in% colnames(new)), all(cells %in% colnames(old)))

grid <- GRanges(new$seqnames, IRanges(new$start, new$end))
mid <- resize(grid, width = 1, fix = "center")

# old segments -> grid bins by midpoint
old_gr <- GRanges(old$seqnames, IRanges(old$start, old$end))
hit_old <- findOverlaps(mid, old_gr, select = "first")

ref <- read.table(file.path(sc_dir, "K562_ref_ploidy_hg38.bed"),
                  col.names = c("chr", "start", "end", "cn"))
# bed is 0-based half-open
ref_gr <- GRanges(ref$chr, IRanges(ref$start + 1, ref$end))
hit_ref <- findOverlaps(mid, ref_gr, select = "first")

gc <- letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, grid), "GC", as.prob = TRUE)[, 1] /
  (1 - letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, grid), "N", as.prob = TRUE)[, 1])

keep <- !is.na(hit_ref) & !is.na(hit_old) & new$seqnames != "chrY" & is.finite(gc)
gc_cut <- quantile(gc[keep], 1 - gc_drop_frac)
low_gc <- keep & gc < gc_cut
cat(sprintf("%d grid bins; %d covered by the reference and both call sets (chrY dropped)\n",
            length(grid), sum(keep)))
cat(sprintf("GC filter: drop bins with GC >= %.3f, leaving %d bins\n", gc_cut, sum(low_gc)))

ref_cn <- ref$cn[hit_ref]
cn <- list(
  old = as.matrix(old[hit_old, cells]),
  new = as.matrix(new[, cells])
)

# ---- per-cell agreement ------------------------------------------------------

score <- function(m, bins) {
  x <- m[bins, , drop = FALSE]
  r <- ref_cn[bins]
  data.frame(
    cell = colnames(x),
    pearson = apply(x, 2, function(v) if (sd(v) > 0) cor(v, r) else NA),
    frac_exact = colMeans(x == r),
    ploidy = colMeans(x)
  )
}

per_cell <- bind_rows(lapply(names(cn), function(run) {
  bind_rows(
    mutate(score(cn[[run]], keep), bins = "All bins"),
    mutate(score(cn[[run]], low_gc), bins = "High-GC bins removed")
  ) %>% mutate(run = run)
})) %>%
  left_join(select(meta, cell = X, pop), by = "cell")

# dominant ploidy per run, as in K562_cnv_event_tree.R
dominant <- per_cell %>%
  filter(bins == "All bins") %>%
  group_by(run) %>%
  reframe(cell = cell[ploidy >= 2.5 & ploidy < 3.0])
in_both <- dominant %>% count(cell) %>% filter(n == 2) %>% pull(cell)

summary_tbl <- bind_rows(
  per_cell %>% mutate(cells = "All 1000"),
  per_cell %>% filter(cell %in% in_both) %>% mutate(cells = "Dominant ploidy in both runs")
) %>%
  group_by(cells, bins, run) %>%
  summarise(n_cells = n(),
            pearson_median = median(pearson, na.rm = TRUE),
            pearson_q25 = quantile(pearson, 0.25, na.rm = TRUE),
            pearson_q75 = quantile(pearson, 0.75, na.rm = TRUE),
            frac_exact_median = median(frac_exact),
            .groups = "drop")

# paired test: does GC correction raise each cell's agreement?
paired <- per_cell %>%
  select(cell, bins, run, pearson) %>%
  pivot_wider(names_from = run, values_from = pearson) %>%
  group_by(bins) %>%
  summarise(n = sum(!is.na(old) & !is.na(new)),
            frac_cells_improved = mean(new > old, na.rm = TRUE),
            median_change = median(new - old, na.rm = TRUE),
            wilcoxon_p = wilcox.test(new, old, paired = TRUE)$p.value,
            .groups = "drop")

write.csv(per_cell, file.path(output_dir, "per_cell_agreement.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(output_dir, "agreement_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(output_dir, "paired_change.csv"), row.names = FALSE)
print(as.data.frame(summary_tbl), digits = 3)
print(as.data.frame(paired), digits = 3)

# ---- plot ----------------------------------------------------------------------

# each cell's reference correlation before (x) and after (y) GC correction
per_cell %>%
  select(cell, bins, run, pearson) %>%
  pivot_wider(names_from = run, values_from = pearson) %>%
  ggplot(aes(old, new)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~bins) +
  coord_equal(xlim = c(-0.2, 1), ylim = c(-0.2, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  xlab("Pearson r with K562 reference, uncorrected") +
  ylab("Pearson r, GC-corrected")
ggsave(file.path(output_dir, "reference_correlation_old_vs_new.svg"),
       bg = "transparent", height = 3.25, width = 5.5)

cat("\nWrote results to:", output_dir, "\n")
