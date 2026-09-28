# Apparent genome doublings in the K562 single cells are a GC-bias artefact.
#
# Input: scripts/utils/aneufinder_gc_toggle.R, which called copy number on the
# same 1000 cells, 1 Mb bins and blacklist with and without GC correction
# (AneuFinder edivisive, identical settings and seeds).
#
# AneuFinder scales each cell's relative segment depths to integer copy numbers
# by choosing the multiplier in [1.5, 6] that puts them closest to whole
# numbers. A genome and its doubling fit equally well, so systematic
# non-integer depth - which GC bias produces - decides the choice.
#
# Per-cell GC bias: |Spearman rho| between the uncorrected bin counts, divided
# by the K562 reference copy number, and bin GC content.
#
# Panels: (A) ploidy with vs without correction per cell; (B) the two ploidy
# distributions; (C) doubled and mis-scaled cells by GC-bias quintile.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
af_dir <- file.path(project_root, "results/K562_tree/aneufinder_gc/output")
toggle_dir <- file.path(af_dir, "gc_toggle")
output_dir <- file.path(project_root, "paper_figures/output/K562_gc_toggle_ploidy/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ploidy <- read.csv(file.path(toggle_dir, "ploidy_by_arm.csv"))

# ---- per-cell GC bias from the uncorrected bins ------------------------------

ref <- read.table(file.path(sc_dir, "K562_ref_ploidy_hg38.bed"), col.names = c("chr", "s", "e", "cn"))
ref_gr <- GRanges(ref$chr, IRanges(ref$s + 1, ref$e))
gc_bias <- sapply(ploidy$cell, function(cell) {
  f <- list.files(file.path(af_dir, "binned-GC"), paste0("^", cell, "\\.bam_"), full.names = TRUE)
  b <- get(load(f)); if (is(b, "GRangesList")) b <- b[[1]]    # carries the bin GC column
  r <- get(load(list.files(file.path(af_dir, "binned"), paste0("^", cell, "\\.bam_"), full.names = TRUE)))
  if (is(r, "GRangesList")) r <- r[[1]]
  h <- findOverlaps(resize(granges(r), 1, fix = "center"), ref_gr, select = "first")
  keep <- !is.na(h) & as.character(seqnames(r)) != "chrY" & r$counts > 0 & b$GC > 0
  abs(cor(r$counts[keep] / ref$cn[h][keep], b$GC[keep], method = "spearman"))
})
ploidy$gc_bias <- gc_bias[ploidy$cell]

# ---- scale classes -----------------------------------------------------------

main <- median(ploidy$ploidy_with_gc[ploidy$ploidy_with_gc >= 2.5 & ploidy$ploidy_with_gc < 3])
scale_class <- function(p) cut(p / main, c(0, 0.75, 1.25, 1.75, 2.5, 9),
                               labels = c("~0.5x", "~1x", "~1.5x", "~2x", ">2.5x"))
ploidy <- ploidy %>%
  mutate(class_without = scale_class(ploidy_without_gc),
         class_with = scale_class(ploidy_with_gc),
         gc_quintile = ntile(gc_bias, 5))
write.csv(ploidy, file.path(output_dir, "ploidy_with_without_gc.csv"), row.names = FALSE)

cat(sprintf("main-clone ploidy (with GC correction): %.2f\n", main))
cat("scale relative to the main clone, rows = without GC correction, columns = with:\n")
print(addmargins(table(without = ploidy$class_without, with = ploidy$class_with)))

by_q <- ploidy %>%
  group_by(gc_quintile) %>%
  summarise(cells = n(),
            gc_bias_range = sprintf("%.2f-%.2f", min(gc_bias), max(gc_bias)),
            doubled_without = mean(class_without %in% c("~2x", ">2.5x")),
            doubled_with = mean(class_with %in% c("~2x", ">2.5x")),
            misscaled_without = mean(class_without != "~1x"),
            misscaled_with = mean(class_with != "~1x"),
            .groups = "drop")
write.csv(by_q, file.path(output_dir, "misscaling_by_gc_bias.csv"), row.names = FALSE)
print(as.data.frame(by_q), digits = 3)

doubled <- ploidy$class_without %in% c("~2x", ">2.5x")
cat(sprintf("\ncells doubled without correction: %d; of these, %d are ~1x with correction\n",
            sum(doubled), sum(doubled & ploidy$class_with == "~1x")))
cat(sprintf("median GC bias: doubled-without-correction %.3f, ~1x-without-correction %.3f (Wilcoxon p = %.2g)\n",
            median(ploidy$gc_bias[doubled]), median(ploidy$gc_bias[ploidy$class_without == "~1x"]),
            wilcox.test(ploidy$gc_bias[doubled], ploidy$gc_bias[ploidy$class_without == "~1x"])$p.value))
cat(sprintf("Spearman, |ploidy without / main - 1| vs GC bias: %.3f\n",
            cor(abs(ploidy$ploidy_without_gc / main - 1), ploidy$gc_bias, method = "spearman")))

# ---- figure ------------------------------------------------------------------

guides_at <- main * c(0.5, 1, 1.5, 2)
pA <- ggplot(ploidy, aes(ploidy_without_gc, ploidy_with_gc, fill = gc_bias)) +
  geom_hline(yintercept = main, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = guides_at, color = "grey85", linewidth = 0.3, linetype = "dashed") +
  geom_point(shape = 21, size = 1.2, stroke = 0.1, alpha = 0.85) +
  scale_fill_gradient(low = "#f0f0f0", high = "#08306b", limits = c(0, 1), name = "Raw GC bias") +
  annotate("text", x = guides_at, y = max(ploidy$ploidy_with_gc) * 1.02,
           label = c("0.5x", "1x", "1.5x", "2x"), size = 2.2, color = "grey40", vjust = 0) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  xlab("Ploidy, no GC correction") + ylab("Ploidy, GC-corrected")

pB <- ploidy %>%
  pivot_longer(c(ploidy_without_gc, ploidy_with_gc), names_to = "arm", values_to = "p") %>%
  mutate(arm = factor(arm, c("ploidy_without_gc", "ploidy_with_gc"),
                      c("No GC correction", "GC-corrected"))) %>%
  ggplot(aes(p)) +
  geom_histogram(binwidth = 0.1, fill = "grey30") +
  geom_vline(xintercept = guides_at, color = "grey70", linewidth = 0.3, linetype = "dashed") +
  facet_wrap(~arm, ncol = 1, scales = "free_y") +
  theme_classic() + theme(strip.background = element_blank()) +
  xlab("Ploidy") + ylab("Cells")

pC <- by_q %>%
  select(gc_quintile, doubled_without, doubled_with, misscaled_without, misscaled_with) %>%
  pivot_longer(-gc_quintile, names_to = c("measure", "arm"), names_sep = "_") %>%
  mutate(measure = recode(measure, doubled = "Doubled (~2x or more)", misscaled = "Any wrong scale"),
         arm = factor(arm, c("without", "with"), c("No GC correction", "GC-corrected"))) %>%
  ggplot(aes(factor(gc_quintile), value, fill = arm)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("No GC correction" = "#cb181d", "GC-corrected" = "#08306b"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~measure, scales = "free_y") +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "bottom") +
  xlab("Raw GC-bias quintile (1 = least)") + ylab("Cells")

fig <- (pA | pB) / pC + plot_annotation(tag_levels = "A")
ggsave(file.path(output_dir, "gc_toggle_ploidy.svg"), fig, bg = "transparent", height = 7, width = 8)
ggsave(file.path(output_dir, "gc_toggle_ploidy.png"), fig, bg = "white", height = 7, width = 8, dpi = 200)

cat("\nWrote results to:", output_dir, "\n")
