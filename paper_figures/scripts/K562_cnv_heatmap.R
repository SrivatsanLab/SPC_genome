# Genome-wide copy-number heatmap and cell-cell similarity for the 1000
# sc_PolE_novaseq K562 cells, from the GC-corrected AneuFinder calls.
#
# Heatmap: cells x 1 Mb bins (chr1-22, X), absolute copy number on a diverging
# scale centred on 3 (K562 is near-triploid). Rows are ordered by hierarchical
# clustering (Ward.D2 on Manhattan distance of ploidy-rescaled copy number over
# unmasked bins). Side strips: population, ploidy group, and the cell's Pearson
# r with the K562 reference karyotype; a top strip marks repeat-masked bins.
#
# Similarity: for each pair of cells, the fraction of unmasked bins at the same
# copy number after rescaling each cell to the median ploidy, so a uniform
# scale offset does not count as a difference.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

sc_dir <- Sys.getenv("K562_SC_PROJECT",
                     "/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq")
result_csv <- file.path(project_root, "results/K562_tree/aneufinder_gc/output/result.csv")
mask_csv <- file.path(project_root, "results/K562_tree/masks/bin_mask_fraction.csv")
ref_csv <- file.path(project_root, "paper_figures/output/K562_cnv_vs_reference/per_cell_agreement.csv")
output_dir <- file.path(project_root, "paper_figures/output/K562_cnv_heatmap/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mask_min <- 0.25
pop_colors <- c("P_0" = "#260091", "P_1" = "#1E90FF", "P_2" = "#FFDB58",
                "P_3" = "#FF9D71", "P_4" = "#FF1B5E", "P_5" = "#E3E6E6")
ploidy_colors <- c("<2.5" = "#5b8fd6", "2.5-3.0" = "#d9d9d9", ">3.0" = "#d6604d")

# ---- data ------------------------------------------------------------------

cnv <- read.csv(result_csv, row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(sc_dir, "results/full_meta.csv"))
cells <- meta$X
bins <- cnv[, c("seqnames", "start", "end")]
cn <- as.matrix(cnv[, cells])

keep_chr <- bins$seqnames %in% c(paste0("chr", 1:22), "chrX")
bins <- bins[keep_chr, ]
cn <- cn[keep_chr, ]

mf <- read.csv(mask_csv)
masked <- mf$mask_frac[match(paste(bins$seqnames, bins$start),
                             paste(mf$seqnames, mf$start))] >= mask_min
modal <- apply(cn, 1, function(x) as.integer(names(which.max(table(x)))))
usable <- !masked & modal > 0

ploidy <- colMeans(cn)
ploidy_group <- cut(ploidy, c(0, 2.5, 3.0, 99), labels = names(ploidy_colors), right = FALSE)
ref <- read.csv(ref_csv) %>% filter(run == "new", bins == "All bins")
ref_r <- setNames(ref$pearson, ref$cell)[cells]
pop <- setNames(meta$pop, meta$X)[cells]

# rescale to the median ploidy, as in K562_cnv_event_tree.R
cn_norm <- round(sweep(cn, 2, median(ploidy) / ploidy, "*"))

# ---- ordering and similarity -------------------------------------------------

x <- t(cn_norm[usable, ])
hc <- hclust(dist(x, method = "manhattan"), method = "ward.D2")
ord <- hc$order

# fraction of usable bins at identical copy number, for every pair of cells
same <- matrix(0, length(cells), length(cells), dimnames = list(cells, cells))
for (b in seq_len(ncol(x))) {
  v <- x[, b]
  same <- same + outer(v, v, "==")
}
same <- same / ncol(x)

lt <- lower.tri(same)
grp <- as.character(ploidy_group)
pair_type <- function(a, b) ifelse(a == b, "same", "different")
pairs <- data.frame(
  similarity = same[lt],
  pop = pair_type(pop[row(same)[lt]], pop[col(same)[lt]]),
  both_main = grp[row(same)[lt]] == "2.5-3.0" & grp[col(same)[lt]] == "2.5-3.0",
  any_off = grp[row(same)[lt]] != "2.5-3.0" | grp[col(same)[lt]] != "2.5-3.0"
)
summ <- bind_rows(
  pairs %>% summarise(set = "all pairs", n = n(), median = median(similarity),
                      q25 = quantile(similarity, 0.25), q75 = quantile(similarity, 0.75)),
  pairs %>% filter(both_main) %>% summarise(set = "both main-peak cells", n = n(), median = median(similarity),
                      q25 = quantile(similarity, 0.25), q75 = quantile(similarity, 0.75)),
  pairs %>% filter(any_off) %>% summarise(set = "at least one off-peak cell", n = n(), median = median(similarity),
                      q25 = quantile(similarity, 0.25), q75 = quantile(similarity, 0.75)),
  pairs %>% filter(both_main, pop == "same") %>% summarise(set = "main peak, same population", n = n(),
                      median = median(similarity), q25 = quantile(similarity, 0.25), q75 = quantile(similarity, 0.75)),
  pairs %>% filter(both_main, pop == "different") %>% summarise(set = "main peak, different population", n = n(),
                      median = median(similarity), q25 = quantile(similarity, 0.25), q75 = quantile(similarity, 0.75))
)
# fraction of usable bins where each cell differs from the per-bin modal state
diff_from_modal <- colMeans(cn_norm[usable, ] != modal[usable])
write.csv(summ, file.path(output_dir, "pairwise_similarity_summary.csv"), row.names = FALSE)
write.csv(data.frame(cell = cells, pop = pop, ploidy = ploidy, ploidy_group = grp,
                     reference_r = ref_r, frac_bins_differ_from_modal = diff_from_modal,
                     cluster_order = match(seq_along(cells), ord)),
          file.path(output_dir, "per_cell_cnv_summary.csv"), row.names = FALSE)
print(as.data.frame(summ), digits = 3)
cat("fraction of usable bins differing from the modal state, by ploidy group:\n")
print(tapply(diff_from_modal, grp, median))

# ---- genome-wide heatmap -----------------------------------------------------

bins$idx <- seq_len(nrow(bins))
chr_bounds <- bins %>% group_by(seqnames) %>%
  summarise(first = min(idx), last = max(idx), .groups = "drop") %>%
  mutate(mid = (first + last) / 2, label = sub("chr", "", seqnames)) %>%
  arrange(first)
row_pos <- setNames(seq_along(ord), cells[ord])

cn_long <- data.frame(
  bin = rep(bins$idx, times = length(cells)),
  cell = rep(cells, each = nrow(bins)),
  cn = pmin(as.vector(cn), 8)
) %>% mutate(row = row_pos[cell])

cn_fill <- scale_fill_gradientn(
  colours = c("#08306b", "#6baed6", "#f0f0f0", "#fc9272", "#cb181d", "#67000d"),
  values = scales::rescale(c(0, 2, 3, 4, 6, 8)), limits = c(0, 8),
  breaks = 0:8, labels = c(0:7, "8+"), name = "Copy number"
)

p_heat <- ggplot(cn_long, aes(bin, row, fill = cn)) +
  geom_raster() +
  cn_fill +
  geom_vline(xintercept = chr_bounds$first[-1] - 0.5, color = "black", linewidth = 0.2) +
  scale_x_continuous(breaks = chr_bounds$mid, labels = chr_bounds$label, expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(), axis.text.x = element_text(size = 6)) +
  xlab(NULL) + ylab(sprintf("%d cells (clustered)", length(cells)))

p_mask <- ggplot(data.frame(bin = bins$idx, masked = masked), aes(bin, 1, fill = masked)) +
  geom_raster() +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey30"),
                    labels = c("FALSE" = "kept", "TRUE" = "repeat-masked"), name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() + theme(legend.position = "right")

strips <- data.frame(row = row_pos[cells], pop = pop, ploidy = grp, ref = ref_r)
strip <- function(var, scale, lab) {
  ggplot(strips, aes(1, row, fill = .data[[var]])) + geom_raster() + scale +
    scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
    theme_void() + labs(title = lab) +
    theme(plot.title = element_text(size = 6, angle = 90, hjust = 0, vjust = 0.5))
}
p_pop <- strip("pop", scale_fill_manual(values = pop_colors, name = "Population"), "Population")
p_ploidy <- strip("ploidy", scale_fill_manual(values = ploidy_colors, name = "Ploidy"), "Ploidy")
p_ref <- strip("ref", scale_fill_gradient(low = "white", high = "black", limits = c(0, 1),
                                          name = "r with K562\nreference"), "Reference r")

heat <- p_mask + p_pop + p_ploidy + p_ref + p_heat +
  plot_layout(design = "
####AAAAAAAAAAAAAAAAAAAAAA
BCD#EEEEEEEEEEEEEEEEEEEEEE
", heights = c(0.02, 1), guides = "collect") &
  theme(legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
ggsave(file.path(output_dir, "cnv_heatmap_all_cells.png"), heat, height = 9, width = 12, dpi = 200, bg = "white")
ggsave(file.path(output_dir, "cnv_heatmap_all_cells.svg"), heat, height = 9, width = 12, bg = "white")

# ---- cell x cell similarity --------------------------------------------------

sim_long <- data.frame(
  i = rep(seq_along(ord), times = length(ord)),
  j = rep(seq_along(ord), each = length(ord)),
  s = as.vector(same[ord, ord])
)
p_sim <- ggplot(sim_long, aes(i, j, fill = s)) +
  geom_raster() +
  scale_fill_gradient(low = "white", high = "#08306b", limits = c(0, 1),
                      name = "Fraction of bins\nat equal copy number") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) +
  coord_equal() +
  theme_void() +
  theme(legend.title = element_text(size = 7), legend.text = element_text(size = 6))
ggsave(file.path(output_dir, "cell_similarity_heatmap.png"), p_sim, height = 6, width = 7, dpi = 200, bg = "white")

cat("\nWrote results to:", output_dir, "\n")
