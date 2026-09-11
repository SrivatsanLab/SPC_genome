# Worm6 haplotype assignment UMAP.
#
# Ported from notebooks/worm6_final_haplotype_assignment.ipynb cell 58, which
# draws the embedding with seaborn. This renders the same points in the
# paper_figures house style.
#
# The embedding itself is computed upstream in the notebook and is NOT saved to
# disk by it, so this script cannot regenerate it. Add one line after the cell
# that computes U (cell 57) and re-run that notebook:
#
#   pd.DataFrame({"UMAP1": U[:, 0], "UMAP2": U[:, 1],
#                 "donor": call_k, "purity": purity[keep]}) \
#     .to_csv("../results/worm6_final/figures/clean_ind_umap_coords.csv",
#             index=False)
#
# then stage that CSV with sync_missing_inputs.sh. Re-running the UMAP from
# scratch is not equivalent: the layout depends on umap-learn / numba / numpy
# versions even with random_state = 1, so it would not match the published
# figure.
#
# Upstream pipeline, for reference:
#   keep    = ~background & donor not in worm02/04/14/19, background = purity < BG_PURITY
#   Xk      = per-core-variant VAF (AD/DP at DP >= 1), column centred, NaN -> 0
#   emb     = TruncatedSVD(n_components = 18, random_state = 1)
#   U       = umap.UMAP(n_neighbors = 30, min_dist = 0.3, random_state = 1)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)

data_dir <- file.path(project_root, "paper_figures/data/worm6_final/")
output_dir <- file.path(project_root, "paper_figures/output/worm6_haplotype_umap/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

coords_path <- file.path(data_dir, "clean_ind_umap_coords.csv")

if (!file.exists(coords_path)) {
  message("Skipping worm6 UMAP: ", coords_path, " not found.\n",
          "  The notebook does not save the embedding; see the header of this ",
          "script for the one line that exports it.")
} else {

  # extended_colors[1:] from the notebook, i.e. the palette minus its leading
  # near-black, assigned to the worms in ascending order and recycled
  extended_colors <- c(
    "#FF8DAF", "#FF1B5E",
    "#FFEDA0", "#FFDB57",
    "#A4C8FE", "#438CFD",
    "#FFCF80", "#FF8C00",
    "#C4F9AE", "#80F15E",
    "#A0FFFF", "#03FFFF",
    "#A99FD4", "#5941A9",
    "#F4BFFE", "#E980FC",
    "#E8E4F0", "#C4BDD6",
    "#E8C89A", "#C17F3E"
  )

  umap_df <- read.csv(coords_path, stringsAsFactors = FALSE)

  worms <- sort(unique(umap_df$donor[startsWith(umap_df$donor, "worm")]))
  others <- sort(unique(umap_df$donor[!startsWith(umap_df$donor, "worm")]))

  pal <- setNames(
    extended_colors[(seq_along(worms) - 1) %% length(extended_colors) + 1],
    worms
  )
  # non-worm donors are greyed out, as in the notebook
  pal <- c(pal, setNames(rep("#BFBFBF", length(others)), others))

  umap_df$donor <- factor(umap_df$donor, levels = c(worms, others))

  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = donor)) +
    geom_point(shape = 21, size = 1.8, colour = "black", stroke = 0.2) +
    scale_fill_manual(values = pal) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.key.height = unit(14, "pt")
    ) +
    guides(fill = guide_legend(override.aes = list(size = 3.5))) +
    xlab("UMAP 1") +
    ylab("UMAP 2")
  ggsave(file.path(output_dir, "clean_ind_umap.svg"),
         bg = "transparent",
         height = 6, width = 7.5)

  # purity version, matching cell 59
  if ("purity" %in% names(umap_df)) {
    ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = purity)) +
      geom_point(shape = 21, size = 1.8, colour = "black", stroke = 0.2) +
      scale_fill_viridis_c(option = "magma", limits = c(0, 1), name = "purity") +
      theme_classic() +
      theme(
        axis.title = element_text(size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      ) +
      xlab("UMAP 1") +
      ylab("UMAP 2")
    ggsave(file.path(output_dir, "clean_ind_umap_purity.svg"),
           bg = "transparent",
           height = 6, width = 7)
  }

  cat("Wrote UMAP to:", output_dir, "\n")
}
