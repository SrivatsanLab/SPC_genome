# Draw the population-level ("bulk") SBS trees from K562_sbs_tree_rebuild.R
# (original 700 sites, pruned of germline/chrM/unplaced sites, and reselected
# on unflagged sites), in the style of K562_consensus_tree.R: vertical
# dendrogram, divergence down the page, tips as filled points in the P_0..P_5
# palette, dashed depth guides, no legend. One 3 x 3 in PDF per tree and all
# three side by side.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
})

tree_dir <- file.path(project_root, "paper_figures/output/K562_sbs_tree_rebuild")
output_dir <- tree_dir

pop_colors <- c("P_0" = "#260091", "P_1" = "#1E90FF", "P_2" = "#FFDB58",
                "P_3" = "#FF9D71", "P_4" = "#FF1B5E", "P_5" = "#E3E6E6")

plot_pop_tree <- function(tree, title = NULL) {
  depth <- max(node.depth.edgelength(tree))
  breaks <- seq(0, ceiling(depth * 10) / 10, by = 0.1)
  tips <- ggtree(tree)$data
  tips <- tips[tips$isTip, ]
  p <- ggtree(tree, layout = "roundrect", size = 1) +
    geom_point(data = tips, aes(x = x, y = y, fill = label),
               shape = 21, size = 6, stroke = 1, color = "black") +
    scale_fill_manual(values = pop_colors) +
    scale_x_reverse(breaks = breaks) +
    coord_flip() +
    theme_tree2() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 10, color = "black"),
          axis.line.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          plot.title = element_text(size = 9))
  if (!is.null(title)) p <- p + ggtitle(title)
  p$layers <- c(geom_vline(xintercept = breaks, color = "grey", linetype = "dashed",
                           linewidth = 0.5, alpha = 0.6), p$layers)
  p
}

sets <- c(original = "original (700 sites)", pruned = "pruned (333 sites)",
          reselected = "reselected (700 clean sites)")
plots <- list()
for (k in names(sets)) {
  tr <- read.tree(file.path(tree_dir, paste0("population_tree_", k, ".newick")))
  plots[[k]] <- plot_pop_tree(tr)
  ggsave(file.path(output_dir, paste0("population_tree_", k, ".pdf")), plots[[k]], width = 3, height = 3)
}
combined <- wrap_plots(lapply(names(sets), function(k) plots[[k]] + ggtitle(sets[[k]])), nrow = 1)
ggsave(file.path(output_dir, "population_trees.pdf"), combined, width = 9, height = 3.4)
ggsave(file.path(output_dir, "population_trees.png"), combined, width = 9, height = 3.4, dpi = 200, bg = "white")
cat("Wrote population tree PDFs to", output_dir, "\n")
