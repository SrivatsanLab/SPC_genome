# K562 grouped bootstrap consensus tree.
#
# Ported from notebooks/K562_tree.ipynb cell 97, which draws the tree with
# baltic + matplotlib. This renders the same tree with ggtree, following the
# style in draw_trees.R (roundrect layout, tips as geom_point coloured by
# population, no legend).
#
# The tree itself is built upstream in the notebook: per-population consensus
# over bootstrap replicates, summarised with `sumtrees --summary-target=mcct`,
# then written out as newick. Internal node labels carry bootstrap support.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(ggtree)
library(ape)

data_dir <- file.path(project_root, "paper_figures/data/K562_tree/")
output_dir <- file.path(project_root, "paper_figures/output/K562_consensus_tree/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# same P_0..P_5 palette as draw_trees.R, plus the WT outgroup
pop_colors <- c(
  "P_0" = "#260091",
  "P_1" = "#1E90FF",
  "P_2" = "#FFDB58",
  "P_3" = "#FF9D71",
  "P_4" = "#FF1B5E",
  "P_5" = "#E3E6E6",
  "WT"  = "black"
)

grid_breaks <- seq(0, 0.6, by = 0.1)

# dendropy writes the population names single-quoted ('P_5'), and ape keeps
# the quotes in the label, so strip them before matching against pop_colors
read_consensus <- function(path) {
  tree <- read.tree(path)
  tree$tip.label <- gsub("^'|'$", "", tree$tip.label)
  tree
}

plot_consensus_tree <- function(tree) {
  tips <- ggtree(tree)$data %>% dplyr::filter(isTip)

  p <- ggtree(tree, layout = "roundrect", size = 1) +
    geom_point(
      data = tips,
      aes(x = x, y = y, fill = label),
      shape = 21, size = 6, stroke = 1.2, color = "black"
    ) +
    scale_fill_manual(values = pop_colors) +
    # Vertical dendrogram: divergence runs down the page, tips along the
    # bottom. This is what layout_dendrogram() does, spelled out so the depth
    # breaks can be set - calling it as well would apply scale_x_reverse twice
    # and detach the tip points from the branches.
    scale_x_reverse(breaks = grid_breaks) +
    coord_flip() +
    theme_tree2() +
    theme(
      legend.position = "none",
      # depth axis, vertical after the flip
      axis.text.y = element_text(size = 16, color = "black"),
      axis.line.y = element_blank(),
      # tip-index axis, meaningless here
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    )

  # prepend so the guides sit behind the branches, matching the notebook's
  # zorder = 0; layers added the usual way would draw on top of the tree
  p$layers <- c(
    geom_vline(xintercept = grid_breaks, color = "grey",
               linetype = "dashed", linewidth = 0.75, alpha = 0.6),
    p$layers
  )
  p
}

consensus <- read_consensus(file.path(data_dir, "grouped_bootstrap_consensus.newick"))

plot_consensus_tree(consensus)
ggsave(file.path(output_dir, "grouped_consensus_tree.pdf"),
       bg = "transparent",
       height = 5, width = 5)

# UPGMA version of the same consensus, for comparison
upgma_path <- file.path(data_dir, "grouped_bootstrap_consensus_upgma.newick")
if (file.exists(upgma_path)) {
  plot_consensus_tree(read_consensus(upgma_path))
  ggsave(file.path(output_dir, "grouped_consensus_tree_upgma.pdf"),
         bg = "transparent",
         height = 5, width = 5)
}

# Bootstrap support sits in the internal node labels; report it rather than
# crowding the figure, matching the notebook (its support scatter is
# commented out in cell 97).
support <- data.frame(
  node = seq_along(consensus$node.label) + length(consensus$tip.label),
  support = consensus$node.label
)
print(support, row.names = FALSE)
write.csv(support, file.path(output_dir, "bootstrap_support.csv"), row.names = FALSE)

cat("\nWrote tree to:", output_dir, "\n")
