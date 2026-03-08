project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(dplyr)
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(tidyverse)
library(treeio)

# Set a global seed for reproducibility
set.seed(50)
setwd(file.path(project_root, "paper_figures/data/SingleCellTrees/"))
# Function definitions ----------------------------------------------------


# Generate a tree and compute depth
generate_tree_data <- function(n_tips, tree_type = "balanced") {
  tree <- stree(n_tips, type = tree_type)
  ggtree_obj <- ggtree(tree)
  tree_data <- ggtree_obj$data %>%
    mutate(depth = node.depth(tree)[node])
  list(tree = tree, tree_data = tree_data)
}

# Randomly assign types (A/B) to internal nodes at a specific depth
assign_types_to_depth_nodes <- function(tree_data, target_depth) {
  tree_data %>%
    dplyr::filter(depth == target_depth, isTip == FALSE) %>%
    mutate(type = sample(c("A", "B"), n(), replace = TRUE))
}

# Use phangorn to get descendants of a node
get_descendants <- function(tree, node) {
  phangorn::Descendants(tree, node, type = "all")
}

# Propagate A/B label from a node to all its descendants
propagate_labels <- function(tree_data, tree, labeled_nodes) {
  tree_data$type <- NA
  for (i in seq_len(nrow(labeled_nodes))) {
    node <- labeled_nodes$node[i]
    label <- labeled_nodes$type[i]
    descendants <- get_descendants(tree, node)
    tree_data$type[tree_data$node %in% c(node, descendants)] <- label
  }
  tree_data
}

# Plot the tree with colored types
plot_colored_tree <- function(tree, tree_data, layout = "roundrect") {
  ggtree(tree, layout = layout) +
    geom_tree(color = NA) +
    geom_point(
      data = tree_data %>% dplyr::filter(!is.na(type)),
      aes(x = x, y = y, color = type), size = 1.5
    ) +
    scale_color_manual(values = c("A" = "#ff1a5e", "B" = "#79AFEA")) +
    theme_tree(legend.position = "none")
}

# Save a plot
save_tree_plot <- function(plot, filename, height = 5, width = 2) {
  ggsave(plot = plot, filename = filename, height = height, width = width)
}

# Make Plots --------------------------------------------------------------

# Output directory
output_dir <- file.path(project_root, "paper_figures/output/draw_trees/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Early split tree (simulate split at depth 64)
early <- generate_tree_data(128, "balanced")
early_depth_nodes <- assign_types_to_depth_nodes(early$tree_data, target_depth = 64)
early_labeled <- propagate_labels(early$tree_data, early$tree, early_depth_nodes)
early_plot <- plot_colored_tree(early$tree, early_labeled)
early_plot
save_tree_plot(early_plot, file.path(output_dir, "early_split_tree.pdf"))

# 2. Late split tree (simulate split at depth 4)
late <- generate_tree_data(128, "balanced")
late_depth_nodes <- assign_types_to_depth_nodes(late$tree_data, target_depth = 4)
late_labeled <- propagate_labels(late$tree_data, late$tree, late_depth_nodes)
late_plot <- plot_colored_tree(late$tree, late_labeled)
late_plot
save_tree_plot(late_plot, file.path(output_dir, "late_split_tree.pdf"))

# 3. Custom labeled x-axis tree (non-colored, custom ticks)
custom_tree <- stree(6, type = "left")
custom_plot <- ggtree(custom_tree, layout = "roundrect", size = 1, color = "white") +
  theme_tree() +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5),
                     labels = c("9", "11", "13", "15", "17", "19")) +
  theme(
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  ) +
  theme_void()

save_tree_plot(custom_plot, file.path(output_dir, "our_tree.pdf"), height = 1.5, width = 5)




# Advanced Schematic ------------------------------------------------------

complex_tree = generate_tree_data(n_tips = 256)

sampled_node_1 = 
  complex_tree[[2]] %>%
  filter(depth == 16) %>%
  sample_n(1) %>%
  pull(node)

descendants_1 = offspring(.data = complex_tree[[1]],.node = sampled_node_1)

sampled_node_2 = 
  complex_tree[[2]] %>%
  filter(depth == 8) %>%
  pull(node) %>%
  intersect(descendants_1) %>%
  sample(1)

descendants_2 = offspring(.data = complex_tree[[1]],.node = sampled_node_2)

sampled_node_3 = 
  complex_tree[[2]] %>%
  filter(depth == 4) %>%
  pull(node) %>%
  intersect(descendants_2) %>%
  sample(1)

descendants_3 = offspring(.data = complex_tree[[1]],.node = sampled_node_3)

custom_plot = ggtree(complex_tree$tree,layout = "roundrect") +
  geom_point(data = complex_tree$tree_data %>%
               filter(node %in% c(sampled_node_1, descendants_1)),
             aes(x = x,
                 y = y),
             color = "#ff1a5e",
             stroke = 0,
             size = 4)  +
  geom_point(data = complex_tree$tree_data %>%
               filter(node %in% c(sampled_node_2, descendants_2)),
             aes(x = x-0.2,
                 y = y),
             color = "#79AFEA",
             stroke = 0,
             size = 4)  +
  geom_point(data = complex_tree$tree_data %>%
               filter(node %in% c(sampled_node_3, descendants_3)),
             aes(x = x-0.4,
                 y = y),
             color = "#FFDB58",
             stroke = 0,
             size = 4)  +
  theme_void()
custom_plot




plot_ILS_descendants <- function(tree_obj,
                                 depth_sequences = list(c(16, 8, 4)),
                                 color_sequences = list(c("#ff1a5e", "#79AFEA", "#FFDB58")),
                                 offset_step = 0.1,
                                 seed = 50,
                                 tree_size = 1,
                                 point_size = 2) {
  set.seed(seed)
  tree <- tree_obj$tree
  tree_data <- tree_obj$tree_data
  tree_data$color = NA
  counter = 0
  p <- ggtree(tree, layout = "roundrect") + theme_void()
  for (seq_idx in seq_along(depth_sequences)) {
    depth_sequence <- depth_sequences[[seq_idx]]
    colors <- color_sequences[[seq_idx]]
    descendants <- NULL
    counter = 0
    
    for (i in seq_along(depth_sequences[[seq_idx]])) {
      depth <- depth_sequences[[seq_idx]][i]
      # Get only the nodes at that depth
      nodes_at_depth <- 
        tree_data %>% 
        filter(depth == depth_sequence[i])

      # If this is the highest depth, pick one node 
      if (i == 1) {
        sampled_node <- 
          nodes_at_depth %>%
          sample_n(1) %>% 
          pull(node)
      } else {
        sampled_node <-
        nodes_at_depth %>%
          pull(node) %>%
          intersect(descendants) %>%
          sample(1)
        
      }
      
      descendants <- offspring(.data = tree, .node = sampled_node)
      highlight_nodes <- c(sampled_node, descendants)
      tree_data =
        tree_data %>%
        mutate(color = ifelse(node %in% highlight_nodes,colors[i],color))

      p <- p + geom_point(
        data = tree_data %>%
          filter(node %in% highlight_nodes) %>%
          mutate(x_offset = x + (i-1) * offset_step),
        aes(x = x_offset, y = y),
        color = colors[i],
        size = point_size,
        stroke = 0)
        
        
      
      print(counter)
    }
  }
  
  return(list(p,tree_data))
}

complex_tree <- generate_tree_data(n_tips = 128)
depths_list <- list(
  c(128,16,4),
  c(64,16,4),
  c(32, 4, 2),
  c(32, 8, 2),
  c(16, 4, 2)
)

colors_list <- list(
  c("red", "blue", "yellow"),
  c("#9A68A4", "#34C924", "pink"),
  c("#8B5E3C", "#A3B18A", "#DAD7CD"),
  c("#E63946", "#457B9D", "#F1FAEE"),
  c("darkgreen", "#F3722C", "black")
  
)

res = plot_ILS_descendants(complex_tree, depths_list, colors_list,tree_size = 0.5,seed = 49,offset_step = 0.2)
custom_plot = res[[1]]
tree_data = res[[2]]
custom_plot
save_tree_plot(custom_plot, file.path(output_dir, "sparse_tree.pdf"))


split_tree_data = 
  tree_data %>%
  filter(depth == 2) %>%
  mutate(line_group = sample(c(1, 2), size = n(), replace = TRUE))

  
# Analyze our single cell tree --------------------------------------------

library(ggtree)

# Read the Newick file
sc_tree <- read.tree("sc_test.newick")

# Plot the tree
p = ggtree(sc_tree,layout = "roundrect")
tree_df = p$data

cell_metadata = 
  read.csv("snp_counts.csv", 
         header = T,
         sep = ",") %>%
  dplyr::rename(label = X)

tree_df = 
  left_join(tree_df, cell_metadata)

# The root node has the itself as a parent
find_self_parent_root <- function(df) {
  df %>%
    filter(parent == node) %>%
    pull(node)
}

root_node = find_self_parent_root(tree_df)

tree_df <- tree_df %>% mutate(parent2 = ifelse(parent == node, NA_integer_, parent))

compute_hops_to_root <- function(df, root_node) {
  # Lookup: node → parent2
  parent_map <- df %>%
    select(node, parent2) %>%
    deframe()
  
  # Count how many steps to reach the specified root
  hops_to_root <- sapply(names(parent_map), function(node) {
    current <- node
    hops <- 0
    visited <- character()
    
    while (!is.na(current) && current != root_node) {
      if (current %in% visited) {
        warning(paste("Cycle detected at node:", current))
        return(NA_integer_)
      }
      visited <- c(visited, current)
      current <- as.character(parent_map[[current]])
      hops <- hops + 1
    }
    
    # Return NA if we never reach the root
    if (is.na(current)) return(NA_integer_)
    return(hops)
  })
  
  df %>%
    mutate(hops_to_root = hops_to_root[as.character(node)])
}

tree_df = compute_hops_to_root(tree_df,root_node = root_node)
tree_df <- tree_df %>%
  mutate(pop_factor = recode(pop,
                             "P_0" = 6,
                             "P_1" = 5,
                             "P_2" = 4,
                             "P_3" = 3,
                             "P_4" = 2,
                             "P_5" = 1
  ))

pop_colors <- c(
  "P_0" = "#260091",
  "P_1" = "#1E90FF",
  "P_2" = "#FFDB58",
  "P_3" = "#FF9D71",
  "P_4" = "#FF1B5E",
  "P_5" = "#E3E6E6"
)

tree_df %>%
  filter(isTip) %>%
  ggplot() +
  geom_violin(aes(x = pop_factor,
                   y = hops_to_root,
                   fill = pop),
              alpha = 0.5,
              color = "black") +
  geom_jitter(aes(x = pop_factor,
                  y = hops_to_root,
                  fill = pop),
              shape = 21,
              color = "black",
              width = 0.25,
              height = 0.75) +
  
  # geom_smooth(aes(x = pop_factor,
  #             y = hops_to_root),
  #         method = "lm",
  #         color = "red") +
  ylab("Depth") +
  xlab("Population")+ 
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13)) +
  scale_fill_manual(values = pop_colors) 
  # scale_color_manual(values = pop_colors)
ggsave(paste(output_dir,"/DepthPerPopulation.pdf",sep = ""),
       height = 2.5,
       width = 2.5)

tips = 
  tree_df %>%
  filter(isTip)
cor.test(tips$pop_factor,tips$hops_to_root)

tree_data %>%
  arrange(-desc(x))

ggtree(sc_tree,layout = "roundrect", size = 0.15) +
  geom_point(data = tree_df %>%
               dplyr::filter(isTip),
             aes(x = x,
                 y = y,
                 color = pop),
             size = .75,
             stroke = 0) +
  geom_vline(xintercept = 0.2,
             color = NA,
             size = 0)+
  scale_color_manual(values = pop_colors) +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "sc_tree.pdf"),
       height = 6,
       width = 2)

# White theme

ggtree(sc_tree,layout = "roundrect", size = 0.15, color = "white") + 
  geom_point(data = tree_df %>%
               dplyr::filter(isTip),
             aes(x = x,
                 y = y,
                 color = pop),
             size = .75,
             stroke = 0) +
  geom_vline(xintercept = 0.2,
             color = NA,
             size = 0)+
  scale_color_manual(values = pop_colors) +
  theme_void() +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position = "none")
  #theme(legend.position = "none")
ggsave(file.path(output_dir, "sc_tree_desktop_legacy.pdf"),
       height = 6,
       width = 2,bg = "transparent")

  # Assume `tree` is a phylo object and `tip_types` is a named vector:
# names(tip_types) == tree$tip.label

# 1. Convert tree_df to phylo
sc_tree <- as.phylo(tree_df)

# 2. Extract tip labels and match them back to tree_df to get types
tip_nodes <- tree_df %>% filter(isTip) %>% filter(!is.na(pop)) 

tip_labels <- tip_nodes$label
names(tip_labels) <- tip_nodes$node

# Ensure type or pop column exists
stopifnot("pop" %in% colnames(tip_nodes))

tip_types <- tip_nodes$pop
names(tip_types) <- tip_nodes$label

# 3. Compute pairwise distances between all tip labels
cophenetic_mat <- cophenetic(sc_tree)

# 4. For each tip, find the closest other tip and compare types
closest_type_check <- map_dfr(rownames(cophenetic_mat), function(tip_label) {
  dists <- cophenetic_mat[tip_label, ]
  dists <- dists[names(dists) != tip_label]  # remove self
  closest_label <- names(which.min(dists))
  
  tibble(
    tip = tip_label,
    closest_tip = closest_label,
    same_type = tip_types[tip_label] == tip_types[closest_label],
    tip_type = tip_types[tip_label],
    closest_type = tip_types[closest_label],
    distance = dists[closest_label]
  )
})

print(closest_type_check)
closest_type_check %>% 
  group_by(tip_type) %>%
  summarise(sum(same_type))


# Try 1 not grouped by type -----------------------------------------------

# Assume `tree` is a phylo object and `tip_types` is a named vector:
# names(tip_types) == tree$tip.label

# 1. Convert tree_df to phylo
sc_tree <- as.phylo(tree_df)

# 2. Extract tip labels and match them back to tree_df to get types
tip_nodes <- tree_df %>% filter(isTip)
tip_labels <- tip_nodes$label
names(tip_labels) <- tip_nodes$node

# Ensure type or pop column exists
stopifnot("pop" %in% colnames(tip_nodes))

tip_types <- tip_nodes$pop
names(tip_types) <- tip_nodes$label

# 3. Compute pairwise distances between all tip labels
cophenetic_mat <- cophenetic(sc_tree)

# 4. Observed value
observed_matches <- map_lgl(names(tip_types), function(tip_label) {
  dists <- cophenetic_mat[tip_label, ]
  dists <- dists[names(dists) != tip_label]
  closest_label <- names(which.min(dists))
  tip_types[tip_label] == tip_types[closest_label]
})
observed_matches
observed_fraction <- mean(observed_matches)

# 5. Null distribution by random shuffling
set.seed(42)  # for reproducibility
n_permutations <- 1000
null_fractions <- replicate(n_permutations, {
  shuffled_types <- sample(tip_types)  # shuffle types
  names(shuffled_types) <- names(tip_types)
  mean(map_lgl(names(tip_types), function(tip_label) {
    dists <- cophenetic_mat[tip_label, ]
    dists <- dists[names(dists) != tip_label]
    closest_label <- names(which.min(dists))
    shuffled_types[tip_label] == shuffled_types[closest_label]
  }))
})

# 6. Report
cat("Observed fraction of same-type nearest neighbors:", round(observed_fraction, 3), "\n")
cat("Null mean ± sd:", round(mean(null_fractions), 3), "±", round(sd(null_fractions), 3), "\n")
cat("P-value (empirical):", mean(null_fractions >= observed_fraction), "\n")

# Convert null values to a data frame
null_df <- data.frame(fraction = null_fractions)

# Create the plot
ggplot(null_df, aes(x = fraction)) +
  geom_histogram(bins = 50, fill = "gray70", color = "black") +
  geom_vline(xintercept = observed_fraction, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Null distribution of same-type matches",
    x = "Fraction of NN's with same pop (null)",
    y = "Count"
  ) +
  theme_minimal()

# Save it
ggsave(file.path(output_dir, "null_distribution_of_same_type.pdf"), width = 4, height = 4)

# Try 2 grouped by type with a null distrbution ---------------------------

library(dplyr)
library(purrr)
library(tibble)
library(pbapply)  # for progress bar
library(ape)

# Step 0: Cophenetic matrix
cophenetic_mat <- cophenetic(sc_tree)
tip_labels <- rownames(cophenetic_mat)

# Function to get closest neighbor match info
get_closest_type_check <- function(tip_types) {
  map_dfr(tip_labels, function(tip_label) {
    dists <- cophenetic_mat[tip_label, ]
    dists <- dists[names(dists) != tip_label]
    closest_label <- names(which.min(dists))
    
    tibble(
      tip = tip_label,
      closest_tip = closest_label,
      same_type = tip_types[tip_label] == tip_types[closest_label],
      tip_type = tip_types[tip_label],
      closest_type = tip_types[closest_label],
      distance = dists[closest_label]
    )
  })
}

# Step 1: Observed data
observed_check <- get_closest_type_check(tip_types)

observed_summary <- observed_check %>%
  group_by(tip_type) %>%
  summarise(
    n = n(),
    frac_same = mean(same_type),
    .groups = "drop"
  )

# Step 2: Null permutations with progress bar
n_perm <- 25
set.seed(42)

null_summaries_long <- pblapply(seq_len(n_perm), function(i) {
  shuffled_types <- sample(tip_types)
  names(shuffled_types) <- names(tip_types)
  
  null_check <- get_closest_type_check(shuffled_types)
  
  null_check %>%
    group_by(tip_type) %>%
    summarise(
      frac_same = mean(same_type),
      .groups = "drop"
    ) %>%
    mutate(iteration = i)
}) %>% bind_rows()

plot_df <- null_summaries_long %>%
  left_join(observed_summary %>% select(tip_type, observed_frac = frac_same), by = "tip_type")

ggplot(plot_df, aes(x = frac_same)) +
  geom_histogram(bins = 50, fill = "gray80", color = "black") +
  geom_vline(aes(xintercept = observed_frac), color = "red", linetype = "dashed") +
  facet_wrap(~ tip_type, scales = "free_y") +
  labs(
    title = "Null distribution of same-type nearest neighbors",
    x = "Fraction of same-type neighbors (null)",
    y = "Count"
  ) +
  theme_minimal()

# Run a Chi Square Test ---------------------------------------------------

# 1. Observed counts
observed_counts <- observed_check %>%
  group_by(tip_type) %>%
  summarise(observed_same_type = sum(same_type), .groups = "drop")

# 2. Expected counts from null
expected_counts <- null_summaries_long %>%
  group_by(tip_type) %>%
  summarise(expected_same_type = mean(frac_same) * n(), .groups = "drop") %>%
  left_join(observed_counts, by = "tip_type")

# 3. Prepare vectors for chi-square
observed_vec <- expected_counts$observed_same_type
expected_vec <- expected_counts$expected_same_type
names(observed_vec) <- expected_counts$tip_type

# 4. Chi-square test
chisq_test <- chisq.test(x = observed_vec, p = expected_vec / sum(expected_vec))

print(chisq_test)

# Null summary: mean and sd of expected same-type counts
expected_stats <- null_summaries_long %>%
  group_by(tip_type) %>%
  summarise(
    expected_mean = mean(frac_same) * n(),
    expected_sd = sd(frac_same) * n(),
    .groups = "drop"
  )
plot_df <- observed_check %>%
  group_by(tip_type) %>%
  summarise(observed_same_type = sum(same_type), .groups = "drop") %>%
  left_join(expected_stats, by = "tip_type") %>%
  pivot_longer(cols = c(observed_same_type, expected_mean), 
               names_to = "type", values_to = "count") %>%
  mutate(type = recode(type,
                       "observed_same_type" = "Observed",
                       "expected_mean" = "Expected"))

# Make the plot for the figure ---------------------------------------------------

plot_df <- plot_df %>%
  mutate(
    ymin = ifelse(type == "Expected", count - expected_stats$expected_sd[match(tip_type, expected_stats$tip_type)], NA),
    ymax = ifelse(type == "Expected", count + expected_stats$expected_sd[match(tip_type, expected_stats$tip_type)], NA)
  )

ggplot(plot_df) +
  geom_point(data = plot_df %>%
               filter(type == "Observed"),
             aes(x = fct_rev(tip_type),
                 y = count,
                 fill = tip_type),
             shape = 21,
             size = 4,
             stroke = 0.5,
             color = "black") +
  geom_point(data = plot_df %>%
               filter(type == "Expected"),
             aes(x = fct_rev(tip_type),
                 y = count),
             fill = 'white',
             shape = 21,
             stroke = 0.5,
             color = "black",
             size = 2.5) +
  geom_errorbar(data = plot_df %>%
                  filter(type == "Expected"),
    aes(x = fct_rev(tip_type),
        y = count, 
        ymin = ymin, 
        ymax = ymax),
    position = position_dodge(width = 0.9),
    width = 0.05) +
  theme(legend.position = "none") +
  theme_classic() +
  xlab("Population")+ 
  ylab("Homotypic Neighbors") +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black"),
        plot.background = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = pop_colors) 
ggsave(file.path(output_dir, "SameNeighbor.pdf"),
       height = 2.5,
       width = 2.5)



# CNV data ----------------------------------------------------------------
library(dplyr)
library(purrr)
library(tibble)
library(pbapply)  # for progress bar
library(ape)
library(ggtree)
library(ggplot2)
setwd(file.path(project_root, "paper_figures/data/Anneufinder/"))


# Read the Newick file
cnv_tree <- read.tree("cnv_sc_test.newick")
sc_tree <- read.tree("../SingleCellTrees/sc_test.newick")

# Plot the tree
p_cnv = ggtree(cnv_tree,layout = "roundrect")
p_sc = ggtree(sc_tree,layout = "roundrect")

tree_df = p_sc$data

cell_metadata = 
  read.csv("full_meta.csv", 
           header = T,
           sep = ",") %>%
  dplyr::rename(label = X)

tree_df = 
  left_join(tree_df, cell_metadata)


ggtree(sc_tree,layout = "roundrect", size = 0.15, color = 'white') +
  geom_point(data = tree_df %>%
               dplyr::filter(isTip),
             aes(x = x,
                 y = y,
                 color = (round(ploidy))),
             size = 1,
             stroke = 0) +
  scale_color_viridis_c(name = "Ploidy") +
  theme_void()
ggsave(file.path(output_dir, "test.pdf"),
       height = 6,
       width = 2)

cell_metadata %>% 
  ggplot() +
  geom_histogram(aes(x = ploidy,
                     fill = as.factor(round(ploidy))),
                 binwidth = 0.1,
                 color = "black",
                 size =0.2) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7),labels = c("1N","2N","3N","4N","5N","6N","7N")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 13, color = "black")) +
  scale_fill_brewer(palette = "Spectral")

ggsave(paste(output_dir,"/ploidy_histogram.pdf",sep = ""),
       height = 2.5,
       width = 2.5)

ggplot() +
  geom_jitter(data = tree_df %>%
               filter(isTip,),
             aes(y =ploidy-.1/50,
                 x = y,
                 fill = factor(round(ploidy))),
             height = 0,
             width = 0,
             shape = 21,
             size = 1.5,
             color = "black",
             stroke = 0.1) +
  scale_fill_brewer(palette = "Spectral") + 
  scale_color_brewer(palette = "Spectral") + 
  theme_classic() +
  scale_y_continuous(breaks = c(1,2,3,4,5,6)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black")) 
ggsave(paste(output_dir,"/ploidy_dist.pdf",sep = ""),
       height = 1.5,
       width = 2.5)

ggtree(sc_tree,layout = "roundrect", size = 0.15) +
  theme_void() 
  ggsave(paste(output_dir,"/small_tree.pdf",sep = ""),
       height = 2.5,
       width = 1)

  
