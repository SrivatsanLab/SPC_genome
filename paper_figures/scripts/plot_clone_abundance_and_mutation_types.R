project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
output_dir = file.path(project_root, "paper_figures/output/plot_clone_abundance_and_mutation_types/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# clonal_cell_count.csv is not in this checkout; skip rather than halt the
# whole script, as bulk_VAF_bottleneck.R does for its missing input
cell_count_path <- file.path(project_root,
                             "paper_figures/data/Clonal_K562/clonal_cell_count.csv")
if (file.exists(cell_count_path)) {
  setwd(file.path(project_root, "paper_figures/data/Clonal_K562/"))

  cell_counts = read.csv(file = "clonal_cell_count.csv")

  cell_counts$Population <- factor(cell_counts$Population, levels = c("P5", "P4", "P3", "P2", "P1", "P0"))

  cell_counts %>%
    ggplot() +
    geom_point(aes(x = Population,
                   y = Cell.Count,
                   fill = Population),
               size = 2,
               shape = 21,
               stroke = 0.5) +
    scale_fill_manual(values = c("P0" = "#260091", "P1" = "#1e90ff",
                                 "P2" = "#ffdb58", "P3" = "#ff9d71",
                                 "P4" = "#ff1b5e", "P5" = "#e3e6e6")) +
    theme_classic() +
    theme(legend.position = "none",
          title = element_blank()) +
    ylim(0,1250000)+
    ylab("Cell Count")

  ggsave(filename = file.path(output_dir, "population_cell_count_legacy.svg"),
         bg = "transparent",
         height = 2.5, width = 2.25)

  ggsave(filename = file.path(output_dir, "population_cell_count.svg"),
         bg = "transparent",
         height = 2, width = 2)
} else {
  message("Skipping population cell count: ", cell_count_path, " not found")
}

setwd(file.path(project_root, "paper_figures/data/Single_cell_bottlenecking_summary_statistics/"))
output_dir = file.path(project_root, "paper_figures/output/plot_clone_abundance_and_mutation_types/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

snps_per_Cell <- read.csv("filt_ns_snp_count.csv") %>%
  mutate(AD_2 = SNP_count_de_novo)

# same P0-P5 palette as the population_cell_count panel above. Names are bound
# to the sorted pop levels, so a mismatch fails loudly instead of mis-colouring.
pop_levels <- sort(unique(as.character(snps_per_Cell$pop)))
stopifnot(length(pop_levels) == 6)
pop_palette <- setNames(
  c("#260091", "#1e90ff", "#ffdb58", "#ff9d71", "#ff1b5e", "#e3e6e6"),
  pop_levels
)

snps_per_Cell %>%
  dplyr::select(-X) %>%
  tidyr::pivot_longer(-pop) %>%
  filter(name == "AD_2") %>%
  # ascending order so the vertical panel runs P0 (left) -> P5 (right)
  mutate(pop = factor(pop, levels = pop_levels)) %>%
  ggplot + 
  geom_boxplot(aes(x = pop,
                   y = value,
                   fill = pop),
               outlier.stroke = 0,
               outlier.size = 1) +
  theme_classic() +
  scale_fill_manual(values = pop_palette) +
  scale_y_log10() +
  xlab("Population") +
  scale_x_discrete(labels = c("P0","P1","P2","P3","P4","P5")) +
  ylab("Variants") +
  theme(legend.position = "none")
ggsave(filename = file.path(output_dir, "variants_per_cell.svg"),
       bg = "transparent",
       height = 3.5, width = 2.25)

# Original Figure 4 panel: the same per-cell counts as a horizontal boxplot in
# the pastel palette (P0 at the bottom). The legacy input was SNP_counts.csv,
# plotting its AD_2 column; SNP_count_de_novo is the reprocessed equivalent.
snps_per_Cell %>%
  mutate(pop = factor(pop, levels = pop_levels)) %>%
  ggplot() +
  geom_boxplot(aes(x = pop,
                   y = AD_2,
                   fill = pop),
               outlier.stroke = 0,
               outlier.size = 1) +
  theme_classic() +
  scale_fill_manual(values = setNames(
    c("#a6d8a5", "#f9dda7", "#e3a9a6", "#deaacd", "#c2bbd9", "#a8c7e6"),
    pop_levels
  )) +
  scale_y_log10() +
  xlab("Population") +
  scale_x_discrete(labels = c("P0","P1","P2","P3","P4","P5")) +
  ylab("Variants") +
  coord_flip() +
  theme(legend.position = "none")
ggsave(filename = file.path(output_dir, "mutations_per_population.svg"),
       bg = "transparent",
       height = 2.5, width = 1.75)


single_cell_clones_data = read.csv("filt_ns_snp_count.csv") %>%
  dplyr::rename(Cell = X)

# SBS per cell as a violin, styled to match the boxplot above: same fills,
# ascending P0 -> P5, theme_classic, log10 y.
single_cell_clones_data %>%
  mutate(pop = factor(pop, levels = pop_levels)) %>%
  ggplot() +
  geom_violin(aes(x = pop,
                  y = SNP_count_de_novo,
                  fill = pop),
              color = "black",
              linewidth = 0.25) +
  scale_fill_manual(values = pop_palette) +
  theme_classic() +
  scale_y_log10() +
  scale_x_discrete(labels = c("P0","P1","P2","P3","P4","P5")) +
  xlab("Population") +
  ylab("Variants") +
  theme(legend.position = "none")

ggsave(filename = file.path(output_dir, "SBS_per_cell_violin.svg"),
       bg = "transparent",
       height = 2, width = 2)

# Original Figure 4 violin: one facet per mutation category. It saved to
# mutations_per_population.png too, overwriting the boxplot, so it gets its own
# name here. The reprocessed table has none of the category columns, so skip
# rather than halt, as with the cell count panel.
category_cols <- c("germline_count", "mito_snps", "SNP_count_de_novo_singletons")
missing_cols <- setdiff(category_cols, colnames(single_cell_clones_data))
if (length(missing_cols) == 0) {
  single_cell_clones_data %>%
    mutate(shared_SNPs = SNP_count_de_novo - SNP_count_de_novo_singletons) %>%
    dplyr::select(Cell, pop, all_of(category_cols), shared_SNPs) %>%
    tidyr::pivot_longer(-c(Cell, pop), names_to = "mutation") %>%
    mutate(mutation = factor(mutation, levels = c(category_cols, "shared_SNPs")),
           # descending so P5 is drawn first and P0 sits on top
           pop = factor(pop, levels = rev(pop_levels))) %>%
    ggplot() +
    geom_violin(aes(x = 1,
                    y = value,
                    fill = pop),
                color = "black",
                linewidth = 0.25) +
    scale_fill_manual(values = pop_palette) +
    facet_wrap(~mutation, scales = "free_y", nrow = 1) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"),
          strip.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 6))

  ggsave(filename = file.path(output_dir, "mutation_categories_per_population.svg"),
         bg = "transparent",
         height = 2, width = 4)
} else {
  message("Skipping mutation category violin: filt_ns_snp_count.csv lacks ",
          paste(missing_cols, collapse = ", "))
}

single_cell_clones_data %>%
  select(Cell, SNP_count_de_novo) %>%
  summarise(mean_sbs = mean(SNP_count_de_novo),
            sd_sbs = sqrt(var(SNP_count_de_novo)))


single_cell_clones_data %>%
  pull(Cell) %>% length()
