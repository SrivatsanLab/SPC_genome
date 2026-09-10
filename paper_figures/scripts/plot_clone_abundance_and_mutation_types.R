project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(ggridges)
setwd(file.path(project_root, "paper_figures/data/Clonal_K562/"))
output_dir = file.path(project_root, "paper_figures/output/plot_clone_abundance_and_mutation_types/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
  #theme_white_text_and_axes

ggsave(filename = file.path(output_dir, "population_cell_count_legacy.png"),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.25)

ggsave(filename = paste(output_dir,"population_cell_count.png",sep = ""),
       dpi= 600,bg = "transparent",
       height = 2, width = 2)

setwd(file.path(project_root, "paper_figures/data/Single_cell_bottlenecking_summary_statistics/"))
output_dir = file.path(project_root, "paper_figures/output/plot_clone_abundance_and_mutation_types/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

snps_per_Cell <- read.csv("snp_counts.csv") %>%
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
ggsave(filename = file.path(output_dir, "variants_per_cell.pdf"),
       bg = "transparent",
       height = 3.5, width = 2.25)
  

single_cell_clones_data = read.csv("snp_counts.csv") %>%
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

ggsave(filename = file.path(output_dir, "SBS_per_cell_violin.pdf"),
       bg = "transparent",
       height = 3.5, width = 2.25)

single_cell_clones_data %>%
  select(Cell, SNP_count_de_novo) %>%
  summarise(mean_sbs = mean(SNP_count_de_novo),
            sd_sbs = sqrt(var(SNP_count_de_novo)))


single_cell_clones_data %>%
  pull(Cell) %>% length()
