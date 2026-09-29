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

snps_per_Cell %>%
  dplyr::select(-X) %>%
  pivot_longer(-pop) %>%
  filter(name == "AD_2") %>%
  ggplot + 
  geom_boxplot(aes(x = pop,
                   y = value,
                   fill = pop),
               outlier.stroke = 0,
               outlier.size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("#a6d8a5","#f9dda7","#e3a9a6","#deaacd","#c2bbd9","#a8c7e6")) +
  scale_y_log10() +
  xlab("Population") +
  scale_x_discrete(labels = c("P0","P1","P2","P3","P4","P5")) +
  ylab("Variants") +
  coord_flip() +
  theme(legend.position = "none")
ggsave(filename = paste(output_dir,"mutations_per_population.png",sep = ""),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 1.75)
  

single_cell_clones_data= read.csv("snp_counts.csv")

single_cell_clones =
  single_cell_clones_data %>%
  mutate(shared_SNPs =  SNP_count_de_novo - SNP_count_de_novo_singletons) %>%
  dplyr::select(-SNP_count_de_novo)
  
single_cell_clones = 
single_cell_clones %>%
  dplyr::rename(Cell = X) %>% 
  dplyr::select(-tot_observed_sites,-SNP_count) %>%
  tidyr::pivot_longer(-c(Cell,pop),names_to = "mutation")

single_cell_clones$mutation <- 
  factor(single_cell_clones$mutation, 
         levels = c("germline_count", "mito_snps", "SNP_count_de_novo_singletons","shared_SNPs"))

pop_levels = rev(paste0("P_", 0:5))

single_cell_clones %>%
  mutate(pop = factor(pop, levels = pop_levels)) %>%
ggplot() +
  geom_violin(aes(x = 1,
                   y = value,
                   fill= pop),
              color = "black",
              size = 0.25) +
  scale_fill_manual(values = c("P_0" = "#260091", "P_1" = "#1e90ff",
                               "P_2" = "#ffdb58", "P_3" = "#ff9d71",
                               "P_4" = "#ff1b5e", "P_5" = "#e3e6e6")) +
  facet_wrap(~mutation,scales = "free_y", nrow = 1)+
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

ggsave(filename = paste(output_dir,"mutations_per_population.png",sep = ""),
       dpi= 600,bg = "transparent",
       height = 2, width = 4)

single_cell_clones_data %>% 
  rename(Cell = X)%>%
  select(Cell, SNP_count_de_novo) %>%
  summarise(mean_sbs = mean(SNP_count_de_novo),
            sd_sbs = sqrt(var(SNP_count_de_novo)))


single_cell_clones_data %>%
  pull(X) %>% length()
