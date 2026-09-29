project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)


library(ggplot2)
library(dplyr)
library(ggridges)
library(scales)  # for number formatting

rm(list = setdiff(ls(), "project_root"))

setwd(file.path(project_root, "paper_figures/data/Single_cell_barnyard_summary_statistics/"))
output_dir = file.path(project_root, "paper_figures/output/plot_genome_barnyard/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


aligned_reads = read.table(file = "barnyard.csv", sep = ",", header = T)
metadata = read.table(file = "barn10_metadata.csv", sep = ",", header = T) 

head(aligned_reads)
head(metadata)

aligned_reads = left_join(aligned_reads, metadata,by = "barcode")

aligned_reads =
  aligned_reads %>%
  mutate(total_reads = mouse_count + human_count)

aligned_reads %>% 
  arrange(-total_reads) %>%
  mutate(rank = row_number()) %>%
  ggplot() +
  geom_line(aes(x = rank,
                y = total_reads)) +
  scale_x_log10() +
  scale_y_log10()

aligned_reads$species = ifelse((is.na(aligned_reads$species)& aligned_reads$total_reads > 2.5e5), "collision",aligned_reads$species)

aligned_reads %>% 
  filter(!is.na(species)) %>%
  ggplot() +
  geom_point(aes(x = mouse_count,
                 y = human_count, 
                 color = species),
             size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = c(1e6, 2e6, 3e6, 4e6, 5e6),
    labels =c("1e6", "2e6", "3e6", "4e6", "5e6"),
    limits = c(0,5e6)) +
  scale_y_continuous(
    breaks = c(1e6, 2e6, 3e6, 4e6, 5e6),
    labels =c("1e6", "2e6", "3e6", "4e6", "5e6"),
    limits = c(0,5e6)) +
  xlab("Mouse Reads") +
  ylab("Human Reads") +
  scale_color_manual(values = c("#9C60E2","#04BC09","#ff1a5e"))
ggsave(filename = paste(output_dir,"genome_barnyard.pdf"),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.5)
  


aligned_reads %>% 
  filter(!is.na(species)) %>%
  ggplot() +
  geom_point(aes(x = mouse_count,
                 y = human_count, 
                 color = species)) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = c( 2e5, 4e5, 6e5, 8e5, 1e6),
    labels =c("2e5", "4e5", "6e5", "8e5", "1e6"),
    limits = c(0,1e6)) +
  scale_y_continuous(
    breaks = c( 2e5, 4e5, 6e5, 8e5, 1e6),
    labels =c("2e5", "4e5", "6e5", "8e5", "1e6"),
    limits = c(0,1e6)) +
  xlab("Mouse Reads") +
  ylab("Human Reads") +
  scale_color_manual(values = c("#9C60E2","#04BC09","#ff1a5e"))
ggsave(filename = paste(output_dir,"genome_barnyard_zoom.pdf"),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.5)


# Plot Gini Coefficients --------------------------------------------------

metadata %>%
  ggplot() +
  geom_boxplot(aes(x = amplification,
                   y = gini,
                   fill = species)) + 
  ylim(0,1) +
  scale_fill_manual(values = c("#04BC09","#ff1a5e")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  xlab("WGA Method") +
  ylab("Gini Coefficient")

ggsave(filename = paste(output_dir,"gini_per_cell.png"),
       dpi= 600,bg = "transparent",
       height = 4, width = 2)
