project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)


library(ggplot2)
library(dplyr)
library(tidyr)
setwd(file.path(project_root, "paper_figures/data/Single_cell_barnyard_summary_statistics/"))
output_dir = file.path(project_root, "paper_figures/output/plot_gini_lorenz/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


data = read.csv("lorenz.csv")

data %>%
  pivot_longer(-X,names_to = "WGA") %>%
  ggplot() +
  geom_abline(linetype = "dotdash") +
  geom_line(aes(x = X,
                y = value,
                color = WGA),
            size = 1.1) +
  theme_classic() +
  theme(#legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  xlab("Cumulative Fraction of Genome") +
  ylab("Cumulative Fraction of Reads") +
  scale_color_manual(values = c("red","blue"))

ggsave(filename = paste(output_dir,"lorenz_curve.png",sep = ""),
       dpi= 600,bg = "transparent",
       height = 3, width = 5)

# NEed to 