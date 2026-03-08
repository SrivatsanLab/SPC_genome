project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)


library(ggplot2)
library(dplyr)
setwd(file.path(project_root, "paper_figures/data/retention_oligomers/"))
output_dir = file.path(project_root, "paper_figures/output/plot_segmented_SPC_oligomers/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files = list.files()[grep(pattern = "mask_summaries_",list.files())]
files

results <- list()  # Initialize an empty list
for (file in files) {
  a <- read.csv(file, header = TRUE) %>% as.data.frame()  # Read and convert to a data frame
  print(dim(a))  # Print the dimensions of the current file's data
  results <- append(results, list(a))  # Append the data frame to the list
}

# Print the length of the results list
print(length(results))

# Combine all data frames in the list and print the dimensions of the resulting data frame
combined_results <- do.call(rbind, results)

combined_results$Time = ifelse(combined_results$Time == "hour1",60,combined_results$Time)
combined_results$Time = ifelse(combined_results$Time == "Day1",1440,combined_results$Time)
combined_results$Time = ifelse(combined_results$Time == "pre",1,combined_results$Time)
combined_results$Time = as.numeric(combined_results$Time)
combined_results$Time %>% unique()

# Visualize the results to get rid of any mis-segmented objects
ggplot(combined_results) +
  geom_point(aes(x = Size, 
                 y = Mask.Sum,
                 color = Design)) +
  scale_y_log10()

# Remove the very large segmented objects in each image
combined_results = 
  combined_results %>%
  filter(Size < 2e6) 

combined_results_mean_0 = 
  combined_results %>%
  filter(Time == 1) %>%
  group_by(Design) %>%
  summarise(mean_time_0 = mean(Mask.Sum))
  
combined_results %>%
  left_join(combined_results_mean_0) %>%
  mutate(Mask.Sum = Mask.Sum / mean_time_0) %>%
  group_by(File.Name,Design,Time) %>%
  mutate(mean_fu = mean(Mask.Sum)) %>%
  ggplot() +
  geom_smooth(aes(x = Time ,
                y = Mask.Sum,
                color = Design),method = "lm") +
  geom_point(aes(x = Time,
                 y = mean_fu,
                 color = Design)) +
  theme_classic() +
  scale_x_continuous(breaks = c(60, 1440, 2880, 4320, 7200,14400),
                     labels = c("1hr","day 1", "day 2", "day 3", "day 5","day 10")) +
  scale_color_manual(values = c("#ff1a5e","#438cfd")) +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10)) +
  scale_y_continuous(limits = c(0, 1.3)) +
  xlab("Time") +
  ylab("Flourescence per SPC (AU)") 
ggsave(filename = paste(output_dir,"line_plot.pdf"),
       height = 3, width = 6)

# Replot the data
ggplot(combined_results) +
  geom_point(aes(x = Size, 
                 y = Mask.Sum,
                 color = as.factor(Time))) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~Design)

ggplot(combined_results) +
  geom_boxplot(aes(x = as.factor(Design), 
                 y = Mask.Sum,
                 fill = as.factor(Time)),
               #outlier.stroke = 0,
               #outlier.size = 0.25,
               linewidth = 0.5,
               color = "black") +
  #ylim(3e4,5e5) +
  theme_classic() +
  xlab("Oligomer") +
  ylab("Cumulative Flourescence") +
  scale_fill_viridis_d() +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10))
ggsave(filename = paste(output_dir,"boxplot.pdf"),
       height = 3, width = 3)

