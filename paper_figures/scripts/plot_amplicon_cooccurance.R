project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(ggridges)

setwd(file.path(project_root, "paper_figures/data/Amplicons/"))
output_dir = file.path(project_root, "paper_figures/output/plot_amplicon_cooccurance/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_white_text_and_axes <- theme(
  text = element_text(color = "white"),                # All text elements white
  axis.text = element_text(color = "white", size = 8),          # Axis text white
  axis.title = element_text(color = "white", size = 10),         # Axis titles white
  axis.ticks = element_line(color = "white"),         # Axis ticks white
  axis.line = element_line(color = "white"),          # Axis lines white
  plot.background = element_blank(), #Black plot background
  panel.grid = element_blank(),
  legend.position = "none",
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.ontop = F,
  strip.background = element_blank(),
  strip.text = element_blank(),# Optional: Remove gridlines
)

files = list.files()[grep(pattern = "counts",list.files())]
results <- list()  # Initialize an empty list
files
for (file in files){
  a <- read.csv(file, header = TRUE) %>% as.data.frame()  # Read and convert to a data frame
  a$name = file
  print(dim(a))  # Print the dimensions of the current file's data
  results <- append(results, list(a))  # Append the data frame to the list
}

# Print the length of the results list
print(length(results))
results = 
  do.call(rbind, results) %>% 
  as.data.frame()

results$type = ifelse(grepl(x = results$name,"non"),"Intact","Fragmented" )

results$mix = ifelse(grepl(x = results$name,"90"),"10:90",ifelse(grepl(x = results$name,"80"),"80:20","100:0"))
results$mix <- factor(results$mix, levels = c("10:90", "80:20", "100:0"))


ggplot() +
  geom_point(data = results %>%
               dplyr::select(-mix) %>%
               filter(type == "Intact"),
             aes(x = Ref_Counts,
                 y = Alt_Counts),
             color = "grey90",
             size = 1) +
  theme_classic() +
  xlab("Reads Allele C (PBMC)") +
  ylab("Reads Allele T (K562)") +
  geom_point(data = results %>%
               filter(type == "Intact"),
             aes(x = Ref_Counts,
                 y = Alt_Counts,
                 color = mix),
             size = 1) +
  facet_wrap(~mix) +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10),
        strip.text = element_blank()) +
  scale_color_manual(values = c("black","black","black"))
ggsave(filename = paste(output_dir,"barnyard_amplicon.png"),
       dpi= 600,bg = "transparent",
       height = 2, width = 4.5)


ggplot() +
  geom_point(data = results %>%
               dplyr::select(-mix) %>%
               filter(type == "Intact"),
             aes(x = Ref_Counts,
                 y = Alt_Counts),
             color = "grey30",
             size = 1) +
  theme_classic() +
  xlab("Reads Allele C (PBMC)") +
  ylab("Reads Allele T (K562)") +
  geom_point(data = results %>%
               filter(type == "Intact"),
             aes(x = Ref_Counts,
                 y = Alt_Counts,
                 color = mix),
             size = 1) +
  facet_wrap(~mix) +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10),
        strip.text = element_blank()) +
  theme_white_text_and_axes +
  scale_color_manual(values = c("white","white","white"))
ggsave(filename = paste(output_dir,"barnyard_amplicon_white.png"),
       dpi= 600,bg = "transparent",
       height = 2, width = 4.5)

             
results %>%
  filter(mix != "100:0") %>%
  mutate(percent_T = Alt_Counts/(Alt_Counts + Ref_Counts)) %>%
  group_by(type) %>%
  filter(percent_T > 0.5) %>%
  summarise(mean_percent_alt = mean(percent_T),
            count = n(),
            mean_alt = mean(Alt_Counts),
            std_alt = sqrt(var(percent_T)))


  

# Het Allele ------------------------------------------------------------------

setwd(file.path(project_root, "paper_figures/data/Amplicons/Het_Allele/"))
output_dir = file.path(project_root, "paper_figures/output/plot_amplicon_cooccurance/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


files = list.files()[grep(pattern = "counts",list.files())]
results <- list()  # Initialize an empty list
files
for (file in files){
  a <- read.csv(file, header = TRUE) %>% as.data.frame()  # Read and convert to a data frame
  a$name = file
  print(dim(a))  # Print the dimensions of the current file's data
  results <- append(results, list(a))  # Append the data frame to the list
}

# Print the length of the results list
print(length(results))
results = 
  do.call(rbind, results) %>% 
  as.data.frame()
results


ggplot() +
  xlab("Reads Allele C (Ref)") +
  ylab("Reads Allele T (Alt)") +
  geom_point(data = results,
             aes(x = Reference_Counts,
                 y = Alternate_Counts),
             size = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10),
        strip.text = element_blank()) 
ggsave(filename = paste(output_dir,"het_amplicon.png"),
       dpi= 600,bg = "transparent",
       height = 2, width = 2)

ggplot() +
  xlab("Reads Allele C (Ref)") +
  ylab("Reads Allele T (Alt)") +
  geom_point(data = results,
             aes(x = Reference_Counts,
                 y = Alternate_Counts),
             color = "white",
             size = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10),
        strip.text = element_blank()) +
  theme_white_text_and_axes 
ggsave(filename = paste(output_dir,"het_amplicon_white.png"),
       dpi= 600,bg = "transparent",
       height = 2, width = 2)



# Percent detection of amplicons ------------------------------------------


setwd(file.path(project_root, "paper_figures/data/Amplicons/Fragmented/"))
output_dir = file.path(project_root, "paper_figures/output/plot_amplicon_cooccurance/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files = list.files()[grep(pattern = "pset",list.files())]
results <- list()  # Initialize an empty list
files
for (file in files){
  a <- read.csv(file, header = TRUE) %>% as.data.frame()  # Read and convert to a data frame
  a$name = file
  print(dim(a))  # Print the dimensions of the current file's data
  results <- append(results, list(a))  # Append the data frame to the list
}

  
# Print the length of the results list
print(length(results))
results_fragmented = 
  do.call(rbind, results) %>% 
  as.data.frame()

results_fragmented = 
  results_fragmented %>%
  group_by(AKT3, FANCC, BRCA1, BIRC6) %>%
  summarise(value = sum(value)) %>%
  mutate(type = "Fragmented")



setwd(file.path(project_root, "paper_figures/data/Amplicons/FullLength/"))
output_dir = file.path(project_root, "paper_figures/output/plot_amplicon_cooccurance/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
rm(results)
rm(files)

files = list.files()[grep(pattern = "pset",list.files())]
results <- list()  # Initialize an empty list
files
for (file in files){
  a <- read.csv(file, header = TRUE) %>% as.data.frame()  # Read and convert to a data frame
  a$name = file
  print(dim(a))  # Print the dimensions of the current file's data
  results <- append(results, list(a))  # Append the data frame to the list
}

results
# Print the length of the results list
print(length(results))
results_full_length = 
  do.call(rbind, results) %>% 
  as.data.frame()

results_full_length = 
  results_full_length %>%
  group_by(AKT3, FANCC, BRCA1, BIRC6) %>%
  summarise(value = sum(value)) %>%
  mutate(type = "Intact")

results = rbind(results_fragmented,results_full_length)
results

results$AKT3 = as.logical(results$AKT3)
results$FANCC = as.logical(results$FANCC)
results$BRCA1 = as.logical(results$BRCA1)
results$BIRC6 = as.logical(results$BIRC6)

head(results)

results %>%
  ungroup() %>%
  mutate(total_detected = AKT3 +FANCC + BRCA1 + BIRC6) %>%
  group_by(total_detected,type) %>%
  summarise(value = sum(value)) %>% 
  group_by(type) %>%
  mutate(percent = 100 * value / sum(value)) %>% 
  ggplot() +
  geom_bar(aes(x = total_detected,
               y = percent),
           color = "black",
           stat = "identity",
           fill = "grey80",
           linewidth = 0.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10),
        strip.text = element_blank()) +
  scale_fill_viridis_d(option = "A") +
  facet_wrap(~type) +
  ylab("Percentage of Cells") +
  xlab("Amplicons Detected Per Cell")
ggsave(filename = paste(output_dir,"barnyard_amplicon_percent_detected.png"),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.75)

  

# Percent detection of FL amplicons ------------------------------------------


setwd(file.path(project_root, "paper_figures/data/Amplicons/FullLength/"))
output_dir = file.path(project_root, "paper_figures/output/plot_amplicon_cooccurance/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files = list.files()[grep(pattern = "pset",list.files())]
results <- list()  # Initialize an empty list
files
for (file in files){
  a <- read.csv(file, header = TRUE) %>% as.data.frame()  # Read and convert to a data frame
  a$name = file
  print(dim(a))  # Print the dimensions of the current file's data
  results <- append(results, list(a))  # Append the data frame to the list
}


# Print the length of the results list
print(length(results))
results = 
  do.call(rbind, results) %>% 
  as.data.frame()

results = 
  results %>%
  group_by(AKT3, FANCC, BRCA1, BIRC6) %>%
  summarise(value = sum(value)) 

results$AKT3 = as.logical(results$AKT3)
results$FANCC = as.logical(results$FANCC)
results$BRCA1 = as.logical(results$BRCA1)
results$BIRC6 = as.logical(results$BIRC6)

head(results)

results %>%
  ungroup() %>%
  mutate(total_detected = AKT3 +FANCC + BRCA1 + BIRC6) %>%
  group_by(total_detected) %>%
  summarise(value = sum(value)) %>%
  ungroup() %>%
  mutate(percent = 100 * value / sum(value),
         type = "Fragmented") %>%
  ggplot() +
  geom_bar(aes(x = type,
               y = percent,
               fill = as.character(total_detected)),
           color = "black",
           stat = "identity") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = "none",
        axis.text = element_text(size = 10),
        strip.text = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_viridis_d(option = "G") +
  ylab("Percent Detected")
ggsave(filename = paste(output_dir,"barnyard_amplicon_percent_detected.png"),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 1.5)
