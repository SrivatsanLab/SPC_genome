project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)


library(flowCore)
library(ggplot2)
library(dplyr)
library(ggridges)
setwd(file.path(project_root, "paper_figures/data/FlowData/"))

# Define the path to your folder with FCS files
folder_path <- "."
output_dir <- file.path(project_root, "paper_figures/output/analyse_flowdata/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# List all .fcs files in the folder
fcs_files <- list.files(path = folder_path, pattern = "\\.fcs$", full.names = TRUE)

# Read all FCS files into a flowSet
fcs_data <- read.flowSet(files = fcs_files, transformation = FALSE)
length(fcs_data)

# Get sample (file) names from the flowSet
file_names <- sampleNames(fcs_data)

parse_metadata <- function(filename) {
  # Determine treatment
  treatment <- case_when(
    str_detect(filename, "nocodazole_treated") ~ "Nocodazole",
    str_detect(filename, "taxol_treated") ~ "Taxol",
    TRUE ~ "No Treatment"
  )
  
  # Remove prefix and suffix to extract sample name
  sample <- filename %>%
    str_replace("^(nocodazole_treated_|taxol_treated_|K562 Cells_)", "") %>%
    str_replace("_\\d+.*\\.fcs$", "") %>%    # Remove timestamp + .fcs
    str_replace("\\.fcs$", "")               # For files without timestamps
  
  sample <- str_replace_all(sample, " ", "_")  # Standardize spacing with underscores
  
  return(list(Treatment = treatment, Sample = sample))
}

# Build combined data frame with metadata
plot_df <- lapply(seq_along(fcs_data), function(i) {
  parsed <- parse_metadata(file_names[i])
  exprs(fcs_data[[i]]) %>%
    as.data.frame() %>%
    mutate(Sample = parsed$Sample,
           Treatment = parsed$Treatment)
}) %>%
  bind_rows()
# Replace all the - with _
colnames(plot_df) <- gsub("-", "_", colnames(plot_df))

# Make a column called supersample for whether the cells have a WT polE or mutatn
plot_df <- 
  plot_df %>%
  mutate(
    supersample = case_when(
      grepl("AAVS", Sample, ignore.case = TRUE) ~ "AAVS",
      grepl("PolE", Sample, ignore.case = TRUE) ~ "PolE",
      TRUE ~ "Other"
    )
  )

colnames(plot_df) <- gsub("-", "_", colnames(plot_df))
plot_df$Sample<- gsub("_", " ", plot_df$Sample)

plot_df[plot_df$DAPI_A < 262144,] %>%
  as.data.frame() %>%
  dplyr::filter(Treatment == "No Treatment") %>%
  ggplot() +
  geom_density_ridges(aes(x = DAPI_H, y = Sample, fill = supersample),stat = "binline",bins = 250) +
  #geom_vline(xintercept = 27500, color = "black", linetype = "dashed") +
  theme_minimal() +
  labs(x = "DAPI-A") +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_brewer(palette = "Set1") 
  #scale_x_continuous(breaks = c(100000,200000)) 
ggsave(paste(output_dir,"no_treatment_dapi_A.pdf",sep =""),
       height = 2,
       width = 4)


plot_df$Sample %>%
  unique()

plot_df[plot_df$DAPI_A > 10000,] %>%
ggplot() + 
  geom_point(aes(x = (FSC_A),
                 y = log10(SSC_A))) +
  facet_wrap(Sample ~ Treatment, scales = "free_y",ncol = 3)

dispersion_stats <- plot_df %>%
  group_by(Sample, Treatment) %>%
  summarise(
    median_DAPI = median(DAPI_A, na.rm = TRUE),
    mad_DAPI = mad(DAPI_A, na.rm = TRUE),
    iqr_DAPI = IQR(DAPI_A, na.rm = TRUE),
    sd_DAPI = sd(DAPI_A, na.rm = TRUE),
    .groups = "drop"
  )

print(dispersion_stats)

dispersion_stats <- dispersion_stats %>%
  mutate(
    supersample = case_when(
      grepl("AAVS", Sample, ignore.case = TRUE) ~ "AAVS",
      grepl("POLE", Sample, ignore.case = TRUE) ~ "POLE",
      TRUE ~ "Other"
    )
  )

ggplot(dispersion_stats, aes(x = 1, y = sd_DAPI, fill = supersample)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(group = supersample),
             position = position_dodge(width = 0.75),
             size = 2, shape = 21, stroke = 0.2, alpha = 0.9) +
  theme_classic() +
  facet_wrap(~Treatment,scales = "free_y") +
  labs(
    x = "Treatment",
    y = "SD (DAPI-A)"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(paste(output_dir,"treatment_spread.pdf",sep =""),
       height = 2,
       width = 4)

