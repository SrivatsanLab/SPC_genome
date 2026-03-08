project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)


library(ggplot2)
library(dplyr)
library(ggridges)
setwd(file.path(project_root, "paper_figures/data/Images/20250111_SPCs_of_different_Sizes/"))
output_dir = file.path(project_root, "paper_figures/output/plot_spc_sizes/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(pattern = "mask_summaries.*\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No mask_summaries CSV files found in paper_figures/data/Images/20250111_SPCs_of_different_Sizes/")
}
results = read.csv(files[1], header = TRUE) %>% as.data.frame()
if (!("Type" %in% names(results))) {
  file_name_col <- NULL
  if ("File Name" %in% names(results)) {
    file_name_col <- "File Name"
  } else if ("File.Name" %in% names(results)) {
    file_name_col <- "File.Name"
  }
  if (!is.null(file_name_col)) {
    results <- results %>%
      mutate(
        Type = case_when(
          grepl("G4", .data[[file_name_col]], ignore.case = TRUE) ~ "G4",
          grepl("C4", .data[[file_name_col]], ignore.case = TRUE) ~ "C4",
          grepl("C2", .data[[file_name_col]], ignore.case = TRUE) ~ "C2",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Type))
  } else {
    stop("Input CSV must contain either a Type column or a 'File Name' column.")
  }
}

results = 
  results %>%
  filter(Size < 3e4)

cutoffs = tibble(Type = c("G4","C4","C2"),
           Upper = c(25000,10000,2000),
           Lower = c(10000,5000,0))
results = 
  left_join(results, cutoffs)

results = results[results$Size < results$Upper,]
results = results[results$Size > results$Lower,]

results %>%
  ggplot() +
  geom_histogram(aes(x = Size),bins = 100) +
  facet_wrap(~Type,scales = "free_x") 

results %>%
  mutate(diameter_um = sqrt(Size/pi) * 2 * 0.8838) %>%
  group_by(Type) %>%
  summarise(mean_size = mean(diameter_um),
            sd_size = sqrt(var(diameter_um)),
            cv_size = var(diameter_um)/mean_size)
