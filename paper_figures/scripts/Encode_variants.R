project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)



library(ggplot2)
library(dplyr)
library(ggridges)
rm(list = setdiff(ls(), "project_root"))
setwd(file.path(project_root, "paper_figures/data/encode_variant_annotation/"))
output_dir = file.path(project_root, "paper_figures/output/encode_variants/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required_files <- c(
  "germline_max_peaks.tsv",
  "somatic_max_peaks.tsv",
  "random_max_peaks.tsv",
  "k562_variants_HMM18.summary",
  "sc_PolE_variants_HMM18.summary",
  "random_variants_HMM18.summary"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    sprintf(
      "Missing required input files in encode_variant_annotation: %s",
      paste(missing_files, collapse = ", ")
    )
  )
}

germline_peaks = 
  read.csv("germline_max_peaks.tsv",
           sep = "\t" ,
           header = T) %>%
  data.frame() %>%
  dplyr::select(peak_bin) %>%
  mutate(type = "Germline")

somatic_peaks = 
  read.csv("somatic_max_peaks.tsv",
           sep = "\t" ,
           header = T) %>%
  data.frame() %>%
  dplyr::select(peak_bin) %>%
  mutate(type = "Somatic")

random_peaks = 
  read.csv("random_max_peaks.tsv",
           sep = "\t" ,
           header = T) %>%
  data.frame() %>%
  dplyr::select(peak_bin) %>%
  mutate(type = "Random Sample")

all_peaks = 
  rbind(germline_peaks, 
        somatic_peaks, 
        random_peaks)

all_peaks = 
  all_peaks %>%
  group_by(type,peak_bin) %>%
  summarise(n = n())

all_peaks$peak_bin <- 
  factor(all_peaks$peak_bin, 
         levels = c("G1", "S1", "S2", "S3", "S4", "G2"))

all_peaks %>%
  group_by(type) %>%
  mutate(total = sum(n)) %>%
  ggplot() +
  geom_bar(aes(x = peak_bin,
               y = 100*n/total,
               fill = type),
           stat = "identity",
           position = "dodge",
           color = "black") +
  theme_classic() +
  ylab("Percent Total") +
  xlab("Cell Cycle Phase") +
  scale_y_continuous(breaks = c(0,5,10,15,20), labels = c("0","5%","10%","15%","20%")) +
  scale_fill_brewer(palette = 6) +
  theme(#legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))

ggsave(paste(output_dir,"cellcycle_phase.png"),
       height = 3, width = 5)

germline_HMM18 = 
  read.csv("k562_variants_HMM18.summary",
           sep = "\t" ,
           header = F) %>%
  data.frame() %>%
  mutate(type = "Germline")


somatic_variants_HMM18 = 
  read.csv("sc_PolE_variants_HMM18.summary",
           sep = "\t" ,
           header = F) %>%
  data.frame() %>%
  mutate(type = "Somatic")

random_variants_HMM18 = 
  read.csv("random_variants_HMM18.summary",
           sep = "\t" ,
           header = F) %>%
  data.frame() %>%
  mutate(type = "Random Sample")

all_annotations = 
  rbind(germline_HMM18, 
        somatic_variants_HMM18, 
        random_variants_HMM18)
colnames(all_annotations) = c("count","annotation","type")

all_peaks$peak_bin <- 
  factor(all_peaks$peak_bin, 
         levels = c("G1", "S1", "S2", "S3", "S4", "G2"))



all_annotations %>%
  group_by(type) %>%
  mutate(total = sum(count)) %>%
  ggplot() +
  geom_bar(aes(x = annotation,
               y = 100*count/total,
               fill = type),
           stat = "identity",
           position = "dodge",
           color = "black") +
  theme_classic() +
  ylab("Percent Total") +
  xlab("HMM18 Encode Annotation") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60), labels = c("0","10%","20%","30%","40%","50%","60%")) +
  scale_fill_brewer(palette = 6) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.margin = unit(c(1, 1, 1, 1), "lines"))

ggsave(paste(output_dir,"encode_annotations.png"),
       height = 3, width = 5)


all_annotations %>%
  group_by(type) %>%
  mutate(percent_total = count/sum(count)) %>%
  filter(annotation == "TxWk")

