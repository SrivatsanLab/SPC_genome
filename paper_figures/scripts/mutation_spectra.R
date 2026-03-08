project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggtree)
library(tidyverse)
library(Biostrings)  # for reverseComplement
set.seed(50)
cosmic_path <- file.path(project_root, "paper_figures/data/external/COSMIC_v3.4_SBS_GRCh38.txt")
if (!file.exists(cosmic_path)) {
  stop(
    paste0(
      "Missing COSMIC signatures file: ", cosmic_path,
      ". Place COSMIC_v3.4_SBS_GRCh38.txt in paper_figures/data/external/."
    )
  )
}


setwd(file.path(project_root, "paper_figures/data/K562_bulk_WGS/"))
output_dir = file.path(project_root, "paper_figures/output/mutation_spectra/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
spectra_diff <- read.csv("spectra_diff.csv", col.names = c("mutation", "value"))
spectra_diff
spectra_diff$density = spectra_diff$value/sum(spectra_diff$value)

# Extract bases
spectra_diff <- spectra_diff %>%
  mutate(
    ref_base = substr(sub("^(...).(...).*$", "\\1", mutation), 2, 2),
    alt_base = substr(sub("^(...).(...).*$", "\\2", mutation), 2, 2),
    base_change = paste0(ref_base, ">", alt_base)
  )

# Reverse complement function for trinucleotides
revcomp_trinuc <- function(trinuc) {
  comp <- c(A="T", T="A", C="G", G="C")
  paste0(rev(comp[strsplit(trinuc, "")[[1]]]), collapse = "")
}

# Convert TTT>TGT → T[T>G]T
to_cosmic_format <- function(tri_mut) {
  ref <- substr(tri_mut, 1, 3)
  alt <- substr(tri_mut, 5, 7)
  pre  <- substr(ref, 1, 1)
  mid  <- substr(ref, 2, 2)
  post <- substr(ref, 3, 3)
  mid_alt <- substr(alt, 2, 2)
  
  # If middle base is purine, reverse complement
  if (mid %in% c("A", "G")) {
    ref_rc <- revcomp_trinuc(ref)
    alt_rc <- revcomp_trinuc(alt)
    pre  <- substr(ref_rc, 1, 1)
    mid  <- substr(ref_rc, 2, 2)
    post <- substr(ref_rc, 3, 3)
    mid_alt <- substr(alt_rc, 2, 2)
  }
  
  paste0(pre, "[", mid, ">", mid_alt, "]", post)
}

# COSMIC 96-order definition
bases <- c("A", "C", "G", "T")
subs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
cosmic_order <- unlist(lapply(subs, function(sub) {
  ref <- substr(sub, 1, 1)
  alt <- substr(sub, 3, 3)
  sapply(bases, function(pre) {
    sapply(bases, function(post) {
      paste0(pre, "[", ref, ">", alt, "]", post)
    })
  })
}))
cosmic_order <- as.character(cosmic_order)

# Apply all transformations
spectra_diff <- spectra_diff %>%
  mutate(
    # Normalize base change to pyrimidine context
    base_change_new = case_when(
      base_change == "A>C" ~ "T>G",
      base_change == "A>G" ~ "T>C",
      base_change == "A>T" ~ "T>A",
      base_change == "G>T" ~ "C>A",
      base_change == "G>C" ~ "C>G",
      base_change == "G>A" ~ "C>T",
      TRUE ~ base_change
    ),
    # Convert to COSMIC-style
    mutation_cosmic = sapply(paste0(mutation, ">", substr(alt_base, 1, 1)), to_cosmic_format),
    # Set factor order
    mutation_cosmic = factor(mutation_cosmic, levels = cosmic_order)
  )
spectra_diff


mutation_colors <- c(
  "C>A" = "#00AEEF",
  "C>G" = "#000000",
  "C>T" = "#EE2E2F",
  "T>A" = "#BFBFBF",
  "T>C" = "#92D050",
  "T>G" = "#E9C3C3"
)

# Step 2: Plot
ggplot(spectra_diff) +
  geom_bar(
    aes(x = mutation_cosmic, y = density, fill = base_change_new),
    stat = "identity",
    color = "black",
    size = 0.25
  ) +
  scale_fill_manual(values = mutation_colors, drop = FALSE) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size = 6,hjust = 0,vjust = 0.15))
ggsave(paste(output_dir,"spectra.pdf",sep = ""),
       height = 4,
       width = 12)


# Decomposition -----------------------------------------------------------
library(deconstructSigs)

m = read.csv(cosmic_path,sep = "\t") %>%
  as_tibble() %>%
  column_to_rownames(var = "Type")

cosmic_raw <- read.csv(cosmic_path, sep = "\t", check.names = FALSE)

# First column is trinucleotide context
contexts <- cosmic_raw[, 1]

# Remaining columns are the signature profiles
cosmic_mat <- as.data.frame(t(cosmic_raw[, -1]))
colnames(cosmic_mat) <- contexts
# signatures_SBS = rownames(cosmic_mat)
cosmic_mat <- as.data.frame(cosmic_mat)

# Step 1: Create a named vector with trinucleotide contexts
signature_vector <- setNames(spectra_diff$density, spectra_diff$mutation_cosmic)

# Step 2: Create a 1-row data frame as input to deconstructSigs
sigs.input.df <- data.frame(t(signature_vector))
rownames(sigs.input.df) <- "Sample1"
colnames(sigs.input.df) = names(signature_vector) 

# Step 3: Load COSMIC reference signatures
data(signatures.cosmic)

# Step 4: Run signature decomposition
fit <- whichSignatures(
  tumor.ref = sigs.input.df,
  signatures.ref = cosmic_mat,
  sample.id = "Sample1",
  contexts.needed = FALSE,  # You already supplied trinucleotide context
  tri.counts.method = 'default'
)

weights_df = 
  data.frame(Signature = names(fit$weights),
             Weight = as.numeric(fit$weights))

SBS_order = weights_df$Signature
  
weights_df <- 
  weights_df %>%
  mutate(Signature = factor(Signature, levels = SBS_order)) %>%
  arrange(Signature)

data.frame(Signature = names(fit$weights),
           Weight = as.numeric(fit$weights)) %>%
  mutate(Signature = str_replace(Signature, "^w\\.", ""),
         Signature = str_replace(Signature, "\\.", " ")) %>%
  mutate(Signature = factor(Signature, levels = paste("Signature", 1:30)),
         Weight = Weight *100)

weights_df %>%
  ggplot() +
  geom_bar(aes(x = Signature,
               y = Weight*100),
           stat = "identity",
           fill = "grey80",
           color = "black") +
  coord_flip() + 
  theme_classic() +
  ylab("Percent Contribution")
ggsave(file.path(output_dir, "spectrum_decomposition.pdf"),
       height = 10,
       width = 4)


# Look at filtered variants from the single cell dataset ------------------

setwd(file.path(project_root, "paper_figures/data/Single_cell_bottlenecking_summary_statistics/"))
output_dir = file.path(project_root, "paper_figures/output/mutation_spectra/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

filtered_variants = 
  read.csv("sc_spectra_no_singletons.csv",
           sep = ",")

average_density = 
filtered_variants %>%
  as.data.frame() %>%
  dplyr::select(-X)

average_density = average_density/rowSums(average_density)

average_density <- colMeans(average_density, na.rm = TRUE)

average_df <- data.frame(
  mutation = names(average_density),
  value = as.numeric(average_density))


average_df$density = average_df$value

# Extract bases
average_df <- average_df %>%
  mutate(
    ref_base = substr(sub("^(...).(...).*$", "\\1", mutation), 2, 2),
    alt_base = substr(sub("^(...).(...).*$", "\\2", mutation), 2, 2),
    base_change = paste0(ref_base, ">", alt_base)
  )

# Apply all transformations
average_df <- average_df %>%
  mutate(
    # Normalize base change to pyrimidine context
    base_change_new = case_when(
      base_change == "A>C" ~ "T>G",
      base_change == "A>G" ~ "T>C",
      base_change == "A>T" ~ "T>A",
      base_change == "G>T" ~ "C>A",
      base_change == "G>C" ~ "C>G",
      base_change == "G>A" ~ "C>T",
      TRUE ~ base_change
    ),
    # Convert to COSMIC-style
    mutation_cosmic = sapply(paste0(mutation, ">", substr(alt_base, 1, 1)), to_cosmic_format),
    # Set factor order
    mutation_cosmic = factor(mutation_cosmic, levels = cosmic_order)
  )

average_df


plot_spectra <-function(x){
  mutation_colors <- c(
    "C>A" = "#00AEEF",
    "C>G" = "#000000",
    "C>T" = "#EE2E2F",
    "T>A" = "#BFBFBF",
    "T>C" = "#92D050",
    "T>G" = "#E9C3C3"
  )
  
  
  ggplot(x) +
    geom_bar(
      aes(x = mutation_cosmic, y = density, fill = base_change_new),
      stat = "identity",
      color = "black",
      size = 0.25
    ) +
    scale_fill_manual(values = mutation_colors, drop = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, size = 6,hjust = 0,vjust = 0.15))
}
plot_spectra(average_df)
ggsave(paste(output_dir,"spectra_filtered.pdf",sep = ""),
       height = 4,
       width = 12)


# Look at single cells ----------------------------------------------------

single_cell_spectrum <- filtered_variants %>%
  sample_n(1) %>%
  select(-X) %>%
  t() %>%
  as.data.frame()

# Add mutation names as rownames to a column
single_cell_spectrum$mutation <- rownames(single_cell_spectrum)
colnames(single_cell_spectrum)[1] <- "value"

# Calculate density
single_cell_spectrum$density <- single_cell_spectrum$value / sum(single_cell_spectrum$value)


# Extract bases
single_cell_spectrum <- 
  single_cell_spectrum %>%
  mutate(
    ref_base = substr(sub("^(...).(...).*$", "\\1", mutation), 2, 2),
    alt_base = substr(sub("^(...).(...).*$", "\\2", mutation), 2, 2),
    base_change = paste0(ref_base, ">", alt_base),
    # Normalize base change to pyrimidine context
    base_change_new = case_when(
      base_change == "A>C" ~ "T>G",
      base_change == "A>G" ~ "T>C",
      base_change == "A>T" ~ "T>A",
      base_change == "G>T" ~ "C>A",
      base_change == "G>C" ~ "C>G",
      base_change == "G>A" ~ "C>T",
      TRUE ~ base_change
    ),
    # Convert to COSMIC-style
    mutation_cosmic = sapply(paste0(mutation, ">", substr(alt_base, 1, 1)), to_cosmic_format),
    # Set factor order
    mutation_cosmic = factor(mutation_cosmic, levels = cosmic_order)
  )

plot_spectra(single_cell_spectrum)


# Sample 10 single cells
ten_cell_spectra <- filtered_variants %>%
  sample_n(10) %>%
  mutate(cell_id = paste0("Cell_", row_number())) %>%
  pivot_longer(-c(X, cell_id), names_to = "mutation", values_to = "value") %>%
  group_by(cell_id) %>%
  mutate(
    density = value / sum(value),
    ref_base = substr(sub("^(...).(...).*$", "\\1", mutation), 2, 2),
    alt_base = substr(sub("^(...).(...).*$", "\\2", mutation), 2, 2),
    base_change = paste0(ref_base, ">", alt_base),
    base_change_new = case_when(
      base_change == "A>C" ~ "T>G",
      base_change == "A>G" ~ "T>C",
      base_change == "A>T" ~ "T>A",
      base_change == "G>T" ~ "C>A",
      base_change == "G>C" ~ "C>G",
      base_change == "G>A" ~ "C>T",
      TRUE ~ base_change
    ),
    # Convert to COSMIC-style
    mutation_cosmic = sapply(paste0(mutation, ">", substr(alt_base, 1, 1)), to_cosmic_format),
    # Set factor order
    mutation_cosmic = factor(mutation_cosmic, levels = cosmic_order)
    
  ) %>%
  ungroup()

# OPTIONAL: define COSMIC-style order if needed
# ten_cell_spectra$mutation_cosmic <- factor(...)

plot_spectra(ten_cell_spectra) + 
  facet_wrap(~cell_id, ncol = 1) +
  theme(panel.background = element_blank(),
        strip.text = element_blank())
ggsave(paste(output_dir,"spectra_filtered_10_cells.pdf",sep = ""),
       height = 7,
       width = 8.5)

# Calculate spread of vectors ---------------------------------------------

library(philentropy)
set.seed(123)

single_cell_spectra <- 
  filtered_variants %>%
  select(-X)

single_cell_spectra <- as.data.frame(t(apply(single_cell_spectra, 1, function(x) x / sum(x))))

null_vector <- colMeans(single_cell_spectra)  # assuming it's a 1000 × 96 matrix

js_divs <- apply(single_cell_spectra, 1, function(vec) {
  distance(rbind(vec, null_vector), method = "jensen-shannon")
})


n_boot <- 1000

mat <- as.matrix(single_cell_spectra)
obs_jsd <- mean(apply(mat, 1, function(row) {
  distance(rbind(row, null_vector), method = "jensen-shannon")
}))

jsd_null_dist <- replicate(n_boot, {
  fake_mat <- t(apply(mat, 1, function(x) sample(null_vector)))  # shuffled version
  mean(apply(fake_mat, 1, function(row) {
    distance(rbind(row, null_vector), method = "jensen-shannon")
  }))
})

rbind(data.frame(type = "Empirical Null",
                 js_div = jsd_null_dist),
      data.frame(type = "Observed",
                 js_div = js_divs)) %>%
  ggplot() +
  geom_violin(aes(x=type,
                  y = js_div,
                  fill = type)) +
  geom_point(aes(x=type,
                  y = js_div,
                  fill = type)) +
  
  theme_classic() +
  ylab("JS Divergence") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1")
ggsave(paste(output_dir,"JS_divergence.pdf",sep = ""),
       height = 2,
       width = 3)


mean(js_divs)
sqrt(var(js_divs))

# Decomposition -----------------------------------------------------------
library(deconstructSigs)
library(tidyverse)
m = read.csv(cosmic_path,sep = "\t") %>%
  as_tibble() %>%
  column_to_rownames(var = "Type")

cosmic_raw <- read.csv(cosmic_path, sep = "\t", check.names = FALSE)

# First column is trinucleotide context
contexts <- cosmic_raw[, 1]

# Remaining columns are the signature profiles
cosmic_mat <- as.data.frame(t(cosmic_raw[, -1]))
colnames(cosmic_mat) <- contexts
# signatures_SBS = rownames(cosmic_mat)
cosmic_mat <- as.data.frame(cosmic_mat)

# Step 1: Create a named vector with trinucleotide contexts
signature_vector <- setNames(average_df$density, average_df$mutation_cosmic)

# Step 2: Create a 1-row data frame as input to deconstructSigs
sigs.input.df <- data.frame(t(signature_vector))
rownames(sigs.input.df) <- "Sample1"
colnames(sigs.input.df) = names(signature_vector) 

# Step 3: Load COSMIC reference signatures
data(signatures.cosmic)

# Step 4: Run signature decomposition
fit <- whichSignatures(
  tumor.ref = sigs.input.df,
  signatures.ref = cosmic_mat,
  sample.id = "Sample1",
  contexts.needed = FALSE,  # You already supplied trinucleotide context
  tri.counts.method = 'default'
)

weights_df = 
  data.frame(Signature = names(fit$weights),
             Weight = as.numeric(fit$weights))

SBS_order = weights_df$Signature

weights_df <- 
  weights_df %>%
  mutate(Signature = factor(Signature, levels = SBS_order)) %>%
  arrange(Signature)

data.frame(Signature = names(fit$weights),
           Weight = as.numeric(fit$weights)) %>%
  mutate(Signature = str_replace(Signature, "^w\\.", ""),
         Signature = str_replace(Signature, "\\.", " ")) %>%
  mutate(Signature = factor(Signature, levels = paste("Signature", 1:30)),
         Weight = Weight *100)

weights_df %>%
  filter(Weight > 0) %>%
  ggplot() +
  geom_bar(aes(x = Signature,
               y = Weight*100),
           stat = "identity",
           fill = "grey80",
           color = "black") +
  coord_flip() + 
  theme_classic() +
  ylab("Percent Contribution")
ggsave(paste(output_dir,"SBS_contribution.pdf",sep = ""),
       height = 2,
       width = 3)
