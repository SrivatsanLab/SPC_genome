# Shared 96-context mutation spectrum plotting.
#
# Extracted verbatim from mutation_spectra.R so other scripts can render
# spectra in the same style. Source it with:
#   source(file.path(project_root, "paper_figures/scripts/spectrum_utils.R"))
#
# Expects `mutation` labels in TCT>TAT form (ref trinucleotide, ">",
# alt trinucleotide).

library(ggplot2)
library(dplyr)

# Reverse complement function for trinucleotides
revcomp_trinuc <- function(trinuc) {
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  paste0(rev(comp[strsplit(trinuc, "")[[1]]]), collapse = "")
}

# Convert TTT>TGT -> T[T>G]T
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

mutation_colors <- c(
  "C>A" = "#00AEEF",
  "C>G" = "#000000",
  "C>T" = "#EE2E2F",
  "T>A" = "#BFBFBF",
  "T>C" = "#92D050",
  "T>G" = "#E9C3C3"
)

# Add the base-change, pyrimidine-normalised base change and COSMIC-ordered
# mutation factor used by plot_spectra(). This is the mutate chain that
# mutation_spectra.R applies to each spectrum before plotting.
annotate_spectrum <- function(df) {
  df %>%
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
      mutation_cosmic = sapply(mutation, to_cosmic_format),
      # Set factor order
      mutation_cosmic = factor(mutation_cosmic, levels = cosmic_order)
    )
}

# outline_width = 0 drops the bar outline, which otherwise swamps the fill
# colour when 96 bars are squeezed into a small (e.g. facetted) panel.
plot_spectra <- function(x, yvar = "density", outline_width = 0.25) {
  ggplot(x) +
    geom_bar(
      aes(x = mutation_cosmic, y = .data[[yvar]], fill = base_change_new),
      stat = "identity",
      color = if (outline_width > 0) "black" else NA,
      linewidth = outline_width
    ) +
    scale_fill_manual(values = mutation_colors, drop = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, size = 6, hjust = 0, vjust = 0.15))
}
