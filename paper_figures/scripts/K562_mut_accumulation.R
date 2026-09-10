# K562 Pol-epsilon (P286R) mutation accumulation: SBS10a accumulation over
# passages, and the 96-context mutation spectra.
#
# Ported from notebooks/K562_mut_accumulation.ipynb (D. Mullane). The variant
# filtering lives upstream in that notebook (cellspec / anndata); this script
# consumes the spectrum tables it writes and renders them with the shared
# plotting code in spectrum_utils.R.
#
# Upstream filtering chain reproduced in the notebook, in order:
#   1. spc.pp.filter_to_snps()                      -- SNVs only
#   2. spc.pp.annotate_contexts(hg38)               -- trinucleotide contexts
#   3. spc.tl.compute_bulk_vaf(target_dp = 100)
#   4. bulk_vaf < 0.15                              -- drop ancestral variants
#   5. per sample: filter_by_coverage(min_depth = 10, max_depth = 100)
#   6. per sample: AD > 1                           -- >1 alt read
#   7. per sample: VAF > 0.15 at target_dp = 100    -- "high_vaf" call
# spectrum.csv is the de novo spectrum over sites passing (1-7);
# spectrum_background.csv is the ancestral spectrum over the bulk_vaf > 0.15
# sites removed by step 4. mutation_accumulation.csv holds the per-lineage
# first-appearance counts per passage.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

source(file.path(project_root, "paper_figures/scripts/spectrum_utils.R"))

data_dir <- file.path(project_root, "paper_figures/data/K562_mut_accumulation/")
output_dir <- file.path(project_root, "paper_figures/output/K562_mut_accumulation/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# SBS10a accumulation -----------------------------------------------------
# Notebook cell 55, in the house style. One passage is 21 days (3 weeks).
# Points and axis text are sized up relative to the other panels so the
# figure stays readable when it is scaled down.

genotype_colors <- c("WT" = "#438CFD", "P286R" = "#ff1a5e")

# row.names = 1 drops the unnamed index column the notebook writes
time_series <-
  read.csv(file.path(data_dir, "mutation_accumulation.csv"), row.names = 1) %>%
  select(lineage, split, passage, genotype, accumulated_mutations, SBS10a) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "P286R")),
    weeks = passage * 3
  )

ggplot(time_series, aes(x = weeks, y = SBS10a / 1e3,
                        color = genotype, fill = genotype)) +
  # se = FALSE: the default OLS band treats all 17 points as independent when
  # they come from 6 lineages sampled over passages, so it is the wrong
  # interval for this design. The notebook drops it in cell 55 too.
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.5) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  scale_x_continuous(breaks = sort(unique(time_series$weeks))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15)
  ) +
  xlab("Time (weeks)") +
  ylab("SBS10a (thousands)")
ggsave(file.path(output_dir, "accumulated_SBS10a.pdf"),
       bg = "transparent",
       height = 3.25, width = 3.5)


# De novo spectra ---------------------------------------------------------
# Per-sample 96-context counts, pooled by construct (AAVS = WT, PolE = P286R).

spectrum_wide <- read.csv(file.path(data_dir, "spectrum.csv"),
                          check.names = FALSE, row.names = 1)

spectrum_by_construct <-
  spectrum_wide %>%
  rownames_to_column("sample") %>%
  mutate(construct = sub("_.*$", "", sample)) %>%
  select(-sample) %>%
  group_by(construct) %>%
  summarise(across(everything(), sum), .groups = "drop")

pooled_spectra <- list()

for (this_construct in spectrum_by_construct$construct) {
  spec <-
    spectrum_by_construct %>%
    filter(construct == this_construct) %>%
    select(-construct) %>%
    pivot_longer(everything(), names_to = "mutation", values_to = "value") %>%
    mutate(density = value / sum(value)) %>%
    annotate_spectrum()

  # counts, in thousands
  plot_spectra(spec, yvar = "value") +
    scale_y_continuous(labels = function(y) y / 1e3)
  ggsave(file.path(output_dir, paste0("spectrum_", this_construct, ".pdf")),
         height = 4, width = 12)

  plot_spectra(spec)
  ggsave(file.path(output_dir, paste0("spectrum_", this_construct, "_density.pdf")),
         height = 4, width = 12)

  write.csv(spec, file.path(output_dir, paste0("spectrum_", this_construct, ".csv")),
            row.names = FALSE)

  pooled_spectra[[this_construct]] <- spec
}


# Per-sample spectra ------------------------------------------------------
# All 17 samples faceted, densities so lineages are comparable.

spectrum_per_sample <-
  spectrum_wide %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "mutation", values_to = "value") %>%
  group_by(sample) %>%
  mutate(density = value / sum(value)) %>%
  ungroup() %>%
  annotate_spectrum()

plot_spectra(spectrum_per_sample, outline_width = 0) +
  facet_wrap(~sample, ncol = 3, scales = "free_y") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(file.path(output_dir, "spectrum_per_sample.pdf"),
       height = 10, width = 12)

write.csv(spectrum_per_sample,
          file.path(output_dir, "spectrum_per_sample.csv"),
          row.names = FALSE)


# Background spectrum -----------------------------------------------------
# Ancestral variants (bulk_vaf > 0.15), i.e. the sites removed by filter
# step 4 above.

background <-
  read.csv(file.path(data_dir, "spectrum_background.csv"),
           col.names = c("mutation", "value")) %>%
  mutate(density = value / sum(value)) %>%
  annotate_spectrum()

# counts, in millions
plot_spectra(background, yvar = "value") +
  scale_y_continuous(labels = function(y) y / 1e6)
ggsave(file.path(output_dir, "spectrum_background.pdf"),
       height = 4, width = 12)

plot_spectra(background)
ggsave(file.path(output_dir, "spectrum_background_density.pdf"),
       height = 4, width = 12)

write.csv(background, file.path(output_dir, "spectrum_background.csv"),
          row.names = FALSE)


# Shared-scale densities --------------------------------------------------
# Background, AAVS and PolE densities as separate panels but on one common
# y axis, so the SBS10a T[C>A]T peak reads against the baseline spectra
# rather than against a per-panel rescaling.

shared_scale_spectra <- list(
  background = background,
  AAVS = pooled_spectra[["AAVS"]],
  PolE = pooled_spectra[["PolE"]]
)

# common limit with a little headroom above the tallest bar in any panel
density_limit <- max(vapply(shared_scale_spectra,
                            function(x) max(x$density), numeric(1))) * 1.05

for (this_spectrum in names(shared_scale_spectra)) {
  # add_context_axis() sets the shared y limit via coord_cartesian, so the
  # annotation drawn below y = 0 is not clipped away
  add_context_axis(plot_spectra(shared_scale_spectra[[this_spectrum]]),
                   ymax = density_limit)
  ggsave(file.path(output_dir,
                   paste0("spectrum_", this_spectrum, "_density_shared_scale.pdf")),
         height = 2.5, width = 8)
}

cat("Shared density y limit:", density_limit, "\n")

cat("\nWrote spectra to:", output_dir, "\n")
