# FACS-style cell-cycle profile of the K562 sc_PolE cells inferred from coverage.
#
# x = inferred DNA content relative to G1 = 1 + fraction of the genome
# replicated: S-phase cells (K562_replication_timing_coverage.R phase call)
# take the replicated fraction from the replication fit
# (scripts/utils/replication_regress_segment.R, mean of the fitted logistic
# p(RT) over bins); non-replicating cells are placed at 1. Coverage is
# normalised per cell, so G1 and G2 cannot be told apart: both sit in the
# left peak, and there is no 4N (x = 2) G2 peak as in a DNA-content stain.
# y = cells: histogram (bin width 0.05) with a kernel density scaled to counts.
# Written without legends or panel text beyond axis titles.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(ggplot2))

fit <- read.csv(file.path(project_root, "results/K562_tree/aneufinder_fragment_gc_rtcorr/output/replication_fit.csv"))
phase <- read.csv(file.path(project_root, "paper_figures/output/K562_replication_timing_coverage/per_cell_phase.csv"))
output_dir <- file.path(project_root, "paper_figures/output/K562_cell_cycle_profile/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- merge(phase[, c("cell", "phase")], fit[, c("cell", "a", "replicated_frac", "r2")], by = "cell")
d$dna_content <- ifelse(d$phase == "S" & !is.na(d$replicated_frac), 1 + d$replicated_frac, 1)
write.csv(d, file.path(output_dir, "inferred_dna_content.csv"), row.names = FALSE)

binw <- 0.05
p <- ggplot(d, aes(dna_content)) +
  geom_histogram(binwidth = binw, center = 1, fill = "grey75", colour = NA) +
  geom_density(aes(y = after_stat(count) * binw), bw = 0.03, colour = "black", linewidth = 0.4) +
  scale_x_continuous(limits = c(0.95, 2.05), breaks = c(1, 1.25, 1.5, 1.75, 2),
                     labels = c("1 (G1/G2)", "1.25", "1.5", "1.75", "2")) +
  labs(x = "Inferred DNA content (G1 = 1)", y = "Cells") +
  theme_classic(base_size = 8)
ggsave(file.path(output_dir, "cell_cycle_profile.pdf"), p, width = 3.5, height = 2.5)
ggsave(file.path(output_dir, "cell_cycle_profile.svg"), p, width = 3.5, height = 2.5)
ggsave(file.path(output_dir, "cell_cycle_profile.png"), p, width = 3.5, height = 2.5, dpi = 200, bg = "white")

cat(sprintf("%d cells: %d at 1 (G1/G2), %d in S; S-phase DNA content quantiles 10/50/90%%: %s\n",
            nrow(d), sum(d$dna_content == 1), sum(d$dna_content > 1),
            paste(round(quantile(d$dna_content[d$dna_content > 1], c(.1, .5, .9)), 2), collapse = " ")))
cat("Outputs written to", output_dir, "\n")
