# K562 Pol-epsilon (P286R) mutation accumulation figures.
#
# Ported from notebooks/K562_mut_accumulation.ipynb (D. Mullane). The variant
# filtering lives upstream in that notebook (cellspec / anndata); this script
# consumes the derived tables it writes and re-renders the panels in the
# paper_figures house style.
#
# Upstream filtering chain reproduced in the notebook, in order:
#   1. spc.pp.filter_to_snps()                      -- SNVs only
#   2. spc.pp.annotate_contexts(hg38)               -- trinucleotide contexts
#   3. spc.tl.compute_bulk_vaf(target_dp = 100)
#   4. bulk_vaf < 0.15                              -- drop ancestral variants
#   5. per sample: filter_by_coverage(min_depth = 10, max_depth = 100)
#   6. per sample: AD > 1                           -- >1 alt read
#   7. per sample: VAF > 0.15 at target_dp = 100    -- "high_vaf" call
# Sites passing (5-7) in a sample are that sample's high-confidence set; the
# per-lineage accumulation is the first-appearance partition of those sets
# across passages P1/P2/P3.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(lme4)
set.seed(50)

data_dir <- file.path(project_root, "paper_figures/data/K562_mut_accumulation/")
output_dir <- file.path(project_root, "paper_figures/output/K562_mut_accumulation/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

genotype_colors <- c("WT" = "#438CFD", "P286R" = "#ff1a5e")

# 21 days per passage (3 weeks); see proliferation section below.
days_per_passage <- 21

stats_log <- character(0)
log_stat <- function(...) {
  line <- paste0(...)
  stats_log <<- c(stats_log, line)
  cat(line, "\n", sep = "")
}
# print a table to the console and keep it in the stats log
log_table <- function(df) {
  lines <- capture.output(print(as.data.frame(df), row.names = FALSE))
  stats_log <<- c(stats_log, lines)
  cat(lines, sep = "\n")
  cat("\n")
}


# Proliferation -----------------------------------------------------------
# Growth curves collected by I. Campbell; counts are transcribed from the
# notebook (cell 6). Seeded at 3e3 cells, counted on days 0/5/7/10.

prolif <- tribble(
  ~sample,             ~clone, ~genotype, ~day_0, ~day_5,  ~day_7,  ~day_10,
  "WT_clone1_rep1",    "1",    "WT",      3000,   1.23e5,  8.27e5,  1.44e6,
  "WT_clone1_rep2",    "1",    "WT",      3000,   1.11e5,  1.30e6,  1.82e6,
  "WT_clone2_rep1",    "2",    "WT",      3000,   3.28e5,  3.69e5,  1.84e6,
  "WT_clone2_rep2",    "2",    "WT",      3000,   1.11e5,  6.98e5,  8.21e5,
  "WT_clone3_rep1",    "3",    "WT",      3000,   3.52e5,  1.71e6,  1.34e6,
  "WT_clone3_rep4",    "3",    "WT",      3000,   5.04e5,  1.52e6,  2.38e6,
  "P286R_clone1_rep1", "1",    "P286R",   3000,   2.70e5,  7.57e5,  7.98e5,
  "P286R_clone1_rep2", "1",    "P286R",   3000,   7.04e4,  3.17e5,  9.21e5,
  "P286R_clone2_rep1", "2",    "P286R",   3000,   1.17e4,  5.86e4,  4.69e5,
  "P286R_clone2_rep2", "2",    "P286R",   3000,   2.93e4,  3.52e4,  8.33e5,
  "P286R_clone3_rep1", "3",    "P286R",   3000,   3.52e4,  7.62e4,  1.51e6,
  "P286R_clone3_rep2", "3",    "P286R",   3000,   7.03e4,  6.54e4,  7.21e5
)

prolif_long <-
  prolif %>%
  pivot_longer(starts_with("day_"), names_to = "day", values_to = "count") %>%
  mutate(
    day = as.numeric(sub("day_", "", day)),
    genotype = factor(genotype, levels = c("WT", "P286R"))
  )

ggplot(prolif_long, aes(x = day, y = count, color = genotype, fill = genotype)) +
  geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, linewidth = 0.6) +
  geom_point(size = 1.4, shape = 21, color = "black", stroke = 0.3) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  scale_y_continuous(labels = function(y) y / 1e6) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  xlab("Time (days)") +
  ylab("Cell count (millions)")
ggsave(file.path(output_dir, "proliferation.png"),
       dpi = 600, bg = "transparent",
       height = 2.5, width = 2.75)

# Log-linear growth rate per replicate: ln(N) = ln(N0) + r*t
growth_rates <-
  prolif_long %>%
  filter(count > 0) %>%
  group_by(sample, clone, genotype) %>%
  summarise(
    fit = list(lm(log(count) ~ day)),
    .groups = "drop"
  ) %>%
  mutate(
    r_per_day = vapply(fit, function(f) coef(f)[["day"]], numeric(1)),
    doublings_per_day = r_per_day / log(2),
    r_squared = vapply(fit, function(f) summary(f)$r.squared, numeric(1))
  ) %>%
  select(-fit)

log_stat("=== Log-linear growth rates (per replicate) ===")
log_table(growth_rates)

doublings_per_passage <-
  growth_rates %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    mean_doublings_per_day = mean(doublings_per_day),
    sd = sd(doublings_per_day),
    sem = sd(doublings_per_day) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    doublings = mean_doublings_per_day * days_per_passage,
    # half-width of the 95% CI, propagated through the linear day -> passage scaling
    doublings_ci95 = qt(0.975, df = n - 1) * sem * days_per_passage
  )

log_stat("")
log_stat("=== Doublings per ", days_per_passage, "-day passage ===")
log_table(doublings_per_passage)

welch <- t.test(r_per_day ~ genotype, data = growth_rates, var.equal = FALSE)
log_stat("")
log_stat(sprintf("Welch's t-test on r_per_day (WT vs P286R): t = %.3f, df = %.2f, p = %.4g",
                 welch$statistic, welch$parameter, welch$p.value))

write.csv(growth_rates, file.path(output_dir, "growth_rates.csv"), row.names = FALSE)
write.csv(doublings_per_passage, file.path(output_dir, "doublings_per_passage.csv"),
          row.names = FALSE)


# Mutation accumulation ---------------------------------------------------
# NOTE: the `accumulated_mutations` / `SBS10a` columns come straight from the
# notebook's first-appearance partition. For the 3-passage lineages that
# partition sets `split_2 = s3 - split_1`, which leaves the shared ancestral
# set (`split_0`) inside `split_2`, so the split_3 total is inflated by exactly
# the split_1 count. `accumulated_split_3_corrected` below reports the value
# with that double-count removed; the plotted values are left as published.

time_series <-
  read.csv(file.path(data_dir, "mutation_accumulation.csv"),
           check.names = FALSE) %>%
  select(lineage, split, passage, genotype,
         accumulated_mutations, accumulated_mutations_norm,
         cumulative_doublings, SBS10a) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "P286R")),
    weeks = passage * days_per_passage / 7
  )

split_1_counts <-
  time_series %>%
  filter(split == "split_1") %>%
  select(lineage, split_1_count = accumulated_mutations)

accumulation_check <-
  time_series %>%
  filter(split == "split_3") %>%
  left_join(split_1_counts, by = "lineage") %>%
  transmute(
    lineage, genotype,
    accumulated_split_3 = accumulated_mutations,
    accumulated_split_3_corrected = accumulated_mutations - split_1_count,
    percent_inflation = 100 * split_1_count / accumulated_mutations
  )

log_stat("")
log_stat("=== split_3 double-count check (upstream set-logic) ===")
log_table(accumulation_check)
write.csv(accumulation_check, file.path(output_dir, "accumulation_split3_check.csv"),
          row.names = FALSE)

week_breaks <- sort(unique(time_series$weeks))

plot_accumulation <- function(df, yvar, ylab, y_scale = 1) {
  ggplot(df, aes(x = weeks, y = .data[[yvar]] / y_scale,
                 color = genotype, fill = genotype)) +
    geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, linewidth = 0.6) +
    geom_point(size = 1.4, shape = 21, color = "black", stroke = 0.3) +
    scale_color_manual(values = genotype_colors) +
    scale_fill_manual(values = genotype_colors) +
    scale_x_continuous(breaks = week_breaks) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    xlab("Time (weeks)") +
    ylab(ylab)
}

plot_accumulation(time_series, "accumulated_mutations_norm",
                  "SBS per doubling")
ggsave(file.path(output_dir, "accumulated_mutations.png"),
       dpi = 600, bg = "transparent",
       height = 2.5, width = 2.75)

plot_accumulation(time_series, "SBS10a",
                  "SBS10a (thousands)", y_scale = 1e3)
ggsave(file.path(output_dir, "accumulated_SBS10a.png"),
       dpi = 600, bg = "transparent",
       height = 2.5, width = 2.75)

plot_accumulation(time_series, "accumulated_mutations",
                  "Accumulated SBS (thousands)", y_scale = 1e3)
ggsave(file.path(output_dir, "accumulated_mutations_raw.png"),
       dpi = 600, bg = "transparent",
       height = 2.5, width = 2.75)


# Mixed models ------------------------------------------------------------
# Random intercept per lineage, matching statsmodels mixedlm(re_formula = "~1").
# Wald z p-values are reported so the numbers line up with statsmodels, which
# does not apply a Satterthwaite correction either.

fit_accumulation_model <- function(df, response) {
  form <- as.formula(paste(response, "~ passage * genotype + (1 | lineage)"))
  model <- lmer(form, data = df, REML = TRUE)
  coefs <- summary(model)$coefficients
  out <- data.frame(
    response = response,
    term = rownames(coefs),
    estimate = coefs[, "Estimate"],
    std_error = coefs[, "Std. Error"],
    statistic = coefs[, "t value"],
    p_value = 2 * pnorm(abs(coefs[, "t value"]), lower.tail = FALSE),
    row.names = NULL
  )
  list(model = model, coefs = out)
}

model_tables <- list()
for (response in c("accumulated_mutations_norm", "SBS10a", "accumulated_mutations")) {
  fit <- fit_accumulation_model(time_series, response)
  model_tables[[response]] <- fit$coefs

  b_passage <- fit$coefs$estimate[fit$coefs$term == "passage"]
  b_inter <- fit$coefs$estimate[fit$coefs$term == "passage:genotypeP286R"]
  se_passage <- fit$coefs$std_error[fit$coefs$term == "passage"]
  se_inter <- fit$coefs$std_error[fit$coefs$term == "passage:genotypeP286R"]
  p_inter <- fit$coefs$p_value[fit$coefs$term == "passage:genotypeP286R"]

  # SE of the summed coefficients, using the real covariance rather than
  # assuming independence (the notebook's markdown assumes zero covariance).
  V <- as.matrix(vcov(fit$model))
  idx <- c("passage", "passage:genotypeP286R")
  se_sum <- sqrt(sum(V[idx, idx]))

  log_stat("")
  log_stat("=== Mixed model: ", response, " ~ passage * genotype + (1 | lineage) ===")
  log_table(fit$coefs)
  log_stat(sprintf("WT rate            = %.1f per passage (SE %.1f)", b_passage, se_passage))
  log_stat(sprintf("P286R excess rate  = %.1f per passage (SE %.1f)", b_inter, se_inter))
  log_stat(sprintf("P286R total rate   = %.1f per passage (SE %.1f)", b_passage + b_inter, se_sum))
  log_stat(sprintf("Interaction p      = %.3g", p_inter))
}

model_coefs <- do.call(rbind, model_tables)
write.csv(model_coefs, file.path(output_dir, "mixed_model_coefficients.csv"),
          row.names = FALSE)


# Per-division SBS10a rate ------------------------------------------------
# mu = M / D, with M = mutations per passage from the mixed model and
# D = doublings per passage from the proliferation assay. The two come from
# independent experiments, so the delta method for a ratio applies.

sbs10a_coefs <- model_tables[["SBS10a"]]
V <- as.matrix(vcov(fit_accumulation_model(time_series, "SBS10a")$model))
idx <- c("passage", "passage:genotypeP286R")

rate_rows <- list()
for (g in c("WT", "P286R")) {
  if (g == "WT") {
    M <- sbs10a_coefs$estimate[sbs10a_coefs$term == "passage"]
    se_M <- sbs10a_coefs$std_error[sbs10a_coefs$term == "passage"]
  } else {
    M <- sum(sbs10a_coefs$estimate[sbs10a_coefs$term %in% idx])
    se_M <- sqrt(sum(V[idx, idx]))
  }
  D <- doublings_per_passage$doublings[doublings_per_passage$genotype == g]
  # convert the reported 95% CI half-width back to an SE
  n_g <- doublings_per_passage$n[doublings_per_passage$genotype == g]
  se_D <- doublings_per_passage$doublings_ci95[doublings_per_passage$genotype == g] /
    qt(0.975, df = n_g - 1)

  mu <- M / D
  se_mu <- abs(mu) * sqrt((se_M / M)^2 + (se_D / D)^2)

  rate_rows[[g]] <- data.frame(
    genotype = g,
    mutations_per_passage = M,
    se_mutations_per_passage = se_M,
    doublings_per_passage = D,
    se_doublings_per_passage = se_D,
    sbs10a_per_division = mu,
    se_sbs10a_per_division = se_mu,
    ci95_low = mu - 1.96 * se_mu,
    ci95_high = mu + 1.96 * se_mu
  )
}
per_division_rate <- do.call(rbind, rate_rows)

log_stat("")
log_stat("=== SBS10a per division (delta method) ===")
log_table(per_division_rate)
write.csv(per_division_rate, file.path(output_dir, "sbs10a_per_division.csv"),
          row.names = FALSE)


# Mutation spectra --------------------------------------------------------
# 96-context spectra rendered with the same COSMIC ordering and palette as
# mutation_spectra.R.

revcomp_trinuc <- function(trinuc) {
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  paste0(rev(comp[strsplit(trinuc, "")[[1]]]), collapse = "")
}

# Convert TTT>TGT -> T[T>G]T
to_cosmic_format <- function(tri_mut) {
  ref <- substr(tri_mut, 1, 3)
  alt <- substr(tri_mut, 5, 7)
  pre <- substr(ref, 1, 1)
  mid <- substr(ref, 2, 2)
  post <- substr(ref, 3, 3)
  mid_alt <- substr(alt, 2, 2)

  # If middle base is purine, reverse complement
  if (mid %in% c("A", "G")) {
    ref_rc <- revcomp_trinuc(ref)
    alt_rc <- revcomp_trinuc(alt)
    pre <- substr(ref_rc, 1, 1)
    mid <- substr(ref_rc, 2, 2)
    post <- substr(ref_rc, 3, 3)
    mid_alt <- substr(alt_rc, 2, 2)
  }

  paste0(pre, "[", mid, ">", mid_alt, "]", post)
}

# COSMIC 96-order definition
bases <- c("A", "C", "G", "T")
subs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
cosmic_order <- as.character(unlist(lapply(subs, function(sub) {
  ref <- substr(sub, 1, 1)
  alt <- substr(sub, 3, 3)
  sapply(bases, function(pre) {
    sapply(bases, function(post) {
      paste0(pre, "[", ref, ">", alt, "]", post)
    })
  })
})))

mutation_colors <- c(
  "C>A" = "#00AEEF",
  "C>G" = "#000000",
  "C>T" = "#EE2E2F",
  "T>A" = "#BFBFBF",
  "T>C" = "#92D050",
  "T>G" = "#E9C3C3"
)

annotate_spectrum <- function(df) {
  df %>%
    mutate(
      ref_base = substr(mutation, 2, 2),
      alt_base = substr(mutation, 6, 6),
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
      mutation_cosmic = factor(vapply(mutation, to_cosmic_format, character(1)),
                               levels = cosmic_order)
    )
}

plot_spectra <- function(x, yvar = "density") {
  ggplot(x) +
    geom_bar(
      aes(x = mutation_cosmic, y = .data[[yvar]], fill = base_change_new),
      stat = "identity",
      color = "black",
      linewidth = 0.25
    ) +
    scale_fill_manual(values = mutation_colors, drop = FALSE) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, size = 6, hjust = 0, vjust = 0.15))
}

# De novo spectra, pooled by construct (AAVS = WT, PolE = P286R)
spectrum_wide <- read.csv(file.path(data_dir, "spectrum.csv"),
                          check.names = FALSE, row.names = 1)

spectrum_by_construct <-
  spectrum_wide %>%
  rownames_to_column("sample") %>%
  mutate(construct = sub("_.*$", "", sample)) %>%
  select(-sample) %>%
  group_by(construct) %>%
  summarise(across(everything(), sum), .groups = "drop")

for (construct in spectrum_by_construct$construct) {
  spec <-
    spectrum_by_construct %>%
    filter(construct == !!construct) %>%
    select(-construct) %>%
    pivot_longer(everything(), names_to = "mutation", values_to = "value") %>%
    mutate(density = value / sum(value)) %>%
    annotate_spectrum()

  plot_spectra(spec, yvar = "value") +
    scale_y_continuous(labels = function(y) y / 1e3)
  ggsave(file.path(output_dir, paste0("spectrum_", construct, ".pdf")),
         height = 4, width = 12)

  plot_spectra(spec, yvar = "density")
  ggsave(file.path(output_dir, paste0("spectrum_", construct, "_density.pdf")),
         height = 4, width = 12)

  write.csv(spec, file.path(output_dir, paste0("spectrum_", construct, ".csv")),
            row.names = FALSE)
}

# Background spectrum (bulk_vaf > 0.15 sites, i.e. the ancestral variants
# removed by filter step 4 above)
background <-
  read.csv(file.path(data_dir, "spectrum_background.csv"),
           col.names = c("mutation", "value")) %>%
  mutate(density = value / sum(value)) %>%
  annotate_spectrum()

plot_spectra(background, yvar = "value") +
  scale_y_continuous(labels = function(y) y / 1e6)
ggsave(file.path(output_dir, "spectrum_background.pdf"),
       height = 4, width = 12)

plot_spectra(background, yvar = "density")
ggsave(file.path(output_dir, "spectrum_background_density.pdf"),
       height = 4, width = 12)


# COSMIC signature activities ---------------------------------------------
# SigProfilerAssignment activities from the notebook, restricted to the ten
# signatures with the largest P286R - WT difference.

sig_activities <- read.csv(file.path(data_dir, "top10_EDT.csv"),
                           check.names = FALSE)
names(sig_activities)[1] <- "sample"

sig_order <- setdiff(names(sig_activities), "sample")

sig_per_sample <-
  sig_activities %>%
  filter(!sample %in% c("WT", "P286R", "difference")) %>%
  mutate(genotype = ifelse(startsWith(sample, "PolE"), "P286R", "WT")) %>%
  pivot_longer(all_of(sig_order), names_to = "signature", values_to = "activity") %>%
  mutate(
    signature = factor(signature, levels = rev(sig_order)),
    genotype = factor(genotype, levels = c("WT", "P286R"))
  )

ggplot(sig_per_sample, aes(x = signature, y = activity / 1e3, fill = genotype)) +
  stat_summary(fun = mean, geom = "bar", color = "black", linewidth = 0.25,
               position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(position = position_dodge(width = 0.8), size = 0.6,
             shape = 21, color = "black", stroke = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = genotype_colors) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  ) +
  xlab(NULL) +
  ylab("Mean SBS (thousands)")
ggsave(file.path(output_dir, "signature_activities.pdf"),
       height = 4, width = 4)

sig_difference <-
  sig_activities %>%
  filter(sample == "difference") %>%
  pivot_longer(all_of(sig_order), names_to = "signature", values_to = "difference") %>%
  mutate(signature = factor(signature, levels = rev(sig_order)))

ggplot(sig_difference) +
  geom_bar(aes(x = signature, y = difference / 1e3),
           stat = "identity", fill = "grey80", color = "black", linewidth = 0.25) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  ) +
  xlab(NULL) +
  ylab("P286R - WT (thousands)")
ggsave(file.path(output_dir, "signature_difference.pdf"),
       height = 4, width = 3.5)

writeLines(stats_log, file.path(output_dir, "stats.txt"))
cat("\nWrote outputs to:", output_dir, "\n")
