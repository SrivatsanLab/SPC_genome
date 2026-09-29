project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
script_dir <- file.path(project_root, "paper_figures", "scripts")

scripts <- c(
  "DNA_SPC_signal.R",
  "Protein_SPC_signal.R",
  "plot_segmented_SPC_oligomers.R",
  "plot_spc_sizes.R",
  "plot_amplicon_cooccurance.R",
  "plot_genome_barnyard.R",
  "plot_gini_lorenz.R",
  "Encode_variants.R",
  "plot_mutations_bulk_WGA.R",
  "plot_clone_abundance_and_mutation_types.R",
  "mutation_spectra.R",
  "analyse_flowdata.R",
  "anneufinder_plot.R",
  "draw_trees.R"
)

results <- lapply(scripts, function(s) {
  path <- file.path(script_dir, s)
  if (!file.exists(path)) {
    return(data.frame(script = s, status = "missing", message = "script file not found"))
  }

  message(sprintf("Running %s", s))
  out <- tryCatch(
    {
      setwd(project_root)
      source(path, local = new.env(parent = globalenv()))
      data.frame(script = s, status = "ok", message = "")
    },
    error = function(e) {
      data.frame(script = s, status = "error", message = conditionMessage(e))
    }
  )
  out
})

summary_df <- do.call(rbind, results)
summary_dir <- file.path(project_root, "paper_figures", "output")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(summary_dir, "run_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

print(summary_df)
cat("\nSummary written to:", summary_path, "\n")

if (any(summary_df$status == "error")) {
  quit(status = 1)
}
