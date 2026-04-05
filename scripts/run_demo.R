get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) == 0) {
    normalizePath(getwd())
  } else {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
  }
}

project_root <- get_project_root()
setwd(project_root)

pkgload::load_all(project_root, export_all = FALSE, helpers = FALSE, quiet = TRUE)

dat <- read_example_studies()
targets <- read_example_targets()

moderators <- c("mean_age", "female_pct", "baseline_risk", "followup_months")

result <- fit_transport_engine(
  dat = dat,
  targets = targets,
  moderators = moderators,
  yi_col = "yi",
  vi_col = "vi",
  study_col = "study",
  target_label_col = "label",
  effect_scale = "log_rr",
  rob_col = "risk_of_bias",
  high_risk_values = c("high")
)

output_dir <- file.path(project_root, "outputs", "demo_frequentist")
save_transport_results(result, output_dir)

cat("Frequentist demo complete.\n")
cat(sprintf("Outputs written to: %s\n\n", output_dir))
print(result$target_estimates)
