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

result <- fit_bayesian_transport(
  dat = dat,
  targets = targets,
  moderators = moderators,
  effect_scale = "log_rr",
  chains = 2,
  iter = 800,
  warmup = 400,
  seed = 20260404,
  refresh = 100
)

output_dir <- file.path(project_root, "outputs", "demo_bayesian")
save_bayesian_results(result, output_dir)

cat("Bayesian demo complete.\n")
cat(sprintf("Outputs written to: %s\n\n", output_dir))
print(result$target_posterior)
