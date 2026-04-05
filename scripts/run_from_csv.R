args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/run_from_csv.R <studies_csv> <targets_csv> <moderator1,moderator2,...> [output_dir]"
    ),
    call. = FALSE
  )
}

get_project_root <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)

  if (length(file_arg) == 0) {
    normalizePath(getwd())
  } else {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."))
  }
}

project_root <- get_project_root()
setwd(project_root)

pkgload::load_all(project_root, export_all = FALSE, helpers = FALSE, quiet = TRUE)

studies_csv <- normalizePath(args[1], mustWork = TRUE)
targets_csv <- normalizePath(args[2], mustWork = TRUE)
moderators <- trimws(strsplit(args[3], ",", fixed = TRUE)[[1]])
output_dir <- if (length(args) >= 4) normalizePath(args[4], winslash = "\\", mustWork = FALSE) else {
  file.path(project_root, "outputs", "custom_run")
}

dat <- utils::read.csv(studies_csv, stringsAsFactors = FALSE)
targets <- utils::read.csv(targets_csv, stringsAsFactors = FALSE)

result <- fit_transport_engine(
  dat = dat,
  targets = targets,
  moderators = moderators,
  yi_col = "yi",
  vi_col = "vi",
  study_col = if ("study" %in% names(dat)) "study" else names(dat)[1],
  target_label_col = if ("label" %in% names(targets)) "label" else names(targets)[1],
  effect_scale = "identity",
  rob_col = if ("risk_of_bias" %in% names(dat)) "risk_of_bias" else NULL
)

save_transport_results(result, output_dir)

cat("Custom run complete.\n")
cat(sprintf("Outputs written to: %s\n", output_dir))
