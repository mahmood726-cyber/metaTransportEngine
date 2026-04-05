test_that("bundled example files are available", {
  expect_true(file.exists(example_studies_file()))
  expect_true(file.exists(example_targets_file()))
  expect_true(file.exists(example_truthcert_schema_file()))
})

test_that("frequentist engine fits and exports results", {
  dat <- read_example_studies()
  targets <- read_example_targets()

  result <- fit_transport_engine(
    dat = dat,
    targets = targets,
    moderators = c("mean_age", "female_pct", "baseline_risk", "followup_months"),
    effect_scale = "log_rr",
    rob_col = "risk_of_bias"
  )

  expect_s3_class(result, "meta_transport_engine")
  expect_equal(nrow(result$target_estimates), nrow(targets))
  expect_true(all(c("transported_mean", "lci", "uci", "lpi", "upi") %in% names(result$target_estimates)))
  expect_true(is.finite(result$pooled_mean))
  expect_s3_class(plot_transport_forest(result), "ggplot")

  output_dir <- file.path(tempdir(), "metaTransportEngine-frequentist")
  save_transport_results(result, output_dir)

  expect_true(file.exists(file.path(output_dir, "target_estimates.csv")))
  expect_true(file.exists(file.path(output_dir, "transport_forest_plot.png")))
  expect_true(file.exists(file.path(output_dir, "summary.json")))
})
