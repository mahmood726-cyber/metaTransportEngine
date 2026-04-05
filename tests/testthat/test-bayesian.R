test_that("bayesian engine fits a minimal smoke-test model", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("rstan")

  dat <- read_example_studies()
  targets <- read_example_targets()

  result <- fit_bayesian_transport(
    dat = dat,
    targets = targets,
    moderators = c("mean_age", "female_pct", "baseline_risk", "followup_months"),
    effect_scale = "log_rr",
    chains = 1,
    iter = 200,
    warmup = 100,
    refresh = 0,
    seed = 1
  )

  expect_s3_class(result, "bayesian_meta_transport")
  expect_equal(nrow(result$target_posterior), nrow(targets))
  expect_true(all(c("mean", "lci", "uci") %in% names(result$target_posterior)))

  output_dir <- file.path(tempdir(), "metaTransportEngine-bayesian")
  save_bayesian_results(result, output_dir)

  expect_true(file.exists(file.path(output_dir, "target_posterior.csv")))
  expect_true(file.exists(file.path(output_dir, "tau_summary.csv")))
})
