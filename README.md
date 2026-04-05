# Meta-Transport Engine

[![R-CMD-check](https://github.com/mahmood726-cyber/metaTransportEngine/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mahmood726-cyber/metaTransportEngine/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/mahmood726-cyber/metaTransportEngine/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/mahmood726-cyber/metaTransportEngine/actions/workflows/pkgdown.yaml)
[![site](https://img.shields.io/badge/docs-pkgdown-1f6feb)](https://mahmood726-cyber.github.io/metaTransportEngine/)

Standalone R package project for aggregate-data meta-analysis with transport to one or more target populations.

The project implements three layers around one common model:

1. Frequentist random-effects meta-regression with `metafor`
2. Bayesian hierarchical transport model in Stan via `rstan`
3. Validators, schema examples, and export/plot helpers

## Core model

For study `s = 1, ..., S`:

- `y_s | theta_s ~ Normal(theta_s, v_s)`
- `theta_s = beta_0 + z_s^T beta + u_s`
- `u_s ~ Normal(0, tau^2)`

For a target modifier profile `z_T`:

- `theta_T = beta_0 + z_T^T beta`

The engine reports:

- the transported mean effect `theta_T`
- a confidence interval for the transported mean
- a prediction interval for a realized target effect using `tau^2`

## Project layout

- `R/common.R`: shared helpers
- `R/frequentist_engine.R`: frequentist engine, validators, plot/export helpers
- `R/bayesian_engine.R`: Bayesian wrapper functions
- `inst/stan/`: Stan models
- `inst/extdata/`: example study CSV, target CSV, and TruthCert-style YAML schema
- `scripts/`: runnable entry points
- `tests/testthat/`: package tests

## Package workflow

Build and install locally:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\R.exe' CMD build .
& 'C:\Program Files\R\R-4.5.2\bin\R.exe' CMD INSTALL metaTransportEngine_0.1.0.tar.gz
```

Run package checks:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\R.exe' CMD check --no-manual metaTransportEngine_0.1.0.tar.gz
```

Build the pkgdown site locally:

```powershell
$env:RSTUDIO_PANDOC = 'C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools'
$env:PATH = $env:RSTUDIO_PANDOC + ';' + $env:PATH
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' -e "pkgdown::build_site('.')"
```

Browse the installed vignette:

```r
browseVignettes("metaTransportEngine")
```

## Quick start

Because `Rscript` is not on `PATH` here, use the full executable:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' scripts/run_demo.R
```

That loads the local package with `pkgload`, reads the bundled example files from `inst/extdata`, writes CSV/JSON/PNG outputs into `outputs/demo_frequentist`, and prints the transported target estimates.

To fit your own files:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' scripts/run_from_csv.R `
  inst/extdata/studies_example.csv `
  inst/extdata/targets_example.csv `
  mean_age,female_pct,baseline_risk,followup_months `
  outputs/custom_run
```

## Bayesian demo

Run the Stan layer separately because compilation and sampling take longer:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' scripts/run_bayesian_demo.R
```

That writes posterior summaries into `outputs/demo_bayesian`. The Bayesian path auto-detects a local Windows R toolchain when `make` is not already on `PATH`.

## Input contract

Study CSV must include:

- a study label column, default `study`
- `yi`: study effect estimate on the analysis scale
- `vi`: within-study sampling variance
- one or more numeric moderator columns

Optional study columns:

- `risk_of_bias`: categorical bias label for sensitivity analysis
- `bias_score`: numeric bias score for the Bayesian bias model

Target CSV must include:

- a target label column, default `label`
- the same moderator columns used in the model

## Notes

- If the analysis scale is `log_rr`, `log_or`, or `log_hr`, the engine also reports exponentiated summaries.
- The engine checks whether the target sits outside the observed study support and computes Mahalanobis-based extrapolation scores.
- With very few studies relative to the number of moderators, meta-regression becomes unstable. The engine warns when the study count is too small for the requested moderator set.
- The package includes a getting-started vignette in `vignettes/getting-started.Rmd`.
- `_pkgdown.yml` and `.github/workflows/pkgdown.yaml` are included so the package can publish a documentation site once it is pushed to GitHub.
