# ============================================================
# harness-guard.R -- before/after for the extrapolation guardrail.
#
# The validation harness (harness.R) measured that under model misspecification,
# extrapolating the transport beyond the studies' modifier support collapses
# coverage of the true effect (-> 0.00) with NO widening of the CI: it fails
# SILENTLY. This harness measures whether the new guardrail (study_data flag on
# predict_transport_target) RAISES a flag on exactly those dangerous cells while
# staying quiet on safe interpolation.
#
#   BEFORE: study_data = NULL  -> no flag column exists (silent, flag rate 0).
#   AFTER : study_data = dat   -> extrapolation_flag emitted; we measure its rate.
#
# Truth-first: coverage numbers are re-measured from seeded sims and reported
# next to the flag rate, so the flag's value is judged against actual coverage.
# Run:  Rscript truth-recovery/harness-guard.R 600
# ============================================================
suppressMessages(library(metafor))
this <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
rdir <- file.path(dirname(this), "..", "R")
source(file.path(rdir, "common.R")); source(file.path(rdir, "frequentist_engine.R"))

args <- commandArgs(TRUE)
NSIM <- if (length(args) >= 1) as.integer(args[1]) else 600
BETA0 <- 0.2; BETA1 <- 0.5; TAU2 <- 0.02
true_effect <- function(xt, curve) BETA0 + BETA1 * xt + curve * (xt^2)

gen <- function(k, curve, seed) {
  set.seed(seed)
  x <- runif(k, 0, 1)
  theta <- BETA0 + BETA1 * x + curve * x^2 + rnorm(k, 0, sqrt(TAU2))
  vi <- runif(k, 0.01, 0.06)
  data.frame(study = paste0("s", seq_len(k)), yi = theta + rnorm(k, 0, sqrt(vi)), vi = vi, x = x)
}

run_cell <- function(k, xt, curve) {
  cov <- c(); flag_before <- c(); flag_after <- c()
  truth <- true_effect(xt, curve)
  for (s in 1:NSIM) {
    dat <- gen(k, curve, 5000 + s)
    fit <- tryCatch(rma.uni(yi = dat$yi, vi = dat$vi, mods = ~ x, data = dat, method = "REML"),
                    error = function(e) NULL)
    if (is.null(fit)) next
    tr <- data.frame(x = xt)
    # BEFORE: no study_data -> no guard column
    pb <- tryCatch(predict_transport_target(fit, tr, "x", "t"), error = function(e) NULL)
    if (is.null(pb)) next
    # AFTER: study_data supplied -> guard column present
    pa <- tryCatch(suppressWarnings(predict_transport_target(fit, tr, "x", "t", study_data = dat)),
                   error = function(e) NULL)
    if (is.null(pa)) next
    cov <- c(cov, pb$lci <= truth & truth <= pb$uci)
    flag_before <- c(flag_before, !is.null(pb$extrapolation_flag) && isTRUE(pb$extrapolation_flag))
    flag_after  <- c(flag_after, isTRUE(pa$extrapolation_flag))
  }
  list(cov = round(mean(cov), 3),
       fb = round(mean(flag_before), 3),
       fa = round(mean(flag_after), 3),
       n = length(cov))
}

cat(sprintf("\n# Extrapolation-guard before/after -- metaTransportEngine  nsim=%d\n", NSIM))
cat("studies' modifier x in [0,1]; true effect = 0.2 + 0.5*x (+ curve*x^2)\n")
cat("flag_before = guard OFF (study_data=NULL); flag_after = guard ON (study_data=dat)\n\n")
cat(sprintf("%-34s | %8s | %10s | %9s\n", "scenario", "coverage", "flag_before", "flag_after"))
cells <- list(
  list(l="interpolation, linear (xt=0.5)",      xt=0.5, c=0),
  list(l="extrapolation, linear (xt=2.0)",      xt=2.0, c=0),
  list(l="extrapolation, linear (xt=3.0)",      xt=3.0, c=0),
  list(l="interpolation, NONLINEAR (xt=0.5)",   xt=0.5, c=0.6),
  list(l="extrapolation, NONLINEAR (xt=2.0)",   xt=2.0, c=0.6),
  list(l="extrapolation, NONLINEAR (xt=3.0)",   xt=3.0, c=0.6)
)
for (cl in cells) {
  r <- run_cell(12, cl$xt, cl$c)
  cat(sprintf("%-34s | %8.3f | %10.3f | %9.3f\n", cl$l, r$cov, r$fb, r$fa))
}
cat("\nGuard works if: extrapolation rows (esp. the coverage~0 NONLINEAR cells) get flag_after~1.0,\n")
cat("while interpolation rows stay near 0 (low false-alarm). flag_before is 0 everywhere (silent).\n")
