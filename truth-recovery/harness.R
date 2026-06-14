# ============================================================
# harness.R -- Truth-recovery yardstick for metaTransportEngine.
#
# The engine transports a meta-analytic effect to a TARGET population by
# meta-regression on effect modifiers and projecting to the target's modifier
# value (predict_transport_target). The honest test: inject a KNOWN effect-
# modification truth, set a target modifier value, and check whether the
# transported CI covers the TRUE target-population effect under
#   (a) INTERPOLATION (target within the studies' modifier range),
#   (b) EXTRAPOLATION (target outside the range), and
#   (c) MODEL MISSPECIFICATION (true modification is NON-linear, engine fits
#       linear) -- the regime where transport is most dangerous.
#
# Uses the app's own predict_transport_target on metafor rma.uni fits, unchanged.
# Truth-first: every number is produced from seeded simulation. R via PowerShell.
# Run:  Rscript truth-recovery/harness.R 600
# ============================================================
suppressMessages(library(metafor))
this <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
rdir <- file.path(dirname(this), "..", "R")
source(file.path(rdir, "common.R")); source(file.path(rdir, "frequentist_engine.R"))

args <- commandArgs(TRUE)
NSIM <- if (length(args) >= 1) as.integer(args[1]) else 600

BETA0 <- 0.2; BETA1 <- 0.5; TAU2 <- 0.02

# true target effect at modifier x_t (linear or +quadratic curvature)
true_effect <- function(xt, curve) BETA0 + BETA1 * xt + curve * (xt^2)

gen <- function(k, xlo, xhi, curve, seed) {
  set.seed(seed)
  x <- runif(k, xlo, xhi)
  theta <- BETA0 + BETA1 * x + curve * (x^2) + rnorm(k, 0, sqrt(TAU2))
  vi <- runif(k, 0.01, 0.06)
  yi <- theta + rnorm(k, 0, sqrt(vi))
  data.frame(study = paste0("s", seq_len(k)), yi = yi, vi = vi, x = x)
}

run_cell <- function(k, xlo, xhi, xt, curve) {
  cover <- c(); width <- c()
  truth <- true_effect(xt, curve)
  for (s in 1:NSIM) {
    dat <- gen(k, xlo, xhi, curve, 5000 + s)
    fit <- tryCatch(rma.uni(yi = dat$yi, vi = dat$vi, mods = ~ x, data = dat, method = "REML"),
                    error = function(e) NULL)
    if (is.null(fit)) next
    target <- data.frame(x = xt)
    pr <- tryCatch(predict_transport_target(fit, target, "x", "t"), error = function(e) NULL)
    if (is.null(pr)) next
    cover <- c(cover, pr$lci <= truth & truth <= pr$uci)
    width <- c(width, pr$uci - pr$lci)
  }
  list(cov = round(mean(cover), 3), width = round(mean(width), 3), n = length(cover))
}

cat(sprintf("\n# Truth-recovery yardstick -- metaTransportEngine  nsim=%d\n", NSIM))
cat("studies' modifier x in [0,1]; true effect = 0.2 + 0.5*x (+ curve*x^2)\n\n")
cat(sprintf("%-34s | %8s | %8s\n", "scenario", "coverage", "CI width"))
cells <- list(
  list(l="interpolation, linear (xt=0.5)",      k=12, lo=0, hi=1, xt=0.5, c=0),
  list(l="extrapolation, linear (xt=2.0)",      k=12, lo=0, hi=1, xt=2.0, c=0),
  list(l="extrapolation, linear (xt=3.0)",      k=12, lo=0, hi=1, xt=3.0, c=0),
  list(l="interpolation, NONLINEAR (xt=0.5)",   k=12, lo=0, hi=1, xt=0.5, c=0.6),
  list(l="extrapolation, NONLINEAR (xt=2.0)",   k=12, lo=0, hi=1, xt=2.0, c=0.6),
  list(l="extrapolation, NONLINEAR (xt=3.0)",   k=12, lo=0, hi=1, xt=3.0, c=0.6)
)
for (cl in cells) {
  r <- run_cell(cl$k, cl$lo, cl$hi, cl$xt, cl$c)
  cat(sprintf("%-34s | %8.3f | %8.3f\n", cl$l, r$cov, r$width))
}
cat("\n(coverage of the TRUE target-population effect; should be ~0.95.\n")
cat(" NONLINEAR rows = engine fits LINEAR but truth is quadratic = model misspecification.)\n")
