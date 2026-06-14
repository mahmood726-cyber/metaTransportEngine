# Rscript truth-recovery/test-truth-recovery.R   (exit 0 = all pass)
# Measured invariants for the metaTransportEngine yardstick. Seeded.
suppressMessages(library(metafor))
this <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
rdir <- file.path(dirname(this), "..", "R")
source(file.path(rdir, "common.R")); source(file.path(rdir, "frequentist_engine.R"))

NSIM <- 400; BETA0 <- 0.2; BETA1 <- 0.5; TAU2 <- 0.02
true_effect <- function(xt, curve) BETA0 + BETA1 * xt + curve * (xt^2)
gen <- function(k, curve, seed) { set.seed(seed); x <- runif(k, 0, 1)
  theta <- BETA0 + BETA1 * x + curve * x^2 + rnorm(k, 0, sqrt(TAU2)); vi <- runif(k, 0.01, 0.06)
  data.frame(study = paste0("s", 1:k), yi = theta + rnorm(k, 0, sqrt(vi)), vi = vi, x = x) }
cell <- function(k, xt, curve) {
  cov <- c(); w <- c(); truth <- true_effect(xt, curve)
  for (s in 1:NSIM) {
    dat <- gen(k, curve, 5000 + s)
    fit <- tryCatch(rma.uni(yi = dat$yi, vi = dat$vi, mods = ~ x, data = dat, method = "REML"), error = function(e) NULL)
    if (is.null(fit)) next
    pr <- tryCatch(predict_transport_target(fit, data.frame(x = xt), "x", "t"), error = function(e) NULL)
    if (is.null(pr)) next
    cov <- c(cov, pr$lci <= truth & truth <= pr$uci); w <- c(w, pr$uci - pr$lci)
  }
  list(cov = mean(cov), w = mean(w))
}

ok <- TRUE
report <- function(name, cond, detail) { cat(sprintf("%-4s %s  %s\n", if (cond) "PASS" else "FAIL", name, detail)); if (!cond) ok <<- FALSE }

li <- cell(12, 0.5, 0)     # linear interpolation
le <- cell(12, 3.0, 0)     # linear extrapolation
report("VALIDATION: correctly-specified linear transport covers the truth (interpolation)", li$cov > 0.92, sprintf("(%.3f)", li$cov))
report("CI widens appropriately on extrapolation (variance propagation works)", le$w > 4 * li$w, sprintf("(width %.2f extrap vs %.2f interp)", le$w, li$w))
report("linear extrapolation keeps reasonable coverage", le$cov > 0.85, sprintf("(%.3f)", le$cov))

ne <- cell(12, 3.0, 0.6)   # NONLINEAR truth, engine fits linear, extrapolated
report("CRITICAL: under misspecification, extrapolation coverage COLLAPSES", ne$cov < 0.2, sprintf("(%.3f -- confident but wrong)", ne$cov))
report("the CI does NOT widen to reflect misspecification (fails silently)", abs(ne$w - le$w) < 0.3 * le$w, sprintf("(misspec width %.2f vs linear %.2f)", ne$w, le$w))

cat(if (ok) "\nAll measured invariants hold.\n" else "\nSOME INVARIANTS FAILED.\n")
quit(status = if (ok) 0 else 1)
