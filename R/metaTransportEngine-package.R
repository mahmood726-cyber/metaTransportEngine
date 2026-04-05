#' metaTransportEngine: Meta-Analysis Transport Models for Target Populations
#'
#' `metaTransportEngine` fits aggregate-data meta-regression models that
#' transport evidence from study populations to one or more target populations.
#' It includes a frequentist `metafor` engine, a Bayesian Stan engine via
#' `rstan`, and diagnostics for overlap, extrapolation, instability, and bias
#' sensitivity.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("estimate", "label", "lcl", "lpi", "type", "ucl", "upi"))
