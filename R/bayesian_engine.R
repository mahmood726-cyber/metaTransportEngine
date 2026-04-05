configure_windows_toolchain <- function() {
  if (.Platform$OS.type != "windows") {
    return(invisible(NULL))
  }

  if (nzchar(Sys.which("make")) && nzchar(Sys.which("g++"))) {
    return(invisible(NULL))
  }

  candidate_roots <- c(
    "C:/RBuildTools/4.5",
    "C:/RBuildTools/4.4",
    "C:/rtools45",
    "C:/rtools44",
    "C:/rtools43"
  )

  path_sep <- .Platform$path.sep
  current_path <- strsplit(Sys.getenv("PATH"), path_sep, fixed = TRUE)[[1]]

  for (root in candidate_roots) {
    tool_paths <- c(file.path(root, "usr", "bin"), file.path(root, "ucrt64", "bin"))

    if (all(dir.exists(tool_paths))) {
      updated_path <- paste(unique(c(tool_paths, current_path)), collapse = path_sep)
      Sys.setenv(PATH = updated_path)

      if (nzchar(Sys.which("make")) && nzchar(Sys.which("g++"))) {
        return(invisible(tool_paths))
      }
    }
  }

  invisible(NULL)
}

resolve_stan_file <- function(model_name = c("meta_transport", "meta_transport_bias"), stan_file = NULL) {
  model_name <- match.arg(model_name)
  stan_file <- stan_file %||% resolve_inst_file("stan", paste0(model_name, ".stan"))

  if (!file.exists(stan_file)) {
    stop(
      sprintf("Stan file not found: %s.", stan_file),
      call. = FALSE
    )
  }

  stan_file
}

make_stan_transport_data <- function(dat,
                                     targets,
                                     moderators,
                                     yi_col = "yi",
                                     vi_col = "vi",
                                     bias_score_col = NULL) {
  assert_required_columns(dat, c(yi_col, vi_col, moderators), "dat")
  assert_required_columns(targets, moderators, "targets")
  assert_numeric_columns(dat, c(yi_col, vi_col, moderators), "dat")
  assert_numeric_columns(targets, moderators, "targets")

  stan_data <- list(
    S = nrow(dat),
    P = length(moderators),
    T = nrow(targets),
    y = dat[[yi_col]],
    v = dat[[vi_col]],
    Z = as.matrix(dat[, moderators, drop = FALSE]),
    ZT = as.matrix(targets[, moderators, drop = FALSE])
  )

  if (!is.null(bias_score_col)) {
    assert_required_columns(dat, bias_score_col, "dat")
    assert_numeric_columns(dat, bias_score_col, "dat")
    stan_data$q <- dat[[bias_score_col]]
  }

  stan_data
}

summarise_parameter_draws <- function(draws_matrix, term_names) {
  if (is.list(draws_matrix) && !is.data.frame(draws_matrix)) {
    draws_matrix <- do.call(cbind, draws_matrix)
  }

  draws_matrix <- tryCatch(
    as.matrix(draws_matrix),
    error = function(e) matrix(as.numeric(draws_matrix), ncol = 1)
  )

  if (is.null(dim(draws_matrix)) || length(dim(draws_matrix)) < 2) {
    draws_matrix <- matrix(as.numeric(draws_matrix), ncol = 1)
  }

  if (ncol(draws_matrix) != length(term_names)) {
    if (length(term_names) == 1) {
      draws_matrix <- matrix(as.numeric(draws_matrix), ncol = 1)
    } else {
      stop("Parameter draw columns do not match the supplied term names.", call. = FALSE)
    }
  }

  rows <- lapply(seq_len(ncol(draws_matrix)), function(i) {
    stats_i <- summarise_draws_vector(draws_matrix[, i])
    data.frame(
      term = term_names[i],
      mean = unname(stats_i["mean"]),
      median = unname(stats_i["median"]),
      sd = unname(stats_i["sd"]),
      lci = unname(stats_i["lci"]),
      uci = unname(stats_i["uci"]),
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}

#' Fit a Bayesian meta-transport model
#'
#' Fits the hierarchical transport model in Stan and returns posterior and
#' predictive summaries for one or more target populations.
#'
#' @param dat Study-level data frame.
#' @param targets Target-population data frame.
#' @param moderators Character vector of moderator column names.
#' @param yi_col Name of the study effect estimate column.
#' @param vi_col Name of the within-study variance column.
#' @param target_label_col Name of the target label column.
#' @param effect_scale Effect scale used for optional back-transformation.
#' @param bias_score_col Optional numeric bias score column for the extended
#'   bias model.
#' @param stan_file Optional path to a Stan file. By default the bundled model
#'   is used.
#' @param chains Number of MCMC chains.
#' @param iter Total iterations per chain.
#' @param warmup Warmup iterations per chain.
#' @param seed Random seed passed to Stan.
#' @param refresh Stan progress refresh interval.
#' @param adapt_delta Stan `adapt_delta` control parameter.
#' @param max_treedepth Stan `max_treedepth` control parameter.
#' @param auto_write Whether `rstan` should cache compiled model objects next
#'   to the Stan file.
#'
#' @return An object of class `"bayesian_meta_transport"`.
#' @export
fit_bayesian_transport <- function(dat,
                                   targets,
                                   moderators,
                                   yi_col = "yi",
                                   vi_col = "vi",
                                   target_label_col = "label",
                                   effect_scale = c("identity", "log_rr", "log_or", "log_hr"),
                                   bias_score_col = NULL,
                                   stan_file = NULL,
                                   chains = 4,
                                   iter = 2000,
                                   warmup = floor(iter / 2),
                                   seed = 1234,
                                   refresh = 0,
                                   adapt_delta = 0.95,
                                   max_treedepth = 12,
                                   auto_write = FALSE) {
  effect_scale <- match.arg(effect_scale)

  if (!target_label_col %in% names(targets)) {
    targets[[target_label_col]] <- paste0("target_", seq_len(nrow(targets)))
  }

  warn_if_thin_data(dat, moderators)
  configure_windows_toolchain()

  model_name <- if (is.null(bias_score_col)) "meta_transport" else "meta_transport_bias"
  stan_file <- resolve_stan_file(model_name, stan_file)
  stan_data <- make_stan_transport_data(dat, targets, moderators, yi_col, vi_col, bias_score_col)

  rstan::rstan_options(auto_write = auto_write)

  sm <- rstan::stan_model(file = stan_file)
  fit <- rstan::sampling(
    object = sm,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    refresh = refresh,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )

  draws <- rstan::extract(fit)

  target_posterior <- summarise_draws_matrix(
    draws$theta_T,
    labels = targets[[target_label_col]],
    effect_scale = effect_scale
  )
  names(target_posterior)[1] <- "label"

  target_predictive <- summarise_draws_matrix(
    draws$theta_T_pred,
    labels = targets[[target_label_col]],
    effect_scale = effect_scale
  )
  names(target_predictive)[1] <- "label"

  modifier_summary <- summarise_parameter_draws(
    cbind(draws$beta0, draws$beta),
    c("(Intercept)", moderators)
  )

  tau_summary <- summarise_parameter_draws(draws$tau, "tau")

  bias_summary <- NULL
  if (!is.null(bias_score_col)) {
    bias_summary <- list(
      delta = summarise_parameter_draws(draws$delta, "delta"),
      sigma_b = summarise_parameter_draws(draws$sigma_b, "sigma_b")
    )
  }

  result <- list(
    fit = fit,
    dat = dat,
    targets = targets,
    spec = list(
      moderators = moderators,
      yi_col = yi_col,
      vi_col = vi_col,
      target_label_col = target_label_col,
      effect_scale = effect_scale,
      bias_score_col = bias_score_col,
      stan_file = stan_file
    ),
    target_posterior = target_posterior,
    target_predictive = target_predictive,
    modifier_summary = modifier_summary,
    tau_summary = tau_summary,
    bias_summary = bias_summary
  )

  class(result) <- "bayesian_meta_transport"
  result
}

#' Save Bayesian transport outputs
#'
#' Writes posterior and predictive summaries to `output_dir`.
#'
#' @param result A fitted `"bayesian_meta_transport"` object.
#' @param output_dir Output directory path.
#'
#' @return Invisibly returns `output_dir`.
#' @export
save_bayesian_results <- function(result, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  utils::write.csv(result$target_posterior, file.path(output_dir, "target_posterior.csv"), row.names = FALSE)
  utils::write.csv(result$target_predictive, file.path(output_dir, "target_predictive.csv"), row.names = FALSE)
  utils::write.csv(result$modifier_summary, file.path(output_dir, "modifier_summary.csv"), row.names = FALSE)
  utils::write.csv(result$tau_summary, file.path(output_dir, "tau_summary.csv"), row.names = FALSE)

  if (!is.null(result$bias_summary)) {
    utils::write.csv(result$bias_summary$delta, file.path(output_dir, "bias_delta_summary.csv"), row.names = FALSE)
    utils::write.csv(result$bias_summary$sigma_b, file.path(output_dir, "bias_sigma_summary.csv"), row.names = FALSE)
  }

  invisible(output_dir)
}

#' Print a Bayesian meta-transport fit
#'
#' @param x A fitted `"bayesian_meta_transport"` object.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.bayesian_meta_transport <- function(x, ...) {
  cat("Bayesian Meta-Transport Model\n")
  cat(sprintf("Targets: %d\n", nrow(x$target_posterior)))
  cat(sprintf("Moderators: %s\n", paste(x$spec$moderators, collapse = ", ")))
  cat(sprintf("Stan file: %s\n", x$spec$stan_file))
  invisible(x)
}
