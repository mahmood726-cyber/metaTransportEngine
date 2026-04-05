build_modifier_formula <- function(moderators) {
  if (length(moderators) == 0) {
    stats::as.formula("~ 1")
  } else {
    stats::as.formula(paste("~", paste(moderators, collapse = " + ")))
  }
}

build_design_matrix <- function(data, moderators) {
  formula <- build_modifier_formula(moderators)
  stats::model.matrix(formula, data = data)
}

warn_if_thin_data <- function(dat, moderators) {
  k <- nrow(dat)
  p <- length(moderators)

  if (k <= p + 1) {
    stop(
      sprintf(
        "Not enough studies for meta-regression: %d studies for %d moderators plus an intercept.",
        k,
        p
      ),
      call. = FALSE
    )
  }

  if (k <= p + 4) {
    warning(
      sprintf(
        "Only %d studies are available for %d moderators. Estimates may be unstable.",
        k,
        p
      ),
      call. = FALSE
    )
  }
}

summarise_coefficients <- function(fit) {
  beta_hat <- as.numeric(stats::coef(fit))
  names(beta_hat) <- names(stats::coef(fit))
  se_hat <- sqrt(diag(stats::vcov(fit)))

  data.frame(
    term = names(beta_hat),
    estimate = beta_hat,
    se = se_hat,
    lci = beta_hat - 1.96 * se_hat,
    uci = beta_hat + 1.96 * se_hat,
    row.names = NULL
  )
}

predict_transport_target <- function(fit, target_row, moderators, target_label, effect_scale = "identity") {
  x_target <- as.numeric(build_design_matrix(target_row, moderators)[1, , drop = TRUE])
  beta_hat <- as.numeric(stats::coef(fit))
  v_beta <- stats::vcov(fit)
  tau2_hat <- fit$tau2 %||% 0

  mu_target <- sum(x_target * beta_hat)
  var_mu_target <- as.numeric(t(x_target) %*% v_beta %*% x_target)
  var_mu_target <- max(var_mu_target, 0)

  se_mu_target <- sqrt(var_mu_target)
  se_pred_target <- sqrt(var_mu_target + tau2_hat)

  out <- data.frame(
    label = target_label,
    transported_mean = mu_target,
    lci = mu_target - 1.96 * se_mu_target,
    uci = mu_target + 1.96 * se_mu_target,
    lpi = mu_target - 1.96 * se_pred_target,
    upi = mu_target + 1.96 * se_pred_target,
    row.names = NULL
  )

  if (effect_scale != "identity") {
    out$transported_mean_bt <- back_transform_values(out$transported_mean, effect_scale)
    out$lci_bt <- back_transform_values(out$lci, effect_scale)
    out$uci_bt <- back_transform_values(out$uci, effect_scale)
    out$lpi_bt <- back_transform_values(out$lpi, effect_scale)
    out$upi_bt <- back_transform_values(out$upi, effect_scale)
  }

  out
}

check_modifier_overlap <- function(dat, targets, moderators, target_label_col = "label") {
  rows <- lapply(seq_len(nrow(targets)), function(i) {
    target_row <- targets[i, , drop = FALSE]

    do.call(
      rbind,
      lapply(moderators, function(mod) {
        min_val <- min(dat[[mod]], na.rm = TRUE)
        max_val <- max(dat[[mod]], na.rm = TRUE)
        target_val <- target_row[[mod]][1]
        outside <- target_val < min_val || target_val > max_val

        data.frame(
          label = target_row[[target_label_col]][1],
          modifier = mod,
          target_value = target_val,
          study_min = min_val,
          study_max = max_val,
          outside_support = outside,
          distance_to_support = if (outside) {
            min(abs(target_val - min_val), abs(target_val - max_val))
          } else {
            0
          },
          row.names = NULL
        )
      })
    )
  })

  do.call(rbind, rows)
}

safe_cov_inverse <- function(x) {
  cov_mat <- stats::cov(x)

  if (is.null(dim(cov_mat))) {
    cov_mat <- matrix(cov_mat, nrow = 1, ncol = 1)
  }

  ridge_scale <- mean(abs(diag(cov_mat)))
  if (!is.finite(ridge_scale) || ridge_scale == 0) {
    ridge_scale <- 1
  }

  cov_mat + diag(ridge_scale * 1e-8, nrow(cov_mat))
}

compute_extrapolation_scores <- function(dat, targets, moderators, target_label_col = "label") {
  z_study <- as.matrix(dat[, moderators, drop = FALSE])
  center <- colMeans(z_study)
  cov_regularised <- safe_cov_inverse(z_study)
  inv_cov <- qr.solve(cov_regularised)
  p <- length(moderators)

  rows <- lapply(seq_len(nrow(targets)), function(i) {
    z_target <- as.numeric(targets[i, moderators, drop = FALSE])
    delta <- z_target - center
    d2 <- as.numeric(t(delta) %*% inv_cov %*% delta)
    q95 <- stats::qchisq(0.95, df = p)
    q99 <- stats::qchisq(0.99, df = p)

    data.frame(
      label = targets[[target_label_col]][i],
      mahalanobis_d2 = d2,
      mahalanobis_d = sqrt(max(d2, 0)),
      threshold_95 = sqrt(q95),
      threshold_99 = sqrt(q99),
      extrapolation_flag = d2 > q95,
      extrapolation_severity = if (d2 > q99) {
        "high"
      } else if (d2 > q95) {
        "moderate"
      } else {
        "low"
      },
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}

compute_transport_gaps <- function(fit, dat, target_estimates, moderators, study_col = "study") {
  x_study <- build_design_matrix(dat, moderators)
  beta_hat <- as.numeric(stats::coef(fit))
  study_pred <- as.numeric(x_study %*% beta_hat)

  rows <- lapply(seq_len(nrow(target_estimates)), function(i) {
    data.frame(
      label = target_estimates$label[i],
      study = dat[[study_col]],
      predicted_study_effect = study_pred,
      transported_mean = target_estimates$transported_mean[i],
      transport_gap = target_estimates$transported_mean[i] - study_pred,
      absolute_gap = abs(target_estimates$transported_mean[i] - study_pred),
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}

leave_one_out_transport <- function(dat, moderators, yi_col = "yi", vi_col = "vi", study_col = "study", method = "REML") {
  full_formula <- build_modifier_formula(moderators)
  full_fit <- metafor::rma.uni(
    yi = dat[[yi_col]],
    vi = dat[[vi_col]],
    mods = full_formula,
    method = method,
    data = dat
  )

  full_beta <- stats::coef(full_fit)
  full_terms <- names(full_beta)

  loo_rows <- lapply(seq_len(nrow(dat)), function(i) {
    reduced <- dat[-i, , drop = FALSE]

    if (nrow(reduced) <= length(moderators) + 1) {
      return(
        data.frame(
          omitted_study = dat[[study_col]][i],
          term = full_terms,
          estimate = NA_real_,
          delta_from_full = NA_real_,
          abs_delta = NA_real_,
          note = "insufficient_studies_after_omission",
          row.names = NULL
        )
      )
    }

    fit_i <- try(
      metafor::rma.uni(
        yi = reduced[[yi_col]],
        vi = reduced[[vi_col]],
        mods = full_formula,
        method = method,
        data = reduced
      ),
      silent = TRUE
    )

    if (inherits(fit_i, "try-error")) {
      return(
        data.frame(
          omitted_study = dat[[study_col]][i],
          term = full_terms,
          estimate = NA_real_,
          delta_from_full = NA_real_,
          abs_delta = NA_real_,
          note = "fit_failed",
          row.names = NULL
        )
      )
    }

    beta_i <- rep(NA_real_, length(full_terms))
    names(beta_i) <- full_terms
    beta_i[names(stats::coef(fit_i))] <- as.numeric(stats::coef(fit_i))

    data.frame(
      omitted_study = dat[[study_col]][i],
      term = full_terms,
      estimate = beta_i,
      delta_from_full = beta_i - as.numeric(full_beta),
      abs_delta = abs(beta_i - as.numeric(full_beta)),
      note = "ok",
      row.names = NULL
    )
  })

  loo_detail <- do.call(rbind, loo_rows)

  loo_summary <- do.call(
    rbind,
    lapply(split(loo_detail, loo_detail$term), function(df_term) {
      valid <- df_term$abs_delta[is.finite(df_term$abs_delta)]
      data.frame(
        term = df_term$term[1],
        max_abs_delta = if (length(valid) == 0) NA_real_ else max(valid),
        median_abs_delta = if (length(valid) == 0) NA_real_ else stats::median(valid),
        fits_failed = sum(df_term$note != "ok"),
        row.names = NULL
      )
    })
  )

  list(detail = loo_detail, summary = loo_summary)
}

risk_of_bias_sensitivity <- function(dat,
                                     targets,
                                     moderators,
                                     yi_col = "yi",
                                     vi_col = "vi",
                                     target_label_col = "label",
                                     rob_col = NULL,
                                     high_risk_values = c("high"),
                                     method = "REML",
                                     effect_scale = "identity") {
  if (is.null(rob_col) || !rob_col %in% names(dat)) {
    return(NULL)
  }

  keep <- !(normalise_risk_values(dat[[rob_col]]) %in% normalise_risk_values(high_risk_values))

  if (sum(keep) <= length(moderators) + 1) {
    return(
      data.frame(
        note = "Not enough studies remain after excluding high-risk studies.",
        row.names = NULL
      )
    )
  }

  filtered <- dat[keep, , drop = FALSE]
  fit <- metafor::rma.uni(
    yi = filtered[[yi_col]],
    vi = filtered[[vi_col]],
    mods = build_modifier_formula(moderators),
    method = method,
    data = filtered
  )

  estimates <- do.call(
    rbind,
    lapply(seq_len(nrow(targets)), function(i) {
      predict_transport_target(
        fit = fit,
        target_row = targets[i, , drop = FALSE],
        moderators = moderators,
        target_label = targets[[target_label_col]][i],
        effect_scale = effect_scale
      )
    })
  )

  estimates$analysis <- "exclude_high_risk"
  estimates$studies_remaining <- nrow(filtered)
  estimates
}

#' Plot transported effect estimates
#'
#' Creates a forest-style plot for study estimates and transported target
#' estimates, including target prediction intervals.
#'
#' @param result A fitted `"meta_transport_engine"` object.
#' @param effect_label Axis label for the effect scale.
#'
#' @return A `ggplot` object.
#' @export
plot_transport_forest <- function(result, effect_label = "Effect") {
  studies_df <- data.frame(
    label = result$dat[[result$spec$study_col]],
    estimate = result$dat[[result$spec$yi_col]],
    lcl = result$dat[[result$spec$yi_col]] - 1.96 * sqrt(result$dat[[result$spec$vi_col]]),
    ucl = result$dat[[result$spec$yi_col]] + 1.96 * sqrt(result$dat[[result$spec$vi_col]]),
    lpi = NA_real_,
    upi = NA_real_,
    type = "Study",
    row.names = NULL
  )

  targets_df <- data.frame(
    label = paste0("Target: ", result$target_estimates$label),
    estimate = result$target_estimates$transported_mean,
    lcl = result$target_estimates$lci,
    ucl = result$target_estimates$uci,
    lpi = result$target_estimates$lpi,
    upi = result$target_estimates$upi,
    type = "Target",
    row.names = NULL
  )

  plot_df <- rbind(studies_df, targets_df)
  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = estimate, y = label, color = type)) +
    ggplot2::geom_vline(xintercept = result$pooled_mean, linetype = "dashed", color = "grey45") +
    ggplot2::geom_segment(
      data = plot_df[plot_df$type == "Target", , drop = FALSE],
      ggplot2::aes(x = lpi, xend = upi, y = label, yend = label),
      inherit.aes = FALSE,
      linewidth = 1.3,
      alpha = 0.55,
      color = "grey55"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = lcl, xend = ucl, y = label, yend = label),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 2.7) +
    ggplot2::scale_color_manual(values = c(Study = "#4C6A92", Target = "#B34E36")) +
    ggplot2::labs(
      x = effect_label,
      y = NULL,
      color = NULL,
      title = "Transported Effect Estimates",
      subtitle = "Dashed vertical line is the baseline pooled mean; grey bars are target prediction intervals."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )
}

#' Save frequentist transport outputs
#'
#' Writes the main result tables and a forest plot to `output_dir`.
#'
#' @param result A fitted `"meta_transport_engine"` object.
#' @param output_dir Output directory path.
#'
#' @return Invisibly returns `output_dir`.
#' @export
save_transport_results <- function(result, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  utils::write.csv(result$target_estimates, file.path(output_dir, "target_estimates.csv"), row.names = FALSE)
  utils::write.csv(result$coefficients, file.path(output_dir, "coefficients.csv"), row.names = FALSE)
  utils::write.csv(result$overlap_checks, file.path(output_dir, "overlap_checks.csv"), row.names = FALSE)
  utils::write.csv(result$extrapolation_scores, file.path(output_dir, "extrapolation_scores.csv"), row.names = FALSE)
  utils::write.csv(result$transport_gaps, file.path(output_dir, "transport_gaps.csv"), row.names = FALSE)
  utils::write.csv(result$leave_one_out$detail, file.path(output_dir, "leave_one_out_detail.csv"), row.names = FALSE)
  utils::write.csv(result$leave_one_out$summary, file.path(output_dir, "leave_one_out_summary.csv"), row.names = FALSE)

  if (!is.null(result$bias_sensitivity)) {
    utils::write.csv(result$bias_sensitivity, file.path(output_dir, "bias_sensitivity.csv"), row.names = FALSE)
  }

  jsonlite::write_json(
    list(
      pooled_mean = result$pooled_mean,
      tau2_baseline = result$baseline_fit$tau2,
      tau2_transport = result$transport_fit$tau2,
      r2_meta = result$heterogeneity$r2_meta
    ),
    path = file.path(output_dir, "summary.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )

  plot_height <- max(6, 0.45 * (nrow(result$dat) + nrow(result$target_estimates)))
  ggplot2::ggsave(
    filename = file.path(output_dir, "transport_forest_plot.png"),
    plot = result$plot,
    width = 10,
    height = plot_height,
    dpi = 300
  )

  invisible(output_dir)
}

#' Fit a frequentist meta-transport engine
#'
#' Fits a baseline random-effects meta-analysis and a transport
#' meta-regression, then produces transported target estimates, overlap checks,
#' extrapolation diagnostics, leave-one-out instability summaries, and optional
#' risk-of-bias sensitivity results.
#'
#' @param dat Study-level data frame.
#' @param targets Target-population data frame.
#' @param moderators Character vector of moderator column names.
#' @param yi_col Name of the study effect estimate column.
#' @param vi_col Name of the within-study variance column.
#' @param study_col Name of the study label column.
#' @param target_label_col Name of the target label column.
#' @param effect_scale Effect scale used for optional back-transformation.
#' @param method Estimation method passed to [metafor::rma.uni()].
#' @param rob_col Optional risk-of-bias column used for sensitivity analysis.
#' @param high_risk_values Values in `rob_col` that should be excluded in the
#'   bias sensitivity run.
#'
#' @return An object of class `"meta_transport_engine"`.
#' @export
fit_transport_engine <- function(dat,
                                 targets,
                                 moderators,
                                 yi_col = "yi",
                                 vi_col = "vi",
                                 study_col = "study",
                                 target_label_col = "label",
                                 effect_scale = c("identity", "log_rr", "log_or", "log_hr"),
                                 method = "REML",
                                 rob_col = NULL,
                                 high_risk_values = c("high")) {
  effect_scale <- match.arg(effect_scale)

  assert_required_columns(dat, c(study_col, yi_col, vi_col, moderators), "dat")
  assert_required_columns(targets, moderators, "targets")
  assert_numeric_columns(dat, c(yi_col, vi_col, moderators), "dat")
  assert_numeric_columns(targets, moderators, "targets")

  if (!target_label_col %in% names(targets)) {
    targets[[target_label_col]] <- paste0("target_", seq_len(nrow(targets)))
  }

  warn_if_thin_data(dat, moderators)

  baseline_fit <- metafor::rma.uni(
    yi = dat[[yi_col]],
    vi = dat[[vi_col]],
    method = method,
    data = dat
  )

  transport_fit <- metafor::rma.uni(
    yi = dat[[yi_col]],
    vi = dat[[vi_col]],
    mods = build_modifier_formula(moderators),
    method = method,
    data = dat
  )

  pooled_mean <- as.numeric(stats::coef(baseline_fit)[1])
  coefficients <- summarise_coefficients(transport_fit)

  target_estimates <- do.call(
    rbind,
    lapply(seq_len(nrow(targets)), function(i) {
      predict_transport_target(
        fit = transport_fit,
        target_row = targets[i, , drop = FALSE],
        moderators = moderators,
        target_label = targets[[target_label_col]][i],
        effect_scale = effect_scale
      )
    })
  )
  target_estimates$delta_vs_pooled <- target_estimates$transported_mean - pooled_mean

  overlap_checks <- check_modifier_overlap(dat, targets, moderators, target_label_col)
  extrapolation_scores <- compute_extrapolation_scores(dat, targets, moderators, target_label_col)
  transport_gaps <- compute_transport_gaps(transport_fit, dat, target_estimates, moderators, study_col)
  loo <- leave_one_out_transport(dat, moderators, yi_col, vi_col, study_col, method)
  bias_sensitivity <- risk_of_bias_sensitivity(
    dat = dat,
    targets = targets,
    moderators = moderators,
    yi_col = yi_col,
    vi_col = vi_col,
    target_label_col = target_label_col,
    rob_col = rob_col,
    high_risk_values = high_risk_values,
    method = method,
    effect_scale = effect_scale
  )

  heterogeneity <- data.frame(
    tau2_baseline = baseline_fit$tau2,
    tau2_transport = transport_fit$tau2,
    r2_meta = if (baseline_fit$tau2 > 0) {
      1 - (transport_fit$tau2 / baseline_fit$tau2)
    } else {
      NA_real_
    },
    row.names = NULL
  )

  result <- list(
    dat = dat,
    targets = targets,
    spec = list(
      moderators = moderators,
      yi_col = yi_col,
      vi_col = vi_col,
      study_col = study_col,
      target_label_col = target_label_col,
      effect_scale = effect_scale,
      method = method
    ),
    baseline_fit = baseline_fit,
    transport_fit = transport_fit,
    pooled_mean = pooled_mean,
    coefficients = coefficients,
    target_estimates = target_estimates,
    overlap_checks = overlap_checks,
    extrapolation_scores = extrapolation_scores,
    transport_gaps = transport_gaps,
    leave_one_out = loo,
    bias_sensitivity = bias_sensitivity,
    heterogeneity = heterogeneity
  )

  result$plot <- plot_transport_forest(result, effect_label = toupper(effect_scale))
  class(result) <- "meta_transport_engine"
  result
}

#' Print a frequentist meta-transport fit
#'
#' @param x A fitted `"meta_transport_engine"` object.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.meta_transport_engine <- function(x, ...) {
  cat("Meta-Transport Engine\n")
  cat(sprintf("Studies: %d\n", nrow(x$dat)))
  cat(sprintf("Targets: %d\n", nrow(x$target_estimates)))
  cat(sprintf("Moderators: %s\n", paste(x$spec$moderators, collapse = ", ")))
  cat(sprintf("Pooled mean: %.4f\n", x$pooled_mean))
  cat(sprintf("Tau^2 baseline: %.4f\n", x$baseline_fit$tau2))
  cat(sprintf("Tau^2 transport: %.4f\n", x$transport_fit$tau2))
  invisible(x)
}
