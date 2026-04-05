resolve_inst_file <- function(subdir, filename, package = "metaTransportEngine") {
  installed_path <- system.file(subdir, filename, package = package)
  if (nzchar(installed_path)) {
    return(installed_path)
  }

  source_path <- file.path("inst", subdir, filename)
  if (file.exists(source_path)) {
    return(normalizePath(source_path))
  }

  stop(
    sprintf("Could not locate bundled file '%s/%s'.", subdir, filename),
    call. = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_required_columns <- function(data, columns, data_name = "data") {
  missing_cols <- setdiff(columns, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        data_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

assert_numeric_columns <- function(data, columns, data_name = "data") {
  non_numeric <- columns[!vapply(data[columns], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      sprintf(
        "%s must contain numeric columns for: %s",
        data_name,
        paste(non_numeric, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

back_transform_values <- function(x, effect_scale = c("identity", "log_rr", "log_or", "log_hr")) {
  effect_scale <- match.arg(effect_scale)
  if (effect_scale == "identity") {
    x
  } else {
    exp(x)
  }
}

normalise_risk_values <- function(x) {
  tolower(trimws(as.character(x)))
}

summarise_draws_vector <- function(draws) {
  stats::setNames(
    c(
      mean = mean(draws),
      median = stats::median(draws),
      sd = stats::sd(draws),
      lci = stats::quantile(draws, 0.025),
      uci = stats::quantile(draws, 0.975)
    ),
    c("mean", "median", "sd", "lci", "uci")
  )
}

summarise_draws_matrix <- function(draws_matrix, labels, effect_scale = "identity") {
  if (is.null(dim(draws_matrix))) {
    draws_matrix <- matrix(draws_matrix, ncol = 1)
  }

  if (ncol(draws_matrix) != length(labels)) {
    stop("Labels do not match the number of draw columns.", call. = FALSE)
  }

  rows <- lapply(seq_len(ncol(draws_matrix)), function(i) {
    summary_i <- summarise_draws_vector(draws_matrix[, i])
    transformed <- back_transform_values(
      c(summary_i["mean"], summary_i["lci"], summary_i["uci"]),
      effect_scale = effect_scale
    )

    data.frame(
      label = labels[i],
      mean = unname(summary_i["mean"]),
      median = unname(summary_i["median"]),
      sd = unname(summary_i["sd"]),
      lci = unname(summary_i["lci"]),
      uci = unname(summary_i["uci"]),
      mean_bt = unname(transformed[1]),
      lci_bt = unname(transformed[2]),
      uci_bt = unname(transformed[3]),
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}
