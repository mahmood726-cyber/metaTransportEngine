#' Locate the bundled example studies CSV
#'
#' Returns the path to the packaged example study-level dataset.
#'
#' @return A single file path.
#' @export
example_studies_file <- function() {
  resolve_inst_file("extdata", "studies_example.csv")
}

#' Locate the bundled example targets CSV
#'
#' Returns the path to the packaged example target-population dataset.
#'
#' @return A single file path.
#' @export
example_targets_file <- function() {
  resolve_inst_file("extdata", "targets_example.csv")
}

#' Locate the bundled TruthCert-style schema example
#'
#' Returns the path to the packaged YAML schema example.
#'
#' @return A single file path.
#' @export
example_truthcert_schema_file <- function() {
  resolve_inst_file("extdata", "truthcert_meta_transport_example.yaml")
}

#' Read the bundled example studies dataset
#'
#' @return A data frame of study-level effect estimates and modifiers.
#' @export
read_example_studies <- function() {
  utils::read.csv(example_studies_file(), stringsAsFactors = FALSE)
}

#' Read the bundled example targets dataset
#'
#' @return A data frame of target modifier profiles.
#' @export
read_example_targets <- function() {
  utils::read.csv(example_targets_file(), stringsAsFactors = FALSE)
}
