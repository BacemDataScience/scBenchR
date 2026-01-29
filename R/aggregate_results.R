#' Aggregate benchmarking results from multiple runs
#'
#' This function aggregates standardized CSV benchmark outputs
#' (Step 7) into a single tidy table suitable for summarization
#' and plotting.
#'
#' @param path Character. Directory containing benchmark CSV files.
#' @return A tibble with aggregated benchmark results.
#' @export
aggregate_results <- function(path) {

  stopifnot(is.character(path), length(path) == 1)

  read_bench_csv(path)
}
