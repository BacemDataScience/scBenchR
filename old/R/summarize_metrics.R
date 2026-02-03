#' Summarize benchmark metrics
#'
#' Computes summary statistics (n, mean, sd, median, min, max) for elapsed time,
#' grouped by dataset/method/step. If optional columns such as `peak_ram_mb` are
#' present, they are summarized as well.
#'
#' @param x A data frame/tibble of benchmark results (e.g., from `aggregate_results()`).
#' @param group_vars Character vector of grouping columns.
#' @return A tibble of summary statistics.
#' @export
summarize_metrics <- function(
    x,
    group_vars = c("dataset_id", "method", "step")
) {
  stopifnot(is.data.frame(x))
  stopifnot(is.character(group_vars))

  missing_g <- setdiff(group_vars, names(x))
  if (length(missing_g) > 0) {
    stop("Missing grouping columns: ", paste(missing_g, collapse = ", "), call. = FALSE)
  }
  if (!("elapsed_sec" %in% names(x))) {
    stop("Missing required column: elapsed_sec", call. = FALSE)
  }

  g <- dplyr::syms(group_vars)

  out <- x |>
    dplyr::group_by(!!!g) |>
    dplyr::summarise(
      n = dplyr::n(),
      elapsed_mean = mean(.data$elapsed_sec, na.rm = TRUE),
      elapsed_sd   = stats::sd(.data$elapsed_sec, na.rm = TRUE),
      elapsed_med  = stats::median(.data$elapsed_sec, na.rm = TRUE),
      elapsed_min  = min(.data$elapsed_sec, na.rm = TRUE),
      elapsed_max  = max(.data$elapsed_sec, na.rm = TRUE),
      .groups = "drop"
    )

  # Optional memory column if present
  if ("peak_ram_mb" %in% names(x)) {
    out_mem <- x |>
      dplyr::group_by(!!!g) |>
      dplyr::summarise(
        ram_mean = mean(.data$peak_ram_mb, na.rm = TRUE),
        ram_sd   = stats::sd(.data$peak_ram_mb, na.rm = TRUE),
        ram_med  = stats::median(.data$peak_ram_mb, na.rm = TRUE),
        ram_min  = min(.data$peak_ram_mb, na.rm = TRUE),
        ram_max  = max(.data$peak_ram_mb, na.rm = TRUE),
        .groups = "drop"
      )
    out <- dplyr::left_join(out, out_mem, by = group_vars)
  }

  tibble::as_tibble(out)
}
