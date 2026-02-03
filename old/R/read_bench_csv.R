#' Read benchmark CSV files from a directory
#'
#' Reads all `.csv` files in `path`, binds them into a single table,
#' and performs minimal validation (non-empty + required columns).
#'
#' @param path Character. Directory containing benchmark CSV files.
#' @return A tibble with concatenated benchmark results.
#' @export
read_bench_csv <- function(path) {

  stopifnot(is.character(path), length(path) == 1)

  if (!dir.exists(path)) {
    stop("Directory does not exist: ", path, call. = FALSE)
  }

  # Avoid regex issues on Windows: list all files, then filter by extension.
  files <- list.files(path, full.names = TRUE)
  files <- files[tolower(tools::file_ext(files)) == "csv"]

  if (length(files) == 0) {
    stop("No .csv files found in: ", path, call. = FALSE)
  }

  out <- dplyr::bind_rows(lapply(files, function(f) {
    x <- readr::read_csv(f, show_col_types = FALSE)
    x$.source_file <- basename(f)  # provenance
    x
  }))

  if (nrow(out) == 0) {
    stop("CSV files were read but produced 0 rows: ", path, call. = FALSE)
  }

  required <- c("dataset_id", "method", "run_id", "step", "elapsed_sec")
  missing <- setdiff(required, names(out))
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  tibble::as_tibble(out)
}
