#' Default base directory for scBenchR runs (Windows)
#' @keywords internal
scbenchr_default_base_dir <- function() "D:/scRNAseq"

#' Create output directory structure
#'
#' @param out_dir Output directory.
#' @return The normalized output directory path (invisibly).
#' @export
scbenchr_init_outdirs <- function(out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "markers"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "bench"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "integrated"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "paper"), showWarnings = FALSE, recursive = TRUE)
  invisible(normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}

#' Write session info for reproducibility
#'
#' @param out_dir Output directory.
#' @param filename Name of the session info file.
#' @export
scbenchr_write_session_info <- function(out_dir, filename = "session_info.txt") {
  scbenchr_init_outdirs(out_dir)
  f <- file.path(out_dir, filename)
  sink(f)
  cat("==== Session Info ====\n")
  print(Sys.info())
  cat("\n==== R Session ====\n")
  print(utils::sessionInfo())
  sink()
  invisible(f)
}
