#' Time an expression (seconds)
#' @keywords internal
scbenchr_time_seconds <- function(expr) {
  t0 <- Sys.time()
  out <- eval.parent(substitute(expr))
  t1 <- Sys.time()
  list(value = out, sec = as.numeric(difftime(t1, t0, units = "secs")))
}

#' Time and return a standardized timing row
#' @keywords internal
scbenchr_time_it <- function(step, expr) {
  t0 <- Sys.time()
  out <- eval.parent(substitute(expr))
  t1 <- Sys.time()
  list(
    value = out,
    timing = data.frame(
      step = step,
      elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")),
      stringsAsFactors = FALSE
    )
  )
}
