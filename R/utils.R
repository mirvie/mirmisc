#' Check that `mirmodels` is installed and error if not.
#'
#' @return `TRUE`, invisibly.
#'
#' @noRd
ensure_mirmodels <- function() {
  if (!requireNamespace("mirmodels", quietly = TRUE)) {
    stop("Ensure that mirmodels is installed.")
  }
  invisible(TRUE)
}
