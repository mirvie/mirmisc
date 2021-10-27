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

#' Get the longest common substring within two strings.
#'
#' @param x,y Strings.
#'
#' @return A string.
#'
#' @noRd
longest_common_substring <- function(x, y) {
  checkmate::assert_string(x)
  checkmate::assert_string(y)
  if (min(nchar(x), nchar(y)) == 0) return(character(0L))
  x <- strex::str_to_vec(x)
  y <- strex::str_to_vec(y)
  paste(qualV::LCS(x, y)$LCS, collapse = "")
}
