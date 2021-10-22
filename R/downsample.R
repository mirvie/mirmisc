#' Downsample a vector of counts.
#'
#' Uniformly downsample a vector of non-negative integers to have a specified
#' sum. That is, keep randomly subtracting 1 from nonzero elements of the vector
#' until it has the desired sum. Each count is equally likely to be taken. That
#' is, an element with value 8 is 4 times more likely to be decremented than an
#' element with value 2.
#'
#' `downsample_count_mat_rows()` and `downsample_count_mat_cols()` just do
#' `downsample_vec()` to all rows and columns of a matrix using [apply()].
#'
#' `downsample_gene_counts()` downsamples on the subset of columns in the
#' data frame that have names in `get_gene_names()`.
#'
#' If `end_sum > sum(vec)`, `vec` is returned unchanged.
#'
#' @param vec A vector of non-negative integers.
#' @param end_sum The number that you would like vector to sum to after
#'   downsampling. This must be less than the initial sum.
#'
#' @return A vector of non-negative integers.
#'
#' @examples
#' downsample_count_vec(1:24, 24)
#' @export
downsample_count_vec <- function(vec, end_sum) {
  checkmate::assert_integerish(vec, lower = 0, min.len = 1, any.missing = FALSE)
  checkmate::assert_count(end_sum)
  vec_sum <- sum(vec)
  n_to_draw <- vec_sum - end_sum
  draws <- 0
  if (n_to_draw > 0) {
    draws <- detrendr::rfromboxes(n_to_draw, vec, vec)
  }
  vec - draws
}

#' @rdname downsample_count_vec
#' @param mat A matrix of non-negative integers.
#' @examples
#' mat <- matrix(sample.int(100, size = 6^2, replace = TRUE), nrow = 6)
#' downsample_count_mat_rows(mat, end_sum = 6)
#' @export
downsample_count_mat_rows <- function(mat, end_sum) {
  t(apply(mat, 1, downsample_count_vec, end_sum = end_sum))
}

#' @rdname downsample_count_vec
#' @param mat A matrix of non-negative integers.
#' @examples
#' downsample_count_mat_cols(mat, end_sum = 6)
#' @export
downsample_count_mat_cols <- function(mat, end_sum) {
  apply(mat, 2, downsample_count_vec, end_sum = end_sum)
}

#' @rdname downsample_count_vec
#' @param df A data frame with gene names as columns.
#' @examples
#' if (rlang::is_installed("mirmodels")) {
#'   ms_data <- mirmodels::get_ms_data(gene_predicate = ~ median(.) > 0)
#'   downsampled_ms <- downsample_gene_counts(ms_data, end_sum = 1e6)
#' }
#' @export
downsample_gene_counts <- function(df, end_sum) {
  checkmate::assert_data_frame(df, col.names = "named")
  df_names <- names(df)
  checkmate::assert_true(any(get_gene_names() %in% names(df)))
  gene_mat <- as.matrix(dplyr::select(df, dplyr::any_of(get_gene_names())))
  downsampled_gene_mat <- downsample_count_mat_rows(gene_mat, end_sum)
  df_without_genes <- dplyr::select(df, -dplyr::any_of(get_gene_names()))
  dplyr::bind_cols(df_without_genes, dplyr::as_tibble(downsampled_gene_mat)) %>%
    dplyr::select(dplyr::all_of(df_names))
}
