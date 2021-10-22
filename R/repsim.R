#' Repeatedly simulate gene counts for cases and controls of a condition.
#'
#' Given a number of gene repetitions, a number of samples, the condition
#' prevalence and the theoretical mean of a gene for cases and controls,
#' repeatedly simulate gene counts for that gene across the condition.
#'
#' Poisson gene counts are assumed.
#'
#' @param n_gene_reps The number of times to repeatedly simulate the gene counts
#'   for that gene.
#' @param n_samples The number of samples in the simulation.
#' @param cond_prevalence A number in (0, 1). The condition prevalence.
#' @param control_mean Theoretical mean of the gene for controls.
#' @param case_mean Theoretical mean of the gene for cases.
#'
#' @return A tibble with `n_gene_reps + 1` columns called `generep_1`,
#'   `generep_2`, . . ., `generep_{n_gene_reps}`, `cond` and `n_samples`
#'   columns.
#'
#' @examples
#' if (rlang::is_installed("mirmodels")) {
#'   rsgc5000 <- repsim_gene_cond(15000, 5000, 1 / 10, 0.01, 0.05)
#'   sde <- mirmodels::cor_de(rsgc5000, "cond", head(names(rsgc5000), -1))
#' }
#' @export
repsim_gene_cond <- function(n_gene_reps, n_samples, cond_prevalence,
                             control_mean, case_mean) {
  checkmate::assert_count(n_gene_reps)
  checkmate::assert_count(n_samples)
  checkmate::assert_number(cond_prevalence, lower = 0, upper = 1)
  checkmate::assert_number(control_mean, lower = 0)
  checkmate::assert_number(case_mean, lower = 0)
  cond_col <- purrr::rbernoulli(n_samples, cond_prevalence)
  if (dplyr::n_distinct(cond_col) == 1) {
    custom_stop(
      "In this simulation, only
       {dplyr::if_else(cond_col[1], 'cases', 'controls')} were simulated.",
      "You need to do the simulation in such a way that both cases and controls
       are present."
    )
  }
  control_df <- stats::rpois(n_gene_reps * sum(!cond_col), control_mean) %>%
    matrix(ncol = n_gene_reps) %>%
    magrittr::set_colnames(paste0("generep_", seq_len(n_gene_reps))) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(cond = FALSE)
  case_df <- stats::rpois(n_gene_reps * sum(cond_col), case_mean) %>%
    matrix(ncol = n_gene_reps) %>%
    magrittr::set_colnames(paste0("generep_", seq_len(n_gene_reps))) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(cond = TRUE)
  dplyr::bind_rows(control_df, case_df) %>%
    dplyr::sample_frac()
}
