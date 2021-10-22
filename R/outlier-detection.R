#' Detect the outlying samples in a cohort.
#'
#' This function wraps [mirmodels::compute_pcas()] and hence uses
#' [rrcov::PcaGrid()] to do robust PCA analysis and detect outliers.
#'
#' Prior to PCA calculation (and outlier detection), a call to
#' [mirmodels::linear_correct()] is made to regress away the effect of the
#' total number of counts on gene expression levels, with care taken to not
#' regress away the effect of gestational age.
#'
#' There's an Easter egg. You can pass a data frame directly as the `cohort`
#' argument and then the function will use that rather than having to call
#' `get_*_data()` to get the data. I advise `get_*_data(log2 = TRUE, tot_counts
#' = TRUE, gene_pred = ~median(.) > 0)`.
#'
#' @param cohort A two-character string, e.g. `"BW"`.
#'
#' @return An object of class `mirvie_cohort_outliers`. This is a data frame
#'   with 5 principal components named `PC1`, `PC2`, . . ., `PC5`. It also has
#'   columns `meta_collectionga`, `mirvie_id` and `outlier` which is a boolean
#'   column where `TRUE` indicates an outlier. This object has attributes
#'   `var_exp` and `loadings`. Read the documentation of
#'   [mirmodels::compute_pcas()] for more on those.
#'
#' @examples
#' if (require("mirmodels")) {
#'   ga_outliers <- get_cohort_outliers("ga")
#'   autoplot(ga_outliers)
#' }
#' @seealso [autoplot.mirvie_cohort_outliers()]
#'
#' @export
get_cohort_outliers <- function(cohort) {
  ensure_mirmodels()
  if (is.data.frame(cohort)) {
    data <- cohort
  } else {
    checkmate::assert_string(cohort)
    cohort <- stringr::str_trim(cohort)
    checkmate::assert_true(nchar(cohort) == 2)
    cohort <- tolower(cohort)
    call <- stringr::str_glue(
      "mirmodels::get_{cohort}_data(",
      "log2 = TRUE, tot_counts = TRUE, ",
      "gene_predicate = ~median(.) > 0)"
    )
    data <- eval(parse(text = call))
  }
  data <- dplyr::mutate(data, tot_counts = log2(tot_counts + 1))
  genes <- dplyr::intersect(names(data), get_gene_names())
  data <- data %>%
    dplyr::select(
      mirvie_id, meta_collectionga, tot_counts,
      dplyr::all_of(genes)
    )
  corrected_data <- mirmodels::linear_correct(
    data,
    correct_cols = genes,
    correct_for_cols = "tot_counts",
    keep_effect_cols = "meta_collectionga",
    robust = TRUE
  )[[1]]
  pca <- mirmodels::compute_pcas(corrected_data,
    subset = genes,
    normalize = TRUE, robust = TRUE
  )
  out <- pca %>%
    dplyr::mutate(outlier = attr(., "outlier")) %>%
    structure(
      var_exp = attr(pca, "var_exp"),
      loadings = attr(pca, "loadings"),
      class = c("mirvie_cohort_outliers", class(.))
    )
  out
}

#' Plot a `mirvie_cohort_outliers` object.
#'
#' A `mirvie_cohort_outliers` object is the output of a call to
#' [get_cohort_outliers()].
#'
#' @param object A `mirvie_cohort_outliers` object.
#' @param pcx An integer between 1 and 5. The principal component that will be
#'   on the x axis.
#' @param pcy An integer between 1 and 5. The principal component that will be
#'   on the y axis.
#' @param plotly A flag. Make the plot interactive (with mirvie ID tooltips)?
#' @param ... Not currently used.
#'
#' @return A [ggplot2::ggplot()] or a [plotly::ggplotly()].
#'
#' @export
autoplot.mirvie_cohort_outliers <- function(object, pcx = 1, pcy = 2,
                                            plotly = interactive(), ...) {
  checkmate::assert_class(object, "mirvie_cohort_outliers")
  pcx <- checkmate::assert_int(pcx, lower = 1, upper = 5, coerce = TRUE)
  pcy <- checkmate::assert_int(pcy, lower = 1, upper = 5, coerce = TRUE)
  pcx <- paste0("PC", pcx)
  pcy <- paste0("PC", pcy)
  object$pcx <- object[[pcx]]
  object$pcy <- object[[pcy]]
  cohort <- object$mirvie_id %>%
    stringr::str_sub(1, 2) %>%
    unique() %>%
    paste(collapse = ",")
  out <- ggplot2::ggplot(
    object,
    ggplot2::aes(pcx, pcy, color = outlier, label = mirvie_id)
  ) +
    ggplot2::geom_point() +
    ggthemes::scale_color_colorblind() +
    ggplot2::xlab(pcx) +
    ggplot2::ylab(pcy) +
    ggplot2::ggtitle(paste(toupper(cohort), "outliers"))
  if (plotly) out <- plotly::ggplotly(out, tooltip = "label")
  out
}
