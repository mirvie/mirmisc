#' @importFrom utils globalVariables
#' @importFrom magrittr '%>%' '%T>%'
#' @importFrom zeallot '%<-%'
#' @importFrom foreach '%dopar%'
NULL

## quiets concerns of R CMD check re: the .'s that appear when using pipes
if (getRversion() >= "2.15.1") {
  globalVariables(
    c(
      ".", ".x", "specificity", "sensitivity", "sens_lower", "sens_upper",
      "tbl", "spec", "sens", "model", "response", "predictor", "AUC", "x", "y",
      "gene", "meta_collectionga", "mirvie_id", "outlier", "tot_counts",
      "GeneName"
    )
  )
}
