roxy_de_return <- function() {
  c(
    "A tibble with 4 columns:",
    "* `gene`: The gene name.",
    "* `log2fc`: The log2 fold-change between case and control.",
    "* `pvalue`: The p-value for that gene being differentially expressed.",
    "* `padj`: Adjusted p-values."
  )
}
