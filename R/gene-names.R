#' Get the names of all of the genes that we use.
#'
#' This includes the Ensembl and colloquial names (so in that sense there is
#' duplication).
#'
#' @return A character vector.
#'
#' @examples
#' get_gene_names()
#' @export
get_gene_names <- function() {
  ensembl_gene_names <- system.file("extdata", "ensembl-gene-names.txt",
    package = "mirmisc"
  ) %>%
    readr::read_lines(progress = FALSE) %>%
    stringr::str_subset(stringr::coll("ENSG"))
  colloquial_gene_names <- system.file("extdata", "ensembl-to-gene-curated.csv",
    package = "mirmisc"
  ) %>%
    readr::read_csv(col_types = readr::cols(), progress = FALSE) %>%
    dplyr::pull(GeneName)
  unique(c(colloquial_gene_names, ensembl_gene_names))
}

#' Fill missing gene columns in one data frame with those from another.
#'
#' The 'gene names' are those returned by [get_gene_names()]. This function
#' takes two data frames `df` and `df_fill_from` and, if there are any columns
#' in `df_fill_from` with gene names as column names which don't exist in `df`,
#' they are copied into `df`.
#'
#' If the gene names in `df` are contiguously located, the copied genes are
#' inserted right after those. Otherwise, they are inserted on the end.
#'
#' @param df A data frame.
#' @param df_fill_from A data frame with the same number of rows as `df`.
#'
#' @return A data frame.
#'
#' @examples
#' if (require("mirmodels")) {
#'   st_data <- get_st_data()
#'   st_data_median0 <- get_st_data(gene_predicate = ~ median(.) == 0)
#'   dim(st_data)
#'   dim(st_data_median0)
#'   dim(df_fill_missing_genes(st_data_median0, st_data))
#' }
#' @export
df_fill_missing_genes <- function(df, df_fill_from) {
  checkmate::assert_data_frame(df, col.names = "unique")
  checkmate::assert_data_frame(df_fill_from, col.names = "unique")
  checkmate::assert_true(nrow(df) == nrow(df_fill_from))
  gene_names <- get_gene_names()
  df_gene_names <- dplyr::intersect(names(df), gene_names)
  df_gene_locs <- match(df_gene_names, names(df))
  contig <- all(diff(df_gene_locs) == 1)
  df_ff_gene_names <- dplyr::intersect(names(df_fill_from), gene_names)
  df_ff_genes_to_fill <- setdiff(df_ff_gene_names, df_gene_names)
  if (length(df_ff_genes_to_fill)) {
    if (contig) {
      if (max(df_gene_locs) == ncol(df)) {
        df <- dplyr::bind_cols(df, df_fill_from[df_ff_genes_to_fill])
      } else {
        df <- dplyr::bind_cols(
          df[seq_len(max(df_gene_locs))],
          df_fill_from[df_ff_genes_to_fill],
          df[seq(max(df_gene_locs) + 1, ncol(df))]
        )
      }
    } else {
      df <- dplyr::bind_cols(df, df_fill_from[df_ff_genes_to_fill])
    }
  }
  df
}

#' Which gene names are also column names?
#'
#' This is just `intersect(names(df), get_gene_names())`.
#'
#' @param df A data frame.
#'
#' @return A character vector.
#'
#' @export
get_df_gene_names <- function(df) {
  checkmate::assert_data_frame(df, col.names = "named")
  dplyr::intersect(names(df), get_gene_names())
}

#' Apply the same function to all columns whose names are gene names.
#'
#' Gene names are elements of [get_gene_names()].
#'
#' @param df A data frame.
#' @param f A function or a formula coercible to a function by
#'   [rlang::as_function()].
#'
#' @return A data frame.
#'
#' @export
mutate_genes <- function(df, f) {
  checkmate::assert_data_frame(df, col.names = "named")
  if (rlang::is_formula(f)) f <- rlang::as_function(f)
  df_gene_names <- get_df_gene_names(df)
  for (gn in df_gene_names) {
    df[[gn]] <- f(df[[gn]])
  }
  df
}

#' Filter genes using a predicate function.
#'
#' Gene names are elements of [get_gene_names()]. Genes for whom the predicate
#' function `f()` evaluates to `FALSE` are dropped.
#'
#' @param df A data frame some of whose column names are gene names.
#' @param f A predicate function or a formula coercible to a function by
#'   [rlang::as_function()].
#'
#' @return A data frame.
#'
#' @export
filter_genes <- function(df, f) {
  checkmate::assert_data_frame(df, col.names = "named")
  if (rlang::is_formula(f)) f <- rlang::as_function(f)
  df_gene_names <- get_df_gene_names(df)
  for (gn in df_gene_names) {
    if (!f(df[[gn]])) df[[gn]] <- NULL
  }
  df
}

#' Convert gene names to/from Ensembl.
#'
#' This function works in a particular way: inputs that don't look like a gene
#' at all are returned as is.
#'
#' @param x A character vector.
#' @param ensembl A string. Either `"from"` or `"to"`. With `"to"` the result is
#'   Ensembl gene names. With `"from"` the result is colloquial gene names.
#'
#' @return A character vector.
#'
#' @export
convert_gene_names <- function(x, ensembl = "from") {
  checkmate::assert_character(x, min.chars = 1, min.len = 1)
  ensembl <- strex::match_arg(ensembl, c("from", "to"), ignore_case = TRUE)
  ensembl_colloquial <- readr::read_csv(
    system.file("extdata", "ensembl-to-gene-curated.csv",
                package = "mirmisc"),
    progress = FALSE, col_types = readr::cols()
  )
  if (ensembl == "from") {
    conversion_vector <- rlang::set_names(ensembl_colloquial$GeneName,
                                          ensembl_colloquial$EnsemblGeneID)
  } else {
    conversion_vector <- rlang::set_names(ensembl_colloquial$EnsemblGeneID,
                                          ensembl_colloquial$GeneName)
  }
  out <- conversion_vector[x]
  dplyr::if_else(is.na(out), x, out)
}
