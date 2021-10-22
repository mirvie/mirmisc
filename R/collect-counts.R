#' Collect the count files of several samples into a single data frame.
#'
#' This function takes a directory path `dir_path` and searches that directory
#' for files whose names end in '_counts.txt'. It reads those files and
#' concatenates them. Each file is assumed to correspond to a single sample
#' whose name is contained in the first part of the file name (the bit before
#' '_counts.txt'). These sample names are used as column names in the output
#' data frame.
#'
#' @param dir_path A character vector. The path to the directory containing
#'   '*_counts.txt' files. To specify several directories, use a list of paths.
#' @param convert_genenames A flag. Convert gene names from Ensembl IDs to more
#'   widely-used names?
#' @param cpm A flag. Convert raw counts to counts-per-million (on a per sample
#'   basis)?
#' @param log2 A flag. Transform the counts using `log2(x + 1)`?
#' @param remove_ercc A flag. Remove ERCC counts from the results?
#' @param remove_rp_4_11 A flag. RP4 and RP11 genes come from a particular donor
#'   when building the genome and are mostly not useful. The default (`TRUE`) is
#'   to remove them.
#' @param remove_metadata A flag. The count files can contain non-gene metadata
#'   (e.g. mapping stats). The default (`TRUE`) is to remove them.
#' @param remove_controls A flag. The count files can contain positive and
#'   negative controls (denoted by having names ending in 'PC_counts.txt' or
#'   'NC_counts.txt'). The default (`TRUE`) is to remove them.
#' @param write A flag or string. Write the results to disk as a tab-separated
#'   file? If `TRUE`, the file will be written to the working directory with
#'   name 'genes_counts.txt', 'genes_cpm.txt', genes_log2.txt' or
#'   'genes_log2_cpm.txt'. To write the file elsewhere, pass the path through
#'   this argument as a string.
#'
#' @return A data frame object.
#'
#' @examples
#' \dontrun{
#' collect_counts("path/to/dir/with/count/files")
#' }
#'
#' @export
collect_counts <- function(dir_path,
                           convert_genenames = TRUE,
                           cpm = FALSE,
                           log2 = FALSE,
                           remove_ercc = TRUE,
                           remove_rp_4_11 = TRUE,
                           remove_metadata = TRUE,
                           remove_controls = TRUE,
                           write = FALSE) {
  checkmate::assert_character(dir_path,
    unique = TRUE,
    min.len = 1, min.chars = 1
  )
  for (dp in dir_path) checkmate::assert_directory_exists(dp)
  checkmate::assert_flag(convert_genenames)
  checkmate::assert_flag(cpm)
  checkmate::assert_flag(log2)
  checkmate::assert_flag(remove_ercc)
  checkmate::assert_flag(remove_rp_4_11)
  checkmate::assert_flag(remove_metadata)
  checkmate::assert_flag(remove_controls)
  checkmate::assert(
    checkmate::check_flag(write),
    checkmate::check_string(write, min.chars = 1)
  )
  file_paths <- fs::dir_ls(dir_path, type = "file", glob = "*_counts.txt")
  if (remove_controls) {
    file_paths <- stringr::str_subset(file_paths, "[PN]C_counts\\.txt$",
      negate = TRUE
    )
  }
  sample_names <- file_paths %>%
    fs::path_file() %>%
    strex::str_before_first(stringr::coll("_counts.txt"))
  counts <- purrr::map2(
    file_paths, sample_names,
    ~ suppressMessages(
      readr::read_tsv(.x, col_names = c("gene", .y), progress = FALSE)
    )
  ) %>%
    purrr::reduce(dplyr::full_join, by = c("gene" = "gene"))
  counts <- clean_ccs(counts,
    convert_genenames = convert_genenames,
    remove_ercc = remove_ercc,
    remove_rp_4_11 = remove_rp_4_11,
    remove_metadata = remove_metadata
  )
  if ((!cpm) && (!log2)) {
    out <- counts
    write_ccs_if_appropriate(out, write, "genes_counts.txt")
  } else if (cpm && (!log2)) {
    out <- dplyr::mutate_if(counts, is.numeric, ~ . / sum(.) * 1e6)
    write_ccs_if_appropriate(out, write, "genes_cpm.txt")
  } else if ((!cpm) && log2) {
    out <- dplyr::mutate_if(counts, is.numeric, ~ log2(. + 1))
    write_ccs_if_appropriate(out, write, "genes_log2.txt")
  } else { # cpm and log2
    out <- dplyr::mutate_if(
      counts, is.numeric,
      ~ log2((. / sum(.) * 1e6) + 1)
    )
    write_ccs_if_appropriate(out, write, "genes_log2_cpm.txt")
  }
  out
}

#' Clean the collected counts.
#'
#' This is a helper function used by [collect_counts()]. It makes sure the
#' count table is integer type (if appropriate) and optionally converts gene
#' names and removes unwanted "genes" such as ERCCs.
#'
#' @param counts Counts data frame collected from several `*_counts.txt` files
#'   that has not yet been cleaned.
#' @inheritParams collect_counts
#'
#' @return A data frame.
#'
#' @noRd
clean_ccs <- function(counts, convert_genenames,
                      remove_ercc, remove_rp_4_11,
                      remove_metadata) {
  checkmate::assert_data_frame(counts, min.rows = 1, min.cols = 2)
  checkmate::assert_subset("gene", names(counts))
  checkmate::assert_flag(convert_genenames)
  checkmate::assert_flag(remove_ercc)
  checkmate::assert_flag(remove_rp_4_11)
  checkmate::assert_flag(remove_metadata)
  if (convert_genenames) {
    curated_gene_conversion_df <- suppressMessages(readr::read_csv(
      system.file("extdata", "ensembl-to-gene-curated.csv",
        package = "mirmisc"
      )
    ))
    joined <- dplyr::left_join(counts, curated_gene_conversion_df,
      by = c("gene" = "EnsemblGeneID")
    )
    counts$gene <- dplyr::if_else(
      is.na(joined$GeneName),
      counts$gene, joined$GeneName
    )
  }
  if (remove_ercc) {
    counts <- dplyr::filter(
      counts,
      stringr::str_detect(gene, "ERCC-", negate = TRUE)
    )
  }
  if (remove_rp_4_11) {
    counts <- dplyr::filter(
      counts,
      stringr::str_detect(gene, "RP4-|RP11-", negate = TRUE)
    )
  }
  if (remove_metadata) {
    counts <- dplyr::filter(
      counts,
      stringr::str_detect(gene, "__|NIST", negate = TRUE)
    )
  }
  for (nm in names(counts)) {
    cts_nm <- counts[[nm]]
    if (is.numeric(cts_nm) && isTRUE(all.equal(cts_nm, floor(cts_nm)))) {
      counts[[nm]] <- as.integer(cts_nm)
    }
  }
  counts
}

#' Write the collected counts to disk as a tsv file.
#'
#' Writing is only done if `write` is a string or `TRUE`.
#'
#' @param df The data frame to write to disk.
#' @param write A string or a flag. If `True`, the file is written with name
#'   equal to the `default` argument. If a string, the file is written with that
#'   name.
#' @param default A string.
#'
#' @return Invisibly, `TRUE` if writing took place and `FALSE` otherwise.
#'
#' @noRd
write_ccs_if_appropriate <- function(df, write, default) {
  checkmate::assert_data_frame(df, min.rows = 1, min.cols = 2)
  checkmate::assert_subset("gene", names(df))
  checkmate::assert(
    checkmate::check_flag(write),
    checkmate::check_string(write, min.chars = 1)
  )
  checkmate::assert_string(default, min.chars = 1)
  if (is.character(write) || write) {
    out_name <- ifelse(is.character(write), write, default)
    suppressMessages(readr::write_tsv(df, out_name))
    invisible(TRUE)
  }
  invisible(FALSE)
}

#' Get the paths to `htseq/` directories for a given cohort.
#'
#' The `*_counts.txt` files live in directories called `htseq/`. This function
#' helps you to find all such directories for a given cohort.
#'
#' @param base_dir A string. The path to a directory that the cohort directories
#'   live under. The cohort directories have name structure `###_XY` where `#`
#'   is a digit and `XY` is the cohort code. For example, `007_RS`.
#' @param cohort_code A string with exactly two characters. E.g. `"RS"`.
#'
#' @return A character vector of paths.
#'
#' @examples
#' \dontrun{
#' get_htseq_paths(cohort_code = "RS")
#' }
#'
#' @export
get_htseq_paths <- function(base_dir = "/mnt/storage/Cohorts", cohort_code) {
  checkmate::assert_string(base_dir)
  checkmate::assert_directory_exists(base_dir)
  checkmate::assert_string(cohort_code)
  checkmate::assert_true(nchar(cohort_code) == 2)
  rgx <- paste0("\\d\\d\\d_", cohort_code, "$")
  cohort_dir <- fs::dir_ls(base_dir, recurse = TRUE, regexp = rgx)
  if (length(cohort_dir) == 0) {
    custom_stop("Cohort directory for cohort code {cohort_code} not found.")
  }
  if (length(cohort_dir) > 1) {
    n_slash <- stringr::str_count(cohort_dir, stringr::coll("/"))
    cohort_dir <- cohort_dir[n_slash == min(n_slash)][[1]]
  }
  fs::dir_ls(cohort_dir, recurse = TRUE, type = "dir", glob = "*/htseq")
}
