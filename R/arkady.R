#' Check the sample descriptions data frame for an R&D experiment.
#'
#' * It must contain a column called `sample_name` with unique values.
#' * It must contain a column called `condition`.
#' * If it contains a column called `id_tech_rep`, then for a given
#'   `id_tech_rep`, the `sample_name`s must have a common substring at least two
#'   characters long.
#'
#' @param sample_descs A data frame of sample descriptions that includes
#'   columns called `sample_name` and `condition`.
#'
#' @return `TRUE` (invisibly) if the check passes. Otherwise, an error is
#'   thrown.
#'
#' @family Arkady
#' @export
check_sample_descs <- function(sample_descs) {
  checkmate::assert_data_frame(sample_descs)
  sample_descs <- janitor::clean_names(sample_descs)
  checkmate::assert_names(names(sample_descs),
                          must.include = c("sample_name", "condition"))
  if ("id_tech_rep" %in% names(sample_descs)) {
    collapsed_id_df <- sample_descs %>%
      dplyr::group_by(.data[["id_tech_rep"]]) %>%
      dplyr::summarise(
        collapsed_id = purrr::reduce(.data[["sample_name"]],
                                     longest_common_substring)
      )
    first_dup <- anyDuplicated(collapsed_id_df$collapsed_id)
    if (first_dup) {
      custom_stop(
        "The collapsed sample names are not unique.",
        "id_tech_rep '{collapsed_id_df[[first_dup, 'id_tech_rep']]}' is the
         offender."
      )
    }
  }
  invisible(TRUE)
}

#' Match sample description sheet names to count file names.
#'
#' Given a data frame of sample descriptions from the lab team and a data frame
#' from a `genes_counts.txt` file from a sequencing run, match the sample names
#' and return a data frame of sample descriptions detailing the matching.
#'
#' Name matching is done by first replacing all whitespace and underscores with
#' hyphens and then asserting that the `sampleName` must a substring of the
#' column name. The matching is done via a pair of columns `id_sample_descs`
#' and `id_genes_counts` in the output.
#'
#' @param sample_descs A data frame of sample descriptions. The `sampleName`
#'   column holds the sample names.
#' @param genes_counts A data frame of gene counts. The first column is `gene`
#'   and the rest are sample names.
#'
#' @return A data frame.
#'
#' @family Arkady
#' @export
reconcile_names <- function(sample_descs, genes_counts) {
  checkmate::assert_data_frame(sample_descs)
  sample_descs <- janitor::clean_names(sample_descs)
  check_sample_descs(sample_descs)
  checkmate::assert_data_frame(genes_counts)
  checkmate::assert_names(names(genes_counts), must.include = "gene")
  genes_counts_cleaned_names <- stringr::str_replace_all(names(genes_counts),
                                                         "[-._ ]+", "_")
  sample_descs %>%
    dplyr::rename(id_sample_descs = "sample_name") %>%
    dplyr::mutate(
      id_sample_descs_clean = stringr::str_replace_all(
        .data[["id_sample_descs"]],
        "[-._ ]+", "_"
      ),
      id_genes_counts_indices = purrr::map(
        .data[["id_sample_descs_clean"]],
        ~stringr::str_which(genes_counts_cleaned_names,
                            stringr::coll(.))
      )
    ) %>%
    dplyr::filter(lengths(id_genes_counts_indices) == 1) %>%
    dplyr::mutate(
      id_genes_counts_indices = unlist(.data[["id_genes_counts_indices"]]),
      id_genes_counts = names(genes_counts)[
        .data[["id_genes_counts_indices"]]
      ]
    ) %>%
    dplyr::select(id_sample_descs, "id_genes_counts", dplyr::everything(),
                  -"id_sample_descs_clean", -"id_genes_counts_indices")
}

#' Join the sample descriptions and genes counts data frames into a wide format.
#'
#' This function calls [reconcile_names()] and [check_sample_descs()]
#' internally.
#'
#' @inheritParams reconcile_names
#'
#' @return A data frame. The gene counts table with columns corresponding to the
#'   desired subset of data.
#'
#' @family Arkady
#' @export
prep_rd_input <- function(sample_descs, genes_counts) {
  sample_descs <- reconcile_names(sample_descs = sample_descs,
                                  genes_counts = genes_counts)
  genes_counts <- data.frame(genes_counts, row.names = "gene") %>%
    t() %>%
    dplyr::as_tibble(rownames = "id_genes_counts")
  dplyr::inner_join(sample_descs, genes_counts, by = "id_genes_counts")
}

#' Collapse (sum) technical replicates in a set of R&D experiments.
#'
#' Technical replicates, identified by the `id_tech_rep` column, are summed into
#' a single sample. The resulting `id_genes_counts` is the longest common
#' substring of the `id_genes_counts` of the input technical replicates.
#'
#' @inheritParams plot_var_mean_ratio
#'
#' @return A dataframe with columns `id_genes_counts`, `condition` and the gene
#'   columns.
#' @family Arkady
#' @export
collapse_tech_reps <- function(prepped_rd_input) {
  checkmate::assert_data_frame(prepped_rd_input)
  checkmate::assert_names(
    names(prepped_rd_input),
    must.include = c("id_genes_counts", "id_tech_rep", "condition")
  )
  prepped_rd_input <- dplyr::select(
    prepped_rd_input,
    "id_genes_counts", "id_tech_rep", "condition",
    dplyr::any_of(get_gene_names())
  )
  prepped_rd_input %>%
    dplyr::group_by(.data[["id_tech_rep"]]) %>%
    dplyr::summarise(
      dplyr::bind_cols(
        id_genes_counts = purrr::reduce(.data[["id_genes_counts"]],
                                        longest_common_substring),
        condition = unique(.data[["condition"]]),
        purrr::map_dfc(
          dplyr::select(dplyr::cur_data(), dplyr::any_of(get_gene_names())),
          sum
        )
      )
    )
}

#' Plot median and other quantile values for variance-to-mean ratios obtained
#' across all genes.
#'
#' Useful in comparison to a reference: if median and other quantile ratios are
#' higher or lower than in the reference, then one can draw a conclusion about a
#' relative variability among the samples in a given experimental condition.
#'
#' @param prepped_rd_input A data frame. The output of a call to
#'   [prep_rd_input()].
#'
#'
#' @return A ggplot.
#'
#' @family Arkady
#' @export
plot_var_mean_ratio <- function(prepped_rd_input) {
  checkmate::assert_data_frame(prepped_rd_input)
  checkmate::assert_names(names(prepped_rd_input), must.include = "condition")
  prepped_rd_input <- janitor::remove_constant(prepped_rd_input)
  plot_df <- prepped_rd_input %>%
    dplyr::group_by(.data[["condition"]]) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(dplyr::any_of(get_gene_names()),
                        names_to = "gene", values_to = "counts") %>%
    dplyr::group_by(.data[["condition"]], .data[["gene"]]) %>%
    dplyr::summarise(median=stats::median(.data[["counts"]]),
                     q25=stats::quantile(.data[["counts"]], 0.25),
                     q75=stats::quantile(.data[["counts"]], 0.75),
                     iqr = q75 - q25,
                     var=stats::var(.data[["counts"]]),
                     mean=mean(.data[["counts"]]),
                     v_m = .data[["var"]] / .data[["mean"]],
                     v_m_cat = if_else(v_m > 1, "greater", "less"),
                     .groups = "drop_last") %>%
    dplyr::summarise(
      median_v_m=stats::median(.data[["v_m"]], na.rm = TRUE),
      q75_v_m=stats::quantile(.data[["v_m"]], 0.75, na.rm=TRUE),
      q90_v2m=stats::quantile(.data[["v_m"]], 0.9, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(-"condition",
                        names_to = "statistic", values_to = "value") %>%
    dplyr::arrange(.data[["condition"]])
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = forcats::fct_reorder(.data[["condition"]],
                                          .data[["value"]]),
                 y = .data[["value"]], fill = .data[["condition"]])
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::facet_wrap("statistic", nrow = 1) +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::labs(x = "condition", y = "Variance to mean ratio")
}

#' Plots comparisons between conditions in an experiment.
#'
#' Actual experimental samples and positive controls, using all genes in the
#' count table or a subset of genes.
#'
#' @inheritParams plot_var_mean_ratio
#' @param gene_subset A character vector of genes in a subset of interest.
#'
#' @return A list of 6 ggplots.
#'
#'
#' @family Arkady
#' @export
plot_metrics_global <- function(prepped_rd_input, gene_subset = NULL) {
  checkmate::assert_data_frame(prepped_rd_input)
  checkmate::assert_names(
    names(prepped_rd_input),
    must.include = c("id_genes_counts", "condition")
  )
  checkmate::assert_character(gene_subset, min.chars = 1, null.ok = TRUE)
  if (!is.null(gene_subset)) {
    genes_outside_subset <- dplyr::setdiff(get_gene_names(), gene_subset)
    prepped_rd_input <- dplyr::select(prepped_rd_input,
                                      -dplyr::any_of(genes_outside_subset))
  }
  prepped_rd_input <- janitor::remove_constant(prepped_rd_input)
  no_controls <- dplyr::filter(
    prepped_rd_input,
    stringr::str_detect(id_genes_counts, "_PC|_NC|_NTC", negate = TRUE)
  )
  pcs_only <- dplyr::filter(
    prepped_rd_input,
    stringr::str_detect(id_genes_counts, "_PC")
  )
  n_missing_plot_df <- dplyr::mutate(
    no_controls,
    n_missing = rowSums(
      as.matrix(
        dplyr::select(no_controls, dplyr::any_of(get_gene_names()))
      ) == 0
    )
  )
  tot_counts_plot_df <- dplyr::mutate(
    no_controls,
    total_counts = rowSums(
      dplyr::select(no_controls, dplyr::any_of(get_gene_names()))
    )
  )
  q75_plot_df <- dplyr::mutate(
    no_controls,
    q75 = apply(
      as.matrix(
        dplyr::select(no_controls, dplyr::any_of(get_gene_names()))
      ),
      1,
      stats::quantile, 0.75
    )
  )
  n_missing_plot_df_pc <- dplyr::mutate(
    pcs_only,
    n_missing = rowSums(
      as.matrix(
        dplyr::select(pcs_only, dplyr::any_of(get_gene_names()))
      ) == 0
    )
  )
  tot_counts_plot_df_pc <- dplyr::mutate(
    pcs_only,
    total_counts = rowSums(
      dplyr::select(pcs_only, dplyr::any_of(get_gene_names()))
    )
  )
  q75_plot_df_pc <- dplyr::mutate(
    pcs_only,
    q75 = apply(
      as.matrix(
        dplyr::select(pcs_only, dplyr::any_of(get_gene_names()))
      ),
      1,
      stats::quantile, 0.75
    )
  )
  out <- list(
    ggstatsplot::ggbetweenstats(n_missing_plot_df,
                                x=condition, y=n_missing,
                                bf.message = FALSE, results.subtitle = FALSE,
                                centrality.plotting = FALSE),
    ggstatsplot::ggbetweenstats(tot_counts_plot_df,
                                x=condition, y=total_counts,
                                bf.message = FALSE, results.subtitle = FALSE,
                                centrality.plotting = FALSE),
    ggstatsplot::ggbetweenstats(q75_plot_df,
                                x=condition, y=q75,
                                bf.message = FALSE, results.subtitle = FALSE,
                                centrality.plotting = FALSE),
    ggplot2::ggplot(n_missing_plot_df_pc,
                    ggplot2::aes(x=condition, y=n_missing)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.5),
    ggplot2::ggplot(tot_counts_plot_df_pc,
                    ggplot2::aes(x=condition, y=total_counts)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.5),
    ggplot2::ggplot(q75_plot_df_pc,
                    ggplot2::aes(x=condition, y=q75)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.5)
  )
  print(patchwork::wrap_plots(out, nrow = 2))
  invisible(out)
}

#' Compute Pearson correlation between all pairs of rows.
#'
#' @param num_mat A numeric matrix with at leas 2 rows and 2 columns.
#'
#' @return A numeric vector.
#'
#' @noRd
pearson_all_row_pairs <- function(num_mat) {
  checkmate::assert_matrix(num_mat, mode = "numeric",
                           min.rows = 2, min.cols = 2)
  pairs <- utils::combn(seq_len(nrow(num_mat)), 2, simplify = FALSE)
  purrr::map_dbl(pairs, ~stats::cor(num_mat[.[1], ], num_mat[.[2], ]))
}

#' Violin plots for within-condition gene Pearson correlation coefficients.
#'
#' @inheritParams plot_metrics_global
#' @param qntl A number between 0 and 1. A quantile-based threshold below which
#'   the data is subset. The purpose of this is to improve the sensitivity of
#'   the reproducibility analysis by focusing on subsets of data with lower
#'   counts.
#'
#' @return A ggplot.
#'
#' @family Arkady
#' @export
plot_cond_paired_pearson <- function(prepped_rd_input,
                                     gene_subset = NULL, qntl=NULL){
  checkmate::assert_data_frame(prepped_rd_input)
  checkmate::assert_names(
    names(prepped_rd_input),
    must.include = c("id_genes_counts", "condition")
  )
  checkmate::assert_character(gene_subset, min.chars = 1, null.ok = TRUE)
  checkmate::assert_number(qntl, lower = 0, upper = 1, null.ok = TRUE)
  if (!is.null(gene_subset)) {
    genes_outside_subset <- dplyr::setdiff(get_gene_names(), gene_subset)
    prepped_rd_input <- dplyr::select(prepped_rd_input,
                                      -dplyr::any_of(genes_outside_subset))
  }
  prepped_rd_input <- janitor::remove_constant(prepped_rd_input)
  if (!is.null(qntl)) {
    gene_meds <- prepped_rd_input %>%
      dplyr::select(dplyr::any_of(get_df_gene_names())) %>%
      purrr::map_dbl(stats::median)
    qntl_thresh <- stats::quantile(row_meds, qntl)
    prepped_rd_input <- filter_genes(prep_rd_input,
                                     ~stats::median(.) <= qntl_thresh)
  }
  plot_df <- prepped_rd_input %>%
    dplyr::group_by(condition) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::summarise(
      dplyr::tibble(
        pearson = pearson_all_row_pairs(
          as.matrix(
            dplyr::select(dplyr::cur_data(), dplyr::any_of(get_gene_names()))
          )
        )
      ),
      .groups = "drop"
    )
  ggstatsplot::ggbetweenstats(plot_df, condition, pearson,
                              bf.message = FALSE, results.subtitle = FALSE,
                              centrality.plotting = FALSE,
                              package = "pals", palette = "alphabet")
}

#' Get the [mirmodels::deseq()] tables from pairs of conditions.
#'
#' Technical replicates in `genes_counts` should be collapsed with
#' [collapse_tech_reps()] prior to running this function.
#'
#' @inheritParams plot_var_mean_ratio
#'
#' @return A dataframe. A series of [mirmodels::deseq()] tables bound together
#'   with [dplyr::bind_rows()]. There are two extra columns `cond1` and `cond2`
#'   to record the conditions being compared.
#'
#' @family Arkady
#' @export
deseq_pairs <- function(prepped_rd_input) {
  checkmate::assert_data_frame(prepped_rd_input)
  checkmate::assert_names(
    names(prepped_rd_input),
    must.include = c("condition")
  )
  prepped_rd_input <- janitor::remove_constant(prepped_rd_input)
  conditions <- prepped_rd_input$condition %>%
    janitor::tabyl() %>%
    dplyr::filter(n > 1) %>%
    .[[1]] %>%
    sort()
  condition_pairs <- utils::combn(conditions, 2, simplify = FALSE)
  dsq_input_dfs <- purrr::map(
    condition_pairs,
    ~prepped_rd_input %>%
      dplyr::filter(condition %in% .env[[".x"]]) %>%
      dplyr::mutate(condition = factor(as.character(condition)))
  )
  gene_names <- get_df_gene_names(prepped_rd_input)
  dsq_outs <- purrr::map(dsq_input_dfs,
                         purrr::possibly(mirmodels::deseq,
                                         otherwise = data.frame()),
                         condition = "condition",
                         genes = gene_names, shrink = TRUE)
  keep <- purrr::map_lgl(dsq_outs, ~nrow(.) > 0)
  purrr::map2_dfr(
    condition_pairs[keep], dsq_outs[keep],
    ~dplyr::bind_cols(cond1 = .x[[1]], cond2 = .x[[2]], .y)
  ) %>%
    dplyr::arrange(padj, cond1, cond2, gene)
}

#' Obtain DESeq object from raw counts wide dataframe.
#'
#' Wide is samples in columns, genes in rows.
#'
#' @param data A wide datafame with raw counts that contains raw counts for
#'   samples of interest (may contain other samples as well, which will be
#'   filtered out).
#' @param design A dataframe with all of the variables specified in `compvar` as
#'   column names. Must have the same number of rows as `data` with samples in
#'   the same order.
#' @param compvar -  a character vector with the name(s) of the variable(s),
#'   which must be identical to the variable name(s) in the design dataframe.
#'
#' @return A DESeq2 object.
#'
#' @examples
#' rs_genes <- arrow::read_feather(
#'   system.file("extdata", "rs_genes_draw_3_4.feather", package = "mirmisc")
#' )
#' rs_meta <- arrow::read_feather(
#'   system.file("extdata", "rs_meta_draw_3_4.feather", package = "mirmisc")
#' )
#' deseq_wide(rs_genes, rs_meta, compvar = "meta_draw")
#'
#' @family Arkady
#' @export
deseq_wide <- function(data, design, compvar) {
  checkmate::assert_data_frame(data, min.cols = 1, min.rows = 1)
  checkmate::assert_data_frame(design, min.cols = 1, nrows = nrow(data))
  checkmate::assert_names(names(data), must.include = "mirvie_id")
  checkmate::assert_names(names(design), must.include = c("mirvie_id", compvar))
  checkmate::assert_character(data$mirvie_id, unique = TRUE)
  checkmate::assert_character(design$mirvie_id, unique = TRUE)
  mirvie_id_intersection <- intersect(data$mirvie_id, design$mirvie_id)
  data <- dplyr::filter(data, mirvie_id %in% mirvie_id_intersection)
  design <- dplyr::filter(design, mirvie_id %in% mirvie_id_intersection)
  df <- data %>%
    filter_genes(~max(.) > 0) %>%
    dplyr::select(dplyr::any_of(get_gene_names())) %>%
    as.data.frame() %>%
    magrittr::set_rownames(data$mirvie_id) %>%
    magrittr::set_colnames(convert_gene_names(names(.))) %>%
    t() %>%
    as.data.frame()
  for (cmpvr in compvar) {
    if (is.character(design[[cmpvr]]) || cmpvr == dplyr::last(compvar)) {
      design[[cmpvr]] <- as.factor(design[[cmpvr]])
    }
  }
  ddsIN <- DESeq2::DESeqDataSetFromMatrix(
    countData = df,
    colData = design,
    design = as.formula(paste("~", paste(compvar, collapse = "+")))
  )
  DESeq2::DESeq(ddsIN)
}


#' Obtain pairwise contrasts with shrunken log2FC from a DESeq  object.
#'
#' @param ddsx A DESeq object, output from [deseq_wide()].
#' @param compr A character vector for pairwise comparisons of interest, eg.,
#'   `c("4", "1")` from draw 4 with samples from draw 1.
#' @param alpha A number. The significance level.
#' @param shrink A flag. `DESeq2::lfcShrink()` the result? Default yes.
#'
#' @return A dataframe.
#'
#' @examples
#' rs_genes <- arrow::read_feather(
#'   system.file("extdata", "rs_genes_draw_3_4.feather", package = "mirmisc")
#' )
#' rs_meta <- arrow::read_feather(
#'   system.file("extdata", "rs_meta_draw_3_4.feather", package = "mirmisc")
#' )
#' ddsx <- deseq_wide(rs_genes, rs_meta, compvar = "meta_draw")
#' deseq_contrasts(ddsx, compr = c("4", "3"))
#'
#' @family Arkady
#' @export
deseq_contrasts <- function(ddsx, compr, alpha = 0.05, shrink = TRUE) {
  out <- DESeq2::results(ddsx,
                         contrast = c("meta_draw",compr[1],compr[2]),
                         alpha = alpha)
  if (shrink) {
    out <- DESeq2::lfcShrink(
      ddsx,
      contrast = c("meta_draw",compr[1],compr[2]), res=out, type="normal"
    )
  }
  as.data.frame(out)
}

#' Prep ranking metric gene input for [fgsea::fgsea()].
#'
#' @param x The dataframe with the contrasts from DESeq, results generated by
#'   [deseq_contrasts()].
#'
#' @return A data frame.
#'
#' @examples
#' rs_genes <- arrow::read_feather(
#'   system.file("extdata", "rs_genes_draw_3_4.feather", package = "mirmisc")
#' )
#' rs_meta <- arrow::read_feather(
#'   system.file("extdata", "rs_meta_draw_3_4.feather", package = "mirmisc")
#' )
#' ddsx <- deseq_wide(rs_genes, rs_meta, compvar = "meta_draw")
#' dsq_contrasts <- deseq_contrasts(ddsx, compr = c("4", "3"))
#' prep_gsea_input_lfc(dsq_contrasts)
#'
#' @family Arkady
#' @export
prep_gsea_input_lfc <- function(x) {
  x %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(diff = .data[["log2FoldChange"]]) %>%
    dplyr::select("gene", "diff") %>%
    tibble::column_to_rownames("gene")
}

#' Prepare gene set input for [fgsea::fgsea()].
#'
#' Categories of interest from
#' http://www.gsea-msigdb.org/gsea/msigdb/index.jsp.
#'
#' @param cetegory A string. Category of interest.
#'
#' @return A list object for fgsea.
#'
#' @examples
#' prep_gsea_input_gene_set("C8")
#'
#' @family Arkady
#' @export
prep_gsea_input_gene_set <- function(category) {
  dX <-  msigdbr::msigdbr(
    species = "Homo sapiens",
    category = category
  )[c("gs_name", "gene_symbol")]

  split(dX$gene_symbol,dX$gs_name)
}

#' Get results from fast Gene Set Enrichment Analysis (GSEA)
#'
#' Uses [fgsea::fgsea()] defaults.
#'
#' @param prepped_gsea_input_lfc Output of [prep_gsea_input_lfc()].
#' @param prepped_gsea_input_gene_set Output of [prep_gsea_input_gene_set()].
#'
#' @return A dataframe.
#'
#' @examples
#' rs_genes <- arrow::read_feather(
#'   system.file("extdata", "rs_genes_draw_3_4.feather", package = "mirmisc")
#' )
#' rs_meta <- arrow::read_feather(
#'   system.file("extdata", "rs_meta_draw_3_4.feather", package = "mirmisc")
#' )
#' ddsx <- deseq_wide(rs_genes, rs_meta, compvar = "meta_draw")
#' dsq_contrasts <- deseq_contrasts(ddsx, compr = c("4", "3"))
#' prepped_gsea_input_lfc <- prep_gsea_input_lfc(dsq_contrasts)
#' prepped_gsea_input_gene_set <- prep_gsea_input_gene_set("C8")
#' fgsea_basic(prepped_gsea_input_lfc, prepped_gsea_input_gene_set)
#'
#' @family Arkady
#' @export
fgsea_basic <- function(prepped_gsea_input_lfc, prepped_gsea_input_gene_set) {
  gene_list <- prepped_gsea_input_lfc %>%
    dplyr::pull(diff) %>%
    purrr::set_names(rownames(prepped_gsea_input_lfc)) %>%
    sort()
  fgsea::fgsea(prepped_gsea_input_gene_set, gene_list)
}

#' Plot top `n` most significantly enriched pathways in the upward or downward
#' direction.
#'
#' @param fgsea_out Output of [fgsea::fgsea()] or [fgsea_basic()].
#' @param padj_thresh Adjusted p-value threshold.
#' @param size_thresh number of genes to show.
#'
#' @return A [ggplot2::ggplot()].
#'
#' @examples
#' rs_genes <- arrow::read_feather(
#'   system.file("extdata", "rs_genes_draw_3_4.feather", package = "mirmisc")
#' )
#' rs_meta <- arrow::read_feather(
#'   system.file("extdata", "rs_meta_draw_3_4.feather", package = "mirmisc")
#' )
#' ddsx <- deseq_wide(rs_genes, rs_meta, compvar = "meta_draw")
#' dsq_contrasts <- deseq_contrasts(ddsx, compr = c("4", "3"))
#' prepped_gsea_input_lfc <- prep_gsea_input_lfc(dsq_contrasts)
#' prepped_gsea_input_gene_set <- prep_gsea_input_gene_set("C8")
#' fgs <- fgsea_basic(prepped_gsea_input_lfc, prepped_gsea_input_gene_set)
#' plot_down(fgs)
#' plot_up(fgs)
#'
#' @family Arkady
#' @export
plot_down <- function(fgsea_out, padj_thresh = 0.05, n = 20) {
  fgsea_out %>%
    dplyr::filter(size >= 5) %>%
    dplyr::mutate(core = lengths(leadingEdge)/size) %>%
    dplyr::filter(padj<padj_thresh & NES<0) %>%
    dplyr::slice_min(padj, n=n) %>%
    ggplot2::ggplot(
      ggplot2::aes(x=core, y=stats::reorder(pathway,core),
                   color=-log10(padj), size=NES)
    ) +
    ggplot2::geom_point() +
    ggplot2::labs(y="Gene sets", x="Fraction of gene set",
                  title=paste("Top", n, "donwregulated pathways")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=7),
                   plot.title = ggplot2::element_text(hjust=0.7),
                   aspect.ratio = 1.7)
}

#' @rdname plot_down
#'
#' @export
plot_up <- function(fgsea_out, padj_thresh = 0.05, n = 20) {
  fgsea_out %>%
    dplyr::filter(size >= 5) %>%
    dplyr::mutate(core = lengths(leadingEdge)/size) %>%
    dplyr::filter(padj<padj_thresh & NES>0) %>%
    dplyr::slice_min(padj, n=n) %>%
    ggplot2::ggplot(
      ggplot2::aes(x=core, y=stats::reorder(pathway,core),
                   color=-log10(padj), size=NES)
    ) +
    ggplot2::geom_point() +
    ggplot2::labs(y="Gene sets", x="Fraction of gene set",
                  title="Top-20 upregulated pathways")+
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=7),
                   plot.title = ggplot2::element_text(hjust=0.7),
                   aspect.ratio = 1.7)
}
