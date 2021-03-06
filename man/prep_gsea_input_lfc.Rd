% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arkady.R
\name{prep_gsea_input_lfc}
\alias{prep_gsea_input_lfc}
\title{Prep ranking metric gene input for \code{\link[fgsea:fgsea]{fgsea::fgsea()}}.}
\usage{
prep_gsea_input_lfc(x)
}
\arguments{
\item{x}{The dataframe with the contrasts from DESeq, results generated by
\code{\link[=deseq_contrasts]{deseq_contrasts()}}.}
}
\value{
A data frame.
}
\description{
Prep ranking metric gene input for \code{\link[fgsea:fgsea]{fgsea::fgsea()}}.
}
\examples{
rs_genes <- arrow::read_feather(
  system.file("extdata", "rs_genes_draw_3_4.feather", package = "mirmisc")
)
rs_meta <- arrow::read_feather(
  system.file("extdata", "rs_meta_draw_3_4.feather", package = "mirmisc")
)
ddsx <- deseq_wide(rs_genes, rs_meta, compvar = "meta_draw")
dsq_contrasts <- deseq_contrasts(ddsx, compr = c("4", "3"))
prep_gsea_input_lfc(dsq_contrasts)

}
\seealso{
Other Arkady: 
\code{\link{check_sample_descs}()},
\code{\link{collapse_tech_reps}()},
\code{\link{deseq_contrasts}()},
\code{\link{deseq_pairs}()},
\code{\link{deseq_wide}()},
\code{\link{fgsea_basic}()},
\code{\link{plot_cond_paired_pearson}()},
\code{\link{plot_down}()},
\code{\link{plot_metrics_global}()},
\code{\link{plot_var_mean_ratio}()},
\code{\link{prep_gsea_input_gene_set}()},
\code{\link{prep_rd_input}()},
\code{\link{reconcile_names}()}
}
\concept{Arkady}
