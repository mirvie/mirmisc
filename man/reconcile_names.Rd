% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arkady.R
\name{reconcile_names}
\alias{reconcile_names}
\title{Match sample description sheet names to count file names.}
\usage{
reconcile_names(sample_descs, genes_counts)
}
\arguments{
\item{sample_descs}{A data frame of sample descriptions. The \code{sampleName}
column holds the sample names.}

\item{genes_counts}{A data frame of gene counts. The first column is \code{gene}
and the rest are sample names.}
}
\value{
A data frame.
}
\description{
Given a data frame of sample descriptions from the lab team and a data frame
from a \code{genes_counts.txt} file from a sequencing run, match the sample names
and return a data frame of sample descriptions detailing the matching.
}
\details{
Name matching is done by first replacing all whitespace and underscores with
hyphens and then asserting that the \code{sampleName} must a substring of the
column name. The matching is done via a pair of columns \code{id_sample_descs}
and \code{id_genes_counts} in the output.
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
\code{\link{prep_gsea_input_lfc}()},
\code{\link{prep_rd_input}()}
}
\concept{Arkady}
