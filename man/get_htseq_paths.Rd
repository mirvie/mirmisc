% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collect-counts.R
\name{get_htseq_paths}
\alias{get_htseq_paths}
\title{Get the paths to \verb{htseq/} directories for a given cohort.}
\usage{
get_htseq_paths(base_dir = "/mnt/storage/Cohorts", cohort_code)
}
\arguments{
\item{base_dir}{A string. The path to a directory that the cohort directories
live under. The cohort directories have name structure \verb{###_XY} where \verb{#}
is a digit and \code{XY} is the cohort code. For example, \verb{007_RS}.}

\item{cohort_code}{A string with exactly two characters. E.g. \code{"RS"}.}
}
\value{
A character vector of paths.
}
\description{
The \verb{*_counts.txt} files live in directories called \verb{htseq/}. This function
helps you to find all such directories for a given cohort.
}
\examples{
\dontrun{
get_htseq_paths(cohort_code = "RS")
}

}
