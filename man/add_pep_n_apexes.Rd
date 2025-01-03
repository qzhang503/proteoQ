% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{add_pep_n_apexes}
\alias{add_pep_n_apexes}
\title{Add column \code{pep_n_apexes}}
\usage{
add_pep_n_apexes(df, key = "pep_seq_mod", max_n_apexes = 2L)
}
\arguments{
\item{df}{A PSM sub-table with only one \code{RAW_File}.}

\item{key}{The column where unique \code{pep_apex_ret}s were counted.}

\item{max_n_apexes}{The maximum number of apexes for consideration.}
}
\description{
Add column \code{pep_n_apexes}
}
