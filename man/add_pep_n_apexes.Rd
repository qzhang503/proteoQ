% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{add_pep_n_apexes}
\alias{add_pep_n_apexes}
\title{Add column \code{pep_n_apexes}}
\usage{
add_pep_n_apexes(df, key = "pep_seq_modz", max_n_apexes = 2L, cleanup = TRUE)
}
\arguments{
\item{df}{A PSM sub-table with only one \code{RAW_File}.}

\item{key}{The column where unique \code{pep_apex_ret}s were counted.}

\item{max_n_apexes}{The maximum number of apexes under the same \code{key}
for consideration. Note that the clean-up does not apply to the lengths of
apexes indicated under column \code{apx_ps}.}

\item{cleanup}{Logical; clean up data or not.}
}
\description{
Add column \code{pep_n_apexes}
}
