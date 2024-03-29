% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{rm_cols_mqpsm}
\alias{rm_cols_mqpsm}
\title{Removes columns from PSM table}
\usage{
rm_cols_mqpsm(df = NULL, group_psm_by = "pep_seq", set_idx = 1L)
}
\arguments{
\item{df}{PSM data}

\item{group_psm_by}{A character string specifying the method in PSM grouping.
At the \code{pep_seq} default, descriptive statistics will be calculated
based on the same \code{pep_seq} groups. At the \code{pep_seq_mod}
alternative, peptides with different variable modifications will be treated
as different species and descriptive statistics will be calculated based on
the same \code{pep_seq_mod} groups.}

\item{set_idx}{Numeric.  The index of a multiplex TMT experiment in metadata
files such as \code{label_scheme.xlsx} and \code{frac_scheme.xlsx}.}
}
\description{
Removes columns from PSM table
}
