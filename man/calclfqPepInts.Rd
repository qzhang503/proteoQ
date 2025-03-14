% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptable.R
\name{calcLFQPepInts}
\alias{calcLFQPepInts}
\title{Calculates \code{pep_tot_int}, \code{pep_razor_int} and
\code{pep_unique_int} in LFQ}
\usage{
calcLFQPepInts(df, filelist, group_psm_by)
}
\arguments{
\item{df}{A data frame.}

\item{filelist}{A list of individual peptide tables.}

\item{group_psm_by}{A character string specifying the method in PSM grouping.
At the \code{pep_seq} default, descriptive statistics will be calculated
based on the same \code{pep_seq} groups. At the \code{pep_seq_mod}
alternative, peptides with different variable modifications will be treated
as different species and descriptive statistics will be calculated based on
the same \code{pep_seq_mod} groups.}
}
\description{
Calculates \code{pep_tot_int}, \code{pep_razor_int} and
\code{pep_unique_int} in LFQ
}
