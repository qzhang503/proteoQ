% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptable.R
\name{calc_tmt_prnnums}
\alias{calc_tmt_prnnums}
\title{Helper: calculates the TMT log2FC and reporter-ion intensity of proteins}
\usage{
calc_tmt_prnnums(
  df,
  use_unique_pep = TRUE,
  id = "prot_acc",
  method_pep_prn = "median",
  use_spec_counts = FALSE
)
}
\arguments{
\item{df}{A data frame.}

\item{use_unique_pep}{Logical. If TRUE, only entries that are \code{TRUE} or
equal to \code{1} under the column \code{pep_isunique} in \code{Peptide.txt}
will be used, for summarizing the \code{log2FC} and the \code{intensity} of
peptides into protein values. The default is to use unique peptides only.
For \code{MaxQuant} data, the levels of uniqueness are according to the
\code{pep_unique_by} in \code{\link{normPSM}}. The argument currently do
nothing to \code{Spectrum Mill} data where both unique and shared peptides
will be kept.}

\item{id}{Always "prot_acc".}

\item{use_spec_counts}{Logical; If TRUE, uses spectrum counts for
quantitation with Mascot or Mzion outputs.}
}
\description{
Helper: calculates the TMT log2FC and reporter-ion intensity of proteins
}
