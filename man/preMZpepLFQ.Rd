% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{preMZpepLFQ}
\alias{preMZpepLFQ}
\title{Pre-processing of Mzion PSM tables}
\usage{
preMZpepLFQ(
  df = NULL,
  group_psm_by = "pep_seq_mod",
  group_pep_by = "prot_acc",
  dat_dir = NULL,
  set_idx = 1L,
  injn_idx = 1L,
  tmt_plex = 1L,
  rm_allna = FALSE
)
}
\arguments{
\item{df}{A data frame.}

\item{group_psm_by}{A character string specifying the method in PSM grouping.
At the \code{pep_seq} default, descriptive statistics will be calculated
based on the same \code{pep_seq} groups. At the \code{pep_seq_mod}
alternative, peptides with different variable modifications will be treated
as different species and descriptive statistics will be calculated based on
the same \code{pep_seq_mod} groups.}

\item{group_pep_by}{A character string specifying the method in peptide
grouping. At the \code{prot_acc} default, descriptive statistics will be
calculated based on the same \code{prot_acc} groups. At the \code{gene}
alternative, proteins with the same gene name but different accession
numbers will be treated as one group.}

\item{dat_dir}{A working directory.}

\item{set_idx}{A set index.}

\item{injn_idx}{An injection index.}

\item{rm_allna}{Logical; if TRUE, removes data rows that are exclusively NA
across ratio columns of \code{log2_R126} etc. The setting also applies to
\code{log2_R000} in LFQ.}

\item{TMT_plex}{The multiplexity of TMT.}
}
\description{
Pre-processing of Mzion PSM tables
}
