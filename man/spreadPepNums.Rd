% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptable.R
\name{spreadPepNums}
\alias{spreadPepNums}
\title{Spreads peptide numbers.}
\usage{
spreadPepNums(df, basenames, group_psm_by, ok_mbr = FALSE)
}
\arguments{
\item{df}{A data frame of peptide table.}

\item{basenames}{Names of peptide table files.}

\item{group_psm_by}{A character string specifying the method in PSM grouping.
At the \code{pep_seq} default, descriptive statistics will be calculated
based on the same \code{pep_seq} groups. At the \code{pep_seq_mod}
alternative, peptides with different variable modifications will be treated
as different species and descriptive statistics will be calculated based on
the same \code{pep_seq_mod} groups.}

\item{ok_mbr}{Logical Capable of MBR or not.}
}
\description{
Spreads fields of numeric values: sd_log2_R, log2_R, log2_R, I, N_I by TMT
sets.
}
\details{
Also works for LFQ as each sample corresponds to a TMT set.

For single SILAC sample, the values of log2Ratios spreads into
\emph{MULTIPLE} columns of heavy, light etc. Despite, log2Ratios remains NA,
just like regular single-sample LFQ. The log2Ratios will be later calculated
with \link{calcLFQPepNums} that are based on intensity values.
}
