% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calcSD_Splex}
\alias{calcSD_Splex}
\title{Calculates CV per TMT_Set and LCMS_injection}
\usage{
calcSD_Splex(df, id, type = "log2_R")
}
\arguments{
\item{df}{The name of a primary data file. By default, it will be determined
automatically after matching the types of data and analysis with an
\code{id} among \code{c("pep_seq", "pep_seq_mod", "prot_acc", "gene")}. A
primary file contains normalized peptide or protein data and is among
\code{c("Peptide.txt", "Peptide_pVal.txt", "Peptide_impNA_pVal.txt",
"Protein.txt", "Protein_pVal.txt", "protein_impNA_pVal.txt")}. For analyses
require the fields of significance p-values, the \code{df} will be one of
\code{c("Peptide_pVal.txt", "Peptide_impNA_pVal.txt", "Protein_pVal.txt",
"protein_impNA_pVal.txt")}.}

\item{id}{Character string; one of \code{pep_seq}, \code{pep_seq_mod},
\code{prot_acc} and \code{gene}.}

\item{type}{Character string; the type of data for SD calculations.}
}
\description{
Calculates CV per TMT_Set and LCMS_injection
}
