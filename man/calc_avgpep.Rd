% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbs.R
\name{calc_avgpep}
\alias{calc_avgpep}
\title{Averaged peptide molecular weight (for proteins)}
\usage{
calc_avgpep(aa_seq = NULL, digits = 4L)
}
\arguments{
\item{aa_seq}{Character string; a peptide sequences with one-letter
representation of amino acids.}

\item{digits}{Integer; the number of decimal places to be used.}
}
\description{
Calculates the average molecular weight of a polypeptide (neutral species).
}
\examples{
\donttest{
calc_avgpep("AAIDWFDGKEFSGNPIK")
}

}
