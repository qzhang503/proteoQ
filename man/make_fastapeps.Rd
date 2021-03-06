% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbs.R
\name{make_fastapeps}
\alias{make_fastapeps}
\title{Make peptide sequences from fastas.}
\usage{
make_fastapeps(fasta_db, max_miss = 2, min_len = 1, max_len = 100)
}
\arguments{
\item{fasta_db}{Fasta database.}

\item{max_miss}{The maximum number of mis-cleavages per peptide sequence.}

\item{min_len}{Integer; the minimum length of peptides. Shorter peptides will
be excluded.}

\item{max_len}{Integer; the maximum length of peptides. Longer peptides will
be excluded.}
}
\description{
Make peptide sequences from fastas.
}
