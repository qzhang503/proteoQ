% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{rm_lowIntPSMs}
\alias{rm_lowIntPSMs}
\title{Remove PSM spikes by a group of \code{pep_seq_modz}}
\usage{
rm_lowIntPSMs(df, yfrac = 0.02, keep_best_score = TRUE)
}
\arguments{
\item{df}{A data frame for a \code{pep_seq_modz} at a given charge state.}

\item{yfrac}{A factor of intensity fraction. Peaks with intensities \eqn{le}
the fraction of the base peak will be ignored.}

\item{keep_best_score}{Logical; if TRUE, keeps the best-score PSM under a
\code{pep_seq_modz}.}
}
\description{
Remove PSM spikes by a group of \code{pep_seq_modz}
}
