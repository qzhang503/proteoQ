% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{aggrTopn}
\alias{aggrTopn}
\title{Sums top_n}
\usage{
aggrTopn(f)
}
\arguments{
\item{f}{A function for data summarization.}
}
\description{
\code{aggrTopn} summarizes \code{log2FC} and \code{intensity} by the
descriptive statistics of \code{c("mean", "median", "sum")}. Note the
difference to \link{TMT_top_n}, which uses mean statistics.
}
\examples{
\donttest{df_num <- aggrTopn(sum)(df, prot_acc, 3, na.rm = TRUE)}
}
