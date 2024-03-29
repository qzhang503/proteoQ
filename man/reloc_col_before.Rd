% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{reloc_col_before}
\alias{reloc_col_before}
\title{Relocates column "to_move" immediately before anchor column "col_after".}
\usage{
reloc_col_before(df = NULL, to_move = "before_anchor", col_after = "anchor")
}
\arguments{
\item{df}{The original data frame.}

\item{to_move}{The column to be moved.}

\item{col_after}{The anchor column to which the \code{to_move} will be moved
before.}
}
\description{
The same as \code{reloc_col}.
}
