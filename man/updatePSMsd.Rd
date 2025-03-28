% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{updatePSMsd}
\alias{updatePSMsd}
\title{Update cross-sample PSM standard deviation}
\usage{
updatePSMsd(dfs, filelist, key = "pep_seq_modz", qt_sd = 0.97, cleanup = FALSE)
}
\arguments{
\item{dfs}{A list of PSM tables.}

\item{filelist}{The file names corresponding to \code{dfs}.}

\item{key}{The key for data grouping.}

\item{qt_sd}{The quantile cut-off in retention time SD values.}

\item{cleanup}{Logical; nullify data of retention outliers or not.}
}
\description{
Update cross-sample PSM standard deviation
}
