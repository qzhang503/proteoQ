% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{msubsetBrukerMGF}
\alias{msubsetBrukerMGF}
\title{Parallel subsetBrukerMGF}
\usage{
msubsetBrukerMGF(
  filepath,
  begin_offset = 5L,
  charge_offset = 5L,
  topn_ms2ions = Inf,
  n_cores = 1L
)
}
\arguments{
\item{filepath}{A file path to MGF.}

\item{begin_offset}{The number of lines before a BEGIN line.}

\item{charge_offset}{The number lines after a BEGIN line to a following
CHARGE line.}

\item{topn_ms2ions}{Top-n MS2 ions to be retained.}

\item{n_cores}{The number of CPU cores.}
}
\description{
Parallel subsetBrukerMGF
}
