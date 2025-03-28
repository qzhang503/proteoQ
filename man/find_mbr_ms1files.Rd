% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptable.R
\name{find_mbr_ms1files}
\alias{find_mbr_ms1files}
\title{Find MBR MS1 files}
\usage{
find_mbr_ms1files(dat_dir, n_files, abort = FALSE)
}
\arguments{
\item{dat_dir}{The working directory}

\item{n_files}{The number of files.}

\item{abort}{Logical; to abort the run or not. TRUE at \link{mergePep};
otherwise, \code{df} is unique by pep_seq_modz and need to be collapsed to
pep_seq_mod.}
}
\description{
Find MBR MS1 files
}
