% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qlfq.R
\name{haddApexRTs_oneset}
\alias{haddApexRTs_oneset}
\title{Helper of adding apex retention times for one PSM table}
\usage{
haddApexRTs_oneset(
  file_psm = NULL,
  files_ms1 = NULL,
  dat_dir = NULL,
  path_ms1 = NULL,
  max_n_apexes = 2L,
  from = 115,
  step = 1e-05,
  rt_size = 240,
  rt_margin = 480,
  max_rt_delta = 240,
  yco = 100,
  min_y = 2e+06,
  n_cores = 1L
)
}
\arguments{
\item{file_psm}{The file names of a PSM table, e.g.,
TMTSet1_LCMSinj1_PSM_N.txt.}

\item{files_ms1}{The file name of an \code{ms1full_.rds}.}

\item{dat_dir}{The working directory.}

\item{path_ms1}{The file path to \code{ms1full_.rds}.}

\item{max_n_apexes}{The maximum number of apexes for consideration.}

\item{from}{A starting mass.}

\item{step}{A step size.}

\item{rt_size}{The width of each LC retention times in seconds.}

\item{rt_margin}{The bracketing margin before and after an LC retention time
window.}

\item{max_rt_delta}{The maximum allowance in retention-time difference
between the identifying MS2 event and apex.}

\item{yco}{The cut-off in intensities.}

\item{min_y}{The minimum peak area.}

\item{n_cores}{The number of parallel processes calling the current utility.}
}
\description{
Inputs are at the levels of one TMTSet_LCMSInj.txt file.
}
