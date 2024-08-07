% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{smbr_psms}
\alias{smbr_psms}
\title{Single MBR.}
\usage{
smbr_psms(dfqi, dfci, cols_mbr, cols_dfq, mbr_ret_tol = 25)
}
\arguments{
\item{dfqi}{An i-th entry of dfq.}

\item{dfci}{An i-th entry of dfc.}

\item{cols_mbr}{Columns for MBR information borrowing.}

\item{cols_dfq}{Columns for copying from dfq.}

\item{mbr_ret_tol}{Retention time tolerance (in seconds) for LFQ-MBR.}
}
\description{
At the same \code{group_psm_by}.
}
