% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qlfq.R
\name{mpepLFQ}
\alias{mpepLFQ}
\title{LFQ of peptides}
\usage{
mpepLFQ(
  basenames,
  are_refs,
  are_smpls,
  xsmat,
  ysmat,
  tsmat,
  ssmat,
  mmat,
  pmat,
  mbr_peps,
  sp_centers,
  pri_sps = 0L,
  new_na_species = ".other",
  rt_tol = 30,
  mbr_ret_tol = 30,
  step = 0.005,
  lfq_mbr = TRUE,
  err_log2r = 0.25
)
}
\arguments{
\item{basenames}{The base names of a peptide table files
(\code{TMTset[i]_LCMSinj[j]_Peptide_N.txt}).}

\item{are_refs}{Logical; are references or not.}

\item{are_smpls}{Logical; are samples or not.}

\item{xsmat}{A matrix of m-over-z values.}

\item{ysmat}{A matrix of apex intensities.}

\item{tsmat}{A matrix of apex retention times.}

\item{ssmat}{A matrix of apex scan numbers.}

\item{mmat}{A logical matrix of MBR status.}

\item{pmat}{A confounding matrix or vectors of peptides shared across
species.}

\item{mbr_peps}{A vector of \code{pep_seq_modz} sequences in the universe.}

\item{sp_centers}{The centers of log2FC for each species; names: species,
values: log2FC.}

\item{pri_sps}{A primary species.}

\item{new_na_species}{A replace value for NA species.}

\item{rt_tol}{Error tolerance in retention times.}

\item{mbr_ret_tol}{The tolerance in MBR retention time in seconds.}

\item{step}{The step size in binning retention times.}

\item{lfq_mbr}{Logical; perform LFQ-MBR or not.}

\item{err_log2r}{Not yet used. Error tolerance in log2FC.}
}
\description{
Peptides are shared by species.
}
