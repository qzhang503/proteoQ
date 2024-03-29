% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{prep_fraction_scheme}
\alias{prep_fraction_scheme}
\title{Loads the information of analyte prefractionation}
\usage{
prep_fraction_scheme(
  dat_dir = NULL,
  expt_smry = "expt_smry.xlsx",
  frac_smry = "frac_smry.xlsx"
)
}
\arguments{
\item{dat_dir}{A character string to the working directory. The default is to
match the value under the global environment.}

\item{expt_smry}{A character string to a \code{.xlsx} file containing the
metadata of TMT or LFQ experiments. The default is \code{expt_smry.xlsx}.}

\item{frac_smry}{A character string to a \code{.xlsx} file containing peptide
fractionation summary. The default is \code{frac_smry.xlsx}.}
}
\description{
Loads the information of analyte prefractionation
}
