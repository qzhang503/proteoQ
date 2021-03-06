% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbs.R
\name{calc_monopep}
\alias{calc_monopep}
\title{Calculates the mono-isotopic mass of a peptide sequence.}
\usage{
calc_monopep(
  aa_seq,
  aa_masses,
  include_insource_nl = FALSE,
  maxn_vmods_per_pep = 5,
  maxn_sites_per_vmod = 3,
  digits = 5
)
}
\arguments{
\item{aa_seq}{Character string; a peptide sequences with one-letter
representation of amino acids.}

\item{aa_masses}{A named list containing the (mono-isotopic) masses of amino
acid residues.}

\item{include_insource_nl}{Logical Logical; if TRUE, includes MS1 precursor
  masses with the losses of neutral species prior to MS2 fragmentation. The
  default is FALSE.

  Note that there is no combination of neutral losses for fixed modifications
  at the precursor levels. Changes from fixed to variable modifications for
  complete combinations.}

\item{maxn_vmods_per_pep}{The maximum number of variable modifications per
peptide.}

\item{maxn_sites_per_vmod}{Integer; the maximum number of combinatorial
variable modifications per site in a per peptide sequence.}

\item{digits}{Integer; the number of decimal places to be used.}
}
\description{
Typically coupled to \link{subpeps_by_vmods} for automatic dispatching of
peptide sequences by sets of fixed and variable modifications. For manual
calculations, uses \link{calc_monopeptide}.
}
