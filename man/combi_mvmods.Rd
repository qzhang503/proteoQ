% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sitecombi.R
\name{combi_mvmods}
\alias{combi_mvmods}
\title{Finds the combinations of variable modifications (multiple sites)}
\usage{
combi_mvmods(
  amods,
  ntmod,
  ctmod,
  famods,
  fntmod,
  fctmod,
  aas,
  aa_masses,
  maxn_vmods_per_pep = 5,
  maxn_sites_per_vmod = 3,
  maxn_vmods_sitescombi_per_pep = 32,
  digits = 5
)
}
\arguments{
\item{amods}{\code{Anywhere} variable modifications.}

\item{ntmod}{The parsed \emph{variable} \code{N-term} modifications from
\code{aa_masses}.}

\item{ctmod}{The parsed \emph{variable} \code{C-term} modifications from
\code{aa_masses}.}

\item{famods}{\code{Anywhere} fixed modifications.}

\item{fntmod}{The parsed \emph{fixed} \code{N-term} modifications from
\code{aa_masses}.}

\item{fctmod}{The parsed \emph{fixed} \code{C-term} modifications from
\code{aa_masses}.}

\item{aas}{\code{aa_seq} split in a sequence of LETTERS.}

\item{aa_masses}{A named list containing the (mono-isotopic) masses of amino
acid residues.}

\item{maxn_vmods_per_pep}{The maximum number of variable modifications per
peptide.}

\item{maxn_sites_per_vmod}{Integer; the maximum number of combinatorial
variable modifications per site in a per peptide sequence.}

\item{maxn_vmods_sitescombi_per_pep}{Integer; the maximum number of
combinatorial variable modifications per peptide sequence.}

\item{digits}{Integer; the number of decimal places to be used.}
}
\value{
Lists by residues in \code{amods}.
}
\description{
For all the \code{Anywhere} modifications specified in \code{amods}.
}
