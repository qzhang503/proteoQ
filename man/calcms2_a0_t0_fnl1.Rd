% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sitecombi.R
\name{calcms2_a0_t0_fnl1}
\alias{calcms2_a0_t0_fnl1}
\title{Helper in calculating MS2 masses.}
\usage{
calcms2_a0_t0_fnl1(
  nl_combi,
  aas,
  aa_masses,
  ntmod,
  ctmod,
  type_ms2ions,
  digits
)
}
\arguments{
\item{nl_combi}{Lists of combinations of neutral losses for corresponding
\code{vmods_combi}. Each list contains a table where each column
corresponds to a set of neutral loss. The first column corresponds to the
combination without NLs.}

\item{aas}{\code{aa_seq} split in a sequence of LETTERS.}

\item{aa_masses}{A named list containing the (mono-isotopic) masses of amino
acid residues.}

\item{ntmod}{The parsed \emph{variable} \code{N-term} modifications from
\code{aa_masses}.}

\item{ctmod}{The parsed \emph{variable} \code{C-term} modifications from
\code{aa_masses}.}

\item{type_ms2ions}{Character; the type of
\href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
and y-ions.}

\item{digits}{Integer; the number of decimal places to be used.}
}
\description{
(6) "amods- tmod- vnl- fnl+".
}
\details{
Note that no need to update the \code{aas} (since `amods-`). The utility runs
column-wisely through the table of \code{nl_combi} (fixed) for corresponding
masses of residues.
}
