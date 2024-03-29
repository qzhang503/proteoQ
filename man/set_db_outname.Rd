% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/go.R
\name{set_db_outname}
\alias{set_db_outname}
\title{Helper to set the output file name for a data base}
\usage{
set_db_outname(filename = NULL, species = "human", signature)
}
\arguments{
\item{filename}{Character string; An output file name. At the \code{NULL}
default, the name will be determined automatically at a given
\code{species}; i.e., \code{msig_hs.rds} for \code{human} data. The file is
saved as a \code{.rds} object for uses with \code{\link{prnGSPA}}.}

\item{species}{Character string; the name of a species for the
\emph{conveninent} preparation of \code{MSig} data bases. The species
available for the convenience feature is in one of \code{c("human", "mouse",
"rat")} with \code{"human"} being the default. The argument is not required
for other species; instead, users will provide values under arguments
\code{ortho_mart} for the lookup of orthologs to human.}

\item{signature}{A character string, i.e. "go" for uses in an output filename.}
}
\description{
Helper to set the output file name for a data base
}
