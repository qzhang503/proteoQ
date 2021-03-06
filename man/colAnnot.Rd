% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{colAnnot}
\alias{colAnnot}
\title{Sets up the column annotation in heat maps}
\usage{
colAnnot(annot_cols = NULL, sample_ids = NULL, annot_colnames = NULL)
}
\arguments{
\item{annot_cols}{A character vector of column keys in \code{expt_smry.xlsx}.
The values under the selected keys will be used to color-code sample IDs on
the top of the indicated plot. The default is NULL without column
annotation.}

\item{sample_ids}{A character vector containing the sample IDs for an ascribing analysis.}

\item{annot_colnames}{A character vector of replacement name(s) to
\code{annot_cols}. The default is NULL without name replacement.}
}
\description{
Sets up the column annotation in heat maps
}
