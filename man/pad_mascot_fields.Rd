% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{pad_mascot_fields}
\alias{pad_mascot_fields}
\title{Pads Mascot PSM exports}
\usage{
pad_mascot_fields(df = NULL)
}
\arguments{
\item{df}{An intermediate PSM table from Mascot export.}
}
\description{
Intensity or ratio columns are handled during channel padding and thus
excluded.
}
