% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptable.R
\name{impPepNA}
\alias{impPepNA}
\title{Imputes NA values in a peptide table}
\usage{
impPepNA(df, fold = 50, is_intensity = TRUE)
}
\arguments{
\item{df}{A intensity data frame (not matrix) of peptides.}

\item{fold}{The fold difference to mean.}

\item{is_intensity}{Logical; is intensity data or not (log2FC).}
}
\description{
Imputes NA values in a peptide table
}
