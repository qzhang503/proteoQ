% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{flatten_pepouts}
\alias{flatten_pepouts}
\title{Flattens the outputs of ion matches.}
\usage{
flatten_pepouts(data, outcol = "matches")
}
\arguments{
\item{data}{Nested list of table from ion matches. Each table contains
columns \code{theo} and \code{expt}.}

\item{outcol}{The output column name.}
}
\description{
Flattens the outputs of ion matches.
}
\examples{
\donttest{
x <- flatten_pepouts(res$matches, "matches")
scores <- flatten_pepouts(scores, "pep_score")
}
}
