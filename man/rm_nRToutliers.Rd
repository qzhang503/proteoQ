% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptable.R
\name{rm_nRToutliers}
\alias{rm_nRToutliers}
\title{Remove retention-time outliers for more than two values}
\usage{
rm_nRToutliers(mx, my, mt, ms, err)
}
\arguments{
\item{mx}{A multi-row matrix of m-over-z values.}

\item{my}{A multi-row matrix of apex intensity values.}

\item{mt}{A multi-row matrix of apex retention times.}

\item{ms}{A multi-row matrix of apex scan numbers.}

\item{err}{The tolerance in retention times.}
}
\description{
Remove retention-time outliers for more than two values
}
