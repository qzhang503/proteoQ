% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_formalArgs}
\alias{check_formalArgs}
\title{Checks conflicts in function arguments.}
\usage{
check_formalArgs(w = anal_prnTrend, f = cmeans, excludes = NULL)
}
\arguments{
\item{w}{Expression; the wrapper function.}

\item{f}{Expression; the function to be wrapped.}

\item{excludes}{Character string; arguments from \code{f} to be excluded for
checking.}
}
\description{
Conflicts between wrapper and function to be wrapped. Note that \code{f}
contains function name but not package name.
}
\details{
if `method` in the formalArgs and "m" in the "..." of `foo`, 
foo(m = 2) without `method` interpreted as `method`.
Call foo(method = kmeans, m =2) to avoid such ambiguity.
}
