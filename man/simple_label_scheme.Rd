% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{simple_label_scheme}
\alias{simple_label_scheme}
\title{Simplifies label schemes from \code{label_scheme_full}}
\usage{
simple_label_scheme(dat_dir, label_scheme_full)
}
\arguments{
\item{dat_dir}{A character string to the working directory. The default is to
match the value under the global environment.}

\item{label_scheme_full}{The label_scheme with the probable inclusion of
different LCMS_inj under the same TMT_Set.}
}
\description{
Removes duplicated sample entries under different LC/MS injections.
}
