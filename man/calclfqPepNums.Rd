% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptable.R
\name{calclfqPepNums}
\alias{calclfqPepNums}
\title{Calculates log2FC of peptides based on LFQ intensity.}
\usage{
calclfqPepNums(df, omit_single_lfq = FALSE)
}
\arguments{
\item{df}{A data frame.}

\item{omit_single_lfq}{Logical; if TRUE, omits LFQ entries with single
measured values across all samples. The default is FALSE.}
}
\description{
With LFQ, values of \code{log2_R000 (...)} and \code{N_log2_R000 (...)} in
\code{df_num} are not yet filled after \link{spreadPepNums}. This utility
calculates the ratio using the values of LFQ intensity.
}
\details{
For SILAC, intensity of heave, light etc. have been spread into multiple
columns and thus the log2Ratios can be calculated like regular LFQ.
}
