% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{spline_coefs}
\alias{spline_coefs}
\title{calculates the median-centered \eqn{x}'s positions by spline points}
\usage{
spline_coefs(
  df,
  label_scheme,
  label_scheme_fit,
  is_prot_lfq = FALSE,
  method_align = "MC",
  ...
)
}
\arguments{
\item{df}{A data frame.}

\item{label_scheme}{The metadata of label scheme.}

\item{label_scheme_fit}{the subset of \code{label_scheme} used in fitting.}

\item{is_prot_lfq}{Logical; is protein LFQ data or not. About half of the
protein intensity values can be missing with LFQ and imputed with small
values. The typically causes a bimodality in protein log2FC distributions
and need to be handled especially at \code{method_align = "MC"}.}

\item{method_align}{Character string indicating the method in aligning
\code{log2FC} across samples. \code{MC}: median-centering; \code{MGKernel}:
the kernel density defined by multiple Gaussian functions
(\code{\link[mixtools]{normalmixEM}}). At the \code{MC} default, the ratio
profiles of each sample will be aligned in that the medians of the
\code{log2FC} are zero. At \code{MGKernel}, the ratio profiles of each
sample will be aligned in that the \code{log2FC} at the maximums of kernel
density are zero.}

\item{...}{Additional parameters, including slice_dots.}
}
\description{
Only used with \code{method_align = MC}. The values of \eqn{x}'s for sample
IDs indicated in \code{label_scheme_fit} will be updated for subsequent
adjustment in median centering. The rest of \eqn{x}'s stays the same,
corresponding to no further changes in center positions.
}
