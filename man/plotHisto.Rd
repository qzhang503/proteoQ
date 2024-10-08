% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histo.R
\name{plotHisto}
\alias{plotHisto}
\title{Plots histograms}
\usage{
plotHisto(
  df = NULL,
  id,
  label_scheme_sub,
  scale_log2r = TRUE,
  complete_cases = FALSE,
  cut_points,
  show_curves = FALSE,
  show_vline = FALSE,
  scale_y = TRUE,
  filepath = NULL,
  filename = NULL,
  theme = NULL,
  ...
)
}
\arguments{
\item{df}{The name of a primary data file. By default, it will be determined
automatically after matching the types of data and analysis with an
\code{id} among \code{c("pep_seq", "pep_seq_mod", "prot_acc", "gene")}. A
primary file contains normalized peptide or protein data and is among
\code{c("Peptide.txt", "Peptide_pVal.txt", "Peptide_impNA_pVal.txt",
"Protein.txt", "Protein_pVal.txt", "protein_impNA_pVal.txt")}. For analyses
require the fields of significance p-values, the \code{df} will be one of
\code{c("Peptide_pVal.txt", "Peptide_impNA_pVal.txt", "Protein_pVal.txt",
"protein_impNA_pVal.txt")}.}

\item{id}{Character string; one of \code{pep_seq}, \code{pep_seq_mod},
\code{prot_acc} and \code{gene}.}

\item{label_scheme_sub}{A data frame. Subset entries from \code{label_scheme}
for selected samples.}

\item{scale_log2r}{Logical; if TRUE, adjusts \code{log2FC} to the same scale
of standard deviation across all samples. The default is TRUE. At
\code{scale_log2r = NA}, the raw \code{log2FC} without normalization will
be used.}

\item{complete_cases}{Logical; if TRUE, only cases that are complete with no
missing values will be used. The default is FALSE.}

\item{cut_points}{A named, numeric vector defines the cut points (knots) in
histograms. The default is \code{cut_points = c(mean_lint = NA)} where the
cut points correspond to the quantile values under column \code{mean_lint}
(mean log10 intensity) of input data. Values of \code{log2FC} will be then
binned from \eqn{-Inf} to \eqn{Inf} according to the cut points. To disable
data binning, set \code{cut_points = Inf} or \code{-Inf}. The binning of
\code{log2FC} can also be achieved through a different numeric column,
e.g., \code{cut_points = c(prot_icover = seq(.25, .75, .25))}. See also
\code{\link{mergePep}} for data alignment with binning.}

\item{show_curves}{Logical; if TRUE, shows the fitted curves. At the TRUE
default, the curve parameters are based on the latest call to
\code{\link{standPep}} or \code{\link{standPrn}} with \code{method_align =
MGKernel}. This feature can inform the effects of data filtration on the
alignment of \code{logFC} profiles. Also see \code{\link{standPep}} and
\code{\link{standPrn}} for more examples.}

\item{show_vline}{Logical; if TRUE, shows the vertical lines at \code{x = 0}.
The default is TRUE.}

\item{scale_y}{Logical; if TRUE, scale data on the \code{y-axis}. The default
is TRUE.}

\item{filepath}{A file path to output results. By default, it will be
determined automatically by the name of the calling function and the value
of \code{id} in the \code{call}.}

\item{filename}{A representative file name to outputs. By default, the
name(s) will be determined automatically. For text files, a typical file
extension is \code{.txt}. For image files, they are typically saved via
\code{\link[ggplot2]{ggsave}} or \code{\link[pheatmap]{pheatmap}} where the
image type will be determined by the extension of the file name.}

\item{theme}{A
\href{https://ggplot2.tidyverse.org/reference/ggtheme.html}{ggplot2} theme,
i.e., theme_bw(), or a custom theme. At the NULL default, a system theme
will be applied.}

\item{...}{\code{filter_}: Variable argument statements for the row
filtration of data against the column keys in \code{Peptide.txt} for
peptides or \code{Protein.txt} for proteins. Each statement contains to a
list of logical expression(s). The \code{lhs} needs to start with
\code{filter_}. The logical condition(s) at the \code{rhs} needs to be
enclosed in \code{exprs} with round parenthesis. \cr \cr For example,
\code{pep_len} is a column key in \code{Peptide.txt}. The statement
\code{filter_peps_at = exprs(pep_len <= 50)} will remove peptide entries
with \code{pep_len > 50}. See also \code{\link{normPSM}}. \cr \cr
Additional parameters for plotting with \code{ggplot2}: \cr \code{xmin},
the minimum \eqn{x} at a log2 scale; the default is -2. \cr \code{xmax},
the maximum \eqn{x} at a log2 scale; the default is +2. \cr \code{xbreaks},
the breaks in \eqn{x}-axis at a log2 scale; the default is 1. \cr
\code{binwidth}, the binwidth of \code{log2FC}; the default is \eqn{(xmax -
xmin)/80}. \cr \code{ncol}, the number of columns; the default is 1. \cr
\code{width}, the width of plot; \cr \code{height}, the height of plot. \cr
\code{scales}, should the scales be fixed across panels; the default is
"fixed" and the alternative is "free".}
}
\description{
Plots histograms
}
