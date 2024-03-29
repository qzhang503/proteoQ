% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcanos.R
\name{fullVolcano}
\alias{fullVolcano}
\title{Volcano plots for all proteins or peptides in a data set}
\usage{
fullVolcano(
  df = NULL,
  id = "gene",
  contrast_groups = NULL,
  theme = NULL,
  fml_nm = NULL,
  filepath = NULL,
  filename = NULL,
  adjP = FALSE,
  topn_labels = 20,
  highlights = NULL,
  grids = NULL,
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

\item{contrast_groups}{The contrast groups defined under a formula at \code{fml_nm}.}

\item{theme}{A
\href{https://ggplot2.tidyverse.org/reference/ggtheme.html}{ggplot2} theme,
i.e., theme_bw(), or a custom theme. At the NULL default, a system theme
will be applied.}

\item{fml_nm}{A character string; the name of \code{fml}.}

\item{filepath}{A file path to output results. By default, it will be
determined automatically by the name of the calling function and the value
of \code{id} in the \code{call}.}

\item{filename}{A representative file name to outputs. By default, the
name(s) will be determined automatically. For text files, a typical file
extension is \code{.txt}. For image files, they are typically saved via
\code{\link[ggplot2]{ggsave}} or \code{\link[pheatmap]{pheatmap}} where the
image type will be determined by the extension of the file name.}

\item{adjP}{Logical; if TRUE, use Benjamini-Hochberg pVals in volcano plots.
The default is FALSE.}

\item{topn_labels}{A non-negative integer; the top-n species for labeling in a
plot. At \code{topn_labels = 0}, no labels of proteins/peptides will be
shown. The default is to label the top-20 species with the lowest p-values.}

\item{highlights}{A list of entries for highlighting. See also \code{filter_}
varargs for the format.}

\item{grids}{An integer or integer vector for subset visualization of
contrasts within a group. For example with a group of three contrasts,
\code{grids = 2:3} will hide the first grid from displaying. At the
NULL default, all available grids will be shown.}

\item{...}{\code{filter_}: Variable argument statements for the row filtration
against data in a primary file linked to \code{df}. See also
\code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
Additional parameters for plotting: \cr \code{xco}, the cut-off lines of
fold changes at position \code{x}; the default is at \eqn{-1.2} and
\eqn{+1.2}. \cr \code{yco}, the cut-off line of \code{pVal} at position
\code{y}; the default is \eqn{0.05}. \cr \code{width}, the width of plot;
\cr \code{height}, the height of plot. \cr \code{nrow}, the number of rows
in a plot. \cr \code{xmin}, the minimum \code{x}. \cr \code{xmax}, the
maximum \code{x}. \cr \code{ymin}, the minimum \code{y}. \cr \code{ymax},
the maximum \code{y}. \cr \code{x_label}, the label on \code{x}. \cr
\code{y_label}, the label on \code{y}.}
}
\description{
Volcano plots for all proteins or peptides in a data set
}
