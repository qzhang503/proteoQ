% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mds.R
\name{plotMDS}
\alias{plotMDS}
\title{Plots MDS}
\usage{
plotMDS(
  df = NULL,
  id = NULL,
  label_scheme_sub = NULL,
  choice = "cmdscale",
  dist_co = log2(1),
  adjEucDist = FALSE,
  method = "euclidean",
  p = 2,
  k = 3,
  dimension = 2,
  folds = 1,
  show_ids = FALSE,
  show_ellipses = FALSE,
  col_group = NULL,
  col_color = NULL,
  col_fill = NULL,
  col_shape = NULL,
  col_size = NULL,
  col_alpha = NULL,
  color_brewer = NULL,
  fill_brewer = NULL,
  size_manual = NULL,
  shape_manual = NULL,
  alpha_manual = NULL,
  scale_log2r = TRUE,
  complete_cases = FALSE,
  filepath = NULL,
  filename = NULL,
  center_features = TRUE,
  scale_features = TRUE,
  theme = NULL,
  anal_type = "MDS",
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

\item{choice}{Character string; the MDS method in \code{c("cmdscale",
"isoMDS")}. The default is "cmdscale".}

\item{dist_co}{Numeric; The cut-off in the absolute distance measured by
\eqn{d = abs(x_i - x_j)}. Data pairs, \eqn{x_i} and \eqn{x_j}, with
corresponding \eqn{d} smaller than \code{dist_co} will be excluded from
distance calculations by \link[stats]{dist}. The default is no distance
cut-off at \eqn{dist_co = log2(1)}.}

\item{adjEucDist}{Logical; if TRUE, adjusts the inter-plex \code{Euclidean}
distance by \eqn{1/sqrt(2)} at \code{method = "euclidean"}. The option
\code{adjEucDist = TRUE} may be suitable when \code{reference samples} from
each TMT plex undergo approximately the same sample handling process as the
samples of interest. For instance, \code{reference samples} were split at
the levels of protein lysates. Typically, \code{adjEucDist = FALSE} if
\code{reference samples} were split near the end of a sample handling
process, for instance, at the stages immediately before or after TMT
labeling. Also see online
\href{https://github.com/qzhang503/proteoQ}{README, section MDS} for a brief
reasoning.}

\item{method}{Character string; the distance measure in one of c("euclidean",
"maximum", "manhattan", "canberra", "binary") for \code{\link[stats]{dist}}.
The default method is "euclidean".}

\item{p}{Numeric; The power of the Minkowski distance in
\code{\link[stats]{dist}}. The default is 2.}

\item{k}{Numeric; The desired dimension for the solution passed to
\code{\link[stats]{cmdscale}}. The default is 3.}

\item{dimension}{Numeric; The desired dimension for pairwise visualization.
The default is 2.}

\item{folds}{Not currently used. Integer; the degree of folding data into
subsets. The default is one without data folding.}

\item{show_ids}{Logical; if TRUE, shows the sample IDs in \code{MDS/PCA}
plots. The default is TRUE.}

\item{show_ellipses}{Logical; if TRUE, shows the ellipses by sample groups
according to \code{col_group}. The default is FALSE.}

\item{col_group}{Character string to a column key in \code{expt_smry.xlsx}.
Samples corresponding to non-empty entries under \code{col_group} will be
used for sample grouping in the indicated analysis. At the NULL default, the
column key \code{Group} will be used. No data annotation by groups will be
performed if the fields under the indicated group column is empty.}

\item{col_color}{Character string to a column key in \code{expt_smry.xlsx}.
Values under which will be used for the \code{color} aesthetics in plots. At
the NULL default, the column key \code{Color} will be used. If NA, bypasses
the aesthetics (a means to bypass the look-up of column \code{Color} and
handle duplication in aesthetics).}

\item{col_fill}{Character string to a column key in \code{expt_smry.xlsx}.
Values under which will be used for the \code{fill} aesthetics in plots. At
the NULL default, the column key \code{Fill} will be used. If NA, bypasses
the aesthetics (a means to bypass the look-up of column \code{Fill} and
handle duplication in aesthetics).}

\item{col_shape}{Character string to a column key in \code{expt_smry.xlsx}.
Values under which will be used for the \code{shape} aesthetics in plots. At
the NULL default, the column key \code{Shape} will be used. If NA, bypasses
the aesthetics (a means to bypass the look-up of column \code{Shape} and
handle duplication in aesthetics).}

\item{col_size}{Character string to a column key in \code{expt_smry.xlsx}.
Values under which will be used for the \code{size} aesthetics in plots. At
the NULL default, the column key \code{Size} will be used. If NA, bypasses
the aesthetics (a means to bypass the look-up of column \code{Size} and
handle duplication in aesthetics).}

\item{col_alpha}{Character string to a column key in \code{expt_smry.xlsx}.
Values under which will be used for the \code{alpha} (transparency)
aesthetics in plots. At the NULL default, the column key \code{Alpha} will
be used. If NA, bypasses the aesthetics (a means to bypass the look-up of
column \code{Alpha} and handle duplication in aesthetics).}

\item{color_brewer}{Character string to the name of a color brewer for use in
\href{https://ggplot2.tidyverse.org/reference/scale_brewer.html}{ggplot2::scale_color_brewer},
 i.e., \code{color_brewer = Set1}. At the NULL default, the setting in
\code{ggplot2} will be used.}

\item{fill_brewer}{Character string to the name of a color brewer for use in
\href{https://ggplot2.tidyverse.org/reference/scale_brewer.html}{ggplot2::scale_fill_brewer},
 i.e., \code{fill_brewer = Spectral}. At the NULL default, the setting in
\code{ggplot2} will be used.}

\item{size_manual}{Numeric vector to the scale of sizes for use in
\href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_size_manual},
 i.e., \code{size_manual = c(8, 12)}. At the NULL default, the setting in
\code{ggplot2} will be used.}

\item{shape_manual}{Numeric vector to the scale of shape IDs for use in
\href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_shape_manual},
 i.e., \code{shape_manual = c(5, 15)}. At the NULL default, the setting in
\code{ggplot2} will be used.}

\item{alpha_manual}{Numeric vector to the scale of transparency of objects for
use in
\href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_alpha_manual}
 , i.e., \code{alpha_manual = c(.5, .9)}. At the NULL default, the setting
in \code{ggplot2} will be used.}

\item{scale_log2r}{Logical; if TRUE, adjusts \code{log2FC} to the same scale
of standard deviation across all samples. The default is TRUE. At
\code{scale_log2r = NA}, the raw \code{log2FC} without normalization will
be used.}

\item{complete_cases}{Logical; if TRUE, only cases that are complete with no
missing values will be used. The default is FALSE.}

\item{filepath}{A file path to output results. By default, it will be
determined automatically by the name of the calling function and the value
of \code{id} in the \code{call}.}

\item{filename}{A representative file name to outputs. By default, the
name(s) will be determined automatically. For text files, a typical file
extension is \code{.txt}. For image files, they are typically saved via
\code{\link[ggplot2]{ggsave}} or \code{\link[pheatmap]{pheatmap}} where the
image type will be determined by the extension of the file name.}

\item{center_features}{Logical; if TRUE, adjusts log2FC to center zero by
features (proteins or peptides). The default is TRUE. Note the difference to
data alignment with \code{method_align} in \code{\link{standPrn}}
or \code{\link{standPep}} where log2FC are aligned by observations
(samples).}

\item{scale_features}{Logical; if TRUE, adjusts log2FC to the same scale of
variance by features (protein or peptide entries). The default is TRUE. Note
the difference to data scaling with \code{scale_log2r} where log2FC are
scaled by observations (samples).}

\item{theme}{A
\href{https://ggplot2.tidyverse.org/reference/ggtheme.html}{ggplot2} theme,
i.e., theme_bw(), or a custom theme. At the NULL default, a system theme
will be applied.}

\item{anal_type}{Character string; the type of analysis that are preset for
method dispatch in function factories. The value will be determined
automatically. Exemplary values include \code{anal_type = c("PCA",
"Corrplot", "EucDist", "GSPA", "Heatmap", "Histogram", "MDS", "Model",
"NMF", "Purge", "Trend", "LDA", ...)}.}

\item{...}{\code{filter_}: Variable argument statements for the row filtration
against data in a primary file linked to \code{df}. See also
\code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
Additional parameters for \code{ggsave}: \cr \code{width}, the width of
plot; \cr \code{height}, the height of plot \cr \code{...}}
}
\description{
Plots MDS
}
