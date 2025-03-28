% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histo.R
\name{pepHist}
\alias{pepHist}
\alias{prnHist}
\title{Histogram visualization}
\usage{
pepHist(
  col_select = NULL,
  scale_log2r = TRUE,
  complete_cases = FALSE,
  cut_points = c(mean_lint = NA),
  show_curves = TRUE,
  show_vline = TRUE,
  scale_y = TRUE,
  df = NULL,
  filepath = NULL,
  filename = NULL,
  theme = NULL,
  ...
)

prnHist(
  col_select = NULL,
  scale_log2r = TRUE,
  complete_cases = FALSE,
  cut_points = c(mean_lint = NA),
  show_curves = TRUE,
  show_vline = TRUE,
  scale_y = TRUE,
  df = NULL,
  filepath = NULL,
  filename = NULL,
  theme = NULL,
  ...
)
}
\arguments{
\item{col_select}{Character string to a column key in \code{expt_smry.xlsx}.
At the \code{NULL} default, the column key of \code{Select} in
\code{expt_smry.xlsx} will be used. In the case of no samples being
specified under \code{Select}, the column key of \code{Sample_ID} will be
used. The non-empty entries under the ascribing column will be used in
indicated analysis.}

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

\item{df}{The name of a primary data file. By default, it will be determined
automatically after matching the types of data and analysis with an
\code{id} among \code{c("pep_seq", "pep_seq_mod", "prot_acc", "gene")}. A
primary file contains normalized peptide or protein data and is among
\code{c("Peptide.txt", "Peptide_pVal.txt", "Peptide_impNA_pVal.txt",
"Protein.txt", "Protein_pVal.txt", "protein_impNA_pVal.txt")}. For analyses
require the fields of significance p-values, the \code{df} will be one of
\code{c("Peptide_pVal.txt", "Peptide_impNA_pVal.txt", "Protein_pVal.txt",
"protein_impNA_pVal.txt")}.}

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
\value{
Histograms of \code{log2FC}; raw histogram data:
  \code{[...]_raw.txt}; fitted data for curves: \code{[...]_fitted.txt}
}
\description{
\code{pepHist} plots the histograms of peptide \code{log2FC}.

\code{prnHist} plots the histograms of protein \code{log2FC}.
}
\details{
In the histograms, the \code{log2FC} under each TMT channel are color-coded by
their contributing reporter-ion or LFQ intensity.
}
\examples{
\donttest{
# ===================================
# Histogram
# ===================================

## !!!require the brief working example in `?load_expts`

## examplary `MGKernel` alignment
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
)

standPrn(
  method_align = MGKernel, 
  n_comp = 2, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
)

## (1) effects of data scaling
# peptide without log2FC scaling
pepHist(scale_log2r = FALSE)

# with scaling
pepHist(scale_log2r = TRUE)

## (2) sample column selection
# sample IDs indicated under column `Select` in `expt_smry.xlsx`
pepHist(col_select = Select, filename = colsel.png)

# protein data for samples under column `W2` in `expt_smry.xlsx`
prnHist(col_select = W2, filename = w2.png)

## (3) row filtration of data
# exclude oxidized methione or deamidated asparagine
pepHist(
  # filter_by = exprs(!grepl("[mn]", pep_seq_mod)),
  filter_by = exprs(not_contain_chars_in("mn", pep_seq_mod)),
  filename = "no_mn.png",
)

# phosphopeptide subset (error message if no matches)
pepHist(
  filter_peps = exprs(contain_chars_in("sty", pep_seq_mod)), 
  scale_y = FALSE, 
  filename = phospho.png,
)

# or use `grepl` directly
pepHist(
  filter_by = exprs(grepl("[sty]", pep_seq_mod)),
  filename = same_phospho.png,
)

## (4) between lead and lag
# leading profiles
pepHist(
  filename = lead.png,
)

# lagging profiles at
#   (1) n_psm >= 10
#   (2) and no methionine oxidation or asparagine deamidation
pepHist(
  filter_peps_by_npsm = exprs(pep_n_psm >= 10),
  filter_peps_by_mn = exprs(not_contain_chars_in("mn", pep_seq_mod)),
  filename = lag.png,
)

## (5) Data binning by `prot_icover`
pepHist(
  cut_points = c(prot_icover = NA),
  filename = prot_icover_coded.png,
)

## (6) custom theme
library(ggplot2)
my_histo_theme <- theme_bw() + theme(
  axis.text.x  = element_text(angle=0, vjust=0.5, size=18),
  axis.ticks.x  = element_blank(), # x-axis ticks
  axis.text.y  = element_text(angle=0, vjust=0.5, size=18),
  axis.title.x = element_text(colour="black", size=24),
  axis.title.y = element_text(colour="black", size=24),
  plot.title = element_text(colour="black", size=24, hjust=.5, vjust=.5),
  
  strip.text.x = element_text(size = 18, colour = "black", angle = 0),
  strip.text.y = element_text(size = 18, colour = "black", angle = 90),
  
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  
  legend.key = element_rect(colour = NA, fill = 'transparent'),
  legend.background = element_rect(colour = NA,  fill = "transparent"),
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=18),
  legend.text.align = 0,
  legend.box = NULL
)

pepHist(
  theme = my_histo_theme,
  filename = my_theme.png,
)

pepHist(
  col_select = BI_1,
  theme = theme_dark(),
  filename = bi1_dark.png,
)


## (7) direct uses of ggplot2
library(ggplot2)
res <- pepHist(filename = default.png)

# names(res)

p <- ggplot() +
  geom_histogram(data = res$raw, aes(x = value, y = ..count.., fill = Int_index),
                 color = "white", alpha = .8, binwidth = .05, size = .1) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  labs(title = "", x = expression("Ratio (" * log[2] * ")"), y = expression("Frequency")) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1),
                     labels = as.character(seq(-2, 2, by = 1))) +
  scale_y_continuous(limits = NULL) + 
  facet_wrap(~ Sample_ID, ncol = 5, scales = "fixed") # + 
  # my_histo_theme

p <- p + 
  geom_line(data = res$fitted, mapping = aes(x = x, y = value, colour = variable), size = .2) +
  scale_colour_manual(values = c("gray", "gray", "gray", "black"), name = "Gaussian",
                      breaks = c(c("G1", "G2", "G3"), paste(c("G1", "G2", "G3"), collapse = " + ")),
                      labels = c("G1", "G2", "G3", "G1 + G2 + G3"))

p <- p + geom_vline(xintercept = 0, size = .25, linetype = "dashed")

ggsave(file.path(dat_dir, "Peptide/Histogram/my_ggplot2.png"), 
       width = 22, height = 48, limitsize = FALSE)

\dontrun{
# sample selection
pepHist(
  col_select = "a_column_key_not_in_`expt_smry.xlsx`",
)

# data filtration
pepHist(
  filter_by = exprs(!grepl("[m]", a_column_key_not_in_data_table)),
)

prnHist(
  lhs_not_start_with_filter_ = exprs(n_psm >= 5),
)  
}
}
}
\seealso{
\emph{Metadata} \cr 
 \code{\link{load_expts}} for metadata preparation and a reduced working example in data normalization \cr

 \emph{Data normalization} \cr 
 \code{\link{normPSM}} for extended examples in PSM data normalization \cr
 \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
 \code{\link{mergePep}} for extended examples in peptide data merging \cr 
 \code{\link{standPep}} for extended examples in peptide data normalization \cr
 \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
 \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
 \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
 \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
 \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
 
 \emph{Variable arguments of `filter_...`} \cr 
 \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
 \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
 \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
 \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
 
 \emph{Missing values} \cr 
 \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
 
 \emph{Informatics} \cr 
 \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
 \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
 \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
 \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
 \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
 \code{\link{prnGSVA}} for gene set variance analysis \cr 
 \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
 \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
 \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA visualization \cr 
 \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr 
 \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
 \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
 \code{\link{prnCorr_logInt}}  for correlation plots \cr 
 \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
 \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
 \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
 \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
 
 \emph{Custom databases} \cr 
 \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
 \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
 \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
 ontology}} \cr 
 \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
 signatures} \cr 
 \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
 
 \emph{Column keys in PSM, peptide and protein outputs} \cr 
 system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
 system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
 system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
}
