#'Histogram visualization
#'
#'\code{proteoHist} plots the histograms of protein or peptide \code{log2FC}.
#'Users should avoid calling the method directly, but instead use the following
#'wrappers.
#'
#'In the histograms, the \code{log2FC} under each TMT channel are color-coded by
#'their contributing reporter-ion intensity.
#'
#'@param id Character string to indicate the type of data. The value will be
#'  determined automatically. Peptide data will be used at \code{id = pep_seq}
#'  or \code{pep_seq_mod}, and protein data will be used at \code{id = prot_acc}
#'  or \code{gene}.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2FC} to the same scale
#'  of standard deviation across all samples. The default is TRUE.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used. The default is FALSE.
#'@param show_curves Logical; if TRUE, shows the fitted curves. At the TRUE
#'  default, the curve parameters are based on the latest call to
#'  \code{\link{standPep}} or \code{\link{standPrn}} with \code{method_align =
#'  MGKernel}. This feature can inform the effects of data filtration on the
#'  alignment of \code{logFC} profiles. Also see \code{\link{standPep}} and
#'  \code{\link{standPrn}} for more examples.
#'@param show_vline Logical; if TRUE, shows the vertical lines at \code{x = 0}.
#'  The default is TRUE.
#'@param scale_y Logical; if TRUE, scale data on the \code{y-axis}. The default
#'  is TRUE.
#'@param df The name of input data file. By default, it will be determined
#'  automatically after matching the data type with an \code{id} among
#'  \code{c("pep_seq", "pep_seq_mod", "prot_acc", "gene")}.
#'@param filepath A file path to output results. By default, it will be
#'  determined automatically by the name of the calling function and the value
#'  of \code{id} in the \code{call}.
#'@param filename A representative file name to outputs. By default, the name(s)
#'  will be determined automatically. For text files, a typical file extension
#'  is \code{.txt}. For image files, they are typically saved via
#'  \code{\link[ggplot2]{ggsave}} or \code{\link[pheatmap]{pheatmap}} where the
#'  image type will be determined by the extension of the file name.
#'@param theme A
#'  \code{\href{https://ggplot2.tidyverse.org/reference/ggtheme.html}{ggplot2}}
#'  theme, i.e., theme_bw(), or a custom theme. At the NULL default, a system
#'  theme will be applied.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  of data against the column keys in \code{Peptide.txt} for peptides or
#'  \code{Protein.txt} for proteins. Each statement contains to a list of
#'  logical expression(s). The \code{lhs} needs to start with \code{filter_}.
#'  The logical condition(s) at the \code{rhs} needs to be enclosed in
#'  \code{exprs} with round parenthesis. \cr \cr For example, \code{pep_len} is
#'  a column key present in \code{Mascot} peptide tables of \code{Peptide.txt}.
#'  The statement \code{filter_peps_at = exprs(pep_len <= 50)} will remove
#'  peptide entries with \code{pep_len > 50}. See also \code{\link{normPSM}}.
#'  \cr \cr Additional parameters for plotting with \code{ggplot2}: \cr
#'  \code{xmin}, the minimum \eqn{x} at a log2 scale; the default is -2. \cr
#'  \code{xmax}, the maximum \eqn{x} at a log2 scale; the default is +2. \cr
#'  \code{xbreaks}, the breaks in \eqn{x}-axis at a log2 scale; the default is
#'  1. \cr \code{binwidth}, the binwidth of \code{log2FC}; the default is
#'  \eqn{(xmax - xmin)/80}. \cr \code{ncol}, the number of columns; the default
#'  is 1. \cr \code{width}, the width of plot; \cr \code{height}, the height of
#'  plot. \cr \code{scales}, should the scales be fixed across panels; the
#'  default is "fixed" and the alternative is "free".
#'@inheritParams standPep
#'@seealso \code{\link{load_expts}} for a reduced working example in data
#'  normalization \cr
#'
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization
#'  \cr \code{\link{mergePep}} for extended examples in peptide data merging \cr
#'  \code{\link{standPep}} for extended examples in peptide data normalization
#'  \cr \code{\link{Pep2Prn}} for extended examples in peptide to protein
#'  summarization \cr \code{\link{standPrn}} for extended examples in protein
#'  data normalization. \cr \code{\link{purgePSM}} and \code{\link{purgePep}}
#'  for extended examples in data purging \cr \code{\link{pepHist}} and
#'  \code{\link{prnHist}} for extended examples in histogram visualization. \cr
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting
#'  MS file names \cr
#'
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr
#'
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation
#'  \cr \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization
#'  \cr
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnHist_.R
#'
#'@return Histograms of \code{log2FC}
#'@export
proteoHist <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                        col_select = NULL, 
                        scale_log2r = TRUE, complete_cases = FALSE, 
                        show_curves = TRUE, show_vline = TRUE, scale_y = TRUE, 
                        df = NULL, filepath = NULL, filename = NULL, theme = NULL, ...) {
  
  old_opt <- options(max.print = 99999, warn = 0)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 2000000, warn = 1)
  
  id <- rlang::enexpr(id)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
  
  stopifnot(rlang::is_logical(show_curves), 
            rlang::is_logical(show_vline), 
            rlang::is_logical(scale_y))
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  dots <- rlang::enexprs(...)
  if (!is.null(dots$impute_na)) {
    dots$impute_na <- NULL
    rlang::warn("No NA imputation with histograms.")
  }
  
  info_anal(id = !!id, col_select = !!col_select, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = FALSE,
            df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "Histogram")(show_curves = show_curves,
                                     show_vline = show_vline, scale_y = scale_y, 
                                     theme = theme, !!!dots)
}


prnHist_wrap <- function (...) {
  check_dots(c("id", "anal_type"), ...)
  id <- match_call_arg(normPSM, group_pep_by)
  proteoHist(id = !!id, ...)
}



#'Visualization of MDS plots
#'
#'\code{proteoMDS} visualizes the results from multidimensional scaling (MDS).
#'Users should avoid calling the method directly, but instead use the following
#'wrappers.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}}, followed by a metric
#'(\code{\link[stats]{cmdscale}}) or non-metric (\code{\link[MASS]{isoMDS}})
#'MDS. The default is metric MDS with the input dissimilarities being euclidean
#'distances.
#'
#'@inheritParams  proteoHist
#'@inheritParams proteoHM
#'@param  col_group Not used.
#'@param  col_color Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{color} aesthetics in plots. At
#'  the NULL default, the column key \code{Color} will be used.
#'@param  col_fill Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{fill} aesthetics in plots. At
#'  the NULL default, the column key \code{Fill} will be used.
#'@param  col_shape Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{shape} aesthetics in plots. At
#'  the NULL default, the column key \code{Shape} will be used.
#'@param  col_size Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{size} aesthetics in plots. At
#'  the NULL default, the column key \code{Size} will be used.
#'@param  col_alpha Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{alpha} (transparency)
#'  aesthetics in plots. At the NULL default, the column key \code{Alpha} will
#'  be used.
#'@param color_brewer Character string to the name of a color brewer for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_brewer.html}{ggplot2::scale_color_brewer},
#'   i.e., \code{color_brewer = Set1}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param fill_brewer Character string to the name of a color brewer for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_brewer.html}{ggplot2::scale_fill_brewer},
#'   i.e., \code{fill_brewer = Spectral}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param size_manual Numeric vector to the scale of sizes for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_size_manual},
#'   i.e., \code{size_manual = c(8, 12)}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param shape_manual Numeric vector to the scale of shape IDs for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_shape_manual},
#'   i.e., \code{shape_manual = c(5, 15)}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param alpha_manual Numeric vector to the scale of transparency of objects for
#'  use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_alpha_manual}
#'   , i.e., \code{alpha_manual = c(.5, .9)}. At the NULL default, the setting
#'  in \code{ggplot2} will be used.
#'@param adjEucDist Logical; if TRUE, adjusts the inter-plex Euclidean distance
#'  by \eqn{1/sqrt(2)}. The option \code{adjEucDist = TRUE} may be suitable when
#'  \code{reference samples} from each TMT plex undergo approximately the same
#'  sample handling process as the samples of interest. For instance,
#'  \code{reference samples} were split at the levels of protein lysates.
#'  Typically, \code{adjEucDist = FALSE} if \code{reference samples} were split
#'  near the end of a sample handling process, for instance, at the stages
#'  immediately before or after TMT labeling.
#'@param classical Logical. Metric MDS will be performed at TRUE and non-metric
#'  MDS at FALSE (see also \code{\link[stats]{cmdscale}} and
#'  \code{\link[MASS]{isoMDS}}). The default is TRUE.
#'@param method Character string; the distance measure in
#'  \code{\link[stats]{dist}}. The default method is "euclidean".
#'@param p Numeric; The power of the Minkowski distance in
#'  \code{\link[stats]{dist}}. The default is 2.
#'@param k Numeric; The desired dimension for the solution passed to
#'  \code{\link[stats]{cmdscale}}. The default is 3.
#'@param show_ids Logical; if TRUE, shows the sample IDs in \code{MDS/PCA}
#'  plots. The default is TRUE.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  peptide or protein data; also see \code{\link{normPSM}}. \cr
#'  \code{arrange_}: Logical expression(s) for the row order of data; also see
#'  \code{\link{prnHM}}. \cr \cr Additional parameters for \code{ggsave}: \cr
#'  \code{width}, the width of plot; \cr \code{height}, the height of plot \cr
#'  \code{...}
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data
#'  normalization \cr \code{\link{normPSM}} for extended examples in PSM data
#'  normalization \cr \code{\link{PSM2Pep}} for extended examples in PSM to
#'  peptide summarization \cr \code{\link{mergePep}} for extended examples in
#'  peptide data merging \cr \code{\link{standPep}} for extended examples in
#'  peptide data normalization \cr \code{\link{Pep2Prn}} for extended examples
#'  in peptide to protein summarization \cr \code{\link{standPrn}} for extended
#'  examples in protein data normalization. \cr \code{\link{pepHist}} and
#'  \code{\link{prnHist}} for extended examples in histogram visualization. \cr
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in
#'  data purging \cr \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr \code{\link{pepImp}} and
#'  \code{\link{prnImp}} for missing value imputation \cr \code{\link{pepSig}}
#'  and \code{\link{prnSig}} for significance tests \cr \code{\link{pepHM}} and
#'  \code{\link{prnHM}} for heat map visualization \cr \code{\link{pepMDS}} and
#'  \code{\link{prnMDS}} for MDS visualization \cr \code{\link{pepPCA}} and
#'  \code{\link{prnPcA}} for PCA visualization \cr
#'@example inst/extdata/examples/prnMDS_.R
#'
#'@return MDS plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoMDS <- function (id = gene,
                       col_select = NULL, col_color = NULL, col_fill = NULL,
                       col_shape = NULL, col_size = NULL, col_alpha = NULL,
                       color_brewer = NULL, fill_brewer = NULL, 
                       size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                       scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                       adjEucDist = FALSE, classical = TRUE, 
                       method = "euclidean", p = 2, k = 3, 
                       show_ids = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                       theme = NULL, ...) {
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- rlang::enexpr(id)
  col_select <- rlang::enexpr(col_select)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)
  
  color_brewer <- rlang::enexpr(color_brewer)
  fill_brewer <- rlang::enexpr(fill_brewer)
  size_manual <- rlang::enexpr(size_manual)
  shape_manual <- rlang::enexpr(shape_manual)
  alpha_manual <- rlang::enexpr(alpha_manual)
  
  method <- rlang::as_string(rlang::enexpr(method))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, method = method, 
                               p = p, k = k, show_ids = show_ids,
                               theme = theme, ...)
}



#'Visualization of PCA plots
#'
#'\code{proteoPCA} visualizes the results from principal component analysis
#'(PCA). Users should avoid calling the method directly, but instead use the
#'following wrappers.
#'
#'\code{log2FC} are used in PCA (\code{\link[stats]{prcomp}}).
#'
#'@inheritParams prnHist
#'@inheritParams proteoHM
#'@inheritParams prnMDS
#'@param complete_cases Logical; always TRUE for PCA. 
#'@param type Character string indicating the type of PCA. At the \code{type =
#'  obs} default, the components are by observations; at \code{type = feats},
#'  the components are by features.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr \code{arrange_}: Logical
#'  expression(s) for the row order of data; also see \code{\link{prnHM}}. \cr
#'  \cr Additional parameters for \code{ggsave}: \cr \code{width}, the width of
#'  plot; \cr \code{height}, the height of plot \cr \code{...}
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data
#'  normalization \cr \code{\link{normPSM}} for extended examples in PSM data
#'  normalization \cr \code{\link{PSM2Pep}} for extended examples in PSM to
#'  peptide summarization \cr \code{\link{mergePep}} for extended examples in
#'  peptide data merging \cr \code{\link{standPep}} for extended examples in
#'  peptide data normalization \cr \code{\link{Pep2Prn}} for extended examples
#'  in peptide to protein summarization \cr \code{\link{standPrn}} for extended
#'  examples in protein data normalization. \cr \code{\link{pepHist}} and
#'  \code{\link{prnHist}} for extended examples in histogram visualization. \cr
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in
#'  data purging \cr \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr \code{\link{pepImp}} and
#'  \code{\link{prnImp}} for missing value imputation \cr \code{\link{pepSig}}
#'  and \code{\link{prnSig}} for significance tests \cr \code{\link{pepHM}} and
#'  \code{\link{prnHM}} for heat map visualization \cr \code{\link{pepMDS}} and
#'  \code{\link{prnMDS}} for MDS visualization \cr
#'@example inst/extdata/examples/prnPCA_.R
#'
#'@return PCA plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoPCA <- function (id = gene, type = "obs", 
                       col_select = NULL, col_color = NULL, 
                       col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                       color_brewer = NULL, fill_brewer = NULL, 
                       size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                       scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                       show_ids = TRUE, 
                       df = NULL, filepath = NULL, filename = NULL, 
                       theme = NULL, ...) {
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- rlang::enexpr(id)
  col_select <- rlang::enexpr(col_select)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)
  
  color_brewer <- rlang::enexpr(color_brewer)
  fill_brewer <- rlang::enexpr(fill_brewer)
  size_manual <- rlang::enexpr(size_manual)
  shape_manual <- rlang::enexpr(shape_manual)
  alpha_manual <- rlang::enexpr(alpha_manual)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  type <- rlang::as_string(rlang::enexpr(type))
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "PCA")(type = type, show_ids = show_ids, theme = theme, ...)
}


#'Visualization of the Euclidean distance matrix
#'
#'\code{proteoEucDist} visualizes the heat map of Euclidean distances. Users
#'should avoid calling the method directly, but instead use the following
#'wrappers.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}} for heat map visualization.
#'
#'@inheritParams prnHist
#'@inheritParams prnMDS
#'@inheritParams proteoHM
#'@param annot_cols A character vector of column keys in \code{expt_smry.xlsx}.
#'  The values under the selected keys will be used to color-code sample IDs on
#'  the top of the indicated plot. The default is NULL without column
#'  annotation.
#'@param annot_colnames A character vector of replacement name(s) to
#'  \code{annot_cols}. The default is NULL without name replacement.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  of data against the column keys in \code{Peptide.txt}, \code{Protein.txt}
#'  etc. see also \code{\link{normPSM}} \cr \cr \code{arrange_}: Logical
#'  expression(s) for the row ordering of data; see also \code{\link{prnHM}}.
#'  \cr \cr Additional parameters for plotting: \cr \code{width}, the width of
#'  plot \cr \code{height}, the height of plot \cr \cr Additional arguments for
#'  \code{\link[pheatmap]{pheatmap}} \cr Note
#'  arguments disabled from \code{pheatmap}: \cr \code{annotation_col}; instead
#'  use keys indicated in \code{annot_cols} \cr \code{annotation_row}; instead
#'  use keys indicated in \code{annot_rows}
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data
#'  normalization \cr \code{\link{normPSM}} for extended examples in PSM data
#'  normalization \cr \code{\link{PSM2Pep}} for extended examples in PSM to
#'  peptide summarization \cr \code{\link{mergePep}} for extended examples in
#'  peptide data merging \cr \code{\link{standPep}} for extended examples in
#'  peptide data normalization \cr \code{\link{Pep2Prn}} for extended examples
#'  in peptide to protein summarization \cr \code{\link{standPrn}} for extended
#'  examples in protein data normalization. \cr \code{\link{pepHist}} and
#'  \code{\link{prnHist}} for extended examples in histogram visualization. \cr
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in
#'  data purging \cr \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr \code{\link{pepImp}} and
#'  \code{\link{prnImp}} for missing value imputation \cr \code{\link{pepSig}}
#'  and \code{\link{prnSig}} for significance tests \cr \code{\link{pepHM}} and
#'  \code{\link{prnHM}} for heat map visualization \cr \code{\link{pepMDS}} and
#'  \code{\link{prnMDS}} for MDS visualization \cr \code{\link{pepPCA}} and
#'  \code{\link{prnPcA}} for PCA visualization \cr
#'@example inst/extdata/examples/prnEucDist_.R
#'@return Heat map visualization of distance matrices.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoEucDist <- function (id = gene, col_select = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                           adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, 
                           df = NULL, filepath = NULL, filename = NULL, ...) {
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- rlang::enexpr(id)
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = NULL, col_fill = NULL, 
            col_shape = NULL, col_size = NULL, col_alpha = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "EucDist")(adjEucDist = adjEucDist, 
                                   annot_cols = annot_cols, annot_colnames = annot_colnames, ...)
}



#'Visualization of heat maps
#'
#'\code{proteoHM} visualizes the heat maps of protein or peptide \code{log2FC}.
#'Users should avoid calling the method directly, but instead use the following
#'wrappers.
#'
#'Data rows without non-missing pairs will result in NA distances in inter-row
#'dissimilarities (\code{\link[stats]{dist}}). At \code{complet_cases = TRUE},
#'the data subset that are complete without missing values will be used. At
#'\code{impute_na = TRUE}, all data rows will be used with NA imputation (see
#'\code{\link{prnImp}}). At the default of \code{complet_cases = FALSE} and
#'\code{impute_na = FALSE}, NA distances will be arbitrarily replaced with the
#'mean value of the row-disance matrix for hierarchical row clustering.
#'
#'Similar to data rows, NA distances in data columns will be replaced with the
#'mean value of the column-distance matrix.
#'
#'To avoid memory failure, row aggregation using the \code{kmeans_k} option
#'(\code{\link[pheatmap]{pheatmap}}) may be considered for large data sets.
#'
#'
#'@inheritParams  prnEucDist
#'@param  col_benchmark Not used.
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used. The default is FALSE.
#'@param annot_rows A character vector of column keys that can be found from
#'  input files of \code{Peptide.txt}, \code{Protein.txt} et al. The values
#'  under the selected keys will be used to color-code peptides or proteins on
#'  the side of the indicated plot. The default is NULL without row annotation.
#'@param xmin  Numeric; the minimum x at a log2 scale. The default is -1.
#'@param xmax  Numeric; the maximum  x at a log2 scale. The default is 1.
#'@param xmargin  Numeric; the margin in heat scales. The default is 0.1.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  of data against the column keys in \code{Peptide.txt}, \code{Protein.txt}
#'  etc. Each statement contains to a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key present in \code{Mascot} peptide
#'  tables of \code{Peptide.txt}. The statement \code{filter_peps_at =
#'  exprs(pep_len <= 50)} will remove peptide entries with \code{pep_len > 50}.
#'  See also \code{\link{pepHist}}, \code{\link{normPSM}}. \cr \cr
#'  \code{arrange_}: Logical expression(s) for the row ordering of data. The
#'  \code{lhs} needs to start with \code{arrange_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{arrange_peps_by = exprs(gene, prot_n_pep)} will arrange
#'  entries by \code{gene}, then by \code{prot_n_pep}. \cr \cr Additional
#'  parameters for plotting: \cr \code{width}, the width of plot \cr
#'  \code{height}, the height of plot \cr \cr Additional arguments for
#'  \code{\link[pheatmap]{pheatmap}}, i.e., \code{cluster_rows}... \cr \cr Note
#'  arguments disabled for \code{pheatmap}: \cr \code{annotation_col}; instead
#'  use keys indicated in \code{annot_cols} \cr \code{annotation_row}; instead
#'  use keys indicated in \code{annot_rows}
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data
#'  normalization \cr
#'
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization
#'  \cr \code{\link{mergePep}} for extended examples in peptide data merging \cr
#'  \code{\link{standPep}} for extended examples in peptide data normalization
#'  \cr \code{\link{Pep2Prn}} for extended examples in peptide to protein
#'  summarization \cr \code{\link{standPrn}} for extended examples in protein
#'  data normalization. \cr \code{\link{purgePSM}} and \code{\link{purgePep}}
#'  for extended examples in data purging \cr \code{\link{pepHist}} and
#'  \code{\link{prnHist}} for extended examples in histogram visualization. \cr
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting
#'  MS file names \cr
#'
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr
#'
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation
#'  \cr \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization
#'  \cr
#'
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein
#'  significance pVals \cr \code{\link{gspaMap}} for mapping GSPA to volcano
#'  plot visualization \cr \code{\link{prnGSPAHM}} for heat map and network
#'  visualization of GSPA results \cr \code{\link{prnGSVA}} for gene set
#'  variance analysis \cr \code{\link{prnGSEA}} for data preparation for online
#'  GSEA. \cr
#'
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr
#'  \code{\link{pepPCA}} and \code{\link{prnPcA}} for PCA visualization \cr
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'  \code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for
#'  correlation plots \cr
#'
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for protein
#'  trend analysis and visualization \cr \code{\link{anal_pepNMF}},
#'  \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'  \code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for protein NMF
#'  analysis and visualization \cr
#'
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@example inst/extdata/examples/prnHM_.R
#'
#'@return Heat maps and optional sub trees.
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoHM <- function (id = gene, col_select = NULL, col_benchmark = NULL,
                      scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                      df = NULL, filepath = NULL, filename = NULL,
                      annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                      xmin = -1, xmax = 1, xmargin = 0.1, ...) {
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- rlang::enexpr(id)
  col_select <- rlang::enexpr(col_select)
  col_benchmark <- rlang::enexpr(col_benchmark)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, complete_cases = complete_cases,impute_na = impute_na, 
            df = !!df, filepath = !!filepath,
            filename = !!filename, anal_type = "Heatmap")(xmin = xmin, xmax = xmax, xmargin = xmargin,
                                                          annot_cols = annot_cols, 
                                                          annot_colnames = annot_colnames, 
                                                          annot_rows = annot_rows, ...)
}


#'Correlation Plots
#'
#'\code{proteoCorr} plots Pearson correlation for both \code{logFC} and
#'\code{intensity} data. Users should avoid call the method directly, but
#'instead use the following wrappers.
#'
#'The function matches the current \code{id} to the grouping argument in the
#'latest \code{call} to \code{\link{normPSM}}. See also \code{\link{prnHist}}
#'for details.
#'
#'@inheritParams proteoHist
#'@inheritParams proteoMDS
#'@param  col_order Character string to a column key in \code{expt_smry.xlsx}.
#'  Numeric values under which will be used for the left-to-right arrangement of
#'  samples in graphic outputs or top-to-bottom arrangement in text outputs. At
#'  the NULL default, the column key \code{Order} will be used. If values under
#'  column \code{Order} are left blank, samples will be ordered by their names.
#'@param data_select The subset of data to be selected. The value will be
#'  determined automatically. At default, \code{logFC} will be used; at
#'  \code{logInt}, intensity with \code{log10} transformation will be used.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr Additional parameters for
#'  plotting: \cr \code{width}, the width of plot \cr \code{height}, the height
#'  of plot \cr \code{xmin}, the minimum \eqn{x} of logFC or intensity \cr
#'  \code{xmax}, the maximum \eqn{x} of logFC data or intensity data \cr
#'  \code{xbreaks}, the breaks on \eqn{x} axis; the same breaks will be applied
#'  to \eqn{y} axis.
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data
#'  normalization \cr
#'
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization
#'  \cr \code{\link{mergePep}} for extended examples in peptide data merging \cr
#'  \code{\link{standPep}} for extended examples in peptide data normalization
#'  \cr \code{\link{Pep2Prn}} for extended examples in peptide to protein
#'  summarization \cr \code{\link{standPrn}} for extended examples in protein
#'  data normalization. \cr \code{\link{purgePSM}} and \code{\link{purgePep}}
#'  for extended examples in data purging \cr \code{\link{pepHist}} and
#'  \code{\link{prnHist}} for extended examples in histogram visualization. \cr
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting
#'  MS file names \cr
#'
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr
#'
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation
#'  \cr \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization
#'  \cr
#'
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein
#'  significance pVals \cr \code{\link{gspaMap}} for mapping GSPA to volcano
#'  plot visualization \cr \code{\link{prnGSPAHM}} for heat map and network
#'  visualization of GSPA results \cr \code{\link{prnGSVA}} for gene set
#'  variance analysis \cr \code{\link{prnGSEA}} for data preparation for online
#'  GSEA. \cr
#'
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr
#'  \code{\link{pepPCA}} and \code{\link{prnPcA}} for PCA visualization \cr
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'  \code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for
#'  correlation plots \cr
#'
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend
#'  analysis and visualization \cr \code{\link{anal_pepNMF}},
#'  \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'  \code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for NMF
#'  analysis and visualization \cr
#'
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@example inst/extdata/examples/prnCorr_.R
#'
#'@return Correlation plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoCorr <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                        data_select = c("logFC", "logInt"), 
                        col_select = NULL, col_order = NULL, 
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                        df = NULL, filepath = NULL, filename = NULL, ...) {
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- rlang::enexpr(id)
  if (length(id) != 1) id <- rlang::expr(gene)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "Corrplot")(data_select, ...)
}


#'Significance tests
#'
#'\code{proteoSigtest} performs significance tests aganist peptide or protein
#'\code{log2FC}. Users should avoid calling the method directly, but instead use
#'the following wrappers.
#'
#'In general, special characters of \code{+} or \code{-} should be avoided from
#'contrast terms. Occasionally, such as in biological studies, it may be
#'convenient to use \code{A+B} to denote a condition of combined treatment of
#'\code{A} and \code{B} . In the case, one can put the term(s) containing
#'\code{+} or \code{-} into a pair of pointy brackets. The syntax in the
#'following hypothetical example will compare the effects of \code{A}, \code{B},
#'\code{A+B} and the average of \code{A} and \code{B} to control \code{C}:
#'
#'\code{prnSig(fml = ~ Term["A - C", "B - C", "<A + B> - C", "(A + B)/2 - C"])}
#'
#'Note that \code{<A + B>} stands for one sample and \code{(A + B)} has two
#'samples in it.
#'
#'@inheritParams  prnHist
#'@inheritParams  prnHM
#'@param filename A file name to output results. The default is
#'  \code{Peptide_pVals.txt} for peptides and \code{Protein_pVals} for proteins.
#'@param method Character string; the method of linear modeling. The default is
#'  \code{limma}. At \code{method = lm}, the \code{lm()} in base R will be used
#'  for models without random effects and the \code{\link[lmerTest]{lmer}} will
#'  be used for models with random effects.
#'@param var_cutoff Numeric; the cut-off in the variances of \code{log2FC}.
#'  Entries with variances smaller than the threshold will be removed from
#'  linear modeling. The default is 1E-3.
#'@param pval_cutoff Numeric; the cut-off in significance \code{pVal}. Entries
#'  with \code{pVals} smaller than the threshold will be removed from multiple
#'  test corrections. The default is at \code{1} to include all entries.
#'@param logFC_cutoff Numeric; the cut-off in \code{log2FC}. Entries with
#'  absolute \code{log2FC} smaller than the threshold will be removed from
#'  multiple test corrections. The default is at \code{log2(1)} to include all
#'  entries.
#'@param ... User-defined formulae for linear modeling. The syntax starts with a
#'  tilde, followed by the name of an available column key in
#'  \code{expt_smry.xlsx} and square brackets. The contrast groups are then
#'  quoted with one to multiple contrast groups separated by commas. The default
#'  column key is \code{Term} in `expt_smry.xlsx`: \cr \code{~ Term["A - C", "B
#'  - C"]}. \cr Additive random effects are indicated by \code{+ (1|col_key_1) +
#'  (1|col_key_2)}... Currently only a syntax of single contrast are supported
#'  for uses with random effects: \cr \code{~ Term["A - C"] + (1|col_key_1) +
#'  (1|col_key_2)} \cr \cr \code{filter_}: Logical expression(s) for the row
#'  filtration of data; also see \code{\link{normPSM}}.
#'@return The primary output is
#'  \code{~\\dat_dir\\Peptide\\Model\\Peptide_pVals.txt} for peptide data or
#'  \code{~\\dat_dir\\Protein\\Model\\Protein_pVals.txt} for protein data. At
#'  \code{impute_na = TRUE}, the corresponding outputs are
#'  \code{Peptide_impNA_pvals.txt} or \code{Protein_impNA_pvals.txt}.
#'
#'@example inst/extdata/examples/prnSig_.R
#'@seealso \code{\link{load_expts}} for a reduced working example in data normalization \cr
#'
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
#'  \code{\link{mergePep}} for extended examples in peptide data merging \cr 
#'  \code{\link{standPep}} for extended examples in peptide data normalization \cr
#'  \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
#'  \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
#'  
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPcA}} for PCA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for protein trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for protein NMF analysis and visualization \cr 
#'  
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoSigtest <- function (df = NULL, id = gene, filepath = NULL, filename = NULL,
                           scale_log2r = TRUE, impute_na = TRUE, complete_cases = FALSE, method = "limma",
                           var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), ...) {
  
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% c(dots) %>% save_call("pepSig")
    } else if (id %in% c("prot_acc", "gene")) {
      load(file.path(dat_dir, "Calls\\prnSig_formulas.rda"))
      dots <- my_union(dots, prnSig_formulas)
      mget(names(formals()), current_env()) %>% c(dots) %>% save_call("prnSig")
    }
    , add = TRUE
  )
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  dots <- rlang::enexprs(...)
  
  id <- rlang::enexpr(id)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  method <- rlang::as_string(rlang::enexpr(method))
  
  reload_expts()
  
  if (!impute_na & method != "limma") {
    impute_na <- TRUE
    warning("Coerce `impute_na = ", impute_na, "` at method = ", method, call. = FALSE)
  }
  
  # Sample selection criteria:
  #   !is_reference under "Reference" ->
  #   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
  info_anal(df = !!df, id = !!id, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            filepath = !!filepath, filename = !!filename, 
            anal_type = "Model")(method = method, var_cutoff, pval_cutoff, logFC_cutoff, ...)
}



#'Volcano plot visualization
#'
#'\code{proteoVolcano} visualizes the volcano plots of peptide or protein data,
#'or protein subgroups under the same gene sets. Users should avoid call the
#'method directly, but instead use the following wrappers.
#'
#'@inheritParams  prnEucDist
#'@inheritParams  prnHM
#'@inheritParams  prnGSPA
#'@param scale_log2r Not used for full volcano plots by \code{pepVol} or 
#'\code{prnVol}; matched to the value in the latest call to \code{\link{prnGSPA}} 
#'for \code{gspaMap}. 
#'@param impute_na Logical. At the NULL default, the TRUE or FALSE will match
#'  the choice in \code{\link{pepSig}} for peptide and \code{\link{prnSig}} for
#'  protein data.
#'@param adjP Logical; if TRUE, use Benjamini-Hochberg pVals in volcano plots.
#'  The default is FALSE.
#'@param show_labels Logical; if TRUE, shows the labels of top twenty entries.
#'  The default is TRUE.
#'@param show_sig Character string indicating the type of significance values to
#'  be shown with \code{\link{gspaMap}}. The default is \code{"none"}.
#'  Additional choices are from \code{c("pVal", "qVal")} where \code{pVal} or
#'  \code{qVal} will be shown, respectively, in the facet grid of the plots.
#'@param pval_cutoff Numeric value or vector for uses with
#'  \code{\link{gspaMap}}. \code{Gene sets} with enrichment \code{pVals} less
#'  significant than the threshold will be excluded from volcano plot
#'  visualization. The default signficance is 0.05 for all formulae matched to
#'  or specified in argument \code{fml_nms}. Formula-specific threshold is
#'  allowed by supplying a vector of cut-off values.
#'@param logFC_cutoff Numeric value or vector for uses with
#'  \code{\link{gspaMap}}. \code{Gene sets} with absolute enrichment
#'  \code{log2FC} less than the threshold will be excluded from volcano plot
#'  visualization. The default magnitude is \code{log2(1.2) } for all formulae
#'  matched to or specified in argument \code{fml_nms}. Formula-specific
#'  threshold is allowed by supplying a vector of absolute values in
#'  \code{log2FC}.
#'@param gset_nms Character vector containing the name(s) of gene sets for uses
#'  with \code{\link{gspaMap}}. By default, the names will match those used in
#'  \code{\link{prnGSPA}}.
#'@inheritParams prnSig
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr \cr Additional parameters for
#'  plotting: \cr \code{xco}, the cut-off lines of fold changes at position
#'  \code{x}; the default is at \eqn{-1.2} and \eqn{+1.2}. \cr \code{yco}, the
#'  cut-off line of \code{pVal} at position \code{y}; the default is \eqn{0.05}.
#'  \cr \code{width}, the width of plot; \cr \code{height}, the height of plot.
#'  \cr \code{nrow}, the number of rows in a plot.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnVol_.R
#'@seealso \code{\link{load_expts}} for a reduced working example in data normalization \cr
#'
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
#'  \code{\link{mergePep}} for extended examples in peptide data merging \cr 
#'  \code{\link{standPep}} for extended examples in peptide data normalization \cr
#'  \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
#'  \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
#'  
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPcA}} for PCA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'  
#'@export
proteoVolcano <- function (id = "gene", anal_type = "Volcano", df = NULL, scale_log2r = TRUE,
                           filepath = NULL, filename = NULL, fml_nms = NULL, 
                           impute_na = NULL, adjP = FALSE, show_labels = TRUE, 
                           pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                           show_sig = "none", gset_nms = c("go_sets", "kegg_sets"), 
                           ...) {
  
  old_opt <- options(
    scipen = 0, 
    warn = 0, 
    max.print = 99999
  )
  
  options(
    scipen = 999,
    warn = 1,
    max.print = 2000000
  )
  
  on.exit(options(old_opt), add = TRUE)
  
  err_msg_1 <- "Unrecognized 'id'; needs to be one of \"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\" or \"term\""
  err_msg_2 <- "Unrecognized 'anal_type'; needs to be one of \"Volcano\" or \"GSPA\""
  err_msg_3 <- "Volcano plots of peptides not available for GSPA."
  err_msg_4 <- "GSPA results not found. Perform prnGSPA() first."
  err_msg_5 <- "Peptide_pVals.txt not found at impute_na = FALSE. Perform pepSig(impute_na = FALSE) first."
  err_msg_6 <- "Peptide_impNA_pVals.txt not found at impute_na = TRUE. Perform both pepImp() and pepSig(impute_na = TRUE) first."
  err_msg_7 <- "Protein_pVals.txt not found at impute_na = FALSE. Perform prnSig(impute_na = FALSE) first."
  err_msg_8 <- "Protein_impNA_pVals.txt not found at impute_na = TRUE. Perform both prnImp() and prnSig(impute_na = TRUE) first."
  
  # scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- rlang::as_string(rlang::enexpr(id))
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)	
  show_sig <- rlang::as_string(rlang::enexpr(show_sig))
  if (is.null(impute_na)) {
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      impute_na <- match_call_arg(pepSig, impute_na)
    } else if (id %in% c("prot_acc", "gene")) {
      impute_na <- match_call_arg(prnSig, impute_na)
    }
  }
  
  if (impute_na) {
    mscale_log2r <- match_call_arg(prnSig_impTRUE, scale_log2r)
  } else {
    mscale_log2r <- match_call_arg(prnSig_impFALSE, scale_log2r)
  }
  
  if (scale_log2r != mscale_log2r) {
    warning("scale_log2r = ", mscale_log2r, " after matching to `sigTest`.", call. = FALSE)
  }
  scale_log2r <- mscale_log2r
  rm(mscale_log2r)
  
  stopifnot(is_logical(scale_log2r))
  stopifnot(is_logical(impute_na))
  stopifnot(is_logical(adjP))
  stopifnot(is_logical(show_labels))
  stopifnot(is_double(pval_cutoff))
  stopifnot(is_double(logFC_cutoff))
  
  gset_nms <- local({
    file <- file.path(dat_dir, "Calls\\anal_prnGSPA.rda")
    
    if (file.exists(file)) {
      gset_nms <- gset_nms %>% 
        .[. %in% match_call_arg(anal_prnGSPA, gset_nms)]
      
      if (is.null(gset_nms)) stop ("Unknown gene sets.")
    }
    
    return(gset_nms)
  })
  
  if (! id %in% c("prot_acc", "gene", "pep_seq", "pep_seq_mod", "term")) stop(err_msg_1, call. = FALSE)
  if (! anal_type %in% c("Volcano", "GSPA", "mapGSPA")) stop(err_msg_2, call. = FALSE)
  
  if (id %in% c("prot_acc", "gene", "term")) {
    data_type <- "Protein"
  } else if (id %in% c("pep_seq", "pep_seq_mod")) {
    data_type <- "Peptide"
  }
  
  if (is.null(filepath)) {
    filepath = file.path(dat_dir, data_type, anal_type)
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (is.null(filename)) {
    fn_prefix <- paste(data_type, anal_type, sep = "_")
    fn_suffix <- "png"
  } else {
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
  }
  
  fn_prefix <- fn_prefix %>% 
    ifelse(impute_na, paste0(., "_impNA"), .) 	
  
  filename <- paste0(fn_prefix, ".", fn_suffix)
  
  if (is.null(df)) {
    if (anal_type %in% c("Volcano", "mapGSPA")) {
      if (id %in% c("pep_seq", "pep_seq_mod")) {
        fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
        fn_imp_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_impNA_pvals.txt")
        src_path <- ifelse(impute_na, fn_imp_p, fn_p)
        
        if (!file.exists(src_path)) {
          if (!impute_na) stop(err_msg_5, call. = FALSE) else stop(err_msg_6, call. = FALSE) 
        }
      } else if (id %in% c("prot_acc", "gene")) {
        fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
        fn_imp_p <- file.path(dat_dir, "Protein\\Model", "Protein_impNA_pvals.txt")
        src_path <- ifelse(impute_na, fn_imp_p, fn_p)
        
        if (!file.exists(src_path)) {
          if (!impute_na) stop(err_msg_7, call. = FALSE) else stop(err_msg_8, call. = FALSE) 
        }
      }
    } else if (anal_type %in% c("GSPA")) {
      if (id %in% c("pep_seq", "pep_seq_mod")) {
        stop(err_msg_3, call. = FALSE)
      } else if (id %in% c("prot_acc", "gene", "term")) {
        fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
        fn_imp_p <- file.path(dat_dir, "Protein\\Model", "Protein_impNA_pvals.txt")
        src_path <- ifelse(impute_na, fn_imp_p, fn_p)
        
        if (!file.exists(src_path)) {
          if (!impute_na) stop(err_msg_7, call. = FALSE) else stop(err_msg_8, call. = FALSE) 
        }
      }
    }
    
    df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE,
                            sep = "\t", comment.char = "#"), error = function(e) NA)
    
    if (!is.null(dim(df))) {
      message(paste("File loaded:", gsub("\\\\", "/", src_path)))
    } else {
      stop(paste("Non-existed file or directory:", gsub("\\\\", "/", src_path)))
    }
    
    df <- df %>%
      `rownames<-`(.[, id]) %>%
      dplyr::select(-grep("I[0-9]{3}|log2_R[0-9]{3}|^FC", names(.)))
  } else {
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      fn_raw <- file.path(dat_dir, "Peptide\\Model", df)
    } else if (id %in% c("prot_acc", "gene")) {
      fn_raw <- file.path(dat_dir, "Protein\\Model", df)
    }
    
    df <- tryCatch(read.csv(fn_raw, check.names = FALSE, header = TRUE, sep = "\t",
                            comment.char = "#"), error = function(e) NA)
    
    if (!is.null(dim(df))) {
      message(paste("File loaded:", gsub("\\\\", "/", fn_raw)))
    } else {
      stop(paste("Non-existed file or directory:", gsub("\\\\", "/", fn_raw)))
    }
  }
  
  df <- df %>% rm_pval_whitespace()
  
  species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
  if (!is_empty(species)) load_dbs(gset_nms = gset_nms, species = species)
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls, fml_nms, dots)
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  if (!adjP) {
    df <- df %>%
      dplyr::select(-contains("adjP"))
  } else {
    df <- df %>%
      dplyr::select(-contains("pVal")) %>%
      `names<-`(gsub("adjP", "pVal", names(.)))
  }
  
  plotVolcano(df, !!id, filepath, filename, adjP, show_labels, anal_type,
              pval_cutoff, logFC_cutoff, show_sig, gset_nms, 
              scale_log2r, impute_na, 
              !!!dots)
}










