#' Plots histograms
#' 
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr purrr rlang mixtools ggplot2 RColorBrewer
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
plotHisto <- function (df = NULL, id, label_scheme_sub, scale_log2r, complete_cases, 
                       show_curves, show_vline, scale_y, filepath = NULL, filename, 
                       theme, ...) {

  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  if (scale_log2r) {
    fn_par <- file.path(filepath, "MGKernel_params_Z.txt")
  } else {
    fn_par <- file.path(filepath, "MGKernel_params_N.txt")
  }
  
  if (file.exists(fn_par)) {
    params <- read.csv(fn_par, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")
    params <- within(params, {mean = mean - x})
  } else {
    params <- NULL
    show_curves <- FALSE
  }

  stopifnot(nrow(label_scheme_sub) > 0)
  stopifnot(rlang::is_logical(scale_log2r))
  stopifnot(rlang::is_logical(show_curves))
  stopifnot(rlang::is_logical(show_vline))
  stopifnot(rlang::is_logical(scale_y))

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

	xmin <- eval(dots$xmin, env = caller_env()) # `xmin = -1` is `language`
	xmax <- eval(dots$xmax, env = caller_env()) 
	xbreaks <- eval(dots$xbreaks, env = caller_env())
	binwidth <- eval(dots$binwidth, env = caller_env())
	alpha <- eval(dots$alpha, env = caller_env())
	ncol <- eval(dots$ncol, env = caller_env())
	width <- eval(dots$width, env = caller_env())
	height <- eval(dots$height, env = caller_env())
	scales <- eval(dots$scales, env = caller_env())

	if (is.null(xmin)) xmin <- -2
	if (is.null(xmax)) xmax <- 2
	if (is.null(xbreaks)) xbreaks <- 1
	if (is.null(binwidth)) binwidth <- (xmax - xmin)/80
	if (is.null(alpha)) alpha <- .8
	if (is.null(ncol)) ncol <- 5
	if (is.null(width)) width <- 4 * ncol + 2
	if (is.null(height)) height <- length(label_scheme_sub$Sample_ID) * 4 / ncol
	if (is.null(scales)) scales <- "fixed"
  
	dots <- dots %>% 
    .[! names(.) %in% c("xmin", "xmax", "xbreaks", "binwidth", "ncol", "alpha", 
                        "width", "height", "scales")]

	filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
	dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
	
	if (scale_y) {
	  df <- df %>% 
	    filters_in_call(!!!filter_dots) %>% 
	    arrangers_in_call(!!!arrange_dots)
	  
	  stopifnot(nrow(df) > 0)
	}
	
	by = (xmax - xmin)/200
	nrow <- nrow(df)
	x_label <- expression("Ratio ("*log[2]*")")
	NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")

	if (!is.null(params)) {
		n_comp <- max(params$Component)
		nm_comps <- paste0("G", 1:n_comp)
		nm_full <- c(nm_comps, paste(nm_comps, collapse = " + "))

		# offset by the percentage of non-NA values
		perc_nna <- df %>%
			dplyr::select(grep(paste0(NorZ_ratios, "[0-9]{3}"), names(.))) %>%
			# `names<-`(gsub(".*_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>%
		  `names<-`(gsub(".*_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
			lapply(function(x) sum(!is.na(x)) / length(x) * nrow * binwidth)

		perc_nna <- perc_nna[names(perc_nna) %in% label_scheme_sub$Sample_ID]

		# density profiles
		fit <- params %>%
			dplyr::filter(.$Sample_ID %in% label_scheme_sub$Sample_ID) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) %>%
			split(.$Sample_ID) %>%
			lapply(sumdnorm, xmin, xmax, by = by)

		fit <- fit %>%
			purrr::map(~ .[, grep("^G[0-9]{1}|^Sum", names(.))]) %>%
			purrr::map2(perc_nna, `*`) %>%
			do.call(rbind, .) %>%
			dplyr::bind_cols(fit %>% do.call(rbind, .) %>% dplyr::select(c("x", "Sample_ID"))) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) %>%
			dplyr::rename(!!sym(paste(nm_comps, collapse = " + ")) := Sum) %>%
			tidyr::gather(key = variable, value = value, -x, -Sample_ID) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
			dplyr::mutate(variable = factor(variable, levels = nm_full)) %>%
			dplyr::arrange(Sample_ID, variable) %>%
			dplyr::filter(Sample_ID %in% label_scheme_sub$Sample_ID)

		myPalette <- c(rep("gray", n_comp), "black")
	} else {
		n_comp <- NA
		nm_comps <- NA
		myPalette <- NA
		nm_full <- NA
		perc_nna <- NA
		fit <- NA
	}


	proteoq_histo_theme <- theme_bw() + theme(
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
	
	if (is.null(theme)) theme <- proteoq_histo_theme

	seq <- c(-Inf, seq(4, 7, .5), Inf)

	if (!scale_y) {
	  df <- df %>% 
	  filters_in_call(!!!filter_dots) %>% 
	  arrangers_in_call(!!!arrange_dots)
	}

	stopifnot(nrow(df) > 0)
	
	df_melt <- df %>%
		dplyr::select(grep(paste0(NorZ_ratios, "[0-9]{3}", "|^N_I[0-9]{3}"), names(.))) %>%
		dplyr::filter(rowSums(!is.na(.[, grepl("[IR][0-9]{3}", names(.))])) > 0) %>%
		dplyr::select(which(not_all_zero(.))) %>%
		dplyr::select(which(colSums(!is.na(.)) > 0)) %>%
		dplyr::mutate(Int_index = log10(rowMeans(.[, grepl("^N_I[0-9]{3}", names(.))], na.rm = TRUE))) %>%
		dplyr::mutate_at(.vars = "Int_index", cut, seq, labels = seq[1:(length(seq)-1)]) %>%
		dplyr::select(-grep("^N_I[0-9]{3}", names(.))) %>%
		# `names<-`(gsub(".*log2_R[0-9]{3}.*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
	  `names<-`(gsub("^[NZ]{1}_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
		tidyr::gather(key = Sample_ID, value = value, -Int_index) %>%
		dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
		dplyr::arrange(Sample_ID) %>%
		dplyr::filter(!is.na(value), !is.na(Int_index)) %>%
		dplyr::filter(Sample_ID %in% label_scheme_sub$Sample_ID) %>% 
	  dplyr::mutate(value = setHMlims(value, xmin, xmax))

	p <- ggplot() +
		geom_histogram(data = df_melt, aes(x = value, y = ..count.., fill = Int_index),
		               color = "white", alpha = alpha, binwidth = binwidth, size = .1) +
		scale_fill_brewer(palette = "Spectral", direction = -1) +
		labs(title = "", x = x_label, y = expression("Frequency")) +
		scale_x_continuous(limits = c(xmin, xmax), breaks = seq(xmin, xmax, by = xbreaks),
		                   labels = as.character(seq(xmin, xmax, by = xbreaks))) +
		facet_wrap(~ Sample_ID, ncol = ncol, scales = scales) + 
	  theme

	if (show_curves) p <- p + geom_line(data = fit, mapping = aes(x = x, y = value, colour = variable),
	                                   size = .2) +
											scale_colour_manual(values = myPalette, name = "Gaussian",
											                    breaks = c(nm_comps, paste(nm_comps, collapse = " + ")),
											                    labels = nm_full)

	if (show_vline) p <- p + geom_vline(xintercept = 0, size = .25, linetype = "dashed")

	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), 
	             width = width, height = height, limitsize = FALSE, dots)
	do.call(ggsave, gg_args)
}


#'Histogram visualization
#'
#'\code{pepHist} plots the histograms of peptide \code{log2FC}.
#'
#'@rdname prnHist
#'@import purrr
#'@export
pepHist <- function (col_select = NULL, scale_log2r = TRUE, complete_cases = FALSE, 
                     show_curves = TRUE, show_vline = TRUE, scale_y = TRUE, 
                     df = NULL, filepath = NULL, filename = NULL, theme = NULL, ...) {
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  stopifnot(rlang::is_logical(show_curves), 
            rlang::is_logical(show_vline), 
            rlang::is_logical(scale_y))
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  dots <- rlang::enexprs(...)
  if (!is.null(dots$impute_na)) {
    dots$impute_na <- NULL
    rlang::warn("No NA imputation with histograms.")
  }

  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = FALSE,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Histogram")(show_curves = show_curves,
                                     show_vline = show_vline, scale_y = scale_y, 
                                     theme = theme, !!!dots)
}


#'Histogram visualization
#'
#'\code{prnHist} plots the histograms of protein \code{log2FC}.
#'
#'In the histograms, the \code{log2FC} under each TMT channel are color-coded by
#'their contributing reporter-ion intensity.
#'
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
#'@param df The name of a primary data file. By default, it will be determined
#'  automatically after matching the types of data and analysis with an
#'  \code{id} among \code{c("pep_seq", "pep_seq_mod", "prot_acc", "gene")}. A
#'  primary file contains normalized peptide or protein data and is among
#'  \code{c("Peptide.txt", "Peptide_pVal.txt", "Peptide_impNA_pVal.txt",
#'  "Protein.txt", "Protein_pVal.txt", "protein_impNA_pVal.txt")}. For analyses
#'  require the fields of significance p-values, the \code{df} will be one of
#'  \code{c("Peptide_pVal.txt", "Peptide_impNA_pVal.txt", "Protein_pVal.txt",
#'  "protein_impNA_pVal.txt")}.
#'@param filepath A file path to output results. By default, it will be
#'  determined automatically by the name of the calling function and the value
#'  of \code{id} in the \code{call}.
#'@param filename A representative file name to outputs. By default, the name(s)
#'  will be determined automatically. For text files, a typical file extension
#'  is \code{.txt}. For image files, they are typically saved via
#'  \code{\link[ggplot2]{ggsave}} or \code{\link[pheatmap]{pheatmap}} where the
#'  image type will be determined by the extension of the file name.
#'@param theme A
#'  \href{https://ggplot2.tidyverse.org/reference/ggtheme.html}{ggplot2}
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
#'@seealso 
#'  \emph{Metadata} \cr 
#'  \code{\link{load_expts}} for metadata preparation and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr 
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
#'  \emph{Variable arguments of `filter_...`} \cr 
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \emph{Missing values} \cr 
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  
#'  \emph{Informatics} \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{prepEntrez}} for lookups between UniProt accessions and Entrez IDs \cr
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnHist_.R
#'
#'@return Histograms of \code{log2FC}
#'@export
prnHist <- function (col_select = NULL, scale_log2r = TRUE, complete_cases = FALSE, 
                     show_curves = TRUE, show_vline = TRUE, scale_y = TRUE, 
                     df = NULL, filepath = NULL, filename = NULL, theme = NULL, ...) {
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))

  stopifnot(rlang::is_logical(show_curves), 
            rlang::is_logical(show_vline), 
            rlang::is_logical(scale_y))
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  dots <- rlang::enexprs(...)
  if (!is.null(dots$impute_na)) {
    dots$impute_na <- NULL
    rlang::warn("No NA imputation with histograms.")
  }

  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = FALSE,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Histogram")(show_curves = show_curves,
                                     show_vline = show_vline, scale_y = scale_y, 
                                     theme = theme, !!!dots)
}

