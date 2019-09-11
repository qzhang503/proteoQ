#' Plots histograms
#'
#' @import dplyr purrr rlang mixtools ggplot2 RColorBrewer
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
plotHisto <- function (df = NULL, id, label_scheme_sub, params, scale_log2r, pep_pattern, 
                       show_curves, show_vline, filepath = NULL, filename, ...) {

  stopifnot(nrow(label_scheme_sub) > 0)
  stopifnot(rlang::is_logical(scale_log2r))
  stopifnot(!grepl("[^A-z_]", pep_pattern))
  stopifnot(rlang::is_logical(show_curves))
  stopifnot(rlang::is_logical(show_vline))

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

	xmin <- eval(dots$xmin, env = caller_env()) # `xmin = -1` is `language`
	xmax <- eval(dots$xmax, env = caller_env()) # `xmax = +1` is `language`
	x_breaks <- eval(dots$x_breaks, env = caller_env())
	binwidth <- eval(dots$binwidth, env = caller_env())
	alpha <- eval(dots$alpha, env = caller_env())
	ncol <- eval(dots$ncol, env = caller_env())
	width <- eval(dots$width, env = caller_env())
	height <- eval(dots$height, env = caller_env())

	if (is.null(xmin)) xmin <- -2
	if (is.null(xmax)) xmax <- 2
	if (is.null(x_breaks)) x_breaks <- 1
	if (is.null(binwidth)) binwidth <- (xmax - xmin)/80
	if (is.null(alpha)) alpha <- .8
	if (is.null(ncol)) ncol <- 5
	if (is.null(width)) width <- 4 * ncol + 2
	if (is.null(height)) height <- length(label_scheme_sub$Sample_ID) * 4 / ncol
	
	nm_idx <- names(dots) %in% c("xmin", "xmax", "x_breaks", "binwidth",
	                             "ncol", "alpha", "width", "height")
	dots[nm_idx] <- NULL

	lang_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	dots <- dots %>% .[! . %in% lang_dots]
	df <- df %>% filters_in_call(!!!lang_dots)
	
	by = (xmax - xmin)/200
	nrow <- nrow(df)
	x_label <- expression("Ratio ("*log[2]*")")
	NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")

	if (!is.null(params)) {
		n_comp <- max(params$Component)
		nm_comps <- paste0("G", 1:n_comp)
		nm_full <- c(nm_comps, paste(nm_comps, collapse = " + "))

		# Offset by the percentage of non-NA values
		perc_nna <- df %>%
			dplyr::select(grep(paste0(NorZ_ratios, "[0-9]{3}"), names(.))) %>%
			`names<-`(gsub(".*_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>%
			lapply(function(x) sum(!is.na(x)) / length(x) * nrow * binwidth)

		perc_nna <- perc_nna[names(perc_nna) %in% label_scheme_sub$Sample_ID]

		# Density profiles
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


	my_theme <- theme_bw() + theme(
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

	seq <- c(-Inf, seq(4, 7, .5), Inf)

	if (pep_pattern != "zzz") {
	  df <- df %>% 
	    dplyr::filter(grepl(paste0("[", pep_pattern, "]"), !!rlang::sym(id)))
	  stopifnot(nrow(df) > 0)
	}
	
	df_melt <- df %>%
		dplyr::select(grep(paste0(NorZ_ratios, "[0-9]{3}", "|^N_I[0-9]{3}"), names(.))) %>%
		dplyr::filter(rowSums(!is.na(.[, grepl("[IR][0-9]{3}", names(.))])) > 0) %>%
		dplyr::select(which(not_all_zero(.))) %>%
		dplyr::select(which(colSums(!is.na(.)) > 0)) %>%
		dplyr::mutate(Int_index = log10(rowMeans(.[, grepl("^N_I[0-9]{3}", names(.))], na.rm = TRUE))) %>%
		dplyr::mutate_at(.vars = "Int_index", cut, seq, labels = seq[1:(length(seq)-1)]) %>%
		dplyr::select(-grep("^N_I[0-9]{3}", names(.))) %>%
		`names<-`(gsub(".*log2_R[0-9]{3}.*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
		tidyr::gather(key = Sample_ID, value = value, -Int_index) %>%
		dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
		dplyr::arrange(Sample_ID) %>%
		dplyr::filter(!is.na(value), !is.na(Int_index)) %>%
		dplyr::filter(Sample_ID %in% label_scheme_sub$Sample_ID)

	p <- ggplot() +
		geom_histogram(data = df_melt, aes(x = value, y = ..count.., fill = Int_index),
		               color = "white", alpha = alpha, binwidth = binwidth, size = .1) +
		scale_fill_brewer(palette = "Spectral", direction = -1) +
		labs(title = "", x = x_label, y = expression("Frequency")) +
		scale_x_continuous(limits = c(xmin, xmax), breaks = seq(xmin, xmax, by = x_breaks),
		                   labels = as.character(seq(xmin, xmax, by = x_breaks))) +
		my_theme +
		facet_wrap(~ Sample_ID, ncol = ncol, scales = "free")

	if(show_curves) p <- p + geom_line(data = fit, mapping = aes(x = x, y = value, colour = variable),
	                                   size = .2) +
											scale_colour_manual(values = myPalette, name = "Gaussian",
											                    breaks = c(nm_comps, paste(nm_comps, collapse = " + ")),
											                    labels = nm_full)

	if(show_vline) p <- p + geom_vline(xintercept = 0, size = .25, linetype = "dashed")

	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), 
	             width = width, height = height, limitsize = FALSE, dots)
	do.call(ggsave, gg_args)
}



#'Histogram visualization
#'
#'\code{proteoHist} plots the histograms of protein or peptide \code{log2FC}.
#'Users should avoid call the method directly, but instead use the following
#'wrappers.
#'
#'In the histograms, the \code{log2FC} under each TMT channel are color-coded by
#'their contributing reporter-ion intensity.
#'
#'The function matches the current \code{id} to the grouping argument in the
#'latest \code{call} to \code{\link{normPSM}} or \code{\link{normPep}}.  For
#'example, if \code{normPSM(group_psm_by = pep_seq, ...)} was called earlier,
#'the setting of \code{id = pep_seq_mod} in the current call will be matched to
#'\code{id = pep_seq}. Similarly, if \code{normPep(group_pep_by = gene, ...)}
#'was employed, the setting of \code{id = prot_acc} in the current call will be
#'matched to \code{id = gene}.
#'
#'@param id Character string to indicate the type of data. The value will be
#'  determined automatically by the program. Peptide data will be used at
#'  \code{id = pep_seq} or \code{pep_seq_mod}, and protein data will be used at
#'  \code{id = prot_acc} or \code{gene}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  Samples corresponding to non-empty entries under the column key will be
#'  included in the indicated analysis. At the NULL default, the column key will
#'  be \code{Select}.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2FC} to the same scale
#'  of standard deviation for all samples.
#'@param pep_pattern Character string containing one-letter representation of
#'  amino acids. At the "zzz" default, all peptides will be used. Letters in the
#'  character string are case sensitive. For example, \code{pep_pattern = "y"}
#'  will extract peptides with tyrosine phosphorylation. To extract peptides
#'  with N-terminal acetylation, use \code{pep_pattern = "_"}. The parameter
#'  provides a means for high-level subsetting of peptide entries in a data set.
#'  In general, one can use the \code{filter-in-function} feature described in
#'  \code{?normPSM} to subset data.
#'
#'@param show_curves Logical; if TRUE, shows the fitted curves. The curve
#'  parameters are based on the latest call to \code{normPep} or \code{normPrn}.
#'  This feature can inform the effects of data filtration on the alignment of
#'  \code{logFC} profiles.
#'@param show_vline Logical; if TRUE, shows the vertical lines at \code{x = 0}.
#'@param df The name of input data file. By default, it will be determined
#'  automatically by the value of \code{id}.
#'@param filepath A file path to output results. By default, it will be
#'  determined automatically by the name of the calling function and the value
#'  of \code{id} in the \code{call}.
#'@param filename A representative file name to output image(s). By default, it
#'  will be determined automatically by the name of the current \code{call}. The
#'  image(s) are saved via \code{\link[ggplot2]{ggsave}} where the image type
#'  will be determined by the extension of the file name. A \code{.png} format
#'  will be used at default or an unrecognized file extension.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr Additional parameters for
#'  plotting: \cr \code{xmin}, the minimum \eqn{x} at a log2 scale; the default
#'  is -2. \cr \code{xmax}, the maximum \eqn{x} at a log2 scale; the default is
#'  +2. \cr \code{x_breaks}, the breaks in \eqn{x}-axis at a log2 scale; the
#'  default is 1. \cr \code{binwidth}, the binwidth of \code{log2FC}; the
#'  default is \eqn{(xmax - xmin)/80}. \cr \code{ncol}, the number of columns;
#'  the default is 1. \cr \code{width}, the width of plot; \cr \code{height},
#'  the height of plot.
#'@return The histograms of \code{log2FC} under
#'  \code{~\\dat_dir\\Peptide\\Histogram} or
#'  \code{~\\dat_dir\\Protein\\Histogram}.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoHist <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                        col_select = NULL, scale_log2r = FALSE, pep_pattern = "zzz", 
                        show_curves = TRUE, show_vline = TRUE, 
                        df = NULL, filepath = NULL, filename = NULL, ...) {

  id <- rlang::enexpr(id)
	if(length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))

	col_select <- rlang::enexpr(col_select)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	pep_pattern <- rlang::as_string(rlang::enexpr(pep_pattern))
	
	reload_expts()

	info_anal(id = !!id, col_select = !!col_select, scale_log2r = scale_log2r, impute_na = FALSE,
	          df = !!df, filepath = !!filepath, filename = !!filename,
	          anal_type = "Histogram")(pep_pattern = pep_pattern, show_curves = show_curves,
	                                   show_vline = show_vline, ...)
}


#'Visualizes the histograms of peptide \code{log2FC}
#'
#'\code{pepHist} is a wrapper of \code{\link{proteoHist}} for peptide data
#'
#'@rdname proteoHist
#'
#' @examples
#' ## visualization of normalization
#' # without scaling
#' pepHist(scale_log2r = FALSE)
#'
#' # with scaling
#' pepHist(scale_log2r = TRUE)
#'
#' ## samples for use are indicated under a column, i.e. `Select`, in `expt_smry.xlsx`
#' # peptides
#' pepHist(col_select = Select)
#'
#' # proteins
#' prnHist(col_select = Select)
#'
#' ## data filtration
#' # phosphopeptide subset
#' pepHist(
#'   pep_pattern = "sty",
#'   filename = "pepHist_fil_by_sty.png",
#' )
#'
#' # recommended way to extract phosphopeptides
#' pepHist(
#'   filter_by = exprs(grepl("[sty]", pep_seq_mod)),
#'   filename = "pepHist_fil_by_sty.png",
#' )
#'
#' # exclude oxidized methione or deamidated asparagine
#' pepHist(
#'   filter_by = exprs(!grepl("[mn]", pep_seq_mod)),
#'   filename = "pepHist_fil_no_mn.png",
#' )
#'
#' ## between lead and lag
#' # leading logFC profiles of peptides
#' pepHist(
#'   filename = "pepHist.png",
#' )
#'
#' # lagging logFC profiles of peptides at
#' #   (1) n_psm >= 10
#' #   (2) and no methionine oxidation or asparagine deamidation
#' # may exclude sample(s) with considerable offset(s) 
#' #   between the 'bootstrapped' lead and lag from further analysis
#' pepHist(
#'   filter_by_npsm = exprs(n_psm >= 10),
#'   filter_by_mn = exprs(!grepl("[mn]", pep_seq_mod)),
#'   filename = "pepHist_filtered.png",
#' )
#'
#' \dontrun{
#' # sample selection
#' pepHist(
#'   col_select = "a_column_key_not_in_`expt_smry.xlsx`",
#' )
#'
#' # data filtration
#' pepHist(
#'   filter_by = exprs(!grepl("[m]", a_column_key_not_in_data_table)),
#' )
#'
#' prnHist(
#'   lhs_not_start_with_filter_ = exprs(n_psm >= 5),
#' )
#'
#' }
#'
#'@import purrr
#'@export
pepHist <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\Histogram\\log"), recursive = TRUE, showWarnings = FALSE)

  quietly_log <- purrr::quietly(proteoHist)(id = pep_seq, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\Histogram\\log","pepHist_log.csv"), append = TRUE)  
}


#'Visualizes the histograms of protein \code{log2FC}
#'
#'\code{prnHist} is a wrapper of \code{\link{proteoHist}} for protein data
#'
#'@rdname proteoHist
#'
#'@import purrr
#'@export
prnHist <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Histogram\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoHist)(id = gene, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Histogram\\log","prnHist_log.csv"), append = TRUE)
}
