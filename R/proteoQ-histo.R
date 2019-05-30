#' Plots histograms
#'
#' @import dplyr purrr rlang mixtools ggplot2 RColorBrewer
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
plotHisto <- function (df = NULL, label_scheme_sub, params, scale_log2r, show_curves,
                       show_vline, filepath = NULL, filename, ...) {

	dots <- rlang::enexprs(...)

	xmin <- eval(dots$xmin, env = caller_env())
	xmax <- eval(dots$xmax, env = caller_env())
	x_breaks <- eval(dots$x_breaks, env = caller_env())
	binwidth <- eval(dots$binwidth, env = caller_env())
	alpha <- eval(dots$alpha, env = caller_env())
	ncol <- eval(dots$ncol, env = caller_env())
	width <- eval(dots$width, env = caller_env())
	height <- eval(dots$height, env = caller_env())

	if(is.null(xmin)) xmin <- -2
	if(is.null(xmax)) xmax <- 2
	if(is.null(x_breaks)) x_breaks <- 1
	if(is.null(binwidth)) binwidth <- (xmax - xmin)/80
	if(is.null(alpha)) alpha <- .8
	if(is.null(ncol)) ncol <- 1
	if(is.null(width)) width <- 4 * ncol + 2
	if(is.null(height)) height <- length(label_scheme_sub$Sample_ID) * 4 / ncol
	
	nm_idx <- names(dots) %in% c("xmin", "xmax", "x_breaks", "binwidth",
	                             "ncol", "alpha", "width", "height")
	dots[nm_idx] <- NULL

	by = (xmax - xmin)/200
	nrow <- nrow(df)
	x_label <- expression("Ratio ("*log[2]*")")
	NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")

	if(!is.null(params)) {
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

	filename <- gg_imgname(filename)
	ggsave(file.path(filepath, filename), p, width = width, height = height, limitsize = FALSE, units = "in")
}



#'Visualizes histograms
#'
#'\code{proteoHist} prepares the histogram visualization of \code{log2-ratios}
#'for proteins or peptides data.
#'
#'In the histograms, the \code{log2-ratios} under each TMT channel are
#'color-coded by their contributing reporter-ion intensity.
#'
#'The function matches the current \code{id} to those in the latest \code{call}
#'to \code{\link{normPep}} or \code{\link{normPrn}}.  For example, if
#'\code{pep_seq} was used in \code{\link{normPep()}}, the current \code{id =
#'pep_seq_mod} will be matched to \code{id = pep_seq}.
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at \code{id = pep_seq} or \code{pep_seq_mod}, and protein data at
#'  \code{id = prot_acc} or \code{gene}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param show_curves Logical; if TRUE, shows the fitted curves.
#'@param show_vline Logical; if TRUE, shows the vertical lines at \code{x = 0}.
#'@param  new_fit Not currently used.
#'@param df The file name of input data. By default, it will be determined by
#'  the value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the name of the current function \code{call} and the value of
#'  \code{id}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}. The images are
#'  saved via \code{\link[ggplot2]{ggsave}} and the image type will be
#'  determined by the extension of file names with the default being \code{.png}.
#'@param ... Additional parameters for plotting: \cr \code{xmin}, the minimum x
#'  at a log2 scale; the default is -2. \cr \code{xmax}, the maximum x at a log2
#'  scale; the default is +2. \cr \code{x_breaks}, the breaks in x-axis at a
#'  log2 scale; the default is 1. \cr \code{binwidth}, the binwidth of
#'  \code{log2-ratios}; the default is (xmax - xmin)/80. \cr \code{ncol}, the
#'  number of columns; the default is 1. \cr \code{width}, the width of plot; \cr \code{height},
#'  the height of plot.
#'@return The histograms of \code{log2-ratios} under a "\code{Histogram}" sub
#'  directory.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoHist <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), col_select = NULL,
                        scale_log2r = FALSE, show_curves = TRUE, show_vline = TRUE, new_fit = FALSE,
                        df = NULL, filepath = NULL, filename = NULL, ...) {

  id <- rlang::enexpr(id)
	if(length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))

	col_select <- rlang::enexpr(col_select)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()

	info_anal(id = !!id, col_select = !!col_select, scale_log2r = scale_log2r, impute_na = FALSE,
	          df = !!df, filepath = !!filepath, filename = !!filename,
	          anal_type = "Histogram")(new_fit = new_fit, show_curves = show_curves,
	                                   show_vline = show_vline, ...)
}





#'Visualizes the histograms of peptide \code{log2FC}
#'
#'\code{pepHist} is a wrapper function of \code{\link{proteoHist}}
#'
#'@rdname proteoHist
#'
#' @examples
#' # without scaling normalization
#' pepHist(
#'   scale_log2r = FALSE,
#'   xmin = -1,
#'   xmax = 1,
#'   ncol = 5
#' )
#'
#' # with scaling normalization
#' pepHist(
#'   scale_log2r = TRUE,
#'   xmin = -1,
#'   xmax = 1,
#'   ncol = 5
#' )
#'
#' \dontrun{
#' pepHist(
#'   col_select = "a_non_existed_column_key"
#' )
#' }
#'
#'@export
pepHist <- function (...) {
	proteoHist(id = pep_seq, ...)
}


#'Visualizes the histograms of protein \code{log2FC}
#'
#'\code{prnHist} is a wrapper function of \code{\link{proteoHist}}
#'
#'@rdname proteoHist
#'
#' @examples
#' # no scaling normalization
#' prnHist(
#'   ncol = 10
#' )
#'
#' # scaling normalization
#' prnHist(
#'   scale_log2r = TRUE,
#'   ncol = 10
#' )
#'
#'
#'@export
prnHist <- function (...) {
	proteoHist(id = gene, ...)
}
