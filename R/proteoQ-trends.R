#' Plots trends
#'
#' @import Mfuzz plyr dplyr rlang Biobase ggplot2 RColorBrewer
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
plotTrend <- function (df, id, col_group, col_order, label_scheme_sub, n_clust,
                       complete_cases, scale_log2r, filepath, filename, ...) {

	mestimate <- function (eset) {
		N <- dim(Biobase::exprs(eset))[[1]]
		D <- dim(Biobase::exprs(eset))[[2]]
		m.sj <- 1 + (1418/N + 22.05) * D^(-2) + (12.33/N + 0.243) *
			D^(-0.0406 * log(N) - 0.1134)
		return(m.sj)
	}

	mfuzz <- function (eset, centers, m, ...) {
		cl <- e1071::cmeans(Biobase::exprs(eset), centers = centers,
		                    method = "cmeans", m = m, ...)
	}


	stopifnot(nrow(label_scheme_sub) > 0)

	id <- rlang::as_string(rlang::enexpr(id))

	dots <- rlang::enexprs(...)

	col_group <- rlang::enexpr(col_group)
	col_order <- rlang::enexpr(col_order)

	ymin <- eval(dots$ymin, env = caller_env())
	ymax <- eval(dots$ymax, env = caller_env())
	y_breaks <- eval(dots$y_breaks, env = caller_env())
	ncol <- dots$ncol
	nrow <- dots$nrow
	width <- dots$width
	height <- dots$height

	if(is.null(ymin)) ymin <- -2
	if(is.null(ymax)) ymax <- 2
	if(is.null(y_breaks)) y_breaks <- 1
	if(is.null(ncol)) ncol <- 1
	if(is.null(nrow)) nrow <- 2

	dots$ymin <- NULL
	dots$ymax <- NULL
	dots$y_breaks <- NULL
	dots$ncol <- NULL
	dots$nrow <- NULL

	x_label <- expression("Ratio ("*log[2]*")")
	fn_prx <- gsub("\\..*$", "", filename)
	fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)

	df_mean <- t(df) %>%
		data.frame(check.names = FALSE) %>%
		tibble::rownames_to_column("Sample_ID") %>%
		dplyr::left_join(label_scheme_sub, by = "Sample_ID") %>%
		dplyr::filter(!is.na(!!col_group)) %>%
		dplyr::mutate(Group := !!col_group) %>%
		dplyr::group_by(Group) %>%
		dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

	if(rlang::as_string(col_order) %in% names(df_mean)) {
		df_mean <- df_mean %>%
			dplyr::arrange(!!col_order) %>%
			dplyr::rename(Order := !!col_order) %>%
			dplyr::select(-Order, -TMT_Set)
	}

	df_mean <- t(df_mean[, -1]) %>%
		data.frame(check.names = FALSE) %>%
		`colnames<-`(df_mean$Group) %>%
		tibble::rownames_to_column(id) %>%
		dplyr::filter(complete.cases(.[, !grepl(id, names(.))])) %>%
		tibble::column_to_rownames(id)

	if(!file.exists(file.path(filepath, paste0(fn_prx, "-",n_clust,"_Clusters.csv")))) {
  	df_fuzzy = new('ExpressionSet', exprs = as.matrix(df_mean))
  	m1 <- mestimate(df_fuzzy)

  	cl <- mfuzz(df_fuzzy, c = n_clust, m = m1)
  	O <- Mfuzz::overlap(cl)

  	df_op <- data.frame(cl_cluster = cl$cluster) %>%
  		tibble::rownames_to_column() %>%
  		dplyr::rename(!!id := rowname) %>%
  		dplyr::left_join(df %>% tibble::rownames_to_column(id), by = id)

  	write.csv(df_op, file.path(filepath, paste0(fn_prx, "-",n_clust,"_Clusters.csv")),
  	          row.names = FALSE)
	} else {
	  df_op <- read.csv(file.path(filepath, paste0(fn_prx, "-",n_clust,"_Clusters.csv")),
	                 check.names = FALSE, header = TRUE, comment.char = "#")
	}

	Levels <- names(df_mean)
	df_mean <- df_mean %>%
		tibble::rownames_to_column(id) %>%
		left_join(df_op %>% dplyr::select(id, "cl_cluster"), by = id) %>%
		tidyr::gather(key = variable, value = value, -id, -cl_cluster) %>%
		dplyr::mutate(variable = factor(variable, levels = Levels)) %>%
		dplyr::arrange(variable)
	rm(Levels)

	if(is.null(width)) width <- n_clust * 8 / nrow
	if(is.null(height)) height <- 8 * nrow

	dots$width <- NULL
	dots$height <- NULL

	my_theme <- theme_bw() + theme(
		axis.text.x  = element_text(angle=60, vjust=0.5, size=24),
		axis.ticks.x  = element_blank(), # x-axis ticks
		axis.text.y  = element_text(angle=0, vjust=0.5, size=24),
		axis.title.x = element_text(colour="black", size=24),
		axis.title.y = element_text(colour="black", size=24),
		plot.title = element_text(face="bold", colour="black",
		                          size=20, hjust=.5, vjust=.5),
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),
		panel.background = element_rect(fill = '#0868ac', colour = 'red'),

		strip.text.x = element_text(size = 24, colour = "black", angle = 0),
		strip.text.y = element_text(size = 24, colour = "black", angle = 90),

		legend.key = element_rect(colour = NA, fill = 'transparent'),
		legend.background = element_rect(colour = NA,  fill = "transparent"),
		legend.position = "none",
		legend.title = element_text(colour="black", size=18),
		legend.text = element_text(colour="black", size=18),
		legend.text.align = 0,
		legend.box = NULL
	)

	p <- ggplot(data = df_mean,
	            mapping = aes(x = variable, y = value, group = !!rlang::sym(id))) +
		geom_line(colour = "white", alpha = .25) +
	  # stat_density2d(aes(fill = ..density..),
	  #                geom = "raster", alpha = .75, contour = FALSE) +
		# scale_fill_distiller(limits = c(0,.5), palette = "Blues", direction = -1,
		#                      na.value = brewer.pal(n = 9, name = "Blues")[1]) +
	  scale_y_continuous(limits = c(ymin, ymax), breaks = c(ymin, 0, ymax)) +
	  labs(title = "", x = "", y = x_label) +
		my_theme
	p <- p + facet_wrap(~ cl_cluster, nrow = nrow, labeller = label_value)

	ggsave(file.path(filepath, paste0(fn_prx, "-", n_clust, "_Clusters.", fn_suffix)),
	       p, width = width, height = height, units = "in")

}


#'Clustering by trends
#'
#'Analyzes and visualizes the trend clustering of peptide or protein
#'\code{log2-ratios}
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at \code{id = pep_seq} or \code{pep_seq_mod}, and protein data at
#'  \code{id = prot_acc} or \code{gene}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param  col_order Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Order}. Data will be visualized by the order
#'  specified by the key.
#'@param impute_na Logical; if TRUE, imputes missing values.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used for visualization.
#'@param n_clust Numeric; the number of clusters that data will be divided into.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the names of the current functional \code{call}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... Additional parameters for plotting: \cr \code{ymin}, the minimum y;
#'  \cr \code{ymax}, the maximum y; \cr \code{y_breaks}, the breaks in y-axis;
#'  \cr \code{ncol}, the number of columns; \cr \code{nrow}, the number of rows;
#'  \cr \code{width}, the width of plot; \cr \code{height}, the height of plot.
#'@return Trend plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoTrend <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"),
                         col_select = NULL, col_order = NULL,
												impute_na = FALSE, complete_cases = FALSE,
												n_clust = 6, scale_log2r = TRUE, df = NULL,
												filepath = NULL, filename = NULL, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

	id <- rlang::enexpr(id)
	if(id == rlang::expr(c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))) {
		id <- "gene"
	} else {
		id <- rlang::as_string(id)
		stopifnot(id %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
	}

	col_select <- rlang::enexpr(col_select)
	col_order <- rlang::enexpr(col_order)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()

	if(!impute_na) complete_cases <- TRUE

	info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
	          scale_log2r = scale_log2r, impute_na = impute_na,
						df = !!df, filepath = !!filepath, filename = !!filename,
						anal_type = "Trend")(n_clust = n_clust, complete_cases = complete_cases, ...)
}


#'Trend analysis
#'
#'Trend analysis of protein \code{log2-ratios}
#'
#' @examples
#' prnTrend(
#'   scale_log2r = TRUE,
#'   col_order = Order,
#'   n_clust = 6
#' )
#'
#'@seealso \code{\link{proteoTrend}} for parameters
#'@export
prnTrend <- function (...) {
	proteoTrend(id = gene, ...)
}
