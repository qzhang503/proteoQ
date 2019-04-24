#' Make Trend lines from fuzzy c-mean clustering
#' 
#' @import Mfuzz plyr dplyr rlang Biobase ggplot2 RColorBrewer
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
plotTrend <- function (df, id, col_group, col_order, label_scheme_sub, n_clust, complete_cases, scale_log2r, filepath, filename, ...) {
	
	mestimate <- function (eset) {
		N <- dim(Biobase::exprs(eset))[[1]]
		D <- dim(Biobase::exprs(eset))[[2]]
		m.sj <- 1 + (1418/N + 22.05) * D^(-2) + (12.33/N + 0.243) * 
			D^(-0.0406 * log(N) - 0.1134)
		return(m.sj)
	}

	mfuzz <- function (eset, centers, m, ...) {
		cl <- e1071::cmeans(Biobase::exprs(eset), centers = centers, method = "cmeans", m = m, ...)
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
	
	# the mean of each sample group
	# -------------------------------------
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

	df_fuzzy = new('ExpressionSet', exprs = as.matrix(df_mean))
	m1 <- mestimate(df_fuzzy)

	# overlaps of clusters
	cl <- mfuzz(df_fuzzy, c = n_clust, m = m1)
	O <- Mfuzz::overlap(cl) # X11(); mfuzz.plot(df_fuzzy, cl=cl, mfrow=c(2,5)); Ptemp <- overlap.plot(cl,over=O,thres=0.05)

	df_op <- data.frame(cl_cluster = cl$cluster) %>% 
		tibble::rownames_to_column() %>% 
		dplyr::rename(!!id := rowname) %>% 
		dplyr::left_join(df %>% tibble::rownames_to_column(id), by = id)

	write.csv(df_op, file.path(filepath, paste0(fn_prx, "-",n_clust,"_Clusters.csv")), row.names = FALSE) 

	
	# plots
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
	
	# --------------------------------------------------------------------------------------------
	# myPalette <- rev(c(brewer.pal(n = 9, name = "Blues"), "#FFFFFF"))
	# myPalette <- colorRampPalette(c("white","royalblue","yellow"))
	# plot(rep(1,50),col=(myPalette(50)), pch=19,cex=2)			
	# myPalette <- c("0"="#FFFFFF", "0.02"="#08306B", "0.04"="#08306B", "0.06"="#08519C", "0.08"="#2171B5", "0.1"="#4292C6", "0.12"="#6BAED6", "0.14"="#9ECAE1", "0.16"="#C6DBEF", "0.18"="#DEEBF7", "0.2"="#F7FBFF")
		
	my_theme <- theme_bw() + theme(
		axis.text.x  = element_text(angle=60, vjust=0.5, size=24),   # rotate x-axis label and define label size
		axis.ticks.x  = element_blank(), # x-axis ticks
		axis.text.y  = element_text(angle=0, vjust=0.5, size=24),   # rotate y-axis label and define label size
		axis.title.x = element_text(colour="black", size=24),  # x-axis title size
		axis.title.y = element_text(colour="black", size=24),  # y-axis title size
		plot.title = element_text(face="bold", colour="black", size=20, hjust=.5, vjust=.5),   # Main title size
		
		panel.grid.major.x = element_blank(), 
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),
		
		panel.background = element_rect(fill = '#0868ac', colour = 'red'), 
		
		strip.text.x = element_text(size = 24, colour = "black", angle = 0), 
		strip.text.y = element_text(size = 24, colour = "black", angle = 90),   # Facet title size
		
		legend.key = element_rect(colour = NA, fill = 'transparent'), 
		legend.background = element_rect(colour = NA,  fill = "transparent"),
		# legend.position = c(0.15, 0.78),  
		legend.position = "none",  # to remove all legends
		# legend.title = element_blank(),
		legend.title = element_text(colour="black", size=18),
		legend.text = element_text(colour="black", size=18),
		legend.text.align = 0, 
		legend.box = NULL
	)
	
	p <- ggplot(data = df_mean, mapping = aes(x = variable, y = value, group = !!rlang::sym(id))) + 
		# stat_density2d(aes(fill = ..density..), geom = "raster", alpha = .75, contour = FALSE) + 
		geom_line(colour = "white", alpha = .25) +
		# scale_fill_distiller(limits = c(0,.5), palette = "Blues", direction = -1, na.value = brewer.pal(n = 9, name = "Blues")[1]) + 
		labs(title = "", x = "", y = x_label) + 
		my_theme
	p <- p + facet_wrap(~ cl_cluster, nrow = nrow, labeller = label_value)
	
	ggsave(file.path(filepath, paste0(fn_prx, "-", n_clust, "_Clusters.", fn_suffix)), p, width = width, height = height, units = "in")

}


#' Plots the trend of \code{log2-ratios} after fuzzy c-mean classification
#'
#' \code{proteoTrend} produces histograms.
#'
#' reads the data from either "\code{~\\Direcotry\\Peptide\\Peptide All.txt}" at
#' \code{id = pep_seq_mod}, or "\code{~\\Direcotry\\Protein\\Protein All by
#' Accession.txt}" at \code{id = prot_acc} or
#' "\code{~\\Direcotry\\Protein\\Protein All by gene.txt}" at \code{id = gene}.
#'
#' @param id The name of a unique identifier (see \code{\link[proteoQ]{MDS}}).
#' @param scale_log2r Logical; if TRUE, rescales \code{log2-ratios} to the same
#'   scale of standard deviation for all samples.
#' @param annot_kinases Logical; if TRUE, annotates proteins being kinases or
#'   not.
#' @return Images stored under the file folders that are associated to
#'   \code{id}, \code{anal_type} and \code{annot_kinases}.
#'
#' @examples
#' MA(
#' 	id = gene,
#' 	scale_log2r = scale_log2r,
#' 	annot_kinases = annot_kinases,
#' )
#'
#' \dontrun{
#' }
#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
#' @export
proteoTrend <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), col_select = NULL, col_order = NULL, 
												impute_na = FALSE, complete_cases = FALSE, 
												n_clust = 6, scale_log2r = FALSE, df = NULL, filepath = NULL, filename = NULL, ...) {
	
	id <- rlang::enexpr(id)
	if(id == rlang::expr(c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))) {
		id <- "gene"
	} else {
		id <- rlang::as_string(id)
		stopifnot(id %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
	}
	
	col_select <- rlang::enexpr(col_select)
	col_order <- rlang::enexpr(col_order)
	
	if(!impute_na) complete_cases <- TRUE

	info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order, scale_log2r = scale_log2r, impute_na = impute_na, 
						df = df, filepath = filepath, filename = filename, 
						anal_type = "Trend")(n_clust = n_clust, complete_cases = complete_cases, ...)
}


#'Trend Analysis of Protein \code{log2-ratios}
#'@seealso \code{\link{proteoTrend}} for parameters
#'@export
prnTrend <- function (...) {
	proteoTrend(id = gene, ...)
}

