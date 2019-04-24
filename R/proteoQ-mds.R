#' Make MDS plots
#' 
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotMDS <- function (df_mds, col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
										label_scheme_sub = label_scheme_sub, filepath, filename, show_ids, ...) {
	
	dots <- rlang::exprs(...)
	
	# col_select <- rlang::enexpr(col_select)
	# col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color) # %>% rlang::as_string()
	col_fill <- rlang::enexpr(col_fill) # %>% rlang::as_string()
	col_shape <- rlang::enexpr(col_shape) # %>% rlang::as_string()
	col_size <- rlang::enexpr(col_size) # %>% rlang::as_string()
	col_alpha <- rlang::enexpr(col_alpha) # %>% rlang::as_string()
	
	my_theme <- theme_bw() +
			theme(
				 axis.text.x  = element_text(angle=0, vjust=0.5, size=20),
				 axis.text.y  = element_text(angle=0, vjust=0.5, size=20),   
				 axis.title.x = element_text(colour="black", size=20),
				 axis.title.y = element_text(colour="black", size=20),
				 plot.title = element_text(face="bold", colour="black", size=20, hjust=0.5, vjust=0.5),

				 panel.grid.major.x = element_blank(), 
				 panel.grid.minor.x = element_blank(),
				 panel.grid.major.y = element_blank(),
				 panel.grid.minor.y = element_blank(), 

				 legend.key = element_rect(colour = NA, fill = 'transparent'), 
				 legend.background = element_rect(colour = NA,  fill = "transparent"),
				 legend.title = element_blank(),
				 legend.text = element_text(colour="black", size=14),
				 legend.text.align = 0, 
				 legend.box = NULL
			)

	myPalette <- c(rep("blue", 7), rep("red", 7))

	df_mds <- df_mds %>% 
		tibble::rownames_to_column("Sample_ID") %>% 
		dplyr::left_join(label_scheme_sub) %>% 
		dplyr::rename(Color := !!col_color, Fill := !!col_fill, Shape := !!col_shape, Size := !!col_size, Alpha := !!col_alpha) %>% 
		dplyr::mutate_at(vars(one_of("Color", "Fill", "Shape", "Size", "Alpha")), ~ as.factor(.)) %>% 
		dplyr::select(which(not_all_NA(.)))	

	mapping <- ggplot2::aes(x = Coordinate.1, y = Coordinate.2, 
			colour = Color, fill = Fill, shape = Shape, 
			size = Size, alpha = Alpha, stroke = NA)
	
	nms <- names(df_mds) %>% 
		.[. %in% c("Color", "Fill", "Shape", "Size", "Alpha")] 
	
	aes_nms <- mapping %>% 
		purrr::map(., rlang::quo_name) %>% 
		.[. %in% nms]
	
	mapping_var <- mapping %>% 
		.[names(.) %in% c("x", "y", names(aes_nms))]
	
	mapping_fix <- mapping %>% 
		.[!names(.) %in% c("x", "y", names(aes_nms))] %>% 
		purrr::map(., rlang::quo_name)
	
	rm(mapping, aes_nms, nms)
	
	# nms <- names(df_mds) %>% .[. %in% c("Color", "Shape", "Size", "Alpha")] 
	# mapping = ggplot2::aes(x = Coordinate.1, y = Coordinate.2, colour = Color, fill = Fill, shape = Shape, size = Size, alpha = Alpha)
	# aes_nms <- lapply(mapping, rlang::quo_name) %>% .[. %in% nms]
	# mapping_var <- mapping %>% .[names(.) %in% c("x", "y", names(aes_nms))]
	# mapping_fix <- mapping %>% .[!names(.) %in% c("x", "y", names(aes_nms))] %>% purrr::map(., rlang::quo_name)	
	# rm(mapping, aes_nms, nms)

	# fix_args <- list(color = "red", fill = NA, shape = 21, size = 4, alpha = 0.9) %>% .[names(.) %in% names(mapping_fix)]
	fix_args <- list(color = "white", fill = NA, shape = 21, size = 4, alpha = 0.9, stroke = .02) %>% .[names(.) %in% names(mapping_fix)]

	p <- ggplot() + rlang::eval_tidy(rlang::quo(geom_point(data = df_mds, mapping = mapping_var, !!!fix_args)))

	p <- p + 
		scale_fill_brewer(palette = "Set1") + 
		scale_color_brewer(palette = "Set1") + 
		# scale_size_area() + 
		# scale_size_continuous(range = c(1, 8), guide=FALSE)+
		# scale_alpha_continuous(range = c(.1, .5), guide=FALSE)+
		# stat_ellipse(data=df_mds, mapping=aes(x=Coordinate.1, y=Coordinate.2, group=MDS_Shape), color="darkgray") + 
		# scale_colour_gradientn(colours = rainbow(7), guide=guide_colourbar(reverse = TRUE), limits=c(1,25), breaks=c(1,5,10,25)) + 
		# scale_x_continuous(limits=c(-.6, .8), breaks = c(-.5, 0, .5), labels=c("-0.5", "0", "0.5")) + 
		# stat_smooth(method=lm, se=FALSE, size=.5) + 
		# coord_cartesian(xlim=c(-30, 30), ylim=c(-30, 30)) + 
		labs(title = "", x = expression("Coordinate 1"), y = expression("Coordinate 2")) + 
		# scale_colour_manual(values=myPalette, name="", breaks=c("shift.10", "shift.20", "shift.30", "shift.40", "shift.50"), labels=c("r = 1.1", "r = 1.2", "r = 1.3", "r = 1.4", "r = 1.5")) + 
		# scale_x_continuous(name="Coordinate 1", limits=c(-.1, .1), breaks = c(-.1, 0, .1), labels=c("-0.1", "0", "0.1")) + 
		# scale_y_continuous(name="Coordinate 2", limits=c(-.1, .1), breaks = c(-.1, 0, .1), labels=c("-0.1", "0", "0.1")) + 
		coord_fixed() + 
		# scale_colour_discrete(name="") +
	my_theme
	
	if (show_ids) 
	  p <- p + geom_text(data = df_mds, mapping = aes(x = Coordinate.1, y = Coordinate.2, label = df_mds$Sample_ID), color = "gray", size = 3)

	ggsave(file.path(filepath, filename), p, width = 8, height = 6, dpi = 600, units = "in")

}


#' Make EucDist plots
#' 
#' @import dplyr ggplot2 rlang pheatmap
#' @importFrom magrittr %>%
plotEucDist <- function (D, annot_cols, filepath, filename, ...) {
	
	dots <- rlang::enexprs(...)
	dots[grepl("^mds_", names(dots))] <- NULL

	D_matrix <- as.matrix(D)
	
	# set up color palette
	n_color <- 500
	xmin <- 0
	xmax <- ceiling(max(D_matrix))
	x_margin <- xmax/10
	color_breaks <- c(seq(xmin, x_margin, length = n_color/2)[1 : (n_color/2-1)], seq(x_margin, xmax, length = n_color/2)[2 : (n_color/2)])	

	if (is.null(dots$color)) mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color) else mypalette <- eval(dots$color, env = caller_env())

	# set up annotation columns
	if (is.null(annot_cols)) annotation_col <- NA else 
		annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = attr(D, "Labels")) 
	
	# set up the color palettes of the annotation columns
	if (is.null(dots$annotation_colors)) {
		annotation_colors <- setHMColor(annotation_col)
	} else if (is.na(dots$annotation_colors)) {
		annotation_colors <- NA
	} else {
		annotation_colors <- eval(dots$annotation_colors, env = caller_env())
	}

	# parameter(s) in "dots" being disabled 
	nm_idx <- names(dots) %in% c("mat", "filename", "annotation_col", "color", "annotation_colors", "breaks")
	dots[nm_idx] <- NULL

	n_TMT_sets <- n_TMT_sets(label_scheme)
	max_width <- 77
	
	if(is.null(dots$width)) dots$width <- pmin(10*n_TMT_sets*1.2, max_width)

	if(dots$width >= max_width) {
		cat("The width for the graphic device is", dots$width, "inches or more.\n")
		stop("Please consider a a smaller `cellwidth`.")
	}
	
	my_pheatmap(
		mat = D_matrix, 
		filename = file.path(filepath, filename), 
		annotation_col = annotation_col, 
		# annotation_col = NULL, 
		color = mypalette, 
		annotation_colors = annotation_colors, 
		breaks = color_breaks, 
		!!!dots
	)
	
}


#' Calculate MDS scores
#' 
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scoreMDS <- function (df, label_scheme_sub, scale_log2r, adjEucDist = FALSE, classical, ...) {

	dots <- rlang::exprs(...)

	M <- cor(df, use="pairwise.complete.obs")
	D <- dist(t(df), method = "euclidean", diag = TRUE, upper = TRUE)
	D_matrix <- as.matrix(D)

	if (adjEucDist) { ## adjust Euclidean distance for samples from two different TMT sets
		annotation_col <- colAnnot(annot_cols = c("TMT_Set"), sample_ids = attr(D, "Labels")) 
		
		for (i in 1:ncol(D_matrix)) {
			for (j in 1:ncol(D_matrix)) {
				if(annotation_col$TMT_Set[i] != annotation_col$TMT_Set[j]) D_matrix[i, j] <- D_matrix[i, j]/sqrt(2)
			}
		}
	}

	k <- 3
	if (!classical) {
		df_mds <- isoMDS(D_matrix, k = k)
		df_mds <- data.frame(df_mds$points)
	} else {
		df_mds <- data.frame(cmdscale(D_matrix, k = k))
	}
	
	nms <- lapply(dots, expr_name) %>% 
		.[. %in% names(label_scheme_sub)]
		
	lookup <- nms %>% 
		data.frame(check.names = FALSE)

	label_scheme_mds <- label_scheme_sub[, names(label_scheme_sub) %in% nms] %>% 
		`names<-`(names(lookup)) %>% 
		dplyr::bind_cols(label_scheme_sub[, "Sample_ID", drop = FALSE]) %>% 
		dplyr::select(which(not_all_NA(.)))
	
	dots[names(dots) %in% names(lookup)] <- NULL
	
	df_mds <- df_mds %>% 
		`colnames<-`(paste("Coordinate", 1:k, sep = ".")) %>% 
		tibble::rownames_to_column("Sample_ID") %>% 
		dplyr::left_join(label_scheme_mds, by = "Sample_ID") %>% 
		dplyr::select(which(not_all_zero(.))) %>% 
		tibble::column_to_rownames(var = "Sample_ID")	
	
	# -----------------
	# transposed data for pca by sample IDs
	df_t <- df[complete.cases(df), ] %>% 
		t() %>% 
		data.frame(check.names = FALSE) %>% 
		bind_cols(label_scheme_sub)
	
	pr_out <- df_t %>% 
		dplyr::select(which(colnames(.) %in% rownames(df))) %>% 
		prcomp(scale = scale_log2r)
		
	rownames(pr_out$x) <- label_scheme_sub$Sample_ID
	
	# Proportion of Variance
	prop_var <- summary(pr_out)$importance[2, ] %>% scales::percent()	
	
	df_pca <- pr_out$x %>% 
		data.frame(check.names = FALSE) %>% 
		tibble::rownames_to_column("Sample_ID") %>% 
		dplyr::left_join(label_scheme_mds, by = "Sample_ID") %>% 
		dplyr::select(which(not_all_zero(.))) %>% 
		tibble::column_to_rownames(var = "Sample_ID")	%>% 
		`names<-`(gsub("^PC", "Coordinate\\.", names(.))) 
		
	# pca for biplot
	pr_bi <- df[complete.cases(df), ] %>% prcomp(scale = scale_log2r) 
	
	df_pr_bi <- pr_bi$x %>% 
		data.frame(check.names = FALSE) %>% 
		`names<-`(gsub("^PC", "Coordinate\\.", names(.))) 

	# Proportion of Variance
	prop_var_bi <- summary(pr_bi)$importance[2, ] %>% scales::percent()
	
	# -----------------

	return(list("MDS" = df_mds, "D" = D, "PCA" = df_pca, "prop_var" = prop_var, "pr_bi" = df_pr_bi, "prop_var_bi" = prop_var_bi))
}


#' Make PCA plots
#' 
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotPCA <- function (df, col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
										label_scheme_sub = label_scheme_sub, prop_var, pr_bi, prop_var_bi, filepath, filename, show_ids, ...) { # Plot PCA

	dots <- rlang::enexprs(...)
	
	# col_select <- rlang::enexpr(col_select)
	# col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color) # %>% rlang::as_string()
	col_fill <- rlang::enexpr(col_fill) # %>% rlang::as_string()
	col_shape <- rlang::enexpr(col_shape) # %>% rlang::as_string()
	col_size <- rlang::enexpr(col_size) # %>% rlang::as_string()
	col_alpha <- rlang::enexpr(col_alpha) # %>% rlang::as_string()
	
	
	my_theme <- theme_bw() +
	theme(
				 axis.text.x  = element_text(angle=0, vjust=0.5, size=20),   # rotate x-axis label and define label size
				 # axis.ticks.x  = element_blank(), # x-axis ticks
				 axis.text.y  = element_text(angle=0, vjust=0.5, size=20),   # rotate y-axis label and define label size
				 axis.title.x = element_text(colour="black", size=20),  # x-axis title size
				 axis.title.y = element_text(colour="black", size=20),  # y-axis title size
				 plot.title = element_text(face="bold", colour="black", size=20, hjust=0.5, vjust=0.5),   # main title size

				 panel.grid.major.x = element_blank(), 
				 panel.grid.minor.x = element_blank(),
				 panel.grid.major.y = element_blank(),
				 panel.grid.minor.y = element_blank(), 

				 legend.key = element_rect(colour = NA, fill = 'transparent'), 
				 legend.background = element_rect(colour = NA,  fill = "transparent"),
				 # legend.position = c(0.85, 0.85),  
				 # legend.position = "none",  # to remove all legends
				 legend.title = element_blank(),
				 # legend.title = element_text(colour="black", size=16), 
				 legend.text = element_text(colour="black", size=14),
				 legend.text.align = 0, 
				 legend.box = NULL
	)

	myPalette <- c(rep("blue", 7), rep("red", 7))
	
	df <- df %>% 
		tibble::rownames_to_column("Sample_ID") %>% 
		dplyr::left_join(label_scheme_sub) %>% 
		dplyr::rename(Color := !!col_color, Fill := !!col_fill, Shape := !!col_shape, Size := !!col_size, Alpha := !!col_alpha) %>% 
		dplyr::mutate_at(vars(one_of("Color", "Fill", "Shape", "Size", "Alpha")), ~ as.factor(.)) %>% 
		dplyr::select(which(not_all_NA(.)))
		
	# ---------------
	mapping <- ggplot2::aes(x = Coordinate.1, y = Coordinate.2, 
			colour = Color, fill = Fill, shape = Shape, 
			size = Size, alpha = Alpha, stroke = NA)
	
	nms <- names(df) %>% 
		.[. %in% c("Color", "Fill", "Shape", "Size", "Alpha")] 
	
	aes_nms <- mapping %>% 
		purrr::map(., rlang::quo_name) %>% 
		.[. %in% nms]
	
	mapping_var <- mapping %>% 
		.[names(.) %in% c("x", "y", names(aes_nms))]
	
	mapping_fix <- mapping %>% 
		.[!names(.) %in% c("x", "y", names(aes_nms))] %>% 
		purrr::map(., rlang::quo_name)
	
	rm(mapping, aes_nms, nms)

	fix_args <- list(color = "white", fill = NA, shape = 21, size = 4, alpha = 0.9, stroke = .02) %>% .[names(.) %in% names(mapping_fix)]

	p <- ggplot() + rlang::eval_tidy(rlang::quo(geom_point(data = df, mapping = mapping_var, !!!fix_args)))

	p <- p + 
		scale_fill_brewer(palette = "Set1") + 
		scale_color_brewer(palette = "Set1") + 
		labs(title = "", x = paste0("PC1 (", prop_var[1], ")"), y = paste0("PC2 (", prop_var[2], ")")) + 
		coord_fixed() + 
	my_theme
		
	# ---------------
	if (show_ids) p <- p + geom_text(data = df, mapping = aes(x = Coordinate.1, y = Coordinate.2, label = df$Sample_ID), color = "gray", size = 3)
	ggsave(file.path(filepath, filename), p, width = 8, height = 6, dpi = 600, units = "in")
	
	bi_plot <- FALSE
	if (bi_plot) {
		# bi-plot
		p_bi <- ggplot() + 
			geom_point(data = pr_bi, mapping = aes(x = Coordinate.1, y = Coordinate.2), alpha = .9, size = 1) + 
			scale_fill_brewer(palette = "Set1") + 
			scale_color_brewer(palette = "Set1") + 
			labs(title = "", x = paste0("PC1 (", prop_var_bi[1], ")"), y = paste0("PC2 (", prop_var_bi[2], ")")) + 
			coord_fixed() + 
		my_theme
		
		if (show_ids) p_bi <- p_bi + geom_text(data = pr_bi, mapping = aes(x = Coordinate.1, y = Coordinate.2, label = rownames(pr_bi)), color = "gray", size = 3)
		
		fn_prx <- gsub("\\..*$", "", filename)
		fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename) 
		
		ggsave(file.path(filepath, paste0(fn_prx, "bi_plot", ".", fn_suffix)), p_bi, width = 8, height = 6, dpi = 600, units = "in")
	}
	
}


#'MDS Plots
#'
#'\code{proteoMDS} visualizes the results from multidimensional scaling (MDS)
#'and principal component analysis (PCA).
#'
#'An Euclidean distance matrix of \code{log2-ratios} is returned by
#'\code{\link[base]{dist}}, followed by a metric (\code{\link[stats]{cmdscale}})
#'or non-metric (\code{\link[MASS]{isoMDS}}) MDS. The default is metric MDS with
#'the input dissimilarities being euclidean distances.
#'
#'\code{log2-ratios} are used in PCA (\code{\link[stats]{prcomp}}).
#'
#'@param id Character string to indicate the type of data (a column ID). Peptide
#'  data will be used at column \code{id = "pep_seq"} or \code{"pep_seq_mod"},
#'  and protein data at \code{id = "prot_acc"} or \code{"gene"}.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param adjEucDist Logical; if TRUE, adjusts the inter-plex Euclidean distance
#'  by \eqn{1/sqrt(2)}. The option \code{adjEucDist = TRUE} may be suitable when
#'  \code{reference samples} in each TMT plex undergo approximately the same
#'  sample handling process as the rest of the samples. For instance,
#'  \code{reference samples} were split at the levels of protein lysates.
#'  Typically, \code{adjEucDist = FALSE} if \code{reference samples} were split
#'  near the end of a sample handling process, for instance, at the stages
#'  immediately before or after TMT labeling.
#'@param classical Logical; performs metric MDS at TRUE and non-metric MDS at
#'  FALSE.
#'@param show_ids Logical; if TRUE, shows the sample IDs in MDS and PCA plots.
#'@param annot_cols Column names available in \code{expt_smry.csv}. Values
#'  under the selected columns will be used to color-code sample IDs on the top
#'  of \code{EucDist} plots.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the name of the current function \code{call} and the value of
#'  \code{id}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... More parameters for plotting:
#'
#'  \code{mds_color}, character string to a column ID in \code{Label
#'  scheme.csv}; values under which will be used to color-code data points in
#'  \code{MDS/PCA} plots.
#'
#'  \code{mds_shape}, character string to a column ID in \code{Label
#'  scheme.csv}; values under which will be used to shape-code data points in
#'  \code{MDS/PCA} plots.
#'
#'  \code{mds_alpha}, character string to a column ID in \code{Label
#'  scheme.csv}; values under which will be used to control the transparency of
#'  data points in \code{MDS/PCA} plots.
#'@return The plots of MDS, PCA and Euclidean distance.
#'
#' @examples
#' head(label_scheme)
#'
#' proteoMDS(
#'   id = pep_seq_mod,
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#' )
#'
#' proteoMDS(
#'   id = pep_seq_mod,
#'   scale_log2r = TRUE,
#'
#'   # parameters for MDS
#'   adjEucDist = FALSE,
#'
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "MDS_Group"),
#'
#'   display_numbers = TRUE,
#'   number_color = "grey30",
#'   number_format = "%.2f",
#'
#'   clustering_distance_rows = "euclidean",
#'   clustering_distance_cols = "euclidean",
#'
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 16,
#'
#'   cluster_rows = TRUE,
#'   show_rownames = TRUE,
#'   show_colnames = TRUE,
#'   border_color = "grey60",
#'   cellwidth = 50,
#'   cellheight = 50
#' )
#' 
#' \dontrun{
#' proteoMDS(
#'   id = pep_seq_mod,
#'   col_color = "column_key_not_in_label_scheme", 
#'   col_shape = "also_a_missing_column_key", 
#'   annot_cols = c("bad_column_key", "yet_another_bad_column_key")
#' )
#' }
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoMDS <- function (id = gene, 
											col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL, 
											col_shape = NULL, col_size = NULL, col_alpha = NULL, 
											scale_log2r = FALSE, adjEucDist = FALSE, classical = TRUE, show_ids = TRUE, 
                      annot_cols = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {
	
	id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)

	info_anal(id = !!id, 
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill, 
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
		scale_log2r = scale_log2r, impute_na = FALSE, df = df, filepath = filepath, filename = filename, 
		anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, show_ids = show_ids, annot_cols = annot_cols, 
		...
	)
}


#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoPCA <- function (id = gene, 
											col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL, 
											col_shape = NULL, col_size = NULL, col_alpha = NULL, 
											scale_log2r = FALSE, show_ids = TRUE, 
                      annot_cols = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {
	
	id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)

	info_anal(id = !!id, 
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill, 
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
		scale_log2r = scale_log2r, impute_na = FALSE, df = df, filepath = filepath, filename = filename, 
		anal_type = "PCA")(show_ids = show_ids, annot_cols = annot_cols, 
		...
	)
}


#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoEucDist <- function (id = gene, 
											col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL, 
											col_shape = NULL, col_size = NULL, col_alpha = NULL, 
											scale_log2r = FALSE, adjEucDist = FALSE, 
                      annot_cols = NULL, 
											df = NULL, filepath = NULL, filename = NULL, ...) {
	
	id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)

	info_anal(id = !!id, 
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill, 
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
		scale_log2r = scale_log2r, impute_na = FALSE, df = df, filepath = filepath, filename = filename, 
		anal_type = "EucDist")(adjEucDist = adjEucDist, 
		annot_cols = annot_cols, 
		...
	)
}




#'MDS plots of Protein Data
#'
#' @examples
#' prnMDS(
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   col_color = Color, 
#'   col_shape = Shape
#' )
#'
#' prnMDS(
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   col_color = Color, 
#'   col_shape = Shape, 
#'   
#'   display_numbers = TRUE,
#'   number_color = "grey30",
#'   number_format = "%.2f",
#'   
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 16,
#'   
#'   cluster_rows = TRUE,
#'   show_rownames = TRUE,
#'   show_colnames = TRUE,
#'   border_color = "grey60",
#'   cellwidth = 50,
#'   cellheight = 50
#' )
#'
#' \dontrun{
#' prnMDS(
#'   col_color = "column_key_not_in_label_scheme", 
#'   col_shape = "also_a_missing_column_key", 
#'   annot_cols = c("bad_column_key", "yet_another_bad_column_key")
#' )
#' }
#'
#'@seealso \code{\link{proteoMDS}} for parameters
#'@export
prnMDS <- function (...) {
	proteoMDS(id = gene, ...)
}

#'MDS Plots for Peptide Data
#'
#' @examples
#' pepMDS(
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   col_color = Color, 
#'   col_shape = Shape
#' )
#'
#' pepMDS(
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   col_color = Color, 
#'   col_shape = Shape, 
#'   
#'   display_numbers = TRUE,
#'   number_color = "grey30",
#'   number_format = "%.2f",
#'   
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 16,
#'   
#'   cluster_rows = TRUE,
#'   show_rownames = TRUE,
#'   show_colnames = TRUE,
#'   border_color = "grey60",
#'   cellwidth = 50,
#'   cellheight = 50
#' )
#' 
#' \dontrun{
#' pepMDS(
#'   col_color = "column_key_not_in_label_scheme", 
#'   col_shape = "also_a_missing_column_key", 
#'   annot_cols = c("bad_column_key", "yet_another_bad_column_key")
#' )
#' }
#'
#'@seealso \code{\link{proteoMDS}} for parameters
#'@export
pepMDS <- function (...) {
	proteoMDS(id = pep_seq, ...)
}


#'@export
prnPCA <- function (...) {
	proteoPCA(id = gene, ...)
}


#'@export
pepPCA <- function (...) {
	proteoPCA(id = pep_seq, ...)
}


#'@export
prnEucDist <- function (...) {
	proteoEucDist(id = gene, ...)
}


#'@export
pepEucDist <- function (...) {
	proteoEucDist(id = pep_seq, ...)
}

