#' Plots MDS
#'
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotMDS <- function (df, col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL,
                     col_alpha = NULL, label_scheme_sub = label_scheme_sub, filepath, filename,
                     show_ids, ...) {

  dots <- rlang::exprs(...)

	col_fill <- rlang::enexpr(col_fill)
	col_color <- rlang::enexpr(col_color)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	
	map_color <- map_fill <- map_shape <- map_size <- map_alpha <- NA

	if (col_color != rlang::expr(Color) | !rlang::as_string(sym(col_color)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_color))), "X")
	if (col_fill != rlang::expr(Fill)  | !rlang::as_string(sym(col_fill)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_fill))), "X")
	if (col_shape != rlang::expr(Shape) | !rlang::as_string(sym(col_shape)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_shape))), "X")
	if (col_size != rlang::expr(Size) | !rlang::as_string(sym(col_size)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_size))), "X")
	if (col_alpha != rlang::expr(Alpha) | !rlang::as_string(sym(col_alpha)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_alpha))), "X")
	
	if (!is.na(map_color)) col_color <- NULL
	if (!is.na(map_fill)) col_fill <- NULL
	if (!is.na(map_shape)) col_shape <- NULL
	if (!is.na(map_size)) col_size <- NULL
	if (!is.na(map_alpha)) col_alpha <- NULL
	
	rm(map_color, map_fill, map_shape, map_size, map_alpha)

	my_theme <- theme_bw() + theme(
		 axis.text.x  = element_text(angle=0, vjust=0.5, size=16),
		 axis.text.y  = element_text(angle=0, vjust=0.5, size=16),
		 axis.title.x = element_text(colour="black", size=18),
		 axis.title.y = element_text(colour="black", size=18),
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
	
	mapping <- ggplot2::aes(x = Coordinate.1, y = Coordinate.2,
	                        colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
	                        size = !!col_size, alpha = !!col_alpha)

	idx <- purrr::map(mapping, `[[`, 1) %>% 
	  purrr::map_lgl(is.null)
	
	mapping_var <- mapping[!idx]
	mapping_fix <- mapping[idx]

	fix_args <- list(colour = "darkgray", fill = NA, shape = 21, size = 4, alpha = 0.9) %>% 
	  .[names(.) %in% names(mapping_fix)] %>% 
	  .[!is.na(.)]
	fix_args$stroke <- 0.02

	# aes_qns <- mapping %>% purrr::map(rlang::quo_name) %>% purrr::flatten_chr(.) %>% .[. %in% names(df)]
	p <- ggplot() + rlang::eval_tidy(rlang::quo(geom_point(data = df, mapping = mapping_var, !!!fix_args)))

	p <- p +
		# scale_fill_brewer(palette = "Set1") +
		# scale_color_brewer(palette = "Set1") +
		labs(title = "", x = expression("Coordinate 1"), y = expression("Coordinate 2")) +
		coord_fixed() +
	my_theme

	if (show_ids)
	  p <- p +
	    geom_text(data = df, mapping = aes(x = Coordinate.1, y = Coordinate.2,
	                                           label = df$Sample_ID), color = "gray", size = 3)

	ggsave(file.path(filepath, filename), p, ...)
}


#' Plots EucDist
#'
#' @import dplyr ggplot2 rlang pheatmap
#' @importFrom magrittr %>%
plotEucDist <- function (D, annot_cols, filepath, filename, annot_colnames, ...) {
	dots <- rlang::enexprs(...)

	D_matrix <- as.matrix(D)

	n_color <- 500
	xmin <- 0
	xmax <- ceiling(max(D_matrix))
	x_margin <- xmax/10
	color_breaks <- c(seq(xmin, x_margin, length = n_color/2)[1 : (n_color/2-1)],
	                  seq(x_margin, xmax, length = n_color/2)[2 : (n_color/2)])

	if (is.null(dots$color)) {
	  mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
	} else {
	  mypalette <- eval(dots$color, env = caller_env())
	}

	if (is.null(annot_cols)) annotation_col <- NA else
		annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = attr(D, "Labels"))
	
	if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
	  colnames(annotation_col) <- annot_colnames
	}

	if (is.null(dots$annotation_colors)) {
		annotation_colors <- setHMColor(annotation_col)
	} else if (is.na(dots$annotation_colors)) {
		annotation_colors <- NA
	} else {
		annotation_colors <- eval(dots$annotation_colors, env = caller_env())
	}

	nm_idx <- names(dots) %in% c("mat", "filename", "annotation_col", "color",
	                             "annotation_colors", "breaks")
	dots[nm_idx] <- NULL

	load(file = file.path(dat_dir, "label_scheme.Rdata"))
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


#' Plots PCA
#'
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotPCA <- function (df, col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL,
                     col_alpha = NULL, label_scheme_sub = label_scheme_sub, prop_var, pr_bi,
                     prop_var_bi, filepath, filename, show_ids, ...) {

	dots <- rlang::enexprs(...)

	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	
	map_color <- map_fill <- map_shape <- map_size <- map_alpha <- NA
	
	if (col_color != rlang::expr(Color) | !rlang::as_string(sym(col_color)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_color))), "X")
	if (col_fill != rlang::expr(Fill)  | !rlang::as_string(sym(col_fill)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_fill))), "X")
	if (col_shape != rlang::expr(Shape) | !rlang::as_string(sym(col_shape)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_shape))), "X")
	if (col_size != rlang::expr(Size) | !rlang::as_string(sym(col_size)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_size))), "X")
	if (col_alpha != rlang::expr(Alpha) | !rlang::as_string(sym(col_alpha)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_alpha))), "X")
	
	if (!is.na(map_color)) col_color <- NULL
	if (!is.na(map_fill)) col_fill <- NULL
	if (!is.na(map_shape)) col_shape <- NULL
	if (!is.na(map_size)) col_size <- NULL
	if (!is.na(map_alpha)) col_alpha <- NULL
	
	rm(map_color, map_fill, map_shape, map_size, map_alpha)
	
	my_theme <- theme_bw() + theme(
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

	mapping <- ggplot2::aes(x = Coordinate.1, y = Coordinate.2,
	                        colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
	                        size = !!col_size, alpha = !!col_alpha)

	idx <- purrr::map(mapping, `[[`, 1) %>% 
	  purrr::map_lgl(is.null)
	
	mapping_var <- mapping[!idx]
	mapping_fix <- mapping[idx]
	
	fix_args <- list(colour = "darkgray", fill = NA, shape = 21, size = 4, alpha = 0.9) %>% 
	  .[names(.) %in% names(mapping_fix)] %>% 
	  .[!is.na(.)]
	fix_args$stroke <- 0.02

	p <- ggplot() +
	  rlang::eval_tidy(rlang::quo(geom_point(data = df, mapping = mapping_var, !!!fix_args)))

	p <- p +
		# scale_fill_brewer(palette = "Set1") +
		# scale_color_brewer(palette = "Set1") +
		labs(title = "", x = paste0("PC1 (", prop_var[1], ")"), y = paste0("PC2 (", prop_var[2], ")")) +
		coord_fixed() +
	my_theme

	if (show_ids)
	  p <- p +
	    geom_text(data = df,
	              mapping = aes(x = Coordinate.1, y = Coordinate.2, label = df$Sample_ID),
	              color = "gray", size = 3)

	ggsave(file.path(filepath, filename), p, ...)

	bi_plot <- FALSE
	if (bi_plot) {
		p_bi <- ggplot() +
			geom_point(data = pr_bi, mapping = aes(x = Coordinate.1, y = Coordinate.2),
			           alpha = .9, size = 1) +
			scale_fill_brewer(palette = "Set1") +
			scale_color_brewer(palette = "Set1") +
			labs(title = "", x = paste0("PC1 (", prop_var_bi[1], ")"),
			     y = paste0("PC2 (", prop_var_bi[2], ")")) +
			coord_fixed() +
		my_theme

		if (show_ids)
		  p_bi <- p_bi + geom_text(data = pr_bi,
		                           mapping = aes(x = Coordinate.1, y = Coordinate.2,
		                                         label = rownames(pr_bi)), color = "gray", size = 3)

		fn_prx <- gsub("\\..*$", "", filename)
		fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)

		ggsave(file.path(filepath, paste0(fn_prx, "bi_plot", ".", fn_suffix)), p_bi, ...)
	}
}


#' Scores MDS and PCA
#'
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scoreMDS <- function (df, label_scheme_sub, scale_log2r, adjEucDist = FALSE, classical, ...) {

	dots <- rlang::exprs(...)

	M <- cor(df, use="pairwise.complete.obs")
	D <- dist(t(df), method = "euclidean", diag = TRUE, upper = TRUE)
	D_matrix <- as.matrix(D)

	# adjust Euclidean distance for samples from two different TMT sets
	if (adjEucDist) {
		annotation_col <- colAnnot(annot_cols = c("TMT_Set"), sample_ids = attr(D, "Labels"))

		for (i in 1:ncol(D_matrix)) {
			for (j in 1:ncol(D_matrix)) {
				if (annotation_col$TMT_Set[i] != annotation_col$TMT_Set[j])
				  D_matrix[i, j] <- D_matrix[i, j]/sqrt(2)
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

	# Proportion of variance
	prop_var_bi <- summary(pr_bi)$importance[2, ] %>% scales::percent()

	return(list("MDS" = df_mds, "D" = D, "PCA" = df_pca, "prop_var" = prop_var,
	            "pr_bi" = df_pr_bi, "prop_var_bi" = prop_var_bi))
}


#'Visualizes MDS plots
#'
#'\code{proteoMDS} visualizes the results from multidimensional scaling (MDS).
#'
#'An Euclidean distance matrix of \code{log2-ratios} is returned by
#'\code{\link[base]{dist}}, followed by a metric (\code{\link[stats]{cmdscale}})
#'or non-metric (\code{\link[MASS]{isoMDS}}) MDS. The default is metric MDS with
#'the input dissimilarities being euclidean distances.
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at column \code{id = "pep_seq"} or \code{"pep_seq_mod"}, and protein
#'  data at \code{id = "prot_acc"} or \code{"gene"}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param  col_group Not currently used.
#'@param  col_color Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Color}. Values under which will be used for the
#'  \code{color} aesthetics in \code{ggplot2}.
#'@param  col_fill Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Fill}. Values under which will be used for the
#'  \code{fill} aesthetics in \code{ggplot2}.
#'@param  col_shape Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Shape}. Values under which will be used for the
#'  \code{shape} aesthetics in \code{ggplot2}.
#'@param  col_size Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Size}. Values under which will be used for the
#'  \code{size} aesthetics in \code{ggplot2}.
#'@param  col_alpha Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Alpha}. Values under which will be used for the
#'  \code{alpha} aesthetics in \code{ggplot2}.
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
#'@param show_ids Logical; if TRUE, shows the sample IDs in \code{MDS/PCA}
#'  plots.
#'@param annot_cols Not currently used.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the name of the current function \code{call} and the value of
#'  \code{id}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... Parameters inherited from \code{ggsave}: \cr \code{width}, the
#'  width of plot; \cr \code{height}, the height of plot \cr \code{...}
#'
#'@return MDS plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoMDS <- function (id = gene,
											col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL,
											col_shape = NULL, col_size = NULL, col_alpha = NULL,
											scale_log2r = TRUE, adjEucDist = FALSE, classical = TRUE, show_ids = TRUE,
                      annot_cols = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

	id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	
	reload_expts()

	info_anal(id = !!id,
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill,
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
		scale_log2r = scale_log2r, impute_na = FALSE, df = df, filepath = filepath, filename = filename,
		anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, show_ids = show_ids,
		                   annot_cols = annot_cols, ...)
}


#'Visualizes PCA plots
#'
#'\code{protePCA} visualizes the results from principal component analysis
#'(PCA).
#'
#'\code{log2-ratios} are used in PCA (\code{\link[stats]{prcomp}}).
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at column \code{id = "pep_seq"} or \code{"pep_seq_mod"}, and protein
#'  data at \code{id = "prot_acc"} or \code{"gene"}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param  col_group Not currently used.
#'@param  col_color Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Color}. Values under which will be used for the
#'  \code{color} aesthetics in \code{ggplot2}.
#'@param  col_fill Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Fill}. Values under which will be used for the
#'  \code{fill} aesthetics in \code{ggplot2}.
#'@param  col_shape Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Shape}. Values under which will be used for the
#'  \code{shape} aesthetics in \code{ggplot2}.
#'@param  col_size Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Size}. Values under which will be used for the
#'  \code{size} aesthetics in \code{ggplot2}.
#'@param  col_alpha Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Alpha}. Values under which will be used for the
#'  \code{alpha} aesthetics in \code{ggplot2}.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param show_ids Logical; if TRUE, shows the sample IDs in \code{MDS/PCA}
#'  plots.
#'@param annot_cols Not currently used.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the name of the current function \code{call} and the value of
#'  \code{id}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... Parameters inherited from \code{ggsave}: \cr \code{width}, the
#'  width of plot; \cr \code{height}, the height of plot \cr \code{...}
#'
#'@return PCA plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoPCA <- function (id = gene,
											col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL,
											col_shape = NULL, col_size = NULL, col_alpha = NULL,
											scale_log2r = TRUE, show_ids = TRUE,
                      annot_cols = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

  id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	
	reload_expts()

	info_anal(id = !!id,
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill,
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
		scale_log2r = scale_log2r, impute_na = FALSE, df = df, filepath = filepath, filename = filename,
		anal_type = "PCA")(show_ids = show_ids, annot_cols = annot_cols, ...)
}


#'Visualizes Euclidean distances
#'
#'\code{proteoEucDist} visualizes the Euclidean distances of protein or peptide
#'data.
#'
#'An Euclidean distance matrix of \code{log2-ratios} is returned by
#'\code{\link[base]{dist}}.
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at column \code{id = "pep_seq"} or \code{"pep_seq_mod"}, and protein
#'  data at \code{id = "prot_acc"} or \code{"gene"}.
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
#'@param annot_cols Column names available in \code{expt_smry.xlsx}. Values
#'  under the selected columns will be used to color-code sample IDs on the top
#'  of \code{EucDist} plots.
#'@param annot_colnames Character vector of replacement name(s) to
#'  \code{annot_cols}.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the name of the current function \code{call} and the value of
#'  \code{id}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... Parameters inherited from \code{\link[pheatmap]{pheatmap}}
#'
#'@return Plots of distance matrices.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoEucDist <- function (id = gene,
											col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL,
											col_shape = NULL, col_size = NULL, col_alpha = NULL,
											scale_log2r = TRUE, adjEucDist = FALSE,
                      annot_cols = NULL, annot_colnames = NULL, 
											df = NULL, filepath = NULL, filename = NULL, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

  id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	
	reload_expts()
	
	# if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
	#   load(file = file.path(dat_dir, "label_scheme.Rdata"))
	#   
	#   label_scheme <- label_scheme %>% 
	#     dplyr::select(annot_cols) %>% 
	#     `colnames<-`(annot_colnames) %>% 
	#     dplyr::select(-which(names(.) %in% annot_cols)) %>% 
	#     dplyr::bind_cols(label_scheme, .)
	#   
	#    save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))
	# }

	info_anal(id = !!id,
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill,
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
		scale_log2r = scale_log2r, impute_na = FALSE, df = df, filepath = filepath, filename = filename,
		anal_type = "EucDist")(adjEucDist = adjEucDist,
		annot_cols = annot_cols, annot_colnames = annot_colnames, ...
	)
}


#'MDS plots
#'
#'MDS plots of protein data
#'
#' @examples
#' prnMDS(
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#' )
#'
#' \dontrun{
#' prnMDS(
#'   col_color = "column_key_not_existed",
#'   col_shape = "also_a_missing_column_key"
#' )
#' }
#'
#'@seealso \code{\link{proteoMDS}} for parameters
#'@export
prnMDS <- function (...) {
	proteoMDS(id = gene, ...)
}


#'MDS plots
#'
#'MDS plots of peptide data
#'
#' @examples
#' pepMDS(
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#' )
#'
#' \dontrun{
#' pepMDS(
#'   col_color = "column_key_not_existed",
#'   col_shape = "also_a_missing_column_key"
#' )
#' }
#'
#'@seealso \code{\link{proteoMDS}} for parameters
#'@export
pepMDS <- function (...) {
	proteoMDS(id = pep_seq, ...)
}


#'PCA plots
#'
#'PCA plots of protein data
#'
#' @examples
#' prnPCA(
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#' )
#'
#' \dontrun{
#' prnPCA(
#'   col_color = "column_key_not_existed",
#'   col_shape = "also_a_missing_column_key"
#' )
#' }
#'
#'@seealso \code{\link{proteoPCA}} for parameters
#'@export
prnPCA <- function (...) {
	proteoPCA(id = gene, ...)
}


#'PCA plots
#'
#'PCA plots of peptide data
#'
#' @examples
#' pepPCA(
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#' )
#'
#' \dontrun{
#' pepPCA(
#'   col_color = "column_key_not_existed",
#'   col_shape = "also_a_missing_column_key"
#' )
#' }
#'
#'@seealso \code{\link{proteoPCA}} for parameters
#'@export
pepPCA <- function (...) {
	proteoPCA(id = pep_seq, ...)
}


#'Distance plots
#'
#'Euclidean distance of protein data
#'
#' @examples
#' prnEucDist(
#'   scale_log2r = TRUE,
#'   adjEucDist = FALSE,
#'   show_ids = FALSE,
#'   annot_cols = c("Peptide_Yield", "Group"),
#'   display_numbers = TRUE,
#'   number_color = "grey30",
#'   number_format = "%.2f",
#'   clustering_distance_rows = "euclidean",
#'   clustering_distance_cols = "euclidean",
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 8,
#'   cluster_rows = TRUE,
#'   show_rownames = TRUE,
#'   show_colnames = TRUE,
#'   border_color = "grey60",
#'   cellwidth = 24,
#'   cellheight = 24
#' )
#'
#'
#' \dontrun{
#' prnEucDist(
#'   col_color = "column_key_not_existed",
#'   col_shape = "also_a_missing_column_key",
#'   annot_cols = c("bad_column_key", "yet_another_bad_column_key")
#' )
#' }
#'
#'@seealso \code{\link{proteoEucDist}} for parameters
#'@export
prnEucDist <- function (...) {
	proteoEucDist(id = gene, ...)
}


#'Distance plots
#'
#'Euclidean distance of peptide data
#'
#' @examples
#' pepEucDist(
#'   scale_log2r = TRUE,
#'   adjEucDist = FALSE,
#'   show_ids = FALSE,
#'   annot_cols = c("Peptide_Yield", "Group"),
#'   display_numbers = TRUE,
#'   number_color = "grey30",
#'   number_format = "%.2f",
#'   clustering_distance_rows = "euclidean",
#'   clustering_distance_cols = "euclidean",
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 8,
#'   cluster_rows = TRUE,
#'   show_rownames = TRUE,
#'   show_colnames = TRUE,
#'   border_color = "grey60",
#'   cellwidth = 24,
#'   cellheight = 24
#' )
#'
#'
#' \dontrun{
#' pepEucDist(
#'   col_color = "column_key_not_existed",
#'   col_shape = "also_a_missing_column_key",
#'   annot_cols = c("bad_column_key", "yet_another_bad_column_key")
#' )
#' }
#'
#'@seealso \code{\link{proteoEucDist}} for parameters
#'@export
pepEucDist <- function (...) {
	proteoEucDist(id = pep_seq, ...)
}
