#' Plots MDS
#'
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotMDS <- function (df, col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                     color_brewer = NULL, fill_brewer = NULL, 
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                     label_scheme_sub = label_scheme_sub, filepath, filename,
                     show_ids, ...) {

  dots <- rlang::enexprs(...)

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
	
	color_brewer <- rlang::enexpr(color_brewer)
	fill_brewer <- rlang::enexpr(fill_brewer)
	if (!is.null(color_brewer)) color_brewer <- rlang::as_string(color_brewer)
	if (!is.null(fill_brewer)) fill_brewer <- rlang::as_string(fill_brewer)
	
	size_manual <- eval_bare(size_manual, env = caller_env())
	shape_manual <- eval_bare(shape_manual, env = caller_env())
	alpha_manual <- eval_bare(alpha_manual, env = caller_env())
	
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

	p <- ggplot() + rlang::eval_tidy(rlang::quo(geom_point(data = df, mapping = mapping_var, !!!fix_args)))

	if (!is.null(fill_brewer)) p <- p + scale_color_brewer(palette = fill_brewer)
	if (!is.null(color_brewer)) p <- p + scale_color_brewer(palette = color_brewer)
	
	if ((!is.null(col_size)) & (!is.null(size_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_size]])) == length(size_manual))
	  p <- p + scale_size_manual(values = size_manual)
	}
	
	if ((!is.null(col_shape)) & (!is.null(shape_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_shape]])) == length(shape_manual))
	  p <- p + scale_shape_manual(values = shape_manual)
	}
	
	if ((!is.null(col_alpha)) & (!is.null(alpha_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_alpha]])) == length(alpha_manual))
	  p <- p + scale_shape_manual(values = alpha_manual)
	}
	
	p <- p +
		labs(title = "", x = expression("Coordinate 1"), y = expression("Coordinate 2")) +
		coord_fixed() +
	my_theme

	if (show_ids) {
	  p <- p +
	    geom_text(data = df, 
	              mapping = aes(x = Coordinate.1, y = Coordinate.2, 
	                            label = df$Sample_ID), color = "gray", size = 3)	  
	}

	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), dots)
	do.call(ggsave, gg_args)
}


#' Plots EucDist
#'
#' @import dplyr ggplot2 rlang pheatmap
#' @importFrom magrittr %>%
plotEucDist <- function (D, annot_cols, annot_colnames, filepath, filename, ...) {
	dots <- rlang::enexprs(...)

	stopifnot(is.matrix(D))

	n_color <- 500
	xmin <- 0
	xmax <- ceiling(max(as.vector(D)))
	x_margin <- xmax/10
	color_breaks <- c(seq(xmin, x_margin, length = n_color/2)[1 : (n_color/2-1)],
	                  seq(x_margin, xmax, length = n_color/2)[2 : (n_color/2)])

	if (is.null(dots$color)) {
	  mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
	} else {
	  mypalette <- eval(dots$color, env = caller_env())
	}

	if (is.null(annot_cols)) {
	  annotation_col <- NA
	} else {
	  annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = rownames(D))
	  idx <- which(annot_cols %in% colnames(annotation_col))
	  
	  annot_cols <- annot_cols[idx]
	  annot_colnames <- annot_colnames[idx]
	}

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

	load(file = file.path(dat_dir, "label_scheme.rda"))
	n_TMT_sets <- n_TMT_sets(label_scheme)
	max_width <- 77

	if (is.null(dots$width)) dots$width <- pmin(10*n_TMT_sets*1.2, max_width)

	if (dots$width >= max_width) {
		cat("The width for the graphic device is", dots$width, "inches or more.\n")
		stop("Please consider a a smaller `cellwidth`.")
	}

	filename <- gg_imgname(filename)
	
	my_pheatmap(
		mat = D,
		filename = file.path(filepath, filename),
		annotation_col = annotation_col,
		annotation_row = NA, 
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
plotPCA <- function (df, col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                     color_brewer = NULL, fill_brewer = NULL, 
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                     label_scheme_sub = label_scheme_sub, prop_var, 
                     filepath, filename, show_ids, ...) {

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
	
	color_brewer <- rlang::enexpr(color_brewer)
	fill_brewer <- rlang::enexpr(fill_brewer)
	if (!is.null(color_brewer)) color_brewer <- rlang::as_string(color_brewer)
	if (!is.null(fill_brewer)) fill_brewer <- rlang::as_string(fill_brewer)
	
	size_manual <- eval_bare(size_manual, env = caller_env())
	shape_manual <- eval_bare(shape_manual, env = caller_env())
	alpha_manual <- eval_bare(alpha_manual, env = caller_env())
	
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

	if (!is.null(fill_brewer)) p <- p + scale_color_brewer(palette = fill_brewer)
	if (!is.null(color_brewer)) p <- p + scale_color_brewer(palette = color_brewer)
	
	if ((!is.null(col_size)) & (!is.null(size_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_size]])) == length(size_manual))
	  p <- p + scale_size_manual(values = size_manual)
	}
	
	if ((!is.null(col_shape)) & (!is.null(shape_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_shape]])) == length(shape_manual))
	  p <- p + scale_shape_manual(values = shape_manual)
	}
	
	if ((!is.null(col_alpha)) & (!is.null(alpha_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_alpha]])) == length(alpha_manual))
	  p <- p + scale_shape_manual(values = alpha_manual)
	}
	
	p <- p +
		labs(title = "", x = paste0("PC1 (", prop_var[1], ")"), y = paste0("PC2 (", prop_var[2], ")")) +
		coord_fixed() +
	my_theme

	if (show_ids)
	  p <- p +
	    geom_text(data = df,
	              mapping = aes(x = Coordinate.1, y = Coordinate.2, label = df$Sample_ID),
	              color = "gray", size = 3)

	filename <- gg_imgname(filename)
	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), dots)
	do.call(ggsave, gg_args)
	# ggsave(file.path(filepath, filename), p, ...)
}


#' Scores MDS
#'
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scoreMDS <- function (df, id, label_scheme_sub, anal_type, scale_log2r, 
                      adjEucDist = FALSE, classical, k = 3, ...) {

	dots <- rlang::enexprs(...)
	id <- rlang::as_string(rlang::enexpr(id))
	
	stopifnot(nrow(df) > 50)

	df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
	             sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
	  .$log2R

	D <- dist(t(df), method = "euclidean", diag = TRUE, upper = TRUE)
	if (anyNA(D)) stop("Distance cannot be calculated for one more sample pairs.")
	D <- as.matrix(D)

	if (adjEucDist) {
	  D <- local({
  	  annotation_col <- colAnnot(annot_cols = c("TMT_Set"), sample_ids = attr(D, "Labels"))
  
  		for (i in 1:ncol(D)) {
  			for (j in 1:ncol(D)) {
  				if (annotation_col$TMT_Set[i] != annotation_col$TMT_Set[j])
  				  D[i, j] <- D[i, j]/sqrt(2)
  			}
  		}
  		
  	  return(D)
		})
	}

	if (!classical) {
		df_mds <- data.frame(isoMDS(D, k = k)$points)
	} else {
		df_mds <- data.frame(cmdscale(D, k = k))
	}
	
	df_mds <- df_mds %>%
		`colnames<-`(paste("Coordinate", 1:k, sep = ".")) %>%
		tibble::rownames_to_column("Sample_ID") %>%
		dplyr::select(which(not_all_zero(.))) %>%
		tibble::column_to_rownames(var = "Sample_ID")
}


#' Scores PCA
#'
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scorePCA <- function (df, id, label_scheme_sub, anal_type, scale_log2r, type, ...) {
  
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  
  stopifnot(nrow(df) > 50)
  
  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R 

  df <- df[complete.cases(df), ]
  
  if (type == "obs") {
    df_t <- df %>% 
      t() %>%
      data.frame(check.names = FALSE) %>%
      bind_cols(label_scheme_sub)
    
    pr_out <- df_t %>%
      dplyr::select(which(colnames(.) %in% rownames(df))) %>%
      prcomp(scale = scale_log2r)
    
    rownames(pr_out$x) <- label_scheme_sub$Sample_ID
    
    prop_var <- summary(pr_out)$importance[2, ] %>% scales::percent()
    
    df_pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      `names<-`(gsub("^PC", "Coordinate\\.", names(.)))    
  } else if (type == "feats") {
    pr_out <- df %>% prcomp(scale = scale_log2r)
    
    df_pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      `names<-`(gsub("^PC", "Coordinate\\.", names(.)))
    
    prop_var <- summary(pr_out)$importance[2, ] %>% scales::percent()    
  } else {
    stop("Unkown `type` for PCA.", call. = FALSE)
  }

  return(list("PCA" = df_pca, "prop_var" = prop_var))
}


#' Scores Euclidean distance
#'
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scoreEucDist <- function (df, id, label_scheme_sub, anal_type, scale_log2r, adjEucDist = FALSE, ...) {

  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  
  stopifnot(nrow(df) > 50)
  
  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R
  
  D <- dist(t(df), method = "euclidean", diag = TRUE, upper = TRUE)
  if (anyNA(D)) stop("Distance cannot be calculated for one more sample pairs.")
  
  D <- as.matrix(D)
  
  if (adjEucDist) {
    D <- local({
      annotation_col <- colAnnot(annot_cols = c("TMT_Set"), colnames(D))
      
      for (i in 1:ncol(D)) {
        for (j in 1:ncol(D)) {
          if (annotation_col$TMT_Set[i] != annotation_col$TMT_Set[j])
            D[i, j] <- D[i, j]/sqrt(2)
        }
      }
      
      return(D)
    })
  }
  
  return(D)
}


#'Visualization of MDS plots
#'
#'\code{proteoMDS} visualizes the results from multidimensional scaling (MDS).
#'Users should avoid call the method directly, but instead use the following
#'wrappers.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}}, followed by a metric
#'(\code{\link[stats]{cmdscale}}) or non-metric (\code{\link[MASS]{isoMDS}})
#'MDS. The default is metric MDS with the input dissimilarities being euclidean
#'distances.
#'
#'The function matches the current \code{id} to the grouping argument in the
#'latest \code{call} to \code{\link{normPSM}}. See also \code{\link{prnHist}}
#'for details.
#'
#'@inheritParams  proteoHist
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
#'   i.e., \code{color_brewer = Set1}.
#'@param fill_brewer Character string to the name of a color brewer for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_brewer.html}{ggplot2::scale_fill_brewer},
#'   i.e., \code{fill_brewer = Spectral}.
#'@param size_manual Numeric vector to the scale of sizes for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_size_manual},
#'   i.e., \code{size_manual = c(8, 12)}.
#'@param shape_manual Numeric vector to the scale of shape IDs for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_shape_manual},
#'   i.e., \code{shape_manual = c(5, 15)}. .
#'@param alpha_manual Numeric vector to the scale of transparency of objects for
#'  use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_alpha_manual}
#'   , i.e., \code{alpha_manual = c(.5, .9)}.
#'@param adjEucDist Logical; if TRUE, adjusts the inter-plex Euclidean distance
#'  by \eqn{1/sqrt(2)}. The option \code{adjEucDist = TRUE} may be suitable when
#'  \code{reference samples} from each TMT plex undergo approximately the same
#'  sample handling process as the samples of interest. For instance,
#'  \code{reference samples} were split at the levels of protein lysates.
#'  Typically, \code{adjEucDist = FALSE} if \code{reference samples} were split
#'  near the end of a sample handling process, for instance, at the stages
#'  immediately before or after TMT labeling.
#'@param classical Logical; performs metric MDS at TRUE and non-metric MDS at
#'  FALSE (see also \code{\link[stats]{cmdscale}} and
#'  \code{\link[MASS]{isoMDS}}).
#'@param k Numeric; The desired dimension for the solution passed to
#'  \code{\link[stats]{cmdscale}}.
#'@param show_ids Logical; if TRUE, shows the sample IDs in \code{MDS/PCA}
#'  plots.
#'@param annot_cols Not used.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr Additional parameters for
#'  \code{ggsave}: \cr \code{width}, the width of plot; \cr \code{height}, the
#'  height of plot \cr \code{...}
#'
#'@return MDS plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoMDS <- function (id = gene,
											col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL,
											col_shape = NULL, col_size = NULL, col_alpha = NULL,
											color_brewer = NULL, fill_brewer = NULL, 
											size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
											scale_log2r = TRUE, adjEucDist = FALSE, classical = TRUE, k = 3, 
											show_ids = TRUE, annot_cols = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

	id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
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
	
	reload_expts()

	info_anal(id = !!id,
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill,
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
		color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
		size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
		scale_log2r = scale_log2r, impute_na = FALSE, df = !!df, filepath = !!filepath, filename = !!filename,
		anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, k = k, show_ids = show_ids,
		                   annot_cols = annot_cols, ...)
}


#'Visualization of PCA plots
#'
#'\code{proteoPCA} visualizes the results from principal component analysis
#'(PCA). Users should avoid call the method directly, but instead use the
#'following wrappers.
#'
#'\code{log2FC} are used in PCA (\code{\link[stats]{prcomp}}).
#'
#'The function matches the current \code{id} to the grouping argument in the
#'latest \code{call} to \code{\link{normPSM}}. See also \code{\link{prnHist}}
#'for details.
#'
#'@inheritParams proteoMDS
#'@param type Character string indicating the type of PCA. At the \code{type =
#'  obs} default, the components are by observations; at \code{type = feats},
#'  the components are by features.
#'
#'@return PCA plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoPCA <- function (id = gene, type = "obs", 
                       col_select = NULL, col_group = NULL, col_color = NULL, 
                       col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                       color_brewer = NULL, fill_brewer = NULL, 
                       size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                       scale_log2r = TRUE, show_ids = TRUE, annot_cols = NULL, 
                       df = NULL, filepath = NULL, filename = NULL, ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
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
	type <- rlang::enexpr(type)	
	
	reload_expts()

	info_anal(id = !!id,
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill,
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
		color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
		size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
		scale_log2r = scale_log2r, impute_na = FALSE, df = !!df, filepath = !!filepath, filename = !!filename,
		anal_type = "PCA")(type = !!type, show_ids = show_ids, annot_cols = annot_cols, ...)
}


#'Visualization of the Euclidean distance matrix
#'
#'\code{proteoEucDist} visualizes the heat map of Euclidean distances. Users
#'should avoid call the method directly, but instead use the following wrappers.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}} for heat map visualization.
#'
#'The function matches the current \code{id} to the grouping argument in the
#'latest \code{call} to \code{\link{normPSM}}. See also \code{\link{prnHist}}
#'for details.
#'
#'@inheritParams proteoMDS
#'@param annot_cols A character vector of column keys in \code{expt_smry.xlsx}.
#'  The values under the selected keys will be used to color-code sample IDs on
#'  the top of the indicated plot.
#'@param annot_colnames A character vector of replacement name(s) to
#'  \code{annot_cols}.
#'@param ... Parameters for \code{\link[pheatmap]{pheatmap}}
#'
#'@return Heat map visualization of distance matrices.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoEucDist <- function (id = gene, col_select = NULL, col_group = NULL, col_color = NULL, col_fill = NULL, 
                           col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                           scale_log2r = TRUE, adjEucDist = FALSE, 
                           annot_cols = NULL, annot_colnames = NULL, 
                           df = NULL, filepath = NULL, filename = NULL, ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()
	
	info_anal(id = !!id,
		col_select = !!col_select, col_group = !!col_group, col_color = !!col_color, col_fill = !!col_fill,
		col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
		scale_log2r = scale_log2r, impute_na = FALSE, df = !!df, filepath = !!filepath, filename = !!filename,
		anal_type = "EucDist")(adjEucDist = adjEucDist,
		annot_cols = annot_cols, annot_colnames = annot_colnames, ...)
}


#'MDS plots
#'
#'\code{pepMDS} is a wrapper of \code{\link{proteoMDS}} for peptide data
#'
#'@rdname proteoMDS
#'
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/pepseqmod_min.R
#'@example inst/extdata/examples/normPep_min.R
#'@example inst/extdata/examples/normPrn_min.R
#'@example inst/extdata/examples/imputeNA_examples.R
#' 
#'@import purrr
#'@export
pepMDS <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\MDS\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_pepid()
  
  quietly_log <- purrr::quietly(proteoMDS)(id = !!id, ...)
  purrr::walk(quietly_log[-1], write, 
              file.path(dat_dir, "Peptide\\MDS\\log\\pepMDS_log.csv"), append = TRUE)
}


#'MDS plots
#'
#'\code{prnMDS} is a wrapper of \code{\link{proteoMDS}} for protein data
#'
#'@rdname proteoMDS
#'@seealso \code{\link{load_expts}} for PSM, peptide and protein data
#'  preparation, \code{\link{pepImp}} for NA value imputation and
#'  \code{\link{pepSig}} for linear modelings.
#' @examples
#' # ===================================
#' # MDS
#' # ===================================
#' scale_log2r <- TRUE
#' 
#' # peptide
#' pepMDS(
#'   scale_log2r = TRUE,
#'   col_select = Select, 
#'   filter_by_npsm = exprs(pep_n_psm >= 10),
#'   filename = "pepMDS_filtered.png",
#' )
#'
#' # protein
#' prnMDS(
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#'   filter_by_npep = exprs(prot_n_pep >= 5),
#'   filename = "prnMDS_filtered.png",
#' )
#'
#' # custom palette
#' prnMDS(
#'   scale_log2r = TRUE,
#'   col_shape = Shape,
#'   color_brewer = Set1,
#'   show_ids = TRUE,
#'   filename = "my_palette.png",
#' )
#'
#' \dontrun{
#' prnMDS(
#'   col_color = "column_key_not_existed",
#'   col_shape = "another_missing_column_key"
#' )
#' }
#'
#'@import purrr
#'@export
prnMDS <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\MDS\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_normPSM_protid()
  
  quietly_log <- purrr::quietly(proteoMDS)(id = !!id, ...)
  purrr::walk(quietly_log[-1], write, 
              file.path(dat_dir, "Protein\\MDS\\log\\prnMDS_log.csv"), append = TRUE)
}


#'PCA plots
#'
#'\code{pepPCA} is a wrapper of \code{\link{proteoPCA}} for peptide data.
#'
#'@rdname proteoPCA
#'
#'@seealso \code{\link{load_expts}} for PSM, peptide and protein data
#'  preparation, \code{\link{pepImp}} for NA value imputation and
#'  \code{\link{pepSig}} for linear modelings.
#'  
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/pepseqmod_min.R
#'@example inst/extdata/examples/normPep_min.R
#'@example inst/extdata/examples/normPrn_min.R
#'@example inst/extdata/examples/imputeNA_examples.R
#'
#'@import purrr
#'@export
pepPCA <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\PCA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_pepid()
  
  quietly_log <- purrr::quietly(proteoPCA)(id = !!id, ...)
  purrr::walk(quietly_log[-1], write, 
              file.path(dat_dir, "Peptide\\PCA\\log\\pepPCA_log.csv"), append = TRUE)
}


#'PCA plots
#'
#'\code{prnPCA} is a wrapper of \code{\link{proteoPCA}} for protein data
#'
#'@rdname proteoPCA
#'
#' @examples
#' # ===================================
#' # PCA
#' # ===================================
#' scale_log2r <- TRUE
#' 
#' # peptide
#' pepPCA(
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#'   filter_by_npsm = exprs(pep_n_psm >= 10),
#'   filename = "pepPCA_filtered.png",
#' )
#' 
#' # protein
#' prnPCA(
#'   scale_log2r = TRUE,
#'   col_color = Color,
#'   col_shape = Shape,
#'   show_ids = TRUE,
#'   filter_by_npep = exprs(prot_n_pep >= 5),
#'   filename = "prnPCA_filtered.png",
#' )
#'
#' # by features
#' prnPCA(
#'   type = feats,
#'   scale_log2r = TRUE,
#'   filename = "prnPCA_by_feats.png",
#' )
#'
#' \dontrun{
#' prnPCA(
#'   col_color = "column_key_not_existed",
#'   col_shape = "another_missing_column_key"
#' )
#' }
#'
#'@import purrr
#'@export
prnPCA <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\PCA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_protid()
  
  quietly_log <- purrr::quietly(proteoPCA)(id = !!id, ...)
  purrr::walk(quietly_log[-1], write, 
              file.path(dat_dir, "Protein\\PCA\\log\\prnPCA_log.csv"), append = TRUE)
}


#'Distance plots
#'
#'\code{pepEucDist} is a wrapper of \code{\link{proteoEucDist}} for peptide data
#'
#'@rdname proteoEucDist
#'
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/pepseqmod_min.R
#'@example inst/extdata/examples/normPep_min.R
#'@example inst/extdata/examples/normPrn_min.R
#'@example inst/extdata/examples/imputeNA_examples.R
#'
#'@import purrr
#'@export
pepEucDist <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\EucDist\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_pepid()
  
  quietly_log <- purrr::quietly(proteoEucDist)(id = !!id, ...)
  purrr::walk(quietly_log[-1], write, 
              file.path(dat_dir, "Peptide\\EucDist\\log\\pepEucDist_log.csv"), append = TRUE)
}


#'Distance plots
#'
#'\code{prnEucDist} is a wrapper of \code{\link{proteoEucDist}} for protein data
#'
#'@rdname proteoEucDist
#'@seealso \code{\link{load_expts}} for PSM, peptide and protein data
#'  preparation, \code{\link{pepImp}} for NA value imputation and
#'  \code{\link{pepSig}} for linear modelings.
#' @examples
#' # ===================================
#' # Euclidean distance
#' # ===================================
#' scale_log2r <- TRUE
#' 
#' # peptide
#' pepEucDist(
#'   annot_cols = c("Peptide_Yield", "Group"),
#'   width = 16,
#'   height = 12,
#' )
#' 
#' # protein
#' prnEucDist(
#'   scale_log2r = TRUE,
#'   annot_cols = c("Peptide_Yield", "Group"),
#'   annot_colnames = c("New_Yield", "New_Group"),
#'
#'   # parameters for `pheatmap`
#'   display_numbers = TRUE,
#'   number_color = "grey30",
#'   number_format = "%.1f",
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
#'   cellheight = 24,
#'   filter_by_npep = exprs(prot_n_pep >= 5),
#'   filename = "prnEucDist_filtered.png",
#' )
#'
#'
#' \dontrun{
#' prnEucDist(
#'   col_color = "column_key_not_existed",
#'   col_shape = "another_missing_column_key",
#'   annot_cols = c("bad_column_key", "yet_another_bad_column_key")
#' )
#' }
#'
#'@import purrr
#'@export
prnEucDist <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\EucDist\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_normPSM_protid()
  
  quietly_log <- purrr::quietly(proteoEucDist)(id = !!id, ...)
  purrr::walk(quietly_log[-1], write, 
              file.path(dat_dir, "Protein\\EucDist\\log\\prnEucDist_log.csv"), append = TRUE)
}
