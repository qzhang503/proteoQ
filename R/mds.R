#' Plots MDS
#' 
#' @inheritParams prnMDS
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotMDS <- function (df = NULL, id = NULL, label_scheme_sub = NULL, 
                     adjEucDist = FALSE, classical = TRUE, method = "euclidean", p = 2, 
                     k = 3, show_ids = FALSE, 
                     col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, 
                     col_alpha = NULL, 
                     color_brewer = NULL, fill_brewer = NULL, 
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                     scale_log2r = TRUE, complete_cases = FALSE,
                     filepath = NULL, filename = NULL, theme = NULL, anal_type = "MDS", 
                     ...) {
  
  stopifnot(vapply(c(scale_log2r, complete_cases, adjEucDist, classical, show_ids), 
                   rlang::is_logical, logical(1)))
  stopifnot(vapply(c(p, k), is.numeric, logical(1)))
  stopifnot(nrow(label_scheme_sub) > 0)
  
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
  
  anal_dots <- dots %>% .[names(.) %in% c("d", "k", "eig", "add", "x.ret", "list.", # cmdscale
                                          "y", "maxit", "trace", "tol", # isoMDS
                                          "x", "method", "diag", "upper", "p", "m")] # dist
  dots <- dots %>% .[! . %in% anal_dots]

  if (!is.null(anal_dots$x)) warning("Argument `x` in `dist()` automated.", call. = FALSE)
  if (!is.null(anal_dots$diag)) warning("Argument `diag` in `dist()` automated.", call. = FALSE)
  if (!is.null(anal_dots$upper)) warning("Argument `upper` in `dist()` automated.", call. = FALSE)
  if (!is.null(anal_dots$m)) warning("Argument `m` in `dist()` not used.", call. = FALSE)
  
  if (!is.null(anal_dots$d)) warning("Distance object `d` automated.", call. = FALSE)
  if (!is.null(anal_dots$eig)) warning("Argument `eig` in `cmdscale()` automated.", call. = FALSE)
  if (!is.null(anal_dots$x.ret)) warning("Argument `x.ret` in `cmdscale()` automated.", call. = FALSE)
  if (!is.null(anal_dots$list.)) warning("Argument `list.` in `cmdscale()` automated.", call. = FALSE)
  
  if (!is.null(anal_dots$y)) warning("Argument `y` in `isoMDS()` automated.", call. = FALSE)

  # note that `method`, `k`, `p` already in main arguments
  anal_dots <- anal_dots %>% .[! names(.) %in% c("x", "method", "diag", "upper", "p", "m")]
  anal_dots <- anal_dots %>% .[! names(.) %in% c("d", "k", "eig", "x.ret", "list.")] # "add", 
  anal_dots <- anal_dots %>% .[! names(.) %in% c("y")] # "maxit", "trace", "tol"

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    scoreMDS(
      id = !!id, 
      label_scheme_sub = label_scheme_sub, 
      anal_type = anal_type, 
      scale_log2r = scale_log2r, 
      adjEucDist = adjEucDist, 
      classical = classical, 
      method = method, 
      p = p, 
      k = k, 
      !!!anal_dots) %>% 
    cmbn_meta(label_scheme_sub) %T>% 
    readr::write_tsv(file.path(filepath, paste0(fn_prefix, "_res.txt")))

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
	
	proteoq_mds_theme <- theme_bw() + theme(
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
	
	if (is.null(theme)) theme <- proteoq_mds_theme

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
	  theme

	if (show_ids) {
	  p <- p +
	    geom_text(data = df, 
	              mapping = aes(x = Coordinate.1, y = Coordinate.2, 
	                            label = df$Sample_ID), color = "gray", size = 3)	  
	}

	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), dots)
	do.call(ggsave, gg_args)
	
	invisible(df)
}


#' Plots EucDist
#'
#' @inheritParams prnEucDist
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 rlang pheatmap
#' @importFrom magrittr %>%
plotEucDist <- function (df = NULL, id = NULL, label_scheme_sub = NULL, adjEucDist = FALSE, 
                         scale_log2r = TRUE, complete_cases = FALSE, 
                         annot_cols, annot_colnames, 
                         filepath = NULL, filename = NULL, anal_type = "EucDist", 
                         ...) {
  stopifnot(vapply(c(scale_log2r, complete_cases, adjEucDist), 
                   rlang::is_logical, logical(1)))
  stopifnot(nrow(label_scheme_sub) > 0)
  
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]

  D <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    scoreEucDist(
      id = !!id, 
      label_scheme_sub = label_scheme_sub, 
      anal_type = anal_type, 
      scale_log2r = scale_log2r, 
      adjEucDist = adjEucDist, ) 

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

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
	n_TMT_sets <- n_TMT_sets(label_scheme_sub)
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
#' @inheritParams prnPCA
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotPCA <- function (df = NULL, id = NULL, label_scheme_sub = NULL, type = "obs", show_ids = TRUE, 
                     col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                     color_brewer = NULL, fill_brewer = NULL, 
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                     # prop_var, 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     filepath = NULL, filename = NULL, theme = NULL, 
                     anal_type = "PCA", ...) {

  stopifnot(vapply(c(scale_log2r, complete_cases, impute_na, show_ids), 
                   rlang::is_logical, logical(1)))
  stopifnot(nrow(label_scheme_sub) > 0)
  
  complete_cases <- to_complete_cases(complete_cases = complete_cases, impute_na = impute_na)
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
  
  # no need to handle "scale": prnPCA(scale = ...) parsed as prnPCA(scale_log2r = ...)
  anal_dots <- dots %>% .[names(.) %in% c("x", "retx", "center", "scale.", "tol", "rank.")]
  fml_dots <- dots[purrr::map_lgl(dots, is_formula)]
  dots <- dots %>% .[! . %in% c(anal_dots, fml_dots)]
  
  if (!is.null(anal_dots$scale.)) {
    anal_dots$scale. <- NULL
    warning("Argument `scale.` disabled; instead `scale_log2r` will be used for `scale.`.", 
            call. = FALSE)
  }
  
  if (!is.null(anal_dots$x)) {
    anal_dots$x <- NULL
    warning("Argument `x` in `prcomp()` automated.", call. = FALSE)
  }
  
  if (!purrr::is_empty(fml_dots)) {
    fml_dots <- NULL
    warning("The method for class 'formula' is not yet available in proteoQ.", call. = FALSE)
  }

  res <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    scorePCA(
      id = !!id, 
      label_scheme_sub = label_scheme_sub, 
      anal_type = anal_type, 
      scale_log2r = scale_log2r, 
      type = type, 
      !!!anal_dots)
      
  df <- res$PCA
  df <- df %>% cmbn_meta(label_scheme_sub)
  prop_var <- res$prop_var %>% 
    gsub("%", "", .) %>% 
    as.numeric() %>% 
    paste0("%")
  rm(res)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

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
	
	proteoq_pca_theme <- theme_bw() + theme(
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
	if (is.null(theme)) theme <- proteoq_pca_theme

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
	
	if (is.null(anal_dots$center)) {
	  x_lable = paste0("PC1 (", prop_var[1], ")")
	  y_lable = paste0("PC2 (", prop_var[2], ")")
	} else {
	  if (anal_dots$center) {
	    x_lable = paste0("PC1 (", prop_var[1], ")")
	    y_lable = paste0("PC2 (", prop_var[2], ")")
	  } else {
	    x_lable = "PC1"
	    y_lable = "PC2"
	  }
	}

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
		labs(title = "", x = x_lable, y = y_lable) +
		coord_fixed() +
	theme

	if (show_ids) {
	  p <- p +
	    geom_text(data = df,
	              mapping = aes(x = Coordinate.1, y = Coordinate.2, label = df$Sample_ID),
	              color = "gray", size = 3)
	}

	filename <- gg_imgname(filename)
	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), dots)
	do.call(ggsave, gg_args)
	
	invisible(list(pca = df, var = prop_var))
}


#' Scores MDS
#'
#' @inheritParams prnMDS
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scoreMDS <- function (df, id, label_scheme_sub, anal_type, scale_log2r, 
                      adjEucDist = FALSE, classical, method = "euclidean", p = 2, k = 3, ...) {

	dots <- rlang::enexprs(...)
	id <- rlang::as_string(rlang::enexpr(id))
	
	stopifnot(nrow(df) > 50)

	df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
	             sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
	  .$log2R
	
	center = FALSE # not matter to `dist` calculation
	if (center) df <- sweep(df, 1, rowMeans(df, na.rm = TRUE), "-")

	D <- dist(x = t(df), method = method, diag = TRUE, upper = TRUE, p = p)
	if (anyNA(D)) stop("Distance cannot be calculated between one more sample pairs.\n", 
	                   "Check entries under the column corresponding to `col_select` in metadata.", 
	                   call. = FALSE)
	
	D <- as.matrix(D)

	if (adjEucDist && method == "euclidean") {
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
		isomds_dots <- dots %>% .[names(.) %in% c("y", "maxit", "trace", "tol")]
		isomds_dots$y <- NULL
		
		df_mds <- rlang::expr(MASS::isoMDS(d = !!D, k = !!k, p = !!p, !!!isomds_dots)) %>% 
		  rlang::eval_bare(env = caller_env()) %>% 
		  .$points %>% 
		  data.frame()
	} else {
		cmdscale_dots <- dots %>% .[names(.) %in% c("eig", "add", "x.ret", "list.")]
		cmdscale_dots$list. <- TRUE
		
		df_mds <- rlang::expr(stats::cmdscale(d = !!D, k = !!k, !!!cmdscale_dots)) %>% 
		  rlang::eval_bare(env = caller_env()) %>% 
		  .$points %>% 
		  data.frame()
	}
	
	df_mds <- df_mds %>%
		`colnames<-`(paste("Coordinate", 1:k, sep = ".")) %>%
		tibble::rownames_to_column("Sample_ID") %>%
		dplyr::select(which(not_all_zero(.))) %>%
		tibble::column_to_rownames(var = "Sample_ID")
}


#' Scores PCA
#'
#' @inheritParams prnPCA
#' @inheritParams info_anal
#' @inheritParams gspaTest
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

  if (type == "obs") {
    df_t <- df %>% 
      t() %>%
      data.frame(check.names = FALSE) %>%
      bind_cols(label_scheme_sub)
    
    pr_out <- local({
      tempdata <- df_t %>% dplyr::select(which(colnames(.) %in% rownames(df)))
      pca_call <- rlang::expr(stats::prcomp(x = !!tempdata, scale. = !!scale_log2r, !!!dots))
      rlang::eval_bare(pca_call, env = caller_env())
    })
    
    rownames(pr_out$x) <- label_scheme_sub$Sample_ID
    
    prop_var <- summary(pr_out)$importance[2, ] %>% round(., digits = 3) %>% scales::percent()
    
    df_pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      `names<-`(gsub("^PC", "Coordinate\\.", names(.)))    
  } else if (type == "feats") {
    pr_out <- local({
      pca_call <- rlang::expr(stats::prcomp(x = !!df, scale. = !!scale_log2r, !!!dots))
      rlang::eval_bare(pca_call, env = caller_env())
    })
    
    df_pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      `names<-`(gsub("^PC", "Coordinate\\.", names(.)))
    
    prop_var <- summary(pr_out)$importance[2, ] %>% round(., digits = 3) %>% scales::percent()    
  } else {
    stop("Unkown `type` for PCA.", call. = FALSE)
  }

  return(list("PCA" = df_pca, "prop_var" = prop_var))
}


#' Scores Euclidean distance
#'
#' @inheritParams prnEucDist
#' @inheritParams info_anal
#' @inheritParams gspaTest
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


#'MDS plots
#'
#'\code{pepMDS} visualizes the multidimensional scaling (MDS) of peptide \code{log2FC}.
#'
#'@rdname prnMDS
#'
#'@import purrr
#'@export
pepMDS <- function (col_select = NULL, col_color = NULL, col_fill = NULL,
                    col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjEucDist = FALSE, classical = TRUE, 
                    method = "euclidean", p = 2, k = 3, 
                    show_ids = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
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
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, method = method, 
                               p = p, k = k, show_ids = show_ids,
                               theme = theme, ...)
  
}


#'MDS plots
#'
#'\code{prnMDS} visualizes the multidimensional scaling (MDS) of protein
#'\code{log2FC}.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}}, followed by a metric
#'(\code{\link[stats]{cmdscale}}) or non-metric (\code{\link[MASS]{isoMDS}})
#'MDS. The default is metric MDS with the input dissimilarities being euclidean
#'distances.
#'
#'@inheritParams  prnHist
#'@inheritParams prnHM
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
#'@param adjEucDist Logical; if TRUE, adjusts the inter-plex \code{Euclidean}
#'  distance by \eqn{1/sqrt(2)} at \code{method = "euclidean"}. The option
#'  \code{adjEucDist = TRUE} may be suitable when \code{reference samples} from
#'  each TMT plex undergo approximately the same sample handling process as the
#'  samples of interest. For instance, \code{reference samples} were split at
#'  the levels of protein lysates. Typically, \code{adjEucDist = FALSE} if
#'  \code{reference samples} were split near the end of a sample handling
#'  process, for instance, at the stages immediately before or after TMT
#'  labeling. Also see online
#'  \href{https://github.com/qzhang503/proteoQ}{README, section MDS} for a brief
#'  reasoning.
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
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional parameters for
#'  \code{ggsave}: \cr \code{width}, the width of plot; \cr \code{height}, the
#'  height of plot \cr \code{...}
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
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
#'@example inst/extdata/examples/prnMDS_.R
#'
#'@return MDS plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnMDS <- function (col_select = NULL, col_color = NULL, col_fill = NULL,
                    col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjEucDist = FALSE, classical = TRUE, 
                    method = "euclidean", p = 2, k = 3, 
                    show_ids = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
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
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, method = method, 
                               p = p, k = k, show_ids = show_ids,
                               theme = theme, ...)
}


#'PCA plots
#'
#'\code{prnPCA} visualizes the principal component analysis (PCA) for peptide
#'data.
#'
#'@rdname prnPCA
#'
#'@import purrr
#'@export
pepPCA <- function (col_select = NULL, col_color = NULL, 
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    show_ids = TRUE, type = c("obs", "feats"), 
                    df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  
  message("=== The feature of 'formula prcomp` is not yet available in proteoQ. ===\n", 
          "=== ONLY arguments with `Default prcomp` will be used. ===\n")
  
  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  type <- rlang::enexpr(type)
  if (type == rlang::expr(c("obs", "feats"))) {
    type <- "obs"
  } else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("obs", "feats"))
  }
  
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
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "PCA")(type = type, show_ids = show_ids, theme = theme, ...)
}


#'PCA plots
#'
#'\code{prnPCA} visualizes the principal component analysis (PCA) for protein
#'data.
#'
#'\code{log2FC} are used in PCA (\code{\link[stats]{prcomp}}).
#'
#'@inheritParams prnHist
#'@inheritParams prnHM
#'@inheritParams prnMDS
#'@param complete_cases Logical; always TRUE for PCA.
#'@param type Character string indicating the type of PCA. At the \code{type =
#'  obs} default, the components are by observations; at \code{type = feats},
#'  the components are by features.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional parameters for
#'  \code{prcomp}: \cr \code{center}, data centering; \cr \code{tol}, the
#'  tolerance \cr \code{...} \cr \cr Additional parameters for
#'  \code{ggsave}: \cr \code{width}, the width of plot; \cr \code{height}, the
#'  height of plot \cr \code{...}
#'
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
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
#'@example inst/extdata/examples/prnPCA_.R
#'
#'@return PCA plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnPCA <- function (col_select = NULL, col_color = NULL, 
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    show_ids = TRUE, type = c("obs", "feats"), 
                    df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  
  message("=== The feature of 'formula prcomp` is not yet available in proteoQ. ===\n", 
          "=== ONLY arguments with `Default prcomp` will be used. ===\n")

  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  type <- rlang::enexpr(type)
  if (type == rlang::expr(c("obs", "feats"))) {
    type <- "obs"
  } else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("obs", "feats"))
  }

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
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "PCA")(type = type, show_ids = show_ids, theme = theme, ...)
}


#'Distance plots
#'
#'\code{pepEucDist} visualizes the heat map of Euclidean distances for peptide
#'data.
#'
#'@rdname prnEucDist
#'
#'@import purrr
#'@export
pepEucDist <- function (col_select = NULL, 
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                        adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, 
                        df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "col_group", "col_color", "col_fill", 
               "col_shape", "col_size", "col_alpha", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = NULL, col_fill = NULL, 
            col_shape = NULL, col_size = NULL, col_alpha = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "EucDist")(adjEucDist = adjEucDist, 
                                   annot_cols = annot_cols, annot_colnames = annot_colnames, ...)
}


#'Distance plots
#'
#'\code{prnEucDist} visualizes the heat map of Euclidean distances for protein
#'data.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}} for heat map visualization. 
#'
#'@inheritParams prnHist
#'@inheritParams prnMDS
#'@inheritParams prnHM
#'@param annot_cols A character vector of column keys in \code{expt_smry.xlsx}.
#'  The values under the selected keys will be used to color-code sample IDs on
#'  the top of the indicated plot. The default is NULL without column
#'  annotation.
#'@param annot_colnames A character vector of replacement name(s) to
#'  \code{annot_cols}. The default is NULL without name replacement.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional parameters for
#'  plotting: \cr \code{width}, the width of plot \cr \code{height}, the height
#'  of plot \cr 
#'  \cr Additional arguments for \code{\link[pheatmap]{pheatmap}}: \cr 
#'  \code{cluster_rows, clustering_method, clustering_distance_rows}... \cr 
#'  \cr Notes about \code{pheatmap}:
#'  \cr \code{annotation_col} disabled; instead use keys indicated in \code{annot_cols}
#'  \cr \code{annotation_row} disabled; instead use keys indicated in \code{annot_rows}
#'
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
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
#'@example inst/extdata/examples/prnEucDist_.R
#'@return Heat map visualization of distance matrices.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnEucDist <- function (col_select = NULL, 
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                        adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, 
                        df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "col_group", "col_color", "col_fill", 
               "col_shape", "col_size", "col_alpha", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = NULL, col_fill = NULL, 
            col_shape = NULL, col_size = NULL, col_alpha = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "EucDist")(adjEucDist = adjEucDist, 
                                   annot_cols = annot_cols, annot_colnames = annot_colnames, ...)
}
