#' Geom plot for ggpairs
#'
#' @param data Input data.
#' @param mapping The mapping for ggplot2.
#' @param params Additional parameters for geom_point.
#' @inheritParams prnPCA
#' @import dplyr ggplot2 
#' @importFrom magrittr %>% %T>% %$% %<>% 
geom_lower_text <- function(data, mapping, params, show_ids, ...) 
{
  mapping_xy <- mapping %>% .[names(.) %in% c("x", "y")]

  p <- ggplot() +
    rlang::eval_tidy(rlang::quo(geom_point(data = data, mapping = mapping, !!!params)))

  if (show_ids) {
    stopifnot("Label" %in% names(data))

    p <- p + rlang::eval_tidy(rlang::quo(geom_text(data = data, mapping = mapping_xy,
                                                   label = data$Label,
                                                   color = "gray", ...)))
  }

  p
}


#' Plots MDS
#'
#' @inheritParams prnMDS
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotMDS <- function (df = NULL, id = NULL, label_scheme_sub = NULL,
                     choice = "cmdscale", 
                     dist_co = log2(1), adjEucDist = FALSE, 
                     method = "euclidean", p = 2,
                     k = 3, dimension = 2, folds = 1, show_ids = FALSE, 
                     show_ellipses = FALSE,
                     col_group = NULL,
                     col_color = NULL, col_fill = NULL, col_shape = NULL, 
                     col_size = NULL,
                     col_alpha = NULL,
                     color_brewer = NULL, fill_brewer = NULL,
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                     scale_log2r = TRUE, complete_cases = FALSE,
                     filepath = NULL, filename = NULL, 
                     center_features = TRUE, scale_features = TRUE,
                     theme = NULL, anal_type = "MDS",
                     ...) 
{
  stopifnot(vapply(c(scale_log2r, complete_cases, adjEucDist, 
                     show_ids, 
                     center_features, scale_features),
                   rlang::is_logical, logical(1L)))
  
  stopifnot(vapply(c(p, k), is.numeric, logical(1L)))
  
  if (!nrow(label_scheme_sub))
    stop("Empty metadata.")

  col_group <- rlang::enexpr(col_group)
  col_fill <- rlang::enexpr(col_fill)
  col_color <- rlang::enexpr(col_color)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)

  if (complete_cases) 
    df <- my_complete_cases(df, scale_log2r, label_scheme_sub)

  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  arrange_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter_dots, arrange_dots)]

  anal_dots <- dots %>% 
    .[names(.) %in% c("d", "k", "eig", "add", "x.ret", "list.", # cmdscale
                      "y", "maxit", "trace", "tol", # isoMDS
                      "x", "method", "diag", "upper", "p", "m")] # dist
  
  dots <- dots %>% 
    .[! . %in% anal_dots]
  
  # (1) overwriting args: to this <- from `cmdscale`
  # ... NA ...

  # (2) excluded formalArgs: 
  if (!is.null(anal_dots[["x"]])) 
    warning("Argument `x` in `dist()` automated.", call. = FALSE)
                                     
  if (!is.null(anal_dots[["diag"]])) 
    warning("Argument `diag` in `dist()` automated.", call. = FALSE)
  
  if (!is.null(anal_dots[["upper"]])) 
    warning("Argument `upper` in `dist()` automated.", call. = FALSE)
                                         
  if (!is.null(anal_dots[["m"]])) 
    warning("Argument `m` in `dist()` not used.", call. = FALSE)

  if (!is.null(anal_dots[["d"]])) 
    warning("Distance object `d` automated.", call. = FALSE)
  
  if (!is.null(anal_dots[["eig"]])) 
    warning("Argument `eig` in `cmdscale()` automated.", call. = FALSE)
  
  if (!is.null(anal_dots[["x.ret"]])) 
    warning("Argument `x.ret` in `cmdscale()` automated.", call. = FALSE)
  
  if (!is.null(anal_dots[["list."]])) 
    warning("Argument `list.` in `cmdscale()` automated.", call. = FALSE)

  if (!is.null(anal_dots[["y"]])) 
    warning("Argument `y` in `isoMDS()` automated.", call. = FALSE)

  # note that `method`, `k`, `p` already in main arguments
  anal_dots <- anal_dots %>% 
    .[! names(.) %in% c("x", "method", "diag", "upper", "p", "m")]
  anal_dots <- anal_dots %>% 
    .[! names(.) %in% c("d", "k", "eig", "x.ret", "list.")] # "add",
  anal_dots <- anal_dots %>% 
    .[! names(.) %in% c("y")] # "maxit", "trace", "tol"

  # (3) conversion: expr to character string
  # ... NA ...
  
  # (4) overwriting from this -> to `cmdscale` defaults
  # ... NA ...
  
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
      center_features = center_features,
      scale_features = scale_features,
      dist_co = dist_co,
      adjEucDist = adjEucDist,
      choice = choice, 
      method = method,
      p = p,
      k = k,
      col_group = !!col_group,
      folds = folds,
      out_file = file.path(filepath, paste0(fn_prefix, "_res.txt")),
      !!!anal_dots)

  # key `Label` used in `geom_lower_text()`
  df$Label <- if ("Sample_ID" %in% names(df)) 
    df$Sample_ID
  else if (id %in% names(df)) 
    df$df[[id]]
  else 
    df[, 1, drop = FALSE]

	map_color <- map_fill <- map_shape <- map_size <- map_alpha <- NA
	nms <- names(df)

	if (col_color != rlang::expr(Color) || !rlang::as_string(sym(col_color)) %in% nms)
	  assign(paste0("map_", tolower(rlang::as_string(col_color))), "X")
	
	if (col_fill != rlang::expr(Fill)  || !rlang::as_string(sym(col_fill)) %in% nms)
	  assign(paste0("map_", tolower(rlang::as_string(col_fill))), "X")
	
	if (col_shape != rlang::expr(Shape) || !rlang::as_string(sym(col_shape)) %in% nms)
	  assign(paste0("map_", tolower(rlang::as_string(col_shape))), "X")
	
	if (col_size != rlang::expr(Size) || !rlang::as_string(sym(col_size)) %in% nms)
	  assign(paste0("map_", tolower(rlang::as_string(col_size))), "X")
	
	if (col_alpha != rlang::expr(Alpha) || !rlang::as_string(sym(col_alpha)) %in% nms)
	  assign(paste0("map_", tolower(rlang::as_string(col_alpha))), "X")

	if (!is.na(map_color)) col_color <- NULL
	if (!is.na(map_fill)) col_fill <- NULL
	if (!is.na(map_shape)) col_shape <- NULL
	if (!is.na(map_size)) col_size <- NULL
	if (!is.na(map_alpha)) col_alpha <- NULL

	rm(list = c("map_color", "map_fill", "map_shape", "map_size", "map_alpha", "nms"))
	suppressWarnings(rm(list = c("map_.")))

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
		 plot.title = element_text(face="bold", colour="black", 
		                           size=20, hjust=0.5, vjust=0.5),

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

	# --- check dimension ---
	if (dimension < 2L) {
	  warning("The `dimension` increased from ", dimension, 
	          " to a minimum of 2.")
	  dimension <- 2L
	}

	ranges <- seq_len(dimension)
	cols <- names(df) %>% .[. %in% paste0("Coordinate.", ranges)]

	max_dim <- names(df) %>% .[grepl("^Coordinate\\.[0-9]+", .)] %>% length()
	
	if (dimension > max_dim) {
	  warning("The `dimension` decreased from ", dimension, 
	          " to a maximum of ", max_dim, ".")
	  dimension <- max_dim
	}
	
	rm(list = c("max_dim"))

	# --- set up aes ---
	if ((!is.null(col_color)) && rlang::as_string(col_color) == ".") 
	  col_color <- NULL
	if ((!is.null(col_fill)) && rlang::as_string(col_fill) == ".") 
	  col_fill <- NULL
	if ((!is.null(col_shape)) && rlang::as_string(col_shape) == ".") 
	  col_shape <- NULL
	if ((!is.null(col_size)) && rlang::as_string(col_size) == ".") 
	  col_size <- NULL
	if ((!is.null(col_alpha)) && rlang::as_string(col_alpha) == ".") 
	  col_alpha <- NULL
	
	if (dimension > 2L) {
	  mapping <- ggplot2::aes(colour = !!col_color, 
	                          fill = !!col_fill, 
	                          shape = !!col_shape,
	                          size = !!col_size, 
	                          alpha = !!col_alpha)
	} 
	else {
	  mapping <- ggplot2::aes(x = Coordinate.1, 
	                          y = Coordinate.2,
	                          colour = !!col_color, 
	                          fill = !!col_fill, 
	                          shape = !!col_shape,
	                          size = !!col_size, 
	                          alpha = !!col_alpha)
	}

	idxes <- purrr::map(mapping, `[[`, 1) %>% purrr::map_lgl(is.null)

	mapping_var <- mapping[!idxes]
	mapping_fix <- mapping[idxes]
	
	local({
	  nms <- names(mapping_var)
	  not_xy <- which(!nms %in% c("x", "y"))
	  
	  vars <- mapping_var[not_xy]
	  
	  if (length(vars)) {
	    for (var in vars) {
	      col <- quo_name(var)
	      
	      if (anyNA(df[[col]])) {
	        warning("NA/incomplete aesthetics in column `", col, "`.\n")
	      }
	    }
	  }
	})

	fix_args <- list(colour = "darkgray", 
	                 fill = NA, 
	                 shape = 21, 
	                 size = 4, 
	                 alpha = 0.9) %>%
	  .[names(.) %in% names(mapping_fix)] %>%
	  .[!is.na(.)]
	
	fix_args$stroke <- 0.02

	# --- set up axis labels ---
	col_labs <- cols %>% gsub("\\.", " ", .)

	# --- plots ---
	if (dimension > 2L) {
	  p <- GGally::ggpairs(df,
	                       axisLabels = "internal",
	                       columns = cols,
	                       mapping = mapping_var,
	                       columnLabels = col_labs,
	                       labeller = label_wrap_gen(10),
	                       title = "",
	                       lower = list(continuous = wrap(geom_lower_text,
	                                                      params = fix_args,
	                                                      show_ids = show_ids,
	                                                      size = 3)),
	                       upper = "blank")

	  p <- p + theme

	  if (!is.null(fill_brewer)) {
	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_fill_brewer(palette = fill_brewer)
	      }
	    }
	  }

	  if (!is.null(color_brewer)) {
	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_color_brewer(palette = color_brewer)
	      }
	    }
	  }

	  if ((!is.null(col_size)) && (!is.null(size_manual))) {
	    check_aes_length(label_scheme_sub, col_size, "size_manual", size_manual)

	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_size_manual(values = size_manual)
	      }
	    }
	  }

	  if ((!is.null(col_shape)) && (!is.null(shape_manual))) {
	    check_aes_length(label_scheme_sub, col_shape, "shape_manual", shape_manual)

	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_shape_manual(values = shape_manual)
	      }
	    }
	  }

	  if ((!is.null(col_alpha)) && (!is.null(alpha_manual))) {
	    check_aes_length(label_scheme_sub, col_alpha, "alpha_manual", alpha_manual)

	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_shape_manual(values = alpha_manual)
	      }
	    }
	  }

	  if (show_ellipses) {
	    if (anyNA(label_scheme_sub[[col_group]])) {
	      warning("(Partial) NA aesthetics under column `", col_group, 
	              "` in expt_smry.xlsx",
	              call. = FALSE)
	    }

	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + stat_ellipse(
	          data = df,
	          aes(x = !!rlang::sym(paste("Coordinate", y, sep = ".")),
	              y = !!rlang::sym(paste("Coordinate", x, sep = ".")),
	              fill = !!rlang::sym(col_group)),
	          geom = "polygon",
	          alpha = .4,
	          show.legend = FALSE,
	        )
	      }
	    }
	  }

	}
	else {
	  p <- ggplot() +
	    rlang::eval_tidy(rlang::quo(geom_point(data = df, 
	                                           mapping = mapping_var, 
	                                           !!!fix_args))) +
	    coord_fixed()

	  check_ggplot_aes(p)
	  
	  if (show_ellipses) {
	    if (anyNA(label_scheme_sub[[col_group]])) {
	      warning("(Partial) NA aesthetics under column `", col_group, 
	              "` in expt_smry.xlsx")
	    }

	    p <- p + stat_ellipse(
	      data = df,
	      aes(x = Coordinate.1, y = Coordinate.2, fill = !!rlang::sym(col_group)),
	      geom = "polygon",
	      alpha = .4,
	      show.legend = FALSE,
	    )
	  }

	  if (show_ids) {
	    p <- p +
	      geom_text(data = df,
	                mapping = aes(x = Coordinate.1, y = Coordinate.2, label = Sample_ID),
	                color = "gray", size = 3)
	  }

	  p <- p +
	    labs(title = "", x = col_labs[1], y = col_labs[2]) + theme

  	if (!is.null(fill_brewer)) p <- p + scale_fill_brewer(palette = fill_brewer)
  	if (!is.null(color_brewer)) p <- p + scale_color_brewer(palette = color_brewer)

  	if ((!is.null(col_size)) && (!is.null(size_manual))) {
  	  check_aes_length(label_scheme_sub, col_size, "size_manual", size_manual)
  	  p <- p + scale_size_manual(values = size_manual)
  	}

  	if ((!is.null(col_shape)) && (!is.null(shape_manual))) {
  	  check_aes_length(label_scheme_sub, col_shape, "shape_manual", shape_manual)
  	  p <- p + scale_shape_manual(values = shape_manual)
  	}

  	if ((!is.null(col_alpha)) && (!is.null(alpha_manual))) {
  	  check_aes_length(label_scheme_sub, col_alpha, "alpha_manual", alpha_manual)
  	  p <- p + scale_shape_manual(values = alpha_manual)
  	}
	}

	ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot"))
	rlang::eval_tidy(rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(filename)),
	                                   plot = p, 
	                                   !!!ggsave_dots)))

	invisible(df)
}


#' Plots EucDist
#'
#' @inheritParams prnEucDist
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 pheatmap
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotEucDist <- function (df = NULL, id = NULL, label_scheme_sub = NULL, adjEucDist = FALSE,
                         scale_log2r = TRUE, complete_cases = FALSE,
                         annot_cols, annot_colnames,
                         filepath = NULL, filename = NULL, anal_type = "EucDist",
                         ...) 
{
  stopifnot(vapply(c(scale_log2r, complete_cases, adjEucDist),
                   rlang::is_logical, logical(1L)))
  
  if (!nrow(label_scheme_sub))
    stop("Empty metadata.")

  if (complete_cases) 
    df <- my_complete_cases(df, scale_log2r, label_scheme_sub)

  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  arrange_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter_dots, arrange_dots)]

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

	mypalette <- if (is.null(dots$color)) {
	  colorRampPalette(c("blue", "white", "red"))(n_color)
	} 
	else {
	  eval(dots$color, envir = rlang::caller_env())
	}

	if (is.null(annot_cols)) {
	  annotation_col <- NA
	} 
	else {
	  annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = rownames(D))
	  idx <- which(annot_cols %in% colnames(annotation_col))

	  annot_cols <- annot_cols[idx]
	  annot_colnames <- annot_colnames[idx]
	}

	if (!is.null(annot_colnames) && length(annot_colnames) == length(annot_cols)) {
	  colnames(annotation_col) <- annot_colnames
	}

	annotation_colors <- if (is.null(dots$annotation_colors)) {
		setHMColor(annotation_col)
	} 
	else if (is.na(dots$annotation_colors)) {
		NA
	} 
	else {
		eval(dots$annotation_colors, envir = rlang::caller_env())
	}

	nm_idx <- names(dots) %in% c("mat", "filename", "annotation_col", "color",
	                             "annotation_colors", "breaks")
	dots[nm_idx] <- NULL

	n_TMT_sets <- n_TMT_sets(label_scheme_sub)
	max_width <- 77

	if (is.null(dots$width)) 
	  dots$width <- min(10*n_TMT_sets*1.2, max_width)

	if (dots$width >= max_width) {
		message("The width for the graphic device is", dots$width, "inches or more.")
		stop("Please consider a a smaller `cellwidth`.", call. = FALSE)
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


#' Scores MDS
#'
#' @param out_file A file path object to an output file.
#' @inheritParams prnMDS
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr 
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>% %T>% %$% %<>% 
scoreMDS <- function (df, id, label_scheme_sub, anal_type, scale_log2r,
                      center_features, scale_features,
                      dist_co = log2(1), adjEucDist = FALSE, 
                      choice = "cmdscale", 
                      method = "euclidean",
                      p = 2L, k = 3L, col_group, folds, out_file, ...) 
{
	dots <- rlang::enexprs(...)
	id <- rlang::as_string(rlang::enexpr(id))
	col_group <- rlang::enexpr(col_group)

	if (nrow(df) <= 50L)
	  stop("Need more than 50 entries for MDS.")

	df <- prepDM(df = df, id = !!id, 
	             scale_log2r = scale_log2r,
	             sub_grp = label_scheme_sub$Sample_ID, 
	             anal_type = anal_type, 
	             rm_allna = TRUE) %>%
	  .$log2R

	nms <- names(df)
	n_rows <- nrow(df)

	label_scheme_sub <- label_scheme_sub %>% 
	  dplyr::filter(Sample_ID %in% nms)

	res <- prep_folded_tdata(df, folds, label_scheme_sub, !!col_group)
	df_t <- res$df_t
	ls_sub <- res$ls_sub
	rm(list = "res")

	if (rlang::as_string(col_group) %in% names(df_t)) {
	  df_t <- df_t %>% dplyr::select(-!!rlang::sym(col_group))
	}

	stopifnot(vapply(df_t, is.numeric, logical(1L)))
	
	# `center_features` not affect `dist` but `scale` calculation
	if (scale_features) {
	  df_t <- df_t %>% scale(center = center_features, scale = TRUE)
	} 
	else if (center_features) {
	  message("Distance measures not affected by data centering.")
	  df_t <- df_t %>% sweep(., 2, colMeans(., na.rm = TRUE), "-")
	}
  
	if (dist_co > 0) {
	  D <- local({
	    len <- nrow(df_t)
	    D <- matrix(ncol = len, nrow = len)
	    colnames(D) <- rownames(D) <- rownames(df_t)
	    
	    for (i in seq_len(len)) {
	      for (j in 1:i) {
	        x_i <- df_t[i, ]
	        x_j <- df_t[j, ]
	        oks <- (abs(x_i - x_j) < dist_co)
	        
	        x_i[oks] <- x_j[oks] <- NA
	        x <- rbind(x_i, x_j)
	        
	        D[j, i] <- D[i, j] <- 
	          dist(x = x, method = method, diag = TRUE, upper = TRUE, p = p) %>% 
	          `[`(1)
	      }
	      
	      D[i, i] <- 0
	    }
	    
	    invisible(D)
	  })
	} 
	else if (identical(dist_co, 0)) {
	  D <- as.matrix(dist(x = df_t, method = method, diag = TRUE, upper = TRUE, p = p))
	} 
	else {
	  stop("`dist_co` cannot be negative.", call. = FALSE)
	}
	
	if (anyNA(D)) {
	  stop("Distance cannot be calculated between one or more sample pairs.\n", 
	       "Check entries under the column corresponding to `col_select` in metadata.", 
	       call. = FALSE)
	}

	if (adjEucDist && method == "euclidean") {
	  D <- local({
	    annotation_col <- colAnnot(annot_cols = c("TMT_Set"), sample_ids = nms) %>%
	      .[rep(seq_len(nrow(.)), folds), , drop = FALSE] %>%
	      `rownames<-`(paste(rep(nms, folds), rep(seq_len(folds), each = length(nms)), 
	                         sep = "."))

	    for (i in 1:ncol(D)) {
	      for (j in 1:ncol(D)) {
	        if (annotation_col$TMT_Set[i] != annotation_col$TMT_Set[j])
	          D[i, j] <- D[i, j]/sqrt(2)
	      }
	    }

	    D
	  })
	}
  
	if (choice == "cmdscale") {
	  cmdscale_dots <- dots %>% .[names(.) %in% c("eig", "add", "x.ret", "list.")]
	  cmdscale_dots$list. <- TRUE
	  
	  df_mds <- rlang::expr(stats::cmdscale(d = !!D, k = !!k, !!!cmdscale_dots)) %>%
	    rlang::eval_bare(env = caller_env()) %>%
	    .$points %>%
	    data.frame()
	} 
	else if (choice == "isoMDS") {
	  isomds_dots <- dots %>% .[names(.) %in% c("y", "maxit", "trace", "tol")]
	  isomds_dots$y <- NULL
	  
	  df_mds <- rlang::expr(MASS::isoMDS(d = !!D, k = !!k, p = !!p, !!!isomds_dots)) %>%
	    rlang::eval_bare(env = caller_env()) %>%
	    .$points %>%
	    data.frame()
	} 
	else {
	  stop("Unknown choice = ", choice, call. = FALSE)
	}
	
	df_mds %>%
	  `colnames<-`(paste("Coordinate", 1:k, sep = ".")) %>%
	  tibble::rownames_to_column("Sample_ID") %>%
	  dplyr::select(which(not_all_zero(.))) %>%
	  `rownames<-`(NULL) %>% 
	  tibble::column_to_rownames(var = "Sample_ID") %>%
	  cmbn_meta(ls_sub) %T>%
	  readr::write_tsv(out_file)
}


#' Prepares folded, transposed data
#'
#' @param df Input data frame
#' @param label_scheme_sub Metadata for data subset
#' @inheritParams prnPCA
#' @inheritParams anal_pepNMF
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
prep_folded_tdata <- function (df, folds, label_scheme_sub, col_group) 
{
  nms <- names(df)
  n_rows <- nrow(df)

  # not used
  col_group = rlang::enexpr(col_group)

  if (folds == 1) {
    nms_feat <- rownames(df)
    nms_smpl <- nms
    ls_sub <- label_scheme_sub
    df_t <- df %>% t()
  } 
  else {
    n_feats <- floor(n_rows/folds)
    nms_feat <- paste("x", seq_len(n_feats), sep = ".")
    nms_smpl <- paste(rep(nms, folds), rep(seq_len(folds), each = length(nms)), sep = ".")
    
    ls_sub <- label_scheme_sub %>%
      .[rep(seq_len(nrow(.)), folds), , drop = FALSE] %>%
      dplyr::mutate(Sample_ID = nms_smpl)

    df_t <- purrr::map(seq_len(folds), ~ {
      df[sample(n_rows, n_feats, replace = FALSE), ] %>%
        `rownames<-`(nms_feat) %>%
        t()
    }) %>% do.call(rbind, .)
  }

  # need col_group column for LDA
  df_t <- df_t %>%
    data.frame(check.names = FALSE) %>%

    dplyr::mutate(Sample_ID = rep(nms, folds)) %>%
    dplyr::left_join(label_scheme_sub %>% 
                       dplyr::select(Sample_ID, !!rlang::sym(col_group)), 
                     by = "Sample_ID") %>%
    dplyr::select(-Sample_ID) %>%

    `rownames<-`(nms_smpl)

  invisible(list(df_t = df_t, ls_sub = ls_sub))
}


#' Scores Euclidean distance
#'
#' @inheritParams prnEucDist
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr 
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>% %T>% %$% %<>% 
scoreEucDist <- function (df, id, label_scheme_sub, anal_type, 
                          scale_log2r, adjEucDist = FALSE, ...) 
{
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))

  if (nrow(df) <= 50L)
    stop("Need more than 50 entries for distance calculations.")
  
  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r,
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type, 
               rm_allna = TRUE) %>%
    .$log2R

  D <- dist(t(df), method = "euclidean", diag = TRUE, upper = TRUE)
  
  if (anyNA(D)) 
    stop("Distance cannot be calculated for one more sample pairs.")

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

      D
    })
  }

  invisible(D)
}


#'MDS plots
#'
#'\code{pepMDS} visualizes the multidimensional scaling (MDS) of peptide \code{log2FC}.
#'
#'@rdname prnMDS
#'
#'@import purrr
#'@export
pepMDS <- function (col_select = NULL, col_group = NULL, col_color = NULL, 
                    col_fill = NULL,
                    col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL,
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                    choice = c("cmdscale", "isoMDS"), 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                    dist_co = log2(1), adjEucDist = FALSE, 
                    method = "euclidean", p = 2, k = 3, dimension = 2, folds = 1,
                    center_features = TRUE, scale_features = TRUE,
                    show_ids = TRUE, show_ellipses = FALSE,
                    df = NULL, filepath = NULL, filename = NULL,
                    theme = NULL, ...) 
{
  old_opts <- options()
  options(warn = 1, warnPartialMatchArgs = TRUE)
  on.exit(options(old_opts), add = TRUE)
  
  check_dots(c("id", "df2", "anal_type"), ...)
  check_depreciated_args(list(c("classical", "choice")), ...)
  
  purrr::walk2(formals(pepMDS)[["choice"]] %>% eval(), 
               list(c("d", "k"), c("d", "k", "p")), 
               ~ check_formalArgs(pepMDS, !!.x, .y))
  
  check_formalArgs(pepMDS, dist, c("method", "p"))
  
  choice <- rlang::enexpr(choice)
  choice <- if (length(choice) > 1L) "cmdscale" else rlang::as_string(choice)

  id <- match_call_arg(normPSM, group_psm_by)
  
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

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

  method <- rlang::enexpr(method)
  method <- if (length(method) > 1L) "euclidean" else rlang::as_string(method)

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()

  info_anal(id = !!id,
            col_select = !!col_select, 
            col_group = !!col_group,
            col_color = !!col_color, 
            col_fill = !!col_fill,
            col_shape = !!col_shape, 
            col_size = !!col_size, 
            col_alpha = !!col_alpha,
            color_brewer = !!color_brewer, 
            fill_brewer = !!fill_brewer,
            size_manual = !!size_manual, 
            shape_manual = !!shape_manual, 
            alpha_manual = !!alpha_manual,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "MDS")(choice = choice, 
                               dist_co = dist_co, 
                               adjEucDist = adjEucDist, 
                               method = method,
                               p = p, 
                               k = k, 
                               dimension = dimension, 
                               folds = folds, 
                               show_ids = show_ids,
                               show_ellipses = show_ellipses, 
                               center_features = center_features,
                               scale_features = scale_features,
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
#'distances. Note that the \code{center_features} alone will not affect the
#'results of \code{\link[stats]{dist}}; it together with \code{scale_features}
#'will be passed to \code{\link[base]{scale}}.
#'
#'@inheritParams  prnHist
#'@inheritParams prnHM
#'@inheritParams anal_prnNMF
#'@param col_color Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{color} aesthetics in plots. At
#'  the NULL default, the column key \code{Color} will be used. If NA, bypasses
#'  the aesthetics (a means to bypass the look-up of column \code{Color} and
#'  handle duplication in aesthetics).
#'@param  col_fill Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{fill} aesthetics in plots. At
#'  the NULL default, the column key \code{Fill} will be used. If NA, bypasses
#'  the aesthetics (a means to bypass the look-up of column \code{Fill} and
#'  handle duplication in aesthetics).
#'@param col_shape Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{shape} aesthetics in plots. At
#'  the NULL default, the column key \code{Shape} will be used. If NA, bypasses
#'  the aesthetics (a means to bypass the look-up of column \code{Shape} and
#'  handle duplication in aesthetics).
#'@param col_size Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{size} aesthetics in plots. At
#'  the NULL default, the column key \code{Size} will be used. If NA, bypasses
#'  the aesthetics (a means to bypass the look-up of column \code{Size} and
#'  handle duplication in aesthetics).
#'@param col_alpha Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{alpha} (transparency)
#'  aesthetics in plots. At the NULL default, the column key \code{Alpha} will
#'  be used. If NA, bypasses the aesthetics (a means to bypass the look-up of
#'  column \code{Alpha} and handle duplication in aesthetics).
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
#'@param method Character string; the distance measure in one of c("euclidean",
#'  "maximum", "manhattan", "canberra", "binary") for \code{\link[stats]{dist}}.
#'  The default method is "euclidean".
#'@param p Numeric; The power of the Minkowski distance in
#'  \code{\link[stats]{dist}}. The default is 2.
#'@param k Numeric; The desired dimension for the solution passed to
#'  \code{\link[stats]{cmdscale}}. The default is 3.
#'@param dimension Numeric; The desired dimension for pairwise visualization.
#'  The default is 2.
#'@param dist_co Numeric; The cut-off in the absolute distance measured by
#'  \eqn{d = abs(x_i - x_j)}. Data pairs, \eqn{x_i} and \eqn{x_j}, with
#'  corresponding \eqn{d} smaller than \code{dist_co} will be excluded from
#'  distance calculations by \link[stats]{dist}. The default is no distance
#'  cut-off at \eqn{dist_co = log2(1)}.
#'@param show_ids Logical; if TRUE, shows the sample IDs in \code{MDS/PCA}
#'  plots. The default is TRUE.
#'@param show_ellipses Logical; if TRUE, shows the ellipses by sample groups
#'  according to \code{col_group}. The default is FALSE.
#'@param choice Character string; the MDS method in \code{c("cmdscale",
#'  "isoMDS")}. The default is "cmdscale".
#'@inheritParams prnPCA
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  Additional parameters for \code{ggsave}: \cr \code{width}, the width of
#'  plot; \cr \code{height}, the height of plot \cr \code{...}
#'@seealso \emph{Metadata} \cr \code{\link{load_expts}} for metadata preparation
#'  and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr \code{\link{normPSM}} for extended examples in
#'  PSM data normalization \cr \code{\link{PSM2Pep}} for extended examples in
#'  PSM to peptide summarization \cr \code{\link{mergePep}} for extended
#'  examples in peptide data merging \cr \code{\link{standPep}} for extended
#'  examples in peptide data normalization \cr \code{\link{Pep2Prn}} for
#'  extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization.
#'  \cr \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples
#'  in data purging \cr \code{\link{pepHist}} and \code{\link{prnHist}} for
#'  extended examples in histogram visualization. \cr \code{\link{extract_raws}}
#'  and \code{\link{extract_psm_raws}} for extracting MS file names \cr
#'
#'  \emph{Variable arguments of `filter_...`} \cr \code{\link{contain_str}},
#'  \code{\link{contain_chars_in}}, \code{\link{not_contain_str}},
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}},
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings
#'  \cr
#'
#'  \emph{Missing values} \cr \code{\link{pepImp}} and \code{\link{prnImp}} for
#'  missing value imputation \cr
#'
#'  \emph{Informatics} \cr \code{\link{pepSig}} and \code{\link{prnSig}} for
#'  significance tests \cr \code{\link{pepVol}} and \code{\link{prnVol}} for
#'  volcano plot visualization \cr \code{\link{prnGSPA}} for gene set enrichment
#'  analysis by protein significance pVals \cr \code{\link{gspaMap}} for mapping
#'  GSPA to volcano plot visualization \cr \code{\link{prnGSPAHM}} for heat map
#'  and network visualization of GSPA results \cr \code{\link{prnGSVA}} for gene
#'  set variance analysis \cr \code{\link{prnGSEA}} for data preparation for
#'  online GSEA. \cr \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS
#'  visualization \cr \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA
#'  visualization \cr \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA
#'  visualization \cr \code{\link{pepHM}} and \code{\link{prnHM}} for heat map
#'  visualization \cr \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'  \code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for
#'  correlation plots \cr \code{\link{anal_prnTrend}} and
#'  \code{\link{plot_prnTrend}} for trend analysis and visualization \cr
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}},
#'  \code{\link{plot_pepNMFCon}}, \code{\link{plot_prnNMFCon}},
#'  \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr
#'
#'  \emph{Custom databases} \cr \code{\link{Uni2Entrez}} for lookups between
#'  UniProt accessions and Entrez IDs \cr \code{\link{Ref2Entrez}} for lookups
#'  among RefSeq accessions, gene names and Entrez IDs \cr \code{\link{prepGO}}
#'  for
#'  \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'   ontology}} \cr \code{\link{prepMSig}} for
#'  \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular
#'   signatures} \cr \code{\link{prepString}} and \code{\link{anal_prnString}}
#'  for STRING-DB \cr
#'
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@example inst/extdata/examples/prnMDS_.R
#'
#'@return MDS plots.
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
prnMDS <- function (col_select = NULL, col_group = NULL, col_color = NULL, 
                    col_fill = NULL,
                    col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL,
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                    choice = c("cmdscale", "isoMDS"), 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                    dist_co = log2(1), adjEucDist = FALSE, 
                    method = "euclidean", p = 2, k = 3, dimension = 2, folds = 1,
                    center_features = TRUE, scale_features = TRUE,
                    show_ids = TRUE, show_ellipses = FALSE,
                    df = NULL, filepath = NULL, filename = NULL,
                    theme = NULL, ...) 
{
  old_opts <- options()
  options(warn = 1, warnPartialMatchArgs = TRUE)
  on.exit(options(old_opts), add = TRUE)
  
  check_dots(c("id", "df2", "anal_type"), ...)
  check_depreciated_args(list(c("classical", "choice")), ...)
  
  purrr::walk2(formals(prnMDS)[["choice"]] %>% eval(), 
               list(c("d", "k"), c("d", "k", "p")), 
              ~ check_formalArgs(prnMDS, !!.x, .y))
  
  check_formalArgs(prnMDS, dist, c("method", "p"))
  
  choice <- rlang::enexpr(choice)
  choice <- if (length(choice) > 1L) "cmdscale" else rlang::as_string(choice)

  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

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

  method <- rlang::enexpr(method)
  method <- if (length(method) > 1L) "euclidean" else rlang::as_string(method)

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()

  info_anal(id = !!id,
            col_select = !!col_select, 
            col_group = !!col_group,
            col_color = !!col_color, 
            col_fill = !!col_fill,
            col_shape = !!col_shape, 
            col_size = !!col_size, 
            col_alpha = !!col_alpha,
            color_brewer = !!color_brewer, 
            fill_brewer = !!fill_brewer,
            size_manual = !!size_manual, 
            shape_manual = !!shape_manual, 
            alpha_manual = !!alpha_manual,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "MDS")(choice = choice, 
                               dist_co = dist_co, 
                               adjEucDist = adjEucDist, 
                               method = method,
                               p = p, 
                               k = k, 
                               dimension = dimension, 
                               folds = folds, 
                               show_ids = show_ids,
                               show_ellipses = show_ellipses, 
                               center_features = center_features,
                               scale_features = scale_features,
                               theme = theme, ...)
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
                        df = NULL, filepath = NULL, filename = NULL, ...) 
{
  check_dots(c("id", "col_group", "col_color", "col_fill",
               "col_shape", "col_size", "col_alpha", "anal_type", "df2"), ...)
  check_formalArgs(prnEucDist, dist)

  id <- match_call_arg(normPSM, group_psm_by)
  
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()

  info_anal(id = !!id,
            col_select = !!col_select, 
            col_group = NULL, 
            col_color = NULL, 
            col_fill = NULL,
            col_shape = NULL, 
            col_size = NULL, 
            col_alpha = NULL,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "EucDist")(adjEucDist = adjEucDist,
                                   annot_cols = annot_cols, 
                                   annot_colnames = annot_colnames, ...)
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
#'  \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr
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
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@example inst/extdata/examples/prnEucDist_.R
#'@return Heat map visualization of distance matrices.
#'
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
prnEucDist <- function (col_select = NULL,
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                        adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL,
                        df = NULL, filepath = NULL, filename = NULL, ...) 
{
  check_dots(c("id", "col_group", "col_color", "col_fill",
               "col_shape", "col_size", "col_alpha", "anal_type", "df2"), ...)
  check_formalArgs(prnEucDist, dist)

  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()

  info_anal(id = !!id,
            col_select = !!col_select, 
            col_group = NULL, 
            col_color = NULL, 
            col_fill = NULL,
            col_shape = NULL, 
            col_size = NULL, 
            col_alpha = NULL,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "EucDist")(adjEucDist = adjEucDist,
                                   annot_cols = annot_cols, 
                                   annot_colnames = annot_colnames, ...)
}

