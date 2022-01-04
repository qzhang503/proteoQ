#' Plots LDA
#' 
#' @inheritParams prnLDA
#' @inheritParams prnPCA
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotLDA <- function (df = NULL, id = NULL, label_scheme_sub = NULL, 
                     choice = "lda", 
                     method = "moment",
                     type = "obs",
                     dimension = 2, folds = 1,
                     show_ids = TRUE, show_ellipses = FALSE,
                     col_group = NULL,
                     col_color = NULL, col_fill = NULL, col_shape = NULL,
                     col_size = NULL, col_alpha = NULL,
                     color_brewer = NULL, fill_brewer = NULL,
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                     filepath = NULL, filename = NULL,
                     center_features = TRUE, scale_features = TRUE,
                     theme = NULL,
                     anal_type = "PCA", 
                     
                     formula = NULL, 
                     data = NULL, 
                     x = NULL, 
                     grouping = NULL, 
                     prior = NULL, 
                     subset = NULL, 
                     CV = NULL, 
                     na.action = NULL, 
                     nu = NULL, 
                     ...) 
{
  stopifnot(vapply(c(scale_log2r, complete_cases, impute_na, show_ids,
                     center_features, scale_features),
                   rlang::is_logical, logical(1)))
  stopifnot(nrow(label_scheme_sub) > 0)

  col_group <- rlang::enexpr(col_group)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)

  complete_cases <- to_complete_cases(complete_cases = complete_cases, 
                                      impute_na = impute_na)
  if (complete_cases) 
    df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

  if (show_ellipses && type == "feats") {
    show_ellipses <- FALSE
    warning("No ellipses at `type = feats`.", 
            call. = FALSE)
  }

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

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  
  if (choice == "lda") {
    dummies <- c("formula", "data", "x", "grouping", "prior", 
                 "subset", "CV", "na.action", "nu")
    
    msgs <- c(
      "`formula` in `lda()` not used.",
      "`data` in `lda()` not used.",
      "`x` in `lda()` automated.",
      paste0("`grouping` in `lda()` disabled; \n", 
             "instead, use `col_group` to indicate the grouping variable."), 
      "`prior` in `lda()` automated.",
      "`subset` in reduced-rank `lda()` not used.",
      "`CV` in `lda()` automated.",
      "`na.action` in `lda()` automated.",
      "`nu` in `lda()` not used."
    )
    
    stopifnot(length(dummies) == length(msgs))
    
    check_formalArgs(prnLDA, lda, dummies)
    check_formalArgs(pepLDA, lda, dummies)
    
    purrr::walk2(dummies, msgs, ~ {
      if (!is.null(get(.x, envir = rlang::env_parent(), inherits = FALSE))) {
        warning(.y, call. = FALSE)
        assign(.x, NULL, envir = rlang::env_parent(), inherits = FALSE)
      } 
    })
    
    # (4) `anal_dots`
    anal_dots <- dots %>% 
      .[names(.) %in% c("formula", "data", "x", "grouping", "prior", "tol", 
                        "subset", "na.action", "method", "CV", "nu")] %>% 
      .[! names(.) %in% dummies]
    
    dots <- dots %>% 
      .[! . %in% anal_dots]
  }
  

  res <- df %>%
    filters_in_call(!!!filter_dots) %>%
    arrangers_in_call(!!!arrange_dots) %>%
    scoreLDA(
      id = !!id,
      label_scheme_sub = label_scheme_sub,
      anal_type = anal_type,
      scale_log2r = scale_log2r,
      center_features = center_features,
      scale_features = scale_features,
      choice = choice, 
      method = method,
      type = type,
      col_group = !!col_group,
      folds = folds,
      out_file = file.path(filepath, paste0(fn_prefix, "_res.txt")),
      !!!anal_dots
      )

  df <- res$x

  df$Label <- if ("Sample_ID" %in% names(df))
    df$Sample_ID
  else if (id %in% names(df))
    df[[id]]
  else 
    df$df[, 1, drop = FALSE]

	map_color <- map_fill <- map_shape <- map_size <- map_alpha <- NA

	if (col_color != rlang::expr(Color) || !rlang::as_string(sym(col_color)) %in% names(df))
	  assign(paste0("map_", tolower(rlang::as_string(col_color))), "X")
	if (col_fill != rlang::expr(Fill)  || !rlang::as_string(sym(col_fill)) %in% names(df))
	  assign(paste0("map_", tolower(rlang::as_string(col_fill))), "X")
	if (col_shape != rlang::expr(Shape) || !rlang::as_string(sym(col_shape)) %in% names(df))
	  assign(paste0("map_", tolower(rlang::as_string(col_shape))), "X")
	if (col_size != rlang::expr(Size) || !rlang::as_string(sym(col_size)) %in% names(df))
	  assign(paste0("map_", tolower(rlang::as_string(col_size))), "X")
	if (col_alpha != rlang::expr(Alpha) || !rlang::as_string(sym(col_alpha)) %in% names(df))
	  assign(paste0("map_", tolower(rlang::as_string(col_alpha))), "X")

	if (!is.na(map_color)) col_color <- NULL
	if (!is.na(map_fill)) col_fill <- NULL
	if (!is.na(map_shape)) col_shape <- NULL
	if (!is.na(map_size)) col_size <- NULL
	if (!is.na(map_alpha)) col_alpha <- NULL

	rm(list = c("map_color", "map_fill", "map_shape", "map_size", "map_alpha"))
	suppressWarnings(rm(list = c("map_.")))

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
	if (is.null(theme)) theme <- proteoq_pca_theme

	# --- check dimension ---
	if (dimension < 2) {
	  warning("The `dimension` increased from ", dimension, " to a minimum of 2.", 
	          call. = FALSE)
	  dimension <- 2
	}

	ranges <- seq_len(dimension)
	cols <- names(df) %>% .[. %in% paste0("LD", ranges)]

	max_dim <- names(df) %>% .[grepl("^LD[0-9]+", .)] %>% length()
	if (dimension > max_dim) {
	  warning("The `dimension` decreased from ", dimension, 
	          " to a maximum of ", max_dim, ".",
	          call. = FALSE)
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
	
	if (dimension > 2) {
	  mapping <- ggplot2::aes(colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
	                          size = !!col_size, alpha = !!col_alpha)
	} else {
	  mapping <- ggplot2::aes(x = LD1, y = LD2,
	                          colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
	                          size = !!col_size, alpha = !!col_alpha)
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
	        warning("NA/incomplete aesthetics in column `", col, "`.\n", 
	                call. = FALSE)
	      }
	    }
	  }
	})
	
	fix_args <- list(colour = "darkgray", fill = NA, shape = 21, size = 4, alpha = 0.9) %>%
	  .[names(.) %in% names(mapping_fix)] %>%
	  .[!is.na(.)]
	fix_args$stroke <- 0.02

	# --- set up axis labels ---
	col_labs <- cols

	# --- plots ---
	if (dimension > 2) {
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

	  if ((!is.null(col_size)) & (!is.null(size_manual))) {
	    stopifnot(length(unique(label_scheme_sub[[col_size]])) == length(size_manual))
	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_size_manual(values = size_manual)
	      }
	    }
	  }

	  if ((!is.null(col_shape)) & (!is.null(shape_manual))) {
	    stopifnot(length(unique(label_scheme_sub[[col_shape]])) == length(shape_manual))
	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_shape_manual(values = shape_manual)
	      }
	    }
	  }

	  if ((!is.null(col_alpha)) & (!is.null(alpha_manual))) {
	    stopifnot(length(unique(label_scheme_sub[[col_alpha]])) == length(alpha_manual))
	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + scale_shape_manual(values = alpha_manual)
	      }
	    }
	  }

	  if (show_ellipses) {
	    if (anyNA(label_scheme_sub[[col_group]])) {
	      warning("(Partial) NA aesthetics under column `", col_group, "` in expt_smry.xlsx",
	              call. = FALSE)
	    }

	    for (x in 2:dimension) {
	      for (y in 1:(x-1)) {
	        p[x, y] <- p[x, y] + stat_ellipse(
	          data = df,
	          aes(x = !!rlang::sym(paste0("LD", y)),
	              y = !!rlang::sym(paste0("LD", x)),
	              fill = !!rlang::sym(col_group)),
	          geom = "polygon",
	          alpha = .4,
	          show.legend = FALSE,
	        )
	      }
	    }
	  }
	} else if (dimension == 2) {
	  p <- ggplot() +
	    rlang::eval_tidy(rlang::quo(geom_point(data = df, 
	                                           mapping = mapping_var, 
	                                           !!!fix_args))) +
	    coord_fixed()
	  
	  check_ggplot_aes(p)

	  if (show_ellipses) {
	    if (anyNA(label_scheme_sub[[col_group]])) {
	      warning("(Partial) NA aesthetics under column `", col_group, 
	              "` in expt_smry.xlsx",
	              call. = FALSE)
	    }

	    p <- p + ggplot2::stat_ellipse(
	      data = df,
	      aes(x = LD1, y = LD2, fill = !!rlang::sym(col_group)),
	      geom = "polygon",
	      alpha = .4,
	      show.legend = FALSE,
	    )
	  }

	  if (show_ids) {
	    p <- p +
	      geom_text(data = df,
	                mapping = aes(x = LD1, y = LD2, label = Sample_ID),
	                color = "gray", size = 3)
	  }

	  p <- p +
	    labs(title = "", x = col_labs[1], y = col_labs[2]) + theme

  	if (!is.null(fill_brewer)) p <- p + scale_fill_brewer(palette = fill_brewer)
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
	} else if (dimension == 1) {
	  p <- ggplot(data = df, aes(x = Label, y = LD1)) +
	    geom_bar(stat="identity") +
	    theme(axis.text.x  = element_text(angle=60, vjust=0.5, size=8))
	}

	ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot"))
	rlang::eval_tidy(rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(filename)),
	                                   plot = p, 
	                                   !!!ggsave_dots)))

	invisible(res)
}


#' Scores LDA
#'
#' @inheritParams prnPCA
#' @inheritParams scoreMDS
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr 
#' @importFrom MASS lda
#' @importFrom magrittr %>% %T>% %$% %<>% 
scoreLDA <- function (df, id, label_scheme_sub, anal_type, scale_log2r,
                      center_features, scale_features,
                      choice = "lda", method = "moment",
                      type, col_group,
                      folds, out_file, ...) {
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  col_group <- rlang::enexpr(col_group)

  if (rlang::as_string(col_group) == "Select") {
    stop("Specify a `col_group` column other than `Select.`", call. = FALSE)
  }

  df_orig <- df

  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r,
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type, 
               rm_allna = TRUE) %>%
    .$log2R

  nms <- names(df)
  n_rows <- nrow(df)

  if (n_rows <= 50) {
    stop("Need 50 or more data rows for LDA.", call. = FALSE)
  }

  label_scheme_sub <- label_scheme_sub %>% dplyr::filter(Sample_ID %in% nms)

  if (type == "obs") {
    res <- prep_folded_tdata(df, folds, label_scheme_sub, !!col_group)
    df_t <- res$df_t
    ls_sub <- res$ls_sub
    rm(res)

    lda_out <- local({
      message("Class labels according to the column `", col_group, 
              "` in metadata file.")

      dm_t <- df_t %>% dplyr::select(-!!rlang::sym(col_group))
      if (scale_features) {
        dm_t <- dm_t %>% scale(center = center_features, scale = TRUE)
      } else if (center_features) {
        dm_t <- dm_t %>% sweep(., 2, colMeans(., na.rm = TRUE), "-")
      }

      grps <- df_t[[col_group]]

      if (choice == "lda") {
        rlang::expr(MASS::lda(x = !!dm_t, grouping = !!grps, !!!dots)) %>%
          rlang::eval_bare(env = caller_env()) %>%
          predict(df_t %>% dplyr::select(-!!rlang::sym(col_group)))
      }
    })

    lda_out$x <- lda_out$x %>%
      data.frame(check.names = FALSE) %>%
      cmbn_meta(ls_sub) %T>%
      readr::write_tsv(out_file)

  } else if (type == "feats") {
    stop("LDA by features (proteins/peptides) not currently available.", call. = FALSE)

    if (folds > 1) {
      message("Coerce to `folds = 1` at `type = feats`.")
    }
  } else {
    stop("Unkown `type` for LDA.", call. = FALSE)
  }

  if (names(lda_out$x) %>% .[grepl("^LD\\d+$", .)] %>% length() == 1) {
    warning("2-D LDA plots not available with only two classes under column `",
            col_group, "` in expt_smry.xlsx.", call. = FALSE)
  }

  invisible(lda_out)
}


#'LDA plots
#'
#'\code{pepLDA} visualizes the linear discriminant analysis (LDA) of peptide \code{log2FC}.
#'
#'@rdname prnLDA
#'
#'@import purrr
#'@export
pepLDA <- function (col_select = NULL, col_group = NULL, col_color = NULL,
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL,
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                    choice = c("lda"), 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                    center_features = TRUE, scale_features = TRUE,
                    show_ids = TRUE, show_ellipses = FALSE,
                    type = c("obs", "feats"), 
                    method = c("moment", "mle", "mve"), 
                    dimension = 2, folds = 1,
                    df = NULL, filepath = NULL, filename = NULL,
                    theme = NULL, 
                    
                    formula = NULL, data = NULL, 
                    x = NULL, grouping = NULL, 
                    prior = NULL, subset = NULL, CV = NULL, 
                    na.action = NULL, nu = NULL, 
                    ...) {

  old_opts <- options()
  options(warn = 1, warnPartialMatchArgs = TRUE)
  on.exit(options(old_opts), add = TRUE)
  
  check_dots(c("id", "df2", "anal_type"), ...)

  choice <- rlang::enexpr(choice)
  if (length(choice) > 1) {
    choice <- "lda"
  } else {
    choice <- rlang::as_string(choice)
  }
  
  method <- rlang::enexpr(method)
  if (length(method) > 1) {
    method <- "moment"
  } else {
    method <- rlang::as_string(method)
  }
  
  type <- rlang::enexpr(type)
  if (length(type) > 1) {
    type <- "obs"
  } else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("obs", "feats"))
  }  
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1)

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
            anal_type = "LDA")(choice = choice, 
                               method = method,
                               type = type, 
                               dimension = dimension, 
                               folds = folds, 
                               show_ids = show_ids,
                               show_ellipses = show_ellipses,
                               center_features = center_features,
                               scale_features = scale_features,
                               theme = theme, 
                               
                               formula = formula, 
                               data = data, 
                               x = x, 
                               grouping = grouping, 
                               prior = prior, 
                               subset = subset, 
                               CV = CV, 
                               na.action = na.action, 
                               nu = nu, 
                               ...)
}


#'LDA plots
#'
#'\code{prnLDA} visualizes the linear discriminant analysis (LDA) of protein
#'\code{log2FC}.
#'
#'The utility is a wrapper of \code{\link[MASS]{lda}}.
#'
#'@param choice Character string; the LDA method in one of \code{c("lda")}. The
#'  default is "lda".
#'@inheritParams  prnHist
#'@inheritParams prnHM
#'@inheritParams anal_prnNMF
#'@inheritParams prnPCA
#'@param formula Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param data Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param x Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param grouping Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param prior Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param subset Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param CV Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param na.action Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
#'@param nu Dummy argument to avoid incurring the corresponding argument in
#'  a pre-existed function by partial argument matches.
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
#'@example inst/extdata/examples/prnLDA_.R
#'
#'@return LDA plots.
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
prnLDA <- function (col_select = NULL, col_group = NULL, col_color = NULL,
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL,
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                    choice = c("lda"), 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                    center_features = TRUE, scale_features = TRUE,
                    show_ids = TRUE, show_ellipses = FALSE,
                    type = c("obs", "feats"), 
                    method = c("moment", "mle", "mve"), 
                    dimension = 2, folds = 1,
                    df = NULL, filepath = NULL, filename = NULL,
                    theme = NULL, 
                    
                    formula = NULL, data = NULL, 
                    x = NULL, grouping = NULL, 
                    prior = NULL, subset = NULL, CV = NULL, 
                    na.action = NULL, nu = NULL, 
                    ...) {

  old_opts <- options()
  options(warn = 1, warnPartialMatchArgs = TRUE)
  on.exit(options(old_opts), add = TRUE)
  
  check_dots(c("id", "df2", "anal_type"), ...)

  choice <- rlang::enexpr(choice)
  if (length(choice) > 1) {
    choice <- "lda"
  } else {
    choice <- rlang::as_string(choice)
  }
  
  method <- rlang::enexpr(method)
  if (length(method) > 1) {
    method <- "moment"
  } else {
    method <- rlang::as_string(method)
  }
  
  type <- rlang::enexpr(type)
  if (length(type) > 1) {
    type <- "obs"
  } else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("obs", "feats"))
  }  
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1)

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
            anal_type = "LDA")(choice = choice, 
                               method = method,
                               type = type, 
                               dimension = dimension, 
                               folds = folds, 
                               show_ids = show_ids,
                               show_ellipses = show_ellipses,
                               center_features = center_features,
                               scale_features = scale_features,
                               theme = theme, 
                               
                               formula = formula, 
                               data = data, 
                               x = x, 
                               grouping = grouping, 
                               prior = prior, 
                               subset = subset, 
                               CV = CV, 
                               na.action = na.action, 
                               nu = nu, 
                               ...)
}


