#' Function factories for informatic analysis
#'
#' \code{info_anal} produces functions for selected informatic analysis.
#'
#' @param anal_type Character string; the type of analysis that are preset for
#'   method dispatch in function factories. The value will be determined
#'   automatically. Examplary values include \code{anal_type = c("PCA",
#'   "Corrplot", "EucDist", "GSPA", "Heatmap", "Histogram", "MDS", "Model",
#'   "NMF", "Purge", "Trend", ...)}.
#' @return a function to the given \code{anal_type}.
#'
#' @import dplyr rlang ggplot2 pheatmap openxlsx
#' @importFrom magrittr %>%
info_anal <- function (id = gene, col_select = NULL, col_group = NULL, col_order = NULL,
                       col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                       color_brewer = NULL, fill_brewer = NULL, 
                       size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                       col_benchmark = NULL, scale_log2r = TRUE,
                       impute_na = FALSE, df = NULL, filepath = NULL, filename = NULL,
                       anal_type = c("Corrplot", "Heatmap", "Histogram", "MA", "MDS", "Model",
                                     "NMF", "Trend")) {
  
  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")

  old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)
	
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_order <- rlang::enexpr(col_order)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	col_benchmark <- rlang::enexpr(col_benchmark)

	col_select <- ifelse(is.null(col_select), rlang::expr(Select), rlang::sym(col_select))
	col_group <- ifelse(is.null(col_group), rlang::expr(Group), rlang::sym(col_group))
	col_order <- ifelse(is.null(col_order), rlang::expr(Order), rlang::sym(col_order))
	col_color <- ifelse(is.null(col_color), rlang::expr(Color), rlang::sym(col_color))
	col_fill <- ifelse(is.null(col_fill), rlang::expr(Fill), rlang::sym(col_fill))
	col_shape <- ifelse(is.null(col_shape), rlang::expr(Shape), rlang::sym(col_shape))
	col_size <- ifelse(is.null(col_size), rlang::expr(Size), rlang::sym(col_size))
	col_alpha <- ifelse(is.null(col_alpha), rlang::expr(Alpha), rlang::sym(col_alpha))
	col_benchmark <- ifelse(is.null(col_benchmark), rlang::expr(Benchmark), rlang::sym(col_benchmark))
	
	if (col_select == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_group == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_order == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_color == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_fill == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_shape == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_size == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_alpha == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_benchmark == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	
	color_brewer <- rlang::enexpr(color_brewer)
	fill_brewer <- rlang::enexpr(fill_brewer)
	if (!is.null(color_brewer)) color_brewer <- rlang::as_string(color_brewer)
	if (!is.null(fill_brewer)) fill_brewer <- rlang::as_string(fill_brewer)
	
	size_manual <- rlang::enexpr(size_manual)
	shape_manual <- rlang::enexpr(shape_manual)
	alpha_manual <- rlang::enexpr(alpha_manual)

	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)

	load(file = file.path(dat_dir, "label_scheme.rda"))
	
	if (is.null(label_scheme[[col_select]])) {
		stop("Column \'", rlang::as_string(col_select), "\' not found.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
		stop("No samples under column \'", rlang::as_string(col_select), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_group]])) {
		col_group <- rlang::expr(Select)
		warning("Column \'", rlang::as_string(col_group), "\' not found; use column \'Select\'.", call. = FALSE)
			
	} else if (sum(!is.na(label_scheme[[col_group]])) == 0) {
		col_group <- rlang::expr(Select)
		warning("No samples under \'", rlang::as_string(col_group), "\'; use column \'Select\'.", call. = FALSE)
			
	}

	if (is.null(label_scheme[[col_order]])) {
		warning("Column \'", rlang::as_string(col_order), "\' not found; arranged by the alphebatics.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_order]])) == 0) {
	  warning("No samples under column \'", rlang::as_string(col_order), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_color]])) {
		warning("Column \'", rlang::as_string(col_color), "\' not found.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_color]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_color), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_fill]])) {
		warning("Column \'", rlang::as_string(col_fill), "\' not found.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_fill]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_fill), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_shape]])) {
		warning("Column \'", rlang::as_string(col_shape), "\' not found.")
	} else if(sum(!is.na(label_scheme[[col_shape]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_shape), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_size]])) {
		warning("Column \'", rlang::as_string(col_size), "\' not found.")
	} else if (sum(!is.na(label_scheme[[col_size]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_size), "\'.", call. = FALSE)
	}

	if(is.null(label_scheme[[col_alpha]])) {
		warning("Column \'", rlang::as_string(col_alpha), "\' not found.")
	} else if(sum(!is.na(label_scheme[[col_alpha]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_alpha), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_benchmark]])) {
		warning("Column \'", rlang::as_string(col_benchmark), "\' not found.")
	} else if (sum(!is.na(label_scheme[[col_benchmark]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_benchmark), "\'.", call. = FALSE)
	}

	id <- rlang::as_string(rlang::enexpr(id))
	if (length(id) != 1)
	  stop("\'id\' must be one of \'pep_seq\', \'pep_seq_mod\', \'prot_acc\' or \'gene\'")

	if (id %in% c("prot_acc", "gene")) {
	  # cat(paste0("id = \"", id, "\"", " after parameter matching to normPrn()\n"))
		data_type <- "Protein"
	} else if (id %in% c("pep_seq", "pep_seq_mod")) {
	  # cat(paste0("id = \"", id, "\"", " after parameter matching to normPep()\n"))
		data_type <- "Peptide"
	} else {
		stop("Unrecognized 'id';
		     needs to be in c(\"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\")", call. = TRUE)
	}

	anal_type <- rlang::as_string(rlang::enexpr(anal_type))
	
	if (is.null(filepath)) {
		filepath = file.path(dat_dir, data_type, anal_type)
		dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
	}

	if (is.null(filename)) {
		fn_prefix <- paste(data_type, anal_type, sep = "_")
		fn_prefix <- paste0(fn_prefix, "_", ifelse(scale_log2r, "Z", "N"))
		fn_prefix <- fn_prefix %>% ifelse(impute_na, paste0(., "_impNA"), .)
		fn_suffix <- "png"
	} else {
		fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
		fn_prefix <- gsub("\\.[^.]*$", "", filename)
	}
	
	err_msg2 <- "not found. \n Run functions `normPSM()`, `normPep()` and `normPrn()` first."
	err_msg3 <- "not found. \nImpute NA values with `pepImp()` or `prnImp()` or set `impute_na = FALSE`."
	err_msg4 <- "not found at impute_na = TRUE. \nRun `prnSig(impute_na = TRUE)` first."
	err_msg5 <- "not found at impute_na = FALSE. \nRun `prnSig(impute_na = FALSE)` first."
	
	if (is.null(df)) {
	  if (id %in% c("pep_seq", "pep_seq_mod")) {
	    fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
	    fn_imp_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_impNA_pVals.txt")
	    fn_raw <- file.path(dat_dir, "Peptide", "Peptide.txt")
	    fn_imp <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
	  } else if (id %in% c("prot_acc", "gene")) {
	    fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
	    fn_imp_p <- file.path(dat_dir, "Protein\\Model", "Protein_impNA_pVals.txt")
	    fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")
	    fn_imp <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
	  } else stop("Unknown `id`.", call. = FALSE)
	  
	  if (anal_type %in% c("Histogram", "Corrplot", "MA")) { # never impute_na
	    if (file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, err_msg2), call. = FALSE)
	  } else if (anal_type %in% c("Model")) { # optional impute_na but no pVals
	    if (impute_na) {
	      if (file.exists(fn_imp)) src_path <- fn_imp else stop(paste(fn_imp, err_msg3), call. = FALSE)
	    } else {
	      if (file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, err_msg2), call. = FALSE)
	    }
	  } else if (anal_type %in% c("GSPA", "GSEA")) { # always use data with pVals
	    if (impute_na) {
	      if (file.exists(fn_imp_p)) src_path <- fn_imp_p else stop(paste(fn_imp_p, err_msg4), call. = FALSE)
	    } else {
	      if (file.exists(fn_p)) src_path <- fn_p else stop(paste(fn_p, err_msg5), call. = FALSE)
	    }
	  } else if (anal_type %in% c("Heatmap", "MDS", "PCA", "EucDist", "Trend", "NMF", "GSVA")) { # optional impute_na and possible p_vals
	    if (impute_na) {
	      if (file.exists(fn_imp_p)) src_path <- fn_imp_p else if (file.exists(fn_imp)) src_path <- fn_imp else 
	        stop(paste(fn_imp, err_msg3), call. = FALSE)
	    } else {
	      if (file.exists(fn_p)) src_path <- fn_p else if (file.exists(fn_raw)) src_path <- fn_raw else 
	        stop(paste(fn_raw, err_msg2), call. = FALSE)
	    }
	  } 
	  
	  df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
	                          comment.char = "#"), error = function(e) NA)
	  
	  if (!is.null(dim(df))) {
	    message(paste("File loaded:", gsub("\\\\", "/", src_path)))
	  } else {
	    stop(paste("Non-existed file or directory:", gsub("\\\\", "/", src_path)))
	  }
	} else {
	  if (id %in% c("pep_seq", "pep_seq_mod")) {
	    fn_raw <- file.path(dat_dir, "Peptide", df)
	  } else if (id %in% c("prot_acc", "gene")) {
	    fn_raw <- file.path(dat_dir, "Protein", df)
	  }
	  
	  df <- tryCatch(read.csv(fn_raw, check.names = FALSE, header = TRUE, sep = "\t",
	                          comment.char = "#"), error = function(e) NA)
	  
	  if (!is.null(dim(df))) {
	    message(paste("File loaded:", gsub("\\\\", "/", fn_raw)))
	  } else {
	    stop(paste("Non-existed file or directory:", gsub("\\\\", "/", fn_raw)))
	  }
	}
	
	df <- df %>% rm_pval_whitespace()
	
	if (anal_type %in% c("Model", "GSPA", "GSVA", "GSEA")) {
		label_scheme_sub <- label_scheme %>% # to be subset by the "key" in "formulae"
			dplyr::filter(!grepl("^Empty\\.[0-9]+", .$Sample_ID), !Reference)
	} else {
		label_scheme_sub <- label_scheme %>%
			dplyr::select(Sample_ID, TMT_Set, !!col_select, !!col_group, !!col_order, !!col_color,
			              !!col_fill, !!col_shape, !!col_size, !!col_alpha, !!col_benchmark) %>%
			dplyr::filter(!is.na(!!col_select))
	}

	if (nrow(label_scheme_sub) == 0)
	  stop(paste0("No samples or conditions were defined for \"", anal_type, "\""))
	
	if (anal_type == "MDS") {
		function(complete_cases = FALSE, 
		         adjEucDist = FALSE, classical = TRUE, k = 3, show_ids = TRUE, annot_cols = NULL, ...) {
		  dots <- rlang::enexprs(...)
		  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
		  dots <- dots %>% .[! . %in% filter_dots]

		  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

		  df %>% 
		    filters_in_call(!!!filter_dots) %>% 
		    scoreMDS(
		      id = !!id, 
		      label_scheme_sub = label_scheme_sub, 
		      anal_type = anal_type, 
		      scale_log2r = scale_log2r, 
		      adjEucDist = adjEucDist, 
		      classical = classical, 
		      k = k, 
		      !!!filter_dots) %>% 
		    cmbn_meta(label_scheme_sub) %T>% 
		    write.csv(file.path(filepath, paste0("res_", fn_prefix, ".csv")), row.names = FALSE) %>% 
		    plotMDS(col_color = !!col_color, 
		            col_fill = !!col_fill, 
		            col_shape = !!col_shape, 
		            col_size = !!col_size, 
		            col_alpha = !!col_alpha, 
		            color_brewer = !!color_brewer,
		            fill_brewer = !!fill_brewer, 
		            size_manual = size_manual,
		            shape_manual = shape_manual,
		            alpha_manual = alpha_manual, 
		            label_scheme_sub = label_scheme_sub, 
		            filepath = filepath, 
		            filename = paste0(fn_prefix, ".", fn_suffix), 
		            complete_cases = complete_cases, 
		            show_ids = show_ids, 
		            !!!dots)
		}
	} else if (anal_type == "PCA") {
		function(complete_cases = FALSE, type = obs, show_ids = TRUE, annot_cols = NULL, ...) {
		  dots <- rlang::enexprs(...)
		  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
		  dots <- dots %>% .[! . %in% filter_dots]
		  
		  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

		  type <- rlang::as_string(rlang::enexpr(type))

		  df_pca <- df %>% 
		    filters_in_call(!!!filter_dots) %>% 
		    scorePCA(
		      id = !!id, 
		      label_scheme_sub = label_scheme_sub, 
		      anal_type = anal_type, 
		      scale_log2r = scale_log2r, 
		      type = type, 
		      !!!filter_dots)
		  
		  df_pca$PCA %>% 
		    cmbn_meta(label_scheme_sub) %>% 
		    plotPCA(col_color = !!col_color, 
		            col_fill = !!col_fill, 
		            col_shape = !!col_shape, 
		            col_size = !!col_size, 
		            col_alpha = !!col_alpha, 
		            color_brewer = !!color_brewer,
		            fill_brewer = !!fill_brewer, 
		            size_manual = size_manual,
		            shape_manual = shape_manual,
		            alpha_manual = alpha_manual, 
		            label_scheme_sub = label_scheme_sub, 
		            prop_var = df_pca$prop_var, 
		            filepath = filepath, 
		            filename = paste0(fn_prefix, ".", fn_suffix), 
		            show_ids = show_ids, 
		            !!!dots)
		}
	} else if (anal_type == "EucDist") {
		function(complete_cases = FALSE, adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, ...) {
		  dots <- rlang::enexprs(...)
		  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
		  dots <- dots %>% .[! . %in% filter_dots]

		  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

		  df %>% 
		    filters_in_call(!!!filter_dots) %>% 
		    scoreEucDist(
		      id = !!id, 
		      label_scheme_sub = label_scheme_sub, 
		      anal_type = anal_type, 
		      scale_log2r = scale_log2r, 
		      adjEucDist = adjEucDist, 
		      !!!filter_dots) %>% 
		    plotEucDist(annot_cols, 
		                annot_colnames, 
		                filepath, 
		                filename = paste0(fn_prefix, ".", fn_suffix), 
		                !!!dots)
		}
	} else if (anal_type == "Heatmap") {
		function(complete_cases = FALSE, xmin = -1, xmax = 1, xmargin = 0.1,
		         annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, ...) {

		  plotHM(df = df, 
		         id = !!id, 
		         scale_log2r = scale_log2r, 
		         col_benchmark = !!col_benchmark,
			       label_scheme_sub = label_scheme_sub,
			       filepath = filepath, 
			       filename = paste0(fn_prefix, ".", fn_suffix),
			       complete_cases = complete_cases, 
			       annot_cols = annot_cols, 
			       annot_colnames = annot_colnames, 
			       annot_rows = annot_rows, 
			       xmin = xmin, 
			       xmax = xmax, 
			       xmargin = xmargin,
			       ...)
		}
	} else if (anal_type == "Histogram") {
		function(show_curves = TRUE, show_vline = TRUE, scale_y = TRUE, ...) {
			if (scale_log2r) {
				fn_par <- file.path(filepath, "MGKernel_params_Z.txt")
			} else {
			  fn_par <- file.path(filepath, "MGKernel_params_N.txt")
			}

			if (file.exists(fn_par)) {
				params <- read.csv(fn_par, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")
				params <- within(params, {mean = mean - x})
			} else {
				params <- NULL
				show_curves <- FALSE
			}

			plotHisto(df = df, 
			          id = !!id, 
			          label_scheme_sub = label_scheme_sub, 
			          params = params,
			          scale_log2r = scale_log2r, 
			          show_curves = show_curves, 
			          show_vline = show_vline, 
			          scale_y = scale_y, 
			          filepath = filepath, 
			          filename = paste0(fn_prefix, ".", fn_suffix), 
			          ...)
		}
	} else if (anal_type == "Corrplot") {
		function(data_select = "logFC", ...) {
		  plotCorr(df = df, 
		           id = !!id, 
		           anal_type = anal_type, 
		           data_select = data_select, 
		           col_select = !!col_select, 
		           col_order = !!col_order,
		           label_scheme_sub = label_scheme_sub, 
		           scale_log2r = scale_log2r, 
		           filepath = filepath, 
		           filename = paste0(fn_prefix, "_", data_select, ".", fn_suffix), 
		           ...)
		}
	} else if (anal_type == "Trend") {
		function(n_clust = NULL, complete_cases = FALSE, task = anal, ...) {
		  fn_prefix <- fn_prefix %>% 
		    ifelse(impute_na, paste0(., "_impNA"), .) %>% 
		    paste0(., "_n", n_clust)

		  dots <- rlang::enexprs(...)
		  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
		  dots <- dots %>% .[! . %in% filter_dots]

		  df <- df %>% 
		    filters_in_call(!!!filter_dots) %>% 
		    prepDM(id = !!id, scale_log2r = scale_log2r, 
		           sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
		    .$log2R
		  
		  switch(rlang::as_string(rlang::enexpr(task)), 
		         anal = trendTest(df = df, id = !!id, 
		                          col_group = !!col_group, 
		                          col_order = !!col_order,
		                          label_scheme_sub = label_scheme_sub, 
		                          n_clust = n_clust,
		                          complete_cases = complete_cases, 
		                          scale_log2r = scale_log2r,
		                          filepath = filepath, 
		                          filename = paste0(fn_prefix, ".csv"), 
		                          !!!dots), 
		         plot = plotTrend(id = !!id, 
		                          col_group = !!col_group, 
		                          col_order = !!col_order, 
		                          label_scheme_sub = label_scheme_sub, 
		                          n_clust = n_clust, 
		                          scale_log2r = scale_log2r,
		                          filepath = filepath, 
		                          filename = paste0(fn_prefix, ".", fn_suffix), 
		                          !!!dots), 
		         stop("Invalid `task`.", Call. = TRUE)
		  )
		}
	} else if (anal_type == "NMF") {
		function(r = NULL, nrun = 200, complete_cases = FALSE, task = anal, ...) {
		  fn_prefix <- fn_prefix %>% 
		    ifelse(impute_na, paste0(., "_impNA"), .) %>% 
		    paste0(., "_r", r)

		  switch(rlang::as_string(rlang::enexpr(task)), 
		         anal = nmfTest(df = df, 
		                        id = !!id, 
		                        r = r, 
		                        nrun = nrun, 
		                        col_group = !!col_group, 
		                        label_scheme_sub = label_scheme_sub,
		                        anal_type = anal_type, 
		                        scale_log2r = scale_log2r, 
		                        filepath = filepath, 
		                        filename = paste0(fn_prefix, ".csv"), 
		                        complete_cases = complete_cases, 
		                        ...), 
		         plotcon = plotNMFCon(id = !!id, 
		                              r = r, 
		                              label_scheme_sub = label_scheme_sub, 
		                              filepath = filepath, 
		                              filename = paste0(fn_prefix, ".", fn_suffix), 
		                              scale_log2r = scale_log2r, 
		                              ...), 
		         plotcoef = plotNMFCoef(id = !!id, 
		                                r = r, 
		                                label_scheme_sub = label_scheme_sub, 
		                                filepath = filepath, 
		                                filename = paste0(fn_prefix, ".", fn_suffix), 
		                                scale_log2r = scale_log2r, 
		                                ...), 
		         plotmeta = plotNMFmeta(df = df, id = !!id, r = r, 
		                                label_scheme_sub = label_scheme_sub, 
		                                anal_type = anal_type, 
		                                scale_log2r = scale_log2r, 
		                                filepath = filepath, 
		                                filename = paste0(fn_prefix, ".", fn_suffix), 
		                                ...), 
		         stop("Invalid `task`.", Call. = TRUE)
		  )
		}
	} else if (anal_type == "Model") {
		function(complete_cases = FALSE, method = "limma", var_cutoff = 1E-3, pval_cutoff = 1,
		         logFC_cutoff = log2(1), ...) {
		  
		  method <- rlang::enexpr(method)

		  dots <- rlang::enexprs(...)
		  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
		  dots <- dots %>% .[! . %in% filter_dots]
		  
		  fn_prefix2 <- ifelse(impute_na, "_impNA_pVals.txt", "_pVals.txt")

		  df_op <- df %>% 
		    filters_in_call(!!!filter_dots) %>% 
		    prepDM(id = !!id, 
		           scale_log2r = scale_log2r, 
		           sub_grp = label_scheme_sub$Sample_ID, 
		           anal_type = anal_type) %>% 
		    .$log2R %>% 
		    sigTest(id = !!id, 
		            label_scheme_sub = label_scheme_sub,
		            filepath = filepath, 
		            filename = paste0(fn_prefix, ".", fn_suffix),
		            complete_cases = complete_cases, 
		            method = !!method, 
		            var_cutoff = var_cutoff,
		            pval_cutoff = pval_cutoff, 
		            logFC_cutoff = logFC_cutoff, 
		            !!!dots)
		  
		  df_op <- df_op %>%
		    tibble::rownames_to_column(id) %>% 
		    dplyr::mutate(!!id := forcats::fct_explicit_na(!!rlang::sym(id))) %>% 
		    dplyr::right_join(df, ., by = id) %T>% 
		    write.table(file.path(filepath, paste0(data_type, fn_prefix2)), sep = "\t",
		                col.names = TRUE, row.names = FALSE)
		    
			wb <- createWorkbook("proteoQ")
			addWorksheet(wb, sheetName = "Results")
			openxlsx::writeData(wb, sheet = "Results", df_op)
			saveWorkbook(wb, file = file.path(filepath, paste0(data_type, "_pVals.xlsx")), overwrite = TRUE) 
		}
	} else if (anal_type == "GSPA") {
		function(complete_cases = FALSE, gset_nms = "go_sets", var_cutoff = .5,
		         pval_cutoff = 1E-2, logFC_cutoff = log2(1.1), gspval_cutoff = 1E-2, 
		         min_size = 10, max_size = Inf, min_greedy_size = 1, task = anal, ...) {

		  # "id" only for tibbling rownames
		  switch(rlang::as_string(rlang::enexpr(task)), 
		         anal = gspaTest(df = df, 
		                         id = !!id, 
		                         label_scheme_sub = label_scheme_sub,
		                         filepath = filepath, 
		                         filename = paste0(fn_prefix, ".csv"),
		                         complete_cases = complete_cases, 
		                         gset_nms = gset_nms, 
		                         var_cutoff = var_cutoff,
		                         pval_cutoff = pval_cutoff, 
		                         logFC_cutoff = logFC_cutoff, 
		                         gspval_cutoff = gspval_cutoff,
		                         min_size = min_size, 
		                         max_size = max_size, 
		                         min_greedy_size = min_greedy_size, 
		                         scale_log2r = scale_log2r, 
		                         anal_type = anal_type, 
		                         ...), 
		         plothm = gspaHM(filepath = filepath, 
		                         filename = paste0(fn_prefix, ".", fn_suffix), 
		                         ...), 
		         stop("Invalid `task`.", Call. = TRUE)
		  )
		}

	} else if(anal_type == "GSVA") {
	  function(complete_cases = FALSE, lm_method = "limma", gset_nms = "go_sets", var_cutoff = .5,
	           pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), ...) {
	    
	    dots <- rlang::enexprs(...)
	    filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	    dots <- dots %>% .[! . %in% filter_dots]
	    
	    species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
	    gsets <- load_dbs(gset_nms, species)
	    
	    stopifnot(length(gsets) > 0)

	    df_log2r <- df %>% 
	      filters_in_call(!!!filter_dots) %>% 
	      prepDM(id = id, 
	             scale_log2r = scale_log2r, 
	             sub_grp = label_scheme_sub$Sample_ID, 
	             anal_type = anal_type) %>% 
	      .$log2R %>% 
	      gsvaTest(id = !!id, 
	               label_scheme_sub = label_scheme_sub, 
	               filepath = filepath, 
	               filename = paste0(fn_prefix, ".csv"), 
	               complete_cases = complete_cases, 
	               lm_method = lm_method, 
	               gsets = gsets, 
	               var_cutoff = var_cutoff, 
	               pval_cutoff = pval_cutoff, 
	               logFC_cutoff = logFC_cutoff, !!!dots)
	  }
	} else if (anal_type == "GSEA") {
	  function(complete_cases = FALSE, gset_nms = "go_sets", var_cutoff = .5,
	           pval_cutoff = 1E-2, logFC_cutoff = log2(1.1), gspval_cutoff = 1E-2, 
	           min_size = 10, max_size = Inf, min_greedy_size = 1, task = anal, ...) {

	    switch(rlang::as_string(rlang::enexpr(task)), 
	           anal = gspaTest(df = df, 
	                           id = !!id, 
	                           label_scheme_sub = label_scheme_sub,
	                           filepath = filepath, 
	                           filename = paste0(fn_prefix, ".csv"),
	                           complete_cases = complete_cases, 
	                           gset_nms = gset_nms, 
	                           var_cutoff = var_cutoff,
	                           pval_cutoff = pval_cutoff, 
	                           logFC_cutoff = logFC_cutoff, 
	                           gspval_cutoff = gspval_cutoff,
	                           min_size = min_size, 
	                           max_size = max_size, 
	                           min_greedy_size = 1, 
	                           scale_log2r = scale_log2r, 
	                           anal_type = anal_type, 
	                           ...), 
	           plothm = gseaHM(filepath = filepath, 
	                           filename = paste0(fn_prefix, ".", fn_suffix), 
	                           ...), 
	           stop("Invalid `task`.", Call. = TRUE)
	    )
	  }
	  
	} 
	
}


