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
                       size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, col_benchmark = NULL, 
                       scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                       df = NULL, filepath = NULL, filename = NULL,
                       anal_type = c("Corrplot", "Heatmap", "Histogram", "MA", "MDS", "Model",
                                     "NMF", "Trend")) {

  stopifnot(rlang::is_logical(scale_log2r), 
            rlang::is_logical(complete_cases),
            rlang::is_logical(impute_na))
  
  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")
  warn_msg1 <- "Coerce `complete_cases = TRUE` at `impute_na = FALSE`."
  
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
	  # warning("No samples under column \'", rlang::as_string(col_order), "\'.", call. = FALSE)
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
		data_type <- "Protein"
	} else if (id %in% c("pep_seq", "pep_seq_mod")) {
		data_type <- "Peptide"
	} else {
		stop("Unrecognized 'id'; needs to be in c(\"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\")", 
		     call. = TRUE)
	}

	anal_type <- rlang::as_string(rlang::enexpr(anal_type))
	
	if (is.null(filepath)) {
		if (grepl("Trend", anal_type)) {
		  filepath = file.path(dat_dir, data_type, "Trend")
		} else if (grepl("NMF", anal_type)) {
		  filepath = file.path(dat_dir, data_type, "NMF")
		} else if (grepl("GSPA", anal_type)) {
		  filepath = file.path(dat_dir, data_type, "GSPA")
		} else {
		  filepath = file.path(dat_dir, data_type, anal_type)
		}
		dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
	} else {
	  stop("Use default `filepath`.", call. = FALSE)
	}

	if (is.null(filename)) {
		fn_prefix <- paste(data_type, anal_type, sep = "_")
		fn_prefix <- paste0(fn_prefix, "_", ifelse(scale_log2r, "Z", "N"))
		fn_prefix <- fn_prefix %>% ifelse(impute_na, paste0(., "_impNA"), .)

		fn_suffix <- ifelse(anal_type == "Model", "txt", "png")
	} else {
	  if (length(filename) > 1) stop("Do not provide multiple file names.", call. = FALSE)
	  fn_prefix <- gsub("\\.[^.]*$", "", filename)
	  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
		if (fn_prefix == fn_suffix) stop("No '.' in the file name.", call. = FALSE)
	  
	  if (anal_type %in% c("Trend", "NMF", "GSPA")) {
	    fn_prefix <- paste(fn_prefix, data_type, anal_type, sep = "_")
	    fn_prefix <- paste0(fn_prefix, "_", ifelse(scale_log2r, "Z", "N"))
	    fn_prefix <- fn_prefix %>% ifelse(impute_na, paste0(., "_impNA"), .)
	  }
	}

	use_pri_data <- c("MDS", "PCA", "EucDist", "Heatmap", "Histogram", "Corrplot", 
	                  "Model", "Trend", "NMF", "NMF_meta", "GSPA", "GSVA", "GSEA", "String")
	use_sec_data <- c("Trend_line", "NMF_con", "NMF_coef", "GSPA_hm")

	if (anal_type %in% use_pri_data) {
	  df <- find_pri_df(anal_type = !!anal_type, id = !!id, impute_na = impute_na)
	} else if (anal_type %in% use_sec_data) {
	  df <- find_sec_df(df = !!df, anal_type = !!anal_type, id = !!id)
	} else {
	  stop("Unknown `anal_type`", call. = FALSE)
	}

	if (!is.null(dim(df))) df <- df %>% rm_pval_whitespace()

	if (anal_type %in% c("Model", "GSPA", "GSVA", "GSEA")) {
		label_scheme_sub <- label_scheme %>% # to be subset by "formulae"
			dplyr::filter(!grepl("^Empty\\.[0-9]+", .$Sample_ID), !Reference)
	} else {
		label_scheme_sub <- label_scheme %>%
			dplyr::select(Sample_ID, TMT_Set, !!col_select, !!col_group, !!col_order, !!col_color,
			              !!col_fill, !!col_shape, !!col_size, !!col_alpha, !!col_benchmark) %>%
			dplyr::filter(!is.na(!!col_select))
	}

	if (nrow(label_scheme_sub) == 0)
	  stop(paste0("No samples or conditions were defined for \"", anal_type, "\""))
	
	force(scale_log2r)
	force(impute_na)
	force(complete_cases)
	
	message("scale_log2r = ", scale_log2r)
	message("impute_na = ", impute_na)
	message("complete_cases = ", complete_cases)
	
	## rules for anal_ functions
	# `impute_na` determines the `src_path` for `df`
	# `scale_log2r` determines `N_log2R` or `Z_log2R` columns in `df`
	# `complete_cases` subsets data rows in `df` for samples in `label_scheme_sub`
	
	## special-cases
	# `GSPA`: use `pVals` instead of `N_log2R` or `Z_log2R`; 
	#    therefore `scale_log2r` only for output file names
	#    `complete_cases` for `entrez` and `pVals` so actually has no effect
	# `Model` and `GSVA`: `complete_cases` depends on lm formulae

	## rules for plot_ functions
	# `scale_log2r` for matching '_N' or '_Z' in input filenames from prior anal_ functions
	# `impute_na` for matching '[_NZ]_impNA' or '[_NZ]' in input filenames

	if (anal_type == "MDS") {
		function(adjEucDist = FALSE, classical = TRUE, k = 3, show_ids = TRUE, annot_cols = NULL, ...) {
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
		            show_ids = show_ids, 
		            !!!dots)
		}
	} else if (anal_type == "PCA") {
		function(type = obs, show_ids = TRUE, annot_cols = NULL, ...) {
		  dots <- rlang::enexprs(...)
		  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
		  dots <- dots %>% .[! . %in% filter_dots]
		  
		  # remove `warn_msg1` after moving the following codes to `plotPCA`
		  
		  if (!(impute_na || complete_cases)) {
		    complete_cases <- TRUE
		    rlang::warn(warn_msg1)
		  }
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
		function(adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, ...) {
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
		function(xmin = -1, xmax = 1, xmargin = 0.1,
		         annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, ...) {
		  # if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
		  
		  plotHM(df = df, 
		         id = !!id, 
		         col_benchmark = !!col_benchmark,
			       label_scheme_sub = label_scheme_sub,
			       filepath = filepath, 
			       filename = paste0(fn_prefix, ".", fn_suffix),
			       scale_log2r = scale_log2r, 
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
		  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
		  
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
		  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
		  
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
	} else if (anal_type == "Model") {
	  function(method = "limma", var_cutoff = 1E-3, pval_cutoff = 1, logFC_cutoff = log2(1), ...) {
	    method <- rlang::enexpr(method)
	    
	    dots <- rlang::enexprs(...)
	    filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	    dots <- dots %>% .[! . %in% filter_dots]
	    
	    fn_prefix2 <- ifelse(impute_na, "_impNA_pVals.txt", "_pVals.txt")
	    
	    # `complete_cases` depends on lm contrasts
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
	    
	    # record the `scale_log2r` status; otherwise, need to indicate it in a way
	    # for example, `_N` or `_Z` in file names
	    dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)
	    local({
	      if (data_type == "Peptide") type <- "pep" else if (data_type == "Protein") type <- "prn"
	      file <- paste0(type, "Sig_imp", ifelse(impute_na, "TRUE", "FALSE"), ".rda")
	      
	      call_pars <- c(scale_log2r = scale_log2r, 
	                     complete_cases = complete_cases, 
	                     impute_na = impute_na) %>% 
	        as.list()
	      
	      save(call_pars, file = file.path(dat_dir, "Calls", file))
	    })
	    
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
	} else if (anal_type == "Trend") {
		function(n_clust = NULL, ...) {
		  analTrend(df = df, 
		            id = !!id, 
                col_group = !!col_group, 
                col_order = !!col_order,
                label_scheme_sub = label_scheme_sub, 
                n_clust = n_clust,
		            scale_log2r = scale_log2r,
		            complete_cases = complete_cases, 
		            impute_na = impute_na,
                filepath = filepath, 
                filename = paste0(fn_prefix %>% paste0(., "_nclust", n_clust), ".txt"), 
		            anal_type = anal_type,
                ...)
		}
	} else if (anal_type == "Trend_line") {
	  function(n_clust = NULL, theme = NULL, ...) {
	    plotTrend(id = !!id, 
	              col_group = !!col_group, 
	              col_order = !!col_order, 
	              label_scheme_sub = label_scheme_sub, 
	              n_clust = n_clust, 
	              scale_log2r = scale_log2r,
	              complete_cases = complete_cases, 
	              impute_na = impute_na,
	              filepath = filepath, 
	              filename = paste0(fn_prefix, ".", fn_suffix), 
	              theme = theme, 
	              ...)
	  }
	} else if (anal_type == "NMF") {
		function(rank = NULL, nrun = 50, seed = NULL, ...) {
		  analNMF(df = df, 
		          id = !!id, 
		          rank = rank, 
		          nrun = nrun, 
		          seed = seed, 
		          col_group = !!col_group, 
		          label_scheme_sub = label_scheme_sub,
		          scale_log2r = scale_log2r,
		          complete_cases = complete_cases, 
		          impute_na = impute_na,
		          filepath = filepath, 
		          filename = paste0(fn_prefix %>% paste0(., "_rank", rank), ".txt"), 
		          anal_type = anal_type,
		          ...)
		}
	} else if (anal_type == "NMF_con") {
	  function(rank = NULL, ...) {
	    plotNMFCon(id = !!id, 
	               rank = rank, 
	               label_scheme_sub = label_scheme_sub, 
	               scale_log2r = scale_log2r, 
	               complete_cases = complete_cases, 
	               impute_na = impute_na,
	               filepath = filepath, 
	               filename = paste0(fn_prefix, ".", fn_suffix), 
	               ...)
	  }
	} else if (anal_type == "NMF_coef") {
	  function(rank = NULL, ...) {
	    plotNMFCoef(id = !!id, 
	                rank = rank, 
	                label_scheme_sub = label_scheme_sub, 
	                scale_log2r = scale_log2r,
	                complete_cases = complete_cases,
	                impute_na = impute_na,                
	                filepath = filepath, 
	                filename = paste0(fn_prefix, ".", fn_suffix), 
	                ...)
	  }
	} else if (anal_type == "NMF_meta") {
	  function(rank = NULL, ...) {
	    plotNMFmeta(df = df, 
	                id = !!id, 
	                rank = rank, 
	                label_scheme_sub = label_scheme_sub, 
	                scale_log2r = scale_log2r,
	                complete_cases = complete_cases, 
	                impute_na = impute_na, 
	                filepath = filepath, 
	                filename = paste0(fn_prefix, ".", fn_suffix), 
	                anal_type = anal_type,
	                ...)
	  }
	} else if (anal_type == "GSPA") {
		function(gset_nms = "go_sets", var_cutoff = .5,
		         pval_cutoff = 1E-2, logFC_cutoff = log2(1.1), gspval_cutoff = 1E-2, 
		         min_size = 10, max_size = Inf, min_greedy_size = 1, ...) {

		  # currently Protein[_impNA]_pVals.txt does not contain `scale_log2_r` info in the file name;
		  #   therefore, `scale_log2r` will be matched to those in `prnSig()` and indicated in 
		  #   file names `_N[_impNA].txt` or `_Z[_impNA].txt` to inform user the `scale_log_r` status
		  # `complete_cases` is for entrez IDs and pVals, so typically has no effect
		  # "id" only for tibbling rownames
		  gspaTest(df = df, 
		           id = !!id, 
		           label_scheme_sub = label_scheme_sub, 
		           scale_log2r = scale_log2r, 
		           complete_cases = complete_cases, 
		           impute_na = impute_na, 
		           filepath = filepath, 
		           filename = paste0(fn_prefix, ".txt"),
		           gset_nms = gset_nms, 
		           var_cutoff = var_cutoff,
		           pval_cutoff = pval_cutoff, 
		           logFC_cutoff = logFC_cutoff, 
		           gspval_cutoff = gspval_cutoff,
		           min_size = min_size, 
		           max_size = max_size, 
		           min_greedy_size = min_greedy_size, 
		           anal_type = anal_type, 
		           ...)
		}
	} else if (anal_type == "GSPA_hm") {
	  function(...) {
	    
	    # `scale_log2r` to match file names from `gspaTest()`
	    # `impute_na` to match file names from `gspaTest()`
	    # `complete_cases` has no effects
	    gspaHM(scale_log2r = scale_log2r, 
	           complete_cases = complete_cases, 
	           impute_na = impute_na, 
	           filepath = filepath, 
	           filename = paste0(fn_prefix, ".", fn_suffix), 
	           ...)
	  }
	} else if(anal_type == "GSVA") {
	  function(lm_method = "limma", gset_nms = "go_sets", var_cutoff = .5,
	           pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), ...) {
      gsvaTest(df = df, 
               id = !!id, 
               label_scheme_sub = label_scheme_sub, 
               filepath = filepath, 
               filename = paste0(fn_prefix, ".txt"), 
               scale_log2r = scale_log2r, 
               complete_cases = complete_cases, 
               impute_na = impute_na, 
               gset_nms = gset_nms,
               lm_method = lm_method, 
               gsets = gsets, 
               var_cutoff = var_cutoff, 
               pval_cutoff = pval_cutoff, 
               logFC_cutoff = logFC_cutoff, 
               anal_type = anal_type, 
               ...)
	  }
	} else if (anal_type == "GSEA") {
	  function(gset_nms = "go_sets", var_cutoff = .5,
	           pval_cutoff = 1E-2, logFC_cutoff = log2(1.1), gspval_cutoff = 1E-2, 
	           min_size = 10, max_size = Inf, min_greedy_size = 1, task = anal, ...) {
	    
	    # DO NOT go through lm; therefore, `complete_cases` can be applied to 
	    #   all samples in label_scheme_sub
	    if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

	    switch(rlang::as_string(rlang::enexpr(task)), 
	           anal = gspaTest(df = df, 
	                           id = !!id, 
	                           label_scheme_sub = label_scheme_sub,
	                           filepath = filepath, 
	                           filename = paste0(fn_prefix, ".txt"),
	                           complete_cases = complete_cases, # redundant
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
	           stop("Invalid `task`.", Call. = TRUE)
	    )
	  }
	} else if (anal_type == "String") {
	  function(db_path = "~\\proteoQ\\dbs\\string", score_cutoff = .7, adjP = FALSE, task = anal, ...) {
	    dots <- rlang::enexprs(...)
	    filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	    dots <- dots %>% .[! . %in% filter_dots]
	    
	    if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
	    
	    df <- df %>% filters_in_call(!!!filter_dots)

	    # `scale_log2r` not used; both `_N` and `_Z` columns will be kept
	    switch(rlang::as_string(rlang::enexpr(task)), 
	           anal = stringTest(df = df, 
	                             id = !!id, 
	                             col_group = !!col_group, 
	                             col_order = !!col_order,
	                             label_scheme_sub = label_scheme_sub, 
	                             db_path = db_path,
	                             score_cutoff = score_cutoff,
	                             adjP = adjP, 
	                             scale_log2r = scale_log2r, 
	                             filepath = filepath, 
	                             filename = paste0(fn_prefix, ".csv"), 
	                             !!!dots), 
	           stop("Invalid `task`.", Call. = TRUE)
	    )
	  }
	} 
	
}



#' Helper for finding input `df`
#'
find_pri_df <- function (anal_type = "Model", df = NULL, id = "gene", impute_na = FALSE, ...) {
  err_msg2 <- "not found. \n Run functions PSM, peptide and protein normalization first."
  err_msg3 <- "not found. \nImpute NA values with `pepImp()` or `prnImp()` or set `impute_na = FALSE`."
  err_msg4 <- "not found at impute_na = TRUE. \nRun `prnSig(impute_na = TRUE)` first."
  err_msg5 <- "not found at impute_na = FALSE. \nRun `prnSig(impute_na = FALSE)` first."

  anal_type <- rlang::as_string(rlang::enexpr(anal_type))
  id <- rlang::as_string(rlang::enexpr(id))

  df <- rlang::enexpr(df)

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
    } else {
      stop("Unknown `id`.", call. = FALSE)
    }
    
    if (anal_type %in% c("Histogram", "MA")) { # never impute_na
      if (file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, err_msg2), call. = FALSE)
    } else if (anal_type %in% c("Model")) { # optional impute_na but no pVals
      if (impute_na) {
        if (file.exists(fn_imp)) src_path <- fn_imp else stop(paste(fn_imp, err_msg3), call. = FALSE)
      } else {
        if (file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, err_msg2), call. = FALSE)
      }
    } else if (anal_type %in% c("GSPA", "GSEA", "String")) { # always use data with pVals
      if (impute_na) {
        if (file.exists(fn_imp_p)) src_path <- fn_imp_p else stop(paste(fn_imp_p, err_msg4), call. = FALSE)
      } else {
        if (file.exists(fn_p)) src_path <- fn_p else stop(paste(fn_p, err_msg5), call. = FALSE)
      }
    } else if (anal_type %in% c("Heatmap", "MDS", "PCA", "EucDist", "Trend", "NMF", "NMF_meta", 
                                "GSVA", "Corrplot")) { # optional impute_na and possible p_vals
      if (impute_na) {
        if (file.exists(fn_imp_p)) src_path <- fn_imp_p else if (file.exists(fn_imp)) src_path <- fn_imp else 
          stop(paste(fn_imp, err_msg3), call. = FALSE)
      } else {
        if (file.exists(fn_p)) src_path <- fn_p else if (file.exists(fn_raw)) src_path <- fn_raw else 
          stop(paste(fn_raw, err_msg2), call. = FALSE)
      }
    } 
  } else {
    df <- rlang::as_string(df)
    
    if (anal_type == "Model") stop("Use default file name.", call. = FALSE)
    
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      src_path <- file.path(dat_dir, "Peptide\\Model", df)
      if (!file.exists(src_path)) src_path <- file.path(dat_dir, "Peptide", df)
    } else if (id %in% c("prot_acc", "gene")) {
      src_path <- file.path(dat_dir, "Protein\\Model", df)
      if (!file.exists(src_path)) src_path <- file.path(dat_dir, "Protein", df)
    }
  }
  
  df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
                          comment.char = "#"), error = function(e) NA)

  if (is.null(dim(df))) stop(src_path, " not found.", call. = FALSE)
  
  message(paste("File loaded:", gsub("\\\\", "/", src_path)))
  
  return(df)
}


#' Helper for finding input `df`
#'
find_sec_df <- function (df = NULL, anal_type = NULL, id = NULL, ...) {
  df <- rlang::enexpr(df)
  anal_type <- rlang::enexpr(anal_type)
  id <- rlang::enexpr(id)

  if (is.null(df) || is.null(anal_type) || is.null(id)) return (NULL)
  
  df <- rlang::as_string(df)
  anal_type <- rlang::as_string(anal_type)
  id <- rlang::as_string(id)

  new_anal_type <- anal_type %>% 
    gsub("^Trend_.*", "Trend", .) %>% 
    gsub("^NMF_.*", "NMF", .) %>% 
    gsub("^GSPA_.*", "GSPA", .)
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    src_path <- file.path(dat_dir, "Peptide", new_anal_type, df)
  } else if (id %in% c("prot_acc", "gene")) {
    src_path <- file.path(dat_dir, "Protein", new_anal_type, df)
  } else {
    stop("Unknown `id`", call. = FALSE)
  }
  
  df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
                          comment.char = "#"), error = function(e) NA)
  
  if (is.null(dim(df))) stop(src_path, " not found.", call. = FALSE)
  
  return(df)
}