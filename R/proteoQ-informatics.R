#' Functional factories for informatic analysis
#'
#' \code{info_anal} produces functions for selected informatic analysis.
#'
#' @param anal_type The type of analysis; \code{anal_type = c("MDS", "PCA",
#'   "EucDist", "MA", "Heatmap", "Histogram", "Corrplot", "Trend", "NMF",
#'   "Model")}.
#' @return a function to the given \code{anal_type}.
#'
#' @import dplyr rlang ggplot2 pheatmap
#' @importFrom magrittr %>%
info_anal <- function (id = gene, col_select = NULL, col_group = NULL, col_order = NULL,
											col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL,
											col_alpha = NULL, col_benchmark = NULL, scale_log2r = TRUE,
											impute_na = FALSE, df = NULL, filepath = NULL, filename = NULL,
                      anal_type = c("Corrplot", "Heatmap", "Histogram", "MA", "MDS", "Model",
                                    "NMF", "Trend")) {
  
  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")

  old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)
	
	anal_type <- rlang::as_string(rlang::enexpr(anal_type))
	
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
	
	if(col_select == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_group == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_order == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_color == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_fill == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_shape == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_size == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_alpha == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if(col_benchmark == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)

	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	
	if(is.null(label_scheme[[col_select]])) {
		stop("Column \'", rlang::as_string(col_select), "\' does not exist.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_select]])) == 0) {
		stop("No samples were selected under column \'", rlang::as_string(col_select), "\'.",
		     call. = FALSE)
	}

	if(is.null(label_scheme[[col_group]])) {
		col_group <- rlang::expr(Select)
		warning("Column \'", rlang::as_string(col_group), "\' does not exist.
			Use column \'Select\' instead.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_group]])) == 0) {
		col_group <- rlang::expr(Select)
		warning("No samples were specified under column \'", rlang::as_string(col_group), "\'.
			Use column \'Select\' instead.", call. = FALSE)
	}

	if(is.null(label_scheme[[col_order]])) {
		warning("Column \'", rlang::as_string(col_order), "\' does not exist.
			Samples will be arranged by the alphebatic order.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_order]])) == 0) {
	  warning("No samples were specified under column \'", rlang::as_string(col_order), "\'.",
	          call. = FALSE)
	}

	if(is.null(label_scheme[[col_color]])) {
		warning("Column \'", rlang::as_string(col_color), "\' does not exist.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_color]])) == 0) {
		warning("No samples were specified under column \'", rlang::as_string(col_color), "\'.",
		        call. = FALSE)
	}

	if(is.null(label_scheme[[col_fill]])) {
		warning("Column \'", rlang::as_string(col_fill), "\' does not exist.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_fill]])) == 0) {
		warning("No samples were specified under column \'", rlang::as_string(col_fill), "\'.",
		        call. = FALSE)
	}

	if(is.null(label_scheme[[col_shape]])) {
		warning("Column \'", rlang::as_string(col_shape), "\' does not exist.")
	} else if(sum(!is.na(label_scheme[[col_shape]])) == 0) {
		warning("No samples were specified under column \'", rlang::as_string(col_shape), "\'.",
		        call. = FALSE)
	}

	if(is.null(label_scheme[[col_size]])) {
		warning("Column \'", rlang::as_string(col_size), "\' does not exist.")
	} else if(sum(!is.na(label_scheme[[col_size]])) == 0) {
		warning("No samples were specified under column \'", rlang::as_string(col_size), "\'.",
		        call. = FALSE)
	}

	if(is.null(label_scheme[[col_alpha]])) {
		warning("Column \'", rlang::as_string(col_alpha), "\' does not exist.")
	} else if(sum(!is.na(label_scheme[[col_alpha]])) == 0) {
		warning("No samples were specified under column \'", rlang::as_string(col_alpha), "\'.",
		        call. = FALSE)
	}

	if(is.null(label_scheme[[col_benchmark]])) {
		warning("Column \'", rlang::as_string(col_benchmark), "\' does not exist.")
	} else if(sum(!is.na(label_scheme[[col_benchmark]])) == 0) {
		warning("No samples were specified under column \'", rlang::as_string(col_benchmark), "\'.",
		        call. = FALSE)
	}

	anal_type <- rlang::as_string(rlang::enexpr(anal_type))

	id <- rlang::as_string(rlang::enexpr(id))
	if(length(id) != 1)
	  stop("\'id\' must be one of \'pep_seq\', \'pep_seq_mod\', \'prot_acc\' or \'gene\'")

	cat(paste0("id = \"", id, "\"", " by the current call\n"))
	id <- match_identifier(id)

	if(id %in% c("prot_acc", "gene")) {
	  cat(paste0("id = \"", id, "\"", " after parameter matching to normPrn()\n"))
		data_type <- "Protein"
	} else if(id %in% c("pep_seq", "pep_seq_mod")) {
	  cat(paste0("id = \"", id, "\"", " after parameter matching to normPep()\n"))
		data_type <- "Peptide"
	} else {
		stop("Unrecognized 'id';
		     needs to be in c(\"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\")", call. = TRUE)
	}

	if(is.null(filepath)) {
		filepath = file.path(dat_dir, data_type, anal_type)
		dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
	}

	if(is.null(filename)) {
		fn_prx <- paste(data_type, anal_type, sep = "_")
		fn_prx <- paste0(fn_prx, "_", ifelse(scale_log2r, "Z", "N"))
		fn_suffix <- "png"
	} else {
		fn_prx <- gsub("\\..*$", "", filename)
		fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
	}

	if(is.null(df)) {
		err_msg <- "File doesn't exist"

		if(id %in% c("pep_seq", "pep_seq_mod")) {
			fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
			fn_imp <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
			fn_raw <- file.path(dat_dir, "Peptide", "Peptide.txt")
		} else if(id %in% c("prot_acc", "gene")) {
			fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
			fn_imp <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
			fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")
		}

		if(anal_type %in% c("Histogram", "Corrplot", "MDS", "PCA", "EucDist", "MA")) {
		  if(file.exists(fn_raw)) src_path <- fn_raw else
		    stop(paste(fn_raw, "not found. \n Run normPSM(), normPep() and normPrn() first"),
		         call. = TRUE)
		} else if(anal_type %in% c("Heatmap", "Trend", "NMF")) {
		  if(impute_na) {
		    if(file.exists(fn_imp)) src_path <- fn_imp else
		      stop(paste(fn_imp, "not found. \nRun imputeNA() first."), call. = TRUE)
		  } else {
		    if(file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, "not found."),
		                                                         call. = TRUE)
		  }
		} else if(anal_type %in% c("Model", "ESGAGE", "GSVA")) {
			if(impute_na) {
		    if(file.exists(fn_imp)) src_path <- fn_imp else
		      stop(paste(fn_imp, "not found. \nRun imputeNA() first."), call. = TRUE)
		  } else {
		    if(file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, "not found."),
		                                                         call. = TRUE)
		  }
		} 

		df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
		                        comment.char = "#"), error = function(e) NA)

		if(!is.null(dim(df))) {
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
	  
	  if(!is.null(dim(df))) {
	    message(paste("File loaded:", gsub("\\\\", "/", fn_raw)))
	  } else {
	    stop(paste("Non-existed file or directory:", gsub("\\\\", "/", fn_raw)))
	  }
	}

	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	if(anal_type %in% c("Model", "ESGAGE", "GSVA")) {
		label_scheme_sub <- label_scheme %>% # to be subset by the "key" in "formulas"
			dplyr::filter(!grepl("^Empty\\.[0-9]+", .$Sample_ID), !Reference)
	} else {
		label_scheme_sub <- label_scheme %>%
			dplyr::select(Sample_ID, TMT_Set, !!col_select, !!col_group, !!col_order, !!col_color,
			              !!col_fill, !!col_shape, !!col_size, !!col_alpha, !!col_benchmark) %>%
			dplyr::filter(!is.na(!!col_select))
	}

	if (nrow(label_scheme_sub) == 0)
	  stop(paste0("No samples or conditions were defined for \"", anal_type, "\""))

	dfw <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r,
	              sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type)


	if(!any(names(df) == "kin_attr")) {
		cat("Columns of kinase annoation not found.\n")
		cat("Set 'annot_kinases = FALSE'.\n")
		annot_kinases <- FALSE
	}	else {
		cat("Columns of kinase annoation found.\n")
		cat("Set 'annot_kinases = TRUE'.\n")
		annot_kinases <- TRUE
	}

	if(annot_kinases) {
		dir.create(file.path(filepath, "Kinases"), recursive = TRUE, showWarnings = FALSE)

		dfw_kinase <- df %>%
		  dplyr::filter(kin_attr)

		if (nrow(dfw_kinase) > 0) {
		  dfw_kinase <- prepDM(dfw_kinase, !!id, scale_log2r,
		                       sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type)
		} else {
		  stop("No proteins were annotated being kinase. Please check the `Accession_Type`
		       and/or `Species` in `expt_smry.xlxs`.", call. = FALSE)
		}
		# prepDM(!!id, scale_log2r = scale_log2r, sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type)
	}

	force(scale_log2r)
	force(annot_kinases)
	force(anal_type)

	force(col_select)
	force(col_group)
	force(col_order)
	force(col_color)
	force(col_fill)
	force(col_shape)
	force(col_size)
	force(col_alpha)
	force(col_benchmark)

	if(anal_type == "MDS") {
		function(adjEucDist = adjEucDist, classical = classical, show_ids = show_ids,
		         annot_cols = NULL, ...) {

		  df_mds <- scoreMDS(df = dfw$log2R,
				label_scheme_sub = label_scheme_sub,
				scale_log2r = scale_log2r,
				adjEucDist,
				classical, ...
			)
		  
			df_mds$MDS %>% 
			  cmbn_meta(label_scheme_sub) %>% 
			  plotMDS(col_color = !!col_color, 
			          col_fill = !!col_fill, col_shape = !!col_shape, 
			          col_size = !!col_size, col_alpha = !!col_alpha, 
			          label_scheme_sub = label_scheme_sub, 
			          filepath = filepath, 
			          filename = paste0(fn_prx, ".", fn_suffix), 
			          show_ids = show_ids, ...)

			if(annot_kinases) {
				df_kinase_mds <- scoreMDS(df = dfw_kinase$log2R,
					label_scheme_sub = label_scheme_sub,
				  scale_log2r = scale_log2r, adjEucDist, classical, ...
				)

				df_kinase_mds$MDS %>% 
				  cmbn_meta(label_scheme_sub) %>% 
				  plotMDS(col_color = !!col_color, 
				          col_fill = !!col_fill, col_shape = !!col_shape, 
				          col_size = !!col_size, col_alpha = !!col_alpha, 
				          label_scheme_sub = label_scheme_sub, 
				          filepath = file.path(filepath, "Kinases"), 
				          filename = paste0(fn_prx, "_Kinases.", fn_suffix), 
				          show_ids = show_ids, ...)
			}
		}
	} else if(anal_type == "PCA") {
		function(show_ids = show_ids, annot_cols = NULL, ...) {
			df_mds <- scoreMDS(df = dfw$log2R,
				label_scheme_sub = label_scheme_sub,
			  scale_log2r = scale_log2r,
				adjEucDist = FALSE,
				classical = TRUE,
				...
			)
			
			df_mds$PCA %>% 
			  cmbn_meta(label_scheme_sub) %>% 
			  plotPCA(col_color = !!col_color, 
			          col_fill = !!col_fill, 
			          col_shape = !!col_shape, 
			          col_size = !!col_size, 
			          col_alpha = !!col_alpha, 
			          label_scheme_sub = label_scheme_sub, 
			          prop_var = df_mds$prop_var, 
			          pr_bi = df_mds$pr_bi, 
			          prop_var_bi = df_mds$prop_var_bi, 
			          filepath = filepath, 
			          filename = paste0(fn_prx, ".", fn_suffix), 
			          show_ids = show_ids, ...)

			if(annot_kinases) {
				df_kinase_mds <- scoreMDS(df = dfw_kinase$log2R,
					label_scheme_sub = label_scheme_sub,
				  scale_log2r = scale_log2r,
					adjEucDist = FALSE,
					classical = TRUE, ...
				)

				df_kinase_mds$PCA %>% 
				  cmbn_meta(label_scheme_sub) %>% 
				  plotPCA(col_color = !!col_color, 
				          col_fill = !!col_fill, 
				          col_shape = !!col_shape, 
				          col_size = !!col_size, 
				          col_alpha = !!col_alpha, 
				          label_scheme_sub = label_scheme_sub, 
				          prop_var = df_kinase_mds$prop_var, 
				          pr_bi = df_kinase_mds$pr_bi, 
				          prop_var_bi = df_kinase_mds$prop_var_bi, 
				          filepath = file.path(filepath, "Kinases"), 
				          filename = paste0(fn_prx, "_Kinases.", fn_suffix), 
				          show_ids = show_ids, ...)
			}
		}
	} else if(anal_type == "EucDist") {
		function(adjEucDist = adjEucDist, annot_cols = NULL, annot_colnames = NULL, ...) {
			df_mds <- scoreMDS(df = dfw$log2R,
				label_scheme_sub = label_scheme_sub,
			  scale_log2r = scale_log2r,
				adjEucDist = adjEucDist,
				classical = TRUE, ...
			)

			df_mds$D %>% plotEucDist(annot_cols, 
			                         annot_colnames, 
			                         filepath = filepath, 
			                         filename = paste0(fn_prx, ".", fn_suffix), ...
			)

			if(annot_kinases) {
				df_kinase_mds <- scoreMDS(df = dfw_kinase$log2R,
					label_scheme_sub = label_scheme_sub,
				  scale_log2r = scale_log2r,
					adjEucDist = adjEucDist,
					classical = TRUE, ...
				)

				df_kinase_mds$D %>% plotEucDist(annot_cols, 
				                                annot_colnames, 
				                                filepath = file.path(filepath, "Kinases"), 
				                                filename = paste0(fn_prx, "_Kinases.", fn_suffix), ...)
			}
		}
	} else if(anal_type == "MA") {
		function(...) {

			plotMA(df = dfw, col_group = !!col_group, label_scheme_sub = label_scheme_sub,
			       filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix), ...)

			if(annot_kinases) {
				plotMA(df = dfw_kinase, col_group = !!col_group, label_scheme_sub = label_scheme_sub,
				       filepath = file.path(filepath, "Kinases"),
				       filename = paste0(fn_prx, "_Kinases.", fn_suffix), ...)
			}
		}
	} else if(anal_type == "Heatmap") {
		function(complete_cases = complete_cases, xmin = xmin, xmax = xmax, x_margin = x_margin,
		         annot_cols = annot_cols, annot_colnames = annot_colnames, ...) {
			plotHM(df = df, id = !!id, scale_log2r = scale_log2r, col_benchmark = !!col_benchmark,
			       label_scheme_sub = label_scheme_sub,
						filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix),
						complete_cases = complete_cases, xmin = xmin, xmax = xmax, x_margin = x_margin,
						annot_cols, annot_colnames, ...)

			if(annot_kinases) {
				plotKinHM(df = df[df$kin_attr, ], id = !!id, scale_log2r = scale_log2r,
				          col_benchmark = !!col_benchmark, label_scheme_sub = label_scheme_sub,
								filepath = file.path(filepath, "Kinases"),
								filename = paste0(fn_prx, "_Kinases.", fn_suffix),
			 					complete_cases = complete_cases, xmin = xmin, xmax = xmax, x_margin = x_margin,
								annot_cols, annot_colnames, ...)
			}
		}
	} else if(anal_type == "Histogram") {
		function(new_fit = new_fit, show_curves = show_curves, show_vline = show_vline, ...) {

			if (scale_log2r) {
				fn <- file.path(filepath, "MGKernel_params_Z.txt")
			} else {
				fn <- file.path(filepath, "MGKernel_params_N.txt")
			}

			if (file.exists(fn)) {
				params <- read.csv(fn, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")
				params <- within(params, {mean = mean - x})
			} else {
				params <- NULL
				show_curves <- FALSE
			}

			plotHisto(df = df, label_scheme_sub = label_scheme_sub, params = params,
			          scale_log2r = scale_log2r, show_curves = show_curves, show_vline = show_vline,
			          filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix), ...)

			if(annot_kinases) {
				if(is.null(dim(df[df$kin_attr, ]))) stop("Kinase annoation not found.", call. = TRUE)

				plotHisto(df = df[df$kin_attr, ], label_scheme_sub = label_scheme_sub, params = params,
				          scale_log2r = scale_log2r, show_curves = FALSE, show_vline = show_vline,
				          filepath = file.path(filepath, "Kinases"),
				          filename = paste0(fn_prx, "_Kinases.", fn_suffix), ...)
			}
		}
	} else if (anal_type == "Corrplot") {
		function(use_log10 = use_log10, min_int = min_int, max_int = max_int, min_log2r = min_log2r,
		         max_log2r = max_log2r, ...) {
			if(ncol(dfw$log2R) > 44) stop("Maximal number of samples for correlation plots is 44.", call. = TRUE)

			plotCorr(df = dfw, col_select = !!col_select, col_order = !!col_order,
			         label_scheme_sub = label_scheme_sub, use_log10 = use_log10,
			         scale_log2r = scale_log2r, min_int = min_int, max_int = max_int,
			         min_log2r = min_log2r, max_log2r = max_log2r,
							 filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix), ...)

			if(annot_kinases) {
				plotCorr(df = dfw_kinase, col_select = !!col_select, col_order = !!col_order,
				         label_scheme_sub = label_scheme_sub, use_log10 = use_log10,
				         scale_log2r = scale_log2r, min_int = min_int, max_int = max_int,
				         min_log2r = min_log2r, max_log2r = max_log2r,
								 filepath = file.path(filepath, "Kinases"),
								 filename = paste0(fn_prx, "_Kinases.", fn_suffix), ...)
			}
		}
	} else if (anal_type == "Trend") {
		function(n_clust = n_clust, complete_cases = complete_cases, task = !!task, ...) {
		  switch(rlang::as_string(rlang::enexpr(task)), 
		         anal = trendTest(df = dfw$log2R, id = !!id, 
		                          col_group = !!col_group, col_order = !!col_order,
		                          label_scheme_sub = label_scheme_sub, n_clust = n_clust,
		                          complete_cases = complete_cases, scale_log2r = scale_log2r,
		                          filepath = filepath, 
		                          filename = paste0(fn_prx, "_n", n_clust, ".csv"), ...), 
		         plot = plotTrend(id = !!id, 
		                          col_group = !!col_group, 
		                          col_order = !!col_order, 
		                          label_scheme_sub = label_scheme_sub, 
		                          n_clust = n_clust, 
		                          filepath = filepath, 
		                          in_nm = paste0(fn_prx, "_n", n_clust, ".csv"), 
		                          out_nm = paste0(fn_prx, "_n", n_clust, ".png"), ...)
		  )
		}
	} else if (anal_type == "NMF") {
		function(r = r, nrun = nrun, complete_cases = complete_cases, task = !!task, ...) {
		  switch(rlang::as_string(rlang::enexpr(task)), 
		         anal = nmfTest(df = dfw$log2R, id = !!id, r = r, nrun = nrun, 
		                        col_group = !!col_group, label_scheme_sub = label_scheme_sub,
		                        filepath = filepath, filename = paste0(fn_prx, "_r", r, ".csv"),
		                        complete_cases = complete_cases, ...), 
		         plotcon = plotNMFCon(id = !!id, 
		                               r = r, 
		                               label_scheme_sub = label_scheme_sub, 
		                               filepath = filepath, 
		                               in_nm = paste0(fn_prx, "_r", r, ".rda"), 
		                               out_nm = paste0(fn_prx, "_r", r, "_consensus.png"), ...), 
		         plotcoef = plotNMFCoef(id = !!id, 
		                              r = r, 
		                              label_scheme_sub = label_scheme_sub, 
		                              filepath = filepath, 
		                              in_nm = paste0(fn_prx, "_r", r, ".rda"), 
		                              out_nm = paste0(fn_prx, "_r", r, "_coef.png"), ...), 
		         plotmeta = plotNMFmeta(df = dfw$log2R, id = !!id, r = r, 
		                                label_scheme_sub = label_scheme_sub, 
		                                filepath = filepath, 
		                                in_nm = paste0(fn_prx, "_r", r, ".rda"), 
		                                out_nm = paste0(fn_prx, "_r", r, "_metagene.png"), ...)
		  )

		}
	} else if(anal_type == "ESGAGE") {
		function(complete_cases = FALSE, method = "limma", gset_nm = "go_sets", var_cutoff = 1E-3,
		         pval_cutoff = .05, logFC_cutoff = log2(1.1), ...) {
			df_op <- gageTest(df = dfw$log2R, id = !!id, label_scheme_sub = label_scheme_sub,
			                  filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix),
			                  complete_cases, method = method, gset_nm = gset_nm, var_cutoff = var_cutoff,
			                  pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, ...)
		}
	} else if(anal_type == "GSVA") {
		function(complete_cases = FALSE, method = "limma", gset_nm = "go_sets", var_cutoff = .5,
		         pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), mx.diff = TRUE, ...) {

			# "id" only for tibbling rownames
			df_op <- gsvaTest(df = dfw$log2R, id = !!id, label_scheme_sub = label_scheme_sub,
				filepath, filename = paste0(fn_prx, ".", fn_suffix),
				complete_cases = complete_cases, method = method,
				gset_nm = gset_nm, var_cutoff = var_cutoff,
				pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, mx.diff = mx.diff, ...)
		}

	} else if (anal_type == "Model") {
		function(complete_cases = FALSE, method = "limma", var_cutoff = 1E-3, pval_cutoff = 1,
		         logFC_cutoff = log2(1), ...) {
		  method <- rlang::enexpr(method)
		  
		  df_op <- sigTest(df = dfw$log2R, id = !!id, label_scheme_sub = label_scheme_sub,
			                 filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix),
			                 complete_cases = complete_cases, method = !!method, var_cutoff = var_cutoff,
			                 pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, ...) %>%
								tibble::rownames_to_column(id) %>%
								dplyr::right_join(df, ., by = id)

			write.table(df_op, file.path(filepath, paste0(data_type, "_pVals.txt")), sep = "\t",
			            col.names = TRUE, row.names = FALSE)
		}
	}
}


