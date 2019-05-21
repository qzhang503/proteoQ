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
		} else if(anal_type %in% c("Trend", "NMF")) {
		  if(file.exists(fn_imp)) src_path <- fn_imp else
		    stop(paste(fn_imp, "not found. \nRun imputeNA() first."), call. = TRUE)
		} else if(anal_type %in% c("Heatmap")) {
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
			          filename = paste0(fn_prx, ".png"), 
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
			          filename = paste0(fn_prx, ".png"), 
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
		function(adjEucDist = adjEucDist, annot_cols = NULL, annot_colnames, ...) {
			df_mds <- scoreMDS(df = dfw$log2R,
				label_scheme_sub = label_scheme_sub,
			  scale_log2r = scale_log2r,
				adjEucDist = adjEucDist,
				classical = TRUE, ...
			)

			df_mds$D %>% plotEucDist(annot_cols, 
			                         annot_colnames, 
			                         filepath = filepath, 
			                         filename = paste0(fn_prx, ".png"), ...
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
		         annot_cols = annot_cols, ...) {

			plotHM(df = df, id = !!id, scale_log2r = scale_log2r, col_benchmark = !!col_benchmark,
			       label_scheme_sub = label_scheme_sub,
						filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix),
						complete_cases = complete_cases, xmin = xmin, xmax = xmax, x_margin = x_margin,
						annot_cols, ...)

			if(annot_kinases) {
				plotKinHM(df = df[df$kin_attr, ], id = !!id, scale_log2r = scale_log2r,
				          col_benchmark = !!col_benchmark, label_scheme_sub = label_scheme_sub,
								filepath = file.path(filepath, "Kinases"),
								filename = paste0(fn_prx, "_Kinases.", fn_suffix),
			 					complete_cases = complete_cases, xmin = xmin, xmax = xmax, x_margin = x_margin,
								annot_cols, ...)
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
		function(n_clust = n_clust, complete_cases = complete_cases, ...) {

			plotTrend(df = dfw$log2R, id = !!id, col_group = !!col_group, col_order = !!col_order,
			          label_scheme_sub = label_scheme_sub, n_clust = n_clust,
			          complete_cases = complete_cases, scale_log2r = scale_log2r,
			          filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix), ...)

			if(annot_kinases) {
				plotTrend(df = dfw_kinase$log2R, col_group = !!col_group, col_order = !!col_order, id = !!id,
				          label_scheme_sub = label_scheme_sub, n_clust = n_clust,
				          complete_cases = complete_cases, scale_log2r = scale_log2r,
				          filepath = file.path(filepath, "Kinases"),
				          filename = paste0(fn_prx, "_Kinases", ".png"), ...)
			}
		}

	} else if (anal_type == "NMF") {
		function(r = r, nrun = nrun,
						complete_cases = complete_cases,
						xmin = -1, xmax = 1, x_margin = .1, annot_cols = annot_cols,
						width_consensus = width_consensus, height_consensus = height_consensus,
						width_coefmap = width_coefmap, height_coefmap = height_coefmap, ...) {

			plotNMF(df = dfw$log2R, id = !!id, col_group = !!col_group,
			        label_scheme_sub = label_scheme_sub,
							filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix),
			        r = r, nrun = nrun, complete_cases = complete_cases,
			        xmin = xmin, xmax = xmax, x_margin = x_margin, annot_cols = annot_cols,
			        width_consensus = width_consensus, height_consensus = height_consensus,
			        width_coefmap = width_coefmap, height_coefmap = height_coefmap, ...)
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

			df_op <- sigTest(df = dfw$log2R, id = !!id, label_scheme_sub = label_scheme_sub,
			                 filepath = filepath, filename = paste0(fn_prx, ".", fn_suffix),
			                 complete_cases = complete_cases, method = method, var_cutoff = var_cutoff,
			                 pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, ...) %>%
								tibble::rownames_to_column(id) %>%
								dplyr::right_join(df, ., by = id)

			write.table(df_op, file.path(filepath, paste0(data_type, "_pVals.txt")), sep = "\t",
			            col.names = TRUE, row.names = FALSE)
		}

	}

}


#' Prepares data for analysis
#'
#' \code{prepDM} prepares a minimal data frame for subsequent analysis.
#'
#' @param df A data frame containing only numeric values.
#' @param id The name of unqiue identifiers.
#' @param scale_log2r Logical; if TRUE, rescales \code{log2-ratios} to the same
#'   scale.
#' @param sub_grp Numeric.  A list of sample IDs that will be used in subsequent
#'   analysis.
#' @return A data frame tailored for subsequent analysis.
#'
#' @examples
#' tempData <- prepDM(df, entrez, scale_log2r, sub_grp = label_scheme_sub$Sample_ID)
#'
#' \dontrun{
#' }
#' @import dplyr
#' @importFrom magrittr %>%
prepDM <- function(df, id, scale_log2r, sub_grp, type = "ratio", anal_type) {
	stopifnot(nrow(df) > 0)

  load(file = file.path(dat_dir, "label_scheme.Rdata"))

  id <- rlang::as_string(rlang::enexpr(id))
	if(anal_type %in% c("ESGAGE", "GSVA")) id <- "entrez"

	NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")

	# data filtration dominated by log2R, not Intensity
	df <- df %>%
			dplyr::filter(!duplicated(!!rlang::sym(id)),
			              !is.na(!!rlang::sym(id)),
			              rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
			reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id)

	Levels <- sub_grp %>%
		as.character(.) %>%
		.[!grepl("^Empty\\.[0-9]+", .)]

	dfR <- df %>%
			dplyr::select(grep(NorZ_ratios, names(.))) %>%
			`colnames<-`(label_scheme$Sample_ID) %>%
			dplyr::select(which(names(.) %in% sub_grp)) %>%
			dplyr::select(which(not_all_zero(.))) %>% # reference will drop with single reference
			dplyr::select(Levels[Levels %in% names(.)]) # ensure the same order

	# dominated by log2R, no need to filter all-NA intensity rows
	dfI <- df %>%
			dplyr::select(grep("^N_I[0-9]{3}", names(.))) %>%
			`colnames<-`(label_scheme$Sample_ID) %>%
			dplyr::select(which(names(.) %in% sub_grp)) %>%
			dplyr::select(which(not_all_zero(.))) %>%
			dplyr::select(which(names(.) %in% names(dfR))) %>%
			dplyr::select(Levels[Levels %in% names(.)])

	tempI <- dfI
	rownames(tempI) <- paste(rownames(tempI), "Intensity", sep = "@")

	tempR <- dfR
	rownames(tempR) <- paste(rownames(tempR), "log2R", sep = "@")

	return(list(log2R = dfR, Intensity = dfI, IR = rbind(tempR, tempI)))
}


#' A wrapper around pheatmap
#'
#' @import dplyr rlang pheatmap
#' @importFrom magrittr %>%
my_pheatmap <- function(mat, filename, annotation_col, color, annotation_colors, breaks, ...) {
	mat <- rlang::enexpr(mat)
	filename <- rlang::enexpr(filename)
	annotation_col <- rlang::enexpr(annotation_col)
	color <- rlang::enexpr(color)
	annotation_colors <- rlang::enexpr(annotation_colors)
	breaks <- rlang::enexpr(breaks)

	dots <- rlang::enexprs(...)

	dots$mat <- NULL
	dots$filename <- NULL
	dots$annotation_col <- NULL
	dots$color <- NULL
	dots$annotation_colors <- NULL
	dots$breaks <- NULL

	ph_call <- rlang::expr(
	  pheatmap(mat = !!mat, filename = !!filename, annotation_col = !!annotation_col, color = !!color,
		annotation_colors = !!annotation_colors, breaks = !!breaks, !!!dots))

	rlang::expr_print(ph_call)
	rlang::eval_bare(ph_call, env = caller_env())
}


#' Makes heat maps
#'
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
plotHM <- function(df, id, scale_log2r, col_benchmark, label_scheme_sub, filepath, filename,
                   complete_cases, xmin = -1, xmax = 1, x_margin = .1, annot_cols = NULL, ...) {

	id <- rlang::as_string(rlang::enexpr(id))
	col_benchmark <- rlang::as_string(rlang::enexpr(col_benchmark))

	dots <- rlang::enexprs(...)

	if (is.null(dots$cluster_rows)) {
	  cluster_rows <- TRUE
	} else {
	  cluster_rows <- dots$cluster_rows
	}

	if(is.null(dots$clustering_distance_rows)) {
	  clustering_distance_rows <- "euclidean"
	} else {
	  clustering_distance_rows <- dots$clustering_distance_rows
	}

	# parameter(s) in "dots" being disabled for pheatmap()
	dots$annotation_col <- NULL
	nm_idx <- names(dots) %in% c("df", "id", "scale_log2r", "annot_kinases",
	                             "complete_cases", "label_scheme_sub", "filepath", "filename")
	dots[nm_idx] <- NULL

	fn_prx <- gsub("\\..*$", "", filename)
	fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)

	dir.create(file.path(filepath, "Subtrees"), recursive = TRUE, showWarnings = FALSE)

	x_label <- expression("Ratio ("*log[2]*")")
	NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
	NorZ_ratios_to_ctrl <- paste("toCtrl", NorZ_ratios, sep = "_")

	n_color <- 500
	if (is.null(dots$breaks)) {
		color_breaks <- c(seq(from = xmin, -x_margin, length = n_color/2)[1:(n_color/2-1)],
											seq(-x_margin, x_margin, length = 3),
											seq(x_margin, xmax, length = n_color/2)[2:(n_color/2)])
	} else if (is.na(dots$breaks)) {
		color_breaks <- NA
	} else {
		color_breaks <- eval(dots$breaks, env = caller_env())
	}

	if (is.null(dots$color)) {
		mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
	} else if (is.na(dots$color)) {
		mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
	} else {
		mypalette <- eval(dots$color, env = caller_env())
	}

	# mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(n_color)
	# mypalette <- colorRampPalette(c("darkorange3", "white", "darkblue"))(n_color)
	# image(matrix(1:n_color, nrow = n_color, ncol = 10),
	#   col = mypalette, xaxt = "n", yaxt = "n", useRaster = TRUE)

	acc_type <- load(file = file.path(dat_dir, "label_scheme.Rdata")) %>% find_acctype()

	if(acc_type == "refseq_acc") {
		df <- df %>% mutate(Species = gsub(".*\\[(.*)\\].*", "\\1", .$prot_desc))
	} else if (acc_type == "uniprot_id") {
		df <- df %>% mutate(Species = gsub(".*\\s+OS\\=(.*)\\s+GN.*", "\\1", .$prot_desc))
	}

	sample_ids <- label_scheme_sub$Sample_ID

	df <- df %>%
		dplyr::mutate_at(vars(grep("^pVal|^adjP", names(.))), as.numeric) %>%
		dplyr::mutate(Mean_log10Int = log10(rowMeans(.[, grepl("^I[0-9]{3}", names(.))],
		                                             na.rm = TRUE))) %>%
		dplyr::mutate_at(vars(grep("log2_R[0-9]{3}", names(.))), ~setHMlims(., xmin, xmax)) %>%
		dplyr::filter(!duplicated(!!rlang::sym(id)),
		              !is.na(!!rlang::sym(id)),
		              rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
		reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id)

	dfR <- df %>%
			dplyr::select(grep(NorZ_ratios, names(.))) %>%
			`colnames<-`(label_scheme$Sample_ID) %>%
			dplyr::select(which(names(.) %in% sample_ids)) %>%
			dplyr::select(which(not_all_zero(.))) %>%
			dplyr::select(as.character(sample_ids)) # ensure the same order

	df <- df %>%
		dplyr::select(-grep("log2_R[0-9]{3}", names(.))) %>%
		dplyr::bind_cols(., dfR) %>%
		`rownames<-`(.[[id]])

	if(is.null(annot_cols)) {
		annotation_col <- NA
	} else {
		annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)
	}

	if (is.null(dots$annotation_colors)) {
		annotation_colors <- setHMColor(annotation_col)
	} else if (is.na(dots$annotation_colors)) {
		annotation_colors <- NA
	} else {
		annotation_colors <- eval(dots$annotation_colors, env = caller_env())
	}

	if (complete_cases) {
		df_hm <- df %>%
			dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
	} else {
		df_hm <- df
	}

	df_hm <- df_hm %>%
			`rownames<-`(.[[id]])	%>%
			dplyr::select(which(names(.) %in% sample_ids))

	if(cluster_rows) {
		d <- dist(df_hm, method = clustering_distance_rows)
		d[is.na(d)] <- .5 * max(d, na.rm = TRUE)
		h <- hclust(d)
		dots$cluster_rows <- h

		rm(d, h)
	} else {
		dots$cluster_rows <- FALSE
	}

	p <- my_pheatmap(
		mat = df_hm,
		filename = file.path(filepath, filename),
		annotation_col = annotation_col,
		color = mypalette,
		annotation_colors = annotation_colors,
		breaks = color_breaks,
		!!!dots
	)

	# when cutree_rows is not NA
	cutree_rows <- eval(dots$cutree_rows, env = caller_env())

	if(!is.null(cutree_rows) & cluster_rows) {
		if(is.numeric(cutree_rows)) {
			Cluster <- data.frame(Cluster = cutree(p$tree_row, k = cutree_rows)) %>%
				dplyr::mutate(!!id := rownames(.)) %>%
				dplyr::left_join(df, by = id)

			write.csv(Cluster,
			          file.path(filepath, "Subtrees", paste0(fn_prx, " n-", cutree_rows, "_subtrees.csv")),
			          row.names = FALSE)

			Cluster <- Cluster %>%
				tibble::column_to_rownames(var = id)

			for(cluster_id in unique(Cluster$Cluster)) {
				df_sub <- Cluster[Cluster$Cluster == cluster_id, names(Cluster) %in% sample_ids]

				if(complete_cases) {
					df_sub <- df_sub %>%
						tibble::rownames_to_column(id) %>%
						dplyr::filter(complete.cases(.[, names(.) %in% sample_ids])) %>%
						tibble::column_to_rownames(id)
				}

				if(cluster_rows) {
					d_sub <- dist(df_sub, method = clustering_distance_rows)
					d_sub[is.na(d_sub)] <- .5 * max(d_sub, na.rm = TRUE)

					if (nrow(df_sub) < 2) {
					  h_sub <- FALSE
					} else {
					  h_sub <- hclust(d_sub)
					}
				}

				pheatmap(
					mat = df_sub,
					main = paste("Cluster", cluster_id),
					cluster_rows = h_sub,
					show_rownames = TRUE,
					show_colnames = TRUE,
					annotation_col = annotation_col,
					color = mypalette,
					breaks = color_breaks,
					cellwidth = 14,
					fontsize_row = ceiling(200/nrow(Cluster %>% dplyr::filter(Cluster == cluster_id))),
					annotation_colors = annotation_colors,
					filename = file.path(filepath, "Subtrees",
					                     paste0("Subtree_", cutree_rows, "-", cluster_id, ".png"))
				)

				rm(d_sub, h_sub)
			}
		}
	}

	# to the benchmark group
	if(sum(!is.na(label_scheme_sub[[col_benchmark]])) > 0) {
		dfc <- ratio_toCtrl(df, !!rlang::sym(id), label_scheme_sub, nm_ctrl = !!col_benchmark)

		if (complete_cases) {
			dfc_hm <- dfc %>%
				dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
		} else {
			dfc_hm <- dfc %>%
				dplyr::filter(rowSums(!is.na(.[, names(.) %in% sample_ids])) > 0)
		}

		dfc_hm <- dfc_hm %>%
				`rownames<-`(.[[id]])	%>%
				dplyr::select(which(names(.) %in% sample_ids))

		if(cluster_rows) {
			d <- dist(dfc_hm, method = clustering_distance_rows)
			d[is.na(d)] <- .5  * max(d, na.rm = TRUE)
			h <- hclust(d)
			dots$cluster_rows <- h

			rm(d, h)
		} else {
			dots$cluster_rows <- FALSE
		}

		p <- my_pheatmap(
			mat = dfc_hm,
			filename = file.path(filepath, paste0(fn_prx, "_toCtrl.", fn_suffix)),
			annotation_col = annotation_col,
			color = mypalette,
			annotation_colors = annotation_colors,
			breaks = color_breaks,
			!!!dots
		)
	}
}


#'Plots heat maps
#'
#'\code{proteoHM} produces the heat map visualization of \code{log2-ratios} for
#'proteins or peptides data.
#'
#'Data columns with complete missing values will be removed prior to
#'hierarchical column clustering. Data rows without non-missing pairs will
#'result in NA distances in inter-row dissimilarities
#'(\code{\link[stats]{dist}}). At \code{complet_cases = TRUE}, the data subset
#'that are complete without missing values will be used. At \code{impute_na =
#'TRUE}, all data rows will be used with NA imputation (\code{\link{imputeNA}}).
#'At the default of \code{complet_cases = FALSE} and \code{impute_na = FALSE},
#'NA distances will be arbitrarily replaced with the mean value in the
#'row-disance matrix for hierarchical row clustering.
#'
#'To avoid memory failure, row aggregation using the \code{kmeans_k} option
#'(\code{\link[pheatmap]{pheatmap}}) may be considered for large data sets.
#'
#'The function matches the current \code{id} to those in the latest \code{calls}
#'to \code{\link{normPep}} or \code{\link{normPrn}}.  For example, if
#'\code{pep_seq} was used in \code{\link{normPep}}, the current \code{id =
#'pep_seq_mod} will be matched to \code{id = pep_seq}.
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at \code{id = pep_seq} or \code{pep_seq_mod}, and protein data at
#'  \code{id = prot_acc} or \code{gene}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param  col_benchmark Character string to a column key in
#'  \code{expt_smry.xlsx}. The default key is \code{Benchmark}. Samples
#'  corresponding to non-empty entries under \code{col_benchmark} will be used
#'  as benchmarks in the indicated analysis.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param impute_na Logical; if TRUE, imputes missing values.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used for visualization.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the name of the current function \code{call} and the value of
#'  \code{id}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... More parameters for plotting: \cr \code{xmin}, the minimum x; \cr
#'  \code{xmax}, the maximum x; \cr \code{x_margin}, the margin in heat
#'  scales;\cr \code{annot_cols}, the column keys in \code{expt_smry.xlsx} for
#'  use in the color coding of samples;\cr \code{width}, the width of plot; \cr
#'  \code{height}, the height of plot; \cr additional arguments inherited from
#'  \code{\link[pheatmap]{pheatmap}}.
#'@return Heat map images.
#'
#' @examples
#'proteoHM(
#'  id = prot_acc,
#'  scale_log2r = TRUE,
#'  xmin = -1,
#'  xmax = 1,
#'  x_margin = 0.1,
#'  annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'  cluster_rows = TRUE,
#'  cutree_rows = 6,
#'  show_rownames = FALSE,
#'  show_colnames = TRUE,
#'  fontsize_row = 3,
#'  cellwidth = 14,
#'  width = 10,
#'  height = 12
#')
#'
#'proteoHM(
#'  id = prot_acc,
#'  scale_log2r = TRUE,
#'  xmin = -1,
#'  xmax = 1,
#'  x_margin = 0.1,
#'  annot_cols = c("Group"),
#'  cluster_rows = TRUE,
#'  clustering_distance_rows  = "maximum",
#'  cutree_rows = 6,
#'  show_rownames = FALSE,
#'  show_colnames = TRUE,
#'  fontsize_row = 3,
#'  cellwidth = 14,
#'  width = 10,
#'  height = 12
#')
#'
#'proteoHM(
#'  id = prot_acc,
#'  scale_log2r = TRUE,
#'  xmin = -1,
#'  xmax = 1,
#'  x_margin = 0.1,
#'  annot_cols = c("Group"),
#'  cluster_rows = FALSE, # no row clustering
#'  clustering_distance_rows  = "maximum",
#'  cutree_rows = 6, # will overrule tree cutting at 'cluster_rows = FALSE'
#'  show_rownames = FALSE,
#'  show_colnames = TRUE,
#'  fontsize_row = 3,
#'  cellwidth = 14,
#'  width = 10,
#'  height = 12
#')
#'
#' \dontrun{
#' }
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoHM <- function (id = gene, col_select = NULL, col_benchmark = NULL,
                      scale_log2r = TRUE,impute_na = FALSE, complete_cases = FALSE,
											df = NULL, filepath = NULL, filename = NULL,
											xmin = -1, xmax = 1, x_margin = 0.1, annot_cols = NULL, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

  id <- rlang::enexpr(id)
	col_select <- rlang::enexpr(col_select)
	col_benchmark <- rlang::enexpr(col_benchmark)
	
	reload_expts()

	info_anal(id = !!id, col_select = !!col_select, col_benchmark = !!col_benchmark,
	          scale_log2r = scale_log2r, impute_na = impute_na, df = df, filepath = filepath,
	          filename = filename, anal_type = "Heatmap")(complete_cases = complete_cases,
	                                                      xmin = xmin, xmax = xmax,
	                                                      x_margin = x_margin,
	                                                      annot_cols = annot_cols, ...)
}


#'Visualizes peptide heat maps
#'@seealso \code{\link{proteoHM}} for parameters
#'@export
pepHM <- function (...) {
	proteoHM(id = pep_seq, ...)
}


#'Visualizes protein heat maps
#'@seealso \code{\link{proteoHM}} for parameters
#'@export
prnHM <- function (...) {
	proteoHM(id = gene, ...)
}


#' Plots kinase heat maps
#'
#' Specialized for kinase heat maps
#'
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
#' @export
plotKinHM <- function(id, scale_log2r, col_benchmark, label_scheme_sub, df, filepath, filename,
                      complete_cases, xmin = -1, xmax = 1, x_margin = .1, annot_cols = NULL, ...) {

  id <- rlang::as_string(rlang::enexpr(id))
	col_benchmark <- rlang::as_string(rlang::enexpr(col_benchmark))

	dots <- rlang::enexprs(...)

	# parameter(s) in "dots" disabled for pheatmap()
	dots$annotation_col <- NULL
	nm_idx <- names(dots) %in% c("df", "id", "scale_log2r", "annot_kinases", "complete_cases",
	                             "label_scheme_sub", "filepath", "filename")
	dots[nm_idx] <- NULL

	fn_prx <- gsub("\\..*$", "", filename)
	fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)

	dir.create(file.path(filepath, "Subtrees"), recursive = TRUE, showWarnings = FALSE)

	x_label <- expression("Ratio ("*log[2]*")")
	NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
	NorZ_ratios_to_ctrl <- paste("toCtrl", NorZ_ratios, sep = "_")

	n_color <- 500
	if (is.null(dots$breaks)) {
		color_breaks <- c(seq(from = xmin, -x_margin, length = n_color/2)[1:(n_color/2-1)],
											seq(-x_margin, x_margin, length = 3),
											seq(x_margin, xmax, length = n_color/2)[2:(n_color/2)])
	} else if (is.na(dots$breaks)) {
		color_breaks <- NA
	} else {
		color_breaks <- eval(dots$breaks, env = caller_env())
	}

	if (is.null(dots$color)) {
		mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
	} else if (is.na(dots$color)) {
		mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
	} else {
		mypalette <- eval(dots$color, env = caller_env())
	}

	acc_type <- load(file = file.path(dat_dir, "label_scheme.Rdata")) %>% find_acctype()

	if(acc_type == "refseq_acc") {
		df <- df %>% mutate(Species = gsub(".*\\[(.*)\\].*", "\\1", .$prot_desc))
	} else if (acc_type == "uniprot_id") {
		df <- df %>% mutate(Species = gsub(".*\\s+OS\\=(.*)\\s+GN.*", "\\1", .$prot_desc))
	}

	sample_ids <- label_scheme_sub$Sample_ID

	kin_levels <- c("TK", "TKL", "STE", "CK1", "AGC", "CAMK", "CMGC", "Atypical", "Other", "RGC",
	                "Unclassified", "FAM20", "Lipid")

	df <- df %>%
		dplyr::mutate_at(vars(grep("^pVal|^adjP", names(.))), as.numeric) %>%
		dplyr::mutate(Mean_log10Int = log10(rowMeans(.[, grepl("^I[0-9]{3}", names(.))],
		                                             na.rm = TRUE))) %>%
		dplyr::mutate_at(vars(grep("log2_R[0-9]{3}", names(.))), ~setHMlims(., xmin, xmax)) %>%
		dplyr::filter(!duplicated(!!rlang::sym(id)),
		              !is.na(!!rlang::sym(id)),
		              rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
		reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id) %>%
		dplyr::mutate(kin_class = factor(kin_class, levels = kin_levels)) %>%
		dplyr::arrange(kin_class, !!rlang::sym(id))

	dfR <- df %>%
			dplyr::select(grep(NorZ_ratios, names(.))) %>%
			`colnames<-`(label_scheme$Sample_ID) %>%
			dplyr::select(which(names(.) %in% sample_ids)) %>%
			dplyr::select(which(not_all_zero(.))) %>%
			dplyr::select(as.character(sample_ids)) # ensure the same order

	df <- df %>%
		dplyr::select(-grep("log2_R[0-9]{3}", names(.))) %>%
		dplyr::bind_cols(., dfR) %>%
		`rownames<-`(.[[id]])

	annotation_row <- df %>%
			dplyr::select(id, kin_class) %>%
			dplyr::mutate(kin_class = factor(kin_class, levels =
			                                   c("TK", "TKL", "STE", "CK1", "AGC", "CAMK", "CMGC",
			                                     "Atypical", "Other", "RGC", "Unclassified",
			                                     "FAM20", "Lipid"))) %>%
			tibble::column_to_rownames(id)

	if (is.null(annot_cols)) {
	  annotation_col <- NA
	} else {
	  annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)
	}

	if (is.null(dots$annotation_colors)) {
		annotation_colors <- setHMColor(annotation_col)
	} else if (is.na(dots$annotation_colors)) {
		annotation_colors <- NA
	} else {
		annotation_colors <- eval(dots$annotation_colors, env = caller_env())
	}

	# complete_cases <- FALSE
	if (complete_cases) {
		df_hm <- df %>%
			dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
	} else {
		df_hm <- df %>%
			dplyr::mutate_if(is.numeric, funs(replace(., is.na(.), 0)))
	}

	df_hm <- df_hm %>%
			`rownames<-`(.[[id]])	%>%
			dplyr::select(which(names(.) %in% sample_ids))

	p <- my_pheatmap(
		mat = df_hm,
		filename = file.path(filepath, filename),
		annotation_col = annotation_col,
		annotation_row = annotation_row,
		color = mypalette,
		annotation_colors = annotation_colors,
		breaks = color_breaks,
		!!!dots)

	# to the benchmark group
	if(sum(!is.na(label_scheme_sub[[col_benchmark]])) > 0) {
		dfc <- ratio_toCtrl(df, !!id, label_scheme_sub, nm_ctrl = !!col_benchmark)

		if (complete_cases) {
			dfc_hm <- dfc %>%
				dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
		} else {
			dfc_hm <- dfc %>%
				dplyr::filter(rowSums(!is.na(.[, names(.) %in% sample_ids])) > 0) %>%
				dplyr::mutate_if(is.numeric, funs(replace(., is.na(.), 0)))
		}

		dfc_hm <- dfc_hm %>%
				`rownames<-`(.[[id]])	%>%
				dplyr::select(which(names(.) %in% sample_ids))

		p <- my_pheatmap(
			mat = dfc_hm,
			filename = file.path(filepath, paste0(fn_prx, "_toCtrl.", fn_suffix)),
			annotation_col = annotation_col,
			annotation_row = annotation_row,
			color = mypalette,
			annotation_colors = annotation_colors,
			breaks = color_breaks,
			!!!dots)
	}
}


#' Visualizes kinase heat maps
#'
#' Specialized for kinase heat maps
#'
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
#' @export
proteoKinHM <- function (id = gene, col_select = NULL, col_benchmark = NULL, scale_log2r = TRUE,
                         df = NULL, filepath = NULL, filename = NULL, complete_cases = FALSE,
                         impute_na = FALSE, anal_type = "Heatmap", xmin = -1, xmax = 1,
                         x_margin = 0.1, annot_cols = NULL, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

  col_select <- rlang::enexpr(col_select)
	col_benchmark <- rlang::enexpr(col_benchmark)

	if(is.null(col_select)) {
		col_select <- rlang::expr(Select)
		if(sum(!is.na(label_scheme$Select)) == 0)
			stop("Column \'", rlang::as_string(col_select),
			     "\' is either empty or missing from the \'label scheme\' file.\n",
				"\tEnter sample names under the \'", rlang::as_string(col_select),
				"\' column for informatic analysis.", call. = TRUE)
	} else if(col_select == rlang::expr(Sample_ID)) {
		stop(paste0("\'Sample_ID\' is a reserved column name.\n",
			"\tChoose or add a different column key to indicate samples for informatic analysis."),
			call. = TRUE)
	} else if(is.null(label_scheme[[col_select]])) {
		stop("Column \'", rlang::as_string(col_select), "\' is empty.\n",
			"\tEnter sample names under the \'", rlang::as_string(col_select),
			"\' column for informatic analysis.", call. = TRUE)
	}

	if(is.null(col_benchmark)) {
		col_benchmark <- rlang::expr(Benchmark)
		if(sum(!is.na(label_scheme$Benchmark)) == 0)
			warning("Default column name \'", rlang::as_string(col_benchmark),
			        "\' is empty or missing.", call. = TRUE)
	} else if(col_benchmark == rlang::expr(Sample_ID)) {
		stop(paste0("Column \'Sample_ID\' is reserved; choose a different key.\n"))
	} else if(is.null(label_scheme[[col_benchmark]])) {
		stop("Column name \'", rlang::as_string(col_benchmark), "\' is missing.", call. = TRUE)
	}

	id <- rlang::as_string(rlang::enexpr(id))
	id <- match_identifier(id)

	dots <- enexprs(...)

	if(!id %in% c("pep_seq_mod", "prot_acc", "gene"))
		stop("Unrecognized 'id'; needs to be one of pep_seq_mod, prot_acc or gene", call. = TRUE)

	if (id %in% c("prot_acc", "gene")) {
		data_type <- "Protein"
	} else if (id %in% c("pep_seq_mod")) {
		data_type <- "Peptide"
	}

	if(is.null(filepath)) {
		filepath = file.path(dat_dir, data_type, anal_type)
		dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
	}

	if(is.null(filename)) {
		fn_prx <- paste(data_type, anal_type, sep = "_")
		fn_suffix <- "png"
	} else {
		fn_prx <- gsub("\\..*$", "", filename)
		fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
	}

	if(is.null(df)) {
		err_msg <- "File doesn't exist"

		if (id %in% c("pep_seq", "pep_seq_mod")) {
			fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
			fn_imp <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
			fn_raw <- file.path(dat_dir, "Peptide", "Peptide.txt")
		} else if (id %in% c("prot_acc", "gene")) {
			fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
			fn_imp <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
			fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")
		}

		if(file.exists(fn_p)) {
			src_path <- fn_p
		} else {
			src_path <- ifelse(impute_na, fn_imp, fn_raw)
		}

		df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
		                        comment.char = "#"), error = function(e) NA)

		if(!is.null(dim(df))) {
			message(paste("File loaded:", gsub("\\\\", "/", src_path)))
		} else {
			stop(paste("No such file or directory:", gsub("\\\\", "/", src_path)))
		}
	}

	reload_expts()

	label_scheme_sub <- label_scheme %>%
		dplyr::select(Sample_ID, TMT_Set, !!col_select) %>%
		dplyr::filter(!is.na(!!col_select))

	dfw <- prepDM(df, !!rlang::sym(id), scale_log2r, label_scheme_sub$Sample_ID, anal_type = anal_type)

	dir.create(file.path(filepath, "Kinases"), recursive = TRUE, showWarnings = FALSE)

	dfw_kinase <- df %>%
	  dplyr::filter(kin_attr) %>%
	  prepDM(!!rlang::sym(id), scale_log2r, label_scheme_sub$Sample_ID, anal_type = anal_type)

	plotKinHM(
		df = df[df$kin_attr, ],
		id = !!rlang::sym(id),
		scale_log2r = scale_log2r,
		col_benchmark = !!col_benchmark,
		label_scheme_sub = label_scheme_sub,
		filepath = file.path(filepath, "Kinases"),
		filename = paste0(fn_prx, "_Kinases.", fn_suffix),
		complete_cases = complete_cases,
		xmin = xmin,
		xmax = xmax,
		x_margin = x_margin,
		annot_cols = annot_cols,
		!!!dots)
}
