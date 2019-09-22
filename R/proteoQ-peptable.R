#' Make new column names
#'
#' \code{newColnames} match names to Sample_ID in label_scheme
#'
#' @param i Integer; the index of TMT experiment 
#' @param x Data frame; log2FC data
#' @param label_scheme Experiment summary
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
newColnames <- function(i, x, label_scheme) {
  label_scheme_sub <- label_scheme %>%
    dplyr::filter(TMT_Set == i)
  
  cols <- grep(paste0("[RI][0-9]{3}[NC]*_", i, "$"), names(x))
  nm_channel <- gsub(paste0("([RI][0-9]{3}[NC]*)_", i, "$"), "\\1", names(x)[cols])
  names(x)[cols] <- paste0(nm_channel, " (", as.character(label_scheme_sub$Sample_ID), ")")
  
  cols <- grep("[RI][0-9]{3}.*\\s+\\(.*\\)$", names(x))
  
  # cols with new names go first
  if (length(cols) < ncol(x)) x <- dplyr::bind_cols(x[, cols], x[, -cols, drop = FALSE])
  
  return(x)
}


#' Peptide reports for individual TMT experiments
#'
#' \code{normPep_Splex} prepares peptide data for each TMT experiment at
#' different LC/MS injections.
#'
#' @inheritParams normPep
#' @return Results in \code{.txt} files for each of TMT experiments and LC/MS
#'   injections.
#'
#' @examples
#' 	normPep_Splex(
#' 		id = pep_seq_mod,
#' 		method_psm_pep = median
#' 	)
#'
#' \dontrun{
#' }
#' @import stringr dplyr purrr rlang  magrittr
normPep_Splex <- function (id = "pep_seq_mod", method_psm_pep = "median", group_pep_by = "prot_acc") {
	on.exit(message("Generation of individual peptide tables by RAW filenames --- Completed."),
	        add = TRUE)

	calcPepide <- function(df, label_scheme, id, method_psm_pep, set_idx, injn_idx) {
	  id <- rlang::as_string(rlang::enexpr(id))
	  
	  channelInfo <- label_scheme %>%
	    dplyr::filter(TMT_Set == set_idx) %>%
	    channelInfo(set_idx)
	  
	  df <- df[rowSums(!is.na(df[, grepl("^N_log2_R[0-9]{3}", names(df)), drop = FALSE])) > 0, ] %>%
	    dplyr::arrange(!!rlang::sym(id), prot_acc) %>%
	    # dplyr::filter(pep_isunique == 1) %>%
	    dplyr::select(-grep("^R[0-9]{3}", names(.))) %>%
	    dplyr::mutate(pep_scan_title = gsub("\\\\", "~~", pep_scan_title)) %>%
	    dplyr::mutate(pep_scan_title = gsub("^File.*~~", "", pep_scan_title))
	  
	  # summarise log2FC and intensity from the same `set_idx` but different LCMS injections
	  if (method_psm_pep == "mean") {
	    df_num <- aggrNums(mean)(df, !!rlang::sym(id), na.rm = TRUE)
	  } else if (method_psm_pep == "top.3") {
	    df_num <- TMT_top_n(df, !!rlang::sym(id), na.rm = TRUE)
	  } else if (method_psm_pep == "weighted.mean") {
	    df_num <- TMT_wt_mean(df, !!rlang::sym(id), na.rm = TRUE)
	  } else {
	    df_num <- aggrNums(median)(df, !!rlang::sym(id), na.rm = TRUE)
	  }

	  df_score <- df %>% 
	    dplyr::select(!!rlang::sym(id), pep_score) %>% 
	    dplyr::group_by(!!rlang::sym(id)) %>% 
	    dplyr::summarise(pep_score = max(pep_score, na.rm = TRUE))
	  
	  df_expect <- df %>% 
	    dplyr::select(!!rlang::sym(id), pep_expect) %>% 
	    dplyr::group_by(!!rlang::sym(id)) %>% 
	    dplyr::summarise(pep_expect = min(pep_expect, na.rm = TRUE))
	  
	  df_rank <- df %>% 
	    dplyr::select(!!rlang::sym(id), pep_rank) %>% 
	    dplyr::group_by(!!rlang::sym(id)) %>% 
	    dplyr::summarise(pep_rank = min(pep_rank, na.rm = TRUE))
	  
	  df_isbold <- df %>% 
	    dplyr::select(!!rlang::sym(id), pep_isbold) %>% 
	    dplyr::group_by(!!rlang::sym(id)) %>% 
	    dplyr::summarise(pep_isbold = max(pep_isbold, na.rm = TRUE))
	  
	  df_exp_mr <- df %>% 
	    dplyr::select(!!rlang::sym(id), pep_exp_mr) %>% 
	    dplyr::group_by(!!rlang::sym(id)) %>% 
	    dplyr::summarise(pep_exp_mr = median(pep_exp_mr, na.rm = TRUE))
	  
	  df_exp_z <- df %>% 
	    dplyr::select(!!rlang::sym(id), pep_exp_z) %>% 
	    dplyr::group_by(!!rlang::sym(id)) %>% 
	    dplyr::summarise(pep_exp_z = median(pep_exp_z, na.rm = TRUE))

	  df_first <- df %>% 
	    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% # already in `df_num`
	    dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score",
	                                         "prot_matches", "prot_matches_sig", "prot_sequences", "prot_sequences_sig", 
	                                         "pep_score", "pep_expect", "pep_rank", "pep_isbold", 
	                                         "pep_exp_mr", "pep_exp_z", "pep_exp_mz", "pep_delta", "pep_var_mod", 
	                                         "pep_var_mod_pos", "pep_scan_title", "pep_res_before", "pep_res_after", 
	                                         "raw_file", "pep_query", "pep_summed_mod_pos", "pep_local_mod_pos"))) %>% 
	    dplyr::filter(!duplicated(!!rlang::sym(id)))
	  
	  df_first <- dplyr::bind_cols(
	    df_first %>% dplyr::select(-pep_calc_mr), 
	    df_first %>% dplyr::select(pep_calc_mr), 
	  )
	  
	  if(id == "pep_seq_mod") {
	    df_first <- df_first %>% dplyr::select(-pep_seq)
	  } else {
	    df_first <- df_first %>% dplyr::select(-pep_seq_mod)
	  }
	  
	  df <- list(df_first, df_exp_mr, df_exp_z, df_score, df_expect, df_rank, df_isbold, 
	             df_num) %>%
	    purrr::reduce(left_join, by = id) %>%
	    data.frame(check.names = FALSE)
	  
	  df <- dplyr::bind_cols(
	    df %>% dplyr::select(grep("^pep_", names(.))), 
	    df %>% dplyr::select(-grep("^pep_", names(.))), 
	  )

	  df <- cbind.data.frame(df[, !grepl("[I|log2_R][0-9]{3}", names(df)), drop = FALSE],
	                         df[, grepl("[log2_R][0-9]{3}", names(df)), drop = FALSE],
	                         df[, grepl("[I][0-9]{3}", names(df)), drop = FALSE]) %>%
	    dplyr::mutate_at(.vars = grep("[I|log2_R][0-9]{3}", names(.)),
	                     list(~ replace(.x, is.infinite(.x), NA) ))
	  
	  # median-centered across TMT channels under the same multiplex experiment
	  col_r <- grepl("^log2_R[0-9]{3}", names(df))
	  cf <- apply(df[, col_r, drop = FALSE], 2, median, na.rm = TRUE)
	  df <- cbind(df[, -grep("^N_log2_R[0-9]{3}", names(df))],
	              sweep(df[, col_r], 2, cf, "-") %>%
	                `colnames<-`(paste("N", colnames(.), sep="_")))
	  
	  col_int <- grepl("^I[0-9]{3}", names(df))
	  df  <- cbind(df[, -grep("^N_I[0-9]{3}", names(df))],
	               sweep(df[, col_int, drop=FALSE], 2, 2^cf, "/") %>%
	                 `colnames<-`(paste("N", colnames(.), sep="_")))
	  
	  rm(cf, col_r, col_int)
	  
	  df <- df %>% 
	    calcSD_Splex(group_pep_by) %>% 
	    `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
	    dplyr::right_join(df, by = group_pep_by) %>% 
	    na_zeroIntensity() %>% 
	    dplyr::mutate(TMT_Set = set_idx)

	  return(df)
	}
	
	
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 2000000)

	id <- rlang::as_string(rlang::enexpr(id))
	method_psm_pep <- rlang::as_string(rlang::enexpr(method_psm_pep))

	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
		reorder_files(n_TMT_sets(label_scheme_full))

	purrr::walk(as.list(filelist), ~ {
		fn_prx <- gsub("_PSM_N.txt", "", .x, fixed = TRUE)

		set_idx <- as.integer(gsub(".*TMTset(\\d+)_.*", "\\1", fn_prx))
		injn_idx <- as.integer(gsub(".*LCMSinj(\\d+).*", "\\1", fn_prx))

		df <- read.csv(file.path(dat_dir, "PSM", .x), check.names = FALSE, header = TRUE,
		               sep = "\t", comment.char = "#") %>% 
		  dplyr::select(-grep("^sd_log2_R", names(.))) %>% 
		  calcPepide(label_scheme = label_scheme, id = !!id, method_psm_pep = method_psm_pep, 
		             set_idx = set_idx, injn_idx = injn_idx)
		
		df <- dplyr::bind_cols(
		  df %>% dplyr::select(grep("^pep_", names(.))), 
		  df %>% dplyr::select(-grep("^pep_", names(.))), 
		)
		
		df <- dplyr::bind_cols(
		  df %>% dplyr::select(grep("^prot_", names(.))),
		  df %>% dplyr::select(-grep("^prot_", names(.))),
		)

		write.table(df, file.path(dat_dir, "Peptide", paste0(fn_prx, "_Peptide_N", ".txt")),
		            sep = "\t", col.names = TRUE, row.names = FALSE)
	})
}


#'Reports peptide results
#'
#'\code{normPep} summarises
#'\code{\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSMs}}
#'into peptides and normalizes the data across
#'\code{\href{https://en.wikipedia.org/wiki/Tandem_mass_tag}{TMT}} experiments
#'and \code{LC/MS} injections.
#'
#'In the primary output file, "\code{Peptide.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'\code{log2_R...} with median-centering alignment. Values under columns
#'\code{Z_log2_R...} are \code{N_log2_R...} with scaling normalization. Values
#'under columns \code{I...} are \code{reporter-ion intensity} before
#'normalization. Values under columns \code{N_I...} are normalized \code{I...}.
#'Values under columns \code{sd_log2_R...} are the standard deviation of the
#'\code{log2FC} of proteins from ascribing peptides.
#'
#'Also see \code{\link{normPrn}} for more description of the column keys in the
#'output.
#'
#'@param id Character string; the variable for summarisation of \code{PSMs} into
#'  peptides. The option \code{id = pep_seq} corresponds to the summarisation by
#'  the primary sequences of peptides. The option \code{id = pep_seq_mod}
#'  corresponds to the summarisation by both the primary sequences and variable
#'  modifications of peptides. \code{PSMs} data with the same value in
#'  \code{pep_seq} or \code{pep_seq_mod} will be summarised into a single entry
#'  of peptide. The value of \code{id} will match automatically to the value of
#'  \code{group_psm_by} in \code{normPSM}.
#'@param group_pep_by A character string for the grouping of peptide entries. At
#'  the \code{prot_acc} default, descriptive statistics will be calculated based
#'  on the same \code{prot_acc} groups. At \code{group_pep_by = gene}, proteins
#'  with the same gene name but different accession numbers will be treated as
#'  one group. The value of \code{group_pep_by} will match automatically to the
#'  value of \code{group_pep_by} in \code{normPSM}.
#'@param method_psm_pep Character string; the method to summarise the
#'  \code{log2FC} and the \code{intensity} of \code{PSMs} by peptide entries.
#'  The descriptive statistics includes \code{c("mean", "median", "top.3",
#'  "weighted.mean")}. The \code{log10-intensity} of reporter ions at the
#'  \code{PSMs} levels will be the weight when summarising \code{log2FC} with
#'  \code{"top.3"} or \code{"weighted.mean"}.
#'@param method_align Character string or a list of gene symbols; the method to
#'  align the \code{log2FC} of peptide/protein entries across samples.
#'  \code{MC}: median-centering; \code{MGKernel}: the kernal density defined by
#'  multiple Gaussian functions (\code{\link[mixtools]{normalmixEM}}). At
#'  \code{method_align = "MC"}, the ratio profiles of each sample will be
#'  aligned in that the medians of the \code{log2FC} are zero. At
#'  \code{method_align = "MGKernel"}, the \code{log2FC} will be aligned in that
#'  the maximums of kernel density are zero. It is also possible to align the
#'  \code{log2FC} to the median of a list of user-supplied genes:
#'  \code{method_align = c("ACTB", "GAPDH", ...)}.
#'@param range_log2r Numeric vector at length two; the range of the
#'  \code{log2FC} of peptide/protein entries for use in the scaling
#'  normalization of standard deviation across samples. The default is between
#'  the 10th and the 90th quantiles.
#'@param range_int Numeric vector at length two; the range of the
#'  \code{intensity} of reporter ions for use in the scaling normalization of
#'  standard deviation across samples. The default is between the 5th and the
#'  95th quantiles.
#'@param n_comp Integer; the number of Gaussian components to be used with
#'  \code{method_align = "MGKernel"}. A typical value is 2 or 3. The variable
#'  \code{n_comp} overwrites the augument \code{k} in
#'  \code{\link[mixtools]{normalmixEM}}.
#'@param seed Integer; a seed for reproducible fitting at \code{method_align =
#'  MGKernel}.
#'@param col_refit Character string to a column key in \code{expt_smry.xlsx}.
#'  Samples corresponding to non-empty entries under \code{col_refit} will be
#'  used in the refit of \code{log2FC} using multiple Gaussian kernels. The
#'  density estimates from an earlier analyis will be kept for the remaining
#'  samples. At the NULL default, the column key \code{Sample_ID} will be used,
#'  which results in the refit of the \code{log2FC} for all samples.
#'@param cache Logical; if TRUE, use cache.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr Additional parameters for
#'  \code{\link[mixtools]{normalmixEM}}:\cr \code{maxit}, the maximum number of
#'  iterations allowed; \cr \code{epsilon}, tolerance limit for declaring
#'  algorithm convergence.
#'@inheritParams normPSM
#'@inheritParams mixtools::normalmixEM
#'@seealso \code{\link{normPSM}} for PSM an \code{\link{normPrn}} for proteins.
#'
#'@return The primary output is in \code{~\\dat_dir\\Peptide\\Peptide.txt}.
#'
#' @examples
#' \dontrun{
#' # peptides results with examplary `filter_...`
#' normPep(
#'   method_psm_pep = median,
#'   method_align = MGKernel,
#'   range_log2r = c(5, 95),
#'   range_int = c(5, 95),
#'   n_comp = 3,
#'   seed = 749662,
#'   maxit = 200,
#'   epsilon = 1e-05,
#'
#'   filter_by = exprs(pep_n_psm >= 2),
#'   # filter_by_sp = exprs(species == "human"),
#' )
#'
#' # examplary peptide purging; n: the number of peptides
#' purgePep(max_cv = .5, min_n = 2)
#'
#' }
#'@import stringr dplyr tidyr purrr data.table rlang
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@importFrom plyr ddply
#'@export
normPep <- function (id = c("pep_seq", "pep_seq_mod"),
										method_psm_pep = c("median", "mean", "weighted.mean", "top.3", "mqpep"), 
										group_pep_by = c("prot_acc", "gene"), 
										method_align = c("MC", "MGKernel"), range_log2r = c(10, 90),
										range_int = c(5, 95), n_comp = NULL, seed = NULL, 
										annot_kinases = FALSE, 
										col_refit = NULL, cache = TRUE, 
										plot_log2FC_cv = TRUE, ...) {

	dir.create(file.path(dat_dir, "Peptide\\cache"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Peptide\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Peptide\\log2FC_cv\\raw"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Peptide\\log2FC_cv\\purged"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Peptide\\log"), recursive = TRUE, showWarnings = FALSE)
	
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 2000000)
	
	reload_expts()

	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	# depreciated; instead matched to normPSM()
	id <- rlang::enexpr(id)
	if (id == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
	  id <- "pep_seq_mod"
	} else {
	  id <- rlang::as_string(id)
	  stopifnot(id %in% c("pep_seq", "pep_seq_mod"))
	}
	id <- match_normPSM_pepid()
	
	# depreciated; instead matched to normPSM()
	group_pep_by <- rlang::enexpr(group_pep_by)
	if (group_pep_by == rlang::expr(c("prot_acc", "gene"))) {
	  group_pep_by <- "prot_acc"
	} else {
	  group_pep_by <- rlang::as_string(group_pep_by)
	  stopifnot(group_pep_by %in% c("prot_acc", "gene"))
	}
	group_pep_by <- match_normPSM_protid()
	
	col_refit <- rlang::enexpr(col_refit)
	col_refit <- ifelse(is.null(col_refit), rlang::expr(Sample_ID), rlang::sym(col_refit))
	
	if (is.null(label_scheme[[col_refit]])) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("Column \'", rlang::as_string(col_refit), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_refit]])) == 0) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("No samples were specified under column \'", rlang::as_string(col_refit), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
	}
	
	method_psm_pep <- rlang::enexpr(method_psm_pep)
	if (method_psm_pep == rlang::expr(c("median", "mean", "weighted.mean", "top.3", "mqpep"))) {
		method_psm_pep <- "median"
	} else {
		method_psm_pep <- rlang::as_string(method_psm_pep)
		stopifnot(method_psm_pep %in% c("mean", "top.3", "median", "weighted.mean", "mqpep"))
	}

	method_align <- rlang::enexpr(method_align)
	if (method_align == rlang::expr(c("MC", "MGKernel"))) {
		method_align <- "MC"
	} else {
		method_align <- rlang::as_string(method_align)
		if(!method_align %in% c("MC", "MGKernel"))
			warning("Assume the value of 'method_align' is a list of housekeeping proteins")
	}
	
	stopifnot(length(range_log2r) == 2)
	stopifnot(length(range_int) == 2)
	
	if (range_log2r[2] <= 1) range_log2r <- range_log2r * 100
	if (range_int[2] <= 1) range_int <- range_int * 100
	
	stopifnot(range_log2r[1] < range_log2r[2] & 
	            range_log2r[1] >= 0 & range_log2r[2] <= 100)
	
	stopifnot(range_int[1] < range_int[2] & 
	            range_int[1] >= 0 & range_int[2] <= 100)
	
	if (is.null(n_comp)) n_comp <- if(nrow(df) > 3000) 3L else 2L
	n_comp <- n_comp %>% as.integer()
	stopifnot(n_comp >= 2)
	
	dots <- rlang::enexprs(...)
	lang_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	dots <- dots %>% .[! . %in% lang_dots]

	mget(names(formals()), rlang::current_env()) %>% save_call("normPep")

	if (!(cache & file.exists(file.path(dat_dir, "Peptide", "Peptide.txt")))) {
	  normPep_Splex(id = !!id, method_psm_pep = method_psm_pep, group_pep_by = group_pep_by)

		df <- do.call(rbind,
			lapply(list.files(path = file.path(dat_dir, "Peptide"), 
			                  pattern = paste0("TMTset[0-9]+_LCMSinj[0-9]+_Peptide_N\\.txt$"),
			                  full.names = TRUE), read.csv, check.names = FALSE, header = TRUE,
			       sep = "\t", comment.char = "#")) %>%
			dplyr::mutate(TMT_Set = factor(TMT_Set)) %>%
			dplyr::arrange(TMT_Set) 

		cat("Available column keys in PSM tables for data filtration: \n")
		cat(paste0(names(df), "\n"))
		df <- df %>% filters_in_call(!!!lang_dots)

		dup_peps <- df %>%
			dplyr::select(!!rlang::sym(id), prot_acc) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise(N = n_distinct(prot_acc)) %>%
			dplyr::filter(N > 1)

		if (nrow(dup_peps) > 0) {
		  df <- df %>% dplyr::filter(! (!!rlang::sym(id) %in% dup_peps[[id]]))
		  write.csv(dup_peps, file.path(dat_dir, "Peptide", "dbl_dipping_peptides.csv"), row.names = FALSE)
		}
		
		write.csv(df, file.path(dat_dir, "Peptide\\cache", "unambi_peptides.csv"), row.names = FALSE)
		rm(dup_peps)

		# median summarisation of data from the same TMT experiment at different LCMS injections
		df_num <- df %>% 
		  dplyr::select(!!rlang::sym(id), 
		                TMT_Set, 
		                grep("^sd_log2_R[0-9]{3}", names(.)), 
		                grep("^log2_R[0-9]{3}", names(.)), 
		                grep("^N_log2_R[0-9]{3}", names(.)), 
		                grep("^Z_log2_R[0-9]{3}", names(.)), 
		                grep("^I[0-9]{3}", names(.)), 
		                grep("^N_I[0-9]{3}", names(.))) %>% 
		  dplyr::group_by(!!rlang::sym(id), TMT_Set) %>%
		  dplyr::summarise_all(~ median(.x, na.rm = TRUE))

		df_num <- df_num %>%
			dplyr::arrange(TMT_Set) %>%
			tidyr::gather(grep("R[0-9]{3}|I[0-9]{3}", names(.)), key = ID, value = value) %>%
			tidyr::unite(ID, ID, TMT_Set)

		# define the levels of TMT channels;
		# otherwise, the order of channels will flip between N(itrogen) and C(arbon)
		Levels <- unique(df_num$ID)
		df_num <- df_num %>%
			dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
			tidyr::spread(ID, value)
		rm(Levels)

		for (set_idx in seq_len(n_TMT_sets(label_scheme))) {
		  df_num <- newColnames(set_idx, df_num, label_scheme)
		}
		
		df_num <- df_num %>% 
		  dplyr::select(!!rlang::sym(id), grep("[RI][0-9]{3}[NC]*", names(.))) %>% 
		  dplyr::arrange(!!rlang::sym(id)) %T>%
		  write.csv(file.path(dat_dir, "Peptide\\cache", "pep_num.csv"), row.names = FALSE)

		# summed `pep_n_psm`, `prot_n_psm` and `prot_n_pep` after data merging
		pep_n_psm <- df %>%
		  dplyr::select(!!rlang::sym(id), pep_n_psm) %>%
		  dplyr::group_by(!!rlang::sym(id)) %>%
		  dplyr::summarise(pep_n_psm = sum(pep_n_psm)) %>% 
		  dplyr::arrange(!!rlang::sym(id))
		
		prot_n_psm <- df %>% 
		  dplyr::select(pep_n_psm, !!rlang::sym(group_pep_by)) %>% 
		  dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
		  dplyr::summarise(prot_n_psm = sum(pep_n_psm))

		prot_n_pep <- df %>% 
		  dplyr::select(!!rlang::sym(id), !!rlang::sym(group_pep_by)) %>% 
		  dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
		  dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
		  dplyr::summarise(prot_n_pep = n())

		df_first <- df %>% 
		  dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
		  dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
		  dplyr::select(-pep_n_psm, -prot_n_psm, -prot_n_pep, -TMT_Set) %>% # remove old values from single `TMT_Set`
		  dplyr::arrange(!!rlang::sym(id))

		df <- list(pep_n_psm, df_first, df_num) %>%
			purrr::reduce(left_join, by = id)
		
		df <- list(df, prot_n_psm, prot_n_pep) %>%
		  purrr::reduce(left_join, by = group_pep_by)
		
		df <- dplyr::bind_cols(
		  df %>% select(grep("^pep_", names(.))), 
		  df %>% select(-grep("^pep_", names(.))), 
		)
		
		df <- dplyr::bind_cols(
		  df %>% select(grep("^prot_", names(.))), 
		  df %>% select(-grep("^prot_", names(.))), 
		)

		df <- df %>% 
		  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
		  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
		  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
		  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
		  dplyr::mutate_at(vars(grep("sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3))

		df <- df %>%
			dplyr::filter(!duplicated(.[[id]])) %>% 
		  dplyr::filter(rowSums(!is.na(.[, grepl("N_log2_R", names(.))])) > 0) %>% 
		  dplyr::arrange(!!rlang::sym(id)) %T>%
		  write.csv(file.path(dat_dir, "Peptide\\cache", "Peptide_no_norm.csv"), row.names = FALSE)		  
	} else {
		df <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"),
			check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>%
			filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0) %>% 
		  dplyr::arrange(!!rlang::sym(id))
	}
	
	quietly_out <- purrr::quietly(normMulGau)(
	  df = df,
	  method_align = method_align,
	  n_comp = n_comp,
	  seed = seed,
	  range_log2r = range_log2r,
	  range_int = range_int,
	  filepath = file.path(dat_dir, "Peptide\\Histogram"),
	  col_refit = col_refit, 
	  !!!dots
	)
	
	purrr::walk(quietly_out[-1], write, 
	            file.path(dat_dir, "Peptide\\log","pep_MulGau_log.csv"), append = TRUE)
	
	df <- quietly_out$result

	df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-  
	  df[, grepl("^Z_log2_R[0-9]{3}", names(df))] %>%
	  dplyr::mutate_if(is.logical, as.numeric) %>%
	  round(., digits = 3)
	
	df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-  
	  df[, grepl("^N_log2_R[0-9]{3}", names(df))] %>%
	  dplyr::mutate_if(is.logical, as.numeric) %>%
	  round(., digits = 3)
	
	df[, grepl("I[0-9]{3}", names(df))] <-  
	  df[, grepl("I[0-9]{3}", names(df))] %>%
		dplyr::mutate_if(is.logical, as.numeric) %>%
		round(., digits = 0)
	
	suppressWarnings(
	  df <- df %>% 
	    dplyr::select(-one_of(
	      "m/z", "Scan number", "Scan index", "Length", "Deamidation (N) Probabilities", 
	      "Oxidation (M) Probabilities", "Deamidation (N) Score Diffs", "Oxidation (M) Score Diffs", 
	      "Acetyl (Protein N-term)", "Deamidation (N)", "Glu->pyro-Glu", "Oxidation (M)", 
	      "Charge", "Fragmentation", "Mass analyzer", "Type", "Scan event number", 
	      "Isotope index", "Mass", "Mass error [ppm]", "Mass error [Da]", 
	      "Simple mass error [ppm]", "Retention time", "Delta score", "Score diff", 
	      "Localization prob", "Combinatorics", "PIF", 
	      "Fraction of total spectrum", "Base peak fraction", "Precursor Full ScanNumber", 
	      "Precursor Intensity", "Precursor Apex Fraction", "Precursor Apex Offset", 
	      "Precursor Apex Offset Time", "Matches", "Intensities", "Mass Deviations [Da]", 
	      "Mass Deviations [ppm]", "Masses", "Number of Matches", "Intensity coverage", "Peak coverage", 
	      "Neutral loss level", "ETD identification type", "Reverse", "All scores", "All sequences", 
	      "All modified sequences", "Reporter PIF", "Reporter fraction", "ID", "Protein group IDs", 
	      "Peptide ID", "Mod. peptide ID", "Evidence ID", "Deamidation (N) site IDs", 
	      "Oxidation (M) site IDs"
	    )) %>% 
	    dplyr::select(-grep("^Reporter mass deviation", names(.)))
	)

	write.table(df, file.path(dat_dir, "Peptide", "Peptide.txt"), sep = "\t",
	            col.names = TRUE, row.names = FALSE)
	
	if (plot_log2FC_cv & TMT_plex(label_scheme) > 0) {
	  df %>% 
	    dplyr::select(group_pep_by, grep("^sd_log2_R[0-9]{3}", names(.))) %>% 
	    dplyr::filter(!duplicated(.[[group_pep_by]])) %>% 
	    dplyr::filter(rowSums(!is.na(.[grep("^sd_log2_R[0-9]{3}", names(.))])) > 0) %>% 
	    sd_violin(!!group_pep_by, file.path(dat_dir, "Peptide\\log2FC_cv\\raw", "Peptide_sd.png"), 
	              8 * n_TMT_sets(label_scheme), 8)
	}
	
}




#' Median-centering normalization
#' 
#' @import dplyr purrr rlang  magrittr
logfcPep <- function(df, label_scheme, set_idx) {
  label_scheme_sub <- label_scheme %>% 
    .[.$TMT_Set == set_idx & .$LCMS_Injection == 1, ]
  
  channelInfo <- channelInfo(label_scheme_sub, set_idx)
  
  col_sample <- grep("^I[0-9]{3}", names(df))
  
  if (length(channelInfo$refChannels) > 0) {
    ref_index <- channelInfo$refChannels
  } else {
    ref_index <- channelInfo$labeledChannels
  }
  
  df <- sweep(df[, col_sample, drop = FALSE], 1,
              rowMeans(df[, col_sample[ref_index], drop = FALSE], na.rm = TRUE), "/") %>%
    log2(.) %>%
    `colnames<-`(gsub("I", "log2_R", names(.)))	%>%
    cbind(df, .) %>%
    dplyr::mutate_at(.vars = grep("[I|R][0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA))
  
  col_log2Ratio <- grepl("^log2_R[0-9]{3}", names(df))
  cf <- apply(df[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  
  df <- sweep(df[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(df, .)
  
  df <- sweep(df[, grepl("^I[0-9]{3}", names(df)), drop = FALSE], 2, 2^cf, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(df, .)
  
  df <- df %>%
    na_zeroIntensity()
}


