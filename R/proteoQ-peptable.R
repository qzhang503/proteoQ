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
#' \code{normPep_Splex} prepares peptide data for each TMT experiment at one or
#' multiple LC/MS injections.
#'
#' @inheritParams normPep
#' @return Results in \code{.txt} files for each of TMT experiments and LC/MS
#'   injections.
#'
#' @examples
#' 	normPep_Splex(
#' 		id = pep_seq_mod,
#' 		method_psm_pep = median, 
#' 		group_pep_by = prot_acc, 
#' 	)
#'
#' \dontrun{
#' }
#' @import stringr dplyr purrr rlang  magrittr
normPep_Splex <- function (id = "pep_seq_mod", method_psm_pep = "median", group_pep_by = "prot_acc") {
	on.exit(message("Generation of individual peptide tables by TMT experiments --- Completed."),
	        add = TRUE)

	calcPepide <- function(df, label_scheme, id, method_psm_pep, set_idx, injn_idx) {
	  stopifnot("prot_acc" %in% names(df))
	  
	  id <- rlang::as_string(rlang::enexpr(id))
	  
	  channelInfo <- label_scheme %>%
	    dplyr::filter(TMT_Set == set_idx) %>%
	    channelInfo(set_idx)

	  df <- df[rowSums(!is.na(df[, grepl("^N_log2_R[0-9]{3}", names(df)), drop = FALSE])) > 0, ] %>%
	    dplyr::arrange(!!rlang::sym(id), prot_acc) %>%
	    dplyr::select(-grep("^R[0-9]{3}", names(.)))
	  
	  df <- df %>% 
	    dplyr::select(-which(names(.) %in% c(
	      "prot_hit_num", "prot_family_member", "prot_score", 
	      "prot_matches", "prot_sequences", 
	      "pep_var_mod", "pep_var_mod_pos", "pep_scan_title", 
	      "pep_res_before", "pep_res_after", 
	      "raw_file", "pep_query", "pep_summed_mod_pos", "pep_local_mod_pos")))
	  
	  df <- local({
	    col_start <- which(names(df) == "Modifications") + 1
	    col_end <- which(names(df) == "Charge") - 1

	    if (!(is_empty(col_start) | is_empty(col_end))) {
	      df <- df %>% dplyr::select(-(col_start : col_end))
	    }

	    return(df)
	  })
	  
    df <- df %>% 
      dplyr::select(-grep("\\s{1}Probabilities$", names(.))) %>% 
      dplyr::select(-grep("\\s{1}Score\\s{1}Diffs$", names(.))) %>% 
      dplyr::select(-which(names(.) %in% c(
        "Scan number", "Scan index", 
        "Deamidation (N) Probabilities", "Oxidation (M) Probabilities", 
        "Deamidation (N) Score Diffs", "Oxidation (M) Score Diffs", 
        "Acetyl (Protein N-term)", "Deamidation (N)", "Gln->pyro-Glu", "Oxidation (M)", 
        "Modifications", 
        "Fragmentation", "Mass analyzer", "Type", 
        "Scan event number", "Isotope index", "Simple mass error [ppm]", 
        "Retention time", "Delta score", "Score diff", 
        "Localization prob", "Precursor Full ScanNumber", 
        "Precursor Apex Offset", "Precursor Apex Offset Time", 
        "Matches", "Intensities", 
        "Mass Deviations [Da]", "Mass Deviations [ppm]", 
        "Masses", "Number of Matches", 
        "Neutral loss level", "ETD identification type", 
        "Reverse", "All scores", 
        "All sequences", "All modified sequences", 
        "Reporter PIF", "Reporter fraction", 
        "ID", "Protein group IDs", 
        "Peptide ID", "Mod. peptide ID", "Evidence ID", 
        "Length"))) %>% 
      dplyr::select(-grep("site IDs$", names(.)))

	  df <- df %>% 
	    dplyr::select(-which(names(.) %in% c(
	      "number", 
	      "variableSites", "nterm", "previous_aa", "sequence", "next_aa", 
	      "cys", "searchCycle", "L/H", "accession_numbers", "entry_name", 
	      "matched_parent_mass"))) 
	  
	  if ("pep_scan_title" %in% names(df)) {
	    df <- df %>% 
	    dplyr::mutate(pep_scan_title = gsub("\\\\", "~~", pep_scan_title)) %>%
	    dplyr::mutate(pep_scan_title = gsub("^File.*~~", "", pep_scan_title))	      
	  }

	  # summarise log2FC and intensity from the same `set_idx` at one or multiple LCMS injections
	  if (method_psm_pep == "mean") {
	    df_num <- aggrNums(mean)(df, !!rlang::sym(id), na.rm = TRUE)
	  } else if (method_psm_pep == "top.3") {
	    df_num <- TMT_top_n(df, !!rlang::sym(id), na.rm = TRUE)
	  } else if (method_psm_pep == "weighted.mean") {
	    df_num <- TMT_wt_mean(df, !!rlang::sym(id), na.rm = TRUE)
	  } else {
	    df_num <- aggrNums(median)(df, !!rlang::sym(id), na.rm = TRUE)
	  }
	  
	  df_first <- df %>% 
	    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
	    med_summarise_keys(id)

	  df <- list(df_first, df_num) %>%
	    purrr::reduce(left_join, by = id)
	  
	  if (id == "pep_seq_mod") {
	    df <- df %>% dplyr::select(-pep_seq)
	  } else {
	    df <- df %>% dplyr::select(-pep_seq_mod)
	  }

	  df <- cbind.data.frame(df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df)), drop = FALSE],
	                         df[, grepl("I[0-9]{3}", names(df)), drop = FALSE], 
	                         df[, grepl("log2_R[0-9]{3}", names(df)), drop = FALSE]) %>%
	    dplyr::mutate_at(.vars = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)),
	                     list(~ replace(.x, is.infinite(.x), NA)))
	  
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
#'aligned \code{log2_R...} according to \code{method_align} without scaling
#'normalization. Values under columns \code{Z_log2_R...} are \code{N_log2_R...}
#'with additional scaling normalization. Values under columns \code{I...} are
#'\code{reporter-ion intensity} before normalization. Values under columns
#'\code{N_I...} are normalized \code{I...}. Values under columns
#'\code{sd_log2_R...} are the standard deviation of the \code{log2FC} of
#'proteins from ascribing peptides.
#'
#'In general, median statistics is applied when summarising numeric values from
#'PSMs to peptides. One exception is \code{pep_expect} where geometric mean is
#'used.
#'
#'Also see \code{\link{normPrn}} for more description of the column keys in the
#'output.
#'
#'The peptide counts in individual peptide tables,
#'\code{TMTset1_LCMSinj1_Peptide_N.txt} et al., may be fewer than the entries
#'indicated under the \code{prot_n_pep} column after the peptide
#'removals/cleanups using \code{purgePSM}. Values under columns
#'\code{N_log2_R...} are intermediate reports by median-centering
#'\code{log2_R...} without scaling normalization.
#'
#'@param id Depreciated: character string; the variable for summarisation of
#'  \code{PSMs} into peptides. The option \code{id = pep_seq} corresponds to the
#'  summarisation by the primary sequences of peptides. The option \code{id =
#'  pep_seq_mod} corresponds to the summarisation by both the primary sequences
#'  and variable modifications of peptides. \code{PSMs} data with the same value
#'  in \code{pep_seq} or \code{pep_seq_mod} will be summarised into a single
#'  entry of peptide. NB: the value of \code{id} will match automatically to the
#'  value of \code{group_psm_by} in \code{normPSM}.
#'@param group_pep_by Depreciated: a character string for the grouping of
#'  peptide entries. At the \code{prot_acc} default, descriptive statistics will
#'  be calculated based on the same \code{prot_acc} groups. At
#'  \code{group_pep_by = gene}, proteins with the same gene name but different
#'  accession numbers will be treated as one group. NB: the value of
#'  \code{group_pep_by} will match automatically to the value of
#'  \code{group_pep_by} in \code{normPSM}.
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
#'  samples. At the \code{NULL} default, the column key of \code{Sample_ID} will
#'  be used, which results in the refit of the \code{log2FC} for all samples.
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
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/normPSM_examples.R
#'@example inst/extdata/examples/normPep_examples.R
#'@example inst/extdata/examples/purgePep_examples.R
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
	  # leave the Splex results under the `Peptide` folder, 
	  # so users can find the column keys for row filtration
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
		  dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
		  med_summarise_keys(id) %>% 
		  dplyr::select(-pep_n_psm, -prot_n_psm, -prot_n_pep, -TMT_Set) %>% # remove old values from single `TMT_Set`
		  dplyr::arrange(!!rlang::sym(id))  

		df <- list(pep_n_psm, df_first, df_num) %>%
			purrr::reduce(left_join, by = id)
		
		df <- list(df, prot_n_psm, prot_n_pep) %>%
		  purrr::reduce(left_join, by = group_pep_by)
		
		if (("pep_seq_mod" %in% names(df)) & (match_normPSM_par("use_lowercase_aa") %>% as.logical())) {
		  df <- df %>% 
		    dplyr::mutate(pep_mod_protnt = ifelse(grepl("^[A-z\\-]\\.~", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_protntac = ifelse(grepl("^[A-z\\-]\\._", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_pepnt = ifelse(grepl("^[A-z\\-]\\.[_~]?\\^", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_m = ifelse(grepl("m", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_n = ifelse(grepl("n", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_sty = ifelse(grepl("[sty]", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_pepct = ifelse(grepl("[\\^]{1}[_~]?\\.[A-z\\-]{1}$", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_protctam = ifelse(grepl("_{1}\\.[A-z\\-]{1}$", pep_seq_mod), TRUE, FALSE)) %>% 
		    dplyr::mutate(pep_mod_protct = ifelse(grepl("~{1}\\.[A-z\\-]{1}$", pep_seq_mod), TRUE, FALSE))
		}
		
		df <- dplyr::bind_cols(
		  df %>% select(grep("^pep_", names(.))), 
		  df %>% select(-grep("^pep_", names(.))), 
		)
		
		df <- dplyr::bind_cols(
		  df %>% select(grep("^prot_", names(.))), 
		  df %>% select(-grep("^prot_", names(.))), 
		)
		
		df <- dplyr::bind_cols(
		  df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
		  df %>% dplyr::select(grep("^I[0-9]{3}[NC]*", names(.))), 
		  df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]*", names(.))), 
		  df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
		  df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]*", names(.))), 
		  df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), 
		  df %>% dplyr::select(grep("^Z_log2_R[0-9]{3}[NC]*", names(.))), 
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
			filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0) 
		
		if (! id %in% names(df)) {
		  try(unlink(file.path(dat_dir, "Peptide\\Peptide.txt")))
		  try(unlink(file.path(dat_dir, "Protein\\Protein.txt")))
		  
		  stop("Column `", id, "` not found in `Peptide.txt`.", 
		       "\nDeleted the older `Peptide.txt` and try `normPep(...)` again.", 
		       call. = FALSE)
		}
		
		df <- df %>% 
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

	write.table(df, file.path(dat_dir, "Peptide", "Peptide.txt"), sep = "\t",
	            col.names = TRUE, row.names = FALSE)
	
	if (plot_log2FC_cv & TMT_plex(label_scheme) > 0) {
	  sd_violin(df, !!group_pep_by, 
	            file.path(dat_dir, "Peptide\\log2FC_cv\\raw", "Peptide_sd.png"), 
	            width = 8 * n_TMT_sets(label_scheme), height = 8, type = "log2_R", adjSD = FALSE, is_psm = FALSE)
	}
	
}


#' Summary of peptide keys by mean or geomean
#' 
#' @import dplyr purrr rlang  magrittr
med_summarise_keys <- function(df, id) {
  # Mascot keys or others converted to Mascot keys
  mascot_median_keys <- c("pep_score", "pep_rank", "pep_isbold", "pep_exp_mr", "pep_delta", 
                          "pep_exp_mz", "pep_exp_z")
  mascot_geomean_keys <- c("pep_expect")
  
  df_mascot_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mascot_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mascot_geomean <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mascot_geomean_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ my_geomean(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(mascot_median_keys, mascot_geomean_keys)))
  
  # MaxQuant keys
  df_mq_rptr_mass_dev <- df %>% 
    dplyr::select(!!rlang::sym(id), grep("^Reporter mass deviation", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep("^Reporter mass deviation", names(.)))	  
  
  mq_median_keys <- c(
    "Score", 
    "Charge", "Mass", "PIF", "Fraction of total spectrum", "Mass error [ppm]", 
    "Mass error [Da]", "Base peak fraction", "Precursor Intensity", 
    "Precursor Apex Fraction", "Intensity coverage", "Peak coverage", 
    "Combinatorics"
  )
  mq_geomean_keys <- c("PEP")
  
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mq_geomean <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_geomean_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ my_geomean(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(mq_median_keys, mq_geomean_keys)))
  
  # Spectrum Mill keys
  sm_median_keys <- c(
    "deltaForwardReverseScore", "percent_scored_peak_intensity", "totalIntensity", 
    "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
    "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
    "delta_parent_mass", "delta_parent_mass_ppm")
  
  df_sm_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% sm_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% sm_median_keys))
  
  # put together
  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id)))
  
  df <- list(df_first, 
             df_mascot_med, df_mascot_geomean, 
             df_mq_rptr_mass_dev, df_mq_med, df_mq_geomean, 
             df_sm_med) %>%
    purrr::reduce(left_join, by = id) %>%
    data.frame(check.names = FALSE)
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
}


#' Median-centering normalization
#' 
#' # not currently used?
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


