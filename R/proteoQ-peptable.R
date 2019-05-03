#' Peptide results for individual TMT experiments
#'
#' \code{normPep_Splex} prepares peptide data for each TMT experiment at
#' different LC/MS injections.
#'
#' @param id The variable for summarising PSMs into peptides. PSMs with the same
#'   value in \code{id} will be summarised into a single entry of peptide.
#' @param method_psm_pep The method for summarising \code{log2-ratios} from PSMs
#'   to peptides. \code{log10-intensities} will be the weights if
#'   \code{method_psm_pep = c("top.3", "weighted.mean")}.
#' @return Results in \code{.txt} files for each of TMT experiments and LC/MS
#'   injections.
#'
#' @examples
#' 	normPep_Splex(
#' 		id = pep_seq_mod,
#' 		method_psm_pep = "median"
#' 	)
#'
#' \dontrun{
#' }
#' @import stringr dplyr purrr rlang
#' @importFrom magrittr %>%
normPep_Splex <- function (id = "pep_seq_mod", method_psm_pep = "median") {

	dir.create(file.path(dat_dir, "Peptide\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	on.exit(message("Generation of individual peptide tables by RAW filenames --- Completed."),
	        add = TRUE)

	calcPepide <- function(df, label_scheme, id, method_psm_pep,
	                       set_idx, injn_idx, range_log2r = c(0, 100), range_int = c(30, 80)) {

	  # Not used: injn_idx
		# Not used: range_log2r
		# Not used: range_int

		id <- rlang::as_string(rlang::enexpr(id))

		channelInfo <- label_scheme %>%
			dplyr::filter(TMT_Set == set_idx) %>%
			channelInfo(set_idx)

		df <- df[rowSums(!is.na(df[, grepl("^N_log2_R[0-9]{3}", names(df)), drop = FALSE])) > 0, ] %>%
				dplyr::arrange(!!rlang::sym(id), prot_acc) %>%
				dplyr::filter(pep_isunique == 1) %>%
				dplyr::select(-grep("^R[0-9]{3}", names(.))) %>%
				dplyr::mutate(pep_scan_title = gsub("\\\\", "~~", pep_scan_title)) %>%
				dplyr::mutate(pep_scan_title = gsub("^File.*~~", "", pep_scan_title))

		# summarise log2-ratios and intensity from the same TMT experiment at different LCMS injections
		if (method_psm_pep == "mean") {
			df_num <- aggrNums(mean)(df, !!rlang::sym(id), na.rm = TRUE)
		} else if (method_psm_pep == "top.3") {
			df_num <- TMT_top_n(df, !!rlang::sym(id), na.rm = TRUE)
		} else if (method_psm_pep == "weighted.mean") {
			df_num <- TMT_wt_mean(df, !!rlang::sym(id), na.rm = TRUE)
		} else {
			df_num <- aggrNums(median)(df, !!rlang::sym(id), na.rm = TRUE)
		}

		df_psm <- df %>%
				dplyr::select(!!rlang::sym(id)) %>%
				dplyr::group_by(!!rlang::sym(id)) %>%
				dplyr::summarise(n_psm = n())


		# summary of categories using all strings
		df_str_a <- df %>%
      # already processed
		  dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
		  # will not be used
      dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score",
                                           "prot_matches", "prot_matches_sig", "prot_sequences",
                                           "prot_sequences_sig"))) %>%
      # redundant data available in "pep_seq_mod"; "pep_var_mod_pos" but kept for convenience
		  dplyr::select(-which(names(.) %in% c("pep_res_before", "pep_res_after",
		                                       "pep_summed_mod_pos", "pep_local_mod_pos"))) %>%
		  # will be processed later
      dplyr::select(-which(names(.) %in% c("prot_acc", "pep_start", "pep_end", "pep_miss",
                                           "pep_var_mod", "pep_var_mod_pos", "prot_desc",
                                           "prot_mass")))

		# will also be processed later
		if(id == "pep_seq_mod") df_str_a <- df_str_a %>% dplyr::select(-pep_seq)

		df_str_a <- df_str_a %>%
				dplyr::group_by(!!rlang::sym(id)) %>%
				dplyr::summarise_all(toString, na.rm = TRUE)

		# summary of categories using only the first string
		nm_a <- names(df_str_a)[-which(names(df_str_a) == id)]
		df_str_b <- df %>%
				dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% # already processed
				dplyr::select(-which(names(.) %in% nm_a)) %>% # already processed under "df_str_a"
				dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score",
				                                     "prot_matches", "prot_matches_sig", "prot_sequences",
				                                     "prot_sequences_sig"))) %>% # will not be used
				dplyr::select(-which(names(.) %in% c("pep_res_before", "pep_res_after", "pep_summed_mod_pos",
				                                     "pep_local_mod_pos"))) %>% # redundant kept for convenience
				dplyr::filter(!duplicated(!!id))

		df <- list(df_psm, df_str_b, df_str_a, df_num) %>%
				purrr::reduce(left_join, by = id) %>%
				data.frame(check.names = FALSE)

		df <- cbind.data.frame(df[, !grepl("[I|log2_R][0-9]{3}", names(df)), drop = FALSE],
					df[, grepl("[log2_R][0-9]{3}", names(df)), drop = FALSE],
					df[, grepl("[I][0-9]{3}", names(df)), drop = FALSE]) %>%
					dplyr::mutate_at(.vars = grep("[I|log2_R][0-9]{3}", names(.)),
					                 list(~replace(., is.infinite(.), NA) ))

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

		df <- na_zeroIntensity(df) %>% mutate(TMT_Set = set_idx)
	}


	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 2000000)

	id <- rlang::as_string(rlang::enexpr(id))

	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	filelist = list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
		reorder_files(n_TMT_sets(label_scheme_full))

	purrr::walk(as.list(filelist), function (x) {
		fn_prx <- gsub("_PSM_N.txt", "", x, fixed = TRUE)

		set_idx <- as.integer(gsub(".*TMTset(\\d+)_.*", "\\1", fn_prx))
		injn_idx <- as.integer(gsub(".*LCMSinj(\\d+).*", "\\1", fn_prx))

		df <- read.csv(file.path(dat_dir, "PSM", x), check.names = FALSE, header = TRUE,
		               sep = "\t", comment.char = "#") %>%
		  calcPepide(label_scheme = label_scheme, id = !!id, method_psm_pep = method_psm_pep,
		             set_idx = set_idx, injn_idx = injn_idx, range_int = range_int, range_log2r = range_log2r)

		if(grepl("|", df$prot_acc[1], fixed = TRUE)) {
			temp <- strsplit(as.character(df$prot_acc), '|', fixed = TRUE)
			temp <- plyr::ldply(temp, rbind)
			if(ncol(temp) > 1) df$prot_acc <- temp[, 2]
			rm(temp)
		}

		# df <- replace_na_genes(df, find_acctype(label_scheme))

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
#'See \code{\link{normPrn}} for the description of column keys.
#'
#'@param id The variable for summarisation of \code{PSMs} into peptides. The
#'  option \code{id = pep_seq} corresponds to the summarisation by the primary
#'  sequences of peptides. The option \code{id = pep_seq_mod} corresponds to the
#'  summarisation by both the primary sequences and variable modifications of
#'  peptides. \code{PSMs} data with the same value in \code{pep_seq} or
#'  \code{pep_seq_mod} will be summarised into a single entry of peptide.
#'@param method_psm_pep The method to summarise the \code{log2-ratios} and the
#'  \code{intensity} of \code{PSMs} by peptide entries. The descriptive
#'  statistics includes \code{c("mean", "median", "top.3", "weighted.mean")}.
#'
#'  The \code{log10-intensity} of reporter ions at the \code{PSMs} levels will
#'  be the weight when summarising \code{log2-ratios} with \code{"top.3"} or
#'  \code{"weighted.mean"}.
#'@param method_align The method to align the \code{log2-ratios} of peptides
#'  entries across samples. \code{MC}: median-centering; \code{MGKernel}: the
#'  kernal density defined by multiple Gaussian functions
#'  (\code{\link[mixtools]{normalmixEM}}).
#'
#'  At \code{method_align = "MC"}, the ratio profiles of each sample will be
#'  aligned in that the medians of the \code{log2-ratios} are zero. At
#'  \code{method_align = "MGKernel"}, the \code{log2-ratios} will be aligned in
#'  that the maximums of kernel density are zero. It is also possible to align
#'  the \code{log2-ratios} to the median of a list of user-supplied genes:
#'  \code{method_align = c("ACTB", "GAPDH", ...)}.
#'@param range_log2r The range of the \code{log2-ratios} of peptide entries for
#'  use in the scaling normalization of standard deviation across samples. The
#'  default is between the 20th and the 90th quantiles.
#'@param range_int The range of the \code{intensity} of reporter ions for use in
#'  the scaling normalization of standard deviation across samples. The default
#'  is between the 5th and the 95th quantiles.
#'@param n_comp The number of Gaussian components to be used with
#'  \code{method_align = "MGKernel"}. A typical value is 2 or 3. The variable
#'  \code{n_comp} overwrites the augument \code{k} in
#'  \code{\link[mixtools]{normalmixEM}}.
#'@param seed Integer; a seed for reproducible fitting at \code{method_align =
#'  MGKernel}.
#'@param annot_kinases Logical; if TRUE, annotates kinase attributes of
#'  proteins.
#'@param ... Additional parameters from \code{\link[mixtools]{normalmixEM}}:\cr
#'  \code{maxit}, the maximum number of iterations allowed; \cr \code{epsilon},
#'  tolerance limit for declaring algorithm convergence.
#'@inheritParams mixtools::normalmixEM
#'@seealso \code{\link{normPSM}} for PSM and \code{\link{normPrn}} for proteins.
#'
#'@return The primary output in \code{C:\\my_direcotry\\Peptide\\Peptide.txt}.
#'
#' @examples
#' 	normPep(
#' 		id = "pep_seq_mod",
#' 		method_psm_pep = "median",
#' 		method_align = "MGKernel",
#' 		range_log2r = c(10, 95),
#' 		range_int = c(5, 95),
#' 		n_comp = 3,
#' 		seed = 1234,
#' 		annot_kinases = TRUE,
#' 		maxit = 200,
#' 		epsilon = 1e-05
#' 	)
#'
#' \dontrun{
#' }
#'@import stringr dplyr tidyr purrr data.table rlang
#'@importFrom plyr ddply
#'@importFrom magrittr %>%
#'@export
normPep <- function (id = c("pep_seq", "pep_seq_mod"),
										method_psm_pep = c("median", "mean", "weighted.mean", "top.3"),
										method_align = c("MC", "MGKernel"), range_log2r = c(20, 90),
										range_int = c(5, 95), n_comp = NULL, seed = NULL,
										annot_kinases = FALSE, ...) {

	dir.create(file.path(dat_dir, "Peptide\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Peptide\\cache"), recursive = TRUE, showWarnings = FALSE)

	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 2000000)

	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	newColnames <- function(i, x) {
		label_scheme_sub <- label_scheme %>%
			dplyr::filter(TMT_Set == i)

		cols <- grep(paste0("[RI][0-9]{3}[NC]*_", i, "$"), names(x))
		nm_channel <- gsub(paste0("([RI][0-9]{3}[NC]*)_", i, "$"), "\\1", names(x)[cols])
		names(x)[cols] <- paste0(nm_channel, " (", as.character(label_scheme_sub$Sample_ID), ")")

		cols <- grep("[RI][0-9]{3}.*\\s+\\(.*\\)$", names(x))
		if (is.data.table(x)) { # for data.table
			if (length(cols) < ncol(x)) x <- cbind(x[, ..cols], x[, -..cols])
		} else { # for data.frame
			if (length(cols) < ncol(x)) x <- dplyr::bind_cols(x[, cols], x[, -cols, drop = FALSE])
		}

		return(x)
	}

	id <- rlang::enexpr(id)
	if(id == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
		id <- "pep_seq_mod"
	} else {
		id <- rlang::as_string(id)
		stopifnot(id %in% c("pep_seq", "pep_seq_mod"))
	}

	method_psm_pep <- rlang::enexpr(method_psm_pep)
	if(method_psm_pep == rlang::expr(c("median", "mean", "weighted.mean", "top.3"))) {
		method_psm_pep <- "median"
	} else {
		method_psm_pep <- rlang::as_string(method_psm_pep)
		stopifnot(method_psm_pep %in% c("mean", "top.3", "median", "weighted.mean"))
	}

	method_align <- rlang::enexpr(method_align)
	if(method_align == rlang::expr(c("MC", "MGKernel"))) {
		method_align <- "MC"
	} else {
		method_align <- rlang::as_string(method_align)
		if(!method_align %in% c("MC", "MGKernel"))
			warning("Assume the value of 'method_align' is a list of housekeeping proteins")
	}

	dots <- rlang::enexprs(...)

	mget(names(formals()), rlang::current_env()) %>% save_call("normPep")

	if(!file.exists(file.path(dat_dir, "Peptide", "Peptide.txt"))) {
		# normalize peptide data by each RAW file
		normPep_Splex(
			id = !!id,
			method_psm_pep = method_psm_pep
		)

		df <- do.call(rbind,
			lapply(list.files(path = file.path(dat_dir, "Peptide"), pattern = paste0("TMTset"),
			                  full.names = TRUE), read.csv, check.names = FALSE, header = TRUE,
			       sep = "\t", comment.char = "#")) %>%
			dplyr::mutate(TMT_Set = factor(TMT_Set)) %>%
			dplyr::arrange(TMT_Set)

		# remove peptides that have been assigned to more than one protein accession
		dup_peps <- df %>%
			dplyr::select(pep_seq, prot_acc) %>%
			dplyr::group_by(pep_seq) %>%
			dplyr::summarise(N = n_distinct(prot_acc)) %>%
			dplyr::filter(N > 1)

		if (nrow(dup_peps) > 0) {
			df <- df %>% dplyr::filter(!pep_seq %in% dup_peps$pep_seq)

			write.csv(df %>% dplyr::filter(pep_seq %in% dup_peps$pep_seq),
				file.path(dat_dir, "Peptide", "pep_w_ambiguity_in_prn_ids.csv"),
				sep = "\t", row.names = FALSE)
		}
		rm(dup_peps)

		# median summarisation of data from the same TMT experiment at different LCMS injections
		df_num <- df %>%
			dplyr::select(!!rlang::sym(id), TMT_Set, grep(
			  "^log2_R[0-9]{3}|^I[0-9]{3}|^N_log2_R[0-9]{3}|^N_I[0-9]{3}|^Z_log2_R[0-9]{3}", names(.))) %>%
			dplyr::group_by(!!rlang::sym(id), TMT_Set) %>%
			dplyr::summarise_all(~ median(., na.rm = TRUE))

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

		# faster when using data.table to transform "df_num" to wide format
		# df_num <- dcast(setDT(df_num), eval(as.name(Identifier)) ~ TMT_Set,
		#   value.var = names(df_num)[!names(df_num) %in% c(Identifier, "TMT_Set")]) %>%
		# 	dplyr::mutate(!! rlang::sym(id) := Identifier) %>%
		# 	dplyr::select(-Identifier) %>%
		# 	setDT(.)


		for (set_idx in seq_len(n_TMT_sets(label_scheme_full)))
		  df_num <- df_num %>% newColnames(set_idx, .)

		write.csv(df_num, file.path(dat_dir, "Peptide\\cache", "pep_num.csv"), row.names = FALSE)

		# calculate the number of PSM for each peptide
		df_psm <- df %>%
				dplyr::select(!!rlang::sym(id), n_psm) %>%
				dplyr::group_by(!!rlang::sym(id)) %>%
				dplyr::summarise(n_psm = sum(n_psm))

		# categories for summary using all strings
		max_char <- 500
		df_str_a <- df %>%
			dplyr::select(-n_psm, -TMT_Set,
			              -grep("^log2_R[0-9]{3}|^I[0-9]{3}|^N_log2_R[0-9]{3}|^N_I[0-9]{3}|^Z_log2_R[0-9]{3}",
			                    names(.))) %>% # already processed
			dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score",
			                                     "prot_matches", "prot_matches_sig", "prot_sequences",
			                                     "prot_sequences_sig"))) %>% # will not be used
			dplyr::select(-which(names(.) %in% c("pep_res_before", "pep_res_after", "pep_summed_mod_pos",
			                                     "pep_local_mod_pos"))) %>% # redundant but kept for convenience
			dplyr::select(-which(names(.) %in% c("pep_start", "pep_end", "pep_miss", "prot_acc", "prot_desc",
			                                     "prot_mass", "pep_var_mod"))) # will be processed later

		# will also be processed later
		if(id == "pep_seq_mod") df_str_a <- df_str_a %>% dplyr::select(-pep_seq)

		df_str_a <- df_str_a %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise_all(toString, na.rm = TRUE) %>%
			dplyr::mutate(pep_scan_title =
			                ifelse(nchar(pep_scan_title) > max_char,
			                       paste0(str_sub(pep_scan_title, 1, max_char), "..."), pep_scan_title))

		write.csv(df_str_a, file.path(dat_dir, "Peptide\\cache", "pep_stra.csv"), row.names = FALSE)

		# categories for summary using the first string
		nm_a <- names(df_str_a)[-which(names(df_str_a) == id)]
		df_str_b <- df %>%
			dplyr::select(-n_psm, -TMT_Set,
			              -grep("^log2_R[0-9]{3}|^I[0-9]{3}|^N_log2_R[0-9]{3}|^N_I[0-9]{3}|^Z_log2_R[0-9]{3}",
			                    names(.))) %>% # already processed
			dplyr::select(-which(names(.) %in% nm_a)) %>% # already processed under "df_str_a"
			dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score",
			                                     "prot_matches", "prot_matches_sig", "prot_sequences",
			                                     "prot_sequences_sig"))) %>% # will not be used
			dplyr::select(-which(names(.) %in% c("pep_res_before", "pep_res_after", "pep_summed_mod_pos",
			                                     "pep_local_mod_pos"))) %>% # redundant
			dplyr::filter(!duplicated(!!rlang::sym(id)))

		df <- list(df_str_b, df_psm, df_str_a, df_num) %>%
			purrr::reduce(left_join, by = id) %>%
			data.frame(check.names = FALSE)

		df[, grepl("I[0-9]{3}", names(df))] <- df[, grepl("I[0-9]{3}", names(df))] %>%
			dplyr::mutate_if(is.integer, as.numeric) %>%
			dplyr::mutate_if(is.logical, as.numeric) %>%
			round(., digits = 0)

		df[, grepl("log2_R[0-9]{3}", names(df))] <- df[, grepl("log2_R[0-9]{3}", names(df))] %>%
			dplyr::mutate_if(is.integer, as.numeric) %>%
			dplyr::mutate_if(is.logical, as.numeric) %>%
			round(., digits = 3)

		df <- cbind.data.frame(
			df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df))],
			df[, grep("^I[0-9]{3}", names(df))],
			df[, grep("^N_I[0-9]{3}", names(df))],
			df[, grep("^log2_R[0-9]{3}", names(df))],
			df[, grep("^N_log2_R[0-9]{3}", names(df))],
			df[, grep("^Z_log2_R[0-9]{3}", names(df))]
		)

		# annotatation including "refseq_acc" for kinase lookup
		acc_type <- find_acctype(label_scheme)
		df <- annotPrn(df, acc_type)
		if(annot_kinases) df <- annotKin(df, acc_type)

		# remove duplicated entries after the mapping
		df <- df %>%
			dplyr::filter(!duplicated(.[[id]]))

		df <- reorderCols(df, endColIndex = grep("I[0-9]{3}|R[0-9]{3}", names(df)), col_to_rn = id)

		df <- df[rowSums(!is.na(df[, grepl("N_log2_R", names(df))])) > 0, ]

		# df <- replace_na_genes(df, acc_type)

		write.csv(df, file.path(dat_dir, "Peptide\\cache", "Peptide_no_norm.csv"), row.names = FALSE)
	} else {
		df <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"),
			check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>%
			filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
	}

	if(is.null(n_comp)) n_comp <- if(nrow(df) > 3000) 3L else 2L

	df <- normMulGau(
		df = df,
		method_align = method_align,
		n_comp = n_comp,
		seed = seed,
		range_log2r = range_log2r,
		range_int = range_int,
		filepath = file.path(dat_dir, "Peptide\\Histogram"),
		!!!dots
	)

	df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-  df[, grepl("^Z_log2_R[0-9]{3}", names(df))] %>%
		dplyr::mutate_if(is.logical, as.numeric) %>%
		round(., digits = 3)

	write.table(df, file.path(dat_dir, "Peptide", "Peptide.txt"), sep = "\t",
	            col.names = TRUE, row.names = FALSE)

	invisible(df)
}
