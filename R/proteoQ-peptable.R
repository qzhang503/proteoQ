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
#' @import stringr dplyr purrr rlang
#' @importFrom magrittr %>%
normPep_Splex <- function (id = "pep_seq_mod", method_psm_pep = "median") {

	dir.create(file.path(dat_dir, "Peptide\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	on.exit(message("Generation of individual peptide tables by RAW filenames --- Completed."),
	        add = TRUE)

	calcPepide <- function(df, label_scheme, id, method_psm_pep, set_idx, injn_idx) {
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
		
		df_sd <- df %>% 
		  dplyr::select(!!rlang::sym(id), grep("^log2_R", names(.))) %>%
		  dplyr::group_by(!!rlang::sym(id)) %>%
		  dplyr::summarise_at(vars(starts_with("log2_R")), ~ sd(.x, na.rm = TRUE)) %>% 
		  dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
		  `names_pos<-`(2:ncol(.), paste0("sd_", names(.)[2:ncol(.)]))
		
		dir.create(file.path(dat_dir, "Peptide\\SD"), recursive = TRUE, showWarnings = FALSE)
		
		# violin plots of SD
		if (TMT_plex(label_scheme) > 0) {
		  TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()
		  Levels <- TMT_levels %>% gsub("^TMT-", "", .)
		  
		  df_sd2 <- df_sd %>%
		    `names<-`(gsub("sd_log2_R", "", names(.))) %>% 
		    tidyr::gather(key = !!rlang::sym(id), value = "SD") %>%
		    dplyr::rename(Channel := !!rlang::sym(id)) %>% 
		    dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
		    dplyr::filter(!is.na(SD))
		  
		  Width <- 8
		  Height <- 8
		  filename <- paste0("SD_set_", set_idx, "_inj_", injn_idx)
		  
		  p <- ggplot() +
		    geom_violin(df_sd2, mapping = aes(x = Channel, y = SD, fill = Channel), size = .25) +
		    geom_boxplot(df_sd2, mapping = aes(x = Channel, y = SD), width = 0.1, lwd = .2, fill = "white") +
		    stat_summary(df_sd2, mapping = aes(x = Channel, y = SD), fun.y = "mean", geom = "point",
		                 shape=23, size=2, fill="white", alpha=.5) +
		    labs(title = expression("Peptide"), x = expression("Channel"), y = expression("SD ("*log[2]*"FC)")) +
		    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
		    theme_psm_violin
		  
		  ggsave(file.path(dat_dir, "Peptide\\SD", paste0(filename, ".png")), p, width = Width, height = Height, units = "in")
		  rm(TMT_levels, Levels, Width, Height, filename, p, df_sd2)
		}
		
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

		df_psm <- df %>%
				dplyr::select(!!rlang::sym(id)) %>%
				dplyr::group_by(!!rlang::sym(id)) %>%
				dplyr::summarise(n_psm = n())
		
		df_score <- df %>% 
		  dplyr::select(!!rlang::sym(id), pep_score) %>% 
		  dplyr::group_by(!!rlang::sym(id)) %>% 
		  dplyr::summarise(pep_score = max(pep_score, na.rm = TRUE))
		
		df_expect <- df %>% 
		  dplyr::select(!!rlang::sym(id), pep_expect) %>% 
		  dplyr::group_by(!!rlang::sym(id)) %>% 
		  dplyr::summarise(pep_expect = min(pep_expect, na.rm = TRUE))

		df <- df %>% 
		  dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% # already in `df_num`
		  dplyr::select(-pep_score, -pep_expect) %>% # in `df_score` and `df_expect`
		  dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score",
		                                       "prot_matches", "prot_matches_sig", "prot_sequences",
		                                       "prot_sequences_sig", "raw_file", "pep_query", 
		                                       "pep_rank", "pep_isbold", "pep_isunique", 
		                                       "pep_exp_mr", "pep_exp_z", "pep_calc_mr", 
		                                       "pep_exp_mz", "pep_delta", "pep_miss", "pep_var_mod", 
		                                       "pep_var_mod_pos", "pep_scan_title", "pep_res_before", 
		                                       "pep_res_after", "pep_summed_mod_pos", 
		                                       "pep_local_mod_pos")))
		
		df_first <- df %>% 
		  dplyr::filter(!duplicated(!!rlang::sym(id)))
		
		if(id == "pep_seq_mod") {
		  df_first <- df_first %>% dplyr::select(-pep_seq)
		} else {
		  df_first <- df_first %>% dplyr::select(-pep_seq_mod)
		}

		df <- list(df_psm, df_sd, df_first, df_score, df_expect, df_num) %>%
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
	method_psm_pep <- rlang::as_string(rlang::enexpr(method_psm_pep))

	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
		reorder_files(n_TMT_sets(label_scheme_full))

	purrr::walk(as.list(filelist), function (x) {
		fn_prx <- gsub("_PSM_N.txt", "", x, fixed = TRUE)

		set_idx <- as.integer(gsub(".*TMTset(\\d+)_.*", "\\1", fn_prx))
		injn_idx <- as.integer(gsub(".*LCMSinj(\\d+).*", "\\1", fn_prx))

		df <- read.csv(file.path(dat_dir, "PSM", x), check.names = FALSE, header = TRUE,
		               sep = "\t", comment.char = "#") 
		
		df <- df %>%
		  calcPepide(label_scheme = label_scheme, id = !!id, method_psm_pep = method_psm_pep,
		             set_idx = set_idx, injn_idx = injn_idx)

		if(grepl("|", df$prot_acc[1], fixed = TRUE)) {
			temp <- strsplit(as.character(df$prot_acc), '|', fixed = TRUE)
			temp <- plyr::ldply(temp, rbind)
			if(ncol(temp) > 1) df$prot_acc <- temp[, 2]
			rm(temp)
		}

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
#'See \code{\link{normPrn}} for the description of column keys in the output.
#'
#'@param id Character string; the variable for summarisation of \code{PSMs} into
#'  peptides. The option \code{id = pep_seq} corresponds to the summarisation by
#'  the primary sequences of peptides. The option \code{id = pep_seq_mod}
#'  corresponds to the summarisation by both the primary sequences and variable
#'  modifications of peptides. \code{PSMs} data with the same value in
#'  \code{pep_seq} or \code{pep_seq_mod} will be summarised into a single entry
#'  of peptide.
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
#'  the 20th and the 90th quantiles.
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
#'@param annot_kinases Logical; if TRUE, annotates kinase attributes of
#'  proteins.
#'@param col_refit Character string to a column key in \code{expt_smry.xlsx}.
#'  Samples corresponding to non-empty entries under \code{col_refit} will be
#'  used in the refit of \code{log2FC} using multiple Gaussian kernels. The
#'  density estimates from an earlier analyis will be kept for the remaining
#'  samples. At the NULL default, the column key \code{Sample_ID} will be used,
#'  which results in the refit of the \code{log2FC} for all samples.
#'@param cache Logical; if TRUE, use cache.
#'@param ... Additional parameters for \code{\link[mixtools]{normalmixEM}}:\cr
#'  \code{maxit}, the maximum number of iterations allowed; \cr \code{epsilon},
#'  tolerance limit for declaring algorithm convergence.
#'@inheritParams mixtools::normalmixEM
#'@seealso \code{\link{normPSM}} for PSM an \code{\link{normPrn}} for proteins.
#'
#'@return The primary output in \code{~\\dat_dir\\Peptide\\Peptide.txt}.
#'
#' @examples
#' \dontrun{
#'normPep(
#'  id = pep_seq_mod,
#'  method_psm_pep = median,
#'  method_align = MGKernel,
#'  range_log2r = c(10, 95),
#'  range_int = c(5, 95),
#'  n_comp = 3,
#'  seed = 1234,
#'  maxit = 200,
#'  epsilon = 1e-05
#')
#' }
#'@import stringr dplyr tidyr purrr data.table rlang
#'@importFrom plyr ddply
#'@importFrom magrittr %>%
#'@export
normPep <- function (id = c("pep_seq", "pep_seq_mod"),
										method_psm_pep = c("median", "mean", "weighted.mean", "top.3"),
										method_align = c("MC", "MGKernel"), range_log2r = c(20, 90),
										range_int = c(5, 95), n_comp = NULL, seed = NULL,
										annot_kinases = FALSE, col_refit = NULL, cache = TRUE, ...) {

	dir.create(file.path(dat_dir, "Peptide\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Peptide\\cache"), recursive = TRUE, showWarnings = FALSE)

	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 2000000)
	
	reload_expts()

	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	newColnames <- function(i, x) {
		label_scheme_sub <- label_scheme %>%
			dplyr::filter(TMT_Set == i)

		cols <- grep(paste0("[RI][0-9]{3}[NC]*_", i, "$"), names(x))
		nm_channel <- gsub(paste0("([RI][0-9]{3}[NC]*)_", i, "$"), "\\1", names(x)[cols])
		names(x)[cols] <- paste0(nm_channel, " (", as.character(label_scheme_sub$Sample_ID), ")")

		cols <- grep("[RI][0-9]{3}.*\\s+\\(.*\\)$", names(x))
		if (is.data.table(x)) {
			if (length(cols) < ncol(x)) x <- cbind(x[, ..cols], x[, -..cols])
		} else {
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

	col_refit <- rlang::enexpr(col_refit)
	col_refit <- ifelse(is.null(col_refit), rlang::expr(Sample_ID), rlang::sym(col_refit))
	
	if(is.null(label_scheme[[col_refit]])) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("Column \'", rlang::as_string(col_refit), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_refit]])) == 0) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("No samples were specified under column \'", rlang::as_string(col_refit), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
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
	
	stopifnot(length(range_log2r) == 2)
	stopifnot(length(range_int) == 2)
	
	if (range_log2r[2] <= 1) range_log2r <- range_log2r * 100
	if (range_int[2] <= 1) range_int <- range_int * 100
	
	stopifnot(range_log2r[1] < range_log2r[2] & 
	            range_log2r[1] >= 0 & range_log2r[2] <= 100)
	
	stopifnot(range_int[1] < range_int[2] & 
	            range_int[1] >= 0 & range_int[2] <= 100)
	
	if(is.null(n_comp)) n_comp <- if(nrow(df) > 3000) 3L else 2L
	n_comp <- n_comp %>% as.integer()
	stopifnot(n_comp >= 2)

	dots <- rlang::enexprs(...)

	mget(names(formals()), rlang::current_env()) %>% save_call("normPep")

	if(!(cache & file.exists(file.path(dat_dir, "Peptide", "Peptide.txt")))) {
		# normalize peptide data within each RAW file
	  normPep_Splex(id = !!id, method_psm_pep = method_psm_pep)
	  
		df <- do.call(rbind,
			lapply(list.files(path = file.path(dat_dir, "Peptide"), pattern = paste0("TMTset"),
			                  full.names = TRUE), read.csv, check.names = FALSE, header = TRUE,
			       sep = "\t", comment.char = "#")) %>%
			dplyr::mutate(TMT_Set = factor(TMT_Set)) %>%
			dplyr::arrange(TMT_Set)

		# remove peptides that have been assigned to more than one protein accession
		dup_peps <- df %>%
			dplyr::select(!!rlang::sym(id), prot_acc) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise(N = n_distinct(prot_acc)) %>%
			dplyr::filter(N > 1)

		if (nrow(dup_peps) > 0) {
			df <- df %>% dplyr::filter(! (!!rlang::sym(id) %in% dup_peps[[id]]))

			write.csv(dup_peps,
				file.path(dat_dir, "Peptide", "dbl_dipping_peptides.csv"),
				row.names = FALSE)
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

		df_first <- df %>% 
		  dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
		  dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
		  dplyr::select(-n_psm, -TMT_Set)

		df <- list(df_psm, df_first, df_num) %>%
			purrr::reduce(left_join, by = id) %>%
			data.frame(check.names = FALSE)
		
		df <- df %>% 
		  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
		  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
		  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
		  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3))
		
		df <- cbind.data.frame(
			df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df))],
			df[, grep("^I[0-9]{3}", names(df))],
			df[, grep("^N_I[0-9]{3}", names(df))],
			df[, grep("^log2_R[0-9]{3}", names(df))],
			df[, grep("^N_log2_R[0-9]{3}", names(df))],
			df[, grep("^Z_log2_R[0-9]{3}", names(df))]
		)

		# annotation including historical lookup of kinases by "refseq_acc"
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

	df <- normMulGau(
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

	invisible(df)
}
