#'Reports Protein Results
#'
#'\code{normPrn} summarises peptides into proteins and performs data
#'normalization.
#'
#'The column keys are defined in
#'\code{\href{http://www.matrixscience.com/help/csv_headers.html}{Matrix
#'Science}} with the following additions:
#'
#'\tabular{ll}{ \strong{Header}   \tab \strong{Descrption}\cr entrez \tab
#'\code{\href{https://www.ncbi.nlm.nih.gov/gene}{Protein Entrez ID}} \cr gene
#'\tab Protein gene name \cr kin_attr       \tab Optional; the attribute of
#'proteins being kinases \cr kin_class      \tab Optional; the classes of
#'kinases, e.g., TK, TKL... \cr kin_order      \tab Optional; the order of
#'"kin_class" using the \code{\href{http://kinase.com/human/kinome/}{kinase tree
#'diagram}} \cr n_pep          \tab The number of significant peptide sequences
#'matched to a proposed protein \cr n_psm \tab The number of significant MS/MS
#'queries matched to a peptide sequence \cr pep_seq_mod    \tab One-letter
#'representation of peptide sequences: amino acid residue(s) with variable
#'modifications are in the lower cases; the acetylations of protein N-terminals
#'indicated by '_' and the flanking residues on the N- or C-terminal side of
#'peptides separated by '.', e.g. "-._mAsGVAVSDGVIK.V"\cr prot_cover     \tab
#'Protein sequence coverage calculated from all available raw files   \cr
#'raw_file \tab MS filename(s) where peptides or proteins were identified \cr
#'refseq_acc \tab \code{\href{https://www.ncbi.nlm.nih.gov/refseq/}{Protein
#'RefSeq accession number}} \cr uniprot_acc    \tab
#'\code{\href{https://www.uniprot.org/help/accession_numbers}{Protein UniProt
#'accession number}} \cr uniprot_id     \tab
#'\code{\href{https://www.uniprot.org/help/entry_name}{Protein UniProt entry
#'name}} \cr I... (...)     \tab Reporter-ion intensity calculated from the
#'descriptive statistics in method_\code{psm_pep} for indicated samples \cr
#'N_I... (...)   \tab Normalized \code{I... (...)}. The calibration factors for
#'median centerring of \code{log2-ratios} are used to scale the reporter-ion
#'intensity \cr log2_R (...)   \tab \code{log2-ratios} relative to reference
#'materials for indicated samples \cr N_log2_R (...) \tab Aligned \code{log2_R
#'(...)} according to method_align without scaling normalization \cr Z_log2_R
#'(...) \tab \code{N_log2_R (...)} with scaling normalization \cr }
#'
#'@param id The variable to summarise peptides into proteins. The option
#'  \code{"prot_acc"} corresponds to the summarisation by the accession numbers
#'  or entry names of proteins The option \code{"gene"} corresponds to the
#'  summarisation by the gene names of proteins.
#'
#'  At \code{id = gene}, data under the same gene name but different acccesssion
#'  numbers or entry names will be summarised into one entry.
#'@param method_pep_prn The method to summarise the \code{log2-ratios} and the
#'  \code{intensity} of peptides by protein entries. The descriptive statistics
#'  includes \code{c("mean", "median", "top.3", "weighted.mean")}.
#'
#'  The representative \code{log10-intensity} of reporter ions at the peptide
#'  levels (from \code{\link{normPep}}) will be the weigth when summarising
#'  log2-ratios with "top.3" or "weighted.mean".
#'@param method_align The method to align the \code{log2-ratios} of proteins
#'  entries across samples. \code{MC}: median-centering; \code{MGKernel}: the
#'  kernal density defined by multiple Gaussian functions
#'  (\code{\link[mixtools]{normalmixEM}}).
#'
#'  At \code{method_align = "MC"}, the ratio profiles of each sample will be
#'  aligned in that the medians of the \code{log2-ratios} are zero. At
#'  \code{method_align = "MGKernel"}, the \code{log2-ratios} will be aligned in
#'  that the maximums of kernel density are zero. It is also possible to align
#'  the \code{log2-ratios} to the median of a list of user-supplied genes:
#'  \code{method_align = c("ACTB", "GAPDH", ...)}. At \code{method_align =
#'  NULL}, the protein summary by \code{method_pep_prn} will be taken without
#'  furthur normalization.
#'@param range_log2r The range of the \code{log2-ratios} of protein entries for
#'  use in the scaling normalization of standard deviation across samples. The
#'  default is between the 20th and the 90th quantiles.
#'@param range_int The range of the representative reporter-ion \code{intensity}
#'  of proteins for use in the scaling normalization of standard deviation
#'  across samples. The default is between the 5th and the 95th quantiles.
#'@param n_comp The number of Gaussian components to be used with
#'  \code{method_align = "MGKernel"}. A typical value is 2 or 3. The variable
#'  \code{n_comp} overwrites the augument \code{k} in
#'  \code{\link[mixtools]{normalmixEM}}.
#'@param seed Integer; a seed for reproducible fitting at \code{method_align =
#'  MGKernel}.
#'@param fasta The file name with prepended directory path to the fasta database
#'  being used in MS/MS ion search. If not specified, a pre-compiled system file
#'  will be used.
#'@param ... Additional parameters from \code{\link[mixtools]{normalmixEM}}:\cr
#'  \code{maxit}, the maximum number of iterations allowed; \cr \code{epsilon},
#'  tolerance limit for declaring algorithm convergence.
#'@inheritParams mixtools::normalmixEM
#'@family aggregate functions
#'@seealso \code{\link{normPSM}} for PSMs and \code{\link{normPep}} for
#'  peptides.
#'@return The primary output in "\code{C:\\my_directory\\Protein\\Protein.txt}".
#'
#' @examples
#' 	normPrn(
#' 		id = "prot_acc",
#' 		method_pep_prn = "median",
#' 		method_align = "MGKernel",
#' 		range_log2r = c(20, 90),
#' 		range_int = c(5, 95),
#' 		n_comp = 3,
#' 		maxit = 200,
#' 		epsilon = 1e-05
#' 	)
#'
#' 	normPrn(
#' 		id = "prot_acc",
#' 		fasta = "C:\\Results\\DB\\Refseq\\RefSeq_HM_Frozen_20130727.fasta"
#' 	)
#'
#' \dontrun{
#' }
#'@import stringr dplyr tidyr purrr data.table rlang mixtools
#'@importFrom magrittr %>%
#'@export
normPrn <- function (id = c("prot_acc", "gene"), 
                     method_pep_prn = c("median", "mean", "weighted.mean", "top.3"), 
                     method_align = c("MC", "MGKernel"), range_log2r = c(20, 90), 
                     range_int = c(5, 95), n_comp = NULL, seed = NULL, fasta = NULL, ...) {

	dir.create(file.path(dat_dir, "Protein\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Protein\\cache"), recursive = TRUE, showWarnings = FALSE)
	
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)
	
	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)
	
	# on.exit(message("Generation of a combined protein table --- Completed."), add = TRUE)

	options(max.print = 2000000) 
	
	method_pep_prn <- rlang::enexpr(method_pep_prn)
	if(method_pep_prn == rlang::expr(c("median", "mean", "weighted.mean", "top.3"))) {
		method_pep_prn <- "median"
	} else {
		method_pep_prn <- rlang::as_string(method_pep_prn)
		stopifnot(method_pep_prn %in% c("mean", "top.3", "median", "weighted.mean"))
	}
	
	method_align <- rlang::enexpr(method_align)
	if(method_align == rlang::expr(c("MC", "MGKernel"))) {
		method_align <- "MC"
	} else {
		method_align <- rlang::as_string(method_align)
		if(! method_align %in% c("MC", "MGKernel")) 
			warning("Assume the value of 'method_align' is a list of housekeeping proteins")
	}
	
	dots <- rlang::enexprs(...)
	
	id <- rlang::enexpr(id)
	if(id == rlang::expr(c("prot_acc", "gene"))) {
		id <- "gene"
	} else {
		id <- rlang::as_string(id)
	}
	
	mget(names(formals()), rlang::current_env()) %>% save_call("normPrn")
	stopifnot(id %in% c("prot_acc", "gene"))

	if(id == "gene") {
		gn_rollup <- TRUE
		id <- "prot_acc"
	} else {
		gn_rollup <- FALSE
	}
	
	stopifnot(id == "prot_acc")
	
	if(!file.exists(file.path(dat_dir, "Protein\\cache", "Protein_no_norm.csv"))) {
		df <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"), check.names = FALSE, 
										header = TRUE, sep = "\t", comment.char = "#") %>% 
			filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)

		# summarise data from the same TMT experiment at different LCMS injections
		df_num <- df %>% 
				dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
				dplyr::group_by(!!rlang::sym(id))

		df_num <- switch(method_pep_prn, 
					mean = aggrNums(mean)(df_num, !!rlang::sym(id), na.rm = TRUE), 
					top.3 = TMT_top_n(df_num, !!rlang::sym(id), na.rm = TRUE), 
					weighted.mean = TMT_wt_mean(df_num, !!rlang::sym(id), na.rm = TRUE), 
					median = aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE), 
					aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE)
				)
				
		write.csv(df_num, file.path(dat_dir, "Protein\\cache", "prn_num.csv"), row.names = FALSE)

		# the number of PSM
		df_psm <- df %>% 
				dplyr::select(!!rlang::sym(id), n_psm) %>% 
				dplyr::group_by(!!rlang::sym(id)) %>% 
				dplyr::summarise(n_psm = sum(n_psm))
		write.csv(df_psm, file.path(dat_dir, "Protein\\cache", "prn_npsm.csv"), row.names = FALSE)
		
		# the number of unique peptides
		df_seq <- df %>% 
				dplyr::select(!!rlang::sym(id), pep_seq) %>% 
				dplyr::filter(!duplicated(pep_seq)) %>% 
				dplyr::count(!!rlang::sym(id)) %>% # number of peptides for a protein
				dplyr::rename(n_pep = n)
		write.csv(df_seq, file.path(dat_dir, "Protein\\cache", "prn_npep.csv"), row.names = FALSE)
	
		# summarise categories using all strings
		max_char <- 200
		df_str_a <- df %>% 
				dplyr::select(-n_psm, 
				              -grep("^log2_R[0-9]{3}|^I[0-9]{3}|^N_log2_R[0-9]{3}|^N_I[0-9]{3}|^Z_log2_R[0-9]{3}", 
				                    names(.))) %>% # already processed 
				dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score", 
				                                     "prot_matches", "prot_matches_sig", "prot_sequences", 
				                                     "prot_sequences_sig"))) %>% # will not be used
				dplyr::select(-which(names(.) %in% c("pep_seq", "pep_var_mod", 
				                                     "pep_var_mod_pos"))) %>% # will not be used in protein report
				dplyr::select(-which(names(.) %in% c("prot_desc", "prot_mass", "uniprot_acc", "uniprot_id", 
				                                     "gene", "entrez", "kin_attr", "kin_class", 
				                                     "kin_order"))) %>% # will be processed later
				dplyr::group_by(!!rlang::sym(id)) %>% 
				dplyr::summarise_all(toString, na.rm = TRUE) %>% 
				dplyr::mutate_at(vars(which(names(.) != id)), 
				                 function (x) {ifelse(nchar(x) > max_char, 
				                                      paste0(str_sub(x, 1, max_char), "..."), x)})
		
		write.csv(df_str_a, file.path(dat_dir, "Protein\\cache", "prn_str_a.csv"), row.names = FALSE)
		
		# summarise categories using only the first string
		nm_a <- names(df_str_a)[-which(names(df_str_a) == id)]
		df_str_b <- df %>% 
				dplyr::select(-n_psm, 
				              -grep("^log2_R[0-9]{3}|^I[0-9]{3}|^N_log2_R[0-9]{3}|^N_I[0-9]{3}|^Z_log2_R[0-9]{3}", 
				                    names(.))) %>% # already processed 
				dplyr::select(-which(names(.) %in% nm_a)) %>% # already processed under "df_str_a"
				dplyr::select(-which(names(.) %in% c("prot_hit_num", "prot_family_member", "prot_score", 
				                                     "prot_matches", "prot_matches_sig", "prot_sequences", 
				                                     "prot_sequences_sig"))) %>% # will not be used
				dplyr::select(-which(names(.) %in% c("pep_seq", "pep_var_mod", "pep_var_mod_pos"))) %>% # redundant
				dplyr::filter(!duplicated(!!rlang::sym(id)))
		
		write.csv(df_str_b, file.path(dat_dir, "Protein\\cache", "prn_str_b.csv"), row.names = FALSE)
	
		df <- list(df_str_b, df_seq, df_psm, df_str_a, df_num) %>% 
				purrr::reduce(left_join, by = id) %>% 
				data.frame(check.names = FALSE)

		rm(df_num, df_psm, df_str_a, nm_a, df_str_b)

		df[, grepl("log2_R[0-9]{3}", names(df)) & !sapply(df, is.logical)] <- 
				df[, grepl("log2_R[0-9]{3}", names(df)) & !sapply(df, is.logical)] %>% 
				dplyr::mutate_if(is.integer, as.numeric) %>% 
				round(., digits = 3)

		df[, grepl("I[0-9]{3}", names(df)) & !sapply(df, is.logical)] <- 
				df[, grepl("I[0-9]{3}", names(df)) & !sapply(df, is.logical)] %>% 
				dplyr::mutate_if(is.integer, as.numeric) %>% 
				round(., digits = 0)
	
		df <- cbind.data.frame(
				df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df))], 
				df[, grep("^I[0-9]{3}", names(df))], 
				df[, grep("^N_I[0-9]{3}", names(df))], 
				df[, grep("^log2_R[0-9]{3}", names(df))], 
				df[, grep("^N_log2_R[0-9]{3}", names(df))], 
				df[, grep("^Z_log2_R[0-9]{3}", names(df))])

		acc_type <- load(file = file.path(dat_dir, "label_scheme.Rdata")) %>% find_acctype()
		
		# replace NA genes with prot_acc; some genes are represented by protein acessions
		df <- replace_na_genes(df, acc_type)
		df$gene <- as.factor(df$gene)

		df <- df[rowSums(!is.na(df[, grepl("N_log2_R", names(df))])) > 0, ] 
	
		# add prot_cover
		psm_data <- do.call(rbind, 
		                    lapply(list.files(path = file.path(dat_dir, "PSM"), 
		                                      pattern = "Clean\\.txt$", full.names = TRUE), 
		                           read.csv, 
		                           check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")) %>% 
							dplyr::select(prot_acc, pep_seq, pep_start, pep_end) %>% 
							dplyr::filter(!duplicated(pep_seq))	%>% 
							annotPrn(acc_type) 
		
		write.csv(psm_data, file.path(dat_dir, "Protein\\cache", "Protein_coverage.csv"), row.names = FALSE)
		
		df <- psm_data %>% 
			calc_cover(id = !!id, fasta = fasta) %>% 
			dplyr::right_join(df, by = id)
		
		write.csv(df, file.path(dat_dir, "Protein\\cache", "Protein_no_norm.csv"), row.names = FALSE)
	} else {
		df <- read.csv(file.path(dat_dir, "Protein\\cache", "Protein_no_norm.csv"), 
										check.names = FALSE, header = TRUE, comment.char = "#") %>% 
			dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
	}
	
	if(is.null(n_comp)) n_comp <- if(nrow(df) > 3000) 3L else 2L
	
	df <- normMulGau(
		df = df, 
		method_align = method_align, 
		n_comp = n_comp, 
		seed = seed, 
		range_log2r = range_log2r, 
		range_int = range_int, 
		filepath = file.path(dat_dir, "Protein\\Histogram"), 
		!!!dots
	)

	df[, grepl("log2_R[0-9]{3}", names(df))] <- df[, grepl("log2_R[0-9]{3}", names(df))] %>% 
		dplyr::mutate_if(is.integer, as.numeric) %>% 
		dplyr::mutate_if(is.logical, as.numeric) %>% 
		round(., digits = 3)
	
	if(gn_rollup) {
		dfa <- df %>% 
			dplyr::select(gene, grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
			dplyr::filter(!is.na(gene)) %>% # proteins with empty "gene" but "prot_acc" are still kept 
			dplyr::group_by(gene) %>% 
			dplyr::summarise_all(list(~median(., na.rm = TRUE)))
		
		dfb <- df %>% 
			dplyr::select(-prot_cover, -grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
			dplyr::filter(!is.na(gene)) 
		
		dfc <- df %>% 
			dplyr::select(gene, prot_cover) %>% 
			dplyr::filter(!is.na(gene), !is.na(prot_cover)) %>% 
			dplyr::group_by(gene) %>% 
			dplyr::mutate(prot_cover = as.numeric(sub("%", "", prot_cover))) %>% 
			dplyr::summarise_all(~max(., na.rm = TRUE)) %>% 
      dplyr::mutate(prot_cover = paste0(prot_cover, "%"))

		# the number of unique gene in psm_data may be shorter than those in 
		# dfa and dfb for the reason of empty "genes" in dfa 
		df <- list(dfc, dfb, dfa) %>% 
			purrr::reduce(right_join, by = "gene") %>% 
			dplyr::filter(!is.na(gene), !duplicated(gene))
	}
	
	df <- df %>% dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0)
	
	write.table(df, file.path(dat_dir, "Protein", "Protein.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
	
	invisible(df)
}

