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
#'raw_file \tab MS file name(s) where peptides or proteins were identified \cr
#'refseq_acc \tab \code{\href{https://www.ncbi.nlm.nih.gov/refseq/}{Protein
#'RefSeq accession number}} \cr uniprot_acc    \tab
#'\code{\href{https://www.uniprot.org/help/accession_numbers}{Protein UniProt
#'accession number}} \cr uniprot_id     \tab
#'\code{\href{https://www.uniprot.org/help/entry_name}{Protein UniProt entry
#'name}} \cr I... (...)     \tab Reporter-ion intensity calculated from the
#'descriptive statistics in method_\code{psm_pep} for indicated samples \cr
#'N_I... (...)   \tab Normalized \code{I... (...)}. The calibration factors for
#'median centerring of \code{log2FC} are used to scale the reporter-ion
#'intensity \cr log2_R (...)   \tab \code{log2FC} relative to reference
#'materials for indicated samples \cr N_log2_R (...) \tab Aligned \code{log2_R
#'(...)} according to method_align without scaling normalization \cr Z_log2_R
#'(...) \tab \code{N_log2_R (...)} with scaling normalization \cr }
#'
#'@param id The variable to summarise peptides into proteins. The option
#'  \code{"prot_acc"} corresponds to the summarisation by the accession numbers
#'  or entry names of proteins The option \code{"gene"} corresponds to the
#'  summarisation by the gene names of proteins. At \code{id = gene}, data under
#'  the same gene name but different acccesssion numbers or entry names will be
#'  summarised into one entry.
#'@param method_pep_prn The method to summarise the \code{log2FC} and the
#'  \code{intensity} of peptides by protein entries. The descriptive statistics
#'  includes \code{c("mean", "median", "top.3", "weighted.mean")}. The
#'  representative \code{log10-intensity} of reporter ions at the peptide levels
#'  (from \code{\link{normPep}}) will be the weigth when summarising \code{log2FC}
#'  with "top.3" or "weighted.mean".
#'@inheritParams normPep
#'@param fasta The file name with prepended directory path to the fasta database
#'  being used in MS/MS ion search. If not specified, a pre-compiled system file
#'  will be used.
#'@inheritParams mixtools::normalmixEM
#'@family aggregate functions
#'@seealso \code{\link{normPSM}} for PSMs and \code{\link{normPep}} for
#'  peptides.
#'@return The primary output in "\code{~\\dat_dir\\Protein\\Protein.txt}".
#'
#' @examples
#' 	normPrn(
#' 		id = gene,
#' 		method_pep_prn = median,
#' 		method_align = MGKernel,
#' 		range_log2r = c(5, 95),
#' 		range_int = c(5, 95),
#' 		n_comp = 2,
#' 		seed = 749662,
#' 		fasta = "C:\\Results\\DB\\Refseq\\RefSeq_HM_Frozen_20130727.fasta",
#' 		maxit = 200,
#' 		epsilon = 1e-05
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
                     range_int = c(5, 95), n_comp = NULL, seed = NULL, fasta = NULL, 
                     col_refit = NULL, ...) {

	dir.create(file.path(dat_dir, "Protein\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Protein\\cache"), recursive = TRUE, showWarnings = FALSE)
	
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)
	
	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)
	
	# on.exit(message("Generation of a combined protein table --- Completed."), add = TRUE)

	options(max.print = 2000000) 
	
	reload_expts()
	
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
	
	id <- rlang::enexpr(id)
	if(id == rlang::expr(c("prot_acc", "gene"))) {
		id <- "gene"
	} else {
		id <- rlang::as_string(id)
	}
	
	col_refit <- rlang::enexpr(col_refit)
	col_refit <- ifelse(is.null(col_refit), rlang::expr(Sample_ID), rlang::sym(col_refit))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	
	if(is.null(label_scheme[[col_refit]])) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("Column \'", rlang::as_string(col_refit), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_refit]])) == 0) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("No samples were specified under column \'", rlang::as_string(col_refit), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
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
	
	if(!is.null(fasta)) {
	  fasta <- rlang::as_string(rlang::enexpr(fasta))
	  stopifnot(file.exists(fasta))
	}
	
	if(!file.exists(file.path(dat_dir, "Protein", "Protein.txt"))) {
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
		                 aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE))
				
		write.csv(df_num, file.path(dat_dir, "Protein\\cache", "prn_num.csv"), row.names = FALSE)

		df_psm <- df %>% 
				dplyr::select(!!rlang::sym(id), n_psm) %>% 
				dplyr::group_by(!!rlang::sym(id)) %>% 
				dplyr::summarise(n_psm = sum(n_psm))
		write.csv(df_psm, file.path(dat_dir, "Protein\\cache", "prn_npsm.csv"), row.names = FALSE)
		
		pep_id <- match_identifier("pep_seq")
		
		df_seq <- df %>% 
				dplyr::select(!!rlang::sym(id), !!rlang::sym(pep_id)) %>% 
				dplyr::filter(!duplicated(!!rlang::sym(pep_id))) %>% 
				dplyr::count(!!rlang::sym(id)) %>% # number of peptides for a protein
				dplyr::rename(n_pep = n)
		write.csv(df_seq, file.path(dat_dir, "Protein\\cache", "prn_npep.csv"), row.names = FALSE)
	
		df_first <- df %>% 
		  dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
		  dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
		  dplyr::select(-n_psm)	%>% 
		  dplyr::select(-pep_start, -pep_end, -pep_score, -pep_expect) %>% 
		  dplyr::select(-!!rlang::sym(pep_id))

		df <- list(df_seq, df_psm, df_first, df_num) %>% 
				purrr::reduce(left_join, by = id) %>% 
				data.frame(check.names = FALSE)

		rm(df_num, df_psm, df_seq, df_first)

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
		df <- do.call(rbind, 
		                    lapply(list.files(path = file.path(dat_dir, "PSM"), 
		                                      pattern = "Clean\\.txt$", full.names = TRUE), 
		                           read.csv, 
		                           check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")) %>% 
							dplyr::select(prot_acc, prot_desc, pep_seq, pep_start, pep_end) %>% 
							dplyr::filter(!duplicated(pep_seq))	%>% 
							annotPrn(acc_type) %>% 
		  calc_cover(id = !!id, fasta = fasta) %>% 
		  dplyr::right_join(df, by = id)
		
		write.csv(df, file.path(dat_dir, "Protein\\cache", "Protein_no_norm.csv"), row.names = FALSE)
		
		if (gn_rollup) {
		  # a slow step
		  dfa <- df %>% 
		    dplyr::select(gene, grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
		    dplyr::filter(!is.na(gene)) %>% # proteins with empty "gene" but "prot_acc" will be kept 
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
		  write.csv(df, file.path(dat_dir, "Protein\\cache", "prn_rollup.csv"), row.names = FALSE)
		}
	} else {
	  df <- read.csv(file.path(dat_dir, "Protein", "Protein.txt"), sep = "\t", 
	                 check.names = FALSE, header = TRUE, comment.char = "#") %>% 
	    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
	}
	
	df <- normMulGau(
		df = df, 
		method_align = method_align, 
		n_comp = n_comp, 
		seed = seed, 
		range_log2r = range_log2r, 
		range_int = range_int, 
		filepath = file.path(dat_dir, "Protein\\Histogram"), 
		col_refit = col_refit, 
		!!!dots
	)

	df <- df %>% dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0)
	
	df <- df %>% 
	  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
	  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
	  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
	  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3))
	
	
	write.table(df, file.path(dat_dir, "Protein", "Protein.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
	
	invisible(df)
}
