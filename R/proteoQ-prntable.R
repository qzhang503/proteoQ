#'Reports Protein Results
#'
#'\code{normPrn} summarises peptides into proteins and performs data
#'normalization.
#'
#'For the analysis of Mascot results, the column keys are defined in
#'\code{\href{http://www.matrixscience.com/help/csv_headers.html}{Matrix
#'Science}} with the following additions:
#'
#'\tabular{ll}{ \strong{Header}   \tab \strong{Descrption}\cr entrez \tab
#'\code{\href{https://www.ncbi.nlm.nih.gov/gene}{Protein Entrez ID}} \cr gene
#'\tab Protein gene name \cr kin_attr       \tab Optional; the attribute of
#'proteins being kinases \cr kin_class      \tab Optional; the classes of
#'kinases, e.g., TK, TKL... \cr kin_order      \tab Optional; the order of
#'"kin_class" from the \code{\href{http://kinase.com/human/kinome/}{kinase tree
#'diagram}} \cr prot_n_pep \tab The number of significant peptide sequences
#'matched to a proposed protein. \cr prot_n_psm \tab The number of significant
#'MS/MS queries matched to a proposed protein. \cr pep_n_psm \tab The number of
#'significant MS/MS queries matched to a proposed peptide sequence. \cr
#'pep_seq_mod    \tab One-letter representation of peptide sequences: amino acid
#'residue(s) with variable modifications are in the lower cases; the
#'acetylations of protein N-terminals indicated by '_' and the flanking residues
#'on the N- or C-terminal side of peptides separated by '.', e.g.
#'"-._mAsGVAVSDGVIK.V"\cr prot_cover     \tab Protein sequence coverage
#'calculated from all available raw files   \cr raw_file \tab MS file name(s)
#'where peptides or proteins were identified \cr refseq_acc \tab
#'\code{\href{https://www.ncbi.nlm.nih.gov/refseq/}{Protein RefSeq accession
#'number}} \cr uniprot_acc    \tab
#'\code{\href{https://www.uniprot.org/help/accession_numbers}{Protein UniProt
#'accession number}} \cr uniprot_id     \tab
#'\code{\href{https://www.uniprot.org/help/entry_name}{Protein UniProt entry
#'name}} \cr I... (...)     \tab Reporter-ion intensity calculated from the
#'descriptive statistics in \code{mothod_psm_pep} or \code{mothod_pep_prn} for
#'indicated samples \cr N_I... (...)   \tab Normalized \code{I... (...)}. The
#'calibration factors for the alignment of \code{log2FC} are used to scale the
#'reporter-ion intensity \cr log2_R (...)   \tab \code{log2FC} relative to
#'reference materials for indicated samples \cr N_log2_R (...) \tab Aligned
#'\code{log2_R (...)} according to \code{method_align} without scaling
#'normalization \cr Z_log2_R (...) \tab \code{N_log2_R (...)} with scaling
#'normalization \cr }
#'
#'@param id Character string; the variable to summarise peptides into proteins.
#'  The option \code{prot_acc} corresponds to the summarisation by the accession
#'  numbers or entry names of proteins. The option \code{gene} corresponds to
#'  the summarisation by the gene names of proteins. At \code{id = gene}, data
#'  under the same gene name but different acccesssion numbers or entry names
#'  will be summarised into one entry. Currently, the value of \code{id} will
#'  match automatically to the value of \code{group_pep_by} in \code{normPep}.
#'@param method_pep_prn Character string; the method to summarise the
#'  \code{log2FC} and the \code{intensity} of peptides by protein entries. The
#'  descriptive statistics includes \code{c("mean", "median", "top.3",
#'  "weighted.mean")}. The representative \code{log10-intensity} of reporter
#'  ions at the peptide levels (from \code{\link{normPep}}) will be the weigth
#'  when summarising \code{log2FC} with \code{top.3} or \code{weighted.mean}.
#'@param use_unique_pep Logical; if TRUE, only peptides that are unique by
#'  protein groups or families will be used to summarise the \code{log2FC} and
#'  the \code{intensity} of peptides by protein entries.
#'@inheritParams normPep
#'@inheritParams mixtools::normalmixEM
#'@family aggregate functions
#'@seealso \code{\link{normPSM}} for PSMs and \code{\link{normPep}} for
#'  peptides.
#'@return The primary output is in "\code{~\\dat_dir\\Protein\\Protein.txt}".
#'
#' @examples
#' \dontrun{
#' # proteins results with examplary `filter_...`
#' normPrn(
#'   method_pep_prn = median,
#'   method_align = MGKernel,
#'   range_log2r = c(20, 95),
#'   range_int = c(5, 95),
#'   n_comp = 2,
#'   seed = 749662,
#'   maxit = 200,
#'   epsilon = 1e-05,
#'
#'   filter_by = exprs(prot_n_psm >= 5, prot_n_pep >=2),
#' )
#'
#' }
#'@import stringr dplyr tidyr purrr data.table rlang mixtools
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@export
normPrn <- function (id = c("prot_acc", "gene"), 
                     method_pep_prn = c("median", "mean", "weighted.mean", "top.3"), 
                     use_unique_pep = TRUE, method_align = c("MC", "MGKernel"), 
                     range_log2r = c(10, 90), range_int = c(5, 95), n_comp = NULL, seed = NULL, 
                     col_refit = NULL, cache = TRUE, ...) {

	dir.create(file.path(dat_dir, "Protein\\Histogram"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Protein\\cache"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(dat_dir, "Protein\\log"), recursive = TRUE, showWarnings = FALSE)
	
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)
	
	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 2000000) 
	
	reload_expts()
	
	method_pep_prn <- rlang::enexpr(method_pep_prn)
	if (method_pep_prn == rlang::expr(c("median", "mean", "weighted.mean", "top.3"))) {
		method_pep_prn <- "median"
	} else {
		method_pep_prn <- rlang::as_string(method_pep_prn)
		stopifnot(method_pep_prn %in% c("mean", "top.3", "median", "weighted.mean"))
	}
	
	method_align <- rlang::enexpr(method_align)
	if (method_align == rlang::expr(c("MC", "MGKernel"))) {
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
	
	if (is.null(n_comp)) n_comp <- if(nrow(df) > 3000) 3L else 2L
	n_comp <- n_comp %>% as.integer()
	stopifnot(n_comp >= 2)
	
	dots <- rlang::enexprs(...)
	lang_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	dots <- dots %>% .[! . %in% lang_dots]
	
	# depreciated
	id <- rlang::enexpr(id)
	if (id == rlang::expr(c("prot_acc", "gene"))) {
		id <- "gene"
	} else {
		id <- rlang::as_string(id)
	}
	
	id <- match_normPSM_protid()
	pep_id <- match_normPSM_pepid()

	col_refit <- rlang::enexpr(col_refit)
	col_refit <- ifelse(is.null(col_refit), rlang::expr(Sample_ID), rlang::sym(col_refit))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	
	if (is.null(label_scheme[[col_refit]])) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("Column \'", rlang::as_string(col_refit), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_refit]])) == 0) {
	  col_refit <- rlang::expr(Sample_ID)
	  warning("No samples were specified under column \'", rlang::as_string(col_refit), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
	}
	
	mget(names(formals()), rlang::current_env()) %>% save_call("normPrn")
	stopifnot(id %in% c("prot_acc", "gene"))

	if (id == "gene") {
		gn_rollup <- TRUE
		id <- "prot_acc"
	} else {
		gn_rollup <- FALSE
	}
	
	stopifnot(id == "prot_acc")
	
	if (!(cache & file.exists(file.path(dat_dir, "Protein", "Protein.txt")))) {
	  fasta <- seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
	                              seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
	  
	  df <- read.csv(file.path(dat_dir, "Peptide\\Peptide.txt"), check.names = FALSE, 
		               header = TRUE, sep = "\t", comment.char = "#") %>% 
			dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
	  
	  cat("Available column keys for data filtration: \n")
	  cat(paste0(names(df), "\n"))
	  df <- df %>% filters_in_call(!!!lang_dots)
	  
	  if (all(c("pep_start", "pep_end") %in% names(df))) {
	    if ("prot_cover" %in% names(df)) {
	      warning("`prot_cover` recalculated after data merging.\n")
	      df$prot_cover <- NULL
	    } 
	    
	    df <- df %>% 
	      dplyr::select(prot_acc, acc_type, prot_desc, !!rlang::sym(pep_id), pep_start, pep_end) %>% 
	      calc_cover(id = !!id, fasta = fasta) %>% 
	      dplyr::right_join(df, by = id)		  
	  } else {
	    df$prot_cover <- NA
	  }
		
		if (use_unique_pep & "pep_isunique" %in% names(df)) 
		  df <- df %>% dplyr::filter(pep_isunique == 1)

		df_num <- df %>% 
				dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
				dplyr::group_by(!!rlang::sym(id))

		df_num <- switch(method_pep_prn, 
		                 mean = aggrNums(mean)(df_num, !!rlang::sym(id), na.rm = TRUE), 
		                 top.3 = TMT_top_n(df_num, !!rlang::sym(id), na.rm = TRUE), 
		                 weighted.mean = TMT_wt_mean(df_num, !!rlang::sym(id), na.rm = TRUE), 
		                 median = aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE), 
		                 aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE))

		df_first <- df %>% 
		  dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
		  dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
		  dplyr::select(-grep("^pep_", names(.)))
		
		df <- list(df_first, df_num) %>% 
		  purrr::reduce(left_join, by = id) %>% 
		  data.frame(check.names = FALSE)

		rm(df_num, df_first)

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
		
		df <- df %>% 
		  .[rowSums(!is.na(.[, grepl("N_log2_R", names(.))])) > 0, ]

		if (gn_rollup) {
		  dfa <- df %>% 
		    dplyr::select(gene, grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
		    dplyr::filter(!is.na(gene)) %>% # proteins with empty "gene" but "prot_acc" will be kept 
		    dplyr::group_by(gene) %>% 
		    dplyr::summarise_all(list(~ median(.x, na.rm = TRUE)))
		  
		  dfb <- df %>% 
		    dplyr::select(-prot_cover, -grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
		    dplyr::filter(!is.na(gene)) %>% 
		    dplyr::filter(!duplicated(.$gene))
		  
		  dfc <- df %>% 
		    dplyr::select(gene, prot_cover) %>% 
		    dplyr::filter(!is.na(gene), !is.na(prot_cover)) %>% 
		    dplyr::group_by(gene) %>% 
		    dplyr::mutate(prot_cover = as.numeric(sub("%", "", prot_cover))) %>% 
		    dplyr::summarise_all(~ max(.x, na.rm = TRUE)) %>% 
		    dplyr::mutate(prot_cover = paste0(prot_cover, "%"))
		  
		  # the number of unique gene in psm_data may be shorter than those in 
		  # dfa and dfb for the reason of empty "genes" in dfa 
		  df <- list(dfc, dfb, dfa) %>% 
		    purrr::reduce(right_join, by = "gene") %>% 
		    dplyr::filter(!is.na(gene), !duplicated(gene))
		}
	} else {
	  df <- read.csv(file.path(dat_dir, "Protein", "Protein.txt"), sep = "\t", 
	                 check.names = FALSE, header = TRUE, comment.char = "#") %>% 
	    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
	}
	
	quietly_out <- purrr::quietly(normMulGau)(
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
	
	purrr::walk(quietly_out[-1], write, 
	            file.path(dat_dir, "Protein\\log","prn_MulGau_log.csv"), append = TRUE)
	
	df <- quietly_out$result %>% 
	  dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0) %>% 
	  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
	  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
	  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
	  dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %T>% 
	  write.table(., file.path(dat_dir, "Protein", "Protein.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
}
