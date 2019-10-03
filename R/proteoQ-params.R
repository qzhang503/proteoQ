#' Processes the metadata of TMT experiments
#'
#' @import dplyr tidyr purrr openxlsx
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
prep_label_scheme <- function(dat_dir, filename) {

	my_channels <- function (x) {
		x <- as.character(x)
		pos <- !grepl("^TMT", x)
		x[pos] <- paste0("TMT-", x[pos])

		return(x)
	}

	replace_na_smpls <- function(x, prefix) {
	  i <- is.na(x)
	  replace(as.character(x), i, paste(prefix, seq_len(sum(i)), sep = "."))
	}
	
	not_trival <- function (x) {
	  ok <- (!is.na(x)) & (x != FALSE) & (x != 0)
	}
	
	
	if (is.null(dat_dir)) dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv),
	                                         error = function(e) 1)

	if (dat_dir == 1) stop("Set up the working directory first.")

	if (!file.exists(file.path(dat_dir, filename)))
	  stop(filename, " not found under '", dat_dir, "'.")

	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	if (fn_suffix %in% c("xls", "xlsx")) {
		label_scheme_full <- readxl::read_excel(file.path(dat_dir, filename), sheet = "Setup") %>%
												dplyr::filter(rowSums(!is.na(.)) > 0)
	} else if (fn_suffix == "csv") {
		label_scheme_full <- read.csv(file.path(dat_dir, filename), check.names = TRUE,
		                              header = TRUE, comment.char = "#", na.strings = c("", "NA")) %>%
												dplyr::filter(rowSums(!is.na(.)) > 0)
	} else {
		stop(filename, " needs to be '.xls' or '.xlsx'.")
	}
	
	label_scheme_full <- label_scheme_full %>% 
	  dplyr::mutate(Sample_ID = ifelse(grepl("^Empty\\.[0-9]+", Sample_ID), NA, Sample_ID))
	
	check_tmt126 <- label_scheme_full %>% 
	  dplyr::filter(TMT_Channel == "TMT-126") %>% 
	  dplyr::filter(is.na(TMT_Set)|is.na(LCMS_Injection))
	
	if (nrow(check_tmt126) > 0) {
	  stop("`TMT_Set` and/or `LCMS_Injection` indeces corresponding to `TMT-126` in `expt_smry.xlsx` cannot be empty.", 
	       call. = FALSE)
	}
	rm(check_tmt126)

	must_have <- c("TMT_Channel", "TMT_Set", "LCMS_Injection", "RAW_File",
								"Sample_ID", "Reference")

	missing_cols <- must_have[!must_have %in% names(label_scheme_full)]
	if (length(missing_cols) > 0) {
		purrr::walk(missing_cols, ~ cat(paste0("\'", ., "\' must be present in \'", filename, "\'\n")))
		stop("Not all required columns are present in \'", filename, "\'", call. = TRUE)
	}

	default_names <- c("Select", "Group", "Order", "Fill",  "Color", "Shape", "Size", "Alpha", "Peptide_Yield")

	purrr::walk(as.list(default_names), ~ {
		if (!.x %in% names(label_scheme_full)) {
		  message("Column \'", .x, "\' added to \'", filename, "\'")
			label_scheme_full[[.x]] <<- NA
		}
	}, label_scheme_full)

	# a case of label-free data
	if (dplyr::n_distinct(label_scheme_full$TMT_Channel) == 1) label_scheme_full$TMT_Channel <- NA

	TMT_plex <- TMT_plex(label_scheme_full)
	TMT_levels <- TMT_levels(TMT_plex)

	label_scheme_full <- label_scheme_full %>% 
	  dplyr::mutate_at(vars(c("TMT_Channel")), ~ my_channels(.x)) %>% 
	  dplyr::filter(rowSums(is.na(.)) < ncol(.)) %>% 
	  dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File, ignore.case = TRUE)) %>% 
		dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File, ignore.case = TRUE)) %>% # Bruker
	  dplyr::mutate_at(vars(c("Reference")), ~ not_trival(.x)) %>%
	  dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ as.numeric(.x)) %>%
	  dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ round(.x, digits = 2)) %>%
	  tidyr::fill(one_of("TMT_Set", "LCMS_Injection", "RAW_File")) %>%
	  dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
	  dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)

	# add IDs to unused TMT channels
	label_scheme_empty <- label_scheme_full %>%
		dplyr::select(TMT_Channel, TMT_Set, Sample_ID) %>%
		tidyr::unite(key, TMT_Channel, TMT_Set, remove = TRUE) %>%
		dplyr::filter(!duplicated(key)) %>%
		dplyr::mutate(Sample_ID = replace_na_smpls(Sample_ID, "Empty"))

	label_scheme_full <- label_scheme_full %>%
		tidyr::unite(key, TMT_Channel, TMT_Set, remove = FALSE) %>%
		dplyr::select(-Sample_ID) %>%
		dplyr::left_join(label_scheme_empty, by = "key") %>%
		dplyr::select(-key) %>%
		dplyr::mutate(Sample_ID = factor(Sample_ID))
	  
	label_scheme_full <- dplyr::bind_cols(
	  label_scheme_full %>% dplyr::select(Sample_ID), 
	  label_scheme_full %>% dplyr::select(-Sample_ID))

	rm(label_scheme_empty)

	# check the completeness of TMT_Channel
	check_tmt <- label_scheme_full %>%
		tidyr::complete(TMT_Set, LCMS_Injection, TMT_Channel) %>%
		dplyr::filter(is.na(RAW_File)) %>%
		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
		dplyr::mutate(n = n()) %>%
		dplyr::filter(n != TMT_plex)

	if (nrow(check_tmt) > 0) {
		check_tmt %>%
			dplyr::select(TMT_Set, LCMS_Injection, TMT_Channel) %>%
			print()

		stop(paste(check_tmt$TMT_Channel, "are missing at set", check_tmt$TMT_Set, "injection",
		           check_tmt$LCMS_Injection, "\n"))
	}

	# check the uniqueness of RAW_File per TMT_Set and LCMS_Injection
	check_fn <- label_scheme_full %>%
		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
		dplyr::summarise(count = n_distinct(RAW_File)) %>%
		dplyr::filter(count > 1) %>%
		dplyr::select(-count)

	if (nrow(check_fn) > 0) {
		check_fn %>% print()
		stop("More than one RAW filename in the above combination of TMT sets and LCMS injections.")
	}

	# check the uniqueness of Sample_ID
	check_smpls <- label_scheme_full %>%
		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
		dplyr::summarise(count = n_distinct(Sample_ID)) %>%
		dplyr::filter(count != TMT_plex)

	if (nrow(check_smpls) > 0 & TMT_plex > 0) {
		check_smpls %>% print()
		stop(paste("Need", TMT_plex,
		           "unique samples in the above combination of TMT sets and LCMS injections." ))
	}

	save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.Rdata"))

	wb <- openxlsx::loadWorkbook(file.path(dat_dir, filename))
	openxlsx::writeData(wb, sheet = "Setup", label_scheme_full)
	openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
	
	simple_label_scheme(dat_dir, label_scheme_full)
}


#' Loads the information of analyte prefractionation
#'
#' @import dplyr purrr tidyr openxlsx
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
prep_fraction_scheme <- function(dat_dir, filename) {
	if (is.null(dat_dir))
	  dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)

	if (dat_dir == 1) stop("Set up the working directory first.")

	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	if (file.exists(file.path(dat_dir, filename))) {
		if (fn_suffix %in% c("xls", "xlsx")) {
			fraction_scheme <- readxl::read_excel(file.path(dat_dir, filename), sheet = "Fractions") %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else if (fn_suffix == "csv") {
			fraction_scheme <- read.csv(file.path(dat_dir, filename), check.names = TRUE, header = TRUE,
			                            comment.char = "#", na.strings = c("", "NA"))  %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else {
			stop(filename, " needs to be in a file format of '.xls' or '.xlsx'.")
		}
	  
	  fraction_scheme <- fraction_scheme %>% 
	    tidyr::fill(TMT_Set, LCMS_Injection) %>% 
	    dplyr::group_by(TMT_Set, LCMS_Injection) %>% 
	    dplyr::mutate(Fraction = row_number()) 

	  wb <- openxlsx::loadWorkbook(file.path(dat_dir, filename))
	  openxlsx::writeData(wb, sheet = "Fractions", fraction_scheme)
	  openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
 	} else {
 	  assign(".auto_frac_smry", TRUE, envir = .GlobalEnv)
 	  
 	  # warning: data in a auto-generated `frac_smry.xlsx` will be incorrect 
 	  #   if they were based on wrong information from `expt_smry.xlsx`
 	  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
 	  
 	  # in case forget to enter RAW_File names
 	  if (anyNA(label_scheme_full$RAW_File)) stop("Enter RAW file names in the experimental summary file")

		fraction_scheme <- label_scheme_full %>%
			dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
			dplyr::filter(!duplicated(RAW_File)) %>%
			dplyr::group_by(TMT_Set, LCMS_Injection) %>%
			dplyr::mutate(Fraction = row_number())

		wb <- openxlsx::createWorkbook()
		openxlsx::addWorksheet(wb, sheetName = "Fractions")
		openxlsx::writeData(wb, sheet = "Fractions", fraction_scheme)
		openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
 	}
	
	save(fraction_scheme, file = file.path(dat_dir, "fraction_scheme.Rdata"))
}


#'Loads species-specific Databases
#'
#'A function loads a set of precompiled tables in protein accessions, \code{GO}
#'data, \code{KEGG} pathways and
#'\code{\href{http://software.broadinstitute.org/gsea/msigdb}{molecular
#'signatures}}.
#'
#'@seealso \code{\link{load_expts}} for supported species.
#'
#' @examples
#' load_dbs("human")
#'
#'@import dplyr rlang
#'@importFrom magrittr %>%
#'@export
load_dbs <- function (gset_nms = "go_sets", species = "human") {
  allowed <- c("go_sets", "kegg_sets", "c2_msig")
  stopifnot(all(gset_nms %in% allowed))
  
	if (is.null(species)) stop("Unrecognized species.")

  abbr_sp <- purrr::map_chr(species, sp_lookup)
  filelist <- map(abbr_sp, ~ paste0(gset_nms, "_", .x)) %>% unlist()
  
  # added abbr_sp to make GO terms species specific 
  # $`hs_GO:0000018 regulation of DNA recombination`
  data(package = "proteoQ", list = filelist)
  gsets <- purrr::map(filelist, ~ try(get(.x))) %>% do.call(`c`, .)
  
  stopifnot(length(gsets) > 0)
  try(rm(list = filelist, envir = .GlobalEnv))
  
  if (length(gsets) > 0) names(gsets) <- gsub("/", "-", names(gsets))
  
  assign("gsets", gsets, envir = .GlobalEnv)
} 



#'Load TMT experiments
#'
#'A function processes \code{.xlsx} files containing the metadata of TMT
#'experiments
#'
#'@section \code{expt_smry.xlsx}: The \code{expt_smry.xlsx} files should be
#'  located immediately under the file folder defined by the character vector
#'  \code{dat_dir}. The tab containing the metadata of TMT experiments should be
#'  named \code{Setup}. The \code{Excel} spread sheet therein is comprised of
#'  three tiers of fields: (1) essential, (2) optional default and (3) optional
#'  open. The \code{essential} columns contain the mandatory information of TMT
#'  experiments. The \code{optional default} columns serve as the fields for
#'  default lookups in sample selection, grouping, ordering, aesthetics, etc.
#'  The \code{optional open} fields allow users to define their own analysis,
#'  aesthetics, etc.
#'
#'  \tabular{ll}{ \strong{Essential column}   \tab \strong{Descrption}\cr
#'  Sample_ID \tab Unique sample IDs \cr TMT_Channel \tab TMT channel names:
#'  \code{126}, \code{127N}, \code{127C} et al. \cr TMT_Set \tab TMT experiment
#'  indeces  \cr LCMS_Injection   \tab LC/MS injection indeces under a
#'  \code{TMT_Set} \cr RAW_File \tab MS data filenames originated by
#'  \code{Xcalibur} with or without the \code{.raw} extension \cr Reference \tab
#'  Labels indicating reference samples in TMT experiments \cr }
#'
#'  \code{Sample_ID}: values should be unique for entries at a unique
#'  combination of \code{TMT_Channel} and \code{TMT_Set}, or left blank for
#'  unused entries. Samples with the same indeces of \code{TMT_Channel} and
#'  \code{TMT_Set} but different indeces of \code{LCMS_Injection} should have
#'  the same value in \code{Sample_ID}. No white space or special characters are
#'  allowed.
#'
#'  \code{RAW_File}: \code{OS} file names may be altered by MS users and thus
#'  different to those recorded in \code{Xcalibur}. The original names by
#'  \code{Xcalibur} should be used. For analysis with off-line fractionations of
#'  peptides before LC/MS, the \code{RAW_File} column should be left blank. The
#'  correspondence between the fractions and \code{RAW_File} names should be
#'  specified in a separate file, for example, \code{frac_smry.xlsx}. The
#'  utility \code{extract_raws(raw_dir)} can be used to extract the names of raw
#'  MS files under the directory indicated by `raw_dir`. More details can be
#'  found via \code{?extract_raws}. The utility \code{extract_psm_raws(dat_dir)}
#'  extracts the names of raw MS files that are actually present in
#'  \code{Mascot} PSM outputs.
#'
#'  \code{Reference}: reference entrie(s) are indicated with non-void string(s).
#'
#'  \tabular{ll}{ \strong{Optional default column}   \tab \strong{Descrption}\cr
#'  Select \tab Samples to be selected for indicated analysis \cr Group \tab
#'  Aesthetic labels annotating the prior knowledge of sample groups, e.g.,
#'  Ctrl_T1, Ctrl_T2, Disease_T1, Disease_T2, ...\cr Order \tab Numeric labels
#'  specifying the order of sample \code{groups} \cr Fill \tab Aesthetic labels
#'  for sample annotation by filled color\cr Color \tab Aesthetic labels for
#'  sample annotation by edge color\cr Shape \tab Aesthetic labels for sample
#'  annotation by shape\cr Size \tab Aesthetic labels for sample annotation by
#'  size \cr Alpha \tab Aesthetic labels for sample annotation by transparency
#'  \cr Benchmark \tab Depreciated; indicators of benchmark sample (groups) for
#'  use in heat map visualizations.\cr}
#'
#'  \tabular{ll}{ \strong{Exemplary optional open column}   \tab
#'  \strong{Descrption}\cr Term \tab Categorial terms for statistical modeling.
#'  \cr Duplicate \tab Indicators of duplicated samples for corrections in
#'  statistical significance \cr Peptide_Yield \tab Yields of peptides in sample
#'  handling \cr}
#'
#'
#'@section \code{frac_smry.xlsx}: It is not necessary to prepare a
#'  \code{frac_smry.xlsx} file if no peptide fractionations were performed in
#'  TMT experiments. 
#'
#'  \tabular{ll}{ \strong{Column}   \tab \strong{Descrption}\cr TMT_Set \tab
#'  v.s.  \cr LCMS_Injection   \tab v.s. \cr Fraction \tab Fraction indeces
#'  under a \code{TMT_Set} \cr RAW_File \tab v.s. }
#'
#'@param dat_dir A character string to the working directory.
#'@param expt_smry A character string to the \code{.xlsx} file containing the
#'  metadata of TMT experiments. The default is "expt_smry.xlsx".
#'@param frac_smry A character string to the \code{.xlsx} file containing
#'  peptide fractionation summary. The default is "frac_smry.xlsx".
#'
#' @examples
#' # An examplary "expt_smry.xlsx"
#' system.file("extdata", "expt_smry.xlsx", package = "proteoQ")
#'
#' # An examplary "frac_smry.xlsx"
#' system.file("extdata", "frac_smry.xlsx", package = "proteoQ")
#'
#' \dontrun{
#' # set up a working directory
#' dir.create("C:\\The\\Mascot\\Example", recursive = TRUE, showWarnings = FALSE)
#' dat_dir <- c("C:\\The\\Mascot\\Example")
#'
#' # copy fasta (do it once)
#' library(proteoQDA)
#' copy_refseq_hs("~\\proteoQ\\dbs\\fasta\\refseq")
#' copy_refseq_mm("~\\proteoQ\\dbs\\fasta\\refseq")
#'
#' # copy Mascot PSM data
#' cptac_csv_1(dat_dir)
#'
#' # copy "expt_smry.xlsx" and "frac_smry.xlsx"
#' cptac_expt_1(dat_dir)
#' cptac_frac_1(dat_dir)
#'
#' # load experiments
#' library(proteoQ)
#' load_expts()
#'
#'
#' # process PSMs with in-function filtration of data by `filter_`
#' normPSM(
#'   group_psm_by = pep_seq,
#'   group_pep_by = gene,
#'   fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
#'             "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
#'   rptr_intco = 3000,
#'   rm_craps = TRUE,
#'   rm_krts = FALSE,
#'   rm_outliers = FALSE,
#'   annot_kinases = TRUE,
#'   plot_rptr_int = TRUE,
#'   plot_log2FC_cv = TRUE,
#'
#'   filter_peps = exprs(pep_expect <= .1),
#'   filter_by_more = exprs(pep_rank == 1, pep_exp_z > 1),
#' )
#'
#' # exemplary PSM purging; n: the number of PSMs under a peptide
#' purgePSM(max_cv = .5, min_n = 2)
#'
#'
#' # peptides results with exemplary `filter_...`
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
#'   # filter_by = exprs(pep_n_psm >= 2),
#'   # filter_by_sp = exprs(species == "human"),
#' )
#'
#' # exemplary peptide purging; n: the number of peptides under a protein
#' purgePep(max_cv = .5, min_n = 2)
#'
#'
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
#'   # filter_by = exprs(prot_n_psm >= 5, prot_n_pep >=2),
#' )
#'
#'
#' # validation steps...
#' ?pepHist
#'
#' }
#'
#'@import dplyr rlang fs
#'@importFrom magrittr %>%
#'@export
load_expts <- function (dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx") {
  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))

  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  
  if (!fs::dir_exists(dat_dir)) {
    stop(dat_dir, " not existed.", call. = FALSE)
  }

	mget(names(formals()), rlang::current_env()) %>% save_call("load_expts")

  prep_label_scheme(dat_dir, expt_smry)
  prep_fraction_scheme(dat_dir, frac_smry)
}


#' Reload the "expt_smry.xlsx" and "frac_smry.xlsx"
#'
#' @import rlang
#' @importFrom magrittr %>%
#' @importFrom fs file_info
reload_expts <- function() {
  expt_smry <- match_expt()
  frac_smry <- match_frac()
  
  fi_xlsx <- fs::file_info(file.path(dat_dir, expt_smry))$change_time
  
  if (is.na(fi_xlsx)) stop("Time stamp of `expt_smry.xlsx` not available.")
  
  fi_rda <- fs::file_info(file.path(dat_dir, "label_scheme.Rdata"))$change_time
  if (fi_xlsx > fi_rda) {
    load_expts(dat_dir = dat_dir, expt_smry = !!expt_smry, frac_smry = !!frac_smry)
  }
}


#' Extracts the channel information in TMT experiments
#'
#' A function returns the indeces of TMT channels that are associated to
#' reference(s), sample(s) and probable unused void(s).
#'
#' @param label_scheme The data frame returned by \code{\link{load_expts}}.
#' @param set_idx Numeric.  The index of a multiplex TMT experiment.
#' @return Three lists of indeces: \code{refChannels}, reference channels(s);
#'   \code{emptyChannels}, empty channel(s) that were not used for sample
#'   labeling; \code{labeledChannels}, non-empty channels including both
#'   reference(s) and sample(s).
#'
#' @importFrom dplyr select filter
channelInfo <- function (label_scheme, set_idx) {
	stopifnot(length(set_idx) == 1)

	label_scheme_sub <- label_scheme %>%
	  dplyr::filter(!duplicated(Sample_ID), TMT_Set == set_idx)

	ref <- label_scheme_sub$Reference

	empty_channel_sub <- is.na(label_scheme_sub$Sample_ID) |
	  grepl("^Empty|^Outlier", label_scheme_sub$Sample_ID, ignore.case = TRUE)

	label_scheme_sub <- !empty_channel_sub

	out <- list(
		refChannels = ref,
		emptyChannels = empty_channel_sub,
		labeledChannels = label_scheme_sub
	)

	lapply(out, which)
}


#' Finds the number of multiplex TMT experiments
#'
#' \code{n_TMT_sets} returns the number of multiplex TMT experiments.
n_TMT_sets <- function (label_scheme_full) {
	length(unique(label_scheme_full$TMT_Set))
}

#' Finds the multiplxity of TMT labels
#'
#' \code{TMT_plex} returns the multiplxity of TMT labels.
TMT_plex <- function (label_scheme_full) {
	nlevels(as.factor(label_scheme_full$TMT_Channel))
}

#' Finds the factor levels of TMT labels
#'
#' \code{TMT_levels} returns the factor levels of TMT labels.
TMT_levels <- function (TMT_plex) {
	if(TMT_plex == 10) {
		TMT_levels <- c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", "TMT-129N",
		                "TMT-129C", "TMT-130N", "TMT-130C", "TMT-131")
	} else if(TMT_plex == 11) {
		TMT_levels <- c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", "TMT-129N",
		                "TMT-129C", "TMT-130N", "TMT-130C", "TMT-131N", "TMT-131C")
	} else if(TMT_plex == 6) {
		TMT_levels <- c("TMT-126", "TMT-127", "TMT-128", "TMT-129", "TMT-130", "TMT-131")
	} else if(TMT_plex == 1) {
		TMT_levels <- c("TMT-126")
	} else if(TMT_plex == 0) {
		TMT_levels <- NULL
	}
}

#' Simplifies label schemes from \code{label_scheme_full}
#'
#' Removes duplicated sample entries under different LC/MS injections.
simple_label_scheme <- function (dat_dir, label_scheme_full) {
	TMT_plex <- TMT_plex(label_scheme_full)
	TMT_levels <- TMT_levels(TMT_plex)

	label_scheme <- label_scheme_full %>%
		dplyr::filter(!duplicated(Sample_ID), !is.na(Sample_ID)) %>%
		dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
		dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)
	
	if (nrow(label_scheme) <(TMT_plex * n_TMT_sets(label_scheme))) {
	  stop("Duplicated sample ID(s) in `expt_smry.xlsx`", call. = FALSE)
	}

	save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))
}

#' Checks the uniqueness of sample IDs in \code{label_scheme_full}
#'
#' \code{check_label_scheme} will stop the analysis if the number of unique
#' samples are less than expected.
check_label_scheme <- function (label_scheme_full) {
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	TMT_plex <- TMT_plex(label_scheme)
	if(!is.null(TMT_plex)) {
		if((nlevels(as.factor(label_scheme$Sample_ID))) < 
		   (TMT_plex * nlevels(as.factor(label_scheme$TMT_Set))))
			stop("Not enough observations in unique 'Sample_ID'")
	}
}


#' Determine the protein accession type from PSM tables
#'
#' Find the protein accession from a non-cRAP entry and parse it.
#'
#' \code{parse_acc} parse the protein accession.
parse_acc <- function(df) {
  stopifnot("prot_acc" %in% names(df))
  
  prn_acc <- df %>%
    dplyr::filter(!grepl("^REV__", prot_acc)) %>% 
    dplyr::filter(!grepl("^CON__", prot_acc)) %>% 
    dplyr::select(prot_acc) %>%
    unlist %>%
    .[1]

  if (grepl("_[A-Z]+", prn_acc)) {
    acc_type <- "uniprot_id"
  } else if (grepl("^NP_[0-9]+", prn_acc)) {
    acc_type <- "refseq_acc"
  } else if (grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", prn_acc)) {
    acc_type <- "uniprot_acc"
  } else {
    stop("Unknown protein accession", call. = FALSE)
  }
  
  return(acc_type)
}


#' Match scale_log2r
#'
#' \code{match_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in caller environment.
match_scale_log2r <- function(scale_log2r) {
  stopifnot(rlang::is_logical(scale_log2r))

  global_var <-tryCatch(global_var <-get("scale_log2r", envir = .GlobalEnv),
                        error = function(e) "e")
  if(global_var != "e" & is.logical(global_var)) scale_log2r <- global_var

  return(scale_log2r)
}


#' Match to a global logical variable
#'
#' @examples
#' scale_log2r <- TRUE
#' foo <- function(scale_log2r = FALSE) {
#'   match_logi_gv("scale_log2r", scale_log2r)
#' }
#' foo()
match_logi_gv <- function(var, val) {
  gvar <-tryCatch(gvar <-get(var, envir = .GlobalEnv), error = function(e) "e")
  
  if(gvar != "e") {
    stopifnot(rlang::is_logical(gvar))
    return(gvar)
  } else {
    return(val)
  }
}


#' Find mismatches in RAW file names
#'
#' \code{check_raws} finds mismatched RAW files between expt_smry.xlsx and
#' PSM outputs.
check_raws <- function(df) {
  stopifnot ("RAW_File" %in% names(df))
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  load(file = file.path(dat_dir, "fraction_scheme.Rdata"))

  ## program-generated frac_smry.xlsx may be based on wrong information from expt_smry.xlsx
  ls_raws <- label_scheme_full$RAW_File %>% unique()
  fs_raws <- fraction_scheme$RAW_File %>% unique()
  if (!(all(is.na(ls_raws)) | all(ls_raws %in% fs_raws))) {
    pars <- read.csv(file.path(dat_dir, "Calls", "load_expts.txt"), 
                     check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")

    fn_frac <- pars %>% 
      dplyr::filter(var == "frac_smry") %>% 
      dplyr::select("value.1") %>% 
      unlist() %>% 
      as.character()
    
    unlink(file.path(dat_dir, fn_frac))
    prep_fraction_scheme(dat_dir, fn_frac)
    load(file = file.path(dat_dir, "fraction_scheme.Rdata"))
  }

  tmtinj_raw <- fraction_scheme %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>%
    dplyr::select(-Fraction) %>%
    dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
		dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File)) # Bruker
  
  ms_raws <- df$RAW_File %>% unique()
  label_scheme_raws <- tmtinj_raw$RAW_File %>% unique()
  
  missing_ms_raws <- ms_raws %>% .[! . %in% label_scheme_raws]
  wrong_label_scheme_raws <- label_scheme_raws[! label_scheme_raws %in% ms_raws]
  
  if(!purrr::is_empty(missing_ms_raws) | !purrr::is_empty(wrong_label_scheme_raws)) {
    cat("RAW MS file name(s) missing from the `experimental summary` and/or `fraction summary` files:\n")
    cat(paste0(missing_ms_raws, "\n"))
    
    cat("RAW MS files in `experimental summary` and/or `fraction summary` files not found in PSM data:\n")
    cat(paste0("\t", wrong_label_scheme_raws, "\n"))
    
    stop(paste("Check file names under the RAW_File column in the experimental summary file."),
         call. = FALSE)
  }

  return(tmtinj_raw)
}

