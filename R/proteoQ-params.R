#' Processes the metadata of experiments
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
prep_label_scheme <- function(dat_dir, filename) {

	my_channels <- function (x) {
		x <- as.character(x)

		pos <- !grepl("^TMT", x)
		x[pos] <- paste0("TMT-", x[pos])

		return(x)
	}


	if(is.null(dat_dir)) dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv),
	                                         error = function(e) 1)

	if(dat_dir == 1) stop("Set up the working directory first.")

	if(!file.exists(file.path(dat_dir, filename)))
	  stop(filename, " not found under '", dat_dir, "'.")

	fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)

	if(fn_suffix %in% c("xls", "xlsx")) {
		label_scheme_full <- readxl::read_excel(file.path(dat_dir, filename), sheet = "Setup") %>%
												dplyr::filter(rowSums(!is.na(.)) > 0)
	} else if(fn_suffix == "csv") {
		label_scheme_full <- read.csv(file.path(dat_dir, filename), check.names = TRUE,
		                              header = TRUE, comment.char = "#", na.strings = c("", "NA")) %>%
												dplyr::filter(rowSums(!is.na(.)) > 0)
	} else {
		stop(filename, " needs to be in a file format of '.csv', '.xls' or '.xlsx'.")
	}

	must_have <- c("TMT_Channel", "TMT_Set", "LCMS_Injection", "RAW_File",
								"Sample_ID", "Reference")

	missing_cols <- must_have[!must_have %in% names(label_scheme_full)]
	if(length(missing_cols) > 0) {
		purrr::walk(missing_cols, ~ cat(paste0("\'", ., "\' must be present in \'", filename, "\'\n")))
		stop("Not all required columns are present in \'", filename, "\'", call. = TRUE)
	}

	default_names <- c("Reference", "Select", "Group", "Order", "Fill",  "Color",
										"Shape", "Size", "Alpha", "Term", "Duplicate", "Benchmark")

	purrr::walk(default_names, ~ {
		if(!.x %in% names(label_scheme_full)) {
			message("Column \'", .x, "\' added to \'", filename, "\'")
			label_scheme_full[[.x]] <<- NA
		}
	})

	# a case of label-free data
	if(dplyr::n_distinct(label_scheme_full$TMT_Channel) == 1) label_scheme_full$TMT_Channel <- NA

	TMT_plex <- TMT_plex(label_scheme_full)
	TMT_levels <- TMT_levels(TMT_plex)

	label_scheme_full <- label_scheme_full %>%
		dplyr::mutate_at(vars(c("TMT_Channel")), ~ my_channels(.)) %>%
		dplyr::filter(rowSums(is.na(.)) < ncol(.)) %>%
		dplyr::mutate(Reference = !is.na(Reference),
		              RAW_File = gsub("\\.raw$", "", RAW_File, ignore.case = TRUE)) %>%
		dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ as.numeric(.)) %>%
		dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ round(., digits = 2)) %>%
		# dplyr::mutate_at(vars(one_of("DAT_File")), ~ gsub("\\.csv$", "", ., ignore.case = TRUE)) %>%
		# dplyr::mutate_at(vars(one_of("Windows_File")), ~ gsub("\\.raw$", "", ., ignore.case = TRUE)) %>%
		tidyr::fill(TMT_Set, LCMS_Injection, RAW_File) %>%
		dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
		dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)

	# replace the ids of NA samples
	replace_na_smpls <- function(x, prefix) {
		i <- is.na(x)
		replace(as.character(x), i, paste(prefix, seq_len(sum(i)), sep = "."))
	}

	label_scheme_temp <- label_scheme_full %>%
		dplyr::select(TMT_Channel, TMT_Set, Sample_ID) %>%
		tidyr::unite(key, TMT_Channel, TMT_Set, remove = TRUE) %>%
		dplyr::filter(!duplicated(key)) %>%
		dplyr::mutate(Sample_ID = replace_na_smpls(Sample_ID, "Empty"))

	label_scheme_full <- label_scheme_full %>%
		tidyr::unite(key, TMT_Channel, TMT_Set, remove = FALSE) %>%
		dplyr::select(-Sample_ID) %>%
		dplyr::left_join(label_scheme_temp, by = "key") %>%
		dplyr::select(-key) %>%
		dplyr::mutate(Sample_ID = factor(Sample_ID))

	rm(label_scheme_temp)

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
		           check_tmt$LCMS_Injection))
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
	simple_label_scheme(dat_dir, label_scheme_full)

	# load(file = file.path(dat_dir, "label_scheme_full.Rdata"), envir =  .GlobalEnv)
	# load(file = file.path(dat_dir, "label_scheme.Rdata"), envir =  .GlobalEnv)
}


#' Loads the information of analyte prefractionation
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
prep_fraction_scheme <- function(dat_dir, filename) {
	if(is.null(dat_dir))
	  dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)

	if(dat_dir == 1) stop("Set up the working directory first.")

	fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
	fn_prefix <- gsub("\\..*$", "", filename)

	if(file.exists(file.path(dat_dir, filename))) {
		if(fn_suffix %in% c("xls", "xlsx")) {
			fraction_scheme <- readxl::read_excel(file.path(dat_dir, filename), sheet = "Fractions") %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else if(fn_suffix == "csv") {
			fraction_scheme <- read.csv(file.path(dat_dir, filename), check.names = TRUE, header = TRUE,
			                            comment.char = "#", na.strings = c("", "NA"))  %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else {
			stop(filename, " needs to be in a file format of '.csv', '.xls' or '.xlsx'.")
		}
 	} else {
		# load(file = file.path(dat_dir, "label_scheme_full.Rdata"), envir =  .GlobalEnv)
		load(file = file.path(dat_dir, "label_scheme_full.Rdata"))

		fraction_scheme <- label_scheme_full %>%
			dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
			dplyr::filter(!duplicated(RAW_File)) %>%
			dplyr::group_by(TMT_Set, LCMS_Injection) %>%
			dplyr::mutate(Fraction = row_number())
	}

	write.csv(fraction_scheme, file.path(dat_dir, paste0(fn_prefix, ".csv")), row.names = FALSE)
	save(fraction_scheme, file = file.path(dat_dir, "fraction_scheme.Rdata"))
	# load(file = file.path(dat_dir, "fraction_scheme.Rdata"), envir = .GlobalEnv)
}


#'Loads species-specific Databases
#'
#'A function loads a set of precompiled tables in protein accessions, \code{GO}
#'data, \code{KEGG} pathways and
#'\code{\href{http://software.broadinstitute.org/gsea/msigdb}{molecular
#'signatures}}.
#'
#'@seealso \code{\link{setup_expts}} for supported species.
#'
#' @examples
#' load_dbs()
#'
#'@export
#'@import dplyr
#'@importFrom magrittr %>%
load_dbs <- function (dat_dir, expt_smry = "expt_smry.xlsx") {
	species <- find_species(label_scheme)

	if(is.null(species)) stop("'species' in '", expt_smry,
	                          "' not recognized. It needs to be one of 'human', 'mouse', 'rat' or 'pdx'.")

	if (species %in% c("homo sapiens", "human")) {
		data(package = "proteoQ", prn_annot_hs)
		data(package = "proteoQ", go_sets_hs)
		data(package = "proteoQ", kegg_sets_hs)
		data(package = "proteoQ", c2_msig_hs)

		prn_annot <- prn_annot_hs
		go_sets <- go_sets_hs
		kegg_sets <- kegg_sets_hs
		c2_sets <- c2_msig_hs

		kegg_species <- "hsa"
		go_species <- "Human"

		rm(prn_annot_hs, go_sets_hs, kegg_sets_hs, c2_msig_hs, envir = .GlobalEnv)
	} else if (species %in% c("mus musculus", "mouse")) {
		data(package = "proteoQ", prn_annot_mm)
		data(package = "proteoQ", go_sets_mm)
		data(package = "proteoQ", kegg_sets_mm)

		prn_annot <- prn_annot_mm
		go_sets <- go_sets_mm
		kegg_sets <- kegg_sets_mm
		c2_sets <- NULL

		kegg_species <- "mmu"
		go_species <- "Mouse"

		rm(prn_annot_mm, go_sets_mm, kegg_sets_mm, envir = .GlobalEnv)
	} else if (species %in% c("rattus norvegicus", "rat")) {
		data(package = "proteoQ", prn_annot_rn)
		data(package = "proteoQ", go_sets_rn)
		data(package = "proteoQ", kegg_sets_rn)

		prn_annot <- prn_annot_rn
		go_sets <- go_sets_rn
		kegg_sets <- kegg_sets_rn
		c2_sets <- NULL

		kegg_species <- "rno"
		go_species <- "Rat"

		rm(prn_annot_rn, go_sets_rn, kegg_sets_rn, envir = .GlobalEnv)
	} else if (species %in% c("drosophila", "fly")) {
		# load(file = ".\\data\\prn_annot_dm.RData")
		data(package = "proteoQ", prn_annot_dm)
		go_sets_dm <- NULL
		kegg_sets_dm <- NULL

		prn_annot <- prn_annot_dm
		go_sets <- go_sets_dm
		kegg_sets <- kegg_sets_dm
		c2_sets <- NULL

		kegg_species <- "dme"
		go_species <- "Fly"

		rm(prn_annot_dm, envir = .GlobalEnv)
	} else if(species %in% c("arabidopsis")) {
		prn_annot_at <- NULL
		go_sets_at <- NULL
		kegg_sets_at <- NULL
		c2_sets <- NULL

		prn_annot <- prn_annot_at
		go_sets <- go_sets_at
		kegg_sets <- kegg_sets_at

		kegg_species <- "ath"
		go_species <- "Anopheles"
	} else if (species %in% c("human and mouse", "pdx")) {
		data(package = "proteoQ", prn_annot_hs)
		data(package = "proteoQ", prn_annot_mm)
		prn_annot_hsmm <- rbind(prn_annot_hs, prn_annot_mm)

		data(package = "proteoQ", go_sets_hs)
		data(package = "proteoQ", go_sets_mm)

		names(go_sets_hs) <- paste("Hs", names(go_sets_hs), sep = ".")
		names(go_sets_mm) <- paste("Mm", names(go_sets_mm), sep = ".")
		go_sets_hsmm <- c(go_sets_hs, go_sets_mm)

		data(package = "proteoQ", kegg_sets_hs)
		data(package = "proteoQ", kegg_sets_mm)

		kegg_sets_hsmm <- c(kegg_sets_hs, kegg_sets_mm)

		data(package = "proteoQ", c2_msig_hs)
		c2_sets <- c2_msig_hs

		prn_annot <- prn_annot_hsmm
		go_sets <- go_sets_hsmm
		kegg_sets <- kegg_sets_hsmm

		kegg_species <- "hsammu"
		go_species <- "Human and Mouse"

		rm(prn_annot_hs, go_sets_hs, kegg_sets_hs, c2_msig_hs, envir = .GlobalEnv)
		rm(prn_annot_mm, go_sets_mm, kegg_sets_mm, envir = .GlobalEnv)
	}

	if(!is.na(go_species)) names(go_sets) <- gsub("/", "-", names(go_sets))
	if(!is.na(kegg_species)) names(kegg_sets) <- gsub("/", "~", names(kegg_sets))

	assign("dbs",
	       list(prn_annot = prn_annot, go_sets = go_sets, kegg_sets = kegg_sets, c2_msig = c2_sets),
	       envir = .GlobalEnv)
}


#'Load experiments
#'
#'A function processes \code{.xlsx} files containing the metadata of TMT
#'experiments
#'
#'@section \code{expt_smry.xlsx}: The \code{.xlsx} files should be located
#'  immediately under the file folder defined by the character vector of
#'  \code{dat_dir}.
#'
#'  In the \code{expt_smry.xlsx} file, values under the \code{Sample_ID} column
#'  should be unique for entries with unique combination of \code{TMT_Channel}
#'  and \code{TMT_Set}, or left blank for used entries. Samples with the same
#'  indeces of \code{TMT_Channel} and \code{TMT_Set} but different indeces of
#'  \code{LCMS_Injection} should have the same value in \code{Sample_ID}.
#'
#'  Under the \code{Reference} column, reference entrie(s) are indicated with
#'  non-void string(s).
#'
#'  Users can add or modify the \code{optional column keys} in the \code{.xlsx}
#'  file for suitable analysis (see \code{\link{prnSig}}... for the
#'  customization of column keys).
#'
#'  \tabular{ll}{ \strong{Reserved column key}   \tab \strong{Descrption}\cr
#'  TMT_Channel \tab TMT channel names: \code{126}, \code{127N}, \code{127C} et
#'  al. \cr TMT_Set \tab TMT experiment indeces  \cr LCMS_Injection   \tab LC/MS
#'  injection indeces under a \code{TMT_Set} \cr RAW_File \tab MS data filenames
#'  originated by \code{Xcalibur} with or without the \code{.raw} extension
#'  (\strong{Notes}: (1) \code{OS} filenames can be different to those recorded
#'  in \code{Xcalibur} and the later should be used. (2) For analysis with
#'  offline fractionations of peptides before LC/MS, the \code{RAW_File} column
#'  should be left blank. The correspondence between the fractions and
#'  \code{RAW_File} names should be specified in a separate file such as
#'  \code{fraction_scheme.xlsx}.)\cr Sample_ID \tab Unique sample IDs
#'  (\strong{Note}: no space in name; leave the sample IDs blank for unused TMT
#'  channels.) \cr Reference \tab Labels indicating reference samples in TMT
#'  experiments (\strong{Note}: leave the labels blank for non-reference TMT
#'  channels.) \cr \cr }
#'
#'  \tabular{ll}{ \strong{Optional column key}   \tab \strong{Descrption}\cr
#'  Select \tab Selected \code{Sample_IDs} for indicated analysis of "MDS",
#'  "Heat map"... If the values underneath are left blank, all samples including
#'  references will be used. \cr Group \tab Aesthetic labels annotating the
#'  prior knowledge of sample groups, e.g., Ctrl_T1, Ctrl_T2, Disease_T1,
#'  Disease_T2, ...\cr Order \tab Numeric labels specifying the order of sample
#'  \code{groups} (\strong{Note}: leave the labels blank for unsupervised
#'  analysis.)\cr Fill \tab Aesthetic labels for sample annotation by filled
#'  color\cr Color \tab Aesthetic labels for sample annotation by edge color\cr
#'  Shape \tab Aesthetic labels for sample annotation by shape\cr Size \tab
#'  Aesthetic labels for sample annotation by size \cr Alpha \tab Aesthetic
#'  labels for sample annotation by transparency \cr Term \tab Categorial terms
#'  for statistical modeling. Users can add columns, e.g. "Time", "Dose",
#'  "Time_and_Dose"..., for additional groups of contrasts. \cr Duplicate \tab
#'  Indicators of duplicated samples for corrections in statistical significance
#'  \cr Benchmark \tab Indicators of benchmark sample (groups) for use in heat
#'  map visualizations.\cr ... \tab Customed columns of metadata for aesthetic
#'  mapping to indicated analysis}
#'
#'@section \code{frac_smry.xlsx}: It is not necessary to prepare a
#'  "frac_smry.xlsx" file if no peptide fractionations were performed in TMT
#'  experiments.
#'
#'  \tabular{ll}{ \strong{Column key}   \tab \strong{Descrption}\cr TMT_Set \tab
#'  v.s.  \cr LCMS_Injection   \tab v.s. \cr Fraction \tab Fraction indeces
#'  under a \code{TMT_Set} \cr RAW_File \tab v.s. }
#'
#'@param dat_dir A character string to the working directory.
#'@param expt_smry A character string to the file of TMT experimental summary.
#'  The default is "expt_smry.xlsx".
#'@param frac_smry A character string to the file of peptide fractionation
#'  summary. The default is "frac_smry.xlsx".
#'
#' @examples
#' # An examplary "expt_smry.xlsx"
#' system.file("extdata", "expt_smry.xlsx", package = "proteoQ")
#'
#' # An examplary "frac_smry.xlsx"
#' system.file("extdata", "frac_smry.xlsx", package = "proteoQ")
#'
#' \dontrun{
#' dat_dir <- c("C:\\my_direcotry")
#' load_expts()
#' }
#'
#'@export
#'@import dplyr
#'@importFrom magrittr %>%
load_expts <- function (dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx") {
  if(is.null(dat_dir)) dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv),
                                           error = function(e) 1)

  if(dat_dir == 1) stop("Set up the working directory first.", call. = FALSE)

  if(!file.exists(file.path(dat_dir, "acctype_sp.txt"))) {
    acctype_sp <- data.frame(Accession_Type = "uniprot_id", Species = "human")
    # stop("The information of protein accession type and species not availble.\n Please run `normPSM()` first.")
  } else {
    acctype_sp <- read.csv(file.path(dat_dir, "acctype_sp.txt"), check.names = FALSE,
                           sep = "\t", header = TRUE, comment.char = "#")

  }

  prep_label_scheme(dat_dir, expt_smry)
  # prep_fraction_scheme(dat_dir = dat_dir, filename = frac_smry)

  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))

  label_scheme_full <- label_scheme_full %>%
    dplyr::mutate(Accession_Type = acctype_sp[1, "Accession_Type"], Species = acctype_sp[1, "Species"])
  save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.Rdata"))

  label_scheme <- label_scheme %>%
    dplyr::mutate(Accession_Type = acctype_sp[1, "Accession_Type"], Species = acctype_sp[1, "Species"])
  save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))

  load_dbs(dat_dir = dat_dir, expt_smry = expt_smry)
}


#' Extracts the channel information in TMT experiments
#'
#' A function returns the indeces of TMT channels that are associated to
#' reference(s), sample(s) and probable unused void(s).
#'
#' @param label_scheme The data frame returned by \code{\link{setup_expts}}.
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
		if((nlevels(as.factor(label_scheme$Sample_ID))) < (TMT_plex * nlevels(as.factor(label_scheme$TMT_Set))))
			stop("Not enough observations in unique 'Sample_ID'")
	}
}

#' Finds the Accession_Type
#'
#' \code{find_acctype} finds the \code{Accession_Type} that will be used in the
#' annotation of protein IDs.
find_acctype <- function (label_scheme_full) {
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	label_scheme %>%
		dplyr::filter(!is.na(Accession_Type), dplyr::row_number() == 1) %>%
		dplyr::select(Accession_Type) %>%
		unlist() %>%
		as.character() %>%
		tolower()
}

#' Finds the species
#'
#' \code{find_species} find the species that will be used in analysis such as GSEA.
find_species <- function (label_scheme_full) {
	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	label_scheme %>%
		dplyr::filter(!is.na(Species), dplyr::row_number() == 1) %>%
		dplyr::select(Species) %>%
		unlist() %>%
		as.character() %>%
		tolower()
}


#' Determine the protein accession type from PSM tables
#'
#' Find the protein accession from a non-cRAP entry and parse it.
#'
#' \code{parse_acc} parse the protein accession.
parse_acc_temp <- function(df, prx) {
  if (grepl("[0-9]{1}::", prx)) { # multiple databases
    prn_acc <- df %>%
      dplyr::filter(grepl(paste0("^", prx), prot_acc)) %>%
      dplyr::mutate(prot_acc = gsub(prx, "", prot_acc))
  } else { # single database
    prn_acc <- df
  }

  prn_acc <- prn_acc %>%
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
}


#' Determine the protein accession type from PSM tables
#'
#' Find the protein accession from a non-cRAP entry and parse it.
#'
#' \code{parse_acc} parse the protein accession.
parse_acc <- function(df) {
  prn_acc <- df %>%
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
}


#' Match scale_log2r
#'
#' \code{match_scale_log2r} matches `scale_log2r` from the caller environment.
match_scale_log2r <- function(scale_log2r) {
  stopifnot(is.logical(scale_log2r))

  global_var <-tryCatch(global_var <-get("scale_log2r", envir = .GlobalEnv),
                        error = function(e) "e")
  if(global_var != "e" & is.logical(global_var)) scale_log2r <- global_var

  return(scale_log2r)
}


#' Match to a global logical variable
#'
#' @examples
#' foo <- function(scale_log2r = FALSE) {
#'   match_logi_gv(scale_log2r)
#' }
#' foo()
match_logi_gv <- function(var) {
  nm <- rlang::as_string(rlang::enexpr(var))
  gvar <-tryCatch(gvar <-get(nm, envir = .GlobalEnv), error = function(e) "e")
  if(gvar != "e" & is.logical(gvar)) return(gvar) else return(var)
}
