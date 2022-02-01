#' Processes the metadata of TMT experiments
#'
#' @inheritParams load_expts
#' @inheritParams prnHist
#' @import dplyr tidyr purrr openxlsx
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom readxl read_excel
prep_label_scheme <- function(dat_dir = NULL, filename = "expt_smry.xlsx") 
{
	my_channels <- function (x) {
		x <- as.character(x)
		pos <- !grepl("^TMT", x)
		x[pos] <- paste0("TMT-", x[pos])

		invisible(x)
	}

	replace_na_smpls <- function(x, prefix) {
	  i <- is.na(x)
	  replace(as.character(x), i, paste(prefix, seq_len(sum(i)), sep = "."))
	}
	
	not_trival <- function (x) {
	  ok <- (!is.na(x)) & (x != FALSE) & (x != 0)
	}
	
	rm_empty_expts <- function (label_scheme_full, TMT_plex) {
	  label_scheme_full %>% 
	    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = FALSE) %>% 
	    dplyr::group_by(TMT_inj) %>% 
	    dplyr::mutate(.count = sum(is.na(Sample_ID))) %>% 
	    dplyr::filter(.count != TMT_plex) %>% 
	    dplyr::ungroup() %>% 
	    dplyr::select(-TMT_inj, -.count)
	}
	
	
	if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()

	if (!file.exists(file.path(dat_dir, filename)))
	  stop(filename, " not found under '", dat_dir, "'.", call. = FALSE)

	fn_suffix <- gsub("^.*\\.([^\\.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^\\.]*$", "", filename)

	if (fn_suffix %in% c("xls", "xlsx")) {
		label_scheme_full <- 
		  readxl::read_excel(file.path(dat_dir, filename), sheet = "Setup") %>% 
		  dplyr::filter(rowSums(!is.na(.)) > 0)
	} else if (fn_suffix == "csv") {
		label_scheme_full <- 
		  read.csv(file.path(dat_dir, filename), check.names = TRUE, 
		           header = TRUE, comment.char = "#", na.strings = c("", "NA")) %>% 
		  dplyr::filter(rowSums(!is.na(.)) > 0)
	} else {
		stop(filename, " needs to be '.xls' or '.xlsx'.")
	}

	local({
	  cols_noname <- grepl("^\\.{3}", names(label_scheme_full))
	  nms <- colnames(label_scheme_full)[cols_noname]
	  
	  if (length(nms)) 
	    warning("Missing column name in `", filename, "`: ", paste(nms, collapse = ", "), 
	            call. = FALSE)
	})
	
	local({
  	must_have <- c("TMT_Channel", "TMT_Set", "LCMS_Injection", "RAW_File",
  								"Sample_ID", "Reference")
  	
  	missing_cols <- must_have %>% .[! . %in% names(label_scheme_full)]
  	
  	if (length(missing_cols)) {
  	  purrr::walk(missing_cols, 
  	              ~ message(paste0("\'", ., "\' must be present in \'", 
  	                               filename, "\'\n")))
  	  stop("Not all required columns are present in \'", filename, "\'", 
  	       call. = FALSE)
  	}
	})

	default_names <- c("Select", "Group", "Order", "Fill",  "Color", 
	                   "Shape", "Size", "Alpha", "Peptide_Yield")

	purrr::walk(as.list(default_names), ~ {
		if (!.x %in% names(label_scheme_full)) {
		  message("Column \'", .x, "\' added to \'", filename, "\'")
			label_scheme_full[[.x]] <<- NA
		}
	}, label_scheme_full)

	label_scheme_full <- label_scheme_full %>% 
	  dplyr::mutate(Sample_ID = ifelse(grepl("^Empty\\.[0-9]+", Sample_ID), 
	                                   NA, Sample_ID))

	local({
  	check_tmt126 <- label_scheme_full %>% 
  	  dplyr::filter(TMT_Channel == "TMT-126" || TMT_Channel == "126") %>% 
  	  dplyr::filter(is.na(TMT_Set) || is.na(LCMS_Injection))
  	
  	if (nrow(check_tmt126) > 0) {
  	  stop("The indexes of `TMT_Set` and/or `LCMS_Injection` corresponding to \n", 
  	       "  the `TMT-126` rows in `expt_smry.xlsx` cannot be empty.", 
  	       call. = FALSE)
  	}
	})

	# --- standardize TMT channel names
	label_scheme_full <- local({
	  tmt_pairs <- tibble::tibble(
	    tmt127 = c("127", "127N"), 
	    tmt128 = c("128", "128N"), 
	    tmt129 = c("129", "129N"),
	    tmt130 = c("130", "130N"), 
	    tmt131 = c("131", "131N"),
	  )
	  
	  unique_tmt <- unique(label_scheme_full$TMT_Channel) %>% 
	    gsub("^TMT-", "", .)
	  
	  for (tmt_pair in tmt_pairs) {
	    tmt_pair <- unlist(tmt_pair)
	    
	    if (tmt_pair[1] %in% unique_tmt && tmt_pair[2] %in% unique_tmt) {
	      label_scheme_full <- label_scheme_full %>% 
	        dplyr::mutate(TMT_Channel = gsub(paste0(tmt_pair[1], "$"), 
	                                         tmt_pair[2], 
	                                         TMT_Channel))
	    }
	  }
	  
	  label_scheme_full <- label_scheme_full %>% 
	    dplyr::mutate(TMT_Channel = gsub("134$", "134N", TMT_Channel)) %>% 
	    dplyr::mutate(TMT_Channel = gsub("126N$", "126", TMT_Channel))
	  
	  return(label_scheme_full)
	})
	
	TMT_plex <- TMT_plex(label_scheme_full)
	
	if (!TMT_plex %in% c(0, 6, 10, 11, 16, 18)) {
	  stop("Apparent `TMT_plex = ", TMT_plex, 
	       "` is not one of c(0, 6, 10, 11, 16, 18).\n", 
	       "For non-TMT experiments (LFQ or simple qualitative processsing), \n", 
	       "  leave the fields under `TMT_Channel` blank.\n", 
	       "For other TMT-plexes, ", 
	       "pad the `TMT_Channel` upstream to an available multiplicity: \n", 
	       "  e.g. 1 -> 6, 8 -> 10 (regular TMT) or 9 -> 18, 15 -> 18 (TMTpro).", 
	       call. = FALSE)
	}

	if (TMT_plex == 10) {
	  label_scheme_full <- label_scheme_full %>% 
	    dplyr::mutate(TMT_Channel = gsub("126N$", "126", TMT_Channel)) %>% 
	    dplyr::mutate(TMT_Channel = gsub("131N$", "131", TMT_Channel))
	} else if (TMT_plex == 6) {
	  label_scheme_full <- label_scheme_full %>% 
	    dplyr::mutate(TMT_Channel = gsub("(12[6-9]{1})N$", "\\1", TMT_Channel)) %>% 
	    dplyr::mutate(TMT_Channel = gsub("(13[0-1]{1})N$", "\\1", TMT_Channel)) 
	}
	
	# --- auto fill
	TMT_levels <- TMT_levels(TMT_plex)
	
	label_scheme_full <- label_scheme_full %>% 
	  dplyr::filter(rowSums(is.na(.)) < ncol(.)) %>% 
	  dplyr::mutate_at(vars(c("TMT_Channel")), 
	                   ~ if (is.null(TMT_levels)) .x else my_channels(.x)) %>% 
	  dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File, ignore.case = TRUE)) %>% 
		dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File, ignore.case = TRUE)) %>% # Bruker
	  dplyr::mutate_at(vars(c("Reference")), ~ not_trival(.x)) %>%
	  dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ as.numeric(.x)) %>%
	  dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ round(.x, digits = 2)) %>%
	  tidyr::fill(one_of("TMT_Set", "LCMS_Injection", "RAW_File")) %>%
	  
	  tidyr::complete(TMT_Set, LCMS_Injection, TMT_Channel) %>%
	  { if (is.null(TMT_levels)) . else rm_empty_expts(., TMT_plex) } %>% 
	  tidyr::fill(one_of("TMT_Set", "LCMS_Injection", "RAW_File")) %>% 
	  dplyr::mutate_at(vars(c("Reference")), ~ not_trival(.x)) %>%
	  
	  dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
	  dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)
	
	if (!(is.null(TMT_levels) || rlang::is_integerish(label_scheme_full$TMT_Set))) {
	  stop("Values under `expt_smry.xlsx::TMT_Set` need to be integers.", 
	       call. = FALSE)
	}
	
	if (!is.null(TMT_levels)) {
	  if (!rlang::is_integerish(label_scheme_full$LCMS_Injection)) {
	    stop("Values under `expt_smry.xlsx::LCMS_Injection` need to be integers.", 
	         call. = FALSE)
	  }
	} else {
	  if (all(is.na(label_scheme_full$LCMS_Injection))) {
	    label_scheme_full$LCMS_Injection <- 1
	  }
	}

	# Check Sample_ID
	if (!is.null(TMT_levels)) {
  	# add IDs to unused TMT channels
  	label_scheme_full <- local({
    	
  	  # the first pass - all possible Empty.1, Empty.2 ... Sample_IDs
  	  label_scheme_empty <- local({
  	    label_scheme_empty <- label_scheme_full %>% 
  	      dplyr::filter(is.na(Sample_ID))
  	    
  	    empty_smpls <- label_scheme_empty %>%
  	      tidyr::unite(key, TMT_Channel, TMT_Set, remove = FALSE) %>%
  	      dplyr::filter(!duplicated(key)) %>%
  	      dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel) %>% 
  	      dplyr::mutate(Sample_ID_new = replace_na_smpls(Sample_ID, "Empty")) %>% 
  	      dplyr::select(TMT_Set, LCMS_Injection, TMT_Channel, Sample_ID_new) %>% 
  	      tidyr::unite(key, TMT_Channel, TMT_Set, LCMS_Injection, remove = TRUE)
  	    
  	    label_scheme_empty <- label_scheme_empty %>% 
  	      tidyr::unite(key, TMT_Channel, TMT_Set, LCMS_Injection, remove = FALSE) %>% 
  	      dplyr::left_join(empty_smpls, id = "key") %>% 
  	      dplyr::mutate(Sample_ID = Sample_ID_new) %>% 
  	      dplyr::select(-Sample_ID_new) %>% 
  	      dplyr::select(-key) %>% 
  	      tidyr::unite(key_2, TMT_Channel, TMT_Set, remove = FALSE)    
  	  })
  	  
  	  # the second pass - different Empty Sample_IDs at different LCMS_Injections
  	  label_scheme_empty <- local({
  	    empty_smpls_2 <- label_scheme_empty %>% dplyr::filter(is.na(Sample_ID))
  	    
  	    empty_smpls_2c <- label_scheme_empty %>% 
  	      dplyr::filter(!is.na(Sample_ID), key_2 %in% empty_smpls_2$key_2) %>% 
  	      dplyr::select(key_2, Sample_ID)
  	    
  	    empty_smpls_2 <- empty_smpls_2 %>% 
  	      dplyr::rename(Sample_ID_old = Sample_ID) %>% 
  	      dplyr::left_join(empty_smpls_2c, by = "key_2") %>% 
  	      dplyr::mutate(Sample_ID_old = Sample_ID) %>% 
  	      dplyr::select(-Sample_ID) %>% 
  	      dplyr::rename(Sample_ID = Sample_ID_old)
  	    
  	    label_scheme_empty <- label_scheme_empty %>% 
  	      dplyr::filter(!is.na(Sample_ID)) %>% 
  	      dplyr::bind_rows(empty_smpls_2) %>% 
  	      dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel) %>% 
  	      dplyr::select(-key_2)
  	  })
  	  
  	  label_scheme_full <- label_scheme_full %>% 
  	    dplyr::filter(!is.na(Sample_ID)) %>% 
  	    dplyr::bind_rows(label_scheme_empty) %>% 
  	    dplyr::mutate(Sample_ID = factor(Sample_ID)) %>% 
  	    dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)
  	  
    	label_scheme_full <- dplyr::bind_cols(
    	  label_scheme_full %>% dplyr::select(Sample_ID), 
    	  label_scheme_full %>% dplyr::select(-Sample_ID))
    	
    	# clean up of unused TMT channels
    	# Once Sample_ID is empty, Reference <- FALSE and opt_nms <- NA
    	label_scheme_full <- local({
      	opt_nms <- names(label_scheme_full) %>% 
      	  .[! . %in% c("TMT_Channel", "TMT_Set", "LCMS_Injection", 
      	               "RAW_File", "Sample_ID", "Reference")]
      	
      	empty_rows <- label_scheme_full %>% 
      	  dplyr::filter(grepl("^Empty\\.[0-9]+", Sample_ID))
      	
      	if (nrow(empty_rows) > 0 ) {
        	empty_rows[, opt_nms] <- NA
        	empty_rows[, "Reference"] <- FALSE  	  
      	}
    
      	label_scheme_full <- label_scheme_full %>% 
      	  dplyr::filter(! Sample_ID %in% empty_rows$Sample_ID) %>% 
      	  dplyr::bind_rows(empty_rows) %>% 
      	  dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)    	  
    	})    	
  	})
  	
  	# rows may be deleted later: check the completeness of TMT_Channel
  	# (RAW_File is either all NA or all filled)
  	local({
    	check_tmt <- label_scheme_full %>%
    		tidyr::complete(TMT_Set, LCMS_Injection, TMT_Channel) %>%
    		dplyr::filter(is.na(RAW_File)) %>%
    		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
    		dplyr::mutate(n = n()) %>%
    		dplyr::filter(n != TMT_plex)
    
    	if (nrow(check_tmt)) {
    	  message("=======================================")
    	  message("=== Mismatched TMT channel & plexes ===")
    		check_tmt %>%
    			dplyr::select(TMT_Set, LCMS_Injection, TMT_Channel) %>%
    	    data.frame(check.names = FALSE) %>% 
    			print()
    		message("=======================================")
    		
    		stop("Some parsed `TMT_plex` are different to `", TMT_plex, "`.", 
    		     "\n(E.g., Use `TMT-131`, instead of `TMT-131N`, ", 
    		     "for 10-plex experiment(s).)", 
    		     call. = FALSE)
    	}
  	})
  	
  	# check the uniqueness of RAW_File per TMT_Set and LCMS_Injection
  	local({
    	check_fn <- label_scheme_full %>%
    		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
    		dplyr::summarise(count = n_distinct(RAW_File)) %>%
    		dplyr::filter(count > 1) %>%
    		dplyr::select(-count)
    
    	if (nrow(check_fn) > 0) {
    		check_fn %>% print()
    		stop("More than one unique RAW filename in the above combination of ", 
    		     "TMT sets and LCMS injections.", 
    		     call. = FALSE)
    	}	  
  	})
  
  	# check the counts of Sample_ID under each TMT_Set and LCMS_Injection
  	# should equal to the TMT_plex
  	local({
    	check_smpls <- label_scheme_full %>%
    		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
    		dplyr::summarise(count = n_distinct(Sample_ID)) %>%
    		dplyr::filter(count != TMT_plex)
    
    	if (nrow(check_smpls) > 0 && TMT_plex > 0) {
    		check_smpls %>% print()
    	  stop("Need ", TMT_plex, " unique sample IDs ", 
    	       "in the above combination of TMT sets and LCMS injections.", 
    	       call. = FALSE)
    	}	  
  	})
  	
  	# check the uniqueness of Sample_ID 
  	# under the same TMT_Set but different LCMS_Injection; 
  	# OK one is `sample_1` and the other is `Empty.xxx`
  	local({
  	  dups <- label_scheme_full %>% 
  	    dplyr::filter(!grepl("^Empty\\.[0-9]+$", Sample_ID)) %>% 
  	    dplyr::group_by(TMT_Set, TMT_Channel) %>% 
  	    dplyr::summarise(count = n_distinct(Sample_ID)) %>% 
  	    dplyr::filter(count > 1) 
  	  
  	  if (nrow(dups) > 1L) {
  	    stop("Fix the above duplication(s) in sample IDs.\n", 
  	         print(data.frame(dups)),
  	         call. = FALSE)
  	  }
  	})
  	
  	if (!all(unique(label_scheme_full$TMT_Channel) %in% TMT_levels)) {
  	  stop("Not all `TMT_Channel` in `expt_smry.xlsx` in \n", 
  	       purrr::reduce(TMT_levels, paste, sep = ", "), 
  	       call. = FALSE)
  	}
  	
  	# Re-check TMT_Channel (after my_channels() for "126" to become "TMT-126" etc.)
  	# i.e., 16-plex TMT may become 10-plex after Sample_ID removals
  	# but TMT_Channel may include 134 etc. (make sure is compatible with mix-plexes)
  	local({
  	  ls_channels <- unique(label_scheme_full$TMT_Channel)
  	  wrongs <- ls_channels %>% .[! . %in% TMT_levels]
  	  
  	  if (!purrr::is_empty(wrongs)) {
  	    stop("Channels not belong to TMT-", TMT_plex, ":\n", 
  	         purrr::reduce(wrongs, paste, sep = ", "),
  	         call. = FALSE)
  	  }	  
  	})
  	
  	wb <- openxlsx::loadWorkbook(file.path(dat_dir, filename))
  	openxlsx::writeData(wb, sheet = "Setup", label_scheme_full)
  	openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
	} else {
	  label_scheme_full <- label_scheme_full %>% 
	    dplyr::filter(!is.na(Sample_ID)) 
	  
	  # check the uniqueness of Sample_ID for LFQ
	  local({
	    dups <- label_scheme_full %>% 
	      tidyr::unite(key, LCMS_Injection, Sample_ID, remove = FALSE) %>% 
	      dplyr::filter(duplicated(key))
	    
	    if (nrow(dups)) {
	      dups %>%
	        dplyr::select(LCMS_Injection, Sample_ID) %>% 
	        print()
	      
	      stop("Non-unique RAW_File in the above combination of ", 
	           "LCMS_Injection and Sample_ID.", 
	           call. = FALSE)
	    }
	  })
	  
	  label_scheme_full <- local({
  	  label_scheme_full <- label_scheme_full%>% 
  	    dplyr::mutate(Sample_ID. = Sample_ID) %>% 
  	    tidyr::complete(LCMS_Injection, Sample_ID) %>% 
  	    dplyr::mutate(Sample_ID = Sample_ID., 
  	                  Sample_ID = replace_na_smpls(Sample_ID, "Empty")) %>% 
  	    dplyr::select(-Sample_ID.) %>% 
  	    dplyr::arrange(LCMS_Injection) %>% 
  	    dplyr::group_by(LCMS_Injection) %>% 
  	    dplyr::mutate(TMT_Set = row_number()) %>% 
  	    dplyr::ungroup()
  	  
  	  if (!all(is.na(label_scheme_full$RAW_File))) {
  	    label_scheme_full <- label_scheme_full %>% 
  	      dplyr::filter(!is.na(RAW_File))
  	    
  	    # check the unique of Raw_file for LFQ
  	    dup_raws <- label_scheme_full %>% 
  	      dplyr::filter(duplicated(RAW_File)) %>% 
  	      dplyr::select(RAW_File) %>% 
  	      unlist()
  	    
  	    if (!purrr::is_empty(dup_raws)) {
  	      stop("Duplicaed entries under `RAW_File`:\n", 
  	           purrr::reduce(dup_raws, paste, sep = ", "), 
  	           call. = FALSE)
  	    }
  	  }
  	  
  	  return(label_scheme_full)
	  })

	  label_scheme_full <- dplyr::bind_cols(
	    label_scheme_full %>% dplyr::select(Sample_ID), 
	    label_scheme_full %>% dplyr::select(-Sample_ID))
	  
	  wb <- openxlsx::loadWorkbook(file.path(dat_dir, filename))
	  openxlsx::renameWorksheet(wb, "Setup", "Setup_old")
	  openxlsx::addWorksheet(wb, sheetName = "Setup")
	  openxlsx::removeWorksheet(wb, "Setup_old")
	  openxlsx::writeData(wb, sheet = "Setup", label_scheme_full)
	  openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
	}
	
	save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.rda"))
	simple_label_scheme(dat_dir, label_scheme_full)
}


#' Loads the information of analyte prefractionation
#'
#' @inheritParams load_expts
#' @inheritParams prnHist
#' @import dplyr purrr tidyr openxlsx
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom readxl read_excel
prep_fraction_scheme <- function(dat_dir = NULL, filename = "frac_smry.xlsx") 
{
  old_opts <- options()
  on.exit(options(old_opts), add = TRUE)
  options(warning.length = 5000L)
  
  tbl_mascot <- c(
    "\nMascot (similarly for MaxQuant and Spectrum Mill):\n",
    "-----------------------------------------------------------------------\n", 
    "     TMT_Set    | LCMS_Injection  |     RAW_File     |     PSM_File    \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "        1       |        1        | TMT1_Inj1_F1.raw |    F000001.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "                |                 | TMT1_Inj1_F2.raw |    F000001.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "       ...      |       ...       |        ...       |        ...      \n",
    "----------------|-----------------|------------------|-----------------\n",
    "        2       |        1        | TMT2_Inj1_F1.raw |    F000002.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n",
    "                |                 | TMT2_Inj1_F2.raw |    F000002.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "       ...      |       ...       |        ...       |        ...      \n",
    "-----------------------------------------------------------------------\n"
  )
  
  tbl_mq <- c(
    "\nExamplary MaxQuant:\n",
    "-----------------------------------------------------------------------\n", 
    "     TMT_Set    | LCMS_Injection  |     RAW_File     |     PSM_File    \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "        1       |        1        | dup_msfile_1.raw |   msms_xxx.txt  \n", 
    "----------------|-----------------|------------------|-----------------\n",
    "        2       |        1        | dup_msfile_1.raw |   msms_yyy.txt  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "      ...       |       ...       |        ...       |        ...      \n",
    "-----------------------------------------------------------------------\n"
  )
  
  tbl_lfq <- c(
    "\nExamplary LFQ:\n",
    "-------------------------------------------------------------------------------\n", 
    "    Sample_ID   |    TMT_Set    | LCMS_Injection |    RAW_File    |  Fraction  \n", 
    "----------------|---------------|----------------|----------------|------------\n", 
    "       s1       |  (automated)  |       1        |  msfile_1.raw  |      1     \n", 
    "----------------|---------------|----------------|----------------|------------\n", 
    "       s1       |               |       1        |  msfile_2.raw  |      2     \n", 
    "----------------|---------------|----------------|----------------|------------\n",
    "       ...      |               |      ...       |       ...      |     ...    \n",
    "-------------------------------------------------------------------------------\n"
  )
  
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  load(file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex2()
  TMT_levels <- TMT_levels(TMT_plex)

	fn_suffix <- gsub("^.*\\.([^\\.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^\\.]*$", "", filename)

	if (file.exists(file.path(dat_dir, filename))) {
		if (fn_suffix %in% c("xls", "xlsx")) {
			fraction_scheme <- readxl::read_excel(file.path(dat_dir, filename), 
			                                      sheet = "Fractions") %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else if (fn_suffix == "csv") {
			fraction_scheme <- read.csv(file.path(dat_dir, filename), 
			                            check.names = TRUE, 
			                            header = TRUE, comment.char = "#", 
			                            na.strings = c("", "NA"))  %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else {
			stop(filename, " needs to be in a file format of '.xls' or '.xlsx'.", 
			     call. = FALSE)
		}
	  
	  if (all(is.na(fraction_scheme$LCMS_Injection))) {
	    fraction_scheme$LCMS_Injection <- 1
	  }
	  
	  fraction_scheme <- fraction_scheme %>% 
	    dplyr::filter(!is.na(RAW_File)) %>% 
	    dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File, ignore.case = TRUE)) %>% 
	    dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File, ignore.case = TRUE))
	  
	  if (any(duplicated(fraction_scheme$RAW_File)) && 
	      is.null(fraction_scheme[["PSM_File"]])) {
	    stop("\nDuplicated `RAW_File` names in `", filename, "`:\n", 
	         "This can occur with one RAW file corresponding to ", 
	         "multiple PSM files \n", 
	         "(e.g. searching the same RAW file with different parameter sets).\n", 
	         "To distinguish, add PSM file names to column `PSM_File`:\n", 
	         tbl_mascot, 
	         call. = FALSE)
	  }	
	  
	  if (!is.null(fraction_scheme[["PSM_File"]])) {
	    fraction_scheme <- fraction_scheme %>% tidyr::fill(PSM_File) 

	    local({
	      meta_psm_files <- unique(fraction_scheme$PSM_File) 
	      
	      if (!is.null(meta_psm_files)) {
	        not_oks <- 
	          purrr::map(meta_psm_files, 
	                     ~ list.files(dat_dir, pattern = 
	                                    paste0(.x, "\\.csv$|\\.txt$|\\.ssv$|\\.tsv$"))) %>% 
	          purrr::map_lgl(purrr::is_empty)
	        
	        if (any(not_oks)) {
	          stop("PSM file (names) in `", filename, "` not found under ", dat_dir, ":\n", 
	               purrr::reduce(meta_psm_files[not_oks], paste, sep = "\n"), 
	               call. = FALSE)
	        }
	      }
	    })

	    fraction_scheme <- fraction_scheme %>% 
	      # make explicit file extension to avoid ".dot.csv" etc.
	      dplyr::mutate(PSM_File = gsub("\\.csv$|\\.txt$|\\.ssv$|\\.tsv$", "", PSM_File))

	    local({
	      raw_psm <- fraction_scheme %>% 
	        tidyr::unite(RAW_File, RAW_File, PSM_File, sep = "@") 
	      
	      if (any(duplicated(raw_psm$RAW_File))) {
	        stop("The combination of `RAW_File` and `PSM_File` is not unique in `", 
	             filename, "`.",
	             call. = FALSE)
	      }
	    })
	  } 

	  local({
	    expt_raws <- label_scheme_full %>% 
	      dplyr::select(RAW_File) %>% 
	      filter(!duplicated(RAW_File)) %>% 
	      dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
	      dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File)) %>% 
	      unlist()
	    
	    frac_raws <- fraction_scheme %>% 
	      dplyr::select(RAW_File) %>% 
	      filter(!duplicated(RAW_File)) %>% 
	      dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
	      dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File)) %>% 
	      unlist()
	    
	    # no prefractionation
	    if (!all(is.na(expt_raws))) {
	      not_oks <- frac_raws %>% .[! . %in% expt_raws]
	      
	      if (!purrr::is_empty(not_oks)) {
	        stop(
	          "Remove superfluous file(s) in `", filename, 
	          "` not in `expt_smry.xlsx`:\n", 
	          purrr::reduce(not_oks, paste, sep = ", "), 
	          call. = FALSE
	        )
	      }
	      
	      not_oks_2 <- expt_raws %>% .[! . %in% frac_raws]
	      
	      if (!purrr::is_empty(not_oks_2)) {
	        stop("File(s) in `expt_smry` not in `", filename, "`:\n", 
	             purrr::reduce(not_oks_2, paste, sep = ", "), 
	             call. = FALSE)
	      }
	      
	    }
	  })

	  if (!(is.null(TMT_levels) || rlang::is_integerish(fraction_scheme$TMT_Set))) {
	    stop("Values under `frac_smry.xlsx::TMT_Set` need to be integers.", 
	         call. = FALSE)
	  }
	  
	  if (!is.null(TMT_levels)) {
	    stopifnot("TMT_Set" %in% names(fraction_scheme))
	    
	    fraction_scheme <- fraction_scheme %>% 
	      tidyr::fill(TMT_Set, LCMS_Injection) %>% 
	      dplyr::group_by(TMT_Set, LCMS_Injection) %>% 
	      dplyr::mutate(Fraction = row_number()) %>% 
	      dplyr::arrange(TMT_Set, LCMS_Injection, Fraction)
	  } else {
	    # Although one-to-one between `Sample_ID` and `TMT_Set` in LFQ, 
	    # the use of `TMT_Set` may be counterintuitive and 
	    # users may leave `TMT_Set` column blank. Therefore, column `Sample_ID` enforced.
	    if (!"Sample_ID" %in% names(fraction_scheme)) {
	      stop("Need column `Sample_ID` in `frac_smry.xlsx` for LFQ.\n", tbl_lfq, 
	           call. = FALSE)
	    }
	    
	    fraction_scheme <- fraction_scheme %>% 
	      tidyr::fill(Sample_ID, LCMS_Injection) %>% 
	      dplyr::group_by(Sample_ID, LCMS_Injection) %>% 
	      dplyr::mutate(Fraction = row_number()) %>% 
	      dplyr::arrange(Sample_ID, LCMS_Injection, Fraction) %>% 
	      dplyr::ungroup()
	    
	    # compare Sample_ID(s) before joining
	    local({
	      smpls_fracsmry <- unique(fraction_scheme$Sample_ID)
	      smpls_exptsmry <- unique(label_scheme_full$Sample_ID)
	      
	      missings <- smpls_fracsmry %>% 
	        .[! . %in% smpls_exptsmry]
	      
	      if (!purrr::is_empty(missings)) {
	        stop("Sample IDs in `", filename, "` not found in `expt_smry.xlsx`: \n", 
	             purrr::reduce(missings, paste, sep = ", "), 
	             call. = FALSE)
	      }
	      
	      missings2 <- smpls_exptsmry %>% 
	        .[! . %in% smpls_fracsmry]
	      
	      if (!purrr::is_empty(missings2)) {
	        stop("Sample IDs in `expt_smry.xlsx` not found in `", filename, "`: \n", 
	             purrr::reduce(missings2, paste, sep = ", "), 
	             call. = FALSE)
	      }
	    })
	    
	    fraction_scheme <- fraction_scheme %>% 
	      dplyr::select(-TMT_Set) %>% 
	      dplyr::left_join(label_scheme_full %>% 
	                         dplyr::select(Sample_ID, TMT_Set) %>% 
	                         dplyr::filter(!duplicated(Sample_ID)), 
	                       by = "Sample_ID") %>% 
	      dplyr::select(names(fraction_scheme))
	  }
	  
	  if (!rlang::is_integerish(fraction_scheme$LCMS_Injection)) {
	    stop("Values under `frac_smry.xlsx::LCMS_Injection` need to be integers.", 
	         call. = FALSE)
	  }
	  
	  if (!rlang::is_integerish(fraction_scheme$Fraction)) {
	    stop("Values under `frac_smry.xlsx::Fraction` need to be integers.", 
	         call. = FALSE)
	  }
	  
	  wb <- openxlsx::createWorkbook()
	  openxlsx::addWorksheet(wb, sheetName = "Fractions")
	  openxlsx::writeData(wb, sheet = "Fractions", fraction_scheme)
	  openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
 	} else {
 	  # warning: data in a auto-generated `frac_smry.xlsx` will be incorrect 
 	  #   if they were based on wrong information from `expt_smry.xlsx`

 	  # in case forget to enter RAW_File names
 	  if (anyNA(label_scheme_full$RAW_File)) 
 	    stop("File name(s) of RAW MS data not found under the column ", 
 	         "`RAW_File` in `expt_smry.xlsx`.", 
 	         call. = FALSE)

		# the same RAW file can go into different searches
 	  # e.g. the same RAW but different TMT_Set
 	  if (!is.null(TMT_levels)) {
 	    fraction_scheme <- label_scheme_full %>%
 	      # tidyr::unite(tmt_lcms_raw, c("TMT_Set", "LCMS_Injection", "RAW_File"), remove = FALSE) %>% 
 	      tidyr::unite(tmt_lcms, c("TMT_Set", "LCMS_Injection"), remove = FALSE) %>% 
 	      dplyr::filter(!duplicated(tmt_lcms)) %>%
 	      dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
 	      dplyr::group_by(TMT_Set, LCMS_Injection) %>%
 	      dplyr::mutate(Fraction = row_number())
 	  } else {
 	    fraction_scheme <- label_scheme_full %>%
 	      tidyr::unite(tmt_lcms, c("TMT_Set", "LCMS_Injection"), remove = FALSE) %>% 
 	      dplyr::filter(!duplicated(tmt_lcms)) %>%
 	      dplyr::select(Sample_ID, TMT_Set, LCMS_Injection, RAW_File) %>%
 	      dplyr::group_by(TMT_Set, LCMS_Injection) %>%
 	      dplyr::mutate(Fraction = row_number())
 	  }
		
		if (any(duplicated(fraction_scheme$RAW_File)) && 
		    is.null(fraction_scheme[["PSM_File"]])) {
		  stop("\nDuplicated `RAW_File` names detected during the auto-generation of `", 
		       filename, "`.\n",
		       "This could occur with one RAW file corresponding to multiple PSM files \n", 
		       "(e.g. searching the same RAW file with different parameter sets).\n", 
		       
		       "To distinguish, add PSM file names to column `PSM_File`:\n", 
		       tbl_mascot, 
		       call. = FALSE)
		}	
		
		wb <- openxlsx::createWorkbook()
		openxlsx::addWorksheet(wb, sheetName = "Fractions")
		openxlsx::writeData(wb, sheet = "Fractions", fraction_scheme)
		openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
 	}
	
	# LFQ: Sample_ID in .xlsx, but not .rda
	# TMT: no Sample_ID in either
	fraction_scheme <- 
	  suppressWarnings(dplyr::select(fraction_scheme, -one_of("Sample_ID")))

	save(fraction_scheme, file = file.path(dat_dir, "fraction_scheme.rda"))
}


#'Loads species-specific Databases
#'
#'A function loads a set of precompiled gene sets of 
#'\href{http://current.geneontology.org/products/pages/downloads.html}{GO}
#'and
#'\href{http://software.broadinstitute.org/gsea/msigdb}{molecular signatures}.
#'@seealso \code{\link{load_expts}} for supported species.
#'
#' @examples
#' \donttest{load_dbs("go_sets", "human")}
#'
#'@param species Character string; the name of a species. 
#'@inheritParams prnGSPA
#'@import dplyr
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
load_dbs <- function (gset_nms = NULL, species = NULL) 
{
  if (is.null(gset_nms)) 
    stop("`gset_nms` cannot be NULL.", call. = FALSE)
  if (is.null(species)) 
    stop("`species` cannot be NULL.", call. = FALSE)
  
  defaults <- c("go_sets", "kegg_sets", "c2_msig", "kinsub")
  sys_defs <- gset_nms %>% .[. %in% defaults]
  not_sys_defs <- gset_nms %>% .[! . %in% defaults]

  if (!purrr::is_empty(sys_defs)) {
    abbr_sp <- purrr::map_chr(species, sp_lookup)
    
    if (purrr::is_empty(abbr_sp)) {
      warning("Unknown species.", call. = FALSE)
    }
    
    filelist <- purrr::map(abbr_sp, ~ paste0(sys_defs, "_", .x)) %>% 
      unlist()

    suppressWarnings(data(package = "proteoQ", list = filelist, 
                          envir = environment()))
    gsets <- purrr::map(filelist, ~ try(get(.x), silent = TRUE)) %>% 
      do.call(`c`, .)
    # suppressWarnings(rm(list = filelist, envir = .GlobalEnv))
    
    if (length(gsets) > 0) names(gsets) <- gsub("/", "-", names(gsets))      
  } else {
    gsets <- NULL
  }
  
  if (!purrr::is_empty(not_sys_defs)) {
    if (!all(grepl("\\.rds$", not_sys_defs))) {
      stop("Custom gene set files indicated by `gset_nms` must ", 
           "end with extension `.rds`.", 
           call. = FALSE)
    }

    not_oks <- not_sys_defs %>% .[!file.exists(not_sys_defs)]
    if (!purrr::is_empty(not_oks)) {
      stop("File not found: \n", purrr::reduce(not_oks, paste, sep = ", \n"), 
           call. = FALSE)
    }
    
    gsets2 <- purrr::map(not_sys_defs, readRDS) %>% do.call(`c`, .)
    
    if (length(gsets2))  {
      names(gsets2) <- gsub("/", "-", names(gsets2))
    } else {
      stop("Empty data file in: \n", purrr::reduce(not_sys_defs, paste, sep = ", \n"), 
           call. = FALSE)
    }
  } else {
    gsets2 <- NULL
  }
  
  gsets <- c(gsets, gsets2) %>% .[!duplicated(.)]
  
  if (!length(gsets)) {
    warning("Zero entries in `gsets`.", call. = FALSE)
  }

  assign("gsets", gsets, envir = .GlobalEnv)
} 


#'Load TMT or LFQ experiments
#'
#'\code{load_expts} processes \code{.xlsx} or \code{.csv} files containing the
#'metadata of TMT or LFQ experiments. For simplicity, \code{.xlsx} will be
#'assumed in the document.
#'
#'@section \code{expt_smry.xlsx}: The \code{expt_smry.xlsx} files should be
#'  placed immediately under the file folder defined by \code{dat_dir}. The tab
#'  containing the metadata of TMT or LFQ experiments should be named
#'  \code{Setup}. The \code{Excel} spread sheet therein is comprised of three
#'  tiers of fields: (1) essential, (2) optional default and (3) optional open.
#'  The \code{essential} columns contain the mandatory information of the
#'  experiments. The \code{optional default} columns serve as the fields for
#'  default lookups in sample selection, grouping, ordering, aesthetics, etc.
#'  The \code{optional open} fields allow users to define their own analysis,
#'  aesthetics, etc.
#'
#'  \tabular{ll}{ \strong{Essential column}   \tab \strong{Descrption}\cr
#'  Sample_ID \tab Unique sample IDs \cr TMT_Channel \tab TMT channel names:
#'  \code{126}, \code{127N}, \code{127C} etc. (left void for LFQ) \cr TMT_Set
#'  \tab TMT experiment indexes 1, 2, 3, ... (auto-filled for LFQ) \cr
#'  LCMS_Injection   \tab LC/MS injection indexes 1, 2, 3, ... under a
#'  \code{TMT_Set} \cr RAW_File \tab MS data file names originated by \code{MS}
#'  software(s) \cr Reference \tab Labels indicating reference samples in TMT or
#'  LFQ experiments \cr }
#'
#'  \code{Sample_ID}: values should be unique for entries at a unique
#'  combination of \code{TMT_Channel} and \code{TMT_Set}, or voided for unused
#'  entries. Samples with the same indexes of \code{TMT_Channel} and
#'  \code{TMT_Set} but different indexes of \code{LCMS_Injection} should have
#'  the same value in \code{Sample_ID}. No white space or special characters are
#'  allowed. See also posts for
#'  \href{https://proteoq.netlify.app/post/sample-exlusion-from-metadata}{sample
#'  exclusion}.
#'
#'  \code{RAW_File}: (a) for analysis with off-line fractionation of peptides
#'  before LC/MS, values under the \code{RAW_File} column should be left void.
#'  Instead, the correspondence between the fraction numbers and \code{RAW_File}
#'  names should be specified in a separate file, for example,
#'  \code{frac_smry.xlsx}. (2) For analysis without off-line fractionation, it
#'  is recommended as well to leave the field under the \code{RAW_File} column
#'  blank and instead enter the MS file names in \code{frac_smry.xlsx}.
#'
#'  The set of \code{RAW_File} names in metadata needs to be identifiable in PSM
#'  data. Impalpable mismatches might occur when \code{OS} file names were
#'  altered by MS users and thus different to those recorded internally in MS
#'  data for parsing by search engine(s). In the case, machine-generated MS file
#'  names should be used. In addition, MS files may occasionally have no
#'  contributions to PSM findings. In the case, users will be prompted to remove
#'  these MS file names.
#'
#'  Utilities \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} may
#'  aid matching MS file names between metadata and PSM data. Utility
#'  \code{\link{extract_raws}} extracts the names of MS files under a file
#'  folder. Utility \code{\link{extract_psm_raws}} extracts the names of MS
#'  files that are available in PSM data.
#'
#'  \code{Reference}: reference entry(entries) are indicated with non-void
#'  string(s).
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
#'  \cr \cr}
#'
#'  \tabular{ll}{ \strong{Exemplary optional open column}   \tab
#'  \strong{Descrption}\cr Term \tab Categorical terms for statistical modeling.
#'  \cr Peptide_Yield \tab Yields of peptides in sample handling \cr}
#'
#'
#'@section \code{frac_smry.xlsx}: \tabular{ll}{ \strong{Column}   \tab
#'  \strong{Descrption}\cr Sample_ID \tab Unique sample IDs (only required with
#'  LFQ) \cr TMT_Set \tab TMT experiment indexes (auto-filled for LFQ) \cr
#'  LCMS_Injection \tab LC/MS injection indexes \cr Fraction \tab Fraction
#'  indexes under a \code{TMT_Set} \cr RAW_File \tab MS data file names \cr
#'  PSM_File \tab Names of PSM files. Required only when one \code{RAW_File} can
#'  be linked to multiple PSM files (e.g. F012345.csv and F012346.csv both from
#'  ms_1.raw). }
#'
#'@family normalization functions
#'@seealso \emph{Data normalization} \cr \code{\link{normPSM}} for extended
#'  examples in PSM data normalization \cr \code{\link{PSM2Pep}} for extended
#'  examples in PSM to peptide summarization \cr \code{\link{mergePep}} for
#'  extended examples in peptide data merging \cr \code{\link{standPep}} for
#'  extended examples in peptide data normalization \cr \code{\link{Pep2Prn}}
#'  for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization.
#'  \cr \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples
#'  in data purging \cr \code{\link{pepHist}} and \code{\link{prnHist}} for
#'  extended examples in histogram visualization. \cr \code{\link{extract_raws}}
#'  and \code{\link{extract_psm_raws}} for extracting MS file names \cr
#'
#'@family data row filtration
#'@seealso \emph{User-friendly utilities for variable arguments of `filter_...`}
#'  \cr \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr
#'
#'@family missing value imputation
#'@seealso \emph{Missing values} \cr \code{\link{pepImp}} and
#'  \code{\link{prnImp}} for missing value imputation \cr
#'
#'@family basic informatics
#'@seealso \emph{Informatics} \cr \code{\link{pepSig}} and \code{\link{prnSig}}
#'  for significance tests \cr \code{\link{pepVol}} and \code{\link{prnVol}} for
#'  volcano plot visualization \cr \code{\link{prnGSPA}} for gene set enrichment
#'  analysis by protein significance pVals \cr \code{\link{gspaMap}} for mapping
#'  GSPA to volcano plot visualization \cr \code{\link{prnGSPAHM}} for heat map
#'  and network visualization of GSPA results \cr \code{\link{prnGSVA}} for gene
#'  set variance analysis \cr \code{\link{prnGSEA}} for data preparation for
#'  online GSEA. \cr \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS
#'  visualization \cr \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA
#'  visualization \cr \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA
#'  visualization \cr \code{\link{pepHM}} and \code{\link{prnHM}} for heat map
#'  visualization \cr \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'  \code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for
#'  correlation plots \cr \code{\link{anal_prnTrend}} and
#'  \code{\link{plot_prnTrend}} for trend analysis and visualization \cr
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}},
#'  \code{\link{plot_pepNMFCon}}, \code{\link{plot_prnNMFCon}},
#'  \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr
#'
#'  \emph{Custom databases} \cr \code{\link{Uni2Entrez}} for lookups between
#'  UniProt accessions and Entrez IDs \cr \code{\link{Ref2Entrez}} for lookups
#'  among RefSeq accessions, gene names and Entrez IDs \cr \code{\link{prepGO}}
#'  for
#'  \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'   ontology}} \cr \code{\link{prepMSig}} for
#'  \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular
#'   signatures} \cr \code{\link{prepString}} and \code{\link{anal_prnString}}
#'  for STRING-DB \cr
#'
#'  \emph{Workflow scripts} \cr # TMT \cr system.file("extdata",
#'  "workflow_tmt_base.R", package = "proteoQ") \cr system.file("extdata",
#'  "workflow_tmt_ext.R", package = "proteoQ") \cr
#'
#'  # LFQ \cr system.file("extdata", "workflow_lfq_base.R", package = "proteoQ")
#'  \cr
#'
#'  \emph{Metadata files} \cr # TMT, no fractionation --- OK without
#'  `frac_smry.xlsx` \cr # (a. no references) \cr system.file("extdata",
#'  "expt_smry_no_prefrac.xlsx", package = "proteoQDA") \cr # (b. W2 and W16
#'  references) \cr system.file("extdata",
#'  "expt_smry_no_prefrac_ref_w2_w16.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation \cr # (a. no references) \cr system.file("extdata",
#'  "expt_smry_gtmt.xlsx", package = "proteoQDA") \cr system.file("extdata",
#'  "frac_smry_gtmt.xlsx", package = "proteoQDA") \cr
#'
#'  # (b. W2 references) \cr system.file("extdata", "expt_smry_ref_w2.xlsx",
#'  package = "proteoQDA") \cr system.file("extdata", "frac_smry_gtmt.xlsx",
#'  package = "proteoQDA") \cr
#'
#'  # (c. W2 and W16 references) \cr system.file("extdata",
#'  "expt_smry_ref_w2_w16.xlsx", package = "proteoQDA") \cr
#'  system.file("extdata", "frac_smry_gtmt.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation (global + phospho) \cr system.file("extdata",
#'  "expt_smry_tmt_cmbn.xlsx", package = "proteoQDA") \cr system.file("extdata",
#'  "frac_smry_tmt_cmbn.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation, one MS to multiple PSM files \cr
#'  system.file("extdata", "expt_smry_psmfiles.xlsx", package = "proteoQDA") \cr
#'  system.file("extdata", "frac_smry_psmfiles.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation, mixed-plexes \cr # (column PSM_File needed; as with
#'  this example, \cr #  mixed-plexes results are actually from the same MS
#'  files \cr #  but searched separately at 6- and 10-plex settings!) \cr
#'  system.file("extdata", "expt_smry_mixplexes.xlsx", package = "proteoQDA")
#'  \cr system.file("extdata", "frac_smry_mixplexes.xlsx", package =
#'  "proteoQDA") \cr
#'
#'  # LFQ, prefractionation \cr system.file("extdata", "expt_smry_plfq.xlsx",
#'  package = "proteoQDA") \cr system.file("extdata", "frac_smry_plfq.xlsx",
#'  package = "proteoQDA") \cr
#'
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@family Monoisotopic masses of peptides
#'@seealso \emph{MS1 peptide masses} \cr \code{\link[proteoM]{calc_pepmasses}}
#'  for mono-isotopic masses of peptides from fasta databases \cr
#'  \code{\link[proteoM]{calc_monopeptide}} for mono-isotopic masses of peptides
#'  from individual sequences \cr \code{\link[proteoM]{parse_unimod}} for
#'  parsing \href{https://www.unimod.org/}{Unimod} fixed modifications, variable
#'  modifications and neutral losses. \cr \code{\link[proteoM]{find_unimod}} for
#'  finding a Unimod
#'
#'
#'
#'@param dat_dir A character string to the working directory. The default is to
#'  match the value under the global environment.
#'@param expt_smry A character string to a \code{.xlsx} file containing the
#'  metadata of TMT or LFQ experiments. The default is \code{expt_smry.xlsx}.
#'@param frac_smry A character string to a \code{.xlsx} file containing peptide
#'  fractionation summary. The default is \code{frac_smry.xlsx}.
#'
#'@example inst/extdata/examples/load_expts_.R
#'
#'@import dplyr fs
#'@rawNamespace import(rlang, except = c(list_along, invoke, flatten_raw,
#'  modify, as_function, flatten_dbl, flatten_lgl, flatten_int, flatten_chr,
#'  splice, flatten, prepend, "%@%", ":="))
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
load_expts <- function (dat_dir = NULL, expt_smry = "expt_smry.xlsx", 
                        frac_smry = "frac_smry.xlsx") 
{
  on.exit(mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
            save_call("load_expts"))
  
  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))

  dat_dir <- set_dat_dir(dat_dir)

  prep_label_scheme(dat_dir, expt_smry)
  prep_fraction_scheme(dat_dir, frac_smry)
  
  local({
    session_info <- sessionInfo()
    dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)
    try(save(session_info, file = file.path(dat_dir, "Calls", "proteoQ.rda")))
  })
}


#' Reload the "expt_smry.xlsx" and "frac_smry.xlsx"
#'
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom fs file_info
reload_expts <- function() 
{
  dat_dir <- get_gl_dat_dir()
  
  expt_smry <- match_call_arg(load_expts, expt_smry)
  frac_smry <- match_call_arg(load_expts, frac_smry)
  
  fi_xlsx <- fs::file_info(file.path(dat_dir, expt_smry))$change_time
  
  if (is.na(fi_xlsx)) 
    stop("Time stamp of ", expt_smry, " not available.")
  
  fi_rda <- 
    fs::file_info(file.path(dat_dir, "label_scheme_full.rda"))$change_time
  
  if (fi_xlsx > fi_rda) 
    load_expts(dat_dir = dat_dir, 
               expt_smry = !!expt_smry, 
               frac_smry = !!frac_smry)
}


#' Extracts the channel information in TMT experiments
#'
#' A function returns the indexes of TMT channels that are associated to
#' reference(s), sample(s) and probable unused void(s).
#'
#' @param set_idx Numeric.  The index of a multiplex TMT experiment in metadata
#'   files such as \code{label_scheme.xlsx} and \code{frac_scheme.xlsx}.
#' @param injn_idx Numeric. The index of \code{LCMS_Inj} in metadata files such
#'   as \code{label_scheme.xlsx} and \code{frac_scheme.xlsx}.
#' @inheritParams load_expts
#' @return Three lists of indexes: \code{refChannels}, reference channels(s);
#'   \code{emptyChannels}, empty channel(s) that were not used for sample
#'   labeling; \code{labeledChannels}, non-empty channels including both
#'   reference(s) and sample(s).
#' 
#' @importFrom dplyr select filter
channelInfo <- function (dat_dir = NULL, set_idx = NULL, injn_idx = 1) 
{
	stopifnot(length(set_idx) == 1L)
  
  if (is.null(set_idx)) 
    stop("Need to specify `set_idx`.", call. = FALSE)
  if (is.null(dat_dir)) 
    stop("Need to specify `dat_dir`.", call. = FALSE)

  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  
  label_scheme_sub <- label_scheme_full %>%
    dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx)
  
  if (nrow(label_scheme_sub) == 0) 
    stop("No samples at TMT_Set ", set_idx, " and LCMS_Injection ", injn_idx, ".", 
         call. = FALSE)

	ref <- label_scheme_sub$Reference

	empty_channels <- is.na(label_scheme_sub$Sample_ID) |
	  grepl("^Empty\\.[0-9]+$", label_scheme_sub$Sample_ID, ignore.case = TRUE)

	labeled_channels <- !empty_channels

	out <- list(
		refChannels = ref,
		emptyChannels = empty_channels,
		labeledChannels = labeled_channels
	)

	lapply(out, which)
}


#' Finds the number of multiplex TMT experiments
#'
#' @param label_scheme_full The label_scheme with the probable inclusion of
#'   different LCMS_inj under the same TMT_Set.
n_TMT_sets <- function (label_scheme_full) 
  length(unique(label_scheme_full$TMT_Set))


#' Finds the maximum number of LCMS injections under the same TMT_Set.
#'
#' @inheritParams n_TMT_sets
n_LCMS <- function (label_scheme_full) 
{
  label_scheme_full %>% 
    dplyr::select(TMT_Set, LCMS_Injection) %>% 
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = FALSE) %>% 
    dplyr::filter(!duplicated(TMT_inj)) %>% 
    dplyr::group_by(TMT_Set) %>% 
    dplyr::mutate(n_LCMS = n()) %>% 
    dplyr::select(TMT_Set, n_LCMS) %>% 
    dplyr::filter(!duplicated(TMT_Set)) %>% 
    dplyr::ungroup() 
}


#' Finds the multiplexity of TMT labels
#'
#' \code{TMT_plex} returns the multiplexity of TMT labels.
#' @inheritParams n_TMT_sets
TMT_plex <- function (label_scheme_full) 
  nlevels(as.factor(label_scheme_full$TMT_Channel))

#' Finds the multiplexity of TMT labels from previously made
#' label_scheme_full.rda
#'
#' \code{TMT_plex} returns the multiplexity of TMT labels.
TMT_plex2 <- function () 
{
  load(file.path(get_gl_dat_dir(), "label_scheme_full.rda"))
  TMT_plex(label_scheme_full)
}


#' Finds the factor levels of TMT labels
#'
#' \code{TMT_levels} returns the factor levels of TMT labels.
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
TMT_levels <- function (TMT_plex) 
{
  TMT_levels <- if (TMT_plex == 18) 
    c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
                    "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
                    "TMT-131N", "TMT-131C", "TMT-132N", "TMT-132C", 
                    "TMT-133N", "TMT-133C", "TMT-134N", "TMT-134C", "TMT-135N")
  else if (TMT_plex == 16) 
	  c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
	                  "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
	                  "TMT-131N", "TMT-131C", "TMT-132N", "TMT-132C", 
	                  "TMT-133N", "TMT-133C", "TMT-134N")
	else if (TMT_plex == 11) 
	  c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
	                  "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
	                  "TMT-131N", "TMT-131C")
	else if (TMT_plex == 10) 
	  c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
	                  "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", "TMT-131")
	else if (TMT_plex == 6) 
	  c("TMT-126", "TMT-127", "TMT-128", "TMT-129", "TMT-130", "TMT-131")
	else if (TMT_plex == 1) 
	  c("TMT-126")
	else if (TMT_plex == 0) 
	  NULL
}


#' Simplifies label schemes from \code{label_scheme_full}
#'
#' Removes duplicated sample entries under different LC/MS injections.
#' @inheritParams load_expts
#' @inheritParams n_TMT_sets
simple_label_scheme <- function (dat_dir, label_scheme_full) 
{
	TMT_plex <- TMT_plex(label_scheme_full)
	TMT_levels <- TMT_levels(TMT_plex)

	# OK Sample_ID of "sample_x" and "Empty.1" at the same TMT_Set and TMT_Channel,
	# but different LCMS_Injection
	label_scheme <- label_scheme_full %>%
	  dplyr::filter(!duplicated(Sample_ID), !is.na(Sample_ID)) %>%
	  dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
	  dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)
	
	multi_dips <- label_scheme %>% 
	  tidyr::unite(key, TMT_Channel, TMT_Set, remove = FALSE) %>% 
	  dplyr::group_by(key) %>% 
	  dplyr::mutate(n = n()) %>% 
	  dplyr::filter(n > 1) %>% 
	  dplyr::ungroup()
	
	label_scheme <- label_scheme %>% 
	  dplyr::filter(! .data$Sample_ID %in% multi_dips$Sample_ID)
	
	label_scheme <- multi_dips %>% 
	  dplyr::filter(!grepl("^Empty\\.[0-9]+$", Sample_ID)) %>% 
	  dplyr::select(names(label_scheme)) %>% 
	  dplyr::bind_rows(label_scheme) %>% 
	  dplyr::mutate(LCMS_Injection = 1) %>% 
	  dplyr::arrange(TMT_Set, TMT_Channel)

	if (nrow(label_scheme) < (TMT_plex * n_TMT_sets(label_scheme))) 
	  stop("Duplicated sample ID(s) in `expt_smry.xlsx`", call. = FALSE)

	save(label_scheme, file = file.path(dat_dir, "label_scheme.rda"))
	
	invisible(label_scheme)
}


#' Find mismatches in RAW file names
#'
#' \code{check_raws} finds mismatched RAW files between expt_smry.xlsx and
#' PSM outputs.
#' @param df A data frame containing the PSM table from database searches.
check_raws <- function(df) 
{
  stopifnot ("RAW_File" %in% names(df))
  
  dat_dir <- get_gl_dat_dir()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))

  ## program-generated frac_smry.xlsx may be based on wrong information from expt_smry.xlsx
  local({
    ls_raws <- label_scheme_full$RAW_File %>% unique()
    fs_raws <- fraction_scheme$RAW_File %>% unique()
    
    if (!(all(is.na(ls_raws)) || all(ls_raws %in% fs_raws))) {
      load(file.path(dat_dir, "Calls", "load_expts.rda"))
      fn_frac <- call_pars$frac_smry
      unlink(file.path(dat_dir, fn_frac))
      prep_fraction_scheme(dat_dir, fn_frac)
      load(file = file.path(dat_dir, "fraction_scheme.rda"))
    }
  })

  local({
    ls_tmt <- unique(label_scheme_full$TMT_Set)
    fs_tmt <- unique(fraction_scheme$TMT_Set)
    extra_fs_tmt <- fs_tmt %>% .[! . %in% ls_tmt]
    extra_ls_tmt <- ls_tmt %>% .[! . %in% fs_tmt]

    if (length(extra_fs_tmt)) 
      stop("TMT_Set ", purrr::reduce(extra_fs_tmt, paste, sep = ", "), 
           " in fraction scheme not found in label scheme.", 
           call. = FALSE)
    
    if (length(extra_ls_tmt)) 
      stop("TMT_Set ", purrr::reduce(extra_ls_tmt, paste, sep = ", "), 
           " in label scheme not found in fraction scheme.", 
           call. = FALSE)
  })
  
  local({
    ls_tmtinj <- label_scheme_full %>% 
      tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
      dplyr::select(TMT_inj) %>% 
      unique() %>% 
      unlist()
    
    fs_tmtinj <- fraction_scheme %>% 
      tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
      dplyr::select(TMT_inj) %>% 
      unique() %>% 
      unlist()
    
    extra_fs_tmt <- fs_tmtinj %>% .[! . %in% ls_tmtinj]
    extra_ls_tmt <- ls_tmtinj %>% .[! . %in% fs_tmtinj]
    
    if (length(extra_fs_tmt)) 
      stop("The combination of TMT_Set & LCMS ", 
           purrr::reduce(extra_fs_tmt, paste, sep = ", "), 
           " in fraction scheme not found in label scheme.", 
           call. = FALSE)
    
    if (length(extra_ls_tmt)) 
      stop("The combination of TMT_Set & LCMS ", 
           purrr::reduce(extra_fs_tmt, paste, sep = ", "), 
           " in label scheme not found in fraction scheme.", 
           call. = FALSE)
  })
  
  tmtinj_raw <- fraction_scheme %>% 
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>%
    dplyr::select(-Fraction) %>%
    dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
		dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File)) # Bruker
  
  ms_raws <- df$RAW_File %>% unique()
  label_scheme_raws <- tmtinj_raw$RAW_File %>% unique()
  
  missing_ms_raws <- ms_raws %>% .[! . %in% label_scheme_raws]
  wrong_label_scheme_raws <- label_scheme_raws %>% .[! . %in% ms_raws]
  
  # --- rm `missing_ms_raws` from `df`
  if (length(missing_ms_raws)) {
    df <- df %>% dplyr::filter(! RAW_File %in% missing_ms_raws)
    
    warning("RAW file (names) not in metadata and ", 
            "corresponding entries removed from PSM data:\n", 
            purrr::reduce(missing_ms_raws, paste, sep = "\n"), 
            call. = FALSE)
  }
  
  if (length(wrong_label_scheme_raws)) 
    stop("\n=========================================================================\n", 
         "RAW file name(s) in metadata have no corresponding entries in PSM data:\n", 
         "(Hint: check the possibility that MS file(s) have ", 
         "no PSM contributions.)\n\n", 
         purrr::reduce(wrong_label_scheme_raws, paste, sep = "\n"), 
         "\n=========================================================================\n", 
         call. = FALSE)
  
  ## RAW_File may not be unique if the same RAW goes into different DAT files (searches)
  # df <- df %>% dplyr::left_join(tmtinj_raw, id = "RAW_File")
  invisible (list(lookup = tmtinj_raw, df = df))
}


#' Find the TMT-plex from Mascot PSM exports
#' @param df A PSM export from Mascot
#' @param pep_seq_rows The row(s) contain character string "pep_seq".
find_masct_tmtplex <- function(df, pep_seq_rows = NULL) 
{
  if (is.null(pep_seq_rows)) 
    pep_seq_rows <- grep("pep_seq", df)
  
  if (purrr::is_empty(pep_seq_rows)) 
    stop("The row of `pep_seq` not found.", call. = FALSE)
  
  first_line <- df[pep_seq_rows[1] + 1]
  
  mascot_tmtplex <- if (grepl("\"135C\"", first_line, fixed = TRUE) || 
      grepl("\"135\"", first_line, fixed = TRUE)) 
    18 
  else if (grepl("\"134N\"", first_line, fixed = TRUE) || 
      grepl("\"134\"", first_line, fixed = TRUE)) 
    16
  else if (grepl("\"131C\"", first_line, fixed = TRUE)) 
    11
  else if (grepl("\"131\"", first_line, fixed = TRUE) && 
             grepl("\"130C\"", first_line, fixed = TRUE)) 
    10
  else if (grepl("\"131\"", first_line, fixed = TRUE)) 
    6
  else 
    0
}

