#' Matches the current id to the id in normPep or normPrn
#'
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_identifier <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene")) {
  if (id %in% c("prot_acc", "gene")) {
    fn_pars <- "normPrn.txt"
  } else if (id %in% c("pep_seq", "pep_seq_mod")) {
    fn_pars <- "normPep.txt"
  }
  
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if (!is.null(dim(call_pars))) {
    id <- call_pars %>%
      dplyr::filter(var == "id") %>%
      dplyr::select("value.1") %>%
      unlist() %>%
      as.character()
  }
  
  return(id)
}


#' Add Z_log2_R, sd_N_log2_R and sd_Z_log2_R
calc_more_psm_sd <- function (df, group_psm_by, range_log2r, range_int, set_idx, injn_idx) {
  # SD columns for "N_log2_R"
  df <- df %>% 
    calcSD_Splex(group_psm_by, "N_log2_R") %>% 
    `names<-`(gsub("^N_log2_R", "sd_N_log2_R", names(.))) %>% 
    dplyr::right_join(df, by = group_psm_by)
  
  # "Z_log2_R" columns and their SD
  cf_SD <- calc_sd_fcts_psm(df, range_log2r, range_int, set_idx, injn_idx)
  
  cf_x <- df %>%
    dplyr::select(matches("^N_log2_R[0-9]{3}")) %>%
    `colnames<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.))) %>%
    dplyr::summarise_all(~ median(.x, na.rm = TRUE)) 
  
  df <- mapply(normSD, df[,grepl("^N_log2_R[0-9]{3}", names(df))],
               center = cf_x, SD = cf_SD$fct, SIMPLIFY = FALSE) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
    cbind(df, .)
  
  df_z <- df %>% dplyr::select(grep("^Z_log2_R", names(.)))
  nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
  df_z[, nan_cols] <- 0
  df[, grep("^Z_log2_R", names(df))] <- df_z
  rm(df_z, nan_cols)
  
  df <- df %>% 
    calcSD_Splex(group_psm_by, "Z_log2_R") %>% 
    `names<-`(gsub("^Z_log2_R", "sd_Z_log2_R", names(.))) %>% 
    dplyr::right_join(df, by = group_psm_by)
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
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


#' Matches the current pep_id to the pep_id in normPSM
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_normPSM_pepid <- function (id = c("pep_seq")) {
  stopifnot(id %in% c("pep_seq", "pep_seq_mod"))
  
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", "normPSM.txt"), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if (!is.null(dim(call_pars))) {
    id <- call_pars %>%
      dplyr::filter(var == "group_psm_by") %>%
      dplyr::select("value.1") %>%
      unlist() %>%
      as.character()
  }
  
  return(id)
}


#' Matches the current prot_id to the prot_id in normPSM
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_normPSM_protid <- function (id = c("prot_acc")) {
  stopifnot(id %in% c("prot_acc", "gene"))
  
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", "normPSM.txt"), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if (!is.null(dim(call_pars))) {
    id <- call_pars %>%
      dplyr::filter(var == "group_pep_by") %>%
      dplyr::select("value.1") %>%
      unlist() %>%
      as.character()
  }
  
  return(id)
}


#' Saves the arguments in a function call
#'
#' @param call_pars Language.
#' @param fn The name of function being saved.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
save_call <- function(call_pars, fn) {
	dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)

	call_pars[names(call_pars) == "..."] <- NULL
	save(call_pars, file = file.path(dat_dir, "Calls", paste0(fn, ".rda")))

	purrr::map(call_pars , as.character) %>%
	  do.call(purrr::quietly(rbind), .) %>% 
	  `[[`(1) %>% 
		data.frame() %>%
		`colnames<-`(paste("value", 1:ncol(.), sep = ".")) %>%
		tibble::rownames_to_column("var") %>%
		write.table(file.path(dat_dir, "Calls", paste0(fn, ".txt")), sep = '\t',
		            col.names = TRUE, row.names = FALSE)
}

#' Find the setting of current par in normPSM
#' @param par
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_normPSM_par <- function (par = "use_lowercase_aa") {
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", "normPSM.txt"), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if (!is.null(dim(call_pars))) {
    res <- call_pars %>%
      dplyr::filter(var == par) %>%
      dplyr::select("value.1") %>%
      unlist() %>%
      as.character()
  }
  
  return(res)
}


#' Matches the file name containing the summary of TMT experiment
#'
#' The default file name is \code{expt_smry.xlsx}. The \code{match_expt} matches
#' the file name provided by users.
#'
#' @param fn_pars.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_expt <- function (fn_pars = "load_expts.txt") {
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
	                               header = TRUE, sep = "\t", comment.char = "#"),
	                      error = function(e) NA)

	if(!is.null(dim(call_pars))) {
		expt_smry <- call_pars %>%
			dplyr::filter(var == "expt_smry") %>%
			dplyr::select("value.1") %>%
			unlist() %>%
			as.character()
	} else {
		expt_smry <- "expt_smry.xlsx"
	}

	return(expt_smry)
}


#' Matches the file name containing the information of analyte fractionation
#'
#' The default file name is \code{frac_smry.xlsx}. The \code{match_frac} matches
#' the file name provided by users.
#'
#' @param fn_pars.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_frac <- function (fn_pars = "load_expts.txt") {
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if(!is.null(dim(call_pars))) {
    frac_smry <- call_pars %>%
      dplyr::filter(var == "frac_smry") %>%
      dplyr::select("value.1") %>%
      unlist() %>%
      as.character()
  } else {
    frac_smry <- "frac_smry.xlsx"
  }
  
  return(frac_smry)
}


#' Matches the name of fasta files to those in normPSM
#'
#' @param fasta.
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
match_fasta <- function () {
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", "normPSM.txt"), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if (!is.null(dim(call_pars))) {
    fasta <- call_pars %>%
      dplyr::filter(var == "fasta") %>%
      dplyr::select(-var) %>%
      unlist() %>%
      as.character()
    
    if (is_empty(fasta)) 
      stop("No fasta file names(s) found from the call to `normPSM`.", call. = FALSE)
  } else {
    stop("No saved parameters found from the call to `normPSM`.", call. = FALSE)
  }
  
  return(fasta)
}


#' Matches the name of GSPA result file
#'
#' @param anal_type
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
match_gspa_filename <- function (anal_type = "GSPA", subdir = NULL) {
  stopifnot(!is.null(subdir))
  
  if (anal_type == "GSPA") {
    call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", "prnGSPA.txt"), check.names = FALSE, 
                                   header = TRUE, sep = "\t", comment.char = "#"), 
                          error = function(e) NA)    
  }

  if (!is.null(dim(call_pars))) {
    filename <- call_pars %>%
      dplyr::filter(var == "filename") %>% 
      dplyr::select(-var) %>%
      unlist() %>%
      as.character() %>% 
      unique()

    if (is_empty(filename)) {
      filename <- list.files(path = file.path(dat_dir, "Protein\\GSPA", subdir), 
                             pattern = "^Protein_GSPA_.*\\.csv$", 
                             full.names = FALSE)
      
      if (length(filename) > 1) stop("More than one result file found under `", subdir, "`", call. = FALSE)
      
      if (is_empty(filename)) stop("No result file found under `", sub_dir, "`", call. = FALSE)
    }
  } else {
    stop("No saved parameters found from the call to `prnGSPA`.", call. = FALSE)
  }
  
  return(filename)
}


#' Matches impute_na
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_sigTest_imputena <- function (id = c("prot_acc")) {
  stopifnot(is_string(id))
  stopifnot(id %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls\\pepSig.txt"), check.names = FALSE,
                                   header = TRUE, sep = "\t", comment.char = "#"),
                          error = function(e) NA)
  } else if (id %in% c("prot_acc", "gene")) {
    call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls\\prnSig.txt"), check.names = FALSE,
                                   header = TRUE, sep = "\t", comment.char = "#"),
                          error = function(e) NA)
  }
  
  if (!is.null(dim(call_pars))) {
    impute_na <- call_pars %>%
      dplyr::filter(var == "impute_na") %>%
      dplyr::select("value.1") %>%
      unlist() %>%
      as.character() %>% 
      as.logical()
  }
  
  return(impute_na)
}

#' Matches the file name containing the summary of TMT experiment
#'
#' The default file name is \code{expt_smry.xlsx}. The \code{match_expt} matches
#' the file name provided by users.
#'
#' @param fn_pars.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_expt <- function (fn_pars = "load_expts.rda") {
  file <- file.path(dat_dir, "Calls\\load_expts.rda")
  stopifnot(file.exists(file))
  
  load(file = file)
  if (is.null(call_pars$expt_smry)) call_pars$expt_smry <- "expt_smry.xlsx"
  call_pars$expt_smry
}


#' Matches the file name containing the information of analyte fractionation
#'
#' The default file name is \code{frac_smry.xlsx}. The \code{match_frac} matches
#' the file name provided by users.
#'
#' @param fn_pars.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_frac <- function (fn_pars = "load_expts.rda") {
  file <- file.path(dat_dir, "Calls\\load_expts.rda")
  stopifnot(file.exists(file))
  
  load(file = file)
  if (is.null(call_pars$frac_smry)) call_pars$frac_smry <- "frac_smry.xlsx"
  call_pars$frac_smry
}

#' Matches the name of fasta files to those in normPSM
#'
#' @param fasta.
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
match_fasta <- function () {
  file <- file.path(dat_dir, "Calls\\normPSM.rda")
  if (!file.exists(file)) stop("Run `normPSM` first.", call. = FALSE)

  load(file = file)
  if (is.null(call_pars$fasta)) 
    stop("No fasta file names(s) found from the call to `normPSM`.", call. = FALSE)
  
  call_pars$fasta
}

#' Find the setting of current par in normPSM
#' @param par
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_normPSM_par <- function (par = "use_lowercase_aa") {
  file <- file.path(dat_dir, "Calls\\normPSM.rda")
  stopifnot(file.exists(file))
  
  load(file = file)
  stopifnot(par %in% names(call_pars))
  call_pars[[par]]
}


#' Matches the current pep_id to the pep_id in normPSM
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_normPSM_pepid <- function (id = c("pep_seq")) {
  file <- file.path(dat_dir, "Calls\\normPSM.rda")
  stopifnot(id %in% c("pep_seq", "pep_seq_mod"), file.exists(file))

  load(file = file)
  stopifnot("group_psm_by" %in% names(call_pars))
  call_pars$group_psm_by
}


#' Matches the current prot_id to the prot_id in normPSM
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_normPSM_protid <- function (id = c("prot_acc")) {
  file <- file.path(dat_dir, "Calls\\normPSM.rda")
  stopifnot(id %in% c("prot_acc", "gene"), file.exists(file))
  
  load(file = file)
  stopifnot("group_pep_by" %in% names(call_pars))
  call_pars$group_pep_by
}





