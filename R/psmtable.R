#' Extract RAW MS file names
#'
#' Extract a list of \code{RAW} file names that can be passed to metadata file
#' of \code{expt_smry.xlsx} and/or \code{frac_smry.xlsx}
#'
#' @param raw_dir A character string to the directory of MS data. 
#' @inheritParams load_expts
#' @examples
#' \dontrun{
#' # Supposed that RAW MS files are stored under "~/my_raw"
#' extract_raws("~/my_raw")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom tools md5sum
#' @export
extract_raws <- function(raw_dir = NULL, dat_dir = NULL) {
  if (is.null(raw_dir)) stop("`raw_dir` cannot be `NULL`.", 
                             call. = FALSE)
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()

  fns <- names(tools::md5sum(dir(raw_dir, pattern = "\\.raw$", 
                                 full.names = FALSE)))
  
  data.frame(Index = seq_along(fns), RAW_File = fns) %T>% 
    write.table(file.path(dat_dir, "raw_list.txt"), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE)
		
	message("List of `.raw` file names saved in ", 
	        file.path(dat_dir, "raw_list.txt"))
}


#' Extract the names of RAW MS files from PSM data
#'
#' \code{extract_psm_raws} extracts the MS file names from the PSM data under
#' the working directory.
#'
#' @param type Character string indicating the type of PSM.
#' @inheritParams load_expts
#' @examples
#' \donttest{
#' extract_psm_raws(mascot)
#' 
#' extract_psm_raws(maxquant)
#' 
#' extract_psm_raws(spectrum_mill)
#' }
#'
#' @import dplyr tidyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @export
extract_psm_raws <- function(type = c("mascot", "maxquant", "spectrum_mill", "msfragger"), 
                             dat_dir = NULL) {
  batchPSMheader_2 <- function(filelist) {
    df <- readLines(file.path(dat_dir, filelist))
    
    pep_seq_rows <- grep("pep_seq", df)
    
    TMT_plex <- find_masct_tmtplex(df, pep_seq_rows)
    
    unassign_hits_row <- grep("Peptide matches not assigned to protein hits", df)
    if (! purrr::is_empty(unassign_hits_row)) {
      psm_end_row <- unassign_hits_row - 2
    } else if (length(pep_seq_rows) > 1) {
      psm_end_row <- pep_seq_rows[2] - 4
    } else {
      psm_end_row <- length(df)
    }
    
    df <- df[pep_seq_rows[1] : psm_end_row]
    df <- gsub("\"---\"", -1, df, fixed = TRUE)
    df <- gsub("\"###\"", -1, df, fixed = TRUE)
    
    df[1] <- paste0(df[1], paste(rep(",", TMT_plex * 4 -2), collapse = ''))
    
    output_prefix <- gsub("\\.csv$", "", filelist)
    writeLines(df, file.path(dat_dir, "PSM/temp", 
                             paste0(output_prefix, "_hdr_rm.csv")))
    
    return(TMT_plex)
  }
    
  find_mascot_psmraws <-function(filelist) {
    dir.create(file.path(dat_dir, "PSM/temp"), 
               recursive = TRUE, showWarnings = FALSE)
    
    TMT_plex <- purrr::map(filelist, batchPSMheader_2)

    purrr::walk2(gsub("\\.csv$", "_hdr_rm.csv", filelist), TMT_plex, ~ {
      df <- read.delim(file.path(dat_dir, "PSM/temp", .x), 
                       sep = ',', check.names = FALSE, 
                 header = TRUE, stringsAsFactors = FALSE, 
                 quote = "\"",fill = TRUE , skip = 0)
      
      TMT_plex <- .y
      
      r_start <- which(names(df) == "pep_scan_title") + 1
      int_end <- ncol(df)
      if(int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
      
      if (TMT_plex == 16) {
        col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                       "R130N", "R130C", "R131N", "R131C", 
                       "R132N", "R132C", "R133N", "R133C", "R134N")
        col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                     "I130N", "I130C", "I131N", "I131C", 
                     "I132N", "I132C", "I133N", "I133C", "I134N")
      } else if (TMT_plex == 11) {
        col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                       "R130N", "R130C", "R131N", "R131C")
        col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                     "I130N", "I130C", "I131N", "I131C")
      } else if (TMT_plex == 10) {
        col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                       "R130N", "R130C", "R131")
        col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                     "I130N", "I130C", "I131")
      } else if(TMT_plex == 6) {
        col_ratio <- c("R127", "R128", "R129", "R130", "R131")
        col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
      } else {
        col_ratio <- NULL
        col_int <- NULL
      }
      
      df <- df[, "pep_scan_title", drop = FALSE] %>% 
        dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', 
                                      .$pep_scan_title)) %>% 
        dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', 
                                      .$RAW_File)) 
      
      raws <- unique(df$RAW_File)

      if (!purrr::is_empty(raws)) {
        out_nm <- paste0("psm_raws_", 
                         .x %>% gsub("\\.[^.]*$", "", .) %>% gsub("_hdr_rm$", "", .), 
                         ".txt")
        data.frame(RAW_File = raws) %>% 
          write.table(file.path(dat_dir, out_nm), sep = "\t", 
                      col.names = TRUE, row.names = FALSE, quote = FALSE)
        message("MS file names stored in ", file.path(dat_dir, out_nm))
      }
    })

    unlink(file.path(dat_dir, "PSM/temp"), recursive = TRUE, force = TRUE)
  } 
  
  find_mq_psmraws <- function(filelist) {
    purrr::walk(filelist, ~ {
      df <- read.csv(file.path(dat_dir, .x), 
                     check.names = FALSE, header = TRUE, 
                     sep = "\t", comment.char = "#") 
      
      raws <- unique(df$`Raw file`)
      
      if (!purrr::is_empty(raws)) {
        out_nm <- paste0("psm_raws_", gsub("\\.[^.]*$", "", .x), ".txt")
        data.frame(RAW_File = raws) %>% 
          write.table(file.path(dat_dir, out_nm), sep = "\t", col.names = TRUE, 
                      row.names = FALSE, quote = FALSE)
        message("MS file names stored in ", file.path(dat_dir, out_nm))
      }
    })
  }
  
  find_sm_psmraws <- function(filelist) {
    purrr::walk(filelist, ~ {
      df <- readr::read_delim(file.path(dat_dir, .x), delim = ";") %>% 
        dplyr::mutate(filename = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", filename))
      
      raws <- unique(df$filename)
      if (!purrr::is_empty(raws)) {
        out_nm <- paste0("psm_raws_", gsub("\\.[^.]*$", "", .x), ".txt")
        data.frame(RAW_File = raws) %>% 
          write.table(file.path(dat_dir, out_nm), sep = "\t", col.names = TRUE, 
                      row.names = FALSE, quote = FALSE)
        message("MS file names stored in ", file.path(dat_dir, out_nm))
      }
    })
  }
  
  find_mf_psmraws <- function(filelist) {
    purrr::walk(filelist, ~ {
      df <- readr::read_delim(file.path(dat_dir, .x), delim = "\t") %>% 
        # maybe only two fields for none ion-mobility data???
        dplyr::mutate(Spectrum = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", Spectrum))
      
      raws <- unique(df$Spectrum)
      if (!purrr::is_empty(raws)) {
        out_nm <- paste0("psm_raws_", gsub("\\.[^.]*$", "", .x), ".txt")
        data.frame(RAW_File = raws) %>% 
          write.table(file.path(dat_dir, out_nm), sep = "\t", col.names = TRUE, 
                      row.names = FALSE, quote = FALSE)
        message("MS file names stored in ", file.path(dat_dir, out_nm))
      }
    })
  }
  
  
  type <- rlang::enexpr(type)
  if (type == rlang::expr(c("mascot", "maxquant", "spectrum_mill", "msfragger"))) {
    type <- "mascot"
  } else {
    type <- rlang::as_string(type)
  }
  
  type <- tolower(type)
  if (type == "ms") type <- "mascot"
  if (type == "mq") type <- "maxquant"
  if (type == "sm") type <- "spectrum_mill"
  if (type == "mf") type <- "msfragger"

  if (is.null(dat_dir)) {
    dat_dir <- get_gl_dat_dir()
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }

  pattern <- switch(type, 
    mascot = "^F[0-9]+\\.csv$", 
    maxquant = "^msms.*\\.txt$",
    spectrum_mill = "^PSMexport.*\\.ssv$", 
    msfragger = "^psm.*\\.tsv$", 
    stop("Data type needs to be one of `mascot`, `maxquant`,", 
         " `spectrum_mill` or `msfragger`.", 
         call. = FALSE)
  )
  
  filelist <- list.files(path = file.path(dat_dir), pattern = pattern)
  if (purrr::is_empty(filelist)) 
    stop("No ", toupper(type), " files(s) under ", dat_dir, 
         call. = FALSE)

  switch (type,
    mascot = find_mascot_psmraws(filelist),
    maxquant = find_mq_psmraws(filelist),
    spectrum_mill = find_sm_psmraws(filelist),
    msfragger = find_mf_psmraws(filelist),
  )
}


#' Removes PSM headers
#'
#' \code{rmPSMHeaders} removes the header of PSM from
#' \href{https://http://www.matrixscience.com/}{Mascot} outputs. It also
#' removes the spacer columns in the fields of ratio and intensity values.
#'
#' @return Intermediate PSM table(s).
#'
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
rmPSMHeaders <- function () {
	old_opts <- options(warn = 0)
	options(warn = 1)
	on.exit(options(old_opts), add = TRUE)
	
	on.exit(message("Remove PSM headers --- Completed."), 
	        add = TRUE)

	dat_dir <- get_gl_dat_dir()
	
	filelist <- list.files(path = file.path(dat_dir), 
	                       pattern = "^F[0-9]+\\.csv$")

	if (purrr::is_empty(filelist)) {
	  stop("No PSM files(s) with `.csv` extension under ", dat_dir, 
	       call. = FALSE)
	}

	load(file = file.path(dat_dir, "label_scheme_full.rda"))
	TMT_plex <- TMT_plex(label_scheme_full)

	batchPSMheader <- function(filelist, TMT_plex) {
		data_all <- readLines(file.path(dat_dir, filelist))
		output_prefix <- gsub("\\.csv$", "", filelist)

		pep_seq_rows <- grep("pep_seq", data_all)
		
		local({
  		data_header <- data_all[1 : (pep_seq_rows[1] - 1)]
  		data_header <- gsub("\"", "", data_header, fixed = TRUE)
  
  		write.table(data_header, file.path(dat_dir, "PSM/cache",
  		            paste0(output_prefix, "_header", ".txt")),
  		            sep = "\t", col.names = FALSE, row.names = FALSE)
		})

		if (length(pep_seq_rows) > 1) {
		  local({
  		  data_queries <- data_all[pep_seq_rows[2] : length(data_all)]
  		  data_queries <- gsub("\"", "", data_queries, fixed = TRUE)
  		  writeLines(data_queries, 
  		             file.path(dat_dir, "PSM/cache", 
  		                       paste0(output_prefix, "_queries.csv")))
  		  prep_queries()
		  })
		}

		unassign_hits_row <- 
		  grep("Peptide matches not assigned to protein hits", data_all)
		
		if (! purrr::is_empty(unassign_hits_row)) {
		  psm_end_row <- unassign_hits_row - 2
		} else if (length(pep_seq_rows) > 1) {
		  psm_end_row <- pep_seq_rows[2] - 4
		} else {
		  psm_end_row <- length(data_all)
		}
	
		local({
		  mascot_tmtplex <- find_masct_tmtplex(data_all)

		  if (mascot_tmtplex != TMT_plex) {
		    warning("Mascot PSMs suggest a TMT ", mascot_tmtplex, "-plex, ", 
		            "which is different to the ", TMT_plex, "-plex in `expt_smry.xlsx`.", 
		            call. = FALSE)
		  }
		})
		
		data_psm <- data_all[pep_seq_rows[1] : psm_end_row]
		data_psm <- gsub("\"---\"", -1, data_psm, fixed = TRUE)
		data_psm <- gsub("\"###\"", -1, data_psm, fixed = TRUE)
		
		# --- for simulated data
		data_psm <- gsub("---", -1, data_psm, fixed = TRUE)
		data_psm <- gsub("###", -1, data_psm, fixed = TRUE)
		
		if (TMT_plex > 0) {
		  data_psm[1] <- paste0(data_psm[1], paste(rep(",", TMT_plex * 4 -2), 
		                                           collapse = ''))
		}
		
		writeLines(data_psm, 
		           file.path(dat_dir, "PSM/cache", 
		                     paste0(output_prefix, "_hdr_rm.csv")))
		rm(data_psm)
	}

	purrr::walk(filelist, batchPSMheader, TMT_plex)
}


#' Prepare query files from extended Mascot exports
#' 
prep_queries <- function() {
  
  prep_queries_sub <- function(filelist, dat_dir) {
    fn_prefix <- gsub("\\.[^.]*$", "", filelist)
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filelist)
    filepath <- file.path(dat_dir, "PSM/cache")
  
    df <- readLines(file.path(dat_dir, "PSM/cache", filelist))
  
    nms <- df[1] %>% stringr::str_split(",") %>% unlist()
    
    ind_a <- local({
      ind <- which(nms == "NumVals")
      
      if (purrr::is_empty(ind)) {
        ind <- which(nms == "StringIons1") - 1
        if (purrr::is_empty(ind)) {
          ind <- which(nms == "pep_rank") - 1
          if (purrr::is_empty(ind)) {
            stop("Column `StringIons1` or `pep_rank` not found.", 
                 call. = FALSE)
          }
        } 
      }
      
      return(ind)
    })
    
    nms_a <- nms[1:ind_a]
    nms_c <- nms %>% 
      .[(ind_a + 1): length(.)] %>% 
      .[!grepl("^StringIons", .)]
    
    vals <- df[-1] %>% stringr::str_split(",")
    
    # before msms `StringIons`
    vals_a <- vals %>% 
      purrr::map(`[`, 1:ind_a) %>% 
      do.call(rbind, .) %>% 
      data.frame() %>% 
      `colnames<-`(nms_a) 
    
    vals <- vals %>% purrr::map(~ .x[-c(1:ind_a)])
  
    # ion table
    try(
      if ("StringIons1" %in% nms) {
        val_msms <- vals %>% 
          purrr::map(~ .x[grepl(":", .x)]) %>% 
          purrr::map(~ .x[!grepl("^File:", .x)]) %T>% 
          saveRDS(file.path(filepath, paste0(fn_prefix, "_msms.rds")))         
      }
    )

    ## results after str_split:
    # (1) with non_empty "StringIons1", "StringIons2" and "StringIons3"
    # e.g. "670.380430:8541", "710.409550:7887",... , "", ""
    # (2) with empty "StringIons1", "StringIons2" and "StringIons3"
    # "", "", ""
    # so entries without msms "StringIons" will have one more column
    
    # vals_c: information after `msms` e.g. phopho location probability
    local({
      vals_c <- vals %>% 
        purrr::map(~ .x[!grepl("[0-9]+\\.{0,1}[0-9]*:", .x)]) %>% 
        purrr::map(~ .x[-length(.x)])
      
      min_inds <- vals_c %>% 
        purrr::map(~ which(nchar(.x) != 0)) %>% 
        purrr::map(`[`, 1) %>% 
        purrr::map(~ replace(.x, is.na(.x), 0))
      
      vals_c <- purrr::map2(vals_c, min_inds, ~ {
        .x[.y : length(.x)]
      }) 
      
      vals_c <- vals_c %>% 
        plyr::ldply(rbind) %>% 
        `colnames<-`(nms_c)

      vals_ac <- dplyr::bind_cols(vals_a, vals_c)
      vals_ac[vals_ac == ""] <- NA
      
      saveRDS(vals_ac, file.path(filepath, paste0(fn_prefix, "_conf.rds")))
    })
  }
  
  
  dat_dir <- get_gl_dat_dir()
  
  filepath <- file.path(dat_dir, "PSM/cache")
  filelist = list.files(path = filepath, pattern = "^F[0-9]+_queries.csv$")
  
  if (length(filelist) > 0) {
    message(paste("Process PSM query files under", filepath))
    purrr::walk(filelist, prep_queries_sub, dat_dir)
  }
}


#' Add column "pep_var_mod_conf"
#' 
#' @param df A data frame of PSMs
#' @inheritParams load_expts
add_mod_conf <- function(df, dat_dir) {
  dat_file <- unique(df$dat_file)
  stopifnot(length(dat_file) == 1)
  
  filepath <- file.path(dat_dir, "PSM/cache", paste0(dat_file, "_queries_conf.rds"))
  
  if (! file.exists(filepath)) return(df)
  
  # (1) assume Mascot does not report multiple "pep_var_mod_pos" 
  #   under the same c("pep_query", "pep_seq")
  # (2) no need of `dat_file` being part of the `uniq_by` 
  #   as `df` has been split under a give `dat_file`
  uniq_by <- c("pep_query", "pep_seq")
  
  queries <- readRDS(filepath) %>% 
    dplyr::filter(!is.na(pep_rank), !is.na(pep_var_mod_conf)) %>% 
    dplyr::rename(pep_query = query_number) %>% 
    tidyr::unite(uniq_id, uniq_by, sep = ".", remove = FALSE) 

  if (! "pep_var_mod_conf" %in% names(queries)) return(df)
  
  df <- df %>% 
    tidyr::unite(uniq_id, uniq_by, sep = ".", remove = FALSE) %>% 
    dplyr::left_join(queries[, c("uniq_id", "pep_var_mod_conf")], 
                     by = "uniq_id") %>% 
    dplyr::mutate(.n = row_number()) %>% 
    dplyr::mutate(pep_var_mod_conf = 
                    as.numeric(sub("%", "", pep_var_mod_conf))) %>% 
    dplyr::mutate(pep_var_mod_conf = pep_var_mod_conf/100) %>% 
    dplyr::arrange(-pep_var_mod_conf) %>% 
    dplyr::group_by(uniq_id) 
  
  df_first <- df %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup(uniq_id) %>% 
    dplyr::arrange(.n) 
  
  df_second <- df %>% 
    dplyr::filter(row_number() == 2) %>% 
    dplyr::ungroup(uniq_id) %>% 
    dplyr::arrange(.n) %>% 
    dplyr::select(uniq_id, pep_var_mod_conf) %>%
    dplyr::rename(pep_var_mod_conf_2 = pep_var_mod_conf)
  
  ## possible `pep_var_mod_conf` values are all NA; 
  #  thus `pep_locprob` & `pep_locdiff` both be NA 
  # query_number                  pep_seq              pep_var_mod_pos pep_var_mod_conf
  # 1       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400000000000000050002000.0             <NA>
  # 2       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400000000000000500002000.0             <NA>
  # 3       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400000000005000000002000.0             <NA>
  # 4       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400050000000000000002000.0             <NA>
  
  df <- df_first %>% 
    dplyr::left_join(df_second, by = "uniq_id") %>% 
    dplyr::mutate(pep_locdiff = pep_var_mod_conf - pep_var_mod_conf_2) %>% 
    dplyr::select(-pep_var_mod_conf_2) %>% 
    dplyr::rename(pep_locprob = pep_var_mod_conf) %>% 
    dplyr::arrange(.n) %>% 
    dplyr::select(-uniq_id, -.n)
}


#' Add the `pep_seq_mod` field to Mascot PSMs
#' 
#' @inheritParams locate_outliers
#' @inheritParams splitPSM
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_mascot_pepseqmod <- function(df, use_lowercase_aa, purge_phosphodata) {
  dat_dir <- get_gl_dat_dir()
  
  stopifnot("dat_file" %in% names(df))
  
  dat_id <- df$dat_file %>% unique()
  
  dat_file <- file.path(dat_dir, "PSM/cache", paste0(dat_id, "_header.txt"))
  stopifnot(length(dat_id)== 1, file.exists(dat_file))
  
  df_header <- readLines(dat_file)

  fixed_mods <- df_header[((grep("Fixed modifications", df_header))[1] + 3) : 
                            ((grep("Variable modifications", df_header))[1] - 2)] %>%
    gsub("\"", "", ., fixed = TRUE) %>%
    data.frame() %>%
    tidyr::separate(".", sep = ",", c("Mascot_abbr", "Description", "Delta_mass"))
  
  var_mods <- df_header[((grep("Variable modifications", df_header))[1] + 3) :
                          ((grep("Search Parameters", df_header))[1] - 2)] %>%
    gsub("\"", "", ., fixed = TRUE) %>%
    data.frame() %>%
    tidyr::separate(".", sep = ",", extra = "drop", 
                    c("Mascot_abbr", "Description", "Delta_mass")) %>%
    dplyr::mutate(Filename = gsub("[\\\\\\/\\:\\*\\?\\'\\<\\>\\|]", ".", Description))
  
  if (is.null(df$pep_seq)) stop("Column `pep_seq` not found.", call. = FALSE)

  if (nrow(var_mods) == 0) {
    df$pep_seq_mod <- df$pep_seq
    return(df)
  } 
  
  is_phospho_expt <- any(grepl("Phospho\\s{1}\\(", var_mods$Description))
  if (is_phospho_expt && purge_phosphodata) {
    if ("pep_var_mod" %in% names(df)) {
      df <- df %>% dplyr::filter(grepl("Phospho\\s{1}\\(", pep_var_mod))
    } else {
      warning("Missing column `pep_var_mod` for non-phosphopeptide removals.", 
              call. = FALSE)
    }
  }
  rm(is_phospho_expt)

  if (!use_lowercase_aa) {
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste0(pep_seq, "[", pep_var_mod_pos, "]"))
  } else {
    df$pep_seq_mod <- df$pep_seq

    # (1) non terminal modifications
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(!grepl("N-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("C-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            pos_matrix  <- gregexpr(mod, df_sub$pep_var_mod_pos) %>%
              plyr::ldply(., rbind) %>%
              # "-2" for the two characters, "0." ..., in 'pep_var_mod_pos'
              purrr::map(function(x) {x - 2}) %>%
              data.frame(check.names = FALSE)
            
            for (k in 1:ncol(pos_matrix)) {
              rows <- !is.na(pos_matrix[, k])
              locales <- pos_matrix[rows, k]
              
              lowers <- substr(df_sub$pep_seq_mod[rows], locales, locales) %>% 
                tolower()
              
              substr(df_sub$pep_seq_mod[rows], locales, locales) <- lowers
            }
            
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (2-1) add "_" to sequences from protein N-terminal acetylation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Acetyl (Protein N-term)", Description, fixed = TRUE))
      
      nrow <- nrow(mod_tbl)
      stopifnot(nrow <= 1)
      
      if (nrow == 1) {
        mod <- mod_tbl$Mascot_abbr[1]
        
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)

        df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
        df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
        
        if (nrow(df_sub) > 0) {
          df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0("_", pep_seq_mod))
          df <- rbind(df_rest, df_sub)
        }        
      }

      return(df)
    })
    
    # (2-2) add "_" to sequences from protein C-terminal amidation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Amidated (Protein C-term)", Description, fixed = TRUE))

      nrow <- nrow(mod_tbl)
      stopifnot(nrow <= 1)
      
      if (nrow == 1) {
        mod <- mod_tbl$Mascot_abbr[1]
        
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
        
        df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
        df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
        
        if (nrow(df_sub) > 0) {
          df_sub <- df_sub %>% 
            dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "_"))
          df <- rbind(df_rest, df_sub)
        }        
      }

      return(df)
    })
    
    # (3-1) "~" for "(Protein N-term)" other than acetylation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Protein N-term", Description, fixed = TRUE)) %>%
        dplyr::filter(!grepl("Acetyl (Protein N-term)", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = paste0("~", pep_seq_mod))
            
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (3-2) "~" for "(Protein C-term)" other than amidation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Protein C-term", Description, fixed = TRUE)) %>%
        dplyr::filter(!grepl("Amidated (Protein C-term)", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "~"))
            
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (4-1) "^" for peptide "(N-term)"  
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("N-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Protein N-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = 
                              gsub("(^[_~]{0,1})(.)", 
                                   paste0("\\1", "^", "\\2"), 
                                   pep_seq_mod)) 
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (4-2) "^" peptide "(C-term)" 
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("C-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Protein C-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
          
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = gsub("(.)([_~]{0,1}$)", 
                                               paste0("\\1", "^", "\\2"), 
                                               pep_seq_mod)) 
            
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })

    # (5) paste "pep_res_before" and "pep_res_after"
    df <- df %>%
      dplyr::mutate(pep_seq = 
                      paste(pep_res_before, pep_seq, pep_res_after, 
                            sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = 
                      paste(pep_res_before, pep_seq_mod, pep_res_after, 
                            sep = "."))
  }
  
  # `pep_seq` may not be in fasta
  df %>% 
    dplyr::mutate(
      pep_seq = gsub("^NA.", "?.", pep_seq), 
      pep_seq = gsub(".NA$", ".?", pep_seq), 
      pep_seq_mod = gsub("^NA.", "?.", pep_seq_mod), 
      pep_seq_mod = gsub(".NA$", ".?", pep_seq_mod)
    )
}


#' Parallel split of PSM tables
#' 
#' @param df A list of data frames.
#' @param nm Names of the \code{df}'s.
#' 
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams TMT_levels
psm_msplit <- function (df, nm, fn_lookup, dat_dir, plot_rptr_int, TMT_plex) {
  df <- df %>% dplyr::select(-TMT_inj)
  
  out_fn <- fn_lookup %>%
    dplyr::filter(TMT_inj == nm) %>%
    dplyr::select(filename) %>%
    unique() %>%
    unlist() %>%
    paste0(., ".csv")
  
  df <- df %>% 
    dplyr::rename(raw_file = RAW_File) %T>% 
    write.csv(file.path(dat_dir, "PSM/cache", out_fn), row.names = FALSE)
  
  if (plot_rptr_int) {
    if (TMT_plex > 0) {
      df_int <- df %>% .[, grepl("^I[0-9]{3}", names(.))]
    } else {
      # change "Precursor Intensity", "Intensity" to "I000" later
      df_int <- df %>% 
        dplyr::select(which(names(.) %in% 
                              c("Precursor Intensity", "Intensity", "I000"))) %>% 
        `colnames<-`("I000")
    }
    
    rptr_violin(df = df_int, 
                filepath = file.path(dat_dir, "PSM/rprt_int/raw", 
                                     gsub("\\.csv", "\\.png", out_fn)), 
                width = 8, height = 8)
  }
  
  invisible(out_fn)
}


#' Find peptides that are unique by protein entries
#' 
#' @param df PSM data.
#' @param pep_id the column key of peptide sequences.
#' @param prot_id The column key of protein accessions.
find_literal_unique_peps <- function(df, pep_id = "pep_seq", prot_id = "prot_acc") {
  stopifnot(pep_id %in% names(df), prot_id %in% names(df))
  
  uniq_by <- c(pep_id, prot_id)
  
  peps_shared_by_prots <- df %>% 
    tidyr::unite(uniq_id, uniq_by, sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(uniq_id)) %>% 
    tidyr::separate(uniq_id, into = c(pep_id, prot_id), 
                    sep = "@", remove = TRUE) %>% 
    dplyr::group_by(!!rlang::sym(pep_id)) %>% 
    dplyr::summarise(pep_n_prot = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(pep_n_prot >= 2) %>% 
    .[[pep_id]]
  
  invisible(peps_shared_by_prots)
}


#' Splits PSM tables
#'
#' \code{splitPSM} splits the PSM outputs after \code{rmPSMHeaders()}. It
#' separates PSM data by TMT/LFQ experiment and LC/MS injection.
#'
#' Arguments \code{group_psm_by} and \code{group_pep_by} are used to update
#' \code{prot_matches_sig} and \code{prot_sequences_sig} after data merge.
#'
#' @param fasta Character string(s) to the name(s) of fasta file(s) with
#'   prepended directory path. The \code{fasta} database(s) need to match those
#'   used in MS/MS ion search. There is no default and users need to provide the
#'   correct file path(s) and name(s).
#' @param entrez Character string(s) to the name(s) of entrez file(s) with
#'   prepended directory path. At the \code{NULL} default, a convenience lookup
#'   is available for species among \code{c("human", "mouse", "rat")}. For other
#'   species, users need to provide the file path(s) and name(s) for the lookup
#'   table(s). See also \code{\link{Uni2Entrez}} and \code{\link{Ref2Entrez}}
#'   for preparing custom entrez files.
#' @param rm_craps Logical; if TRUE,
#'   \href{https://www.thegpm.org/crap/}{cRAP} proteins will be removed.
#'   The default is FALSE.
#' @param rm_krts Logical; if TRUE, keratin entries will be removed. The default
#'   is FALSE.
#' @param annot_kinases Logical; if TRUE, proteins of human or mouse origins
#'   will be annotated with their kinase attributes. The default is FALSE.
#' @param rptr_intco Numeric; the threshold of reporter ion intensity. The
#'   default is 0 without cut-offs. The filter is also applicable to the
#'   precursor intensity with MaxQuant LFQ.
#' @param plot_rptr_int Logical; if TRUE, the distributions of reporter-ion
#'   intensities will be plotted. The default is TRUE. The argument is also
#'   applicable to the precursor intensity with MaxQuant LFQ.
#' @param use_lowercase_aa Logical; if TRUE, modifications in amino acid
#'   residues will be abbreviated with lower-case and/or \code{^_~}. See the
#'   table below for details. The default is TRUE. 
#' @param purge_phosphodata Logical; if TRUE and phosphorylation present as
#'   variable modification(s), entries without phosphorylation will be removed.
#'   The default is TRUE.
#' @param parallel Logical; if TRUE, performs parallel computation. The default
#'   is TRUE.
#' @inheritParams annotPSM
#' @import dplyr tidyr stringr
#' @importFrom magrittr %>% %T>% %$% %<>% 
splitPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                     fasta = NULL, entrez = NULL, pep_unique_by = "group", 
                     rm_craps = FALSE, rm_krts = FALSE, purge_phosphodata = TRUE, 
                     annot_kinases = FALSE, plot_rptr_int = TRUE, rptr_intco = 0, 
                     use_lowercase_aa = TRUE, parallel = TRUE, ...) {
  
  on.exit(message("Split PSMs by TMT experiments and LCMS series --- Completed."), 
          add = TRUE)

  dat_dir <- get_gl_dat_dir()
	load(file = file.path(dat_dir, "label_scheme_full.rda"))
	load(file = file.path(dat_dir, "fraction_scheme.rda"))

	TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist = list.files(path = file.path(dat_dir, "PSM/cache"),
                        pattern = "^F[0-9]+_hdr_rm.csv$")

	if (length(filelist) == 0) {
	  stop(paste("Missing intermediate PSM file(s) under: ", 
	             file.path(dat_dir, "PSM/cache")), 
	       call. = FALSE)
	}

  df <- purrr::map(filelist, pad_mascot_channels) %>% 
    pad_mascot_fields() %>% 
    do.call(rbind, .)
  
  # --- remove spacer columns ---
  df <- local({
    r_start <- which(names(df) == "")[1]
    int_end <- ncol(df) - 1
    if (int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
    
    col_int <- find_int_cols(TMT_plex)
    col_ratio <- find_ratio_cols(TMT_plex)
    
    if (TMT_plex > 0) {
      colnames(df)[r_start:(r_start+TMT_plex-2)] <- col_ratio
      colnames(df)[(r_start+TMT_plex-1):(r_start+TMT_plex+TMT_plex-2)] <- col_int
    }
    
    return(df)
  })
  
  # --- remove rows of emPAI etc.
  # partial columns of emPAI etc. form new rows with "empty" prot_acc
  df <- df %>% dplyr::filter(prot_acc != "")
    
  # --- check the grouping of proteins ---
  if (! "prot_family_member" %in% names(df)) {
    warning("PSMs were exported from Mascot without protein grouping.", 
            call. = FALSE)
    
    df <- df %>% 
      dplyr::mutate(prot_family_member = 0) %>% 
      reloc_col_after("prot_family_member", "prot_hit_num")
  }

  df <- local({
    df <- df %>% 
      dplyr::mutate(pep_isbold = ifelse(pep_isbold == 1, TRUE, FALSE), 
                    pep_isunique = ifelse(pep_isunique == 1, TRUE, FALSE)) %>% 
      add_entry_ids("pep_seq", "pep_index") %>% 
      add_entry_ids("prot_acc", "prot_index")
    
    peps_shared_by_prots <- df %>% 
      find_literal_unique_peps("pep_seq", "prot_acc")
    
    df <- df %>% 
      dplyr::mutate(pep_literal_unique = 
                      ifelse(pep_seq %in% peps_shared_by_prots, FALSE, TRUE), 
                    pep_razor_unique = 
                      pep_isunique) %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
    
    if (pep_unique_by == "group") {
      df <- df
    } else if (pep_unique_by == "protein") {
      df <- df %>% 
        dplyr::mutate(pep_isunique = pep_literal_unique)
    } else if (pep_unique_by == "none") {
      df <- df %>% 
        dplyr::mutate(pep_isunique = TRUE)
    } else {
      df <- df
    }
  })
  
  # convenience crap removals where the fasta names ended with "|"
  if (rm_craps) {
    df <- df %>% dplyr::filter(!grepl("\\|.*\\|$", prot_acc))
  }

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `F[...].csv` ", 
          "or `PSM/cache/[...]_hdr_rm.csv` for `filter_` varargs.")
  
  # (1) the same `pep_query` can be assigned to different `pep_seq` at different `pep_rank`
  # (2) the same combination of `pep_query`, `pep_seq` and `pep_var_mod_pos` 
  #     (positional difference) 
  # can be assigned to different `prot_acc` 
  
  # remove redundant peptides under each `dat_file`
  # dbl-dipping peptides across `dat_file` will be handled in `mergePep`
  df <- local({
    uniq_by <- c("pep_query", "pep_seq", "pep_var_mod_pos")
    
    if (length(unique(df$dat_file)) >= 1) {
      uniq_by <- c(uniq_by, "dat_file")
    }
    
    df <- df %>% 
      tidyr::unite(uniq_id, uniq_by, sep = ".", remove = FALSE) %>% 
      dplyr::mutate(.n = row_number()) %>% 
      dplyr::arrange(-prot_matches, -pep_isbold, -prot_mass) %>% 
      dplyr::filter(!duplicated(uniq_id)) %>% 
      dplyr::arrange(.n) %>% 
      dplyr::select(-uniq_id, -.n)
    
    ## before
    # prot_acc   I126 I127N   I127C
    # 1  2::HBB1_MOUSE  29410 20920   27430
    # 2  2::HBE_MOUSE  29410 20920   27430
    # 3  2::HBE_MOUSE 137000 38330 1115000
    
    ## after
    # prot_acc   I126 I127N   I127C
    # 1  2::HBB1_MOUSE  29410 20920   27430
    # 2  2::HBE_MOUSE 137000 38330 1115000
    
    # remove subset proteins
    if ("prot_family_member" %in% names(df)) {
      df <- df %>% dplyr::filter(!is.na(prot_family_member))
    }
    
    # `pep_seq` changed from such as MENGQSTAAK to K.MENGQSTAAK.L
    df <- df %>% 
      split(., .$dat_file, drop = TRUE) %>% 
      purrr::map(add_mod_conf, dat_dir) %>% 
      purrr::map(add_mascot_pepseqmod, use_lowercase_aa, purge_phosphodata) %>% 
      dplyr::bind_rows() 
    
    df <- df %>% 
      dplyr::mutate(prot_acc_orig = prot_acc) %>% 
      dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc)) %>% 
      annotPrn(fasta, entrez) %>% 
      dplyr::mutate(prot_acc = prot_acc_orig) %>% 
      dplyr::select(-prot_acc_orig)
  })
  
  df <- local({
    prot_matches_sig <- df %>%
      dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
      dplyr::summarise(prot_matches_sig_new = n())
    
    prot_sequences_sig <- df %>%
      dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
      dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
      dplyr::summarise(prot_sequences_sig_new = n())
    
    df <- list(df, prot_matches_sig, prot_sequences_sig) %>% 
      purrr::reduce(dplyr::left_join, by = group_pep_by) %>% 
      dplyr::mutate(prot_matches_sig = prot_matches_sig_new, 
                    prot_sequences_sig = prot_sequences_sig_new) %>%
      dplyr::select(-prot_matches_sig_new, -prot_sequences_sig_new)
  })
  
  # `filter_dots` placed after the calculation of prot_matches_sig & prot_sequences_sig: 
  #    to reflect counts before `filter_dots` so not to confuse users by any mismatches.
  # remove prefix 1:: etc. here in case users try varargs against the prefix
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc))

  acc_types <- unique(df$acc_type) %>% .[!is.na(.)]

  # re-apply craps after annotation
  if (rm_craps) {
    craps <- load_craps(acc_types)
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% 
      dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) {
    df <- df %>% 
      split(., .$acc_type, drop = TRUE) %>% 
      .[acc_types] %>% 
      purrr::imap( ~ annotKin(.x, .y)) %>% 
      do.call(rbind, .)
  } else {
    df <- df %>% 
      dplyr::mutate(kin_attr = FALSE, kin_class = NA, kin_order = NA)
  }

  # easier to just recalculate
  df$pep_miss <- NULL
  df <- df %>% 
    annotPeppos() %>% 
    dplyr::mutate(pep_seq_bare = gsub("^.*\\.([^\\.]+)\\..*", "\\1", pep_seq)) %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq_bare)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq_bare), 
                                    stringr::str_count(pep_seq_bare, "[KR]") - 1,
                                    stringr::str_count(pep_seq_bare, "[KR]"))) %>% 
    add_prot_icover() %>% 
    dplyr::select(-pep_seq_bare) %>% 
    reloc_col_after("pep_seq_mod", "pep_seq") %>% 
    reloc_col_before("pep_res_before", "pep_seq") %>% 
    reloc_col_after("pep_res_after", "pep_seq_mod") %>% 
    reloc_col_after("pep_start", "pep_res_after") %>% 
    reloc_col_after("pep_end", "pep_start") %>% 
    reloc_col_after("pep_len", "pep_end") %>% 
    reloc_col_after("pep_miss", "pep_len")
  
  if (!("prot_cover" %in% names(df) && length(filelist) == 1)) {
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by))
  } 
  
  # --  column arrangements --- 
  df <- df %>% 
    reloc_col_after("prot_index", "prot_cover") 
  
  df <- df %>% 
    order_psm_cols(c("gene", "fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                     "other_acc", "entrez", "species", "acc_type", 
                     "kin_attr", "kin_class", "kin_order", "is_tryptic"))
  
  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    # Shared peptides will be removed in this step by `rowSums` if 
    #   checked 'Unique peptide only' during Mascot PSM export
    df <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), 
                       as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                       ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                       ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0)
  }
  
  df <- df %>%
    dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', 
                                  .$pep_scan_title)) %>% 
    dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', 
                                  .$RAW_File)) %>% # Bruker
    dplyr::mutate(prot_acc = gsub("\\d::", "", 
                                  .$prot_acc)) %>%
    dplyr::arrange(RAW_File, pep_seq, prot_acc) %>% 
    reloc_col_before("dat_file", "RAW_File")
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))),
    df %>% dplyr::select(grep("^pep_", names(.))),
    df %>% dplyr::select(-grep("^prot_|^pep_", names(.))),
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("^[IR][0-9]{3}[NC]{0,1}", names(.))),
    df %>% dplyr::select(grep("^[IR][0-9]{3}[NC]{0,1}", names(.))),
  )

  # --- compare experimental and user-provided PSM files ---
  res <- check_raws(df)
  df <- res$df
  tmtinj_raw_map <- res$lookup
  rm(res)
  
  # --- split by TMT and LCMS ---
  if (! "PSM_File" %in% names(tmtinj_raw_map)) {
    df_split <- df %>%
      dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>% 
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>% 
      split(., .$TMT_inj, drop = TRUE)
  } else {
    tmtinj_raw_map <- tmtinj_raw_map %>% 
      dplyr::mutate(PSM_File = gsub("\\.csv$", "", PSM_File)) %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "PSM_File"), sep = "@") 
    
    df_split <- df %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "dat_file"), sep = "@", remove = FALSE) %>% 
      dplyr::left_join(tmtinj_raw_map, by = "RAW_File2") %>% 
      dplyr::select(-RAW_File2) %>% 
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>% 
      split(., .$TMT_inj, drop = TRUE)
  }

  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    warning("TMT sets and LC/MS injections not have corresponindg PSM files:\n", 
            call. = FALSE)
    message(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` from experiment summary file."), 
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), 
                                   sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")

  df_split <- df_split %>% .[names(.) %>% sort_tmt_lcms()]
  
  if (parallel) {
    n_cores <- pmin(parallel::detectCores(), length(filelist))
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

    suppressWarnings(
      silent_out <- parallel::clusterMap(
        cl, psm_msplit, df_split, names(df_split), 
        MoreArgs = list(fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
      )
    )
    
    parallel::stopCluster(cl)
  } else {
    purrr::walk2(df_split, names(df_split), psm_msplit, 
                 fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
  }
}


#' Pad Mascot TMT channels to the highest plex
#' 
#' @param df An intermediate PSM table from Mascot export
pad_mascot_channels <- function(file) {
  find_mascot_plex <- function (df) {
    to_126 <- grep("/126", df[1, ], fixed = TRUE)
    
    if (purrr::is_empty(to_126)) {
      return (0)
    } else {
      return (length(to_126) + 1)
    }
  }
  
  make_mascot_ratios <- function(TMT_plex, nrow) {
    if (TMT_plex == 16) {
      col_keys <- c("127N/126", "127C/126", "128N/126", "128C/126", "129N/126", "129C/126",
                    "130N/126", "130C/126", "131N/126", "131C/126", 
                    "132N/126", "132C/126", "133N/126", "133C/126", "134N/126")
    } else if (TMT_plex == 11) {
      col_keys <- c("127N/126", "127C/126", "128N/126", "128C/126", "129N/126", "129C/126",
                    "130N/126", "130C/126", "131N/126", "131C/126")
    } else if (TMT_plex == 10) {
      col_keys <- c("127N/126", "127C/126", "128N/126", "128C/126", "129N/126", "129C/126",
                    "130N/126", "130C/126", "131/126")
    } else if(TMT_plex == 6) {
      col_keys <- c("127/126", "128/126", "129/126", "130/126", "131/126")
    } else {
      col_keys <- NULL
    }
    
    col_keys <- rep(col_keys, each = 2) 
    col_keys[which(seq_along(col_keys)%%2 == 0)] <- NA
    
    suppressMessages(purrr::map(col_keys, rep, each = nrow) %>% 
      dplyr::bind_cols() %>% 
      `colnames<-`(""))
  }
  
  make_mascot_intensities <- function(TMT_plex, nrow) {
    if (TMT_plex == 16) {
      col_keys <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                    "130N", "130C", "131N", "131C", 
                    "132N", "132C", "133N", "133C", "134N")
    } else if (TMT_plex == 11) {
      col_keys <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                    "130N", "130C", "131N", "131C")
    } else if (TMT_plex == 10) {
      col_keys <- c("126", "127N", "127C", "128N", "128C", "129N", "129C",
                    "130N", "130C", "131")
    } else if(TMT_plex == 6) {
      col_keys <- c("126", "127", "128", "129", "130", "131")
    } else {
      col_keys <- NULL
    }
    
    col_keys <- rep(col_keys, each = 2) 
    col_keys[which(seq_along(col_keys)%%2 == 0)] <- -1
    
    suppressMessages(purrr::map(col_keys, rep, each = nrow) %>% 
      `names<-`(paste0("nm_", seq_along(.))) %>% 
      dplyr::bind_cols() %>% 
      `colnames<-`(""))
  }

  find_padding_pos <- function (this_plex, TMT_plex) {
    if (this_plex == 6) {
      if (TMT_plex == 10) {
        pos <- c(3, 5, 7, 9)
      } else if (TMT_plex == 11) {
        pos <- c(3, 5, 7, 9, 11)
      } else if (TMT_plex == 16) {
        pos <- c(3, 5, 7, 9, 11:16)
      } else {
        stop("TMT_plex not one of c(10, 11, 16).", call. = FALSE)
      }
    } else if (this_plex == 10) {
      if (TMT_plex == 11) {
        pos <- 11
      } else if (TMT_plex == 16) {
        pos <- c(11:16)
      } else {
        stop("TMT_plex not one of c(11, 16).", call. = FALSE)
      }
    } else if (this_plex == 11) {
      if (TMT_plex == 16) {
        pos <- c(12:16)
      } else {
        stop("TMT_plex not one of c(16).", call. = FALSE)
      }
    }
    
    invisible(pos)
  }
  
  
  dat_dir <- get_gl_dat_dir()
  
  base_name <- file %>% gsub("_hdr_rm\\.csv", "", .)
  
  df <- read.delim(file.path(dat_dir, "PSM/cache", file), 
                   sep = ',', check.names = FALSE, 
                   header = TRUE, stringsAsFactors = FALSE, 
                   quote = "\"", fill = TRUE , skip = 0)

  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  
  this_plex <- find_mascot_plex(df)
  TMT_plex <- TMT_plex(label_scheme_full)
  stopifnot(this_plex <= TMT_plex, this_plex >= 0)
  
  if (this_plex > 0 && this_plex < TMT_plex) {
    pos <- find_padding_pos(this_plex, TMT_plex)
    
    nas <- data.frame(rep(NA, nrow(df)))
    
    cols_both <- which(names(df) %>% stringr::str_count() == 0)

    df_ratio <- local({
      cols_ratio <- cols_both %>% .[1:((this_plex-1)*2)]
      df_ratio_raw <- df[, cols_ratio]
      df_ratio <- df_ratio_raw[, seq_along(cols_ratio) %% 2 %>% as.logical() %>% `!`()]
      df_ratio <- cbind(nas, df_ratio)
      
      for (idx in seq_along(pos)) {
        df_ratio <- suppressMessages(add_cols_at(df_ratio, nas, pos[idx] - 1))
      }
      
      holders_ratio <- make_mascot_ratios(TMT_plex, nrow(df))
      holders_ratio <- cbind(rep("126/126", nrow(df)), nas, holders_ratio)
      
      for (idx in seq_len(TMT_plex)) {
        holders_ratio[, 2*idx] <- df_ratio[, idx]
      }
      
      holders_ratio <- holders_ratio[, -c(1:2)]
      names(holders_ratio) <- ""
      
      return(holders_ratio)
    })

    df_int <- local({
      cols_int <- cols_both %>% .[((this_plex-1)*2+1):(this_plex*4-2)]
      df_int_raw <- df[, cols_int]
      df_int <- df_int_raw[, seq_along(cols_int) %% 2 %>% as.logical() %>% `!`()]
      
      for (idx in seq_along(pos)) {
        df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))
      }
      
      holders_int <- make_mascot_intensities(TMT_plex, nrow(df))
      
      for (idx in seq_len(TMT_plex)) {
        holders_int[, 2*idx] <- df_int[, idx]
      }
      
      names(holders_int) <- ""
      return(holders_int)
    })
    
    df <- suppressMessages(dplyr::bind_cols(
      df %>% .[, which(names(.) != "")],
      df_ratio,
      df_int,
    ))
    
    cols_no_nm <- (ncol(df) - ncol(df_ratio) - ncol(df_int) + 1) : ncol(df)
    colnames(df)[cols_no_nm] <- ""
  }
  
  df$dat_file <- base_name

  invisible(df)
}


#' Pad Mascot PSM exports 
#' 
#' Intensity or ratio columns are excluded
#' @param df An intermediate PSM table from Mascot export
pad_mascot_fields <- function(df) {
  ncols <- purrr::map_dbl(df, ncol)
  
  if (length(unique(ncols)) > 1) {
    nms_union <- local({
      nms <- df %>% 
        purrr::map(~ names(.x) %>% .[. != ""] %>% .[. != "dat_file"]) 
      
      nms_union <- purrr::reduce(nms, union) 
      
      c(nms_union[nms_union != "pep_scan_title"], 
        nms_union[nms_union == "pep_scan_title"])
    })
    
    df <- purrr::map(df, ~ {
      nms <- names(.x) %>% 
        .[. != ""] %>% 
        .[. != "dat_file"]
      
      missing_nms <- setdiff(nms_union, nms)
      
      if (!purrr::is_empty(missing_nms)) {
        for (nm in missing_nms) .x[[nm]] <- NA
      } 
      
      dplyr::bind_cols(
        .x[nms_union], 
        .x %>% 
          .[, ! names(.) %in% nms_union] %>% 
          .[, ! names(.) == "dat_file"], 
        .x["dat_file"],
      )
    })
  }
  
  return(df)
}


#' Parallel PSM cleanup
#' 
#' @param file A PSM file (name).
#' @inheritParams load_expts
#' @inheritParams cleanupPSM
#' @inheritParams TMT_levels
psm_mcleanup <- function(file, rm_outliers, group_psm_by, dat_dir, TMT_plex) {
  df <- read.csv(file.path(dat_dir, "PSM/cache", file), 
                 check.names = FALSE, header = TRUE, comment.char = "#")
  
  stopifnot(group_psm_by %in% names(df))
  
  if (TMT_plex == 0) {
    df$psm_index <- NULL
    fn <- paste0(gsub(".csv", "", file), "_Clean.txt")
    write.table(df, file.path(dat_dir, "PSM/cache", fn), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
    message(file, "processed (no PSM cleanup for MS1-based LFQ).\n")
    
    return(fn)
  }
  
  # remove all "-1" ratio rows
  df <- local({
    N <- sum(grepl("^R[0-9]{3}", names(df)))
    df <- df %>%
      dplyr::mutate(n = rowSums(.[, grep("^R[0-9]{3}", names(.))] == -1, na.rm = TRUE)) %>%
      dplyr::filter(n != N) %>%
      dplyr::select(-n)
  })
  
  set_idx <- file %>% 
    gsub("TMTset(\\d+)_.*", "\\1", .) %>% 
    as.integer()
  
  injn_idx <- file %>% 
    gsub("^TMTset\\d+_LCMSinj(\\d+)\\.csv$", "\\1", .) %>% 
    as.integer()
  
  channelInfo <- channelInfo(dat_dir = dat_dir, set_idx = set_idx, injn_idx = injn_idx)

  # add column "R126"
  df <- local({
    pos_af <- min(grep("^R[0-9]{3}", names(df)))
    df$R126 <- 1
    df <- cbind.data.frame(df[, 1:(pos_af-1)], 
                           R126 = df$R126, 
                           df[, (pos_af):(ncol(df)-1)]) %>%
      dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, 
                       ~ replace(.x, , NA)) %>%
      dplyr::filter(rowSums(!is.na(.[, grep("^I[0-9]{3}", names(.) )])) > 0) %>%
      dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, 
                       ~ replace(.x, , 0)) %>%
      dplyr::mutate_at(.vars = which(names(.)=="R126")-1+channelInfo$emptyChannels, 
                       ~ replace(.x, , NA)) %>%
      dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) # "> 1" not "0"    
  })

  if (rm_outliers) {
    df <- local({
      dfw_split <- df %>%
        dplyr::select(grep("^I[0-9]{3}", names(.))) %>%
        dplyr::mutate(RM = rowMeans(.[, grep("^I[0-9]{3}", 
                                             names(.))[channelInfo$labeledChannels]],
                                    na.rm = TRUE)) %>%
        dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ log2(.x/RM)) %>%
        dplyr::select(-c("RM")) %>%
        `colnames<-`(gsub("I", "X", names(.))) %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         ~ replace(.x, is.infinite(.x), NA)) %>%
        dplyr::bind_cols(df[, c("psm_index", group_psm_by)], .) %>%
        split(., .[[group_psm_by]], drop = TRUE)
      
      range_colRatios <- grep("^X[0-9]{3}", names(dfw_split[[1]]))
      
      dfw_split <- do.call("rbind", lapply(dfw_split, locate_outliers, range_colRatios)) %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         ~ replace(.x, is.infinite(.x), NA)) %>%
        tidyr::unite(pep_seq_i, !!group_psm_by, psm_index, sep = ":") %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         ~ replace(.x, !is.na(.x), 1))
      
      # column order not charged other than the removal of "psm_index"
      df <- df %>% 
        tidyr::unite(pep_seq_i, !!group_psm_by, psm_index, sep = ":") %>%
        dplyr::left_join(dfw_split, by = "pep_seq_i") %>%
        tidyr::separate(pep_seq_i, into = c(group_psm_by, "psm_index"), 
                        sep = ":", remove = TRUE) %>%
        dplyr::select(-psm_index)
    })
    
    df[, grepl("^I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()

    df <- dplyr::bind_cols(
      df[, !grepl("^R[0-9]{3}|^I[0-9]{3}|^X[0-9]{3}", names(df)), drop = FALSE],
      df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df)), drop = FALSE], 
    ) %>% # "> 1" as "R126 == 1"
      dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) %>% 
      dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 + channelInfo$emptyChannels, 
                       ~ replace(.x, , 0))
  } else {
    df <- df %>% dplyr::select(-psm_index)
    
    df <- dplyr::bind_cols(
      df[, !grepl("^R[0-9]{3}|^I[0-9]{3}", names(df)), drop = FALSE],
      df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df)), drop = FALSE], 
    ) %>% 
      dplyr::mutate_at(.vars = which(names(.)=="I126") - 1 +
                         channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
      dplyr::mutate_at(.vars = which(names(.)=="R126") - 1 +
                         channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
      dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1)
  }
  
  fn <- paste0(gsub(".csv", "", file), "_Clean.txt")
  write.table(df, file.path(dat_dir, "PSM/cache", fn), 
              sep = "\t", col.names = TRUE, row.names = FALSE)
  message(file, "processed\n")
  
  invisible(fn)
}


#' Cleans Up PSM results
#'
#' \code{cleanupPSM} removes PSM outliers after \code{splitPSM},
#' \code{splitPSM_mq} or \code{splitPSM_sm}. The outlier removals will be
#' assessed at the basis of per peptide per TMT channel.
#'
#' Dixon's method will be used when \eqn{2 < n \le 25}; Rosner's method will be
#' used when \eqn{n > 25}.
#'
#' @param rm_outliers Logical; if TRUE, PSM outlier removals will be performed
#'   for peptides with more than two identifying PSMs. Dixon's method will be
#'   used when \eqn{2 < n \le 25} and Rosner's method will be used when \eqn{n >
#'   25}. The default is FALSE.
#' @inheritParams splitPSM
#'
#' @import dplyr tidyr
#' @importFrom stringr str_split
cleanupPSM <- function(rm_outliers = FALSE, group_psm_by = "pep_seq", parallel = TRUE) {
  dat_dir <- get_gl_dat_dir()

  load(file = file.path(dat_dir, "label_scheme_full.rda"))
	TMT_plex <- TMT_plex(label_scheme_full)
	
	filelist <- list.files(path = file.path(dat_dir, "PSM/cache"),
	                      pattern = "^TMT.*LCMS.*\\.csv$") %>% 
	  reorder_files()

	if (parallel) {
	  n_cores <- pmin(parallel::detectCores(), length(filelist))
	  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
	  
	  parallel::clusterExport(cl, list("add_cols_at"), 
	                          envir = env_where("add_cols_at"))
	  parallel::clusterExport(cl, list("reloc_col_before"), 
	                          envir = env_where("reloc_col_before"))
	  
	  suppressWarnings(
	    silent_out <- parallel::clusterApply(
	      cl, filelist, psm_mcleanup, rm_outliers, group_psm_by, dat_dir, TMT_plex)
	  )
	  rm(silent_out)

  	parallel::stopCluster(cl)	  
	} else {
	  purrr::walk(filelist, psm_mcleanup, rm_outliers, group_psm_by, dat_dir, TMT_plex)
	}
}


#'Median-centering normalization of PSM data
#'
#'\code{mcPSM} adds fields of \code{log2_R, N_log2_R and N_I} to PSM tables.
#'
#'@param df A data frame containing the PSM table from database searches.
#'@inheritParams channelInfo
#'@inheritParams annotPSM
#'@inheritParams calcPeptide
#'@import dplyr tidyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>%
mcPSM <- function(df, set_idx, injn_idx, mc_psm_by, group_psm_by, group_pep_by) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  channelInfo <- channelInfo(dat_dir = dat_dir, set_idx = set_idx, injn_idx = injn_idx)
  
  label_scheme_sub <- label_scheme_full %>% 
    dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx) 

  dfw <- df[rowSums(!is.na(df[, grepl("^R[0-9]{3}", names(df)), drop = FALSE])) > 1, ] %>%
    dplyr::arrange(pep_seq, prot_acc) %>%
    dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 +
                       channelInfo$emptyChannels, ~ replace(.x, , NaN))
  
  col_sample <- grep("^I[0-9]{3}", names(dfw))
  
  if (length(channelInfo$refChannels) > 0) {
    ref_index <- channelInfo$refChannels
  } else {
    ref_index <- channelInfo$labeledChannels
  }
  
  dfw <- sweep(dfw[, col_sample], 1,
               rowMeans(dfw[, col_sample[ref_index], drop = FALSE], na.rm = TRUE), "/") %>%
    log2(.) %>% 
    dplyr::mutate_all(~ replace(.x, is.infinite(.), NA)) %>% 
    `colnames<-`(gsub("I", "log2_R", names(.)))	%>%
    cbind(dfw, .) %>%
    dplyr::mutate_at(.vars = grep("[I|R][0-9]{3}[NC]{0,1}", names(.)), 
                     ~ replace(.x, is.infinite(.), NA))
  
  col_log2Ratio <- grepl("^log2_R[0-9]{3}[NC]{0,1}", names(dfw))
  
  if (mc_psm_by == "peptide") {
    cf <- dfw %>% 
      dplyr::select(group_psm_by, grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
      dplyr::group_by(!!rlang::sym(group_psm_by)) %>% 
      dplyr::summarise_all(~ median(.x, na.rm = TRUE)) %>% 
      dplyr::summarise_at(.vars = grep("^log2_R[0-9]{3}[NC]{0,1}", names(.)), 
                          ~ median(.x, na.rm = TRUE)) %>% 
      unlist()
  } else if (mc_psm_by == "protein") {
    cf <- dfw %>% 
      dplyr::select(group_pep_by, grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>% 
      dplyr::summarise_all(~ median(.x, na.rm = TRUE)) %>% 
      dplyr::summarise_at(.vars = grep("^log2_R[0-9]{3}[NC]{0,1}", names(.)), 
                          ~ median(.x, na.rm = TRUE)) %>% 
      unlist()
  } else if (mc_psm_by == "psm") {
    cf <- apply(dfw[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  } else {
    warning("Unknown setting in `mc_psm_by`; the default will be used", call. = FALSE)
  }

  dfw <- sweep(dfw[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- sweep(dfw[, grepl("^I[0-9]{3}", names(dfw)), drop = FALSE], 2, 2^cf, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  # not yet SD scaling at PSM levels
  run_scripts <- FALSE
  if (run_scripts) {
    sd_coefs <- calc_sd_fcts_psm(df = dfw, range_log2r = c(5, 95), range_int = c(5, 95), 
                                 set_idx = set_idx, injn_idx = injn_idx)
    
    nm_log2r_n <- names(dfw) %>% .[grepl("^N_log2_R[0-9]{3}[NC]{0,1}", .)]
    nm_int_n <- names(dfw) %>% .[grepl("^N_I[0-9]{3}[NC]{0,1}", .)]
    
    df_z <- mapply(normSD, dfw[, nm_log2r_n, drop = FALSE], 
                   center = 0, SD = sd_coefs$fct, SIMPLIFY = FALSE) %>%
      data.frame(check.names = FALSE) %>%
      `rownames<-`(rownames(df)) 
    
    ref_cols <- label_scheme_sub %>% 
      dplyr::filter(Reference) %>% 
      dplyr::select(TMT_Channel) %>% 
      unlist() %>% 
      gsub("^TMT-", "N_log2_R", .)
    df_z[, ref_cols] <- 0
    rm(ref_cols)
    
    dfw[, nm_log2r_n] <- df_z
  }
  
  dfw <- dfw %>%
    reorderCols(endColIndex = grep("[RI][0-9]{3}", names(.)), col_to_rn = "pep_seq_mod") %>%
    na_zeroIntensity()
}


#'Annotates PSM results
#'
#'\code{annotPSM} adds fields of annotation to Mascot PSM tables after
#'\code{rmPSMHeaders}, \code{splitPSM} and \code{cleanupPSM}.
#'
#'@param group_psm_by A character string specifying the method in PSM grouping.
#'  At the \code{pep_seq} default, descriptive statistics will be calculated
#'  based on the same \code{pep_seq} groups. At the \code{pep_seq_mod}
#'  alternative, peptides with different variable modifications will be treated
#'  as different species and descriptive statistics will be calculated based on
#'  the same \code{pep_seq_mod} groups.
#'@param group_pep_by A character string specifying the method in peptide
#'  grouping. At the \code{prot_acc} default, descriptive statistics will be
#'  calculated based on the same \code{prot_acc} groups. At the \code{gene}
#'  alternative, proteins with the same gene name but different accession
#'  numbers will be treated as one group.
#'@param mc_psm_by A character string specifying the method in the median
#'  centering of PSM \code{log2FC} across samples. At the \code{peptide}
#'  default, the median description of PSMs (grouped by \code{pep_seq} or
#'  \code{pep_seq_mod} according to \code{group_psm_by}) will be first
#'  calculated and the offsets to zero (of logarithmic 2) will be used for the
#'  centering of PSMs across samples. At \code{mc_psm_by = protein}, the median
#'  description of PSMs (grouped by \code{prot_acc} or \code{gene} according to
#'  \code{group_pep_by}) will be calculated and the corresponding offsets to
#'  zero will be applied. At the \code{mc_psm_by = psm}, PSMs will be median
#'  centered without grouping.
#'@param plot_log2FC_cv Logical; if TRUE, the distributions of the CV of peptide
#'  \code{log2FC} will be plotted. The default is TRUE.
#'@param ... Not currently used.
#'@inheritParams load_expts
#'@inheritParams splitPSM
#'@import dplyr tidyr purrr ggplot2 RColorBrewer
#'@importFrom stringr str_split
#'@importFrom tidyr gather
#'@importFrom magrittr %>% %T>% %$% %<>% 
annotPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                     mc_psm_by = "peptide", 
                     fasta = NULL, expt_smry = "expt_smry.xlsx", 
                     plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, ...) {

  dat_dir <- get_gl_dat_dir()
  
  hd_fn <- list.files(path = file.path(dat_dir, "PSM/cache"),
                      pattern = "^F\\d+_header.txt$")
  assign("df_header", readLines(file.path(dat_dir, "PSM/cache", hd_fn[1])))
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM/cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files()
  
  set_indexes <- gsub("^TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort() 
  
  for (set_idx in set_indexes) {
    sublist <- filelist[grep(paste0("^TMTset", set_idx, "_"), 
                             filelist, 
                             ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename = 
                           do.call('rbind', strsplit(as.character(sublist), 
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    # --- LCMS injections under the same TMT experiment ---
    for (idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM/cache", sublist[idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#")

      if (TMT_plex > 0) {
        df <- local({
          injn_idx <- sublist[idx] %>% 
            gsub("^TMTset\\d+_LCMSinj(\\d+)_Clean.txt$", "\\1", .) %>% 
            as.integer()
          
          df <- mcPSM(df, set_idx, injn_idx, mc_psm_by, group_psm_by, group_pep_by)
        })
      }
      
      df <- df %>% 
        dplyr::mutate(!!group_psm_by := as.character(!!rlang::sym(group_psm_by)))
      
      df <- df %>% 
        calcSD_Splex(id = group_psm_by, type = "log2_R") %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
        dplyr::right_join(df, by = group_psm_by) %>% 
        dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 0)) %>% 
        dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 4)) %>% 
        na_genes_by_acc() %>% 
        reloc_col_before("pep_seq", "pep_res_after") %>% 
        reloc_col_after("pep_seq_mod", "pep_seq")

      pep_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by)) %>%
        dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
        dplyr::summarise(pep_n_psm = n())
      
      prot_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_psm = n())
      
      prot_n_pep <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_pep = n())
      
      df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)

      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(dplyr::left_join, by = group_pep_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))), 
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
      ) 

      df <- dplyr::bind_cols(
        df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^I[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))),
      ) %>% 
        dplyr::select(-which(names(.) %in% c("pep_index", "prot_index"))) %T>% 
        write.table(file.path(dat_dir, "PSM", paste0(out_fn[idx, 1], ".txt")),
                    sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int && TMT_plex > 0) {
        df_int <- df %>% .[, grepl("^N_I[0-9]{3}", names(.))]
        rptr_violin(df = df_int, 
                    filepath = 
                      file.path(dat_dir, "PSM/rprt_int/mc", 
                                paste0(gsub("_PSM_N", "", out_fn[idx, 1]), 
                                       "_rprt.png")), 
                    width = 8, height = 8)
      }

      if (plot_log2FC_cv && TMT_plex > 0) {
        sd_violin(df = df, id = !!group_psm_by, 
                  filepath = 
                    file.path(dat_dir, "PSM/log2FC_cv/raw", 
                              paste0(gsub("_PSM_N", "", out_fn[idx, 1]), "_sd.png")), 
                  width = 8, height = 8, type = "log2_R", 
                  adjSD = FALSE, is_psm = TRUE, ...)
      }
    }
  }
}


#'Standardization of PSM
#'
#'\code{normPSM} standardizes
#'\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSM}
#'results from database search engines.
#'
#'In each primary output file, "\code{...PSM_N.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the average
#'intensity of \code{reference(s)} within each multiplex TMT set, or to the
#'row-mean intensity within each plex if no \code{reference(s)} are present.
#'Values under columns \code{N_log2_R...} are \code{log2_R...} with
#'median-centering alignment. Values under columns \code{I...} are raw
#'\code{reporter-ion intensity} from database searches. Values under columns
#'\code{N_I...} are normalized \code{reporter-ion intensity}. Values under
#'columns \code{sd_log2_R...} are the standard deviation of the \code{log2FC} of
#'peptides from ascribing PSMs. Character strings under \code{pep_seq_mod}
#'denote peptide sequences with applicable variable modifications.
#'
#'\cr \strong{Nomenclature of \code{pep_seq_mod}}:
#'
#'\tabular{ll}{ \emph{Variable modification}   \tab \emph{Abbreviation}\cr
#'Non-terminal \tab A letter from upper to lower case and the flanking residues
#'on the N- or C-terminal side of the peptide separated by a dot, e.g.,
#'\code{-.mtFPEADILLK.S} \cr N-term \tab A hat to the left of a peptide
#'sequence, e.g., \code{K.^QDGTHVVEAVDATHIGK.L} \cr C-term \tab A hat to the
#'right of a peptide sequence, e.g., \code{K.DAYYNLCLPQRPnMI^.-} \cr Acetyl
#'(Protein N-term) \tab A underscore to the left of a peptide sequence, e.g.,
#'\code{-._mAsGVAVSDGVIK.V}. \cr Amidated (Protein C-term) \tab A underscore to
#'the right of a peptide sequence, e.g., \code{K.DAYYNLCLPQRPnMI_.-}. \cr Other
#'(Protein N-term) \tab A tilde to the left of a peptide sequence, e.g.,
#'\code{-.~mAsGVAVSDGVIK.V} \cr Other (Protein C-term) \tab An tilde to the
#'right of a peptide sequence, e.g. \code{K.DAYYNLCLPQRPnMI~.-} \cr }
#'
#'@section \code{Mascot}: End users will export \code{PSM} data from
#'  \href{https://http://www.matrixscience.com/}{Mascot} at a \code{.csv}
#'  format and store them under the file folder indicated by \code{dat_dir}. The
#'  header information should be included during the \code{.csv} export. The
#'  format of file name(s) should be defaulted by
#'  \href{https://http://www.matrixscience.com/}{Mascot}: starting with
#'  the letter \code{'F'}, followed by digits without space and ended with a
#'  \code{'.csv'} extension \code{(e.g., F004453.csv)}.
#'
#'@section \code{MaxQuant}: End users will copy over \code{msms.txt} file(s)
#'  from \href{https://www.maxquant.org/}{MaxQuant} to the \code{dat_dir}
#'  directory. In the case of multiple \code{msms.txt} files for processing, the
#'  file names need to be compiled in that they all start with \code{'msms'} and
#'  end with a \code{'.txt'} extension.
#'
#'@section \code{Spectrum Mill}: End users will copy over \code{PSMexport.1.ssv}
#'  file(s) from
#'  \href{https://www.agilent.com/en/products/software-informatics/masshunter-suite/masshunter-for-life-science-research/spectrum-mill}{Spectrum Mill} 
#'  to the \code{dat_dir} directory. In the case of multiple
#'  \code{PSMexport} files for processing, the file names need to be compiled in
#'  that they all start with \code{'PSMexport'} and end with a \code{'.ssv'}
#'  extension.
#'
#'@param ... \code{filter_}: Variable argument statements for the filtration of
#'  data rows. Each statement contains to a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_expect} is a column key present in \code{Mascot} PSM
#'  exports and \code{filter_psms_at = exprs(pep_expect <= 0.1)} will remove PSM
#'  entries with \code{pep_expect > 0.1}.
#'@inheritParams load_expts
#'@inheritParams splitPSM
#'@inheritParams splitPSM_mq
#'@inheritParams cleanupPSM
#'@inheritParams annotPSM
#'@seealso 
#'  \emph{Metadata} \cr 
#'  \code{\link{load_expts}} for metadata preparation and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr 
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
#'  \code{\link{mergePep}} for extended examples in peptide data merging \cr 
#'  \code{\link{standPep}} for extended examples in peptide data normalization \cr
#'  \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
#'  \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
#'  
#'  \emph{Variable arguments of `filter_...`} \cr 
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \emph{Missing values} \cr 
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  
#'  \emph{Informatics} \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA visualization \cr 
#'  \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@section \code{Variable arguments and data files}: Variable argument (vararg)
#'  statements of \code{filter_} and \code{arrange_} are available in
#'  \code{proteoQ} for flexible filtration and ordering of data rows, via
#'  functions at users' interface. To take advantage of the feature, users need
#'  to be aware of the column keys in input files. As indicated by their names,
#'  \code{filter_} and \code{filter2_} perform row filtration against column
#'  keys from a primary data file, \code{df}, and secondary data file(s),
#'  \code{df2}, respectively. The same correspondence is applicable for
#'  \code{arrange_} and \code{arrange2_} varargs. \cr \cr Users will typically
#'  employ either primary or secondary vararg statements, but not both. In the
#'  more extreme case of \code{gspaMap(...)}, it links \code{\link{prnGSPA}}
#'  findings in \code{df2} to the significance \code{pVals} and abundance fold
#'  changes in \code{df} for volcano plot visualizations by gene sets. The table
#'  below summarizes the \code{df} and the \code{df2} for varargs in
#'  \code{proteoQ}.
#'
#'  \tabular{lllll}{ 
#'  \strong{Utility} \tab \strong{Vararg_} \tab \strong{df} \tab \strong{Vararg2_} \tab \strong{df2} \cr
#'  normPSM \tab filter_ \tab Mascot, \code{F[...].csv}; MaxQuant, \code{msms[...].txt}; 
#'  SM, \code{PSMexport[...].ssv} \tab NA \tab NA \cr 
#'  PSM2Pep \tab NA \tab NA \tab NA \tab NA \cr 
#'  mergePep \tab filter_ \tab \code{TMTset1_LCMSinj1_Peptide_N.txt} \tab NA \tab NA \cr 
#'  standPep \tab slice_ \tab \code{Peptide.txt} \tab NA \tab NA \cr 
#'  Pep2Prn \tab filter_ \tab \code{Peptide.txt} \tab NA \tab NA \cr 
#'  standPrn \tab slice_\tab \code{Protein.txt} \tab NA \tab NA \cr 
#'  pepHist \tab filter_\tab \code{Peptide.tx}t \tab NA \tab NA \cr 
#'  prnHist \tab filter_\tab \code{Protein.txt} \tab NA \tab NA \cr 
#'  pepSig \tab filter_\tab \code{Peptide[_impNA].txt} \tab NA \tab NA \cr 
#'  prnSig \tab filter_\tab \code{Protein[_impNA].txt} \tab NA \tab NA \cr 
#'  pepMDS \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnMDS \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepPCA \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnPCA \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepLDA \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnLDA \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepEucDist \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnEucDist \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepCorr_logFC \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnCorr_logFC \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepHM \tab filter_, arrange_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnHM \tab filter_, arrange_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  
#'  anal_prnTrend \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  plot_prnTrend \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_Trend_{NZ}[_impNA][...].txt} \cr 
#'  
#'  anal_pepNMF \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  anal_prnNMF \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  plot_pepNMFCon \tab NA \tab NA \tab filter2_\tab \code{[...]Peptide_NMF[...]_consensus.txt} \cr 
#'  plot_prnNMFCon \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_NMF[...]_consensus.txt} \cr 
#'  plot_pepNMFCoef \tab NA \tab NA \tab filter2_\tab \code{[...]Peptide_NMF[...]_coef.txt} \cr 
#'  plot_prnNMFCoef \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_NMF[...]_coef.txt} \cr 
#'  plot_metaNMF \tab filter_, arrange_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  
#'  prnGSPA \tab filter_\tab \code{Protein[_impNA]_pVals.txt} \tab NA \tab NA \cr 
#'  prnGSPAHM \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_GSPA_{NZ}[_impNA]_essmap.txt} \cr 
#'  gspaMap \tab filter_\tab \code{Protein[_impNA]_pVal.txt} \tab filter2_\tab \code{[...]Protein_GSPA_{NZ}[_impNA].txt} \cr 
#'  
#'  anal_prnString \tab filter_\tab \code{Protein[_impNA][_pVals].tx}t \tab NA \tab NA \cr 
#'  }
#'
#'@return Outputs are under the directory of \code{PSM} sub to \code{dat_dir}.
#'  Primary results are in \code{TMTset1_LCMSinj1_PSM_N.txt,
#'  TMTset2_LCMSinj1_PSM_N.txt, ...} The indexes of TMT experiment and LC/MS
#'  injection are indicated in the file names.
#'@example inst/extdata/examples/normPSM_.R
#'@import dplyr purrr ggplot2 RColorBrewer
#'@importFrom stringr str_split
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
normPSM <- function(group_psm_by = c("pep_seq", "pep_seq_mod"), 
                    group_pep_by = c("prot_acc", "gene"), 
                    mc_psm_by = c("peptide", "protein", "psm"), 
                    dat_dir = NULL, expt_smry = "expt_smry.xlsx", 
                    frac_smry = "frac_smry.xlsx", 
                    fasta = NULL, entrez = NULL, 
                    pep_unique_by = c("group", "protein", "none"), 
                    corrected_int = TRUE, rm_reverses = TRUE, 
                    rptr_intco = 0, rm_craps = FALSE, 
                    rm_krts = FALSE, rm_outliers = FALSE, 
                    annot_kinases = FALSE, plot_rptr_int = TRUE, 
                    plot_log2FC_cv = TRUE, 
                    use_lowercase_aa = TRUE, purge_phosphodata = TRUE, 
                    parallel = TRUE, ...) {
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
            c(dots) %>% 
            save_call("normPSM"), 
          add = TRUE)
  
  dots <- rlang::enexprs(...)

  if (is.null(dat_dir)) {
    dat_dir <- get_gl_dat_dir()
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
    message("Variable `dat_dir` added to the Global Environment.")
  }

  if (is.null(fasta)) {
    stop("Path(s) to fasta file(s) cannot be empty.", 
         call. = FALSE)
  }
  
  group_psm_by <- rlang::enexpr(group_psm_by)
  if (group_psm_by == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
    group_psm_by <- "pep_seq"
  } else {
    group_psm_by <- rlang::as_string(group_psm_by)
    stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"), 
              length(group_psm_by) == 1)
  }
  
  group_pep_by <- rlang::enexpr(group_pep_by)
  if (group_pep_by == rlang::expr(c("prot_acc", "gene"))) {
    group_pep_by <- "prot_acc"
  } else {
    group_pep_by <- rlang::as_string(group_pep_by)
    stopifnot(group_pep_by %in% c("prot_acc", "gene"), length(group_pep_by) == 1)
  }

  dir.create(file.path(dat_dir, "PSM/cache"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/rprt_int/raw"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/rprt_int/mc"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/log2FC_cv/raw"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/log2FC_cv/purged"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/individual_mods"), 
             recursive = TRUE, showWarnings = FALSE)
  
  if (!purrr::is_empty(list.files(path = file.path(dat_dir), 
                                  pattern = "^F[0-9]+\\.csv$"))) {
    type <- "mascot"
    message("Mascot results found.")
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), 
                                         pattern = "^msms.*\\.txt$"))) {
    type <- "mq"
    message("MaxQuant results found.")
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), 
                                         pattern = "^PSMexport.*\\.ssv$"))) {
    type <- "sm"
    message("Spectrum Mill results found.")
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), 
                                         pattern = "^psm.*\\.tsv$"))) {
    type <- "mf"
    message("MSFragger results found.")
  } else {
    stop("Unknow data type or missing data files.", 
         call. = FALSE)
	}

  pep_unique_by <- rlang::enexpr(pep_unique_by)
  if (pep_unique_by == rlang::expr(c("group", "protein", "none"))) {
    pep_unique_by <- "group"
  } else {
    pep_unique_by <- rlang::as_string(pep_unique_by)
    stopifnot(pep_unique_by %in% c("group", "protein", "none"), 
              length(pep_unique_by) == 1)
  }
  
  mc_psm_by <- rlang::enexpr(mc_psm_by)
  if (mc_psm_by == rlang::expr(c("peptide", "protein", "psm"))) {
    mc_psm_by <- "peptide"
  } else {
    mc_psm_by <- rlang::as_string(mc_psm_by)
    stopifnot(mc_psm_by %in% c("peptide", "protein", "psm"), 
              length(mc_psm_by) == 1)
  }

  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))
  
  stopifnot(vapply(c(rptr_intco), is.numeric, logical(1)))
  stopifnot(vapply(c(corrected_int, 
                     rm_reverses, 
                     rm_craps, 
                     rm_krts, 
                     rm_outliers, 
                     annot_kinases, 
                     plot_rptr_int, 
                     plot_log2FC_cv, 
                     use_lowercase_aa, 
                     purge_phosphodata), 
                   rlang::is_logical, logical(1)))
  
  reload_expts()

  if (type == "mascot") {
    rmPSMHeaders()
    splitPSM(group_psm_by = group_psm_by, 
             group_pep_by = group_pep_by, 
             fasta = fasta, 
             entrez = entrez, 
             pep_unique_by = pep_unique_by, 
             rm_craps = rm_craps, 
             rm_krts = rm_krts, 
             purge_phosphodata = purge_phosphodata, 
             annot_kinases = annot_kinases, 
             plot_rptr_int = plot_rptr_int, 
             rptr_intco = rptr_intco, 
             use_lowercase_aa = use_lowercase_aa, 
             parallel = parallel, ...)
    cleanupPSM(rm_outliers = rm_outliers, 
               group_psm_by = group_psm_by, 
               parallel = parallel)
    annotPSM(group_psm_by = group_psm_by, 
             group_pep_by = group_pep_by, 
             mc_psm_by = mc_psm_by, 
             fasta = fasta, 
             expt_smry = expt_smry, 
             plot_rptr_int = plot_rptr_int, 
             plot_log2FC_cv = plot_log2FC_cv, ...)
  } else if (type == "mq") {
    splitPSM_mq(group_psm_by = group_psm_by, 
                group_pep_by = group_pep_by, 
                fasta = fasta, 
                entrez = entrez, 
                pep_unique_by = pep_unique_by, 
                corrected_int = corrected_int, 
                rm_craps = rm_craps, 
                rm_krts = rm_krts, 
                rm_reverses = rm_reverses, 
                purge_phosphodata = purge_phosphodata, 
                annot_kinases = annot_kinases, 
                plot_rptr_int = plot_rptr_int, 
                rptr_intco = rptr_intco, 
                use_lowercase_aa = use_lowercase_aa, 
                parallel = parallel, ...)
    cleanupPSM(rm_outliers = rm_outliers, 
               group_psm_by = group_psm_by, 
               parallel = parallel)
		annotPSM_mq(group_psm_by = group_psm_by, 
		            group_pep_by = group_pep_by, 
		            mc_psm_by = mc_psm_by, 
		            fasta = fasta, 
		            expt_smry = expt_smry, 
		            plot_rptr_int = plot_rptr_int, 
		            plot_log2FC_cv = plot_log2FC_cv, 
		            ...)
  } else if (type == "sm") {
    splitPSM_sm(group_psm_by = group_psm_by, 
                group_pep_by = group_pep_by, 
                fasta = fasta, 
                entrez = entrez, 
                pep_unique_by = pep_unique_by, 
                rm_craps = rm_craps, 
                rm_krts = rm_krts, 
                purge_phosphodata = purge_phosphodata, 
                annot_kinases = annot_kinases, 
                plot_rptr_int = plot_rptr_int, 
                rptr_intco = rptr_intco, 
                use_lowercase_aa = use_lowercase_aa, 
                parallel = parallel, ...)
    cleanupPSM(rm_outliers = rm_outliers, 
               group_psm_by = group_psm_by, 
               parallel = parallel)
    annotPSM_sm(group_psm_by = group_psm_by, 
                group_pep_by = group_pep_by, 
                mc_psm_by = mc_psm_by, 
                fasta = fasta, 
                expt_smry = expt_smry, 
                plot_rptr_int = plot_rptr_int, 
                plot_log2FC_cv = plot_log2FC_cv, 
                ...)
  } else if (type == "mf") {
    splitPSM_mf(group_psm_by = group_psm_by, 
                group_pep_by = group_pep_by, 
                fasta = fasta, 
                entrez = entrez, 
                pep_unique_by = pep_unique_by, 
                rm_craps = rm_craps, 
                rm_krts = rm_krts,
                purge_phosphodata = purge_phosphodata, 
                annot_kinases = annot_kinases, 
                plot_rptr_int = plot_rptr_int, 
                rptr_intco = rptr_intco, 
                use_lowercase_aa = use_lowercase_aa, 
                parallel = parallel, 
                ...)
    cleanupPSM(rm_outliers = rm_outliers, 
               group_psm_by = group_psm_by, 
               parallel = parallel)
    annotPSM_mf(group_psm_by = group_psm_by, 
                group_pep_by = group_pep_by, 
                mc_psm_by = mc_psm_by, 
                fasta = fasta, 
                expt_smry = expt_smry, 
                plot_rptr_int = plot_rptr_int, 
                plot_log2FC_cv = plot_log2FC_cv, 
                ...)
  } 
}


#' Remove columns from PSM table
#' 
#' @param df PSM data
#' @inheritParams annotPSM
#' @inheritParams channelInfo
rm_cols_mqpsm <- function(df, group_psm_by, set_idx) {
  df <- local({
    uniq_by <- c("pep_seq_mod", "Charge", "Retention time", "raw_file")
    
    df %>% 
      tidyr::unite(uniq_id, uniq_by, sep = "@", remove = FALSE) %>% 
      dplyr::group_by(uniq_id) %>% 
      dplyr::arrange(-pep_expect) %>% 
      dplyr::mutate(.n = row_number()) %>% 
      dplyr::filter(.n == 1) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-uniq_id, -.n) %>% 
      dplyr::mutate(TMT_Set = set_idx)
  })
  
  df <- local({
    if ("Proteins" %in% names(df)) {
      if ("TMT_Set" %in% names(df)) {
        df <- df %>% reloc_col_before("Proteins", "TMT_Set")
      } else {
        df <- df %>% reloc_col_after_last("Proteins")
      }
    }
    
    if ("Protein group IDs" %in% names(df)) {
      if ("TMT_Set" %in% names(df)) {
        df <- df %>% reloc_col_before("Protein group IDs", "TMT_Set")
      } else {
        df <- df %>% reloc_col_after_last("Protein group IDs")
      }
    }
    
    col_start <- which(names(df) == "Scan number")
    col_end <- which(names(df) == "Retention time")
    
    if (!(purrr::is_empty(col_start) || purrr::is_empty(col_end))) {
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
      "ID", # "Protein group IDs", 
      "Peptide ID", "Mod. peptide ID", "Evidence ID", 
      "Length", 
      "Precursor full scan number", "Precursor apex fraction", 
      "Precursor apex offset", "Precursor apex offset time", 
      "Number of matches"))) %>% 
    dplyr::select(-grep("site IDs$", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("m/z", "Charge", 
                                         "Mass", "Mass error [ppm]", "Mass error [Da]", 
                                         "Combinatorics", "PIF", "Fraction of total spectrum", 
                                         "Base peak fraction", "Precursor Intensity", 
                                         "Precursor Apex Fraction", 
                                         "Intensity coverage", "Peak coverage"))) %>% 
    dplyr::select(-grep("^Reporter mass deviation ", names(.))) 
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c("raw_file")))
  
  if (group_psm_by == "pep_seq_mod") {
    df <- df %>% dplyr::select(-pep_seq)
  } else {
    df <- df %>% dplyr::select(-pep_seq_mod)
  }
}


#' Calculate peptide data for individual TMT experiments
#' 
#' Argument \code{injn_idx} not currently used.
#' 
#' @inheritParams mcPSM
#' @inheritParams PSM2Pep
#' @inheritParams annotPSM
#' @inheritParams channelInfo
#' @inheritParams locate_outliers
#' @inheritParams TMT_levels
calcPeptide <- function(df, group_psm_by, method_psm_pep, 
                       group_pep_by, dat_dir, set_idx, injn_idx, TMT_plex) {
  
  stopifnot("prot_acc" %in% names(df))

  channelInfo <- channelInfo(dat_dir = dat_dir, set_idx = set_idx, injn_idx = injn_idx)
  
  if (TMT_plex > 0) {
    df <- df[rowSums(!is.na(df[, grepl("^N_log2_R[0-9]{3}", names(df)), drop = FALSE])) > 0, ]
  } 
  
  df <- df %>%
    dplyr::arrange(!!rlang::sym(group_psm_by), prot_acc) %>%
    dplyr::select(-grep("^R[0-9]{3}", names(.)))
  
  # --- Mascot ---
  is_mascot <- if ("pep_var_mod" %in% names(df)) TRUE else FALSE
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "prot_hit_num", "prot_family_member", "prot_score", 
      "prot_matches", "prot_sequences", 
      "pep_var_mod", "pep_var_mod_pos", "pep_scan_title", 
      "raw_file", "pep_query", "pep_summed_mod_pos", "pep_local_mod_pos", 
      "pep_rank", "pep_isbold", "pep_exp_mz", "pep_exp_mr", "pep_exp_z", 
      "pep_delta", "pep_calc_mr")))
  
  if ("pep_scan_title" %in% names(df)) {
    df <- df %>% 
      dplyr::mutate(pep_scan_title = gsub("\\\\", "~~", pep_scan_title)) %>%
      dplyr::mutate(pep_scan_title = gsub("^File.*~~", "", pep_scan_title))
  }
  
  # --- MaxQuant ---
  is_mq <- if (any(c("Scan number", "Scan Number") %in% names(df))) TRUE else FALSE
  
  df <- local({
    col_start <- which(names(df) == "Scan number")
    col_end <- which(names(df) == "Retention time")
    
    if (!(purrr::is_empty(col_start) || purrr::is_empty(col_end))) {
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
      "ID", # "Protein group IDs", 
      "Peptide ID", "Mod. peptide ID", "Evidence ID", 
      "Length"))) %>% 
    dplyr::select(-grep("site IDs$", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("m/z", "Charge", "Mass", "Mass error [ppm]", 
                                         "Mass error [Da]", 
                                         "Combinatorics", "PIF", "Fraction of total spectrum", 
                                         "Base peak fraction", "Precursor Intensity", 
                                         "Precursor Apex Fraction", 
                                         "Intensity coverage", "Peak coverage"))) %>% 
    dplyr::select(-grep("^Reporter mass deviation ", names(.)))
  
  # --- Spectrum Mill ---
  is_sm <- if ("matched_parent_mass" %in% names(df)) TRUE else FALSE
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "number", "modifications", 
      "variableSites", "nterm", "previous_aa", "sequence", "next_aa", 
      "cys", "searchCycle", "L/H", "accession_numbers", "entry_name", 
      "matched_parent_mass", "parent_charge", 
      "deltaForwardReverseScore", "percent_scored_peak_intensity", 
      "totalIntensity", "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
      "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
      "delta_parent_mass", "delta_parent_mass_ppm"))) 
  
  # --- MSFragger ---
  is_mf <- if ("Number of Enzymatic Termini" %in% names(df)) TRUE else FALSE
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "Spectrum", "Spectrum File", "Ion Mobility", "raw_file", "Peptide Length", 
      "Charge", "Retention", "Observed Mass", "Calibrated Observed Mass", 
      "Observed M/Z", "Calibrated Observed M/Z", "Calculated Peptide Mass", 
      "Calculated M/Z", "Delta Mass", "Number of Missed Cleavages", 
      "Number of Enzymatic Termini", 
      "Ion Mobility", "Assigned Modifications", "Observed Modifications", 
      "Entry Name", "Protein Description")))
  
  
  # summarizes log2FC and intensity from the same `set_idx` at one or multiple LCMS series
  df_num <- switch(method_psm_pep, 
                   mean = 
                     aggrNums(mean)(df, !!rlang::sym(group_psm_by), na.rm = TRUE), 
                   median = 
                     aggrNums(median)(df, !!rlang::sym(group_psm_by), na.rm = TRUE), 
                   weighted_mean = 
                     tmt_wtmean(df, !!rlang::sym(group_psm_by), na.rm = TRUE), 
                   top_3_mean = 
                     TMT_top_n(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   lfq_max = 
                     aggrTopn(sum)(df, !!rlang::sym(group_psm_by), 1, na.rm = TRUE), 
                   lfq_top_2_sum = 
                     aggrTopn(sum)(df, !!rlang::sym(group_psm_by), 2, na.rm = TRUE), 
                   lfq_top_3_sum = 
                     aggrTopn(sum)(df, !!rlang::sym(group_psm_by), 3, na.rm = TRUE), 
                   lfq_all = 
                     aggrTopn(sum)(df, !!rlang::sym(group_psm_by), Inf, na.rm = TRUE), 
                   aggrNums(median)(df, !!rlang::sym(group_psm_by), na.rm = TRUE))
  
  if (TMT_plex == 0) {
    df_num <- df_num %>% 
      dplyr::mutate(log2_R000 = NA, N_log2_R000 = NA)
  }
  
  df_first <- df %>% 
    dplyr::select(-which(names(.) %in% c("pep_tot_int", 
                                         "pep_unique_int", 
                                         "pep_razor_int"))) %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    med_summarise_keys(group_psm_by) 
  
  df <- list(df_first, df_num) %>%
    purrr::reduce(left_join, by = group_psm_by) %>% 
    reloc_col_after("pep_score", "pep_razor_unique") %>% 
    reloc_col_after("pep_expect", "pep_score") 

  if (is_mf) {
    df <- local({
      percent_na <- sum(is.na(df$gene))/length(unique(df$gene))
      
      if (group_pep_by == "gene" && percent_na < .5) {
        col_map <- "Mapped Genes"
      } else {
        col_map <- "Mapped Proteins"
      }
      
      df <- df %>% 
        dplyr::mutate(pep_tot_int = I000, 
                      pep_unique_int = ifelse(`Is Unique` & is.na(.[[col_map]]), I000, 0), 
                      pep_razor_int = ifelse(`Is Unique`, I000, 0))
    })
  } else {
    # df <- df %>% dplyr::mutate(pep_tot_int = NA, pep_unique_int = NA, pep_razor_int = NA)
  }
  
  if (group_psm_by == "pep_seq_mod") {
    df <- df %>% dplyr::select(-pep_seq)
  } else {
    df <- df %>% dplyr::select(-pep_seq_mod)
  }
  
  df <- cbind.data.frame(df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df)), drop = FALSE],
                         df[, grepl("I[0-9]{3}", names(df)), drop = FALSE], 
                         df[, grepl("log2_R[0-9]{3}", names(df)), drop = FALSE]) %>%
    dplyr::mutate_at(.vars = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)),
                     list(~ replace(.x, is.infinite(.x), NA)))
  
  if (TMT_plex > 0) {
    df <- local({
      col_r <- grepl("^log2_R[0-9]{3}", names(df))
      cf <- apply(df[, col_r, drop = FALSE], 2, median, na.rm = TRUE)
      
      df <- cbind(df[, -grep("^N_log2_R[0-9]{3}", names(df))],
                  sweep(df[, col_r, drop = FALSE], 2, cf, "-") %>%
                    `colnames<-`(paste("N", colnames(.), sep="_")))
      
      col_int <- grepl("^I[0-9]{3}", names(df))
      cbind(df[, -grep("^N_I[0-9]{3}", names(df))],
                   sweep(df[, col_int, drop = FALSE], 2, 2^cf, "/") %>%
                     `colnames<-`(paste("N", colnames(.), sep = "_")))
    })
  }

  df <- df %>% 
    dplyr::mutate(!!group_pep_by := as.character(!!rlang::sym(group_pep_by)))
  
  df <- df %>% 
    calcSD_Splex(group_pep_by) %>% 
    `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
    dplyr::right_join(df, by = group_pep_by) %>% 
    na_zeroIntensity() %>% 
    dplyr::mutate(TMT_Set = set_idx)
  
  return(df)
}


#' Helper of PSM2Pep
#' 
#' @param file The name of a PSM file.
#' @inheritParams load_expts
#' @inheritParams calcPeptide
#' @importFrom magrittr %>% %T>% %$% %<>% 
psm_to_pep <- function (file, dat_dir, label_scheme_full, 
                        group_psm_by, group_pep_by, method_psm_pep) {
  fn_prx <- gsub("_PSM_N.txt", "", file, fixed = TRUE)
  set_idx <- as.integer(gsub(".*TMTset(\\d+)_.*", "\\1", fn_prx))
  injn_idx <- as.integer(gsub(".*LCMSinj(\\d+).*", "\\1", fn_prx))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  df <- read.csv(file.path(dat_dir, "PSM", file), check.names = FALSE, header = TRUE,
                 sep = "\t", comment.char = "#") %>% 
    dplyr::select(-grep("^sd_log2_R", names(.))) 
  
  # special handling for MaxQuant; no `Precursor Intensity` for .d files
  is_mq_lfq <- local({
    # inconsistent upper or lower cases
    cols <- c("Retention time", "Retention Time")
    is_mq <- if (any(cols %in% names(df))) TRUE else FALSE
    is_mq && (TMT_plex == 0)
  })

  if (is_mq_lfq) {
    message("Precursor intensity/area not available at PSM levels.\n",  
            "(LFQ) intensity from MaxQuant `peptides[...].txt` will be used.")
    
    df <- rm_cols_mqpsm(df, group_psm_by, set_idx)
  } else {
    df <- df %>% 
        calcPeptide(group_psm_by, method_psm_pep, group_pep_by, dat_dir, 
                    set_idx, injn_idx, TMT_plex)
  }
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]{0,1}", names(.))), 
    df %>% dplyr::select(grep("^I[0-9]{3}[NC]{0,1}", names(.))), 
    df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]{0,1}", names(.))), 
    df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]{0,1}", names(.))), 
    df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))), 
    df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))),
  ) %T>% 
    write.table(file.path(dat_dir, "Peptide", paste0(fn_prx, "_Peptide_N.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
}


#'Interim peptide tables
#'
#'\code{PSM2Pep} summarizes
#'\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSMs} to
#'peptides by individual TMT (or LFQ) experiments and LC/MS series.
#'
#'In general, fields other than \code{log2FC} and \code{intensity} are
#'summarized with median statistics. One exception is with \code{pep_expect} in
#'Mascot or \code{PEP} in MaxQuant where geometric mean is applied.
#'
#'@param method_psm_pep Character string; the method to summarize the
#'  \code{log2FC} and the \code{intensity} of \code{PSMs} by peptide entries.
#'  The descriptive statistics includes \code{c("mean", "median",
#'  "weighted_mean", "top_3_mean", "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum",
#'  "lfq_all")} with \code{median} being the default for TMT and
#'  \code{lfq_top_2_sum} for LFQ. The \code{log10-intensity} of reporter (or
#'  LFQ) ions at the \code{PSMs} levels will be the weight when summarizing
#'  \code{log2FC} with various \code{"top_n"} statistics or
#'  \code{"weighted_mean"}.
#'@inheritParams splitPSM
#'@param ... Not currently used.
#'@seealso 
#'  \emph{Metadata} \cr 
#'  \code{\link{load_expts}} for metadata preparation and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr 
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
#'  \code{\link{mergePep}} for extended examples in peptide data merging \cr 
#'  \code{\link{standPep}} for extended examples in peptide data normalization \cr
#'  \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
#'  \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
#'  
#'  \emph{Variable arguments of `filter_...`} \cr 
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \emph{Missing values} \cr 
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  
#'  \emph{Informatics} \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA visualization \cr 
#'  \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'@family custom database preparation
#'@seealso 
#'  \emph{Custom databases} \cr 
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@return Tables under \code{PSM} folder for each TMT experiment and LC/MS
#'  series: \code{TMTset1_LCMSinj1_PSM_N.txt}, \code{TMTset1_LCMSinj2_PSM_N.txt} etc.
#'
#'@example inst/extdata/examples/PSM2Pep_.R
#'@import stringr dplyr purrr 
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
PSM2Pep <- function (method_psm_pep = c("median", "mean", "weighted_mean", "top_3_mean", 
                                        "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum", 
                                        "lfq_all"), 
                     parallel = TRUE, ...) {
  
  dat_dir <- get_gl_dat_dir()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  old_opts <- options()
  on.exit(options(old_opts), add = TRUE)
  options(warn = 1)
  
  on.exit(mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
            c(dots) %>% save_call("PSM2Pep"), 
          add = TRUE)
  
  dots <- rlang::enexprs(...)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  method_psm_pep <- rlang::enexpr(method_psm_pep)
  if (TMT_plex > 0) {
    if (length(method_psm_pep) > 1) {
      method_psm_pep <- "median"
    } else {
      method_psm_pep <- rlang::as_string(method_psm_pep)
    }
  } else {
    if (length(method_psm_pep) > 1) {
      method_psm_pep <- "lfq_top_2_sum"
    } else {
      method_psm_pep <- rlang::as_string(method_psm_pep)
    }
  }
  
  if (method_psm_pep == "top.3") {
    stop("Method `top.3` depreciated; instead use `top_3_mean`.", 
         call. = FALSE)
  } else if (method_psm_pep == "weighted.mean") {
    stop("Method `weighted.mean` depreciated; instead use `weighted_mean`.", 
         call. = FALSE)
  }
  
  stopifnot(method_psm_pep %in% c("median", "mean", "weighted_mean", "top_3_mean", 
                                  "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum", 
                                  "lfq_all"), 
            length(method_psm_pep) == 1)
  
  dir.create(file.path(dat_dir, "Peptide/cache"), 
             recursive = TRUE, showWarnings = FALSE)
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), 
                         pattern = "*_PSM_N\\.txt$") %>%
    reorder_files()
  
  if (parallel) {
    n_cores <- pmin(parallel::detectCores(), length(filelist))
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(cl, list("add_cols_at"), 
                            envir = env_where("add_cols_at"))
    
    parallel::clusterExport(cl, list("reloc_col_before"), 
                            envir = env_where("reloc_col_before"))
    
    suppressWarnings(
      silent_out <- parallel::clusterApply(
        cl, filelist, psm_to_pep, dat_dir, label_scheme_full, 
        group_psm_by, group_pep_by, method_psm_pep
      )
    )
    
    parallel::stopCluster(cl)
  } else {
    purrr::walk(as.list(filelist), psm_to_pep, 
                dat_dir, label_scheme_full, 
                group_psm_by, group_pep_by, method_psm_pep)
  }
}


#' To lower case
#' 
#' @param x A character string of amino acid sequence.
#' @param ch A tag before the letter to conversion to its lower case.
my_tolower <- function(x, ch = "^") {
  locales <- gregexpr(ch, x) %>% .[[1]] %>% `+`(., 1)
  lowers <- map(locales, ~ substr(x, .x, .x)) %>% tolower()
  
  for (i in seq_along(lowers)) {
    substr(x, locales[i], locales[i]) <- lowers[i]
  }
  
  x <- gsub(ch, "", x)
  
  return(x)
}


#' Add the `pep_seq_mod` field to MaxQuant PSMs
#'
#' @inheritParams splitPSM
#' @inheritParams locate_outliers
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_maxquant_pepseqmod <- function(df, use_lowercase_aa) {
  if (!use_lowercase_aa) {
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, 
                                    sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, 
                                        sep = "."))
  } else {
    # (1) all non-terminal modifications: M(Oxidation (M)) -> m ...
    # assume "(Acetyl (Protein N-term))" goes before peptide N-term modifications etc.
    # note: residues before and after not yet included; only added dots to both ends.
    
    # "_(Acetyl (Protein N-term))(Gln->pyro-Glu)QAAAAQGSN(Deamidation (N))GPVK_" -> 
    # ".(Acetyl (Protein N-term))(Gln->pyro-Glu)QAAAAQGSnGPVK."
    
    # "(Hex(1)HexNAc(1) (ST))", "(Hex(5) HexNAc(4) NeuAc(2) Sodium (N))"
    # No cases like: X(YYY (zzz) xxx) but X(YYY (zzz) xxx (kkk)) ends with "))"
    
    df <- df %>% 
      tidyr::separate("pep_seq_mod", c("nt", "interior", "ct"), sep = "_") %>% 
      dplyr::mutate(interior = gsub("([A-Z]){1}\\([^\\(\\)]*\\)", 
                                    paste0("@", "\\1"), interior)) %>% 
      dplyr::mutate(interior = gsub("([A-Z]{1})\\(.*\\s+?\\(+?.*\\){2}", 
                                    paste0("@", "\\1"), interior)) %>% 
      dplyr::mutate_at(vars("interior"), ~ purrr::map_chr(.x, my_tolower, "@")) %>% 
      tidyr::unite(pep_seq_mod, nt, interior, ct, sep = ".", remove = TRUE)
    
    # (2-1) add "_" to sequences from protein N-terminal acetylation
    # ".(Acetyl (Protein N-term))(Gln->pyro-Glu)QAAAAQGSnGPVK." -> 
    # "._(Gln->pyro-Glu)QAAAAQGSnGPVK."
    df <- local({
      n_ac <- df %>% 
        dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% 
        dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))

      if (nrow(n_ac) > 0) {
        n_ac <- n_ac %>% 
          dplyr::mutate(pep_seq_mod = 
                          gsub("(Acetyl (Protein N-term))", "_", pep_seq_mod, fixed = TRUE))
        df <- rbind(rest, n_ac)
      }
      
      return(df)
    })
    
    # (2-2) add "_" to sequences from protein C-terminal amidation
    # ".AAASNGPVK(Xxx->Yyy)(Amidated (Protein C-term))." -> 
    # ".AAASNGPVK(Xxx->Yyy)_."
    df <- local({
      c_am <- df %>% 
        dplyr::filter(grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% 
        dplyr::filter(!grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      if (nrow(c_am) > 0) {
        c_am <- c_am %>% 
          dplyr::mutate(pep_seq_mod = gsub("(Amidated (Protein C-term))", "_", 
                                           pep_seq_mod, fixed = TRUE))
        df <- rbind(rest, c_am)
      }
      
      return(df)
    })
    
    # (3-1) "~" for "(Protein N-term)" other than acetylation
    # assume no dual (My (Protein N-term)) + "Acetyl": "._(My (Protein N-term))AAASSLTK."    
    # ".(My (Protein N-term))(Gln->pyro-Glu)QAAAAQGSnGPVK." -> ".~(Gln->pyro-Glu)QAAAAQGSnGPVK."

    df <- local({
      n_ac <- df %>% 
        dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      other_n <- df %>% 
        dplyr::filter(grepl("Protein N-term", Modifications, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% 
        dplyr::filter(!grepl("Protein N-term", Modifications, fixed = TRUE))

      if (nrow(other_n) > 0) {
        other_n <- other_n %>% 
          dplyr::mutate(pep_seq_mod = 
                          gsub("^\\.\\(.*\\s+?\\(+?.*\\){2}", ".~", pep_seq_mod))
        df <- rbind(rest, n_ac, other_n)
      }
      
      return(df)
    })
    
    # (3-2) "~" for "(Protein C-term)" other than amidation
    # ".AGALAPGPL(Yyy->Xxx)(Other (Protein C-term))." -> ".AGALAPGPL(Yyy->Xxx)~."
    df <- local({
      c_am <- df %>% 
        dplyr::filter(grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      other_c <- df %>% 
        dplyr::filter(grepl("Protein C-term", Modifications, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% 
        dplyr::filter(!grepl("Protein C-term", Modifications, fixed = TRUE))
      
      if (nrow(other_c) > 0) {
        other_c <- other_c %>% 
          dplyr::mutate(pep_seq_mod = 
                          gsub("^(.*)\\(.*\\s+?\\(+?.*\\){2}\\.$", paste0("\\1", "~."),
                                           pep_seq_mod))
        df <- rbind(rest, c_am, other_c)
      }
      
      return(df)
    })
    
    # (4-1) "^" for peptide "(N-term)" modification
    # "._(Carbamyl (N-term))AAAAGALAPGPLPDLAAR." -> "._^AAAAGALAPGPLPDLAAR."
    df <- local({
      nt <- df %>% 
        dplyr::filter(grepl("(N-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% 
        dplyr::filter(!grepl("(N-term)", Modifications, fixed = TRUE))
      
      if (nrow(nt) > 0) {
        nt <- nt %>% 
          dplyr::mutate(pep_seq_mod = gsub("(^\\.[_~]{0,1})\\(.*\\s+?\\(+?.*\\){2}", 
                                           paste0("\\1", "^"), pep_seq_mod))
        df <- rbind(rest, nt)
      }
      
      return(df)
    })
    
    # (4-2) "^" for peptide "(C-term)" modification
    # ".AAAANLCPGQDR(My (C-term))_." -> ".AAAANLCPGQDR^_."
    df <- local({
      ct <- df %>% 
        dplyr::filter(grepl("(C-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% 
        dplyr::filter(!grepl("(C-term)", Modifications, fixed = TRUE))
      
      if (nrow(ct) > 0) {
        ct <- ct %>% 
          dplyr::mutate(pep_seq_mod = gsub("^(.*)\\(.*\\s+?\\(+?.*\\){2}([_~]{0,1}\\.$)", 
                                           paste0("\\1", "^", "\\2"), pep_seq_mod)) 
        df <- rbind(rest, ct)
      }
      
      return(df)
    })
    
    # (5) remove "." at both ends 
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("\\.", "", pep_seq_mod))

    # (6) other N- or C-terminal modifications better but not named with "N-term" or "C-term": 
    # i.e. (Gln->pyro-Glu)
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      gsub("(^[_~]{0,1})\\([^\\(\\)]*\\)", 
                           paste0("\\1", "^"), 
                           pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = 
                      gsub("(^[_~]{0,1})\\(.*\\s+?\\(+?.*\\){2}", 
                           paste0("\\1", "^"), 
                           pep_seq_mod)) %>% 
      
      dplyr::mutate(pep_seq_mod = 
                      gsub("\\([^\\(\\)]*\\)([_~]{0,1}$)", 
                           paste0("^", "\\1"), 
                           pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = 
                      gsub("\\(.*\\s+?\\(+?.*\\){2}([_~]{0,1}$)", 
                           paste0("^", "\\1"), 
                           pep_seq_mod)) %>% 
      
      dplyr::mutate(pep_seq_mod = 
                      paste(pep_res_before, pep_seq_mod, pep_res_after, 
                            sep = ".")) %>% 
      dplyr::mutate(pep_seq = 
                      paste(pep_res_before, pep_seq, pep_res_after, 
                            sep = "."))
  }
  
  # `pep_seq` may not be in fasta
  df %>% 
    dplyr::mutate(
      pep_seq = gsub("^NA.", "?.", pep_seq), 
      pep_seq = gsub(".NA$", ".?", pep_seq), 
      pep_seq_mod = gsub("^NA.", "?.", pep_seq_mod), 
      pep_seq_mod = gsub(".NA$", ".?", pep_seq_mod)
    )
}


#' Add the `pep_seq_mod` field to MSFragger PSMs
#'
#' @inheritParams splitPSM
#' @inheritParams locate_outliers
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_msfragger_pepseqmod <- function(df, use_lowercase_aa) {
  if (!use_lowercase_aa) {
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, 
                                    sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, 
                                        sep = "."))
  } else {
    # (1) all non-terminal modifications
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("([A-Z]){1}\\[[^\\(\\)]*\\]", 
                                       paste0("@", "\\1"), pep_seq_mod)) %>% 
      dplyr::mutate_at(vars("pep_seq_mod"), ~ purrr::map_chr(.x, my_tolower, "@")) %>% 
      dplyr::mutate(.n = row_number())

    # (2-1) add "_" to sequences from protein N-terminal acetylation
    df <- local({
      df_sub <- df %>% 
        dplyr::filter(pep_start <= 2)
      
      df_rest <- df %>% 
        dplyr::filter(! .n %in% df_sub$.n)
      
      df_sub <- df_sub %>% 
        dplyr::mutate(pep_seq_mod = gsub("^n\\[43\\]", "_", pep_seq_mod))
      
      dplyr::bind_rows(df_sub, df_rest)
    })

    # (2-2) add "_" to sequences from protein C-terminal amidation
    df <- local({
      df_sub <- df %>% dplyr::filter(pep_end == prot_len)
      df_rest <- df %>% dplyr::filter(! .n %in% df_sub$.n)

      df_sub <- df_sub %>% 
        dplyr::mutate(pep_seq_mod = gsub("c\\[17\\]$", "_", pep_seq_mod))
      
      dplyr::bind_rows(df_sub, df_rest) 
    })
    
    # (3-1) "~" for "(Protein N-term)" other than acetylation
    df <- local({
      df_sub <- df %>% dplyr::filter(pep_start <= 2)
      df_rest <- df %>% dplyr::filter(! .n %in% df_sub$.n)
      
      df_sub <- df_sub %>% 
        dplyr::mutate(pep_seq_mod = gsub("^n\\[.*\\]", "~", pep_seq_mod))
      
      dplyr::bind_rows(df_sub, df_rest)
    })

    # (3-2) "~" for "(Protein C-term)" other than amidation
    df <- local({
      df_sub <- df %>% dplyr::filter(pep_end == prot_len)
      df_rest <- df %>% dplyr::filter(! .n %in% df_sub$.n)
      
      df_sub <- df_sub %>% 
        dplyr::mutate(pep_seq_mod = gsub("c\\[.*\\]$", "~", pep_seq_mod))
      
      dplyr::bind_rows(df_sub, df_rest) 
    })
    
    # (4-1) "^" for peptide "(N-term)" modification
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("n\\[.*\\].*?", "^", pep_seq_mod))
    
    # (4-2) "^" for peptide "(C-term)" modification
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("c\\[.*\\]$", "^", pep_seq_mod))
    
    # (5) cleanup
    df <- df %>% dplyr::select(-.n) 
    
    # (6) paste "pep_res_before" and "pep_res_after"
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      paste(pep_res_before, pep_seq_mod, pep_res_after, sep = ".")) %>% 
      dplyr::mutate(pep_seq = 
                      paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) 
  }
  
  # `pep_seq` may not be in fasta
  df %>% 
    dplyr::mutate(
      pep_seq = gsub("^NA.", "?.", pep_seq), 
      pep_seq = gsub(".NA$", ".?", pep_seq), 
      pep_seq_mod = gsub("^NA.", "?.", pep_seq_mod), 
      pep_seq_mod = gsub(".NA$", ".?", pep_seq_mod)
    )
}


#' Pad columns to a placeholder data frame.
#' 
#' @param df The original data frame.
#' @param df2 The data frame to be inserted.
#' @param idx The index of \code{df} column for \code{df2} to be inserted
#'   (after).
add_cols_at <- function (df, df2, idx) {
  stopifnot(idx >= 0)
  
  if (idx == 0) {
    bf <- NULL
  } else {
    bf <- df[, seq_len(idx), drop = FALSE]
  }
  
  if ((idx + 1) <= ncol(df)) {
    af <- df[, (idx + 1) : ncol(df), drop = FALSE]
  } else {
    af <- NULL
  }
  
  dplyr::bind_cols(
    bf,
    df2,
    af,
  )
}


#' Replace columns in the original PSM table.
#' 
#' The column index(es) need to be continuous.
#' 
#' @param df The original data frame.
#' @param df2 The data columns to replace those in \code{df}.
#' @param idxs The sequences of column indexes in \code{df}. Note that
#'   \code{idxs} need to be a continuous sequences.
replace_cols_at <- function (df, df2, idxs) {
  ncol <- ncol(df)
  stopifnot(all(idxs >= 1), all(idxs <= ncol))
  
  idxs <- sort(idxs)
  stopifnot(all.equal(idxs - idxs[1] + 1, seq_along(idxs)))
  
  if (idxs[1] >= 2) {
    bf <- df[, 1:(idxs[1]-1), drop = FALSE]
  } else {
    bf <- NULL
  }
  
  if (idxs[length(idxs)] < ncol(df)) {
    af <- df[, (idxs[length(idxs)]+1):ncol(df), drop = FALSE]
  } else {
    af <- NULL
  }
  
  dplyr::bind_cols(
    bf,
    df2,
    af
  )
}


#' Relocate column "m/z" to be immediately before column "Mass".
#' 
#' Not currently used.
#' 
#' @param df The original data frame.
#' @param from The column to be moved from.
#' @param to The column to which the \code{from} will be moved before.
reloc_col <- function (df, from = "m/z", to = "Mass") {
  df0 <- df
  
  df2 <- suppressWarnings(df %>% dplyr::select(one_of(from)))
  if (ncol(df2) == 0) return(df0)
  
  df <- df %>% dplyr::select(-one_of(from))
  
  idx <- which(names(df) == to) - 1
  if (purrr::is_empty(idx)) return(df0)
  
  df <- add_cols_at(df, df2, idx)
  
  return(df)
}


#' Relocate column "to_move" immediately after column "col_before".
#' 
#' @param df The original data frame.
#' @param to_move The column to be moved.
#' @param col_before The anchor column to which the \code{to_move} will be moved after.
reloc_col_after <- function (df, to_move = "after_anchor", col_before = "anchor") {
  if (!(to_move %in% names(df) && col_before %in% names(df))) return(df)
  
  if (to_move == col_before) return(df)
  
  df2 <- df %>% dplyr::select(one_of(to_move))
  df <- df %>% dplyr::select(-one_of(to_move))
  
  idx <- which(names(df) == col_before)
  df <- add_cols_at(df, df2, idx)
  
  return(df)
}


#' Relocate column "to_move" immediately after the last column.
#' 
#' @inheritParams reloc_col_after
reloc_col_after_last <- function (df, to_move = "after_anchor") {
  col_last <- names(df)[ncol(df)]
  reloc_col_after(df, to_move, col_last)
}


#' Relocate column "to_move" immediately after the first column.
#' 
#' @inheritParams reloc_col_after
reloc_col_after_first <- function (df, to_move = "after_anchor") {
  col_first <- names(df)[1]
  reloc_col_after(df, to_move, col_first)
}


#' Relocate column "to_move" immediately before anchor column "col_after".
#' 
#' The same as \code{reloc_col}.
#' 
#' @param df The original data frame.
#' @param to_move The column to be moved.
#' @param col_after The anchor column to which the \code{to_move} will be moved before.
reloc_col_before <- function (df, to_move = "before_anchor", col_after = "anchor") {
  if (!(to_move %in% names(df) && col_after %in% names(df))) return(df)

  df2 <- df %>% dplyr::select(one_of(to_move))
  df <- df %>% dplyr::select(-one_of(to_move))
  
  idx <- which(names(df) == col_after)
  
  df <- add_cols_at(df, df2, idx - 1)
  
  return(df)
}


#' Relocate column "to_move" immediately before the last column.
#' 
#' @inheritParams reloc_col_after
reloc_col_before_last <- function (df, to_move = "after_anchor") {
  col_last <- names(df)[ncol(df)]
  reloc_col_before(df, to_move, col_last)
}


#' Relocate column "to_move" immediately before the first column.
#' 
#' @inheritParams reloc_col_after
reloc_col_before_first <- function (df, to_move = "after_anchor") {
  col_first <- names(df)[1]
  reloc_col_before(df, to_move, col_first)
}


#' Order PSM columns.
#' 
#' Only for certain columns.
#' @param df A PSM data frame.
order_psm_cols <- function (df, cols) {
  purrr::walk(cols, ~ {
    if (is.null(df[[.x]])) df[[.x]] <- NA
    df <<- df
  }, df)
  
  dplyr::bind_cols(
    df %>% dplyr::select(-cols), 
    df %>% dplyr::select(cols),
  )
}


#' Pad MaxQuant TMT channels to the highest plex. 
#' 
#' @param file A file name of PSM table from MaxQuant export
#' @inheritParams splitPSM
pad_mq_channels <- function(file, fasta, entrez) {
  dat_dir <- get_gl_dat_dir()
  
  base_name <- file %>% gsub("\\.txt$", "", .)
  
  df <- read.csv(file.path(dat_dir, file), 
                 check.names = FALSE, header = TRUE, 
                 sep = "\t", comment.char = "#")
  
  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  
  # --- A patch for inconsistency in MaxQuant msms.txt ---
  # TMT: no `Gene Names` and `Protein Names` columns
  # LFQ: 
  #  (1) UniProt identifiers: with `Gene Names` and `Protein Names` columns
  #  (2) RefSeq identifiers: without `Gene Names` and `Protein Names` columns
  if (!("Gene Names" %in% names(df) && "Protein Names" %in% names(df))) {
    stopifnot("Proteins" %in% names(df))
    
    df$"Gene names" <- df$"Gene Names" <- df$"Protein names" <- df$"Protein Names" <- NULL

    df <- local({
      df <- df %>% 
        dplyr::mutate(prot_acc = gsub(";.*$", "", Proteins))
      
      tempdata <- suppressWarnings(
        df %>% 
          dplyr::select(prot_acc) %>% 
          dplyr::filter(!duplicated(prot_acc)) %>% 
          annotPrn(fasta, entrez) %>% 
          dplyr::select(prot_acc, gene, prot_desc) %>% 
          dplyr::rename(`Gene Names` = "gene", "Protein Names" = "prot_desc")
      )
      
      df <- df %>% 
        dplyr::left_join(tempdata, by = "prot_acc") %>% 
        dplyr::select(-prot_acc) 
      
      col <- which(names(df) == "Proteins")
      
      dplyr::bind_cols(
        df[, 1:col],
        df[, c("Gene Names", "Protein Names")],
        df %>% 
          dplyr::select(-c("Gene Names", "Protein Names")) %>% 
          dplyr::select((col+1):ncol(.))
      )
    })
  }
  
  if (! "PSM_File" %in% names(fraction_scheme)) {
    label_scheme_sub <- local({
      raw_files <- unique(df$`Raw file`) %>% 
        gsub("\\.raw$", "", .) %>% 
        gsub("\\.d$", "", .)
      
      tmt_sets <- fraction_scheme %>% 
        dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File), 
                      RAW_File = gsub("\\.d$", "", RAW_File)) %>% 
        dplyr::filter(RAW_File %in% raw_files, !duplicated(TMT_Set)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(TMT_Set) %>% 
        unlist()
  
      label_scheme_full %>% dplyr::filter(TMT_Set %in% tmt_sets)
    })
  } else {
    label_scheme_sub <- local({
      tmt_sets <- fraction_scheme %>% 
          dplyr::mutate(PSM_File = gsub("\\.txt$", "", PSM_File)) %>% 
          dplyr::filter(PSM_File == base_name, !duplicated(TMT_Set)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(TMT_Set) %>% 
        unlist()
    
      label_scheme_full %>% dplyr::filter(TMT_Set %in% tmt_sets)
    })
  }

  nas <- data.frame(rep(NA, nrow(df)))
  sample_ids <- as.character(label_scheme_sub$Sample_ID)

  str_int1 <- "^Reporter intensity [0-9]+"
  str_int2 <- "^Reporter intensity corrected [0-9]+"
  str_dev <- "^Reporter mass deviation \\[mDa\\] [0-9]+"
  
  df_int <- df %>% dplyr::select(grep(str_int1, names(.)))
  df_int2 <- df %>% dplyr::select(grep(str_int2, names(.)))
  df_dev <- df %>% dplyr::select(grep(str_dev, names(.)))
  
  this_plex <- ncol(df_int)
  TMT_plex <- TMT_plex(label_scheme_full)
  stopifnot(this_plex <= TMT_plex, this_plex >= 0)
  
  # Empty.xxx can be due to either channel padding or removals
  if (this_plex > 0 && this_plex < TMT_plex) {
    if (this_plex == 6) {
      if (TMT_plex == 10) {
        pos <- c(3, 5, 7, 9)
      } else if (TMT_plex == 11) {
        pos <- c(3, 5, 7, 9, 11)
      } else if (TMT_plex == 16) {
        pos <- c(3, 5, 7, 9, 11:16)
      } else {
        stop("TMT_plex not one of c(10, 11, 16).", call. = FALSE)
      }
    } else if (this_plex == 10) {
      if (TMT_plex == 11) {
        pos <- 11
      } else if (TMT_plex == 16) {
        pos <- c(11:16)
      } else {
        stop("TMT_plex not one of c(11, 16).", call. = FALSE)
      }
    } else if (this_plex == 11) {
      if (TMT_plex == 16) {
        pos <- c(12:16)
      } else {
        stop("TMT_plex not one of c(16).", call. = FALSE)
      }
    }
    
    for (idx in seq_along(pos)) {
      df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))
      df_int2 <- suppressMessages(add_cols_at(df_int2, nas, pos[idx] - 1))
      df_dev <- suppressMessages(add_cols_at(df_dev, nas, pos[idx] - 1))
    }
    rm(idx)
    
    len <- length(sample_ids)
    if (ncol(df_int) == len && ncol(df_int2) == len && ncol(df_dev) == len) {
      names(df_int) <- paste("Reporter intensity", seq_len(len))
      df <- replace_cols_at(df, df_int, grep(str_int1, names(df)))
      
      names(df_int2) <- paste("Reporter intensity corrected", seq_len(len))
      df <- replace_cols_at(df, df_int2, grep(str_int2, names(df)))
      
      names(df_dev) <- paste("Reporter mass deviation [mDa]", seq_len(len))
      df <- replace_cols_at(df, df_dev, grep(str_dev, names(df)))
    }
  }

  # the same `raw_file` may be at different `dat_file`s
  df$dat_file <- base_name
  
  invisible(df)
}


#'Splits PSM tables
#'
#'\code{splitPSM_mq} splits the PSM outputs by TMT experiment and LC/MS
#'injection.
#'
#'Different to \code{splitPSM} used in \code{Mascot} processes,
#'\code{pep_seq_mod}, \code{prot_n_psm}, \code{prot_n_pep} and \code{pep_n_psm}
#'are calculated later in \code{annotPSM_mq}. This is suitable mainly because
#'there is no columns like \code{prot_matches_sig} and \code{prot_sequences_sig}
#'need to be updated after data merging in \code{splitPSM_mq}.
#'
#'@param pep_unique_by A character string for annotating the uniqueness of
#'  peptides in \code{MaxQuant} PSMs. At the \code{group} default, the
#'  uniqueness of peptides is by groups with the collapses of same-set or
#'  sub-set proteins. At a more stringent criterion of \code{protein}, the
#'  uniqueness of peptides is by protein entries without grouping. On the other
#'  extreme of choice \code{none}, all peptides are treated as unique. A new
#'  column of \code{pep_isunique} with corresponding logical TRUE or FALSE will
#'  be added to the PSM reports. Note that the choice of \code{none} is only for
#'  convenience, as the same can be achieved by setting \code{use_unique_pep =
#'  FALSE} in \link{Pep2Prn}.
#'@param corrected_int A logical argument for uses with \code{MaxQuant} data. At
#'  the TRUE default, values under columns "Reporter intensity corrected..." in
#'  \code{MaxQuant} PSM results (\code{msms.txt}) will be used. Otherwise,
#'  "Reporter intensity" values without corrections will be used.
#'@param rm_reverses A logical argument for uses with \code{MaxQuant} data. At
#'  the TRUE default, \code{Reverse} entries will be removed.
#'@inheritParams splitPSM
#'@import dplyr tidyr
#'@importFrom stringr str_split
#'@importFrom magrittr %>% %T>% %$% %<>% equals
splitPSM_mq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        fasta = NULL, entrez = NULL, 
                        pep_unique_by = "group", corrected_int = TRUE, 
                        rm_craps = FALSE, rm_krts = FALSE, 
                        rm_reverses = TRUE, purge_phosphodata = TRUE, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, 
                        rptr_intco = 0, use_lowercase_aa = TRUE, 
                        parallel = TRUE, ...) {

  on.exit(message("Split PSMss by sample IDs and LCMS series --- Completed."), 
          add = TRUE)
  
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$")
  
  if (purrr::is_empty(filelist)) {
    stop(paste("No PSM files of `msms[...].txt` under", file.path(dat_dir), ".",
               "\nMake sure that the names of PSM files start with `msms`."), 
         call. = FALSE)
  }

  df <- purrr::map(filelist, pad_mq_channels, fasta, entrez)
  df <- suppressWarnings(dplyr::bind_rows(df))
  # nms_mq <- names(df) %>% .[! . %in% c("dat_file")]

  df <- local({
    stopifnot(c("Protein group IDs", "Proteins") %in% names(df))
    
    df <- df %>% 
      dplyr::mutate(pep_literal_unique = ifelse(grepl(";", Proteins), FALSE, TRUE), 
                    pep_razor_unique = ifelse(grepl(";", `Protein group IDs`), FALSE, TRUE))
    
    if (pep_unique_by == "group") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    } else if (pep_unique_by == "protein") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_literal_unique)
    } else if (pep_unique_by == "none") {
      df <- df %>% dplyr::mutate(pep_isunique = TRUE)
    } else {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    }
    
    df <- df %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
    
    # An unhandled case with MaxQuant: "single protein but multiple groups"
    #   unique at protein but not at group
    # Sequence	Proteins	Gene Names	Protein group IDs	
    # EAEDSLRR	NP_004530	NARS1	4870;13627
  })

  # MaxQuant exception and not necessarily "Reverse" entries
  df <- df %>% dplyr::filter(as.character(.$Proteins) > 0)

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `msms[...].txt` for `filter_` varargs.")
  
  df <- df %>% filters_in_call(!!!filter_dots) 
  if ("id" %in% names(df)) df <- df %>% dplyr::rename(ID = id)

  if (corrected_int) {
    df <- df %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}\\d+$", names(.)))
  } else {
    df <- df %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}corrected\\s{1}\\d+$", names(.)))
  }
  
  if (rm_craps) {
    df <- df %>% dplyr::filter(!grepl("^CON_", .[["Proteins"]]))
  }
  
  if (rm_reverses) {
    df <- df %>% dplyr::filter(is.na(Reverse) || Reverse != "+")
  }

  df <- df %>% 
    `names_pos<-`(grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+$", names(.)), 
                  gsub("TMT-", "I", as.character(TMT_levels)))
  
  col_int <- find_int_cols(TMT_plex)
  
  df <- local({
    phos_idx <- grep("Phospho (STY) Probabilities", names(df), fixed = TRUE)
    
    if (length(phos_idx) > 1) phos_idx <- phos_idx[1]

    if (!purrr::is_empty(phos_idx)) {
      if (purge_phosphodata) {
        df <- df %>% 
          dplyr::filter(!nchar(as.character(.[["Phospho (STY) Probabilities"]])) == 0) 
      }

      phos <- df %>% 
        dplyr::select(phos_idx) %>% 
        purrr::map(~ stringr::str_extract_all(.x, "\\([^()]+\\)")) %>% 
        .[[1]] %>% 
        purrr::map(~ substring(.x, 2, nchar(.x) - 1)) %>% 
        purrr::map(as.numeric) %>% 
        purrr::map(sort, decreasing = TRUE)
      
      # may contain numeric(0) from non phospho entries
      phos_max <- suppressWarnings(phos %>% purrr::map(`[`, 1)) %>% unlist()
      sec_max <- suppressWarnings(phos %>% purrr::map(`[`, 2)) %>% unlist()
      
      df <- dplyr::bind_cols(df, 
                             pep_phospho_locprob = phos_max, 
                             pep_phospho_locdiff = phos_max - sec_max) %>% 
        dplyr::mutate(is_locprob_one = magrittr::equals(1, pep_phospho_locprob)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locdiff"), 
                         ~ replace(.x, is_locprob_one, 1)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locprob"), 
                         ~ replace(.x, is.infinite(.x), NA)) %>% 
        dplyr::select(-is_locprob_one)
    }
    
    return(df)
  })
  
  # make available `prot_acc` & `pep_seq` 
  if (TMT_plex > 0) {
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), 
                       ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
  } else {
    # (1) `Modified sequence` not yet available in MaxQuant `peptides.txt` 
    #     for matching with those in `msms.ttx`
    # (2) `Precusor intensity` not available in `msms.txt`
    
    warning("Currently with MaxQuant LFQ,\n", 
            "(1) `group_psm_by` will only be `pep_seq`.\n", 
            # "(2) `rptr_intco` has no effects.", 
            call. = FALSE)
    
    if (group_psm_by == "pep_seq_mod") {
      stop("Currently only `group_psm_by = pep_seq` to match queries in `peptides.txt`.", 
           call. = FALSE)
    }
  }
  
  df <- df %>% 
    dplyr::rename(
      pep_seq = Sequence, 
      pep_score = Score,
      pep_expect = PEP, 
      RAW_File = `Raw file`, # for consistency with metadata
    ) %>% 
    dplyr::mutate(prot_acc = gsub("\\;.*", "", Proteins)) %>% 
    reloc_col_after("pep_expect", "pep_score")
  
  df <- df %>% 
    add_entry_ids("pep_seq", "pep_index") %>% 
    add_entry_ids("prot_acc", "prot_index")

  # cannot to (re)calculate `pep_exp_mr`: fixed modifications not contained in msms.txt
  # nor TMT tags as one may only specify a fix TMT to e.g. N-term but not lysine
  # be nice to have fixed mod info and peak area to make msms.txt more self-contained

  df <- df %>% annotPrn(fasta, entrez) %>% 
    dplyr::select(-which(names(.) %in% c("Gene Names", "Gene names", 
                                         "Protein Names", "Protein Names")))
  
  acc_types <- unique(df$acc_type) %>% .[!is.na(.)]
  
  if (rm_craps) {
    craps <- load_craps(acc_types)
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) {
    df <- df %>% 
      split(., .$acc_type, drop = TRUE) %>% 
      .[acc_types] %>% 
      purrr::imap( ~ annotKin(.x, .y)) %>% 
      do.call(rbind, .)
  } else {
    df <- df %>% 
      dplyr::mutate(kin_attr = FALSE, kin_class = NA, kin_order = NA)
  }
  
  # make available `pep_res_before` and `pep_res_after`
  df <- df %>% 
    annotPeppos() %>% 
    reloc_col_before("pep_seq", "pep_res_after") %>% 
    reloc_col_before("pep_res_before", "pep_seq")
  
  df <- df %>% 
    dplyr::mutate(pep_seq_bare = gsub("^.*\\.([^\\.]+)\\..*", "\\1", pep_seq)) %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq_bare)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq_bare), 
                                    stringr::str_count(pep_seq_bare, "[KR]") - 1,
                                    stringr::str_count(pep_seq_bare, "[KR]"))) %>% 
    add_prot_icover() %>% 
    dplyr::select(-pep_seq_bare) %>% 
    dplyr::select(-which(names(.) %in% c("Length", "Missed cleavages", "Missed Clevages"))) 
  
  if (!("prot_cover" %in% names(df) && length(filelist) == 1)) {
    df <- df %>% calc_cover(id = !!rlang::sym(group_pep_by))
  }
  
  # --  column arrangements --- 
  df <- df %>% reloc_col_after("prot_index", "prot_cover") 
  
  df <- df %>% 
    order_psm_cols(c("gene", "fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                     "other_acc", "entrez", "species", "acc_type", 
                     "kin_attr", "kin_class", "kin_order", "is_tryptic"))

  df <- cbind.data.frame(
    df[, grepl("^[a-z]", names(df))], 
    df[, grepl("^[A-Z]", names(df)) & !grepl("^[IR][0-9]{3}[NC]{0,1}", names(df))], 
    df[, grepl("^[R][0-9]{3}[NC]{0,1}", names(df))], 
    df[, grepl("^[I][0-9]{3}[NC]{0,1}", names(df))]
  ) %>% reloc_col_before("m/z", "Mass")
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))),
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
  )

  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), 
                       as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                       ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                       ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) 
  } else {
    df <- df %>%
      dplyr::mutate_at(.vars = grep("^Precursor Intensity$|^Intensity$", names(.)), 
                       as.numeric) %>% 
      dplyr::mutate_at(.vars = grep("^Precursor Intensity$|^Intensity$", names(.)), 
                       ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      # temporarily disabled as values not filled for Bruker's .d files
      # dplyr::filter(!is.na(.[grep("^Precursor Intensity$|^Intensity$", names(.))])) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  # make available `pep_seq_mod` 
  stopifnot("Modified sequence" %in% names(df))
  
  df <- df %>% 
    dplyr::rename(pep_seq_mod = `Modified sequence`) %>%
    add_maxquant_pepseqmod(use_lowercase_aa) %>% 
    reloc_col_before("dat_file", "RAW_File") 
  
  df <- df %>% 
    reloc_col_after("pep_seq_mod", "pep_seq") %>% 
    reloc_col_before("pep_res_before", "pep_seq") %>% 
    reloc_col_after("pep_res_after", "pep_seq_mod") %>% 
    reloc_col_after("pep_start", "pep_res_after") %>% 
    reloc_col_after("pep_end", "pep_start") %>% 
    reloc_col_after("pep_len", "pep_end") %>% 
    reloc_col_after("pep_miss", "pep_len") %>% 
    reloc_col_before("pep_isunique", "pep_res_before") %>% 
    reloc_col_before("pep_literal_unique", "pep_res_before") %>% 
    reloc_col_before("pep_razor_unique", "pep_res_before") %>% 
    reloc_col_before("pep_score", "pep_res_before") %>% 
    reloc_col_before("pep_expect", "pep_res_before")
  
  # --- compare experimental and user-provided PSM files ---
  res <- check_raws(df)
  df <- res$df
  tmtinj_raw_map <- res$lookup
  rm(res)
  
  # --- split by TMT and LCMS ---
  if (! "PSM_File" %in% names(tmtinj_raw_map)) {
    df_split <- df %>%
      dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>%
      split(., .$TMT_inj, drop = TRUE)
  } else {
    tmtinj_raw_map <- tmtinj_raw_map %>% 
      dplyr::mutate(PSM_File = gsub("\\.txt$", "", PSM_File)) %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "PSM_File"), sep = "@") 
    
    df_split <- df %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "dat_file"), 
                   sep = "@", remove = FALSE) %>% 
      dplyr::left_join(tmtinj_raw_map, by = "RAW_File2") %>% 
      dplyr::select(-RAW_File2) %>% 
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>%
      split(., .$TMT_inj, drop = TRUE)
  }
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    warning("TMT sets and LC/MS injections not have corresponindg PSM files:\n", 
            call. = FALSE)
    message(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` from experiment summary file."), 
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")

  df_split <- df_split %>% .[names(.) %>% sort_tmt_lcms()]
  
  if (parallel) {
    n_cores <- pmin(parallel::detectCores(), length(filelist))
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

    suppressWarnings(
      silent_out <- parallel::clusterMap(
        cl, psm_msplit, df_split, names(df_split), 
        MoreArgs = list(fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
      )
    )
    
    parallel::stopCluster(cl)
  } else {
    purrr::walk2(df_split, names(df_split), psm_msplit, 
                 fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
  }
}


#' Annotates MaxQuant PSM results
#'
#' \code{annotPSM_mq} adds fields of annotation to MaxQuant PSM tables.
#'
#' @inheritParams load_expts
#' @inheritParams annotPSM
#' @inheritParams normPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @param ... Not currently used.
annotPSM_mq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", mc_psm_by = "peptide", 
                        fasta = NULL, expt_smry = "expt_smry.xlsx", 
                        plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                        ...) {
  
  dat_dir <- get_gl_dat_dir()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM/cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files()
  
  set_indexes <- gsub("^TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort() 
  
  for (set_idx in set_indexes) {
    sublist <- filelist[grep(paste0("^TMTset", set_idx, "_"), 
                             filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(
      Filename = do.call('rbind', strsplit(as.character(sublist), '.txt', fixed = TRUE))
    ) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    # --- LCMS injections under the same TMT experiment --- 
    for (idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM/cache", sublist[idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") 

      df <- df %>% 
        dplyr::mutate(pep_start_discrepancy = 
                        ifelse(grepl("Protein N-term", Modifications) & (pep_start > 2), 
                               TRUE, FALSE)) %>% 
        dplyr::filter(!pep_start_discrepancy) %>% 
        dplyr::select(-pep_start_discrepancy)
      
      if (TMT_plex > 0) {
        df <- local({
          injn_idx <- sublist[idx] %>% 
            gsub("^TMTset\\d+_LCMSinj(\\d+)_Clean.txt$", "\\1", .) %>% 
            as.integer()
          
          df <- mcPSM(df, set_idx, injn_idx, mc_psm_by, group_psm_by, group_pep_by)
        })
      }

			df <- df %>% 
			  calcSD_Splex(group_psm_by) %>% 
			  `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
			  dplyr::right_join(df, by = group_psm_by) %>% 
			  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
			                   ~ round(.x, digits = 0)) %>% 
			  dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), 
			                   ~ round(.x, digits = 3)) %>% 
			  dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), 
			                   ~ round(.x, digits = 3)) %>% 
			  dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
			                   ~ round(.x, digits = 4)) %>% 
			  na_genes_by_acc() 

			pep_n_psm <- df %>%
			  dplyr::select(!!rlang::sym(group_psm_by)) %>%
			  dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
			  dplyr::summarise(pep_n_psm = n())
			
			prot_n_psm <- df %>% 
			  dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
			  dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
			  dplyr::summarise(prot_n_psm = n())
			
			prot_n_pep <- df %>%
			  dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
			  dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
			  dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
			  dplyr::summarise(prot_n_pep = n())
			
			df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)

			df <- list(df, prot_n_psm, prot_n_pep) %>%
			  purrr::reduce(left_join, by = group_pep_by)
			
			df <- dplyr::bind_cols(
			  df %>% dplyr::select(grep("^prot_", names(.))),
			  df %>% dplyr::select(grep("^pep_", names(.))), 
			  df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
			) %>% 
			  reloc_col_before("pep_seq", "pep_res_after") %>% 
			  reloc_col_after("pep_seq_mod", "pep_seq") 

			df <- dplyr::bind_cols(
			  df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]{0,1}", names(.))), 
			  df %>% dplyr::select(grep("^I[0-9]{3}[NC]{0,1}", names(.))), 
			  df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]{0,1}", names(.))), 
			  df %>% dplyr::select(grep("^R[0-9]{3}[NC]{0,1}", names(.))), 
			  df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]{0,1}", names(.))), 
			  df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))), 
			  df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))),
			) %>% 
			  dplyr::select(-which(names(.) %in% c("pep_index", "prot_index"))) %T>% 
			  write.table(file.path(dat_dir, "PSM", paste0(out_fn[idx, 1], ".txt")),
			                   sep = "\t", col.names = TRUE, row.names = FALSE)
      
			if (plot_rptr_int) {
			  if (TMT_plex > 0) {
			    df_int <- df %>% .[, grepl("^N_I[0-9]{3}", names(.))]
			  } else {
			    df_int <- df %>% 
			      dplyr::select("Precursor Intensity") %>% 
			      dplyr::rename(N_I000 = "Precursor Intensity")
			  }
			  
			  rptr_violin(
			    df = df_int, 
			    filepath = file.path(dat_dir, "PSM/rprt_int/mc", 
			                         paste0(gsub("_PSM_N", "", out_fn[idx, 1]), "_rprt.png")), 
			    width = 8, height = 8)
			}
      
      if (plot_log2FC_cv && TMT_plex > 0) {
        sd_violin(
          df = df, id = !!group_psm_by, 
          filepath = file.path(dat_dir, "PSM/log2FC_cv/raw", 
                               paste0(gsub("_PSM_N", "", out_fn[idx, 1]), "_sd.png")), 
          width = 8, height = 8, type = "log2_R", 
          adjSD = FALSE, is_psm = TRUE, , ...)
      }
    }
  }
}


#' Pad Spectrum Mill TMT channels to the highest plex
#' @param file A file name of PSM table from MaxQuant export
pad_sm_channels <- function(file) {
  dat_dir <- get_gl_dat_dir()
  
  base_name <- file %>% gsub("\\.ssv$", "", .)
  
  df <- suppressWarnings(readr::read_delim(file.path(dat_dir, file), delim = ";", 
                                           col_types = cols(filename = col_character())))
  
  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  
  if (! "PSM_File" %in% names(fraction_scheme)) {
    label_scheme_sub <- local({
      raw_files <- df$filename %>% 
        gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", .) %>% 
        unique() %>% 
        gsub("\\.raw$", "", .) %>% 
        gsub("\\.d$", "", .)
      
      tmt_sets <- fraction_scheme %>% 
        dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File), 
                      RAW_File = gsub("\\.d$", "", RAW_File)) %>% 
        dplyr::filter(RAW_File %in% raw_files, !duplicated(TMT_Set)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(TMT_Set) %>% 
        unlist()
      
      label_scheme_full %>% dplyr::filter(TMT_Set %in% tmt_sets)
    })
  } else {
    label_scheme_sub <- local({
      tmt_sets <- fraction_scheme %>% 
        dplyr::mutate(PSM_File = gsub("\\.ssv$", "", PSM_File)) %>% 
        dplyr::filter(PSM_File == base_name, !duplicated(TMT_Set)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(TMT_Set) %>% 
        unlist()
      
      label_scheme_full %>% dplyr::filter(TMT_Set %in% tmt_sets)
    })
  }

  nas <- data.frame(rep(NA, nrow(df)))
  sample_ids <- as.character(label_scheme_sub$Sample_ID)
  
  str_ratio <- "^TMT_[0-9]{3}[NC]{0,1}_[0-9]{3}[NC]{0,1}"
  str_int <- "^TMT_[0-9]{3}[NC]{0,1}$"
  
  df_int <- df %>% dplyr::select(grep(str_int, names(.)))
  
  ref <- names(df) %>% 
    .[grepl(str_ratio, .)] %>% 
    gsub(".*_([0-9]{3}[NC]{0,1})$", "\\1", .) %>% 
    unique()
  
  stopifnot(length(ref) == 1)
  
  this_plex <- ncol(df_int)
  TMT_plex <- TMT_plex(label_scheme_full)
  stopifnot(this_plex <= TMT_plex, this_plex >= 0)
  
  if (this_plex > 0 && this_plex < TMT_plex) {
    if (TMT_plex == 16) {
      keys_int <- paste0("TMT_", c("126", "127N", "127C", "128N", "128C", "129N", "129C", 
                                   "130N", "130C", "131N", "131C", "132N", "132C", 
                                   "133N", "133C", "134N"))
    } else if (TMT_plex == 11) {
      keys_int <- paste0("TMT_", c("126", "127N", "127C", "128N", "128C", "129N", "129C", 
                                   "130N", "130C", "131N", "131C"))
    } else if (TMT_plex == 10) {
      keys_int <- paste0("TMT_", c("126", "127N", "127C", "128N", "128C", "129N", "129C", 
                                   "130N", "130C", "131"))
    } else if(TMT_plex == 6) {
      keys_int <- paste0("TMT_", c("126", "127", "128", "129", "130", "131"))
    } else {
      keys_int <- NULL
    }
    keys_ratio <- paste0(keys_int, "_", ref) %>% .[!grepl(paste0("_", ref, "_", ref), .)]
    
    if (this_plex == 6) {
      if (TMT_plex == 10) {
        pos <- c(3, 5, 7, 9)
      } else if (TMT_plex == 11) {
        pos <- c(3, 5, 7, 9, 11)
      } else if (TMT_plex == 16) {
        pos <- c(3, 5, 7, 9, 11:16)
      } else {
        stop("TMT_plex not one of c(10, 11, 16).", call. = FALSE)
      }
    } else if (this_plex == 10) {
      if (TMT_plex == 11) {
        pos <- 11
      } else if (TMT_plex == 16) {
        pos <- c(11:16)
      } else {
        stop("TMT_plex not one of c(11, 16).", call. = FALSE)
      }
    } else if (this_plex == 11) {
      if (TMT_plex == 16) {
        pos <- c(12:16)
      } else {
        stop("TMT_plex not one of c(16).", call. = FALSE)
      }
    }
    
    for (idx in seq_along(pos)) {
      df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))
    }
    rm(idx)
    
    if (ncol(df_int) == length(sample_ids)) {
      names(df_int) <- keys_int
      df <- replace_cols_at(df, df_int, grep(str_int, names(df)))
    }
    
    df_ratio <- local({
      df_int <- as.data.frame(df_int)
      
      sweep(df_int, 1, df_int[, paste0("TMT_", ref)], "/") %>% 
        `colnames<-`(paste0(names(.), "_", ref)) %>% 
        dplyr::select(-grep(paste0("^TMT_", ref, "_", ref, "$"), names(.))) %>% 
        dplyr::mutate_all(~ replace(.x, is.infinite(.x), NA))    
    })
    
    if ((ncol(df_ratio) + 1) == TMT_plex) {
      names(df_ratio) <- keys_ratio
      df <- replace_cols_at(df, df_ratio, grep(str_ratio, names(df)))
    }
  }
  
  df$dat_file <- base_name
  
  return(df)
}


#' Add the `pep_seq_mod` field to MaxQuant PSMs
#'
#' @inheritParams splitPSM
#' @inheritParams locate_outliers
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_sm_pepseqmod <- function(df, use_lowercase_aa) {
  
  if (!use_lowercase_aa) {
    df <- df %>% 
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste0(pep_seq, "[", variableSites, "]"))
  } else {
    # (1) mods under column `variableSites`
    df <- local({
      if (!is_empty(which(names(df) == "variableSites"))) {
        df_sub <- df %>% dplyr::filter(!is.na(variableSites))
        
        if (nrow(df_sub) > 0) {
          pos_matrix <- df_sub %>% 
            dplyr::select(variableSites) %>% 
            dplyr::mutate(variableSites = as.character(variableSites)) %$% 
            stringr::str_split(.$variableSites, " ") %>% 
            plyr::ldply(., rbind) %>% 
            `names<-`(paste0("mod_", names(.))) %>% 
            purrr::map(~ gsub("[A-z]", "", .)) %>% 
            dplyr::bind_cols() %>% 
            dplyr::mutate_at(.vars = grep("^mod_", names(.)), as.numeric) %>% 
            dplyr::bind_cols(df_sub %>% dplyr::select(pep_seq, pep_start)) %>% 
            dplyr::mutate_at(.vars = grep("^mod_", names(.)), ~ {.x + 1 - pep_start}) %>% 
            dplyr::select(-pep_seq, -pep_start) %>% 
            data.frame(check.names = FALSE)
          
          for (k in seq_along(pos_matrix)) {
            rows <- !is.na(pos_matrix[, k])
            locales <- pos_matrix[rows, k]
            
            lowers <- substr(df_sub$pep_seq_mod[rows], locales, locales) %>% tolower()
            substr(df_sub$pep_seq_mod[rows], locales, locales) <- lowers
          }
          
          df <- dplyr::bind_rows(
            df %>% dplyr::filter(is.na(variableSites)), 
            df_sub
          )            
        }
      }
      
      return(df)
    })
    
    # (2-1) protein n-term acetylation
    if (!purrr::is_empty(which(names(df) == "nterm"))) {
      df_sub <- df %>% dplyr::filter(nterm == "Acetyl")
      
      if (nrow(df_sub) > 0) {
        df_sub$pep_seq_mod <- paste0("_", df_sub$pep_seq_mod)
        
        df <- dplyr::bind_rows(
          df %>% dplyr::filter(nterm != "Acetyl"), 
          df_sub
        )
      }        
    }
    
    # (2-2) protein C-terminal amidation
    if (!purrr::is_empty(which(names(df) == "cterm"))) {
      # hypothetical; no yet known how it will be named in SM
      df_sub <- df %>% dplyr::filter(grepl("^Amidate", cterm))
      
      if (nrow(df_sub) > 0) {
        df_sub$pep_seq_mod <- paste0(df_sub$pep_seq_mod, "_")
        
        df <- dplyr::bind_rows(
          df %>% dplyr::filter(!grepl("^Amidate", cterm)), 
          df_sub
        )
      }
    }
    
    # (3-1) other protien n-term not yet defined in SM
    # ...
    
    # (3-2) other protien c-term not yet defined in SM
    # ...
    
    # (4-1) peptide "(N-term)" modification: 
    #   only "Pyroglutamic acid (N-termQ)" and under column `variableSites`;
    #   as a result, becomes lower-case `q` after step (1)
    
    # (4-2) peptide "(C-term)" modification not yet defined in SM
    # ...
    
    # paste "pep_res_before" and "pep_res_after"
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, 
                                    sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, 
                                        sep = "."))        
  }
  
  # `pep_seq` may not be in fasta
  df %>% 
    dplyr::mutate(
      pep_seq = gsub("^NA.", "?.", pep_seq), 
      pep_seq = gsub(".NA$", ".?", pep_seq), 
      pep_seq_mod = gsub("^NA.", "?.", pep_seq_mod), 
      pep_seq_mod = gsub(".NA$", ".?", pep_seq_mod)
    )
}


#' Splits PSM tables from \code{Spectrum Mill}
#'
#' \code{splitPSM_sm} splits the PSM outputs from \code{Spectrum Mill}. It
#' separates PSM data by TMT experiment and LC/MS injection.
#'
#' @inheritParams splitPSM
#' @inheritParams splitPSM_mq
#' @inheritParams normPSM
#' @import dplyr tidyr readr
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
splitPSM_sm <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        fasta = NULL, entrez = NULL, pep_unique_by = "group", 
                        rm_craps = FALSE, rm_krts = FALSE, purge_phosphodata = TRUE, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, 
                        rptr_intco = 0, 
                        use_lowercase_aa = TRUE, parallel = TRUE, ...) {
  
  on.exit(message("Split PSMs by sample IDs and LCMS series --- Completed."), 
          add = TRUE)
  
  if (is.null(fasta)) stop("FASTA file(s) not provided.", call. = FALSE)
  
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$")
  
  if (rlang::is_empty(filelist)) 
    stop(paste("No PSM files were found under", file.path(dat_dir), 
               "\nCheck that the names of PSM files start with `PSMexport`."), 
         call. = FALSE)
  
  df <- purrr::map(filelist, pad_sm_channels)
  df <- suppressWarnings(dplyr::bind_rows(df))
  nms_sm <- names(df) %>% .[! . %in% c("dat_file")]
  
  df <- df %>% 
    add_entry_ids("sequence", "pep_index") %>% 
    add_entry_ids("accession_number", "prot_index")
  
  # in addition to true unique
  # if shared sequences with accession_numbers not in primary accessions -> razor?

  df <- local({
    stopifnot(c("sequence", "accession_numbers") %in% names(df))
    
    peps_shared_by_prots <- df %>% 
      find_literal_unique_peps("sequence", "accession_numbers")
    
    df <- df %>% 
      dplyr::mutate(pep_isunique = NA, # not yet available
                    pep_literal_unique = ifelse(sequence %in% peps_shared_by_prots, 
                                                FALSE, TRUE), 
                    pep_razor_unique = pep_isunique) %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
    
    if (pep_unique_by == "group") {
      df <- df
    } else if (pep_unique_by == "protein") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_literal_unique)
    } else if (pep_unique_by == "none") {
      df <- df %>% dplyr::mutate(pep_isunique = TRUE)
    } else {
      df <- df
    }
  })

  is_phospho_expt <- any(grepl("Phosphorylated", df$modifications))
  if (is_phospho_expt && purge_phosphodata) {
    if ("modifications" %in% names(df)) {
      df <- df %>% dplyr::filter(grepl("Phosphorylated", modifications))
    } else {
      warning("Missing column `modifications` for non-phosphopeptide removals.", 
              call. = FALSE)
    }
  }
  rm(is_phospho_expt)
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `PSMexport[...].ssv` for `filter_` varargs.")
  
  # make available `prot_acc`
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate(pep_expect = NA) %>% 
    dplyr::mutate(pep_isunique = NA) %>% 
    dplyr::rename(prot_acc = accession_number, 
                  pep_score = score) %>% 
    annotPrn(fasta, entrez)  
  
  col_int <- find_int_cols(TMT_plex)
  
  if (TMT_plex > 0) {
    df <- df %>% 
      `names<-`(gsub("^TMT_([0-9]{3}[NC]?)", "I\\1", names(.))) %>% 
      dplyr::select(-grep("^I[0-9]{3}[NC]?_[0-9]{3}[NC]?$", names(.))) %>% 
      as.data.frame()
    
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), 
                       ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
    
    # make available `pep_seq`
    df <- df %>% 
      dplyr::rename(RAW_File = `filename`) %>% 
      dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", RAW_File)) %>% 
      dplyr::mutate(pep_seq = toupper(sequence)) 
  }
  
  acc_types <- unique(df$acc_type) %>% .[!is.na(.)]
  
  if (rm_craps) {
    craps <- load_craps(acc_types)
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) {
    df <- df %>% 
      split(., .$acc_type, drop = TRUE) %>% 
      .[acc_types] %>% 
      purrr::imap( ~ annotKin(.x, .y)) %>% 
      do.call(rbind, .)
  } else {
    df <- df %>% 
      dplyr::mutate(kin_attr = FALSE, kin_class = NA, kin_order = NA)
  }
  
  # make available `pep_res_before` and `pep_res_after`
  df <- df %>% 
    annotPeppos() %>% 
    dplyr::mutate(pep_seq_bare = gsub("^.*\\.([^\\.]+)\\..*", "\\1", pep_seq)) %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq_bare)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq_bare), 
                                    stringr::str_count(pep_seq_bare, "[KR]") - 1,
                                    stringr::str_count(pep_seq_bare, "[KR]"))) %>% 
    add_prot_icover() %>% 
    dplyr::select(-pep_seq_bare)

  # M._sequence.c; -._sequence.c; n.sequence.c; -.sequence.c
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df <- df %>% calc_cover(id = !!rlang::sym(group_pep_by))
  }
  
  df <- df %>% dplyr::mutate(pep_seq = toupper(pep_seq))

  # --  column arrangements --- 
  df <- df %>% reloc_col_after("prot_index", "prot_cover") 
  
  df <- df %>% 
    order_psm_cols(c("gene", "fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                     "other_acc", "entrez", "species", "acc_type", 
                     "kin_attr", "kin_class", "kin_order", "is_tryptic"))
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))),
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_|^pep_", names(.))), 
    df %>% 
      dplyr::select(-grep("^prot_|^pep_|[RI]{1}[0-9]{3}[NC]{0,1}", names(.))) %>% 
      dplyr::select(-which(names(.) %in% nms_sm)), 
    df %>% 
      dplyr::select(-grep("^prot_|^pep_|[RI]{1}[0-9]{3}[NC]{0,1}", names(.))) %>% 
      dplyr::select(which(names(.) %in% nms_sm)), 
    df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
  )

  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), 
                       as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                       ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                       ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) # %>% 
      # dplyr::select(which(not_all_zero(.))) %>% # don't: empty channels are all NA too
      # dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    # don't yet know how LFQ looks like with SM
    df <- df %>% dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  # make available `pep_seq_mod` 
  df <- df %>%
    dplyr::mutate(pep_seq_mod = pep_seq) %>% 
    dplyr::mutate(pep_seq_mod = as.character(pep_seq_mod)) %>% 
    dplyr::select(which(names(.) == "pep_seq_mod"),
                  which(names(.) != "pep_seq_mod")) %>% 
    add_sm_pepseqmod(use_lowercase_aa) %>% 
    reloc_col_before("dat_file", "RAW_File") %>% 
    reloc_col_after("pep_seq_mod", "pep_seq") %>% 
    reloc_col_before("pep_res_before", "pep_seq") %>% 
    reloc_col_after("pep_res_after", "pep_seq_mod") %>% 
    reloc_col_after("pep_start", "pep_res_after") %>% 
    reloc_col_after("pep_end", "pep_start") %>% 
    reloc_col_after("pep_len", "pep_end") %>% 
    reloc_col_after("pep_miss", "pep_len") %>% 
    reloc_col_before("pep_isunique", "pep_res_before") %>% 
    reloc_col_before("pep_literal_unique", "pep_res_before") %>% 
    reloc_col_before("pep_razor_unique", "pep_res_before") %>% 
    reloc_col_before("pep_score", "pep_res_before") %>% 
    reloc_col_before("pep_expect", "pep_res_before") %>% 
    reloc_col_after("dat_file", "is_tryptic") %>% 
    reloc_col_after("RAW_File", "dat_file")

  # --- compare experimental and user-provided PSM files ---
  res <- check_raws(df)
  df <- res$df
  tmtinj_raw_map <- res$lookup
  rm(res)
  
  # --- split by TMT and LCMS ---
  if (! "PSM_File" %in% names(tmtinj_raw_map)) {
    df_split <- df %>%
      dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>%
      split(., .$TMT_inj, drop = TRUE)
  } else {
    tmtinj_raw_map <- tmtinj_raw_map %>% 
      dplyr::mutate(PSM_File = gsub("\\.ssv$", "", PSM_File)) %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "PSM_File"), sep = "@") 
    
    df_split <- df %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "dat_file"), sep = "@", remove = FALSE) %>% 
      dplyr::left_join(tmtinj_raw_map, by = "RAW_File2") %>% 
      dplyr::select(-RAW_File2) %>% 
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>% 
      split(., .$TMT_inj, drop = TRUE)
  }
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    warning("TMT sets and LC/MS injections not have corresponindg PSM files:\n", 
            call. = FALSE)
    message(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` from experiment summary file."), 
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), 
                                   sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")
  
  df_split <- df_split %>% .[names(.) %>% sort_tmt_lcms()]
  
  if (parallel) {
    n_cores <- pmin(parallel::detectCores(), length(filelist))
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))

    suppressWarnings(
      silent_out <- parallel::clusterMap(
        cl, psm_msplit, df_split, names(df_split), 
        MoreArgs = list(fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
      )
    )
    
    parallel::stopCluster(cl)
  } else {
    purrr::walk2(df_split, names(df_split), psm_msplit, 
                 fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
  }
}


#' Annotates \code{Spectrum Mill} PSMs following \code{splitPSM_sm} and
#' \code{cleanupPSM}
#'
#' \code{annotPSM_sm} adds fields of annotation to PSM tables.
#'
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams annotPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @param ... Not currently used.
annotPSM_sm <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        mc_psm_by = "peptide", 
                        fasta = NULL, expt_smry = "expt_smry.xlsx", 
                        plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                        ...) {
  
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM/cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files()
  
  set_indexes <- gsub("^TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort() 
  
  for(set_idx in set_indexes) {
    sublist <- filelist[grep(paste0("^TMTset", set_idx, "_"), filelist, 
                             ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename = 
                           do.call('rbind', strsplit(as.character(sublist),
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    for (idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM/cache", sublist[idx]),
                     check.names = FALSE, header = TRUE, 
                     sep = "\t", comment.char = "#") 
      
      # median centering
      if (TMT_plex > 0) {
        df <- local({
          injn_idx <- sublist[idx] %>% 
            gsub("^TMTset\\d+_LCMSinj(\\d+)_Clean.txt$", "\\1", .) %>% 
            as.integer()
          
          df <- mcPSM(df, set_idx, injn_idx, mc_psm_by, group_psm_by, group_pep_by)
        })
      }
      
      # add SD columns
      df <- df %>% 
        calcSD_Splex(group_psm_by) %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
        dplyr::right_join(df, by = group_psm_by) %>% 
        dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 0)) %>% 
        dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 4)) %>% 
        na_genes_by_acc() %>% 
        reloc_col_before("pep_seq", "pep_res_after") %>% 
        reloc_col_after("pep_seq_mod", "pep_seq")
      
      # add column `pep_n_psm` etc.
      pep_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by)) %>%
        dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
        dplyr::summarise(pep_n_psm = n())
      
      prot_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_psm = n())
      
      prot_n_pep <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_pep = n())
      
      df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)

      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))),
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
      ) 
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^I[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))),
      ) %>% 
        dplyr::select(-which(names(.) %in% c("pep_index", "prot_index"))) %T>% 
        write.table(file.path(dat_dir, "PSM", paste0(out_fn[idx, 1], ".txt")),
                         sep = "\t", col.names = TRUE, row.names = FALSE)

      if (plot_rptr_int & TMT_plex > 0) {
        df_int <- df %>% 
          .[, grepl("^N_I[0-9]{3}", names(.))]
        
        rptr_violin(
          df = df_int, 
          filepath = file.path(dat_dir, "PSM/rprt_int/mc", 
                               paste0(gsub("_PSM_N", "", out_fn[idx, 1]), "_rprt.png")), 
          width = 8, height = 8)
      }
      
      if (plot_log2FC_cv & TMT_plex > 0) {
        sd_violin(
          df = df, 
          id = !!group_psm_by, 
          filepath = file.path(dat_dir, "PSM/log2FC_cv/raw", 
                               paste0(gsub("_PSM_N", "", out_fn[idx, 1]), "_sd.png")), 
          width = 8, height = 8, type = "log2_R", 
          adjSD = FALSE, is_psm = TRUE, , ...)
      }
    }
  }
}


#' Pad MSFragger TMT channels to the highest plex. 
#' 
#' @param file A file name of PSM table from MaxQuant export
#' @inheritParams splitPSM
pad_mf_channels <- function(file, fasta, entrez) {
  dat_dir <- get_gl_dat_dir()
  
  base_name <- file %>% gsub("\\.tsv$", "", .)
  
  df <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, file), 
                    col_types = cols(`Is Unique` = col_logical()))
  )

  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  
  stopifnot("Spectrum" %in% names(df))
  
  df <- df %>% 
    dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", Spectrum))
  
  if (! "PSM_File" %in% names(fraction_scheme)) {
    label_scheme_sub <- local({
      raw_files <- unique(df$RAW_File) %>% 
        gsub("\\.raw$", "", .) %>% 
        gsub("\\.d$", "", .)
      
      tmt_sets <- fraction_scheme %>% 
        dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File), 
                      RAW_File = gsub("\\.d$", "", RAW_File)) %>% 
        dplyr::filter(RAW_File %in% raw_files, !duplicated(TMT_Set)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(TMT_Set) %>% 
        unlist()
      
      label_scheme_full %>% dplyr::filter(TMT_Set %in% tmt_sets)
    })
  } else {
    label_scheme_sub <- local({
      tmt_sets <- fraction_scheme %>% 
        dplyr::mutate(PSM_File = gsub("\\.tsv$", "", PSM_File)) %>% 
        dplyr::filter(PSM_File == base_name, !duplicated(TMT_Set)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(TMT_Set) %>% 
        unlist()
      
      label_scheme_full %>% dplyr::filter(TMT_Set %in% tmt_sets)
    })
  }
  
  nas <- data.frame(rep(NA, nrow(df)))
  sample_ids <- as.character(label_scheme_sub$Sample_ID)
  
  str_int1 <- "^Reporter intensity [0-9]+"
  str_int2 <- "^Reporter intensity corrected [0-9]+"
  str_dev <- "^Reporter mass deviation \\[mDa\\] [0-9]+"
  
  df_int <- df %>% dplyr::select(grep(str_int1, names(.)))
  df_int2 <- df %>% dplyr::select(grep(str_int2, names(.)))
  df_dev <- df %>% dplyr::select(grep(str_dev, names(.)))
  
  this_plex <- ncol(df_int)
  TMT_plex <- TMT_plex(label_scheme_full)
  stopifnot(this_plex <= TMT_plex, this_plex >= 0)
  
  # Empty.xxx can be due to either channel padding or removals
  if (this_plex > 0 && this_plex < TMT_plex) {
    if (this_plex == 6) {
      if (TMT_plex == 10) {
        pos <- c(3, 5, 7, 9)
      } else if (TMT_plex == 11) {
        pos <- c(3, 5, 7, 9, 11)
      } else if (TMT_plex == 16) {
        pos <- c(3, 5, 7, 9, 11:16)
      } else {
        stop("TMT_plex not one of c(10, 11, 16).", call. = FALSE)
      }
    } else if (this_plex == 10) {
      if (TMT_plex == 11) {
        pos <- 11
      } else if (TMT_plex == 16) {
        pos <- c(11:16)
      } else {
        stop("TMT_plex not one of c(11, 16).", call. = FALSE)
      }
    } else if (this_plex == 11) {
      if (TMT_plex == 16) {
        pos <- c(12:16)
      } else {
        stop("TMT_plex not one of c(16).", call. = FALSE)
      }
    }
    
    for (idx in seq_along(pos)) {
      df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))
      df_int2 <- suppressMessages(add_cols_at(df_int2, nas, pos[idx] - 1))
      df_dev <- suppressMessages(add_cols_at(df_dev, nas, pos[idx] - 1))
    }
    rm(idx)
    
    len <- length(sample_ids)
    if (ncol(df_int) == len && ncol(df_int2) == len && ncol(df_dev) == len) {
      names(df_int) <- paste("Reporter intensity", seq_len(len))
      df <- replace_cols_at(df, df_int, grep(str_int1, names(df)))
      
      names(df_int2) <- paste("Reporter intensity corrected", seq_len(len))
      df <- replace_cols_at(df, df_int2, grep(str_int2, names(df)))
      
      names(df_dev) <- paste("Reporter mass deviation [mDa]", seq_len(len))
      df <- replace_cols_at(df, df_dev, grep(str_dev, names(df)))
    }
  }
  
  # the same `raw_file` may be at different `dat_file`s
  df$dat_file <- base_name
  
  invisible(df)
}



#'Splits PSM tables
#'
#'\code{splitPSM_mf} splits the PSM outputs by TMT experiment and LC/MS
#'injection.
#'@inheritParams splitPSM
#'@import dplyr tidyr
#'@importFrom stringr str_split
#'@importFrom magrittr %>% %T>% %$% %<>% equals
splitPSM_mf <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        fasta = NULL, entrez = NULL, 
                        pep_unique_by = "group", 
                        rm_craps = FALSE, rm_krts = FALSE,
                        purge_phosphodata = TRUE, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, 
                        rptr_intco = 1000, use_lowercase_aa = TRUE, 
                        parallel = TRUE, ...) {
  
  on.exit(message("Split PSMs by sample IDs and LCMS series --- Completed."), 
          add = TRUE)
  
  # the same pep_seq_mod may have different
  # (1) `Charge`, (2) `Retention` and (3) `Delta Mass`
  
  # use top-3 intensities for the same pep_seq_mod with 
  # (1) isotope effect of `Delta Mass` ignored
  # (2) `Retention` difference ignored for now
  # (3) `Charge` difference ignored; may later be distinguished
  
  # numeric column `top_n_used` after top-3; the higher the more ambiguous
  # logical column `multi_charge_states`
  # argument for aligning charge states `align_lfq_cs` 
  #   in Pep2PSM -> pep_seq_mod@2, pep_seq_mod@3
  # ... -> Pep2Prn top-3 pep_seq_mod@2 pep_seq_mod@3 etc. 
  # may be column `pep_seq_mod_more` K.peptidek.-@2 for charge 
  
  # `Total Intensity`,	`Unique Intensity`,	`Razor Intensity`
  # `Unique Intensity`: NA under `Mapped Genes`,	`Mapped Proteins`

  # filename: 3 or 2 additional fields after Raw file???
  # assume it is always 3 for now.
  
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^psm.*\\.tsv$")
  
  if (purrr::is_empty(filelist)) {
    stop(paste("No PSM files of `psm[...].tsv` under", file.path(dat_dir), ".",
               "\nMake sure that the names of PSM files start with `psm`."), 
         call. = FALSE)
  }
  
  # not yet for TMT
  df <- purrr::map(filelist, pad_mf_channels)
  df <- suppressWarnings(dplyr::bind_rows(df))
  
  df <- df %>% 
    add_entry_ids("Peptide", "pep_index") %>% 
    add_entry_ids("Protein", "prot_index")
  
  df <- local({
    stopifnot(c("Is Unique", "Mapped Proteins") %in% names(df))

    df <- df %>% 
      dplyr::mutate(pep_literal_unique = ifelse(`Is Unique` & is.na(`Mapped Proteins`), 
                                                TRUE, FALSE), 
                    pep_razor_unique = `Is Unique`)

    if (pep_unique_by == "group") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    } else if (pep_unique_by == "protein") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_literal_unique)
    } else if (pep_unique_by == "none") {
      df <- df %>% dplyr::mutate(pep_isunique = TRUE)
    } else {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    }
    
    df <- df %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
  })

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `psm[...].tsv` for `filter_` varargs.")
  
  df <- df %>% filters_in_call(!!!filter_dots) 

  # proteins:
  # (1) prot_tot_int: Gene == Ckap5 | Mapped Genes == Ckap5
  # (2) prot_razor_int: Gene == Ckap5
  # (3) prot_uniq_int: Gene == Ckap5 & void Mapped Genes
  # top3 for each
  
  # only with LFQ
  df <- df %>% 
    dplyr::mutate(pep_tot_int = Intensity) %>% 
    dplyr::mutate(pep_unique_int = ifelse(
      `Is Unique` & is.na(`Mapped Proteins`), Intensity, 0)) %>% 
    dplyr::mutate(pep_razor_int = ifelse(`Is Unique`, Intensity, 0))

  # not yet for TMT
  run_scripts <- FALSE
  if (run_scripts) {
    df <- df %>% 
      `names_pos<-`(grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+$", names(.)), 
                    gsub("TMT-", "I", as.character(TMT_levels)))
  }

  col_int <- find_int_cols(TMT_plex)
  
  df <- local({
    phos_idx <- grep("Phospho (STY) Probabilities", names(df), fixed = TRUE)
    
    if (length(phos_idx) > 1) phos_idx <- phos_idx[1]
    
    if (!purrr::is_empty(phos_idx)) {
      if (purge_phosphodata) {
        df <- df %>% 
          dplyr::filter(!nchar(as.character(.[["Phospho (STY) Probabilities"]])) == 0) 
      }
      
      phos <- df %>% 
        dplyr::select(phos_idx) %>% 
        purrr::map(~ str_extract_all(.x, "\\([^()]+\\)")) %>% 
        .[[1]] %>% 
        purrr::map(~ substring(.x, 2, nchar(.x) - 1)) %>% 
        purrr::map(as.numeric) %>% 
        purrr::map(sort, decreasing = TRUE)
      
      # may contain numeric(0) from non phospho entries
      phos_max <- suppressWarnings(phos %>% purrr::map(`[`, 1)) %>% unlist()
      sec_max <- suppressWarnings(phos %>% purrr::map(`[`, 2)) %>% unlist()
      
      df <- dplyr::bind_cols(df, 
                             pep_phospho_locprob = phos_max, 
                             pep_phospho_locdiff = phos_max - sec_max) %>% 
        dplyr::mutate(is_locprob_one = magrittr::equals(1, pep_phospho_locprob)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locdiff"), 
                         ~ replace(.x, is_locprob_one, 1)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locprob"), 
                         ~ replace(.x, is.infinite(.x), NA)) %>% 
        dplyr::select(-is_locprob_one)
    }
    
    return(df)
  })
  
  # make available `prot_acc` & `pep_seq` 
  if (TMT_plex > 0) {
    stop("TMT workflows not yet available with MSFragger.", call. = FALSE)
    
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), 
                       ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
  } else {
    df <- df %>% 
      dplyr::rename(I000 = Intensity) %>% 
      # dplyr::filter(I000 > 0) %>% 
      dplyr::mutate(R000 = I000/I000, 
                    R000 = ifelse(is.infinite(R000), NA, R000))
  }
  
  # column `RAW_File` made during pad_mf_channels
  df <- df %>% 
    dplyr::rename(
      pep_seq = Peptide, 
      prot_acc = `Protein ID`, 
      pep_expect = Expectation, 
      pep_score = Hyperscore, 
    ) %>% 
    annotPrn(fasta, entrez)

  acc_types <- unique(df$acc_type) %>% .[!is.na(.)]
  
  if (rm_craps) {
    craps <- load_craps(acc_types)
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% 
      dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) {
    df <- df %>% 
      split(., .$acc_type, drop = TRUE) %>% 
      .[acc_types] %>% 
      purrr::imap( ~ annotKin(.x, .y)) %>% 
      do.call(rbind, .)
  } else {
    df <- df %>% 
      dplyr::mutate(kin_attr = FALSE, kin_class = NA, kin_order = NA)
  }
  
  # make available `pep_res_before` and `pep_res_after`
  df <- df %>% 
    annotPeppos() %>% 
    dplyr::mutate(pep_seq_bare = gsub("^.*\\.([^\\.]+)\\..*", "\\1", pep_seq)) %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq_bare)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq_bare), 
                                    stringr::str_count(pep_seq_bare, "[KR]") - 1,
                                    stringr::str_count(pep_seq_bare, "[KR]"))) %>% 
    add_prot_icover() %>% 
    dplyr::select(-pep_seq_bare)
  
  if (!("prot_cover" %in% names(df) && length(filelist) == 1)) {
    df <- df %>% calc_cover(id = !!rlang::sym(group_pep_by))
  }
  
  # --  column arrangements --- 
  df <- df %>% reloc_col_after("prot_index", "prot_cover") 
  
  df <- df %>% 
    order_psm_cols(c("gene", "fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                     "other_acc", "entrez", "species", "acc_type", 
                     "kin_attr", "kin_class", "kin_order", "is_tryptic"))

  df <- cbind.data.frame(
    df[, grepl("^[a-z]", names(df))], 
    df[, grepl("^[A-Z]", names(df)) & !grepl("^[IR][0-9]{3}[NC]{0,1}", names(df))], 
    df[, grepl("^[R][0-9]{3}[NC]{0,1}", names(df))], 
    df[, grepl("^[I][0-9]{3}[NC]{0,1}", names(df))]
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))),
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
  )
  
  df <- df %>%
    dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), 
                     as.numeric) %>%
    dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                     ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
    dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
    dplyr::arrange(RAW_File, pep_seq, prot_acc) 
  
  # make available `pep_seq_mod` 
  stopifnot("Modified Peptide" %in% names(df))
  
  df <- df %>% dplyr::mutate(`Modified Peptide` = 
                               ifelse(is.na(`Modified Peptide`), 
                                      pep_seq, 
                                      `Modified Peptide`))
  
  df <- df %>% 
    dplyr::rename(pep_seq_mod = `Modified Peptide`) %>% 
    add_msfragger_pepseqmod(use_lowercase_aa) %>% 
    reloc_col_before("dat_file", "RAW_File") 
  
  df <- df %>% 
    reloc_col_after("pep_seq_mod", "pep_seq") %>% 
    reloc_col_before("pep_res_before", "pep_seq") %>% 
    reloc_col_after("pep_res_after", "pep_seq_mod") %>% 
    reloc_col_after("pep_start", "pep_res_after") %>% 
    reloc_col_after("pep_end", "pep_start") %>% 
    reloc_col_after("pep_len", "pep_end") %>% 
    reloc_col_after("pep_miss", "pep_len") %>% 
    reloc_col_before("pep_isunique", "pep_res_before") %>% 
    reloc_col_before("pep_literal_unique", "pep_res_before") %>% 
    reloc_col_before("pep_razor_unique", "pep_res_before") %>% 
    reloc_col_before("pep_score", "pep_res_before") %>% 
    reloc_col_before("pep_expect", "pep_res_before")

  # --- compare experimental and user-provided PSM files ---
  res <- check_raws(df)
  df <- res$df
  tmtinj_raw_map <- res$lookup
  rm(res)
  
  # --- split by TMT and LCMS ---
  if (! "PSM_File" %in% names(tmtinj_raw_map)) {
    df_split <- df %>%
      dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>%
      split(., .$TMT_inj, drop = TRUE)
  } else {
    tmtinj_raw_map <- tmtinj_raw_map %>% 
      dplyr::mutate(PSM_File = gsub("\\.tsv$", "", PSM_File)) %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "PSM_File"), sep = "@") 
    
    df_split <- df %>% 
      tidyr::unite(RAW_File2, c("RAW_File", "dat_file"), 
                   sep = "@", remove = FALSE) %>% 
      dplyr::left_join(tmtinj_raw_map, by = "RAW_File2") %>% 
      dplyr::select(-RAW_File2) %>% 
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>%
      dplyr::select(-dat_file) %>% 
      split(., .$TMT_inj, drop = TRUE)
  }
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    warning("TMT sets and LC/MS injections not have corresponindg PSM files:\n", 
            call. = FALSE)
    message(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    stop("Remove mismatched `TMT_Set` and/or `LC/MS` ",
         "from experiment summary file.", 
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), 
                                   sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")
  
  df_split <- df_split %>% .[names(.) %>% sort_tmt_lcms()]
  
  if (parallel) {
    n_cores <- pmin(parallel::detectCores(), length(filelist))
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    suppressWarnings(
      silent_out <- parallel::clusterMap(
        cl, psm_msplit, df_split, names(df_split), 
        MoreArgs = list(fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
      )
    )
    
    parallel::stopCluster(cl)
  } else {
    purrr::walk2(df_split, names(df_split), psm_msplit, 
                 fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
  }
}


theme_psm_violin <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 24),
    axis.ticks.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, vjust = 0.5, size = 24),
    axis.title.x = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    plot.title = element_text(colour = "black", size  =24, hjust = .5, vjust = .5),
    
    strip.text.x = element_text(size = 18, colour = "black", angle = 0),
    strip.text.y = element_text(size = 18, colour = "black", angle = 90),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    legend.key = element_rect(colour = NA, fill = 'transparent'),
    legend.background = element_rect(colour = NA,  fill = "transparent"),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 18),
    legend.text.align = 0,
    legend.box = NULL
  )







#' Annotates MSFragger PSM results
#'
#' \code{annotPSM_mf} adds fields of annotation to MSFragger PSM tables.
#'
#' @inheritParams load_expts
#' @inheritParams annotPSM
#' @inheritParams normPSM
#' @param ... Not currently used.
annotPSM_mf <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        mc_psm_by = "peptide", 
                        fasta = NULL, expt_smry = "expt_smry.xlsx", 
                        plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                        ...) {
  
  dat_dir <- get_gl_dat_dir()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM/cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files()
  
  set_indexes <- gsub("^TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort() 
  
  for (set_idx in set_indexes) {
    sublist <- filelist[grep(paste0("^TMTset", set_idx, "_"), 
                             filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(
      Filename = do.call('rbind', strsplit(as.character(sublist), '.txt', fixed = TRUE))
    ) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    # --- LCMS injections under the same TMT experiment --- 
    for (idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM/cache", sublist[idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") 
      
      if (TMT_plex > 0) {
        df <- local({
          injn_idx <- sublist[idx] %>% 
            gsub("^TMTset\\d+_LCMSinj(\\d+)_Clean.txt$", "\\1", .) %>% 
            as.integer()
          
          df <- mcPSM(df, set_idx, injn_idx, mc_psm_by, group_psm_by, group_pep_by)
        })
      } else {
        df <- df %>% 
          dplyr::mutate(N_I000 = I000, R000 = NA, log2_R000 = NA, N_log2_R000 = NA)
      }
      
      df <- df %>% 
        calcSD_Splex(group_psm_by) %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
        dplyr::right_join(df, by = group_psm_by) %>% 
        dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 0)) %>% 
        dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
                         ~ round(.x, digits = 4)) %>% 
        na_genes_by_acc() 
      
      pep_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by)) %>%
        dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
        dplyr::summarise(pep_n_psm = n())
      
      prot_n_psm <- df %>% 
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_psm = n())
      
      prot_n_pep <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_pep = n())
      
      df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)
      
      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))),
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
      ) %>% 
        reloc_col_before("pep_seq", "pep_res_after") %>% 
        reloc_col_after("pep_seq_mod", "pep_seq") 

      df <- dplyr::bind_cols(
        df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^I[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))), 
        df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))),
      ) %>% 
        dplyr::select(-which(names(.) %in% c("pep_index", "prot_index"))) %T>% 
        write.table(file.path(dat_dir, "PSM", paste0(out_fn[idx, 1], ".txt")),
                         sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int) {
        if (TMT_plex > 0) {
          df_int <- df %>% .[, grepl("^N_I[0-9]{3}", names(.))]
        } else {
          df_int <- df %>% 
            dplyr::select("I000") %>% 
            dplyr::rename(N_I000 = "I000")
        }
        
        rptr_violin(
          df = df_int, 
          filepath = file.path(dat_dir, "PSM/rprt_int/mc", 
                               paste0(gsub("_PSM_N", "", out_fn[idx, 1]), "_rprt.png")), 
          width = 8, height = 8)
      }
      
      if (plot_log2FC_cv && TMT_plex > 0) {
        sd_violin(
          df = df, id = !!group_psm_by, 
          filepath = file.path(dat_dir, "PSM/log2FC_cv/raw", 
                               paste0(gsub("_PSM_N", "", out_fn[idx, 1]), "_sd.png")), 
          width = 8, height = 8, type = "log2_R", 
          adjSD = FALSE, is_psm = TRUE, , ...)
      }
    }
  }
}

