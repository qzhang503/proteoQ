#' Extract RAW MS file names
#'
#' Extract a list of \code{RAW} file names that can be passed to \code{frac_smry.xlsx}
#'
#' @examples
#' \dontrun{
#' # Supposed that RAW MS files are stored under "~\my_raw"
#' extract_raws("~\\my_raw")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
#' @importFrom tools md5sum
#' @export
extract_raws <- function(raw_dir) {
  dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
  if (dat_dir == 1) 
    stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)

  fns <- names(tools::md5sum(dir(raw_dir, pattern = "\\.raw$", full.names = FALSE)))
  data.frame(Index = seq_along(fns), RAW_File = fns) %T>% 
    write.table(file.path(dat_dir, "raw_list.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
		
	message("RAW MS file names in ", file.path(dat_dir, "raw_list.txt"))
}


#' Extract RAW file names from Mascot PSM outputs
#'
#' \code{extract_psm_raws} extracts the RAW file names from the PSM data under
#' the current working directory.
#'
#' @param type Character string indicating the type of PSM.
#' @inheritParams load_expts
#' @examples
#' extract_psm_raws(mascot)
#' 
#' extract_psm_raws(maxquant)
#' 
#' extract_psm_raws(spectrum_mill)
#'
#' @import dplyr tidyr rlang
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
extract_psm_raws <- function(type = c("mascot", "maxquant", "spectrum_mill"), dat_dir = NULL) {
  batchPSMheader_2 <- function(filelist, TMT_plex) {
    df <- readLines(file.path(dat_dir, filelist))
    
    pep_seq_rows <- grep("pep_seq", df)
    
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
    writeLines(df, file.path(dat_dir, "PSM\\temp", paste0(output_prefix, "_hdr_rm.csv")))
  }
    
  find_mascot_psmraws <-function() {
    dir.create(file.path(dat_dir, "PSM\\temp"), recursive = TRUE, showWarnings = FALSE)
    
    purrr::walk(filelist, batchPSMheader_2, TMT_plex)
    
    df <- purrr::map(gsub("\\.csv$", "_hdr_rm.csv", filelist), ~ {
      read.delim(file.path(dat_dir, "PSM\\temp", .x), sep = ',', check.names = FALSE, 
                 header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    }) %>% 
      do.call(rbind, .)
    
    r_start <- which(names(df) == "pep_scan_title") + 1
    int_end <- ncol(df)
    if(int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
    
    if (TMT_plex == 11) {
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
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title))
    
    raws <- unique(df$RAW_File)
    
    if (!purrr::is_empty(raws)) {
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, "mascot_psm_raws.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      message("MS file names in ", file.path(dat_dir, "mascot_psm_raws.txt"))
    }
    
    unlink(file.path(dat_dir, "PSM\\temp"), recursive = TRUE, force = TRUE)
  } 
  
  find_sm_psmraws <- function() {
    df <- purrr::map(file.path(dat_dir, filelist), readr::read_delim, delim = ";") %>% 
      dplyr::bind_rows() %>% 
      dplyr::rename(RAW_File = `filename`) %>% 
      dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", RAW_File))
    
    raws <- unique(df$RAW_File)
    
    if (!purrr::is_empty(raws)) {
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, "sm_psm_raws.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      message("MS file names in ", file.path(dat_dir, "sm_psm_raws.txt"))
    }
  }
  
  find_mq_psmraws <- function() {
    df <- purrr::map(file.path(dat_dir, filelist), read.csv, 
                     check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
      dplyr::bind_rows() %>% 
      dplyr::rename(RAW_File = `Raw file`) 

    raws <- unique(df$RAW_File)
    
    if (!purrr::is_empty(raws)) {
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, "mq_psm_raws.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      message("MS file names in ", file.path(dat_dir, "mq_psm_raws.txt"))
    }
  }
  
  type <- rlang::enexpr(type)
  if (type == rlang::expr(c("mascot", "maxquant", "spectrum_mill"))) {
    type <- "mascot"
  } else {
    type <- rlang::as_string(type)
  }
  
  type <- tolower(type)
  if (type == "ms") type <- "mascot"
  if (type == "mq") type <- "maxquant"
  if (type == "sm") type <- "spectrum_mill"

  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Assign the working directory to variable `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }

  pattern <- switch(type, 
    mascot = "^F[0-9]+\\.csv$", 
    maxquant = "^msms.*\\.txt$",
    spectrum_mill = "^PSMexport.*\\.ssv$", 
    stop("Data type needs to be one of `mascot`, `maxquant` or `spectrum_mill`.", Call. = FALSE)
  )
  
  filelist <- list.files(path = file.path(dat_dir), pattern = pattern)
  if (purrr::is_empty(filelist)) stop("No ", toupper(type), " files(s) under ", dat_dir, call. = FALSE)

  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme)

  switch (type,
    mascot = find_mascot_psmraws(),
    maxquant = find_mq_psmraws(),
    spectrum_mill = find_sm_psmraws(),
  )
}


#' Removes PSM headers
#'
#' \code{rmPSMHeaders} removes the header of PSM from
#' \code{\href{https://http://www.matrixscience.com/}{Mascot}} outputs. It also
#' removes the spacer columns in the fields of ratio and intensity values.
#'
#' @return Intermediate PSM table(s).
#'
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>%
rmPSMHeaders <- function () {
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)
	options(max.print = 5000000)

	on.exit(message("Remove PSM headers --- Completed."), add = TRUE)

	filelist = list.files(path = file.path(dat_dir), pattern = "^F[0-9]+\\.csv$")

	if (purrr::is_empty(filelist))
	  stop("No PSM files(s) with `.csv` extension under ", dat_dir, call. = FALSE)

	load(file = file.path(dat_dir, "label_scheme.rda"))
	TMT_plex <- TMT_plex(label_scheme)

	batchPSMheader <- function(filelist, TMT_plex) {
		data_all <- readLines(file.path(dat_dir, filelist))

		pep_seq_rows <- grep("pep_seq", data_all)
		data_header <- data_all[1 : (pep_seq_rows[1] - 1)]
		data_header <- gsub("\"", "", data_header, fixed = TRUE)

		output_prefix <- gsub("\\.csv$", "", filelist)
		write.table(data_header, file.path(dat_dir, "PSM\\cache",
		            paste0(output_prefix, "_header", ".txt")),
		            sep = "\t", col.names = FALSE, row.names = FALSE)
		rm(data_header)
		
		if (length(pep_seq_rows) > 1) {
		  data_queries <- data_all[pep_seq_rows[2] : length(data_all)]
		  data_queries <- gsub("\"", "", data_queries, fixed = TRUE)
		  writeLines(data_queries, file.path(dat_dir, "PSM\\cache", paste0(output_prefix, "_queries.csv")))
		  rm(data_queries)
		}

		unassign_hits_row <- grep("Peptide matches not assigned to protein hits", data_all)
		if (! purrr::is_empty(unassign_hits_row)) {
		  psm_end_row <- unassign_hits_row - 2
		} else if (length(pep_seq_rows) > 1) {
		  psm_end_row <- pep_seq_rows[2] - 4
		} else {
		  psm_end_row <- length(data_all)
		}
	
		data_psm <- data_all[pep_seq_rows[1] : psm_end_row]
		data_psm <- gsub("\"---\"", -1, data_psm, fixed = TRUE)
		data_psm <- gsub("\"###\"", -1, data_psm, fixed = TRUE)
		
		if (TMT_plex > 0) data_psm[1] <- paste0(data_psm[1], paste(rep(",", TMT_plex * 4 -2), collapse = ''))
		
		writeLines(data_psm, file.path(dat_dir, "PSM\\cache", paste0(output_prefix, "_hdr_rm.csv")))
		rm(data_psm)
	}

	purrr::walk(filelist, batchPSMheader, TMT_plex)
}


#' Add the `pep_seq_mod` field to Mascot PSMs
#'
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
add_mascot_pepseqmod <- function(df, use_lowercase_aa) {
  dat_id <- df$dat_file %>% unique()
  dat_file <- file.path(dat_dir, "PSM\\cache", paste0(dat_id, "_header.txt"))
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
    tidyr::separate(".", sep = ",", extra = "drop", c("Mascot_abbr", "Description", "Delta_mass")) %>%
    dplyr::mutate(Filename = gsub("[\\\\\\/\\:\\*\\?\\'\\<\\>\\|]", ".", Description))
  
  if (is.null(df$pep_seq)) stop("column `pep_seq` not found.")

  if (nrow(var_mods) == 0) {
    df$pep_seq_mod <- df$pep_seq
    return(df)
  } 

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
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
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
              
              lowers <- substr(df_sub$pep_seq_mod[rows], locales, locales) %>% tolower()
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
        
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)

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
          df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "_"))
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
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0("~", pep_seq_mod))
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
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "~"))
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
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = gsub("(^[_~]{0,1})(.)", paste0("\\1", "^", "\\2"), pep_seq_mod)) 
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
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
          
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = gsub("(.)([_~]{0,1}$)", paste0("\\1", "^", "\\2"), pep_seq_mod)) 
            
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })

    # (5) paste "pep_res_before" and "pep_res_after"
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = "."))
  }

  purrr::walk2(var_mods$Description, var_mods$Filename, ~ {
    try(
      df %>% 
        dplyr::filter(grepl(.x, pep_var_mod, fixed = TRUE)) %>% 
        write.table(file.path(dat_dir, "PSM\\individual_mods", paste0(.y, ".txt")), 
                    sep = "\t", col.names = TRUE, row.names = FALSE)
    )
  })

  return(df)
}


#' Splits PSM tables
#'
#' \code{splitPSM} splits the PSM outputs after \code{rmPSMHeaders()}. It
#' separates PSM data by TMT experiment and LC/MS injection.
#'
#' Arguments \code{group_psm_by} and \code{group_pep_by} are used to update
#' \code{prot_matches_sig} and \code{prot_sequences_sig} after data merge.
#'
#' @param fasta Character string(s) to the name(s) of fasta file(s) with
#'   prepended directory path. The \code{fasta} database(s) need to match those
#'   used in MS/MS ion search. There is no default and users need to provide the
#'   correct file path(s) and name(s).
#' @param rm_craps Logical; if TRUE,
#'   \code{\href{https://www.thegpm.org/crap/}{cRAP}} proteins will be removed.
#'   The default is FALSE.
#' @param rm_krts Logical; if TRUE, keratin entries will be removed. The default
#'   is FALSE.
#' @param annot_kinases Logical; if TRUE, proteins of human or mouse origins
#'   will be annotated with their kinase attributes. The default is FALSE.
#' @param rptr_intco Numeric; the threshold of reporter ion intensity. The
#'   default is 1,000.
#' @param plot_rptr_int Logical; if TRUE, the distributions of reporter-ion
#'   intensities will be plotted. The default is TRUE.
#' @param use_lowercase_aa Logical; if TRUE, modifications in amino acid
#'   residues will be abbreviated with lower-case and/or \code{^_~}. See the
#'   table below for details. The default is TRUE. 
#' @import dplyr tidyr seqinr stringr
#' @importFrom magrittr %>%
splitPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, 
                     rm_craps = FALSE, rm_krts = FALSE, rptr_intco = 1000, 
                     annot_kinases = FALSE, plot_rptr_int = TRUE, use_lowercase_aa = TRUE, ...) {

	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)
	on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)

	load(file = file.path(dat_dir, "label_scheme_full.rda"))
	load(file = file.path(dat_dir, "label_scheme.rda"))
	load(file = file.path(dat_dir, "fraction_scheme.rda"))

	TMT_plex <- TMT_plex(label_scheme_full)

  filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
                        pattern = "^F[0-9]{6}\\_hdr_rm.csv$")

	if (length(filelist) == 0) stop(paste("No PSM files under", file.path(dat_dir, "PSM")))

  df <- purrr::map(filelist, ~ {
    data <- read.delim(file.path(dat_dir, "PSM\\cache", .x), sep = ',', check.names = FALSE, 
                       header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    data$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    return(data)
  }) %>% 
    do.call(rbind, .)
  
  r_start <- which(names(df) == "pep_scan_title") + 1
  int_end <- ncol(df) - 1
	if (int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]

  if (TMT_plex == 11) {
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

	if (TMT_plex > 0) {
		colnames(df)[r_start:(r_start+TMT_plex-2)] <- col_ratio
		colnames(df)[(r_start+TMT_plex-1):(r_start+TMT_plex+TMT_plex-2)] <- col_int
		rm(r_start, int_end, col_ratio, col_int)
	}
  
  # convenience craps removals where their uniprot afasta names ended with "|"
  if (rm_craps) df <- df %>% dplyr::filter(!grepl("\\|.*\\|$", prot_acc))

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  # message("Primary column keys in `PSM/cache/[...]_hdr_rm.csv` for `filter_` varargs.")
  message("Primary column keys in `F[...].csv` or `PSM/cache/[...]_hdr_rm.csv`  for `filter_` varargs.")

  # note pep_seq: from such as MENGQSTAAK to K.MENGQSTAAK.L
  df <- df %>% 
    dplyr::mutate(pep_len = str_length(pep_seq)) %>% 
    split(., .$dat_file, drop = TRUE) %>% 
    purrr::map(add_mascot_pepseqmod, use_lowercase_aa) %>% 
    bind_rows() %>% 
    dplyr::select(-dat_file)

  df <- df %>% 
    dplyr::mutate(prot_acc_orig = prot_acc) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc)) %>% 
    annotPrn(fasta) %>% 
    dplyr::mutate(prot_acc = prot_acc_orig) %>% 
    dplyr::select(-prot_acc_orig)
  
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
    purrr::reduce(left_join, by = group_pep_by) %>% 
    dplyr::mutate(prot_matches_sig = prot_matches_sig_new, 
                  prot_sequences_sig = prot_sequences_sig_new) %>%
    dplyr::select(-prot_matches_sig_new, -prot_sequences_sig_new)
  
  rm(prot_matches_sig, prot_sequences_sig)
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc))

  # re-apply craps after annotation
  # 'acc_type' will be NA for entries not found in fasta
  acc_type <- unique(df$acc_type) %>% .[!is.na(.)]

  stopifnot(length(acc_type) == 1)
  
  if (rm_craps) {
    data(package = "proteoQ", prn_annot_crap)

    craps <- prn_annot_crap %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::select(acc_type) %>% 
      unlist()
    
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }

  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) df <- annotKin(df, acc_type)
  
  # `pep_start`, `pep_end` and `gene` will be used for protein percent coverage calculation
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  } 
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(-grep("^prot_", names(.))), 
  )

  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>%
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>%
      # a special case of redundant entries from Mascot
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>% 
			dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` from experimental summary file."),
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

	for (i in seq_along(df_split)) {
		df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)

		out_fn <- fn_lookup %>%
			dplyr::filter(TMT_inj == names(df_split)[i]) %>%
			dplyr::select(filename) %>%
			unique() %>%
			unlist() %>%
			paste0(., ".csv")

		df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
		
		write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
		
		if (plot_rptr_int & TMT_plex > 0) {
		  df_int <- df_split[[i]] %>% 
		    .[, grepl("^I[0-9]{3}", names(.))]
		  
		  rptr_violin(df = df_int, filepath = file.path(dat_dir, "PSM\\rprt_int\\raw", gsub("\\.csv", "\\.png", out_fn)), 
		              width = 8, height = 8)
		}
	}
}


#' Locates the positions of outliersf
#' 
#' @param df A data frame containing the PSM table from Mascot.
#' @param range_colRatios The range of columns.
#' @return A data frame.
#' @examples
#' locate_outliers(df, 2:3)
locate_outliers <- function (df, range_colRatios) {
	for(col_index in range_colRatios) {
		counts <- colSums(!is.na(df[col_index]))
		if(counts > 25) df[, col_index] <- Rosner_outliers(df[, col_index])
		else if(counts > 2) df[, col_index] <- Dixon_outliers(df[, col_index])
	}

	return(df)
}


#' Outlier removals with Rosner's method
Rosner_outliers <- function(x) {

	if (length(unique(x)) < 5) return(x)
	# up to 9-number outliers; may get warnings with NA being an outlier
	gofOutlier_obj <- rosnerTest(as.numeric(x), 9)

	if(gofOutlier_obj$n.outliers > 0) {
		Index <- with(gofOutlier_obj$all.stat, Obs.Num[Outlier == TRUE])
		x[Index] <- NA
	}

	return(x)
}


#' Outlier removals with Dixon's method
Dixon_outliers <- function(x) {
	# x = c(0.0000000, 0.0000000, 1.0271542, 0.0000000, 0.2080097)
	# x = c(0.0000000, 0.0000000, NA, 0.0000000, 0.2080097)
	# x = c(0.0000000, 0.0000000, 0.0000000, 0.2080097)
	# x = c(NA, NA, NA, 0.2080097)

	newx <- x[!is.na(x)]
	len_newx <- length(newx)
	uni_newx <- length(unique(newx))

	if (len_newx > 2 & uni_newx > 1) {
		gofOutlier_obj <- dixon.test(as.numeric(x), type = 0)

		while(gofOutlier_obj$p.value < 0.05) {
			if (grepl("^highest", gofOutlier_obj$alternative)) x[which.max(x)] <- NA else x[which.min(x)] <- NA

			newx <- x[!is.na(x)]
			len_newx <- length(newx)
			uni_newx <- length(unique(newx))

			if (len_newx > 2 & uni_newx > 1)
			  gofOutlier_obj <- dixon.test(as.numeric(x), type = 0) else gofOutlier_obj$p.value <- 1
		}
	}

	return (x)
}


#' Outlier removals with Grubbs's method
Grubbs_outliers <- function(x, type = 10) {
	newx <- x[!is.na(x)]
	len_newx <- length(newx)
	uni_newx <- length(unique(newx))

	if (len_newx > 2 & uni_newx > 1) {
		gofOutlier_obj <- grubbs.test(as.numeric(x, type))

		while(gofOutlier_obj$p.value < 0.05) {
			if (grepl("^highest", gofOutlier_obj$alternative)) x[which.max(x)] <- NA else x[which.min(x)] <- NA

			newx <- x[!is.na(x)]
			len_newx <- length(newx)
			uni_newx <- length(unique(newx))

			if (len_newx > 2 & uni_newx > 1)
			  gofOutlier_obj <- grubbs.test(as.numeric(x), type = type) else gofOutlier_obj$p.value <- 1
		}
	}

	return (x)
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
#'
#' @import dplyr tidyr
#' @importFrom stringr str_split
#' @importFrom outliers dixon.test
#' @importFrom EnvStats rosnerTest
cleanupPSM <- function(rm_outliers = FALSE) {
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 5000000)

	load(file = file.path(dat_dir, "label_scheme.rda"))
	load(file = file.path(dat_dir, "label_scheme_full.rda"))
	TMT_plex <- TMT_plex(label_scheme_full)

	filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
	                      pattern = "^TMT.*LCMS.*\\.csv$")

	for (i in seq_along(filelist)) {
		df <- read.csv(file.path(dat_dir, "PSM\\cache", filelist[i]), check.names = FALSE,
		               header = TRUE, comment.char = "#")

		if (TMT_plex == 0) {
			# lable-free data
		  # re-save ".csv" as ".txt"
			fn <- paste0(gsub(".csv", "", filelist[i]), "_Clean.txt")
			write.table(df, file.path(dat_dir, "PSM\\cache", fn), sep = "\t", col.names = TRUE,
			            row.names = FALSE)
			cat(filelist[i], "processed\n")

			next
		}

		# remove all "-1" ratio rows
		N <- sum(grepl("^R[0-9]{3}", names(df)))

		df <- df %>%
			dplyr::mutate(n = rowSums(.[, grep("^R[0-9]{3}", names(.))] == -1)) %>%
			dplyr::filter(n != N) %>%
			dplyr::select(-n)

		channelInfo <- channelInfo(label_scheme, set_idx =
		                  as.integer(gsub("TMTset(\\d+)_.*", "\\1", filelist[i])))

		# add a column of "R126"
		pos_af <- min(grep("^R1[0-9]{2}", names(df)))

		df$R126 <- 1
		df <- cbind.data.frame(df[, 1:(pos_af-1)], R126 = df$R126, df[, (pos_af):(ncol(df)-1)]) %>%
				dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^I[0-9]{3}", names(.) )])) > 0) %>%
				dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
				dplyr::mutate_at(.vars = which(names(.)=="R126")-1+channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) # note that "> 1" not "0"

		if (rm_outliers) {
			dfw_split <- df %>%
				dplyr::select(grep("^I[0-9]{3}", names(.))) %>%
				dplyr::mutate(RM = rowMeans(.[, grep("^I[0-9]{3}", names(.))[channelInfo$labeledChannels]],
				                            na.rm = TRUE)) %>%
				dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ log2(.x/RM)) %>%
				dplyr::select(-c("RM")) %>%
				`colnames<-`(gsub("I", "X", names(.))) %>%
				dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA)) %>%
				dplyr::bind_cols(df[, c("psm_index", "pep_seq")], .) %>%
				split(., .$pep_seq, drop = TRUE)

			range_colRatios <- grep("^X[0-9]{3}", names(dfw_split[[1]]))

			dfw_split <- do.call("rbind", lapply(dfw_split, locate_outliers, range_colRatios)) %>%
					dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA)) %>%
					tidyr::unite(pep_seq_i, pep_seq, psm_index, sep = ":") %>%
					dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, !is.na(.x), 1))

			df <- df %>%
					tidyr::unite(pep_seq_i, pep_seq, psm_index, sep = ":") %>%
					dplyr::left_join(., dfw_split, by = "pep_seq_i") %>%
					tidyr::separate(pep_seq_i, into = c("pep_seq", "psm_index"), sep = ":", remove = TRUE) %>%
					dplyr::select(-c("psm_index"))

			rm(dfw_split, range_colRatios)

			df[, grepl("^I[0-9]{3}", names(df))] <-
			  purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
			              as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
			  dplyr::bind_rows()

			df[, grepl("^R[0-9]{3}", names(df))] <-
			  purrr::map2(as.list(df[, grepl("^R[0-9]{3}", names(df))]),
			              as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
			  dplyr::bind_rows()

			df <- cbind.data.frame(raw_file = df[, c("raw_file")],
					df[, !grepl("^R[0-9]{3}|^I[0-9]{3}|^X[0-9]{3}|^raw_file$", names(df))],
					df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df))]) %>%
					dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) %>% # "> 1" as "R126 == 1"
					dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 + channelInfo$emptyChannels, ~ replace(.x, , 0))
		} else {
			df <- cbind.data.frame(raw_file = df[, c("raw_file")],
				df[, !grepl("^R[0-9]{3}|^I[0-9]{3}|^psm_index$|^raw_file$", names(df))],
				df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df))]) %>%
			  dplyr::mutate_at(.vars = which(names(.)=="I126") - 1 +
				                   channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
				dplyr::mutate_at(.vars = which(names(.)=="R126") - 1 +
				                   channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1)

		}

		fn <- paste0(gsub(".csv", "", filelist[i]), "_Clean.txt")
		write.table(df, file.path(dat_dir, "PSM\\cache", fn), sep = "\t", col.names = TRUE,
		            row.names = FALSE)
		cat(filelist[i], "processed\n")
	}

}


#'Median-centering normalization of PSM data
#'
#'\code{mcPSM} adds fields of \code{log2_R, N_log2_R and N_I} to PSM tables.
#'
#'@import dplyr tidyr purrr
#'@importFrom magrittr %>%
mcPSM <- function(df, set_idx) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  label_scheme_sub <- label_scheme[label_scheme$TMT_Set == set_idx & 
                                     label_scheme$LCMS_Injection == 1, ]
  
  channelInfo <- channelInfo(label_scheme_sub, set_idx)
  
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
    `colnames<-`(gsub("I", "log2_R", names(.)))	%>%
    cbind(dfw, .) %>%
    dplyr::mutate_at(.vars = grep("[I|R][0-9]{3}", names(.)), ~ replace(.x, is.infinite(.), NA))
  
  col_log2Ratio <- grepl("^log2_R[0-9]{3}", names(dfw))
  cf <- apply(dfw[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  
  dfw <- sweep(dfw[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- sweep(dfw[, grepl("^I[0-9]{3}", names(dfw)), drop = FALSE], 2, 2^cf, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- dfw %>%
    reorderCols(endColIndex = grep("[RI][0-9]{3}", names(dfw)), col_to_rn = "pep_seq_mod") %>%
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
#'@param plot_log2FC_cv Logical; if TRUE, the distributions of the CV of peptide
#'  \code{log2FC} will be plotted. The default is TRUE.
#'@inheritParams load_expts
#'@inheritParams splitPSM
#'@import dplyr tidyr purrr ggplot2 RColorBrewer
#'@importFrom stringr str_split
#'@importFrom tidyr gather
#'@importFrom magrittr %>%
annotPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                     fasta = NULL, expt_smry = "expt_smry.xlsx", 
                     plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 5000000)
  
  hd_fn <- list.files(path = file.path(dat_dir, "PSM\\cache"),
                      pattern = "^F\\d+_header.txt$")
  assign("df_header", readLines(file.path(dat_dir, "PSM\\cache", hd_fn[1])))
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(n_TMT_sets)
  
  for (set_idx in seq_len(n_TMT_sets)) {
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename =
                           do.call('rbind', strsplit(as.character(sublist),
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    # LCMS injections under the same TMT experiment
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#")

      acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
      stopifnot(length(acc_type) == 1)
      
      species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()

      if (TMT_plex > 0) df <- mcPSM(df, set_idx)
      
      df <- df %>% 
        calcSD_Splex(id = group_psm_by, type = "log2_R") %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
        dplyr::right_join(df, by = group_psm_by) %>% 
        dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
        dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 4))
      
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
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^pep_", names(.))), 
      )
      
      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))), 
        df %>% dplyr::select(-grep("^prot_", names(.))), 
      )
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
      )
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int & TMT_plex > 0) {
        df_int <- df %>% .[, grepl("^N_I[0-9]{3}", names(.))]
        rptr_violin(df = df_int, 
                    filepath = file.path(dat_dir, "PSM\\rprt_int\\mc",
                                         paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_rprt.png")), 
                    width = 8, height = 8)
      }

      if (plot_log2FC_cv & TMT_plex > 0) {
        sd_violin(df = df, id = !!group_psm_by, 
                  filepath = file.path(dat_dir, "PSM\\log2FC_cv\\raw", 
                                       paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_sd.png")), 
                  width = 8, height = 8, type = "log2_R", adjSD = FALSE, is_psm = TRUE)
      }
      
    }
    
  }
}


#'Standardization of PSM results
#'
#'\code{normPSM} standarizes
#'\code{\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSM}}
#'results from \code{\href{https://en.wikipedia.org/wiki/Tandem_mass_tag}{TMT}}
#'experiments.
#'
#'In each primary output file, "\code{...PSM_N.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'\code{log2_R...} with median-centering alignment. Values under columns
#'\code{I...} are raw \code{reporter-ion intensity} from database searches.
#'Values under columns \code{N_I...} are normalized \code{reporter-ion
#'intensity}. Values under columns \code{sd_log2_R...} are the standard
#'deviation of the \code{log2FC} of peptides from ascribing PSMs. Character
#'strings under \code{pep_seq_mod} denote peptide sequences with applicable
#'variable modifications.
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
#'  \code{\href{https://http://www.matrixscience.com/}{Mascot}} at a \code{.csv}
#'  format and store them under the file folder indicated by \code{dat_dir}. The
#'  header information should be included during the \code{.csv} export. The
#'  file name(s) should be defaulted by
#'  \code{\href{https://http://www.matrixscience.com/}{Mascot}}: starting with
#'  the letter \code{'F'}, followed by digits without space and ended with a
#'  \code{'.csv'} extension \code{(e.g., F004453.csv)}.
#'
#'  See \code{\link{normPrn}} for the description of column keys in the output.
#'
#'@section \code{MaxQuant}: End users will copy over \code{msms.txt} file(s)
#'  from \code{\href{https://www.maxquant.org/}{MaxQuant}} to the \code{dat_dir}
#'  directory. In the case of multiple \code{msms.txt} files for processing, the
#'  file names need to be compiled in that they all start with \code{'msms'} and
#'  end with a \code{'.txt'} extension.
#'
#'@section \code{Spectrum Mill}: End users will copy over \code{PSMexport.1.ssv}
#'  file(s) from
#'  \code{\href{https://www.agilent.com/en/products/software-informatics/masshunter-suite/masshunter-for-life-science-research/spectrum-mill}{Spectrum
#'   Mill}} to the \code{dat_dir} directory. In the case of multiple
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
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@section \code{Variable arguments and data files}: Variable argument (vararg)
#'  statements of \code{filter_} and \code{arrange_} are available in
#'  \code{proteoQ} for flexible filtration and ordering of data rows, via
#'  functions at users' interface. To take advantage of the feature, users need
#'  to be aware of the column keys in input files. As indicated by their names,
#'  \code{filter_} and \code{filter2_} perform row filtration against column
#'  keys from a primary data file, \code{df}, and secondary data file(s),
#'  \code{df2}, respectively. The same correspondance is applicable for
#'  \code{arrange_} and \code{arrange2_} varargs. \cr \cr Users will typically
#'  employ either primary or secondary vararg statements, but not both. In the
#'  more extreme case of \code{gspaMap(...)}, it links \code{\link{prnGSPA}}
#'  findings in \code{df2} to the significance \code{pVals} and abundance fold
#'  changes in \code{df} for volcano plot visualizaitons by gene sets. The table
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
#'  TMTset2_LCMSinj1_PSM_N.txt, ...} The indeces of TMT experiment and LC/MS
#'  injection are indicated in the file names.
#'@example inst/extdata/examples/normPSM_.R
#'@import rlang dplyr purrr ggplot2 RColorBrewer
#'@importFrom stringr str_split
#'@importFrom magrittr %>%
#'@export
normPSM <- function(group_psm_by = c("pep_seq", "pep_seq_mod"), group_pep_by = c("prot_acc", "gene"), 
                    dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx", 
                    fasta = NULL, pep_unique_by = "group", corrected_int = TRUE, rm_reverses = TRUE, 
                    rptr_intco = 1000, rm_craps = FALSE, rm_krts = FALSE, rm_outliers = FALSE, 
                    annot_kinases = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                    use_lowercase_aa = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  options(max.print = 2000000)
  on.exit(options(old_opt), add = TRUE)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("normPSM"), add = TRUE)
  
  dots <- rlang::enexprs(...)

  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to the variable first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
    message("Variable `dat_dir` added to the Global Environment.")
  }
  
  if (is.null(fasta)) stop("Path(s) to fasta file(s) not found.", call. = FALSE)
  
  group_psm_by <- rlang::enexpr(group_psm_by)
  if (group_psm_by == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
    group_psm_by <- "pep_seq"
  } else {
    group_psm_by <- rlang::as_string(group_psm_by)
    stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  }
  
  group_pep_by <- rlang::enexpr(group_pep_by)
  if (group_pep_by == rlang::expr(c("prot_acc", "gene"))) {
    group_pep_by <- "prot_acc"
  } else {
    group_pep_by <- rlang::as_string(group_pep_by)
    stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  }

  dir.create(file.path(dat_dir, "PSM\\cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\rprt_int\\raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\rprt_int\\mc"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\purged"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\individual_mods"), recursive = TRUE, showWarnings = FALSE)
  
  if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^F[0-9]+\\.csv$"))) {
    type <- "mascot"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$"))) {
    type <- "mq"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$"))) {
    type <- "sm"
  } else {
    stop("Unknow data type or missing data files.", call. = FALSE)
	}

  pep_unique_by <- rlang::as_string(rlang::enexpr(pep_unique_by))

  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))
  reload_expts()
  
  if (type == "mascot") {
    rmPSMHeaders()
    splitPSM(group_psm_by, group_pep_by, fasta, rm_craps, rm_krts, rptr_intco, 
             annot_kinases, plot_rptr_int, use_lowercase_aa, ...)
    cleanupPSM(rm_outliers)
    annotPSM(group_psm_by, group_pep_by, fasta, expt_smry, plot_rptr_int, plot_log2FC_cv, ...)
  } else if (type == "mq") {
    splitPSM_mq(group_psm_by, group_pep_by, fasta, pep_unique_by, corrected_int, 
                rptr_intco, rm_craps, rm_reverses, annot_kinases, plot_rptr_int, ...)
    cleanupPSM(rm_outliers)
		annotPSM_mq(group_psm_by, group_pep_by, fasta, expt_smry, rm_krts, 
		            plot_rptr_int, plot_log2FC_cv, use_lowercase_aa, ...)
  } else if (type == "sm") {
    splitPSM_sm(group_psm_by, group_pep_by, fasta, rm_craps, rm_krts, rptr_intco, 
                annot_kinases, plot_rptr_int, ...)
    cleanupPSM(rm_outliers)
    annotPSM_sm(group_psm_by, group_pep_by, fasta, expt_smry, rm_krts, plot_rptr_int, 
                plot_log2FC_cv, use_lowercase_aa, ...)
  }
}



#' Calculate peptide data for individual TMT experiments
#' Argument injn_idx not used
calcPepide <- function(df, label_scheme, id, method_psm_pep, group_pep_by, set_idx, injn_idx) {
  stopifnot("prot_acc" %in% names(df))
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  channelInfo <- label_scheme %>%
    dplyr::filter(TMT_Set == set_idx) %>%
    channelInfo(set_idx)
  
  df <- df[rowSums(!is.na(df[, grepl("^N_log2_R[0-9]{3}", names(df)), drop = FALSE])) > 0, ] %>%
    dplyr::arrange(!!rlang::sym(id), prot_acc) %>%
    dplyr::select(-grep("^R[0-9]{3}", names(.)))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "prot_hit_num", "prot_family_member", "prot_score", 
      "prot_matches", "prot_sequences", 
      "pep_var_mod", "pep_var_mod_pos", "pep_scan_title", 
      "pep_res_before", "pep_res_after", 
      "raw_file", "pep_query", "pep_summed_mod_pos", "pep_local_mod_pos")))
  
  df <- local({
    col_start <- which(names(df) == "Modifications") + 1
    col_end <- which(names(df) == "Charge") - 1
    
    if (!(is_empty(col_start) | is_empty(col_end))) {
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
      "ID", "Protein group IDs", 
      "Peptide ID", "Mod. peptide ID", "Evidence ID", 
      "Length"))) %>% 
    dplyr::select(-grep("site IDs$", names(.)))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "number", "modifications", 
      "variableSites", "nterm", "previous_aa", "sequence", "next_aa", 
      "cys", "searchCycle", "L/H", "accession_numbers", "entry_name", 
      "matched_parent_mass"))) 
  
  if ("pep_scan_title" %in% names(df)) {
    df <- df %>% 
      dplyr::mutate(pep_scan_title = gsub("\\\\", "~~", pep_scan_title)) %>%
      dplyr::mutate(pep_scan_title = gsub("^File.*~~", "", pep_scan_title))
  }
  
  # summarise log2FC and intensity from the same `set_idx` at one or multiple LCMS injections
  if (method_psm_pep == "mean") {
    df_num <- aggrNums(mean)(df, !!rlang::sym(id), na.rm = TRUE)
  } else if (method_psm_pep == "top.3") {
    df_num <- TMT_top_n(df, !!rlang::sym(id), na.rm = TRUE)
  } else if (method_psm_pep == "weighted.mean") {
    df_num <- TMT_wt_mean(df, !!rlang::sym(id), na.rm = TRUE)
  } else {
    df_num <- aggrNums(median)(df, !!rlang::sym(id), na.rm = TRUE)
  }
  
  df_first <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    med_summarise_keys(id)
  
  df <- list(df_first, df_num) %>%
    purrr::reduce(left_join, by = id)
  
  if (id == "pep_seq_mod") {
    df <- df %>% dplyr::select(-pep_seq)
  } else {
    df <- df %>% dplyr::select(-pep_seq_mod)
  }
  
  df <- cbind.data.frame(df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df)), drop = FALSE],
                         df[, grepl("I[0-9]{3}", names(df)), drop = FALSE], 
                         df[, grepl("log2_R[0-9]{3}", names(df)), drop = FALSE]) %>%
    dplyr::mutate_at(.vars = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)),
                     list(~ replace(.x, is.infinite(.x), NA)))
  
  df <- local({
    col_r <- grepl("^log2_R[0-9]{3}", names(df))
    cf <- apply(df[, col_r, drop = FALSE], 2, median, na.rm = TRUE)
    df <- cbind(df[, -grep("^N_log2_R[0-9]{3}", names(df))],
                sweep(df[, col_r], 2, cf, "-") %>%
                  `colnames<-`(paste("N", colnames(.), sep="_")))
    
    col_int <- grepl("^I[0-9]{3}", names(df))
    df  <- cbind(df[, -grep("^N_I[0-9]{3}", names(df))],
                 sweep(df[, col_int, drop=FALSE], 2, 2^cf, "/") %>%
                   `colnames<-`(paste("N", colnames(.), sep="_")))
    
    return(df)
  })
  
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


#'Interim peptide tables
#'
#'\code{PSM2Pep} summarises
#'\code{\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSMs}} to
#'peptides by individual TMT experiments and LC/MS series.
#'
#'In general, fields other than \code{log2FC} and \code{intensity} are
#'summarized with median statistics. One exception is with \code{pep_expect} in
#'Mascot or \code{PEP} in MaxQuant where geometric mean is applied.
#'
#'@param method_psm_pep Character string; the method to summarise the
#'  \code{log2FC} and the \code{intensity} of \code{PSMs} by peptide entries.
#'  The descriptive statistics includes \code{c("mean", "median", "top.3",
#'  "weighted.mean")} with \code{median} being the default. The
#'  \code{log10-intensity} of reporter ions at the \code{PSMs} levels will be
#'  the weight when summarising \code{log2FC} with \code{"top.3"} or
#'  \code{"weighted.mean"}.
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
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'  
#'@return Tables under \code{PSM} folder for each TMT experiment and LC/MS
#'  series: \code{TMTset1_LCMSinj1_PSM_N.txt}, \code{TMTset1_LCMSinj2_PSM_N.txt}...
#'
#'@example inst/extdata/examples/PSM2Pep_.R
#'@import stringr dplyr purrr rlang  magrittr
#'@export
PSM2Pep <- function (method_psm_pep = c("median", "mean", "weighted.mean", "top.3"), ...) {
  old_opt <- options(max.print = 99999, warn = 0)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 2000000, warn = 1)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("PSM2Pep"), add = TRUE)
  
  dots <- rlang::enexprs(...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)

  method_psm_pep <- rlang::enexpr(method_psm_pep)
  if (method_psm_pep == rlang::expr(c("median", "mean", "weighted.mean", "top.3"))) {
    method_psm_pep <- "median"
  } else {
    method_psm_pep <- rlang::as_string(method_psm_pep)
    stopifnot(method_psm_pep %in% c("median", "mean", "weighted.mean", "top.3"))
  }
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  dir.create(file.path(dat_dir, "Peptide\\cache"), recursive = TRUE, showWarnings = FALSE)
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
    reorder_files(n_TMT_sets(label_scheme_full))
  
  purrr::walk(as.list(filelist), ~ {
    fn_prx <- gsub("_PSM_N.txt", "", .x, fixed = TRUE)
    set_idx <- as.integer(gsub(".*TMTset(\\d+)_.*", "\\1", fn_prx))
    injn_idx <- as.integer(gsub(".*LCMSinj(\\d+).*", "\\1", fn_prx))
    fn_pep <- file.path(dat_dir, "Peptide", paste0("TMTset", set_idx, "_LCMSinj", injn_idx, "_Peptide_N.txt"))
    
    df <- read.csv(file.path(dat_dir, "PSM", .x), check.names = FALSE, header = TRUE,
                   sep = "\t", comment.char = "#") %>% 
      dplyr::select(-grep("^sd_log2_R", names(.))) %>% 
      calcPepide(label_scheme = label_scheme, id = !!id, method_psm_pep = method_psm_pep, 
                 group_pep_by = group_pep_by, set_idx = set_idx, injn_idx = injn_idx)
    
    df <- dplyr::bind_cols(
      df %>% dplyr::select(grep("^pep_", names(.))), 
      df %>% dplyr::select(-grep("^pep_", names(.))), 
    )
    
    df <- dplyr::bind_cols(
      df %>% dplyr::select(grep("^prot_", names(.))),
      df %>% dplyr::select(-grep("^prot_", names(.))),
    )
    
    write.table(df, file.path(dat_dir, "Peptide", paste0(fn_prx, "_Peptide_N.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE)		
  })
}


#' Add the `pep_seq_mod` field to MaxQuant PSMs
#'
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
add_maxquant_pepseqmod <- function(df, use_lowercase_aa) {

  my_tolower <- function(x, ch = "^") {
    locales <- gregexpr(ch, x) %>% .[[1]] %>% `+`(., 1)
    lowers <- map(locales, ~ substr(x, .x, .x)) %>% tolower()
    
    for (i in seq_along(lowers)) {
      substr(x, locales[i], locales[i]) <- lowers[i]
    }
    
    x <- gsub(ch, "", x)
    
    return(x)
  }
  
  if (!use_lowercase_aa) {
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = "."))
  } else {
    # (1) all non-terminal modifications: M(ox) -> m ...
    df <- df %>% 
      tidyr::separate("pep_seq_mod", c("nt", "interior", "ct"), sep = "_") %>% 
      dplyr::mutate(interior = gsub("([A-Z]){1}\\([^\\(\\)]*\\)", paste0("@", "\\1"), interior)) %>% 
      dplyr::mutate_at(vars("interior"), ~ map_chr(.x, my_tolower, "@")) %>% 
      tidyr::unite(pep_seq_mod, nt, interior, ct, sep = ".", remove = TRUE)

    # (2) phospho: pS -> s, pT -> t, pY -> y
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("pS", "s", pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = gsub("pT", "t", pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = gsub("pY", "y", pep_seq_mod)) 

    # (3-1) add "_" to sequences from protein N-terminal acetylation
    df <- local({
      n_ac <- df %>% dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      if (nrow(n_ac) > 0) {
        n_ac <- n_ac %>% dplyr::mutate(pep_seq_mod = gsub("^\\.\\(ac\\)", "_", pep_seq_mod))
        df <- rbind(rest, n_ac)
      }
      
      return(df)
    })
    
    # (3-2) add "_" to sequences from protein C-terminal amidation
    df <- local({
      c_am <- df %>% dplyr::filter(grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      if (nrow(c_am) > 0) {
        c_am <- c_am %>% dplyr::mutate(pep_seq_mod = gsub("\\.\\(am\\)$", "_", pep_seq_mod))
        df <- rbind(rest, c_am)
      }
      
      return(df)
    })
    
    # (4-1) "~" for "(Protein N-term)" other than acetylation
    df <- local({
      n_ac <- df %>% dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      other_n <- df %>% 
        dplyr::filter(grepl("Protein N-term", Modifications, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% dplyr::filter(!grepl("Protein N-term", Modifications, fixed = TRUE))
      
      if (nrow(other_n) > 0) {
        other_n <- other_n %>% dplyr::mutate(pep_seq_mod = gsub("^\\.\\([^\\(\\)]*\\)", "~", pep_seq_mod))
        df <- rbind(rest, n_ac, other_n)
      }
      
      return(df)
    })    
    
    # (4-2) "~" for "(Protein C-term)" other than amidation
    df <- local({
      c_am <- df %>% dplyr::filter(grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      other_c <- df %>% 
        dplyr::filter(grepl("Protein C-term", Modifications, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% dplyr::filter(!grepl("Protein C-term", Modifications, fixed = TRUE))
      
      if (nrow(other_c) > 0) {
        other_c <- other_c %>% dplyr::mutate(pep_seq_mod = gsub("\\.\\([^\\(\\)]*\\)$", "~", pep_seq_mod))
        df <- rbind(rest, c_am, other_c)
      }
      
      return(df)
    })
    
    # (5-1) "^" peptide "(N-term)" modification
    df <- local({
      nt <- df %>% dplyr::filter(grepl("(N-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("(N-term)", Modifications, fixed = TRUE))
      
      if (nrow(nt) > 0) {
        nt <- nt %>% 
          dplyr::mutate(pep_seq_mod = gsub("^\\.([_~]{0,1})\\([^\\(\\)]*\\)", paste0("\\1", "^"), pep_seq_mod)) 
        
        df <- rbind(rest, nt)
      }
      
      return(df)
    })
    
    # (5-2) "^" peptide "(C-term)" modification
    df <- local({
      ct <- df %>% dplyr::filter(grepl("(C-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("(C-term)", Modifications, fixed = TRUE))
      
      
      if (nrow(ct) > 0) {
        ct <- ct %>% 
          dplyr::mutate(pep_seq_mod = gsub("\\.\\([^\\(\\)]*\\)([_~]{0,1}$)", paste0("^", "\\1"), pep_seq_mod)) 
        
        df <- rbind(rest, ct)
      }
      
      return(df)
    })
    
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("\\.", "", pep_seq_mod))
    
    # (6) other N- or C-terminal modifications better but not named with "N-term" or "C-term": 
    #     (py)C -> c, (gl)Q -> q
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("(^[_~]{0,1})\\([^\\(\\)]*\\)", paste0("\\1", "^"), pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = gsub("\\([^\\(\\)]*\\)([_~]{0,1}$)", paste0("^", "\\1"), pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = ".")) %>% 
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = "."))
  }
  
  return(df)
}


#'Splits PSM tables
#'
#'\code{splitPSM_mq} splits the PSM outputs after \code{rmPSMHeaders()}. It
#'separates PSM data by TMT experiment and LC/MS injection.
#'
#'Different to \code{splitPSM} used in \code{Mascot} processes,
#'\code{pep_seq_mod}, \code{prot_n_psm}, \code{prot_n_pep} and \code{pep_n_psm}
#'are calculated later in \code{annotPSM_mq}. This is suitable mainly because
#'there is no columns like \code{prot_matches_sig} and \code{prot_sequences_sig}
#'need to be updated after data merging in \code{splitPSM_mq}.
#'
#'@param pep_unique_by A character string for annotating the uniqueness of
#'  peptides in \code{MaxQuant} PSMs. At the \code{group} default, the
#'  uniqueness of peptides is by protein groups. At a more strigent criterion of
#'  \code{protein}, the uniqueness of peptides is by protein entries. A new
#'  column of \code{pep_isunique} with corresponding logical TRUE or FALSE will
#'  be added to the PSM reports.
#'@param corrected_int A logical argument for uses with \code{MaxQuant} data. At
#'  the TRUE default, values under columns "Reporter intensity corrected..." in
#'  \code{MaxQuant} PSM results (\code{msms.txt}) will be used. Otherwise,
#'  "Reporter intensity" values without corrections will be used.
#'@param rm_reverses A logical argument for uses with \code{MaxQuant} data. At
#'  the TRUE default, \code{Reverse} entries will be removed.
#'@inheritParams splitPSM
#'@import dplyr tidyr
#'@importFrom stringr str_split
#'@importFrom magrittr %>%
splitPSM_mq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, 
                        pep_unique_by = "group", corrected_int = TRUE, rptr_intco = 1000, 
                        rm_craps = FALSE, rm_reverses = TRUE, annot_kinases = FALSE, 
                        plot_rptr_int = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$")
  
  if (rlang::is_empty(filelist)) {
    stop(paste("No PSM files were found under", file.path(dat_dir), 
               "\nCheck that the names of PSM files start with `msms`."), call. = FALSE)
  }

  df <- purrr::map(file.path(dat_dir, filelist), read.csv, 
                   check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::bind_rows() 
  
  # exception: empty string under `Proteins`
  df <- df %>% dplyr::filter(as.character(.$Proteins) > 0)

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `msms.txt` for `filter_` varargs.")
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::rename(ID = id)
  
  if (pep_unique_by == "group") {
    df <- df %>% 
      dplyr::mutate(pep_isunique = ifelse(grepl(";", `Protein group IDs`), 0, 1))
  } else if (pep_unique_by == "protein") {
    df <- df %>% 
      dplyr::mutate(pep_isunique = ifelse(grepl(";", Proteins), 0, 1))
  }
  
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
    df <- df %>% dplyr::filter(.$Reverse != "+")
  }
  
  df <- df %>% 
    `names_pos<-`(grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+$", names(.)), 
                  gsub("TMT-", "I", as.character(TMT_levels)))
  
  if (TMT_plex == 11) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_int <- NULL
  }  
  
  df <- local({
    phos_idx <- grep("Phospho (STY) Probabilities", names(df), fixed = TRUE)

    if (!is_empty(phos_idx)) {
      phos <- df %>% 
        dplyr::select(phos_idx) %>% 
        purrr::map(~ str_extract_all(.x, "\\([^()]+\\)")) %>% 
        .[[1]] %>% 
        purrr::map(~ substring(.x, 2, nchar(.x) - 1)) %>% 
        purrr::map(as.numeric) 
      
      phos_max <- phos %>% 
        purrr::map(max, na.rm = TRUE) %>% 
        unlist()

      phos_sort <- phos %>% purrr::map(sort, decreasing = TRUE)
      maxs <- phos_sort %>% purrr::map(`[`, 1) %>% unlist()
      sec_maxs <- phos_sort %>% purrr::map(`[`, 2) %>% unlist()
      
      df <- bind_cols(df, pep_phospho_locprob = phos_max, pep_phospho_locdiff = maxs - sec_maxs) %>% 
        dplyr::mutate(is_locprob_one = equals(1, pep_phospho_locprob)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locdiff"), ~ replace(.x, is_locprob_one, 1)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locprob"), ~ replace(.x, is.infinite(.x), NA)) %>% 
        dplyr::select(-is_locprob_one)
    }
    
    return(df)
  })
  
  if (TMT_plex > 0) {
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) %>% 
      dplyr::rename(
        pep_seq = Sequence, 
        prot_acc = Proteins, 
        RAW_File = `Raw file`, 
      ) 

    df <- df %>% 
      dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc))
  }
  
  acc_type <- parse_acc(df)
  stopifnot(length(acc_type) == 1)
  
  df <- df %>% annotPrn(fasta)
  if (annot_kinases) df <- df %>% annotKin(acc_type)
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  }
  
  df <- cbind.data.frame(
    df[, grepl("^[a-z]", names(df))], 
    df[, grepl("^[A-Z]", names(df)) & !grepl("^[IR][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[R][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[I][0-9]{3}[NC]*", names(df))]
  )
  
  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>%
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` under the experimental summary file."),
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

  for (i in seq_along(df_split)) {
    df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)
    
    out_fn <- fn_lookup %>%
      dplyr::filter(TMT_inj == names(df_split)[i]) %>%
      dplyr::select(filename) %>%
      unique() %>%
      unlist() %>%
      paste0(., ".csv")
    
    df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
    
    write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
    
    if (plot_rptr_int & TMT_plex > 0) {
      df_int <- df_split[[i]] %>% 
        .[, grepl("^I[0-9]{3}", names(.))]
      
      rptr_violin(df = df_int, filepath = file.path(dat_dir, "PSM\\rprt_int\\raw", gsub("\\.csv", "\\.png", out_fn)), 
                  width = 8, height = 8)
    }
  }
}


#' Annotates MaxQuant PSM results
#'
#' \code{annotPSM_mq} adds fields of annotation to MaxQuant PSM tables.
#'
#'@inheritParams load_expts
#'@inheritParams annotPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
annotPSM_mq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, expt_smry = "expt_smry.xlsx", 
                        rm_krts = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, use_lowercase_aa = TRUE, ...) {
	
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 5000000)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(., n_TMT_sets)
  
  for (set_idx in seq_len(n_TMT_sets)) {
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename = 
                           do.call('rbind', strsplit(as.character(sublist), '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    # LCMS injections under the same TMT experiment
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") %>%
        dplyr::rename(pep_seq_mod = `Modified sequence`) %>%
        dplyr::select(which(names(.) == "pep_seq_mod"),
                      which(names(.) != "pep_seq_mod"))
      
      acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
      species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()

      if (rm_krts) {
        df <- df %>% 
          dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
      }

      df <- df %>% add_maxquant_pepseqmod(use_lowercase_aa)
      
      df <- df %>% 
        dplyr::mutate(pep_start_discrepancy = ifelse(grepl("Protein N-term", Modifications) & (pep_start > 2), TRUE, FALSE)) %>% 
        dplyr::filter(!pep_start_discrepancy) %>% 
        dplyr::select(-pep_start_discrepancy)
      
      if (TMT_plex > 0) df <- mcPSM(df, set_idx)
			
			df <- df %>% 
			  calcSD_Splex(group_psm_by) %>% 
			  `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
			  dplyr::right_join(df, by = group_psm_by) %>% 
			  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
			  dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
			  dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
			  dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 4))
			
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

			df <- dplyr::bind_cols(
			  df %>% dplyr::select(grep("^pep_", names(.))), 
			  df %>% dplyr::select(-grep("^pep_", names(.))), 
			)
			
			df <- list(df, prot_n_psm, prot_n_pep) %>%
			  purrr::reduce(left_join, by = group_pep_by)
			
			df <- dplyr::bind_cols(
			  df %>% dplyr::select(grep("^prot_", names(.))), 
			  df %>% dplyr::select(-grep("^prot_", names(.))), 
			)
			
			df <- dplyr::bind_cols(
			  df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
			  df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
			  df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
			)
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int & TMT_plex > 0) {
        df_int <- df %>% .[, grepl("^N_I[0-9]{3}", names(.))]
        rptr_violin(df = df_int, 
                    filepath = file.path(dat_dir, "PSM\\rprt_int\\mc", 
                                         paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_rprt.png")), 
                    width = 8, height = 8)
      }
      
      if (plot_log2FC_cv & TMT_plex > 0) {
        sd_violin(df = df, id = !!group_psm_by, 
                  filepath = file.path(dat_dir, "PSM\\log2FC_cv\\raw", 
                                       paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_sd.png")), 
                  width = 8, height = 8, type = "log2_R", adjSD = FALSE, is_psm = TRUE)
      }
    }
    
  }
}



#' Splits PSM tables from \code{Spectrum Mill}
#'
#' \code{splitPSM_sm} splits the PSM outputs from \code{Spectrum Mill}. It
#' separates PSM data by TMT experiment and LC/MS injection.
#'
#' @inheritParams splitPSM
#' @inheritParams splitPSM_mq
#' @import dplyr tidyr readr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
splitPSM_sm <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, 
                        rm_craps = FALSE, rm_krts = FALSE, rptr_intco = 1000, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)
  
  if (is.null(fasta)) stop("FASTA file(s) not provided.", call. = FALSE)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$")
  
  if (rlang::is_empty(filelist)) 
    stop(paste("No PSM files were found under", file.path(dat_dir), 
               "\nCheck that the names of PSM files start with `PSMexport`."), call. = FALSE)
  
  df <- purrr::map(file.path(dat_dir, filelist), readr::read_delim, delim = ";") %>% 
    dplyr::bind_rows() 
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `PSMexport[...].ssv` for `filter_` varargs.")
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::rename(prot_acc = accession_number) %>% 
    annotPrn(fasta)  
  
  if (TMT_plex == 11) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_int <- NULL
  }  
  
  if (TMT_plex > 0) {
    df <- df %>% 
      `names<-`(gsub("^TMT_([0-9]{3}[NC]?)", "I\\1", names(.))) %>% 
      dplyr::select(-grep("^I[0-9]{3}[NC]?_[0-9]{3}[NC]?$", names(.))) %>% 
      as.data.frame()
    
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
    
    df <- df %>% 
      dplyr::rename(RAW_File = `filename`) %>% 
      dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", RAW_File)) %>% 
      dplyr::mutate(pep_seq = toupper(sequence)) %>% 
      dplyr::mutate(pep_miss = ifelse(.$next_aa == "(-)", str_count(pep_seq, "[KR]"), str_count(pep_seq, "[KR]") - 1))
  }
  
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)
  
  if (rm_craps) {
    data(package = "proteoQ", prn_annot_crap)
    
    craps <- prn_annot_crap %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::select(acc_type) %>% 
      unlist()
    
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) df <- annotKin(df, acc_type)
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  # M._sequence.c; -._sequence.c; n.sequence.c; -.sequence.c
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  }
  
  df <- df %>% 
    dplyr::mutate(pep_seq = toupper(pep_seq))
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(-grep("^prot_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
  )
  
  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>% 
      # dplyr::select(which(not_all_zero(.))) %>% # don't: empty channels are all NA too
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>% 
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` under the experimental summary file."),
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
  
  for (i in seq_along(df_split)) {
    df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)
    
    out_fn <- fn_lookup %>%
      dplyr::filter(TMT_inj == names(df_split)[i]) %>%
      dplyr::select(filename) %>%
      unique() %>%
      unlist() %>%
      paste0(., ".csv")
    
    df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
    
    write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
    
    if (plot_rptr_int & TMT_plex > 0) {
      df_int <- df_split[[i]] %>% 
        .[, grepl("^I[0-9]{3}", names(.))]
      
      rptr_violin(df = df_int, filepath = file.path(dat_dir, "PSM\\rprt_int\\raw", gsub("\\.csv", "\\.png", out_fn)), 
                  width = 8, height = 8)
    }
  }
}


#' Annotates \code{Spectrum Mill} PSMs following \code{splitPSM_sm} and
#' \code{cleanupPSM}
#'
#' \code{annotPSM_sm} adds fields of annotation to PSM tables.
#'
#' @inheritParams load_expts
#' @inheritParams annotPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
annotPSM_sm <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, expt_smry = "expt_smry.xlsx", 
                        rm_krts = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, use_lowercase_aa = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 5000000)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(., n_TMT_sets)
  
  for(set_idx in seq_len(n_TMT_sets)){
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename = 
                           do.call('rbind', strsplit(as.character(sublist),
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") %>%
        dplyr::mutate(pep_seq_mod = pep_seq) %>% 
        dplyr::mutate(pep_seq_mod = as.character(pep_seq_mod)) %>% 
        dplyr::select(which(names(.) == "pep_seq_mod"),
                      which(names(.) != "pep_seq_mod"))
      
      acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
      species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
      
      if (rm_krts) {
        df <- df %>% 
          dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
      }
      
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
        if (!is_empty(which(names(df) == "nterm"))) {
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
        if (!is_empty(which(names(df) == "cterm"))) {
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
          dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
          dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = "."))        
      }

      # median centering
      if (TMT_plex > 0) df <- mcPSM(df, set_idx)
      
      # add SD columns
      df <- df %>% 
        calcSD_Splex(group_psm_by) %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
        dplyr::right_join(df, by = group_psm_by) %>% 
        dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
        dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 4))
      
      # add column `pep_n_psm` et al.
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
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^pep_", names(.))), 
      )
      
      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))), 
        df %>% dplyr::select(-grep("^prot_", names(.))), 
      )
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
      )
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int & TMT_plex > 0) {
        df_int <- df %>% 
          .[, grepl("^N_I[0-9]{3}", names(.))]
        
        rptr_violin(df = df_int, 
                    filepath = file.path(dat_dir, "PSM\\rprt_int\\mc", 
                                         paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_rprt.png")), 
                    width = 8, height = 8)
      }
      
      if (plot_log2FC_cv & TMT_plex > 0) {
        sd_violin(df = df, id = !!group_psm_by, 
                  filepath = file.path(dat_dir, "PSM\\log2FC_cv\\raw", 
                                       paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_sd.png")), 
                  width = 8, height = 8, type = "log2_R", adjSD = FALSE, is_psm = TRUE)
      }
      
    }
    
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


