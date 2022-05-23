#' Extracts RAW MS file names
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
extract_raws <- function(raw_dir = NULL, dat_dir = NULL) 
{
  if (is.null(raw_dir)) 
    stop("\"raw_dir\" cannot be NULL.", call. = FALSE)

  if (is.null(dat_dir)) 
    dat_dir <- get_gl_dat_dir()
  
  fns <- dir(raw_dir, pattern = "\\.raw$", full.names = FALSE)

  data.frame(Index = seq_along(fns), RAW_File = fns) %T>% 
    write.table(file.path(dat_dir, "raw_list.txt"), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE)
		
	message("List of \".raw\" file names stored in ", 
	        file.path(dat_dir, "raw_list.txt"))
}


#' Finds the names of Mascot PSM files.
#' 
#' @param filelist A list of PSM files.
#' @inheritParams load_expts
find_mascot_psmraws <-function(filelist = NULL, dat_dir = NULL) 
{
  dir.create(file.path(dat_dir, "PSM/cache"), 
             recursive = TRUE, 
             showWarnings = FALSE)
  
  purrr::walk(filelist, batchPSMheader, 
              dat_dir, 
              TMT_plex = 16, # 18 later?
              proc_extdata = FALSE, 
              warns = FALSE)
  
  purrr::walk(gsub("\\.csv$", "_hdr_rm.csv", filelist), ~ {
    df <- read.delim(file.path(dat_dir, "PSM/cache", .x), 
                     sep = ',', 
                     check.names = FALSE, 
                     header = TRUE, 
                     stringsAsFactors = FALSE, 
                     quote = "\"", 
                     fill = TRUE , 
                     skip = 0)
    
    l <- df$pep_scan_title[1]
    
    if (grepl("File:~", l, fixed = TRUE)) { # MSConvert
      raws <- df[, "pep_scan_title"] %>% 
        gsub("^[^ ]+ File:\\~(.*)\\.(raw|d)\\~.*", "\\1", .) %>% 
        unique()
    } 
    else if (grepl("File: ~", l, fixed = TRUE)) { # PD
      raws <- df[, "pep_scan_title"] %>% 
        gsub("\\\\", "/", .) %>% 
        gsub('^.*/(.*)\\.(raw|d)\\~.*', '\\1', .) %>% 
        unique()
    }
    else if (grepl(", File:", l, fixed = TRUE)) { # custom from Bruker MGF, add "~" later...
      raws <- df[, "pep_scan_title"] %>% 
        gsub("\\\\", "/", .) %>% 
        gsub('^.*, File:(.*)\\.(raw|d)\\, .*', '\\1', .) %>% 
        unique()
    }
    else {
      stop("Mascot RAW_File name in \"pep_scan_title\" is not in the format of ", 
           "\"MSConvert\" or \"Proteome DisCoverer\".\n", 
           "Contact the developer about the new format.")
    }
    
    if (length(raws)) {
      out_nm <- paste0("msraws_in_", 
                       .x %>% 
                         gsub("\\.[^.]*$", "", .) %>% 
                         gsub("_hdr_rm$", "", .), 
                       ".txt")
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, out_nm), 
                    sep = "\t", 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE)
      
      message("MS file names stored in ", file.path(dat_dir, out_nm))
    }
  })
} 


#' Finds the names of MaxQuant PSM files.
#' 
#' @param filelist A list of PSM files.
#' @inheritParams load_expts
find_mq_psmraws <- function(filelist = NULL, dat_dir = NULL) 
{
  purrr::walk(filelist, ~ {
    df <- suppressWarnings(
      readr::read_tsv(file.path(dat_dir, .x), 
                      col_types = cols(`Raw file` = col_character())))
    
    raws <- unique(df$`Raw file`)
    
    if (length(raws)) {
      out_nm <- paste0("msraws_in_", gsub("\\.[^.]*$", "", .x), ".txt")
      
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, out_nm), 
                    sep = "\t", 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE)
      
      message("MS file names stored in ", file.path(dat_dir, out_nm))
    }
  })
}


#' Finds the names of Spectrum Mill PSM files.
#' 
#' @param filelist A list of PSM files.
#' @inheritParams load_expts
find_sm_psmraws <- function(filelist = NULL, dat_dir = NULL) 
{
  purrr::walk(filelist, ~ {
    df <- suppressWarnings(
      readr::read_delim(file.path(dat_dir, .x), 
                        delim = ";", 
                        col_types = cols(filename = col_character()))) %>% 
      dplyr::mutate(filename = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", filename))
    
    raws <- unique(df$filename)
    
    if (length(raws)) {
      out_nm <- paste0("msraws_in_", gsub("\\.[^.]*$", "", .x), ".txt")
      
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, out_nm), 
                    sep = "\t", 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE)
      
      message("MS file names stored in ", file.path(dat_dir, out_nm))
    }
  })
}


#' Finds the names of MSFragger PSM files.
#' 
#' @param filelist A list of PSM files.
#' @inheritParams load_expts
find_mf_psmraws <- function(filelist = NULL, dat_dir = NULL) 
{
  purrr::walk(filelist, ~ {
    
    df <- suppressWarnings(
      readr::read_delim(file.path(dat_dir, .x), 
                        delim = "\t", 
                        col_types = cols(Spectrum = col_character())) 
    ) %>% 
      dplyr::mutate(Spectrum = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", Spectrum))
    
    raws <- unique(df$Spectrum)
    
    if (length(raws)) {
      out_nm <- paste0("msraws_in_", gsub("\\.[^.]*$", "", .x), ".txt")
      
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, out_nm), 
                    sep = "\t", 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE)
      
      message("MS file names stored in ", file.path(dat_dir, out_nm))
    }
  })
}


#' Finds the names of proteoQ PSM files.
#' 
#' @param filelist A list of PSM files.
#' @inheritParams load_expts
find_pq_psmraws <- function(filelist = NULL, dat_dir = NULL) 
{
  purrr::walk(filelist, ~ {
    
    df <- suppressWarnings(
      readr::read_tsv(file.path(dat_dir, .x), 
                      col_types = cols(prot_acc = col_character())) 
    )
    
    raws <- unique(df$raw_file)
    
    if (length(raws)) {
      out_nm <- paste0("msraws_in_", gsub("\\.[^.]*$", "", .x), ".txt")
      
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, out_nm), 
                    sep = "\t", 
                    col.names = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE)
      
      message("MS file names stored in ", file.path(dat_dir, out_nm))
    }
  })
}


#' Extracts the names of RAW MS files from PSM data
#'
#' \code{extract_psm_raws} extracts the MS file names from the PSM data under
#' the working directory.
#'
#' @inheritParams load_expts
#' @examples
#' \donttest{
#' extract_psm_raws()
#' }
#'
#' @import dplyr tidyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @export
extract_psm_raws <- function(dat_dir = NULL) 
{
  if (is.null(dat_dir)) 
    dat_dir <- get_gl_dat_dir()
  else 
    assign("dat_dir", dat_dir, envir = .GlobalEnv)

  type <- find_search_engine(dat_dir)

  pattern <- switch(type, 
    mascot = "^F[0-9]+.*\\.csv$", 
    mq = "^msms.*\\.txt$",
    sm = "^PSMexport.*\\.ssv$", 
    mf = "^psm.*\\.tsv$", 
    pq = "^psm[QC]{1}.*\\.txt$", 
    stop("Data type needs to be one of \"mascot\", \"maxquant\",", 
         " \"spectrum_mill\", \"msfragger\" or \"proteoM\".", 
         call. = FALSE)
  )
  
  filelist <- list.files(path = file.path(dat_dir), pattern = pattern)
  
  if (!length(filelist)) 
    stop("No ", toupper(type), " files(s) under ", dat_dir, call. = FALSE)

  if (type == "pq") {
    fileQ <- filelist %>% 
      .[grepl("^psmQ.*\\.txt$", .)]

    if (length(fileQ)) 
      filelist <- fileQ
  }

  switch (type,
    mascot = find_mascot_psmraws(filelist, dat_dir),
    mq = find_mq_psmraws(filelist, dat_dir),
    sm = find_sm_psmraws(filelist, dat_dir),
    mf = find_mf_psmraws(filelist, dat_dir),
    pq = find_pq_psmraws(filelist, dat_dir), 
    stop("Unknown \"type\".", call. = FALSE)
  )
}


#' Sets names for Mascot ratio and intensity columns. 
#' 
#' Including spacer columns.
#' 
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
set_mascot_colnms <- function (TMT_plex = 10L) 
{
  cols_ratio <- find_ratio_cols(TMT_plex)
  cols_int <- find_int_cols(TMT_plex)
  
  dummy_ratios <- paste0("dummy_", cols_ratio)
  dummy_ints <- paste0("dummy_", cols_int)
  
  nms_ratio <- purrr::map2(dummy_ratios, cols_ratio, ~ {
    c(.x, .y)
  }) %>% 
    unlist()
  
  nms_int <- purrr::map2(dummy_ints, cols_int, ~ {
    c(.x, .y)
  }) %>% 
    unlist()
  
  nms <- paste(
    c(nms_ratio, nms_int),
    collapse = ','
  )
}


#' Processes Mascot PSMs
#' 
#' @param filename The name of a PSM file.
#' @param proc_extdata Logical; if TRUE, processes extended PSM outputs;
#'   otherwise, skips them. Setting the value to FALSE allows bypassing
#'   time-consuming steps.
#' @param warns Logical; if TRUE, show warning message(s).
#' @inheritParams load_expts
#' @inheritParams TMT_levels
batchPSMheader <- function(filename = NULL, dat_dir = NULL, TMT_plex = 10L, 
                           proc_extdata = TRUE, warns = TRUE) 
{
  output_prefix <- gsub("\\.csv$", "", filename)
  
  data_all <- readLines(file.path(dat_dir, filename))
  
  fixed_mod_row <- grep("Fixed modifications", data_all)
  var_mod_row <- grep("Variable modifications", data_all)
  pep_seq_rows <- grep("pep_seq", data_all)
  unassign_hits_row <- grep("Peptide matches not assigned to protein hits", data_all)
  
  err_msg <- c(
    "\nCheck at least the following items when exporting PSMs:\n", 
    # "[x] - required\n", 
    # "[s] - suggested\n", 
    # "[o] - optional\n", 
    
    "\n=== Search Information ===\n", 
    "[x] Modification deltas\n", 
    
    # "\n=== Protein Hit Information ===\n", 
    # "[x] Mass (Da)\n", 
    # "[x] Number of queries matched\n\n", 
    
    "\n=== Peptide Match Information ===\n", 
    "[x] Sequence\n", 
    "[x] Variable Modifications\n", 
    "[x] Query title\n", 
    "[x] Peptide quantitation (for TMT)\n\n", 
    
    "Check more applicable items to get more out of your data.\n"
  )
  
  local({
    keys <- list(`Fixed modifications` = fixed_mod_row, 
                 `Variable modifications` = var_mod_row)

    oks <- !purrr::map_lgl(list(fixed_mod_row, var_mod_row), 
                           purrr::is_empty)
    
    if (!all(oks)) {
      warning("\nNot all of \"", 
              paste(keys, collapse = ", "), 
              "\"\n are present in the header of ", filename, ".\n",
              err_msg,
              call. = FALSE)
    }
  })

  if (!length(pep_seq_rows)) {
    stop("Column \"pep_seq\" not found in ", filename, ".\n",
         "(1) \"Sequence\", \"Variable Modifications\", \"Query title\" ", 
         "are required during PSM exports.\n", 
         "(2) \"Peptide quantitation\" is further required for TMT quantitation.", 
         call. = FALSE)
  }
  
  local({
    keys <- c("pep_seq", "pep_var_mod_pos", "pep_scan_title", 
              "prot_acc", "pep_isbold")

    oks <- purrr::map_lgl(keys, ~ grepl(.x, data_all[pep_seq_rows])[1]) %>% 
      `names<-`(keys) 
    
    if (!all(oks)) {
      stop("\nNot all of \"", 
           paste(keys, collapse = ", "), "\"\n are present in ", 
           filename, ".\n", 
           err_msg, 
           call. = FALSE)
    }
  })

  if (length(pep_seq_rows) > 2L) {
    stop("Unhandled input formats: more than two \"pep_seq\" found in rows ", 
         pep_seq_rows, " in ", filename, ".", 
         call. = FALSE)
  }

  # https://stackoverflow.com/questions/18893390/splitting-on-comma-outside-quotes
  pattern <- ",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)"
  
  local({
    data_header <- data_all[1 : (pep_seq_rows[1] - 1)]
    data_header <- gsub("\"", "", data_header, fixed = TRUE)
    
    write.table(data_header, file.path(dat_dir, "PSM/cache",
                                       paste0(output_prefix, "_header", ".txt")),
                sep = "\t", col.names = FALSE, row.names = FALSE)
  })
  
  # --- "Query Level Information" ---
  if (length(pep_seq_rows) > 1L) {
    # seq(), comp(), tag(), etc.: 
    # "StringTitle","Scan number range","Retention time range","qualifiers",
    
    # [x] Query level Search Parameters (not yet used): 
    # "Peptide Mass Tolerance","Peptide Mass Tolerance Units",
    # "Variable modifications","Instrument type",
    
    # [x] MS/MS Peak lists: 
    # "StringIons1","StringIons2","StringIons3"
    
    if (proc_extdata) {
      message("Processing \"Query Level Information\": ", filename)
      
      queries <- stringr::str_split(data_all[pep_seq_rows[2]:length(data_all)], 
                                    pattern) %>% 
        list_to_dataframe() %>% 
        `names<-`(.[1, ]) %>% 
        `names<-`(gsub("\"", "", names(.))) %>% 
        `[`(-1, ) %>% 
        dplyr::select(-ncol(.)) 
      
      purrr::walk(c("query_number", "moverz", "charge", "intensity"), ~ {
        if (! .x %in% names(queries)) {
          stop("Column ", .x, " not found in ", filename, 
               call. = FALSE)
        }
      }, queries)
      
      purrr::walk(c("qualifiers", "charge", "pep_var_mod", "pep_scan_title", 
                    "Scan number range", "Retention time range", 
                    "StringTitle", "StringIons1", "StringIons2", 
                    "StringIons3"), ~ {
                      if (.x %in% names(queries)) {
                        queries[[.x]] <- gsub("\"", "", queries[[.x]])
                        queries[[.x]][stringr::str_length(queries[[.x]]) == 0] <- NA
                        queries <<- queries
                      }
                    }, queries)
      
      purrr::walk(c("query_number", "moverz", "intensity", "pep_index", 
                    "pep_rank", "pep_exp_mz", "pep_exp_mr", "pep_exp_z", 
                    "pep_calc_mr", "pep_delta", "pep_start", "pep_end", 
                    "pep_miss", "pep_score", "pep_homol", "pep_ident", 
                    "pep_expect", 
                    "Scan number range", 
                    "Retention time range", "TotalIonsIntensity", 
                    "NumVals"), ~ {
                      if (.x %in% names(queries)) {
                        queries[[.x]] <<- as.numeric(queries[[.x]])
                      }
                    }, queries)
      
      purrr::walk(c("pep_isbold", "pep_isunique"), ~ {
        if (.x %in% names(queries)) {
          queries[[.x]] <<- as.logical(queries[[.x]])
        }
      }, queries)
      
      # note that pep_var_mod_conf: 81.72% ...
      if ("pep_var_mod_conf" %in% names(queries)) {
        queries <- queries %>% 
          dplyr::mutate(pep_var_mod_conf = 
                          as.numeric(sub("%", "", pep_var_mod_conf))) %>% 
          dplyr::mutate(pep_var_mod_conf = pep_var_mod_conf/100)
      } 
      else {
        queries <- queries %>% 
          dplyr::mutate(pep_var_mod_conf = NA)
      }
      
      # assume the same intensity for the same pep_query at different pep_seq
      # pep_query	moverz	intensity	pep_expect	pep_seq	pep_var_mod_pos	pep_var_mod_conf
      # 20107	669.37537	251692.1094	0.0058	FNNVEAGK	0.00200000.0	50.00%
      # 20107	NA	NA	0.0058	FNNVEAGK	0.02000000.0	50.00%
      # 20107	NA	NA	0.45	GFGDVGEAK		
      # 20107	NA	NA	0.45	EAFQEQK		
      
      ##  --- down-fill values --- 
      # Don't tidyr::fill before dplyr::group as 
      # TRUE NA (not measurable) may occur to "intensity"
      # 
      # rows <- queries %>% 
      #   dplyr::filter(!is.na(moverz), is.na(intensity)) %>% 
      #   .[["query_number"]] %>% 
      #   unique()
      
      queries <- queries %>% 
        dplyr::arrange(query_number, intensity) %>% 
        dplyr::group_by(query_number)
      
      # same query, same c("moverz", "charge", "intensity")
      purrr::walk(c("moverz", "charge", "intensity"), ~ {
        if (.x %in% names(queries)) {
          queries <<- queries %>% tidyr::fill(!!rlang::sym(.x)) 
        }
      }, queries)
      
      queries <- queries %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(intensity = round(intensity, digits = 0))
      
      queries %>% 
        saveRDS(file.path(dat_dir, "PSM/cache", 
                          paste0(gsub("\\.[^.]*$", "", filename), "_queryInfo.rds")))
      
      rm(list = c("queries"))
    }
    
    data_all <- data_all[1:(pep_seq_rows[2] - 4)]
  }

  if (length(unassign_hits_row)) {
    if (proc_extdata) {
      message("Processing \"Unassigned queries\": ", filename)
      
      stringr::str_split(data_all[(unassign_hits_row+2):length(data_all)], 
                         pattern) %>% 
        list_to_dataframe() %>% 
        dplyr::select(-ncol(.)) %>% 
        saveRDS(file.path(dat_dir, "PSM/cache", 
                          paste0(gsub("\\.[^.]*$", "", filename), "_unassigned.rds")))
    }

    data_all <- data_all[1:(unassign_hits_row - 2)]
  } 
  
  # use "this_plex", not TMT_plex by metadata
  # (not yet channel or column padding)
  this_plex <- find_mascot_tmtplex(data_all)
  
  if (this_plex != TMT_plex && warns) {
    warning("Mascot PSMs suggest a TMT ", this_plex, "-plex, ", 
            "which is different to the ", TMT_plex, "-plex in \"expt_smry.xlsx\".", 
            call. = FALSE)
  }
  
  # !!! KEEP quotation marks !!!
  # Do NOT gsub("\"", "", data_psm)
  # e.g., pep_scan_title may contain "Last name, First name" in the file path.

  data_psm <- data_all[pep_seq_rows[1]:length(data_all)]
  data_psm <- gsub("\"---\"", "\"-1\"", data_psm, fixed = TRUE)
  data_psm <- gsub("\"###\"", "\"-1\"", data_psm, fixed = TRUE)
  
  # --- for simulated data
  data_psm <- gsub("---", "\"-1\"", data_psm, fixed = TRUE)
  data_psm <- gsub("###", "\"-1\"", data_psm, fixed = TRUE)
  
  if (this_plex) {
    data_psm[1] <- paste0(data_psm[1], ",", set_mascot_colnms(this_plex))
  } 
  
  # [x] emPAI
  # [x] Protein quantitation
  # -> partial columns to the right
  # ..."emPAI",0.80,"Quantitation summary for protein","127N/126","1.139"...
  
  data_psm <- local({
    # exclude "Protein quantitation"
    for(i in seq_along(data_psm)) {
      data_psm[i] <- gsub(",\"Quantitation summary for protein\".*", 
                          "", data_psm[i])
    }
    
    empai <- data_psm %>% .[grepl("\"emPAI\"", .)] 
    
    if (!length(empai)) 
      return(data_psm)
    
    empai <- c(data_psm[1], empai) %>% 
      stringr::str_split(pattern) %>% 
      list_to_dataframe() %>% 
      `names<-`(.[1, ]) %>% 
      `names<-`(gsub("\"", "", names(.), fixed = TRUE)) %>% 
      `[`(-1, ) 
    
    cols <- c(which(names(empai) == "prot_acc"), 
              which(empai[1, ] == "\"emPAI\"") + 1)

    stopifnot(length(cols) == 2L)
    
    empai <- empai[, cols] %>% 
      `names<-`(c("prot_acc", "prot_empai")) %>% 
      dplyr::mutate(prot_acc = gsub("[1-9]+::", "", prot_acc)) %>% 
      dplyr::mutate(prot_acc = gsub("\"", "", prot_acc, fixed = TRUE), 
                    prot_empai = as.numeric(prot_empai)) %T>% 
      saveRDS(file.path(dat_dir, "PSM/cache", 
                        paste0(gsub("\\.[^.]*$", "", filename), "_empai.rds")))
                    
    for(i in seq_along(data_psm)) {
      data_psm[i] <- gsub(",\"emPAI\".*", "", data_psm[i])
    }
    
    invisible(data_psm)
  })

  writeLines(data_psm, 
             file.path(dat_dir, "PSM/cache", paste0(output_prefix, "_hdr_rm.csv")))
}


#' Removes PSM headers
#'
#' \code{rmPSMHeaders} removes the header of PSM from
#' \href{https://http://www.matrixscience.com/}{Mascot} outputs. It also
#' removes the spacer columns in the fields of ratio and intensity values.
#' 
#' @inheritParams normPSM
#' @return Intermediate PSM table(s).
#'
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
rmPSMHeaders <- function(parallel = TRUE) 
{
  old_opts <- options(warn = 0)
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        message("Remove PSM headers --- Completed.")
      }
    }, 
    add = TRUE
  )
  
  dat_dir <- get_gl_dat_dir()
  
  filelist <- list.files(path = file.path(dat_dir), 
                         pattern = "^F[0-9]+.*\\.csv$")
  
  if (!length(filelist)) {
    stop("No PSM files(s) with \".csv\" extension under ", dat_dir, 
         call. = FALSE)
  }
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  n_files <- length(filelist)
  
  if (parallel && (n_files > 1)) {
    n_cores <- min(parallel::detectCores(), n_files)
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterMap(cl, batchPSMheader, filelist, 
                         MoreArgs = list(dat_dir, TMT_plex))
    
    parallel::stopCluster(cl)
  } 
  else {
    purrr::walk(filelist, batchPSMheader, dat_dir, TMT_plex)
  }
  
  .saveCall <- TRUE
  
  invisible(NULL)
}


#' Adds column "pep_var_mod_conf"
#'
#' Currently only for Mascot.
#'
#' Not yet supports Mascot timsTOF format and thus \code{pep_scan_range} is
#' \code{NA_integer_} not \code{NA_character_}.
#'
#' @param df A data frame of PSMs
#' @inheritParams load_expts
add_mod_conf <- function(df = NULL, dat_dir = NULL) 
{
  dat_file <- unique(df$dat_file)
  
  if (length(dat_file) != 1L)
    stop("The number of unique \"dat_file\"s is not one.")

  file <- file.path(dat_dir, "PSM/cache", paste0(dat_file, "_queryInfo.rds"))
  
  if (!file.exists(file)) {
    df <- df %>% 
      dplyr::mutate(pep_tot_int = NA_real_, 
                    pep_scan_range = NA_integer_, 
                    pep_ret_range = NA_integer_, 
                    pep_ms2_sumint = NA_real_,
                    pep_n_ions = NA_integer_, 
                    pep_locprob = NA_real_,
                    pep_locdiff = NA_real_)
    
    return(df)
  }
  
  # (1) no need of "raw_file" being part of "uniq_id" as "pep_query" is 
  #     by each "dat_file" and irrespective to "raw_file"
  # 
  # (nevertheless, the same `pep_query` may have multiple `pep_seq`s and 
  #  each `pep_seq` may have multiple `pep_var_mod_pos`s)
  # 
  # (2) no need of `dat_file` being part of the `uniq_by` 
  #   as `df` has been split by `dat_file`
  # 
  # (3) when left_join(df, queries) by uniq_id = c("pep_query", "pep_seq"), 
  # the same `uniq_id` but different "pep_var_mod_pos" forms different rows; 
  # different rows under the same `uniq_id` are used for `df_first` and `df_second`.

  # Don't: 
  # uniq_by <- c("pep_query", "pep_seq", "pep_var_mod_pos")

  uniq_by <- c("pep_query", "pep_seq")
  
  queries <- readRDS(file) %>% 
    dplyr::rename(pep_query = query_number) %>% 
    tidyr::unite(uniq_id, uniq_by, sep = ".", remove = FALSE) 
  
  purrr::walk2(c("intensity", "Scan number range", 
                 "Retention time range", 
                 "TotalIonsIntensity", "NumVals", 
                 "StringIons1", "StringIons2", "StringIons3"), 
               c("pep_tot_int", "pep_scan_range", 
                 "pep_ret_range", "pep_ms2_sumint", 
                 "pep_n_ions", 
                 "pep_ions_first", "pep_ions_second", "pep_ions_third"), 
               ~ {
                 if (.x %in% names(queries)) 
                   queries <<- queries %>% 
                     dplyr::rename(!!.y := !!.x)
               }, queries)
  
  ## === CPTAC phosphopeptides, fractions F5 and F10 ===
  # 
  # Senario 1: one `uniq_id`, multiple `pep_var_mod_conf`
  # 
  # In `df`:
  # uniq_id           pep_var_mod_pos
  # 48940.SGEGEVSGLMR 0.50000000000.0
  # 
  # In `queries`: 
  # uniq_id           pep_tot_int pep_var_mod_pos pep_var_mod_conf
  # 48940.SGEGEVSGLMR    21739660 0.50000000000.0           1.00  
  # 48940.SGEGEVSGLMR    21739660 0.00000050000.0           0.0002
  # 
  # `df2` after left_join by `uniq_id`: 
  # uniq_id           pep_var_mod  pep_var_mod_pos pep_var_mod_conf
  # 48940.SGEGEVSGLMR Phospho (ST) 0.50000000000.0           1.00
  # 48940.SGEGEVSGLMR Phospho (ST) 0.50000000000.0           0.0002 <-- 
  # 
  # the same `pep_var_mod_pos` but different `pep_var_mod_conf`; 
  # namely, the second `pep_var_mod_pos` etc. do not link to `pep_var_mod_conf`.
  # to avoid clumsiness, one solution is to order data by decreasing 
  #   `pep_var_mod_conf` and only keeps the fist and shows the delta.
  # this are what `df_first` and `df_second` for.
  # 
  # 
  # Senario 2: one `uniq_id`, multiple `pep_var_mod_conf` and `prot_acc`
  # 
  # "183371.VFPGSTTEDYNLIVIER"
  # 
  # `df`:
  # prot_acc                  uniq_id       pep_var_mod_pos
  # TIF1B_HUMAN 183371.VFPGSTTEDYNLIVIER 0.00000500000000000.0
  # TIF1B_MOUSE 183371.VFPGSTTEDYNLIVIER 0.00000500000000000.0
  # 
  # `queries` (`prot_acc` not available)
  # uniq_id                  pep_tot_int pep_var_mod_pos       pep_var_mod_conf
  # 183371.VFPGSTTEDYNLIVIER     1714404 0.00000500000000000.0           0.493 
  # 183371.VFPGSTTEDYNLIVIER     1714404 0.00005000000000000.0           0.493 
  # 183371.VFPGSTTEDYNLIVIER     1714404 0.00000050000000000.0           0.0131
  # 183371.VFPGSTTEDYNLIVIER     1714404 0.00000000060000000.0           0     
  # 
  # `df2` after left_join by `uniq_id`
  # prot_acc    uniq_id                  pep_var_mod  pep_var_mod_pos       pep_var_mod_conf
  # TIF1B_HUMAN 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0.493 
  # TIF1B_HUMAN 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0.493 
  # TIF1B_MOUSE 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0.493 
  # TIF1B_MOUSE 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0.493 
  # TIF1B_HUMAN 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0.0131
  # TIF1B_MOUSE 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0.0131
  # TIF1B_HUMAN 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0     
  # TIF1B_MOUSE 183371.VFPGSTTEDYNLIVIER Phospho (ST) 0.00000500000000000.0           0     
  # 
  # Again, `pep_var_mod_pos` values are not exactly correct but OK if only keep the first 
  # in each group of `uniq_id`. 
  # Next calculates the delta; however, rows corresponding to `TIF1B_MOUSE` will be gone. 
  # TIF1B_MOUSE` is considered double-dipping and be later removed when applying parsimony. 
  # However, for better self-containedness, 
  #  join by uniq_id2 <- c("prot_acc", uniq_by) instead of uniq_id. 
  # 
  # Note that queries doesn't have column `prot_acc` and a two-stage joining 
  # with `uniq_id` followed by `uniq_id2` is required.
  

  # two-stage joining and subsetting
  # (1) df + queries by uniq_id
  # (2) grouped by uniq_id2
  uniq_by2 <- c("prot_acc", uniq_by)
  
  df2 <- df %>% 
    tidyr::unite(uniq_id, uniq_by, sep = ".", remove = FALSE) %>% 
    dplyr::left_join(queries %>% 
                       # be explicit on intended columns
                       dplyr::select(which(names(.) %in% 
                                             c("uniq_id", 
                                               "pep_tot_int", 
                                               "pep_var_mod_conf", 
                                               "pep_scan_range", 
                                               "pep_ret_range", 
                                               "pep_ms2_sumint", 
                                               "pep_n_ions", 
                                               "pep_ions_first", 
                                               "pep_ions_second", 
                                               "pep_ions_third"))), 
                     by = "uniq_id") %>% 
    dplyr::mutate(.n = row_number()) %>% 
    dplyr::arrange(-pep_var_mod_conf) %>% 
    tidyr::unite(uniq_id2, uniq_by2, sep = ".", remove = FALSE) %>% 
    dplyr::group_by(uniq_id2)
  
  df_first <- df2 %>% 
    dplyr::filter(row_number() == 1L) %>% 
    dplyr::ungroup(uniq_id2) %>% 
    dplyr::arrange(.n) 
  
  df_second <- df2 %>% 
    dplyr::filter(row_number() == 2L) %>% 
    dplyr::ungroup(uniq_id2) %>% 
    dplyr::arrange(.n) %>% 
    dplyr::select(uniq_id2, pep_var_mod_conf) %>%
    dplyr::rename(pep_var_mod_conf_2 = pep_var_mod_conf)
  
  ## possible `pep_var_mod_conf` values are all NA; 
  #  thus `pep_locprob` & `pep_locdiff` both be NA 
  # query_number                  pep_seq              pep_var_mod_pos pep_var_mod_conf
  # 1       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400000000000000050002000.0             <NA>
  # 2       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400000000000000500002000.0             <NA>
  # 3       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400000000005000000002000.0             <NA>
  # 4       312209 MDPLSELQDDLTLDDTSQALNQLK 1.400050000000000000002000.0             <NA>
  
  df_first %>% 
    dplyr::left_join(df_second, by = "uniq_id2") %>% 
    dplyr::mutate(pep_locdiff = pep_var_mod_conf - pep_var_mod_conf_2) %>% 
    dplyr::select(-pep_var_mod_conf_2) %>% 
    dplyr::rename(pep_locprob = pep_var_mod_conf) %>% 
    dplyr::arrange(.n) %>% 
    dplyr::select(-c("uniq_id", "uniq_id2", ".n"))
}


#' Adds the \code{pep_seq_mod} field to Mascot PSMs
#' 
#' @inheritParams locate_outliers
#' @inheritParams splitPSM
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_mascot_pepseqmod <- function(df = NULL, use_lowercase_aa = TRUE, 
                                 purge_phosphodata = TRUE) 
{
  dat_dir <- get_gl_dat_dir()
  
  if (!"dat_file" %in% names(df))
    stop("Column \"dat_file\" not found.")

  dat_id <- unique(df$dat_file)
  dat_file <- file.path(dat_dir, "PSM/cache", paste0(dat_id, "_header.txt"))
  
  if (length(dat_id) != 1L)
    stop("The number of \"dat_file\" is not one.")
  
  if (!file.exists(dat_file))
    stop(dat_file, " not found.")

  df_header <- readLines(dat_file)

  fixed_mods <- local({
    zero_out <- tibble::tibble(Mascot_abbr = character(),
                               Description = character(), 
                               Delta_mass = character()) 

    spacer <- "\"\""
    rows_spacer <- which(df_header == spacer)
    
    row_fixmod <- grep("Fixed modifications", df_header)[1]
    
    if (is.na(row_fixmod)) 
      return(zero_out)
    
    row_fixmod_next <- row_fixmod + 1L
    
    if (df_header[row_fixmod_next] != spacer) 
      return(zero_out)

    row_next_spacer <- rows_spacer[which(rows_spacer > row_fixmod_next)][1]

    df_header[(row_fixmod_next + 2L) : (row_next_spacer - 1L)] %>%
      gsub("\"", "", ., fixed = TRUE) %>%
      data.frame() %>%
      tidyr::separate(".", sep = ",", c("Mascot_abbr", "Description", "Delta_mass"))
  })
  
  var_mods <- local({
    zero_out <- tibble::tibble(Mascot_abbr = character(),
                               Description = character(), 
                               Delta_mass = character(), 
                               Filename = character()) 
    
    spacer <- "\"\""
    rows_spacer <- which(df_header == spacer)
    
    row_varmod <- grep("Variable modifications", df_header)[1]
    
    if (is.na(row_varmod)) 
      return(zero_out)
    
    row_varmod_next <- row_varmod + 1L
    
    if (df_header[row_varmod_next] != spacer) 
      return(zero_out)
    
    row_next_spacer <- rows_spacer[which(rows_spacer > row_varmod_next)][1]

    df_header[(row_varmod_next + 2) : (row_next_spacer - 1)] %>%
      gsub("\"", "", ., fixed = TRUE) %>%
      data.frame() %>%
      tidyr::separate(".", sep = ",", extra = "drop", 
                      c("Mascot_abbr", "Description", "Delta_mass")) %>%
      dplyr::mutate(Filename = gsub("[\\\\\\/\\:\\*\\?\\'\\<\\>\\|]", ".", Description))
  })
  
  if (is.null(df$pep_seq)) {
    stop("Column \"pep_seq\" not found.", call. = FALSE)
  }

  df <- local({
    is_phospho_expt <- any(grepl("Phospho\\s{1}\\(", var_mods$Description))
    
    if (is_phospho_expt && purge_phosphodata) {
      if ("pep_var_mod" %in% names(df)) {
        df <- df %>% dplyr::filter(grepl("Phospho\\s{1}\\(", pep_var_mod))
      } 
      else {
        warning("Missing column \"pep_var_mod\" for non-phosphopeptide removals.", 
                call. = FALSE)
      }
    }
    
    df
  })
  
  if (!nrow(var_mods)) {
    df$pep_seq_mod <- df$pep_seq
    
    warning("No variable modifications available in the header of \"", 
            dat_id, ".csv\".", 
            call. = FALSE)
            
    return(df)
  } 

  if (!use_lowercase_aa) {
    df <- df %>%
      dplyr::mutate(pep_seq_mod = paste0(pep_seq, "[", pep_var_mod_pos, "]"))
  } 
  else {
    df$pep_seq_mod <- df$pep_seq
    
    # (1) non terminal modifications
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(!grepl("N-term", Description, fixed = TRUE), 
                      !grepl("C-term", Description, fixed = TRUE)) 
      
      if (nrow(mod_tbl)) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub)) {
            # "-2" for the two characters, "0." ..., in 'pep_var_mod_pos'
            pos_matrix  <- gregexpr(mod, df_sub$pep_var_mod_pos) %>% 
              list_to_dataframe() %>% 
              purrr::map(~ {.x - 2}) %>%
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

      df
    })
    
    # (2-1) add "_" to sequences from protein N-terminal acetylation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Acetyl (Protein N-term)", Description, fixed = TRUE))
      
      nrow <- nrow(mod_tbl)
      stopifnot(nrow <= 1L)
      
      if (nrow == 1L) {
        mod <- mod_tbl$Mascot_abbr[1]
        
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)

        df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
        df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
        
        if (nrow(df_sub)) {
          df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0("_", pep_seq_mod))
          df <- rbind(df_rest, df_sub)
        }        
      }

      df
    })
    
    # (2-2) add "_" to sequences from protein C-terminal amidation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Amidated (Protein C-term)", Description, fixed = TRUE))

      nrow <- nrow(mod_tbl)
      stopifnot(nrow <= 1L)
      
      if (nrow == 1L) {
        mod <- mod_tbl$Mascot_abbr[1]
        
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
        
        df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
        df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
        
        if (nrow(df_sub)) {
          df_sub <- df_sub %>% 
            dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "_"))
          df <- rbind(df_rest, df_sub)
        }        
      }

      df
    })
    
    # (3-1) "~" for "(Protein N-term)" other than acetylation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Protein N-term", Description, fixed = TRUE)) %>%
        dplyr::filter(!grepl("Acetyl (Protein N-term)", Description, fixed = TRUE))
      
      if (nrow(mod_tbl)) {
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub)) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = paste0("~", pep_seq_mod))
            
            df <- rbind(df_rest, df_sub)
          }
        } 
      }

      df
    })
    
    # (3-2) "~" for "(Protein C-term)" other than amidation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Protein C-term", Description, fixed = TRUE)) %>%
        dplyr::filter(!grepl("Amidated (Protein C-term)", Description, fixed = TRUE))
      
      if (nrow(mod_tbl)) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub)) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "~"))
            
            df <- rbind(df_rest, df_sub)
          }
        } 
      }

      df
    })
    
    # (4-1) "^" for peptide "(N-term)"  
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("N-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Protein N-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl)) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub)) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = 
                              gsub("(^[_~]{0,1})(.)", 
                                   paste0("\\1", "^", "\\2"), 
                                   pep_seq_mod)) 
            df <- rbind(df_rest, df_sub)
          }
        } 
      }

      df
    })
    
    # (4-2) "^" peptide "(C-term)" 
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("C-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Protein C-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl)) {
        var_mods <<- var_mods %>% 
          dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
          
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% 
            dplyr::filter(grepl(mod, pep_var_mod_pos))
          
          df_rest <- df %>% 
            dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub)) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = gsub("(.)([_~]{0,1}$)", 
                                               paste0("\\1", "^", "\\2"), 
                                               pep_seq_mod)) 
            
            df <- rbind(df_rest, df_sub)
          }
        } 
      }

      df
    })
  }
  
  invisible(df)
}


#' Adds column "prot_empai"
#' 
#' @inheritParams add_mod_conf
add_empai <- function(df = NULL, dat_dir = NULL) 
{
  dat_file <- unique(df$dat_file)
  
  if (length(dat_file) != 1L)
    stop("The number of unique \"dat_file\"s is not one.")
  
  file <- file.path(dat_dir, "PSM/cache", paste0(dat_file, "_empai.rds"))
  
  if (!file.exists(file)) {
    return(df %>% dplyr::mutate(prot_empai = NA))
  }
  
  empai <- readRDS(file) %>% dplyr::filter(!duplicated(prot_acc))
  
  df %>% 
    dplyr::mutate(.prot_acc = gsub("[1-9]+::", "", prot_acc)) %>% 
    dplyr::left_join(empai, by = c(".prot_acc" = "prot_acc")) %>% 
    dplyr::select(-.prot_acc)
}


#' Adds columns \code{pep_n_psm}, \code{prot_n_psm} and \code{prot_n_pep}.
#'
#' @param df A data frame.
#' @param uniq_by A vector of column keys in \code{df} defining the levels of
#'   uniqueness in PSM entries.
#' @inheritParams normPSM
add_quality_cols <- function(df = NULL, group_psm_by = "pep_seq", 
                             group_pep_by = "prot_acc", uniq_by = NULL) 
{
  # PSMs might be duplicated: one PSM -> multiple rows (by NLs or pep_seq_mod)
  if (!is.null(uniq_by)) 
    df <- unique(df, by = uniq_by)
  
  group_psm_by <- rlang::enexpr(group_psm_by)
  group_pep_by <- rlang::enexpr(group_pep_by)
  
  pep_n_psm <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by)) %>%
    dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
    dplyr::summarise(pep_n_psm = n())
  
  prot_n_psm <- df %>%
    dplyr::select(!!rlang::sym(group_pep_by)) %>%
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_n_psm = n())
  
  prot_n_pep <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    unique() %>% 
    ## primary and same-set proteins (group_pep_by) can share the same peptides
    ## (group_psm_by). the primary or same-set proteins will be removed with the
    ## following filter. Values of prot_n_pep will be NA for the removed proteins 
    ## when joining back with df.
    # dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_n_pep = n())
  
  df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)
  
  df <- list(df, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = group_pep_by)
}


#' Parallel splits of PSM tables
#' 
#' @param df A list of data frames.
#' @param nm Names of the \code{df}'s.
#' @param fn_lookup A lookup of file names.
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams TMT_levels
psm_msplit <- function(df = NULL, nm = NULL, fn_lookup = NULL, dat_dir = NULL, 
                       plot_rptr_int = TRUE, TMT_plex = 10L) 
{
  df <- df %>% dplyr::select(-TMT_inj)
  
  out_fn <- fn_lookup %>%
    dplyr::filter(TMT_inj == nm) %>%
    dplyr::select(filename) %>%
    unique() %>%
    unlist() %>%
    paste0(".csv")
  
  df <- df %>% 
    dplyr::rename(raw_file = RAW_File) %T>% 
    readr::write_csv(file.path(dat_dir, "PSM/cache", out_fn))

  if (plot_rptr_int) {
    try(
      df_int <- df %>% 
        dplyr::select(grep("^I[0-9]{3}", names(.))) %>% 
        rptr_violin(filepath = file.path(dat_dir, "PSM/rprt_int/raw", 
                                         gsub("\\.csv", "\\.png", out_fn)), 
                    width = 8, height = 8)
    )
  }
  
  invisible(out_fn)
}


#' Finds peptides that are shared by proteins.
#'
#' Only for PSM tables where the same \code{pep_seq}  duplicated in rows with
#' different \code{prot_acc}'s . This for now means Mascot with same-set and
#' sub-set proteins enabled. Otherwise, e.g. when applying to other search
#' engines without such redundancy, every peptide would turn out to be unique.
#'
#' Output rows are unique at the levels of \code{c("prot_acc", "pep_seq")}. By
#' default, sequences of \code{Protein} terminals will be treated separately.
#'
#' @param df PSM data.
#' @param pep_id the column key of peptide sequences.
#' @param prot_id The column key of protein accessions.
find_shared_prots <- function(df = NULL, pep_id = "pep_seq", prot_id = "prot_acc") 
{
  # At pep_id == "pep_seq", 
  #   not yet to differentiate interior versus Protein terminal sequences.
  # E.g. "Acetyl (Protein N-term)" versus "Acetyl (N-term)":
  #   where the Protein qualifier in the former is provided by users.
  # The same MS query can, however, often give matches under either specification.
  # In other words, the "Protein" qualifier can remain arbitrary.
  
  local({
    nms <- names(df)
    
    lapply(c(pep_id, prot_id), function (x) {
      if (! x %in% nms) stop("Column \"", x, "\" not found.")
    })
  })
  
  if (pep_id == "pep_seq" && "pep_seq_mod" %in% names(df)) {
    df <- df %>% 
      dplyr::mutate(row_id. = row_number(), 
                    pep_isprotterm = ifelse(grepl("[_~]", pep_seq_mod), TRUE, FALSE))

    pt <- df %>% 
      dplyr::filter(pep_isprotterm) %>% 
      hfind_shared_prots(pep_id, prot_id) %>% 
      dplyr::mutate(pep_isprotterm = 1L)
    
    npt <- df %>% 
      dplyr::filter(!pep_isprotterm) %>% 
      hfind_shared_prots(pep_id, prot_id) %>% 
      dplyr::mutate(pep_isprotterm = 0L)
    
    cols <- c("pep_seq_term", 
              paste0("shared_", prot_id, "s"), 
              paste0("pep_n_", prot_id, "s"))

    x <- dplyr::bind_rows(pt, npt) %>% 
      dplyr::mutate(pep_seq_term = paste(pep_seq, pep_isprotterm, sep = ".")) %>% 
      dplyr::select(which(names(.) %in% cols)) %>% 
      dplyr::filter(!duplicated(pep_seq_term))
    
    rm(list = c("pt", "npt"))
    
    df <- df %>% 
      dplyr::mutate(pep_isprotterm = as.integer(pep_isprotterm)) %>% 
      dplyr::mutate(pep_seq_term = paste(pep_seq, pep_isprotterm, sep = ".")) %>% 
      dplyr::select(-pep_isprotterm)
    
    # "pep_seq_term" needed for Mascot
    ans <- df %>% 
      dplyr::left_join(x, by = "pep_seq_term") %>% 
      dplyr::arrange(row_id.) %>% 
      dplyr::select(-row_id.)
  }
  # should not incur
  else {
    cols <- c("pep_seq", 
              paste0("shared_", prot_id, "s"), 
              paste0("pep_n_", prot_id, "s"))
    
    x <- df %>% 
      hfind_shared_prots(pep_id, prot_id) %>% 
      dplyr::select(which(names(.) %in% cols)) %>% 
      dplyr::filter(!duplicated(pep_seq))
    
    ans <- df %>% 
      dplyr::left_join(x, by = "pep_seq")
  }
  
  invisible(ans)
}


#' Helper of find_shared_prots.
#' 
#' @inheritParams find_shared_prots
hfind_shared_prots <- function(df, pep_id = "pep_seq", prot_id = "prot_acc") 
{
  uniq_by <- c(pep_id, prot_id)
  rows <- !duplicated(df[, uniq_by])
  
  entries <- df[rows, uniq_by] %>% 
    dplyr::mutate(!!prot_id := gsub("[1-9]+::", "", !!rlang::sym(prot_id))) %>% 
    dplyr::group_by(!!rlang::sym(pep_id))
  
  # Keep the redundancy at `pep_id: 
  # the same peptides under different proteins -> different rows
  # pep_x  prot_a  a, b
  # pep_x  prot_b  b, a
  
  prot_maps <- entries %>% 
    dplyr::summarise(!!paste0("shared_", prot_id, "s") := 
                       paste(!!rlang::sym(prot_id), collapse  = ", "))
  
  pep_n_prot <- entries %>% 
    dplyr::summarise(!!paste0("pep_n_", prot_id, "s") :=  n())

  list(entries, prot_maps, pep_n_prot) %>% 
    purrr::reduce(dplyr::left_join, by = pep_id) %>% 
    dplyr::ungroup()
}


#' Find peptides that are shared by genes.
#' 
#' For MaxQuant and MSFragger.
#' 
#' @param df PSM data.
#' @param key the column key of mapped proteins.
#' @param sep Character string; the separator of the mapped proteins.
#' @inheritParams splitPSM
add_shared_genes <- function (df = NULL, key = "Proteins", sep = ";", 
                                 fasta = NULL, entrez = NULL) 
{
  stopifnot(all(c("pep_seq", key) %in% names(df)))

  sep_new = ", "
  
  if (sep != sep_new) {
    df <- df %>% 
      dplyr::mutate(shared_prot_accs = gsub(sep, sep_new, !!rlang::sym(key)))
  }

  # don't use preexisted acc_lookuap.rda
  # as many `mapped proteins` in `Proteins` not in the df[["prot_acc"]].
  
  # instead compile the new list of prot_accs from `Proteins`,
  # for new `acc_lookup`

  pep_seqs <- df[["shared_prot_accs"]] %>% 
    stringr::str_split(sep_new) %>% 
    list_to_dataframe() %>% 
    dplyr::bind_cols(
      df %>% dplyr::select(pep_seq), 
    ) %>% 
    # dplyr::mutate_at(.vars = which(names(.) != "pep_seq"), 
    #                  ~ gsub("\\.[0-9]*$", "", .x)) %>% 
    tidyr::gather("id.", "prot_acc", -pep_seq) %>% 
    dplyr::select(-id.) %>% 
    dplyr::filter(!is.na(prot_acc)) %>% 
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(pep_prot.)) %>% 
    dplyr::select(-pep_prot.)

  acc_lookup <- pep_seqs %>% 
    dplyr::filter(!duplicated(prot_acc)) %>% 
    parse_fasta(fasta, entrez, warns = FALSE) 

  out <- dplyr::left_join(pep_seqs, 
                   acc_lookup[, c("prot_acc", "gene")], 
                   by = "prot_acc") %>% 
    dplyr::select(-prot_acc) %>% 
    dplyr::filter(!is.na(gene)) %>% 
    tidyr::unite(pep_gene., pep_seq, gene, sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(pep_gene.)) %>% 
    dplyr::select(-pep_gene.) %>% 
    dplyr::group_by(pep_seq) %>% 
    dplyr::summarise(!!paste0("shared_", "gene", "s") := 
                       paste(gene, collapse  = ", ")) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(df, ., by = "pep_seq")
}


#' Finds the shared protein accessions when using MSFragger.
#'
#' Protein accessions will be used if the mapped accessions are present under
#' \code{prot_acc}; otherwise, \code{fasta_name} will be used. This allows quick
#' distinction whether mapped proteins are same-set/sub-set to the primary
#' entries or not.
#' 
#' @param df PSM data.
add_shared_prot_accs_mf <- function (df) 
{
  stopifnot(all(c("pep_seq", "prot_acc", "Mapped Proteins") %in% names(df)))

  ok <- tryCatch(load(file = file.path(dat_dir, "acc_lookup.rda")),
                 error = function(e) "e")
  
  if (ok != "acc_lookup") {
    stop("\"acc_lookup.rda\" not found under ", dat_dir, ".", 
         call. = FALSE)
  }
  
  # no need of `Mapped Proteins` to be part of tidyr::unite
  # (razors are those in fasta but in prot_acc; 
  # namely, `Mapped Proteins` consider all possible fasta entries)
  
  df_acc <- df %>% 
    dplyr::select(c("pep_seq", "prot_acc", "Mapped Proteins")) %>% 
    tidyr::unite("id.", c("pep_seq", "prot_acc"), sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(id.)) %>% 
    dplyr::select(-id.)
  
  # may be later remove "^rev_" and "^REV_" entries...
  
  df_accs <- df_acc[["Mapped Proteins"]] %>% 
    stringr::str_split(", ") %>% 
    list_to_dataframe()

  df_accs2 <- purrr::imap(df_accs, ~ {
    data.frame(fasta_name = .x) %>% 
      dplyr::left_join(acc_lookup[, c("prot_acc", "fasta_name")], 
                       by = "fasta_name") %>% 
      dplyr::mutate(prot_acc. = ifelse(is.na(prot_acc), fasta_name, prot_acc)) %>% 
      dplyr::select(prot_acc.) %>% 
      dplyr::rename(!!.y := prot_acc.)
  }) %>% 
    dplyr::bind_cols(
      df_acc %>% dplyr::select(pep_seq, prot_acc)
    )
  
  df_shared <- df_accs2 %>% 
    tidyr::gather("id.", "maps.", -prot_acc, -pep_seq) %>% 
    dplyr::select(-id.) %>% 
    tidyr::unite(uniq_id, c("pep_seq", "prot_acc", "maps."), sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(uniq_id)) %>% 
    tidyr::separate(uniq_id, into = c("pep_seq", "prot_acc", "maps."), 
                    sep = "@", remove = TRUE) %>% 
    dplyr::arrange(pep_seq, prot_acc, maps., na.last = FALSE) %>% 
    dplyr::mutate(maps. = replace(maps., maps. == "NA", NA_character_), 
                  maps. = ifelse(is.na(maps.), prot_acc, maps.)) %>% 
    dplyr::group_by(pep_seq) %>%
    dplyr::summarise(prot_acc = dplyr::first(prot_acc), 
                     shared_prot_accs = paste(maps., collapse = ", ")) %>% 
    tidyr::unite(pep_prot, c("pep_seq", "prot_acc"), sep = "@") %>% 
    dplyr::left_join(df %>% 
                       tidyr::unite(pep_prot, c("pep_seq", "prot_acc"), 
                                    sep = "@", 
                                    remove = FALSE), 
                     ., 
                     by = "pep_prot") %>% 
    dplyr::select(-pep_prot)
}


#' Find peptides that are shared by genes.
#'
#' For Spectrum Mill with duplicates, such as P06748|Q61937|P06748|Q61937,
#' removals under \code{accession_numbers}.
#'
#' @param df PSM data.
#' @param key the column key of mapped proteins.
#' @param sep Character string; the separator of the mapped proteins.
#' @inheritParams splitPSM
add_shared_sm_genes <- function (df = NULL, key = "Proteins", sep = ";", 
                                 fasta = NULL, entrez = NULL) 
{
  stopifnot(all(c("pep_seq", key) %in% names(df)))
  
  sep_new = ", "
  
  if (sep != sep_new) {
    df <- df %>% 
      dplyr::mutate(shared_prot_accs = gsub(sep, sep_new, !!rlang::sym(key)))
  }
  
  # don't use preexisted acc_lookuap.rda
  # as many `mapped proteins` in `Proteins` not in the df[["prot_acc"]].
  
  # instead compile the new list of prot_accs from `Proteins`,
  # for new `acc_lookup`
  
  pep_seqs <- df[["shared_prot_accs"]] %>% 
    stringr::str_split(sep_new) %>% 
    list_to_dataframe() %>% 
    dplyr::bind_cols(
      df %>% dplyr::select(pep_seq), 
    ) %>% 
    tidyr::gather("id.", "prot_acc", -pep_seq) %>% 
    dplyr::select(-id.) %>% 
    dplyr::filter(!is.na(prot_acc)) %>% 
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(pep_prot.)) %>% 
    dplyr::select(-pep_prot.)
  
  # need to update `shared_prot_accs` in `df` 
  # (`accession_numbers` redundancy: P06748, Q61937, P06748, Q61937)
  
  df <- pep_seqs %>% 
    dplyr::group_by(pep_seq) %>% 
    dplyr::summarise(shared_prot_accs. = paste(prot_acc, collapse  = ", ")) %>% 
    dplyr::left_join(df, ., by = "pep_seq") %>% 
    dplyr::mutate(shared_prot_accs = shared_prot_accs.) %>% 
    dplyr::select(-shared_prot_accs.)

  acc_lookup <- pep_seqs %>% 
    dplyr::filter(!duplicated(prot_acc)) %>% 
    parse_fasta(fasta, entrez, warns = FALSE) 
  
  out <- dplyr::left_join(pep_seqs, 
                          acc_lookup[, c("prot_acc", "gene")], 
                          by = "prot_acc") %>% 
    dplyr::select(-prot_acc) %>% 
    dplyr::filter(!is.na(gene)) %>% 
    tidyr::unite(pep_gene., pep_seq, gene, sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(pep_gene.)) %>% 
    dplyr::select(-pep_gene.) %>% 
    dplyr::group_by(pep_seq) %>% 
    dplyr::summarise(!!paste0("shared_", "gene", "s") := 
                       paste(gene, collapse  = ", ")) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(df, ., by = "pep_seq")
}


#' Processes PSMs.
#'
#' Common routines in PSM processing across search engines. 
#'
#' @param df PSM data.
#' @inheritParams normPSM
procPSMs <- function (df = NULL, scale_rptr_int = FALSE, 
                      rptr_intco = 0, rptr_intrange = c(0, 100), 
                      rm_craps = FALSE, rm_krts = FALSE, rm_allna = FALSE, 
                      annot_kinases = FALSE, 
                      plot_rptr_int = TRUE, parallel = TRUE) 
{
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  
  if (TMT_plex && scale_rptr_int) {
    df <- local({
      pattern <- "^I[0-9]{3}[NC]{0,1}"
      
      df <- df %>% 
        dplyr::mutate_at(.vars = grep(pattern, names(.)), 
                         ~ ifelse(.x == -1L, NA_real_, .x)) %>% 
        dplyr::mutate(.sumint = rowSums(.[grepl(pattern, names(.))], na.rm = TRUE)) %>% 
        dplyr::mutate(.pep_tot_int = ifelse(is.na(pep_tot_int), 
                                            .sumint, pep_tot_int)) %>% 
        dplyr::mutate_at(.vars = grep(pattern, names(.)), 
                         ~ .x * .pep_tot_int / .sumint) %>% 
        dplyr::select(-.sumint, -.pep_tot_int) 
    })
  } 
  
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
      split(.$acc_type, drop = TRUE) %>% 
      .[acc_types] %>% 
      purrr::imap( ~ annotKin(.x, .y)) %>% 
      do.call(rbind, .)
  } 
  else {
    df <- df %>% 
      dplyr::mutate(kin_attr = FALSE, kin_class = NA, kin_order = NA)
  }
  
  df <- df %>% order_mascot_psm_cols()
  
  # (1) Shared peptides will be removed in this step if 
  #   checked 'Unique peptide only' during Mascot PSM export and rm_allNA = TRUE;
  # (2) `> 1L` to bypass non-quantitative or LFQ workflows;
  # (3) MaxQuant Bruker no `Precursor Intensity` being filled;
  
  df <- df %>%
    dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), 
                     as.numeric) %>%
    dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), 
                     ~ ifelse(.x <= rptr_intco, NA_real_, .x)) %>% 
    dplyr::mutate_at(vars(grep("^I[0-9]{3}[NC]{0,1}", names(.))), ~ {
      x <- .x
      qts <- quantile(x, probs = rptr_intrange/100, na.rm = TRUE)
      x <- ifelse((x < qts[1]) | (x > qts[2]), NA_real_, .x)
    }) %>% 
    dplyr::arrange(RAW_File, pep_seq, prot_acc) 
  
  if (length(grep("^I[0-9]{3}", names(df))) > 1L && rm_allna) {
    df <- df %>% 
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}[NC]{0,1}", names(.))])) > 0L, 
                    rowSums(!is.na(.[grep("^I[0-9]{3}[NC]{0,1}", names(.))])) > 0L)
  }
  
  if (!nrow(df)) {
    stop("No non-trivial reporter-ion intensities/ratios available.\n", 
         "In the case of Mascot ", 
         "Have you forgot to check \"Peptide quantitation\" when exporting PSMs?\n",
         "See also https://proteoq.netlify.app/post/exporting-mascot-psms/.", 
         call. = FALSE)
  }
  
  # compare experimental and user-provided PSM files
  res <- check_raws(df)
  df <- res$df
  tmtinj_raw_map <- res$lookup
  rm(list = c("res"))

  # split by TMT and LCMS
  if (! "PSM_File" %in% names(tmtinj_raw_map)) {
    df_split <- df %>%
      dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>% 
      dplyr::group_by(TMT_inj) %>%
      dplyr::mutate(psm_index = row_number()) %>%
      data.frame(check.names = FALSE) %>% 
      split(.$TMT_inj, drop = TRUE)
  } 
  else {
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
      split(.$TMT_inj, drop = TRUE)
  }
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  
  if (length(missing_tmtinj)) {
    cat("\n")
    warning("TMT sets and LC/MS injections not have corresponindg PSM files:\n", 
            call. = FALSE)
    message(paste0("  TMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop("Remove mismatched \"TMT_Set\" and/or \"LCMS_Injection\" ", 
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
  
  n_files <- length(df_split)
  
  if (parallel && (n_files > 1L)) {
    n_cores <- min(parallel::detectCores(), n_files)
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    suppressWarnings(
      silent_out <- parallel::clusterMap(
        cl, psm_msplit, df_split, names(df_split), 
        MoreArgs = list(fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
      )
    )
    
    parallel::stopCluster(cl)
  } 
  else {
    purrr::walk2(df_split, names(df_split), psm_msplit, 
                 fn_lookup, dat_dir, plot_rptr_int, TMT_plex)
  }
}


#' Finds mismatches in RAW file names
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
  
  # automated frac_smry.xlsx may be based on wrong information 
  # from expt_smry.xlsx (e.g. Sample_ID removals)
  local({
    ls_raws <-  unique(label_scheme_full$RAW_File)
    fs_raws <-  unique(fraction_scheme$RAW_File)
    
    if (!(all(is.na(ls_raws)) || all(ls_raws %in% fs_raws))) {
      expt_smry <- match_call_arg(load_expts, expt_smry)
      frac_smry <- match_call_arg(load_expts, frac_smry)
      unlink(file.path(dat_dir, frac_smry))
      prep_fraction_scheme(dat_dir, expt_smry, frac_smry)

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
      stop("The combination of TMT_Set & LCMS_Injection ", 
           paste(extra_fs_tmt, collapse = ", "), 
           " in fraction scheme not found in label scheme.", 
           call. = FALSE)
    
    if (length(extra_ls_tmt)) 
      stop("The combination of TMT_Set & LCMS ", 
           paste(extra_fs_tmt, collapse = ", "), 
           " in label scheme not found in fraction scheme.", 
           call. = FALSE)
  })
  
  tmtinj_raw <- fraction_scheme %>% 
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>%
    dplyr::select(-Fraction) %>%
    dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
    dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File)) # Bruker
  
  ms_raws <- unique(df$RAW_File)
  
  # (sometime due to parsing error, i.e., of unknown formats)
  if (all(is.na(ms_raws)))
    stop("All values of \"RAW_File\" are NA.")
  
  label_scheme_raws <- unique(tmtinj_raw$RAW_File)
  
  # ms_raws <- df$RAW_File %>% unique()
  # label_scheme_raws <- tmtinj_raw$RAW_File %>% unique()
  
  missing_ms_raws <- ms_raws %>% .[! . %in% label_scheme_raws]
  wrong_label_scheme_raws <- label_scheme_raws %>% .[! . %in% ms_raws]
  
  # --- rm `missing_ms_raws` from `df`
  if (length(missing_ms_raws)) {
    df <- df %>% dplyr::filter(! RAW_File %in% missing_ms_raws)
    
    warning("RAW file (names) not in metadata and ", 
            "corresponding entries removed from PSM data:\n", 
            paste(missing_ms_raws, collapse = "\n"), 
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
#'@param pep_unique_by A character string for annotating the uniqueness of
#'  peptides. At the \code{group} default, the uniqueness of peptides is by
#'  groups with the collapses of same-set or sub-set proteins. At a more
#'  stringent criterion of \code{protein}, the uniqueness of peptides is by
#'  protein entries without grouping. On the other extreme of choice
#'  \code{none}, all peptides are treated as unique. A new column of
#'  \code{pep_isunique} with corresponding logical TRUE or FALSE will be added
#'  to the PSM reports. Note that the choice of \code{none} is only for
#'  convenience, as the same may be achieved by setting \code{use_unique_pep =
#'  FALSE} in \link{Pep2Prn}.
#' @param scale_rptr_int Logical; if TRUE, scales (up) MS2 reporter-ion
#'  intensities by MS1 precursor intensity: \eqn{I_{MS1}*(I_{x}/\sum I_{MS2})}.
#'  \eqn{I_{MS1}}, MS1 precursor intensity; \eqn{I_{MS2}}, MS2 reporter-ion
#'  intensity; \eqn{I_{x}}, MS2 reporter-ion intensity under TMT channel
#'  \eqn{x}. Note that the scaling will not affect \code{log2FC}.
#' @param rm_craps Logical; if TRUE,
#'   \href{https://www.thegpm.org/crap/}{cRAP} proteins will be removed.
#'   The default is FALSE.
#' @param rm_krts Logical; if TRUE, keratin entries will be removed. The default
#'   is FALSE.
#' @param annot_kinases Logical; if TRUE, proteins of human or mouse origins
#'   will be annotated with their kinase attributes. The default is FALSE.
#' @param rptr_intco Numeric; the threshold of reporter-ion intensity (TMT:
#'   \code{I126} etc.; LFQ: \code{I000}) being considered non-trivial. The
#'   default is 0 without cut-offs. The data nullification will not be applied
#'   synchronously to the precursor intensity (\code{pep_tot_int}) under the
#'   same PSM query. To guard against odds such as higher MS2 reporter-ion
#'   intensities than their contributing MS1 precursor intensity, employs for
#'   example \code{filter_... = rlang::exprs(pep_tot_int >= my_ms1_cutoff)}
#'   during \link{PSM2Pep}. The rule of thumb is that \code{pep_tot_int} is a
#'   single column; thus the corresponding data filtration against it may be
#'   readily achieved without introducing new arguments. By contrast,
#'   \code{rptr_intco} applies to a set of columns, \code{I126} etc.; it might
#'   be slightly more involved/laborious when applying suitable statements of
#'   \code{filter_} varargs.
#' @param rptr_intrange Numeric vector at length two. The argument specifies the
#'   range of reporter-ion intensities (TMT: \code{I126} etc.; LFQ: \code{I000})
#'   being considered non-trivial. The default is between 0 and 100 percentile
#'   without cut-offs. While argument \code{rptr_intco} employs a universal
#'   cut-off across samples by absolute values, \code{range_int} provides an
#'   alternative means of sample-specific thresholding of intensities by
#'   percentiles. The data nullification will not be applied synchronously to
#'   the precursor intensity under the same PSM query.
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
#' @inheritParams normPSM
#' @import dplyr tidyr stringr
#' @importFrom magrittr %>% %T>% %$% %<>% 
splitPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                     fasta = NULL, entrez = NULL, pep_unique_by = "group", 
                     scale_rptr_int = FALSE, 
                     rm_craps = FALSE, rm_krts = FALSE, rm_allna = FALSE, 
                     purge_phosphodata = TRUE, 
                     annot_kinases = FALSE, plot_rptr_int = TRUE, 
                     rptr_intco = 0, rptr_intrange = c(0, 100), 
                     use_lowercase_aa = TRUE, parallel = TRUE, ...) 
{
  # --- Outlines ---
  # (1.1) row filtration, column padding and psm file combinations
  # (Mascot: not yet available I000 or I126 etc., 
  # as ms1 intensities need to be imported from query data)
  # 
  # (1.2) remove spacer columns
  # 
  # (1.3) clean up prot_acc (1::, 2::)
  # convenience crap removals where the fasta names ended with "|"
  # (this can affect the uniqueness of peptides: 
  # by excluding the assessment of peptide sharedness with cRAPs)
  # 
  # (1.4) add query data & pep_seq_mod
  # dependency: pep_seq_mod for `group_psm_by = ...`
  # (precursor intensity and empai added if available)
  # (`pep_seq` changed from such as MENGQSTAAK to K.MENGQSTAAK.L)
  # 
  # 
  # (2.1) annotate proteins
  # dependency: `gene` for prot_n_pep etc.
  # (NA genes are replaced with accessions)
  # (not yet prot_cover and prot_icover as they require pep_start and pep_end)
  # 
  # (2.2) compile "preferred" columns
  # 
  # (2.3) find the shared prot_accs and genes for each peptide
  # 
  # (2.4) remove non-essential proteins after shared_prot_accs, shared_genes
  # 
  # (2.5) find the levels of uniqueness for peptides
  # uniqueness currently by prot_acc, not gene
  # 
  # (2.6) add peptide properties and prot_cover, prot_icover
  # dependency: prot_cover and prot_icover need pep_start and pep_end
  # 
  # (2.7) apply parsimony
  # (a) chimeric: 
  #     the same `pep_query` can be assigned to different `pep_seq` at different `pep_rank`
  # (b) positional difference:
  #     the same combination of `pep_query`, `pep_seq` and `pep_var_mod_pos` 
  #     can be assigned to different `prot_acc` 
  # 
  # (c) remove redundant peptides under each `dat_file`
  #     dbl-dipping peptides across `dat_file` will be handled in `mergePep`
  # 
  # (2.8) add columns pep_n_psm, prot_n_psm, prot_n_pep
  # dependency: after the removals of redundant PSMs, peptides
  # 
  # (2.9) other fields
  # for compatibility, updates after data merging etc.
  # --- End ---
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        message("Split PSMs by TMT experiments and LCMS series --- Completed.")
      }
    }, 
    add = TRUE
  )
  
  dat_dir <- get_gl_dat_dir()
	load(file = file.path(dat_dir, "label_scheme_full.rda"))
	load(file = file.path(dat_dir, "fraction_scheme.rda"))

	TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist = list.files(path = file.path(dat_dir, "PSM/cache"),
                        pattern = "^F[0-9]+.*_hdr_rm.csv$")

	if (!length(filelist)) {
	  stop(paste("Missing intermediate PSM file(s) under: ", 
	             file.path(dat_dir, "PSM/cache")), 
	       call. = FALSE)
	}
  
  # (1.1) row filtration, column padding and psm file combinations
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in \"F[...].csv\" for \"filter_\" varargs.")
  
  df <- purrr::map(filelist, pad_mascot_channels, !!!filter_dots) %>% 
    pad_mascot_fields() %>% 
    do.call(rbind, .)

  stopifnot(all(c("RAW_File", "dat_file") %in% names(df)))

  # (1.2) remove spacer columns
  df <- local({
    r_start <- grep("^dummy_", names(df))[1]

    if (is.na(r_start)) {
      # placeholder for the "same symmetry" as TMT
      df <- df %>% dplyr::mutate(R000 = 1, I000 = NA_real_) 
    } 
    else {
      int_end <- max(grep("^I[0-9]{3}[NC]{0,1}", names(df)))
      if (int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
    }

    df
  })
  
  # (1.3) clean up prot_acc
  df <- df %>% 
    dplyr::filter(prot_acc != "") %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]+::", "", prot_acc)) %>% 
    { if (rm_craps) dplyr::filter(., !grepl("\\|.*\\|$", prot_acc)) else . } 
  
  # (1.4) add query data & pep_seq_mod
  df <- df %>% 
    split(.$dat_file, drop = TRUE) %>% 
    purrr::map(add_mod_conf, dat_dir) %>% 
    purrr::map(add_empai, dat_dir) %>% 
    purrr::map(add_mascot_pepseqmod, use_lowercase_aa, purge_phosphodata) %>% 
    dplyr::bind_rows() 
  
  stopifnot("pep_tot_int" %in% names(df), 
            "pep_seq_mod" %in% names(df))
  
  if (!TMT_plex) {
    df <- df %>% 
      dplyr::mutate(I000 = pep_tot_int) %>% 
      dplyr::mutate(R000 = I000/I000, 
                    R000 = ifelse(is.infinite(R000), NA, R000)) 
  }
  
  stopifnot(any(grepl("^I[0-9]{3}[NC]{0,1}", names(df))))
  
  # (1.4.2) phospho
  if (all(c("pep_var_mod", "pep_locprob", "pep_locdiff") %in% names(df))) {
    df <- df %>% 
      dplyr::mutate(pep_phospho_locprob = 
                      ifelse(grepl("Phospho", pep_var_mod), pep_locprob, NA_real_)) %>% 
      dplyr::mutate(pep_phospho_locdiff = 
                      ifelse(grepl("Phospho", pep_var_mod), pep_locdiff, NA_real_))
  }
  
  # (2.1) annotate proteins
  df <- df %>% 
    annotPrn(fasta, entrez) %>%  
    { if (!"gene" %in% names(.)) dplyr::mutate(., gene = prot_acc) else . } %>% 
    dplyr::mutate(gene = ifelse(is.na(gene), prot_acc, gene))

  # (2.2) compile "preferred" columns
  df <- local({
    if (! "prot_hit_num" %in% names(df)) {
      df <- df %>% 
        add_entry_ids("prot_acc", "prot_hit_num") %>% 
        dplyr::arrange(prot_hit_num) %>% 
        reloc_col_before_first("prot_hit_num")
    } 
    else {
      df <- df %>% 
        dplyr::arrange(prot_hit_num)
    }
    
    if (! "prot_family_member" %in% names(df)) {
      warning("PSMs were exported from Mascot without protein grouping.", 
              call. = FALSE)
      
      df <- df %>% 
        dplyr::mutate(prot_family_member = 0L) %>% 
        reloc_col_after("prot_family_member", "prot_hit_num")
    }
    
    if (! "pep_isbold" %in% names(df)) {
      df$pep_isbold <- TRUE
    } 
    else {
      df <- df %>% 
        dplyr::mutate(pep_isbold = ifelse(pep_isbold == 1, TRUE, FALSE))
    }
    
    if (! "pep_isunique" %in% names(df)) {
      df$pep_isunique <- TRUE
    } 
    else {
      df <- df %>% 
        dplyr::mutate(pep_isunique = ifelse(pep_isunique == 1, TRUE, FALSE))
    }
    
    if ("pep_index" %in% names(df)) df$pep_index <- NULL
    if ("prot_index" %in% names(df)) df$prot_index <- NULL
    
    df <- df %>% 
      add_entry_ids("pep_seq", "pep_index") %>% 
      add_entry_ids("prot_acc", "prot_index")

    invisible(df)
  })
  
  # (2.3) finds the shared prot_accs and genes for each peptide
  df <- df %>% 
    find_shared_prots("pep_seq", "prot_acc") %>% 
    find_shared_prots("pep_seq", "gene")
  
  # (2.4) remove non-essential proteins after shared_prot_accs, shared_genes
  # run this asap; otherwise some primary peptides may be removed by other filters
  if ("prot_family_member" %in% names(df)) 
    df <- dplyr::filter(df, !is.na(prot_family_member))
  
  # (2.5) find the uniqueness of peptides
  df <- local({
    key <- if ("pep_seq_term" %in% names(df)) "pep_seq_term" else "pep_seq"
    
    cols <- c(key, 
              paste0("shared_", group_pep_by, "s"), 
              paste0("pep_n_", group_pep_by, "s"))
    
    shared_peps <- df %>% 
      dplyr::select(which(names(.) %in% cols)) %>% 
      dplyr::filter(!!rlang::sym(paste0("pep_n_", group_pep_by, "s")) >= 2L) %>% 
      .[[key]] %>% 
      unique() %>% 
      sort()
    
    df <- df %>% 
      dplyr::mutate(pep_literal_unique = 
                      ifelse(!!rlang::sym(key) %in% shared_peps, FALSE, TRUE), 
                    pep_razor_unique = pep_isunique) %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
    
    if (pep_unique_by == "group") 
      df <- df
    else if (pep_unique_by == "protein") 
      df <- dplyr::mutate(df, pep_isunique = pep_literal_unique)
    else if (pep_unique_by == "none") 
      df <- dplyr::mutate(df, pep_isunique = TRUE)
    else 
      df <- df
    
    df <- df %>% 
      dplyr::select(-which(names(.) %in% c("pep_seq_term", "pep_n_prot_accs", 
                                           "pep_n_genes")))
  })
  
  # 0, not NA since we known... even pep_tot_int is NA
  df <- df %>% 
    dplyr::mutate(pep_unique_int = 
                    ifelse(pep_literal_unique, pep_tot_int, 0)) %>% 
    dplyr::mutate(pep_razor_int = 
                    ifelse(pep_razor_unique, pep_tot_int, 0)) %>% 
    reloc_col_after("pep_unique_int", "pep_tot_int") %>% 
    reloc_col_after("pep_razor_int", "pep_unique_int")
  
  stopifnot(all(c("pep_isunique", "pep_literal_unique", "pep_razor_unique", 
                  "pep_tot_int", "pep_unique_int", "pep_razor_int") %in% 
                  names(df)))
  
  # (2.6) add peptide properties and prot_cover, prot_icover
  df$pep_miss <- NULL
  
  df <- df %>% 
    annotPeppos() %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq), 
                                    stringr::str_count(pep_seq, "[KR]") - 1,
                                    stringr::str_count(pep_seq, "[KR]"))) %>% 
    add_prot_icover(id = group_pep_by) %>% 
    { if (!("prot_cover" %in% names(.) && length(filelist) == 1L)) 
        calc_cover(., id = !!rlang::sym(group_pep_by)) 
      else 
        dplyr::mutate(., prot_cover = prot_cover/100) } 

  # (2.7) apply parsimony
  # allow PSM duplicates that are different at `pep_var_mod_pos` if there are any;
  # mostly not with current Mascot but possible later redundancy of 
  #   one query -> multiple `pep_var_mod_pos`
  
  # uniq_by <- c("RAW_File", "pep_query", "pep_seq")
  uniq_by <- c("RAW_File", "pep_query", "pep_seq", "pep_var_mod_pos")

  if (length(unique(df$dat_file)) >= 1L) 
    uniq_by <- c(uniq_by, "dat_file")
  
  # redundancy at prot_acc also removed (kept the heaviest one)
  df <- df %>% dplyr::arrange(pep_rank, -pep_isbold, -prot_mass)
  rows <- !duplicated(df[, uniq_by])
  df <- df[rows, ]
  rm(list = c("rows"))
  
  df <- try(
    dplyr::arrange(df, prot_hit_num, prot_family_member, pep_start, pep_end)
  )
  
  ## before
  #      prot_acc   I126 I127N   I127C
  # 1  HBB1_MOUSE  29410 20920   27430
  # 2  HBE_MOUSE   29410 20920   27430
  # 3  HBE_MOUSE  137000 38330 1115000
  
  ## after
  #      prot_acc   I126 I127N   I127C
  # 1  HBB1_MOUSE  29410 20920   27430
  # 2  HBE_MOUSE  137000 38330 1115000
  
  # remove subset proteins
  ## one example: I <-> L
  # (No prot_family_member)
  # prot_hit_num	prot_family_member	prot_acc	pep_rank	pep_isbold	pep_isunique	pep_seq
  # 1	            1	                  P04114	1	        1	          1	            AAIQALR
  # 1		                              Q9BY43	1	        0	          1	            AALQALR

  # (2.8) adds columns pep_n_psm, prot_n_psm, prot_n_pep
  # prot_n_psm and pep_n_psm may be inflated depends on the PSM redundancy 
  df <- df %>% add_quality_cols(!!group_psm_by, !!group_pep_by, uniq_by = NULL)
  
  # (2.9) update original Mascot fields
  # (e.g. updated counts after data merging)

  .saveCall <- TRUE
  
  invisible(df)
}


#' Finds the plexes according to the Mascot outputs.
#' 
#' @param df A data frame of Mascot outputs with headers being removed.
find_mascot_plex <- function (df = NULL) 
{
  to_126 <- grep("/126", df[1, ], fixed = TRUE)
  len <- length(to_126)
  
  if (len) len + 1L else 0L
}


#' Defines the ratio columns of Mascot outputs.
#' 
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
#' @param nrow Positive integer; the number of rows in the data.
make_mascot_ratios <- function(TMT_plex = 10L, nrow = 1L) 
{
  if (TMT_plex == 18) {
    col_keys <- c("127N/126", "127C/126", "128N/126", "128C/126", "129N/126", 
                  "129C/126", "130N/126", "130C/126", "131N/126", "131C/126", 
                  "132N/126", "132C/126", "133N/126", "133C/126", "134N/126", 
                  "134C/126", "135N/126")
  }
  else if (TMT_plex == 16) {
    col_keys <- c("127N/126", "127C/126", "128N/126", "128C/126", "129N/126", 
                  "129C/126", "130N/126", "130C/126", "131N/126", "131C/126", 
                  "132N/126", "132C/126", "133N/126", "133C/126", "134N/126")
  } 
  else if (TMT_plex == 11) {
    col_keys <- c("127N/126", "127C/126", "128N/126", "128C/126", "129N/126", 
                  "129C/126", "130N/126", "130C/126", "131N/126", "131C/126")
  } 
  else if (TMT_plex == 10) {
    col_keys <- c("127N/126", "127C/126", "128N/126", "128C/126", "129N/126", 
                  "129C/126", "130N/126", "130C/126", "131/126")
  } 
  else if(TMT_plex == 6) {
    col_keys <- c("127/126", "128/126", "129/126", "130/126", "131/126")
  } 
  else {
    col_keys <- NULL
  }
  
  col_keys <- rep(col_keys, each = 2) 
  col_keys[which(seq_along(col_keys)%%2 == 0)] <- NA
  
  suppressMessages(purrr::map(col_keys, rep, each = nrow) %>% 
                     dplyr::bind_cols() %>% 
                     `colnames<-`(""))
}


#' Defines the intensity columns of Mascot outputs.
#' 
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
#' @param nrow Positive integer; the number of rows in the data.
make_mascot_intensities <- function(TMT_plex = 10L, nrow = 1L) 
{
  if (TMT_plex == 18) {
    col_keys <- c("126", "127N", "127C", "128N", "128C", "129N", 
                  "129C", "130N", "130C", "131N", "131C", 
                  "132N", "132C", "133N", "133C", "134N", 
                  "134C", "135N")
  }
  else if (TMT_plex == 16) {
    col_keys <- c("126", "127N", "127C", "128N", "128C", "129N", 
                  "129C", "130N", "130C", "131N", "131C", 
                  "132N", "132C", "133N", "133C", "134N")
  } 
  else if (TMT_plex == 11) {
    col_keys <- c("126", "127N", "127C", "128N", "128C", "129N", 
                  "129C", "130N", "130C", "131N", "131C")
  } 
  else if (TMT_plex == 10) {
    col_keys <- c("126", "127N", "127C", "128N", "128C", 
                  "129N", "129C", "130N", "130C", "131")
  } 
  else if(TMT_plex == 6) {
    col_keys <- c("126", "127", "128", "129", "130", "131")
  } 
  else {
    col_keys <- NULL
  }
  
  col_keys <- rep(col_keys, each = 2) 
  col_keys[which(seq_along(col_keys)%%2 == 0)] <- -1
  
  suppressMessages(purrr::map(col_keys, rep, each = nrow) %>% 
                     `names<-`(paste0("nm_", seq_along(.))) %>% 
                     dplyr::bind_cols() %>% 
                     `colnames<-`(""))
}

#' Pads Mascot TMT channels to the highest plex
#' 
#' @param file An intermediate PSM table.
#' @param ... filter_dots.
pad_mascot_channels <- function(file = NULL, ...) 
{
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dat_dir <- get_gl_dat_dir()
  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  base_name <- gsub("_hdr_rm\\.csv$", "", file)
  
  df <- read.delim(file.path(dat_dir, "PSM/cache", file), 
                   sep = ',', check.names = FALSE, 
                   header = TRUE, stringsAsFactors = FALSE, 
                   quote = "\"", fill = TRUE , skip = 0) %>% 
    filters_in_call(!!!filter_dots)
  
  # parses `RAW_File`
  l <- df$pep_scan_title[1]
  
  if (grepl("File:~", l, fixed = TRUE)) { # MSConvert
    df <- df %>% 
      dplyr::mutate(RAW_File = pep_scan_title) %>% 
      dplyr::mutate(RAW_File = gsub("^[^ ]+ File:\\~(.*)\\.raw\\~.*", 
                                    "\\1", RAW_File)) %>% 
      dplyr::mutate(RAW_File = gsub("^[^ ]+ File:\\~(.*)\\.d\\~.*", 
                                    "\\1", RAW_File))
  } 
  else if (grepl("File: ~", l, fixed = TRUE)) { # PD
    df <- df %>% 
      dplyr::mutate(RAW_File = gsub("\\\\", "/", pep_scan_title), 
                    RAW_File = gsub('^.*/([^/]*)\\.raw\\~.*', '\\1', RAW_File), 
                    RAW_File = gsub('^.*/([^/]*)\\.d\\~.*', '\\1', RAW_File)) 
  }
  else if (grepl(", File:", l, fixed = TRUE)) { # custom from Bruker MGF, add "~" later...
    df <- df %>% 
      dplyr::mutate(RAW_File = gsub("\\\\", "/", pep_scan_title), 
                    RAW_File = gsub('^.*, File:(.*)\\.(raw|d)\\, .*', '\\1', RAW_File)) 
  }
  else {
    stop("Mascot RAW_File name in \"pep_scan_title\" is not in the format of ", 
         "\"MSConvert\" or \"Proteome DisCoverer\".\n", 
         "Contact the developer about the new format.")
  }
  
  rm(list = "l")
  
  # with some refseq_acc
  df <- df %>%
    dplyr::mutate(prot_acc = gsub("\\.[0-9]*$", "", prot_acc)) 
  
  # compares plexes between Mascot outputs and those in metadata
  this_plex <- find_mascot_plex(df)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  if (this_plex < 0L)
    stop("The multiplexity of data is less than 0 (i.e., TMT > 0, LFQ = 0).")

  if (this_plex > TMT_plex) {
    stop("\nPSM multiplexity is \"", this_plex, "\" with ", file, ".\n",
         "Metadata multiplexity is however smaller at \"", TMT_plex, "\".\n", 
         "Don't know which channels to exclude from PSM file.", 
         call. = FALSE)
  }
  
  if (TMT_plex == 0L) {
    if (!file.exists(file.path(dat_dir, "PSM/cache/", 
                               paste0(base_name, "_queryInfo.rds")))) {
      warning("Query file not found for LFQ.\n", 
              "Make sure that \"Query Level Information\" -> \"Raw peptide match data\" ", 
              "was checked \nwhen exporting Mascot PSMs.", 
              call. = FALSE)
    }
  }
  
  if ((this_plex > 0L) && (this_plex < TMT_plex)) {
    nms <- set_mascot_colnms(TMT_plex) %>% stringr::str_split(",", simplify = TRUE)
    nms_ratio <- nms %>% .[grepl("R[0-9]{3}[NC]{0,1}$", .)]
    nms_int <- nms %>% .[grepl("I[0-9]{3}[NC]{0,1}$", .)]
    
    pos <- find_padding_pos(this_plex, TMT_plex)
    nas <- data.frame(R126 = rep(NA, nrow(df)))

    df_ratio <- local({
      df_ratio <- df %>% dplyr::select(grep("^R[0-9]{3}[NC]{0,1}$", names(.)))
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
      names(holders_ratio) <- nms_ratio
      
      holders_ratio
    })

    df_int <- local({
      df_int <- df %>% dplyr::select(grep("^I[0-9]{3}[NC]{0,1}$", names(.)))
      
      for (idx in seq_along(pos)) {
        df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))
      }
      
      holders_int <- make_mascot_intensities(TMT_plex, nrow(df))
      
      for (idx in seq_len(TMT_plex)) {
        holders_int[, 2*idx] <- df_int[, idx]
      }
      
      names(holders_int) <- nms_int
      
      holders_int
    })
    
    df <- dplyr::bind_cols(
      df %>% select(-grep("[IR]{1}[0-9]{3}[NC]{0,1}$", names(.))), 
      df_ratio,
      df_int,
    )
  }
  
  df$dat_file <- base_name

  invisible(df)
}


#' Pads Mascot PSM exports
#'
#' Intensity or ratio columns are handled during channel padding and thus
#' excluded.
#'
#' @param df An intermediate PSM table from Mascot export.
pad_mascot_fields <- function(df = NULL) 
{
  local({
    list_a <- list.files(path = file.path(dat_dir), 
                         pattern = "^F[0-9]+.*\\.csv$")
    
    list_b <- list.files(path = file.path(dat_dir, "PSM/cache"), 
                         pattern = "^F[0-9]+.*\\_hdr_rm\\.csv$")
    
    len_a <- length(list_a)
    len_b <- length(list_b)
    
    if (len_a != len_b) {
      stop("\nThe number of PSM files under ", file.path(dat_dir), 
           " is ", len_a, ".\n", 
           "The number of PSM files under ", file.path(dat_dir, "PSM/cache"), 
           " is ", len_b, ".", 
           call. = FALSE)
    }
  })
  
  ncols <- purrr::map_dbl(df, ncol)
  
  nms <- lapply(df, function (x) {
    names(x) %>% 
      .[!grepl("[IR]{1}[0-9]{3}[NC]{0,1}$", .)] %>% 
      .[. != "dat_file"]
  })
  
  nms_union <- purrr::reduce(nms, union) 
  
  nms_union <- c(nms_union[nms_union != "pep_scan_title"], 
                 nms_union[nms_union == "pep_scan_title"])
  
  snms_union <- sort(nms_union)
  
  ok_nms <- purrr::map_lgl(nms, function (x) {
    identical(sort(x), snms_union)
  })
  
  if (!all(ok_nms)) {
    warning("Inequal numbers or names of columns deteted: \n", 
            paste(ncols, collapse = "\n"), 
            call. = FALSE)
    
    for (i in seq_along(df)) {
      nms_i <- nms[[i]]
      df_i <- df[[i]]
      
      missing_nms <- setdiff(nms_union, nms_i)
      
      if (length(missing_nms)) {
        for (nm in missing_nms) df_i[[nm]] <- NA
      } 
      
      df[[i]] <- dplyr::bind_cols(
        df_i[nms_union], 
        df_i %>% 
          .[, ! names(.) %in% nms_union, drop = FALSE] %>% 
          .[, ! names(.) == "dat_file", drop = FALSE], 
        df_i["dat_file"],
      )
    }
  }
  
  invisible(df)
}


#' Parallel PSM cleanup
#' 
#' @param file A PSM file (name).
#' @inheritParams load_expts
#' @inheritParams cleanupPSM
#' @inheritParams TMT_levels
psm_mcleanup <- function(file = NULL, rm_outliers = FALSE, 
                         group_psm_by = "pep_seq", dat_dir = NULL, 
                         TMT_plex = 10L, rm_allna = FALSE) 
{
  df <- suppressWarnings(
    readr::read_csv(file.path(dat_dir, "PSM/cache", file), 
                    col_types = get_col_types(), 
                    show_col_types = FALSE)
  )
  
  stopifnot(group_psm_by %in% names(df))
  
  # e.g. TMT_plex is 10 but actually one sample and nine empties
  if (TMT_plex == 0L || sum(grepl("^I[0-9]{3}[NC]{0,1}", names(df))) == 1) {
    fn <- paste0(gsub(".csv", "", file), "_Clean.txt")
    df$psm_index <- NULL
    readr::write_tsv(df, file.path(dat_dir, "PSM/cache", fn))

    message(file, " processed (no PSM cleanup for MS1-based LFQ).")
    
    return(fn)
  } 
  
  # remove all "-1" ratio rows
  df <- local({
    N <- sum(grepl("^R[0-9]{3}", names(df)))
    
    df <- df %>%
      dplyr::mutate(.n = rowSums(.[grep("^R[0-9]{3}", names(.))] == -1, na.rm = TRUE)) %>% 
      { if (rm_allna) dplyr::filter(., .n != N) else . } %>% 
      dplyr::select(-.n)
  })

  set_idx <- file %>% 
    gsub("TMTset(\\d+)_.*", "\\1", .) %>% 
    as.integer()
  
  injn_idx <- file %>% 
    gsub("^TMTset\\d+_LCMSinj(\\d+)\\.csv$", "\\1", .) %>% 
    as.integer()
  
  channelInfo <- channelInfo(dat_dir = dat_dir, 
                             set_idx = set_idx, 
                             injn_idx = injn_idx)

  # add column "R126"
  df <- local({
    pos_af <- min(grep("^R[0-9]{3}", names(df)))
    df$R126 <- 1
    
    df <- cbind.data.frame(df[, 1:(pos_af-1)], 
                           R126 = df$R126, 
                           df[, (pos_af):(ncol(df)-1)]) %>%
      dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, 
                       ~ replace(.x, , NA)) %>% 
      { if (rm_allna) dplyr::filter(., rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) 
        else . } %>% 
      dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, 
                       ~ replace(.x, , 0)) %>%
      dplyr::mutate_at(.vars = which(names(.)=="R126")-1+channelInfo$emptyChannels, 
                       ~ replace(.x, , NA)) %>%
      # "> 1" not "0" 
      { if (rm_allna) dplyr::filter(., rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 1) 
        else . }
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
      
      dfw_split <- do.call("rbind", 
                           lapply(dfw_split, locate_outliers, range_colRatios)) %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         ~ replace(.x, is.infinite(.x), NA)) %>%
        tidyr::unite(pep_seq_i, !!group_psm_by, psm_index, sep = ":") %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         ~ replace(.x, !is.na(.x), 1))
      
      # column order NOT altered other than the removal of "psm_index"
      
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
      { if (rm_allna) dplyr::filter(., rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 1) 
        else . } %>% 
      dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 + channelInfo$emptyChannels, 
                       ~ replace(.x, , 0))
  } 
  else {
    df <- df %>% dplyr::select(-psm_index)
    
    df <- dplyr::bind_cols(
      df[, !grepl("^R[0-9]{3}|^I[0-9]{3}", names(df)), drop = FALSE],
      df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df)), drop = FALSE], 
    ) %>% 
      dplyr::mutate_at(.vars = which(names(.)=="I126") - 1 +
                         channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
      dplyr::mutate_at(.vars = which(names(.)=="R126") - 1 +
                         channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
      { if (rm_allna) dplyr::filter(., rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 1) 
        else . }
  }
  
  fn <- paste0(gsub("\\.csv$", "", file), "_Clean.txt")

  readr::write_tsv(df, file.path(dat_dir, "PSM/cache", fn))
  
  message(file, " processed.")
  
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
cleanupPSM <- function(rm_outliers = FALSE, group_psm_by = "pep_seq", 
                       rm_allna = FALSE, parallel = TRUE) 
{
  dat_dir <- get_gl_dat_dir()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(path = file.path(dat_dir, "PSM/cache"),
                         pattern = "^TMT.*LCMS.*\\.csv$") %>% 
    reorder_files()
  
  n_files <- length(filelist)
  
  if (parallel && (n_files > 1L)) {
    n_cores <- min(parallel::detectCores(), n_files)
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(
      cl,
      c("add_cols_at", 
        "reloc_col_before"), 
      envir = environment(proteoQ:::add_cols_at)
    )
    
    suppressWarnings(
      silent_out <- parallel::clusterApply(
        cl, filelist, psm_mcleanup, 
        rm_outliers, 
        group_psm_by, 
        dat_dir, 
        TMT_plex, 
        rm_allna)
    )
    
    rm(list = "silent_out")
    
    parallel::stopCluster(cl)
  } 
  else {
    purrr::walk(filelist, psm_mcleanup, 
                rm_outliers, 
                group_psm_by, 
                dat_dir, 
                TMT_plex, 
                rm_allna)
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
mcPSM <- function(df = NULL, set_idx = 1L, injn_idx = 1L, mc_psm_by = "peptide", 
                  group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                  rm_allna = FALSE) 
{
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  
  channelInfo <- channelInfo(dat_dir = dat_dir, 
                             set_idx = set_idx, 
                             injn_idx = injn_idx)
  
  label_scheme_sub <- label_scheme_full %>% 
    dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx) 
  
  dfw <- df %>% 
    { if (rm_allna) .[rowSums(!is.na(.[grepl("^R[0-9]{3}[NC]{0,1}", names(.))])) > 1, ] 
      else . } %>%
    dplyr::arrange(pep_seq, prot_acc) 
  
  if (length(channelInfo$emptyChannels)) {
    dfw <- dfw %>%
      dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 +
                         channelInfo$emptyChannels, ~ replace(.x, , NaN))
  }

  col_sample <- grep("^I[0-9]{3}[NC]{0,1}", names(dfw))
  
  if (length(channelInfo$refChannels)) {
    ref_index <- channelInfo$refChannels
  } 
  else {
    ref_index <- channelInfo$labeledChannels
  }
  
  dfw <- sweep(dfw[, col_sample], 1,
               rowMeans(dfw[col_sample[ref_index]], na.rm = TRUE), "/") %>%
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
  } 
  else if (mc_psm_by == "protein") {
    cf <- dfw %>% 
      dplyr::select(group_pep_by, grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>% 
      dplyr::summarise_all(~ median(.x, na.rm = TRUE)) %>% 
      dplyr::summarise_at(.vars = grep("^log2_R[0-9]{3}[NC]{0,1}", names(.)), 
                          ~ median(.x, na.rm = TRUE)) %>% 
      unlist()
  } 
  else if (mc_psm_by == "psm") {
    cf <- apply(dfw[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  } 
  else {
    warning("Unknown setting in \"mc_psm_by\"; the default will be used", 
            call. = FALSE)
  }

  dfw <- sweep(dfw[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- sweep(dfw[, grepl("^I[0-9]{3}", names(dfw)), drop = FALSE], 2, 2^cf, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  # not yet SD scaling at PSM levels
  if (FALSE) {
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
  
  dfw %>%
    reorderCols(endColIndex = grep("[RI][0-9]{3}", names(.)), 
                col_to_rn = "pep_seq_mod") %>%
    na_zeroIntensity()
}


#' Annotates PSM results
#'
#' \code{annotPSM} adds fields of annotation to Mascot PSM tables after
#' \code{rmPSMHeaders}, \code{splitPSM} and \code{cleanupPSM}.
#'
#' @param group_psm_by A character string specifying the method in PSM grouping.
#'   At the \code{pep_seq} default, descriptive statistics will be calculated
#'   based on the same \code{pep_seq} groups. At the \code{pep_seq_mod}
#'   alternative, peptides with different variable modifications will be treated
#'   as different species and descriptive statistics will be calculated based on
#'   the same \code{pep_seq_mod} groups.
#' @param group_pep_by A character string specifying the method in peptide
#'   grouping. At the \code{prot_acc} default, descriptive statistics will be
#'   calculated based on the same \code{prot_acc} groups. At the \code{gene}
#'   alternative, proteins with the same gene name but different accession
#'   numbers will be treated as one group.
#' @param mc_psm_by A character string specifying the method in the median
#'   centering of PSM \code{log2FC} across samples. At the \code{peptide}
#'   default, the median description of PSMs (grouped by \code{pep_seq} or
#'   \code{pep_seq_mod} according to \code{group_psm_by}) will be first
#'   calculated and the offsets to zero (of logarithmic 2) will be used for the
#'   centering of PSMs across samples. At \code{mc_psm_by = protein}, the median
#'   description of PSMs (grouped by \code{prot_acc} or \code{gene} according to
#'   \code{group_pep_by}) will be calculated and the corresponding offsets to
#'   zero will be applied. At the \code{mc_psm_by = psm}, PSMs will be median
#'   centered without grouping.
#' @param plot_log2FC_cv Logical; if TRUE, the distributions of the CV of
#'   peptide \code{log2FC} will be plotted. The default is TRUE.
#' @param ... Not currently used.
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams normPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %T>% %$% %<>%
annotPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                     mc_psm_by = "peptide", 
                     fasta = NULL, expt_smry = "expt_smry.xlsx", 
                     plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                     rm_allna = FALSE, type_sd = "log2_R", ...) 
{
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

    out_fn <- sublist %>% 
      purrr::map_chr(~ gsub("_Clean\\.txt$", "_PSM_N.txt", .x))
    
    # --- LCMS injections under the same TMT experiment ---
    
    for (idx in seq_along(sublist)) {
      df <- suppressWarnings(
        readr::read_tsv(file.path(dat_dir, "PSM/cache", sublist[idx]), 
                        col_types = get_col_types(), 
                        show_col_types = FALSE)
      )

      df <- df %>% 
        add_pep_retsd(group_psm_by) %>% 
        add_n_pepexpz(group_psm_by)

      # e.g. TMT 10-plex but 9 are empties
      if (TMT_plex && (sum(grepl("^I[0-9]{3}[NC]{0,1}", names(df))) > 1L)) {
        injn_idx <- sublist[idx] %>% 
          gsub("^TMTset\\d+_LCMSinj(\\d+)_Clean.txt$", "\\1", .) %>% 
          as.integer()
        
        df <- mcPSM(df = df, 
                    set_idx = set_idx, 
                    injn_idx = injn_idx, 
                    mc_psm_by = mc_psm_by, 
                    group_psm_by = group_psm_by, 
                    group_pep_by = group_pep_by, 
                    rm_allna = rm_allna)
        
        rm(list = "injn_idx")
      } 
      else {
        df <- df %>% 
          dplyr::mutate(N_I000 = I000, 
                        R000 = NA_real_, 
                        log2_R000 = NA_real_, 
                        N_log2_R000 = NA_real_)
      }
      
      df <- df %>% 
        dplyr::mutate(!!group_psm_by := as.character(!!rlang::sym(group_psm_by))) %>% 
        calcSD_Splex(id = group_psm_by, type = type_sd) %>% 
        `names<-`(gsub(paste0("^", type_sd), "sd_log2_R", names(.))) %>% 
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
        readr::write_tsv(file.path(dat_dir, "PSM", out_fn[idx]))

      if (plot_rptr_int) {
        try(
          df_int <- df %>% 
            dplyr::select(grep("^I[0-9]{3}", names(.))) %>% 
            rptr_violin(filepath = 
                          file.path(dat_dir, "PSM/rprt_int/mc", 
                                    gsub("_PSM_N\\.txt", "_rprt.png", out_fn[idx])), 
                        width = 8, height = 8)
        )
      }

      if (plot_log2FC_cv && TMT_plex) {
        try(
          df %>% 
            sd_violin(id = !!group_psm_by, 
                      filepath = 
                        file.path(dat_dir, "PSM/log2FC_cv/raw", 
                                  gsub("_PSM_N\\.txt", "_sd.png", out_fn[idx])), 
                      width = 8, height = 8, type = "log2_R", 
                      adjSD = FALSE, is_psm = TRUE, ...)
        )
      }
      
    }
  }
}


#' Standardization of PSM
#'
#' \code{normPSM} standardizes
#' \href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSM} results
#' from database search engines.
#'
#' In each primary output file, "\code{...PSM_N.txt}", values under columns
#' \code{log2_R...} are logarithmic ratios at base 2 in relative to the average
#' intensity of \code{reference(s)} within each multiplex TMT set, or to the
#' row-mean intensity within each plex if no \code{reference(s)} are present.
#' Values under columns \code{N_log2_R...} are \code{log2_R...} with
#' median-centering alignment. Values under columns \code{I...} are raw
#' \code{reporter-ion intensity} from database searches. Values under columns
#' \code{N_I...} are normalized \code{reporter-ion intensity}. Values under
#' columns \code{sd_log2_R...} are the standard deviation of the \code{log2FC}
#' of peptides from ascribing PSMs. Character strings under \code{pep_seq_mod}
#' denote peptide sequences with applicable variable modifications.
#' 
#' \cr \strong{Nomenclature of \code{pep_seq_mod}}:
#' 
#' \tabular{ll}{ \emph{Variable modification}   \tab \emph{Abbreviation}\cr
#' Non-terminal \tab A letter from upper to lower case, e.g., \code{mtFPEADILLK}
#' \cr N-term \tab A hat to the left of a peptide sequence, e.g.,
#' \code{^QDGTHVVEAVDATHIGK} \cr C-term \tab A hat to the right of a peptide
#' sequence, e.g., \code{DAYYNLCLPQRPnMI^} \cr Acetyl (Protein N-term) \tab A
#' underscore to the left of a peptide sequence, e.g., \code{_mAsGVAVSDGVIK}.
#' \cr Amidated (Protein C-term) \tab A underscore to the right of a peptide
#' sequence, e.g., \code{DAYYNLCLPQRPnMI_}. \cr Other (Protein N-term) \tab A
#' tilde to the left of a peptide sequence, e.g., \code{~mAsGVAVSDGVIK} \cr
#' Other (Protein C-term) \tab An tilde to the right of a peptide sequence, e.g.
#' \code{DAYYNLCLPQRPnMI~} \cr }
#'
#' @section \code{Mascot}: Users will export \code{PSM} data from
#'   \href{https://http://www.matrixscience.com/}{Mascot} at a \code{.csv}
#'   format and store them under the file folder indicated by \code{dat_dir}.
#'   The header information should be included during the \code{.csv} export.
#'   The file name(s) should start with the letter \code{'F'} and ended with a
#'   \code{'.csv'} extension \code{(e.g., F004452.csv, F004453_this.csv etc.)}.
#'
#' @section \code{MaxQuant}: Users will copy over \code{msms.txt} file(s) from
#'   \href{https://www.maxquant.org/}{MaxQuant} to the \code{dat_dir} directory.
#'   The file name(s) should start with \code{'msms'} and end with a
#'   \code{'.txt'} extension \code{(e.g., msms.txt, msms_this.txt etc.)}.
#'
#' @section \code{MSFragger}: Users will copy over \code{psm.tsv} file(s) from
#'   \href{http://msfragger.nesvilab.org/}{MSFragger} to the \code{dat_dir}
#'   directory. The file name(s) should start with \code{'psm'} and end with a
#'   \code{'.tsv'} extension \code{(e.g., psm.tsv, psm_this.tsv etc.)}.
#'
#' @section \code{Spectrum Mill}: Users will copy over \code{PSMexport.1.ssv}
#'   file(s) from
#'   \href{https://www.agilent.com/en/products/software-informatics/masshunter-suite/masshunter-for-life-science-research/spectrum-mill}{Spectrum
#'    Mill} to the \code{dat_dir} directory. The file name(s) should start with
#'   \code{'PSMexport'} and end with a \code{'.ssv'} extension \code{(e.g.,
#'   PSMexport.ssv, PSMexport_this.ssv etc.)}.
#' 
#' @param rm_allna Logical; if TRUE, removes data rows that are exclusively NA
#'   across ratio columns of \code{log2_R126} etc. The setting also applies to
#'   \code{log2_R000} in LFQ.
#' @param type_sd Character string; the type of log2Ratios for SD calculations.
#'   The value is one \code{log2_R}, \code{N_log2_R} or \code{Z_log2_R}.
#' @param ... \code{filter_}: Variable argument statements for the filtration of
#'   data rows. Each statement contains to a list of logical expression(s). The
#'   \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'   the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'   For example, \code{pep_expect} is a column key present in \code{Mascot} PSM
#'   exports and \code{filter_psms_at = exprs(pep_expect <= 0.1)} will remove
#'   PSM entries with \code{pep_expect > 0.1}.
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams splitPSM_mq
#' @inheritParams splitPSM_mf
#' @inheritParams cleanupPSM
#' @inheritParams annotPSM
#' @seealso
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
#'  \emph{Variable arguments of \code{filter_...}} \cr 
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
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'  
#' @section \code{Variable arguments and data files}: Variable argument (vararg)
#'   statements of \code{filter_} and \code{arrange_} are available in
#'   \code{proteoQ} for flexible filtration and ordering of data rows, via
#'   functions at users' interface. To take advantage of the feature, users need
#'   to be aware of the column keys in input files. As indicated by their names,
#'   \code{filter_} and \code{filter2_} perform row filtration against column
#'   keys from a primary data file, \code{df}, and secondary data file(s),
#'   \code{df2}, respectively. The same correspondence is applicable for
#'   \code{arrange_} and \code{arrange2_} varargs. \cr \cr Users will typically
#'   employ either primary or secondary vararg statements, but not both. In the
#'   more extreme case of \code{gspaMap(...)}, it links \code{\link{prnGSPA}}
#'   findings in \code{df2} to the significance \code{pVals} and abundance fold
#'   changes in \code{df} for volcano plot visualizations by gene sets. The
#'   table below summarizes the \code{df} and the \code{df2} for varargs in
#'   \code{proteoQ}.
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
#' @return Outputs are interim and final PSM tables under the directory of
#'   \code{PSM} sub to \code{dat_dir}. Primary results are in
#'   \emph{standardized} PSM tables of \code{TMTset1_LCMSinj1_PSM_N.txt,
#'   TMTset2_LCMSinj1_PSM_N.txt, etc.} The indexes of TMT experiment and LC/MS
#'   injection are indicated in the file names.
#' @example inst/extdata/examples/normPSM_.R
#' @import dplyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @export
normPSM <- function(dat_dir = NULL, 
                    expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx", 
                    fasta = NULL, entrez = NULL, 
                    group_psm_by = c("pep_seq", "pep_seq_mod"), 
                    group_pep_by = c("prot_acc", "gene"), 
                    pep_unique_by = c("group", "protein", "none"), 
                    mc_psm_by = c("peptide", "protein", "psm"), 
                    scale_rptr_int = FALSE, 
                    rptr_intco = 0, rptr_intrange = c(0, 100), 
                    rm_craps = FALSE, rm_krts = FALSE, rm_outliers = FALSE, 
                    rm_allna = FALSE, type_sd = c("log2_R", "N_log2_R", "Z_log2_R"), 
                    purge_phosphodata = TRUE, 
                    annot_kinases = FALSE, 
                    plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                    use_lowercase_aa = TRUE, 
                    corrected_int = TRUE, rm_reverses = TRUE, 
                    parallel = TRUE, ...) 
{
  old_opts <- options()
  options(warn = 1L)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = environment())) {
      if (.savecall) {
        mget(names(formals()), envir = environment(), inherits = FALSE) %>% 
          c(dots) %>% 
          save_call("normPSM")
      }
    } 
    else {
      warning("The current call is not saved.", call. = TRUE)
    }, 
    add = TRUE
  )
  
  # ---
  dots <- rlang::enexprs(...)

  if (is.null(dat_dir)) {
    dat_dir <- get_gl_dat_dir()
  } 
  else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
    message("Variable \"dat_dir\" added to the Global Environment.")
  }

  if (is.null(fasta)) 
    stop("Path(s) to fasta file(s) cannot be empty.", call. = FALSE)
  
  # ---
  group_psm_by <- rlang::enexpr(group_psm_by)
  oks <- eval(formals()[["group_psm_by"]])
  
  group_psm_by <- if (length(group_psm_by) > 1L) 
    oks[[1]]
  else 
    rlang::as_string(group_psm_by)
  
  if (!group_psm_by %in% oks)
    stop("\"group_psm_by\" is not one of ", paste(oks, collapse = ", "), 
         call. = FALSE)
  
  if (length(group_psm_by) != 1L) 
    stop("Length of \"group_psm_by\" is not one.", call. = FALSE)

  rm(list = c("oks"))

  # ---
  group_pep_by <- rlang::enexpr(group_pep_by)
  
  oks <- eval(formals()[["group_pep_by"]])
  
  group_pep_by <- if (length(group_pep_by) > 1L) 
    oks[[1]]
  else 
    group_pep_by <- rlang::as_string(group_pep_by)
  
  if (!group_pep_by %in% oks)
    stop("\"group_pep_by\" is not one of ", paste(oks, collapse = ", "), 
         call. = FALSE)
  
  if (length(group_pep_by) != 1L) 
    stop("Length of \"group_pep_by\" is not one.", call. = FALSE)
  
  rm(list = c("oks"))
  
  # ---
  type_sd <- rlang::enexpr(type_sd)
  
  if (length(type_sd) > 1L) 
    type_sd <- "log2_R"
  else
    type_sd <- rlang::as_string(type_sd)
  
  stopifnot(type_sd %in% c("log2_R", "N_log2_R", "Z_log2_R"), 
            length(type_sd) == 1L)

  # ---
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
  
  type <- find_search_engine(dat_dir)
  
  # ---
  pep_unique_by <- rlang::enexpr(pep_unique_by)
  
  oks <- eval(formals()[["pep_unique_by"]])
  
  pep_unique_by <- if (length(pep_unique_by) > 1L) 
    oks[[1]]
  else 
    rlang::as_string(pep_unique_by)
  
  stopifnot(pep_unique_by %in% oks, length(pep_unique_by) == 1L)
  
  rm(list = c("oks"))
  
  # ---
  mc_psm_by <- rlang::enexpr(mc_psm_by)
  
  oks <- eval(formals()[["mc_psm_by"]])
  
  mc_psm_by <- if (length(mc_psm_by) > 1L) 
    oks[[1]]
  else 
    rlang::as_string(mc_psm_by)
  
  if (!mc_psm_by %in% oks)
    stop("\"mc_psm_by\" is not one of ", paste(oks, collapse = ", "), 
         call. = FALSE)
  
  if (length(mc_psm_by) != 1L) 
    stop("Length of \"mc_psm_by\" is not one.", call. = FALSE)
  
  rm(list = c("oks"))

  # ---
  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))
  
  stopifnot(vapply(c(rptr_intco), is.numeric, logical(1)))
  
  stopifnot(vapply(c(corrected_int, 
                     rm_reverses, 
                     rm_craps, 
                     rm_krts, 
                     rm_outliers, 
                     rm_allna, 
                     annot_kinases, 
                     plot_rptr_int, 
                     plot_log2FC_cv, 
                     use_lowercase_aa, 
                     purge_phosphodata, 
                     scale_rptr_int), 
                   rlang::is_logical, logical(1)))
  
  reload_expts()
  
  if (type == "pq") {
    df <- splitPSM_pq(group_psm_by = group_psm_by, 
                      group_pep_by = group_pep_by, 
                      fasta = fasta, 
                      entrez = entrez, 
                      pep_unique_by = pep_unique_by, 
                      scale_rptr_int = scale_rptr_int, 
                      rm_craps = rm_craps, 
                      rm_krts = rm_krts, 
                      rm_allna = rm_allna, 
                      purge_phosphodata = purge_phosphodata, 
                      annot_kinases = annot_kinases, 
                      plot_rptr_int = plot_rptr_int, 
                      rptr_intco = rptr_intco, 
                      rptr_intrange = rptr_intrange, 
                      use_lowercase_aa = use_lowercase_aa, 
                      parallel = parallel, ...)
  } 
  else if (type == "mascot") {
    rmPSMHeaders(parallel = parallel)
    
    df <- splitPSM(group_psm_by = group_psm_by, 
                   group_pep_by = group_pep_by, 
                   fasta = fasta, 
                   entrez = entrez, 
                   pep_unique_by = pep_unique_by, 
                   scale_rptr_int = scale_rptr_int, 
                   rm_craps = rm_craps, 
                   rm_krts = rm_krts, 
                   rm_allna = rm_allna, 
                   purge_phosphodata = purge_phosphodata, 
                   annot_kinases = annot_kinases, 
                   plot_rptr_int = plot_rptr_int, 
                   rptr_intco = rptr_intco, 
                   rptr_intrange = rptr_intrange, 
                   use_lowercase_aa = use_lowercase_aa, 
                   parallel = parallel, ...)
  } 
  else if (type == "mq") {
    df <- splitPSM_mq(group_psm_by = group_psm_by, 
                      group_pep_by = group_pep_by, 
                      fasta = fasta, 
                      entrez = entrez, 
                      pep_unique_by = pep_unique_by, 
                      scale_rptr_int = scale_rptr_int, 
                      corrected_int = corrected_int, 
                      rm_craps = rm_craps, 
                      rm_krts = rm_krts, 
                      rm_allna = rm_allna, 
                      rm_reverses = rm_reverses, 
                      purge_phosphodata = purge_phosphodata, 
                      annot_kinases = annot_kinases, 
                      plot_rptr_int = plot_rptr_int, 
                      rptr_intco = rptr_intco, 
                      rptr_intrange = rptr_intrange, 
                      use_lowercase_aa = use_lowercase_aa, 
                      parallel = parallel, ...)
  } 
  else if (type == "sm") {
    df <- splitPSM_sm(group_psm_by = group_psm_by, 
                      group_pep_by = group_pep_by, 
                      fasta = fasta, 
                      entrez = entrez, 
                      pep_unique_by = pep_unique_by, 
                      scale_rptr_int = scale_rptr_int, 
                      rm_craps = rm_craps, 
                      rm_krts = rm_krts, 
                      rm_allna = rm_allna, 
                      purge_phosphodata = purge_phosphodata, 
                      annot_kinases = annot_kinases, 
                      plot_rptr_int = plot_rptr_int, 
                      rptr_intco = rptr_intco, 
                      rptr_intrange = rptr_intrange, 
                      use_lowercase_aa = use_lowercase_aa, 
                      parallel = parallel, ...)
  } 
  else if (type == "mf") {
    df <- splitPSM_mf(group_psm_by = group_psm_by, 
                      group_pep_by = group_pep_by, 
                      fasta = fasta, 
                      entrez = entrez, 
                      pep_unique_by = pep_unique_by, 
                      scale_rptr_int = scale_rptr_int, 
                      rm_craps = rm_craps, 
                      rm_krts = rm_krts, 
                      rm_allna = rm_allna, 
                      purge_phosphodata = purge_phosphodata, 
                      annot_kinases = annot_kinases, 
                      plot_rptr_int = plot_rptr_int, 
                      rptr_intco = rptr_intco, 
                      rptr_intrange = rptr_intrange, 
                      use_lowercase_aa = use_lowercase_aa, 
                      parallel = parallel, ...)
  }

  procPSMs(df = df, 
           scale_rptr_int = scale_rptr_int, 
           rptr_intco = rptr_intco, 
           rptr_intrange = rptr_intrange, 
           rm_craps = rm_craps, 
           rm_krts = rm_krts, 
           rm_allna = rm_allna, 
           annot_kinases = annot_kinases, 
           plot_rptr_int = plot_rptr_int, 
           parallel = parallel)
  
  cleanupPSM(rm_outliers = rm_outliers, 
             group_psm_by = group_psm_by, 
             rm_allna = rm_allna, 
             parallel = parallel)
  
  annotPSM(group_psm_by = group_psm_by, 
           group_pep_by = group_pep_by, 
           mc_psm_by = mc_psm_by, 
           fasta = fasta, 
           expt_smry = expt_smry, 
           plot_rptr_int = plot_rptr_int, 
           plot_log2FC_cv = plot_log2FC_cv, 
           rm_allna = rm_allna, 
           type_sd = type_sd, 
           ...)
  
  .savecall <- TRUE
}


#' Removes columns from PSM table
#' 
#' @param df PSM data
#' @inheritParams annotPSM
#' @inheritParams channelInfo
rm_cols_mqpsm <- function(df = NULL, group_psm_by = "pep_seq", set_idx = 1L) 
{
  df <- local({
    uniq_by <- c("pep_seq_mod", "Charge", "pep_ret_range", "raw_file")
    
    df %>% 
      tidyr::unite(uniq_id, uniq_by, sep = "@", remove = FALSE) %>% 
      dplyr::group_by(uniq_id) %>% 
      dplyr::arrange(-pep_expect) %>% 
      dplyr::filter(row_number() == 1) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-uniq_id) %>% 
      dplyr::mutate(TMT_Set = set_idx)
  })
  
  df <- local({
    if ("Proteins" %in% names(df)) {
      if ("TMT_Set" %in% names(df)) {
        df <- df %>% reloc_col_before("Proteins", "TMT_Set")
      } 
      else {
        df <- df %>% reloc_col_after_last("Proteins")
      }
    }
    
    if ("Protein Group Ids" %in% names(df)) {
      if ("TMT_Set" %in% names(df)) {
        df <- df %>% reloc_col_before("Protein Group Ids", "TMT_Set")
      } 
      else {
        df <- df %>% reloc_col_after_last("Protein Group Ids")
      }
    }
    
    col_start <- which(names(df) == "Scan Number")
    col_end <- which(names(df) == "Retention Time")
    
    if (length(col_start) && length(col_end)) {
      df <- df %>% dplyr::select(-(col_start : col_end))
    }
    
    df
  })
  
  df <- df %>% 
    dplyr::select(-grep("\\s{1}Probabilities$", names(.))) %>% 
    dplyr::select(-grep("\\s{1}Score\\s{1}Diffs$", names(.))) %>% 
    dplyr::select(-which(names(.) %in% stringr::str_to_title(
      c(
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
        "ID", # "Protein Group Ids", 
        "Peptide ID", "Mod. peptide ID", "Evidence ID", 
        "Length", 
        "Precursor full scan number", "Precursor apex fraction", 
        "Precursor apex offset", "Precursor apex offset time", 
        "Number of matches")
    ))) %>% 
    dplyr::select(-grep(str_to_title("site IDs$"), names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("m/z", stringr::str_to_title("PIF")))) %>% 
    dplyr::select(-which(names(.) %in% str_to_title(
      c("Charge", 
        "Mass", "Mass error [ppm]", "Mass error [Da]", 
        "Combinatorics", 
        "Fraction of total spectrum", 
        "Base peak fraction", 
        "Precursor Apex Fraction", 
        "Intensity coverage", "Peak coverage")
    ))) %>% 
    dplyr::select(-grep(str_to_title("^Reporter mass deviation "), names(.))) 
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c("raw_file")))
  
  if (group_psm_by == "pep_seq_mod") {
    df <- df %>% dplyr::select(-pep_seq)
  } else {
    df <- df %>% dplyr::select(-pep_seq_mod)
  }
}


#' Calculates peptide data for individual TMT experiments
#'
#' Argument \code{injn_idx} does not currently used.
#'
#' @param lfq_ret_tol The tolerance of retention time (in seconds) for the
#'   aggregation of LFQ data.
#' @inheritParams mcPSM
#' @inheritParams PSM2Pep
#' @inheritParams annotPSM
#' @inheritParams channelInfo
#' @inheritParams locate_outliers
#' @inheritParams TMT_levels
calcPeptide <- function(df = NULL, group_psm_by = "pep_seq", 
                        method_psm_pep = "median", group_pep_by = "prot_acc", 
                        dat_dir = NULL, set_idx = 1L, injn_idx = 1L, 
                        TMT_plex = 10L, lfq_ret_tol = 60L, rm_allna = FALSE, 
                        type_sd = "log2_R") 
{
  nms <- names(df)
  
  if (!"prot_acc" %in% nms)
    stop("Column \"prot_acc\" not found.")
  
  if (length(group_psm_by) != 1L)
    stop("Length of \"group_psm_by\" is not one.")
  
  if (!group_psm_by %in% nms)
    stop("Column \"", group_psm_by, "\" not found.")
  
  if (!any(grepl("log2_R[0-9]{3}|I[0-9]{3}", nms)))
    stop("Columns of log2 ratios or intensity not found.")
  
  group_psm_by0 <- group_psm_by
  
  group_ids <- if ("pep_group" %in% nms) unique(df$pep_group) else NULL

  if (length(group_ids) >= 2L) {
    group_psm_by <- c(group_psm_by, "pep_group")
    save(group_ids, file = file.path(dat_dir, "group_ids.rds"))
  }
  
  channelInfo <- channelInfo(dat_dir = dat_dir, 
                             set_idx = set_idx, 
                             injn_idx = injn_idx)
  
  # (currently no rm_allna for LFQ: Maxquant timsTOF LFQ intensities are all NA)
  if (TMT_plex && rm_allna) {
    df <- df %>% 
      dplyr::filter(rowSums(!is.na(.[grepl("^N_log2_R[0-9]{3}", names(.))])) > 0) 
    
    if (!nrow(df)) {
      stop("Fields 'N_log2_R' are all NA at TMT_plex = ", TMT_plex, ".\n", 
           "Is this an error-tolerant search ", 
           "without intensity/ratio values in the psm table?", 
           "May consider replacing 'NA' with '1' to proceed.", 
           call. = FALSE)
    }
  }
  
  df <- df %>%
    dplyr::arrange_at(c(group_psm_by, "prot_acc")) %>% 
    dplyr::select(-grep("^R[0-9]{3}", names(.)))
  
  # --- Mascot ---
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "prot_hit_num", "prot_family_member", "prot_score", 
      "prot_matches", "prot_sequences", 
      "pep_var_mod", "pep_var_mod_pos", "pep_scan_title", 
      "raw_file", "pep_query", "pep_summed_mod_pos", "pep_local_mod_pos", 
      "pep_rank", "pep_isbold", "pep_exp_mz", "pep_exp_mr", "pep_exp_z", 
      "pep_delta", "pep_calc_mr", "pep_homol", "pep_ident", 
      "pep_num_match", "pep_scan_range", 
      "pep_ms2_sumint", "pep_n_ions", 
      "pep_scan_num", "pep_mod_group", "pep_fmod", "pep_vmod", "pep_ivmod", 
      "pep_n_ms2", "pep_rank_nl", 
      "pep_ms2_moverzs", "pep_ms2_ints", 
      "pep_ms2_theos", "pep_ms2_theos2", 
      "pep_ms2_exptints", "pep_ms2_exptints2", 
      "pep_n_matches", "pep_n_matches2", "pep_ms2_deltas", 
      "pep_ms2_ideltas", "pep_ms2_deltas2", "pep_ms2_ideltas2", 
      "pep_ms2_deltas_mean", "pep_ms2_deltas_sd", 
      "pep_ions_first", "pep_ions_second", "pep_ions_third")))
  
  # --- MaxQuant ---
  # ("Charge" also in MSFragger)
  df <- local({
    col_start <- which(names(df) == "Scan Number")
    col_end <- which(names(df) == "Retention Time")
    
    if (length(col_start) && length(col_end)) {
      df <- df %>% dplyr::select(-(col_start : col_end))
    }
    
    df
  })
  
  df <- df %>% 
    dplyr::select(-grep("\\s{1}Probabilities$", names(.))) %>% 
    dplyr::select(-grep("\\s{1}Score\\s{1}Diffs$", names(.))) %>% 
    dplyr::select(-which(names(.) %in% stringr::str_to_title(
      c(
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
        "ID", # "Protein Group Ids", 
        "Peptide ID", "Mod. peptide ID", "Evidence ID", 
        "Length")
    ))) %>% 
    dplyr::select(-grep(stringr::str_to_title("site IDs$"), names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("m/z", "PIF", "Pif"))) %>% 
    dplyr::select(-which(names(.) %in% stringr::str_to_title(
      c("Charge", "Mass", "Mass error [ppm]", 
        "Mass error [Da]", 
        "Combinatorics", 
        "Fraction of total spectrum", 
        "Base peak fraction", "Precursor Intensity", 
        "Precursor Apex Fraction", 
        "Intensity coverage", "Peak coverage")
    ))) %>% 
    dplyr::select(-grep(stringr::str_to_title("^Reporter mass deviation "), 
                        names(.)))
  
  # --- Spectrum Mill ---
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "number", "modifications", 
      "variableSites", "nterm", "previous_aa", "sequence", "next_aa", 
      "cys", "searchCycle", "L/H", "accession_numbers", "entry_name", 
      "matched_parent_mass", "parent_charge", 
      "deltaForwardReverseScore", "percent_scored_peak_intensity", 
      "totalIntensity", "precursorAveragineChiSquared", 
      "precursorIsolationPurityPercent", 
      "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
      "delta_parent_mass", "delta_parent_mass_ppm"))) 
  
  # --- MSFragger ---
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "Spectrum", "Spectrum File", "Ion Mobility", "raw_file", "Peptide Length", 
      "Charge", "Retention", "Observed Mass", "Calibrated Observed Mass", 
      "Observed M/Z", "Calibrated Observed M/Z", "Calculated Peptide Mass", 
      "Calculated M/Z", "Delta Mass", "Number of Missed Cleavages", 
      "Number of Enzymatic Termini", 
      "Ion Mobility", "Assigned Modifications", "Observed Modifications", 
      "Entry Name", "Protein Description")))
  
  
  # summarizes log2FC and intensity from the same `set_idx` 
  #   at one or multiple LCMS series
  if (!TMT_plex) {
    cols <- grep("^I[0-9]{3}[NC]{0,1}", names(df))
    
    if (length(cols))
      ok_intensity <- TRUE
    else {
      warning("Columns of precursor intensities not found.")
      ok_intensity <- FALSE
    }
    
    if ("pep_ret_range" %in% names(df))
      ok_ret <- TRUE
    else {
      warning("\"pep_ret_range\" not found.")
      ok_ret <- FALSE
    }
    
    if (!all(is.na(df[["pep_ret_range"]]))) 
      ok_ret2 <- TRUE
    else {
      warning("Values of \"pep_ret_range\" are all NA.")
      ok_ret2 <- FALSE
    }
    
    ok_lfq_ret <- ok_intensity && ok_ret && ok_ret2

    if (ok_lfq_ret) {
      # subset by top_n intensity
      df <- local({
        n <- switch(method_psm_pep, 
                    lfq_max = 1, 
                    lfq_top_2_sum = 2, 
                    lfq_top_3_sum = 3, 
                    lfq_all = Inf, 
                    2)
        
        # but since I000 are all NA with MaxQuant timsTOF 
        # -> would have to replace NA with 0
        
        df <- df %>% 
          dplyr::mutate(I000 = ifelse(is.na(I000), 0, I000)) %>% 
          dplyr::group_by_at(group_psm_by)

        if (!is.infinite(n)) {
          df <- df  %>% dplyr::top_n(n = n, wt = I000)
        }
        
        df <- df %>%
          dplyr::arrange(-I000) %>% 
          dplyr::ungroup()
      })
      
      # subset by retention time
      df <- local({
        uniq_by <- c(group_psm_by, "pep_ret_range")

        # 1.1 no yet differentiation of neutral losses at different pep_ivmod
        # 1.2 no yet differentiation of different charge states
        # 1.2 no differentiation by raw_file
        # 1.3 no differentiation by dat_file
        # 2.1 the same pep_seq at different prot_acc makes no difference
        
        df_uniq <- unique(df[, c(group_psm_by, "pep_ret_range")])

        # (group_psm_by0 not group_psm_by)
        df_split <- split(df_uniq, df_uniq[[group_psm_by0]])

        rows <- unlist(lapply(df_split, function (x) nrow(x) == 1L))
        dfu <- dplyr::bind_rows(df_split[rows])
        df_split <- df_split[!rows]
        
        dfm <- lapply(df_split, function (x) {
          a <- x$pep_ret_range
          d <- c(0, abs(a[2:length(a)] - a[1]))
          rows <- ifelse(d <= lfq_ret_tol, TRUE, FALSE)
          x[rows, ]
        }) %>% 
          dplyr::bind_rows()
        
        oks <- dplyr::bind_rows(dfu, dfm) %>% 
          tidyr::unite(uniq_id, uniq_by, sep = ".", remove = TRUE) %>% 
          dplyr::mutate(keep. = TRUE)
        
        df <- df %>% 
          tidyr::unite(uniq_id, uniq_by, sep = ".", remove = FALSE) %>% 
          dplyr::left_join(oks, by = "uniq_id") %>% 
          dplyr::filter(keep.) %>% 
          dplyr::select(-c("keep.", "uniq_id"))
      })
      
      df_num <- aggrLFQs(sum)(df, group_psm_by, na.rm = TRUE) %>% 
        dplyr::mutate(log2_R000 = NA, N_log2_R000 = NA)
    }
  }
  
  # works even without precursor intensity (in MaxQuant)
  if (!exists("df_num", envir = environment())) {
    df_num <- switch(method_psm_pep, 
                     mean = 
                       aggrNums(mean)(df, !!rlang::sym(group_psm_by), na.rm = TRUE), 
                     median = 
                       aggrNums(median)(df, !!rlang::sym(group_psm_by), na.rm = TRUE), 
                     weighted_mean = 
                       tmt_wtmean(df, !!rlang::sym(group_psm_by), na.rm = TRUE), 
                     top_3_mean = 
                       TMT_top_n(df_num, !!rlang::sym(group_psm_by), na.rm = TRUE), 
                     lfq_max = 
                       aggrTopn(sum)(df, !!rlang::sym(group_psm_by), 1, na.rm = TRUE), 
                     lfq_top_2_sum = 
                       aggrTopn(sum)(df, !!rlang::sym(group_psm_by), 2, na.rm = TRUE), 
                     lfq_top_3_sum = 
                       aggrTopn(sum)(df, !!rlang::sym(group_psm_by), 3, na.rm = TRUE), 
                     lfq_all = 
                       aggrTopn(sum)(df, !!rlang::sym(group_psm_by), Inf, na.rm = TRUE), 
                     aggrNums(median)(df, !!rlang::sym(group_psm_by), na.rm = TRUE))
  }

  # `pep_unique_int` and `pep_razor_int` recalculated  
  #   from the newly derived `pep_tot_int` (always by sum over `group_psm_by`). 
  # `pep_tot_int` is different to the summary of I000, I126 etc. 
  #   which is according to `method_psm_pep`
  
  # if were to calculate the statistics of `pep_ret_range` using the 
  # PSM with most intensity `pep_tot_int` or `pep_unique_int`,
  # ARRANGE DATA HERE and then summarize...
  
  # form PSM to peptide, median statistics of `pep_n_nl` and `pep_n_exp_z` 
  #   is the same as "the first" statistics (at a small cost of computation, 
  #   but allow the reuse of med_summarise_keys in mergePep)
  df_first <- df %>% 
    dplyr::select(-which(names(.) %in% c("pep_unique_int", 
                                         "pep_razor_int"))) %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    med_summarise_keys(group_psm_by) 
  
  df <- local({
    colnm_before <- find_preceding_colnm(df, "pep_tot_int")
    
    list(df_first, df_num) %>%
      purrr::reduce(dplyr::left_join, by = group_psm_by) %>% 
      reloc_col_after("pep_score", "pep_razor_unique") %>% 
      reloc_col_after("pep_expect", "pep_score") %>% 
      reloc_col_after("pep_ret_sd", "pep_ret_range") %>% 
      reloc_col_after("pep_tot_int", colnm_before)
  })
  
  df <- local({
    stopifnot(all(c("pep_tot_int", "shared_prot_accs", "shared_genes") %in% 
                    names(df)))
    
    if (group_pep_by == "gene")
      col_map <- "shared_genes"
    else
      col_map <- "shared_prot_accs"

    # floating `pep_isunique` can be either `pep_razor_unique` or `pep_literal_unique`
    # so don't use `pep_isunique`
    df %>% 
      dplyr::mutate(pep_unique_int = ifelse(pep_literal_unique, pep_tot_int, 0), 
                    pep_razor_int = ifelse(pep_razor_unique, pep_tot_int, 0))
  })
  
  df <- df %>% 
    reloc_col_after("pep_unique_int", "pep_tot_int") %>% 
    reloc_col_after("pep_razor_int", "pep_unique_int")
  
  df <- cbind.data.frame(df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df)), 
                            drop = FALSE],
                         df[, grepl("I[0-9]{3}", names(df)), 
                            drop = FALSE], 
                         df[, grepl("log2_R[0-9]{3}", names(df)), 
                            drop = FALSE]) %>%
    dplyr::mutate_at(.vars = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)),
                     list(~ replace(.x, is.infinite(.x), NA_real_)))
  
  if (TMT_plex) {
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
  
  df <- local({
    res <- df %>% 
      calcSD_Splex(id = group_pep_by, type = type_sd) %>% 
      `names<-`(gsub(paste0("^", type_sd), "sd_log2_R", names(.)))
    
    res %>% 
      dplyr::right_join(df, by = group_pep_by) %>% 
      na_zeroIntensity() %>% 
      dplyr::mutate(TMT_Set = set_idx)
  })

  invisible(df)
}


#' Helper of PSM2Pep
#' 
#' @param file The name of a PSM file.
#' @param ... filter_dots.
#' @inheritParams PSM2Pep
#' @inheritParams load_expts
#' @inheritParams n_TMT_sets
#' @inheritParams normPSM
#' @importFrom magrittr %>% %T>% %$% %<>% 
psm_to_pep <- function (file = NULL, dat_dir = NULL, label_scheme_full = NULL, 
                        group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        method_psm_pep = "median", lfq_ret_tol = 60L, 
                        rm_allna = FALSE, type_sd = "log2_R", ...) 
{
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  fn_prx <- gsub("_PSM_N.txt", "", file, fixed = TRUE)
  set_idx <- as.integer(gsub(".*TMTset(\\d+)_.*", "\\1", fn_prx))
  injn_idx <- as.integer(gsub(".*LCMSinj(\\d+).*", "\\1", fn_prx))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  df <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, "PSM", file), 
                    col_types = get_col_types(), 
                    show_col_types = FALSE)
  )

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::select(-grep("^sd_log2_R", names(.)))
  
  # special handling for MaxQuant; no `Precursor Intensity` for .d files
  is_mq_lfq <- if (find_search_engine(dat_dir) == "mq" && 
                   TMT_plex == 0L && 
                   all(is.nan(df[["I000"]]))) TRUE else FALSE
  
  if (is_mq_lfq) {
    message("Precursor intensity not available at PSM levels (.d files).\n",  
            "(LFQ) intensity from MaxQuant \"peptides[...].txt\" will be used.")
    
    df <- rm_cols_mqpsm(df, group_psm_by, set_idx)
  } 
  else {
    df <- calcPeptide(df = df, group_psm_by = group_psm_by, 
                      method_psm_pep = method_psm_pep, 
                      group_pep_by = group_pep_by, dat_dir = dat_dir, 
                      set_idx = set_idx, injn_idx = injn_idx, 
                      TMT_plex = TMT_plex, 
                      lfq_ret_tol = lfq_ret_tol, rm_allna = rm_allna, 
                      type_sd = type_sd)
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
    readr::write_tsv(file.path(dat_dir, "Peptide", paste0(fn_prx, "_Peptide_N.txt")))
}


#' Interim peptide tables
#' 
#' \code{PSM2Pep} summarizes
#' \href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSMs} to
# 'peptides by individual TMT (or LFQ) experiments and LC/MS series.
#' 
#' In general, fields other than \code{log2FC} and \code{intensity} are
#' summarized with median statistics. One exception is with \code{pep_expect} in
#' Mascot or \code{PEP} in MaxQuant where geometric mean is applied.
#' 
#' @param method_psm_pep Character string; the method to summarize the
#'   \code{log2FC} and the \code{intensity} of \code{PSMs} by peptide entries.
#'   The descriptive statistics includes \code{c("mean", "median",
#'   "weighted_mean", "top_3_mean", "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum",
#'   "lfq_all")} with \code{median} being the default for TMT and
#'   \code{lfq_top_2_sum} for LFQ. The \code{log10-intensity} of reporter (or
#'   LFQ) ions at the \code{PSMs} levels will be the weight when summarizing
#'   \code{log2FC} with various \code{"top_n"} statistics or
#'   \code{"weighted_mean"}.
#' @param lfq_ret_tol The tolerance of retention time (in seconds) for the
#'   aggregation of LFQ data.
#' @inheritParams normPSM
#' @inheritParams calcPeptide
#' @param ... \code{filter_}: Variable argument statements for the filtration of
#'   data rows. See also \code{\link{normPSM}}.
#' @seealso 
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
#'  \emph{Variable arguments of \code{filter_...}} \cr 
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
#' @family custom database preparation
#' @seealso 
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
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'  
#' @return Tables under \code{Peptide} folder for each TMT experiment and LC/MS
#'   series: \code{TMTset1_LCMSinj1_Peptide_N.txt} etc.
#'
#' @example inst/extdata/examples/PSM2Pep_.R
#' @import stringr dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @export
PSM2Pep <- function(method_psm_pep = c("median", "mean", "weighted_mean", 
                                       "top_3_mean", "lfq_max", "lfq_top_2_sum", 
                                       "lfq_top_3_sum", "lfq_all"), 
                    lfq_ret_tol = 60L, rm_allna = FALSE, 
                    type_sd = c("log2_R", "N_log2_R", "Z_log2_R"), ...) 
{
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  old_opts <- options()
  on.exit(options(old_opts), add = TRUE)
  options(warn = 1)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
          c(dots) %>% save_call("PSM2Pep")
      }
    }, 
    add = TRUE
  )
  
  dots <- rlang::enexprs(...)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  # ---
  method_psm_pep <- rlang::enexpr(method_psm_pep)
  
  if (TMT_plex) {
    if (length(method_psm_pep) > 1L) 
      method_psm_pep <- "median"
    else
      method_psm_pep <- rlang::as_string(method_psm_pep)
  } 
  else {
    if (length(method_psm_pep) > 1L)
      method_psm_pep <- "lfq_top_2_sum"
    else 
      method_psm_pep <- rlang::as_string(method_psm_pep)
  }
  
  if (method_psm_pep == "top.3") {
    stop("Method \"top.3\" depreciated; instead use \"top_3_mean\".", 
         call. = FALSE)
  } 
  else if (method_psm_pep == "weighted.mean") {
    stop("Method \"weighted.mean\" depreciated; instead use \"weighted_mean\".", 
         call. = FALSE)
  }
  
  stopifnot(method_psm_pep %in% c("median", "mean", "weighted_mean", "top_3_mean", 
                                  "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum", 
                                  "lfq_all"), 
            length(method_psm_pep) == 1L)
  
  # ---
  type_sd <- rlang::enexpr(type_sd)
  
  if (length(type_sd) > 1L) 
    type_sd <- "log2_R"
  else
    type_sd <- rlang::as_string(type_sd)
  
  stopifnot(type_sd %in% c("log2_R", "N_log2_R", "Z_log2_R"), 
            length(type_sd) == 1L)
  
  # ---
  stopifnot(vapply(c(lfq_ret_tol), is.numeric, logical(1)))
  stopifnot(vapply(c(rm_allna), rlang::is_logical, logical(1)))
  
  dir.create(file.path(dat_dir, "Peptide/cache"), 
             recursive = TRUE, showWarnings = FALSE)
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), 
                         pattern = "_PSM_N\\.txt$") %>%
    reorder_files()
  
  if (!length(filelist)) 
    stop("Files of \"_PSM_N.txt\" not found.")
  
  message("Primary column keys in \"PSM/TMTset1_LCMSinj1_PSM_N.txt\" etc. ", 
          "for \"filter_\" varargs.")
  
  n_files <- length(filelist)
  
  if (n_files > 1L) {
    n_cores <- min(parallel::detectCores(), n_files)
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(cl, list("add_cols_at", "reloc_col_before"), 
                            envir = environment(proteoQ:::add_cols_at))
    
    parallel::clusterExport(cl, list("exprs"), 
                            envir = environment(rlang::exprs))
    
    suppressWarnings(
      silent_out <- parallel::clusterApply(
        cl, filelist, psm_to_pep, 
        dat_dir = dat_dir, label_scheme_full = label_scheme_full, 
        group_psm_by = group_psm_by, group_pep_by = group_pep_by, 
        method_psm_pep = method_psm_pep, lfq_ret_tol = lfq_ret_tol, 
        rm_allna = rm_allna, type_sd = type_sd, ...)
    )
    
    parallel::stopCluster(cl)
  } 
  else {
    psm_to_pep(filelist, dat_dir = dat_dir, label_scheme_full = label_scheme_full, 
               group_psm_by = group_psm_by, group_pep_by = group_pep_by, 
               method_psm_pep = method_psm_pep, lfq_ret_tol = lfq_ret_tol, 
               rm_allna = rm_allna, type_sd = type_sd, ...)
  }
  
  .saveCall <- TRUE
  
  invisible(NULL)
}


#' To lower cases
#' 
#' @param x A character string of amino acid sequence.
#' @param ch A tag before the letter to conversion to its lower case.
my_tolower <- function(x = "", ch = "^") 
{
  locales <- gregexpr(ch, x) %>% .[[1]] %>% `+`(., 1)
  lowers <- purrr::map(locales, ~ substr(x, .x, .x)) %>% tolower()
  
  for (i in seq_along(lowers)) 
    substr(x, locales[i], locales[i]) <- lowers[i]

  gsub(ch, "", x)
}


#' Adds the \code{pep_seq_mod} field to MaxQuant PSMs
#'
#' @inheritParams splitPSM
#' @inheritParams locate_outliers
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_maxquant_pepseqmod <- function(df = NULL, use_lowercase_aa = TRUE) 
{
  if (!use_lowercase_aa) 
    return(df)
  
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
      dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, 
                          fixed = TRUE))
    
    rest <- df %>% 
      dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, 
                           fixed = TRUE))
    
    if (nrow(n_ac)) {
      n_ac <- n_ac %>% 
        dplyr::mutate(pep_seq_mod = 
                        gsub("(Acetyl (Protein N-term))", "_", pep_seq_mod, 
                             fixed = TRUE))
      df <- rbind(rest, n_ac)
    }
    
    df
  })
  
  # (2-2) add "_" to sequences from protein C-terminal amidation
  # ".AAASNGPVK(Xxx->Yyy)(Amidated (Protein C-term))." -> 
  # ".AAASNGPVK(Xxx->Yyy)_."
  
  df <- local({
    c_am <- df %>% 
      dplyr::filter(grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
    
    rest <- df %>% 
      dplyr::filter(!grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
    
    if (nrow(c_am)) {
      c_am <- c_am %>% 
        dplyr::mutate(pep_seq_mod = gsub("(Amidated (Protein C-term))", "_", 
                                         pep_seq_mod, fixed = TRUE))
      df <- rbind(rest, c_am)
    }
    
    df
  })
  
  # (3-1) "~" for "(Protein N-term)" other than acetylation
  # assume no dual (My (Protein N-term)) + "Acetyl": "._(My (Protein N-term))AAASSLTK."    
  # ".(My (Protein N-term))(Gln->pyro-Glu)QAAAAQGSnGPVK." -> 
  # ".~(Gln->pyro-Glu)QAAAAQGSnGPVK."
  
  df <- local({
    n_ac <- df %>% 
      dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
    
    other_n <- df %>% 
      dplyr::filter(grepl("Protein N-term", Modifications, fixed = TRUE)) %>% 
      dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
    
    rest <- df %>% 
      dplyr::filter(!grepl("Protein N-term", Modifications, fixed = TRUE))
    
    if (nrow(other_n)) {
      other_n <- other_n %>% 
        dplyr::mutate(pep_seq_mod = 
                        gsub("^\\.\\(.*\\s+?\\(+?.*\\){2}", ".~", pep_seq_mod))
      df <- rbind(rest, n_ac, other_n)
    }
    
    df
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
    
    if (nrow(other_c)) {
      other_c <- other_c %>% 
        dplyr::mutate(pep_seq_mod = 
                        gsub("^(.*)\\(.*\\s+?\\(+?.*\\){2}\\.$", paste0("\\1", "~."),
                             pep_seq_mod))
      df <- rbind(rest, c_am, other_c)
    }
    
    df
  })
  
  # (4-1) "^" for peptide "(N-term)" modification
  # "._(Carbamyl (N-term))AAAAGALAPGPLPDLAAR." -> "._^AAAAGALAPGPLPDLAAR."
  
  df <- local({
    nt <- df %>% 
      dplyr::filter(grepl("(N-term)", Modifications, fixed = TRUE))
    
    rest <- df %>% 
      dplyr::filter(!grepl("(N-term)", Modifications, fixed = TRUE))
    
    if (nrow(nt)) {
      nt <- nt %>% 
        dplyr::mutate(pep_seq_mod = gsub("(^\\.[_~]{0,1})\\(.*\\s+?\\(+?.*\\){2}", 
                                         paste0("\\1", "^"), pep_seq_mod))
      df <- rbind(rest, nt)
    }
    
    df
  })
  
  # (4-2) "^" for peptide "(C-term)" modification
  # ".AAAANLCPGQDR(My (C-term))_." -> ".AAAANLCPGQDR^_."
  
  df <- local({
    ct <- df %>% 
      dplyr::filter(grepl("(C-term)", Modifications, fixed = TRUE))
    
    rest <- df %>% 
      dplyr::filter(!grepl("(C-term)", Modifications, fixed = TRUE))
    
    if (nrow(ct)) {
      ct <- ct %>% 
        dplyr::mutate(pep_seq_mod = gsub("^(.*)\\(.*\\s+?\\(+?.*\\){2}([_~]{0,1}\\.$)", 
                                         paste0("\\1", "^", "\\2"), pep_seq_mod)) 
      df <- rbind(rest, ct)
    }
    
    df
  })
  
  # (5) remove "." at both ends
  
  df <- df %>% 
    dplyr::mutate(pep_seq_mod = gsub("\\.", "", pep_seq_mod))
  
  # (6) other N- or C-terminal modifications better 
  # but not named with "N-term" or "C-term": 
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
                         pep_seq_mod))
}


#' Adds the \code{pep_seq_mod} field to MSFragger PSMs
#'
#' @inheritParams splitPSM
#' @inheritParams locate_outliers
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_msfragger_pepseqmod <- function(df = NULL, use_lowercase_aa = TRUE) 
{
  if (!use_lowercase_aa) return(df)
  
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
}


#' Pads columns to a placeholder data frame.
#' 
#' @param df The original data frame.
#' @param df2 The data frame to be inserted.
#' @param idx The index of \code{df} column for \code{df2} to be inserted
#'   (after).
add_cols_at <- function(df = NULL, df2 = NULL, idx = 0L) 
{
  stopifnot(idx >= 0L)
  
  if (idx == 0L) {
    bf <- NULL
  } 
  else {
    bf <- df[, seq_len(idx), drop = FALSE]
  }
  
  if ((idx + 1) <= ncol(df)) {
    af <- df[, (idx + 1) : ncol(df), drop = FALSE]
  } 
  else {
    af <- NULL
  }
  
  dplyr::bind_cols(
    bf,
    df2,
    af,
  )
}


#' Replaces columns in the original PSM table.
#' 
#' The column index(es) need to be continuous.
#' 
#' @param df The original data frame.
#' @param df2 The data columns to replace those in \code{df}.
#' @param idxs The sequences of column indexes in \code{df}. Note that
#'   \code{idxs} need to be a continuous sequences.
replace_cols_at <- function(df = NULL, df2 = NULL, idxs = 1L) 
{
  ncol <- ncol(df)
  stopifnot(all(idxs >= 1L), all(idxs <= ncol))
  
  idxs <- sort(idxs)
  stopifnot(all.equal(idxs - idxs[1] + 1, seq_along(idxs)))
  
  if (idxs[1] >= 2L) {
    bf <- df[, 1:(idxs[1]-1), drop = FALSE]
  } 
  else {
    bf <- NULL
  }
  
  if (idxs[length(idxs)] < ncol(df)) {
    af <- df[, (idxs[length(idxs)]+1):ncol(df), drop = FALSE]
  } 
  else {
    af <- NULL
  }
  
  dplyr::bind_cols(
    bf,
    df2,
    af
  )
}


#' Relocates column "m/z" to be immediately before column "Mass".
#' 
#' Not currently used.
#' 
#' @param df The original data frame.
#' @param from The column to be moved from.
#' @param to The column to which the \code{from} will be moved before.
reloc_col <- function (df = NULL, from = "m/z", to = "Mass") 
{
  df0 <- df
  df2 <- suppressWarnings(df %>% dplyr::select(one_of(from)))
  
  if (!ncol(df2)) 
    return(df0)
  
  df <- df %>% dplyr::select(-one_of(from))
  
  idx <- which(names(df) == to) - 1
  
  if (!length(idx)) 
    return(df0)
  
  df <- add_cols_at(df, df2, idx)
  
  return(df)
}


#' Relocates column "to_move" immediately after column "col_before".
#'
#' @param df The original data frame.
#' @param to_move The column to be moved.
#' @param col_before The anchor column to which the \code{to_move} will be moved
#'   after.
reloc_col_after <- function(df = NULL, to_move = "after_anchor", 
                            col_before = "anchor") 
{
  if (!(to_move %in% names(df) && col_before %in% names(df))) return(df)
  
  if (to_move == col_before) 
    return(df)
  
  df2 <- df %>% dplyr::select(one_of(to_move))
  df <- df %>% dplyr::select(-one_of(to_move))
  
  idx <- which(names(df) == col_before)
  
  add_cols_at(df, df2, idx)
}


#' Relocates column "to_move" immediately after the last column.
#' 
#' @inheritParams reloc_col_after
reloc_col_after_last <- function (df = NULL, to_move = "after_anchor") 
{
  col_last <- names(df)[ncol(df)]
  reloc_col_after(df, to_move, col_last)
}


#' Relocates column "to_move" immediately after the first column.
#' 
#' @inheritParams reloc_col_after
reloc_col_after_first <- function(df = NULL, to_move = "after_anchor") 
{
  col_first <- names(df)[1]
  reloc_col_after(df, to_move, col_first)
}


#' Relocates column "to_move" immediately before anchor column "col_after".
#'
#' The same as \code{reloc_col}.
#'
#' @param df The original data frame.
#' @param to_move The column to be moved.
#' @param col_after The anchor column to which the \code{to_move} will be moved
#'   before.
reloc_col_before <- function(df = NULL, to_move = "before_anchor", 
                             col_after = "anchor") 
{
  if (!(to_move %in% names(df) && col_after %in% names(df))) 
    return(df)

  df2 <- df %>% dplyr::select(one_of(to_move))
  df <- df %>% dplyr::select(-one_of(to_move))
  
  idx <- which(names(df) == col_after)
  
  add_cols_at(df, df2, idx - 1)
}


#' Relocates column "to_move" immediately before the last column.
#' 
#' @inheritParams reloc_col_after
reloc_col_before_last <- function(df = NULL, to_move = "after_anchor") 
{
  col_last <- names(df)[ncol(df)]
  reloc_col_before(df, to_move, col_last)
}


#' Relocates column "to_move" immediately before the first column.
#' 
#' @inheritParams reloc_col_after
reloc_col_before_first <- function (df = NULL, to_move = "after_anchor") 
{
  col_first <- names(df)[1]
  reloc_col_before(df, to_move, col_first)
}


#' Helper: finds the column name before \code{to_move}.
#' 
#' To keep columns at the same order after descriptive summary.
#' 
#' @inheritParams reloc_col_after
find_preceding_colnm <- function(df = NULL, to_move = NULL) 
{
  if (!to_move %in% names(df)) {
    stop("Column ", to_move, " not found.", 
         call. = FALSE)
  }
  
  ind_bf <- which(names(df) == to_move) - 1
  
  if (ind_bf == 0) {
    names(df)[1]
  } 
  else {
    names(df)[ind_bf]
  }
}


#' Orders PSM columns.
#'
#' Only for certain columns.
#' @param df A PSM data frame.
#' @param cols A vector; the names of columns to be moved to the front of
#'   \code{df}.
order_psm_cols <- function(df = NULL, cols = NULL) 
{
  purrr::walk(cols, ~ {
    if (is.null(df[[.x]])) df[[.x]] <- NA
    df <<- df
  }, df)
  
  dplyr::bind_cols(
    df %>% dplyr::select(-cols), 
    df %>% dplyr::select(cols),
  )
}


#' Orders Mascot PSM columns
#' 
#' Columns indicated by \code{psm_cols} will be moved to the front, and columns
#' of ratios and ratios and intensities to the end. 
#' 
#' @param df A data frame of PSM.
#' @param psm_cols A character string of column names.
#' @param rm_na_cols Logical; if TRUE, remove columns of all NAs.
order_mascot_psm_cols <- function(df, psm_cols = NULL, rm_na_cols = FALSE) 
{
  if (is.null(psm_cols)) {
    psm_cols <- c("prot_hit_num", "prot_family_member", "prot_acc", 
                  "prot_desc", "prot_score", "prot_mass", 
                  "prot_matches",	"prot_matches_sig",	
                  "prot_sequences", "prot_sequences_sig",	
                  "prot_n_psm", "prot_n_pep", 
                  "prot_len",	"prot_pi", "prot_tax_str", "prot_tax_id", 
                  "prot_seq",	"prot_empai",	"prot_icover",	"prot_cover",
                  "prot_index",	
                  
                  "pep_query",	"pep_rank",	"pep_n_psm", "pep_isbold",	
                  "pep_isunique",	"pep_literal_unique",	"pep_razor_unique",	
                  "pep_tot_int", "pep_unique_int",	"pep_razor_int", 
                  "pep_exp_mz",	"pep_exp_mr",	"pep_exp_z",	"pep_calc_mr", 
                  "pep_delta",	"pep_score", "pep_homol",	"pep_ident",
                  "pep_expect",	"pep_res_before", "pep_seq",
                  "pep_seq_mod",	"pep_res_after",	"pep_start", 
                  "pep_end",	"pep_len",	"pep_miss",	"pep_istryptic",	
                  "pep_frame", "pep_var_mod",	"pep_var_mod_pos", 
                  "pep_summed_mod_pos",	"pep_local_mod_pos",
                  "pep_num_match",	"pep_scan_title",	"pep_index",	
                  "pep_scan_range",	"pep_ret_range",
                  "pep_ms2_sumint",	"pep_n_ions",	"pep_locprob",	"pep_locdiff",
                  "pep_phospho_locprob", "pep_phospho_locdiff", 
                  "pep_ions_first", "pep_ions_second", "pep_ions_third", 
                  
                  "gene",	"fasta_name",	"uniprot_acc",	"uniprot_id",
                  "refseq_acc",	"other_acc",	"entrez",	"species",
                  "acc_type",	"shared_prot_accs",	"shared_genes",	"kin_attr",
                  "kin_class",	"kin_order",	

                  "dat_file",	"psm_index")
  }
  
  df <- df %>% 
    order_psm_cols(psm_cols) 
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(psm_cols), 
    df %>% dplyr::select(-psm_cols), 
  )
    
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^prot_|^pep_", names(.))),
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("^[IR]{1}[0-9]{3}[NC]{0,1}", names(.))),
    df %>% dplyr::select(grep("^[IR]{1}[0-9]{3}[NC]{0,1}", names(.))),
  )
  
  df <- df %>% 
    { if (rm_na_cols) dplyr::select(., not_all_NA(.)) else . }
}


#' Pads MaxQuant TMT channels to the highest plex. 
#' 
#' @param ... filter_dots.
#' @inheritParams splitPSM
#' @inheritParams splitPSM_mq
#' @inheritParams pad_mascot_channels
pad_mq_channels <- function(file = NULL, fasta = NULL, entrez = NULL, 
                            corrected_int = TRUE, ...) 
{
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dat_dir <- get_gl_dat_dir()
  base_name <- gsub("\\.txt$", "", file)
  
  # (QE, timsTOF: different cases in column keys)
  df <- read.csv(file.path(dat_dir, file), 
                 check.names = FALSE, header = TRUE, 
                 sep = "\t", comment.char = "#") %>% 
    filters_in_call(!!!filter_dots)

  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  
  # interim solution; wait till the column keys become consistent with MQ
  df <- df %>% 
    { if ("Precursor intensity" %in% names(.)) 
      dplyr::rename(., `Precursor Intensity` = `Precursor intensity`) else . } %>% 
    { if ("Precursor apex fraction" %in% names(.)) 
      dplyr::rename(., `Precursor Apex Fraction` = `Precursor apex fraction`) else . }
  
  if (! "Precursor Intensity" %in% names(df)) {
    stop("Column \"Precursor Intensity\" not found.\n", 
         "May need to change the case: \"Precursor intensity -> Precursor Intensity\".",
         call. = FALSE)
  } 
  else if (! "Precursor Apex Fraction" %in% names(df)) {
    stop("Column \"Precursor Apex Fraction\" not found.\n", 
         "May need to change the case: \"Precursor apex fraction -> Precursor Apex Fraction\".",
         call. = FALSE)
  } 
  else {
    df <- df %>% 
      dplyr::mutate(`Precursor Intensity` = `Precursor Intensity`/`Precursor Apex Fraction`)
  }

  # --- A patch for inconsistency in MaxQuant msms.txt columns ---
  # TMT: no `Gene Names` and `Protein Names` columns for RefSeq etc.
  # LFQ: 
  #  (1) UniProt identifiers: with `Gene Names` and `Protein Names` columns
  #  (2) RefSeq identifiers: without `Gene Names` and `Protein Names` columns

  if (!all(c("Gene Names", "Protein Names") %in% stringr::str_to_title(names(df))) ) {
    stopifnot("Proteins" %in% names(df))
    
    df$"Gene names" <- df$"Gene Names" <- df$"Protein names" <- df$"Protein Names" <- NULL

    # some proteins were not annotated...
    df <- local({
      df <- df %>% 
        dplyr::mutate(Proteins = gsub("\\.[0-9]*", "", Proteins)) %>% 
        dplyr::mutate(prot_acc = gsub(";.*$", "", Proteins))
      
      old_opts <- options()
      options(warn = -1)
      
      tempdata <- suppressMessages(
        df %>% 
          dplyr::select(prot_acc) %>% 
          dplyr::filter(!duplicated(prot_acc)) %>% 
          annotPrn(fasta, entrez) %>% 
          dplyr::select(prot_acc, gene, prot_desc) %>% 
          dplyr::rename(`Gene Names` = "gene", "Protein Names" = "prot_desc")
      )
      
      options(old_opts)
      
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
  
  TMT_plex <- TMT_plex(label_scheme_full)
  
  if (TMT_plex == 0L) {
    # no "Corrected" intensity for LFQ
    df <- df %>% 
      dplyr::mutate(dat_file = base_name, 
                    I000 = `Precursor Intensity`)

    return(df)
  }
  
  ## TMT only (no LFQ) from this point on
  
  if (! "PSM_File" %in% names(fraction_scheme)) {
    
    # raw_files and tmt_sets may be later used for error messages
    
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
    
    label_scheme_sub <- label_scheme_full %>% 
      dplyr::filter(TMT_Set %in% tmt_sets)
  } 
  else {
    raw_files <- NULL
    
    tmt_sets <- fraction_scheme %>% 
      dplyr::mutate(PSM_File = gsub("\\.txt$", "", PSM_File)) %>% 
      dplyr::filter(PSM_File == base_name, !duplicated(TMT_Set)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(TMT_Set) %>% 
      unlist()
    
    label_scheme_sub <- label_scheme_full %>% 
      dplyr::filter(TMT_Set %in% tmt_sets)
  }

  nas <- data.frame(rep(NA, nrow(df)))
  sample_ids <- as.character(label_scheme_sub$Sample_ID)

  str_int1 <- "^Reporter Intensity [0-9]+"
  str_int2 <- "^Reporter Intensity Corrected [0-9]+"
  str_dev <- "^Reporter Mass Deviation \\[Mda\\] [0-9]+"

  df_int <- df %>% dplyr::select(grep(str_int1, stringr::str_to_title(names(.))))
  df_int2 <- df %>% dplyr::select(grep(str_int2, stringr::str_to_title(names(.))))
  df_dev <- df %>% dplyr::select(grep(str_dev, stringr::str_to_title(names(.))))
  
  this_plex <- ncol(df_int)
  TMT_plex <- TMT_plex(label_scheme_full)
  stopifnot(this_plex <= TMT_plex, this_plex >= 0)
  
  # Empty.xxx can be due to either channel padding or removals
  if ((this_plex > 0) && (this_plex < TMT_plex)) {
    pos <- find_padding_pos(this_plex, TMT_plex)

    for (idx in seq_along(pos)) {
      df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))
      df_int2 <- suppressMessages(add_cols_at(df_int2, nas, pos[idx] - 1))
      df_dev <- suppressMessages(add_cols_at(df_dev, nas, pos[idx] - 1))
    }
    
    rm(list = c("idx"))

    len <- length(sample_ids)
    
    if (ncol(df_int) == len && ncol(df_int2) == len && ncol(df_dev) == len) {
      names(df_int) <- paste("Reporter intensity", seq_len(len))
      df <- replace_cols_at(df, df_int, 
                            grep(str_int1, stringr::str_to_title(names(df))))
      
      names(df_int2) <- paste("Reporter intensity corrected", seq_len(len))
      df <- replace_cols_at(df, df_int2, 
                            grep(str_int2, stringr::str_to_title(names(df))))
      
      names(df_dev) <- paste("Reporter mass deviation [mDa]", seq_len(len))
      df <- replace_cols_at(df, df_dev, 
                            grep(str_dev, stringr::str_to_title(names(df))))
    }
  }

  # add I126 etc. and remove `Reporter Intensity ...`
  
  if (!corrected_int) {
    df <- df %>% 
      dplyr::bind_cols(df_int %>% `names<-`(find_int_cols(TMT_plex)))
  } 
  else {
    df <- df %>% 
      dplyr::bind_cols(df_int2 %>% `names<-`(find_int_cols(TMT_plex)))
  }

  df <- df %>% 
    dplyr::select(-grep("^Reporter Intensity [0-9]+", 
                        stringr::str_to_title(names(.)))) %>% 
    dplyr::select(-grep("^Reporter Intensity Corrected [0-9]+", 
                        stringr::str_to_title(names(.))))

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
#'are calculated later in \code{annotPSM}. This is suitable mainly because there
#'is no columns like \code{prot_matches_sig} and \code{prot_sequences_sig} need
#'to be updated after data merging in \code{splitPSM_mq}.
#'
#'@param corrected_int A logical argument for uses with \code{MaxQuant} TMT. At
#'  the TRUE default, values under columns "Reporter intensity corrected..." in
#'  \code{MaxQuant} PSM results (\code{msms.txt}) will be used. Otherwise,
#'  "Reporter intensity" values without corrections will be used.
#'@param rm_reverses A logical argument for uses with \code{MaxQuant} TMT and
#'  LFQ. At the TRUE default, \code{Reverse} entries will be removed.
#'@inheritParams splitPSM
#'@import dplyr tidyr stringr
#'@importFrom magrittr %>% %T>% %$% %<>% equals
splitPSM_mq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        fasta = NULL, entrez = NULL, 
                        pep_unique_by = "group", 
                        scale_rptr_int = FALSE, corrected_int = TRUE, 
                        rm_craps = FALSE, rm_krts = FALSE, rm_allna = FALSE, 
                        rm_reverses = TRUE, purge_phosphodata = TRUE, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, 
                        rptr_intco = 0, rptr_intrange = c(0, 100), 
                        use_lowercase_aa = TRUE, 
                        parallel = TRUE, ...) 
{
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        message("Split PSMs by TMT experiments and LCMS series --- Completed.")
      }
    }, 
    add = TRUE
  )

  old_opts <- options()
  on.exit(options(old_opts), add = TRUE)
  options(warn = 1, warning.length = 5000L)
  
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$")
  
  if (!length(filelist)) {
    stop("No PSM files of \"msms[...].txt\" under", file.path(dat_dir), ".",
         "\nMake sure that the names of PSM files start with \"msms\".", 
         call. = FALSE)
  }
  
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in \"msms[...].txt\" for \"filter_\" varargs.")
  
  # (1.1) row filtration, column padding and psm file combinations
  # (make also `RAW_File`, I000 or I126 etc.; 
  # not yet `RAW_File`: need to first change MaxQuant names to the title case)
  
  df <- purrr::map(filelist, pad_mq_channels, fasta, entrez, 
                   corrected_int, !!!filter_dots)
  df <- suppressWarnings(dplyr::bind_rows(df))
  
  if (!"dat_file" %in% names(df)) 
    stop("Column \"dat_file\" not found.")
  
  if (!any(grepl("^I[0-9]{3}[NC]{0,1}", names(df))))
    stop("Intensity columns not found.")

  # (1.2.1) handle inconsistent column keys 
  # (after `filter_dots`)
  df <- local({
    nms_old <- names(df)

    df_int <- df %>% 
      dplyr::select(grep("^I[0-9]{3}[NC]{0,1}$", names(.)))
    
    df <- df %>% 
      dplyr::select(-grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %>% 
      `names<-`(stringr::str_to_title(names(.))) %>% 
      dplyr::rename(dat_file = Dat_file, 
                    pep_expect = Pep, 
                    pep_exp_mz = `M/Z`, 
                    pep_seq = Sequence, 
                    pep_score = Score,
                    RAW_File = `Raw File`) %>% 
      dplyr::bind_cols(df_int)

    data.frame(original_name = nms_old, title_case_name = names(df)) %>% 
      readr::write_excel_csv(file.path(dat_dir, "PSM/cache/psm_colkeys_lookup.csv"))
    
    warning("\n================================================================\n", 
            "Column keys in \"msms[...].txt\" ", 
            "may be changed to a \"Title Case\":\n",
            "( ) \"Protein Names\" <-> \"Protein Names\"... \n", 
            "(x) \"Protein names\" -> \"Protein Names\"... \n", 
            "(Lookups in \"~/PSM/cache/MQ_colkey_lookup.csv\".)\n", 
            "================================================================\n",
            call. = FALSE)
    
    invisible(df)
  })
  
  stopifnot("RAW_File" %in% names(df))
  
  # (1.2.2) add ratio columns
  if (TMT_plex && (sum(grepl("^I[0-9]{3}[NC]{0,1}", names(df))) > 1L)) {
    df <- sweep(df[, find_int_cols(TMT_plex), drop = FALSE], 
                1, df[["I126"]], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), 
                       ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
  } 
  else {
    # (1) `Modified sequence` not yet available in MaxQuant `peptides.txt` 
    #     for matching with those in `msms.ttx`
    # (2) `Precursor intensity` not available in `msms.txt`
    
    local({
      if (all(is.nan(df[["I000"]]))) {
        if (group_psm_by == "pep_seq_mod") {
          stop("Currently with MaxQuant timsTOF, only \"group_psm_by = pep_seq\" ", 
               "to match queries in `peptides.txt`.", 
               call. = FALSE)
        }
        
        filelist <- list.files(path = file.path(dat_dir), 
                               pattern = "^peptides.*\\.txt$")
        
        if (!file.exists(file.path(dat_dir, filelist))) {
          stop("All NA values under \"Precursor Intensity\" in \"msms[...].txt\".\n", 
               "To proceed, provide \"peptides[...].txt\" for intensity back-filling.", 
               call. = FALSE)
        }
      }
    })
    
    df <- df %>% 
      reloc_col_after_last("I000") %>% 
      dplyr::mutate(R000 = I000/I000, 
                    R000 = ifelse(is.infinite(R000), NA_real_, R000)) 
  }
  
  # (1.3) clean up prot_acc
  # zero-char exception in `Proteins` and not necessarily "Reverse" entries
  df <- df %>% 
    dplyr::filter(as.character(.[["Proteins"]]) > 0) %>% 
    { if (rm_craps) dplyr::filter(., !grepl("^CON_", .[["Proteins"]])) else . } %>% 
    { if (rm_reverses) dplyr::filter(., is.na(Reverse) || Reverse != "+") else . } %>% 
    dplyr::mutate(prot_acc = gsub("\\;.*", "", Proteins)) %>% 
    dplyr::filter(prot_acc != "") %>% 
    { if (rm_craps) dplyr::filter(., !grepl("\\|.*\\|$", prot_acc)) else . } 

  # (1.4.1) add pep_seq_mod
  stopifnot("Modified Sequence" %in% names(df))
  
  df <- df %>% 
    dplyr::rename(pep_seq_mod = `Modified Sequence`) %>%
    add_maxquant_pepseqmod(use_lowercase_aa) 
  
  # (1.4.2) phospho
  df <- local({
    nms_phospho <- stringr::str_to_title("Phospho (STY) Probabilities")
    
    phos_idx <- grep(nms_phospho, names(df), fixed = TRUE)
    
    if (length(phos_idx) >= 2L) 
      phos_idx <- phos_idx[1]
    
    if (length(phos_idx) >= 1L) {
      if (purge_phosphodata) {
        df <- df %>% dplyr::filter(!nchar(as.character(.[[nms_phospho]])) == 0) 
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
    
    df
  })
  
  # (2.1a) annotate proteins
  df <- df %>% 
    dplyr::mutate(prot_acc = gsub("\\.[0-9]*$", "", prot_acc)) %>% 
    annotPrn(fasta, entrez) %>%  
    { if (!"gene" %in% names(.)) dplyr::mutate(., gene = prot_acc) else .} %>% 
    dplyr::mutate(gene = ifelse(is.na(gene), prot_acc, gene)) %>% 
    dplyr::select(-which(names(.) %in% c("Gene Names", "Gene names", 
                                         "Protein Names", "Protein names")))
  
  # (2.1b) add `pep_res_before` and `pep_res_after`
  df <- df %>% 
    annotPeppos() %>% 
    reloc_col_before("pep_seq", "pep_res_after") %>% 
    reloc_col_before("pep_res_before", "pep_seq")

  # (2.2) compile "preferred" columns
  df <- local({
    df <- df %>% 
      { if("Scan Number" %in% names(.)) . 
        else dplyr::mutate(., `Scan Number` = NA) } %>% 
      { if("Retention Time" %in% names(.)) . 
        else dplyr::mutate(., `Retention Time` = NA) } %>% 
      { if("Intensities" %in% names(.)) . 
        else dplyr::mutate(., `Intensities` = NA) } %>% 
      { if("Number Of Matches" %in% names(.)) . 
        else dplyr::mutate(., `Number Of Matches` = NA) } %>% 
      { if("Localization Prob" %in% names(.)) . 
        else dplyr::mutate(., `Localization Prob` = NA) } 
    
    df <- df %>% 
      dplyr::mutate(pep_scan_range = `Scan Number`,
                    pep_ret_range = `Retention Time`, 
                    pep_ms2_sumint = .[, "Intensities"] %>% 
                      stringr::str_split(., pattern = ";") %>% 
                      purrr::map_dbl(~ .x %>% as.numeric() %>% sum(na.rm = TRUE)) ,
                    pep_n_ions = `Number Of Matches`,
                    pep_locprob = `Localization Prob`, 
                    pep_locdiff = NA) %>% 
      # `Scan Number` and `Retention Time` kept as anchors 
      # for column removals during `psm_to_pep`
      dplyr::select(-c("Number Of Matches", "Localization Prob"))
    
    if ("pep_index" %in% names(df)) df$pep_index <- NULL
    if ("prot_index" %in% names(df)) df$prot_index <- NULL
    df <- df %>% 
      add_entry_ids("pep_seq", "pep_index") %>% 
      add_entry_ids("prot_acc", "prot_index")
  })

  # (2.4) find the shared prot_accs and genes for each peptide
  message("\nParsing shared proteins...\n")
  
  df <- df %>% 
    dplyr::mutate(Proteins = gsub("\\.[0-9]*", "", Proteins)) %>% 
    add_shared_genes(key = "Proteins", sep = ";", fasta, entrez)

  # (2.5) find the uniqueness of peptides
  df <- local({
    stopifnot(all(c("Protein Group Ids", "Proteins") %in% names(df)))
    
    # Not possible (provided one prot_acc one gene):
    # pep1 | acc1  | gn1, gn2
    # 
    # Unique by protein  <=> unique by gene
    # df %>% filter(!grepl(";", `Proteins`)) %>% filter(grepl(",", `shared_genes`))
    # df %>% filter(grepl(",", `shared_genes`)) %>% filter(!grepl(";", `Proteins`))
    # 
    # NOT unique by prot_acc group -> can still be unique by gene group
    # acc_1, acc_2 group are in C3 group 
    # 
    # Unique by prot_acc group  -> unique by gene group
    # (by definition, gn2 must be razor as no evidence of unique presence of acc2 
    # -> no evidence of gn2)
    # (1) grp1 | acc1, (razor acc2) | gn1 
    # (2) grp1 | acc1, (razor acc2) | gn1 gn2
    # 
    # Could simply use `Protein Group Ids` as indicators of razor genes. 
    # However, there are a few discrepant examples in MQ 
    # With SINGLE msms.txt file: 
    # shouldn't be any below, but there are a few MQ cases of one gene multiple groups 
    # gn_grp_maps %>% dplyr::filter(!grepl(";", gene_groups), grepl(", ", gene_groups))
    # gene   grpid_gene
    # POLR2M 4094
    # POLR2M 1286
    
    ## new calculations (not yet for multiple msms.txt): 
    # likely to have different `Protein Group Ids` 
    #   for the same gene at different msms.txt files
    # MAP4K3  3431;5721, 3415, 3392, 3617
    # Elf2 9951;13287, 10127;13596
    
    if (FALSE) {
      gn_grp_maps <- local({
        uniq_by <- c("Protein Group Ids", "gene")
        
        gn_grp_maps <- df %>% 
          dplyr::select(uniq_by) %>% 
          tidyr::unite(uniq_id, uniq_by, sep = "@", remove = FALSE) %>% 
          dplyr::filter(!duplicated(uniq_id)) %>% 
          dplyr::select(-uniq_id) %>% 
          dplyr::group_by(gene) 
        
        # keep the first longest
        gn_grp_maps <- gn_grp_maps %>% 
          dplyr::mutate(n_commas = stringr::str_count(`Protein Group Ids`, ";")) %>% 
          dplyr::arrange(-n_commas) %>% 
          dplyr::filter(row_number() == 1)
        
        gn_grp_maps <- gn_grp_maps %>% 
          dplyr::summarise(gene_groups = paste(`Protein Group Ids`, collapse  = ", ")) %>% 
          dplyr::mutate(gene_groups = gsub(";", ", ", gene_groups))
        
        # non-redundant map
        # (would be redundant at multiple msms.txt(s))
        gn_grp_maps <- gn_grp_maps[["gene_groups"]] %>% 
          stringr::str_split(", ") %>% 
          list_to_dataframe() %>% 
          dplyr::bind_cols(
            gn_grp_maps %>% dplyr::select("gene"), 
          ) %>% 
          tidyr::gather("id.", "gene_groups", -gene) %>% 
          dplyr::select(-id.) %>% 
          dplyr::filter(!is.na(gene_groups)) %>% 
          tidyr::unite(gn_grp., gene, gene_groups, sep = "@", remove = FALSE) %>% 
          dplyr::filter(!duplicated(gn_grp.)) %>% 
          dplyr::select(-gn_grp.)
      })
      
      # ??? redundancy
      df <- df %>% 
        dplyr::left_join(gn_grp_maps, by = "gene") 
    }
    
    if (FALSE) {
      df <- local({
        df_grpids <- df %>% 
          split(.$dat_file, drop = TRUE) %>% 
          purrr::imap(~ {
            # length(unique(.x$gene)); 7932
            
            # MQ: the same gene in the same dat_file, 
            #   indexes in group_ids can still not unique
            # ABLIM1  457;5779
            # ABLIM1  457;5779;1649
            x_s <- .x %>% 
              dplyr::select(gene, `Protein Group Ids`) %>% 
              dplyr::filter(grepl(";", `Protein Group Ids`)) %>% 
              tidyr::unite(uniq_id, c("gene", "Protein Group Ids"), sep = "@", 
                           remove = FALSE) %>% 
              dplyr::filter(!duplicated(uniq_id)) %>% 
              dplyr::select(-uniq_id) %>% 
              dplyr::group_by(gene) %>% 
              dplyr::filter(row_number() == 1)
            # length(unique(x_s$gene)); 2465
            
            # MQ: the same gene in the same dat_file, can be both unique and shared
            # THADA 9717               
            # THADA 9717;13107
            
            # perhaps in MQ, group IDs calculated based on invidual RAW files

            x_u <- .x %>% 
              dplyr::select(gene, `Protein Group Ids`) %>% 
              dplyr::filter(!grepl(";", `Protein Group Ids`)) %>% 
              dplyr::filter(! .data$gene %in% x_s$gene) %>% 
              tidyr::unite(uniq_id, c("gene", "Protein Group Ids"), sep = "@", 
                           remove = FALSE) %>% 
              dplyr::filter(!duplicated(uniq_id)) %>% 
              dplyr::select(-uniq_id) %>% 
              dplyr::group_by(gene) %>% 
              dplyr::filter(row_number() == 1)
            # length(unique(x_u$gene)); 5467
            
            x <- dplyr::bind_rows(x_u, x_s) %>% 
              dplyr::mutate(dat_file = .y)
            
            # given that x_s calculated earlier than x_u;
            # a gene found shared once is considered shared.
          }) %>% 
          dplyr::bind_rows() %>% 
          tidyr::unite(gn_dat., c("gene", "dat_file"), sep = "@", 
                       remove = TRUE)
        
        # updata `Protein Group Ids`
        df <- df %>% 
          tidyr::unite(gn_dat., c("gene", "dat_file"), sep = "@", 
                       remove = FALSE) %>% 
          dplyr::select(-`Protein Group Ids`) %>% 
          dplyr::left_join(df_grpids, by = "gn_dat.") %>% 
          dplyr::select(-gn_dat.)
      })
    }
    
    
    df <- df %>% 
      dplyr::mutate(prot_acc_groups = gsub(";", ", ", `Protein Group Ids`)) %>% 
      dplyr::mutate(gene_groups = prot_acc_groups)
    
    if (group_pep_by == "prot_acc") {
      df <- df %>% 
        dplyr::mutate(pep_literal_unique = ifelse(grepl(", ", shared_prot_accs), 
                                                  FALSE, TRUE)) %>% 
        dplyr::mutate(pep_razor_unique = ifelse(grepl(", ", prot_acc_groups), 
                                                FALSE, TRUE))
    } 
    else if (group_pep_by == "gene") {
      df <- df %>% 
        dplyr::mutate(pep_literal_unique = ifelse(grepl(", ", shared_genes), 
                                                  FALSE, TRUE)) %>% 
        dplyr::mutate(pep_razor_unique = ifelse(grepl(", ", gene_groups), 
                                                FALSE, TRUE))
    } 
    else {
      stop("\"group_pep_by\" is not one of \"prot_acc\" or \"gene\".", 
           call. = FALSE)
    }
    
    df <- df %>% 
      dplyr::mutate(pep_razor_unique = ifelse(pep_literal_unique, 
                                              pep_literal_unique, 
                                              pep_razor_unique))

    # (keeps the original columns for examination) 
    df <- df %>% 
      dplyr::select(-c("prot_acc_groups", "gene_groups"))
    
    if (pep_unique_by == "group") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    } 
    else if (pep_unique_by == "protein") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_literal_unique)
    } 
    else if (pep_unique_by == "none") {
      df <- df %>% dplyr::mutate(pep_isunique = TRUE)
    } 
    else {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    }
    
    df <- df %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
    
    # An unknown MQ case: "single protein but multiple groups"
    #   unique at protein but not at group
    # Sequence	Proteins	Gene Names	Protein Group Ids	
    # EAEDSLRR	NP_004530	NARS1	4870;13627
  })
  
  df <- local({
    colnm_int <- "Precursor Intensity"
    
    df %>% 
      dplyr::mutate(pep_tot_int = !!rlang::sym(colnm_int)) %>% 
      dplyr::mutate(pep_unique_int = 
                      ifelse(pep_literal_unique, pep_tot_int, 0)) %>% 
      dplyr::mutate(pep_razor_int = 
                      ifelse(pep_razor_unique, pep_tot_int, 0)) %>% 
      dplyr::select(-!!rlang::sym(colnm_int)) %>% 
      reloc_col_after("pep_unique_int", "pep_tot_int") %>% 
      reloc_col_after("pep_razor_int", "pep_unique_int")
  })
  
  stopifnot(all(c("pep_isunique", "pep_literal_unique", "pep_razor_unique", 
                  "pep_tot_int", "pep_unique_int", "pep_razor_int") %in% 
                  names(df)))
  
  # (2.6) add peptide properties and prot_cover, prot_icover
  df <- df %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq), 
                                    stringr::str_count(pep_seq, "[KR]") - 1,
                                    stringr::str_count(pep_seq, "[KR]"))) %>% 
    add_prot_icover(id = group_pep_by) %>% 
    { if (!("prot_cover" %in% names(.) && length(filelist) == 1)) 
        calc_cover(., id = !!rlang::sym(group_pep_by)) 
      else . } %>% 
    dplyr::select(-which(names(.) %in% c("Length", "Missed cleavages", "Missed Clevages"))) 
  
  # (2.3) add columns pep_n_psm, prot_n_psm, prot_n_pep
  df <- df %>% add_quality_cols(!!group_psm_by, !!group_pep_by)
  
  .saveCall <- TRUE
  
  invisible(df)
}


#' Pads Spectrum Mill TMT channels to the highest plex
#' 
#' @param ... filter_dots.
#' @inheritParams pad_mascot_channels
pad_sm_channels <- function(file = NULL, ...) 
{
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dat_dir <- get_gl_dat_dir()
  
  base_name <- file %>% gsub("\\.ssv$", "", .)
  
  df <- suppressWarnings(
    readr::read_delim(file.path(dat_dir, file), delim = ";", 
                      col_types = cols(filename = col_character()))
  ) %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate_at(.vars = which(names(.) %in% c("accession_number", 
                                                   "accession_numbers")), 
                     ~ gsub("\\.[0-9]$", "", .x))
  
  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  
  if (TMT_plex == 0L) {
    df <- df %>% 
      dplyr::mutate(dat_file = base_name, 
                    I000 = `totalIntensity`)
    
    return(df)
  }
  
  if (! "PSM_File" %in% names(fraction_scheme)) {
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
    
    label_scheme_sub <- label_scheme_full %>% 
      dplyr::filter(TMT_Set %in% tmt_sets)
  } else {
    raw_files <- NULL
    
    tmt_sets <- fraction_scheme %>% 
      dplyr::mutate(PSM_File = gsub("\\.ssv$", "", PSM_File)) %>% 
      dplyr::filter(PSM_File == base_name, !duplicated(TMT_Set)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(TMT_Set) %>% 
      unlist()
    
    label_scheme_sub <- label_scheme_full %>% 
      dplyr::filter(TMT_Set %in% tmt_sets)
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
  
  stopifnot(length(ref) == 1L)
  
  this_plex <- ncol(df_int)
  TMT_plex <- TMT_plex(label_scheme_full)
  stopifnot(this_plex <= TMT_plex, this_plex >= 0L)
  
  if ((this_plex > 0L) && (this_plex < TMT_plex)) {
    if (TMT_plex == 18) {
      keys_int <- paste0("TMT_", c("126", "127N", "127C", "128N", "128C", "129N", 
                                   "129C", "130N", "130C", "131N", "131C", "132N", 
                                   "132C", "133N", "133C", "134N", 
                                   "134C", "135N"))
    } 
    else if (TMT_plex == 16) {
      keys_int <- paste0("TMT_", c("126", "127N", "127C", "128N", "128C", "129N", 
                                   "129C", "130N", "130C", "131N", "131C", "132N", 
                                   "132C", "133N", "133C", "134N"))
    } 
    else if (TMT_plex == 11) {
      keys_int <- paste0("TMT_", c("126", "127N", "127C", "128N", "128C", "129N", 
                                   "129C", "130N", "130C", "131N", "131C"))
    } 
    else if (TMT_plex == 10) {
      keys_int <- paste0("TMT_", c("126", "127N", "127C", "128N", "128C", "129N", 
                                   "129C", "130N", "130C", "131"))
    } 
    else if(TMT_plex == 6) {
      keys_int <- paste0("TMT_", c("126", "127", "128", "129", "130", "131"))
    } 
    else {
      keys_int <- NULL
    }
    
    keys_ratio <- paste0(keys_int, "_", ref) %>% 
      .[!grepl(paste0("_", ref, "_", ref), .)]
    
    pos <- find_padding_pos(this_plex, TMT_plex)

    for (idx in seq_along(pos)) {
      df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))
    }
    
    rm(list = "idx")

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
  
  df <- df %>% 
    dplyr::bind_cols(df_int %>% `names<-`(find_int_cols(TMT_plex)))
  
  df <- df %>% 
    dplyr::select(-grep("^TMT_[0-9]{3}[NC]{0,1}$", names(.))) %>% 
    dplyr::select(-grep("^TMT_[0-9]{3}[NC]{0,1}_[0-9]{3}[NC]{0,1}$", names(.)))

  df$dat_file <- base_name
  
  return(df)
}


#' Adds the \code{pep_seq_mod} field to MaxQuant PSMs
#'
#' @inheritParams splitPSM
#' @inheritParams locate_outliers
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_sm_pepseqmod <- function(df = NULL, use_lowercase_aa = TRUE) 
{
  if (!use_lowercase_aa) return(df)
  
  if (! "pep_start" %in% names(df)) {
    stop("Need \"pep_start\" for \"pep_seq_mod\".", call. = FALSE)
  }
  
  # (1) mods under column `variableSites`
  df <- local({
    if (length(which(names(df) == "variableSites"))) {
      df_sub <- df %>% dplyr::filter(!is.na(variableSites))
      
      if (nrow(df_sub)) {
        pos_matrix <- df_sub %>% 
          dplyr::select(variableSites) %>% 
          dplyr::mutate(variableSites = as.character(variableSites)) %$% 
          stringr::str_split(.$variableSites, " ") %>% 
          list_to_dataframe() %>% 
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
    
    df
  })
  
  # (2-1) protein n-term acetylation
  if (length(which(names(df) == "nterm"))) {
    df_sub <- df %>% dplyr::filter(nterm == "Acetyl")
    
    if (nrow(df_sub)) {
      df_sub$pep_seq_mod <- paste0("_", df_sub$pep_seq_mod)
      
      df <- dplyr::bind_rows(
        df %>% dplyr::filter(nterm != "Acetyl"), 
        df_sub
      )
    }        
  }
  
  # (2-2) protein C-terminal amidation
  if (length(which(names(df) == "cterm"))) {
    # hypothetical; no yet known how it will be named in SM
    df_sub <- df %>% dplyr::filter(grepl("^Amidate", cterm))
    
    if (nrow(df_sub)) {
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
  
  invisible(df)
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
                        scale_rptr_int = FALSE, 
                        rm_craps = FALSE, rm_krts = FALSE, rm_allna = FALSE, 
                        purge_phosphodata = TRUE, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, 
                        rptr_intco = 0, rptr_intrange = c(0, 100), 
                        use_lowercase_aa = TRUE, parallel = TRUE, ...) 
{
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        message("Split PSMs by TMT experiments and LCMS series --- Completed.")
      }
    }, 
    add = TRUE
  )
  
  if (is.null(fasta)) 
    stop("FASTA file(s) not provided.", call. = FALSE)
  
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$")
  
  if (!length(filelist)) 
    stop("No PSM files were found under", file.path(dat_dir), 
         "\nCheck that the names of PSM files start with \"PSMexport\".", 
         call. = FALSE)
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in \"PSMexport[...].ssv\" for \"filter_\" varargs.")
  
  # (1.1) row filtration, column padding and psm file combinations
  # (make also `RAW_File`, I000 or I126 etc.)
  df <- purrr::map(filelist, pad_sm_channels, !!!filter_dots)
  df <- suppressWarnings(dplyr::bind_rows(df))
  
  stopifnot(c("sequence", "accession_numbers") %in% names(df))
  stopifnot("dat_file" %in% names(df))
  stopifnot(any(grepl("^I[0-9]{3}[NC]{0,1}", names(df))))

  nms_sm <- names(df) %>% .[! . %in% c("dat_file")]
  
  # (1.2.1) handle inconsistent column keys 
  # ...

  # (1.2.2) add ratio columns
  if ((TMT_plex > 0L) && (sum(grepl("^I[0-9]{3}[NC]{0,1}", names(df))) > 1L)) {
    df <- sweep(df[, find_int_cols(TMT_plex), drop = FALSE], 
                1, df[["I126"]], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), 
                       ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
  } 
  else {
    df <- df %>% 
      reloc_col_after_last("I000") %>% 
      dplyr::mutate(R000 = I000/I000, 
                    R000 = ifelse(is.infinite(R000), NA, R000)) 
  }
  
  # (1.3) clean up prot_acc
  # (no craps removal by grepl("\\|.*\\|$", accession_number); ending "|" not present)
  df <- df %>% 
    dplyr::mutate(pep_seq = toupper(sequence), 
                  prot_acc = accession_number) %>% 
    dplyr::filter(prot_acc != "") %>% 
    { if (rm_craps && "Protein" %in% names(.)) 
      dplyr::filter(., !grepl("\\|.*\\|$", Protein)) else . } 
  
  # (2.1a) annotate proteins
  df <- df %>% 
    annotPrn(fasta, entrez) %>%  
    { if (!"gene" %in% names(.)) dplyr::mutate(., gene = prot_acc) else .} %>% 
    dplyr::mutate(gene = ifelse(is.na(gene), prot_acc, gene))
  
  # (2.1b) add pep_res_before and pep_res_after
  # (needed before adding pep_seq_mod)
  df <- df %>% 
    annotPeppos() %>% 
    reloc_col_before("pep_seq", "pep_res_after") %>% 
    reloc_col_before("pep_res_before", "pep_seq")
  
  # (1.4.1) add pep_seq_mod
  # (need `pep_start`)
  df <- df %>%
    dplyr::mutate(pep_seq_mod = as.character(pep_seq)) %>% 
    add_sm_pepseqmod(use_lowercase_aa)
  
  # (1.4.2) phospho
  df <- local({
    is_phospho_expt <- any(grepl("Phosphorylated", df$modifications))
    
    if (is_phospho_expt && purge_phosphodata) {
      if ("modifications" %in% names(df)) {
        df <- df %>% dplyr::filter(grepl("Phosphorylated", modifications))
      } 
      else {
        warning("Missing column \"modifications\" for non-phosphopeptide removals.", 
                call. = FALSE)
      }
    }
    
    df
  })
  
  # (2.2) compile "preferred" columns
  df <- local({
    df <- df %>% 
      dplyr::mutate(pep_scan_range = NA,
                    pep_ret_range = NA, 
                    pep_ms2_sumint = NA,
                    pep_n_ions = NA,
                    pep_locprob = NA, 
                    pep_locdiff = NA) 
    
    df <- df %>% 
      dplyr::rename(RAW_File = `filename`) %>% 
      dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", RAW_File)) %>% 
      dplyr::mutate(pep_expect = NA) %>% 
      dplyr::mutate(pep_isunique = NA) %>% 
      dplyr::rename(pep_score = score) 

    if ("pep_index" %in% names(df)) df$pep_index <- NULL
    if ("prot_index" %in% names(df)) df$prot_index <- NULL
    
    df <- df %>% 
      add_entry_ids("pep_seq", "pep_index") %>% 
      add_entry_ids("prot_acc", "prot_index")
  })
  
  # (2.4) find the shared prot_accs and genes for each peptide
  df <- df %>% 
    add_shared_sm_genes(key = "accession_numbers", sep = "\\|", fasta, entrez)
  
  # (2.5) find the uniqueness of peptides
  df <- local({
    warning("No protein grouping yet available; ", 
            "all peptide sequences are treated as razor unique.", 
            call. = FALSE)
    
    df <- df %>% 
      dplyr::mutate(pep_literal_unique = ifelse(grepl(", ", shared_prot_accs), 
                                                FALSE, TRUE), 
                    pep_razor_unique = TRUE)
    
    if (pep_unique_by == "group") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    } 
    else if (pep_unique_by == "protein") {
      df <- df %>% dplyr::mutate(pep_isunique = pep_literal_unique)
    } 
    else if (pep_unique_by == "none") {
      df <- df %>% dplyr::mutate(pep_isunique = TRUE)
    } 
    else {
      df <- df %>% dplyr::mutate(pep_isunique = pep_razor_unique)
    }
    
    df %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
  })
  
  df <- df %>% 
    dplyr::mutate(pep_tot_int = NA, 
                  pep_scan_range = NA, 
                  pep_ret_range = NA, 
                  pep_ms2_sumint = NA,
                  pep_n_ions = NA, 
                  pep_locprob = NA,
                  pep_locdiff = NA) %>% 
    dplyr::mutate(pep_tot_int = totalIntensity) %>% 
    dplyr::mutate(pep_unique_int = 
                    ifelse(pep_literal_unique, pep_tot_int, 0)) %>% 
    dplyr::mutate(pep_razor_int = 
                    ifelse(pep_razor_unique, pep_tot_int, 0)) %>% 
    reloc_col_after("pep_unique_int", "pep_tot_int") %>% 
    reloc_col_after("pep_razor_int", "pep_unique_int")
  
  stopifnot(all(c("pep_isunique", "pep_literal_unique", "pep_razor_unique", 
                  "pep_tot_int", "pep_unique_int", "pep_razor_int") %in% 
                  names(df)))
  
  # (2.6) add peptide properties and prot_cover, prot_icover
  # M._sequence.c; -._sequence.c; n.sequence.c; -.sequence.c
  df <- df %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq), 
                                    stringr::str_count(pep_seq, "[KR]") - 1,
                                    stringr::str_count(pep_seq, "[KR]"))) %>% 
    add_prot_icover(id = group_pep_by) %>% 
    { if (!("prot_cover" %in% names(.) && length(filelist) == 1)) 
        calc_cover(., id = !!rlang::sym(group_pep_by)) 
      else . }
  
  # (2.3) add columns pep_n_psm, prot_n_psm, prot_n_pep
  df <- df %>% add_quality_cols(!!group_psm_by, !!group_pep_by)
  
  .saveCall <- TRUE
  
  invisible(df)
}


#' Pads MSFragger TMT channels to the highest plex. 
#' 
#' @param ... filter_dots.
#' @inheritParams pad_mascot_channels
pad_mf_channels <- function(file = NULL, ...) 
{
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]

  dat_dir <- get_gl_dat_dir()
  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  base_name <- gsub("\\.tsv$", "", file)
  
  df <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, file), 
                    col_types = cols(`Is Unique` = col_logical()))
  ) %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate_at(.vars = which(names(.) %in% c("Protein", "Protein ID",	
                                                   "Entry Name", "Mapped Genes", 
                                                   "Mapped Proteins")), 
                     ~ gsub("\\.[0-9]$", "", .x))

  nms_df <- names(df)
  
  if (!"Spectrum" %in% nms_df)
    stop("Column 'Spectrum' not found.")

  # !!! "Intensity" in both TMT and LFQ
  # "Purity" only with TMT
  if ("Purity" %in% nms_df) {
    df_int <- local({
      df_int <- df[(which(nms_df == "Purity") + 1):ncol(df)]
      
      nums <- purrr::map_lgl(df_int, is.numeric)
      csum <- cumsum(nums)
      cols <- diff(csum, 1)
      
      df_int[, names(cols[cols == 1]), drop = FALSE]
    })
  } 
  else if ("Intensity" %in% nms_df) {
    return(dplyr::mutate(df, dat_file = base_name, I000 = Intensity))
  } 
  else {
    stop("Neither column 'Purity' or 'Intensity' were found.")
  }
  
  rm(list = "nms_df")
  
  # (Assumed MSFragger still uses MGF from MSConvert)
  df <- df %>% 
    dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", Spectrum))

  ## TMT only (no LFQ) from this point on
  
  if (!"PSM_File" %in% names(fraction_scheme)) {
    # raw_files and tmt_sets may be later used for error messages
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
    
    label_scheme_sub <- label_scheme_full %>% 
      dplyr::filter(TMT_Set %in% tmt_sets)
  } 
  else {
    raw_files <- NULL
    
    tmt_sets <- fraction_scheme %>% 
      dplyr::mutate(PSM_File = gsub("\\.tsv$", "", PSM_File)) %>% 
      dplyr::filter(PSM_File == base_name, !duplicated(TMT_Set)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(TMT_Set) %>% 
      unlist()
    
    label_scheme_sub <- label_scheme_full %>% 
      dplyr::filter(TMT_Set %in% tmt_sets)
  }
  
  sample_ids <- as.character(label_scheme_sub$Sample_ID)
  this_plex <- ncol(df_int)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  stopifnot(this_plex <= TMT_plex, this_plex >= 0L)

  # Empty.xxx can be due to either channel padding or removals
  
  if (this_plex && (this_plex < TMT_plex)) {
    pos <- find_padding_pos(this_plex, TMT_plex)
    nas <- data.frame(rep(NA_real_, nrow(df)))
    
    for (idx in seq_along(pos))
      df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1))

    rm(list = c("idx", "nas"))
  }
  
  # new column names
  df <- local({
    n_ints <- ncol(df_int)
    n_samples <- length(sample_ids)
    
    stopifnot(n_ints >= 1L, n_samples >= n_ints)
    
    if (n_ints < n_samples) {
      if (!((n_samples %% n_ints) == 0)) {
        stop("Number of intensity columns: ", n_ints, 
             " is not a multiple of number of samples: ", n_samples, "\n",
             call. = FALSE)
      } 
      else {
        warning("In TMT plex(es) ", 
                paste(tmt_sets, collapse = ", "), ": \n", 
                "number of intensity columns: ", n_ints, "\n",
                "number of samples: ", n_samples, ".\n",
                "Assume a merged search at folds: ", n_samples/n_ints, ".", 
                call. = FALSE)
      }
    }
    
    nms_old <- names(df_int)
    names(df_int) <- find_int_cols(TMT_plex)
    
    # OK that sample IDs from MSFragger may be named differently to 
    # those in label_scheme

    df <- df %>% 
      dplyr::select(-which(names(.) %in% nms_old), 
                    -which(names(.) %in% names(df_int))) %>% 
      dplyr::bind_cols(df_int)
  })

  # the same `raw_file` may be at different `dat_file`s
  
  df$dat_file <- base_name
  
  invisible(df)
}


#'Splits MSFragger PSM tables
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
                        scale_rptr_int = FALSE, 
                        rm_craps = FALSE, rm_krts = FALSE, rm_allna = FALSE, 
                        purge_phosphodata = TRUE, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, 
                        rptr_intco = 0, rptr_intrange = c(0, 100), 
                        use_lowercase_aa = TRUE, 
                        parallel = TRUE, ...) 
{
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        message("Split PSMs by TMT experiments and LCMS series --- Completed.")
      }
    }, 
    add = TRUE
  )
  
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
  # assume it is always 3 for now (as by MSConvert).
  
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^psm.*\\.tsv$")
  
  if (!length(filelist)) {
    stop("No PSM files of \"psm[...].tsv\" under", file.path(dat_dir), ".",
         "\nMake sure that the names of PSM files start with \"psm\".", 
         call. = FALSE)
  }
  
  message("Primary column keys in \"psm[...].tsv\" for \"filter_\" varargs.")
  
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dots <- dots %>% .[! . %in% filter_dots]

  # (1.1) row filtration, column padding and psm file combinations
  # (make also `RAW_File`, I000 or I126 etc.)
  df <- purrr::map(filelist, pad_mf_channels, !!!filter_dots)
  df <- suppressWarnings(dplyr::bind_rows(df))
  
  stopifnot(c("Is Unique", # razor uniqueness
              "Mapped Proteins", # literal and razor uniqueness
              "Intensity", # precursor intensity
              "Modified Peptide" # pep_seq_mod
  ) %in% names(df))
  
  # (1.2.1) handle inconsistent column keys 
  # (Not applicable)
  
  # (1.2.2) add ratio columns
  if (TMT_plex && (sum(grepl("^I[0-9]{3}[NC]{0,1}", names(df))) > 1)) {
    df <- sweep(df[, find_int_cols(TMT_plex), drop = FALSE], 
                1, df[["I126"]], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), 
                       ~ replace(.x, is.infinite(.x), NA_real_)) %>% 
      dplyr::bind_cols(df, .) 
  } 
  else {
    df <- df %>% 
      reloc_col_after_last("I000") %>% 
      dplyr::mutate(R000 = I000/I000, 
                    R000 = ifelse(is.infinite(R000), NA_real_, R000)) 
  }
  
  # (1.3) clean up prot_acc
  # (no craps removal by grepl("\\|.*\\|$", Protein); ending "|" not present)
  df <- df %>% 
    dplyr::rename(pep_seq = Peptide, 
                  prot_acc = `Protein ID`, 
                  ) %>% 
    dplyr::filter(prot_acc != "") %>% 
    { if (rm_craps && "Protein" %in% names(.)) 
      dplyr::filter(., !grepl("\\|.*\\|$", Protein)) else . } 
  
  # (2.1a) annotate proteins
  df <- df %>% 
    annotPrn(fasta, entrez) %>%  
    { if (!"gene" %in% names(.)) dplyr::mutate(., gene = prot_acc) else .} %>% 
    dplyr::mutate(gene = ifelse(is.na(gene), prot_acc, gene))

  # (2.1b) add pep_res_before and pep_res_after
  # (needed before adding pep_seq_mod)
  df <- df %>% 
    annotPeppos() %>% 
    reloc_col_before("pep_seq", "pep_res_after") %>% 
    reloc_col_before("pep_res_before", "pep_seq")
  
  # (1.4.1) add pep_seq_mod
  df <- df %>% 
    dplyr::mutate(`Modified Peptide` = ifelse(is.na(`Modified Peptide`), 
                                              pep_seq, `Modified Peptide`)) %>% 
    dplyr::rename(pep_seq_mod = `Modified Peptide`) %>% 
    add_msfragger_pepseqmod(use_lowercase_aa) 
  
  # (1.4.2) phospho
  df <- local({
    # assumed use_lowercase_aa = TRUE; otherwise on lower-case sty?
    is_phospho_expt <- any(grepl("[sty]", df$pep_seq_mod))
    
    if (is_phospho_expt && purge_phosphodata)
      df <- dplyr::filter(df, grepl("[sty]", pep_seq_mod))
    
    if ("STY:79.9663" %in% names(df)) {
      prs <- stringr::str_extract_all(df[["STY:79.9663"]], "\\([^()]+\\)")
      prs <- lapply(prs, function (x) as.numeric(substring(x, 2L, nchar(x) - 1L)))
      prs <- lapply(prs, sort, decreasing = TRUE)
      
      df$pep_phospho_locprob <- unlist(lapply(prs, `[`, 1), recursive = FALSE)
      
      df$pep_phospho_locdiff <- lapply(prs, function (x) {
        if (length(x) <= 1L) 0 else x[1] - x[2]
      }) %>% 
        unlist(recursive = FALSE)
      
      rm(list = "prs")
    }

    invisible(df)
  })
  
  # (2.2) compile "preferred" columns
  df <- local({
    df <- df %>% 
      { if("Retention" %in% names(.)) . else dplyr::mutate(., Retention = NA_real_) } %>% 
      dplyr::mutate(pep_scan_num = NA_character_,
                    pep_ret_range = Retention, 
                    pep_ms2_sumint = NA_real_,
                    pep_n_ions = NA_integer_,
                    pep_locprob = NA_real_, 
                    pep_locdiff = NA_real_) %>% 
      dplyr::select(-which(names(.) == "Retention"))
    
    df <- df %>% 
      mutate(pep_scan_num = gsub("[^\\.]+?\\.(.*)", "\\1", Spectrum)) %>% 
      mutate(pep_scan_num = gsub("([^\\.]+?)\\..*", "\\1", pep_scan_num)) %>% 
      mutate(pep_scan_num = gsub("^0+", "", pep_scan_num)) # may be leading zeros
    
    df <- df %>% 
      dplyr::mutate(
        pep_expect = Expectation, 
        pep_score = Hyperscore, 
        
        pep_exp_mz = `Observed M/Z`, 
        pep_exp_mr = `Observed Mass`, 
        pep_exp_z = Charge,
        pep_n_exp_z = NA_integer_, 
        pep_calc_mr = `Calculated M/Z`, 
        pep_delta = `Delta Mass`, 
      ) %>% 
      add_n_pepexpz(group_psm_by)

    if ("pep_index" %in% names(df)) 
      df$pep_index <- NULL
    
    if ("prot_index" %in% names(df)) 
      df$prot_index <- NULL
    
    df <- df %>% 
      add_entry_ids("pep_seq", "pep_index") %>% 
      add_entry_ids("prot_acc", "prot_index")
  })
  
  # (2.4) find the shared prot_accs and genes for each peptide
  # (the original uniqueness of peptides by MSFragger may not holder after 
  # the joining of PSM files)
  df <- add_shared_prot_accs_mf(df)

  if (all(is.na(df$Gene))) {
    df <- add_shared_genes(df, key = "shared_prot_accs", sep = ", ", fasta, entrez)
  } 
  else {
    df <- dplyr::mutate(df, shared_genes = 
                          ifelse(is.na(`Mapped Genes`), 
                                 gene, 
                                 paste(gene, `Mapped Genes`, sep = ", ")))
  }

  # (2.5) find the uniqueness of peptides
  df <- local({
    if (group_pep_by == "prot_acc") {
      df <- df %>% 
        dplyr::mutate(pep_literal_unique = 
                        ifelse(`Is Unique` & !grepl(", ", shared_prot_accs), TRUE, FALSE), 
                      pep_razor_unique = `Is Unique`)
    } 
    else if (group_pep_by == "gene") {
      df <- df %>% 
        dplyr::mutate(pep_literal_unique = 
                        ifelse(`Is Unique` & !grepl(", ", shared_genes), TRUE, FALSE), 
                      pep_razor_unique = `Is Unique`)
    } 
    else {
      stop("\"group_pep_by\" is not one of \"prot_acc\" or \"gene\".")
    }
    
    if (pep_unique_by == "group") {
      df <- dplyr::mutate(df, pep_isunique = pep_razor_unique)
    } 
    else if (pep_unique_by == "protein") {
      df <- dplyr::mutate(df, pep_isunique = pep_literal_unique)
    } 
    else if (pep_unique_by == "none") {
      df <- dplyr::mutate(df, pep_isunique = TRUE)
    } 
    else {
      df <- dplyr::mutate(df, pep_isunique = pep_razor_unique)
    }
    
    df <- df %>% 
      reloc_col_after("pep_literal_unique", "pep_isunique") %>% 
      reloc_col_after("pep_razor_unique", "pep_literal_unique")
  })

  # proteins:
  # (1) prot_tot_int: Gene == Ckap5 | Mapped Genes == Ckap5
  # (2) prot_razor_int: Gene == Ckap5
  # (3) prot_uniq_int: Gene == Ckap5 & void Mapped Genes
  # top3 for each
  
  df <- df %>% 
    dplyr::mutate(pep_tot_int = Intensity) %>% 
    dplyr::mutate(pep_unique_int = ifelse(pep_literal_unique, Intensity, 0)) %>% 
    dplyr::mutate(pep_razor_int = ifelse(pep_razor_unique, Intensity, 0)) %>% 
    dplyr::select(-Intensity)

  # (2.6) add peptide properties and prot_cover, prot_icover
  df <- df %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq), 
                                    stringr::str_count(pep_seq, "[KR]") - 1,
                                    stringr::str_count(pep_seq, "[KR]"))) %>% 
    add_prot_icover(id = group_pep_by) %>% 
    { if (!("prot_cover" %in% names(.) && length(filelist) == 1L)) 
        calc_cover(., id = !!rlang::sym(group_pep_by)) 
      else . } 
  
  # (2.3) add columns pep_n_psm, prot_n_psm, prot_n_pep
  df <- add_quality_cols(df, !!group_psm_by, !!group_pep_by)
  
  .saveCall <- TRUE
  
  invisible(df)
}


#' Finds the padding positions of TMT channels.
#'
#' @param this_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc. before
#'   padding.
#' @param Numeric; the maximum multiplexity of TMT indicated in metadata of
#'   \code{label_scheme.xlsx}.
find_padding_pos <- function (this_plex = 10L, TMT_plex = 10L) 
{
  if (!((this_plex > 0L) && (this_plex < TMT_plex)))
    return(NULL)

  if (this_plex == 6L) {
    if (TMT_plex == 10L)
      pos <- c(3, 5, 7, 9)
    else if (TMT_plex == 11L) 
      pos <- c(3, 5, 7, 9, 11)
    else if (TMT_plex == 16L)
      pos <- c(3, 5, 7, 9, 11:16)
    else if (TMT_plex == 18L)
      pos <- c(3, 5, 7, 9, 11:18)
    else
      stop("TMT_plex is not one of c(10, 11, 16, 18).", call. = FALSE)
  } 
  else if (this_plex == 10L) {
    if (TMT_plex == 11L)
      pos <- 11
    else if (TMT_plex == 16L)
      pos <- c(11:16)
    else if (TMT_plex == 18L)
      pos <- c(11:18)
    else
      stop("TMT_plex is not one of c(11, 16, 18).", call. = FALSE)
  } 
  else if (this_plex == 11L) {
    if (TMT_plex == 16L)
      pos <- c(12:16)
    else if (TMT_plex == 18L)
      pos <- c(12:18)
    else
      stop("TMT_plex is not one of c(16, 18).", call. = FALSE)
  } 
  else if (this_plex == 16L) {
    if (TMT_plex == 18L)
      pos <- c(17:18)
    else
      stop("TMT_plex is not one of c(18).", call. = FALSE)
  }
  else {
    stop("TMT_plex is not one of c(6, 10, 11, 16, 18).", call. = FALSE)
  }
  
  invisible(pos)
}


#' Pads TMT channels to the highest multiplex.
#' 
#' For proteoM.
#' 
#' @param file An intermediate PSM table.
#' @param ... filter_dots.
pad_tmt_channels <- function(file = NULL, ...) 
{
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dat_dir <- get_gl_dat_dir()
  load(file.path(dat_dir, "label_scheme_full.rda"))
  load(file.path(dat_dir, "fraction_scheme.rda"))
  
  base_name <- gsub("\\.txt$", "", file)
  
  df <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, file), 
                    col_types = cols(
                      prot_acc = col_character(), 
                      prot_issig = col_logical(), 
                      prot_isess = col_logical(),
                      prot_tier = col_integer(), 
                      prot_hit_num = col_integer(), 
                      prot_family_member = col_integer(), 
                      prot_es = col_number(), 
                      prot_es_co = col_number(), 
                      pep_seq = col_character(), 
                      pep_n_ms2 = col_integer(), 
                      pep_scan_title = col_character(), 
                      pep_exp_mz = col_number(),
                      pep_exp_mr = col_number(), 
                      pep_exp_z = col_character(), 
                      pep_calc_mr = col_number(), 
                      pep_delta = col_number(),
                      pep_tot_int = col_number(), 
                      pep_ret_range = col_number(), 
                      pep_scan_num = col_character(), # timsTOF
                      pep_mod_group = col_integer(), 
                      pep_frame = col_integer(), 
                      pep_fmod = col_character(),
                      pep_vmod = col_character(),
                      pep_isdecoy = col_logical(),
                      pep_ivmod = col_character(),
                      pep_len = col_integer(), 
                      pep_issig = col_logical(),
                      pep_score = col_double(),
                      pep_score_co = col_double(),
                      pep_rank = col_integer(), 
                      pep_locprob = col_double(),
                      pep_locdiff = col_double(),
                      pep_rank_nl = col_integer(), 
                      pep_literal_unique = col_logical(),
                      pep_razor_unique = col_logical(),
                      raw_file = col_character(), ), 
                    show_col_types = FALSE)
  ) %>% 
    filters_in_call(!!!filter_dots)

  if ("raw_file" %in% names(df)) {
    df <- df %>% 
      dplyr::rename(RAW_File = raw_file) %>% 
      dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
      dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File))
    
    df <- df %>%
      dplyr::mutate(pep_scan_title = gsub("\\\\", "/", pep_scan_title))
  } 
  else {
    df <- df %>%
      dplyr::mutate(pep_scan_title = gsub("\\\\", "/", pep_scan_title)) %>% 
      dplyr::mutate(RAW_File = pep_scan_title) %>% 
      dplyr::mutate(RAW_File = gsub("^.*/([^/]*)\\.raw[\\\"]{0,1}; .*", "\\1", 
                                    RAW_File)) %>% 
      dplyr::mutate(RAW_File = gsub("^.*/([^/]*)\\.d[\\\"]{0,1}; .*", "\\1", 
                                    RAW_File))
  }
  
  # (with some refseq_acc)
  df <- df %>%
    dplyr::mutate(prot_acc = gsub("\\.[0-9]*$", "", prot_acc)) 
  
  this_plex <- sum(grepl("^I[0-9]{3}[NC]{0,1}$", names(df)))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  stopifnot(this_plex >= 0L)
  
  if (this_plex == 0L) {
    return(dplyr::mutate(df, dat_file = base_name, I000 = pep_tot_int))
  }
  
  ## TMT only (not LFQ) from this point on
  
  if (this_plex > TMT_plex) {
    stop("\nPSM multiplexity is \"", this_plex, "\" with ", file, ".\n",
         "Metadata multiplexity is however smaller at \"", TMT_plex, "\".\n", 
         "Don't know which channels to exclude from PSM file.", 
         call. = FALSE)
  }
  
  # (Empty.xxx can be due to either channel padding or removals)
  if ((this_plex > 0L) && (this_plex < TMT_plex)) {
    df <- local({
      df_int <- df[, grepl("^I[0-9]{3}[NC]{0,1}$", names(df))]
      pos <- find_padding_pos(this_plex, TMT_plex)
      nas <- data.frame(rep(NA, nrow(df)))
      
      for (idx in seq_along(pos)) {
        df_int <- suppressMessages(add_cols_at(df_int, nas, pos[idx] - 1L))
      }
      
      names(df_int) <- find_int_cols(TMT_plex)

      df <- dplyr::bind_cols(
        df %>% select(-grep("^I[0-9]{3}[NC]{0,1}$", names(.))), 
        df_int)
    })
  } 

  df$dat_file <- base_name
  
  invisible(df)
}


#' Pads PSM exports.
#'
#' Intensity or ratio columns are handled during channel padding and thus
#' excluded.
#'
#' @param dfs Intermediate PSM tables.
pad_psm_fields <- function(dfs = NULL) 
{
  ncols <- purrr::map_dbl(dfs, ncol)
  
  nms <- lapply(dfs, function (x) {
    names(x) %>% 
      .[!grepl("[IR]{1}[0-9]{3}[NC]{0,1}$", .)] %>% 
      .[. != "dat_file"]
  })
  
  nms_union <- purrr::reduce(nms, union)
  
  nms_union <- c(nms_union[nms_union != "pep_scan_title"], 
                 nms_union[nms_union == "pep_scan_title"])
  
  snms_union <- sort(nms_union)
  
  ok_nms <- purrr::map_lgl(nms, function (x) {
    identical(sort(x), snms_union)
  })
  
  if (!all(ok_nms)) {
    warning("Inequal numbers or names of columns deteted: \n", 
            paste(ncols, collapse = "\n"), 
            call. = FALSE)
    
    for (i in seq_along(dfs)) {
      nms_i <- nms[[i]]
      df_i <- dfs[[i]]
      
      missing_nms <- setdiff(nms_union, nms_i)
      
      if (length(missing_nms)) {
        for (nm in missing_nms) df_i[[nm]] <- NA
      } 
      
      dfs[[i]] <- dplyr::bind_cols(
        df_i[nms_union], 
        df_i %>% 
          .[, ! names(.) %in% nms_union, drop = FALSE] %>% 
          .[, ! names(.) == "dat_file", drop = FALSE], 
        df_i["dat_file"], )
    }
  }
  
  invisible(dfs)
}


#' Adds the \code{pep_seq_mod} field to PSMs.
#' 
#' For proteoM.
#' 
#' @inheritParams locate_outliers
#' @inheritParams splitPSM
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>% %T>% %$% %<>% 
add_pepseqmod <- function(df, use_lowercase_aa = TRUE, purge_phosphodata = TRUE) 
{
  col_nms <- names(df)
  
  lapply(c("pep_seq", "pep_ivmod"), function (x) {
    if (! x %in% col_nms)
      stop("Column \"", x, "\" not found.", call. = FALSE)
  })
  
  if ("pep_seq_mod" %in% col_nms) {
    warning("Recompile column \"pep_seq_mod\".")
  }

  df <- df %>% dplyr::mutate(pep_seq_mod = pep_seq)
  
  fmods <- unique(df$pep_fmod)
  
  vmods <- unique(df$pep_vmod) %>% 
    .[! . == ""] %>% 
    .[!is.na(.)]

  df <- local({
    pat <- "Phospho\\s{1}\\("
    
    is_phospho_expt <- any(grepl(pat, c(fmods, vmods)))
    
    if (is_phospho_expt && purge_phosphodata) {
      df <- df %>% dplyr::filter(grepl(pat, pep_vmod) | grepl(pat, pep_fmod))
    }
    
    df
  })
  
  if (!length(vmods)) {
    warning("No applicable variable modifications", call. = FALSE)
    return(df)
  }
  
  if (use_lowercase_aa) {
    # (1) non terminal modifications
    pos_matrix  <- gregexpr("[1-f]", df$pep_ivmod) %>%
      list_to_dataframe() %>% 
      data.frame(check.names = FALSE)
    
    pos_matrix[pos_matrix < 0L] <- NA
    
    for (k in 1:ncol(pos_matrix)) {
      rows <- !is.na(pos_matrix[, k])
      locales <- pos_matrix[rows, k]
      lowers <- tolower(substr(df$pep_seq_mod[rows], locales, locales))
      substr(df$pep_seq_mod[rows], locales, locales) <- lowers
    }
    
    # (2-1) add "_" to sequences from protein N-terminal acetylation
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      ifelse(grepl("Acetyl (Protein N-term)", pep_vmod, fixed = TRUE), 
                             paste0("_", pep_seq_mod), 
                             pep_seq_mod))
    
    # (2-2) add "_" to sequences from protein C-terminal amidation
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      ifelse(grepl("Amidated (Protein C-term)", pep_vmod, fixed = TRUE), 
                             paste0(pep_seq_mod, "_"), 
                             pep_seq_mod))
    
    # (3-1) "~" for "(Protein N-term)" other than acetylation
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      ifelse(grepl("Protein N-term", pep_vmod) & 
                               !grepl("Acetyl (Protein N-term)", pep_vmod, fixed = TRUE), 
                             paste0("~", pep_seq_mod), 
                             pep_seq_mod))
    
    # (3-2) "~" for "(Protein C-term)" other than amidation
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      ifelse(grepl("Protein C-term", pep_vmod) & 
                               !grepl("Amidated (Protein C-term)", pep_vmod, fixed = TRUE), 
                             paste0(pep_seq_mod, "~"), 
                             pep_seq_mod))
    
    # (4-1) "^" for peptide "(N-term)"
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      ifelse(grepl("N-term", pep_vmod) & 
                               !grepl("Protein N-term", pep_vmod, fixed = TRUE), 
                             gsub("(^[_~]{0,1})", paste0("\\1", "^"), pep_seq_mod), 
                             pep_seq_mod))
    
    # (4-2) "^" peptide "(C-term)"
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = 
                      ifelse(grepl("C-term", pep_vmod) & 
                               !grepl("Protein C-term", pep_vmod, fixed = TRUE), 
                             gsub("([_~]{0,1}$)", paste0("^", "\\1"), pep_seq_mod), 
                             pep_seq_mod))
  } 
  else {
    df <- df %>%
      dplyr::mutate(pep_seq_mod = paste0(pep_seq, "[", pep_ivmod, "]"))
  }
  
  invisible(df)
}


#' Splits PSM tables from \link{matchMS}.
#' 
#' @inheritParams splitPSM
splitPSM_pq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                        fasta = NULL, entrez = NULL, pep_unique_by = "group", 
                        scale_rptr_int = FALSE, 
                        rm_craps = FALSE, rm_krts = FALSE, rm_allna = FALSE, 
                        purge_phosphodata = TRUE, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, 
                        rptr_intco = 0, rptr_intrange = c(0, 100), 
                        use_lowercase_aa = TRUE, parallel = TRUE, ...) 
{
  # --- Outlines ---
  # (1.1) row filtration, column padding and psm file combinations
  # 
  # (1.2) add ratio columns
  # 
  # (1.3) clean up prot_acc (1::, 2::)
  # convenience crap removals where the fasta names ended with "|"
  # (this can affect the uniqueness of peptides: 
  # by excluding the assessment of peptide sharedness with cRAPs)
  # 
  # (1.4) add pep_seq_mod
  # 
  # (2.1) annotate proteins
  # dependency: `gene` for prot_n_pep etc.
  # (NA genes are replaced with accessions)
  # (not yet prot_cover and prot_icover as they require pep_start and pep_end)
  # 
  # (2.2) compile "preferred" columns
  # 
  # (2.3) find the shared prot_accs and genes for each peptide
  # 
  # (2.4) remove non-essential proteins after shared_prot_accs, shared_genes
  # 
  # (2.5) find the levels of uniqueness for peptides
  # uniqueness currently by prot_acc, not gene
  # 
  # (2.6) add peptide properties and prot_cover, prot_icover
  # dependency: prot_cover and prot_icover need pep_start and pep_end
  # 
  # (2.7) apply parsimony
  # (a) chimeric: 
  #     the same `pep_query` can be assigned to different `pep_seq` at different `pep_rank`
  # (b) positional difference:
  #     the same combination of `pep_query`, `pep_seq` and `pep_var_mod_pos` 
  #     can be assigned to different `prot_acc` 
  # 
  # (c) remove redundant peptides under each `dat_file`
  #     dbl-dipping peptides across `dat_file` will be handled in `mergePep`
  # 
  # (2.8) add columns pep_n_psm, prot_n_psm, prot_n_pep
  # dependency: after the removal of duplicated PSM entries (under the same 
  #   dat_file + RAW_File + pep_scan_num + pep_seq)
  # 
  # (2.x) other fields
  # for compatibility, updates after data merging etc.
  # --- End ---
  
  on.exit(
    if (exists(".savecall", envir = environment())) {
      if (.savecall) 
        message("Split PSMs by TMT experiments and LCMS series --- Completed.")
    }, 
    add = TRUE
  )
  
  # ---
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  filelist <- find_psmQ_files(dat_dir)

  # (1.1) row filtration, column padding and psm file combinations
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in \"psmQ[...].txt\" for \"filter_\" varargs.", "\n")
  message("Found PSM files: \n  ", paste(filelist, collapse = ", \n  "), "\n")
  
  df <- filelist %>% 
    purrr::map(pad_tmt_channels, !!!filter_dots) %>% 
    pad_psm_fields() %>% 
    dplyr::bind_rows()
  
  col_nms <- names(df)
  
  lapply(c("RAW_File", "dat_file", 
           "pep_tot_int", "pep_vmod", "pep_ivmod", 
           "prot_hit_num", "prot_family_member"), 
         function (x) {
           if (! x %in% col_nms) stop("Column \"", x, "\" not found.")
         })
  
  # (has different meaning in Mascot)
  if ("pep_frame" %in% col_nms) 
    df$pep_frame <- NULL
  
  # (back-compatibility)
  if ("pep_ret_time" %in% col_nms)
    df <- dplyr::rename(df, pep_ret_range = pep_ret_time)
  
  if ("pep_ms2_n" %in% col_nms)
    df <- dplyr::rename(df, pep_n_ms2 = pep_ms2_n)
  
  # maybe new set cover of pep_seq's by prot_acc's here at length(filelist) > 1L; 
  # add_prot_acc(df) to avoid dependency w.r.t. terminal vs interior sequences,
  #   but need path info to access the cached lookups

  # (1.2) add ratio columns
  if (TMT_plex) {
    if (! "I126" %in% col_nms)
      stop("Column \"I126\" not found.", call. = FALSE)

    if (! TMT_plex %in% c(6L, 10L, 11L, 16L, 18L)) {
      stop("For TMT-plexes other than 6, 10, 11, 16 and 18, \n", 
           "  pad channels to the applicble higher plexes (e.g., 8 -> 10 or 11).", 
           call. = FALSE)
    }

    df <- sweep(df[, find_int_cols(TMT_plex), drop = FALSE], 
                1, df[["I126"]], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}[NC]{0,1}", names(.))), 
                       ~ replace(.x, is.infinite(.x), NA_real_)) %>% 
      dplyr::bind_cols(df, .) 
  } 
  else {
    if (! "I000" %in% col_nms)
      stop("Column `I000` not found.", call. = FALSE)

    df <- df %>% 
      reloc_col_after_last("I000") %>% 
      dplyr::mutate(R000 = I000/I000, 
                    R000 = ifelse(is.infinite(R000), NA_real_, R000)) 
  }
  
  rm(list = "col_nms")
  col_nms2 <- names(df)
  
  if (!any(grepl("^I[0-9]{3}[NC]{0,1}", col_nms2)))
    stop("Intensity column(s) \"I[...]\" not found.", call. = FALSE)
  
  if (!any(grepl("^R[0-9]{3}[NC]{0,1}", col_nms2)))
    stop("Ratio column(s) \"R[...]\" not found.", call. = FALSE)
  
  # (1.3) clean up prot_acc
  df <- df %>% 
    { if (rm_craps) dplyr::filter(., !grepl("\\|$", prot_acc)) else . } 
  
  # (1.4) add pep_seq_mod
  # (removals NL indicators at the end of `pep_ivmod`; 
  df$pep_ivmod <- with(df, gsub(" .*", "", pep_ivmod))

  df <- df %>% 
    split(.$dat_file, drop = TRUE) %>% 
    purrr::map(add_pepseqmod, use_lowercase_aa, purge_phosphodata) %>% 
    dplyr::bind_rows() 
  
  # (1.4.2) phospho
  # (may be deleted later)
  if (all(c("pep_vmod", "pep_locprob", "pep_locdiff") %in% col_nms2)) {
    df <- df %>% 
      dplyr::mutate(pep_phospho_locprob = 
                      ifelse(grepl("Phospho", pep_vmod), 
                             pep_locprob, 
                             NA_real_)) %>% 
      dplyr::mutate(pep_phospho_locdiff = 
                      ifelse(grepl("Phospho", pep_vmod), 
                             pep_locdiff, 
                             NA_real_))
  }
  
  rm(list = c("col_nms2"))
  
  # (2.1) annotate proteins
  df <- df %>% 
    annotPrn(fasta, entrez) %>%  
    { if (! "gene" %in% names(.)) dplyr::mutate(., gene = prot_acc) else . } %>% 
    dplyr::mutate(gene = ifelse(is.na(gene), prot_acc, gene))

  # (2.2) compile "preferred" columns
  df <- df %>% 
    dplyr::arrange(prot_hit_num, prot_family_member)
  
  # (2.3) find the shared prot_accs and genes for each peptide
  df <- df %>% 
    find_shared_prots("pep_seq", "prot_acc") %>% 
    find_shared_prots("pep_seq", "gene") %>% 
    dplyr::select(-which(names(.) %in% c("pep_seq_term", "pep_n_prot_accs", 
                                         "pep_n_genes")))
  
  # (2.4) remove non-essential proteins after shared_prot_accs, shared_genes
  df <- df %>% 
    dplyr::filter(!is.na(prot_family_member))

  # (2.5) find the uniqueness of peptides
  df <- if (pep_unique_by == "group") 
    dplyr::mutate(df, pep_isunique = pep_razor_unique)
  else if (pep_unique_by == "protein") 
    dplyr::mutate(df, pep_isunique = pep_literal_unique)
  else if (pep_unique_by == "none") 
    dplyr::mutate(df, pep_isunique = TRUE)
  else
    df
  
  # 0, not NA, since we known... even pep_tot_int is NA
  df <- df %>% 
    dplyr::mutate(pep_unique_int = 
                    ifelse(pep_literal_unique, pep_tot_int, 0)) %>% 
    dplyr::mutate(pep_razor_int = 
                    ifelse(pep_razor_unique, pep_tot_int, 0)) %>% 
    reloc_col_after("pep_unique_int", "pep_tot_int") %>% 
    reloc_col_after("pep_razor_int", "pep_unique_int")
  
  stopifnot(all(c("pep_isunique", "pep_literal_unique", "pep_razor_unique", 
                  "pep_tot_int", "pep_unique_int", "pep_razor_int") %in% 
                  names(df)))
  
  # (2.6) add peptide properties and prot_cover, prot_icover
  df <- df %>% 
    annotPeppos() %>% 
    dplyr::mutate(pep_len = stringr::str_length(pep_seq)) %>% 
    dplyr::mutate(pep_miss = ifelse(grepl("[KR]$", pep_seq), 
                                    stringr::str_count(pep_seq, "[KR]") - 1L,
                                    stringr::str_count(pep_seq, "[KR]"))) %>% 
    dplyr::mutate(pep_istryptic = as.logical(pep_istryptic), 
                  pep_miss = as.integer(pep_miss)) %>% 
    add_prot_icover(id = group_pep_by) %>% 
    calc_cover(id = !!rlang::sym(group_pep_by)) 

  # (2.7) apply parsimony
  # redundant prot_acc's removed with the heaviest being kept
  df <- df %>% 
    dplyr::arrange(pep_rank, -prot_mass)

  uniq_by <- c("RAW_File", "pep_scan_num", "pep_seq")
  
  if (length(unique(df$dat_file))) 
    uniq_by <- c("dat_file", uniq_by)
  
  # (often `pep_rank_nl == 1` after this step)
  rows <- !duplicated(df[, uniq_by])
  df <- df[rows, ]
  rm(list = c("rows"))
  
  # only for the demonstration of finer redundancy being kept at 
  #   "pep_seq_mod" and/or "neuloss"
  if (FALSE) {
    uniq_by2 <- c("RAW_File", "pep_scan_num", "pep_seq", "pep_ivmod")
    
    if (length(unique(df$dat_file))) 
      uniq_by2 <- c("dat_file", uniq_by2)

    ## handle ties in NL
    # 
    # the same `pep_ivmod` at different NLs with ties in `pep_rank` 
    #   -> keeps the best NL at EACH rank
    # 000050050, NL(1), 1; 000050050, NL(2), 1; 000050050, NL(3), 1
    # 000050050, NL(1), 2; 000050050, NL(2), 2
    #   -> pep_rank == 1: 000050050, NL(1), 1
    #   -> pep_rank == 2: 000050050, NL(1), 2
    df <- df %>% 
      dplyr::arrange(-pep_score) %>% 
      dplyr::group_by_at(c(uniq_by2, "pep_rank")) %>% 
      dplyr::filter(row_number() == 1L) %>% 
      dplyr::ungroup()
    
    ## keep the best pep_ivmod
    # (since the same raw, scan, pep_seq, & pep_ivmod, non-best NL are trivial)
    # 000050050, NL(1), 1
    # 000050050, NL(1), 2
    #   -> 000050050, NL(1), 1
    df <- df %>% 
      dplyr::group_by_at(c(uniq_by2)) %>% 
      dplyr::filter(row_number() == 1L)
    
    ## handle ties in pep_seq_mod
    if (psm_redundancy_at == "pep_seq_mod") {
      uniq_by3 <- uniq_by2[uniq_by2 != "pep_ivmod"]
      
      # keeps the first pep_seq_mod at each rank
      # pep_seq, 000050050, 1; pep_seq, 000050500, 1, pep_seq, 000055000, 1
      # pep_seq, 000050060, 2; pep_seq, 000050600, 2
      #   -> pep_seq, 000050050, 1
      #   -> pep_seq, 000050060, 2
      df <- df %>% 
        dplyr::group_by_at(c(uniq_by3, "pep_rank")) %>% 
        # dplyr::arrange(-pep_score) %>% 
        dplyr::filter(row_number() == 1L) %>% 
        dplyr::ungroup() 
      
      # keeps the best `pep_seq_mod` under the same scan
      # (since the same raw, scan & pep_seq, non-best pep_ivmod are trivial)
      # pep_seq, 000050050, 1
      # pep_seq, 000050060, 2
      #   -> pep_seq, 000050050, 1
      df <- df %>% 
        dplyr::group_by_at(c(uniq_by3)) %>% 
        # dplyr::arrange(-pep_score) %>% 
        dplyr::filter(row_number() == 1L) 
    }
  }

  df <- df %>% 
    dplyr::arrange(prot_hit_num, prot_family_member, pep_start, pep_end)
  
  # (2.8) add columns pep_n_psm, prot_n_psm, prot_n_pep
  # (after the non-redundant PSM entries)
  df <- df %>% 
    add_quality_cols(!!group_psm_by, !!group_pep_by, uniq_by)

  .saveCall <- TRUE
  
  invisible(df)
}


#' Adds the standard deviation in the retention times of peptides.
#' 
#' @param df A data frame.
#' @param use_unique Logical; if TRUE, filter data by uniqueness.
#' @inheritParams normPSM
add_pep_retsd <- function (df, group_psm_by = "pep_seq", use_unique = TRUE) 
{
  if (!group_psm_by %in% names(df)) {
    # warning("Column \"", group_psm_by, "\" not found.", call. = FALSE)
    return(df)
  }
  
  if (!"pep_ret_range" %in% names(df)) {
    # warning("Column \"pep_ret_range\" not available for assessing ", 
    #         "the standard deviations of peptide retention times.", 
    #         call. = FALSE)
    return (df)
  }
  
  if (all(is.na(df$pep_ret_range))) {
    # warning("All NA values under \"pep_ret_range\".", call. = FALSE)
    df <- df %>% dplyr::mutate(pep_ret_sd = NA_real_)
  }
  else {
    pep_ret <- calc_pep_retsd(df, group_psm_by, use_unique)

    df <- df %>% 
      # in case `pep_ret_sd` already in `df`
      dplyr::select(-which(names(df) == "pep_ret_sd")) %>% 
      dplyr::left_join(pep_ret, by = group_psm_by) 
  }
  
  df <- df %>% 
    reloc_col_after("pep_ret_sd", "pep_ret_range") %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
}


#' Calculates the standard deviation in the retention times of peptides.
#'
#' Sets \code{use_unique = TRUE} when used with normPSM where the same
#' \code{group_psm_by} can be redundant by \code{prot_acc}. Sets
#' \code{use_unique = FALSE} when used with mergePep where the same
#' \code{group_psm_by} from different samples and LCMS series are all taken into
#' account.
#'
#' @inheritParams normPSM
#' @inheritParams add_pep_retsd
calc_pep_retsd <- function (df, group_psm_by = "pep_seq", use_unique = TRUE) 
{
  pep_ret <- df[, c(group_psm_by, "pep_ret_range")] %>% 
    { if (use_unique) unique(.) else . } %>% 
    dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
    dplyr::summarise(pep_ret_sd = sd(pep_ret_range, na.rm = TRUE)) %>% 
    dplyr::mutate(pep_ret_sd = round(pep_ret_sd, digits = 2L)) %>% 
    dplyr::mutate(pep_ret_sd = ifelse(is.na(pep_ret_sd), 0, pep_ret_sd))
}


#' Calculates the number of unique charge states of peptides.
#' 
#' @param df A data frame.
#' @param use_unique Logical; if TRUE, filter data by uniqueness.
#' @inheritParams normPSM
add_n_pepexpz <- function (df, group_psm_by = "pep_seq", use_unique = TRUE) 
{
  if (!group_psm_by %in% names(df)) {
    # warning("Column \"", group_psm_by, "\" not found.", call. = FALSE)
    return(df)
  }
  
  if (!"pep_exp_z" %in% names(df)) {
    # warning("Column \"pep_exp_z\" not available for assessing ", 
    #         "the number of unique peptide charge states.", 
    #         call. = FALSE)
    return (df)
  }
  
  if (all(is.na(df$pep_exp_z))) {
    # warning("All NA values under `pep_exp_z`.", call. = FALSE)
    df <- df %>% dplyr::mutate(pep_n_exp_z = NA_integer_)
  }
  else {
    pep_z <- df[, c(group_psm_by, "pep_exp_z")] %>% 
      { if (use_unique) unique(.) else . } %>% 
      dplyr::select(group_psm_by) %>% 
      dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
      dplyr::summarise(pep_n_exp_z = n())
    
    df <- df %>% 
      # in case `pep_n_exp_z` already in `df`
      dplyr::select(-which(names(df) == "pep_n_exp_z")) %>% 
      dplyr::left_join(pep_z, by = group_psm_by) 
  }
  
  df <- df %>% 
    reloc_col_after("pep_n_exp_z", "pep_exp_z") %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
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





