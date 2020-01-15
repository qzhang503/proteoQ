#' TMT labeling efficiency
#' @import dplyr rlang ggplot2
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams cleanupPSM
#' @inheritParams annotPSM
#' @inheritParams normPSM
#' @export
labEffPSM <- function(group_psm_by = c("pep_seq", "pep_seq_mod"), group_pep_by = c("prot_acc", "gene"), 
                      dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx", 
                      fasta = NULL, pep_unique_by = "group", corrected_int = TRUE, rm_reverses = TRUE, 
                      rptr_intco = 1000, rm_craps = FALSE, rm_krts = FALSE, rm_outliers = FALSE, 
                      annot_kinases = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                      use_lowercase_aa = TRUE, ...) {
  
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
  
  rmPSMHeaders()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
                        pattern = "^F[0-9]{6}\\_hdr_rm.csv$")
  
  if (length(filelist) == 0) stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
  df <- purrr::map(filelist, ~ {
    df <- read.delim(file.path(dat_dir, "PSM\\cache", .x), sep = ',', check.names = FALSE, 
                     header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    df$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    
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
    
    dat_id <- df$dat_file %>% unique()
    dat_file <- file.path(dat_dir, "PSM\\cache", paste0(dat_id, "_header.txt"))
    stopifnot(length(dat_id)== 1, file.exists(dat_file))
    
    df_header <- readLines(dat_file)
    
    tmt_line <- df_header %>% 
      .[grep("Fixed modifications", .)] %>% 
      .[grep("TMT[0-9]{1}plex", .)]
    
    if (rlang::is_empty(tmt_line)) {
      tmt_type <- "neither"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(K\\)", tmt_line) & grepl("TMT[0-9]{1}plex\\s\\(N-term\\)", tmt_line)) {
      tmt_type <- "both"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(K\\)", tmt_line)) {
      tmt_type <- "k_only"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(N-term\\)", tmt_line)) {
      tmt_type <- "nt_only"
    } 
    
    if (tmt_type == "nt_only") {
      df <- df %>% 
        dplyr::filter(grepl("K", .$pep_seq))
    }
    
    df <- df %>% dplyr::mutate(tmt_type = tmt_type)
    
    return(df)
  }) %>% 
    do.call(rbind, .)
  
  # convenience craps removals where their uniprot afasta names ended with "|"
  if (rm_craps) df <- df %>% dplyr::filter(!grepl("\\|.*\\|$", prot_acc))
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  # note pep_seq: from such as MENGQSTAAK to K.MENGQSTAAK.L
  df <- df %>% 
    dplyr::mutate(pep_len = str_length(pep_seq)) %>% 
    split(., .$dat_file, drop = TRUE) %>% 
    purrr::map(add_mascot_pepseqmod, use_lowercase_aa) %>% 
    bind_rows() # %>% 
  # dplyr::select(-dat_file)
  
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
      # dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      # dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>%
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) # %>%
    # dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>% 
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  df_split <- df_split %>% split(., .$RAW_File, drop = TRUE)
  
  dat_files <- df %>% 
    dplyr::select(dat_file, tmt_type) %>% 
    dplyr::filter(!duplicated(dat_file))
  
  lab_eff <- purrr::map2_dbl(df_split, names(df_split), ~ {
    res <- .x %>% 
      dplyr::group_by(dat_file) %>% 
      dplyr::summarize(N = n()) %>% 
      dplyr::right_join(dat_files, by = "dat_file") %>% 
      dplyr::arrange(-N, dat_file) %T>% 
      write.csv(file.path(dat_dir, "PSM\\cache", paste0(.y, ".csv")), row.names = FALSE)
    
    both <- res %>% 
      dplyr::filter(tmt_type == "both") %>% 
      dplyr::mutate(N = 2 *N) %>% 
      .[["N"]]
    
    k_only <- res %>% 
      dplyr::filter(tmt_type == "k_only") %>% 
      .[["N"]]
    
    nt_only <- res %>% 
      dplyr::filter(tmt_type == "nt_only") %>% 
      .[["N"]]
    
    neither <- res %>% 
      dplyr::filter(tmt_type == "neither") %>% 
      dplyr::mutate(N = 2 *N) %>% 
      .[["N"]]
    
    1 - (k_only + nt_only + neither) / (both + k_only + nt_only)
  })
}



