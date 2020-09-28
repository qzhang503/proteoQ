#' Make new column names
#'
#' \code{newColnames} match names to Sample_ID in label_scheme
#'
#' @param i Integer; the index of TMT experiment 
#' @param x Data frame; log2FC data
#' @param label_scheme Experiment summary
#' @import dplyr purrr forcats
#' @importFrom magrittr %>% %T>% %$% %<>% 
newColnames <- function(i, x, label_scheme, pattern = "[RI][0-9]{3}[NC]{0,1}") {
  label_scheme_sub <- label_scheme %>%
    dplyr::filter(TMT_Set == i) 
  
  cols <- grep(paste0(pattern, "_", i, "$"), names(x))
  nm_channel <- gsub(paste0("(", pattern, ")_", i, "$"), "\\1", names(x)[cols])
  names(x)[cols] <- paste0(nm_channel, " (", as.character(label_scheme_sub$Sample_ID), ")")
  
  # update the column indexes with TMT_Set being replaced with Sample_ID
  cols <- grep(paste0(pattern, "\\s+\\(.*\\)$"), names(x))
  
  # cols with the updated names go first
  if (length(cols) < ncol(x)) x <- dplyr::bind_cols(x[, cols], x[, -cols, drop = FALSE])
  
  return(x)
}


#' Check the single presence of MaxQuant peptides[...].txt.
#' 
#' @inheritParams load_expts
use_mq_peptable <- function (dat_dir, label_scheme_full) {
  filelist <- list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  
  if (length(filelist) == 0) {
    message("Primary column keys in `Peptide/TMTset1_LCMSinj1_Peptide_N.txt` etc. ", 
            "for `filter_` varargs.")
    return(FALSE)
  } else if (length(filelist) > 1) {
    stop("Only single MaxQuant `peptides.txt` allowed.", call. = FALSE)
  } else {
    if (!file.exists(file.path(dat_dir, filelist))) {
      stop("No MaxQuant LFQ file `peptides[...].txt` udner", dat_dir, call. = FALSE)
    } else {
      message("MaxQuant file `", filelist, "` found.")
      return(TRUE)
    }
  }
}


#' Check the single presence of MaxQuant proteinGroups[...].txt.
#' 
#' @inheritParams load_expts
single_mq_prntable <- function (dat_dir) {
  filelist <- list.files(path = file.path(dat_dir), pattern = "^proteinGroups.*\\.txt$")
  
  if (length(filelist) > 1) {
    stop("Only single MaxQuant `proteinGroups.txt` allowed.", call. = FALSE)
  }
  
  if (!file.exists(file.path(dat_dir, filelist))) {
    stop("No MaxQuant LFQ file `proteinGroups[...].txt` udner", dat_dir, call. = FALSE)
  }
  
  message("MaxQuant file `", filelist, "` found.")  
  
  invisible(filelist)
}


#' Compare sample IDS between MaxQuant results and expt_smry.xlsx 
#' 
#' @param df A data frame of MaxQuant peptides.txt or proteinGroups.txt.
#' @inheritParams n_TMT_sets
check_mq_df <- function (df, label_scheme) {
  mq_nms <- names(df) %>% 
    .[grepl("^LFQ intensity ", .)] %>% 
    gsub("^LFQ intensity ", "", .)
  
  ls_nms <- label_scheme$Sample_ID %>% 
    .[!grepl("^Empty\\.[0-9]+", .)] %>% 
    unique() 
  
  if (!all(mq_nms %in% ls_nms)) {
    missing_nms <- mq_nms %>% .[! . %in% ls_nms]
    
    warning("Sample ID(s) in MaxQuant `peptides.txt` not found in metadata:\n", 
            purrr::reduce(missing_nms, paste, sep = ", "), 
            call. = FALSE)

    # the same ID occurs in "Identification type", "Experiment", "Intensity", "LFQ intensity"
    purrr::walk(missing_nms, ~ {
      df[, grep(paste0(" ", .x, "$"), names(df))] <- NULL
      df <<- df
    }, df)
  }
  
  invisible(df)
}


#' Extracts MaxQuant intensity values and calculates log2FC.
#' 
#' @param df A data frame of MaxQuant results.
extract_mq_ints <- function (df) {
  calc_mq_log2r <- function (df, type, refChannels) {
    type <- paste0("^", type, " ")
    col_smpls <- grep(type, names(df))
    
    if (length(refChannels) > 0) {
      col_refs <- paste0(type, refChannels) %>% 
        purrr::map_dbl(~ grep(.x, names(df)))
    } else {
      col_refs <- grep(type, names(df)) 
    }
    
    if (type == "^Intensity ") {
      prefix <- "log2_R000"
    } else if (type == "^LFQ intensity ") {
      prefix <- "N_log2_R000"
    } else {
      stop("`type` needs to be either `Intensity` or `LFQ intensity`.", call. = FALSE)
    }
    
    sweep(df[, col_smpls], 1,
          rowMeans(df[, col_refs, drop = FALSE], na.rm = TRUE), "/") %>%
      log2(.)  %>% 
      dplyr::mutate_all(~ replace(.x, is.infinite(.), NA)) %>% 
      `colnames<-`(gsub(paste0(type, "(.*)$"), paste0(prefix, " \\(", "\\1", "\\)"), names(.))) %>%
      cbind(df, .)
  }
  
  load(file.path(dat_dir, "label_scheme.rda"))
  
  refChannels <- label_scheme %>% 
    dplyr::filter(Reference) %>% 
    dplyr::select(Sample_ID) %>% 
    unlist()
  
  df <- df %>% 
    calc_mq_log2r("Intensity", refChannels) %>% 
    calc_mq_log2r("LFQ intensity", refChannels) %>% 
    dplyr::mutate_at(.vars = grep("log2_R000\\s", names(.)), ~ replace(.x, is.infinite(.), NA)) 
  
  log2sd <- df %>% 
    dplyr::select(grep("Intensity\\s", names(.))) %>% 
    `names<-`(gsub("^Intensity\\s(.*)", "sd_log2_R000 \\(\\1\\)", names(.))) %>% 
    dplyr::mutate_all(~ replace(.x, !is.na(.x), NA))
  
  df <- dplyr::bind_cols(df, log2sd)
  
  df <- df %>% 
    `names<-`(gsub("^Intensity (.*)$", paste0("I000 \\(", "\\1", "\\)"), names(.))) %>% 
    `names<-`(gsub("^LFQ intensity (.*)$", paste0("N_I000 \\(", "\\1", "\\)"), names(.)))
  
  if (purrr::is_empty(grep("log2_R000", names(df)))) {
    stop("No `log2_R000...` columns available.\n",
         "Probably inconsistent sample IDs between metadata and MaxQuant `peptides.txt`.", 
         call. = FALSE)
  }
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[IR]{1}000 \\(", names(.))),
    df %>% dplyr::select(grep("^I000 \\(", names(.))),
    df %>% dplyr::select(grep("^N_I000 \\(", names(.))),
    df %>% dplyr::select(grep("^sd_log2_R000 \\(", names(.))),
    df %>% dplyr::select(grep("^log2_R000 \\(", names(.))),
    df %>% dplyr::select(grep("^N_log2_R000 \\(", names(.))),  
  )
}


#' Handling of MaxQuant peptide.txt
#'
#' Fill back temporarily intensity values that are not filled in LFQ msms.txt.
#' @inheritParams n_TMT_sets
pep_mq_lfq <- function(label_scheme, omit_single_lfq) {
  dat_dir <- get_gl_dat_dir()
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  
  if (purrr::is_empty(filelist)) {
    stop(paste("No MaxQuant LFQ file of `peptides[...].txt` under", 
               file.path(dat_dir), ".\n",
               "Make sure that the name of file starts with `peptides`."), 
         call. = FALSE)
  }
  
  df <- read.csv(file.path(dat_dir, filelist), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(not_allzero_rows(.[grep("^LFQ intensity ", names(.))])) 
  
  if (omit_single_lfq) {
    df <- df %>% na_single_lfq("^LFQ intensity ")
  }

  ## handle inconsistency in MaxQuant column keys
  # (1) "Gene names" vs "Gene Names" etc.
  # (2) presence or absence of "Gene Names", "Protein Names" etc.
  if (!("Gene Names" %in% names(df) && "Protein Names" %in% names(df))) {
    stopifnot("Proteins" %in% names(df))
    
    df$"Gene names" <- df$"Gene Names" <- df$"Protein names" <- df$"Protein Names" <- NULL
    
    fasta <- match_call_arg(normPSM, fasta)
    entrez <- match_call_arg(normPSM, entrez)
    
    df <- local({
      df <- df %>% dplyr::mutate(prot_acc = gsub(";.*$", "", Proteins))
      
      tempdata <- df %>% 
        dplyr::select(prot_acc) %>% 
        dplyr::filter(!duplicated(prot_acc)) %>% 
        annotPrn(fasta, entrez) %>% 
        dplyr::select(prot_acc, gene, prot_desc) %>% 
        dplyr::rename(`Gene Names` = "gene", "Protein Names" = "prot_desc")
      
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
  
  df <- df %>% check_mq_df(label_scheme) %>% 
    dplyr::rename(
      pep_seq = Sequence, 
      prot_acc = Proteins, 
    ) %>% 
    dplyr::mutate(gene = gsub("\\;.*", "", `Gene Names`)) %>% 
    dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc))

  # `Modified sequence` not available in `peptides.txt`
  # group_psm_by = "pep_seq" only in `normPSM()`
  if (group_psm_by == "pep_seq_mod") {
    if ("Modified sequence" %in% names(df)) {
      use_lowercase_aa <- match_call_arg(normPSM, use_lowercase_aa)
      
      df <- df %>% 
        dplyr::rename(pep_seq_mod = `Modified sequence`) %>% 
        dplyr::select(which(names(.) == group_psm_by),
                      which(names(.) != group_psm_by)) %>% 
        add_maxquant_pepseqmod(use_lowercase_aa = FALSE)
    } else {
      stop("Column `Modified sequence` not found in MaxQuant `peptides.txt`.\n", 
           "Rerun `normPSM(group_psm_by = pep_seq, ...)`.\n", 
           call. = FALSE)
    }
  } else {
    df <- df %>% 
      dplyr::mutate(pep_seq = paste(`Amino acid before`, pep_seq, `Amino acid after`, sep = ".")) 
  }
  
  df <- local({
    df_vals <- df %>% 
      dplyr::select(group_psm_by, grep("^Intensity\\s|^LFQ\\sintensity\\s", names(.))) %>% 
      extract_mq_ints() %>% 
      dplyr::select(-grep("sd_log2_R000", names(.)))
    
    df_sds <- df %>% 
      dplyr::left_join(df_vals, by = group_psm_by) %>% 
      calcSD_Splex(group_pep_by) %>% 
      `names<-`(gsub("^log2_R", "sd_log2_R", names(.)))
    
    df %>% 
      dplyr::left_join(df_vals, by = group_psm_by) %>% 
      dplyr::left_join(df_sds, by = group_pep_by) %>% 
      na_zeroIntensity() 
  })

  df %>% 
    dplyr::select(
      group_psm_by, 
      grep("^sd_log2_R000", names(.)), 
      grep("^log2_R000", names(.)), 
      grep("^N_log2_R000", names(.)), 
      grep("^I000", names(.)), 
      grep("^N_I000", names(.)), 
    )
}


#' Total, razor and unique intensity of peptides 
#' 
#' Temporary handling using MaxQuant peptide.txt.
pep_mq_lfq2 <- function() {
  dat_dir <- get_gl_dat_dir()
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  corrected_int <- match_call_arg(normPSM, corrected_int)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  
  if (purrr::is_empty(filelist)) {
    stop(paste("No MaxQuant LFQ file of `peptides[...].txt` under", 
               file.path(dat_dir), ".\n",
               "Make sure that the name of file starts with `peptides`."), 
         call. = FALSE)
  }
  
  df <- read.csv(file.path(dat_dir, filelist), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(not_allzero_rows(.[grep("^LFQ intensity ", names(.))])) 
  
  stopifnot("Proteins" %in% names(df), 
            "Sequence" %in% names(df), 
            "Unique (Groups)" %in% names(df), 
            "Unique (Proteins)" %in% names(df), 
            "Amino acid before" %in% names(df), 
            "Amino acid after" %in% names(df))
  
  df <- df %>% dplyr::mutate(prot_acc = gsub(";.*$", "", Proteins))
  
  prefix <- if (corrected_int) "LFQ intensity" else "^Intensity"
  
  df_slim <- df %>% 
    dplyr::mutate(Sequence = paste(`Amino acid before`, Sequence, `Amino acid after`, 
                                   sep = ".")) %>% 
    dplyr::select(Sequence, `Unique (Groups)`, `Unique (Proteins)`, 
                  grep(prefix, names(.)))
  
  df_tot <- df_slim %>% 
    `names<-`(gsub(paste0(prefix, " (.*)"), 
                   paste0("pep_tot_int \\(", "\\1", "\\)"), 
                   names(.))) %>% 
    dplyr::rename(pep_seq = Sequence, 
                  uniq_grp = `Unique (Groups)`, 
                  uniq_prot = `Unique (Proteins)`) %>% 
    dplyr::mutate(uniq_grp = ifelse(uniq_grp == "yes", TRUE, FALSE)) %>% 
    dplyr::mutate(uniq_prot = ifelse(uniq_prot == "yes", TRUE, FALSE))
  
  df_razor <- df_tot %>% 
    dplyr::mutate_at(grep("^pep_tot_int", names(.)), ~ ifelse(uniq_grp, .x, 0)) %>% 
    `names<-`(gsub("pep_tot_int", "pep_razor_int", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("pep_seq", "uniq_grp", "uniq_prot")))
  
  df_uniq <- df_tot %>% 
    dplyr::mutate_at(grep("^pep_tot_int", names(.)), ~ ifelse(uniq_prot, .x, 0)) %>% 
    `names<-`(gsub("pep_tot_int", "pep_unique_int", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("pep_seq", "uniq_grp", "uniq_prot")))
  
  dplyr::bind_cols(df_tot, df_razor, df_uniq) %>% 
    dplyr::select(-which(names(.) %in% c("uniq_grp", "uniq_prot")))
}


#' load MaxQuant protein table
#' 
#' @param label_scheme Experiment summary
prn_mq_lfq <- function(label_scheme) {
  dat_dir <- get_gl_dat_dir()
  group_pep_by <- "prot_acc"
  
  filelist <- single_mq_prntable(dat_dir)
  
  df <- read.csv(file.path(dat_dir, filelist), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(not_allzero_rows(.[grep("^LFQ intensity ", names(.))]))
  
  df <- check_mq_df(df, label_scheme)
  
  df <- df %>% 
    dplyr::rename(
      prot_acc = "Majority protein IDs", 
      gene = "Gene names", 
    ) %>% 
    dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc), 
                  gene = gsub("\\;.*", "", gene))
  
  df <- df %>% 
    dplyr::select(group_pep_by, grep("^Intensity |^LFQ intensity ", names(.))) %>% 
    extract_mq_ints()
  
  # (1) calculate Z_log2_R
  sd_coefs <- df %>% calc_sd_fcts(range_log2r = c(5, 95), range_int = c(5, 95), label_scheme)
  
  x_vals <- df %>%
    dplyr::select(grep("^N_log2_R[0-9]{3}", names(.))) %>% 
    `colnames<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
    dplyr::summarise_all(~ median(.x, na.rm = TRUE)) %>%
    unlist() %>%
    data.frame(x = .) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
    dplyr::arrange(Sample_ID) %>% 
    dplyr::mutate(x = ifelse(is.na(x), NA, 0))
  
  df <- update_df(df, label_scheme, x_vals, sd_coefs)
  
  df %>% 
    dplyr::select(
      group_pep_by, 
      grep("^I", names(.)), 
      grep("^N_I", names(.)), 
      grep("^sd_log2_R", names(.)), 
      grep("^log2_R", names(.)), 
      grep("^N_log2_R", names(.)), 
      grep("^Z_log2_R", names(.)), 
    )
}


#' Helper: calculates LFQ log2FC.
#' 
#' @param df A data frame.
#' @param type The type of intensity data.
#' @param refChannels The reference channels.
calc_lfq_log2r <- function (df, type, refChannels) {
  stopifnot(type %in% c("I000", "N_I000"))
  
  new_type <- paste0("^", type, " ")
  prefix <- gsub("I000", "log2_R000", new_type)
  no_hat <- gsub("\\^", "", prefix)
  
  df <- df %>% dplyr::select(-grep(prefix, names(.)))

  col_smpls <- grep(new_type, names(df))
  
  if (purrr::is_empty(col_smpls)) {
    stop("No sample columns start with ", type, ".", call. = FALSE)
  }
  
  if (length(refChannels) > 0) {
    col_refs <- which(names(df) %in% paste0(type, " (", refChannels, ")"))
  } else {
    col_refs <- col_smpls
  }
  
  sweep(df[, col_smpls], 1,
        rowMeans(df[, col_refs, drop = FALSE], na.rm = TRUE), "/") %>%
    log2(.)  %>% 
    dplyr::mutate_all(~ replace(.x, is.infinite(.), NA)) %>% 
    `colnames<-`(gsub(paste0(new_type, "(.*)$"), paste0(no_hat, "\\1"), names(.))) %>% 
    cbind(df, .)
}


#' Trivializes data rows with single LFQ intensity
#' 
#' @param df A data frame.
#' @param pattern The pattern of intensity fields.
na_single_lfq <- function (df, pattern = "^I000 ") {
  df_lfq <- df[, grep(pattern, names(df))]
  not_single_zero <- (rowSums(df_lfq > 0, na.rm = TRUE) > 1) 
  not_single_zero[!not_single_zero] <- NA
  
  df_lfq[] <- purrr::map(df_lfq, `*`, not_single_zero)
  df[, grep(pattern, names(df))] <- df_lfq
  
  invisible(df)
}


#' Calculates log2FC of peptides based on LFQ intensity.
#' 
#' With LFQ, values of \code{log2_R000 (...)} and \code{N_log2_R000 (...)} in
#' \code{df_num} are not yet filled after \link{calc_tmt_nums}. This utility
#' calculates the ratio using the values of LFQ intensity.
#' 
#' @param df A data frame.
calc_lfq_pepnums <- function (df, omit_single_lfq) {
  load(file.path(dat_dir, "label_scheme.rda"))
  
  if (omit_single_lfq) {
    df <- df %>% 
      na_single_lfq("^I000 ") %>% 
      na_single_lfq("^N_I000 ")
  }
  
  refChannels <- label_scheme %>% 
    dplyr::filter(Reference) %>% 
    dplyr::select(Sample_ID) %>% 
    unlist()
  
  df <- df %>% 
    dplyr::mutate_if(is.numeric, ~ replace(.x, is.nan(.x), NA_real_)) %>% 
    calc_lfq_log2r("I000", refChannels) %>% 
    calc_lfq_log2r("N_I000", refChannels) %>% 
    dplyr::mutate_at(.vars = grep("log2_R000\\s", names(.)), 
                     ~ replace(.x, is.infinite(.), NA)) 

  if (purrr::is_empty(grep("log2_R000", names(df)))) {
    stop("No `log2_R000...` columns available.\n",
         "Probably inconsistent sample IDs between metadata and MaxQuant `peptides.txt`.", 
         call. = FALSE)
  }
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[IR]{1}000 \\(", names(.))),
    df %>% dplyr::select(grep("^I000 \\(", names(.))),
    df %>% dplyr::select(grep("^N_I000 \\(", names(.))),
    df %>% dplyr::select(grep("^sd_log2_R000 \\(", names(.))),
    df %>% dplyr::select(grep("^log2_R000 \\(", names(.))),
    df %>% dplyr::select(grep("^N_log2_R000 \\(", names(.))),  
  )
}


#' Calculates \code{pep_tot_int}, \code{pep_razor_int} and
#' \code{pep_unique_int} in LFQ
#'
#' @param df A data frame.
calc_lfq_pepnums2 <- function (df, filelist, group_psm_by) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  df_num <- df %>% 
    dplyr::select(!!rlang::sym(group_psm_by), 
                  TMT_Set,    
                  pep_tot_int, 
                  pep_unique_int,
                  pep_razor_int) %>% 
    dplyr::group_by(!!rlang::sym(group_psm_by), TMT_Set) %>%
    dplyr::arrange(TMT_Set) %>%
    tidyr::gather(grep("^pep_.*_int$", names(.)), key = ID, value = value) %>%
    tidyr::unite(ID, ID, TMT_Set)
  
  Levels <- unique(df_num$ID)
  df_num <- df_num %>%
    dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
    tidyr::spread(ID, value)
  rm(Levels)
  
  set_indexes <- gsub("^.*TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort()
  
  for (set_idx in set_indexes) {
    df_num <- newColnames(set_idx, df_num, label_scheme, "pep_.*_int")
  }
  
  df_num <- df_num %>% 
    dplyr::select(!!rlang::sym(group_psm_by), grep("^pep_.*_int ", names(.))) %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
}


#' Calculates numeric values
#' 
#' Fields of numeric values: sd_log2_R, log2_R, log2_R, I, N_I.
#'
#' @param df A data frame of PSM.
#' @param filelist A list of PSM files.
#' @inheritParams splitPSM
#' @inheritParams annotPSM
calc_tmt_nums <- function (df, filelist, group_psm_by, parallel) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))

  df_num <- df %>% 
    dplyr::select(!!rlang::sym(group_psm_by), 
                  TMT_Set, 
                  grep("^sd_log2_R[0-9]{3}", names(.)), 
                  grep("^log2_R[0-9]{3}", names(.)), 
                  grep("^N_log2_R[0-9]{3}", names(.)), 
                  # grep("^Z_log2_R[0-9]{3}", names(.)), 
                  grep("^I[0-9]{3}", names(.)), 
                  grep("^N_I[0-9]{3}", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(group_psm_by), TMT_Set)
  
  tbl_lcms <- n_LCMS(label_scheme_full)
  
  if (any(tbl_lcms$n_LCMS > 1)) {
    tbl_n <- tbl_lcms %>% dplyr::filter(n_LCMS > 1)
    df_n <- df_num %>% dplyr::filter(TMT_Set %in% tbl_n$TMT_Set)
    df_1 <- df_num %>% dplyr::filter(! TMT_Set %in% tbl_n$TMT_Set)
    
    if (parallel) {
      nms <- names(df_n)
      stopifnot(nms[1] == group_psm_by, nms[2] == "TMT_Set")
      
      n_cores <- parallel::detectCores()
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      df_num <- suppressWarnings(
        parallel::clusterApply(
          cl, 3:ncol(df_n), function(i) {
            df_n %>% 
              dplyr::select(c(1:2, i)) %>% 
              dplyr::summarise_all(~ median(.x, na.rm = TRUE))
          }
        )
      ) %>% 
        purrr::map(~ tidyr::unite(.x, id, group_psm_by, TMT_Set, sep = "@"))  %>% 
        purrr::reduce(dplyr::left_join, by = "id") %>% 
        tidyr::separate(id, into = c(group_psm_by, "TMT_Set"), sep = "@", remove = TRUE) %>% 
        dplyr::bind_rows(df_1) # the same order of columns ensured
      
      parallel::stopCluster(cl)
    } else {
      df_num <- df_n %>% 
        dplyr::summarise_all(~ median(.x, na.rm = TRUE)) %>% 
        dplyr::bind_rows(df_1)
    }
  } 
  
  df_num <- df_num %>%
    dplyr::arrange(TMT_Set) %>%
    tidyr::gather(grep("R[0-9]{3}|I[0-9]{3}", names(.)), key = ID, value = value) %>%
    tidyr::unite(ID, ID, TMT_Set)
  
  # define the levels of TMT channels;
  # otherwise, the order of channels will flip between N(itrogen) and C(arbon)
  Levels <- unique(df_num$ID)
  df_num <- df_num %>%
    dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
    tidyr::spread(ID, value)
  rm(Levels)
  
  set_indexes <- gsub("^.*TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort()
  
  for (set_idx in set_indexes) {
    df_num <- newColnames(set_idx, df_num, label_scheme)
  }
  
  df_num <- df_num %>% 
    dplyr::select(!!rlang::sym(group_psm_by), grep("[RI][0-9]{3}[NC]*", names(.))) %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
}


#' combined peptide reports across multiple TMT experiments
#' 
#' median summarization of data from the same TMT experiment at different LCMS injections
#' summed \code{pep_n_psm}, \code{prot_n_psm}, and \code{prot_n_pep} after data merging
#' no Z_log2_R yet available
#'   use \code{col_select = expr(Sample_ID)} not \code{col_select} to get all Z_log2_R
#'   why: users may specify \code{col_select} only partial to Sample_ID entries
#'
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @inheritParams splitPSM
#' @inheritParams mergePep
normPep_Mplex <- function (group_psm_by = "pep_seq_mod", group_pep_by = "prot_acc", 
                           use_duppeps = TRUE, cut_points = Inf, 
                           omit_single_lfq = TRUE, use_mq_pep = FALSE, 
                           parallel = TRUE, ...) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)

  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  filelist <- list.files(path = file.path(dat_dir, "Peptide"), 
                         pattern = paste0("TMTset[0-9]+_LCMSinj[0-9]+_Peptide_N\\.txt$"), 
                         full.names = TRUE)
  
  if (purrr::is_empty(filelist)) {
    stop("No individual peptide tables available; run `PSM2Pep()` first.", 
         call. = FALSE)
  }
  
  df <- do.call(rbind, 
                purrr::map(filelist, read.csv, check.names = FALSE, header = TRUE, 
                           sep = "\t", comment.char = "#")) %>% 
    dplyr::select(-which(names(.) %in% c("dat_file"))) %>% 
    dplyr::mutate(TMT_Set = factor(TMT_Set)) %>%
    dplyr::arrange(TMT_Set) 
  
  df <- df %>% filters_in_call(!!!filter_dots)
  
  if ("gene" %in% names(df)) {
    df <- df %>% dplyr::mutate(gene = forcats::fct_explicit_na(gene))
  }
  
  df <- df %>% assign_duppeps(group_psm_by, group_pep_by, use_duppeps)
  
  if (!use_mq_pep) {
    df_num <- calc_tmt_nums(df, filelist, group_psm_by, parallel)
    df_num2 <- df %>% dplyr::select(group_psm_by)
    
    if (TMT_plex == 0) {
      df_num <- calc_lfq_pepnums(df_num, omit_single_lfq)
      df_num2 <- calc_lfq_pepnums2(df, filelist, group_psm_by)
    } 
  } else {
    # temporarily back-fill from a peptide table
    df_num <- pep_mq_lfq(label_scheme, omit_single_lfq)
    df_num2 <- pep_mq_lfq2()
  }
  
  df_num <- dplyr::left_join(df_num2, df_num, by = group_psm_by)
  rm(df_num2)
  
  df_num <- local({
    count_nna <- df_num %>% 
      dplyr::select(grep("N_log2_R[0-9]{3}[NC]{0,1}", names(.)))%>% 
      dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Ref\\.[0-9]+\\)$", names(.))) %>% 
      dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Empty\\.[0-9]+\\)$", names(.))) %>% 
      is.na() %>% 
      `!`() %>% 
      rowSums()
    
    dplyr::bind_cols(count_nna = count_nna, df_num)
  })

  pep_n_psm <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), pep_n_psm) %>%
    dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
    dplyr::summarise(pep_n_psm = sum(pep_n_psm)) %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
  
  prot_n_psm <- df %>% 
    dplyr::select(pep_n_psm, !!rlang::sym(group_pep_by)) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_n_psm = sum(pep_n_psm))

  prot_n_pep <- df %>% 
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>% 
    dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_n_pep = n())
  
  df_first <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    med_summarise_keys(group_psm_by) %>% 
    dplyr::select(-which(names(.) %in% c("pep_n_psm", "prot_n_psm", 
                                         "prot_n_pep", "TMT_Set"))) %>% 
    dplyr::select(-which(names(.) %in% c("pep_tot_int", "pep_unique_int", 
                                         "pep_razor_int"))) %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))

  df <- list(pep_n_psm, df_first, df_num) %>%
    purrr::reduce(dplyr::left_join, by = group_psm_by) %>% 
    reloc_col_before(group_psm_by, "pep_res_after")
  
  df <- list(df, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = group_pep_by)
  
  # --- update sd_log2R000 (...) ---
  if (TMT_plex == 0 && !use_mq_pep) {
    df <- local({
      df <- df %>% 
        dplyr::select(-grep("^sd_log2_R000 ", names(.)))

      df_sds <- df %>% 
        calcSD_Splex(group_pep_by) %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.)))
      
      dplyr::right_join(df_sds, df, by = group_pep_by) %>% 
        na_zeroIntensity()
    })
  }
  
  # --- add varmod columns ---
  if (("pep_seq_mod" %in% names(df)) && (match_call_arg(normPSM, use_lowercase_aa))) {
    df <- df %>% 
      dplyr::mutate(pep_mod_protnt = ifelse(grepl("^[A-z\\-]\\.~", pep_seq_mod), 
                                            TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_protntac = ifelse(grepl("^[A-z\\-]\\._", pep_seq_mod), 
                                              TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_pepnt = ifelse(grepl("^[A-z\\-]\\.[_~]?\\^", pep_seq_mod), 
                                           TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_m = ifelse(grepl("m", pep_seq_mod), 
                                       TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_n = ifelse(grepl("n", pep_seq_mod), 
                                       TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_sty = ifelse(grepl("[sty]", pep_seq_mod), 
                                         TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_pepct = ifelse(grepl("[\\^]{1}[_~]?\\.[A-z\\-]{1}$", pep_seq_mod), 
                                           TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_protctam = ifelse(grepl("_{1}\\.[A-z\\-]{1}$", pep_seq_mod), 
                                              TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_protct = ifelse(grepl("~{1}\\.[A-z\\-]{1}$", pep_seq_mod), 
                                            TRUE, FALSE))
  } else {
    df <- df %>% 
      dplyr::mutate(pep_mod_protnt = NA, 
                    pep_mod_protntac = NA, 
                    pep_mod_pepnt = NA, 
                    pep_mod_m = NA, 
                    pep_mod_n = NA, 
                    pep_mod_sty = NA, 
                    pep_mod_pepct = NA, 
                    pep_mod_protctam = NA, 
                    pep_mod_protct = NA)
  }
  
  # --- tidy up ---
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))),
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^prot_|^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^Z_log2_R[0-9]{3}[NC]*", names(.))), 
  )
  
  df <- df %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
                     as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
                     ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), 
                     as.numeric) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), 
                     ~ round(.x, digits = 3)) %>% 
    dplyr::mutate_at(vars(grep("sd_log2_R[0-9]{3}[NC]*", names(.))), 
                     ~ round(.x, digits = 3))
  
  df <- df %>%
    dplyr::filter(!duplicated(.[[group_psm_by]])) %>% 
    dplyr::filter(rowSums(!is.na(.[, grepl("N_log2_R", names(.))])) > 0) %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
  
  # a placeholder so no need to handle the exception of 
  # no `Z_log2_R` columns before the first `normMulGau`
  if (purrr::is_empty(grep("^Z_log2_R[0-9]{3}[NC]{0,1}", names(df)))) {
    df <- df %>% 
      dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
      `names<-`(gsub("^N_log2_R", "Z_log2_R", names(.))) %>% 
      dplyr::bind_cols(df, .)
  }
  
  df <- normMulGau(
    df = df,
    method_align = "MC",
    n_comp = 1L,
    range_log2r = c(0, 100),
    range_int = c(0, 100),
    filepath = file.path(dat_dir, "Peptide/Histogram"),
    col_select = rlang::expr(Sample_ID), 
    cut_points = cut_points, 
  )
}


#' Summary of peptide keys by mean or geomean
#' 
#' @inheritParams info_anal
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
med_summarise_keys <- function(df, id) {
  ## --- Mascot ---
  mascot_median_keys <- c("pep_score", "pep_rank", "pep_isbold", "pep_exp_mr", "pep_delta", 
                          "pep_exp_mz", "pep_exp_z", "pep_locprob", "pep_locdiff")
  mascot_geomean_keys <- c("pep_expect")
  
  df_mascot_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mascot_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mascot_geomean <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mascot_geomean_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ my_geomean(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(mascot_median_keys, mascot_geomean_keys)))
  
  ## --- MaxQuant ---
  df_mq_rptr_mass_dev <- df %>% 
    dplyr::select(!!rlang::sym(id), grep("^Reporter mass deviation", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep("^Reporter mass deviation", names(.)))	  
  
  mq_median_keys <- c(
    "Score", 
    "Charge", "Mass", "PIF", "Fraction of total spectrum", "Mass error [ppm]", 
    "Mass error [Da]", "Base peak fraction", "Precursor Intensity", 
    "Precursor Apex Fraction", "Intensity coverage", "Peak coverage", 
    "Combinatorics"
  )
  mq_geomean_keys <- c("PEP")
  
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mq_geomean <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_geomean_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ my_geomean(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(mq_median_keys, mq_geomean_keys)))
  
  if ("Length" %in% names(df)) {
    df <- df %>% dplyr::mutate(pep_len = Length)
  }
  if ("Missed cleavages" %in% names(df)) {
    df <- df %>% dplyr::mutate(pep_miss = `Missed cleavages`)
  }
  if ("Missed Cleavages" %in% names(df)) {
    df <- df %>% dplyr::mutate(pep_miss = `Missed cleavages`)
  }
  
  df <- suppressWarnings(
    df %>% 
      dplyr::select(-one_of("Scan number", "Scan index", "Length", 
                            "Missed cleavages", "Missed Cleavages"))
  )
  
  if (all(c("Modifications", "Retention time") %in% names(df))) {
    df <- df %>% dplyr::select(!`Modifications`:`Retention time`)
  }

  df <- suppressWarnings(
    df %>% dplyr::select(-one_of("Delta score", "Score diff", "Localization prob", 
                                 "Precursor full scan number", "Precursor apex fraction",
                                 "Precursor apex offset", "Precursor apex offset time", 
                                 "Matches", 
                                 "Mass deviations [ppm]", "Masses", "Number of matches", 
                                 "Neutral loss level", "Intensities", "Mass deviations [Da]", 
                                 "ETD identification type", "Reverse", "All scores", 
                                 "All sequences", 
                                 "All modified sequences", 
                                 "Reporter PIF", "Reporter fraction", 
                                 "ID", # "Protein group IDs", 
                                 "Peptide ID", "Mod. peptide ID", 
                                 "Evidence ID"))
  ) %>% 
    dplyr::select(-grep("site\\s+IDs$", names(.)))
  
  ## --- Spectrum Mill ---
  sm_median_keys <- c(
    "score", "parent_charge", 
    "deltaForwardReverseScore", "percent_scored_peak_intensity", "totalIntensity", 
    "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
    "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
    "delta_parent_mass", "delta_parent_mass_ppm")
  sm_geomean_keys <- NA
  
  df_sm_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% sm_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% sm_median_keys))
  
  ## --- MSFragger ---
  mf_median_keys <- c("Nextscore", "PeptideProphet Probability")
  mf_geomean_keys <- NA
  
  df_mf_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mf_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% mf_median_keys)) 

  ## --- put together ---
  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id)))
  
  df <- list(df_first, 
             df_mascot_med, df_mascot_geomean, 
             df_mq_rptr_mass_dev, df_mq_med, df_mq_geomean, 
             df_sm_med, df_mf_med) %>%
    purrr::reduce(dplyr::left_join, by = id) %>%
    data.frame(check.names = FALSE)
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  ) 
}


#' load prior Peptide.txt
#' 
#' @inheritParams info_anal
load_prior <- function(filename, id) {
  dat_dir <- get_gl_dat_dir()
  
  stopifnot(file.exists(filename))
  
  df <- read.csv(filename, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0) 
  
  if (! id %in% names(df)) {
    try(unlink(file.path(dat_dir, "Peptide/Peptide.txt")))
    try(unlink(file.path(dat_dir, "Protein/Protein.txt")))
    stop("`Peptide.txt` deleted as column `", id, "` not available.", call. = FALSE)
  }
  
  df <- df %>% dplyr::arrange(!!rlang::sym(id))
}


#' format numeric columns
#' 
#' @inheritParams info_anal
fmt_num_cols <- function (df) {
  df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-  
    df[, grepl("^Z_log2_R[0-9]{3}", names(df))] %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    round(., digits = 3)
  
  df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-  
    df[, grepl("^N_log2_R[0-9]{3}", names(df))] %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    round(., digits = 3)
  
  df[, grepl("I[0-9]{3}", names(df))] <-  
    df[, grepl("I[0-9]{3}", names(df))] %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    round(., digits = 0)
  
  return(df)
}


#'Merge peptide table(s) into one
#'
#'\code{mergePep} merges individual peptide table(s),
#'\code{TMTset1_LCMSinj1_Peptide_N.txt, TMTset1_LCMSinj2_Peptide_N.txt} etc.,
#'into one interim \code{Peptide.txt}. The \code{log2FC} values in the interim
#'result are centered with the medians at zero (median centering). The utility
#'is typically applied after the conversion of PSMs to peptides via
#'\code{\link{PSM2Pep}} and is required even with a experiment at one multiplex
#'TMT and one LC/MS series.
#'
#'In the interim output file, "\code{Peptide.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means within
#'each plex if no \code{reference(s)} are present. Values under columns
#'\code{N_log2_R...} are median-centered \code{log2_R...} without scaling
#'normalization. Values under columns \code{Z_log2_R...} are \code{N_log2_R...}
#'with additional scaling normalization. Values under columns \code{I...} are
#'reporter-ion or LFQ intensity before normalization. Values under columns
#'\code{N_I...} are normalized \code{I...}. Values under columns
#'\code{sd_log2_R...} are the standard deviation of the \code{log2FC} of
#'proteins from ascribing peptides.
#'
#'Description of the column keys in the output: \cr 
#'\code{system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ")}
#'
#'The peptide counts in individual peptide tables,
#'\code{TMTset1_LCMSinj1_Peptide_N.txt} etc., may be fewer than the entries
#'indicated under the \code{prot_n_pep} column after the peptide
#'removals/cleanups using \code{purgePSM}.
#'
#'@param cut_points A numeric vector defines the cut points (knots) for the
#'  median-centering of \code{log2FC} by sections. The values of the knots
#'  indicate the summarized \code{log10(Intentisy)} of ions. For
#'  example, at \code{cut_points = seq(4, 7, .5)}, values of \code{log2FC} will
#'  be binned from \eqn{-Inf} to \eqn{Inf} according to the cut points at the
#'  reporter-ion (or LFQ) intensity of \eqn{10^4, 10^4.5, ... 10^7}. The default is
#'  \code{cut_points = Inf}, or equivalently \code{-Inf}, where the
#'  \code{log2FC} under each sample will be median-centered all together. See
#'  also \code{\link{prnHist}} for data binning and intensity-coded histograms.
#'@param use_duppeps Logical; if TRUE, re-assigns double/multiple dipping
#'  peptide sequences to the most likely proteins by majority votes.
#'@param omit_single_lfq Logical for MaxQuant LFQ; if TRUE, omits LFQ entries
#'  with single measured values across all samples. The default is TRUE.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  of data against the column keys in individual peptide tables of
#'  \code{TMTset1_LCMSinj1_Peptide_N.txt, TMTset1_LCMSinj2_Peptide_N.txt}, etc.
#'  \cr \cr The variable argument statements should be in the following format:
#'  each statement contains to a list of logical expression(s). The \code{lhs}
#'  needs to start with \code{filter_}. The logical condition(s) at the
#'  \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis. For
#'  example, \code{pep_len} is a column key present in \code{Mascot} peptide
#'  tables of \code{TMTset1_LCMSinj1_Peptide_N.txt},
#'  \code{TMTset1_LCMSinj2_Peptide_N.txt} etc. The statement
#'  \code{filter_peps_at = exprs(pep_len <= 50)} will remove peptide entries
#'  with \code{pep_len > 50}. See also \code{\link{normPSM}}.
#'@inheritParams normPSM
#'@inheritParams splitPSM
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
#'@return The primary output is in \code{.../Peptide/Peptide.txt}.
#'
#'@example inst/extdata/examples/mergePep_.R
#'@import stringr dplyr tidyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@importFrom plyr ddply
#'@export
mergePep <- function (plot_log2FC_cv = TRUE, use_duppeps = TRUE, cut_points = Inf, 
                      omit_single_lfq = TRUE, parallel = TRUE, ...) {
  dat_dir <- get_gl_dat_dir()  
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("mergePep"), 
          add = TRUE)
  
  stopifnot(vapply(c(plot_log2FC_cv), rlang::is_logical, logical(1)))

  reload_expts()

  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))

  dir.create(file.path(dat_dir, "Peptide/cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/Histogram"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/log2FC_cv/raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/log2FC_cv/purged"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/log"), recursive = TRUE, showWarnings = FALSE)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  filename <- file.path(dat_dir, "Peptide/Peptide.txt")
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  
  use_mq_pep <- use_mq_peptable(dat_dir, label_scheme_full)

  df <- normPep_Mplex(group_psm_by = group_psm_by, 
                      group_pep_by = group_pep_by, 
                      use_duppeps = use_duppeps, 
                      cut_points = cut_points, 
                      omit_single_lfq = omit_single_lfq,
                      parallel = parallel, 
                      use_mq_pep = use_mq_pep, 
                      !!!filter_dots) %T>% 
    write.table(filename, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  if (plot_log2FC_cv) {
    quiet_out <- purrr::quietly(sd_violin)(
      df = df, 
      id = !!group_pep_by, 
      filepath = file.path(dat_dir, "Peptide/log2FC_cv/raw/Peptide_sd.png"), 
      width = 8 * n_TMT_sets(label_scheme), 
      height = 8, 
      type = "log2_R", 
      adjSD = FALSE, 
      is_psm = FALSE
    )
  }
}


#'Standardize peptide results
#'
#'\code{standPep} standardizes peptide results from \code{\link{mergePep}} with
#'additional, stand-alone choices in data alignment. The utility is typically
#'applied after the assembly of peptide data via \code{\link{mergePep}}. It
#'further supports iterative normalization against data under selected sample
#'columns, data rows or both.
#'
#'
#'In the primary output file, "\code{Peptide.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means within
#'each plex if no \code{reference(s)} are present. Values under columns
#'\code{N_log2_R...} are aligned \code{log2_R...} according to
#'\code{method_align} without scaling normalization. Values under columns
#'\code{Z_log2_R...} are \code{N_log2_R...} with additional scaling
#'normalization. Values under columns \code{I...} are reporter-ion or LFQ
#'intensity before normalization. Values under columns \code{N_I...} are
#'normalized \code{I...}. Values under columns \code{sd_log2_R...} are the
#'standard deviation of the \code{log2FC} of proteins from ascribing peptides.
#'
#'In general, median statistics is applied when summarizing numeric peptide data
#'from different LCMS series. One exception is \code{pep_expect} where geometric
#'mean is used.
#'
#'Description of the column keys in the inputs and outputs: \cr
#'\code{system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ")}
#'
#'@param method_align Character string indicating the method in aligning
#'  \code{log2FC} across samples. \code{MC}: median-centering; \code{MGKernel}:
#'  the kernel density defined by multiple Gaussian functions
#'  (\code{\link[mixtools]{normalmixEM}}). At the \code{MC} default, the ratio
#'  profiles of each sample will be aligned in that the medians of the
#'  \code{log2FC} are zero. At \code{MGKernel}, the ratio profiles of each
#'  sample will be aligned in that the \code{log2FC} at the maximums of kernel
#'  density are zero.
#'@param col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  At the \code{NULL} default, the column key of \code{Select} in
#'  \code{expt_smry.xlsx} will be used. In the case of no samples being
#'  specified under \code{Select}, the column key of \code{Sample_ID} will be
#'  used. The non-empty entries under the ascribing column will be used in
#'  indicated analysis.
#'@param range_log2r Numeric vector at length two. The argument specifies the
#'  range of the \code{log2FC} for use in the scaling normalization of standard
#'  deviation across samples. The default is between the 10th and the 90th
#'  quantiles.
#'@param range_int Numeric vector at length two. The argument specifies the
#'  range of the \code{intensity} of reporter ions for use in the scaling
#'  normalization of standard deviation across samples. The default is between
#'  the 5th and the 95th quantiles.
#'@param n_comp Integer; the number of Gaussian components to be used with
#'  \code{method_align = MGKernel}. A typical value is 2 or 3. The variable
#'  \code{n_comp} overwrites the argument \code{k} in
#'  \code{\link[mixtools]{normalmixEM}}.
#'@param seed Integer; a seed for reproducible fitting at \code{method_align =
#'  MGKernel}.
#'@param ... \code{slice_}: variable argument statements for the identification
#'  of row subsets. The partial data will be taken for parameterizing the
#'  alignment of \code{log2FC} across samples. The full data set will be updated
#'  subsequently with the newly derived parameters. Note that there is no data
#'  entry removals from the complete data set with the \code{slice_} procedure.
#'  \cr \cr The variable argument statements should be in the following format:
#'  each of the statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{slice_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key present in \code{Peptide.txt}.
#'  The \code{slice_peps_at = exprs(pep_len >= 10, pep_len <= 50)} will extract
#'  peptide entries with the number of amino acid residues betwen 10 and 50 for
#'  \code{log2FC} alignment. Shorter or longer peptide sequences will remain in
#'  \code{Peptide.txt} but not used in the parameterization. See also
#'  \code{\link{normPSM}} for the variable arguments of \code{filter_}. \cr \cr
#'  Additional parameters from \code{\link[mixtools]{normalmixEM}}, i.e., \cr
#'  \code{maxit}, the maximum number of iterations allowed; \cr \code{epsilon},
#'  tolerance limit for declaring algorithm convergence.
#'@inheritParams normPSM
#'@seealso \emph{Metadata} \cr \code{\link{load_expts}} for metadata preparation
#'and a reduced working example in data normalization \cr
#'
#'\emph{Data normalization} \cr \code{\link{normPSM}} for extended examples in
#'PSM data normalization \cr \code{\link{PSM2Pep}} for extended examples in PSM
#'to peptide summarization \cr \code{\link{mergePep}} for extended examples in
#'peptide data merging \cr \code{\link{standPep}} for extended examples in
#'peptide data normalization \cr \code{\link{Pep2Prn}} for extended examples in
#'peptide to protein summarization \cr \code{\link{standPrn}} for extended
#'examples in protein data normalization. \cr \code{\link{purgePSM}} and
#'\code{\link{purgePep}} for extended examples in data purging \cr
#'\code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in
#'histogram visualization. \cr \code{\link{extract_raws}} and
#'\code{\link{extract_psm_raws}} for extracting MS file names \cr
#'
#'\emph{Variable arguments of `filter_...`} \cr \code{\link{contain_str}},
#'\code{\link{contain_chars_in}}, \code{\link{not_contain_str}},
#'\code{\link{not_contain_chars_in}}, \code{\link{start_with_str}},
#'\code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and
#'\code{\link{ends_with_chars_in}} for data subsetting by character strings \cr
#'
#'\emph{Missing values} \cr \code{\link{pepImp}} and \code{\link{prnImp}} for
#'missing value imputation \cr
#'
#'\emph{Informatics} \cr \code{\link{pepSig}} and \code{\link{prnSig}} for
#'significance tests \cr \code{\link{pepVol}} and \code{\link{prnVol}} for
#'volcano plot visualization \cr \code{\link{prnGSPA}} for gene set enrichment
#'analysis by protein significance pVals \cr \code{\link{gspaMap}} for mapping
#'GSPA to volcano plot visualization \cr \code{\link{prnGSPAHM}} for heat map
#'and network visualization of GSPA results \cr \code{\link{prnGSVA}} for gene
#'set variance analysis \cr \code{\link{prnGSEA}} for data preparation for
#'online GSEA. \cr \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS
#'visualization \cr \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA
#'visualization \cr \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA
#'visualization \cr \code{\link{pepHM}} and \code{\link{prnHM}} for heat map
#'visualization \cr \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'\code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for correlation
#'plots \cr \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for
#'trend analysis and visualization \cr \code{\link{anal_pepNMF}},
#'\code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'\code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'\code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for NMF analysis
#'and visualization \cr
#'
#'\emph{Custom databases} \cr \code{\link{Uni2Entrez}} for lookups between
#'UniProt accessions and Entrez IDs \cr \code{\link{Ref2Entrez}} for lookups
#'among RefSeq accessions, gene names and Entrez IDs \cr \code{\link{prepGO}}
#'for
#'\code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'ontology}} \cr \code{\link{prepMSig}} for
#'\href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular
#'signatures} \cr \code{\link{prepString}} and \code{\link{anal_prnString}} for
#'STRING-DB \cr
#'
#'\emph{Column keys in PSM, peptide and protein outputs} \cr 
#'system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output is in \code{.../Peptide/Peptide.txt}.
#'
#'@example inst/extdata/examples/normPep_.R
#'
#'@import stringr dplyr tidyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@importFrom plyr ddply
#'@export
standPep <- function (method_align = c("MC", "MGKernel"), col_select = NULL, range_log2r = c(10, 90), 
                      range_int = c(5, 95), n_comp = NULL, seed = NULL, plot_log2FC_cv = FALSE, ...) {

  dat_dir <- get_gl_dat_dir()
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("standPep"), add = TRUE)

  reload_expts()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  ok_existing_params(file.path(dat_dir, "Peptide/Histogram/MGKernel_params_N.txt"))

  filename <- file.path(dat_dir, "Peptide/Peptide.txt")
  if (!file.exists(filename)) {
    stop(filename, " not found; run `mergePep(...)` first.", call. = FALSE)
  }
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  method_align <- rlang::enexpr(method_align)
  if (method_align == rlang::expr(c("MC", "MGKernel"))) {
    method_align <- "MC"
  } else {
    method_align <- rlang::as_string(method_align)
    stopifnot(method_align %in% c("MC", "MGKernel"), 
              length(method_align) == 1)
  }
  
  col_select <- rlang::enexpr(col_select)
  col_select <- ifelse(is.null(col_select), rlang::expr(Sample_ID), rlang::sym(col_select))
  
  if (is.null(label_scheme[[col_select]])) {
    col_select <- rlang::expr(Sample_ID)
    warning("Column \'", rlang::as_string(col_select), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    col_select <- rlang::expr(Sample_ID)
    warning("No samples were specified under column \'", rlang::as_string(col_select), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  }
  
  local({
    vars <- list(range_int, range_log2r)
    stopifnot(vapply(vars, function (x) length(x) == 2, logical(1)))
    
    ranges <- vapply(vars, range, c(min = 0, max = 100))
    stopifnot(ranges[2, ] > ranges[1, ], 
              ranges[1, ] > 0,
              ranges[2, ] < 100) 
  })

  if (range_log2r[2] <= 1) range_log2r <- range_log2r * 100
  if (range_int[2] <= 1) range_int <- range_int * 100
  
  stopifnot(vapply(c(plot_log2FC_cv), rlang::is_logical, logical(1)))
  
  dots <- rlang::enexprs(...)
  
  message("Primary column keys in `Peptide/Peptide.txt` for `slice_` varargs.")

  df <- load_prior(filename, group_psm_by) %>% 
    normMulGau(
      df = .,
      method_align = method_align,
      n_comp = n_comp,
      seed = seed,
      range_log2r = range_log2r,
      range_int = range_int,
      filepath = file.path(dat_dir, "Peptide/Histogram"),
      col_select = col_select, 
      cut_points = Inf, 
      !!!dots, 
    ) %>% 
    fmt_num_cols() %T>% 
    write.table(file.path(dat_dir, "Peptide", "Peptide.txt"), 
                sep = "\t", col.names = TRUE, row.names = FALSE)

  if (plot_log2FC_cv & TMT_plex(label_scheme) > 0) {
    sd_violin(df = df, id = !!group_pep_by, 
              filepath = file.path(dat_dir, "Peptide/log2FC_cv/raw", "Peptide_sd.png"), 
              width = 8 * n_TMT_sets(label_scheme), height = 8, 
              type = "log2_R", adjSD = FALSE, is_psm = FALSE)
  }
}



#'Interim protein data
#'
#'\code{Pep2Prn} summarizes \code{Peptide.txt} to an interim protein report in
#'\code{Protein.txt}.
#'
#'Fields other than \code{log2FC} and \code{intensity} are summarized with
#'median statistics.
#'
#' @param method_pep_prn Character string; the method to summarize the
#'   \code{log2FC} and the \code{intensity} of peptides by protein entries. The
#'   descriptive statistics includes \code{c("mean", "median", "weighted_mean",
#'   "top_3_mean", "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum", "lfq_all")} with
#'   \code{median} being the default for TMT and \code{lfq_top_3_sum} for LFQ.
#'   The representative \code{log10-intensity} of reporter (or LFQ) ions at the
#'   peptide levels will be the weight when summarizing \code{log2FC} with
#'   various \code{"top_n"} statistics or \code{"weighted_mean"}.
#' @param use_unique_pep Logical. If TRUE, only entries that are \code{TRUE} or
#'   equal to \code{1} under the column \code{pep_isunique} in
#'   \code{Peptide.txt} will be used, for summarizing the \code{log2FC} and the
#'   \code{intensity} of peptides into protein values. The default is to use
#'   unique peptides only. For \code{MaxQuant} data, the levels of uniqueness
#'   are according to the \code{pep_unique_by} in \code{\link{normPSM}}. The
#'   argument currently do nothing to \code{Spectrum Mill} data where both
#'   unique and shared peptides will be kept.
#' @param ... \code{filter_}: Variable argument statements for the filtration of
#'   data rows. Each statement contains a list of logical expression(s). The
#'   \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'   the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'   For example, \code{pep_len} is a column key in \code{Peptide.txt}. The
#'   statement of \code{filter_peps_at = exprs(pep_len <= 50)} will remove
#'   peptide entries with \code{pep_len > 50}.
#' @inheritParams mergePep
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
#'@return The primary output in "\code{.../Protein/Protein.txt}".
#'
#'@example inst/extdata/examples/Pep2Prn_.R
#'@import stringr dplyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
Pep2Prn <- function (method_pep_prn = c("median", "mean", "weighted_mean", "top_3_mean", 
                                        "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum", "lfq_all"), 
                     use_unique_pep = TRUE, cut_points = Inf, ...) {
  
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Protein/Histogram"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/log"), recursive = TRUE, showWarnings = FALSE)
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("Pep2Prn"), add = TRUE)
  
  reload_expts()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  method_pep_prn <- rlang::enexpr(method_pep_prn)
  if (TMT_plex > 0) {
    if (length(method_pep_prn) > 1) {
      method_pep_prn <- "median"
    } else {
      method_pep_prn <- rlang::as_string(method_pep_prn)
    }
  } else {
    if (length(method_pep_prn) > 1) {
      method_pep_prn <- "lfq_top_3_sum"
    } else {
      method_pep_prn <- rlang::as_string(method_pep_prn)
    }
  }
  
  if (method_pep_prn == "top.3") {
    stop("Method `top.3` depreciated; instead use `top_3_mean`.", 
         call. = FALSE)
  } else if (method_pep_prn == "weighted.mean") {
    stop("Method `weighted.mean` depreciated; instead use `weighted_mean`.", 
         call. = FALSE)
  }
  
  stopifnot(method_pep_prn %in% c("median", "mean", 
                                  "weighted_mean", "top_3_mean", 
                                  "lfq_max", "lfq_top_2_sum", 
                                  "lfq_top_3_sum", "lfq_all"), 
            length(method_pep_prn) == 1)

  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  stopifnot(group_pep_by %in% c("prot_acc", "gene"), length(group_pep_by) == 1)
  
  stopifnot(vapply(c(use_unique_pep), rlang::is_logical, logical(1)))
  
  if (group_pep_by == "gene") {
    gn_rollup <- TRUE
    group_pep_by <- "prot_acc"
  } else {
    gn_rollup <- FALSE
  }
  
  stopifnot(group_pep_by == "prot_acc") 
  
  message("Primary column keys in `Peptide/Peptide.txt` ", 
          "for `filter_` varargs.")
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  df <- pep_to_prn(!!group_pep_by, 
                   method_pep_prn, 
                   use_unique_pep, 
                   gn_rollup, 
                   !!!filter_dots) 
  
  df <- normMulGau(
    df = df,
    method_align = "MC",
    n_comp = 1L,
    range_log2r = c(0, 100),
    range_int = c(0, 100),
    filepath = file.path(dat_dir, "Protein/Histogram"),
    col_select = rlang::expr(Sample_ID), 
    cut_points = cut_points, 
  ) 
  
  df <- df %>% 
    dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
                     as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), 
                     ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), 
                     as.numeric) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), 
                     ~ round(.x, digits = 3)) %T>% 
    write.table(file.path(dat_dir, "Protein/Protein.txt"), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
}


#' Helper: maps peptides to possible proteins
#' 
#' @param df A data frame.
#' @inheritParams normPSM
map_peps_prots <- function (df, group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                            type = "mf") {
  dat_dir <- get_gl_dat_dir()
  
  df <- df %>% dplyr::filter(!duplicated(.[[group_psm_by]]))
  
  if (type == "mf") {
    prot_key <- "Protein"
    mapped_prot_key <- "Mapped Proteins"
    gene_key <- "Gene"
    mapped_gene_key <- "Mapped Genes"
    sep = ", "
  } else if (type == "mq") {
    prot_key <- NULL
    mapped_prot_key <- "Proteins"
    gene_key <- NULL
    mapped_gene_key <- "gene"
    sep = ";"
  }
  
  if (group_pep_by == "gene") {
    x <- df %>% dplyr::rename(pri_col = gene_key, sec_col = mapped_gene_key)
  } else {
    x <- df %>% dplyr::rename(pri_col = prot_key, sec_col = mapped_prot_key)
  }
  
  if (type == "mf") {
    x <- x %>% dplyr::select(pri_col, sec_col)
  } else if (type == "mq") {
    x <- x %>% 
      dplyr::mutate(pri_col = NA) %>% 
      dplyr::select(pri_col, sec_col)
  }
  
  stopifnot(nrow(x) == nrow(df))
  
  x <- x %>% 
    dplyr::select(sec_col) %$% 
    stringr::str_split(.$sec_col, sep) %>% 
    purrr::map2(as.list(x$pri_col), ., ~ c(.x, .y) %>% .[!is.na(.)] %>% unique()) %>% 
    `names<-`(df[[group_psm_by]]) 
  
  pep_prot_map <- x 
  save(pep_prot_map, file = file.path(dat_dir, "pep_prot_map.rda"))
  
  invisible(x)
}


#' Helper: finds the row indexes of a protein within a family
#' 
#' Not currently used. The rows include primary (razor) and secondary (mapped)
#' proteins.
#' 
#' @param df A data frame.
#' @inheritParams normPSM
find_prot_family_rows <- function (df, group_psm_by, group_pep_by) {
  dat_dir <- get_gl_dat_dir()
  
  pep_prot_map <- map_peps_prots(df, group_psm_by, group_pep_by) 
  
  prots <- unique(df[[group_pep_by]])
  
  message("Find peptides unique to or shared by proteins; please wait...")
  
  prot_family_rows <- purrr::map(prots, ~ {
    prot <- .x
    prot_rows <- purrr::map_lgl(pep_prot_map, ~ prot %in% .x) %>% 
      which()
  }) %>% 
    `names<-`(prots)
  
  save(prot_family_rows, file = file.path(dat_dir, "prot_family_rows.rda"))
  invisible(prot_family_rows)
}


#' Calculates the \code{total}, \code{razor} and \code{unique} intensities of
#' proteins
#' 
#' @param df A data frame.
#' @param n Numeric; the top \code{n} entries to be used. The default is 3. All
#'   entries will be used at \eqn{n = Inf}.
#' @inheritParams normPSM
calc_lfq_prnnums <- function (df, use_unique_pep, group_psm_by, 
                              group_pep_by, method_pep_prn, type = "mf") {
  my_sum_n <- function (x, n = 3, ...) {
    if (n < 1) stop("`n` need to be a positive integer.", call. = FALSE)
    if (is.infinite(n)) n <- length(x)
    x %>% sort(decreasing = TRUE) %>% .[1:n] %>% sum(...)
  }
  
  
  load(file.path(dat_dir, "label_scheme.rda"))
  
  refChannels <- label_scheme %>% 
    dplyr::filter(Reference) %>% 
    dplyr::select(Sample_ID) %>% 
    unlist()
  
  pep_prot_map <- df %>% 
    map_peps_prots(group_psm_by, group_pep_by, type) %>% 
    plyr::ldply(rbind) %>% 
    `names_pos<-`(1, group_psm_by) %>% 
    tidyr::gather(-group_psm_by, key = n, value = !!rlang::sym(group_pep_by)) %>% 
    dplyr::select(-n) %>% 
    dplyr::filter(!is.na(!!rlang::sym(group_pep_by))) %>% 
    dplyr::rename(prot_map = group_pep_by)
  
  if (grepl("^lfq_top_", method_pep_prn)) {
    n <- gsub("^lfq_top_(\\d+)_[A-z]+", "\\1", method_pep_prn) %>% as.integer()
  } else if (method_pep_prn == "lfq_max") {
    n <- 1
  } else if (method_pep_prn == "lfq_all") {
    n <- Inf
  } else {
    n <- 3
  }
  
  # --- top_n ---
  # 1. summarizes top_n  of tot, razor and unique 
  # 2. selects one of them according to `pep_unique_by`
  three_lfqints <- local({
    # currently if top_n, not to filter by `use_unique_pep`
    #   to keep them the same as MSFragger calculations
    # instead let `pep_unique_by` decides which one to use
    df <- df %>% dplyr::left_join(pep_prot_map, by = group_psm_by) 

    prot_tot_ints <- df %>% 
      dplyr::select(prot_map, grep("^pep_tot_int\\s\\(", names(.))) %>% 
      dplyr::group_by(prot_map) %>% 
      dplyr::summarise_all(~ my_sum_n(.x, n = n, na.rm = TRUE))
    
    prot_razor_ints <- df %>% 
      dplyr::select(prot_map, grep("^pep_razor_int\\s\\(", names(.))) %>% 
      dplyr::group_by(prot_map) %>% 
      dplyr::summarise_all(~ my_sum_n(.x, n = n, na.rm = TRUE))
    
    prot_unique_ints <- df %>% 
      dplyr::select(prot_map, grep("^pep_unique_int\\s\\(", names(.))) %>% 
      dplyr::group_by(prot_map) %>% 
      dplyr::summarise_all(~ my_sum_n(.x, n = n, na.rm = TRUE))
    
    list(prot_tot_ints, prot_razor_ints, prot_unique_ints) %>% 
      purrr::reduce(dplyr::left_join, by = "prot_map") %>% 
      `names<-`(gsub("^pep_", "prot_", names(.))) 
  })

  pep_unique_by <- match_call_arg(normPSM, pep_unique_by)
  if (use_unique_pep) {
    if (pep_unique_by == "group") {
      dfw <- three_lfqints %>% 
        dplyr::select(grep("^prot_razor_int\\s\\(", names(.)))
    } else if (pep_unique_by == "protein") {
      dfw <- three_lfqints %>% 
        dplyr::select(grep("^prot_unique_int\\s\\(", names(.)))
    } else {
      stop("`pep_unique_by` need to be `group` or `protein`.", 
           call. = FALSE)
    }
  } else {
    dfw <- three_lfqints %>% 
      dplyr::select(grep("^prot_tot_int\\s\\(", names(.)))
  }
  
  dfw <- dfw %>% 
    `names<-`(gsub("^prot_.*_int", "I000", names(.))) %>% 
    calc_lfq_log2r("I000", refChannels) 
  
  # --- median centering ---
  dfw <- local({
    col_log2Ratio <- grepl("^log2_R[0-9]{3}[NC]{0,1}", names(dfw))
    cf <- apply(dfw[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
    
    dfw <- sweep(dfw[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
      `colnames<-`(paste("N", names(.), sep="_"))	%>%
      cbind(dfw, .)
    
    dfw <- sweep(dfw[, grepl("^I[0-9]{3}", names(dfw)), drop = FALSE], 2, 2^cf, "/") %>%
      `colnames<-`(paste("N", names(.), sep="_"))	%>%
      cbind(dfw, .)
    
    dfw <- dfw %>% dplyr::mutate(!!group_pep_by := three_lfqints[["prot_map"]])
    
    if (purrr::is_empty(grep("log2_R000", names(dfw)))) {
      stop("No `log2_R000...` columns available.\n", call. = FALSE)
    }
    
    dfw <- dfw %>% 
      dplyr::mutate_at(.vars = grep("log2_R000\\s", names(.)), 
                       ~ replace(.x, is.infinite(.), NA))
    
    if (purrr::is_empty(grep("^Z_log2_R[0-9]{3}[NC]{0,1}", names(dfw)))) {
      dfw <- dfw %>% 
        dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
        `names<-`(gsub("^N_log2_R", "Z_log2_R", names(.))) %>% 
        dplyr::bind_cols(dfw, .)
    }
    
    dplyr::bind_cols(
      dfw %>% dplyr::select(group_pep_by),
      dfw %>% dplyr::select(grep("^I000 \\(", names(.))),
      dfw %>% dplyr::select(grep("^N_I000 \\(", names(.))),
      dfw %>% dplyr::select(grep("^log2_R000 \\(", names(.))),
      dfw %>% dplyr::select(grep("^N_log2_R000 \\(", names(.))),
      dfw %>% dplyr::select(grep("^Z_log2_R000 \\(", names(.))),
    )
  })
}


#' Helper: calculates the TMT log2FC and reporter-ion intensity of proteins
#' 
#' @param id Always "prot_acc".
#' @inheritParams calc_lfq_prnnums
#' @inheritParams Pep2Prn
calc_tmt_prnnums <- function (df, use_unique_pep, id = "prot_acc", method_pep_prn) {
  if (use_unique_pep && "pep_isunique" %in% names(df)) {
    df <- df %>% dplyr::filter(pep_isunique == 1)
  }
  
  df_num <- df %>% 
    dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    dplyr::select(-grep("^sd_log2_R[0-9]{3}", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id))
  
  df_num <- switch(method_pep_prn, 
                   mean = aggrNums(mean)(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   median = aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE),
                   top_3_mean = TMT_top_n(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   weighted_mean = tmt_wtmean(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE))
  
  run_scripts <- FALSE
  if (run_scripts) {
    df_num <- local({
      prn_tot_int <- df_num %>% 
        dplyr::select(grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
        `names<-`(gsub("^log2_R[0-9]{3}[NC]{0,1}", "prot_tot_int", names(.))) 
      prn_tot_int[] <- NA
      
      prn_razor_int <- prn_tot_int %>% 
        `names<-`(gsub("prot_tot_in", "prn_razor_int", names(.)))
      
      prn_unique_int <- prn_tot_int %>% 
        `names<-`(gsub("prot_tot_in", "prn_unique_int", names(.)))
      
      prn_ints <- dplyr::bind_cols(
        prn_tot_int,
        prn_razor_int,
        prn_unique_int, 
        df_num,
      )
    })
  }
  
  invisible(df_num)
}


#' Helper of Pep2Prn
#' 
#' @param gn_rollup Logical; if TRUE, rolls up protein accessions to gene names.
#' @inheritParams info_anal
#' @inheritParams Pep2Prn
pep_to_prn <- function(id, method_pep_prn, use_unique_pep, gn_rollup, ...) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  # `id` is always "prot_acc"; `group_pep_by` could be "gene"
  id <- rlang::as_string(rlang::enexpr(id))
  
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  df <- read.csv(file.path(dat_dir, "Peptide/Peptide.txt"), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
  
  if (! "pep_isunique" %in% names(df)) {
    df$pep_isunique <- TRUE
    warning("Column `pep_isunique` created and TRUE values assumed.", 
            call. = FALSE)
  } else if (all(is.na(df$pep_isunique))) {
    df$pep_isunique <- TRUE
    warning("Values of `pep_isunique` are all NA and coerced to TRUE.", 
            call. = FALSE)
  }

  df <- df %>% filters_in_call(!!!filter_dots)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  df <- local({
    df_shared <- df %>% 
      dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by), 
                    pep_n_psm, prot_n_psm, prot_n_pep, pep_isunique) %>% 
      dplyr::filter(!pep_isunique)
    
    prot_n_sharepeps <- df_shared %>% 
      dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>% 
      dplyr::summarise(prot_n_sharepeps = n())
    
    prot_n_sharepsms <- df_shared %>% 
      dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by), pep_n_psm) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>% 
      dplyr::summarise(prot_n_sharepsms = sum(pep_n_psm))
    
    df <- list(df, prot_n_sharepeps, prot_n_sharepsms) %>% 
      purrr::reduce(dplyr::left_join, by = group_pep_by) %>% 
      dplyr::mutate(prot_n_sharepeps = replace(prot_n_sharepeps, is.na(prot_n_sharepeps), 0), 
                    prot_n_sharepsms = replace(prot_n_sharepsms, is.na(prot_n_sharepsms), 0)) %>% 
      dplyr::mutate(prot_n_uniqpep = prot_n_pep - prot_n_sharepeps, 
                    prot_n_uniqpsm = prot_n_psm - prot_n_sharepsms) %>% 
      dplyr::select(-prot_n_sharepeps, -prot_n_sharepsms)
    
    df %>% 
      ins_cols_after(which(names(.) == "prot_n_psm"), which(names(.) == "prot_n_uniqpsm")) %>% 
      ins_cols_after(which(names(.) == "prot_n_pep"), which(names(.) == "prot_n_uniqpep"))
  })
  
  # temporarily for special handling MaxQuant
  if ("Mapped Proteins" %in% names(df)) {
    type = "mf"
  } else if ("Proteins" %in% names(df)) {
    type = "mq"
  }
  
  if (grepl("^lfq_", method_pep_prn)) {
    df_num <- calc_lfq_prnnums(df, use_unique_pep, group_psm_by, 
                               "fasta_name", method_pep_prn, type) %>% 
      dplyr::left_join(df[, c("fasta_name", id)] %>% 
                         dplyr::filter(!duplicated(!!rlang::sym(id))),
                       by = "fasta_name") %>% 
      dplyr::select(-fasta_name) %>% 
      dplyr::filter(!is.na(!!rlang::sym(id)))
  } else {
    df_num <- calc_tmt_prnnums(df, use_unique_pep, id, method_pep_prn)
  }
  
  df_num <- local({
    count_nna <- df_num %>% 
      dplyr::select(grep("N_log2_R[0-9]{3}[NC]{0,1}", names(.)))%>% 
      dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Ref\\.[0-9]+\\)$", 
                          names(.))) %>% 
      dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Empty\\.[0-9]+\\)$", 
                          names(.))) %>% 
      is.na() %>% 
      `!`() %>% 
      rowSums()
    
    dplyr::bind_cols(count_nna = count_nna, df_num)
  })
  
  df <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("is_tryptic", "count_nna"))) %>% 
    dplyr::select(-grep("^Reporter mass deviation", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("m/z", "Charge", "Mass", "Mass error [ppm]", 
                                         "Mass error [Da]", "Score", "Combinatorics", 
                                         "PIF", "Fraction of total spectrum", 
                                         "Base peak fraction", 
                                         "Precursor Intensity", "Precursor intensity", 
                                         "Precursor Apex Fraction", 
                                         "Intensity coverage", "Intensity Coverage", 
                                         "Peak coverage", "Peak Coverage", "PEP", 
                                         "Proteins", "Protein group IDs"))) %>% 
    dplyr::select(-which(names(.) %in% c("Is Unique", "Protein", "Gene", "Mapped Genes", 
                                         "Mapped Proteins", "Nextscore", 
                                         "PeptideProphet Probability")))
  
  mq_median_keys <- NULL
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  df <- df %>% dplyr::select(-which(names(.) %in% mq_median_keys))
  rm(mq_median_keys)
  
  sm_median_keys <- c(
    "deltaForwardReverseScore", "percent_scored_peak_intensity", "totalIntensity", 
    "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
    "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
    "matched_parent_mass", "delta_parent_mass", "delta_parent_mass_ppm")
  df_sm_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% sm_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  df <- df %>% dplyr::select(-which(names(.) %in% sm_median_keys))
  rm(sm_median_keys)
  
  mf_median_keys <- NULL
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mf_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  df <- df %>% dplyr::select(-which(names(.) %in% mf_median_keys))
  rm(mf_median_keys)

  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
    dplyr::select(-grep("^pep_", names(.)))
  
  df <- list(df_first, 
             df_mq_med, 
             df_sm_med, 
             df_num) %>% 
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
    df$count_nna <- NULL
    
    dfa <- local({
      dfa <- df %>% 
        dplyr::select(gene, grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
        dplyr::filter(!is.na(gene)) %>% 
        dplyr::group_by(gene) %>% 
        dplyr::summarise_all(list(~ median(.x, na.rm = TRUE)))
      
      count_nna <- dfa %>% 
        dplyr::select(grep("N_log2_R[0-9]{3}[NC]{0,1}", 
                           names(.)))%>% 
        dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Ref\\.[0-9]+\\)$", 
                            names(.))) %>% 
        dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Empty\\.[0-9]+\\)$", 
                            names(.))) %>% 
        is.na() %>% 
        `!`() %>% 
        rowSums() 
      
      dplyr::bind_cols(count_nna = count_nna, dfa)
    })

    dfb <- df %>% 
      dplyr::select(-prot_cover, -grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
      dplyr::filter(!is.na(gene)) %>% 
      dplyr::filter(!duplicated(.$gene))
    
    dfc <- df %>% 
      dplyr::select(gene, prot_cover) %>% 
      dplyr::filter(!is.na(gene), !is.na(prot_cover)) %>% 
      dplyr::filter(prot_cover != "NA%") %>% 
      dplyr::group_by(gene) %>% 
      dplyr::mutate(prot_cover = as.numeric(sub("%", "", prot_cover))) %>% 
      dplyr::summarise_all(~ max(.x, na.rm = TRUE)) %>% 
      dplyr::mutate(prot_cover = paste0(prot_cover, "%"))
    
    df <- list(dfc, dfb, dfa) %>% 
      purrr::reduce(right_join, by = "gene") %>% 
      dplyr::filter(!is.na(gene), !duplicated(gene))
    
    df <- dplyr::bind_cols(
      df %>% dplyr::select(prot_acc), 
      df %>% dplyr::select(-prot_acc), 
    ) %>% 
      reloc_col_before("gene", "fasta_name") %>% 
      reloc_col_before("prot_cover", "prot_n_psm")
  } 
  
  return(df)
}


#' Assign duplicated peptides to a leading protein
#' @param df A PSM data frame
#' @inheritParams mergePep
#' @inheritParams annotPSM
assign_duppeps <- function(df, group_psm_by, group_pep_by, use_duppeps = TRUE) {
  # Scenario: 
  # In `dat_file_1`, peptide_x assigned to Prn_MOUSE against "human + mouse" databases.
  # In `dat_file_2` the same peptide_x assigned to PRN_HUMAN against "human only" database.
  # When combining, `dat_file_1` and `dat_file_2`, all the peptide entries will be 
  #   re-assigned to the protein id with the greater `prot_n_pep`.

  dat_dir <- get_gl_dat_dir()
  
  dup_peps <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
    dplyr::summarise(N = n_distinct(!!rlang::sym(group_pep_by))) %>%
    dplyr::filter(N > 1)
  
  if (nrow(dup_peps) > 0) {
    if (use_duppeps) {
      # to ensure the same order in column names during replacement
      col_nms <- suppressWarnings(
        df %>% 
          dplyr::select(grep("^prot_", names(.)), 
                        one_of(c("gene", "acc_type", "entrez", "species", 
                                 "kin_attr", "kin_class", "kin_order"))) %>% 
          dplyr::select(-one_of(c("prot_n_psm", "prot_n_pep", "prot_cover", 
                                  "prot_matches_sig", "prot_sequences_sig"))) %>% 
          names()
      )
      
      df_dups <- purrr::map(as.character(unique(dup_peps[[group_psm_by]])), ~ {
        df_sub <- df %>% 
          dplyr::filter(!!rlang::sym(group_psm_by) == .x) %>% 
          dplyr::arrange(-prot_n_pep, -prot_n_psm, -prot_mass)
  
        # if (group_pep_by prot_acc) use the first prot_acc and also the first gene
        # if (group_pep_by gene) use the first gene and also the first prot_acc
        
        cols_replace <- df_sub %>% 
          dplyr::select(col_nms) %>% 
          dplyr::slice(1)

        df_sub2 <- df_sub %>% dplyr::slice(-1)
        df_sub2[, col_nms] <- cols_replace
        
        df_sub <- dplyr::bind_rows(df_sub %>% dplyr::slice(1), df_sub2)
      }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::select(names(df)) 
      
      df <- dplyr::bind_rows(
        df %>% dplyr::filter(! (!!rlang::sym(group_psm_by) %in% df_dups[[group_psm_by]])), 
        df_dups) 
      
      # update `dup_peps`; should be empty
      dup_peps_af <- df %>% 
        dplyr::filter(!!rlang::sym(group_psm_by) %in% dup_peps[[group_psm_by]]) %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
        dplyr::summarise(N = n_distinct(!!rlang::sym(group_pep_by))) %>%
        dplyr::filter(N > 1)
      
      if (nrow(dup_peps_af) > 0) {
        write.csv(dup_peps_af, file.path(dat_dir, "Peptide/dbl_dipping_peptides.csv"), 
                  row.names = FALSE)
        df <- df %>% dplyr::filter(! (!!rlang::sym(group_psm_by) %in% dup_peps_af[[group_psm_by]]))
      }
    } else {
      write.csv(dup_peps, file.path(dat_dir, "Peptide/dbl_dipping_peptides.csv"), 
                row.names = FALSE)
      df <- df %>% dplyr::filter(! (!!rlang::sym(group_psm_by) %in% dup_peps[[group_psm_by]]))
    }
  }
  
  return(df)
}


