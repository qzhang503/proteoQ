#' Make new column names
#'
#' \code{newColnames} match names to Sample_ID in label_scheme
#'
#' @param i Integer; the index of TMT experiment 
#' @param x Data frame; log2FC data
#' @param label_scheme Experiment summary
#' @import dplyr purrr rlang forcats
#' @importFrom magrittr %>% %T>% %$% %<>% 
newColnames <- function(i, x, label_scheme) {
  label_scheme_sub <- label_scheme %>%
    dplyr::filter(TMT_Set == i) # %>% 
    # dplyr::arrange(TMT_Channel)
  
  cols <- grep(paste0("[RI][0-9]{3}[NC]{0,1}_", i, "$"), names(x))
  nm_channel <- gsub(paste0("([RI][0-9]{3}[NC]{0,1})_", i, "$"), "\\1", names(x)[cols])
  names(x)[cols] <- paste0(nm_channel, " (", as.character(label_scheme_sub$Sample_ID), ")")
  
  # update the column indexes with TMT_Set being replaced with Sample_ID
  cols <- grep("[RI][0-9]{3}[NC]{0,1}\\s+\\(.*\\)$", names(x))
  
  # cols with the updated names go first
  if (length(cols) < ncol(x)) x <- dplyr::bind_cols(x[, cols], x[, -cols, drop = FALSE])
  
  return(x)
}


#' Check the single presence of MaxQuant peptides[...].txt.
#' 
#' @inheritParams load_expts
single_mq_peptable <- function (dat_dir) {
  filelist <- list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  
  if (length(filelist) > 1) {
    stop("Only single MaxQuant `peptides.txt` allowed.", call. = FALSE)
  }
  
  if (!file.exists(file.path(dat_dir, filelist))) {
    stop("No MaxQuant LFQ file `peptides[...].txt` udner", dat_dir, call. = FALSE)
  }
  
  message("MaxQuant file `", filelist, "` found.")  
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
#' 
check_mq_df <- function (df, label_scheme) {
  mq_nms <- names(df) %>% 
    .[grepl("^LFQ intensity ", .)] %>% 
    gsub("^LFQ intensity ", "", .)
  
  ls_nms <- label_scheme$Sample_ID %>% 
    .[!grepl("^Empty\\.[0-9]+", .)] %>% 
    unique() 
  
  if (!all(mq_nms %in% ls_nms)) {
    missing_nms <- mq_nms %>% .[! . %in% ls_nms]

    # the same ID occurs in "Identification type", "Experiment", "Intensity", "LFQ intensity"
    df <- df %>% 
      dplyr::mutate_at(.vars = grep(paste0(" ", missing_nms, "$"), names(.)), ~ {.x <- NULL})
  }
  
  invisible(df)
}


#' Extracts MaxQuant intensity values and calculates log2FC.
#' 
#' @param df A data frame of MaxQuant results.
extract_mq_ints <- function (df) {
  calc_log2r <- function (df, type, refChannels) {
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
      `colnames<-`(gsub(paste0(type, "(.*)$"), paste0(prefix, " \\(", "\\1", "\\)"), names(.))) %>%
      cbind(df, .)
  }
  
  
  load(file.path(dat_dir, "label_scheme.rda"))
  
  refChannels <- label_scheme %>% 
    dplyr::filter(Reference) %>% 
    dplyr::select(Sample_ID) %>% 
    unlist()
  
  df <- df %>% 
    calc_log2r("Intensity", refChannels) %>% 
    calc_log2r("LFQ intensity", refChannels) %>% 
    dplyr::mutate_at(.vars = grep("log2_R000\\s", names(.)), ~ replace(.x, is.infinite(.), NA)) 
  
  log2sd <- df %>% 
    dplyr::select(grep("Intensity\\s", names(.))) %>% 
    `names<-`(gsub("^Intensity\\s(.*)", "sd_log2_R000 \\(\\1\\)", names(.))) %>% 
    dplyr::mutate_all(~ replace(.x, !is.na(.x), NA))
  
  df <- dplyr::bind_cols(df, log2sd)
  
  df <- df %>% 
    `names<-`(gsub("^Intensity (.*)$", paste0("I000 \\(", "\\1", "\\)"), names(.))) %>% 
    `names<-`(gsub("^LFQ intensity (.*)$", paste0("N_I000 \\(", "\\1", "\\)"), names(.)))
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[IR]{1}000 \\(", names(.))),
    df %>% dplyr::select(grep("^I000 \\(", names(.))),
    df %>% dplyr::select(grep("^N_I000 \\(", names(.))),
    df %>% dplyr::select(grep("^sd_log2_R000 \\(", names(.))),
    df %>% dplyr::select(grep("^log2_R000 \\(", names(.))),
    df %>% dplyr::select(grep("^N_log2_R000 \\(", names(.))),  
  )
}





extract_mq_ints_orig <- function (df) {
  df <- df %>% 
    dplyr::mutate(Mean_Int = rowMeans(.[, grepl("^Intensity\\s", names(.))], na.rm = TRUE))
  
  log2r <- df %>% 
    dplyr::select(grep("^Intensity\\s", names(.)), Mean_Int) %>% 
    dplyr::mutate_at(.vars = grep("^Intensity\\s", names(.)), 
                     ~ log2(.x/Mean_Int)) %>% 
    dplyr::select(-Mean_Int) %>% 
    `names<-`(gsub("Intensity\\s(.*)", "log2_R000 \\(\\1\\)", names(.))) %>% 
    dplyr::mutate_at(.vars = grep("^log2_R000\\s", names(.)), 
                     ~ replace(.x, is.infinite(.x), NA)) 
  
  df <- df %>% 
    dplyr::mutate(Mean_Int = rowMeans(.[, grepl("^LFQ intensity\\s", names(.))], na.rm = TRUE))
  
  n_log2r <- df %>% 
    dplyr::select(grep("^LFQ intensity\\s", names(.)), Mean_Int) %>% 
    dplyr::mutate_at(.vars = grep("^LFQ intensity\\s", names(.)), 
                     ~ log2(.x/Mean_Int)) %>% 
    dplyr::select(-Mean_Int) %>% 
    `names<-`(gsub("LFQ intensity\\s(.*)", "N_log2_R000 \\(\\1\\)", names(.))) %>% 
    dplyr::mutate_at(.vars = grep("^N_log2_R000\\s", names(.)), 
                     ~ replace(.x, is.infinite(.x), NA)) 
  
  log2sd <- log2r %>% 
    `names<-`(paste0("sd_", names(.))) %>% 
    dplyr::mutate_all(~ replace(.x, !is.na(.x), NA))
  
  df <- list(df, log2sd, log2r, n_log2r) %>% 
    dplyr::bind_cols() %>% 
    dplyr::select(-Mean_Int)
  
  df <- df %>% 
    `names<-`(gsub("^Intensity (.*)$", paste0("I000 \\(", "\\1", "\\)"), names(.))) %>% 
    `names<-`(gsub("^LFQ intensity (.*)$", paste0("N_I000 \\(", "\\1", "\\)"), names(.)))
}

#' load MaxQuant peptide table
#' 
#' @inheritParams n_TMT_sets
pep_mq_lfq <- function(label_scheme) {
  dat_dir <- get_gl_dat_dir()
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  
  if (purrr::is_empty(filelist)) {
    stop(paste("No MaxQuant LFQ file of `peptides[...].txt` under", file.path(dat_dir), ".\n",
               "Make sure that the name of file starts with `peptides`."), 
         call. = FALSE)
  }

  df <- read.csv(file.path(dat_dir, filelist), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(not_allzero_rows(.[grep("^LFQ intensity ", names(.))])) 
  
  df <- check_mq_df(df, label_scheme)

  # note: empty `Gene names` will not be coerced to `prot_acc`
  df <- df %>% 
    dplyr::rename(
      pep_seq = Sequence, 
      prot_acc = Proteins, 
    ) %>% 
    dplyr::mutate(gene = gsub("\\;.*", "", `Gene names`)) %>% 
    dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc))
  
  # `Modified sequence` not available in `peptides.txt` at version 1.6.15
  # group_psm_by = "pep_seq" only in normPSM
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
      dplyr::select(group_psm_by, grep("^Intensity |^LFQ intensity ", names(.))) %>% 
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

  # calculate log2FC to reference(s)

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
                           use_duppeps = TRUE, parallel = TRUE, ...) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)

  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  
  filelist <- list.files(path = file.path(dat_dir, "Peptide"), 
                         pattern = paste0("TMTset[0-9]+_LCMSinj[0-9]+_Peptide_N\\.txt$"), 
                         full.names = TRUE)
  
  if (purrr::is_empty(filelist)) stop("No individual peptide tables; run `PSM2Pep()` first.")

  df <- do.call(rbind, 
                purrr::map(filelist, read.csv, check.names = FALSE, header = TRUE, 
                           sep = "\t", comment.char = "#")) %>%
    dplyr::mutate(TMT_Set = factor(TMT_Set)) %>%
    dplyr::arrange(TMT_Set) 
  
  df <- df %>% filters_in_call(!!!filter_dots)
  
  if ("gene" %in% names(df)) 
    df <- df %>% dplyr::mutate(gene = forcats::fct_explicit_na(gene))

  df <- df %>% assign_duppeps(group_psm_by, group_pep_by, use_duppeps)
  
  if (TMT_plex > 0) {
    df_num <- calc_tmt_nums(df, filelist, group_psm_by, parallel)
  } else {
    df_num <- pep_mq_lfq(label_scheme)
  }

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
    dplyr::select(-pep_n_psm, -prot_n_psm, -prot_n_pep, -TMT_Set) %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))  

  df <- list(pep_n_psm, df_first, df_num) %>%
    purrr::reduce(dplyr::left_join, by = group_psm_by)
  
  df <- list(df, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = group_pep_by)
  
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
  }
  
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
  
  df <- normMulGau(
    df = df,
    method_align = "MC",
    n_comp = 1L,
    range_log2r = c(0, 100),
    range_int = c(0, 100),
    filepath = file.path(dat_dir, "Peptide/Histogram"),
    col_select = rlang::expr(Sample_ID), 
  )
}


#' Summary of peptide keys by mean or geomean
#' 
#' @inheritParams info_anal
#' @import dplyr purrr rlang
#' @importFrom magrittr %>% %T>% %$% %<>% 
med_summarise_keys <- function(df, id) {
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
  
  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id)))
  
  df <- list(df_first, 
             df_mascot_med, df_mascot_geomean, 
             df_mq_rptr_mass_dev, df_mq_med, df_mq_geomean, 
             df_sm_med) %>%
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
#'\code{\link{PSM2Pep}} and is required even for a experiment with one multiplex
#'TMT and one LCMS injection.
#'
#'In the interim output file, "\code{Peptide.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'median-centered \code{log2_R...} without scaling normalization. Values under
#'columns \code{Z_log2_R...} are \code{N_log2_R...} with additional scaling
#'normalization. Values under columns \code{I...} are \code{reporter-ion
#'intensity} before normalization. Values under columns \code{N_I...} are
#'normalized \code{I...}. Values under columns \code{sd_log2_R...} are the
#'standard deviation of the \code{log2FC} of proteins from ascribing peptides.
#'
#'Description of the column keys in the output: \cr \code{system.file("extdata",
#'"mascot_peptide_keys.txt", package = "proteoQ")} \cr
#'\code{system.file("extdata", "maxquant_peptide_keys.txt", package =
#'"proteoQ")}
#'
#'The peptide counts in individual peptide tables,
#'\code{TMTset1_LCMSinj1_Peptide_N.txt} etc., may be fewer than the entries
#'indicated under the \code{prot_n_pep} column after the peptide
#'removals/cleanups using \code{purgePSM}.
#'
#' @param use_duppeps Logical; if TRUE, re-assigns double/multiple dipping
#'   peptide sequences to the most likely proteins by majority votes.
#' @param ... \code{filter_}: Variable argument statements for the row filtration
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
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@return The primary output is in \code{.../Peptide/Peptide.txt}.
#'
#'@example inst/extdata/examples/mergePep_.R
#'@import stringr dplyr tidyr purrr rlang
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@importFrom plyr ddply
#'@export
mergePep <- function (plot_log2FC_cv = TRUE, use_duppeps = TRUE, parallel = TRUE, ...) {
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

  if (TMT_plex(label_scheme_full) > 0) {
    message("Primary column keys in `Peptide/TMTset1_LCMSinj1_Peptide_N.txt` etc. for `filter_` varargs.")
  } else {
    single_mq_peptable(dat_dir)
  }
  
  df <- normPep_Mplex(group_psm_by = group_psm_by, 
                      group_pep_by = group_pep_by, 
                      use_duppeps = use_duppeps, 
                      parallel = parallel, 
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
#'additional choices in data alignment. The utility is typically applied after
#'the assembly of peptide data via \code{\link{mergePep}}. It further supports
#'iterative normalization against data under selected sample columns, data rows
#'or both.
#'
#'
#'In the primary output file, "\code{Peptide.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'aligned \code{log2_R...} according to \code{method_align} without scaling
#'normalization. Values under columns \code{Z_log2_R...} are \code{N_log2_R...}
#'with additional scaling normalization. Values under columns \code{I...} are
#'\code{reporter-ion intensity} before normalization. Values under columns
#'\code{N_I...} are normalized \code{I...}. Values under columns
#'\code{sd_log2_R...} are the standard deviation of the \code{log2FC} of
#'proteins from ascribing peptides.
#'
#'In general, median statistics is applied when summarizing numeric peptide data
#'from different LCMS series. One exception is \code{pep_expect} with Mascot
#'workflow where geometric mean is used.
#'
#'Description of the column keys in the inputs and outputs: \cr
#'\code{system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ")}
#'\cr \code{system.file("extdata", "maxquant_peptide_keys.txt", package =
#'"proteoQ")}
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
#'@param seed Integer; a seed setting a starting point for reproducible
#'  analyses.
#'@param ... \code{slice_}: variable argument statements for the identification
#'  of row subsets. The partial data will be taken for parameterizing the
#'  alignment of \code{log2FC} across samples. The full data set will be updated
#'  subsequently with the newly derived paramters. Note that there is no data
#'  entry removals from the complete data set with the \code{slice_} procedure.
#'  \cr \cr The variable argument statements should be in the following format:
#'  each of the statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{slice_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key present in \code{Peptide.txt}
#'  with \code{Mascot} workflows. The \code{slice_peps_at = exprs(pep_len >= 10,
#'  pep_len <= 50)} will extract peptide entries with the number of amino acid
#'  residues betwen 10 and 50 for \code{log2FC} alignment. Shorter or longer
#'  peptide sequences will remain in \code{Peptide.txt} but not used in the
#'  parameterization. See also \code{\link{normPSM}} for the variable arguments
#'  of \code{filter_}. \cr \cr Additional parameters from
#'  \code{\link[mixtools]{normalmixEM}}, i.e., \cr \code{maxit}, the maximum
#'  number of iterations allowed; \cr \code{epsilon}, tolerance limit for
#'  declaring algorithm convergence.
#'@inheritParams normPSM
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
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output is in \code{.../Peptide/Peptide.txt}.
#'
#'@example inst/extdata/examples/normPep_.R
#'
#'@import stringr dplyr tidyr purrr rlang
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
  if (!file.exists(filename)) stop(filename, " not found; run `mergePep(...)` first.", call. = FALSE)
  
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
#'@param method_pep_prn Character string; the method to summarize the
#'  \code{log2FC} and the \code{intensity} of peptides by protein entries. The
#'  descriptive statistics includes \code{c("mean", "median", "top.3",
#'  "weighted.mean")} with \code{median} being the default. The representative
#'  \code{log10-intensity} of reporter ions at the peptide levels (from
#'  \code{\link{standPep}}) will be the weight when summarizing \code{log2FC}
#'  with \code{top.3} or \code{weighted.mean}.
#'@param use_unique_pep Logical. If TRUE, only entries that are \code{TRUE} or
#'  equal to \code{1} under the column \code{pep_isunique} in \code{Peptide.txt}
#'  will be used, for summarizing the \code{log2FC} and the \code{intensity} of
#'  peptides into protein values. The default is to use unique peptides only.
#'  For \code{MaxQuant} data, the levels of uniqueness are according to the
#'  \code{pep_unique_by} in \code{\link{normPSM}}. The argument currently do
#'  nothing to \code{Spectrum Mill} data where both unique and shared peptides
#'  will be kept.
#'@param ... \code{filter_}: Variable argument statements for the filtration of
#'  data rows. Each statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key present in \code{Peptide.txt}
#'  with \code{Mascot} workflows. The statement of \code{filter_peps_at =
#'  exprs(pep_len <= 50)} will remove peptide entries with \code{pep_len > 50}.
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
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@return The primary output in "\code{.../Protein/Protein.txt}".
#'
#'@example inst/extdata/examples/Pep2Prn_.R
#'@import stringr dplyr purrr rlang
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
Pep2Prn <- function (method_pep_prn = c("median", "mean", "weighted.mean", "top.3"), 
                     use_unique_pep = TRUE, ...) {
  
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Protein/Histogram"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/log"), recursive = TRUE, showWarnings = FALSE)
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("Pep2Prn"), add = TRUE)
  
  reload_expts()
  
  method_pep_prn <- rlang::enexpr(method_pep_prn)
  if (method_pep_prn == rlang::expr(c("median", "mean", "weighted.mean", "top.3"))) {
    method_pep_prn <- "median"
  } else {
    method_pep_prn <- rlang::as_string(method_pep_prn)
    stopifnot(method_pep_prn %in% c("mean", "top.3", "median", "weighted.mean"), 
              length(method_pep_prn) == 1)
  }
  
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
  
  message("Primary column keys in `Peptide/Peptide.txt` for `filter_` varargs.")
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  
  df <- pep_to_prn(!!group_pep_by, method_pep_prn, use_unique_pep, gn_rollup, !!!filter_dots) 
  
  df <- normMulGau(
    df = df,
    method_align = "MC",
    n_comp = 1L,
    range_log2r = c(0, 100),
    range_int = c(0, 100),
    filepath = file.path(dat_dir, "Protein/Histogram"),
    col_select = rlang::expr(Sample_ID), 
  ) 
  
  df <- df %>% 
    dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %T>% 
    write.table(file.path(dat_dir, "Protein/Protein.txt"), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
}


#' Helper of Pep2Prn
#' 
#' @param gn_rollup Logical; if TRUE, rolls up protein accessions to gene names.
#' @inheritParams info_anal
#' @inheritParams Pep2Prn
#' @rawNamespace import(seqinr, except = c(consensus, count, zscore))
pep_to_prn <- function(id, method_pep_prn, use_unique_pep, gn_rollup, ...) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme.rda"))
  id <- rlang::as_string(rlang::enexpr(id))
  
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  fn_fasta <- file.path(dat_dir, "my_project.fasta")
  stopifnot(file.exists(fn_fasta))
  fasta <- seqinr::read.fasta(fn_fasta, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
  
  df <- read.csv(file.path(dat_dir, "Peptide/Peptide.txt"), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
  
  if (! "pep_isunique" %in% names(df)) {
    df$pep_isunique <- TRUE
    warning("Column `pep_isunique` created and all TRUE values assumed.", call. = FALSE)
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
  
  if (use_unique_pep && "pep_isunique" %in% names(df)) df <- df %>% dplyr::filter(pep_isunique == 1)
  
  df_num <- df %>% 
    dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id))

  df_num <- switch(method_pep_prn, 
                   mean = aggrNums(mean)(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   top.3 = TMT_top_n(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   weighted.mean = tmt_wtmean(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   median = aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.)))
  
  df_mq_rptr_mass_dev <- df %>% 
    dplyr::select(!!rlang::sym(id), grep("^Reporter mass deviation", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep("^Reporter mass deviation", names(.)))	  
  
  mq_median_keys <- c(
    "Score", "Missed cleavages", "PEP", 
    "Charge", "Mass", "PIF", "Fraction of total spectrum", "Mass error [ppm]", 
    "Mass error [Da]", "Base peak fraction", "Precursor Intensity", 
    "Precursor Apex Fraction", "Intensity coverage", "Peak coverage", 
    "Combinatorics"
  )
  
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% mq_median_keys))		
  
  sm_median_keys <- c(
    "deltaForwardReverseScore", "percent_scored_peak_intensity", "totalIntensity", 
    "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
    "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
    "matched_parent_mass", "delta_parent_mass", "delta_parent_mass_ppm")
  
  df_sm_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% sm_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% sm_median_keys))
  
  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
    dplyr::select(-grep("^pep_", names(.)))    
  
  df <- list(df_first, 
             df_mq_rptr_mass_dev, df_mq_med, 
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
    dfa <- df %>% 
      dplyr::select(gene, grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
      dplyr::filter(!is.na(gene)) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::summarise_all(list(~ median(.x, na.rm = TRUE)))
    
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
        write.csv(dup_peps_af, file.path(dat_dir, "Peptide/dbl_dipping_peptides.csv"), row.names = FALSE)
        df <- df %>% dplyr::filter(! (!!rlang::sym(group_psm_by) %in% dup_peps_af[[group_psm_by]]))
      }
    } else {
      write.csv(dup_peps, file.path(dat_dir, "Peptide/dbl_dipping_peptides.csv"), row.names = FALSE)
      df <- df %>% dplyr::filter(! (!!rlang::sym(group_psm_by) %in% dup_peps[[group_psm_by]]))
    }
  }
  
  return(df)
}


