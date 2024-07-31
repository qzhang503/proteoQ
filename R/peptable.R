#' Make new column names
#'
#' According to the Sample_ID(s) in label_scheme.
#'
#' @param TMT_Set Integer; the index of TMT experiment.
#' @param df Data frame; log2FC data.
#' @param label_scheme Experiment summary.
#' @param pattern Regex pattern for capturing intensity and ratio columns.
#' @param group_ids A character vector; for example, c("heavy", "light") that
#'   can be found from column \code{pep_group} in \code{df}; \code{group_ids} is
#'   NULL if column \code{pep_group} is not in \code{df}.
#' @import dplyr purrr forcats
#' @importFrom magrittr %>% %T>% %$% %<>% not
newColnames <- function(TMT_Set = 1L, df, label_scheme, 
                        pattern = "[RI][0-9]{3}[NC]{0,1}", group_ids = NULL) 
{
  label_scheme_sub <- label_scheme %>%
    dplyr::filter(.data$TMT_Set == .env$TMT_Set)
  
  if (is.null(group_ids)) {
    df <- hsubColnames(group_ids, df, pattern, TMT_Set, label_scheme_sub)
  }
  else {
    # no lapply, need side effects
    for (id in group_ids) {
      df <- hsubColnames(id, df, pattern, TMT_Set, label_scheme_sub)
    }
  }

  invisible(df)
}


#' Helper of \link{newColnames}.
#'
#' @param id A character string; for example, "heavy" or "light" that can be
#'   found from column \code{pep_group} in \code{df}. The value of \code{id} is
#'   NULL if column \code{pep_group} is not in \code{df}.
#' @param label_scheme_sub Experiment summary at \code{TMT_Set}.
#' @inheritParams newColnames
hsubColnames <- function (id = NULL, df, pattern = "[RI][0-9]{3}[NC]{0,1}", 
                          TMT_Set = 1L, label_scheme_sub) 
{
  nms <- names(df)
  sids <- as.character(label_scheme_sub$Sample_ID)
  backref <- paste0("(", pattern, ")")
  
  if (is.null(id)) {
    cols <- grep(paste0(pattern, "_", TMT_Set, "$"), nms)
    pat <- paste0(paste(backref, TMT_Set, sep = "_"), "$")
    bares <- gsub(pat, "\\1", nms[cols])
    names(df)[cols] <- paste0(bares, " (", sids, ")")
  }
  else {
    cols <- grep(paste0(paste(pattern, TMT_Set, id, sep = "_"), "$"), nms)
    pat <- paste0(paste(backref, TMT_Set, id, sep = "_"), "$")
    bares <- gsub(pat, "\\1", nms[cols])
    names(df)[cols] <- paste0(bares, " (", sids, " [", id, "]", ")")
  }
  
  invisible(df)
}


#' Check the single presence of MaxQuant peptides[...].txt.
#' 
#' @param dat_dir The working directory.
#' @param label_scheme_full The metadata
#' @inheritParams normPSM
use_mq_peptable <- function (dat_dir, label_scheme_full, 
                             group_psm_by = "pep_seq") 
{
  # only works for "pep_seq"
  filelist <- if (group_psm_by == "pep_seq_mod") {
    list.files(path = file.path(dat_dir), 
               pattern = "^modificationSpecificPeptides.*\\.txt$")
  }
  else if (group_psm_by == "pep_seq") {
    list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  }
  
  if (!length(filelist)) {
    return(FALSE)
  }
  else if (length(filelist) > 1L) {
    stop("Only single MaxQuant `peptides.txt` allowed.", call. = FALSE)
  } 
  else {
    if (!file.exists(file.path(dat_dir, filelist))) {
      stop("No MaxQuant LFQ file `peptides[...].txt` udner", dat_dir)
    } 
    else {
      message("MaxQuant file `", filelist, "` found.")
      return(TRUE)
    }
  }
}


#' Check the single presence of MSFragger peptide table
#' 
#' @param dat_dir The working directory.
#' @param label_scheme_full The metadata
#' @inheritParams normPSM
use_mf_peptable <- function (dat_dir, label_scheme_full, 
                             group_psm_by = "pep_seq_mod") 
{
  filelist <- if (group_psm_by == "pep_seq_mod") {
    list.files(path = file.path(dat_dir), 
               pattern = "^combined_modified_peptide.*\\.tsv$")
  }
  else if (group_psm_by == "pep_seq") {
    list.files(path = file.path(dat_dir), 
               pattern = "^combined_peptide.*\\.tsv$")
  }

  if (!length(filelist)) {
    return(FALSE)
  } 
  else if (length(filelist) > 1L) {
    stop("Only single MSFragger `combined_modified_peptide.tsv` allowed.")
  } 
  else {
    if (!file.exists(file.path(dat_dir, filelist))) {
      stop("No MSFragger `combined_modified_peptide.tsv` udner", dat_dir)
    } 
    else {
      message("MSFragger file `", filelist, "` found.")
      return(TRUE)
    }
  }
}


#' Check the single presence of MSFragger combined_protein[...].tsv
#' 
#' @param dat_dir The working directory.
use_mf_prottable <- function (dat_dir) 
{
  filelist <- list.files(path = file.path(dat_dir), 
                         pattern = "^combined_protein.*\\.tsv$")

  if (!length(filelist)) {
    return(FALSE)
  } 
  else if (length(filelist) > 1L) {
    stop("Only single MSFragger `combined_protein.tsv` allowed.")
  } 
  else {
    if (!file.exists(file.path(dat_dir, filelist))) {
      stop("No MSFragger `combined_protein.tsv` udner", dat_dir)
    } 
    else {
      message("MSFragger file `", filelist, "` found.")
      return(TRUE)
    }
  }
}


#' Check the single presence of MaxQuant combined_protein[...].tsv
#' 
#' @param dat_dir The working directory.
use_mq_prottable <- function (dat_dir) 
{
  filelist <- list.files(path = file.path(dat_dir), 
                         pattern = "^proteinGroups*\\.txt$")
  
  if (!length(filelist)) {
    return(FALSE)
  } 
  else if (length(filelist) > 1L) {
    stop("Only single MaxQuant `proteinGroups.txt` allowed.")
  } 
  else {
    if (!file.exists(file.path(dat_dir, filelist))) {
      stop("No MaxQuant `proteinGroups.txt` udner", dat_dir)
    } 
    else {
      message("MaxQuant file `", filelist, "` found.")
      return(TRUE)
    }
  }
}


#' Adds MSFragger intensity and ratio columns
#'
#' @param dat_dir The working directory.
fillMFprnnums <- function (dat_dir) 
{
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  
  refChannels <- label_scheme |>
    dplyr::filter(Reference) |>
    dplyr::select(Sample_ID) |> 
    unlist()
  
  sids <- label_scheme$Sample_ID
  
  df <- readr::read_tsv(file.path(dat_dir, "combined_protein.tsv"))
  df <- df[, c("Protein ID", paste0(sids, " Intensity"))] |>
    dplyr::rename(prot_acc = `Protein ID`)
  names(df)[2:ncol(df)] <- paste0("I000 (", sids, ")")
  df <- df |> calc_lfq_log2r("I000", refChannels)
  nms <- names(df)
  df <- df |>
    dplyr::mutate_at(.vars = grep("log2_R000\\s", nms), 
                     function (x) replace(x, is.nan(x), NA_real_))

  df_n <- df[, 2:ncol(df), drop = FALSE]
  names(df_n) <- paste0("N_", names(df_n))
  
  df_z <- df_n[, grepl("^N_log2_R000 \\(", names(df_n))]
  names(df_z) <- gsub("^N_log2", "Z_log2", names(df_z))
  
  df <- dplyr::bind_cols(df, df_n, df_z)
  
  nms <- names(df)
  df <- dplyr::bind_cols(
    df |> dplyr::select("prot_acc"),
    df |> dplyr::select(grep("^I000 \\(", nms)),
    df |> dplyr::select(grep("^N_I000 \\(", nms)),
    df |> dplyr::select(grep("^log2_R000 \\(", nms)),
    df |> dplyr::select(grep("^N_log2_R000 \\(", nms)),
    df |> dplyr::select(grep("^Z_log2_R000 \\(", nms)), )
}


#' Adds MSFragger intensity and ratio columns
#'
#' @param dat_dir The working directory.
fillMQprnnums <- function (dat_dir) 
{
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  
  refChannels <- label_scheme |>
    dplyr::filter(Reference) |>
    dplyr::select(Sample_ID) |> 
    unlist()
  
  sids <- label_scheme$Sample_ID
  
  df <- readr::read_tsv(file.path(dat_dir, "proteinGroups.txt"))
  df <- df[, c("Majority protein IDs", paste0("LFQ intensity ", sids))] |>
    dplyr::rename(prot_acc = `Majority protein IDs`) |>
    dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc)) |>
    dplyr::filter(prot_acc != "")

  names(df)[2:ncol(df)] <- paste0("I000 (", sids, ")")
  df <- df |> calc_lfq_log2r("I000", refChannels)
  nms <- names(df)
  df <- df |>
    dplyr::mutate_at(.vars = grep("log2_R000\\s", nms), 
                     function (x) replace(x, is.nan(x), NA_real_))
  
  df_n <- df[, 2:ncol(df), drop = FALSE]
  names(df_n) <- paste0("N_", names(df_n))
  
  df_z <- df_n[, grepl("^N_log2_R000 \\(", names(df_n))]
  names(df_z) <- gsub("^N_log2", "Z_log2", names(df_z))
  
  df <- dplyr::bind_cols(df, df_n, df_z)
  
  nms <- names(df)
  df <- dplyr::bind_cols(
    df |> dplyr::select("prot_acc"),
    df |> dplyr::select(grep("^I000 \\(", nms)),
    df |> dplyr::select(grep("^N_I000 \\(", nms)),
    df |> dplyr::select(grep("^log2_R000 \\(", nms)),
    df |> dplyr::select(grep("^N_log2_R000 \\(", nms)),
    df |> dplyr::select(grep("^Z_log2_R000 \\(", nms)))
}


#' Check the single presence of MaxQuant proteinGroups[...].txt.
#' 
#' @param dat_dir The working directory.
single_mq_prntable <- function (dat_dir) 
{
  filelist <- list.files(path = file.path(dat_dir), 
                         pattern = "^proteinGroups.*\\.txt$")
  
  if (length(filelist) > 1) {
    stop("Only single MaxQuant `proteinGroups.txt` allowed.")
  }
  
  if (!file.exists(file.path(dat_dir, filelist))) {
    stop("No MaxQuant LFQ file `proteinGroups[...].txt` udner", dat_dir)
  }
  
  message("MaxQuant file `", filelist, "` found.")  
  
  invisible(filelist)
}


#' Compare sample IDS between MaxQuant results and expt_smry.xlsx 
#' 
#' @param df A data frame of MaxQuant peptides.txt or proteinGroups.txt.
#' @param label_scheme Experiment summary
check_mq_df <- function (df, label_scheme) 
{
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

    # the same ID occurs in "Identification type", "Experiment", 
    # "Intensity", "LFQ intensity"
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
extract_mq_ints <- function (df) 
{
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
      stop("`type` needs to be either `Intensity` or `LFQ intensity`.", 
           call. = FALSE)
    }
    
    sweep(df[, col_smpls, drop = FALSE], 1,
          rowMeans(df[, col_refs, drop = FALSE], na.rm = TRUE), "/") %>%
      log2(.)  %>% 
      dplyr::mutate_all(~ replace(.x, is.infinite(.), NA)) %>% 
      `colnames<-`(gsub(paste0(type, "(.*)$"), 
                        paste0(prefix, " \\(", "\\1", "\\)"), 
                        names(.))) %>%
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
    dplyr::mutate_at(.vars = grep("log2_R000\\s", names(.)), 
                     ~ replace(.x, is.infinite(.), NA)) 
  
  log2sd <- df %>% 
    dplyr::select(grep("Intensity\\s", names(.))) %>% 
    `names<-`(gsub("^Intensity\\s(.*)", 
                   "sd_log2_R000 \\(\\1\\)", 
                   names(.))) %>% 
    dplyr::mutate_all(~ replace(.x, !is.na(.x), NA))
  
  df <- dplyr::bind_cols(df, log2sd)
  
  df <- df %>% 
    `names<-`(gsub("^Intensity (.*)$", 
                   paste0("I000 \\(", "\\1", "\\)"), 
                   names(.))) %>% 
    `names<-`(gsub("^LFQ intensity (.*)$", 
                   paste0("N_I000 \\(", "\\1", "\\)"), 
                   names(.)))
  
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
#' Fill back intensity values that are not filled in LFQ msms.txt.
#' 
#' @param label_scheme Experiment summary
#' @inheritParams n_TMT_sets
#' @inheritParams mergePep
pep_mq_lfq <- function(label_scheme, omit_single_lfq) 
{
  dat_dir <- get_gl_dat_dir()
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)

  filelist <- 
    list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  
  if (!length(filelist)) {
    stop(paste("No MaxQuant LFQ file of `peptides[...].txt` under", 
               file.path(dat_dir), ".\n",
               "Make sure that the name of file starts with `peptides`."))
  }
  
  df <- local({
    df <- read.csv(file.path(dat_dir, filelist), check.names = FALSE, 
                   header = TRUE, sep = "\t", comment.char = "#")
    
    # Leading razor protein
    
    df <- df %>% 
      dplyr::mutate(Proteins = gsub("\\.[0-9]*", "", Proteins), 
                    `Leading razor protein` = 
                      gsub("\\.[0-9]*", "", `Leading razor protein`))

    df <- local({
      if (purrr::is_empty(grep("^LFQ intensity |^Intensity ", names(df)))) {
        nas <- data.frame(rep(NA, nrow(df)))
        sample_ids <- as.character(label_scheme$Sample_ID)
        
        df_int <- purrr::map(sample_ids, ~ {
          nas %>% `colnames<-`(paste("Intensity", .x))
        }) %>% 
          dplyr::bind_cols()
        
        df_int2 <- purrr::map(sample_ids, ~ {
          nas %>% `colnames<-`(paste("LFQ intensity", .x))
        }) %>% 
          dplyr::bind_cols()
        
        df <- dplyr::bind_cols(df, df_int, df_int2)
      } else if (purrr::is_empty(grep("^LFQ intensity ", names(df)))) {
        warning("Columns `LFQ intensity` not found in `", filelist, "`; ", 
                "columns `Intensity` used instead.")
        
        df_int <- df %>% 
          dplyr::select(grep("^Intensity ", names(.))) %>% 
          `names<-`(gsub("^Intensity ", "LFQ intensity ", names(.)))
        
        df <- dplyr::bind_cols(df, df_int)
      } else if (purrr::is_empty(grep("^Intensity ", names(df)))) {
        warning("Columns `Intensity` not found in `", filelist, "`; ", 
                "columns `LFQ intensity` used instead.")
        
        df_int <- df %>% 
          dplyr::select(grep("^LFQ intensity ", names(.))) %>% 
          `names<-`(gsub("^LFQ intensity ", "Intensity ", names(.)))
        
        df <- dplyr::bind_cols(df, df_int)
      }
      
      invisible(df)
    })
    
    # MaxQuant peptide.txt may contains extra_sample_ids not in label_scheme
    # e.g. one may delete a Sample_ID from label_scheme, 
    #   the corresponding Sample_ID will be removed from compiled PSM tables.
    # However, the same needs to be done 
    #  when backfilling intensity data from peptide.txt.
    # Otherwise would cause columns of `pep_razor_int (extra_sample_ids)` etc. 
    #   in peptide table.
    # This will cause an error with Pep2Prn(method_pep_prn = lfq_...), 
    # which involves columns of `pep_razor_int (extra_sample_ids)` etc. 
    # for the identification of top_n.
    df <- local({
      extra_sample_ids <- names(df) %>% 
        .[grep("^Intensity\\s.*", .)] %>% 
        gsub("^Intensity\\s(.*)$", "\\1", .) %>% 
        .[! . %in% label_scheme$Sample_ID]
      
      if (!purrr::is_empty(extra_sample_ids)) {
        warning("\nSample IDs in `peptide.txt` not found", 
                " in `label_scheme.xlsx` and removed: \n", 
                purrr::reduce(extra_sample_ids, paste, sep = ", "), 
                "\n\n=======================================================================", 
                "\nWith the temporary data backfilling of `msms.txt` <- `peptide.txt`, ", 
                "\nit is currently not possible to tell the above mismatches is by either \n",
                "(a) intended sample removals in `label_scheme.xlsx` or \n",
                "(b) inadvertence in naming samples differently between `label_scheme.xlsx` ", 
                "and MaXquant searches.\n\n",
                
                "Fow now, identical sample IDs between `label_scheme.xlsx` ", 
                "and `peptide.txt` are required ", 
                "(only with the backfilling procedures).\n", 
                "With the temporary requirement of identical sample IDs being fulfilled, ", 
                "users may then perform sample exclusions via `label_scheme.xlsx`.",
                "\n=======================================================================\n", 
                call. = FALSE)
        
        purrr::walk(extra_sample_ids, ~ {
          smpl <- paste0(" ", .x, "$")
          df <<- df %>% dplyr::select(-grep(smpl, names(.)))
          
          message("Mismatched sample ID removed from `peptide.txt`: ", .x)
        })
        
        if (purrr::is_empty(grep("^LFQ intensity |^Intensity ", names(df)))) {
          stop("No samples matched to the IDs in `label_scheme.xlsx`.", 
               call. = FALSE)
        }
      }
      
      invisible(df)
    })
    
    df <- df %>% 
      dplyr::filter(not_allzero_rows(.[grep("^LFQ intensity ", names(.))])) 
  })

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
      df <- df %>% 
        dplyr::mutate(prot_acc = gsub(";.*$", "", Proteins))
      
      tempdata <- df %>% 
        dplyr::select(prot_acc) %>% 
        dplyr::filter(!duplicated(prot_acc)) %>% 
        annotPrn(fasta, entrez) %>% 
        dplyr::select(prot_acc, gene, prot_desc) %>% 
        dplyr::rename(`Gene Names` = "gene", 
                      "Protein Names" = "prot_desc")
      
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
    # df <- df %>% 
    #   dplyr::mutate(pep_seq = paste(`Amino acid before`, pep_seq, `Amino acid after`, 
    #                                 sep = ".")) 
  }
  
  df <- local({
    df_vals <- df %>% 
      dplyr::select(group_psm_by, 
                    grep("^Intensity\\s|^LFQ\\sintensity\\s", names(.))) %>% 
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
#' 
#' @param label_scheme Experiment summary.
pep_mq_lfq2 <- function(label_scheme) 
{
  dat_dir <- get_gl_dat_dir()
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  corrected_int <- match_call_arg(normPSM, corrected_int)
  
  filelist <- list.files(path = file.path(dat_dir), 
                         pattern = "^peptides.*\\.txt$")
  
  if (purrr::is_empty(filelist)) {
    stop(paste("No MaxQuant LFQ file of `peptides[...].txt` under", 
               file.path(dat_dir), ".\n",
               "Make sure that the name of file starts with `peptides`."), 
         call. = FALSE)
  }
  
  # identical codes to those in pep_mq_lfq
  # since temporary exception handling so no makes of function
  df <- local({
    df <- read.csv(file.path(dat_dir, filelist), check.names = FALSE, 
                   header = TRUE, sep = "\t", comment.char = "#")
    
    df <- local({
      if (purrr::is_empty(grep("^LFQ intensity |^Intensity ", names(df)))) {
        nas <- data.frame(rep(NA, nrow(df)))
        sample_ids <- as.character(label_scheme$Sample_ID)
        
        df_int <- purrr::map(sample_ids, ~ {
          nas %>% `colnames<-`(paste("Intensity", .x))
        }) %>% 
          dplyr::bind_cols()
        
        df_int2 <- purrr::map(sample_ids, ~ {
          nas %>% `colnames<-`(paste("LFQ intensity", .x))
        }) %>% 
          dplyr::bind_cols()
        
        df <- dplyr::bind_cols(df, df_int, df_int2)
      } else if (purrr::is_empty(grep("^LFQ intensity ", names(df)))) {
        warning("Columns `LFQ intensity` not found in `", filelist, "`; ", 
                "columns `Intensity` used instead.",
                call. = FALSE)
        
        df_int <- df %>% 
          dplyr::select(grep("^Intensity ", names(.))) %>% 
          `names<-`(gsub("^Intensity ", "LFQ intensity ", names(.)))
        
        df <- dplyr::bind_cols(df, df_int)
      } else if (purrr::is_empty(grep("^Intensity ", names(df)))) {
        warning("Columns `Intensity` not found in `", filelist, "`; ", 
                "columns `LFQ intensity` used instead.",
                call. = FALSE)
        
        df_int <- df %>% 
          dplyr::select(grep("^LFQ intensity ", names(.))) %>% 
          `names<-`(gsub("^LFQ intensity ", "Intensity ", names(.)))
        
        df <- dplyr::bind_cols(df, df_int)
      }

      invisible(df)
    })
    
    # MaxQuant peptide.txt may contains extra_sample_ids not in label_scheme
    # e.g. one may delete a Sample_ID from label_scheme, 
    #   the corresponding Sample_ID will be removed from compiled PSM tables.
    # However, the same needs to be done 
    #  when backfilling intensity data from peptide.txt.
    # Otherwise would cause columns of `pep_razor_int (extra_sample_ids)` etc. 
    #   in peptide table.
    # This will cause an error with Pep2Prn(method_pep_prn = lfq_...), 
    # which involves columns of `pep_razor_int (extra_sample_ids)` etc. 
    # for the identification of top_n.
    df <- local({
      extra_sample_ids <- names(df) %>% 
        .[grep("^Intensity\\s.*", .)] %>% 
        gsub("^Intensity\\s(.*)$", "\\1", .) %>% 
        .[! . %in% label_scheme$Sample_ID]
      
      df <- df %>% 
        dplyr::select(-grep(paste0(" ", extra_sample_ids, "$"), names(.)))
    })
    
    df <- df %>% 
      dplyr::filter(not_allzero_rows(.[grep("^LFQ intensity ", names(.))])) 
  })
  
  
  stopifnot("Proteins" %in% names(df), 
            "Sequence" %in% names(df), 
            "Unique (Groups)" %in% names(df), 
            "Unique (Proteins)" %in% names(df), 
            "Amino acid before" %in% names(df), 
            "Amino acid after" %in% names(df))
  
  df <- df %>% dplyr::mutate(prot_acc = gsub(";.*$", "", Proteins))
  
  prefix <- if (corrected_int) "LFQ intensity" else "^Intensity"
  
  df_slim <- local({
    sample_ids <- paste0("LFQ intensity ", label_scheme$Sample_ID)
    cols <- c("Sequence", "Unique (Groups)", "Unique (Proteins)")
    
    df %>% 
      dplyr::select(which(names(.) %in% cols), 
                    grep(prefix, names(.))) %>% 
      dplyr::select(c(cols, sample_ids))
  })
  
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


#' Handling of MSFragger peptide table
#'
#' Fill back intensity values.
#' 
#' @param label_scheme Experiment summary
#' @inheritParams n_TMT_sets
#' @inheritParams normPSM
pep_mf_lfq <- function(label_scheme, omit_single_lfq, group_psm_by, 
                       use_lowercase_aa = FALSE) 
{
  dat_dir <- get_gl_dat_dir()
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  filelist <- if (group_psm_by == "pep_seq_mod") {
    list.files(path = file.path(dat_dir), pattern = "^combined_modified_peptide.*\\.tsv$")
  }
  else if (group_psm_by == "pep_seq") {
    list.files(path = file.path(dat_dir), pattern = "^combined_peptide.*\\.tsv$")
  }
  else {
    stop("Unhandled condition.")
  }

  n_files <- length(filelist)
  if (!n_files) {
    stop(paste("No MSFragger LFQ file of `combined_peptide[...].tsv` under", 
               file.path(dat_dir), ".\n",
               "Make sure that the name of file starts with `combined_peptide`."))
  }
  
  if (n_files > 1L) {
    stop(paste("Only allowed one MSFragger LFQ file of `combined_peptide[...].tsv` under", 
               file.path(dat_dir), ".\n"))
  }
  
  df <- readr::read_tsv(file.path(dat_dir, filelist))
  pat  <- " Intensity"
  pat2 <- " MaxLFQ Intensity"
  col_nms <- names(df)
  
  col_ints  <- grep(pat, col_nms)
  col_ints2 <- grep(pat2, col_nms)
  col_ints  <- col_ints[!col_ints %in% col_ints2]
  sids <- gsub(pat, "", col_nms[col_ints])
  names(df)[col_ints]  <- paste0("I000 (", sids, ")")
  
  # can be empty for unknown reason
  # no_maxlfq <- length(col_ints2) == 0L
  n_ints2 <- length(col_ints2)
  if (n_ints2) {
    names(df)[col_ints2] <- paste0("N_I000 (", sids, ")")
  }

  df <- df |>
    dplyr::select(which(col_nms %in% c("Peptide Sequence", "Modified Sequence", 
                                       "Prev AA", "Next AA", "Start", "End")), 
                  col_ints, col_ints2) |>
    dplyr::rename(pep_seq = `Peptide Sequence`, 
                  pep_seq_mod = `Modified Sequence`, 
                  pep_res_before = `Prev AA`, 
                  pep_res_after = `Next AA`, 
                  pep_start = Start, 
                  pep_end = End)

  if (group_psm_by == "pep_seq_mod") {
    df <- df |> dplyr::select(-"pep_seq")
    
    df <- df |>
      reloc_col_before(group_psm_by, "pep_res_after") |>
      reloc_col_before("pep_res_before", group_psm_by) |>
      add_msfragger_pepseqmod(use_lowercase_aa)
  }
  else if (group_psm_by == "pep_seq") {
    df <- df |> dplyr::select(-"pep_seq_mod")
  }

  if (omit_single_lfq) {
    df <- na_single_lfq(df, "^N_I000 ")
  }

  df <- df |>
    dplyr::select(-c("pep_res_before", "pep_res_after", "pep_start", "pep_end"))
  
  df_empty <- matrix(ncol = length(sids), nrow = nrow(df)) |>
    data.frame() |>
    tibble::tibble()
  
  # special handling
  if (!n_ints2) {
    df_n_ints <- df[, grepl("^I000 ", names(df))]
    names(df_n_ints) <- paste0("N_I000 (", sids, ")")
    df <- dplyr::bind_cols(df, df_n_ints)
  }
  
  df_n_log2r <- df_log2r <- df_sd <- df_empty
  names(df_sd) <- paste0("sd_log2_R000 (", sids, ")")
  names(df_log2r) <- paste0("log2_R000 (", sids, ")")
  names(df_n_log2r) <- paste0("N_log2_R000 (", sids, ")")
  
  df <- dplyr::bind_cols(df, df_sd, df_log2r, df_n_log2r)
  
  ## log2 ratios
  refChannels <- label_scheme |>
    dplyr::filter(Reference) |>
    dplyr::select(Sample_ID) |>
    unlist()
  
  df <- df |>
    calc_lfq_log2r("I000", refChannels) |>
    calc_lfq_log2r("N_I000", refChannels) |>
    # NaN for log2 if zero-intensity under reference channels: log(0/0)
    dplyr::mutate_if(is.numeric, function (x) replace(x, is.nan(x), NA_real_))
  
  df <- df |>
    dplyr::mutate_at(grep("log2_R000\\s", names(df)), 
                     function (x) replace(x, is.infinite(x), NA_real_)) 
  
  df
}


#' Helper: calculates LFQ log2FC.
#'
#' For rows with single intensity values: original intensities kept but
#' \code{log2_R000} coerced from 0 to NA. This will keep single-value rows away
#' from data alignment that are based on \code{log2_R000}. Otherwise, the
#' alignment may be trapped to the spike at \code{log2_R000 = 0}.
#'
#' @param df A data frame.
#' @param type The type of intensity data.
#' @param refChannels The reference channels.
calc_lfq_log2r <- function (df, type, refChannels) 
{
  if (!type %in% c("I000", "N_I000")) {
    stop("The value of `type` need to be `I000` or `N_I000`.")
  }

  new_type <- paste0("^", type, " ")
  prefix <- gsub("I000", "log2_R000", new_type)
  no_hat <- gsub("\\^", "", prefix)
  
  df <- df[, !grepl(prefix, names(df))]
  col_smpls <- grep(new_type, names(df))
  
  if (!length(col_smpls)) {
    stop("No sample columns start with ", type, ".")
  }

  col_refs <- if (length(refChannels)) 
    which(names(df) %in% paste0(type, " (", refChannels, ")"))
  else 
    col_smpls
  
  dfx <- df[, col_smpls, drop = FALSE]
  rows_sv <- (rowSums(!is.na(dfx)) == 1L)
  
  ans <- sweep(dfx, 1,
        rowMeans(df[, col_refs, drop = FALSE], na.rm = TRUE), "/") %>%
    log2(.)  %>% 
    dplyr::mutate_all(~ replace(.x, is.infinite(.), NA_real_)) %>% 
    `colnames<-`(gsub(paste0(new_type, "(.*)$"), 
                      paste0(no_hat, "\\1"), 
                      names(.)))
  
  ans[rows_sv, ] <- NA_real_
  
  cbind(df, ans)
}


#' Trivializes data rows with single LFQ intensity
#' 
#' @param df A data frame.
#' @param pattern The pattern of intensity fields.
na_single_lfq <- function (df, pattern = "^I000 ") 
{
  cols <- grep(pattern, names(df))
  df_lfq <- df[, cols, drop = FALSE]

  # slightly more flexible than (rowSums(!is.na(df_lfq)) > 1L)
  # in case that upstream NA's are replaced 0's
  # (e.g. MaxQuant's timsTOF PSMs are all NA values)
  not_single_zero <- (rowSums(df_lfq > 0, na.rm = TRUE) > 1L) 
  not_single_zero[!not_single_zero] <- NA # NA logical
  
  df_lfq[] <- lapply(df_lfq, `*`, not_single_zero)
  df[, cols] <- df_lfq
  
  invisible(df)
}


#' Calculates log2FC of peptides based on LFQ intensity.
#'
#' With LFQ, values of \code{log2_R000 (...)} and \code{N_log2_R000 (...)} in
#' \code{df_num} are not yet filled after \link{spreadPepNums}. This utility
#' calculates the ratio using the values of LFQ intensity.
#'
#' For SILAC, intensity of heave, light etc. have been spread into multiple
#' columns and thus the log2Ratios can be calculated like regular LFQ.
#'
#' @param df A data frame.
#' @param min_int Minimal Y intensity to replace 0 values.
calcLFQPepNums <- function (df, min_int = 1E3, imp_vals = FALSE) 
{
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  refChannels  <- label_scheme$Sample_ID[label_scheme[["Reference"]]]
  
  nms <- names(df)
  cols_both <- nms[grepl("^I000 \\(", nms)]

  if (length(refChannels)) {
    cols_ref  <- paste0("I000 (", refChannels, ")")
    cols_smpl <- cols_both[!cols_both %in% cols_ref]
  }
  else {
    cols_smpl <- cols_ref <- cols_both
  }

  # (0) no possible to have all NA rows for both references and samples
  if (FALSE) {
    row_sums <- df[, c(cols_ref, cols_smpl), drop = FALSE] |>
      is.na() |>
      rowSums()
    oks <- row_sums < length(c(cols_ref, cols_smpl))
    df <- df[oks, ]
    rm(list = c("row_sums", "oks"))
  }
  
  # (1) handle missing reference
  df <- local({
    ref_rowsums <- df[, cols_ref, drop = FALSE] |>
      is.na() |>
      rowSums()
    empties <- ref_rowsums == length(cols_ref)
    df0 <- df[empties, ] # works at no empties
    df1 <- df[!empties, ] # some references are NA
    ys  <- rowMeans(df0[, cols_smpl, drop = FALSE], na.rm = TRUE)

    df0[, cols_ref] <- lapply(cols_ref, function (col) {
      df0[, col] <- ys
    })
    
    df0[, paste0("N_", cols_ref)] <- df0[, cols_ref, drop = FALSE]
    df <- dplyr::bind_rows(df0, df1)
  })

  # 2. With references
  if (imp_vals && length(cols_smpl) > 1L) {
    df <- local({
      smpl_rowsums <- rowSums(is.na(df[, cols_smpl, drop = FALSE]))
      empties <- smpl_rowsums == length(cols_smpl)
      df0 <- df[empties, ]
      df1 <- df[!empties, ]
      
      nms  <- names(df)
      cols <- which(nms %in% cols_smpl)
      class(df1) <- "data.frame" # remove `tbl_df` and `tbl`
      # df1[, cols] <- impPepNA_byrow(df1[, cols, drop = FALSE])
      df1[, cols] <- impPepNA(df1[, cols, drop = FALSE])
      
      df1[, paste0("N_", cols_smpl)] <- df1[, cols_smpl, drop = FALSE]
      df <- dplyr::bind_rows(df0, df1)
    })
  }

  if (FALSE) {
    df <- df |>
      na_single_lfq("^I000 ") |>
      na_single_lfq("^N_I000 ")
  }

  df <- df |>
    # dplyr::mutate_at(grep("I[0-9]{3}[NC]{0,1}", names(df)), 
    #                  function (x) { x[x < min_int] <- min_int; x }) |>
    calc_lfq_log2r("I000", refChannels) |>
    calc_lfq_log2r("N_I000", refChannels) |>
    # NaN for log2 if zero-intensity under reference channels: log(0/0)
    dplyr::mutate_if(is.numeric, function (x) replace(x, is.nan(x), NA_real_))

  df <- df |>
    dplyr::mutate_at(grep("log2_R000\\s", names(df)), 
                     function (x) replace(x, is.infinite(x), NA_real_)) 

  if (!length(grep("log2_R000", names(df)))) {
    stop("No `log2_R000...` columns available.\n",
         "Probably inconsistent sample IDs between metadata and ", 
         "MaxQuant `peptides.txt`.")
  }
  
  nms <- names(df)
  df <- dplyr::bind_cols(
    df |> dplyr::select(-grep("[IR]{1}000 \\(", nms)),
    df |> dplyr::select(grep("^I000 \\(", nms)),
    df |> dplyr::select(grep("^N_I000 \\(", nms)),
    df |> dplyr::select(grep("^sd_log2_R000 \\(", nms)),
    df |> dplyr::select(grep("^log2_R000 \\(", nms)),
    df |> dplyr::select(grep("^N_log2_R000 \\(", nms)))
}


#' Calculates \code{pep_tot_int}, \code{pep_razor_int} and
#' \code{pep_unique_int} in LFQ
#' 
#' @param df A data frame.
#' @param filelist A list of individual peptide tables.
#' @inheritParams normPSM 
calcLFQPepInts <- function (df, filelist, group_psm_by) 
{
  dat_dir <- get_gl_dat_dir()
  label_scheme <- load_ls_group(dat_dir, label_scheme, prefer_group = FALSE)
  
  cols_grp  <- c(group_psm_by, "TMT_Set")
  cols_grp2 <- c("TMT_Set")
  nms <- names(df)
  
  lapply(cols_grp, 
         function (x) if (!x %in% nms) stop("Column `", x, "` not found."))
  
  group_ids <- if ("pep_group" %in% nms) unique(df$pep_group) else NULL
  is_mulgrps <- length(group_ids) >= 2L
  
  if (is_mulgrps) {
    cols_grp <- c(cols_grp, "pep_group")
    cols_grp2 <- c(cols_grp2, "pep_group")
  }

  df_num <- 
    df[, c(cols_grp, "pep_tot_int", "pep_unique_int", "pep_razor_int")] |>
    dplyr::group_by_at(cols_grp) |>
    dplyr::arrange(TMT_Set)
  nms <- names(df_num)
  df_num <- df_num |>
    tidyr::gather(grep("^pep_.*_int$", nms), key = ID, value = value) |>
    tidyr::unite(ID, ID, cols_grp2, sep = "_")
  
  ## pep_tot_int_1_light, pep_unique_int_1_light, pep_razor_int_1_light
  type_levels <- c("pep_tot_int", "pep_unique_int", "pep_razor_int")
  tmt_levels <- NULL
  pat <- "pep_.*_int"
  fct_cols <- c("type", "set", "group")
  
  lapply(c("type_levels"), function (x) {
    if (is.factor(df_num[[x]])) stop("`", x, "` cannot be factor.")
  })
  
  df_num <- df_num %>% 
    dplyr::mutate(type = gsub(paste0("^(", pat, ")", "_.*"), "\\1", ID), 
                  set = gsub(paste0("^", pat, "_([0-9]+).*"), "\\1", ID), ) |>
    dplyr::mutate(type = factor(type, levels = type_levels), 
                  set = as.integer(set), )
  
  lapply(c("type", "set"), function (x) {
    if (any(is.na(df_num[[x]]))) stop("Unexpected NA under column `", x, "`.")
  })

  df_num <- df_num |>
    dplyr::mutate(group = gsub(paste0("^", pat, "_[0-9]+_(.*)$"), "\\1", ID), 
                  group = factor(group, levels = group_ids)) |>
    dplyr::arrange_at(fct_cols)
  df_num <- df_num[, !names(df_num) %in% fct_cols]
  id_levels <- unique(df_num$ID)
  
  df_num <- df_num %>%
    dplyr::group_by(pep_seq_mod, ID) |>
    dplyr::summarise_at("value", sum, na.rm = TRUE) |>
    dplyr::mutate(ID = factor(ID, levels = id_levels)) |>
    tidyr::pivot_wider(names_from = ID, values_from = value)
  
  set_indexes <- gsub("^.*TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort()
  
  for (set_idx in set_indexes) {
    df_num <- newColnames(TMT_Set = set_idx, 
                          df = df_num, 
                          label_scheme = label_scheme, 
                          pattern = "pep_.*_int", 
                          group_ids = group_ids)
  }
  
  if (is_mulgrps) {
    label_scheme_group <- 
      load_ls_group(dat_dir, label_scheme, prefer_group = TRUE)
    
    df_num <- pad_grp_samples(
      df = df_num, 
      sids = as.character(label_scheme_group$Sample_ID), 
      tmt_levels = tmt_levels, 
      type_levels = type_levels)
  }
  
  df_num %>% 
    dplyr::select(!!rlang::sym(group_psm_by), 
                  grep("^pep_.*_int ", names(.))) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
}


#' Spreads peptide numbers.
#'
#' Spreads fields of numeric values: sd_log2_R, log2_R, log2_R, I, N_I by TMT
#' sets.
#'
#' Also works for LFQ as each sample corresponds to a TMT set.
#'
#' For single SILAC sample, the values of log2Ratios spreads into
#' \emph{MULTIPLE} columns of heavy, light etc. Despite, log2Ratios remains NA,
#' just like regular single-sample LFQ. The log2Ratios will be later calculated
#' with \link{calcLFQPepNums} that are based on intensity values.
#'
#' @param df A data frame of peptide table.
#' @param basenames Names of peptide table files.
#' @param ok_mbr Logical Capable of MBR or not. 
#' @inheritParams splitPSM_ma
#' @inheritParams annotPSM
spreadPepNums <- function (df, basenames, group_psm_by, ok_mbr = FALSE) 
{
  dat_dir <- get_gl_dat_dir()
  
  label_scheme <- 
    load_ls_group(dat_dir, label_scheme, prefer_group = FALSE)
  label_scheme_full <- 
    load_ls_group(dat_dir, label_scheme_full, prefer_group = FALSE)
  tmt_plexes <- TMT_plex(label_scheme)
  tmt_levels <- TMT_levels(tmt_plexes)
  tmt_levels <- if (tmt_plexes) gsub("^TMT-", "", tmt_levels) else "000"
  
  cols_grp  <- c(group_psm_by, "TMT_Set")
  cols_grp2 <- c("TMT_Set")
  nms <- names(df)
  
  lapply(cols_grp, function (x) 
    if (! x %in% nms) stop("Column `", x, "` not found."))

  group_ids  <- if ("pep_group" %in% nms) unique(df$pep_group) else NULL
  is_mulgrps <- length(group_ids) >= 2L

  if (is_mulgrps) {
    cols_grp  <- c(cols_grp, "pep_group")
    cols_grp2 <- c(cols_grp2, "pep_group")
    
    label_scheme_group <- rep_ls_groups(group_ids)
    label_scheme_full_group <- rep_ls_groups(group_ids)
    
    save(label_scheme_group, 
         file = file.path(dat_dir, "label_scheme_group.rda"))
    save(label_scheme_full_group, 
         file = file.path(dat_dir, "label_scheme_full_group.rda"))
    
    write_excel_wb(label_scheme_full_group, "Setup", dat_dir, 
                   "label_scheme_full_group.xlsx")
    write_excel_wb(label_scheme_group, "Setup", dat_dir, 
                   "label_scheme_group.xlsx")
  }
  else {
    # save(label_scheme, file = file.path(dat_dir, "label_scheme_group.rda"))
    # save(label_scheme_full, file = file.path(dat_dir, "label_scheme_group.rda"))
  }

  ## Numeric fields
  pat_sd_log2_R <- "^sd_log2_R[0-9]{3}[NC]{0,1}"
  pat_log2_R <- "^log2_R[0-9]{3}[NC]{0,1}"
  pat_N_log2_R <- "^N_log2_R[0-9]{3}[NC]{0,1}"
  # pat_Z_log2_R <- "^Z_log2_R[0-9]{3}[NC]{0,1}"
  pat_I <- "^I[0-9]{3}[NC]{0,1}"
  pat_N_I <- "^N_I[0-9]{3}[NC]{0,1}"
  
  # aggregation different LCMS injections under the same TMT_Set
  nms <- names(df)
  
  if (tmt_plexes) {
    df_num <- df |>
      dplyr::select(cols_grp, 
                    grep(pat_sd_log2_R, nms), 
                    grep(pat_log2_R, nms), 
                    grep(pat_N_log2_R, nms), 
                    # grep(pat_Z_log2_R, nms), 
                    grep(pat_I, nms), 
                    grep(pat_N_I,nms)) |>
      dplyr::group_by_at(cols_grp) |>
      aggrNumLCMS(group_psm_by, label_scheme_full) |>
      dplyr::mutate(TMT_Set = as.integer(TMT_Set))
  }
  else {
    if (ok_mbr) {
      df_mbr <- lapply(basenames, function (file) {
        fi <- paste0("MBRpeps_", file)
        df <- readr::read_tsv(file.path(dat_dir, "Peptide/cache", fi))
        idx_tmtset <- as.integer(gsub("^.*TMTset(\\d+).*", "\\1", file))
        df$TMT_Set <- idx_tmtset
        df
      }) |>
        dplyr::bind_rows() |>
        dplyr::rename(I000 = pep_tot_int)
    }
    else {
      df_mbr <- df
    }
    
    df_num <- df_mbr |>
      dplyr::select(c(cols_grp, "I000")) |>
      dplyr::group_by_at(cols_grp) |>
      aggrNumLCMS(group_psm_by, label_scheme_full) |>
      dplyr::ungroup() |>
      dplyr::mutate(TMT_Set = as.integer(TMT_Set)) |>
      dplyr::mutate(N_I000 = I000, sd_log2_R000 = NA_real_, 
                    log2_R000 = NA_real_, N_log2_R000 = NA_real_)
  }
  
  # if (anyDuplicated(df_num$pep_seq_mod)) {
  #   stop("Developer: duplicated pep_seq_mod detected.")
  # }

  ## Format: ..._log2_R126_1_[base], ..._I126_1_[base]
  nms <- names(df_num)
  df_num <- df_num |>
    tidyr::gather(grep("R[0-9]{3}[NC]{0,1}|I[0-9]{3}[NC]{0,1}", nms), 
                  key = ID, value = value) |>
    dplyr::arrange_at(cols_grp2) |>
    tidyr::unite(ID, c(ID, cols_grp2))
  
  type_levels <- c("I", "N_I", "sd_log2_R", "log2_R", "N_log2_R")

  ## define the levels of TMT channels;
  #  otherwise, the order of channels will flip between N(itrogen) and C(arbon)
  lapply(c("type_levels", "tmt_levels"), function (x) {
    if (is.factor(df_num[[x]])) stop("`", x, "` cannot be factor.")
  })
  
  df_num <- df_num |>
    dplyr::mutate(type = gsub("^(.*[RI]{1})[0-9]{3}[NC]{0,1}_.*", "\\1", ID), 
                  set = gsub("^.*[RI]{1}[0-9]{3}[NC]{0,1}_([0-9]+).*", "\\1", ID), 
                  channel = gsub("^.*[RI]{1}([0-9]{3}[NC]{0,1})_.*", "\\1", ID)) |>
    dplyr::mutate(type = factor(type, levels = type_levels), 
                  set = as.integer(set), 
                  channel = factor(channel, levels = tmt_levels)) 
  
  lapply(c("type", "channel", "set"), function (x) {
    if (any(is.na(df_num[[x]]))) stop("Unexpected NA under column `", x, "`.")
  })
  
  fct_cols <- c("type", "set", "channel", "group")
  
  df_num  <- df_num |>
    dplyr::mutate(
      group = gsub("^.*[RI]{1}[0-9]{3}[NC]{0,1}_[0-9]+_(.*)$", "\\1", ID), 
      group = factor(group, levels = group_ids)) |>
    dplyr::arrange_at(fct_cols) 
  df_num <- df_num[, !names(df_num) %in% fct_cols]

  id_levels <- unique(df_num$ID)
  
  df_num <- df_num |>
    dplyr::mutate(ID = factor(ID, levels = id_levels)) |>
    tidyr::pivot_wider(names_from = ID, values_from = value)

  set_indexes <- gsub("^.*TMTset(\\d+).*", "\\1", basenames) |>
    unique() |>
    as.integer() |>
    sort()
  
  if (is.factor(group_ids)) {
    stop("`group_ids` cannot be factor.")
  }

  for (set_idx in set_indexes) {
    df_num <- newColnames(TMT_Set = set_idx, 
                          df = df_num, 
                          label_scheme = label_scheme, 
                          pattern = "[RI][0-9]{3}[NC]{0,1}", 
                          group_ids = group_ids)
  }
  
  nms <- names(df_num)
  df_num <- df_num %>% 
    dplyr::select(!!rlang::sym(group_psm_by), 
                  grep("[RI][0-9]{3}[NC]{0,1}", nms)) |>
    dplyr::ungroup() |>
    dplyr::arrange(!!rlang::sym(group_psm_by))
  
  if (is_mulgrps) {
    df_num <- pad_grp_samples(
      df = df_num, 
      sids = as.character(label_scheme_group$Sample_ID), 
      tmt_levels = tmt_levels, 
      type_levels = type_levels)
  }

  invisible(df_num)
}


#' Aggregates data from the same TMT_Set at different LCMS_Injection.
#' 
#' @param df A data frame.
#' @param label_scheme_full Metadata.
#' @param cols Columns for data aggregation.
aggrNumLCMS <- function (df, group_psm_by = "pep_seq_mod", label_scheme_full, 
                         cols = NULL)
{
  tbl_lcms <- n_LCMS(label_scheme_full)
  
  if (!any(tbl_lcms$n_LCMS > 1L)) {
    return(df)
  }
  
  tb_n <- dplyr::filter(tbl_lcms, n_LCMS > 1L)
  rows <- df$TMT_Set %in% tb_n$TMT_Set
  df_n <- df[rows, ]
  df_1 <- df[!rows, ]
  
  nms_n <- names(df_n)
  
  if (is.null(cols)) {
    cols <- grep("log2_R[0-9]{3}[NC]{0,1}|I[0-9]{3}[NC]{0,1}", nms_n)
    len  <- length(cols)
  }
  else {
    cols <- which(cols %in% nms_n)
    len  <- length(cols)
    
    if (!len) {
      stop("Not all column found for data summary.")
    }
  }
  
  col1 <- which(nms_n == group_psm_by)
  col2 <- which(nms_n == "TMT_Set")
  col3 <- which(nms_n == "pep_group")
  cols123 <- c(col1, col2, col3)
  
  df_n <- df_n |> dplyr::group_by_at(cols123)
  
  if (len == 1L) {
    df <- df_n[, c(cols123, cols)] |> 
      dplyr::summarise_all(function (x) median(x, na.rm = TRUE)) |>
      dplyr::ungroup()
    return(df)
  }
  
  # slower; later to summary only over the original intensity columns
  if (FALSE) {
    df_n <- df_n |> 
      dplyr::group_by_at(cols123) |>
      dplyr::group_split() # 40094
    
    nrows <- unlist(lapply(df_n, nrow))
    df_n0 <- df_n[nrows == 1L] # 13240
    df_n0 <- dplyr::bind_rows(df_n0)
    df_n1 <- df_n[nrows > 1L] # 26854
    
    df_n1_keys <- lapply(df_n1, function (x) x[1, -cols]) |>
      dplyr::bind_rows() # 26854
    
    data <- lapply(df_n1, function (x) 
      dplyr::summarise_at(x, .vars = cols, median, na.rm = TRUE))
    data <- dplyr::bind_rows(data)
    df_n1 <- dplyr::bind_cols(df_n1_keys, data)

    if (FALSE) {
      n_cores <- max(min(parallel::detectCores() - 1L), 4L)
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      
      data <- parallel::clusterApply(
        cl, parallel::clusterSplit(cl, df_n1), 
        function (dfs) {
          lapply(dfs, function (df) 
            dplyr::summarise_at(df, .vars = cols, median, na.rm = TRUE)) |>
            dplyr::bind_rows()
        }
      )
      parallel::stopCluster(cl)
      
      data <- dplyr::bind_rows(data)
      df_n1 <- dplyr::bind_cols(df_n1_keys, data)
      df_n <- dplyr::bind_rows(df_n0, df_n1)
      df <- dplyr::bind_rows(df_1, df_n)
    }
  }
  
  n_cores <- min(parallel::detectCores() - 1L, len)
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  df <- suppressWarnings(
    parallel::clusterApply(
      cl, cols, function(col) {
        df_n[, c(cols123, col)] |> 
          dplyr::summarise_all(function (x) median(x, na.rm = TRUE)) |>
          dplyr::ungroup()
      }
    )
  )
  parallel::stopCluster(cl)

  keys <- df[[1]][, cols123, drop = FALSE]
  xs <- lapply(df, function (x) x[, -cols123, drop = FALSE])
  
  df <- dplyr::bind_cols(keys, xs) |>
    dplyr::bind_rows(df_1)
  
  invisible(df)
}


#' Pads sample groups.
#'
#' For uses with group searches (with non-trivial values under
#' \code{pep_group}). For example, heavy may be missing in complete under one
#' sample and light missing in another.
#'
#' @param df A data frame.
#' @param sids The Sample_IDs from label_scheme.
#' @param tmt_levels The levels of TMT: 000, 126, 127 etc. without the
#'   \code{TMT-} prefix.
#' @param type_levels The levels of five types.
#' 
#' @examples 
#' \donttest{
#' # 10-plex
#' tmt_levels <- TMT_levels(10)
#' tmt_levels <- gsub("^TMT-", "", tmt_levels)
#' 
#' group_ids <- c("light", "heavy")
#' 
#' sids <- LETTERS[1:10]
#' len <- length(sids)
#' sids <- rep(sids, each = length(group_ids))
#' 
#' group_ids <- rep(group_ids, len)
#' sids <- paste0(sids, " [", group_ids, "]")
#' 
#' rm(list = c("len", "group_ids"))
#' 
#' # 1-plex
#' #' tmt_levels <- TMT_levels(1)
#' tmt_levels <- gsub("^TMT-", "", tmt_levels)
#' 
#' group_ids <- c("light", "heavy")
#' 
#' sids <- LETTERS[1]
#' len <- length(sids)
#' sids <- rep(sids, each = length(group_ids))
#' 
#' group_ids <- rep(group_ids, len)
#' sids <- paste0(sids, " [", group_ids, "]")
#' 
#' rm(list = c("len", "group_ids"))
#' }
pad_grp_samples <- function (df, sids, tmt_levels = NULL, 
                             type_levels = c("I", "N_I", "sd_log2_R", "log2_R", "N_log2_R")) 
{
  len_tp <- length(type_levels)
  len_si <- length(sids)
  len_tmt <- length(tmt_levels)

  univ <- unlist(lapply(type_levels, function (x) paste0(x, tmt_levels)))
  
  # pept_tot_int etc.: len_tmt == 0L and is.null(tmt_levels)
  if (len_tmt <= 1L) {
    univ <- unlist(lapply(univ, function (x) paste0(x, " (", sids, ")")))
  }
  else {
    n_grps <- len_si/len_tmt
    univ <- rep(univ, each = n_grps)
    univ <- paste0(univ, " (", rep(sids, len_tp), ")")
  }
  
  more_nms <- univ[!univ %in% names(df)]
  len <- length(more_nms)
  
  if (len) {
    for (i in 1:len) {
      more_nms_i <- more_nms[[i]]
      df[[more_nms_i]] <- NA_real_
    }
    
    nms <- names(df)
    
    # stopifnot(all(univ %in% nms))
    
    df <- dplyr::bind_cols(
      df[, setdiff(nms, univ), drop = FALSE], 
      df[, univ, drop = FALSE])
  }
  
  invisible(df)
}


#' Replicates new label_scheme with new Sample_IDs by group_ids.
#' 
#' The new Sample_ID: "Sample [heavy]", "Sample [light]" etc.
#' 
#' @inheritParams newColnames
rep_ls_groups <- function (group_ids) 
{
  label_scheme <- load_ls_group(dat_dir, label_scheme, prefer_group = FALSE)

  nrow <- nrow(label_scheme)
  sids <- label_scheme$Sample_ID
  n_grps <- length(group_ids)
  row_ids <- rep(1:nrow(label_scheme), each = n_grps)
  sids_grps <- paste0(rep(sids, each = n_grps), " [", group_ids, "]")
  
  label_scheme$Sample_ID_Orig <- label_scheme$Sample_ID
  label_scheme <- label_scheme[row_ids, ]
  label_scheme$Sample_ID <- sids_grps
  label_scheme$Sample_ID_Group <- rep(group_ids, nrow)

  invisible(label_scheme)
}


#' Combines peptide reports across multiple experiments.
#'
#' Median summary of data from the same TMT or LFQ experiment at different LCMS
#' injections summed \code{pep_n_psm}, \code{prot_n_psm}, and \code{prot_n_pep}
#' after data merging no Z_log2_R yet available use \code{col_select =
#' expr(Sample_ID)} not \code{col_select} to get all Z_log2_R why: users may
#' specify \code{col_select} only partial to Sample_ID entries.
#'
#' @param engine The name of search engine.
#' @param use_mq_pep Logical; if TRUE, uses the peptides.txt from MaxQuant.
#' @param use_mf_pep Logical; if TRUE, uses the peptides.txt from MSFragger.
#' @param max_mbr_fold The maximum absolute fold change in MBR.
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @inheritParams splitPSM_ma
#' @inheritParams mergePep
normPep_Mplex <- function (group_psm_by = "pep_seq_mod", 
                           group_pep_by = "prot_acc", 
                           engine = "mz", lfq_mbr = TRUE, 
                           use_duppeps = TRUE, duppeps_repair = "denovo", 
                           cut_points = Inf, omit_single_lfq = FALSE, 
                           use_mq_pep = FALSE, use_mf_pep = FALSE, 
                           rm_allna = FALSE, mbr_ret_tol = 25, max_mbr_fold = 4, 
                           ret_sd_tol = Inf, rm_ret_outliers = FALSE, ...) 
{
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  pat <- paste0("TMTset[0-9]+_LCMSinj[0-9]+_Peptide_N\\.txt$")
  
  filelist  <- list.files(path = file.path(dat_dir, "Peptide"), pattern = pat, 
                          full.names = TRUE)
  basenames <- basename(filelist)
  n_files   <- length(filelist)
  
  if (!length(filelist)) {
    stop("No individual peptide tables available; run `PSM2Pep()` first.")
  }
  
  ok_mbr <- if (TMT_plex || n_files == 1L || engine != "mz") FALSE else TRUE
  ok_mbr <- if (ok_mbr && lfq_mbr) TRUE else FALSE
  
  ###
  # add row means before imputing...
  # separate handling of single value rows?
  
  # 1. Missing references:
  # 1.1 All missing
  # 1.2 row means before imputation
  # 1.3 reference becomes row mean
  # 
  # 2. With references
  # 2.1 All samples missing, both reference and and samples NA
  # 2.2 Impute without reference
  ###
  
  if (ok_mbr) {
    ms1files   <- list.files(dat_dir, pattern = "^ms1full_.*\\.rds$")
    n_ms1files <- length(ms1files)

    if (!n_ms1files) {
      # stop: otherwise df is unique by pep_seq_modz and need to be collapsed
      # to pep_seq_mod
      stop("MS1 peak lists of `^ms1full_[...].rds` not found for MBR.\n", 
              "Please copy them from the Mzion folder ", 
              "to your working directory.")
      ok_mbr <- FALSE
    }
    else if (n_ms1files != n_files) {
      stop("The number of MS1 `^ms1full_[...].rds` files is different ", 
              "to the number of peptide ", pat, " files.")
      ok_mbr <- FALSE
    }
    
    if (ok_mbr) {
      hpeptideMBR(ms1files = ms1files, filelist = filelist, dat_dir = dat_dir, 
                  mbr_ret_tol = mbr_ret_tol, max_mbr_fold = max_mbr_fold, 
                  step = 1e-5)
    }
  }

  if (ok_mbr) {
    df <- lapply(filelist, addMBRpeps, dat_dir = dat_dir, 
                 group_psm_by = group_psm_by)
  }
  else {
    df <- suppressWarnings(
      lapply(filelist, function (x) {
        df <- readr::read_tsv(x, col_types = get_col_types(), 
                              show_col_types = FALSE)

        # MaxQuant all groups contain only single ID and become integer
        if ("Protein Group Ids" %in% names(df)) {
          df <- df |>
            dplyr::mutate(`Protein Group Ids` = as.character(`Protein Group Ids`))
        }
        
        df
      })
    )
  }
  
  df <- dplyr::bind_rows(df)
  df$pep_tot_int <- NULL

  df <- df[, !names(df) %in% c("dat_file")] |>
    dplyr::arrange(TMT_Set) |>
    filters_in_call(!!!filter_dots)

  if ("gene" %in% names(df)) {
    df <- dplyr::mutate(df, gene = forcats::fct_na_value_to_level(gene))
  }
  
  df <- 
    assign_duppeps(df, group_psm_by, group_pep_by, use_duppeps, duppeps_repair)
  
  if (!(use_mq_pep || use_mf_pep)) {
    # make `I000 (SID_1)`, `I000 (SID_2)`, ...
    df_num <- spreadPepNums(df, basenames = basenames, 
                            group_psm_by = group_psm_by, ok_mbr = ok_mbr)

    # (updates metadata in case of "group searches")
    # label_scheme <- load_ls_group(dat_dir, label_scheme)
    # label_scheme_full <- load_ls_group(dat_dir, label_scheme_full)

    if (!TMT_plex) {
      df_num <- calcLFQPepNums(df_num)
    }
  } 
  else {
    use_lowercase_aa <- match_call_arg("normPSM", "use_lowercase_aa")
    
    # back-fill from a MaxQuant peptide table
    if (use_mq_pep) {
      df_num <- pep_mq_lfq(label_scheme, omit_single_lfq)
    }
    else if (use_mf_pep) {
      df_num <- pep_mf_lfq(label_scheme, omit_single_lfq, 
                           group_psm_by = group_psm_by, 
                           use_lowercase_aa = use_lowercase_aa)
    }
    else {
      stop("Unhandled condition.")
    }
  }
  
  df_num <- df_num %>% 
    dplyr::mutate(mean_lint = 
                    log10(rowMeans(.[, grepl("^N_I[0-9]{3}[NC]{0,1}", names(.)), 
                                     drop = FALSE], 
                                   na.rm = TRUE)), 
                  mean_lint = round(mean_lint, digits = 2))
  
  count_nna <- df_num %>% 
    dplyr::select(grep("N_log2_R[0-9]{3}[NC]{0,1}", 
                       names(.))) %>% 
    dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Ref\\.[0-9]+\\)$", 
                        names(.))) %>% 
    dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Empty\\.[0-9]+\\)$", 
                        names(.))) %>% 
    is.na() %>% 
    magrittr::not() %>% 
    rowSums()
  
  df_num <- dplyr::bind_cols(count_nna = count_nna, df_num) |>
    reloc_col_before("mean_lint", "count_nna")

  df$pep_seq_modz <- NULL
  
  # Not for Mzion...
  if ("pep_ret_range" %in% names(df) && !all(is.na(df$pep_ret_range))) {
    if (rm_ret_outliers) {
      df <- df %>% dplyr::mutate(id. = row_number())
      
      oks <- local({
        df <- df[, c(group_psm_by, "pep_ret_range", "id.")]
        col <- which(names(df) == "pep_ret_range")
        df <- df %>% split(.[[group_psm_by]], drop = TRUE)
        
        rows <- unlist(lapply(df, function (x) nrow(x) <= 2L))
        df0 <- dplyr::bind_rows(df[rows])
        df1 <- df[!rows] 
        
        if (length(df1)) {
          df1 <- lapply(df1, locate_outliers, col) %>% 
            dplyr::bind_rows() %>% 
            dplyr::filter(!is.na(pep_ret_range))
        }
        
        c(df0$id., df1$id.)
      })
      
      if (length(oks)) 
        df <- df %>% dplyr::filter(id. %in% oks)
      
      df <- df %>% dplyr::select(-id.)
      rm(list = "oks")
    }
    
    pep_ret_sd <- df %>% 
      calc_pep_retsd(group_psm_by, use_unique = FALSE) %>% 
      dplyr::arrange(!!rlang::sym(group_psm_by))
    
    if (!is.infinite(ret_sd_tol)) {
      message("Removal of `", group_psm_by, "` entries at ", 
              "retention time SD >= ", ret_sd_tol, ".")
      
      # (all entries under the same peptide will be removed)
      df <- local({
        peps <- pep_ret_sd %>% 
          dplyr::filter(pep_ret_sd <= ret_sd_tol) %>% 
          .[[group_psm_by]]
        
        df %>% dplyr::filter(!!rlang::sym(group_psm_by) %in% peps)
      })
    }
  }
  else {
    if (rm_ret_outliers) {
      message("Peptide retention times not available for data filtration ", 
              "by `rm_ret_outliers`.")
    }
    
    if (!is.infinite(ret_sd_tol)) {
      message("Peptide retention times not available for data filtration ", 
              "by `ret_sd_tol`.")
    }

    pep_ret_sd <- df %>%
      dplyr::select(group_psm_by) %>% 
      unique() %>% 
      dplyr::mutate(pep_ret_sd = NA_real_) %>% 
      dplyr::arrange(!!rlang::sym(group_psm_by))
  }
  
  # up to this point pep_tot_int == I000;
  # following the summary statistics, pep_tot_int is the sum across all samples
  df_first <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("prot_matches_sig", "prot_sequences_sig", 
                                         "pep_ret_sd", "pep_exp_mz", "pep_exp_z", 
                                         "TMT_Set"))) %>% 
    # `pep_tot_int` will be summed over samples; 
    # the field is different to `pep_tot_int (SID)`
    med_summarise_keys(group_psm_by) %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by)) %>% 
    reloc_col_after("pep_expect", "pep_score") %>% 
    reloc_col_after("pep_phospho_locprob", "pep_locdiff") %>% 
    reloc_col_after("pep_phospho_locdiff", "pep_phospho_locprob") %>% 
    reloc_col_after("pep_n_nl", "pep_rank_nl")

  # colnm_before <- find_preceding_colnm(df, "pep_tot_int")
  
  df <- list(df_first, pep_ret_sd, df_num) |>
    purrr::reduce(dplyr::left_join, by = group_psm_by) |>
    reloc_col_before(group_psm_by, "pep_res_after") |>
    reloc_col_after("pep_ret_sd", "pep_ret_range")
  
  if (separate_lfq_cols <- TRUE) {
    df <- df[, !grepl("^pep_.*_int \\(", names(df))]
  }

  # --- update sd_log2R000 (...) ---
  if (!(TMT_plex || use_mq_pep || use_mf_pep)) {
    df <- local({
      df <- df[, !grepl("^sd_log2_R000 ", names(df))]

      df_sds <- df %>% 
        calcSD_Splex(group_pep_by) %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.)))
      
      dplyr::right_join(df_sds, df, by = group_pep_by) |>
        na_zeroIntensity()
    })
  }
  
  # --- add varmod columns ---
  if (("pep_seq_mod" %in% names(df)) && 
      (match_call_arg(normPSM, use_lowercase_aa))) {
    df <- df |>
      dplyr::mutate(pep_mod_protnt = ifelse(grepl("^~", pep_seq_mod), 
                                            TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_protntac = ifelse(grepl("^_", pep_seq_mod), 
                                              TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_pepnt = ifelse(grepl("^[_~]?\\^", pep_seq_mod), 
                                           TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_m = ifelse(grepl("m", pep_seq_mod), 
                                       TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_n = ifelse(grepl("n", pep_seq_mod), 
                                       TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_sty = ifelse(grepl("[sty]", pep_seq_mod), 
                                         TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_pepct = ifelse(grepl("[\\^]{1}[_~]?$",
                                                 pep_seq_mod), 
                                           TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_protctam = ifelse(grepl("_{1}$", pep_seq_mod), 
                                              TRUE, FALSE)) |>
      dplyr::mutate(pep_mod_protct = ifelse(grepl("~{1}$", pep_seq_mod), 
                                            TRUE, FALSE))
  } 
  else {
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
  nms <- names(df)
  df <- dplyr::bind_cols(
    df |> dplyr::select(grep("^prot_", nms)),
    df |> dplyr::select(grep("^pep_", nms)), 
    df |> dplyr::select(-grep("^prot_|^pep_", nms)), 
  )
  
  nms <- names(df)
  df <- dplyr::bind_cols(
    df |> dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", nms)), 
    df |> dplyr::select(grep("^I[0-9]{3}[NC]*", nms)), 
    df |> dplyr::select(grep("^N_I[0-9]{3}[NC]*", nms)), 
    df |> dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]*", nms)), 
    df |> dplyr::select(grep("^log2_R[0-9]{3}[NC]*", nms)), 
    df |> dplyr::select(grep("^N_log2_R[0-9]{3}[NC]*", nms)), 
    df |> dplyr::select(grep("^Z_log2_R[0-9]{3}[NC]*", nms)), 
  )

  df <- df %>%
    dplyr::filter(!duplicated(.[[group_psm_by]])) %>% 
    { if (TMT_plex && rm_allna) 
      .[rowSums(!is.na(.[grepl("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))])) > 0, ] 
      else . } %>% 
    dplyr::arrange(!!rlang::sym(group_psm_by))
  
  # a placeholder so no need to handle the exception of 
  # no `Z_log2_R` columns before the first `normMulGau`
  if (!length(grep("^Z_log2_R[0-9]{3}[NC]{0,1}", names(df)))) {
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
    is_prot_lfq = FALSE,
  ) %>% 
    fmt_num_cols() %T>% 
    write.table(file.path(dat_dir, "Peptide/Peptide.txt"), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
}


#' Helper of \link{peptideMBR}
#'
#' @param ms1files The names of MS1 data files.
#' @param filelist The name of peptide files.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param max_mbr_fold The maximum absolute fold change in MBR.
#' @param step The mass error in \code{ppm / 1e6}.
#' @param dat_dir The working directory.
hpeptideMBR <- function (ms1files, filelist, dat_dir, mbr_ret_tol = 25, 
                         max_mbr_fold = 4, step = 1e-5)
{
  # (1) grouped split of ms1files in correspondence to filelist (peptide tables)
  b_nms <- basename(filelist)
  
  ms1files <- local({
    load(file.path(dat_dir, "fraction_scheme.rda"))
    
    fraction_scheme <- fraction_scheme |>
      tidyr::unite(tmt_lcms, c("TMT_Set", "LCMS_Injection"), remove = FALSE)
    
    meta <- data.frame(
      filelist = filelist, 
      basename = b_nms, 
      TMT_Set = gsub("^TMTset(\\d+).*", "\\1", b_nms), 
      LCMS_Injection = 
        gsub("^TMTset\\d+_LCMSinj(\\d+)_Peptide_N.txt$", "\\1", b_nms)) |>
      tidyr::unite(tmt_lcms, c("TMT_Set", "LCMS_Injection"), remove = TRUE) |>
      dplyr::left_join(fraction_scheme, by = "tmt_lcms")
    
    readr::write_tsv(meta, file.path(dat_dir, "Peptide/cache/mbr_meta.txt"))
    
    # raws under each peptide table
    raws <- lapply(filelist, function (file) {
      b_nm <- basename(file)
      rows <- meta$basename == b_nm
      raws <- unique(meta$RAW_File[rows])
      raws <- gsub("\\.raw$", "", raws)
      raws <- gsub("\\.d$", "", raws)
    })
    
    # group ms1files by raws
    ms1raws <- gsub("^ms1full_(.*)\\.rds$", "\\1", ms1files)
    ms1raws <- gsub("\\.raw$", "", ms1raws)
    ms1raws <- gsub("\\.d$", "", ms1raws)
    ms1raws <- lapply(raws, function (x) ms1raws[ms1raws %in% x])
    
    # groups of ms1raws in accordance to peptide tables
    ms1file1 <- ms1files[[1]]
    
    suffix <- if (grepl("\\.raw\\.rds", ".raw.rds", ms1file1)) 
      ".raw.rds"
    else if (grepl("\\.d\\.rds", ".raw.rds", ms1file1))
      ".d.rds"
    else
      ".raw.rds"
    
    ms1files <- lapply(ms1raws, function (x) paste0("ms1full_", x, suffix))
  })
  
  # (2) MBR candidates
  dfs <- suppressWarnings(
    lapply(filelist, function (x) {
      df <- readr::read_tsv(x, col_types = get_col_types(), 
                            show_col_types = FALSE)
      # df <- df[!duplicated(df$pep_seq_modz), ] # just in case
    }))
  pep_seq_modzs <- lapply(dfs, `[[`, "pep_seq_modz")
  pep_apex_rets <- lapply(dfs, `[[`, "pep_apex_ret")
  pep_exp_mzs   <- lapply(dfs, `[[`, "pep_exp_mz")
  pep_exp_ints  <- lapply(dfs, `[[`, "pep_tot_int")

  len <- length(filelist)
  if (len < 2L) {
    warning("Requires at least two files for MBR.")
    return(NULL)
  }
  
  mbr_peps <- mbr_rets <- mbr_mzs <- mbr_ys <- vector("list", len)
  
  for (i in seq_along(filelist)) {
    peps <- unlist(pep_seq_modzs[-i])
    rets <- unlist(pep_apex_rets[-i])
    mzs  <- unlist(pep_exp_mzs[-i])
    ys   <- unlist(pep_exp_ints[-i])
    oks  <- which(!duplicated(peps)) # uses the first among duplicates
    peps <- peps[oks]
    rets <- rets[oks]
    mzs  <- mzs[oks]
    ys   <- ys[oks]
    
    pepsi <- pep_seq_modzs[[i]]
    missings <- which(!peps %in% pepsi)
    mbr_peps[[i]] <- peps[missings]
    mbr_rets[[i]] <- rets[missings]
    mbr_mzs[[i]]  <- mzs[missings]
    mbr_ys[[i]]   <- ys[missings]
  }
  rm(list = c("pep_seq_modzs", "pep_apex_rets", "pep_exp_mzs", 
              "peps", "rets", "mzs", "oks", "missings"))
  
  # (3) MBR
  if ((n_cores <- min(parallel::detectCores(), len)) > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    out <- parallel::clusterMap(
      cl, peptideMBR, 
      ms1files = ms1files, mbr_peps = mbr_peps, mbr_rets = mbr_rets, 
      mbr_mzs = mbr_mzs, mbr_ys = mbr_ys, base_name = b_nms, 
      MoreArgs = list(dat_dir = dat_dir, step = step, mbr_ret_tol = mbr_ret_tol, 
                      max_mbr_fold = max_mbr_fold))
    
    parallel::stopCluster(cl)
  }
  else {
    out <- mapply(
      peptideMBR, ms1files = ms1files, mbr_peps = mbr_peps, mbr_rets = mbr_rets, 
      mbr_mzs = mbr_mzs, mbr_ys = mbr_ys, base_name = b_nms, 
      MoreArgs = list(dat_dir = dat_dir, step = step, mbr_ret_tol = mbr_ret_tol, 
                      max_mbr_fold = max_mbr_fold), 
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  
  out
}


#' MBR of peptides
#'
#' @param ms1files The names of MS1 data files corresponding to a peptide table,
#'   e.g., \code{TMTset1_LCMSinj1_Peptide_N.txt}.
#' @param mbr_peps Sequences of peptides for MBR.
#' @param mbr_rets Values of retention times for MBR.
#' @param mbr_mzs Values of m-over-zs for MBR.
#' @param mbr_ys Intensity values for MBR.
#' @param base_name The basename of a peptide table.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param max_mbr_fold The maximum absolute fold change in MBR.
#' @param step The mass error in \code{ppm / 1e6}.
#' @param dat_dir The working data directory.
peptideMBR <- function (ms1files, mbr_peps, mbr_rets, mbr_mzs, mbr_ys, 
                        base_name, dat_dir, mbr_ret_tol = 25, max_mbr_fold = 4, 
                        step = 1e-5)
{
  df <- lapply(ms1files, function (file) {
    qs::qread(file.path(dat_dir, file))
  }) |>
    dplyr::bind_rows()
  
  tss <- df$ret_time
  xss <- df$msx_moverzs
  yss <- df$msx_ints
  sss <- df$scan_num
  
  rm(list = c("df"))
  
  yints <- rets <- scans <- vector("numeric", length(mbr_peps))

  for (i in seq_along(mbr_peps)) {
    mbr_pep <- mbr_peps[[i]]
    mbr_ret <- mbr_rets[[i]]
    mbr_mz  <- mbr_mzs[[i]]
    mbr_y   <- mbr_ys[[i]]
    
    rng <- tss >= (mbr_ret - 90) & tss <= (mbr_ret + 90)
    rng <- .Internal(which(rng))
    tsj <- tss[rng]
    xsj <- xss[rng]
    ysj <- yss[rng]
    ssj <- sss[rng]
    
    ysub <- mapply(function (xs, ys) ys[abs(xs - mbr_mz) / mbr_mz <= step], 
                   xsj, ysj, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    lens <- lengths(ysub)
    ysub[lens == 0L] <- 0
    rows2 <- .Internal(which(lens > 1L))
    
    if (length(rows2)) {
      ysub[rows2] <- lapply(ysub[rows2], `[[`, 1L)
    }
    
    ysub <- .Internal(unlist(ysub, recursive = FALSE, use.names = FALSE))
    ans <- find_mbr_int(ys = ysub, ts = tsj, ss = ssj, mbr_ret = mbr_ret, 
                        mbr_y = mbr_y, mbr_ret_tol = mbr_ret_tol, 
                        max_mbr_fold = max_mbr_fold, n_dia_scans = 4)
    yints[[i]] <- ans$y
    rets[[i]]  <- ans$t
    scans[[i]] <- ans$s
  }
  
  out <- data.frame(
    pep_seq_modz = mbr_peps, 
    pep_apex_ret = rets, 
    pep_apex_scan = scans, 
    pep_exp_mz = mbr_mzs, 
    pep_tot_int = yints)

  readr::write_tsv(
    out, file.path(dat_dir, "Peptide/cache", paste0("MBR_", base_name)))
  
  out
}


#' Find the MBR intensity
#'
#' @param ys A vector of intensity values.
#' @param ts A vector of time values.
#' @param ss A vector of scan numbers.
#' @param mbr_ret The reference MBR retention time.
#' @param mbr_y The reference MBR intensity.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param max_mbr_fold The maximum absolute fold change in MBR.
#' @param n_dia_scans The maximum number of zero-intensity scans for gap filling
#'   in peak tracing.
find_mbr_int <- function (ys, ts, ss, mbr_ret, mbr_y, mbr_ret_tol = 25, 
                          max_mbr_fold = 4, n_dia_scans = 4)
{
  gates <- mzion:::find_lc_gates(ys = ys, ts = ts, n_dia_scans = n_dia_scans)
  apexs <- gates$apex
  lenp  <- length(apexs)
  
  if (lenp == 0L) {
    return(list(y = 0, t = 0, s = 0))
  }
  
  yints <- gates$yints
  
  upr <- mbr_y * max_mbr_fold
  lwr <- mbr_y / max_mbr_fold
  
  if (lenp == 1L) {
    tval <- ts[[apexs]]
    
    if (abs(tval - mbr_ret) <= mbr_ret_tol) {
      if (yints <= upr && yints >= lwr) {
        return(list(y = yints, t = tval, s = ss[[apexs]]))
      }
      else {
        # return(list(y = mbr_y, t = tval, s = ss[[apexs]]))
        return(list(y = 0, t = 0, s = 0))
      }
    }
    else {
      return(list(y = 0, t = 0, s = 0))
    }
  }
  
  tvals <- ts[apexs]
  scans <- ss[apexs]
  tdiff <- abs(tvals - mbr_ret)
  ydiff <- abs(yints - mbr_y)
  idxt  <- .Internal(which.min(tdiff))
  # idxy  <- which.max(yints)
  idxy  <- .Internal(which.min(ydiff))
  
  if (tdiff[idxy] <= mbr_ret_tol) {
    yi <- yints[idxy]
    
    if (yi <= upr && yi >= lwr) {
      return(list(y = yi, t = tvals[idxy], s = scans[idxy]))
    }
    else {
      return(list(y = 0, t = 0, s = 0))
    }
  }
  
  if (tdiff[idxt] <= mbr_ret_tol) {
    yi <- yints[idxt]
    
    if (yi <= upr && yi >= lwr) {
      return(list(y = yi, t = tvals[idxt], s = scans[idxt]))
    }
    else {
      return(list(y = 0, t = 0, s = 0))
    }
  }
  
  return(list(y = 0, t = 0, s = 0))
}


#' Add MBR peptide intensities
#' 
#' @param file A file name with prepending path to a peptide table.
#' @param dat_dir A working directory.
#' @inheritParams normPSM
addMBRpeps <- function (file, dat_dir, group_psm_by = "pep_seq_mod")
{
  base_name <- basename(file)
  mbr_file  <- file.path(dat_dir, "Peptide/cache", paste0("MBR_", base_name))
  
  cols <- c("pep_seq_modz", "pep_seq_mod", "pep_seq", "pep_exp_mz", "pep_exp_z", 
            "pep_apex_ret", "pep_apex_scan", "pep_tot_int", "pep_score")
  
  pep_tbl <- suppressWarnings(
    readr::read_tsv(file, col_types = get_col_types(), show_col_types = FALSE))
  
  df <- readr::read_tsv(mbr_file) |>
    dplyr::mutate(pep_score = 0) |> 
    tidyr::separate(pep_seq_modz, into = c("pep_seq_mod", "pep_exp_z"), 
                    sep = "@", remove = FALSE) |>
    dplyr::mutate(pep_seq = gsub("^[\\^_~]", "", pep_seq_mod), 
                  pep_seq = gsub("[\\^_~]$", "", pep_seq), 
                  pep_seq = toupper(pep_seq))
  
  if (!all(cols %in% names(df))) {
    stop("Developer: missing columns.")
  }
  
  df <- dplyr::bind_rows(pep_tbl[, cols], df[, cols]) |>
    dplyr::filter(pep_apex_ret > 0) |>
    dplyr::mutate(pep_exp_z = as.integer(pep_exp_z))

  # (2) add `pep_ret_sd`
  df <- df |> 
    dplyr::group_by_at("pep_seq_mod") |> 
    dplyr::mutate(N = dplyr::n())
  
  rows <- df$N == 1L
  df0 <- df[rows, ]
  df <- df[!rows, ]
  df0$N <- df$N <- NULL
  
  df0$pep_ret_sd <- 0
  df <- calc_pep_retsd(df, group_psm_by = "pep_seq_mod") |>
    dplyr::left_join(df, by = "pep_seq_mod")
  df <- df[, names(df0)]
  
  # (3) pep_seq_modz -> pep_seq_mod
  mbr <- groupMZPSM0(df, group_psm_by = group_psm_by)
  mbr <- dplyr::bind_rows(df0[, names(mbr)], mbr)
  out_name <- file.path(dat_dir, "Peptide/cache", paste0("MBRpeps_", base_name))
  readr::write_tsv(mbr, out_name)

  # (4) side effect to return the Peptide.txt
  pep_tbl
}


#' Summary of peptide keys by mean or geomean
#'
#' @param df A data frame.
#' @param ids A vector of column keys being the identifier of the level of
#'   uniqueness. For example, it may be "pep_seq_mod" for non-SILAC or
#'   c("pep_seq_mod", "pep_group") for SILAC with the conversion from PSMs to
#'   peptides.
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
med_summarise_keys <- function(df, ids) 
{
  ## --- Mascot ---
  mascot_median_keys <- c("pep_score", "pep_rank", "pep_isbold", 
                          "pep_exp_mr", "pep_delta", 
                          "pep_exp_mz", "pep_exp_z", 
                          "pep_locprob", "pep_locdiff", 
                          "pep_ret_range", "pep_n_nl", "pep_n_exp_z", 
                          # "pep_score_co", 
                          "pep_n_nl")
  mascot_geomean_keys <- c("pep_expect")
  mascot_sum_keys <- c("pep_tot_int")
  
  df_mascot_med <- df %>% 
    dplyr::select(ids, which(names(.) %in% mascot_median_keys)) %>% 
    dplyr::group_by_at(ids) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mascot_geomean <- df %>% 
    dplyr::select(ids, which(names(.) %in% mascot_geomean_keys)) %>% 
    dplyr::group_by_at(ids) %>% 
    dplyr::summarise_all(~ my_geomean(.x, na.rm = TRUE))
  
  df_all_sum <- df %>% 
    dplyr::select(ids, which(names(.) %in% mascot_sum_keys)) %>% 
    dplyr::group_by_at(ids) %>% 
    dplyr::summarise_all(~ sum(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(mascot_median_keys, 
                                         mascot_geomean_keys, 
                                         mascot_sum_keys)))
  
  ## --- MaxQuant ---
  df_mq_rptr_mass_dev <- df %>% 
    dplyr::select(ids, grep(str_to_title("^Reporter mass deviation"), names(.))) %>% 
    dplyr::group_by_at(ids) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep(str_to_title("^Reporter mass deviation"), names(.)))
  
  mq_median_keys <- str_to_title(c(
    "Score", 
    "Charge", "Mass", "PIF", "Fraction of total spectrum", "Mass error [ppm]", 
    "Mass error [Da]", "Base peak fraction", # "Precursor Intensity", 
    "Precursor Apex Fraction", "Intensity coverage", "Peak coverage", 
    "Combinatorics"
  ))
  mq_geomean_keys <- c("PEP")
  
  df_mq_med <- df %>% 
    dplyr::select(ids, which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by_at(ids) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mq_geomean <- df %>% 
    dplyr::select(ids, which(names(.) %in% mq_geomean_keys)) %>% 
    dplyr::group_by_at(ids) %>% 
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
    df <- df %>% dplyr::mutate(pep_miss = `Missed Cleavages`)
  }
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% str_to_title(
      c("Scan number", "Scan index", "Length", "Missed Cleavages")
    )))
  
  if (all(c("Modifications", "Retention Time") %in% names(df))) {
    df <- df %>% dplyr::select(!`Modifications`:`Retention Time`)
  }
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% str_to_title(
      c("Delta score", "Score diff", "Localization prob", 
        "Precursor full scan number", "Precursor apex fraction",
        "Precursor apex offset", "Precursor apex offset time", 
        "Matches", 
        "Mass deviations [ppm]", "Masses", "Number of matches", 
        "Neutral loss level", "Intensities", 
        "Mass deviations [Da]", 
        "ETD identification type", "Reverse", "All scores", 
        "All sequences", 
        "All modified sequences", 
        "Reporter PIF", "Reporter fraction", 
        "ID", # "Protein group IDs", 
        "Peptide ID", "Mod. peptide ID", 
        "Evidence ID", "Diagnostic Peak Phospho (Sty) Y")
    ))) %>% 
    dplyr::select(-grep("Site IDs$", names(.)))
  
  ## --- Spectrum Mill ---
  sm_median_keys <- c(
    "score", "parent_charge", 
    "deltaForwardReverseScore", "percent_scored_peak_intensity", "totalIntensity", 
    "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
    "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
    "delta_parent_mass", "delta_parent_mass_ppm")
  sm_geomean_keys <- NA
  
  df_sm_med <- df %>% 
    dplyr::select(ids, which(names(.) %in% sm_median_keys)) %>% 
    dplyr::group_by_at(ids) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% sm_median_keys))
  
  ## --- MSFragger ---
  mf_median_keys <- c("Nextscore", "PeptideProphet Probability")
  mf_geomean_keys <- NA
  
  df_mf_med <- df %>% 
    dplyr::select(ids, which(names(.) %in% mf_median_keys)) %>% 
    dplyr::group_by_at(ids) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% mf_median_keys)) 
  
  ## --- put together ---
  df_first <- df %>% 
    tidyr::unite(ids., ids, sep = ".", remove = FALSE) %>% 
    dplyr::filter(!duplicated(ids.)) %>% 
    dplyr::select(-ids.)
  
  df <- list(df_first, 
             df_mascot_med, df_mascot_geomean, df_all_sum, 
             df_mq_rptr_mass_dev, df_mq_med, df_mq_geomean, 
             df_sm_med, df_mf_med) %>%
    purrr::reduce(dplyr::left_join, by = ids) %>%
    data.frame(check.names = FALSE)
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  ) 
}


#' Loads prior Peptide.txt
#' 
#' @inheritParams info_anal
load_prior <- function(filename, id) 
{
  dat_dir <- get_gl_dat_dir()
  
  stopifnot(file.exists(filename))
  
  df <- read.csv(filename, check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") 

  if (! id %in% names(df)) {
    try(unlink(file.path(dat_dir, "Peptide/Peptide.txt")))
    try(unlink(file.path(dat_dir, "Protein/Protein.txt")))
    stop("`Peptide.txt` deleted as column `", id, "` not available.", 
         call. = FALSE)
  }
  
  df <- df %>% dplyr::arrange(!!rlang::sym(id))
}


#' Formats numeric columns
#' 
#' @inheritParams info_anal
fmt_num_cols <- function (df) 
{
  df %>% 
    dplyr::mutate_at(vars(grep("[IR][0-9]{3}[NC]{0,1}", names(.))), 
                     as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]{0,1}", names(.))), 
                     ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("R[0-9]{3}[NC]{0,1}", names(.))), 
                     ~ round(.x, digits = 3)) 
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
#'Description of the column keys in the output: \cr \code{system.file("extdata",
#'"peptide_keys.txt", package = "proteoQ")}
#'
#'The peptide counts in individual peptide tables,
#'\code{TMTset1_LCMSinj1_Peptide_N.txt} etc., may be fewer than the entries
#'indicated under the \code{prot_n_pep} column after the peptide
#'removals/cleanups using \code{purgePSM}.
#'
#'@param duppeps_repair Not currently used (or only with \code{majority}).
#'  Character string; the method of reparing double-dipping peptide sequences
#'  upon data pooling.
#'
#'  For instance, the same sequence of PEPTIDE may be assigned to protein
#'  accession PROT_ACC1 in data set 1 and PROT_ACC2 in data set 2. At the
#'  \code{denovo} default, the peptide to protein association will be
#'  re-established freshly. At the \code{majority} alternative, a majority rule
#'  will be applied for the re-assignments.
#'@param cut_points A named, numeric vector defines the cut points (knots) for
#'  the median-centering of \code{log2FC} by sections. For example, at
#'  \code{cut_points = c(mean_lint = seq(4, 7, .5))}, \code{log2FC} will be
#'  binned according to the intervals of \eqn{-Inf, 4, 4.5, ..., 7, Inf} under
#'  column \code{mean_lint} (mean log10 intensity) in the input data. The
#'  default is \code{cut_points = Inf}, or equivalently \code{-Inf}, where the
#'  \code{log2FC} under each sample will be median-centered as one piece. See
#'  also \code{\link{prnHist}} for data binning in histogram visualization.
#'@param use_duppeps Logical; if TRUE, re-assigns double/multiple dipping
#'  peptide sequences to the most likely proteins by majority votes.
#'@param ret_sd_tol Numeric; the tolerance in the variance of retention time
#'  (w.r.t. measures in seconds). The thresholding applies to TMT data. The
#'  default is \code{Inf}. Depends on the setting of LCMS gradients, a setting
#'  of, e.g., 150 might be suitable.
#'@param rm_ret_outliers Logical; if TRUE, removes peptide entries with outlying
#'  retention times across samples and/or LCMS series.
#'@param omit_single_lfq Logical; if TRUE, omits LFQ entries with single
#'  measured values across all samples. The default is FALSE.
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
#'@inheritParams splitPSM_ma
#'
#'@seealso \emph{Metadata} \cr \code{\link{load_expts}} for metadata preparation
#'  and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr \code{\link{normPSM}} for extended examples in
#'  PSM data normalization \cr \code{\link{PSM2Pep}} for extended examples in
#'  PSM to peptide summarization \cr \code{\link{mergePep}} for extended
#'  examples in peptide data merging \cr \code{\link{standPep}} for extended
#'  examples in peptide data normalization \cr \code{\link{Pep2Prn}} for
#'  extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization.
#'  \cr \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples
#'  in data purging \cr \code{\link{pepHist}} and \code{\link{prnHist}} for
#'  extended examples in histogram visualization. \cr \code{\link{extract_raws}}
#'  and \code{\link{extract_psm_raws}} for extracting MS file names \cr
#'
#'  \emph{Variable arguments of `filter_...`} \cr \code{\link{contain_str}},
#'  \code{\link{contain_chars_in}}, \code{\link{not_contain_str}},
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}},
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings
#'  \cr
#'
#'  \emph{Missing values} \cr \code{\link{pepImp}} and \code{\link{prnImp}} for
#'  missing value imputation \cr
#'
#'  \emph{Informatics} \cr \code{\link{pepSig}} and \code{\link{prnSig}} for
#'  significance tests \cr \code{\link{pepVol}} and \code{\link{prnVol}} for
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
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output is in \code{.../Peptide/Peptide.txt}.
#'
#'@example inst/extdata/examples/mergePep_.R
#'@import stringr dplyr tidyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
mergePep <- function (plot_log2FC_cv = TRUE, use_duppeps = TRUE, 
                      duppeps_repair = c("majority", "denovo"), 
                      cut_points = Inf, rm_allna = FALSE, 
                      omit_single_lfq = FALSE, ret_sd_tol = Inf, 
                      rm_ret_outliers = FALSE, ...) 
{
  dat_dir <- get_gl_dat_dir()
  engine <- find_search_engine(dat_dir)
  
  old_opts <- options()
  options(warn = 1L)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = environment())) {
      if (.savecall) {
        mget(names(formals()), envir = environment(), inherits = FALSE) %>% 
          c(dots) %>% 
          save_call("mergePep")
      }
    }, add = TRUE)

  # ---
  duppeps_repair <- "majority"
  duppeps_repair <- rlang::enexpr(duppeps_repair)
  oks <- eval(formals()[["duppeps_repair"]])
  
  duppeps_repair <- if (length(duppeps_repair) > 1L) 
    oks[[1]]
  else 
    rlang::as_string(duppeps_repair)
  
  stopifnot(duppeps_repair %in% oks, length(duppeps_repair) == 1L)
  
  rm(list = c("oks"))
  # ---
  
  stopifnot(vapply(c(plot_log2FC_cv, use_duppeps, rm_allna, omit_single_lfq, 
                     rm_ret_outliers), 
                   rlang::is_logical, logical(1)))
  stopifnot(cut_points >= 0, ret_sd_tol > 0)

  reload_expts()

  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))

  dir.create(file.path(dat_dir, "Peptide/cache"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/Histogram"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/log2FC_cv/raw"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/log2FC_cv/purged"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  mbr_ret_tol  <- match_call_arg(normPSM, mbr_ret_tol)

  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  # later apply also the back-filling with MSFragger data...
  message("Primary column keys in `Peptide/TMTset1_LCMSinj1_Peptide_N.txt` etc. ", 
          "for `filter_` varargs.")
  use_mq_pep <- use_mq_peptable(dat_dir, label_scheme_full)
  use_mf_pep <- use_mf_peptable(dat_dir, label_scheme_full, 
                                group_psm_by = group_psm_by)
  
  if (use_mq_pep && use_mf_pep) {
    stop("Both MaxQuant and MSFragger peptide tables found.", 
         "Use either, not both.")
  }

  df <- normPep_Mplex(group_psm_by = group_psm_by, 
                      group_pep_by = group_pep_by, 
                      engine = engine, 
                      lfq_mbr = match_call_arg(normPSM, lfq_mbr),
                      use_duppeps = use_duppeps, 
                      duppeps_repair = duppeps_repair, 
                      cut_points = cut_points, 
                      omit_single_lfq = omit_single_lfq,
                      use_mq_pep = use_mq_pep, 
                      use_mf_pep = use_mf_pep, 
                      rm_allna = rm_allna, 
                      mbr_ret_tol = mbr_ret_tol, 
                      max_mbr_fold = 4, 
                      ret_sd_tol = ret_sd_tol, 
                      rm_ret_outliers = rm_ret_outliers, 
                      !!!filter_dots) 
  
  if (plot_log2FC_cv) {
    quiet_out <- purrr::quietly(sd_violin)(
      df = df, 
      id = !!group_pep_by, 
      filepath = file.path(dat_dir, "Peptide/log2FC_cv/raw/Peptide_sd.png"), 
      width = 8 * n_TMT_sets(label_scheme), 
      height = 8, 
      type = "log2_R", 
      adjSD = FALSE, 
      is_psm = FALSE, 
      ...
    )
  }
  
  .saveCall <- TRUE
  
  invisible(NULL)
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
#'\code{system.file("extdata", "peptide_keys.txt", package = "proteoQ")}
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
#'  range of the \code{intensity} of reporter ions (including \code{I000}) for
#'  use in the scaling normalization of standard deviation across samples. The
#'  default is between the 5th and the 95th quantiles.
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
#'system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output is in \code{.../Peptide/Peptide.txt}.
#'
#'@example inst/extdata/examples/normPep_.R
#'
#'@import stringr dplyr tidyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
standPep <- function (method_align = c("MC", "MGKernel"), col_select = NULL, 
                      range_log2r = c(10, 90), 
                      range_int = c(5, 95), n_comp = NULL, seed = NULL, 
                      plot_log2FC_cv = FALSE, ...) 
{
  dat_dir <- get_gl_dat_dir()
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        mget(names(formals()), envir = rlang::current_env(), 
             inherits = FALSE) %>% 
          c(dots) %>% 
          save_call("standPep")
      }
    }, 
    add = TRUE
  )
  
  stopifnot(vapply(c(plot_log2FC_cv), rlang::is_logical, logical(1)))
  
  reload_expts()
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  method_align <- rlang::enexpr(method_align)
  if (method_align == rlang::expr(c("MC", "MGKernel"))) {
    method_align <- "MC"
  } 
  else {
    method_align <- rlang::as_string(method_align)
    stopifnot(method_align %in% c("MC", "MGKernel"), 
              length(method_align) == 1L)
  }
  
  range_log2r <- prep_range(range_log2r)
  range_int <- prep_range(range_int)
  
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  label_scheme_full <- load_ls_group(dat_dir, label_scheme_full)
  ok_existing_params(file.path(dat_dir, "Peptide/Histogram/MGKernel_params_N.txt"))

  col_select <- rlang::enexpr(col_select)
  col_select <- if (is.null(col_select)) 
    rlang::expr(Sample_ID) 
  else 
    rlang::sym(col_select)
  col_select <- parse_col_select(rlang::as_string(col_select), label_scheme)

  dots <- rlang::enexprs(...)
  
  filename <- file.path(dat_dir, "Peptide/Peptide.txt")
  
  if (!file.exists(filename)) {
    stop(filename, " not found; run `mergePep(...)` first.", call. = FALSE)
  }

  message("Primary column keys in `Peptide/Peptide.txt` for `slice_` varargs.")

  df <- load_prior(filename, group_psm_by) 
  
  if (sum(grepl("^log2_R[0-9]{3}[NC]{0,1}", names(df))) <= 1) {
    stop("Need more than one sample for `standPep` or `standPrn`.\n", 
         "Skip this module for qualitative analysis.", 
         call. = FALSE)
  }
  
  df <- normMulGau(
    df = df,
    method_align = method_align,
    n_comp = n_comp,
    seed = seed,
    range_log2r = range_log2r,
    range_int = range_int,
    filepath = file.path(dat_dir, "Peptide/Histogram"),
    col_select = col_select, 
    cut_points = Inf, 
    is_prot_lfq = FALSE, 
    !!!dots, 
  ) %>% 
    fmt_num_cols() %T>% 
    write.table(file.path(dat_dir, "Peptide", "Peptide.txt"), 
                sep = "\t", col.names = TRUE, row.names = FALSE)

  if (plot_log2FC_cv & TMT_plex(label_scheme)) {
    sd_violin(df = df, id = !!group_pep_by, 
              filepath = file.path(dat_dir, "Peptide/log2FC_cv/raw", "Peptide_sd.png"), 
              width = 8 * n_TMT_sets(label_scheme), height = 8, 
              type = "log2_R", adjSD = FALSE, is_psm = FALSE)
  }
  
  .saveCall <- TRUE
  
  invisible(NULL)
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
#'  descriptive statistics includes \code{c("mean", "median", "weighted_mean",
#'  "top_3_mean", "lfq_max", "lfq_top_2_sum", "lfq_top_3_sum", "lfq_all")} with
#'  \code{median} being the default for TMT and \code{lfq_top_3_sum} for LFQ.
#'  The representative \code{log10-intensity} of reporter (or LFQ) ions at the
#'  peptide levels will be the weight when summarizing \code{log2FC} with
#'  various \code{"top_n"} statistics or \code{"weighted_mean"}.
#'@param use_unique_pep Logical. If TRUE, only entries that are \code{TRUE} or
#'  equal to \code{1} under the column \code{pep_isunique} in \code{Peptide.txt}
#'  will be used, for summarizing the \code{log2FC} and the \code{intensity} of
#'  peptides into protein values. The default is to use unique peptides only.
#'  For \code{MaxQuant} data, the levels of uniqueness are according to the
#'  \code{pep_unique_by} in \code{\link{normPSM}}. The argument currently do
#'  nothing to \code{Spectrum Mill} data where both unique and shared peptides
#'  will be kept.
#' @param impute_prot_na Logical; impute NA values of protein log2FC or not.
#' @param mc Logical. At the TRUE default, performs median-centering of
#'  \code{log2FC} after the peptide-to-protein aggregation. Otherwise, the
#'  summarized \code{log2FC} values will be left as they are.
#'@param ... \code{filter_}: Variable argument statements for the filtration of
#'  data rows. Each statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key in \code{Peptide.txt}. The
#'  statement of \code{filter_peps_at = exprs(pep_len <= 50)} will remove
#'  peptide entries with \code{pep_len > 50}.
#'@inheritParams mergePep
#'@inheritParams normPSM
#'@seealso \emph{Metadata} \cr \code{\link{load_expts}} for metadata preparation
#'  and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr \code{\link{normPSM}} for extended examples in
#'  PSM data normalization \cr \code{\link{PSM2Pep}} for extended examples in
#'  PSM to peptide summarization \cr \code{\link{mergePep}} for extended
#'  examples in peptide data merging \cr \code{\link{standPep}} for extended
#'  examples in peptide data normalization \cr \code{\link{Pep2Prn}} for
#'  extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization.
#'  \cr \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples
#'  in data purging \cr \code{\link{pepHist}} and \code{\link{prnHist}} for
#'  extended examples in histogram visualization. \cr \code{\link{extract_raws}}
#'  and \code{\link{extract_psm_raws}} for extracting MS file names \cr
#'
#'  \emph{Variable arguments of `filter_...`} \cr \code{\link{contain_str}},
#'  \code{\link{contain_chars_in}}, \code{\link{not_contain_str}},
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}},
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings
#'  \cr
#'
#'  \emph{Missing values} \cr \code{\link{pepImp}} and \code{\link{prnImp}} for
#'  missing value imputation \cr
#'
#'  \emph{Informatics} \cr \code{\link{pepSig}} and \code{\link{prnSig}} for
#'  significance tests \cr \code{\link{pepVol}} and \code{\link{prnVol}} for
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
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output in "\code{.../Protein/Protein.txt}".
#'
#'@example inst/extdata/examples/Pep2Prn_.R
#'@import stringr dplyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
Pep2Prn <- function (method_pep_prn = 
                       c("median", "mean", "weighted_mean", "lfq_top_3_sum", 
                         "lfq_all", "lfq_top_2_sum", "top_3_mean", "lfq_max"), 
                     impute_prot_na = TRUE, use_unique_pep = TRUE, 
                     cut_points = Inf, rm_outliers = FALSE, rm_allna = FALSE, 
                     mc = TRUE, ...) 
{
  dat_dir <- get_gl_dat_dir()
  
  dir.create(file.path(dat_dir, "Protein/Histogram"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/cache"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) |>
          c(dots) |>
          save_call("Pep2Prn")
      }
    }, add = TRUE)

  reload_expts()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  method_pep_prn <- rlang::enexpr(method_pep_prn)
  
  if (TMT_plex) {
    if (length(method_pep_prn) > 1L) 
      method_pep_prn <- "median"
    else 
      method_pep_prn <- rlang::as_string(method_pep_prn)
  } 
  else {
    if (length(method_pep_prn) > 1L) 
      method_pep_prn <- "lfq_all"
    else 
      method_pep_prn <- rlang::as_string(method_pep_prn)
  }
  
  if (method_pep_prn == "top.3") {
    stop("Method `top.3` depreciated; instead use `top_3_mean`.")
  } 
  else if (method_pep_prn == "weighted.mean") {
    stop("Method `weighted.mean` depreciated; instead use `weighted_mean`.")
  }
  
  stopifnot(method_pep_prn %in% c("median", "mean", 
                                  "weighted_mean", "top_3_mean", 
                                  "lfq_max", "lfq_top_2_sum", 
                                  "lfq_top_3_sum", "lfq_all"), 
            length(method_pep_prn) == 1L)

  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(group_pep_by %in% c("prot_acc", "gene"), length(group_pep_by) == 1)
  stopifnot(vapply(c(use_unique_pep, mc), rlang::is_logical, logical(1)))
  
  gn_rollup <- if (group_pep_by == "gene") TRUE else FALSE

  message("Column keys in `Peptide/Peptide.txt` ", "for \"filter_\" varargs.")

  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  engine <- find_search_engine(dat_dir)
  use_mq_prot <- use_mq_prottable(dat_dir) && engine == "mq"
  use_mf_prot <- use_mf_prottable(dat_dir) && engine == "mf"

  if (use_mq_prot && use_mf_prot) {
    stop("Both MaxQuant and MSFragger protein tables found.", 
         "Use either, not both.")
  }
  
  is_prot_lfq <- grepl("^lfq_", method_pep_prn)
  
  df <- pep_to_prn(id = prot_acc, 
                   method_pep_prn = method_pep_prn, 
                   use_unique_pep = use_unique_pep, 
                   gn_rollup = gn_rollup, 
                   rm_outliers = rm_outliers, 
                   rm_allna = rm_allna, 
                   is_prot_lfq = is_prot_lfq, 
                   engine = engine, 
                   impute_prot_na = impute_prot_na, 
                   use_mq_prot = use_mq_prot, 
                   use_mf_prot = use_mf_prot, 
                   !!!filter_dots) 
  
  if (mc) {
    df <- normMulGau(
      df = df,
      method_align = "MC",
      n_comp = 1L,
      range_log2r = c(0, 100),
      range_int = c(0, 100),
      filepath = file.path(dat_dir, "Protein/Histogram"),
      col_select = rlang::expr(Sample_ID), 
      is_prot_lfq = is_prot_lfq, 
      cut_points = cut_points) 
  }
  else {
    warning("No data centering performed.", call. = FALSE)
  }

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
  
  .saveCall <- TRUE
  
  invisible(df)
}


#' Helper: maps peptides to possible proteins
#' 
#' @param df A data frame.
#' @inheritParams normPSM
map_peps_prots <- function (df, group_psm_by = "pep_seq", 
                            group_pep_by = "prot_acc") 
{
  dat_dir <- get_gl_dat_dir()
  
  df <- df |> dplyr::filter(!duplicated(group_psm_by))

  if (group_pep_by == "gene") {
    pep_prot_map <- df |> dplyr::select(shared_genes)
  }
  else if (group_pep_by == "prot_acc") {
    pep_prot_map <- df |> dplyr::select(shared_prot_accs)
  }
  else {
    stop("`group_pep_by` needs to be either `prot_acc` or `gene`.")
  }

  pep_prot_map <- unlist(pep_prot_map, recursive = FALSE, use.names = FALSE) |>
    stringr::str_split(", ") |>
    setNames(df[[group_psm_by]])
  
  save(pep_prot_map, file = file.path(dat_dir, "pep_prot_map.rda"))
  
  invisible(pep_prot_map)
}


#' Helper: finds the row indexes of a protein within a family
#' 
#' Not currently used. The rows include primary (razor) and secondary (mapped)
#' proteins.
#' 
#' @param df A data frame.
#' @inheritParams normPSM
find_prot_family_rows <- function (df, group_psm_by, group_pep_by) 
{
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


#' Sums top-n.
#' 
#' @param x A numeric vector.
#' @param n An positive integer.
my_sum_n <- function (x, n = 3, ...) 
{
  if (n < 1L) 
    stop("\"n\" need to be a positive integer.")
  
  if (is.infinite(n)) 
    n <- length(x)
  
  # need to convert x to integers if to partial sort
  sum(sort(x, decreasing = TRUE)[1:n], ...)
}


#' Adds the \code{total}, \code{razor} and \code{unique} intensities of
#' proteins.
#'
#' @param df A data frame.
#' @param impute_prot_na Logical; impute NA values of protein log2FC or not.
#' @inheritParams normPSM
#' @inheritParams Pep2Prn
calcLFQprnnums <- function (df, use_unique_pep = TRUE, 
                            group_psm_by = "pep_seq", group_pep_by = "gene", 
                            method_pep_prn = "lfq_top_3_sum", 
                            impute_prot_na = TRUE)
{
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  refChannels  <- label_scheme$Sample_ID[label_scheme[["Reference"]]]

  if (use_unique_pep) {
    pep_unique_by <- match_call_arg(normPSM, pep_unique_by)
    
    if (pep_unique_by == "group") {
      df <- df |> dplyr::filter(pep_razor_unique)
    }
    else if (pep_unique_by == "protein") {
      df <- df |> dplyr::filter(pep_literal_unique)
    }
    else {
      stop("`pep_unique_by` need to be `group` or `protein`.")
    }
  }
  
  prot_ints <- df |>
    dplyr::select(group_pep_by, grep("^I[0-9]{3}[NC]{0,1}", names(df)))
  prot_ints <- switch(
    method_pep_prn, 
    mean = aggrNums(mean)(prot_ints, !!rlang::sym(group_pep_by), na.rm = TRUE), 
    median = aggrNums(median)(prot_ints, !!rlang::sym(group_pep_by), na.rm = TRUE),
    top_3_mean = TMT_top_n(prot_ints, !!rlang::sym(group_pep_by), na.rm = TRUE), 
    weighted_mean = tmt_wtmean(prot_ints, !!rlang::sym(group_pep_by), na.rm = TRUE), 
    aggrNums(median)(prot_ints, !!rlang::sym(group_pep_by), na.rm = TRUE))
  
  prot_log2rs <- df |>
    dplyr::select(group_pep_by, grep("^log2_R[0-9]{3}[NC]{0,1}", names(df)))
  prot_log2rs <- 
    aggrNums(median)(prot_log2rs, !!rlang::sym(group_pep_by), na.rm = TRUE)
  
  if (impute_prot_na) {
    cols <- grep("^log2_R[0-9]{3}[NC]{0,1}", names(prot_log2rs))
    
    if (length(cols) > 1L) {
      class(prot_log2rs) <- "data.frame"
      prot_log2rs[, cols] <- 
        impPepNA(prot_log2rs[, cols, drop = FALSE], is_intensity = FALSE)
    }
  }
  
  # e.g. expand the same pep_seq_mod at different prot_acc's -> multiple rows
  if (FALSE) {
    pep_prot_map <- df %>% 
      map_peps_prots(group_psm_by, group_pep_by) %>% 
      list_to_dataframe() %>% 
      `names_pos<-`(1, group_psm_by) %>% 
      tidyr::gather(-group_psm_by, key = n, value = !!rlang::sym(group_pep_by)) %>% 
      dplyr::select(-n) %>% 
      dplyr::filter(!is.na(!!rlang::sym(group_pep_by))) %>% 
      dplyr::rename(prot_map = group_pep_by)

    df <- df |> dplyr::left_join(pep_prot_map, by = group_psm_by)
  }

  # --- median centering ---
  dfw <- prot_ints |> dplyr::left_join(prot_log2rs, by = group_pep_by)
  cols_log2Ratio <- grepl("^log2_R[0-9]{3}[NC]{0,1}", names(dfw))
  
  # cfs <- apply(dfw[, cols_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  cfs <- lapply(dfw[, cols_log2Ratio, drop = FALSE], function (x) {
    median(x[x >= -2 & x <= 2], na.rm = TRUE)
  }) |>
    unlist()
  
  dfw <- sweep(dfw[, cols_log2Ratio, drop = FALSE], 2, cfs, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    dplyr::bind_cols(dfw, .)
  
  dfw <- sweep(dfw[, grepl("^I[0-9]{3}", names(dfw)), drop = FALSE], 
               2, 2^cfs, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    dplyr::bind_cols(dfw, .)
  
  if (!length(grep("log2_R000", names(dfw)))) {
    stop("No \"log2_R000...\" columns available.\n")
  }
  
  dfw <- dfw %>% 
    dplyr::mutate_at(.vars = grep("log2_R000\\s", names(.)), 
                     function (x) replace(x, is.infinite(x) | is.nan(x), 
                                          NA_real_))
  
  if (!length(grep("^Z_log2_R[0-9]{3}[NC]{0,1}", names(dfw)))) {
    dfw <- dfw %>% 
      dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
      `names<-`(gsub("^N_log2_R", "Z_log2_R", names(.))) %>% 
      dplyr::bind_cols(dfw, .)
  }
  
  nms <- names(dfw)
  dfw <- dplyr::bind_cols(
    dfw |> dplyr::select(group_pep_by),
    dfw |> dplyr::select(grep("^I000 \\(", nms)),
    dfw |> dplyr::select(grep("^N_I000 \\(", nms)),
    dfw |> dplyr::select(grep("^log2_R000 \\(", nms)),
    dfw |> dplyr::select(grep("^N_log2_R000 \\(", nms)),
    dfw |> dplyr::select(grep("^Z_log2_R000 \\(", nms)), )
}


#' Helper: calculates the TMT log2FC and reporter-ion intensity of proteins
#' 
#' @param id Always "prot_acc".
#' @inheritParams calcLFQprnnums
#' @inheritParams Pep2Prn
calc_tmt_prnnums <- function (df, use_unique_pep = TRUE, id = "prot_acc", 
                              method_pep_prn = "median") 
{
  if (use_unique_pep && "pep_isunique" %in% names(df)) {
    df <- df %>% dplyr::filter(pep_isunique)
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
  
  invisible(df_num)
}


#' Helper of Pep2Prn
#' 
#' @param gn_rollup Logical; if TRUE, rolls up protein accessions to gene names.
#' @param engine The name of search engine.
#' @param use_mq_prot Logical; use MQ protein table or not.
#' @param use_mf_prot Logical; use MSFragger protein table or not.
#' @inheritParams info_anal
#' @inheritParams Pep2Prn
pep_to_prn <- function(id = "prot_acc", method_pep_prn = "median", 
                       use_unique_pep = TRUE, gn_rollup = TRUE, 
                       rm_outliers = FALSE, rm_allna = FALSE, 
                       is_prot_lfq = FALSE, engine = "mz", impute_prot_na = TRUE, 
                       use_mq_prot = FALSE, use_mf_prot = FALSE, ...)
{
  dat_dir <- get_gl_dat_dir()
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  TMT_plex <- TMT_plex(label_scheme)

  # `id` is always "prot_acc"; `group_pep_by` could be "gene"
  id <- rlang::as_string(rlang::enexpr(id))
  
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  df <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, "Peptide/Peptide.txt"), 
                    col_types = get_col_types(), 
                    show_col_types = FALSE))

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    { if (TMT_plex && rm_allna) 
      .[rowSums(!is.na(.[grepl("^log2_R[0-9]{3}[NC]{0,1}", names(.))])) > 0, ] 
      else . } 
  
  if (rm_outliers) {
    df <- local({
      df$pep_index <- seq_along(1:nrow(df))
      
      dfw_split <- df %>% 
        dplyr::select(!!rlang::sym(id), 
                      pep_index, 
                      grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))) %>% 
        dplyr::group_by(!!rlang::sym(id)) %>%
        `colnames<-`(gsub("^log2_R", "X", names(.))) %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         ~ replace(.x, is.infinite(.x), NA)) %>% 
        data.frame(check.names = FALSE) %>% 
        split(.[[id]], drop = TRUE)
      
      range_colRatios <- grep("^X[0-9]{3}", names(dfw_split[[1]]))
      
      dfw_split <- dfw_split %>% 
        purrr::map(locate_outliers, range_colRatios) %>% 
        dplyr::bind_rows() %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         function (x) replace(x, is.infinite(x), NA_real_)) %>% 
        tidyr::unite(prot_acc_i, !!rlang::sym(id), pep_index, sep = ":") %>%
        dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), 
                         function (x) replace(x, !is.na(x), 1))
      
      df <- df %>% 
        tidyr::unite(prot_acc_i, !!rlang::sym(id), pep_index, sep = ":") %>%
        dplyr::left_join(dfw_split, by = "prot_acc_i") %>%
        tidyr::separate(prot_acc_i, into = c(id, "pep_index"), 
                        sep = ":", remove = TRUE) %>%
        dplyr::select(-pep_index)
    })
    
    df[, grepl("^I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^Z_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df <- df %>% 
      { if (rm_allna) 
        .[rowSums(!is.na(.[grepl("^log2_R[0-9]{3}[NC]{0,1}", names(.))])) > 0, ] 
        else . } %>% 
      dplyr::select(-grep("^X[0-9]{3}[NC]{0,1}", names(.)))
  } 

  if (! "pep_isunique" %in% names(df)) {
    df$pep_isunique <- TRUE
    warning("Column \"pep_isunique\" created and TRUE values assumed.")
  }
  else if (all(is.na(df$pep_isunique))) {
    df$pep_isunique <- TRUE
    warning("Values of \"pep_isunique\" are all NA and coerced to TRUE.")
  }
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  # add `prot_n_uniqpep` and `prot_n_uniqpsm`
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
      dplyr::select(!!rlang::sym(group_pep_by), pep_n_psm) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>% 
      dplyr::summarise(prot_n_sharepsms = sum(pep_n_psm))

    df <- list(df, prot_n_sharepeps, prot_n_sharepsms) %>% 
      purrr::reduce(dplyr::left_join, by = group_pep_by) %>% 
      dplyr::mutate(prot_n_sharepeps = replace(prot_n_sharepeps, 
                                               is.na(prot_n_sharepeps), 0), 
                    prot_n_sharepsms = replace(prot_n_sharepsms, 
                                               is.na(prot_n_sharepsms), 0)) %>% 
      dplyr::mutate(prot_n_uniqpep = prot_n_pep - prot_n_sharepeps, 
                    prot_n_uniqpsm = prot_n_psm - prot_n_sharepsms) %>% 
      dplyr::select(-prot_n_sharepeps, -prot_n_sharepsms)
    
    df <- df %>% 
      ins_cols_after(which(names(.) == "prot_n_psm"), 
                     which(names(.) == "prot_n_uniqpsm")) %>% 
      ins_cols_after(which(names(.) == "prot_n_pep"), 
                     which(names(.) == "prot_n_uniqpep"))
  })
  
  # first by `id = prot_acc` and later optional roll-up to `gene`
  if (is_prot_lfq) {
    if (TMT_plex) {
      stop("\"method_pep_prn = lfq_[...]\" only for LFQ.")
    }

    if (engine == "mz") {
      df_num <- calcLFQprnnums(
        df, use_unique_pep = use_unique_pep, group_psm_by = group_psm_by, 
        group_pep_by = id, method_pep_prn = method_pep_prn, 
        impute_prot_na = impute_prot_na)
    }
    else if (engine == "mf") {
      df_num <- fillMFprnnums(dat_dir)
    }
    else if (engine == "mq") {
      df_num <- fillMQprnnums(dat_dir)
    }
    
    df_num <- df_num |> dplyr::filter(!is.na(!!rlang::sym(id)))
  } 
  else {
    df_num <- calc_tmt_prnnums(df, use_unique_pep, id, method_pep_prn)
  }
  
  df_num <- df_num %>% 
    dplyr::mutate(mean_lint = 
                    log10(rowMeans(.[, grepl("^N_I[0-9]{3}[NC]{0,1}", names(.)), 
                                     drop = FALSE], na.rm = TRUE)), 
                  mean_lint = round(mean_lint, digits = 2L))
  
  nms <- names(df_num)
  count_nna <- df_num |>
    dplyr::select(grep("N_log2_R[0-9]{3}[NC]{0,1}", nms)) |>
    dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Ref\\.[0-9]+\\)$", 
                        nms)) |>
    dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Empty\\.[0-9]+\\)$", 
                        nms)) |>
    is.na() |>
    magrittr::not() |>
    rowSums()
  
  df_num <- dplyr::bind_cols(count_nna = count_nna, df_num) |>
    reloc_col_before("mean_lint", "count_nna")
  
  df <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("pep_istryptic", "pep_semitryptic", 
                                         "mean_lint", "count_nna", 
                                         "shared_prot_accs", "shared_genes"))) %>% 
    dplyr::select(-grep("^Reporter mass deviation", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c("m/z", "PIF", "PEP"))) %>% 
    dplyr::select(-which(names(.) %in% stringr::str_to_title(
      c("Charge", "Mass", "Mass error [ppm]", 
        "Mass error [Da]", "Score", "Combinatorics", 
        "Fraction of total spectrum", 
        "Base peak fraction", 
        "Precursor Intensity", "Precursor intensity", 
        "Precursor Apex Fraction", 
        "Intensity coverage", "Intensity Coverage", 
        "Peak coverage", "Peak Coverage",  
        "Proteins", "Protein group IDs")
    ))) %>% 
    dplyr::select(-which(names(.) %in% c("Is Unique", "Protein", "Gene", 
                                         "Mapped Genes", 
                                         "Mapped Proteins", "Nextscore", 
                                         "PeptideProphet Probability")))

  mq_median_keys <- NULL
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  df <- df %>% dplyr::select(-which(names(.) %in% mq_median_keys))
  rm(list = "mq_median_keys")
  
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
  rm(list = "sm_median_keys")
  
  mf_median_keys <- NULL
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mf_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  df <- df %>% dplyr::select(-which(names(.) %in% mf_median_keys))
  rm(list = "mf_median_keys")

  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
    dplyr::select(-grep("^pep_", names(.)))
  
  df <- list(df_first, 
             df_mq_med, 
             df_sm_med, 
             df_num) %>% 
    purrr::reduce(left_join, by = id) %>% 
    data.frame(check.names = FALSE)
  
  nms <- names(df)
  cols <- grepl("log2_R[0-9]{3}", nms) & !sapply(df, is.logical)
  df[, cols] <- df[, cols] |>
    dplyr::mutate_if(is.integer, as.numeric) |>
    round(digits = 3L)

  cols <- grepl("I[0-9]{3}", nms) & !sapply(df, is.logical)
  df[, cols] <- df[, cols] |>
    dplyr::mutate_if(is.integer, as.numeric) |>
    round(digits = 0L)
  
  df <- dplyr::bind_cols(
    df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", nms), drop = FALSE], 
    df[, grep("^I[0-9]{3}", nms), drop = FALSE], 
    df[, grep("^N_I[0-9]{3}", nms), drop = FALSE], 
    df[, grep("^log2_R[0-9]{3}", nms), drop = FALSE], 
    df[, grep("^N_log2_R[0-9]{3}", nms), drop = FALSE], 
    df[, grep("^Z_log2_R[0-9]{3}", nms), drop = FALSE])
  
  if (rm_allna) {
    df <- df %>% 
      .[rowSums(!is.na(.[grepl("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))])) > 0, ]
  }

  if (gn_rollup) {
    # don't move: for keeping track of the original column names
    nms <- names(df)
    
    df$mean_lint <- df$count_nna <- NULL

    dfa <- local({
      dfa <- df %>% 
        dplyr::select(gene, grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
        dplyr::filter(!is.na(gene)) %>% 
        dplyr::group_by(gene) %>% 
        dplyr::summarise_all(list(function (x) median(x, na.rm = TRUE)))
      
      dfa <- dfa %>% 
        dplyr::mutate(
          mean_lint = log10(rowMeans(.[, grepl("^N_I[0-9]{3}[NC]{0,1}", names(.)), 
                                       drop = FALSE], na.rm = TRUE)), 
          mean_lint = round(mean_lint, digits = 2L))

      count_nna <- dfa %>% 
        dplyr::select(grep("N_log2_R[0-9]{3}[NC]{0,1}", 
                           names(.)))%>% 
        dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Ref\\.[0-9]+\\)$", 
                            names(.))) %>% 
        dplyr::select(-grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s\\(Empty\\.[0-9]+\\)$", 
                            names(.))) %>% 
        is.na() %>% 
        magrittr::not() %>% 
        rowSums() 
      
      dfa <- dplyr::bind_cols(count_nna = count_nna, dfa) %>% 
        reloc_col_before("mean_lint", "count_nna")
    })

    dfb <- df %>% 
      dplyr::select(-prot_cover, 
                    -grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
      dplyr::filter(!is.na(gene)) %>% 
      dplyr::filter(!duplicated(gene))
    
    dfc <- suppressWarnings(
      df %>% 
        dplyr::select(gene, prot_cover) %>% 
        dplyr::filter(!is.na(gene)) %>% 
        dplyr::group_by(gene) %>% 
        dplyr::summarise_all(~ max(.x, na.rm = TRUE)))

    df <- list(dfc, dfb, dfa) %>% 
      purrr::reduce(right_join, by = "gene") %>% 
      dplyr::filter(!is.na(gene), !duplicated(gene)) %>% 
      dplyr::select(nms)
  }
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(-grep("^prot_", names(.))))
}


#' Assign duplicated peptides to a leading protein
#' @param df A PSM data frame
#' @inheritParams mergePep
#' @inheritParams annotPSM
assign_duppeps <- function(df, group_psm_by = "pep_seq", 
                           group_pep_by = "pep_seq_mod", use_duppeps = TRUE, 
                           duppeps_repair = "denovo") 
{
  # Scenario: 
  # In `dat_file_1`, peptide_x assigned to Prn_MOUSE against "human + mouse" databases.
  # In `dat_file_2` the same peptide_x assigned to PRN_HUMAN against "human only" database.
  # When combining, `dat_file_1` and `dat_file_2`, all the peptide entries will be 
  #   re-assigned to the protein id with the greater `prot_n_pep`.
  
  message("Assigning multiple-dipped peptide sequences.\n")

  dat_dir <- get_gl_dat_dir()
  
  dup_peps <- df |>
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) |>
    dplyr::group_by(!!rlang::sym(group_psm_by)) |>
    dplyr::summarise(N = dplyr::n_distinct(!!rlang::sym(group_pep_by))) |>
    dplyr::filter(N > 1L)
  
  if (nrow(dup_peps)) {
    if (use_duppeps) {
      if (duppeps_repair == "denovo") {
        df <- local({
          # grps <- readRDS(file.path(dat_dir, "grps.rds"))
          grps <- mzion:::groupProts(unique(df[, c("prot_acc", "pep_seq")]), 
                                     out_path = dat_dir)
          sets <- readRDS(file.path(dat_dir, "prot_pep_setcover.rds"))
          ids <- with(sets, paste0(prot_acc, ".", pep_seq))
          
          # e.g. "prot_hit_num", "prot_family_member" may be not in df
          col_nms <- names(df)

          grps <- grps %>% 
            .[, names(.) %in% col_nms] %>% 
            tidyr::unite(prot_pep, prot_acc, pep_seq, sep = ".", remove = FALSE) %>% 
            dplyr::filter(prot_pep %in% ids) %>% 
            dplyr::select(-prot_pep)

          updated_nms <- names(grps) %>% .[! . == "pep_seq"]

          ans <- df %>% 
            dplyr::select(-which(names(.) %in% updated_nms)) %>% 
            dplyr::right_join(grps, by = "pep_seq") %>% 
            dplyr::select(col_nms)

          # keep the original pep_n_psm
          # update prot_n_psm, prot_n_pep and pep_n_psm later...
          # update other ^prot_ fields
          # and more...
          
          ans
        })
      }
      else if (duppeps_repair == "majority") {
        df <- local({
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
          
          rows <- df[[group_psm_by]] %in% dup_peps[[group_psm_by]]
          
          dups <- df[rows, ]
          unis <- df[!rows, ]
          
          ans <- lapply(split(dups, dups[[group_psm_by]]), 
                        replace_by_rowone, col_nms) %>% 
            dplyr::bind_rows() %>% 
            dplyr::select(names(df))
          
          dplyr::bind_rows(unis, ans)
        })
      }
      else {
        stop("Invalide choice of `duppeps_repair`." )
      }

      if (FALSE) {
        # update `dup_peps`; should be empty
        dup_peps_af <- df %>% 
          dplyr::filter(!!rlang::sym(group_psm_by) %in% dup_peps[[group_psm_by]]) %>%
          dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
          dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
          dplyr::summarise(N = n_distinct(!!rlang::sym(group_pep_by))) %>%
          dplyr::filter(N > 1)
        
        if (nrow(dup_peps_af)) {
          write.csv(dup_peps_af, file.path(dat_dir, "Peptide/dbl_dipping_peptides.csv"), 
                    row.names = FALSE)
          
          df <- df %>% 
            dplyr::filter(! (!!rlang::sym(group_psm_by) %in% dup_peps_af[[group_psm_by]]))
        }
      }
      
    } 
    else {
      write.csv(dup_peps, file.path(dat_dir, "Peptide/dbl_dipping_peptides.csv"), 
                row.names = FALSE)
      
      df <- df %>% 
        dplyr::filter(! (!!rlang::sym(group_psm_by) %in% dup_peps[[group_psm_by]]))
    }
  }
  else {
    # save the pep-prot map...
  }
  
  invisible(df)
}


#' Replaces values by the first row.
#'
#' @param df A data frame.
#' @param col_nms The column names under which the values in \code{df} will be
#'   replaced with those in the first row.
replace_by_rowone <- function (df, col_nms) 
{
  stopifnot(all(c("prot_n_pep", "prot_n_psm", "prot_mass") %in% names(df)))
  
  df <- df %>% 
    dplyr::arrange(-prot_n_pep, -prot_n_psm, -prot_mass)
  
  # if (group_pep_by == "prot_acc") use the first prot_acc and also the first gene
  # if (group_pep_by == "gene") use the first gene and also the first prot_acc
  
  cols_replace <- df[1, col_nms]
  
  df2 <- df[-1, ]
  df2[, col_nms] <- cols_replace
  
  dplyr::bind_rows(df[1, ], df2)
}


#' Imputes NA values in a peptide table
#' 
#' Not used.
#' 
#' @param df A intensity data frame (not matrix) of peptides.
#' @param fold The fold difference to mean.
impPepNA_byrow <- function (df, fold = 20)
{
  set.seed(1422)
  ncol <- ncol(df)
  
  for (i in seq_len(nrow(df))) {
    dfi <- df[i, ]
    nas <- is.na(dfi)
    nna <- sum(nas)
    
    if (nna == 0) {
      next
    }
    
    nok <- ncol - nna
    mv   <- sum(dfi[!nas]) / nok / fold
    errs <- 2^rnorm(nna, 0, .5)
    dfi[nas] <- mv * errs
    df[i, ] <- dfi
  }

  df
}


#' Imputes NA values in a peptide table
#'
#' @param df A intensity data frame (not matrix) of peptides.
#' @param fold The fold difference to mean.
#' @param is_intensity Logical; is intensity data or not (log2FC).
impPepNA <- function (df, fold = 50, is_intensity = TRUE)
{
  # if (!is.data.frame(df)) stop("Input `df` needs to be a data frame.")
  set.seed(1422)
  
  nas <- is.na(df)
  
  if (sum(nas) == 0) {
    return(df)
  }
  
  if (is_intensity) {
    mv   <- mean(df[!nas]) / fold
    errs <- 2^rnorm(sum(nas), 0, .5)
  }
  else {
    mv   <- mean(df[!nas])
    errs <- 2^rnorm(sum(nas), -2, .5)
  }

  nvs  <- mv * errs
  df[nas] <- nvs

  df
}

