#'Rmove SD outliers
#'
#'\code{purgeData} removes entries with SD > cv_cutoff
#'
#'@import dplyr purrr rlang
#'@importFrom magrittr %>%
purgeData <- function(df, id, label_scheme, cv_cutoff = NULL, nseq_cutoff = 1, ...) {
  id <- rlang::as_string(rlang::enexpr(id))
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  stopifnot(nseq_cutoff > 0 & nseq_cutoff%%1 == 0)
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
      reorder_files(n_TMT_sets(label_scheme_full))
    
    df_sd <- purrr::map(as.list(filelist), ~ {
      df <- read.csv(file.path(dat_dir, "PSM", .x), check.names = FALSE, header = TRUE,
                     sep = "\t", comment.char = "#") %>% 
        dplyr::mutate(TMT_Set = as.integer(gsub("^TMTset(\\d+)_LCMSinj.*", "\\1", .x)))
      
      return(df)
    }) %>% do.call("rbind", .)
    
    df_sd <- df_sd %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::select(!!rlang::sym(id), TMT_Set, grep("^log2_R[0-9]{3}", names(.))) %>%
      dplyr::group_by(!!rlang::sym(id), TMT_Set) %>%
      dplyr::summarise_at(vars(starts_with("log2_R")), ~ sd(.x, na.rm = TRUE)) %>% 
      dplyr::arrange(TMT_Set) %T>%
      write.csv(file.path(dat_dir, "Peptide\\cache", "pep_sd.csv"), row.names = FALSE) %>% 
      tidyr::gather(grep("log2_R[0-9]{3}", names(.)), key = ID, value = value) %>%
      tidyr::unite(ID, ID, TMT_Set)
    
    Levels <- unique(df_sd$ID)
    df_sd <- df_sd %>%
      dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
      tidyr::spread(ID, value) %>% 
      dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3))
    
    for (set_idx in seq_len(n_TMT_sets(label_scheme))) {
      df_sd <- newColnames(set_idx, df_sd, label_scheme)
    }
  } else if (id %in% c("prot_acc", "gene")) {
    df_sd <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"), check.names = FALSE, 
                      header = TRUE, sep = "\t", comment.char = "#") %>% 
      filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0)
    
    df_sd <- df_sd %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::select(!!rlang::sym(id), grep("^log2_R[0-9]{3}", names(.))) %>%
      dplyr::group_by(!!rlang::sym(id)) %>%
      dplyr::summarise_at(vars(starts_with("log2_R")), ~ sd(.x, na.rm = TRUE)) %T>% 
      write.csv(file.path(dat_dir, "Protein\\cache", "prn_sd.csv"), row.names = FALSE) %>% 
      tidyr::gather(grep("log2_R[0-9]{3}", names(.)), key = ID, value = value)
    
    Levels <- unique(df_sd$ID)
    df_sd <- df_sd %>%
      dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
      tidyr::spread(ID, value) %>% 
      dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3))
  }
  
  # violin plots
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    filepath <- file.path(dat_dir, "Peptide\\Purge")
    filename <- "Peptide_SD.png"
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  } else if (id %in% c("prot_acc", "gene")) {
    filepath <- file.path(dat_dir, "Protein\\Purge")
    filename <- "Protein_SD.png"
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  }
  sd_violin_full(df_sd, !!id, label_scheme, filepath, filename)

  # purging
  if (!is.null(cv_cutoff)) {
    stopifnot(is.numeric(cv_cutoff))

    df_sd <- df_sd %>% 
      dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, is.na(.x), 0)) %>% 
      dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, .x > cv_cutoff, NA)) %>% 
      dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, !is.na(.x), 1)) %>% 
      dplyr::filter(!!rlang::sym(id) %in% df[[id]]) %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      tibble::column_to_rownames(id) %>% 
      `names<-`(paste0("sd_", names(.))) %>% 
      tibble::rownames_to_column(id)
    
    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd, by = id)
  
    df[, grepl("^log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^Z_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df <- df %>% 
      dplyr::select(-grep("^sd_log2_R", names(.)))
  }

  if (id %in% c("pep_seq", "pep_seq_mod")) {
    df <- df %>% 
      dplyr::filter(n_psm >= nseq_cutoff)
  } else if (id %in% c("prot_acc", "gene")) {
    df <- df %>% 
      dplyr::filter(n_pep >= nseq_cutoff)
  }
  
  df <- df %>% 
    dplyr::filter(rowSums(!is.na(.[, grep("^log2_R[0-9]{3}", names(.) )])) > 0)

  return(df)
}


#'Purge ata
#'
#'\code{proteoPurge} provides additional cleanup of protein or peptide data.
#'
#'The function matches the current \code{id} to those in the latest \code{call}
#'to \code{\link{normPep}} or \code{\link{normPrn}}.  For example, if
#'\code{pep_seq} was used in \code{\link{normPep()}}, the current \code{id =
#'pep_seq_mod} will be matched to \code{id = pep_seq}.
#'
#'@inheritParams proteoHist
#'@param cv_cutoff Numeric; the cut-off of CV. The CVs are from the ascribing
#'  PSMs for peptide data or ascribing peptides for protein data.
#'@param nseq_cutoff Positive integer; for peptide data, peptide entries with the
#'  number of identifying PSMs smaller than the cutoff will be removed; for
#'  protein data, protein entries with the number of identifying peptides
#'  smaller than the cutoff will be removed.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoPurge <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                         cv_cutoff = NULL, nseq_cutoff = 1, 
                         df = NULL, filepath = NULL, filename = NULL, ...) {
  
  id <- rlang::enexpr(id)
  if(length(id) != 1) id <- rlang::expr(gene)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "Purge")(cv_cutoff = cv_cutoff, nseq_cutoff = nseq_cutoff, ...)
}


#'Purge peptide data
#'
#'\code{purPep} is a wrapper of \code{\link{proteoPurge}} for peptide data
#'
#'@rdname proteoPurge
#'
#' @examples
#' \dontrun{
#' purPep(
#'   cv_cutoff = 0.4,
#'   nseq_cutoff = 1,
#' )
#' }
#'
#'@export
purPep <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  proteoPurge(id = pep_seq, ...)
}


#'Purge protein data
#'
#'\code{purPrn} is a wrapper of \code{\link{proteoPurge}} for protein data
#'
#'@rdname proteoPurge
#'
#' @examples
#' \dontrun{
#' purPrn(
#'   cv_cutoff = 0.4,
#'   nseq_cutoff = 1,
#' )
#' }
#'
#'@export
purPrn <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  proteoPurge(id = gene, ...)
}


