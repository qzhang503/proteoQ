#' Violin plots of SDs
#'
#' \code{sd_violin_full} visualizes the SD distribution of SD
#'
#' @import dplyr purrr rlang ggplot2
#' @importFrom magrittr %>%
sd_violin_full <- function(df_sd, id, label_scheme, filepath, filename) {
  id <- rlang::as_string(rlang::enexpr(id))
  
  Levels <- names(df_sd) %>% 
    .[grepl("^log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
    gsub("^log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", .)  
  
  df_sd <- df_sd %>%
    `names<-`(gsub("^log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = !!rlang::sym(id), value = "SD") %>%
    dplyr::rename(Channel := !!rlang::sym(id)) %>% 
    dplyr::ungroup(Channel) %>% 
    dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
    dplyr::filter(!is.na(SD))
  
  p <- ggplot() +
    geom_violin(df_sd, mapping = aes(x = Channel, y = SD, fill = Channel), size = .25) +
    geom_boxplot(df_sd, mapping = aes(x = Channel, y = SD), width = 0.1, lwd = .2, fill = "white") +
    stat_summary(df_sd, mapping = aes(x = Channel, y = SD), fun.y = "mean", geom = "point",
                 shape=23, size=2, fill="white", alpha=.5) +
    labs(title = expression(""), x = expression("Channel"), y = expression("SD ("*log[2]*"FC)")) +
    scale_y_continuous(limits = c(0, .6), breaks = seq(0, .6, .2)) +
    theme_psm_violin
  
  try(ggsave(file.path(filepath, filename), p, width = 7* n_TMT_sets(label_scheme), height = 7, units = "in"))
}


#'Rmove SD outliers
#'
#'\code{purgeData} removes entries with SD > cv_cutoff
#'
#'@import dplyr purrr rlang
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
purgeData <- function(df, id, cv_cutoff = NULL, nseq_cutoff = 1, ...) {
  id <- rlang::as_string(rlang::enexpr(id))
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  stopifnot(nseq_cutoff > 0 & nseq_cutoff%%1 == 0)
  
  df <- df %>% dplyr::arrange(!!rlang::sym(id))
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
      reorder_files(n_TMT_sets(label_scheme_full))
    
    if (is_empty(filelist)) {
      stop("PSM files not found; may be you are starting with MaxQuant peptide outputs.", call. = FALSE)
    }
    
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
      dplyr::arrange(TMT_Set) %>% 
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
    
    df_sd <- df_sd %>% 
      dplyr::select(!!rlang::sym(id), grep("log2_R[0-9]{3}[NC]*", names(.))) %>% 
      dplyr::arrange(!!rlang::sym(id)) %T>%
      write.csv(file.path(dat_dir, "Peptide\\cache", "pep_sd.csv"), row.names = FALSE)
  } else if (id %in% c("prot_acc", "gene")) {
    df_sd <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"), check.names = FALSE, 
                      header = TRUE, sep = "\t", comment.char = "#") %>% 
      filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0) %>% 
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
      dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
      dplyr::arrange(!!rlang::sym(id))
  }
  
  # purging
  if (!is.null(cv_cutoff)) {
    stopifnot(is.numeric(cv_cutoff))
    
    df_sd <- df_sd %>% 
      dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, .x > cv_cutoff, NA)) 
    
    df_sd_lgl <- df_sd%>% 
      dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, !is.na(.x), 1)) %>% 
      dplyr::filter(!!rlang::sym(id) %in% df[[id]]) %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      tibble::column_to_rownames(id) %>% 
      `names<-`(paste0("sd_", names(.))) %>% 
      tibble::rownames_to_column(id)
    
    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd_lgl, by = id)
    
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
    filepath <- file.path(dat_dir, "Peptide\\Purge")
    filename <- "Peptide_SD.png"
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
    df <- df %>% dplyr::filter(n_psm >= nseq_cutoff)
  } else if (id %in% c("prot_acc", "gene")) {
    filepath <- file.path(dat_dir, "Protein\\Purge")
    filename <- "Protein_SD.png"
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
    df <- df %>% dplyr::filter(n_pep >= nseq_cutoff)
  }
  
  df <- df %>% 
    dplyr::filter(rowSums(!is.na(.[, grep("^log2_R[0-9]{3}", names(.) )])) > 0)
  
  # violin plots
  df_sd <- df_sd %>% 
    dplyr::filter(!!rlang::sym(id) %in% df[[id]])
  
  sd_violin_full(df_sd, !!id, label_scheme, filepath, filename)

  return(df)
}


#'Purge data
#'
#'\code{proteoPurge} provides additional cleanup of protein or peptide data.
#'
#'The function matches the current \code{id} to those in the latest \code{call}
#'to \code{\link{normPep}} or \code{\link{normPrn}}.  For example, if
#'\code{pep_seq} was used in \code{\link{normPep}}, the current \code{id =
#'pep_seq_mod} will be matched to \code{id = pep_seq}.
#'
#'@inheritParams proteoHist
#'@param cv_cutoff Numeric; the cut-off of CV. The CVs are from the ascribing
#'  PSMs for peptide data or ascribing peptides for protein data.
#'@param nseq_cutoff Positive integer. When calling from \code{purPep}, peptide
#'  entries in \code{Peptide.txt} with the number of identifying PSMs smaller
#'  than \code{nseq_cutoff} will be replaced with NA. When calling from
#'  \code{purPrn}, protein entries in \code{Protein.txt} with the number of
#'  identifying peptides smaller than \code{nseq_cutoff} will be replaced with
#'  NA.
#'@param ... Additional parameters for plotting yet to be defined.
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
  
  info_anal(id = !!id, 
            df = !!df, filepath = !!filepath, filename = !!filename,
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
#'@import purrr
#'@export
purPep <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\Purge\\log"), recursive = TRUE, showWarnings = FALSE)
  quietly_log <- purrr::quietly(proteoPurge)(id = pep_seq, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\Purge\\log","pepPurge_log.csv"), append = TRUE)
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
#'@import purrr
#'@export
purPrn <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Purge\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoPurge)(id = gene, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Purge\\log","prnPurge_log.csv"), append = TRUE)
}






#' Filter data groups by a CV cut-off
#'
#' \code{purge_by_cv} replaces the data entries at \code{group CV > max_cv} to
#' NA.
#'
#' @inheritParams proteoHist
#' @param max_cv Numeric; the cut-off in CV. The CV are from the ascribing PSMs
#'   for the filtration of peptide entries or ascribing peptides for the
#'   filtration of protein entries.
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
purge_by_cv <- function (df, id, max_cv) {
  if (!is.null(max_cv)) {
    stopifnot(is.numeric(max_cv))
    
    df_sd_lgl <- df %>% 
      dplyr::select(id, grep("^sd_log2_R[0-9]{3}[NC]*", names(.))) %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
                       ~ replace(.x, .x > max_cv, NA)) %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
                       ~ replace(.x, !is.na(.x), 1)) %>% 
      `names<-`(gsub("^sd_log2_R", "lgl_log2_sd", names(.))) %>% 
      dplyr::filter(!duplicated(.[[id]]))
    
    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd_lgl, by = id)
    
    df[, grepl("^log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^Z_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^sd_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df <- df %>% 
      dplyr::select(-grep("^lgl_log2_sd", names(.))) %>% 
      dplyr::filter(rowSums(!is.na(.[, grep("^log2_R[0-9]{3}", names(.) )])) > 0)
  }
  
  return(df)
}


#' Filter data groups by a minimal number of observations (n_obs)
#'
#' \code{purge_by_n} replaces the data entries at \code{group n_obs < min_n} to
#' NA.
#'
#' @inheritParams proteoHist
#' @param min_n Positive integer. When calling from \code{purgePSM}, peptide
#'   entries in PSM tables with the number of identifying PSMs smaller than
#'   \code{min_n} will be replaced with NA. When calling from \code{purgePep},
#'   protein entries in peptide tables with the number of identifying peptides
#'   smaller than \code{min_n} will be replaced with NA.
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
purge_by_n <- function (df, id, min_n) {
  kept <- df %>%
    dplyr::select(!!rlang::sym(id)) %>%
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= min_n) %>% 
    dplyr::select(id) %>% 
    unlist() %>% 
    as.character()

  df %>% 
    dplyr::filter(.[[id]] %in% kept) 
}


#'Purge PSM data
#'
#'\code{purgePSM} removes \code{peptide} entries from PSM tables by selection
#'criteria.
#'
#'The CV of peptides are calculated from contributing PSMs at the basis of per
#'TMT experiment per LCMS run.
#'
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@param ... Additional parameters for plotting: \cr \code{ymax}, the maximum
#'  \eqn{y} at a log2 scale; the default is +0.6. \cr \code{y_breaks}, the
#'  breaks in \eqn{y}-axis at a log2 scale; the default is 0.2. \cr
#'  \code{width}, the width of plot. \cr \code{height}, the height of plot.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@export
purgePSM <- function (dat_dir = NULL, max_cv = NULL, min_n = 1, ...) {

  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  
  group_psm_by <- match_normPSM_pepid()
  group_pep_by <- match_normPSM_protid()
  
  stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  stopifnot(min_n > 0 & min_n%%1 == 0)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
    reorder_files(n_TMT_sets(label_scheme_full))
  
  dir.create(file.path(dat_dir, "PSM\\Copy"), recursive = TRUE, showWarnings = FALSE)

  for (fn in filelist) {
    file.copy(file.path(dat_dir, "PSM", fn), file.path(dat_dir, "PSM\\Copy", fn))
    
    df <- read.csv(file.path(dat_dir, "PSM", fn), check.names = FALSE, 
                   header = TRUE, sep = "\t", comment.char = "#")
    
    df <- df %>% 
      purge_by_cv(group_psm_by, max_cv) %>% 
      purge_by_n(group_psm_by, min_n) %>% 
      dplyr::filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0)
    
    if ((!is.null(max_cv)) | (min_n > 1)) {
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
        dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_pep = n())
      
      df <- df %>% 
        dplyr::left_join(pep_n_psm, by = group_psm_by) %>% 
        dplyr::mutate(pep_n_psm.x = pep_n_psm.y) %>% 
        dplyr::rename("pep_n_psm" = "pep_n_psm.x") %>% 
        dplyr::select(-pep_n_psm.y)
  
      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by) %>% 
        dplyr::mutate(prot_n_psm.x = prot_n_psm.y, prot_n_pep.x = prot_n_pep.y) %>% 
        dplyr::rename(prot_n_psm = prot_n_psm.x, prot_n_pep = prot_n_pep.x) %>% 
        dplyr::select(-prot_n_psm.y, -prot_n_pep.y) %T>% 
        write.table(file.path(dat_dir, "PSM", fn), sep = "\t", col.names = TRUE, row.names = FALSE)      
    }

    df %>% 
      dplyr::select(group_psm_by, grep("^sd_log2_R[0-9]{3}", names(.))) %>% 
      dplyr::filter(!duplicated(.[[group_psm_by]])) %>% 
      dplyr::filter(rowSums(!is.na(.[grep("^sd_log2_R[0-9]{3}", names(.))])) > 0) %>%
      sd_violin(!!group_psm_by, 
                file.path(dat_dir, "PSM\\log2FC_cv\\purged", gsub("_PSM_N.txt", "_sd.png", fn)), ...)
  }
  
  mget(names(formals()), rlang::current_env()) %>% save_call("purPSM")
}


#'Purge peptide data
#'
#'\code{purgePep} removes \code{protein} entries from peptide table(s) by
#'selection criteria.
#'
#'The CV of proteins are calculated from contributing peptides at the basis of
#'per TMT experiment per LCMS run.
#'
#'@inheritParams purgePSM
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@export
purgePep <- function (dat_dir = NULL, max_cv = NULL, min_n = 1, ...) {
  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  
  group_psm_by <- match_normPSM_pepid()
  group_pep_by <- match_normPSM_protid()

  stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  stopifnot(min_n > 0 & min_n%%1 == 0)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  dir.create(file.path(dat_dir, "Peptide\\Copy"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(dat_dir, "Peptide\\Peptide.txt"), file.path(dat_dir, "Peptide\\Copy\\Peptide.txt"))
  
  fn <- file.path(dat_dir, "Peptide", "Peptide.txt")
  df <- read.csv(fn, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    purge_by_cv(group_pep_by, max_cv) %>% 
    purge_by_n(group_pep_by, min_n) %>% 
    dplyr::filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0) 
  
  if ((!is.null(max_cv)) | (min_n > 1)) {
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
  
    df <- df %>% 
      dplyr::left_join(pep_n_psm, by = group_psm_by) %>% 
      dplyr::mutate(pep_n_psm.x = pep_n_psm.y) %>% 
      dplyr::rename("pep_n_psm" = "pep_n_psm.x") %>% 
      dplyr::select(-pep_n_psm.y)
    
    df <- list(df, prot_n_psm, prot_n_pep) %>%
      purrr::reduce(left_join, by = group_pep_by) %>% 
      dplyr::mutate(prot_n_psm.x = prot_n_psm.y, prot_n_pep.x = prot_n_pep.y) %>% 
      dplyr::rename(prot_n_psm = prot_n_psm.x, prot_n_pep = prot_n_pep.x) %>% 
      dplyr::select(-prot_n_psm.y, -prot_n_pep.y) %T>% 
      write.table(fn, sep = "\t", col.names = TRUE, row.names = FALSE)    
  }
  
  dots <- rlang::enexprs(...)
  width <- eval(dots$width, env = caller_env())
  height <- eval(dots$height, env = caller_env())
  
  if (is.null(width)) width <- 8 * n_TMT_sets(label_scheme)
  if (is.null(height)) height <- 8
  
  dots <- dots %>% .[! names(.) %in% c("width", "height")]

  df %>% 
    dplyr::select(group_pep_by, grep("^sd_log2_R[0-9]{3}", names(.))) %>% 
    dplyr::filter(!duplicated(.[[group_pep_by]])) %>% 
    dplyr::filter(rowSums(!is.na(.[grep("^sd_log2_R[0-9]{3}", names(.))])) > 0) %>% 
    sd_violin(!!group_pep_by, file.path(dat_dir, "Peptide\\log2FC_cv\\purged", "Peptide_sd.png"), 
              width = width, height = height, type = "log2_R", !!!dots)

  mget(names(formals()), rlang::current_env()) %>% save_call("purPSM")
}
