#' Update data groups of `log2_R` etc. by logical matrix `lgl_log2_sd`
sd_lgl_cleanup <- function (df) {
  df[, grepl("^log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^N_log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^Z_log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^I[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^N_I[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^N_I[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^sd_log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df <- df %>% 
    dplyr::select(-grep("^lgl_log2_sd", names(.))) %>% 
    dplyr::filter(rowSums(!is.na(.[, grep("^log2_R[0-9]{3}", names(.) )])) > 0)
}


#' Filter data groups by a CV cut-off
#'
#' \code{purge_by_cv} replaces the data entries at \code{group CV > max_cv} to
#' NA.
#'
#' @inheritParams proteoHist
#' @param max_cv Numeric; the cut-off in maximum CV. Values above the threshold
#'   will be replaced with NA.
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
purge_by_cv <- function (df, id, max_cv, keep_ohw = TRUE) {
  if (!is.null(max_cv)) {
    stopifnot(is.numeric(max_cv))

    df_sd_lgl <- df %>% 
      dplyr::select(id, grep("^sd_log2_R[0-9]{3}[NC]*", names(.)))
    
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      # remove duplicated PSMs under the same id
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::filter(!duplicated(.[[id]]))
    } else if (id %in% c("prot_acc", "gene")) {
      # no duplicated id, but protein `sd` under each channel can be a mix of non-NA (one value) and NA
      # this is due to the `wide` joinining of data when forming `Peptide.txt`
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::group_by(!!rlang::sym(id)) %>% 
        dplyr::summarise_all(~ dplyr::first(na.omit(.x)))
    }
    
    if (keep_ohw) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, is.na(.x), -1E-3))
    }
    
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, .x > max_cv, NA)) %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, !is.na(.x), 1)) %>% 
      `names<-`(gsub("^sd_log2_R", "lgl_log2_sd", names(.)))

    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd_lgl, by = id) %>% 
      sd_lgl_cleanup()
  }
  
  return(df)
}


#' Filter data groups by quantiles of CV
#'
#' \code{purge_by_qt} replaces the data entries at \code{CV > quantile} with NA.
#'
#' @inheritParams proteoHist
#' @param pt_cv Numeric between 0 and 1; the percentile of CV. Values above the
#'   percentile threshold will be replaced with NA.
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
purge_by_qt <- function(df, id, pt_cv = NULL, keep_ohw = TRUE) {
  if (!is.null(pt_cv)) {
    stopifnot(is.numeric(pt_cv))
    if (pt_cv > 1) pt_cv <- pt_cv / 100
    
    df_sd_lgl <- df %>% 
      dplyr::select(id, grep("^sd_log2_R[0-9]{3}[NC]*", names(.)))
    
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::filter(!duplicated(.[[id]]))
    } else if (id %in% c("prot_acc", "gene")) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::group_by(!!rlang::sym(id)) %>% 
        dplyr::summarise_all(~ dplyr::first(na.omit(.x)))
    }
    
    qts <- df_sd_lgl %>% 
      .[, grepl("^sd_log2_R[0-9]{3}", names(.))] %>% 
      purrr::map_dbl(~ quantile(.x, probs = pt_cv, na.rm = TRUE))
    
    if (keep_ohw) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, is.na(.x), -1E-3))
    }
    
    df_sd_lgl[, grepl("^sd_log2_R[0-9]{3}[NC]*", names(df_sd_lgl))] <- 
      purrr::map2(as.list(df_sd_lgl[, grepl("^sd_log2_R[0-9]{3}[NC]*", names(df_sd_lgl))]), 
                  as.list(qts), ~ {
                    .x[.x > .y] <- NA
                    return(.x)
                  }) %>% 
      dplyr::bind_cols()
    
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, !is.na(.x), 1)) %>% 
      `names<-`(gsub("^sd_log2_R", "lgl_log2_sd", names(.)))
    
    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd_lgl, by = id) %>% 
      sd_lgl_cleanup()
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
#'TMT experiment per series of LC/MS. Note that greater CV may be encountered
#'for samples that are more different to reference material(s).
#'
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@inheritParams purge_by_qt
#'@param adjSD Logical; if TRUE, adjust the standard deviation in relative to
#'  the width of ratio profiles.
#'@param keep_ohw Logical; if TRUE, keep one-hit-wonders with unknown CV. The
#'  default is TRUE.
#'@param ... Additional parameters for plotting: \cr \code{ymax}, the maximum
#'  \eqn{y} at a log2 scale; the default is +0.6. \cr \code{ybreaks}, the breaks
#'  in \eqn{y}-axis at a log2 scale; the default is 0.2. \cr \code{width}, the
#'  width of plot. \cr \code{height}, the height of plot. \cr \code{flip_coord},
#'  logical; if TRUE, flip \code{x} and \code{y} axis.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/normPSM_examples.R
#'@example inst/extdata/examples/purgePSM_examples.R
#'@seealso \code{\link{purgePep}} for peptide data purging.
#'@export
purgePSM <- function (dat_dir = NULL, pt_cv = NULL, max_cv = NULL, adjSD = FALSE, keep_ohw = TRUE, ...) {
  
  dots <- rlang::enexprs(...)
  lang_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% lang_dots]
  
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
  # stopifnot(min_n > 0 & min_n%%1 == 0)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
    reorder_files(n_TMT_sets(label_scheme_full))
  
  dir.create(file.path(dat_dir, "PSM\\Copy"), recursive = TRUE, showWarnings = FALSE)

  for (fn in filelist) {
    file.copy(file.path(dat_dir, "PSM", fn), file.path(dat_dir, "PSM\\Copy", fn))
    
    df <- read.csv(file.path(dat_dir, "PSM", fn), check.names = FALSE, 
                   header = TRUE, sep = "\t", comment.char = "#") %>% 
      # filters_in_call(!!!lang_dots) %>% 
      purge_by_qt(group_psm_by, pt_cv, keep_ohw) %>% 
      purge_by_cv(group_psm_by, max_cv, keep_ohw) %>% 
      dplyr::filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0)
    
    cdns_changed <- any(!is.null(pt_cv), !is.null(max_cv))
    
    if (cdns_changed) {
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
    
    if (adjSD) {
      dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\purged_adj"), recursive = TRUE, showWarnings = FALSE)
      filepath <- file.path(dat_dir, "PSM\\log2FC_cv\\purged_adj", gsub("_PSM_N.txt", "_sd.png", fn))
    } else {
      filepath <- file.path(dat_dir, "PSM\\log2FC_cv\\purged", gsub("_PSM_N.txt", "_sd.png", fn))
    }    

    sd_violin(df, !!group_psm_by, filepath, type = "log2_R", adjSD = adjSD, is_psm = TRUE, ...)
  }
  
  mget(names(formals()), rlang::current_env()) %>% save_call("purPSM")
}


#'Purge peptide data
#'
#'\code{purgePep} removes \code{protein} entries from peptide table(s) by
#'selection criteria.
#'
#'The CV of proteins under each sample are first calculated from contributing
#'peptides. In the event of multiple sereis of LC/MS injections, the CV of the
#'same protein from different LC/MS will be summarised by median statistics.
#'
#'The data nullification will be applied column-wisely for all available
#'samples. Argument \code{col_select} is merely used to subsettng samples for
#'violin plot visualization.
#'
#'@inheritParams purgePSM
#'@inheritParams proteoHist
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@inheritParams purge_by_qt
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/normPSM_examples.R
#'@example inst/extdata/examples/normPep_examples.R
#'@example inst/extdata/examples/purgePep_examples.R
#'@seealso \code{\link{purgePSM}} for PSM data purging.
#'@export
purgePep <- function (dat_dir = NULL, pt_cv = NULL, max_cv = NULL, adjSD = FALSE, keep_ohw = TRUE, 
                      col_select = NULL, col_order = NULL, filename = NULL, ...) {
  
  dots <- rlang::enexprs(...)
  lang_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% lang_dots]

  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  filename <- rlang::enexpr(filename)
  if (is.null(filename)) {
    filename <- "Peptide_sd.png"
    fn_prefix <- "Peptide_sd"
    fn_suffix <- "png"
  } else {
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
  }
  
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
  # stopifnot(min_n > 0 & min_n%%1 == 0)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))

  dir.create(file.path(dat_dir, "Peptide\\Copy"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(dat_dir, "Peptide\\Peptide.txt"), file.path(dat_dir, "Peptide\\Copy\\Peptide.txt"))
  
  fn <- file.path(dat_dir, "Peptide", "Peptide.txt")
  df <- read.csv(fn, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    # filters_in_call(!!!lang_dots) %>% 
    purge_by_qt(group_pep_by, pt_cv, keep_ohw) %>% 
    purge_by_cv(group_pep_by, max_cv, keep_ohw) %>% 
    dplyr::filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0) 

  cdns_changed <- any(!is.null(pt_cv), !is.null(max_cv))
  
  if (cdns_changed) {
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
  
  width <- eval(dots$width, env = caller_env())
  height <- eval(dots$height, env = caller_env())
  
  if (is.null(width)) width <- 8 * n_TMT_sets(label_scheme)
  if (is.null(height)) height <- 8
  dots <- dots %>% .[! names(.) %in% c("width", "height")]
  
  if (adjSD) {
    dir.create(file.path(dat_dir, "Peptide\\log2FC_cv\\purged_adj"), recursive = TRUE, showWarnings = FALSE)
    filepath <- file.path(dat_dir, "Peptide\\log2FC_cv\\purged_adj", filename)
  } else {
    filepath <- file.path(dat_dir, "Peptide\\log2FC_cv\\purged", filename)
  }
  
  sd_violin(df, !!group_pep_by, filepath, width = width, height = height, 
            type = "log2_R", adjSD = adjSD, is_psm = FALSE, 
            col_select = !!col_select, col_order = !!col_order, !!!dots)

  mget(names(formals()), rlang::current_env()) %>% save_call("purPSM")
}





