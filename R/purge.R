#' Updates data groups of `log2_R_` and `I_` etc. by logical matrix
#' `lgl_log2_...`.
#'
#' @inheritParams info_anal
#' @param type Type of logical matrix.
lgl_cleanup <- function (df, type = "sd", rm_allna = FALSE) 
{
  pat <- "[0-9]{3}[NC]{0,1}"
  fields <- paste0("^lgl_log2_", type, pat)
  pat_log2_R <- paste0("^log2_R", pat)
  pat_N_log2_R <- paste0("^N_log2_R", pat)
  pat_Z_log2_R <- paste0("^Z_log2_R", pat)
  pat_I <- paste0("^I", pat)
  pat_N_I <- paste0("^N_I", pat)
  pat_sd <- paste0("^sd_log2_R", pat)

  nms <- names(df)
  
  df[, grepl(pat_log2_R, nms)] <-
    purrr::map2(as.list(df[, grepl(pat_log2_R, nms)]),
                as.list(df[, grepl(fields, nms)]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl(pat_N_log2_R, nms)] <-
    purrr::map2(as.list(df[, grepl(pat_N_log2_R, nms)]),
                as.list(df[, grepl(fields, nms)]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl(pat_Z_log2_R, nms)] <-
    purrr::map2(as.list(df[, grepl(pat_Z_log2_R, nms)]),
                as.list(df[, grepl(fields, nms)]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl(pat_I, nms)] <-
    purrr::map2(as.list(df[, grepl(pat_I, nms)]),
                as.list(df[, grepl(fields, nms)]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl(pat_N_I, nms)] <-
    purrr::map2(as.list(df[, grepl(pat_N_I, nms)]),
                as.list(df[, grepl(fields, nms)]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl(pat_sd, nms)] <-
    purrr::map2(as.list(df[, grepl(pat_sd, nms)]),
                as.list(df[, grepl(fields, nms)]), `*`) %>%
    dplyr::bind_cols()
  
  df <- df %>% 
    dplyr::select(-grep(fields, names(.))) %>% 
    { if (rm_allna) .[rowSums(!is.na(.[grepl(pat_log2_R, names(.))])) > 0L, ] 
      else .} 
}


#' Filters data groups by a CV cut-off
#'
#' \code{purge_by_cv} replaces the data entries at \code{group CV > max_cv} to
#' NA.
#'
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @inheritParams purgePSM
#' @param max_cv Numeric; the cut-off in maximum CV. Values above the threshold
#'   will be replaced with NA. The default is NULL with no data trimming by max
#'   CV.
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
purge_by_cv <- function (df, id, max_cv = NULL, min_n = 1L, rm_allna = FALSE) 
{
  if (is.null(max_cv)) 
    return(df)

  if (!is.numeric(max_cv)) 
    stop("Need numeric value for `max_cv`.")

  if (max_cv > 1) 
    max_cv <- max_cv / 100
  
  pat <- "^sd_log2_R[0-9]{3}[NC]{0,1}"
  
  df_sd_lgl <- df %>% 
    dplyr::select(id, grep(pat, names(.)))

  # remove duplicated PSMs under the same id
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::filter(!duplicated(.[[id]]))
  } 
  # no duplicated id, but protein `sd` under each channel can be 
  #   a mix of non-NA (one value) and NA
  # this is due to the `wide` joining of data when forming `Peptide.txt`
  else if (id %in% c("prot_acc", "gene")) {
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::group_by(!!rlang::sym(id)) %>% 
      dplyr::summarise_all(~ dplyr::first(na.omit(.x)))
  }
  
  df_sd_lgl <- df_sd_lgl %>% 
    dplyr::mutate_at(vars(grep(pat, names(.))), 
                     ~ replace(.x, .x > max_cv, NA_real_)) %>% 
    dplyr::mutate_at(vars(grep(pat, names(.))), 
                     ~ replace(.x, !is.na(.x), 1)) %>% 
    `names<-`(gsub("^sd_log2_R", "lgl_log2_sd", names(.)))
  
  df <- df %>% 
    dplyr::arrange(!!rlang::sym(id)) %>% 
    dplyr::left_join(df_sd_lgl, by = id) %>% 
    lgl_cleanup(rm_allna = rm_allna)
}


#' Filters data groups by quantiles of CV
#'
#' \code{purge_by_qt} replaces the data entries at \code{CV > quantile} with NA.
#'
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @inheritParams purgePSM
#' @param pt_cv Numeric between 0 and 1; the percentile of CV. Values above the
#'   percentile threshold will be replaced with NA. The default is NULL with no
#'   data trimming by CV percentile. The precedence in data purging is
#'   \code{pt_cv} \eqn{\ge} \code{max_cv} \eqn{\ge} \code{min_n}.
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
purge_by_qt <- function(df, id, pt_cv = NULL, min_n = 1L, rm_allna = FALSE) 
{
  if (is.null(pt_cv)) 
    return(df)
  
  if (!is.numeric(pt_cv))
    stop("Need numeric value for `pt_cv`.")
  
  if (pt_cv > 1) 
    pt_cv <- pt_cv / 100
  
  pat <- "^sd_log2_R[0-9]{3}[NC]{0,1}"
  
  df_sd_lgl <- df %>% 
    dplyr::select(id, grep(pat, names(.)))
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::filter(!duplicated(.[[id]]))
  } 
  else if (id %in% c("prot_acc", "gene")) {
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::group_by(!!rlang::sym(id)) %>% 
      dplyr::summarise_all(~ dplyr::first(na.omit(.x)))
  }
  
  qts <- df_sd_lgl %>% 
    .[, grepl(pat, names(.))] %>% 
    purrr::map_dbl(~ quantile(.x, probs = pt_cv, na.rm = TRUE))
  
  df_sd_lgl[, grepl(pat, names(df_sd_lgl))] <- 
    purrr::map2(as.list(df_sd_lgl[, grepl(pat, names(df_sd_lgl))]), 
                as.list(qts), ~ {
                  .x[.x > .y] <- NA_real_
                  return(.x)
                }) %>% 
    dplyr::bind_cols()
  
  df_sd_lgl <- df_sd_lgl %>% 
    dplyr::mutate_at(vars(grep(pat, names(.))), 
                     ~ replace(.x, !is.na(.x), 1)) %>% 
    `names<-`(gsub("^sd_log2_R", "lgl_log2_sd", names(.)))
  
  df <- df %>% 
    dplyr::arrange(!!rlang::sym(id)) %>% 
    dplyr::left_join(df_sd_lgl, by = id) %>% 
    lgl_cleanup(rm_allna = rm_allna)
}


#' Filters data groups by a minimal number of observations (n_obs)
#'
#' \code{purge_by_n} replaces the data entries at \code{group n_obs < min_n} to
#' NA.
#'
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @param min_n Positive integer. When calling from \code{purgePSM}, peptide
#'   entries in PSM tables with the number of identifying PSMs smaller than
#'   \code{min_n} will be replaced with NA. When calling from \code{purgePep},
#'   protein entries in peptide tables with the number of identifying peptides
#'   smaller than \code{min_n} will be replaced with NA. 
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
purge_by_n <- function (df, id, min_n = 1L, rm_allna = FALSE) 
{
  if (!is.numeric(min_n))
    stop("Need numeric value for `min_n`.")

  if (min_n == 1L) 
    return(df)
  
  pat <- "[0-9]{3}[NC]{0,1}"
  pat_log2_R <- paste0("^log2_R", pat)
  pat_lgl <- paste0("^lgl_log2_R", pat)

  df_lgl <- df %>% 
    dplyr::select(id, grep(pat_log2_R, names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ sum(!is.na(.x))) %>% 
    `names<-`(gsub("^log2_R", "lgl_log2_R", names(.)))
  
  nms <- names(df_lgl)
  
  df_lgl[, grepl(pat_lgl, nms)] <- 
    purrr::map(as.list(df_lgl[, grepl(pat_lgl, nms)]), 
               ~ {
                 .x[.x < min_n] <- NA_real_
                 return(.x)
               }, min_n) %>% 
    dplyr::bind_cols() %>% 
    dplyr::mutate_at(vars(grep(pat_lgl, names(.))), 
                     ~ replace(.x, !is.na(.x), 1))
  
  df <- df %>% 
    dplyr::arrange(!!rlang::sym(id)) %>% 
    dplyr::left_join(df_lgl, by = id) %>% 
    lgl_cleanup(type = "R", rm_allna = rm_allna)
}


#' Helper of \link{purgePSM}
#' 
#' @param file A list of file (names).
#' @inheritParams purgePSM
#' @inheritParams annotPSM
#' @inheritParams prnHist
psm_mpurge <- function (file, dat_dir = NULL, group_psm_by = "pep_seq", 
                        group_pep_by = "prot_acc", pt_cv = NULL, max_cv = NULL, 
                        min_n = 1L, rm_allna = FALSE, theme, ...) 
{
  dots <- rlang::enexprs(...)

  file.copy(file.path(dat_dir, "PSM", file), file.path(dat_dir, "PSM/Copy", file))

  df <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, "PSM", file), col_types = get_col_types(), 
                    show_col_types = FALSE)
  )
  
  pat_log2_R <- "^log2_R[0-9]{3}[NC]{0,1}"

  df <- df %>% 
    purge_by_qt(id = group_psm_by, pt_cv = pt_cv, 
                min_n = min_n, rm_allna = rm_allna) %>% 
    purge_by_cv(id = group_psm_by, max_cv = max_cv, 
                min_n = min_n, rm_allna = rm_allna) %>% 
    purge_by_n(id = group_psm_by, min_n = min_n, rm_allna = rm_allna) %>% 
    { if (rm_allna) .[rowSums(!is.na(.[grepl(pat_log2_R, names(.))])) > 0L, ] 
      else . }

  # update
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
    purrr::reduce(dplyr::left_join, by = group_pep_by) %>% 
    dplyr::mutate(prot_n_psm.x = prot_n_psm.y, prot_n_pep.x = prot_n_pep.y) %>% 
    dplyr::rename(prot_n_psm = prot_n_psm.x, prot_n_pep = prot_n_pep.x) %>% 
    dplyr::select(-prot_n_psm.y, -prot_n_pep.y) %T>% 
    readr::write_tsv(file.path(dat_dir, "PSM", file))

  # plot
  width <- eval(dots$width, envir = rlang::caller_env())
  height <- eval(dots$height, envir = rlang::caller_env())
  
  if (is.null(width)) width <- 8 
  if (is.null(height)) height <- 8
  
  dots <- dots %>% .[! names(.) %in% c("width", "height")]
  
  filepath <- file.path(dat_dir, "PSM/log2FC_cv/purged", 
                        gsub("_PSM_N.txt", "_sd.png", file))
  
  quiet_out <- purrr::quietly(sd_violin)(df = df, 
                                         id = !!group_psm_by, 
                                         filepath = filepath, 
                                         width = width, 
                                         height = height, 
                                         type = "log2_R", 
                                         adjSD = FALSE, 
                                         is_psm = TRUE, 
                                         col_select = NULL, 
                                         col_order = NULL, 
                                         theme = theme, 
                                         !!!dots)
}


#'Purges PSM data
#'
#'\code{purgePSM} removes \code{peptide} entries from PSM tables by quality
#'criteria. It further plots the distributions of the CV of \code{log2FC} by TMT
#'experiments and LC/MS series. The utility will have no effect against LFQ data
#'using MS1 peak area/intensity.
#'
#'The CV of peptides are calculated from contributing PSMs at the basis of per
#'TMT experiment per series of LC/MS. Note that greater CV may be encountered
#'for samples that are more different to reference material(s).
#'
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@inheritParams purge_by_qt
#'@inheritParams plot_prnTrend
#'@param adjSD Not currently used. If TRUE, adjust the standard deviation in
#'  relative to the width of ratio profiles.
#'@param ... Additional parameters for plotting: \cr \code{ymax}, the maximum
#'  \eqn{y} at a log2 scale. \cr \code{ybreaks}, the breaks in \eqn{y}-axis at a
#'  log2 scale. \cr \code{width}, the width of plot. \cr \code{height}, the
#'  height of plot. \cr \code{flip_coord}, logical; if TRUE, flip \code{x} and
#'  \code{y} axis.
#'@import dplyr ggplot2
#'@importFrom rlang exprs expr
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@example inst/extdata/examples/purgePSM_.R
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
#'@export
purgePSM <- function (dat_dir = NULL, pt_cv = NULL, max_cv = NULL, adjSD = FALSE, 
                      min_n = 1L, rm_allna = FALSE, theme= NULL, ...) 
{
  on.exit(
    mget(names(formals()), envir = environment(), inherits = FALSE) %>% 
      c(rlang::enexprs(...)) %>% 
      save_call(fun), add = TRUE)
  
  this_call <- match.call()
  fun <- as.character(this_call[[1]])
  
  check_dots(c("id", "anal_type", "filename", "filepath"), ...)
  
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dots <- dots[! dots %in% filter_dots]
  
  if (length(filter_dots)) {
    warning("`filter_` not applicable with `", fun, "`")
  }
  
  if (is.null(dat_dir)) {
    dat_dir <- get_gl_dat_dir()
  }
  else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  
  dir.create(file.path(dat_dir, "PSM/log2FC_cv/purged"), 
             recursive = TRUE, showWarnings = FALSE)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(
    group_psm_by %in% c("pep_seq", "pep_seq_mod"), 
    group_pep_by %in% c("prot_acc", "gene"), 
    length(group_psm_by) == 1L, 
    length(group_pep_by) == 1L
  )
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))

  if (!(tmt_plex <- TMT_plex(label_scheme_full))) {
    stop("No PSM purging for LFQ using MS1 peak area.", call. = FALSE)
  }
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), 
                         pattern = "_PSM_N\\.txt$", full.names = TRUE) |>
    reorder_files2()
  basenames  <- attr(filelist, "basenames")
  set_idxes  <- attr(filelist, "set_idxes")
  injn_idxes <- attr(filelist, "injn_idxes")

  if (!(n_files <- length(filelist))) {
    stop("Files of \"_PSM_N.txt\" not found.")
  }
  
  dir.create(
    file.path(dat_dir, "PSM/Copy"), recursive = TRUE, showWarnings = FALSE)
  
  if (n_files > 3L) {
    n_cores <- min(parallel::detectCores(), n_files)
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(cl, list("purge_by_qt", "purge_by_cv", "purge_by_n"), 
                            envir = environment(proteoQ:::purge_by_qt))
    parallel::clusterExport(cl, "dat_dir")
    
    args <- c(list(dat_dir = dat_dir,  
                   group_psm_by = group_psm_by, 
                   group_pep_by = group_pep_by, 
                   pt_cv = pt_cv, 
                   max_cv = max_cv, 
                   min_n = min_n, 
                   rm_allna = rm_allna, 
                   theme = theme), 
              dots)
    
    suppressMessages(
      parallel::clusterMap(
        cl, psm_mpurge, basenames,  
        MoreArgs = args)
    )
    
    parallel::stopCluster(cl)
  } 
  else {
    purrr::walk(basenames, psm_mpurge, 
                dat_dir = dat_dir, 
                group_psm_by = group_psm_by, 
                group_pep_by = group_pep_by, 
                pt_cv = pt_cv, 
                max_cv = max_cv, 
                min_n = min_n, 
                rm_allna = rm_allna, 
                theme = theme, 
                !!!dots)
  }
}


#'Purges peptide data
#'
#'\code{purgePep} removes \code{protein} entries from \code{Peptide.txt} by
#'quality criteria. It further plots the distributions of the CV of
#'\code{log2FC}.
#'
#'The CV of proteins under each sample are first calculated from contributing
#'peptides. In the event of multiple series of LC/MS injections, the CV of the
#'same protein from different LC/MS will be summarized by median statistics.
#'
#'The data nullification will be applied column-wisely for all available
#'samples. Argument \code{col_select} is merely used to subset samples for the
#'visualization of \code{log2FC} distributions.
#'
#'@inheritParams purgePSM
#'@inheritParams prnHist
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@inheritParams purge_by_qt
#'@inheritParams plot_prnTrend
#'@import dplyr ggplot2
#'@importFrom rlang exprs expr
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@example inst/extdata/examples/purgePep_.R
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
#'@export
purgePep <- function (dat_dir = NULL, pt_cv = NULL, max_cv = NULL, 
                      adjSD = FALSE, min_n = 1L, rm_allna = FALSE, 
                      col_select = NULL, col_order = NULL, filename = NULL, 
                      theme= NULL, ...) 
{
  on.exit(
    mget(names(formals()), envir = environment(), inherits = FALSE) %>% 
      c(rlang::enexprs(...)) %>% 
      save_call(fun), 
    add = TRUE
  )
  
  this_call <- match.call()
  fun <- as.character(this_call[[1]])
  
  check_dots(c("id", "anal_type", "filename", "filepath"), ...)
  
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% filter_dots]
  
  if (length(filter_dots)) {
    warning("`filter_` not applicable with `", fun, "`", call. = FALSE)
  }

  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  filename <- rlang::enexpr(filename)
  
  if (is.null(filename)) {
    filename <- "Peptide_sd.png"
    fn_prefix <- "Peptide_sd"
    fn_suffix <- "png"
  } 
  else {
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
  }
  
  if (is.null(dat_dir)) 
    dat_dir <- get_gl_dat_dir()
  else 
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  
  dir.create(file.path(dat_dir, "Peptide/log2FC_cv/purged"), 
             recursive = TRUE, showWarnings = FALSE)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(
    group_psm_by %in% c("pep_seq", "pep_seq_mod"), 
    group_pep_by %in% c("prot_acc", "gene"), 
    length(group_psm_by) == 1L, 
    length(group_pep_by) == 1L
  )

  load(file = file.path(dat_dir, "label_scheme_full.rda"))

  dir.create(file.path(dat_dir, "Peptide/Copy"), 
             recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(dat_dir, "Peptide/Peptide.txt"), 
            file.path(dat_dir, "Peptide/Copy/Peptide.txt"))
  
  fn <- file.path(dat_dir, "Peptide", "Peptide.txt")
  
  df <- suppressWarnings(
    readr::read_tsv(fn, col_types = get_col_types(), show_col_types = FALSE))
  
  pat_log2_R <- "^log2_R[0-9]{3}[NC]{0,1}"

  df <- df %>% 
    purge_by_qt(group_pep_by, pt_cv, min_n, rm_allna) %>% 
    purge_by_cv(group_pep_by, max_cv, min_n, rm_allna) %>% 
    purge_by_n(group_pep_by, min_n, rm_allna) %>% 
    { if (rm_allna) .[rowSums(!is.na(.[grepl(pat_log2_R, names(.))])) > 0L, ] 
      else . }

  # update
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
    purrr::reduce(dplyr::left_join, by = group_pep_by) %>% 
    dplyr::mutate(prot_n_psm.x = prot_n_psm.y, prot_n_pep.x = prot_n_pep.y) %>% 
    dplyr::rename(prot_n_psm = prot_n_psm.x, prot_n_pep = prot_n_pep.x) %>% 
    dplyr::select(-prot_n_psm.y, -prot_n_pep.y) %T>% 
    write.table(fn, sep = "\t", col.names = TRUE, row.names = FALSE)

  # plot
  width <- eval(dots$width, envir = rlang::caller_env())
  height <- eval(dots$height, envir = rlang::caller_env())
  
  if (is.null(width)) width <- 8 * n_TMT_sets(label_scheme_full)
  if (is.null(height)) height <- 8
  
  dots <- dots %>% .[! names(.) %in% c("width", "height")]
  
  filepath <- file.path(dat_dir, "Peptide/log2FC_cv/purged", filename)
  
  quiet_out <- purrr::quietly(sd_violin)(df = df, 
                                         id = !!group_pep_by, 
                                         filepath = filepath, 
                                         width = width, 
                                         height = height, 
                                         type = "log2_R", 
                                         adjSD = FALSE, 
                                         is_psm = FALSE, 
                                         col_select = !!col_select, 
                                         col_order = !!col_order, 
                                         theme = theme, 
                                         !!!dots)
}


