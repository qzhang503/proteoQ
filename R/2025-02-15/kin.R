#' Annotates kinase-substrate interactions
#' 
#' @inheritParams KinSubTest 
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>%
annot_KinSub <- function(df, db_nms, match_orgs, id, filepath = NULL, 
                         filename = NULL, ...) 
{
  dat_dir <- get_gl_dat_dir()
  
  if (!all(file.exists(db_nms))) {
    stop("Missing PSP file(s): \n", 
         purrr::reduce(db_nms %>% .[!file.exists(.)], paste, sep = "\n"), 
         "\nCheck the file path or download the database from 'www.phosphosite.org'.", 
         call. = FALSE)
  }
  
  if (all(is.na(df$entrez))) {
    stop("Entrez IDs not available in primary 'df' for annotation.\n", 
         "Consider utilities 'Uni2Entrez' or 'Ref2Entrez'for adding entrez IDs.", 
         call. = FALSE)
  }
  
  nms_df <- names(df)
  
  # --- add GENE_ID to the PSP table ---
  ok <- tryCatch(load(file = file.path(dat_dir, "acc_lookup.rda")),
                 error = function(e) "e")
  
  if (ok != "acc_lookup") {
    stop("`acc_lookup.rda` not found under ", dat_dir, ".", 
         call. = FALSE)
  }
  
  rm(list = "ok")

  acc_lookup <- acc_lookup %>% 
    dplyr::filter(!is.na(entrez), !is.na(gene), !duplicated(entrez)) %>% 
    dplyr::select(gene, entrez) 
  
  dbs <- readr::read_tsv(file.path(db_nms), 
                         col_types = cols(GENE = col_character(), 
                                          KIN_ORGANISM = col_character(),
                                          SUB_GENE_ID = col_character(), 
                                          SUB_GENE = col_character(),
                                          SUB_ORGANISM = col_character())) %>% 
    { if (match_orgs) dplyr::filter(., KIN_ORGANISM == SUB_ORGANISM) else . } %>% 
    dplyr::left_join(acc_lookup, by = c("GENE" = "gene")) %>% 
    dplyr::filter(!is.na(entrez)) %>% 
    dplyr::rename(GENE_ID = entrez) %>% 
    reloc_col_after("GENE_ID", "GENE") 
  
  # --- keep kinase and substrate proteins --- 
  universe <- c(dbs$GENE_ID, dbs$SUB_GENE_ID) %>% unique()

  df <- df %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id))) %>% 
    dplyr::filter(entrez %in% universe) %>% 
    dplyr::left_join(dbs, by = c("entrez" = "GENE_ID")) # %>% 
    # dplyr::filter(!is.na(SUB_GENE_ID)) # kinase may bot be other kinase's substrate

  df <- dplyr::bind_cols(
    df %>% dplyr::select(-which(names(.) %in% nms_df)), 
    df %>% dplyr::select(which(names(.) %in% nms_df))
  ) %T>% 
    readr::write_delim(file.path(filepath, filename), delim = find_delim(filename))
}


#' Analysis of kinase-substrate interactions
#'
#' The argument \code{scale_log2r} is not used in that both `_N` and `_Z`
#' columns from primary \code{df} will be kept. 
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnCorr_logFC
#' @inheritParams anal_KinSub
#' @import dplyr purrr fs
#' @importFrom magrittr %>% %T>% %$% %<>% 
KinSubTest <- function(df = NULL, id = gene, label_scheme_sub = NULL, 
                       db_nms = NULL, match_orgs = TRUE, 
                       scale_log2r = TRUE, complete_cases = FALSE, 
                       filepath = NULL, filename = NULL, ...) 
{
  dat_dir <- get_gl_dat_dir()
  id <- rlang::as_string(rlang::enexpr(id))
  
  stopifnot(id %in% names(df), nrow(label_scheme_sub) > 0L)

  if (complete_cases)
    df <- my_complete_cases(df, scale_log2r, label_scheme_sub)

  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  arrange_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter_dots, arrange_dots)]
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    rm_pval_whitespace()

  annot_KinSub(df, db_nms, match_orgs, id, filepath, filename)
}


#' PSP outputs of kinase-substrate interactions
#'
#' \code{anal_KinSub} adds the data of \href{https://www.phosphosite.org/}{PSP}
#' kinase-substrate interactions to peptide or protein results.
#'
#' OUtputs under folder \code{KinSub}.
#'
#' @inheritParams anal_prnString
#' @param db_nms  Character string(s) to the name(s) of PSP database(s) with
#'   prepended directory path(s). Users need to download the kinase-substrate
#'   table, e.g. \code{Kinase_Substrate_Dataset.txt} directly from the PSP
#'   website and supply the corresponding file path(s) and name(s). Currently
#'   assume single database file.
#' @param scale_log2r Not currently used. Values before and after scaling will
#'   be both reported.
#' @param filename The name of a output file.
#' @param type The type of data for annotation. The default is \code{peptide}.
#' @param match_orgs Logical; if TRUE, matches the organism between kinases and
#'   their acting substrates. The default is TRUE.
#' @param ... \code{filter_}: Variable argument statements for the row
#'   filtration against data in a primary file linked to \code{df}. See also
#'   \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'   \code{arrange_}: Variable argument statements for the row ordering against
#'   data in a primary file linked to \code{df}. See also \code{\link{prnHM}}
#'   for the format of \code{arrange_} statements.
#' @example inst/extdata/examples/KinSub_.R
#' @import dplyr purrr
#' @export
anal_KinSub <- function (db_nms = "~/proteoQ/dbs/psp/Kinase_Substrate_Dataset.txt", 
                         type = c("peptide", "protein"), match_orgs = TRUE, 
                         scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                         df = NULL, filepath = NULL, filename = NULL, ...) 
{
  on.exit(
    if (rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call("pepKinSub")
    } 
    else if (rlang::as_string(id) %in% c("prot_acc", "gene")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call("prnKinSub")
    }, 
    add = TRUE
  )
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  type <- rlang::enexpr(type)
  
  if (length(type) > 1) {
    type <- "peptide"
  } 
  else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("peptide", "protein"), 
              length(type) == 1L)
  }

  id <- switch(type, 
               peptide = match_call_arg(normPSM, group_psm_by), 
               protein = match_call_arg(normPSM, group_pep_by), 
               stop("`type` need to be one of `peptide` or `protein`.", 
                    call. = FALSE))
  
  dir.create(file.path(get_gl_dat_dir(), stringr::str_to_title(type), "KinSub/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(vapply(c(scale_log2r, complete_cases, impute_na), 
                   rlang::is_logical, logical(1L)))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            scale_log2r = TRUE, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            filepath = !!filepath, 
            filename = !!filename, 
            anal_type = "KinSub")(db_nms = db_nms, match_orgs = match_orgs, ...)
}


