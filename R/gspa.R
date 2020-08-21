#'GSPA of peptide data
#'
#'\code{pepGSPA} performs the analysis of Gene Set Probability Asymmetricity
#'(GSPA) against peptide \code{log2FC} data.
#'
#'@rdname prnGSPA
#'
#'@import purrr rlang dplyr
#'@export
pepGSPA <- function (gset_nms = c("go_sets", "c2_msig"), method = "mean", 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), 
                     min_size = 10, max_size = Inf, 
                     min_delta = 4, min_greedy_size = 1, 
                     use_adjP = FALSE, 
                     fml_nms = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {
  
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call(paste0("anal", "_pepGSPA"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call(paste0("anal", "_prnGSPA"))
    }
    , add = TRUE
  )  
  
  check_dots(c("id", "anal_type", "var_cutoff"), ...)
  check_gset_nms(gset_nms)
  
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Peptide/GSPA/log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), length(id) == 1)
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  method <- rlang::enexpr(method)
  if (rlang::is_call(method)) {
    method <- eval(method, env = caller_env())
  } else {
    method <- rlang::as_string(method)
  }
  stopifnot(all(method %in% c("mean", "limma")))
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSPA")
  
  reload_expts()
  
  info_anal(df = !!df, id = !!id, filepath = !!filepath, filename = !!filename, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            anal_type = "GSPA")(gset_nms = gset_nms, 
                                var_cutoff = 1000, 
                                pval_cutoff = pval_cutoff, 
                                logFC_cutoff = logFC_cutoff, 
                                gspval_cutoff = gspval_cutoff, 
                                gslogFC_cutoff = gslogFC_cutoff,
                                min_size = min_size, 
                                max_size = max_size, 
                                min_delta = min_delta, 
                                min_greedy_size = min_greedy_size, 
                                use_adjP = use_adjP, 
                                method = method, 
                                !!!dots)
}



#'GSPA of protein data
#'
#'\code{prnGSPA} performs the analysis of Gene Set Probability Asymmetricity
#'(GSPA) against protein \code{log2FC} data.
#'
#'The significance \code{pVals} of individual proteins are first obtained from
#'\code{\link{prnSig}}, followed by log10 transformation and separation into up-
#'or down-expressed groups for each gene set. At the default of \code{method =
#'mean}, the geometric mean of \code{pVals}, \eqn{P}, are each calculated for
#'the groups of up or down regulated proteins, with a penalty-like term
#'
#'\deqn{-log10(P)=(\sum_{i=1}^{n}-log10(p)+m)/(n+m)}{-log10(P)=(-log10(p_1*p_2*...*p_n)+m)/(n+m)}
#'
#'where \eqn{n} and \eqn{m} are the numbers of entries with \eqn{p} values
#'\eqn{\le} or \eqn{>} a significance cut-off, respectively. The quotient of the
#'two \eqn{P} values is then used to represent the significance of gene set
#'enrichment. The arguments \code{pval_cutoff} and \code{logFC_cutoff} are used
#'to discriminate low influence genes. Additional subsetting of data via the
#'\code{vararg} approach of \code{filter_} is feasible. At \code{method =
#'limma}, moderated t-tests are performed against \code{-log10(pVals)} between
#'the up and the down groups via \code{\link[limma]{eBayes}}.
#'
#'@inheritParams anal_pepNMF
#'@inheritParams prnHist
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param gset_nms Character string or vector containing the shorthanded name(s),
#'  full file path(s) or both to gene sets for enrichment analysis. For species
#'  among \code{"human", "mouse", "rat"}, the default of \code{c("go_sets",
#'  "c2_msig")} will utilize terms from both gene ontology (\code{GO}) and
#'  molecular signatures (\code{MSig}). Custom data bases of \code{GO} and
#'  curated \code{MSig}, and/or additional species are also supported. See also
#'  \code{\link{prepGO}} for the preparation of custom \code{GO} and
#'  \code{\link{prepMSig}} for the preparation of custom \code{MSig}.
#'@param method Character string or vector; the method to assess the p-values of
#'  GSPA. The default is \code{mean} and the alternative is \code{limma}. See
#'  also section \code{Details} for the calculations.
#'@param pval_cutoff Numeric value or vector; the cut-off in protein
#'  significance \code{pVal}. Entries with \code{pVals} less significant than
#'  the threshold will be excluded from enrichment analysis. The default is 0.05
#'  for all formulas matched to or specified in argument \code{fml_nms}.
#'  Formula-specific threshold is allowed by supplying a vector of cut-off
#'  values.
#'@param logFC_cutoff Numeric value or vector; the cut-off in protein
#'  \code{log2FC}. Entries with absolute \code{log2FC} smaller than the
#'  threshold will be excluded from enrichment analysis. The default magnitude
#'  is \code{log2(1.2)} for all formulas matched to or specified in argument
#'  \code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
#'  of absolute values in \code{log2FC}.
#'@param gspval_cutoff Numeric value or vector; the cut-off in gene-set
#'  significance \code{pVal}. Only enrichment terms with \code{pVals} more
#'  significant than the threshold will be reported. The default is 0.05 for all
#'  formulas matched to or specified in argument \code{fml_nms}.
#'  Formula-specific threshold is allowed by supplying a vector of cut-off
#'  values.
#'@param gslogFC_cutoff Numeric value or vector; the cut-off in gene-set
#'  enrichment fold change. Only enrichment terms with absolute fold change
#'  greater than the threshold will be reported. The default magnitude is
#'  \code{log2(1.2)} for all formulas matched to or specified in argument
#'  \code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
#'  of absolute values in \code{log2FC}.
#'@param min_size Numeric value or vector; minimum number of protein entries for
#'  consideration in gene set tests. The number is after data filtration by
#'  \code{pval_cutoff}, \code{logFC_cutoff} or varargs expressions under
#'  \code{filter_}. The default is 10 for all formulas matched to or specified
#'  in argument \code{fml_nms}. Formula-specific threshold is allowed by
#'  supplying a vector of sizes.
#'@param max_size Numeric value or vector; maximum number of protein entries for
#'  consideration in gene set tests. The number is after data filtration by
#'  \code{pval_cutoff}, \code{logFC_cutoff} or varargs expressions under
#'  \code{filter_}. The default in infinite for all formulas matched to or
#'  specified in argument \code{fml_nms}. Formula-specific threshold is allowed
#'  by supplying a vector of sizes.
#'@param min_delta Numeric value or vector; the minimum count difference between
#'  the up- and the down-expressed group of proteins for consideration in gene
#'  set tests. For example at \code{min_delta = 4}, a gene set will 6
#'  upregulated proteins and 2 down-expressed proteins, or vice versa, will be
#'  assessed. The number is after data filtration by \code{pval_cutoff},
#'  \code{logFC_cutoff} or varargs expressions under \code{filter_}. The default
#'  is 4 for all formulas matched to or specified in argument \code{fml_nms}.
#'  Formula-specific threshold is allowed by supplying a vector of sizes.
#'@param min_greedy_size Numeric value or vector; minimum number of unique
#'  protein entries for a gene set to be considered essential. The default in
#'  \code{1} for all formulas matched to or specified in argument
#'  \code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
#'  of sizes.
#'@param use_adjP Logical; if TRUE, use Benjamini-Hochberg pVals. The default is
#'  FALSE.
#'@param fml_nms Character string or vector; the formula name(s). By default,
#'  the formula(s) will match those used in \code{\link{pepSig}} or
#'  \code{\link{prnSig}}.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration
#'  against data in a primary file of \code{/Model/Protein[_impNA]_pVals.txt}.
#'  See also \code{\link{normPSM}} for the format of \code{filter_} statements.
#'  \cr \cr \code{arrange_}: Variable argument statements for the row ordering
#'  against data in a primary file of \code{/Model/Protein[_impNA]_pVals.txt}.
#'  See also \code{\link{prnHM}} for the format of \code{arrange_} statements.
#'@import dplyr rlang ggplot2 
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'
#'@example inst/extdata/examples/prnGSPA_.R
#'
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
#'@export
prnGSPA <- function (gset_nms = c("go_sets", "c2_msig"), method = "mean", 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), 
                     min_size = 10, max_size = Inf, 
                     min_delta = 4, min_greedy_size = 1, 
                     use_adjP = FALSE, 
                     fml_nms = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {

  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call(paste0("anal", "_pepGSPA"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call(paste0("anal", "_prnGSPA"))
    }
    , add = TRUE
  )  
  
  check_dots(c("id", "anal_type", "var_cutoff"), ...)
  check_gset_nms(gset_nms)
  
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Protein/GSPA/log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)

  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  method <- rlang::enexpr(method)
  if (rlang::is_call(method)) {
    method <- eval(method, env = caller_env())
  } else {
    method <- rlang::as_string(method)
  }
  stopifnot(all(method %in% c("mean", "limma")))
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSPA")

  reload_expts()
  
  info_anal(df = !!df, id = !!id, filepath = !!filepath, filename = !!filename, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            anal_type = "GSPA")(gset_nms = gset_nms, 
                                var_cutoff = 1000, 
                                pval_cutoff = pval_cutoff, 
                                logFC_cutoff = logFC_cutoff, 
                                gspval_cutoff = gspval_cutoff, 
                                gslogFC_cutoff = gslogFC_cutoff,
                                min_size = min_size, 
                                max_size = max_size, 
                                min_delta = min_delta, 
                                min_greedy_size = min_greedy_size, 
                                use_adjP = use_adjP, 
                                method = method, 
                                !!!dots)
}


#' Perform GSPA tests
#'
#' logFC_cutoff Numeric A threshold for the subset of data before the calculation
#' of adjusted pvals
#'
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams prnGSPA
#' @inheritParams prnGSVA
#' @inheritParams info_anal
#' @param id Currently only "entrez".
#' @param label_scheme_sub A data frame. Subset entries from \code{label_scheme}
#'   for selected samples.
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %T>% %$% %<>% 
gspaTest <- function(df = NULL, id = "entrez", label_scheme_sub = NULL, 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                     filepath = NULL, filename = NULL, 
                     gset_nms = "go_sets", var_cutoff = 0.5, 
                     pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1), 
                     min_size = 6, max_size = Inf, min_delta = 4, min_greedy_size = 1, 
                     use_adjP = FALSE, 
                     method = "mean", anal_type = "GSPA", ...) {

  # `GSPA`: use `pVals` instead of `N_log2R` or `Z_log2R`; 
  # currently Protein[_impNA]_pVals.txt does not contain `scale_log2_r` info in the file name
  #   and cannot tell `pVals` are from `N_log2R` or `Z_log2R`;  
  #   therefore, `scale_log2r` will be matched to those in `prnSig()` and indicated in 
  #   the names of out files as `_N[_impNA].txt` or `_Z[_impNA].txt` to inform user the `scale_log_r` status
  # `complete_cases` is for entrez IDs, pVals, log2FC
  # "id" for tibbling rownames
  
  stopifnot(vapply(c(pval_cutoff, logFC_cutoff, min_size, max_size, min_delta, min_greedy_size), 
                   rlang::is_double, logical(1)))
  stopifnot(nrow(label_scheme_sub) > 0)
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots = rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]

  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  if (nrow(df) == 0) stop("No data available after row filtration.", call. = FALSE)
  
  if (purrr::is_empty(fmls)) stop("Formula(s) of contrasts not available.", call. = FALSE)

  species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
  gsets <- load_dbs(gset_nms = gset_nms, species = species)
  
  stopifnot(length(gsets) > 0)

  fml_nms <- names(df) %>% 
    .[grepl("pVal\\s*\\(", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique() %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls %>% .[names(.) %in% fml_nms]
  fml_nms <- fml_nms %>% .[map_dbl(., ~ which(.x == names(fmls)))]

  if (purrr::is_empty(fml_nms)) {
    stop("No formula matached; compare the formula name(s) with those in `prnSig(..)`")
  }

  col_ind <- purrr::map(fml_nms, ~ grepl(.x, names(df))) %>%
    `names<-`(paste0("nm_", seq_along(.))) %>% 
    dplyr::bind_cols() %>%
    rowSums() %>%
    `>`(0)
  
  stopifnot(sum(col_ind) > 0)

  switch (anal_type,
          GSPA = purrr::pwalk(list(fmls, 
                                   fml_nms, 
                                   pval_cutoff, 
                                   logFC_cutoff, 
                                   gspval_cutoff, 
                                   gslogFC_cutoff, 
                                   min_size, 
                                   max_size, 
                                   min_delta, 
                                   min_greedy_size, 
                                   method
                              ), 
                              fml_gspa, 
                              df = df, 
                              col_ind = col_ind, 
                              id = !!id, 
                              gsets = gsets, 
                              label_scheme_sub = label_scheme_sub, 
                              complete_cases = complete_cases, 
                              scale_log2r = scale_log2r, 
                              filepath = filepath, 
                              filename = filename, 
                              use_adjP = use_adjP, 
                              !!!filter_dots), 
          GSEA = purrr::pwalk(list(fmls, 
                                   fml_nms, 
                                   var_cutoff, 
                                   pval_cutoff, 
                                   logFC_cutoff, 
                                   gspval_cutoff, 
                                   gslogFC_cutoff, 
                                   min_size, 
                                   max_size
                              ), 
                              fml_gsea, 
                              df = df, 
                              col_ind = col_ind, 
                              id = !!id, 
                              gsets = gsets, 
                              label_scheme_sub = label_scheme_sub, 
                              complete_cases = complete_cases, 
                              scale_log2r = scale_log2r, 
                              filepath = filepath, 
                              filename = filename, 
                              !!!dots), 
    stop("Unhandled task.")
  )
  
  invisible(gset_nms)
}


#' A helper function for GSPA
#' 
#' @param fml A character string; the formula used in \link{prnSig}.
#' @param fml_nm A character string; the name of \code{fml}.
#' @param col_ind Numeric vector; the indexes of columns for the ascribed \code{fml_nm}.
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams gspaTest
#' @inheritParams gsVolcano
#' @import purrr dplyr rlang
#' @importFrom magrittr %>% %T>% %$% %<>% 
fml_gspa <- function (fml, fml_nm, pval_cutoff, logFC_cutoff, gspval_cutoff, gslogFC_cutoff, 
                      min_size, max_size, min_delta, min_greedy_size, method, 
                      df, col_ind, id, gsets, label_scheme_sub, complete_cases, scale_log2r, 
                      filepath, filename, use_adjP = FALSE, ...) {

  on.exit(
    pars <- mget(names(formals()), current_env()) %>% 
      .[! names(.) %in% c("gsets", "df", "label_scheme_sub", "col_ind")] %>% 
      c(rlang::enexprs(...)) %>% 
      save_call(paste0(fn_prefix, "@", fml_nm))
  )
  
  dir.create(file.path(filepath, fml_nm), recursive = TRUE, showWarnings = FALSE)
  id <- rlang::as_string(rlang::enexpr(id))
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  df <- df %>% prep_gspa(id = "entrez", fml_nm = fml_nm, col_ind = col_ind, 
                         pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, 
                         use_adjP = use_adjP) 
  
  if (complete_cases) {
    df <- df[complete.cases(df), ]
    if (!purrr::is_empty(names(df))) 
      message("Complete cases against columns: ", purrr::reduce(names(df), paste, sep = ", "))
  }
  
  df <- df %>% 
    dplyr::mutate(p_val = -log10(p_val)) %>% 
    dplyr::mutate(valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>% 
    dplyr::mutate(valence = factor(valence, levels = c("neg", "pos"))) %>% 
    tidyr::complete(entrez, contrast, valence)
  
  # do not re-`arrange` df, df_sub for limma after this point
  df <- df %>% dplyr::arrange(entrez, contrast, valence)

  gsets <- gsets %>% 
    .[!grepl("molecular_function$", names(.))] %>% 
    .[!grepl("cellular_component$", names(.))] %>% 
    .[!grepl("biological_process$", names(.))] %>%     
    .[purrr::map_lgl(., ~ length(.x) >= min_size)] %>% 
    .[purrr::map_lgl(., ~ length(.x) <= max_size)]
  
  res_pass <- local({
    contrast_groups <- unique(df$contrast) %>% as.character()
    
    contrast_table <- tibble(contrast = rep(contrast_groups, 2), 
                             valence = rep(c("neg", "pos"), length(contrast_groups)), 
                             p_val = rep(-log10(0.05), 2 * length(contrast_groups))) %>% 
      dplyr::mutate(contrast = factor(contrast, levels = levels(df$contrast)), 
                    valence = factor(valence, levels(df$valence)))
    
    res <- switch(method, 
                  limma = purrr::map(gsets, 
                                     gspa_summary_limma, 
                                     df = df, 
                                     min_size = min_size, 
                                     min_delta = min_delta, 
                                     gspval_cutoff = gspval_cutoff, 
                                     gslogFC_cutoff = gslogFC_cutoff), 
                  mean = purrr::map(gsets, 
                                    gspa_summary_mean, 
                                    df = df, 
                                    min_size = min_size, 
                                    min_delta = min_delta, 
                                    gspval_cutoff = gspval_cutoff, 
                                    gslogFC_cutoff = gslogFC_cutoff), 
                  rlang::abort("Unknown `method`."))

    idx <- purrr::map_dbl(res, is.null)
    
    res <- res[!idx] %>% 
      do.call(rbind, .) %>% 
      tibble::rownames_to_column("term") %>% 
      dplyr::mutate(term = factor(term)) %>% 
      dplyr::mutate_at(.vars = grep("^pVal\\s+\\(", names(.)), ~ abs(.x)) %>% 
      dplyr::group_by(term) %>% 
      dplyr::arrange(term) %>% 
      data.frame(check.names = FALSE) 
    
    pass <- purrr::map(res[paste0("pVal (", contrast_groups, ")")], ~ .x <= gspval_cutoff) %>%
      dplyr::bind_cols() %>%
      rowSums() %>%
      `>`(0)
    
    res_pass <- res[pass, ] %>% 
      dplyr::mutate_at(vars(grep("^pVal", names(.))), ~ replace(.x, .x < 1E-50, 1E-50))
  })
  
  if (nrow(res_pass) == 0) return(NULL)
  
  if (!requireNamespace("RcppGreedySetCover", quietly = TRUE)) {
    stop("\n============================================================================", 
         "\nNeed install package \"RcppGreedySetCover\" needed for this function to work.",
         "\n============================================================================",
         call. = FALSE)
  }

  out <- local({
    adjp <- res_pass %>% 
      dplyr::select(grep("^pVal\\s+\\(", names(.))) %>% 
      purrr::map(~ p.adjust(., "BH")) %>% 
      data.frame(check.names = FALSE) %>% 
      `names<-`(gsub("pVal", "q_val", colnames(.))) %>% 
      `names<-`(gsub("^q_val\\s+\\((.*)\\)", "\\1", names(.))) %>% 
      tidyr::gather(key = contrast, value = q_val) %>% 
      dplyr::select(-contrast)
    
    log2fc <- res_pass %>% 
      dplyr::select(grep("^log2Ratio\\s+\\(", names(.))) %>% 
      `names<-`(gsub("^log2Ratio\\s+\\((.*)\\)", "\\1", names(.))) %>% 
      tidyr::gather(key = contrast, value = log2fc) %>% 
      dplyr::select(-contrast)
    
    pval <- res_pass %>% 
      dplyr::select(-grep("^log2Ratio\\s+\\(", names(.))) %>% 
      `names<-`(gsub("^pVal\\s+\\((.*)\\)", "\\1", names(.))) %>% 
      tidyr::gather(key = contrast, value = p_val, -term)
    
    dplyr::bind_cols(pval, adjp, log2fc) %>% 
      dplyr::mutate_at(.vars = grep("^p_val$|^adjP$", names(.)), format, 
                       scientific = TRUE, digits = 2) %>%
      dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2)
  })

  sig_sets <- purrr::map2(gsets, names(gsets), ~ {
    tibble::tibble(id = as.character(.x)) %>% 
      dplyr::mutate(set = .y)
  }) %>% 
    dplyr::bind_rows() %>% 
    `names<-` (c("id", "term")) %>% 
    dplyr::filter(id %in% unique(df$entrez)) %>% 
    dplyr::select(c("term", "id")) %>% 
    dplyr::filter(term %in% res_pass$term)
  
  out <- out %>% dplyr::mutate(term = as.character(term))
  
  out <- sig_sets %>% 
    dplyr::group_by(term) %>% 
    dplyr::summarise(size = n()) %>% # the number of entries matched to input `df`
    dplyr::right_join(out, by = "term")

  res_greedy <- sig_sets %>% 
    RcppGreedySetCover::greedySetCover(FALSE) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_resgreedy.txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE) %>% 
    dplyr::group_by(term) %>% 
    dplyr::summarise(ess_size = n()) %>% 
    dplyr::filter(ess_size >= min_greedy_size)
  
  sig_sets <- res_greedy %>% 
    dplyr::right_join(out, by = "term") %>% 
    dplyr::mutate(is_essential = ifelse(!is.na(ess_size), TRUE, FALSE)) %>% 
    dplyr::select(term, is_essential, size, ess_size, contrast, p_val, q_val, log2fc) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, ".txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE) %>% 
    dplyr::filter(is_essential) %>% 
    dplyr::filter(!duplicated(term)) %>% 
    dplyr::select(-contrast, -p_val, -q_val, -log2fc) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_essmeta.txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE) %>% 
    dplyr::select(term, ess_size) %>% 
    dplyr::right_join(sig_sets, by = "term")

  map_essential(sig_sets) %>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_essmap.txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE)	
}


#' check the size of a gene set
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
ok_min_size <- function (df, min_delta, gspval_cutoff = 1, gslogFC_cutoff = 0) {
  data <- df %>% 
    dplyr::group_by(contrast, valence) %>% 
    dplyr::filter(!is.na(log2Ratio), !is.na(p_val)) %>% 
    dplyr::mutate(n = n()) %>% 
    tidyr::unite(con_val, contrast, valence, sep = ".", remove = FALSE) 
  
  delta_p <- data %>% 
    dplyr::summarise(p_val = mean(p_val, na.rm = TRUE)) %>% 
    dplyr::mutate(p_val = ifelse(is.nan(p_val), 0, p_val)) %>% 
    tidyr::spread(key = contrast, value = p_val) %>% 
    dplyr::select(-valence)
  delta_p[is.na(delta_p)] <- 0
  if (nrow(delta_p) == 2) delta_p <- delta_p[2, ] - delta_p[1, ]

  delta_fc <- data %>% 
    dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) %>% 
    dplyr::mutate(log2Ratio = ifelse(is.nan(log2Ratio), 0, log2Ratio)) %>% 
    tidyr::spread(key = contrast, value = log2Ratio) %>% 
    dplyr::select(-valence)
  delta_fc[is.na(delta_fc)] <- 0
  if (nrow(delta_fc) == 2) delta_fc <- delta_fc[2, ] + delta_fc[1, ]
  
  delta_n <- data %>% 
    dplyr::summarise(n = mean(n, na.rm = TRUE)) %>% # empty levels will be kept
    tidyr::spread(key = contrast, value = n) %>% 
    dplyr::select(-valence)
  delta_n[is.na(delta_n)] <- 0
  if (nrow(delta_n) == 2) delta_n <- delta_n[2, ] - delta_n[1, ]

  oks <- sign(delta_p) == sign(delta_n) & 
    abs(delta_p) >= -log10(gspval_cutoff) & 
    abs(delta_n) >= min_delta & 
    abs(delta_fc) >= gslogFC_cutoff
  
  if (length(oks) < nlevels(df$contrast)) {
    levs <- levels(df$contrast)
    missing_nms <- levs[!levs %in% colnames(oks)]
    oks <- dplyr::bind_cols(!!missing_nms := FALSE, oks)
    oks <- oks[, levs, drop = FALSE]
  }
  
  return(oks)  
}


#' gspa pVal calculations using geomean
#' 
#' @param gset Character string; a gene set indicated by \code{gset_nms}.
#' @inheritParams prnHist
#' @inheritParams prnGSPA
gspa_summary_mean <- function(gset, df, min_size = 10, min_delta = 4, 
                              gspval_cutoff = 0.05, gslogFC_cutoff = log2(1)) {
  df <- df %>% dplyr::filter(.[["entrez"]] %in% gset)
  
  if (length(unique(df$entrez)) < min_size) return(NULL)
  
  df <- df %>% dplyr::group_by(contrast, valence)

  ok_delta_n <- ok_min_size(df = df, 
                            min_delta = min_delta, 
                            gspval_cutoff = gspval_cutoff, 
                            gslogFC_cutoff = gslogFC_cutoff)
  
  # penalty term
  dfw <- df %>% dplyr::mutate(p_val = ifelse(is.na(p_val), 1, p_val))

  delta_p <- dfw %>% 
    dplyr::summarise(p_val = mean(p_val, na.rm = TRUE)) %>% 
    dplyr::mutate(p_val = replace(p_val, is.nan(p_val), 1)) %>% 
    tidyr::spread(key = contrast, value = p_val) %>% 
    dplyr::select(-valence) %>% 
    `colnames<-`(paste0("pVal (", colnames(.), ")")) 
  if (nrow(delta_p) == 2) delta_p <- delta_p[2, ] - delta_p[1, ]
  delta_p <- abs(delta_p)
  delta_p <- 1/10^delta_p
  delta_p[] <- purrr::map2(as.list(delta_p), as.list(ok_delta_n), ~ {if (.y) .x else 1}) 

  delta_fc <- dfw %>% 
    dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) %>% 
    dplyr::mutate(log2Ratio = replace(log2Ratio, is.nan(log2Ratio), 0)) %>% 
    tidyr::spread(key = contrast, value = log2Ratio) %>% 
    dplyr::select(-valence) %>% 
    `colnames<-`(paste0("log2Ratio (", colnames(.), ")")) %>% 
    abs(.) 
  if (nrow(delta_fc) == 2) delta_fc <- delta_fc[2, ] - delta_fc[1, ]
  
  return(dplyr::bind_cols(delta_p, delta_fc))
}  


#' gspa pVal calculations using limma
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams gspa_summary_mean
gspa_summary_limma <- function(gset, df, min_size = 10, min_delta = 4, 
                               gspval_cutoff = 0.05, gslogFC_cutoff = 1.2) {
  df_sub <- df %>% dplyr::filter(.[["entrez"]] %in% gset) 
  if (length(unique(df_sub$entrez)) < min_size) return(NULL)
  
  lm_gspa(df_sub, min_delta, gspval_cutoff, gslogFC_cutoff)
}  


#' significance tests of pVals between the up and the down groups
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
lm_gspa <- function(df, min_delta, gspval_cutoff, gslogFC_cutoff) {
  delta_fc <- df %>% 
    dplyr::group_by(contrast, valence) %>% 
    dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) %>% 
    dplyr::mutate(log2Ratio = replace(log2Ratio, is.nan(log2Ratio), 0)) %>% 
    tidyr::spread(key = contrast, value = log2Ratio) %>% 
    dplyr::select(-valence) 
  
  nms <- colnames(delta_fc)
  delta_fc <- delta_fc%>% 
    `colnames<-`(paste0("log2Ratio (", colnames(.), ")")) %>% 
    abs(.) 
  if (nrow(delta_fc) == 2) delta_fc <- delta_fc[2, ] - delta_fc[1, ]
  
  ok_delta_n <- ok_min_size(df = df, 
                            min_delta = min_delta, 
                            gspval_cutoff = gspval_cutoff, 
                            gslogFC_cutoff = gslogFC_cutoff)

  data <- df %>% 
    dplyr::group_by(contrast, valence) %>% 
    dplyr::select(-log2Ratio) %>% 
    tidyr::unite(sample, entrez, contrast, valence, sep = ".", remove = TRUE) %>% 
    tidyr::spread(key = sample, value = p_val) 
  
  label_scheme_gspa <- tibble(Sample_ID = names(data), Group = Sample_ID) %>% 
    dplyr::mutate(Group = gsub("^\\d+\\.", "", Group)) %>% 
    dplyr::mutate(Group = factor(Group))
  
  design <- model.matrix(~0+label_scheme_gspa$Group) %>% 
    `colnames<-`(levels(label_scheme_gspa$Group)) %>% 
    data.frame(check.names = TRUE)

  neg_cols <- colnames(design)[seq_along(design) %% 2 %>% as.logical()]
  pos_cols <- colnames(design)[!(seq_along(design) %% 2)]
  
  contrs <- paste(pos_cols, neg_cols, sep = "-")
  
  contr_mat <- makeContrasts(contrasts = contrs, levels = data.frame(design)) %>% 
    `colnames<-`(contrs) %>% 
    `rownames<-`(colnames(design))
  
  fit <- data %>%
    lmFit(design = design) %>%
    contrasts.fit(contr_mat) %>%
    eBayes()
  
  delta_p <- fit$p.value %>%
    replace(., is.na(.), 1) %>% 
    data.frame(check.names = FALSE) %>%
    `colnames<-`(paste0("pVal (", nms, ")"))

  delta_p[] <- purrr::map2(as.list(delta_p), as.list(ok_delta_n), ~ {if (.y) .x else 1}) 

  return(dplyr::bind_cols(delta_p, delta_fc))
}


#' A helper function for GSPA
#'
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams fml_gspa
#' 
#' @import purrr dplyr rlang
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom readr read_tsv
prep_gspa <- function(df, id, fml_nm, col_ind, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                      use_adjP = FALSE) {
  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- df %>%
    dplyr::select(grep(fml_nm, names(.), fixed = TRUE)) %>%
    `colnames<-`(gsub(paste0(fml_nm, "."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .) %>% 
    rm_pval_whitespace() %>% 
    dplyr::select(id, grep("^pVal|^adjP|^log2Ratio", names(.))) %>% 
    dplyr::mutate(!!id := as.character(.[[id]]))
  
  if (use_adjP) {
    df <- df %>% 
      dplyr::select(-grep("pVal", names(.))) %>% 
      `names<-`(gsub("adjP", "pVal", names(.)))
  } else {
    df <- df %>% 
      dplyr::select(-grep("adjP", names(.)))
  }
  
  contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
    gsub("^log2Ratio\\s+\\(|\\)$", "", .)
  
  pvals <- df %>% 
    dplyr::select(-grep("^log2Ratio\\s+\\(", names(.))) %>% 
    `names<-`(gsub("^pVal\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = contrast, value = p_val, -id)
  
  log2rs <- df %>% 
    dplyr::select(-grep("^pVal\\s+\\(", names(.))) %>% 
    `names<-`(gsub("^log2Ratio\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = contrast, value = log2Ratio, -id) %>% 
    dplyr::select(-id, -contrast)
  
  df <- dplyr::bind_cols(pvals, log2rs) %>% 
    dplyr::filter(p_val <= pval_cutoff) %>% # NA will be removed
    dplyr::filter(abs(log2Ratio) >= logFC_cutoff) %>% 
    dplyr::filter(!is.na(id)) %>% 
    dplyr::mutate(contrast = factor(contrast, levels = contrast_groups)) %>%
    dplyr::arrange(contrast) 

  return(df)
}


#' A helper function for mapping between gene sets and essential gene sets
#'
#' @param sig_sets A data frame containing the gene sets that are significant
#'   under given criteria.
#' @import purrr dplyr rlang 
#' @importFrom magrittr %>% %T>% %$% %<>% 
map_essential <- function (sig_sets) {
  ess_terms <- sig_sets %>% 
    dplyr::filter(!is.na(ess_size), !duplicated(term)) %>% 
    dplyr::select(term) %>% 
    unlist

  sig_greedy_sets <- sig_sets %>% dplyr::filter(term %in% ess_terms)

  res <- purrr::map(unique(sig_sets$term), ~ {
    curr_sig_set <- sig_sets %>% 
      dplyr::filter(term %in% .x)
    
    sig_greedy_sets %>% 
      dplyr::filter(id %in% curr_sig_set$id) %>% 
      dplyr::group_by(term) %>%
      dplyr::summarise(size = n()) %>% 
      dplyr::left_join(curr_sig_set %>% 
                         dplyr::select(c("term", "ess_size")) %>% 
                         dplyr::filter(!duplicated(term)), 
                       by = "term") %>% 
      dplyr::mutate(fraction = size/nrow(curr_sig_set)) %>% 
      dplyr::mutate(curr_sig_set = .x)
  }, sig_sets, sig_greedy_sets) %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(curr_sig_set) %>%
    dplyr::rename(ess_term = "term", term = "curr_sig_set") %>%
    dplyr::select(c("term", "ess_term", "size", "ess_size", "fraction")) %>% 
    dplyr::mutate(fraction = round(fraction, digits = 3)) %>% 
    # dplyr::ungroup(ess_term) %>% 
    dplyr::mutate(distance = 1 - fraction) %>% 
    dplyr::mutate(term = factor(term, levels = unique(as.character(term)))) %>% 
    dplyr::mutate(ess_term = factor(ess_term, levels = levels(term))) %>%
    dplyr::mutate(idx = as.numeric(term), ess_idx = as.numeric(ess_term))
}


#'Heat map visualization of GSPA results
#'
#'\code{prnGSPAHM} visualizes distance heat maps and networks between essential
#'and all gene sets.
#'
#'The list of gene sets and the associative quality metrics of \code{size} and
#'\code{ess_size} are assessed after data filtration with the criteria specified
#'by arguments \code{pval_cutoff} and \code{logFC_cutoff}, as well as optional
#'varargs of \code{filter_}.
#'
#'@section \code{Protein_GSPA_[...].txt}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Description}\cr term \tab a gene
#'  set term \cr is_essential \tab a logical indicator of gene set essentiality
#'  \cr size \tab the number of IDs under a \code{term} \cr ess_size \tab the
#'  number of IDs that can be found under a corresponding essential set \cr
#'  contrast \tab a contrast of sample groups \cr p_val \tab significance p
#'  values \cr q_val \tab \code{p_val} with \code{BH} adjustment of multiple
#'  tests \cr log2fc \tab the fold change of a gene set at logarithmic base of 2
#'  \cr }
#'
#'@section \code{Protein_GSPA_[...]essmap.txt}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Descrption}\cr term \tab a gene
#'  set term \cr ess_term \tab an essential gene set term \cr size \tab the
#'  number of IDs under a \code{term} with matches to an \code{ess_term} \cr
#'  ess_size \tab the number of essential IDs under a \code{term} with matches
#'  to an \code{ess_term} \cr fraction \tab a fraction of matches in IDs between
#'  a \code{term} and a \code{ess_term} \cr distance \tab 1 - \code{fraction}
#'  \cr idx \tab a numeric index of \code{term} \cr ess_idx \tab a numeric index
#'  of \code{ess_term} \cr }
#'
#'@inheritParams plot_prnTrend
#'@inheritParams  prnEucDist
#'@param impute_na Logical; at TRUE, input files with \code{_impNA[...].txt} in
#'  name will be loaded. Otherwise, files without \code{_impNA} in name will be
#'  taken. An error will be thrown if no files are matched under given
#'  conditions. The default is FALSE.
#'@param fml_nms Character string or vector; the formula name(s). By default,
#'  the formula(s) will match those used in \code{\link{pepSig}} or
#'  \code{\link{prnSig}}.
#'@param annot_cols A character vector of column keys that can be found in
#'  \code{_essmap.txt}. The values under the selected keys will be used to
#'  color-code enrichment terms on the top of heat maps. The default is NULL
#'  without column annotation.
#'@param annot_rows A character vector of column keys that can be found from
#'  \code{_essmeta.txt} . The values under the selected keys will be used to
#'  color-code essential terms on the side of heat maps. The default is NULL
#'  without row annotation.
#'@param ... \code{filter2_}: Variable argument statements for the row
#'  filtration against data in secondary file(s) of \code{_essmap.txt}. Each
#'  statement contains to a list of logical expression(s). The \code{lhs} needs
#'  to start with \code{filter2_}. The logical condition(s) at the \code{rhs}
#'  needs to be enclosed in \code{exprs} with round parenthesis. For example,
#'  \code{distance} is a column key in \code{Protein_GSPA_Z_essmap.txt}. The
#'  statement \code{filter2_ = exprs(distance <= .95),} will remove entries with
#'  \code{distance > 0.95}. See also \code{\link{normPSM}} for the format of
#'  \code{filter_} statements against primary data. \cr \cr \code{arrange2_}:
#'  Variable argument statements for the row ordering against data in secondary
#'  file(s) of \code{_essmap.txt}. The \code{lhs} needs to start with
#'  \code{arrange2_}. The expression(s) at the \code{rhs} needs to be
#'  enclosed in \code{exprs} with round parenthesis. For example,
#'  \code{distance} and \code{size} are column keys in
#'  \code{Protein_GSPA_Z_essmap.txt}. The statement \code{arrange2_ =
#'  exprs(distance, size),} will order entries by \code{distance}, then by
#'  \code{size}. See also \code{\link{prnHM}} for the format of \code{arrange_}
#'  statements against primary data. \cr \cr Additional arguments for
#'  \code{\link[pheatmap]{pheatmap}}, i.e., \code{fontsize }... \cr \cr Note
#'  arguments disabled from \code{pheatmap}: \cr \code{annotation_col}; instead
#'  use keys indicated in \code{annot_cols} \cr \code{annotation_row}; instead
#'  use keys indicated in \code{annot_rows}
#'@import purrr
#'
#'@example inst/extdata/examples/prnGSPAHM_.R
#'
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
#'@export
prnGSPAHM <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, fml_nms = NULL, 
                       annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                       df2 = NULL, filename = NULL, ...) {
  check_dots(c("id", "anal_type", "df", "filepath"), ...)
  
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Protein/GSPA/log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)

  df2 <- rlang::enexpr(df2)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  annot_rows <- rlang::enexpr(annot_rows) 
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSPA")
  
  reload_expts()

  info_anal(id = !!id, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, df2 = !!df2, filepath = NULL, filename = !!filename, 
            anal_type = "GSPA_hm")(annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, 
                                   annot_rows = !!annot_rows, !!!dots)
}


#' Plot distance heat map of GSPA
#' 
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams gspaMap
gspaHM <- function(scale_log2r, complete_cases, impute_na, df2, filepath, filename, ...) {
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  if (purrr::is_empty(fmls))
    stop("No formula(s) of contrasts available.", call. = TRUE)
  
  fml_nms <- names(fmls)
  if (length(fml_nms) > 0) {
    purrr::walk(fml_nms, byfml_gspahm, df2, filepath, filename, scale_log2r, impute_na, !!!dots)
  }
}


#' A helper function for gspaHM by formula names
#' 
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import purrr dplyr rlang pheatmap 
#' @importFrom magrittr %>% %T>% %$% %<>% 
byfml_gspahm <- function (fml_nm, df2, filepath, filename, scale_log2r, impute_na, ...) {
  ins <- list.files(path = file.path(filepath, fml_nm), pattern = "_essmap\\.txt$")
  
  if (purrr::is_empty(ins)) {
    message("No GSPA results at ", fml_nm)
    return(NULL)
  }
  
  if (is.null(df2)) {
    ins <- ins %>% 
      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} %>% 
      {if (scale_log2r) .[grepl("_GSPA_Z", .)] else .[grepl("_GSPA_N", .)]}
  
    if (rlang::is_empty(ins)) {
      stop("No GSPA inputs correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
           call. = FALSE)
    }
  } else {
    local({
      non_exists <- df2 %>% .[! . %in% ins]
      if (!purrr::is_empty(non_exists)) {
        stop("Missing _essmap file(s): ", purrr::reduce(non_exists, paste, sep = ", "), 
             call. = FALSE)
      }
    })
    
    if (purrr::is_empty(df2)) stop("Input file(s) not found.", call. = FALSE)
    ins <- ins %>% .[. %in% df2]    
  }

  meta_ins <- list.files(path = file.path(filepath, fml_nm), pattern = "_essmeta\\.txt$") 
  meta_ins <- local({
    required <- ins %>% gsub("_essmap\\.txt$", "_essmeta.txt", .)
    meta_ins <- meta_ins %>% .[. %in% required]
  })
  
  purrr::walk2(ins, meta_ins, byfile_gspahm, fml_nm, filepath, filename, scale_log2r, impute_na, ...)
}


#' A helper function for gspaHM by input file names
#' 
#' @param ess_in A list of file names for input data
#' @param meta_in A list of file names for input metadata
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams fml_gspa
#' @import purrr dplyr rlang pheatmap 
#' @importFrom magrittr %>% %T>% %$% %<>% 
byfile_gspahm <- function (ess_in, meta_in, fml_nm, filepath, filename, scale_log2r, impute_na, ...) {
  custom_prefix <- gsub("(.*_{0,1})Protein_GSPA.*", "\\1", ess_in)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) # %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  filename <- paste0(custom_prefix, fn_prefix, ".", fn_suffix)

  dots <- rlang::enexprs(...)
  if (!purrr::is_empty(dots)) {
    if (any(grepl("^filter_", names(dots)))) {
      stop("Primary `filter_` depreciated; use secondary `filter2_`.")
    }      
  }
  filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
  arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
  select2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select2_", names(.))]
  dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots, select2_dots)]
  
  all_by_greedy <- tryCatch(readr::read_tsv(file.path(filepath, fml_nm, ess_in), 
                                            col_types = cols(ess_term = col_factor(), 
                                                             size = col_double(), 
                                                             ess_size = col_double(), 
                                                             fraction = col_double(),
                                                             distance = col_double(), 
                                                             idx = col_double(),
                                                             ess_idx = col_double())), 
                            error = function(e) NA) %>% 
    filters_in_call(!!!filter2_dots) %>% 
    arrangers_in_call(!!!arrange2_dots)
  
  rm(filter2_dots, arrange2_dots, select2_dots)
  
  if (nrow(all_by_greedy) == 0) stop("No GSPA terms available after data filtration.")
  
  if (max(all_by_greedy$distance) == 0) {
    warning("Identical, all-zero distance detected; 
         try lower the criteria in data filtrations or 
         rerun `prnGSPA` at more relaxed `gspval_cutoff` threshold.")
    
    return(NULL)
  }

  if (!is.null(dim(all_by_greedy))) {
    message(paste("File loaded:", ess_in, "at", fml_nm))
  } else {
    stop("Essential GSPA not found.")
  }
  
  ess_vs_all <- all_by_greedy %>% 
    dplyr::select(-which(names(.) %in% c("size", "ess_size", "idx", "ess_idx", "distance"))) %>% 
    tidyr::spread(term, fraction) %>% 
    dplyr::mutate_at(vars(which(names(.) != "ess_term")), ~ {1 - .x}) %>% 
    tibble::column_to_rownames("ess_term")
  
  d_row <- dist(ess_vs_all)
  d_col <- dist(t(ess_vs_all))
  max_d_row <- max(d_row, na.rm = TRUE)
  max_d_col <- max(d_col, na.rm = TRUE)
  
  d_row[is.na(d_row)] <- 1.2 * max_d_row
  d_col[is.na(d_col)] <- 1.2 * max_d_col
  
  if (is.infinite(max_d_row)) {
    cluster_rows <- FALSE
  } else {
    cluster_rows <- hclust(d_row)
  }
  
  if (is.infinite(max_d_col)) {
    cluster_cols <- FALSE
  } else {
    cluster_cols <- hclust(d_col)
  }
  
  if (is.null(dots$color)) {
    mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  } else {
    mypalette <- eval(dots$color, env = caller_env())
  }
  
  if (is.null(dots$annot_cols)) {
    annot_cols <- NULL
  } else {
    annot_cols <- eval(dots$annot_cols, env = caller_env()) 
  }
  
  if (is.null(dots$annot_colnames)) {
    annot_colnames <- NULL
  } else {
    annot_colnames <- eval(dots$annot_colnames, env = caller_env()) 
  }
  
  if (is.null(dots$annot_rows)) {
    annot_rows <- NULL
  } else {
    annot_rows <- eval(dots$annot_rows, env = caller_env()) 
  }
  
  if (is.null(annot_cols)) {
    annotation_col <- NA
  } else {
    annotation_col <- gspa_colAnnot(annot_cols = annot_cols, 
                                    df = all_by_greedy, 
                                    sample_ids = all_by_greedy$term)
  }
  
  if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
    colnames(annotation_col) <- annot_colnames
  }
  
  if (is.null(annot_rows)) {
    annotation_row <- NA
  } else {
    ess_meta <- tryCatch(readr::read_tsv(file.path(filepath, fml_nm, meta_in), 
                                         col_types = cols(term = col_factor())), 
                         error = function(e) NA)
    
    annot_rows <- annot_rows %>% .[. %in% names(ess_meta)]
    
    if (is_empty(annot_rows)) {
      annotation_row <- NA
    } else {
      annotation_row <- ess_meta %>% 
        dplyr::rename(ess_term = term) %>% 
        dplyr::select("ess_term", annot_rows) %>% 
        dplyr::mutate(ess_term = factor(ess_term, levels = levels(all_by_greedy$ess_term))) %>% 
        tibble::column_to_rownames("ess_term")
    }
  }
  
  if (is.null(dots$annotation_colors)) {
    annotation_colors <- setHMColor(annotation_col)
  } else if (is.na(dots$annotation_colors)) {
    annotation_colors <- NA
  } else {
    annotation_colors <- eval(dots$annotation_colors, env = caller_env())
  }
  
  nrow <- nrow(ess_vs_all)
  ncol <- ncol(ess_vs_all)
  max_width <- 44
  max_height <- 44
  
  if (is.null(dots$width)) {
    if (ncol <= 150) {
      fontsize <- 12
      fontsize_col <- cellwidth <- 12
      show_colnames <- TRUE
      width <- pmax(ncol/1.2, 8)
    } else {
      fontsize <- fontsize_col <- 1
      cellwidth <- NA
      show_colnames <- FALSE
      width <- max_width
    }
  } else {
    width <- dots$width
    if (is.null(dots$fontsize)) fontsize <- 1 else fontsize <- dots$fontsize
    if (is.null(dots$fontsize_col)) fontsize_col <- 1 else fontsize_col <- dots$fontsize_col
    if (is.null(dots$cellwidth)) cellwidth <- NA else cellwidth <- dots$cellwidth
    if (is.null(dots$show_colnames)) show_colnames <- FALSE else show_colnames <- dots$show_colnames
  }
  
  if (is.null(dots$height)) {
    if (nrow <= 150) {
      fontsize <- 12 
      fontsize_row <- cellheight <- 12
      show_rownames <- TRUE
      height <- pmax(nrow/1.2, 8)
    } else {
      fontsize <- fontsize_row <- 1
      cellheight <- NA
      show_rownames <- FALSE
      height <- max_height
    }
  } else {
    height <- dots$height
    if (is.null(dots$fontsize)) fontsize <- 1 else fontsize <- dots$fontsize
    if (is.null(dots$fontsize_row)) fontsize_row <- 1 else fontsize_row <- dots$fontsize_row
    if (is.null(dots$cellheight)) cellheight <- NA else cellheight <- dots$cellheight
    if (is.null(dots$show_rownames)) show_rownames <- FALSE else show_rownames <- dots$show_rownames
  }
  
  if ((!is.na(width)) & (width > max_width)) {
    warning("The plot width is set to ", max_width, call. = FALSE)
    width <- max_width
    height <- pmin(max_height, width * nrow / ncol)
  } 
  
  if ((!is.na(height)) & (height > max_height)) {
    warning("The plot height is set to ", max_height, call. = FALSE)
    height <- max_height
    width <- pmin(max_width, height * ncol / nrow)
  }
  
  dots <- dots %>% 
    .[!names(.) %in% c("annot_cols", "annot_colnames", "annot_rows", 
                       "mat", "filename", "annotation_col", "annotation_row", 
                       "color", "annotation_colors", "breaks", 
                       "cluster_rows", "cluster_cols", 
                       "width", "height", "fontsize_row", "fontsize_col", 
                       "cellheight", "cellwidth")]
  
  ph <- my_pheatmap(
    mat = ess_vs_all,
    filename = file.path(filepath, fml_nm, filename),
    annotation_col = annotation_col,
    annotation_row = annotation_row, 
    color = mypalette,
    annotation_colors = annotation_colors, 
    breaks = NA, # not used
    cluster_rows = cluster_rows, 
    cluster_cols = cluster_cols,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    cellheight = cellheight,
    cellwidth = cellwidth,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    width = width, 
    height = height, 
    !!!dots,
  )
  
  # networks
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop("\n====================================================================", 
         "\nNeed install package \"networkD3\" needed for this function to work.",
         "\n====================================================================",
         call. = FALSE)
  }
  
  cluster <- data.frame(cluster = cutree(ph$tree_col, h = max_d_col)) %>% 
    tibble::rownames_to_column("term")
  
  all_by_greedy <- local({
    essterms <- all_by_greedy$ess_term %>% unique() %>% as.character()
    terms <- all_by_greedy$term %>% unique() %>% .[! . %in% essterms] %>% as.character()
    universe <- union(essterms, terms)
    
    all_by_greedy <- all_by_greedy %>% 
      dplyr::select(-idx, -ess_idx) %>% 
      dplyr::mutate(term = factor(term, levels = universe)) %>% 
      dplyr::mutate(ess_term = factor(ess_term, levels = universe)) %>%
      dplyr::mutate(source = as.numeric(term), target = as.numeric(ess_term))
    
    min_target <- min(all_by_greedy$target, na.rm = TRUE)
    
    all_by_greedy %>% 
      dplyr::mutate(source = source - min_target, target = target - min_target) %>% 
      dplyr::mutate(term = as.character(term)) %>% 
      dplyr::left_join(cluster, by = "term") %>% 
      dplyr::mutate(term = factor(term, levels(ess_term)))
  })
  
  my_nodes <- local({
    my_nodes <- all_by_greedy %>% 
      dplyr::arrange(-size) %>% 
      dplyr::select(c("term", "cluster", "size")) %>% 
      dplyr::filter(!duplicated(term)) %>% 
      dplyr::arrange(cluster) %>% 
      dplyr::arrange(term)
    
    # an `ess_term` may be missing from `term` after data row filtration
    temp_nodes <- all_by_greedy %>% 
      dplyr::select(c("ess_term", "cluster", "size")) %>% 
      dplyr::filter(!duplicated(ess_term)) %>% 
      dplyr::filter(! ess_term %in% my_nodes$term) %>% 
      dplyr::arrange(cluster) %>% 
      dplyr::arrange(ess_term) %>% 
      dplyr::rename(term = ess_term)
    
    my_nodes <- bind_rows(my_nodes, temp_nodes) %>% 
      dplyr::arrange(cluster) %>% 
      dplyr::arrange(term)
  }) %>% 
    data.frame(check.names = FALSE)
  
  my_links <- all_by_greedy %>% 
    dplyr::arrange(term) %>% 
    dplyr::select(source, target, fraction) %>% 
    dplyr::mutate(fraction = fraction * 10) %>% 
    dplyr::distinct() %>% 
    data.frame(check.names = FALSE)
  
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  networkD3::forceNetwork(Links = my_links, Nodes = my_nodes, Source = "source",
                          Target = "target", Value = "fraction", NodeID = "term", Nodesize = "size", 
                          Group = "cluster", opacity = 0.8, zoom = TRUE) %>% 
    networkD3::saveNetwork(file = file.path(filepath, fml_nm, paste0(fn_prefix, ".html")))
}


#' Sets up the column annotation in heat maps
#' 
#' @param sample_ids A character vecotr containing the sample IDs for an ascribing analysis. 
#' @inheritParams prnEucDist
#' @import dplyr rlang
#' @importFrom magrittr %>% %T>% %$% %<>% 
gspa_colAnnot <- function (annot_cols = NULL, df, sample_ids, annot_colnames = NULL) {
  if (is.null(annot_cols)) return(NA)
  
  exists <- annot_cols %in% names(df)
  
  if (sum(!exists) > 0) {
    warning(paste0("Column '", annot_cols[!exists], "'",
                   " not found in GSPA inputs and will be skipped."))
    annot_cols <- annot_cols[exists]
  }
  
  if (length(annot_cols) == 0) return(NA)
  
  x <- df %>%
    dplyr::filter(term %in% sample_ids, !duplicated(term)) %>%
    dplyr::select(annot_cols, term) %>%
    dplyr::select(which(not_all_NA(.))) %>%
    data.frame(check.names = FALSE) %>%
    `rownames<-`(.[["term"]])
  
  if (any(duplicated(x[["term"]]))) stop("Duplicated terms found\n")
  
  if (!"term" %in% annot_cols) x <- x %>% dplyr::select(-term)
  
  if (ncol(x) == 0) return(NA)
  
  return(x)
}


