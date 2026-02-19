#'GSPA of peptide data
#'
#'\code{pepGSPA} performs the analysis of Gene Set Probability Asymmetricity
#'(GSPA) against peptide \code{log2FC} data.
#'
#'@rdname prnGSPA
#'
#'@import purrr dplyr
#'@export
pepGSPA <- function (gset_nms = c("go_sets", "c2_msig", "kinsub"), 
                     method = "mean", scale_log2r = TRUE, 
                     complete_cases = FALSE, impute_na = FALSE, 
                     pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     gsscore_cutoff = 10.0, gslogFC_cutoff = log2(1.2), 
                     gspval_cutoff = 1E-2, 
                     min_size = 10L, max_size = .Machine$integer.max, 
                     min_delta = 5L, min_greedy_size = 1L, 
                     use_adjP = FALSE, 
                     fml_nms = NULL, df = NULL, id_gspa = "entrez", 
                     filepath = NULL, filename = NULL, ...) 
{
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) |>
        c(rlang::enexprs(...)) |>
        save_call(paste0("anal", "_pepGSPA"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) |>
        c(rlang::enexprs(...)) |>
        save_call(paste0("anal", "_prnGSPA"))
    }, 
    add = TRUE
  )  
  
  check_dots(c("id", "anal_type", "var_cutoff"), ...)
  check_gset_nms(gset_nms)
  
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Peptide/GSPA/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  id <- tryCatch(
    match_call_arg(normPSM, group_psm_by), 
    error = function(e) NA)
  
  if (is.na(id)) {
    id <- tryCatch(
      match_call_arg(makePepDIANN, group_psm_by), 
      error = function(e) NA)
  }
  
  if (is.na(id)) {
    id <- "pep_seq_mod"
  }

  stopifnot(
    rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), length(id) == 1L)
  
  id_gspa <- rlang::as_string(rlang::enexpr(id_gspa))

  scale_log2r <- 
    match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  method <- rlang::enexpr(method)
  method <- if (rlang::is_call(method))
    eval(method, envir = rlang::caller_env())
  else
    rlang::as_string(method)

  stopifnot(all(method %in% c("mean", "limma")))
  
  dots <- rlang::enexprs(...)
  
  dots <- concat_fml_dots(
    fmls = dots %>% .[grepl("^\\s*~", .)], 
    fml_nms = fml_nms, 
    dots = dots %>% .[!grepl("^\\s*~", .)], 
    anal_type = "GSPA"
  )
  
  reload_expts()
  
  info_anal(df = !!df, 
            id = !!id, 
            id_gspa = id_gspa, 
            filepath = !!filepath, 
            filename = !!filename, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            anal_type = "GSPA")(
              gset_nms = gset_nms, 
              var_cutoff = 1000, 
              pval_cutoff = pval_cutoff, 
              logFC_cutoff = logFC_cutoff, 
              gsscore_cutoff = gsscore_cutoff, 
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
#'  full file path(s), or both, to gene sets for enrichment analysis. For
#'  species among \code{"human", "mouse", "rat"}, the default of
#'  \code{c("go_sets", "c2_msig", "kinsub")} will utilize terms from gene
#'  ontology (\code{GO}), molecular signatures (\code{MSig}) and
#'  kinase-substrate network (\code{PSP Kinase-Substrate}). Custom \code{GO},
#'  \code{MSig} and other data bases at given species are also supported. See
#'  also: \code{\link{prepGO}} for the preparation of custom \code{GO};
#'  \code{\link{prepMSig}} for the preparation of custom \code{MSig}. For other
#'  custom data bases, follow the same format of list as \code{GO} or
#'  \code{MSig}.
#'@param id_gspa Character string; the ID type used in GSPA. The default is
#'  \code{entrez}.
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
#'@param gspval_cutoff Depreciated. Numeric value or vector; the cut-off in
#'  gene-set significance \code{pVal}. Only enrichment terms with \code{pVals}
#'  more significant than the threshold will be reported. The default is 0.05
#'  for all formulas matched to or specified in argument \code{fml_nms}.
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
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>%
#'
#'@example inst/extdata/examples/prnGSPA_.R
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
#'  analysis by protein significance pVals \cr \code{\link{prnGSPAMap}} for
#'  mapping GSPA to volcano plot visualization \cr \code{\link{prnGSPAHM}} for
#'  heat map and network visualization of GSPA results \cr \code{\link{prnGSVA}}
#'  for gene set variance analysis \cr \code{\link{prnGSEA}} for data
#'  preparation for online GSEA. \cr \code{\link{pepMDS}} and
#'  \code{\link{prnMDS}} for MDS visualization \cr \code{\link{pepPCA}} and
#'  \code{\link{prnPCA}} for PCA visualization \cr \code{\link{pepLDA}} and
#'  \code{\link{prnLDA}} for LDA visualization \cr \code{\link{pepHM}} and
#'  \code{\link{prnHM}} for heat map visualization \cr
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
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
#'  among RefSeq accessions, gene names and Entrez IDs \cr
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'  ontology}} \cr
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular
#'  signatures} \cr
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
prnGSPA <- function (
    gset_nms = c("go_sets", "c2_msig", "kinsub"), method = "mean", 
    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
    pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
    gsscore_cutoff = 10.0, gslogFC_cutoff = log2(1.2), 
    gspval_cutoff = 1E-2, min_size = 10L, 
    max_size = .Machine$integer.max, min_delta = 5L, min_greedy_size = 1L, 
    use_adjP = FALSE, fml_nms = NULL, df = NULL, id_gspa = "entrez", 
    filepath = NULL, filename = NULL, ...) 
{
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) |>
        c(rlang::enexprs(...)) |>
        save_call(paste0("anal", "_pepGSPA"))
    } 
    else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) |>
        c(rlang::enexprs(...)) |>
        save_call(paste0("anal", "_prnGSPA"))
    }, 
    add = TRUE
  )  
  
  check_dots(c("id", "anal_type", "var_cutoff"), ...)
  check_gset_nms(gset_nms)
  
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Protein/GSPA/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  id <- tryCatch(
    match_call_arg(normPSM, group_pep_by), 
    error = function(e) NA)
  
  if (is.na(id)) {
    id <- tryCatch(
      match_call_arg(makeProtDIANN, group_pep_by), 
      error = function(e) NA)
  }
  
  if (is.na(id)) {
    id <- "gene"
  }
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1L)

  scale_log2r <- match_prnSig_scale_log2r(
    scale_log2r = scale_log2r, impute_na = impute_na)

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  id_gspa <- rlang::as_string(rlang::enexpr(id_gspa))

  method <- rlang::enexpr(method)
  
  method <- if (rlang::is_call(method))
    eval(method, envir = rlang::caller_env())
  else
    rlang::as_string(method)

  stopifnot(all(method %in% c("mean", "limma")))
  
  dots <- rlang::enexprs(...)
  
  dots <- concat_fml_dots(
    fmls = dots[grepl("^\\s*~", dots)], 
    fml_nms = fml_nms, 
    dots = dots[!grepl("^\\s*~", dots)], 
    anal_type = "GSPA"
  )

  reload_expts()
  
  info_anal(df = !!df, 
            id = !!id, 
            id_gspa = id_gspa, 
            filepath = !!filepath, 
            filename = !!filename, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            anal_type = "GSPA")(
              gset_nms = gset_nms, 
              var_cutoff = 1000, 
              pval_cutoff = pval_cutoff, 
              logFC_cutoff = logFC_cutoff, 
              gsscore_cutoff = gsscore_cutoff, 
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
#' @import limma stringr purrr tidyr dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% is_greater_than not 
gspaTest <- function(df = NULL, id = "entrez", id_gspa = "entrez", 
                     label_scheme_sub = NULL, 
                     scale_log2r = TRUE, complete_cases = FALSE, 
                     impute_na = FALSE,filepath = NULL, filename = NULL, 
                     gset_nms = "go_sets", var_cutoff = 0.5, 
                     pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     gsscore_cutoff = 10.0, gslogFC_cutoff = log2(1), 
                     gspval_cutoff = 1E-2, 
                     min_size = 6L, max_size = .Machine$integer.max, 
                     min_delta = 5L, min_greedy_size = 1L, use_adjP = FALSE, 
                     method = "mean", anal_type = "GSPA", ...) 
{
  # `GSPA`: use `pVals` instead of `N_log2R` or `Z_log2R`; 
  # Protein[_impNA]_pVals.txt does not contain `scale_log2_r` info in the file name, 
  #   and cannot tell `pVals` are from `N_log2R` or `Z_log2R`;  
  #   therefore, `scale_log2r` will be matched to those in `prnSig()` and indicated in 
  #   the names of out files as `_N[_impNA].txt` or `_Z[_impNA].txt` 
  #     to inform user the `scale_log_r` status
  # `complete_cases` is for entrez IDs, pVals, log2FC
  # "id" for tibbling rownames
  
  stopifnot(vapply(c(pval_cutoff, 
                     logFC_cutoff, 
                     min_size, 
                     max_size, 
                     min_delta,
                     min_greedy_size), 
                   rlang::is_double, logical(1L)))
  
  if (!nrow(label_scheme_sub)) {
    stop("Empty metadata.")
  }

  if (!"entrez" %in% names(df)) {
    stop("Column \"entrez\" not found.")
  }

  if (all(is.na(df[["entrez"]]))) {
    stop("All values are NA under the column `entrez` in the input data.")
  }

  id   <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  fmls <- dots[grepl("^\\s*~", dots)]
  dots <- dots[!names(dots) %in% names(fmls)]

  lang_dots    <- dots[unlist(lapply(dots, is.language))]
  filter_dots  <- lang_dots[grepl("^filter_", names(lang_dots))]
  arrange_dots <- lang_dots[grepl("^arrange_", names(lang_dots))]
  select_dots  <- lang_dots[grepl("^select_", names(lang_dots))]
  dots         <- dots[!dots %in% c(filter_dots, arrange_dots, select_dots)]

  df <- df |>
    filters_in_call(!!!filter_dots) |>
    arrangers_in_call(!!!arrange_dots)
  
  if (!nrow(df)) {
    stop("No data available after row filtration.")
  }

  if (!length(fmls)) {
    stop("Formula(s) of contrasts not available.")
  }

  species <- unique(df$species) %>% .[!is.na(.)] %>% as.character()
  gsets   <- load_dbs(gset_nms = gset_nms, species = species)
  
  fml_nms <- names(df) %>% 
    .[grepl("pVal\\s*\\(", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique() %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls[names(fmls) %in% fml_nms]
  fml_nms <- fml_nms[fml_nms %in% names(fmls)]

  if (!length(fml_nms)) {
    stop("No formula matached; ", 
         "compare the formula name(s) with those in \"prnSig(..)\"")
  }

  # Note: wrong number of arguments to '>' at devtools building of a package
  # (https://github.com/Mouse-Imaging-Centre/RMINC/issues/226)
  col_ind <- fml_nms %>% 
    lapply(function (x) grepl(paste0("^", x, "\\."), names(df))) %>%
    `names<-`(paste0("nm_", seq_along(.))) %>% 
    dplyr::bind_cols() %>%
    rowSums() %>%
    magrittr::is_greater_than(0)
  
  if (!sum(col_ind)) {
    stop("No matching columns found for the given formula(s).")
  }

  switch(anal_type,
         GSPA = purrr::pwalk(
           list(fml = fmls, 
                fml_nm = fml_nms, 
                pval_cutoff = pval_cutoff, 
                logFC_cutoff = logFC_cutoff, 
                gsscore_cutoff = gsscore_cutoff, 
                gslogFC_cutoff = gslogFC_cutoff, 
                gspval_cutoff = gspval_cutoff, 
                min_size = min_size, 
                max_size = max_size, 
                min_delta = min_delta, 
                min_greedy_size = min_greedy_size, 
                method = method), 
         fml_gspa, 
         df = df, 
         col_ind = col_ind, 
         id = !!id, 
         id_gspa = id_gspa,
         gsets = gsets, 
         label_scheme_sub = label_scheme_sub, 
         complete_cases = complete_cases, 
         scale_log2r = scale_log2r, 
         filepath = filepath, 
         filename = filename, 
         use_adjP = use_adjP, 
         !!!filter_dots), 
         
         GSEA = purrr::pwalk(
           list(fml = fmls, 
                fml_nm = fml_nms, 
                var_cutoff = var_cutoff, 
                pval_cutoff = pval_cutoff, 
                logFC_cutoff = logFC_cutoff, 
                # gsscore_cutoff = gsscore_cutoff, 
                gspval_cutoff = gspval_cutoff, 
                gslogFC_cutoff = gslogFC_cutoff, 
                min_size = min_size, 
                max_size = max_size), 
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


#' Helper of GSPA.
#' 
#' @param fml A character string; the formula used in \link{prnSig}.
#' @param fml_nm A character string; the name of \code{fml}.
#' @param col_ind Numeric vector; the indexes of columns for the ascribed \code{fml_nm}.
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams gspaTest
#' @inheritParams gsVolcano
#' @import purrr dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
fml_gspa <- function (fml, fml_nm, 
                      pval_cutoff = .05, logFC_cutoff = log2(1.2), 
                      gsscore_cutoff = 10.0, gslogFC_cutoff = log2(1.2), 
                      gspval_cutoff = 0.01, min_size = 10L, 
                      max_size = .Machine$integer.max, min_delta = 5L, 
                      min_greedy_size = 1L, method = "mean", 
                      df, col_ind, id = "gene", id_gspa = "entrez", 
                      gsets = "go_sets", label_scheme_sub, 
                      complete_cases = FALSE, scale_log2r = TRUE, 
                      filepath, filename, use_adjP = FALSE, ...) 
{
  on.exit(
    pars <- mget(names(formals()), envir = rlang::current_env(), 
                 inherits = FALSE) %>% 
      .[! names(.) %in% c("gsets", "df", "label_scheme_sub", "col_ind")] %>% 
      c(rlang::enexprs(...)) %>% 
      save_call(paste0(fn_prefix, "@", fml_nm))
  )
  
  if (!length(gsets)) {
    stop("No gene sets available.")
  }

  dir.create(file.path(filepath, fml_nm), recursive = TRUE, 
             showWarnings = FALSE)

  id <- rlang::as_string(rlang::enexpr(id))
  fn_prefix <- tools::file_path_sans_ext(filename)
  
  # Column 'contrast' is factor and sorted by contrast_groups
  df <- prep_gspa(df, id = id_gspa, fml_nm = fml_nm, col_ind = col_ind, 
                  pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, 
                  use_adjP = use_adjP)
  contrast_groups <- attr(df, "contrast_groups")
  n_contrs <- length(contrast_groups)
  
  if (complete_cases) {
    df <- df[complete.cases(df), ]
    
    if (!nrow(df)) {
      stop("No data row aviable after complete casing.")
    }

    if (length(nms <- names(df))) {
      message("Complete cases against columns: ", paste(nms, collapse = ", "))
    }

    rm(list = "nms")
  }

  df <- df |>
    dplyr::mutate(p_val = -log10(p_val)) |>
    dplyr::mutate(valence = ifelse(log2Ratio >= 0, "pos", "neg"), 
                  valence = factor(valence, levels = c("neg", "pos"))) |>
    tidyr::complete(entrez, contrast, valence)
  
  # do not re-`arrange` df, df_sub for limma after this point
  df <- df |> 
    dplyr::arrange(entrez, contrast, valence)
  
  gsets <- gsets[!grepl("molecular_function$", names(gsets))]
  gsets <- gsets[!grepl("cellular_component$", names(gsets))]
  gsets <- gsets[!grepl("biological_process$", names(gsets))]
  lens  <- lengths(gsets)
  gsets <- gsets[lens >= min_size & lens <= max_size]

  res <- switch(method, 
                limma = purrr::map(gsets, 
                                   gspa_summary_limma, 
                                   df = df, 
                                   min_size = min_size, 
                                   min_delta = min_delta, 
                                   gsscore_cutoff = gsscore_cutoff, 
                                   gspval_cutoff = gspval_cutoff, 
                                   gslogFC_cutoff = gslogFC_cutoff), 
                mean = purrr::map(gsets, 
                                  gspa_summary_mean, 
                                  df = df, 
                                  min_size = min_size, 
                                  min_delta = min_delta, 
                                  gsscore_cutoff = gsscore_cutoff, 
                                  gspval_cutoff = gspval_cutoff, 
                                  gslogFC_cutoff = gslogFC_cutoff), 
                rlang::abort("Unknown `method`."))
  
  res <- dplyr::bind_rows(res, .id = "term")
  
  if (!nrow(res)) {
    stop("No GSPA results available.")
  }
  
  res <- res |>
    dplyr::select(dplyr::one_of(c("term", "contrast", "log2Ratio", "score")))

  uids <- unique(df$entrez)
  sig_sets <- purrr::imap_dfr(gsets, ~ tibble::tibble(id = .x, term = .y)) |>
    dplyr::filter(id %in% uids) |>
    dplyr::select(c("term", "id"))
  
  res <- sig_sets |>
    dplyr::group_by(term) %>% 
    # the number of entries matched to input `df`
    dplyr::summarise(size = dplyr::n()) |>
    dplyr::right_join(res, by = "term")
  
  res_greedy <- greedysetcover(sig_sets) %T>% 
    write.table(file.path(filepath, fml_nm, 
                          paste0(fn_prefix, "_resgreedy.txt")), 
                sep = "\t", col.names = TRUE, 
                row.names = FALSE, quote = FALSE) %>% 
    dplyr::group_by(term) %>% 
    dplyr::summarise(ess_size = dplyr::n()) %>% 
    dplyr::filter(ess_size >= min_greedy_size)
  
  tempdata <- res_greedy %>% 
    dplyr::right_join(res, by = "term") %>% 
    dplyr::mutate(is_essential = ifelse(!is.na(ess_size), TRUE, FALSE)) %>% 
    dplyr::select(term, is_essential, size, ess_size, 
                  contrast, score, log2Ratio) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, ".txt")), 
                sep = "\t", col.names = TRUE, 
                row.names = FALSE, quote = FALSE) %>% 
    dplyr::filter(is_essential) %>% 
    dplyr::filter(!duplicated(term)) %>% 
    dplyr::select(-contrast, -score, -log2Ratio) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_essmeta.txt")), 
                sep = "\t", col.names = TRUE, 
                row.names = FALSE, quote = FALSE) %>% 
    dplyr::select(term, ess_size)
  
  sig_sets <- tempdata %>% 
    dplyr::right_join(sig_sets, by = "term")

  map_essential(sig_sets) %>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_essmap.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE, 
                quote = FALSE)

  invisible(sig_sets)
}


#' checks the size of a gene set
#'
#' @param n_ignore The number of entries to be ignored on both the left and the
#'   right side of volcano to handle, e.g., one-hit-wonders on a side.
#' @param max_low_n The maximum number of entries in the lower side of a volcano
#'   plot.
#' @param p_x A trivialized p-value for replacing NA values.
#' @param FUN A summary function.
#' @inheritParams prnHist
#' @inheritParams prnGSPA
ok_min_size <- function (df, min_delta = 4L, max_low_n = 3L, gspval_cutoff = 1.0, 
                         gsscore_cutoff = 10.0, gslogFC_cutoff = log2(1.2), 
                         n_ignore = 1L, p_x = 1.0, FUN = mean) 
{
  # data "already" ordered by entrez, contrast and then valence: "neg", "pos"
  dfw <- df |>
    # dplyr::arrange(entrez, contrast, valence) |>
    # dplyr::group_by(contrast, valence) |>
    dplyr::filter(!is.na(log2Ratio), !is.na(p_val)) |>
    dplyr::mutate(n = dplyr::n())

  if (FALSE) {
    # contrast can drop and number of columns reduced
    dfw <- dfw |>
      dplyr::mutate(n = n - n_ignore) |>
      dplyr::filter(n > 0L)
  }
  
  # "unique" can change order: 
  # unique(dfw[, c("contrast", "valence", "n")]) |> dplyr::arrange(contrast, valence)
  delta_n <- dfw |>
    dplyr::summarise(n = mean(n, na.rm = TRUE)) |>
    # dplyr::slice_head(n = 1L) |> 
    # dplyr::select(c("contrast", "valence", "n")) |>
    # dplyr::ungroup() |> 
    tidyr::pivot_wider(
      names_from = contrast, 
      values_from = n, 
      values_fill = 0L) |>
    dplyr::arrange(valence) |>
    dplyr::select(-valence)

  if (nrow(delta_n) == 1L) {
    ok_min_n <- rep_len(TRUE, ncol(delta_n))
    names(ok_min_n) <- colnames(delta_n)
  } else {
    min_n <- sapply(delta_n, min)
    max_n <- sapply(delta_n, max)
    ok_min_n <- (min_n <= max_low_n) | 
      (max_n / min_n >= 4. & min_n <= 2 * max_low_n)
  }

  delta_p <- dfw |>
    dplyr::summarise(p_val = sum(p_val, na.rm = TRUE)) |>
    tidyr::pivot_wider(
      names_from = contrast, 
      values_from = p_val, 
      values_fill = p_x) |>
    dplyr::arrange(valence) |>
    dplyr::select(-valence)

  delta_fc <- dfw |>
    dplyr::summarise(log2Ratio = median(log2Ratio, na.rm = TRUE)) |>
    tidyr::pivot_wider(
      names_from = contrast, 
      values_from = log2Ratio, 
      values_fill = 0.0) |>
    dplyr::arrange(valence) |>
    dplyr::select(-valence)

  # delta_sc <- abs(delta_p  * log2(delta_n + 1L)/2) # with mean stat
  delta_sc <- delta_p  / log2(delta_n + 2L) # with sum stat

  if (nrow(delta_p) == 2L) {
    delta_p1  <- delta_p[2, ]  - delta_p[1, ]
  } else {
    # all pos or all neg
    delta_p1  <- delta_p[1, ]
  }
  
  if (nrow(delta_fc) == 2L) {
    delta_fc1 <- delta_fc[2, ] - delta_fc[1, ]
  } else {
    delta_fc1 <- delta_fc[1, ]
  }
  
  if (nrow(delta_n) == 2L) {
    delta_n1  <- delta_n[2, ]  - delta_n[1, ]
  } else {
    delta_n1  <- delta_n[1, ]
  }
  
  if (nrow(delta_sc) == 2L) {
    delta_sc1 <- delta_sc[2, ] - delta_sc[1, ]
  } else {
    delta_sc1 <- delta_sc[1, ]
  }

  passes <- sign(delta_sc1) == sign(delta_n1) & 
    # (abs(delta_p1) >= -log10(gspval_cutoff) | abs(delta_sc1) >= gsscore_cutoff) & 
    abs(delta_sc1) >= gsscore_cutoff & 
    ok_min_n & 
    abs(delta_n1) >= min_delta & 
    abs(delta_fc1) >= gslogFC_cutoff
  
  ## pad missing contrast(s)
  levs <- levels(df$contrast)
  if (length(passes) < length(levs)) {
    missing_nms <- levs[!levs %in% colnames(passes)]
    
    if (length(missing_nms)) {
      for (i in seq_along(missing_nms)) {
        nmi <- missing_nms[i]
        passes    <- dplyr::bind_cols(!!nmi := FALSE, passes)
        delta_p1  <- dplyr::bind_cols(!!nmi := -log10(p_x), delta_p1)
        delta_fc1 <- dplyr::bind_cols(!!nmi := 0.0, delta_fc1)
        delta_n1  <- dplyr::bind_cols(!!nmi := 0L, delta_n1)
        delta_sc1 <- dplyr::bind_cols(!!nmi := 0.0, delta_sc1)
      }
    }
    
    passes   <- passes[, levs, drop = FALSE]
    delta_p1  <- delta_p1[, levs, drop = FALSE]
    delta_fc1 <- delta_fc1[, levs, drop = FALSE]
    delta_n1  <- delta_n1[, levs, drop = FALSE]
    delta_sc1 <- delta_sc1[, levs, drop = FALSE]
  }
  
  # Outputs
  ans <- dplyr::bind_cols(
    contrast = colnames(passes), 
    pass     = as.vector(t(passes)), 
    pVal     = abs(as.vector(t(delta_p1))), 
    log2Ratio = as.vector(t(delta_fc1)), 
    # n  = abs(as.vector(t(delta_n1))), 
    score = abs(as.vector(t(delta_sc1))), )

  ans <- ans |>
    dplyr::mutate(contrast = factor(contrast, levels = levels(df$contrast))) |>
    dplyr::arrange(contrast)

  ans |> dplyr::filter(pass)
}


#' gspa pVal calculations using geomean
#' 
#' @param gset Character string; a gene set indicated by \code{gset_nms}.
#' @param p_x A trivalized p-value for replacing NA values.
#' @inheritParams prnHist
#' @inheritParams prnGSPA
gspa_summary_mean <- function(gset, df, min_size = 10, min_delta = 4L, 
                              gsscore_cutoff = 10., gspval_cutoff = 0.05, 
                              gslogFC_cutoff = log2(1.0), p_x = 1.0) 
{
  df <- df |>
    dplyr::filter(entrez %in% gset) |>
    dplyr::group_by(contrast, valence)
  
  if (length(unique(df$entrez)) < min_size) {
    return(NULL)
  }

  ok_min_size(
    df = df, min_delta = min_delta, gspval_cutoff = gspval_cutoff, 
    gsscore_cutoff = gsscore_cutoff, gslogFC_cutoff = gslogFC_cutoff)
} 


#' GSPA pVal calculations using limma
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams gspa_summary_mean
gspa_summary_limma <- function(gset, df, min_size = 10, min_delta = 4, 
                               gsscore_cutoff = 10., gspval_cutoff = 0.05, 
                               gslogFC_cutoff = 1.2) 
{
  df <- df %>% 
    dplyr::filter(entrez %in% gset)
  
  if (!nrow(df))
    stop("No entrez IDs matched to the provided gene sets.")

  if (length(unique(df$entrez)) < min_size) 
    return(NULL)
  
  lm_gspa(df, min_delta, gspval_cutoff, gslogFC_cutoff)
}  


#' significance tests of pVals between the up and the down groups
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
lm_gspa <- function(df, min_delta, gsscore_cutoff = 10.0, 
                    gspval_cutoff = .01, gslogFC_cutoff = 1.2) 
{
  delta_fc <- df |>
    dplyr::group_by(contrast, valence) |>
    dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) |>
    dplyr::mutate(log2Ratio = replace(log2Ratio, is.nan(log2Ratio), 0)) |>
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = log2Ratio
    ) |>
    dplyr::select(-valence)

  nms <- colnames(delta_fc)
  
  delta_fc <- delta_fc %>% 
    `colnames<-`(paste0("log2Ratio (", colnames(.), ")")) %>% 
    abs()
  
  if (nrow(delta_fc) == 2L) 
    delta_fc <- delta_fc[2, ] - delta_fc[1, ]
  
  ok_delta_n <- ok_min_size(df = df, 
                            min_delta = min_delta, 
                            gspval_cutoff = gspval_cutoff, 
                            gsscore_cutoff = gsscore_cutoff, 
                            gslogFC_cutoff = gslogFC_cutoff)

  data <- df %>% 
    dplyr::group_by(contrast, valence) %>% 
    dplyr::select(-log2Ratio) %>% 
    tidyr::unite(sample, entrez, contrast, valence, sep = ".", 
                 remove = TRUE) %>% 
    tidyr::spread(key = sample, value = p_val) 
  
  label_scheme_gspa <- tibble::tibble(Sample_ID = names(data), 
                                      Group = Sample_ID) %>% 
    dplyr::mutate(Group = gsub("^\\d+\\.", "", Group)) %>% 
    dplyr::mutate(Group = factor(Group))
  
  design <- model.matrix(~0+label_scheme_gspa$Group) %>% 
    `colnames<-`(levels(label_scheme_gspa$Group)) %>% 
    data.frame(check.names = TRUE)

  neg_cols <- colnames(design)[as.logical(seq_along(design) %% 2)]
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

  delta_p[] <- purrr::map2(as.list(delta_p), as.list(ok_delta_n), 
                           ~ {if (.y) .x else 1}) 

  invisible(dplyr::bind_cols(delta_p, delta_fc))
}


#' Helper of GSPA
#'
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams fml_gspa
#' 
#' @import purrr dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom readr read_tsv
prep_gspa <- function(df = NULL, id = NULL, fml_nm = NULL, 
                      col_ind = 0L, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                      use_adjP = FALSE) 
{
  df <- df %>%
    dplyr::select(matches(paste0("^", fml_nm, "\\."))) %>%
    `colnames<-`(gsub(paste0("^", fml_nm, "\\."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .) %>% 
    rm_pval_whitespace() %>% 
    dplyr::select(id, grep("^pVal|^adjP|^log2Ratio", names(.))) %>% 
    dplyr::mutate(!!id := as.character(.[[id]]))
  
  if (use_adjP) {
    df <- df |>
      dplyr::select(-dplyr::matches("pVal")) %>% 
      `names<-`(gsub("adjP", "pVal", names(.)))
  } 
  else {
    df <- df |>
      dplyr::select(-dplyr::matches("adjP"))
  }
  
  contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
    gsub("^log2Ratio\\s+\\(|\\)$", "", .)
  
  pvals <- df |>
    dplyr::select(-dplyr::matches("^log2Ratio\\s+\\(")) |>
    dplyr::rename_with(
      function (x) gsub("^pVal\\s+\\((.*)\\)$", "\\1", x)
    ) |>
    tidyr::pivot_longer(
      cols = -id,
      names_to = "contrast",
      values_to = "p_val"
    )
  
  log2rs <- df |>
    dplyr::select(-dplyr::matches("^pVal\\s+\\(")) |>
    dplyr::rename_with(
      function (x) gsub("^log2Ratio\\s+\\((.*)\\)$", "\\1", x)
    ) |>
    tidyr::pivot_longer(
      cols = -id,
      names_to = "contrast",
      values_to = "log2Ratio"
    ) |>
    dplyr::select(-id, -contrast)

  df <- dplyr::bind_cols(pvals, log2rs) |>
    dplyr::filter(p_val <= pval_cutoff, 
                  abs(log2Ratio) >= logFC_cutoff) |> # NA will be removed
    dplyr::filter(!is.na(id)) |>
    dplyr::mutate(contrast = factor(contrast, levels = contrast_groups)) |>
    dplyr::arrange(contrast)
  
  # e.g., id == "pep_seq_mod" for phosphoproteomics
  if (id != "entrez") {
    df <- df |>
      dplyr::rename(entrez = !!id)
  }
  
  attr(df, "contrast_groups") <- contrast_groups

  df
}


#' Helper for mapping between gene sets and essential gene sets
#'
#' @param sig_sets A data frame containing the gene sets that are significant
#'   under given criteria.
#' @import purrr dplyr  
#' @importFrom magrittr %>% %T>% %$% %<>% 
map_essential <- function (sig_sets) 
{
  ess_terms <- sig_sets %>% 
    dplyr::filter(!is.na(ess_size), !duplicated(term)) %>% 
    dplyr::select(term) %>% 
    unlist

  sig_greedy_sets <- sig_sets %>% 
    dplyr::filter(term %in% ess_terms)

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


#' Heat map visualization of GSPA results
#'
#' \code{prnGSPAHM} visualizes distance heat maps and networks between essential
#' and all gene sets.
#'
#' The list of gene sets and the associative quality metrics of \code{size} and
#' \code{ess_size} are assessed after data filtration with the criteria
#' specified by arguments \code{pval_cutoff} and \code{logFC_cutoff}, as well as
#' optional varargs of \code{filter_}.
#'
#' @section \code{Protein_GSPA_[...].txt}:
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
#' @section \code{Protein_GSPA_[...]essmap.txt}:
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
#' @inheritParams plot_prnTrend
#' @inheritParams  prnEucDist
#' @param impute_na Logical; at TRUE, input files with \code{_impNA[...].txt} in
#'   name will be loaded. Otherwise, files without \code{_impNA} in name will be
#'   taken. An error will be thrown if no files are matched under given
#'   conditions. The default is FALSE.
#' @param fml_nms Character string or vector; the formula name(s). By default,
#'   the formula(s) will match those used in \code{\link{pepSig}} or
#'   \code{\link{prnSig}}.
#' @param annot_cols A character vector of column keys that can be found in
#'   \code{_essmap.txt}. The values under the selected keys will be used to
#'   color-code enrichment terms on the top of heat maps. The default is NULL
#'   without column annotation.
#' @param annot_rows A character vector of column keys that can be found from
#'   \code{_essmeta.txt} . The values under the selected keys will be used to
#'   color-code essential terms on the side of heat maps. The default is NULL
#'   without row annotation.
#' @param ... \code{filter2_}: Variable argument statements for the row
#'   filtration against data in secondary file(s) of \code{_essmap.txt}. Each
#'   statement contains to a list of logical expression(s). The \code{lhs} needs
#'   to start with \code{filter2_}. The logical condition(s) at the \code{rhs}
#'   needs to be enclosed in \code{exprs} with round parenthesis. For example,
#'   \code{distance} is a column key in \code{Protein_GSPA_Z_essmap.txt}. The
#'   statement \code{filter2_ = exprs(distance <= .95),} will remove entries
#'   with \code{distance > 0.95}. See also \code{\link{normPSM}} for the format
#'   of \code{filter_} statements against primary data. \cr \cr
#'   \code{arrange2_}: Variable argument statements for the row ordering against
#'   data in secondary file(s) of \code{_essmap.txt}. The \code{lhs} needs to
#'   start with \code{arrange2_}. The expression(s) at the \code{rhs} needs to
#'   be enclosed in \code{exprs} with round parenthesis. For example,
#'   \code{distance} and \code{size} are column keys in
#'   \code{Protein_GSPA_Z_essmap.txt}. The statement \code{arrange2_ =
#'   exprs(distance, size),} will order entries by \code{distance}, then by
#'   \code{size}. See also \code{\link{prnHM}} for the format of \code{arrange_}
#'   statements against primary data. \cr \cr Additional arguments for
#'   \code{\link[pheatmap]{pheatmap}}, i.e., \code{fontsize }... \cr \cr Note
#'   arguments disabled from \code{pheatmap}: \cr \code{annotation_col}; instead
#'   use keys indicated in \code{annot_cols} \cr \code{annotation_row}; instead
#'   use keys indicated in \code{annot_rows}
#' @import purrr
#' 
#' @example inst/extdata/examples/prnGSPAHM_.R
#'
#' @seealso 
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
#'  \code{\link{prnGSPAMap}} for mapping GSPA to volcano plot visualization \cr 
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
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'  
#' @export
prnGSPAHM <- function (scale_log2r = TRUE, complete_cases = FALSE, 
                       impute_na = FALSE, fml_nms = NULL, 
                       annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                       df2 = NULL, filename = NULL, ...) 
{
  old_opts <- options()
  on.exit(options(old_opts), add = TRUE)
  options(warn = 1)
  
  check_dots(c("id", "anal_type", "df", "filepath"), ...)
  
  dat_dir <- get_gl_dat_dir()
  
  dir.create(file.path(dat_dir, "Protein/GSPA/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  # id <- match_call_arg(normPSM, group_pep_by)
  id <- tryCatch(
    match_call_arg(normPSM, group_pep_by), 
    error = function(e) NA)
  
  if (is.na(id)) {
    id <- tryCatch(
      match_call_arg(makeProtDIANN, group_pep_by), 
      error = function(e) NA)
  }
  
  if (is.na(id)) {
    id <- "gene"
  }

  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1L)

  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, 
                                          impute_na = impute_na)

  df2 <- rlang::enexpr(df2)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  annot_rows <- rlang::enexpr(annot_rows) 
  
  dots <- rlang::enexprs(...)
  
  dots <- concat_fml_dots(
    fmls = dots %>% .[grepl("^\\s*~", .)], 
    fml_nms = fml_nms, 
    dots = dots %>% .[!grepl("^\\s*~", .)], 
    anal_type = "GSPA"
  )

  reload_expts()

  info_anal(id = !!id, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            df = NULL, 
            df2 = !!df2, 
            filepath = NULL, 
            filename = !!filename, 
            anal_type = "GSPA_hm")(annot_cols = !!annot_cols, 
                                   annot_colnames = !!annot_colnames, 
                                   annot_rows = !!annot_rows, 
                                   !!!dots)
}


#' Plots distance heat map of GSPA
#' 
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams prnGSPAMap
gspaHM <- function(scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                   df2, filepath, filename, ...) 
{
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  if (!length(fmls))
    stop("No formula(s) of contrasts available.")
  
  fml_nms <- names(fmls)
  
  if (length(fml_nms)) {
    purrr::walk(fml_nms, byfml_gspahm, df2, filepath, filename, 
                scale_log2r, impute_na, !!!dots)
  }
}


#' A helper function for prnGSPAMap by formula names
#' 
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams prnGSPAMap
#' @inheritParams fml_gspa
#' @import purrr dplyr pheatmap 
#' @importFrom magrittr %>% %T>% %$% %<>% 
byfml_gspahm <- function (fml_nm, df2, filepath, filename, scale_log2r = TRUE, 
                          impute_na = FALSE, ...) 
{
  ins <- list.files(path = file.path(filepath, fml_nm), 
                    pattern = "GSPA_[ONZ]_essmap\\.txt$")
  
  if (!length(ins)) {
    message("No GSPA results at ", fml_nm)
    return(NULL)
  }
  
  if (is.null(df2)) {
    ins <- ins %>% 
      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} 
    
    ins <- if (is.na(scale_log2r))
      ins[grepl("_GSPA_O_", ins)]
    else if (scale_log2r)
      ins[grepl("_GSPA_Z_", ins)]
    else
      ins[grepl("_GSPA_N_", ins)]

    if (!length(ins)) {
      stop("No GSPA inputs correspond to impute_na = ", impute_na, 
           ", scale_log2r = ", scale_log2r)
    }
  } 
  else {
    local({
      non_exists <- df2 %>% .[! . %in% ins]
      
      if (length(non_exists))
        stop("Missing _essmap file(s): ", paste(non_exists, collapse = ", "))
    })
    
    if (!length(df2)) 
      stop("Input file(s) not found.")
    
    ins <- ins %>% .[. %in% df2]    
  }

  meta_ins <- list.files(path = file.path(filepath, fml_nm), 
                         pattern = "_essmeta\\.txt$")
  
  meta_ins <- local({
    required <- gsub("_essmap\\.txt$", "_essmeta.txt", ins)
    meta_ins %>% .[. %in% required]
  })
  
  purrr::walk2(ins, meta_ins, byfile_gspahm, fml_nm, filepath, 
               filename, scale_log2r, impute_na, ...)
}


#' A helper function for prnGSPAMap by input file names
#' 
#' @param ess_in A list of file names for input data
#' @param meta_in A list of file names for input metadata
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams fml_gspa
#' @import purrr dplyr pheatmap 
#' @importFrom magrittr %>% %T>% %$% %<>% 
byfile_gspahm <- function (ess_in, meta_in, fml_nm, filepath, filename, 
                           scale_log2r, impute_na, ...) 
{
  custom_prefix <- gsub("(.*_{0,1})Protein_GSPA.*", "\\1", ess_in)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) # %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  filename <- paste0(custom_prefix, fn_prefix, ".", fn_suffix)

  dots <- rlang::enexprs(...)
  
  if (length(dots)) {
    if (any(grepl("^filter_", names(dots))))
      stop("Primary `filter_` depreciated; use secondary `filter2_`.")
  }
  
  filter2_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter2_", names(.))]
  
  arrange2_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange2_", names(.))]
  
  select2_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^select2_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter2_dots, arrange2_dots, select2_dots)]
  
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
  
  rm(list = c("filter2_dots", "arrange2_dots", "select2_dots"))

  if (!nrow(all_by_greedy))
    stop("No GSPA terms available after data filtration.")

  if (max(all_by_greedy$distance) == 0) {
    warning("Identical, all-zero distance detected; ", 
            "try lower the criteria in data filtrations or ", 
            "rerun `prnGSPA` at more relaxed `gspval_cutoff` threshold.", 
            call. = FALSE)
    
    return(NULL)
  }

  if (!is.null(dim(all_by_greedy)))
    message(paste("File loaded:", ess_in, "at", fml_nm))
  else
    stop("Essential GSPA not found.", call. = FALSE)

  ess_vs_all <- all_by_greedy %>% 
    dplyr::select(-which(names(.) %in% c("size", "ess_size", 
                                         "idx", "ess_idx", "distance"))) %>% 
    tidyr::spread(term, fraction) %>% 
    dplyr::mutate_at(vars(which(names(.) != "ess_term")), ~ {1 - .x}) %>% 
    `rownames<-`(NULL) %>% 
    tibble::column_to_rownames("ess_term")
  
  d_row <- dist(ess_vs_all)
  d_col <- dist(t(ess_vs_all))
  max_d_row <- max(d_row, na.rm = TRUE)
  max_d_col <- max(d_col, na.rm = TRUE)
  
  d_row[is.na(d_row)] <- 1.2 * max_d_row
  d_col[is.na(d_col)] <- 1.2 * max_d_col
  
  cluster_rows <- if (is.infinite(max_d_row)) FALSE else hclust(d_row)
  cluster_cols <- if (is.infinite(max_d_col)) FALSE else hclust(d_col)
  
  mypalette <- if (is.null(dots$color))
    colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  else
    eval(dots$color, envir = rlang::caller_env())

  annot_cols <- if (is.null(dots$annot_cols))
    NULL
  else
    eval(dots$annot_cols, envir = rlang::caller_env()) 

  annot_colnames <- if (is.null(dots$annot_colnames))
    NULL
  else
    eval(dots$annot_colnames, envir = rlang::caller_env()) 

  annot_rows <- if (is.null(dots$annot_rows))
    NULL
  else
    eval(dots$annot_rows, envir = rlang::caller_env()) 

  annotation_col <- if (is.null(annot_cols))
    NA
  else
    gspa_colAnnot(annot_cols = annot_cols, 
                  df = all_by_greedy, 
                  sample_ids = all_by_greedy$term)

  if (!is.null(annot_colnames) && length(annot_colnames) == length(annot_cols))
    colnames(annotation_col) <- annot_colnames

  if (is.null(annot_rows))
    annotation_row <- NA
  else {
    ess_meta <- tryCatch(readr::read_tsv(file.path(filepath, fml_nm, meta_in), 
                                         col_types = cols(term = col_factor())), 
                         error = function(e) NA)
    
    annot_rows <- annot_rows %>% .[. %in% names(ess_meta)]
    
    if (!length(annot_rows))
      annotation_row <- NA
    else {
      annotation_row <- ess_meta %>% 
        dplyr::rename(ess_term = term) %>% 
        dplyr::select("ess_term", annot_rows) %>% 
        dplyr::mutate(ess_term = 
                        factor(ess_term, levels = levels(all_by_greedy$ess_term))) %>% 
        `rownames<-`(NULL) %>% 
        tibble::column_to_rownames("ess_term")
    }
  }
  
  annotation_colors <- if (is.null(dots$annotation_colors))
    setHMColor(annotation_col)
  else if (is.na(dots$annotation_colors))
    NA
  else
    eval(dots$annotation_colors, envir = rlang::caller_env())

  nrow <- nrow(ess_vs_all)
  ncol <- ncol(ess_vs_all)
  max_width <- 44
  max_height <- 44
  
  if (is.null(dots$width)) {
    if (ncol <= 150L) {
      fontsize <- 12
      fontsize_col <- cellwidth <- 12
      show_colnames <- TRUE
      width <- pmax(ncol/1.2, 8)
    } 
    else {
      fontsize <- fontsize_col <- 1
      cellwidth <- NA
      show_colnames <- FALSE
      width <- max_width
    }
  } 
  else {
    width <- dots$width
    fontsize <- if (is.null(dots$fontsize)) 1 else dots$fontsize
    fontsize_col <- if (is.null(dots$fontsize_col)) 1 else dots$fontsize_col
    cellwidth <- if (is.null(dots$cellwidth)) NA else dots$cellwidth
    show_colnames <- if (is.null(dots$show_colnames)) FALSE else dots$show_colnames
  }
  
  if (is.null(dots$height)) {
    if (nrow <= 150) {
      fontsize <- 12 
      fontsize_row <- cellheight <- 12
      show_rownames <- TRUE
      height <- pmax(nrow/1.2, 8)
    } 
    else {
      fontsize <- fontsize_row <- 1
      cellheight <- NA
      show_rownames <- FALSE
      height <- max_height
    }
  } 
  else {
    height <- dots$height
    fontsize <- if (is.null(dots$fontsize)) 1 else dots$fontsize
    fontsize_row <- if (is.null(dots$fontsize_row))  1 else dots$fontsize_row
    cellheight <- if (is.null(dots$cellheight)) NA else dots$cellheight
    show_rownames <- if (is.null(dots$show_rownames)) FALSE else dots$show_rownames
  }
  
  if ((!is.na(width)) && (width > max_width)) {
    warning("The plot width is set to ", max_width, call. = FALSE)
    width <- max_width
    height <- min(max_height, width * nrow / ncol)
  } 
  
  if ((!is.na(height)) && (height > max_height)) {
    warning("The plot height is set to ", max_height, call. = FALSE)
    height <- max_height
    width <- min(max_width, height * ncol / nrow)
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
         "\nPackage \"networkD3\" needed for this function to work.",
         "\n====================================================================",
         call. = FALSE)
  }
  
  cluster <- data.frame(cluster = cutree(ph$tree_col, h = max_d_col)) %>% 
    tibble::rownames_to_column("term")
  
  all_by_greedy <- local({
    essterms <- unique(all_by_greedy$ess_term) %>% as.character()
    terms <- unique(all_by_greedy$term) %>% .[! . %in% essterms] %>% as.character()
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
                          Target = "target", Value = "fraction", NodeID = "term", 
                          Nodesize = "size", 
                          Group = "cluster", opacity = 0.8, zoom = TRUE) %>% 
    networkD3::saveNetwork(file = file.path(filepath, fml_nm, paste0(fn_prefix, ".html")))
}


#' Sets up the column annotation in heat maps
#' 
#' @param sample_ids A character vector containing the sample IDs for an ascribing analysis. 
#' @inheritParams prnEucDist
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
gspa_colAnnot <- function (annot_cols = NULL, df, sample_ids, annot_colnames = NULL) 
{
  if (is.null(annot_cols)) 
    return(NA)
  
  exists <- annot_cols %in% names(df)
  
  if (sum(!exists) > 0) {
    warning(paste0("Column '", annot_cols[!exists], "'",
                   " not found in GSPA inputs and will be skipped."))
    annot_cols <- annot_cols[exists]
  }
  
  if (!length(annot_cols))
    return(NA)
  
  x <- df %>%
    dplyr::filter(term %in% sample_ids, !duplicated(term)) %>%
    dplyr::select(annot_cols, term) %>%
    dplyr::select(which(not_all_NA(.))) %>%
    data.frame(check.names = FALSE) %>%
    `rownames<-`(.[["term"]])
  
  if (any(duplicated(x[["term"]]))) 
    stop("Duplicated terms found\n")
  
  if (!"term" %in% annot_cols)
    x <- dplyr::select(x, -term)

  if (!ncol(x))
    return(NA)
  
  invisible(x)
}




##########################################
# From mzion
##########################################

#' Greedy set cover.
#'
#' @param df A two-column data frame.
greedysetcover <- function (df) 
{
  ## assume: (1) two columns 
  # stopifnot(ncol(df) == 2L)
  
  len <- length(unique(df[[1]]))
  
  if (len == 1L)
    return(df)

  nms <- colnames(df)
  colnames(df) <- c("s", "a")
  
  # (2) no duplicated entries (for speeds)
  # df <- df %>%
  #   tidyr::unite(sa, s, a, sep = "@", remove = FALSE) %>%
  #   dplyr::filter(!duplicated(sa)) %>%
  #   dplyr::select(-sa)
  
  # ---
  cts <- df %>%
    dplyr::group_by(s) %>%
    dplyr::summarise(n = n()) 
  
  df <- dplyr::left_join(cts, df, by = "s") %>%
    dplyr::arrange(-n)
  
  sets <- NULL
  
  while(nrow(df)) {
    s <- df[[1, "s"]]
    sa <- df[df$s == s, c("s", "a")]
    
    sets <- rbind(sets, sa)
    
    # may consider partial sorting
    df <- df %>%
      dplyr::filter(! a %in% sa[["a"]]) %>%
      dplyr::group_by(s) %>%
      dplyr::mutate(n = n()) %>%
      dplyr::arrange(-n)
  }
  
  colnames(sets) <- nms
  
  invisible(sets)
}

