#'GSVA of protein data
#'
#'\code{prnGSVA} performs the GSVA aganist protein \code{log2FC}
#'
#'The formula(s) of contrast(s) used in \code{\link{pepSig}} will be taken by
#'default.
#'
#'
#'@inheritParams proteoSigtest
#'@inheritParams  proteoEucDist
#'@inheritParams  proteoHM
#'@inheritParams  info_anal
#'@inheritParams proteoGSPA
#'@param impute_na Logical. At the NULL default, the TRUE or FALSE will match
#'  the choice in \code{\link{pepSig}} for peptide and \code{\link{prnSig}} for
#'  protein data.
#'@param lm_method Character string indicating the linear modeling method for
#'  significance assessment of GSVA enrichment scores. The default is
#'  \code{limma}. At \code{method = lm}, the \code{lm()} in base R will be used
#'  for models without random effects and the \code{\link[lmerTest]{lmer}} will
#'  be used for models with random effects.
#'@param var_cutoff Numeric; the cut-off in the variances of protein log2FC.
#'  Entries with variances smaller than the threshold will be removed from GSVA.
#'  The default is 0.5.
#'@param pval_cutoff Numeric; the cut-off in enrichment pVals. Terms with
#'  enrichment pVals smaller than the threshold will be removed from multiple
#'  test corrections. The default is 1e-04.
#'@param logFC_cutoff Numeric; the cut-off in enrichment log2FC. Terms with
#'  absolute log2FC smaller than the threshold will be removed from multiple
#'  test corrections. The default is at log2(1.1).
#'@param filepath Use system default.
#'@param filename Use system default.
#'@param ... Arguments for \code{\link{GSVA::gsva}}
#'
#'@example inst/extdata/examples/prnGSVA_.R
#'@seealso \code{\link{load_expts}} for a reduced working example in data normalization \cr
#'
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
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPcA}} for PCA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'  
#'@import dplyr rlang ggplot2 GSVA
#'@importFrom magrittr %>%
#'@export
prnGSVA <- function (df = NULL, filepath = NULL, filename = NULL, 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = NULL, lm_method = "limma", 
                     gset_nms = "go_sets", var_cutoff = .5, pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), ...) {

  on.exit(
    if (rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(dots) %>% 
        save_call("pepGSVA")
    } else if (rlang::as_string(id) %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(dots) %>% 
        save_call("prnGSVA")
    }
    , add = TRUE
  )
  
  dir.create(file.path(dat_dir, "Protein\\GSVA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  id <- match_call_arg(normPSM, group_pep_by)
  
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	lm_method <- rlang::as_string(rlang::enexpr(lm_method))

	if (is.null(impute_na)) {
	  if (id %in% c("pep_seq", "pep_seq_mod")) {
	    impute_na <- match_call_arg(pepSig, impute_na)
	  } else if (id %in% c("prot_acc", "gene")) {
	    impute_na <- match_call_arg(prnSig, impute_na)
	  }
	}
	
	stopifnot(rlang::is_logical(scale_log2r), 
	          rlang::is_logical(impute_na), 
	          rlang::is_logical(complete_cases))

	dots <- rlang::enexprs(...)
	fmls <- dots %>% .[grepl("^\\s*~", .)]
	dots <- dots[!names(dots) %in% names(fmls)]

	if (purrr::is_empty(fmls)) {
	  fml_file <-  file.path(dat_dir, "Calls\\prnSig_formulas.rda")
	  if (file.exists(fml_file)) {
	    load(file = fml_file)
	    dots <- c(dots, prnSig_formulas)
	  } else {
	    stop("Run `prnSig()` first.")
	  }
	} else {
	  match_fmls(fmls)
	  dots <- c(dots, fmls)
	}
	
	# load(file = file.path(dat_dir, "Calls\\prnSig_formulas.rda"))
	# dots <- c(rlang::enexprs(...), prnSig_formulas)
		
	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
	info_anal(df = !!df, id = !!id, filepath = !!filepath, filename = !!filename, 
	          scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
	          anal_type = "GSVA")(lm_method, gset_nms, var_cutoff, pval_cutoff, logFC_cutoff, !!!dots)
}


#' Perform GSVA tests
#'
#' logFC_cutoff subsets data for adjusted pvals
#'
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$%
#' @importFrom outliers grubbs.test
#' @importFrom broom.mixed tidy
gsvaTest <- function(df = NULL, id = "entrez", label_scheme_sub = NULL, filepath = NULL, filename = NULL, 
                     complete_cases = FALSE, lm_method = "limma", gsets = NULL, 
                     var_cutoff = .5, pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), ...) {

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  if (purrr::is_empty(fmls)) stop("Formula(s) of contrasts not available.", call. = FALSE)
  
  stopifnot(length(gsets) > 0)
  
  res_es <- df %>% 
    filterData(var_cutoff = var_cutoff) %>% 
    as.matrix() %>% 
    GSVA::gsva(gsets, !!!dots) %>% 
    data.frame(check.names = FALSE) %>% 
    tibble::rownames_to_column("term") %T>%
    write.table(file.path(filepath, paste0(fn_prefix, "_ES.txt")), sep = "\t", col.names = TRUE, row.names = FALSE)    

  quietly_log <- 
    purrr::map(fmls, ~ purrr::quietly(model_onechannel)(
      res_es %>% tibble::column_to_rownames("term"), id = !!id,
      .x, label_scheme_sub, complete_cases, method = lm_method,
      var_cutoff, pval_cutoff, logFC_cutoff
    )) 
  
  out_path <- file.path(dat_dir, "Protein\\GSVA\\log\\prnGSVA_log.txt")
  
  purrr::map(quietly_log, ~ {
    .x[[1]] <- NULL
    return(.x)
  }) %>% 
    reduce(., `c`) %>% 
    purrr::walk(., write, out_path, append = TRUE)
  
  res <- purrr::map(quietly_log, `[[`, 1) %>%
    do.call("cbind", .) %>% 
    tibble::rownames_to_column("term") %>% 
    rm_pval_whitespace()
  
  rm(quietly_log)
  
  kept_rows <- res %>%
    tibble::column_to_rownames("term") %>%
    dplyr::select(grep("pVal", names(.))) %>%
    dplyr::mutate(Kept = rowSums(!is.na (.)) > 0) %>%
    dplyr::select(Kept) %>%
    unlist()
  
  qvals <- res[kept_rows, ] %>%
    `rownames<-`(.$term) %>%
    dplyr::select(grep("pVal", names(.))) %>%
    my_padj(pval_cutoff) %>%
    tibble::rownames_to_column("term")
  
  res <- res %>%
    tibble::column_to_rownames("term") %>%
    dplyr::select(-which(names(.) %in% names(qvals))) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::right_join(qvals, by = "term") %T>% 
    write.table(file.path(filepath, paste0(gsub("\\..*$", "", gsub("\\..*$", "", filename)), "_pVals.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
}


