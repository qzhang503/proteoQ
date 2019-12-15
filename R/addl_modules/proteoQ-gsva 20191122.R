#'GSVA of protein data
#'
#'\code{prnGSVA} performs the GSVA aganist protein \code{log2FC}
#'
#'The formula(s) of contrast(s) used in \code{\link{prnSig}} will be taken by
#'default.
#'
#'
#'@inheritParams  proteoEucDist
#'@inheritParams  proteoHM
#'@inheritParams  info_anal
#'@inheritParams proteoSigtest
#'@param lm_method Character string indicating the linear modeling method for
#'  significance assessment of GSVA enrichment scores
#'@param gset_nm Character vector of the name(s) of gene sets. The currently
#'  available data sets include c("go_sets", "kegg_sets", "c2_msig").
#'@param filepath Use system default.
#'@param filename Use system default.
#'@param Parameters for \code{\link[GSVA]{gsva}}
#'
#'@import dplyr rlang ggplot2 GSVA
#'@importFrom magrittr %>%
#' @examples
#'prnGSVA(
#'  scale_log2r = TRUE, 
#'  impute_na = FALSE, 
#'  min.sz = 10, 
#'  verbose = FALSE, 
#'  parallel.sz = 0, 
#'  mx.diff = TRUE, 
#'  gset_nms = c("go_sets", "kegg_sets"), 
#')
#'
#'@export
prnGSVA <- function (id = gene,
										scale_log2r = TRUE, df = NULL, filepath = NULL, filename = NULL,
										impute_na = TRUE, complete_cases = FALSE, lm_method = "limma",
										gset_nms = "go_sets", var_cutoff = .5, pval_cutoff = 1E-4,
										logFC_cutoff = log2(1.1), ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  id <- rlang::enexpr(id)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	lm_method <- rlang::as_string(rlang::enexpr(lm_method))

	stopifnot(rlang::is_logical(scale_log2r))
	stopifnot(rlang::is_logical(impute_na))
	stopifnot(rlang::is_logical(complete_cases))
	
	dots = rlang::enexprs(...)
	fmls <- dots %>% .[grepl("~", .)]
	dots <- dots[!names(dots) %in% names(fmls)]
	
	if(purrr::is_empty(fmls)) {
	  fml_file <-  file.path(dat_dir, "Calls\\prnSig_formulas.Rdata")
	  if(file.exists(fml_file)) {
	    load(file = fml_file)
	    dots <- c(dots, prnSig_formulas)
	  } else {
	    stop("Run `prnSig()` first.")
	  }
	} else {
	  match_fmls(fmls)
	  dots <- c(dots, fmls)
	}
	
	# load(file = file.path(dat_dir, "Calls\\prnSig_formulas.Rdata"))
	# dots <- c(rlang::enexprs(...), prnSig_formulas)
		
	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
	info_anal(df = !!df, id = !!id, scale_log2r = scale_log2r,
					filepath = !!filepath, filename = !!filename,
					impute_na = impute_na,
					anal_type = "GSVA")(complete_cases, lm_method, gset_nms, var_cutoff, pval_cutoff,
					                    logFC_cutoff, !!!dots)
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
                     complete_cases = FALSE, lm_method = "limma", gsets = NULL, var_cutoff = .5, pval_cutoff = 1E-4, 
                     logFC_cutoff = log2(1.1), ...) {
  
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
  
  run_scripts <- FALSE
  if (run_scripts) {
    is_null <- purrr::map_lgl(gsets, ~ is.null(.x))
    purrr::walk2(is_null, names(is_null), ~ {
      if (.x) warning("Gene set: `", .y, "` not found", call. = FALSE)
    })
    gsets[is_null] <- NULL    
  }
  
  stopifnot(length(gsets) > 0)
  
  res_es <- df %>% 
    filterData(var_cutoff = var_cutoff) %>% 
    as.matrix() %>% 
    GSVA::gsva(gsets, !!!dots) %>% 
    data.frame(check.names = FALSE) %>% 
    tibble::rownames_to_column("term") %T>%
    write.table(file.path(filepath, paste0(fn_prefix, "_ES.txt")), sep = "\t", col.names = TRUE, row.names = FALSE)    

  res <- purrr::map(fmls,
             ~ model_onechannel(res_es %>% tibble::column_to_rownames("term"), id = !!id,
                                .x, label_scheme_sub, complete_cases, method = lm_method,
                                var_cutoff, pval_cutoff, logFC_cutoff)) %>%
    do.call("cbind", .) %>%
    tibble::rownames_to_column("term") %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.character) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ gsub("\\s*", "", .x) ) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.numeric)
  
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








