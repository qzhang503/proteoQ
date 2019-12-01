#'GSPA of protein data
#'
#'\code{prnGSPA} performs the analysis of Gene Set Probability Asymmetricity
#'(GSPA) aganist protein \code{log2FC} data. Users should avoid call the method
#'directly, but instead use the following wrappers.
#'
#'The significance \code{pVals} of proteins are first obtained from
#'\code{\link{prnSig}}, which involves linear modelings using
#'\code{\link[limma]{eBayes}} or \code{\link[lmerTest]{lmer}}. The geometric
#'mean of \code{pVals} are then each calculated for the groups of up or down
#'regulated proteins. The quotient of the two \code{pVals} is used to prepsent
#'the significance of gene set enrichment. The arguments \code{pval_cutoff} and
#'\code{logFC_cutoff} are used to filter out low influence genes. Additional
#'subsetting of data via the \code{vararg} approach of \code{filter_} is
#'feasible.
#'
#'The formula(e) of contrast(s) used in \code{\link{prnSig}} will be taken by
#'default.
#'
#'@inheritParams proteoSigtest
#'@inheritParams  proteoEucDist
#'@inheritParams  proteoHM
#'@inheritParams  info_anal
#'@param id Character string indicating the type of data. The value will be
#'  determined automatically.
#'@param filename A file name to output results. By default, it will be
#'  determined automatically by the name of the calling function and the value
#'  of id in the call.
#'@param gset_nms Character string or vector containing the name(s) of gene sets
#'  for enrichment analysis. The possible values are in \code{c("go_sets",
#'  "kegg_sets")}. Note that the currently supported species are human, mouse
#'  and rat.
#'@param var_cutoff Not currently used.
#'@param pval_cutoff Numeric value or vector; the cut-off in protein
#'  significance \code{pVal}. Entries with \code{pVals} less significant than
#'  the threshold will be ignored during enrichment analysis.
#'@param logFC_cutoff Numeric value or vector; the cut-off in protein
#'  \code{log2FC}. Entries with absolute \code{log2FC} smaller than the
#'  threshold will be ignored during enrichment analysis.
#'@param min_size Numeric value or vector; minimum number of protein entries for
#'  consideration in gene set tests. The number is the sum of up or
#'  down-expressed proteins after data filtration by \code{pval_cutoff},
#'  \code{logFC_cutoff} or varargs expressions under \code{filter_}.
#'@param max_size Numeric value or vector; maximum number of protein entries for
#'  consideration in gene set tests. The number is the sum of up or
#'  down-expressed proteins after data filtration by \code{pval_cutoff},
#'  \code{logFC_cutoff} or varargs expressions under \code{filter_}. The default
#'  in infinite. For samples that are vastly different between contrasts,
#'  essential gene sets may be blackholed with trivial terms such as cell parts,
#'  molecular functions et al. In cases like this, a smaller \code{max_size} may
#'  be considered.
#'@param min_greedy_size Numeric value or vector; minimum number of unique
#'  protein entries for a gene set to be considered essential.
#'@param gspval_cutoff Numeric value or vector; the cut-off in gene-set
#'  significance \code{pVal}. Only enrichment terms with \code{pVals} more
#'  significant than the threshold will be reported.
#'@param fml_nms Character string or vector; the forumula name(s). By default,
#'  the names match to those used in \code{\link{pepSig}} or
#'  \code{\link{prnSig}}.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}.
#'@import dplyr rlang ggplot2 networkD3
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/pepseqmod_min.R
#'@example inst/extdata/examples/normPep_min.R
#'@example inst/extdata/examples/normPrn_min.R
#'@example inst/extdata/examples/imputeNA_examples.R
#'@example inst/extdata/examples/gspa_examples.R
#'
#'@export
proteoGSPA <- function (id = gene, scale_log2r = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                        impute_na = TRUE, complete_cases = FALSE, gset_nms = "go_sets", 
                        var_cutoff = .5, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                        gspval_cutoff = 5E-2, min_size = 10, max_size = Inf, min_greedy_size = 1, 
                        fml_nms = NULL, task = "anal", ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  id <- rlang::enexpr(id)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	task <- rlang::enexpr(task)

	stopifnot(is_logical(scale_log2r), is_logical(impute_na), is_logical(complete_cases))

	dots <- rlang::enexprs(...)
	fmls <- dots %>% .[grepl("^\\s*~", .)]
	dots <- dots[!names(dots) %in% names(fmls)]
	dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSPA")
	
	curr_call <- mget(names(formals()), rlang::current_env()) %>% 
	  .[names(.) != "..."] %>% 
	  c(dots) 
	
	if (task == "anal") {
	  curr_call %>% save_call("prnGSPA")
	} else if (task == "plothm") {
	  curr_call %>% save_call("prnGSPAHM")
	}

	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]

	info_anal(df = !!df, id = !!id, scale_log2r = scale_log2r, 
	          filepath = !!filepath, filename = !!filename, impute_na = impute_na, 
	          anal_type = "GSPA")(complete_cases = complete_cases, gset_nms = gset_nms, 
	                              var_cutoff = var_cutoff, pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, 
	                              gspval_cutoff = gspval_cutoff, 
	                              min_size = min_size, max_size = max_size, min_greedy_size = min_greedy_size, 
	                              task = !!task, !!!dots)
}

#'GSPA
#'
#'\code{prnGSPA} is a wrapper of \code{\link{proteoGSPA}} for enrichment analysis
#'
#'@rdname proteoGSPA
#'
#'@import purrr
#'@export
prnGSPA <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\GSPA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoGSPA)(id = gene, task = anal, ...)
  quietly_log$result <- NULL
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\GSPA\\log","prnGSPA_log.csv"), append = TRUE)
}


#' Perform GSPA tests
#'
#' logFC_cutoff setset of data before the calculation of adjusted pvals
#'
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$% %T>%
#' @importFrom outliers grubbs.test
#' @importFrom broom.mixed tidy
gspaTest <- function(df = NULL, id = "entrez", label_scheme_sub = NULL, 
                     filepath = NULL, filename = NULL, complete_cases = FALSE,
                     gset_nms = "go_sets", var_cutoff = .5, pval_cutoff = 5E-2,
                     logFC_cutoff = log2(1.2), gspval_cutoff = 5E-2, 
                     min_size = 6, max_size = Inf, min_greedy_size = 1, 
                     scale_log2r = TRUE, anal_type = "GSPA", ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  dots = rlang::enexprs(...)
  
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]

  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]

  if (purrr::is_empty(fmls)) stop("Formula(s) of contrasts not available.", call. = FALSE)

  species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
  gsets <- load_dbs(gset_nms, species)
  
  stopifnot(length(gsets) > 0)

  fml_nms <- names(df) %>% 
    .[grepl("pVal\\s*\\(", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique() %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls %>% .[names(.) %in% fml_nms]
  fml_nms <- fml_nms %>% .[map_dbl(., ~ which(.x == names(fmls)))]

  if (is_empty(fml_nms)) stop("No formula matached; compare the formula name(s) with those in `prnSig(..)`")

  col_ind <- purrr::map(fml_nms, ~ grepl(.x, names(df))) %>%
    dplyr::bind_cols() %>%
    rowSums() %>%
    `>`(0)
  
  stopifnot(sum(col_ind) > 0)
  
  cat("Column keys available for data filtration are in `Protein\\Model\\Protein_pVals.txt`.\n")

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  # var_cutoff only for GSEA
  # min_greedy_size and gspval_cutoff only for GSPA
  if (anal_type == "GSPA") {
    purrr::pmap(list(fmls, fml_nms, pval_cutoff, logFC_cutoff, gspval_cutoff, min_size, max_size, min_greedy_size), 
                fml_gspa, df = df, col_ind = col_ind, id = !!id, gsets = gsets, label_scheme_sub = label_scheme_sub, 
                complete_cases = complete_cases, scale_log2r = scale_log2r, 
                filepath = filepath, filename = filename, !!!dots)
  } else if (anal_type == "GSEA") {
    # need scale_log2r for var_cutoff
    # need formula for pval_cutoff and logFC_cutoff
    purrr::pmap(list(fmls, fml_nms, var_cutoff, pval_cutoff, logFC_cutoff, gspval_cutoff, min_size, max_size), 
                fml_gsea, df, col_ind, id = !!id, gsets, label_scheme_sub, complete_cases, scale_log2r, 
                filepath, filename, !!!dots)
    
  } else {
    stop("Unhandled task.")
  }
}


#' A helper function for GSPA
#'
#' @import purrr dplyr rlang
#' @importFrom magrittr %>% %T>%
fml_gspa <- function (fml, fml_nm, pval_cutoff, logFC_cutoff, gspval_cutoff, min_size, max_size, min_greedy_size, 
                      df, col_ind, id, gsets, label_scheme_sub, complete_cases, scale_log2r, 
                      filepath, filename, ...) {

  gapa_summary <- function(gsets, df, min_size = 1) {
    df_sub <- df %>% 
      dplyr::filter(.[["entrez"]] %in% gsets) 
    
    if (nrow(df_sub) >= min_size) {
      delta_p <- df_sub %>% 
        dplyr::group_by(contrast, valence) %>% 
        dplyr::summarise(p_val = mean(p_val, na.rm = TRUE)) %>% 
        tidyr::spread(key = contrast, value = p_val) %>% 
        dplyr::select(-valence) %>% 
        `colnames<-`(paste0("pVal (", colnames(.), ")"))
      
      delta_p <- delta_p[2, ] - delta_p[1, ]
      
      delta_fc <- df_sub %>% 
        dplyr::group_by(contrast, valence) %>% 
        dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) %>% 
        tidyr::spread(key = contrast, value = log2Ratio) %>% 
        dplyr::select(-valence) %>% 
        `colnames<-`(paste0("log2Ratio (", colnames(.), ")")) %>% 
        abs(.) # may leave as an option later for dual significance
      
      delta_fc <- delta_fc[2, ] - delta_fc[1, ]
      
      return(dplyr::bind_cols(delta_p, delta_fc))
    } else {
      return(NULL)
    }
  }  

  dir.create(file.path(filepath, fml_nm), recursive = TRUE, showWarnings = FALSE)
  id <- rlang::as_string(rlang::enexpr(id))
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  # penaltize with NA imputation to 0 to increases the number of entries in descriptive "mean" calculation
  df <- df %>% prep_gspa(id = "entrez", fml_nm = fml_nm, col_ind = col_ind, 
                         pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff) %>% 
    dplyr::mutate(p_val = -log10(p_val)) %>% 
    dplyr::mutate(valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>% 
    dplyr::mutate(valence = factor(valence, levels = c("neg", "pos"))) %>% 
    tidyr::complete(entrez, contrast, valence) %>% 
    dplyr::mutate(p_val = ifelse(is.na(p_val), 0, p_val)) # %>% 
    # dplyr::mutate(log2Ratio = ifelse(is.na(log2Ratio), 0, log2Ratio))

  gsets <- gsets %>% 
    .[!grepl("molecular_function$", names(.))] %>% 
    .[!grepl("cellular_component$", names(.))] %>% 
    .[!grepl("biological_process$", names(.))] %>%     
    .[purrr::map_lgl(., ~ length(.x) >= min_size)] %>% 
    .[purrr::map_lgl(., ~ length(.x) <= max_size)]
  
  res_pass <- local({
    contrast_groups <- unique(df$contrast) %>% as.character()
    
    res <- purrr::map(gsets, gapa_summary, df, min_size * length(contrast_groups) * 2) 
    idx <- purrr::map_dbl(res, is.null)
    
    res <- res[!idx] %>% 
      do.call(rbind, .) %>% 
      tibble::rownames_to_column("term") %>% 
      dplyr::mutate(term = factor(term)) %>% 
      dplyr::mutate_at(.vars = grep("^pVal\\s+\\(", names(.)), ~ abs(.x)) %>% 
      dplyr::group_by(term) %>% 
      dplyr::arrange(term) %>% 
      data.frame(check.names = FALSE) %>% 
      dplyr::mutate_at(vars(grep("^pVal\\s+\\(", names(.))), ~ 1/10^.x)
    
    pass <- purrr::map(res[paste0("pVal (", contrast_groups, ")")], ~ .x <= gspval_cutoff) %>%
      dplyr::bind_cols() %>%
      rowSums() %>%
      `>`(0)
    
    res_pass <- res[pass, ]
  })
  
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
      dplyr::mutate_at(.vars = grep("^p_val$|^adjP$", names(.)), format, scientific = TRUE, digits = 2) %>%
      dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2)
  })

  sig_sets <- purrr::map2(gsets, names(gsets), ~ {
    .x %>% 
      tibble() %>% 
      dplyr::mutate(set = .y)
  }) %>% 
    bind_rows() %>% 
    `names<-` (c("id", "term")) %>% 
    dplyr::filter(id %in% unique(df$entrez)) %>% 
    dplyr::select(c("term", "id")) %>% 
    dplyr::filter(term %in% res_pass$term)
  
  out <- sig_sets %>% 
    dplyr::group_by(term) %>% 
    dplyr::summarise(size = n()) %>% # the number of entries matched to input `df`
    dplyr::right_join(out, by = "term")

  res_greedy <- sig_sets %>% 
    RcppGreedySetCover::greedySetCover(FALSE) %T>% 
    write.csv(file.path(filepath, fml_nm, paste0("resgreedy_", fn_prefix, ".csv")), row.names = FALSE) %>% 
    dplyr::group_by(term) %>% 
    dplyr::summarise(ess_size = n()) %>% 
    dplyr::filter(ess_size >= min_greedy_size)
  
  sig_sets <- res_greedy %>% 
    dplyr::right_join(out, by = "term") %>% 
    dplyr::mutate(is_essential = ifelse(!is.na(ess_size), TRUE, FALSE)) %>% 
    dplyr::select(term, is_essential, size, ess_size, contrast, p_val, q_val, log2fc) %T>% 
    write.csv(file.path(filepath, fml_nm, paste0(fn_prefix, ".csv")), row.names = FALSE) %>% 
    dplyr::filter(is_essential) %>% 
    dplyr::filter(!duplicated(term)) %>% 
    dplyr::select(-contrast, -p_val, -q_val, -log2fc) %T>% 
    write.csv(file.path(filepath, fml_nm, paste0("essmeta_", fn_prefix, ".csv")), row.names = FALSE) %>% 
    dplyr::select(term, ess_size) %>% 
    dplyr::right_join(sig_sets, by = "term")

  map_essential(sig_sets) %>% 
    write.csv(file.path(filepath, fml_nm, paste0("essmap_", fn_prefix, ".csv")), row.names = FALSE)
}


#' A helper function for GSPA
#'
#' @import purrr dplyr rlang
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
prep_gspa <- function(df, id, fml_nm, col_ind, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2)) {
  id <- as_string(enexpr(id))
  
  df <- df %>%
    dplyr::select(grep(fml_nm, names(.), fixed = TRUE)) %>%
    `colnames<-`(gsub(paste0(fml_nm, "."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .) %>% 
    rm_pval_whitespace() %>% 
    dplyr::select(id, grep("^pVal|^log2Ratio", names(.))) %>% 
    dplyr::mutate(!!id := as.character(.[[id]]))
  
  contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
    gsub("^log2Ratio\\s+\\(|\\)$", "", .)
  
  pvals <- df %>% 
    dplyr::select(-grep("^log2Ratio\\s+\\(", names(.))) %>% 
    `names<-`(gsub("^pVal\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = contrast, value = p_val, -id) # %>% 
    # tidyr::unite(key, id, contrast, sep = "_", remove = FALSE) 
  
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


#' A helper function for mapping btw gene sets and essential gene sets
#'
#' @import purrr dplyr rlang RcppGreedySetCover
#' @importFrom magrittr %>%
map_essential <- function (sig_sets) {
  ess_terms <- sig_sets %>% 
    dplyr::filter(!is.na(ess_size), !duplicated(term)) %>% 
    dplyr::select(term) %>% 
    unlist

  sig_greedy_sets <- sig_sets %>% dplyr::filter(term %in% ess_terms)

  purrr::map(unique(sig_sets$term), ~ {
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
    dplyr::ungroup(ess_term) %>% 
    dplyr::mutate(distance = 1 - fraction) %>% 
    dplyr::mutate(term = factor(term, levels = unique(as.character(term)))) %>% 
    dplyr::mutate(ess_term = factor(ess_term, levels = levels(term))) %>%
    dplyr::mutate(idx = as.numeric(term), ess_idx = as.numeric(ess_term))
}


#'Heat map visualization of GSPA results
#'
#'\code{prnGSPAHM} visualizes distance heat maps between essential and all gene
#'sets.
#'
#'The list of gene sets and the associative quality metrics of \code{size} and
#'\code{ess_size} are assessed after data filtration with the criteria specified
#'by arguments \code{pval_cutoff} and \code{logFC_cutoff}, as well as optional
#'varargs of \code{filter_}.
#'
#'@section \code{Protein_GSPA_...csv}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Descrption}\cr term \tab a gene
#'  set term \cr is_essential \tab a logicial indicator of gene set essentiality
#'  \cr size \tab the number of IDs under a \code{term} \cr ess_size \tab the
#'  number of IDs that can be found under a corresponding essential set \cr
#'  contrast \tab a contrast of sample groups \cr p_val \tab significance p
#'  values \cr q_val \tab \code{p_val} with \code{BH} adjustment of multiple
#'  tests \cr log2fc \tab the fold change of a gene set at logarithmic base of 2
#'  \cr }
#'
#'@section \code{essmap_Protein_GSPA_...csv}:
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
#'@inheritParams  proteoEucDist
#'@inheritParams proteoHM
#'@param annot_cols A character vector of column keys that can be found in
#'  \code{essmap_.*.csv}. The values under the selected keys will be used to
#'  color-code enrichment terms on the top of heat maps.
#'@param annot_rows A character vector of column keys that can be found from
#'  \code{essmeta_.*.csv} . The values under the selected keys will be used to
#'  color-code essential terms on the side of heat maps.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr \cr Additional arguments for
#'  \code{\link[pheatmap]{pheatmap}}, i.e., \code{fontsize }... \cr Note
#'  arguments disabled from \code{pheatmap}: \cr \code{annotation_col}; instead
#'  use keys indicated in \code{annot_cols} \cr \code{annotation_row}; instead
#'  use keys indicated in \code{annot_rows}
#'@import purrr
#'
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/pepseqmod_min.R
#'@example inst/extdata/examples/normPep_min.R
#'@example inst/extdata/examples/normPrn_min.R
#'@example inst/extdata/examples/imputeNA_examples.R
#'@example inst/extdata/examples/sigtest_min.R
#'@example inst/extdata/examples/prnGSPAHM_examples.R
#'
#'@seealso \code{\link{prnGSPA}} for gene set enrichment analysis.
#'@export
prnGSPAHM <- function (annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, ...) {
  err_msg <- "Duplicated arguments in `annot_cols`, `annot_colnames` or `annot_rows`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)
  
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  annot_rows <- rlang::enexpr(annot_rows)
  
  dir.create(file.path(dat_dir, "Protein\\GSPA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoGSPA)(id = gene, task = plothm, 
                                            annot_cols = !!annot_cols, 
                                            annot_colnames = !!annot_colnames, 
                                            annot_rows = !!annot_rows, 
                                            ...)
  purrr::walk(quietly_log, write, file.path(dat_dir, "Protein\\GSPA\\log","plot_hmGSPA_log.csv"), 
              append = TRUE)  
}


#'Plot distance heat map of GSPA
gspaHM <- function(filepath, filename, ...) {
  dots <- rlang::enexprs(...)
  
  fmls <- dots %>% .[grepl("~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  if (purrr::is_empty(fmls))
    stop("No formula(s) of contrasts available.", call. = TRUE)
  
  fml_nms <- names(fmls)
  if (length(fml_nms) > 0) {
    purrr::walk(fml_nms, fml_gspahm, filepath, filename, !!!dots)
  }
}


#' A helper function for gspaHM
#'
#' @import purrr dplyr rlang pheatmap networkD3
#' @importFrom magrittr %>%
fml_gspahm <- function (fml_nm, filepath, filename, ...) {
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  ins <- list.files(path = file.path(filepath, fml_nm), pattern = "^essmap_.*\\.csv$")
  stopifnot(length(ins) == 1)

  all_by_greedy <- tryCatch(read.csv(file.path(filepath, fml_nm, ins), check.names = FALSE, header = TRUE, 
                                     comment.char = "#"), error = function(e) NA) %>% 
    filters_in_call(!!!filter_dots)
  
  if (nrow(all_by_greedy) == 0) stop("No GSPA terms available after data filtration.")
  
  if (max(all_by_greedy$distance) == 0) stop("Identical, all-zero distance detected; 
                                             try lower the criteria in data filtrations
                                             or rerun `prnGSPA` at more relaxed `gspval_cutoff` threshold.")

  if (!is.null(dim(all_by_greedy))) {
    message(paste("Essential GSPA loaded."))
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
    annotation_col <- gspa_colAnnot(annot_cols = annot_cols, df = all_by_greedy, sample_ids = all_by_greedy$term)
  }

  if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
    colnames(annotation_col) <- annot_colnames
  }

  if (is.null(annot_rows)) {
    annotation_row <- NA
  } else {
    meta_ins <- list.files(path = file.path(filepath, fml_nm), pattern = "^essmeta_.*\\.csv$")
    stopifnot(length(meta_ins) == 1)
    
    ess_meta <- tryCatch(read.csv(file.path(filepath, fml_nm, meta_ins), check.names = FALSE, header = TRUE, 
                                  comment.char = "#"), error = function(e) NA)
    
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
    warning("The plot width is set to", max_width)
    width <- max_width
    height <- pmin(max_height, width * nrow / ncol)
  } 
  
  if ((!is.na(height)) & (height > max_height)) {
    warning("The plot height is set to", max_height)
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
  })

  my_links <- all_by_greedy %>% 
    dplyr::arrange(term) %>% 
    dplyr::select(source, target, fraction) %>% 
    dplyr::mutate(fraction = fraction * 10) %>% 
    dplyr::distinct() 

  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  networkD3::forceNetwork(Links = my_links, Nodes = my_nodes, Source = "source",
                          Target = "target", Value = "fraction", NodeID = "term", Nodesize = "size", 
                          Group = "cluster", opacity = 0.8, zoom = TRUE) %>% 
    networkD3::saveNetwork(file = file.path(filepath, fml_nm, paste0(fn_prefix, ".html")))
}


#' Sets up the column annotation in heat maps
#'
#' @import dplyr rlang
#' @importFrom magrittr %>%
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




