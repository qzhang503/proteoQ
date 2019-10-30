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
#'\code{logFC_cutoff} are used to filter out low influence genes.
#'
#'The formula(s) of contrast(s) used in \code{\link{prnSig}} will be taken by
#'default.
#'
#'\code{prnGSPAHM} visualizes the distance heat map between essential and all
#'gene sets. The distance is \eqn{1 - f} where \eqn{f} is the fraction overlap
#'for ids under a gene set.
#'
#'@section \code{Protein_GSPA_...csv}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Descrption}\cr term \tab the name
#'  of a gene set \cr is_essential \tab a logicial indicator of gene set
#'  essentiality \cr count \tab the number of entrez IDs with significance
#'  \code{pVal} <= \code{pval_cutoff} under a \code{term} \cr contrast \tab a
#'  contrast of sample groups \cr p_val \tab significance p values \cr q_val
#'  \tab \code{p_val} with \code{BH} adjustment for multiple tests \cr log2fc
#'  \tab the fold change of a gene set at logarithmic base of 2 \cr }
#'
#'@section \code{essmap_Protein_GSPA_...csv}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Descrption}\cr ess_term \tab the
#'  name of an essential gene set \cr term \tab the name of a gene set \cr
#'  ess_count \tab the number of entrez IDs under a \code{term} with matches to
#'  a correspoding \code{ess_term} \cr fraction \tab a fraction of matches in
#'  entrez IDs between a \code{term} and a \code{ess_term}; for example at
#'  \code{fraction = 1} between a \code{term}, \eqn{x}, and an \code{ess_term},
#'  \eqn{y}, all of the gene IDs under \eqn{x} are present in \eqn{y} \cr
#'  distance \tab 1 - \code{fraction} \cr source \tab a numeric index of
#'  essential gene sets \cr target \tab a numeric index of all gene sets \cr }
#'
#'@inheritParams proteoSigtest
#'@inheritParams  proteoEucDist
#'@inheritParams  proteoHM
#'@inheritParams  info_anal
#'@param gset_nms Character vector containing the name(s) of gene sets for
#'  enrichment analysis. The possible values are in \code{c("go_sets",
#'  "kegg_sets", "c2_msig")}.
#'@param filepath Use system default.
#'@param filename Use system default.
#'@param var_cutoff Not currently used.
#'@param pval_cutoff The cut-off in protein significance \code{pVal}. Entries
#'  with \code{pVals} smaller than the threshold will be ignored.
#'@param logFC_cutoff The cut-off in protein \code{log2FC}. Entries with
#'  \code{log2FC} smaller than the threshold will be ignored.
#'@param min_size Minimum number of protein entries for consideration of a gene
#'  set test. The number is the sum of both up and down-expressed proteins after
#'  data filtration by `pval_cutoff`, `logFC_cutoff` or varargs `filter_`.
#'@param gspval_cutoff The cut-off in significance \code{pVal} of gene sets.
#'  Gene sets with \code{pVals} smaller than the threshold will be ignored.
#'@param annot_cols A character vector of column keys that can be found in
#'  \code{essmap_Protein_GSPA_...csv} output files. The values under the
#'  selected keys will be used to color-code sample IDs on the top of heat maps
#'  by \code{prnGSPAHM}.
#'@param annot_colnames A character vector of replacement name(s) to
#'  \code{annot_cols}.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}.
#'@import dplyr rlang ggplot2 networkD3
#'@importFrom magrittr %>%
#' @examples
#' # enrichment analysis
#' prnGSPA(
#'   scale_log2r = TRUE,
#'   impute_na = FALSE,
#'   pval_cutoff = 5E-2,
#'   gspval_cutoff = 5E-3,
#'   gset_nms = c("go_sets", "kegg_sets"),
#'
#'   filter_by_npep = exprs(prot_n_pep >= 2),
#'   # `filename(s)` will be automated,
#' )
#'
#'
#' # distance heat map and network of GSPA terms
#' # a `term` is subset to an `ess_term` if the distance is zero
#' # `source` is a column key in `essmap_Protein_GSPA...csv`
#' prnGSPAHM(
#'   filter_by = exprs(distance <= .2),
#'   annot_cols = "source",
#'   annot_colnames = "Essential set index",
#'
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 8,
#'   border_color = "grey60",
#'   cellwidth = 24,
#'   cellheight = 24,
#'   filename = show_redundancy_at_small_dist.png,
#' )
#'
#' # human terms only
#' prnGSPAHM(
#'   filter_num = exprs(distance <= .8),
#'   filter_sp = exprs(start_with_str("hs", term)),
#'   annot_cols = "source",
#'   annot_colnames = "Essential set index",
#'
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 8,
#'   border_color = "grey60",
#'   cellwidth = 24,
#'   cellheight = 24,
#'   filename = show_connectivity_at_large_dist.png,
#' )
#'
#' # custom color palette
#' prnGSPAHM(
#'   annot_cols = c("source", "ess_count"),
#'   annot_colnames = c("Essential set index", "Count"),
#'   filter_by = exprs(ess_count >= 8, distance <= .80),
#'
#'   fontsize = 16,
#'   fontsize_row = 20,
#'   fontsize_col = 20,
#'   fontsize_number = 8,
#'   border_color = "grey60",
#'   cellwidth = 24,
#'   cellheight = 24,
#'   color = colorRampPalette(c("blue", "white", "red"))(100),
#'   filename = "custom_colors.png"
#' )
#'
#' \dontrun{
#' }
#'
#'@export
proteoGSPA <- function (id = gene, scale_log2r = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                        impute_na = TRUE, complete_cases = FALSE, gset_nms = "go_sets", 
                        var_cutoff = .5, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                        gspval_cutoff = 1E-2, min_size = 10, task = "anal", ...) {

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
	
	if (purrr::is_empty(fmls)) {
	  fml_file <-  file.path(dat_dir, "Calls\\prnSig_formulas.Rdata")
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
	
	mget(names(formals()), rlang::current_env()) %>% save_call("prnGSPA")
	
	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]

	info_anal(df = !!df, id = !!id, scale_log2r = scale_log2r, 
	          filepath = !!filepath, filename = !!filename, impute_na = impute_na, 
	          anal_type = "GSPA")(complete_cases, gset_nms, 
	                              var_cutoff, pval_cutoff, logFC_cutoff, gspval_cutoff, min_size, 
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
gspaTest <- function(df, id = "entrez", label_scheme_sub, filepath, filename, complete_cases = FALSE,
                     gset_nms = "go_sets", var_cutoff = .5, pval_cutoff = 5E-2,
                     logFC_cutoff = log2(1.2), gspval_cutoff = 5E-2, min_size = 6, ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  dots = rlang::enexprs(...)
  
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]

  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]

  if (purrr::is_empty(fmls)) stop("Formula(s) of contrasts not available.", call. = FALSE)

  species <- df$species %>% unique() %>% as.character()
  gsets <- load_dbs(gset_nms, species)
  
  formulas <- names(df) %>% 
    .[grepl("pVal", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique()
  
  col_ind <- purrr::map(formulas, ~ grepl(.x, names(df))) %>%
    dplyr::bind_cols() %>%
    rowSums() %>%
    `>`(0)
  
  cat("Column keys available for data filtration are in `Protein\\Model\\Protein_pVals.txt`.\n")

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)

  if (length(formulas) > 0) purrr::map(formulas, fml_gspa, df = df, col_ind = col_ind, 
                                       id = !!id, gsets = gsets, pval_cutoff, logFC_cutoff, 
                                       gspval_cutoff = gspval_cutoff, min_size = min_size, 
                                       filepath = filepath, filename = filename, 
                                       !!!dots)
}


#' A helper function for GSPA
#'
#' @import purrr dplyr rlang
#' @importFrom magrittr %>%
fml_gspa <- function (df, formula, col_ind, id, gsets, pval_cutoff, logFC_cutoff, gspval_cutoff, min_size, 
                      filepath, filename, ...) {
  
  gapa_summary <- function(gsets, df, min_size = 1) {
    df_sub <- df %>% 
      dplyr::filter(.[["entrez"]] %in% gsets) 
    
    if (nrow(df_sub) > min_size) {
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

  dir.create(file.path(filepath, formula), recursive = TRUE, showWarnings = FALSE)
  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- prep_gspa(df = df, formula = formula, col_ind = col_ind, 
                  pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff) 
  
  df <- df %>% 
    dplyr::mutate(p_val = -log10(p_val)) %>% 
    dplyr::mutate(valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>% 
    dplyr::mutate(valence = factor(valence, levels = c("neg", "pos"))) 
  
  # imputation of NA to 0 increases the number of entries in descriptive "mean" calculation, 
  # like a penalty function
  df <- df %>% 
    tidyr::complete(entrez, contrast, valence) %>% 
    dplyr::mutate(p_val = ifelse(is.na(p_val), 0, p_val)) # %>% 
    # dplyr::mutate(log2Ratio = ifelse(is.na(log2Ratio), 0, log2Ratio))
  
  contrast_groups <- unique(df$contrast) %>% as.character()
  gsets <- gsets %>% .[purrr::map_lgl(., ~ length(.x) > min_size)]
  
  res <- purrr::map(gsets, gapa_summary, df, min_size * length(contrast_groups) * 2) 
  idx <- purrr::map_dbl(res, is.null)

  res <- res[!idx] %>% 
    do.call(rbind, .) %>% 
    tibble::rownames_to_column("term") %>% 
    dplyr::mutate(term = factor(term)) %>% 
    dplyr::mutate_at(.vars = grep("pVal\\s+", names(.)), ~ abs(.x)) %>% 
    dplyr::group_by(term) %>% 
    dplyr::arrange(term) %>% 
    data.frame(check.names = FALSE) %>% 
    dplyr::mutate_at(vars(grep("^pVal\\s+\\(", names(.))), ~ 1/10^.x)

  pass <- purrr::map(res[paste0("pVal (", contrast_groups, ")")], ~ .x <= gspval_cutoff) %>%
    dplyr::bind_cols() %>%
    rowSums() %>%
    `>`(0)
  
  res_pass <- res[pass, ]
  
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

  out <- dplyr::bind_cols(pval, adjp, log2fc) %>% 
    dplyr::mutate_at(.vars = grep("^p_val$|^adjP$", names(.)), format, scientific = TRUE, digits = 2) %>%
    dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2)
  
  # sig_sets: column "1" - significant gsets; column "2" - entrez ids can be found from `df`
  gsets <- gsets %>% 
    .[!grepl("molecular_function$", names(.))] %>% 
    .[!grepl("cellular_component$", names(.))] %>% 
    .[!grepl("biological_process$", names(.))]
  
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
    dplyr::summarise(count = n()) %>% # the number of entries matched to input `df`
    dplyr::right_join(., out, by = "term")

  sig_greedy_terms <- sig_sets %>% 
    RcppGreedySetCover::greedySetCover(FALSE) %>% 
    dplyr::select(term) %>% 
    dplyr::filter(!duplicated(term))
  
  fn_prefix <- paste0(gsub("\\.[^.]*$", "", filename), "_", formula)  
  
  out <- sig_greedy_terms %>% 
    dplyr::mutate(is_essential = TRUE) %>% 
    dplyr::right_join(., out, by = "term") %T>% 
    write.csv(file.path(filepath, formula, paste0(fn_prefix, ".csv")), row.names = FALSE)

  map_essential(sig_sets, unique(sig_greedy_terms$term)) %>% 
    write.csv(file.path(filepath, formula, paste0("essmap_", fn_prefix, ".csv")), row.names = FALSE)
}


#' A helper function for GSPA
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
prep_gspa <- function(df, formula, col_ind, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2)) {
  df <- df %>%
    dplyr::select(grep(formula, names(.), fixed = TRUE)) %>%
    `colnames<-`(gsub(paste0(formula, "."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .) %>% 
    rm_pval_whitespace %>% 
    dplyr::select(grep("^entrez$|^pVal|^log2Ratio", names(.))) %>% 
    dplyr::mutate(entrez = as.character(entrez)) 
  
  contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
    gsub("^log2Ratio\\s+\\(|\\)$", "", .)

  pvals <- df %>% 
    dplyr::select(-grep("^log2Ratio\\s+\\(", names(.))) %>% 
    `names<-`(gsub("^pVal\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = contrast, value = p_val, -entrez) # %>% 
    # tidyr::unite(key, entrez, contrast, sep = "_", remove = FALSE) 
  
  log2rs <- df %>% 
    dplyr::select(-grep("^pVal\\s+\\(", names(.))) %>% 
    `names<-`(gsub("^log2Ratio\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = contrast, value = log2Ratio, -entrez) %>% 
    dplyr::select(-entrez, -contrast)
  
  df <- dplyr::bind_cols(pvals, log2rs) %>% 
    dplyr::filter(p_val <= pval_cutoff) %>% 
    dplyr::filter(abs(log2Ratio) >= logFC_cutoff) %>% 
    dplyr::filter(!is.na(entrez)) %>% 
    dplyr::mutate(contrast = factor(contrast, levels = contrast_groups)) %>%
    dplyr::arrange(contrast) # %>% 

  return(df)
}


#' A helper function for mapping btw gene sets and essential gene sets
#'
#' @import purrr dplyr rlang RcppGreedySetCover
#' @importFrom magrittr %>%
map_essential <- function (sig_sets, sig_greedy_terms) {
  # sig_sets with terms that can be found in sig_greedy_terms
  sig_greedy_sets <- sig_sets %>% dplyr::filter(term %in% sig_greedy_terms)

  # distribution of sig_sets under sig_greedy_sets
  purrr::map(unique(sig_sets$term), ~ {
    curr_sig_set <- sig_sets %>% dplyr::filter(term %in% .x)
    
    sig_greedy_sets %>% 
      dplyr::filter(id %in% curr_sig_set$id) %>% 
      dplyr::group_by(term) %>%
      dplyr::summarise(ess_count = n()) %>% 
      dplyr::mutate(fraction = ess_count/nrow(curr_sig_set)) %>% 
      dplyr::mutate(curr_sig_set = .x)
  }, sig_sets, sig_greedy_sets) %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(curr_sig_set) %>%
    dplyr::rename(ess_term = "term", term = "curr_sig_set") %>%
    dplyr::select(c("ess_term", "term", "ess_count", "fraction")) %>% 
    dplyr::mutate(fraction = round(fraction, digits = 3)) %>% 
    dplyr::ungroup(ess_term) %>% 
    dplyr::mutate(distance = 1 - fraction) %>% 
    dplyr::mutate(term = factor(term, levels = unique(as.character(term)))) %>% 
    dplyr::mutate(ess_term = factor(ess_term, levels = levels(term))) %>%
    dplyr::mutate(source = as.numeric(ess_term), target = as.numeric(term)) 
}


#'GSPA
#'
#'\code{prnGSPAHM} is a wrapper of \code{\link{proteoGSPA}} for heat map and
#'network visualization
#'
#'@rdname proteoGSPA
#'
#'@import purrr
#'@export
prnGSPAHM <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  err_msg <- "Don't call the function with arguments `annot_cols` and/or `annot_colnames`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)
  
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  dir.create(file.path(dat_dir, "Protein\\GSPA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoGSPA)(id = gene, task = plothm, 
                                           annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, 
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
  
  formulas <- names(fmls)
  if (length(formulas) > 0) {
    purrr::walk(formulas, fml_gspahm, filepath, filename, !!!dots)
  }
}


#' A helper function for gspaHM
#'
#' @import purrr dplyr rlang
#' @importFrom magrittr %>%
fml_gspahm <- function (formula, filepath, filename, ...) {
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  ins <- list.files(path = file.path(filepath, formula), pattern = "^essmap_Protein_GSPA_.*\\.csv$")
  stopifnot(length(ins) == 1)

  all_by_greedy <- tryCatch(read.csv(file.path(filepath, formula, ins), check.names = FALSE, header = TRUE, 
                                     comment.char = "#"), error = function(e) NA)
  
  all_by_greedy <- all_by_greedy %>% 
   filters_in_call(!!!filter_dots)
  
  if (nrow(all_by_greedy) == 0) stop("No GSPA terms available after data filtration.")
  
  if (!is.null(dim(all_by_greedy))) {
    message(paste("Essential GSPA loaded."))
  } else {
    stop("Essential GSPA not found.")
  }

  ess_vs_all <- all_by_greedy %>% 
    dplyr::select(-which(names(.) %in% c("ess_count", "source", "target", "distance"))) %>% 
    tidyr::spread(term, fraction) %>% 
    dplyr::mutate_at(vars(which(names(.) != "ess_term")), ~ {1 - .x}) %>% 
    tibble::column_to_rownames("ess_term")
  
  d_row <- dist(ess_vs_all)
  d_col <- dist(t(ess_vs_all))
  max_d_row <- max(d_row, na.rm = TRUE)
  max_d_col <- max(d_col, na.rm = TRUE)
  d_row[is.na(d_row)] <- 1.2 * max_d_row
  d_col[is.na(d_col)] <- 1.2 * max_d_col

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
    annotation_row <- ess_vs_all %>% 
      tibble::rownames_to_column("ess_term") %>% 
      dplyr::filter(!duplicated(ess_term)) %>% 
      dplyr::select(ess_term)
  }  
  
  if (is.null(dots$annotation_colors)) {
    annotation_colors <- setHMColor(annotation_col)
  } else if (is.na(dots$annotation_colors)) {
    annotation_colors <- NA
  } else {
    annotation_colors <- eval(dots$annotation_colors, env = caller_env())
  }

  dots <- dots %>% 
    .[!names(.) %in% c("annot_cols", "annot_colnames", "annot_rows", 
                       "mat", "filename", "annotation_col", "annotation_row", 
                       "color", "annotation_colors", "breaks", 
                       "cluster_rows", "cluster_cols")]
                       
  ph <- my_pheatmap(
    mat = ess_vs_all,
    filename = file.path(filepath, formula, filename),
    annotation_col = annotation_col,
    annotation_row = NA, # not used
    color = mypalette,
    annotation_colors = annotation_colors, 
    breaks = NA, # not used
    cluster_rows =  hclust(d_row), 
    cluster_cols = hclust(d_col),
    !!!dots,
  )

  # networks
  cluster <- data.frame(cluster = cutree(ph$tree_col, h = max_d_col)) %>% 
    tibble::rownames_to_column("term")

  all_by_greedy <- all_by_greedy %>% 
    dplyr::mutate(term = factor(term, levels = unique(as.character(term)))) %>% 
    dplyr::mutate(ess_term = factor(ess_term, levels = levels(term))) %>%
    dplyr::mutate(source = as.numeric(ess_term), target = as.numeric(term)) 

  min_target <- min(all_by_greedy$target, na.rm = TRUE)
  all_by_greedy <- all_by_greedy %>% 
    dplyr::mutate(source = source - min_target, target = target - min_target)
  rm(min_target)
  
  all_by_greedy <- all_by_greedy %>% 
    dplyr::left_join(cluster, by = "term")
  
  my_nodes <- all_by_greedy %>% 
    dplyr::arrange(-ess_count) %>% 
    dplyr::select(c("term", "cluster", "ess_count")) %>% 
    dplyr::filter(!duplicated(term)) %>% 
    dplyr::arrange(cluster) %>% 
    dplyr::mutate(term = factor(term, levels = as.character(term))) %>% 
    dplyr::arrange(term)
  
  my_links <- all_by_greedy %>% 
    dplyr::mutate(term = factor(term, levels = levels(my_nodes$term))) %>% 
    dplyr::mutate(ess_term = factor(ess_term, levels = levels(my_nodes$term))) %>% 
    dplyr::mutate(source = as.numeric(term), target = as.numeric(ess_term)) %>% 
    dplyr::arrange(term) %>% 
    dplyr::mutate(source = source - 1, target = target -1) %>% 
    dplyr::select(source, target, fraction) %>% 
    dplyr::mutate(fraction = fraction * 10) %>% 
    dplyr::distinct() 

  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  networkD3::forceNetwork(Links = my_links, Nodes = my_nodes, Source = "source",
                          Target = "target", Value = "fraction", NodeID = "term", Nodesize = "ess_count", 
                          Group = "cluster", opacity = 0.8, zoom = TRUE) %>% 
    networkD3::saveNetwork(file = file.path(filepath, formula, paste0(fn_prefix, ".html")))
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




