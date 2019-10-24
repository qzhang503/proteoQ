#'GSPA of protein data
#'
#'\code{prnGSPA} performs the analysis of Gene Set Probability Asymmetricity
#'(GSPA) aganist protein \code{log2FC} data
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
#'@param pval_cutoff The cut-off in significance \code{pVal}. Protein entries
#'  with \code{pVals} smaller than the threshold will be ignored.
#'@param logFC_cutoff The cut-off in \code{log2FC}. Protein entries with
#'  \code{log2FC} smaller than the threshold will be ignored.
#'@param min_size Minimum number of protein entries for consideration of a gene
#'  set test.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#' @examples
#' prnGSPA(
#'   scale_log2r = TRUE,
#'   impute_na = FALSE,
#'   pval_cutoff = 5E-2,
#'   gset_nms = c("go_sets", "kegg_sets"),
#'
#'   filter_by_npep = exprs(prot_n_pep >= 2),
#'   # `filename(s)` will be automated,
#' )
#'
#' \dontrun{
#' }
#'
#'@export
proteoGSPA <- function (id = gene, scale_log2r = TRUE, df = NULL, filepath = NULL, 
										filename = NULL, impute_na = TRUE, complete_cases = FALSE, 
										gset_nms = "go_sets", 
										var_cutoff = .5, pval_cutoff = 1E-2, logFC_cutoff = log2(1.2), 
										min_size = 10, 
										task = "anal", ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  id <- rlang::enexpr(id)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	task <- rlang::enexpr(task)

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
	
	mget(names(formals()), rlang::current_env()) %>% save_call("prnGSPA")
	
	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]

	info_anal(df = !!df, id = !!id, scale_log2r = scale_log2r,
					filepath = !!filepath, filename = !!filename,
					impute_na = impute_na,
					anal_type = "GSPA")(complete_cases, gset_nms, 
					                    var_cutoff, pval_cutoff, logFC_cutoff, min_size, 
					                    task = !!task, !!!dots)
}

#'GSPA
#'
#'\code{prnGSPA} is a wrapper of \code{\link{proteoGSPA}} 
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
#' logFC_cutoff subsets data for adjusted pvals
#'
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$% %T>%
#' @importFrom outliers grubbs.test
#' @importFrom broom.mixed tidy
gspaTest <- function(df, id = "entrez", label_scheme_sub, filepath, filename, complete_cases = FALSE,
                     gset_nms = "go_sets", var_cutoff = .5, pval_cutoff = 5E-2,
                     logFC_cutoff = log2(1.2), min_size = 6, ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  dots = rlang::enexprs(...)
  
  fmls <- dots %>% .[grepl("~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]

  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]

  if(purrr::is_empty(fmls))
    stop("No formula(s) of contrasts available.", call. = TRUE)
  
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
  
  # cat("Available column keys for data filtration: \n")
  # cat(paste0(names(df), "\n"))
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)

  if (length(formulas) > 0) purrr::map(formulas, fml_gspa, df = df, col_ind = col_ind, 
                                       id = !!id, gsets = gsets, pval_cutoff, logFC_cutoff, 
                                       filepath = filepath, filename = filename, min_size = min_size, 
                                       !!!dots)
}


#' A helper function for GSPA
#'
#' @import purrr dplyr rlang
#' @importFrom magrittr %>%
fml_gspa <- function (df, formula, col_ind, id, gsets, pval_cutoff, logFC_cutoff, filepath, filename, min_size, ...) {
  
  gapa_summary <- function(gsets, df, min_size = 1) {
    df_sub <- df %>% 
      dplyr::filter(.[["entrez"]] %in% gsets) 
    
    if(nrow(df_sub) > min_size) {
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
  # which is like a penalty function
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

  pass <- purrr::map(res[paste0("pVal (", contrast_groups, ")")], ~ .x <= pval_cutoff) %>%
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
  
  fn_prefix <- paste0(gsub("\\.[^.]*$", "", filename), "_", formula)
  out_nm <- paste0(fn_prefix, ".csv")

  out <- dplyr::bind_cols(pval, adjp, log2fc) %>% 
    dplyr::mutate_at(.vars = grep("^p_val$|^adjP$", names(.)), format, scientific = TRUE, digits = 2) %>%
    dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2)
  
  # find essential sets
  gsets <- gsets %>% 
    .[!grepl("molecular_function$", names(.))] %>% 
    .[!grepl("cellular_component$", names(.))] %>% 
    .[!grepl("biological_process$", names(.))]
  
  all_sets <- purrr::map2(gsets, names(gsets), ~ {
    .x %>% 
      tibble() %>% 
      dplyr::mutate(set = .y)
  }) %>% 
    bind_rows() %>% 
    `names<-` (c("id", "term")) %>% 
    dplyr::filter(id %in% unique(df$entrez)) %>% 
    dplyr::select(c("term", "id")) %>% 
    dplyr::filter(term %in% res_pass$term)
  
  # one id only assign to one term without overlap
  # 969 significant ids in 124 terms
  # non-redundant mapping of ids to greedy terms
  # no redundant ids
  greedy_sets <- all_sets %>% 
    RcppGreedySetCover::greedySetCover(FALSE) %>% 
    dplyr::select(term) %>% 
    dplyr::filter(!duplicated(term)) %>% 
    dplyr::mutate(is_essential = TRUE)
  
  out <- greedy_sets %>% dplyr::right_join(., out, by = "term") %T>% 
    write.csv(file.path(filepath, formula, out_nm), row.names = FALSE)
  
  map_essential(all_sets, greedy_sets, out) %>% 
    write.csv(file.path(filepath, formula, paste0("essmap_", fn_prefix, ".csv")), row.names = FALSE)
  
  invisible(out)

  # dplyr::bind_cols(pval, adjp, log2fc) %>% 
  #   dplyr::mutate_at(.vars = grep("^p_val$|^adjP$", names(.)), format, scientific = TRUE, digits = 2) %>%
  #   dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2) %>%
  #   write.csv(file.path(filepath, formula, out_nm), row.names = FALSE)
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




#' A helper function
#'
#' @import purrr dplyr rlang RcppGreedySetCover
#' @importFrom magrittr %>%
map_essential <- function (all_sets, greedy_sets, out) {
  # 969 significant ids in 124 terms
  # with redundancy between ids and terms
  # complete map between greedy terms and ids
  # the same id can be found in multiple greedy terms
  greedy_sets_all <- all_sets %>% 
    dplyr::filter(term %in% unique(greedy_sets$term))

  all_by_greedy <- purrr::map(unique(out$term), ~ {
    curr_set <- all_sets %>% 
      dplyr::filter(term %in% .x)
    
    # distribution of curr_set under greedy_sets_all
    greedy_sets_all %>% 
      dplyr::filter(id %in% curr_set$id) %>% # limit the complete map to ids found in curr_set
      dplyr::group_by(term) %>%
      dplyr::summarise(n = n()) %>% # number of curr_set ids under greedy sets
      dplyr::mutate(percent = n/nrow(curr_set)) %>% 
      dplyr::mutate(curr_set = .x)
  }, all_sets, greedy_sets_all) %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(curr_set) %>%
    dplyr::rename(ess_term = "term", count = "n") %>%
    dplyr::rename(term = "curr_set") %>% 
    dplyr::select(c("ess_term", "term", "count", "percent")) %>% 
    dplyr::mutate(percent = round(percent, digits = 3)) %>% 
    dplyr::ungroup(ess_term)

  all_by_greedy <- all_by_greedy %>% 
    dplyr::filter(count > 10, percent >= .33) %>% 
    dplyr::mutate(term = factor(term, levels = unique(as.character(term)))) %>% 
    dplyr::mutate(ess_term = factor(ess_term, levels = levels(term))) %>%
    dplyr::mutate(source = as.numeric(ess_term), target = as.numeric(term)) 
}


#' Greedy set cover of GSPA
#'
#' @import purrr dplyr rlang RcppGreedySetCover
#' @importFrom magrittr %>%
essential_gsets_interim <- function (out, df, filepath, filename, formula, all_sets, greedy_sets) {
  # dir.create(file.path(filepath, formula, "essential_sets"), recursive = TRUE, showWarnings = FALSE)
  # fn_gsc_prefix <- paste0(gsub("\\.[^.]*$", "", filename), "_", formula)
  # out_nm_gsc <- paste0(fn_gsc_prefix, "_greedysets.csv")

  # instead: load file "essmap_Protein_GSPA_Z_W2_bat.csv"
  all_by_greedy <- map_essential(all_sets, greedy_sets, out)
  
  # heat map
  ess_vs_all <- all_by_greedy %>% 
    dplyr::select(-which(names(.) %in% c("count", "source", "target"))) %>% 
    tidyr::spread(term, percent) %>% 
    dplyr::mutate_at(vars(which(names(.) != "ess_term")), ~ {1 - .x}) %>% 
    tibble::column_to_rownames("ess_term")
  
  d_row <- dist(ess_vs_all)
  d_col <- dist(t(ess_vs_all))
  max_d_row <- max(d_row, na.rm = TRUE)
  max_d_col <- max(d_col, na.rm = TRUE)
  d_row[is.na(d_row)] <- 1.2 * max_d_row
  d_col[is.na(d_col)] <- 1.2 * max_d_col
  
  ph <- pheatmap(
    ess_vs_all,
    cluster_rows =  hclust(d_row), 
    cluster_cols = hclust(d_col),
    treeheight_col = 1.1 * max_d_col, 
    treeheight_row = 1.1 * max_d_row, 
    show_rownames = FALSE, 
    show_colnames = FALSE,
  )

  cluster <- data.frame(cluster = cutree(ph$tree_col, h = max_d_col)) %>% 
    tibble::rownames_to_column("term")

  # all_by_greedy %>% 
  #   tidyr::spread(ess_set, percent) %>% 
  #   write.csv(file.path(filepath, formula, "essential_sets", "gs_greedy_map.csv"), row.names = FALSE)
  
  # for java
  min_target <- min(all_by_greedy$target, na.rm = TRUE)
  all_by_greedy <- all_by_greedy %>% 
    dplyr::mutate(source = source - min_target, target = target - min_target)
  rm(min_target)
  
  all_by_greedy <- all_by_greedy %>% 
    dplyr::left_join(cluster, by = "term")

  my_links <- all_by_greedy %>% 
    dplyr::select(c("source", "target", "percent")) %>% 
    dplyr::mutate(percent = percent * 10)
  my_nodes <- all_by_greedy %>% dplyr::select(c("term", "cluster", "count"))

  networkD3::forceNetwork(Links = my_links, Nodes = my_nodes, Source = "source",
               Target = "target", Value = "percent", NodeID = "term", Nodesize = "count", 
               Group = "cluster", opacity = 0.8, zoom = TRUE) %>% 
    networkD3::saveNetwork(file = 'GSPA.html')
}









plot_GSPAHM <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  err_msg <- "Don't call the function with arguments `annot_cols` and/or `annot_colnames`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)
  
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  dir.create(file.path(dat_dir, "Protein\\GSPA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  proteoGSPA(id = gene, task = plothm, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
  browser()
  
  quietly_log <- purrr::quietly(proteoGSPA)(id = gene, task = plothm, 
                                           annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, 
                                           ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\GSPA\\log","plot_hmGSPA_log.csv"), append = TRUE)  
}



gspaHM <- function(id, label_scheme_sub, anal_type, scale_log2r, 
                        filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  fmls <- dots %>% .[grepl("~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  if(purrr::is_empty(fmls))
    stop("No formula(s) of contrasts available.", call. = TRUE)
  
  # need formula subdir
  # dots also contain annot_cols and annot_colnames
  formulas <- names(fmls)
  
  
  # df <- df %>% 
  #   filters_in_call(!!!filter_dots) %>% 
  #   arrangers_in_call(!!!arrange_dots)
  
  if (length(formulas) > 0) purrr::map(formulas, fml_gspahm, 
                                      id = !!id, filepath = filepath, filename = filename, 
                                      !!!dots)

  # ins <- list.files(path = filepath, pattern = "_r\\d+\\.rda$")

}



fml_gspahm <- function (formula, id, filepath, filename, ...) {

  # instead: load file "essmap_Protein_GSPA_Z_W2_bat.csv"
  # all_by_greedy <- map_essential(all_sets, greedy_sets, out)
  
  dir.create(file.path(filepath, formula), recursive = TRUE, showWarnings = FALSE)
  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- prep_gspa(df = df, formula = formula, col_ind = col_ind, 
                  pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff) 
  
  df <- df %>% 
    dplyr::mutate(p_val = -log10(p_val)) %>% 
    dplyr::mutate(valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>% 
    dplyr::mutate(valence = factor(valence, levels = c("neg", "pos"))) 
  
  # imputation of NA to 0 increases the number of entries in descriptive "mean" calculation, 
  # which is like a penalty function
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
  
  pass <- purrr::map(res[paste0("pVal (", contrast_groups, ")")], ~ .x <= pval_cutoff) %>%
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
  
  fn_prefix <- paste0(gsub("\\.[^.]*$", "", filename), "_", formula)
  out_nm <- paste0(fn_prefix, ".csv")
  
  out <- dplyr::bind_cols(pval, adjp, log2fc) %>% 
    dplyr::mutate_at(.vars = grep("^p_val$|^adjP$", names(.)), format, scientific = TRUE, digits = 2) %>%
    dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2)
  
  # find essential sets
  gsets <- gsets %>% 
    .[!grepl("molecular_function$", names(.))] %>% 
    .[!grepl("cellular_component$", names(.))] %>% 
    .[!grepl("biological_process$", names(.))]
  
  all_sets <- purrr::map2(gsets, names(gsets), ~ {
    .x %>% 
      tibble() %>% 
      dplyr::mutate(set = .y)
  }) %>% 
    bind_rows() %>% 
    `names<-` (c("id", "term")) %>% 
    dplyr::filter(id %in% unique(df$entrez)) %>% 
    dplyr::select(c("term", "id")) %>% 
    dplyr::filter(term %in% res_pass$term)
  
  # one id only assign to one term without overlap
  # 969 significant ids in 124 terms
  # non-redundant mapping of ids to greedy terms
  # no redundant ids
  greedy_sets <- all_sets %>% 
    RcppGreedySetCover::greedySetCover(FALSE) %>% 
    dplyr::select(term) %>% 
    dplyr::filter(!duplicated(term)) %>% 
    dplyr::mutate(is_essential = TRUE)
  
  out <- greedy_sets %>% dplyr::right_join(., out, by = "term") %T>% 
    write.csv(file.path(filepath, formula, out_nm), row.names = FALSE)
  
  map_essential(all_sets, greedy_sets, out) %>% 
    write.csv(file.path(filepath, formula, paste0("essmap_", fn_prefix, ".csv")), row.names = FALSE)
  
  invisible(out)
}



