#'Make GSEA gct
make_gct <- function(df, filepath, fn_prefix) {
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  
  df <- df %>% 
    tibble::rownames_to_column("NAME") %>% 
    dplyr::mutate(Description = NA) 
  
  df <- bind_cols(
    df %>% dplyr::select(NAME, Description), 
    df %>% dplyr::select(-NAME, -Description), 
  )
  
  hr_l1 <- "#1.2\n"
  hr_l2 <- paste0(c(nrow(df), ncol(df) - 2), collapse = "\t") %>% 
    paste0("\n")
  
  hr_nms <- paste0(names(df), collapse = "\t") %>% 
    paste0("\n")
  
  header <- paste0(hr_l1, hr_l2, hr_nms)
  file <- file.path(filepath, paste0(fn_prefix, ".gct"))
  cat(header, file = file)
  df %>% write.table(file, append = TRUE, sep = '\t', na = "", 
                     col.names = FALSE, row.names = FALSE, quote = FALSE)
}


#'Make GSEA cls 
make_cls <- function(df, nms, filepath, fn_prefix) {
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  
  cls_l1 <- paste0(c(length(nms), length(unique(nms)), "1"), collapse = "\t") %>% 
    paste0("\n")
  
  grps <- nms %>% 
    as.character() %>% 
    unique() %>% 
    paste0(collapse = "\t") %>% 
    paste0("#", ., "\n")
  
  smpls <- nms %>% 
    as.character() %>% 
    paste0(collapse = "\t") %>% 
    paste0("\n")
  
  cls <- paste0(cls_l1, grps, smpls)
  fn_cls <- file.path(filepath, paste0(fn_prefix, ".cls"))
  cat(cls, file = fn_cls)  
}


#'GSEA of protein data
#'
#'\code{proteoGSEA} prepares data for the analysis of
#'\code{\href{http://software.broadinstitute.org/gsea/index.jsp}{GSEA}} aganist
#'protein \code{log2FC} data. Users should avoid calling the method directly,
#'but instead use the following wrappers.
#'
#'The arguments \code{var_cutoff}, \code{pval_cutoff} and \code{logFC_cutoff}
#'are used to filter out low influence genes. Additional subsetting of data via
#'the \code{vararg} approach of \code{filter_} is feasible.
#'
#'The outputs include \code{Protein_GSEA.gct} and \code{protein_GSEA.cls} for
#'samples indicated in file \code{Protein_pVals.txt} or
#'\code{Protein_impNA_pVals.txt}. These outputs can be used with online
#'\code{\href{http://software.broadinstitute.org/gsea/index.jsp}{GSEA}}.
#'
#'The current GSEA may not support the comparisons between two grouped
#'conditions, i.e.,  (grpA + grpB) versus (grpC + grpD). The \code{prnGSEA}
#'utility further breaks the input data into pairs of groups according to the
#'formulae and contrasts defined in \code{pepSig} or \code{prnSig}. The
#'phenotype labels are then reformed in reflection of the original group names,
#'weights and directions, i.e., \code{0.5xgrpA&0.5xgrpB	-0.5xgrpC&-0.5xgrpD}.
#'The corresponding \code{.gct} and \code{.cls} files can be used with the
#'online or the github version of R-GSEA.
#'
#'@inheritParams proteoGSPA
#'@inheritParams proteoSigtest
#'@inheritParams  proteoEucDist
#'@inheritParams  proteoHM
#'@inheritParams  info_anal
#'@param impute_na Logical. At the NULL default, the TRUE or FALSE will match
#'  the choice in \code{\link{pepSig}} for peptide and \code{\link{prnSig}} for
#'  protein data.
#'@param gset_nms Not currently used.
#'@param var_cutoff Numeric value or vector; the cut-off in the variance of
#'  protein \code{log2FC} across samples. Entries with \code{variances} less
#'  than the threshold will be removed for enrichment analysis. The default is
#'  0.5.
#'@param min_size Not currently used.
#'@param max_size Not currently used.
#'@param min_greedy_size Not currently used.
#'@param gspval_cutoff Not currently used.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}.
#'@import dplyr rlang ggplot2 networkD3
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnGSEA_.R
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data normalization \cr
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'@export
proteoGSEA <- function (id = gene, scale_log2r = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                        impute_na = NULL, complete_cases = FALSE, gset_nms = "go_sets", 
                        var_cutoff = .5, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                        gspval_cutoff = 5E-2, min_size = 10, max_size = Inf, min_greedy_size = 1, 
                        fml_nms = NULL, task = "anal", ...) {
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- rlang::enexpr(id)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  task <- rlang::enexpr(task)
  if (is.null(impute_na)) impute_na <- match_sigTest_imputena(as_string(id))

  stopifnot(is_logical(scale_log2r), is_logical(impute_na), is_logical(complete_cases))
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSEA")

  curr_call <- mget(names(formals()), rlang::current_env()) %>% 
    .[names(.) != "..."] %>% 
    c(dots) 
  
  if (task == "anal") {
    curr_call %>% save_call("prnGSEA")
  } else if (task == "plothm") {
    curr_call %>% save_call("prnGSEAHM")
  }
  
  # Sample selection criteria:
  #   !is_reference under "Reference"
  #   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
  info_anal(df = !!df, id = !!id, scale_log2r = scale_log2r, 
            filepath = !!filepath, filename = !!filename, impute_na = impute_na, 
            anal_type = "GSEA")(complete_cases = complete_cases, gset_nms = gset_nms, 
                                var_cutoff = var_cutoff, pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, 
                                gspval_cutoff = gspval_cutoff, 
                                min_size = min_size, max_size = max_size, min_greedy_size = min_greedy_size, 
                                task = !!task, !!!dots)
}


#'Protein GSEA
#'
#'\code{prnGSEA} is a wrapper of \code{\link{proteoGSEA}} for gene set
#'enrichment analysis.
#'
#'@rdname proteoGSEA
#'
#'@import purrr
#'@export
prnGSEA <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\GSEA\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_normPSM_protid()
  
  quietly_log <- purrr::quietly(proteoGSEA)(id = !!id, task = anal, ...)
  quietly_log$result <- NULL
  purrr::walk(quietly_log, write, file.path(dat_dir, "Protein\\GSEA\\log","prnGSEA_log.csv"), append = TRUE)
}


#'Protein GSEA by formula(e) in `pepSig` or `prnSig`
fml_gsea <- function (fml, fml_nm, var_cutoff, pval_cutoff, 
                      logFC_cutoff, gspval_cutoff, min_size, max_size, 
                      df, col_ind, id, gsets, label_scheme_sub, complete_cases, scale_log2r, 
                      filepath, filename, ...) {
  
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  df <- local({
    ids <- df %>% 
      prep_gspa(id = !!id, fml_nm = fml_nm, col_ind = col_ind, 
                pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff) %>% 
      dplyr::select(id) %>% 
      unlist() %>% 
      unique()
    
    df <- df %>% dplyr::filter(!!sym(id) %in% ids)
    stopifnot(nrow(df) > 0)

    ids <- df %>% 
      `rownames<-`(.[[id]]) %>% 
      filterData(cols = grep(NorZ_ratios, names(.)), var_cutoff = var_cutoff) %>% 
      tibble::rownames_to_column(id) %>% 
      dplyr::select(id) %>% 
      unlist() %>% 
      unique()
    
    df <- df %>% dplyr::filter(!!sym(id) %in% ids)
    
    if (complete_cases) {
      rows <- df %>% 
        .[, grep(NorZ_ratios, names(.)), drop = FALSE] %>% 
        `names<-`(gsub(paste0(NorZ_ratios, "[0-9]{3}[NC]*\\s+\\((.*)\\)$"), "\\1", names(.))) %>% 
        dplyr::select(as.character(label_scheme_sub$Sample_ID)) %>% 
        complete.cases()
      
      df <- df[rows, ]      
    }
    
    return(df)
  })
  
  if (id != "gene") {
    df <- df %>% gn_rollup(grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)))
    id <- "gene"
  }
  
  df_log2r <- df %>% 
    prepDM(id = !!id, 
           scale_log2r = scale_log2r, 
           sub_grp = label_scheme_sub$Sample_ID, 
           anal_type = "GSEA") %>% 
    .$log2R
  
  make_gct(df = df_log2r, filepath = filepath, fn_prefix = "protein_GSEA")
  make_cls(df = df_log2r, nms = label_scheme_sub$Group, filepath = filepath, fn_prefix = "protein_GSEA")
  
  fml_ops <- prepFml(fml, label_scheme_sub, ...)
  contr_mat <- fml_ops$contr_mat
  design <- fml_ops$design
  key_col <- fml_ops$key_col
  random_vars <- fml_ops$random_vars
  label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub
  
  df <- df %>% 
    `rownames<-`(NULL) %>% 
    tibble::column_to_rownames(id) %>% 
    `names<-`(gsub(paste0(NorZ_ratios, "[0-9]{3}[NC]*\\s+\\((.*)\\)$"), "\\1", names(.))) %>% 
    dplyr::select(as.character(label_scheme_sub_sub$Sample_ID))

  nms <- purrr::map(1:ncol(contr_mat), ~ {
    rows <- contr_mat[, .x, drop = FALSE]
    
    rows_pos <- rows[rows > 0, , drop = FALSE] 
    rows_pos[, 1] <- rows_pos[, 1] %>% purrr::map_dbl(~ round(.x, digits = 2))
    nms_pos <- paste0(rows_pos, "x", rownames(rows_pos)) %>% 
      reduce(paste, sep = "&")
    
    rows_neg <- rows[rows < 0, , drop = FALSE]
    rows_neg[, 1] <- rows_neg[, 1] %>% purrr::map_dbl(~ round(.x, digits = 2))
    nms_neg <- paste0(rows_neg, "x", rownames(rows_neg)) %>% 
      reduce(paste, sep = "&")
    
    data.frame(pos = nms_pos, neg = nms_neg)
  }) %>% 
    do.call(rbind, .) %>% 
    `rownames<-`(colnames(contr_mat))
  
  nms <- nms %>% 
    purrr::map(~ gsub("[\\\\\\/\\:\\*\\?\\'\\<\\>\\|]", ".", .x)) %>% 
    bind_rows()

  contr_identity <- ifelse(contr_mat == 0, 0, 1)
  
  n_col <- ncol(design)
  df_sub <- rep(list(NULL), n_col)
  
  n_col2 <- ncol(contr_identity)
  pos <- neg <- both <- rep(list(NULL), n_col2)
  
  for (i in 1:n_col) df_sub[[i]] <- df[, which(design[, i] == 1), drop = FALSE]
  
  for (i in 1:n_col2) {
    filepath_i <- file.path(filepath, fml_nm)

    both[[i]] <- map2(df_sub, contr_identity[, i], ~ .x * .y) %>% 
      do.call(cbind, .) %>% 
      .[, names(df), drop = FALSE]
    
    idx_pos <- design[, contr_mat[, i] > 0, drop = FALSE] %>% 
      rowSums()
    pos[[i]] <- both[[i]][, which(idx_pos == 1), drop = FALSE]
    
    idx_neg <- design[, contr_mat[, i] < 0, drop = FALSE] %>% 
      rowSums()
    neg[[i]] <- both[[i]][, which(idx_neg == 1), drop = FALSE]
    
    df_i <- cbind(pos[[i]], neg[[i]])
    out_nm <- paste0("fml-", fml_nm, "_contr-", i)
    make_gct(df_i, filepath_i, out_nm)

    nms_pos_i <- rep(nms$pos[i], ncol(pos[[i]])) %>% as.character()
    nms_neg_i <- rep(nms$neg[i], ncol(neg[[i]])) %>% as.character()
    nms_both <- c(nms_pos_i, nms_neg_i)
    make_cls(df = df_i, nms = nms_both, filepath = filepath_i, fn_prefix = out_nm)
    
    run_scripts <- FALSE
    if (run_scripts) {
      dir_out <- file.path(filepath_i, out_nm, "res")
      dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
      
      res <- try (
        GSEA::GSEA(
          input.ds = file.path(filepath_i, paste0(out_nm, ".gct")),
          input.cls = file.path(filepath_i, paste0(out_nm, ".cls")),
          input.chip = "NOCHIP", 
          gs.db = file.path("~\\proteoQ\\dbs\\gsea\\msig\\c2_cgp_hs.gmt"),
          collapse.dataset = FALSE,
          collapse.mode = "NOCOLLAPSE",
          output.directory = dir_out,
        ), silent = TRUE)
    }
    
  }    

}


