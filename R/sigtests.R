#' Row Variance
#' 
#' @param x A data frame.
#' @param na.rm The same as in \code{mean}.
rowVars <- function (x, na.rm = TRUE) 
{
  n <- rowSums(!is.na(x))
  n[n <= 1] <- NA_real_
  
  rowSums((x - rowMeans(x, na.rm = na.rm))^2, na.rm = na.rm)/(n - 1)
}


#' Filters rows by variance quantiles
#' 
#' @inheritParams info_anal
#' @inheritParams gn_rollup
#' @inheritParams prnSig
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
filterData <- function (df, cols = NULL, var_cutoff = 1E-3) 
{
  if (is.null(cols)) 
    cols <- 1:ncol(df)
  
  if (length(cols) > 1L) {
    df <- df %>% 
      dplyr::select(cols) %>% 
      dplyr::mutate(Variance = rowVars(.), 
                    rowname = rownames(.)) 
    
    Quantile <- quantile(df$Variance, probs = var_cutoff, na.rm = TRUE)
    
    df <- df %>%
      dplyr::filter(Variance >= pmax(Quantile, 1E-3)) %>%
      tibble::remove_rownames() %>% 
      tibble::column_to_rownames() %>% 
      dplyr::select(-Variance)
  }
  
  invisible(df)
}


#' Prepare formulas
#' 
#' @inheritParams gspaTest
#' @inheritParams model_onechannel
#' @importFrom magrittr %>% %T>% %$% %<>% 
prepFml <- function(formula, label_scheme_sub, ...) 
{
  # formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]  + (1|TMT_Set) + (1|Duplicate)
  # formula = ~ Term["Ner-V", "Ner_PLUS_PD-PD", "(Ner_PLUS_PD-PD)-(Ner-V)"]
  # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|TMT_Set)
  # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|Duplicate)
  # formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]
  # formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|TMT_Set)
  # formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|Duplicate)
  # formula = ~ Term["Ner-PD", "V-PD", "(Ner_PLUS_PD-PD)"] + (1|Duplicate)
  # formula = log2Ratio ~ Term
  # formula = ~ Term[~V] # no interaction terms
  # formula = ~ Term
  
  dots <- rlang::enexprs(...)
  
  fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
  len <- length(fml)
  
  key_col <- fml[len] %>% gsub("(.*)\\[\\s*\\~*.*\\].*", "\\1", .)
  
  label_scheme_sub <- label_scheme_sub %>% 
    dplyr::filter(!is.na(!!sym(key_col)))
  
  # formula = ~ Term[ ~ V]
  base <- if (grepl("\\[\\s*\\~\\s*", fml[len])) 
    gsub(".*\\[\\s*\\~\\s*(.*)\\]", "\\1", fml[len])
  else 
    NULL
  
  if (!is.null(base)) { # formula = ~ Term[~V]
    new_levels <- label_scheme_sub[[key_col]] %>% levels()
    new_levels <- c(new_levels[new_levels == base], 
                    new_levels[new_levels != base])
    label_scheme_sub <- label_scheme_sub %>%
      dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col), levels = new_levels))
  } 
  else if (!grepl("\\[", fml[len])) { # formula = log2Ratio ~ Term
    new_levels <- label_scheme_sub[[key_col]] %>% levels() # leveled by the alphabetic order
  }
  else { # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]
    new_levels <- NULL
  }
  
  if (!is.null(new_levels)) {
    contrs <- paste(new_levels[-1], new_levels[1], sep = "-")
    elements <- new_levels
  } 
  else {
    contrs <- fml[len] %>%
      gsub(".*\\[(.*)\\].*", "\\1", .) %>%
      gsub("\\\"", "", .) %>%
      str_split(",\\s*", simplify = TRUE) %>%
      as.character()
    
    new_contrs <- fml[len] %>%
      gsub("^.*\\[(.*)\\].*", "\\1", .) %>% # may have random terms at the end
      gsub("\\\"", "", .) %>% 
      str_split(",\\s*", simplify = TRUE) %>% 
      gsub("\\s+", "", .) %>% 
      gsub("<([^>]*?)\\+([^>]*?)>", "<\\1.plus.\\2>", .) %>% 
      gsub("<([^>]*?)\\-([^>]*?)>", "<\\1.minus.\\2>", .) %>% 
      gsub("[ <>]+", "", .)
    
    new_elements <- new_contrs %>%
      gsub("/[0-9]", "", .) %>% # (A+B+C)/3-D
      gsub("[\\(\\)]", "", .) %>%
      str_split("[\\+\\-]\\s*", simplify = TRUE) %>%
      as.character() %>%
      unique() %>% 
      .[. != ""]
    
    elements <- new_elements %>% 
      gsub(".plus.", "+", ., fixed = TRUE) %>% 
      gsub(".minus.", "-", ., fixed = TRUE)
    
    message("\ncontrs: ", contrs %>% as.character, "\n")
    message("new_contrs: ", new_contrs %>% as.character, "\n")
    message("elements: ", paste(elements, collapse = ", "), "\n")
    message("new_elements: ", paste(new_elements, collapse = ", "), "\n\n")
  }
  
  label_scheme_sub_sub <- label_scheme_sub %>%
    dplyr::filter(!!sym(key_col) %in% elements) %>%
    dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col)))
  
  if (!nrow(label_scheme_sub_sub))
    stop("No samples were found for formula ", formula, 
         "\nCheck the terms under column ", key_col)
  
  design <- model.matrix(~0+label_scheme_sub_sub[[key_col]]) %>%
    `colnames<-`(levels(label_scheme_sub_sub[[key_col]]))
  
  new_design_nms <- colnames(design) %>% 
    gsub("+", ".plus.", ., fixed = TRUE) %>% 
    gsub("-", ".minus.", ., fixed = TRUE)
  
  new_design <- design %>% 
    `colnames<-`(new_design_nms)
  
  contr_mat <- makeContrasts(contrasts = new_contrs, 
                             levels = data.frame(new_design)) %>% 
    `colnames<-`(contrs) %>% 
    `rownames<-`(colnames(design))
  
  rm(new_design_nms, new_design)
  
  random_vars <- fml[len] %>%
    gsub("\\[.*\\]+?", "", .) %>%
    paste("~", .) %>%
    as.formula() %>%
    terms.formula(.) %>%
    attr(., "term.labels") %>%
    .[grepl("\\|", .)] %>%
    gsub("1\\s*\\|\\s*(.*)", "\\1", .)
  
  message("random_vars: ", as.character(random_vars), "\n\n")
  
  list(design = design, 
       contr_mat = contr_mat, 
       key_col = key_col, 
       random_vars = random_vars,
       elements = elements, 
       label_scheme_sub_sub = label_scheme_sub_sub)
}


#' Adjusted pVal
#' 
#' @param df_pval A data frame containing pVals
#' @inheritParams prnSig 
#' @importFrom magrittr %>% %T>% %$% %<>% 
my_padj <- function(df_pval, pval_cutoff) 
{
  df_pval %>%
    purrr::map(~ .x <= pval_cutoff)	%>%
    purrr::map(~ ifelse(!.x, NA, .x)) %>%
    purrr::map2(as.list(df_pval), `*`) %>%
    purrr::map(~ p.adjust(.x, "BH")) %>%
    dplyr::bind_cols() %>%
    `names<-`(gsub("pVal", "adjP", colnames(.))) %>%
    dplyr::mutate(rowname = rownames(df_pval)) %>%
    dplyr::bind_cols(df_pval, .) %>%
    dplyr::mutate_at(.vars = grep("pVal\\s+", names(.)), format, 
                     scientific = TRUE, digits = 2) %>%
    dplyr::mutate_at(.vars = grep("adjP\\s+", names(.)), format, 
                     scientific = TRUE, digits = 2) %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames()
}


#' Model summary
#' 
#' @param pvals A data frame of pVal
#' @param log2rs A data frame of log2FC
#' @inheritParams prnSig
#' @importFrom magrittr %>% %T>% %$% %<>% 
lm_summary <- function(pvals, log2rs, pval_cutoff, logFC_cutoff, 
                       padj_method = "BH") 
{
  nms <- rownames(pvals)
  
  pass_pvals <- purrr::map(pvals, ~ .x <= pval_cutoff)
  pass_fcs <- purrr::map(log2rs, ~ abs(.x) >= log2(logFC_cutoff))
  
  pass_both <- purrr::map2(pass_pvals, pass_fcs, `&`) %>% 
    purrr::map(~ ifelse(!.x, NA, .x))
  
  res_padj <- pvals %>%
    purrr::map2(pass_both, `*`) %>%
    purrr::map(~ p.adjust(.x, padj_method)) %>%
    data.frame(check.names = FALSE) %>%
    `names<-`(gsub("pVal", "adjP", colnames(.))) %>%
    `rownames<-`(nms) %>%
    tibble::rownames_to_column() %>% 
    dplyr::bind_cols(pvals, .) %>%
    dplyr::mutate_at(.vars = grep("pVal\\s+", names(.)), format, 
                     scientific = TRUE, digits = 2L) %>%
    dplyr::mutate_at(.vars = grep("adjP\\s+", names(.)), format, 
                     scientific = TRUE, digits = 2L) %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames()
  
  log2rs <- log2rs %>%
    to_linfc() %>%
    `colnames<-`(gsub("log2Ratio", "FC", names(.))) %>%
    dplyr::bind_cols(log2rs, .) %>%
    dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2L) %>%
    `rownames<-`(nms)
  
  cbind.data.frame(res_padj, log2rs)
}


#' Factorial model formula for interaction terms. All factors under one channel
#' in `label_scheme_sub`
#' 
#' @param formula Language; the formula in linear modeling.
#' @param dfI Intensity data frame.
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnSig
#' @importFrom MASS ginv
model_onechannel <- function (dfR = NULL, dfI = NULL, 
                              id, formula, label_scheme_sub, 
                              complete_cases = FALSE, impute_group_na = TRUE, 
                              method = "limma", padj_method = "BH", 
                              var_cutoff = 1E-3, pval_cutoff = 1.00, 
                              logFC_cutoff = log2(1), ...) 
{
  # formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]  + (1|TMT_Set) + (1|Duplicate)
  # formula = ~ Term["Ner-V", "Ner_PLUS_PD-PD", "(Ner_PLUS_PD-PD)-(Ner-V)"]
  # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|TMT_Set)
  # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|Duplicate)
  # formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]
  # formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|TMT_Set)
  # formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|Duplicate)
  # formula = ~ Term["Ner-PD", "V-PD", "(Ner_PLUS_PD-PD)"] + (1|Duplicate)
  # formula = log2Ratio ~ Term
  # formula = ~ Term[~V] # no interaction terms
  # formula = ~ Term
  
  options(warn = 1L)
  
  dots <- rlang::enexprs(...)
  # lmFit_dots <- dots %>% .[. %in% c("method")]
  # eBayes_dots <- dots %>% .[. %in% c("proportion")]
  
  fml_ops     <- prepFml(formula, label_scheme_sub, ...)
  elements    <- fml_ops$elements
  contr_mat   <- fml_ops$contr_mat
  design      <- fml_ops$design
  key_col     <- fml_ops$key_col
  random_vars <- fml_ops$random_vars
  label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub
  
  local({
    col_nms <- colnames(contr_mat)
    if (length(dups <- which(duplicated(col_nms)))) {
      stop("Duplicated contrasts found: ", 
           paste(col_nms[dups], collapse = ", "))
    }
    
    ss <- (colSums(design) == 1)
    if (length(bads <- ss[ss])) 
      warning("Single sample condition: ", paste0(names(bads), collapse = ", "), 
              " under `", formula, "`.")
  })
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  # keep the name list as rows may drop in filtration
  df_nms <- tibble::tibble(!!id := rownames(dfR))
  
  if (complete_cases) {
    dfR <- dfR[complete.cases(dfR), ]
  }

  sids <- label_scheme_sub_sub[["Sample_ID"]]
  dfR  <- filterData(dfR, var_cutoff = var_cutoff)
  dfR  <- dfR[, sids, drop = FALSE]
  dfI  <- dfI[, sids, drop = FALSE]
  dfI  <- dfI[rownames(dfI) %in% rownames(dfR), ]

  stopifnot(identical(colnames(dfI), colnames(dfR)))
  stopifnot(identical(nrow(dfI), nrow(dfR)))
  
  ###
  # Haven't yet test within-group imputation at mixed effect modeling...
  ###
  
  if (impute_group_na && !is.null(dfI)) {
    ans <- base_sigtest_y(
      dfR = dfR, dfI = dfI, elements = elements, key_col = key_col, 
      label_scheme_sub_sub = label_scheme_sub_sub, seed = 1234L)
    dfsR <- ans[["log2R"]]
    dfsI <- ans[["Intensity"]]
    ys_base <- ans[["baseline"]]
    rm(list = "ans")

    mapply(function (x, y) {
      stopifnot(identical(rownames(x), rownames(y)))
      stopifnot(identical(colnames(x), colnames(y)))
    }, dfsR, dfsI)

    ans2 <- impute_baseline_ints(
      dfsR = dfsR, dfsI = dfsI, ys_base = ys_base, 
      sample_ids = sids)
    
    dfR <- dplyr::bind_cols(lapply(ans2, `[[`, "log2R"))
    dfR <- dfR[, sids, drop = FALSE]
    # dfI <- dplyr::bind_cols(lapply(ans2, `[[`, "Intensity"))
    # dfI <- dfI[, sids, drop = FALSE]
    rm(list = "ans2")
  }

  local({
    ncol <- ncol(dfR)
    
    if (ncol < 4L) 
      warning(formula, ": the total number of samples is ", ncol, ".\n", 
              "May need more samples for statistical tests.")
  })
  
  if (length(random_vars)) {
    if (length(random_vars) > 1L) {
      warning("Uses only the first random variable: ", random_vars[1])
    }
    
    design_random <- label_scheme_sub_sub[[random_vars[1]]] 
    corfit <- duplicateCorrelation(dfR, design = design, block = design_random)
    
    fit <- suppressWarnings(
      dfR %>%
        lmFit(design = design, 
              block = design_random, 
              correlation = corfit$consensus) %>%
        contrasts.fit(contr_mat) %>%
        eBayes()
    )
  } 
  else {
    fit <- suppressWarnings(
      lmFit(dfR, design = design) |>
        contrasts.fit(contr_mat) |>
        eBayes()
    )
  }
  
  print(design)
  print(contr_mat)
  
  # limma
  log2rs <- fit$coefficients %>%
    data.frame(check.names = FALSE) %>%
    `names<-`(paste0("log2Ratio (", names(.), ")"))
  
  pvals <- fit$p.value %>%
    data.frame(check.names = FALSE) %>%
    `names<-`(paste0("pVal (", names(.), ")"))
  
  res_lm <- lm_summary(pvals, log2rs, pval_cutoff, logFC_cutoff, padj_method)
  
  if (method %in% c("lmer", "lme", "lm")) {
    fml_rhs <- gsub("\\[.*\\]", "", formula) %>% .[length(.)]
    new_formula <- as.formula(paste("log2Ratio", "~", fml_rhs))
    
    contr_mat_lm <- t(contr_mat) %>%
      MASS::ginv() %>%
      `colnames<-`(colnames(contr_mat)) %>%
      `rownames<-`(rownames(contr_mat))
    
    names(dimnames(contr_mat_lm)) <- c("Levels", "Contrasts")
    contr_mat_lm <- list(Cdn = contr_mat_lm) %>% `names<-`(key_col)
    
    smpl_levels <- names(dfR)
    contr_levels <- attributes(contr_mat_lm[[key_col]])$dimnames$Contrasts
    
    df_lm <- dfR %>%
      tibble::rownames_to_column(id) %>%
      tidyr::gather(-id, key = Sample_ID, value = log2Ratio) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = smpl_levels)) %>%
      dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", key_col, random_vars)], 
                       by = "Sample_ID") %>%
      dplyr::select(which(not_all_NA(.))) %>%
      dplyr::group_by(!!rlang::sym(id)) %>%
      tidyr::nest()
    
    if (length(random_vars)) {
      if (!requireNamespace("broom.mixed", quietly = TRUE)) {
        stop("\n============================================================", 
             "\nNeed package \"broom.mixed\" for this function to work.",
             "\n============================================================")
      }
      
      if (!requireNamespace("lmerTest", quietly = TRUE)) {
        stop("\n============================================================", 
             "\nNeed package \"lmerTest\" for this function to work.",
             "\n============================================================")
      }
      
      res_lm <- df_lm %>%
        dplyr::mutate(
          model = purrr::map(data, 
                             ~ lmerTest::lmer(data = .x, formula = new_formula, 
                                              contrasts = contr_mat_lm))) %>%
        dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>%
        tidyr::unnest(glance, keep_empty = TRUE) %>% 
        dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>%
        dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>%
        dplyr::mutate(term = gsub(key_col, "", term)) %>% 
        dplyr::select(-data, -model) %>% 
        dplyr::mutate(term = factor(term, levels = contr_levels)) %>%
        tidyr::spread(term , p.value) %>%
        `rownames<-`(NULL) %>% 
        tibble::column_to_rownames(id) %>%
        `names<-`(paste0("pVal (", names(.), ")")) %>%
        lm_summary(log2rs, pval_cutoff, logFC_cutoff, padj_method)
    } 
    else {
      if (!requireNamespace("broom", quietly = TRUE)) {
        stop("\n============================================================", 
             "\nNeed package \"broom\" for this function to work.",
             "\n============================================================")
      }
      
      res_lm <- df_lm %>%
        dplyr::mutate(model = purrr::map(data, ~ lm(data = .x, formula = new_formula,
                                                    contrasts = contr_mat_lm))) %>%
        dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>%
        tidyr::unnest(glance, keep_empty = TRUE) %>%	
        dplyr::filter(!grepl("Intercept", term)) %>%
        dplyr::select(-c("std.error", "estimate", "statistic")) %>% 
        dplyr::mutate(term = gsub(key_col, "", term)) %>% 
        dplyr::select(-data, -model) %>% 
        dplyr::mutate(term = factor(term, levels = contr_levels))	%>% 
        tidyr::spread(term , p.value)	%>% 
        `rownames<-`(NULL) %>% 
        tibble::column_to_rownames(id) %>%
        `names<-`(paste0("pVal (", names(.), ")"))  %>% 
        lm_summary(log2rs, pval_cutoff, logFC_cutoff, padj_method)
    }
  }
  
  bads <- df_nms[!df_nms[[id]] %in% rownames(res_lm), ] |>
    tibble::column_to_rownames(var = id)
  
  if (FALSE) {
    df_op <- res_lm %>% 
      tibble::rownames_to_column(id) %>%
      dplyr::right_join(df_nms, by = id) %>%
      `rownames<-`(NULL) %>% 
      tibble::column_to_rownames(var = id)
  }
  
  df_op <- dplyr::bind_rows(res_lm, bads)
  
  df_op <- df_op %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("^(pVal \\(|adjP \\()"),
        function (x) tidyr::replace_na(as.numeric(x), 1)
      )
    )
  
  df_op <- df_op %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("^log2Ratio \\(", ),
        function (x) tidyr::replace_na(x, 0.0)
      )
    )
  
  df_op <- df_op %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("^FC \\(", ),
        function (x) tidyr::replace_na(x, 1.0)
      )
    )
}


#' Significance tests
#' 
#' @param data_type The type of data being either \code{Peptide} or \code{Protein}.
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnSig
#' @import limma stringr purrr tidyr dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
sigTest <- function(df, id, label_scheme_sub, 
                    scale_log2r = TRUE, complete_cases = FALSE, 
                    impute_na = FALSE, impute_group_na = TRUE, 
                    rm_allna = FALSE, method_replace_na, 
                    filepath, filename, 
                    method, padj_method, var_cutoff, pval_cutoff, logFC_cutoff, 
                    data_type, anal_type, ...) 
{
  dat_dir <- get_gl_dat_dir()
  
  stopifnot(vapply(c(var_cutoff, pval_cutoff, logFC_cutoff), is.numeric, 
                   logical(1L)))
  
  id     <- rlang::as_string(rlang::enexpr(id))
  method <- rlang::as_string(rlang::enexpr(method))
  dots   <- rlang::enexprs(...)
  
  lang_dots    <- dots[unlist(lapply(dots, is.language))]
  filter_dots  <- lang_dots[grepl("^filter_", names(lang_dots))]
  arrange_dots <- lang_dots[grepl("^arrange_", names(lang_dots))]
  dots         <- dots[!dots %in% c(filter_dots, arrange_dots)]
  non_fml_dots <- dots[!purrr::map_lgl(dots, is_formula)]
  dots         <- dots[purrr::map_lgl(dots, is_formula)]
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    pepSig_formulas <- dots
    save(pepSig_formulas, file = file.path(dat_dir, "Calls/pepSig_formulas.rda"))
    rm(list = "pepSig_formulas")
  } 
  else if (id %in% c("prot_acc", "gene")) {
    if (length(dots)) {
      prnSig_formulas <- dots
    } 
    else {
      prnSig_formulas <- dots <- concat_fml_dots()
    }
    
    save(prnSig_formulas, file = file.path(dat_dir, "Calls", "prnSig_formulas.rda"))
    rm(list = "prnSig_formulas")
  }	
  
  fn_prefix2 <- if (impute_na) "_impNA_pVals.txt" else "_pVals.txt"
  
  tempdata <- df |>
    filters_in_call(!!!filter_dots) |>
    arrangers_in_call(!!!arrange_dots) |>
    prepDM(id = !!id, 
           scale_log2r = scale_log2r, 
           sub_grp = label_scheme_sub$Sample_ID, 
           anal_type = anal_type, 
           rm_allna = rm_allna)
  dfR <- tempdata[["log2R"]]
  dfI <- tempdata[["Intensity"]]
  
  # in case of all-NA sample columns being removed
  label_scheme_sub <- label_scheme_sub |> 
    dplyr::filter(Sample_ID %in% colnames(dfR))
  
  dfR <- local({
    if ((!impute_na) && method_replace_na == "min") {
      row_mins <- apply(dfR, 1, FUN = min, na.rm = TRUE)
      row_sds  <- apply(dfR, 1, FUN = sd,  na.rm = TRUE)
      row_sds[is.na(row_sds)] <- 5E-3
      n_nas <- rowSums(is.na(dfR))
      
      for (i in seq_along(row_mins)) {
        n_i <- n_nas[[i]]
        
        if (n_i > 0) {
          vals <- rnorm(n_i, mean = row_mins[[i]] - 1, sd = row_sds[[i]]/5)
          x <- dfR[i, ]
          x[is.na(x)] <- vals
          dfR[i, ] <- x
        }
      }
    }
    
    dfR
  })

  # `complete_cases` depends on lm contrasts
  df_op <- lapply(dots, function (formula) model_onechannel(
    dfR = dfR, dfI = dfI, id = !!id, formula = formula, 
    label_scheme_sub = label_scheme_sub, complete_cases = complete_cases, 
    impute_group_na = impute_group_na, method = method, 
    padj_method = padj_method, var_cutoff = var_cutoff, 
    pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff, 
    !!!non_fml_dots))
  df_op <- do.call("cbind", df_op)
  
  # record the `scale_log2r` status; otherwise, need to indicate it in a way
  # for example, `_N` or `_Z` in file names
  local({
    dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)	  
    
    type <- if (data_type == "Peptide")
      "pep"
    else if (data_type == "Protein")
      "prn"
    else
      stop("`data_type` needs to be either `Peptide` or `Protein`.")
    
    file <- paste0(type, "Sig_imp", ifelse(impute_na, "TRUE", "FALSE"), ".rda")
    
    call_pars <- c(scale_log2r = scale_log2r, 
                   complete_cases = complete_cases, 
                   impute_na = impute_na) %>% 
      as.list()
    
    save(call_pars, file = file.path(dat_dir, "Calls", file))
  })
  
  suppressWarnings(
    df_op <- df_op %>%
      tibble::rownames_to_column(id) %>% 
      dplyr::mutate(!!id := forcats::fct_explicit_na(!!rlang::sym(id))) %>% 
      dplyr::right_join(df, ., by = id) %T>% 
      write.table(file.path(filepath, paste0(data_type, fn_prefix2)), sep = "\t",
                  col.names = TRUE, row.names = FALSE)	  
  )
  
  wb <- createWorkbook("proteoQ")
  addWorksheet(wb, sheetName = "Results")
  openxlsx::writeData(wb, sheet = "Results", df_op)
  saveWorkbook(wb, file = file.path(filepath, paste0(data_type, "_pVals.xlsx")), 
               overwrite = TRUE) 
  
  invisible(df_op)
}


#' Impute baseline intensities.
#'
#' The values of log2Ratios are derived accordingly
#'
#' @param dfsR List of log2Ratio data frames. Data are grouped according to the
#'   formulas in \link{prnSig}.
#' @param dfsI The corresponding intensity data frames.
#' @param ys_base The baseline intensities. The length is equal to the number of
#'   columns in \code{dfR}.
#' @param sample_ids The sample IDs for joining data at the original order.
#' @param max_fold The maximum fold change allowed.
#' @importFrom magrittr %>%
impute_baseline_ints <- function (dfsR, dfsI, ys_base, sample_ids, 
                                  max_fold = 100)
{
  grps <- names(dfsR)

  ans  <- mapply(function (dfR, dfI, g, b) {
    # Complementary columns
    oks  <- !grps %in% g
    dfRc <- dfsR[oks]
    dfIc <- dfsI[oks]
    dfRc <- dplyr::bind_cols(dfRc)
    dfIc <- dplyr::bind_cols(dfIc)
    
    if (sum(nas <- rowSums(is.na(dfR)) == ncol(dfR))) {
      ybars <- rowMeans(dfIc, na.rm = TRUE)
      
      for (i in seq_along(dfI)) {
        dfI[nas, i] <- b[[i]]
      }
      
      dfR[nas, ] <- dfR_emp <- log2(dfI[nas, ] / ybars[nas])
    }
    
    # Lower bounds
    min_log2_fold <- -log2(max_fold)
    dfR_emp[dfR_emp <= min_log2_fold] <- min_log2_fold
    dfR[nas, ] <- dfR_emp
    
    list(log2R = dfR, Intensity = dfI)
  }, dfsR, dfsI, grps, ys_base, 
  SIMPLIFY = FALSE)
  
  ans
}


#' Calculate baseline intensities.
#'
#' @param dfR The complete data frame of log2Ratios.
#' @param dfI The corresponding intensity data frame.
#' @param ys_base The baseline intensities. The length is equal to the number of
#'   columns in \code{dfR}.
#' @param elements The elements of groups in contrast fits.
#' @param key_col The key column in contrast fits, e.g., column \code{Term}.
#' @param label_scheme_sub_sub The metadata corresponding \code{dfR}.
#' @param seed A seed for reproducible random number generations.
#' @importFrom magrittr %>%
base_sigtest_y <- function(dfR, dfI, elements = NULL, key_col = NULL, 
                           label_scheme_sub_sub, seed = NULL) 
{
  sids <- colnames(dfR)
  rnms <- rownames(dfR)
  
  gl_min <- max(
    min(dfI, na.rm = TRUE),  # * 5,
    quantile(dfI, .0001, na.rm = TRUE)
  )
  log10_gl_min <- log10(gl_min)
  
  ans <- lapply(elements, function (element) {
    sids_sub <- label_scheme_sub_sub |>
      dplyr::filter(!!rlang::sym(key_col) %in% element) |>
      dplyr::select(Sample_ID) |>
      unlist(use.names = FALSE)
    
    oks <- sids %in% sids_sub
    dfR_sub <- dfR[, oks, drop = FALSE]
    dfI_sub <- dfI[, oks, drop = FALSE]
    
    list(log2R = dfR_sub, Intensity = dfI_sub)
  })
  names(ans) <- elements
  dfsR <- lapply(ans, `[[`, "log2R")
  dfsI <- lapply(ans, `[[`, "Intensity")
  n_samples <-lengths(dfsR)
  rm(list = "ans")
  
  ## (1) Generate random, baseline intensities for each contrast group
  slope <- mapply(function (x, y) {
    coef(lm(log10(x) ~ rowMeans(log10(y), na.rm = TRUE)))[[2]]
  }, lapply(dfsI, rowVars), dfsI) |>
    mean()
  ys_base <- lapply(n_samples, gen_randoms, mu = log10_gl_min, sigma = slope, 
                    seed = seed)
  ys_base <- lapply(ys_base, function (x) 10^x + gl_min)
  
  list(baseline = ys_base, log2R = dfsR, Intensity = dfsI)
}


#' Generate random intensity values.
#' 
#' @param n The number of random values.
#' @param mu The mean.
#' @param sigma The standard deviation.
#' @param seed A seed for reproducible random number generations.
#' @importFrom magrittr %>%
gen_randoms <- function (n = 3L, mu = 0, sigma = 1.8, seed = NULL)
{
  if (!is.null(seed)) set.seed(seed)
  
  x <- NULL
  
  while(length(x) < n) {
    x <- c(x, rnorm(n - length(x), mean = mu, sd = sigma))
    x <- x[abs(x - mu) <= 2 * sigma]
  }
  
  x
}


#' Significance tests of peptide/protein log2FC
#' 
#' \code{pepSig} performs significance tests against peptide log2FC. 
#'
#' @rdname prnSig
#' 
#' @import purrr
#' @export
pepSig <- function (scale_log2r = TRUE, impute_na = FALSE, 
                    impute_group_na = TRUE, complete_cases = FALSE, 
                    rm_allna = FALSE, method = c("limma", "lm"), 
                    padj_method = "BH", method_replace_na = c("none", "min"), 
                    var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), 
                    df = NULL, filepath = NULL, filename = NULL, ...) 
{
  on.exit({
    mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
      c(rlang::enexprs(...)) %>% 
      save_call("pepSig")
  }, add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- tryCatch(
    match_call_arg(normPSM, group_psm_by), error = function(e) "pep_seq_mod")

  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  method <- rlang::enexpr(method)
  
  method <- if (method == rlang::expr(c("limma", "lm")))
    "limma"
  else
    rlang::as_string(method)
  
  # replace NA other than mice imputation
  # (only for sigTest)
  method_replace_na <- rlang::enexpr(method_replace_na)
  
  method_replace_na <- if (method_replace_na == rlang::expr(c("none", "min")))
    "none"
  else
    rlang::as_string(method_replace_na)
  
  stopifnot(method_replace_na %in% c("none", "min"), 
            length(method_replace_na) == 1L)
  
  stopifnot(method %in% c("limma", "lm"), length(method) == 1L)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  padj_method <- rlang::as_string(rlang::enexpr(padj_method))

  reload_expts()
  
  if ((!impute_na) && (method != "limma")) {
    impute_na <- TRUE
    warning("Coerce `impute_na = ", impute_na, "` at method = ", method)
  }
  
  stopifnot(vapply(c(scale_log2r, impute_na, complete_cases, rm_allna), 
                   rlang::is_logical, logical(1L)))
  
  stopifnot(vapply(c(var_cutoff, pval_cutoff, logFC_cutoff), 
                   is.numeric, logical(1L)))
  
  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            impute_group_na = impute_group_na, 
            method_replace_na = method_replace_na, 
            filepath = !!filepath, 
            filename = !!filename, 
            anal_type = "Model")(method = method, 
                                 padj_method = padj_method, 
                                 var_cutoff = var_cutoff, 
                                 pval_cutoff = pval_cutoff, 
                                 logFC_cutoff = logFC_cutoff, 
                                 rm_allna = rm_allna, 
                                 ...)
}


#'Significance tests of peptide/protein log2FC.
#'
#'\code{prnSig} performs significance tests against protein log2FC.
#'
#'In general, special characters of \code{+} or \code{-} should be avoided from
#'contrast terms. Occasionally, such as in biological studies, it may be
#'convenient to use \code{A+B} to denote a condition of combined treatment of
#'\code{A} and \code{B} . In the case, one can put the term(s) containing
#'\code{+} or \code{-} into a pair of pointy brackets. The syntax in the
#'following hypothetical example will compare the effects of \code{A}, \code{B},
#'\code{A+B} and the average of \code{A} and \code{B} to control \code{C}:
#'
#'\code{prnSig(fml = ~ Term["A - C", "B - C", "<A + B> - C", "(A + B)/2 - C"])}
#'
#'Note that \code{<A + B>} stands for one sample and \code{(A + B)} has two
#'samples in it.
#'
#'@inheritParams prnHist
#'@inheritParams prnHM
#'@inheritParams normPSM
#'@param filename A file name to output results. The default is
#'  \code{Peptide_pVals.txt} for peptides and \code{Protein_pVals} for proteins.
#'@param method Character string; the method of linear modeling. The default is
#'  \code{limma}. At \code{method = lm}, the \code{lm()} in base R will be used
#'  for models without random effects and the \code{\link[lmerTest]{lmer}} will
#'  be used for models with random effects.
#'@param impute_group_na Logical; if TRUE, impute (intensity) values that are
#'  exclusively NA under a sample group. The sample grouping is defined under
#'  \code{label_scheme}.
#'@param padj_method Character string; the method of multiple-test corrections
#'  for uses with \link[stats]{p.adjust}. The default is "BH". See
#'  ?p.adjust.methods for additional choices.
#'@param var_cutoff Numeric; the cut-off in the variances of \code{log2FC}.
#'  Entries with variances smaller than the threshold will be excluded from
#'  linear modeling. The default is 1E-3.
#'@param pval_cutoff Numeric; the cut-off in significance \code{pVal}. Entries
#'  with \code{pVals} smaller than the threshold will be excluded from multiple
#'  test corrections. The default is at \code{1} to include all entries.
#'@param logFC_cutoff Numeric; the cut-off in \code{log2FC}. Entries with
#'  absolute \code{log2FC} smaller than the threshold will be removed from
#'  multiple test corrections. The default is at \code{log2(1)} to include all
#'  entries.
#'@param method_replace_na The method to replace NA values by rows. The default
#'  is \code{none} by doing nothing. At \code{method_replace_na = min}, the row
#'  minimums will be used. The argument is only a device to assess \code{pVals},
#'  e.g., by handling the circumstance of all NA values under one group and
#'  non-trivial values under another. The setting of \code{min} might be useful
#'  at the experimenters' discretion of ascribing \code{NA} values to the lack
#'  of signals. The argument has no effects with \code{impute_na = TRUE}.
#'@param ... User-defined formulas for linear modeling. The syntax starts with a
#'  tilde, followed by the name of an available column key in
#'  \code{expt_smry.xlsx} and square brackets. The contrast groups are then
#'  quoted with one to multiple contrast groups separated by commas. The
#'   default column key is \code{Term} in \code{expt_smry.xlsx}: \cr \code{~
#'   Term["A - C", "B - C"]}. \cr \cr Additive random effects are indicated by
#'  \code{+ (1|col_key_1) + (1|col_key_2)}... Currently only a syntax of single
#'   contrast are supported for uses with random effects: \cr \code{~ Term["A -
#'   C"] + (1|col_key_1) + (1|col_key_2)} \cr \cr \code{filter_}: Logical
#'  expression(s) for the row filtration against data in a primary file of
#'  \code{Peptide[_impNA].txt} or \code{Protein[_impNA].txt}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements.
#'@return The primary output is \code{.../Peptide/Model/Peptide_pVals.txt} for
#'  peptide data or \code{.../Protein/Model/Protein_pVals.txt} for protein data.
#'  At \code{impute_na = TRUE}, the corresponding outputs are
#'  \code{Peptide_impNA_pvals.txt} or \code{Protein_impNA_pvals.txt}.
#'
#'@example inst/extdata/examples/prnSig_.R
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
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>%
#'
#'@export
prnSig <- function (scale_log2r = TRUE, impute_na = FALSE, 
                    impute_group_na = TRUE, complete_cases = FALSE, 
                    rm_allna = FALSE, method = c("limma", "lm"), 
                    padj_method = "BH", method_replace_na = c("none", "min"), 
                    var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), 
                    df = NULL, filepath = NULL, filename = NULL, ...) 
{
  on.exit({
    load(file.path(get_gl_dat_dir(), "Calls/prnSig_formulas.rda"))
    dots <- my_union(rlang::enexprs(...), prnSig_formulas)
    mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
      c(dots) %>% 
      save_call("prnSig")
  }, add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- tryCatch(
    match_call_arg(normPSM, group_pep_by), error = function(e) "gene")
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  method <- rlang::enexpr(method)
  method <- if (method == rlang::expr(c("limma", "lm"))) {
    "limma"
  } else {
    rlang::as_string(method)
  }
  stopifnot(method %in% c("limma", "lm"), length(method) == 1L)
  
  # replace NA other than mice imputation
  # (only for sigTest)
  method_replace_na <- rlang::enexpr(method_replace_na)
  method_replace_na <- if (method_replace_na == rlang::expr(c("none", "min"))) {
    "none"
  } else {
    rlang::as_string(method_replace_na)
  }

  stopifnot(method_replace_na %in% c("none", "min"), 
            length(method_replace_na) == 1L)

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  if ((!impute_na) && (method != "limma")) {
    impute_na <- TRUE
    warning("Coerce `impute_na = ", impute_na, "` at method = ", method)
  }
  
  stopifnot(vapply(c(scale_log2r, impute_na, complete_cases, rm_allna), 
                   rlang::is_logical, logical(1L)))
  
  stopifnot(vapply(c(var_cutoff, pval_cutoff, logFC_cutoff), 
                   is.numeric, logical(1L)))
  
  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            impute_group_na = impute_group_na, 
            method_replace_na = method_replace_na, 
            filepath = !!filepath, 
            filename = !!filename, 
            anal_type = "Model")(method = method, 
                                 padj_method = padj_method, 
                                 var_cutoff = var_cutoff, 
                                 pval_cutoff = pval_cutoff, 
                                 logFC_cutoff = logFC_cutoff, 
                                 rm_allna = rm_allna, 
                                 ...)
}


