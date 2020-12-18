#' Make GSEA gct
#' 
#' @param fn_prefix The base name of a file.
#' @inheritParams info_anal
make_gct <- function(df, filepath, fn_prefix) {
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  
  df <- df %>% 
    tibble::rownames_to_column("NAME") %>% 
    dplyr::mutate(Description = NA) 
  
  df <- dplyr::bind_cols(
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
#' 
#' @param nms The names of sample groups.
#' @inheritParams info_anal
#' @inheritParams make_gct
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


#'Protein GSEA
#'
#'\code{prnGSEA} prepares data for the analysis of
#'\href{http://software.broadinstitute.org/gsea/index.jsp}{GSEA} against
#'protein \code{log2FC} data. 
#'
#'The arguments \code{var_cutoff}, \code{pval_cutoff} and \code{logFC_cutoff}
#'are used to filter out low influence genes. Additional subsetting of data via
#'the \code{vararg} approach of \code{filter_} is feasible.
#'
#'The outputs include \code{Protein_GSEA.gct} and \code{protein_GSEA.cls} for
#'samples indicated in file \code{Protein_pVals.txt} or
#'\code{Protein_impNA_pVals.txt}. These outputs can be used with online
#'\href{http://software.broadinstitute.org/gsea/index.jsp}{GSEA}.
#'
#'The current GSEA may not support the comparisons between two grouped
#'conditions, i.e.,  (grpA + grpB) versus (grpC + grpD). The \code{prnGSEA}
#'utility further breaks the input data into pairs of groups according to the
#'formulas and contrasts defined in \code{pepSig} or \code{prnSig}. The
#'phenotype labels are then reformed in reflection of the original group names,
#'weights and directions, i.e., \code{0.5xgrpA&0.5xgrpB	-0.5xgrpC&-0.5xgrpD}.
#'The corresponding \code{.gct} and \code{.cls} files can be used with the
#'online or the github version of R-GSEA.
#'
#'@inheritParams prnGSPA
#'@param gset_nms Not currently used (to be chosen by users during online GSEA).
#'@param var_cutoff Numeric value or vector; the cut-off in the variance of
#'  protein \code{log2FC} across samples. Entries with \code{variances} less
#'  than the threshold will be removed for enrichment analysis. The default is
#'  0.5.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. 
#'@import dplyr ggplot2 
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'
#'@example inst/extdata/examples/prnGSEA_.R
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
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'  
#'@export
prnGSEA <- function (gset_nms = "go_sets", 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     df = NULL, filepath = NULL, filename = NULL, 
                     var_cutoff = .5, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     fml_nms = NULL, ...) {

  on.exit(
    if (rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(dots) %>% 
        save_call("pepGSEA")
    } else if (rlang::as_string(id) %in% c("prot_acc", "gene")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(dots) %>% 
        save_call("prnGSEA")
    }, 
    add = TRUE
  )
  
  check_dots(c("id", "anal_type", "df2"), ...)

  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Protein/GSEA/log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  stopifnot(vapply(c(scale_log2r, complete_cases, impute_na), rlang::is_logical, logical(1)))
  stopifnot(vapply(c(var_cutoff, pval_cutoff, logFC_cutoff), rlang::is_double, logical(1)))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  dots <- rlang::enexprs(...)
  dots <- concat_fml_dots(
    fmls = dots %>% .[grepl("^\\s*~", .)], 
    fml_nms = fml_nms, 
    dots = dots %>% .[!grepl("^\\s*~", .)], 
    anal_type = "GSEA"
  )

  reload_expts()
  
  # Sample selection criteria:
  #   !is_reference under "Reference"
  #   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
  
  message("`prnGSEA` compiles data for online GSEA and argument `gset_nms` not currently used.")
  
  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            filepath = !!filepath, 
            filename = !!filename, 
            anal_type = "GSEA")(gset_nms = "go_sets", 
                                var_cutoff = var_cutoff, 
                                pval_cutoff = pval_cutoff, 
                                logFC_cutoff = logFC_cutoff, 
                                !!!dots)
}


#'Protein GSEA by formula(s) in `pepSig` or `prnSig`
#'
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams fml_gspa
#' @inheritParams gsVolcano
fml_gsea <- function (fml, fml_nm, var_cutoff, pval_cutoff, 
                      logFC_cutoff, gspval_cutoff, gslogFC_cutoff, 
                      min_size, max_size, 
                      df, col_ind, id, gsets, label_scheme_sub, 
                      complete_cases, scale_log2r, 
                      filepath, filename, ...) {
  
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  df <- local({
    ids <- df %>% 
      prep_gspa(id = !!id, 
                fml_nm = fml_nm, 
                col_ind = col_ind, 
                pval_cutoff = pval_cutoff, 
                logFC_cutoff = logFC_cutoff, 
                use_adjP = FALSE) %>% 
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
    
    # No `lm` so `complete_cases` applied to all samples in label_scheme_sub
    if (complete_cases) {
      df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
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
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  make_gct(df = df_log2r, 
           filepath = filepath, 
           fn_prefix = fn_prefix)
  make_cls(df = df_log2r, 
           nms = label_scheme_sub$Group, 
           filepath = filepath, 
           fn_prefix = fn_prefix)
  
  fml_ops <- suppressMessages(prepFml(fml, label_scheme_sub, ...))
  contr_mat <- fml_ops$contr_mat
  design <- fml_ops$design
  key_col <- fml_ops$key_col
  random_vars <- fml_ops$random_vars
  label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub
  
  df <- df %>% 
    `rownames<-`(NULL) %>% 
    tibble::column_to_rownames(id) %>% 
    `names<-`(gsub(paste0(NorZ_ratios, "[0-9]{3}[NC]*\\s+\\((.*)\\)$"), "\\1", 
                   names(.))) %>% 
    dplyr::select(as.character(label_scheme_sub_sub$Sample_ID))

  nms <- purrr::map(1:ncol(contr_mat), ~ {
    rows <- contr_mat[, .x, drop = FALSE]
    
    rows_pos <- rows[rows > 0, , drop = FALSE] 
    rows_pos[, 1] <- rows_pos[, 1] %>% purrr::map_dbl(~ round(.x, digits = 2))
    nms_pos <- paste0(rows_pos, "x", rownames(rows_pos)) %>% 
      purrr::reduce(paste, sep = "&")
    
    rows_neg <- rows[rows < 0, , drop = FALSE]
    rows_neg[, 1] <- rows_neg[, 1] %>% purrr::map_dbl(~ round(.x, digits = 2))
    nms_neg <- paste0(rows_neg, "x", rownames(rows_neg)) %>% 
      purrr::reduce(paste, sep = "&")
    
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
    make_gct(df_i, filepath_i, paste0(fn_prefix, "_", out_nm))

    nms_pos_i <- rep(nms$pos[i], ncol(pos[[i]])) %>% as.character()
    nms_neg_i <- rep(nms$neg[i], ncol(neg[[i]])) %>% as.character()
    nms_both <- c(nms_pos_i, nms_neg_i)
    make_cls(df = df_i, 
             nms = nms_both, 
             filepath = filepath_i, 
             fn_prefix = paste0(fn_prefix, "_", out_nm))
  }    
}


