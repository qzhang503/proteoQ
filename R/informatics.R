#' Function factories for informatic analysis
#'
#' \code{info_anal} produces functions for selected informatic analysis.
#'
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param anal_type Character string; the type of analysis that are preset for
#'   method dispatch in function factories. The value will be determined
#'   automatically. Exemplary values include \code{anal_type = c("PCA",
#'   "Corrplot", "EucDist", "GSPA", "Heatmap", "Histogram", "MDS", "Model",
#'   "NMF", "Purge", "Trend", "LDA", ...)}.
#' @param method_replace_na A method of NA value replacements.
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams prnMDS
#' @inheritParams anal_pepNMF
#' @inheritParams prnGSPAMap
#' @inheritParams prnCorr_logFC
#'
#' @return A function to the given \code{anal_type}.
#' @import dplyr ggplot2 pheatmap openxlsx
#' @importFrom magrittr %>% %T>% %$% %<>% 
info_anal <- function (id = gene, id_gspa = "entrez", 
                       col_select = NULL, col_group = NULL, col_order = NULL,
                       col_color = NULL, col_fill = NULL, col_shape = NULL, 
                       col_size = NULL, col_alpha = NULL,
                       color_brewer = NULL, fill_brewer = NULL,
                       size_manual = NULL, shape_manual = NULL, 
                       alpha_manual = NULL, col_benchmark = NULL,
                       scale_log2r = TRUE, complete_cases = FALSE, 
                       impute_na = FALSE, impute_group_na = FALSE, 
                       method_replace_na = "none", 
                       df = NULL, df2 = NULL, filepath = NULL, filename = NULL,
                       anal_type = c("Corrplot", "Heatmap", "Histogram", 
                                     "MA", "MDS", "Model",
                                     "NMF", "Trend")) 
{
  scipen <- if (anal_type %in% c("MDS", "Volcano", "mapGSPA")) 999 else 0
  
  old_opts <- options()
  
  options(
    scipen = scipen,
    warn = 1L
  )
  
  on.exit(options(old_opts), add = TRUE)
  
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  
  stopifnot(vapply(c(scale_log2r, complete_cases, impute_na), 
                   rlang::is_logical, logical(1L)))
  
  err_msg1  <- paste0("Column key \'Sample_ID\' is reserved. ", 
                      "Choose a different key value.")
  warn_msg1 <- "Coerce `complete_cases = TRUE` at `impute_na = FALSE`."
  
  col_select <- rlang::enexpr(col_select)
  col_group  <- rlang::enexpr(col_group)
  col_order  <- rlang::enexpr(col_order)
  col_color  <- rlang::enexpr(col_color)
  col_fill   <- rlang::enexpr(col_fill)
  col_shape  <- rlang::enexpr(col_shape)
  col_size   <- rlang::enexpr(col_size)
  col_alpha  <- rlang::enexpr(col_alpha)
  col_benchmark <- rlang::enexpr(col_benchmark)
  
  col_select <- if (is.null(col_select)) {
    rlang::expr(Select)
  } else {
    rlang::sym(col_select)
  }

  col_group <- if (is.null(col_group)) {
    rlang::expr(Group)
  } else {
    rlang::sym(col_group)
  }

  col_order <- if (is.null(col_order)) {
    rlang::expr(Order)
  } else {
    rlang::sym(col_order)
  }

  col_color <- if (is.null(col_color)) {
    rlang::expr(Color)
  } else if(suppressWarnings(is.na(col_color))) {
    rlang::sym(".")
  } else {
    rlang::sym(col_color)
  }

  col_fill <- if (is.null(col_fill)) {
    rlang::expr(Fill)
  } else if(suppressWarnings(is.na(col_fill))) {
    rlang::sym(".")
  } else {
    rlang::sym(col_fill)
  }

  col_shape <- if (is.null(col_shape)) {
    rlang::expr(Shape)
  } else if(suppressWarnings(is.na(col_shape))) {
    rlang::sym(".")
  } else {
    rlang::sym(col_shape)
  }

  col_size <- if (is.null(col_size)) {
    rlang::expr(Size)
  } else if(suppressWarnings(is.na(col_size))) {
    rlang::sym(".")
  } else {
    rlang::sym(col_size)
  }

  col_alpha <- if (is.null(col_alpha)) {
    rlang::expr(Alpha)
  } else if(suppressWarnings(is.na(col_alpha))) {
    rlang::sym(".")
  } else {
    rlang::sym(col_alpha)
  }

  col_benchmark <- if (is.null(col_benchmark)) {
    rlang::expr(Benchmark)
  } else {
    rlang::sym(col_benchmark)
  }

  if (col_select == rlang::expr(Sample_ID)) stop(err_msg1)
  if (col_group == rlang::expr(Sample_ID))  stop(err_msg1)
  if (col_order == rlang::expr(Sample_ID))  stop(err_msg1)
  if (col_color == rlang::expr(Sample_ID))  stop(err_msg1)
  if (col_fill == rlang::expr(Sample_ID))   stop(err_msg1)
  if (col_shape == rlang::expr(Sample_ID))  stop(err_msg1)
  if (col_size == rlang::expr(Sample_ID))   stop(err_msg1)
  if (col_alpha == rlang::expr(Sample_ID))  stop(err_msg1)
  if (col_benchmark == rlang::expr(Sample_ID)) stop(err_msg1)
  
  color_brewer <- rlang::enexpr(color_brewer)
  fill_brewer  <- rlang::enexpr(fill_brewer)
  if (!is.null(color_brewer)) color_brewer <- rlang::as_string(color_brewer)
  if (!is.null(fill_brewer))  fill_brewer  <- rlang::as_string(fill_brewer)
  
  size_manual  <- rlang::enexpr(size_manual)
  shape_manual <- rlang::enexpr(shape_manual)
  alpha_manual <- rlang::enexpr(alpha_manual)
  
  df  <- rlang::enexpr(df)
  df2 <- rlang::enexpr(df2)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  dat_dir <- get_gl_dat_dir()
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  
  if (is.null(label_scheme[[col_select]])) {
    stop("Column \'", rlang::as_string(col_select), "\' not found.")
  } else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    stop("No samples under column \'", rlang::as_string(col_select), "\'.")
  }

  if (is.null(label_scheme[[col_group]])) {
    warning("Column \"", rlang::as_string(col_group), 
            "\" not found; use column \"Select\".")
    col_group <- rlang::expr(Select)
  } else if (sum(!is.na(label_scheme[[col_group]])) == 0) {
    warning("No samples under \"", rlang::as_string(col_group), 
            "\"; use column \"Select\".")
    col_group <- rlang::expr(Select)
  }
  
  if (is.null(label_scheme[[col_order]])) {
    warning("Column \"", rlang::as_string(col_order), 
            "\" not found; arranged by the alphebatics.")
  } else if (sum(!is.na(label_scheme[[col_order]])) == 0) {
    # warning("No samples under column \"", rlang::as_string(col_order), "\".")
  }
  
  if (is.null(label_scheme[[col_color]]) && rlang::as_string(col_color) != ".") {
    warning("Column \'", rlang::as_string(col_color), 
            "\' not found.", 
            call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_color]])) == 0) {
    # warning("No samples under column \'", rlang::as_string(col_color), "\'.")
  }
  
  if (is.null(label_scheme[[col_fill]]) && rlang::as_string(col_fill) != ".") {
    warning("Column \"", rlang::as_string(col_fill), "\" not found.")
  } else if(sum(!is.na(label_scheme[[col_fill]])) == 0) {
    # warning("No samples under column \"", rlang::as_string(col_fill), "\".")
  }
  
  if (is.null(label_scheme[[col_shape]]) && rlang::as_string(col_shape) != ".") {
    warning("Column \"", rlang::as_string(col_shape), "\" not found.")
  } else if(sum(!is.na(label_scheme[[col_shape]])) == 0) {
    # warning("No samples under column \"", rlang::as_string(col_shape), "\".")
  }
  
  if (is.null(label_scheme[[col_size]]) && rlang::as_string(col_size) != ".") {
    warning("Column \"", rlang::as_string(col_size), "\" not found.")
  } else if (sum(!is.na(label_scheme[[col_size]])) == 0) {
    # warning("No samples under column \"", rlang::as_string(col_size), "\".")
  }
  
  if(is.null(label_scheme[[col_alpha]]) && rlang::as_string(col_alpha) != ".") {
    warning("Column \"", rlang::as_string(col_alpha), "\" not found.")
  } else if(sum(!is.na(label_scheme[[col_alpha]])) == 0) {
    # warning("No samples under column \"", rlang::as_string(col_alpha), "\".")
  }
  
  if (is.null(label_scheme[[col_benchmark]])) {
    warning("Column \"", rlang::as_string(col_benchmark), "\" not found.")
  } else if (sum(!is.na(label_scheme[[col_benchmark]])) == 0) {
    # warning("No samples under column \"", rlang::as_string(col_benchmark), "\".")
  }
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  if (length(id) != 1L) {
    stop("'id' must be one of 'pep_seq', 'pep_seq_mod', 'prot_acc' or 'gene'.")
  }

  if (id %in% c("prot_acc", "gene")) {
    data_type <- "Protein"
  } else if (id %in% c("pep_seq", "pep_seq_mod")) {
    data_type <- "Peptide"
  } else {
    stop("Unrecognized 'id'; ", 
         "needs to be in c(\"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\")")
  }

  anal_type <- rlang::as_string(rlang::enexpr(anal_type))
  
  if (is.null(filepath)) {
    filepath <- if (grepl("Trend", anal_type)) {
      file.path(dat_dir, data_type, "Trend")
    } else if (grepl("NMF", anal_type)) {
      file.path(dat_dir, data_type, "NMF")
    } else if (grepl("GSPA", anal_type)) {
      file.path(dat_dir, data_type, "GSPA")
    } else {
      file.path(dat_dir, data_type, anal_type)
    }

    dir.create(
      file.path(filepath, "log"), recursive = TRUE, showWarnings = FALSE)
  } else {
    stop("Use default `filepath`.")
  }
  
  if (is.null(filename)) {
    fn_prefix <- paste(data_type, anal_type, sep = "_")
    
    fn_prefix <- local({
      s_type <- if (is.na(scale_log2r)) "O" else if (scale_log2r) "Z" else "N"
      paste0(fn_prefix, "_", s_type)
    })
    
    ## ifelse handles scale_log2r = NA
    # fn_prefix <- paste0(fn_prefix, "_", ifelse(scale_log2r, "Z", "N"))
    fn_prefix <- fn_prefix %>% ifelse(impute_na, paste0(., "_impNA"), .)
    fn_suffix <- if (anal_type %in% c("Model", "KinSub")) "txt" else "png"
  } else {
    if (length(filename) > 1L) {
      stop("Do not provide multiple file names.")
    }
    
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    
    if (fn_prefix == fn_suffix) {
      stop("No '.' in the file name.")
    }

    if (anal_type %in% c("Trend", "NMF", "GSPA")) {
      fn_prefix <- paste(fn_prefix, data_type, anal_type, sep = "_")
      
      fn_prefix <- local({
        s_type <- if (is.na(scale_log2r)) "O" else if (scale_log2r) "Z" else "N"
        paste0(fn_prefix, "_", s_type)
      })
      
      # fn_prefix <- paste0(fn_prefix, "_", ifelse(scale_log2r, "Z", "N"))
      fn_prefix <- fn_prefix %>% ifelse(impute_na, paste0(., "_impNA"), .)
    }
  }
  
  use_pri_data <- c("MDS", "PCA", "EucDist", "Heatmap", 
                    "Histogram", "Corrplot",
                    "Model", "Volcano", "Trend", "NMF", "NMF_meta", 
                    "GSPA", "mapGSPA",
                    "GSVA", "GSEA", "String", "LDA", 
                    "KinSub")
  
  use_sec_data <- c("Trend_line", "NMF_con", "NMF_coef", 
                    "NMF_meta", "GSPA_hm", "mapGSPA")
  
  if (anal_type %in% use_pri_data) {
    df <- find_pri_df(anal_type = !!anal_type, id = !!id, impute_na = impute_na)
    
    if (!is.null(dim(df))) {
      df <- rm_pval_whitespace(df)
    }
  } else {
    df <- NULL
  }
  
  if (anal_type %in% use_sec_data) {
    df2 <- rlang::eval_bare(df2, env = current_env())
    vararg_secmsg(id = !!id, anal_type = !!anal_type)
  } else {
    df2 <- NULL
  }
  
  if (anal_type %in% c("Model", "GSPA", "GSVA", "GSEA")) {
    label_scheme_sub <- label_scheme %>% # to be subset by "formulas"
      dplyr::filter(!grepl("^Empty\\.[0-9]+", .$Sample_ID), !Reference)
  } else {
    label_scheme_sub <- label_scheme %>%
      dplyr::filter(!is.na(!!col_select))
  }
  
  if (!nrow(label_scheme_sub)) {
    stop(paste0("Samples or conditions not defined for \"", anal_type, "\""))
  }

  message("\nscale_log2r = ", scale_log2r)
  message("impute_na = ", impute_na)
  message("complete_cases = ", complete_cases)
  
  ## primary functions:
  # `impute_na` determines the `src_path` for `df`
  # `scale_log2r` determines `N_log2R` or `Z_log2R` columns in `df`
  # `complete_cases` subsets data rows in `df` for samples in `label_scheme_sub`
  #   (`Model` and `GSVA`: `complete_cases` further subjects to lm formulae)
  
  ## Secondary analysis
  # `scale_log2r` for matching '_N' or '_Z' in input filenames from 
  #    the corresponding primary function
  #  (special case of `GSPA`: match to the value in `prnSig`)
  # `impute_na` for matching '[_NZ]_impNA' or '[_NZ]' in input filenames
  # `complete_cases` subsets data rows in primary outputs
  
  if (anal_type == "MDS") {
    function(choice = "cmdscale", dist_co = log2(1), 
             adjEucDist = FALSE, 
             method = "euclidean", 
             p = 2, k = 3, dimension = 2, folds = 1,
             show_ids = TRUE, show_ellipses = FALsE, 
             center_features = TRUE, scale_features = TRUE,
             theme = NULL, ...) {
      plotMDS(df = df,
              id = !!id,
              label_scheme_sub = label_scheme_sub,
              choice = choice,
              dist_co = dist_co, 
              adjEucDist = adjEucDist,
              method = method,
              p = p,
              k = k,
              dimension = dimension,
              folds = folds,
              show_ids = show_ids,
              show_ellipses = show_ellipses,
              col_group = !!col_group,
              col_color = !!col_color,
              col_fill = !!col_fill,
              col_shape = !!col_shape,
              col_size = !!col_size,
              col_alpha = !!col_alpha,
              color_brewer = !!color_brewer,
              fill_brewer = !!fill_brewer,
              size_manual = size_manual,
              shape_manual = shape_manual,
              alpha_manual = alpha_manual,
              scale_log2r = scale_log2r,
              complete_cases = complete_cases,
              filepath = filepath,
              filename = paste0(fn_prefix, ".", fn_suffix),
              center_features = center_features,
              scale_features = scale_features,
              theme = theme,
              anal_type = anal_type,
              ...)
    }
  } 
  else if (anal_type == "PCA") {
    function(choice = "prcomp", type = "obs", dimension = 2, folds = 1,
             show_ids = TRUE, show_ellipses = FALsE, 
             center_features = TRUE, scale_features = TRUE,
             theme = NULL, ...) {
      plotPCA(df = df,
              id = !!id,
              label_scheme_sub = label_scheme_sub,
              choice = choice,
              type = type,
              dimension = dimension,
              folds = folds,
              show_ids = show_ids,
              show_ellipses = show_ellipses,
              col_group = !!col_group,
              col_color = !!col_color,
              col_fill = !!col_fill,
              col_shape = !!col_shape,
              col_size = !!col_size,
              col_alpha = !!col_alpha,
              color_brewer = !!color_brewer,
              fill_brewer = !!fill_brewer,
              size_manual = size_manual,
              shape_manual = shape_manual,
              alpha_manual = alpha_manual,
              scale_log2r = scale_log2r,
              complete_cases = complete_cases,
              impute_na = impute_na,
              filepath = filepath,
              filename = paste0(fn_prefix, ".", fn_suffix),
              center_features = center_features,
              scale_features = scale_features,
              theme = theme,
              anal_type = anal_type,
              ...)
    }
  } 
  else if (anal_type == "EucDist") {
    function(adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, ...) {
      plotEucDist(df = df,
                  id = !!id,
                  label_scheme_sub = label_scheme_sub,
                  adjEucDist = adjEucDist,
                  scale_log2r = scale_log2r,
                  complete_cases = complete_cases,
                  annot_cols = annot_cols,
                  annot_colnames = annot_colnames,
                  filepath = filepath,
                  filename = paste0(fn_prefix, ".", fn_suffix),
                  anal_type = anal_type,
                  ...)
    }
  } 
  else if (anal_type == "Heatmap") {
    function(xmin = -1, xmax = 1, xmargin = 0.1,
             annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL,
             p_dist_rows = 2, p_dist_cols = 2,
             hc_method_rows = "complete", hc_method_cols = "complete", 
             
             x = NULL, 
             p = NULL, 
             method = NULL, 
             diag = NULL, 
             upper = NULL, 
             annotation_col = NULL, 
             annotation_row = NULL, 
             clustering_method = NULL, 
             rm_allna = TRUE, 
             ...) {
      plotHM(df = df,
             id = !!id,
             col_select = !!col_select, 
             col_order  = !!col_order,
             col_benchmark = !!col_benchmark,
             label_scheme_sub = label_scheme_sub,
             filepath = filepath,
             filename = paste0(fn_prefix, ".", fn_suffix),
             scale_log2r = scale_log2r,
             complete_cases = complete_cases,
             annot_cols = annot_cols,
             annot_colnames = annot_colnames,
             annot_rows = annot_rows,
             xmin = xmin,
             xmax = xmax,
             xmargin = xmargin,
             p_dist_rows = p_dist_rows,
             p_dist_cols = p_dist_cols,
             hc_method_rows = hc_method_rows,
             hc_method_cols = hc_method_cols,
             
             x = x, 
             p = p, 
             method = method, 
             diag = diag, 
             upper = upper, 
             annotation_col = annotation_col, 
             annotation_row = annotation_row, 
             clustering_method = clustering_method, 
             rm_allna = rm_allna, 
             ...)
    }
  } 
  else if (anal_type == "Histogram") {
    function(cut_points = NA, show_curves = TRUE, show_vline = TRUE, 
             scale_y = TRUE, theme = NULL, ...) {
      plotHisto(df = df,
                id = !!id,
                label_scheme_sub = label_scheme_sub,
                scale_log2r = scale_log2r,
                complete_cases = complete_cases,
                cut_points = cut_points,
                show_curves = show_curves,
                show_vline = show_vline,
                scale_y = scale_y,
                filepath = filepath,
                filename = paste0(fn_prefix, ".", fn_suffix),
                theme = theme,
                ...)
    }
  } 
  else if (anal_type == "Corrplot") {
    function(data_select = "logFC", 
             cor_method = "pearson", digits = 2L, ...) {
      plotCorr(df = df,
               id = !!id,
               anal_type = anal_type,
               data_select = data_select,
               col_select = !!col_select,
               col_order = !!col_order,
               label_scheme_sub = label_scheme_sub,
               scale_log2r = scale_log2r,
               complete_cases = complete_cases,
               filepath = filepath,
               filename = paste0(fn_prefix, "_", data_select, ".", fn_suffix),
               cor_method = cor_method, 
               digits = digits,
               ...)
    }
  } 
  else if (anal_type == "Model") {
    function(method = "limma", padj_method = "BH", 
             var_cutoff = 1E-3, pval_cutoff = 1, logFC_cutoff = log2(1), 
             rm_allna = FALSE, ...) {
      sigTest(df = df,
              id = !!id,
              label_scheme_sub = label_scheme_sub,
              scale_log2r = scale_log2r,
              complete_cases = complete_cases,
              impute_na = impute_na,
              impute_group_na = impute_group_na, 
              rm_allna = rm_allna, 
              method_replace_na = method_replace_na, 
              filepath = filepath,
              filename = paste0(fn_prefix, ".", fn_suffix),
              method = !!method,
              padj_method = padj_method, 
              var_cutoff = var_cutoff,
              pval_cutoff = pval_cutoff,
              logFC_cutoff = logFC_cutoff,
              data_type = data_type,
              anal_type = anal_type,
              ...)
    }
  } 
  else if (anal_type == "Volcano") {
    function(fml_nms = NULL, adjP = FALSE, topn_labels = 20, 
             theme = NULL, highlights = NULL, grids = NULL, ...) {
      plotVolcano(df = df,
                  df2 = NULL,
                  id = !!id,
                  id_gspa = id_gspa, 
                  adjP = adjP,
                  topn_labels = topn_labels, 
                  anal_type = anal_type,
                  gspval_cutoff = 1,
                  gslogFC_cutoff = 0,
                  topn_gsets = 0,
                  show_sig = "none",
                  fml_nms = fml_nms,
                  gset_nms = NULL,
                  scale_log2r = scale_log2r,
                  complete_cases = complete_cases,
                  impute_na = impute_na,
                  filepath = filepath,
                  filename = paste0(fn_prefix, ".", fn_suffix),
                  highlights = highlights, 
                  grids = grids, 
                  theme = theme,
                  ...)
    }
  } 
  else if (anal_type == "mapGSPA") {
    function(fml_nms = NULL, adjP = FALSE, topn_labels = 20, 
             gspval_cutoff = 0.05, gslogFC_cutoff = log2(1.2), 
             topn_gsets = Inf, show_sig = "none", 
             gset_nms = "go_sets",
             theme = NULL, ...) {
      plotVolcano(df = df,
                  df2 = df2,
                  id = !!id,
                  id_gspa = id_gspa, 
                  adjP = adjP,
                  topn_labels = topn_labels, 
                  anal_type = anal_type,
                  gspval_cutoff = gspval_cutoff,
                  gslogFC_cutoff = gslogFC_cutoff,
                  topn_gsets = topn_gsets,
                  show_sig = show_sig,
                  fml_nms = fml_nms,
                  gset_nms = gset_nms,
                  scale_log2r = scale_log2r,
                  complete_cases = complete_cases,
                  impute_na = impute_na,
                  filepath = filepath,
                  filename = paste0(fn_prefix, ".", fn_suffix),
                  theme = theme,
                  ...)
    }
  } 
  else if (anal_type == "Trend") {
    function(choice = "cmeans", n_clust = NULL, ...) {
      analTrend(df = df,
                id = !!id,
                col_group = !!col_group,
                col_order = !!col_order,
                label_scheme_sub = label_scheme_sub,
                choice = choice,
                n_clust = n_clust,
                scale_log2r = scale_log2r,
                complete_cases = complete_cases,
                impute_na = impute_na,
                filepath = filepath,
                filename = paste0(fn_prefix %>% paste0(., "_nclust", n_clust), ".txt"),
                anal_type = anal_type,
                ...)
    }
  } 
  else if (anal_type == "Trend_line") {
    function(n_clust = NULL, theme = NULL, ...) {
      plotTrend(df2 = df2,
                id = !!id,
                col_group = !!col_group,
                col_order = !!col_order,
                label_scheme_sub = label_scheme_sub,
                n_clust = n_clust,
                scale_log2r = scale_log2r,
                complete_cases = complete_cases,
                impute_na = impute_na,
                filepath = filepath,
                filename = paste0(fn_prefix, ".", fn_suffix),
                theme = theme,
                ...)
    }
  } 
  else if (anal_type == "NMF") {
    function(rank = NULL, nrun = 50, seed = NULL, ...) {
      analNMF(df = df,
              id = !!id,
              rank = rank,
              nrun = nrun,
              seed = seed,
              col_group = !!col_group,
              label_scheme_sub = label_scheme_sub,
              scale_log2r = scale_log2r,
              complete_cases = complete_cases,
              impute_na = impute_na,
              filepath = filepath,
              filename = paste0(fn_prefix %>% paste0(., "_rank", rank), ".txt"),
              anal_type = anal_type,
              ...)
    }
  } 
  else if (anal_type == "NMF_con") {
    function(rank = NULL, ...) {
      plotNMFCon(df2 = df2,
                 id = !!id,
                 rank = rank,
                 label_scheme_sub = label_scheme_sub,
                 scale_log2r = scale_log2r,
                 complete_cases = complete_cases,
                 impute_na = impute_na,
                 filepath = filepath,
                 filename = paste0(fn_prefix, ".", fn_suffix),
                 ...)
    }
  } 
  else if (anal_type == "NMF_coef") {
    function(rank = NULL, ...) {
      plotNMFCoef(df2 = df2,
                  id = !!id,
                  rank = rank,
                  label_scheme_sub = label_scheme_sub,
                  scale_log2r = scale_log2r,
                  complete_cases = complete_cases,
                  impute_na = impute_na,
                  filepath = filepath,
                  filename = paste0(fn_prefix, ".", fn_suffix),
                  ...)
    }
  } 
  else if (anal_type == "NMF_meta") {
    function(rank = NULL, ...) {
      plotNMFmeta(df = df,
                  df2 = df2,
                  id = !!id,
                  rank = rank,
                  label_scheme_sub = label_scheme_sub,
                  scale_log2r = scale_log2r,
                  complete_cases = complete_cases,
                  impute_na = impute_na,
                  filepath = filepath,
                  filename = paste0(fn_prefix, ".", fn_suffix),
                  anal_type = anal_type,
                  ...)
    }
  } 
  else if (anal_type == "GSPA") {
    function(gset_nms = "go_sets", var_cutoff = .5,
             pval_cutoff = 1E-2, logFC_cutoff = log2(1.1),
             gspval_cutoff = 1E-2, gslogFC_cutoff = log2(1),
             min_size = 10, max_size = Inf, min_delta = 4, min_greedy_size = 1,
             use_adjP = FALSE,
             method = "mean",
             ...) {
      gspaTest(df = df,
               id = !!id,
               id_gspa = id_gspa, 
               label_scheme_sub = label_scheme_sub,
               scale_log2r = scale_log2r,
               complete_cases = complete_cases,
               impute_na = impute_na,
               filepath = filepath,
               filename = paste0(fn_prefix, ".txt"),
               gset_nms = gset_nms,
               var_cutoff = var_cutoff,
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
               anal_type = anal_type,
               ...)
    }
  } 
  else if (anal_type == "GSPA_hm") {
    function(...) {
      gspaHM(df2 = df2,
             scale_log2r = scale_log2r,
             complete_cases = complete_cases,
             impute_na = impute_na,
             filepath = filepath,
             filename = paste0(fn_prefix, ".", fn_suffix),
             ...)
    }
  } 
  else if(anal_type == "GSVA") {
    function(lm_method = "limma", padj_method = "BH", 
             gset_nms = "go_sets", var_cutoff = .5,
             pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), ...) {
      gsvaTest(df = df,
               id = !!id,
               label_scheme_sub = label_scheme_sub,
               filepath = filepath,
               filename = paste0(fn_prefix, ".txt"),
               scale_log2r = scale_log2r,
               complete_cases = complete_cases,
               impute_na = impute_na,
               gset_nms = gset_nms,
               lm_method = lm_method, 
               padj_method = padj_method,
               var_cutoff = var_cutoff,
               pval_cutoff = pval_cutoff,
               logFC_cutoff = logFC_cutoff,
               anal_type = anal_type,
               ...)
    }
  } 
  else if (anal_type == "GSEA") {
    function(gset_nms = "go_sets", var_cutoff = 0.5, pval_cutoff = 1E-2, 
             logFC_cutoff = log2(1.1), ...) {
      gspaTest(df = df,
               id = !!id,
               label_scheme_sub = label_scheme_sub,
               scale_log2r = scale_log2r,
               complete_cases = complete_cases,
               impute_na = impute_na,
               filepath = filepath,
               filename = paste0(fn_prefix, ".txt"),
               gset_nms = gset_nms,
               var_cutoff = var_cutoff,
               pval_cutoff = pval_cutoff,
               logFC_cutoff = logFC_cutoff,
               gspval_cutoff = 0.05, # dummy
               gslogFC_cutoff = log2(1),
               min_size = 10,
               max_size = Inf,
               min_delta = 1,
               min_greedy_size = 1,
               use_adjP = FALSE,
               method = method, # dummy end
               anal_type = anal_type,
               ...)
    }
  } 
  else if (anal_type == "String") {
    function(db_nms = NULL, score_cutoff = .7, ...) {
      stringTest(df = df,
                 id = !!id,
                 label_scheme_sub = label_scheme_sub,
                 db_nms = db_nms,
                 score_cutoff = score_cutoff,
                 scale_log2r = scale_log2r,
                 complete_cases = complete_cases,
                 filepath = filepath,
                 filename = paste0(fn_prefix, ".csv"),
                 ...)
    }
  } 
  else if (anal_type == "LDA") {
    function(choice = "lda", method = "moment", 
             type = "obs", dimension = 2, folds = 1, show_ids = TRUE,
             show_ellipses = FALsE, 
             center_features = TRUE, scale_features = TRUE,
             theme = NULL, 
             
             formula = NULL, 
             data = NULL, 
             x = NULL, 
             grouping = NULL, 
             prior = NULL, 
             subset = NULL, 
             CV = NULL, 
             na.action = NULL, 
             nu = NULL, 
             ...) {
      plotLDA(df = df,
              id = !!id,
              label_scheme_sub = label_scheme_sub,
              choice = choice, 
              method = method,
              type = type,
              dimension = dimension,
              folds = folds,
              show_ids = show_ids,
              show_ellipses = show_ellipses,
              col_group = !!col_group,
              col_color = !!col_color,
              col_fill = !!col_fill,
              col_shape = !!col_shape,
              col_size = !!col_size,
              col_alpha = !!col_alpha,
              color_brewer = !!color_brewer,
              fill_brewer = !!fill_brewer,
              size_manual = size_manual,
              shape_manual = shape_manual,
              alpha_manual = alpha_manual,
              scale_log2r = scale_log2r,
              complete_cases = complete_cases,
              impute_na = impute_na,
              filepath = filepath,
              filename = paste0(fn_prefix, ".", fn_suffix),
              center_features = center_features,
              scale_features = scale_features,
              theme = theme,
              anal_type = anal_type,
              
              formula = formula, 
              data = data, 
              x = x, 
              grouping = grouping, 
              prior = prior, 
              subset = subset, 
              CV = CV, 
              na.action = na.action, 
              nu = nu, 
              ...)
    }
  } 
  else if (anal_type == "KinSub") {
    function(db_nms = NULL, match_orgs = TRUE, ...) {
      KinSubTest(df = df,
                 id = !!id,
                 label_scheme_sub = label_scheme_sub,
                 db_nms = db_nms,
                 match_orgs = match_orgs, 
                 scale_log2r = scale_log2r,
                 complete_cases = complete_cases,
                 filepath = filepath,
                 filename = paste(fn_prefix, fn_suffix, sep = "."),
                 ...)
    }
  }
}



#' helper for finding input `df`
#'
#' @param ... Not currently used.
#' @inheritParams info_anal
find_pri_df <- function (anal_type = "Model", df = NULL, 
                         id = "gene", impute_na = FALSE, ...) 
{
  dat_dir <- get_gl_dat_dir()
  
  err_msg2 <- paste0("not found. \n", 
                     "Run functions of PSM, peptide and protein normalization first.")
  err_msg3 <- paste0("not found. \n", 
                     "Impute NA values with ", 
                     "`pepImp()` or `prnImp()` or set `impute_na = FALSE`.")
  err_msg4 <- "not found at impute_na = TRUE. \nRun `prnSig(impute_na = TRUE)` first."
  err_msg5 <- "not found at impute_na = FALSE. \nRun `prnSig(impute_na = FALSE)` first."

  anal_type <- rlang::as_string(rlang::enexpr(anal_type))
  id <- rlang::as_string(rlang::enexpr(id))

  df <- rlang::enexpr(df)

  if (is.null(df)) {
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      fn_p <- file.path(dat_dir, "Peptide/Model", "Peptide_pVals.txt")
      fn_imp_p <- file.path(dat_dir, "Peptide/Model", "Peptide_impNA_pVals.txt")
      fn_raw <- file.path(dat_dir, "Peptide", "Peptide.txt")
      fn_imp <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
    } else if (id %in% c("prot_acc", "gene")) {
      fn_p <- file.path(dat_dir, "Protein/Model", "Protein_pVals.txt")
      fn_imp_p <- file.path(dat_dir, "Protein/Model", "Protein_impNA_pVals.txt")
      fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")
      fn_imp <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
    } else {
      stop("Unknown `id`.", call. = FALSE)
    }

    if (anal_type %in% c("Histogram")) { # never impute_na and no pVals
      if (file.exists(fn_raw)) 
        src_path <- fn_raw
      else 
        stop(paste(fn_raw, err_msg2), call. = FALSE)

      if (id %in% c("pep_seq", "pep_seq_mod")) 
        message("Primary column keys in `Peptide/Peptide.txt` ", 
                "for `filter_` varargs.")
      else if (id %in% c("prot_acc", "gene")) 
        message("Primary column keys in `Protein/Protein.txt` ", 
                "for `filter_` varargs.")
    } else if (anal_type %in% c("Model")) { # optional impute_na but no pVals
      if (impute_na) {
        if (file.exists(fn_imp)) 
          src_path <- fn_imp
        else 
          stop(paste(fn_imp, err_msg3), call. = FALSE)
      } else {
        if (file.exists(fn_raw)) 
          src_path <- fn_raw
        else 
          stop(paste(fn_raw, err_msg2), call. = FALSE)
      }

      if (id %in% c("pep_seq", "pep_seq_mod")) 
        message("Primary column keys in `Peptide/Peptide[_impNA].txt` ", 
                "for `filter_` varargs.")
      else if (id %in% c("prot_acc", "gene")) 
        message("Primary column keys in `Protein/Protein[_impNA].txt` ", 
                "for `filter_` varargs.")
    } else if (anal_type %in% c("Volcano", "GSPA", "mapGSPA", "GSEA")) { # data with pVals
      if (impute_na) {
        if (file.exists(fn_imp_p)) 
          src_path <- fn_imp_p
        else 
          stop(paste(fn_imp_p, err_msg4), call. = FALSE)
      } else {
        if (file.exists(fn_p)) 
          src_path <- fn_p
        else 
          stop(paste(fn_p, err_msg5), call. = FALSE)
      }

      if (id %in% c("pep_seq", "pep_seq_mod")) 
        message("Primary column keys in `Model/Peptide[_impNA]_pVals.txt` ", 
                "for `filter_` varargs.")
      else if (id %in% c("prot_acc", "gene")) 
        message("Primary column keys in `Model/Protein[_impNA]_pVals.txt` ", 
                "for `filter_` varargs.")
    } else if (anal_type %in% c("Heatmap", "MDS", "PCA", "EucDist", "Trend", 
                                "NMF", "NMF_meta",
                                "GSVA", "Corrplot", "String", "LDA", 
                                "KinSub")) { # optional impute_na and possible pVals
      if (impute_na) {
        if (file.exists(fn_imp_p)) 
          src_path <- fn_imp_p
        else if (file.exists(fn_imp)) 
          src_path <- fn_imp
        else 
          stop(paste(fn_imp, err_msg3), call. = FALSE)
      } else {
        if (file.exists(fn_p)) 
          src_path <- fn_p
        else if (file.exists(fn_raw)) 
          src_path <- fn_raw
        else 
          stop(paste(fn_raw, err_msg2), call. = FALSE)
      }

      if (id %in% c("pep_seq", "pep_seq_mod")) 
        message("Primary column keys in `Model/Peptide[_impNA_pVals].txt` ", 
                "for `filter_` varargs.")
      else if (id %in% c("prot_acc", "gene")) 
        message("Primary column keys in `Model/Protein[_impNA_pVals].txt` ", 
                "for `filter_` varargs.")
    }
  } else {
    df <- rlang::as_string(df)

    if (anal_type == "Model") 
      stop("Use default file name.", call. = FALSE)

    if (id %in% c("pep_seq", "pep_seq_mod")) {
      src_path <- file.path(dat_dir, "Peptide/Model", df)
      
      if (!file.exists(src_path)) 
        src_path <- file.path(dat_dir, "Peptide", df)
    } else if (id %in% c("prot_acc", "gene")) {
      src_path <- file.path(dat_dir, "Protein/Model", df)
      
      if (!file.exists(src_path)) 
        src_path <- file.path(dat_dir, "Protein", df)
    }
  }

  df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
                          comment.char = "#"), error = function(e) NA)

  if (is.null(dim(df))) 
    stop(src_path, " not found.", call. = FALSE)

  df <- df %>% reorderCols2()

  message(paste("Primary file loaded:", src_path))

  invisible(df)
}


#' Helper for finding input \code{df} (not currently used).
#'
#' @param ... Not currently used.
#' @inheritParams info_anal
find_sec_df <- function (df = NULL, anal_type = NULL, id = NULL, ...) 
{
  dat_dir <- get_gl_dat_dir()
  
  df <- rlang::enexpr(df)
  anal_type <- rlang::enexpr(anal_type)
  id <- rlang::enexpr(id)

  if (is.null(df) || is.null(anal_type) || is.null(id)) 
    return (NULL)

  df <- rlang::as_string(df)
  anal_type <- rlang::as_string(anal_type)
  id <- rlang::as_string(id)

  new_anal_type <- anal_type %>%
    gsub("^Trend_.*", "Trend", .) %>%
    gsub("^NMF_.*", "NMF", .) %>%
    gsub("^GSPA_.*", "GSPA", .)

  if (id %in% c("pep_seq", "pep_seq_mod")) {
    src_path <- file.path(dat_dir, "Peptide", new_anal_type, df)
  } else if (id %in% c("prot_acc", "gene")) {
    src_path <- file.path(dat_dir, "Protein", new_anal_type, df)
  } else {
    stop("Unknown `id`", call. = FALSE)
  }

  df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
                          comment.char = "#"), error = function(e) NA)

  if (is.null(dim(df))) stop(src_path, " not found.", call. = FALSE)

  invisible(df)
}


#' Helper for finding input \code{df}.
#'
#' @param ... Not currently used.
#' @inheritParams info_anal
vararg_secmsg <- function (id = NULL, anal_type = NULL, ...) 
{
  id <- rlang::as_string(rlang::enexpr(id))
  anal_type <- rlang::as_string(rlang::enexpr(anal_type))

  if (id %in% c("pep_seq", "pep_seq_mod")) {
    if (anal_type == "Trend_line") {
      message("Secondary column keys in `Trend/[...]Peptide_Trend_{NZ}[_impNA][...].txt` ", 
              "for `filter2_` varargs.")
    } else if (anal_type == "NMF_con") {
      message("Secondary column keys in `NMF/[...]Peptide_NMF[...]_consensus.txt` ", 
              "for `filter2_` varargs.")
    } else if (anal_type == "NMF_coef") {
      message("Secondary column keys in `NMF/[...]Peptide_NMF[...]_coef.txt` ", 
              "for `filter2_` varargs.")
    }
  } else if (id %in% c("prot_acc", "gene")) {
    if (anal_type == "mapGSPA") {
      message("Secondary column keys in `GSPA/[...]Protein_GSPA_{NZ}[_impNA].txt` ", 
              "for `filter2_` varargs.")
    } else if (anal_type == "GSPA_hm") {
      message("Secondary column keys in `GSPA/[...]Protein_GSPA_{NZ}_essmap.txt` ", 
              "for `filter2_` varargs.")
      message("Column keys in `GSPA/[...]Protein_GSPA_{NZ}_essmeta.txt` ", 
              "for heat map annotation.")
    } else if (anal_type == "Trend_line") {
      message("Secondary column keys in `Trend/[...]Protein_Trend_{NZ}[_impNA...].txt` ", 
              "for `filter2_` varargs.")
    } else if (anal_type == "NMF_con") {
      message("Secondary column keys in `NMF/[...]Protein_NMF[...]_consensus.txt` ", 
              "for `filter2_` varargs.")
    } else if (anal_type == "NMF_coef") {
      message("Secondary column keys in `NMF/[...]Protein_NMF[...]_coef.txt` ", 
              "for `filter2_` varargs.")
    }
  } else {
    stop("Unknown `id`", call. = FALSE)
  }
}

