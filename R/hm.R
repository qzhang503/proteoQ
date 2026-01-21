#' A wrapper of pheatmap
#'
#' @param mat The same as in \link[pheatmap]{pheatmap}.
#' @param annotation_col The same as in \link[pheatmap]{pheatmap}.
#' @param annotation_row The same as in \link[pheatmap]{pheatmap}.
#' @param color The same as in \link[pheatmap]{pheatmap}.
#' @param annotation_colors The same as in \link[pheatmap]{pheatmap}.
#' @param breaks The same as in \link[pheatmap]{pheatmap}.
#' @param filename The output filename.
#' @param ... Additional arguments for \link[pheatmap]{pheatmap}.
#' 
#' @import dplyr pheatmap
#' @importFrom magrittr %>% %T>% %$% %<>% 
my_pheatmap <- function(mat, filename, annotation_col, annotation_row, 
                        color, annotation_colors, breaks, ...) 
{
  mat <- rlang::enexpr(mat)
  filename <- rlang::enexpr(filename)
  annotation_col <- rlang::enexpr(annotation_col)
  annotation_row <- rlang::enexpr(annotation_row)
  color <- rlang::enexpr(color)
  annotation_colors <- rlang::enexpr(annotation_colors)
  breaks <- rlang::enexpr(breaks)
  
  dots <- rlang::enexprs(...) %>% 
    .[! names(.) %in% c("mat", "filename", "annotation_col", "annotation_row", 
                        "color", "annotation_colors", "breaks")]

  ph_call <- rlang::expr(
    pheatmap(mat = !!mat, 
             filename = !!filename, 
             annotation_col = !!annotation_col, 
             annotation_row = !!annotation_row, 
             color = !!color,
             annotation_colors = !!annotation_colors, 
             breaks = !!breaks, 
             !!!dots))
  
  rlang::eval_bare(ph_call, env = caller_env())
}


#' Makes heat maps
#' 
#' @inheritParams prnHM
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import stringr dplyr ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotHM <- function(df, id, col_select, col_order, col_benchmark, label_scheme_sub, 
                   filepath, filename, scale_log2r, complete_cases, 
                   annot_cols = NULL, annot_colnames = NULL, 
                   annot_rows = annot_rows, 
                   xmin = -1, xmax = 1, xmargin = .1, 
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
                   ...) 
{
  # (1) `x`, `p` etc. defined as NULL @param
  dummies <- c("x", "diag", "upper", "method", "p", 
               "annotation_col", "annotation_row", "clustering_method")
  msgs <- c(
    "`x` in `dist()` automated.",
    "`diag` in `dist()` automated.",
    "`upper` in `dist()` automated.",
    paste0("`method` in `dist()` replaced with ", 
           "`clustering_distance_rows` and `clustering_distance_cols` ", 
           "in `pheatmap()`."), 
    paste0("`p` in `dist()` replaced with ", 
           "`p_dist_rows` and `p_dist_cols`."), 
    paste0("`annotation_col` in `pheatmap()` disabled; \n", 
           "instead, use `annot_cols`."), 
    paste0("`annotation_row` in `pheatmap()` disabled; \n", 
           "instead, use `annot_rows`."), 
    paste0("`clustering_method` in `pheatmap()` split into ", 
           "`hc_method_rows` and `hc_method_cols`.")
  )
  
  stopifnot(length(dummies) == length(msgs))
  
  # (2) checking (for developer only)
  check_formalArgs(prnHM, dist, dummies)
  check_formalArgs(pepHM, dist, dummies)
  
  # (3) values back to default
  purrr::walk2(dummies, msgs, ~ {
    if (!is.null(get(.x, envir = rlang::env_parent(), inherits = FALSE))) {
      warning(.y)
      assign(.x, NULL, envir = rlang::env_parent(), inherits = FALSE)
    } 
  })
  
  stopifnot(vapply(c(xmin, xmax, xmargin, p_dist_rows, p_dist_cols), 
                   is.numeric, logical(1L)))
  stopifnot(vapply(c(hc_method_rows, hc_method_cols), 
                   rlang::is_string, logical(1L)))
  stopifnot(xmin < xmax, xmargin >= 0, xmargin <= abs(xmax))
  stopifnot(p_dist_rows > 0, p_dist_cols > 0)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  dir.create(file.path(filepath, "Subtrees", fn_prefix), 
             recursive = TRUE, showWarnings = FALSE)
  
  id <- rlang::as_string(rlang::enexpr(id))
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  col_benchmark <- rlang::as_string(rlang::enexpr(col_benchmark))
  
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  arrange_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange_", names(.))]
  
  select_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^select_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter_dots, arrange_dots, select_dots)]

  # needed defaults before calling pheatmap
  cluster_rows <- if (is.null(dots$cluster_rows)) TRUE else dots$cluster_rows
  cluster_cols <- if (is.null(dots$cluster_cols)) TRUE else dots$cluster_cols

  clustering_distance_rows <- if (is.null(dots$clustering_distance_rows)) {
    "euclidean"
  }
  else {
    dots$clustering_distance_rows
  }

  clustering_distance_cols <- if (is.null(dots$clustering_distance_cols)) {
    "euclidean"
  }
  else {
    dots$clustering_distance_cols
  }

  if (!is.null(dots$clustering_method)) {
    dots$clustering_method <- NULL
    
    warning("Argument `clustering_method` disabled; 
            use `hc_method_rows` and `hc_method_cols` instead.")    
  }

  n_color <- 500
  
  if (is.null(dots$breaks)) {
    color_breaks <- c(seq(xmin, -xmargin, length.out = n_color/2)[1:(n_color/2-1)],
                      seq(-xmargin, xmargin, length.out = 3),
                      seq(xmargin, xmax, length.out = n_color/2)[2:(n_color/2)])
    
    color_breaks <- unique(color_breaks[color_breaks >= xmin])
  } 
  else if (is.na(dots$breaks)) {
    color_breaks <- NA
  } 
  else {
    color_breaks <- eval(dots$breaks, envir = rlang::caller_env())
  }
  
  mypalette <- if (is.null(dots$color)) {
    grDevices::colorRampPalette(c("#1f78b4", "white", "#e31a1c"))(n_color)
  }
  else if (is.na(dots$color)) {
    grDevices::colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  }
  else {
    eval(dots$color, envir = rlang::caller_env())
  }

  x_label <- expression("Ratio ("*log[2]*")")
  NorZ_ratios <- find_NorZ(scale_log2r)
  NorZ_ratios_to_ctrl <- paste("toCtrl", NorZ_ratios, sep = "_")
  
  dat_dir <- get_gl_dat_dir()
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  sample_ids <- label_scheme_sub$Sample_ID
  
  ###
  # may check "NA" and replace with NA
  ###
  
  pattern <- 
    "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("
  
  df <- df %>%
    dplyr::mutate_at(vars(grep("^pVal|^adjP", names(.))), as.numeric) %>%
    dplyr::mutate(Mean_log10Int = log10(rowMeans(.[, grepl("^I[0-9]{3}", names(.))],
                                                 na.rm = TRUE)))
  
  df <- df %>%
    # dplyr::mutate_at(vars(grep("log2_R[0-9]{3}", names(.))), ~ setHMlims(.x, xmin, xmax)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>% 
    reorderCols(endColIndex = grep(pattern, names(.)), col_to_rn = id) 
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  if (!nrow(df)) {
    stop("Zero data rows available after data filtration.")
  }

  dfR <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sample_ids)) %>%
    dplyr::select(as.character(sample_ids)) # ensure the same order
  
  # Add pep_start and pep_end
  add_pep_range <- dots[["add_pep_range"]]
  
  if (id %in% c("pep_seq_mod", "pep_seq") && "pep_start" %in% names(df) && 
      "pep_end" %in% names(df) && isTRUE(add_pep_range)) {
    pep_ids <- paste0(rownames(df), " (", df$pep_start, ":", df$pep_end, ")")
    rownames(df) <- df[[id]] <- pep_ids
  } else {
    pep_ids <- df[[id]]
  }

  df <- df %>%
    dplyr::select(-grep("log2_R[0-9]{3}", names(.))) %>%
    dplyr::bind_cols(., dfR) %>% 
    dplyr::filter(!is.na(.[[id]])) %>% 
    `rownames<-`(.[[id]])
  
  rm(list = c("dfR"))

  if (!is.null(dots$annotation_row)) {
    dots$annotation_row <- NULL
    warning("Argument `annotation_row` disabled; use `annot_rows` instead.")
  }
  
  if (!is.null(dots$annotation_col)) {
    dots$annotation_col <- NULL
    warning("Argument `annotation_col` disabled; use `annot_cols` instead.")
  }

  annotation_col <- if (is.null(annot_cols)) {
    NA
  }
  else {
    colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)
  }

  if ((!is.null(annot_colnames)) && (length(annot_colnames) == length(annot_cols))) {
    colnames(annotation_col) <- annot_colnames
  }

  annotation_row <- if (is.null(annot_rows)) {
    NA
  }
  else {
    dplyr::select(df, annot_rows)
  }

  annotation_colors <- if (is.null(dots$annotation_colors)) {
    setHMColor(annotation_col)
  }
  else if (suppressWarnings(is.na(dots$annotation_colors))) {
    NA
  }
  else {
    eval(dots$annotation_colors, envir = rlang::caller_env())
  }

  if (complete_cases) {
    df_hm <- dplyr::filter(df, complete.cases(.[, names(.) %in% sample_ids]))
  } 
  else {
    df_hm <- df
  }

  df_hm <- df_hm %>%
    `rownames<-`(.[[id]])	%>%
    dplyr::select(which(names(.) %in% sample_ids)) %>% 
    { if (rm_allna) .[rowSums(!is.na(.)) > 0L, ] else . } 
  
  if (!nrow(df_hm)) {
    stop("Zero data rows after removing all-NA rows.")
  }

  # sample orders
  if (cluster_cols && (!is.null(col_order))) {
    message("No column-ordering of samples at `cluster_cols = TRUE`.")
  }

  if ((!cluster_cols) && (!is.null(col_order)) && 
      length(unique(label_scheme_sub[[col_order]])) > 1L) {
    
    df_hm <- local({
      plot_orders <- label_scheme_sub %>%
        dplyr::select(Sample_ID, !!col_order) %>%
        dplyr::filter(!is.na(!!col_order)) %>%
        unique() %>%
        dplyr::arrange(!!col_order)
      
      if (nrow(plot_orders) != length(sample_ids)) {
        stop("The number of entries under `", rlang::as_string(col_order), 
             "` is different to the number of selected samples.")
      }
      
      df_hm <- df_hm[, as.character(plot_orders$Sample_ID), drop = FALSE]
    })
  }

  if (cluster_rows) {
    d <- stats::dist(df_hm, method = clustering_distance_rows, p = p_dist_rows)
    d[is.na(d)] <- .5 * max(d, na.rm = TRUE)

    h <- tryCatch(
      hclust(d, hc_method_rows), 
      error = function(e) 1L
    )
    
    if (class(h) != "hclust" && h == 1L) {
      warning("Row clustering cannot be performed.")
      h <- FALSE
    }

    dots$cluster_rows <- h
    rm(list = c("d", "h"))
  } 
  else {
    dots$cluster_rows <- FALSE
    h <- FALSE
  }
  
  if (cluster_cols) {
    d_cols <- stats::dist(t(df_hm), method = clustering_distance_cols, 
                          p = p_dist_cols)
    d_cols[is.na(d_cols)] <- .5 * max(d_cols, na.rm = TRUE)

    h_cols <- tryCatch(
      hclust(d_cols, hc_method_cols), 
      error = function(e) 1
    )
    
    if ((class(h_cols) != "hclust") && (h_cols == 1)) {
      warning("Column clustering cannot be performed.")
      h_cols <- FALSE
    }
    
    dots$cluster_cols <- h_cols
    # rm(list = c("d_cols", "h_cols")) # h_cols also for subtrees
  } 
  else {
    dots$cluster_cols <- FALSE
    h_cols <- FALSE
  }
  
  filename <- gg_imgname(filename)
  
  # forms `annotation_col` and `annotation_row` from `annot_col` and `annot_row` 
  dots <- dots %>% 
    .[! names(.) %in% c("mat", "filename", "annotation_col", "annotation_row", 
                        "clustering_distance_rows", "clustering_distance_cols", 
                        "clustering_method", 
                        "color", "annotation_colors", "breaks")]
  
  # setHMlims after hclust
  df_hm <- setHMlims(df_hm, xmin, xmax)
  
  # references under expt_smry::Sample_ID may be included
  p <- my_pheatmap(
    mat = df_hm,
    filename = file.path(filepath, filename),
    annotation_col = annotation_col,
    annotation_row = annotation_row, 
    color = mypalette,
    annotation_colors = annotation_colors,
    breaks = color_breaks,
    !!!dots
  )
  
  # subtrees
  cutree_rows <- eval(dots$cutree_rows, envir = rlang::caller_env())
  df <- df %>% dplyr::mutate(!!id := as.character(!!rlang::sym(id)))
  
  if ((!is.null(cutree_rows)) && cluster_rows) {
    if (is.numeric(cutree_rows) && (nrow(df) >= cutree_rows)) {
      Cluster <- data.frame(Cluster = cutree(p$tree_row, k = cutree_rows)) %>%
        dplyr::mutate(!!id := rownames(.)) %>%
        dplyr::left_join(df, by = id) %T>% 
        readr::write_tsv(file.path(filepath, "Subtrees", fn_prefix, 
                                   paste0(fn_prefix, " n-", cutree_rows, "_subtrees.txt")))
      
      Cluster <- Cluster %>%
        `rownames<-`(NULL) %>% 
        tibble::column_to_rownames(var = id)
      
      for (cluster_id in unique(Cluster$Cluster)) {
        df_sub <- Cluster[Cluster$Cluster == cluster_id, names(Cluster) %in% sample_ids]
        
        if (complete_cases) {
          df_sub <- df_sub %>%
            tibble::rownames_to_column(id) %>%
            dplyr::filter(complete.cases(.[, names(.) %in% sample_ids])) %>%
            `rownames<-`(NULL) %>% 
            tibble::column_to_rownames(id)
        }

        if (cluster_rows) {
          nrow <- nrow(df_sub)
          d_sub <- stats::dist(df_sub, 
                               method = clustering_distance_rows, 
                               p = p_dist_rows)
          
          max_d_row <- suppressWarnings(max(d_sub, na.rm = TRUE))
          
          if ((!length(d_sub)) || is.infinite(max_d_row)) {
            h_sub <- FALSE
          }
          else {
            d_sub[is.na(d_sub)] <- .5 * max_d_row
            
            if (nrow <= 2L) {
              h_sub <- FALSE
            } 
            else {
              h_sub <- tryCatch(
                hclust(d_sub, hc_method_rows),
                error = function(e) 1L
              )
              
              if (class(h_sub) != "hclust" && h_sub == 1L) {
                warning("No row clustering for subtree: ", cluster_id)
                h_sub <- FALSE
              }
            }
          }
        }
        
        if (cluster_cols) {
          t_df_sub <- t(df_sub)
          
          nrow_trans <- nrow(t_df_sub)
          d_sub_col <- stats::dist(t_df_sub, 
                                   method = clustering_distance_cols, 
                                   p = p_dist_cols)
          
          max_d_col <- suppressWarnings(max(d_sub_col, na.rm = TRUE))
          
          if (!length(d_sub_col)) 
            next
          
          if (is.infinite(max_d_col)) {
            v_sub <- FALSE
          } 
          else {
            d_sub_col[is.na(d_sub_col)] <- .5 * max_d_col
  
            if (nrow_trans <= 2L) {
              v_sub <- FALSE
            } 
            else {
              v_sub <- tryCatch(
                hclust(d_sub_col, hc_method_cols),
                error = function(e) 1
              )
              
              if ((class(v_sub) != "hclust") && (v_sub == 1)) {
                warning("No column clustering for subtree: ", cluster_id)
                v_sub <- FALSE
              }
            }  
          }
        }
        
        if (nrow <= 150L) {
          cellheight <- 5
          fontsize_row <- 5
          show_rownames <- TRUE
          height <- NA
        } 
        else {
          cellheight <- NA
          fontsize_row <- NA
          show_rownames <- FALSE
          height <- NA
        }
        
        try(
          pheatmap(
            mat = df_sub,
            main = paste("Cluster", cluster_id),
            cluster_rows = h_sub,
            # cluster_cols = v_sub, 
            cluster_cols = h_cols, # use whole-tree hierarchy
            show_rownames = show_rownames,
            show_colnames = TRUE,
            annotation_col = annotation_col,
            annotation_row = annotation_row, 
            color = mypalette,
            breaks = color_breaks,
            cellwidth = 14,
            cellheight = cellheight,
            fontsize_row = fontsize_row,
            height = height, 
            annotation_colors = annotation_colors,
            filename = file.path(filepath, "Subtrees", fn_prefix, 
                                 paste0("Subtree_", cutree_rows, "-", cluster_id, ".",
                                        fn_suffix))
          )
        )

        rm(list = c("d_sub", "h_sub"))
      }
    }
  }
}


#'Visualization of heat maps
#'
#'\code{pepHM} applies \code{\link[stats]{dist}} and \code{\link[stats]{hclust}} 
#' for the visualization of the heat maps of peptide \code{log2FC} via 
#' \code{\link[pheatmap]{pheatmap}}.
#'
#'@rdname prnHM
#'
#'@import purrr
#'@export
pepHM <- function (col_select = NULL, col_order = NULL, col_benchmark = NULL,
                   scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                   rm_allna = TRUE, 
                   df = NULL, filepath = NULL, filename = NULL,
                   annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                   xmin = -1, xmax = 1, xmargin = 0.1, 
                   hc_method_rows = "complete", hc_method_cols = "complete", 
                   p_dist_rows = 2, p_dist_cols = 2, 
                   
                   x = NULL, p = NULL, method = NULL, 
                   diag = NULL, upper = NULL, 
                   annotation_col = NULL, annotation_row = NULL, 
                   clustering_method = NULL, ...) 
{
  old_opts <- options()
  options(warn = 1, warnPartialMatchArgs = TRUE)
  on.exit(options(old_opts), add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)

  id <- tryCatch(
    match_call_arg(normPSM, group_psm_by), error = function(e) "pep_seq_mod")

  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  col_benchmark <- rlang::enexpr(col_benchmark)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  hc_method_rows <- rlang::as_string(rlang::enexpr(hc_method_rows))
  hc_method_cols <- rlang::as_string(rlang::enexpr(hc_method_cols))

  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_order = !!col_order, 
            col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename, 
            anal_type = "Heatmap")(xmin = xmin, xmax = xmax, xmargin = xmargin, 
                                   annot_cols = annot_cols, 
                                   annot_colnames = annot_colnames, 
                                   annot_rows = annot_rows, 
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



#'Visualization of heat maps
#'
#'\code{prnHM} applies \code{\link[stats]{dist}} and \code{\link[stats]{hclust}} 
#' for the visualization of the heat maps of protein \code{log2FC} via 
#' \code{\link[pheatmap]{pheatmap}}.
#'
#'Data rows without non-missing pairs will result in NA distances in inter-row
#'dissimilarities (\code{\link[stats]{dist}}). At \code{complet_cases = TRUE},
#'the data subset that are complete without missing values will be used. At
#'\code{impute_na = TRUE}, all data rows will be used with NA imputation (see
#'\code{\link{prnImp}}). At the default of \code{complet_cases = FALSE} and
#'\code{impute_na = FALSE}, NA distances will be arbitrarily replaced with the
#'mean value of the row-distance matrix for hierarchical row clustering.
#'
#'Similar to data rows, NA distances in data columns will be replaced with the
#'mean value of the column-distance matrix.
#'
#'To avoid memory failure, row aggregation using the \code{kmeans_k} option
#'(\code{\link[pheatmap]{pheatmap}}) may be considered for large data sets.
#'
#'
#'@inheritParams  prnEucDist
#'@inheritParams normPSM
#'@inheritParams prnCorr_logFC
#'@param hc_method_rows A character string; the same agglomeration method for 
#'\code{\link[stats]{hclust}} of data rows. The default is \code{complete}. 
#'@param hc_method_cols A character string; similar to \code{hc_method_rows} 
#'but for column data.
#'@param  col_benchmark Not used.
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used. The default is FALSE.
#'@param annot_rows A character vector of column keys that can be found from
#'  input files of \code{Peptide.txt}, \code{Protein.txt} etc. The values
#'  under the selected keys will be used to color-code peptides or proteins on
#'  the side of the indicated plot. The default is NULL without row annotation.
#'@param xmin  Numeric; the minimum x at a log2 scale. The default is -1.
#'@param xmax  Numeric; the maximum  x at a log2 scale. The default is 1.
#'@param xmargin  Numeric; the margin in heat scales. The default is 0.1.
#'@param p_dist_rows Numeric; the power of the Minkowski distance in the measures 
#'  of row \code{\link[stats]{dist}} at \code{clustering_distance_rows = "minkowski"}. 
#'  The default is 2.
#'@param p_dist_cols Numeric; similar to \code{p_dist_rows} but for column data.
#'@param x Dummy argument to avoid incurring the corresponding argument in
#'  \link[stats]{dist} by partial argument matches.
#'@param p Dummy argument to avoid incurring the corresponding argument in
#'  \link[stats]{dist} by partial argument matches.
#'@param method Dummy argument to avoid incurring the corresponding argument in
#'  \link[stats]{dist} by partial argument matches.
#'@param diag Dummy argument to avoid incurring the corresponding argument in
#'  \link[stats]{dist} by partial argument matches.
#'@param upper Dummy argument to avoid incurring the corresponding argument in
#'  \link[stats]{dist} by partial argument matches.
#'@param annotation_col Dummy argument to avoid incurring the corresponding
#'  argument in \link[pheatmap]{pheatmap}.
#'@param annotation_row Dummy argument to avoid incurring the corresponding
#'  argument in \link[pheatmap]{pheatmap}.
#'@param clustering_method Dummy argument to avoid incurring the corresponding
#'  argument in \link[pheatmap]{pheatmap}.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. Each statement contains
#'  to a list of logical expression(s). The \code{lhs} needs to start with
#'  \code{filter_}. The logical condition(s) at the \code{rhs} needs to be
#'  enclosed in \code{exprs} with round parenthesis. For example, \code{pep_len}
#'  is a column key in \code{Peptide.txt}. The statement \code{filter_peps_at =
#'  exprs(pep_len <= 50)} will remove peptide entries with \code{pep_len > 50}.
#'  See also \code{\link{pepHist}}, \code{\link{normPSM}}. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. The \code{lhs} needs to start
#'  with \code{arrange_}. The expression(s) at the \code{rhs} needs to be
#'  enclosed in \code{exprs} with round parenthesis. For example,
#'  \code{arrange_peps_by = exprs(gene, prot_n_pep)} will arrange entries by
#'  \code{gene}, then by \code{prot_n_pep}. \cr \cr Additional parameters for
#'  plotting: \cr \code{width}, the width of plot \cr \code{height}, the height
#'  of plot \cr \cr Additional arguments for \code{\link[pheatmap]{pheatmap}}:
#'  \cr \code{cluster_rows, clustering_method, clustering_distance_rows}... \cr
#'  \cr Notes about \code{pheatmap}: \cr \code{annotation_col} disabled; instead
#'  use keys indicated in \code{annot_cols} \cr \code{annotation_row} disabled;
#'  instead use keys indicated in \code{annot_rows} \cr \code{clustering_method}
#'  breaks into \code{hc_method_rows} for row data and \code{hc_method_cols} for
#'  column data \cr \code{clustering_distance_rows = "minkowski"} allowed
#'  together with the powder of \code{p_dist_rows} and/or \code{p_dist_cols}
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
#'@example inst/extdata/examples/prnHM_.R
#'
#'@return Heat maps and optional sub trees.
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
prnHM <- function (col_select = NULL, col_order = NULL, col_benchmark = NULL,
                   scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                   rm_allna = TRUE, 
                   df = NULL, filepath = NULL, filename = NULL, 
                   annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                   xmin = -1, xmax = 1, xmargin = 0.1, 
                   hc_method_rows = "complete", hc_method_cols = "complete", 
                   p_dist_rows = 2, p_dist_cols = 2, 
                   
                   x = NULL, p = NULL, method = NULL, 
                   diag = NULL, upper = NULL, 
                   annotation_col = NULL, annotation_row = NULL, 
                   clustering_method = NULL, ...) 
{
  ## incorrect match of `x` to `xmargin` if `x` from dot-dot-dot
  ## correct match if `x` defined explictly in `prnHM`
  # prnHM(xmin = -1, xmax = 1, x = df)
  
  ## otherwise need all arguments starting with "x" in `prnHM`
  # prnHM(xmin = -1, xmax = 1, xmargin = .1, x = df)

  old_opts <- options()
  options(warn = 1, warnPartialMatchArgs = TRUE)
  on.exit(options(old_opts), add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)

  id <- tryCatch(
    match_call_arg(normPSM, group_pep_by), error = function(e) "gene")

  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  col_benchmark <- rlang::enexpr(col_benchmark)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  hc_method_rows <- rlang::as_string(rlang::enexpr(hc_method_rows))
  hc_method_cols <- rlang::as_string(rlang::enexpr(hc_method_cols))

  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_order = !!col_order, 
            col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename, 
            anal_type = "Heatmap")(xmin = xmin, 
                                   xmax = xmax, 
                                   xmargin = xmargin, 
                                   annot_cols = annot_cols, 
                                   annot_colnames = annot_colnames, 
                                   annot_rows = annot_rows, 
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


