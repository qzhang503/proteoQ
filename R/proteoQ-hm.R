#' A wrapper of pheatmap
#'
#' @import dplyr rlang pheatmap
#' @importFrom magrittr %>%
my_pheatmap <- function(mat, filename, annotation_col, annotation_row, color, annotation_colors, breaks, ...) {
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
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
plotHM <- function(df, id, col_benchmark, label_scheme_sub, filepath, filename,
                   scale_log2r, complete_cases, 
                   annot_cols = NULL, annot_colnames = NULL, annot_rows = annot_rows, 
                   xmin = -1, xmax = 1, xmargin = .1, ...) {
                   
                   
  dir.create(file.path(filepath, "Subtrees"), recursive = TRUE, showWarnings = FALSE)
  
  id <- rlang::as_string(rlang::enexpr(id))
  col_benchmark <- rlang::as_string(rlang::enexpr(col_benchmark))
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]

  # needed defaults before calling `pheatmap`
  if (is.null(dots$cluster_rows)) {
    cluster_rows <- TRUE
  } else {
    cluster_rows <- dots$cluster_rows
  }
  
  if (is.null(dots$cluster_cols)) {
    cluster_cols <- TRUE
  } else {
    cluster_cols <- dots$cluster_cols
  }
  
  if (is.null(dots$clustering_distance_rows)) {
    clustering_distance_rows <- "euclidean"
  } else {
    clustering_distance_rows <- dots$clustering_distance_rows
  }
  
  if (is.null(dots$clustering_distance_cols)) {
    clustering_distance_cols <- "euclidean"
  } else {
    clustering_distance_cols <- dots$clustering_distance_cols
  }
  
  n_color <- 500
  if (is.null(dots$breaks)) {
    color_breaks <- c(seq(from = xmin, -xmargin, length = n_color/2)[1:(n_color/2-1)],
                      seq(-xmargin, xmargin, length = 3),
                      seq(xmargin, xmax, length = n_color/2)[2:(n_color/2)])
  } else if (is.na(dots$breaks)) {
    color_breaks <- NA
  } else {
    color_breaks <- eval(dots$breaks, env = caller_env())
  }
  
  if (is.null(dots$color)) {
    mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
  } else if (is.na(dots$color)) {
    mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  } else {
    mypalette <- eval(dots$color, env = caller_env())
  }
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  x_label <- expression("Ratio ("*log[2]*")")
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  NorZ_ratios_to_ctrl <- paste("toCtrl", NorZ_ratios, sep = "_")
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)

  sample_ids <- label_scheme_sub$Sample_ID
  
  df <- df %>%
    dplyr::mutate_at(vars(grep("^pVal|^adjP", names(.))), as.numeric) %>%
    dplyr::mutate(Mean_log10Int = log10(rowMeans(.[, grepl("^I[0-9]{3}", names(.))],
                                                 na.rm = TRUE))) %>%
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}", names(.))), ~ setHMlims(.x, xmin, xmax)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>% 
    reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id) 
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  if (nrow(df) == 0) stop("No rows available after data filtratin.", call. = FALSE)

  dfR <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sample_ids)) %>%
    # dplyr::select(which(not_all_zero(.))) %>%
    dplyr::select(as.character(sample_ids)) # ensure the same order

  df <- df %>%
    dplyr::select(-grep("log2_R[0-9]{3}", names(.))) %>%
    dplyr::bind_cols(., dfR) %>% 
    dplyr::filter(!is.na(.[[id]])) %>% 
    `rownames<-`(.[[id]])
  
  # generate `annotation_col` from keys in `annot_col`, et al.
  if (is.null(annot_cols)) {
    annotation_col <- NA
  } else {
    annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)
  }
  
  if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
    colnames(annotation_col) <- annot_colnames
  }
  
  if (is.null(annot_rows)) {
    annotation_row <- NA
  } else {
    annotation_row <- df %>% dplyr::select(annot_rows)
  }
  
  # column annotations
  if (is.null(dots$annotation_colors)) {
    annotation_colors <- setHMColor(annotation_col)
  } else if (is.na(dots$annotation_colors)) {
    annotation_colors <- NA
  } else {
    annotation_colors <- eval(dots$annotation_colors, env = caller_env())
  }
  
  if (complete_cases) {
    df_hm <- df %>%
      dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
  } else {
    df_hm <- df
  }

  df_hm <- df_hm %>%
    `rownames<-`(.[[id]])	%>%
    dplyr::select(which(names(.) %in% sample_ids))
  
  if (cluster_rows) {
    d <- dist(df_hm, method = clustering_distance_rows)
    d[is.na(d)] <- .5 * max(d, na.rm = TRUE)
    h <- hclust(d)
    dots$cluster_rows <- h
    rm(d, h)
  } else {
    dots$cluster_rows <- FALSE
  }
  
  if (cluster_cols) {
    d_cols <- dist(t(df_hm), method = clustering_distance_cols)
    d_cols[is.na(d_cols)] <- .5 * max(d_cols, na.rm = TRUE)
    h_cols <- hclust(d_cols)
    dots$cluster_cols <- h_cols
    # rm(d_cols, h_cols) # h_cols also for subtrees
  } else {
    dots$cluster_cols <- FALSE
  }
  
  filename <- gg_imgname(filename)
  
  # form `annotation_col` and `annotation_row` from `annot_col` and `annot_row` 
  dots <- dots %>% 
    .[! names(.) %in% c("mat", "filename", "annotation_col", "annotation_row", 
                        "clustering_distance_rows", "clustering_distance_cols", 
                        "color", "annotation_colors", "breaks")]
  
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
  cutree_rows <- eval(dots$cutree_rows, env = caller_env())
  df <- df %>% dplyr::mutate(!!id := as.character(!!rlang::sym(id)))
  
  if (!is.null(cutree_rows) & cluster_rows) {
    if (is.numeric(cutree_rows)) {
      Cluster <- data.frame(Cluster = cutree(p$tree_row, k = cutree_rows)) %>%
        dplyr::mutate(!!id := rownames(.)) %>%
        dplyr::left_join(df, by = id) %T>% 
        write.csv(file.path(filepath, "Subtrees", paste0(fn_prefix, " n-", cutree_rows, "_subtrees.csv")), 
                  row.names = FALSE)

      Cluster <- Cluster %>%
        tibble::column_to_rownames(var = id)
      
      for (cluster_id in unique(Cluster$Cluster)) {
        df_sub <- Cluster[Cluster$Cluster == cluster_id, names(Cluster) %in% sample_ids]

        if (complete_cases) {
          df_sub <- df_sub %>%
            tibble::rownames_to_column(id) %>%
            dplyr::filter(complete.cases(.[, names(.) %in% sample_ids])) %>%
            tibble::column_to_rownames(id)
        }

        if (cluster_rows) {
          d_sub <- dist(df_sub, method = clustering_distance_rows)
          if (length(d_sub) == 0) {
            h_sub <- FALSE
            nrow <- nrow(df_sub)
          } else {
            d_sub[is.na(d_sub)] <- .5 * max(d_sub, na.rm = TRUE)
            
            nrow <- nrow(df_sub)
            if (nrow <= 2) {
              h_sub <- FALSE
            } else {
              h_sub <- hclust(d_sub)
            }    
          }
        }
        
        if (cluster_cols) {
          d_sub_col <- dist(t(df_sub), method = clustering_distance_cols)
          if (length(d_sub_col) == 0) next
          
          d_sub_col[is.na(d_sub_col)] <- .5 * max(d_sub_col, na.rm = TRUE)
          nrow_trans <- nrow(t(df_sub))
          
          if (nrow_trans <= 2) {
            v_sub <- FALSE
          } else {
            v_sub <- hclust(d_sub_col)
          }
        }
        
        if (nrow <= 150) {
          cellheight <- 5
          fontsize_row <- 5
          show_rownames <- TRUE
          height <- NA
        } else {
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
            cluster_cols = h_cols, 
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
            filename = file.path(filepath, "Subtrees",
                                 paste0("Subtree_", cutree_rows, "-", cluster_id, ".", fn_suffix))
          )
        )

        rm(d_sub, h_sub)
      }
    }
  }
  
  # to the benchmark group
  # not currently used
  if (sum(!is.na(label_scheme_sub[[col_benchmark]])) > 0) {
    dfc <- ratio_toCtrl(df, !!rlang::sym(id), label_scheme_sub, nm_ctrl = !!col_benchmark)
    
    if (complete_cases) {
      dfc_hm <- dfc %>%
        dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
    } else {
      dfc_hm <- dfc %>%
        dplyr::filter(rowSums(!is.na(.[, names(.) %in% sample_ids])) > 0)
    }

    dfc_hm <- dfc_hm %>%
      `rownames<-`(.[[id]])	%>%
      dplyr::select(which(names(.) %in% sample_ids))
    
    if (cluster_rows) {
      d <- dist(dfc_hm, method = clustering_distance_rows)
      d[is.na(d)] <- .5  * max(d, na.rm = TRUE)
      h <- hclust(d)
      dots$cluster_rows <- h
      
      rm(d, h)
    } else {
      dots$cluster_rows <- FALSE
    }
    
    p <- my_pheatmap(
      mat = dfc_hm,
      filename = file.path(filepath, paste0(fn_prefix, "_toCtrl.", fn_suffix)),
      annotation_col = annotation_col,
      annotation_row = annotation_row, 
      color = mypalette,
      annotation_colors = annotation_colors,
      breaks = color_breaks,
      !!!dots
    )
  }
}


#'Visualization of heat maps
#'
#'\code{pepHM} visualizes the heat maps of peptide \code{log2FC}.
#'
#'@rdname prnHM
#'
#'@import purrr
#'@export
pepHM <- function (col_select = NULL, col_benchmark = NULL,
                   scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                   df = NULL, filepath = NULL, filename = NULL,
                   annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                   xmin = -1, xmax = 1, xmargin = 0.1, ...) {
  check_dots(c("id", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(rlang::is_double(xmin), 
            rlang::is_double(xmax), 
            rlang::is_double(xmargin))
  
  col_select <- rlang::enexpr(col_select)
  col_benchmark <- rlang::enexpr(col_benchmark)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, complete_cases = complete_cases,impute_na = impute_na, 
            df = !!df, filepath = !!filepath,
            filename = !!filename, anal_type = "Heatmap")(xmin = xmin, xmax = xmax, xmargin = xmargin,
                                                          annot_cols = annot_cols, 
                                                          annot_colnames = annot_colnames, 
                                                          annot_rows = annot_rows, ...)
}



#'Visualization of heat maps
#'
#'\code{prnHM} visualizes the heat maps of protein \code{log2FC}.
#'
#'Data rows without non-missing pairs will result in NA distances in inter-row
#'dissimilarities (\code{\link[stats]{dist}}). At \code{complet_cases = TRUE},
#'the data subset that are complete without missing values will be used. At
#'\code{impute_na = TRUE}, all data rows will be used with NA imputation (see
#'\code{\link{prnImp}}). At the default of \code{complet_cases = FALSE} and
#'\code{impute_na = FALSE}, NA distances will be arbitrarily replaced with the
#'mean value of the row-disance matrix for hierarchical row clustering.
#'
#'Similar to data rows, NA distances in data columns will be replaced with the
#'mean value of the column-distance matrix.
#'
#'To avoid memory failure, row aggregation using the \code{kmeans_k} option
#'(\code{\link[pheatmap]{pheatmap}}) may be considered for large data sets.
#'
#'
#'@inheritParams  prnEucDist
#'@param  col_benchmark Not used.
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used. The default is FALSE.
#'@param annot_rows A character vector of column keys that can be found from
#'  input files of \code{Peptide.txt}, \code{Protein.txt} et al. The values
#'  under the selected keys will be used to color-code peptides or proteins on
#'  the side of the indicated plot. The default is NULL without row annotation.
#'@param xmin  Numeric; the minimum x at a log2 scale. The default is -1.
#'@param xmax  Numeric; the maximum  x at a log2 scale. The default is 1.
#'@param xmargin  Numeric; the margin in heat scales. The default is 0.1.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  of data against the column keys in \code{Peptide.txt}, \code{Protein.txt}
#'  etc. Each statement contains to a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key present in \code{Mascot} peptide
#'  tables of \code{Peptide.txt}. The statement \code{filter_peps_at =
#'  exprs(pep_len <= 50)} will remove peptide entries with \code{pep_len > 50}.
#'  See also \code{\link{pepHist}}, \code{\link{normPSM}}. \cr \cr
#'  \code{arrange_}: Logical expression(s) for the row ordering of data. The
#'  \code{lhs} needs to start with \code{arrange_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{arrange_peps_by = exprs(gene, prot_n_pep)} will arrange
#'  entries by \code{gene}, then by \code{prot_n_pep}. \cr \cr Additional
#'  parameters for plotting: \cr \code{width}, the width of plot \cr
#'  \code{height}, the height of plot \cr \cr Additional arguments for
#'  \code{\link[pheatmap]{pheatmap}}, i.e., \code{cluster_rows}... \cr \cr Note
#'  arguments disabled for \code{pheatmap}: \cr \code{annotation_col}; instead
#'  use keys indicated in \code{annot_cols} \cr \code{annotation_row}; instead
#'  use keys indicated in \code{annot_rows}
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data
#'  normalization \cr
#'
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization
#'  \cr \code{\link{mergePep}} for extended examples in peptide data merging \cr
#'  \code{\link{standPep}} for extended examples in peptide data normalization
#'  \cr \code{\link{Pep2Prn}} for extended examples in peptide to protein
#'  summarization \cr \code{\link{standPrn}} for extended examples in protein
#'  data normalization. \cr \code{\link{purgePSM}} and \code{\link{purgePep}}
#'  for extended examples in data purging \cr \code{\link{pepHist}} and
#'  \code{\link{prnHist}} for extended examples in histogram visualization. \cr
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting
#'  MS file names \cr
#'
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr
#'
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation
#'  \cr \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization
#'  \cr
#'
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein
#'  significance pVals \cr \code{\link{gspaMap}} for mapping GSPA to volcano
#'  plot visualization \cr \code{\link{prnGSPAHM}} for heat map and network
#'  visualization of GSPA results \cr \code{\link{prnGSVA}} for gene set
#'  variance analysis \cr \code{\link{prnGSEA}} for data preparation for online
#'  GSEA. \cr
#'
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr
#'  \code{\link{pepPCA}} and \code{\link{prnPcA}} for PCA visualization \cr
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'  \code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for
#'  correlation plots \cr
#'
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for protein
#'  trend analysis and visualization \cr \code{\link{anal_pepNMF}},
#'  \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'  \code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for protein NMF
#'  analysis and visualization \cr
#'
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@example inst/extdata/examples/prnHM_.R
#'
#'@return Heat maps and optional sub trees.
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnHM <- function (col_select = NULL, col_benchmark = NULL,
                   scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                   df = NULL, filepath = NULL, filename = NULL,
                   annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                   xmin = -1, xmax = 1, xmargin = 0.1, ...) {
  check_dots(c("id", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(rlang::is_double(xmin), 
            rlang::is_double(xmax), 
            rlang::is_double(xmargin))
  
  col_select <- rlang::enexpr(col_select)
  col_benchmark <- rlang::enexpr(col_benchmark)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, complete_cases = complete_cases,impute_na = impute_na, 
            df = !!df, filepath = !!filepath,
            filename = !!filename, anal_type = "Heatmap")(xmin = xmin, xmax = xmax, xmargin = xmargin,
                                                          annot_cols = annot_cols, 
                                                          annot_colnames = annot_colnames, 
                                                          annot_rows = annot_rows, ...)
}


