#' Prepares data for analysis
#'
#' \code{prepDM} prepares a minimal adequate data frame for subsequent analysis.
#'
#' @param df A data frame containing only numeric values.
#' @param id The name of unqiue identifiers.
#' @param scale_log2r Logical; if TRUE, rescales \code{log2FC} to the same
#'   scale.
#' @param sub_grp Numeric.  A list of sample IDs that will be used in subsequent
#'   analysis.
#' @return A data frame tailored for subsequent analysis.
#'
#' @examples
#' tempData <- prepDM(df, entrez, scale_log2r, sub_grp = label_scheme_sub$Sample_ID)
#'
#' \dontrun{
#' }
#' @import dplyr
#' @importFrom magrittr %>%
prepDM <- function(df, id, scale_log2r, sub_grp, type = "ratio", anal_type) {
  stopifnot(nrow(df) > 0)
  
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  id <- rlang::as_string(rlang::enexpr(id))
  if(anal_type %in% c("ESGAGE", "GSVA")) id <- "entrez"
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  # data filtration dominated by log2R, not Intensity
  df <- df %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
    reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id)
  
  Levels <- sub_grp %>%
    as.character(.) %>%
    .[!grepl("^Empty\\.[0-9]+", .)]
  
  dfR <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sub_grp)) %>%
    dplyr::select(which(not_all_zero(.))) %>% # reference will drop with single reference
    dplyr::select(Levels[Levels %in% names(.)]) # ensure the same order
  
  # dominated by log2R, no need to filter all-NA intensity rows
  dfI <- df %>%
    dplyr::select(grep("^N_I[0-9]{3}", names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sub_grp)) %>%
    dplyr::select(which(not_all_zero(.))) %>%
    dplyr::select(which(names(.) %in% names(dfR))) %>%
    dplyr::select(Levels[Levels %in% names(.)])
  
  tempI <- dfI
  rownames(tempI) <- paste(rownames(tempI), "Intensity", sep = "@")
  
  tempR <- dfR
  rownames(tempR) <- paste(rownames(tempR), "log2R", sep = "@")
  
  return(list(log2R = dfR, Intensity = dfI, IR = rbind(tempR, tempI)))
}


#' A wrapper around pheatmap
#'
#' @import dplyr rlang pheatmap
#' @importFrom magrittr %>%
my_pheatmap <- function(mat, filename, annotation_col, color, annotation_colors, breaks, ...) {
  mat <- rlang::enexpr(mat)
  filename <- rlang::enexpr(filename)
  annotation_col <- rlang::enexpr(annotation_col)
  color <- rlang::enexpr(color)
  annotation_colors <- rlang::enexpr(annotation_colors)
  breaks <- rlang::enexpr(breaks)
  
  dots <- rlang::enexprs(...)
  
  dots$mat <- NULL
  dots$filename <- NULL
  dots$annotation_col <- NULL
  dots$color <- NULL
  dots$annotation_colors <- NULL
  dots$breaks <- NULL
  
  ph_call <- rlang::expr(
    pheatmap(mat = !!mat, filename = !!filename, annotation_col = !!annotation_col, color = !!color,
             annotation_colors = !!annotation_colors, breaks = !!breaks, !!!dots))
  
  rlang::expr_print(ph_call)
  rlang::eval_bare(ph_call, env = caller_env())
}


#' Makes heat maps
#'
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
plotHM <- function(df, id, scale_log2r, col_benchmark, label_scheme_sub, filepath, filename,
                   complete_cases, xmin = -1, xmax = 1, x_margin = .1, 
                   annot_cols = NULL, annot_colnames = NULL, ...) {
  
  id <- rlang::as_string(rlang::enexpr(id))
  col_benchmark <- rlang::as_string(rlang::enexpr(col_benchmark))
  
  dots <- rlang::enexprs(...)
  
  if (is.null(dots$cluster_rows)) {
    cluster_rows <- TRUE
  } else {
    cluster_rows <- dots$cluster_rows
  }
  
  if(is.null(dots$clustering_distance_rows)) {
    clustering_distance_rows <- "euclidean"
  } else {
    clustering_distance_rows <- dots$clustering_distance_rows
  }
  
  # parameter(s) in "dots" being disabled for pheatmap()
  dots$annotation_col <- NULL
  nm_idx <- names(dots) %in% c("df", "id", "scale_log2r", "annot_kinases",
                               "complete_cases", "label_scheme_sub", "filepath", "filename")
  dots[nm_idx] <- NULL
  
  fn_prx <- gsub("\\..*$", "", filename)
  fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
  
  dir.create(file.path(filepath, "Subtrees"), recursive = TRUE, showWarnings = FALSE)
  
  x_label <- expression("Ratio ("*log[2]*")")
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  NorZ_ratios_to_ctrl <- paste("toCtrl", NorZ_ratios, sep = "_")
  
  n_color <- 500
  if (is.null(dots$breaks)) {
    color_breaks <- c(seq(from = xmin, -x_margin, length = n_color/2)[1:(n_color/2-1)],
                      seq(-x_margin, x_margin, length = 3),
                      seq(x_margin, xmax, length = n_color/2)[2:(n_color/2)])
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
  
  # mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(n_color)
  # mypalette <- colorRampPalette(c("darkorange3", "white", "darkblue"))(n_color)
  # image(matrix(1:n_color, nrow = n_color, ncol = 10),
  #   col = mypalette, xaxt = "n", yaxt = "n", useRaster = TRUE)
  
  acc_type <- load(file = file.path(dat_dir, "label_scheme.Rdata")) %>% find_acctype()
  
  if(acc_type == "refseq_acc") {
    df <- df %>% mutate(Species = gsub(".*\\[(.*)\\].*", "\\1", .$prot_desc))
  } else if (acc_type == "uniprot_id") {
    df <- df %>% mutate(Species = gsub(".*\\s+OS\\=(.*)\\s+GN.*", "\\1", .$prot_desc))
  }
  
  sample_ids <- label_scheme_sub$Sample_ID
  
  df <- df %>%
    dplyr::mutate_at(vars(grep("^pVal|^adjP", names(.))), as.numeric) %>%
    dplyr::mutate(Mean_log10Int = log10(rowMeans(.[, grepl("^I[0-9]{3}", names(.))],
                                                 na.rm = TRUE))) %>%
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}", names(.))), ~setHMlims(., xmin, xmax)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
    reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id)
  
  dfR <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sample_ids)) %>%
    dplyr::select(which(not_all_zero(.))) %>%
    dplyr::select(as.character(sample_ids)) # ensure the same order
  
  df <- df %>%
    dplyr::select(-grep("log2_R[0-9]{3}", names(.))) %>%
    dplyr::bind_cols(., dfR) %>%
    `rownames<-`(.[[id]])
  
  if(is.null(annot_cols)) {
    annotation_col <- NA
  } else {
    annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)
  }
  
  if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
    colnames(annotation_col) <- annot_colnames
  }
  
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
  
  if(cluster_rows) {
    d <- dist(df_hm, method = clustering_distance_rows)
    d[is.na(d)] <- .5 * max(d, na.rm = TRUE)
    h <- hclust(d)
    dots$cluster_rows <- h
    
    rm(d, h)
  } else {
    dots$cluster_rows <- FALSE
  }
  
  filename <- gg_imgname(filename)
  
  p <- my_pheatmap(
    mat = df_hm,
    filename = file.path(filepath, filename),
    annotation_col = annotation_col,
    color = mypalette,
    annotation_colors = annotation_colors,
    breaks = color_breaks,
    !!!dots
  )
  
  # when cutree_rows is not NA
  cutree_rows <- eval(dots$cutree_rows, env = caller_env())
  
  if(!is.null(cutree_rows) & cluster_rows) {
    if(is.numeric(cutree_rows)) {
      Cluster <- data.frame(Cluster = cutree(p$tree_row, k = cutree_rows)) %>%
        dplyr::mutate(!!id := rownames(.)) %>%
        dplyr::left_join(df, by = id)
      
      write.csv(Cluster,
                file.path(filepath, "Subtrees", paste0(fn_prx, " n-", cutree_rows, "_subtrees.csv")),
                row.names = FALSE)
      
      Cluster <- Cluster %>%
        tibble::column_to_rownames(var = id)
      
      for(cluster_id in unique(Cluster$Cluster)) {
        df_sub <- Cluster[Cluster$Cluster == cluster_id, names(Cluster) %in% sample_ids]
        
        if(complete_cases) {
          df_sub <- df_sub %>%
            tibble::rownames_to_column(id) %>%
            dplyr::filter(complete.cases(.[, names(.) %in% sample_ids])) %>%
            tibble::column_to_rownames(id)
        }
        
        if(cluster_rows) {
          d_sub <- dist(df_sub, method = clustering_distance_rows)
          d_sub[is.na(d_sub)] <- .5 * max(d_sub, na.rm = TRUE)
          
          if (nrow(df_sub) < 2) {
            h_sub <- FALSE
          } else {
            h_sub <- hclust(d_sub)
          }
        }
        
        pheatmap(
          mat = df_sub,
          main = paste("Cluster", cluster_id),
          cluster_rows = h_sub,
          show_rownames = TRUE,
          show_colnames = TRUE,
          annotation_col = annotation_col,
          color = mypalette,
          breaks = color_breaks,
          cellwidth = 14,
          fontsize_row = ceiling(200/nrow(Cluster %>% dplyr::filter(Cluster == cluster_id))),
          annotation_colors = annotation_colors,
          filename = file.path(filepath, "Subtrees",
                               paste0("Subtree_", cutree_rows, "-", cluster_id, ".", fn_suffix))
        )
        
        rm(d_sub, h_sub)
      }
    }
  }
  
  # to the benchmark group
  if(sum(!is.na(label_scheme_sub[[col_benchmark]])) > 0) {
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
    
    if(cluster_rows) {
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
      filename = file.path(filepath, paste0(fn_prx, "_toCtrl.", fn_suffix)),
      annotation_col = annotation_col,
      color = mypalette,
      annotation_colors = annotation_colors,
      breaks = color_breaks,
      !!!dots
    )
  }
}


#'Plots heat maps
#'
#'\code{proteoHM} produces the heat map visualization of \code{log2FC} for
#'proteins or peptides data.
#'
#'Data columns with complete missing values will be removed prior to
#'hierarchical column clustering. Data rows without non-missing pairs will
#'result in NA distances in inter-row dissimilarities
#'(\code{\link[stats]{dist}}). At \code{complet_cases = TRUE}, the data subset
#'that are complete without missing values will be used. At \code{impute_na =
#'TRUE}, all data rows will be used with NA imputation (\code{\link{imputeNA}}).
#'At the default of \code{complet_cases = FALSE} and \code{impute_na = FALSE},
#'NA distances will be arbitrarily replaced with the mean value in the
#'row-disance matrix for hierarchical row clustering.
#'
#'To avoid memory failure, row aggregation using the \code{kmeans_k} option
#'(\code{\link[pheatmap]{pheatmap}}) may be considered for large data sets.
#'
#'The function matches the current \code{id} to those in the latest \code{calls}
#'to \code{\link{normPep}} or \code{\link{normPrn}}.  For example, if
#'\code{pep_seq} was used in \code{\link{normPep}}, the current \code{id =
#'pep_seq_mod} will be matched to \code{id = pep_seq}.
#'
#'@inheritParams  proteoEucDist
#'@param  col_benchmark Not used.
#'@param impute_na Logical; if TRUE, imputes missing values.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used for visualization.
#'@param ... More parameters for plotting: \cr \code{xmin}, the minimum \eqn{x}
#'  at a log2 scale; the default is -1 \cr \code{xmax}, the maximum \eqn{x} at a
#'  log2 scale; the default is +1 \cr \code{x_margin}, the margin in heat
#'  scales; the default is 0.1 \cr \code{width}, the width of plot \cr
#'  \code{height}, the height of plot \cr additional arguments inherited from
#'  \code{\link[pheatmap]{pheatmap}}.
#'@return Heat map images.
#'
#' @examples
#'prnHM(
#'  scale_log2r = TRUE,
#'  xmin = -1,
#'  xmax = 1,
#'  x_margin = 0.1,
#'  annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'  annot_colnames = c("Peptide_Yield", "Batch", "Phenotype"),
#'  cluster_rows = TRUE,
#'  cutree_rows = 6,
#'  show_rownames = FALSE,
#'  show_colnames = TRUE,
#'  fontsize_row = 3,
#'  cellwidth = 14,
#'  width = 10,
#'  height = 12
#')
#'
#'prnHM(
#'  scale_log2r = TRUE,
#'  annot_cols = c("Group"),
#'  cluster_rows = TRUE,
#'  clustering_distance_rows  = "maximum",
#'  cutree_rows = 6,
#'  show_rownames = FALSE,
#'  show_colnames = TRUE,
#'  fontsize_row = 3,
#'  cellwidth = 14,
#'  width = 10,
#'  height = 12
#')
#'
#'prnHM(
#'  scale_log2r = TRUE,
#'  annot_cols = c("Group"),
#'  cluster_rows = FALSE, # no row clustering
#'  clustering_distance_rows  = "maximum",
#'  cutree_rows = 6, # will be overruled at 'cluster_rows = FALSE'
#'  show_rownames = FALSE,
#'  show_colnames = TRUE,
#'  fontsize_row = 3,
#'  cellwidth = 14,
#'  width = 10,
#'  height = 12
#')
#'
#' \dontrun{
#' }
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoHM <- function (id = gene, col_select = NULL, col_benchmark = NULL,
                      scale_log2r = TRUE,impute_na = FALSE, complete_cases = FALSE,
                      df = NULL, filepath = NULL, filename = NULL,
                      xmin = -1, xmax = 1, x_margin = 0.1, 
                      annot_cols = NULL, annot_colnames = NULL, ...) {
  
  # scale_log2r <- match_logi_gv(scale_log2r)
  
  id <- rlang::enexpr(id)
  col_select <- rlang::enexpr(col_select)
  col_benchmark <- rlang::enexpr(col_benchmark)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, impute_na = impute_na, df = !!df, filepath = !!filepath,
            filename = !!filename, anal_type = "Heatmap")(complete_cases = complete_cases,
                                                          xmin = xmin, xmax = xmax,
                                                          x_margin = x_margin,
                                                          annot_cols = annot_cols, 
                                                          annot_colnames = annot_colnames, ...)
}


#'Visualizes peptide heat maps
#'
#'\code{pepHM} is a wrapper function of \code{\link{proteoHM}} for peptide data
#'
#'@rdname proteoHM
#'
#'@export
pepHM <- function (...) {
  proteoHM(id = pep_seq, ...)
}


#'Visualizes protein heat maps
#'
#'\code{prnHM} is a wrapper function of \code{\link{proteoHM}} for protein data
#'
#'@rdname proteoHM
#'
#'@export
prnHM <- function (...) {
  proteoHM(id = gene, ...)
}


#' Plots kinase heat maps
#'
#' Specialized for kinase heat maps
#'
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
#' @export
plotKinHM <- function(id, scale_log2r, col_benchmark, label_scheme_sub, df, filepath, filename,
                      complete_cases, xmin = -1, xmax = 1, x_margin = .1, 
                      annot_cols = NULL, annot_colnames = NULL, ...) {
  
  id <- rlang::as_string(rlang::enexpr(id))
  col_benchmark <- rlang::as_string(rlang::enexpr(col_benchmark))
  
  dots <- rlang::enexprs(...)
  
  # parameter(s) in "dots" disabled for pheatmap()
  dots$annotation_col <- NULL
  nm_idx <- names(dots) %in% c("df", "id", "scale_log2r", "annot_kinases", "complete_cases",
                               "label_scheme_sub", "filepath", "filename")
  dots[nm_idx] <- NULL
  
  fn_prx <- gsub("\\..*$", "", filename)
  fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
  
  dir.create(file.path(filepath, "Subtrees"), recursive = TRUE, showWarnings = FALSE)
  
  x_label <- expression("Ratio ("*log[2]*")")
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  NorZ_ratios_to_ctrl <- paste("toCtrl", NorZ_ratios, sep = "_")
  
  n_color <- 500
  if (is.null(dots$breaks)) {
    color_breaks <- c(seq(from = xmin, -x_margin, length = n_color/2)[1:(n_color/2-1)],
                      seq(-x_margin, x_margin, length = 3),
                      seq(x_margin, xmax, length = n_color/2)[2:(n_color/2)])
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
  
  acc_type <- load(file = file.path(dat_dir, "label_scheme.Rdata")) %>% find_acctype()
  
  if(acc_type == "refseq_acc") {
    df <- df %>% mutate(Species = gsub(".*\\[(.*)\\].*", "\\1", .$prot_desc))
  } else if (acc_type == "uniprot_id") {
    df <- df %>% mutate(Species = gsub(".*\\s+OS\\=(.*)\\s+GN.*", "\\1", .$prot_desc))
  }
  
  sample_ids <- label_scheme_sub$Sample_ID
  
  kin_levels <- c("TK", "TKL", "STE", "CK1", "AGC", "CAMK", "CMGC", "Atypical", "Other", "RGC",
                  "Unclassified", "FAM20", "Lipid")
  
  df <- df %>%
    dplyr::mutate_at(vars(grep("^pVal|^adjP", names(.))), as.numeric) %>%
    dplyr::mutate(Mean_log10Int = log10(rowMeans(.[, grepl("^I[0-9]{3}", names(.))],
                                                 na.rm = TRUE))) %>%
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}", names(.))), ~setHMlims(., xmin, xmax)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
    reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id) %>%
    dplyr::mutate(kin_class = factor(kin_class, levels = kin_levels)) %>%
    dplyr::arrange(kin_class, !!rlang::sym(id))
  
  dfR <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sample_ids)) %>%
    dplyr::select(which(not_all_zero(.))) %>%
    dplyr::select(as.character(sample_ids)) # ensure the same order
  
  df <- df %>%
    dplyr::select(-grep("log2_R[0-9]{3}", names(.))) %>%
    dplyr::bind_cols(., dfR) %>%
    `rownames<-`(.[[id]])
  
  annotation_row <- df %>%
    dplyr::select(id, kin_class) %>%
    dplyr::mutate(kin_class = factor(kin_class, levels =
                                       c("TK", "TKL", "STE", "CK1", "AGC", "CAMK", "CMGC",
                                         "Atypical", "Other", "RGC", "Unclassified",
                                         "FAM20", "Lipid"))) %>%
    tibble::column_to_rownames(id)
  
  if (is.null(annot_cols)) {
    annotation_col <- NA
  } else {
    annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)
  }
  
  if (is.null(dots$annotation_colors)) {
    annotation_colors <- setHMColor(annotation_col)
  } else if (is.na(dots$annotation_colors)) {
    annotation_colors <- NA
  } else {
    annotation_colors <- eval(dots$annotation_colors, env = caller_env())
  }
  
  # complete_cases <- FALSE
  if (complete_cases) {
    df_hm <- df %>%
      dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
  } else {
    df_hm <- df %>%
      dplyr::mutate_if(is.numeric, funs(replace(., is.na(.), 0)))
  }
  
  df_hm <- df_hm %>%
    `rownames<-`(.[[id]])	%>%
    dplyr::select(which(names(.) %in% sample_ids))
  
  p <- my_pheatmap(
    mat = df_hm,
    filename = file.path(filepath, filename),
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    color = mypalette,
    annotation_colors = annotation_colors,
    breaks = color_breaks,
    !!!dots)
  
  # to the benchmark group
  if(sum(!is.na(label_scheme_sub[[col_benchmark]])) > 0) {
    dfc <- ratio_toCtrl(df, !!id, label_scheme_sub, nm_ctrl = !!col_benchmark)
    
    if (complete_cases) {
      dfc_hm <- dfc %>%
        dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
    } else {
      dfc_hm <- dfc %>%
        dplyr::filter(rowSums(!is.na(.[, names(.) %in% sample_ids])) > 0) %>%
        dplyr::mutate_if(is.numeric, funs(replace(., is.na(.), 0)))
    }
    
    dfc_hm <- dfc_hm %>%
      `rownames<-`(.[[id]])	%>%
      dplyr::select(which(names(.) %in% sample_ids))
    
    p <- my_pheatmap(
      mat = dfc_hm,
      filename = file.path(filepath, paste0(fn_prx, "_toCtrl.", fn_suffix)),
      annotation_col = annotation_col,
      annotation_row = annotation_row,
      color = mypalette,
      annotation_colors = annotation_colors,
      breaks = color_breaks,
      !!!dots)
  }
}


#' Visualizes kinase heat maps
#'
#' Specialized for kinase heat maps
#'
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
#' @export
proteoKinHM <- function (id = gene, col_select = NULL, col_benchmark = NULL, scale_log2r = TRUE,
                         df = NULL, filepath = NULL, filename = NULL, complete_cases = FALSE,
                         impute_na = FALSE, anal_type = "Heatmap", xmin = -1, xmax = 1,
                         x_margin = 0.1, annot_cols = NULL, ...) {
  
  # scale_log2r <- match_logi_gv(scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_benchmark <- rlang::enexpr(col_benchmark)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  if(is.null(col_select)) {
    col_select <- rlang::expr(Select)
    if(sum(!is.na(label_scheme$Select)) == 0)
      stop("Column \'", rlang::as_string(col_select),
           "\' is either empty or missing from the \'label scheme\' file.\n",
           "\tEnter sample names under the \'", rlang::as_string(col_select),
           "\' column for informatic analysis.", call. = TRUE)
  } else if(col_select == rlang::expr(Sample_ID)) {
    stop(paste0("\'Sample_ID\' is a reserved column name.\n",
                "\tChoose or add a different column key to indicate samples for informatic analysis."),
         call. = TRUE)
  } else if(is.null(label_scheme[[col_select]])) {
    stop("Column \'", rlang::as_string(col_select), "\' is empty.\n",
         "\tEnter sample names under the \'", rlang::as_string(col_select),
         "\' column for informatic analysis.", call. = TRUE)
  }
  
  if(is.null(col_benchmark)) {
    col_benchmark <- rlang::expr(Benchmark)
    if(sum(!is.na(label_scheme$Benchmark)) == 0)
      warning("Default column name \'", rlang::as_string(col_benchmark),
              "\' is empty or missing.", call. = TRUE)
  } else if(col_benchmark == rlang::expr(Sample_ID)) {
    stop(paste0("Column \'Sample_ID\' is reserved; choose a different key.\n"))
  } else if(is.null(label_scheme[[col_benchmark]])) {
    stop("Column name \'", rlang::as_string(col_benchmark), "\' is missing.", call. = TRUE)
  }
  
  id <- rlang::as_string(rlang::enexpr(id))
  id <- match_identifier(id)
  
  dots <- enexprs(...)
  
  if(!id %in% c("pep_seq_mod", "prot_acc", "gene"))
    stop("Unrecognized 'id'; needs to be one of pep_seq_mod, prot_acc or gene", call. = TRUE)
  
  if (id %in% c("prot_acc", "gene")) {
    data_type <- "Protein"
  } else if (id %in% c("pep_seq_mod")) {
    data_type <- "Peptide"
  }
  
  if(is.null(filepath)) {
    filepath = file.path(dat_dir, data_type, anal_type)
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  }
  
  if(is.null(filename)) {
    fn_prx <- paste(data_type, anal_type, sep = "_")
    fn_suffix <- "png"
  } else {
    fn_prx <- gsub("\\..*$", "", filename)
    fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
  }
  
  if(is.null(df)) {
    err_msg <- "File doesn't exist"
    
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
      fn_imp <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
      fn_raw <- file.path(dat_dir, "Peptide", "Peptide.txt")
    } else if (id %in% c("prot_acc", "gene")) {
      fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
      fn_imp <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
      fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")
    }
    
    if(file.exists(fn_p)) {
      src_path <- fn_p
    } else {
      src_path <- ifelse(impute_na, fn_imp, fn_raw)
    }
    
    df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
                            comment.char = "#"), error = function(e) NA)
    
    if(!is.null(dim(df))) {
      message(paste("File loaded:", gsub("\\\\", "/", src_path)))
    } else {
      stop(paste("No such file or directory:", gsub("\\\\", "/", src_path)))
    }
  } else {
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      fn_raw <- file.path(dat_dir, "Peptide", df)
    } else if (id %in% c("prot_acc", "gene")) {
      fn_raw <- file.path(dat_dir, "Protein", df)
    }
    
    df <- tryCatch(read.csv(fn_raw, check.names = FALSE, header = TRUE, sep = "\t",
                            comment.char = "#"), error = function(e) NA)
    
    if(!is.null(dim(df))) {
      message(paste("File loaded:", gsub("\\\\", "/", fn_raw)))
    } else {
      stop(paste("Non-existed file or directory:", gsub("\\\\", "/", fn_raw)))
    }
  }
  
  reload_expts()
  
  label_scheme_sub <- label_scheme %>%
    dplyr::select(Sample_ID, TMT_Set, !!col_select) %>%
    dplyr::filter(!is.na(!!col_select))
  
  dfw <- prepDM(df, !!rlang::sym(id), scale_log2r, label_scheme_sub$Sample_ID, anal_type = anal_type)
  
  dir.create(file.path(filepath, "Kinases"), recursive = TRUE, showWarnings = FALSE)
  
  dfw_kinase <- df %>%
    dplyr::filter(kin_attr) %>%
    prepDM(!!rlang::sym(id), scale_log2r, label_scheme_sub$Sample_ID, anal_type = anal_type)
  
  plotKinHM(
    df = df[df$kin_attr, ],
    id = !!rlang::sym(id),
    scale_log2r = scale_log2r,
    col_benchmark = !!col_benchmark,
    label_scheme_sub = label_scheme_sub,
    filepath = file.path(filepath, "Kinases"),
    filename = paste0(fn_prx, "_Kinases.", fn_suffix),
    complete_cases = complete_cases,
    xmin = xmin,
    xmax = xmax,
    x_margin = x_margin,
    annot_cols = annot_cols,
    !!!dots)
}
