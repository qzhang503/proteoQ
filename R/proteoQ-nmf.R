#' NMF analysis
#'
#' @import dplyr purrr rlang Biobase
#' @importFrom magrittr %>%
#' @importFrom NMF nmf
analNMF <- function(df, id, rank, nrun, seed, col_group, label_scheme_sub, 
                    scale_log2r, complete_cases, impute_na, 
                    filepath, filename, anal_type, ...) {

  stopifnot(nrow(label_scheme_sub) > 0)
  
  complete_cases <- to_complete_cases(complete_cases = complete_cases, impute_na = impute_na)
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    prepDM(id = !!id, scale_log2r = scale_log2r, 
           sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R

  col_group <- rlang::enexpr(col_group) # optional phenotypic information

	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% unique()
	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	stopifnot(length(fn_suffix) == 1)
	
  if (is.null(rank)) {
    rank <- c(4:8)
    nrun <- 5
    fn_prefix <- paste0(fn_prefix, rank)
  } else {
    stopifnot(all(rank >= 2) & all(rank %% 1 == 0))
    stopifnot(all(nrun >= 1) & all(nrun %% 1 == 0))
  }
  
  exprs_data <- data.matrix(2^df)
  
  pData <- label_scheme_sub %>%
    dplyr::filter(!is.na(!!col_group)) %>%
    dplyr::select(Sample_ID, !!col_group) %>%
    dplyr::rename(Group := !!col_group) %>%
    tibble::column_to_rownames("Sample_ID")
  
  metadata <- data.frame(labelDescription = c("Case/control status"), row.names = c("Group"))
  phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)
  experimentData <- new("MIAME", name = "Pierre Fermat", lab = "", contact = "",
                        title = "", abstract = "", url = "", other = list(notes = ""))
  exampleSet <- ExpressionSet(assayData = exprs_data, phenoData = phenoData,
                              experimentData = experimentData)

  purrr::walk(fn_prefix, ~ {
    rank <- gsub("^.*_rank(\\d+).*", "\\1", .x) %>% as.numeric()
    if (!is.null(seed)) set.seed(seed) else set.seed(sample(.Random.seed, 1))
    args <- c(list(x = exampleSet, rank = rank, nrun = nrun, seed = seed), dots)
    res_nmf <- do.call(NMF::nmf, args) 
    save(res_nmf, file = file.path(filepath, paste0(.x, ".rda")))
    write.table(res_nmf@consensus, file.path(filepath, paste0(.x, "_consensus.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(coef(res_nmf) %>% as.matrix, file.path(filepath, paste0(.x, "_coef.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  })
}


#' Plots consensus results from NMF analysis
#'
#' @import NMF dplyr purrr rlang cluster Biobase
#' @importFrom magrittr %>%
plotNMFCon <- function(id, rank, label_scheme_sub, scale_log2r, complete_cases, impute_na, 
                       filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
  
  ins <- list.files(path = filepath, pattern = "_rank\\d+\\.rda$") 
  if (impute_na) ins <- ins %>% .[grepl("_impNA", .)] else ins <- ins %>% .[!grepl("_impNA", .)]
  if (scale_log2r) ins <- ins %>% .[grepl("_NMF_Z", .)] else ins <- ins %>% .[grepl("_NMF_N", .)]

  if (is.null(rank)) {
    filelist <- ins
  } else {
    stopifnot(all(rank >= 2) & all(rank %% 1 == 0))
    
    filelist <- local({
      possible <- ins %>% 
        gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .) %>% 
        as.numeric() %>% 
        `names<-`(ins)
      
      r2 <- rank %>% .[. %in% possible]
      
      filelist <- possible %>% 
        .[. %in% r2] %>% 
        names(.)
    })
  } 

  if (id %in% c("pep_seq", "pep_seq_mod")) {
    custom_prefix <- purrr::map_chr(filelist, ~ {
      gsub("(.*_{0,1})Peptide_NMF.*", "\\1", .x)
    })
  } else if (id %in% c("prot_acc", "gene")) {
    custom_prefix <- purrr::map_chr(filelist, ~ {
      gsub("(.*_{0,1})Protein_NMF.*", "\\1", .x)
    })
  } else {
    stop("Unknown id = ", id, call. = FALSE)
  }

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  if (purrr::is_empty(filelist)) 
    stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
         " at rank = ", paste0(rank, collapse = ", "), call. = FALSE)

  purrr::walk2(filelist, custom_prefix, ~ {
    rank <- gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .x) %>% 
      as.numeric()
    
    out_nm <- paste0(.y, fn_prefix, "_rank", rank, ".", fn_suffix)

    load(file = file.path(file.path(filepath, .x)))
    
    D_matrix <- res_nmf@consensus
    if (complete_cases) D_matrix <- D_matrix %>% .[complete.cases(.), ]
    
    D_matrix <- data.frame(D_matrix, check.names = FALSE) %>% 
      filters_in_call(!!!filter_dots) %>% 
      arrangers_in_call(!!!arrange_dots) %>% 
      dplyr::select(which(names(.) %in% sample_ids)) %>% 
      tibble::rownames_to_column() %>% 
      dplyr::filter(rowname %in% sample_ids) %>% 
      tibble::column_to_rownames() 

    n_color <- 50
    xmin <- 0
    xmax <- ceiling(max(D_matrix))
    xmargin <- (xmax - xmin)/2
    color_breaks <- c(seq(xmin, xmargin, length = n_color/2)[1 : (n_color/2-1)],
                      seq(xmargin, xmax, length = n_color/2)[2 : (n_color/2)])
      
    if (is.null(dots$units)) {
      units <- "in"
    } else {
      units <- dots$units
    }
    
    if (is.null(dots$res)) {
      res <- 300
    } else {
      res <- dots$res
    }
    
    if (is.null(dots$width)) {
      width <- 1.35 * ncol(D_matrix)
    } else {
      width <- dots$width
    }
    
    if (width > 40 & units == "in") {
      warning("The plot width is reduced to 40")
      width <- 40
    }
    
    if (is.null(dots$height)) {
      height <- 1.35 * ncol(D_matrix)
    } else {
      height <- dots$height
    }
    
    if (height > 40 & units == "in") {
      warning("The plot height is reduced to 40")
      height <- 40
    }
    
    if (is.null(dots$color)) {
      mypalette <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(n_color)
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
    
    if (is.null(annot_cols)) annotation_col <- NA else
      annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = colnames(D_matrix))
    
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
    
    clus <- cluster::silhouette(res_nmf)
    attr(clus, "Ordered") <- NULL
    attr(clus, "call") <- NULL
    attr(clus, "class") <- NULL
    clus <- data.frame(clus, check.names = FALSE)
    clus <- clus %>% .[rownames(.) %in% label_scheme_sub$Sample_ID, ]
    
    annotation_col <- annotation_col %>% 
      tibble::rownames_to_column() %>% 
      dplyr::bind_cols(clus) %>% 
      dplyr::select(-c("neighbor", "sil_width")) %>% 
      dplyr::rename(silhouette = cluster) %>% 
      dplyr::mutate(silhouette = factor(silhouette)) %>% 
      tibble::column_to_rownames()
    
    p <- my_pheatmap(
      mat = D_matrix,
      filename = file.path(filepath, out_nm),
      annotation_col = annotation_col,
      annotation_row = NA, 
      color = mypalette,
      annotation_colors = annotation_colors,
      breaks = color_breaks,
      !!!dots
    )

  }, complete_cases = complete_cases)
}


#' Plots coef results from NMF analysis
#'
#' @import NMF dplyr purrr rlang cluster Biobase
#' @importFrom magrittr %>%
plotNMFCoef <- function(id, rank, label_scheme_sub, scale_log2r, complete_cases, impute_na, 
                        filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
  
  ins <- list.files(path = filepath, pattern = "_rank\\d+\\.rda$")
  if (impute_na) ins <- ins %>% .[grepl("_impNA", .)] else ins <- ins %>% .[!grepl("_impNA", .)]
  if (scale_log2r) ins <- ins %>% .[grepl("_NMF_Z", .)] else ins <- ins %>% .[grepl("_NMF_N", .)]

  if (is.null(rank)) {
    filelist <- ins
  } else {
    stopifnot(all(rank >= 2) & all(rank %% 1 == 0))
    
    filelist <- local({
      possible <- ins %>% 
        gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .) %>% 
        as.numeric() %>% 
        `names<-`(ins)
      
      r2 <- rank %>% .[. %in% possible]
      
      filelist <- possible %>% 
        .[. %in% r2] %>% 
        names(.)
    })
  }
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    custom_prefix <- purrr::map_chr(filelist, ~ {
      gsub("(.*_{0,1})Peptide_NMF.*", "\\1", .x)
    })
  } else if (id %in% c("prot_acc", "gene")) {
    custom_prefix <- purrr::map_chr(filelist, ~ {
      gsub("(.*_{0,1})Protein_NMF.*", "\\1", .x)
    })
  } else {
    stop("Unknown id = ", id, call. = FALSE)
  }

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  if (purrr::is_empty(filelist)) 
    stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
         " at rank = ", paste0(rank, collapse = ", "), call. = FALSE)

  purrr::walk2(filelist, custom_prefix, ~ {
    rank <- gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .x) %>% 
      as.numeric()
    
    out_nm <- paste0(.y, fn_prefix, "_rank", rank, ".", fn_suffix)
    
    src_path <- file.path(filepath, .x)
    load(file = file.path(src_path))
    
    D_matrix <- coef(res_nmf) 
    if (complete_cases) D_matrix <- D_matrix %>% .[complete.cases(.), ]
    rownames(D_matrix) <- seq_along(1:nrow(D_matrix))
    
    D_matrix <- data.frame(D_matrix, check.names = FALSE) %>% 
      filters_in_call(!!!filter_dots) %>% 
      arrangers_in_call(!!!arrange_dots) %>% 
      dplyr::select(which(names(.) %in% sample_ids))

    n_color <- 50
    xmin <- 0
    xmax <- max(D_matrix)
    xmargin <- xmax/2
    color_breaks <- c(seq(xmin, xmargin, length = n_color/2)[1 : (n_color/2-1)],
                      seq(xmargin, xmax, length = n_color/2)[2 : (n_color/2)])
    
    if (is.null(dots$width)) {
      width <- 1.35 * ncol(D_matrix)
    } else {
      width <- dots$width
    }
    
    if (is.null(dots$height)) {
      height <- 1.35 * ncol(D_matrix)
    } else {
      height <- dots$height
    }
    
    if (is.null(dots$color)) {
      mypalette <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(n_color)
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
    
    annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = colnames(D_matrix))
    
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
    
    clus <- cluster::silhouette(res_nmf)
    attr(clus, "Ordered") <- NULL
    attr(clus, "call") <- NULL
    attr(clus, "class") <- NULL
    clus <- data.frame(clus, check.names = FALSE)
    
    annotation_col <- annotation_col %>% 
      tibble::rownames_to_column() %>% 
      dplyr::bind_cols(clus) %>% 
      dplyr::select(-c("neighbor", "sil_width")) %>% 
      dplyr::rename(silhouette = cluster) %>% 
      dplyr::mutate(silhouette = factor(silhouette)) %>% 
      tibble::column_to_rownames()
    
    p <- my_pheatmap(
      mat = D_matrix,
      filename = file.path(filepath, out_nm),
      annotation_col = annotation_col,
      annotation_row = NA, 
      color = mypalette,
      annotation_colors = annotation_colors,
      breaks = color_breaks,
      !!!dots
    )

  })
}


#' Plots coef results from NMF analysis
#'
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMFmeta <- function(df, id, rank, label_scheme_sub, scale_log2r, complete_cases, impute_na, 
                        filepath, filename, anal_type, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  
  complete_cases <- to_complete_cases(complete_cases = complete_cases, impute_na = impute_na)
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)    
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]  
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    prepDM(id = !!id, scale_log2r = scale_log2r, 
           sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R
  
  ins <- list.files(path = filepath, pattern = "_rank\\d+\\.rda$")
  if (impute_na) ins <- ins %>% .[grepl("_impNA", .)] else ins <- ins %>% .[!grepl("_impNA", .)]
  if (scale_log2r) ins <- ins %>% .[grepl("_NMF_Z", .)] else ins <- ins %>% .[grepl("_NMF_N", .)]
  
  if (is.null(rank)) {
    filelist <- ins
  } else {
    stopifnot(all(rank >= 2) & all(rank %% 1 == 0))
    
    filelist <- local({
      possible <- ins %>% 
        gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .) %>% 
        as.numeric() %>% 
        `names<-`(ins)
      
      r2 <- rank %>% .[. %in% possible]
      
      filelist <- possible %>% 
        .[. %in% r2] %>% 
        names(.)
    })
  } 
  
  if (purrr::is_empty(filelist)) 
    stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
         " at rank = ", paste0(rank, collapse = ", "), call. = FALSE)

  if (id %in% c("pep_seq", "pep_seq_mod")) {
    custom_prefix <- purrr::map_chr(filelist, ~ {
      gsub("(.*_{0,1})Peptide_NMF.*", "\\1", .x)
    })
  } else if (id %in% c("prot_acc", "gene")) {
    custom_prefix <- purrr::map_chr(filelist, ~ {
      gsub("(.*_{0,1})Protein_NMF.*", "\\1", .x)
    })
  } else {
    stop("Unknown id = ", id, call. = FALSE)
  }
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)  
  
  purrr::walk2(filelist, custom_prefix, ~ {
    rank <- gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .x) %>% as.numeric()
    dir.create(file.path(filepath, rank), recursive = TRUE, showWarnings = FALSE)
    
    out_nm <- paste0(.y, fn_prefix, "_rank", rank, ".", fn_suffix)

    src_path <- file.path(filepath, .x)
    load(file = file.path(src_path))
    
    V_hat <- NMF::fitted(res_nmf)
    s <- NMF::extractFeatures(res_nmf)

    if (is.null(dots$xmin)) {
      xmin <- -1
    } else {
      xmin <- eval(dots$xmin, env = caller_env())
    }
    
    if (is.null(dots$xmax)) {
      xmax <- 1
    } else {
      xmax <- eval(dots$xmax, env = caller_env())
    }
    
    if (is.null(dots$xmargin)) {
      xmargin <- .1
    } else {
      xmargin <- eval(dots$xmargin, env = caller_env())
    }
    
    n_color <- 500
    color_breaks <- c(seq(xmin, xmargin, length = n_color/2)[1 : (n_color/2-1)],
                      seq(xmargin, xmax, length = n_color/2)[2 : (n_color/2)])
    
    if (is.null(dots$color)) {
      mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
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
    
    annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = colnames(df))
    
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
        
    for (i in seq_len(rank)) {
      df_sub <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ]
      
      nrow <- nrow(df_sub)
      
      if (nrow > 0) {
        fn_sub <- paste0(fn_prefix, "_", i, ".", fn_suffix)
        width <- ncol(df_sub) * 2 + 2
        
        if (nrow > 300) {
          height <- width * 1.5
          dots$show_rownames <- FALSE
        } else {
          cellheight <- 5
          height <- cellheight * nrow + 8
          fontsize_row <- 5
          
          dots$cellheight <- cellheight
          dots$fontsize_row <- fontsize_row
        }
        
        if (nrow <= 150) dots$show_rownames <- TRUE
        
        dots <- dots %>% 
          .[!names(.) %in% c("annot_cols", "annot_colnames", "annot_rows", 
                             "mat", "filename", "annotation_col", "annotation_row", 
                             "color", "annotation_colors", "breaks")]
        
        p <- my_pheatmap(
          mat = df_sub,
          filename = file.path(filepath, rank, fn_sub),
          annotation_col = annotation_col,
          annotation_row = NA, 
          color = mypalette,
          annotation_colors = annotation_colors,
          breaks = color_breaks,
          !!!dots
        )
        
        df_op <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ] %>%
          tibble::rownames_to_column(id)
        
        write.table(df_op, file.path(filepath, rank, paste0(fn_prefix, "_", i, ".txt")), 
                    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      }
    }

  })
}


#'NMF Classification
#'
#'\code{anal_pepNMF} performs the NMF classification of peptide \code{log2FC}.
#'The function is a wrapper of \code{\link[NMF]{nmf}}.
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}.
#'
#'@inheritParams anal_prnTrend
#'@inheritParams  prnEucDist
#'@inheritParams  prnHM
#'@inheritParams  info_anal
#'@inheritParams standPep
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is TRUE.
#'@param col_group Character string to a column key in \code{expt_smry.xlsx}.
#'  Samples corresponding to non-empty entries under \code{col_group} will be
#'  used for sample grouping in the indicated analysis. At the NULL default, the
#'  column key \code{Group} will be used. No data annotation by groups will be
#'  performed if the fields under the indicated group column is empty.
#'@param rank Numeric vector; the factorization rank(s) in
#'  \code{\link[NMF]{nmf}}. The default is c(4:8)
#'@param nrun Numeric; the number of runs in \code{\link[NMF]{nmf}}. The default
#'  is 50.
#'@param filepath Use system default.
#'@param filename A representative file name to outputs. By default, it will be
#'  determined automatically by the name of the current call.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data against the column keys in \code{Peptide.txt}, \code{Protein.txt} etc.;
#'  also see \code{\link{normPSM}}. \cr \cr \code{arrange_}: Logical
#'  expression(s) for the row ordering of data; also see \code{\link{prnHM}}.
#'  \cr \cr Additional arguments for \code{\link[NMF]{nmf}}.
#'@return NMF classification of \code{log2FC} data.
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'@example inst/extdata/examples/prnNMF_.R
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
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend
#'  analysis and visualization \cr \code{\link{anal_pepNMF}},
#'  \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'  \code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for NMF
#'  analysis and visualization \cr
#'
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@export
anal_pepNMF <- function (col_select = NULL, col_group = NULL, 
                         scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                         df = NULL, filepath = NULL, filename = NULL, 
                         rank = NULL, nrun = if (length(rank) > 1) 50 else 1, seed = NULL, ...) {
  on.exit(
    mget(names(formals()), current_env()) %>% 
      c(enexprs(...)) %>% 
      save_call(paste0("anal", "_pepNMF"))
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type"), ...)
  
  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_group <- rlang::enexpr(col_group)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = !!col_group, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, filepath = !!filepath, filename = !!filename, 
            anal_type = "NMF")(rank = rank, nrun = nrun, seed = seed, ...)
}


#'NMF Classification
#'
#'\code{anal_prnNMF} performs the NMF classification of protein \code{log2FC}.
#'The function is a wrapper of \code{\link[NMF]{nmf}}.
#'
#'@rdname anal_pepNMF
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'
#'@export
anal_prnNMF <- function (col_select = NULL, col_group = NULL, 
                         scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                         df = NULL, filepath = NULL, filename = NULL, 
                         rank = NULL, nrun = if (length(rank) > 1) 50 else 1, seed = NULL, ...) {
  on.exit(
    mget(names(formals()), current_env()) %>% 
          c(enexprs(...)) %>% 
          save_call(paste0("anal", "_prnNMF"))
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  col_group <- rlang::enexpr(col_group)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  dots <- rlang::enexprs(...)
  if (!is.null(dots$method)) dots$method <- rlang::as_string(dots$method)
  
  info_anal(id = !!id, col_select = !!col_select, col_group = !!col_group, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, filepath = !!filepath, filename = !!filename, 
            anal_type = "NMF")(rank = rank, nrun = nrun, seed = seed, !!!dots)
}


#'NMF plots
#'
#'\code{plot_pepNMFCon} plots the consensus heat maps from the NMF
#'classification of peptide \code{log2FC}.
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}.
#'
#'@param rank Numeric vector; the factorization rank(s) in
#'  \code{\link[NMF]{nmf}}. At the NULL default, all available ranks from the
#'  results of \code{\link{anal_pepNMF}} or \code{\link{anal_pepNMF}} will be
#'  used.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data against the column keys in \code{_NMF[...]_consensus.txt} for consensus
#'  plots or \code{_NMF[...]_coef.txt} for coefficient plots; also see
#'  \code{\link{normPSM}}. \cr \cr \code{arrange_}: Logical expression(s) for
#'  the row ordering of data; also see \code{\link{prnHM}}. \cr \cr Additional
#'  arguments for \code{\link[pheatmap]{pheatmap}}
#'@inheritParams prnHist
#'@inheritParams plot_prnTrend
#'@inheritParams  prnEucDist
#'@return Concensus heat maps from NMF classification.
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'@example inst/extdata/examples/prnNMF_.R
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
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend
#'  analysis and visualization \cr \code{\link{anal_pepNMF}},
#'  \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'  \code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for NMF
#'  analysis and visualization \cr
#'
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@export
plot_pepNMFCon <- function (col_select = NULL, 
                            scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                            filename = NULL, 
                            annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "col_group", "df", "filepath"), ...)

  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)  
  
  reload_expts()

  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, filepath = NULL, filename = !!filename, 
            anal_type = "NMF_con")(rank = rank, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'NMF plots
#'
#'\code{plot_prnNMFCon} plots the consensus heat maps from the NMF
#'classification of protein \code{log2FC}.
#'
#'@rdname plot_pepNMFCon
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'@export
plot_prnNMFCon <- function (col_select = NULL, 
                            scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                            filename = NULL, 
                            annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "col_group", "df", "filepath"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)  
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, filepath = NULL, filename = !!filename, 
            anal_type = "NMF_con")(rank = rank, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'NMF plots
#'
#'\code{plot_pepNMFCoef} plots the coefficient heat maps from the NMF
#'classification of peptide \code{log2FC}.
#'
#'@rdname plot_pepNMFCon
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'@export
plot_pepNMFCoef <- function (col_select = NULL, 
                             scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE, 
                             filename = NULL, 
                             annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "col_group", "df", "filepath"), ...)
  
  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)  
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, filepath = NULL, filename = !!filename, 
            anal_type = "NMF_coef")(rank = rank, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'NMF plots
#'
#'\code{plot_prnNMFCoef} plots the coefficient heat maps from the NMF
#'classification of protein \code{log2FC}.
#'
#'@rdname plot_pepNMFCon
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'@export
plot_prnNMFCoef <- function (col_select = NULL, 
                             scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                             filename = NULL, 
                             annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "col_group", "df", "filepath"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)  
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, filepath = NULL, filename = !!filename, 
            anal_type = "NMF_coef")(rank = rank, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'Heat maps of metagenes from NMF
#'
#'\code{plot_metaNMF} is a wrapper of \code{\link[pheatmap]{pheatmap}} for the
#'visualization of the metagene heat maps from NMF
#'
#'@inheritParams anal_pepNMF
#'@inheritParams  prnEucDist
#'@param rank Numeric vector; the factorization rank(s) in
#'  \code{\link[NMF]{nmf}}. At the NULL default, all available ranks from the
#'  results of \code{\link{anal_pepNMF}} or \code{\link{anal_pepNMF}} will be
#'  used.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data against the column keys in \code{Peptide.txt}, \code{Protein.txt} etc.;
#'  also see \code{\link{normPSM}}. \cr \cr \code{arrange_}: Logical
#'  expression(s) for the row ordering of data; also see \code{\link{prnHM}}.
#'  \cr \cr Additional arguments for \code{\link[pheatmap]{pheatmap}}.
#'@import purrr
#'@example inst/extdata/examples/prnNMF_.R
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
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend
#'  analysis and visualization \cr \code{\link{anal_pepNMF}},
#'  \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'  \code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for NMF
#'  analysis and visualization \cr
#'
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@export
plot_metaNMF <- function (col_select = NULL, 
                          scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                          df = NULL, filepath = NULL, filename = NULL, 
                          rank = NULL, annot_cols = NULL, annot_colnames = NULL, ...) {
  check_dots(c("id", "anal_type", "col_group"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)  
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, filepath = !!filepath, filename = !!filename, 
            anal_type = "NMF_meta")(rank = rank, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}

