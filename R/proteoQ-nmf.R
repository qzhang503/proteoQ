#' NMF analysis
#'
#' @import dplyr rlang Biobase
#' @importFrom magrittr %>%
#' @importFrom NMF nmf
nmfTest <- function(df, id, r, nrun, col_group, label_scheme_sub, filepath, filename,
                    complete_cases, ...) {
  
  stopifnot(nrow(label_scheme_sub) > 0)
  
  sample_ids <- label_scheme_sub$Sample_ID

  col_group <- rlang::enexpr(col_group) # optional phenotypic information
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  fn_prx <- gsub("\\..*$", "", filename)
  fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
  
  if (complete_cases) df <- df[complete.cases(df), ]
  
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

  # featureNames(exampleSet)[1:5]
  # phenoData(exampleSet)
  # mat <- Biobase::exprs(exampleSet)
  
  res_nmf <- NMF::nmf(exampleSet, r, nrun = nrun)

  save(res_nmf, file = file.path(filepath, paste0(fn_prx, ".rda")))
  write.csv(res_nmf@consensus, file.path(filepath, paste0(fn_prx, "_consensus.csv")), row.names = FALSE)
  write.csv(coef(res_nmf) %>% as.matrix, file.path(filepath, paste0(fn_prx, "_coef.csv")), row.names = FALSE)
}



#' Plots consensus results from NMF analysis
#'
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMFCon <- function(id, r, label_scheme_sub, filepath, in_nm, out_nm, 
                       # annot_cols = NULL, annot_colnames = NULL, 
                       ...) {

  stopifnot(nrow(label_scheme_sub) > 0)
  
  sample_ids <- label_scheme_sub$Sample_ID
  
  src_path <- file.path(filepath, in_nm)
  load(file = file.path(src_path))
  
  fn_prx <- gsub("\\..*$", "", in_nm)

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

  D_matrix <- res_nmf@consensus
    
  n_color <- 500
  xmin <- 0
  xmax <- ceiling(max(D_matrix))
  x_margin <- xmax/10
  color_breaks <- c(seq(xmin, x_margin, length = n_color/2)[1 : (n_color/2-1)],
                    seq(x_margin, xmax, length = n_color/2)[2 : (n_color/2)])
  
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
    
  png(file.path(filepath, paste0(fn_prx, "_consensus.png")), 
      width = width, height = height, units="in", res = 300)
    consensusmap(res_nmf, annCol = annotation_col, 
                 annColor=list(Type = 'Spectral', basis = 'Set3',consensus = 'YlOrRd:50'),
                 tracks = c("basis:"), main = '', sub = '')
  dev.off()
  
  # my_pheatmap(
  #   mat = D_matrix,
  #   filename = file.path(filepath, out_nm), 
  #   annotation_col = annotation_col,
  #   color = mypalette,
  #   annotation_colors = annotation_colors,
  #   breaks = color_breaks,
  #   !!!dots
  # )
}


#' Plots coef results from NMF analysis
#'
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMFCoef <- function(id, r, label_scheme_sub, filepath, in_nm, out_nm, 
                       # annot_cols = NULL, annot_colnames = NULL, 
                       ...) {
  
  stopifnot(nrow(label_scheme_sub) > 0)
	
  sample_ids <- label_scheme_sub$Sample_ID
	
  load(file = file.path(filepath, in_nm))
	
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  D_matrix <- coef(res_nmf) %>% as.matrix

  n_color <- 500
  xmin <- 0
  xmax <- ceiling(max(D_matrix))
  x_margin <- xmax/10
  color_breaks <- c(seq(xmin, x_margin, length = n_color/2)[1 : (n_color/2-1)],
                    seq(x_margin, xmax, length = n_color/2)[2 : (n_color/2)])
  
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
  
  png(file.path(filepath, paste0(gsub("\\..*$", "", in_nm), "_coef.png")), 
      width = width, height = height, units="in", res = 300)
    coefmap(res_nmf, annCol = annotation_col, 
            annColor = list(Type = 'Spectral', basis = 'Set3', consensus = 'YlOrRd:50'),
            tracks = c("basis:"))
  dev.off()
}


#' Plots coef results from NMF analysis
#'
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMFmeta <- function(df, id, r, label_scheme_sub, filepath, in_nm, out_nm, 
                        # complete_cases, 
                        ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
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
  
  if (is.null(dots$x_margin)) {
    x_margin <- .1
  } else {
    x_margin <- eval(dots$x_margin, env = caller_env())
  }
  
  n_color <- 500
  color_breaks <- c(seq(xmin, x_margin, length = n_color/2)[1 : (n_color/2-1)],
                    seq(x_margin, xmax, length = n_color/2)[2 : (n_color/2)])
  
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
  
  fn_prx <- gsub("\\..*$", "", in_nm)
  fn_suffix <- gsub(".*\\.(.*)$", "\\1", in_nm)
  
  # if (complete_cases) df <- df[complete.cases(df), ]

  # metagene-specific features
  load(file = file.path(filepath, in_nm))
  V_hat <- NMF::fitted(res_nmf)
  s <- NMF::extractFeatures(res_nmf)
  # s <- featureScore(res_nmf)
  # w <- basis(res_nmf) # get matrix W
  # h <- coef(res_nmf) # get matrix H
  
  for (i in seq_len(r)) {
    df_sub <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ]
    
    nrow <- nrow(df_sub)
    
    if (nrow > 0) {
      fn_sub <- paste0(fn_prx, "_metagene", i, ".png")
      width <- ncol(df_sub) * 2 + 2
      
      if (nrow > 300) {
        height <- width * 1.5
        dots$show_rownames = FALSE
      } else {
        cellheight <- 5
        height <- cellheight * nrow + 8
        fontsize_row <- 5
        
        dots$cellheight <- cellheight
        dots$fontsize_row <- fontsize_row
      }
      
      p <- my_pheatmap(
        mat = df_sub,
        filename = file.path(filepath, fn_sub),
        annotation_col = annotation_col,
        color = mypalette,
        annotation_colors = annotation_colors,
        breaks = color_breaks,
        !!!dots
      )
      
      df_op <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ] %>%
        tibble::rownames_to_column(id)
      
      write.csv(df_op, file.path(filepath, paste0(fn_prx, "_metagene", i, ".csv")),
                row.names = FALSE)
    }
  }
  
}


#'NMF analysis
#'
#'Analyzes and visualizes the NMF clustering of peptide or protein \code{log2FC}
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}
#'
#'@inheritParams  proteoEucDist
#'@inheritParams  proteoHM
#'@inheritParams  info_anal
#'@param col_group Character string to a column key in \code{expt_smry.xlsx}.
#'  Samples corresponding to non-empty entries under \code{col_group} will be
#'  used for sample grouping in the indicated analysis. At the NULL default, the
#'  column key \code{Group} will be used.
#'@param r Numeric; the factorization rank (\code{\link[NMF]{nmf}}).
#'@param nrun Numeric; the number of runs to perform (\code{\link[NMF]{nmf}}).
#'@param task Character string; a switch for different tasks in a functional
#'  factory.
#'@param ... In \code{anal_} functions: additional arguments inherited from
#'  \code{\link[NMF]{nmf}}; in \code{plot_} functions: additional arguments
#'  inherited from \code{proteoEucDist}.
#'@return NMF classification and visualization of \code{log2FC}.
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoNMF <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                       col_select = NULL, col_group = NULL, scale_log2r = TRUE, impute_na = TRUE, 
                       complete_cases = FALSE, df = NULL, filepath = NULL, filename = NULL, 
                       anal_type = "NMF", task = "anal", r = 4, nrun = 200, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

  id <- rlang::enexpr(id)
	if (length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
	
	stopifnot(rlang::is_logical(scale_log2r))
	stopifnot(rlang::is_logical(impute_na))
	stopifnot(rlang::is_logical(complete_cases))
	stopifnot(rlang::as_string(rlang::enexpr(anal_type)) == "NMF")
	stopifnot(r >= 2 & r %% 1 == 0)
	stopifnot(nrun >= 1 & nrun %% 1 == 0)

	task <- rlang::enexpr(task)	
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()

	if (!impute_na) complete_cases <- TRUE

	info_anal(id = !!id, col_select = !!col_select, col_group = !!col_group, scale_log2r = scale_log2r, 
	          impute_na = impute_na, df = !!df, filepath = !!filepath, filename = !!filename, 
	          anal_type = "NMF")(r = r, nrun = nrun, complete_cases = complete_cases, 
	                             task = !!task, ...)
}


#'Peptide NMF Classification
#'
#'\code{anal_pepNMF} is a wrapper function of \code{\link{proteoNMF}} for the
#'NMF analysis of peptide data
#'
#'@rdname proteoNMF
#'@examples
#' library(NMF)
#' anal_pepNMF(
#'   scale_log2r = TRUE,
#'   col_group = Group, # optional a priori knowledge of sample groups
#'   r = 6,
#'   nrun = 200
#' )
#'
#'@export
anal_pepNMF <- function (...) {
	proteoNMF(id = pep_seq, task = anal, ...)
}


#'Protein NMF Classification
#'
#'\code{anal_prnNMF} is a wrapper function of \code{\link{proteoNMF}} for the
#'NMF analysis of protein data
#'
#'@rdname proteoNMF
#' @examples
#' library(NMF)
#' anal_prnNMF(
#'   scale_log2r = TRUE,
#'   col_group = Group, # optional a priori knowledge of sample groups
#'   r = 6,
#'   nrun = 200
#' )
#'
#'@export
anal_prnNMF <- function (...) {
	proteoNMF(id = gene, task = anal, ...)
}


#'NMF consensus
#'
#'\code{plot_prnNMFCon} is a wrapper function of \code{\link{proteoNMF}} for the
#'visualization of the consensus heat map of protein data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#'plot_prnNMFCon(
#'   r = 6, 
#'   annot_cols = c("Color", "Alpha", "Shape"), 
#'   annot_colnames = c("Lab", "Batch", "WHIM"), 
#'   width = 10, 
#'   height = 10
#')
#'
#'@export
plot_prnNMFCon <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  proteoNMF(id = gene, task = plotcon, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'NMF coefficients
#'
#'\code{plot_prnNMFCoef} is a wrapper function of \code{\link{proteoNMF}} for the
#'visualization of the coefficient heat map of protein data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#'plot_prnNMFCoef(
#'   r = 6, 
#'   annot_cols = c("Color", "Alpha", "Shape"), 
#'   annot_colnames = c("Lab", "Batch", "WHIM"), 
#'   width = 10, 
#'   height = 10
#')
#'
#'@export
plot_prnNMFCoef <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)

  proteoNMF(id = gene, task = plotcoef, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'NMF coefficients
#'
#'\code{plot_metaNMF} is a wrapper function of \code{\link{proteoNMF}} for the
#'visualization of the metagene heat maps of protein data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#'plot_metaNMF(
#'   r = 6, 
#'   annot_cols = c("Color", "Alpha", "Shape"), 
#'   annot_colnames = c("Lab", "Batch", "WHIM"), 
#'   
#'   fontsize = 8, 
#'   fontsize_col = 5
#')
#'
#'@export
plot_metaNMF <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)

  proteoNMF(id = gene, task = plotmeta, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


