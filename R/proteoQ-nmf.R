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
plotNMFmeta <- function(df, id, r, label_scheme_sub, filepath, in_nm, out_nm, ...) {

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
  
  # metagene-specific features
  load(file = file.path(filepath, in_nm))
  V_hat <- fitted(res_nmf)
  s <- extractFeatures(res_nmf)
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
#'Analyzes and visualizes the NMF clustering of peptide or protein
#'\code{log2-ratios}
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at \code{id = pep_seq} or \code{pep_seq_mod}, and protein data at
#'  \code{id = prot_acc} or \code{gene}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param impute_na Logical; if TRUE, imputes missing values.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used for visualization.
#'@param r Numeric; the factorization rank (\code{\link[NMF]{nmf}}).
#'@param nrun Numeric; the number of runs to perform (\code{\link[NMF]{nmf}}).
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the names of the current functional \code{call}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... Additional parameters for plotting: \cr \code{xmin}, the minimum x; \cr
#'  \code{xmax}, the maximum x; \cr \code{x_margin}, the margin in heat
#'  scales;\cr \code{annot_cols}, the column keys in \code{expt_smry.xlsx} for
#'  use in the color coding of samples;\cr \code{width_}, the width of plot; \cr
#'  \code{height_}, the height of plot; \cr additional arguments inherited from
#'  \code{\link[pheatmap]{pheatmap}}.
#'@return NMF plots.
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoNMF <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                       col_select = NULL, scale_log2r = TRUE, impute_na = TRUE, complete_cases = FALSE, 
                       df = NULL, filepath = NULL, filename = NULL, 
                       anal_type = "NMF", task = "anal", r = 4, nrun = 200, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

  id <- rlang::enexpr(id)
	if(length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
	
	# anal_type <- rlang::enexpr(anal_type)
	task <- rlang::enexpr(task)	

	col_select <- rlang::enexpr(col_select)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()

	if(!impute_na) complete_cases <- TRUE

	info_anal(id = !!id, col_select = !!col_select, scale_log2r = scale_log2r, impute_na = impute_na,
					df = !!df, filepath = !!filepath, filename = !!filename,
					anal_type = "NMF")(r = r, nrun = nrun, complete_cases = complete_cases,
					task = !!task, ...)
}


#'Peptide NMF Classification
#'
#'NMF Classification of peptide \code{log2-ratios}
#'
#' @examples
#' pepNMF(
#'   scale_log2r = TRUE,
#'   xmin = -1,
#'   xmax = 1,
#'   x_margin = 0.1,
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   r = 6,
#'   nrun = 200
#' )
#'
#'@seealso \code{\link{proteoNMF}} for parameters
#'@export
anal_pepNMF <- function (...) {
	proteoNMF(id = pep_seq, task = anal, ...)
}


#'Protein NMF Classification
#'
#'NMF Classification of protein \code{log2-ratios}
#'
#' @examples
#' prnNMF(
#'   scale_log2r = TRUE,
#'   xmin = -1,
#'   xmax = 1,
#'   x_margin = 0.1,
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   r = 6,
#'   nrun = 200
#' )
#'
#'
#'@seealso \code{\link{proteoNMF}} for parameters
#'@export
anal_prnNMF <- function (...) {
	proteoNMF(id = gene, task = anal, ...)
}


#'NMF consensus
#'
#'Visualizes the consensus heat map
#'
#' @examples
#' plot_prnNMFCon(
#'   scale_log2r = TRUE,
#'   col_order = Order,
#'   n_clust = 6
#' )
#'
#'@seealso \code{\link{proteoTrend}} for parameters
#'@export
plot_prnNMFCon <- function (...) {
  proteoNMF(id = gene, task = plotcon, ...)
}


#'NMF coefficients
#'
#'Visualizes the coefficient heat map
#'
#' @examples
#' plot_prnNMFCoef(
#'   scale_log2r = TRUE,
#'   col_order = Order,
#'   n_clust = 6
#' )
#'
#'@seealso \code{\link{proteoTrend}} for parameters
#'@export
plot_prnNMFCoef <- function (...) {
  proteoNMF(id = gene, task = plotcoef, ...)
}


#'NMF coefficients
#'
#'Visualizes the coefficient heat map
#'
#' @examples
#' plot_metaNMF(
#'   scale_log2r = TRUE,
#'   col_order = Order,
#'   n_clust = 6
#' )
#'
#'@seealso \code{\link{proteoTrend}} for parameters
#'@export
plot_metaNMF <- function (...) {
  proteoNMF(id = gene, task = plotmeta, ...)
}


