#' NMF analysis
#'
#' @import dplyr purrr rlang Biobase
#' @importFrom magrittr %>%
#' @importFrom NMF nmf
nmfTest <- function(df, id, r, nrun, col_group, label_scheme_sub, anal_type, scale_log2r, 
                    filepath, filename, complete_cases, ...) {
                    
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  df <- df %>% filters_in_call(!!!filter_dots)

  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R  
  
  col_group <- rlang::enexpr(col_group) # optional phenotypic information

	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% unique()
	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	stopifnot(length(fn_suffix) == 1)
	
  if (is.null(r)) {
    r <- c(4:8)
    nrun <- 5
    fn_prefix <- paste0(fn_prefix, r)
  } else {
    stopifnot(all(r >= 2) & all(r %% 1 == 0))
    stopifnot(all(nrun >= 1) & all(nrun %% 1 == 0))
  }
  
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

  purrr::walk(fn_prefix, ~ {
    r <- gsub(".*_r(\\d+)$", "\\1", .x) %>% as.numeric()
    
    res_nmf <- NMF::nmf(exampleSet, r, nrun = nrun)
    save(res_nmf, file = file.path(filepath, paste0(.x, ".rda")))
    write.csv(res_nmf@consensus, 
              file.path(filepath, paste0(.x, "_consensus.csv")), row.names = FALSE)
    write.csv(coef(res_nmf) %>% as.matrix, 
              file.path(filepath, paste0(.x, "_coef.csv")), row.names = FALSE)
  })
}


#' Plots consensus results from NMF analysis
#'
#' @import NMF dplyr purrr rlang Biobase
#' @importFrom magrittr %>%
plotNMFCon <- function(id, r, label_scheme_sub, filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  ins <- list.files(path = filepath, pattern = "_r\\d+\\.rda$")
 
  if (is.null(r)) {
    filelist <- ins
  } else {
    stopifnot(all(r >= 2) & all(r %% 1 == 0))
    
    ins_prefix <- gsub("\\.[^.]*$", "", ins)
    
    possible_prefix <- ins_prefix %>% 
      gsub("(.*)_r\\d+$", "\\1", .) %>% 
      unique() %>% 
      paste0("_r", r)
    
    ok_prefix <- ins_prefix %>% 
      .[. %in% possible_prefix]
    rm(ins_prefix, possible_prefix)
    
    filelist <- purrr::map(ok_prefix, ~ list.files(path = filepath, pattern = paste0(.x, "\\.rda$"))) %>% 
      unlist()
  } 
  
  if (purrr::is_empty(filelist)) 
    stop("Missing NMF results under ", filepath, 
         "\nCheck the setting in `scale_log2r` for a probable mismatch.", 
         call. = FALSE)
  
  # both filename extension and prefix ignored
  # only png for now
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  purrr::walk(filelist, ~ {
    out_nm <- paste0(gsub("\\.rda$", "", .x), "_consensus.png")
    
    load(file = file.path(file.path(filepath, .x)))
    
    D_matrix <- res_nmf@consensus
      
    n_color <- 500
    xmin <- 0
    xmax <- ceiling(max(D_matrix))
    xmargin <- xmax/10
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
    
    png(file.path(filepath, out_nm), width = width, height = height, units = units, res = res)
    consensusmap(res_nmf, annCol = annotation_col, 
                 annColor = list(Type = 'Spectral', basis = 'Set3',consensus = 'YlOrRd:50'),
                 tracks = c("basis:"), main = '', sub = '')
    dev.off()
    
  })
}


#' Plots coef results from NMF analysis
#'
#' @import NMF dplyr purrr rlang Biobase
#' @importFrom magrittr %>%
plotNMFCoef <- function(id, r, label_scheme_sub, filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  ins <- list.files(path = filepath, pattern = "_r\\d+\\.rda$")
  
  if (is.null(r)) {
    filelist <- ins
  } else {
    stopifnot(all(r >= 2) & all(r %% 1 == 0))
    
    ins_prefix <- gsub("\\.[^.]*$", "", ins)
    
    possible_prefix <- ins_prefix %>% 
      gsub("(.*)_r\\d+$", "\\1", .) %>% 
      unique() %>% 
      paste0("_r", r)
    
    ok_prefix <- ins_prefix %>% 
      .[. %in% possible_prefix]
    rm(ins_prefix, possible_prefix)
    
    filelist <- purrr::map(ok_prefix, ~ list.files(path = filepath, pattern = paste0(.x, "\\.rda$"))) %>%
      unlist()
  } 

  if (purrr::is_empty(filelist)) 
    stop("Missing NMF results under ", filepath, 
         "\nCheck the setting in `scale_log2r` for a probable mismatch.", 
         call. = FALSE)
  
  # both filename extension and prefix will be ignored
  # only png for now
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  purrr::walk(filelist, ~ {
    # out_nm <- paste0(gsub("\\.rda$", "", .x), "_coef.", fn_suffix)
    out_nm <- paste0(gsub("\\.rda$", "", .x), "_coef.png")
    src_path <- file.path(filepath, .x)
    load(file = file.path(src_path))
    
    D_matrix <- coef(res_nmf) %>% as.matrix
    
    n_color <- 500
    xmin <- 0
    xmax <- ceiling(max(D_matrix))
    xmargin <- xmax/10
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
    
    png(file.path(filepath, out_nm), width = width, height = height, units = "in", res = 300)
    coefmap(res_nmf, annCol = annotation_col, 
            annColor = list(Type = 'Spectral', basis = 'Set3', consensus = 'YlOrRd:50'),
            tracks = c("basis:"))
    dev.off()
  })
}


#' Plots coef results from NMF analysis
#'
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMFmeta <- function(df, id, r, label_scheme_sub, anal_type, scale_log2r, 
                        filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R  
  
  ins <- list.files(path = filepath, pattern = "_r\\d+\\.rda$")
  
  if (is.null(r)) {
    filelist <- ins
  } else {
    stopifnot(all(r >= 2) & all(r %% 1 == 0))
    
    ins_prefix <- gsub("\\.[^.]*$", "", ins)
    
    possible_prefix <- ins_prefix %>% 
      gsub("(.*)_r\\d+$", "\\1", .) %>% 
      unique() %>% 
      paste0("_r", r)
    
    ok_prefix <- ins_prefix %>% 
      .[. %in% possible_prefix]
    rm(ins, ins_prefix, possible_prefix)
    
    filelist <- purrr::map(ok_prefix, ~ list.files(path = filepath, pattern = paste0(.x, "\\.rda$"))) %>%
      unlist()
  } 
  
  if (purrr::is_empty(filelist)) 
    stop("Missing NMF results under ", filepath, 
         "\nCheck the setting in `scale_log2r` for a probable mismatch.", 
         call. = FALSE)

  # both filename extension and prefix will be ignored
  # only png for now
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  # fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  purrr::walk(filelist, ~ {
    fn_prefix <- gsub("\\.rda$", "", .x)
    r <- gsub(".*_r(\\d+)\\.rda$", "\\1", .x) %>% as.numeric()
    dir.create(file.path(filepath, r), recursive = TRUE, showWarnings = FALSE)

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
        
    for (i in seq_len(r)) {
      df_sub <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ]
      
      nrow <- nrow(df_sub)
      
      if (nrow > 0) {
        fn_sub <- paste0(fn_prefix, "_metagene", i, ".", fn_suffix)
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
        
        # probably no need to handle `cluster_rows` and `cluster_cols` for metagenes
        
        p <- my_pheatmap(
          mat = df_sub,
          filename = file.path(filepath, r, fn_sub),
          annotation_col = annotation_col,
          annotation_row = NA, 
          color = mypalette,
          annotation_colors = annotation_colors,
          breaks = color_breaks,
          !!!dots
        )
        
        df_op <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ] %>%
          tibble::rownames_to_column(id)
        
        write.csv(df_op, file.path(filepath, r, paste0(fn_prefix, "_metagene", i, ".csv")),
                  row.names = FALSE)
      }
    }

  })
}


#'NMF analysis
#'
#'\code{proteoNMF} analyzes and visualizes the NMF clustering of peptide or
#'protein \code{log2FC}. Users should avoid calling the method directly, but
#'instead use the following wrappers.
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
#'  column key \code{Group} will be used. No data annotation by groups will be
#'  performed if the fields under the indicated group column is empty.
#'@param r Numeric vector; the factorization rank(s) in \code{\link[NMF]{nmf}}.
#'  The default is c(4:8)
#'@param nrun Numeric; the number of runs in \code{\link[NMF]{nmf}}. The default
#'  is 200.
#'@param task Character string; a signature for task dispatching in a function
#'  factory. The value will be determined automatically.
#'@param filepath Use system default.
#'@param filename Use system default.
#'@param ... In \code{anal_} functions: additional arguments are for
#'  \code{\link[NMF]{nmf}}; in \code{plot_} functions: \code{width},
#'  \code{height}; in \code{plot_metaNMF} functions: additional arguments are
#'  for \code{pheatmap}.
#'@return NMF classification and visualization of \code{log2FC}.
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@example inst/extdata/examples/prnNMF_.R
#'
#'@export
proteoNMF <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                       col_select = NULL, col_group = NULL, scale_log2r = TRUE, impute_na = TRUE, 
                       complete_cases = FALSE, df = NULL, filepath = NULL, filename = NULL, 
                       task = "anal", r = NULL, nrun = 200, ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  id <- rlang::enexpr(id)
	if (length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
	
	stopifnot(rlang::is_logical(scale_log2r))
	stopifnot(rlang::is_logical(impute_na))
	stopifnot(rlang::is_logical(complete_cases))

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
#'\code{anal_pepNMF} is a wrapper of \code{\link{proteoNMF}} for the NMF
#'analysis of peptide data
#'
#'@rdname proteoNMF
#'@examples
#'
#' # ===================================
#' # NMF
#' # ===================================
#' scale_log2r <- TRUE
#' 
#' library(NMF)
#'
#' # peptide NMF at two different r(ank)s
#' anal_pepNMF(
#'   scale_log2r = TRUE,
#'   col_group = Group, # optional a priori knowledge of sample groups
#'   r = c(6, 8),
#'   nrun = 200,
#'   filter_by_npsm = exprs(pep_n_psm >= 2),
#' )
#'
#'@import purrr
#'@export
anal_pepNMF <- function (...) {
  err_msg <- "Don't call the function with arguments `id` and/or `task`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id", "task"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_pepid()
  
  quietly_log <- purrr::quietly(proteoNMF)(id = !!id, task = anal, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\NMF\\log","anal_pepNMF_log.csv"), append = TRUE)  
}


#'Protein NMF Classification
#'
#'\code{anal_prnNMF} is a wrapper of \code{\link{proteoNMF}} for the NMF
#'analysis of protein data
#'
#'@rdname proteoNMF
#' @examples
#' # protein NMF over a range of ranks
#' library(NMF)
#' 
#' anal_prnNMF(
#'   impute_na = FALSE,
#'   scale_log2r = TRUE,
#'   col_group = Group,
#'   r = c(5:8),
#'   nrun = 200, 
#'   filter_by_npep = exprs(prot_n_pep >= 2),
#' )
#'
#'@import purrr
#'@export
anal_prnNMF <- function (...) {
  err_msg <- "Don't call the function with arguments `id` and/or `task`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id", "task"))) stop(err_msg)

  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_protid()

  quietly_log <- purrr::quietly(proteoNMF)(id = !!id, task = anal, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\NMF\\log\\anal_prnNMF_log.csv"), append = TRUE)  
}


#'NMF consensus
#'
#'\code{plot_pepNMFCon} is a wrapper of \code{\link{proteoNMF}} for the
#'visualization of the consensus heat map of peptide data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#'
#' # peptide consensus heat maps at specific ranks
#' plot_pepNMFCon(
#'   r = c(5, 6),
#'   annot_cols = c("Color", "Alpha", "Shape"),
#'   annot_colnames = c("Lab", "Batch", "WHIM"),
#'   width = 10,
#'   height = 10,
#' )
#'
#' # peptide consensus heat maps at all available ranks
#' plot_pepNMFCon(
#'   impute_na = FALSE,
#'   annot_cols = c("Color", "Alpha", "Shape"),
#'   annot_colnames = c("Lab", "Batch", "WHIM"),
#'   width = 10,
#'   height = 10,
#' )
#'
#'@import purrr
#'@export
plot_pepNMFCon <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  err_msg <- "Don't call the function with arguments `annot_cols` and/or `annot_colnames`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)
  
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_pepid()
  
  quietly_log <- purrr::quietly(proteoNMF)(id = !!id, task = plotcon, 
                                           annot_cols = !!annot_cols, 
                                           annot_colnames = !!annot_colnames, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\NMF\\log","plot_pepNMFCon_log.csv"), append = TRUE)
}


#'NMF consensus
#'
#'\code{plot_prnNMFCon} is a wrapper of \code{\link{proteoNMF}} for the
#'visualization of the consensus heat map of protein data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#'
#' # protein consensus heat maps at specific ranks
#' plot_prnNMFCon(
#'   r = c(7:8),
#'   annot_cols = c("Color", "Alpha", "Shape"),
#'   annot_colnames = c("Lab", "Batch", "WHIM"),
#'   width = 10,
#'   height = 10,
#' )
#'
#' # protein consensus heat maps at all available ranks
#' plot_prnNMFCon(
#'   impute_na = FALSE,
#'   annot_cols = c("Color", "Alpha", "Shape"),
#'   annot_colnames = c("Lab", "Batch", "WHIM"),
#'   width = 10,
#'   height = 10
#' )
#'
#'@import purrr
#'@export
plot_prnNMFCon <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  err_msg <- "Don't call the function with arguments `annot_cols` and/or `annot_colnames`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)

  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_protid()
  
  quietly_log <- purrr::quietly(proteoNMF)(id = !!id, task = plotcon, 
                                           annot_cols = !!annot_cols, 
                                           annot_colnames = !!annot_colnames, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\NMF\\log\\plot_prnNMFCon_log.csv"), append = TRUE)
}


#'NMF coefficients
#'
#'\code{plot_pepNMFCoef} is a wrapper of \code{\link{proteoNMF}} for the
#'visualization of the coefficient heat map of peptide data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#'
#' # peptide coefficient heat maps at all ranks
#' plot_pepNMFCoef(
#'   annot_cols = c("Color", "Alpha", "Shape"),
#'   annot_colnames = c("Lab", "Batch", "WHIM"),
#'   width = 10,
#'   height = 10
#' )
#'
#'@import purrr
#'@export
plot_pepNMFCoef <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  err_msg <- "Don't call the function with arguments `annot_cols` and/or `annot_colnames`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)
  
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_pepid()
  
  quietly_log <- purrr::quietly(proteoNMF)(id = !!id, task = plotcoef, 
                                           annot_cols = !!annot_cols, 
                                           annot_colnames = !!annot_colnames, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\NMF\\log","plot_pepNMFCoef_log.csv"), append = TRUE)
}


#'NMF coefficients
#'
#'\code{plot_prnNMFCoef} is a wrapper of \code{\link{proteoNMF}} for the
#'visualization of the coefficient heat map of protein data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#'
#' # protein coefficient heat maps at all ranks
#' plot_prnNMFCoef(
#'   annot_cols = c("Color", "Alpha", "Shape"),
#'   annot_colnames = c("Lab", "Batch", "WHIM"),
#'   width = 10,
#'   height = 10
#' )
#'
#'@import purrr
#'@export
plot_prnNMFCoef <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  err_msg <- "Don't call the function with arguments `annot_cols` and/or `annot_colnames`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)

  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_protid()
  
  quietly_log <- purrr::quietly(proteoNMF)(id = !!id, task = plotcoef, 
                                           annot_cols = !!annot_cols, 
                                           annot_colnames = !!annot_colnames, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\NMF\\log\\plot_prnNMFCoef_log.csv"), append = TRUE)
}


#'NMF coefficients
#'
#'\code{plot_metaNMF} is a wrapper of \code{\link{proteoNMF}} for the
#'visualization of the metagene heat maps of protein data
#'
#'@rdname proteoNMF
#'@inheritParams  proteoEucDist
#' @examples
#' # metagenes heat maps at all available ranks
#' # additional arguments for `pheatmap`
#' plot_metaNMF(
#'   annot_cols = c("Color", "Alpha", "Shape"),
#'   annot_colnames = c("Lab", "Batch", "WHIM"),
#'
#'   fontsize = 8,
#'   fontsize_col = 5
#' )
#'
#'@import purrr
#'@export
plot_metaNMF <- function (annot_cols = NULL, annot_colnames = NULL, ...) {
  err_msg <- "Don't call the function with arguments `annot_cols` and/or `annot_colnames`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("annot_cols", "annot_colnames"))) stop(err_msg)

  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_protid()
  
  quietly_log <- purrr::quietly(proteoNMF)(id = !!id, task = plotmeta, 
                                           annot_cols = !!annot_cols, 
                                           annot_colnames = !!annot_colnames, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\NMF\\log\\plot_metaNMF_log.csv"), append = TRUE)  
}


