#' Trend analysis
#'
#' @import Mfuzz dplyr purrr rlang Biobase
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>%
trendTest <- function (df, id, col_group, col_order, label_scheme_sub, n_clust,
                       complete_cases, scale_log2r, filepath, filename, ...) {

	mestimate <- function (eset) {
		N <- dim(Biobase::exprs(eset))[[1]]
		D <- dim(Biobase::exprs(eset))[[2]]
		m.sj <- 1 + (1418/N + 22.05) * D^(-2) + (12.33/N + 0.243) *
			D^(-0.0406 * log(N) - 0.1134)
		return(m.sj)
	}

	mfuzz <- function (eset, centers, m, ...) {
		cl <- e1071::cmeans(Biobase::exprs(eset), centers = centers,
		                    method = "cmeans", m = m, ...)
	}

	
	stopifnot(nrow(label_scheme_sub) > 0)

	id <- rlang::as_string(rlang::enexpr(id))
	dots <- rlang::enexprs(...)

	col_group <- rlang::enexpr(col_group)
	col_order <- rlang::enexpr(col_order)
	
	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^.]*$", "", filename)
	
	if (is.null(n_clust)) {
	  n_clust <- c(5:6)
	  fn_prefix <- paste0(fn_prefix, n_clust)
	} else {
	  stopifnot(all(n_clust >= 2) & all(n_clust %% 1 == 0))
	}
	
	if (complete_cases) df <- df[complete.cases(df), ]
	
	df_mean <- t(df) %>%
		data.frame(check.names = FALSE) %>%
		tibble::rownames_to_column("Sample_ID") %>%
		dplyr::left_join(label_scheme_sub, by = "Sample_ID") %>%
		dplyr::filter(!is.na(!!col_group)) %>%
		dplyr::mutate(Group := !!col_group) %>%
		dplyr::group_by(Group) %>%
		dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

	if (rlang::as_string(col_order) %in% names(df_mean)) {
		df_mean <- df_mean %>%
			dplyr::arrange(!!col_order) %>%
			dplyr::rename(Order := !!col_order) %>%
			dplyr::select(-Order, -TMT_Set)
	}

	df_mean <- t(df_mean[, -1]) %>%
		data.frame(check.names = FALSE) %>%
		`colnames<-`(df_mean$Group) %>%
		tibble::rownames_to_column(id) %>%
		dplyr::filter(complete.cases(.[, !grepl(id, names(.))])) %>%
		tibble::column_to_rownames(id)

	purrr::walk(fn_prefix, ~ {
	  n_clust <- gsub(".*_n(\\d+)$", "\\1", .x) %>% as.numeric()
	  filename <- paste0(.x, ".csv")

	  df_fuzzy = new('ExpressionSet', exprs = as.matrix(df_mean))
	  m1 <- mestimate(df_fuzzy)
	  cl <- mfuzz(df_fuzzy, c = n_clust, m = m1)
	  O <- Mfuzz::overlap(cl)
	  
	  res_cl <- data.frame(cl_cluster = cl$cluster) %>%
	    tibble::rownames_to_column() %>%
	    dplyr::rename(!!id := rowname)
	  
	  Levels <- names(df_mean)
	  
	  df_mean %>%
	    tibble::rownames_to_column(id) %>% 
	    left_join(res_cl, by = id) %>% 
	    tidyr::gather(key = variable, value = value, -id, -cl_cluster) %>%
	    dplyr::mutate(variable = factor(variable, levels = Levels)) %>%
	    dplyr::arrange(variable) %>% 
	    write.csv(file.path(filepath, filename), row.names = FALSE)	
	})
}


#' Plots trends
#'
#' @import dplyr rlang purrr ggplot2 RColorBrewer
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>%
plotTrend <- function(id, col_group, col_order, label_scheme_sub, n_clust, filepath, filename, scale_log2r, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  # filename extension will be used but prefix ignored
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  impute_na <- ifelse(all(grepl("_impNA", fn_prefix)), TRUE, FALSE)
  ins <- list.files(path = filepath, pattern = "_n\\d+\\.csv$")
  if (impute_na) ins <- ins %>% .[grepl("_impNA", .)] else ins <- ins %>% .[!grepl("_impNA", .)]
  if (scale_log2r) ins <- ins %>% .[grepl("_Trend_Z", .)] else ins <- ins %>% .[grepl("_Trend_N", .)]
  
  # also possible some files with id = prot_acc and some with id = gene
  # as one might run normPSM(group_pep_by = prot_acc) -> anal_prnTrend 
  # then rerun normPSM(group_pep_by = gene) -> anal_prnTrend
  # save pars for anal_prnTrend, then compare id

	if (is.null(n_clust)) {
	  filelist <- ins
	} else {
	  stopifnot(all(n_clust >= 2) & all(n_clust %% 1 == 0))
	  
	  ins_prefix <- gsub("\\.[^.]*$", "", ins)
	  
	  possible_prefix <- ins_prefix %>% 
	    gsub("(.*)_n\\d+$", "\\1", .) %>% 
	    unique() %>% 
	    paste0("_n", n_clust)
	  
	  ok_prefix <- ins_prefix %>% 
	    .[. %in% possible_prefix]
	  rm(ins_prefix, possible_prefix)
	  
	  filelist <- purrr::map(ok_prefix, ~ list.files(path = filepath, pattern = paste0(.x, "\\.csv$"))) %>% 
	    unlist()
	}
	
  if(purrr::is_empty(filelist)) 
    stop("Missing trend results under ", filepath, 
         "\nCheck the setting in `scale_log2r` for a probable mismatch.", call. = FALSE)

  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  
  my_theme <- theme_bw() + theme(
    axis.text.x  = element_text(angle=60, vjust=0.5, size=24),
    axis.ticks.x  = element_blank(), 
    axis.text.y  = element_text(angle=0, vjust=0.5, size=24),
    axis.title.x = element_text(colour="black", size=24),
    axis.title.y = element_text(colour="black", size=24),
    plot.title = element_text(face="bold", colour="black",
                              size=20, hjust=.5, vjust=.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = '#0868ac', colour = 'red'),
    
    strip.text.x = element_text(size = 24, colour = "black", angle = 0),
    strip.text.y = element_text(size = 24, colour = "black", angle = 90),
    
    plot.margin = unit(c(5.5, 55, 5.5, 5.5), "points"), 
    
    legend.key = element_rect(colour = NA, fill = 'transparent'),
    legend.background = element_rect(colour = NA,  fill = "transparent"),
    legend.position = "none",
    legend.title = element_text(colour="black", size=18),
    legend.text = element_text(colour="black", size=18),
    legend.text.align = 0,
    legend.box = NULL
  )
  

  purrr::walk(filelist, ~ {
    out_nm <- paste0(gsub("\\.csv$", "", .x), ".", fn_suffix)
    src_path <- file.path(filepath, .x)
    df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, 
                            comment.char = "#"), error = function(e) NA)
    
    if (!is.null(dim(df))) {
      message(paste("File loaded:", gsub("\\\\", "/", src_path)))
    } else {
      stop(paste("Non-existed file or directory:", gsub("\\\\", "/", src_path)))
    }
    
    Levels <- label_scheme_sub %>% 
      dplyr::arrange(!!col_order) %>% 
      dplyr::select(!!col_group) %>% 
      unique() %>% 
      unlist()
    
    df <- df %>% 
      dplyr::mutate(variable = factor(variable, levels = Levels))
    
    ymin <- eval(dots$ymin, env = caller_env())
    ymax <- eval(dots$ymax, env = caller_env())
    ybreaks <- eval(dots$ybreaks, env = caller_env())
    ncol <- dots$ncol
    nrow <- dots$nrow
    width <- dots$width
    height <- dots$height
    
    if (is.null(ymin)) ymin <- -2
    if (is.null(ymax)) ymax <- 2
    if (is.null(ybreaks)) ybreaks <- 1
    if (is.null(ncol)) ncol <- 1
    if (is.null(nrow)) nrow <- 2
    
    dots$ymin <- NULL
    dots$ymax <- NULL
    dots$ybreaks <- NULL
    dots$ncol <- NULL
    dots$nrow <- NULL
    
    x_label <- expression("Ratio ("*log[2]*")")
    
    n_clust <- length(unique(df$cl_cluster))
    if (is.null(width)) width <- n_clust * 8 / nrow + 2
    if (is.null(height)) height <- 8 * nrow
    
    dots$width <- NULL
    dots$height <- NULL
    
    p <- ggplot(data = df,
                mapping = aes(x = variable, y = value, group = !!rlang::sym(id))) +
      geom_line(colour = "white", alpha = .25) +
      scale_y_continuous(limits = c(ymin, ymax), breaks = c(ymin, 0, ymax)) +
      labs(title = "", x = "", y = x_label) +
      my_theme
    p <- p + facet_wrap(~ cl_cluster, nrow = nrow, labeller = label_value)
    
    ggsave(file.path(filepath, out_nm),
           p, width = width, height = height, units = "in")
  })
}


#'Clustering of data trends
#'
#'\code{proteoTrend} analyzes and visualizes the trends of protein
#'\code{log2FC}. Users should avoid calling the method directly, but instead use
#'the following wrappers.
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}
#'
#'@inheritParams proteoNMF
#'@inheritParams proteoEucDist
#'@inheritParams proteoCorr
#'@inheritParams info_anal
#'@param n_clust Numeric vector; the number(s) of clusters that data will be
#'  divided into. The default is c(5:6).
#'@param filepath Use system default.
#'@param filename Use system default.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr \cr Additional parameters for use
#'  in \code{plot_} functions: \cr \code{ymin}, the minimum y at \code{log2}
#'  scale; \cr \code{ymax}, the maximum y at \code{log2} scale; \cr
#'  \code{ybreaks}, the breaks in y-axis at \code{log2} scale; \cr \code{nrow},
#'  the number of rows; \cr \code{width}, the width of plot; \cr \code{height},
#'  the height of plot.
#'@return Trend classification and visualization of \code{log2FC}.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnTrend_.R
#'
#'@seealso \code{\link{load_expts}} for a reduced working example in data normalization \cr
#'
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
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPcA}} for PCA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for protein trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for protein NMF analysis and visualization \cr 
#'  
#'  \code{\link{dl_stringdbs}} and \code{\link{getStringDB}} for STRING-DB
#'  
#'@export
proteoTrend <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                         task = "anal", 
                         col_select = NULL, col_group = NULL, col_order = NULL, 
                         impute_na = FALSE, complete_cases = FALSE, 
                         n_clust = NULL, scale_log2r = TRUE, df = NULL, 
                         filepath = NULL, filename = NULL, ...) {

  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0(task, "_pepTrend"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0(task, "_prnTrend"))
    }
    , add = TRUE
  )
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

	id <- rlang::enexpr(id)
	if(id == rlang::expr(c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))) {
		id <- "gene"
	} else {
		id <- rlang::as_string(id)
		stopifnot(id %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
	}
	
	stopifnot(rlang::is_logical(scale_log2r))
	stopifnot(rlang::is_logical(impute_na))
	stopifnot(rlang::is_logical(complete_cases))

	task <- rlang::enexpr(task)
	
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_order <- rlang::enexpr(col_order)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()

	if(!impute_na) complete_cases <- TRUE

	info_anal(id = !!id, col_select = !!col_select, col_group = !!col_group, col_order = !!col_order,
	          scale_log2r = scale_log2r, impute_na = impute_na,
						df = !!df, filepath = !!filepath, filename = !!filename,
						anal_type = "Trend")(n_clust = n_clust, complete_cases = complete_cases, 
						                         task = !!task, ...)
}


#'Trend analysis of protein \code{log2FC}
#'
#'\code{anal_prnTrend} is a wrapper of \code{\link{proteoTrend}} for the trend
#'analysis of protein data
#'
#'@rdname proteoTrend
#'
#'@import purrr
#'@export
anal_prnTrend <- function (...) {
  err_msg <- "Don't call the function with arguments `id`, `anal_type` or `task`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id", "anal_type", "task"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Trend\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)
  
  quietly_log <- purrr::quietly(proteoTrend)(id = !!id, task = anal, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Trend\\log\\anal_prnTrend_log.csv"), append = TRUE)
}


#'Trend visualization
#'
#'\code{plot_prnTrend} is a wrapper of \code{\link{proteoTrend}} for the
#'visualization of the trends in protein data
#'
#'@rdname proteoTrend
#'
#'@import purrr
#'@export
plot_prnTrend <- function (...) {
  err_msg <- "Don't call the function with arguments `id`, `anal_type` or `task`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id", "anal_type", "task"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Trend\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)

  quietly_log <- purrr::quietly(proteoTrend)(id = !!id, task = plot, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Trend\\log\\plot_prnTrend_log.csv"), append = TRUE)
}

