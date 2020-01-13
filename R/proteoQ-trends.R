#' Trend analysis
#'
#' @import Mfuzz dplyr purrr rlang Biobase
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>%
analTrend <- function (df, id, col_group, col_order, label_scheme_sub, n_clust,
                       filepath, filename, ...) {

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
	  n_clust <- local({
	    nrow <- nrow(df)
	    if (nrow >= 5000) n_clust <- c(5:6) else if (nrow >= 2000) n_clust <- c(3:5) else n_clust <- 3 
	  })
	  fn_prefix <- paste0(fn_prefix, n_clust)
	} else {
	  stopifnot(all(n_clust >= 2) & all(n_clust %% 1 == 0))
	}
	
	suppressWarnings(
	  df_mean <- t(df) %>%
	    data.frame(check.names = FALSE) %>%
	    tibble::rownames_to_column("Sample_ID") %>%
	    dplyr::left_join(label_scheme_sub, by = "Sample_ID") %>%
	    dplyr::filter(!is.na(!!col_group)) %>%
	    dplyr::mutate(Group := !!col_group) %>%
	    dplyr::group_by(Group) %>%
	    dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
	)

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
	  n_clust <- gsub("^.*_nclust(\\d+).*", "\\1", .x) %>% as.numeric()
	  filename <- paste0(.x, ".txt")

	  df_fuzzy = new('ExpressionSet', exprs = as.matrix(df_mean))
	  m1 <- mestimate(df_fuzzy)
	  cl <- mfuzz(df_fuzzy, c = n_clust, m = m1)
	  O <- Mfuzz::overlap(cl)
	  
	  res_cl <- data.frame(cluster = cl$cluster) %>%
	    tibble::rownames_to_column() %>%
	    dplyr::rename(!!id := rowname)
	  
	  Levels <- names(df_mean)
	  
	  df_mean %>%
	    tibble::rownames_to_column(id) %>% 
	    left_join(res_cl, by = id) %>% 
	    tidyr::gather(key = variable, value = value, -id, -cluster) %>%
	    dplyr::mutate(variable = factor(variable, levels = Levels)) %>%
	    dplyr::arrange(variable) %>% 
	    dplyr::rename(id = !!rlang::sym(id), group = variable, log2FC = value) %>% 
	    write.table(file.path(filepath, filename), sep = "\t", 
	                col.names = TRUE, row.names = FALSE, quote = FALSE)	
	})
}


#' Plots trends
#'
#' @import dplyr rlang purrr ggplot2 RColorBrewer
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>%
plotTrend <- function(id, col_group, col_order, label_scheme_sub, n_clust, 
                      scale_log2r, complete_cases, impute_na, 
                      filepath, filename, ...) {
  
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

  ins <- list.files(path = filepath, pattern = "Trend_[NZ]{1}.*nclust\\d+.*\\.txt$")

  if (impute_na) ins <- ins %>% .[grepl("_impNA", .)] else ins <- ins %>% .[!grepl("_impNA", .)]
  if (scale_log2r) ins <- ins %>% .[grepl("_Trend_Z", .)] else ins <- ins %>% .[grepl("_Trend_N", .)]
  
	if (is.null(n_clust)) {
	  filelist <- ins
	} else {
	  stopifnot(all(n_clust >= 2) & all(n_clust %% 1 == 0))
	  
	  filelist <- local({
  	  possible <- ins %>% 
  	    gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", .) %>% 
  	    as.numeric() %>% 
  	    `names<-`(ins)
  	  
  	  n_clust2 <- n_clust %>% .[. %in% possible]
  	  
  	  filelist <- possible %>% 
  	    .[. %in% n_clust2] %>% 
  	    names(.)	    
	  })
	}
	
  custom_prefix <- purrr::map_chr(filelist, ~ {
    gsub("(.*_{0,1})Protein_Trend.*", "\\1", .x)
  })

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  if (purrr::is_empty(filelist)) 
    stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
         " at n_clust = ", paste0(n_clust, collapse = ","), call. = FALSE)

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

  purrr::walk2(filelist, custom_prefix, ~ {
    n <- gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", .x) %>% 
      as.numeric()
    
    fn_suffix <- "png" # png for now
    out_nm <- paste0(.y, fn_prefix, "_nclust", n, ".", fn_suffix)
    src_path <- file.path(filepath, .x)

    df <- tryCatch(
      readr::read_tsv(src_path, col_types = cols(id = col_character(), group = col_factor())), 
      error = function(e) NA
    )

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
    
    filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
    dots <- dots %>% .[! . %in% filter_dots]
    
    df <- df %>% 
      filters_in_call(!!!filter_dots) %>% 
      dplyr::mutate(group = factor(group, levels = Levels))

    if (complete_cases) df <- df %>% .[complete.cases(.), ]
    
    ymin <- eval(dots$ymin, env = caller_env())
    ymax <- eval(dots$ymax, env = caller_env())
    ybreaks <- eval(dots$ybreaks, env = caller_env())
    ncol <- dots$ncol
    nrow <- dots$nrow
    width <- dots$width
    height <- dots$height
    units <- dots$units
    
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
    
    n_clust <- length(unique(df$cluster))
    if (is.null(width)) width <- n_clust * 8 / nrow + 2
    if (is.null(height)) height <- 8 * nrow
    if (is.null(units)) units <- "in"
    
    dots$width <- NULL
    dots$height <- NULL
    dots$units <- NULL
    
    p <- ggplot(data = df,
                mapping = aes(x = group, y = log2FC, group = id)) +
      geom_line(colour = "white", alpha = .25) + 
      # coord_cartesian(ylim = c(ymin, ymax)) + 
      scale_y_continuous(limits = c(ymin, ymax), breaks = c(ymin, 0, ymax)) +
      labs(title = "", x = "", y = x_label) +
      my_theme
    p <- p + facet_wrap(~ cluster, nrow = nrow, labeller = label_value)
    
    # ggsave(filename = file.path(filepath, out_nm),
    #        plot = p, width = width, height = height, units = units)
    
    gg_args <- c(filename = file.path(filepath, gg_imgname(out_nm)), 
                 width = width, height = height, units = units, dots) %>% 
      do.call(ggsave, .)
  }, complete_cases = complete_cases)
}



#'Trend analysis of protein data
#'
#'\code{anal_prnTrend} applies the soft clustering algorithm in
#'\code{\link[e1071]{cmeans}} for the trend analysis of protein \code{log2FC}.
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}
#'
#'@inheritParams proteoNMF
#'@inheritParams proteoEucDist
#'@inheritParams proteoCorr
#'@inheritParams info_anal
#'@param n_clust Numeric vector; the number(s) of clusters that data will be
#'  divided into. At the NULL default, it will be determined by the number of
#'  data entries.
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param filepath Use system default.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}.
#'@return Trend classification of \code{log2FC}.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnTrend_.R
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
#'@export
anal_prnTrend <- function (col_select = NULL, col_group = NULL, col_order = NULL, 
                           scale_log2r = TRUE, impute_na = FALSE, complete_cases = FALSE, 
                           df = NULL, filepath = NULL, filename = NULL, n_clust = NULL, ...) {
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_pepTrend"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_prnTrend"))
    }
    , add = TRUE
  )  
  
  err_msg <- "Don't call the function with arguments `id` or `anal_type`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id", "anal_type"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Trend\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(rlang::is_logical(scale_log2r), 
            rlang::is_logical(impute_na), 
            rlang::is_logical(complete_cases))
  
  col_select <- rlang::enexpr(col_select)
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = !!col_group, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "Trend")(n_clust = n_clust, ...)
}


#'Visualization of trend results
#'
#'\code{plot_prnTrend} plots the trends of protein expressions from
#'\code{\link{anal_prnTrend}}.
#'
#'The function reads \code{Protein_Trend_[...].txt} files under the
#'\code{...\\Protein\\Trend} directory.
#'
#'@section \code{Protein_Trend_[...].txt}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Descrption}\cr id \tab a gene name
#'  or an acession number for protein data \cr cluster \tab a cluster ID
#'  assigned to an \code{id} \cr group \tab a name of the sample group for a
#'  \code{id} \cr log2FC \tab the mean \code{log2FC} of an \code{id} under a
#'  \code{group} at a given \code{cluster} \cr }
#'
#'@inheritParams proteoNMF
#'@inheritParams proteoCorr
#'@inheritParams  proteoEucDist
#'@inheritParams proteoHM
#'@inheritParams anal_prnTrend
#'@param scale_log2r Logical; at the TRUE default, files with \code{_Z[...].txt}
#'  in name will be used. Otherwise, files with \code{_N[...].txt} in name will
#'  be taken. An error will be shown if no files are matched under given
#'  conditions.
#'@param impute_na Logical; at the TRUE default, files with
#'  \code{_impNA[...].txt} in name will be loaded. Otherwise, files without
#'  \code{_impNA} in name will be taken. An error will be shown if no files are
#'  matched under given conditions.
#'@param n_clust Numeric vector; the cluster ID(s) corresponding to
#'  \code{\link{anal_prnTrend}} for visualization. At the NULL default, all
#'  available cluster IDs will be used.
#'
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data in \code{Protein_Trend_[...].txt}; also see \code{\link{normPSM}}. \cr
#'  \cr Additional parameters for use in \code{plot_} functions: \cr
#'  \code{ymin}, the minimum y at \code{log2} scale; \cr \code{ymax}, the
#'  maximum y at \code{log2} scale; \cr \code{ybreaks}, the breaks in y-axis at
#'  \code{log2} scale; \cr \code{nrow}, the number of rows; \cr \code{width},
#'  the width of plot; \cr \code{height}, the height of plot.
#'@import purrr
#'
#'@example inst/extdata/examples/prnTrend_.R
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
#'@export
plot_prnTrend <- function (col_group = NULL, col_order = NULL, 
                           scale_log2r = TRUE, impute_na = FALSE, complete_cases = FALSE, 
                           filename = NULL, n_clust = NULL, ...) {
  
  err_msg <- "Duplicated arguments in `id` or `anal_type`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id", "anal_type"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Trend\\log"), recursive = TRUE, showWarnings = FALSE)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(rlang::is_logical(scale_log2r), 
            rlang::is_logical(impute_na), 
            rlang::is_logical(complete_cases))
  
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = Select, col_group = !!col_group, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = NULL, filepath = NULL, filename = !!filename,
            anal_type = "Trend_line")(n_clust = n_clust, ...)
}


