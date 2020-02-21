#' Trend analysis
#'
#' @import dplyr purrr rlang Biobase
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>%
analTrend <- function (df, id, col_group, col_order, label_scheme_sub, n_clust,
                       scale_log2r, complete_cases, impute_na, 
                       filepath, filename, anal_type, ...) {

	stopifnot(nrow(label_scheme_sub) > 0)

	complete_cases <- to_complete_cases(complete_cases = complete_cases, impute_na = impute_na)
	if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
	
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

	if (is.null(dots$m)) {
	  dots$m <- local({
	    N <- nrow(df_mean)
	    D <- ncol(df_mean)
	    m <- 1 + (1418/N + 22.05) * D^(-2) + (12.33/N + 0.243) *
	      D^(-0.0406 * log(N) - 0.1134)
	  })
	}
	
	purrr::walk(fn_prefix, ~ {
	  n_clust <- gsub("^.*_nclust(\\d+).*", "\\1", .x) %>% as.numeric()
	  filename <- paste0(.x, ".txt")

	  args <- c(list(x = as.matrix(df_mean), centers = n_clust), dots)
	  cl <- do.call(cmeans, args) 

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
                      df2 = NULL, filepath, filename, theme, ...) {
  
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

  # find input df2 ---------------------------
  ins <- list.files(path = filepath, pattern = "Trend_[NZ]{1}.*nclust\\d+.*\\.txt$")
  if (purrr::is_empty(ins)) stop("No inputs under ", filepath, call. = FALSE)
  
  if (is.null(df2)) {
    if (impute_na) ins <- ins %>% .[grepl("_impNA", .)] else ins <- ins %>% .[!grepl("_impNA", .)]
    if (scale_log2r) ins <- ins %>% .[grepl("_Trend_Z", .)] else ins <- ins %>% .[grepl("_Trend_N", .)]
    
    if (purrr::is_empty(ins)) 
      stop("No inputs correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, call. = FALSE)
    
  	if (is.null(n_clust)) {
  	  df2 <- ins
  	} else {
  	  stopifnot(all(n_clust >= 2) & all(n_clust %% 1 == 0))
  	  
  	  df2 <- local({
    	  possibles <- ins %>% 
    	    gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", .) %>% 
    	    as.numeric() %>% 
    	    `names<-`(ins)
    	  
    	  n_clust2 <- n_clust %>% .[. %in% possibles]
    	  
    	  df2 <- possibles %>% 
    	    .[. %in% n_clust2] %>% 
    	    names(.)	    
  	  })
  	  
  	  if (purrr::is_empty(df2)) 
  	    stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
  	         " at n_clust = ", paste0(n_clust, collapse = ","), call. = FALSE)
  	}    
  } else {
    local({
      non_exists <- df2 %>% .[! . %in% ins]
      if (!purrr::is_empty(non_exists)) {
        stop("Missing trend file(s): ", purrr::reduce(non_exists, paste, sep = ", "), call. = FALSE)
      }
    })
    if (purrr::is_empty(df2)) stop("File(s) not found under ", filepath, call. = FALSE)
  }

  # prepare output filename ---------------------------	
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Peptide_Trend.*", "\\1", .x)
    })
  } else if (id %in% c("prot_acc", "gene")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Protein_Trend.*", "\\1", .x)
    })
  } else {
    stop("Unknown id = ", id, call. = FALSE)
  }
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  # plot data ---------------------------
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)

  proteoq_trend_theme <- theme_bw() + theme(
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
  
  if (is.null(theme)) theme <- proteoq_trend_theme

  purrr::walk2(df2, custom_prefix, ~ {
    n <- gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", .x) %>% 
      as.numeric()
    
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
    
    if (!purrr::is_empty(dots)) {
      if (any(grepl("^filter_", names(dots)))) {
        stop("Primary `filter_` depreciated; use secondary `filter2_`.")
      }      
    }

    filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
    arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
    dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots)]
    
    df <- df %>% 
      dplyr::filter(group %in% Levels) %>% 
      filters_in_call(!!!filter2_dots) %>% 
      arrangers_in_call(!!!arrange2_dots) %>% 
      dplyr::mutate(group = factor(group, levels = Levels))
    rm(filter2_dots, arrange2_dots)

    if (complete_cases) df <- df %>% .[complete.cases(.), ]
    
    ymin <- eval(dots$ymin, env = caller_env())
    ymax <- eval(dots$ymax, env = caller_env())
    ybreaks <- eval(dots$ybreaks, env = caller_env())
    ncol <- dots$ncol
    nrow <- dots$nrow
    width <- dots$width
    height <- dots$height
    units <- dots$units
    color <- dots$color
    alpha <- dots$alpha
    
    if (is.null(ymin)) ymin <- -2
    if (is.null(ymax)) ymax <- 2
    if (is.null(ybreaks)) ybreaks <- 1
    if (is.null(ncol)) ncol <- 1
    if (is.null(nrow)) nrow <- 2
    if (is.null(color)) color <- "#f0f0f0"
    if (is.null(alpha)) alpha <- .25
    
    dots$ymin <- NULL
    dots$ymax <- NULL
    dots$ybreaks <- NULL
    dots$ncol <- NULL
    dots$nrow <- NULL
    dots$color <- NULL
    dots$alpha <- NULL
    
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
      geom_line(colour = color, alpha = alpha) + 
      # coord_cartesian(ylim = c(ymin, ymax)) + 
      scale_y_continuous(limits = c(ymin, ymax), breaks = c(ymin, 0, ymax)) +
      labs(title = "", x = "", y = x_label) +
      theme
    p <- p + facet_wrap(~ cluster, nrow = nrow, labeller = label_value)
    
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
#'@inheritParams prnCorr_logFC
#'@inheritParams anal_prnNMF
#'@param n_clust Numeric vector; the number(s) of clusters that data will be
#'  divided into. At the NULL default, it will be determined by the number of
#'  data entries. The \code{n_clust} overwrites the augument \code{centers} in
#'  \code{\link[e1071]{cmeans}}.
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param filepath Use system default.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering aganist
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional arguments for
#'  \code{\link[e1071]{cmeans}} by noting that: \cr \code{centers} is replaced
#'  with \code{n_clust} \cr \code{m} is according to Schwaemmle and Jensen if
#'  not provided; \cr \code{x} is disabled with input data being determined
#'  automatically
#'@return Fuzzy c-mean classification of \code{log2FC}.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnTrend_.R
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
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@export
#'@references \code{Schwaemmle and Jensen, Bioinformatics,Vol. 26 (22),
#'  2841-2848, 2010. \cr J. C. Bezdek (1981). Pattern recognition with fuzzy
#'  objective function algorithms. New York: Plenum. \cr }
anal_prnTrend <- function (col_select = NULL, col_group = NULL, col_order = NULL, n_clust = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                           df = NULL, filepath = NULL, filename = NULL, ...) {
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
  
  check_dots(c("id", "df2", "anal_type"), ...)
  
  err_msg1 <- "Do not use argument `x`; input data will be determined automatically.\n"
  err_msg2 <- "Do not use argument `centers`; instead use `n_clust`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("x"))) stop(err_msg1, call. = FALSE)
  if (any(names(rlang::enexprs(...)) %in% c("centers"))) stop(err_msg2, call. = FALSE)

  dir.create(file.path(dat_dir, "Protein\\Trend\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = !!col_group, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
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
#'@inheritParams anal_prnNMF
#'@inheritParams prnCorr_logFC
#'@inheritParams prnHist
#'@param df2 Character vector or string; the name(s) of secondary data file(s).
#'  An informatic task, i.e. \code{anal_prnTrend(...)} against a primary
#'  \code{df} generates secondary files such as
#'  \code{Protein_Trend_Z_nclust6.txt} etc. See also the \code{Auguments}
#'  section or \code{\link{prnHist}} for the description of a primary \code{df}.
#'@param scale_log2r Logical; at the TRUE default, input files with
#'  \code{_Z[...].txt} in name will be used. Otherwise, files with
#'  \code{_N[...].txt} in name will be taken. An error will be thrown if no
#'  files are matched under given conditions.
#'@param impute_na Logical; at the TRUE default, input files with
#'  \code{_impNA[...].txt} in name will be loaded. Otherwise, files without
#'  \code{_impNA} in name will be taken. An error will be thrown if no files are
#'  matched under given conditions.
#'@param n_clust Numeric vector; the cluster ID(s) corresponding to
#'  \code{\link{anal_prnTrend}} for visualization. At the NULL default, all
#'  available cluster IDs will be used.
#'@param ...  \code{filter2_}: Variable argument statements for the row
#'  filtration aganist data in secondary file(s) of
#'  \code{[...]Protein_Trend_[...].txt}. See also \code{\link{prnGSPAHM}} for
#'  the format of \code{filter2_} statements. \cr \cr \code{arrange2_}: Variable
#'  argument statements for the row ordering aganist data in secondary file(s)
#'  of \code{[...]Protein_Trend_[...].txt}. See also \code{\link{prnGSPAHM}} for
#'  the format of \code{arrange2_} statements. \cr \cr Additional parameters for
#'  use in \code{plot_} functions: \cr \code{ymin}, the minimum y at \code{log2}
#'  scale; \cr \code{ymax}, the maximum y at \code{log2} scale; \cr
#'  \code{ybreaks}, the breaks in y-axis at \code{log2} scale; \cr \code{nrow},
#'  the number of rows; \cr \code{width}, the width of plot; \cr \code{height},
#'  the height of plot; \cr \code{color}, the color of lines; \cr \code{alpha},
#'  the transparency of lines.
#'@import purrr
#'
#'@example inst/extdata/examples/prnTrend_.R
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
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB
#'
#'@export
plot_prnTrend <- function (col_select = NULL, col_order = NULL, n_clust = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                           df2 = NULL, filename = NULL, theme = NULL, ...) {
  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\Trend\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  filename <- rlang::enexpr(filename)
  df2 <- rlang::enexpr(df2)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = NULL, df2 = !!df2, filepath = NULL, filename = !!filename,
            anal_type = "Trend_line")(n_clust = n_clust, theme = theme, ...)
}

