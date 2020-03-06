#' Correlation plots
#' 
#' @param data_select The type of data to be selected, for example, logFC or logInt.
#' @inheritParams prnCorr_logFC
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import stringr dplyr ggplot2 GGally rlang
#' @importFrom magrittr %>%
plotCorr <- function (df = NULL, id, anal_type, data_select, col_select = NULL, col_order = NULL,
                      label_scheme_sub = label_scheme_sub, 
                      scale_log2r = scale_log2r, complete_cases = complete_cases, 
                      filepath = filepath, filename = filename, ...) {

  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  xmin <- eval(dots$xmin, env = caller_env()) # `xmin = -1` is `language`
  xmax <- eval(dots$xmax, env = caller_env()) # `xmax = +1` is `language`
  xbreaks <- eval(dots$xbreaks, env = caller_env())
  width <- eval(dots$width, env = caller_env())
  height <- eval(dots$height, env = caller_env())

  dots <- dots %>% 
    .[! names(.) %in% c("xmin", "xmax", "xbreaks", "width", "height")]
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)

	col_select <- rlang::enexpr(col_select)
	col_order <- rlang::enexpr(col_order)

	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
	             sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) 
	
	if (data_select == "logFC") {
	  df <- df %>% .$log2R
	  y_label <- x_label <- expression("Ratio ("*log[2]*")")
	  if (is.null(xmin)) xmin <- -2
	  if (is.null(xmax)) xmax <- 2
	  if (is.null(xbreaks)) xbreaks <- 1
	} else if (data_select == "logInt") {
	  df <- df %>% .$Intensity %>% log10()
	  y_label <- x_label <- expression("Intensity ("*log[10]*")")
	  if (is.null(xmin)) xmin <- 3.5
	  if (is.null(xmax)) xmax <- 6
	  if (is.null(xbreaks)) xbreaks <- 1
	} else {
	  stop("`data_select` nees to be either`logFC` or `logInt`.", call. = FALSE)
	}
	
	if (is.null(width)) width <- 1.4 * length(label_scheme_sub$Sample_ID)
	if (is.null(height)) height <- width
	if (ncol(df) > 44) stop("Maximal number of samples for correlation plots is 44.", call. = FALSE)

	if (dplyr::n_distinct(label_scheme_sub[[col_order]]) == 1) {
		df <- df[, order(names(df))]
	} else {
	  corrplot_orders <- label_scheme_sub %>%
	    dplyr::select(Sample_ID, !!col_select, !!col_order) %>%
	    dplyr::filter(!is.na(!!col_order)) %>%
	    unique(.) %>%
	    dplyr::arrange(!!col_order)
	  
	  df <- df[, as.character(corrplot_orders$Sample_ID), drop = FALSE]
	}

	plot_corr_sub(df = df, xlab = x_label, ylab = y_label,
	              filename = filename, filepath = filepath,
	              xmin = xmin, xmax = xmax, xbreaks = xbreaks, width = width, height = height, !!!dots)
}


#' Make correlation plots
#' 
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param xmin minimal x.
#' @param xmax maximal x.
#' @param xbreaks breaks on x-axis.
#' @param width plot width
#' @param height plot height
#' @param ... additional arguments for ggsave.
#' @inheritParams info_anal
#' 
#' @import stringr dplyr ggplot2 GGally purrr rlang
#' @importFrom magrittr %>%
plot_corr_sub <- function (df, xlab, ylab, filename, filepath, 
                           xmin, xmax, xbreaks, width, height, ...) {
                           
  my_fn <- function(data, mapping, method = "lm", ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(alpha = 0.3, size = .1) +
      geom_smooth(method = method, ...)
    p
  }

  lm_with_cor <- function(data, mapping, ..., method = "pearson") {
    x <- data[[deparse(mapping$x)]]
    y <- data[[deparse(mapping$y)]]
    cor <- cor(x, y, method = method)

    ggally_smooth_lm(data, mapping, ...) +
      ggplot2::geom_label(
        data = data.frame(
          x = min(x, na.rm = TRUE),
          y = max(y, na.rm = TRUE),
          lab = round(cor, digits = 3)
        ),
        mapping = ggplot2::aes(x = x, y = y, label = lab, color = NULL),
        hjust = 0, vjust = 1,
        size = 5, fontface = "bold"
      )
  }

  panel_cor = function(x, y, digits = 2, prefix = "", cex.cor, ...){
    usr = par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = abs(cor(x, y, use = 'complete.obs'))
    txt = format(c(r, 0.123456789),digits = digits)[1]
    txt = paste(prefix, txt, sep = '')
    if(missing(cex.cor)) cex.cor = 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor*(1+r)/2)
    bg = "transparent"
  }

  panel_hist = function(x, ...){
    usr = par('usr')
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h = hist(x, plot = FALSE)
    breaks = h$breaks
    nB = length(breaks)
    y=h$counts
    y=y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = 'white', ...)
    bg = "transparent"
  }

  panel_lm = function(x, y, col = par('col'), bg = NA, pch = par('pch'),
                      cex = 1, col.smooth = 'black', ...){
    points(x, y, pch = pch, col = "red",  bg = bg, cex = cex)
    # abline(stats::lm(y~x), col=col.smooth,...)
    bg = "transparent"
  }

  my_lower <- function(data, mapping, method = "lm", ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(size = .02, alpha = .5) +
      geom_abline(alpha = .5, linetype = "dashed", color = "gray", size = 1) +
      geom_smooth(size = 1, method = method, ...)
    p
  }

  my_lower_no_sm <- function(data, mapping, method = "lm", ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(size = .02, alpha = .5) +
      geom_abline(alpha = .5, linetype = "dashed", color = "gray", size = 1)
    p
  }

  my_diag <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_density(fill = "#fec44f", size = .02, alpha = .5, adjust = 3)
    p
  }

  my_custom_cor <- function(data, mapping, color = I("grey50"), sizeRange = c(1, 4), ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)

    ct <- cor.test(x, y)
    sig <- symnum(
      ct$p.value, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )

    r <- unname(ct$estimate)
    rt <- format(r, digits=2)[1]

    cex <- max(sizeRange)

    percent_of_range <- function(percent, range) {
      percent * diff(range) + min(range, na.rm = TRUE)
    }

    ggally_text(
      label = as.character(rt),
      mapping = aes(),
      xP = 0.5, yP = 0.5,
      size = I(percent_of_range(cex * abs(r), sizeRange)),
      color = color,
      ...
    ) +
      ## add the sig stars
      # geom_text(
      #   aes_string(
      #     x = 0.8,
      #     y = 0.8
      #   ),
      #   label = sig,
      #   size = I(cex),
      #   color = color,
      #   ...
      # ) +
    theme_classic() +
      theme(
        panel.background = element_rect(
          color = color,
          linetype = "longdash"
        ),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()
      )
  }

  my_theme <- theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 0, vjust = 0.5, size = 10),
      axis.text.y  = element_text(angle = 0, vjust = 0.5, size = 10),
      axis.title.x = element_text(colour = "black", size = 16),
      axis.title.y = element_text(colour = "black", size = 16),
      plot.title = element_text(face = "bold", colour = "black", size = 20,
                                hjust = 0.5, vjust = 0.5),

      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),

      legend.key = element_rect(colour = NA, fill = 'transparent'),
      legend.background = element_rect(colour = NA,  fill = "transparent"),
      legend.position = "none",
      legend.title = element_text(colour = "black", size = 16),
      legend.text = element_text(colour = "black", size = 14),
      legend.text.align = 0,
      legend.box = NULL
    )

  dots <- rlang::enexprs(...)

  ncol <- ncol(df)
  
  if (is.null(dots$dpi)) {
    dpi <- 300
  } else {
    dpi <- eval(dots$dpi, env = caller_env())
    dots$dpi <- NULL
  }

  p1 <- ggpairs(df, columnLabels = as.character(names(df)), labeller = label_wrap_gen(10),
                title = "", xlab = xlab, ylab = ylab, lower = list(continuous = my_lower_no_sm),
                upper = list(continuous = my_custom_cor, digits = 2))
  p2 <- ggcorr(df, label = TRUE, label_round = 2)

  g2 <- ggplotGrob(p2)
  colors <- g2$grobs[[6]]$children[[3]]$gp$fill

  idx <- 1
  for (k1 in 1:(ncol-1)) { # row
    for (k2 in (k1+1):ncol) { # column
      plt <- getPlot(p1, k1, k2) +
        theme(panel.background = element_rect(fill = colors[idx], color = "white"),
              panel.grid.major = element_line(color = colors[idx]))
      p1 <- putPlot(p1, plt, k1, k2)
      idx <- idx + 1
    }
  }

  idx <- 1
  for (k1 in 1:(ncol-1)) { # column
    for (k2 in (k1+1):ncol) { # row
      plt <- getPlot(p1, k2, k1) +
        theme_bw() +
          theme(panel.background = element_rect(colour = NA,  fill = "transparent"),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()
        )
      p1 <- putPlot(p1, plt, k2, k1)
      idx <- idx + 1
    }
  }

  idx <- 1
  for (k1 in 1:ncol) { # row
    plt <- getPlot(p1, k1, k1) +
      theme_bw() +
        theme(panel.background = element_rect(colour = NA,  fill = "transparent"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()
      )
    p1 <- putPlot(p1, plt, k1, k1)
    idx <- idx + 1
  }

  for (x in 2:ncol) {
    for (y in 1:(x-1)) {
      p1[x, y] <- p1[x, y] +
        scale_x_continuous(limits = c(xmin-.2, xmax+.2), breaks = c(xmin, 0, xmax)) +
        scale_y_continuous(limits = c(xmin-.2, xmax+.2), breaks = c(xmin, 0, xmax))
    }
  }

  for (x in 1:ncol)
    p1[x, x] <- p1[x, x] +
      scale_x_continuous(limits = c(xmin-.2, xmax+.2), breaks = c(xmin, 0, xmax))

  p1 <- p1 +
    theme(plot.title = element_text(face = "bold", colour = "black", size = 20,
                                    hjust = 0.5, vjust = 0.5))

  quietly_log <- purrr::quietly(ggsave)(file.path(filepath, gg_imgname(filename)), 
                                        plot = p1, width = width, height = height, dpi = dpi, units = "in")
}


#'Correlation plots
#'
#'\code{pepCorr_logFC} plots Pearson correlation for peptide \code{logFC}. 
#'data.
#'
#'@rdname prnCorr_logFC
#'
#'@import purrr
#'@export
pepCorr_logFC <- function (col_select = NULL, col_order = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                           df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Corrplot")(data_select = "logFC", ...)
}


#'Correlation plots
#'
#'\code{pepCorr_logInt} plots Pearson correlation of the \code{log10} intensity
#'of reporter ions for peptide data.
#'
#'@rdname prnCorr_logFC
#'
#'@import purrr
#'@export
pepCorr_logInt <- function (col_select = NULL, col_order = NULL, 
                            scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                            df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Corrplot")(data_select = "logInt", ...)
}


#'Correlation plots
#'
#'\code{prnCorr_logFC} plots Pearson correlation for protein \code{logFC}. 
#'
#'The function matches the current \code{id} to the grouping argument in the
#'latest \code{call} to \code{\link{normPSM}}. See also \code{\link{prnHist}}
#'for details.
#'
#'@inheritParams prnHist
#'@inheritParams prnMDS
#'@param col_order Character string to a column key in \code{expt_smry.xlsx}.
#'  Numeric values under which will be used for the left-to-right arrangement of
#'  samples in graphic outputs or top-to-bottom arrangement in text outputs. At
#'  the NULL default, the column key \code{Order} will be used. If values under
#'  column \code{Order} are left blank, samples will be ordered by their names.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  Additional parameters for
#'  plotting: \cr \code{width}, the width of plot \cr \code{height}, the height
#'  of plot \cr \code{xmin}, the minimum \eqn{x} of logFC or intensity \cr
#'  \code{xmax}, the maximum \eqn{x} of logFC data or intensity data \cr
#'  \code{xbreaks}, the breaks on \eqn{x} axis; the same breaks will be applied
#'  to \eqn{y} axis.
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@example inst/extdata/examples/prnCorr_.R
#'
#'@return Correlation plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnCorr_logFC <- function (col_select = NULL, col_order = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                           df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Corrplot")(data_select = "logFC", ...)
}


#'Correlation Plots
#'
#'\code{prnCorr_logInt} plots Pearson correlation of the \code{log10} intensity
#'of reporter ions for protein data.
#'
#'
#'@rdname prnCorr_logFC
#'
#'@import purrr
#'@export
prnCorr_logInt <- function (col_select = NULL, col_order = NULL, 
                            scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                            df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Corrplot")(data_select = "logInt", ...)
}


#' Entrez IDs of human molecular signatures.
#'
#' A dataset containing human entrez IDs by the gene sets of molecular signatures.
#'
#' @format A list at a length of 5501 
#' \describe{
#'   \item{hs_...}{human entries}
#' }
#' @source \url{https://www.gsea-msigdb.org/gsea/index.jsp}
"c2_msig_hs"


#' Entrez IDs of mouse molecular signatures.
#'
#' A dataset containing mouse entrez IDs by the gene sets of molecular signatures.
#'
#' @format A list at a length of 5501 
#' \describe{
#'   \item{mm_...}{mouse entries}
#' }
#' @source \url{https://www.gsea-msigdb.org/gsea/index.jsp}
"c2_msig_mm"


#' Entrez IDs of rat molecular signatures.
#'
#' A dataset containing rat entrez IDs by the gene sets of molecular signatures.
#'
#' @format A list at a length of 5501 
#' \describe{
#'   \item{rn_...}{rat entries}
#' }
#' @source \url{https://www.gsea-msigdb.org/gsea/index.jsp}
"c2_msig_rn"


#' Entrez IDs of human gene-ontology annotations.
#'
#' A dataset containing human entrez IDs by the gene sets of GO.
#'
#' @format A list at a length of 18320 
#' \describe{
#'   \item{hs_...}{human entries}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"go_sets_hs"


#' Entrez IDs of mouse gene-ontology annotations.
#'
#' A dataset containing mouse entrez IDs by the gene sets of GO.
#'
#' @format A list at a length of 18364 
#' \describe{
#'   \item{mm_...}{mouse entries}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"go_sets_mm"


#' Entrez IDs of rat gene-ontology annotations.
#'
#' A dataset containing rat entrez IDs by the gene sets of GO.
#'
#' @format A list at a length of 18726 
#' \describe{
#'   \item{rn_...}{rat entries}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"go_sets_rn"


#' Entrez IDs of human KEGG annotations.
#'
#' A dataset containing human entrez IDs by the gene sets of KEGG pathways.
#'
#' @format A list at a length of 229 
#' \describe{
#'   \item{hsa...}{human entries}
#' }
#' @source \url{https://www.genome.jp/kegg/pathway.html}
"kegg_sets_hs"


#' Entrez IDs of mouse KEGG annotations.
#'
#' A dataset containing mouse entrez IDs by the gene sets of KEGG pathways.
#'
#' @format A list at a length of 225 
#' \describe{
#'   \item{mmu...}{mouse entries}
#' }
#' @source \url{https://www.genome.jp/kegg/pathway.html}
"kegg_sets_mm"


#' Entrez IDs of rat KEGG annotations.
#'
#' A dataset containing rat entrez IDs by the gene sets of KEGG pathways.
#'
#' @format A list at a length of 225 
#' \describe{
#'   \item{mmu...}{rat entries}
#' }
#' @source \url{https://www.genome.jp/kegg/pathway.html}
"kegg_sets_rn"


#' Lookups of human or mouse kinases.
#'
#' A dataset containing the \code{refseq}, \code{uniprot}, \code{gene names} of
#' human or mouse kinases.
#'
#' @format A data frame with 2050 rows and 9 variables:
#' \describe{
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{gene names}
#'   \item{kin_attr}{the attribute being a kinase or not}
#'   \item{kin_class}{the class of a kinase; for example: TK, tyrosine kinase}
#'   \item{kin_order}{the number order of a kinase class by the kinase tree}
#'   \item{uniprot_acc}{UniProt accession number}
#'   \item{uniprot_id}{UniProt ID}
#'   \item{prot_desc}{protein description}
#'   \item{entrez}{Entrez ID}
#' }
#' @source \url{http://kinase.com/human/kinome/phylogeny.html}
"kinase_lookup"


#' BioMart-Ensembl databases for human proteins.
#' 
#' A S4 object from \link[biomaRt]{useMart} at "dataset = hsapiens_gene_ensembl". 
#' @source \url{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}
"mart_hs"


#' BioMart-Ensembl for mouse proteins.
#' 
#' A S4 object from \link[biomaRt]{useMart} at "dataset = mmusculus_gene_ensembl". 
#' @source \url{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}
"mart_mm"


#' BioMart-Ensembl for rat proteins.
#' 
#' A S4 object from \link[biomaRt]{useMart} at "dataset = rnorvegicus_gene_ensembl". 
#' @source \url{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}
"mart_rn"


#' Lookups of cRAP proteins. 
#'
#' A dataset containing the IDs of cRAP proteins. 
#'
#' @format A data frame with 195 rows and 9 variables:
#' \describe{
#'   \item{uniprot_acc}{UniProt accession number}
#'   \item{uniprot_id}{UniProt ID}
#'   \item{status}{the review status according to UniProt}
#'   \item{prot_desc}{protein description}
#'   \item{gene}{gene names}
#'   \item{organism}{the organism attribute according to UniProt}
#'   \item{length}{the number of amino acid residues under a proposed protein}
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{entrez}{Entrez ID (not currently used)}
#' }
#' @source \url{https://www.thegpm.org/crap/}
"prn_annot_crap"


#' Lookups among RefSeq accessions, Entrez IDs and gene names for human
#' proteins.
#'
#' A dataset containing human \code{refseq} accessions, \code{entrez} IDs and
#' \code{gene} names.
#'
#' @format A data frame with 31833 rows and 4 variables:
#' \describe{
#'   \item{entrez}{Entrez ID}
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{gene names}
#' }
"refseq_entrez_hs"


#' Lookups among RefSeq accessions, Entrez IDs and gene names for mouse
#' proteins.
#'
#' A dataset containing mouse \code{refseq} accessions, \code{entrez} IDs and
#' \code{gene} names.
#'
#' @format A data frame with 24906 rows and 4 variables: \describe{
#'   \item{entrez}{Entrez ID} 
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{gene names} }
"refseq_entrez_mm"


#' Lookups among RefSeq accessions, Entrez IDs and gene names for rat
#' proteins.
#'
#' A dataset containing mouse \code{refseq} accessions, \code{entrez} IDs and
#' \code{gene} names.
#'
#' @format A data frame with 158959 rows and 3 variables: \describe{
#'   \item{entrez}{Entrez ID} 
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{gene names} }
"refseq_entrez_rn"


#' Lookups between UniProt accessions and Entrez IDs for human proteins.
#'
#' A dataset containing human \code{entrez} IDs and \code{UniProt} accessions.
#'
#' @format A data frame with 33194 rows and 2 variables:
#' \describe{
#'   \item{entrez}{Entrez ID}
#'   \item{uniprot_acc}{UniProt accession number}
#' }
"uniprot_entrez_hs"


#' Lookups between UniProt accessions and Entrez IDs for mouse proteins.
#'
#' A dataset containing mouse \code{entrez} IDs and \code{UniProt} accessions.
#'
#' @format A data frame with 33231 rows and 2 variables:
#' \describe{
#'   \item{entrez}{Entrez ID}
#'   \item{uniprot_acc}{UniProt accession number}
#' }
"uniprot_entrez_mm"


#' Lookups between UniProt accessions and Entrez IDs for rat proteins.
#'
#' A dataset containing rat \code{entrez} IDs and \code{UniProt} accessions.
#'
#' @format A data frame with 21197 rows and 2 variables:
#' \describe{
#'   \item{entrez}{Entrez ID}
#'   \item{uniprot_acc}{UniProt accession number}
#' }
"uniprot_entrez_rn"


#' A list of UniProt taxon IDs and Species names.
#'
#' The dataset is based on UniProt.ws::availableUniprotSpecies().
#'
#' @format A data frame with 10818 rows and 2 variables:
#' \describe{
#'   \item{taxid}{UniProt taxon ID}
#'   \item{organism}{UniProt Species name}
#' }
"uniprot_species"

#' Helper to create `db_path`
#' 
#' @inheritParams prepGO
create_db_path <- function (db_path) {
  if (!fs::dir_exists(db_path)) {
    new_db_path <- fs::path_expand_r(db_path)
    new_db_path2 <- fs::path_expand(db_path)
    
    if (fs::dir_exists(new_db_path)) {
      db_path <- new_db_path
    } else if (fs::dir_exists(new_db_path2)) {
      db_path <- new_db_path2
    } else {
      dir.create(file.path(db_path, "cache"), recursive = TRUE, showWarnings = FALSE)
      db_path <- new_db_path
    }
  }
}


#' Helper to save `obo` without header
#' 
#' @param fn_obo filename according to \code{obo_url}
#' @inheritParams prepGO
proc_obo <- function(db_path, fn_obo, type = c("biological_process", "cellular_component", "molecular_function")) {
  filepath <- file.path(db_path, "cache", fn_obo)
  if (!file.exists(filepath)) stop("File not found ", filepath, ".", call. = FALSE)

  suppressWarnings(df <- readLines(filepath))
  first_row <- grep("\\[Term\\]", df)[1]
  last_row <- grep("\\[Typedef\\]", df)[1] - 2
  df <- df[first_row:last_row]
  
  go_ids <- df %>% .[grepl("^id:\\s{1}", .)] %>% gsub("^id:\\s{1}", "", .)
  go_nms <- df %>% .[grepl("^name:\\s{1}", .)] %>% gsub("^name:\\s{1}", "", .)
  go_type <- df %>% .[grepl("^namespace:\\s{1}", .)] %>% gsub("^namespace:\\s{1}", "", .)
  
  df <- tibble::tibble(go_id = go_ids, go_name = go_nms, go_space = go_type) %>% 
    dplyr::filter(go_space %in% type) %>% 
    dplyr::select(-go_space) 
}


#' Helper to save `gaf` without header
#' 
#' @param fn_gaf filename according to \code{gaf_url}
#' @inheritParams prepGO
proc_gaf <- function(db_path, fn_gaf) {
  filepath <- file.path(db_path, "cache", fn_gaf)
  if (!file.exists(filepath)) stop("File not found ", filepath, ".", call. = FALSE)

  con <- gzfile(path.expand(filepath))
  
  suppressWarnings(df <- readLines(con))
  first_row <- grep("UniProtKB", df)[1]
  last_row <- length(df)
  df <- df[first_row:last_row]
  
  out_nm <- gsub("\\.gz$", "_hdr_rm.txt", fn_gaf)
  writeLines(df, file.path(db_path, out_nm))
  close(con)
  
  df <- readr::read_tsv(file.path(db_path, out_nm), col_names = FALSE, 
                        col_types = cols(
                          X1 = col_character(),
                          X2 = col_character(),
                          X3 = col_character(),
                          X4 = col_character(),
                          X5 = col_character(),
                          X6 = col_character(),
                          X7 = col_character(),
                          X8 = col_character(),
                          X9 = col_character(),
                          X10 = col_character(),
                          X11 = col_character(),
                          X12 = col_character(),
                          X13 = col_character(),
                          X14 = col_character(),
                          X15 = col_character(),
                          X16 = col_character(),
                          X17 = col_character()
                        )) 
  
  unlink(file.path(db_path, out_nm))
  
  df <- df %>% 
    dplyr::select(X2:X3, X5) %>% 
    `colnames<-`(c("uniprot_acc", "gene", "go_id"))
}


#' Helper to map `SYMBOL` to `ENTREZID` 
#' 
#' @param keys Identifier such as gene names.
#' @param from the type of \code{keys}
#' @param to the type of target IDs
#' @inheritParams prepGO
#' @import dplyr purrr AnnotationDbi
annot_from_to <- function(abbr_species = "Hs", keys = NULL, from = "SYMBOL", to = "ENTREZID") {
  if (all(is.null(keys))) stop("Argument `keys` cannot be NULL", call. = FALSE)

  pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
  
  if (!requireNamespace(pkg_nm, quietly = TRUE)) {
    stop("Run `BiocManager::install(\"", pkg_nm, "\")` first.", call. = FALSE)
  }

  x <- tryCatch(
    get(pkg_nm), 
    error = function(e) 1
  )
  
  if (!is.object(x)) {
    if (x == 1) stop("Did you forget to run `library(", pkg_nm, ")`?", call. = FALSE)
  }

  accessions <- purrr::quietly(AnnotationDbi::select)(
    get(pkg_nm), 
    keys = keys,
    columns = to,
    keytype = from
  )$result
}


#' Helper to get a complete list of `ENTREZID` from `egUNIPROT`
#' 
#' @inheritParams prepGO
#' @inheritParams annot_from_to
get_full_entrez <- function(species, from = "egUNIPROT") {
  abbr_species <- sp_lookup_Ul(species)
  
  pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
  pkg_nm_from <- paste("org", abbr_species, from, sep = ".")
  
  if (!requireNamespace(pkg_nm, quietly = TRUE)) {
    stop("Run `BiocManager::install(\"", pkg_nm, "\")` first.", call. = FALSE)
  }
  
  pkg_nm_from %>% get() %>% mappedkeys()
}


#' Helper to find human orthologs
#' 
#' @inheritParams prepMSig
#' @examples \donttest{res <- find_human_orthologs(mouse)}
find_human_orthologs <- function(species, ortho_mart) {
  if (species == "human") stop("Ortholog `species` needs to be different to `human`.", call. = FALSE)
  
  out_nm <- paste0("ortho_hs", sp_lookup(species))

  data(package = "proteoQ", mart_hs)
  
  if (species == "mouse") {
    data(package = "proteoQ", mart_mm)
    martL <- mart_mm
  } else if (species == "rat") {
    data(package = "proteoQ", mart_rn)
    martL <- mart_rn
  } else {
    martL <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ortho_mart)
  }
  
  res <- biomaRt::getLDS(
    attributes = c("entrezgene_id"),
    filters = "entrezgene_id", 
    values = get_full_entrez(species = "human", from = "egUNIPROT"), 
    mart = mart_hs,
    attributesL = c("entrezgene_id"), 
    martL = martL
  ) %>% 
    `colnames<-`(c("human", species)) 
  
  assign(out_nm, res)
  save(list = out_nm, file = file.path(dat_dir, paste0(out_nm, ".rda")))

  invisible(res)
}


#' Helper to look up GO species 
#' 
#' @inheritParams prepMSig
sp_lookup_go <- function(species) {
  switch(species, 
         human = "human",
         mgi = "mouse",
         rgd = "rat",
         fb = "fly", 
         cow = "cow",
         dog = "dog", 
         crap = "crap")
}


#' Helper to find the two-letter character string of abbreviated species
#' 
#' @import stringr
#' @inheritParams prepMSig
find_abbr_species <- function(species = "human", abbr_species = NULL) {
  species <- rlang::as_string(rlang::enexpr(species))
  abbr_species <- rlang::enexpr(abbr_species)
  
  if (is.null(abbr_species) || species %in% c("human", "mouse", "rat")) {
    abbr_species <- switch(species, 
           human = "Hs",
           mouse = "Mm",
           rat = "Rn",
           stop("`species` not in one of `human`, `mouse` or `rat`.", 
                "\nThus users need to provide a two-letter abbreviation of `abbr_species`, ", 
                "\ni.e., `abbr_species = Ce` for later uses with `org.Ce.eg.db` annotation.", 
                call. = FALSE)
    )
  } else {
    abbr_species <- rlang::as_string(abbr_species)
    
    if (stringr::str_length(abbr_species) != 2) {
      stop("The number of characters needs to be `2` for `abbr_species`.", call. = FALSE)
    }
    
    if (abbr_species != stringr::str_to_title(abbr_species)) {
      stop("`abbr_species` needs to be in Title case, i.e., `Xx`.", call. = FALSE)
    }    
  }

  return(abbr_species)
}


#' Helper to set the output file name for a data base
#' 
#' @param signature A character string, i.e. "go" for uses in an output filename.
#' @inheritParams prepMSig
set_db_outname <- function(filename = NULL, species = "human", signature) {
  filename <- rlang::enexpr(filename)
  
  if (is.null(filename)) {
    abbr_species_lwr <- switch(species, 
                               human = "hs", 
                               mouse = "mm", 
                               rat = "rn", 
                               stop("`species` not in one of `human`, `mouse` or `rat`.", 
                                    "\nThus users need to provide a `filename`.", 
                                    call. = FALSE))
    
    filename <- paste0(signature, "_", abbr_species_lwr, ".rds")
  } else {
    filename <- rlang::as_string(filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    if (fn_prefix == fn_suffix) stop("No '.' to separate a basename and an extension.", call. = FALSE)
    if (fn_suffix != "rds") stop("File extension must be `.rds`.", call. = FALSE)
    filename <- paste0(fn_prefix, ".rds")
  }
}


#' Helper to download `gmt`
#' 
#' @inheritParams prepMSig
dl_msig <- function(msig_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.all.v7.0.entrez.gmt", 
                    db_path = "~\\proteoQ\\dbs\\msig", overwrite = FALSE) {
  db_path <- create_db_path(db_path)
  
  if (!grepl("\\.entrez\\.gmt$", msig_url)) {
    stop("Use a link to a `.entrez.gmt` file; not `.symbols.gmt`.", call. = FALSE)
  }
  
  fn_msig <- msig_url %>% gsub("^.*/(.*)$", "\\1", .)
  
  if ((!file.exists(file.path(db_path, "cache", fn_msig))) | overwrite)  {
    downloader::download(msig_url, file.path(db_path, "cache", fn_msig), mode = "wb")
  }
  
  return(fn_msig)
}


#' Helper to save `gmt` 
#' 
#' @param fn_gmt filename of downloaded gmt results.
#' @inheritParams prepMSig
proc_gmt <- function(species, abbr_species, ortho_mart, fn_gmt, db_path, filename) {
  filepath <- file.path(db_path, "cache", fn_gmt)
  if (!file.exists(filepath)) stop("File not found ", filepath, ".", call. = FALSE)
  
  df <- suppressWarnings(readr::read_tsv(filepath, col_names = FALSE)) %>% 
    `names_pos<-`(1, "term")
  
  df <- local({
    cols_entrez <- purrr::map_lgl(df, is.numeric) %>% which()
    
    df[, c(1, cols_entrez)] %>% 
      tidyr::gather("col", "entrez", -term) %>% 
      dplyr::select(-col)  %>% 
      dplyr::filter(!is.na(entrez)) 
  })
  
  if (species != "human") {
    df <- local({
      orthos <- find_human_orthologs(species, ortho_mart) %>% 
        `names<-`(c("from", "to"))
      
      df <- df %>% 
        dplyr::left_join(orthos, by = c("entrez" = "from")) %>% 
        dplyr::filter(!is.na(to))  %>% 
        dplyr::select(-entrez) %>% 
        dplyr::rename(entrez = to)      
    })
  }
  
  gsets <- df %>% 
    split(., .$term, drop = TRUE) %>% 
    purrr::map(`[[`, 2) %>% 
    `names<-`(paste(tolower(abbr_species), names(.), sep = "_"))
  
  saveRDS(gsets, file.path(db_path, filename))
}


#'Download and prepare gene ontology
#'
#'\code{prepGO} downloads and prepares data bases of
#'\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'ontology} (GO) for enrichment analysis by gene sets.
#'
#'@import rlang dplyr magrittr purrr fs readr downloader org.Hs.eg.db
#'  org.Mm.eg.db org.Rn.eg.db 
#'@inheritParams dl_stringdbs
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of GO. The species available for the
#'  convenience feature is in one of \code{c("human", "mouse", "rat")} with
#'  \code{"human"} being the default. The argument is not required for other
#'  species; instead, users will provide values under arguments
#'  \code{abbr_species}, \code{gaf_url} and \code{obo_url}.
#'@param abbr_species Two-letter character string; the abbreviated name of
#'  species used with
#'  \href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Xx.eg.db}.
#'   The value of \code{abbr_species} will be determined automatically if the
#'  species is in one of \code{c("human", "mouse", "rat")}. Otherwise, for
#'  example, users need to provide \code{abbr_species = Ce} for fetching the
#'  \code{org.Ce.eg.db} package in the name space of proteoQ. For
#'  \href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'  ontology} and \href{https://www.gsea-msigdb.org/gsea/index.jsp}{Molecular
#'  Signatures}, the argument is further applied to differentiate the same
#'  biological terms under different species.
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~\\proteoQ\\dbs\\go"}.
#'@param gaf_url A URL to
#'  \href{http://current.geneontology.org/products/pages/downloads.html}{GO
#'  Annotation File} (GAF). A valid web address is required for species other
#'  than \code{c("human", "mouse", "rat")}. At the NULL default and the species
#'  in one of \code{c("human", "mouse", "rat")}, the link will be determined
#'  automatically; note that users can overwrite the default GAF by providing
#'  their own URL.
#'@param obo_url A URL link to GO terms in an OBO format. At the NULL default,
#'  the web address will be determined automatically. Users can overwrite the
#'  default OBO by providing their own URL.
#'@param filename Character string; An output file name. At the NULL default,
#'  the name will be determined automatically at a given \code{species}; i.e.,
#'  \code{go_hs.rds} for human data. The file is saved as a \code{.rds} object
#'  for uses with \code{\link{prnGSPA}}.
#'@param type Character vector. The name space in gene ontology to be included.
#'  The default is to include all in \code{c("biological_process",
#'  "cellular_component", "molecular_function")}. In the example of \code{type =
#'  c("biological_process", "cellular_component")}, terms of
#'  \code{molecular_function} will be excluded.
#' @examples
#' \donttest{
#' library(proteoQ)
#'
#' # `human` and `mouse` with a default OBO;
#' # outputs under `db_path`
#' prepGO(human)
#' prepGO(mouse)
#'
#' # `mouse` with a slim OBO
#' prepGO(
#'   species = mouse,
#'   obo_url = "http://current.geneontology.org/ontology/subsets/goslim_mouse.obo",
#'   filename = mm_slim.rds,
#' )
#'
#' # `worm` not available for default GO preparation
#' if (!requireNamespace("BiocManager", quietly = TRUE))
#'     install.packages("BiocManager")
#' BiocManager::install("org.Ce.eg.db")
#'
#' library(org.Ce.eg.db)
#'
#' prepGO(
#'   # species = worm,
#'   abbr_species = Ce,
#'   gaf_url = "http://current.geneontology.org/annotations/wb.gaf.gz",
#'   obo_url = "http://purl.obolibrary.org/obo/go/go-basic.obo",
#' )
#' }
#' 
#' \dontrun{
#' gsets <- readRDS(file.path("~\\proteoQ\\dbs\\go", "mm_slim.rds"))
#' }
#'@export
prepGO <- function(species = "human", abbr_species = NULL, gaf_url = NULL, obo_url = NULL, 
                   db_path = "~\\proteoQ\\dbs\\go", 
                   type = c("biological_process", "cellular_component", "molecular_function"), 
                   filename = NULL, overwrite = FALSE) {
  
  old_opt <- options(warn = 0)
  options(warn = 1)
  on.exit(options(old_opt), add = TRUE)
  
  species <- rlang::as_string(rlang::enexpr(species))

  db_path <- create_db_path(db_path)

  if (is.null(gaf_url)) {
    gaf_url <- switch(species, 
                      human = c("goa_human.gaf.gz" = "http://current.geneontology.org/annotations/goa_human.gaf.gz"), 
                      mouse = c("goa_mouse.gaf.gz" = "http://current.geneontology.org/annotations/mgi.gaf.gz"), 
                      rat = c("goa_rat.gaf.gz" = "http://current.geneontology.org/annotations/rgd.gaf.gz"), 
                      stop("`species` need to be one of `human`, `mouse` or `rat` for an auto lookup of GAF files.", 
                           call. = FALSE))
    fn_gaf <- names(gaf_url)
  } else {
    fn_gaf <- gaf_url %>% gsub("^.*/(.*)$", "\\1", .)
    
    species <- local({
      species_go <- fn_gaf %>% 
        gsub("([^\\.]*)[\\.].*", "\\1", .) %>% 
        gsub("^goa_", "", .) %>% 
        sp_lookup_go()
      
      if (is.null(species_go)) {
        species_go <- "unknown"
      } else if (species_go != species) {
        message("The species is `", species_go, "`.")
      }
      
      species <- species_go      
    })
  }
  
  if (is.null(obo_url)) {
    obo_url <- c("go-basic.obo" = "http://purl.obolibrary.org/obo/go/go-basic.obo")
    fn_obo <- names(obo_url)
  } else {
    fn_obo <- obo_url %>% gsub("^.*/(.*)$", "\\1", .)
  }

  if ((!file.exists(file.path(db_path, "cache", fn_gaf))) | overwrite)  {
    downloader::download(gaf_url, file.path(db_path, "cache", fn_gaf), mode = "wb")
  }
  
  if ((!file.exists(file.path(db_path, "cache", fn_obo))) | overwrite) {
    downloader::download(obo_url, file.path(db_path, "cache", fn_obo), mode = "wb")
  }
  
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, "go")

  df <- local({
    df_gaf <- proc_gaf(db_path, fn_gaf)
    df_obo <- proc_obo(db_path, fn_obo, type)
    
    dplyr::left_join(df_gaf, df_obo, by = "go_id") %>% 
      dplyr::mutate(go_name = paste(go_id, go_name)) %>% 
      dplyr::select(go_name, gene)
  })
  
  accessions <- annot_from_to(abbr_species = abbr_species, 
                              keys = unique(df$gene), from = "SYMBOL", to = "ENTREZID") %>% 
    dplyr::rename(gene = SYMBOL, entrez = ENTREZID) %>% 
    dplyr::filter(!is.na(entrez), !is.na(gene))  %>% 
    dplyr::filter(!duplicated(gene))
  
  gsets <- accessions %>% 
    dplyr::left_join(df, by = "gene") %>% 
    dplyr::select(-gene) %>% 
    split(., .$go_name, drop = TRUE) %>% 
    purrr::map(`[[`, 1) %>% 
    `names<-`(paste(tolower(abbr_species), names(.), sep = "_"))
  
  saveRDS(gsets, file.path(db_path, filename))
}


#'Preparation of Molecular Signatures
#'
#'\code{prepMSig} downloads and prepares data bases of
#'\href{https://www.gsea-msigdb.org/gsea/index.jsp}{Molecular Signatures}
#'(MSig) for enrichment analysis by gene sets.
#'
#'@import rlang dplyr magrittr purrr fs readr downloader biomaRt org.Hs.eg.db
#'@inheritParams dl_stringdbs
#'@inheritParams prepGO
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of \code{MSig} data bases. The species
#'  available for the convenience feature is in one of \code{c("human", "mouse",
#'  "rat")} with \code{"human"} being the default. The argument is not required
#'  for other species; instead, users will provide values under arguments
#'  \code{ortho_mart} for the lookup of orthologs to human.
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~\\proteoQ\\dbs\\msig"}.
#'@param msig_url A URL to
#'  \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{MSig
#'  }. At the \code{NULL} default, a \code{c2.all.v[...].entrez.gmt} data will
#'  be used for all species. A valid web address is required for a custom data
#'  base, and for simplicity, only files with entrez IDs will be handled.
#'@param ortho_mart Character string; a dataset name from
#'  \code{\link[biomaRt]{useMart}} and/or \code{\link[biomaRt]{listDatasets}}
#'  for the lookup of orthologs to \code{human} genes. For species in
#'  \code{c("human", "mouse", "rat")}, the value will be determined
#'  automatically unless otherwise specified.
#'@param filename Character string; An output file name. At the \code{NULL}
#'  default, the name will be determined automatically at a given
#'  \code{species}; i.e., \code{msig_hs.rds} for \code{human} data. The file is
#'  saved as a \code{.rds} object for uses with \code{\link{prnGSPA}}.
#' @examples
#' \donttest{
#' library(proteoQ)
#'
#' ## the default `MSig` is `c2.all`
#' # `human`; outputs under `db_path`
#' prepMSig()
#' head(readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_hs.rds")))
#'
#' # `mouse`
#' prepMSig(species = mouse, filename = msig_mm.rds)
#' head(readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_mm.rds")))
#'
#' # `rat`
#' prepMSig(species = rat, filename = msig_rn.rds)
#' head(readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_rn.rds")))
#'
#' # `dog`; need `ortho_mart` for species other than `human`, `mouse` and `rat`
#' # (try `?biomaRt::useMart` for a list of marts)
#' prepMSig(
#'   # species = dog,
#'   abbr_species = Cf, 
#'   ortho_mart = cfamiliaris_gene_ensembl,
#'   filename = msig_cf.rds,
#' )
#'
#' msig_cf <- readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_cf.rds"))
#'
#' # also `dog`
#' prepMSig(
#'   species = my_dog,
#'   abbr_species = Cf, 
#'   ortho_mart = cfamiliaris_gene_ensembl,
#'   filename = msig_cf2.rds,
#' )
#'
#' msig_cf2 <- readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_cf2.rds"))
#'
#'identical(msig_cf, msig_cf2)
#'
#' ## use an `MSig`other than the default of `c2.all`
#' prepMSig(
#'   msig_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.cgp.v7.0.entrez.gmt",
#'   species = human,
#'   filename = c2_cgp_hs.rds,
#' )
#'
#' prepMSig(
#'   msig_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.cgp.v7.0.entrez.gmt",
#'   species = dog,
#'   ortho_mart = cfamiliaris_gene_ensembl,
#'   filename = c2_cgp_cf.rds,
#' )
#' }
#'
#'@export
prepMSig <- function(msig_url = NULL, 
                     species = "human", 
                     abbr_species = NULL, 
                     ortho_mart = switch(species, 
                                         mouse = "mmusculus_gene_ensembl", 
                                         rat = "rnorvegicus_gene_ensembl", 
                                         human = "zzz",
                                         "unknown"), 
                     db_path = "~\\proteoQ\\dbs\\msig", filename = NULL, overwrite = FALSE) {

  old_opt <- options(warn = 0)
  options(warn = 1)
  on.exit(options(old_opt), add = TRUE)
  
  if (is.null(msig_url)) {
    msig_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.all.v7.0.entrez.gmt"
  } else {
    msig_url <- rlang::as_string(rlang::enexpr(msig_url))
  }
  
  fn_gmt <- dl_msig(msig_url, db_path, overwrite)
  
  species <- rlang::as_string(rlang::enexpr(species))
  
  ortho_mart <- local({
    ok <- tryCatch(force(ortho_mart), error = function(e) 1)
    
    if (ok == 1) {
      ortho_mart <- rlang::as_string(rlang::enexpr(ortho_mart))
    } else {
      ortho_mart <- ok
    }
  })

  if (ortho_mart == "unknown") stop("Specify the value of `ortho_mart` for species `", species, "`.", 
                                    call. = FALSE)
  
  if (species == "human" && ortho_mart != "zzz") species <- "unknown"
  
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, "msig")
  proc_gmt(species, abbr_species, ortho_mart, fn_gmt, db_path, filename)
}


#' Map UniProt or Refseq accessions to Entrez IDs
#'
#' @param from Character string; the type of accession keys in c("UNIPROT", "REFSEQ").
#' @inheritParams prepMSig
#' @inheritParams annot_from_to
#' @import dplyr purrr tidyr plyr reshape2 org.Hs.eg.db org.Mm.eg.db
#'   org.Rn.eg.db
#' @export
map_to_entrez <- function(species = "human", abbr_species = NULL, from = "UNIPROT", 
                          filename = NULL, db_path = "~\\proteoQ\\dbs\\entrez", overwrite = FALSE) {
  db_path <- create_db_path(db_path)
  species <- rlang::as_string(rlang::enexpr(species))
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  from <- rlang::as_string(rlang::enexpr(from))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, paste(tolower(from), "entrez", sep = "_" ))
  
  if ((!file.exists(file.path(db_path, filename))) || overwrite)  {
    pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
    if (!requireNamespace(pkg_nm, quietly = TRUE)) {
      stop("Run `BiocManager::install(\"", pkg_nm, "\")` first, 
           then `library(", pkg_nm, ")`", call. = FALSE)
    }
    
    new_from <- paste0("eg", from)
    x <- tryCatch(
      get(paste("org", abbr_species, new_from, sep = ".")),
      error = function(e) 1
    )
    
    if (!is.object(x)) {
      if (x == 1) stop("Did you forget to run `library(", pkg_nm, ")`?", call. = FALSE)
    }
    
    entrez_ids <- mappedkeys(x) 
    
    accessions <- as.list(x[entrez_ids]) %>% 
      plyr::ldply(., rbind) %>% 
      `names_pos<-`(., 1, c("entrez")) %>% 
      `names_pos<-`(., 2:ncol(.), paste(new_from, 1:(length(.)-1), sep = ".")) %>% 
      dplyr::mutate_at(.vars = grep("^eg", names(.)), ~ as.character(.x)) %>% 
      tidyr::gather("variable", "value", -entrez) %>% 
      dplyr::filter(!is.na(entrez), !is.na(value)) %>% 
      dplyr::select(-c("variable")) %>% 
      dplyr::mutate(species = species)
    
    if (from == "UNIPROT") {
      accessions <- accessions %>% dplyr::rename(uniprot_acc = value)
    } else if (from == "REFSEQ") {
      accessions <- accessions %>% dplyr::rename(refseq_acc = value)
      
      accessions <- annot_from_to(abbr_species, unique(accessions$refseq_acc), 
                                  "REFSEQ", "SYMBOL") %>% 
        dplyr::rename(refseq_acc = REFSEQ, gene = SYMBOL) %>% 
        dplyr::filter(!is.na(refseq_acc), !is.na(gene))  %>% 
        # dplyr::filter(!duplicated(gene)) %>% 
        dplyr::left_join(accessions, by = "refseq_acc")
    } else {
      stop("Variable `from` needs to be either `UNIPROT` or `REFSEQ`.", call. = FALSE)
    }
    
    saveRDS(accessions, file.path(db_path, filename))

    invisible(accessions)
  } else {
    invisible(NULL)
  }
}


#'Map UniProt accessions to Entrez IDs
#'
#'\code{Uni2Entrez} prepares lookup tables between UniProt accessions and
#'Entrez IDs for uses with \link{normPSM} and downstream gene-set analysis such
#'as \link{prnGSPA}. The utility is optional for \code{human}, \code{mouse} and
#'\code{rat} data. It is \strong{required} for other species with \link{prnGSPA}
#'in users' workflows. It can also be used to update and overrule the lookups
#'for \code{human}, \code{mouse} and \code{rat} that are defaulted by
#'\code{proteoQ}.
#'
#'@param species Character string; the name of a species. 
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~\\proteoQ\\dbs\\entrez"}.
#'@inheritParams prepMSig
#'@import dplyr purrr tidyr plyr reshape2 org.Hs.eg.db org.Mm.eg.db org.Rn.eg.db
#'@example inst/extdata/examples/prepEntrez_.R
#'@export
Uni2Entrez <- function(species = "human", abbr_species = NULL, filename = NULL, 
                       db_path = "~\\proteoQ\\dbs\\entrez", overwrite = FALSE) {
  map_to_entrez(!!rlang::enexpr(species), 
                !!rlang::enexpr(abbr_species), 
                "UNIPROT", 
                !!rlang::enexpr(filename), 
                db_path, 
                overwrite)
}




#'Map RefSeq accessions to Entrez IDs and gene names
#'
#'\code{Ref2Entrez} prepares lookup tables between RefSeq accessions and
#'Entrez IDs and gene names for uses with \link{normPSM} and downstream gene-set
#'analysis such as \link{prnGSPA}. The utility is optional for \code{human} and
#'\code{mouse} data. It is \strong{required} for other species with
#'\link{prnGSPA} in users' workflows. It can also be used to update and overrule
#'the lookups for \code{human} and \code{mouse} that are defaulted by
#'\code{proteoQ}.
#'
#'@rdname Uni2Entrez
#'@import dplyr purrr tidyr plyr reshape2 org.Hs.eg.db org.Mm.eg.db org.Rn.eg.db
#'@example inst/extdata/examples/prepEntrez_.R
#'@export
Ref2Entrez <- function(species = "human", abbr_species = NULL, filename = NULL, 
                       db_path = "~\\proteoQ\\dbs\\entrez", overwrite = FALSE) {
  map_to_entrez(!!rlang::enexpr(species), 
                !!rlang::enexpr(abbr_species), 
                "REFSEQ", 
                !!rlang::enexpr(filename), 
                db_path, 
                overwrite)
}




#' Map uniprot or refseq to entrez (not currently used)
#'
#' @param os_name An organism name by UniProt.
#' @inheritParams Uni2Entrez
#' @inheritParams annot_from_to
#' @import dplyr purrr tidyr plyr reshape2 org.Hs.eg.db org.Mm.eg.db
#'   org.Rn.eg.db
map_to_entrez_os_name <- function(species = "human", abbr_species = NULL, 
                                 os_name = "Homo sapiens", 
                                 from = "UNIPROT", 
                                 filename = NULL, db_path = "~\\proteoQ\\dbs\\entrez", overwrite = FALSE) {
  # the value of species will be used under column `species` in Protein.txt etc
  # add column species in the rds output
  
  db_path <- create_db_path(db_path)
  species <- rlang::as_string(rlang::enexpr(species))
  os_name <- rlang::as_string(rlang::enexpr(os_name))
  
  create_os_lookup(species, os_name, overwrite)
  
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, paste(tolower(from), "entrez", sep = "_" ))
  
  if ((!file.exists(file.path(db_path, filename))) || overwrite)  {
    pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
    if (!requireNamespace(pkg_nm, quietly = TRUE)) {
      stop("Run `BiocManager::install(\"", pkg_nm, "\")` first, 
           then `library(", pkg_nm, ")`", call. = FALSE)
    }
    
    new_from <- paste0("eg", from)
    x <- tryCatch(
      get(paste("org", abbr_species, new_from, sep = ".")),
      error = function(e) 1
    )
    
    if (!is.object(x)) {
      if (x == 1) stop("Did you forget to run `library(", pkg_nm, ")`?", call. = FALSE)
    }
    
    entrez_ids <- mappedkeys(x) 
    
    accessions <- as.list(x[entrez_ids]) %>% 
      plyr::ldply(., rbind) %>% 
      `names_pos<-`(., 1, c("entrez")) %>% 
      `names_pos<-`(., 2:ncol(.), paste(new_from, 1:(length(.)-1), sep = ".")) %>% 
      dplyr::mutate_at(.vars = grep("^eg", names(.)), ~ as.character(.x)) %>% 
      tidyr::gather("variable", "value", -entrez) %>% 
      dplyr::filter(!is.na(entrez), !is.na(value)) %>% 
      dplyr::select(-c("variable")) # %>% 
    # dplyr::mutate(species = species)
    
    if (from == "UNIPROT") {
      accessions <- accessions %>% dplyr::rename(uniprot_acc = value)
    } else if (from == "REFSEQ") {
      accessions <- accessions %>% dplyr::rename(refseq_acc = value)
    } else {
      stop("Variable `from` needs to be either `UNIPROT` or `REFSEQ`.", call. = FALSE)
    }
    
    saveRDS(accessions, file.path(db_path, filename))
    
    invisible(accessions)
  } else {
    invisible(NULL)
  }
}



#' create a lookup table and remove duplicated entries
#' 
#' @inheritParams map_to_entrez_os_name
create_os_lookup <- function(species, os_name, overwrite = FALSE) {
  my_lookup <- c(
    "Homo sapiens" = "human",
    "Mus musculus" = "mouse",
    "Rattus norvegicus" = "rat"
  )

  # check if conflict between species and os_name
  convert_default_os <- function (species, os_name) {
    if (species %in% my_lookup) {
      lookup_nm <- my_lookup %>% .[. %in% species] %>% names()
      if (!purrr::is_empty(lookup_nm) && lookup_nm != os_name) {
        os_name <- lookup_nm
      }
    }
    
    return(os_name)
  }
  
  os_name <- convert_default_os(species, os_name)
  curr_lookup <- uniprot_entrez_lookup <- setNames(species, os_name)

  if (os_name == "Homo sapiens" && species != "human") 
    stop("`os_name = Homo sapiens` is reversed for `species = human`.", call. = FALSE)
  
  if (os_name == "Mus musculus" && species != "mouse") 
    stop("`os_name = Mus musculus` is reversed for `species = mouse`.", call. = FALSE)
  
  if (os_name == "Rattus norvegicus" && species != "rat") 
    stop("`os_name = Rattus norvegicus` is reversed for `species = rat`.", call. = FALSE)  
  
  # check if the same os_name but diferent species
  file <- file.path(dat_dir, "uniprot_entrez_lookup.rda")
  
  if (file.exists(file)) {
    load(file)
    
    if (overwrite) {
      uniprot_entrez_lookup <- uniprot_entrez_lookup %>% 
        .[! names(.) == os_name] %>% 
        .[! . == species]
    } else {
      old_sp <- uniprot_entrez_lookup %>% .[names(.) == os_name]
      curr_sp <- curr_lookup

      if (!purrr::is_empty(old_sp) && (curr_sp != old_sp)) {
        stop("`", names(curr_lookup), "` was previously linked to `species = ", old_sp, "`.", 
             "\nTo overwrite the value of `os_name`, set `overwrite = TRUE`.", 
             call. = FALSE)
      }
      
      old_nm <- uniprot_entrez_lookup %>% .[. == species] %>% names()
      curr_nm <- os_name

      if (!purrr::is_empty(old_nm) && (curr_nm != old_nm)  && !overwrite) {
        stop("`", species, "` was previously linked to `os_name = ", old_nm, "`.", 
             "\nTo overwrite, set `overwrite = TRUE`.", 
             call. = FALSE)
      }
    }

    uniprot_entrez_lookup <- c(curr_lookup, uniprot_entrez_lookup) %>% 
      .[!duplicated(.)]
  }
  
  save(uniprot_entrez_lookup, file = file)
}


#' Make GSEA gct
#' 
#' @param fn_prefix The base name of a file.
#' @inheritParams info_anal
make_gct <- function(df, filepath, fn_prefix) {
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  
  df <- df %>% 
    tibble::rownames_to_column("NAME") %>% 
    dplyr::mutate(Description = NA) 
  
  df <- bind_cols(
    df %>% dplyr::select(NAME, Description), 
    df %>% dplyr::select(-NAME, -Description), 
  )
  
  hr_l1 <- "#1.2\n"
  hr_l2 <- paste0(c(nrow(df), ncol(df) - 2), collapse = "\t") %>% 
    paste0("\n")
  
  hr_nms <- paste0(names(df), collapse = "\t") %>% 
    paste0("\n")
  
  header <- paste0(hr_l1, hr_l2, hr_nms)
  file <- file.path(filepath, paste0(fn_prefix, ".gct"))
  cat(header, file = file)
  df %>% write.table(file, append = TRUE, sep = '\t', na = "", 
                     col.names = FALSE, row.names = FALSE, quote = FALSE)
}


#'Make GSEA cls 
#' 
#' @param nms The names of sample groups.
#' @inheritParams info_anal
#' @inheritParams make_gct
make_cls <- function(df, nms, filepath, fn_prefix) {
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  
  cls_l1 <- paste0(c(length(nms), length(unique(nms)), "1"), collapse = "\t") %>% 
    paste0("\n")
  
  grps <- nms %>% 
    as.character() %>% 
    unique() %>% 
    paste0(collapse = "\t") %>% 
    paste0("#", ., "\n")
  
  smpls <- nms %>% 
    as.character() %>% 
    paste0(collapse = "\t") %>% 
    paste0("\n")
  
  cls <- paste0(cls_l1, grps, smpls)
  fn_cls <- file.path(filepath, paste0(fn_prefix, ".cls"))
  cat(cls, file = fn_cls)  
}


#'Protein GSEA
#'
#'\code{prnGSEA} prepares data for the analysis of
#'\href{http://software.broadinstitute.org/gsea/index.jsp}{GSEA} against
#'protein \code{log2FC} data. 
#'
#'The arguments \code{var_cutoff}, \code{pval_cutoff} and \code{logFC_cutoff}
#'are used to filter out low influence genes. Additional subsetting of data via
#'the \code{vararg} approach of \code{filter_} is feasible.
#'
#'The outputs include \code{Protein_GSEA.gct} and \code{protein_GSEA.cls} for
#'samples indicated in file \code{Protein_pVals.txt} or
#'\code{Protein_impNA_pVals.txt}. These outputs can be used with online
#'\href{http://software.broadinstitute.org/gsea/index.jsp}{GSEA}.
#'
#'The current GSEA may not support the comparisons between two grouped
#'conditions, i.e.,  (grpA + grpB) versus (grpC + grpD). The \code{prnGSEA}
#'utility further breaks the input data into pairs of groups according to the
#'formulas and contrasts defined in \code{pepSig} or \code{prnSig}. The
#'phenotype labels are then reformed in reflection of the original group names,
#'weights and directions, i.e., \code{0.5xgrpA&0.5xgrpB	-0.5xgrpC&-0.5xgrpD}.
#'The corresponding \code{.gct} and \code{.cls} files can be used with the
#'online or the github version of R-GSEA.
#'
#'@inheritParams prnGSPA
#'@param gset_nms Not currently used (to be chosen by users during online GSEA).
#'@param var_cutoff Numeric value or vector; the cut-off in the variance of
#'  protein \code{log2FC} across samples. Entries with \code{variances} less
#'  than the threshold will be removed for enrichment analysis. The default is
#'  0.5.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. 
#'@import dplyr rlang ggplot2 networkD3
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnGSEA_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@export
prnGSEA <- function (gset_nms = "go_sets", 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     df = NULL, filepath = NULL, filename = NULL, 
                     var_cutoff = .5, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     fml_nms = NULL, ...) {

  on.exit(
    if (rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(dots) %>% 
        save_call("pepGSEA")
    } else if (rlang::as_string(id) %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(dots) %>% 
        save_call("prnGSEA")
    }
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\GSEA\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSEA")
  
  reload_expts()
  
  # Sample selection criteria:
  #   !is_reference under "Reference"
  #   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
  info_anal(df = !!df, df2 = NULL, id = !!id, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            filepath = !!filepath, filename = !!filename, 
            anal_type = "GSEA")(gset_nms = gset_nms, 
                                var_cutoff = var_cutoff, 
                                pval_cutoff = pval_cutoff, 
                                logFC_cutoff = logFC_cutoff, 
                                !!!dots)
}


#'Protein GSEA by formula(s) in `pepSig` or `prnSig`
#'
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams fml_gspa
#' @inheritParams gsVolcano
fml_gsea <- function (fml, fml_nm, var_cutoff, pval_cutoff, 
                      logFC_cutoff, gspval_cutoff, gslogFC_cutoff, min_size, max_size, 
                      df, col_ind, id, gsets, label_scheme_sub, complete_cases, scale_log2r, 
                      filepath, filename, ...) {
  
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  df <- local({
    ids <- df %>% 
      prep_gspa(id = !!id, fml_nm = fml_nm, col_ind = col_ind, 
                pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff) %>% 
      dplyr::select(id) %>% 
      unlist() %>% 
      unique()
    
    df <- df %>% dplyr::filter(!!sym(id) %in% ids)
    stopifnot(nrow(df) > 0)

    ids <- df %>% 
      `rownames<-`(.[[id]]) %>% 
      filterData(cols = grep(NorZ_ratios, names(.)), var_cutoff = var_cutoff) %>% 
      tibble::rownames_to_column(id) %>% 
      dplyr::select(id) %>% 
      unlist() %>% 
      unique()
    
    df <- df %>% dplyr::filter(!!sym(id) %in% ids)
    
    # No `lm` so `complete_cases` applied to all samples in label_scheme_sub
    if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

    return(df)
  })
  
  if (id != "gene") {
    df <- df %>% gn_rollup(grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)))
    id <- "gene"
  }
  
  df_log2r <- df %>% 
    prepDM(id = !!id, 
           scale_log2r = scale_log2r, 
           sub_grp = label_scheme_sub$Sample_ID, 
           anal_type = "GSEA") %>% 
    .$log2R
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  make_gct(df = df_log2r, filepath = filepath, fn_prefix = fn_prefix)
  make_cls(df = df_log2r, nms = label_scheme_sub$Group, filepath = filepath, fn_prefix = fn_prefix)
  
  fml_ops <- suppressMessages(prepFml(fml, label_scheme_sub, ...))
  contr_mat <- fml_ops$contr_mat
  design <- fml_ops$design
  key_col <- fml_ops$key_col
  random_vars <- fml_ops$random_vars
  label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub
  
  df <- df %>% 
    `rownames<-`(NULL) %>% 
    tibble::column_to_rownames(id) %>% 
    `names<-`(gsub(paste0(NorZ_ratios, "[0-9]{3}[NC]*\\s+\\((.*)\\)$"), "\\1", names(.))) %>% 
    dplyr::select(as.character(label_scheme_sub_sub$Sample_ID))

  nms <- purrr::map(1:ncol(contr_mat), ~ {
    rows <- contr_mat[, .x, drop = FALSE]
    
    rows_pos <- rows[rows > 0, , drop = FALSE] 
    rows_pos[, 1] <- rows_pos[, 1] %>% purrr::map_dbl(~ round(.x, digits = 2))
    nms_pos <- paste0(rows_pos, "x", rownames(rows_pos)) %>% 
      reduce(paste, sep = "&")
    
    rows_neg <- rows[rows < 0, , drop = FALSE]
    rows_neg[, 1] <- rows_neg[, 1] %>% purrr::map_dbl(~ round(.x, digits = 2))
    nms_neg <- paste0(rows_neg, "x", rownames(rows_neg)) %>% 
      reduce(paste, sep = "&")
    
    data.frame(pos = nms_pos, neg = nms_neg)
  }) %>% 
    do.call(rbind, .) %>% 
    `rownames<-`(colnames(contr_mat))
  
  nms <- nms %>% 
    purrr::map(~ gsub("[\\\\\\/\\:\\*\\?\\'\\<\\>\\|]", ".", .x)) %>% 
    bind_rows()

  contr_identity <- ifelse(contr_mat == 0, 0, 1)
  
  n_col <- ncol(design)
  df_sub <- rep(list(NULL), n_col)
  
  n_col2 <- ncol(contr_identity)
  pos <- neg <- both <- rep(list(NULL), n_col2)
  
  for (i in 1:n_col) df_sub[[i]] <- df[, which(design[, i] == 1), drop = FALSE]
  
  for (i in 1:n_col2) {
    filepath_i <- file.path(filepath, fml_nm)

    both[[i]] <- map2(df_sub, contr_identity[, i], ~ .x * .y) %>% 
      do.call(cbind, .) %>% 
      .[, names(df), drop = FALSE]
    
    idx_pos <- design[, contr_mat[, i] > 0, drop = FALSE] %>% 
      rowSums()
    pos[[i]] <- both[[i]][, which(idx_pos == 1), drop = FALSE]
    
    idx_neg <- design[, contr_mat[, i] < 0, drop = FALSE] %>% 
      rowSums()
    neg[[i]] <- both[[i]][, which(idx_neg == 1), drop = FALSE]
    
    df_i <- cbind(pos[[i]], neg[[i]])
    out_nm <- paste0("fml-", fml_nm, "_contr-", i)
    make_gct(df_i, filepath_i, paste0(fn_prefix, "_", out_nm))

    nms_pos_i <- rep(nms$pos[i], ncol(pos[[i]])) %>% as.character()
    nms_neg_i <- rep(nms$neg[i], ncol(neg[[i]])) %>% as.character()
    nms_both <- c(nms_pos_i, nms_neg_i)
    make_cls(df = df_i, nms = nms_both, filepath = filepath_i, fn_prefix = paste0(fn_prefix, "_", out_nm))
  }    
}


#'GSPA of protein data
#'
#'\code{prnGSPA} performs the analysis of Gene Set Probability Asymmetricity
#'(GSPA) against protein \code{log2FC} data.
#'
#'The significance \code{pVals} of individual proteins are first obtained from
#'\code{\link{prnSig}}, followed by log10 transformation and separation into up-
#'or down-expressed groups for each gene set. At the default of \code{method =
#'mean}, the geometric mean of \code{pVals}, \eqn{P}, are each calculated for
#'the groups of up or down regulated proteins, with a penalty-like term
#'
#'\deqn{-log10(P)=(\sum_{i=1}^{n}-log10(p)+m)/(n+m)}{-log10(P)=(-log10(p_1*p_2*...*p_n)+m)/(n+m)}
#'
#'where \eqn{n} and \eqn{m} are the numbers of entries with \eqn{p} values
#'\eqn{\le} or \eqn{>} a significance cut-off, respectively. The quotient of the
#'two \eqn{P} values is then used to represent the significance of gene set
#'enrichment. The arguments \code{pval_cutoff} and \code{logFC_cutoff} are used
#'to discriminate low influence genes. Additional subsetting of data via the
#'\code{vararg} approach of \code{filter_} is feasible. At \code{method =
#'limma}, moderated t-tests are performed between the up and the down groups via
#'\code{\link[limma]{eBayes}}.
#'
#'@inheritParams anal_pepNMF
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param gset_nms Character string or vector containing the shorthanded name(s),
#'  full file path(s) or both to gene sets for enrichment analysis. For species
#'  among \code{"human", "mouse", "rat"}, the default of \code{c("go_sets",
#'  "c2_msig")} will utilize terms from both gene ontology (\code{GO}) and
#'  molecular signatures (\code{MSig}). Custom data bases of \code{GO} and
#'  curated \code{MSig}, and/or additional species are also supported. See also
#'  \code{\link{prepGO}} for the preparation of custom \code{GO} and
#'  \code{\link{prepMSig}} for the preparation of custom \code{MSig}.
#'@param method Character string; the method to assess the p-values of GSPA. The
#'  default is \code{mean}. See also section \code{Details} for the
#'  calculations.
#'@param pval_cutoff Numeric value or vector; the cut-off in protein
#'  significance \code{pVal}. Entries with \code{pVals} less significant than
#'  the threshold will be excluded from enrichment analysis. The default is 0.05
#'  for all formulas matched to or specified in argument \code{fml_nms}.
#'  Formula-specific threshold is allowed by supplying a vector of cut-off
#'  values.
#'@param logFC_cutoff Numeric value or vector; the cut-off in protein
#'  \code{log2FC}. Entries with absolute \code{log2FC} smaller than the
#'  threshold will be excluded from enrichment analysis. The default magnitude
#'  is \code{log2(1.2)} for all formulas matched to or specified in argument
#'  \code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
#'  of absolute values in \code{log2FC}.
#'@param gspval_cutoff Numeric value or vector; the cut-off in gene-set
#'  significance \code{pVal}. Only enrichment terms with \code{pVals} more
#'  significant than the threshold will be reported. The default is 0.05 for all
#'  formulas matched to or specified in argument \code{fml_nms}.
#'  Formula-specific threshold is allowed by supplying a vector of cut-off
#'  values.
#'@param gslogFC_cutoff Numeric value or vector; the cut-off in gene-set
#'  enrichment fold change. Only enrichment terms with absolute fold change
#'  greater than the threshold will be reported. The default magnitude is
#'  \code{log2(1.2)} for all formulas matched to or specified in argument
#'  \code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
#'  of absolute values in \code{log2FC}.
#'@param min_size Numeric value or vector; minimum number of protein entries for
#'  consideration in gene set tests. The number is after data filtration by
#'  \code{pval_cutoff}, \code{logFC_cutoff} or varargs expressions under
#'  \code{filter_}. The default is 10 for all formulas matched to or specified
#'  in argument \code{fml_nms}. Formula-specific threshold is allowed by
#'  supplying a vector of sizes.
#'@param max_size Numeric value or vector; maximum number of protein entries for
#'  consideration in gene set tests. The number is after data filtration by
#'  \code{pval_cutoff}, \code{logFC_cutoff} or varargs expressions under
#'  \code{filter_}. The default in infinite for all formulas matched to or
#'  specified in argument \code{fml_nms}. Formula-specific threshold is allowed
#'  by supplying a vector of sizes.
#'@param min_delta Numeric value or vector; the minimum count difference between
#'  the up- and the down-expressed group of proteins for consideration in gene
#'  set tests. For example at \code{min_delta = 4}, a gene set will 6
#'  upregulated proteins and 2 down-expressed proteins, or vice versa, will be
#'  assessed. The number is after data filtration by \code{pval_cutoff},
#'  \code{logFC_cutoff} or varargs expressions under \code{filter_}. The default
#'  is 4 for all formulas matched to or specified in argument \code{fml_nms}.
#'  Formula-specific threshold is allowed by supplying a vector of sizes.
#'@param min_greedy_size Numeric value or vector; minimum number of unique
#'  protein entries for a gene set to be considered essential. The default in
#'  \code{1} for all formulas matched to or specified in argument
#'  \code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
#'  of sizes.
#'@param fml_nms Character string or vector; the formula name(s). By default,
#'  the formula(s) will match those used in \code{\link{pepSig}} or
#'  \code{\link{prnSig}}.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration
#'  against data in a primary file of \code{\\Model\\Protein[_impNA]_pVals.txt}.
#'  See also \code{\link{normPSM}} for the format of \code{filter_} statements.
#'  \cr \cr \code{arrange_}: Variable argument statements for the row ordering
#'  against data in a primary file of \code{\\Model\\Protein[_impNA]_pVals.txt}.
#'  See also \code{\link{prnHM}} for the format of \code{arrange_} statements.
#'@import dplyr rlang ggplot2 networkD3
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnGSPA_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
prnGSPA <- function (gset_nms = c("go_sets", "c2_msig"), method = c("mean","limma"), 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), 
                     min_size = 10, max_size = Inf, 
                     min_delta = 4, min_greedy_size = 1, 
                     fml_nms = NULL, df = NULL, filepath = NULL, filename = NULL, ...) {
  
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_pepGSPA"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_prnGSPA"))
    }
    , add = TRUE
  )  
  
  check_dots(c("id", "anal_type", "var_cutoff"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\GSPA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  method <- rlang::enexpr(method)
  if (method == rlang::expr(c("mean","limma"))) {
    method <- "mean"
  } else {
    method <- rlang::as_string(method)
    stopifnot(method %in% c("limma", "mean"))
  }
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSPA")
  
  reload_expts()
  
  info_anal(df = !!df, id = !!id, filepath = !!filepath, filename = !!filename, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            anal_type = "GSPA")(gset_nms = gset_nms, 
                                var_cutoff = 1000, 
                                pval_cutoff = pval_cutoff, 
                                logFC_cutoff = logFC_cutoff, 
                                gspval_cutoff = gspval_cutoff, 
                                gslogFC_cutoff = gslogFC_cutoff,
                                min_size = min_size, 
                                max_size = max_size, 
                                min_delta = min_delta, 
                                min_greedy_size = min_greedy_size, 
                                method = method, 
                                !!!dots)
}


#' Perform GSPA tests
#'
#' logFC_cutoff Numeric A threshold for the subset of data before the calculation
#' of adjusted pvals
#'
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams prnGSPA
#' @inheritParams prnGSVA
#' @inheritParams info_anal
#' @param id Currently only "entrez".
#' @param label_scheme_sub A data frame. Subset entries from \code{label_scheme}
#'   for selected samples.
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$% %T>%
#' @importFrom outliers grubbs.test
#' @importFrom broom.mixed tidy
gspaTest <- function(df = NULL, id = "entrez", label_scheme_sub = NULL, 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                     filepath = NULL, filename = NULL, 
                     gset_nms = "go_sets", var_cutoff = 0.5, 
                     pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                     gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1), 
                     min_size = 6, max_size = Inf, min_delta = 4, min_greedy_size = 1, 
                     method = "mean", anal_type = "GSPA", ...) {

  # `GSPA`: use `pVals` instead of `N_log2R` or `Z_log2R`; 
  # currently Protein[_impNA]_pVals.txt does not contain `scale_log2_r` info in the file name
  #   and cannot tell `pVals` are from `N_log2R` or `Z_log2R`;  
  #   therefore, `scale_log2r` will be matched to those in `prnSig()` and indicated in 
  #   the names of out files as `_N[_impNA].txt` or `_Z[_impNA].txt` to inform user the `scale_log_r` status
  # `complete_cases` is for entrez IDs, pVals, log2FC
  # "id" for tibbling rownames
  
  stopifnot(nrow(label_scheme_sub) > 0)
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots = rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]

  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  if (nrow(df) == 0) stop("No data available after row filtration.", call. = FALSE)
  
  if (purrr::is_empty(fmls)) stop("Formula(s) of contrasts not available.", call. = FALSE)

  species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
  gsets <- load_dbs(gset_nms = gset_nms, species = species)
  
  stopifnot(length(gsets) > 0)

  fml_nms <- names(df) %>% 
    .[grepl("pVal\\s*\\(", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique() %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls %>% .[names(.) %in% fml_nms]
  fml_nms <- fml_nms %>% .[map_dbl(., ~ which(.x == names(fmls)))]

  if (rlang::is_empty(fml_nms)) {
    stop("No formula matached; compare the formula name(s) with those in `prnSig(..)`")
  }

  col_ind <- purrr::map(fml_nms, ~ grepl(.x, names(df))) %>%
    dplyr::bind_cols() %>%
    rowSums() %>%
    `>`(0)
  
  stopifnot(sum(col_ind) > 0)

  switch (anal_type,
          GSPA = purrr::pwalk(list(fmls, 
                                   fml_nms, 
                                   pval_cutoff, 
                                   logFC_cutoff, 
                                   gspval_cutoff, 
                                   gslogFC_cutoff, 
                                   min_size, 
                                   max_size, 
                                   min_delta, 
                                   min_greedy_size
                              ), 
                              fml_gspa, 
                              df = df, 
                              col_ind = col_ind, 
                              id = !!id, 
                              gsets = gsets, 
                              label_scheme_sub = label_scheme_sub, 
                              complete_cases = complete_cases, 
                              scale_log2r = scale_log2r, 
                              filepath = filepath, 
                              filename = filename, 
                              method = method, 
                              !!!dots), 
          GSEA = purrr::pwalk(list(fmls, 
                                   fml_nms, 
                                   var_cutoff, 
                                   pval_cutoff, 
                                   logFC_cutoff, 
                                   gspval_cutoff, 
                                   gslogFC_cutoff, 
                                   min_size, 
                                   max_size
                              ), 
                              fml_gsea, 
                              df = df, 
                              col_ind = col_ind, 
                              id = !!id, 
                              gsets = gsets, 
                              label_scheme_sub = label_scheme_sub, 
                              complete_cases = complete_cases, 
                              scale_log2r = scale_log2r, 
                              filepath = filepath, 
                              filename = filename, 
                              !!!dots), 
    stop("Unhandled task.")
  )
}


#' A helper function for GSPA
#' 
#' @param fml A character string; the formula used in \link{prnSig}.
#' @param fml_nm A character string; the name of \code{fml}.
#' @param col_ind Numeric vector; the indexes of columns for the ascribed \code{fml_nm}.
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams gspaTest
#' @inheritParams gsVolcano
#' @import purrr dplyr rlang
#' @importFrom magrittr %>% %T>%
fml_gspa <- function (fml, fml_nm, pval_cutoff, logFC_cutoff, gspval_cutoff, gslogFC_cutoff, 
                      min_size, max_size, min_delta, min_greedy_size, 
                      df, col_ind, id, gsets, label_scheme_sub, complete_cases, scale_log2r, 
                      filepath, filename, method, ...) {

  dir.create(file.path(filepath, fml_nm), recursive = TRUE, showWarnings = FALSE)
  id <- rlang::as_string(rlang::enexpr(id))
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  df <- df %>% prep_gspa(id = "entrez", fml_nm = fml_nm, col_ind = col_ind, 
                         pval_cutoff = pval_cutoff, logFC_cutoff = logFC_cutoff) 
  
  if (complete_cases) {
    df <- df[complete.cases(df), ]
    message("Complete cases against: ", reduce(names(df), paste, sep = ", "))
  }
  
  df <- df %>% 
    dplyr::mutate(p_val = -log10(p_val)) %>% 
    dplyr::mutate(valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>% 
    dplyr::mutate(valence = factor(valence, levels = c("neg", "pos"))) %>% 
    tidyr::complete(entrez, contrast, valence)
  
  # no re-`arrange` df, df_sub for limma after this point
  df <- df %>% dplyr::arrange(entrez, contrast, valence)

  gsets <- gsets %>% 
    .[!grepl("molecular_function$", names(.))] %>% 
    .[!grepl("cellular_component$", names(.))] %>% 
    .[!grepl("biological_process$", names(.))] %>%     
    .[purrr::map_lgl(., ~ length(.x) >= min_size)] %>% 
    .[purrr::map_lgl(., ~ length(.x) <= max_size)]
  
  res_pass <- local({
    contrast_groups <- unique(df$contrast) %>% as.character()
    
    contrast_table <- tibble(contrast = rep(contrast_groups, 2), 
                             valence = rep(c("neg", "pos"), length(contrast_groups)), 
                             p_val = rep(-log10(0.05), 2 * length(contrast_groups))) %>% 
      dplyr::mutate(contrast = factor(contrast, levels = levels(df$contrast)), 
                    valence = factor(valence, levels(df$valence)))
    
    res <- switch(method, 
                  limma = purrr::map(gsets, 
                                     gspa_summary_limma, 
                                     df = df, 
                                     min_size = min_size, 
                                     min_delta = min_delta), 
                  mean = purrr::map(gsets, 
                                    gspa_summary_mean, 
                                    df = df, 
                                    min_size = min_size, 
                                    min_delta = min_delta, 
                                    gspval_cutoff = gspval_cutoff, 
                                    gslogFC_cutoff = gslogFC_cutoff), 
                  rlang::abort("Unknown `method`."))

    idx <- purrr::map_dbl(res, is.null)
    
    res <- res[!idx] %>% 
      do.call(rbind, .) %>% 
      tibble::rownames_to_column("term") %>% 
      dplyr::mutate(term = factor(term)) %>% 
      dplyr::mutate_at(.vars = grep("^pVal\\s+\\(", names(.)), ~ abs(.x)) %>% 
      dplyr::group_by(term) %>% 
      dplyr::arrange(term) %>% 
      data.frame(check.names = FALSE) 
    
    pass <- purrr::map(res[paste0("pVal (", contrast_groups, ")")], ~ .x <= gspval_cutoff) %>%
      dplyr::bind_cols() %>%
      rowSums() %>%
      `>`(0)
    
    res_pass <- res[pass, ] %>% 
      dplyr::mutate_at(vars(grep("^pVal", names(.))), ~ replace(.x, .x < 1E-50, 1E-50))
  })
  
  if (nrow(res_pass) == 0) return(NULL)
  
  out <- local({
    adjp <- res_pass %>% 
      dplyr::select(grep("^pVal\\s+\\(", names(.))) %>% 
      purrr::map(~ p.adjust(., "BH")) %>% 
      data.frame(check.names = FALSE) %>% 
      `names<-`(gsub("pVal", "q_val", colnames(.))) %>% 
      `names<-`(gsub("^q_val\\s+\\((.*)\\)", "\\1", names(.))) %>% 
      tidyr::gather(key = contrast, value = q_val) %>% 
      dplyr::select(-contrast)
    
    log2fc <- res_pass %>% 
      dplyr::select(grep("^log2Ratio\\s+\\(", names(.))) %>% 
      `names<-`(gsub("^log2Ratio\\s+\\((.*)\\)", "\\1", names(.))) %>% 
      tidyr::gather(key = contrast, value = log2fc) %>% 
      dplyr::select(-contrast)
    
    pval <- res_pass %>% 
      dplyr::select(-grep("^log2Ratio\\s+\\(", names(.))) %>% 
      `names<-`(gsub("^pVal\\s+\\((.*)\\)", "\\1", names(.))) %>% 
      tidyr::gather(key = contrast, value = p_val, -term)
    
    dplyr::bind_cols(pval, adjp, log2fc) %>% 
      dplyr::mutate_at(.vars = grep("^p_val$|^adjP$", names(.)), format, scientific = TRUE, digits = 2) %>%
      dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2)
  })

  sig_sets <- purrr::map2(gsets, names(gsets), ~ {
    tibble::tibble(id = as.character(.x)) %>% 
      dplyr::mutate(set = .y)
  }) %>% 
    dplyr::bind_rows() %>% 
    `names<-` (c("id", "term")) %>% 
    dplyr::filter(id %in% unique(df$entrez)) %>% 
    dplyr::select(c("term", "id")) %>% 
    dplyr::filter(term %in% res_pass$term)
  
  out <- out %>% dplyr::mutate(term = as.character(term))
  
  out <- sig_sets %>% 
    dplyr::group_by(term) %>% 
    dplyr::summarise(size = n()) %>% # the number of entries matched to input `df`
    dplyr::right_join(out, by = "term")

  res_greedy <- sig_sets %>% 
    RcppGreedySetCover::greedySetCover(FALSE) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_resgreedy.txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE) %>% 
    dplyr::group_by(term) %>% 
    dplyr::summarise(ess_size = n()) %>% 
    dplyr::filter(ess_size >= min_greedy_size)
  
  sig_sets <- res_greedy %>% 
    dplyr::right_join(out, by = "term") %>% 
    dplyr::mutate(is_essential = ifelse(!is.na(ess_size), TRUE, FALSE)) %>% 
    dplyr::select(term, is_essential, size, ess_size, contrast, p_val, q_val, log2fc) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, ".txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE) %>% 
    dplyr::filter(is_essential) %>% 
    dplyr::filter(!duplicated(term)) %>% 
    dplyr::select(-contrast, -p_val, -q_val, -log2fc) %T>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_essmeta.txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE) %>% 
    dplyr::select(term, ess_size) %>% 
    dplyr::right_join(sig_sets, by = "term")

  map_essential(sig_sets) %>% 
    write.table(file.path(filepath, fml_nm, paste0(fn_prefix, "_essmap.txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE)	
}


#' check the size of a gene set
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
ok_min_size <- function (df, min_delta, gspval_cutoff, gslogFC_cutoff) {
  data <- df %>% 
    dplyr::group_by(contrast, valence) %>% 
    dplyr::filter(!is.na(log2Ratio), !is.na(p_val)) %>% 
    dplyr::mutate(n = n()) %>% 
    tidyr::unite(con_val, contrast, valence, sep = ".", remove = FALSE) 
  
  delta_p <- data %>% 
    dplyr::summarise(p_val = mean(p_val, na.rm = TRUE)) %>% 
    dplyr::mutate(p_val = ifelse(is.nan(p_val), 0, p_val)) %>% 
    tidyr::spread(key = contrast, value = p_val) %>% 
    dplyr::ungroup(valence) %>% 
    dplyr::select(-valence)
  delta_p[is.na(delta_p)] <- 0
  if (nrow(delta_p) == 2) delta_p <- delta_p[2, ] - delta_p[1, ]

  delta_fc <- data %>% 
    dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) %>% 
    dplyr::mutate(log2Ratio = ifelse(is.nan(log2Ratio), 0, log2Ratio)) %>% 
    tidyr::spread(key = contrast, value = log2Ratio) %>% 
    dplyr::ungroup(valence) %>% 
    dplyr::select(-valence)
  delta_fc[is.na(delta_fc)] <- 0
  if (nrow(delta_fc) == 2) delta_fc <- delta_fc[2, ] + delta_fc[1, ]
  
  delta_n <- data %>% 
    dplyr::summarise(n = mean(n, na.rm = TRUE)) %>% # empty levels will be kept
    tidyr::spread(key = contrast, value = n) %>% 
    dplyr::ungroup(valence) %>% 
    dplyr::select(-valence)
  delta_n[is.na(delta_n)] <- 0
  if (nrow(delta_n) == 2) delta_n <- delta_n[2, ] - delta_n[1, ]

  sign(delta_p) == sign(delta_n) & 
    abs(delta_p) >= -log10(gspval_cutoff) & 
    abs(delta_n) >= min_delta & 
    abs(delta_fc) >= gslogFC_cutoff
    
}


#' gspa pVal calculations using geomean
#' 
#' @param gset Character string; a gene set indicated by \code{gset_nms}.
#' @inheritParams prnHist
#' @inheritParams prnGSPA
gspa_summary_mean <- function(gset, df, min_size = 10, min_delta = 4, 
                              gspval_cutoff = 0.05, gslogFC_cutoff = log2(1)) {
  df <- df %>% dplyr::filter(.[["entrez"]] %in% gset)
  
  if (length(unique(df$entrez)) < min_size) return(NULL)
  
  df <- df %>% dplyr::group_by(contrast, valence)

  ok_delta_n <- ok_min_size(df = df, min_delta = min_delta, gspval_cutoff = gspval_cutoff, gslogFC_cutoff = gslogFC_cutoff)

  # penalty term
  dfw <- df %>% dplyr::mutate(p_val = ifelse(is.na(p_val), 1, p_val))

  delta_p <- dfw %>% 
    dplyr::summarise(p_val = mean(p_val, na.rm = TRUE)) %>% 
    dplyr::mutate(p_val = replace(p_val, is.nan(p_val), 1)) %>% 
    tidyr::spread(key = contrast, value = p_val) %>% 
    dplyr::select(-valence) %>% 
    `colnames<-`(paste0("pVal (", colnames(.), ")")) 
  if (nrow(delta_p) == 2) delta_p <- delta_p[2, ] - delta_p[1, ]
  delta_p <- abs(delta_p)
  delta_p <- 1/10^delta_p
  delta_p[] <- purrr::map2(as.list(delta_p), as.list(ok_delta_n), ~ {if (.y) .x else 1}) 

  delta_fc <- dfw %>% 
    dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) %>% 
    dplyr::mutate(log2Ratio = replace(log2Ratio, is.nan(log2Ratio), 0)) %>% 
    tidyr::spread(key = contrast, value = log2Ratio) %>% 
    dplyr::select(-valence) %>% 
    `colnames<-`(paste0("log2Ratio (", colnames(.), ")")) %>% 
    abs(.) 
  if (nrow(delta_fc) == 2) delta_fc <- delta_fc[2, ] - delta_fc[1, ]
  
  return(dplyr::bind_cols(delta_p, delta_fc))
}  


#' gspa pVal calculations using limma
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams gspa_summary_mean
gspa_summary_limma <- function(gset, df, min_size = 10, min_delta = 4) {
  df_sub <- df %>% dplyr::filter(.[["entrez"]] %in% gset) 
  if (length(unique(df_sub$entrez)) < min_size) return(NULL)
  
  lm_gspa(df_sub, min_delta)
}  


#' significance tests of pVals between the up and the down groups
#' 
#' @inheritParams prnHist
#' @inheritParams prnGSPA
lm_gspa <- function(df, min_delta) {
  delta_fc <- df %>% 
    dplyr::group_by(contrast, valence) %>% 
    dplyr::summarise(log2Ratio = mean(log2Ratio, na.rm = TRUE)) %>% 
    dplyr::mutate(log2Ratio = replace(log2Ratio, is.nan(log2Ratio), 0)) %>% 
    tidyr::spread(key = contrast, value = log2Ratio) %>% 
    dplyr::select(-valence) 
  
  nms <- colnames(delta_fc)
  delta_fc <- delta_fc%>% 
    `colnames<-`(paste0("log2Ratio (", colnames(.), ")")) %>% 
    abs(.) 
  if (nrow(delta_fc) == 2) delta_fc <- delta_fc[2, ] - delta_fc[1, ]
  
  ok_delta_n <- ok_min_size(df = df, min_delta = min_delta, gspval_cutoff = 1)

  data <- df %>% 
    dplyr::group_by(contrast, valence) %>% 
    dplyr::select(-log2Ratio) %>% 
    tidyr::unite(sample, entrez, contrast, valence, sep = ".", remove = TRUE) %>% 
    tidyr::spread(key = sample, value = p_val) 
  
  label_scheme_gspa <- tibble(Sample_ID = names(data), Group = Sample_ID) %>% 
    dplyr::mutate(Group = gsub("^\\d+\\.", "", Group)) %>% 
    dplyr::mutate(Group = factor(Group))
  
  design <- model.matrix(~0+label_scheme_gspa$Group) %>% 
    `colnames<-`(levels(label_scheme_gspa$Group)) %>% 
    data.frame(check.names = TRUE)

  neg_cols <- colnames(design)[seq_along(design) %% 2 %>% as.logical()]
  pos_cols <- colnames(design)[!(seq_along(design) %% 2)]
  
  contrs <- paste(pos_cols, neg_cols, sep = "-")
  
  contr_mat <- makeContrasts(contrasts = contrs, levels = data.frame(design)) %>% 
    `colnames<-`(contrs) %>% 
    `rownames<-`(colnames(design))
  
  fit <- data %>%
    lmFit(design = design) %>%
    contrasts.fit(contr_mat) %>%
    eBayes()
  
  delta_p <- fit$p.value %>%
    replace(., is.na(.), 1) %>% 
    data.frame(check.names = FALSE) %>%
    `colnames<-`(paste0("pVal (", nms, ")"))

  delta_p[] <- purrr::map2(as.list(delta_p), as.list(ok_delta_n), ~ {if (.y) .x else 1}) 

  return(dplyr::bind_cols(delta_p, delta_fc))
}


#' A helper function for GSPA
#'
#' @inheritParams prnHist
#' @inheritParams prnGSPA
#' @inheritParams fml_gspa
#' 
#' @import purrr dplyr rlang
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
prep_gspa <- function(df, id, fml_nm, col_ind, pval_cutoff = 5E-2, logFC_cutoff = log2(1.2)) {
  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- df %>%
    dplyr::select(grep(fml_nm, names(.), fixed = TRUE)) %>%
    `colnames<-`(gsub(paste0(fml_nm, "."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .) %>% 
    rm_pval_whitespace() %>% 
    dplyr::select(id, grep("^pVal|^log2Ratio", names(.))) %>% 
    dplyr::mutate(!!id := as.character(.[[id]]))
  
  contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
    gsub("^log2Ratio\\s+\\(|\\)$", "", .)
  
  pvals <- df %>% 
    dplyr::select(-grep("^log2Ratio\\s+\\(", names(.))) %>% 
    `names<-`(gsub("^pVal\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = contrast, value = p_val, -id)
  
  log2rs <- df %>% 
    dplyr::select(-grep("^pVal\\s+\\(", names(.))) %>% 
    `names<-`(gsub("^log2Ratio\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = contrast, value = log2Ratio, -id) %>% 
    dplyr::select(-id, -contrast)
  
  df <- dplyr::bind_cols(pvals, log2rs) %>% 
    dplyr::filter(p_val <= pval_cutoff) %>% # NA will be removed
    dplyr::filter(abs(log2Ratio) >= logFC_cutoff) %>% 
    dplyr::filter(!is.na(id)) %>% 
    dplyr::mutate(contrast = factor(contrast, levels = contrast_groups)) %>%
    dplyr::arrange(contrast) 

  return(df)
}


#' A helper function for mapping between gene sets and essential gene sets
#'
#' @param sig_sets A data frame containing the gene sets that are significant
#'   under given criteria.
#' @import purrr dplyr rlang RcppGreedySetCover
#' @importFrom magrittr %>%
map_essential <- function (sig_sets) {
  ess_terms <- sig_sets %>% 
    dplyr::filter(!is.na(ess_size), !duplicated(term)) %>% 
    dplyr::select(term) %>% 
    unlist

  sig_greedy_sets <- sig_sets %>% dplyr::filter(term %in% ess_terms)

  res <- purrr::map(unique(sig_sets$term), ~ {
    curr_sig_set <- sig_sets %>% 
      dplyr::filter(term %in% .x)
    
    sig_greedy_sets %>% 
      dplyr::filter(id %in% curr_sig_set$id) %>% 
      dplyr::group_by(term) %>%
      dplyr::summarise(size = n()) %>% 
      dplyr::left_join(curr_sig_set %>% 
                         dplyr::select(c("term", "ess_size")) %>% 
                         dplyr::filter(!duplicated(term)), 
                       by = "term") %>% 
      dplyr::mutate(fraction = size/nrow(curr_sig_set)) %>% 
      dplyr::mutate(curr_sig_set = .x)
  }, sig_sets, sig_greedy_sets) %>% 
    dplyr::bind_rows() %>% 
    dplyr::group_by(curr_sig_set) %>%
    dplyr::rename(ess_term = "term", term = "curr_sig_set") %>%
    dplyr::select(c("term", "ess_term", "size", "ess_size", "fraction")) %>% 
    dplyr::mutate(fraction = round(fraction, digits = 3)) %>% 
    dplyr::ungroup(ess_term) %>% 
    dplyr::mutate(distance = 1 - fraction) %>% 
    dplyr::mutate(term = factor(term, levels = unique(as.character(term)))) %>% 
    dplyr::mutate(ess_term = factor(ess_term, levels = levels(term))) %>%
    dplyr::mutate(idx = as.numeric(term), ess_idx = as.numeric(ess_term))
}


#'Heat map visualization of GSPA results
#'
#'\code{prnGSPAHM} visualizes distance heat maps and networks between essential
#'and all gene sets.
#'
#'The list of gene sets and the associative quality metrics of \code{size} and
#'\code{ess_size} are assessed after data filtration with the criteria specified
#'by arguments \code{pval_cutoff} and \code{logFC_cutoff}, as well as optional
#'varargs of \code{filter_}.
#'
#'@section \code{Protein_GSPA_[...].txt}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Description}\cr term \tab a gene
#'  set term \cr is_essential \tab a logical indicator of gene set essentiality
#'  \cr size \tab the number of IDs under a \code{term} \cr ess_size \tab the
#'  number of IDs that can be found under a corresponding essential set \cr
#'  contrast \tab a contrast of sample groups \cr p_val \tab significance p
#'  values \cr q_val \tab \code{p_val} with \code{BH} adjustment of multiple
#'  tests \cr log2fc \tab the fold change of a gene set at logarithmic base of 2
#'  \cr }
#'
#'@section \code{Protein_GSPA_[...]essmap.txt}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Descrption}\cr term \tab a gene
#'  set term \cr ess_term \tab an essential gene set term \cr size \tab the
#'  number of IDs under a \code{term} with matches to an \code{ess_term} \cr
#'  ess_size \tab the number of essential IDs under a \code{term} with matches
#'  to an \code{ess_term} \cr fraction \tab a fraction of matches in IDs between
#'  a \code{term} and a \code{ess_term} \cr distance \tab 1 - \code{fraction}
#'  \cr idx \tab a numeric index of \code{term} \cr ess_idx \tab a numeric index
#'  of \code{ess_term} \cr }
#'
#'@inheritParams plot_prnTrend
#'@inheritParams  prnEucDist
#'@param impute_na Logical; at TRUE, input files with \code{_impNA[...].txt} in
#'  name will be loaded. Otherwise, files without \code{_impNA} in name will be
#'  taken. An error will be thrown if no files are matched under given
#'  conditions. The default is FALSE.
#'@param fml_nms Character string or vector; the formula name(s). By default,
#'  the formula(s) will match those used in \code{\link{pepSig}} or
#'  \code{\link{prnSig}}.
#'@param annot_cols A character vector of column keys that can be found in
#'  \code{_essmap.txt}. The values under the selected keys will be used to
#'  color-code enrichment terms on the top of heat maps. The default is NULL
#'  without column annotation.
#'@param annot_rows A character vector of column keys that can be found from
#'  \code{_essmeta.txt} . The values under the selected keys will be used to
#'  color-code essential terms on the side of heat maps. The default is NULL
#'  without row annotation.
#'@param ... \code{filter2_}: Variable argument statements for the row
#'  filtration against data in secondary file(s) of \code{_essmap.txt}. Each
#'  statement contains to a list of logical expression(s). The \code{lhs} needs
#'  to start with \code{filter2_}. The logical condition(s) at the \code{rhs}
#'  needs to be enclosed in \code{exprs} with round parenthesis. For example,
#'  \code{distance} is a column key in \code{Protein_GSPA_Z_essmap.txt}. The
#'  statement \code{filter2_ = exprs(distance <= .95),} will remove entries with
#'  \code{distance > 0.95}. See also \code{\link{normPSM}} for the format of
#'  \code{filter_} statements against primary data. \cr \cr \code{arrange2_}:
#'  Variable argument statements for the row ordering against data in secondary
#'  file(s) of \code{_essmap.txt}. The \code{lhs} needs to start with
#'  \code{arrange2_}. The expression(s) at the \code{rhs} needs to be
#'  enclosed in \code{exprs} with round parenthesis. For example,
#'  \code{distance} and \code{size} are column keys in
#'  \code{Protein_GSPA_Z_essmap.txt}. The statement \code{arrange2_ =
#'  exprs(distance, size),} will order entries by \code{distance}, then by
#'  \code{size}. See also \code{\link{prnHM}} for the format of \code{arrange_}
#'  statements against primary data. \cr \cr Additional arguments for
#'  \code{\link[pheatmap]{pheatmap}}, i.e., \code{fontsize }... \cr \cr Note
#'  arguments disabled from \code{pheatmap}: \cr \code{annotation_col}; instead
#'  use keys indicated in \code{annot_cols} \cr \code{annotation_row}; instead
#'  use keys indicated in \code{annot_rows}
#'@import purrr
#'
#'@example inst/extdata/examples/prnGSPAHM_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
prnGSPAHM <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, fml_nms = NULL, 
                       annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                       df2 = NULL, filename = NULL, ...) {
  check_dots(c("id", "anal_type", "df", "filepath"), ...)
  dir.create(file.path(dat_dir, "Protein\\GSPA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)
  
  df2 <- rlang::enexpr(df2)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  annot_rows <- rlang::enexpr(annot_rows) 
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots, anal_type = "GSPA")
  
  reload_expts()

  info_anal(id = !!id, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, df2 = !!df2, filepath = NULL, filename = !!filename, 
            anal_type = "GSPA_hm")(annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, 
                                   annot_rows = !!annot_rows, !!!dots)
}


#' Plot distance heat map of GSPA
#' 
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams gspaMap
gspaHM <- function(scale_log2r, complete_cases, impute_na, df2, filepath, filename, ...) {
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  if (purrr::is_empty(fmls))
    stop("No formula(s) of contrasts available.", call. = TRUE)
  
  fml_nms <- names(fmls)
  if (length(fml_nms) > 0) {
    purrr::walk(fml_nms, byfml_gspahm, df2, filepath, filename, scale_log2r, impute_na, !!!dots)
  }
}


#' A helper function for gspaHM by formula names
#' 
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import purrr dplyr rlang pheatmap networkD3
#' @importFrom magrittr %>%
byfml_gspahm <- function (fml_nm, df2, filepath, filename, scale_log2r, impute_na, ...) {
  ins <- list.files(path = file.path(filepath, fml_nm), pattern = "_essmap\\.txt$")
  
  if (purrr::is_empty(ins)) {
    message("No GSPA results at ", fml_nm)
    return(NULL)
  }
  
  if (is.null(df2)) {
    ins <- ins %>% 
      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} %>% 
      {if (scale_log2r) .[grepl("_GSPA_Z", .)] else .[grepl("_GSPA_N", .)]}
  
    if (rlang::is_empty(ins)) {
      stop("No GSPA inputs correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
           call. = FALSE)
    }
  } else {
    local({
      non_exists <- df2 %>% .[! . %in% ins]
      if (!purrr::is_empty(non_exists)) {
        stop("Missing _essmap file(s): ", purrr::reduce(non_exists, paste, sep = ", "), call. = FALSE)
      }
    })
    
    if (purrr::is_empty(df2)) stop("Input file(s) not found.", call. = FALSE)
    ins <- ins %>% .[. %in% df2]    
  }

  meta_ins <- list.files(path = file.path(filepath, fml_nm), pattern = "_essmeta\\.txt$") 
  meta_ins <- local({
    required <- ins %>% gsub("_essmap\\.txt$", "_essmeta.txt", .)
    # non_exists <- meta_ins %>% .[! . %in% required]
    # if (!purrr::is_empty(non_exists)) {
    #   stop("Missing _essmeta file(s): ", purrr::reduce(non_exists, paste, sep = ", "), call. = FALSE)
    # }
    meta_ins <- meta_ins %>% .[. %in% required]
  })
  
  purrr::walk2(ins, meta_ins, byfile_gspahm, fml_nm, filepath, filename, scale_log2r, impute_na, ...)
}


#' A helper function for gspaHM by input file names
#' 
#' @param ess_in A list of file names for input data
#' @param meta_in A list of file names for input metadata
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams fml_gspa
#' @import purrr dplyr rlang pheatmap networkD3
#' @importFrom magrittr %>%
byfile_gspahm <- function (ess_in, meta_in, fml_nm, filepath, filename, scale_log2r, impute_na, ...) {
  custom_prefix <- gsub("(.*_{0,1})Protein_GSPA.*", "\\1", ess_in)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) # %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  filename <- paste0(custom_prefix, fn_prefix, ".", fn_suffix)

  dots <- rlang::enexprs(...)
  if (!purrr::is_empty(dots)) {
    if (any(grepl("^filter_", names(dots)))) {
      stop("Primary `filter_` depreciated; use secondary `filter2_`.")
    }      
  }
  filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
  arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
  select2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select2_", names(.))]
  dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots, select2_dots)]
  
  all_by_greedy <- tryCatch(readr::read_tsv(file.path(filepath, fml_nm, ess_in), 
                                            col_types = cols(ess_term = col_factor(), 
                                                             size = col_double(), 
                                                             ess_size = col_double(), 
                                                             fraction = col_double(),
                                                             distance = col_double(), 
                                                             idx = col_double(),
                                                             ess_idx = col_double())), 
                            error = function(e) NA) %>% 
    filters_in_call(!!!filter2_dots) %>% 
    arrangers_in_call(!!!arrange2_dots)
  
  rm(filter2_dots, arrange2_dots, select2_dots)
  
  if (nrow(all_by_greedy) == 0) stop("No GSPA terms available after data filtration.")
  
  if (max(all_by_greedy$distance) == 0) stop("Identical, all-zero distance detected; 
                                             try lower the criteria in data filtrations
                                             or rerun `prnGSPA` at more relaxed `gspval_cutoff` threshold.")
  
  if (!is.null(dim(all_by_greedy))) {
    message(paste("File loaded:", ess_in, "at", fml_nm))
  } else {
    stop("Essential GSPA not found.")
  }
  
  ess_vs_all <- all_by_greedy %>% 
    dplyr::select(-which(names(.) %in% c("size", "ess_size", "idx", "ess_idx", "distance"))) %>% 
    tidyr::spread(term, fraction) %>% 
    dplyr::mutate_at(vars(which(names(.) != "ess_term")), ~ {1 - .x}) %>% 
    tibble::column_to_rownames("ess_term")
  
  d_row <- dist(ess_vs_all)
  d_col <- dist(t(ess_vs_all))
  max_d_row <- max(d_row, na.rm = TRUE)
  max_d_col <- max(d_col, na.rm = TRUE)
  
  d_row[is.na(d_row)] <- 1.2 * max_d_row
  d_col[is.na(d_col)] <- 1.2 * max_d_col
  
  if (is.infinite(max_d_row)) {
    cluster_rows <- FALSE
  } else {
    cluster_rows <- hclust(d_row)
  }
  
  if (is.infinite(max_d_col)) {
    cluster_cols <- FALSE
  } else {
    cluster_cols <- hclust(d_col)
  }
  
  if (is.null(dots$color)) {
    mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
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
  
  if (is.null(dots$annot_rows)) {
    annot_rows <- NULL
  } else {
    annot_rows <- eval(dots$annot_rows, env = caller_env()) 
  }
  
  if (is.null(annot_cols)) {
    annotation_col <- NA
  } else {
    annotation_col <- gspa_colAnnot(annot_cols = annot_cols, 
                                    df = all_by_greedy, 
                                    sample_ids = all_by_greedy$term)
  }
  
  if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
    colnames(annotation_col) <- annot_colnames
  }
  
  if (is.null(annot_rows)) {
    annotation_row <- NA
  } else {
    ess_meta <- tryCatch(readr::read_tsv(file.path(filepath, fml_nm, meta_in), 
                                         col_types = cols(term = col_factor())), 
                         error = function(e) NA)
    
    annot_rows <- annot_rows %>% .[. %in% names(ess_meta)]
    
    if (is_empty(annot_rows)) {
      annotation_row <- NA
    } else {
      annotation_row <- ess_meta %>% 
        dplyr::rename(ess_term = term) %>% 
        dplyr::select("ess_term", annot_rows) %>% 
        dplyr::mutate(ess_term = factor(ess_term, levels = levels(all_by_greedy$ess_term))) %>% 
        tibble::column_to_rownames("ess_term")
    }
  }
  
  if (is.null(dots$annotation_colors)) {
    annotation_colors <- setHMColor(annotation_col)
  } else if (is.na(dots$annotation_colors)) {
    annotation_colors <- NA
  } else {
    annotation_colors <- eval(dots$annotation_colors, env = caller_env())
  }
  
  nrow <- nrow(ess_vs_all)
  ncol <- ncol(ess_vs_all)
  max_width <- 44
  max_height <- 44
  
  if (is.null(dots$width)) {
    if (ncol <= 150) {
      fontsize <- 12
      fontsize_col <- cellwidth <- 12
      show_colnames <- TRUE
      width <- pmax(ncol/1.2, 8)
    } else {
      fontsize <- fontsize_col <- 1
      cellwidth <- NA
      show_colnames <- FALSE
      width <- max_width
    }
  } else {
    width <- dots$width
    if (is.null(dots$fontsize)) fontsize <- 1 else fontsize <- dots$fontsize
    if (is.null(dots$fontsize_col)) fontsize_col <- 1 else fontsize_col <- dots$fontsize_col
    if (is.null(dots$cellwidth)) cellwidth <- NA else cellwidth <- dots$cellwidth
    if (is.null(dots$show_colnames)) show_colnames <- FALSE else show_colnames <- dots$show_colnames
  }
  
  if (is.null(dots$height)) {
    if (nrow <= 150) {
      fontsize <- 12 
      fontsize_row <- cellheight <- 12
      show_rownames <- TRUE
      height <- pmax(nrow/1.2, 8)
    } else {
      fontsize <- fontsize_row <- 1
      cellheight <- NA
      show_rownames <- FALSE
      height <- max_height
    }
  } else {
    height <- dots$height
    if (is.null(dots$fontsize)) fontsize <- 1 else fontsize <- dots$fontsize
    if (is.null(dots$fontsize_row)) fontsize_row <- 1 else fontsize_row <- dots$fontsize_row
    if (is.null(dots$cellheight)) cellheight <- NA else cellheight <- dots$cellheight
    if (is.null(dots$show_rownames)) show_rownames <- FALSE else show_rownames <- dots$show_rownames
  }
  
  if ((!is.na(width)) & (width > max_width)) {
    warning("The plot width is set to ", max_width, call. = FALSE)
    width <- max_width
    height <- pmin(max_height, width * nrow / ncol)
  } 
  
  if ((!is.na(height)) & (height > max_height)) {
    warning("The plot height is set to ", max_height, call. = FALSE)
    height <- max_height
    width <- pmin(max_width, height * ncol / nrow)
  }
  
  dots <- dots %>% 
    .[!names(.) %in% c("annot_cols", "annot_colnames", "annot_rows", 
                       "mat", "filename", "annotation_col", "annotation_row", 
                       "color", "annotation_colors", "breaks", 
                       "cluster_rows", "cluster_cols", 
                       "width", "height", "fontsize_row", "fontsize_col", 
                       "cellheight", "cellwidth")]
  
  ph <- my_pheatmap(
    mat = ess_vs_all,
    filename = file.path(filepath, fml_nm, filename),
    annotation_col = annotation_col,
    annotation_row = annotation_row, 
    color = mypalette,
    annotation_colors = annotation_colors, 
    breaks = NA, # not used
    cluster_rows = cluster_rows, 
    cluster_cols = cluster_cols,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    cellheight = cellheight,
    cellwidth = cellwidth,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    width = width, 
    height = height, 
    !!!dots,
  )
  
  
  # networks
  cluster <- data.frame(cluster = cutree(ph$tree_col, h = max_d_col)) %>% 
    tibble::rownames_to_column("term")
  
  all_by_greedy <- local({
    essterms <- all_by_greedy$ess_term %>% unique() %>% as.character()
    terms <- all_by_greedy$term %>% unique() %>% .[! . %in% essterms] %>% as.character()
    universe <- union(essterms, terms)
    
    all_by_greedy <- all_by_greedy %>% 
      dplyr::select(-idx, -ess_idx) %>% 
      dplyr::mutate(term = factor(term, levels = universe)) %>% 
      dplyr::mutate(ess_term = factor(ess_term, levels = universe)) %>%
      dplyr::mutate(source = as.numeric(term), target = as.numeric(ess_term))
    
    min_target <- min(all_by_greedy$target, na.rm = TRUE)
    
    all_by_greedy %>% 
      dplyr::mutate(source = source - min_target, target = target - min_target) %>% 
      dplyr::mutate(term = as.character(term)) %>% 
      dplyr::left_join(cluster, by = "term") %>% 
      dplyr::mutate(term = factor(term, levels(ess_term)))
  })
  
  my_nodes <- local({
    my_nodes <- all_by_greedy %>% 
      dplyr::arrange(-size) %>% 
      dplyr::select(c("term", "cluster", "size")) %>% 
      dplyr::filter(!duplicated(term)) %>% 
      dplyr::arrange(cluster) %>% 
      dplyr::arrange(term)
    
    # an `ess_term` may be missing from `term` after data row filtration
    temp_nodes <- all_by_greedy %>% 
      dplyr::select(c("ess_term", "cluster", "size")) %>% 
      dplyr::filter(!duplicated(ess_term)) %>% 
      dplyr::filter(! ess_term %in% my_nodes$term) %>% 
      dplyr::arrange(cluster) %>% 
      dplyr::arrange(ess_term) %>% 
      dplyr::rename(term = ess_term)
    
    my_nodes <- bind_rows(my_nodes, temp_nodes) %>% 
      dplyr::arrange(cluster) %>% 
      dplyr::arrange(term)
  }) %>% 
    data.frame(check.names = FALSE)
  
  my_links <- all_by_greedy %>% 
    dplyr::arrange(term) %>% 
    dplyr::select(source, target, fraction) %>% 
    dplyr::mutate(fraction = fraction * 10) %>% 
    dplyr::distinct() %>% 
    data.frame(check.names = FALSE)
  
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  networkD3::forceNetwork(Links = my_links, Nodes = my_nodes, Source = "source",
                          Target = "target", Value = "fraction", NodeID = "term", Nodesize = "size", 
                          Group = "cluster", opacity = 0.8, zoom = TRUE) %>% 
    networkD3::saveNetwork(file = file.path(filepath, fml_nm, paste0(fn_prefix, ".html")))
}


#' Sets up the column annotation in heat maps
#' 
#' @param sample_ids A character vecotr containing the sample IDs for an ascribing analysis. 
#' @inheritParams prnEucDist
#' @import dplyr rlang
#' @importFrom magrittr %>%
gspa_colAnnot <- function (annot_cols = NULL, df, sample_ids, annot_colnames = NULL) {
  if (is.null(annot_cols)) return(NA)
  
  exists <- annot_cols %in% names(df)
  
  if (sum(!exists) > 0) {
    warning(paste0("Column '", annot_cols[!exists], "'",
                   " not found in GSPA inputs and will be skipped."))
    annot_cols <- annot_cols[exists]
  }
  
  if (length(annot_cols) == 0) return(NA)
  
  x <- df %>%
    dplyr::filter(term %in% sample_ids, !duplicated(term)) %>%
    dplyr::select(annot_cols, term) %>%
    dplyr::select(which(not_all_NA(.))) %>%
    data.frame(check.names = FALSE) %>%
    `rownames<-`(.[["term"]])
  
  if (any(duplicated(x[["term"]]))) stop("Duplicated terms found\n")
  
  if (!"term" %in% annot_cols) x <- x %>% dplyr::select(-term)
  
  if (ncol(x) == 0) return(NA)
  
  return(x)
}


#'GSVA of protein data
#'
#'\code{prnGSVA} performs the GSVA against protein \code{log2FC}. It is a
#'wrapper of \code{\link[GSVA]{gsva}}.
#'
#'The formula(s) of contrast(s) used in \code{\link{pepSig}} will be taken by
#'default.
#'
#'
#'@inheritParams prnGSPA
#'@param lm_method Character string indicating the linear modeling method for
#'  significance assessment of GSVA enrichment scores. The default is
#'  \code{limma}. At \code{method = lm}, the \code{lm()} in base R will be used
#'  for models without random effects and the \code{\link[lmerTest]{lmer}} will
#'  be used for models with random effects.
#'@param var_cutoff Numeric; the cut-off in the variances of protein log2FC.
#'  Entries with variances smaller than the threshold will be removed from GSVA.
#'  The default is 0.5.
#'@param pval_cutoff Numeric; the cut-off in enrichment pVals. Terms with
#'  enrichment pVals smaller than the threshold will be removed from multiple
#'  test corrections. The default is 1e-04.
#'@param logFC_cutoff Numeric; the cut-off in enrichment log2FC. Terms with
#'  absolute log2FC smaller than the threshold will be removed from multiple
#'  test corrections. The default is at log2(1.1).
#'@param ... \code{filter_}: Logical expression(s) for the row filtration
#'  against data in a primary file of \code{\\Model\\Protein[_impNA]_pVals.txt}.
#'  See also \code{\link{normPSM}} for the format of \code{filter_} statements.
#'  \cr \cr \code{arrange_}: Variable argument statements for the row ordering
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{prnHM}} for the format of \code{arrange_} statements. \cr \cr
#'  Additional arguments for \code{\link[GSVA]{gsva}}
#'@example inst/extdata/examples/prnGSVA_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@import dplyr rlang ggplot2 GSVA
#'@importFrom magrittr %>%
#'@export
prnGSVA <- function (gset_nms = c("go_sets", "c2_msig"), 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     df = NULL, filepath = NULL, filename = NULL, 
                     var_cutoff = .5, pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), 
                     lm_method = "limma", ...) {

  on.exit(
    if (rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(dots) %>% 
        save_call("pepGSVA")
    } else if (rlang::as_string(id) %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(dots) %>% 
        save_call("prnGSVA")
    }
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  err_msg1 <- "Argument `expr`, `gset.idx.list` and `annotation` will be determined automatically.\n"
  if (any(names(rlang::enexprs(...)) %in% c("expr", "gset.idx.list", "annotation"))) 
    stop(err_msg1, call. = FALSE)

  dir.create(file.path(dat_dir, "Protein\\GSVA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	lm_method <- rlang::as_string(rlang::enexpr(lm_method))

	dots <- rlang::enexprs(...)
	fmls <- dots %>% .[grepl("^\\s*~", .)]
	dots <- dots[!names(dots) %in% names(fmls)]

	if (rlang::is_empty(fmls)) {
	  fml_file <-  file.path(dat_dir, "Calls\\prnSig_formulas.rda")
	  if (file.exists(fml_file)) {
	    load(file = fml_file)
	    dots <- c(dots, prnSig_formulas)
	  } else {
	    stop("Run `prnSig()` first.")
	  }
	} else {
	  match_fmls(fmls)
	  dots <- c(dots, fmls)
	}
	
	reload_expts()
	
	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
	info_anal(df = !!df, df2 = NULL, id = !!id, filepath = !!filepath, filename = !!filename, 
	          scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
	          anal_type = "GSVA")(lm_method, gset_nms, var_cutoff, pval_cutoff, logFC_cutoff, !!!dots)
}


#' Perform GSVA tests
#'
#' logFC_cutoff subsets data for adjusted pvals
#' 
#' @inheritParams prnGSVA
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$%
#' @importFrom outliers grubbs.test
#' @importFrom broom.mixed tidy
gsvaTest <- function(df = NULL, id = "entrez", label_scheme_sub = NULL, 
                     filepath = NULL, filename = NULL, 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                     gset_nms = "go_sets", lm_method = "limma", 
                     var_cutoff = .5, pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), 
                     anal_type = "GSVA", ...) {
  
  stopifnot(nrow(label_scheme_sub) > 0)
  
  species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
  gsets <- load_dbs(gset_nms = gset_nms, species = species)
  stopifnot(length(gsets) > 0)
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  if (purrr::is_empty(fmls)) stop("Formula(s) of contrasts not available.", call. = FALSE)  
  
  # `complete_cases` depends on lm contrasts
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    prepDM(id = id, 
           scale_log2r = scale_log2r, 
           sub_grp = label_scheme_sub$Sample_ID, 
           anal_type = anal_type) %>% 
    .$log2R 

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  res_es <- df %>% 
    filterData(var_cutoff = var_cutoff) %>% 
    as.matrix() %>% 
    GSVA::gsva(gsets, !!!dots) %>% 
    data.frame(check.names = FALSE) %>% 
    tibble::rownames_to_column("term") %T>%
    write.table(file.path(filepath, paste0(fn_prefix, "_ES.txt")), sep = "\t", col.names = TRUE, row.names = FALSE)    

  quietly_log <- 
    purrr::map(fmls, ~ purrr::quietly(model_onechannel)(
      res_es %>% tibble::column_to_rownames("term"), id = !!id,
      .x, label_scheme_sub, complete_cases, method = lm_method,
      var_cutoff, pval_cutoff, logFC_cutoff
    )) 
  
  out_path <- file.path(dat_dir, "Protein\\GSVA\\log\\prnGSVA_log.txt")
  
  purrr::map(quietly_log, ~ {
    .x[[1]] <- NULL
    return(.x)
  }) %>% 
    reduce(., `c`) %>% 
    purrr::walk(., write, out_path, append = TRUE)
  
  res <- purrr::map(quietly_log, `[[`, 1) %>%
    do.call("cbind", .) %>% 
    tibble::rownames_to_column("term") %>% 
    rm_pval_whitespace()
  
  rm(quietly_log)
  
  kept_rows <- res %>%
    tibble::column_to_rownames("term") %>%
    dplyr::select(grep("pVal", names(.))) %>%
    dplyr::mutate(Kept = rowSums(!is.na (.)) > 0) %>%
    dplyr::select(Kept) %>%
    unlist()
  
  qvals <- res[kept_rows, ] %>%
    `rownames<-`(.$term) %>%
    dplyr::select(grep("pVal", names(.))) %>%
    my_padj(pval_cutoff) %>%
    tibble::rownames_to_column("term")
  
  res <- res %>%
    tibble::column_to_rownames("term") %>%
    dplyr::select(-which(names(.) %in% names(qvals))) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::right_join(qvals, by = "term") %T>% 
    write.table(file.path(filepath, paste0(gsub("\\..*$", "", gsub("\\..*$", "", filename)), "_pVals.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
}


#' Plots histograms
#' 
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr purrr rlang mixtools ggplot2 RColorBrewer
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
plotHisto <- function (df = NULL, id, label_scheme_sub, scale_log2r, complete_cases, 
                       show_curves, show_vline, scale_y, filepath = NULL, filename, 
                       theme, ...) {

  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  if (scale_log2r) {
    fn_par <- file.path(filepath, "MGKernel_params_Z.txt")
  } else {
    fn_par <- file.path(filepath, "MGKernel_params_N.txt")
  }
  
  if (file.exists(fn_par)) {
    params <- read.csv(fn_par, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")
    params <- within(params, {mean = mean - x})
  } else {
    params <- NULL
    show_curves <- FALSE
  }

  stopifnot(nrow(label_scheme_sub) > 0)
  stopifnot(rlang::is_logical(scale_log2r))
  stopifnot(rlang::is_logical(show_curves))
  stopifnot(rlang::is_logical(show_vline))
  stopifnot(rlang::is_logical(scale_y))

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

	xmin <- eval(dots$xmin, env = caller_env()) # `xmin = -1` is `language`
	xmax <- eval(dots$xmax, env = caller_env()) 
	xbreaks <- eval(dots$xbreaks, env = caller_env())
	binwidth <- eval(dots$binwidth, env = caller_env())
	alpha <- eval(dots$alpha, env = caller_env())
	ncol <- eval(dots$ncol, env = caller_env())
	width <- eval(dots$width, env = caller_env())
	height <- eval(dots$height, env = caller_env())
	scales <- eval(dots$scales, env = caller_env())

	if (is.null(xmin)) xmin <- -2
	if (is.null(xmax)) xmax <- 2
	if (is.null(xbreaks)) xbreaks <- 1
	if (is.null(binwidth)) binwidth <- (xmax - xmin)/80
	if (is.null(alpha)) alpha <- .8
	if (is.null(ncol)) ncol <- 5
	if (is.null(width)) width <- 4 * ncol + 2
	if (is.null(height)) height <- length(label_scheme_sub$Sample_ID) * 4 / ncol
	if (is.null(scales)) scales <- "fixed"
  
	dots <- dots %>% 
    .[! names(.) %in% c("xmin", "xmax", "xbreaks", "binwidth", "ncol", "alpha", 
                        "width", "height", "scales")]

	filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
	dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
	
	if (scale_y) {
	  df <- df %>% 
	    filters_in_call(!!!filter_dots) %>% 
	    arrangers_in_call(!!!arrange_dots)
	  
	  stopifnot(nrow(df) > 0)
	}
	
	by = (xmax - xmin)/200
	nrow <- nrow(df)
	x_label <- expression("Ratio ("*log[2]*")")
	NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")

	if (!is.null(params)) {
		n_comp <- max(params$Component)
		nm_comps <- paste0("G", 1:n_comp)
		nm_full <- c(nm_comps, paste(nm_comps, collapse = " + "))

		# offset by the percentage of non-NA values
		perc_nna <- df %>%
			dplyr::select(grep(paste0(NorZ_ratios, "[0-9]{3}"), names(.))) %>%
			# `names<-`(gsub(".*_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>%
		  `names<-`(gsub(".*_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
			lapply(function(x) sum(!is.na(x)) / length(x) * nrow * binwidth)

		perc_nna <- perc_nna[names(perc_nna) %in% label_scheme_sub$Sample_ID]

		# density profiles
		fit <- params %>%
			dplyr::filter(.$Sample_ID %in% label_scheme_sub$Sample_ID) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) %>%
			split(.$Sample_ID) %>%
			lapply(sumdnorm, xmin, xmax, by = by)

		fit <- fit %>%
			purrr::map(~ .[, grep("^G[0-9]{1}|^Sum", names(.))]) %>%
			purrr::map2(perc_nna, `*`) %>%
			do.call(rbind, .) %>%
			dplyr::bind_cols(fit %>% do.call(rbind, .) %>% dplyr::select(c("x", "Sample_ID"))) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) %>%
			dplyr::rename(!!sym(paste(nm_comps, collapse = " + ")) := Sum) %>%
			tidyr::gather(key = variable, value = value, -x, -Sample_ID) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
			dplyr::mutate(variable = factor(variable, levels = nm_full)) %>%
			dplyr::arrange(Sample_ID, variable) %>%
			dplyr::filter(Sample_ID %in% label_scheme_sub$Sample_ID)

		myPalette <- c(rep("gray", n_comp), "black")
	} else {
		n_comp <- NA
		nm_comps <- NA
		myPalette <- NA
		nm_full <- NA
		perc_nna <- NA
		fit <- NA
	}


	proteoq_histo_theme <- theme_bw() + theme(
		axis.text.x  = element_text(angle=0, vjust=0.5, size=18),
		axis.ticks.x  = element_blank(), # x-axis ticks
		axis.text.y  = element_text(angle=0, vjust=0.5, size=18),
		axis.title.x = element_text(colour="black", size=24),
		axis.title.y = element_text(colour="black", size=24),
		plot.title = element_text(colour="black", size=24, hjust=.5, vjust=.5),

		strip.text.x = element_text(size = 18, colour = "black", angle = 0),
		strip.text.y = element_text(size = 18, colour = "black", angle = 90),

		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),

		legend.key = element_rect(colour = NA, fill = 'transparent'),
		legend.background = element_rect(colour = NA,  fill = "transparent"),
		legend.title = element_blank(),
		legend.text = element_text(colour="black", size=18),
		legend.text.align = 0,
		legend.box = NULL
	)
	
	if (is.null(theme)) theme <- proteoq_histo_theme

	seq <- c(-Inf, seq(4, 7, .5), Inf)

	if (!scale_y) {
	  df <- df %>% 
	  filters_in_call(!!!filter_dots) %>% 
	  arrangers_in_call(!!!arrange_dots)
	}

	stopifnot(nrow(df) > 0)
	
	df_melt <- df %>%
		dplyr::select(grep(paste0(NorZ_ratios, "[0-9]{3}", "|^N_I[0-9]{3}"), names(.))) %>%
		dplyr::filter(rowSums(!is.na(.[, grepl("[IR][0-9]{3}", names(.))])) > 0) %>%
		dplyr::select(which(not_all_zero(.))) %>%
		dplyr::select(which(colSums(!is.na(.)) > 0)) %>%
		dplyr::mutate(Int_index = log10(rowMeans(.[, grepl("^N_I[0-9]{3}", names(.))], na.rm = TRUE))) %>%
		dplyr::mutate_at(.vars = "Int_index", cut, seq, labels = seq[1:(length(seq)-1)]) %>%
		dplyr::select(-grep("^N_I[0-9]{3}", names(.))) %>%
		# `names<-`(gsub(".*log2_R[0-9]{3}.*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
	  `names<-`(gsub("^[NZ]{1}_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
		tidyr::gather(key = Sample_ID, value = value, -Int_index) %>%
		dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme_sub$Sample_ID)) %>%
		dplyr::arrange(Sample_ID) %>%
		dplyr::filter(!is.na(value), !is.na(Int_index)) %>%
		dplyr::filter(Sample_ID %in% label_scheme_sub$Sample_ID) %>% 
	  dplyr::mutate(value = setHMlims(value, xmin, xmax))

	p <- ggplot() +
		geom_histogram(data = df_melt, aes(x = value, y = ..count.., fill = Int_index),
		               color = "white", alpha = alpha, binwidth = binwidth, size = .1) +
		scale_fill_brewer(palette = "Spectral", direction = -1) +
		labs(title = "", x = x_label, y = expression("Frequency")) +
		scale_x_continuous(limits = c(xmin, xmax), breaks = seq(xmin, xmax, by = xbreaks),
		                   labels = as.character(seq(xmin, xmax, by = xbreaks))) +
		facet_wrap(~ Sample_ID, ncol = ncol, scales = scales) + 
	  theme

	if (show_curves) p <- p + geom_line(data = fit, mapping = aes(x = x, y = value, colour = variable),
	                                   size = .2) +
											scale_colour_manual(values = myPalette, name = "Gaussian",
											                    breaks = c(nm_comps, paste(nm_comps, collapse = " + ")),
											                    labels = nm_full)

	if (show_vline) p <- p + geom_vline(xintercept = 0, size = .25, linetype = "dashed")

	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), 
	             width = width, height = height, limitsize = FALSE, dots)
	do.call(ggsave, gg_args)
}


#'Histogram visualization
#'
#'\code{pepHist} plots the histograms of peptide \code{log2FC}.
#'
#'@rdname prnHist
#'@import purrr
#'@export
pepHist <- function (col_select = NULL, scale_log2r = TRUE, complete_cases = FALSE, 
                     show_curves = TRUE, show_vline = TRUE, scale_y = TRUE, 
                     df = NULL, filepath = NULL, filename = NULL, theme = NULL, ...) {
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  stopifnot(rlang::is_logical(show_curves), 
            rlang::is_logical(show_vline), 
            rlang::is_logical(scale_y))
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  dots <- rlang::enexprs(...)
  if (!is.null(dots$impute_na)) {
    dots$impute_na <- NULL
    rlang::warn("No NA imputation with histograms.")
  }

  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = FALSE,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Histogram")(show_curves = show_curves,
                                     show_vline = show_vline, scale_y = scale_y, 
                                     theme = theme, !!!dots)
}


#'Histogram visualization
#'
#'\code{prnHist} plots the histograms of protein \code{log2FC}.
#'
#'In the histograms, the \code{log2FC} under each TMT channel are color-coded by
#'their contributing reporter-ion intensity.
#'
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2FC} to the same scale
#'  of standard deviation across all samples. The default is TRUE.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used. The default is FALSE.
#'@param show_curves Logical; if TRUE, shows the fitted curves. At the TRUE
#'  default, the curve parameters are based on the latest call to
#'  \code{\link{standPep}} or \code{\link{standPrn}} with \code{method_align =
#'  MGKernel}. This feature can inform the effects of data filtration on the
#'  alignment of \code{logFC} profiles. Also see \code{\link{standPep}} and
#'  \code{\link{standPrn}} for more examples.
#'@param show_vline Logical; if TRUE, shows the vertical lines at \code{x = 0}.
#'  The default is TRUE.
#'@param scale_y Logical; if TRUE, scale data on the \code{y-axis}. The default
#'  is TRUE.
#'@param df The name of a primary data file. By default, it will be determined
#'  automatically after matching the types of data and analysis with an
#'  \code{id} among \code{c("pep_seq", "pep_seq_mod", "prot_acc", "gene")}. A
#'  primary file contains normalized peptide or protein data and is among
#'  \code{c("Peptide.txt", "Peptide_pVal.txt", "Peptide_impNA_pVal.txt",
#'  "Protein.txt", "Protein_pVal.txt", "protein_impNA_pVal.txt")}. For analyses
#'  require the fields of significance p-values, the \code{df} will be one of
#'  \code{c("Peptide_pVal.txt", "Peptide_impNA_pVal.txt", "Protein_pVal.txt",
#'  "protein_impNA_pVal.txt")}.
#'@param filepath A file path to output results. By default, it will be
#'  determined automatically by the name of the calling function and the value
#'  of \code{id} in the \code{call}.
#'@param filename A representative file name to outputs. By default, the name(s)
#'  will be determined automatically. For text files, a typical file extension
#'  is \code{.txt}. For image files, they are typically saved via
#'  \code{\link[ggplot2]{ggsave}} or \code{\link[pheatmap]{pheatmap}} where the
#'  image type will be determined by the extension of the file name.
#'@param theme A
#'  \href{https://ggplot2.tidyverse.org/reference/ggtheme.html}{ggplot2}
#'  theme, i.e., theme_bw(), or a custom theme. At the NULL default, a system
#'  theme will be applied.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  of data against the column keys in \code{Peptide.txt} for peptides or
#'  \code{Protein.txt} for proteins. Each statement contains to a list of
#'  logical expression(s). The \code{lhs} needs to start with \code{filter_}.
#'  The logical condition(s) at the \code{rhs} needs to be enclosed in
#'  \code{exprs} with round parenthesis. \cr \cr For example, \code{pep_len} is
#'  a column key present in \code{Mascot} peptide tables of \code{Peptide.txt}.
#'  The statement \code{filter_peps_at = exprs(pep_len <= 50)} will remove
#'  peptide entries with \code{pep_len > 50}. See also \code{\link{normPSM}}.
#'  \cr \cr Additional parameters for plotting with \code{ggplot2}: \cr
#'  \code{xmin}, the minimum \eqn{x} at a log2 scale; the default is -2. \cr
#'  \code{xmax}, the maximum \eqn{x} at a log2 scale; the default is +2. \cr
#'  \code{xbreaks}, the breaks in \eqn{x}-axis at a log2 scale; the default is
#'  1. \cr \code{binwidth}, the binwidth of \code{log2FC}; the default is
#'  \eqn{(xmax - xmin)/80}. \cr \code{ncol}, the number of columns; the default
#'  is 1. \cr \code{width}, the width of plot; \cr \code{height}, the height of
#'  plot. \cr \code{scales}, should the scales be fixed across panels; the
#'  default is "fixed" and the alternative is "free".
#'@inheritParams standPep
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnHist_.R
#'
#'@return Histograms of \code{log2FC}
#'@export
prnHist <- function (col_select = NULL, scale_log2r = TRUE, complete_cases = FALSE, 
                     show_curves = TRUE, show_vline = TRUE, scale_y = TRUE, 
                     df = NULL, filepath = NULL, filename = NULL, theme = NULL, ...) {
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))

  stopifnot(rlang::is_logical(show_curves), 
            rlang::is_logical(show_vline), 
            rlang::is_logical(scale_y))
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  dots <- rlang::enexprs(...)
  if (!is.null(dots$impute_na)) {
    dots$impute_na <- NULL
    rlang::warn("No NA imputation with histograms.")
  }

  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = FALSE,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "Histogram")(show_curves = show_curves,
                                     show_vline = show_vline, scale_y = scale_y, 
                                     theme = theme, !!!dots)
}

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
#' @inheritParams prnHM
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import stringr dplyr rlang ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
plotHM <- function(df, id, col_benchmark, label_scheme_sub, filepath, filename,
                   scale_log2r, complete_cases, 
                   annot_cols = NULL, annot_colnames = NULL, annot_rows = annot_rows, 
                   xmin = -1, xmax = 1, xmargin = .1, 
                   p_dist_rows = 2, p_dist_cols = 2, 
                   hc_method_rows = "complete", hc_method_cols = "complete", ...) {

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  dir.create(file.path(filepath, "Subtrees", fn_prefix), recursive = TRUE, showWarnings = FALSE)
  
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
  
  if (!is.null(dots$clustering_method)) {
    dots$clustering_method <- NULL
    warning("Argument `clustering_method` disabled; 
            use `hc_method_rows` and `hc_method_cols` instead.", 
            call. = FALSE)    
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
  
  x_label <- expression("Ratio ("*log[2]*")")
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  NorZ_ratios_to_ctrl <- paste("toCtrl", NorZ_ratios, sep = "_")
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)

  sample_ids <- label_scheme_sub$Sample_ID
  
  pattern <- "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("
  
  df <- df %>%
    dplyr::mutate_at(vars(grep("^pVal|^adjP", names(.))), as.numeric) %>%
    dplyr::mutate(Mean_log10Int = log10(rowMeans(.[, grepl("^I[0-9]{3}", names(.))],
                                                 na.rm = TRUE))) %>%
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}", names(.))), ~ setHMlims(.x, xmin, xmax)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>% 
    reorderCols(endColIndex = grep(pattern, names(.)), col_to_rn = id) 
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  if (nrow(df) == 0) stop("Zero data rows available after data filtration.", call. = FALSE)

  dfR <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sample_ids)) %>%
    dplyr::select(as.character(sample_ids)) # ensure the same order

  df <- df %>%
    dplyr::select(-grep("log2_R[0-9]{3}", names(.))) %>%
    dplyr::bind_cols(., dfR) %>% 
    dplyr::filter(!is.na(.[[id]])) %>% 
    `rownames<-`(.[[id]])
  
  if (!is.null(dots$annotation_row)) {
    dots$annotation_row <- NULL
    warning("Argument `annotation_row` disabled; use `annot_rows` instead.", call. = FALSE)
  }
  
  if (!is.null(dots$annotation_col)) {
    dots$annotation_col <- NULL
    warning("Argument `annotation_col` disabled; use `annot_cols` instead.", call. = FALSE)
  }

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
    d <- dist(df_hm, method = clustering_distance_rows, p = p_dist_rows)
    d[is.na(d)] <- .5 * max(d, na.rm = TRUE)

    h <- tryCatch(
      hclust(d, hc_method_rows), 
      error = function(e) 1
    )
    
    if (class(h) != "hclust" && h == 1) {
      warning("Row clustering cannot be performed.", call. = FALSE)
      h <- FALSE
    }

    dots$cluster_rows <- h
    rm(d, h)
  } else {
    dots$cluster_rows <- FALSE
  }
  
  if (cluster_cols) {
    d_cols <- dist(t(df_hm), method = clustering_distance_cols, p = p_dist_cols)
    d_cols[is.na(d_cols)] <- .5 * max(d_cols, na.rm = TRUE)

    h_cols <- tryCatch(
      hclust(d_cols, hc_method_cols), 
      error = function(e) 1
    )
    
    if (class(h_cols) != "hclust" && h_cols == 1) {
      warning("Column clustering cannot be performed.", call. = FALSE)
      h_cols <- FALSE
    }
    
    dots$cluster_cols <- h_cols
    # rm(d_cols, h_cols) # h_cols also for subtrees
  } else {
    dots$cluster_cols <- FALSE
  }
  
  filename <- gg_imgname(filename)
  
  # forms `annotation_col` and `annotation_row` from `annot_col` and `annot_row` 
  dots <- dots %>% 
    .[! names(.) %in% c("mat", "filename", "annotation_col", "annotation_row", 
                        "clustering_distance_rows", "clustering_distance_cols", 
                        "clustering_method", 
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
    if (is.numeric(cutree_rows) && nrow(df) >= cutree_rows) {
      Cluster <- data.frame(Cluster = cutree(p$tree_row, k = cutree_rows)) %>%
        dplyr::mutate(!!id := rownames(.)) %>%
        dplyr::left_join(df, by = id) %T>% 
        write.csv(file.path(filepath, "Subtrees", fn_prefix, 
                            paste0(fn_prefix, " n-", cutree_rows, "_subtrees.csv")), row.names = FALSE)

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
          nrow <- nrow(df_sub)
          d_sub <- dist(df_sub, method = clustering_distance_rows, p = p_dist_rows)
          max_d_row <- suppressWarnings(max(d_sub, na.rm = TRUE))
          
          if (length(d_sub) == 0 || is.infinite(max_d_row)) {
            h_sub <- FALSE
          } else {
            d_sub[is.na(d_sub)] <- .5 * max_d_row
            
            if (nrow <= 2) {
              h_sub <- FALSE
            } else {
              h_sub <- tryCatch(
                hclust(d_sub, hc_method_rows),
                error = function(e) 1
              )
              
              if (class(h_sub) != "hclust" && h_sub == 1) {
                warning("No row clustering for subtree: ", cluster_id, call. = FALSE)
                h_sub <- FALSE
              }
            }    
          }
        }
        
        if (cluster_cols) {
          t_df_sub <- t(df_sub)
          
          nrow_trans <- nrow(t_df_sub)
          d_sub_col <- dist(t_df_sub, method = clustering_distance_cols, p = p_dist_cols)
          max_d_col <- suppressWarnings(max(d_sub_col, na.rm = TRUE))
          
          if (length(d_sub_col) == 0) next
          
          if (is.infinite(max_d_col)) {
            v_sub <- FALSE
          } else {
            d_sub_col[is.na(d_sub_col)] <- .5 * max_d_col
  
            if (nrow_trans <= 2) {
              v_sub <- FALSE
            } else {
              v_sub <- tryCatch(
                hclust(d_sub_col, hc_method_cols),
                error = function(e) 1
              )
              
              if (class(v_sub) != "hclust" && v_sub == 1) {
                warning("No column clustering for subtree: ", cluster_id, call. = FALSE)
                v_sub <- FALSE
              }
            }            
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
                                 paste0("Subtree_", cutree_rows, "-", cluster_id, ".", fn_suffix))
          )
        )

        rm(d_sub, h_sub)
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
pepHM <- function (col_select = NULL, col_benchmark = NULL,
                   scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                   df = NULL, filepath = NULL, filename = NULL,
                   annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                   xmin = -1, xmax = 1, xmargin = 0.1, 
                   hc_method_rows = "complete", hc_method_cols = "complete", 
                   p_dist_rows = 2, p_dist_cols = 2, ...) {
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
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
  hc_method_rows <- rlang::as_string(rlang::enexpr(hc_method_rows))
  hc_method_cols <- rlang::as_string(rlang::enexpr(hc_method_cols))

  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, complete_cases = complete_cases,impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename, 
            anal_type = "Heatmap")(xmin = xmin, xmax = xmax, xmargin = xmargin, 
                                   annot_cols = annot_cols, 
                                   annot_colnames = annot_colnames, 
                                   annot_rows = annot_rows, 
                                   p_dist_rows = p_dist_rows, 
                                   p_dist_cols = p_dist_cols, 
                                   hc_method_rows = hc_method_rows, 
                                   hc_method_cols = hc_method_cols, ...)
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
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. Each statement contains
#'  to a list of logical expression(s). The \code{lhs} needs to start with
#'  \code{filter_}. The logical condition(s) at the \code{rhs} needs to be
#'  enclosed in \code{exprs} with round parenthesis. For example, \code{pep_len}
#'  is a column key present in \code{Mascot} peptide tables of
#'  \code{Peptide.txt}. The statement \code{filter_peps_at = exprs(pep_len <=
#'  50)} will remove peptide entries with \code{pep_len > 50}. See also
#'  \code{\link{pepHist}}, \code{\link{normPSM}}. \cr \cr \code{arrange_}:
#'  Variable argument statements for the row ordering against data in a primary
#'  file linked to \code{df}. The \code{lhs} needs to start with
#'  \code{arrange_}. The expression(s) at the \code{rhs} needs to be enclosed in
#'  \code{exprs} with round parenthesis. For example, \code{arrange_peps_by =
#'  exprs(gene, prot_n_pep)} will arrange entries by \code{gene}, then by
#'  \code{prot_n_pep}. \cr \cr Additional parameters for plotting: \cr
#'  \code{width}, the width of plot \cr \code{height}, the height of plot \cr
#'  \cr Additional arguments for \code{\link[pheatmap]{pheatmap}}: \cr 
#'  \code{cluster_rows, clustering_method, clustering_distance_rows}... \cr 
#'  \cr Notes about \code{pheatmap}:
#'  \cr \code{annotation_col} disabled; instead use keys indicated in \code{annot_cols}
#'  \cr \code{annotation_row} disabled; instead use keys indicated in \code{annot_rows}
#'  \cr \code{clustering_method} breaks into \code{hc_method_rows} for row data 
#'  and \code{hc_method_cols} for column data
#'  \cr \code{clustering_distance_rows = "minkowski"} allowed at the powder of \code{p_dist_rows} 
#'  and/or \code{p_dist_cols}
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
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
                   xmin = -1, xmax = 1, xmargin = 0.1, 
                   hc_method_rows = "complete", hc_method_cols = "complete", 
                   p_dist_rows = 2, p_dist_cols = 2, ...) {

  check_dots(c("id", "anal_type", "df2"), ...)
  
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
  hc_method_rows <- rlang::as_string(rlang::enexpr(hc_method_rows))
  hc_method_cols <- rlang::as_string(rlang::enexpr(hc_method_cols))

  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_benchmark = !!col_benchmark,
            scale_log2r = scale_log2r, complete_cases = complete_cases,impute_na = impute_na, 
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
                                   hc_method_cols = hc_method_cols, ...)
}


#' Function factories for informatic analysis
#'
#' \code{info_anal} produces functions for selected informatic analysis.
#'
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param anal_type Character string; the type of analysis that are preset for
#'   method dispatch in function factories. The value will be determined
#'   automatically. Exemplary values include \code{anal_type = c("PCA",
#'   "Corrplot", "EucDist", "GSPA", "Heatmap", "Histogram", "MDS", "Model",
#'   "NMF", "Purge", "Trend", ...)}.
#' @inheritParams prnHist
#' @inheritParams prnHM
#' @inheritParams prnMDS
#' @inheritParams anal_pepNMF
#' @inheritParams gspaMap
#' @inheritParams prnCorr_logFC
#' 
#' @return A function to the given \code{anal_type}.
#' @import dplyr rlang ggplot2 pheatmap openxlsx
#' @importFrom magrittr %>%
info_anal <- function (id = gene, col_select = NULL, col_group = NULL, col_order = NULL,
                       col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                       color_brewer = NULL, fill_brewer = NULL, 
                       size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, col_benchmark = NULL, 
                       scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                       df = NULL, df2 = NULL, filepath = NULL, filename = NULL,
                       anal_type = c("Corrplot", "Heatmap", "Histogram", "MA", "MDS", "Model",
                                     "NMF", "Trend")) {
  
  if (anal_type %in% c("MDS", "Volcano", "mapGSPA")) scipen = 999 else scipen = 0
  
  old_opt <- options(
    scipen = scipen, 
    warn = 0, 
    max.print = 99999
  )
  
  options(
    scipen = scipen,
    warn = 1,
    max.print = 2000000
  )
  
  on.exit(options(old_opt), add = TRUE)
  
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  
  stopifnot(rlang::is_logical(scale_log2r), 
            rlang::is_logical(complete_cases),
            rlang::is_logical(impute_na))
  
  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")
  warn_msg1 <- "Coerce `complete_cases = TRUE` at `impute_na = FALSE`."
  
	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	col_order <- rlang::enexpr(col_order)
	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)
	col_benchmark <- rlang::enexpr(col_benchmark)

	col_select <- ifelse(is.null(col_select), rlang::expr(Select), rlang::sym(col_select))
	col_group <- ifelse(is.null(col_group), rlang::expr(Group), rlang::sym(col_group))
	col_order <- ifelse(is.null(col_order), rlang::expr(Order), rlang::sym(col_order))
	col_color <- ifelse(is.null(col_color), rlang::expr(Color), rlang::sym(col_color))
	col_fill <- ifelse(is.null(col_fill), rlang::expr(Fill), rlang::sym(col_fill))
	col_shape <- ifelse(is.null(col_shape), rlang::expr(Shape), rlang::sym(col_shape))
	col_size <- ifelse(is.null(col_size), rlang::expr(Size), rlang::sym(col_size))
	col_alpha <- ifelse(is.null(col_alpha), rlang::expr(Alpha), rlang::sym(col_alpha))
	col_benchmark <- ifelse(is.null(col_benchmark), rlang::expr(Benchmark), rlang::sym(col_benchmark))
	
	if (col_select == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_group == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_order == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_color == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_fill == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_shape == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_size == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_alpha == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	if (col_benchmark == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
	
	color_brewer <- rlang::enexpr(color_brewer)
	fill_brewer <- rlang::enexpr(fill_brewer)
	if (!is.null(color_brewer)) color_brewer <- rlang::as_string(color_brewer)
	if (!is.null(fill_brewer)) fill_brewer <- rlang::as_string(fill_brewer)
	
	size_manual <- rlang::enexpr(size_manual)
	shape_manual <- rlang::enexpr(shape_manual)
	alpha_manual <- rlang::enexpr(alpha_manual)

	df <- rlang::enexpr(df)
	df2 <- rlang::enexpr(df2)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	load(file = file.path(dat_dir, "label_scheme.rda"))
	
	if (is.null(label_scheme[[col_select]])) {
		stop("Column \'", rlang::as_string(col_select), "\' not found.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
		stop("No samples under column \'", rlang::as_string(col_select), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_group]])) {
		col_group <- rlang::expr(Select)
		warning("Column \'", rlang::as_string(col_group), "\' not found; use column \'Select\'.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_group]])) == 0) {
		col_group <- rlang::expr(Select)
		warning("No samples under \'", rlang::as_string(col_group), "\'; use column \'Select\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_order]])) {
		warning("Column \'", rlang::as_string(col_order), "\' not found; arranged by the alphebatics.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_order]])) == 0) {
	  # warning("No samples under column \'", rlang::as_string(col_order), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_color]])) {
		warning("Column \'", rlang::as_string(col_color), "\' not found.", call. = FALSE)
	} else if (sum(!is.na(label_scheme[[col_color]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_color), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_fill]])) {
		warning("Column \'", rlang::as_string(col_fill), "\' not found.", call. = FALSE)
	} else if(sum(!is.na(label_scheme[[col_fill]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_fill), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_shape]])) {
		warning("Column \'", rlang::as_string(col_shape), "\' not found.")
	} else if(sum(!is.na(label_scheme[[col_shape]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_shape), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_size]])) {
		warning("Column \'", rlang::as_string(col_size), "\' not found.")
	} else if (sum(!is.na(label_scheme[[col_size]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_size), "\'.", call. = FALSE)
	}

	if(is.null(label_scheme[[col_alpha]])) {
		warning("Column \'", rlang::as_string(col_alpha), "\' not found.")
	} else if(sum(!is.na(label_scheme[[col_alpha]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_alpha), "\'.", call. = FALSE)
	}

	if (is.null(label_scheme[[col_benchmark]])) {
		warning("Column \'", rlang::as_string(col_benchmark), "\' not found.")
	} else if (sum(!is.na(label_scheme[[col_benchmark]])) == 0) {
		# warning("No samples under column \'", rlang::as_string(col_benchmark), "\'.", call. = FALSE)
	}

	id <- rlang::as_string(rlang::enexpr(id))
	if (length(id) != 1)
	  stop("\'id\' must be one of \'pep_seq\', \'pep_seq_mod\', \'prot_acc\' or \'gene\'")

	if (id %in% c("prot_acc", "gene")) {
		data_type <- "Protein"
	} else if (id %in% c("pep_seq", "pep_seq_mod")) {
		data_type <- "Peptide"
	} else {
		stop("Unrecognized 'id'; needs to be in c(\"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\")", 
		     call. = TRUE)
	}

	anal_type <- rlang::as_string(rlang::enexpr(anal_type))
	
	if (is.null(filepath)) {
		if (grepl("Trend", anal_type)) {
		  filepath = file.path(dat_dir, data_type, "Trend")
		} else if (grepl("NMF", anal_type)) {
		  filepath = file.path(dat_dir, data_type, "NMF")
		} else if (grepl("GSPA", anal_type)) {
		  filepath = file.path(dat_dir, data_type, "GSPA")
		} else {
		  filepath = file.path(dat_dir, data_type, anal_type)
		}
		dir.create(file.path(filepath, "log"), recursive = TRUE, showWarnings = FALSE)
	} else {
	  stop("Use default `filepath`.", call. = FALSE)
	}

	if (is.null(filename)) {
		fn_prefix <- paste(data_type, anal_type, sep = "_")
		fn_prefix <- paste0(fn_prefix, "_", ifelse(scale_log2r, "Z", "N"))
		fn_prefix <- fn_prefix %>% ifelse(impute_na, paste0(., "_impNA"), .)

		fn_suffix <- ifelse(anal_type == "Model", "txt", "png")
	} else {
	  if (length(filename) > 1) stop("Do not provide multiple file names.", call. = FALSE)
	  fn_prefix <- gsub("\\.[^.]*$", "", filename)
	  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
		if (fn_prefix == fn_suffix) stop("No '.' in the file name.", call. = FALSE)
	  
	  if (anal_type %in% c("Trend", "NMF", "GSPA")) {
	    fn_prefix <- paste(fn_prefix, data_type, anal_type, sep = "_")
	    fn_prefix <- paste0(fn_prefix, "_", ifelse(scale_log2r, "Z", "N"))
	    fn_prefix <- fn_prefix %>% ifelse(impute_na, paste0(., "_impNA"), .)
	  }
	}

	use_pri_data <- c("MDS", "PCA", "EucDist", "Heatmap", "Histogram", "Corrplot", 
	                  "Model", "Volcano", "Trend", "NMF", "NMF_meta", "GSPA", "mapGSPA", 
	                  "GSVA", "GSEA", "String")
	use_sec_data <- c("Trend_line", "NMF_con", "NMF_coef", "NMF_meta", "GSPA_hm", "mapGSPA")

	if (anal_type %in% use_pri_data) {
	  df <- find_pri_df(anal_type = !!anal_type, id = !!id, impute_na = impute_na)
	  if (!is.null(dim(df))) df <- df %>% rm_pval_whitespace()
	} else {
	  df <- NULL
	}
	
	if (anal_type %in% use_sec_data) {
	  df2 <- rlang::eval_bare(df2, env = current_env())
	  vararg_secmsg(id = !!id, anal_type = !!anal_type)
	} else {
	  df2 <- NULL
	}

	if (anal_type %in% c("Model", "GSPA", "GSVA", "GSEA")) {
		label_scheme_sub <- label_scheme %>% # to be subset by "formulas"
			dplyr::filter(!grepl("^Empty\\.[0-9]+", .$Sample_ID), !Reference)
	} else {
		label_scheme_sub <- label_scheme %>%
			dplyr::select(Sample_ID, TMT_Set, !!col_select, !!col_group, !!col_order, !!col_color,
			              !!col_fill, !!col_shape, !!col_size, !!col_alpha, !!col_benchmark) %>%
			dplyr::filter(!is.na(!!col_select))
	}

	if (nrow(label_scheme_sub) == 0)
	  stop(paste0("No samples or conditions were defined for \"", anal_type, "\""))
	
	force(scale_log2r)
	force(impute_na)
	force(complete_cases)
	
	message("scale_log2r = ", scale_log2r)
	message("impute_na = ", impute_na)
	message("complete_cases = ", complete_cases)
	
	## primary functions: 
	# `impute_na` determines the `src_path` for `df`
	# `scale_log2r` determines `N_log2R` or `Z_log2R` columns in `df`
	# `complete_cases` subsets data rows in `df` for samples in `label_scheme_sub`
	#   (`Model` and `GSVA`: `complete_cases` further subjects to lm formulae)
	
	## Secondary analysis
	# `scale_log2r` for matching '_N' or '_Z' in input filenames from the corresponding primary function
	#   (special case of `GSPA`: match to the value in `prnSig`)
	# `impute_na` for matching '[_NZ]_impNA' or '[_NZ]' in input filenames
	# `complete_cases` subsets data rows in primary outputs

	if (anal_type == "MDS") {
		function(adjEucDist = FALSE, classical = TRUE, method = "euclidean", p = 2, 
		         k = 3, show_ids = TRUE, theme = NULL, ...) {
		  plotMDS(df = df, 
		          id = !!id,
		          label_scheme_sub = label_scheme_sub, 
		          adjEucDist = adjEucDist, 
		          classical = classical, 
		          method = method,
		          p = p, 
		          k = k,
		          show_ids = show_ids, 
		          col_color = !!col_color, 
		          col_fill = !!col_fill, 
		          col_shape = !!col_shape, 
		          col_size = !!col_size, 
		          col_alpha = !!col_alpha, 
		          color_brewer = !!color_brewer,
		          fill_brewer = !!fill_brewer, 
		          size_manual = size_manual,
		          shape_manual = shape_manual,
		          alpha_manual = alpha_manual, 
		          scale_log2r = scale_log2r,
		          complete_cases = complete_cases, 
		          filepath = filepath, 
		          filename = paste0(fn_prefix, ".", fn_suffix), 
		          theme = theme,
		          anal_type = anal_type, 
		          ...)
		}
	} else if (anal_type == "PCA") {
		function(type = "obs", show_ids = TRUE, theme = NULL, ...) {
		  plotPCA(df = df, 
		          id = !!id,
		          label_scheme_sub = label_scheme_sub, 
		          type = type, 
		          show_ids = show_ids, 
		          col_color = !!col_color, 
		          col_fill = !!col_fill, 
		          col_shape = !!col_shape, 
		          col_size = !!col_size, 
		          col_alpha = !!col_alpha, 
		          color_brewer = !!color_brewer,
		          fill_brewer = !!fill_brewer, 
		          size_manual = size_manual,
		          shape_manual = shape_manual,
		          alpha_manual = alpha_manual, 
		          # prop_var = df_pca$prop_var, 
		          scale_log2r = scale_log2r,
		          complete_cases = complete_cases, 
		          impute_na = impute_na, 
		          filepath = filepath, 
		          filename = paste0(fn_prefix, ".", fn_suffix), 
		          theme = theme,
		          anal_type = anal_type, 
		          ...)
		}
	} else if (anal_type == "EucDist") {
		function(adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, ...) {
		  plotEucDist(df = df, 
		              id = !!id,
    		          label_scheme_sub = label_scheme_sub, 
    		          adjEucDist = adjEucDist, 
    		          scale_log2r = scale_log2r,
    		          complete_cases = complete_cases, 
		              annot_cols = annot_cols, 
	                annot_colnames = annot_colnames, 
    		          filepath = filepath, 
    		          filename = paste0(fn_prefix, ".", fn_suffix), 
    		          anal_type = anal_type, 
    		          ...)
		}
	} else if (anal_type == "Heatmap") {
		function(xmin = -1, xmax = 1, xmargin = 0.1,
		         annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
		         p_dist_rows = 2, p_dist_cols = 2, 
		         hc_method_rows = "complete", hc_method_cols = "complete", ...) {
		  plotHM(df = df, 
		         id = !!id, 
		         col_benchmark = !!col_benchmark,
			       label_scheme_sub = label_scheme_sub,
			       filepath = filepath, 
			       filename = paste0(fn_prefix, ".", fn_suffix),
			       scale_log2r = scale_log2r, 
			       complete_cases = complete_cases, 
			       annot_cols = annot_cols, 
			       annot_colnames = annot_colnames, 
			       annot_rows = annot_rows, 
			       xmin = xmin, 
			       xmax = xmax, 
			       xmargin = xmargin,
			       p_dist_rows = p_dist_rows, 
			       p_dist_cols = p_dist_cols, 
			       hc_method_rows = hc_method_rows, 
			       hc_method_cols = hc_method_cols, 
			       ...)
		}
	} else if (anal_type == "Histogram") {
		function(show_curves = TRUE, show_vline = TRUE, scale_y = TRUE, theme = NULL, ...) {
			plotHisto(df = df, 
			          id = !!id, 
			          label_scheme_sub = label_scheme_sub, 
			          scale_log2r = scale_log2r, 
			          complete_cases = complete_cases, 
			          show_curves = show_curves, 
			          show_vline = show_vline, 
			          scale_y = scale_y, 
			          filepath = filepath, 
			          filename = paste0(fn_prefix, ".", fn_suffix), 
			          theme = theme,
			          ...)
		}
	} else if (anal_type == "Corrplot") {
		function(data_select = "logFC", ...) {
		  plotCorr(df = df, 
		           id = !!id, 
		           anal_type = anal_type, 
		           data_select = data_select, 
		           col_select = !!col_select, 
		           col_order = !!col_order,
		           label_scheme_sub = label_scheme_sub, 
		           scale_log2r = scale_log2r, 
		           complete_cases = complete_cases, 
		           filepath = filepath, 
		           filename = paste0(fn_prefix, "_", data_select, ".", fn_suffix), 
		           ...)
		}
	} else if (anal_type == "Model") {
	  function(method = "limma", var_cutoff = 1E-3, pval_cutoff = 1, logFC_cutoff = log2(1), ...) {
	    sigTest(df = df, 
	            id = !!id, 
	            label_scheme_sub = label_scheme_sub,
	            scale_log2r = scale_log2r,
	            complete_cases = complete_cases, 
	            impute_na = impute_na,
	            filepath = filepath, 
	            filename = paste0(fn_prefix, ".", fn_suffix), 
	            method = !!method, 
	            var_cutoff = var_cutoff,
	            pval_cutoff = pval_cutoff, 
	            logFC_cutoff = logFC_cutoff, 
	            data_type = data_type, 
	            anal_type = anal_type,
	            ...)
	  }
	} else if (anal_type == "Volcano") {
	  function(fml_nms = NULL, adjP = FALSE, show_labels = TRUE, theme = NULL, ...) {
	    plotVolcano(df = df, 
	                df2 = NULL, 
	                id = !!id, 
	                adjP = adjP, 
	                show_labels = show_labels, 
	                anal_type = anal_type,
	                gspval_cutoff = 1, 
	                gslogFC_cutoff = 0, 
	                topn = 0, 
	                show_sig = "none", 
	                fml_nms = fml_nms,
	                gset_nms = NULL, 
	                scale_log2r = scale_log2r, 
	                complete_cases = complete_cases, 
	                impute_na = impute_na, 
	                filepath = filepath, 
	                filename = paste0(fn_prefix, ".", fn_suffix), 
	                theme = theme, 
	                ...)
	  }
	} else if (anal_type == "mapGSPA") {
	  function(fml_nms = NULL, adjP = FALSE, show_labels = TRUE, gspval_cutoff = 0.05, 
	           gslogFC_cutoff = log2(1.2), topn = 100, show_sig = "none", gset_nms = "go_sets", 
	           theme = NULL, ...) {
	    plotVolcano(df = df, 
	                df2 = df2,
	                id = !!id, 
	                adjP = adjP, 
	                show_labels = show_labels, 
	                anal_type = anal_type, 
	                gspval_cutoff = gspval_cutoff, 
	                gslogFC_cutoff = gslogFC_cutoff, 
	                topn = topn, 
	                show_sig = show_sig, 
	                fml_nms = fml_nms, 
	                gset_nms = gset_nms, 
	                scale_log2r = scale_log2r, 
	                complete_cases = complete_cases, 
	                impute_na = impute_na, 
	                filepath = filepath, 
	                filename = paste0(fn_prefix, ".", fn_suffix), 
	                theme = theme, 
	                ...)
	  }
	} else if (anal_type == "Trend") {
		function(n_clust = NULL, ...) {
		  analTrend(df = df, 
		            id = !!id, 
                col_group = !!col_group, 
                col_order = !!col_order,
                label_scheme_sub = label_scheme_sub, 
                n_clust = n_clust,
		            scale_log2r = scale_log2r,
		            complete_cases = complete_cases, 
		            impute_na = impute_na,
                filepath = filepath, 
                filename = paste0(fn_prefix %>% paste0(., "_nclust", n_clust), ".txt"), 
		            anal_type = anal_type,
                ...)
		}
	} else if (anal_type == "Trend_line") {
	  function(n_clust = NULL, theme = NULL, ...) {
	    plotTrend(df2 = df2,
	              id = !!id, 
	              col_group = !!col_group, 
	              col_order = !!col_order, 
	              label_scheme_sub = label_scheme_sub, 
	              n_clust = n_clust, 
	              scale_log2r = scale_log2r,
	              complete_cases = complete_cases, 
	              impute_na = impute_na,
	              filepath = filepath, 
	              filename = paste0(fn_prefix, ".", fn_suffix), 
	              theme = theme, 
	              ...)
	  }
	} else if (anal_type == "NMF") {
		function(rank = NULL, nrun = 50, seed = NULL, ...) {
		  analNMF(df = df, 
		          id = !!id, 
		          rank = rank, 
		          nrun = nrun, 
		          seed = seed, 
		          col_group = !!col_group, 
		          label_scheme_sub = label_scheme_sub,
		          scale_log2r = scale_log2r,
		          complete_cases = complete_cases, 
		          impute_na = impute_na,
		          filepath = filepath, 
		          filename = paste0(fn_prefix %>% paste0(., "_rank", rank), ".txt"), 
		          anal_type = anal_type,
		          ...)
		}
	} else if (anal_type == "NMF_con") {
	  function(rank = NULL, ...) {
	    plotNMFCon(df2 = df2, 
	               id = !!id, 
	               rank = rank, 
	               label_scheme_sub = label_scheme_sub, 
	               scale_log2r = scale_log2r, 
	               complete_cases = complete_cases, 
	               impute_na = impute_na,
	               filepath = filepath, 
	               filename = paste0(fn_prefix, ".", fn_suffix), 
	               ...)
	  }
	} else if (anal_type == "NMF_coef") {
	  function(rank = NULL, ...) {
	    plotNMFCoef(df2 = df2, 
	                id = !!id, 
	                rank = rank, 
	                label_scheme_sub = label_scheme_sub, 
	                scale_log2r = scale_log2r,
	                complete_cases = complete_cases,
	                impute_na = impute_na, 
	                filepath = filepath, 
	                filename = paste0(fn_prefix, ".", fn_suffix), 
	                ...)
	  }
	} else if (anal_type == "NMF_meta") {
	  function(rank = NULL, ...) {
	    plotNMFmeta(df = df, 
	                df2 = df2,
	                id = !!id, 
	                rank = rank, 
	                label_scheme_sub = label_scheme_sub, 
	                scale_log2r = scale_log2r,
	                complete_cases = complete_cases, 
	                impute_na = impute_na,
	                filepath = filepath, 
	                filename = paste0(fn_prefix, ".", fn_suffix), 
	                anal_type = anal_type,
	                ...)
	  }
	} else if (anal_type == "GSPA") {
		function(gset_nms = "go_sets", var_cutoff = .5,
		         pval_cutoff = 1E-2, logFC_cutoff = log2(1.1), 
		         gspval_cutoff = 1E-2, gslogFC_cutoff = log2(1), 
		         min_size = 10, max_size = Inf, min_delta = 4, min_greedy_size = 1, 
		         method = "mean", 
		         ...) {
		  gspaTest(df = df, 
		           id = !!id, 
		           label_scheme_sub = label_scheme_sub, 
		           scale_log2r = scale_log2r, 
		           complete_cases = complete_cases, 
		           impute_na = impute_na, 
		           filepath = filepath, 
		           filename = paste0(fn_prefix, ".txt"),
		           gset_nms = gset_nms, 
		           var_cutoff = var_cutoff,
		           pval_cutoff = pval_cutoff, 
		           logFC_cutoff = logFC_cutoff, 
		           gspval_cutoff = gspval_cutoff, 
		           gslogFC_cutoff = gslogFC_cutoff, 
		           min_size = min_size, 
		           max_size = max_size, 
		           min_delta = min_delta,
		           min_greedy_size = min_greedy_size, 
		           method = method, 
		           anal_type = anal_type, 
		           ...)
		}
	} else if (anal_type == "GSPA_hm") {
	  function(...) {
	    gspaHM(df2 = df2,
	           scale_log2r = scale_log2r, 
	           complete_cases = complete_cases, 
	           impute_na = impute_na, 
	           filepath = filepath, 
	           filename = paste0(fn_prefix, ".", fn_suffix), 
	           ...)
	  }
	} else if(anal_type == "GSVA") {
	  function(lm_method = "limma", gset_nms = "go_sets", var_cutoff = .5,
	           pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), ...) {
      gsvaTest(df = df, 
               id = !!id, 
               label_scheme_sub = label_scheme_sub, 
               filepath = filepath, 
               filename = paste0(fn_prefix, ".txt"), 
               scale_log2r = scale_log2r, 
               complete_cases = complete_cases, 
               impute_na = impute_na, 
               gset_nms = gset_nms,
               lm_method = lm_method, 
               var_cutoff = var_cutoff, 
               pval_cutoff = pval_cutoff, 
               logFC_cutoff = logFC_cutoff, 
               anal_type = anal_type, 
               ...)
	  }
	} else if (anal_type == "GSEA") {
	  function(gset_nms = "go_sets", var_cutoff = 0.5, pval_cutoff = 1E-2, logFC_cutoff = log2(1.1), ...) {
	    gspaTest(df = df, 
	            id = !!id, 
	            label_scheme_sub = label_scheme_sub,
	            scale_log2r = scale_log2r, 
	            complete_cases = complete_cases, 
	            impute_na = impute_na,
	            filepath = filepath, 
	            filename = paste0(fn_prefix, ".txt"),
	            gset_nms = gset_nms, 
	            var_cutoff = var_cutoff,
	            pval_cutoff = pval_cutoff, 
	            logFC_cutoff = logFC_cutoff, 
	            gspval_cutoff = 1E-2, # dummy start
	            gslogFC_cutoff = log2(1), 
	            min_size = 10, 
	            max_size = Inf, 
	            min_delta = 1, 
	            min_greedy_size = 1, 
	            method = method, # dummy end
	            anal_type = anal_type, 
	            ...)
	  }
	} else if (anal_type == "String") {
	  function(db_path = "~\\proteoQ\\dbs\\string", score_cutoff = .7, adjP = FALSE, ...) {
	    stringTest(df = df, 
	               id = !!id, 
	               col_group = !!col_group, 
	               col_order = !!col_order,
	               label_scheme_sub = label_scheme_sub, 
	               db_path = db_path,
	               score_cutoff = score_cutoff,
	               adjP = adjP, 
	               scale_log2r = scale_log2r, 
	               complete_cases = complete_cases, 
	               filepath = filepath, 
	               filename = paste0(fn_prefix, ".csv"), 
	               ...)
	  }
	} 
	
}



#' helper for finding input `df`
#' 
#' @param ... Not currently used.
#' @inheritParams info_anal
find_pri_df <- function (anal_type = "Model", df = NULL, id = "gene", impute_na = FALSE, ...) {
  err_msg2 <- "not found. \n Run functions PSM, peptide and protein normalization first."
  err_msg3 <- "not found. \nImpute NA values with `pepImp()` or `prnImp()` or set `impute_na = FALSE`."
  err_msg4 <- "not found at impute_na = TRUE. \nRun `prnSig(impute_na = TRUE)` first."
  err_msg5 <- "not found at impute_na = FALSE. \nRun `prnSig(impute_na = FALSE)` first."

  anal_type <- rlang::as_string(rlang::enexpr(anal_type))
  id <- rlang::as_string(rlang::enexpr(id))

  df <- rlang::enexpr(df)

  if (is.null(df)) {
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
      fn_imp_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_impNA_pVals.txt")
      fn_raw <- file.path(dat_dir, "Peptide", "Peptide.txt")
      fn_imp <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
    } else if (id %in% c("prot_acc", "gene")) {
      fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
      fn_imp_p <- file.path(dat_dir, "Protein\\Model", "Protein_impNA_pVals.txt")
      fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")
      fn_imp <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
    } else {
      stop("Unknown `id`.", call. = FALSE)
    }
    
    if (anal_type %in% c("Histogram", "MA")) { # never impute_na and no pVals
      if (file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, err_msg2), call. = FALSE)
      
      if (id %in% c("pep_seq", "pep_seq_mod")) {
        message("Primary column keys in `Peptide/Peptide.txt` for `filter_` varargs.")
      } else if (id %in% c("prot_acc", "gene")) {
        message("Primary column keys in `Protein/Protein.txt` for `filter_` varargs.")
      } 
    } else if (anal_type %in% c("Model")) { # optional impute_na but no pVals
      if (impute_na) {
        if (file.exists(fn_imp)) src_path <- fn_imp else stop(paste(fn_imp, err_msg3), call. = FALSE)
      } else {
        if (file.exists(fn_raw)) src_path <- fn_raw else stop(paste(fn_raw, err_msg2), call. = FALSE)
      }
      
      if (id %in% c("pep_seq", "pep_seq_mod")) {
        message("Primary column keys in `Peptide/Peptide[_impNA].txt` for `filter_` varargs.")
      } else if (id %in% c("prot_acc", "gene")) {
        message("Primary column keys in `Protein/Protein[_impNA].txt` for `filter_` varargs.")
      }
    } else if (anal_type %in% c("Volcano", "GSPA", "mapGSPA", "GSEA", "String")) { # always use data with pVals
      if (impute_na) {
        if (file.exists(fn_imp_p)) src_path <- fn_imp_p else stop(paste(fn_imp_p, err_msg4), call. = FALSE)
      } else {
        if (file.exists(fn_p)) src_path <- fn_p else stop(paste(fn_p, err_msg5), call. = FALSE)
      }
      
      if (id %in% c("pep_seq", "pep_seq_mod")) {
        message("Primary column keys in `Model/Peptide[_impNA]_pVals.txt` for `filter_` varargs.")
      } else if (id %in% c("prot_acc", "gene")) {
        message("Primary column keys in `Model/Protein[_impNA]_pVals.txt` for `filter_` varargs.")
      }
    } else if (anal_type %in% c("Heatmap", "MDS", "PCA", "EucDist", "Trend", "NMF", "NMF_meta", 
                                "GSVA", "Corrplot")) { # optional impute_na and possible p_vals
      if (impute_na) {
        if (file.exists(fn_imp_p)) src_path <- fn_imp_p else if (file.exists(fn_imp)) src_path <- fn_imp else 
          stop(paste(fn_imp, err_msg3), call. = FALSE)
      } else {
        if (file.exists(fn_p)) src_path <- fn_p else if (file.exists(fn_raw)) src_path <- fn_raw else 
          stop(paste(fn_raw, err_msg2), call. = FALSE)
      }
      
      if (id %in% c("pep_seq", "pep_seq_mod")) {
        message("Primary column keys in `Model/Peptide[_impNA_pVals].txt` for `filter_` varargs.")
      } else if (id %in% c("prot_acc", "gene")) {
        message("Primary column keys in `Model/Protein[_impNA_pVals].txt` for `filter_` varargs.")
      }
    } 
  } else {
    df <- rlang::as_string(df)
    
    if (anal_type == "Model") stop("Use default file name.", call. = FALSE)
    
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      src_path <- file.path(dat_dir, "Peptide\\Model", df)
      if (!file.exists(src_path)) src_path <- file.path(dat_dir, "Peptide", df)
    } else if (id %in% c("prot_acc", "gene")) {
      src_path <- file.path(dat_dir, "Protein\\Model", df)
      if (!file.exists(src_path)) src_path <- file.path(dat_dir, "Protein", df)
    }
  }
  
  df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
                          comment.char = "#"), error = function(e) NA)

  if (is.null(dim(df))) stop(src_path, " not found.", call. = FALSE)
  
  df <- df %>% reorderCols2()

  message(paste("Primary file loaded:", gsub("\\\\", "/", src_path)))
  
  return(df)
}


#' helper for finding input `df` (not currently used)
#' 
#' @param ... Not currently used.
#' @inheritParams info_anal
find_sec_df <- function (df = NULL, anal_type = NULL, id = NULL, ...) {
  df <- rlang::enexpr(df)
  anal_type <- rlang::enexpr(anal_type)
  id <- rlang::enexpr(id)

  if (is.null(df) || is.null(anal_type) || is.null(id)) return (NULL)
  
  df <- rlang::as_string(df)
  anal_type <- rlang::as_string(anal_type)
  id <- rlang::as_string(id)

  new_anal_type <- anal_type %>% 
    gsub("^Trend_.*", "Trend", .) %>% 
    gsub("^NMF_.*", "NMF", .) %>% 
    gsub("^GSPA_.*", "GSPA", .)
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    src_path <- file.path(dat_dir, "Peptide", new_anal_type, df)
  } else if (id %in% c("prot_acc", "gene")) {
    src_path <- file.path(dat_dir, "Protein", new_anal_type, df)
  } else {
    stop("Unknown `id`", call. = FALSE)
  }
  
  df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
                          comment.char = "#"), error = function(e) NA)
  
  if (is.null(dim(df))) stop(src_path, " not found.", call. = FALSE)
  
  return(df)
}


#' helper for finding input `df`
#' 
#' @param ... Not currently used.
#' @inheritParams info_anal
vararg_secmsg <- function (id = NULL, anal_type = NULL, ...) {
  id <- rlang::as_string(rlang::enexpr(id))
  anal_type <- rlang::as_string(rlang::enexpr(anal_type))

  if (id %in% c("pep_seq", "pep_seq_mod")) {
    if (anal_type == "Trend_line") {
      message("Secondary column keys in `Trend/[...]Peptide_Trend_{NZ}[_impNA][...].txt` for `filter2_` varargs.")
    } else if (anal_type == "NMF_con") {
      message("Secondary column keys in `NMF/[...]Peptide_NMF[...]_consensus.txt` for `filter2_` varargs.")
    } else if (anal_type == "NMF_coef") {
      message("Secondary column keys in `NMF/[...]Peptide_NMF[...]_coef.txt` for `filter2_` varargs.")
    } 
  } else if (id %in% c("prot_acc", "gene")) {
    if (anal_type == "mapGSPA") {
      message("Secondary column keys in `GSPA/[...]Protein_GSPA_{NZ}[_impNA].txt` for `filter2_` varargs.")
    } else if (anal_type == "GSPA_hm") {
      message("Secondary column keys in `GSPA/[...]Protein_GSPA_{NZ}_essmap.txt` for `filter2_` varargs.")
      message("Column keys in `GSPA/[...]Protein_GSPA_{NZ}_essmeta.txt` for heat map annotation.")
    } else if (anal_type == "Trend_line") {
      message("Secondary column keys in `Trend/[...]Protein_Trend_{NZ}[_impNA...].txt` for `filter2_` varargs.")
    } else if (anal_type == "NMF_con") {
      message("Secondary column keys in `NMF/[...]Protein_NMF[...]_consensus.txt` for `filter2_` varargs.")
    } else if (anal_type == "NMF_coef") {
      message("Secondary column keys in `NMF/[...]Protein_NMF[...]_coef.txt` for `filter2_` varargs.")
    } 
  } else {
    stop("Unknown `id`", call. = FALSE)
  }
}

#' Plots MDS
#' 
#' @inheritParams prnMDS
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotMDS <- function (df = NULL, id = NULL, label_scheme_sub = NULL, 
                     adjEucDist = FALSE, classical = TRUE, method = "euclidean", p = 2, 
                     k = 3, show_ids = FALSE, 
                     col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                     color_brewer = NULL, fill_brewer = NULL, 
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                     scale_log2r = TRUE, complete_cases = FALSE,
                     filepath = NULL, filename = NULL, theme = NULL, anal_type = "MDS", 
                     ...) {

  stopifnot(nrow(label_scheme_sub) > 0)
  
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    scoreMDS(
      id = !!id, 
      label_scheme_sub = label_scheme_sub, 
      anal_type = anal_type, 
      scale_log2r = scale_log2r, 
      adjEucDist = adjEucDist, 
      classical = classical, 
      method = method, 
      p = p, 
      k = k, ) %>% 
    cmbn_meta(label_scheme_sub) %T>% 
    write.table(file.path(filepath, paste0(fn_prefix, "_res.txt")), sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE)	

	col_fill <- rlang::enexpr(col_fill)
	col_color <- rlang::enexpr(col_color)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)

	map_color <- map_fill <- map_shape <- map_size <- map_alpha <- NA

	if (col_color != rlang::expr(Color) | !rlang::as_string(sym(col_color)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_color))), "X")
	if (col_fill != rlang::expr(Fill)  | !rlang::as_string(sym(col_fill)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_fill))), "X")
	if (col_shape != rlang::expr(Shape) | !rlang::as_string(sym(col_shape)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_shape))), "X")
	if (col_size != rlang::expr(Size) | !rlang::as_string(sym(col_size)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_size))), "X")
	if (col_alpha != rlang::expr(Alpha) | !rlang::as_string(sym(col_alpha)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_alpha))), "X")
	
	if (!is.na(map_color)) col_color <- NULL
	if (!is.na(map_fill)) col_fill <- NULL
	if (!is.na(map_shape)) col_shape <- NULL
	if (!is.na(map_size)) col_size <- NULL
	if (!is.na(map_alpha)) col_alpha <- NULL
	
	rm(map_color, map_fill, map_shape, map_size, map_alpha)
	
	color_brewer <- rlang::enexpr(color_brewer)
	fill_brewer <- rlang::enexpr(fill_brewer)
	if (!is.null(color_brewer)) color_brewer <- rlang::as_string(color_brewer)
	if (!is.null(fill_brewer)) fill_brewer <- rlang::as_string(fill_brewer)
	
	size_manual <- eval_bare(size_manual, env = caller_env())
	shape_manual <- eval_bare(shape_manual, env = caller_env())
	alpha_manual <- eval_bare(alpha_manual, env = caller_env())
	
	proteoq_mds_theme <- theme_bw() + theme(
		 axis.text.x  = element_text(angle=0, vjust=0.5, size=16),
		 axis.text.y  = element_text(angle=0, vjust=0.5, size=16),
		 axis.title.x = element_text(colour="black", size=18),
		 axis.title.y = element_text(colour="black", size=18),
		 plot.title = element_text(face="bold", colour="black", size=20, hjust=0.5, vjust=0.5),

		 panel.grid.major.x = element_blank(),
		 panel.grid.minor.x = element_blank(),
		 panel.grid.major.y = element_blank(),
		 panel.grid.minor.y = element_blank(),

		 legend.key = element_rect(colour = NA, fill = 'transparent'),
		 legend.background = element_rect(colour = NA,  fill = "transparent"),
		 legend.title = element_blank(),
		 legend.text = element_text(colour="black", size=14),
		 legend.text.align = 0,
		 legend.box = NULL
	)
	
	if (is.null(theme)) theme <- proteoq_mds_theme

	mapping <- ggplot2::aes(x = Coordinate.1, y = Coordinate.2,
	                        colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
	                        size = !!col_size, alpha = !!col_alpha)

	idx <- purrr::map(mapping, `[[`, 1) %>% 
	  purrr::map_lgl(is.null)
	
	mapping_var <- mapping[!idx]
	mapping_fix <- mapping[idx]

	fix_args <- list(colour = "darkgray", fill = NA, shape = 21, size = 4, alpha = 0.9) %>% 
	  .[names(.) %in% names(mapping_fix)] %>% 
	  .[!is.na(.)]
	fix_args$stroke <- 0.02

	p <- ggplot() + rlang::eval_tidy(rlang::quo(geom_point(data = df, mapping = mapping_var, !!!fix_args)))

	if (!is.null(fill_brewer)) p <- p + scale_color_brewer(palette = fill_brewer)
	if (!is.null(color_brewer)) p <- p + scale_color_brewer(palette = color_brewer)
	
	if ((!is.null(col_size)) & (!is.null(size_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_size]])) == length(size_manual))
	  p <- p + scale_size_manual(values = size_manual)
	}
	
	if ((!is.null(col_shape)) & (!is.null(shape_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_shape]])) == length(shape_manual))
	  p <- p + scale_shape_manual(values = shape_manual)
	}
	
	if ((!is.null(col_alpha)) & (!is.null(alpha_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_alpha]])) == length(alpha_manual))
	  p <- p + scale_shape_manual(values = alpha_manual)
	}
	
	p <- p +
		labs(title = "", x = expression("Coordinate 1"), y = expression("Coordinate 2")) +
		coord_fixed() +
	  theme

	if (show_ids) {
	  p <- p +
	    geom_text(data = df, 
	              mapping = aes(x = Coordinate.1, y = Coordinate.2, 
	                            label = df$Sample_ID), color = "gray", size = 3)	  
	}

	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), dots)
	do.call(ggsave, gg_args)
}


#' Plots EucDist
#'
#' @inheritParams prnEucDist
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 rlang pheatmap
#' @importFrom magrittr %>%
plotEucDist <- function (df = NULL, id = NULL, label_scheme_sub = NULL, adjEucDist = FALSE, 
                         scale_log2r = TRUE, complete_cases = FALSE, 
                         annot_cols, annot_colnames, 
                         filepath = NULL, filename = NULL, anal_type = "EucDist", 
                         ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)

  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]

  D <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    scoreEucDist(
      id = !!id, 
      label_scheme_sub = label_scheme_sub, 
      anal_type = anal_type, 
      scale_log2r = scale_log2r, 
      adjEucDist = adjEucDist, ) 

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

	stopifnot(is.matrix(D))

	n_color <- 500
	xmin <- 0
	xmax <- ceiling(max(as.vector(D)))
	x_margin <- xmax/10
	color_breaks <- c(seq(xmin, x_margin, length = n_color/2)[1 : (n_color/2-1)],
	                  seq(x_margin, xmax, length = n_color/2)[2 : (n_color/2)])

	if (is.null(dots$color)) {
	  mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
	} else {
	  mypalette <- eval(dots$color, env = caller_env())
	}

	if (is.null(annot_cols)) {
	  annotation_col <- NA
	} else {
	  annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = rownames(D))
	  idx <- which(annot_cols %in% colnames(annotation_col))
	  
	  annot_cols <- annot_cols[idx]
	  annot_colnames <- annot_colnames[idx]
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

	nm_idx <- names(dots) %in% c("mat", "filename", "annotation_col", "color",
	                             "annotation_colors", "breaks")
	dots[nm_idx] <- NULL

	load(file = file.path(dat_dir, "label_scheme.rda"))
	n_TMT_sets <- n_TMT_sets(label_scheme_sub)
	max_width <- 77

	if (is.null(dots$width)) dots$width <- pmin(10*n_TMT_sets*1.2, max_width)

	if (dots$width >= max_width) {
		cat("The width for the graphic device is", dots$width, "inches or more.\n")
		stop("Please consider a a smaller `cellwidth`.")
	}

	filename <- gg_imgname(filename)
	
	my_pheatmap(
		mat = D,
		filename = file.path(filepath, filename),
		annotation_col = annotation_col,
		annotation_row = NA, 
		color = mypalette,
		annotation_colors = annotation_colors,
		breaks = color_breaks,
		!!!dots
	)
}


#' Plots PCA
#'
#' @inheritParams prnPCA
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 rlang
#' @importFrom magrittr %>%
plotPCA <- function (df = NULL, id = NULL, label_scheme_sub = NULL, type = "obs", show_ids = TRUE, 
                     col_color = NULL, col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                     color_brewer = NULL, fill_brewer = NULL, 
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                     # prop_var, 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     filepath = NULL, filename = NULL, theme = NULL, 
                     anal_type = "PCA", ...) {

  stopifnot(nrow(label_scheme_sub) > 0)
  
  complete_cases <- to_complete_cases(complete_cases = complete_cases, impute_na = impute_na)
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]

  res <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    scorePCA(
      id = !!id, 
      label_scheme_sub = label_scheme_sub, 
      anal_type = anal_type, 
      scale_log2r = scale_log2r, 
      type = type, )
      
  df <- res$PCA
  df <- df %>% cmbn_meta(label_scheme_sub)
  prop_var <- res$prop_var
  rm(res)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

	col_color <- rlang::enexpr(col_color)
	col_fill <- rlang::enexpr(col_fill)
	col_shape <- rlang::enexpr(col_shape)
	col_size <- rlang::enexpr(col_size)
	col_alpha <- rlang::enexpr(col_alpha)

	map_color <- map_fill <- map_shape <- map_size <- map_alpha <- NA
	
	if (col_color != rlang::expr(Color) | !rlang::as_string(sym(col_color)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_color))), "X")
	if (col_fill != rlang::expr(Fill)  | !rlang::as_string(sym(col_fill)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_fill))), "X")
	if (col_shape != rlang::expr(Shape) | !rlang::as_string(sym(col_shape)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_shape))), "X")
	if (col_size != rlang::expr(Size) | !rlang::as_string(sym(col_size)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_size))), "X")
	if (col_alpha != rlang::expr(Alpha) | !rlang::as_string(sym(col_alpha)) %in% names(df)) 
	  assign(paste0("map_", tolower(rlang::as_string(col_alpha))), "X")
	
	if (!is.na(map_color)) col_color <- NULL
	if (!is.na(map_fill)) col_fill <- NULL
	if (!is.na(map_shape)) col_shape <- NULL
	if (!is.na(map_size)) col_size <- NULL
	if (!is.na(map_alpha)) col_alpha <- NULL
	
	rm(map_color, map_fill, map_shape, map_size, map_alpha)
	
	color_brewer <- rlang::enexpr(color_brewer)
	fill_brewer <- rlang::enexpr(fill_brewer)
	if (!is.null(color_brewer)) color_brewer <- rlang::as_string(color_brewer)
	if (!is.null(fill_brewer)) fill_brewer <- rlang::as_string(fill_brewer)
	
	size_manual <- eval_bare(size_manual, env = caller_env())
	shape_manual <- eval_bare(shape_manual, env = caller_env())
	alpha_manual <- eval_bare(alpha_manual, env = caller_env())
	
	proteoq_pca_theme <- theme_bw() + theme(
		 axis.text.x  = element_text(angle=0, vjust=0.5, size=20),
		 axis.text.y  = element_text(angle=0, vjust=0.5, size=20),
		 axis.title.x = element_text(colour="black", size=20),
		 axis.title.y = element_text(colour="black", size=20),
		 plot.title = element_text(face="bold", colour="black", size=20, hjust=0.5, vjust=0.5),

		 panel.grid.major.x = element_blank(),
		 panel.grid.minor.x = element_blank(),
		 panel.grid.major.y = element_blank(),
		 panel.grid.minor.y = element_blank(),

		 legend.key = element_rect(colour = NA, fill = 'transparent'),
		 legend.background = element_rect(colour = NA,  fill = "transparent"),
		 legend.title = element_blank(),
		 legend.text = element_text(colour="black", size=14),
		 legend.text.align = 0,
		 legend.box = NULL
	)
	if (is.null(theme)) theme <- proteoq_pca_theme

	mapping <- ggplot2::aes(x = Coordinate.1, y = Coordinate.2,
	                        colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
	                        size = !!col_size, alpha = !!col_alpha)

	idx <- purrr::map(mapping, `[[`, 1) %>% 
	  purrr::map_lgl(is.null)
	
	mapping_var <- mapping[!idx]
	mapping_fix <- mapping[idx]
	
	fix_args <- list(colour = "darkgray", fill = NA, shape = 21, size = 4, alpha = 0.9) %>% 
	  .[names(.) %in% names(mapping_fix)] %>% 
	  .[!is.na(.)]
	fix_args$stroke <- 0.02

	p <- ggplot() +
	  rlang::eval_tidy(rlang::quo(geom_point(data = df, mapping = mapping_var, !!!fix_args)))

	if (!is.null(fill_brewer)) p <- p + scale_color_brewer(palette = fill_brewer)
	if (!is.null(color_brewer)) p <- p + scale_color_brewer(palette = color_brewer)
	
	if ((!is.null(col_size)) & (!is.null(size_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_size]])) == length(size_manual))
	  p <- p + scale_size_manual(values = size_manual)
	}
	
	if ((!is.null(col_shape)) & (!is.null(shape_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_shape]])) == length(shape_manual))
	  p <- p + scale_shape_manual(values = shape_manual)
	}
	
	if ((!is.null(col_alpha)) & (!is.null(alpha_manual))) {
	  stopifnot(length(unique(label_scheme_sub[[col_alpha]])) == length(alpha_manual))
	  p <- p + scale_shape_manual(values = alpha_manual)
	}
	
	p <- p +
		labs(title = "", x = paste0("PC1 (", prop_var[1], ")"), y = paste0("PC2 (", prop_var[2], ")")) +
		coord_fixed() +
	theme

	if (show_ids) {
	  p <- p +
	    geom_text(data = df,
	              mapping = aes(x = Coordinate.1, y = Coordinate.2, label = df$Sample_ID),
	              color = "gray", size = 3)
	}

	filename <- gg_imgname(filename)
	gg_args <- c(filename = file.path(filepath, gg_imgname(filename)), dots)
	do.call(ggsave, gg_args)
}


#' Scores MDS
#'
#' @inheritParams prnMDS
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scoreMDS <- function (df, id, label_scheme_sub, anal_type, scale_log2r, 
                      adjEucDist = FALSE, classical, method = "euclidean", p = 2, k = 3, ...) {

	dots <- rlang::enexprs(...)
	id <- rlang::as_string(rlang::enexpr(id))
	
	stopifnot(nrow(df) > 50)

	df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
	             sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
	  .$log2R

	D <- dist(t(df), method = method, p = p, diag = TRUE, upper = TRUE)
	if (anyNA(D)) stop("Distance cannot be calculated for one more sample pairs.")
	D <- as.matrix(D)

	if (adjEucDist && method == "euclidean") {
	  D <- local({
  	  annotation_col <- colAnnot(annot_cols = c("TMT_Set"), sample_ids = attr(D, "Labels"))
  
  		for (i in 1:ncol(D)) {
  			for (j in 1:ncol(D)) {
  				if (annotation_col$TMT_Set[i] != annotation_col$TMT_Set[j])
  				  D[i, j] <- D[i, j]/sqrt(2)
  			}
  		}
  		
  	  return(D)
		})
	}

	if (!classical) {
		df_mds <- data.frame(isoMDS(D, k = k)$points)
	} else {
		df_mds <- data.frame(cmdscale(D, k = k))
	}
	
	df_mds <- df_mds %>%
		`colnames<-`(paste("Coordinate", 1:k, sep = ".")) %>%
		tibble::rownames_to_column("Sample_ID") %>%
		dplyr::select(which(not_all_zero(.))) %>%
		tibble::column_to_rownames(var = "Sample_ID")
}


#' Scores PCA
#'
#' @inheritParams prnPCA
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scorePCA <- function (df, id, label_scheme_sub, anal_type, scale_log2r, type, ...) {
  
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  
  stopifnot(nrow(df) > 50)
  
  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R 

  if (type == "obs") {
    df_t <- df %>% 
      t() %>%
      data.frame(check.names = FALSE) %>%
      bind_cols(label_scheme_sub)
    
    pr_out <- df_t %>%
      dplyr::select(which(colnames(.) %in% rownames(df))) %>%
      prcomp(scale = scale_log2r)
    
    rownames(pr_out$x) <- label_scheme_sub$Sample_ID
    
    prop_var <- summary(pr_out)$importance[2, ] %>% round(., digits = 3) %>% scales::percent()
    
    df_pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      `names<-`(gsub("^PC", "Coordinate\\.", names(.)))    
  } else if (type == "feats") {
    pr_out <- df %>% prcomp(scale = scale_log2r)
    
    df_pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      `names<-`(gsub("^PC", "Coordinate\\.", names(.)))
    
    prop_var <- summary(pr_out)$importance[2, ] %>% round(., digits = 3) %>% scales::percent()    
  } else {
    stop("Unkown `type` for PCA.", call. = FALSE)
  }

  return(list("PCA" = df_pca, "prop_var" = prop_var))
}


#' Scores Euclidean distance
#'
#' @inheritParams prnEucDist
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr rlang
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>%
scoreEucDist <- function (df, id, label_scheme_sub, anal_type, scale_log2r, adjEucDist = FALSE, ...) {

  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  
  stopifnot(nrow(df) > 50)
  
  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r, 
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>% 
    .$log2R
  
  D <- dist(t(df), method = "euclidean", diag = TRUE, upper = TRUE)
  if (anyNA(D)) stop("Distance cannot be calculated for one more sample pairs.")
  
  D <- as.matrix(D)
  
  if (adjEucDist) {
    D <- local({
      annotation_col <- colAnnot(annot_cols = c("TMT_Set"), colnames(D))
      
      for (i in 1:ncol(D)) {
        for (j in 1:ncol(D)) {
          if (annotation_col$TMT_Set[i] != annotation_col$TMT_Set[j])
            D[i, j] <- D[i, j]/sqrt(2)
        }
      }
      
      return(D)
    })
  }
  
  return(D)
}


#'MDS plots
#'
#'\code{pepMDS} visualizes the multidimensional scaling (MDS) of peptide \code{log2FC}.
#'
#'@rdname prnMDS
#'
#'@import purrr
#'@export
pepMDS <- function (col_select = NULL, col_color = NULL, col_fill = NULL,
                    col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjEucDist = FALSE, classical = TRUE, 
                    method = "euclidean", p = 2, k = 3, 
                    show_ids = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(rlang::is_logical(adjEucDist), 
            rlang::is_logical(classical), 
            rlang::is_logical(show_ids), 
            rlang::is_double(p), 
            rlang::is_double(k))
  
  col_select <- rlang::enexpr(col_select)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)
  
  color_brewer <- rlang::enexpr(color_brewer)
  fill_brewer <- rlang::enexpr(fill_brewer)
  size_manual <- rlang::enexpr(size_manual)
  shape_manual <- rlang::enexpr(shape_manual)
  alpha_manual <- rlang::enexpr(alpha_manual)
  
  method <- rlang::as_string(rlang::enexpr(method))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, method = method, 
                               p = p, k = k, show_ids = show_ids,
                               theme = theme, ...)
  
}


#'MDS plots
#'
#'\code{prnMDS} visualizes the multidimensional scaling (MDS) of protein
#'\code{log2FC}.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}}, followed by a metric
#'(\code{\link[stats]{cmdscale}}) or non-metric (\code{\link[MASS]{isoMDS}})
#'MDS. The default is metric MDS with the input dissimilarities being euclidean
#'distances.
#'
#'@inheritParams  prnHist
#'@inheritParams prnHM
#'@param  col_color Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{color} aesthetics in plots. At
#'  the NULL default, the column key \code{Color} will be used.
#'@param  col_fill Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{fill} aesthetics in plots. At
#'  the NULL default, the column key \code{Fill} will be used.
#'@param  col_shape Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{shape} aesthetics in plots. At
#'  the NULL default, the column key \code{Shape} will be used.
#'@param  col_size Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{size} aesthetics in plots. At
#'  the NULL default, the column key \code{Size} will be used.
#'@param  col_alpha Character string to a column key in \code{expt_smry.xlsx}.
#'  Values under which will be used for the \code{alpha} (transparency)
#'  aesthetics in plots. At the NULL default, the column key \code{Alpha} will
#'  be used.
#'@param color_brewer Character string to the name of a color brewer for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_brewer.html}{ggplot2::scale_color_brewer},
#'   i.e., \code{color_brewer = Set1}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param fill_brewer Character string to the name of a color brewer for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_brewer.html}{ggplot2::scale_fill_brewer},
#'   i.e., \code{fill_brewer = Spectral}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param size_manual Numeric vector to the scale of sizes for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_size_manual},
#'   i.e., \code{size_manual = c(8, 12)}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param shape_manual Numeric vector to the scale of shape IDs for use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_shape_manual},
#'   i.e., \code{shape_manual = c(5, 15)}. At the NULL default, the setting in
#'  \code{ggplot2} will be used.
#'@param alpha_manual Numeric vector to the scale of transparency of objects for
#'  use in
#'  \href{https://ggplot2.tidyverse.org/reference/scale_manual.html}{ggplot2::scale_alpha_manual}
#'   , i.e., \code{alpha_manual = c(.5, .9)}. At the NULL default, the setting
#'  in \code{ggplot2} will be used.
#'@param adjEucDist Logical; if TRUE, adjusts the inter-plex \code{Euclidean}
#'  distance by \eqn{1/sqrt(2)} at \code{method = "euclidean"}. The option
#'  \code{adjEucDist = TRUE} may be suitable when \code{reference samples} from
#'  each TMT plex undergo approximately the same sample handling process as the
#'  samples of interest. For instance, \code{reference samples} were split at
#'  the levels of protein lysates. Typically, \code{adjEucDist = FALSE} if
#'  \code{reference samples} were split near the end of a sample handling
#'  process, for instance, at the stages immediately before or after TMT
#'  labeling. Also see online
#'  \href{https://github.com/qzhang503/proteoQ}{README, section MDS} for a brief
#'  reasoning.
#'@param classical Logical. Metric MDS will be performed at TRUE and non-metric
#'  MDS at FALSE (see also \code{\link[stats]{cmdscale}} and
#'  \code{\link[MASS]{isoMDS}}). The default is TRUE.
#'@param method Character string; the distance measure in
#'  \code{\link[stats]{dist}}. The default method is "euclidean".
#'@param p Numeric; The power of the Minkowski distance in
#'  \code{\link[stats]{dist}}. The default is 2.
#'@param k Numeric; The desired dimension for the solution passed to
#'  \code{\link[stats]{cmdscale}}. The default is 3.
#'@param show_ids Logical; if TRUE, shows the sample IDs in \code{MDS/PCA}
#'  plots. The default is TRUE.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional parameters for
#'  \code{ggsave}: \cr \code{width}, the width of plot; \cr \code{height}, the
#'  height of plot \cr \code{...}
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@example inst/extdata/examples/prnMDS_.R
#'
#'@return MDS plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnMDS <- function (col_select = NULL, col_color = NULL, col_fill = NULL,
                    col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjEucDist = FALSE, classical = TRUE, 
                    method = "euclidean", p = 2, k = 3, 
                    show_ids = TRUE, df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(rlang::is_logical(adjEucDist), 
            rlang::is_logical(classical), 
            rlang::is_logical(show_ids), 
            rlang::is_double(p), 
            rlang::is_double(k))

  col_select <- rlang::enexpr(col_select)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)
  
  color_brewer <- rlang::enexpr(color_brewer)
  fill_brewer <- rlang::enexpr(fill_brewer)
  size_manual <- rlang::enexpr(size_manual)
  shape_manual <- rlang::enexpr(shape_manual)
  alpha_manual <- rlang::enexpr(alpha_manual)
  
  method <- rlang::as_string(rlang::enexpr(method))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "MDS")(adjEucDist = adjEucDist, classical = classical, method = method, 
                               p = p, k = k, show_ids = show_ids,
                               theme = theme, ...)
}


#'PCA plots
#'
#'\code{prnPCA} visualizes the principal component analysis (PCA) for peptide
#'data.
#'
#'@rdname prnPCA
#'
#'@import purrr
#'@export
pepPCA <- function (col_select = NULL, col_color = NULL, 
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    show_ids = TRUE, type = "obs", 
                    df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  type <- rlang::as_string(rlang::enexpr(type))
  
  stopifnot(rlang::is_logical(show_ids))
  
  col_select <- rlang::enexpr(col_select)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)
  
  color_brewer <- rlang::enexpr(color_brewer)
  fill_brewer <- rlang::enexpr(fill_brewer)
  size_manual <- rlang::enexpr(size_manual)
  shape_manual <- rlang::enexpr(shape_manual)
  alpha_manual <- rlang::enexpr(alpha_manual)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "PCA")(type = type, show_ids = show_ids, theme = theme, ...)
}


#'PCA plots
#'
#'\code{prnPCA} visualizes the principal component analysis (PCA) for protein
#'data.
#'
#'\code{log2FC} are used in PCA (\code{\link[stats]{prcomp}}).
#'
#'@inheritParams prnHist
#'@inheritParams prnHM
#'@inheritParams prnMDS
#'@param complete_cases Logical; always TRUE for PCA.
#'@param type Character string indicating the type of PCA. At the \code{type =
#'  obs} default, the components are by observations; at \code{type = feats},
#'  the components are by features.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional parameters for
#'  \code{ggsave}: \cr \code{width}, the width of plot; \cr \code{height}, the
#'  height of plot \cr \code{...}
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@example inst/extdata/examples/prnPCA_.R
#'
#'@return PCA plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnPCA <- function (col_select = NULL, col_color = NULL, 
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL, 
                    color_brewer = NULL, fill_brewer = NULL, 
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL, 
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    show_ids = TRUE, type = "obs", 
                    df = NULL, filepath = NULL, filename = NULL, 
                    theme = NULL, ...) {
  check_dots(c("id", "col_group", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  type <- rlang::as_string(rlang::enexpr(type))
  
  stopifnot(rlang::is_logical(show_ids))

  col_select <- rlang::enexpr(col_select)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)
  
  color_brewer <- rlang::enexpr(color_brewer)
  fill_brewer <- rlang::enexpr(fill_brewer)
  size_manual <- rlang::enexpr(size_manual)
  shape_manual <- rlang::enexpr(shape_manual)
  alpha_manual <- rlang::enexpr(alpha_manual)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha, 
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer, 
            size_manual = !!size_manual, shape_manual = !!shape_manual, alpha_manual = !!alpha_manual, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "PCA")(type = type, show_ids = show_ids, theme = theme, ...)
}


#'Distance plots
#'
#'\code{pepEucDist} visualizes the heat map of Euclidean distances for peptide
#'data.
#'
#'@rdname prnEucDist
#'
#'@import purrr
#'@export
pepEucDist <- function (col_select = NULL, 
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                        adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, 
                        df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "col_group", "col_color", "col_fill", 
               "col_shape", "col_size", "col_alpha", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(rlang::is_logical(adjEucDist))
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = NULL, col_fill = NULL, 
            col_shape = NULL, col_size = NULL, col_alpha = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "EucDist")(adjEucDist = adjEucDist, 
                                   annot_cols = annot_cols, annot_colnames = annot_colnames, ...)
}


#'Distance plots
#'
#'\code{prnEucDist} visualizes the heat map of Euclidean distances for protein
#'data.
#'
#'An Euclidean distance matrix of \code{log2FC} is returned by
#'\code{\link[stats]{dist}} for heat map visualization. 
#'
#'@inheritParams prnHist
#'@inheritParams prnMDS
#'@inheritParams prnHM
#'@param annot_cols A character vector of column keys in \code{expt_smry.xlsx}.
#'  The values under the selected keys will be used to color-code sample IDs on
#'  the top of the indicated plot. The default is NULL without column
#'  annotation.
#'@param annot_colnames A character vector of replacement name(s) to
#'  \code{annot_cols}. The default is NULL without name replacement.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional parameters for
#'  plotting: \cr \code{width}, the width of plot \cr \code{height}, the height
#'  of plot \cr 
#'  \cr Additional arguments for \code{\link[pheatmap]{pheatmap}}: \cr 
#'  \code{cluster_rows, clustering_method, clustering_distance_rows}... \cr 
#'  \cr Notes about \code{pheatmap}:
#'  \cr \code{annotation_col} disabled; instead use keys indicated in \code{annot_cols}
#'  \cr \code{annotation_row} disabled; instead use keys indicated in \code{annot_rows}
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@example inst/extdata/examples/prnEucDist_.R
#'@return Heat map visualization of distance matrices.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnEucDist <- function (col_select = NULL, 
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                        adjEucDist = FALSE, annot_cols = NULL, annot_colnames = NULL, 
                        df = NULL, filepath = NULL, filename = NULL, ...) {
  check_dots(c("id", "col_group", "col_color", "col_fill", 
               "col_shape", "col_size", "col_alpha", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  stopifnot(rlang::is_logical(adjEucDist))
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id,
            col_select = !!col_select, col_group = NULL, col_color = NULL, col_fill = NULL, 
            col_shape = NULL, col_size = NULL, col_alpha = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "EucDist")(adjEucDist = adjEucDist, 
                                   annot_cols = annot_cols, annot_colnames = annot_colnames, ...)
}
#' NMF analysis
#'
#' @inheritParams anal_prnNMF
#' @inheritParams info_anal
#' @inheritParams gspaTest
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
    
    res_nmf@consensus %>% 
      data.frame(check.names = FALSE) %>% 
      tibble::rownames_to_column("Sample_ID") %>% 
      readr::write_tsv(file.path(filepath, paste0(.x, "_consensus.txt")))
    
    coef(res_nmf) %>% 
      data.frame(check.names = FALSE) %>% 
      tibble::rownames_to_column("Cluster") %>% 
      readr::write_tsv(file.path(filepath, paste0(.x, "_coef.txt")))
    
    run_scripts <- FALSE
    if (run_scripts) {
      write.table(res_nmf@consensus, file.path(filepath, paste0(.x, "_consensus.txt")), 
                  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      write.table(coef(res_nmf) %>% as.matrix, file.path(filepath, paste0(.x, "_coef.txt")), 
                  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)      
    }

  })
}


#' Plots consensus results from NMF analysis
#'
#' @inheritParams anal_prnNMF
#' @inheritParams plot_prnNMFCoef
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import NMF dplyr purrr rlang cluster Biobase
#' @importFrom magrittr %>%
plotNMFCon <- function(id, rank, label_scheme_sub, scale_log2r, complete_cases, impute_na, 
                       df2, filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  
  dots <- rlang::enexprs(...)
  if (!purrr::is_empty(dots)) {
    if (any(grepl("^filter_", names(dots)))) {
      stop("Primary `filter_` depreciated; use secondary `filter2_`.")
    }      
  }
  filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
  arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
  dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots)]
  
  # find input df2 ---------------------------
  ins <- list.files(path = filepath, pattern = "NMF_[NZ]{1}.*_rank\\d+\\.rda$")
  if (purrr::is_empty(ins)) stop("No inputs under ", filepath, call. = FALSE)
  
  if (is.null(df2)) {
    ins <- ins %>% 
      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} %>% 
      {if (scale_log2r) .[grepl("_NMF_Z", .)] else .[grepl("_NMF_N", .)]}

    if (purrr::is_empty(ins)) 
      stop("No inputs correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, call. = FALSE)
    
    if (is.null(rank)) {
      df2 <- ins
    } else {
      stopifnot(all(rank >= 2) & all(rank %% 1 == 0))
      
      df2 <- local({
        possibles <- ins %>% 
          gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .) %>% 
          as.numeric() %>% 
          `names<-`(ins)
        
        r2 <- rank %>% .[. %in% possibles]
        
        df2 <- possibles %>% 
          .[. %in% r2] %>% 
          names(.)
      })
    }
    
    if (purrr::is_empty(df2)) 
      stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
           " at rank = ", paste0(rank, collapse = ", "), call. = FALSE)    
  } else {
    df2 <- local({
      possibles <- gsub("_consensus\\.txt$", ".rda", df2)
      non_exists <- possibles %>% .[! . %in% ins]
      if (!purrr::is_empty(non_exists)) {
        stop("Missing consensus file(s): ", purrr::reduce(non_exists, paste, sep = ", "), call. = FALSE)
      }
      df2 <- possibles %>% .[. %in% ins]
    })
    if (purrr::is_empty(df2)) stop("File(s) not found under ", filepath, call. = FALSE)
  }

  # prepare output filename ---------------------------
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Peptide_NMF.*", "\\1", .x)
    })
  } else if (id %in% c("prot_acc", "gene")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Protein_NMF.*", "\\1", .x)
    })
  } else {
    stop("Unknown id = ", id, call. = FALSE)
  }

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  # plot data ---------------------------
  purrr::walk2(df2, custom_prefix, ~ {
    load(file = file.path(filepath, .x))
    D_matrix <- res_nmf@consensus

    if (!is.null(dim(D_matrix))) {
      message(paste("File loaded:", gsub("\\\\", "/", file.path(filepath, .x))))
    } else {
      stop(paste("Non-existed file or directory:", gsub("\\\\", "/", file.path(filepath, .x))))
    }

    rank <- gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .x) %>% as.numeric()
    out_nm <- paste0(.y, fn_prefix, "_rank", rank, ".", fn_suffix)

    if (complete_cases) D_matrix <- D_matrix %>% .[complete.cases(.), ]
    
    D_matrix <- data.frame(D_matrix, check.names = FALSE) %>% 
      filters_in_call(!!!filter2_dots) %>% 
      arrangers_in_call(!!!arrange2_dots) %>% 
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
#' @inheritParams anal_prnNMF
#' @inheritParams plot_prnNMFCoef
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import NMF dplyr purrr rlang cluster Biobase
#' @importFrom magrittr %>%
plotNMFCoef <- function(id, rank, label_scheme_sub, scale_log2r, complete_cases, impute_na, 
                        df2, filepath, filename, ...) {
  stopifnot(nrow(label_scheme_sub) > 0)
  sample_ids <- label_scheme_sub$Sample_ID
  id <- rlang::as_string(rlang::enexpr(id))
  
  dots <- rlang::enexprs(...)
  if (!purrr::is_empty(dots)) {
    if (any(grepl("^filter_", names(dots)))) {
      stop("Primary `filter_` depreciated; use secondary `filter2_`.")
    }      
  }
  filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
  arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
  dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots)]
  
  # find input df2 ---------------------------
  ins <- list.files(path = filepath, pattern = "NMF_[NZ]{1}.*_rank\\d+\\.rda$")
  if (purrr::is_empty(ins)) stop("No inputs under ", filepath, call. = FALSE)
  
  if (is.null(df2)) {
    ins <- ins %>% 
      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} %>% 
      {if (scale_log2r) .[grepl("_NMF_Z", .)] else .[grepl("_NMF_N", .)]}

    if (purrr::is_empty(ins)) 
      stop("No inputs correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, call. = FALSE)
    
    if (is.null(rank)) {
      df2 <- ins
    } else {
      stopifnot(all(rank >= 2) & all(rank %% 1 == 0))
      
      df2 <- local({
        possibles <- ins %>% 
          gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .) %>% 
          as.numeric() %>% 
          `names<-`(ins)
        
        r2 <- rank %>% .[. %in% possibles]
        
        df2 <- possibles %>% 
          .[. %in% r2] %>% 
          names(.)
      })
    }
    
    if (purrr::is_empty(df2)) 
      stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
           " at rank = ", paste0(rank, collapse = ", "), call. = FALSE)    
  } else {
    df2 <- local({
      possibles <- gsub("_coef\\.txt$", ".rda", df2)
      non_exists <- possibles %>% .[! . %in% ins]
      if (!purrr::is_empty(non_exists)) {
        stop("Missing coefficient file(s): ", purrr::reduce(non_exists, paste, sep = ", "), call. = FALSE)
      }
      df2 <- possibles %>% .[. %in% ins]      
    })
    if (purrr::is_empty(df2)) stop("File(s) not found under ", filepath, call. = FALSE)
  }

  # prepare output filename ---------------------------
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Peptide_NMF.*", "\\1", .x)
    })
  } else if (id %in% c("prot_acc", "gene")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Protein_NMF.*", "\\1", .x)
    })
  } else {
    stop("Unknown id = ", id, call. = FALSE)
  }

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  # plot data ---------------------------
  purrr::walk2(df2, custom_prefix, ~ {
    load(file = file.path(filepath, .x))
    D_matrix <- coef(res_nmf)
    
    if (!is.null(dim(D_matrix))) {
      message(paste("File loaded:", gsub("\\\\", "/", file.path(filepath, .x))))
    } else {
      stop(paste("Non-existed file or directory:", gsub("\\\\", "/", file.path(filepath, .x))))
    }

    rank <- gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .x) %>% as.numeric()
    out_nm <- paste0(.y, fn_prefix, "_rank", rank, ".", fn_suffix)

    if (complete_cases) D_matrix <- D_matrix %>% .[complete.cases(.), ]
    rownames(D_matrix) <- seq_along(1:nrow(D_matrix))
    
    D_matrix <- data.frame(D_matrix, check.names = FALSE) %>% 
      filters_in_call(!!!filter2_dots) %>% 
      arrangers_in_call(!!!arrange2_dots) %>% 
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
#' @inheritParams anal_prnNMF
#' @inheritParams plot_prnNMFCoef
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMFmeta <- function(id, rank, label_scheme_sub, scale_log2r, complete_cases, impute_na, 
                        df, df2, filepath, filename, anal_type, ...) {
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
  
  # find input df2 ---------------------------
  ins <- list.files(path = filepath, pattern = "_rank\\d+\\.rda$")
  if (purrr::is_empty(ins)) stop("No inputs under ", filepath, call. = FALSE)
  
  if (is.null(df2)) {
    ins <- ins %>% 
      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} %>% 
      {if (scale_log2r) .[grepl("_NMF_Z", .)] else .[grepl("_NMF_N", .)]}
  
    if (purrr::is_empty(ins)) 
      stop("No inputs correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, call. = FALSE)
    
    if (is.null(rank)) {
      df2 <- ins
    } else {
      stopifnot(all(rank >= 2) & all(rank %% 1 == 0))
      
      df2 <- local({
        possibles <- ins %>% 
          gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .) %>% 
          as.numeric() %>% 
          `names<-`(ins)
        
        r2 <- rank %>% .[. %in% possibles]
        
        df2 <- possibles %>% 
          .[. %in% r2] %>% 
          names(.)
      })
    }
    
    if (purrr::is_empty(df2)) 
      stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
           " at rank = ", paste0(rank, collapse = ", "), call. = FALSE)    
  } else {
    df2 <- local({
      possibles <- gsub("_metagene\\.txt$", ".rda", df2)
      non_exists <- possibles %>% .[! . %in% ins]
      if (!purrr::is_empty(non_exists)) {
        stop("Missing metagene file(s): ", purrr::reduce(non_exists, paste, sep = ", "), call. = FALSE)
      }
      df2 <- possibles %>% .[. %in% ins]
    })
    if (purrr::is_empty(df2)) stop("File(s) not found under ", filepath, call. = FALSE)
  }
  
  # prepare output filename ---------------------------
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Peptide_NMF.*", "\\1", .x)
    })
  } else if (id %in% c("prot_acc", "gene")) {
    custom_prefix <- purrr::map_chr(df2, ~ {
      gsub("(.*_{0,1})Protein_NMF.*", "\\1", .x)
    })
  } else {
    stop("Unknown id = ", id, call. = FALSE)
  }
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)  
  
  # plot data ---------------------------
  purrr::walk2(df2, custom_prefix, ~ {
    rank <- gsub(".*_rank(\\d+)[^\\d]*\\.rda$", "\\1", .x) %>% as.numeric()
    dir.create(file.path(filepath, .y, rank), recursive = TRUE, showWarnings = FALSE)
    
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
          filename = file.path(filepath, .y, rank, fn_sub),
          annotation_col = annotation_col,
          annotation_row = NA, 
          color = mypalette,
          annotation_colors = annotation_colors,
          breaks = color_breaks,
          !!!dots
        )
        
        df_op <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ] %>%
          tibble::rownames_to_column(id)
        
        write.table(df_op, file.path(filepath, .y, rank, paste0(fn_prefix, "_", i, ".txt")), 
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
#'@param ... \code{filter_}: Logical expression(s) for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional arguments for
#'  \code{\link[NMF]{nmf}}.
#'@return NMF classification of \code{log2FC} data.
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'@example inst/extdata/examples/prnNMF_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
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
  
  check_dots(c("id", "df2", "anal_type"), ...)
  
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
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename, 
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
  
  check_dots(c("id", "df2", "anal_type"), ...)
  
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
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename, 
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
#'@param ...  \code{filter2_}: Variable argument statements for the row
#'  filtration against data in secondary file(s) of
#'  \code{_NMF[...]_consensus.txt} for consensus plots or
#'  \code{_NMF[...]_coef.txt} for coefficient plots. See also
#'  \code{\link{prnGSPAHM}} for the format of \code{filter2_} statements. \cr
#'  \cr \code{arrange2_}: Variable argument statements for the row ordering
#'  against data in secondary file(s) of \code{_NMF[...]_consensus.txt} for
#'  consensus plots or \code{_NMF[...]_coef.txt} for coefficient plots. See also
#'  \code{\link{prnGSPAHM}} for the format of \code{arrange2_} statements. \cr
#'  \cr Additional arguments for \code{\link[pheatmap]{pheatmap}}
#'@inheritParams prnHist
#'@inheritParams plot_prnTrend
#'@inheritParams  prnEucDist
#'@return Consensus heat maps from NMF classification.
#'@import NMF dplyr rlang readr ggplot2
#'@importFrom magrittr %>%
#'@example inst/extdata/examples/prnNMF_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
plot_pepNMFCon <- function (col_select = NULL, 
                            scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                            df2 = NULL, filename = NULL, 
                            annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)

  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  df2 <- rlang::enexpr(df2)
  
  reload_expts()

  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, df2 = !!df2, filepath = NULL, filename = !!filename, 
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
                            df2 = NULL, filename = NULL, 
                            annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  df2 <- rlang::enexpr(df2)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, df2 = !!df2, filepath = NULL, filename = !!filename, 
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
                             df2 = NULL, filename = NULL, 
                             annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)
  
  dir.create(file.path(dat_dir, "Peptide\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)  
  df2 <- rlang::enexpr(df2)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, df2 = !!df2, filepath = NULL, filename = !!filename, 
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
                             df2 = NULL, filename = NULL, 
                             annot_cols = NULL, annot_colnames = NULL, rank = NULL, ...) {
  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  filename <- rlang::enexpr(filename)
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)
  df2 <- rlang::enexpr(df2)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = NULL, df2 = !!df2, filepath = NULL, filename = !!filename, 
            anal_type = "NMF_coef")(rank = rank, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'Heat maps of metagenes from NMF
#'
#'\code{plot_metaNMF} is a wrapper of \code{\link[pheatmap]{pheatmap}} for the
#'visualization of the metagene heat maps from NMF
#'
#'@inheritParams anal_pepNMF
#'@inheritParams plot_prnTrend
#'@inheritParams  prnEucDist
#'@param rank Numeric vector; the factorization rank(s) in
#'  \code{\link[NMF]{nmf}}. At the NULL default, all available ranks from the
#'  results of \code{\link{anal_pepNMF}} or \code{\link{anal_pepNMF}} will be
#'  used.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{prnHM}}. No \code{filter2_} available for corresponding
#'  secondary file(s). \cr \cr \code{arrange_}: Variable argument statements for
#'  the row ordering against data in a primary file linked to \code{df}. See
#'  also \code{\link{prnHM}} for the format of \code{arrange_} statements. No
#'  \code{arrange2_} available for corresponding secondary file(s). \cr \cr
#'  Additional arguments for \code{\link[pheatmap]{pheatmap}}.
#'
#'@import purrr
#'@example inst/extdata/examples/prnNMF_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
plot_metaNMF <- function (col_select = NULL, 
                          scale_log2r = TRUE, complete_cases = FALSE, impute_na = TRUE,  
                          df = NULL, df2 = NULL, filepath = NULL, filename = NULL, 
                          rank = NULL, annot_cols = NULL, annot_colnames = NULL, ...) {
  check_dots(c("id", "anal_type", "col_group"), ...)
  
  dir.create(file.path(dat_dir, "Protein\\NMF\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))  
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  df <- rlang::enexpr(df)
  df2 <- rlang::enexpr(df2)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  annot_cols <- rlang::enexpr(annot_cols)
  annot_colnames <- rlang::enexpr(annot_colnames)  
  
  reload_expts()
  
  info_anal(id = !!id, col_select = !!col_select, col_group = NULL, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            df = !!df, df2 = !!df2, filepath = !!filepath, filename = !!filename, 
            anal_type = "NMF_meta")(rank = rank, annot_cols = !!annot_cols, annot_colnames = !!annot_colnames, ...)
}


#'Data normalization
#'
#'\code{normMulGau} normalizes \code{log2FC} under the assumption of multi
#'Gaussian kernels.
#'
#'@param df An input data frame
#'@inheritParams prnHist
#'@inheritParams mixtools::normalmixEM
#'@inheritParams standPep
#'@return A data frame.
#'
#'@import dplyr purrr rlang mixtools
#'@importFrom magrittr %>%
normMulGau <- function(df, method_align, n_comp, seed = NULL, range_log2r, range_int, filepath, 
                       col_select = NULL, ...) {

  # check n_comp
  find_n_comp <- function (n_comp, method_align) {
    if (is.null(n_comp)) {
      if (method_align == "MGKernel") {
        n_comp <- ifelse(nrow(df) > 3000, 3L, 2L)
        n_comp <- n_comp %>% as.integer()
      } else if (method_align == "MC") {
        n_comp <- 1L
      } else {
        n_comp <- 1L
      }
    } else {
      stopifnot(n_comp >= 1)
      if (method_align == "MGKernel") stopifnot(n_comp >= 2)
    }
    
    return(n_comp)
  }
  
  
  # find the x at max density
  my_which_max <- function (params, label_scheme) {
    fit <- params %>%
      split(., .$Sample_ID) %>%
      lapply(sumdnorm, xmin = -2, xmax = 2, by = 2/400) %>%
      do.call(rbind, .) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
      dplyr::arrange(Sample_ID)
    
    cf_x <- fit %>%
      dplyr::group_by(Sample_ID) %>% 
      dplyr::mutate(Max = max(Sum, na.rm = TRUE)) %>% 
      dplyr::mutate(Max = ifelse(is.infinite(Max), NA, Max)) %>% 
      dplyr::filter(Sum == Max) %>%
      dplyr::mutate(x = mean(x, na.rm = TRUE)) %>% # tie-breaking
      dplyr::filter(!duplicated(x)) %>% 
      dplyr::select(-Max)
    
    cf_empty <- fit %>%
      dplyr::filter(! Sample_ID %in% cf_x$Sample_ID)
    
    if (nrow(cf_empty) > 0) {
      cf_empty <- cf_empty %>%
        dplyr::group_by(Sample_ID) %>%
        dplyr::filter(!duplicated(Sum)) %>%
        dplyr::mutate(x = 0)
    }
    
    # add back the empty samples
    cf_x <- dplyr::bind_rows(cf_x, cf_empty) %>% 
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
      dplyr::arrange(Sample_ID)
  }
  
 
  # SD for all non-trival samples
  calc_sd_fcts <- function (df, range_log2r, range_int, label_scheme) {
    label_scheme_sd <- label_scheme %>%
      dplyr::filter(!Reference, !grepl("^Empty\\.", Sample_ID))	%>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = (.$Sample_ID)))
    
    SD <- df %>%
      dplyr::select(grep("^N_log2_R|^N_I", names(.))) %>%
      dblTrim(., range_log2r, range_int) %>%
      `names<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.)))

    cf_SD <- SD/mean(SD %>% .[names(.) %in% label_scheme_sd$Sample_ID], na.rm = TRUE)
    cf_SD <- cbind.data.frame(fct = cf_SD, SD) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
      dplyr::arrange(Sample_ID)
  }
  
  
  # sample IDS for used in the current fitting
  find_fit_nms <- function(nm_a, nm_b) {
    ind <- purrr::map(nm_b, ~ grepl(.x, nm_a, fixed = TRUE)) %>% 
      purrr::reduce(`|`)
    nm_a <- nm_a[ind]
  }
  
  
  # compare to prior n_comp value
  ok_file_ncomp <- function(filepath, filename, n_comp) {
    if (file.exists(file.path(filepath, filename))) {
      params <- read.table(file.path(filepath, filename), 
                           check.names = FALSE, header = TRUE, comment.char = "#")
      
      n_comp == dplyr::n_distinct(params$Component)
    } else {
      return(FALSE)
    }
  }
  

  # update df
  update_df <- function (df, label_scheme_fit, cf_x_fit, sd_coefs_fit) {
    nm_log2r_n <- names(df) %>% 
      .[grepl("^N_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
      find_fit_nms(label_scheme_fit$Sample_ID)
    
    nm_int_n <- names(df) %>% 
      .[grepl("^N_I[0-9]{3}[NC]*\\s+\\(", .)] %>% 
      find_fit_nms(label_scheme_fit$Sample_ID)
    
    nm_log2r_z <- names(df) %>% 
      .[grepl("^Z_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
      find_fit_nms(label_scheme_fit$Sample_ID)  
    
    df_z <- mapply(normSD, df[, nm_log2r_n, drop = FALSE], 
                   center = cf_x_fit$x, SD = sd_coefs_fit$fct, SIMPLIFY = FALSE) %>%
      data.frame(check.names = FALSE) %>%
      `colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
      `rownames<-`(rownames(df))    
    
    nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
    df_z[, nan_cols] <- 0
    rm(nan_cols)
    
    if (is_empty(nm_log2r_z)) {
      df <- cbind(df, df_z)
      
      nm_log2r_z <- names(df) %>% 
        .[grepl("^Z_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
        find_fit_nms(label_scheme_fit$Sample_ID)
    } else {
      df[, nm_log2r_z] <- df_z
    }
    
    # aligned log2FC and intensity after the calculation of "df_z"
    df[, nm_log2r_n] <- sweep(df[, nm_log2r_n, drop = FALSE], 2, cf_x_fit$x, "-")
    df[, nm_int_n] <- sweep(df[, nm_int_n, drop = FALSE], 2, 2^cf_x_fit$x, "/")    
    
    env_bind(caller_env(), nm_log2r_z = nm_log2r_z)
    
    return(df)
  }
  
  
  plot_foo <- function () {
    # Pass arguments by row
    args <- params %>% dplyr::select(c("lambda", "mean", "sd")) %>% pmap(list)
    
    # Define the function, i.e., normal density function with the weight of lambda
    wt_dnorm <- function(x, lambda, mean, sd) lambda * length(x) * dnorm(x, mean, sd)
    
    # The wrapper of stat_function(); easier to pass "x"
    stat_dnorm <- function(x, args) stat_function(aes(x), fun = wt_dnorm, n = 100, args = args, size = .2)
    
    # Pass the list of "args" to the wrapper function "stat_dnorm()" for plotting
    # ggplot() + lapply(args, stat_dnorm, x = seq(-2, 2, 0.1))
  }
  
  
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)

	dots <- rlang::enexprs(...)
	slice_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^slice_", names(.))]
	nonslice_dots <- dots %>% .[! . %in% slice_dots]
	n_comp <- find_n_comp(n_comp, method_align)

	if (!purrr::is_empty(nonslice_dots)) {
	  data.frame(nonslice_dots) %>%
			bind_cols(n_comp = n_comp) %>%
			write.table(., file = file.path(filepath, "normalmixEM_pars.txt"),
			            sep = "\t", col.names = TRUE, row.names = FALSE)
	}

	# if different `n_comp` between two `method_align = MGKernel`, force `col_select` to all samples
	# if `n_comp` is given but with `method_align = MC`, ignore difference in n_comp
	
	ok_N_ncomp <- ok_file_ncomp(filepath, "MGKernel_params_N.txt", n_comp)
	ok_Z_ncomp <- ok_file_ncomp(filepath, "MGKernel_params_Z.txt", n_comp)
	
	if (method_align == "MGKernel") {
  	if ((!ok_N_ncomp) & (col_select != rlang::expr(sample_ID))) {
  	  col_select <- rlang::expr(Sample_ID)
  	}
  
  	if ((!ok_Z_ncomp) & (col_select != rlang::expr(sample_ID))) {
  	  col_select <- rlang::expr(Sample_ID)
  	}	  
	}

	load(file = file.path(dat_dir, "label_scheme.rda"))
	label_scheme_fit <- label_scheme %>% .[!is.na(.[[col_select]]), ]
	
	nm_log2r_n <- names(df) %>% 
	  .[grepl("^N_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
	  find_fit_nms(label_scheme_fit$Sample_ID)

	nm_int_n <- names(df) %>% 
	  .[grepl("^N_I[0-9]{3}[NC]*\\s+\\(", .)] %>% 
	  find_fit_nms(label_scheme_fit$Sample_ID)

	if (method_align == "MGKernel") {
	  message("method_align = ", method_align)
	  message("n_comp = ", n_comp)
	  message("col_select = ", col_select)

	  params_sub <- df %>% 
	    filters_in_call(!!!slice_dots) %>% 
	    dplyr::select(nm_log2r_n) %>% 
	    # `names<-`(gsub("^N_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>% 
	    `names<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
	    fitKernelDensity(n_comp = n_comp, seed = seed, !!!nonslice_dots) %>% 
	    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
	    dplyr::arrange(Sample_ID, Component)

    if (!ok_N_ncomp) {
      # earlierly forced `col_select = Sample_ID` if detected different `n_comp` 
      # so if not `ok`, `params_sub` must be for all samples
      params <- params_sub
    } else {
      params <- read.table(file.path(filepath, "MGKernel_params_N.txt"), 
                           check.names = FALSE, header = TRUE, comment.char = "#") %>% 
	      dplyr::select(names(params_sub)) %>% 
	      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
	      dplyr::arrange(Sample_ID, Component)
      
      if (anyNA(params$Sample_ID)) {
        try(unlink(file.path(dat_dir, "Peptide\\Histogram\\MGKernel_params_[NZ].txt")))
        try(unlink(file.path(dat_dir, "Peptide\\Peptide.txt")))
        try(unlink(file.path(dat_dir, "Protein\\Histogram\\MGKernel_params_[NZ].txt")))
        try(unlink(file.path(dat_dir, "Protein\\Protein.txt")))
        stop("`MGKernel_params` files contain Sample_ID(s) not in `expt_smry.xlsx` and are removed; 
             start the normalization over again.")
      }
	    
	    rows <- params$Sample_ID %in% params_sub$Sample_ID
	    params[rows, ] <- params_sub
    }
    
	  # profile widths based on all sample columns and data rows
	  sd_coefs <- calc_sd_fcts(df, range_log2r, range_int, label_scheme)

		cf_x <- my_which_max(params, label_scheme) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) 

		list(params, cf_x, sd_coefs) %>%
			purrr::reduce(left_join, by = "Sample_ID") %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) %>%
			write.table(., file = file.path(filepath, "MGKernel_params_N.txt"), sep = "\t",
			            col.names = TRUE, row.names = FALSE)
		
		sd_coefs_fit <- sd_coefs %>% dplyr::filter(Sample_ID %in% label_scheme_fit$Sample_ID)
		cf_x_fit <- cf_x %>% dplyr::filter(Sample_ID %in% label_scheme_fit$Sample_ID)
		df <- update_df(df, label_scheme_fit, cf_x_fit, sd_coefs_fit)
		
		# separate fits of Z_log2_R for updating curve parameters only
		if (!ok_Z_ncomp) {
		  params_z <- df %>% 
		    filters_in_call(!!!slice_dots) %>% 
		    dplyr::select(nm_log2r_z) %>% 
		    # `names<-`(gsub("^Z_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>% 
		    `names<-`(gsub("^Z_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
		    fitKernelDensity(n_comp = n_comp, seed = seed, !!!nonslice_dots) %>% 
		    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
		    dplyr::arrange(Sample_ID, Component) %>% 
		    dplyr::mutate(x = 0)
		} else {
		  params_z_sub <- df %>% 
		    filters_in_call(!!!slice_dots) %>% 
		    dplyr::select(nm_log2r_z) %>% 
		    `names<-`(gsub("^Z_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
		    fitKernelDensity(n_comp = n_comp, seed, !!!nonslice_dots) %>% 
		    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
		    dplyr::arrange(Sample_ID, Component)
		  
		  params_z <- read.table(file.path(filepath, "MGKernel_params_Z.txt"), 
		                         check.names = FALSE, header = TRUE, comment.char = "#") %>% 
		    dplyr::select(names(params_z_sub)) %>% 
		    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
		    dplyr::arrange(Sample_ID, Component)
		  
		  rows_z <- params_z$Sample_ID %in% params_z_sub$Sample_ID
		  params_z[rows_z, ] <- params_z_sub
		  
		  params_z$x <- 0
		}	
		
		write.table(params_z, file = file.path(filepath, "MGKernel_params_Z.txt"),
		            sep = "\t", col.names = TRUE, row.names = FALSE)

	} else if (method_align == "MC") {
	  message("method_align = ", method_align)
	  message("n_comp = NULL")
	  message("col_select = ", col_select)
	  
	  # profile widths based on all sample columns and data rows
	  sd_coefs <- df %>% calc_sd_fcts(range_log2r, range_int, label_scheme)
	  
	  # initialization: NA for Empty smpls; 0 for the rest
	  x_vals <- df %>%
	    dplyr::select(matches("^N_log2_R[0-9]{3}")) %>%
	    # `colnames<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.))) %>%
	    `colnames<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
	    dplyr::summarise_all(funs(median(., na.rm = TRUE))) %>%
	    unlist() %>%
	    data.frame(x = .) %>%
	    tibble::rownames_to_column("Sample_ID") %>%
	    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
	    dplyr::arrange(Sample_ID) %>% 
	    dplyr::mutate(x = ifelse(is.na(x), NA, 0))
	  
	  x_vals <- local({
	    x_vals_fit <- df %>%
	      filters_in_call(!!!slice_dots) %>% 
	      dplyr::select(matches("^N_log2_R[0-9]{3}")) %>%
	      # `colnames<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.))) %>%
	      `colnames<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
	      dplyr::select(which(names(.) %in% label_scheme_fit$Sample_ID)) %>% 
	      dplyr::summarise_all(funs(median(., na.rm = TRUE))) 
	    
	    if (any(is.na(x_vals_fit))) {
	      data.frame(Sample_ID = names(x_vals_fit)[is.na(x_vals_fit)], `mean(x)` = NA) %>% 
	        print()
	    }

	    x_vals_fit <- x_vals_fit %>%
	      unlist() %>%
	      data.frame(x = .) %>%
	      tibble::rownames_to_column("Sample_ID") %>%
	      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
	      dplyr::arrange(Sample_ID)
	    
	    rows <- x_vals$Sample_ID %in% x_vals_fit$Sample_ID
	    x_vals[rows, ] <- x_vals_fit
	    
	    return(x_vals)
	  })
	  
	  df <- update_df(df, label_scheme, x_vals, sd_coefs)
	} else { 
		# depreciated
	  if(!all(method_align %in% df[["gene"]])) {
			stop(paste(setdiff(method_align, df[["gene"]]),
			           "in 'method_align' not found in gene names."), call. = TRUE)
		}

		cf_x <- df %>%
			dplyr::filter(.[["gene"]] %in% method_align) %>%
			dplyr::select(matches("^N_log2_R[0-9]{3}")) %>%
			`colnames<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
			dplyr::summarise_all(funs(median(., na.rm = TRUE))) %>%
			unlist() %>%
			data.frame(x = .) %>%
			tibble::rownames_to_column("Sample_ID") %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID)

		if (sum(!is.na(cf_x$x)) == 0) {
			stop(paste("Protein ID/Alignment method", method_align, "not found.") , call. = TRUE)
		} else {
			df <- mapply(normSD, df[,grepl("^N_log2_R[0-9]{3}", names(df))],
			             center = cf_x$x, SD = cf_SD$fct, SIMPLIFY = FALSE) %>% # scale to the same SD
					data.frame(check.names = FALSE) %>%
					`colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
					cbind(df, .)

			df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
			  sweep(df[, grepl("^N_log2_R[0-9]{3}", names(df))], 2, cf_x$x, "-")
			df[, grepl("^N_I[0-9]{3}", names(df))] <-
			  sweep(df[, grepl("^N_I[0-9]{3}", names(df))], 2, 2^cf_x$x, "/")
		}
	}

	return(df)
}


#' Data normalization
#'
#' \code{dblTrim} doubly trims the \code{log2FC} and reporter-ion intensity by
#' the given ranges.
#' 
#' @param type_r Character string for recognizing the columns of \code{log2FC}.
#' @param type_int Character string for recognizing the columns of \code{intensity}.
#' @inheritParams info_anal
#' @inheritParams standPep
#' 
#' @return A data frame.
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
dblTrim <- function(df, range_log2r, range_int, type_r = "N_log2_R", type_int = "N_I") {
	df_trim <- df
	
	type_r <- paste0("^", type_r, "[0-9]{3}")
	type_int <- paste0("^", type_int, "[0-9]{3}")
	
	# trim by log2-ratios
	col_r <- grepl(type_r, names(df_trim))
	df_trim[, col_r] <- lapply(df_trim[, col_r], function (x) {
			q_ratio <- quantile(x, probs = range_log2r/100, na.rm = TRUE)
			x[x < q_ratio[1] | x > q_ratio[2]] <- NA
			return(x)
		}
	)

	# trim by intensity
	col_int <- grepl(type_int, names(df_trim))
	df_trim[, col_int] <- lapply(df_trim[, col_int], function (x) {
			q_intensity <- quantile(x, probs = range_int/100, na.rm = TRUE)
			x[x < q_intensity[1] | x > q_intensity[2]] <- NA
			return(x)
		}
	)

	# doubly trim
	df_trim[!is.na(df_trim)] <- 1

	df_trim <- mapply(`*`, df_trim[, grepl(type_r, names(df_trim))],
	                  df_trim[, grepl(type_int, names(df_trim))], SIMPLIFY = FALSE) %>%
		data.frame(check.names = FALSE)

	df_trim[] <- mapply(`*`, df[, grepl(type_r, names(df))] , df_trim, SIMPLIFY = FALSE)

	sapply(df_trim, sd, na.rm = TRUE)
}


#' Data normalization
#'
#' \code{sumdnorm} calculates summed density from \code{normMulGau}.
#' 
#' @param x A numeric vector.
#' @param xmin the miminal x values.
#' @param xmax the maximal x values.
#' @param by the step length.
#' @return A data frame.
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
sumdnorm <- function (x, xmin = -4, xmax = 4, by = xmax/200) {
	wt_dnorm <- function (x, lambda, mean, sd) lambda * dnorm(x, mean = mean, sd = sd)

	args <- purrr::pmap(x[, names(x) %in% c("lambda", "mean", "sd")], list) %>%
		`names<-`(x$Sample_ID)

	nm_comps <- paste0("G", seq_len(length(args)))

	Seq <- seq(xmin, xmax, by = by)

	lapply(args, function(args) rlang::eval_tidy(rlang::quo(wt_dnorm(Seq, !!! args)))) %>%
		do.call(cbind, .) %>%
		data.frame(check.names = FALSE) %>%
		`colnames<-`(nm_comps) %>%
		dplyr::mutate(Sum = rowSums(.)) %>%
		dplyr::mutate(x = Seq) %>% #
		dplyr::mutate(Sample_ID = names(args)[1])
}


#' Data normalization
#'
#' \code{normSD} normalizes the SD of \code{log2FC}.
#' 
#' @param x A numeric vector.
#' @param center The position of \code{x} to be centered at.
#' @param SD The standard deviation that data will be scaled to.
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
#' @export
normSD <- function (x, center = 0, SD = 1) {
	if (sum(is.na(x)) == length(x)) {
		x
	} else if ((sum(is.na(x)) + sum(x == 0, na.rm = TRUE)) == length(x)) {
		x[1:length(x)] <- NaN
	} else if (all(x == 0)) {
		x
	} else {
		x <- (x - center) / SD
	}
	return (x)
}


#' Data normalization
#'
#' \code{fitKernelDensity} calculates the fitted density of \code{log2FC}.
#' 
#' @param df An input data frame.
#' @inheritParams standPep
#' @return A data frame.
#'
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
fitKernelDensity <- function (df, n_comp = 3, seed = NULL, ...) {

	nmix_params <- function (x, n_comp = 3, seed = seed, ...) {
		dots <- rlang::enexprs(...)
		x <- rlang::enexpr(x)
		
		if (!is.null(dots$k)) {
			cat(paste("k =", dots$k, "replaced by", paste("n_comp =", n_comp, "\n")))
			dots$k <- NULL
		}

		if (sum(is.na(x)) == length(x) |
					(sum(is.na(x)) + sum(x == 0, na.rm = TRUE)) == length(x) |
					all(x == 0)) {
			df_par <- data.frame(Component = c(1:n_comp), lambda = rep(NA, n_comp),
									mean = rep(NA, n_comp), sd = rep(NA, n_comp))
		} else {
			x <- x[!is.na(x)]

			stopifnot(n_comp > 1)

			mixEM_call <- rlang::quo(mixtools::normalmixEM(!!x, k = !!n_comp, !!!dots))
			if (!is.null(seed)) set.seed(seed) else set.seed(sample(.Random.seed, 1))
			
			quietly_out <- purrr::quietly(rlang::eval_tidy)(mixEM_call, caller_env())
			x_k2 <- quietly_out$result
			df_par <- data.frame(Component = 1:n_comp, lambda = x_k2$lambda, mean = x_k2$mu, sd = x_k2$sigma)
		}

		return(df_par)
	}

	dots <- rlang::enexprs(...)

	min_n <- local({
	  ok_nan_cols <- df %>% 
      purrr::map_lgl(not_all_nan, na.rm = TRUE)
	  
	  min_n <- df %>% 
      .[, ok_nan_cols, drop = FALSE] %>% 
      .[, not_all_NA(.), drop = FALSE] %>% 
      purrr::map_dbl(., ~ (!is.na(.x)) %>% sum()) %>% 
      min()
	})

	if (min_n < 50) stop("Too few data points for fitting with multiple Gaussian functions.")

	lapply(df, nmix_params, n_comp, seed = seed, !!!dots) %>%
		do.call(rbind, .) %>%
		dplyr::mutate(Sample_ID = rownames(.)) %>%
		dplyr::mutate(Sample_ID = gsub("(.*)\\.\\d+$", "\\1", Sample_ID)) %>%
		dplyr::mutate(Channel = rep(1:(nrow(.)/n_comp), each = n_comp)) %>%
		dplyr::arrange(Channel, -lambda) %>%
		dplyr::mutate(Component = rep(1:n_comp, nrow(.)/n_comp)) %>%
		dplyr::mutate(Height = .$lambda * dnorm(.$mean, mean = .$mean, sd = .$sd))
}

#' Processes the metadata of TMT experiments
#'
#' @inheritParams load_expts
#' @inheritParams prnHist
#' @import dplyr tidyr purrr openxlsx
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
prep_label_scheme <- function(dat_dir, filename) {

	my_channels <- function (x) {
		x <- as.character(x)
		pos <- !grepl("^TMT", x)
		x[pos] <- paste0("TMT-", x[pos])

		return(x)
	}

	replace_na_smpls <- function(x, prefix) {
	  i <- is.na(x)
	  replace(as.character(x), i, paste(prefix, seq_len(sum(i)), sep = "."))
	}
	
	not_trival <- function (x) {
	  ok <- (!is.na(x)) & (x != FALSE) & (x != 0)
	}
	
	
	if (is.null(dat_dir)) dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv),
	                                         error = function(e) 1)

	if (dat_dir == 1) stop("Set up the working directory first.", call. = FALSE)

	if (!file.exists(file.path(dat_dir, filename)))
	  stop(filename, " not found under '", dat_dir, "'.")

	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	if (fn_suffix %in% c("xls", "xlsx")) {
		label_scheme_full <- readxl::read_excel(file.path(dat_dir, filename), sheet = "Setup") %>%
												dplyr::filter(rowSums(!is.na(.)) > 0)
	} else if (fn_suffix == "csv") {
		label_scheme_full <- read.csv(file.path(dat_dir, filename), check.names = TRUE,
		                              header = TRUE, comment.char = "#", na.strings = c("", "NA")) %>%
												dplyr::filter(rowSums(!is.na(.)) > 0)
	} else {
		stop(filename, " needs to be '.xls' or '.xlsx'.")
	}
	
	label_scheme_full <- label_scheme_full %>% 
	  dplyr::mutate(Sample_ID = ifelse(grepl("^Empty\\.[0-9]+", Sample_ID), NA, Sample_ID))
	
	check_tmt126 <- label_scheme_full %>% 
	  dplyr::filter(TMT_Channel == "TMT-126") %>% 
	  dplyr::filter(is.na(TMT_Set)|is.na(LCMS_Injection))
	
	if (nrow(check_tmt126) > 0) {
	  stop("`TMT_Set` and/or `LCMS_Injection` indexes corresponding to `TMT-126` in `expt_smry.xlsx` cannot be empty.", 
	       call. = FALSE)
	}
	rm(check_tmt126)

	must_have <- c("TMT_Channel", "TMT_Set", "LCMS_Injection", "RAW_File",
								"Sample_ID", "Reference")

	missing_cols <- must_have[!must_have %in% names(label_scheme_full)]
	if (length(missing_cols) > 0) {
		purrr::walk(missing_cols, ~ cat(paste0("\'", ., "\' must be present in \'", filename, "\'\n")))
		stop("Not all required columns are present in \'", filename, "\'", call. = TRUE)
	}

	default_names <- c("Select", "Group", "Order", "Fill",  "Color", "Shape", "Size", "Alpha", "Peptide_Yield")

	purrr::walk(as.list(default_names), ~ {
		if (!.x %in% names(label_scheme_full)) {
		  message("Column \'", .x, "\' added to \'", filename, "\'")
			label_scheme_full[[.x]] <<- NA
		}
	}, label_scheme_full)

	# a case of label-free data
	if (dplyr::n_distinct(label_scheme_full$TMT_Channel) == 1) label_scheme_full$TMT_Channel <- NA

	TMT_plex <- TMT_plex(label_scheme_full)
	TMT_levels <- TMT_levels(TMT_plex)

	label_scheme_full <- label_scheme_full %>% 
	  dplyr::mutate_at(vars(c("TMT_Channel")), ~ my_channels(.x)) %>% 
	  dplyr::filter(rowSums(is.na(.)) < ncol(.)) %>% 
	  dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File, ignore.case = TRUE)) %>% 
		dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File, ignore.case = TRUE)) %>% # Bruker
	  dplyr::mutate_at(vars(c("Reference")), ~ not_trival(.x)) %>%
	  dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ as.numeric(.x)) %>%
	  dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ round(.x, digits = 2)) %>%
	  tidyr::fill(one_of("TMT_Set", "LCMS_Injection", "RAW_File")) %>%
	  dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
	  dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)

	# add IDs to unused TMT channels
	label_scheme_empty <- label_scheme_full %>%
		dplyr::select(TMT_Channel, TMT_Set, Sample_ID) %>%
		tidyr::unite(key, TMT_Channel, TMT_Set, remove = TRUE) %>%
		dplyr::filter(!duplicated(key)) %>%
		dplyr::mutate(Sample_ID = replace_na_smpls(Sample_ID, "Empty"))

	label_scheme_full <- label_scheme_full %>%
		tidyr::unite(key, TMT_Channel, TMT_Set, remove = FALSE) %>%
		dplyr::select(-Sample_ID) %>%
		dplyr::left_join(label_scheme_empty, by = "key") %>%
		dplyr::select(-key) %>%
		dplyr::mutate(Sample_ID = factor(Sample_ID))
	  
	label_scheme_full <- dplyr::bind_cols(
	  label_scheme_full %>% dplyr::select(Sample_ID), 
	  label_scheme_full %>% dplyr::select(-Sample_ID))

	rm(label_scheme_empty)

	# check the completeness of TMT_Channel
	check_tmt <- label_scheme_full %>%
		tidyr::complete(TMT_Set, LCMS_Injection, TMT_Channel) %>%
		dplyr::filter(is.na(RAW_File)) %>%
		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
		dplyr::mutate(n = n()) %>%
		dplyr::filter(n != TMT_plex)

	if (nrow(check_tmt) > 0) {
		check_tmt %>%
			dplyr::select(TMT_Set, LCMS_Injection, TMT_Channel) %>%
			print()
	  
	  stop("`", check_tmt$TMT_Channel, "` not found under set ", check_tmt$TMT_Set, 
	          " injection ", check_tmt$LCMS_Injection, ".", 
	          "\n(Use `TMT-131`, instead of `TMT-131N`, for 10-plex experiment(s).)", 
	          call. = FALSE)
	}

	# check the uniqueness of RAW_File per TMT_Set and LCMS_Injection
	check_fn <- label_scheme_full %>%
		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
		dplyr::summarise(count = n_distinct(RAW_File)) %>%
		dplyr::filter(count > 1) %>%
		dplyr::select(-count)

	if (nrow(check_fn) > 0) {
		check_fn %>% print()
		stop("More than one RAW filename in the above combination of TMT sets and LCMS injections.")
	}

	# check the uniqueness of Sample_ID
	check_smpls <- label_scheme_full %>%
		dplyr::group_by(TMT_Set, LCMS_Injection) %>%
		dplyr::summarise(count = n_distinct(Sample_ID)) %>%
		dplyr::filter(count != TMT_plex)

	if (nrow(check_smpls) > 0 & TMT_plex > 0) {
		check_smpls %>% print()
		stop(paste("Need", TMT_plex,
		           "unique samples in the above combination of TMT sets and LCMS injections." ))
	}

	save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.rda"))

	wb <- openxlsx::loadWorkbook(file.path(dat_dir, filename))
	openxlsx::writeData(wb, sheet = "Setup", label_scheme_full)
	openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
	
	simple_label_scheme(dat_dir, label_scheme_full)
}


#' Loads the information of analyte prefractionation
#'
#' @inheritParams load_expts
#' @inheritParams prnHist
#' @import dplyr purrr tidyr openxlsx
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
prep_fraction_scheme <- function(dat_dir, filename) {
	if (is.null(dat_dir))
	  dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)

	if (dat_dir == 1) stop("Set up the working directory first.")

	fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	if (file.exists(file.path(dat_dir, filename))) {
		if (fn_suffix %in% c("xls", "xlsx")) {
			fraction_scheme <- readxl::read_excel(file.path(dat_dir, filename), sheet = "Fractions") %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else if (fn_suffix == "csv") {
			fraction_scheme <- read.csv(file.path(dat_dir, filename), check.names = TRUE, header = TRUE,
			                            comment.char = "#", na.strings = c("", "NA"))  %>%
			  dplyr::filter(rowSums(!is.na(.)) > 0)
		} else {
			stop(filename, " needs to be in a file format of '.xls' or '.xlsx'.")
		}
	  
	  fraction_scheme <- fraction_scheme %>% 
	    tidyr::fill(TMT_Set, LCMS_Injection) %>% 
	    dplyr::group_by(TMT_Set, LCMS_Injection) %>% 
	    dplyr::mutate(Fraction = row_number()) 

	  wb <- openxlsx::loadWorkbook(file.path(dat_dir, filename))
	  openxlsx::writeData(wb, sheet = "Fractions", fraction_scheme)
	  openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
 	} else {
 	  assign(".auto_frac_smry", TRUE, envir = .GlobalEnv)
 	  
 	  # warning: data in a auto-generated `frac_smry.xlsx` will be incorrect 
 	  #   if they were based on wrong information from `expt_smry.xlsx`
 	  load(file = file.path(dat_dir, "label_scheme_full.rda"))
 	  
 	  # in case forget to enter RAW_File names
 	  if (anyNA(label_scheme_full$RAW_File)) stop("Enter RAW file names in the experimental summary file")

		fraction_scheme <- label_scheme_full %>%
			dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
			dplyr::filter(!duplicated(RAW_File)) %>%
			dplyr::group_by(TMT_Set, LCMS_Injection) %>%
			dplyr::mutate(Fraction = row_number())

		wb <- openxlsx::createWorkbook()
		openxlsx::addWorksheet(wb, sheetName = "Fractions")
		openxlsx::writeData(wb, sheet = "Fractions", fraction_scheme)
		openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
 	}
	
	save(fraction_scheme, file = file.path(dat_dir, "fraction_scheme.rda"))
}


#'Loads species-specific Databases
#'
#'A function loads a set of precompiled gene sets of 
#'\href{http://current.geneontology.org/products/pages/downloads.html}{GO}
#'and
#'\href{http://software.broadinstitute.org/gsea/msigdb}{molecular signatures}.
#'@seealso \code{\link{load_expts}} for supported species.
#'
#' @examples
#' \donttest{load_dbs("go_sets", "human")}
#'
#'@param species Character string; the name of a species. 
#'@inheritParams prnGSPA
#'@import dplyr rlang
#'@importFrom magrittr %>%
#'@export
load_dbs <- function (gset_nms = NULL, species = NULL) {
  if (is.null(gset_nms)) stop("`gset_nms` cannot be NULL.", call. = FALSE)
  if (is.null(species)) stop("`species` cannot be NULL.", call. = FALSE)
  
  defaults <- c("go_sets", "kegg_sets", "c2_msig")
  sys_defs <- gset_nms %>% .[. %in% defaults]
  not_sys_defs <- gset_nms %>% .[! . %in% defaults]

  if (!purrr::is_empty(sys_defs)) {
    abbr_sp <- purrr::map_chr(species, sp_lookup)
    filelist <- map(abbr_sp, ~ paste0(sys_defs, "_", .x)) %>% unlist()
    
    data(package = "proteoQ", list = filelist)
    gsets <- purrr::map(filelist, ~ try(get(.x))) %>% do.call(`c`, .)
    
    try(rm(list = filelist, envir = .GlobalEnv))
    
    if (length(gsets) > 0) names(gsets) <- gsub("/", "-", names(gsets))      
  } else {
    gsets <- NULL
  }
  
  if (!purrr::is_empty(not_sys_defs)) {
    if (!all(grepl("\\.rds$", not_sys_defs))) {
      stop("Custom gene set files indicated by `gset_nms` must end with the `.rds` extension.", call. = FALSE)
    }

    not_oks <- not_sys_defs %>% .[!file.exists(not_sys_defs)]
    if (!purrr::is_empty(not_oks)) {
      stop("File not found: \n", reduce(not_oks, paste, sep = ", \n"), call. = FALSE)
    }
    
    gsets2 <- purrr::map(not_sys_defs, readRDS) %>% do.call(`c`, .)
    
    if (length(gsets2) > 0)  {
      names(gsets2) <- gsub("/", "-", names(gsets2))
    } else {
      stop("Empty data file in: \n", reduce(not_sys_defs, paste, sep = ", \n"), call. = FALSE)
    }
  } else {
    gsets2 <- NULL
  }
  
  gsets <- c(gsets, gsets2) %>% .[!duplicated(.)]
  stopifnot(length(gsets) > 0)
  
  assign("gsets", gsets, envir = .GlobalEnv)
} 


#'Load TMT experiments
#'
#'\code{load_expts} processes \code{.xlsx} files containing the metadata of TMT
#'experiments
#'
#'@section \code{expt_smry.xlsx}: The \code{expt_smry.xlsx} files should be
#'  located immediately under the file folder defined by \code{dat_dir}. The tab
#'  containing the metadata of TMT experiments should be named \code{Setup}. The
#'  \code{Excel} spread sheet therein is comprised of three tiers of fields: (1)
#'  essential, (2) optional default and (3) optional open. The \code{essential}
#'  columns contain the mandatory information of TMT experiments. The
#'  \code{optional default} columns serve as the fields for default lookups in
#'  sample selection, grouping, ordering, aesthetics, etc. The \code{optional
#'  open} fields allow users to define their own analysis, aesthetics, etc.
#'
#'  \tabular{ll}{ \strong{Essential column}   \tab \strong{Descrption}\cr
#'  Sample_ID \tab Unique sample IDs \cr TMT_Channel \tab TMT channel names:
#'  \code{126}, \code{127N}, \code{127C} et al. \cr TMT_Set \tab TMT experiment
#'  indexes 1, 2, 3, ... \cr LCMS_Injection   \tab LC/MS injection indexes 1, 2,
#'  3, ... under a \code{TMT_Set} \cr RAW_File \tab MS data file names
#'  originated by \code{Xcalibur} with or without the \code{.raw} extension \cr
#'  Reference \tab Labels indicating reference samples in TMT experiments \cr }
#'
#'  \code{Sample_ID}: values should be unique for entries at a unique
#'  combination of \code{TMT_Channel} and \code{TMT_Set}, or left blank for
#'  unused entries. Samples with the same indexes of \code{TMT_Channel} and
#'  \code{TMT_Set} but different indexes of \code{LCMS_Injection} should have
#'  the same value in \code{Sample_ID}. No white space or special characters are
#'  allowed.
#'
#'  \code{RAW_File}: for analysis with off-line fractionation of peptides
#'  before LC/MS, the \code{RAW_File} column should be left blank. Instead, the
#'  correspondence between the fraction numbers and \code{RAW_File} names should
#'  be specified in a separate file, for example, \code{frac_smry.xlsx}. For
#'  analysis without off-line fractionation, it is recommended as well to leave
#'  the field under the \code{RAW_File} column blank and instead enter the MS
#'  file names in \code{frac_smry.xlsx}.
#'
#'  The set of \code{RAW_File} names in \code{frac_smry.xlsx} needs to be
#'  identical to those in PSM data. Note that \code{OS} file names may be
#'  altered by MS users and thus different to those recorded in \code{Xcalibur}.
#'  The original names by \code{Xcalibur} should be used. MS files may
#'  occasionally have no contributions to PSM findings. These MS file names
#'  should be removed from \code{frac_smry.xlsx}.
#'
#'  Utilities \code{extract_raws()} and \code{extract_psm_raws()} may aid
#'  matching MS file names between \code{frac_smry.xlsx} and PSM data. Utility
#'  \code{extract_raws()} extracts the list of MS file names under a file
#'  folder. For help, try \code{?extract_raws}. Utility
#'  \code{extract_psm_raws()} extracts the list of MS file names that are
#'  actually present in PSM data. For help, try \code{?extract_psm_raws}.
#'
#'  \code{Reference}: reference entry(entries) are indicated with non-void string(s).
#'
#'  \tabular{ll}{ \strong{Optional default column}   \tab \strong{Descrption}\cr
#'  Select \tab Samples to be selected for indicated analysis \cr Group \tab
#'  Aesthetic labels annotating the prior knowledge of sample groups, e.g.,
#'  Ctrl_T1, Ctrl_T2, Disease_T1, Disease_T2, ...\cr Order \tab Numeric labels
#'  specifying the order of sample \code{groups} \cr Fill \tab Aesthetic labels
#'  for sample annotation by filled color\cr Color \tab Aesthetic labels for
#'  sample annotation by edge color\cr Shape \tab Aesthetic labels for sample
#'  annotation by shape\cr Size \tab Aesthetic labels for sample annotation by
#'  size \cr Alpha \tab Aesthetic labels for sample annotation by transparency
#'  \cr \cr}
#'
#'  \tabular{ll}{ \strong{Exemplary optional open column}   \tab
#'  \strong{Descrption}\cr Term \tab Categorical terms for statistical modeling.
#'  \cr Duplicate \tab Indicators of duplicated samples for corrections in
#'  statistical significance \cr Peptide_Yield \tab Yields of peptides in sample
#'  handling \cr}
#'
#'
#'@section \code{frac_smry.xlsx}: \tabular{ll}{ \strong{Column}   \tab
#'  \strong{Descrption}\cr TMT_Set \tab v.s.  \cr LCMS_Injection   \tab v.s. \cr
#'  Fraction \tab Fraction indexes under a \code{TMT_Set} \cr RAW_File \tab v.s.
#'  }
#'  
#'@family normalization functions
#'@seealso 
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
#'@family data row filtration
#'@seealso 
#'  \emph{Variable arguments of `filter_...`} \cr 
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'@family missing value imputation
#'@seealso 
#'  \emph{Missing values} \cr 
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  
#'@family basic informatics
#'@seealso 
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@param dat_dir A character string to the working directory. The default is to
#'  match the value under the global environment.
#'@param expt_smry A character string to the \code{.xlsx} file containing the
#'  metadata of TMT experiments. The default is \code{expt_smry.xlsx}.
#'@param frac_smry A character string to the \code{.xlsx} file containing
#'  peptide fractionation summary. The default is \code{frac_smry.xlsx}.
#'
#'@example inst/extdata/examples/load_expts_.R
#'
#'@import dplyr rlang fs
#'@importFrom magrittr %>%
#'@export
load_expts <- function (dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx") {
  on.exit(mget(names(formals()), rlang::current_env()) %>% save_call("load_expts"))
  
  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))

  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  
  if (!fs::dir_exists(dat_dir)) {
    new_dat_dir <- fs::path_expand_r(dat_dir)
    new_dat_dir2 <- fs::path_expand(dat_dir)
    
    if (fs::dir_exists(new_dat_dir)) {
      dat_dir <- new_dat_dir
      assign("dat_dir", dat_dir, envir = .GlobalEnv)
      cat("dat_dir <- \"", dat_dir, "\"", sep = "")
    } else if (fs::dir_exists(new_dat_dir2)) {
      dat_dir <- new_dat_dir2
      assign("dat_dir", dat_dir, envir = .GlobalEnv)
      cat("dat_dir <- \"", dat_dir, "\"", sep = "")
    } else {
      stop(dat_dir, " not existed.", call. = FALSE)
    }
    
    rm(new_dat_dir, new_dat_dir2)
  }

  prep_label_scheme(dat_dir, expt_smry)
  prep_fraction_scheme(dat_dir, frac_smry)
}


#' Reload the "expt_smry.xlsx" and "frac_smry.xlsx"
#'
#' @import rlang
#' @importFrom magrittr %>%
#' @importFrom fs file_info
reload_expts <- function() {
  expt_smry <- match_call_arg(load_expts, expt_smry)
  frac_smry <- match_call_arg(load_expts, frac_smry)
  
  fi_xlsx <- fs::file_info(file.path(dat_dir, expt_smry))$change_time
  if (is.na(fi_xlsx)) stop("Time stamp of ", expt_smry, " not available.")
  
  fi_rda <- fs::file_info(file.path(dat_dir, "label_scheme.rda"))$change_time
  if (fi_xlsx > fi_rda) {
    load_expts(dat_dir = dat_dir, expt_smry = !!expt_smry, frac_smry = !!frac_smry)
  }
}


#' Extracts the channel information in TMT experiments
#'
#' A function returns the indexes of TMT channels that are associated to
#' reference(s), sample(s) and probable unused void(s).
#'
#' @param label_scheme The data frame returned by \code{\link{load_expts}}.
#' @param set_idx Numeric.  The index of a multiplex TMT experiment.
#' @return Three lists of indexes: \code{refChannels}, reference channels(s);
#'   \code{emptyChannels}, empty channel(s) that were not used for sample
#'   labeling; \code{labeledChannels}, non-empty channels including both
#'   reference(s) and sample(s).
#'
#' @importFrom dplyr select filter
channelInfo <- function (label_scheme, set_idx) {
	stopifnot(length(set_idx) == 1)

	label_scheme_sub <- label_scheme %>%
	  dplyr::filter(!duplicated(Sample_ID), TMT_Set == set_idx)

	ref <- label_scheme_sub$Reference

	empty_channel_sub <- is.na(label_scheme_sub$Sample_ID) |
	  grepl("^Empty|^Outlier", label_scheme_sub$Sample_ID, ignore.case = TRUE)

	label_scheme_sub <- !empty_channel_sub

	out <- list(
		refChannels = ref,
		emptyChannels = empty_channel_sub,
		labeledChannels = label_scheme_sub
	)

	lapply(out, which)
}


#' Finds the number of multiplex TMT experiments
#'
#' @param label_scheme_full The label_scheme with the probable inclusion of
#'   different LCMS_inj under the same TMT_Set.
n_TMT_sets <- function (label_scheme_full) {
	length(unique(label_scheme_full$TMT_Set))
}


#' Finds the multiplexity of TMT labels
#'
#' \code{TMT_plex} returns the multiplexity of TMT labels.
#' @inheritParams check_label_scheme
TMT_plex <- function (label_scheme_full) {
	nlevels(as.factor(label_scheme_full$TMT_Channel))
}


#' Finds the factor levels of TMT labels
#'
#' \code{TMT_levels} returns the factor levels of TMT labels.
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
TMT_levels <- function (TMT_plex) {
	if (TMT_plex == 16) {
	  TMT_levels <- c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
	                  "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
	                  "TMT-131N", "TMT-131C", "TMT-132N", "TMT-132C", 
	                  "TMT-133N", "TMT-133C", "TMT-134N")
	} else if (TMT_plex == 11) {
	  TMT_levels <- c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
	                  "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
	                  "TMT-131N", "TMT-131C")
	} else if (TMT_plex == 10) {
	  TMT_levels <- c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
	                  "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", "TMT-131")
	} else if (TMT_plex == 6) {
	  TMT_levels <- c("TMT-126", "TMT-127", "TMT-128", "TMT-129", "TMT-130", "TMT-131")
	} else if (TMT_plex == 1) {
	  TMT_levels <- c("TMT-126")
	} else if (TMT_plex == 0) {
	  TMT_levels <- NULL
	}
}


#' Simplifies label schemes from \code{label_scheme_full}
#'
#' Removes duplicated sample entries under different LC/MS injections.
#' @inheritParams load_expts
#' @inheritParams check_label_scheme
simple_label_scheme <- function (dat_dir, label_scheme_full) {
	TMT_plex <- TMT_plex(label_scheme_full)
	TMT_levels <- TMT_levels(TMT_plex)

	label_scheme <- label_scheme_full %>%
		dplyr::filter(!duplicated(Sample_ID), !is.na(Sample_ID)) %>%
		dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
		dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)
	
	if (nrow(label_scheme) <(TMT_plex * n_TMT_sets(label_scheme))) {
	  stop("Duplicated sample ID(s) in `expt_smry.xlsx`", call. = FALSE)
	}

	save(label_scheme, file = file.path(dat_dir, "label_scheme.rda"))
}


#' Checks the uniqueness of sample IDs in \code{label_scheme_full}
#'
#' \code{check_label_scheme} will stop the analysis if the number of unique
#' samples are less than expected.
#' @param label_scheme_full The data frame returned by \code{\link{load_expts}},
#'   including multiple LCMS series.
check_label_scheme <- function (label_scheme_full) {
	load(file = file.path(dat_dir, "label_scheme.rda"))

	TMT_plex <- TMT_plex(label_scheme)
	if(!is.null(TMT_plex)) {
		if((nlevels(as.factor(label_scheme$Sample_ID))) < 
		   (TMT_plex * nlevels(as.factor(label_scheme$TMT_Set))))
			stop("Not enough observations in unique 'Sample_ID'")
	}
}


#' Find mismatches in RAW file names
#'
#' \code{check_raws} finds mismatched RAW files between expt_smry.xlsx and
#' PSM outputs.
#' @param df A data frame containing the PSM table from database searches.
check_raws <- function(df) {
  stopifnot ("RAW_File" %in% names(df))
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))

  ## program-generated frac_smry.xlsx may be based on wrong information from expt_smry.xlsx
  ls_raws <- label_scheme_full$RAW_File %>% unique()
  fs_raws <- fraction_scheme$RAW_File %>% unique()
  if (!(all(is.na(ls_raws)) | all(ls_raws %in% fs_raws))) {
    load(file.path(dat_dir, "Calls", "load_expts.rda"))
    fn_frac <- call_pars$frac_smry
    unlink(file.path(dat_dir, fn_frac))
    prep_fraction_scheme(dat_dir, fn_frac)
    load(file = file.path(dat_dir, "fraction_scheme.rda"))
  }
  
  tmtinj_raw <- fraction_scheme %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>%
    dplyr::select(-Fraction) %>%
    dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
		dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File)) # Bruker
  
  ms_raws <- df$RAW_File %>% unique()
  label_scheme_raws <- tmtinj_raw$RAW_File %>% unique()
  
  missing_ms_raws <- ms_raws %>% .[! . %in% label_scheme_raws]
  wrong_label_scheme_raws <- label_scheme_raws[! label_scheme_raws %in% ms_raws]
  
  if(!purrr::is_empty(missing_ms_raws) | !purrr::is_empty(wrong_label_scheme_raws)) {
    cat("Required RAW MS file name(s) not found from the `expt_smry.xlsx` and/or `frac_smry.xlsx`:\n")
    cat(paste0(missing_ms_raws, "\n"))
    
    cat("RAW MS files in `expt_smry.xlsx` and/or `frac_smry.xlsx` but not present in PSM data:\n")
    cat(paste0("\t", wrong_label_scheme_raws, "\n"))
    
    stop("Check file names under the `RAW_File` column in `expt_smry.xlsx` and/or `frac_smry.xlsx`.", 
         call. = FALSE)
  }

  return(tmtinj_raw)
}

#' Make new column names
#'
#' \code{newColnames} match names to Sample_ID in label_scheme
#'
#' @param i Integer; the index of TMT experiment 
#' @param x Data frame; log2FC data
#' @param label_scheme Experiment summary
#' @import dplyr purrr rlang forcats
#' @importFrom magrittr %>%
newColnames <- function(i, x, label_scheme) {
  label_scheme_sub <- label_scheme %>%
    dplyr::filter(TMT_Set == i)
  
  cols <- grep(paste0("[RI][0-9]{3}[NC]*_", i, "$"), names(x))
  nm_channel <- gsub(paste0("([RI][0-9]{3}[NC]*)_", i, "$"), "\\1", names(x)[cols])
  names(x)[cols] <- paste0(nm_channel, " (", as.character(label_scheme_sub$Sample_ID), ")")
  
  cols <- grep("[RI][0-9]{3}.*\\s+\\(.*\\)$", names(x))
  
  # cols with new names go first
  if (length(cols) < ncol(x)) x <- dplyr::bind_cols(x[, cols], x[, -cols, drop = FALSE])
  
  return(x)
}


#' combined peptide reports across multiple TMT experiments
#' 
#' median summarization of data from the same TMT experiment at different LCMS injections
#' summed \code{pep_n_psm}, \code{prot_n_psm}, and \code{prot_n_pep} after data merging
#' no Z_log2_R yet available
#'   use \code{col_select = expr(Sample_ID)} not \code{col_select} to get all Z_log2_R
#'   why: users may specify \code{col_select} only partial to Sample_ID entries
#'
#' @inheritParams info_anal
#' @inheritParams normPSM
normPep_Mplex <- function (id = "pep_seq_mod", group_pep_by = "prot_acc", ...) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  id <- rlang::as_string(rlang::enexpr(id))
  
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  
  filelist <- list.files(path = file.path(dat_dir, "Peptide"), 
                         pattern = paste0("TMTset[0-9]+_LCMSinj[0-9]+_Peptide_N\\.txt$"), full.names = TRUE)
  
  if (is_empty(filelist)) stop("No individual peptide tables; run `PSM2Pep()` first.")

  df <- do.call(rbind, 
                purrr::map(filelist, read.csv, check.names = FALSE, header = TRUE, 
                           sep = "\t", comment.char = "#")) %>%
    dplyr::mutate(TMT_Set = factor(TMT_Set)) %>%
    dplyr::arrange(TMT_Set) 
  
  df <- df %>% filters_in_call(!!!filter_dots)
  
  if ("gene" %in% names(df)) 
    df <- df %>% dplyr::mutate(gene = forcats::fct_explicit_na(gene))

  df <- local({
    dup_peps <- df %>%
      dplyr::select(!!rlang::sym(id), prot_acc) %>%
      dplyr::group_by(!!rlang::sym(id)) %>%
      dplyr::summarise(N = n_distinct(prot_acc)) %>%
      dplyr::filter(N > 1)
    
    if (nrow(dup_peps) > 0) {
      write.csv(dup_peps, file.path(dat_dir, "Peptide", "dbl_dipping_peptides.csv"), row.names = FALSE)
      df <- df %>% dplyr::filter(! (!!rlang::sym(id) %in% dup_peps[[id]]))
    }
    
    write.csv(df, file.path(dat_dir, "Peptide\\cache", "unambi_peptides.csv"), row.names = FALSE)
    return(df)
  })
  
  df_num <- df %>% 
    dplyr::select(!!rlang::sym(id), 
                  TMT_Set, 
                  grep("^sd_log2_R[0-9]{3}", names(.)), 
                  grep("^log2_R[0-9]{3}", names(.)), 
                  grep("^N_log2_R[0-9]{3}", names(.)), 
                  grep("^Z_log2_R[0-9]{3}", names(.)), 
                  grep("^I[0-9]{3}", names(.)), 
                  grep("^N_I[0-9]{3}", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id), TMT_Set) %>%
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_num <- df_num %>%
    dplyr::arrange(TMT_Set) %>%
    tidyr::gather(grep("R[0-9]{3}|I[0-9]{3}", names(.)), key = ID, value = value) %>%
    tidyr::unite(ID, ID, TMT_Set)
  
  # define the levels of TMT channels;
  # otherwise, the order of channels will flip between N(itrogen) and C(arbon)
  Levels <- unique(df_num$ID)
  df_num <- df_num %>%
    dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
    tidyr::spread(ID, value)
  rm(Levels)
  
  for (set_idx in seq_len(n_TMT_sets(label_scheme))) {
    df_num <- newColnames(set_idx, df_num, label_scheme)
  }
  
  df_num <- df_num %>% 
    dplyr::select(!!rlang::sym(id), grep("[RI][0-9]{3}[NC]*", names(.))) %>% 
    dplyr::arrange(!!rlang::sym(id)) %T>%
    write.csv(file.path(dat_dir, "Peptide\\cache", "pep_num.csv"), row.names = FALSE)
  
  pep_n_psm <- df %>%
    dplyr::select(!!rlang::sym(id), pep_n_psm) %>%
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarise(pep_n_psm = sum(pep_n_psm)) %>% 
    dplyr::arrange(!!rlang::sym(id))
  
  prot_n_psm <- df %>% 
    dplyr::select(pep_n_psm, !!rlang::sym(group_pep_by)) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_n_psm = sum(pep_n_psm))

  prot_n_pep <- df %>% 
    dplyr::select(!!rlang::sym(id), !!rlang::sym(group_pep_by)) %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_n_pep = n())
  
  df_first <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    med_summarise_keys(id) %>% 
    dplyr::select(-pep_n_psm, -prot_n_psm, -prot_n_pep, -TMT_Set) %>% # remove old values from single `TMT_Set`
    dplyr::arrange(!!rlang::sym(id))  

  df <- list(pep_n_psm, df_first, df_num) %>%
    purrr::reduce(left_join, by = id)
  
  df <- list(df, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(left_join, by = group_pep_by)
  
  if (("pep_seq_mod" %in% names(df)) & (match_call_arg(normPSM, use_lowercase_aa))) {
    df <- df %>% 
      dplyr::mutate(pep_mod_protnt = ifelse(grepl("^[A-z\\-]\\.~", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_protntac = ifelse(grepl("^[A-z\\-]\\._", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_pepnt = ifelse(grepl("^[A-z\\-]\\.[_~]?\\^", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_m = ifelse(grepl("m", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_n = ifelse(grepl("n", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_sty = ifelse(grepl("[sty]", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_pepct = ifelse(grepl("[\\^]{1}[_~]?\\.[A-z\\-]{1}$", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_protctam = ifelse(grepl("_{1}\\.[A-z\\-]{1}$", pep_seq_mod), TRUE, FALSE)) %>% 
      dplyr::mutate(pep_mod_protct = ifelse(grepl("~{1}\\.[A-z\\-]{1}$", pep_seq_mod), TRUE, FALSE))
  }
  
  df <- dplyr::bind_cols(
    df %>% select(grep("^pep_", names(.))), 
    df %>% select(-grep("^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% select(grep("^prot_", names(.))), 
    df %>% select(-grep("^prot_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^N_I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^log2_R[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("^Z_log2_R[0-9]{3}[NC]*", names(.))), 
  )
  
  df <- df %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
    dplyr::mutate_at(vars(grep("sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3))
  
  df <- df %>%
    dplyr::filter(!duplicated(.[[id]])) %>% 
    dplyr::filter(rowSums(!is.na(.[, grepl("N_log2_R", names(.))])) > 0) %>% 
    dplyr::arrange(!!rlang::sym(id)) %T>%
    write.csv(file.path(dat_dir, "Peptide\\cache", "Peptide_no_norm.csv"), row.names = FALSE)	
  
  df <- normMulGau(
    df = df,
    method_align = "MC",
    n_comp = 1L,
    range_log2r = c(0, 100),
    range_int = c(0, 100),
    filepath = file.path(dat_dir, "Peptide\\Histogram"),
    col_select = rlang::expr(Sample_ID), 
  )
}


#' Summary of peptide keys by mean or geomean
#' 
#' @inheritParams info_anal
#' @import dplyr purrr rlang  magrittr
med_summarise_keys <- function(df, id) {
  mascot_median_keys <- c("pep_score", "pep_rank", "pep_isbold", "pep_exp_mr", "pep_delta", 
                          "pep_exp_mz", "pep_exp_z")
  mascot_geomean_keys <- c("pep_expect")
  
  df_mascot_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mascot_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mascot_geomean <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mascot_geomean_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ my_geomean(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(mascot_median_keys, mascot_geomean_keys)))
  
  df_mq_rptr_mass_dev <- df %>% 
    dplyr::select(!!rlang::sym(id), grep("^Reporter mass deviation", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep("^Reporter mass deviation", names(.)))	  
  
  mq_median_keys <- c(
    "Score", 
    "Charge", "Mass", "PIF", "Fraction of total spectrum", "Mass error [ppm]", 
    "Mass error [Da]", "Base peak fraction", "Precursor Intensity", 
    "Precursor Apex Fraction", "Intensity coverage", "Peak coverage", 
    "Combinatorics"
  )
  mq_geomean_keys <- c("PEP")
  
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df_mq_geomean <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_geomean_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ my_geomean(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(mq_median_keys, mq_geomean_keys)))
  
  sm_median_keys <- c(
    "score", "parent_charge", 
    "deltaForwardReverseScore", "percent_scored_peak_intensity", "totalIntensity", 
    "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
    "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
    "delta_parent_mass", "delta_parent_mass_ppm")
  sm_geomean_keys <- NA
  
  df_sm_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% sm_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% sm_median_keys))
  
  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id)))
  
  df <- list(df_first, 
             df_mascot_med, df_mascot_geomean, 
             df_mq_rptr_mass_dev, df_mq_med, df_mq_geomean, 
             df_sm_med) %>%
    purrr::reduce(left_join, by = id) %>%
    data.frame(check.names = FALSE)
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
}


#' load prior Peptide.txt
#' @inheritParams info_anal
load_prior <- function(filename, id) {
  stopifnot(file.exists(filename))
  
  df <- read.csv(filename, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0) 
  
  if (! id %in% names(df)) {
    try(unlink(file.path(dat_dir, "Peptide\\Peptide.txt")))
    try(unlink(file.path(dat_dir, "Protein\\Protein.txt")))
    stop("`Peptide.txt` deleted as column `", id, "` not available.", call. = FALSE)
  }
  
  df <- df %>% dplyr::arrange(!!rlang::sym(id))
}


#' format numeric columns
#' @inheritParams info_anal
fmt_num_cols <- function (df) {
  df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-  
    df[, grepl("^Z_log2_R[0-9]{3}", names(df))] %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    round(., digits = 3)
  
  df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-  
    df[, grepl("^N_log2_R[0-9]{3}", names(df))] %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    round(., digits = 3)
  
  df[, grepl("I[0-9]{3}", names(df))] <-  
    df[, grepl("I[0-9]{3}", names(df))] %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    round(., digits = 0)
  
  return(df)
}


#'Merge peptide table(s) into one
#'
#'\code{mergePep} merges individual peptide table(s),
#'\code{TMTset1_LCMSinj1_Peptide_N.txt, TMTset1_LCMSinj2_Peptide_N.txt} etc.,
#'into one interim \code{Peptide.txt}. The \code{log2FC} values in the interim
#'result are centered with the medians at zero (median centering). The utility
#'is typically applied after the conversion of PSMs to peptides via
#'\code{\link{PSM2Pep}} and is required even for a experiment with one multiplex
#'TMT and one LCMS injection.
#'
#'In the interim output file, "\code{Peptide.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'median-centered \code{log2_R...} without scaling normalization. Values under
#'columns \code{Z_log2_R...} are \code{N_log2_R...} with additional scaling
#'normalization. Values under columns \code{I...} are \code{reporter-ion
#'intensity} before normalization. Values under columns \code{N_I...} are
#'normalized \code{I...}. Values under columns \code{sd_log2_R...} are the
#'standard deviation of the \code{log2FC} of proteins from ascribing peptides.
#'
#'Description of the column keys in the output: \cr \code{system.file("extdata",
#'"mascot_peptide_keys.txt", package = "proteoQ")} \cr
#'\code{system.file("extdata", "maxquant_peptide_keys.txt", package =
#'"proteoQ")}
#'
#'The peptide counts in individual peptide tables,
#'\code{TMTset1_LCMSinj1_Peptide_N.txt} etc., may be fewer than the entries
#'indicated under the \code{prot_n_pep} column after the peptide
#'removals/cleanups using \code{purgePSM}.
#'
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  of data against the column keys in individual peptide tables of
#'  \code{TMTset1_LCMSinj1_Peptide_N.txt, TMTset1_LCMSinj2_Peptide_N.txt}, etc.
#'  \cr \cr The variable argument statements should be in the following format:
#'  each statement contains to a list of logical expression(s). The \code{lhs}
#'  needs to start with \code{filter_}. The logical condition(s) at the
#'  \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis. For
#'  example, \code{pep_len} is a column key present in \code{Mascot} peptide
#'  tables of \code{TMTset1_LCMSinj1_Peptide_N.txt},
#'  \code{TMTset1_LCMSinj2_Peptide_N.txt} etc. The statement
#'  \code{filter_peps_at = exprs(pep_len <= 50)} will remove peptide entries
#'  with \code{pep_len > 50}. See also \code{\link{normPSM}}.
#'@inheritParams normPSM
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@return The primary output is in \code{...\\Peptide\\Peptide.txt}.
#'
#'@example inst/extdata/examples/mergePep_.R
#'@import stringr dplyr tidyr purrr data.table rlang
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@importFrom plyr ddply
#'@export
mergePep <- function (plot_log2FC_cv = TRUE, ...) {
  dir.create(file.path(dat_dir, "Peptide\\cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide\\Histogram"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide\\log2FC_cv\\raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide\\log2FC_cv\\purged"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Peptide\\log"), recursive = TRUE, showWarnings = FALSE)
  
  old_opt <- options(max.print = 99999, warn = 0)
  options(max.print = 2000000, warn = 1)
  on.exit(options(old_opt), add = TRUE)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("mergePep"), add = TRUE)

  reload_expts()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  id <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  filename <- file.path(dat_dir, "Peptide\\Peptide.txt")
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]

  message("Primary column keys in `Peptide/TMTset1_LCMSinj1_Peptide_N.txt` etc. for `filter_` varargs.")
  
  df <- normPep_Mplex(!!id, group_pep_by, !!!filter_dots) %T>% 
    write.table(filename, sep = "\t", col.names = TRUE, row.names = FALSE)

  if (plot_log2FC_cv & TMT_plex(label_scheme) > 0) {
    quiet_out <- purrr::quietly(sd_violin)(df = df, id = !!group_pep_by, 
                                           filepath = file.path(dat_dir, "Peptide\\log2FC_cv\\raw", "Peptide_sd.png"), 
                                           width = 8 * n_TMT_sets(label_scheme), height = 8, 
                                           type = "log2_R", adjSD = FALSE, is_psm = FALSE)
  }
}


#'Standardize peptide results
#'
#'\code{standPep} standardizes peptide results from \code{\link{mergePep}} with
#'additional choices in data alignment. The utility is typically applied after
#'the assembly of peptide data via \code{\link{mergePep}}. It further supports
#'iterative normalization against data under selected sample columns, data rows
#'or both.
#'
#'
#'In the primary output file, "\code{Peptide.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'aligned \code{log2_R...} according to \code{method_align} without scaling
#'normalization. Values under columns \code{Z_log2_R...} are \code{N_log2_R...}
#'with additional scaling normalization. Values under columns \code{I...} are
#'\code{reporter-ion intensity} before normalization. Values under columns
#'\code{N_I...} are normalized \code{I...}. Values under columns
#'\code{sd_log2_R...} are the standard deviation of the \code{log2FC} of
#'proteins from ascribing peptides.
#'
#'In general, median statistics is applied when summarizing numeric peptide data
#'from different LCMS series. One exception is \code{pep_expect} with Mascot
#'workflow where geometric mean is used.
#'
#'Description of the column keys in the inputs and outputs: \cr
#'\code{system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ")}
#'\cr \code{system.file("extdata", "maxquant_peptide_keys.txt", package =
#'"proteoQ")}
#'
#'@param method_align Character string indicating the method in aligning
#'  \code{log2FC} across samples. \code{MC}: median-centering; \code{MGKernel}:
#'  the kernel density defined by multiple Gaussian functions
#'  (\code{\link[mixtools]{normalmixEM}}). At the \code{MC} default, the ratio
#'  profiles of each sample will be aligned in that the medians of the
#'  \code{log2FC} are zero. At \code{MGKernel}, the ratio profiles of each
#'  sample will be aligned in that the \code{log2FC} at the maximums of kernel
#'  density are zero.
#'@param col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  At the \code{NULL} default, the column key of \code{Select} in
#'  \code{expt_smry.xlsx} will be used. In the case of no samples being
#'  specified under \code{Select}, the column key of \code{Sample_ID} will be
#'  used. The non-empty entries under the ascribing column will be used in
#'  indicated analysis.
#'@param range_log2r Numeric vector at length two. The argument specifies the
#'  range of the \code{log2FC} for use in the scaling normalization of standard
#'  deviation across samples. The default is between the 10th and the 90th
#'  quantiles.
#'@param range_int Numeric vector at length two. The argument specifies the
#'  range of the \code{intensity} of reporter ions for use in the scaling
#'  normalization of standard deviation across samples. The default is between
#'  the 5th and the 95th quantiles.
#'@param n_comp Integer; the number of Gaussian components to be used with
#'  \code{method_align = MGKernel}. A typical value is 2 or 3. The variable
#'  \code{n_comp} overwrites the argument \code{k} in
#'  \code{\link[mixtools]{normalmixEM}}.
#'@param seed Integer; a seed setting a starting point for reproducible
#'  analyses.
#'@param ... \code{slice_}: variable argument statements for the identification
#'  of row subsets. The partial data will be taken for parameterizing the
#'  alignment of \code{log2FC} across samples. The full data set will be updated
#'  subsequently with the newly derived paramters. Note that there is no data
#'  entry removals from the complete data set with the \code{slice_} procedure.
#'  \cr \cr The variable argument statements should be in the following format:
#'  each of the statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{slice_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key present in \code{Peptide.txt}
#'  with \code{Mascot} workflows. The \code{slice_peps_at = exprs(pep_len >= 10,
#'  pep_len <= 50)} will extract peptide entries with the number of amino acid
#'  residues betwen 10 and 50 for \code{log2FC} alignment. Shorter or longer
#'  peptide sequences will remain in \code{Peptide.txt} but not used in the
#'  parameterization. See also \code{\link{normPSM}} for the variable arguments
#'  of \code{filter_}. \cr \cr Additional parameters from
#'  \code{\link[mixtools]{normalmixEM}}, i.e., \cr \code{maxit}, the maximum
#'  number of iterations allowed; \cr \code{epsilon}, tolerance limit for
#'  declaring algorithm convergence.
#'@inheritParams normPSM
#'@inheritParams mixtools::normalmixEM
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output is in \code{...\\Peptide\\Peptide.txt}.
#'
#'@example inst/extdata/examples/normPep_.R
#'
#'@import stringr dplyr tidyr purrr data.table rlang
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@importFrom plyr ddply
#'@export
standPep <- function (method_align = c("MC", "MGKernel"), col_select = NULL, range_log2r = c(10, 90), 
                    range_int = c(5, 95), n_comp = NULL, seed = NULL, plot_log2FC_cv = FALSE, ...) {

  old_opt <- options(max.print = 99999)
  options(max.print = 2000000)
  on.exit(options(old_opt), add = TRUE)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("standPep"), add = TRUE)
  
  reload_expts()
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))

  filename <- file.path(dat_dir, "Peptide\\Peptide.txt")
  if (!file.exists(filename)) stop(filename, " not found; run `mergePep(...)` first", call. = FALSE)

  id <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  method_align <- rlang::enexpr(method_align)
  if (method_align == rlang::expr(c("MC", "MGKernel"))) {
    method_align <- "MC"
  } else {
    method_align <- rlang::as_string(method_align)
  }
  
  col_select <- rlang::enexpr(col_select)
  col_select <- ifelse(is.null(col_select), rlang::expr(Sample_ID), rlang::sym(col_select))
  
  if (is.null(label_scheme[[col_select]])) {
    col_select <- rlang::expr(Sample_ID)
    warning("Column \'", rlang::as_string(col_select), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    col_select <- rlang::expr(Sample_ID)
    warning("No samples were specified under column \'", rlang::as_string(col_select), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  }
  
  stopifnot(length(range_log2r) == 2)
  stopifnot(length(range_int) == 2)
  
  if (range_log2r[2] <= 1) range_log2r <- range_log2r * 100
  if (range_int[2] <= 1) range_int <- range_int * 100
  
  stopifnot(range_log2r[1] < range_log2r[2] & range_log2r[1] >= 0 & range_log2r[2] <= 100)
  stopifnot(range_int[1] < range_int[2] & range_int[1] >= 0 & range_int[2] <= 100)
  
  dots <- rlang::enexprs(...)
  
  message("Primary column keys in `Peptide/Peptide.txt` for `slice_` varargs.")

  df <- load_prior(filename, id) %>% 
    normMulGau(
      df = .,
      method_align = method_align,
      n_comp = n_comp,
      seed = seed,
      range_log2r = range_log2r,
      range_int = range_int,
      filepath = file.path(dat_dir, "Peptide\\Histogram"),
      col_select = col_select, 
      !!!dots, 
    ) %>% 
    fmt_num_cols() %T>% 
    write.table(file.path(dat_dir, "Peptide", "Peptide.txt"), 
                sep = "\t", col.names = TRUE, row.names = FALSE)

  if (plot_log2FC_cv & TMT_plex(label_scheme) > 0) {
    sd_violin(df = df, id = !!group_pep_by, 
              filepath = file.path(dat_dir, "Peptide\\log2FC_cv\\raw", "Peptide_sd.png"), 
              width = 8 * n_TMT_sets(label_scheme), height = 8, 
              type = "log2_R", adjSD = FALSE, is_psm = FALSE)
  }
}



#'Interim protein data
#'
#'\code{Pep2Prn} summarizes \code{Peptide.txt} to an interim protein report in
#'\code{Protein.txt}.
#'
#'Fields other than \code{log2FC} and \code{intensity} are summarized with
#'median statistics.
#'
#'@param method_pep_prn Character string; the method to summarize the
#'  \code{log2FC} and the \code{intensity} of peptides by protein entries. The
#'  descriptive statistics includes \code{c("mean", "median", "top.3",
#'  "weighted.mean")} with \code{median} being the default. The representative
#'  \code{log10-intensity} of reporter ions at the peptide levels (from
#'  \code{\link{standPep}}) will be the weight when summarizing \code{log2FC}
#'  with \code{top.3} or \code{weighted.mean}.
#'@param use_unique_pep Logical. If TRUE, only entries that are \code{TRUE} or
#'  equal to \code{1} under the column \code{pep_isunique} in \code{Peptide.txt}
#'  will be used, for summarizing the \code{log2FC} and the \code{intensity} of
#'  peptides into protein values. The default is to use unique peptides only.
#'  For \code{MaxQuant} data, the levels of uniqueness are according to the
#'  \code{pep_unique_by} in \code{\link{normPSM}}. The argument currently do
#'  nothing to \code{Spectrum Mill} data where both unique and shared peptides
#'  will be kept.
#'@param ... \code{filter_}: Variable argument statements for the filtration of
#'  data rows. Each statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_len} is a column key present in \code{Peptide.txt}
#'  with \code{Mascot} workflows. The statement of \code{filter_peps_at =
#'  exprs(pep_len <= 50)} will remove peptide entries with \code{pep_len > 50}.
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@return The primary output in "\code{...\\Protein\\Protein.txt}".
#'
#'@example inst/extdata/examples/Pep2Prn_.R
#'@import stringr dplyr purrr rlang  magrittr
#'@export
Pep2Prn <- function (method_pep_prn = c("median", "mean", "weighted.mean", "top.3"), 
                     use_unique_pep = TRUE, ...) {
  
  dir.create(file.path(dat_dir, "Protein\\Histogram"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein\\cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein\\log"), recursive = TRUE, showWarnings = FALSE)
  
  old_opt <- options(max.print = 99999)
  options(max.print = 2000000)
  on.exit(options(old_opt), add = TRUE)
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("Pep2Prn"), add = TRUE)
  
  reload_expts()
  
  method_pep_prn <- rlang::enexpr(method_pep_prn)
  if (method_pep_prn == rlang::expr(c("median", "mean", "weighted.mean", "top.3"))) {
    method_pep_prn <- "median"
  } else {
    method_pep_prn <- rlang::as_string(method_pep_prn)
    stopifnot(method_pep_prn %in% c("mean", "top.3", "median", "weighted.mean"))
  }
  
  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(id %in% c("prot_acc", "gene"))
  
  if (id == "gene") {
    gn_rollup <- TRUE
    id <- "prot_acc"
  } else {
    gn_rollup <- FALSE
  }
  
  stopifnot(id == "prot_acc") 
  
  message("Primary column keys in `Peptide/Peptide.txt` for `filter_` varargs.")
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  
  df <- pep_to_prn(!!id, method_pep_prn, use_unique_pep, gn_rollup, !!!filter_dots) 
  
  to_mc <- TRUE
  if (to_mc) {
    df <- normMulGau(
      df = df,
      method_align = "MC",
      n_comp = 1L,
      range_log2r = c(0, 100),
      range_int = c(0, 100),
      filepath = file.path(dat_dir, "Protein\\Histogram"),
      col_select = rlang::expr(Sample_ID), 
    )    
  }

  df <- df %>% 
    dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %T>% 
    write.table(., file.path(dat_dir, "Protein\\Protein.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
}


#' Helper of Pep2Prn
#' 
#' @param gn_rollup Logical; if TRUE, rolls up protein accessions to gene names.
#' @inheritParams info_anal
#' @inheritParams Pep2Prn
pep_to_prn <- function(id, method_pep_prn, use_unique_pep, gn_rollup, ...) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  id <- rlang::as_string(rlang::enexpr(id))
  
  filter_dots <- rlang::enexprs(...) %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  fn_fasta <- file.path(dat_dir, "my_project.fasta")
  stopifnot(file.exists(fn_fasta))
  fasta <- seqinr::read.fasta(fn_fasta, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
  
  df <- read.csv(file.path(dat_dir, "Peptide\\Peptide.txt"), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
  
  df <- df %>% filters_in_call(!!!filter_dots)
  
  if (use_unique_pep & "pep_isunique" %in% names(df)) df <- df %>% dplyr::filter(pep_isunique == 1)
  
  df_num <- df %>% 
    dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id))
  
  df_num <- switch(method_pep_prn, 
                   mean = aggrNums(mean)(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   top.3 = TMT_top_n(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   weighted.mean = TMT_wt_mean(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   median = aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE), 
                   aggrNums(median)(df_num, !!rlang::sym(id), na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.)))
  
  df_mq_rptr_mass_dev <- df %>% 
    dplyr::select(!!rlang::sym(id), grep("^Reporter mass deviation", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-grep("^Reporter mass deviation", names(.)))	  
  
  mq_median_keys <- c(
    "Score", "Missed cleavages", "PEP", 
    "Charge", "Mass", "PIF", "Fraction of total spectrum", "Mass error [ppm]", 
    "Mass error [Da]", "Base peak fraction", "Precursor Intensity", 
    "Precursor Apex Fraction", "Intensity coverage", "Peak coverage", 
    "Combinatorics"
  )
  
  df_mq_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% mq_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% mq_median_keys))		
  
  sm_median_keys <- c(
    "deltaForwardReverseScore", "percent_scored_peak_intensity", "totalIntensity", 
    "precursorAveragineChiSquared", "precursorIsolationPurityPercent", 
    "precursorIsolationIntensity", "ratioReporterIonToPrecursor", 
    "matched_parent_mass", "delta_parent_mass", "delta_parent_mass_ppm")
  
  df_sm_med <- df %>% 
    dplyr::select(!!rlang::sym(id), which(names(.) %in% sm_median_keys)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% sm_median_keys))
  
  df_first <- df %>% 
    dplyr::filter(!duplicated(!!rlang::sym(id))) %>% 
    dplyr::select(-grep("^pep_", names(.)))    
  
  df <- list(df_first, 
             df_mq_rptr_mass_dev, df_mq_med, 
             df_sm_med, 
             df_num) %>% 
    purrr::reduce(left_join, by = id) %>% 
    data.frame(check.names = FALSE)
  
  rm(df_num, df_first)
  
  df[, grepl("log2_R[0-9]{3}", names(df)) & !sapply(df, is.logical)] <- 
    df[, grepl("log2_R[0-9]{3}", names(df)) & !sapply(df, is.logical)] %>% 
    dplyr::mutate_if(is.integer, as.numeric) %>% 
    round(., digits = 3)
  
  df[, grepl("I[0-9]{3}", names(df)) & !sapply(df, is.logical)] <- 
    df[, grepl("I[0-9]{3}", names(df)) & !sapply(df, is.logical)] %>% 
    dplyr::mutate_if(is.integer, as.numeric) %>% 
    round(., digits = 0)
  
  df <- cbind.data.frame(
    df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df))], 
    df[, grep("^I[0-9]{3}", names(df))], 
    df[, grep("^N_I[0-9]{3}", names(df))], 
    df[, grep("^log2_R[0-9]{3}", names(df))], 
    df[, grep("^N_log2_R[0-9]{3}", names(df))], 
    df[, grep("^Z_log2_R[0-9]{3}", names(df))])
  
  df <- df %>% 
    .[rowSums(!is.na(.[, grepl("N_log2_R", names(.))])) > 0, ]
  
  if (gn_rollup) {
    dfa <- df %>% 
      dplyr::select(gene, grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
      dplyr::filter(!is.na(gene)) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::summarise_all(list(~ median(.x, na.rm = TRUE)))
    
    dfb <- df %>% 
      dplyr::select(-prot_cover, -grep("I[0-9]{3}|log2_R[0-9]{3}", names(.))) %>% 
      dplyr::filter(!is.na(gene)) %>% 
      dplyr::filter(!duplicated(.$gene))
    
    dfc <- df %>% 
      dplyr::select(gene, prot_cover) %>% 
      dplyr::filter(!is.na(gene), !is.na(prot_cover)) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::mutate(prot_cover = as.numeric(sub("%", "", prot_cover))) %>% 
      dplyr::summarise_all(~ max(.x, na.rm = TRUE)) %>% 
      dplyr::mutate(prot_cover = paste0(prot_cover, "%"))
    
    df <- list(dfc, dfb, dfa) %>% 
      purrr::reduce(right_join, by = "gene") %>% 
      dplyr::filter(!is.na(gene), !duplicated(gene))
  }
  
  return(df)
}



#'Standardize protein results
#'
#'\code{standPrn} standardizes protein results from \code{\link{Pep2Prn}} with
#'additional choices in data alignment. The utility further supports iterative
#'normalization against data under selected sample columns, data rows or both.
#'
#'In the primary output file, "\code{Protein.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'aligned \code{log2_R...} according to \code{method_align} without scaling
#'normalization. Values under columns \code{Z_log2_R...} are \code{N_log2_R...}
#'with additional scaling normalization. Values under columns \code{I...} are
#'\code{reporter-ion intensity} before normalization. Values under columns
#'\code{N_I...} are normalized \code{I...}. Values under columns
#'\code{sd_log2_R...} are the standard deviation of the \code{log2FC} of
#'proteins from ascribing peptides.
#'
#'@param cache Not currently used.
#'@param ... \code{slice_}: Variable argument statements for the identification
#'  of row subsets. The partial data will be taken for parameterizing the
#'  alignment of \code{log2FC} across samples. The full data set will be updated
#'  subsequently with the newly derived parameters. Note that there is no data
#'  entry removals from the complete data set with the \code{slice_} procedure.
#'  \cr \cr The variable argument statements should be in the following format:
#'  each of the statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{slice_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{prot_n_pep} is a column key present in
#'  \code{Protein.txt}. The \code{slice_prns_at = exprs(prot_n_pep >= 5)} will
#'  extract protein entries with five or more identifying peptide sequences for
#'  \code{log2FC} alignment. Protein entries with less than five identifying
#'  sequences will remain in \code{Protein.txt} but not used in the
#'  parameterization. See also \code{\link{normPSM}} for the variable arguments
#'  of \code{filter_}. \cr \cr Additional parameters from
#'  \code{\link[mixtools]{normalmixEM}}, i.e., \cr \code{maxit}, the maximum
#'  number of iterations allowed; \cr \code{epsilon}, tolerance limit for
#'  declaring algorithm convergence.
#'@inheritParams normPSM
#'@inheritParams standPep
#'@inheritParams mixtools::normalmixEM
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output is in \code{...\\Protein\\Protein.txt}.
#'
#'@example inst/extdata/examples/normPrn_.R
#'@import stringr dplyr tidyr purrr data.table rlang
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@importFrom plyr ddply
#'@export
standPrn <- function (method_align = c("MC", "MGKernel"), 
                    range_log2r = c(10, 90), range_int = c(5, 95), n_comp = NULL, seed = NULL, 
                    col_select = NULL, cache = TRUE, ...) {
  
  dir.create(file.path(dat_dir, "Protein\\Histogram"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein\\cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein\\log"), recursive = TRUE, showWarnings = FALSE)
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 2000000)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("standPrn"), add = TRUE)
  
  reload_expts()
  
  method_align <- rlang::enexpr(method_align)
  if (method_align == rlang::expr(c("MC", "MGKernel"))) {
    method_align <- "MC"
  } else {
    method_align <- rlang::as_string(method_align)
    if(! method_align %in% c("MC", "MGKernel")) 
      warning("Assume the value of 'method_align' is a list of housekeeping proteins")
  }
  
  stopifnot(length(range_log2r) == 2)
  stopifnot(length(range_int) == 2)
  
  if (range_log2r[2] <= 1) range_log2r <- range_log2r * 100
  if (range_int[2] <= 1) range_int <- range_int * 100
  
  stopifnot(range_log2r[1] < range_log2r[2] & 
              range_log2r[1] >= 0 & range_log2r[2] <= 100)
  
  stopifnot(range_int[1] < range_int[2] & 
              range_int[1] >= 0 & range_int[2] <= 100)
  
  id <- match_call_arg(normPSM, group_pep_by)
  pep_id <- match_call_arg(normPSM, group_psm_by)
  
  col_select <- rlang::enexpr(col_select)
  col_select <- ifelse(is.null(col_select), rlang::expr(Sample_ID), rlang::sym(col_select))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  if (is.null(label_scheme[[col_select]])) {
    col_select <- rlang::expr(Sample_ID)
    warning("Column \'", rlang::as_string(col_select), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    col_select <- rlang::expr(Sample_ID)
    warning("No samples were specified under column \'", rlang::as_string(col_select), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  }
  
  dots <- rlang::enexprs(...)
  
  filename <- file.path(dat_dir, "Protein\\Protein.txt")
  
  if (!file.exists(filename)) {
    stop(filename, " not found; run `Pep2Prn` first.")
  }
  
  df <- read.csv(filename, sep = "\t", check.names = FALSE, header = TRUE, comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)

  message("Primary column keys in `Protein/Protein.txt` for `slice_` varargs.")
  
  df <- normMulGau(
    df = df, 
    method_align = method_align, 
    n_comp = n_comp, 
    seed = seed, 
    range_log2r = range_log2r, 
    range_int = range_int, 
    filepath = file.path(dat_dir, "Protein\\Histogram"), 
    col_select = col_select, 
    !!!dots,
  )
  
  df <- df %>% 
    dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %T>% 
    write.table(., file.path(dat_dir, "Protein", "Protein.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
}

#' Extract RAW MS file names
#'
#' Extract a list of \code{RAW} file names that can be passed to \code{frac_smry.xlsx}
#'
#' @param raw_dir A character string to the directory of MS data. 
#' @examples
#' \dontrun{
#' # Supposed that RAW MS files are stored under "~\my_raw"
#' extract_raws("~\\my_raw")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
#' @importFrom tools md5sum
#' @export
extract_raws <- function(raw_dir = NULL) {
  dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
  if (dat_dir == 1) 
    stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  
  if (is.null(raw_dir)) 
    stop("`raw_dir` cannot be `NULL`.", call. = FALSE)

  fns <- names(tools::md5sum(dir(raw_dir, pattern = "\\.raw$", full.names = FALSE)))
  data.frame(Index = seq_along(fns), RAW_File = fns) %T>% 
    write.table(file.path(dat_dir, "raw_list.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
		
	message("RAW MS file names stored in ", file.path(dat_dir, "raw_list.txt"))
}


#' Extract RAW file names from Mascot PSM outputs
#'
#' \code{extract_psm_raws} extracts the RAW file names from the PSM data under
#' the current working directory.
#'
#' @param type Character string indicating the type of PSM.
#' @inheritParams load_expts
#' @examples
#' \donttest{
#' extract_psm_raws(mascot)
#' 
#' extract_psm_raws(maxquant)
#' 
#' extract_psm_raws(spectrum_mill)
#' }
#'
#' @import dplyr tidyr rlang
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
extract_psm_raws <- function(type = c("mascot", "maxquant", "spectrum_mill"), dat_dir = NULL) {
  batchPSMheader_2 <- function(filelist, TMT_plex) {
    df <- readLines(file.path(dat_dir, filelist))
    
    pep_seq_rows <- grep("pep_seq", df)
    
    unassign_hits_row <- grep("Peptide matches not assigned to protein hits", df)
    if (! purrr::is_empty(unassign_hits_row)) {
      psm_end_row <- unassign_hits_row - 2
    } else if (length(pep_seq_rows) > 1) {
      psm_end_row <- pep_seq_rows[2] - 4
    } else {
      psm_end_row <- length(df)
    }
    
    df <- df[pep_seq_rows[1] : psm_end_row]
    df <- gsub("\"---\"", -1, df, fixed = TRUE)
    df <- gsub("\"###\"", -1, df, fixed = TRUE)
    
    df[1] <- paste0(df[1], paste(rep(",", TMT_plex * 4 -2), collapse = ''))
    
    output_prefix <- gsub("\\.csv$", "", filelist)
    writeLines(df, file.path(dat_dir, "PSM\\temp", paste0(output_prefix, "_hdr_rm.csv")))
  }
    
  find_mascot_psmraws <-function() {
    dir.create(file.path(dat_dir, "PSM\\temp"), recursive = TRUE, showWarnings = FALSE)
    
    purrr::walk(filelist, batchPSMheader_2, TMT_plex)
    
    df <- purrr::map(gsub("\\.csv$", "_hdr_rm.csv", filelist), ~ {
      read.delim(file.path(dat_dir, "PSM\\temp", .x), sep = ',', check.names = FALSE, 
                 header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    }) %>% 
      do.call(rbind, .)
    
    r_start <- which(names(df) == "pep_scan_title") + 1
    int_end <- ncol(df)
    if(int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
    
    if (TMT_plex == 16) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131N", "R131C", 
                     "R132N", "R132C", "R133N", "R133C", "R134N")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131N", "I131C", 
                   "I132N", "I132C", "I133N", "I133C", "I134N")
    } else if (TMT_plex == 11) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131N", "R131C")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131N", "I131C")
    } else if (TMT_plex == 10) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131")
    } else if(TMT_plex == 6) {
      col_ratio <- c("R127", "R128", "R129", "R130", "R131")
      col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
    } else {
      col_ratio <- NULL
      col_int <- NULL
    }
    
    df <- df[, "pep_scan_title", drop = FALSE] %>% 
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title))
    
    raws <- unique(df$RAW_File)
    
    if (!purrr::is_empty(raws)) {
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, "mascot_psm_raws.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      message("MS file names stored in ", file.path(dat_dir, "mascot_psm_raws.txt"))
    }
    
    unlink(file.path(dat_dir, "PSM\\temp"), recursive = TRUE, force = TRUE)
  } 
  
  find_sm_psmraws <- function() {
    df <- purrr::map(file.path(dat_dir, filelist), readr::read_delim, delim = ";") %>% 
      dplyr::bind_rows() %>% 
      dplyr::rename(RAW_File = `filename`) %>% 
      dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", RAW_File))
    
    raws <- unique(df$RAW_File)
    
    if (!purrr::is_empty(raws)) {
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, "sm_psm_raws.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      message("MS file names stored in ", file.path(dat_dir, "sm_psm_raws.txt"))
    }
  }
  
  find_mq_psmraws <- function() {
    df <- purrr::map(file.path(dat_dir, filelist), read.csv, 
                     check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
      dplyr::bind_rows() %>% 
      dplyr::rename(RAW_File = `Raw file`) 

    raws <- unique(df$RAW_File)
    
    if (!purrr::is_empty(raws)) {
      data.frame(RAW_File = raws) %>% 
        write.table(file.path(dat_dir, "mq_psm_raws.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      message("MS file names stored in ", file.path(dat_dir, "mq_psm_raws.txt"))
    }
  }
  
  type <- rlang::enexpr(type)
  if (type == rlang::expr(c("mascot", "maxquant", "spectrum_mill"))) {
    type <- "mascot"
  } else {
    type <- rlang::as_string(type)
  }
  
  type <- tolower(type)
  if (type == "ms") type <- "mascot"
  if (type == "mq") type <- "maxquant"
  if (type == "sm") type <- "spectrum_mill"

  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Assign the working directory to variable `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }

  pattern <- switch(type, 
    mascot = "^F[0-9]+\\.csv$", 
    maxquant = "^msms.*\\.txt$",
    spectrum_mill = "^PSMexport.*\\.ssv$", 
    stop("Data type needs to be one of `mascot`, `maxquant` or `spectrum_mill`.", Call. = FALSE)
  )
  
  filelist <- list.files(path = file.path(dat_dir), pattern = pattern)
  if (purrr::is_empty(filelist)) stop("No ", toupper(type), " files(s) under ", dat_dir, call. = FALSE)

  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex(label_scheme)

  switch (type,
    mascot = find_mascot_psmraws(),
    maxquant = find_mq_psmraws(),
    spectrum_mill = find_sm_psmraws(),
  )
}


#' Removes PSM headers
#'
#' \code{rmPSMHeaders} removes the header of PSM from
#' \href{https://http://www.matrixscience.com/}{Mascot} outputs. It also
#' removes the spacer columns in the fields of ratio and intensity values.
#'
#' @return Intermediate PSM table(s).
#'
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>%
rmPSMHeaders <- function () {
	old_opt <- options(max.print = 99999, warn = 0)
	options(max.print = 5000000, warn = 1)
	on.exit(options(old_opt), add = TRUE)

	on.exit(message("Remove PSM headers --- Completed."), add = TRUE)

	filelist = list.files(path = file.path(dat_dir), pattern = "^F[0-9]+\\.csv$")

	if (purrr::is_empty(filelist))
	  stop("No PSM files(s) with `.csv` extension under ", dat_dir, call. = FALSE)

	load(file = file.path(dat_dir, "label_scheme.rda"))
	TMT_plex <- TMT_plex(label_scheme)

	batchPSMheader <- function(filelist, TMT_plex) {
		data_all <- readLines(file.path(dat_dir, filelist))

		pep_seq_rows <- grep("pep_seq", data_all)
		data_header <- data_all[1 : (pep_seq_rows[1] - 1)]
		data_header <- gsub("\"", "", data_header, fixed = TRUE)

		output_prefix <- gsub("\\.csv$", "", filelist)
		write.table(data_header, file.path(dat_dir, "PSM\\cache",
		            paste0(output_prefix, "_header", ".txt")),
		            sep = "\t", col.names = FALSE, row.names = FALSE)
		rm(data_header)
		
		if (length(pep_seq_rows) > 1) {
		  data_queries <- data_all[pep_seq_rows[2] : length(data_all)]
		  data_queries <- gsub("\"", "", data_queries, fixed = TRUE)
		  writeLines(data_queries, file.path(dat_dir, "PSM\\cache", paste0(output_prefix, "_queries.csv")))
		  rm(data_queries)
		}

		unassign_hits_row <- grep("Peptide matches not assigned to protein hits", data_all)
		if (! purrr::is_empty(unassign_hits_row)) {
		  psm_end_row <- unassign_hits_row - 2
		} else if (length(pep_seq_rows) > 1) {
		  psm_end_row <- pep_seq_rows[2] - 4
		} else {
		  psm_end_row <- length(data_all)
		}
	
		local({
		  first_line <- data_all[pep_seq_rows[1]+1]
		  if (grepl("\"134\"", first_line, fixed = TRUE)) {
		    mascot_tmtplex <- 16
		  } else if (grepl("\"131C\"", first_line, fixed = TRUE)) {
		    mascot_tmtplex <- 11
		  } else if (grepl("\"131\"", first_line, fixed = TRUE) && 
		             grepl("\"130C\"", first_line, fixed = TRUE)) {
		    mascot_tmtplex <- 10
		  } else if (grepl("\"131\"", first_line, fixed = TRUE)) {
		    mascot_tmtplex <- 6
		  } else {
		    mascot_tmtplex <- 0
		  }
		  
		  if (mascot_tmtplex != TMT_plex) {
		    warning("Mascot PSMs suggest a TMT ", mascot_tmtplex, "-plex, ", 
		            "which is different to the ", TMT_plex, "-plex in `expt_smry.xlsx`.", 
		            call. = FALSE)
		  }
		})
		
		data_psm <- data_all[pep_seq_rows[1] : psm_end_row]
		data_psm <- gsub("\"---\"", -1, data_psm, fixed = TRUE)
		data_psm <- gsub("\"###\"", -1, data_psm, fixed = TRUE)
		
		# --- for simulated data
		data_psm <- gsub("---", -1, data_psm, fixed = TRUE)
		data_psm <- gsub("###", -1, data_psm, fixed = TRUE)
		
		if (TMT_plex > 0) data_psm[1] <- paste0(data_psm[1], paste(rep(",", TMT_plex * 4 -2), collapse = ''))
		
		writeLines(data_psm, file.path(dat_dir, "PSM\\cache", paste0(output_prefix, "_hdr_rm.csv")))
		rm(data_psm)
	}

	purrr::walk(filelist, batchPSMheader, TMT_plex)
}


#' Add the `pep_seq_mod` field to Mascot PSMs
#' 
#' @inheritParams locate_outliers
#' @inheritParams splitPSM
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
add_mascot_pepseqmod <- function(df, use_lowercase_aa) {
  dat_id <- df$dat_file %>% unique()
  dat_file <- file.path(dat_dir, "PSM\\cache", paste0(dat_id, "_header.txt"))
  stopifnot(length(dat_id)== 1, file.exists(dat_file))
  
  df_header <- readLines(dat_file)

  fixed_mods <- df_header[((grep("Fixed modifications", df_header))[1] + 3) : 
                            ((grep("Variable modifications", df_header))[1] - 2)] %>%
    gsub("\"", "", ., fixed = TRUE) %>%
    data.frame() %>%
    tidyr::separate(".", sep = ",", c("Mascot_abbr", "Description", "Delta_mass"))
  
  var_mods <- df_header[((grep("Variable modifications", df_header))[1] + 3) :
                          ((grep("Search Parameters", df_header))[1] - 2)] %>%
    gsub("\"", "", ., fixed = TRUE) %>%
    data.frame() %>%
    tidyr::separate(".", sep = ",", extra = "drop", c("Mascot_abbr", "Description", "Delta_mass")) %>%
    dplyr::mutate(Filename = gsub("[\\\\\\/\\:\\*\\?\\'\\<\\>\\|]", ".", Description))
  
  if (is.null(df$pep_seq)) stop("column `pep_seq` not found.")

  if (nrow(var_mods) == 0) {
    df$pep_seq_mod <- df$pep_seq
    return(df)
  } 

  if (!use_lowercase_aa) {
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste0(pep_seq, "[", pep_var_mod_pos, "]"))
  } else {
    df$pep_seq_mod <- df$pep_seq
    
    # (1) non terminal modifications
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(!grepl("N-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("C-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            pos_matrix  <- gregexpr(mod, df_sub$pep_var_mod_pos) %>%
              plyr::ldply(., rbind) %>%
              # "-2" for the two characters, "0." ..., in 'pep_var_mod_pos'
              purrr::map(function(x) {x - 2}) %>%
              data.frame(check.names = FALSE)
            
            for (k in 1:ncol(pos_matrix)) {
              rows <- !is.na(pos_matrix[, k])
              locales <- pos_matrix[rows, k]
              
              lowers <- substr(df_sub$pep_seq_mod[rows], locales, locales) %>% tolower()
              substr(df_sub$pep_seq_mod[rows], locales, locales) <- lowers
            }
            
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (2-1) add "_" to sequences from protein N-terminal acetylation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Acetyl (Protein N-term)", Description, fixed = TRUE))
      
      nrow <- nrow(mod_tbl)
      stopifnot(nrow <= 1)
      
      if (nrow == 1) {
        mod <- mod_tbl$Mascot_abbr[1]
        
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)

        df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
        df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
        
        if (nrow(df_sub) > 0) {
          df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0("_", pep_seq_mod))
          df <- rbind(df_rest, df_sub)
        }        
      }

      return(df)
    })
    
    # (2-2) add "_" to sequences from protein C-terminal amidation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Amidated (Protein C-term)", Description, fixed = TRUE))

      nrow <- nrow(mod_tbl)
      stopifnot(nrow <= 1)
      
      if (nrow == 1) {
        mod <- mod_tbl$Mascot_abbr[1]
        
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
        
        df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
        df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
        
        if (nrow(df_sub) > 0) {
          df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "_"))
          df <- rbind(df_rest, df_sub)
        }        
      }

      return(df)
    })
    
    # (3-1) "~" for "(Protein N-term)" other than acetylation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Protein N-term", Description, fixed = TRUE)) %>%
        dplyr::filter(!grepl("Acetyl (Protein N-term)", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0("~", pep_seq_mod))
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (3-2) "~" for "(Protein C-term)" other than amidation
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("Protein C-term", Description, fixed = TRUE)) %>%
        dplyr::filter(!grepl("Amidated (Protein C-term)", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% dplyr::mutate(pep_seq_mod = paste0(pep_seq_mod, "~"))
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (4-1) "^" for peptide "(N-term)"  
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("N-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Protein N-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
  
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = gsub("(^[_~]{0,1})(.)", paste0("\\1", "^", "\\2"), pep_seq_mod)) 
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })
    
    # (4-2) "^" peptide "(C-term)" 
    df <- local({
      mod_tbl <- var_mods %>% 
        dplyr::filter(grepl("C-term", Description, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Protein C-term", Description, fixed = TRUE))
      
      if (nrow(mod_tbl) > 0) {
        var_mods <<- var_mods %>% dplyr::filter(! Mascot_abbr %in% mod_tbl$Mascot_abbr)
          
        for (mod in mod_tbl$Mascot_abbr) {
          df_sub <- df %>% dplyr::filter(grepl(mod, pep_var_mod_pos))
          df_rest <- df %>% dplyr::filter(!grepl(mod, pep_var_mod_pos))
          
          if (nrow(df_sub) > 0) {
            df_sub <- df_sub %>% 
              dplyr::mutate(pep_seq_mod = gsub("(.)([_~]{0,1}$)", paste0("\\1", "^", "\\2"), pep_seq_mod)) 
            
            df <- rbind(df_rest, df_sub)
          }
        }        
      }

      return(df)
    })

    # (5) paste "pep_res_before" and "pep_res_after"
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = "."))
  }

  purrr::walk2(var_mods$Description, var_mods$Filename, ~ {
    try(
      df %>% 
        dplyr::filter(grepl(.x, pep_var_mod, fixed = TRUE)) %>% 
        write.table(file.path(dat_dir, "PSM\\individual_mods", paste0(.y, ".txt")), 
                    sep = "\t", col.names = TRUE, row.names = FALSE)
    )
  })

  return(df)
}


#' Splits PSM tables
#'
#' \code{splitPSM} splits the PSM outputs after \code{rmPSMHeaders()}. It
#' separates PSM data by TMT experiment and LC/MS injection.
#'
#' Arguments \code{group_psm_by} and \code{group_pep_by} are used to update
#' \code{prot_matches_sig} and \code{prot_sequences_sig} after data merge.
#'
#' @param fasta Character string(s) to the name(s) of fasta file(s) with
#'   prepended directory path. The \code{fasta} database(s) need to match those
#'   used in MS/MS ion search. There is no default and users need to provide the
#'   correct file path(s) and name(s).
#' @param entrez Character string(s) to the name(s) of entrez file(s) with
#'   prepended directory path. At the \code{NULL} default, a convenience lookup
#'   is available for species among \code{c("human", "mouse", "rat")}. For other
#'   species, users need to provide the file path(s) and name(s) for the lookup
#'   table(s). See also \code{\link{Uni2Entrez}} and \code{\link{Ref2Entrez}}
#'   for preparing custom entrez files.
#' @param rm_craps Logical; if TRUE,
#'   \href{https://www.thegpm.org/crap/}{cRAP} proteins will be removed.
#'   The default is FALSE.
#' @param rm_krts Logical; if TRUE, keratin entries will be removed. The default
#'   is FALSE.
#' @param annot_kinases Logical; if TRUE, proteins of human or mouse origins
#'   will be annotated with their kinase attributes. The default is FALSE.
#' @param rptr_intco Numeric; the threshold of reporter ion intensity. The
#'   default is 1,000.
#' @param plot_rptr_int Logical; if TRUE, the distributions of reporter-ion
#'   intensities will be plotted. The default is TRUE.
#' @param use_lowercase_aa Logical; if TRUE, modifications in amino acid
#'   residues will be abbreviated with lower-case and/or \code{^_~}. See the
#'   table below for details. The default is TRUE. 
#' @inheritParams annotPSM
#' @import dplyr tidyr seqinr stringr
#' @importFrom magrittr %>%
splitPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, entrez = NULL, 
                     rm_craps = FALSE, rm_krts = FALSE, rptr_intco = 1000, 
                     annot_kinases = FALSE, plot_rptr_int = TRUE, use_lowercase_aa = TRUE, ...) {

	old_opt <- options(max.print = 99999, warn = 0)
	options(max.print = 2000000, warn = 1)
	on.exit(options(old_opt), add = TRUE)
	on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)

	load(file = file.path(dat_dir, "label_scheme_full.rda"))
	load(file = file.path(dat_dir, "label_scheme.rda"))
	load(file = file.path(dat_dir, "fraction_scheme.rda"))

	TMT_plex <- TMT_plex(label_scheme_full)

  filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
                        pattern = "^F[0-9]{6}\\_hdr_rm.csv$")

	if (length(filelist) == 0) stop(paste("No intermediate PSM files under", file.path(dat_dir, "PSM//cache")))

  df <- purrr::map(filelist, ~ {
    data <- read.delim(file.path(dat_dir, "PSM\\cache", .x), sep = ',', check.names = FALSE, 
                       header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    data$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    return(data)
  }) %>% 
    do.call(rbind, .)
  
  r_start <- which(names(df) == "pep_scan_title") + 1
  int_end <- ncol(df) - 1
	if (int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
  
  if (TMT_plex == 16) {
    col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                   "R130N", "R130C", "R131N", "R131C", 
                   "R132N", "R132C", "R133N", "R133C", "R134N")
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C", 
                 "I132N", "I132C", "I133N", "I133C", "I134N")
  } else if (TMT_plex == 11) {
		col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
		               "R130N", "R130C", "R131N", "R131C")
		col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
		             "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
		col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
		               "R130N", "R130C", "R131")
		col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
		             "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
		col_ratio <- c("R127", "R128", "R129", "R130", "R131")
		col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
		col_ratio <- NULL
		col_int <- NULL
	}

	if (TMT_plex > 0) {
		colnames(df)[r_start:(r_start+TMT_plex-2)] <- col_ratio
		colnames(df)[(r_start+TMT_plex-1):(r_start+TMT_plex+TMT_plex-2)] <- col_int
		rm(r_start, int_end, col_ratio, col_int)
	}
  
  # convenience crap removals where their uniprot fasta names ended with "|"
  if (rm_craps) df <- df %>% dplyr::filter(!grepl("\\|.*\\|$", prot_acc))

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `F[...].csv` or `PSM/cache/[...]_hdr_rm.csv`  for `filter_` varargs.")

  # note that `pep_seq` changed from such as MENGQSTAAK to K.MENGQSTAAK.L
  df <- df %>% 
    dplyr::mutate(pep_len = str_length(pep_seq)) %>% 
    split(., .$dat_file, drop = TRUE) %>% 
    purrr::map(add_mascot_pepseqmod, use_lowercase_aa) %>% 
    bind_rows() %>% 
    dplyr::select(-dat_file)

  df <- df %>% 
    dplyr::mutate(prot_acc_orig = prot_acc) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc)) %>% 
    annotPrn(fasta, entrez) %>% 
    dplyr::mutate(prot_acc = prot_acc_orig) %>% 
    dplyr::select(-prot_acc_orig)
  
  prot_matches_sig <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_matches_sig_new = n())
  
  prot_sequences_sig <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_sequences_sig_new = n())
  
  df <- list(df, prot_matches_sig, prot_sequences_sig) %>% 
    purrr::reduce(left_join, by = group_pep_by) %>% 
    dplyr::mutate(prot_matches_sig = prot_matches_sig_new, 
                  prot_sequences_sig = prot_sequences_sig_new) %>%
    dplyr::select(-prot_matches_sig_new, -prot_sequences_sig_new)
  
  rm(prot_matches_sig, prot_sequences_sig)
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc))

  # re-apply craps after annotation
  # 'acc_type' will be NA for entries not found in fasta
  acc_type <- unique(df$acc_type) %>% .[!is.na(.)]

  stopifnot(length(acc_type) == 1)
  
  if (rm_craps) {
    data(package = "proteoQ", prn_annot_crap)

    craps <- prn_annot_crap %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::select(acc_type) %>% 
      unlist()
    
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }

  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) df <- annotKin(df, acc_type)
  
  # `pep_start`, `pep_end` and `gene` will be used for protein percent coverage calculation
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  } 
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(-grep("^prot_", names(.))), 
  )

  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>%
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>%
      # a special case of redundant entries from Mascot
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>% 
			dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` from experimental summary file."),
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")

	for (i in seq_along(df_split)) {
		df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)

		out_fn <- fn_lookup %>%
			dplyr::filter(TMT_inj == names(df_split)[i]) %>%
			dplyr::select(filename) %>%
			unique() %>%
			unlist() %>%
			paste0(., ".csv")

		df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
		
		write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
		
		if (plot_rptr_int & TMT_plex > 0) {
		  df_int <- df_split[[i]] %>% 
		    .[, grepl("^I[0-9]{3}", names(.))]
		  
		  rptr_violin(df = df_int, filepath = file.path(dat_dir, "PSM\\rprt_int\\raw", gsub("\\.csv", "\\.png", out_fn)), 
		              width = 8, height = 8)
		}
	}
}


#' Locates the positions of outliers
#' 
#' @param df A data frame containing the PSM table from database searches.
#' @param range_colRatios The range of columns.
#' @return A data frame.
#' @examples \donttest{locate_outliers(df, 2:3)}
locate_outliers <- function (df, range_colRatios) {
	for(col_index in range_colRatios) {
		counts <- colSums(!is.na(df[col_index]))
		if(counts > 25) df[, col_index] <- Rosner_outliers(df[, col_index])
		else if(counts > 2) df[, col_index] <- Dixon_outliers(df[, col_index])
	}

	return(df)
}


#' Outlier removals with Rosner's method
#' @param x A matrix or data.frame.
Rosner_outliers <- function(x) {

	if (length(unique(x)) < 5) return(x)
	# up to 9-number outliers; may get warnings with NA being an outlier
	gofOutlier_obj <- rosnerTest(as.numeric(x), 9)

	if(gofOutlier_obj$n.outliers > 0) {
		Index <- with(gofOutlier_obj$all.stat, Obs.Num[Outlier == TRUE])
		x[Index] <- NA
	}

	return(x)
}


#' Outlier removals with Dixon's method
#' @inheritParams Rosner_outliers
Dixon_outliers <- function(x) {
	# x = c(0.0000000, 0.0000000, 1.0271542, 0.0000000, 0.2080097)
	# x = c(0.0000000, 0.0000000, NA, 0.0000000, 0.2080097)
	# x = c(0.0000000, 0.0000000, 0.0000000, 0.2080097)
	# x = c(NA, NA, NA, 0.2080097)

	newx <- x[!is.na(x)]
	len_newx <- length(newx)
	uni_newx <- length(unique(newx))

	if (len_newx > 2 & uni_newx > 1) {
		gofOutlier_obj <- dixon.test(as.numeric(x), type = 0)

		while(gofOutlier_obj$p.value < 0.05) {
			if (grepl("^highest", gofOutlier_obj$alternative)) x[which.max(x)] <- NA else x[which.min(x)] <- NA

			newx <- x[!is.na(x)]
			len_newx <- length(newx)
			uni_newx <- length(unique(newx))

			if (len_newx > 2 & uni_newx > 1)
			  gofOutlier_obj <- dixon.test(as.numeric(x), type = 0) else gofOutlier_obj$p.value <- 1
		}
	}

	return (x)
}


#' Outlier removals with Grubbs's method
#' @param type Type for grubbs.test.
#' @inheritParams Rosner_outliers
Grubbs_outliers <- function(x, type = 10) {
	newx <- x[!is.na(x)]
	len_newx <- length(newx)
	uni_newx <- length(unique(newx))

	if (len_newx > 2 & uni_newx > 1) {
		gofOutlier_obj <- grubbs.test(as.numeric(x, type))

		while(gofOutlier_obj$p.value < 0.05) {
			if (grepl("^highest", gofOutlier_obj$alternative)) x[which.max(x)] <- NA else x[which.min(x)] <- NA

			newx <- x[!is.na(x)]
			len_newx <- length(newx)
			uni_newx <- length(unique(newx))

			if (len_newx > 2 & uni_newx > 1)
			  gofOutlier_obj <- grubbs.test(as.numeric(x), type = type) else gofOutlier_obj$p.value <- 1
		}
	}

	return (x)
}


#' Cleans Up PSM results
#'
#' \code{cleanupPSM} removes PSM outliers after \code{splitPSM},
#' \code{splitPSM_mq} or \code{splitPSM_sm}. The outlier removals will be
#' assessed at the basis of per peptide per TMT channel.
#'
#' Dixon's method will be used when \eqn{2 < n \le 25}; Rosner's method will be
#' used when \eqn{n > 25}.
#'
#' @param rm_outliers Logical; if TRUE, PSM outlier removals will be performed
#'   for peptides with more than two identifying PSMs. Dixon's method will be
#'   used when \eqn{2 < n \le 25} and Rosner's method will be used when \eqn{n >
#'   25}. The default is FALSE.
#'
#' @import dplyr tidyr
#' @importFrom stringr str_split
#' @importFrom outliers dixon.test
#' @importFrom EnvStats rosnerTest
cleanupPSM <- function(rm_outliers = FALSE) {
	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 5000000)

	load(file = file.path(dat_dir, "label_scheme.rda"))
	load(file = file.path(dat_dir, "label_scheme_full.rda"))
	TMT_plex <- TMT_plex(label_scheme_full)

	filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
	                      pattern = "^TMT.*LCMS.*\\.csv$")

	for (i in seq_along(filelist)) {
		df <- read.csv(file.path(dat_dir, "PSM\\cache", filelist[i]), check.names = FALSE,
		               header = TRUE, comment.char = "#")

		if (TMT_plex == 0) {
			# lable-free data
		  # re-save ".csv" as ".txt"
			fn <- paste0(gsub(".csv", "", filelist[i]), "_Clean.txt")
			write.table(df, file.path(dat_dir, "PSM\\cache", fn), sep = "\t", col.names = TRUE,
			            row.names = FALSE)
			cat(filelist[i], "processed\n")

			next
		}

		# remove all "-1" ratio rows
		df <- local({
		  N <- sum(grepl("^R[0-9]{3}", names(df)))
  		df <- df %>%
  		  dplyr::mutate(n = rowSums(.[, grep("^R[0-9]{3}", names(.))] == -1, na.rm = TRUE)) %>%
  			dplyr::filter(n != N) %>%
  			dplyr::select(-n)
		})

		channelInfo <- channelInfo(label_scheme, set_idx =
		                  as.integer(gsub("TMTset(\\d+)_.*", "\\1", filelist[i])))

		# add a column of "R126"
		pos_af <- min(grep("^R1[0-9]{2}", names(df)))

		df$R126 <- 1
		df <- cbind.data.frame(df[, 1:(pos_af-1)], R126 = df$R126, df[, (pos_af):(ncol(df)-1)]) %>%
				dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^I[0-9]{3}", names(.) )])) > 0) %>%
				dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
				dplyr::mutate_at(.vars = which(names(.)=="R126")-1+channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) # note that "> 1" not "0"

		if (rm_outliers) {
			dfw_split <- df %>%
				dplyr::select(grep("^I[0-9]{3}", names(.))) %>%
				dplyr::mutate(RM = rowMeans(.[, grep("^I[0-9]{3}", names(.))[channelInfo$labeledChannels]],
				                            na.rm = TRUE)) %>%
				dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ log2(.x/RM)) %>%
				dplyr::select(-c("RM")) %>%
				`colnames<-`(gsub("I", "X", names(.))) %>%
				dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA)) %>%
				dplyr::bind_cols(df[, c("psm_index", "pep_seq")], .) %>%
				split(., .$pep_seq, drop = TRUE)

			range_colRatios <- grep("^X[0-9]{3}", names(dfw_split[[1]]))

			dfw_split <- do.call("rbind", lapply(dfw_split, locate_outliers, range_colRatios)) %>%
					dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA)) %>%
					tidyr::unite(pep_seq_i, pep_seq, psm_index, sep = ":") %>%
					dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, !is.na(.x), 1))

			df <- df %>%
					tidyr::unite(pep_seq_i, pep_seq, psm_index, sep = ":") %>%
					dplyr::left_join(., dfw_split, by = "pep_seq_i") %>%
					tidyr::separate(pep_seq_i, into = c("pep_seq", "psm_index"), sep = ":", remove = TRUE) %>%
					dplyr::select(-c("psm_index"))

			rm(dfw_split, range_colRatios)

			df[, grepl("^I[0-9]{3}", names(df))] <-
			  purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
			              as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
			  dplyr::bind_rows()

			df[, grepl("^R[0-9]{3}", names(df))] <-
			  purrr::map2(as.list(df[, grepl("^R[0-9]{3}", names(df))]),
			              as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
			  dplyr::bind_rows()

			df <- cbind.data.frame(raw_file = df[, c("raw_file")],
					df[, !grepl("^R[0-9]{3}|^I[0-9]{3}|^X[0-9]{3}|^raw_file$", names(df))],
					df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df))]) %>%
					dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) %>% # "> 1" as "R126 == 1"
					dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 + channelInfo$emptyChannels, ~ replace(.x, , 0))
		} else {
			df <- cbind.data.frame(raw_file = df[, c("raw_file")],
				df[, !grepl("^R[0-9]{3}|^I[0-9]{3}|^psm_index$|^raw_file$", names(df))],
				df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df))]) %>%
			  dplyr::mutate_at(.vars = which(names(.)=="I126") - 1 +
				                   channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
				dplyr::mutate_at(.vars = which(names(.)=="R126") - 1 +
				                   channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1)

		}

		fn <- paste0(gsub(".csv", "", filelist[i]), "_Clean.txt")
		write.table(df, file.path(dat_dir, "PSM\\cache", fn), sep = "\t", col.names = TRUE,
		            row.names = FALSE)
		cat(filelist[i], "processed\n")
	}

}


#'Median-centering normalization of PSM data
#'
#'\code{mcPSM} adds fields of \code{log2_R, N_log2_R and N_I} to PSM tables.
#'
#'@param df A data frame containing the PSM table from database searches.
#'@inheritParams channelInfo
#'@import dplyr tidyr purrr
#'@importFrom magrittr %>%
mcPSM <- function(df, set_idx) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  label_scheme_sub <- label_scheme[label_scheme$TMT_Set == set_idx & 
                                     label_scheme$LCMS_Injection == 1, ]
  
  channelInfo <- channelInfo(label_scheme_sub, set_idx)
  
  dfw <- df[rowSums(!is.na(df[, grepl("^R[0-9]{3}", names(df)), drop = FALSE])) > 1, ] %>%
    dplyr::arrange(pep_seq, prot_acc) %>%
    dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 +
                       channelInfo$emptyChannels, ~ replace(.x, , NaN))
  
  col_sample <- grep("^I[0-9]{3}", names(dfw))
  
  if (length(channelInfo$refChannels) > 0) {
    ref_index <- channelInfo$refChannels
  } else {
    ref_index <- channelInfo$labeledChannels
  }
  
  dfw <- sweep(dfw[, col_sample], 1,
               rowMeans(dfw[, col_sample[ref_index], drop = FALSE], na.rm = TRUE), "/") %>%
    log2(.) %>%
    `colnames<-`(gsub("I", "log2_R", names(.)))	%>%
    cbind(dfw, .) %>%
    dplyr::mutate_at(.vars = grep("[I|R][0-9]{3}", names(.)), ~ replace(.x, is.infinite(.), NA))
  
  col_log2Ratio <- grepl("^log2_R[0-9]{3}", names(dfw))
  cf <- apply(dfw[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  
  dfw <- sweep(dfw[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- sweep(dfw[, grepl("^I[0-9]{3}", names(dfw)), drop = FALSE], 2, 2^cf, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- dfw %>%
    reorderCols(endColIndex = grep("[RI][0-9]{3}", names(.)), col_to_rn = "pep_seq_mod") %>%
    na_zeroIntensity()
}


#'Annotates PSM results
#'
#'\code{annotPSM} adds fields of annotation to Mascot PSM tables after
#'\code{rmPSMHeaders}, \code{splitPSM} and \code{cleanupPSM}.
#'
#'@param group_psm_by A character string specifying the method in PSM grouping.
#'  At the \code{pep_seq} default, descriptive statistics will be calculated
#'  based on the same \code{pep_seq} groups. At the \code{pep_seq_mod}
#'  alternative, peptides with different variable modifications will be treated
#'  as different species and descriptive statistics will be calculated based on
#'  the same \code{pep_seq_mod} groups.
#'@param group_pep_by A character string specifying the method in peptide
#'  grouping. At the \code{prot_acc} default, descriptive statistics will be
#'  calculated based on the same \code{prot_acc} groups. At the \code{gene}
#'  alternative, proteins with the same gene name but different accession
#'  numbers will be treated as one group.
#'@param plot_log2FC_cv Logical; if TRUE, the distributions of the CV of peptide
#'  \code{log2FC} will be plotted. The default is TRUE.
#'@param ... Not currently used.
#'@inheritParams load_expts
#'@inheritParams splitPSM
#'@import dplyr tidyr purrr ggplot2 RColorBrewer
#'@importFrom stringr str_split
#'@importFrom tidyr gather
#'@importFrom magrittr %>%
annotPSM <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", 
                     fasta = NULL, expt_smry = "expt_smry.xlsx", 
                     plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 5000000)
  
  hd_fn <- list.files(path = file.path(dat_dir, "PSM\\cache"),
                      pattern = "^F\\d+_header.txt$")
  assign("df_header", readLines(file.path(dat_dir, "PSM\\cache", hd_fn[1])))
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(n_TMT_sets)
  
  for (set_idx in seq_len(n_TMT_sets)) {
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename =
                           do.call('rbind', strsplit(as.character(sublist),
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    # LCMS injections under the same TMT experiment
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#")

      acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
      stopifnot(length(acc_type) == 1)
      
      species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()

      if (TMT_plex > 0) df <- mcPSM(df, set_idx)
      
      df <- df %>% 
        calcSD_Splex(id = group_psm_by, type = "log2_R") %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
        dplyr::right_join(df, by = group_psm_by) %>% 
        dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
        dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 4))
      
      pep_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by)) %>%
        dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
        dplyr::summarise(pep_n_psm = n())
      
      prot_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_psm = n())
      
      prot_n_pep <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_pep = n())
      
      df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^pep_", names(.))), 
      )
      
      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))), 
        df %>% dplyr::select(-grep("^prot_", names(.))), 
      )
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
      )
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int & TMT_plex > 0) {
        df_int <- df %>% .[, grepl("^N_I[0-9]{3}", names(.))]
        rptr_violin(df = df_int, 
                    filepath = file.path(dat_dir, "PSM\\rprt_int\\mc",
                                         paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_rprt.png")), 
                    width = 8, height = 8)
      }

      if (plot_log2FC_cv & TMT_plex > 0) {
        sd_violin(df = df, id = !!group_psm_by, 
                  filepath = file.path(dat_dir, "PSM\\log2FC_cv\\raw", 
                                       paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_sd.png")), 
                  width = 8, height = 8, type = "log2_R", adjSD = FALSE, is_psm = TRUE)
      }
      
    }
    
  }
}


#'Standardization of PSM results
#'
#'\code{normPSM} standardizes
#'\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSM}
#'results from \href{https://en.wikipedia.org/wiki/Tandem_mass_tag}{TMT}
#'experiments.
#'
#'In each primary output file, "\code{...PSM_N.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means if no
#'\code{reference(s)} are present. Values under columns \code{N_log2_R...} are
#'\code{log2_R...} with median-centering alignment. Values under columns
#'\code{I...} are raw \code{reporter-ion intensity} from database searches.
#'Values under columns \code{N_I...} are normalized \code{reporter-ion
#'intensity}. Values under columns \code{sd_log2_R...} are the standard
#'deviation of the \code{log2FC} of peptides from ascribing PSMs. Character
#'strings under \code{pep_seq_mod} denote peptide sequences with applicable
#'variable modifications.
#'
#'\cr \strong{Nomenclature of \code{pep_seq_mod}}:
#'
#'\tabular{ll}{ \emph{Variable modification}   \tab \emph{Abbreviation}\cr
#'Non-terminal \tab A letter from upper to lower case and the flanking residues
#'on the N- or C-terminal side of the peptide separated by a dot, e.g.,
#'\code{-.mtFPEADILLK.S} \cr N-term \tab A hat to the left of a peptide
#'sequence, e.g., \code{K.^QDGTHVVEAVDATHIGK.L} \cr C-term \tab A hat to the
#'right of a peptide sequence, e.g., \code{K.DAYYNLCLPQRPnMI^.-} \cr Acetyl
#'(Protein N-term) \tab A underscore to the left of a peptide sequence, e.g.,
#'\code{-._mAsGVAVSDGVIK.V}. \cr Amidated (Protein C-term) \tab A underscore to
#'the right of a peptide sequence, e.g., \code{K.DAYYNLCLPQRPnMI_.-}. \cr Other
#'(Protein N-term) \tab A tilde to the left of a peptide sequence, e.g.,
#'\code{-.~mAsGVAVSDGVIK.V} \cr Other (Protein C-term) \tab An tilde to the
#'right of a peptide sequence, e.g. \code{K.DAYYNLCLPQRPnMI~.-} \cr }
#'
#'@section \code{Mascot}: End users will export \code{PSM} data from
#'  \href{https://http://www.matrixscience.com/}{Mascot} at a \code{.csv}
#'  format and store them under the file folder indicated by \code{dat_dir}. The
#'  header information should be included during the \code{.csv} export. The
#'  file name(s) should be defaulted by
#'  \href{https://http://www.matrixscience.com/}{Mascot}: starting with
#'  the letter \code{'F'}, followed by digits without space and ended with a
#'  \code{'.csv'} extension \code{(e.g., F004453.csv)}.
#'
#'@section \code{MaxQuant}: End users will copy over \code{msms.txt} file(s)
#'  from \href{https://www.maxquant.org/}{MaxQuant} to the \code{dat_dir}
#'  directory. In the case of multiple \code{msms.txt} files for processing, the
#'  file names need to be compiled in that they all start with \code{'msms'} and
#'  end with a \code{'.txt'} extension.
#'
#'@section \code{Spectrum Mill}: End users will copy over \code{PSMexport.1.ssv}
#'  file(s) from
#'  \href{https://www.agilent.com/en/products/software-informatics/masshunter-suite/masshunter-for-life-science-research/spectrum-mill}{Spectrum Mill} 
#'  to the \code{dat_dir} directory. In the case of multiple
#'  \code{PSMexport} files for processing, the file names need to be compiled in
#'  that they all start with \code{'PSMexport'} and end with a \code{'.ssv'}
#'  extension.
#'
#'@param ... \code{filter_}: Variable argument statements for the filtration of
#'  data rows. Each statement contains to a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{filter_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{pep_expect} is a column key present in \code{Mascot} PSM
#'  exports and \code{filter_psms_at = exprs(pep_expect <= 0.1)} will remove PSM
#'  entries with \code{pep_expect > 0.1}.
#'@inheritParams load_expts
#'@inheritParams splitPSM
#'@inheritParams splitPSM_mq
#'@inheritParams cleanupPSM
#'@inheritParams annotPSM
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@section \code{Variable arguments and data files}: Variable argument (vararg)
#'  statements of \code{filter_} and \code{arrange_} are available in
#'  \code{proteoQ} for flexible filtration and ordering of data rows, via
#'  functions at users' interface. To take advantage of the feature, users need
#'  to be aware of the column keys in input files. As indicated by their names,
#'  \code{filter_} and \code{filter2_} perform row filtration against column
#'  keys from a primary data file, \code{df}, and secondary data file(s),
#'  \code{df2}, respectively. The same correspondence is applicable for
#'  \code{arrange_} and \code{arrange2_} varargs. \cr \cr Users will typically
#'  employ either primary or secondary vararg statements, but not both. In the
#'  more extreme case of \code{gspaMap(...)}, it links \code{\link{prnGSPA}}
#'  findings in \code{df2} to the significance \code{pVals} and abundance fold
#'  changes in \code{df} for volcano plot visualizations by gene sets. The table
#'  below summarizes the \code{df} and the \code{df2} for varargs in
#'  \code{proteoQ}.
#'
#'  \tabular{lllll}{ 
#'  \strong{Utility} \tab \strong{Vararg_} \tab \strong{df} \tab \strong{Vararg2_} \tab \strong{df2} \cr
#'  normPSM \tab filter_ \tab Mascot, \code{F[...].csv}; MaxQuant, \code{msms[...].txt}; 
#'  SM, \code{PSMexport[...].ssv} \tab NA \tab NA \cr 
#'  PSM2Pep \tab NA \tab NA \tab NA \tab NA \cr 
#'  mergePep \tab filter_ \tab \code{TMTset1_LCMSinj1_Peptide_N.txt} \tab NA \tab NA \cr 
#'  standPep \tab slice_ \tab \code{Peptide.txt} \tab NA \tab NA \cr 
#'  Pep2Prn \tab filter_ \tab \code{Peptide.txt} \tab NA \tab NA \cr 
#'  standPrn \tab slice_\tab \code{Protein.txt} \tab NA \tab NA \cr 
#'  pepHist \tab filter_\tab \code{Peptide.tx}t \tab NA \tab NA \cr 
#'  prnHist \tab filter_\tab \code{Protein.txt} \tab NA \tab NA \cr 
#'  pepSig \tab filter_\tab \code{Peptide[_impNA].txt} \tab NA \tab NA \cr 
#'  prnSig \tab filter_\tab \code{Protein[_impNA].txt} \tab NA \tab NA \cr 
#'  pepMDS \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnMDS \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepPCA \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnPCA \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepEucDist \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnEucDist \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepCorr_logFC \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnCorr_logFC \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  pepHM \tab filter_, arrange_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  prnHM \tab filter_, arrange_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  
#'  anal_prnTrend \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  plot_prnTrend \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_Trend_{NZ}[_impNA][...].txt} \cr 
#'  
#'  anal_pepNMF \tab filter_\tab \code{Peptide[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  anal_prnNMF \tab filter_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  plot_pepNMFCon \tab NA \tab NA \tab filter2_\tab \code{[...]Peptide_NMF[...]_consensus.txt} \cr 
#'  plot_prnNMFCon \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_NMF[...]_consensus.txt} \cr 
#'  plot_pepNMFCoef \tab NA \tab NA \tab filter2_\tab \code{[...]Peptide_NMF[...]_coef.txt} \cr 
#'  plot_prnNMFCoef \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_NMF[...]_coef.txt} \cr 
#'  plot_metaNMF \tab filter_, arrange_\tab \code{Protein[_impNA][_pVal].txt} \tab NA \tab NA \cr 
#'  
#'  prnGSPA \tab filter_\tab \code{Protein[_impNA]_pVals.txt} \tab NA \tab NA \cr 
#'  prnGSPAHM \tab NA \tab NA \tab filter2_\tab \code{[...]Protein_GSPA_{NZ}[_impNA]_essmap.txt} \cr 
#'  gspaMap \tab filter_\tab \code{Protein[_impNA]_pVal.txt} \tab filter2_\tab \code{[...]Protein_GSPA_{NZ}[_impNA].txt} \cr 
#'  
#'  anal_prnString \tab filter_\tab \code{Protein[_impNA][_pVals].tx}t \tab NA \tab NA \cr 
#'  }
#'
#'@return Outputs are under the directory of \code{PSM} sub to \code{dat_dir}.
#'  Primary results are in \code{TMTset1_LCMSinj1_PSM_N.txt,
#'  TMTset2_LCMSinj1_PSM_N.txt, ...} The indexes of TMT experiment and LC/MS
#'  injection are indicated in the file names.
#'@example inst/extdata/examples/normPSM_.R
#'@import rlang dplyr purrr ggplot2 RColorBrewer
#'@importFrom stringr str_split
#'@importFrom magrittr %>%
#'@export
normPSM <- function(group_psm_by = c("pep_seq", "pep_seq_mod"), group_pep_by = c("prot_acc", "gene"), 
                    dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx", 
                    fasta = NULL, entrez = NULL, pep_unique_by = "group", corrected_int = TRUE, rm_reverses = TRUE, 
                    rptr_intco = 1000, rm_craps = FALSE, rm_krts = FALSE, rm_outliers = FALSE, 
                    annot_kinases = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                    use_lowercase_aa = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  options(max.print = 2000000)
  on.exit(options(old_opt), add = TRUE)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("normPSM"), add = TRUE)
  
  dots <- rlang::enexprs(...)

  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to the variable first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
    message("Variable `dat_dir` added to the Global Environment.")
  }
  
  if (is.null(fasta)) stop("Path(s) to fasta file(s) not found.", call. = FALSE)
  
  group_psm_by <- rlang::enexpr(group_psm_by)
  if (group_psm_by == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
    group_psm_by <- "pep_seq"
  } else {
    group_psm_by <- rlang::as_string(group_psm_by)
    stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  }
  
  group_pep_by <- rlang::enexpr(group_pep_by)
  if (group_pep_by == rlang::expr(c("prot_acc", "gene"))) {
    group_pep_by <- "prot_acc"
  } else {
    group_pep_by <- rlang::as_string(group_pep_by)
    stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  }

  dir.create(file.path(dat_dir, "PSM\\cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\rprt_int\\raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\rprt_int\\mc"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\purged"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\individual_mods"), recursive = TRUE, showWarnings = FALSE)
  
  if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^F[0-9]+\\.csv$"))) {
    type <- "mascot"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$"))) {
    type <- "mq"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$"))) {
    type <- "sm"
  } else {
    stop("Unknow data type or missing data files.", call. = FALSE)
	}

  pep_unique_by <- rlang::as_string(rlang::enexpr(pep_unique_by))

  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))
  reload_expts()
  
  if (type == "mascot") {
    rmPSMHeaders()
    splitPSM(group_psm_by, group_pep_by, fasta, entrez, rm_craps, rm_krts, rptr_intco, 
             annot_kinases, plot_rptr_int, use_lowercase_aa, ...)
    cleanupPSM(rm_outliers)
    annotPSM(group_psm_by, group_pep_by, fasta, expt_smry, plot_rptr_int, plot_log2FC_cv, ...)
  } else if (type == "mq") {
    splitPSM_mq(group_psm_by, group_pep_by, fasta, entrez, pep_unique_by, corrected_int, 
                rptr_intco, rm_craps, rm_reverses, annot_kinases, plot_rptr_int, ...)
    cleanupPSM(rm_outliers)
		annotPSM_mq(group_psm_by, group_pep_by, fasta, expt_smry, rm_krts, 
		            plot_rptr_int, plot_log2FC_cv, use_lowercase_aa, ...)
  } else if (type == "sm") {
    splitPSM_sm(group_psm_by, group_pep_by, fasta, entrez, rm_craps, rm_krts, rptr_intco, 
                annot_kinases, plot_rptr_int, ...)
    cleanupPSM(rm_outliers)
    annotPSM_sm(group_psm_by, group_pep_by, fasta, expt_smry, rm_krts, plot_rptr_int, 
                plot_log2FC_cv, use_lowercase_aa, ...)
  }
}



#' Calculate peptide data for individual TMT experiments
#' 
#' Argument \code{injn_idx} not currently used.
#' 
#' @param id The same as \code{group_psm_by} in \code{\link{normPSM}}.
#' @param injn_idx Numeric; an index of \code{LCMS_Inj} in metadata \code{label_scheme.xlsx}.
#' @inheritParams mcPSM
#' @inheritParams PSM2Pep
#' @inheritParams annotPSM
#' @inheritParams channelInfo
#' @inheritParams locate_outliers
calcPepide <- function(df, label_scheme, id, method_psm_pep, group_pep_by, set_idx, injn_idx) {
  stopifnot("prot_acc" %in% names(df))
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  channelInfo <- label_scheme %>%
    dplyr::filter(TMT_Set == set_idx) %>%
    channelInfo(set_idx)
  
  df <- df[rowSums(!is.na(df[, grepl("^N_log2_R[0-9]{3}", names(df)), drop = FALSE])) > 0, ] %>%
    dplyr::arrange(!!rlang::sym(id), prot_acc) %>%
    dplyr::select(-grep("^R[0-9]{3}", names(.)))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "prot_hit_num", "prot_family_member", "prot_score", 
      "prot_matches", "prot_sequences", 
      "pep_var_mod", "pep_var_mod_pos", "pep_scan_title", 
      "pep_res_before", "pep_res_after", 
      "raw_file", "pep_query", "pep_summed_mod_pos", "pep_local_mod_pos")))
  
  df <- local({
    col_start <- which(names(df) == "Modifications") + 1
    col_end <- which(names(df) == "Charge") - 1
    
    if (!(is_empty(col_start) | is_empty(col_end))) {
      df <- df %>% dplyr::select(-(col_start : col_end))
    }
    
    return(df)
  })
  
  df <- df %>% 
    dplyr::select(-grep("\\s{1}Probabilities$", names(.))) %>% 
    dplyr::select(-grep("\\s{1}Score\\s{1}Diffs$", names(.))) %>% 
    dplyr::select(-which(names(.) %in% c(
      "Scan number", "Scan index", 
      "Deamidation (N) Probabilities", "Oxidation (M) Probabilities", 
      "Deamidation (N) Score Diffs", "Oxidation (M) Score Diffs", 
      "Acetyl (Protein N-term)", "Deamidation (N)", "Gln->pyro-Glu", "Oxidation (M)", 
      "Modifications", 
      "Fragmentation", "Mass analyzer", "Type", 
      "Scan event number", "Isotope index", "Simple mass error [ppm]", 
      "Retention time", "Delta score", "Score diff", 
      "Localization prob", "Precursor Full ScanNumber", 
      "Precursor Apex Offset", "Precursor Apex Offset Time", 
      "Matches", "Intensities", 
      "Mass Deviations [Da]", "Mass Deviations [ppm]", 
      "Masses", "Number of Matches", 
      "Neutral loss level", "ETD identification type", 
      "Reverse", "All scores", 
      "All sequences", "All modified sequences", 
      "Reporter PIF", "Reporter fraction", 
      "ID", "Protein group IDs", 
      "Peptide ID", "Mod. peptide ID", "Evidence ID", 
      "Length"))) %>% 
    dplyr::select(-grep("site IDs$", names(.)))
  
  df <- df %>% 
    dplyr::select(-which(names(.) %in% c(
      "number", "modifications", 
      "variableSites", "nterm", "previous_aa", "sequence", "next_aa", 
      "cys", "searchCycle", "L/H", "accession_numbers", "entry_name", 
      "matched_parent_mass"))) 
  
  if ("pep_scan_title" %in% names(df)) {
    df <- df %>% 
      dplyr::mutate(pep_scan_title = gsub("\\\\", "~~", pep_scan_title)) %>%
      dplyr::mutate(pep_scan_title = gsub("^File.*~~", "", pep_scan_title))
  }
  
  # summarise log2FC and intensity from the same `set_idx` at one or multiple LCMS injections
  if (method_psm_pep == "mean") {
    df_num <- aggrNums(mean)(df, !!rlang::sym(id), na.rm = TRUE)
  } else if (method_psm_pep == "top.3") {
    df_num <- TMT_top_n(df, !!rlang::sym(id), na.rm = TRUE)
  } else if (method_psm_pep == "weighted.mean") {
    df_num <- TMT_wt_mean(df, !!rlang::sym(id), na.rm = TRUE)
  } else {
    df_num <- aggrNums(median)(df, !!rlang::sym(id), na.rm = TRUE)
  }
  
  df_first <- df %>% 
    dplyr::select(-grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>% 
    med_summarise_keys(id)
  
  df <- list(df_first, df_num) %>%
    purrr::reduce(left_join, by = id)
  
  if (id == "pep_seq_mod") {
    df <- df %>% dplyr::select(-pep_seq)
  } else {
    df <- df %>% dplyr::select(-pep_seq_mod)
  }
  
  df <- cbind.data.frame(df[, !grepl("I[0-9]{3}|log2_R[0-9]{3}", names(df)), drop = FALSE],
                         df[, grepl("I[0-9]{3}", names(df)), drop = FALSE], 
                         df[, grepl("log2_R[0-9]{3}", names(df)), drop = FALSE]) %>%
    dplyr::mutate_at(.vars = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)),
                     list(~ replace(.x, is.infinite(.x), NA)))
  
  df <- local({
    col_r <- grepl("^log2_R[0-9]{3}", names(df))
    cf <- apply(df[, col_r, drop = FALSE], 2, median, na.rm = TRUE)
    df <- cbind(df[, -grep("^N_log2_R[0-9]{3}", names(df))],
                sweep(df[, col_r], 2, cf, "-") %>%
                  `colnames<-`(paste("N", colnames(.), sep="_")))
    
    col_int <- grepl("^I[0-9]{3}", names(df))
    df  <- cbind(df[, -grep("^N_I[0-9]{3}", names(df))],
                 sweep(df[, col_int, drop=FALSE], 2, 2^cf, "/") %>%
                   `colnames<-`(paste("N", colnames(.), sep="_")))
    
    return(df)
  })
  
  df <- df %>% 
    dplyr::mutate(!!group_pep_by := as.character(!!rlang::sym(group_pep_by)))
  
  df <- df %>% 
    calcSD_Splex(group_pep_by) %>% 
    `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
    dplyr::right_join(df, by = group_pep_by) %>% 
    na_zeroIntensity() %>% 
    dplyr::mutate(TMT_Set = set_idx)
  
  return(df)
}


#'Interim peptide tables
#'
#'\code{PSM2Pep} summarizes
#'\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSMs} to
#'peptides by individual TMT experiments and LC/MS series.
#'
#'In general, fields other than \code{log2FC} and \code{intensity} are
#'summarized with median statistics. One exception is with \code{pep_expect} in
#'Mascot or \code{PEP} in MaxQuant where geometric mean is applied.
#'
#'@param method_psm_pep Character string; the method to summarize the
#'  \code{log2FC} and the \code{intensity} of \code{PSMs} by peptide entries.
#'  The descriptive statistics includes \code{c("mean", "median", "top.3",
#'  "weighted.mean")} with \code{median} being the default. The
#'  \code{log10-intensity} of reporter ions at the \code{PSMs} levels will be
#'  the weight when summarizing \code{log2FC} with \code{"top.3"} or
#'  \code{"weighted.mean"}.
#'@param ... Not currently used.
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
#'@family custom database preparation
#'@seealso 
#'  \emph{Custom databases} \cr 
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@return Tables under \code{PSM} folder for each TMT experiment and LC/MS
#'  series: \code{TMTset1_LCMSinj1_PSM_N.txt}, \code{TMTset1_LCMSinj2_PSM_N.txt}...
#'
#'@example inst/extdata/examples/PSM2Pep_.R
#'@import stringr dplyr purrr rlang  magrittr
#'@export
PSM2Pep <- function (method_psm_pep = c("median", "mean", "weighted.mean", "top.3"), ...) {
  old_opt <- options(max.print = 99999, warn = 0)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 2000000, warn = 1)
  
  on.exit(mget(names(formals()), current_env()) %>% c(dots) %>% save_call("PSM2Pep"), add = TRUE)
  
  dots <- rlang::enexprs(...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)

  method_psm_pep <- rlang::enexpr(method_psm_pep)
  if (method_psm_pep == rlang::expr(c("median", "mean", "weighted.mean", "top.3"))) {
    method_psm_pep <- "median"
  } else {
    method_psm_pep <- rlang::as_string(method_psm_pep)
    stopifnot(method_psm_pep %in% c("median", "mean", "weighted.mean", "top.3"))
  }
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  dir.create(file.path(dat_dir, "Peptide\\cache"), recursive = TRUE, showWarnings = FALSE)
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
    reorder_files(n_TMT_sets(label_scheme_full))
  
  purrr::walk(as.list(filelist), ~ {
    fn_prx <- gsub("_PSM_N.txt", "", .x, fixed = TRUE)
    set_idx <- as.integer(gsub(".*TMTset(\\d+)_.*", "\\1", fn_prx))
    injn_idx <- as.integer(gsub(".*LCMSinj(\\d+).*", "\\1", fn_prx))
    fn_pep <- file.path(dat_dir, "Peptide", paste0("TMTset", set_idx, "_LCMSinj", injn_idx, "_Peptide_N.txt"))
    
    df <- read.csv(file.path(dat_dir, "PSM", .x), check.names = FALSE, header = TRUE,
                   sep = "\t", comment.char = "#") %>% 
      dplyr::select(-grep("^sd_log2_R", names(.))) %>% 
      calcPepide(label_scheme = label_scheme, id = !!id, method_psm_pep = method_psm_pep, 
                 group_pep_by = group_pep_by, set_idx = set_idx, injn_idx = injn_idx)
    
    df <- dplyr::bind_cols(
      df %>% dplyr::select(grep("^pep_", names(.))), 
      df %>% dplyr::select(-grep("^pep_", names(.))), 
    )
    
    df <- dplyr::bind_cols(
      df %>% dplyr::select(grep("^prot_", names(.))),
      df %>% dplyr::select(-grep("^prot_", names(.))),
    )
    
    write.table(df, file.path(dat_dir, "Peptide", paste0(fn_prx, "_Peptide_N.txt")), 
                sep = "\t", col.names = TRUE, row.names = FALSE)		
  })
}


#' Add the `pep_seq_mod` field to MaxQuant PSMs
#'
#' @inheritParams splitPSM
#' @inheritParams locate_outliers
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
add_maxquant_pepseqmod <- function(df, use_lowercase_aa) {

  my_tolower <- function(x, ch = "^") {
    locales <- gregexpr(ch, x) %>% .[[1]] %>% `+`(., 1)
    lowers <- map(locales, ~ substr(x, .x, .x)) %>% tolower()
    
    for (i in seq_along(lowers)) {
      substr(x, locales[i], locales[i]) <- lowers[i]
    }
    
    x <- gsub(ch, "", x)
    
    return(x)
  }
  
  if (!use_lowercase_aa) {
    df <- df %>%
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = "."))
  } else {
    # (1) all non-terminal modifications: M(ox) -> m ...
    df <- df %>% 
      tidyr::separate("pep_seq_mod", c("nt", "interior", "ct"), sep = "_") %>% 
      dplyr::mutate(interior = gsub("([A-Z]){1}\\([^\\(\\)]*\\)", paste0("@", "\\1"), interior)) %>% 
      dplyr::mutate_at(vars("interior"), ~ map_chr(.x, my_tolower, "@")) %>% 
      tidyr::unite(pep_seq_mod, nt, interior, ct, sep = ".", remove = TRUE)

    # (2) phospho: pS -> s, pT -> t, pY -> y
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("pS", "s", pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = gsub("pT", "t", pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = gsub("pY", "y", pep_seq_mod)) 

    # (3-1) add "_" to sequences from protein N-terminal acetylation
    df <- local({
      n_ac <- df %>% dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      if (nrow(n_ac) > 0) {
        n_ac <- n_ac %>% dplyr::mutate(pep_seq_mod = gsub("^\\.\\(ac\\)", "_", pep_seq_mod))
        df <- rbind(rest, n_ac)
      }
      
      return(df)
    })
    
    # (3-2) add "_" to sequences from protein C-terminal amidation
    df <- local({
      c_am <- df %>% dplyr::filter(grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      if (nrow(c_am) > 0) {
        c_am <- c_am %>% dplyr::mutate(pep_seq_mod = gsub("\\.\\(am\\)$", "_", pep_seq_mod))
        df <- rbind(rest, c_am)
      }
      
      return(df)
    })
    
    # (4-1) "~" for "(Protein N-term)" other than acetylation
    df <- local({
      n_ac <- df %>% dplyr::filter(grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      other_n <- df %>% 
        dplyr::filter(grepl("Protein N-term", Modifications, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Acetyl (Protein N-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% dplyr::filter(!grepl("Protein N-term", Modifications, fixed = TRUE))
      
      if (nrow(other_n) > 0) {
        other_n <- other_n %>% dplyr::mutate(pep_seq_mod = gsub("^\\.\\([^\\(\\)]*\\)", "~", pep_seq_mod))
        df <- rbind(rest, n_ac, other_n)
      }
      
      return(df)
    })    
    
    # (4-2) "~" for "(Protein C-term)" other than amidation
    df <- local({
      c_am <- df %>% dplyr::filter(grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      other_c <- df %>% 
        dplyr::filter(grepl("Protein C-term", Modifications, fixed = TRUE)) %>% 
        dplyr::filter(!grepl("Amidated (Protein C-term)", Modifications, fixed = TRUE))
      
      rest <- df %>% dplyr::filter(!grepl("Protein C-term", Modifications, fixed = TRUE))
      
      if (nrow(other_c) > 0) {
        other_c <- other_c %>% dplyr::mutate(pep_seq_mod = gsub("\\.\\([^\\(\\)]*\\)$", "~", pep_seq_mod))
        df <- rbind(rest, c_am, other_c)
      }
      
      return(df)
    })
    
    # (5-1) "^" peptide "(N-term)" modification
    df <- local({
      nt <- df %>% dplyr::filter(grepl("(N-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("(N-term)", Modifications, fixed = TRUE))
      
      if (nrow(nt) > 0) {
        nt <- nt %>% 
          dplyr::mutate(pep_seq_mod = gsub("^\\.([_~]{0,1})\\([^\\(\\)]*\\)", paste0("\\1", "^"), pep_seq_mod)) 
        
        df <- rbind(rest, nt)
      }
      
      return(df)
    })
    
    # (5-2) "^" peptide "(C-term)" modification
    df <- local({
      ct <- df %>% dplyr::filter(grepl("(C-term)", Modifications, fixed = TRUE))
      rest <- df %>% dplyr::filter(!grepl("(C-term)", Modifications, fixed = TRUE))
      
      
      if (nrow(ct) > 0) {
        ct <- ct %>% 
          dplyr::mutate(pep_seq_mod = gsub("\\.\\([^\\(\\)]*\\)([_~]{0,1}$)", paste0("^", "\\1"), pep_seq_mod)) 
        
        df <- rbind(rest, ct)
      }
      
      return(df)
    })
    
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("\\.", "", pep_seq_mod))
    
    # (6) other N- or C-terminal modifications better but not named with "N-term" or "C-term": 
    #     (py)C -> c, (gl)Q -> q
    df <- df %>% 
      dplyr::mutate(pep_seq_mod = gsub("(^[_~]{0,1})\\([^\\(\\)]*\\)", paste0("\\1", "^"), pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = gsub("\\([^\\(\\)]*\\)([_~]{0,1}$)", paste0("^", "\\1"), pep_seq_mod)) %>% 
      dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = ".")) %>% 
      dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = "."))
  }
  
  return(df)
}


#'Splits PSM tables
#'
#'\code{splitPSM_mq} splits the PSM outputs after \code{rmPSMHeaders()}. It
#'separates PSM data by TMT experiment and LC/MS injection.
#'
#'Different to \code{splitPSM} used in \code{Mascot} processes,
#'\code{pep_seq_mod}, \code{prot_n_psm}, \code{prot_n_pep} and \code{pep_n_psm}
#'are calculated later in \code{annotPSM_mq}. This is suitable mainly because
#'there is no columns like \code{prot_matches_sig} and \code{prot_sequences_sig}
#'need to be updated after data merging in \code{splitPSM_mq}.
#'
#'@param pep_unique_by A character string for annotating the uniqueness of
#'  peptides in \code{MaxQuant} PSMs. At the \code{group} default, the
#'  uniqueness of peptides is by protein groups. At a more stringent criterion of
#'  \code{protein}, the uniqueness of peptides is by protein entries. A new
#'  column of \code{pep_isunique} with corresponding logical TRUE or FALSE will
#'  be added to the PSM reports.
#'@param corrected_int A logical argument for uses with \code{MaxQuant} data. At
#'  the TRUE default, values under columns "Reporter intensity corrected..." in
#'  \code{MaxQuant} PSM results (\code{msms.txt}) will be used. Otherwise,
#'  "Reporter intensity" values without corrections will be used.
#'@param rm_reverses A logical argument for uses with \code{MaxQuant} data. At
#'  the TRUE default, \code{Reverse} entries will be removed.
#'@inheritParams splitPSM
#'@import dplyr tidyr
#'@importFrom stringr str_split
#'@importFrom magrittr %>%
splitPSM_mq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, entrez = NULL, 
                        pep_unique_by = "group", corrected_int = TRUE, rptr_intco = 1000, 
                        rm_craps = FALSE, rm_reverses = TRUE, annot_kinases = FALSE, 
                        plot_rptr_int = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$")
  
  if (rlang::is_empty(filelist)) {
    stop(paste("No PSM files were found under", file.path(dat_dir), 
               "\nCheck that the names of PSM files start with `msms`."), call. = FALSE)
  }

  df <- purrr::map(file.path(dat_dir, filelist), read.csv, 
                   check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::bind_rows() 
  
  # exception: empty string under `Proteins`
  df <- df %>% dplyr::filter(as.character(.$Proteins) > 0)

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `msms.txt` for `filter_` varargs.")
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::rename(ID = id)
  
  if (pep_unique_by == "group") {
    df <- df %>% 
      dplyr::mutate(pep_isunique = ifelse(grepl(";", `Protein group IDs`), 0, 1))
  } else if (pep_unique_by == "protein") {
    df <- df %>% 
      dplyr::mutate(pep_isunique = ifelse(grepl(";", Proteins), 0, 1))
  }
  
  if (corrected_int) {
    df <- df %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}\\d+$", names(.)))
  } else {
    df <- df %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}corrected\\s{1}\\d+$", names(.)))
  }
  
  if (rm_craps) {
    df <- df %>% dplyr::filter(!grepl("^CON_", .[["Proteins"]]))
  }
  
  if (rm_reverses) {
    df <- df %>% dplyr::filter(.$Reverse != "+")
  }
  
  df <- df %>% 
    `names_pos<-`(grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+$", names(.)), 
                  gsub("TMT-", "I", as.character(TMT_levels)))
  
  if (TMT_plex == 16) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C", 
                 "I132N", "I132C", "I133N", "I133C", "I134N")
  } else if (TMT_plex == 11) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_int <- NULL
  }
  
  df <- local({
    phos_idx <- grep("Phospho (STY) Probabilities", names(df), fixed = TRUE)

    if (!is_empty(phos_idx)) {
      phos <- df %>% 
        dplyr::select(phos_idx) %>% 
        purrr::map(~ str_extract_all(.x, "\\([^()]+\\)")) %>% 
        .[[1]] %>% 
        purrr::map(~ substring(.x, 2, nchar(.x) - 1)) %>% 
        purrr::map(as.numeric) 
      
      phos_max <- phos %>% 
        purrr::map(max, na.rm = TRUE) %>% 
        unlist()

      phos_sort <- phos %>% purrr::map(sort, decreasing = TRUE)
      maxs <- phos_sort %>% purrr::map(`[`, 1) %>% unlist()
      sec_maxs <- phos_sort %>% purrr::map(`[`, 2) %>% unlist()
      
      df <- bind_cols(df, pep_phospho_locprob = phos_max, pep_phospho_locdiff = maxs - sec_maxs) %>% 
        dplyr::mutate(is_locprob_one = equals(1, pep_phospho_locprob)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locdiff"), ~ replace(.x, is_locprob_one, 1)) %>% 
        dplyr::mutate_at(vars("pep_phospho_locprob"), ~ replace(.x, is.infinite(.x), NA)) %>% 
        dplyr::select(-is_locprob_one)
    }
    
    return(df)
  })
  
  if (TMT_plex > 0) {
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) %>% 
      dplyr::rename(
        pep_seq = Sequence, 
        prot_acc = Proteins, 
        RAW_File = `Raw file`, 
      ) 

    df <- df %>% 
      dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc))
  }
  
  acc_type <- parse_acc(df)
  stopifnot(length(acc_type) == 1)
  
  df <- df %>% annotPrn(fasta, entrez)
  if (annot_kinases) df <- df %>% annotKin(acc_type)
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  }
  
  df <- cbind.data.frame(
    df[, grepl("^[a-z]", names(df))], 
    df[, grepl("^[A-Z]", names(df)) & !grepl("^[IR][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[R][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[I][0-9]{3}[NC]*", names(df))]
  )
  
  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>%
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` under the experimental summary file."),
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")

  for (i in seq_along(df_split)) {
    df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)
    
    out_fn <- fn_lookup %>%
      dplyr::filter(TMT_inj == names(df_split)[i]) %>%
      dplyr::select(filename) %>%
      unique() %>%
      unlist() %>%
      paste0(., ".csv")
    
    df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
    
    write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
    
    if (plot_rptr_int & TMT_plex > 0) {
      df_int <- df_split[[i]] %>% 
        .[, grepl("^I[0-9]{3}", names(.))]
      
      rptr_violin(df = df_int, filepath = file.path(dat_dir, "PSM\\rprt_int\\raw", gsub("\\.csv", "\\.png", out_fn)), 
                  width = 8, height = 8)
    }
  }
}


#' Annotates MaxQuant PSM results
#'
#' \code{annotPSM_mq} adds fields of annotation to MaxQuant PSM tables.
#'
#' @inheritParams load_expts
#' @inheritParams annotPSM
#' @inheritParams normPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @param ... Not currently used.
annotPSM_mq <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, expt_smry = "expt_smry.xlsx", 
                        rm_krts = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, use_lowercase_aa = TRUE, ...) {
	
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 5000000)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(., n_TMT_sets)
  
  for (set_idx in seq_len(n_TMT_sets)) {
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename = 
                           do.call('rbind', strsplit(as.character(sublist), '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    # LCMS injections under the same TMT experiment
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") %>%
        dplyr::rename(pep_seq_mod = `Modified sequence`) %>%
        dplyr::select(which(names(.) == "pep_seq_mod"),
                      which(names(.) != "pep_seq_mod"))
      
      acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
      species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()

      if (rm_krts) {
        df <- df %>% 
          dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
      }

      df <- df %>% add_maxquant_pepseqmod(use_lowercase_aa)
      
      df <- df %>% 
        dplyr::mutate(pep_start_discrepancy = ifelse(grepl("Protein N-term", Modifications) & (pep_start > 2), TRUE, FALSE)) %>% 
        dplyr::filter(!pep_start_discrepancy) %>% 
        dplyr::select(-pep_start_discrepancy)
      
      if (TMT_plex > 0) df <- mcPSM(df, set_idx)
			
			df <- df %>% 
			  calcSD_Splex(group_psm_by) %>% 
			  `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
			  dplyr::right_join(df, by = group_psm_by) %>% 
			  dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
			  dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
			  dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
			  dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 4))
			
			pep_n_psm <- df %>%
			  dplyr::select(!!rlang::sym(group_psm_by)) %>%
			  dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
			  dplyr::summarise(pep_n_psm = n())
			
			prot_n_psm <- df %>%
			  dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
			  dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
			  dplyr::summarise(prot_n_psm = n())
			
			prot_n_pep <- df %>%
			  dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
			  dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
			  dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
			  dplyr::summarise(prot_n_pep = n())
			
			df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)

			df <- dplyr::bind_cols(
			  df %>% dplyr::select(grep("^pep_", names(.))), 
			  df %>% dplyr::select(-grep("^pep_", names(.))), 
			)
			
			df <- list(df, prot_n_psm, prot_n_pep) %>%
			  purrr::reduce(left_join, by = group_pep_by)
			
			df <- dplyr::bind_cols(
			  df %>% dplyr::select(grep("^prot_", names(.))), 
			  df %>% dplyr::select(-grep("^prot_", names(.))), 
			)
			
			df <- dplyr::bind_cols(
			  df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
			  df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
			  df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
			)
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int & TMT_plex > 0) {
        df_int <- df %>% .[, grepl("^N_I[0-9]{3}", names(.))]
        rptr_violin(df = df_int, 
                    filepath = file.path(dat_dir, "PSM\\rprt_int\\mc", 
                                         paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_rprt.png")), 
                    width = 8, height = 8)
      }
      
      if (plot_log2FC_cv & TMT_plex > 0) {
        sd_violin(df = df, id = !!group_psm_by, 
                  filepath = file.path(dat_dir, "PSM\\log2FC_cv\\raw", 
                                       paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_sd.png")), 
                  width = 8, height = 8, type = "log2_R", adjSD = FALSE, is_psm = TRUE)
      }
    }
    
  }
}



#' Splits PSM tables from \code{Spectrum Mill}
#'
#' \code{splitPSM_sm} splits the PSM outputs from \code{Spectrum Mill}. It
#' separates PSM data by TMT experiment and LC/MS injection.
#'
#' @inheritParams splitPSM
#' @inheritParams splitPSM_mq
#' @inheritParams normPSM
#' @import dplyr tidyr readr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
splitPSM_sm <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, entrez = NULL, 
                        rm_craps = FALSE, rm_krts = FALSE, rptr_intco = 1000, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)
  
  if (is.null(fasta)) stop("FASTA file(s) not provided.", call. = FALSE)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$")
  
  if (rlang::is_empty(filelist)) 
    stop(paste("No PSM files were found under", file.path(dat_dir), 
               "\nCheck that the names of PSM files start with `PSMexport`."), call. = FALSE)
  
  df <- purrr::map(file.path(dat_dir, filelist), readr::read_delim, delim = ";") %>% 
    dplyr::bind_rows() 
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  message("Primary column keys in `PSMexport[...].ssv` for `filter_` varargs.")
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::rename(prot_acc = accession_number) %>% 
    annotPrn(fasta, entrez)  
  
  if (TMT_plex == 16) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C", 
                 "I132N", "I132C", "I133N", "I133C", "I134N")
  } else if (TMT_plex == 11) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_int <- NULL
  }
  
  if (TMT_plex > 0) {
    df <- df %>% 
      `names<-`(gsub("^TMT_([0-9]{3}[NC]?)", "I\\1", names(.))) %>% 
      dplyr::select(-grep("^I[0-9]{3}[NC]?_[0-9]{3}[NC]?$", names(.))) %>% 
      as.data.frame()
    
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
    
    df <- df %>% 
      dplyr::rename(RAW_File = `filename`) %>% 
      dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", RAW_File)) %>% 
      dplyr::mutate(pep_seq = toupper(sequence)) %>% 
      dplyr::mutate(pep_miss = ifelse(.$next_aa == "(-)", str_count(pep_seq, "[KR]"), str_count(pep_seq, "[KR]") - 1))
  }
  
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)
  
  if (rm_craps) {
    data(package = "proteoQ", prn_annot_crap)
    
    craps <- prn_annot_crap %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::select(acc_type) %>% 
      unlist()
    
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) df <- annotKin(df, acc_type)
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  # M._sequence.c; -._sequence.c; n.sequence.c; -.sequence.c
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  }
  
  df <- df %>% 
    dplyr::mutate(pep_seq = toupper(pep_seq))
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(-grep("^prot_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
  )
  
  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>% 
      # dplyr::select(which(not_all_zero(.))) %>% # don't: empty channels are all NA too
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>% 
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` under the experimental summary file."),
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")
  
  for (i in seq_along(df_split)) {
    df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)
    
    out_fn <- fn_lookup %>%
      dplyr::filter(TMT_inj == names(df_split)[i]) %>%
      dplyr::select(filename) %>%
      unique() %>%
      unlist() %>%
      paste0(., ".csv")
    
    df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
    
    write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
    
    if (plot_rptr_int & TMT_plex > 0) {
      df_int <- df_split[[i]] %>% 
        .[, grepl("^I[0-9]{3}", names(.))]
      
      rptr_violin(df = df_int, filepath = file.path(dat_dir, "PSM\\rprt_int\\raw", gsub("\\.csv", "\\.png", out_fn)), 
                  width = 8, height = 8)
    }
  }
}


#' Annotates \code{Spectrum Mill} PSMs following \code{splitPSM_sm} and
#' \code{cleanupPSM}
#'
#' \code{annotPSM_sm} adds fields of annotation to PSM tables.
#'
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams annotPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @param ... Not currently used.
annotPSM_sm <- function(group_psm_by = "pep_seq", group_pep_by = "prot_acc", fasta = NULL, expt_smry = "expt_smry.xlsx", 
                        rm_krts = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, use_lowercase_aa = TRUE, ...) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  options(max.print = 5000000)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(., n_TMT_sets)
  
  for(set_idx in seq_len(n_TMT_sets)){
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename = 
                           do.call('rbind', strsplit(as.character(sublist),
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") %>%
        dplyr::mutate(pep_seq_mod = pep_seq) %>% 
        dplyr::mutate(pep_seq_mod = as.character(pep_seq_mod)) %>% 
        dplyr::select(which(names(.) == "pep_seq_mod"),
                      which(names(.) != "pep_seq_mod"))
      
      acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
      species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
      
      if (rm_krts) {
        df <- df %>% 
          dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
      }
      
      if (!use_lowercase_aa) {
        df <- df %>% 
          dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
          dplyr::mutate(pep_seq_mod = paste0(pep_seq, "[", variableSites, "]"))
      } else {
        # (1) mods under column `variableSites`
        df <- local({
          if (!is_empty(which(names(df) == "variableSites"))) {
            df_sub <- df %>% dplyr::filter(!is.na(variableSites))
            
            if (nrow(df_sub) > 0) {
              pos_matrix <- df_sub %>% 
                dplyr::select(variableSites) %>% 
                dplyr::mutate(variableSites = as.character(variableSites)) %$% 
                stringr::str_split(.$variableSites, " ") %>% 
                plyr::ldply(., rbind) %>% 
                `names<-`(paste0("mod_", names(.))) %>% 
                purrr::map(~ gsub("[A-z]", "", .)) %>% 
                dplyr::bind_cols() %>% 
                dplyr::mutate_at(.vars = grep("^mod_", names(.)), as.numeric) %>% 
                dplyr::bind_cols(df_sub %>% dplyr::select(pep_seq, pep_start)) %>% 
                dplyr::mutate_at(.vars = grep("^mod_", names(.)), ~ {.x + 1 - pep_start}) %>% 
                dplyr::select(-pep_seq, -pep_start) %>% 
                data.frame(check.names = FALSE)
              
              for (k in seq_along(pos_matrix)) {
                rows <- !is.na(pos_matrix[, k])
                locales <- pos_matrix[rows, k]
                
                lowers <- substr(df_sub$pep_seq_mod[rows], locales, locales) %>% tolower()
                substr(df_sub$pep_seq_mod[rows], locales, locales) <- lowers
              }
              
              df <- dplyr::bind_rows(
                df %>% dplyr::filter(is.na(variableSites)), 
                df_sub
              )            
            }
          }
  
          return(df)
        })
  
        # (2-1) protein n-term acetylation
        if (!is_empty(which(names(df) == "nterm"))) {
          df_sub <- df %>% dplyr::filter(nterm == "Acetyl")
    
          if (nrow(df_sub) > 0) {
            df_sub$pep_seq_mod <- paste0("_", df_sub$pep_seq_mod)
    
            df <- dplyr::bind_rows(
              df %>% dplyr::filter(nterm != "Acetyl"), 
              df_sub
            )
          }        
        }
  
        # (2-2) protein C-terminal amidation
        if (!is_empty(which(names(df) == "cterm"))) {
          # hypothetical; no yet known how it will be named in SM
          df_sub <- df %>% dplyr::filter(grepl("^Amidate", cterm))
          
          if (nrow(df_sub) > 0) {
            df_sub$pep_seq_mod <- paste0(df_sub$pep_seq_mod, "_")
            
            df <- dplyr::bind_rows(
              df %>% dplyr::filter(!grepl("^Amidate", cterm)), 
              df_sub
            )
          }
        }
        
        # (3-1) other protien n-term not yet defined in SM
        # ...
        
        # (3-2) other protien c-term not yet defined in SM
        # ...
        
        # (4-1) peptide "(N-term)" modification: 
        #   only "Pyroglutamic acid (N-termQ)" and under column `variableSites`;
        #   as a result, becomes lower-case `q` after step (1)
        
        # (4-2) peptide "(C-term)" modification not yet defined in SM
        # ...
  
        # paste "pep_res_before" and "pep_res_after"
        df <- df %>%
          dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
          dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = "."))        
      }

      # median centering
      if (TMT_plex > 0) df <- mcPSM(df, set_idx)
      
      # add SD columns
      df <- df %>% 
        calcSD_Splex(group_psm_by) %>% 
        `names<-`(gsub("^log2_R", "sd_log2_R", names(.))) %>% 
        dplyr::right_join(df, by = group_psm_by) %>% 
        dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
        dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^N_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 4))
      
      # add column `pep_n_psm` et al.
      pep_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by)) %>%
        dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
        dplyr::summarise(pep_n_psm = n())
      
      prot_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_psm = n())
      
      prot_n_pep <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_pep = n())
      
      df <- df %>% dplyr::left_join(pep_n_psm, by = group_psm_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^pep_", names(.))), 
      )
      
      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by)
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))), 
        df %>% dplyr::select(-grep("^prot_", names(.))), 
      )
      
      df <- dplyr::bind_cols(
        df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
        df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
      )
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int & TMT_plex > 0) {
        df_int <- df %>% 
          .[, grepl("^N_I[0-9]{3}", names(.))]
        
        rptr_violin(df = df_int, 
                    filepath = file.path(dat_dir, "PSM\\rprt_int\\mc", 
                                         paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_rprt.png")), 
                    width = 8, height = 8)
      }
      
      if (plot_log2FC_cv & TMT_plex > 0) {
        sd_violin(df = df, id = !!group_psm_by, 
                  filepath = file.path(dat_dir, "PSM\\log2FC_cv\\raw", 
                                       paste0(gsub("_PSM_N", "", out_fn[injn_idx, 1]), "_sd.png")), 
                  width = 8, height = 8, type = "log2_R", adjSD = FALSE, is_psm = TRUE)
      }
      
    }
    
  }
}




theme_psm_violin <- theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 24),
    axis.ticks.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, vjust = 0.5, size = 24),
    axis.title.x = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    plot.title = element_text(colour = "black", size  =24, hjust = .5, vjust = .5),
    
    strip.text.x = element_text(size = 18, colour = "black", angle = 0),
    strip.text.y = element_text(size = 18, colour = "black", angle = 90),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    legend.key = element_rect(colour = NA, fill = 'transparent'),
    legend.background = element_rect(colour = NA,  fill = "transparent"),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 18),
    legend.text.align = 0,
    legend.box = NULL
  )


#' Update data groups of `log2_R` etc. by logical matrix `lgl_log2_sd`
#' 
#' @inheritParams info_anal
sd_lgl_cleanup <- function (df) {
  df[, grepl("^log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^N_log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^Z_log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^I[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^N_I[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^N_I[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df[, grepl("^sd_log2_R[0-9]{3}", names(df))] <-
    purrr::map2(as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]),
                as.list(df[, grepl("^lgl_log2_sd[0-9]{3}", names(df))]), `*`) %>%
    dplyr::bind_cols()
  
  df <- df %>% 
    dplyr::select(-grep("^lgl_log2_sd", names(.))) %>% 
    dplyr::filter(rowSums(!is.na(.[, grep("^log2_R[0-9]{3}", names(.) )])) > 0)
}


#' Filter data groups by a CV cut-off
#'
#' \code{purge_by_cv} replaces the data entries at \code{group CV > max_cv} to
#' NA.
#'
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @inheritParams purgePSM
#' @param max_cv Numeric; the cut-off in maximum CV. Values above the threshold
#'   will be replaced with NA. The default is NULL with no data trimming by max
#'   CV.
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
purge_by_cv <- function (df, id, max_cv, keep_ohw = TRUE) {
  if (!is.null(max_cv)) {
    stopifnot(is.numeric(max_cv))

    df_sd_lgl <- df %>% 
      dplyr::select(id, grep("^sd_log2_R[0-9]{3}[NC]*", names(.)))
    
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      # remove duplicated PSMs under the same id
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::filter(!duplicated(.[[id]]))
    } else if (id %in% c("prot_acc", "gene")) {
      # no duplicated id, but protein `sd` under each channel can be a mix of non-NA (one value) and NA
      # this is due to the `wide` joinining of data when forming `Peptide.txt`
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::group_by(!!rlang::sym(id)) %>% 
        dplyr::summarise_all(~ dplyr::first(na.omit(.x)))
    }
    
    if (keep_ohw) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, is.na(.x), -1E-3))
    }
    
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, .x > max_cv, NA)) %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, !is.na(.x), 1)) %>% 
      `names<-`(gsub("^sd_log2_R", "lgl_log2_sd", names(.)))

    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd_lgl, by = id) %>% 
      sd_lgl_cleanup()
  }
  
  return(df)
}


#' Filter data groups by quantiles of CV
#'
#' \code{purge_by_qt} replaces the data entries at \code{CV > quantile} with NA.
#'
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @inheritParams purgePSM
#' @param pt_cv Numeric between 0 and 1; the percentile of CV. Values above the
#'   percentile threshold will be replaced with NA. The default is NULL with no
#'   data trimming by CV percentile.
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
purge_by_qt <- function(df, id, pt_cv = NULL, keep_ohw = TRUE) {
  if (!is.null(pt_cv)) {
    stopifnot(is.numeric(pt_cv))
    if (pt_cv > 1) pt_cv <- pt_cv / 100
    
    df_sd_lgl <- df %>% 
      dplyr::select(id, grep("^sd_log2_R[0-9]{3}[NC]*", names(.)))
    
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::filter(!duplicated(.[[id]]))
    } else if (id %in% c("prot_acc", "gene")) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::group_by(!!rlang::sym(id)) %>% 
        dplyr::summarise_all(~ dplyr::first(na.omit(.x)))
    }
    
    qts <- df_sd_lgl %>% 
      .[, grepl("^sd_log2_R[0-9]{3}", names(.))] %>% 
      purrr::map_dbl(~ quantile(.x, probs = pt_cv, na.rm = TRUE))
    
    if (keep_ohw) {
      df_sd_lgl <- df_sd_lgl %>% 
        dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, is.na(.x), -1E-3))
    }
    
    df_sd_lgl[, grepl("^sd_log2_R[0-9]{3}[NC]*", names(df_sd_lgl))] <- 
      purrr::map2(as.list(df_sd_lgl[, grepl("^sd_log2_R[0-9]{3}[NC]*", names(df_sd_lgl))]), 
                  as.list(qts), ~ {
                    .x[.x > .y] <- NA
                    return(.x)
                  }) %>% 
      dplyr::bind_cols()
    
    df_sd_lgl <- df_sd_lgl %>% 
      dplyr::mutate_at(vars(grep("^sd_log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, !is.na(.x), 1)) %>% 
      `names<-`(gsub("^sd_log2_R", "lgl_log2_sd", names(.)))
    
    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd_lgl, by = id) %>% 
      sd_lgl_cleanup()
  }
  
  return(df)
}


#' Filter data groups by a minimal number of observations (n_obs)
#'
#' \code{purge_by_n} replaces the data entries at \code{group n_obs < min_n} to
#' NA. 
#'
#' @inheritParams prnHist
#' @inheritParams info_anal
#' @param min_n Positive integer. When calling from \code{purgePSM}, peptide
#'   entries in PSM tables with the number of identifying PSMs smaller than
#'   \code{min_n} will be replaced with NA. When calling from \code{purgePep},
#'   protein entries in peptide tables with the number of identifying peptides
#'   smaller than \code{min_n} will be replaced with NA.
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
purge_by_n <- function (df, id, min_n) {
  kept <- df %>%
    dplyr::select(!!rlang::sym(id)) %>%
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= min_n) %>% 
    dplyr::select(id) %>% 
    unlist() %>% 
    as.character()

  df %>% 
    dplyr::filter(.[[id]] %in% kept) 
}


#'Purge PSM data
#'
#'\code{purgePSM} removes \code{peptide} entries from PSM tables by selection
#'criteria. It further plots the distributions of \code{log2FC} by TMT
#'experiments and LC/MS series.
#'
#'The CV of peptides are calculated from contributing PSMs at the basis of per
#'TMT experiment per series of LC/MS. Note that greater CV may be encountered
#'for samples that are more different to reference material(s).
#'
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@inheritParams purge_by_qt
#'@inheritParams plot_prnTrend
#'@param adjSD Not currently used. If TRUE, adjust the standard
#'  deviation in relative to the width of ratio profiles.
#'@param keep_ohw Logical; if TRUE, keep one-hit-wonders with unknown CV. The
#'  default is TRUE.
#'@param ... Additional parameters for plotting: \cr \code{ymax}, the maximum
#'  \eqn{y} at a log2 scale. \cr \code{ybreaks}, the breaks in \eqn{y}-axis at a
#'  log2 scale. \cr \code{width}, the width of plot. \cr \code{height}, the
#'  height of plot. \cr \code{flip_coord}, logical; if TRUE, flip \code{x} and
#'  \code{y} axis.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@example inst/extdata/examples/purgePSM_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
purgePSM <- function (dat_dir = NULL, pt_cv = NULL, max_cv = NULL, adjSD = FALSE, 
                      keep_ohw = TRUE, theme= NULL, ...) {
  on.exit(
    mget(names(formals()), rlang::current_env()) %>% c(enexprs(...)) %>% save_call("purPSM")
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type", "filename", "filepath"), ...)

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  if (!rlang::is_empty(filter_dots)) warning("No data filtration by `filter_` varargs.", call. = FALSE)
  
  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\purged"), recursive = TRUE, showWarnings = FALSE)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  # stopifnot(min_n > 0 & min_n%%1 == 0)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
    reorder_files(n_TMT_sets(label_scheme_full))
  
  dir.create(file.path(dat_dir, "PSM\\Copy"), recursive = TRUE, showWarnings = FALSE)

  for (fn in filelist) {
    file.copy(file.path(dat_dir, "PSM", fn), file.path(dat_dir, "PSM\\Copy", fn))
    
    df <- read.csv(file.path(dat_dir, "PSM", fn), check.names = FALSE, 
                   header = TRUE, sep = "\t", comment.char = "#") %>% 
      # filters_in_call(!!!filter_dots) %>% 
      purge_by_qt(group_psm_by, pt_cv, keep_ohw) %>% 
      purge_by_cv(group_psm_by, max_cv, keep_ohw) %>% 
      dplyr::filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0)
    
    cdns_changed <- any(!is.null(pt_cv), !is.null(max_cv))
    
    if (cdns_changed) {
      pep_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by)) %>%
        dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
        dplyr::summarise(pep_n_psm = n())
      
      prot_n_psm <- df %>%
        dplyr::select(!!rlang::sym(group_pep_by)) %>%
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_psm = n())
      
      prot_n_pep <- df %>%
        dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
        dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
        dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
        dplyr::summarise(prot_n_pep = n())
      
      df <- df %>% 
        dplyr::left_join(pep_n_psm, by = group_psm_by) %>% 
        dplyr::mutate(pep_n_psm.x = pep_n_psm.y) %>% 
        dplyr::rename("pep_n_psm" = "pep_n_psm.x") %>% 
        dplyr::select(-pep_n_psm.y)
  
      df <- list(df, prot_n_psm, prot_n_pep) %>%
        purrr::reduce(left_join, by = group_pep_by) %>% 
        dplyr::mutate(prot_n_psm.x = prot_n_psm.y, prot_n_pep.x = prot_n_pep.y) %>% 
        dplyr::rename(prot_n_psm = prot_n_psm.x, prot_n_pep = prot_n_pep.x) %>% 
        dplyr::select(-prot_n_psm.y, -prot_n_pep.y) %T>% 
        write.table(file.path(dat_dir, "PSM", fn), sep = "\t", col.names = TRUE, row.names = FALSE)      
    }
    
    adjSD <- FALSE
    if (adjSD) {
      dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\purged_adj"), recursive = TRUE, showWarnings = FALSE)
      filepath <- file.path(dat_dir, "PSM\\log2FC_cv\\purged_adj", gsub("_PSM_N.txt", "_sd.png", fn))
    } else {
      filepath <- file.path(dat_dir, "PSM\\log2FC_cv\\purged", gsub("_PSM_N.txt", "_sd.png", fn))
    }
    
    width <- eval(dots$width, env = caller_env())
    height <- eval(dots$height, env = caller_env())
    
    if (is.null(width)) width <- 8 
    if (is.null(height)) height <- 8
    dots <- dots %>% .[! names(.) %in% c("width", "height")]

    quiet_out <- purrr::quietly(sd_violin)(df = df, id = !!group_psm_by, filepath = filepath, 
                                           width = width, height = height, type = "log2_R", adjSD = adjSD, 
                                           is_psm = TRUE, col_select = NULL, col_order = NULL, 
                                           theme = theme, !!!dots)
  }
}


#'Purge peptide data
#'
#'\code{purgePep} removes \code{protein} entries from \code{Peptide.txt} by
#'selection criteria. It further plots the distributions of \code{log2FC}.
#'
#'The CV of proteins under each sample are first calculated from contributing
#'peptides. In the event of multiple series of LC/MS injections, the CV of the
#'same protein from different LC/MS will be summarized by median statistics.
#'
#'The data nullification will be applied column-wisely for all available
#'samples. Argument \code{col_select} is merely used to subsetting samples for
#'the visualization of \code{log2FC} distributions.
#'
#'@inheritParams purgePSM
#'@inheritParams prnHist
#'@inheritParams normPSM
#'@inheritParams purge_by_cv
#'@inheritParams purge_by_n
#'@inheritParams purge_by_qt
#'@inheritParams plot_prnTrend
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
#'@example inst/extdata/examples/purgePep_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'  
#'@export
purgePep <- function (dat_dir = NULL, pt_cv = NULL, max_cv = NULL, adjSD = FALSE, keep_ohw = TRUE, 
                      col_select = NULL, col_order = NULL, filename = NULL, theme= NULL, ...) {
  on.exit(
    mget(names(formals()), rlang::current_env()) %>% c(enexprs(...)) %>% save_call("purPep")
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type", "filename", "filepath"), ...)
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  if (!rlang::is_empty(filter_dots)) warning("No data filtration by `filter_` varargs.", call. = FALSE)

  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  filename <- rlang::enexpr(filename)
  if (is.null(filename)) {
    filename <- "Peptide_sd.png"
    fn_prefix <- "Peptide_sd"
    fn_suffix <- "png"
  } else {
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
  }
  
  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  dir.create(file.path(dat_dir, "Peptide\\log2FC_cv\\purged"), recursive = TRUE, showWarnings = FALSE)
  
  group_psm_by <- match_call_arg(normPSM, group_psm_by)
  group_pep_by <- match_call_arg(normPSM, group_pep_by)

  stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  # stopifnot(min_n > 0 & min_n%%1 == 0)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))

  dir.create(file.path(dat_dir, "Peptide\\Copy"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(dat_dir, "Peptide\\Peptide.txt"), file.path(dat_dir, "Peptide\\Copy\\Peptide.txt"))
  
  fn <- file.path(dat_dir, "Peptide", "Peptide.txt")
  df <- read.csv(fn, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    # filters_in_call(!!!filter_dots) %>% 
    purge_by_qt(group_pep_by, pt_cv, keep_ohw) %>% 
    purge_by_cv(group_pep_by, max_cv, keep_ohw) %>% 
    dplyr::filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0) 

  cdns_changed <- any(!is.null(pt_cv), !is.null(max_cv))
  
  if (cdns_changed) {
    pep_n_psm <- df %>%
      dplyr::select(!!rlang::sym(group_psm_by), pep_n_psm) %>%
      dplyr::group_by(!!rlang::sym(group_psm_by)) %>%
      dplyr::summarise(pep_n_psm = sum(pep_n_psm)) %>% 
      dplyr::arrange(!!rlang::sym(group_psm_by))
    
    prot_n_psm <- df %>% 
      dplyr::select(pep_n_psm, !!rlang::sym(group_pep_by)) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
      dplyr::summarise(prot_n_psm = sum(pep_n_psm))
    
    prot_n_pep <- df %>% 
      dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>% 
      dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
      dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
      dplyr::summarise(prot_n_pep = n())
  
    df <- df %>% 
      dplyr::left_join(pep_n_psm, by = group_psm_by) %>% 
      dplyr::mutate(pep_n_psm.x = pep_n_psm.y) %>% 
      dplyr::rename("pep_n_psm" = "pep_n_psm.x") %>% 
      dplyr::select(-pep_n_psm.y)
    
    df <- list(df, prot_n_psm, prot_n_pep) %>%
      purrr::reduce(left_join, by = group_pep_by) %>% 
      dplyr::mutate(prot_n_psm.x = prot_n_psm.y, prot_n_pep.x = prot_n_pep.y) %>% 
      dplyr::rename(prot_n_psm = prot_n_psm.x, prot_n_pep = prot_n_pep.x) %>% 
      dplyr::select(-prot_n_psm.y, -prot_n_pep.y) %T>% 
      write.table(fn, sep = "\t", col.names = TRUE, row.names = FALSE)    
  }
  
  width <- eval(dots$width, env = caller_env())
  height <- eval(dots$height, env = caller_env())
  
  if (is.null(width)) width <- 8 * n_TMT_sets(label_scheme)
  if (is.null(height)) height <- 8
  dots <- dots %>% .[! names(.) %in% c("width", "height")]
  
  adjSD <- FALSE
  if (adjSD) {
    dir.create(file.path(dat_dir, "Peptide\\log2FC_cv\\purged_adj"), recursive = TRUE, showWarnings = FALSE)
    filepath <- file.path(dat_dir, "Peptide\\log2FC_cv\\purged_adj", filename)
  } else {
    filepath <- file.path(dat_dir, "Peptide\\log2FC_cv\\purged", filename)
  }
  
  quiet_out <- purrr::quietly(sd_violin)(df = df, id = !!group_pep_by, filepath = filepath, 
                                         width = width, height = height, type = "log2_R", 
                                         adjSD = adjSD, is_psm = FALSE, 
                                         col_select = !!col_select, col_order = !!col_order, 
                                         theme = theme, !!!dots)
}





#' Row Variance
#' 
#' @param x A data frame.
#' @param na.rm The same as in \code{mean}.
#' @importFrom magrittr %>%
rowVars <- function (x, na.rm = TRUE) {
		sqr <- function(x) x * x
		n <- rowSums(!is.na(x))
		n[n <= 1] <- NA
		return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}


#' Filter rows by variance quantiles
#' 
#' @inheritParams info_anal
#' @inheritParams gn_rollup
#' @inheritParams prnSig
#' @import dplyr rlang
#' @importFrom magrittr %>%
filterData <- function (df, cols = NULL, var_cutoff = 1E-3) {
  if (is.null(cols)) cols <- 1:ncol(df)
  
  if (length(cols) > 1) {
    df <- df %>% 
      dplyr::select(cols) %>% 
      dplyr::mutate(Variance = rowVars(.)) %>%
      dplyr::mutate(rowname = rownames(df))
    
    Quantile <- quantile(df$Variance, probs = var_cutoff, na.rm = TRUE)
    
    df <- df %>%
      dplyr::filter(Variance >= pmax(Quantile, 1E-3)) %>%
      dplyr::select(-Variance) %>%
      tibble::column_to_rownames()
  }
  
  return(df)
}

#' Prepare formulas
#' 
#' @inheritParams gspaTest
#' @inheritParams model_onechannel
#' @importFrom magrittr %>%
prepFml <- function(formula, label_scheme_sub, ...) {

	# formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]  + (1|TMT_Set) + (1|Duplicate)
	# formula = ~ Term["Ner-V", "Ner_PLUS_PD-PD", "(Ner_PLUS_PD-PD)-(Ner-V)"]
	# formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|TMT_Set)
	# formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|Duplicate)
	# formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]
	# formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|TMT_Set)
	# formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|Duplicate)
	# formula = ~ Term["Ner-PD", "V-PD", "(Ner_PLUS_PD-PD)"] + (1|Duplicate)
	# formula = log2Ratio ~ Term
	# formula = ~ Term[~V] # no interaction terms
	# formula = ~ Term
  
  dots <- rlang::enexprs(...)

	fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
	len <- length(fml)

	key_col <- fml[len] %>% gsub("(.*)\\[\\s*\\~*.*\\].*", "\\1", .)

	label_scheme_sub <- label_scheme_sub %>% dplyr::filter(!is.na(!!sym(key_col)))

	if (grepl("\\[\\s*\\~\\s*", fml[len])) { # formula = ~ Term[ ~ V]
		base <- fml[len] %>% gsub(".*\\[\\s*\\~\\s*(.*)\\]", "\\1", .)
	} else {
		base <- NULL
	}

	if (!is.null(base)) { # formula = ~ Term[~V]
		new_levels <- label_scheme_sub[[key_col]] %>% levels()
		new_levels <- c(new_levels[new_levels == base], new_levels[new_levels != base])
		label_scheme_sub <- label_scheme_sub %>%
			dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col), levels = new_levels))
	} else if (!grepl("\\[", fml[len])) { # formula = log2Ratio ~ Term
		new_levels <- label_scheme_sub[[key_col]] %>% levels() # leveled by the alphebatic order
	} else { # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]
		new_levels <- NULL
	}

	if (!is.null(new_levels)) {
		contrs <- paste(new_levels[-1], new_levels[1], sep = "-")
		elements <- new_levels
	} else {
	  contrs <- fml[len] %>%
	    gsub(".*\\[(.*)\\].*", "\\1", .) %>%
	    gsub("\\\"", "", .) %>%
	    str_split(",\\s*", simplify = TRUE) %>%
	    as.character()

	  new_contrs <- fml[len] %>%
	    gsub("^.*\\[(.*)\\].*", "\\1", .) %>% # may have random terms at the end
	    gsub("\\\"", "", .) %>% 
	    str_split(",\\s*", simplify = TRUE) %>% 
	    gsub("\\s+", "", .) %>% 
	    gsub("<([^>]*?)\\+([^>]*?)>", "<\\1.plus.\\2>", .) %>% 
	    gsub("<([^>]*?)\\-([^>]*?)>", "<\\1.minus.\\2>", .) %>% 
	    gsub("[ <>]+", "", .)

	  new_elements <- new_contrs %>%
	    gsub("/[0-9]", "", .) %>% # (A+B+C)/3-D
	    gsub("[\\(\\)]", "", .) %>%
	    str_split("[\\+\\-]\\s*", simplify = TRUE) %>%
	    as.character() %>%
	    unique() %>% 
	    .[. != ""]
	  
	  elements <- new_elements %>% 
	    gsub(".plus.", "+", ., fixed = TRUE) %>% 
	    gsub(".minus.", "-", ., fixed = TRUE)
	  
	  message("\ncontrs: ", contrs %>% as.character, "\n")
	  message("new_contrs: ", new_contrs %>% as.character, "\n")
	  message("elements: ", elements %>% as.character, "\n")
	  message("new_elements: ", new_elements %>% as.character, "\n\n")		
	}

	label_scheme_sub_sub <- label_scheme_sub %>%
		dplyr::filter(!!sym(key_col) %in% elements) %>%
		dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col)))
	
	if (nrow(label_scheme_sub_sub) == 0) {
	  stop("No samples were found for formula ", formula, 
	       "\nCheck the terms under column ", key_col, call. = FALSE)
	}

	design <- model.matrix(~0+label_scheme_sub_sub[[key_col]]) %>%
		`colnames<-`(levels(label_scheme_sub_sub[[key_col]]))

	new_design_nms <- colnames(design) %>% 
	  gsub("+", ".plus.", ., fixed = TRUE) %>% 
	  gsub("-", ".minus.", ., fixed = TRUE)
	
	new_design <- design %>% 
	  `colnames<-`(new_design_nms)
	
	contr_mat <- makeContrasts(contrasts = new_contrs, levels = data.frame(new_design)) %>% 
	  `colnames<-`(contrs) %>% 
	  `rownames<-`(colnames(design))
	
	rm(new_design_nms, new_design)
	
	random_vars <- fml[len] %>%
		gsub("\\[.*\\]+?", "", .) %>%
		paste("~", .) %>%
		as.formula() %>%
		terms.formula(.) %>%
		attr(., "term.labels") %>%
		.[grepl("\\|", .)] %>%
		gsub("1\\s*\\|\\s*(.*)", "\\1", .)

	message("random_vars: ", random_vars %>% as.character, "\n\n")
	
	return(list(design = design, contr_mat = contr_mat, key_col = key_col, random_vars = random_vars,
							label_scheme_sub_sub = label_scheme_sub_sub))
}


#' Adjusted pVal
#' 
#' @param df_pval A data frame containing pVals
#' @inheritParams prnSig 
#' @importFrom magrittr %>%
my_padj <- function(df_pval, pval_cutoff) {
	df_pval %>%
		purrr::map(~ .x <= pval_cutoff)	%>%
		purrr::map(~ ifelse(!.x, NA, .x)) %>%
		purrr::map2(as.list(df_pval), `*`) %>%
		purrr::map(~ p.adjust(.x, "BH")) %>%
		dplyr::bind_cols() %>%
		`names<-`(gsub("pVal", "adjP", colnames(.))) %>%
		dplyr::mutate(rowname = rownames(df_pval)) %>%
		bind_cols(df_pval, .) %>%
		mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		tibble::column_to_rownames()
}


#' Model summary
#' 
#' @param pvals A data frame of pVal
#' @param log2rs A data frame of log2FC
#' @inheritParams prnSig
#' @importFrom magrittr %>% %$%
lm_summary <- function(pvals, log2rs, pval_cutoff, logFC_cutoff) {
	nms <- rownames(pvals)

	pass_pvals <- pvals %>% purrr::map(~ .x <= pval_cutoff)
	pass_fcs <- log2rs %>% purrr::map(~ abs(.x) >= log2(logFC_cutoff))
	pass_both <- purrr::map2(pass_pvals, pass_fcs, `&`) %>% 
	  purrr::map(~ ifelse(!.x, NA, .x))

	res_padj <- pvals %>%
		purrr::map2(pass_both, `*`) %>%
		purrr::map(~ p.adjust(.x, "BH")) %>%
		data.frame(check.names = FALSE) %>%
		`names<-`(gsub("pVal", "adjP", colnames(.))) %>%
		`rownames<-`(nms) %>%
		tibble::rownames_to_column() %>%
		dplyr::bind_cols(pvals, .) %>%
		mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		tibble::column_to_rownames()

	log2rs <- log2rs %>%
		to_linfc() %>%
		`colnames<-`(gsub("log2Ratio", "FC", names(.))) %>%
		dplyr::bind_cols(log2rs, .) %>%
		dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2) %>%
		`rownames<-`(nms)

	cbind.data.frame(res_padj, log2rs)
}


#' Factorial model formula for interaction terms. All factors under one channel
#' in `label_scheme_sub`
#' 
#' @param formula Language; the formula in linear modeling.
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnSig
#' @importFrom MASS ginv
model_onechannel <- function (df, id, formula, label_scheme_sub, complete_cases, method, var_cutoff, 
                              pval_cutoff, logFC_cutoff, ...) {

	# formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]  + (1|TMT_Set) + (1|Duplicate)
	# formula = ~ Term["Ner-V", "Ner_PLUS_PD-PD", "(Ner_PLUS_PD-PD)-(Ner-V)"]
	# formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|TMT_Set)
	# formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|Duplicate)
	# formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]
	# formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|TMT_Set)
	# formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|Duplicate)
	# formula = ~ Term["Ner-PD", "V-PD", "(Ner_PLUS_PD-PD)"] + (1|Duplicate)
	# formula = log2Ratio ~ Term
	# formula = ~ Term[~V] # no interaction terms
	# formula = ~ Term

	id <- rlang::as_string(rlang::enexpr(id))

	fml_ops <- prepFml(formula, label_scheme_sub, ...)
	contr_mat <- fml_ops$contr_mat
	design <- fml_ops$design
	key_col <- fml_ops$key_col
	random_vars <- fml_ops$random_vars
	label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub

	# keep the name list as rows may drop in filtration
	df_nms <- df %>%
		tibble::rownames_to_column(id) %>%
		dplyr::select(id)

	if (complete_cases) df <- df[complete.cases(df), ]
	
	df <- df %>% 
	  filterData(var_cutoff = var_cutoff) %>% 
	  dplyr::select(as.character(label_scheme_sub_sub$Sample_ID))

	if (length(random_vars) > 0) {
		# only use the first random variable
	  design_random <- label_scheme_sub_sub[[random_vars[1]]] 
		corfit <- duplicateCorrelation(df, design = design, block = design_random)
		fit <- df %>%
			lmFit(design = design, block = design_random, correlation = corfit$consensus) %>%
			contrasts.fit(contr_mat) %>%
			eBayes()
	} else {
		fit <- df %>%
			lmFit(design = design) %>%
			contrasts.fit(contr_mat) %>%
			eBayes()
	}

	print(design)
	print(contr_mat)
	
	# limma
	log2rs <- fit$coefficients %>%
		data.frame(check.names = FALSE) %>%
		`names<-`(paste0("log2Ratio (", names(.), ")"))

	pvals <- fit$p.value %>%
		data.frame(check.names = FALSE) %>%
		`names<-`(paste0("pVal (", names(.), ")"))

	res_lm <- lm_summary(pvals, log2rs, pval_cutoff, logFC_cutoff)

	if (method %in% c("lmer", "lme", "lm")) {
		fml_rhs <- gsub("\\[.*\\]", "", formula) %>% .[length(.)]
		new_formula <- as.formula(paste("log2Ratio", "~", fml_rhs))

		contr_mat_lm <- t(contr_mat) %>%
			MASS::ginv() %>%
			`colnames<-`(colnames(contr_mat)) %>%
			`rownames<-`(rownames(contr_mat))
		names(dimnames(contr_mat_lm)) <- c("Levels", "Contrasts")
		contr_mat_lm <- list(Cdn = contr_mat_lm) %>% `names<-`(key_col)

		smpl_levels <- names(df)
		contr_levels <- attributes(contr_mat_lm[[key_col]])$dimnames$Contrasts
		
		df_lm <- df %>%
			tibble::rownames_to_column(id) %>%
			tidyr::gather(-id, key = Sample_ID, value = log2Ratio) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = smpl_levels)) %>%
			dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", key_col, random_vars)], by = "Sample_ID") %>%
			dplyr::select(which(not_all_NA(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			tidyr::nest()

		if (!purrr::is_empty(random_vars)) {
			res_lm <- df_lm %>%
				dplyr::mutate(
				  model = purrr::map(
				    data, ~ lmerTest::lmer(data = .x, formula = new_formula, contrasts = contr_mat_lm))) %>%
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>%
			  tidyr::unnest(glance, keep_empty = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>%
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>%
				dplyr::mutate(term = gsub(key_col, "", term)) %>% 
			  dplyr::select(-data, -model) %>% 
				dplyr::mutate(term = factor(term, levels = contr_levels)) %>%
				tidyr::spread(term , p.value) %>%
				tibble::column_to_rownames(id) %>%
				`names<-`(paste0("pVal (", names(.), ")")) %>%
				lm_summary(log2rs, pval_cutoff, logFC_cutoff)
		} else {
			res_lm <- df_lm %>%
				dplyr::mutate(model = purrr::map(data, ~ lm(data = .x, formula = new_formula,
				                                            contrasts = contr_mat_lm))) %>%
				dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>%
			  tidyr::unnest(glance, keep_empty = TRUE) %>%	
				dplyr::filter(!grepl("Intercept", term)) %>%
				dplyr::select(-c("std.error", "estimate", "statistic")) %>% 
				dplyr::mutate(term = gsub(key_col, "", term)) %>% 
			  dplyr::select(-data, -model) %>% 
				dplyr::mutate(term = factor(term, levels = contr_levels))	%>% 
				tidyr::spread(term , p.value)	%>% 
				tibble::column_to_rownames(id) %>%
				`names<-`(paste0("pVal (", names(.), ")"))  %>% 
			  lm_summary(log2rs, pval_cutoff, logFC_cutoff)
		}
	}

  df_op <- res_lm %>% 
    tibble::rownames_to_column(id) %>%
    dplyr::right_join(df_nms, by = id) %>%
    tibble::column_to_rownames(var = id)	
}


#' Perform significance tests
#' 
#' @param data_type The type of data being either \code{Peptide} or \code{Protein}.
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnSig
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$%
#' @importFrom broom.mixed tidy
sigTest <- function(df, id, label_scheme_sub, 
                    scale_log2r, complete_cases, impute_na, 
                    filepath, filename, 
										method, var_cutoff, pval_cutoff, logFC_cutoff, 
										data_type, anal_type, ...) {

	id <- rlang::as_string(rlang::enexpr(id))
	method <- rlang::as_string(rlang::enexpr(method))

	dots <- rlang::enexprs(...)
	filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
	dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
	
	non_fml_dots <- dots[!map_lgl(dots, is_formula)]
	dots <- dots[map_lgl(dots, is_formula)]
	
	if (id %in% c("pep_seq", "pep_seq_mod")) {
	  pepSig_formulas <- dots
	  save(pepSig_formulas, file = file.path(dat_dir, "Calls\\pepSig_formulas.rda"))
	  rm(pepSig_formulas)
	} else if (id %in% c("prot_acc", "gene")) {
	  if (rlang::is_empty(dots)) {
	    prnSig_formulas <- dots <- concat_fml_dots()
	  } else {
	    prnSig_formulas <- dots
	  }
	  save(prnSig_formulas, file = file.path(dat_dir, "Calls", "prnSig_formulas.rda"))
	  rm(prnSig_formulas)
	}	

	fn_prefix2 <- ifelse(impute_na, "_impNA_pVals.txt", "_pVals.txt")
	
	quietly_log <- local({
  	dfw <- df %>% 
  	  filters_in_call(!!!filter_dots) %>% 
  	  arrangers_in_call(!!!arrange_dots) %>% 
  	  prepDM(id = !!id, 
  	         scale_log2r = scale_log2r, 
  	         sub_grp = label_scheme_sub$Sample_ID, 
  	         anal_type = anal_type) %>% 
  	  .$log2R
  	
  	# `complete_cases` depends on lm contrasts
	  purrr::map(dots, ~ purrr::quietly(model_onechannel)
	             (dfw, !!id, .x, label_scheme_sub, 
	               complete_cases, method, var_cutoff, 
	               pval_cutoff, logFC_cutoff, !!!non_fml_dots)) 
	})

	if (id %in% c("pep_seq", "pep_seq_mod")) {
	  out_path <- file.path(dat_dir, "Peptide\\Model\\log\\pepSig_log.txt")
	} else if (id %in% c("prot_acc", "gene")) {
	  out_path <- file.path(dat_dir, "Protein\\Model\\log\\prnSig_log.txt")
	}
	
	purrr::map(quietly_log, ~ {
	  .x[[1]] <- NULL
	  return(.x)
	}) %>% 
	  reduce(., `c`) %>% 
	  purrr::walk(., write, out_path, append = TRUE)

	df_op <- purrr::map(quietly_log, `[[`, 1) %>%
	  do.call("cbind", .)

	local({
  	# record the `scale_log2r` status; otherwise, need to indicate it in a way
  	# for example, `_N` or `_Z` in file names
  	dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)	  
	  
	  if (data_type == "Peptide") type <- "pep" else if (data_type == "Protein") type <- "prn"
	  file <- paste0(type, "Sig_imp", ifelse(impute_na, "TRUE", "FALSE"), ".rda")
	  
	  call_pars <- c(scale_log2r = scale_log2r, 
	                 complete_cases = complete_cases, 
	                 impute_na = impute_na) %>% 
	    as.list()
	  
	  save(call_pars, file = file.path(dat_dir, "Calls", file))
	})
	
	suppressWarnings(
  	df_op <- df_op %>%
  	  tibble::rownames_to_column(id) %>% 
  	  dplyr::mutate(!!id := forcats::fct_explicit_na(!!rlang::sym(id))) %>% 
  	  dplyr::right_join(df, ., by = id) %T>% 
  	  write.table(file.path(filepath, paste0(data_type, fn_prefix2)), sep = "\t",
  	              col.names = TRUE, row.names = FALSE)	  
	)

	wb <- createWorkbook("proteoQ")
	addWorksheet(wb, sheetName = "Results")
	openxlsx::writeData(wb, sheet = "Results", df_op)
	saveWorkbook(wb, file = file.path(filepath, paste0(data_type, "_pVals.xlsx")), overwrite = TRUE) 

	invisible(df_op)
}


#'Significance tests of peptide \code{log2FC}
#'
#'\code{pepSig} performs significance tests peptide \code{log2FC}. 
#'
#'@rdname prnSig
#'
#'@import purrr
#'@export
pepSig <- function (scale_log2r = TRUE, impute_na = TRUE, complete_cases = FALSE, 
                    method = "limma",
                    var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), 
                    df = NULL, filepath = NULL, filename = NULL, ...) {
  on.exit({
    mget(names(formals()), current_env()) %>% c(rlang::enexprs(...)) %>% save_call("pepSig")
  }, add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  method <- rlang::as_string(rlang::enexpr(method))
  
  stopifnot(rlang::is_double(var_cutoff), 
            rlang::is_double(pval_cutoff), 
            rlang::is_double(logFC_cutoff))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()
  
  if (!impute_na & method != "limma") {
    impute_na <- TRUE
    warning("Coerce `impute_na = ", impute_na, "` at method = ", method, call. = FALSE)
  }
  
  info_anal(df = !!df, df2 = NULL, id = !!id, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            filepath = !!filepath, filename = !!filename, 
            anal_type = "Model")(method = method, var_cutoff, pval_cutoff, logFC_cutoff, ...)
}


#'Significance tests
#'
#'\code{prnSig} performs significance tests protein \code{log2FC}.
#'
#'In general, special characters of \code{+} or \code{-} should be avoided from
#'contrast terms. Occasionally, such as in biological studies, it may be
#'convenient to use \code{A+B} to denote a condition of combined treatment of
#'\code{A} and \code{B} . In the case, one can put the term(s) containing
#'\code{+} or \code{-} into a pair of pointy brackets. The syntax in the
#'following hypothetical example will compare the effects of \code{A}, \code{B},
#'\code{A+B} and the average of \code{A} and \code{B} to control \code{C}:
#'
#'\code{prnSig(fml = ~ Term["A - C", "B - C", "<A + B> - C", "(A + B)/2 - C"])}
#'
#'Note that \code{<A + B>} stands for one sample and \code{(A + B)} has two
#'samples in it.
#'
#'@inheritParams  prnHist
#'@inheritParams  prnHM
#'@param filename A file name to output results. The default is
#'  \code{Peptide_pVals.txt} for peptides and \code{Protein_pVals} for proteins.
#'@param method Character string; the method of linear modeling. The default is
#'  \code{limma}. At \code{method = lm}, the \code{lm()} in base R will be used
#'  for models without random effects and the \code{\link[lmerTest]{lmer}} will
#'  be used for models with random effects.
#'@param var_cutoff Numeric; the cut-off in the variances of \code{log2FC}.
#'  Entries with variances smaller than the threshold will be excluded from
#'  linear modeling. The default is 1E-3.
#'@param pval_cutoff Numeric; the cut-off in significance \code{pVal}. Entries
#'  with \code{pVals} smaller than the threshold will be excluded from multiple
#'  test corrections. The default is at \code{1} to include all entries.
#'@param logFC_cutoff Numeric; the cut-off in \code{log2FC}. Entries with
#'  absolute \code{log2FC} smaller than the threshold will be removed from
#'  multiple test corrections. The default is at \code{log2(1)} to include all
#'  entries.
#'@param ... User-defined formulas for linear modeling. The syntax starts with a
#'  tilde, followed by the name of an available column key in
#'  \code{expt_smry.xlsx} and square brackets. The contrast groups are then
#'  quoted with one to multiple contrast groups separated by commas. The default
#'  column key is \code{Term} in `expt_smry.xlsx`: \cr \code{~ Term["A - C", "B
#'  - C"]}. \cr Additive random effects are indicated by \code{+ (1|col_key_1) +
#'  (1|col_key_2)}... Currently only a syntax of single contrast are supported
#'  for uses with random effects: \cr \code{~ Term["A - C"] + (1|col_key_1) +
#'  (1|col_key_2)} \cr \cr \code{filter_}: Logical expression(s) for the row
#'  filtration against data in a primary file of \code{Peptide[_impNA].txt} or
#'  \code{Protein[_impNA].txt}. See also \code{\link{normPSM}} for the format of
#'  \code{filter_} statements.
#'@return The primary output is
#'  \code{~\\dat_dir\\Peptide\\Model\\Peptide_pVals.txt} for peptide data or
#'  \code{~\\dat_dir\\Protein\\Model\\Protein_pVals.txt} for protein data. At
#'  \code{impute_na = TRUE}, the corresponding outputs are
#'  \code{Peptide_impNA_pvals.txt} or \code{Protein_impNA_pvals.txt}.
#'
#'@example inst/extdata/examples/prnSig_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnSig <- function (scale_log2r = TRUE, impute_na = TRUE, complete_cases = FALSE, 
                    method = "limma",
                    var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), 
                    df = NULL, filepath = NULL, filename = NULL, ...) {
  on.exit({
    load(file.path(dat_dir, "Calls\\prnSig_formulas.rda"))
    dots <- my_union(rlang::enexprs(...), prnSig_formulas)
    mget(names(formals()), current_env()) %>% c(dots) %>% save_call("prnSig")
  }, add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)

  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  method <- rlang::as_string(rlang::enexpr(method))
    
  stopifnot(rlang::is_double(var_cutoff), 
            rlang::is_double(pval_cutoff), 
            rlang::is_double(logFC_cutoff))

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()
  
  if (!impute_na & method != "limma") {
    impute_na <- TRUE
    warning("Coerce `impute_na = ", impute_na, "` at method = ", method, call. = FALSE)
  }
  
  info_anal(df = !!df, df2 = NULL, id = !!id, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            filepath = !!filepath, filename = !!filename, 
            anal_type = "Model")(method = method, var_cutoff, pval_cutoff, logFC_cutoff, ...)
}

#' Generate simulated UniProt results from RefSeq PSM data 
#' 
#' @inheritParams load_expts
#' @param type The type of PSM data in \code{c("Masoct", "MaxQuant")}.
#' @import rlang magrittr dplyr tidyr stringr
#' @export
simulUniprotPSM <- function(type = "Mascot", dat_dir = NULL) {
  type <- rlang::as_string(rlang::enexpr(type)) %>% 
    tolower()
  
  if (is.null(dat_dir)) dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv),
                                            error = function(e) 1)
  
  if (dat_dir == 1) stop("Set up the working directory first.", call. = FALSE)
  
  switch(type,
         mascot = mascot_refseq2uniprot(dat_dir), 
         maxquant = maxquant_refseq2uniprot(dat_dir), 
         stop("Unknown `type = ", type, "`", call. = FALSE))
}


#' Helper function 
#' 
#' @inheritParams load_expts
mascot_refseq2uniprot <- function(dat_dir) {
  dir.create(file.path(dat_dir, "Mascot_refseq"))
  local({
    filelist <- list.files(path = file.path(dat_dir), pattern = "^F[0-9]+.csv$")
    file.copy(file.path(dat_dir, filelist), file.path(dat_dir, "Mascot_refseq", filelist))
  })
  
  suppressMessages(rmPSMHeaders())
  
  filelist <- list.files(path = file.path(dat_dir, "PSM\\cache"), pattern = "^F[0-9]{6}\\_hdr_rm.csv$")
  filelist_hdr <- list.files(path = file.path(dat_dir, "PSM\\cache"), pattern = "^F[0-9]{6}\\_header.txt$")
  
  purrr::walk2(filelist, filelist_hdr, ~ {
    df <- read.delim(file.path(dat_dir, "PSM\\cache", .x), sep = ',', check.names = FALSE, 
                     header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    df$psm_index <- seq_along(1:nrow(df))
    
    df_acc <- data.frame(prot_acc_orig = df$prot_acc, 
                         database = gsub("(.*)::.*", "\\1", df$prot_acc), 
                         prot_acc = gsub("[1-9]{1}::", "", df$prot_acc), 
                         prot_desc = df$prot_desc, 
                         psm_index = df$psm_index) %>% 
      dplyr::mutate(prot_acc = as.character(prot_acc))
      
    
    df_acc_hs <- df_acc %>% dplyr::filter(grepl("\\[Homo sapiens\\]", prot_desc))
    df_acc_mm <- df_acc %>% dplyr::filter(grepl("\\[Mus musculus\\]", prot_desc))
    df_others <- df_acc %>% dplyr::filter(!grepl("\\[Homo sapiens\\]|\\[Mus musculus\\]", prot_desc))
    
    uniprot_acc_hs <- local({
      accessions <- annot_from_to(abbr_species = "Hs", 
                                  keys = as.character(df_acc_hs$prot_acc), 
                                  from = "REFSEQ", 
                                  to = "UNIPROT") %>% 
        dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
        dplyr::filter(!duplicated(prot_acc))
      
      dplyr::left_join(df_acc_hs, accessions, by = "prot_acc") %>% 
        dplyr::filter(!duplicated(psm_index)) %>% 
        dplyr::mutate(uniprot_acc = paste(database, uniprot_acc, sep = "::")) %>% 
        dplyr::select(uniprot_acc, psm_index) %>% 
        dplyr::rename(prot_acc = uniprot_acc)
    })
    
    uniprot_acc_mm <- local({
      accessions <- annot_from_to(abbr_species = "Mm", 
                                  keys = as.character(df_acc_mm$prot_acc), 
                                  from = "REFSEQ", 
                                  to = "UNIPROT") %>% 
        dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
        dplyr::filter(!duplicated(prot_acc))
      
      dplyr::left_join(df_acc_mm, accessions, by = "prot_acc") %>% 
        dplyr::filter(!duplicated(psm_index)) %>% 
        dplyr::mutate(uniprot_acc = paste(database, uniprot_acc, sep = "::")) %>% 
        dplyr::select(uniprot_acc, psm_index) %>% 
        dplyr::rename(prot_acc = uniprot_acc)
    })
    
    uniprot_acc_others <- df_others %>% 
      dplyr::select(prot_acc_orig, psm_index) %>% 
      dplyr::rename(prot_acc = prot_acc_orig) %>% 
      dplyr::mutate(prot_acc = as.character(prot_acc))
    
    uniprot_acc <- dplyr::bind_rows(
      uniprot_acc_hs,
      uniprot_acc_mm,
      uniprot_acc_others
    ) %>% 
      dplyr::arrange(psm_index) %>% 
      dplyr::select(prot_acc)
    
    df[["prot_acc"]] <- unlist(uniprot_acc)
    df[["psm_index"]] <- NULL
    
    dir.create(file.path(dat_dir, "PSM\\cache\\temp"))
    write.table(df, file = file.path(dat_dir, "PSM\\cache\\temp", .x), sep = ",")
    df <- readLines(file.path(dat_dir, "PSM\\cache\\temp", .x))
    hdr <- readLines(file.path(dat_dir, "PSM\\cache", .y))
    
    write(append(hdr, df), file = file.path(dat_dir, gsub("_hdr_rm", "", .x)))
  })
  
  unlink(file.path(dat_dir, "PSM\\cache\\temp"), recursive = TRUE)
}


#' Helper function 
#' 
#' @inheritParams load_expts
maxquant_refseq2uniprot <- function(dat_dir) {
  dir.create(file.path(dat_dir, "Maxquant_refseq"))
  filelist <- list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$")
  file.copy(file.path(dat_dir, filelist), file.path(dat_dir, "Maxquant_refseq", filelist))
  
  purrr::walk(filelist, ~ {
    df <- read.csv(file.path(dat_dir, .x), sep = "\t", check.names = FALSE, header = TRUE, comment.char = "#")
    df$psm_index <- seq_along(1:nrow(df))
    
    suppressWarnings(
      df_acc <- df$Proteins %>% 
        stringr::str_split(";", simplify = TRUE) %>% 
        `colnames<-`(paste0("prot_acc", 1:ncol(.))) %>% 
        data.frame() %>% 
        dplyr::mutate(psm_index = df$psm_index) %>% 
        tidyr::gather("Index", "prot_acc", -psm_index) %>% 
        dplyr::filter(nchar(prot_acc) != 0) %>% 
        dplyr::select(-Index) %>% 
        dplyr::arrange(psm_index)      
    )

    accessions_hs <- annot_from_to(abbr_species = "Hs", 
                                   keys = as.character(df_acc$prot_acc), 
                                   from = "REFSEQ", 
                                   to = "UNIPROT") %>% 
      dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_acc, ., by = "prot_acc") %>% 
      dplyr::filter(!is.na(uniprot_acc)) %>% 
      dplyr::select(-prot_acc) 
    
    accessions_mm <- annot_from_to(abbr_species = "Mm", 
                                   keys = as.character(df_acc$prot_acc), 
                                   from = "REFSEQ", 
                                   to = "UNIPROT") %>% 
      dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_acc, ., by = "prot_acc") %>% 
      dplyr::filter(!is.na(uniprot_acc)) %>% 
      dplyr::select(-prot_acc) 
    
    accessions <- dplyr::bind_rows(accessions_hs, accessions_mm) %>% 
      dplyr::arrange(psm_index) %>% 
      tidyr::unite(id, c("uniprot_acc", "psm_index"), remove = FALSE) %>% 
      dplyr::filter(!duplicated(id)) %>% 
      dplyr::select(-id) %>% 
      dplyr::group_by(psm_index) %>% 
      dplyr::summarise(uniprot_acc = paste(uniprot_acc, collapse = ";"))
    
    df_others <- df %>% 
      dplyr::filter(!psm_index %in% unique(accessions$psm_index)) %>% 
      dplyr::select(-psm_index) %>% 
      dplyr::mutate(Proteins = as.character(Proteins))
    
    df_hsmm <- df %>% 
      dplyr::filter(psm_index %in% unique(accessions$psm_index)) %>% 
      dplyr::left_join(accessions) %>% 
      dplyr::mutate(Proteins = uniprot_acc) %>% 
      dplyr::select(-uniprot_acc, -psm_index)
    
    write.table(dplyr::bind_rows(df_hsmm, df_others), file = file.path(dat_dir, .x), sep = "\t", row.names = FALSE)
  })
}


#' TMT labeling efficiency
#' 
#' @import dplyr rlang ggplot2
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams cleanupPSM
#' @inheritParams annotPSM
#' @inheritParams normPSM
#' @examples 
#' \donttest{
#' res <- labEffPSM(
#'   fasta = c("~\\proteoQ\\dbs\\fasta\\uniprot\\uniprot_mm_2014_07.fasta"),
#' )
#' }
#' @export
labEffPSM <- function(group_psm_by = c("pep_seq", "pep_seq_mod"), group_pep_by = c("prot_acc", "gene"), 
                      dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx", 
                      fasta = NULL, entrez = NULL, 
                      pep_unique_by = "group", corrected_int = TRUE, rm_reverses = TRUE, 
                      rptr_intco = 1000, rm_craps = FALSE, rm_krts = FALSE, rm_outliers = FALSE, 
                      annot_kinases = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                      use_lowercase_aa = TRUE, ...) {
  
  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to the variable first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
    message("Variable `dat_dir` added to the Global Environment.")
  }
  
  if (is.null(fasta)) stop("Path(s) to fasta file(s) not found.", call. = FALSE)
  
  group_psm_by <- rlang::enexpr(group_psm_by)
  if (group_psm_by == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
    group_psm_by <- "pep_seq"
  } else {
    group_psm_by <- rlang::as_string(group_psm_by)
    stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  }
  
  group_pep_by <- rlang::enexpr(group_pep_by)
  if (group_pep_by == rlang::expr(c("prot_acc", "gene"))) {
    group_pep_by <- "prot_acc"
  } else {
    group_pep_by <- rlang::as_string(group_pep_by)
    stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  }
  
  dir.create(file.path(dat_dir, "PSM\\cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\rprt_int\\raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\rprt_int\\mc"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\log2FC_cv\\purged"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\individual_mods"), recursive = TRUE, showWarnings = FALSE)
  
  if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^F[0-9]+\\.csv$"))) {
    type <- "mascot"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$"))) {
    type <- "mq"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$"))) {
    type <- "sm"
  } else {
    stop("Unknow data type or missing data files.", call. = FALSE)
  }
  
  pep_unique_by <- rlang::as_string(rlang::enexpr(pep_unique_by))
  
  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))
  reload_expts()
  
  rmPSMHeaders()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
                        pattern = "^F[0-9]{6}\\_hdr_rm.csv$")
  
  if (length(filelist) == 0) stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
  df <- purrr::map(filelist, ~ {
    df <- read.delim(file.path(dat_dir, "PSM\\cache", .x), sep = ',', check.names = FALSE, 
                     header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    df$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    
    r_start <- which(names(df) == "pep_scan_title") + 1
    int_end <- ncol(df) - 1
    if (int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
    
    if (TMT_plex == 16) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131N", "R131C", 
                     "R132N", "R132C", "R133N", "R133C", "R134N")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131N", "I131C", 
                   "I132N", "I132C", "I133N", "I133C", "I134N")
    } else if (TMT_plex == 11) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131N", "R131C")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131N", "I131C")
    } else if (TMT_plex == 10) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131")
    } else if(TMT_plex == 6) {
      col_ratio <- c("R127", "R128", "R129", "R130", "R131")
      col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
    } else {
      col_ratio <- NULL
      col_int <- NULL
    }

    if (TMT_plex > 0) {
      colnames(df)[r_start:(r_start+TMT_plex-2)] <- col_ratio
      colnames(df)[(r_start+TMT_plex-1):(r_start+TMT_plex+TMT_plex-2)] <- col_int
      rm(r_start, int_end, col_ratio, col_int)
    }
    
    dat_id <- df$dat_file %>% unique()
    dat_file <- file.path(dat_dir, "PSM\\cache", paste0(dat_id, "_header.txt"))
    stopifnot(length(dat_id)== 1, file.exists(dat_file))
    
    df_header <- readLines(dat_file)
    
    tmt_line <- df_header %>% 
      .[grep("Fixed modifications", .)] %>% 
      .[grep("TMT[0-9]{1}plex", .)]
    
    if (rlang::is_empty(tmt_line)) {
      tmt_type <- "neither"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(K\\)", tmt_line) & grepl("TMT[0-9]{1}plex\\s\\(N-term\\)", tmt_line)) {
      tmt_type <- "both"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(K\\)", tmt_line)) {
      tmt_type <- "k_only"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(N-term\\)", tmt_line)) {
      tmt_type <- "nt_only"
    } 
    
    if (tmt_type == "nt_only") {
      df <- df %>% 
        dplyr::filter(grepl("K", .$pep_seq))
    }
    
    df <- df %>% dplyr::mutate(tmt_type = tmt_type)
    
    return(df)
  }) %>% 
    do.call(rbind, .)
  
  # convenience craps removals where their uniprot afasta names ended with "|"
  if (rm_craps) df <- df %>% dplyr::filter(!grepl("\\|.*\\|$", prot_acc))
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  # note pep_seq: from such as MENGQSTAAK to K.MENGQSTAAK.L
  df <- df %>% 
    dplyr::mutate(pep_len = str_length(pep_seq)) %>% 
    split(., .$dat_file, drop = TRUE) %>% 
    purrr::map(add_mascot_pepseqmod, use_lowercase_aa) %>% 
    bind_rows()
  
  df <- df %>% 
    dplyr::mutate(prot_acc_orig = prot_acc) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc)) %>% 
    annotPrn(fasta, entrez) %>% 
    dplyr::mutate(prot_acc = prot_acc_orig) %>% 
    dplyr::select(-prot_acc_orig)
  
  prot_matches_sig <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_matches_sig_new = n())
  
  prot_sequences_sig <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_sequences_sig_new = n())
  
  df <- list(df, prot_matches_sig, prot_sequences_sig) %>% 
    purrr::reduce(left_join, by = group_pep_by) %>% 
    dplyr::mutate(prot_matches_sig = prot_matches_sig_new, 
                  prot_sequences_sig = prot_sequences_sig_new) %>%
    dplyr::select(-prot_matches_sig_new, -prot_sequences_sig_new)
  
  rm(prot_matches_sig, prot_sequences_sig)
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc))
  
  # re-apply craps after annotation
  # 'acc_type' will be NA for entries not found in fasta
  acc_type <- unique(df$acc_type) %>% .[!is.na(.)]
  
  stopifnot(length(acc_type) == 1)
  
  if (rm_craps) {
    data(package = "proteoQ", prn_annot_crap)
    
    craps <- prn_annot_crap %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::select(acc_type) %>% 
      unlist()
    
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) df <- annotKin(df, acc_type)
  
  # `pep_start`, `pep_end` and `gene` will be used for protein percent coverage calculation
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  } 
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(-grep("^prot_", names(.))), 
  )
  
  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      # dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      # dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>%
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) # %>%
    # dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>% 
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  df_split <- df_split %>% split(., .$RAW_File, drop = TRUE)
  
  dat_files <- df %>% 
    dplyr::select(dat_file, tmt_type) %>% 
    dplyr::filter(!duplicated(dat_file))
  
  lab_eff <- purrr::map2_dbl(df_split, names(df_split), ~ {
    res <- .x %>% 
      dplyr::group_by(dat_file) %>% 
      dplyr::summarize(N = n()) %>% 
      dplyr::right_join(dat_files, by = "dat_file") %>% 
      dplyr::arrange(-N, dat_file) %T>% 
      write.csv(file.path(dat_dir, "PSM\\cache", paste0(.y, ".csv")), row.names = FALSE)
    
    both <- res %>% 
      dplyr::filter(tmt_type == "both") %>% 
      dplyr::mutate(N = 2 *N) %>% 
      .[["N"]]
    
    k_only <- res %>% 
      dplyr::filter(tmt_type == "k_only") %>% 
      .[["N"]]
    
    nt_only <- res %>% 
      dplyr::filter(tmt_type == "nt_only") %>% 
      .[["N"]]
    
    neither <- res %>% 
      dplyr::filter(tmt_type == "neither") %>% 
      dplyr::mutate(N = 2 *N) %>% 
      .[["N"]]
    
    1 - (k_only + nt_only + neither) / (both + k_only + nt_only)
  })
}


#' Makes heat maps
#' 
#' @param df_meta A file name of meta data.
#' @inheritParams prnHM
#' @inheritParams info_anal
#' @inheritParams gspa_colAnnot
#' 
#' @examples
#' \donttest{
#' proteo_hm(
#'   df = Protein_delta.txt, 
#'   id = gene, 
#'   df_meta = hm_meta.xlsx, 
#'   filepath = file.path(dat_dir, "Protein\\Heatmap"), 
#'   filename = "kin_delta.png",
#'   complete_cases = FALSE, 
#'   annot_cols = NULL, 
#'   annot_colnames = NULL, 
#'   annot_rows = c("kin_class"), 
#'   cluster_rows = FALSE, 
#'   xmin = -1, 
#'   xmax = 1, 
#'   xmargin = .1, 
#'   width = 5, 
#'   height = 12,
#'   arrange2_by = exprs(kin_class, gene), 
#' )
#' }
#' 
#' @import stringr dplyr magrittr readr readxl rlang ggplot2 RColorBrewer pheatmap
#' @export
proteo_hm <- function(df = NULL, id = NULL, df_meta = NULL, sample_ids = NULL, 
                      filepath = NULL, filename = NULL, complete_cases = FALSE, 
                      annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                      xmin = -1, xmax = 1, xmargin = .1, ...) {
  
  dir.create(file.path(filepath), recursive = TRUE, showWarnings = FALSE)
  
  id <- rlang::enexpr(id)
  df <- rlang::enexpr(df)
  df_meta <- rlang::enexpr(df_meta)

  if (is.null(df)) stop("Data file `df` cannot be NULL.", call. = FALSE)
  if (is.null(df_meta)) stop("Metadata file `df_meta` cannot be NULL.", call. = FALSE)
  if (is.null(id)) stop("Column key `id` cannot be NULL.", call. = FALSE)
  if (is.null(filepath)) stop("`filepath` cannot be NULL.", call. = FALSE)
  if (is.null(filename)) stop("`filename` cannot be NULL.", call. = FALSE)

  id <- rlang::as_string(id)
  df <- rlang::as_string(df)
  df_meta <- rlang::as_string(df_meta)
  
  message("Use the default sheet name `Sheet1` for metadata Excel.")
  
  df_path <- file.path(filepath, df)
  if (file.exists(df_path)) {
    df <- readr::read_tsv(df_path)
  } else {
    stop("File not found: ", df_path, call. = FALSE)
  }
  
  df_meta_path <- file.path(filepath, df_meta)
  if (file.exists(df_meta_path)) {
    df_meta <- readxl::read_excel(df_meta_path) %>% dplyr::filter(rowSums(!is.na(.)) > 0)
    sample_ids <- df_meta$Sample_ID
  } else {
    if (is.null(sample_ids)) {
      stop("Provide sample IDs under either column `Sample_ID` in the Excel indicated by `df_meta` 
           or in the vector of `sample_ids`.", call. = FALSE)
    }
  }

  dots <- rlang::enexprs(...)
  filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
  arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
  select2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select2_", names(.))]
  dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots, select2_dots)]
  
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

  df <- df %>%
    dplyr::mutate_at(vars(which(names(.) %in% sample_ids)), as.numeric) %>%
    dplyr::mutate_at(vars(which(names(.) %in% sample_ids)), ~ setHMlims(.x, xmin, xmax)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, which(names(.) %in% sample_ids)])) > 0) 
  
  df <- df %>% 
    filters_in_call(!!!filter2_dots) %>% 
    arrangers_in_call(!!!arrange2_dots)
  
  if (nrow(df) == 0) stop("No rows available after data filtratin.", call. = FALSE)
  
  df <- df %>%
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
    annotation_row <- df %>% dplyr::select(annot_rows) %>% data.frame(check.names = FALSE)
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
}

#'Downloads STRING databases
#'
#'@param species Character string; the species. The currently supported species
#'  include \code{human, mouse, rat, fly, cow, dog}. The default is
#'  \code{human}.
#'@param overwrite Logical; if TRUE, overwrite the downloaded database(s). The
#'  default is FALSE.
#'@inheritParams anal_prnString
#'@import rlang dplyr magrittr purrr fs downloader
#'@seealso \code{\link{anal_prnString}} for protein-protein interaction
#'  networks.
#'@export
dl_stringdbs <- function(species = "human", db_path = "~\\proteoQ\\dbs\\string", overwrite = FALSE) {
  species <- rlang::as_string(rlang::enexpr(species))

  if (!fs::dir_exists(db_path)) {
    new_db_path <- fs::path_expand_r(db_path)
    new_db_path2 <- fs::path_expand(db_path)
    
    if (fs::dir_exists(new_db_path)) {
      db_path <- new_db_path
    } else if (fs::dir_exists(new_db_path2)) {
      db_path <- new_db_path2
    } else {
      stop(db_path, " not existed.", call. = FALSE)
    }
    
    rm(new_db_path, new_db_path2)
  }

  urls <- switch(species, 
                 human = c(
                   "9606.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz",
                   "9606.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz", 
                   "9606.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz"
                 ), 
                 mouse = c(
                   "10090.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/10090.protein.links.full.v11.0.txt.gz",
                   "10090.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz",
                   "10090.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/10090.protein.info.v11.0.txt.gz"
                 ), 
                 rat = c(
                   "10116.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/10116.protein.links.full.v11.0.txt.gz",
                   "10116.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/10116.protein.aliases.v11.0.txt.gz", 
                   "10116.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/10116.protein.info.v11.0.txt.gz"
                 ), 
                 fly = c(
                   "7227.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/7227.protein.links.full.v11.0.txt.gz",
                   "7227.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/7227.protein.aliases.v11.0.txt.gz", 
                   "7227.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/7227.protein.info.v11.0.txt.gz"
                 ), 
                 cow = c(
                   "9913.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/9913.protein.links.full.v11.0.txt.gz", 
                   "9913.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/9913.protein.aliases.v11.0.txt.gz",
                   "9913.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/9913.protein.info.v11.0.txt.gz"
                 ), 
                 dog = c(
                   "9612.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/9612.protein.links.full.v11.0.txt.gz", 
                   "9612.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/9612.protein.aliases.v11.0.txt.gz", 
                   "9612.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/9612.protein.info.v11.0.txt.gz"
                 ), 
                 stop("Unknown `species`.", Call. = FALSE)
  )
  
  abbr_species <- sp_lookup(species) 
  taxid <- taxid_lookup(species)
  
  db_path2 <- file.path(db_path, abbr_species)
  dir.create(file.path(db_path2, "zip"), recursive = TRUE, showWarnings = FALSE)

  for(i in seq_along(urls)) {
    url <- urls[[i]]
    fn_zip <- names(urls)[[i]]
    fn_tsv <- gsub("\\.gz$", "", fn_zip)
    filepath <- file.path(db_path2, fn_zip)
    
    if ((!file.exists(filepath)) | overwrite) {
      downloader::download(url, filepath, mode = "wb")
      con <- gzfile(path.expand(filepath))
      
      if (grepl("protein\\.links", filepath)) {
        read.csv(con, sep = " ", check.names = FALSE, header = TRUE, comment.char = "#") %>% 
          write.table(file.path(db_path2, fn_tsv), sep = "\t", col.names = TRUE, row.names = FALSE)
      } else if (grepl("protein\\.info", filepath)) {
        temp <- readLines(con)
        for (idx in seq_along(temp)) {
          temp[idx] <- gsub("^(.*)\t[^\t].*$", "\\1", temp[idx])
        }
        temp[1] <- "protein_external_id\tpreferred_name\tprotein_size"
        writeLines(temp, file.path(db_path2, fn_tsv))
        
        rm(temp, idx)
      } else if (grepl("protein\\.alias", filepath)) {
        temp <- read.csv(con, sep = "\t", check.names = FALSE, header = TRUE, comment.char = "#")
        
        col_nms <- c("string_protein_id", "alias", "source")
        first_row <- names(temp) %>% 
          data.frame() %>% 
          t() %>% 
          `colnames<-`(col_nms)
        
        temp %>% 
          `colnames<-`(col_nms) %>% 
          dplyr::mutate_all(as.character) %>% 
          rbind(first_row, .) %>% 
          write.table(file.path(db_path2, fn_tsv), sep = "\t", col.names = TRUE, row.names = FALSE)
        
        rm(temp)
      }
      
      # close(con)
    }
  }

}


#'Annotates protein STRING ids by species
#'
#'@inheritParams info_anal
#'@inheritParams dl_stringdbs
#'@inheritParams anal_prnString
#'@import rlang dplyr purrr magrittr
annot_stringdb <- function(species, df, db_path, id, score_cutoff, 
                           filepath = NULL, filename = NULL, ...) {
  abbr_species <- sp_lookup(species) 
  taxid <- taxid_lookup(species)
  dl_stringdbs(!!species)
  
  db_path2 <- file.path(db_path, abbr_species)
  filelist_info <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.info", ".*.txt$"))
  filelist_link <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.links", ".*.txt$"))
  filelist_alias <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.aliases", ".*.txt$"))
  
  df <- df %>% dplyr::mutate(!!id := as.character(!!rlang::sym(id)))

  prn_info <- read.csv(file.path(db_path2, filelist_info), sep = "\t", check.names = FALSE, 
                       header = TRUE, comment.char = "#") %>% 
    dplyr::select(protein_external_id, preferred_name) %>% 
    dplyr::rename(!!id := preferred_name) %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id)))
  
  prn_links <- read.csv(file.path(db_path2, filelist_link), sep = "\t", 
                        check.names = FALSE, header = TRUE, comment.char = "#") %>% 
    dplyr::mutate(protein1 = as.character(protein1), protein2 = as.character(protein2))
  
  prn_alias <- read.csv(file.path(db_path2, filelist_alias), sep = "\t", 
                        check.names = FALSE, header = TRUE, comment.char = "#")
  
  prn_alias_sub <- prn_alias %>% 
    dplyr::filter(.$alias %in% df[[id]]) %>% 
    dplyr::filter(!duplicated(alias)) %>% 
    dplyr::select(-source) %>% 
    dplyr::rename(!!id := alias) %>% 
    dplyr::rename(protein_external_id = string_protein_id) %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id)))

  if (id == "gene") {
    string_map <- df %>%
      dplyr::select(id) %>% 
      dplyr::left_join(prn_info, by = id) %>% 
      dplyr::filter(!is.na(protein_external_id))
  } else {
    string_map <- df %>%
      dplyr::select(id) %>% 
      dplyr::left_join(prn_alias_sub, by = id) %>% 
      dplyr::filter(!is.na(protein_external_id))
  }
  string_map <- string_map %>% 
    dplyr::mutate(protein_external_id = as.character(protein_external_id))
  
  prn_links_sub <- prn_links %>% 
    dplyr::filter(protein1 %in% string_map$protein_external_id) %>% 
    dplyr::left_join(string_map, by = c("protein1" = "protein_external_id")) %>% 
    dplyr::rename(node1 = !!rlang::sym(id)) %>% 
    dplyr::left_join(string_map, by = c("protein2" = "protein_external_id")) %>% 
    dplyr::rename(node2 = !!rlang::sym(id))
  
  first_four <- c("node1", "node2", "protein1", "protein2")
  ppi <- dplyr::bind_cols(
    prn_links_sub[, first_four], 
    prn_links_sub[, !names(prn_links_sub) %in% first_four]
  ) %>% 
    dplyr::filter(!is.na(node1), !is.na(node2)) %>% 
    `names_pos<-`(1:4, c("#node1", "node2", "node1_external_id", "node2_external_id")) %>% 
    dplyr::filter(combined_score > score_cutoff)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  def_pre <- "^Protein_String_[NZ]"
  if (grepl(def_pre, fn_prefix)) {
    outnm_ppi <- paste0(fn_prefix, "_", species, "_ppi.tsv")
    outnm_expr <- paste0(fn_prefix, "_", species, "_expr.tsv")
  } else if (grepl(species, fn_prefix)) {
    outnm_ppi <- paste0(fn_prefix, "_ppi.tsv")
    outnm_expr <- paste0(fn_prefix, "_expr.tsv")
  } else {
    outnm_ppi <- paste0(fn_prefix, "_", species, "_ppi.tsv")
    outnm_expr <- paste0(fn_prefix, "_", species, "_expr.tsv")
  }
  
  write.table(ppi, file.path(dat_dir, "Protein\\String", outnm_ppi), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # expression data
  gns <- c(ppi[["#node1"]], ppi[["node2"]]) %>% unique()
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(id), 
    df %>% dplyr::select(grep("pVal|adjP|log2Ratio", names(.))),
    df %>% dplyr::select(-grep("pVal|adjP|log2Ratio", names(.)), -id)
  ) %>% 
    dplyr::filter(!!rlang::sym(id) %in% gns)

  suppressWarnings(df[is.na(df)] <- "") # Cytoscape compatibility

  write.table(df, file.path(dat_dir, "Protein\\String", outnm_expr), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
}


#' String analysis
#'
#' The input contains pVal fields
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnCorr_logFC
#' @inheritParams anal_prnString
#' @import dplyr purrr rlang fs
#' @importFrom magrittr %>%
stringTest <- function(df = NULL, id = gene, col_group = Group, col_order = Order, 
                       label_scheme_sub = NULL, 
                       db_path = "~\\proteoQ\\dbs\\string", score_cutoff = .7, 
                       scale_log2r = TRUE, complete_cases = FALSE, 
                       filepath = NULL, filename = NULL, ...) {
  
  # `scale_log2r` not used; both `_N` and `_Z` columns will be kept
  
  stopifnot(nrow(label_scheme_sub) > 0)
  stopifnot(rlang::is_double(score_cutoff))
  
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  if (score_cutoff <= 1) score_cutoff <- score_cutoff * 1000
  
  dir.create(file.path(dat_dir, "Protein\\String\\cache"), recursive = TRUE, showWarnings = FALSE)
  
  if (!fs::dir_exists(db_path)) {
    new_db_path <- fs::path_expand_r(db_path)
    new_db_path2 <- fs::path_expand(db_path)
    
    if (fs::dir_exists(new_db_path)) {
      db_path <- new_db_path
    } else if (fs::dir_exists(new_db_path2)) {
      db_path <- new_db_path2
    } else {
      stop(db_path, " not existed.", call. = FALSE)
    }
    
    rm(new_db_path, new_db_path2)
  }
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  
  stopifnot(id %in% names(df))
  stopifnot("species" %in% names(df))
  
  species <- unique(df$species) %>% 
    .[!is.na(.)] %>% 
    as.character()
  
  stopifnot(length(species) >= 1)
  
  run_scripts <- FALSE
  if (run_scripts) {
    if (scale_log2r) {
      df <- df %>% 
        dplyr::select(-grep("N_log2_R", names(.)))
    } else {
      df <- df %>% 
        dplyr::select(-grep("Z_log2_R", names(.)))
    }    
  }

  df <- df %>% rm_pval_whitespace()

  purrr::walk(species, annot_stringdb, df, db_path, id, score_cutoff, filepath, filename)
}


#'STRING outputs of protein-protein interactions
#'
#'\code{anal_prnString} prepares the data of both
#'\href{https://string-db.org/}{STRING} protein-protein interactions
#'(ppi) and companion protein expressions.
#'
#'Convenience features, such as data row filtration via \code{filter_} varargs,
#'are available. The ppi file, \code{..._ppi.tsv}, and the expression file,
#'\code{..._expr.tsv}, are also compatible with a third-party
#'\href{https://cytoscape.org/}{Cytoscape}.
#'
#'@inheritParams anal_prnTrend
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~\\proteoQ\\dbs\\string"}.
#'@param score_cutoff Numeric; the threshold in the \code{combined_score} of
#'  protein-protein interaction. The default is 0.7.
#'@param scale_log2r Not currently used. Values before and after scaling will be
#'  both reported.
#'@param filename Use system default. Otherwise, the basename will be prepended
#'  to \code{_[species]_ppi.tsv} for network data and \code{_[species]_expr.tsv}
#'  for expression data.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration
#'  against data in a primary file of \code{\\Model\\Protein[_impNA]_pVals.txt}.
#'  See also \code{\link{normPSM}} for the format of \code{filter_} statements.
#'  \cr \cr \code{arrange_}: Variable argument statements for the row ordering
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{prnHM}} for the format of \code{arrange_} statements.
#'@seealso \code{\link{dl_stringdbs}} for database downloads. \cr
#'  \code{\link{prnSig}} for significance tests \cr
#'@example inst/extdata/examples/getStringDB_.R
#'
#'@import dplyr purrr rlang fs
#'@export
anal_prnString <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                            db_path = "~\\proteoQ\\dbs\\string", score_cutoff = .7, 
                            df = NULL, filepath = NULL, filename = NULL, ...) {
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_pepString"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_prnString"))
    }
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = NULL, col_group = NULL, col_order = NULL,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "String")(db_path = db_path, score_cutoff = score_cutoff, adjP = adjP, 
                                  ...)
}




#' Trend analysis
#' 
#' @inheritParams anal_prnTrend
#' @inheritParams info_anal
#' @inheritParams gspaTest
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
#' @inheritParams plot_prnTrend
#' @inheritParams info_anal
#' @inheritParams gspaTest
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
#'  \code{arrange_}: Variable argument statements for the row ordering against
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
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
#'  filtration against data in secondary file(s) of
#'  \code{[...]Protein_Trend_[...].txt}. See also \code{\link{prnGSPAHM}} for
#'  the format of \code{filter2_} statements. \cr \cr \code{arrange2_}: Variable
#'  argument statements for the row ordering against data in secondary file(s)
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
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

#'Prepares data for analysis
#'
#'\code{prepDM} prepares a minimal adequate data frame for subsequent analysis.
#'
#'@inheritParams  prnEucDist
#'@inheritParams  info_anal
#'@param sub_grp Numeric.  A list of sample IDs that will be used in subsequent
#'  analysis.
#'@param type The type of data, for example ratio or intensity.
#'@return A data frame tailored for subsequent analysis.
#'
#' @examples
#' \donttest{tempData <- prepDM(df, entrez, scale_log2r, label_scheme_sub$Sample_ID)}
#'@import dplyr
#'@importFrom magrittr %>%
prepDM <- function(df, id, scale_log2r, sub_grp, type = "ratio", anal_type) {
  stopifnot(nrow(df) > 0)
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  id <- rlang::as_string(rlang::enexpr(id))
  if (anal_type %in% c("ESGAGE", "GSVA")) id <- "entrez"
  if ((anal_type %in% c("GSEA")) & (id != "gene")) stop("Primary ID is not `gene`.")
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  pattern <- "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("

  df <- df %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
    reorderCols(endColIndex = grep(pattern, names(.)), col_to_rn = id)

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


#' Prefix form of colnames(x)[c(2, 5, ...)] for use in pipes
#'
#' \code{names_pos<-} rename the columns at the indexes of \code{pos}.
#'
#' @param x A data frame.
#' @param pos Numeric.  The index of columns for name change.
#' @param value Characters.  The new column names.
#' @return The data frame with new names.
#'
#' @import dplyr
#' @importFrom magrittr %>%
`names_pos<-` <- function(x, pos, value) {
	names(x)[pos] <- value
	x
}


#' Re-order file names
#'
#' \code{reorder_files} re-orders file names by TMT set numbers then by LCMS
#' injection numbers.
#'
#' @param filelist A list of file names.
#' @param n_TMT_sets the number of multiplex TMT experiments.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
reorder_files <- function(filelist, n_TMT_sets) {
  newlist <- NULL
  for (i in seq_len(n_TMT_sets))
    newlist <- c(newlist, filelist[grep(paste0("set*.",i), filelist, ignore.case = TRUE)])
  return(newlist)
}


#' Re-order columns in a data frame
#'
#' \code{reorderCols} re-orders columns in a data frame.
#'
#' @param df A data frame.
#' @param endColIndex the indexes of columns to be moved to the end of
#'   \code{df}.
#' @param col_to_rn the column identifier where the values under that column
#'   will be used as row names.
#' @import dplyr rlang
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
reorderCols <- function (df, endColIndex, col_to_rn) {
	if (length(endColIndex) == 0) endColIndex <- grep("I[0-9]{3}|log2_R[0-9]{3}", names(df))
	
	df <- cbind(df[, -endColIndex], df[, endColIndex])

	if(sum(duplicated(df[[col_to_rn]])) > 0) {
		df <- df %>%
				dplyr::group_by(!!ensym(col_to_rn)) %>%
				dplyr::mutate(Index = row_number()) %>%
				dplyr::mutate(indexed_names = paste(get(col_to_rn), Index, sep=".")) %>%
				data.frame(check.names = FALSE) %>%
				`rownames<-`(.[, "indexed_names"]) %>%
				dplyr::select(-c("Index", "indexed_names"))
	} else {
		rownames(df) <- df[[col_to_rn]]
	}

	return(df)
}


#' Re-order columns in a data frame
#'
#' \code{reorderCols2} re-orders columns in a data frame.
#'
#' @param df A data frame.
#' @param pattern columns matched the pattern will be moved to the end of
#'   \code{df}.
#' @import dplyr rlang
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
reorderCols2 <- function (df = NULL, pattern = NULL) {
  if (is.null(df)) stop("`df` cannot be `NULL`.", call. = FALSE)
  
  if (is.null(pattern)) 
    pattern <- "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("
  
  endColIndex <- grep(pattern, names(df))
  
  if (length(endColIndex) > 0) 
    df <- dplyr::bind_cols(df[, -endColIndex], df[, endColIndex])
  
  return(df)
}


#' Replace zero intensity with NA
#'
#' \code{na_zeroIntensity} replaces zero intensity with NA to avoid -Inf in
#' log10 transformation.
#'
#' @param df A list of file names.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
na_zeroIntensity <- function (df) {
	ind <- grep("I1[0-9]{2}", names(df))

	temp <- df[, ind]
	temp[temp == 0] <- NA
	df[, ind] <- temp

	return(df)
}

#' Summarizes numeric values
#'
#' \code{aggrNums} summarizes \code{log2FC} and \code{intensity} by the
#' descriptive statistics of \code{c("mean", "median", "weighted.mean",
#' "top.3")}
#'
#' @param f A function for data summarization.
#' @examples \donttest{df_num <- aggrNums(median)(df, prot_acc, na.rm = TRUE)}
#' @import dplyr rlang
#' @importFrom magrittr %>%
aggrNums <- function(f) {
	function (df, id, ...) {
		id <- rlang::as_string(rlang::enexpr(id))
		dots <- rlang::enexprs(...)

		df %>%
			dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise_all(~ f(., !!!dots))
	}
}


#' Calculates weighted mean
#'
#' \code{TMT_wt_mean} calculates the weighted mean of \code{log2FC} and
#' \code{intensity}.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param id The variable to summarize \code{log2FC}.
#' @param ... Additional arguments for \code{weighted.mean}.
#' @import dplyr rlang
#' @importFrom stringr str_length
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
TMT_wt_mean <- function (x, id, ...) {
	id <- rlang::as_string(rlang::enexpr(id))
	dots <- rlang::enexprs(...)

	load(file = file.path(dat_dir, "label_scheme.rda"))
	TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()

	x_R <- x %>%
			dplyr::select(id, grep("^log2_R[0-9]{3}", names(.))) %>%
			tidyr::gather(key = variable, value = log2_R, -id)

	x_I <- x %>%
			dplyr::select(id, grep("^I[0-9]{3}", names(.))) %>%
			tidyr::gather(key = variable, value = Intensity, -id)

	x_wt <- cbind.data.frame(x_I, log2_R = x_R[, c("log2_R")]) %>%
			dplyr::mutate(variable = gsub("^I", "log2_R", variable)) %>%
			dplyr::mutate(TMT = gsub("^log2_R(.*)\\s+.*", "TMT-\\1", variable)) %>%
			dplyr::mutate(TMT = factor(TMT, levels = TMT_levels)) %>%
			dplyr::mutate(variable = factor(variable, levels = unique(variable))) %>%
			dplyr::mutate(n = stringr::str_length(.[[id]])) %>%
			dplyr::filter(n != 0) %>%
			dplyr::group_by(!!rlang::sym(id), variable) %>%
			dplyr::summarise(log2_R = weighted.mean(log2_R, log10(Intensity), !!!dots)) %>%
			tidyr::spread(variable, log2_R)

	x <- x %>%
			dplyr::select(id, grep("(N_log2_R)|(Z_log2_R)|(I)[0-9]{3}", names(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise_all(funs(mean(., !!!dots))) %>%
			dplyr::left_join(x_wt, by = id)

	rm(x_R, x_I, x_wt)

	cbind.data.frame(
		x[, names(x) == id],
		x[, grepl("^I[0-9]{3}", names(x))],
		x[, grepl("^N_I[0-9]{3}", names(x))],
		x[, grepl("^log2_R[0-9]{3}", names(x))],
		x[, grepl("^N_log2_R[0-9]{3}", names(x))],
		x[, grepl("^Z_log2_R[0-9]{3}", names(x))])
}


#' Calculate top.3
#'
#' \code{TMT_wt_mean} calculates the weighted mean of \code{log2FC} and
#' \code{intensity}.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param id The variable to summarize \code{log2FC}.
#' @param ... Additional arguments for \code{mean}.
#' @examples \donttest{df_num <- TMT_top_n(df, prot_acc, na.rm = TRUE)}
#' @import dplyr rlang
#' @importFrom magrittr %>%
TMT_top_n <- function (x, id, ...) {

	id_var <- rlang::enexpr(id)
	id_nm <- rlang::expr_name(id_var)

	dots <- rlang::enexprs(...)
	args <- lapply(dots, expr_name)

	x %>%
		dplyr::select(id_nm, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
		dplyr::mutate(sum_Intensity =rowSums(.[grep("^I1[0-9]{2}", names(.))], na.rm = TRUE)) %>%
		dplyr::group_by(!!id_var) %>%
		dplyr::top_n(n = 3, wt = sum_Intensity) %>%
		dplyr::select(-sum_Intensity) %>%
		dplyr::summarise_all(funs(median(., !!!args)))
}


#' Finds all-zero column(s)
#'
#' \code{not_all_zero} identifies the column indexes with all NA values.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
not_all_zero <- function (x) (colSums(x != 0, na.rm = TRUE) > 0)


#' Finds all-NA column(s)
#'
#' \code{not_all_NA} identifies the column indexes with all NA values.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @import dplyr rlang
#' @importFrom magrittr %>%
not_all_NA <- function (x) (colSums(!is.na(x), na.rm = TRUE) > 0)


#' Finds all-NaN column(s)
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param ... The same in \code{sum}.
not_all_nan <- function(x, ...) {
  sum(is.nan(x), ...) != length(x)
}


#' Finds all-NaN column(s)
#' 
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param ... The same in \code{sum}.
is_all_nan <- function(x, ...) {
  sum(is.nan(x), ...) == length(x)
}


#' Sets up the column annotation in heat maps
#' 
#' @inheritParams prnEucDist
#' @inheritParams gspa_colAnnot
#' @import dplyr rlang
#' @importFrom magrittr %>%
colAnnot <- function (annot_cols = NULL, sample_ids = NULL, annot_colnames = NULL) {
	if (is.null(annot_cols)) return(NA)

  load(file = file.path(dat_dir, "label_scheme.rda"))
	exists <- annot_cols %in% names(label_scheme)

	if (sum(!exists) > 0) {
		warning(paste0("Column '", annot_cols[!exists], "'",
		               " not found in \'label_scheme\' and will be skipped."))
		annot_cols <- annot_cols[exists]
	}

	if (length(annot_cols) == 0) return(NA)

	x <- label_scheme %>%
		dplyr::filter(Sample_ID %in% sample_ids) %>%
		dplyr::select(annot_cols, Sample_ID) 
	
	idx <- !not_all_NA(x)
	if (sum(idx) > 0) stop("No aesthetics defined under column `", names(x)[idx], "`")

	x <- x %>%
		dplyr::filter(!grepl("^Empty\\.", .[["Sample_ID"]]),
		              !is.na(.[["Sample_ID"]])) %>%
		data.frame(check.names = FALSE) %>%
		`rownames<-`(.[["Sample_ID"]])

	if (any(duplicated(x[["Sample_ID"]]))) stop("Duplicated sample IDs found\n")

	if (!"Sample_ID" %in% annot_cols) x <- x %>% dplyr::select(-Sample_ID)

	if (ncol(x) == 0) return(NA)

	if ("TMT_Set" %in% names(x)) {
		x <- x %>%
			tibble::rownames_to_column() %>%
			mutate(TMT_Set, as.factor(TMT_Set)) %>%
			tibble::column_to_rownames(var = "rowname")
	}

	return(x)
}


#' Customizes the colors in column annotation
#' 
#' @param annotation_col The same as in \link[pheatmap]{pheatmap}.
#' @import dplyr rlang
#' @importFrom magrittr %>%
setHMColor <- function (annotation_col) {
	ncol <- ncol(annotation_col)

	if(is.null(ncol)) return (NA)
	if(ncol == 0) return (NA)

	suppressWarnings(
	  annotation_col <- annotation_col %>%
	    dplyr::mutate_at(vars(one_of("TMT_Set")), ~ factor(TMT_Set)) %>%
	    dplyr::mutate_if(is.character, as.factor)  
	)

	palette <- lapply(names(annotation_col), function(x) {
		n <- nlevels(annotation_col[[x]])

		palette <- if(n <= 9 & n >= 3) {
		  brewer.pal(n, name = "Set2")
		} else if(n > 9) {
			colorRampPalette(brewer.pal(n = 7, "Set1"))(n)
		} else if(n == 2) {
		  c("#66C2A5", "#FC8D62")
		} else if(n == 1) {
			# c("#E41A1C")
		  c("#66C2A5")
		} else if(n == 0) {
			colorRampPalette(brewer.pal(n = 9, "YlOrBr"))(100)
		}

		names(palette) <- levels(annotation_col[[x]])

		return(palette)
	})

	names(palette) <- names(annotation_col)

	return(palette)
}


#' Sets the upper and the lower limits in the range of heat maps
#'
#' \code{setHMlims} imputes values beyond the limits to the corresponding
#' limits.
#'
#' @param x A data frame of \code{log2FC}.
#' @param xmin the lower limit.
#' @param xmax the upper limit.
#' @import dplyr rlang
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
setHMlims <- function (x, xmin, xmax) {
	x[x < xmin] <- xmin
	x[x > xmax] <- xmax

	x
}


#' Calculate the log2-ratio to the control group of samples
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @param nm_ctrl The names of samples that belong to the control group.
#' @examples \donttest{ratio_toCtrl(df, "gene", label_scheme_sub, Heatmap_Group)}
ratio_toCtrl <- function(df, id, label_scheme_sub, nm_ctrl) {
	id <- rlang::as_string(rlang::enexpr(id))

	nm_ctrl <- rlang::as_string(rlang::ensym(nm_ctrl))

	x <- df %>%
		dplyr::select(which(names(.) %in% label_scheme_sub$Sample_ID)) %>%
		`rownames<-`(df[[id]])

	col_ratio <- which(names(x) %in% label_scheme_sub$Sample_ID)

	col_ctrl <- which(names(x) %in%
	                    label_scheme_sub[!is.na(label_scheme_sub[[nm_ctrl]]), "Sample_ID"])

	x[, col_ratio] <- sweep(x[, col_ratio], 1,
	                        rowMeans(x[, col_ctrl, drop = FALSE], na.rm = TRUE), "-")

	x <- x %>% tibble::rownames_to_column(id)

	df %>%
		dplyr::select(which(!names(df) %in% label_scheme_sub$Sample_ID)) %>%
		dplyr::left_join(x, by = id)	%>%
		`rownames<-`(.[[id]])
}


#'Imputation of NA values
#'
#'
#'\code{imputeNA} imputes NA \code{log2FC} values in protein or peptide data.
#'Users should avoid call the method directly, but instead use the following
#'wrappers.
#'
#'@param overwrite Logical. If true, overwrite the previous results. The default
#'  is FALSE.
#'@inheritParams prnHist
#'@inheritParams info_anal
#'@param ... Parameters for \code{\link[mice]{mice}}
#'@return \code{Peptide_impNA.txt} for peptide data and \code{Protein_impNA.txt}
#'  for protein data.
#'
#'@import dplyr purrr rlang mice
#'@importFrom magrittr %>%
#'@export
imputeNA <- function (id, overwrite = FALSE, ...) {
	my_mice <- function (data, ...) {
		data <- rlang::enexpr(data)
		dots <- rlang::enexprs(...)

		mice_call <- rlang::expr(mice(data = !!data, !!!dots))
		rlang::eval_bare(mice_call, env = caller_env())
	}

	handleNA <- function (x, ...) {
		ind <- purrr::map(x, is.numeric) %>% unlist()

		if (sum(ind) < ncol(x)) cat(names(ind)[!ind], "skipped.")

		# handle special column names
		nm_orig <- names(x[, ind])

		x[, ind] <- x[, ind] %>%
			data.frame(check.names = TRUE) %>%
			my_mice(...) %>%
			mice::complete(1) %>%
			`colnames<-`(nm_orig)

		return(x)
	}

	id <- rlang::as_string(rlang::enexpr(id))

	if (id %in% c("pep_seq", "pep_seq_mod")) {
		src_path <- file.path(dat_dir, "Peptide", "Peptide.txt")
		filename <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
	} else if (id %in% c("prot_acc", "gene")) {
		src_path <- file.path(dat_dir, "Protein", "Protein.txt")
		filename <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
	}

	if ((!file.exists(filename)) || overwrite) {
		df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
		                        comment.char = "#"), error = function(e) NA)

		if (!is.null(dim(df))) {
			message(paste("File loaded:", gsub("\\\\", "/", src_path)))
		} else {
			stop(paste("File or directory not found:", gsub("\\\\", "/", src_path)))
		}

		df[, grep("N_log2_R", names(df))]  <- df[, grep("N_log2_R", names(df))]  %>% handleNA(...)

		fn_params <- file.path(dat_dir, "Protein\\Histogram", "MGKernel_params_N.txt")
		if (file.exists(fn_params)) {
			cf_SD <- 
				read.csv(fn_params, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>%
				dplyr::filter(!duplicated(Sample_ID)) %>%
				dplyr::select(Sample_ID, fct)

			df[, grep("Z_log2_R", names(df))] <-
				mapply(normSD,
				       df[, grepl("^N_log2_R[0-9]{3}", names(df))], center = 0, SD = cf_SD$fct, SIMPLIFY = FALSE) %>%
				data.frame(check.names = FALSE) %>%
				`colnames<-`(gsub("N_log2", "Z_log2", names(.)))
		} else {
			df[, grep("Z_log2_R", names(df))]  <- df[, grep("Z_log2_R", names(df))] %>% handleNA(...)
		}

		if (any(duplicated(df[[id]]))) {
			if (id == "pep_seq") {
				warning("\`pep_seq\` is not uqique for rownames; use \`pep_seq_mod\` instead.\n", call. = FALSE)
				rownames(df) <- df[["pep_seq_mod"]]
			}
			if (id == "gene") {
				warning("\`gene\` is not uqique for rownames; use \`prot_acc\` instead.\n", call. = FALSE)
				rownames(df) <- df[["prot_acc"]]
			}
		} else {
			rownames(df) <- df[[id]]
		}

		write.table(df, filename, sep = "\t", col.names = TRUE, row.names = FALSE)
	} else {
	  warning("NA imputation has been previously performed!\n", call. = FALSE)
	}
}


#'Impute NA values for peptide data
#'
#'\code{pepImp} is a wrapper of \code{\link{imputeNA}} for peptide data with
#'auto-determination of \code{id}.
#'
#'@rdname imputeNA
#'
#' @examples 
#' \donttest{
#' # ===================================
#' # Peptide NA imputation
#' # ===================================
#' pepImp(
#'   m = 3,
#'   maxit = 3,
#' )
#' }
#' 
#'@export
pepImp <- function (...) {
  check_dots(c("id"), ...)
  id <- match_call_arg(normPSM, group_psm_by)
  imputeNA(id = !!id, ...)
}


#'Impute NA values for protein data
#'
#'\code{prnImp} is a wrapper of \code{\link{imputeNA}} for protein data with
#'auto-determination of \code{id}.
#'
#'@rdname imputeNA
#'
#' @examples
#' \donttest{
#' # ===================================
#' # Protein NA imputation
#' # ===================================
#' prnImp(
#'   m = 5,
#'   maxit = 5,
#' )
#' }
#'
#'@export
prnImp <- function (...) {
  check_dots(c("id"), ...)
  id <- match_call_arg(normPSM, group_pep_by)
  imputeNA(id = !!id, ...)
}


#' Species lookup
#' @inheritParams load_dbs
sp_lookup <- function(species) {
  switch(species, 
         human = "hs",
         mouse = "mm",
         rat = "rn",
         unknown = "unknown"
  )    
}


#' Taxonomy lookup
#' @inheritParams load_dbs
taxid_lookup <- function(species) {
  switch (species,
          human = 9606, 
          mouse = 10090,
          rat = 10116, 
          unknown = 999999
  )
}


#' Reversed taxonomy lookup
#' @inheritParams load_dbs
taxid_lookup_rev <- function(species) {
  switch (species,
          "9606" = "human", 
          "10090" = "mouse",
          "10116" = "rat", 
          "999999" = "unknown"
  )
}


#' Species lookup UpperLower (Ul)
#' @inheritParams load_dbs
sp_lookup_Ul <- function(species) {
  switch(species, 
         human = "Hs",
         mouse = "Mm",
         rat = "Rn",
         unknown = "Unknown"
  )    
}


#' Add Entrez IDs
#' @param acc_lookup A data frame of protein accession lookups
add_entrez <- function (acc_lookup) {
  sp_map <- c(
    human = "hs",
    mouse = "mm",
    rat = "rn"
  )
  
  acc_map <- c(
    uniprot_acc = "uniprot_",
    uniprot_id = "uniprot_", 
    refseq_acc = "refseq_"
  )
  
  acc_type <- unique(acc_lookup$acc_type)
  species <- unique(acc_lookup$species)
  abbr_sp <- sp_map[species]
  abbr_acc <- acc_map[acc_type]
  
  entrez_key <- switch(acc_type, 
    uniprot_acc = "uniprot_acc", 
    uniprot_id = "uniprot_acc", 
    refseq_acc = "refseq_acc", 
    stop("Unknown `accession type`.", Call. = FALSE)
  )

  stopifnot(length(acc_type) == 1)
  
  if (! all(species %in% c("human", "mouse", "rat"))) {
    warning("No default `entrez` lookups available for species other than `human`, `mouse` and `rat`.", 
         "\nTo annotate, provide the file path and name(s) via argument `entrez`.", 
         "\nSee also `?Uni2Entrez` and `?Ref2Entrez` for preparing custom `entrez` lookups.", 
         call. = FALSE)
  }

  # multiple uniprot_acc(s) can share the same entrez id
  if (all(species %in% c("human", "mouse", "rat"))) {
    filelist <- paste0(abbr_acc, "entrez_", abbr_sp)
    data(package = "proteoQ", list = filelist)
    entrez <- purrr::map(filelist, ~ get(.x)) %>% 
      dplyr::bind_rows() %>% 
      dplyr::filter(!duplicated(.[[entrez_key]])) %>% 
      dplyr::mutate(!!entrez_key := as.character(!!rlang::sym(entrez_key)))
  
    if (acc_type == "uniprot_id") {
      acc_lookup <- left_join(acc_lookup, entrez, by = c("uniprot_acc" = entrez_key))
    } else {
      acc_lookup <- left_join(acc_lookup, entrez, by = c("prot_acc" = entrez_key))
    }    
  } else {
    acc_lookup <- acc_lookup %>% dplyr::mutate(entrez = NA)
  }

  invisible(acc_lookup)
}


#' Add custom Entrez IDs
#' 
#' @inheritParams  add_entrez
#' @inheritParams normPSM
#' @inheritParams annotKin
add_custom_entrez <- function(acc_lookup, entrez, acc_type) {
  if (all(file.exists(entrez))) {
    entrez <- purrr::map(entrez, ~ readRDS(.x)) %>% do.call(rbind, .)
    
    if ("species" %in% names(acc_lookup) && "species" %in% names(entrez)) {
      acc_lookup <- acc_lookup %>% dplyr::select(-species)
    }
    
    stopifnot("prot_acc" %in% names(acc_lookup))
    
    if (acc_type == "uniprot_acc") {
      if (!"uniprot_acc" %in% names(entrez)) {
        stop("Argument `entrez` does not link to UniProt database(s).", call. = FALSE)
      }
      acc_lookup <- acc_lookup %>% dplyr::left_join(entrez, by = c("prot_acc" = "uniprot_acc"))
    } else if ("uniprot_acc" %in% names(acc_lookup)) {
      if (!"uniprot_acc" %in% names(entrez)) {
        stop("Argument `entrez` does not link to UniProt database(s).", call. = FALSE)
      }
      acc_lookup <- acc_lookup %>% dplyr::left_join(entrez, by = "uniprot_acc")
    } else if (acc_type == "refseq_acc") {
      if (!"refseq_acc" %in% names(entrez)) {
        stop("Argument `entrez` does not link to RefSeq database(s).", call. = FALSE)
      }
      acc_lookup <- acc_lookup %>% dplyr::left_join(entrez, by = c("prot_acc" = "refseq_acc"))
    } else {
      stop("Unknown protein accession types.", call. = FALSE)
    }
    
    if (all(is.na(acc_lookup$entrez))) {
      warning("No matched UniProt accessions between PSM data and `entrez` database(s).", 
              "\nSpecies in the `entrez` file(s) are probably incorrect.", 
              call. = FALSE)
    } else if ((sum(is.na(acc_lookup$entrez))/nrow(acc_lookup)) > .8) {
      warning("Over 80% UniProt accessions from PSM data do not have corresponding `entrez` IDs.", 
              "\nSpecies in the `entrez` file(s) are probably incorrect or incomplete.", 
              call. = FALSE)
    }
    
  } else {
    stop("Wrong `entrez` file path(s) or name(s).", 
         "\nMake sure that the file type is `.rds`, not `.rda`", 
         call. = FALSE)
  }
  
  invisible(acc_lookup)
} 


#' Determine the protein accession type from PSM tables
#'
#' Find the protein accession from a non-cRAP entry and parse it.
#' 
#'@param df An input data frame.
parse_acc <- function(df) {
  stopifnot("prot_acc" %in% names(df))
  
  prn_acc <- df %>%
    dplyr::filter(!grepl("^REV__", prot_acc)) %>% 
    dplyr::filter(!grepl("^CON__", prot_acc)) %>% 
    dplyr::select(prot_acc) %>%
    unlist %>%
    .[1]

  if (grepl("[[:alnum:]]+_[A-Z]{1,5}$", prn_acc)) {
    acc_type <- "uniprot_id"
  } else if (grepl("^[XN]{1}[MRP]{1}_", prn_acc)) {
    acc_type <- "refseq_acc"
  } else if (grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", prn_acc)) {
    acc_type <- "uniprot_acc"
  } else {
    acc_type <- "unknown"
  }
  
  return(acc_type)
}


#' Parse UniProt FASTA for accession lookups
#'
#' @param df An input data frame.
#' @inheritParams add_entrez
#' @inheritParams normPSM
parse_uniprot_fasta <- function (df, fasta, entrez) {
  
  na_genes_by_acc <- function(acc_lookup) {
    stopifnot("prot_acc" %in% names(acc_lookup), "gene" %in% names(acc_lookup))
    
    na_gene <- (is.na(acc_lookup$gene)) | (stringr::str_length(acc_lookup$gene) == 0)
    acc_lookup$gene <- as.character(acc_lookup$gene)
    acc_lookup$gene[na_gene] <- as.character(acc_lookup$prot_acc[na_gene])
    
    return(acc_lookup)
  }
  
  na_species_by_org <- function(acc_lookup) {
    stopifnot("organism" %in% names(acc_lookup), "species" %in% names(acc_lookup))
    
    na_species <- (is.na(acc_lookup$species)) | (stringr::str_length(acc_lookup$species) == 0)
    acc_lookup$species <- as.character(acc_lookup$species)
    acc_lookup$species[na_species] <- as.character(acc_lookup$organism[na_species])
    
    return(acc_lookup)
  }
  
  
  my_lookup <- c(
    "Homo sapiens" = "human",
    "Mus musculus" = "mouse",
    "Rattus norvegicus" = "rat"
  )

  stopifnot ("prot_acc" %in% names(df))

  acc_type <- parse_acc(df)
  
  if (acc_type == "unknown") {
    warning("Unknown accession; use hypothetical `refseq_acc`.")
    acc_type <- "refseq_acc"
  } else if (! acc_type %in% c("refseq_acc", "uniprot_id", "uniprot_acc")) {
    stop("The type of protein accesion needs to one of \'uniprot_id\', \'uniprot_acc\' or \'refseq_acc\'",
         call. = FALSE)
  }
  
  if (is.null(fasta)) stop("FASTA file(s) are required.", call. = FALSE)
  if (!all(file.exists(fasta))) stop("Wrong FASTA file path(s) or name(s).", call. = FALSE)
  
  fasta <- purrr::map(fasta, ~ {
    seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
  }) %>% do.call(`c`, .) 
  
  if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
    
    # check cRAP...
    
    acc_lookup <- names(fasta) %>% 
      gsub("^..\\|", "", .) %>% 
      stringr::str_split("\\|", simplify = TRUE) %>% 
      data.frame() %>% 
      `colnames<-`(c("uniprot_acc", "uniprot_id")) %>% 
      dplyr::bind_cols(data.frame(fasta_name = names(fasta)), .) %>% 
      dplyr::filter(.[[acc_type]] %in% unique(df$prot_acc)) %>% 
      dplyr::filter(!duplicated(.[[acc_type]]))
  } else if (acc_type == "refseq_acc") {
    acc_lookup <- tibble::tibble(fasta_name = names(fasta), refseq_acc = fasta_name) %>% 
      dplyr::filter(.[[acc_type]] %in% unique(df$prot_acc)) %>% 
      dplyr::filter(!duplicated(.[[acc_type]]))
  }
  
  fasta <- fasta %>% 
    .[names(.) %in% unique(acc_lookup$fasta_name)] %>% 
    .[! duplicated(names(.))]
  
  if (length(fasta) == 0) {
    stop("No fasta entries match protein accessions; probably wrong fasta file.", 
         call. = FALSE)
  } else {
    write.fasta(sequences = fasta, names = seqinr::getAnnot(fasta), nbchar = 80, 
                file.out = file.path(dat_dir, "my_project.fasta"))
  }
  
  if (acc_type == "uniprot_acc") {
    names(fasta) <- gsub("^..\\|(.*)\\|.*$", "\\1", names(fasta))
  } else if (acc_type == "uniprot_id") {
    names(fasta) <- gsub("^.*\\|(.*)$", "\\1", names(fasta))
  } else if (acc_type == "refseq_acc") {
    names(fasta) <- gsub("\\.[^\\.]*$", "", names(fasta))
  }
  
  acc_lookup <- local({
    fasta_smry <- dplyr::bind_cols(
      prot_acc = names(fasta), 
      prot_desc = seqinr::getAnnot(fasta) %>% 
        purrr::map(., `[[`, 1) %>% 
        unlist(), 
      prot_mass = purrr::map_dbl(fasta, ~ {seqinr::getSequence(.x) %>% seqinr::pmw()}), 
      prot_len = getLength(fasta)
    ) %>% 
      dplyr::filter(!duplicated(.$prot_acc)) %>% 
      dplyr::mutate(acc_type = acc_type) %>% 
      dplyr::mutate(prot_mass = round(prot_mass, digits = 0))
    
    if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
      fasta_smry <- fasta_smry %>% 
        dplyr::mutate(prot_desc = gsub("^.*\\|.*\\s+?", "", prot_desc))
    } else if (acc_type == "refseq_acc") {
      fasta_smry <- fasta_smry %>% 
        dplyr::mutate(prot_desc = gsub("^.*_.*\\s+?", "", prot_desc))
    }
    
    acc_lookup <- acc_lookup %>% 
      dplyr::mutate(!!acc_type := as.character(!!rlang::sym(acc_type)))
    
    if (acc_type != "uniprot_acc" && "uniprot_acc" %in% names(acc_lookup)) {
      acc_lookup <- acc_lookup %>% 
        dplyr::mutate(uniprot_acc = as.character(uniprot_acc))
    }
    
    acc_lookup <- fasta_smry %>% 
      dplyr::left_join(acc_lookup, by = c("prot_acc" = acc_type))
  })
  
  if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
    genes <- local({
      genes <- acc_lookup %>% dplyr::select(prot_acc, prot_desc)
      
      na_genes <- genes %>% 
        dplyr::filter(!grepl("GN=", .$prot_desc)) %>% 
        dplyr::mutate(gene = NA)
      
      genes <- genes %>% 
        dplyr::filter(grepl("GN=", .$prot_desc)) %>% 
        dplyr::mutate(gene = gsub("^.*GN=(\\S+)\\s*.*", "\\1", prot_desc)) %>% 
        dplyr::bind_rows(., na_genes) %>% 
        na_genes_by_acc()
    })
    
    acc_lookup <- local({
      na_org <- genes %>% 
        dplyr::filter(!grepl("OS=", .$prot_desc)) %>% 
        dplyr::mutate(organism = NA)
      
      acc_lookup <- genes %>% 
        dplyr::filter(grepl("OS=", .$prot_desc)) %>% 
        dplyr::mutate(organism = gsub("^.*OS=(.*?)=.*$", "\\1", prot_desc)) %>% 
        dplyr::mutate(organism = gsub("\\s\\S*$", "", organism)) %>% 
        dplyr::bind_rows(., na_org) %>% 
        dplyr::select(-prot_desc) %>% 
        dplyr::right_join(acc_lookup, by = "prot_acc")
      
      acc_lookup <- acc_lookup %>% 
        dplyr::mutate(species = my_lookup[.$organism]) %>% 
        na_species_by_org()
    })
    
    if (is.null(entrez)) {
      acc_lookup <- add_entrez(acc_lookup)
    } else {
      acc_lookup <- add_custom_entrez(acc_lookup, entrez, acc_type)
    }
  } else if (acc_type == "refseq_acc") {
    warning("`refseq_acc` only available for 'human' and/or 'mouse'.", call. = FALSE)
    
    ## df$prot_acc can contain cRAPs
    
    acc_lookup <- acc_lookup %>% 
      dplyr::filter(grepl("^[XN]{1}[MRP]{1}_", prot_acc)) %>% 
      # dplyr::filter(organism %in% uniprot_species$organism) %>% 
      # dplyr::filter(!grepl("B99907|TRYP_PromTArt7|P00974|BPT1_BOVIN", prot_acc)) %>% 
      # dplyr::mutate(prot_acc = gsub("^..\\|(.*)\\|.*", "\\1", prot_acc)) %>% 
      # dplyr::filter(! prot_acc %in% prn_annot_crap$refseq_acc) %>% 
      # dplyr::filter(! prot_acc %in% prn_annot_crap$uniprot_acc) %>% 
      dplyr::mutate(organism = gsub("^.*\\[(.*)\\]\\.*.*", "\\1", prot_desc)) %>% 
      dplyr::mutate(species = my_lookup[.$organism]) %>% 
      na_species_by_org()

    if (is.null(entrez)) {
      acc_lookup <- add_entrez(acc_lookup)
    } else {
      # check the gene column is present...
      acc_lookup <- add_custom_entrez(acc_lookup, entrez, acc_type)
    }
    
    acc_lookup <- acc_lookup %>% na_genes_by_acc()
  }
  
  save(acc_lookup, file = file.path(dat_dir, "acc_lookup.rda"))
  
  return(acc_lookup)
}


#' Adds protein annotation
#'
#' \code{annotPrn} cross-referencing proteins among \code{uniprot_acc},
#' \code{uniprot_id}, \code{refseq} and \code{entrez}.
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @import plyr dplyr purrr rlang seqinr stringr 
#' @importFrom magrittr %>% %$% %T>% 
annotPrn <- function (df, fasta, entrez) {
	acc_lookup <- parse_uniprot_fasta(df, fasta, entrez)
	
	acc_lookup <- dplyr::bind_cols(
	  acc_lookup %>% 
	    dplyr::select(prot_acc), 
	  acc_lookup %>% 
	    dplyr::select(-which(names(.) %in% names(df))) %>% 
	    dplyr::select(-"organism", -"fasta_name"), 
	) 
	
	# (1) multiple uniprot_acc(s) can share the same entrez id
	# prot_acc    gene      acc_type   uniprot_acc species entrez
	# A1AT1_MOUSE Serpina1a uniprot_id P07758      mouse    20703
	# A1AT4_MOUSE Serpina1d uniprot_id Q00897      mouse    20703
		
	# (2) same uniprot_acc can have different entrez ids
	# From	To
	# P02088	100503605
	# P02088	101488143
	# P02088	15129	
	
	# (3) ok that some uniprot_accs(s) have no corresponding entrez id

	df <- df %>% 
	  dplyr::mutate(psm_index = row_number()) %>% 
	  dplyr::left_join(acc_lookup, by = "prot_acc") %>% 
	  dplyr::filter(!duplicated(psm_index)) %>% 
	  dplyr::select(-psm_index)

	return(df)
}


#' Adds kinase annotation
#'
#' @inheritParams info_anal
#' @param acc_type Character string; the type of protein accessions in one of
#'   c("refseq_acc", "uniprot_acc", "uniprot_id")
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
annotKin <- function (df, acc_type) {
	stopifnot ("prot_acc" %in% names(df))
	
  data(package = "proteoQ", kinase_lookup)
  stopifnot(acc_type %in% names(kinase_lookup))

  lookup <- kinase_lookup %>% 
    dplyr::select(acc_type, kin_attr, kin_class, kin_order) %>%
    dplyr::filter(!duplicated(.)) %>% 
    dplyr::filter(!is.na(.[[acc_type]])) %>% 
    dplyr::mutate(!!acc_type := as.character(!!rlang::sym(acc_type)))

  df <- df %>% 
    dplyr::left_join(lookup, by = c("prot_acc" = acc_type))

	df$kin_attr[is.na(df$kin_attr)] <- FALSE

	return(df)
}


#' Saves the arguments in a function call
#'
#' @param call_pars Language.
#' @param fn The name of function being saved.
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
save_call <- function(call_pars, fn) {
	dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)
	call_pars[names(call_pars) == "..."] <- NULL
	save(call_pars, file = file.path(dat_dir, "Calls", paste0(fn, ".rda")))
}


#' Matches formulas to those in calls to pepSig or prnSig
#'
#' @param formulas Language; the formulas in linear modeling.
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_fmls <- function(formulas) {
  fml_file <-  file.path(dat_dir, "Calls\\pepSig_formulas.rda")
  
  if (file.exists(fml_file)) {
    load(file = fml_file)
  } else {
    stop("Run `pepSig()` first.")
  }
  
  fml_chr <- formulas %>%
    as.character() %>%
    gsub("\\s+", "", .)
  
  pepSig_chr <- pepSig_formulas %>%
    purrr::map(~ .[is_call(.)]) %>%
    as.character() %>%
    gsub("\\s+", "", .)
  
  ok <- purrr::map_lgl(fml_chr, ~ . %in% pepSig_chr)
  
  if (!all(ok))
    stop("Formula match failed: ", formulas[[which(!ok)]],
         " not found in the latest call to 'pepSig(...)'.")
}


#' Matches the arg to anal_prnGSPA
#'
#' @param call_rda the name of a rda.
#' @param arg Argument to be matched.
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
match_call_arg <- function (call_rda = "foo", arg = "scale_log2r") {
  call_rda <- rlang::as_string(rlang::enexpr(call_rda))
  arg <- rlang::as_string(rlang::enexpr(arg))
  
  rda <- paste0(call_rda, ".rda")
  file <- file.path(dat_dir, "Calls", rda)
  if (!file.exists(file)) stop(rda, " not found.")
  
  load(file = file)
  if (is.null(call_pars[[arg]])) 
    stop(arg, " not found in the latest call to ", call_rda, call. = FALSE)
  
  call_pars[[arg]]
}


#' Matches the name of GSPA result file (not currently used)
#'
#' @param anal_type Always \code{GSPA}; maybe different value for future uses.
#' @param subdir Character string of sub directory
#' @inheritParams prnHM
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
match_gspa_filename <- function (anal_type = "GSPA", subdir = NULL, scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(!is.null(subdir))

  if (anal_type == "GSPA") {
    file <- file.path(dat_dir, "Calls\\anal_prnGSPA.rda")
    if (!file.exists(file)) stop("Run `prnGSPA` first.", call. = FALSE)
    load(file = file)
  }
  
  filename <- call_pars$filename
  if (rlang::is_empty(filename)) {
    filename <- list.files(path = file.path(dat_dir, "Protein\\GSPA", subdir), 
                           pattern = "^Protein_GSPA_.*\\.txt$", 
                           full.names = FALSE)
    
    if (rlang::is_empty(filename)) stop("No result file found under `", subdir, "`", call. = FALSE)
    
    if (scale_log2r) filename <- filename %>% .[grepl("^Protein_GSPA_Z", .)] else 
      filename <- filename %>% .[grepl("^Protein_GSPA_N", .)]
    
    if (impute_na) filename <- filename %>% .[grepl("Protein_GSPA_[NZ]_impNA\\.txt", .)] else 
      filename <- filename %>% .[grepl("Protein_GSPA_[NZ]\\.txt", .)]

    if (length(filename) > 1) stop("More than one result file found under `", subdir, "`", call. = FALSE)
    if (rlang::is_empty(filename)) 
      stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, call. = FALSE)
  }

  return(filename)
}


#' Matches gset_nms to prnGSPA (not currently used)
#' 
#' @inheritParams prnGSPA
match_gset_nms <- function (gset_nms = NULL) {
  file <- file.path(dat_dir, "Calls\\anal_prnGSPA.rda")
  
  if (is.null(gset_nms)) {
    if (file.exists(file)) {
      gset_nms <- match_call_arg(anal_prnGSPA, gset_nms)
    } else {
      stop("`gset_nms` is NULL and ", file, " not found.", call. = FALSE)
    }
  } else {
    if (file.exists(file)) {
      gset_nms <- gset_nms %>% 
        .[. %in% match_call_arg(anal_prnGSPA, gset_nms)]
    } else {
      warning("File `", file, "`` not available for `gset_nms` matches.", 
              "\nUse the `gset_nms` value as it.", call. = FALSE)
    }
  }
  
  if (is.null(gset_nms)) {
    stop ("The `gset_nms` is NULL after matching parameters to the latest `prnGSPA(...)`.", 
          "\n\tConsider providing explicitly the `gset_nms`: ", 
          "\n\t\t`gset_nms = \"go_sets\"` or ", 
          "\n\t\t`gset_nms = \"~\\\\proteoQ\\\\dbs\\\\go_hs.rds\"` ...",
          call. = FALSE)
  }
  
  if (purrr::is_empty(gset_nms)) {
    stop ("The `gset_nms` is EMPTY after parameter matching.", 
          "\n\tCheck the values of `gset_nms` in the latest call to `prnGSPA(...)`.", 
          call. = FALSE)
  }
  
  return(gset_nms)
}


#' Match scale_log2r
#'
#' \code{match_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in caller environment.
#' 
#' @inheritParams prnHist
match_scale_log2r <- function(scale_log2r) {
  stopifnot(rlang::is_logical(scale_log2r))
  
  global_var <-tryCatch(global_var <-get("scale_log2r", envir = .GlobalEnv),
                        error = function(e) "e")
  if (global_var != "e" & is.logical(global_var)) scale_log2r <- global_var
  
  return(scale_log2r)
}


#' Match to a global logical variable
#' 
#' @param var Character string representation of a variable.
#' @param val The value of \code{var} before matching.
#' @examples
#' \donttest{
#' foo <- function(scale_log2r = FALSE) {
#'   match_logi_gv("scale_log2r", scale_log2r)
#' }
#' 
#' scale_log2r <- TRUE
#' foo()
#' }
match_logi_gv <- function(var, val) {
  oval <- val
  gvar <-tryCatch(gvar <- get(var, envir = .GlobalEnv), error = function(e) "e")
  
  if (gvar != "e") {
    stopifnot(rlang::is_logical(gvar))
    if (gvar != oval) {
      warning("Coerce ", var, " to ", gvar, " after matching to the global setting.", 
              call. = FALSE)
    }
    return(gvar)
  } else {
    return(val)
  }
}


#' Match scale_log2r 
#'
#' \code{match_prnSig_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in the most recent prnSig at a given impute_na status
#' 
#' @inheritParams prnHM
match_prnSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  if (impute_na) {
    mscale_log2r <- match_call_arg(prnSig_impTRUE, scale_log2r)
  } else {
    mscale_log2r <- match_call_arg(prnSig_impFALSE, scale_log2r)
  }
  
  if (scale_log2r != mscale_log2r) {
    warning("scale_log2r = ", mscale_log2r, " after matching to prnSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }
  
  scale_log2r <- mscale_log2r
}


#' Match scale_log2r 
#'
#' \code{match_pepSig_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in the most recent pepSig at a given impute_na status.
#' @inheritParams prnHM
match_pepSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  if (impute_na) {
    mscale_log2r <- match_call_arg(pepSig_impTRUE, scale_log2r)
  } else {
    mscale_log2r <- match_call_arg(pepSig_impFALSE, scale_log2r)
  }
  
  if (scale_log2r != mscale_log2r) {
    warning("scale_log2r = ", mscale_log2r, " after matching to pepSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }
  
  scale_log2r <- mscale_log2r
}


#' Replaces NA genes
#' 
#' @inheritParams info_anal
#' @inheritParams annotKin
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
replace_na_genes <- function(df, acc_type) {
	acc_type <- tolower(acc_type)

	if (acc_type == "refseq_acc") {
		na_gene <- (is.na(df[, c("gene")])) | (str_length(df$gene) == 0)
		df$gene <- as.character(df$gene)
		df$gene[na_gene] <- as.character(df$prot_acc[na_gene])
	} else if (acc_type == "uniprot_id") {
		temp <- data.frame(do.call('rbind', strsplit(as.character(df$prot_desc), 'GN=', fixed = TRUE))) %>%
				dplyr::select(2) %>%
				`colnames<-`("gene") %>%
				dplyr::mutate(gene = gsub("PE\\=.*", "", gene)) %>%
				dplyr::mutate(gene = gsub("\\s+.*", "", gene))

		na_gene <- is.na(df$gene)
		df[na_gene, c("gene")] <- temp[na_gene, c("gene")]
		rm(temp)

		df <- df %>%
			dplyr::mutate(gene = gsub(".*\\|", "", gene))
	}

	return(df)
}


#' Find peptide start and end positions
#'
#' \code{find_pep_pos} finds the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}.
#' 
#' @param prot_acc Protein accession
#' @param pep_seq Peptide sequence
#' @inheritParams normPSM
#' @import dplyr purrr rlang stringr seqinr tidyr
#' @importFrom magrittr %>% %$%
find_pep_pos <- function (prot_acc, pep_seq, fasta) {
  fasta_sub <- fasta %>% .[names(.) == prot_acc]
  pep_seq <- as.character(pep_seq)
  
  if (!rlang::is_empty(fasta_sub)) {
    pep_pos <- stringr::str_locate(fasta_sub, pattern = pep_seq)
    
    pos_bf <- pep_pos[1] - 1
    pos_af <- pep_pos[2] + 1
    
    pep_res_before <- stringr::str_sub(fasta_sub, pos_bf, pos_bf)
    pep_res_after <- stringr::str_sub(fasta_sub, pos_af, pos_af)
    
    # Mascot can alter the original sequence in fasta
    # prot_acc: "XP_003960355", original "QERFCQXK" becomes "QERFCQVK"
    if (any(is.na(c(pep_res_before, pep_res_after)))) {
      pep_pos <- cbind(pep_seq, pep_res_before = NA, start = NA, end = NA, 
                       pep_res_after = NA, prot_acc = prot_acc, is_tryptic = NA)
      
      return(pep_pos)
    }
    
    if (nchar(pep_res_before) == 0) pep_res_before <- "-"
    if (nchar(pep_res_after) == 0) pep_res_after <- "-"
    
    # ADVSLPSMQGDLK|NP_612429: not "E.ADVSLPSMQGDLK.T" but "K.ADVSLPSMQGDLK.T"
    if (pep_res_before %in% c("K", "R", "-")) { # the first match is tryptic
      pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, pep_res_after, prot_acc, is_tryptic = TRUE)
    } else if (pep_res_before == "M" & pep_pos[1] == 2) { # the first match is also tryptic
      pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, pep_res_after, prot_acc, is_tryptic = TRUE)
    } else { # the first match is non-tryptic
      pep_seq_new <- paste0(c("K", "R"), pep_seq)
      pep_pos_new_all <- purrr::map(pep_seq_new, ~ str_locate(fasta_sub, .x))
      ok_pos <- purrr::map_lgl(pep_pos_new_all, ~ !is.na(.x[[1]]))
      
      if (sum(ok_pos) > 0) { # tryptic match existed
        pep_pos_new <- pep_pos_new_all[[which(ok_pos)[1]]]
        
        pos_bf_new <- pep_pos_new[1]
        pos_af_new <- pep_pos_new[2] + 1
        
        pep_res_before_new <- stringr::str_sub(fasta_sub, pos_bf_new, pos_bf_new)
        pep_res_after_new <- stringr::str_sub(fasta_sub, pos_af_new, pos_af_new)
        
        pep_pos_new[1] <- pep_pos_new[1] + 1
        
        pep_pos <- cbind(pep_seq, pep_res_before_new, pep_pos_new, pep_res_after_new, prot_acc, is_tryptic = TRUE)        
      } else { # no tryptic matches
        pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, pep_res_after, prot_acc, is_tryptic = FALSE)
      }
    }
  } else { # no fasta matches
    pep_pos <- cbind(pep_seq, pep_res_before = NA, start = NA, end = NA, 
                     pep_res_after = NA, prot_acc = prot_acc, is_tryptic = FALSE)
  }
}


#' Annotation of peptide positions and adjacent amino acid residues
#'
#' \code{annotPeppos} annotates the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}. It also annotates the
#' preceding and the following AA residues.
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @import dplyr purrr rlang stringr seqinr tidyr
#' @importFrom magrittr %>% %$%
annotPeppos <- function (df, fasta){
  stopifnot(all(c("prot_acc", "pep_seq") %in% names(df)))
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)
  
  load(file = file.path(dat_dir, "label_scheme.rda"))

  # ok cases that same `pep_seq` but different `prot_acc`
  # (K)	MENGQSTAAK	(L) NP_510965
  # (-)	MENGQSTAAK	(L) NP_001129505  
  
  if (! "pep_seq_bare" %in% names(df)) {
    df <- df %>% 
      dplyr::mutate(pep_seq_bare = gsub("^.*\\.([^\\.]+)\\..*", "\\1", pep_seq))
  }

  df <- df %>% 
    dplyr::mutate(prot_acc = gsub("^.*\\|(.*)\\|.*$", "\\1", prot_acc)) %>% 
    dplyr::mutate(pep_prn = paste(pep_seq_bare, prot_acc, sep = "|"))

  df_pep_prn <- df %>% 
    dplyr::filter(!duplicated(pep_prn)) %>% 
    dplyr::select(c("pep_seq_bare", "prot_acc")) 
  
  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- purrr::map(fasta, ~ {
        seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
      }) %>% do.call(`c`, .) %>% 
        `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
        .[names(.) %in% unique(df_pep_prn$prot_acc)]
      
      if (length(fasta) == 0) {
        stop("No fasta entries matched protein accessions; probably wrong fasta file(s).", 
             call. = FALSE)
      }
      
      pep_pos_all <- purrr::map2(as.list(df_pep_prn$prot_acc), as.list(df_pep_prn$pep_seq_bare), 
                                 find_pep_pos, fasta) %>% 
        do.call(rbind, .) %>% 
        `colnames<-`(c("pep_seq_bare", "pep_res_before", "pep_start", "pep_end", "pep_res_after", 
                       "prot_acc", "is_tryptic")) %>% 
        data.frame(check.names = FALSE) %>% 
        tidyr::unite(pep_prn, pep_seq_bare, prot_acc, sep = "|", remove = TRUE)
    } else {
      stop("Wrong FASTA file path or name(s).", call. = FALSE)
    }
  } else {
    stop("FASTA file(s) not provided.")
  }
  
  rm(fasta)
  
  if ("pep_res_before" %in% names(df)) pep_pos_all$pep_res_before <- NULL
  if ("pep_res_after" %in% names(df)) pep_pos_all$pep_res_after <- NULL
  if ("pep_start" %in% names(df)) pep_pos_all$pep_start <- NULL
  if ("pep_end" %in% names(df)) pep_pos_all$pep_end <- NULL
  
  df$pep_seq_bare <- NULL

  df <- df %>% 
    dplyr::left_join(pep_pos_all, by = "pep_prn") %>% 
    dplyr::select(-pep_prn)
}


#' Subset fasta by accession type
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @inheritParams annotKin
#' @import plyr dplyr purrr rlang seqinr
#' @importFrom magrittr %>%
subset_fasta <- function (df, fasta, acc_type) {
  stopifnot("prot_acc" %in% names(df))
  
  if (! acc_type %in% c("refseq_acc", "uniprot_id", "uniprot_acc")) {
    stop("The type of protein accesion needs to one of \'uniprot_id\', \'uniprot_acc\' or \'refseq_acc\'",
         call. = FALSE)
  }
  
  fasta <- purrr::map(fasta, ~ {
    seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
  }) %>% do.call(`c`, .)
  
  if (acc_type == "uniprot_id") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|.*\\|(.*)$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } else if (acc_type == "uniprot_acc") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } else if (acc_type == "refseq_acc") {
    fasta <- fasta %>% 
      .[names(.) %in% unique(df$prot_acc)]
  }    
}


#' Calculates protein percent coverage
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @import plyr dplyr purrr rlang seqinr
#' @importFrom magrittr %>%
calc_cover <- function(df, id, fasta = NULL) {
  stopifnot(all(c("prot_acc", "gene", "pep_start", "pep_end") %in% names(df)))

  if (all(is.factor(df$pep_start))) {
    df$pep_start <- df$pep_start %>% as.character() %>% as.numeric()
  }
    
  if (all(is.factor(df$pep_end))) {
    df$pep_end <- df$pep_end %>% as.character() %>% as.numeric()
  }

  id <- rlang::as_string(rlang::enexpr(id))
  if (id == "gene") {
    gn_rollup <- TRUE
    id <- "prot_acc"
  } else {
    gn_rollup <- FALSE
  }
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "acc_lookup.rda"))
  
  if (length(fasta) == 0) {
    stop("No fasta entries matched the type of protein accession. Check the correctness of fasta file(s).", 
         call. = FALSE)
  }
  
  if (length(fasta) <= 200) {
    warning("Less than 200 entries in fasta matched by protein accession. 
            Make sure the fasta file is correct.")
  }
  
  df_sels <- df %>%
    dplyr::select(prot_acc, pep_start, pep_end) %>%
    dplyr::mutate(index = row_number()) %>% 
    dplyr::left_join(acc_lookup, by = "prot_acc") %>%
    dplyr::filter(!is.na(prot_len), !duplicated(index)) %>% 
    dplyr::select(-index)
  
  if (nrow(df_sels) == 0) stop("Probably incorrect accession types in the fasta file(s).", call. = FALSE)
  
  df_sels <- df_sels %>%
    dplyr::filter(pep_start <= prot_len) %>%
    dplyr::filter(pep_end <= prot_len) %>%
    split(.[["prot_acc"]], drop = TRUE) %>%
    purrr::map(function (.x) {
      len <- .x[1, "prot_len"]
      aa_map <- rep(NA, len)
      for (i in 1:nrow(.x)) aa_map[.x[i, ]$pep_start : .x[i, ]$pep_end] <- TRUE
      sum(aa_map, na.rm = TRUE)/len
    } ) %>%
    do.call("rbind", .) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`("prot_cover") %>%
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::mutate(prot_cover = ifelse(prot_cover > 1, 1, prot_cover)) 
  
  if (gn_rollup) {
    df_sels <- df %>% 
      dplyr::select(prot_acc, gene) %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_sels, by = "prot_acc") %>% 
      dplyr::select(-prot_acc) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::summarise_all(~ max(.x, na.rm = TRUE))
    
    df_sels <- df %>% 
      dplyr::select(prot_acc, gene) %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_sels, by = "gene") %>% 
      dplyr::select(-gene) 
  }
  
  df <- df %>% 
    dplyr::mutate(index = row_number()) %>% 
    dplyr::left_join(df_sels, by = "prot_acc") %>% 
    dplyr::filter(!duplicated(index)) %>% 
    dplyr::select(-index) %>% 
    dplyr::mutate(prot_cover = round(prot_cover * 100, digits = 1)) %>%
    dplyr::mutate(prot_cover = paste0(prot_cover, "%"))

  return(df)
}


#' Converts log2FC to linear fold changes
#' 
#' @inheritParams info_anal
#' @import dplyr purrr
#' @importFrom magrittr %>%
to_linfc <- function(df) {
		nms <- rownames(df)

		df %>%
			purrr::map(~ {ifelse(.x > 0, 2^.x, -1/(2^.x))}) %>%
			data.frame(check.names = FALSE) %>%
			`rownames<-`(nms)
}


#' Remove single-value columns
#' 
#' @inheritParams setHMlims
#' @import dplyr purrr
#' @importFrom magrittr %>%
rm_sglval_cols <- function (x) {
  sgl_val <- x %>% 
    summarise_all(funs(n_distinct(.))) %>% 
    purrr::map(~ .x == 1) %>% 
    purrr::flatten_lgl()
  
  x[, !sgl_val, drop = FALSE]
}


#' Combine data with metadata
#' 
#' @param data A data frame
#' @param metadata Another data frame 
#' @import dplyr purrr
#' @importFrom magrittr %>%
cmbn_meta <- function(data, metadata) {
  suppressWarnings(
    data %>% 
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::left_join(metadata, by = "Sample_ID") %>%
      dplyr::mutate_at(vars(one_of("Color", "Fill", "Shape", "Size", "Alpha")), ~ as.factor(.)) %>%
      dplyr::select(which(not_all_NA(.))) %>% 
      rm_sglval_cols()    
  )
}


#' Check file names for ggsave()
#' @param filename Character string; An output file name.
gg_imgname <- function(filename) {
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  exts <- c("png", "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg") 
  
  if(! fn_suffix %in% exts) {
    warning(paste0("Unrecognized file extenstion: '", fn_suffix, "'. Image will be saved as a '.png'.\n"))
    fn_suffix <- "png"
  }
  
  paste0(fn_prefix, ".", fn_suffix)
}


#' Check file names for ggsave()
#' 
#' @inheritParams info_anal
#' @import dplyr purrr
#' @importFrom magrittr %>%
rm_pval_whitespace <- function(df) {
  df <- df %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.character) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ gsub("\\s*", "", .x) ) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ suppressWarnings(as.numeric(.x)))
}


#' Filter rows
#'
#' @param df a data frame. 
#' @param ... Arguments for \link[dplyr]{filter}
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
filters_in_call <- function (df, ...) {
  dots <- rlang::enexprs(...)
  nms <- names(dots)
  
  for (i in seq_along(dots)) {
    row_exprs <- dots[[nms[i]]] %>% 
      rlang::eval_bare()
    
    if (!rlang::is_list(row_exprs)) row_exprs <- list(row_exprs)
    
    row_vals <- row_exprs %>% 
      purrr::map(eval_tidy, df) %>% 
      purrr::reduce(`&`, .init = 1)
    
    row_vals <- row_vals %>% ifelse(is.na(.), FALSE, .)

    stopifnot(is.logical(row_vals))
    
    df <- df[row_vals, , drop = FALSE]
  }
  
  return(df)
}


#' Arrange rows
#'
#' @param .df a data frame. 
#' @param .na.last The same as \link[base]{order}
#' @param ... Arguments for \link[dplyr]{arrange}
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>% 
arrangers_in_call <- function(.df, ..., .na.last = TRUE) {
  dots <- rlang::enexprs(...)
  nms <- names(dots)
  
  for (i in seq_along(dots)) {
    row_orders <- dots[[nms[i]]] %>% 
      rlang::eval_bare()
    
    if (!rlang::is_list(row_orders)) row_orders <- list(row_orders)
    
    order_call <- rlang::expr(order(!!!row_orders, na.last = !!.na.last))
    ord <- rlang::eval_tidy(order_call, .df)
    stopifnot(length(ord) == nrow(.df))

    .df <- .df[ord, , drop = FALSE]
  }
  
  return(.df)
}


#' Calculate PSM SDs
#' 
#' @inheritParams info_anal
#' @inheritParams standPep
#' @inheritParams channelInfo
#' @inheritParams calcPepide
calc_sd_fcts_psm <- function (df, range_log2r, range_int, set_idx, injn_idx) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  label_scheme <- label_scheme %>% 
    dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx) %>% 
    dplyr::mutate(tmt_nm = gsub("TMT-", "N_log2_R", TMT_Channel))
  
  label_scheme_sd <- label_scheme %>%
    dplyr::filter(!Reference, !grepl("^Empty\\.", Sample_ID))	%>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = (.$Sample_ID)))
  
  SD <- df %>%
    dplyr::select(grep("^N_log2_R|^N_I", names(.))) %>%
    dblTrim(., range_log2r, range_int) %>%
    `names<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.)))
  
  cf_SD <- SD/mean(SD %>% .[names(.) %in% label_scheme_sd$tmt_nm], na.rm = TRUE)
  cf_SD <- cbind.data.frame(fct = cf_SD, SD) %>%
    tibble::rownames_to_column("tmt_nm") %>%
    dplyr::mutate(tmt_nm = factor(tmt_nm, levels = label_scheme$tmt_nm)) %>%
    dplyr::arrange(tmt_nm)
}


#' Calculate CV per TMT_Set and LCMS_injection
#' @inheritParams info_anal
#' @param type Character string; the type of data.
calcSD_Splex <- function (df, id, type = "log2_R") {
  if (type == "log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^log2_R[0-9]{3}", names(.)))
  } else if (type == "N_log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^N_log2_R[0-9]{3}", names(.)))
  } else if (type == "Z_log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^Z_log2_R[0-9]{3}", names(.)))
  }
  
  run_scripts <- FALSE
  if (run_scripts) {
    df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::group_by(!!rlang::sym(id)) %>%
      dplyr::summarise_at(vars(starts_with(type)), ~ sd(.x, na.rm = TRUE))     
  }

  df %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id))) %>% 
    dplyr::arrange(!!rlang::sym(id)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarise_at(vars(starts_with(type)), ~ sd(.x, na.rm = TRUE)) # %>% 
    # dplyr::mutate(!!id := as.factor(!!rlang::sym(id)))
}



#' Violin plots of CV per TMT_Set and LCMS_injection
#' 
#' @param width The width of a plot.
#' @param height The height of a plot.
#' @param is_psm Logical; indicator if the data belong to a PSM table .
#' 
#' @inheritParams info_anal
#' @inheritParams purgePSM
#' @inheritParams prnCorr_logFC
#' @inheritParams calcSD_Splex
#' @import dplyr rlang ggplot2
sd_violin <- function(df = NULL, id = NULL, filepath = NULL, width = NULL, height = NULL, 
                      type = "log2_R", adjSD = FALSE, 
                      is_psm = FALSE, col_select = NULL, col_order = NULL, theme = NULL, ...) {

  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  col_select <- ifelse(is.null(col_select), rlang::expr(Select), rlang::sym(col_select))
  col_order <- ifelse(is.null(col_order), rlang::expr(Order), rlang::sym(col_order))
  
  if (col_select == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
  if (col_order == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  if (is.null(label_scheme[[col_select]])) {
    stop("Column \'", rlang::as_string(col_select), "\' does not exist.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    stop("No samples were selected under column \'", rlang::as_string(col_select), "\'.",
         call. = FALSE)
  }
  
  if (is.null(label_scheme[[col_order]])) {
    warning("Column \'", rlang::as_string(col_order), "\' does not exist.
			Samples will be arranged by the alphebatic order.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_order]])) == 0) {
    warning("No samples were specified under column \'", rlang::as_string(col_order), "\'.",
            call. = FALSE)
  }
  
  label_scheme_sub <- label_scheme %>% 
    dplyr::mutate(new_id = paste0(TMT_Channel, " (", Sample_ID, ")")) %>% 
    dplyr::mutate(new_id = gsub("TMT-", "", new_id)) %>% 
    dplyr::select(Sample_ID, TMT_Set, new_id, !!col_select, !!col_order) %>%
    dplyr::filter(!is.na(!!col_select))

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  if (rlang::is_missing(width)) width <- 8
  if (rlang::is_missing(height)) height <- 8
  
  ymax <- eval(dots$ymax, env = caller_env())
  ybreaks <- eval(dots$ybreaks, env = caller_env())
  
  flip_coord <- eval(dots$flip_coord, env = caller_env())
  if (is.null(flip_coord)) flip_coord <- FALSE
  
  df <- df %>% dplyr::filter(!duplicated(.[[id]]))

  if (type == "log2_R") {
    df_sd <- df %>% dplyr::select(id, grep("^sd_log2_R[0-9]{3}[NC]*", names(.)))
  } else if (type == "N_l og2_R") {
    df_sd <- df %>% dplyr::select(id, grep("^sd_N_log2_R[0-9]{3}[NC]*", names(.)))
  } else if (type == "Z_log2_R") {
    df_sd <- df %>% dplyr::select(id, grep("^sd_Z_log2_R[0-9]{3}[NC]*", names(.)))
  }
  
  # all-NA first removed for finding all-NaN columns
  df_sd <- df_sd %>% 
    dplyr::filter(rowSums(!is.na(.[grep("^.*log2_R[0-9]{3}", names(.))])) > 0)
  
  if (adjSD) {
    SD <- df %>%
      dplyr::select(grep("^log2_R[0-9]{3}|^I[0-9]{3}", names(.))) %>%
      dblTrim(., range_log2r = c(0, 100), range_int = c(0, 100), type_r = "log2_R", type_int = "I")
    
    df_sd[, grep("^.*log2_R", names(df_sd))] <- df_sd[, grep("^.*log2_R", names(df_sd)), drop = FALSE] %>% 
      sweep(., 2, sqrt(SD), "/")
    
    df_z <- df_sd %>% dplyr::select(grep("^.*log2_R[0-9]{3}", names(.)))
    nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
    df_z[, nan_cols] <- 0
    df_sd[, grep("^.*_log2_R[0-9]{3}", names(df_sd))] <- df_z
    
    rm(df_z, nan_cols, SD)
  }

  df_sd <- df_sd %>% 
    `names<-`(gsub("^.*log2_R", "", names(.))) 
  
  if (!is_psm) {
    df_sd <- df_sd %>% dplyr::select(id, which(names(.) %in% label_scheme_sub$new_id))
  }

  Levels <- names(df_sd) %>% .[! . %in% id]

  if (!purrr::is_empty(Levels)) {
    df_sd <- df_sd %>%
      tidyr::gather(key = !!rlang::sym(id), value = "SD") %>%
      dplyr::rename(Channel := !!rlang::sym(id)) %>% 
      dplyr::ungroup(Channel) %>% 
      dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
      dplyr::filter(!is.na(SD))
    
    p <- ggplot() +
      geom_violin(df_sd, mapping = aes(x = Channel, y = SD, fill = Channel), size = .25, draw_quantiles = c(.95, .99)) +
      geom_boxplot(df_sd, mapping = aes(x = Channel, y = SD), width = 0.1, lwd = .2, fill = "white") +
      stat_summary(df_sd, mapping = aes(x = Channel, y = SD), fun.y = "mean", geom = "point",
                   shape=23, size=2, fill="white", alpha=.5) +
      labs(title = expression(""), x = expression("Channel"), y = expression("SD ("*log[2]*"FC)")) 
    
    if (!is.null(ymax)) {
      if (is.null(ybreaks)) {
        ybreaks <- ifelse(ymax > 1, 0.5, ifelse(ymax > 0.5, 0.2, 0.1))
      }
      p <- p + scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, ybreaks))
    }

    if (is.null(theme)) theme <- theme_psm_violin
    p <- p + theme

    if (flip_coord) {
      p <- p + coord_flip()
      width_temp <- width
      width <- height
      height <- width_temp
      rm(width_temp)
    }
    
    dots <- dots %>% .[! names(.) %in% c("width", "height", "in", "limitsize")]
    my_call <- rlang::expr(ggplot2::ggsave(filename = !!filepath, plot = !!p, 
                                           width = !!width, height = !!height, 
                                  units = "in", limitsize = FALSE, !!!dots))
    try(eval(my_call, caller_env()))
  }
}


#' Violin plots of reporter-ion intensity per TMT_Set and LCMS_injection
#' 
#' @inheritParams info_anal
#' @inheritParams sd_violin
rptr_violin <- function(df, filepath, width, height) {
  df_int <- df %>% 
    `names<-`(gsub("^N_I|^I", "", names(.))) 
  
  Levels <- names(df_int)
  
  df_int <- df_int %>%
    tidyr::gather(key = "Channel", value = "Intensity") %>%
    dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
    dplyr::filter(!is.na(Intensity))
  
  mean_int <- df_int %>% 
    dplyr::group_by(Channel) %>% 
    dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
    dplyr::mutate(Intensity = round(Intensity, digit = 1))
  
  p <- ggplot() +
    geom_violin(df_int, mapping = aes(x = Channel, y = log10(Intensity), fill = Channel), size = .25) +
    geom_boxplot(df_int, mapping = aes(x = Channel, y = log10(Intensity)), width = 0.2, lwd = .2, fill = "white") +
    stat_summary(df_int, mapping = aes(x = Channel, y = log10(Intensity)), fun.y = "mean", geom = "point",
                 shape = 23, size = 2, fill = "white", alpha = .5) +
    labs(title = expression("Reporter ions"), x = expression("Channel"), y = expression("Intensity ("*log[10]*")")) + 
    geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
    theme_psm_violin
  
  try(ggplot2::ggsave(filepath, p, width = width, height = height, units = "in"))
}


#' geometric mean
#' 
#' @param x A data frame.
#' @param ... The same in \code{mean}.
my_geomean <- function (x, ...) {
  x <- log10(x) %>% mean(...)
  10^x
}


#' phospho counts
count_phosphopeps <- function() {
  df <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
  
  id <- match_call_arg(normPSM, group_psm_by)
  
  df_phos <- df %>% dplyr::filter(grepl("[sty]", .[[id]]))
  
  n_phos_peps <- nrow(df_phos)
  n_phos_sites <- stringr::str_count(df_phos[[id]], "[sty]") %>% sum()
  
  write.csv(
    data.frame(n_peps = n_phos_peps, n_sites = n_phos_sites), 
    file.path(dat_dir, "Peptide\\cache", "phos_pep_nums.csv"), 
    row.names = FALSE
  )
}


#' peptide mis-cleavage counts
count_pepmiss <- function() {
  dir.create(file.path(dat_dir, "PSM\\cache"), recursive = TRUE, showWarnings = FALSE)
  
  rmPSMHeaders()
  
  filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
                        pattern = "^F[0-9]{6}\\_hdr_rm.csv$")
  
  if (length(filelist) == 0) stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
  df <- purrr::map(filelist, ~ {
    data <- read.delim(file.path(dat_dir, "PSM\\cache", .x), sep = ',', check.names = FALSE, 
                       header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    
    data$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    data <- data %>% 
      dplyr::filter(!duplicated(pep_seq))
    
    tot <- nrow(data)
    mis <- data %>% dplyr::filter(pep_miss > 0) %>% nrow()
    
    tibble::tibble(total = tot, miscleavage = mis, percent = miscleavage/tot)
  }) %>% do.call(rbind, .)
  
  write.csv(df, file.path(dat_dir, "PSM\\cache\\miscleavage_nums.csv"), row.names = FALSE)
}


#' Row filtration helpers
#'
#' \code{contain_str}: contain a literal string; "PEPTIDES" contain_str "TIDE".
#' 
#' @param match A character string containing the pattern for matching.
#' @param vars A character string of the name of a variable. The default is
#'   FALSE.
#' @param ignore.case Logical; if TRUE, ignores case when matching.
#' @examples
#' \donttest{
#' pepHist(
#'   col_select = BI,
#'   scale_log2r = TRUE,
#'   filter_peps = exprs(contain_chars_in("sty", pep_seq_mod)),
#'   scale_y = FALSE,
#'   ncol = 4,
#'   filename = "BI_pSTY_scaley_no.png",
#' )
#' }
#' @export
contain_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(match, vars, fixed = TRUE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{contain_chars_in}: contain some of the characters in a literal string;
#' "PEPTIDES" contain_chars_in "XP".
#' 
#' @rdname contain_str
#' @export
contain_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{not_contain_str}" not contain a literal string; "PEPTIDES"
#' not_contain_str "TED".
#' 
#' @rdname contain_str
#' @export
not_contain_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(match, vars, fixed = TRUE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{not_contain_chars_in}: not contain any of the characters in a literal
#' string; "PEPTIDES" not_contain_chars_in  "CAB".
#' 
#' @rdname contain_str
#' @export
not_contain_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}

#' Row filtration helpers
#'
#' \code{start_with_str}: start with a literal string. "PEPTIDES" start_with_str
#' "PEP".
#' 
#' @rdname contain_str
#' @export
start_with_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("^", match), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{end_with_str}: end with a literal string. "PEPTIDES" end_with_str
#' "TIDES".
#' 
#' @rdname contain_str
#' @export
end_with_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0(match, "$"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{start_with_chars_in}: start with one of the characters in a literal
#' string. "PEPTIDES" start_with_chars_in "XP".
#' 
#' @rdname contain_str
#' @export
start_with_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("^[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{ends_with_chars_in}: end with one of the characters in a literal
#' string. "PEPTIDES" ends_with_chars_in "XS".
#' 
#' @rdname contain_str
#' @export
ends_with_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[", match, "]$"), vars, fixed = FALSE, ignore.case)
}


#' Row filtration helpers
#'
#' \code{rows_are_all}: rows are all
#' @rdname contain_str
rows_are_all <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}


#' Row filtration helpers
#'
#' \code{rows_are_all}: rows are all
#' @rdname contain_str
rows_are_not_all <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}


#' Concatenate formula(s) to varargs of dots
#' 
#' @param fmls A character vector of formula(s)
#' @param dots A character vector of formula(s) in \code{dots}
#' @param fml_nms A character vector containing the names of \code{fmls}.
#' @inheritParams info_anal
concat_fml_dots <- function(fmls = NULL, fml_nms = NULL, dots = NULL, anal_type = "zzz") {
  if ((!is_empty(fmls)) & (anal_type == "GSEA")) return(c(dots, fmls))
  
  if (purrr::is_empty(fmls)) {
    fml_file <-  file.path(dat_dir, "Calls\\pepSig_formulas.rda")
    if (file.exists(fml_file)) {
      load(file = fml_file)
      
      if (!is.null(fml_nms)) {
        stopifnot(all(fml_nms %in% names(pepSig_formulas)))
        # stopifnot(all(names(pepSig_formulas) %in% fml_nms))
        # pepSig_formulas <- pepSig_formulas %>% .[map_dbl(names(.), ~ which(.x == fml_nms))]
        pepSig_formulas <- pepSig_formulas %>% .[names(.) %in% fml_nms]
      }
      
      dots <- c(dots, pepSig_formulas)
    } else {
      stop("Run both `pepSig()` and `prnSig()` first.")
    }
  } else {
    match_fmls(fmls)
    dots <- c(dots, fmls)
  }
  
  return(dots)
}


#' Roll up genes
#' 
#' @param df A data frame
#' @param cols Column indexes
gn_rollup <- function (df, cols) {
  if (! "gene" %in% names(df)) return(df)
  
  dfa <- df %>% 
    dplyr::select(gene, cols) %>% 
    dplyr::filter(!is.na(gene)) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))

  dfb <- df %>% 
    dplyr::select(-cols) %>% 
    dplyr::select(-which(names(.) %in% c("prot_cover"))) %>% 
    dplyr::filter(!is.na(gene)) %>% 
    dplyr::filter(!duplicated(gene))
  
  if ("prot_cover" %in% names(df)) {
    dfc <- df %>% 
      dplyr::select(gene, prot_cover) %>% 
      dplyr::filter(!is.na(gene), !is.na(prot_cover)) %>% 
      dplyr::group_by(gene) %>% 
      dplyr::mutate(prot_cover = as.numeric(sub("%", "", prot_cover))) %>% 
      dplyr::summarise_all(~ max(.x, na.rm = TRUE)) %>% 
      dplyr::mutate(prot_cover = paste0(prot_cover, "%"))
  } else {
    dfc <- df %>% 
      dplyr::select(gene) %>% 
      dplyr::mutate(prot_cover = NA)
  }
  
  df <- list(dfc, dfb, dfa) %>% 
    purrr::reduce(right_join, by = "gene") %>% 
    dplyr::filter(!is.na(gene), !duplicated(gene))
}


#' Compare dot-dot-dot between prior and current
#' 
#' @param call_nm The name of a function call.
#' @param curr_dots The values of the current dots.
#' @param pattern The pattern for comparison.
identical_dots <- function(call_nm, curr_dots, pattern) {
  file <- file.path(dat_dir, "Calls", paste0(call_nm, ".rda"))
  if (!file.exists(file)) return(FALSE)
  
  load(file = file)
  identical(call_pars %>% .[grepl(pattern, names(.))], curr_dots)
}


#' Complete cases among sample IDs in label_scheme_sub, not label_scheme
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
my_complete_cases <- function (df, scale_log2r, label_scheme_sub) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  rows <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% label_scheme_sub$Sample_ID)) %>% 
    complete.cases(.)
  
  if (sum(rows) == 0) stop("None of the cases are complete.", call. = FALSE)
  
  df <- df[rows, ]
}


#' my union of named list 
#' 
#' names will be kept after the unification
#' 
#' @param x A list of values.
#' @param y Another list of values.
my_union <- function (x, y) {
  x %>% 
    .[! names(.) %in% names(y)] %>% 
    c(y)
}


#' find the base names of files (not currently used)
#' 
#' @param filenames A character vector of filenames
find_fn_bases <- function (filenames) {
  gsub("\\.[^.]*$", "", filenames) # %>% .[1] 
}


#' find the extensions of files
#' 
#' @param filename A character string of filename
#' @param type the type of filename extension
find_fn_exts <- function (filename, type = "text") {
  purrr::map_chr(filename, ~ {
    if (!grepl("\\.", .x)) {
      fn_ext <- switch(
        text = "txt",
        graphics = "png",
        stop("Invalid extension in file name(s).", Call. = FALSE)
      )
    } else {
      fn_ext <- gsub("^.*\\.([^.]*)$", "\\1", .x) # %>% .[1]
    }
  })
}


#' check duplicate argments in 'dots'
#' 
#' @param blacklist A character vector of variable names.
#' @param ... A list of arguments for checking.
check_dots <- function (blacklist = NULL, ...) {
  dots <- rlang::enexprs(...)
  dups <- purrr::map_lgl(names(dots), ~ .x %in% blacklist)
  nms <- names(dots[dups])
  
  if (!rlang::is_empty(nms)) {
    stop("Do not use argument(s): ", purrr::reduce(nms, paste, sep = ", "), call. = FALSE)
  }
}


#' check depreciated arguments
#' 
#' @param ... A list of arguments for checking.
#' @inheritParams check_dots
check_depreciated_args <- function (blacklist = NULL, ...) {
  dots <- rlang::enexprs(...)
  old_args <- purrr::map_chr(blacklist, `[[`, 1)
  new_args <- purrr::map_chr(blacklist, `[[`, 2)
  
  depreciated <- purrr::map_lgl(names(dots), ~ .x %in% old_args)
  nms <- names(dots[depreciated])
  
  ind <- which(old_args %in% nms)
  old_args <- old_args %>% .[ind]
  new_args <- new_args %>% .[ind]
  
  if (!rlang::is_empty(nms)) {
    message("Depreciated argument(s): ", purrr::reduce(old_args, paste, sep = ", "))
    stop("Use replacement argument(s): ", purrr::reduce(new_args, paste, sep = ", "), call. = FALSE)
  }
}


#' force 'complete_cases = TRUE' at 'impute_na = FALSE'
#' 
#' @inheritParams prnHM
to_complete_cases <- function (complete_cases = FALSE, impute_na = FALSE) {
  warn_msg1 <- "Coerce `complete_cases = TRUE` at `impute_na = FALSE`."
  if (!(impute_na || complete_cases)) {
    complete_cases <- TRUE
    rlang::warn(warn_msg1)
  }
  return(complete_cases)
}


#' Plot volcanos
#' 
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @import limma stringr purrr dplyr rlang grid gridExtra gtable
#' @importFrom magrittr %>%
plotVolcano <- function(df = NULL, df2 = NULL, id = "gene", 
                        adjP = FALSE, show_labels = TRUE, anal_type = "Volcano", 
                        gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), topn = Inf, 
                        show_sig = "none", fml_nms = NULL, gset_nms = "go_sets", 
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                        filepath = NULL, filename = NULL, theme = NULL, ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- df %>%
    `rownames<-`(.[, id]) %>%
    dplyr::select(-grep("I[0-9]{3}|log2_R[0-9]{3}|^FC", names(.)))
  
  dots <- local({
    dots <- rlang::enexprs(...)
    fmls <- dots %>% .[grepl("^\\s*~", .)]
    dots <- dots[!names(dots) %in% names(fmls)]
    dots <- concat_fml_dots(fmls, fml_nms, dots)    
  })
  
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  rm(filter_dots, arrange_dots, select_dots)
  
  if (!adjP) {
    df <- df %>%
      dplyr::select(-contains("adjP"))
  } else {
    df <- df %>%
      dplyr::select(-contains("pVal")) %>%
      `names<-`(gsub("adjP", "pVal", names(.)))
  }
  
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  fml_nms <- names(df) %>% .[grepl("pVal", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique()  %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls %>% .[names(.) %in% fml_nms]
  fml_nms <- fml_nms %>% .[map_dbl(., ~ which(.x == names(fmls)))]
  
  if (rlang::is_empty(fml_nms)) {
    stop("No formula (names) matched to those in `pepSig(...)` or `prnSig(...)`", call. = FALSE)
  }
  
  purrr::walk(file.path(filepath, fml_nms), ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))
  
  col_ind <- purrr::map(fml_nms, ~ grepl(.x, names(df))) %>%
    dplyr::bind_cols() %>%
    rowSums() %>%
    `>`(0)
  
  species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
  if (!(rlang::is_empty(species) || is.null(gset_nms))) load_dbs(gset_nms = gset_nms, species = species)

  proteoq_volcano_theme <- theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.5, size = 24),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 0, vjust = 0.5, size = 24),
      axis.title.x = element_text(colour = "black", size = 24),
      axis.title.y = element_text(colour="black", size = 24),
      plot.title = element_text(face = "bold", colour = "black", size = 14, hjust = .5, vjust = .5),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      strip.text.x = element_text(size = 16, colour = "black", angle = 0),
      strip.text.y = element_text(size = 16, colour = "black", angle = 90),
      
      legend.key = element_rect(colour = NA, fill = 'transparent'),
      legend.background = element_rect(colour = NA,  fill = "transparent"),
      legend.position = "none",
      legend.title = element_text(colour="black", size = 18),
      legend.text = element_text(colour="black", size = 18),
      legend.text.align = 0,
      legend.box = NULL
    )
  
  if (is.null(theme)) theme <- proteoq_volcano_theme

  if (length(fml_nms) > 0) purrr::pwalk(list(fml_nms, gspval_cutoff, gslogFC_cutoff, topn), 
                                        byfml_volcano, 
                                        df = df, 
                                        df2 = df2, 
                                        col_ind = col_ind, 
                                        id = !!id, 
                                        filepath = filepath, 
                                        filename = filename, 
                                        adjP = adjP, 
                                        show_labels = show_labels, 
                                        anal_type = anal_type, 
                                        show_sig = show_sig, 
                                        gset_nms = gset_nms, 
                                        scale_log2r = scale_log2r, 
                                        complete_cases = complete_cases, 
                                        impute_na = impute_na, 
                                        theme = theme, 
                                        !!!dots)
}


#' Formula specific volcano plots
#' 
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import purrr dplyr rlang
#' @importFrom magrittr %>%
byfml_volcano <- function (fml_nm, gspval_cutoff, gslogFC_cutoff, topn, df, df2, 
                           col_ind, id, 
                           filepath, filename, adjP, show_labels, anal_type, 
                           show_sig, gset_nms, scale_log2r, complete_cases, impute_na, 
                           theme = NULL, ...) {

  pval_complete_cases <- function (df) {
    rows <- df %>% 
      dplyr::select(grep("^pVal\\s+", names(.))) %>% 
      complete.cases()
    
    if (sum(rows) == 0) stop("None of the cases are complete.", call. = FALSE)
    
    df <- df[rows, ]
  }  
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- df %>%
    dplyr::select(grep(fml_nm, names(.), fixed = TRUE)) %>%
    `colnames<-`(gsub(paste0(fml_nm, "."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .) 
  
  if (complete_cases) df <- df %>% pval_complete_cases()

  byfile_plotVolcano(df = df, df2 = df2, id = !!id, fml_nm = fml_nm, filepath = filepath, filename = filename, 
                     adjP = adjP, show_labels = show_labels, anal_type = anal_type, gset_nms = gset_nms, 
                     scale_log2r = scale_log2r, 
                     impute_na = impute_na)(gspval_cutoff = gspval_cutoff, 
                                            gslogFC_cutoff = gslogFC_cutoff, 
                                            topn = topn, 
                                            show_sig = show_sig, 
                                            theme = theme, ...)
}


#' Plot volcanos
#' 
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import limma stringr purrr dplyr rlang grid gridExtra gtable
#' @importFrom magrittr %>%
byfile_plotVolcano <- function(df = NULL, df2 = NULL, id = "gene", fml_nm = NULL, filepath = NULL, filename = NULL, 
                               adjP = FALSE, show_labels = TRUE, anal_type = "Volcano", gset_nms = "go_sets", 
                               scale_log2r, impute_na, theme = NULL, ...) {

  id <- rlang::as_string(rlang::enexpr(id))

	contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
	  gsub("^log2Ratio\\s+\\(|\\)$", "", .)

	if (anal_type %in% c("Volcano")) {
		function(gspval_cutoff = 1, gslogFC_cutoff = 0, topn = Inf, show_sig = "none", theme = theme, ...) {
			rm(gspval_cutoff, gslogFC_cutoff, show_sig)

			fullVolcano(df = df, id = !!id, contrast_groups = contrast_groups,
				theme = theme, fml_nm = fml_nm, filepath = filepath, filename = filename,
				adjP = adjP, show_labels = show_labels, ...)
		}
	} else if (anal_type == "mapGSPA")
	  function(gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), topn = Inf, show_sig = "none", theme = theme, ...) {
	    stopifnot(!is.null(fml_nm))
	    
	    filepath_fml <- file.path(filepath, fml_nm)
	    
	    in_names <- list.files(path = filepath_fml, pattern = "_GSPA_[NZ]{1}.*\\.txt$")
	    if (purrr::is_empty(in_names)) stop("No inputs under ", filepath_fml, call. = FALSE)
	    
	    if (is.null(df2)) {
  	    in_names <- in_names %>% 
  	      .[!grepl("_essmap|_essmeta|_resgreedy", .)] %>% 
  	      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} %>% 
  	      {if (scale_log2r) .[grepl("_GSPA_Z", .)] else .[grepl("_GSPA_N", .)]}

  	    if (rlang::is_empty(in_names)) {
  	      stop("No inputs correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, 
  	           " at fml_nms = ", fml_nm, call. = FALSE)
  	    }
  	    
  	    df2 <- in_names
	    } else {
	      local({
	        if (grepl("_essmap|_essmeta|_resgreedy", df2)) {
	          stop("Do not use `_essmap`, `_essmeta` or `_resgreedy` for `df2`", call. = FALSE)
	        }
	        
	        non_exists <- df2 %>% .[! . %in% in_names]
	        if (!purrr::is_empty(non_exists)) {
	          stop("Missing file(s): ", purrr::reduce(non_exists, paste, sep = ", "), call. = FALSE)
	        }
	      })
	      
	      if (purrr::is_empty(df2)) stop("File(s) not found under ", filepath_fml, call. = FALSE)
	    }
	    
	    # plot data ---------------------------
	    purrr::walk(df2, gsVolcano, df = df, contrast_groups = contrast_groups, 
	                gsea_key = "term", 
	                gsets = gsets, 
	                theme = theme, 
	                fml_nm = fml_nm, 
	                filepath = filepath, filename = filename, 
	                adjP = adjP, show_labels = show_labels, show_sig = show_sig, 
	                gspval_cutoff = gspval_cutoff, gslogFC_cutoff = gslogFC_cutoff, topn = topn, ...)
	  }
}


#' Volcano plots for all proteins or peptides in a data set
#' 
#' @param contrast_groups The contrast groups defined under a formula at \code{fml_nm}.
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import dplyr purrr rlang ggplot2
#' @importFrom magrittr %>%
#' @importFrom limma vennDiagram
fullVolcano <- function(df = NULL, id = "gene", contrast_groups = NULL, theme = NULL,
                        fml_nm = NULL, filepath = NULL, filename = NULL, adjP = FALSE, 
                        show_labels = TRUE, ...) {

  id <- rlang::as_string(rlang::enexpr(id))

	dots <- rlang::enexprs(...)
	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)

	x_label <- expression("Ratio ("*log[2]*")")

	stopifnot(length(contrast_groups) > 0)

	dfw <- do.call(rbind, purrr::map(contrast_groups, ~ {
			df[, grepl(paste0(" (", .x, ")"), names(df), fixed = TRUE)] %>%
				`colnames<-`(gsub("\\s+\\(.*\\)$", "", names(.))) %>%
				mutate(Contrast = .x) %>%
				bind_cols(df[, !grepl("^pVal\\s+|^adjP\\s+|^log2Ratio\\s+", names(df)), drop = FALSE], .)
		} )) %>% 
	  dplyr::mutate(
	    Contrast = factor(Contrast, levels = contrast_groups),
	    # pVal = as.numeric(as.character(.$pVal)), 
	    valence = ifelse(.$log2Ratio > 0, "pos", "neg")
	  ) %>%
	  dplyr::filter(!is.na(pVal))
	  
	dfw_sub <- dfw %>%
		dplyr::filter(pVal < yco & abs(log2Ratio) > log2(xco)) %>%
		dplyr::arrange(Contrast, pVal) %>%
		dplyr::group_by(Contrast) %>%
		dplyr::mutate(Index = row_number()) %>%
		data.frame (check.names = FALSE)

	dfw_sub_top20 <- dfw_sub %>%
		dplyr::group_by(Contrast) %>%
		dplyr::top_n(n = -20, wt = pVal) %>%
		data.frame (check.names = FALSE)

	# tt <- gridExtra::ttheme_minimal(core = list(fg_params=list(cex = .7)), 
	#                                 colhead = list(fg_params=list(cex = .7), parse=TRUE), 
	#                                 rowhead = list(fg_params=list(cex = .7)))
	# nrow <- max(min(nrow(dfw_sub),20), 1)
	# tbl <- tableGrob(dfw_sub_top20[,c("Contrast", "Index","gene")], rows=NULL, col=NULL, theme=tt)
	# tbl$heights <- tbl$heights*.6

	# data table for labels
	dt <- purrr::map(contrast_groups, ~ {
	  to_csv_(dfw_sub_top20 %>%
            dplyr::filter(Contrast == .x) %>%
            dplyr::select(c("Index", id))) %>%
            {if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .} # a zero-entry exception
	  })  %>%
	  do.call(rbind, .) %>%
	  data.frame(Contrast = contrast_groups, id = ., stringsAsFactors = FALSE) %>%
	  dplyr::rename(!!rlang::sym(id) := id) %>%
	  dplyr::mutate(Contrast = factor(Contrast,  levels = contrast_groups))

	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	myPalette <- c("#377EB8", "#E41A1C")
	xmax <- ceiling(pmax(abs(min(dfw$log2Ratio, na.rm = TRUE)), max(dfw$log2Ratio, na.rm = TRUE)))
	ymax <- ceiling(max(-log10(dfw$pVal), na.rm = TRUE)) * 1.1
	dfw_greater <- dfw_sub[dfw_sub$valence == "pos", ]
	dfw_less <- dfw_sub[dfw_sub$valence == "neg", ]

	if(is.null(dots$nrow)) {
		nrow <- ifelse(length(unique(dfw$Contrast)) > 3, 2, 1)
	} else {
		nrow <- dots$nrow
	}

	p <-ggplot() +
		geom_point(data = dfw, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, colour = "gray", shape = 20, alpha = .5) +
		geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, color = myPalette[2], shape = 20, alpha = .8) +
		geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, color = myPalette[1], shape = 20, alpha = .8) +
		geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = log2(xco), linetype = "longdash", size = .5) +
		scale_x_continuous(limits = c(-xmax, xmax)) +
		scale_y_continuous(limits = c(0, ymax)) +
		labs(title = "", x = x_label, y = expression("P-value ("*-log[10]*")")) +
		theme

	if (show_labels) {
		p <- p + geom_text(data = dfw_sub_top20, 
		                   mapping = aes(x = log2Ratio, y = -log10(pVal), label = Index, color = Index),
		                   size = 3, alpha = .5, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05, 
		                   na.rm = TRUE)
		p <- p + facet_wrap(~ Contrast, nrow = nrow, labeller = label_value)
		p <- p + geom_table(data = dt, aes(table = !!rlang::sym(id)), x = -xmax*.85, y = ymax/2)
	} else {
		p <- p + facet_wrap(~ Contrast, nrow = nrow, labeller = label_value)
	}

	if(is.null(dots$width)) {
		width <- ifelse(nrow > 1, 6*length(unique(dfw$Contrast))/nrow + 1, 
		                6*length(unique(dfw$Contrast)) / nrow)
	} else {
		width <- dots$width
	}

	if(is.null(dots$height)) {
		height <- 6*nrow
	} else {
		height <- dots$height
	}

	ggsave(file.path(filepath, fml_nm, filename), p, width = width, height = height, units = "in")

	# Venn
	summ_venn <- function(df, id, contrast_groups) {
		stopifnot(length(contrast_groups) <= 5)

		universe <- sort(unique(df[[id]]))

		Counts <- matrix(0, nrow = length(universe), ncol = length(contrast_groups)) %>%
			`rownames<-`(universe) %>%
			`colnames<-`(contrast_groups)

		for (Group in contrast_groups) {
			ids <- df %>% dplyr::filter(Contrast == Group) %>% dplyr::select(id) %>% unlist()
			Counts[, Group] <- rownames(Counts) %in% ids
		}

		return(Counts)
	}

	plot_venn <- function(Counts, filepath, direction) {
		myPalette <- brewer.pal(n = 9, name = "Set1")
		Width <- 3
		Height <- 3

		fn_prefix <- paste0(direction, "_venn")

		png(file.path(filepath, fml_nm, paste0(fn_prefix, ".png")), width = Width, height = Height, units = "in", res = 300)
			limma::vennDiagram(limma::vennCounts(Counts), circle.col = myPalette, cex = .5)
		dev.off()
		write.table(Counts, file.path(filepath, fml_nm, paste0(fn_prefix, ".txt")), sep = "\t", 
		            col.names = TRUE, row.names = TRUE, quote = FALSE)	
	}

	summ_venn(dfw_greater, id, contrast_groups) %>% plot_venn(filepath, paste0(fn_prefix, "_greater"))
	summ_venn(dfw_less, id, contrast_groups) %>% plot_venn(filepath, paste0(fn_prefix, "_less"))
}


#' Volcano plots of protein \code{log2FC} under given gene sets
#'
#' @param gsea_key Character string; the column key indicating the terms of gene sets.
#' @param gsets The gene sets.
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @inheritParams fullVolcano
#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
gsVolcano <- function(df2 = NULL, df = NULL, contrast_groups = NULL, 
                      gsea_key = "term", gsets = NULL, 
                      theme = NULL, 
                      fml_nm = NULL, 
                      filepath = NULL, filename = NULL, adjP = FALSE, show_labels = TRUE, show_sig = "none", 
                      gspval_cutoff = 1E-6, gslogFC_cutoff = log2(1.2), topn = Inf, ...) {
  
  custom_prefix <- gsub("(.*_{0,1})Protein_GSPA.*", "\\1", df2)
  dir.create(path = file.path(filepath, fml_nm, custom_prefix), recursive = TRUE, showWarnings = FALSE)

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) 
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  filename <- paste0(custom_prefix, fn_prefix, ".", fn_suffix)
  
  dots <- rlang::enexprs(...)

	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)
	x_label <- expression("Ratio ("*log[2]*")")
	y_label <- ifelse(adjP, expression("pVal ("*-log[10]*")"), expression("adjP ("*-log[10]*")"))
	
	gsea_res <- tryCatch(readr::read_tsv(file.path(filepath, fml_nm, df2), 
	                                          col_types = cols(term = col_character(),
	                                                           is_essential = col_logical(),
	                                                           size = col_double(),
	                                                           ess_size = col_double(),
	                                                           contrast = col_character(),
	                                                           p_val = col_double(),
	                                                           q_val = col_double(),
	                                                           log2fc = col_double())), 
	                          error = function(e) NA)
	
	message("Secondary file loaded: ", gsub("\\\\", "/", file.path(filepath, fml_nm, df2)))

	filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
	arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
	select2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select2_", names(.))]
	dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots, select2_dots)]
	
	gsea_res <- gsea_res %>% 
	  dplyr::arrange(p_val, abs(log2fc)) %>% 
	  filters_in_call(!!!filter2_dots) %>% 
	  arrangers_in_call(!!!arrange2_dots)
	rm(filter2_dots, arrange2_dots, select2_dots)

	if (nrow(gsea_res) == 0) stop("No GSPA terms available after data filtration.")

	topn <- pmin(dplyr::n_distinct(gsea_res[[gsea_key]]), topn)
	terms <- gsea_res %>%
	  dplyr::arrange(p_val) %>% 
	  dplyr::slice(1:(topn*length(contrast_groups))) %>% 
	  dplyr::filter(p_val <= gspval_cutoff, abs(log2fc) >= gslogFC_cutoff) %>%
	  dplyr::select(gsea_key) %>%
	  unique() %>%
	  unlist() %>% 
	  .[. %in% names(gsets)]
	
	gsea_res <- gsea_res %>% 
	  dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>% 
	  dplyr::mutate(q_val = format(p_val, scientific = TRUE, digits = 2)) %>% 
	  dplyr::mutate(log2fc = round(log2fc, digits = 2))

	if (length(terms) > 0) {
  	dfw <- do.call(rbind,
  		purrr::map(contrast_groups, ~ {
  			df[, grepl(paste0(" (", .x, ")"), names(df), fixed = TRUE)] %>%
  				`colnames<-`(gsub("\\s+\\(.*\\)$", "", names(.))) %>%
  				mutate(Contrast = .x) %>%
  				bind_cols(df[, !grepl("^pVal\\s+|^adjP\\s+|^log2Ratio\\s+", names(df))], .)
  		} )) %>%
  		dplyr::mutate(Contrast = factor(Contrast, levels = contrast_groups),
  			pVal = as.numeric(pVal),
  			valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>%
  		dplyr::filter(!is.na(pVal))

  	lapply(terms, function(gt) {
  	  # some results may be based on gene sets from an older version, 
  	  # which are no longer available in the current version
  	  
  	  gsets_sub <- gsets %>% .[names(.) == gt]
  		if (length(gsets_sub) == 0) return(NULL)

  		fn <- gsub(":", "~", gsub("/", "or", names(gsets_sub)[[1]]), fixed = TRUE)

  		res_sub <- gsea_res[as.character(gsea_res$term) == gt, ] %>% 
  		  data.frame(check.names = FALSE)

  		dfw_sub <- dfw[as.character(dfw$entrez) %in% gsets_sub[[1]], ]

  		if (nrow(dfw_sub) == 0) return(NULL) 

  		xmax <- ceiling(pmax(abs(min(dfw_sub$log2Ratio)), max(dfw_sub$log2Ratio)))
  		ymax <- ceiling(max(-log10(dfw_sub$pVal))) * 1.1

  		# ensure the same levels between "Levels" and "newLevels"
  		Levels <- levels(dfw_sub$Contrast)

  		dfw_sub <- dfw_sub %>%
  			dplyr::arrange(Contrast, pVal) %>%
  			dplyr::group_by(Contrast) %>%
  			dplyr::mutate(Index = row_number()) %>%
  			data.frame(check.names = FALSE) %>% 
  		  dplyr::mutate(Contrast = as.character(Contrast)) %>% 
  			dplyr::left_join(., res_sub, by = c("Contrast" = "contrast")) %>%
  			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
  			dplyr::arrange(Contrast) %>%
  			# dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>%
  			# dplyr::mutate(p_val = as.numeric(p_val)) %>%
  			# dplyr::mutate(q_val = format(q_val, scientific = TRUE, digits = 2)) %>%
  			# dplyr::mutate(q_val = as.numeric(q_val)) %>%
  			# mutate(sig_level = ifelse(.$q.val > 0.05, "n.s.", ifelse(.$q.val > 0.005, "*", "**"))) %>%
  			# mutate(newContrast = paste0(Contrast, " (", sig_level, ")"))
  			dplyr::mutate(newContrast = Contrast)
				
  		if (show_sig != "none") {
  			if (grepl("^p", show_sig)) {
  				dfw_sub <- dfw_sub %>%
  					dplyr::mutate(newContrast = paste0(Contrast, " (p = ", p_val, ")"))
  			} else if (grepl("^q", show_sig)) {
  				dfw_sub <- dfw_sub %>%
  					dplyr::mutate(newContrast = paste0(Contrast, " (q = ", q_val, ")"))
  			}
  		}

  		newLevels <- unique(dfw_sub$newContrast)

  		dfw_sub <- dfw_sub %>%
  			dplyr::mutate(newContrast = factor(newContrast, levels = newLevels)) %>%
  			dplyr::arrange(newContrast)

  		dfw_sub_top20 <- dfw_sub %>%
  			dplyr::group_by(newContrast) %>%
  			dplyr::top_n(n = -20, wt = pVal)

  		dt <- purrr::map(newLevels, ~ {
  				dfw_sub_top20 %>%
  					dplyr::filter(newContrast == .x) %>%
  					data.frame(check.names = FALSE) %>%
  					dplyr::select(c("Index", "gene")) %>%
  					to_csv_() %>%
  					{if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .}
  			}) %>%
  			do.call(rbind, .) %>%
  			data.frame(newContrast = newLevels, Contrast = Levels, Gene = ., stringsAsFactors = FALSE) %>%
  			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
  			dplyr::mutate(newContrast = factor(newContrast, levels = newLevels))

  		dfw_greater <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio > log2(xco))
  		dfw_less <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio < -log2(xco))

  		dt_pos <- ifelse(nrow(dfw_greater) > nrow(dfw_less), -xmax*.85, xmax*.6) # table position
  		myPalette <- c("#377EB8", "#E41A1C")
  		
  		if(is.null(dots$nrow)) {
  		  nrow <- ifelse(length(unique(dfw_sub$Contrast)) > 3, 2, 1)
  		} else {
  		  nrow <- dots$nrow
  		}

  		p <- ggplot() +
  			geom_point(data = dfw_sub, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, colour = "gray", shape = 20, alpha = .5) +
  			geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[2], shape = 20, alpha = .8) +
  			geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[1], shape = 20, alpha = .8) +
  			geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = log2(xco), linetype = "longdash", size =.5) +
  			labs(title = names(gsets_sub), x = x_label, y = y_label) +
  			scale_x_continuous(limits = c(-xmax, xmax)) +
  			scale_y_continuous(limits = c(0, ymax)) +
  			theme
  		p <- p + facet_wrap(~ newContrast, nrow = nrow, labeller = label_value)

  		if(show_labels)
  			p <- p + geom_text(data = dfw_sub_top20, mapping = aes(x = log2Ratio, y = -log10(pVal),
  					label = Index, color = Index), size = 2, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05) +
  				geom_table(data = dt, aes(table = Gene), x = dt_pos, y = ymax/2)

  		if (nchar(fn) > 50) fn <- paste0(str_sub(fn, 1, 50), "...") # avoid long fns for pdf()
  		
  		if (is.null(dots$width)) {
  		  width <- ifelse(nrow > 1, 6*length(unique(dfw_sub$newContrast))/nrow + 1, 
  		                  6*length(unique(dfw$Contrast)) / nrow)
  		} else {
  		  width <- dots$width
  		}
  		
  		if (is.null(dots$height)) {
  		  height <- 6*nrow
  		} else {
  		  height <- dots$height
  		}
  		
  		ggsave(file.path(filepath, fml_nm, custom_prefix, paste0(fn, ".", fn_suffix)), 
  		       p, width = width, height = height, dpi = 300, units = "in")

  		write.table(dfw_sub, file = file.path(filepath, fml_nm, custom_prefix, paste0(fn, ".txt")), 
  		            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)	
  		            
  	})
	}
}



#'Volcano plots
#'
#'\code{pepVol} visualizes the volcano plots of peptide data.
#'
#'@rdname prnVol
#'
#'@import purrr
#'@export
pepVol <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjP = FALSE, show_labels = TRUE, 
                    df = NULL, filepath = NULL, filename = NULL, 
                    fml_nms = NULL, theme = NULL, ...) {
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"))
  
  scale_log2r <- match_pepSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)	
  
  stopifnot(rlang::is_logical(adjP), 
            rlang::is_logical(show_labels))

  reload_expts()

  info_anal(df = !!df, df2 = NULL, id = !!id, filepath = !!filepath, filename = !!filename, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            anal_type = "Volcano")(fml_nms = fml_nms, 
                                   adjP = adjP, 
                                   show_labels = show_labels, 
                                   theme = theme, 
                                   ...)
}


#'Volcano plots
#'
#'\code{prnVol} visualizes the volcano plots of protein data.
#'
#'@inheritParams prnGSPA
#'@inheritParams prnHist
#'@param adjP Logical; if TRUE, use Benjamini-Hochberg pVals in volcano plots.
#'  The default is FALSE.
#'@param show_labels Logical; if TRUE, shows the labels of top twenty entries.
#'  The default is TRUE.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  Additional parameters for plotting: \cr \code{xco}, the cut-off lines of
#'  fold changes at position \code{x}; the default is at \eqn{-1.2} and
#'  \eqn{+1.2}. \cr \code{yco}, the cut-off line of \code{pVal} at position
#'  \code{y}; the default is \eqn{0.05}. \cr \code{width}, the width of plot;
#'  \cr \code{height}, the height of plot. \cr \code{nrow}, the number of rows
#'  in a plot.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnVol_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
prnVol <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjP = FALSE, show_labels = TRUE, 
                    df = NULL, filepath = NULL, filename = NULL, 
                    fml_nms = NULL, theme = NULL, ...) {
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)	
  
  stopifnot(rlang::is_logical(adjP), 
            rlang::is_logical(show_labels))

  reload_expts()

  info_anal(df = !!df, df2 = NULL, id = !!id, filepath = !!filepath, filename = !!filename, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            anal_type = "Volcano")(fml_nms = fml_nms, 
                                   adjP = adjP, 
                                   show_labels = show_labels, 
                                   theme = theme, 
                                   ...)
}


#'Volcano plots of protein \code{log2FC} under gene sets
#'
#'\code{gspaMap} visualizes the volcano plots of protein subgroups under the
#'same gene sets.
#'
#'@inheritParams prnGSPA
#'@inheritParams prnVol
#'@inheritParams prnHist
#'@inheritParams plot_prnTrend
#'@inheritParams anal_prnTrend
#'@param filename Use system default for each gene set.
#'@param show_sig Character string indicating the type of significance values to
#'  be shown with \code{\link{gspaMap}}. The default is \code{"none"}.
#'  Additional choices are from \code{c("pVal", "qVal")} where \code{pVal} or
#'  \code{qVal} will be shown, respectively, in the facet grid of the plots.
#'@param gspval_cutoff Numeric value or vector for uses with
#'  \code{\link{gspaMap}}. \code{Gene sets} with enrichment \code{pVals} less
#'  significant than the threshold will be excluded from volcano plot
#'  visualization. The default significance is 0.05 for all formulas matched to
#'  or specified in argument \code{fml_nms}. Formula-specific threshold is
#'  allowed by supplying a vector of cut-off values.
#'@param gslogFC_cutoff Numeric value or vector for uses with
#'  \code{\link{gspaMap}}. \code{Gene sets} with absolute enrichment
#'  \code{log2FC} less than the threshold will be excluded from volcano plot
#'  visualization. The default magnitude is \code{log2(1.2) } for all formulas
#'  matched to or specified in argument \code{fml_nms}. Formula-specific
#'  threshold is allowed by supplying a vector of absolute values in
#'  \code{log2FC}.
#'@param topn Numeric value or vector; top entries in gene sets ordered by
#'  increasing \code{pVal} for visualization. The default is to use all
#'  available entries.
#'@param gset_nms Character string or vector containing the shorthanded name(s),
#'  full file path(s) or both to gene sets for enrichment analysis. For species
#'  among \code{"human", "mouse", "rat"}, the default of \code{c("go_sets",
#'  "c2_msig")} will utilize terms from both gene ontology (\code{GO}) and
#'  molecular signatures (\code{MSig}). Custom data bases of \code{GO} and
#'  curated \code{MSig}, and/or additional species are also supported. See also
#'  \code{\link{prepGO}} for the preparation of custom \code{GO} and
#'  \code{\link{prepMSig}} for the preparation of custom \code{MSig}.
#'
#'  Note that it is users' responsibility to ensure that the custom gene sets
#'  contain terms that can be found from the one or multiple preceding analyses
#'  of \code{\link{prnGSPA}}. For simplicity, it is generally applicable to
#'  include \emph{all} of the data bases that have been applied to
#'  \code{\link{prnGSPA}} and in that way no terms will be missed for
#'  visualization. See also \code{\link{prnGSPA}} for examples of custom data
#'  bases.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements and the
#'  association between \code{filter_} and \code{df}. \cr \cr \code{filter2_}:
#'  Variable argument statements for the row filtration against data in
#'  secondary file(s) linked to \code{df2}. See also \code{\link{prnGSPAHM}} for
#'  the format of \code{filter2_}, \code{normPSM} for the association between
#'  \code{filter_} and \code{df}. \cr \cr Additional parameters for plotting:
#'  \cr \code{xco}, the cut-off lines of fold changes at position \code{x}; the
#'  default is at \eqn{-1.2} and \eqn{+1.2}. \cr \code{yco}, the cut-off line of
#'  \code{pVal} at position \code{y}; the default is \eqn{0.05}. \cr
#'  \code{width}, the width of plot; \cr \code{height}, the height of plot. \cr
#'  \code{nrow}, the number of rows in a plot.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/prnVol_.R
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
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{dl_stringdbs}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  # Mascot \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'  
#'  # MaxQuant \cr
#'  system.file("extdata", "maxquant_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "maxquant_protein_keys.txt", package = "proteoQ") \cr
#'
#'@export
gspaMap <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     df = NULL, df2 = NULL, filepath = NULL, filename = NULL, fml_nms = NULL, 
                     adjP = FALSE, show_labels = TRUE, show_sig = "none", 
                     gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), topn = Inf, 
                     gset_nms = c("go_sets", "c2_msig"), theme = NULL, ...) {
  check_dots(c("id", "anal_type"), ...)
  check_depreciated_args(list(c("pval_cutoff", "gspval_cutoff"), c("logFC_cutoff", "gslogFC_cutoff")), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  df2 <- rlang::enexpr(df2)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)	
  show_sig <- rlang::as_string(rlang::enexpr(show_sig))
  
  stopifnot(rlang::is_logical(adjP), 
            rlang::is_logical(show_labels), 
            rlang::is_double(gspval_cutoff), 
            rlang::is_double(gslogFC_cutoff))
  
  reload_expts()

  info_anal(df = !!df, df2 = !!df2, id = !!id, filepath = !!filepath, filename = !!filename, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            anal_type = "mapGSPA")(fml_nms = fml_nms, 
                                   adjP = adjP, 
                                   show_labels = show_labels, 
                                   gspval_cutoff = gspval_cutoff, 
                                   gslogFC_cutoff = gslogFC_cutoff, 
                                   topn = topn, 
                                   show_sig = show_sig,
                                   gset_nms = gset_nms, 
                                   theme = theme, 
                                   ...)
}



GeomTable <- ggproto(
	"GeomTable",
	Geom,
	required_aes = c("x", "y",  "table"),
	default_aes = aes(widthx = 10, widthy = 10, rownames = NA),
	draw_key = draw_key_blank,

	draw_panel = function(data, panel_scales, coord) {
  	if (nrow(data) != 1) {
  		stop(sprintf("only one table per panel allowed, got %s (%s)", 
  		             nrow(data), 
  		             as.character(data)), 
  		     call. = FALSE)
  	}
	  
	  wx = data$widthx / 2
	  wy = data$widthy / 2

  	corners <- data.frame(x = c(data$x - wx, data$x + wx), y = c(data$y - wy, data$y + wy))
  	d <- coord$transform(corners, panel_scales)
  
  	table = read.csv(text = data$table, header = TRUE)
  	if (!is.na(data$rownames)) {
  		rownames(table) <-
  		unlist(strsplit(data$rownames, "|", fixed = TRUE))
  	}
  
  	x_rng <- range(d$x, na.rm = TRUE)
  	y_rng <- range(d$y, na.rm = TRUE)

  	vp <- viewport(x = mean(x_rng), y = mean(y_rng), width = diff(x_rng), height = diff(y_rng), 
  	               just = c("center", "center"))
  	  
  	grob <- tableGrob(table, rows = NULL, col = NULL,
  	                  theme = ttheme_minimal(core = list(fg_params=list(cex = .7)),
  	                                         colhead = list(fg_params=list(cex = .7), parse=TRUE),
  	                                         rowhead = list(fg_params=list(cex = .7))))
  	grob$heights <- grob$heights*.6
  
  	## add a line across the header
  	# grob <- gtable_add_grob(
  	#   grob,
  	#   grobs = segmentsGrob(y1 = unit(0, "npc"), gp = gpar(lwd = 2.0)),
  	#   t = 1,
  	#   b = 1,
  	#   l = 1,
  	#   r = ncol(d) + 1
  	# )
  	editGrob(grob, vp = vp, name = paste(grob$name, facet_id()))
	}
)


facet_id <- local({
	i <- 1
	function() {
  	i <<- i + 1
  	i
	}
})

#' Print table in ggplot2 images
#' 
#' @param mapping The same as ggplot2.
#' @param data same as ggplot2.
#' @param stat The same as ggplot2.
#' @param position same as ggplot2.
#' @param na.rm The same as ggplot2.
#' @param show.legend same as ggplot2.
#' @param inherit.aes The same as ggplot2.
#' @param ... same as ggplot2.
geom_table <- function(mapping = NULL, data = NULL, stat = "identity",
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, ...) {
  layer(geom = GeomTable, mapping = mapping, data = data, stat = stat, position = position,
        show.legend = show.legend, inherit.aes = inherit.aes, params = list(na.rm = na.rm, ...)
  )
}


#' Convert a data frame column to a csv vector
#' 
#' @param x A data frame column.
to_csv_ <- function(x) {
	paste(capture.output(write.csv(x, stdout(), row.names = F)), collapse = "\n")
}


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!")
}
#' proteoQ: A package for processing mass spectrometry data using tandem mass
#' tags (\url{https://en.wikipedia.org/wiki/Tandem_mass_tag}).
#'
#' The proteoQ package provides three categories of functions in data
#' preprocessing, quality assurance assessments and informatic analysis.
#'
#' @section proteoQ functions in data preprocessing: normPSM -> (purgePSM) ->
#'   PSM2Pep -> mergePep -> standPep -> (purgePep) -> standPrn
#' @section proteoQ functions in data QA: pepHist, prnHist, pepCorr_logFC,
#'   pepCorr_logInt, prnCorr_logFC, prnCorr_logInt,
#' @section proteoQ functions in informatic analysis: pepMDS, prnMDS, pepPCA,
#'   prnPCA, pepEucDist, prnEucDist, pepHM, prnHM, anal_pepNMF, anal_prnNMF,
#'   plot_pepNMFCon, plot_prnNMFCon, plot_pepNMFCoef, plot_prnNMFCoef,
#'   plot_metaNMF, anal_prnTrend, plot_prnTrend, pepSig, prnSig, pepVol, prnVol,
#'   prnGSPA, prnGSPAHM, gspaMap, anal_prnString, pepImp, prnImp, prnGSVA,
#'   prnGSEA
#' @section proteoQ functions in custom database preparation: prepGO, prepMSig, 
#' Uni2Entrez, Ref2Entrez
#'
#' @docType package
#' @name proteoQ
NULL
#> NULL
