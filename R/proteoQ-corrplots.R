#' Correlation plots
#'
#' @import stringr dplyr ggplot2 GGally rlang
#' @importFrom magrittr %>%
plotCorr <- function (df = NULL, id, anal_type, data_select, col_select = NULL, col_order = NULL,
                      label_scheme_sub = label_scheme_sub, 
                      scale_log2r = scale_log2r, 
                      filepath = filepath, filename = filename, ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  xmin <- eval(dots$xmin, env = caller_env()) # `xmin = -1` is `language`
  xmax <- eval(dots$xmax, env = caller_env()) # `xmax = +1` is `language`
  x_breaks <- eval(dots$x_breaks, env = caller_env())
  width <- eval(dots$width, env = caller_env())
  height <- eval(dots$height, env = caller_env())

  nm_idx <- names(dots) %in% c("xmin", "xmax", "x_breaks", "width", "height")
                               
  dots[nm_idx] <- NULL
  
  lang_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% lang_dots]
  
  # cat("Available column keys for data filtration: \n")
  # cat(paste0(names(df), "\n"))
  
  df <- df %>% filters_in_call(!!!lang_dots)

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
	  if (is.null(x_breaks)) x_breaks <- 1
	} else if (data_select == "logInt") {
	  df <- df %>% .$Intensity %>% log10()
	  y_label <- x_label <- expression("Intensity ("*log[10]*")")
	  if (is.null(xmin)) xmin <- 3.5
	  if (is.null(xmax)) xmax <- 6
	  if (is.null(x_breaks)) x_breaks <- 1
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
	              xmin = xmin, xmax = xmax, x_breaks = x_breaks, width = width, height = height, !!!dots)
}


#' Make correlation plots
#'
#' @import stringr dplyr ggplot2 GGally purrr rlang
#' @importFrom magrittr %>%
plot_corr_sub <- function (df, xlab, ylab, filename, filepath, 
                           xmin, xmax, x_breaks, width, height, ...) {
                           
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

  ggsave(file.path(filepath, gg_imgname(filename)), plot = p1, width = width, height = height, dpi = dpi, units = "in")
         
}


#'Correlation Plots
#'
#'\code{proteoCorr} plots Pearson correlation for both \code{logFC} and
#'\code{intensities} data. Users should avoid call the method directly, but
#'instead use the following wrappers.
#'
#'The function matches the current \code{id} to the grouping argument in the
#'latest \code{call} to \code{\link{normPSM}} or \code{\link{normPep}}.  For
#'example, if \code{normPSM(group_psm_by = pep_seq, ...)} was called earlier,
#'the setting of \code{id = pep_seq_mod} in the current call will be matched to
#'\code{id = pep_seq}. Similarly, if \code{normPep(group_pep_by = gene, ...)}
#'was employed, the setting of \code{id = prot_acc} in the current call will be
#'matched to \code{id = gene}.
#'
#'@inheritParams proteoHist
#'@param  col_order Character string to a column key in \code{expt_smry.xlsx}.
#'  Numeric values under which will be used for the left-to-right arrangement of
#'  samples in plots. At the NULL default, the column key \code{Order} will be
#'  used.
#'@param data_select The subset of data to be selected. At default, \code{logFC}
#'  will be used; at \code{logInt}, intensity with \code{log10} transformation
#'  will be used.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr Additional parameters for
#'  plotting: \cr \code{width}, the width of plot; \cr \code{height}, the height
#'  of plot; \cr \code{xmin}, the minimum \eqn{x} of logFC or intensity. \cr
#'  \code{xmax}, the maximum \eqn{x} of logFC data or intensity data.
#'@return Correlation plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoCorr <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                        data_select = c("logFC", "logInt"), 
                        col_select = NULL, col_order = NULL, scale_log2r = TRUE, 
                        df = NULL, filepath = NULL, filename = NULL, ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  id <- rlang::enexpr(id)
	if (length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))

	col_select <- rlang::enexpr(col_select)
	col_order <- rlang::enexpr(col_order)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)

	reload_expts()

	info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
	          scale_log2r = scale_log2r, impute_na = FALSE,
						df = !!df, filepath = !!filepath, filename = !!filename,
						anal_type = "Corrplot")(data_select, ...)
						                        
						                        
}


#'Correlation plots of the \code{log2FC} of peptide data
#'
#'\code{pepCorr_logFC} is a wrapper of \code{\link{proteoCorr}} for peptide logFC
#'data.
#'
#'@rdname proteoCorr
#'
#' @examples
#' pepCorr_logFC(
#'  width = 10,
#'  height = 10,
#'  filter_by = exprs(n_psm >= 3),
#'  filename = pepcorr_logfc_npsm3.png,
#' )
#'
#'@import purrr
#'@export
pepCorr_logFC <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\Corrplot\\log"), recursive = TRUE, showWarnings = FALSE)

  quietly_log <- purrr::quietly(proteoCorr)(id = pep_seq, data_select = "logFC", ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\Corrplot\\log","pepCorr_log.csv"), append = TRUE)
}


#'Correlation plots of the \code{log10} intensity of peptide data
#'
#'\code{pepCorr_logInt} is a wrapper of \code{\link{proteoCorr}} for peptide
#'intensity data.
#'
#'@rdname proteoCorr
#'
#' @examples
#' pepCorr_logInt(
#'  width = 10,
#'  height = 10,
#'  filter_by = exprs(n_psm >= 3),
#'  filename = pepcorr_int_npsm3.png,
#' )
#'
#'@import purrr
#'@export
pepCorr_logInt <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\Corrplot\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoCorr)(id = pep_seq, data_select = "logInt", ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\Corrplot\\log","pepCorr_log.csv"), append = TRUE)
}



#'Correlation plots of the \code{log2FC} of protein data
#'
#'\code{prnCorr_logFC} is a wrapper of \code{\link{proteoCorr}} for protein
#'logFC data.
#'
#'@rdname proteoCorr
#'
#' @examples
#' prnCorr_logFC(
#'  width = 10,
#'  height = 10,
#'  filter_npep = exprs(n_pep >= 5),
#'  filename = prncorr_logfc_npep5.png,
#' )
#'
#'@import purrr
#'@export
prnCorr_logFC <- function (...) {
  err_msg <- "Don't call the function with arguments `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Corrplot\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoCorr)(id = gene, data_select = "logFC", ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Corrplot\\log","prnCorr_log.csv"), append = TRUE)
}


#'Correlation plots of the \code{log10} intensity of protein data
#'
#'\code{prnCorr_logInt} is a wrapper of \code{\link{proteoCorr}} for protein
#'intensity data.
#'
#'@rdname proteoCorr
#'
#' @examples
#' prnCorr_logInt(
#'  width = 10,
#'  height = 10,
#'  filter_npep = exprs(n_pep >= 5),
#'  filename = prncorr_int_npep5.png,
#' )
#'
#'@import purrr
#'@export
prnCorr_logInt <- function (...) {
  err_msg <- "Don't call the function with arguments `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Corrplot\\log"), recursive = TRUE, showWarnings = FALSE)

  quietly_log <- purrr::quietly(proteoCorr)(id = gene, data_select = "logInt", ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Corrplot\\log","prnCorr_log.csv"), append = TRUE)
}


