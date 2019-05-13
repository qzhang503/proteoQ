#' Correlation plots
#'
#' @import stringr dplyr ggplot2 GGally rlang
#' @importFrom magrittr %>%
plotCorr <- function (df = NULL, col_select = NULL, col_order = NULL,
                      label_scheme_sub = label_scheme_sub, use_log10 = TRUE,
                      scale_log2r = scale_log2r, min_int = min_int, max_int = max_int,
                      min_log2r = min_log2r, max_log2r = max_log2r,
                      filepath = filepath, filename = filename, ...) {

	dots <- rlang::enexprs(...)

	col_select <- rlang::enexpr(col_select)
	col_order <- rlang::enexpr(col_order)

	fn_prx <- gsub("\\..*$", "", filename)

	if (use_log10) df$Intensity <- log10(df$Intensity)

	if(dplyr::n_distinct(label_scheme_sub$Order) == 1) {
		df$Intensity <- df$Intensity[, order(names(df$Intensity))]
		df$log2R <- df$log2R[, order(names(df$log2R))]
	} else {
		corrplot_orders <- label_scheme_sub %>%
			dplyr::select(!!col_select, !!col_order) %>%
			dplyr::filter(!is.na(!!col_order)) %>%
			unique(.) %>%
			dplyr::arrange(!!col_order)

		df$Intensity <- df$Intensity[, as.character(corrplot_orders[[col_select]]), drop = FALSE]
		df$log2R <- df$log2R[, as.character(corrplot_orders[[col_select]]), drop = FALSE]
	}

	y_label <- x_label <- expression("Ratio ("*log[2]*")")
	y_label_int <- x_label_int <- ifelse(use_log10,
	                                     expression("Intensity ("*log[10]*")"), "Intensity")

	plot_corr_sub(df = df$Intensity, xlab = x_label_int, ylab = y_label_int,
	              filename = paste0(fn_prx, "_Intensity.png"), filepath = filepath,
	              xmin = min_int, xmax = max_int, ...)

	plot_corr_sub(df = df$log2R, xlab = x_label, ylab = y_label,
	              filename = paste0(fn_prx, "_log2Ratio.png"), filepath = filepath,
	              xmin = min_log2r, xmax = max_log2r, ...)
}


#' Make correlation plots
#'
#' @import stringr dplyr ggplot2 GGally purrr rlang
#' @importFrom magrittr %>%
plot_corr_sub <- function (df, xlab, ylab, filename, filepath, xmin, xmax, ...) {

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

  fn_prx <- gsub("\\..*$", "", filename)

  n_TMT_sets <- n_TMT_sets(label_scheme)

  if(is.null(dots$width)) {
    width <- 10*n_TMT_sets*1.2
  } else {
    width <- eval(dots$width, env = caller_env())
    dots$width <- NULL
  }

  if(is.null(dots$height)) {
    height <- 10*n_TMT_sets*1.2
  } else {
    height <- eval(dots$height, env = caller_env())
    dots$height <- NULL
  }

  if(is.null(dots$dpi)) {
    dpi <- 300
  } else {
    dpi <- eval(dots$dpi, env = caller_env())
    dots$dpi <- NULL
  }

  png(file.path(filepath, filename), width = width, height = height, units = "in", res = 300)
    pairs(df, pch = ".", upper.panel = panel_cor, lower.panel = panel.smooth)
  dev.off()

  ncol <- ncol(df)

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

  ggsave(file.path(filepath, paste0(fn_prx, "_gg.png")), p1, width = width, height = height,
         dpi = dpi, units = "in")
}


#'Correlation Plots
#'
#'\code{proteoCorrplot} produces correlation plots for both \code{log2-ratios}
#'and reporter-ion \code{intensities}.
#'
#'The function matches the current \code{id} to those in the latest \code{calls}
#'to \code{\link{normPep}} or \code{\link{normPrn}}.  For example, if
#'\code{pep_seq} was used in \code{\link{normPep()}}, the current \code{id =
#'pep_seq_mod} will be matched to \code{id = pep_seq}.
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at \code{id = pep_seq} or \code{pep_seq_mod}, and protein data at
#'  \code{id = prot_acc} or \code{gene}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param  col_order Not currently used.
#'@param use_log10 Logical; if TRUE, \code{log10} transformation of
#'  \code{intensity}.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param min_int The minimum intensity at \code{log10} scale.
#'@param max_int The maximum intensity at \code{log10} scale.
#'@param min_log2r The minimum \code{log2-ratio} for display.
#'@param max_log2r The maximum \code{log2-ratio} for display.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the names of the current functional \code{call}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... Parameters inherited from \code{ggsave}: \cr \code{width}, the
#'  width of plot; \cr \code{height}, the height of plot \cr \code{...}
#'@return Correlation plots.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoCorrplot <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), col_select = NULL,
                            col_order = NULL, use_log10 = TRUE, scale_log2r = TRUE,
														min_int = 3.5, max_int = 6.5, min_log2r = -2, max_log2r = 2,
														df = NULL, filepath = NULL, filename = NULL, ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

  id <- rlang::enexpr(id)
	if(length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))

	col_select <- rlang::enexpr(col_select)
	col_order <- rlang::enexpr(col_order)
	
	reload_expts()

	info_anal(id = !!id, col_select = !!col_select, col_order = !!col_order,
	          scale_log2r = scale_log2r, impute_na = FALSE,
						df = df, filepath = filepath, filename = filename,
						anal_type = "Corrplot")(use_log10 = use_log10, min_int = min_int,
						                        max_int = max_int, min_log2r = min_log2r,
						                        max_log2r = max_log2r, ...)
}


#'Correlation plots of peptide data
#'
#'@seealso \code{\link{proteoCorrplot}} for parameters
#'
#' @examples
#' pepCorr(
#'   use_log10 = TRUE,
#'   scale_log2r = TRUE,
#'   min_int = 3.5,
#'   max_int = 6.5,
#'   min_log2r = -2,
#'   max_log2r = 2
#' )
#'
#'@export
pepCorr <- function (...) {
	proteoCorrplot(id = pep_seq, ...)
}

#'Correlation plots of protein data
#'
#'@seealso \code{\link{proteoCorrplot}} for parameters
#'
#' @examples
#' prnCorr(
#'   use_log10 = TRUE,
#'   scale_log2r = TRUE,
#'   min_int = 3.5,
#'   max_int = 6.5,
#'   min_log2r = -2,
#'   max_log2r = 2
#' )
#'
#'@export
prnCorr <- function (...) {
	proteoCorrplot(id = gene, ...)
}
