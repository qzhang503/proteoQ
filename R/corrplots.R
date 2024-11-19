#' Correlation plots
#' 
#' @param data_select The type of data to be selected, for example, logFC or logInt.
#' @inheritParams prnCorr_logFC
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import stringr dplyr ggplot2 GGally 
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotCorr <- function (df = NULL, id = NULL, anal_type, data_select, 
                      col_select = NULL, col_order = NULL, label_scheme_sub = NULL, 
                      scale_log2r = TRUE, complete_cases = FALSE, 
                      filepath = NULL, filename = NULL, 
                      cor_method = "pearson", digits = 2L, ...) 
{
  if (complete_cases) {
    df <- my_complete_cases(df, scale_log2r, label_scheme_sub)
  }

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  xmin <- eval(dots$xmin, envir = rlang::caller_env()) 
  xmax <- eval(dots$xmax, envir = rlang::caller_env()) 
  xbreaks <- eval(dots$xbreaks, envir = rlang::caller_env())
  width <- eval(dots$width, envir = rlang::caller_env())
  height <- eval(dots$height, envir = rlang::caller_env())
  
  dots <- dots %>% 
    .[! names(.) %in% c("xmin", "xmax", "xbreaks", "width", "height")]
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  arrange_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter_dots, arrange_dots)]
  
  save_data <- dots[["save_data"]]
  
  if (is.null(save_data)) {
    save_data <- FALSE
  } 
  else {
    save_data <- TRUE
    dots[["save_data"]] <- NULL
  }
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  df <- prepDM(df = df, 
               id = !!id, 
               scale_log2r = scale_log2r, 
               sub_grp = label_scheme_sub$Sample_ID, 
               anal_type = anal_type, 
               rm_allna = TRUE) 
  
  if (data_select == "logFC") {
    df <- df$log2R
    y_label <- x_label <- expression("Ratio ("*log[2]*")")
    if (is.null(xmin)) xmin <- -2
    if (is.null(xmax)) xmax <- 2
    if (is.null(xbreaks)) xbreaks <- 1
  } 
  else if (data_select == "logInt") {
    df <- log10(df$Intensity)
    y_label <- x_label <- expression("Intensity ("*log[10]*")")
    if (is.null(xmin)) xmin <- 3.5
    if (is.null(xmax)) xmax <- 6
    if (is.null(xbreaks)) xbreaks <- 1
  } 
  else {
    stop("`data_select` nees to be either`logFC` or `logInt`.")
  }
  
  label_scheme_sub <- label_scheme_sub %>% 
    dplyr::filter(Sample_ID %in% colnames(df))
  
  if (is.null(width)) {
    width <- 1.4 * length(label_scheme_sub$Sample_ID)
  }

  if (is.null(height)) {
    height <- width
  }

  if (ncol(df) > 44L) {
    stop("Maximum number of samples for correlation plots is 44.")
  }

  if (dplyr::n_distinct(label_scheme_sub[[col_order]]) == 1L) {
    df <- df[, order(names(df))]
  }
  else {
    corrplot_orders <- label_scheme_sub %>%
      dplyr::select(Sample_ID, !!col_select, !!col_order) %>%
      dplyr::filter(!is.na(!!col_order)) %>%
      unique(.) %>%
      dplyr::arrange(!!col_order)
    
    df <- df[, as.character(corrplot_orders$Sample_ID), drop = FALSE]
  }
  
  if (save_data) {
    readr::write_tsv(df, file.path(filepath, paste0(fn_prefix, "_data.txt")))
  }

  plot_corr_sub(df = df, 
                cor_method = cor_method, 
                xlab = x_label, ylab = y_label,
                filename = filename, filepath = filepath,
                xmin = xmin, xmax = xmax, xbreaks = xbreaks, 
                width = width, height = height, digits = digits, !!!dots)
}


#' Custom correlation plots
#' 
#' @param data A data frame
#' @param mapping A mapping aesthetics.
#' @param color A color.
#' @param sizeRange A range of sizes.
#' @param cor_method A correlation method.
#' 
#' @examples 
#' \donttest{
#' my_custom_cor(data, aes(x = col_nm_1, y = col_nm_1))
#' }
my_custom_cor <- function(data, mapping, color = I("grey50"), sizeRange = c(1, 4), 
                          cor_method = "pearson", digits = 2L, ...) 
{
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  x[is.infinite(x)] <- NA_real_
  y[is.infinite(y)] <- NA_real_
  
  ct <- cor.test(x, y, method = cor_method)
  
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  
  r <- unname(ct$estimate)
  rt <- format(r, digits = digits)[1]
  
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
#' @import stringr dplyr ggplot2 GGally purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
plot_corr_sub <- function (df, cor_method = "pearson", 
                           xlab, ylab, filename, filepath, 
                           xmin, xmax, xbreaks, width, height, 
                           digits = 2L, ...) 
{
  # not used
  my_fn <- function(data, mapping, method = "lm", ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(alpha = 0.3, size = .1) +
      geom_smooth(method = method, ...)
    p
  }

  # not used
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

  # not used
  panel_cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
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

  # not used
  panel_hist <- function(x, ...) {
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

  # not used
  panel_lm <- function(x, y, col = par('col'), bg = NA, pch = par('pch'),
                       cex = 1, col.smooth = 'black', ...){
    points(x, y, pch = pch, col = "red",  bg = bg, cex = cex)
    # abline(stats::lm(y~x), col=col.smooth,...)
    bg = "transparent"
  }

  # not used
  my_lower <- function(data, mapping, method = "lm", ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(size = .02, alpha = .5) +
      geom_abline(alpha = .5, linetype = "dashed", color = "gray", size = 1) +
      geom_smooth(size = 1, method = method, ...)
    p
  }

  my_lower_no_sm <- function(data, mapping, method = "lm", ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(size = .02, alpha = .5) +
      geom_abline(alpha = .5, linetype = "dashed", color = "gray", size = 1)
    p
  }

  # not u sed
  my_diag <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_density(fill = "#fec44f", size = .02, alpha = .5, adjust = 3)
    p
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
  
  df <- df %>% dplyr::mutate_all(~ replace(.x, is.infinite(.), NA_real_))
  sids <- as.character(names(df))
  
  p1 <- ggpairs(df, columnLabels = sids, 
                labeller = label_wrap_gen(10),
                title = "", xlab = xlab, ylab = ylab, 
                lower = list(continuous = my_lower_no_sm),
                upper = list(continuous = wrap(my_custom_cor, 
                                               cor_method = cor_method, 
                                               digits = digits)))
  
  p2 <- ggcorr(df, label = TRUE, label_round = 2)
  
  local({
    cors <- df %>% 
      cor(use = "pairwise.complete.obs", method = cor_method) %>% 
      data.frame(check.names = FALSE) %>% 
      tibble::rownames_to_column("Sample_ID")

    filename <- gsub("(.*\\.)[^.]*$", paste0("\\1", "txt"), filename)
    
    readr::write_delim(cors,file.path(filepath, filename), 
                       delim = find_delim(filename))
  })

  g2 <- ggplotGrob(p2)
  colors <- g2$grobs[[6]]$children[[3]]$gp$fill
  
  ###
  c(
    "#F73F1E","#BFD5DC","#BED4DC","#BDD4DC","#BED4DC","#F63C1C","#B5D0D8","#ACCBD5","#FFAF98","#FFB7A3","#AFCCD6","#AACAD4",
    "#FFAD96","#FFB39D","#FA5835","#DBE4E7","#E4E9EA","#A4C7D2","#ADCCD6","#C1D6DD","#C8DAE0","#D5E0E4","#D9E2E6","#ADCBD5",
    "#B7D1D9","#CFDDE2","#D0DEE3","#F73F1E"
  )
  
  ###
  colors2 <- vector("character", length(colors))
  
  # deduce col and row positions from the index of colors
  ncol2 <- ncol - 1L
  cols  <- rep.int(1:ncol2, times = 1:ncol2)
  rows  <- unlist(lapply(seq_len(ncol2), seq_len))
  
  for (i in seq_along(colors)) {
    col <- cols[[i]]
    row <- rows[[i]]
    gap <- if (row == 1L) 0L else (row - 1) * ncol2 - sum(seq_len(row - 2))
    ps2 <- gap + col - row + 1L
    colors2[[ps2]] <- colors[[i]]
  }
  
  colors <- colors2
  ###
  
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

  for (x in 1:ncol) {
    p1[x, x] <- p1[x, x] +
      scale_x_continuous(limits = c(xmin-.2, xmax+.2), breaks = c(xmin, 0, xmax))
  }

  p1 <- p1 +
    theme(plot.title = element_text(face = "bold", colour = "black", size = 20,
                                    hjust = 0.5, vjust = 0.5))

  ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
  
  suppressWarnings(
    rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(filename)),
                      plot = p1, 
                      width = width, 
                      height = height, 
                      !!!ggsave_dots)) %>% 
      rlang::eval_tidy()
  )
  
  invisible(p2$data)
}


#'Correlation plots
#'
#'\code{pepCorr_logFC} plots correlation for peptide \code{logFC}. data.
#'
#'@rdname prnCorr_logFC
#'
#'@import purrr
#'@export
pepCorr_logFC <- function (col_select = NULL, col_order = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, 
                           impute_na = FALSE, df = NULL, filepath = NULL, 
                           filename = NULL, cor_method = "pearson", 
                           digits = 2L, ...) 
{
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  cor_method <- rlang::as_string(rlang::enexpr(cor_method))
  
  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_order = !!col_order,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "Corrplot")(data_select = "logFC", 
                                    cor_method = cor_method, 
                                    digits = digits,
                                    ...)
}


#'Correlation plots
#'
#'\code{pepCorr_logInt} plots correlation of the \code{log10} intensity of ions
#'for peptide data.
#'
#'@rdname prnCorr_logFC
#'
#'@import purrr
#'@export
pepCorr_logInt <- function (col_select = NULL, col_order = NULL, 
                            scale_log2r = TRUE, complete_cases = FALSE, 
                            impute_na = FALSE, df = NULL, filepath = NULL, 
                            filename = NULL, cor_method = "pearson", 
                            digits = 2L, ...) 
{
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  cor_method <- rlang::as_string(rlang::enexpr(cor_method))
  
  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_order = !!col_order,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "Corrplot")(data_select = "logInt", 
                                    cor_method = cor_method, 
                                    digits = digits,
                                    ...)
}


#'Correlation plots
#'
#'\code{prnCorr_logFC} plots correlation for protein \code{logFC}.
#'
#'With TMT experiments, the same polypeptide may be triggered for MS2 any where
#'between the baseline and the apex levels during a peak elution. The direct
#'comparison of reporter-ion intensities between plex-es might have little
#'meaning.
#'
#'@inheritParams prnHist
#'@inheritParams prnMDS
#'@param col_order Character string to a column key in \code{expt_smry.xlsx}.
#'  Numeric values under which will be used for the left-to-right arrangement of
#'  samples in graphic outputs or top-to-bottom arrangement in text outputs. At
#'  the NULL default, the column key \code{Order} will be used. If values under
#'  column \code{Order} are left blank, samples will be ordered by their names.
#'@param cor_method A character string indicating which correlation coefficient
#'  is to be computed. One of \code{"pearson"} (default), \code{"kendall"}, or
#'  \code{"spearman"}.
#'@param digits The number of decimal places in correlation values to be
#'  displayed.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  Additional parameters for plotting: \cr \code{width}, the width of plot \cr
#'  \code{height}, the height of plot \cr \code{xmin}, the minimum \eqn{x} of
#'  logFC or intensity \cr \code{xmax}, the maximum \eqn{x} of logFC data or
#'  intensity data \cr \code{xbreaks}, the breaks on \eqn{x} axis; the same
#'  breaks will be applied to \eqn{y} axis.
#'
#'@seealso \emph{Metadata} \cr \code{\link{load_expts}} for metadata preparation
#'  and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr \code{\link{normPSM}} for extended examples in
#'  PSM data normalization \cr \code{\link{PSM2Pep}} for extended examples in
#'  PSM to peptide summarization \cr \code{\link{mergePep}} for extended
#'  examples in peptide data merging \cr \code{\link{standPep}} for extended
#'  examples in peptide data normalization \cr \code{\link{Pep2Prn}} for
#'  extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization.
#'  \cr \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples
#'  in data purging \cr \code{\link{pepHist}} and \code{\link{prnHist}} for
#'  extended examples in histogram visualization. \cr \code{\link{extract_raws}}
#'  and \code{\link{extract_psm_raws}} for extracting MS file names \cr
#'
#'  \emph{Variable arguments of `filter_...`} \cr \code{\link{contain_str}},
#'  \code{\link{contain_chars_in}}, \code{\link{not_contain_str}},
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}},
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings
#'  \cr
#'
#'  \emph{Missing values} \cr \code{\link{pepImp}} and \code{\link{prnImp}} for
#'  missing value imputation \cr
#'
#'  \emph{Informatics} \cr \code{\link{pepSig}} and \code{\link{prnSig}} for
#'  significance tests \cr \code{\link{pepVol}} and \code{\link{prnVol}} for
#'  volcano plot visualization \cr \code{\link{prnGSPA}} for gene set enrichment
#'  analysis by protein significance pVals \cr \code{\link{gspaMap}} for mapping
#'  GSPA to volcano plot visualization \cr \code{\link{prnGSPAHM}} for heat map
#'  and network visualization of GSPA results \cr \code{\link{prnGSVA}} for gene
#'  set variance analysis \cr \code{\link{prnGSEA}} for data preparation for
#'  online GSEA. \cr \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS
#'  visualization \cr \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA
#'  visualization \cr \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA
#'  visualization \cr \code{\link{pepHM}} and \code{\link{prnHM}} for heat map
#'  visualization \cr \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'  \code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for
#'  correlation plots \cr \code{\link{anal_prnTrend}} and
#'  \code{\link{plot_prnTrend}} for trend analysis and visualization \cr
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}},
#'  \code{\link{plot_pepNMFCon}}, \code{\link{plot_prnNMFCon}},
#'  \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr
#'
#'  \emph{Custom databases} \cr \code{\link{Uni2Entrez}} for lookups between
#'  UniProt accessions and Entrez IDs \cr \code{\link{Ref2Entrez}} for lookups
#'  among RefSeq accessions, gene names and Entrez IDs \cr \code{\link{prepGO}}
#'  for
#'  \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'   ontology}} \cr \code{\link{prepMSig}} for
#'  \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular
#'   signatures} \cr \code{\link{prepString}} and \code{\link{anal_prnString}}
#'  for STRING-DB \cr
#'
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@example inst/extdata/examples/prnCorr_.R
#'
#'@return Correlation plots.
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
prnCorr_logFC <- function (col_select = NULL, col_order = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, 
                           impute_na = FALSE, df = NULL, filepath = NULL, 
                           filename = NULL, cor_method = "pearson", 
                           digits = 2L, ...) 
{
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  cor_method <- rlang::as_string(rlang::enexpr(cor_method))
  
  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_order = !!col_order,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "Corrplot")(data_select = "logFC", 
                                    cor_method = cor_method, 
                                    digits = digits, 
                                    ...)
}


#'Correlation Plots
#'
#'\code{prnCorr_logInt} plots correlation of the \code{log10} intensity of ions
#'for protein data.
#'
#'
#'@rdname prnCorr_logFC
#'
#'@import purrr
#'@export
prnCorr_logInt <- function (col_select = NULL, col_order = NULL, 
                            scale_log2r = TRUE, complete_cases = FALSE, 
                            impute_na = FALSE, df = NULL, filepath = NULL, 
                            filename = NULL, cor_method = "pearson", 
                            digits = 2L, ...) 
{
  check_dots(c("id", "anal_type", "data_select", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  cor_method <- rlang::as_string(rlang::enexpr(cor_method))
  
  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_order = !!col_order,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "Corrplot")(data_select = "logInt", 
                                    cor_method = cor_method, 
                                    digits = digits, 
                                    ...)
}


