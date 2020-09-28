#' Plots PCA
#'
#' @inheritParams prnPCA
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr ggplot2 
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotPCA <- function (df = NULL, id = NULL, label_scheme_sub = NULL, type = "obs",
                     dimension = 2, folds = 1,
                     show_ids = TRUE, show_ellipses = FALSE,
                     col_group = NULL,
                     col_color = NULL, col_fill = NULL, col_shape = NULL,
                     col_size = NULL, col_alpha = NULL,
                     color_brewer = NULL, fill_brewer = NULL,
                     size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                     filepath = NULL, filename = NULL, 
                     center_features = TRUE, scale_features = TRUE,
                     theme = NULL,
                     anal_type = "PCA", ...) {
  
  stopifnot(vapply(c(scale_log2r, complete_cases, impute_na, show_ids, 
                     center_features, scale_features),
                   rlang::is_logical, logical(1)))
  stopifnot(nrow(label_scheme_sub) > 0)
  
  col_group <- rlang::enexpr(col_group)
  col_color <- rlang::enexpr(col_color)
  col_fill <- rlang::enexpr(col_fill)
  col_shape <- rlang::enexpr(col_shape)
  col_size <- rlang::enexpr(col_size)
  col_alpha <- rlang::enexpr(col_alpha)
  
  complete_cases <- to_complete_cases(complete_cases = complete_cases, impute_na = impute_na)
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  if (show_ellipses && type == "feats") {
    show_ellipses <- FALSE
    warning("No ellipses at `type = feats`.", call. = FALSE)
  }
  
  id <- rlang::enexpr(id)
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]
  
  anal_dots <- dots %>% .[names(.) %in% c("x", "retx", "center", "scale.", "tol", "rank.")]
  fml_dots <- dots[purrr::map_lgl(dots, is_formula)]
  dots <- dots %>% .[! . %in% c(anal_dots, fml_dots)]
  
  if (!is.null(anal_dots$scale.)) {
    if (type == "obs") {
      scale_features <- anal_dots$scale.
      warning("Overwrite `scale_features` with `scale.` ; suggest use only `scale_features`.", 
              call. = FALSE)
    } else if (type == "feats") {
      warning("Argument `scale.` not used; data scaling already set by `scale_log2r`.", 
              call. = FALSE)
    }
    anal_dots$scale. <- NULL
  }
  
  if (!is.null(anal_dots$center)) {
    if (type == "obs") {
      center_features <- anal_dots$center
      warning("Overwrite `center_features` with `center` ; suggest use only `center_features`.", 
              call. = FALSE)
    } else if (type == "feats") {
      warning("Argument `center` not used; data already centered with `standPep()` or `standPrn()`.", 
              call. = FALSE)
    }
    anal_dots$center <- NULL
  }
  
  if (!is.null(anal_dots$x)) {
    anal_dots$x <- NULL
    warning("Argument `x` in `prcomp()` automated.", call. = FALSE)
  }
  
  if (!purrr::is_empty(fml_dots)) {
    fml_dots <- NULL
    warning("The method for class 'formula' is not yet available in proteoQ.", call. = FALSE)
  }
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  res <- df %>%
    filters_in_call(!!!filter_dots) %>%
    arrangers_in_call(!!!arrange_dots) %>%
    scorePCA(
      id = !!id,
      label_scheme_sub = label_scheme_sub,
      anal_type = anal_type,
      scale_log2r = scale_log2r,
      center_features = center_features,
      scale_features = scale_features,
      type = type,
      col_group = !!col_group,
      folds = folds,
      out_file = file.path(filepath, paste0(fn_prefix, "_res.txt")),
      !!!anal_dots)
  
  df <- res$pca
  
  # key `Label` used in `geom_lower_text()`
  if ("Sample_ID" %in% names(df)) {
    df$Label <- df$Sample_ID
  } else if (id %in% names(df)) {
    df$Label <- df[[id]]
  } else {
    df$Label <- df[, 1, drop = FALSE]
  }
  
  res$prop_var <- res$prop_var %>%
    gsub("%", "", .) %>%
    as.numeric() %>%
    paste0("%")
  
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
  
  # --- check dimension ---
  if (dimension < 2) {
    warning("The `dimension` increased from ", dimension, " to a minimum of 2.", call. = FALSE)
    dimension <- 2
  }
  
  ranges <- seq_len(dimension)
  cols <- colnames(df) %>% .[. %in% paste0("PC", ranges)]
  
  max_dim <- names(df) %>% .[grepl("^PC[0-9]+", .)] %>% length()
  if (dimension > max_dim) {
    warning("The `dimension` decreased from ", dimension, " to a maximum of ", max_dim, ".",
            call. = FALSE)
    dimension <- max_dim
  }
  rm(max_dim)
  
  # --- set up aes ---
  if (dimension > 2) {
    mapping <- ggplot2::aes(colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
                            size = !!col_size, alpha = !!col_alpha)
  } else {
    mapping <- ggplot2::aes(x = PC1, y = PC2,
                            colour = !!col_color, fill = !!col_fill, shape = !!col_shape,
                            size = !!col_size, alpha = !!col_alpha)
  }
  
  idx <- purrr::map(mapping, `[[`, 1) %>% purrr::map_lgl(is.null)
  
  mapping_var <- mapping[!idx]
  mapping_fix <- mapping[idx]
  
  if (type == "obs") {
    dot_shape <- 21
    dot_size <- 4
    dot_alpha <- .9
    dot_stroke <- 0.02
    text_size <- 3
  } else {
    dot_shape <- 20
    dot_size <- 2
    dot_alpha <- .6
    dot_stroke <- NA
    text_size = 2
  }
  
  fix_args <- list(colour = "darkgray", fill = NA, shape = dot_shape,
                   size = dot_size, alpha = dot_alpha) %>%
    .[names(.) %in% names(mapping_fix)] %>%
    .[!is.na(.)]
  fix_args$stroke <- dot_stroke
  
  # --- set up axis labels ---
  if (is.null(anal_dots$center)) {
    col_labs <-
      purrr::imap_chr(res$prop_var, ~ paste0("PC", .y, " (", .x, ")")) %>%
      .[ranges]
  } else {
    if (anal_dots$center) {
      col_labs <-
        purrr::imap_chr(res$prop_var, ~ paste0("PC", .y, " (", .x, ")")) %>%
        .[ranges]
    } else {
      col_labs <- cols
    }
  }
  
  # --- plots ---
  if (dimension > 2) {
    p <- GGally::ggpairs(df,
                         axisLabels = "internal",
                         columns = cols,
                         mapping = mapping_var,
                         columnLabels = col_labs,
                         labeller = label_wrap_gen(10),
                         title = "",
                         lower = list(continuous = wrap(geom_lower_text,
                                                        params = fix_args,
                                                        show_ids = show_ids,
                                                        size = text_size)),
                         upper = "blank")
    
    p <- p + theme
    
    if (!is.null(fill_brewer)) {
      for (x in 2:dimension) {
        for (y in 1:(x-1)) {
          p[x, y] <- p[x, y] + scale_fill_brewer(palette = fill_brewer)
        }
      }
    }
    
    if (!is.null(color_brewer)) {
      for (x in 2:dimension) {
        for (y in 1:(x-1)) {
          p[x, y] <- p[x, y] + scale_color_brewer(palette = color_brewer)
        }
      }
    }
    
    if ((!is.null(col_size)) & (!is.null(size_manual))) {
      stopifnot(length(unique(label_scheme_sub[[col_size]])) == length(size_manual))
      for (x in 2:dimension) {
        for (y in 1:(x-1)) {
          p[x, y] <- p[x, y] + scale_size_manual(values = size_manual)
        }
      }
    }
    
    if ((!is.null(col_shape)) & (!is.null(shape_manual))) {
      stopifnot(length(unique(label_scheme_sub[[col_shape]])) == length(shape_manual))
      for (x in 2:dimension) {
        for (y in 1:(x-1)) {
          p[x, y] <- p[x, y] + scale_shape_manual(values = shape_manual)
        }
      }
    }
    
    if ((!is.null(col_alpha)) & (!is.null(alpha_manual))) {
      stopifnot(length(unique(label_scheme_sub[[col_alpha]])) == length(alpha_manual))
      for (x in 2:dimension) {
        for (y in 1:(x-1)) {
          p[x, y] <- p[x, y] + scale_shape_manual(values = alpha_manual)
        }
      }
    }
    
    if (show_ellipses) {
      if (anyNA(label_scheme_sub[[col_group]])) {
        warning("(Partial) NA aesthetics under column `", col_group, "` in expt_smry.xlsx",
                call. = FALSE)
      }
      
      for (x in 2:dimension) {
        for (y in 1:(x-1)) {
          p[x, y] <- p[x, y] + stat_ellipse(
            data = df,
            aes(x = !!rlang::sym(paste0("PC", y)),
                y = !!rlang::sym(paste0("PC", x)),
                fill = !!rlang::sym(col_group)),
            geom = "polygon",
            alpha = .4,
            show.legend = FALSE,
          )
        }
      }
    }
  } else {
    p <- ggplot() +
      rlang::eval_tidy(rlang::quo(geom_point(data = df, mapping = mapping_var, !!!fix_args))) +
      coord_fixed()
    
    if (show_ellipses) {
      if (anyNA(label_scheme_sub[[col_group]])) {
        warning("(Partial) NA aesthetics under column `", col_group, "` in expt_smry.xlsx",
                call. = FALSE)
      }
      
      p <- p + ggplot2::stat_ellipse(
        data = df,
        aes(x = PC1, y = PC2, fill = !!rlang::sym(col_group)),
        geom = "polygon",
        alpha = .4,
        show.legend = FALSE,
      )
    }
    
    if (show_ids) {
      p <- p +
        geom_text(data = df,
                  mapping = aes(x = PC1, y = PC2, label = Label),
                  color = "gray", size = text_size)
    }
    
    p <- p +
      labs(title = "", x = col_labs[1], y = col_labs[2]) + theme
    
    if (!is.null(fill_brewer)) p <- p + scale_fill_brewer(palette = fill_brewer)
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
  }
  
  rlang::eval_tidy(rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(filename)),
                                     plot = p, !!!dots)))
  
  invisible(res)
}


#' Scores PCA
#'
#' @inheritParams prnPCA
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams scoreMDS
#' @import dplyr 
#' @importFrom MASS isoMDS
#' @importFrom magrittr %>% %T>% %$% %<>% 
scorePCA <- function (df, id, label_scheme_sub, anal_type, scale_log2r, 
                      center_features, scale_features, 
                      type, col_group, 
                      folds, out_file, ...) {
  dots <- rlang::enexprs(...)
  id <- rlang::as_string(rlang::enexpr(id))
  col_group <- rlang::enexpr(col_group)
  
  if (! purrr::is_empty(dots$rank.)) {
    if (dots$rank. < 2) {
      warning("PCA `rank.` increased from ", dots$rank., " to a minimum of 2.", 
              call. = FALSE)
      dots$rank. <- 2
    }
  }
  
  df_orig <- df
  
  df <- prepDM(df = df, id = !!id, scale_log2r = scale_log2r,
               sub_grp = label_scheme_sub$Sample_ID, anal_type = anal_type) %>%
    .$log2R
  
  nms <- names(df)
  n_rows <- nrow(df)
  
  if (n_rows <= 50) {
    stop("Need 50 or more data rows for PCA.", call. = FALSE)
  }
  
  label_scheme_sub <- label_scheme_sub %>% dplyr::filter(Sample_ID %in% nms)
  
  if (type == "obs") {
    res <- prep_folded_tdata(df, folds, label_scheme_sub, !!col_group)
    df_t <- res$df_t
    ls_sub <- res$ls_sub
    rm(res)
    
    if (rlang::as_string(col_group) %in% names(df_t)) {
      df_t <- df_t %>% dplyr::select(-!!rlang::sym(col_group))
    }
    
    stopifnot(vapply(df_t, is.numeric, logical(1)))
    
    pr_out <- rlang::expr(stats::prcomp(x = !!df_t, 
                                        center = !!center_features, 
                                        scale. = !!scale_features, 
                                        !!!dots)) %>%
      rlang::eval_bare(env = caller_env())
    
    pr_out$pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      cmbn_meta(ls_sub) %T>%
      readr::write_tsv(out_file)
  } else if (type == "feats") {
    if (folds > 1) {
      message("Coerce to `k_fold = 1` at `type = feats`.\n")
    }
    
    stopifnot(vapply(df, is.numeric, logical(1)))
    
    run_scripts <- FALSE
    if (run_scripts) {
      if (scale_features) {
        df <- df %>% 
          t() %>% 
          scale(center = center_features, scale = TRUE) %>% 
          t()
      } else if (center_features) {
        df <- df %>% sweep(., 1, rowMeans(., na.rm = TRUE), "-")
      }
    }
    
    message("\nAt `type = feats` (peptides/protiens in rows and samples in columns), 
            arguments `center_features` and `scale_features` will not be used.
            Instead, data were already aligned with `method_align` in `standPep()` or `standPrn()` 
            and optionally scaled with `scale_log2r`.", 
            "\nSee also https://proteoq.netlify.app/post/wrapping-pca-into-proteoq/\n")

    pr_out <- rlang::expr(stats::prcomp(x = !!df, 
                                        center = FALSE, 
                                        scale. = !!scale_log2r, 
                                        !!!dots)) %>%
      rlang::eval_bare(env = caller_env())
    
    pr_out$pca <- pr_out$x %>%
      data.frame(check.names = FALSE) %>%
      tibble::rownames_to_column(id) %>%
      dplyr::left_join(df_orig, by = id) %T>%
      readr::write_tsv(out_file)
  } else {
    stop("Unkown `type` for PCA.", call. = FALSE)
  }
  
  pr_out$prop_var <- summary(pr_out)$importance[2, ] %>% round(., digits = 3) %>% scales::percent()
  
  invisible(pr_out)
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
pepPCA <- function (col_select = NULL, col_group = NULL, col_color = NULL,
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL,
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                    center_features = TRUE, scale_features = TRUE,
                    show_ids = TRUE, show_ellipses = FALSE,
                    dimension = 2, folds = 1,
                    df = NULL, filepath = NULL, filename = NULL,
                    theme = NULL, type = c("obs", "feats"), ...) {
  
  check_dots(c("id", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  type <- rlang::enexpr(type)
  if (type == rlang::expr(c("obs", "feats"))) {
    type <- "obs"
  } else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("obs", "feats"))
  }
  
  col_select <- rlang::enexpr(col_select)
  col_group <- rlang::enexpr(col_group)
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
            col_select = !!col_select, col_group = !!col_group,
            col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer,
            size_manual = !!size_manual, shape_manual = !!shape_manual, 
            alpha_manual = !!alpha_manual,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "PCA")(type = type, dimension = dimension, folds = folds, 
                               show_ids = show_ids,
                               show_ellipses = show_ellipses, 
                               center_features = center_features,
                               scale_features = scale_features,
                               theme = theme, ...)
}


#'PCA plots
#'
#'\code{prnPCA} visualizes the principal component analysis (PCA) for protein
#'data.
#'
#'The utility is a wrapper of \code{\link[stats]{prcomp}} against \code{log2FC}.
#'The results are then visualized by either \emph{observations} or
#'\emph{features}. See also
#'https://proteoq.netlify.app/post/wrapping-pca-into-proteoq/ for data centering
#'by either observations or features.
#'
#'@inheritParams prnHist
#'@inheritParams prnHM
#'@inheritParams prnMDS
#'@inheritParams anal_pepNMF
#'@param complete_cases Logical; always TRUE for PCA.
#'@param center_features Logical; if TRUE, adjusts log2FC to center zero by
#'  features (proteins or peptides). The default is TRUE. Note the difference to
#'  data alignment with \code{method_align} in \code{\link[proteoQ]{standPrn}}
#'  or \code{\link[proteoQ]{standPep}} where log2FC are aligned by observations
#'  (samples).
#'@param scale_features Logical; if TRUE, adjusts log2FC to the same scale of
#'  variance by features (protein or peptide entries). The default is TRUE. Note
#'  the difference to data scaling with \code{scale_log2r} where log2FC are
#'  scaled by observations (samples).
#'@param type Character string indicating the type of PCA by either
#'  \emph{observations} or \emph{features}. At the \code{type = obs} default,
#'  observations (samples) are in rows and features (peptides or proteins) in
#'  columns for \code{\link[stats]{prcomp}}. The principal components are then
#'  plotted by observations. Alternatively at \code{type = feats}, features
#'  (peptides or proteins) are in rows and observations (samples) are in
#'  columns. The principal components are then plotted by features.
#'@param folds Not currently used. Integer; the degree of folding data into
#'  subsets. The default is one without data folding.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  Arguments passed to \code{\link[stats]{prcomp}}: \code{rank.}, \code{tol}
#'  etc. At \code{type = obs}, argument \code{scale} becomes
#'  \code{scale_features} and \code{center} matches \code{center_features}. At
#'  \code{type = feats}, the setting of \code{scale_log2r} will be applied for
#'  data scaling and data centering be automated by
#'  \code{\link[proteoQ]{stanPep}} or \code{\link[proteoQ]{stanPrn}}. \cr \cr
#'  Additional arguments for \code{ggsave}: \cr \code{width}, the width of plot;
#'  \cr \code{height}, the height of plot \cr \code{...}
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
#'  \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr
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
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "mascot_psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "mascot_protein_keys.txt", package = "proteoQ") \cr
#'
#'@example inst/extdata/examples/prnPCA_.R
#'
#'@return PCA plots.
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
prnPCA <- function (col_select = NULL, col_group = NULL, col_color = NULL,
                    col_fill = NULL, col_shape = NULL, col_size = NULL, col_alpha = NULL,
                    color_brewer = NULL, fill_brewer = NULL,
                    size_manual = NULL, shape_manual = NULL, alpha_manual = NULL,
                    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE,
                    center_features = TRUE, scale_features = TRUE,
                    show_ids = TRUE, show_ellipses = FALSE,
                    dimension = 2, folds = 1,
                    df = NULL, filepath = NULL, filename = NULL,
                    theme = NULL, type = c("obs", "feats"), ...) {

  check_dots(c("id", "df2", "anal_type"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  type <- rlang::enexpr(type)
  if (type == rlang::expr(c("obs", "feats"))) {
    type <- "obs"
  } else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("obs", "feats"))
  }
  
  col_select <- rlang::enexpr(col_select)
  col_group <- rlang::enexpr(col_group)
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
            col_select = !!col_select, col_group = !!col_group,
            col_color = !!col_color, col_fill = !!col_fill,
            col_shape = !!col_shape, col_size = !!col_size, col_alpha = !!col_alpha,
            color_brewer = !!color_brewer, fill_brewer = !!fill_brewer,
            size_manual = !!size_manual, shape_manual = !!shape_manual, 
            alpha_manual = !!alpha_manual,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "PCA")(type = type, dimension = dimension, folds = folds, 
                               show_ids = show_ids,
                               show_ellipses = show_ellipses, 
                               center_features = center_features,
                               scale_features = scale_features,
                               theme = theme, ...)
}

