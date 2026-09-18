#' Bar plots of median subcellular statistics.
#' 
#' @param col_cluster The column key of cluster nodes.
#' @param pat A pattern of input files.
#' @param ext The extension of input files.
#' @param data_type The type of data.
#' @param levels_fraction The levels of subcellular fractions.
#' @param tie_method The method to handle ties.
#' @param panel_ids Panel IDs for plotting.
#' @param qt The quantile of localization score for thresholding subcellular
#'   compartment specificity. The default is to use median.
#' @param make_plot Logical; to bypass the plotting at FALSE.
#' @inheritParams plot_prnTrend
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr purrr ggplot2 RColorBrewer
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotSubcellular <- function(
    id = "gene", dat_dir = NULL, levels_group = NULL, 
    col_group = NULL, col_order = NULL, col_cluster = "cluster", 
    col_fraction = NULL, levels_fraction = NULL, 
    col_subtype = NULL, levels_subtype = NULL, 
    tie_method = "none", 
    label_scheme_sub = NULL, n_clust = NULL, panel_ids = NULL, 
    pat = "Trend_[ONZ]_.*nclust\\d+\\.txt$", ext = "txt", data_type = "Trend", 
    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
    df2 = NULL, filepath = NULL, filename = NULL, qt = .5, 
    make_plot = TRUE, theme = NULL, ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")

  dots <- rlang::enexprs(...)

  # Find input `df2`
  df2 <- find_trend_df2(
    df2 = df2, n_clust = n_clust, pat = pat, ext = ext, 
    scale_log2r = scale_log2r, impute_na = impute_na, filepath = filepath)
  
  # Prepare output file name
  custom_prefix <- if (id %in% c("pep_seq", "pep_seq_mod")) {
    purrr::map_chr(
      df2, ~ sub(paste0("^(.*?)_?Peptide_", data_type, ".*$"), "\\1", .x))
  } else if (id %in% c("prot_acc", "gene")) {
    purrr::map_chr(
      df2, ~ sub(paste0("^(.*?)_?Protein_", data_type, ".*$"), "\\1", .x))
  } else {
    stop("Unknown id = ", id)
  }
  
  fn_prefix <- tools::file_path_sans_ext(filename)
  fn_suffix <- tools::file_ext(filename)
  
  ## Plot data
  if (is.null(theme)) theme <- theme_bw()

  tempsub <- mapply(
    plotSubcellular_sub, 
    df2 = df2, 
    custom_prefix = custom_prefix, 
    MoreArgs = list(
      fn_prefix = fn_prefix, 
      fn_suffix = fn_suffix, 
      df = df, 
      id = id, 
      qt = qt, 
      levels_group = levels_group, 
      col_group = col_group, 
      col_order = col_order, 
      col_cluster = col_cluster,
      col_fraction = col_fraction, 
      levels_fraction = levels_fraction,
      col_subtype = col_subtype, 
      levels_subtype = levels_subtype, 
      tie_method = tie_method, 
      complete_cases = complete_cases, 
      panel_ids = panel_ids, 
      filepath = filepath, 
      label_scheme_sub = label_scheme_sub, 
      data_type = data_type, 
      pat = pat, 
      ext = ext, 
      theme = theme, 
      make_plot = make_plot, 
      dots = dots), 
    SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  df_meds <- lapply(tempsub, `[[`, "df_med")
  # become "svm" for non-Trend data such as svmProbs
  col_fraction <- tempsub[[1]][["col_fraction"]]
  col_subtype  <- tempsub[[1]][["col_subtype"]]
  col_x <- tempsub[[1]][["col_x"]]
  col_y <- tempsub[[1]][["col_y"]]

  ## Classify trends
  tempdata <- lapply(
    df_meds, 
    hclassify_trends, 
    col_subtype = col_subtype, 
    col_fraction = col_fraction, 
    col_group = col_group, 
    col_order = col_order, 
    col_cluster = col_cluster, 
    col_x = col_x, 
    col_y = col_y, 
    tie_method = tie_method)

  res_score_co <- lapply(tempdata, `[[`, "res")
  col_group <- tempdata[[1]][["col_group"]]
  col_order <- tempdata[[1]][["col_order"]]
  col_fraction <- tempdata[[1]][["col_fraction"]]
  col_subtype <- tempdata[[1]][["col_subtype"]]
  col_x <- tempdata[[1]][["col_x"]]
  col_y <- tempdata[[1]][["col_y"]]

  ## Update subcellular locations
  res_score_co <- res_score_co |>
    dplyr::bind_rows(.id = "n_clust") |>
    dplyr::mutate(
      n_clust = 
        as.integer(gsub(paste0(".*_nclust(\\d+)[^\\d]*\\.", ext, "$"), "\\1", 
                        n_clust))) |>
    dplyr::arrange(!!rlang::sym(col_fraction), !!rlang::sym(col_y))
  
  # readr::write_tsv(res_score_co, file.path(filepath, "all_loc_scores.tsv"))
  
  res_score_co <- res_score_co |>
    dplyr::group_by_at(c(col_fraction, col_subtype)) |>
    dplyr::slice(1L) |>
    dplyr::rename(score_co = !!rlang::sym(col_y), ) |>
    dplyr::ungroup() |>
    dplyr::arrange(score_co)

  join_cols <- c(col_fraction, col_subtype)
  names(join_cols) <- c(col_fraction, col_subtype)
  
  dfs <- lapply(df2, function (fn) {
    df <- readr::read_tsv(file.path(filepath, fn)) |>
      dplyr::left_join(
        res_score_co[, c(col_fraction, "score_co", col_subtype)], 
        by = join_cols)

    # Keep the original data structure; do not filter
    # may split by "score_co" if no ties
    df <- lapply(split(df, df[[col_fraction]]), function (x) {
      frac <- x[[col_fraction]][[1]]
      sco  <- x[["score_co"]][[1]]
      
      x <- x |>
        dplyr::mutate(
          sub_location = ifelse(loc_score >= sco, frac, NA_character_))
    }) |>
      dplyr::bind_rows() |>
      dplyr::select(-c("score_co"))
    
    # Not yet for svmProb: need column 'log10Int'
    if (data_type == "Trend") {
      df <- threshold_subcell_by_int(df, id = id, col_fraction = col_fraction)
    }
    
    # Update data
    readr::write_tsv(df, file.path(filepath, fn))
    
    df
  })
  
  invisible(dfs)
}


#' Threshold subcellular locations by intensity.
#' 
#' @param df A data frame.
#' @param id An identifier
#' @param col_fraction The column key indicating subcellular locations.
#' @param fct A multiplification factor in intensity thresholding.
threshold_subcell_by_int <- function (df, id = "gene", col_fraction = NULL, 
                                      fct = 1.0)
{
  if (!"log10Int" %in% names(df)) {
    warning("Column 'log10Int' not found.")
    return(df)
  }
  
  dfs <- split(df, df[[id]])
  
  lapply(dfs, function (dfx) {
    nas <- is.na(dfx[["sub_location"]])
    
    if (all(nas) || !any(nas)) {
      return(dfx)
    }
    
    dfa <- dfx[!nas, ]
    dfb <- dfx[nas, ]
    frs <- dfb[[col_fraction]]
    oks <- dfb[["log10Int"]] >= min(dfa[["log10Int"]], na.rm = TRUE) + log10(fct)
    
    dfb[["sub_location"]][oks] <- frs[oks]
    dfx[["sub_location"]][nas] <- dfb[["sub_location"]]
    
    dfx
  }) |>
    dplyr::bind_rows()
}


#' Plot UMAP of subcellular scores.
#'
#' Against localization scores.
#' 
#' @param dat_dir The current working directory.
#' @param filepath The file path to Trend data.
#' @param filename_ref Optional; a file name with prepending path to a reference
#'   subcellular annotation file, e.g., UniProt and GO intersection.
#' @param id An identifier.
#' @param key_subcellular The key of subcellular compartments from the reference
#'   data.
#' @param key_col The column key where its values are used for UMAP.
#' @param seed A seed for umap.
#' @inheritParams plot_prnTrend
#' @inheritParams plot_prnSubcellular
plot_prnSubcellular_UMAP <- function (
    df2 = NULL, id = "gene", levels_fraction = NULL, levels_subtype = NULL, 
    dat_dir = NULL, filepath = NULL, 
    key_subcellular = "location_slim3", key_col = "loc_score", 
    levels_ref_compartments = c("CP", "NP", "ChA"), 
    scale_log2r = TRUE, impute_na = FALSE, 
    filename_ref = NULL, filename = NULL, seed = 42, ...) 
{
  dots <- rlang::enexprs(...)
  lang_dots     <- dots[unlist(lapply(dots, is.language))]
  filter2_dots  <- lang_dots[grepl("^filter2_", names(lang_dots))]
  arrange2_dots <- lang_dots[grepl("^arrange2_", names(lang_dots))]
  dots <- dots[!dots %in% c(filter2_dots, arrange2_dots)]
  
  if (is.null(dat_dir)) dat_dir  <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  label_scheme_sub <- load_ls_group(dat_dir = dat_dir)
  
  ## Special handling of character vectors
  if (FALSE) {
    levels_fraction <- eval(rlang::enexpr(levels_fraction))
    levels_subtype <- eval(rlang::enexpr(levels_subtype))
  }

  ## Find input `df2`
  df2 <- find_trend_df2(
    df2 = df2, n_clust = NULL, scale_log2r = scale_log2r, 
    impute_na = impute_na, filepath = filepath)
  
  if (!length(df2)) {
    warning("No trend results found.", " Run 'anal_prnTrend' first.")
    return(NULL)
  }
  
  if (is.null(filename)) {
    fn_prefix <- paste0(tools::file_path_sans_ext(df2[[1]]), "_umap")
    fn_suffix <- "png"
    filename  <- paste0(fn_prefix, ".", fn_suffix)
  } else {
    fn_prefix <- tools::file_path_sans_ext(filename)
    fn_suffix <- tools::file_ext(filename)
  }
  
  # A data frame containing either log2FC or intensity values
  df <- readr::read_tsv(file.path(filepath, df2[[1]]))
  
  df <- df |>
    filters_in_call(!!!filter2_dots) |>
    arrangers_in_call(!!!arrange2_dots)
  
  # Obtain the column keys encoding fraction and sub-type information
  col_fraction <- df[["col_fraction"]][[1]]
  col_subtype  <- df[["col_subtype"]][[1]]
  
  check_subcell(
    levels_fraction = levels_fraction, levels_subtype = levels_subtype, 
    col_subtype = col_subtype)
  
  df <- df |>
    dplyr::filter(.data[[col_subtype]] %in% levels_subtype)

  # String `col_fraction` is the key that encodes sample fraction; 
  # Column `sub_location` contains classification results
  
  # Remove is.na(sub_location) can cause extreme values in UMAP
  # gene  Veh_CP TG_CP TG_ChA TG_NP Veh_NP Veh_ChA
  # ALG9     0.2     0      0     0  0.209       0
  if (key_col == "log10Int") {
    df_expr <- hmake_subcell_expr(
      df = df, label_scheme_sub = label_scheme_sub, id = id, key_col = key_col, 
      col_subtype = col_subtype, levels_subtype = levels_subtype, 
      col_fraction = col_fraction, levels_fraction = levels_fraction)
    
    saveRDS(df_expr, file.path(filepath, "df_exprs_int.rds"))
    
    df_expr <- do.call(cbind, df_expr) |>
      as.data.frame() |>
      tibble::rownames_to_column(id)
    
    readr::write_tsv(df_expr, file.path(filepath, "df_expr_int.tsv"))
    saveRDS(df_expr, file.path(filepath, "df_expr_int.rds"))
  } else {
    df_expr <- df[, c(id, key_col, "sub_location", "group")] |>
      dplyr::filter(!is.na(sub_location)) |>
      dplyr::arrange_at(id) |>
      unique() |>
      tidyr::pivot_wider(
        id_cols = id, 
        names_from = c("group"), 
        values_from = key_col,
        values_fill = 0, )
    
    readr::write_tsv(df_expr, file.path(filepath, "df_expr_log2Ratio.tsv"))
    saveRDS(df_expr, file.path(filepath, "df_expr_log2Ratio.rds"))
  }
  
  file_umap <- file.path(filepath, paste0(fn_prefix, "_", key_col, ".tsv"))
  
  if (file.exists(file_umap)) {
    umap_df <- readr::read_tsv(file_umap)
  } else {
    set.seed(seed)
    res_umap <- umap::umap(df_expr[, !colnames(df_expr) == "gene", drop = FALSE])
    umap_df <- res_umap$layout[, 1:2]
    colnames(umap_df) <- paste0("UMAP", 1:2)
    rownames(umap_df) <- df_expr$gene
    umap_df <- tibble::rownames_to_column(data.frame(umap_df), "gene")
    
    readr::write_tsv(umap_df, file_umap)
    
    umap_df
  }
  
  ## Plot of predicted subcellular results
  ok_subtypes <- sort(unique(df[[col_subtype]]))
  
  if (n_subtypes <- length(levels_subtype)) {
    if (length(bads <- levels_subtype[!levels_subtype %in% ok_subtypes])) {
      warning("levels_subtype not found: ", paste0(bads, collapse = ", "))
      levels_subtype <- ok_subtypes
      n_subtypes <- length(levels_subtype)
    }
  }

  ncol   <- dots$ncol
  nrow   <- dots$nrow
  size   <- dots$size
  alpha  <- dots$alpha
  width  <- dots$width
  height <- dots$height
  dpi    <- dots$dpi
  
  if (is.null(width))  width <- n_subtypes * 8 / nrow + 2
  if (is.null(height)) height <- 8 * nrow
  if (is.null(ncol))  ncol <- 2
  if (is.null(nrow))  nrow <- 1
  if (is.null(size))  size <- .01
  if (is.null(alpha)) alpha <- .5
  if (is.null(dpi)) dpi <- 300
  
  if (is.null(ncol)) {
    if (is.null(nrow)) {
      ncol <- 4L
      nrow <- ceiling(n_subtypes / ncol)
    } else {
      ncol <- ceiling(n_subtypes / nrow)
    }
  } else {
    nrow <- ceiling(n_subtypes / ncol)
  }
  
  dots$ncol <- NULL
  dots$nrow <- NULL
  dots$width <- NULL
  dots$height <- NULL
  dots$alpha <- NULL
  dots$size <- NULL
  dots$dpi <- NULL
  
  df <- df |> 
    dplyr::mutate(
      !!col_subtype := 
        factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
      !!col_fraction := 
        factor(!!rlang::sym(col_fraction), levels = levels_fraction),
      sub_location = 
        factor(sub_location, levels = levels_fraction),)

  # Special handling of
  # gene  Veh_CP TG_CP TG_ChA TG_NP Veh_NP Veh_ChA
  # ALG9     0.2     0      0     0  0.209       0
  umap_df <- umap_df[rowSums(
    abs(umap_df[, grepl("^UMAP\\d+", colnames(umap_df)), drop = FALSE]) > 100)
    == 0, ]
  
  df_plot <- df |>
    dplyr::left_join(umap_df, by = "gene") |> 
    dplyr::filter(!is.na(sub_location))
  
  my_colors <- local({
    set1_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                     "#FFFF33", "#A65628", "#F781BF")
    add_colors  <- c("#008080", "#3A1188", "#808000", "#882255", "#00BFFF", 
                     "#99d8c9", "#fee6ce")
    last_color  <- "#999999"
    
    n_comps <- length(levels_ref_compartments)
    
    if (any(is.na(n_comps)) || "Unknown" %in% levels_ref_compartments) {
      rng_colors <- 1:(n_comps - 1L)
    } else {
      rng_colors <- 1:n_comps
    }
    
    my_colors <- c(set1_colors, add_colors)[rng_colors]
    my_colors <- c(my_colors, last_color)
    
    if (any(is.na(my_colors))) {
      warning("More compartments than the number of colors.")
    }
    
    my_colors
  })
  
  p <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = sub_location)) +
    geom_point(shape = 46, alpha = alpha, size = size) +
    scale_color_manual(values = my_colors, drop = FALSE) + 
    guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
    theme_bw() +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Subcellular\nLocation"
    ) +
    # facet_wrap(vars(!!rlang::sym(col_subtype)), nrow = nrow, labeller = label_value) + 
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
  
  if (n_subtypes < 2L) {
    px <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_fraction)), ncol = 3, 
        labeller = label_value, drop = FALSE)
  } else {
    px <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_subtype), !!rlang::sym(col_fraction)), ncol = 3, 
        
        labeller = labeller(.default = label_value, .multi_line = FALSE))
    
    p <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_subtype)), ncol = ncol, 
        labeller = labeller(.default = label_value, .multi_line = FALSE), )
  }
  
  ggsave(plot = px, file.path(filepath, paste0(fn_prefix, "_", key_col, ".", fn_suffix)), 
         width = width, height = height, dpi = dpi) # 8, 5
  
  ggsave(plot = p, file.path(filepath, paste0(fn_prefix, "_", key_col, "_combined.", fn_suffix)), 
         width = height + 1, height = height - 2, dpi = dpi) # 6, 3
  
  if (FALSE) {
    # before classification
    dfx <- df[, c("gene", col_fraction, col_subtype, "purity", "entropy", 
                  key_col, "sub_location")] |>
      unique()
    out_nmx <- paste0(gsub("umap$", "hist", fn_prefix), "_all.", fn_suffix)
    
    # after classification
    dfx <- df_plot
    out_nmx <- paste0(gsub("umap$", "hist", fn_prefix), "_fil.", fn_suffix)

    ggplot(dfx, 
           aes(x = !!rlang::sym(key_col), fill = !!rlang::sym(col_subtype))) +
      geom_histogram(aes(y = after_stat(count)), 
                     color = "white", linewidth = 0.3, binwidth = 0.01, 
                     alpha = .9) +
      scale_x_continuous(
        labels = scales::label_number(accuracy = .2),
        breaks = scales::breaks_width(.4)
      ) +
      facet_wrap(vars(!!rlang::sym(col_fraction)), nrow = 1, 
                 labeller = label_value) + 
      scale_fill_brewer(palette = "Set2") +
      labs(x = "Score", y = "Count") +
      theme_classic() + 
      theme(
        plot.title = element_blank(),
        axis.title.x = element_text(size = 14, margin = margin(t = 14)), 
        axis.title.y = element_text(size = 14, margin = margin(r = 14)),
        strip.text = element_text(size = 12, hjust = 0), # face = "bold", 
        axis.line = element_line(linewidth = 0.3, color = "grey30"),
        axis.text = element_text(color = "grey20"),
        legend.title = element_blank(),
        legend.position = "top", 
        legend.justification = "left",
        panel.spacing.y = unit(1.2, "lines"),
        plot.margin = margin(t = 5, r = 15, b = 5, l = 5, unit = "pt")
      )
    
    ggsave(file.path(filepath, out_nmx), width = width, height = height, 
           dpi = dpi)
  }
  
  ## Map to the reference of UniProt + GO
  # Only at `loc_score`
  if ((!is.null(filename_ref)) && file.exists(filename_ref)) {
    filename2 <- paste0(fn_prefix, "_", key_col, "_ref.", fn_suffix)

    df_ref <- readr::read_tsv(file.path(filename_ref)) |>
      dplyr::filter(source != "HPA")
    
    if (use_intersect <- TRUE) {
      df_ref <- local({
        dfs <- split(df_ref, df_ref$source)
        gns <- lapply(dfs, `[[`, "gene")
        gns <- purrr::reduce(gns, intersect)
        
        df_ref <- df_ref |> 
          dplyr::filter(gene %in% gns)
      })
    }
    
    df_ref <- df_ref |>
      dplyr::rename(sub_location = !!rlang::sym(key_subcellular)) |>
      dplyr::select(c("gene", "sub_location")) |>
      unique() |>
      dplyr::filter(gene %in% unique(df$gene))
    
    df_plot_ref <- df_plot |>
      dplyr::select(-one_of("sub_location")) |>
      dplyr::left_join(df_ref, by = "gene") |>
      dplyr::filter(!is.na(sub_location)) |>
      dplyr::select(dplyr::one_of(
        c("gene", "UMAP1", "UMAP2", "sub_location", col_subtype))) |>
      unique() |>
      dplyr::mutate(
        !!col_subtype := 
          factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
        sub_location = 
          factor(sub_location, levels = levels_fraction),
      )

    p_ref <- ggplot(df_plot_ref, aes(x = UMAP1, y = UMAP2, color = sub_location)) +
      geom_point(shape = 46, alpha = alpha, size = size) +
      scale_color_manual(
        values = c("CP" = "#E41A1C", "NP" = "#377EB8", "ChA" = "#4DAF4A")
      ) +
      # Override the legend point size, shape, and transparency
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
      theme_bw() +
      labs(
        title = NULL,
        x = "UMAP 1",
        y = "UMAP 2",
        color = "Subcellular\nLocation"
      ) +
      facet_wrap(vars(!!rlang::sym(col_subtype[1])), nrow = nrow, 
                 labeller = label_value) + 
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.title = element_blank(),
        panel.grid = element_blank()
      )
    
    ggsave(file.path(filepath, filename2), width = height + 1, height = height - 2, dpi = dpi)
    
    ## Combined plot
    if (combine_plots <- TRUE) {
      df_plot_all <- dplyr::bind_rows(
        df_plot |> dplyr::mutate(Source = "This work"), 
        df_plot_ref |> dplyr::mutate(Source = "UP and GO"), 
      )
      
      p_all <- ggplot(df_plot_all, aes(x = UMAP1, y = UMAP2, color = sub_location)) +
        geom_point(shape = 46, alpha = alpha, size = size) +
        scale_color_manual(
          values = c("CP" = "#E41A1C", "NP" = "#377EB8", "ChA" = "#4DAF4A")
        ) +
        # Override the legend point size, shape, and transparency
        guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
        theme_bw() +
        labs(
          title = NULL,
          x = "UMAP 1",
          y = "UMAP 2",
          color = "Subcellular\nLocation"
        ) +
        facet_wrap(vars(Source), nrow = nrow, 
                   labeller = label_value) + 
        theme(
          plot.title = element_text(face = "bold", size = 14),
          legend.title = element_blank(),
          panel.grid = element_blank()
        )
      
      filename_all <- paste0(fn_prefix, "_", key_col, "_all.", fn_suffix)
      ggsave(file.path(filepath, filename_all), 
             width = height + 1, height = height - 2.1, dpi = dpi)
    }
    
    
  }
}


#' Plots median trends at a given \code{n_clust}.
#'
#' Input data include results from \link{anal_prnTrend} at \code{is_subcellular
#' = TRUE} or the probabilities from \link{svmSubcell}.
#'
#' @param custom_prefix A custom filename prefix.
#' @param fn_prefix A file name prefix.
#' @param fn_suffix A file name suffix.
#' @param filepath An output file path.
#' @param df A data frame.
#' @param label_scheme_sub A metadata subset.
#' @param qt The quantile of localization score for thresholding subcellular
#'   compartment specificity. The default is to use median.
#' @param make_plot Logical; to bypass the plotting at FALSE.
#' @param dots Variable arguments.
#' @inheritParams plotTrend
#' @inheritParams plotSubcellular
#' @importFrom magrittr %>% %T>%
#' @import dplyr ggplot2 RColorBrewer
plotSubcellular_sub <- function (
    df2 = NULL, custom_prefix = NULL, fn_prefix = NULL, fn_suffix = NULL, 
    df = NULL, id = "gene", col_group = NULL, col_order = NULL, 
    col_cluster = "cluster", levels_group = NULL, col_fraction = NULL, 
    levels_fraction = NULL, col_subtype = NULL, levels_subtype = NULL, 
    data_type = "Trend", pat = "Trend_[ONZ]_.*nclust\\d+\\.txt$", ext = "txt", 
    tie_method = "none", qt = .5, complete_cases = FALSE, panel_ids = NULL, 
    filepath = NULL, label_scheme_sub = NULL, theme = NULL, make_plot = TRUE, 
    dots = NULL)
{
  ## (1) Load input data
  col_x <- "group"
  col_y <- "loc_score"
  cl_id <- as.integer(sub(paste0(".*_nclust(\\d+)[^\\d]*\\.", ext), "\\1", df2))

  fn_pre <- paste0(custom_prefix, fn_prefix, "_nclust", cl_id)
  out_nm <- paste0(fn_pre, ".", fn_suffix)

  # [x] Column format from 'anal_prnTrend':
  #   gene, cluster, group, log2FC, log10Int, purity, entropy, loc_score, 
  #   sub_location | 
  #   Fraction, Type	species, col_group, col_order, col_fraction, col_subtype	
  #   
  # [x] Column format from svmSubcell: 
  #   gene, cluster, compartment, purity, entropy, loc_score | Type 

  if (!file.exists(src_path <- file.path(filepath, df2))) {
    stop("File not found: ", src_path)
  }
  
  if (is.null(dim(df <- readr::read_tsv(src_path)))) {
    stop("File contains not data: ", src_path)
  }
  
  message(paste("File loaded:", src_path))
  
  ## (2) Non-Trend data: compile the master X-axis column, 'col_x'
  #  Changed 'col_fraction' to 'svm'!!!
  if (!(data_type == "Trend" || col_x %in% names(df))) {
    levels_fraction <- sort(unique(df[["compartment"]]))
    col_fraction <- "svm"

    # Keep the order of 'paste0' consistent: Fraction_Subtype
    if (length(col_subtype)) {
      if (!length(levels_subtype)) {
        levels_subtype <- sort(unique(df[[col_subtype]]))
      }
      
      df <- df |>
        dplyr::mutate(!!col_subtype := factor(!!rlang::sym(col_subtype), 
                                              levels = levels_subtype)) |>
        dplyr::arrange_at(c(col_fraction, col_subtype)) |> 
        dplyr::mutate(!!col_x := paste0(
          .data[[col_fraction]], "_", .data[[col_subtype]]), )
    } else {
      df <- df |>
        dplyr::mutate(!!col_x := .data[[col_fraction]], )
    }
  }

  ## (3) Define the order the X-axis
  if (!length(levels_group)) {
    if (data_type == "Trend") {
      levels_group <- local({
        col_group <- df[["col_group"]][1]
        col_order <- df[["col_order"]][1]
        
        label_scheme_sub |>
          dplyr::arrange(!!rlang::sym(col_order)) |>
          dplyr::pull(!!rlang::sym(col_group)) |> 
          unique()
      })
    } else {
      if (length(levels_subtype)) {
        grid <- expand.grid(sub = levels_fraction, frac = levels_subtype)
        grid <- grid[with(grid, order(sub, frac)), ]
        levels_group <- paste0(grid[, 1], "_", grid[, 2])
      } else {
        levels_group <- levels_fraction
      }
    }
  }
  
  df <- df |>
    dplyr::mutate(
      !!col_x := factor(!!rlang::sym(col_x), levels = levels_group), 
      !!col_cluster := as.factor(as.integer(!!rlang::sym(col_cluster))), )

  # (x) 'anal_prnTrend.rda' overwritten with the latest call ->
  #     embedded 'col_group' etc. into the trend analysis output
  # (x) Consider to call from corresponding Protein_Trend_Z_nclust[...].rds ...
  
  if (data_type == "Trend") {
    local({
      levs_df  <- levels(df[[col_x]])
      mis_levs <- levs_df[!levs_df %in% levels_group]
      
      if (length(mis_levs)) {
        if (length(mis_levs) > 12L) {
          mis_levs <- c(mis_levs[1:12], "...")
        }

        stop(
          "\n--- Mismatches in data levels ---\n\n", 
          "Levels in `", df2, "`:\n", paste(mis_levs, collapse = ", "), "\n\n", 
          "Levels by `levels_group:\n", paste(levels_group, collapse = ", "), 
          "\n\n", 
          "??? Check for consistency in the setting of ", 
          "`anal_prnTrend(col_group = ...)` ", 
          "and `plot_prnTrend(col_group = ...)` for file `", df2, "`.")
      }
    })

    # col_fraction: for aesthetics: same col_fraction, same color
    if (!length(col_fraction)) col_fraction <- df[["col_fraction"]][1] 
    
    # col_subtype: for separate median statistics
    if (!length(col_subtype)) col_subtype <- df[["col_subtype"]][1]
  }
  
  ## (4) Data preprocessing
  lang_dots     <- dots[unlist(lapply(dots, is.language))]
  filter2_dots  <- lang_dots[grepl("^filter2_", names(lang_dots))]
  arrange2_dots <- lang_dots[grepl("^arrange2_", names(lang_dots))]
  dots <- dots[!dots %in% c(filter2_dots, arrange2_dots)]
  
  if (length(dots) && any(grepl("^filter_", names(dots)))) {
    stop("Primary `filter_` depreciated; use secondary `filter2_`.")
  }
  
  if (complete_cases) {
    df <- df[complete.cases(df), ]
  }
  
  if (data_type == "Trend") {
    df <- df |>
      dplyr::filter(!!rlang::sym(col_x) %in% levels_group) |>
      filters_in_call(!!!filter2_dots) |>
      arrangers_in_call(!!!arrange2_dots) |>
      dplyr::mutate(!!col_x := factor(.data[[col_x]], levels = levels_group))
  }
  
  ## (3.2) Calculate median description of localization scores
  # User defined column `group" provides the finest granularity including 
  #  subcellular, sample type, etc.
  target_cols <- c("purity", "entropy", col_y)
  
  df_med <- lapply(split(df, df[[col_cluster]]), function(dfx) {
    dfx |>
      dplyr::group_by(.data[[col_x]]) |>
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::any_of(target_cols),
          .fns  = ~ as.numeric(quantile(.x, probs = qt, na.rm = TRUE)),
          .names = "{.col}"
        ),
        .groups = "drop"
      )
  }) |>
    dplyr::bind_rows(.id = col_cluster) |>
    dplyr::mutate(
      !!col_cluster := as.integer(.data[[col_cluster]]),
      !!col_cluster := factor(
        .data[[col_cluster]], levels = sort(unique(.data[[col_cluster]]))), )
  
  df_med <- df_med |>
    dplyr::left_join(
      df |>
        dplyr::select(
          dplyr::one_of(c(col_cluster, col_x, col_fraction, col_subtype))) |>
        unique(), 
      by = c(col_cluster, col_x))
  
  readr::write_tsv(
    df_med, file.path(filepath, paste0(fn_pre, "_", data_type, "_med.tsv")))

  ### (4) Visualization
  # n_clust <- length(unique(df[[col_cluster]]))
  
  if (make_plot) {
    try(
      make_subcell_plot(
        df_med = df_med, col_cluster = col_cluster, 
        col_x = col_x, col_y = col_y, col_fraction = col_fraction, 
        panel_ids = panel_ids, cl_id = cl_id, 
        filepath = filepath, out_nm = out_nm, !!!dots)
    )
  }
  
  list(
    df_med = df_med, 
    tie_method = tie_method, 
    col_fraction = col_fraction, 
    col_subtype = col_subtype, 
    col_cluster = col_cluster,
    col_x = col_x, 
    col_y = col_y)
}


#' Helper of making plots of subcellular scores.
#' 
#' @param df_med A data frame of median summarization.
#' @param col_cluster The column key of \code{cluster}.
#' @param col_x The column key of X-axis.
#' @param col_y The column key of Y-axis.
#' @param levels_group The levels of groups.
#' @param col_fraction The column key of subcellular fractions.
#' @param levels_fraction The levels of fractions.
#' @param panel_ids The subset of panel_ids for visualization.
#' @param cl_id A cluster ID.
#' @param filepath An output filepath.
#' @param out_nm An output filename.
#' @param theme The theme of plot.
#' @param ... Additional parameters for aesthetics.
make_subcell_plot <- function (
    df_med, n_clust, col_cluster = "cluster", 
    col_x = "group", col_y = "loc_score", 
    levels_group = NULL, col_fraction = "Fraction", levels_fraction = NULL, 
    panel_ids = NULL, cl_id = NULL, filepath = NULL, out_nm = NULL, 
    theme = NULL, ...)
{
  dots <- rlang::enexprs(...)
  
  ## Subset data by panel_ids and levels_group
  n_panels  <- length(panel_ids)
  
  if (!n_panels) {
    panel_ids <- unique(df_med[[col_cluster]])
    n_panels  <- length(panel_ids)
  }
  
  df_med <- df_med |>
    dplyr::filter(.data[[col_cluster]] %in% panel_ids) |>
    dplyr::mutate(
      !!col_cluster := factor(!!rlang::sym(col_cluster), levels = panel_ids)
    )
  
  if (length(levels_group)) {
    df_med <- df_med |>
      dplyr::filter(.data[[col_x]] %in% levels_group) |>
      dplyr::mutate(
        !!col_x := factor(!!rlang::sym(col_x), levels = levels_group)
      )
  }
 
  ## Set up plot parameters
  ymin    <- eval(dots$ymin, envir = rlang::caller_env())
  ymax    <- eval(dots$ymax, envir = rlang::caller_env())
  ybreaks <- eval(dots$ybreaks, envir = rlang::caller_env())
  ncol    <- dots$ncol
  nrow    <- dots$nrow
  width   <- dots$width
  height  <- dots$height
  color   <- dots$color
  alpha   <- dots$alpha
  
  if (is.null(ymin)) ymin <- 0
  if (is.null(ymax)) ymax <- 1
  if (is.null(ybreaks)) ybreaks <- .2
  if (is.null(ncol)) ncol <- 1
  if (is.null(nrow)) nrow <- 2
  if (is.null(color)) color <- "#2c7fb8"
  if (is.null(alpha)) alpha <- .25
  
  if (is.null(ncol)) {
    if (is.null(nrow)) {
      ncol <- 4L
      nrow <- ceiling(n_panels / ncol)
    } else {
      ncol <- ceiling(n_panels / nrow)
    }
  } else {
    nrow <- ceiling(n_panels / ncol)
  }
  
  dots$ymin <- NULL
  dots$ymax <- NULL
  dots$ybreaks <- NULL
  dots$ncol <- NULL
  dots$nrow <- NULL
  dots$color <- NULL
  dots$alpha <- NULL
  
  
  if (is.null(width)) width <- n_panels * 8 / nrow + 2
  if (is.null(height)) height <- 8 * nrow
  
  dots$width <- NULL
  dots$height <- NULL
  
  # Handle multiple sub_type later...
  # e.g., df[[col_x]]: Veh_CP, TG_CP, df[[col_fraction]]: CP, CP
  # This allow the same fill color for 
  # different df[[col_x]] at the same df[[col_fraction]]
  # if (!length(col_fraction)) col_fraction <- col_x
  
  proteoq_trend_theme <- theme_bw() + 
    theme(
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
      # panel.background = element_rect(fill = '#0868ac', colour = 'red'),
      panel.background = element_rect(fill = "#f5f5f5", colour = NA), 
      
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
  
  if (is.null(theme)) {
    theme <- proteoq_trend_theme
  }

  p <- ggplot(
    df_med, aes(x = !!rlang::sym(col_x), y = !!rlang::sym(col_y), 
                fill = !!rlang::sym(col_fraction))) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(vars(!!rlang::sym(col_cluster)), ncol = ncol, 
               scales = "free_y") +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    scale_fill_viridis_d(option = "plasma", end = 0.85) + 
    scale_y_continuous(
      labels = scales::label_number(accuracy = .1),
      breaks = scales::breaks_width(.2)
    ) +
    labs(
      x = NULL,
      y = "Score",
      fill = col_fraction
    ) +
    # theme_classic(base_size = 11) +
    theme
  
  ggsave_dots <- 
    set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
  
  rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(out_nm)),
                    # plot = p, 
                    width = width, 
                    height = height, 
                    !!!ggsave_dots)) |>
    rlang::eval_tidy()
}


#' Helper of \link{classify_trends}.
#' 
#' @param col_x The column name of X values.
#' @param col_y The column name of Y values for summarization.
#' @param col_cluster The column key of cluster nodes.
#' @param tie_method A method handling ties in localization scores.
#' @param df_med A median description of localization scores by modules.
#' @inheritParams anal_prnTrend
hclassify_trends <- function (
    df_med = NULL, col_subtype = NULL, col_fraction = NULL, 
    col_group = NULL, col_order = NULL, col_cluster = "cluster", 
    col_x = NULL, col_y = NULL, tie_method = "none")
  
{
  if (length(col_subtype) && col_subtype %in% names(df_med)) {
    res <- lapply(split(df_med, df_med[[col_subtype]]), function (dfx) {
      tryCatch(
        classify_trends(df_med = dfx, 
                        tie_method = tie_method, 
                        col_fraction = col_fraction, 
                        col_subtype = col_subtype,
                        col_cluster = col_cluster,
                        col_y = col_y),
        error = function(e) NULL)
    }) |>
      dplyr::bind_rows(.id = col_subtype) |>
      dplyr::ungroup() |>
      dplyr::arrange(!!rlang::sym(col_fraction), !!rlang::sym(col_y))
  } else {
    res <- tryCatch(
      classify_trends(df_med = df_med, 
                      tie_method = tie_method, 
                      col_fraction = col_fraction, 
                      col_subtype = col_subtype, 
                      col_cluster = col_cluster,
                      col_y = col_y),
      error = function(e) NULL)
  }
  
  # Currently use the lowest score across all cell states of `col_subtype`
  if (FALSE) {
    res <- res |>
      dplyr::group_by(!!rlang::sym(col_fraction)) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::select(-dplyr::one_of(col_subtype))
  }
  
  list(res = res, 
       col_group = col_group,
       col_order = col_order,
       col_fraction = col_fraction, 
       col_subtype = col_subtype, 
       col_x = col_x, 
       col_y = col_y)
}


#' Classify trends for one value of \code{col_subtype}.
#' 
#' @param col_y The column name of Y values for summarization.
#' @param col_cluster The column key of cluster nodes.
#' @param tie_method A method handling ties in localization scores.
#' @param df_med A median description of localization scores by modules.
#' @inheritParams anal_prnTrend
classify_trends <- function (df_med = NULL, tie_method = "none", 
                             col_y = "score", col_cluster = "cluster", 
                             col_fraction = NULL, col_subtype = NULL)
{
  subcell_levs <- levels(df_med[[col_fraction]])
  
  ## 1. Assign the primary subcellular fraction for each cluster
  df_pri <- df_med |>
    dplyr::select(dplyr::one_of(
      c(col_cluster, "purity", "entropy", col_y, col_fraction)
    )) |>
    dplyr::arrange(!!rlang::sym(col_cluster), -!!rlang::sym(col_y))

  # Keep the first combination among multiple col_subtype of Ctrl, Treatment...
  combs <- interaction(df_pri[[col_cluster]], df_pri[[col_fraction]], 
                       drop = TRUE)
  df_pri <- lapply(
    split(df_pri, combs), 
    function(x) x[1, , drop = FALSE]
  ) |> 
    dplyr::bind_rows() |>
    dplyr::arrange(!!rlang::sym(col_cluster), -!!rlang::sym(col_y))
  
  # [x] May assign to both fractions if < .01, but not contribute to the thresholds
  # [x] Or count the percent of missing values per group and assign the cluster to 
  #  the most complete groups; must be >= 50% values...
  is_unambiguous <- sapply(split(df_pri, df_pri[[col_cluster]]), function (x) {
    nrow <- nrow(x)
    
    # Unambiguous assignment of the primary fraction for a cluster module
    if (nrow == 1L) {
      return(TRUE)
    }
    
    fracs <- x[[col_fraction]]
    n_fracs <- length(fracs)
    scos <- x[[col_y]]
    # scos[[1]] - scos[[2]] >= .01
    scos[[1]] - scos[[2]] >= .025 * exp(-n_fracs / 3)
  })
  
  # (1.a) Fractions that never peak in localization scores are dropped here
  df_pri <- df_pri |>
    dplyr::group_by(!!rlang::sym(col_cluster))
  
  df_sec <- df_pri |>
    dplyr::slice(-1L) |>
    dplyr::ungroup()
  
  df_fst <- df_pri |>
    dplyr::slice(1L) |>
    dplyr::ungroup()
  
  # (1.b) Clusters that are tied in the scores between the primary and the 
  # secondary fractions are dropped here
  df_fst[["subcellular..."]] <- df_fst[[col_fraction]]
  df_fst[[col_fraction]][!is_unambiguous] <- NA_character_
  
  # Find the probable score cut-offs for each subcellular group;
  df_pri_min <- df_fst |>
    dplyr::arrange(!!rlang::sym(col_fraction), !!rlang::sym(col_y)) |>
    dplyr::group_by(!!rlang::sym(col_fraction)) |>
    dplyr::slice(1L) |>
    dplyr::filter(!is.na(!!rlang::sym(col_fraction)))
  
  if (tie_method != "none") { # "all" or "first"
    ## Add back ties ("1.b" fractions)
    df_tie <- df_fst[!is_unambiguous, ] |>
      dplyr::filter(
        !.data[["subcellular..."]] %in% df_pri_min[[col_fraction]])
    df_tie[[col_fraction]] <- df_tie[["subcellular..."]]
    
    # Use the minimum score as the cut-off (not maximum)
    df_tie_min <- df_tie |>
      dplyr::arrange(!!rlang::sym(col_fraction), !!rlang::sym(col_y)) |>
      dplyr::group_by(!!rlang::sym(col_fraction)) |>
      dplyr::slice(1L) 
    
    df_pri_min <- df_pri_min |>
      dplyr::bind_rows(df_tie_min) |>
      dplyr::arrange(!!rlang::sym(col_fraction))
    df_pri_min[["subcellular..."]] <- NULL
    
    ## Add back fractions "1.a" by finding the nearest neighbors
    mis_levs <- subcell_levs[!subcell_levs %in% df_pri_min[[col_fraction]]]

    if (length(mis_levs)) {
      df_secs <- split(df_sec, df_sec[[col_cluster]])
      
      df_no_apex <- lapply(mis_levs, function (l) {
        mt_rownums_oks <- sapply(df_secs, function (x) {
          which(x[[col_fraction]] == l)
        })
        
        rns_min <- which(mt_rownums_oks == min(mt_rownums_oks, na.rm = TRUE))
        df_secs_min <- df_secs[rns_min]
        mt_rownums_oks_min <- mt_rownums_oks[rns_min]
        
        res_sec <- vector("list", length(rns_min))
        for (i in seq_along(rns_min)) {
          dfx <- df_secs_min[[i]]
          rnx <- mt_rownums_oks_min[[i]]
          res_sec[[i]] <- dfx[rnx, ]
        }
        
        cl_max <- which.max(lapply(res_sec, function (x) (x[[col_y]])))
        res_sec[[cl_max]]
      }) |>
        dplyr::bind_rows()
      
      df_pri_min <- df_pri_min |> 
        dplyr::bind_rows(df_no_apex) |>
        dplyr::arrange(.data[[col_fraction]])
    }
  }
  
  return(df_pri_min)
}


#' Bar plots of median subcellular scores.
#'
#' \code{plot_prnSubcellular} plots the subcellular purity of protein
#' expressions from \code{\link{anal_prnTrend}} or the probabilities from
#' \link{svmSubcell}.
#'
#' @param dat_dir The working directory.
#' @param levels_group The levels of data groups.
#' @param levels_fraction The levels of subcellular fractions.
#' @param col_subtye For uses with subcellular data. The column key in metadata
#'   (\code{expt_smry.xlsx}) linking to sample sub-types, e.g., \code{Control,
#'   Treated}.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param pat The pattern of input filenames.
#' @param ext The extension of input filenames.
#' @param data_type The type of input data.
#' @param qt Not used. The quantile of localization score for thresholding
#'   subcellular compartment specificity. The default is to use median.
#' @param make_plot Logical; to bypass the plotting at FALSE.
#' @inheritParams anal_prnTrend
#' @inheritParams plot_prnTrend
#' @examples \donttest{
#' # c("Veh_CP", "TG_CP", "Veh_NP", "TG_NP", "Veh_ChA", "TG_ChA")
#' levels_group <- load_ls_group(dat_dir = dat_dir) |>
#'   dplyr::arrange_at("Order") |>
#'   dplyr::pull("Group") |>
#'   unique()
#' 
#' levels_fraction <- c("CP", "NP", "ChA")
#' levels_subtype <- c("Veh", "TG")
#' 
#' # Automatic deduction of 'levels_group' or overruled with custom orders
#' plot_prnSubcellular(
#'   # levels_group = levels_group, 
#'   col_fraction = "Fraction", levels_fraction = levels_fraction, 
#'   col_subtype = "Type", levels_subtype = levels_subtype, 
#'   ncol = 4, limitsize = FALSE, 
#' )
#' }
#' @export
plot_prnSubcellular <- function (
    dat_dir = NULL, levels_group = NULL, 
    col_select = NULL, col_order = NULL, 
    col_fraction = NULL, levels_fraction = NULL, 
    col_subtype = NULL, levels_subtype = NULL, tie_method = "none", qt = .5, 
    n_clust = NULL, panel_ids = NULL, 
    pat = "Trend_[ONZ]_.*nclust\\d+\\.txt$", ext = "txt", data_type = "Trend", 
    scale_log2r = TRUE, complete_cases = FALSE, 
    impute_na = FALSE, df2 = NULL, filepath = NULL, filename = NULL, 
    make_plot = TRUE, theme = NULL, ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")

  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  id <- tryCatch(match_call_arg(normPSM, group_pep_by), error = function(e) NA)
  
  if (is.na(id)) {
    id <- tryCatch(
      match_call_arg(makeProtDIANN, group_pep_by), 
      error = function(e) NA)
  }
  
  if (is.na(id)) {
    id <- "gene"
  }
  
  if (length(id) != 1L || !id %in% c("prot_acc", "gene")) {
    stop("Value of 'id' need to be one of 'prot_acc', 'gene'.")
  }
  
  ## Argument with single option and NULL default (quotation marks or not)
  col_select  <- rlang::enexpr(col_select)
  col_order   <- rlang::enexpr(col_order)
  col_fraction <- rlang::enexpr(col_fraction)
  col_subtype <- rlang::enexpr(col_subtype)
  filepath   <- rlang::enexpr(filepath)
  filename   <- rlang::enexpr(filename)
  df2 <- rlang::enexpr(df2)
  
  # NULL or symbol -> length 0 or 1
  if (!is.character(col_select)) col_select <- as.character(col_select)
  if (!is.character(col_order)) col_order <- as.character(col_order)
  if (!is.character(col_fraction)) col_fraction <- as.character(col_fraction)
  if (!is.character(col_subtype)) col_subtype <- as.character(col_subtype)
  if (!is.character(filepath)) filepath <- as.character(filepath)
  if (!is.character(filename)) filename <- as.character(filename)
  if (!is.character(df2)) df2 <- as.character(df2)

  tie_method <- if (n_subcell_levs <- length(levels_fraction)) {
    if (n_subcell_levs <= 3) {
      "none"
    }
    else if (n_subcell_levs <= 6) {
      "first"
    } else {
      "all"
    }
  } else {
    "none"
  }
  
  dir.create(file.path(filepath, "log"), recursive = TRUE, showWarnings = FALSE)
  
  reload_expts()
  
  info_anal(id = id, 
            col_select = col_select, 
            col_group = NULL, 
            col_order = col_order,
            df2 = df2, 
            filepath = filepath, 
            filename = filename,
            
            dat_dir = dat_dir, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = NULL, 
            anal_type = "Subcellular_plot")(
              levels_group = levels_group, 
              col_fraction = col_fraction, 
              levels_fraction = levels_fraction, 
              col_subtype = col_subtype, 
              levels_subtype = levels_subtype, 
              tie_method = tie_method, 
              n_clust = n_clust, 
              pat = pat, ext = ext, data_type = data_type, 
              qt = qt, panel_ids = panel_ids, 
              make_plot = make_plot, 
              theme = theme, ...)
}


#' SVM classification of subcellular locations.
#'
#' @param species Species in one of human, mouse and rat.
#' @param id Identifier.
#' @param col_y The column key in \code{df_ref} for training labels.
#' @param col_compartment The identifier of subcellular compartment.
#' @param dat_dir The current working directory.
#' @param filepath The filepath of inputs and outputs.
#' @param n_clust The number of clusters.
#' @param u_subtypes The unique levels of sample types. A NULL value indicates
#'   only one subtype.
#' @param levels_fraction The levels of subcellular fractions.
#' @param seed A seed for SVM.
#' @inheritParams anal_prnTrend col_subtype
svmSubcell <- function (
    species = "human", id = "gene", dat_dir = NULL, filepath = NULL, 
    label_scheme_sub = NULL, 
    n_clust = NULL, col_compartment = "compartment", col_y = "location", 
    col_group = NULL, 
    col_subtype = NULL, levels_subtype = NULL, u_subtypes = NULL, 
    col_fraction = NULL, 
    levels_fraction = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), 
    seed = 42, ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(label_scheme_sub)) label_scheme_sub <- load_ls_group(dat_dir)
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  ## (1) Obtain expression data and reference data
  tempdata <- prepSubcell(
    species = species, id = id, dat_dir = dat_dir, filepath = filepath, 
    label_scheme_sub = label_scheme_sub, 
    n_clust = n_clust, col_group = col_group, 
    col_subtype = col_subtype, levels_subtype = levels_subtype, 
    u_subtypes = u_subtypes, 
    col_fraction = col_fraction, levels_fraction = levels_fraction)
  
  df_exprs <- tempdata[["df_exprs"]]
  df_ref   <- tempdata[["df_ref"]]

  saveRDS(df_exprs, file.path(filepath, "df_exprs.rds"))
  saveRDS(df_ref, file.path(filepath, "df_ref.rds"))
  rm(list = "tempdata")
  
  if (length(u_subtypes) > 1L) {
    if (!identical(names(df_exprs), u_subtypes)) {
      stop("Developer: mismatches between expression data and sample types.")
    }
  }

  ## (2) Results by sub_types
  res <- vector("list", length(df_exprs))
  for (i in seq_along(res)) {
    res[[i]] <- train_loc_svm(
      df_expr = df_exprs[[i]], sub_type = u_subtypes[[i]], 
      df_ref = df_ref, id = id, col_y = col_y, col_subtype = col_subtype, 
      col_compartment = col_compartment, filepath = filepath, seed = seed)
  }
  
  names(res) <- names(df_exprs)
  
  list(prob = res, expr = df_exprs)
}


#' Train localization SVM.
#'
#' @param df_expr A matrix of expression data, typically normalized intensities
#'   that summed to one across fractions for each protein.
#' @param df_ref A reference data frame.
#' @param id An identifier.
#' @param col_y The column key in \code{df_ref} for training labels.
#' @param col_compartment The identifier of subcellular compartment.
#' @param sub_type A label of data subtype, e.g., "Control".
#' @param seed A seed for SVM.
#' @return Two lists of (1) SVM probablity matrix and (2) expression data.
train_loc_svm <- function (
    df_expr, sub_type = NULL, df_ref, id = "gene", col_y = "location", 
    col_subtype = NULL, col_compartment = "compartment", filepath = NULL, 
    seed = 42L, ...)
{
  Xmat <- df_expr[df_ref[[id]], ]
  Xmat <- Xmat[match(rownames(Xmat), df_ref[[id]]), ]
  
  if (!identical(rownames(Xmat), df_ref[[id]])) {
    stop("Develepor: mismatches in identifiers of ", id, ".")
  }
  
  if (!col_y %in% colnames(df_ref)) {
    stop("The column of Y labels", col_y, " not found in reference data.")
  }
  
  y_factor <- factor(df_ref[[col_y]])
  df_train <- data.frame(Xmat, check.names = FALSE)
  df_train[[col_y]] <- y_factor
  class_names <- levels(y_factor)
  
  ## (2) Training
  set.seed(seed)
  svm_model <- e1071::svm(
    formula = stats::reformulate(".", response = col_y),
    data = df_train,
    kernel = "radial", 
    probability = TRUE, 
    scale = TRUE
  )
  
  svm_preds <- predict(svm_model, df_expr, probability = TRUE)
  prob_mat  <- attr(svm_preds, "probabilities")
  prob_mat  <- prob_mat[, class_names]
  
  if (FALSE) {
    svm_preds_ <- predict(svm_model, df_train, probability = TRUE)
    prob_mat_ <- attr(svm_preds_, "probabilities")
    prob_mat_ <- prob_mat_[, class_names]
    
    data.frame(
      True      = as.vector(table(y_factor)),
      Predicted = as.vector(table(factor(svm_preds_, levels = class_names))),
      row.names = class_names
    )
    
    # Overall training accuracy
    acc <- mean(df_train$location == svm_preds_)
    cat("Training Accuracy:", round(acc * 100, 2), "%\n")
    
    # Confusion matrix
    table(True = df_train$location, Predicted = svm_preds_)
  }
  
  ## Outputs
  df_pred    <- tibble::tibble(!!id := names(svm_preds), svm = svm_preds)
  entropies  <- -rowSums(prob_mat * log10(prob_mat))
  scores_loc <- prob_mat * (1 - entropies)
  entropies  <- tibble::tibble(!!id := names(entropies), entropy = entropies)
  entropies  <- entropies |> dplyr::left_join(df_pred, by = id)
  
  scores_loc <- scores_loc |>
    data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = id) |> 
    tidyr::pivot_longer(
      -c(!!rlang::sym(id)), names_to = col_compartment, values_to = "loc_score")
  
  prob_mat <- prob_mat |>
    data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = id) |> 
    tidyr::pivot_longer(
      -c(!!rlang::sym(id)), names_to = col_compartment, values_to = "purity")
  
  ans <- prob_mat |>
    dplyr::left_join(entropies, by = id) |>
    dplyr::left_join(scores_loc, by = c(id, col_compartment))
}


#' Prepare subcellular data.
#' 
#' @param species Species in one of human, mouse and rat.
#' @param id Identifier.
#' @param dat_dir The current working directory.
#' @param filepath The filepath of inputs and outputs.
#' @param key_col The column key where data to be used.
#' @param levels_fraction The levels of subcellular fractions.
#' @param u_subtypes The unique levels of sample types. A NULL value indicates
#'   only one subtype.
#' @inheritParams anal_prnTrend col_subtype
prepSubcell <- function (
    species = "human", id = "gene", dat_dir = NULL, filepath = NULL, 
    label_scheme_sub = NULL, 
    n_clust = NULL, key_col = "log10Int", col_group = NULL, 
    col_subtype = NULL, levels_subtype = NULL, u_subtypes = NULL, 
    col_fraction = NULL, 
    levels_fraction = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(label_scheme_sub)) label_scheme_sub <- load_ls_group(dat_dir)
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  abbr_species <- switch(
    species, 
    human = "hs",
    mouse = "mm",
    rat = "rn",
    stop("species: ", species, " is not in one of `human`, `mouse` or `rat`.")
  )

  df_ref <- local({
    file_db <- paste0("subcell_", abbr_species)
    data(package = "proteoQ", list = file_db, envir = environment())
    df_ref <- get(file_db)
    df_ref |> dplyr::mutate(location = sub("ERGIC", "ER", location))
  })

  if (length(levels_fraction) && 
      all(levels_fraction %in% c("NP", "CP", "ChA"))) {
    df_ref$location <- df_ref$location_slim3
  }
  
  # Find the file names of secondary 'Trend' inputs
  # (dummies: 'scale_log2r' and 'impute_na')
  scale_log2r <- TRUE
  impute_na   <- FALSE
  
  filenames_df2 <- tryCatch(
    find_trend_df2(
      df2 = NULL, n_clust = n_clust, scale_log2r = scale_log2r, 
      impute_na = impute_na, filepath = filepath), 
    error = function (e) NULL)
  
  if (is.null(filenames_df2)) {
    df <- tryCatch(
      readr::read_tsv(file.path(filepath, "Protein_Intensity_summary.tsv")), 
      error = function (e) NULL)
    
    if (is.null(df)) {
      stop("No inputs under", filepath, "\n", "Run 'anal_prnTrend' first.")
    }
    
    df <- df |>
      tidyr::pivot_longer(-!!id, names_to = "group", values_to = key_col) |>
      dplyr::mutate(!!key_col := log10(!!rlang::sym(key_col)))
  } else {
    df <- readr::read_tsv(file.path(filepath, filenames_df2[[1]]))
  }

  vals_g <- df[["group"]]
  nms_df <- names(df)
  
  if (no_subtype <- !length(col_subtype)) col_subtype <- "sample_type_unknown"
  
  pat_frac <- paste0("(", paste(levels_fraction, collapse = "|"), ")$")
  pat_subtype  <- paste0("^(", paste(levels_subtype, collapse = "|"), ")")
  
  if (!col_subtype %in% nms_df) {
    df[[col_subtype]] <- ifelse(
      grepl(pat_frac, vals_g), sub(pat_frac, "", vals_g), vals_g)
    
    df[[col_subtype]] <- ifelse(
      grepl(pat_subtype, df[[col_subtype]]),
      sub(paste0( pat_subtype, ".*"), "\\1", df[[col_subtype]]),
      NA_character_
    )
    
    # Value: "" if 'no_subtype' 
    u_subtypes <- unique(df[[col_subtype]])
  }
  
  # Should not incur
  # if (!length(col_fraction)) col_fraction <- "fraction..."

  if (!col_fraction %in% nms_df) {
    df[[col_fraction]] <- ifelse(
      grepl(pat_frac, vals_g),
      sub(paste0(".*", pat_frac), "\\1", vals_g),
      NA_character_)
  }

  # rm(list = c("vals_g", "nms_df", "pat_frac", "pat_subtype"))
  
  dfs <- split(df, df[[col_subtype]])
  
  # Match the order of levels_subtype
  if (length(levels_subtype) > 1L) {
    dfs <- dfs[levels_subtype]
  }
  
  # Normalized intensity (percentages)
  sc_lookup <- label_scheme_sub |>
    dplyr::select(dplyr::one_of(col_group, col_fraction, col_subtype)) |>
    unique() |>
    dplyr::rename(Sample_ID := col_group)
  
  df_exprs <- vector("list", length(dfs))
  for (i in seq_along(dfs)) {
    df_expr <- make_subcell_expr(
      df = dfs[[i]], id = id, key_col = key_col, 
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      scale = TRUE)
    
    # Change from values under col_group to values under col_fraction
    colnames(df_expr) <- sc_lookup[[col_fraction]][
      match(colnames(df_expr), sc_lookup[["Sample_ID"]])
    ]
    
    df_expr <- df_expr[, levels_fraction]
    df_exprs[[i]] <- df_expr
  }
  
  # Check for identical rownames
  local({
    if ((n_subtypes <- length(u_subtypes)) > 1L) {
      for (i in (n_subtypes - 1L)) {
        if (!identical(rownames(df_exprs[[i]]), rownames(df_exprs[[i+1]]))) {
          stop("Developer: mismatched rownames.")
        }
      }
    }
  })

  df_ref <- df_ref[df_ref[[id]] %in% rownames(df_exprs[[1]]), ]
  
  if (!nrow(df_ref)) {
    stop("No reference data retained after ID matching. \n", 
         "Check the correctionness of 'species = ", species, "'.")
  }
  
  df_ref <- df_ref[!is.na(df_ref[[id]]), ]
  
  list(
    df_exprs = df_exprs, 
    df_ref = df_ref, 
    abbr_species = abbr_species, 
    filenames_df2 = filenames_df2
  )
}


#' Helper of Preparing the matrix of expression data.
#' 
#' @inheritParams make_subcell_expr
hmake_subcell_expr <- function (
    df = NULL, label_scheme_sub = NULL, id = "gene", key_col = "log10Int", 
    col_subtype = NULL, levels_subtype = NULL, 
    col_fraction = NULL, levels_fraction = NULL)
{
  dfs <- split(df, df[[col_subtype]])
  label_schemes <- split(label_scheme_sub, label_scheme_sub[[col_subtype]])
  
  # Match the order of levels_subtype
  if (length(levels_subtype) > 1L) {
    dfs <- dfs[levels_subtype]
    label_schemes <- label_schemes[levels_subtype]
  }
  
  # Normalized intensity (percentages)
  df_exprs <- vector("list", length(dfs))
  for (i in seq_along(dfs)) {
    df_exprs[[i]] <- make_subcell_expr(
      df = dfs[[i]], 
      id = id, key_col = key_col, 
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      scale = TRUE)
  }
  
  local({
    u_subtypes <- unique(df[["col_subtype"]])
    
    if ((n_subtypes <- length(u_subtypes)) > 1L) {
      for (i in (n_subtypes - 1L)) {
        if (!identical(rownames(df_exprs[[i]]), rownames(df_exprs[[i+1]]))) {
          stop("Developer: mismatched rownames.")
        }
      }
    }
  })
  
  df_exprs
}


#' Prepare the matrix of expression data.
#' 
#' One subtype at a time.
#' 
#' @param df An input data frame.
#' @param id An identifier.
#' @param key_col The column key where values will be used for preparing
#'   expression data.
#' @param key_group The column key where group information will be found for
#'   preparing expression data.
#' @param levels_fraction The levels of subcellular fractions.
#' @param scale Logical; to scale values by rowSums or not.
make_subcell_expr <- function (
    df = NULL, id = "gene", key_col = "log10Int", # label_scheme_sub = NULL, 
    key_group = "group", col_fraction = NULL, 
    levels_fraction = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), 
    scale = TRUE)
{
  if (!key_group %in% colnames(df)) {
    print(head(df))
    stop("Developer: column 'group' not found.")
  }
  
  df_expr <- df[, c(id, key_group, key_col)]
  
  if (key_col == "log10Int") {
    df_expr <- df_expr |> 
      dplyr::mutate(Int = 10^(!!rlang::sym(key_col))) |>
      dplyr::select(-dplyr::one_of(key_col))
    
    key_col <- "Int"
  }
  
  local({
    u_grps <- unique(df_expr[[key_group]])
    n_grps <- length(u_grps)
    n_levs <- length(levels_fraction)
    
    if (n_grps > n_levs) {
      stop("The number of groups is different to the number of ", 
           "fraction levels: \n", key_group, ": ", 
           paste(u_grps, collapse = ", "), "\n",
           "levels_fraction", ": ", 
           paste(levels_fraction, collapse = ", "))
    }
  })

  df_expr <- df_expr |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym(key_group), 
      values_from = !!rlang::sym(key_col))
  
  rns <- df_expr[[id]]
  df_expr <- as.matrix(df_expr[, -which(names(df_expr) %in% id)])
  rownames(df_expr) <- rns
  
  # Normalized intensity
  if (scale) {
    df_expr <- df_expr / rowSums(df_expr)
  }
  
  df_expr
}


#' Cluster SVM probabilities for a given sample subtype.
#'
#' @param df_prob A data frame of SVM probabilities (under column
#'   \code{loc_score}).
#' @param dat_dir The working directory.
#' @param filepath A filepath for inputs and outputs.
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param n_clust The number of clusters.
#' @param choice A clustering method.
cluster_svmprobs <- function (
    df_prob = NULL, sub_type = NULL, dat_dir = NULL, filepath = NULL, 
    id = "gene", n_clust = NULL, 
    choice = c("cmeans", "clara", "kmeans", "pam", "fanny"), 
    ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  choice <- match.arg(choice)
  dots <- rlang::enexprs(...)
  
  # Compartment probabilities
  df_comp <- local({
    df_comp <- df_prob[, c(id, "compartment", "loc_score")] |>
      tidyr::pivot_wider(
        id_cols = id, 
        names_from = c("compartment"), 
        values_from = "loc_score",
        # values_fill = 0, 
      )
    
    ids <- df_comp[[id]]
    df_comp <- as.matrix(df_comp[, colnames(df_comp) != id])
    rownames(df_comp) <- ids
    
    df_comp
  })
  
  ans_dots <- find_trend_m(
    dots, choice = choice, n_clust = n_clust, df_mean_log2r = df_comp)
  dots <- ans_dots[["dots"]]
  n_clust <- ans_dots[["n_clust"]]
  rm(list = c("ans_dots"))
  
  # Add cluster numbers
  out <- lapply(
    paste0("svmProb_Z_nclust", n_clust), makeTrendRes, 
    choice = choice, dots = dots, id = id, df_mean_log2r = df_comp)
  out <- lapply(out, dplyr::left_join, df_prob, by = id)
  names(out) <- n_clust

  if (FALSE) {
    out_nms <- sapply(n_clust, function (cl) {
      out_nm <- paste0("Protein_svmProb_Z_nclust", cl)
      out_nm <- paste(c(out_nm, sub_type), collapse = "_")
      out_nm <- sub("_+$", "", out_nm)
      out_nm <- paste0(out_nm, ".tsv")
    })
  }
  
  out
}


#' A wrapper of classifying subcellular locations.
#'
#' @param dat_dir The working directory.
#' @param filepath A filepath for inputs and outputs.
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param species A species in one of \code{"human", "mouse", "rat"} for
#'   identifying the database of reference markers.
#' @param n_clust The number of clusters.
#' @param choice A clustering method.
#' @param min_prot_n_peps The minimum number of detected peptides under a
#'   protein for consideration.
#' @param make_plot Logical; an option to plot clusters of probabilities.
#' @param levels_group The levels of data group for ordering on the X-axis in
#'   visualization.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param seed A random seed.
#' @param ... Additional arguments.
#' @import ggplot2 RColorBrewer
#'
#' @examples \donttest{
#' dat_dir <- "my/working/directory"
#' 
#' # The names of organelle compartments can be either (1) coincide with or (2) 
#' # different to sample fractions.
#' 
#' ## (1) Sample fractions coincide with organelle compartments
#' levels_fraction <- c("CP", "NP", "ChA")
#' levels_subtype <- c("Veh", "TG")
#'
#' # c("Veh_CP", "TG_CP", "Veh_NP", "TG_NP", "Veh_ChA", "TG_ChA")
#' levels_group <- load_ls_group(dat_dir = dat_dir) |>
#'   dplyr::arrange_at("Order") |>
#'   dplyr::pull("Group") |>
#'   unique()
#'
#' # 'levels_group' can be either users' provided or deduced automatically.
#' classify_subcellular(
#'   dat_dir = dat_dir, id = "gene", species = "human", n_clust = 16:36,
#'   col_group = "Group", levels_group = levels_group,
#'   col_fraction = "Fraction", levels_fraction = levels_fraction,
#'   col_subtype = "Type", levels_subtype = levels_subtype)
#'
#' ## (2) Sample fractions are different to organelle compartments
#' data("subcell_mm")
#' compartments <- names(sort(table(subcell_mm$location), decreasing = TRUE))
#' levels_group <- compartments # Assume one 'subtype'
#'
#' classify_subcellular(
#'   dat_dir = dat_dir, id = "gene", species = "mouse", n_clust = 16:36,
#'   col_group = "Group", # for group summarization of raw expression data
#'   levels_group = levels_group, # for the order on the X-axis in visualization
#'   col_fraction = "Fraction", # for preparation of summarized expression data
#'   levels_fraction = # to replace the colnames of summarized expression data
#'     c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"))
#' }
#' @export
classify_subcellular <- function (
    dat_dir = NULL, filepath = NULL, id = "gene", species = "human", 
    n_clust = 16:36, levels_group = NULL, 
    col_fraction = NULL, levels_fraction = NULL, 
    col_subtype = NULL, levels_subtype = NULL, 
    col_group = NULL, col_order = NULL, 
    min_prot_n_peps = 2L, make_plot = TRUE, seed = 42L, ...) 
{
  dots <- rlang::enexprs(...)
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  label_scheme <- load_ls_group(dat_dir)

  if (TRUE) {
    ## Argument with single option and NULL default (quotation marks or not)
    col_group    <- rlang::enexpr(col_group)
    col_order    <- rlang::enexpr(col_order)
    col_fraction <- rlang::enexpr(col_fraction)
    col_subtype  <- rlang::enexpr(col_subtype)
    
    # At NULL default: length(0); at symbol -> length(1)
    if (!is.character(col_group)) col_group <- as.character(col_group)
    if (!is.character(col_order)) col_order <- as.character(col_order)
    if (!is.character(col_fraction)) col_fraction <- as.character(col_fraction)
    if (!is.character(col_subtype)) col_subtype <- as.character(col_subtype)
    
    ## Argument with single, non-NULL option (with or without quotation mark)
    #  (Never NULL since the default is not NULL)
    id <- rlang::as_string(rlang::enexpr(id))
    species <- rlang::as_string(rlang::enexpr(species))
  }

  if (FALSE) {
    ## Special handling of character vectors
    levels_fraction <- eval(rlang::enexpr(levels_fraction))
    levels_subtype <- eval(rlang::enexpr(levels_subtype))
  }
  
  check_subcell(
    levels_fraction = levels_fraction, levels_subtype = levels_subtype, 
    col_subtype = col_subtype)

  ### (1) Obtain aggregated expression data  
  # Need `!!` since `anal_prnTrend` is a UI function
  ans_trend <- anal_prnTrend(
    col_group = !!col_group, 
    col_order = !!col_order,
    col_fraction = !!col_fraction, 
    col_subtype = !!col_subtype, 
    cluster_data = FALSE, # to skip clustering steps
    n_clust = n_clust, 
    filter_by_npep = exprs(prot_n_pep >= !!min_prot_n_peps), )
  # df_mean_log2r <- ans_trend[["log2R"]]
  # df_mean_int   <- ans_trend[["Intensity"]]

  ### (2) Generate an SVM probability matrix
  no_subtype <- !length(col_subtype)
  u_subtypes <- if (no_subtype) NULL else unique(label_scheme[[col_subtype]])
  
  ans_svm <- svmSubcell(
    species = species, id = id, dat_dir = dat_dir, filepath = filepath, 
    col_group = col_group, 
    col_subtype = col_subtype, levels_subtype = levels_subtype, 
    col_compartment = "compartment", 
    col_fraction = col_fraction, levels_fraction = levels_fraction, 
    u_subtypes = u_subtypes, n_clust = n_clust, seed = seed)
  svm_probs <- ans_svm[["prob"]]
  df_exprs  <- ans_svm[["expr"]]

  ### (3) Classify SVM probabilities
  # ans_cl: list by sub_types
  #  list of data frames by n_clust
  # Not `mapply` since u_subtypes can be NULL
  ans_cl <- vector("list", length(svm_probs))
  names(ans_cl) <- u_subtypes
  
  for (i in seq_along(svm_probs)) {
    ans_cl[[i]] <- cluster_svmprobs(
      df_prob = svm_probs[[i]], sub_type = u_subtypes[[i]], dat_dir = dat_dir, 
      filepath = filepath, id = id, n_clust = n_clust, iter.max = 250)
  }
  
  ## (4) Combine data under the same n_clust across all sub_types
  ans_cl <- local({
    out_nms <- paste0("Protein_svmProb_Z_nclust", n_clust, ".tsv")
    
    ans <- vector("list", length(n_clust))
    
    for (i in seq_along(n_clust)) {
      ans[[i]] <- lapply(ans_cl, `[[`, i) |>
        dplyr::bind_rows(.id = if (no_subtype) NULL else col_subtype)
      
      readr::write_tsv(ans[[i]], file.path(filepath, out_nms[[i]]))
    }
    
    ans
  })
  names(ans_cl) <- n_clust
  
  ans_cl <- ans_cl |>
    dplyr::bind_rows(.id = "n_clust")

  ### (5) Classify subcellular compartments by medium statistics
  plot_prnSubcellular(
    col_order = "colnames", # X-axix ordered by compartment
    levels_group = levels_group, 
    col_subtype = !!col_subtype, 
    levels_subtype = levels_subtype,
    col_fraction = !!col_fraction,
    levels_fraction = levels_fraction,
    q = .5, 
    ncol = 6,
    pat = "svmProb_[ONZ]_.*nclust\\d+.*\\.tsv$", 
    ext = "tsv", 
    data_type = "svmProb", 
    limitsize = FALSE,
    make_plot = make_plot, 
  )
  
  ### (5) Optional: Plot clusters of SVM probabilities
  if (make_plot) {
    plotTrend(
      id = id, 
      dat_dir = dat_dir, col_group = NULL, col_order = NULL, 
      label_scheme_sub = NULL, anal_type = "svmProb", 
      pat = "svmProb_[ONZ]_.*nclust\\d+\\_{0,1}.*\\.tsv$", ext = "tsv",  
      n_clust = n_clust, panel_ids = NULL, show_panel_ids = TRUE, 
      group_data_by = "ratio", scale_log2r = TRUE, 
      complete_cases = FALSE, impute_na = FALSE, 
      df2 = NULL, filepath = filepath, 
      filename = "Protein_svmProb_line_Z.png", theme = NULL, 
      ncol = 6, ymin = 0, ymax = 1, limitsize = FALSE, )
  }
}


#' A wrapper of unsupervised classification of subcellular locations.
#'
#' Typically used with ChESS MS data of fractions \code{CP}, \code{NP} and
#' \code{ChA}.
#' 
#' @param dat_dir The working directory.
#' @param filepath A filepath for inputs and outputs.
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param species A species in one of \code{"human", "mouse", "rat"} for
#'   identifying the database of reference markers.
#' @param n_clust The number of clusters.
#' @param choice A clustering method.
#' @param min_prot_n_peps The minimum number of detected peptides under a
#'   protein for consideration.
#' @param make_plot Logical; an option to plot clusters of probabilities.
#' @param levels_group The levels of data group for ordering on the X-axis in
#'   visualization. The argument overrule the orders by \code{col_order}.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param seed A random seed.
#' @param ... Additional arguments.
#' @import ggplot2 RColorBrewer
#'
#' @examples \donttest{
#' dat_dir <- "my/working/directory"
#'
#' classify_subcellular_unsupervised(
#'   dat_dir = dat_dir, species = "human",
#'   n_clust = 16:36, levels_group = NULL,
#'   col_fraction = "Fraction", levels_fraction = c("CP", "NP", "ChA"),
#'   col_subtype = "Type", levels_subtype = c("Veh", "TG"),
#'   col_group = "Group", col_order = "Order")
#' }
#' @export
classify_subcellular_unsupervised <- function (
    dat_dir = NULL, filepath = NULL, col_select = NULL, 
    id = "gene", species = "human", 
    n_clust = 16:36, levels_group = NULL, 
    col_fraction = NULL, levels_fraction = NULL, 
    col_subtype = NULL, levels_subtype = NULL, 
    col_group = NULL, col_order = NULL, 
    min_prot_n_peps = 2L, make_plot = TRUE, seed = 42L, ...) 
{
  dots <- rlang::enexprs(...)
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  label_scheme <- load_ls_group(dat_dir)
  
  if (TRUE) {
    ## Argument with single option and NULL default (quotation marks or not)
    col_select   <- rlang::enexpr(col_select)
    col_order    <- rlang::enexpr(col_order)
    col_group    <- rlang::enexpr(col_group)
    col_fraction <- rlang::enexpr(col_fraction)
    col_subtype  <- rlang::enexpr(col_subtype)

    # At NULL default: length(0); at symbol -> length(1)
    if (!is.character(col_select)) col_select <- as.character(col_select)
    if (!is.character(col_order)) col_order <- as.character(col_order)
    if (!is.character(col_group)) col_group <- as.character(col_group)
    if (!is.character(col_fraction)) col_fraction <- as.character(col_fraction)
    if (!is.character(col_subtype)) col_subtype <- as.character(col_subtype)
    
    ## Argument with single, non-NULL option (with or without quotation mark)
    #  (Never NULL since the default is not NULL)
    id <- rlang::as_string(rlang::enexpr(id))
    species <- rlang::as_string(rlang::enexpr(species))
  }
  
  if (FALSE) {
    ## Special handling of character vectors
    levels_fraction <- eval(rlang::enexpr(levels_fraction))
    levels_subtype <- eval(rlang::enexpr(levels_subtype))
  }
  
  check_subcell(
    levels_fraction = levels_fraction, levels_subtype = levels_subtype, 
    col_subtype = col_subtype)
  
  # label_scheme_sub <- label_scheme %>% dplyr::filter(!is.na(!!rlang::sym(col_select)))
  
  ### (1) Obtain aggregated expression data  
  # Need `!!` since `anal_prnTrend` is a UI function
  anal_prnTrend(
    col_select = !!col_select, 
    col_group = !!col_group, 
    col_order = !!col_order,
    col_fraction = !!col_fraction, 
    col_subtype = !!col_subtype, 
    n_clust = n_clust, 
    iter.max = 250, 
    filter_by_npep = exprs(prot_n_pep >= !!min_prot_n_peps), )
  
  ### (2) Optional visualzation of trend clusters
  if (make_plot) {
    try(
      plot_prnTrend(
        col_select = !!col_select,
        levels_group = !!levels_group, 
        col_order = !!col_order,
        ymax = 10, 
        ymin = -10,
        limitsize = FALSE,
        ncol = 4, 
      )
    )
  }
  
  ### (3) Classify subcellular
  plot_prnSubcellular(
    col_select = !!col_select, 
    levels_group = levels_group, 
    col_order = !!col_order, 
    col_fraction = !!col_fraction, levels_fraction = levels_fraction, 
    col_subtype = !!col_subtype, levels_subtype = levels_subtype, 
    ncol = 6, limitsize = FALSE, make_plot = make_plot, 
  )
}


#' Check the parameters for subcellular analysis.
#'
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param levels_fraction The levels of subcellular fractions, e.g., CP, NP.
#'   ChA.
#' @param col_subtype The column key in \code{expt_smry.xlsx} that defines
#'   sample subtype of control, treated etc.
check_subcell <- function (levels_fraction = NULL, levels_subtype = NULL, 
                           col_subtype = NULL)
{
  if (!(n_levs <- length(levels_fraction))) {
    stop("'levels_fraction' cannot be empty.")
  } else if (n_levs == 1L) {
    if (levels_fraction == "levels_fraction") {
      stop("Use 'levels_fraction = !!levels_fraction' ", 
           "instead of 'levels_fraction = levels_fraction'.")
    }
    
    stop("Need more than one level for 'levels_fraction'.")
  }
  
  if ((n_subtypes <- length(levels_subtype)) == 1L) {
    if (levels_subtype == "levels_subtype") {
      stop("Use 'levels_subtype = !!levels_subtype' ", 
           "instead of 'levels_subtype = levels_subtype'.")
    }
  }
  
  if (n_subtypes && (no_subtype <- !length(col_subtype))) {
    stop("Defined 'levels_subtype' but not 'col_subtype'.")
  }
}


#' Plot UMAP of sub-cellular scores.
#'
#' Inputs are expression data of normalized intensities.
#'
#' @param dat_dir The current working directory.
#' @param df_expr A matrix of expression data.
#' @param id An identifier.
#' @param n_clust Not used. The number of clusters.
#' @param levels_ref_compartments The factor levels of subcellular
#'   compartments. Set \code{levels_ref_compartments = NULL} for automatic
#'   determination.
#' @param df2 Not used. The names of secondary input files.
#' @param filepath A file path of inputs and outputs.
#' @param filename An output filename.
#' @param ext An extension of output filename.
#' @param wrap_facets Wrap facets or not.
#' @param seed A seed for umap.
#'
#' @examples
#' \dontrun{
#' plot_prnSubcell_UMAP(
#'   dat_dir = dat_dir, df_expr = NULL, id = "gene",
#'   df2 = NULL, filepath = NULL, filename = NULL, ext = "pdf",
#'   wrap_facets = TRUE, seed = 123, nrow = 2, ncol = 7, width = 16, height = 6)
#' }
plot_prnSubcell_UMAP <- function (
    dat_dir = NULL, df_expr = NULL, id = "gene", n_clust = 16:36, 
    levels_ref_compartments = c(
      "Cytosol", "ER", "Mitochondrion", "Nucleus", "Chromatin", "Lysosome", 
      "Plasma_membrane", "Proteasome", "Golgi", "Peroxisome", "Endosome", 
      "Ribosome", "Centrosome", "P_body", "Stress_granule"), 
    col_subtype = NULL, levels_subtype = NULL, 
    df2 = NULL, filepath = NULL, filename = NULL, ext = "png", use_cache = TRUE, 
    wrap_facets = TRUE, seed = 123, ...) 
{
  dots <- rlang::enexprs(...)
  
  if (is.null(dat_dir)) dat_dir  <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  
  ## Find input `df2` of SVM probabilities
  fns_svmprob <- find_trend_df2(
    df2 = df2, n_clust = NULL, 
    pat = "Protein_svmProb_Z_nclust\\d+[_]{0,1}.*\\.tsv", 
    ext = "tsv", scale_log2r = TRUE, impute_na = FALSE, filepath = filepath)
  
  if (!length(fns_svmprob)) {
    warning("No SVM results found.")
    return(NULL)
  }
  
  ## Load classification results
  df <- readr::read_tsv(file.path(filepath, fns_svmprob[[1]]))
  
  if (FALSE) {
    sub_types <- sub(paste0(".*_nclust\\d+[_]{0,1}(.*)\\.", "tsv"), "\\1", fns_svmprob)
    
    # Use for-loop: lapply crashed sometime
    df <- local({
      if (!length(mts <- match(levels_subtype, sub_types))) {
        mts <- seq_along(fns_svmprob)
      }
      
      filepaths <- file.path(filepath, fns_svmprob[mts])
      dfs <- vector("list", length(filepaths))
      for (i in seq_along(filepaths)) {
        dfs[[i]] <- readr::read_tsv(filepaths[[i]])
      }
      
      df <- dplyr::bind_rows(dfs)
    })
  }

  ## Compile subcellular factor levels
  if (!length(levels_ref_compartments)) {
    levels_ref_compartments <- local({
      univ <- unique(df$compartment)
      levs <- names(sort(table(df[["sub_location"]]), decreasing = TRUE))
      levs <- c(levs, univ[!univ %in% levs], "Unknown")
    })
  }
  
  print(paste("The levels of 'levels_ref_compartments': ", 
              paste(levels_ref_compartments, collapse = ", ")))

  df <- df |>
    dplyr::select(dplyr::one_of(id, "sub_location", col_subtype)) |>
    dplyr::filter(!is.na(sub_location)) |>
    unique()
  
  if (length(levels_ref_compartments)) {
    df <- df |>
      dplyr::mutate(
        sub_location = factor(sub_location, levels = levels_ref_compartments))
  }

  if (length(col_subtype) && length(levels_subtype) && 
      col_subtype %in% names(df)) {
    df <- df |>
      dplyr::mutate(
        !!col_subtype := 
          factor(!!rlang::sym(col_subtype), levels = levels_subtype), )
  }
  
  if (is.null(filename)) {
    fn_prefix <- paste0(tools::file_path_sans_ext(fns_svmprob[[1]]), "_umap")
    fn_suffix <- ext
    filename  <- paste0(fn_prefix, ".", fn_suffix)
  } else {
    fn_prefix <- tools::file_path_sans_ext(filename)
    fn_suffix <- tools::file_ext(filename)
  }
  
  ## Load the matrices of normalized intensities
  #  (see also 'prepSubcell' and 'make_subcell_expr')
  if (file.exists(fn_expr <- file.path(filepath, "df_exprs.rds"))) {
    df_exprs  <- readRDS(fn_expr)
    sub_types <- names(df_exprs)
  } else {
    warning("Expression file not found: ", file.path(filepath, "df_exprs.rds"))
    return(NULL)
  }
  rm(list = "fn_expr")
  
  file_umaps <- file.path(filepath, paste0("umap_", sub_types, ".tsv"))
  
  if (use_cache && all(file.exists(file_umaps))) {
    df_umap <- lapply(file_umaps, readr::read_tsv) |>
      dplyr::bind_rows()
  } else {
    set.seed(seed)
    
    df_umap <- mapply(function (df_expr, sub_type, file_umap) {
      res_umap <- umap::umap(df_expr)
      df_umap  <- res_umap$layout[, 1:2]
      colnames(df_umap) <- paste0("UMAP", 1:2)
      rownames(df_umap) <- rownames(df_expr)
      df_umap <- tibble::rownames_to_column(data.frame(df_umap), id)
      
      if (length(col_subtype)) {
        df_umap[[col_subtype]] <- sub_type
      }
      
      readr::write_tsv(df_umap, file_umap)
      
      df_umap
    }, df_exprs, sub_types, file_umaps, SIMPLIFY = FALSE) |>
      dplyr::bind_rows()
  }

  # Join predicted 'sub_location'
  df_plot <- tibble::tibble(df_umap)
  if (length(col_subtype)) {
    df_plot <- df_plot |>
      dplyr::left_join(df[, c(id, "sub_location", col_subtype)], 
                       by = c(id, col_subtype))
  } else {
    df_plot <- df_plot |>
      dplyr::left_join(df[, c(id, "sub_location")], by = id)
  }

  df_plot <- df_plot |>
    dplyr::mutate(
      sub_location = as.character(sub_location), 
      sub_location = ifelse(is.na(sub_location), "Unknown", sub_location), 
      sub_location = factor(sub_location, levels = levels_ref_compartments)) |>
    unique()
  
  my_colors <- local({
    set1_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                     "#FFFF33", "#A65628", "#F781BF")
    add_colors  <- c("#008080", "#3A1188", "#808000", "#882255", "#00BFFF", 
                     "#99d8c9", "#fee6ce")
    last_color  <- "#999999"
    
    n_comps <- length(levels_ref_compartments)
    
    if (any(is.na(n_comps)) || "Unknown" %in% levels_ref_compartments) {
      rng_colors <- 1:(n_comps - 1L)
    } else {
      rng_colors <- 1:n_comps
    }

    my_colors <- c(set1_colors, add_colors)[rng_colors]
    my_colors <- c(my_colors, last_color)
    
    if (any(is.na(my_colors))) {
      warning("More compartments than the number of colors.")
    }
    
    my_colors
  })

  ncol   <- dots$ncol
  nrow   <- dots$nrow
  size   <- dots$size
  alpha  <- dots$alpha
  width  <- dots$width
  height <- dots$height
  dpi    <- dots$dpi
  
  n_subtypes <- length(levels_subtype)
  
  if (is.null(nrow))  nrow <- 1
  if (is.null(width))  width <- n_subtypes * 8 / nrow + 2
  if (is.null(height)) height <- 8 * nrow
  if (is.null(ncol))  ncol <- 2
  
  if (is.null(size))  size <- .01
  if (is.null(alpha)) alpha <- .5
  if (is.null(dpi)) dpi <- 300
  
  if (is.null(ncol)) {
    if (is.null(nrow)) {
      ncol <- 4L
      nrow <- ceiling(n_subtypes / ncol)
    } else {
      ncol <- ceiling(n_subtypes / nrow)
    }
  } else {
    nrow <- ceiling(n_subtypes / ncol)
  }
  
  dots$ncol <- NULL
  dots$nrow <- NULL
  dots$width <- NULL
  dots$height <- NULL
  dots$alpha <- NULL
  dots$size <- NULL
  dots$dpi <- NULL
  
  p <- ggplot() +
    geom_point(
      data = df_plot, 
      aes(x = UMAP1, y = UMAP2, color = sub_location), 
      size = .02, 
      alpha = 0.5
    ) +
    scale_color_manual(values = my_colors, drop = FALSE) + 
    theme_minimal(base_size = 14) +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Location"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"), 
      legend.title = element_blank()
    ) + 
    guides(
      color = guide_legend(override.aes = list(size = 3, alpha = 1))
    )
  
  if (n_subtypes < 2L) {
    px <- p + 
      facet_wrap(
        ~ sub_location, ncol = ncol, labeller = label_value, drop = FALSE)
  } else {
    px <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_subtype), sub_location), 
        ncol = ncol, labeller = label_value)
    
    p <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_subtype)),
        ncol = ncol, labeller = label_value, )
  }
  
  ggsave(plot = px, file.path(filepath, paste0(fn_prefix, ".", fn_suffix)), 
         width = width, height = height, dpi = dpi)
  
  ggsave(plot = p, file.path(filepath, paste0(fn_prefix, "_combined.", fn_suffix)), 
         width = height + 2, height = height, dpi = dpi)
  
  invisible(df_plot)
}


#' Update subcellular protein or peptide tables by contrast.
#' 
#' A contrast can be nested.
#' 
#' @param contr_elements Contrast elements.
#' @param contr_nm Contrast names.
#' @param df A data frame.
update_subell_by_contr <- function (
    contr_elements = NULL, contr_nm = NULL, df = NULL)
{
  nas <- local({
    dfx <- data.frame(df[, contr_elements], check.names = FALSE)
    
    nas <- sapply(1:nrow(dfx), function(i) {
      vals <- dfx[i, ]
      # any(is.na(vals) | vals == 0)
      all(is.na(vals) | vals == 0)
    })
  })
  
  df <- local({
    cnms <- colnames(df)
    # No special characters
    cols <- grepl(paste0("\\(", contr_nm, "\\)$"), cnms)
    
    # Contains special characters such as parentheses
    if (!any(cols)) {
      cols <- grepl(paste0("(", contr_nm, ")"), cnms, fixed = TRUE)
    }
    
    if (!any(cols)){
      warning("Data column not found.")
      return(df)
    }
    
    dfw <- df[, cols, drop = FALSE]
    dfw[nas, ] <- NA_real_
    df[, cols] <- dfw
    
    df
  })
  
  df
}


#' Update subcellular protein or peptide tables by formula.
#' 
#' There can be multiple contrasts under a formula.
#' 
#' @param fml_nm Formula name.
#' @param fml Formula.
#' @param df A data frame.
#' @param label_scheme_sub Metadata.
update_subcell_by_fml <- function (
    fml_nm = NULL, fml = NULL, df = NULL, label_scheme_sub = NULL)
{
  fml_ops     <- prepFml(fml, label_scheme_sub)
  elements    <- fml_ops[["elements"]]
  contr_mat   <- fml_ops[["contr_mat"]]
  sub_contrs  <- fml_ops[["sub_contrs"]] # NULL if no nested structure
  # label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub
  
  if (FALSE && !is.null(sub_contrs)) {
    return(df)
  }
  
  contr_nms   <- colnames(contr_mat)
  
  contrs_elements <- if (is.null(sub_contrs)) {
    strsplit(contr_nms, "-")
  } else {
    lapply(sub_contrs, function (x) unique(unlist(strsplit(x, "-"))))
  }
  
  for (i in seq_along(contr_nms)) {
    contr_nm <- contr_nms[[i]]
    contr_elements <- contrs_elements[[i]]
    
    df <- update_subell_by_contr(
      contr_elements = contr_elements, contr_nm = contr_nm, df = df)
  }
  
  df
}


#' Update subcellular protein or peptide tables.
#'
#' Nullify pVals for proteins found but not belonging to a subcellular fraction.
#' Intensity fields of \code{I[0-9]{3}} and log2Ratio files of
#' \code{[NZ]{0,1}_log2R} remains intact.
#'
#' @param df Automatically determined: a data frame from
#'   \code{Protein_pVals.txt} or \code{Peptide_pVals.txt}.
#' @param dat_dir The working directory.
#' @param data_type The type of data in either \code{Protein} or \code{Peptide}.
#' @param filepath Automatically determined: the fiepath to \code{Protein} or
#'   \code{Peptide} table.
#' @param fml_nms Automatically determined: the formula names used in
#'   \code{prnSig} or \code{pepSig}.
#' @param label_scheme_sub Automatically determined: metadata.
#' @export
update_subcell <- function (
    df = NULL, data_type = c("Protein", "Peptide"), 
    dat_dir = NULL, filepath = NULL, fml_nms = NULL, label_scheme_sub = NULL, 
    ...)
{
  if (is.null(dat_dir)) {
    dat_dir <- get_gl_dat_dir()
  }
  
  data_type <- match.arg(data_type)
  filename  <- paste0(data_type, "_pVals.txt")
  
  group_psm_by <- tryCatch(match_call_arg(normPSM, group_psm_by), error = function(e) NA)
  
  if (is.na(group_psm_by)) {
    group_psm_by <- tryCatch(
      match_call_arg(makePepDIANN, group_psm_by), 
      error = function(e) NA)
  }
  
  if (is.na(group_psm_by)) {
    group_psm_by <- "pep_seq_mod"
  }
  
  group_pep_by <- tryCatch(match_call_arg(normPSM, group_pep_by), error = function(e) NA)
  
  if (is.na(group_pep_by)) {
    group_pep_by <- tryCatch(
      match_call_arg(makePepDIANN, group_psm_by), 
      error = function(e) NA)
  }
  
  if (is.na(group_pep_by)) {
    group_pep_by <- "gene"
  }
  
  if (is.null(filepath)) {
    filepath <- file.path(dat_dir, data_type, "Model")
  }
  
  if (is.null(df)) {
    df <- readr::read_tsv(file.path(filepath, filename))
  }
  
  if (is.null(label_scheme_sub)) {
    label_scheme_sub <- load_ls_group(dat_dir = dat_dir)
  }
  
  dots <- rlang::enexprs(...)
  
  dots <- concat_fml_dots(
    fmls = dots[grepl("^\\s*~", dots)], 
    fml_nms = fml_nms, 
    dots = dots[!grepl("^\\s*~", dots)], )
  
  lang_dots    <- dots[unlist(lapply(dots, is.language))]
  filter_dots  <- lang_dots[grepl("^filter_", names(lang_dots))]
  arrange_dots <- lang_dots[grepl("^arrange_", names(lang_dots))]
  select_dots  <- lang_dots[grepl("^select_", names(lang_dots))]
  dots         <- dots[!dots %in% c(filter_dots, arrange_dots, select_dots)]
  
  df <- df |>
    filters_in_call(!!!filter_dots) |>
    arrangers_in_call(!!!arrange_dots)
  
  rm(list = c("filter_dots", "arrange_dots", "select_dots"))
  
  fmls <- dots[grepl("^\\s*~", dots)]
  dots <- dots[!names(dots) %in% names(fmls)]
  
  fml_nms <- names(df) %>% 
    .[grepl("pVal", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique() %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls[names(fmls) %in% fml_nms]
  fml_nms <- fml_nms[map_dbl(fml_nms, function (x) which(x == names(fmls)))]
  
  if (!length(fml_nms)) {
    stop("No formula (names) matched to those in 'pepSig' or 'prnSig'.")
  }
  
  ## Add subcellular specificity
  nms_df <- colnames(df)
  
  res_subcell <- add_subcell(
    df= df, dat_dir = dat_dir, group_pep_by = group_pep_by, 
    is_subcellular = TRUE, scale_log2r = TRUE, impute_na = FALSE)
  df <- res_subcell[["df"]]
  col_fraction <- res_subcell[["col_fraction"]]
  col_subtype <- res_subcell[["col_subtype"]]
  col_subgroup <- res_subcell[["col_group"]]
  rm(list = "res_subcell")
  
  for (i in seq_along(fml_nms)) {
    fml_nm <- fml_nms[[i]]
    fml <- fmls[[i]]
    
    df <- update_subcell_by_fml(
      fml_nm = fml_nm, fml = fml, df = df, label_scheme_sub = label_scheme_sub)
  }
  
  df <- df[, names(df) %in% nms_df]
  
  path_out <- file.path(dat_dir, data_type, "Model", filename)
  readr::write_tsv(df, path_out)
  
  invisible(df)
}


#' Finds low-fidelity entries in translocation measures.
#' 
#' The low-fidelity entries are excluded from being hits of translocation.
#' 
#' @param dat_dir The current working directory.
#' @param label_scheme_sub Metadata.
#' @param col_group The column key in metadata linked to the information of
#'   sample groups.
#' @param from The \code{from} compartment.
#' @param to The \code{to} compartment.
#' @param col_fraction The column key in metadata to subcellular fractions.
#' @param col_subtype The column key in metadata to smple subtypes.
#' @param levels_fraction The levels of subcellular compartments.
#' @param levels_subtype The levels of sample subtypes or treatment conditions.
find_bad_transloc <- function (
    from = "CP", to = "NP", col_group = "Group", 
    dat_dir = NULL, label_scheme_sub = NULL, 
    col_fraction = "Fraction", levels_fraction = c("CP", "NP", "ChA"), 
    col_subtype = "Type", levels_subtype = c("Veh", "TG")) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  
  dfN  <- make_dfn(dat_dir = dat_dir,, col_group = col_group)
  grps <- colnames(dfN)
  
  ty_1 <- levels_subtype[[1]]
  ty_2 <- levels_subtype[[2]]
  
  sc_lookup <- label_scheme_sub |>
    dplyr::select(dplyr::one_of(col_group, col_fraction, col_subtype)) |>
    unique() |>
    dplyr::rename(Sample_ID := col_group)
  labels_frac <- sc_lookup[[col_fraction]]
  labels_subtype <- sc_lookup[[col_subtype]]

  dfN_from <- dfN[, labels_frac %in% from]
  dfN_to <- dfN[, labels_frac %in% to]
  bads_from <- rowSums(dfN_from > 1) == 0
  bads_to <- rowSums(dfN_to > 1) == 0
  
  dfN_sub_1 <- dfN[, labels_subtype %in% ty_1 & labels_frac %in% c(from, to)]
  dfN_sub_2 <- dfN[, labels_subtype %in% ty_2 & labels_frac %in% c(from, to)]
  
  bads_sub_1 <- rowSums(dfN_sub_1 > 1) == 0
  bads_sub_2 <- rowSums(dfN_sub_2 > 1) == 0
  
  bads <- bads_from | bads_to | bads_sub_1 | bads_sub_2
  
  ids_excl <- rownames(dfN)[bads]
}


#' Find translocation outliers between a pair of subcellular compartments.
#'
#' @param dat_dir The current working directory.
#' @param filepath An output filepath.
#' @param key_col The key column for measures of translocation. The value is in
#'   one of \code{c("purity", "loc_score", "entropy")}.
#' @param id An identifier.
#' @param use_updated_subcell Logical; if TRUE, use clean subcellular metrics of
#'   purity etc. Otherwise, use raws subcellular metrics.
#' @param col_group The column key in metadata linked to the information of
#'   sample groups.
#' @param scale_log2r Not used. Logical; scale log2Ratio values or not.
#' @param impute_na Not used. Loigcal; impute missing values or not.
#' @param from_to The pair of subcellular compartments that defines the scope of
#'   translocation.
#' @param levels_ref_compartments The levels of subcellular compartments.
#' @param colors A color scheme.
#' @param color_outlier The color of outliers.
#' @param levels_subtype The levels of sample subtypes or treatment conditions.
#' @param ... Additional arguments.
#' 
#' @examples \donttest{
#' res <- find_from_to_transloc(
#'   from_to = list(c("CP", "NP"), c("NP", "ChA")), 
#'   levels_subtype = c("Veh", "TG"), transloc_co = .3, width = 9, height = 3.4)
#' }
#' @return A table of \code{corr_subcell_combined.tsv}.
#' @export
find_from_to_transloc <- function (
    from_to = list(c("CP", "NP"), c("NP", "ChA")), 
    key_col = c("purity", "loc_score", "entropy"), 
    levels_ref_compartments = c("CP", "NP", "ChA"), 
    colors = c("CP" = "coral", "NP" = "skyblue", "ChA" = "mediumseagreen"), 
    color_outlier = "red", 
    levels_subtype = NULL, transloc_co = .3, 
    use_updated_subcell = TRUE, 
    col_select = "Select", col_group = "Group", col_order = "Order", 
    id = "gene", dat_dir = NULL, filepath = NULL, 
    scale_log2r = TRUE, impute_na = FALSE, 
    ...)
{
  lens <- lengths(from_to)
  if (!all(lens == 2L)) {
    stop("Lengths of 'from_to' needs to be two for each pair.")
  }
  rm(list = "lens")
  
  key_col <- match.arg(key_col)
  
  if (!length(col_group)) {
    warning("'col_group' not specified; Use the default 'Group'.")
    col_order <- "Group"
  }
  
  if (!length(col_order)) {
    warning("'col_order' not specified; Use the default 'Order'.")
    col_order <- "Order"
  }
  
  if (!length(col_select)) {
    warning("'col_select' not specified; Use the default 'Select'.")
    col_select <- "Select"
  }
  
  dots <- rlang::enexprs(...)
  lang_dots     <- dots[unlist(lapply(dots, is.language))]
  filter2_dots  <- lang_dots[grepl("^filter2_", names(lang_dots))]
  arrange2_dots <- lang_dots[grepl("^arrange2_", names(lang_dots))]
  dots <- dots[!dots %in% c(filter2_dots, arrange2_dots)]
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Transloc")
  if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
  
  label_scheme <- load_ls_group(dat_dir = dat_dir)
  
  label_scheme_sub <- label_scheme %>% 
    dplyr::filter(!is.na(.data[[col_select]]))
  
  filepath_trend <- file.path(dat_dir, "Protein", "Trend")
  
  files_trend <- find_trend_df2(
    scale_log2r = scale_log2r, impute_na = impute_na, filepath = filepath_trend)
  
  if (!length(files_trend)) {
    stop("No trend results found.", " Run 'anal_prnTrend' first.")
  }
  
  ## Update loc_score, purity, entropy and trend results
  file_one <- files_trend[[1]]
  if (use_updated_subcell) {
    df <- readr::read_tsv(file.path(filepath_trend, "subcell_update", file_one))
  } else {
    df <- readr::read_tsv(file.path(filepath_trend, file_one))
  }
  
  ## Plots
  #  Load updated loc_score, purity, entropy and trend results
  col_fraction <- df[["col_fraction"]][[1]]
  col_subtype  <- df[["col_subtype"]][[1]]
  
  if (!length(levels_subtype)) {
    levels_subtype <- unique(df[[col_subtype]])
  }
  
  n_subtype <- length(levels_subtype)
  if (n_subtype != 2L) {
    stop("Proive two values to 'levels_subtype' for translocation analysis.")
  }
  
  cdn_x <- levels_subtype[[1]]
  cdn_y <- levels_subtype[[2]]
  
  df <- df |>
    dplyr::filter(!is.na(sub_location))

  df <- df |>
    dplyr::select(dplyr::one_of(c(id, key_col, col_fraction, col_subtype))) |>
    dplyr::filter(!is.na(.data[[id]])) |> 
    tidyr::pivot_wider(
      id_cols = c(id, col_fraction), 
      names_from = col_subtype, 
      values_from = key_col)
  
  df <- local({
    for (i in seq_along(levels_subtype)) {
      subtype <- levels_subtype[[i]]
      nas <- is.na(df[[subtype]])
      df[[subtype]][nas] <- 0.0
    }
    
    df
  })

  df <- df |>
    dplyr::mutate(
      delta_loc_score  = .data[[cdn_y]] - .data[[cdn_x]], 
      rdelta_loc_score = (.data[[cdn_y]] - .data[[cdn_x]]) / .data[[cdn_x]]) |>
    dplyr::arrange(delta_loc_score)
  
  df <- df |> 
    dplyr::mutate(
      !!col_fraction := factor(!!rlang::sym(col_fraction), 
                               levels = levels_ref_compartments)) |>
    tidyr::drop_na(dplyr::all_of(levels_subtype))
  
  # Plot arguments
  ncol   <- dots$ncol
  nrow   <- dots$nrow
  size   <- dots$size
  alpha  <- dots$alpha
  width  <- dots$width
  height <- dots$height
  dpi    <- dots$dpi
  
  n_subtypes <- length(levels_subtype)
  
  if (is.null(width))  width <- n_subtypes * 8 / nrow + 2
  if (is.null(height)) height <- 8 * nrow
  if (is.null(ncol))  ncol <- 2
  if (is.null(nrow))  nrow <- 1
  if (is.null(size))  size <- .01
  if (is.null(alpha)) alpha <- .5
  if (is.null(dpi)) dpi <- 300
  
  if (is.null(ncol)) {
    if (is.null(nrow)) {
      ncol <- 4L
      nrow <- ceiling(n_subtypes / ncol)
    } else {
      ncol <- ceiling(n_subtypes / nrow)
    }
  } else {
    nrow <- ceiling(n_subtypes / ncol)
  }
  
  dots$ncol <- NULL
  dots$nrow <- NULL
  dots$width <- NULL
  dots$height <- NULL
  dots$alpha <- NULL
  dots$size <- NULL
  dots$dpi <- NULL
  
  # Visualization
  dfs <- ps <- vector("list", length(from_to))
  
  for (i in seq_along(from_to)) {
    pair <- from_to[[i]]
    from <- pair[[1]]
    to <- pair[[2]]
    
    bad_ids <- find_bad_transloc(
      from = from, to = to, col_group = col_group, 
      dat_dir = dat_dir, label_scheme_sub = label_scheme_sub, 
      col_fraction = col_fraction, levels_fraction = levels_ref_compartments, 
      col_subtype = col_subtype, levels_subtype = levels_subtype) 
    
    dfs[[i]] <- dfi <- df |>
      dplyr::filter(.data[[col_fraction]] %in% pair, 
                    !.data[[id]] %in% bad_ids) |>
      dplyr::mutate(is_outlier = 
                      ifelse(abs(delta_loc_score) > transloc_co, TRUE, FALSE))
    
    ps[[i]] <- make_transloc_scatter(
      df = dfi, cdn_x = cdn_x, cdn_y = cdn_y, col_fraction = col_fraction, 
      color_outlier = color_outlier, colors = colors)

    fn_prefix <- paste0("corr_subcell_", from, "_", to)
    out_name <- file.path(filepath, paste0(fn_prefix, ".tsv"))
    readr::write_tsv(dfi, out_name)
    
    ggsave(file.path(filepath, paste0(fn_prefix, ".png")), 
           width = width, height = height)
  }
  
  df_combined <- dfs |>
    dplyr::bind_rows() |>
    unique()
  
  out_name_combined <- file.path(filepath, "corr_subcell_combined.tsv")
  readr::write_tsv(df_combined, out_name_combined)
  
  p <- make_transloc_scatter(
    df = df_combined, cdn_x = cdn_x, cdn_y = cdn_y, col_fraction = col_fraction, 
    color_outlier = color_outlier, colors = colors)
  
  ggsave(file.path(filepath, "corr_subcell_combined.png"), 
         width = width, height = height)
  
  df_combined
}


#' Update subcellular purity, entropy and localization scores.
#' 
#' Also update trend results under sub directory \code{subcell_update}.
#' 
#' @param id An identifier.
#' @param dat_dir The current working directory.
#' @param filepath The filepath to protein \code{Trend} results.
#' @param files The filenames to protein \code{Trend} results.
#' @param col_select The column key in metadata linked to the information of
#'   sample being selected.
#' @param col_group The column key in metadata linked to the information of
#'   sample groups.
#' @param col_oder The column key in metadata linked to the information of
#'   sample orders.
#' @param col_fraction The column key in metadata linked to the information of
#'   sample subcellular fractions.
#' @param col_subtype The column key in metadata linked to the information of
#'   sample subtypes or treatment conditions.
#' @param levels_subtype The levels of sample subtypes.
#' @param scale_log2r Not used. Logical; scale log2Ratio values or not.
#' @param impute_na Not used. Loigcal; impute missing values or not.
#' @param ... Additional arguments.
#' 
#' @examples \donttest{
#' update_subcell_metrics(col_group = "Group", levels_subtype = c("Veh", "TG"))
#' }
#' 
#' @export
update_subcell_metrics <- function (
    id = "gene", dat_dir = NULL, filepath = NULL, files = NULL, 
    col_select = NULL, col_group = NULL, col_order = NULL, 
    # col_fraction = NULL, col_subtype = NULL, 
    levels_subtype = NULL, 
    scale_log2r = TRUE, impute_na = FALSE, ...)
{
  dots <- rlang::enexprs(...)
  lang_dots     <- dots[unlist(lapply(dots, is.language))]
  filter2_dots  <- lang_dots[grepl("^filter2_", names(lang_dots))]
  arrange2_dots <- lang_dots[grepl("^arrange2_", names(lang_dots))]
  dots <- dots[!dots %in% c(filter2_dots, arrange2_dots)]
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  
  if (!length(files)) {
    files <- find_trend_df2(
      scale_log2r = scale_log2r, impute_na = impute_na, filepath = filepath)
  }
  
  df <- readr::read_tsv(file.path(filepath, files[[1]]))
  
  # Retain only proteins with subcellular assignments
  df <- df |>
    dplyr::filter(!is.na(sub_location))
  
  col_fraction <- df[["col_fraction"]][[1]]
  col_subtype  <- df[["col_subtype"]][[1]]
  
  if (!length(levels_subtype)) {
    levels_subtype <- unique(df[[col_subtype]])
  }
  
  if (!length(col_group)) {
    warning("'col_group' not specified; Use the default 'Group'.")
    col_order <- "Group"
  }
  
  if (!length(col_order)) {
    warning("'col_order' not specified; Use the default 'Order'.")
    col_order <- "Order"
  }
  
  if (!length(col_select)) {
    warning("'col_select' not specified; Use the default 'Select'.")
    col_select <- "Select"
  }
  
  label_scheme <- load_ls_group(dat_dir = dat_dir)
  
  label_scheme_sub <- label_scheme %>%
    dplyr::filter(!is.na(.data[[col_select]]))
  
  # Aggregated intensity matrix
  df_mean_int <- df[, c(id, "group", "log10Int")] |>
    dplyr::mutate(log10Int = 10^log10Int) |>
    tidyr::pivot_wider(
      id_cols = id, 
      names_from = c("group"), 
      values_from = "log10Int",
      values_fill = 1, 
    ) |>
    as.data.frame()
  
  rownames(df_mean_int) <- df_mean_int[[id]]
  df_mean_int[[id]] <- NULL
  
  ## Values of purity, entropy and localization score
  grps <- colnames(df_mean_int)
  
  if (!length(n_subtypes <- length(levels_subtype))) {
    levels_subtype <- "unknown"
    n_subtypes <- 1L
    col_subtype <- "sample_type_unknown"
    label_scheme_sub[[col_subtype]] <- levels_subtype
  }
  
  # Subcellular lookup
  sc_lookup <- label_scheme_sub |>
    dplyr::select(dplyr::one_of(col_group, col_fraction, col_subtype)) |>
    unique() |>
    dplyr::rename(Sample_ID := col_group)
  
  purities <- vector("list", n_subtypes)
  for (i in seq_along(levels_subtype)) {
    sids_sub <- sc_lookup |>
      dplyr::filter(.data[[col_subtype]] == levels_subtype[[i]]) |>
      dplyr::pull("Sample_ID")
    
    dfy <- df_mean_int[, grps %in% sids_sub]
    purities[[i]] <- dfy / rowSums(dfy, na.rm = TRUE)
  }
  rm(list = c("dfy", "sids_sub"))
  
  entropies <- lapply(purities, function (x) {
    x[x == 0] <- 1E-4
    -rowSums(x * log10(x))
  })
  
  scores_loc <- mapply(function (p, y) {
    p * (1 - y)
  }, purities, entropies, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  scores_loc <- lapply(scores_loc, function (x) {
    x[x < 1E-3] <- 0.0
    x
  })
  
  purities <- lapply(purities, function (x) {
    x[x < 1E-3] <- 0.0
    x
  })
  
  entropies <- lapply(entropies, function (x) {
    x[x < 1E-3] <- 0.0
    x
  })
  
  purities <- do.call(cbind, purities)
  scores_loc <- do.call(cbind, scores_loc)
  entropies <- do.call(cbind, entropies)
  colnames(entropies) <- levels_subtype
  
  nans <- sapply(purities, is.nan)
  purities[nans] <- NA_real_
  purities <- round(purities, digits = 3L)
  purities <- data.frame(purities, check.names = FALSE)
  purities <- purities[, grps]
  
  nans <- sapply(scores_loc, is.nan)
  scores_loc[nans] <- NA_real_
  scores_loc <- round(scores_loc, digits = 3L)
  scores_loc <- data.frame(scores_loc, check.names = FALSE)
  scores_loc <- scores_loc[, grps]
  
  nans <- sapply(entropies, is.nan)
  entropies[nans] <- NA_real_
  entropies <- round(entropies, digits = 3L)
  entropies <- data.frame(entropies, check.names = FALSE)
  
  purities <- purities |>
    tibble::rownames_to_column(id) |>
    tidyr::pivot_longer(
      cols = -id,
      names_to = "group",
      values_to = "purity"
    )
  
  entropies <- entropies |>
    tibble::rownames_to_column(id) |>
    tidyr::pivot_longer(
      cols = -id,
      names_to = col_subtype,
      values_to = "entropy"
    )
  
  scores_loc <- scores_loc |>
    tibble::rownames_to_column(id) |>
    tidyr::pivot_longer(
      cols = -id,
      names_to = "group",
      values_to = "loc_score"
    )
  
  ## Update purity, loc_score etc. in trend results
  filepath_updated <- file.path(filepath, "subcell_update")
  dir.create(filepath_updated, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_along(files)) {
    fni <- files[[i]]
    dfi <- readr::read_tsv(file.path(filepath, fni))
    
    dfx <- list(
      dplyr::select(dfi, -dplyr::one_of("purity", "entropy", "loc_score")), 
      purities, 
      scores_loc) |>
      purrr::reduce(dplyr::left_join, by = c(id, "group")) |>
      dplyr::left_join(entropies, by = c(id, col_subtype))
    
    dfi <- dfx[, colnames(dfi)]
    readr::write_tsv(dfi, file.path(filepath_updated, fni))
  }
  
  invisible(
    list(
      purity = purities, 
      entropy = entropies, 
      score_loc = scores_loc, 
      sc_lookup = sc_lookup, 
      col_fraction = col_fraction, 
      col_subtype = col_subtype
    )
  )
}


#' Helper of making translocation scatter plots.
#' 
#' @param df A data frame.
#' @param cdn_x The condition of X-axis.
#' @param cdn_y The condition of Y-axis.
#' @param col_fraction The column key in metadata linked to the information of
#'   sample subcellular fractions.
#' @param color_outlier The color for uses with outliers.
#' @param colors The color scheme.
make_transloc_scatter <- function (
    df, cdn_x, cdn_y, col_fraction, color_outlier, colors)
{
  p <- ggplot(df, aes(x = .data[[cdn_x]], y = .data[[cdn_y]], 
                 color = !!rlang::sym(col_fraction))) +
    geom_point(size = 0.8, alpha = 0.4) +
    geom_point(
      data = subset(df, is_outlier),
      color = color_outlier, size = 1.5, alpha = 0.5
    ) +
    geom_text_repel(
      data = subset(df, is_outlier),
      aes(label = gene),
      color = "black",
      size = 1.8,
      max.overlaps = 15,
      segment.size = .1, 
      segment.color = "grey30",
      show.legend = FALSE
    ) +
    # scale_color_brewer(palette = "Set2") +
    scale_color_manual(values = colors) + 
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1), 
                       breaks = scales::breaks_width(0.4)) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1), 
                       breaks = scales::breaks_width(0.4)) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.grid = element_blank(), 
      # legend.text = element_blank(), 
      # legend.title = element_blank(), 
      legend.position = "none"
    )
  
  p <- p + 
    facet_wrap(
      vars(!!rlang::sym(col_fraction)), ncol = 3, 
      labeller = label_value, drop = FALSE)
}


#' Update protein tables of \code{Protein.txt} etc. by subcellular locations.
#' 
#' Not used.
#' 
#' @param id Protein identifier.
#' @param dat_dir The current working directory.
#' @param filepath The filepath to trend data.
#' @param scale_log2r Not used. Logical; scale log2Ratio values or not.
#' @param impute_na Not used. Loigcal; impute missing values or not.
#' @param ... Additional parameters.
update_subcellular <- function (id = "gene", dat_dir = NULL, filepath = NULL, 
                                scale_log2r = TRUE, impute_na = FALSE, ...) 
{
  if (is.null(dat_dir)) {
    dat_dir <- get_gl_dat_dir()
  }
  
  if (is.null(filepath)) {
    filepath <- file.path(dat_dir, "Protein", "trend")
  }
  
  dots <- rlang::enexprs(...)
  label_scheme_sub <- load_ls_group(dat_dir = dat_dir)
  
  ## Find input `df2`
  df2 <- find_trend_df2(
    df2 = NULL, n_clust = NULL, scale_log2r = scale_log2r, 
    impute_na = impute_na, filepath = filepath)
  
  if (!length(df2)) {
    stop("No trend results found.", " Run 'anal_prnTrend' first.")
  }
  
  df_trend <- readr::read_tsv(file.path(filepath, df2[[1]]))
  col_fraction <- df_trend[["col_fraction"]][[1]]
  col_subtype <- df_trend[["col_subtype"]][[1]]
  df_trend <- df_trend[, c(id, "loc_score", "sub_location", col_subtype)]
  
  hupdate_subcellular(
    df_trend = df_trend, dat_dir = dat_dir, id = id, 
    label_scheme_sub = label_scheme_sub, col_fraction = col_fraction, 
    col_subtype = col_subtype) 
}


#' Update subcellular protein and peptide tables.
#' 
#' Not used.
#' 
#' @param df_trend Trend results
#' @param dat_dir The working directory.
#' @param id A protein identifier.
#' @param label_scheme_sub Metadata.
#' @param col_fraction The column key in metadata linking to subcellular
#'   values.
#' @param col_subtye The column key metadata linking to sample sub-types.
hupdate_subcellular <- function (
    df_trend, dat_dir, id = "gene", label_scheme_sub = NULL, 
    col_fraction = NULL, col_subtype = NULL) 
{
  df_prn <- readr::read_tsv(file.path(dat_dir, "Protein", "Protein.txt"))
  df_pep <- readr::read_tsv(file.path(dat_dir, "Peptide", "Peptide.txt"))
  dfp_prn <- readr::read_tsv(file.path(dat_dir, "Protein", "Model", "Protein_pVals.txt"))
  dfp_pep <- readr::read_tsv(file.path(dat_dir, "Peptide", "Model", "Peptide_pVals.txt"))
  
  metas <- split(
    label_scheme_sub, label_scheme_sub[, c(col_fraction, col_subtype)])
  
  for (i in seq_along(metas)) {
    meta <- metas[[i]]
    
    gns <- df_trend |>
      dplyr::filter(
        sub_location == meta[[col_fraction]][[1]], 
        !!rlang::sym(col_subtype) == meta[[col_subtype]][[1]], ) |>
      dplyr::pull(id)
    
    pats <- paste0("\\(", meta[["Sample_ID"]], "\\)$")
    
    cols_prn <- grepl(paste0(pats, collapse = "|"), colnames(df_prn))
    cols_pep <- grepl(paste0(pats, collapse = "|"), colnames(df_pep))
    rows_prn <- !df_prn[[id]] %in% gns
    rows_pep <- !df_pep[[id]] %in% gns
    
    df_prn[rows_prn, cols_prn] <- NA_real_
    df_pep[rows_pep, cols_pep] <- NA_real_
    
    cols_prnp <- grepl(paste0(pats, collapse = "|"), colnames(dfp_prn))
    cols_pepp <- grepl(paste0(pats, collapse = "|"), colnames(dfp_pep))
    rows_prnp <- !dfp_prn[[id]] %in% gns
    rows_pepp <- !dfp_pep[[id]] %in% gns
    
    dfp_prn[rows_prnp, cols_prnp] <- NA_real_
    dfp_pep[rows_pepp, cols_pepp] <- NA_real_
  }
  
  readr::write_tsv(df_prn, file.path(dat_dir, "Protein", "Protein_subcell.txt"))
  readr::write_tsv(df_pep, file.path(dat_dir, "Peptide", "Peptide_subcell.txt"))
  readr::write_tsv(dfp_prn, file.path(dat_dir, "Protein", "Model", "Protein_pVals_subcell.txt"))
  readr::write_tsv(dfp_pep, file.path(dat_dir, "Peptide", "Model", "Peptide_pVals_subcell.txt"))
  
  invisible(df_prn)
}


