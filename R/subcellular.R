#' Subcellular plots.
#' 
#' @param col_cluster The column key of cluster nodes.
#' @param pat A pattern of input files.
#' @param ext The extension of input files.
#' @param data_type The type of data.
#' @param levels_subcellular The levels of subcellular fractions.
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
    id = "gene", dat_dir = NULL, col_group = NULL, col_order = NULL, 
    col_cluster = "cluster", 
    levels_subcellular = NULL, levels_subtype = NULL, tie_method = "none", 
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
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
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
      col_group = col_group, 
      col_order = col_order, 
      col_cluster = col_cluster, 
      levels_subcellular = levels_subcellular,
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
  col_fraction <- tempsub[[1]][["col_fraction"]]
  col_subtype <- tempsub[[1]][["col_subtype"]]
  col_x <- tempsub[[1]][["col_x"]]
  col_y <- tempsub[[1]][["col_y"]]
  rm(list = "tempsub")
  
  ## Classify trends
  tempdata <- lapply(
    df_meds, 
    hclassify_trends, 
    col_subtype = col_subtype, 
    col_fraction = col_fraction, 
    col_group = rlang::enexpr(col_group), 
    col_order = rlang::enexpr(col_order), 
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
  rm(list = "tempdata")
  
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
    
    # Not yet for svmProb: need column log10Int
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
#' @param filename_ref Optional; a file name with prepending path to a reference
#'   subcellular annotation file, e.g., UniProt and GO intersection.
#' @param key_subcellular The key of subcellular compartments from the reference
#'   data.
#' @param key_col The column key where its values are used for UMAP. 
#' @param seed A seed for umap.
#' @inheritParams plot_prnTrend
#' @inheritParams plot_prnSubcellular
#' @export
plot_prnSubcellular_UMAP <- function (
    df2 = NULL, levels_subcellular = NULL, levels_subtype = NULL, 
    key_subcellular = "location_slim3", key_col = "loc_score", 
    filename_ref = NULL, scale_log2r = TRUE, impute_na = FALSE, 
    filename = NULL, seed = 42, ...) 
{
  dots <- rlang::enexprs(...)
  dat_dir  <- get_gl_dat_dir()
  filepath <- file.path(dat_dir, "Protein", "Trend")
  
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
  
  # Obtain the column keys encoding fraction and sub-type information
  col_fraction <- df[["col_fraction"]][[1]]
  col_subtype  <- df[["col_subtype"]][[1]]

  # String `col_fraction` is the key that encodes sample fraction; 
  # Column `sub_location` contains classification results
  df_expr <- df[, c("gene", key_col, "sub_location", "group")] |>
    dplyr::filter(!is.na(sub_location)) |>
    dplyr::arrange(gene) |>
    unique() |>
    tidyr::pivot_wider(
      id_cols = "gene", 
      names_from = c("group"), 
      values_from = key_col,
      values_fill = 0, 
    )
  
  set.seed(seed)
  res_umap <- umap::umap(df_expr[, !colnames(df_expr) == "gene", drop = FALSE])
  
  umap_df <- res_umap$layout[, 1:2]
  colnames(umap_df) <- paste0("UMAP", 1:2)
  rownames(umap_df) <- df_expr$gene
  umap_df <- tibble::rownames_to_column(data.frame(umap_df), "gene")
  readr::write_tsv(umap_df, file.path(filepath, "loc_score_umap.tsv"))

  ## Plot of predicted subcellular results
  ok_subtypes <- sort(unique(df[[col_subtype]]))
  if (is.null(levels_subtype)) {
    levels_subtype <- ok_subtypes
  } else {
    if (length(bads <- levels_subtype[!levels_subtype %in% ok_subtypes])) {
      stop("levels_subtype not found: ", paste0(bads, collapse = ", "))
    }
  }
  
  n_subtypes <- length(levels_subtype)
  
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
        factor(!!rlang::sym(col_fraction), levels = levels_subcellular),
      sub_location = 
        factor(sub_location, levels = levels_subcellular),
    )
  
  df_plot <- df |>
    dplyr::left_join(umap_df, by = "gene") |> 
    dplyr::filter(!is.na(sub_location))

  p <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = sub_location)) +
    geom_point(shape = 46, alpha = alpha, size = size) +
    # Override the legend point size, shape, and transparency
    guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
    theme_bw() +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Subcellular\nLocation"
    ) +
    facet_wrap(
      vars(!!rlang::sym(col_subtype[1])), nrow = nrow, labeller = label_value) + 
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave(file.path(filepath, filename), width = width, height = height, 
         dpi = dpi) # , type = "cairo"
  
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
        # FIX: Adds 15 points of padding to the right side of the plot
        plot.margin = margin(t = 5, r = 15, b = 5, l = 5, unit = "pt")
      )
    
    ggsave(file.path(filepath, out_nmx), width = width, height = height, 
           dpi = dpi)
  }
  
  ## Map to UniProt etc.
  if ((!is.null(filename_ref)) && file.exists(filename_ref)) {
    filename2 <- paste0(fn_prefix, "_ref.", fn_suffix)

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
          factor(sub_location, levels = levels_subcellular),
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
    
    ggsave(file.path(filepath, filename2), width = width, height = height)
  }
}


#' Plots trends at a given \code{n_clust}.
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
    col_cluster = "cluster", data_type = "Trend", 
    pat = "Trend_[ONZ]_.*nclust\\d+\\.txt$", ext = "txt", 
    levels_subcellular = NULL, levels_subtype = NULL, tie_method = "none", 
    qt = .5, complete_cases = FALSE, panel_ids = NULL, filepath = NULL, 
    label_scheme_sub = NULL, theme = NULL, make_plot = TRUE, dots = NULL)
{
  ## (1) Load input data
  cl_id  <- sub(paste0(".*_nclust(\\d+)[^\\d]*\\.", ext), "\\1", df2) |>
    as.numeric()
  out_nm <- paste0(custom_prefix, fn_prefix, "_nclust", cl_id, ".", fn_suffix)
  src_path <- file.path(filepath, df2)
  
  if (!file.exists(src_path)) {
    stop("File not found: ", src_path)
  }
  
  # Warning msg: column 'group' not yet existed with data_type != "Trend"
  df <- readr::read_tsv(
    src_path, 
    col_types = cols(group = col_factor(), !!col_cluster := col_integer()))
  
  if (is.null(dim(df))) {
    stop("File contains not data: ", src_path)
  }
  
  if (!col_cluster %in% names(df)) {
    stop("Column 'cluster' not found in input data.")
  }
  
  message(paste("File loaded:", src_path))
  
  if (is.null(panel_ids)) {
    df <- df |>
      dplyr::mutate(
        !!col_cluster := as.integer(.data[[col_cluster]]), 
        !!col_cluster := factor(.data[[col_cluster]], 
                                levels = sort(unique(.data[[col_cluster]])))
      )
  } else {
    df <- df |>
      dplyr::filter(!!rlang::sym(col_cluster) %in% panel_ids) |>
      dplyr::mutate(!!rlang::sym(col_cluster) := 
                      factor(.data[[col_cluster]], levels = panel_ids))
  }
  
  # (x) 'anal_prnTrend.rda' overwritten with the latest call ->
  #     embedded 'col_group' etc. into the trend analysis output
  # (x) later to call from corresponding Protein_Trend_Z_nclust[...].rds
  col_group <- df[["col_group"]][1]
  col_order <- df[["col_order"]][1]
  col_fraction <- df[["col_fraction"]][1]
  col_subtype <- df[["col_subtype"]][1]
  
  df[["col_group"]] <- df[["col_order"]] <- df[["col_fraction"]] <- 
    df[["col_subtype"]] <- NULL

  ## (2.1) Set up plot parameters
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
      nrow <- ceiling(cl_id / ncol)
    } else {
      ncol <- ceiling(cl_id / nrow)
    }
  } else {
    nrow <- ceiling(cl_id / ncol)
  }
  
  dots$ymin <- NULL
  dots$ymax <- NULL
  dots$ybreaks <- NULL
  dots$ncol <- NULL
  dots$nrow <- NULL
  dots$color <- NULL
  dots$alpha <- NULL
  
  n_clust <- length(unique(df$cluster))
  if (is.null(width)) width <- n_clust * 8 / nrow + 2
  if (is.null(height)) height <- 8 * nrow
  
  dots$width <- NULL
  dots$height <- NULL

  ## (2.2) Data preprocessing
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
  
  ## (3.1) Set up the order of plot on the x-axis
  if (data_type == "Trend") {
    col_y <- "loc_score"
    col_x <- "group"

    df <- local({
      lev_x  <- label_scheme_sub |>
        dplyr::arrange(!!rlang::sym(col_order)) |>
        dplyr::select(!!rlang::sym(col_group)) |>
        unique() |>
        dplyr::pull(col_group)
      
      levs_df  <- levels(df[[col_x]])
      mis_lev_x <- levs_df[!levs_df %in% lev_x]
      
      if (length(mis_lev_x)) {
        if (length(mis_lev_x) > 12L) 
          mis_lev_x <- c(mis_lev_x[1:12], "...")
        
        stop("\n--- Mismatches in data levels ---\n\n", 
             "Levels in `", df2, "`:\n",
             paste(mis_lev_x, collapse = ", "), 
             "\n\n", 
             "Levels by `col_group = ", rlang::as_string(col_group), "`:\n", 
             paste(levs, collapse = ", "), "\n\n", 
             "??? Check for consistency in the setting of ", 
             "`anal_prnTrend(col_group = ...)` ", 
             "and `plot_prnTrend(col_group = ...)` for file `", 
             df2, "`.")
      }
      
      df <- df |>
        dplyr::filter(!!rlang::sym(col_x) %in% lev_x) |>
        filters_in_call(!!!filter2_dots) |>
        arrangers_in_call(!!!arrange2_dots) |>
        dplyr::mutate(!!col_x := factor(.data[[col_x]], levels = lev_x))
    })
  } else {
    col_y <- "loc_score"
    col_x <- "compartment"
    
    levels_subcellular <- sort(unique(df[[col_x]]))

    # To a long form
    df <- df |> 
      # tidyr::pivot_longer(-c(!!id, !!col_cluster), names_to = col_x, values_to = col_y) |>
      dplyr::arrange_at(c(col_cluster, col_x)) |>
      dplyr::mutate(
        !!col_x := factor(.data[[col_x]], levels = levels_subcellular))
  }
  
  # Median description of localization scores
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

  if (data_type == "Trend") {
    df_med <- df_med |>
      dplyr::left_join(
        df |>
          dplyr::select(
            dplyr::one_of(c(col_cluster, col_x, col_fraction, 
                            col_subtype))) |>
          unique(), 
        by = c(col_cluster, col_x))
    
    df_med <- df_med |>
      dplyr::mutate(
        !!col_fraction := 
          factor(!!rlang::sym(col_fraction), levels = levels_subcellular), )
  }
  
  # Handle multiple sub_type later...
  # e.g., df[[col_x]]: Veh_CP, TG_CP, df[[col_fraction]]: CP, CP
  # This allow the same fill color for 
  # different df[[col_x]] at the same df[[col_fraction]]
  if (is.null(col_fraction)) col_fraction <- col_x
  
  if (make_plot) {
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
      theme_classic(base_size = 11) +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.spacing = unit(1, "lines")
      )
    
    ggsave_dots <- 
      set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
    
    rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(out_nm)),
                      plot = p, 
                      width = width, 
                      height = height, 
                      !!!ggsave_dots)) |>
      rlang::eval_tidy()
  }

  list(
    df_med = df_med, 
    tie_method = tie_method, 
    col_fraction = col_fraction, 
    col_subtype = col_subtype, 
    col_cluster = col_cluster,
    col_x = col_x, 
    col_y = col_y
  )
  

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
  if (is.null(col_subtype)) {
    res <- tryCatch(
      classify_trends(df_med = df_med, 
                      tie_method = tie_method, 
                      col_fraction = col_fraction, 
                      col_subtype = col_subtype, 
                      col_cluster = col_cluster,
                      col_y = col_y),
      error = function(e) NULL)
  } else {
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


#' Classify trends.
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


#' Visualization of subcellular results.
#'
#' \code{plot_prnSubcellular} plots the subcellular purity of protein
#' expressions from \code{\link{anal_prnTrend}}.
#' 
#' @param dat_dir The working directory.
#' @param levels_subcellular The levels of subcellular fractions.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param pat The pattern of input filenames.
#' @param ext The extension of input filenames.
#' @param data_type The type of input data.
#' @param qt Not used. The quantile of localization score for thresholding
#'   subcellular compartment specificity. The default is to use median.
#' @param make_plot Logical; to bypass the plotting at FALSE.
#' @inheritParams anal_prnTrend
#' @inheritParams plot_prnTrend
#' @examples
#' if (FALSE) {
#'   cl <- 16; width = 10; height = 7.5;
#'   plot_prnSubcellular(
#'     df2 = paste0("Protein_Trend_Z_nclust", cl, ".txt"),
#'     col_order = Order,
#'     n_clust = cl,
#'     levels_subcellular = c("CP", "NP", "ChA"),
#'     ncol = 4,
#'     width = !!width,
#'     height = !!height,
#'   )
#' }
#' @export
plot_prnSubcellular <- function (
    dat_dir = NULL, 
    col_select = NULL, col_order = NULL, levels_subcellular = NULL, 
    levels_subtype = NULL, tie_method = "none", qt = .5, 
    # levels_col_group = NULL, 
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
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1L)

  ## Argument with single option and NULL default (quotation marks or not)
  col_select <- rlang::enexpr(col_select)
  col_order  <- rlang::enexpr(col_order)
  filepath   <- rlang::enexpr(filepath)
  filename   <- rlang::enexpr(filename)
  df2 <- rlang::enexpr(df2)
  
  # NULL or symbol -> length 0 or 1
  if (!is.character(col_select)) col_select <- as.character(col_select)
  if (!is.character(col_order)) col_order <- as.character(col_order)
  if (!is.character(filepath)) filepath <- as.character(filepath)
  if (!is.character(filename)) filename <- as.character(filename)
  if (!is.character(df2)) df2 <- as.character(df2)

  tie_method <- if (n_subcell_levs <- length(levels_subcellular)) {
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
              levels_subcellular = levels_subcellular, 
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
#' @param col_compartment The identifier of subcellular compartment.
#' @param dat_dir The current working directory.
#' @param filepath The filepath of inputs and outputs.
#' @param n_clust The number of clusters.
#' @param levels_subcellular The levels of subcellular fractions.
#' @param seed A seed for SVm.
#' @inheritParams anal_prnTrend col_subtype 
svmSubcell <- function (
    species = "human", id = "gene", 
    col_subtype = NULL, col_compartment = "compartment", 
    dat_dir = NULL, filepath = NULL, n_clust = NULL, 
    levels_subcellular = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), 
    seed = 42, ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  ## (1) Obtain expression data and reference data
  tempdata <- prepSubcell(
    species = species, id = id, dat_dir = dat_dir, filepath = filepath, 
    n_clust = n_clust, col_subtype = col_subtype, 
    levels_subcellular = levels_subcellular)
  df_expr  <- tempdata$df_expr
  df_ref   <- tempdata$df_ref
  saveRDS(df_expr, file.path(filepath, "df_expr.rds"))
  saveRDS(df_ref, file.path(filepath, "df_ref.rds"))
  rm(list = "tempdata")

  Xmat <- df_expr[df_ref[[id]], ]
  Xmat <- Xmat[match(rownames(Xmat), df_ref[[id]]), ]
  
  stopifnot(identical(rownames(Xmat), df_ref[[id]]))
  
  y_factor <- factor(df_ref$location)
  df_train <- as.data.frame(Xmat)
  df_train$location <- y_factor
  class_names <- levels(y_factor)
  
  ## (2) Training
  set.seed(seed)
  svm_model <- e1071::svm(
    location ~ .,
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
  entropies  <- -rowSums(prob_mat * log10(prob_mat))
  scores_loc <- prob_mat * (1 - entropies)
  entropies  <- tibble::tibble(!!id := names(entropies), entropy = entropies)
  
  scores_loc <- scores_loc |>
    data.frame() |>
    tibble::rownames_to_column(var = id) |> 
    tidyr::pivot_longer(
      -c(!!id), names_to = col_compartment, values_to = "loc_score")
  
  prob_mat <- prob_mat |>
    data.frame() |>
    tibble::rownames_to_column(var = id) |> 
    tidyr::pivot_longer(
      -c(!!id), names_to = col_compartment, values_to = "purity")
  
  ans <- prob_mat |>
    dplyr::left_join(entropies, by = id) |>
    dplyr::left_join(scores_loc, by = c(id, col_compartment))
  
  readr::write_tsv(ans, file.path(filepath, "svm_probs.tsv"))
  
  list(df_loc_score = ans, df_expr = df_expr)
}


#' Prepare subcellular data.
#' 
#' @param species Species in one of human, mouse and rat.
#' @param id Identifier.
#' @param dat_dir The current working directory.
#' @param filepath The filepath of inputs and outputs.
#' @param key_col The column key where data to be used.
#' @param levels_subcellular The levels of subcellular fractions.
#' @inheritParams anal_prnTrend col_subtype 
prepSubcell <- function (
    species = "human", id = "gene", dat_dir = NULL, filepath = NULL, 
    n_clust = NULL, key_col = "log10Int", col_subtype = NULL, 
    levels_subcellular = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
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
    
    df_ref <- df_ref |> 
      dplyr::mutate(location = sub("ERGIC", "ER", location))
  })

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
  
  
  # Check the number of unique values under col_subtype... "sample_type_unknown"
  if (length(col_subtype)) {
    
  }
  
  df_expr <- make_subcell_expr(
    df = df, id = id, key_col = key_col, col_subtype = col_subtype, 
    levels_subcellular = levels_subcellular, 
    scale = TRUE, match_colnames_to_levels_subcellular = TRUE)

  df_ref <- df_ref[df_ref[[id]] %in% rownames(df_expr), ]
  
  if (!nrow(df_ref)) {
    stop("No reference data retained after ID matching. \n", 
         "Check the correctionness of 'species = ", species, "'.")
  }
  
  df_ref <- df_ref[!is.na(df_ref[[id]]), ]
  
  list(df_expr = df_expr, df_ref = df_ref, abbr_species = abbr_species, 
       filenames_df2 = filenames_df2)
}


#' Prepare the matrix of expression data.
#'
#' @param df An input data frame.
#' @param id An identifier.
#' @param key_col The column key where values will be used for preparing
#'   expression data.
#' @param key_group The column key where group information will be found for
#'   preparing expression data.
#' @param levels_subcellular The levels of subcellular fractions.
#' @param scale Logical; to scale values by rowSums or not.
#' @param match_colnames_to_levels_subcellular Llogical; to match the column
#'   names of the expression data with the keys defined in
#'   \code{levels_subcellular} or not.
make_subcell_expr <- function (
    df = NULL, id = "gene", key_col = "log10Int", key_group = "group", 
    col_subtype = NULL, 
    levels_subcellular = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), 
    scale = TRUE, match_colnames_to_levels_subcellular = TRUE)
{
  if (!key_group %in% colnames(df)) {
    print(head(df))
    stop("Developer: column 'group' not found.")
  }
  
  # Then check if one or multiple values under 'col_subtype'
  
  if (!length(col_subtype)) {
    
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
    n_levs <- length(levels_subcellular)
    
    if (n_grps > n_levs) {
      stop("The number of groups is different to the number of ", 
           "fraction levels: \n", key_group, ": ", 
           paste(u_grps, collapse = ", "), "\n",
           "levels_subcellular", ": ", 
           paste(levels_subcellular, collapse = ", "))
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
  
  # Change column names to fraction names
  if (match_colnames_to_levels_subcellular) {
    colnames(df_expr) <- sapply(colnames(df_expr), function(col) {
      matched <- levels_subcellular[sapply(levels_subcellular, function(l) {
        # Handle edge case, e.g., colnames '5K' instead of 'XXX...5k'
        grepl(paste0("(^|[^0-9])", l, "$"), col)
      })]
      
      if (length(matched) == 1L) matched else col
    })
    
    df_expr <- df_expr[, levels_subcellular]
  }
  
  df_expr
}


#' Cluster SVM probabilities.
#'
#' @param dat_dir The working directory.
#' @param filepath A filepath for inputs and outputs.
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param n_clust The number of clusters.
#' @param choice A clustering method.
cluster_svmprobs <- function (
    dat_dir = NULL, filepath = NULL, id = "gene", n_clust = NULL, 
    choice = c("cmeans", "clara", "kmeans", "pam", "fanny"), 
    ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  choice <- match.arg(choice)
  dots <- rlang::enexprs(...)
  
  df_prob <- readr::read_tsv(file.path(filepath, "svm_probs.tsv"))
  
  df_expr <- local({
    df_expr <- df_prob[, c(id, "compartment", "loc_score")] |>
      tidyr::pivot_wider(
        id_cols = id, 
        names_from = c("compartment"), 
        values_from = "loc_score",
        # values_fill = 0, 
      )
    
    ids <- df_expr[[id]]
    df_expr <- as.matrix(df_expr[, colnames(df_expr) != id])
    rownames(df_expr) <- ids
    
    df_expr
  })
  
  ans_dots <- find_trend_m(
    dots, choice = choice, n_clust = n_clust, df_mean_log2r = df_expr)
  dots <- ans_dots[["dots"]]
  n_clust <- ans_dots[["n_clust"]]
  rm(list = c("ans_dots"))
  
  if (FALSE) {
    fns <- paste0("svmProb_Z_nclust", n_clust)
    res <- vector("list", length(n_clust))
    
    for (i in seq_along(n_clust)) {
      res[[i]] <- makeTrendRes(
        fns[[i]], choice = choice, dots = dots, id = id, df_mean_log2r = df_expr)
      readr::write_tsv(
        res[[i]], 
        file.path(filepath, paste0("Protein_svmProb_Z_nclust", i, ".tsv")))
    }
  }
  
  res <- lapply(
    paste0("svmProb_Z_nclust", n_clust), makeTrendRes, 
    choice = choice, dots = dots, id = id, df_mean_log2r = df_expr)
  
  mapply(function (rx, cl) {
    readr::write_tsv(
      dplyr::left_join(rx, df_prob, by = id), 
      file.path(filepath, paste0("Protein_svmProb_Z_nclust", cl, ".tsv")))
  }, res, n_clust)
  
  invisible(res)
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
#' @param plot_prob_clusters Logical; an option to plot clusters of
#'   probabilities.
#' @param seed A random seed.
#' @import ggplot2 RColorBrewer
#' 
#' @examples
#' \dontrun{
#' classify_subcellular(
#'   dat_dir = dat_dir, id = "gene", species = "mouse", n_clust = 16:36, 
#'   col_group = "Group", col_order = "Order", col_fraction = "Fraction", 
#'   levels_subcellular = 
#'     c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"))
#' }
#' 
#' @export
classify_subcellular <- function (
    dat_dir = NULL, filepath = NULL, id = "gene", species = "human", 
    n_clust = 16:36, levels_subcellular = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), 
    col_group = NULL, col_order = NULL, col_fraction = NULL, col_subtype = NULL, 
    min_prot_n_peps = 2L, plot_prob_clusters = TRUE, seed = 42L, ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  dots <- rlang::enexprs(...)

  ### (1) Generate expression data
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
  # Never NULL (since the default is not NULL)
  id <- rlang::as_string(rlang::enexpr(id))
  species <- rlang::as_string(rlang::enexpr(species))
  
  ## Special handling of character vectors
  levels_subcellular <- rlang::enexpr(levels_subcellular)
  levels_subcellular <- eval(levels_subcellular)
  n_levs <- length(levels_subcellular)
  if (n_levs == 1L) {
    if (levels_subcellular == "levels_subcellular") {
      stop("Use 'levels_subcellular = !!levels_subcellular' ", 
           "instead of 'levels_subcellular = levels_subcellular'.")
    }
    
    stop("Need more than one level for 'levels_subcellular'.")
  }
  rm(list = "n_levs")

  # Processed log2FC and intensity data
  # Need `!!` since `anal_prnTrend` is a UI function
  ans_trend <- anal_prnTrend(
    col_group = !!col_group, 
    col_order = !!col_order,
    col_fraction = !!col_fraction, 
    col_subtype = !!col_subtype, 
    cluster_data = FALSE, # to skip clustering steps
    # iter.max = 250, 
    n_clust = n_clust, 
    filter_by_npep = exprs(prot_n_pep >= !!min_prot_n_peps), )
  df_mean_log2r <- ans_trend[["log2R"]]
  df_mean_int   <- ans_trend[["Intensity"]]

  ### (2) Generate an SVM probability matrix
  ans_svm <- svmSubcell(
    species = species, id = id, dat_dir = dat_dir, filepath = filepath, 
    col_subtype = col_subtype, col_compartment = "compartment", 
    levels_subcellular = levels_subcellular, 
    n_clust = n_clust, seed = seed)

  ### (3) Classify SVM probabilities
  cluster_svmprobs(
    dat_dir = dat_dir, filepath = filepath, id = id, n_clust = n_clust, 
    iter.max = 250)
  
  ### (4) Optional: Plot clusters of SVM probabilities
  if (plot_prob_clusters) {
    plotTrend(id = id, 
              dat_dir = dat_dir, col_group = NULL, col_order = NULL, 
              label_scheme_sub = NULL, anal_type = "svmProb", 
              pat = "svmProb_[ONZ]_.*nclust\\d+\\.tsv$", ext = "tsv",  
              n_clust = n_clust, panel_ids = NULL, show_panel_ids = TRUE, 
              group_data_by = "ratio", scale_log2r = TRUE, 
              complete_cases = FALSE, impute_na = FALSE, 
              df2 = NULL, filepath = filepath, 
              filename = "Protein_svmProb_line_Z.png", theme = NULL, 
              ncol = 8, ymin = 0, ymax = 1, limitsize = FALSE, )
  }

  ### (5) Classify subcellular compartments by medium statistics
  plot_prnSubcellular(
    col_order = "colnames",
    q = .5, 
    ncol = 4,
    pat = "svmProb_[ONZ]_.*nclust\\d+\\.tsv$", 
    ext = "tsv", 
    data_type = "svmProb", 
    limitsize = FALSE,
    make_plot = TRUE, 
  )
}


#' Plot UMAP of sub-cellular scores.
#'
#' Inputs are expression data of normalized intensities.
#'
#' @param dat_dir The current working directory.
#' @param df_expr A matrix of expression data.
#' @param id An identifier.
#' @param n_clust Not used. The number of clusters.
#' @param levels_subcell_compartments The factor levels of subcellular
#'   compartments. Set \code{levels_subcell_compartments = NULL} for automatic
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
#' @export
plot_prnSubcell_UMAP <- function (
    dat_dir = NULL, df_expr = NULL, id = "gene", n_clust = 16:36, 
    levels_subcell_compartments = c(
      "Cytosol", "Nucleus", "ER", "PM", "Mitochondrion", "Golgi", "Lysosome",
      "Ribosome", "Stress.granule", "Proteasome", "Peroxisome", "Centrosome", 
      "Endosome", "Unknown"), 
    col_subtype = NULL, df2 = NULL, filepath = NULL, filename = NULL, 
    ext = "png", wrap_facets = TRUE, seed = 123, ...) 
{
  dots <- rlang::enexprs(...)
  if (is.null(dat_dir)) dat_dir  <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  
  ## Find input `df2` of SVM probabilities
  fns_svmprob <- find_trend_df2(
    df2 = df2, n_clust = NULL, pat = "Protein_svmProb_Z_nclust\\d+", 
    ext = "tsv", scale_log2r = TRUE, impute_na = FALSE, filepath = filepath)
  
  if (!length(fns_svmprob)) {
    warning("No SVM probability results found.", " Run 'anal_prnTrend' first.")
    return(NULL)
  }
  
  ## Load classification results
  df <- readr::read_tsv(file.path(filepath, fns_svmprob[[1]]))
  
  ## Compile subcellular factor levels
  possibles <- unique(df$compartment)
  
  df <- df |>
    dplyr::select(dplyr::one_of(id, "sub_location")) |>
    dplyr::filter(!is.na(sub_location)) |>
    unique()
  
  if (!length(levels_subcell_compartments)) {
    levels_subcell_compartments <- names(sort(table(df[["sub_location"]]), decreasing = TRUE))
    levels_missing <- possibles[!possibles %in% levels_subcell_compartments]
    levels_subcell_compartments <- c(levels_subcell_compartments, levels_missing)
    levels_subcell_compartments <- c(levels_subcell_compartments, "Unknown")
  } else {
    levels_missing <- NULL
  }
  
  rm(list = c("possibles", "levels_missing"))
  
  df <- df |>
    dplyr::mutate(
      sub_location = factor(sub_location, levels = levels_subcell_compartments))
  
  if (is.null(filename)) {
    fn_prefix <- paste0(tools::file_path_sans_ext(fns_svmprob[[1]]), "_umap")
    fn_suffix <- ext
    filename  <- paste0(fn_prefix, ".", fn_suffix)
  } else {
    fn_prefix <- tools::file_path_sans_ext(filename)
    fn_suffix <- tools::file_ext(filename)
  }
  
  ## Load the matrix of normalized intensity across fractions
  # Recompile df_expr from scratch
  if (FALSE) {
    df_expr <- local({
      filenames_df2 <- tryCatch(
        find_trend_df2(
          df2 = NULL, n_clust = n_clust, scale_log2r = TRUE, 
          impute_na = FALSE, filepath = filepath), 
        error = function (e) NULL)
      
      if (length(filenames_df2)) {
        df_trend <- readr::read_tsv(file.path(filepath, filenames_df2[[1]]))
        col_fraction <- df_trend[["col_fraction"]][[1]]
        col_subtype <- df_trend[["col_subtype"]][[1]]
        key_col <- "log10Int"
        
        df_trend <- df_trend |>
          dplyr::select(dplyr::one_of(c(id, col_fraction, key_col)))
        
        levels_subcellular <- unique(df_trend[[col_fraction]])
        
        df_expr <- make_subcell_expr(
          df = df_trend, id = id, key_col = key_col, key_group = col_fraction, 
          levels_subcellular = levels_subcellular, col_subtype = col_subtype, 
          scale = TRUE, match_colnames_to_levels_subcellular = TRUE)
      } else {
        df_expr <- NULL
      }
      
      df_expr
    })
    
    if (isnull(df_expr)) {
      if (file.exists(fn_expr <- file.path(filepath, "df_expr.rds"))) {
        df_expr <- readRDS(fn_expr)
      } else {
        warning("Expression file not found: ", file.path(filepath, "df_expr.rds"))
        return(NULL)
      }
    }
  }
  
  if (file.exists(fn_expr <- file.path(filepath, "df_expr.rds"))) {
    df_expr <- readRDS(fn_expr)
  } else {
    warning("Expression file not found: ", file.path(filepath, "df_expr.rds"))
    return(NULL)
  }
  rm(list = "fn_expr")
  
  df_umap <- local({
    set.seed(seed)
    res_umap <- umap::umap(df_expr)
    df_umap  <- res_umap$layout[, 1:2]
    colnames(df_umap) <- paste0("UMAP", 1:2)
    rownames(df_umap) <- rownames(df_expr)
    df_umap <- tibble::rownames_to_column(data.frame(df_umap), id)
    readr::write_tsv(df_umap, file.path(filepath, "umap.tsv"))
    
    df_umap
  })
  
  df_plot <- tibble::tibble(df_umap) |>
    dplyr::left_join(df[, c("gene", "sub_location")], by = id)
  
  df_plot <- df_plot |>
    dplyr::mutate(
      sub_location = as.character(sub_location), 
      sub_location = ifelse(is.na(sub_location), "Unknown", sub_location), 
      sub_location = factor(sub_location, levels = levels_subcell_compartments), ) |>
    unique()
  
  my_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(
    length(levels(df_plot$sub_location)))
  
  ncol   <- dots$ncol
  nrow   <- dots$nrow
  size   <- dots$size
  alpha  <- dots$alpha
  width  <- dots$width
  height <- dots$height
  dpi    <- dots$dpi
  
  n_subtypes <- 1L # temporary
  
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
  
  if (wrap_facets) {
    p <- p + 
      facet_wrap(~ sub_location, ncol = ncol, labeller = label_value, drop = FALSE)
  }
  
  ggsave(file.path(filepath, paste0(fn_prefix, ".", fn_suffix)), 
         width = width, height = height, dpi = dpi)
  
  invisible(df_plot)
}


#' Update protein tables of \code{Protein.txt} etc. by subcellular locations.
#' 
#' Not used.
#' 
#' @param id Identifier.
#' @inheritParams anal_prnTrend
update_subcellular <- function (id = "gene", impute_na = FALSE, 
                                scale_log2r = TRUE, ...) 
{
  dots <- rlang::enexprs(...)
  
  dat_dir  <- get_gl_dat_dir()
  filepath <- file.path(dat_dir, "Protein", "trend")
  
  label_scheme_sub <- load_ls_group(dat_dir)
  
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
  df_trend <- df_trend[, c("gene", "loc_score", "sub_location", col_subtype)]
  
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
#' @param id Identifier.
#' @param label_scheme_sub Metadata.
#' @param col_fraction The column key in metadata linking to subcellular
#'   values.
#' @param col_subtye The column key metadata linking to sample sub-types.
hupdate_subcellular <- function (df_trend, dat_dir, id = "gene", 
                                 label_scheme_sub, 
                                 col_fraction = NULL, col_subtype = NULL) 
{
  path_prn <- file.path(dat_dir, "Protein", "Protein.txt")
  path_pep <- file.path(dat_dir, "Peptide", "Peptide.txt")
  file.copy(path_prn, file.path(dat_dir, "Protein", "Protein_bf_subcell.txt"))
  file.copy(path_pep, file.path(dat_dir, "Peptide", "Peptide_bf_subcell.txt"))
  
  df_prn <- readr::read_tsv(path_prn)
  df_pep <- readr::read_tsv(path_pep)
  
  metas <- split(label_scheme_sub, 
                 label_scheme_sub[, c(col_fraction, col_subtype)])
  
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
  }
  
  readr::write_tsv(df_prn, path_prn)
  readr::write_tsv(df_pep, path_pep)
  
  invisible(df_prn)
}


