#' Subcellular plots.
#' 
#' @param levels_subcellular The levels of subcellular fractions.
#' @param panel_ids Panel IDs for plotting.
#' @inheritParams plot_prnTrend
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr purrr ggplot2 RColorBrewer
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotSubcellular <- function(
    id, col_group = NULL, col_order = NULL, levels_subcellular = NULL, 
    levels_subtype = NULL, label_scheme_sub, n_clust = NULL, panel_ids = NULL, 
    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
    df2 = NULL, filepath, filename, theme, ...) 
{
  if (!nrow(label_scheme_sub)) {
    stop("Empty metadata.")
  }
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

  # Find input `df2`
  df2 <- find_trend_df2(
    df2 = df2, n_clust = n_clust, scale_log2r = scale_log2r, 
    impute_na = impute_na, filepath = filepath)
  
  # Prepare output file name
  custom_prefix <- if (id %in% c("pep_seq", "pep_seq_mod")) {
    purrr::map_chr(df2, ~ gsub("(.*_{0,1})Peptide_Trend.*", "\\1", .x))
  }
  else if (id %in% c("prot_acc", "gene")) {
    purrr::map_chr(df2, ~ gsub("(.*_{0,1})Protein_Trend.*", "\\1", .x))
  }
  else {
    stop("Unknown id = ", id)
  }

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  ## Plot data
  if (is.null(theme)) {
    theme <- theme_bw()
  }
  
  tempdata <- mapply(
    plotSubcellular_sub, 
    df2 = df2, 
    custom_prefix = custom_prefix, 
    MoreArgs = list(
      fn_prefix = fn_prefix, 
      fn_suffix = fn_suffix, 
      df = df, 
      id = id, 
      col_group = rlang::enexpr(col_group), 
      col_order = rlang::enexpr(col_order), 
      levels_subcellular = levels_subcellular,
      levels_subtype = levels_subtype, 
      complete_cases = complete_cases, 
      panel_ids = panel_ids, 
      filepath = filepath, 
      label_scheme_sub = label_scheme_sub, 
      theme = theme, 
      dots = dots), 
    SIMPLIFY = FALSE, USE.NAMES = TRUE
  )
  
  res_score_co <- lapply(tempdata, `[[`, "res")
  col_group <- tempdata[[1]][["col_group"]]
  col_order <- tempdata[[1]][["col_order"]]
  col_subcellular <- tempdata[[1]][["col_subcellular"]]
  col_subtype <- tempdata[[1]][["col_subtype"]]

  ## Update subcellular locations
  res_score_co <- res_score_co |>
    dplyr::bind_rows(.id = "n_clust") |>
    dplyr::mutate(
      n_clust = 
        as.integer(gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", n_clust))) |>
    dplyr::rename(subcellular = subcellular.) |> 
    dplyr::arrange(subcellular, score)
  
  # readr::write_tsv(res_score_co, file.path(filepath, "all_loc_scores.tsv"))
  
  res_score_co <- res_score_co |>
    dplyr::group_by_at(c("subcellular", col_subtype)) |>
    dplyr::slice(1L) |>
    dplyr::rename(score_co = score, ) |>
    dplyr::ungroup() |>
    dplyr::arrange(score_co)

  join_cols <- c("subcellular", col_subtype)
  names(join_cols) <- c(col_subcellular, col_subtype)
  
  dfs <- lapply(df2, function (fn) {
    df <- readr::read_tsv(file.path(filepath, fn)) |>
      dplyr::left_join(
        res_score_co[, c("subcellular", "score_co", col_subtype)], 
        by = join_cols)
    
    # Keep the original data structure; do not filter
    # may split by "score_co" if no ties
    df <- lapply(split(df, df[["group"]]), function (x) {
      frac <- x[[col_subcellular]][[1]]
      sco  <- x[["score_co"]][[1]]
      
      x <- x |>
        dplyr::mutate(
          sub_location = ifelse(loc_score >= sco, frac, NA_character_))
    }) |>
      dplyr::bind_rows() |>
      dplyr::select(-c("score_co"))
    
    if (FALSE) {
      dfx <- df[, c("gene", "sub_location")] |> 
        dplyr::filter(!is.na(sub_location)) |>
        unique() |> 
        dplyr::group_by(gene) |> 
        dplyr::mutate(
          sub_location = factor(sub_location, levels = levels_subcellular), 
          N = dplyr::row_number(), )

      ggplot(dfx, aes(x = sub_location, fill = sub_location)) +
        geom_bar(width = 0.6) + 
        theme_minimal() +
        labs(
          title = NULL,
          x = "Subcellular Location",
          y = "Number of Genes",
          fill = NULL
        ) +
        theme(axis.text.x = element_text(hjust = 1), 
              panel.grid  = element_blank(), ) 
      
      ggsave(file.path(filepath, "prnSubcellular_counts.png"), 
             width = 4, height = 3)

      ggplot(dfx, aes(x = factor(N), fill = factor(N))) +
        geom_bar(width = 0.6) +
        theme_minimal() +
        labs(
          title = NULL,
          x = "Number of Locations per Gene",
          y = "Number of Genes",
          fill = NULL
        ) +
        theme(
          axis.text.x = element_text(),
          panel.grid  = element_blank(), 
          legend.position = "none"
        )
      
      ggsave(file.path(filepath, "subcellular_counts_per_gene.png"), 
             width = 3, height = 3)
    }
    
    # Update data
    readr::write_tsv(df, file.path(filepath, fn))
    
    df
  })

  invisible(dfs)
}


#' Plot UMAP of subcellular scores.
#'
#' @param filename_ref Optional; a file name with prepending path to a reference
#'   subcellular annotation file, e.g., UniProt and GO intersection.
#' @inheritParams plot_prnTrend
#' @inheritParams plot_prnSubcellular
#' @export
plot_prnSubcellular_UMAP <- function (
    df2 = NULL, levels_subcellular = NULL, levels_subtype = NULL, 
    filename_ref = NULL, 
    scale_log2r = TRUE, impute_na = FALSE, filename = NULL, ...) 
{
  dots <- rlang::enexprs(...)
  
  dat_dir  <- get_gl_dat_dir()
  filepath <- file.path(dat_dir, "Protein", "trend")
  
  ## Find input `df2`
  df2 <- find_trend_df2(
    df2 = df2, n_clust = NULL, scale_log2r = scale_log2r, 
    impute_na = impute_na, filepath = filepath)
  
  if (!length(df2)) {
    stop("No trend results found.", " Run 'anal_prnTrend' first.")
  }
  
  if (is.null(filename)) {
    fn_prefix <- paste0(tools::file_path_sans_ext(df2[[1]]), "_umap")
    fn_suffix <- "png"
    filename  <- paste0(fn_prefix, ".", fn_suffix)
  } else {
    fn_prefix <- tools::file_path_sans_ext(filename)
    fn_suffix <- tools::file_ext(filename)
  }
  
  df <- readr::read_tsv(file.path(filepath, df2[[1]]))
  col_subcellular <- df[["col_subcellular"]][[1]]
  col_subtype <- df[["col_subtype"]][[1]]

  # col_subcellular encodes sample fraction; 
  # sub_location: classification results
  df_umap <- df[, c("gene", "loc_score", "sub_location", col_subtype)] |>
    dplyr::filter(!is.na(sub_location)) |>
    dplyr::arrange(gene) |>
    unique() |>
    tidyr::pivot_wider(
      id_cols = "gene", 
      names_from = c("sub_location", col_subtype), 
      values_from = "loc_score",
      values_fill = 0, 
    )
  
  res_umap <- local({
    dfx <- data.frame(df_umap)
    rownames(dfx) <- df_umap$gene
    dfx$gene <- NULL
    
    set.seed(42)
    res <- umap::umap(dfx)
  })
  
  umap_df <- data.frame(
    UMAP1 = res_umap$layout[, 1],
    UMAP2 = res_umap$layout[, 2]
  ) |>
    tibble::rownames_to_column("gene")
  
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
  
  if (is.null(width))  width <- n_subtypes * 8 / nrow + 2
  if (is.null(height)) height <- 8 * nrow

  if (is.null(ncol))  ncol <- 2
  if (is.null(nrow))  nrow <- 1
  if (is.null(size))  size <- .02
  if (is.null(alpha)) alpha <- .5
  
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
  
  df_plot <- df |> 
    dplyr::filter(!is.na(sub_location)) |>
    dplyr::mutate(
      !!col_subtype := 
        factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
      !!col_subcellular := 
        factor(!!rlang::sym(col_subcellular), levels = levels_subcellular),
      sub_location = 
        factor(sub_location, levels = levels_subcellular),
    ) |>
    dplyr::left_join(umap_df, by = "gene")
  
  p <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = sub_location)) +
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
    facet_wrap(
      vars(!!rlang::sym(col_subtype[1])), nrow = nrow, labeller = label_value) + 
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave(file.path(filepath, filename), width = width, height = height)
  
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
      dplyr::rename(sub_location = compartment_2) |>
      dplyr::select(c("gene", "sub_location")) |>
      unique() |>
      dplyr::filter(gene %in% unique(df$gene))
    
    df_plot2 <- df_plot |>
      dplyr::select(-one_of("sub_location")) |>
      dplyr::left_join(df_ref, by = "gene") |>
      dplyr::filter(!is.na(sub_location)) |>
      dplyr::select(
        one_of(c("gene", "UMAP1", "UMAP2", "sub_location", col_subtype))) |>
      unique() |>
      dplyr::mutate(
        !!col_subtype := 
          factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
        sub_location = 
          factor(sub_location, levels = levels_subcellular),
      )

    p2 <- ggplot(df_plot2, aes(x = UMAP1, y = UMAP2, color = sub_location)) +
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
#' @param dots Variable arguments.
#' @inheritParams plotTrend
#' @inheritParams plotSubcellular
#' @importFrom magrittr %>% %T>% 
#' @import dplyr ggplot2 RColorBrewer
plotSubcellular_sub <- function (
    df2 = NULL, custom_prefix = NULL, fn_prefix = NULL, fn_suffix = NULL, 
    df = NULL, id = "gene", col_group = NULL, col_order = NULL, 
    levels_subcellular = NULL, levels_subtype = NULL, 
    complete_cases = FALSE, panel_ids = NULL, filepath = NULL, 
    label_scheme_sub = NULL, theme = NULL, dots = NULL)
{
  cl_id  <- as.numeric(gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", df2))
  out_nm <- paste0(custom_prefix, fn_prefix, "_nclust", cl_id, ".", fn_suffix)
  src_path <- file.path(filepath, df2)
  
  df <- tryCatch(
    readr::read_tsv(
      src_path, 
      col_types = cols(group = col_factor(), cluster = col_integer())), 
    error = function(e) NA)
  
  if (is.null(dim(df))) {
    stop(paste("Non-exist file or directory:", src_path))
  }
  
  message(paste("File loaded:", src_path))
  
  if (is.null(panel_ids)) {
    df <- df |>
      dplyr::mutate(
        cluster = as.integer(cluster),
        cluster = factor(cluster, levels = sort(unique(cluster)))
      )
  } else {
    df <- df |>
      dplyr::filter(cluster %in% panel_ids) |>
      dplyr::mutate(cluster = factor(cluster, levels = panel_ids))
  }
  
  # (x) 'anal_prnTrend.rda' overwritten with the latest call ->
  #     embedded 'col_group' etc. into the trend analysis output
  # (x) later to call from corresponding Protein_Trend_Z_nclust[...].rds
  col_group <- df[["col_group"]][1]
  col_order <- df[["col_order"]][1]
  col_subcellular <- df[["col_subcellular"]][1]
  col_subtype <- df[["col_subtype"]][1]
  
  if (FALSE) {
    col_group <- rlang::as_string(rlang::enexpr(col_group))
    col_order <- rlang::as_string(rlang::enexpr(col_order))
    col_subcellular <- match_call_arg(anal_prnTrend, "col_subcellular")
    col_subcellular <- rlang::as_string(rlang::enexpr(col_subcellular))
    col_subtype <- match_call_arg(anal_prnTrend, "col_subtype")
    col_subtype <- rlang::as_string(rlang::enexpr(col_subtype))
  }
  
  if (FALSE) {
    fn_df2_rds <- paste0(tools::file_path_sans_ext(df2), ".rds")
    fn_sc_lookup <- "subcellular_lookup.tsv"
    if (file.exists(fi_sc <- file.path(filepath, fn_df2_rds))) {
      sc_lookup <- readRDS(fi_sc)[["sc_lookup"]]
    } else if (file.exists(fi_sc <- file.path(filepath, fn_sc_lookup))) {
      sc_lookup <- readr::read_tsv(fi_sc)
    } else {
      sc_lookup <- NULL
    }
    rm(list = c("fn_df2_rds", "fn_sc_lookup"))
  }

  levs <- label_scheme_sub |>
    dplyr::arrange(!!rlang::sym(col_order)) |>
    dplyr::select(!!rlang::sym(col_group)) |>
    unique() |>
    dplyr::pull(col_group)
  
  local({
    levs_df  <- levels(df$group)
    mis_levs <- levs_df[!levs_df %in% levs]
    
    if (length(mis_levs)) {
      if (length(mis_levs) > 12L) 
        mis_levs <- c(mis_levs[1:12], "...")
      
      stop("\n--- Mismatches in data levels ---\n\n", 
           "Levels in `", df2, "`:\n",
           paste(mis_levs, collapse = ", "), 
           "\n\n", 
           "Levels by `col_group = ", rlang::as_string(col_group), "`:\n", 
           paste(levs, collapse = ", "), "\n\n", 
           "??? Check for consistency in the setting of ", 
           "`anal_prnTrend(col_group = ...)` ", 
           "and `plot_prnTrend(col_group = ...)` for file `", 
           df2, "`.")
    }
  })
  
  if (length(dots) && any(grepl("^filter_", names(dots)))) {
    stop("Primary `filter_` depreciated; use secondary `filter2_`.")
  }
  
  lang_dots     <- dots[unlist(lapply(dots, is.language))]
  filter2_dots  <- lang_dots[grepl("^filter2_", names(lang_dots))]
  arrange2_dots <- lang_dots[grepl("^arrange2_", names(lang_dots))]
  dots          <- dots[!dots %in% c(filter2_dots, arrange2_dots)]
  
  df <- df |>
    dplyr::filter(group %in% levs) |>
    filters_in_call(!!!filter2_dots) |>
    arrangers_in_call(!!!arrange2_dots) |>
    dplyr::mutate(group = factor(group, levels = levs))

  if (complete_cases) {
    df <- df[complete.cases(df), ]
  }
  
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
  
  # Median description of localization scores
  # User defined column `group" provides the finest granularity including 
  #  subcellular, sample type, etc.
  df_med <- lapply(split(df, df[["cluster"]]), function (dfx) {
    dfx |>
      dplyr::group_by(group) |>
      dplyr::summarise(
        purity  = median(purity, na.rm = TRUE), 
        entropy = median(entropy, na.rm = TRUE), 
        score   = median(loc_score, na.rm = TRUE), 
      )
  }) |>
    dplyr::bind_rows(.id = "cluster") |>
    dplyr::mutate(
      cluster = as.integer(cluster),
      cluster = factor(cluster, levels = sort(unique(cluster)))
    )
  
  df_med <- df_med |>
    dplyr::left_join(
      unique(df[, c("cluster", "group", col_subcellular, col_subtype)]), 
      by = c("cluster", "group")) |>
    dplyr::mutate(
      !!col_subcellular := 
        factor(!!rlang::sym(col_subcellular), levels = levels_subcellular)
    )
  
  p <- ggplot(
    df_med, aes(x = group, y = score, fill = !!rlang::sym(col_subcellular))) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ cluster, ncol = ncol, scales = "free_y") +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    scale_fill_viridis_d(option = "plasma", end = 0.85) + 
    scale_y_continuous(
      labels = scales::label_number(accuracy = .1),
      breaks = scales::breaks_width(.2)
    ) +
    labs(
      x = NULL,
      y = "Median Localization Score",
      fill = "Subcellular"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
      panel.spacing = unit(1, "lines")
    )

  ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
  
  rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(out_nm)),
                    plot = p, 
                    width = width, 
                    height = height, 
                    !!!ggsave_dots)) |>
    rlang::eval_tidy()
  
  ## Classify trends
  if (FALSE) {
    res <- tryCatch(
      classify_trends(df_med = df_med, col_subcellular = col_subcellular, 
                      col_subtype = col_subtype),
      error = function(e) NULL)
  }

  res <- lapply(split(df_med, df_med[[col_subtype]]), function (dfx) {
    tryCatch(
      classify_trends(df_med = dfx, col_subcellular = col_subcellular, 
                      col_subtype = col_subtype),
      error = function(e) NULL)
  }) |>
    dplyr::bind_rows(.id = col_subtype) |>
    dplyr::ungroup() |>
    dplyr::arrange(subcellular., score)
  
  # Currently use the lowest score across all cell states of `col_subtype`
  if (FALSE) {
    res <- res |>
      dplyr::group_by(subcellular.) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::select(-dplyr::one_of(col_subtype))
  }

  list(res = res, 
       col_group = col_group,
       col_order = col_order,
       col_subcellular = col_subcellular, 
       col_subtype = col_subtype)
}


#' Classify trends.
#' 
#' @param df_med A median description of trend scores.
#' @inheritParams anal_prnTrend
classify_trends <- function (df_med = NULL, col_subcellular, col_subtype)
{
  # Assign the primary subcellular fraction for each cluster
  df_pri <- df_med[, c("cluster", "purity", "entropy", "score", 
                       col_subcellular)] |>
    dplyr::group_by(cluster) |>
    dplyr::arrange(cluster, -score) |>
    dplyr::ungroup()
  
  df_pri <- lapply(
    split(df_pri, interaction(df_pri[["cluster"]], df_pri[[col_subcellular]], 
                              drop = TRUE)), 
    function(x) x[1, , drop = FALSE]
  ) |> 
    dplyr::bind_rows()
  
  df_pri <- df_pri |>
    dplyr::group_by(cluster) |>
    dplyr::arrange(cluster, -score) |>
    dplyr::ungroup()
  
  # [x] May assign to both fractions if < .01, but not contribute to the thresholds
  # [x] Or count the percent of missing values per group and assign the cluster to 
  #  the most complete groups; must be >= 50% values...
  oks <- sapply(split(df_pri, df_pri$cluster), function (x) {
    nrow <- nrow(x)
    fracs <- x[[col_subcellular]]
    scos <- x[["score"]]
    if (nrow >= 2L) scos[[1]] - scos[[2]] >= .01 else TRUE
  })
  
  df_pri <- df_pri |>
    dplyr::group_by(cluster) |>
    # dplyr::slice(if (dplyr::n() >= 2L) 2L else 1L) |>
    dplyr::slice(1L)
  df_pri[[col_subcellular]][!oks] <- NA_character_

  df_pri <- df_pri |>
    dplyr::ungroup() |>
    dplyr::rename(subcellular. := !!col_subcellular)

  # Find the least probable cluster under each subcellular group;
  # or match the value of the next most abundant fraction
  df_pri_min <- df_pri |>
    dplyr::arrange(subcellular., score) |>
    dplyr::group_by(subcellular.) |>
    dplyr::slice(1L) |>
    dplyr::filter(!is.na(subcellular.))
  
  return(df_pri_min)
  
  
  # Elbow, KNN, ... 
  # Or data driven (heat map correlation between clusters)
  score_cos <- df_pri_min[["score"]] * .5
  score_cos <- pmax(score_cos, min(df_pri_min[["score"]]))
  
  ans_interfacial <- list()
  for (i in seq_len(nrow(df_pri_min))) {
    df_pri_min_i <- df_pri_min[i, ]
    score_min <- df_pri_min_i[["score"]]
    subcell_i <- df_pri_min_i$subcellular.
    cl_i <- df_pri_min_i[["cluster"]]
    score_co <- score_cos[[i]]

    # Scores of the current subcellular fraction at each clusters
    dfi <- df_med |>
      dplyr::filter(.data[[col_subcellular]] == subcell_i, 
                    # .data[["cluster"]] != cl_i, 
                    ) |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(
        purity = median(purity, na.rm = TRUE), 
        entropy = median(entropy, na.rm = TRUE), 
        score = median(score, na.rm = TRUE), 
      ) |>
      dplyr::arrange(-score)
    
    # Finds the turning point of 
    # plot(dfi$entropy ~ dfi$purity)
    
    dfi <- dfi |>
      dplyr::left_join(
        df_pri[, c("cluster", "subcellular.")], by = "cluster") |> 
      dplyr::mutate(delta = .data[["score"]] - score_min) 
    
    # plot(dfi$delta)
    
    dfi <- dfi |>
      # dplyr::arrange(subcellular., -delta) |>
      dplyr::filter(
        # purity >= 1 / (length(levels_subcellular) + 2L), 
        # entropy <= .5, 
        score >= .1, 
        delta >= -.1)
    
    ans_interfacial[[subcell_i]]  <- if (nrow(dfi)) {
      dfi[["cluster"]][[which.min(dfi[["delta"]])]]
    } else {
      0L
    }
  }
  
  names(ans_interfacial) <- df_pri_min[["subcellular."]]
  df_pri_min[["cluster_next"]] <- unlist(ans_interfacial, use.names = FALSE)
  
}


#' Visualization of subcellular results.
#' 
#' \code{plot_prnSubcellular} plots the subcellular purity of protein expressions
#' from \code{\link{anal_prnTrend}}.
#' 
#' @param levels_subcellular The levels of subcellular fractions.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
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
#' 
#' @export
plot_prnSubcellular <- function (
    col_select = NULL, col_order = NULL, levels_subcellular = NULL, 
    levels_subtype = NULL, 
    # levels_col_group = NULL, 
    n_clust = NULL, panel_ids = NULL, scale_log2r = TRUE, 
    complete_cases = FALSE, impute_na = FALSE, df2 = NULL, 
    filename = NULL, theme = NULL, ...) 
{
  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)
  
  dir.create(file.path(get_gl_dat_dir(), "Protein/Trend/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  id <- tryCatch(
    match_call_arg(normPSM, group_pep_by), 
    error = function(e) NA)
  
  if (is.na(id)) {
    id <- tryCatch(
      match_call_arg(makeProtDIANN, group_pep_by), 
      error = function(e) NA)
  }
  
  if (is.na(id)) {
    id <- "gene"
  }
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  filename <- rlang::enexpr(filename)
  df2 <- rlang::enexpr(df2)
  
  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_group = NULL, 
            col_order = !!col_order,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = NULL, 
            df2 = !!df2, 
            filepath = NULL, 
            filename = !!filename,
            anal_type = "Subcellular_plot")(
              levels_subcellular = levels_subcellular, 
              levels_subtype = levels_subtype, 
              n_clust = n_clust, panel_ids = panel_ids, theme = theme, ...)
}

