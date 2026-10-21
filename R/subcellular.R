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
    levels_ref_compartments = NULL, tie_method = "none", 
    label_scheme_sub = NULL, n_clust = NULL, panel_ids = NULL, 
    pat = "Trend_[ONZ]_.*nclust\\d+\\.txt$", ext = "txt", data_type = "Trend", 
    scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
    df2 = NULL, filepath = NULL, filename = NULL, 
    type_supervision = "both", qt = .5, 
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

  # Obtain quantile-summarized localization scores by each n_clust
  tempsub <- mapply(
    plotSubcellular_sub, 
    df2 = df2, 
    custom_prefix = custom_prefix, 
    MoreArgs = list(
      fn_prefix = fn_prefix, 
      fn_suffix = fn_suffix, 
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
      tie_method = tie_method, # no used
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
  col_subtype  <- tempsub[[1]][["col_subtype"]]
  col_x <- tempsub[[1]][["col_x"]]
  col_y <- tempsub[[1]][["col_y"]]

  ## Classify trends at each n_clust
  tempdata <- lapply(
    df_meds, 
    hassign_compartment_to_trend_cluster, 
    col_subtype = col_subtype, 
    col_fraction = col_fraction, 
    levels_ref_compartments = levels_ref_compartments, 
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
  # Some compartments may have no corresponding primary clusters. 
  res_score_co <- res_score_co |>
    dplyr::bind_rows(.id = "n_clust") |>
    dplyr::mutate(
      n_clust = 
        as.integer(gsub(paste0(".*_nclust(\\d+)[^\\d]*\\.", ext, "$"), "\\1", 
                        n_clust))) |>
    dplyr::arrange(!!rlang::sym(col_fraction), !!rlang::sym(col_y))
  
  readr::write_tsv(res_score_co, file.path(filepath, "all_loc_scores.tsv"))
  
  res_score_co <- res_score_co |>
    dplyr::group_by_at(c(col_fraction, col_subtype)) |>
    dplyr::slice(1L) |>
    dplyr::rename(score_co = !!rlang::sym(col_y), ) |>
    dplyr::ungroup() |>
    dplyr::arrange(score_co)
  
  message(
    "Leading subcellular compartments with corresponding sample fractions:\n")
  
  print(
    res_score_co |>
      dplyr::select(dplyr::one_of(col_subtype, col_fraction)) |>
      unique()
  )
  
  dfs <- make_sub_location(
    df2 = df2, res_score_co = res_score_co, type_supervision = type_supervision, 
    id = id, dat_dir = dat_dir, filepath = filepath, data_type = data_type, 
    col_fraction = col_fraction, col_subtype = col_subtype)

  invisible(dfs)
}


#' Helper making the column sub_location in trend results.
#' 
#' @param df2 File names of trend results.
#' @param res_score_co Score cut-offs of classes (with leading fractions).
#' @param type_supervision Only for uses with supervised learning: \code{sole},
#'   supervised learning only; \code{both}, both supervised and unsupervised
#'   learning.
#' @param id An identifier.
#' @param dat_dir The current working directory.
#' @param filepath The file path to trend results.
#' @param data_type The type of data, either Trend or svmProb.
#' @param col_fraction The column key to sample fractions.
#' @param col_subtype The column key to sample sub-types.
make_sub_location <- function (
    df2, res_score_co, type_supervision = "both", id = "gene", 
    dat_dir = NULL, filepath = NULL, data_type = "Trend", 
    col_fraction = NULL, col_subtype = NULL)
{
  join_cols <- c(col_fraction, col_subtype)
  names(join_cols) <- c(col_fraction, col_subtype)
  
  is_supervised <- data_type %in% c("svmProb")
  
  dfs <- lapply(df2, function (fn) {
    ## 1. Add unsupervised learning results
    df <- readr::read_tsv(file.path(filepath, fn))
    df[["score_co"]] <- NULL # just in case
    
    df <- df |>
      dplyr::left_join(
        res_score_co |>
          dplyr::select(
            dplyr::all_of(
              c(col_fraction, "score_co", col_subtype, "is_leading_fraction")
            )
          ), 
        by = join_cols)
    
    # Keep the original data structure; do not filter
    # may split by "score_co" if no ties are allowed
    df <- lapply(split(df, df[[col_fraction]]), function (x) {
      frac <- x[[col_fraction]][[1]]
      sco  <- x[["score_co"]][[1]]
      
      # e.g. is.na if without a leading organelle being assigned
      if (isTRUE(sco)) {
        x[["sub_location"]] <- NA_character_
      } else {
        x <- x |>
          dplyr::mutate(
            sub_location = ifelse(loc_score >= sco, frac, NA_character_))
      }
      
      x
    }) |>
      dplyr::bind_rows() |> 
      dplyr::select(-c("score_co"))
    
    ## 2. Apply SVM results for proteins without a leading organelle being assigned
    if (is_supervised) {
      df_fst <- df |>
        dplyr::ungroup() |>
        dplyr::mutate(index. = dplyr::row_number()) |>
        dplyr::arrange(-loc_score) |>
        dplyr::group_by_at(id) |>
        dplyr::mutate(is_first = dplyr::row_number() == 1L) |>
        dplyr::ungroup() |>
        dplyr::arrange(index.) |>
        dplyr::select(-index.)
      df_fst[["svm"]][!df_fst[["is_first"]]] <- NA_character_
      df_fst[["is_first"]] <- NULL
      
      if (type_supervision == "sole") {
        # Overwrite the 'sub_location' findings from unsupervised learning
        df[["sub_location"]] <- df_fst[["svm"]]
        readr::write_tsv(df, file.path(filepath, fn))
        
        return(df)
      }
      
      if (type_supervision == "both") {
        rows <- is.na(df[["sub_location"]]) & !is.na(df_fst[["svm"]])
        df[["sub_location"]][rows] <- df_fst[["svm"]][rows]
        rm(list = c("rows"))
      }
    }
    
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
    id = "gene", col_group = NULL, col_order = NULL, 
    col_cluster = "cluster", levels_group = NULL, col_fraction = NULL, 
    levels_fraction = NULL, col_subtype = NULL, levels_subtype = NULL, 
    data_type = "Trend", pat = "Trend_[ONZ]_.*nclust\\d+\\.txt$", ext = "txt", 
    tie_method = "none", qt = .5, complete_cases = FALSE, panel_ids = NULL, 
    filepath = NULL, label_scheme_sub = NULL, theme = NULL, make_plot = TRUE, 
    dots = NULL)
{
  ## (1) Load input data
  col_x <- col_group
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
  
  # By each panel in a n_clust, e.g., 1:16 in n_clust = 16
  df_med <- lapply(split(df, df[[col_cluster]]), function(dfx) {
    dfx <- dfx |>
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
        df_med = df_med, n_clust = n_clust, col_cluster = col_cluster, 
        levels_group = levels_group, 
        col_x = col_x, col_y = col_y, 
        col_fraction = col_fraction, levels_fraction = levels_fraction, 
        color_fractions = NULL, theme = NULL, 
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
#' @param color_fractions The color scheme of fractions.
#' @param panel_ids The subset of panel_ids for visualization.
#' @param cl_id A cluster ID.
#' @param filepath An output filepath.
#' @param out_nm An output filename.
#' @param theme The theme of plot.
#' @param ... Additional parameters for aesthetics.
make_subcell_plot <- function (
    df_med, n_clust, col_cluster = "cluster", 
    col_x = "Group", col_y = "loc_score", 
    levels_group = NULL, col_fraction = "Fraction", levels_fraction = NULL, 
    color_fractions = NULL, 
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
  
  if (length(levels_fraction)) {
    df_med <- df_med |>
      dplyr::mutate(
        !!col_fraction := 
          factor(!!rlang::sym(col_fraction), levels = levels_fraction)
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
    # scale_fill_viridis_d(option = "plasma", end = 0.85) + 
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
  
  if (length(color_fractions)) {
    p <- p + 
      scale_fill_manual(
        values = color_fractions,
        breaks = c("CP", "NP", "ChA"), 
        drop = FALSE
      )
  } else {
    p <- p +
      scale_fill_viridis_d(option = "plasma", end = 0.85)
  }

  ggsave_dots <- 
    set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
  
  rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(out_nm)),
                    # plot = p, 
                    width = width, 
                    height = height, 
                    !!!ggsave_dots)) |>
    rlang::eval_tidy()
}


#' Helper of \link{assign_compartment_to_trend_cluster}.
#' 
#' @param col_x The column name of X values.
#' @param col_y The column name of Y values for summarization.
#' @param col_cluster The column key of cluster nodes.
#' @param tie_method A method handling ties in localization scores.
#' @param df_med A median description of localization scores by modules.
#' @inheritParams anal_prnTrend
hassign_compartment_to_trend_cluster <- function (
    df_med = NULL, col_subtype = NULL, col_fraction = NULL,
    levels_ref_compartments = NULL, 
    col_group = NULL, col_order = NULL, col_cluster = "cluster", 
    col_x = NULL, col_y = NULL, tie_method = "none")
  
{
  if (length(col_subtype) && col_subtype %in% names(df_med)) {
    res <- lapply(split(df_med, df_med[[col_subtype]]), function (dfx) {
      tryCatch(
        assign_compartment_to_trend_cluster(
          df_med = dfx, 
          tie_method = tie_method, 
          col_fraction = col_fraction, 
          col_subtype = col_subtype,
          col_cluster = col_cluster,
          col_y = col_y, 
          levels_ref_compartments = levels_ref_compartments),
        error = function(e) NULL)
    }) |>
      dplyr::bind_rows(.id = col_subtype) |>
      dplyr::ungroup() |>
      dplyr::arrange(!!rlang::sym(col_fraction), !!rlang::sym(col_y))
  } else {
    res <- tryCatch(
      assign_compartment_to_trend_cluster(
        df_med = df_med, 
        tie_method = tie_method, 
        col_fraction = col_fraction, 
        col_subtype = col_subtype, 
        col_cluster = col_cluster,
        col_y = col_y, 
        levels_ref_compartments = levels_ref_compartments),
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


#' Assign the primary subcellular compartment using median-summarized
#' localization scores.
#' 
#' For one \code{col_subtype}.
#' 
#' No guarantee that each compartment will have a corresponding primary cluster.
#' 
#' @param col_y The column name of Y values for summarization.
#' @param col_cluster The column key of cluster nodes.
#' @param tie_method A method handling ties in localization scores.
#' @param df_med A median description of localization scores by modules.
#' @inheritParams anal_prnTrend
assign_compartment_to_trend_cluster <- function (
    df_med = NULL, tie_method = "none", col_y = "score", col_cluster = "cluster", 
    col_fraction = NULL, col_subtype = NULL, levels_ref_compartments = NULL)
{
  df_med <- df_med |>
    dplyr::mutate(
      !!col_fraction := 
        factor(!!rlang::sym(col_fraction), levels = levels_ref_compartments), )
  
  subcell_levs <- unique(df_med[[col_fraction]])
  subcell_levs <- factor(subcell_levs, levels = levels_ref_compartments)
  
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
    df_pri_min[["is_leading_fraction"]] <- TRUE
    rm(list = c("df_tie", "df_tie_min"))
    
    ## Add back fractions "1.a" by finding the nearest neighbors
    df_sec[["is_leading_fraction"]] <- FALSE
    mis_levs <- subcell_levs[!subcell_levs %in% df_pri_min[[col_fraction]]]
    
    if (length(mis_levs)) {
      # By each panel under a n_clust, e.g., 16
      df_secs <- split(df_sec, df_sec[[col_cluster]])
      
      # By each organelle
      df_no_apex <- lapply(mis_levs, function (l) {
        # The row index of the target organelle 'l' at each panel
        mt_rownums_oks <- sapply(df_secs, function (x) {
          which(x[[col_fraction]] == l)
        })
        
        # The panel that gives the lowest row index for the target organelle 'l'
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
  } else {
    df_pri_min[["subcellular..."]] <- NULL
    df_pri_min[["is_leading_fraction"]] <- TRUE
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
#' @param col_fraction The column key to sample fractions.
#' @param levels_fraction The levels of subcellular fractions.
#' @param col_subtye For uses with subcellular data. The column key in metadata
#'   (\code{expt_smry.xlsx}) linking to sample sub-types, e.g., \code{Control,
#'   Treated}.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param pat The pattern of input filenames.
#' @param ext The extension of input filenames.
#' @param data_type The type of input data.
#' @param tie_method Not used. A method for tie-breaking.
#' @param type_supervision Only for uses with supervised learning: \code{sole},
#'   supervised learning only; \code{both}, both supervised and unsupervised
#'   learning.
#' @param qt The quantile of localization score for thresholding subcellular
#'   compartment specificity. The default is to use median.
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
    dat_dir = NULL, col_group = NULL, levels_group = NULL, 
    col_select = NULL, col_order = NULL, 
    col_fraction = NULL, levels_fraction = NULL, 
    levels_ref_compartments = NULL, 
    col_subtype = NULL, levels_subtype = NULL, tie_method = "none", 
    type_supervision = c("sole", "both"), qt = .5, 
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
  
  id <- tryCatch(
    match_call_arg("normPSM", "group_pep_by"), error = function(e) NA)
  
  if (is.na(id)) {
    id <- tryCatch(
      match_call_arg("makePepDIANN", "group_pep_by"), 
      error = function(e) NA)
  }
  
  if (is.na(id)) {
    id <- "gene"
  }
  
  if (length(id) != 1L || !id %in% c("prot_acc", "gene")) {
    stop("Value of 'id' need to be one of 'prot_acc', 'gene'.")
  }
  
  if (FALSE) {
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
  }

  type_supervision <- match.arg(type_supervision)
  tie_method <- "none"
  
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
              levels_ref_compartments = levels_ref_compartments, 
              col_subtype = col_subtype, 
              levels_subtype = levels_subtype, 
              tie_method = tie_method, 
              n_clust = n_clust, 
              pat = pat, ext = ext, data_type = data_type, 
              type_supervision = type_supervision, qt = qt, 
              panel_ids = panel_ids, make_plot = make_plot, 
              theme = theme, ...)
}


#' SVM classification of subcellular locations.
#'
#' @param species Species in one of human, mouse and rat.
#' @param id Identifier.
#' @param col_y The column key in \code{df_ref} for training labels.
#' @param dat_dir The current working directory.
#' @param filepath The filepath of inputs and outputs.
#' @param n_clust The number of clusters.
#' @param df_ref A data frame of reference markers.
#' @param u_subtypes The unique levels of sample types. A NULL value indicates
#'   only one subtype.
#' @param levels_fraction The levels of subcellular fractions.
#' @param ranges for \code{cost} and \code{gamma} for grid searching with
#'   \code{radial} kernel.
#' @param use_cache Logical; use cache or not.
#' @param seed A seed for SVM.
#' @inheritParams anal_prnTrend col_subtype
svmSubcell <- function (
    species = "human", id = "gene", dat_dir = NULL, filepath = NULL, 
    label_scheme_sub = NULL, df_ref = NULL, 
    n_clust = NULL, col_y = "location", 
    col_group = NULL, levels_group = NULL, 
    col_subtype = NULL, levels_subtype = NULL, u_subtypes = NULL, 
    col_fraction = NULL, 
    levels_fraction = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), 
    ranges = list(cost = c(1, 10, 50, 100), gamma = c(0.01, 0.1, 0.5, 1)), 
    use_cache = FALSE, seed = 42, ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(label_scheme_sub)) label_scheme_sub <- load_ls_group(dat_dir)
  if (is.null(n_clust)) stop("Proive a range for 'n_clust'.")
  
  ## (1) Obtain expression data and reference data
  tempdata <- prepSubcell(
    species = species, id = id, dat_dir = dat_dir, filepath = filepath, 
    label_scheme_sub = label_scheme_sub, df_ref = df_ref, 
    n_clust = n_clust, col_group = col_group, 
    col_subtype = col_subtype, levels_subtype = levels_subtype, 
    u_subtypes = u_subtypes, 
    col_fraction = col_fraction, levels_fraction = levels_fraction)
  df_exprs  <- tempdata[["df_exprs"]]
  df_ref    <- tempdata[["df_ref"]]
  sc_lookup <- tempdata[["sc_lookup"]]

  saveRDS(df_exprs, file.path(filepath, "df_exprs.rds"))
  saveRDS(df_ref, file.path(filepath, "df_ref.rds"))
  rm(list = "tempdata")
  
  if (length(u_subtypes) > 1L && !identical(names(df_exprs), u_subtypes)) {
    stop("Developer: mismatches between expression data and sample types.")
  }

  ## (2) Results by sub-types
  res <- vector("list", length(df_exprs))
  names(res) <- names(df_exprs)
  
  for (i in seq_along(res)) {
    res[[i]] <- train_loc_svm(
      df_expr = df_exprs[[i]], sub_type = u_subtypes[[i]], 
      df_ref = df_ref, id = id, col_y = col_y, col_subtype = col_subtype, 
      col_fraction = col_fraction, filepath = filepath, ranges = ranges, 
      use_cache = use_cache, seed = seed)
  }

  ## 'purity' is SVM probability
  #   gene  Fraction    purity entropy svm   loc_score
  # 1 Aaas  Centrosome  0.155    -76.3 ER        11.9 
  # ...
  # 5 Aaas  ER         45.5      -76.3 ER      3513. 
  
  list(prob = res, expr = df_exprs)
}


#' Train localization SVM.
#'
#' @param df_expr A matrix of expression data, typically normalized intensities
#'   that summed to one across fractions for each protein.
#' @param df_ref A reference data frame.
#' @param id An identifier.
#' @param col_y The column key in \code{df_ref} for training labels.
#' @param col_fraction The identifier of subcellular compartment.
#' @param sub_type A label of data subtype, e.g., "Control".
#' @param ranges for \code{cost} and \code{gamma} for grid searching with
#'   \code{radial} kernel.
#' @param use_cache Logical; use cache or not.
#' @param seed A seed for SVM.
#' @return Two lists of (1) SVM probability matrix and (2) expression data.
train_loc_svm <- function (
    df_expr, sub_type = NULL, df_ref, id = "gene", col_y = "location", 
    col_subtype = NULL, col_fraction = NULL, filepath = NULL, 
    ranges = list(cost = c(1, 10, 50, 100), gamma = c(0.01, 0.1, 0.5, 1)), 
    use_cache = FALSE, seed = 42L, ...)
{
  ## 1. Setup
  Xmat <- df_expr[df_ref[[id]], ]
  Xmat <- Xmat[match(rownames(Xmat), df_ref[[id]]), ]
  
  if (!identical(rownames(Xmat), df_ref[[id]])) {
    stop("Develepor: mismatches in identifiers: ", id, ".")
  }
  
  if (!col_y %in% colnames(df_ref)) {
    stop("The column of Y labels", col_y, " not found in reference data.")
  }
  
  y_factor <- factor(df_ref[[col_y]])
  df_train <- data.frame(Xmat, check.names = FALSE)
  df_train[[col_y]] <- y_factor
  class_names <- levels(y_factor)
  
  ## 2. Training
  if (!file.exists(file_svmp <- file.path(filepath, "tuned_svmp.rds"))) {
    use_cache <- FALSE
  }
  
  if (use_cache) {
    tuned_svmp <- readRDS(file_svmp)
  } else {
    class_counts <- table(df_train[[col_y]])
    weights <- max(class_counts) / class_counts
    
    set.seed(seed)
    tuned_svmp <- e1071::tune(
      e1071::svm,
      train.x = stats::reformulate(".", response = col_y),
      data = df_train,
      kernel = "radial",
      class.weights = weights, 
      ranges = ranges, 
      probability = TRUE, 
      scale = TRUE)
    
    saveRDS(tuned_svmp, file_svmp)
    rm(list = c("class_counts", "weights"))
  }

  svm_model <- tuned_svmp$best.model
  
  # Calculate training set priors (class proportions)
  priors <- table(df_train[[col_y]]) / nrow(df_train)
  
  # Raw model probabilities for all classes
  svm_preds <- predict(svm_model, newdata = df_expr, probability = TRUE)
  raw_probs <- attr(svm_preds, "probabilities")
  
  # Adjusted probability by its corresponding class prior
  prob_mat <- sweep(raw_probs, MARGIN = 2, STATS = priors, FUN = "/")
  prob_mat  <- prob_mat[, class_names]
  
  # Predict the class with the maximum prior-adjusted probability
  svm_final <- colnames(prob_mat)[apply(prob_mat, 1, which.max)]
  svm_final <- factor(svm_final, levels = class_names)
  
  df_pred   <- tibble::tibble(!!id := names(svm_preds), svm = svm_final)

  ## 3. Outputs
  entropies  <- -rowSums(prob_mat * log10(prob_mat))
  scores_loc <- prob_mat * (1 - entropies)
  entropies  <- tibble::tibble(!!id := names(entropies), entropy = entropies)
  entropies  <- entropies |> dplyr::left_join(df_pred, by = id)
  
  scores_loc <- scores_loc |>
    data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = id) |> 
    tidyr::pivot_longer(
      -c(!!rlang::sym(id)), names_to = col_fraction, values_to = "loc_score")
  
  prob_mat <- prob_mat |>
    data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = id) |> 
    tidyr::pivot_longer(
      -c(!!rlang::sym(id)), names_to = col_fraction, values_to = "purity")
  
  ans <- prob_mat |>
    dplyr::left_join(entropies, by = id) |>
    dplyr::left_join(scores_loc, by = c(id, col_fraction))
}


#' Prepare subcellular data.
#' 
#' @param species Species in one of human, mouse and rat.
#' @param id Identifier.
#' @param dat_dir The current working directory.
#' @param filepath The filepath of inputs and outputs.
#' @param df_ref A data frame of reference markers.
#' @param key_col The column key where data to be used.
#' @param levels_fraction The levels of subcellular fractions.
#' @param u_subtypes The unique levels of sample types. A NULL value indicates
#'   only one subtype.
#' @inheritParams anal_prnTrend col_subtype
prepSubcell <- function (
    species = "human", id = "gene", dat_dir = NULL, filepath = NULL, 
    label_scheme_sub = NULL, df_ref = NULL, 
    n_clust = NULL, key_col = "log10Int", col_group = NULL, 
    col_subtype = NULL, levels_subtype = NULL, u_subtypes = NULL, 
    col_fraction = NULL, 
    levels_fraction = 
      c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"), ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(label_scheme_sub)) label_scheme_sub <- load_ls_group(dat_dir)
  if (is.null(n_clust)) stop("Provide a range of 'n_clust'.")
  
  abbr_species <- switch(
    species, 
    human = "hs",
    mouse = "mm",
    rat = "rn",
    stop("species: ", species, " is not in one of `human`, `mouse` or `rat`.")
  )

  if (is.null(df_ref)) {
    df_ref <- local({
      file_db <- paste0("subcell_", abbr_species)
      data(package = "proteoQ", list = file_db, envir = environment())
      df_ref <- get(file_db)
      df_ref |> dplyr::mutate(location = sub("ERGIC", "ER", location))
    })
  }

  if (length(levels_fraction) && 
      all(levels_fraction %in% c("NP", "CP", "ChA"))) {
    df_ref[["location"]] <- df_ref[["location_slim3"]]
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
      tidyr::pivot_longer(-!!id, names_to = col_group, values_to = key_col) |>
      dplyr::mutate(!!key_col := log10(!!rlang::sym(key_col)))
  } else {
    df <- readr::read_tsv(file.path(filepath, filenames_df2[[1]]))
  }

  if (no_subtype <- !length(col_subtype)) {
    label_scheme_sub[[col_subtype <- "sample_type_unknown"]] <- "unknown"
  }
  
  sc_lookup <- label_scheme_sub |>
    dplyr::select(dplyr::one_of(col_group, col_fraction, col_subtype)) |>
    unique() # |> dplyr::rename(Sample_ID := col_group)

  rows <- match(df[[col_group]], sc_lookup[[col_group]])
  
  if (!col_subtype %in% names(df)) {
    df <- dplyr::bind_cols(df, sc_lookup[rows, col_subtype])
    u_subtypes <- unique(df[[col_subtype]]) # Value: "" if 'no_subtype' 
  }
  
  if (!length(levels_subtype)) {
    levels_subtype <- unique(sc_lookup[[col_subtype]])
  }
  
  # Should not incur
  # if (!length(col_fraction)) col_fraction <- "fraction..."

  if (!col_fraction %in% names(df)) {
    df <- dplyr::bind_cols(df, sc_lookup[rows, col_fraction])
  }

  rm(list = c("rows"))
  
  dfs <- split(df, df[[col_subtype]])
  
  # Match the order of levels_subtype
  if (length(levels_subtype) > 1L) {
    dfs <- dfs[levels_subtype]
  }
  
  # Normalized intensity (percentages)
  df_exprs <- vector("list", length(dfs))
  names(df_exprs) <- names(dfs)
  
  for (i in seq_along(dfs)) {
    df_expr <- make_subcell_expr(
      df = dfs[[i]], id = id, key_col = key_col, key_group = col_group, 
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      scale = TRUE)
    
    # Change from values under col_group to values under col_fraction
    local({
      sc_fracs <- sc_lookup[[col_fraction]]

      if (any(bads <- !levels_fraction %in% sc_fracs)) {
        stop("User-provided 'levels_fraction' not found in metadata:\n", 
             levels_fraction[bads])
      }
    })
    
    colnames(df_expr) <- sc_lookup[[col_fraction]][
      match(colnames(df_expr), sc_lookup[[col_group]])
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
    sc_lookup = sc_lookup, 
    abbr_species = abbr_species, 
    filenames_df2 = filenames_df2
  )
}


#' Make UMAP for 10-fraction data.
#' 
#' More generally where sample fractions do not match biological compartments.
#' 
#' @param df_expr Expression data.
#' @param sub_type The subtype of samples.
#' @param file_umap The file name for saving UMAP results.
#' @param dat_dir The current working directory.
#' @param filepath The file path to \code{file}.
#' @param id An identifier.
#' @param col_subtype The column key to \code{sub_type}.
#' @param scale Logical; scale data or not.
#' @param seed A seed.
smake_f10_umap <- function (
    df_expr, sub_type, file_umap = NULL, dat_dir = NULL, filepath = NULL, 
    id = "gene", col_subtype = "Type", scale = TRUE, seed = 42)
{
  if (is.null(file_umap)) {
    file_umap  <- file.path(filepath, paste0("umap_", sub_type, ".tsv"))
  }
  
  ncol <- ncol(df_expr)
  rows <- rowSums(is.na(df_expr)) == ncol | 
    rowSums(sapply(df_expr, is.infinite)) == ncol
  df_expr_b <- df_expr[rows, ]
  df_expr   <- df_expr[!rows, ]
  
  if (scale) {
    df_expr <- df_expr / rowSums(df_expr, na.rm = TRUE)
  }
  
  set.seed(seed)
  res_umap <- umap::umap(df_expr)
  df_umap  <- res_umap$layout[, 1:2]
  df_umap  <- round(df_umap, digits = 2L)
  colnames(df_umap) <- paste0("UMAP", 1:2)
  rownames(df_umap) <- rownames(df_expr)
  
  df_umap <- tibble::rownames_to_column(data.frame(df_umap), id)
  df_umap[[col_subtype]] <- sub_type
  
  readr::write_tsv(df_umap, file_umap)
  
  df_umap
}


#' Helper of Preparing the matrix of expression data.
#' 
#' @inheritParams make_subcell_expr
hmake_subcell_expr <- function (
    df = NULL, label_scheme_sub = NULL, id = "gene", key_col = "log10Int", 
    scale = TRUE, key_group = NULL, 
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
  names(df_exprs) <- names(dfs)
  
  for (i in seq_along(dfs)) {
    df_exprs[[i]] <- make_subcell_expr(
      df = dfs[[i]], 
      id = id, key_col = key_col, 
      key_group = key_group, 
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      scale = scale)
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
    key_group = "Group", col_fraction = NULL, 
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
  
  # Normalized intensity; do not scale if log2Ratios
  if (scale) {
    df_expr <- df_expr / rowSums(df_expr, na.rm = TRUE)
  }
  
  df_expr
}


#' Cluster SVM probabilities for a given sample sub-type using an unsupervised
#' method.
#'
#' @param df_prob A data frame of SVM classification (column \code{svm}),
#'   probabilities (column \code{purity}) and localization scores (column
#'   \code{loc_score}).
#' @param col_score The column key of localization scores.
#' @param dat_dir The working directory.
#' @param filepath A filepath for inputs and outputs.
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param n_clust The number of clusters.
#' @param choice A clustering method.
#' @param use_cache Logical; use cache or not.
cluster_svmprobs <- function (
    df_prob = NULL, sub_type = NULL, col_score = "loc_score", 
    col_fraction = NULL, 
    dat_dir = NULL, filepath = NULL, id = "gene", n_clust = NULL, 
    choice = c("cmeans", "clara", "kmeans", "pam", "fanny"), 
    use_cache = FALSE, ...) 
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("'n_clust' cannot be NULL.")
  
  choice <- match.arg(choice)
  dots <- rlang::enexprs(...)
  
  file_trend_res <- file.path(filepath, "trendRes.rds")
  if (use_cache && !file.exists(file_trend_res)) {
    use_cache <- FALSE
  }
  
  if (use_cache) {
    return(readRDS(file_trend_res))
  }
  
  # Matrix of compartment purity based on SVM probabilities
  Mat <- df_prob[, c(id, col_fraction, col_score)] |>
    tidyr::pivot_wider(
      id_cols = id, 
      names_from = col_fraction, 
      values_from = col_score,
      # values_fill = 0, 
    )
  
  ids <- Mat[[id]]
  Mat <- as.matrix(Mat[, colnames(Mat) != id])
  rownames(Mat) <- ids

  ans_dots <- find_trend_m(
    dots, choice = choice, n_clust = n_clust, df_mean_log2r = Mat)
  dots <- ans_dots[["dots"]]
  n_clust <- ans_dots[["n_clust"]]
  # rm(list = c("ans_dots"))
  
  ##  Unsupervised classification of SVM probabilities at each 'n_clust'
  out <- lapply(
    n_clust, makeTrendRes, 
    df_mean_log2r = Mat, choice = choice, dots = dots, id = id)

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
  
  saveRDS(out, file_trend_res)
  
  out
}


#' Join the group information for subcellular data.
#' 
#' @param ans_cl Clustered trend results.
#' @param filepath A file path.
#' @param label_scheme_sub Metadata.
#' @param col_group The column key to group information.
#' @param levels_group The levels of data group for ordering on the X-axis in
#'   visualization.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param col_fraction The column key to subcellular fractions.
#' @param levels_fraction The levels of subcellular compartments.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param n_clust The number of clusters.
#' @param no_subtype Logical; data contain more than one sub-type or not.
#' @param is_pseudo_fraction Logical; data contain pseudo fractions or not.
join_subcell_group <- function (
    ans_cl, filepath, label_scheme_sub, col_group, levels_group, 
    col_fraction, col_subtype, n_clust = 16:32, 
    no_subtype = FALSE, is_pseudo_fraction = TRUE)
{
  sc_lookup <- label_scheme_sub |>
    dplyr::select(dplyr::one_of(col_group, col_fraction, col_subtype)) |>
    unique()
  
  out_nms <- paste0("Protein_svmProb_Z_nclust", n_clust, ".tsv")
  
  ans <- vector("list", length(n_clust))
  for (i in seq_along(n_clust)) {
    # Combine results from the i-th n_clust for each sub-type
    tempdata <- lapply(ans_cl, `[[`, i) |>
      dplyr::bind_rows(.id = if (no_subtype) NULL else col_subtype)
    
    # Join 'col_group' information:
    # However, 'col_fraction' from ans[[i]] may be pseudo fractions from SVM, 
    # such as compartment Centrosome etc., not sample fractions CP, NP..., 
    # if experimental fractions do not coincide with compartments
    # In that case values under 'col_group' would be all NA
    if (is_pseudo_fraction) {
      tempdata[[col_group]] <- paste0(
        tempdata[[col_subtype]], "_", tempdata[[col_fraction]])
    } else {
      tempdata <- tempdata |>
        dplyr::left_join(sc_lookup, by = c(col_fraction, col_subtype))
    }
    
    tempdata <- tempdata |>
      reloc_col_before_first(col_group) |>
      reloc_col_after(col_fraction, col_subtype)
    
    u_data_grps <- unique(tempdata[[col_group]])
    missing_grp_levs <- u_data_grps[!u_data_grps %in% levels_group]
    
    if (length(missing_grp_levs)) {
      warning(
        "Mismatched 'levels_group' between data and user's specification:\n", 
        paste(missing_grp_levs, collapse = ", "), "\n", 
        "Ignore data arrangement by 'levels_group'.")
    } else {
      tempdata <- tempdata |>
        dplyr::mutate(
          !!col_group := factor(!!rlang::sym(col_group), levels = levels_group))
    }
    
    readr::write_tsv(tempdata, file.path(filepath, out_nms[[i]]))
    
    ans[[i]] <- tempdata
  }
  
  names(ans) <- n_clust
  
  ans <- ans |> dplyr::bind_rows(.id = "n_clust")
  
  invisible(ans)
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
#' @param df_ref A data frame of reference markers.
#' @param choice A clustering method.
#' @param min_prot_n_peps The minimum number of detected peptides under a
#'   protein for consideration.
#' @param make_plot Logical; an option to plot clusters of probabilities.
#' @param levels_group The levels of data group for ordering on the X-axis in
#'   visualization.
#' @param levels_subtype The levels of sample subtypes, e.g., Control, Treated.
#' @param levels_ref_compartments Levels of reference compartments.
#' @param col_fraction The column key to subcellular fractions.
#' @param levels_fraction The levels of subcellular compartments.
#' @param type_supervision \code{sole}, supervised learning only; \code{both},
#'   both supervised and unsupervised learning.
#' @param qt The quantile of localization score for thresholding
#'   subcellular compartment specificity. The default is to use median.
#' @param use_cache Logical; Should be FALSE. if TRUE, redo trend
#'   analysis.
#' @param ranges for \code{cost} and \code{gamma} for grid searching with
#'   \code{radial} kernel.
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
#' levels_ref_compartments <- levels_fraction
#'
#' # c("Veh_CP", "TG_CP", "Veh_NP", "TG_NP", "Veh_ChA", "TG_ChA")
#' levels_group <- load_ls_group(dat_dir = dat_dir) |>
#'   dplyr::arrange_at("Order") |>
#'   dplyr::pull("Group") |>
#'   unique()
#'
#' # 'levels_group' can be either users' provided or deduced automatically.
#' classify_subcellular_supervised(
#'   dat_dir = dat_dir, id = "gene", species = "human", n_clust = 16:36,
#'   col_group = "Group", levels_group = levels_group,
#'   col_fraction = "Fraction", levels_fraction = levels_fraction,
#'   col_subtype = "Type", levels_subtype = levels_subtype,
#'   levels_ref_compartments = levels_ref_compartments)
#'
#' ## (2) Sample fractions are different to organelle compartments
#' data("subcell_mm")
#' compartments <- names(sort(table(subcell_mm$location), decreasing = TRUE))
#' levels_ref_compartments <- compartments
#' levels_group <- compartments # Assume one 'subtype'
#'
#' classify_subcellular_supervised(
#'   dat_dir = dat_dir, id = "gene", species = "mouse", n_clust = 16:36,
#'   col_group = "Group", # for group summarization of raw expression data
#'   levels_group = levels_group, # for the order on the X-axis in visualization
#'   col_fraction = "Fraction", # for preparation of summarized expression data
#'   levels_fraction = # to replace the colnames of summarized expression data
#'     c("1K", "3K", "5K", "9K", "12K", "15K", "30K", "79K", "120K", "SN"),
#'   levels_ref_compartments = levels_ref_compartments)
#' }
#' @export
classify_subcellular_supervised <- function (
    dat_dir = NULL, filepath = NULL, label_scheme = NULL, 
    id = "gene", species = "human", n_clust = 16:36, df_ref = NULL, 
    col_group = NULL, levels_group = NULL, 
    col_fraction = NULL, levels_fraction = NULL, 
    levels_ref_compartments = NULL, 
    col_subtype = NULL, levels_subtype = NULL, 
    col_select = NULL, col_order = NULL, 
    type_supervision = c("both", "sole"), qt = .75, min_prot_n_peps = 2L, 
    ranges = list(cost = c(1, 10, 50, 100), gamma = c(0.01, 0.1, 0.5, 1)), 
    # use_cache = FALSE, 
    make_plot = FALSE, seed = 42L, 
    ...) 
{
  dots <- rlang::enexprs(...)
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("Specifify a range of 'n_clust'.")
  if (is.null(label_scheme)) label_scheme <- load_ls_group(dat_dir)
  
  if (!length(levels_fraction)) {
    stop("Need specification of levels_fraction")
  }
  
  if (!length(levels_ref_compartments)) {
    stop("Need specification of levels_ref_compartments.")
  }
  
  if (qt < .5) {
    warning("The quantile value is ", qt, "\nConsider a value >= 0.5.")
  }
  
  type_supervision <- match.arg(type_supervision)
  is_pseudo_fraction <- !setequal(levels_fraction, levels_ref_compartments)

  check_subcell(
    levels_fraction = levels_fraction, levels_subtype = levels_subtype, 
    col_subtype = col_subtype)
  
  label_scheme_sub <- if (length(col_select)) {
    label_scheme |> dplyr::filter(!is.na(.data[[col_select]]))
  } else {
    label_scheme
  }
  
  ### (1) Obtain aggregated expression data if no prior 'anal_prnTrend' execution
  use_cache <- dots[["use_cache"]]
  if (is.null(use_cache)) use_cache <- FALSE
  dots[["use_cache"]] <- NULL
  
  file_log10int <- "Protein_Intensity_summary.tsv"
  if (use_cache && !file.exists(file.path(filepath, file_log10int))) {
    use_cache <- FALSE
  }
  
  if (!use_cache) {
    anal_prnTrend(
      col_group = col_group, 
      col_order = col_order,
      col_fraction = col_fraction, 
      col_subtype  = col_subtype, 
      cluster_data = FALSE, # to skip clustering steps
      n_clust = n_clust, 
      filter_by_npep = exprs(prot_n_pep >= !!min_prot_n_peps), )
  }
  
  ### (2) Obtain SVM classification and probability matrix (supervised)
  n_subtype  <- length(col_subtype)
  no_subtype <- !n_subtype
  u_subtypes <- if (no_subtype) NULL else unique(label_scheme_sub[[col_subtype]])
  
  ans_svm <- svmSubcell(
    species = species, id = id, dat_dir = dat_dir, filepath = filepath, 
    df_ref = df_ref, label_scheme_sub = label_scheme_sub, 
    col_group = col_group, levels_group = levels_group, 
    col_subtype = col_subtype, levels_subtype = levels_subtype, 
    col_fraction = col_fraction, levels_fraction = levels_fraction, 
    u_subtypes = u_subtypes, n_clust = n_clust, 
    ranges = ranges, use_cache = use_cache, seed = seed)
  
  svm_probs <- ans_svm[["prob"]] # Table: SVM classifications and Probabilities
  # df_exprs  <- ans_svm[["expr"]] # Raw expression data by fractions

  ### (3) Classify SVM probabilities (unsupervised)
  # ans_cl: list by 'sub_types'
  #  list of data frames by 'n_clust'
  # Not `mapply` since 'u_subtypes' can be NULL
  ans_cl <- vector("list", length(svm_probs))
  names(ans_cl) <- u_subtypes
  
  # Cluster SVM probabilities by each sub-type
  #  ans_cl[[i]]: A list contains sub-lists by each 'n_clust'
  for (i in seq_along(svm_probs)) {
    ans_cl[[i]] <- cluster_svmprobs(
      df_prob = svm_probs[[i]], sub_type = u_subtypes[[i]], 
      col_fraction = col_fraction, dat_dir = dat_dir, 
      filepath = filepath, id = id, n_clust = n_clust, use_cache = use_cache, 
      iter.max = 250)
  }
  
  ## (4) Combine data under the same n_clust across all sub_types
  sc_lookup <- label_scheme_sub |>
    dplyr::select(dplyr::one_of(col_group, col_fraction, col_subtype)) |>
    unique()

  join_subcell_group(
    ans_cl = ans_cl, filepath = filepath, label_scheme_sub = label_scheme_sub, 
    col_group = col_group, levels_group = levels_group, 
    col_fraction = col_fraction, col_subtype = col_subtype, n_clust = n_clust, 
    no_subtype = no_subtype, is_pseudo_fraction = is_pseudo_fraction)

  ### (5) Classify subcellular compartments by medium statistics
  plot_prnSubcellular(
    col_order = "colnames", # X-axix ordered by compartment
    levels_group = levels_group, 
    col_subtype = col_subtype, 
    levels_subtype = levels_subtype,
    col_fraction = col_fraction,
    levels_fraction = levels_fraction,
    levels_ref_compartments = levels_ref_compartments,
    type_supervision = type_supervision, 
    qt = qt, 
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
#' @param label_scheme The metadata generated from \link{load_expts}.
#' @param id Character string; one of \code{pep_seq}, \code{pep_seq_mod},
#'   \code{prot_acc} and \code{gene}.
#' @param species A species in one of \code{"human", "mouse", "rat"} for
#'   identifying the database of reference markers.
#' @param n_clust The number of clusters.
#' @param qt The quantile of localization score for thresholding
#'   subcellular compartment specificity. The default is to use median.
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
    dat_dir = NULL, filepath = NULL, label_scheme = NULL, 
    id = "gene", species = "human", n_clust = 16:36, qt = .5, 
    col_group = NULL, levels_group = NULL, 
    col_fraction = NULL, levels_fraction = NULL, 
    col_subtype = NULL, levels_subtype = NULL, 
    col_select = NULL, col_order = NULL, 
    min_prot_n_peps = 2L, make_plot = TRUE, seed = 42L, ...) 
{
  options(warn = 1L)
  
  dots <- rlang::enexprs(...)
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(n_clust)) stop("Specify a range of 'n_clust'.")
  if (is.null(label_scheme)) label_scheme <- load_ls_group(dat_dir)

  check_subcell(
    levels_fraction = levels_fraction, levels_subtype = levels_subtype, 
    col_subtype = col_subtype)
  
  label_scheme_sub <- if (length(col_select)) {
    label_scheme |> dplyr::filter(!is.na(.data[[col_select]]))
  } else {
    label_scheme
  }

  ### (1) Obtain aggregated expression data  
  # min_prot_n_peps <- eval(min_prot_n_peps)
  anal_prnTrend(
    col_select = col_select, 
    col_group  = col_group, 
    col_order  = col_order,
    col_fraction = col_fraction, 
    col_subtype  = col_subtype, 
    n_clust = n_clust, 
    iter.max = 250, 
    filter_by_npep = exprs(prot_n_pep >= !!min_prot_n_peps), )
  
  ### (2) Optional visualization of trend clusters
  if (make_plot) {
    try(
      plot_prnTrend(
        col_select = col_select,
        col_group = col_group, 
        levels_group = levels_group, 
        col_order = col_order,
        ymax = 10, 
        ymin = -10,
        limitsize = FALSE,
        ncol = 4, 
      )
    )
  }
  
  ### (3) Classify subcellular
  plot_prnSubcellular(
    col_select = col_select, 
    col_order = col_order, 
    col_group = col_group, 
    levels_group = levels_group, 
    col_subtype = col_subtype, 
    levels_subtype = levels_subtype,
    col_fraction = col_fraction,
    levels_fraction = levels_fraction,
    levels_ref_compartments = levels_fraction,
    qt = qt, 
    ncol = 6,
    pat = "Trend_[ONZ]_.*nclust\\d+\\.txt$", 
    ext = "txt", 
    data_type = "Trend", 
    limitsize = FALSE,
    make_plot = make_plot, 
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
#' @param col_fraction The column key in metadata to subcellular fractions.
#' @param levels_fraction The levels of subcellular compartments.
#' @param col_group The column key in metadata linked to the information of
#'   sample groups.
#' @param levels_group The levels of data group.
#' @param col_subtype Column key of sample types.
#' @param levels_subtype The levels of sample subtypes or treatment conditions.
#' @param df_umap A matrix UMAP coordinates. Automatically determined at
#'   \code{NULL}.
#' @param highlight_markers Logical; highlight marker proteins or not.
#' @param levels_ref_compartments For uses at \code{highlight_markers = TRUE}.
#'   The factor levels of subcellular compartments for uses in marker
#'   visualization.
#'
#' @param label_scheme The metadata generated from \link{load_expts}.
#' @param id An identifier.
#' @param species The species.
#' @param df_markers A data frame of reference markers. Automatically determined
#'   at \code{NULL}.
#' @param n_clust Not used. The number of clusters.
#' @param key_col The key column where the values being used for UMAP.
#' @param col_pred The prediction column.
#' @param type_fraction The type of sub-cellular fraction. At the value of
#'   \code{matched}, the experimental fractions matches biological compartments,
#'   e.g., CP, NP and ChA. At the value of \code{pseudo}, the experimental
#'   fractions (e.g., F1, F2 etc.) are different from biological compartments
#'   (e.g. Centrosome, ER, etc.).
#'
#'   The information will be used to decide which secondary file (df2) to used
#'   and to join the encoding of sub-cellular locations therein for color
#'   aesthetics.
#'
#'   It is possible that users may exclude some levels from UMAP visualization
#'   and thus a simple set difference assessment may not be adequate. For now,
#'   the developer leave this argument for users own discretion.
#'
#' @param df2 Not used. The names of secondary input files.
#' @param filepath A file path of inputs and outputs. Automatically determined
#'   at \code{NULL}.
#' @param filename An output filename.
#' @param ext An extension of output filename, e.g., pdf or png.
#' @param polygons The vertices of polygons.
#' @param exclude_inside Logical; if TRUE, exclude points inside the polygon
#'   defined by \code{poly_x} and \code{poly_y}. Otherwise, exclude points
#'   outside of the polygon.
#' @param wrap_facets Wrap facets or not.
#' @param width_combined The width of a combined plot.
#' @param height_combined The height of a combined plot.
#' @param width_sep The width of a plot separated by compartments.
#' @param height_sep The height of a plot separated by compartments.
#' @param overwrite_umap_res Logical; rewrite UMAP result files or not.
#' @param use_cache Logical; use cached results or not.
#' @param seed A seed for umap.
#'
#' @examples
#' \dontrun{
#' plot_prnSubcell_UMAP(
#'   dat_dir = dat_dir, id = "gene", species = "human",
#'   df2 = NULL, filepath = NULL, filename = NULL, ext = "pdf",
#'   wrap_facets = TRUE, seed = 123, nrow = 2, ncol = 7, width = 16, height = 6)
#' }
#' @export
plot_prnSubcell_UMAP <- function (
    dat_dir = NULL, label_scheme = NULL, df_umap = NULL, id = "gene", 
    species = "human", df_markers = NULL, highlight_markers = FALSE, 
    n_clust = 16:36, 
    key_col = c("log2FC", "log10Int"), 
    type_fraction = c("matched", "pseudo"), 
    levels_ref_compartments = NULL, 
    col_pred = "sub_location", levels_pred = NULL, 
    col_fraction = NULL, levels_fraction =  NULL, 
    col_subtype = NULL, levels_subtype = NULL, 
    col_group = NULL, levels_group = NULL, 
    df2 = NULL, filepath = NULL, filename = NULL, ext = "png", 
    polygons = NULL, exclude_inside = TRUE, 
    overwrite_umap_res = TRUE, 
    wrap_facets = TRUE, break_width = 5, hide_grids = TRUE, drop = FALSE, 
    width_combined = 5.5, height_combined = 4, 
    width_sep = 14, height_sep = 6, 
    use_cache = TRUE, seed = 42, ...) 
{
  dots <- rlang::enexprs(...)

  if (FALSE) {
    use_cache <- dots[["use_cache"]]
    if (is.null(use_cache)) use_cache <- FALSE
    dots[["use_cache"]] <- NULL
  }

  # For better positioning of marker proteins
  dot_write_tetherome_res <- dots[["write_tetherome_res"]]
  dot_drop <- dots[["drop"]]
  dot_wrap_facet <- dots[["wrap_facet"]]
  dots[["write_tetherome_res"]] <- dots[["drop"]] <- dots[["wrap_facet"]] <- NULL
  
  if (!is.null(dot_write_tetherome_res) && 
      is.logical(is.null(dot_write_tetherome_res)) && 
      !dot_write_tetherome_res) {
    write_tetherome_res <- FALSE
  } else {
    write_tetherome_res <- TRUE
  }
  
  if (!is.null(dot_drop) && is.logical(dot_drop) && dot_drop) {
    drop <- TRUE
  } else {
    drop <- FALSE
  }
  
  if (!is.null(dot_wrap_facet) && is.logical(dot_wrap_facet) && !dot_wrap_facet) {
    wrap_facet <- FALSE
  } else {
    wrap_facet <- TRUE
  }
  
  ### 0. Metadata etc.
  if (is.null(dat_dir)) dat_dir  <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  if (is.null(col_group)) col_group <- "Group"
  if (is.null(col_subtype)) col_subtype <- "Type"
  if (is.null(col_fraction)) col_fraction <- "Fraction"
  
  message("The output path: ", filepath)
  
  type_fraction <- match.arg(type_fraction)
  key_col <- match.arg(key_col)
  col_x   <- col_group
  
  if (key_col == "log2FC") {
    file_expr <- "Protein_log2FC_summary.tsv" 
    scale <- FALSE
  } else if (key_col == "log10Int") {
    file_expr <- "Protein_Intensity_summary.tsv"
    scale <- TRUE
  } else {
    stop("Unknown key_col = ", key_col)
  }
  
  if (is.null(label_scheme)) label_scheme <- load_ls_group(dat_dir = dat_dir)

  label_scheme_sub <- label_scheme |>
    dplyr::filter(.data[[col_subtype]] %in% levels_subtype)
  
  ### 1. Marker proteins
  df_markers <- load_ref_markers(
    df_markers = df_markers, 
    species = species, 
    levels_ref_compartments = levels_ref_compartments)
  
  # Insufficient since some levels_ref_compartments
  #  may be excluded by users for UMAP visualization 
  # is_pseudo_fraction <- !setequal(levels_fraction, levels_ref_compartments)

  ### 2. Find input `df2` of classication probabilities
  files_2 <- if (type_fraction == "pseudo") {
    find_trend_df2(
      df2 = df2, n_clust = NULL, 
      pat = "Protein_svmProb_Z_nclust\\d+[_]{0,1}.*\\.tsv", 
      ext = "tsv", scale_log2r = TRUE, impute_na = FALSE, filepath = filepath)
  } else if (type_fraction == "matched") {
    find_trend_df2(
      df2 = df2, n_clust = NULL, 
      scale_log2r = TRUE, impute_na = FALSE, filepath = filepath)
  } else {
    stop("Unknown 'type_fraction'.")
  }
  
  if (!length(files_2)) {
    warning("No classification results found.")
    return(NULL)
  }

  file_2 <- files_2[[1]]
  rm(list = c("files_2"))
  
  ### 3. Load classification results (column sub_location for color aesthetics)
  # Don't filter out entries without sub_location assignment; 
  # Otherwise, it can cause too many missing values for UMAP
  df <- readr::read_tsv(file.path(filepath, file_2)) |>
    # dplyr::filter(!is.na(.data[["sub_location"]])) |>
    dplyr::left_join(df_markers, by = id)

  cnms <- colnames(df)
  if (col_subtype %in% cnms && length(levels_subtype)) {
    df <- df |> 
      dplyr::filter(.data[[col_subtype]] %in% levels_subtype) |>
      dplyr::mutate(
        !!col_subtype := 
          factor(!!rlang::sym(col_subtype), levels = levels_subtype))
  }
  
  # 'col_fraction' encodes pseudo fraction at type_fraction = "pseudo"
  if (length(col_fraction) && length(levels_fraction)) {
    df <- df |>
      dplyr::filter(.data[[col_fraction]] %in% levels_fraction) |>
      dplyr::mutate(
        !!col_fraction := 
          factor(!!rlang::sym(col_fraction), levels = levels_fraction))
  }
  
  if (col_group %in% cnms && length(levels_group)) {
    df <- df |> 
      dplyr::filter(.data[[col_group]] %in% levels_group) |>
      dplyr::mutate(
        !!col_group := 
          factor(!!rlang::sym(col_group), levels = levels_group))
  }
  
  if (length(col_pred) && length(levels_pred)) {
    df <- df |>
      dplyr::filter(.data[[col_pred]] %in% levels_pred) |>
      dplyr::mutate(
        !!col_pred := 
          factor(!!rlang::sym(col_pred), levels = levels_pred))
  }

  if (is.null(filename)) {
    fn_prefix <- if (type_fraction == "matched") {
      paste0(tools::file_path_sans_ext(file_2), "_umap")
    } else if (type_fraction == "pseudo") {
      paste0(tools::file_path_sans_ext(file_2), "_sumap")
    } else {
      stop("Unknown 'type_fraction'.")
    }

    fn_suffix <- ext
    filename  <- paste0(fn_prefix, ".", fn_suffix)
  } else {
    fn_prefix <- tools::file_path_sans_ext(filename)
    fn_suffix <- tools::file_ext(filename)
  }
  
  files_umap  <- file.path(filepath, paste0("umap_", levels_subtype, ".tsv"))
  ok_files_umap <- all(file.exists(files_umap))

  if (use_cache && !all(file.exists(files_umap))) {
    use_cache <- FALSE
  }
  
  # Better to 'bind_row' data since can always subset data by specifying 'levels_subtype'

  ### Obtain UMAP coordinates
  if (use_cache) {
    df_umap <- lapply(files_umap, readr::read_tsv) |> 
      dplyr::bind_rows() |>
      tibble::tibble()
  } else {
    # Load expression data from disk
    if (file.exists(file.path(filepath, file_expr))) {
      dfs_expr <- sep_data_by_subtype(
        file = file_expr, filepath = filepath, id = id, 
        label_scheme_sub = label_scheme_sub, 
        col_group = col_group, col_fraction = col_fraction, 
        col_subtype = col_subtype, levels_subtype = levels_subtype, 
        scale = scale)
    } else {
      # Obtain expression data directly
      dfs_expr <- hmake_subcell_expr(
        df = df, id = id, key_col = key_col, key_group = col_group, 
        label_scheme_sub = label_scheme_sub, 
        col_subtype = col_subtype, levels_subtype = levels_subtype, 
        col_fraction = col_fraction, levels_fraction = levels_ref_compartments, 
        scale = scale)
    }
    
    if (!identical(names(dfs_expr), levels_subtype)) {
      stop("Developer: check for mismatches.")
    }
    
    print(
      paste("The levels of 'levels_ref_compartments': ", 
            paste(levels_ref_compartments, collapse = ", "))
    )
    
    if (is.null(df_umap)) {
      df_umap <- mapply(
        smake_f10_umap, 
        dfs_expr, levels_subtype, files_umap, 
        MoreArgs = list(
          dat_dir = dat_dir, filepath = filepath, 
          id = id, col_subtype = col_subtype, scale = scale, seed = seed
        ), 
        SIMPLIFY = FALSE
      ) |>
        dplyr::bind_rows()
    } else {
      df_umap <- df_umap
    }
  }
  
  ### Combine df and df_umap
  if (length(levels_subtype)) {
    df_umap <- df_umap |>
      dplyr::filter(.data[[col_subtype]] %in% levels_subtype)
  }
  
  col_pred <- "sub_location"
  col_join <- c(id, col_pred)
  key_join <- id

  if (length(col_subtype)) {
    col_join <- c(col_join, col_subtype)
    key_join <- c(key_join, col_subtype)
  }
  
  # Add 'sub_location' column for color aesthetics
  df <- df_umap |>
    dplyr::left_join(df, by = key_join)
  
  # Need re-factoring since col_subtype is part of the the 'key_join'
  df <- df |> 
    dplyr::mutate(
      !!col_subtype := 
        factor(!!rlang::sym(col_subtype), levels = levels_subtype))
  
  # 'col_fraction' encodes pseudo fractions at 'type_fraction' of 'pseudo'
  df <- df |>
    dplyr::arrange_at(c(col_group, col_fraction, col_subtype))

  df <- df |>
    dplyr::mutate(
      !!col_pred := as.character(!!rlang::sym(col_pred)), 
      !!col_pred := dplyr::coalesce(.data[[col_pred]], "Unknown"), 
      !!col_pred := factor(!!rlang::sym(col_pred), levels = levels_pred)
    ) |>
    # Can have NA with some reference levels removed by users or 
    # 'id' in df_umap but not in df
    dplyr::filter(!is.na(.data[[col_pred]])) |> 
    unique()
  
  n_subtypes <- length(levels_subtype)
  draw_polygon <- lapply(polygons, function (x) !is.null(x))
  
  if (!length(draw_polygon)) {
    draw_polygon <- rep_len(FALSE, n_subtypes)
    names(draw_polygon) <- levels_subtype
    draw_polygon <- draw_polygon[levels_subtype]
  }
  
  find_tethers(
    df = df, id = id, dat_dir = dat_dir, filepath = filepath, 
    col_fraction = col_fraction, levels_fraction = levels_fraction, 
    col_subtype = col_subtype, levels_subtype = levels_subtype, 
    levels_ref_compartments = levels_ref_compartments, 
    polygons = polygons, exclude_inside = exclude_inside, 
    write_tetherome_res = write_tetherome_res)
  
  # Make plots
  df_plot <- plot_umap(
    df = df, polygons = polygons, draw_polygon = draw_polygon, 
    highlight_markers = highlight_markers, 
    id = id, 
    col_pred = col_pred, levels_pred =  levels_pred, 
    col_subtype = col_subtype, levels_subtype = levels_subtype, 
    levels_ref_compartments = levels_ref_compartments, 
    col_fraction = col_fraction, levels_fraction = levels_fraction, 
    exclude_inside = exclude_inside, 
    write_tetherome_res = write_tetherome_res,
    filepath = filepath, fn_prefix = fn_prefix, fn_suffix = fn_suffix, 
    break_width = break_width, hide_grids = hide_grids, 
    wrap_facet = wrap_facet, drop = drop, ext = ext, 
    width_combined = width_combined, height_combined = height_combined, 
    width_sep = width_sep, height_sep = height_sep, 
    !!!dots)
  
  invisible(df_plot)
}


#' Helper of finding tetherome
#' 
#' @inheritParams find_tethers_by_subtype
find_tethers <- function (
    df, id = "gene", dat_dir = NULL, filepath = NULL, 
    col_fraction = "Fraction", levels_fraction = NULL, 
    col_subtype = "Type", levels_subtype = NULL, levels_ref_compartments = NULL, 
    polygons = NULL, exclude_inside = TRUE, write_tetherome_res = TRUE, 
    out_name = "marker_tetherome.tsv")
{
  dfs <- split(df, df[[col_subtype]])
  dfs <- dfs[levels_subtype]
  polygons <- polygons[levels_subtype]
  
  df_tether <- vector("list", length(dfs))
  
  for (i in seq_along(dfs)) {
    df_tether[[i]] <- find_tethers_by_subtype(
      df = dfs[[i]], subtype = levels_subtype[[i]], polygons = polygons[[i]], 
      id = id, levels_ref_compartments = levels_ref_compartments, 
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      col_subtype = col_subtype, levels_subtype = levels_subtype, 
      exclude_inside = exclude_inside)
  }
  
  df_tether <- dplyr::bind_rows(df_tether)
  
  if (isTRUE(write_tetherome_res)) {
    readr::write_tsv(df_tether, file.path(filepath, out_name))
  }
}


#' Find tetherome by subtype.
#' 
#' @param df A data frame.
#' @param subtype A sample subtype.
#' @param polygons The vertices of polygons.
#' @param id An identifier.
#' @param col_subtype Column key of sample types.
#' @param levels_subtype The levels of sub_type.
#' @param col_fraction The column key in metadata to subcellular fractions.
#' @param levels_fraction The levels of subcellular compartments.
#' @param levels_ref_compartments The factor levels of subcellular
#'   compartments. Set \code{levels_ref_compartments = NULL} for automatic
#'   determination.
#' @param exclude_inside Logical; if TRUE, exclude points inside the polygon
#'   defined by \code{poly_x} and \code{poly_y}. Otherwise, exclude points
#'   outside of the polygon.
find_tethers_by_subtype <- function (
    df, subtype, polygons, id = "gene", levels_ref_compartments = NULL, 
    col_fraction = "Fraction", levels_fraction = NULL, col_subtype = "Type", 
    levels_subtype = NULL, exclude_inside = TRUE) 
{
  df_tether <- vector("list", length(levels_fraction))
  
  for (i in seq_along(levels_fraction)) {
    fr <- levels_fraction[[i]]
    pg <- polygons[[fr]]
    
    df <- hremove_contams(
      df = df, id = id, col_fraction = col_fraction, tether_fraction = fr, 
      polygons = pg, exclude_inside = exclude_inside)
    
    # The results for markers, not for the whole proteome
    df_tether[[i]] <- tibble::tibble(
      !!id := attr(df, "tethers"),
      !!col_subtype := subtype, 
      !!col_fraction := fr)
  }
  
  df_tether <- dplyr::bind_rows(df_tether)
  
  # Column not existed if all(is.null(polygons))
  ok_polygons <- !all(is.null(polygons))
  
  if (ok_polygons && id %in% colnames(df_tether)) {
    df_tether <- df_tether |>
      # NA gene if no polygon for a fraction
      # dplyr::filter(!is.na(.data[[id]])) |>
      dplyr::mutate(
        !!col_fraction := 
          factor(!!rlang::sym(col_fraction), levels = levels_fraction))
  } else {
    df_tether <- tibble::tibble(df_tether, !!id := NA_character_, )
  }
  
  df_tether
}


#' Help of plotting UMAP results.
#'
#' @param df A data frame combining UMAP, marker and trend results.
#' @param id An identifier.
#' @param col_subtype Column key of sample types.
#' @param levels_subtype The levels of sub_type.
#' @param levels_ref_compartments The factor levels of subcellular
#'   compartments. Set \code{levels_ref_compartments = NULL} for automatic
#'   determination.
#' @param fn_prefix A file name prefix.
#' @param fn_suffix A file name suffix.
#' @param ... Additional parameters for plot controls.
plot_umap <- function (
    df, polygons = NULL, draw_polygon = FALSE, highlight_markers = FALSE, 
    id = "gene", 
    col_pred = "sub_location", levels_pred = NULL, 
    col_fraction = NULL, levels_fraction = NULL, exclude_inside = TRUE, 
    col_subtype = NULL, levels_subtype = NULL, 
    levels_ref_compartments = NULL, write_tetherome_res = TRUE, 
    filepath = NULL, fn_prefix, fn_suffix, 
    break_width = 5, hide_grids = TRUE, wrap_facet = TRUE, drop = FALSE, 
    width_combined = 6, height_combined = 4, ext = "png", 
    width_sep = 14, height_sep = 6, 
    ...)
{
  dots <- rlang::enexprs(...)

  my_colors <- prep_subcell_colors(levels_pred)
  my_ref_colors <- prep_subcell_colors(levels_ref_compartments)

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
  
  xmin <- eval(dots$xmin)
  xmax <- eval(dots$xmax)
  ymin <- eval(dots$ymin)
  ymax <- eval(dots$ymax)
  
  if (is.null(xmin)) xmin <- -12
  if (is.null(xmax)) xmax <- 12
  if (is.null(ymin)) ymin <- -12
  if (is.null(ymax)) ymax <- 12
  
  dots$xmin <- dots$xmax <- dots$ymin <- dots$ymax <- NULL
  
  if (highlight_markers) {
    plot_sc_markers(
      df = df, polygons = polygons, draw_polygon = draw_polygon, 
      my_colors = my_ref_colors, id = id, col_x = col_group, 
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      col_subtype = col_subtype, levels_subtype = levels_subtype, 
      levels_ref_compartments = levels_ref_compartments, 
      dat_dir = dat_dir, filepath = filepath, ext = ext, 
      break_width = break_width, hide_grids = hide_grids, 
      wrap_facet = wrap_facet, drop = drop, 
      ncol = ncol, nrow = nrow, width = width_sep, height = height_sep, 
      dpi = dpi, 
      !!!dots)
    
    return(NULL)
  }
  
  p <- ggplot() +
    geom_point(
      data = df, 
      aes(x = UMAP1, y = UMAP2, color = !!rlang::sym(col_pred)), 
      shape = 46, size = 1, alpha = 1.0) +
    scale_color_manual(values = my_colors, breaks = levels_pred, 
                       drop = drop)
  
  p <- p +
    guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
    theme_bw() +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2",
      color = col_pred
    )
  
  if (hide_grids) {
    p <- p +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.title = element_blank(),
        panel.grid = element_blank())
  } else {
    p <- p +
      theme(
        plot.title = element_text(face = "bold"), 
        legend.title = element_blank())
  }
  
  p <- p + 
    guides(color = guide_legend(
      override.aes = list(size = 3, shape = 16, alpha = 1)))
  
  if (n_subtypes < 2L) {
    px <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_pred)), 
        ncol = ncol, labeller = label_value, drop = drop)
  } else {
    px <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_subtype), !!rlang::sym(col_pred)), ncol = ncol, 
        labeller = labeller(.default = label_value, .multi_line = FALSE), 
        drop = drop)
    
    p <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_subtype)),
        ncol = ncol, labeller = label_value, )
  }
  
  ggsave(plot = px, 
         file.path(filepath, paste0(fn_prefix, ".", fn_suffix)), 
         width = width_sep, height = height_sep, dpi = dpi)
  
  ggsave(plot = p, 
         file.path(filepath, paste0(fn_prefix, "_combined.", fn_suffix)), 
         width = width_combined, height = height_combined, dpi = dpi)
  
  invisible(df)
}


#' Prepare subcellular color schemes.
#' 
#' @param compartments A set of subcellular compartments.
prep_subcell_colors <- function (compartments)
{
  set1_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                   "#FFFF33", "#A65628", "#F781BF")
  add_colors  <- c("#008080", "#3A1188", "#808000", "#882255", "#00BFFF", 
                   "#99d8c9", "#fee6ce")
  last_color  <- "#999999"
  
  n_comps <- length(compartments)
  
  if (any(is.na(n_comps)) || "Unknown" %in% compartments) {
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
}


#' Helper in loading data of marker proteins.
#' 
#' @param df_markers A data frame of marker proteins.
#' @param species An species name.
#' @param levels_ref_compartments Levels of reference compartments.
load_ref_markers <- function (df_markers = NULL, species = "hs", 
                              levels_ref_compartments)
{
  abbr_species <- switch(
    species, 
    human = "hs",
    mouse = "mm",
    rat = "rn",
    stop("species: ", species, " is not in one of `human`, `mouse` or `rat`.")
  )
  
  if (is.null(df_markers)) {
    db_nm <- paste0("subcell_", abbr_species)
    data(list = db_nm)
    df_markers <- get(db_nm, rlang::current_env())
  }
  
  if (!"location" %in% names(df_markers)) {
    stop("Required column 'location' not found.")
  }
  
  if (FALSE && all(levels_ref_compartments %in% c("CP", "NP", "ChA")) && 
      "location_slim3" %in% names(df_markers) && 
      !all(is.na(df_markers[["location_slim3"]])) && 
      !identical(df_markers[["location"]], df_markers[["location_slim3"]])) {
    df_markers[["location"]] <- df_markers[["location_slim3"]]
  }
  
  if (is.null(levels_ref_compartments)) {
    locs <- df_markers[["location"]]
    univ <- unique(locs)
    levels_ref_compartments <- names(sort(table(locs), decreasing = TRUE))
  }
  
  df_markers <- df_markers |>
    dplyr::mutate(location = factor(location, levels = levels_ref_compartments))
}


#' Helper of separating data by sub-types.
#' 
#' @param df A data frame of expression data.
#' @param id An identifier.
#' @param scale Logical; scale data or not.
#' @inheritParams plot_prnSubcell_UMAP
sep_data_by_subtype <- function (
    file = NULL, df = NULL, id = "gene", filepath = NULL, label_scheme_sub, 
    col_group = NULL, col_fraction = NULL, 
    col_subtype = NULL, levels_subtype = NULL, 
    scale = TRUE) 
{
  if (is.null(df)) {
    df <- readr::read_tsv(file.path(filepath, file))
  }
  
  sc_lookup <- label_scheme_sub |>
    dplyr::select(dplyr::one_of(col_group, col_fraction, col_subtype)) |>
    unique()
  
  nms <- colnames(df)
  
  dfs <- vector("list", length(levels_subtype))
  names(dfs) <- levels_subtype
  
  for (i in seq_along(levels_subtype)) {
    sids_sub <- sc_lookup |>
      dplyr::filter(.data[[col_subtype]] == levels_subtype[[i]]) |>
      dplyr::pull(col_group)
    
    df_expr <- df[, nms %in% sids_sub]
    
    if (scale) {
      df_expr <- df_expr / rowSums(df_expr, na.rm = TRUE)
    }
    
    df_expr <- data.frame(df_expr, check.names = FALSE)
    rownames(df_expr) <- df[[id]]
    
    dfs[[i]] <- df_expr
  }
  
  dfs
}


#' Helper of plotting subcellular markers.
#' 
#' @param df A data frame.
#' @param id An identifier.
#' @param df Background data. 
#' @param col_x The column key of X-axis.
#' @param draw_polygon Logical; draw polygon or not.
#' @param col_subtype Column key of sample types.
#' @param levels_subtype The levels of sub_type.
#' @param col_fraction The column key in metadata to subcellular fractions.
#' @param levels_fraction The levels of subcellular compartments.
#' @param levels_ref_compartments The factor levels of subcellular
#'   compartments. Set \code{levels_ref_compartments = NULL} for automatic
#'   determination.
#' @param polygons The vertices of polygons.
#' @param filepath An output file path.
#' @param ext An file name extension.
#' @param break_width The width of X and Y breaks.
#' @param hide_grids Logical; hide grid lines or not.
#' @param wrap_facet Logical; wrap_facet or not.
#' @param drop Drop unknown levels or not.
#' @param my_colors A color scheme.
#' @param dat_dir A workding directory.
#' @param filepath A file path to inputs and outputs.
plot_sc_markers <- function (
    df, polygons = NULL, draw_polygon = TRUE, my_colors = NULL, 
    id = "gene", col_x = "Group", 
    col_fraction = "Fraction", levels_fraction = NULL, 
    col_subtype = "Type", levels_subtype = NULL, 
    levels_ref_compartments = NULL, 
    dat_dir = NULL, filepath = NULL, ext = "png", 
    ncol = 7L, nrow = 2L, width = 14, height = 3, dpi = 600, 
    break_width = 5L, hide_grids = TRUE, wrap_facet = TRUE, drop = FALSE, 
    ...)
{
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  
  dfs <- split(df, df[[col_subtype]])
  dfs <- dfs[levels_subtype]
  polygons <- polygons[levels_subtype]
  draw_polygon <- draw_polygon[levels_subtype]
  
  for (i in seq_along(levels_subtype)) {
    plot_sc_markers_by_subtype(
      df = dfs[[i]], subtype = levels_subtype[[i]], polygons = polygons[[i]], 
      draw_polygon = draw_polygon[[i]], my_colors = my_colors, 
      id = id, col_x = col_x, 
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      col_subtype = col_subtype, levels_subtype = levels_subtype, 
      levels_ref_compartments = levels_ref_compartments, 
      dat_dir = dat_dir, filepath = filepath, ext = ext, 
      ncol = ncol, nrow = nrow, 
      width = width, height = height, dpi = dpi, 
      break_width = break_width, hide_grids = hide_grids, 
      wrap_facet = wrap_facet, drop = drop, ...)
  }
}


#' Plot subcellular markers by subtype.
#' 
#' @param df A data frame.
#' @param subtype A sample subtype.
#' @param id An identifier.
#' @param df Background data. 
#' @param col_x The column key of X-axis.
#' @param draw_polygon Logical; draw polygon or not.
#' @param col_subtype Column key of sample types.
#' @param levels_subtype The levels of sub_type.
#' @param col_fraction The column key in metadata to subcellular fractions.
#' @param levels_fraction The levels of subcellular compartments.
#' @param levels_ref_compartments The factor levels of subcellular
#'   compartments. Set \code{levels_ref_compartments = NULL} for automatic
#'   determination.
#' @param polygons The vertices of polygons.
#' @param filepath An output file path.
#' @param ext An file name extension.
#' @param break_width The width of X and Y breaks.
#' @param hide_grids Logical; hide grid lines or not.
#' @param wrap_facet Logical; wrap_facet or not.
#' @param drop Drop unknown levels or not.
#' @param my_colors A color scheme.
#' @param dat_dir A workding directory.
#' @param filepath A file path to inputs and outputs.
plot_sc_markers_by_subtype <- function (
    df, subtype, polygons = NULL, draw_polygon = TRUE, my_colors = NULL, 
    id = "gene", col_x = "Group", 
    col_fraction = "Fraction", levels_fraction = NULL, 
    col_subtype = "Type", levels_subtype = NULL, 
    levels_ref_compartments = NULL, 
    dat_dir = NULL, filepath = NULL, ext = "png", 
    ncol = 7L, nrow = 2L, width = 14, height = 3, dpi = 600, 
    break_width = 5L, hide_grids = TRUE, wrap_facet = TRUE, drop = FALSE, ...)
{
  dots <- rlang::enexprs(...)
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  
  draw_polygon <- isTRUE(draw_polygon)
  
  # Expression subset found with marker annotations
  df_loc <- df |>
    dplyr::filter(!is.na(.data[["location"]]))
  
  n_fracs <- length(levels_fraction)

  xmin <- eval(dots$xmin)
  xmax <- eval(dots$xmax)
  ymin <- eval(dots$ymin)
  ymax <- eval(dots$ymax)
  
  if (is.null(xmin)) xmin <- -12
  if (is.null(xmax)) xmax <- 12
  if (is.null(ymin)) ymin <- -12
  if (is.null(ymax)) ymax <- 12
  
  dots$xmin <- dots$xmax <- dots$ymin <- dots$ymax <- NULL
  
  top_layers_1  <- df_loc[["location"]] %in% c("Centrosome", "Chromatin")
  top_layers_2  <- df_loc[["location"]] %in% c("Nucleus")
  df_loc_top_1  <- df_loc[top_layers_1, ]
  df_loc_top_2  <- df_loc[top_layers_2, ]
  df_loc_others <- df_loc[!(top_layers_1 | top_layers_2), ]
  
  p <- ggplot() +
    geom_point(data = df, aes(x = UMAP1, y = UMAP2), color = "gray", 
               shape = 1, size = 1.5, stroke = .2)
  
  p <- p + 
    geom_point(data = df_loc_others, aes(x = UMAP1, y = UMAP2, fill = location), 
               shape = 21, color = "white", size = 1.5, , stroke = .1, alpha = 1) + 
    geom_point(data = df_loc_top_2, aes(x = UMAP1, y = UMAP2, fill = location), 
               shape = 21, color = "white", size = 1.5, , stroke = .1, alpha = 1) + 
    geom_point(data = df_loc_top_1, aes(x = UMAP1, y = UMAP2, fill = location), 
               shape = 21, color = "white", size = 1.5, , stroke = .1, alpha = 1) + 
    scale_fill_manual(values = my_colors, drop = drop)
  
  p <- p + 
    scale_x_continuous(
      limits = c(xmin, xmax), breaks = scales::breaks_width(break_width)) +
    scale_y_continuous(
      limits = c(ymin, ymax), breaks = scales::breaks_width(break_width))
  
  p <- p + 
    theme_minimal(base_size = 14) +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2")
  
  if (hide_grids) {
    p <- p +
      theme(
        strip.text = element_text(
          size = 14,
          color = "black"
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), 
        legend.title = element_blank()
      ) 
  } else {
    p <- p +
      theme(
        strip.text = element_text(
          size = 14,
          color = "black"
        ),
        plot.title = element_text(face = "bold"), 
        legend.title = element_blank()
      ) 
  }
  
  p <- p +
    guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))
  
  if (wrap_facet) {
    p <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_fraction)), 
        ncol = ncol, 
        labeller = labeller(.default = label_value, .multi_line = FALSE))
  }
  
  if (draw_polygon) {
    # k <- 0.08685714; b <- 1.148571
    df_poly <- vector("list", n_fracs)
    names(df_poly) <- levels_fraction
    
    for (i in seq_along(levels_fraction)) {
      fr <- levels_fraction[[i]]
      pg <- polygons[[fr]]
      
      if (length(pg)) {
        stopifnot(length(pg[["x"]]) == length(pg[["y"]]))
      }
      
      df_poly[[i]] <- tibble::tibble(
        x = pg[["x"]], y = pg[["y"]], !!col_fraction := fr)
    }
    
    df_poly <- dplyr::bind_rows(df_poly) |>
      dplyr::mutate(
        !!col_fraction := 
          factor(!!rlang::sym(col_fraction), levels = levels_fraction))
    
    p <- p + 
      geom_polygon(
        df_poly, 
        mapping = aes(x = x, y = y), 
        fill = NA, color = "black", linewidth = .2)
    
    # k <- 10/3; b <- -28/3; 1.5 * k + b
    if (FALSE) {
      p <- p + 
        geom_segment(
          # x = -0.2, y = -10, xend = 1.5, yend = -4.333
          # x = 1.5, y = -4.333, xend = -0.2, yend = -10
          data = data.frame(
            x = poly_x_np_cha[[1]], y = poly_y_np_cha[[1]], 
            xend = poly_x_np_cha[[2]], yend = poly_y_np_cha[[2]]),
          aes(x = x, y = y, xend = xend, yend = yend),
          linetype = "longdash",
          color = "darkgray",
          linewidth = 0.5)
    }
  }
  
  if (draw_polygon) 
    out_name <- paste0("umap_markers_", subtype, "_polygon.", ext)
  else 
    out_name <- paste0("umap_markers_", subtype, ".", ext)
  
  ggsave(filename = file.path(filepath, out_name), 
         width = width, height = height, dpi = dpi)
  
  p
}


#' Update subcellular protein or peptide tables by contrast.
#' 
#' A contrast can be nested.
#' 
#' @param contr_elements Contrast elements.
#' @param contr_nm Contrast names.
#' @param FUN A logical function of \code{all} or \code{any}.
#' @param df A data frame.
update_subell_by_contr <- function (
    contr_elements = NULL, contr_nm = NULL, df = NULL, FUN = all)
{
  # No assigned subcellular locations
  nas <- local({
    dfx <- data.frame(df[, contr_elements], check.names = FALSE)
    
    sapply(1:nrow(dfx), function(i) {
      vals <- dfx[i, ]
      FUN(is.na(vals) | vals == 0)
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
#' @param FUN A logical function of \code{all} or \code{any}.
#' @param label_scheme_sub Metadata.
update_subcell_by_fml <- function (
    fml_nm = NULL, fml = NULL, df = NULL, label_scheme_sub = NULL, FUN = all)
{
  fml_ops    <- prepFml(formula = fml, label_scheme_sub = label_scheme_sub)
  elements   <- fml_ops[["elements"]]
  contr_mat  <- fml_ops[["contr_mat"]]
  sub_contrs <- fml_ops[["sub_contrs"]] # NULL if no nested structure
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
      contr_elements = contr_elements, contr_nm = contr_nm, df = df, FUN = FUN)
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
#' @param rm_tetherome Logical; remove tetherome entries or not by nullifying
#'   pVals.
#' @param dat_dir The working directory.
#' @param data_type The type of data in either \code{Protein} or \code{Peptide}.
#' @param filepath Automatically determined: the fiepath to \code{Protein} or
#'   \code{Peptide} table.
#' @param fml_nms Automatically determined: the formula names used in
#'   \code{prnSig} or \code{pepSig}.
#' @param label_scheme_sub Automatically determined: metadata.
#' @export
update_subcell <- function (
    data_type = c("Protein", "Peptide"), rm_tetherome = FALSE, df = NULL, 
    dat_dir = NULL, filepath = NULL, fml_nms = NULL, label_scheme_sub = NULL, 
    ...)
{
  data_type <- match.arg(data_type)
  filename  <- paste0(data_type, "_pVals.txt")
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, data_type, "Model")

  group_psm_by <- tryCatch(
    match_call_arg("normPSM", "group_psm_by"), error = function(e) NA)
  
  if (is.na(group_psm_by)) {
    group_psm_by <- tryCatch(
      match_call_arg("makePepDIANN", "group_psm_by"), 
      error = function(e) NA)
  }
  
  if (is.na(group_psm_by)) {
    group_psm_by <- "pep_seq_mod"
  }
  
  group_pep_by <- tryCatch(
    match_call_arg("normPSM", "group_pep_by"), error = function(e) NA)
  
  if (is.na(group_pep_by)) {
    group_pep_by <- tryCatch(
      match_call_arg("makePepDIANN", "group_psm_by"), 
      error = function(e) NA)
  }
  
  if (is.na(group_pep_by)) {
    group_pep_by <- "gene"
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
    is_subcellular = TRUE, rm_tetherome = rm_tetherome, 
    update_trend_res = if (data_type == "Protein") TRUE else FALSE, 
    scale_log2r = TRUE, impute_na = FALSE)
  df <- res_subcell[["df"]]
  col_fraction <- res_subcell[["col_fraction"]]
  col_subtype <- res_subcell[["col_subtype"]]
  col_subgroup <- res_subcell[["col_group"]]
  rm(list = "res_subcell")
  
  for (i in seq_along(fml_nms)) {
    fml_nm <- fml_nms[[i]]
    fml <- fmls[[i]]
    
    df <- update_subcell_by_fml(
      fml_nm = fml_nm, fml = fml, df = df, label_scheme_sub = label_scheme_sub, 
      FUN = all)
  }
  
  df <- df[, names(df) %in% nms_df]
  
  path_out <- file.path(dat_dir, data_type, "Model", filename)
  readr::write_tsv(df, path_out)
  
  invisible(df)
}


#' Update pVals and log2Ratios by incorporating subcellular specificity.
#'
#' Nullifying p-values and log2Ratios of protein and peptides
#'
#' @param rm_tetherome Logical; if TRUE, nullify additionally pVals and
#'   log2Ratios of tetherome entries.
#' @export
update_subcell_stat <- function (rm_tetherome = TRUE)
{
  update_subcell(data_type = "Protein", rm_tetherome = rm_tetherome)
  update_subcell(data_type = "Peptide", rm_tetherome = rm_tetherome)
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
    unique() # |> dplyr::rename(Sample_ID := col_group)
  
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
    levels_subtype = NULL, transloc_co = .5, 
    use_updated_subcell = TRUE, 
    col_select = "Select", col_order = "Order", # col_group = "Group", 
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
  col_group    <- df[["col_group"]][[1]]
  col_fraction <- df[["col_fraction"]][[1]]
  col_subtype  <- df[["col_subtype"]][[1]]
  
  if (!length(levels_subtype)) {
    levels_subtype <- unique(df[[col_subtype]])
  }
  
  n_subtype <- length(levels_subtype)
  if (n_subtype != 2L) {
    stop("Proive two values to 'levels_subtype' for ", 
         "pairwised translocation analysis.")
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
  
  invisible(df_combined)
}


#' Update subcellular purity, entropy and localization scores.
#'
#' Also update trend results under sub directory \code{subcell_update}.
#'
#' @param id An identifier.
#' @param dat_dir The current working directory.
#' @param filepath The filepath to protein \code{Trend} results.
#' @param files The filenames to protein \code{Trend} results.
#' @param col_select Not used. The column key in metadata linked to the
#'   information of sample being selected.
#' @param col_group The column key in metadata linked to the information of
#'   sample groups.
#'
#'   Not used since its values can be changed from run to run by users. To
#'   decide, use the first values under the column \code{"col_group"}.
#' 
#' @param col_oder The column key in metadata linked to the information of
#'   sample orders.
#' @param col_fraction The column key in metadata linked to the information of
#'   sample subcellular fractions.
#' @param col_subtype The column key in metadata linked to the information of
#'   sample subtypes or treatment conditions.
#'   
#'   Not used. See also parameter \code{col_group}.
#'   
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
    col_select = NULL, col_order = NULL, 
    # col_group = NULL, col_fraction = NULL, col_subtype = NULL, 
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
  
  col_group    <- df[["col_group"]][[1]]
  col_fraction <- df[["col_fraction"]][[1]]
  col_subtype  <- df[["col_subtype"]][[1]]
  
  if (!length(levels_subtype)) {
    levels_subtype <- unique(df[[col_subtype]])
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
  df_mean_int <- df[, c(id, col_group, "log10Int")] |>
    dplyr::mutate(log10Int = 10^log10Int) |>
    tidyr::pivot_wider(
      id_cols = id, 
      names_from = col_group, 
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
    unique()

  purities <- vector("list", n_subtypes)
  # DON'T: names(purities) <- levels_subtype
  for (i in seq_along(levels_subtype)) {
    sids_sub <- sc_lookup |>
      dplyr::filter(.data[[col_subtype]] == levels_subtype[[i]]) |>
      dplyr::pull(col_group)

    dfy <- df_mean_int[, grps %in% sids_sub]
    purities[[i]] <- dfy / rowSums(dfy, na.rm = TRUE)
  }
  rm(list = c("dfy", "sids_sub"))
  
  entropies  <- lapply(purities, function (x) {
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
  
  purities   <- lapply(purities, function (x) {
    x[x < 1E-3] <- 0.0
    x
  })
  
  entropies  <- lapply(entropies, function (x) {
    x[x < 1E-3] <- 0.0
    x
  })
  
  purities   <- do.call(cbind, purities)
  scores_loc <- do.call(cbind, scores_loc)
  entropies  <- do.call(cbind, entropies)
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
      names_to = col_group,
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
      names_to = col_group,
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
      purrr::reduce(dplyr::left_join, by = c(id, col_group)) |>
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


#' Plot marker proteins onto the scatter of 2-D UMAP.
#'
#' @param id An identifier.
#' @param species The species.
#' @param data_type The type of data.
#' @param df_markers A data frame of marker proteins.
#' @param dat_dir The current working directory.
#' @param filepath The file path.
#' @param label_scheme Metadata.
#' @param levels_group The levels of data group.
#' @param col_subtype Column key of sample types.
#' @param levels_subtype The levels of sub_type.
#' @param col_fraction The column key in metadata to subcellular fractions.
#' @param levels_fraction The levels of subcellular compartments.
#' @param levels_ref_compartments The factor levels of subcellular
#'   compartments. Set \code{levels_ref_compartments = NULL} for automatic
#'   determination.
#' @param polygons The vertices of polygons.
#' @param pat_df2 The pattern of secondary files (df2).
#' @param ext A file name extension.
#' @param break_width The width of breaks on X and Y axis.
#' @param hide_grids Logical; hide grid lines or not.
#' @param ... Additional parameters for plot controls.
#' @param exclude_inside Logical; if TRUE, exclude points inside the polygon
#'   defined by \code{poly_x} and \code{poly_y}. Otherwise, exclude points
#'   outside of the polygon.
#' 
#' @importFrom magrittr %>%
#' @export
plot_subcell_markers <- function (
    id = "gene", species = "human", 
    # data_type = c("expr", "prob"), 
    method = c("unsupervised", "supervised"), 
    df_markers = NULL, dat_dir = NULL, filepath = NULL, 
    label_scheme = NULL, 
    levels_ref_compartments = c(
      "Cytosol", "ER", "Mitochondrion", "Nucleus", "Chromatin", "Lysosome", 
      "Plasma_membrane", "Proteasome", "Golgi", "Peroxisome", "Endosome", 
      "Ribosome", "Centrosome", "P_body", "Stress_granule"), 
    col_group = NULL, levels_group = NULL, 
    col_subtype = NULL,  levels_subtype = NULL, 
    col_fraction = NULL, levels_fraction = NULL, 
    polygons = NULL, exclude_inside = TRUE, 
    ext = "png", break_width = 5, hide_grids = TRUE, 
    ...)
{
  dots <- rlang::enexprs(...)
  
  write_tetherome_res <- TRUE
  drop <- FALSE
  wrap_facet <- TRUE
  
  dot_write_tetherome_res <- dots[["write_tetherome_res"]]
  dot_drop <- dots[["drop"]]
  dot_wrap_facet <- dots[["wrap_facet"]]
  dots[["write_tetherome_res"]] <- dots[["drop"]] <- dots[["wrap_facet"]] <- NULL
  
  if (!is.null(dot_write_tetherome_res) && 
      is.logical(is.null(dot_write_tetherome_res)) && 
      !dot_write_tetherome_res) {
    write_tetherome_res <- FALSE
  }
  
  if (!is.null(dot_drop) && is.logical(dot_drop) && dot_drop) {
    drop <- TRUE
  }
  
  if (!is.null(dot_wrap_facet) && is.logical(dot_wrap_facet) && !dot_wrap_facet) {
    wrap_facet <- FALSE
  }
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  if (is.null(filepath)) filepath <- file.path(dat_dir, "Protein", "Trend")
  
  message("The output path: ", filepath)
  
  ### 1. Marker proteins
  abbr_species <- switch(
    species, 
    human = "hs",
    mouse = "mm",
    rat = "rn",
    stop("species: ", species, " is not in one of `human`, `mouse` or `rat`.")
  )
  
  if (is.null(df_markers)) {
    db_nm <- paste0("subcell_", abbr_species)
    data(list = db_nm)
    df_markers <- get(db_nm, rlang::current_env())
  }
  
  if (!"location" %in% names(df_markers)) {
    stop("Required column 'location' not found.")
  }
  
  if (is.null(levels_ref_compartments)) {
    levels_ref_compartments <- local({
      locs <- df_markers[["location"]]
      univ <- unique(locs)
      levs <- names(sort(table(locs), decreasing = TRUE))
    })
  }
  
  df_markers <- df_markers |>
    dplyr::mutate(location = factor(location, levels = levels_ref_compartments))
  
  ### 
  method <- match.arg(method)
  
  ### 2. Metadata
  if (is.null(label_scheme)) label_scheme <- load_ls_group(dat_dir = dat_dir)
  
  label_scheme_sub <- label_scheme |>
    dplyr::filter(.data[[col_subtype]] %in% levels_subtype)
  
  if (is.null(levels_group)) {
    levels_group <- label_scheme_sub |>
      dplyr::arrange_at("Order") |>
      dplyr::pull("Group") |>
      unique()
  }
  
  ### 3. Expression data
  if (method == "supervised") {
    pat_df2 <- "svmProb_[ONZ]_.*nclust\\d+\\.tsv$"
    ext_df2 <- "tsv"
  } else if (method == "unsupervised") {
    pat_df2 <- "Trend_[ONZ]_.*nclust\\d+\\.txt$"
    ext_df2 <- "txt"
  } else {
    pat_df2 <- "Trend_[ONZ]_.*nclust\\d+\\.txt$"
    ext_df2 <- "txt"
  }
  
  file_df2 <- find_trend_df2(pat = pat_df2, ext = ext_df2)[[1]]
  df <- readr::read_tsv(file.path(dat_dir, "Protein", "Trend", file_df2))
  
  if (method == "unsupervised") {
    col_fraction <- df[["col_fraction"]][[1]]
    col_subtype  <- df[["col_subtype"]][[1]]
    col_group    <- df[["col_group"]][[1]]
    
    key_col <- "log2FC"
    col_x   <- col_group
  } else if (data_type == "supervised") {
    key_col <- "log2FC"
    col_x   <- col_group
  } else {
    col_fraction <- df[["col_fraction"]][[1]]
    col_subtype  <- df[["col_subtype"]][[1]]
    col_group    <- df[["col_group"]][[1]]
    
    key_col <- "log2FC"
    col_x   <- col_group
  }
  

  df <- df |>
    dplyr::filter(.data[[col_x]] %in% levels_group, 
                  .data[[col_subtype]] %in% levels_subtype, 
                  .data[[col_fraction]] %in% levels_fraction) |>
    dplyr::mutate(
      !!col_x := 
        factor(!!rlang::sym(col_x), levels = levels_group), 
      !!col_subtype := 
        factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
      !!col_fraction := 
        factor(!!rlang::sym(col_fraction), levels = levels_fraction), ) |>
    dplyr::arrange_at(c(col_x, col_fraction, col_subtype)) |>
    dplyr::filter(!is.na(.data[["sub_location"]]), ) |>
    dplyr::left_join(df_markers, by = id)
  
  dfs <- split(df, df[[col_subtype]])
  dfs <- dfs[levels_subtype]
  
  ### 5. UMAP data
  files_umap <- local({
    files_umap <- file.path(filepath, paste0("umap_", levels_subtype, ".tsv"))
    
    oks <- file.exists(files_umap)
    
    if (any(bads <- !oks)) {
      stop("UMAP results not found: ", paste(files_umap[bads], collapse = ", "))
    }
    
    files_umap <- files_umap[levels_subtype]
    names(files_umap) <- levels_subtype
    
    if (!length(files_umap)) {
      stop("UMAP results not found. Run 'plot_prnSubcell_UMAP' first.")
    }
    
    files_umap
  })
  
  dfs_umap   <- lapply(files_umap, readr::read_tsv)
  dfs_umap   <- lapply(dfs_umap, function (x) { x[[col_subtype]] <- NULL; x})
  names(dfs_umap) <- levels_subtype
  
  ## Check polygons
  n_subtypes <- length(levels_subtype)
  
  if (!all(names(polygons) %in% levels_subtype)) {
    stop("Mismatched levels between polygons and levels_subtype")
  }
  
  polygons <- polygons[levels_subtype]
  
  mapply(function (ps, subtype) {
    oks <- lapply(lapply(ps, lengths), function (n) {
      n[["x"]] == n[["y"]]
    }) |>
      unlist()
    
    if (any(bads <- !oks)) {
      stop("Mismatched lengths of x and y vertices: ", subtype, " | ", 
           paste(names(bads)[bads == TRUE], collapse = ", "))
    }
  }, polygons, levels_subtype)
  
  draw_polygon <- lapply(polygons, function (x) !is.null(x))
  
  if (!length(draw_polygon)) {
    draw_polygon <- rep_len(FALSE, n_subtypes)
  }

  df_tether <- vector("list", n_subtypes)
  
  
  
  for (i in seq_along(levels_subtype)) {
    df_tether[[i]] <- find_tetherome_by_subtype(
      df = dfs[[i]], 
      df_umap = dfs_umap[[i]], 
      polygons = polygons[[i]], 
      draw_polygon = draw_polygon[[i]], 
      subtype = levels_subtype[[i]], 
      id = id, 
      col_fraction = col_fraction,
      levels_fraction = levels_fraction, 
      col_subtype = col_subtype, 
      levels_subtype = levels_subtype, 
      levels_ref_compartments = levels_ref_compartments, 
      exclude_inside = exclude_inside)
  }
  
  for (i in seq_along(levels_subtype)) {
    plot_subcell_markers_by_subtype(
      df = dfs[[i]], 
      df_umap = dfs_umap[[i]], 
      polygons = polygons[[i]], 
      draw_polygon = draw_polygon[[i]], 
      subtype = levels_subtype[[i]], 
      id = id, col_x = col_x, 
      # df_bg = df, 
      col_fraction = col_fraction,
      levels_fraction = levels_fraction, 
      col_subtype = col_subtype, 
      levels_subtype = levels_subtype, 
      levels_ref_compartments = levels_ref_compartments, 
      filepath = filepath, 
      ext = ext, 
      break_width = break_width, 
      hide_grids = hide_grids, 
      wrap_facet = wrap_facet, 
      drop = drop, 
      !!!dots
    )
  }
  
  df_tether <- dplyr::bind_rows(df_tether)
  
  if (isTRUE(write_tetherome_res)) {
    readr::write_tsv(df_tether, file.path(filepath, "marker_tetherome.tsv"))
  }

  df_tether
}


#' Helper remove tethered mitochondrial, ER proteins.
#' 
#' @param df A data frame.
#' @param col_fraction The column key in metadata to subcellular fractions.
#' @param tether_fraction The fraction where tetherome will be considered.
#' @param polygons The vertices of polygons.
#' @param exclude_inside Logical; if TRUE, exclude points inside the polygon
#'   defined by \code{poly_x} and \code{poly_y}. Otherwise, exclude points
#'   outside of the polygon.
hremove_contams <- function (
    df, id = "gene", col_fraction = "Fraction", tether_fraction = "NP",
    polygons = NULL, exclude_inside = TRUE)
{
  rows <- df[[col_fraction]] == tether_fraction
  
  dfx <- remove_contams(
    df[rows, ], 
    id = id, 
    poly_x = polygons[["x"]], 
    poly_y = polygons[["y"]], 
    exclude_inside = exclude_inside)
  
  tethers <- attr(dfx, "tethers")
  
  df <- dplyr::bind_rows(df[!rows, ], dfx)
  attr(df, "tethers") <- tethers
  
  df
}


#' Remove tethered mitochondrial, ER proteins.
#'
#' @param df A data frame.
#' @param An identifiers.
#' @param poly_x The vertices of X coordinates of a polygon.
#' @param poly_x The vertices of Y coordinates of a polygon.
#' @inheritParams hremove_contams
remove_contams <- function (
    df, id = "gene", poly_x = NULL, poly_y = NULL, exclude_inside = TRUE)
{
  if (is.null(poly_x) || is.null(poly_y)) {
    return(df)
  }
  
  inside_status <- sp::point.in.polygon(
    point.x = df[["UMAP1"]], 
    point.y = df[["UMAP2"]], 
    pol.x   = poly_x, 
    pol.y   = poly_y
  )
  
  ids_inside  <- unique(df[[id]][inside_status > 0])
  ids_outside <- unique(df[[id]][inside_status <= 0])
  
  rows <- df[[id]] %in% ids_inside
  
  df_inside  <- df[rows, ]
  df_outside <- df[!rows, ]
  
  if (exclude_inside) {
    res <- df_outside
    attr(res, "tethers") <- ids_inside
  } else {
    res <- df_inside
    attr(res, "tethers") <- ids_outside
  }
  
  res
}


