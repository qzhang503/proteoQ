#' Prepare Bandle.
#' 
#' @param df_expr_sub An expression data by subtype.
#' @param comartment_2 A set of markers.
prep_bandle <- function (df_expr_sub, location_slim3)
{
  pheno_data <- data.frame(
    FractionNumber = seq_len(ncol(df_expr_sub)),
    row.names = colnames(df_expr_sub)
  )
  
  feature_data <- data.frame(
    protein_id = rownames(df_expr_sub),
    markers = location_slim3,
    stringsAsFactors = FALSE
  )
  
  rownames(feature_data) <- feature_data$protein_id
  MSnbase::MSnSet(exprs = df_expr_sub, fData = feature_data, pData = pheno_data)
}


#' Calculate Bandle.
#' 
#' @param my_data Expression data.
#' @param sub_grp A sub group.
#' @param dat_dir A working directory.
#' @param use_svm Use SVM or not.
calc_bandle <- function(my_data, sub_grp, dat_dir, use_svm = FALSE)
{
  ## BANDLE
  # 1. Train the MAP model parameters
  map_par <- pRoloc::tagmMapTrain(my_data, fcol = "markers")
  
  # 2. Run the prediction using those parameters
  map_res <- pRoloc::tagmPredict(my_data, params = map_par, fcol = "markers")
  
  # Apply a threshold
  tagm_thresholded <- pRoloc::getPredictions(
    map_res, fcol = "tagm.map.allocation", scol = "tagm.map.probability", 
    t = 0)
  
  ans_bayes <- Biobase::fData(
    tagm_thresholded)[, c("protein_id", "tagm.map.allocation", "tagm.map.probability")]
  
  # table(ans_bayes$tagm.map.allocation)
  
  readr::write_tsv(
    ans_bayes, 
    file.path(dat_dir, "Protein", "Bandle", 
              paste0("bandle_pred_", sub_dir, "_", sub_grp, ".tsv")))

  ## SVM
  if (use_svm) {
    # 1. Optimize the SVM parameters (sigma and cost) using a grid search
    set.seed(123)
    # 1: In lapply(X = X, FUN = FUN, ...) :
    #   NaNs found in 'precision' with hyperparameters cost:16 sigma:1.
    # 2: In lapply(X = X, FUN = FUN, ...) :
    #   NaNs found in 'precision' with hyperparameters cost:8 sigma:10.
    svm_param <- pRoloc::svmOptimisation(my_data, fcol = "markers", times = 10)
    
    # Define a constrained hyperparameter grid
    if (FALSE) {
      svm_param <- pRoloc::svmOptimisation(
        my_data, 
        fcol = "markers", 
        times = 10,
        hyperpar = list(
          cost = 2^(0:4),        # e.g., 1, 2, 4, 8, 16
          sigma = 10^(-2:1)      # e.g., 0.01, 0.1, 1, 10
        )
      )
    }
    
    # 2. Run the final classifier using the best parameters
    svm_res <- pRoloc::svmClassification(my_data, svm_param, fcol = "markers")
    # plot2D(svm_res, fcol = "svm")
    
    readr::write_tsv(
      Biobase::fData(svm_res), 
      file.path(dat_dir, "Protein", "Bandle", 
                paste0("svm_pred_", sub_dir, "_", sub_grp, ".tsv")))
  }
  
  if (use_svm) {
    ans <- list(bandle = ans_bayes, svm = svm_res)
  } else {
    ans <- list(bandle = ans_bayes, svm = NULL)
  }
  
  ans
}


#' Analyze Bandle.
#' 
#' @param work_dir A parent working directory. 
#' @param sub_grp A sub group.
#' @param dat_dir A working directory.
#' @param use_svm Use SVM or not.
anal_bandle <- function (work_dir, sub_dir, db_path, dat_dir, 
                         subgrp_prefix = "", species = "human", 
                         col_group = "Group", 
                         levels_fraction = c("CP", "NP", "ChA"), 
                         levels_subtype = c("Veh", "TG"), 
                         use_svm = TRUE, ...) 
{
  ## Bandle
  dir.create(file.path(dat_dir, "Protein", "Bandle"))
  
  df_sc <- readr::read_tsv(
    file.path(db_path, paste0("bandle_", species, "_single_subcellular.tsv"))
  )
  
  df <- readr::read_tsv(
    file.path(dat_dir, "Protein", "Trend", 
              paste0(subgrp_prefix, "Protein_Trend_Z_nclust32.txt")))
  
  df_expr <- df[, c("gene", col_group, "log10Int")] |>
    dplyr::mutate(Int = 10^log10Int) |>
    dplyr::select(-"log10Int") |> 
    tidyr::pivot_wider(names_from = !!rlang::sym(col_group), values_from = Int, 
                       values_fill = 0) |>
    dplyr::left_join(df_sc, by = "gene") |>
    dplyr::mutate(location = ifelse(is.na(location), "unknown", location),
                  location_slim3 = ifelse(is.na(location_slim3), 
                                          "unknown", location_slim3),)
  
  # Convert to data matrix
  genes <- df_expr$gene
  locations <- df_expr$location
  location_slim3 <- df_expr$location_slim3
  df_expr <- as.matrix(
    df_expr[, -which(names(df_expr) %in% c("gene", "location", "location_slim3"))])
  rownames(df_expr) <- genes
  df_expr <- df_expr / rowSums(df_expr)
  
  # Can have NaN if only select a subset of input data
  # df_expr[is.na(df_expr)] <- 0.0
  oks <- rowSums(is.na(df_expr)) == 0L
  df_expr <- df_expr[oks, ]
  genes <- genes[oks]
  locations <- locations[oks]
  location_slim3 <- location_slim3[oks]
  
  # Separate into sub-groups
  n_subtypes <- length(levels_subtype)
  df_exprs <- vector("list", n_subtypes)
  for (i in seq_along(levels_subtype)) {
    df_exprs[[i]] <- 
      df_expr[, grepl(paste0("^", levels_subtype[[i]], "_"), colnames(df_expr))]
  }
  names(df_exprs) <- levels_subtype
  rm(list = "df_expr")
  
  markers <- unique(df_sc$location_slim3)
  bandle_inputs <- lapply(df_exprs, prep_bandle, location_slim3)
  
  res <- vector("list", length(bandle_inputs))
  nms_bandle <- names(bandle_inputs)
  for (i in seq_along(bandle_inputs)) {
    res_i <- calc_bandle(bandle_inputs[[i]], nms_bandle[[i]], dat_dir = dat_dir)
    
    res[[i]] <- if (use_svm) res_i[["svm"]] else res_i[["bandle"]]
  }

  res
}


#' Helper finding the confusion matrix between Bandle and proteoQ.
#'
#' @param df_pq Subcellular classification by proteoQ with unclassified entries
#'   under each of experimental CP, NP and ChA fractions already being removed.
#' @param df_bandle Subcellular classification by Bandle.
#' @param id A column key of identifier.
#' @param col_subtype A column key of sample subtypes, e.g., Veh, TG.
#' @param col_fraction A column key of sample fractions, e.g., CP, NP and ChA.
#' @param bandle_pred_col A column key of Bandle subcellular prediction.
confu_bandle <- function (
    df_pq, df_bandle, id = "gene", col_subtype = "Type", 
    col_fraction = "Fraction", bandle_pred_col = "pred_bandle", 
    remove_tethers = FALSE, poly_x_cp_np, poly_y_cp_np, 
    poly_x_np_cha, poly_y_np_cha) 
{
  df_pqs <- split(df_pq, df_pq[[col_subtype]])
  df_banles <- split(df_bandle, df_bandle[[col_subtype]])
  
  df_inners <- df_pq[, c(id, col_fraction, col_subtype)] |>
    dplyr::inner_join(df_bandle, by = c(id, col_subtype)) 
  df_inners <- split(df_inners, df_inners[[col_subtype]])

  ans_confu <- lapply(
    df_inners, confu_bandle_by_frac, id = id, col_fraction = col_fraction, 
    bandle_pred_col = bandle_pred_col)
  
  mat_confu <- vector("list", length(ans_confu))
  for (i in seq_along(ans_confu)) {
    mat_confu[[i]] <- 
      table(ans_confu[[i]][[col_fraction]], ans_confu[[i]][[bandle_pred_col]])
  }
  
  attr(ans_confu, "mat_confu") <- mat_confu
  
  ans_confu
}


#' Finding the confusion matrix between Bandle and proteoQ by each subtype.
#'
#' @param df_inner Joint data of Bandle subcellular classification by
#'   experimental fraction with proteoQ removals of unclassified entries.
#' @param col_fraction A column key of sample fractions, e.g., CP, NP and ChA.
#' @param bandle_pred_col A column key of Bandle subcellular prediction.
confu_bandle_by_frac <- function (
    df_inner, id = "gene", col_fraction = "Fraction", 
    bandle_pred_col = "pred_bandle") 
{
  df_fracs <- split(df_inner, df_inner[[bandle_pred_col]])
  
  purrr::imap(df_fracs, function(df_frac, frac_name) {
    df_frac <- df_frac |>
      dplyr::filter(.data[[col_fraction]] != frac_name)
  }) |>
    dplyr::bind_rows()
}


#' Helper remove tethered mitochondrial, ER proteins.
#' 
hremove_contams <- function (df, col_fraction = "Fraction", 
                             tether_fraction = "NP", 
                             # poly_x, poly_y, 
                             polygons)
{
  rows <- df[[col_fraction]] == tether_fraction
  
  dfx <- remove_contams(
    df[rows, ], poly_x = polygons[["x"]], poly_y = polygons[["y"]])
  tethers <- attr(dfx, "tethers")
  
  df <- dplyr::bind_rows(df[!rows, ], dfx)
  attr(df, "tethers") <- tethers
  
  df
}


#' Remove tethered mitochondrial, ER proteins.
#' 
remove_contams <- function (df, id = "gene", poly_x, poly_y)
{
  if (is.null(poly_x) || is.null(poly_y)) {
    return(df)
  }
  
  inside_status <- sp::point.in.polygon(
    point.x = df$UMAP1, 
    point.y = df$UMAP2, 
    pol.x   = poly_x, 
    pol.y   = poly_y
  )
  
  ids_inside <- unique(df[[id]][inside_status > 0])
  
  rows <- df[[id]] %in% ids_inside
  dfx  <- df[rows, ]
  df <- df[!rows, ]
  
  if (FALSE) {
    readr::write_tsv(dfx, "~/cp_np_interfacial.tsv")
  }
  
  attr(df, "tethers") <- ids_inside
  
  df
}


#' Plot UMAP of markers with proteoQ annotation.
#' 
plot_marks_by_pq <- function (
    df, df_markers, id = "gene", sub_dir, col_fraction = "Fraction", 
    col_subtype = "Type", key_col = "log2FC", levels_ref_compartments = NULL, 
    my_colors = NULL, remove_tethers = FALSE, 
    poly_x_cp_np, poly_y_cp_np, poly_x_np_cha, poly_y_np_cha, 
    polygons)
{
  df <- df |>
    dplyr::left_join(df_markers, by = id) |>
    dplyr::filter(!is.na(location))
  
  if (remove_tethers) {
    for (fr in levels_fraction) {
      df <- hremove_contams(
        df, col_fraction = col_fraction, tether_fraction = fr, 
        polygons = polygons[[fr]])
    }
  }
  
  ## Combined plot
  n_subtype <- length(unique(df[[col_fraction]]))
  
  p <- ggplot() +
    geom_point(data = df, aes(x = UMAP1, y = UMAP2, color = location), 
               shape = 16, size = 2, alpha = 1) +
    scale_color_manual(values = my_colors, drop = FALSE) + 
    scale_x_continuous(limits = c(-10, 10)) +
    scale_y_continuous(limits = c(-10, 12)) +
    theme_minimal(base_size = 14) +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Location"
    ) +
    theme(
      strip.text = element_text(
        face = "bold",
        size = 14,
        color = "black"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"), 
      legend.title = element_blank()
    ) + 
    guides(
      color = guide_legend(override.aes = list(size = 3, alpha = 1))
    ) + 
    facet_wrap(
      vars(!!rlang::sym(col_fraction)), 
      ncol = n_subtype, 
      labeller = labeller(.default = label_value, .multi_line = FALSE))
  
  if (draw_line <- FALSE) {
    if (FALSE) {
      df_line <- data.frame(
        x = poly_x_cp_np[[1]], y = poly_y_cp_np[[1]], 
        xend = poly_x_cp_np[[2]], yend = poly_y_cp_np[[2]])
    }

    k <- 0.08685714; b <- 1.148571; x <- -2; xend <- 10;
    df_line <- data.frame(
      x = x, y = k * x + b, 
      xend = xend, yend = k * xend + b)
    df_line <- dplyr::bind_rows(
      df_line |> dplyr::mutate(!!col_fraction := "CP"), 
      df_line |> dplyr::mutate(!!col_fraction := "NP")
    )
    df_line <- df_line |>
      dplyr::mutate(
        !!col_fraction := factor(!!rlang::sym(col_fraction), levels = levels_fraction))
    
    p <- p + 
      # geom_hline(yintercept = 1.3, linetype = "longdash", linewidth = .5) +
      geom_segment(
        df_line,
        mapping = aes(x = x, y = y, xend = xend, yend = yend), 
        linetype = "longdash",
        linewidth = 0.5)
    
    # k <- 10/3; b <- -28/3; 1.5 * k + b
    df_line <- data.frame(
      x = poly_x_np_cha[[1]], y = poly_y_np_cha[[1]], 
      xend = poly_x_np_cha[[2]], yend = poly_y_np_cha[[2]])
    df_line <- dplyr::bind_rows(
      df_line |> dplyr::mutate(!!col_fraction := "NP"), 
      df_line |> dplyr::mutate(!!col_fraction := "ChA")
    )
    df_line <- df_line |>
      dplyr::mutate(
        !!col_fraction := factor(!!rlang::sym(col_fraction), levels = levels_fraction))
    
    p + 
      geom_segment(
        df_line,
        mapping = aes(x = x, y = y, xend = xend, yend = yend), 
        linetype = "longdash",
        # color = "darkgrey",
        linewidth = 0.5)
  }
  
  ggsave(
    file.path(
      dat_dir, "Protein", "Bandle", 
      paste0(sub_dir, "_", key_col, "_markers", ".png")), 
    width = 12, height = 6, dpi = 1200)
}


#' Plot Bandle.
#' 
#' @param df_sc A subcellular database.
#' @param work_dir A parent working directory. 
#' @param sub_grp A sub group.
#' @param dat_dir A working directory.
#' @param filename_ref Reference file name of UniProt and GO.
#' @param filename_hpa Reference file name of HPA.
#' @import ggplot2
plot_bandle <- function (id = "gene", dat_dir, dat_dir_pqsvm = NULL, 
                         sub_dir, db_path, df = NULL, 
                         df_sc = NULL, subgrp_prefix = "", 
                         filename_ref = NULL, filename_hpa = NULL, 
                         file_umap = NULL, 
                         FUN = umap::umap, 
                         key_up_go = "location_slim3", 
                         df_markers = NULL, 
                         key_col = c("log2FC", "log10Int"), 
                         pred_type = "bandle", col_subtype = "Type",  
                         col_fraction = "Fraction", 
                         levels_subtype = NULL, levels_group = NULL, 
                         ext = "png", ncol = 3, remove_tethers = FALSE, 
                         use_cache = TRUE, is_hpa = FALSE, 
                         poly_x_cp_np = c(-10, 7.5, 7.5, -10), poly_y_cp_np = c(.28, 1.8, 3.0, 10), 
                         poly_x_np_cha = c(1.5, -.2, 10, 10), poly_y_np_cha = c(-4.333, -10, -10, 10), 
                         
                         polygons = list(
                           Veh = list(
                             # CP = list(x = c(-5, 0, 0, -5), y = c(-5, -5, -2, -2)), 
                             CP = NULL, 
                             # NP = list(x = c(-10, 8, 8, -10), y = c(0, 0, 12, 12)), 
                             NP = list(x = c(5.5, 8.5, .5, -3.5), y = c(.5, 2.5, 11, 11)), 
                             ChA = NULL # ChA = list(x = c(1.5, -.2, 10, 10, 0), y = c(-5, -10, -10, 2, 2))
                           ),
                           TG = list(
                             CP = NULL, 
                             # CP = list(x = c(-7.5, 2.5, 2.5, -7.5), y = c(-7.5, -7.5, -5, 0)), 
                             # NP = list(x = c(-10, 7.5, 7.5, -10), y = c(.28, 1.8, 10, 10)), 
                             NP = list(x = c(1, 6, 2, -4.5, -4.5), y = c(.5, 3, 6, 6, 4)), 
                             ChA = NULL # ChA = list(x = c(1.5, -.2, 10, 10, 5), y = c(-4.333, -10, -10, 10, 10))
                           )
                         ), 
                         
                         width_p = NULL, height_p = NULL, 
                         width_px = 3, height_px = 2)
{
  if (FALSE) {
    polygons = list(
      Veh = list(
        CP = list(x = c(-5, -2, -2, -5), y = c(-5, -5, -2, -2)), 
        NP = list(x = c(-5, 8, 8, -5), y = c(0, 0, 11.2, 11.2)), 
        ChA = list(x = c(1.5, -.2, 10, 10, 0), y = c(-5, -10, -10, 2, 2))),
      TG = list(
        CP = list(x = c(-7.5, 2.5, 2.5, -7.5), y = c(-7.5, -7.5, -5, 0)), 
        NP = list(x = c(-10, 7.5, 7.5, -10), y = c(.28, 1.8, 4, 10)), 
        ChA = list(x = c(1.5, -.2, 10, 10, 5), y = c(-4.333, -10, -10, 10, 10)))
    )
  }

  key_col <- match.arg(key_col)
  col_x <- "Group"
  
  if (is.null(df)) {
    # paste0(subgrp_prefix, find_trend_df2()[[1]]))
    df <- readr::read_tsv(
      file.path(dat_dir, "Protein", "Trend", proteoQ:::find_trend_df2()[[1]])) |>
      dplyr::filter(.data[[col_x]] %in% levels_group)

    df <- df |>
      dplyr::mutate(
        !!col_x := factor(!!rlang::sym(col_x), levels = levels_group), 
        !!col_subtype := factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
        !!col_fraction := factor(!!rlang::sym(col_fraction), levels = levels_fraction), ) |>
      dplyr::arrange_at(c(col_x, col_fraction, col_subtype))
  }

  # Use log2FC
  label_scheme_sub <- load_ls_group(dat_dir = dat_dir) |>
    dplyr::filter(.data[[col_subtype]] %in% levels_subtype)
  
  polygons <- polygons[[levels_subtype]]
  
  if (is.null(levels_group)) {
    levels_group <- label_scheme_sub |>
      dplyr::arrange_at("Order") |>
      dplyr::pull("Group") |>
      unique()
  }

  if (key_col == "log10Int") {
    df_expr <- hmake_subcell_expr(
      df = df, id = id, label_scheme_sub = label_scheme_sub, 
      key_col = key_col, col_subtype = col_subtype, 
      levels_subtype = levels_subtype, levels_fraction = levels_fraction, 
      match_colnames_to_levels_fraction = FALSE)
    df_expr <- as.data.frame(do.call(cbind, df_expr))
    df_expr <- tibble::rownames_to_column(df_expr[, levels_group], id)
  } else {
    df_expr <- df |> 
      dplyr::select(dplyr::all_of(c(id, col_x, key_col))) |> 
      tidyr::pivot_wider(names_from = col_x, values_from = key_col)
  }
  
  if (is.null(file_umap)) {
    file_umap <- file.path(dat_dir, "Protein", "Trend", 
                           paste0(key_col, "_", levels_subtype, "_bandle_umap.tsv"))
  }

  if (use_cache && file.exists(file_umap)) {
    df_umap <- readr::read_tsv(file_umap)
    df_umap[[col_subtype]] <- NULL
  } else {
    df_umap <- local({
      set.seed(42)
      # res_umap <- umap::umap(df_expr[, grepl("_CP$|_NP$|ChA$", names(df_expr))])
      res_umap <- FUN(df_expr[, grepl("_CP$|_NP$|ChA$", names(df_expr))])
      df_umap <- res_umap$layout[, 1:2]
      df_umap <- round(df_umap, digits = 2)
      colnames(df_umap) <- paste0("UMAP", 1:2)
      rownames(df_umap) <- df_expr$gene
      df_umap <- tibble::rownames_to_column(data.frame(df_umap), "gene")
      
      # readr::write_tsv(df_umap, file_umap)
      
      df_umap
    })
  }

  pred_col <- if (pred_type == "bandle") "tagm.map.allocation" else "svm"
  pval_col <- if (pred_type == "bandle") "tagm.map.probability" else "svm.scores"
  bandle_pred_col <- "pred_bandle"
  bandle_pval_col <- "pval_bandle"
  
  df_bandle <- lapply(levels_subtype, function (x) {
    readr::read_tsv(
      file.path(dat_dir, "Protein", "Bandle", 
                paste0(pred_type, "_pred_", sub_dir, "_", x, ".tsv"))) |>
      dplyr::rename(!!id := protein_id) |>
      dplyr::mutate(!!col_subtype := x)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      !!bandle_pred_col := 
        factor(!!rlang::sym(pred_col), levels = levels_fraction), ) |>
    dplyr::rename(!!bandle_pval_col := !!rlang::sym(pval_col)) |>
    dplyr::select(
      dplyr::one_of(c(id, col_subtype, bandle_pred_col, bandle_pval_col)))
  
  df_bandle <- df_umap |> 
    dplyr::left_join(df_bandle, by = id) |>
    dplyr::filter(!is.na(.data[[col_subtype]]))
  
  df_pq <- df[, c(id, col_x, col_fraction, col_subtype, "sub_location")] |>
    dplyr::filter(!is.na(.data[["sub_location"]]))

  df_pq <- df_umap |> dplyr::left_join(df_pq, by = id) |>
    dplyr::filter(!is.na(.data[["sub_location"]]))
  
  file_pq_svm <- file.path(dat_dir_pqsvm, "Protein", "Trend", "Protein_svmProb_Z_nclust16.tsv")
  # file_pq_svm <- file.path(dat_dir_pqsvm, "Protein", "Trend", "umap_Veh.tsv")

  if (isTRUE(file.exists(file_pq_svm))) {
    df_pq_svm <- local({
      dat_dir_pqsvm <- paste0(dat_dir, '_svm')
      
      df_pq_svm <- file_pq_svm |>
        readr::read_tsv() |>
        dplyr::filter(.data[[col_subtype]] %in% levels_subtype) # |>
        # dplyr::rename(!!col_fraction := "compartment")
      
      df_pq_svm <- df_umap |> dplyr::left_join(df_pq_svm, by = id) |>
        dplyr::filter(!is.na(.data[["sub_location"]]))
    })
  } else {
    df_pq_svm <- NULL
  }

  # Marker proteins
  levels_ref_compartments <- local({
    univ <- unique(df_markers[["location"]])
    levs <- names(sort(table(df_markers[["location"]]), decreasing = TRUE))
    levs <- c(levs, univ[!univ %in% levs], "Unknown")
  })
  
  levels_ref_compartments <- c(
    "Cytosol", "ER", "Mitochondrion", "Nucleus", "Chromatin", "Lysosome", 
    "Plasma_membrane", "Proteasome", "Golgi", "Peroxisome", "Endosome", 
    "Ribosome", "Centrosome", "P_body", "Stress_granule"
  )
  
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
  
  df_markers <- df_markers |>
    dplyr::mutate(
      location = factor(location, levels = levels_ref_compartments))
  
  if (remove_tethers) {
    for (fr in levels_fraction) {
      df_pq <- hremove_contams(
        df_pq, col_fraction = col_fraction, tether_fraction = fr, 
        polygons = polygons[[fr]])
    }
  }
  
  ## Confusion matrix
  ans_confu <- confu_bandle(
    df_pq = df_pq, df_bandle = df_bandle, id = id, col_subtype = col_subtype, 
    col_fraction = col_fraction, bandle_pred_col = bandle_pred_col, 
    remove_tethers = remove_tethers, 
    poly_x_cp_np = poly_x_cp_np, poly_y_cp_np = poly_y_cp_np, 
    poly_x_np_cha = poly_x_np_cha, poly_y_np_cha = poly_y_np_cha)
  mat_confu <- attr(ans_confu, "mat_confu")
  names(mat_confu) <- names(ans_confu)
  ans_confu <- do.call(rbind, ans_confu)
  
  local({
    nms_mat <- names(mat_confu)
    
    for (i in seq_along(mat_confu)) {
      mat <- matrix(mat_confu[[i]], nrow = 3, byrow = FALSE, 
                    dimnames = list(Experimental = c("CP", "NP", "ChA"),
                                    Predicted = c("CP", "NP", "ChA")))
      
      df_mat <- as.data.frame(as.table(mat))
      
      ggplot(df_mat, aes(x = Experimental, y = Predicted, fill = Freq)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = Freq), color = "black", size = 5) +
        scale_fill_gradient(low = "#e0f2fe", high = "#0284c7") +
        theme_minimal() +
        labs(title = NULL,
             x = "Experimental",
             y = "Predicted",
             fill = "Count") +
        theme(panel.grid = element_blank(), 
              legend.position = "none", )
      for (ext in c("png", "pdf")) {
        ggsave(file.path(dat_dir, "Protein", "Bandle", 
                         paste0(sub_dir, "_", key_col, "_confusion_", nms_mat[[i]], ".", ext)), 
               width = 3, height = 3, dpi = 600)
      }
    }
  })

  ## Plot markers (has its own `remove_contams`)
  plot_marks_by_pq(
    df = df_pq, df_markers = df_markers, id = id, sub_dir = sub_dir, 
    col_fraction = col_fraction, col_subtype = col_subtype, key_col = key_col, 
    levels_ref_compartments = levels_ref_compartments, my_colors = my_colors, 
    remove_tethers = remove_tethers, 
    poly_x_cp_np = poly_x_cp_np, poly_y_cp_np = poly_y_cp_np, 
    poly_x_np_cha = poly_x_np_cha, poly_y_np_cha = poly_y_np_cha, 
    polygons = polygons)
  
  ## Map to the reference of UniProt + GO
  ok_ref <- (!is.null(filename_ref)) && file.exists(filename_ref)
  
  df_ref <- prep_ref_db(
    filename_ref = filename_ref, id = id, df = df, df_pq = df_pq, 
    key_up_go = key_up_go, col_subtype = col_subtype, 
    levels_subtype = levels_subtype, is_hpa = FALSE)
  
  if (is.null(filename_hpa)) {
    df_hpa <- NULL
  } else {
    df_hpa <- prep_ref_db(
      filename_ref = filename_hpa, id = id, df = df, df_pq = df_pq, 
      key_up_go = key_up_go, col_subtype = col_subtype, 
      levels_subtype = levels_subtype, is_hpa = TRUE)
  }
    
  if (ok_ref) {
    df_plots <- list(df_ref, df_hpa, df_pq, df_pq_svm, df_bandle, ans_confu)
    key_col_pred <- c("location_slim3", "location_slim3", "sub_location", 
                      "sub_location", bandle_pred_col, bandle_pred_col)
    out_nm_type  <- c("Ref", "HPA", "proteoQ", "SVM", pred_type, "confusion")
  } else {
    df_plots <- list(df_pq, df_bandle, ans_confu)
    key_col_pred <- c("sub_location", bandle_pred_col, bandle_pred_col)
    out_nm_type <- c("proteoQ", pred_type, "confusion")
  }
  
  df_plot_all <- vector("list", length(df_plots))

  for (i in seq_along(df_plots)) {
    df_plot <- df_plots[[i]]
    
    if (is.null(df_plot) || !nrow(df_plot)) {
      next
    }
    
    out_nm_i <- out_nm_type[[i]]
    key_col_pred_i <- key_col_pred[[i]]
    
    df_plot <- df_plot |>
      dplyr::filter(!is.na(.data[[col_subtype]]))
    
    n_subtype <- length(unique(df_plot[[col_subtype]]))
    
    df_plot <- df_plot |> 
      dplyr::mutate(!!col_subtype := factor(!!rlang::sym(col_subtype), 
                                            levels = levels_subtype))

    if (length(df_plot[[col_fraction]])) {
      n_frac <- length(unique(df_plot[[col_fraction]]))
    } else {
      n_frac <- length(unique(df_plot[[key_col_pred_i]]))
    }

    if (n_frac < 1L) n_frac <- 1L
    
    df_plot_all[[i]] <- df_plot |> 
      dplyr::mutate(Source = out_nm_i, 
                    Location_ = .data[[key_col_pred_i]])
    
    df_plot <- df_plot |> 
      dplyr::mutate(!!rlang::sym(key_col_pred_i) := 
                      factor(!!rlang::sym(key_col_pred_i), 
                             levels = levels_fraction))

    p <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, 
                                   color = !!rlang::sym(key_col_pred_i))) +
      geom_point(shape = 46, alpha = .5, size = .02) +
      scale_color_manual(
        values = c("CP" = "#E41A1C", "NP" = "#377EB8", "ChA" = "#4DAF4A"),
        breaks = c("CP", "NP", "ChA"), 
        drop = FALSE
      ) +
      guides(color = guide_legend(
        override.aes = list(size = 3, shape = 16, alpha = 1))) +
      theme_bw() +
      labs(
        title = NULL,
        x = "UMAP 1",
        y = "UMAP 2",
        color = "Subcellular\nLocation"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.title = element_blank(),
        panel.grid = element_blank()
      )
    
    if (n_subtype == 1L) {
      px <- p
    } else {
      px <- p + 
        facet_wrap(
          vars(!!rlang::sym(col_subtype)), 
          ncol = n_subtype, 
          labeller = labeller(.default = label_value, .multi_line = FALSE))
    }

    # No Fraction info from Bandle results; use confusion plot instead
    if (!out_nm_i %in% c("bandle", "Ref", "HPA")) {
      p <- p + 
        facet_wrap(
          vars(!!rlang::sym(col_subtype), !!rlang::sym(col_fraction)), 
          ncol = n_frac, 
          labeller = labeller(.default = label_value, .multi_line = FALSE))
    } else {
      p <- p + 
        facet_wrap(
          vars(!!rlang::sym(col_subtype), !!rlang::sym(key_col_pred_i)), 
          ncol = n_frac, 
          labeller = labeller(.default = label_value, .multi_line = FALSE))
    }
    
    # PNG, wrapped by subtype only
    out_path <- file.path(dat_dir, "Protein", "Bandle")
    
    if (FALSE) {
      # poly_x_cp_np <- c(-4, 7.5, 7.5, -4.2); poly_y_cp_np <- c(4.7, -3, 3, 7.2)
      # poly_x_cp_np = c(-4, 7.5, 7.5, -4.6); poly_y_cp_np = c(.8, 1.8, 3, 7.2)
      poly_x_cp_np <- c(-10, 7.5, 7.5, -10); poly_y_cp_np <- c(.28, 1.8, 3.0, 10)
      
      px + 
        geom_polygon(
          data.frame(x = poly_x_cp_np, y = poly_y_cp_np), mapping = aes(x = x, y = y), 
          fill = NA, color = "black", linewidth = .2)
      ggsave(file.path(out_path, 
                       paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, 
                              "_combined", ".", "png")), 
             width = width_px, height = height_px, 
             dpi = 1200)
    }
    
    n_grids <- n_frac * n_subtype
    
    if (is.null(width_px)) {
      width_px <- n_subtype * 2 + 1
    }
    if (is.null(height_px)) {
      height_px <- 2 + .333
    }
    if (is.null(width_p)) {
      width_p <- n_frac * 2 + 1
    }
    if (is.null(height_p)) {
      height_p <- n_subtype * 2 + .4
    }
    
    ggsave(plot = px, 
           file.path(out_path, 
                     paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, 
                            "_combined", ".", "png")), 
           width = width_px, height = height_px, 
           dpi = 1200)
    
    ggsave(plot = p, 
           file.path(
             out_path, 
             paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, ".", "png")), 
           width = width_p, height = height_p, 
           dpi = 1200)
    
    # pdf
    ggsave(plot = px, 
           file.path(out_path, 
                     paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, 
                            "_combined", ".", "pdf")), 
           # width = width_px * 1.25, height = height_px * 1.25, 
           width = width_px, height = height_px, 
           dpi = 1200)
    
    ggsave(plot = p, 
           file.path(
             out_path, 
             paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, ".", "pdf")), 
           # width = width_p, height = height_p, 
           width = width_p, height = height_p, 
           dpi = 1200)
  }
  
  if (FALSE) {
    for (i in seq_along(df_plot_all)) {
      df_plot_all[[i]] <- df_plot_all[[i]] |>
        dplyr::select(dplyr::one_of(c("gene", "UMAP1", "UMAP2", "Type", "Source", 
                                      "Location_")))
    }
  }

  df_plot_all <- df_plot_all |>
    dplyr::bind_rows()
  
  readr::write_tsv(df_plot_all, file.path(out_path, "results_combined.tsv"))
  
  df_plot_all <- df_plot_all |>
    dplyr::filter(.data[["Source"]] %in% c("proteoQ", "Ref", "HPA", "bandle")) |>
    dplyr::mutate(Source = dplyr::recode(
      Source,
      "proteoQ" = "This work",
      "Ref"     = "UP_GO",
      "HPA"     = "HPA",
      "bandle"  = "BANDLE"
    )) |>
    dplyr::mutate(Source = factor(Source, levels = c("This work", "UP_GO", "HPA", "BANDLE")))
  
  # df_plot_all <- df_plot_all |> dplyr::filter(Source != "UP_GO")
  
  p_all <- ggplot(df_plot_all, aes(x = UMAP1, y = UMAP2, 
                                   color = Location_)) +
    geom_point(shape = 16, alpha = .5, size = .25) +
    scale_color_manual(
      values = c("CP" = "#E41A1C", "NP" = "#377EB8", "ChA" = "#4DAF4A"),
      breaks = c("CP", "NP", "ChA"), 
      drop = FALSE
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 3, shape = 16, alpha = 1))) +
    theme_bw() +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Subcellular\nLocation"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
  
  p_all <- p_all + 
    facet_wrap(
      vars(Source), 
      ncol = 2, 
      labeller = labeller(.default = label_value, .multi_line = FALSE))
  
  ggsave(plot = p_all, 
         file.path(
           out_path, 
           "umap_all_sources.png"), 
         # width = width_p, height = height_p, 
         width = 6, height = 5, 
         dpi = 1200)
}


#' Plot combined data
#' 
plot_combined <- function (
    dat_dir, filepath = file.path(dat_dir, "Protein", "Bandle"), 
    filename = "results_combined.tsv", out_name = "umap_all_sources.png", ...)
{
  dots <- rlang::enexprs(...)
  lang_dots     <- dots[unlist(lapply(dots, is.language))]
  filter2_dots  <- lang_dots[grepl("^filter2_", names(lang_dots))]
  arrange2_dots <- lang_dots[grepl("^arrange2_", names(lang_dots))]
  dots <- dots[!dots %in% c(filter2_dots, arrange2_dots)]
  
  df <- readr::read_tsv(file.path(filepath, filename))
  
  df <- df |>
    proteoQ:::filters_in_call(!!!filter2_dots) |>
    proteoQ:::arrangers_in_call(!!!arrange2_dots)

  df <- df |>
    dplyr::mutate(
      Source = dplyr::case_when(
        Source == "Ref" ~ "UP_GO", 
        Source == "proteoQ" ~ "This work", 
        Source == "bandle" ~ "BANDLE", 
        Source == "HPA" ~ "HPA", 
        TRUE ~ "Other"
      )
    )
  
  df <- df |> dplyr::filter(!Source %in% "Other")
  
  df <- local({
    dfs <- split(df, df$Source)
    
    dfs[["This work"]][["location_"]] <- dfs[["This work"]][["sub_location"]]
    dfs[["BANDLE"]][["location_"]] <- dfs[["BANDLE"]][["pred_bandle"]]
    dfs[["UP_GO"]][["location_"]] <- dfs[["UP_GO"]][["location_slim3"]]
    dfs[["HPA"]][["location_"]] <- dfs[["HPA"]][["location_slim3"]]
    
    df <- dfs |> dplyr::bind_rows()
    
    df <- df |> dplyr::mutate(
      Source = factor(Source, levels = c("This work", "BANDLE", "HPA", "UP_GO")),
      location_ = factor(location_, levels = c("CP", "NP", "ChA")), )
  })
  
  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = location_)) +
    geom_point(shape = 46, alpha = .5, size = .02) +
    scale_color_manual(
      values = c("CP" = "#E41A1C", "NP" = "#377EB8", "ChA" = "#4DAF4A"),
      breaks = c("CP", "NP", "ChA")
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 3, shape = 16, alpha = 1))) +
    theme_bw() +
    labs(
      title = NULL,
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Subcellular\nLocation"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_blank(),
      panel.grid = element_blank()
    ) + 
    facet_wrap(~ Source, ncol = 2, labeller = label_value)
  
  ggsave(plot = p, 
         file.path(filepath, "umap_all.pdf"), width = 7, height = 5.8, dpi = 1200)
  
  ggsave(plot = p, 
         file.path(filepath, "umap_all.png"), width = 7, height = 5.8, dpi = 1200)
}



#' Prepare reference data.
#' 
#' @param filename_ref Reference file name of UniProt and GO, or HPA.
#' @param id An identifier.
#' @param df A data frame for subsetting of genes.
#' @param df_pq A data frame of proteoQ results.
#' @param key_up_go The column key of subcellular annotations.
#' @param col_subtype The column of sample sub types.
#' @param levels_subtype The levels of subtypes.
#' @param is_hpa Logical; is the data type HPA.
#' @param use_intercept Logical; use the intercept of UniProt and GO or not.
prep_ref_db <- function (
    filename_ref, id = "gene", df = NULL, df_pq = NULL, 
    key_up_go = "location_slim3", col_subtype = "Type", levels_subtype = NULL, 
    is_hpa = FALSE, use_intersect = TRUE)
{
  df_ref <- readr::read_tsv(file.path(filename_ref))
  
  if ("HPA" %in% colnames(df_ref)) {
    if (is_hpa) {
      df_ref <- df_ref |> dplyr::filter(source == "HPA")
    } else {
      df_ref <- df_ref |> dplyr::filter(source != "HPA")
      
      if (use_intersect) {
        df_ref <- local({
          dfs <- split(df_ref, df_ref$source)
          gns <- lapply(dfs, `[[`, "gene")
          gns <- purrr::reduce(gns, intersect)
          
          df_ref <- df_ref |> 
            dplyr::filter(gene %in% gns)
        })
      }
    }
  }
  
  df_ref <- df_ref[, c("gene", key_up_go)] |>
    unique() |>
    dplyr::filter(gene %in% unique(df$gene))
  
  df_ref <- df_ref |> 
    dplyr::filter(!is.na(.data[[key_up_go]])) |>
    dplyr::left_join(
      df_pq |> dplyr::select(-dplyr::one_of("sub_location")), by = id)
  
  df_ref <- df_ref |>
    dplyr::select(dplyr::one_of(
      c(id, "UMAP1", "UMAP2", key_up_go, col_subtype))) |>
    unique() |>
    dplyr::mutate(
      !!col_subtype := 
        factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
      !!key_up_go := 
        factor(!!rlang::sym(key_up_go), levels = levels_fraction),) |> 
    dplyr::filter(!is.na(!!rlang::sym(col_subtype)))
  
  df_ref
}


#' Prepare HPA data.
#' 
prep_hpa <- function ()
{
  levels_ref_compartments <- c(
    "Cytosol", "ER", "Mitochondrion", "Nucleus", "Chromatin", "Lysosome", 
    "Plasma_membrane", "Proteasome", "Golgi", "Peroxisome", "Endosome", 
    "Ribosome", "Centrosome", "P_body", "Stress_granule"
  )
  
  fraction_map <- c(
    "Chromatin"       = "ChA",
    "Centrosome"      = "NP",
    "Nucleus"         = "NP",
    "Cytosol"         = "CP",
    "Endosome"        = "CP",
    "ER"              = "CP",
    "Golgi"           = "CP",
    "Lysosome"        = "CP",
    "Mitochondrion"   = "CP",
    "P_body"          = "CP",
    "Peroxisome"      = "CP",
    "Plasma_membrane" = "CP",
    "Proteasome"      = "CP",
    "Ribosome"        = "CP",
    "Stress_granule"  = "CP"
  )
  
  db_path <- "D:/QZ/dbs/cc"
  species <- "human"
  
  file_hpa <- "hpa_subcellular_location.tsv"
  
  df_hpa <- readr::read_tsv(file.path(db_path, file_hpa)) |>
    dplyr::rename(gene = "Gene name", 
                  reliability = "Reliability", 
                  location = "Main location")
  
  df_hpa <- df_hpa |>
    dplyr::filter(reliability %in% c("Enhanced", "Approved", "Supported")) |>
    dplyr::select(c("gene", "reliability", "location"))
  
  df_hpa <- df_hpa |>
    tidyr::separate_rows(location, sep = ";") |>
    dplyr::mutate(location = trimws(location)) |> 
    dplyr::filter(!is.na(location))
  
  # unique(df_hpa$location)
  df_hpa <- df_hpa |>
    dplyr:: mutate(
      ref_compartment = dplyr::case_when(
        location %in% c("Cytosol", "Cytoplasmic bodies") ~ "Cytosol",
        location == "Endoplasmic reticulum" ~ "ER",
        location == "Mitochondria" ~ "Mitochondrion",
        location %in% c("Nucleoplasm", "Nuclear speckles", "Nuclear membrane", 
                        "Nuclear bodies", "Nucleoli", "Nucleoli fibrillar center", 
                        "Nucleoli rim") ~ "Nucleus",
        location == "Mitotic chromosome" ~ "Chromatin",
        location == "Lysosomes" ~ "Lysosome",
        location %in% c("Plasma membrane", "Cell Junctions", "Focal adhesion sites") ~ "Plasma_membrane",
        location == "Golgi apparatus" ~ "Golgi",
        location == "Peroxisomes" ~ "Peroxisome",
        location %in% c("Endosomes", "Vesicles") ~ "Endosome",
        location %in% c("Centrosome", "Flagellar centriole", "Centriolar satellite") ~ "Centrosome",
        location == "Aggresome" ~ "Stress_granule",
        
        # Unmapped to a specific ref category? Map them to Cytosol or leave as NA.
        location %in% c("Actin filaments", "Microtubules", "Intermediate filaments", 
                        "Microtubule ends", "Mitotic spindle", "Cytokinetic bridge", 
                        "Midbody", "Midbody ring", "Cleavage furrow", "Kinetochore", 
                        "Rods & Rings", "Lipid droplets", "Basal body", "Primary cilium", 
                        "Primary cilium transition zone", "Primary cilium tip", 
                        "Acrosome", "End piece", "Principal piece", "Equatorial segment", 
                        "Perinuclear theca", "Mid piece", "Connecting piece", 
                        "Annulus", "Calyx") ~ NA_character_, # "Cytosol"
        TRUE ~ NA_character_
      ),
      ref_compartment = factor(ref_compartment, levels = levels_ref_compartments)
    )
  
  df_hpa <- df_hpa |>
    dplyr::rename(location_fine = location, 
                  location = ref_compartment)
  
  df_hpa <- df_hpa[, c("gene", "location")] |>
    unique()
  
  df_hpa <- df_hpa |>
    dplyr::group_by(gene) |>
    dplyr::filter(dplyr::n() == 1) |>
    dplyr::filter(!is.na(location))
  
  df_hpa <- df_hpa |>
    dplyr::mutate(
      location_slim3 = dplyr::recode(location, !!!fraction_map),
      location_slim3 = factor(location_slim3, levels = c("CP", "NP", "ChA"))
    )
  
  readr::write_tsv(df_hpa, file.path(db_path, "hpa_single_compartment.tsv"))
  
  df_hpa
}


