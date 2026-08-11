#' Prepare Bandle.
#' 
#' @param df_expr_sub An expression data by subtype.
#' @param comartment_2 A set of markers.
prep_bandle <- function (df_expr_sub, compartment_2)
{
  pheno_data <- data.frame(
    FractionNumber = seq_len(ncol(df_expr_sub)),
    row.names = colnames(df_expr_sub)
  )
  
  feature_data <- data.frame(
    protein_id = rownames(df_expr_sub),
    markers = compartment_2,
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
    levels_fraction = c("CP", "NP", "ChA"), levels_subtype = c("Veh", "TG"), 
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
  
  df_expr <- df[, c("gene", "group", "log10Int")] |>
    dplyr::mutate(Int = 10^log10Int) |>
    dplyr::select(-"log10Int") |> 
    pivot_wider(names_from = group, values_from = Int, values_fill = 0) |>
    dplyr::left_join(df_sc, by = "gene") |>
    dplyr::mutate(location = ifelse(is.na(location), "unknown", location),
                  compartment_2 = ifelse(is.na(compartment_2), 
                                         "unknown", compartment_2),)
  
  # Convert to data matrix
  genes <- df_expr$gene
  locations <- df_expr$location
  compartment_2 <- df_expr$compartment_2
  df_expr <- as.matrix(
    df_expr[, -which(names(df_expr) %in% c("gene", "location", "compartment_2"))])
  rownames(df_expr) <- genes
  df_expr <- df_expr / rowSums(df_expr)
  
  # Can have NaN if only select a subset of input data
  # df_expr[is.na(df_expr)] <- 0.0
  oks <- rowSums(is.na(df_expr)) == 0L
  df_expr <- df_expr[oks, ]
  genes <- genes[oks]
  locations <- locations[oks]
  compartment_2 <- compartment_2[oks]
  
  # Separate into sub-groups
  n_subtypes <- length(levels_subtype)
  df_exprs <- vector("list", n_subtypes)
  for (i in seq_along(levels_subtype)) {
    df_exprs[[i]] <- 
      df_expr[, grepl(paste0("^", levels_subtype[[i]], "_"), colnames(df_expr))]
  }
  names(df_exprs) <- levels_subtype
  rm(list = "df_expr")
  
  markers <- unique(df_sc$compartment_2)
  bandle_inputs <- lapply(df_exprs, prep_bandle, compartment_2)
  
  res <- mapply(calc_bandle, bandle_inputs, names(bandle_inputs), 
                MoreArgs = list(dat_dir = dat_dir))
  
  res <- if (use_svm) res[["svm"]] else res[["bandle"]]
  
  plot_bandle(dat_dir = dat_dir, sub_dir = sub_dir, db_path = db_path, 
              df = df, df_sc = df_sc, subgrp_prefix = "", 
              pred_type = "bandle", col_subtype = "Type",  
              col_fraction = "Fraction", 
              levels_subtype = levels_subtype)
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
    col_fraction = "Fraction", bandle_pred_col = "pred_bandle") 
{
  df_pqs <- split(df_pq, df_pq[[col_subtype]])
  df_banles <- split(df_bandle, df_bandle[[col_subtype]])
  
  df_inners <- df_pq[, c(id, col_fraction, col_subtype)] |>
    dplyr::inner_join(df_bandle, by = c(id, col_subtype)) %>%
    split(.[[col_subtype]])
  
  ans_confu <- lapply(
    df_inners, confu_bandle_by_frac, col_fraction = col_fraction, 
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
    df_inner, col_fraction = "Fraction", bandle_pred_col = "pred_bandle") 
{
  df_fracs <- split(df_inner, df_inner[[bandle_pred_col]])
  
  purrr::imap(df_fracs, function(df_frac, frac_name) {
    df_frac |>
      dplyr::filter(.data[[col_fraction]] != frac_name)
  }) |>
    dplyr::bind_rows()
}


#' Plot Bandle.
#' 
#' @param df_sc A subcellular database.
#' @param work_dir A parent working directory. 
#' @param sub_grp A sub group.
#' @param dat_dir A working directory.
plot_bandle <- function (id = "gene", dat_dir, sub_dir, db_path, df = NULL, 
                         df_sc = NULL, subgrp_prefix = "", 
                         key_col = c("log2FC", "log10Int"), 
                         pred_type = "bandle", col_subtype = "Type",  
                         col_fraction = "Fraction", 
                         levels_subtype = NULL, levels_group = NULL, 
                         ext = "png")
{
  key_col <- match.arg(key_col)
  col_x <- "group"
  
  if (is.null(df)) {
    # Any one OK
    df <- readr::read_tsv(
      file.path(dat_dir, "Protein", "Trend", 
                paste0(subgrp_prefix, "Protein_Trend_Z_nclust32.txt")))
  }

  # Use log2FC
  if (is.null(levels_group)) {
    levels_group <- load_ls_group(dat_dir = dat_dir) |>
      dplyr::arrange_at("Order") |>
      dplyr::pull("Group") |>
      unique()
  }

  if (key_col == "log10Int") {
    df_expr <- hmake_subcell_expr(
      df = df, id = id, key_col = key_col, col_subtype = col_subtype, 
      levels_subtype = levels_subtype, levels_fraction = levels_fraction, 
      match_colnames_to_levels_fraction = FALSE) %>% 
      do.call(cbind, .) |>
      as.data.frame()
    df_expr <- df_expr[, levels_group] |>
      tibble::rownames_to_column(id)
  } else {
    df_expr <- df[, c(id, col_x, key_col)] |> 
      pivot_wider(names_from = col_x, values_from = key_col)
  }
  
  file_umap <- 
    file.path(dat_dir, "Protein", "Trend", paste0(key_col, "_umap.tsv"))
  
  if (file.exists(file_umap)) {
    umap_df <- readr::read_tsv(file_umap)
  } else {
    umap_df <- local({
      set.seed(42)
      res_umap <- umap::umap(df_expr[, grepl("_CP$|_NP$|ChA$", names(df_expr))])
      umap_df <- res_umap$layout[, 1:2]
      colnames(umap_df) <- paste0("UMAP", 1:2)
      rownames(umap_df) <- df_expr$gene
      umap_df <- tibble::rownames_to_column(data.frame(umap_df), "gene")
      
      readr::write_tsv(
        umap_df, 
        file.path(dat_dir, "Protein", "Trend", paste0(key_col, "_umap.tsv"))
      )
      
      umap_df
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
    dplyr::bind_rows()

  df_bandle <- df_bandle |>
    dplyr::mutate(
      !!bandle_pred_col := factor(!!rlang::sym(pred_col), levels = levels_fraction), ) |>
    dplyr::rename(!!bandle_pval_col := !!rlang::sym(pval_col)) |>
    dplyr::select(dplyr::one_of(c(id, col_subtype, bandle_pred_col, bandle_pval_col)))
  
  df_pq <- df[, c(id, col_x, col_fraction, col_subtype, "sub_location")] |>
    # A cleaner subset by removing unclassified entries
    dplyr::filter(!is.na(.data[["sub_location"]])) |>
    dplyr::mutate(
      !!col_x := factor(!!rlang::sym(col_x), levels = levels_group), 
      !!col_subtype := factor(!!rlang::sym(col_subtype), levels = levels_subtype), 
      !!col_fraction := factor(!!rlang::sym(col_fraction), levels = levels_fraction), ) |>
    dplyr::arrange_at(c(col_x, col_fraction, col_subtype))

  df_bandle <- umap_df |> dplyr::left_join(df_bandle, by = id)
  
  df_pq <- umap_df |> dplyr::left_join(df_pq, by = id) |>
    dplyr::filter(!is.na(.data[["sub_location"]]))
  
  ans_confu <- confu_bandle(
    df_pq = df_pq, df_bandle = df_bandle, id = id, col_subtype = col_subtype, 
    col_fraction = col_fraction, bandle_pred_col = bandle_pred_col)
  mat_confu <- attr(ans_confu, "mat_confu")
  names(mat_confu) <- names(ans_confu)
  ans_confu <- do.call(rbind, ans_confu)
  
  ## Confusiton matrix
  local({
    nms_mat <- names(mat_confu)
    
    for (i in seq_along(mat_confu)) {
      mat <- matrix(mat_confu[[i]], nrow = 3, byrow = FALSE, 
                    dimnames = list(Experimental = c("CP", "NP", "ChA"),
                                    Predicted = c("CP", "NP", "ChA")))
      
      df_mat <- as.data.frame(as.table(mat))
      
      ggplot(df_mat, aes(x = Predicted, y = Experimental, fill = Freq)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = Freq), color = "black", size = 5) +
        scale_fill_gradient(low = "#e0f2fe", high = "#0284c7") +
        theme_minimal() +
        labs(title = NULL,
             x = "Predicted",
             y = "Experimental",
             fill = "Count") +
        theme(panel.grid = element_blank())
      
      for (ext in c("png", "pdf")) {
        ggsave(file.path(dat_dir, "Protein", "Bandle", 
                         paste0(sub_dir, "_", key_col, "_confusion_", nms_mat[[i]], ".", ext)), 
               width = 4, height = 3)
      }
    }
  })

  df_plots <- list(df_pq, df_bandle, ans_confu)
  key_col_pred <- c("sub_location", bandle_pred_col, bandle_pred_col)
  out_nm_type <- c("proteoQ", pred_type, "confusion")
  
  for (i in seq_along(df_plots)) {
    out_nm_i <- out_nm_type[[i]]
    
    p <- ggplot(df_plots[[i]], aes(x = UMAP1, y = UMAP2, 
                                   color = !!rlang::sym(key_col_pred[[i]]))) +
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
      )
    
    px <- p + 
      facet_wrap(
        vars(!!rlang::sym(col_subtype)), 
        ncol = 3, 
        labeller = labeller(.default = label_value, .multi_line = FALSE))
    
    # No Fraction info from Bandle results; use confusion plot instead
    if (out_nm_i != "bandle") {
      p <- p + 
        facet_wrap(
          vars(!!rlang::sym(col_subtype), !!rlang::sym(col_fraction)), ncol = 3, 
          labeller = labeller(.default = label_value, .multi_line = FALSE))
    }

    # png
    ggsave(plot = px, 
           file.path(dat_dir, "Protein", "Bandle", 
                     paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, "_combined", ".", "png")), 
           width = 6, height = 3)
    
    ggsave(plot = p, 
           file.path(
             dat_dir, "Protein", "Bandle", 
             paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, ".", "png")), 
           width = 8, height = 5)
    
    # pdf
    ggsave(plot = px, 
           file.path(dat_dir, "Protein", "Bandle", 
                     paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, "_combined", ".", "pdf")), 
           width = 7.5, height = 3.75, dpi = 600)
    
    ggsave(plot = p, 
           file.path(
             dat_dir, "Protein", "Bandle", 
             paste0(sub_dir, "_", key_col, "_umap_by_", out_nm_i, ".", "pdf")), 
           width = 8, height = 5, dpi = 600)
  }
}


