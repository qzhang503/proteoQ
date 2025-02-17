#' Finds sample IDS for used in the current fitting.
#'
#' @param universe A name list from data table.
#' @param selected A name list from metadata.
#' @param pat A pattern in regular expression.
#'
#' @examples
#' pat <- "^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\("
#'
#' # special character of "+" and pointed brackets
#' universe <- c("N_log2_R000 (<Mix1+Mix2> [base])", "N_log2_R000 (Mix1 [heavy])", 
#'               "N_log2_R000 (Mix2 [heavy])")
#' selected <- c("<Mix1+Mix2> [base]", "Mix2 [heavy]")
#'
#' # special character of "+" and round parenthesis
#' universe <- c("N_log2_R000 ((Mix1+Mix2) [base])", "N_log2_R000 (Mix1 [heavy])", 
#'               "N_log2_R000 (Mix2 [heavy])")
#' selected <- c("(Mix1+Mix2) [base]", "Mix2 [heavy]")
#' 
find_fit_nms <- function(universe, selected, pat = NULL) 
{
  len_u <- length(universe)
  len_s <- length(selected)
  
  if (len_u < len_s) 
    stop("More samples in selected subset than in full set.\n", 
         "  universe = ", paste(universe, collapse = ", "), "\n", 
         "  selected = ", paste(selected, collapse = ", "))
  
  us <- gsub(paste0(pat, "(.*)\\)$"), "\\1", universe)
  ok <- us %in% selected
  us <- us[ok]
  universe <- universe[ok]
  
  ## MaxQuant sample IDs by alphabetic orders and can be different to the 
  #  order of "selected"
  if (!identical(us, selected)) {
    us <- factor(us, levels = selected)
    ord <- order(us)
    us <- us[ord]
    universe <- universe[ord]
  }
  
  invisible(universe)
}


#' Helper of \link{find_fit_nms}.
#' 
#' @param field A field of name prefix.
#' @param pat0 A pattern of regular expression.
#' @param nms A vector of column names with sample IDs in parentheses.
#' @param sub_sids A subset of sample IDs for match with \code{nms}.
hfind_fit_nms <- function (field = "^N_log2_R", nms, sub_sids, 
                           pat0 = "[0-9]{3}[NC]{0,1}\\s+\\(") 
{
  pat <- paste0(field, pat0)
  universe <- nms[grepl(pat, nms)]
  find_fit_nms(universe, sub_sids, pat)
}


#' Calculates standard deviations and the corresponding factors for scaling
#' normalization.
#'
#' The \emph{mean} SD is based on the data from \emph{non-trivial} samples. In
#' the output, \eqn{fct < 1} corresponds to the need of profile broadening in
#' log2Ratios whereas \eqn{fct > 1} corresponds to the need of profile
#' shrinking.
#' 
#' Currently no slice dots for the utility.
#'
#' @param df A data frame.
#' @param label_scheme Experiment summary
#' @inheritParams standPep
#' @return A data frame. Column \code{SD}, the original standard deviations;
#'   \code{fct}, the scaling factors as SD/mean(SD).
calc_sd_fcts <- function (df, range_log2r = c(0, 100), range_int = c(0, 100), 
                          label_scheme) 
{
  sids <- label_scheme$Sample_ID
  
  label_scheme_sd <- label_scheme %>%
    dplyr::filter(!Reference, !grepl("^Empty\\.", Sample_ID))	%>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = (sids)))
  
  non_triv_sids <- label_scheme_sd$Sample_ID
  
  pat <- "^N_log2_R[0-9]{3}[NC]{0,1}|^N_I[0-9]{3}[NC]{0,1}"
  
  SD <- df %>%
    dplyr::select(grep(pat, names(.))) %>%
    dblTrim(range_log2r, range_int) %>%
    `names<-`(replace_NorZ_names(find_NorZ(FALSE), names(.)))
    # `names<-`(gsub("^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\((.*)\\)$", "\\1", names(.)))
  
  non_triv_sds <- SD %>% .[names(.) %in% non_triv_sids]
  coef_sd <- SD/mean(non_triv_sds, na.rm = TRUE)
  
  cbind.data.frame(fct = coef_sd, SD) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>%
    dplyr::arrange(Sample_ID)
}


#' Updates df after normalization.
#' 
#' At a given centering method, the following three take place: 
#' 
#' (1) \code{Z_log2_R} centered and SD adjusted;
#' (2) \code{N_log2_R} centered;
#' (3) \code{N_I} scaled.
#' 
#' @param df A data frame of peptide or protein table.
#' @param label_scheme_fit Experiment summary for samples being selected for
#'   fitting.
#' @param cf_x_fit A data frame containing the \code{x} positions for each
#'   samples indicated in label_scheme_fit.
#' @param sd_coefs_fit The standard deviations for each samples indicated in
#'   label_scheme_fit.
center_df <- function (df, label_scheme_fit, cf_x_fit, sd_coefs_fit) 
{
  ## subset `colnames(df)` according to the Sample_IDs in `label_scheme_fit`
  sub_sids <- label_scheme_fit$Sample_ID
  
  ans <- lapply(c("^N_log2_R", "^N_I", "^Z_log2_R"), hfind_fit_nms, 
                nms = names(df), sub_sids = sub_sids)
  nm_log2r_n <- ans[[1]]
  nm_int_n <- ans[[2]]
  nm_log2r_z <- ans[[3]]
  rm(list = c("ans"))
  
  ## standardize `N_log2_R` to `Z_log2_R` for samples in `sub_sids`
  centers <- cf_x_fit$x
  
  df_z <- mapply(normSD, 
                 df[, nm_log2r_n, drop = FALSE], 
                 center = centers, 
                 SD = sd_coefs_fit$fct, 
                 SIMPLIFY = FALSE) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`(gsub("^N_log2_R", "Z_log2_R", names(.))) %>%
    `rownames<-`(rownames(df))
  
  nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
  df_z[, nan_cols] <- 0
  rm(list = c("nan_cols"))
  
  ## Update `Z_log2_R` in `df`
  df[, nm_log2r_z] <- df_z
  
  if (FALSE) {
    # pre-existed
    if (length(nm_log2r_z)) {
      df[, nm_log2r_z] <- df_z
    }
    # not yet available (should not occurred)
    else {
      df <- cbind(df, df_z)
      nm_log2r_z <- names(df_z)
    }
  }

  ## !!! align `N_log2_R` and `N_I` only AFTER the calculation of "df_z"
  df[, nm_log2r_n] <- sweep(df[, nm_log2r_n, drop = FALSE], 2, centers, "-")
  df[, nm_int_n] <- sweep(df[, nm_int_n, drop = FALSE], 2, 2^centers, "/")    
  
  invisible(list(df = df, nm_log2r_z = nm_log2r_z))
}


#' Adds mean-deviation
#' 
#' @param df A data frame.
#' @param label_scheme_fit A subset of lable_scheme.
#' @param filepath A file path.
add_mean_dev <- function (df, label_scheme_fit, filepath) 
{
  filepath <- gsub("\\\\", "/", filepath)
  
  if (grepl("Peptide/", filepath)) {
    prefix <- "pep_mean_"
  } 
  else if (grepl("Protein/", filepath)) {
    prefix <- "prot_mean_"
  } 
  else {
    return(df)
  }
  
  ## matches the column names in `df`
  ans <- lapply(c("^log2_R", "^N_log2_R", "^Z_log2_R"), hfind_fit_nms, 
                nms = names(df), 
                sub_sids = label_scheme_fit$Sample_ID)
  
  nm_log2r_raw <- ans[[1]]
  nm_log2r_n <- ans[[2]]
  nm_log2r_z <- ans[[3]]
  rm(list = c("ans"))

  ## adds row_mean columns
  df <- df %>% 
    dplyr::mutate(!!paste0(prefix, "raw") := 
                    rowMeans(.[, names(.) %in% nm_log2r_raw, drop = FALSE], 
                             na.rm = TRUE), 
                  !!paste0(prefix, "n") := 
                  rowMeans(.[, names(.) %in% nm_log2r_n, drop = FALSE], 
                           na.rm = TRUE), 
                  !!paste0(prefix, "z") := 
                    rowMeans(.[, names(.) %in% nm_log2r_z, drop = FALSE], 
                             na.rm = TRUE)) 
  
  cols <- grepl(paste0("^", prefix), names(df))
  
  df[, cols] <- df[, cols, drop = FALSE] %>%
    dplyr::mutate_if(is.logical, as.numeric) %>%
    round(digits = 3L)
  
  if (prefix == "pep_mean_") {
    df <- dplyr::bind_cols(
      df %>% dplyr::select(grep("^prot_", names(.))), 
      df %>% dplyr::select(grep("^pep_", names(.))), 
      df %>% dplyr::select(-grep("^pep_|^prot_", names(.)))
    )
  } 
  else {
    df <- local({
      new_nms <- paste0(prefix, c("raw", "n", "z"))
      nms <- names(df)
      
      nm_last <- nms %>% 
        .[! . %in% new_nms] %>% 
        .[grepl(paste0("^", gsub("mean_", "", prefix)), .)] %>% 
        .[length(.)]
      
      idx <- which(nms == nm_last)
      
      dplyr::bind_cols(
        df %>% dplyr::select(1:idx), 
        df %>% dplyr::select(new_nms), 
        df %>% dplyr::select((idx+1):ncol(.)) %>% dplyr::select(-new_nms)
      )    
    })
  }
}


#' Checks \code{n_comp}
#' 
#' @inheritParams normMulGau
find_n_comp <- function (df, n_comp = NULL, method_align = "MC") 
{
  if (is.null(n_comp)) {
    if (method_align == "MGKernel") {
      n_comp <- if (nrow(df) > 3000L) 3L else 2L
    } 
    else if (method_align == "MC") {
      n_comp <- 1L
    } 
    else {
      n_comp <- 1L
    }
  } 
  else {
    if (n_comp < 1L) {
      stop("`n_comp` is not >= 1.")
    }
    
    if (method_align == "MGKernel") {
      if (n_comp < 2L) {
        stop("`n_comp` is not >= 2 at `method_align = MGKernel`.")
      }
    }
  }
  
  invisible(n_comp)
}


#' Finds the x at max density.
#' 
#' @param params Parameters from multiple Gaussians.
#' @param label_scheme The metadata of label scheme.
my_which_max <- function (params, label_scheme) 
{
  sids <- label_scheme$Sample_ID
  
  fit <- params %>%
    split(.$Sample_ID) %>%
    lapply(sumdnorm, xmin = -2, xmax = 2, by = 2/400) %>%
    do.call(rbind, .) %>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>%
    dplyr::arrange(Sample_ID)
  
  cf_x <- suppressWarnings(
    fit %>%
      dplyr::group_by(Sample_ID) %>% 
      dplyr::mutate(Max = max(Sum, na.rm = TRUE)) %>% 
      dplyr::mutate(Max = ifelse(is.infinite(Max), NA_real_, Max)) %>% 
      dplyr::filter(Sum == Max) %>%
      dplyr::mutate(x = mean(x, na.rm = TRUE)) %>% # tie-breaking
      dplyr::filter(!duplicated(x)) %>% 
      dplyr::select(-Max)
  )
  
  cf_empty <- fit %>%
    dplyr::filter(! Sample_ID %in% cf_x$Sample_ID)
  
  if (nrow(cf_empty)) {
    cf_empty <- cf_empty %>%
      dplyr::group_by(Sample_ID) %>%
      dplyr::filter(!duplicated(Sum)) %>%
      dplyr::mutate(x = 0)
  }
  
  cf_x <- dplyr::bind_rows(cf_x, cf_empty) %>% 
    data.frame(check.names = FALSE) %>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>% 
    dplyr::arrange(Sample_ID)
}


#' Compares to prior n_comp value.
#' 
#' Returns TRUE if matched in the value of n_comp; otherwise return FALSE.
#' 
#' @param filepath A file path.
#' @param filename A file name.
#' @inheritParams normMulGau
#' 
#' @return A logical value 
ok_file_ncomp <- function(filepath, filename, n_comp) 
{
  if (!file.exists(file.path(filepath, filename))) 
    return(FALSE)
  
  params <- read.table(file.path(filepath, filename), check.names = FALSE, 
                       header = TRUE, comment.char = "#")
  
  n_comp == dplyr::n_distinct(params$Component)
}


#' calculates the median-centered \eqn{x}'s positions by spline points
#'
#' Only used with \code{method_align = MC}. The values of \eqn{x}'s for sample
#' IDs indicated in \code{label_scheme_fit} will be updated for subsequent
#' adjustment in median centering. The rest of \eqn{x}'s stays the same,
#' corresponding to no further changes in center positions.
#'
#' @param df A data frame.
#' @param label_scheme The metadata of label scheme.
#' @param label_scheme_fit the subset of \code{label_scheme} used in fitting.
#' @param ... Additional parameters, including slice_dots.
#' @inheritParams normMulGau
spline_coefs <- function (df, label_scheme, label_scheme_fit, 
                          is_prot_lfq = FALSE, method_align = "MC", ...) 
{
  dots <- rlang::enexprs(...)
  
  slice_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^slice_", names(.))]
  
  # initialization: NA for Empty samples; 0 for the remaining
  sids <- label_scheme$Sample_ID
  
  pat_pre <- "^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\("
  pat_ref <- paste0(pat_pre, "(.*)\\)$")
  
  x_vals <- label_scheme %>% 
    dplyr::select(Sample_ID) %>% 
    dplyr::mutate(x = ifelse(grepl("^Empty\\.[0-9]+$", Sample_ID), NA_real_, 0), 
                  Sample_ID = factor(Sample_ID, levels = sids)) %>% 
    dplyr::arrange(Sample_ID)

  x_vals_fit <- df %>%
    filters_in_call(!!!slice_dots) %>% 
    dplyr::select(grep(pat_pre, names(.))) %>% 
    `colnames<-`(gsub(pat_ref, "\\1", names(.))) %>%
    dplyr::select(which(names(.) %in% label_scheme_fit$Sample_ID))
  
  if (is_prot_lfq && method_align == "MC") {
    x_vals_fit <- lapply(x_vals_fit, function (vals) {
      qts <- quantile(vals, seq(0, 1, .05), na.rm = TRUE)
      qt_lwr <- qts[qts >= -2][[1]]
      qt_upr <- qts[which(qts <= 2.5) %>% .[length(.)]]

      vals[vals < qt_lwr | vals > qt_upr] <- NA_real_
      median(vals, na.rm = TRUE)
    }) |>
      unlist()
  }
  else {
    x_vals_fit <- x_vals_fit %>% 
      dplyr::summarise_all(function (x) median(x, na.rm = TRUE)) %>%
      unlist()
  }
  
  x_vals_fit <- data.frame(x = x_vals_fit) |>
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) |>
    dplyr::arrange(Sample_ID)
  
  rows <- x_vals$Sample_ID %in% x_vals_fit$Sample_ID
  x_vals[rows, ] <- x_vals_fit
  
  invisible(x_vals)
}


#' Data normalization
#'
#' \code{normMulGau} normalizes \code{log2FC}.
#'
#' When executed with \link{mergePep} or link{Pep2Prn}, the \code{method_align}
#' is always \code{MC}. As a result, peptide or protein data are at first
#' median-centered.
#'
#' It is then up to \link{standPep} or \link{standPep} for alternative choices
#' in \code{method_align}, \code{col_select} etc.
#'
#' @param df An input data frame
#' @param is_prot_lfq Logical; is protein LFQ data or not. About half of the
#'   protein intensity values can be missing with LFQ and imputed with small
#'   values. The typically causes a bimodality in protein log2FC distributions
#'   and need to be handled especially at \code{method_align = "MC"}.
#' @inheritParams prnHist
#' @inheritParams standPep
#' @return A data frame.
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
normMulGau <- function(df, method_align = "MC", n_comp = NULL, seed = NULL, 
                       range_log2r = c(0, 100), range_int = c(0, 100), 
                       filepath = NULL, col_select = NULL, cut_points = Inf, 
                       is_prot_lfq = FALSE, ...) 
{
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  dat_dir <- get_gl_dat_dir()
  
  dots <- rlang::enexprs(...)
  
  slice_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^slice_", names(.))]
  
  nonslice_dots <- dots %>% 
    .[! . %in% slice_dots]
  
  # if different `n_comp` between two `method_align = MGKernel`, 
  #   force `col_select` to `Sample_ID` (all samples); 
  # if `n_comp` is given but with `method_align = MC`, 
  #   ignore difference in `n_comp`
  
  n_comp <- find_n_comp(df, n_comp, method_align)
  # ok_O_ncomp <- ok_file_ncomp(filepath, "MGKernel_params_O.txt", n_comp)
  ok_N_ncomp <- ok_file_ncomp(filepath, "MGKernel_params_N.txt", n_comp)
  ok_Z_ncomp <- ok_file_ncomp(filepath, "MGKernel_params_Z.txt", n_comp)
  
  if (method_align == "MGKernel") {
    if ((!ok_N_ncomp) && (col_select != rlang::expr(sample_ID))) 
      col_select <- rlang::expr(Sample_ID)
    
    if ((!ok_Z_ncomp) && (col_select != rlang::expr(sample_ID))) 
      col_select <- rlang::expr(Sample_ID)
  }
  
  label_scheme <- load_ls_group(dat_dir, "label_scheme")
  label_scheme_fit <- label_scheme %>% .[!is.na(.[[col_select]]), ]
  sids <- label_scheme$Sample_ID
  sub_sids <- label_scheme_fit$Sample_ID
  
  ## finds the subset of column names for fitting
  ans <- lapply(c("^N_log2_R", "^N_I"), hfind_fit_nms, names(df), sub_sids)
  nm_log2r_n <- ans[[1]]
  nm_int_n <- ans[[2]]
  rm(list = c("ans"))
  
  ## fitting
  pat_ref_n <- "^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\((.*)\\)$"
  pat_ref_z <- "^Z_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$"
  
  if (method_align == "MGKernel") {
    message("method_align = ", method_align)
    message("n_comp = ", n_comp)
    message("col_select = ", col_select)
    
    # (1.1) multi-Gaussian params
    params_sub <- df %>% 
      filters_in_call(!!!slice_dots) %>% 
      dplyr::select(nm_log2r_n) %>% 
      `names<-`(gsub(pat_ref_n, "\\1", names(.))) %>% 
      fitKernelDensity(n_comp = n_comp, seed = seed, !!!nonslice_dots) %>% 
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>% 
      dplyr::arrange(Sample_ID, Component)
    
    if (!ok_N_ncomp) {
      # previously forced `col_select = Sample_ID` if detected different `n_comp` 
      # so if not `ok`, `params_sub` must be for all samples
      params <- params_sub
    } 
    else {
      params <- read.table(file.path(filepath, "MGKernel_params_N.txt"), 
                           check.names = FALSE, header = TRUE, 
                           comment.char = "#") %>% 
        dplyr::select(names(params_sub)) %>% 
        dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>%
        dplyr::arrange(Sample_ID, Component)
      
      rows <- params$Sample_ID %in% params_sub$Sample_ID
      params[rows, ] <- params_sub
      rm(list = "rows")
    }
    
    rm(list = "params_sub")
    
    # (1.2) SDs and scaling factors
    sd_coefs  <- calc_sd_fcts(df, range_log2r, range_int, label_scheme)
    bad_coefs <- is.na(sd_coefs$fct)
    
    if (any(bad_coefs)) {
      warning("Scaling of standard deviations failed for ", 
              paste(sd_coefs$Sample_ID[bad_coefs], collapse = ", "), ".\n", 
              "Choose `scale_log2r = FALSE` with your workflow.")
      
      message("May set `range_log2r = c(0, 100)` and `range_int = c(0, 100)`.")
    }
    
    # (1.3) centers
    cf_x <- my_which_max(params, label_scheme) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>%
      dplyr::arrange(Sample_ID) 
    
    # (1.4) putting together
    list(params, cf_x, sd_coefs) %>%
      purrr::reduce(dplyr::left_join, by = "Sample_ID") %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>%
      dplyr::arrange(Sample_ID) %>%
      write.table(file = file.path(filepath, "MGKernel_params_N.txt"), 
                  sep = "\t", col.names = TRUE, row.names = FALSE)
    
    # (2.1) the subset of SDs and scaling factors
    sd_coefs_fit <- sd_coefs %>% dplyr::filter(Sample_ID %in% sub_sids)
    
    # (2.2) the subset of centers
    cf_x_fit <- cf_x %>% dplyr::filter(Sample_ID %in% sub_sids)
    
    # (2.3) data centering: (Z_log2R, N_log2_R) and intensity scaling (N_I)
    # for selected samples
    ans <- center_df(df, label_scheme_fit, cf_x_fit, sd_coefs_fit)
    df <- ans$df
    nm_log2r_z <- ans$nm_log2r_z
    rm(list = "ans")
    
    # (3) separate fits of Z_log2_R for updating curve parameters only
    if (!ok_Z_ncomp) {
      params_z <- df %>% 
        filters_in_call(!!!slice_dots) %>% 
        dplyr::select(nm_log2r_z) %>% 
        `names<-`(gsub(pat_ref_z, "\\1", names(.))) %>% 
        fitKernelDensity(n_comp = n_comp, seed = seed, !!!nonslice_dots) %>% 
        dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>% 
        dplyr::arrange(Sample_ID, Component) %>% 
        dplyr::mutate(x = 0)
    } 
    else {
      params_z_sub <- df %>% 
        filters_in_call(!!!slice_dots) %>% 
        dplyr::select(nm_log2r_z) %>% 
        `names<-`(gsub(pat_ref_z, "\\1", names(.))) %>% 
        fitKernelDensity(n_comp = n_comp, seed, !!!nonslice_dots) %>% 
        dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>% 
        dplyr::arrange(Sample_ID, Component)
      
      params_z <- read.table(file.path(filepath, "MGKernel_params_Z.txt"), 
                             check.names = FALSE, header = TRUE, 
                             comment.char = "#") %>% 
        dplyr::select(names(params_z_sub)) %>% 
        dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sids)) %>%
        dplyr::arrange(Sample_ID, Component)
      
      rows_z <- params_z$Sample_ID %in% params_z_sub$Sample_ID
      params_z[rows_z, ] <- params_z_sub
      rm(list = "rows_z")
      
      params_z$x <- 0
    }
    
    write.table(params_z, file = file.path(filepath, "MGKernel_params_Z.txt"),
                sep = "\t", col.names = TRUE, row.names = FALSE)
  } 
  else if (method_align == "MC") {
    message("method_align = ", method_align)
    message("n_comp = NULL")
    message("col_select = ", col_select)
    
    # (1) when called with `mergePep` or `Pep2Prn`, 
    #       always median-centering against all samples 
    # (2) `mergePep` called before `standPep` and `Pep2Prn` before `standPrn`
    #       -> data are at first median-centered for all samples
    
    # (MC.1) SDs and scaling factors (all samples)
    sd_coefs <- calc_sd_fcts(df, range_log2r, range_int, label_scheme)
    
    # cut_points = c(mean_lint = seq(4, 7, .5))
    # cut_points = c(prot_icover = seq(.25, .75, .25))
    # cut_points = c(prot_icover = Inf)
    # cut_points = c(prot_icover = NULL)
    # cut_points = c(prot_icover = NA)
    # cut_points = Inf
    # cut_points = NULL
    # cut_points = NA
    
    # (MC.2) df's by cut_points (all samples)
    cut_points <- set_cutpoints2(cut_points, df)
    
    if (all(is.infinite(cut_points))) {
      df <- list(df)
    } 
    else {
      df <- df %>% 
        dplyr::mutate(col_cut = !!rlang::sym(names(cut_points)[1])) %>% 
        dplyr::mutate_at(.vars = "col_cut", cut, 
                         breaks = cut_points, 
                         labels = cut_points[1:(length(cut_points)-1)]) %>%
        split(.$col_cut, drop = TRUE) %>% 
        purrr::map(~ .x %>% dplyr::select(-col_cut))
    }
    
    # (MC.3) centers 
    # (a) new x's for selected samples; 
    # (b) 0 for the remaining (no further adjustment)
    x_vals <- df %>% 
      purrr::map(spline_coefs, label_scheme, label_scheme_fit, 
                 is_prot_lfq = is_prot_lfq, method_align = method_align, 
                 !!!slice_dots)

    # (MC.4) data median-centering (selected samples)
    # (label_scheme, not label_scheme_fit as x's are 0's for non-selected)
    ans <- purrr::map2(df, x_vals, ~ center_df(.x, label_scheme, .y, sd_coefs))
    
    df <- lapply(ans, `[[`, "df") %>% 
      dplyr::bind_rows()
    
    rm(list = "ans")
  }
  
  # (1) median deviation based sample IDs in `label_scheme` not `label_scheme_fit`
  #  as sample IDs in `label_scheme_fit` may be only for mixed-bed normalization
  # (2) the `mean` also get updated after normalization against sample subsets 
  
  label_scheme_non_trivial <- label_scheme %>% 
    dplyr::filter(!Reference, !grepl("^Empty\\.[0-9]+", Sample_ID))
  
  df <- df %>% add_mean_dev(label_scheme_non_trivial, filepath)
  
  invisible(df)
}


#' Data normalization
#'
#' \code{dblTrim} doubly trims the \code{log2FC} and reporter-ion intensity by
#' the given ranges.
#' 
#' @param type_r Character string for recognizing the columns of \code{log2FC}.
#' @param type_int Character string for recognizing the columns of \code{intensity}.
#' @inheritParams info_anal
#' @inheritParams standPep
#' 
#' @return A data frame.
#' @import dplyr purrr  
#' @importFrom magrittr %>% %T>% %$% %<>% 
dblTrim <- function(df, range_log2r = c(0, 100), range_int = c(0, 100),  
                   type_r = "N_log2_R",  type_int = "N_I") 
{
  df_trim <- df
  
  type_r <- paste0("^", type_r, "[0-9]{3}")
  type_int <- paste0("^", type_int, "[0-9]{3}")
  
  # trim by log2-ratios
  col_r <- grepl(type_r, names(df_trim))
  
  df_trim[, col_r] <- lapply(df_trim[, col_r], function (x) {
    q_ratio <- quantile(x, probs = range_log2r/100, na.rm = TRUE)
    x[x < q_ratio[1] | x > q_ratio[2]] <- NA_real_
    return(x)
  }
  )
  
  # trim by intensity
  col_int <- grepl(type_int, names(df_trim))
  
  df_trim[, col_int] <- lapply(df_trim[, col_int], function (x) {
    q_intensity <- quantile(x, probs = range_int/100, na.rm = TRUE)
    x[x < q_intensity[1] | x > q_intensity[2]] <- NA_real_
    return(x)
  }
  )
  
  # doubly trim
  df_trim[!is.na(df_trim)] <- 1
  
  df_trim <- mapply(`*`, df_trim[, grepl(type_r, names(df_trim))],
                    df_trim[, grepl(type_int, names(df_trim))], 
                    SIMPLIFY = FALSE) %>%
    data.frame(check.names = FALSE)
  
  df_trim[] <- mapply(`*`, df[, grepl(type_r, names(df))] , df_trim, 
                      SIMPLIFY = FALSE)
  
  sapply(df_trim, sd, na.rm = TRUE)
}


#' Data normalization
#'
#' \code{sumdnorm} calculates summed density from \code{normMulGau}.
#' 
#' @param x A numeric vector.
#' @param xmin the miminal x values.
#' @param xmax the maximal x values.
#' @param by the step length.
#' @return A data frame.
#' @import dplyr purrr  
#' @importFrom magrittr %>% %T>% %$% %<>% 
sumdnorm <- function (x, xmin = -4, xmax = 4, by = xmax/200) 
{
  wt_dnorm <- function (x, lambda, mean, sd) lambda * dnorm(x, mean = mean, sd = sd)
  
  args <- purrr::pmap(x[, names(x) %in% c("lambda", "mean", "sd")], list) %>%
    `names<-`(x$Sample_ID)
  
  nm_comps <- paste0("G", seq_len(length(args)))
  
  Seq <- seq(xmin, xmax, by = by)
  
  lapply(args, function(args) rlang::eval_tidy(rlang::quo(wt_dnorm(Seq, !!! args)))) %>%
    do.call(cbind, .) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`(nm_comps) %>%
    dplyr::mutate(Sum = rowSums(.)) %>%
    dplyr::mutate(x = Seq) %>% #
    dplyr::mutate(Sample_ID = names(args)[1])
}


#' Data normalization
#'
#' \code{normSD} normalizes the SD of \code{log2FC}.
#' 
#' @param x A numeric vector.
#' @param center The position of \code{x} to be centered at.
#' @param SD The standard deviation that data will be scaled to.
#' @import dplyr purrr  
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @export
normSD <- function (x, center = 0, SD = 1) 
{
  if (sum(is.na(x)) == length(x))
    x
  else if ((sum(is.na(x)) + sum(x == 0, na.rm = TRUE)) == length(x)) 
    x[1:length(x)] <- NaN
  else if (all(x == 0)) 
    x
  else 
    x <- (x - center) / SD
  
  invisible(x)
}


#' Data normalization
#'
#' \code{fitKernelDensity} calculates the fitted density of \code{log2FC}.
#' 
#' @param df An input data frame.
#' @inheritParams standPep
#' @return A data frame.
#'
#' @import dplyr purrr  
#' @importFrom mixtools normalmixEM 
#' @importFrom magrittr %>% %T>% %$% %<>% 
fitKernelDensity <- function (df, n_comp = 3L, seed = NULL, ...) 
{
  dots <- rlang::enexprs(...)
  
  ok_nan_cols <- df %>% 
    purrr::map_lgl(not_all_nan, na.rm = TRUE)
  
  min_n <- df %>% 
    .[, ok_nan_cols, drop = FALSE] %>% 
    .[, not_all_NA(.), drop = FALSE] %>% 
    purrr::map_int(~ (!is.na(.x)) %>% sum()) %>% 
    min()
  
  if (min_n < 50L) 
    stop("Too few data points for fitting with multiple Gaussian functions.")
  
  lapply(df, nmix_params, n_comp, seed = seed, !!!dots) %>%
    do.call(rbind, .) %>%
    dplyr::mutate(Sample_ID = rownames(.)) %>%
    dplyr::mutate(Sample_ID = gsub("(.*)\\.\\d+$", "\\1", Sample_ID)) %>%
    dplyr::mutate(Channel = rep(1:(nrow(.)/n_comp), each = n_comp)) %>%
    dplyr::arrange(Channel, -lambda) %>%
    dplyr::mutate(Component = rep(1:n_comp, nrow(.)/n_comp)) %>%
    dplyr::mutate(Height = .$lambda * dnorm(.$mean, mean = .$mean, sd = .$sd))
}


#' Helper of \link{fitKernelDensity}.
#' 
#' @param x A numeric vector of log2Ratios under a Sample_ID.
#' @inheritParams fitKernelDensity
nmix_params <- function (x, n_comp = 3L, seed = seed, ...) 
{
  dots <- rlang::enexprs(...)

  if (!is.null(dots$k)) {
    cat(paste("k =", dots$k, "replaced by", paste("n_comp =", n_comp, "\n")))
    dots$k <- NULL
  }
  
  if (sum(is.na(x)) == length(x) |
      (sum(is.na(x)) + sum(x == 0, na.rm = TRUE)) == length(x) |
      all(x == 0)) {
    df_par <- data.frame(Component = c(1:n_comp), lambda = rep(NA, n_comp),
                         mean = rep(NA, n_comp), sd = rep(NA, n_comp))
  } 
  else {
    x <- x[!is.na(x)]
    
    stopifnot(n_comp > 1L)
    
    mixEM_call <- rlang::quo(mixtools::normalmixEM(!!x, k = !!n_comp, !!!dots))
    
    if (!is.null(seed)) set.seed(seed) else set.seed(sample(.Random.seed, 1))
    
    quietly_out <- purrr::quietly(rlang::eval_tidy)(mixEM_call, caller_env())
    x_k2 <- quietly_out$result
    
    df_par <- data.frame(Component = 1:n_comp, 
                         lambda = x_k2$lambda, 
                         mean = x_k2$mu, 
                         sd = x_k2$sigma)
  }
  
  invisible(df_par)
}


