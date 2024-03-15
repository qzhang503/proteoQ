#' Prepares data for analysis
#'
#' \code{prepDM} prepares a minimal adequate data frame for subsequent analysis.
#'
#' @inheritParams prnEucDist
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @param sub_grp Numeric.  A list of sample IDs that will be used in subsequent
#'   analysis.
#' @param type The type of data, for example ratio or intensity.
#' @return A data frame tailored for subsequent analysis.
#'
#' @examples
#' \donttest{tempData <- prepDM(df, entrez, scale_log2r, label_scheme_sub$Sample_ID)}
#' @import dplyr
#' @importFrom magrittr %>% %T>% %$% %<>%
prepDM <- function(df = NULL, id = "pep_seq", scale_log2r = TRUE, sub_grp = NULL, 
                   type = "ratio", anal_type = NULL, rm_allna = FALSE) 
{
  dat_dir <- get_gl_dat_dir()
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  
  if (!nrow(df))
    stop("Zero row of data after subsetting.")

  id <- rlang::as_string(rlang::enexpr(id))
  
  if (anal_type %in% c("ESGAGE", "GSVA"))
    id <- "entrez"

  if ((anal_type %in% c("GSEA")) && (id != "gene"))
    stop("Primary ID is not `gene`.", call. = FALSE)

  NorZ_ratios <- find_NorZ(scale_log2r)

  pattern <- 
    "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("

  df <- local({
    df <- df %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(!duplicated(!!rlang::sym(id)),
                    !is.na(!!rlang::sym(id)),) 
    
    if (!nrow(df))
      stop("Zero row of data after the removals of duplicated or NA `", id, "`.")

    df <- df %>% 
      dplyr::ungroup() %>% 
      { if (rm_allna) dplyr::filter(., rowSums(!is.na(.[grep(NorZ_ratios, names(.))])) > 0) 
        else . } %>% 
      reorderCols(endColIndex = grep(pattern, names(.)), col_to_rn = id)
    
    if (!nrow(df))
      stop("All intensity are NA or 0 in the input data.")
    
    df
  })

  Levels <- sub_grp %>%
    as.character(.) %>%
    .[!grepl("^Empty\\.[0-9]+", .)]
  
  dfR <- local({
    dfR <- df %>%
      dplyr::select(grep(NorZ_ratios, names(.))) %>%
      `colnames<-`(label_scheme$Sample_ID) %>%
      dplyr::select(which(names(.) %in% sub_grp)) 
    
    # reference will drop with single reference
    non_trivials <- dfR %>% 
      dplyr::select(which(not_all_zero(.))) %>% 
      names()
    
    trivials <- names(dfR) %>% .[! . %in% non_trivials]
    
    if (length(trivials)) {
      warning("Samples with all NA entries being dropped: ", 
              paste(trivials, collapse = ", "))
    }

    dfR <- dfR %>%
      dplyr::select(which(not_all_zero(.))) %>% 
      dplyr::select(Levels[Levels %in% names(.)]) # ensure the same order  
  })

  # dominated by log2R, no need to filter all-NA intensity rows
  dfI <- df %>%
    dplyr::select(grep("^N_I[0-9]{3}", names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sub_grp)) %>%
    dplyr::select(which(not_all_zero(.))) %>%
    dplyr::select(which(names(.) %in% names(dfR))) %>%
    dplyr::select(Levels[Levels %in% names(.)])
  
  tempI <- dfI
  rownames(tempI) <- paste(rownames(tempI), "Intensity", sep = "@")
  
  tempR <- dfR
  rownames(tempR) <- paste(rownames(tempR), "log2R", sep = "@")
  
  invisible(list(log2R = dfR, 
                 Intensity = dfI, 
                 IR = rbind(tempR, tempI)))
}


#' Finds the patterns of log2FC columns.
#'
#' @inheritParams prnHist
find_NorZ <- function (scale_log2r = TRUE)
{
  if (is.na(scale_log2r))
    "^log2_R"
  else if (scale_log2r) 
    "^Z_log2_R"
  else
    "^N_log2_R"
}


#' Replaces column names to get Sample_ID(s).
#' 
#' @param NorZ_ratios A string.
#' @param nms A vector of names.
replace_NorZ_names <- function (NorZ_ratios, nms) 
{
  gsub(paste0(NorZ_ratios, "[0-9]{3}[NC]{0,1}\\s+\\((.*)\\)$"), "\\1", nms)
}

#' Prefix form of colnames(x)[c(2, 5, ...)] for use in pipes
#'
#' \code{names_pos<-} rename the columns at the indexes of \code{pos}.
#'
#' @param x A data frame.
#' @param pos Numeric.  The index of columns for name change.
#' @param value Characters.  The new column names.
#' @return The data frame with new names.
#'
#' @import dplyr
#' @importFrom magrittr %>% %T>% %$% %<>% 
`names_pos<-` <- function(x, pos, value) 
{
  names(x)[pos] <- value
  x
}


#' Re-order file names
#'
#' \code{reorder_files} re-orders file names by TMT set numbers then by LCMS
#' injection numbers.
#'
#' @param filelist A list of file names.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
reorder_files <- function(filelist) 
{
  newlist <- NULL
  
  tmt_sets <- gsub("^TMTset(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort() 
  
  lcms_injs <- gsub("^.*_LCMSinj(\\d+).*", "\\1", filelist) %>% 
    unique() %>% 
    as.integer() %>% 
    sort()
  
  if (anyNA(tmt_sets)) {
    stop("Values under `expt_smry.xlsx::TMT_Set` need to be integers.", 
         call. = FALSE)
  }

  if (anyNA(lcms_injs)) {
    stop("Values under `expt_smry.xlsx::LCMS_Injection` need to be integers.", 
         call. = FALSE)
  }

  for (i in tmt_sets) {
    newlist <- c(newlist, filelist[grep(paste0("TMTset", i, "_"), filelist, 
                                        ignore.case = TRUE)])
  }

  return(newlist)
}


#' Re-order columns in a data frame
#'
#' \code{reorderCols} re-orders columns in a data frame.
#'
#' @param df A data frame.
#' @param endColIndex the indexes of columns to be moved to the end of
#'   \code{df}.
#' @param col_to_rn the column identifier where the values under that column
#'   will be used as row names.
#' @import dplyr 
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
reorderCols <- function(df, endColIndex, col_to_rn) 
{
  if (!length(endColIndex)) {
    endColIndex <- grep("I[0-9]{3}|log2_R[0-9]{3}", names(df))
  }
  
  df <- cbind(df[, -endColIndex], df[, endColIndex])
  
  if (sum(duplicated(df[[col_to_rn]])) > 0) {
    df <- df %>%
      dplyr::group_by(!!ensym(col_to_rn)) %>%
      dplyr::mutate(Index = row_number()) %>%
      dplyr::mutate(indexed_names = paste(get(col_to_rn), Index, sep=".")) %>%
      data.frame(check.names = FALSE) %>%
      `rownames<-`(.[, "indexed_names"]) %>%
      dplyr::select(-c("Index", "indexed_names"))
  } else {
    rownames(df) <- df[[col_to_rn]]
  }
  
  return(df)
}


#' Re-order columns in a data frame
#'
#' \code{reorderCols2} re-orders columns in a data frame.
#'
#' @param df A data frame.
#' @param pattern columns matched the pattern will be moved to the end of
#'   \code{df}.
#' @import dplyr 
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @export
reorderCols2 <- function(df = NULL, pattern = NULL) 
{
  if (is.null(df)) 
    stop("`df` cannot be `NULL`.", call. = FALSE)

  if (is.null(pattern)) {
    pattern <- 
      "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("
  }

  endColIndex <- grep(pattern, names(df))
  
  if (length(endColIndex)) 
    df <- dplyr::bind_cols(df[, -endColIndex], df[, endColIndex])
  
  return(df)
}


#' Re-order columns in a data frame
#'
#' \code{ins_cols_after} re-orders columns in a data frame.
#'
#' @param df A data frame.
#' @param idx_bf A column index for the insertion of columns after.
#' @param idx_ins Column index(es) for the columns to be inserted after
#'   \code{idx_bf}.
#' @import dplyr 
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
ins_cols_after <- function(df = NULL, idx_bf = ncol(df), idx_ins = NULL) 
{
  if (is.null(df)) 
    stop("`df` cannot be `NULL`.", call. = FALSE)
  
  if (is.null(idx_ins)) 
    return(df)
  
  if (idx_bf >= ncol(df)) 
    return(df)
  
  col_nms <- names(df)[idx_ins]
  
  dplyr::bind_cols(
    df[, 1:idx_bf, drop = FALSE], 
    df[, idx_ins, drop = FALSE], 
    df[, (idx_bf + 1):ncol(df), drop = FALSE] %>% dplyr::select(-col_nms),
  )
}


#' Sorts by the indexes of TMT_Set and LCMS_Injection
#' 
#' "8.3"    "8.5"   "11.1"   "11.3"   "99.2"  "99.10" 
#' Not "8.3"    "8.5"   "11.1"   "11.3"   "99.10"  "99.2" 
#' Nor "11.1"   "11.3"  "8.3"    "8.5"   "99.10"  "99.2" 
#' @param nms A list of names of "TMT_Set" and "LCMS_Injection" separated by dot.
sort_tmt_lcms <- function (nms) 
{
  nms %>% 
    purrr::map(~ strsplit(.x, "[.]") %>% unlist()) %>% 
    do.call(rbind, .) %>% 
    data.frame() %>% 
    `colnames<-`(c("TMT_Set", "LCMS_Injection")) %>% 
    dplyr::mutate(TMT_Set = as.numeric(TMT_Set), 
                  LCMS_Injection = as.numeric(LCMS_Injection)) %>% 
    dplyr::arrange(TMT_Set, LCMS_Injection) %>% 
    dplyr::mutate(tmt_inj = paste(TMT_Set, LCMS_Injection, sep = ".")) %>% 
    dplyr::select(tmt_inj) %>% 
    unlist()
}


#' Replaces zero intensity with NA
#'
#' \code{na_zeroIntensity} replaces zero intensity with NA to avoid -Inf in
#' log10 transformation.
#'
#' @param df A list of file names.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
na_zeroIntensity <- function(df) 
{
	ind <- grep("I[0-9]{3}|^LFQ\\sintensity\\s|^Intensity\\s", names(df))

	temp <- df[, ind]
	temp[temp == 0] <- NA
	df[, ind] <- temp

	invisible(df)
}

#' Finds all zero rows
#' 
#' @param x A data frame.
not_allzero_rows <- function(x) (rowSums(x != 0, na.rm = TRUE) > 0)

#' Summarizes numeric values
#'
#' \code{aggrNums} summarizes \code{log2FC} and \code{intensity} by the
#' descriptive statistics of \code{c("mean", "median", "weighted.mean",
#' "top.3")}
#'
#' @param f A function for data summarization.
#' @examples \donttest{df_num <- aggrNums(median)(df, prot_acc, na.rm = TRUE)}
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
aggrNums <- function(f) 
{
  function (df, id, ...) {
    id <- rlang::as_string(rlang::enexpr(id))
    dots <- rlang::enexprs(...)
    
    df %>%
      dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
      dplyr::group_by(!!rlang::sym(id)) %>%
      dplyr::summarise_all(~ f(.x, !!!dots))
  }
}


#' Sums top_n
#'
#' \code{aggrTopn} summarizes \code{log2FC} and \code{intensity} by the
#' descriptive statistics of \code{c("mean", "median", "sum")}. Note the
#' difference to \link{TMT_top_n}, which uses mean statistics.
#'
#' @param f A function for data summarization.
#' @examples \donttest{df_num <- aggrTopn(sum)(df, prot_acc, 3, na.rm = TRUE)}
aggrTopn <- function(f) 
{
  function (df, id, n, ...) {
    id <- rlang::as_string(rlang::enexpr(id))
    dots <- rlang::enexprs(...)
    
    df %>%
      dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
      dplyr::mutate(sum_Intensity = rowSums(.[grep("^I[0-9]{3}", names(.))], 
                                            na.rm = TRUE)) %>%
      dplyr::group_by(!!rlang::sym(id)) %>%
      dplyr::top_n(n = n, wt = sum_Intensity) %>%
      dplyr::select(-sum_Intensity) %>%
      dplyr::summarise_all(~ f(.x, !!!dots))
  }
}


#' Aggregation of LFQ intensity.
#'
#' Summarizes \code{log2FC} and \code{intensity} by the descriptive statistics
#' of \code{c("mean", "median", "sum")}. Note that data are already subset by
#' top_n.
#' 
#' \code{ids} are a list of column keys.
#'
#' @param f A function for data summary.
#' @examples \donttest{df_num <- aggrLFQ(sum)(df, prot_acc, na.rm = TRUE)}
aggrLFQs <- function(f) 
{
  function (df, ids, ...) {
    nms <- names(df)
    
    bads <- ids[!ids %in% nms]
    
    if (length(bads))
      stop("Columns not found `", paste(bads, collapse = ", "), "` not found.")
    
    cols <- grep("log2_R[0-9]{3}|I[0-9]{3}", nms)
    
    if (!length(cols))
      stop("Columns of log2 ratios or intensity not found.")
    
    df %>%
      dplyr::select(ids, cols) |>
      dplyr::group_by_at(ids) |>
      dplyr::summarise_all(function (x) f(x, ...))
  }
}


#' Calculates weighted mean
#'
#' \code{tmt_wtmean} calculates the weighted mean of \code{log2FC} based on
#' \code{intensity}.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param id The variable to summarize \code{log2FC}.
#' @param na.rm Logical; if TRUE, removes NA values.
#' @param ... Additional arguments for \code{weighted.mean}.
#' @import dplyr 
#' @importFrom stringr str_length
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %T>% %$% %<>% 
tmt_wtmean <- function (x, id, na.rm = TRUE, ...) 
{
  my_trim <- function (x, range_int = c(.05, .95), na.rm = TRUE) {
    if (range_int[2] > 1) range_int <- range_int/100
    min_x <- min(x, na.rm = na.rm)
    qs <- quantile(x, probs = range_int, na.rm = na.rm)
    x[x < qs[1] | x > qs[2]] <- min_x
    return(x)
  }
  
  my_wtmean <- function (x) {
    x %>% 
      dplyr::mutate_at(.vars = "Intensity", my_trim, range_int = range_int, 
                       na.rm = TRUE) %>% 
      dplyr::mutate(n = stringr::str_length(.[[id]])) %>%
      dplyr::filter(n != 0) %>%
      dplyr::group_by(!!rlang::sym(id), variable) %>%
      dplyr::summarise(log2_R = weighted.mean(log2_R, log10(Intensity), 
                                              !!!dots)) %>%
      tidyr::spread(variable, log2_R)
  }
  
  
  range_int <- c(.05, .95)
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

  # intensity
  xi <- x %>%
    dplyr::select(id, grep("^I[0-9]{3}[NC]{0,1}", names(.))) %>%
    tidyr::gather(key = variable, value = Intensity, -id)
  
  # mean: `I` and `N_I`
  xi_wt <- x %>%
    dplyr::select(id, grep("^[N]{0,1}[_]{0,1}I[0-9]{3}[NC]{0,1}", names(.))) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarise_all(mean, na.rm = TRUE, !!!dots) 

  # weighted.mean `log2_R`
  xr <- x %>%
    dplyr::select(id, grep("^log2_R[0-9]{3}[NC]{0,1}", names(.))) %>%
    tidyr::gather(key = variable, value = log2_R, -id)
  
  if ("log2_R" %in% names(xr)) {
    xr_wt <- cbind.data.frame(xi, log2_R = xr[, c("log2_R")]) %>% 
      dplyr::mutate(variable = gsub("^I([0-9]{3}[NC]{0,1})", "log2_R\\1", variable), 
                    variable = factor(variable, levels = unique(variable))) %>% 
      my_wtmean()
  } else {
    xr_wt <- xi_wt[, 1, drop = FALSE]
  }

  # weighted.mean `N_log2_R` 
  xnr <- x %>%
    dplyr::select(id, grep("^N_log2_R[0-9]{3}[NC]{0,1}", names(.))) %>%
    tidyr::gather(key = variable, value = log2_R, -id)
  
  if ("log2_R" %in% names(xnr)) {
    xnr_wt <- cbind.data.frame(xi, log2_R = xnr[, c("log2_R")]) %>% 
      dplyr::mutate(variable = gsub("^I([0-9]{3}[NC]{0,1})", "N_log2_R\\1", variable), 
                    variable = factor(variable, levels = unique(variable))) %>% 
      my_wtmean()
  } else {
    xnr_wt <- xi_wt[, 1, drop = FALSE]
  }

  # weighted.mean `Z_log2_R` may not be available for PSM data
  xzr <- x %>%
    dplyr::select(id, grep("^Z_log2_R[0-9]{3}[NC]{0,1}", names(.))) %>%
    tidyr::gather(key = variable, value = log2_R, -id)
  
  if ("log2_R" %in% names(xzr)) {
    xzr_wt <- cbind.data.frame(xi, log2_R = xzr[, c("log2_R")]) %>% 
      dplyr::mutate(variable = gsub("^I([0-9]{3}[NC]{0,1})", "Z_log2_R\\1", variable), 
                    variable = factor(variable, levels = unique(variable))) %>% 
      my_wtmean()    
  } else {
    xzr_wt <- xi_wt[, 1, drop = FALSE]
  }

  x <- list(xi_wt, xr_wt, xnr_wt, xzr_wt) %>% 
    purrr::reduce(dplyr::left_join, by = id)

  dplyr::bind_cols(
    x[, names(x) == id],
    x[, grepl("^I[0-9]{3}[NC]{0,1}", names(x))],
    x[, grepl("^N_I[0-9]{3}[NC]{0,1}", names(x))],
    x[, grepl("^log2_R[0-9]{3}[NC]{0,1}", names(x))],
    x[, grepl("^N_log2_R[0-9]{3}[NC]{0,1}", names(x))],
    x[, grepl("^Z_log2_R[0-9]{3}[NC]{0,1}", names(x))],
  )
}


#' Average top_n
#'
#' \code{TMT_wt_mean} calculates the weighted mean of \code{log2FC} and
#' \code{intensity}. Note the difference to \link{aggrTopn}, which uses sum
#' statistics.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param id The variable to summarize \code{log2FC}.
#' @param ... Additional arguments for \code{mean}.
#' @examples \donttest{df_num <- TMT_top_n(df, prot_acc, na.rm = TRUE)}
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
TMT_top_n <- function (x, id, ...) 
{
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  x %>%
    dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
    dplyr::mutate(sum_Intensity = rowSums(.[grep("^I[0-9]{3}", names(.))], 
                                          na.rm = TRUE)) %>%
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::top_n(n = 3, wt = sum_Intensity) %>%
    dplyr::select(-sum_Intensity) %>%
    dplyr::summarise_all(~ mean(.x, !!!dots))
}


#' Finds not-all-zero column(s)
#'
#' \code{not_all_zero} identifies the column indexes with all NA values.
#'
#' @param df A data frame.
#' @examples \donttest{not_all_zero(mtcars)}
#' @return A logical vector
not_all_zero <- function (df) 
{
  stopifnot(is.data.frame(df))
  colSums(df != 0, na.rm = TRUE) > 0
}


#' Finds not-all-NA column(s)
#'
#' \code{not_all_NA} identifies the column indexes with all NA values.
#'
#' @param df A data frame.
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @examples \donttest{not_all_NA(mtcars)}
#' @return A logical vector
not_all_NA <- function (df) 
{
  stopifnot(is.data.frame(df))
  colSums(!is.na(df), na.rm = TRUE) > 0
}


#' Finds non-all-NaN column(s)
#'
#' @param x A vector.
#' @param ... Additional arguments for \code{sum}.
#' @examples \donttest{purrr::map_lgl(mtcars, not_all_nan, na.rm = TRUE)}
not_all_nan <- function(x, ...) 
{
  stopifnot(is.vector(x))
  sum(is.nan(x), ...) != length(x)
}


#' Finds all-NaN column(s)
#' 
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param ... The same in \code{sum}.
#' @examples \donttest{purrr::map_lgl(mtcars, is_all_nan, na.rm = TRUE)}
is_all_nan <- function(x, ...) 
{
  stopifnot(is.vector(x))
  sum(is.nan(x), ...) == length(x)
}


#' Finds a trivial column
#' 
#' @param df A data frame of \code{log2FC} and \code{intensity}.
#' @param col The index of column.
#' @examples \donttest{is_trivial_col(mtcars, 1)}
is_trivial_col <- function (df, col) 
{
  stopifnot(length(col) == 1L)
  
  x <- df[, col]
  
  x[x == 0] <- NA
  x[is.nan(x)] <- NA
  
  all(is.na(x))
}


#' Finds a trivial vector
#' 
#' @param x A numeric vector
is_trivial_dbl <- function (x) 
{
  stopifnot(is.vector(x))
  
  if (!is.numeric(x)) 
    return(FALSE)
  
  x[x == 0] <- NA
  x[is.nan(x)] <- NA
  
  all(is.na(x))
}


#' Replaces a trivial vector with NA values.
#' 
#' @param x A numeric vector
replace_trivial_with_na <- function (x) 
{
  if (is_trivial_dbl(x)) {
    x <- rep(NA, length(x))
  }

  invisible(x)
}


#' Sets up the column annotation in heat maps
#' 
#' @inheritParams prnEucDist
#' @inheritParams gspa_colAnnot
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
colAnnot <- function (annot_cols = NULL, sample_ids = NULL, annot_colnames = NULL) 
{
  if (is.null(annot_cols)) 
    return(NA)
  
  dat_dir <- get_gl_dat_dir()
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  exists <- annot_cols %in% names(label_scheme)
  
  if (sum(!exists) > 0L) {
    warning(paste0("Column '", annot_cols[!exists], "'",
                   " not found in \"label_scheme\" and will be skipped."), 
            call. = FALSE)
    annot_cols <- annot_cols[exists]
  }
  
  if (!length(annot_cols)) 
    return(NA)
  
  x <- label_scheme %>%
    dplyr::filter(Sample_ID %in% sample_ids) %>%
    dplyr::select(annot_cols, Sample_ID) 
  
  idx <- !not_all_NA(x)
  
  if (sum(idx) > 0L) 
    stop("No aesthetics defined under column(s) `", 
         paste(names(x)[idx], collapse = ", "), "`.")
  
  x <- x %>%
    dplyr::filter(!grepl("^Empty\\.", .[["Sample_ID"]]),
                  !is.na(.[["Sample_ID"]])) %>%
    data.frame(check.names = FALSE) %>%
    `rownames<-`(.[["Sample_ID"]])
  
  if (any(duplicated(x[["Sample_ID"]]))) 
    stop("Duplicated sample IDs found.\n")
  
  if (!"Sample_ID" %in% annot_cols) 
    x <- x %>% dplyr::select(-Sample_ID)
  
  if (!ncol(x)) 
    return(NA)
  
  if ("TMT_Set" %in% names(x)) {
    x <- x %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(TMT_Set = as.factor(TMT_Set)) %>%
      `rownames<-`(NULL) %>% 
      tibble::column_to_rownames(var = "rowname")
  }
  
  x <- x %>% 
    mutate_if(is.logical, factor)
  
  invisible(x)
}


#' Customizes the colors in column annotation
#' 
#' @param annotation_col The same as in \link[pheatmap]{pheatmap}.
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
setHMColor <- function (annotation_col) 
{
  ncol <- ncol(annotation_col)
  
  if (is.null(ncol)) 
    return (NA)
  
  if (!ncol) 
    return (NA)
  
  suppressWarnings(
    annotation_col <- annotation_col %>%
      dplyr::mutate_at(vars(one_of("TMT_Set")), ~ factor(TMT_Set)) %>%
      dplyr::mutate_if(is.character, as.factor)  
  )
  
  palette <- lapply(names(annotation_col), function(x) {
    n <- nlevels(annotation_col[[x]])
    
    palette <- if(n < 9L && n >= 3L) 
      brewer.pal(n, name = "Set2") # Set2: maximum 8 
    else if(n >= 9L) 
      colorRampPalette(brewer.pal(n = 7L, "Set1"))(n)
    else if(n == 2L) 
      c("#66C2A5", "#FC8D62")
    else if(n == 1L) 
      c("#66C2A5")
    else if(n == 0L) 
      colorRampPalette(brewer.pal(n = 9L, "YlOrBr"))(100)
    
    names(palette) <- levels(annotation_col[[x]])
    
    palette
  })
  
  names(palette) <- names(annotation_col)
  
  palette
}


#' Sets the upper and the lower limits in the range of heat maps
#'
#' \code{setHMlims} imputes values beyond the limits to the corresponding
#' limits.
#'
#' @param x A data frame of \code{log2FC}.
#' @param xmin the lower limit.
#' @param xmax the upper limit.
#' @import dplyr 
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
setHMlims <- function (x, xmin, xmax) 
{
  x[x < xmin] <- xmin
  x[x > xmax] <- xmax
  
  invisible(x)
}


#' Calculate the log2-ratio to the control group of samples
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @param nm_ctrl The names of samples that belong to the control group.
#' @examples \donttest{ratio_toCtrl(df, "gene", label_scheme_sub, Heatmap_Group)}
ratio_toCtrl <- function(df, id, label_scheme_sub, nm_ctrl) 
{
  id <- rlang::as_string(rlang::enexpr(id))
  
  nm_ctrl <- rlang::as_string(rlang::ensym(nm_ctrl))
  
  x <- df %>%
    dplyr::select(which(names(.) %in% label_scheme_sub$Sample_ID)) %>%
    `rownames<-`(df[[id]])
  
  col_ratio <- which(names(x) %in% label_scheme_sub$Sample_ID)
  
  col_ctrl <- which(names(x) %in%
                      label_scheme_sub[!is.na(label_scheme_sub[[nm_ctrl]]), "Sample_ID"])
  
  x[, col_ratio] <- sweep(x[, col_ratio], 1,
                          rowMeans(x[, col_ctrl, drop = FALSE], na.rm = TRUE), "-")
  
  x <- x %>% tibble::rownames_to_column(id)
  
  df %>%
    dplyr::select(which(!names(df) %in% label_scheme_sub$Sample_ID)) %>%
    dplyr::left_join(x, by = id)	%>%
    `rownames<-`(.[[id]])
}


#'Imputation of NA values
#'
#'
#'\code{imputeNA} imputes NA \code{log2FC} values in protein or peptide data.
#'Users should avoid call the method directly, but instead use the following
#'wrappers.
#'
#'@param overwrite Logical. If true, overwrite the previous results. The default
#'  is FALSE.
#'@inheritParams prnHist
#'@inheritParams info_anal
#'@param ... Parameters for \code{\link[mice]{mice}}
#'@return \code{Peptide_impNA.txt} for peptide data and \code{Protein_impNA.txt}
#'  for protein data.
#'
#'@import dplyr purrr  
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
imputeNA <- function (id, overwrite = FALSE, ...) 
{
  my_mice <- function (data, ...) {
    data <- rlang::enexpr(data)
    dots <- rlang::enexprs(...)
    
    rlang::eval_bare(rlang::expr(mice::mice(data = !!data, !!!dots)), 
                     env = rlang::caller_env())
  }
  
  handleNA <- function (x, ...) {
    ind <- purrr::map(x, is.numeric) %>% unlist()
    
    if (sum(ind) < ncol(x)) cat(names(ind)[!ind], "skipped.")
    
    # handle special column names
    nm_orig <- names(x[, ind])
    
    x[, ind] <- x[, ind] %>%
      data.frame(check.names = TRUE) %>%
      my_mice(...) %>%
      mice::complete(1) %>%
      `colnames<-`(nm_orig)
    
    x
  }
  
  
  dat_dir <- get_gl_dat_dir()
  
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("\n============================================================", 
         "\nNeed package \"mice\" for this function to work.",
         "\n============================================================",
         call. = FALSE)
  }
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    src_path <- file.path(dat_dir, "Peptide", "Peptide.txt")
    filename <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
  } else if (id %in% c("prot_acc", "gene")) {
    src_path <- file.path(dat_dir, "Protein", "Protein.txt")
    filename <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
  }
  
  if ((!file.exists(filename)) || overwrite) {
    df <- tryCatch(read.csv(src_path, check.names = FALSE, 
                            header = TRUE, sep = "\t",
                            comment.char = "#"), 
                   error = function(e) NA)
    
    if (!is.null(dim(df))) {
      message(paste("File loaded:", src_path))
    } else {
      stop("File or directory not found: ", src_path, 
           call. = FALSE)
    }
    
    df[, grep("N_log2_R", names(df))] <- 
      handleNA(df[, grep("N_log2_R", names(df))], ...)
    
    fn_params <- file.path(dat_dir, "Protein/Histogram", "MGKernel_params_N.txt")
    
    if (file.exists(fn_params)) {
      cf_SD <- read.csv(fn_params, check.names = FALSE, 
                        header = TRUE, sep = "\t", comment.char = "#") %>%
        dplyr::filter(!duplicated(Sample_ID)) %>%
        dplyr::select(Sample_ID, fct)
      
      df[, grep("Z_log2_R", names(df))] <-
        mapply(normSD,
               df[, grepl("^N_log2_R[0-9]{3}", names(df))], 
               center = 0, SD = cf_SD$fct, SIMPLIFY = FALSE) %>%
        data.frame(check.names = FALSE) %>%
        `colnames<-`(gsub("N_log2", "Z_log2", names(.)))
    } else {
      df[, grep("Z_log2_R", names(df))]  <- df[, grep("Z_log2_R", names(df))] %>% 
        handleNA(...)
    }
    
    if (any(duplicated(df[[id]]))) {
      if (id == "pep_seq") {
        warning("\`pep_seq\` is not uqique for rownames; ", 
                "uses \`pep_seq_mod\` instead.\n", 
                call. = FALSE)
        rownames(df) <- df[["pep_seq_mod"]]
      }
      if (id == "gene") {
        warning("\`gene\` is not uqique for rownames; ", 
                "uses \`prot_acc\` instead.\n", 
                call. = FALSE)
        rownames(df) <- df[["prot_acc"]]
      }
    } else {
      rownames(df) <- df[[id]]
    }
    
    write.table(df, filename, sep = "\t", col.names = TRUE, row.names = FALSE)
  } else {
    warning("NA imputation has been previously performed!\n", call. = FALSE)
  }
}


#'Imputes NA values for peptide data
#'
#'\code{pepImp} is a wrapper of \code{\link{imputeNA}} for peptide data with
#'auto-determination of \code{id}.
#'
#'@rdname imputeNA
#'
#' @examples 
#' \donttest{
#' # ===================================
#' # Peptide NA imputation
#' # ===================================
#' pepImp(
#'   m = 3,
#'   maxit = 3,
#' )
#' }
#' 
#'@export
pepImp <- function (...) 
{
  check_dots(c("id"), ...)
  id <- match_call_arg(normPSM, group_psm_by)
  imputeNA(id = !!id, ...)
}


#'Impute NA values for protein data
#'
#'\code{prnImp} is a wrapper of \code{\link{imputeNA}} for protein data with
#'auto-determination of \code{id}.
#'
#'@rdname imputeNA
#'
#' @examples
#' \donttest{
#' # ===================================
#' # Protein NA imputation
#' # ===================================
#' prnImp(
#'   m = 5,
#'   maxit = 5,
#' )
#' }
#'
#'@export
prnImp <- function (...) 
{
  check_dots(c("id"), ...)
  id <- match_call_arg(normPSM, group_pep_by)
  imputeNA(id = !!id, ...)
}


#' Species lookup
#' 
#' @inheritParams load_dbs
sp_lookup <- function(species) 
{
  switch(species, 
         human = "hs",
         mouse = "mm",
         rat = "rn",
         unknown = "unknown",
         species)
}


#' Taxonomy lookup
#' @inheritParams load_dbs
taxid_lookup <- function(species) 
{
  switch (species,
          human = 9606, 
          mouse = 10090,
          rat = 10116, 
          unknown = 999999, 
          999999)
}


#' Reversed taxonomy lookup
#' @inheritParams load_dbs
taxid_lookup_rev <- function(species) 
{
  switch (species,
          "9606" = "human", 
          "10090" = "mouse",
          "10116" = "rat", 
          "unknown")
}


#' Species lookup UpperLower (Ul)
#' 
#' @inheritParams load_dbs
sp_lookup_Ul <- function(species) 
{
  switch(species, 
         human = "Hs",
         mouse = "Mm",
         rat = "Rn",
         unknown = "Unknown",
         "unknown")
}


#' Adds gene IDs for RefSeq
#' 
#' @inheritParams add_entrez
add_refseq_gene <- function (acc_lookup, acc_type) 
{
  sp_map <- c(
    human = "hs",
    mouse = "mm",
    rat = "rn"
  )
  
  acc_map <- c(
    uniprot_acc = "uniprot_",
    uniprot_id = "uniprot_", 
    refseq_acc = "refseq_"
  )
  
  species <- unique(acc_lookup$species) %>% .[!is.na(.)]
  abbr_sp <- sp_map[species]
  abbr_acc <- acc_map[acc_type]
  
  stopifnot(acc_type == "refseq_acc")
  
  entrez_key <- switch(acc_type, 
                       uniprot_acc = "uniprot_acc", 
                       uniprot_id = "uniprot_acc", 
                       refseq_acc = "refseq_acc", 
                       other_acc = "other_acc", 
                       stop("Unknown `accession type`.", 
                            Call. = FALSE))

  # (`species` can be empty)
  if (length(species)  && all(species %in% c("human", "mouse", "rat"))) {
    filelist <- paste0(abbr_acc, "entrez_", abbr_sp)
    data(package = "proteoQ", list = filelist, envir = environment())
    
    entrez_db <- purrr::map(filelist, ~ get(.x)) %>% 
      dplyr::bind_rows() %>% 
      dplyr::filter(!duplicated(.[[entrez_key]])) %>% 
      dplyr::mutate(entrez = as.numeric(entrez),
                    !!entrez_key := as.character(!!rlang::sym(entrez_key))) %>% 
      dplyr::select(entrez_key, "gene")

    acc_lookup <- dplyr::left_join(acc_lookup, entrez_db, 
                                   by = c("prot_acc" = entrez_key)) 
  } else {
    acc_lookup <- acc_lookup %>% 
      dplyr::mutate(gene = NA_character_) 
  }
  
  invisible(acc_lookup)
}


#' Adds Entrez IDs
#' 
#' @param acc_lookup A data frame of protein accession lookups
#' @param acc_type The type of protein accessions
#' @inheritParams parse_fasta
add_entrez <- function (acc_lookup, acc_type, warns = TRUE) 
{
  sp_map <- c(
    human = "hs",
    mouse = "mm",
    rat = "rn"
  )
  
  acc_map <- c(
    uniprot_acc = "uniprot_",
    uniprot_id = "uniprot_", 
    refseq_acc = "refseq_"
  )
  
  def_sps <- c("human", "mouse", "rat")
  species <- unique(acc_lookup$species) %>% .[!is.na(.)]
  all_ok_species <- all(species %in% def_sps)
  
  if (length(species) && (!all_ok_species)) {
    cat("\n")
    
    if (warns) {
      warning("No default `entrez` lookups for species other than: ", 
              paste(def_sps, collapse = ", "), ".\n",
              "To annotate, provide the file path and name(s) via argument `entrez`.\n", 
              "See also `?Uni2Entrez` and `?Ref2Entrez` for custom `entrez` lookups.", 
              call. = FALSE)
      
      if (!all_ok_species) {
        cat("\n")
        
        warning("At `entrez = NULL`, no entrez annoation of ", 
                "non-default species: \n", 
                paste(species[!species %in% def_sps], collapse = ", "), 
                call. = FALSE)
      }
    }

    species <- species %>% .[. %in% def_sps]
  }
  
  if (!length(species)) {
    warning("No default species found for entrez annotation at ", 
            acc_type, ".",  
            call. = FALSE)
    
    acc_lookup <- acc_lookup %>% dplyr::mutate(entrez = NA_integer_)
    
    return(acc_lookup)
  } 
  
  abbr_sp <- sp_map[species] %>% .[!is.na(.)]
  abbr_acc <- acc_map[acc_type]
  
  entrez_key <- switch(acc_type, 
                       uniprot_acc = "uniprot_acc", 
                       uniprot_id = "uniprot_acc", 
                       refseq_acc = "refseq_acc", 
                       other_acc = "other_acc", 
                       stop("Unknown `accession type`.", Call. = FALSE))

  # multiple uniprot_acc(s) can share the same entrez id
  filelist <- paste0(abbr_acc, "entrez_", abbr_sp)
  
  # four columns: 
  # 1. uniprot_acc/refseq_acc/other_acc  2. gene  3. entrez  4. species
  data(package = "proteoQ", list = filelist, envir = environment())
  
  db <- purrr::map(filelist, ~ get(.x)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(!duplicated(.[[entrez_key]])) %>% 
    dplyr::mutate(entrez = as.integer(entrez),
                  !!entrez_key := as.character(!!rlang::sym(entrez_key))) %>% 
    dplyr::select(entrez_key, "entrez")
  
  acc_lookup <- dplyr::left_join(acc_lookup, db, by = entrez_key)
}


#' Adds custom Entrez IDs and optional overwriting species
#' 
#' @inheritParams  add_entrez
#' @inheritParams normPSM
#' @inheritParams annotKin
add_custom_entrez <- function(acc_lookup, entrez, acc_type) 
{
  if (!all(file.exists(entrez))) 
    stop("Wrong `entrez` file path(s) or name(s).", 
         "\nMake sure that the file type is `.rds`, not `.rda`", 
         call. = FALSE)
  
  entrez_key <- switch(acc_type, 
                       uniprot_acc = "uniprot_acc", 
                       uniprot_id = "uniprot_acc", 
                       refseq_acc = "refseq_acc", 
                       other_acc = "other_acc", 
                       stop("Unknown `accession type`.", 
                            Call. = FALSE))

  cols_ess <- c(entrez_key, "entrez")
  cols_out <- c(entrez_key, "entrez", "species")
  
  db <- purrr::map(entrez, ~ {
    db <- readRDS(.x)
    
    if (!all(cols_ess %in% names(db))) 
      stop("Need columns '", 
           purrr::reduce(cols_ess, paste, sep = ", "), 
           "' in ", .x, ".", 
           call. = FALSE)
    
    if (!length(db)) 
      stop("Empty file: ", .x, call. = FALSE)
    
    db %>% 
      dplyr::select(which(names(.) %in% cols_out)) %>% 
      { if ("species" %in% names(.)) . else .[["species"]] <- NA_character_ } %>% 
      dplyr::mutate(entrez = as.integer(entrez)) %>% 
      dplyr::select(cols_out)
  }) %>% 
    dplyr::bind_rows() 
  
  # higher priority for custom "species" in `db`
  if (all(is.na(db$species))) 
    db[["species"]] <- NULL
  else 
    acc_lookup[["species"]] <- NULL
  
  acc_lookup %>% dplyr::left_join(db, by = entrez_key)
} 


#' Determines the protein accession type from PSM tables
#'
#' Find the protein accession from a non-cRAP entry and parse it.
#' 
#'@param df An input data frame.
parse_acc <- function(df) 
{
  stopifnot("prot_acc" %in% names(df))
  
  prot_accs <- unique(df$prot_acc)
  
  pat_uni_acc <- 
    "^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  pat_uni_id <- 
    "^[[:alnum:]]+_[A-Z]{1,10}$"
  pat_ref_acc <- 
    "^[XNY]{1}[MRP]{1}_"
  pat_uni_both <- "^..\\|[[:alnum:]]+\\|[^_]+_[A-Z]{1,10}$"
  
  acc_types <- purrr::map_chr(prot_accs, ~ {
    if (grepl(pat_uni_id, .x)) 
      "uniprot_id"
    else if (grepl(pat_ref_acc, .x)) 
      "refseq_acc"
    else if (grepl(pat_uni_acc, .x)) 
      "uniprot_acc"
    else 
      "other_acc"
  }, prot_accs)
  
  df %>% 
    dplyr::left_join(data.frame(prot_acc = prot_accs, acc_type = acc_types), 
                     by = "prot_acc")
}


#' Parses FASTA for accession lookups
#'
#' NA genes are replaced with prot_acc
#' 
#' @param df An input data frame.
#' @param warns Logical; if TRUE, show warning message(s).
#' @inheritParams add_entrez
#' @inheritParams normPSM
#' @seealso \code{\link{read_fasta}} for the definition of fasta_name(s).
#' @return A lookup table, \code{acc_lookup}.
parse_fasta <- function (df, fasta, entrez, warns = TRUE) 
{
  my_lookup <- c(
    "Homo sapiens" = "human",
    "Mus musculus" = "mouse",
    "Rattus norvegicus" = "rat")

  pat_uni_acc <- "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  pat_uni_id <- "^[[:alnum:]]+_[A-Z]{1,10}$"
  pat_ref_acc <- "^[XNY]{1}[MRP]{1}_"
  
  dat_dir <- get_gl_dat_dir()
  
  # =======================================================================================
  # Two assumptions:
  # (a) `acc_type` parsed from `prot_acc` with rules of UniProt accession etc.,
  #   and must be in one of c("uniprot_acc", "uniprot_id", "refseq_acc", "other_acc") 
  #   with NA `acc_type` being coerced to "other_acc" (see `na_acc_type_to_other()`)
  # 
  # (b) In acc_lookup, `fasta_name` <- `names(fasta)`
  #   (see also `read_fasta()`)
  # =======================================================================================
  # Normally handled `acc_type`s: 
  # (1) at `acc_type` %in% c("uniprot_acc", "uniprot_id", "refseq_acc")
  #   read_fasta(">([^ ]+?) .*", ...) 
  #     -> fasta_name = "sp|Q9H553|ALG2_HUMAN" 
  #       -> "Q9H553" at acc_type == uniprot_acc
  #       -> "ALG2_HUMAN" at acc_type == uniprot_id
  #     -> NP_001254479 titin isoform... 
  #       -> "NP_001254479" at acc_type == refseq_acc
  #    
  # Two special cases: 
  # (a) if the proteoQ fasta does not include all the entries in the database-search fasta, 
  #     such `df$prot_acc` will not be in `acc_lookup` and 
  #     annotated as NA under columns fasta_name, species etc. after data joining;
  #     (possible coercion: NA fasta_name <- prot_acc)
  # 
  # (b) only at "acc_type == other_acc", 
  #     if different parsing in making `df$prot_acc` and `names(fasta)`, 
  #     such `df$prot_acc` will be annotated as NA in species etc.
  # 
  #     when making `acc_lookup`, organism etc. set to `NA_character_`, not "other_acc";
  #     this make consistent NA organism after after data joining for `df$prot_acc`s not 
  #     found in `acc_lookup`. 
  # =======================================================================================
  # Protein attributes:
  #   `fasta_name` kept in `df` for matching with `fasta_db`
  # =======================================================================================
  
  
  # add column `acc_type`; note that `df` not saved
  if (!nrow(df))
    warning("Zero row of input data.")
  
  df <- df %>% 
    parse_acc() %>% 
    dplyr::filter(nchar(prot_acc) > 0L)

  # load fasta_db
  fasta_db <- load_fasta(fasta) %>% 
    `names<-`(gsub("\\.[0-9]*$", "", names(.)))
  
  if (!length(fasta_db)) 
    warning("FASTA databases are empty.")

  # accession lookups by acc_type's
  df_splits <- df %>% 
    split(.$acc_type, drop = TRUE)
  
  acc_lookup <- lapply(df_splits, hparse_fasta, 
                       fasta_db, pat_uni_acc, pat_uni_id, pat_ref_acc, 
                       fasta, warns) %>% 
    dplyr::bind_rows()

  # --- subset and annotate fasta_db with entries in acc_lookup ---
  fasta_db <- local({
    fasta_db <- fasta_db %>% 
      .[names(.) %in% unique(acc_lookup$fasta_name)] %>% 
      .[! duplicated(names(.))]
    
    if (!length(fasta_db)) 
      stop("No fasta entries match protein accessions; probably wrong fasta file.", 
           call. = FALSE)

    fasta_db <- lapply(unique(df$acc_type), add_prot_desc, 
                       acc_lookup = acc_lookup, 
                       fasta_db = fasta_db) %>% 
      do.call(`c`, .) %>% 
      .[!duplicated(names(.))]
    
    save(fasta_db, file = file.path(dat_dir, "fasta_db.rda"))
    
    fasta_db
  }) 

  # -- add columns prot_desc, prot_mass, prot_len ---
  acc_lookup <- local({
    more <- dplyr::bind_cols(
      prot_acc = unname(purrr::map_chr(fasta_db, attr, "prot_acc")), 
      prot_desc = unname(purrr::map_chr(fasta_db, attr, "prot_desc")), 
      prot_mass = unname(purrr::map_dbl(fasta_db, ~ calc_avgpep(.x))), 
      prot_len = unname(nchar(fasta_db)), 
    ) %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::mutate(prot_mass = round(prot_mass, digits = 0))
    
    dplyr::left_join(acc_lookup, more, by = "prot_acc")
  })
  
  # -- add columns gene, organism, species, entrez ---
  df_splits <- df %>% 
    split(.$acc_type, drop = TRUE)
  
  acc_lookup <- purrr::map(df_splits, ~ {
    acc_type <- .x[["acc_type"]][1]
    
    if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
      genes <- local({
        genes <- acc_lookup %>% dplyr::select(prot_acc, prot_desc)
        
        na_genes <- genes %>% 
          dplyr::filter(!grepl("GN=", .$prot_desc)) %>% 
          dplyr::mutate(gene = NA_character_) 
        
        genes <- genes %>% 
          dplyr::filter(grepl("GN=", .$prot_desc)) %>% 
          dplyr::mutate(gene = gsub("^.*GN=(\\S+)\\s*.*", "\\1", prot_desc)) %>% 
          dplyr::bind_rows(., na_genes) %>% 
          na_genes_by_acc()
      })
      
      acc_lookup <- local({
        na_org <- genes %>% 
          dplyr::filter(!grepl("OS=", .$prot_desc)) %>% 
          dplyr::mutate(organism = NA_character_)
        
        # adds gene and organism
        acc_lookup <- genes %>% 
          dplyr::filter(grepl("OS=", .$prot_desc)) %>% 
          dplyr::mutate(organism = gsub("^.*OS=(.*?)=.*$", "\\1", prot_desc)) %>% 
          dplyr::mutate(organism = gsub("\\s\\S*$", "", organism)) %>% 
          dplyr::bind_rows(., na_org) %>% 
          dplyr::select(-prot_desc) %>% 
          dplyr::right_join(acc_lookup, by = "prot_acc")

        # adds species
        acc_lookup <- acc_lookup %>% 
          dplyr::mutate(species = my_lookup[.$organism]) %>% 
          na_species_by_org() 
      })
      
      if (!nrow(acc_lookup)) {
        warning("\nNo matches in protein accessions at `acc_type = ", acc_type, "`.\n", 
                "Fix the FASTA before running into errors", 
                call. = FALSE)
      }

      # adds entrez
      acc_lookup <- if (is.null(entrez)) 
        add_entrez(acc_lookup, acc_type, warns)
      else 
        add_custom_entrez(acc_lookup, entrez, acc_type)
      
      acc_lookup <- acc_lookup %>% 
        dplyr::filter(!is.na(.[["uniprot_acc"]]))
    } 
    else if (acc_type == "refseq_acc") {
      # adds organism and species
      acc_lookup <- acc_lookup %>% 
        dplyr::filter(grepl("^[XNY]{1}[MRP]{1}_", prot_acc)) %>% 
        dplyr::mutate(organism = gsub("^.*\\[(.*)\\]\\.*.*", "\\1", prot_desc)) %>% 
        dplyr::mutate(species = my_lookup[.$organism]) %>% 
        na_species_by_org()
      
      if (!nrow(acc_lookup)) {
        warning("\nNo matches in protein accessions at `acc_type = ", acc_type, "`.\n", 
                "Fix the FASTA before running into an error such as ", 
                "no species in gene name annoation.", 
                call. = FALSE)
      }
      
      # adds entrez
      acc_lookup <- if (is.null(entrez)) 
        add_entrez(acc_lookup, acc_type, warns)
      else 
        add_custom_entrez(acc_lookup, entrez, acc_type)
      
      # separately column gene
      acc_lookup <- add_refseq_gene(acc_lookup, acc_type)
      
      acc_lookup <- acc_lookup %>% 
        dplyr::filter(!is.na(.[["refseq_acc"]]))
    } 
    else {
      acc_lookup <- acc_lookup %>% 
        dplyr::filter(acc_type == "other_acc") %>% 
        dplyr::mutate(gene = NA_character_, 
                      organism = NA_character_, 
                      species = NA_character_, 
                      entrez = NA_integer_) 
    }
    
    acc_lookup <- acc_lookup %>% 
      dplyr::select(c("prot_acc", "gene", "organism", "prot_desc", 
                      "prot_mass", "prot_len", "fasta_name", "uniprot_acc", 
                      "uniprot_id", "refseq_acc", "other_acc", "acc_type", 
                      "entrez", "species")) %>% 
      na_genes_by_acc()
  }) %>% 
    dplyr::bind_rows() %>% 
    reloc_col_after("acc_type", "species")

  save(acc_lookup, file = file.path(dat_dir, "acc_lookup.rda"))
  
  invisible(acc_lookup)
}


#' Helper of \link{parse_fasta} by accession type.
#'
#' NA genes are replaced with prot_acc.
#' 
#' @param df_sub An input data frame with a single type of protein accession.
#' @param fasta A fasta database.
#' @inheritParams parse_fasta
hparse_fasta <- function (df_sub, fasta_db, pat_uni_acc, pat_uni_id, pat_ref_acc, 
                          fasta, warns) 
{
  acc_type <- df_sub[["acc_type"]][1]
  
  if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
    acc_lookup <- names(fasta_db) %>% 
      gsub("^..\\|", "", .) %>% 
      stringr::str_split("\\|", simplify = TRUE) %>% 
      data.frame(check.names = FALSE) %>% 
      dplyr::select(1:2) %>% 
      `names<-`(c("uniprot_acc", "uniprot_id")) %>% 
      dplyr::bind_cols(data.frame(fasta_name = names(fasta_db)), .) %>% 
      dplyr::filter(grepl(pat_uni_acc, uniprot_acc), 
                    grepl(pat_uni_id, uniprot_id)) %>%     
      dplyr::filter(.[[acc_type]] %in% unique(df_sub$prot_acc)) %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::mutate(refseq_acc = NA_character_, 
                    other_acc = NA_character_, 
                    acc_type = acc_type, 
                    prot_acc = !!rlang::sym(acc_type))
  } 
  else if (acc_type == "refseq_acc") {
    unique_accs <- unique(df_sub$prot_acc)
    unique_accs <- gsub("\\.[0-9]+$", "", unique_accs)
    
    acc_lookup <- tibble::tibble(fasta_name = names(fasta_db), 
                                 refseq_acc = names(fasta_db)) %>% 
      dplyr::filter(grepl(pat_ref_acc, refseq_acc)) %>% 
      dplyr::filter(.[[acc_type]] %in% unique_accs) %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::mutate(uniprot_acc = NA_character_, 
                    uniprot_id = NA_character_, 
                    other_acc = NA_character_, 
                    acc_type = acc_type, 
                    prot_acc = !!rlang::sym(acc_type)) %>% 
      dplyr::select(c("fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                      "other_acc", "acc_type", "prot_acc"))
  } 
  else if (acc_type == "other_acc") {
    acc_lookup <- tibble::tibble(fasta_name = names(fasta_db), 
                                 other_acc = names(fasta_db)) %>% 
      dplyr::filter(other_acc %in% unique(df_sub$prot_acc)) %>% 
      dplyr::filter(!duplicated(other_acc)) %>% 
      dplyr::mutate(uniprot_acc = NA_character_, 
                    uniprot_id = NA_character_, 
                    refseq_acc = NA_character_, 
                    acc_type = acc_type, 
                    prot_acc = other_acc) %>% 
      dplyr::select(c("fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                      "other_acc", "acc_type", "prot_acc")) 
  }
  
  n_fasta <- length(fasta_db)
  n_lookup <- nrow(acc_lookup)
  perc <- n_lookup/n_fasta 
  
  message("At acc_type = \"", acc_type, "\":\n", 
          " the number of entries in FASTA: ", n_fasta, "\n", 
          " the number of entries matched: ", n_lookup, "\n", 
          " percent annotated: ", format(round(perc, 3), nsmall = 3))
  
  if (acc_type != "other_acc" && perc < 1e-5 && warns) {
    warning("\n\n============================================================\n", 
            "PLEASE READ:\n\n", 
            "Almost NO protein accessions being annotated.\n", 
            # paste(fasta, collapse = "\n"), " is `", perc, "`.\n", 
            "The fasta database(s) are probably incorrect.", 
            "\n============================================================\n", 
            call. = FALSE)
  } 
  else if (acc_type != "other_acc" && perc < .1 && warns) {
    warning("\n------------------------------------------------------------", 
            "\nDouble check the choice of fasta(s) for probable mismatches.\n", 
            "\nA common source of mistakes: ", 
            "\n choose Uniprot with proteoQ", 
            "\n but used Refseq in database searches, or vice versa.", 
            "\n------------------------------------------------------------", 
            call. = FALSE)
  }

  invisible(acc_lookup)
}


#' Helper of \link{parse_fasta}.
#' 
#' Coerces NA values in genes by the value of accessions.
#' 
#' @param acc_type The type of protein accession.
na_genes_by_acc <- function(acc_lookup) 
{
  if (nrow(acc_lookup) && 
      "prot_acc" %in% names(acc_lookup) && 
      "gene" %in% names(acc_lookup)) 
  {
    na_genes <- (is.na(acc_lookup$gene)) | (stringr::str_length(acc_lookup$gene) == 0)
    acc_lookup$gene <- as.character(acc_lookup$gene)
    acc_lookup$gene[na_genes] <- as.character(acc_lookup$prot_acc[na_genes])
  }
  
  acc_lookup
}


#' Helper of \link{parse_fasta}.
#' 
#' Coerces NA values in species by the value of organism.
#' 
#' @param acc_type The type of protein accession.
na_species_by_org <- function(acc_lookup) 
{
  if (nrow(acc_lookup) && 
      "organism" %in% names(acc_lookup) && 
      "species" %in% names(acc_lookup)) 
  {
    na_species <- 
      (is.na(acc_lookup$species)) | (stringr::str_length(acc_lookup$species) == 0)
    acc_lookup$species <- as.character(acc_lookup$species)
    acc_lookup$species[na_species] <- as.character(acc_lookup$organism[na_species])
  }
  
  acc_lookup
}


#' Helper of \link{parse_fasta}.
#' 
#' Adds attributes to fasta entries.
#' 
#' @param pat_acc A regular expression of protein accession.
#' @inheritParams add_prot_desc
add_fasta_attrs <- function (pat_acc, acc_type, acc_lookup, fasta_db) 
{
  fasta_names_sub <- acc_lookup %>% 
    dplyr::filter(!is.na(!!rlang::sym(acc_type))) %>% 
    .[["fasta_name"]] %>% 
    unique()
  
  fasta_db_sub <- fasta_db %>% 
    .[names(.) %in% fasta_names_sub]
  
  # add attributes
  prot_acc <- if(is.null(pat_acc))
    names(fasta_db_sub)
  else
    gsub(pat_acc, "\\1", names(fasta_db_sub))
  
  fasta_db_sub <- purrr::map2(fasta_db_sub, prot_acc, function (fa, pr) {
    attr(fa, "prot_acc") <- pr
    fa
  })
  
  prot_desc <- fasta_db_sub %>%
    purrr::map(attr, "header") %>% 
    gsub("^.*\\|.*\\s+?", "", .)
  
  fasta_db_sub <- purrr::map2(fasta_db_sub, prot_desc, function (fa, pr) {
    attr(fa, "prot_desc") <- pr
    fa
  })
  
  invisible(fasta_db_sub)
}


#' Helper of \link{parse_fasta}.
#' 
#' Adds \code{prot_desc} by \code{acc_type}.
#' 
#' @param acc_type The type of protein accession.
#' @param acc_lookup A look-up table of protein accessions.
#' @param fasta_db A database of fasta.
add_prot_desc <- function (acc_type, acc_lookup, fasta_db) 
{
  pat_acc <- if (acc_type == "uniprot_acc")
    "^..\\|(.*)\\|.*$"
  else if (acc_type == "uniprot_id")
    "^.*\\|(.*)$"
  else if (acc_type == "refseq_acc")
    "^(.*)\\.[0-9]*$"
  else
    NULL
  
  add_fasta_attrs(pat_acc, acc_type, acc_lookup, fasta_db) 
}


#' Adds protein annotation
#'
#' \code{annotPrn} cross-referencing proteins among \code{uniprot_acc},
#' \code{uniprot_id}, \code{refseq} and \code{entrez}. 
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @import dplyr purrr stringr 
#' @importFrom magrittr %>% %$% %T>% 
annotPrn <- function (df, fasta, entrez) 
{
  na_fasta_name_by_prot_acc <- function(df) 
  {
    if (nrow(df) && ("fasta_name" %in% names(df))) {
      na_fasta_name <- 
        (is.na(df$fasta_name)) | (stringr::str_length(df$fasta_name) == 0)
      
      df$fasta_name[na_fasta_name] <- df$prot_acc[na_fasta_name]
    }
    
    invisible(df)
  }
  
  na_acc_type_to_other <- function(df) 
  {
    if (nrow(df) && "acc_type" %in% names(df)) {
      na_acc_type <- 
        (is.na(df$acc_type)) | (stringr::str_length(df$acc_type) == 0)
      
      df$acc_type[na_acc_type] <- "other_acc"
    }
    
    invisible(df)
  }
  
  # currently with MSGF+ outputs at UniProt DB
  uniboth <- grepl("^..\\|[[:alnum:]]+\\|[^_]+_[A-Z]{1,10}$", df[["prot_acc"]])
  
  if (sum(uniboth)) {
    df[["prot_acc"]][uniboth] <- 
      gsub("^..\\|([[:alnum:]]+)\\|[^_]+_[A-Z]{1,10}$", "\\1", 
           df[["prot_acc"]][uniboth])
  }

  acc_lookup <- parse_fasta(df, fasta, entrez)
  
  acc_lookup <- dplyr::bind_cols(
    acc_lookup %>% 
      dplyr::select(prot_acc), 
    acc_lookup %>% 
      dplyr::select(-which(names(.) %in% names(df))) %>% 
      dplyr::select(-"organism"))
  
  # (1) multiple uniprot_acc(s) can share the same entrez id
  # prot_acc    gene      acc_type   uniprot_acc species entrez
  # A1AT1_MOUSE Serpina1a uniprot_id P07758      mouse    20703
  # A1AT4_MOUSE Serpina1d uniprot_id Q00897      mouse    20703
  
  # (2) same uniprot_acc can have different entrez ids
  # From	To
  # P02088	100503605
  # P02088	101488143
  # P02088	15129	
  
  # (3) ok that some uniprot_accs(s) have no corresponding entrez id
  
  # NA values after joining, for `prot_cc` without corresponding entries in fasta
  # except for the coercion to "other_acc" for `acc_type`
  # [1] "gene"        "prot_len"    "fasta_name"  "uniprot_acc" "uniprot_id" 
  # [7] "refseq_acc"  "other_acc"   "entrez"      "species"     "acc_type"   
  
  df <- df %>% 
    dplyr::mutate(psm_index = row_number()) %>% 
    dplyr::left_join(acc_lookup, by = "prot_acc") %>% 
    dplyr::filter(!duplicated(psm_index)) %>% 
    na_acc_type_to_other() %>% 
    # na_fasta_name_by_prot_acc() %>% 
    dplyr::select(-psm_index) %>% 
    # dplyr::ungroup() %>% 
    reloc_col_after("prot_mass", "prot_acc")
  
  invisible(df)
}


#' Adds kinase annotation
#'
#' @inheritParams info_anal
#' @param acc_type Character string; the type of protein accessions in one of
#'   c("refseq_acc", "uniprot_acc", "uniprot_id")
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
annotKin <- function (df, acc_type) 
{
  if (is.null(df)) 
    return(df)
  
  stopifnot ("prot_acc" %in% names(df))
  
  data(package = "proteoQ", kinase_lookup, envir = environment())
  
  if (!acc_type %in% names(kinase_lookup)) {
    df <- df %>% 
      dplyr::mutate(kin_attr = FALSE, 
                    kin_class = NA_character_, 
                    kin_order = NA_integer_)
    
    return(df)
  }
  
  lookup <- kinase_lookup %>% 
    dplyr::select(acc_type, kin_attr, kin_class, kin_order) %>%
    dplyr::filter(!duplicated(.)) %>% 
    dplyr::filter(!is.na(.[[acc_type]])) %>% 
    dplyr::mutate(!!acc_type := as.character(!!rlang::sym(acc_type)))
  
  df <- df %>% 
    dplyr::left_join(lookup, by = c("prot_acc" = acc_type))
  
  df$kin_attr[is.na(df$kin_attr)] <- FALSE
  
  invisible(df)
}


#' Saves the arguments in a function call
#'
#' @param call_pars Language.
#' @param fn The name of function being saved.
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
save_call <- function(call_pars, fn) 
{
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)
	call_pars[names(call_pars) == "..."] <- NULL
	save(call_pars, file = file.path(dat_dir, "Calls", paste0(fn, ".rda")))
}


#' Matches formulas to those in calls to pepSig or prnSig
#'
#' @param formulas Language; the formulas in linear modeling.
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
match_fmls <- function(formulas) 
{
  dat_dir <- get_gl_dat_dir()
  fml_file <-  file.path(dat_dir, "Calls/pepSig_formulas.rda")
  
  if (file.exists(fml_file)) 
    load(file = fml_file)
  else 
    stop("Run `pepSig()` first.")
  
  fml_chr <- formulas %>%
    as.character() %>%
    gsub("\\s+", "", .)
  
  pepSig_chr <- pepSig_formulas %>%
    purrr::map(~ .[is_call(.)]) %>%
    as.character() %>%
    gsub("\\s+", "", .)
  
  ok <- purrr::map_lgl(fml_chr, ~ . %in% pepSig_chr)
  
  if (!all(ok))
    stop("Formula match failed: ", formulas[[which(!ok)]],
         " not found in the latest call to 'pepSig(...)'.", 
         call. = FALSE)
}


#' Matches the arg to anal_prnGSPA
#'
#' @param call_rda the name of a rda.
#' @param arg Argument to be matched.
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
match_call_arg <- function (call_rda = "foo", arg = "scale_log2r") 
{
  dat_dir <- get_gl_dat_dir()
  
  call_rda <- rlang::as_string(rlang::enexpr(call_rda))
  arg <- rlang::as_string(rlang::enexpr(arg))
  
  rda <- paste0(call_rda, ".rda")
  file <- file.path(dat_dir, "Calls", rda)
  
  if (!file.exists(file)) 
    stop(rda, " not found.")
  
  load(file = file)
  
  if (!arg %in% names(call_pars)) 
    stop(arg, " not found in the latest call to ", call_rda, 
         call. = FALSE)
  
  call_pars[[arg]]
}


#' Matches gset_nms to prnGSPA.
#' 
#' Not currently used.
#' 
#' @inheritParams prnGSPA
match_gset_nms <- function (gset_nms = NULL) 
{
  dat_dir <- get_gl_dat_dir()
  file <- file.path(dat_dir, "Calls/anal_prnGSPA.rda")
  
  if (is.null(gset_nms)) {
    if (file.exists(file)) {
      gset_nms <- match_call_arg(anal_prnGSPA, gset_nms)
    } else {
      stop("`gset_nms` is NULL and ", file, " not found.", 
           call. = FALSE)
    }
  } 
  else {
    if (file.exists(file)) {
      # gset_nms <- gset_nms %>% .[. %in% match_call_arg(anal_prnGSPA, gset_nms)]
      gset_nms <- gset_nms %>% .[. %in% match_call_arg(anal_prnGSPA, .)]
    } else {
      warning("File `", file, "`` not available for `gset_nms` matches.", 
              "\nUse the `gset_nms` value as it.", 
              call. = FALSE)
    }
  }
  
  if (is.null(gset_nms)) {
    stop ("The `gset_nms` is NULL after matching parameters to the latest `prnGSPA(...)`.", 
          "\n\tConsider providing explicitly the `gset_nms`: ", 
          "\n\t\t`gset_nms = \"go_sets\"` or ", 
          "\n\t\t`gset_nms = \"~/proteoQ/dbs/go_hs.rds\"` ...",
          call. = FALSE)
  }
  
  if (!length(gset_nms)) {
    stop ("The `gset_nms` is EMPTY after parameter matching.", 
          "\n\tCheck the values of `gset_nms` in the latest call to `prnGSPA(...)`.", 
          call. = FALSE)
  }
  
  return(gset_nms)
}


#' A match to a global logical variable
#' 
#' @param var Character string representation of a variable.
#' @param val The value of \code{var} before matching.
#' @examples
#' \donttest{
#' foo <- function(scale_log2r = FALSE) {
#'   match_logi_gv("scale_log2r", scale_log2r)
#' }
#' 
#' scale_log2r <- TRUE
#' foo()
#' }
match_logi_gv <- function(var, val) 
{
  oval <- val

  ok <- exists(var, envir = .GlobalEnv)
  gvar <- if (ok) get(var, envir = .GlobalEnv) else NULL

  if (is.null(gvar))
    val
  else {
    if (!is.logical(gvar))
      stop("Argument \"", var, "\" is not logical.")

    if ((!is.na(gvar)) && (gvar != oval)) {
      warning("Coerce ", var, " to ", gvar, 
              " after matching to the global setting.", 
              call. = FALSE)
    }

    gvar
  }
}


#' Matches scale_log2r 
#'
#' \code{match_prnSig_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in the most recent prnSig at a given impute_na status
#' 
#' @inheritParams prnHM
match_prnSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) 
{
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  mscale_log2r <- if (impute_na) 
    match_call_arg(prnSig_impTRUE, scale_log2r)
  else 
    match_call_arg(prnSig_impFALSE, scale_log2r)
  
  if ((!is.na(mscale_log2r)) && (scale_log2r != mscale_log2r)) {
    warning("scale_log2r = ", mscale_log2r, 
            " after matching to prnSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }

  mscale_log2r
}


#' Matches scale_log2r 
#'
#' \code{match_pepSig_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in the most recent pepSig at a given impute_na status.
#' @inheritParams prnHM
match_pepSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) 
{
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  mscale_log2r <- if (impute_na) 
    match_call_arg(pepSig_impTRUE, scale_log2r)
  else 
    match_call_arg(pepSig_impFALSE, scale_log2r)
  
  if ((!is.na(mscale_log2r)) && (scale_log2r != mscale_log2r)) {
    warning("scale_log2r = ", mscale_log2r, 
            " after matching to pepSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }

  mscale_log2r
}


#' Replaces NA genes (not currently used)
#' 
#' @inheritParams info_anal
#' @inheritParams annotKin
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
replace_na_genes <- function(df, acc_type) 
{
  acc_type <- tolower(acc_type)
  
  if (acc_type == "refseq_acc") {
    na_gene <- (is.na(df[, c("gene")])) | (str_length(df$gene) == 0)
    df$gene <- as.character(df$gene)
    df$gene[na_gene] <- as.character(df$prot_acc[na_gene])
  } 
  else if (acc_type == "uniprot_id") {
    temp <- data.frame(do.call('rbind', 
                               strsplit(as.character(df$prot_desc), 
                                        'GN=', 
                                        fixed = TRUE))) %>%
      dplyr::select(2) %>%
      `colnames<-`("gene") %>%
      dplyr::mutate(gene = gsub("PE\\=.*", "", gene)) %>%
      dplyr::mutate(gene = gsub("\\s+.*", "", gene))
    
    na_gene <- is.na(df$gene)
    df[na_gene, c("gene")] <- temp[na_gene, c("gene")]
    
    rm(list = c("temp"))
    
    df <- df %>%
      dplyr::mutate(gene = gsub(".*\\|", "", gene))
  }
  
  invisible(df)
}


#' Replaces NA genes
#' 
#' @param df A data frame.
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
na_genes_by_acc <- function(df) 
{
  stopifnot("prot_acc" %in% names(df), "gene" %in% names(df))
  
  na_genes <- (is.na(df$gene)) | (stringr::str_length(df$gene) == 0)
  df$gene <- as.character(df$gene)
  df$gene[na_genes] <- as.character(df$prot_acc[na_genes])
  
  invisible(df)
}


#' Finds peptide start and end positions
#'
#' \code{find_pep_pos} finds the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}.
#' 
#' @param fasta_name The fasta name
#' @param pep_seq Peptide sequence
#' @param fasta The database of fasta
#' @import dplyr purrr stringr tidyr
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' 
#' @examples
#' \donttest{
#' fasta_name <- "NP_612429"
#' pep_seq <- "ADVSLPSMQGDLK"
#' fasta <- read_fasta("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta")
#' x <- find_pep_pos(fasta_name, pep_seq, fasta)
#' 
#' pep_seq <- "ADVSLPSMQGDL"
#' x <- find_pep_pos(fasta_name, pep_seq, fasta)
#' }
find_pep_pos <- function (fasta_name, pep_seq, fasta) 
{
  if (is.na(fasta_name)) {
    out <- cbind(pep_seq = pep_seq, pep_res_before = NA_character_, 
                 start = NA_integer_, end = NA_integer_, 
                 pep_res_after = NA_character_, fasta_name = fasta_name, 
                 pep_istryptic = NA, pep_semitryptic = NA_character_)
    
    return(out)
  }

  this_fasta <- fasta[names(fasta) == fasta_name]
  len_fasta <- length(this_fasta)
  pep_seq <- as.character(pep_seq)
  
  if (len_fasta == 1L) {
    pep_pos_all <- stringr::str_locate_all(this_fasta, pattern = pep_seq)
    pep_pos_all <- pep_pos_all[[1]]
    
    # can have no matches if the fasta is different to what are used in 
    #  database searches; i.e., the fasta sequence may be altered.

    n_rows <- nrow(pep_pos_all)
    
    if (!n_rows) {
      warning(pep_seq, " not found in ", fasta_name, ".\n", 
              "Probably different fastas between database searches and proteoQ.", 
              call. = FALSE)
      
      out <- cbind(pep_seq = pep_seq, pep_res_before = NA_character_, 
                   start = NA_integer_, end = NA_integer_, 
                   pep_res_after = NA_character_, fasta_name = fasta_name, 
                   pep_istryptic = NA, pep_semitryptic = NA_character_)

      return(out)
    }
    
    cts <- nts <- logical(n_rows)

    for (i in seq_len(n_rows)) {
      pep_pos_i <- pep_pos_all[i, ]
      
      pos_bf <- pep_pos_i[1] - 1L
      pos_af <- pep_pos_i[2] + 1L
      
      pep_res_before <- stringr::str_sub(this_fasta, pos_bf, pos_bf)
      pep_res_after <- stringr::str_sub(this_fasta, pos_af, pos_af)
      
      # Mascot specialty of "X" residues (see also aa_residues.rda)
      # prot_acc: "XP_003960355", original "QERFCQXK" becomes "QERFCQVK"
      
      if (is.na(pep_res_before) || is.na(pep_res_after)) {
        out <- cbind(pep_seq = pep_seq, pep_res_before = NA_character_, 
                     start = NA_integer_, end = NA_integer_, 
                     pep_res_after = NA_character_, fasta_name = fasta_name, 
                     pep_istryptic = NA, pep_semitryptic = NA_character_)
        
        return(out)
      }
      
      if (nchar(pep_res_before) == 0L) pep_res_before <- "-"
      if (nchar(pep_res_after) == 0L) pep_res_after <- "-"

      # N-term M: 
      #   pep_res_before == "M" && pep_pos_i[1] == 2L
      # Not N-term M: 
      #   pep_res_before != "M" || pep_pos_i[1] != 2L
      
      ok_ct <- grepl("[KR]$", pep_seq) || pep_res_after == "-"
      is_nt <- pep_res_before %in% c("K", "R", "-")
      is_nm <- pep_res_before == "M" && pep_pos_i[1] == 2L
      ok_nt <- is_nt || is_nm
      
      if (ok_ct && ok_nt) {
        out <- cbind(pep_seq = pep_seq, 
                     pep_res_before = pep_res_before, 
                     start = pep_pos_i[[1]], end = pep_pos_i[[2]],  
                     pep_res_after = pep_res_after, fasta_name = fasta_name, 
                     pep_istryptic = TRUE, pep_semitryptic = "both")
        
        return(out)
      }
      else {
        cts[[i]] <- ok_ct
        nts[[i]] <- ok_nt
      }
    }
    
    ## no full-tryptic matches
    semi_tryp <- "none"
    
    for (k in seq_along(nts)) {
      if (nts[[k]]) {
        semi_tryp <- "nterm"
        break
      }
      else if (cts[[k]]) {
        semi_tryp <- "cterm"
        break
      }
    }
    
    pep_pos_i <- pep_pos_all[k, ]

    out <- cbind(pep_seq = pep_seq, 
                 pep_res_before = pep_res_before, 
                 start = pep_pos_i[[1]], end = pep_pos_i[[2]], 
                 pep_res_after = pep_res_after, fasta_name = fasta_name, 
                 pep_istryptic = FALSE, pep_semitryptic = semi_tryp)
  } 
  else  if (len_fasta == 0L) {
    out <- cbind(pep_seq = pep_seq, 
                 pep_res_before = NA_character_, 
                 start = NA_integer_, 
                 end = NA_integer_, 
                 pep_res_after = NA_character_, 
                 fasta_name = fasta_name, 
                 pep_istryptic = NA, pep_semitryptic = NA_character_)
  }
  else {
    stop("More than one matched fasta_name.")
  }

  invisible(out)
}


#' Annotation of peptide positions and adjacent amino acid residues
#'
#' \code{annotPeppos} annotates the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}. It also annotates the
#' preceding and the following AA residues.
#' 
#' @param df A data frame
#' @import dplyr purrr stringr tidyr
#' @importFrom magrittr %>% %T>% %$% %<>% 
annotPeppos <- function (df) 
{
  stopifnot(all(c("fasta_name", "pep_seq") %in% names(df)))
  
  file <- file.path(dat_dir, "fasta_db.rda")
  
  if (!file.exists(file)) 
    stop("File not found: ", file, call. = FALSE)
  
  load(file = file)

  # ok cases that same `pep_seq` but different `prot_acc`
  # (K)	MENGQSTAAK	(L) NP_510965
  # (-)	MENGQSTAAK	(L) NP_001129505
  
  df <- df %>% 
    dplyr::mutate(pep_prn = paste(pep_seq, fasta_name, sep = "@"))
  
  df_pep_prn <- df %>% 
    dplyr::filter(!duplicated(pep_prn)) %>% 
    dplyr::select(c("pep_seq", "fasta_name")) 
  
  fasta_db_sub <- fasta_db %>% 
    .[names(.) %in% unique(df_pep_prn$fasta_name)]

  pep_pos_all <- suppressWarnings(
    purrr::map2(as.list(df_pep_prn$fasta_name), 
                as.list(df_pep_prn$pep_seq), 
                find_pep_pos, 
                fasta_db_sub)
  ) %>% 
    do.call(rbind, .) %>% 
    `colnames<-`(c("pep_seq", "pep_res_before", "pep_start", 
                   "pep_end", "pep_res_after", 
                   "fasta_name", "pep_istryptic", "pep_semitryptic")) %>% 
    data.frame(check.names = FALSE) %>% 
    tidyr::unite(pep_prn, pep_seq, fasta_name, sep = "@", remove = TRUE) %>% 
    dplyr::mutate(pep_start = as.integer(pep_start), 
                  pep_end = as.integer(pep_end))
  
  df$pep_res_before <- NULL
  df$pep_res_after <- NULL
  df$pep_start <- NULL
  df$pep_end <- NULL
  nms <- names(df)
  
  if ("pep_res_before" %in% nms) pep_pos_all$pep_res_before <- NULL
  if ("pep_res_after" %in% nms) pep_pos_all$pep_res_after <- NULL
  if ("pep_start" %in% nms) pep_pos_all$pep_start <- NULL
  if ("pep_end" %in% nms) pep_pos_all$pep_end <- NULL
  
  df %>% 
    dplyr::left_join(pep_pos_all, by = "pep_prn") %>% 
    dplyr::select(-pep_prn) %>% 
    reloc_col_before("pep_end", "pep_res_before") %>% 
    reloc_col_before("pep_start", "pep_end")
}


#' Subsets fasta by accession type
#' 
#' Not yet used.
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @inheritParams annotKin
#' @import dplyr purrr  
#' @importFrom magrittr %>% %T>% %$% %<>% 
subset_fasta <- function (df, fasta, acc_type) 
{
  stopifnot("prot_acc" %in% names(df))
  
  if (! acc_type %in% c("refseq_acc", "uniprot_id", "uniprot_acc")) 
    stop("The type of protein accesion needs to one of ", 
         "\'uniprot_id\', \'uniprot_acc\' or \'refseq_acc\'.",
         call. = FALSE)
  
  fasta <- purrr::map(fasta, ~ read_fasta(.x)) %>% do.call(`c`, .)
  
  if (acc_type == "uniprot_id") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|.*\\|(.*)$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } 
  else if (acc_type == "uniprot_acc") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } 
  else if (acc_type == "refseq_acc") {
    fasta <- fasta %>% 
      .[names(.) %in% unique(df$prot_acc)]
  }    
}


#' Finds peptide abundance indexes.
#'
#' In relative to tryptic peptides.
#'
#' @param max_len The maximum length of peptide sequence for consideration.
#' @inheritParams info_anal
#' @inheritParams find_shared_prots
add_prot_icover <- function (df, id = "gene", pep_id = "pep_seq", 
                             max_len = 500L) 
{
  dat_dir <- get_gl_dat_dir()
  
  ok <- tryCatch(load(file = file.path(dat_dir, "fasta_db.rda")),
                 error = function(e) "e")
  
  if (ok != "fasta_db") 
    stop("`fasta_db.rda` not found under ", dat_dir, ".", 
         call. = FALSE)
  
  fasta_db <- fasta_db %>% .[!duplicated(names(.))]
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  if (id == "gene") {
    gn_rollup <- TRUE
    id <- "prot_acc"
  } 
  else {
    gn_rollup <- FALSE
  }
  
  min_len <- min(nchar(df[[pep_id]])) # - 1L
  max_len2 <- max_len # - 1L
  
  zero_npeps <- df %>% 
    dplyr::mutate(pep_istryptic = as.logical(pep_istryptic)) %>% 
    { if ("pep_istryptic" %in% names(df)) dplyr::filter(., pep_istryptic) else . } %>% 
    dplyr::filter(pep_len <= max_len, pep_miss == 0L) %>% 
    dplyr::select(c(pep_id, "fasta_name")) %>%
    dplyr::filter(!duplicated(!!rlang::sym(pep_id))) %>% 
    dplyr::group_by(fasta_name) %>%
    dplyr::summarise(prot_n_pep0 = n())
  
  max_npeps <- fasta_db %>% 
    purrr::map_int(~ {
      .x <- gsub("([KR]{1})", paste0("\\1", "@"), .x)
      x <- stringr::str_split(.x, "@")
      
      x <- lapply(x, function (p) {
        p1 <- p[[1]]
        if (grepl("^M", p1)) c(gsub("^M", "", p1), p) else p
      })

      x %>% purrr::map_int(~ {
        len <- stringr::str_length(.x)
        sum(len >= min_len & len <= max_len2)
      })
    }) %>% 
    tibble::tibble(
      fasta_name = names(.), 
      prot_n_pepi = .,
    ) 

  max_npeps <- max_npeps %>% 
    dplyr::left_join(zero_npeps, by = "fasta_name") %>% 
    dplyr::mutate(prot_n_pep0 = ifelse(is.na(prot_n_pep0), 1L, prot_n_pep0)) %>% 
    dplyr::mutate(prot_icover = prot_n_pep0/prot_n_pepi, 
                  prot_icover = round(prot_icover, digits = 3L)) %>% 
    dplyr::select(-c("prot_n_pepi", "prot_n_pep0"))
  
  if (gn_rollup) {
    df <- local({
      tempdata <- df %>% 
        dplyr::select(fasta_name, gene) %>% 
        dplyr::filter(!duplicated(fasta_name)) %>% 
        dplyr::left_join(max_npeps, by = "fasta_name") %>% 
        dplyr::select(-fasta_name) %>% 
        dplyr::filter(!is.na(gene)) %>% 
        dplyr::group_by(gene)
      
      suppressWarnings(
        tempdata %>% 
          dplyr::summarise_all(~ max(.x, na.rm = TRUE))
      ) %>% 
        dplyr::left_join(df, ., by = "gene")
    })
  } 
  else {
    df <- df %>% 
      dplyr::left_join(max_npeps, by = "fasta_name")
  }
  
  invisible(df)
}


#' Calculates protein percent coverage
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @import dplyr purrr  
#' @importFrom magrittr %>% %T>% %$% %<>% 
calc_cover <- function(df, id) 
{
  dat_dir <- get_gl_dat_dir()
  
  stopifnot(all(c("prot_acc", "gene", "pep_start", "pep_end") %in% names(df)))
  
  df$prot_cover <- NULL
  
  load(file = file.path(dat_dir, "fasta_db.rda"))

  if (all(is.factor(df$pep_start))) {
    df$pep_start <-as.numeric(as.character(df$pep_start))
  }
    
  if (all(is.factor(df$pep_end))) {
    df$pep_end <- as.numeric(as.character(df$pep_end))
  }

  id <- rlang::as_string(rlang::enexpr(id))
  
  if (id == "gene") {
    gn_rollup <- TRUE
    id <- "prot_acc"
  } 
  else {
    gn_rollup <- FALSE
  }
  
  load(file = file.path(dat_dir, "acc_lookup.rda"))
  
  if (!length(fasta_db)) 
    stop("No fasta entries match protein accessions.", 
         "Check the correctness of fasta file(s).", 
         call. = FALSE)
  
  if (length(fasta_db) <= 200L) 
    warning("Less than 200 entries in fasta matched by protein accessions.\n", 
            "Make sure the fasta file is correct.", 
            call. = FALSE)
  
  df_sels <- df %>%
    dplyr::select(prot_acc, pep_start, pep_end) %>%
    unique() %>% 
    dplyr::left_join(acc_lookup, by = "prot_acc") %>%
    dplyr::filter(!is.na(prot_len)) %>% 
    dplyr::filter(pep_start <= prot_len, 
                  pep_end <= prot_len)
  
  if (!nrow(df_sels)) 
    stop("Probably incorrect accession types in the fasta file(s).", 
         call. = FALSE)
  
  df_sels <- df_sels %>%
    split(.[["prot_acc"]], drop = TRUE) %>%
    purrr::map(function (.x) {
      len <- .x[1, "prot_len"]
      aa_map <- rep(NA, len)
      
      for (i in 1:nrow(.x)) {
        aa_map[.x[i, ]$pep_start : .x[i, ]$pep_end] <- TRUE
      }
      
      sum(aa_map, na.rm = TRUE)/len
    } ) %>%
    do.call("rbind", .) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`("prot_cover") %>%
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::mutate(prot_cover = ifelse(prot_cover > 1, 1, prot_cover))
  
  if (gn_rollup) {
    df_sels <- local({
      df_sub <- df %>% 
        dplyr::select(prot_acc, gene) %>% 
        dplyr::filter(!duplicated(prot_acc))
      
      df_sels <- df_sub %>% 
        dplyr::left_join(df_sels, by = "prot_acc") %>% 
        dplyr::select(-prot_acc) %>% 
        dplyr::filter(!is.na(gene)) %>% 
        dplyr::group_by(gene) 
      
      df_sels <- suppressWarnings(
        df_sels %>% dplyr::summarise_all(~ max(.x, na.rm = TRUE))
      )
      
      df_sels <- df_sub %>% 
        dplyr::left_join(df_sels, by = "gene") %>% 
        dplyr::select(-gene)
    })
  }
  
  df <- df %>% 
    dplyr::mutate(index = row_number()) %>% 
    dplyr::left_join(df_sels, by = "prot_acc") %>% 
    dplyr::filter(!duplicated(index)) %>% 
    dplyr::select(-index) 
  
  if ("prot_icover" %in% names(df)) {
    df <- df %>% 
      dplyr::mutate(prot_icover = 
                      ifelse(!is.na(prot_icover), prot_icover, prot_cover), 
                    prot_icover = round(prot_icover, digits = 3L))
  }
  
  df <- df %>% 
    dplyr::mutate(prot_cover = round(prot_cover, digits = 3L)) 

  invisible(df)
}


#' Converts log2FC to linear fold changes
#' 
#' @inheritParams info_anal
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
to_linfc <- function(df) 
{
  nms <- rownames(df)
  
  df %>%
    purrr::map(~ {ifelse(.x > 0, 2^.x, -1/(2^.x))}) %>%
    data.frame(check.names = FALSE) %>%
    `rownames<-`(nms)
}


#' Removes single-value columns
#' 
#' @inheritParams setHMlims
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
rm_sglval_cols <- function (x) 
{
  sgl_val <- x %>% 
    dplyr::summarise_all(~ n_distinct(.x)) %>% 
    purrr::map(~ .x == 1) %>% 
    purrr::flatten_lgl()
  
  cols <- names(x[sgl_val])
  
  if (!purrr::is_empty(cols)) {
    warning("Columns at a single-factor level for aesthetics excluded: \n", 
            purrr::reduce(cols, paste, sep = ", "), ".\n",
            call. = FALSE)
  }

  x[, !sgl_val, drop = FALSE]
}


#' Combines data with metadata
#' 
#' @param data A data frame
#' @param metadata Another data frame 
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
cmbn_meta <- function(data, metadata) 
{
  data <- data %>% 
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::left_join(metadata, by = "Sample_ID") %>%
    dplyr::mutate_at(vars(one_of("Color", "Fill", "Shape", "Size", "Alpha")), 
                     ~ as.factor(.)) 
  
  col_nms <- data %>% 
    .[!not_all_NA(.)] %>% 
    names()
  
  if (!purrr::is_empty(col_nms)) 
    warning("columns with all-NA values for aesthetics excluded: \n", 
            purrr::reduce(col_nms, paste, sep = ", "), ".\n", 
            call. = FALSE)
  
  data %>%
    dplyr::select(which(not_all_NA(.))) %>% 
    rm_sglval_cols()
}


#' Checks file names for ggsave
#' 
#' @param filename Character string; An output file name.
gg_imgname <- function(filename) 
{
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  exts <- c("png", "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg") 
  
  if(! fn_suffix %in% exts) {
    warning("Unrecognized file extenstion: '", fn_suffix, 
            "'. Image will be saved as a '.png'.", 
            call. = FALSE)
    
    fn_suffix <- "png"
  }
  
  paste0(fn_prefix, ".", fn_suffix)
}


#' Checks file names for ggsave
#' 
#' @inheritParams info_anal
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
rm_pval_whitespace <- function(df) 
{
  df <- df %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.character) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ gsub("\\s*", "", .x) ) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ suppressWarnings(as.numeric(.x)))
}


#' Filters rows
#'
#' @param df a data frame. 
#' @param ... Arguments for \link[dplyr]{filter}
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
filters_in_call <- function (df, ...) 
{
  dots <- rlang::enexprs(...)
  nms <- names(dots)
  
  for (i in seq_along(dots)) {
    row_exprs <- dots[[nms[i]]] %>% 
      rlang::eval_bare()
    
    if (!rlang::is_list(row_exprs)) row_exprs <- list(row_exprs)
    
    row_exprs <- purrr::map(row_exprs, ~ {
      is_char <- is.character(.x)
      if (is_char) .x <- rlang::sym(.x)
      
      return(.x)
    })
    
    row_vals <- row_exprs %>% 
      purrr::map(eval_tidy, df) %>% 
      purrr::reduce(`&`, .init = 1)

    row_vals <- row_vals %>% ifelse(is.na(.), FALSE, .)

    stopifnot(is.logical(row_vals))
    
    if (sum(row_vals) == 0) 
      stop("Zero row of data available after `filter_` vararg(s).\n", 
           "Examine both (a) the values under the selected column(s) and ", 
           "(b) the supplied logical condition(s).\n", 
           "For example, values are not all NAs or FALSE under the column(s).", 
           call. = FALSE)
    
    df <- df[row_vals, , drop = FALSE]
  }
  
  return(df)
}


#' Arranges rows
#'
#' @param .df a data frame. 
#' @param .na.last The same as \link[base]{order}
#' @param ... Arguments for \link[dplyr]{arrange}
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
arrangers_in_call <- function(.df, ..., .na.last = TRUE) 
{
  dots <- rlang::enexprs(...)
  nms <- names(dots)
  
  for (i in seq_along(dots)) {
    row_orders <- dots[[nms[i]]] %>% 
      rlang::eval_bare()
    
    if (!rlang::is_list(row_orders)) 
      row_orders <- list(row_orders)
    
    row_orders <- purrr::map(row_orders, ~ {
      is_char <- is.character(.x)
      if (is_char) .x <- rlang::sym(.x)
      
      .x
    })
    
    order_call <- rlang::expr(order(!!!row_orders, na.last = !!.na.last))
    ord <- rlang::eval_tidy(order_call, .df)
    stopifnot(length(ord) == nrow(.df))

    .df <- .df[ord, , drop = FALSE]
  }
  
  return(.df)
}


#' Calculates PSM SDs
#' 
#' The standard deviations are based on samples under each TMT set and LCMS.
#' 
#' @inheritParams info_anal
#' @inheritParams standPep
#' @inheritParams channelInfo
#' @inheritParams calcPeptide
calc_sd_fcts_psm <- function (df, range_log2r = c(5, 95), range_int = c(5, 95), 
                              set_idx, injn_idx) 
{
  dat_dir <- get_gl_dat_dir()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  
  label_scheme <- label_scheme_full %>% 
    dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx) %>% 
    dplyr::mutate(tmt_nm = gsub("TMT-", "N_log2_R", TMT_Channel))
  
  label_scheme_sd <- label_scheme %>%
    dplyr::filter(!Reference, !grepl("^Empty\\.", Sample_ID))	%>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = (.$Sample_ID)))
  
  SD <- df %>%
    dplyr::select(grep("^N_log2_R|^N_I", names(.))) %>%
    dblTrim(., range_log2r, range_int) %>%
    `names<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.)))
  
  cf_SD <- SD/mean(SD %>% .[names(.) %in% label_scheme_sd$tmt_nm], na.rm = TRUE)
  
  cf_SD <- cbind.data.frame(fct = cf_SD, SD) %>%
    tibble::rownames_to_column("tmt_nm") %>%
    dplyr::mutate(tmt_nm = factor(tmt_nm, levels = label_scheme$tmt_nm)) %>%
    dplyr::arrange(tmt_nm)
}


#' Calculates CV per TMT_Set and LCMS_injection
#' 
#' @inheritParams info_anal
#' @param type Character string; the type of data for SD calculations.
calcSD_Splex <- function (df, id, type = "log2_R") 
{
  if (type == "log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^log2_R[0-9]{3}", names(.)))
  } 
  else if (type == "N_log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^N_log2_R[0-9]{3}", names(.)))
  } 
  else if (type == "Z_log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^Z_log2_R[0-9]{3}", names(.)))
  }
  
  df %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id))) %>% 
    dplyr::arrange(!!rlang::sym(id)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarise_at(vars(starts_with(type)), ~ sd(.x, na.rm = TRUE)) %>% 
    dplyr::mutate_at(vars(starts_with(type)), ~ replace_trivial_with_na(.x))
}



#' Violin plots of CV per TMT_Set and LCMS_injection
#' 
#' @param width The width of a plot.
#' @param height The height of a plot.
#' @param is_psm Logical; indicator if the data belong to a PSM table .
#' 
#' @inheritParams info_anal
#' @inheritParams purgePSM
#' @inheritParams prnCorr_logFC
#' @inheritParams calcSD_Splex
#' @import dplyr ggplot2
sd_violin <- function(df = NULL, id = NULL, filepath = NULL, 
                      width = NULL, height = NULL, 
                      type = "log2_R", adjSD = FALSE, 
                      is_psm = FALSE, col_select = NULL, col_order = NULL, 
                      theme = NULL, ...) 
{
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filepath)
  fn_prefix <- gsub("\\.[^.]*$", "", filepath)
  
  df <- df %>% 
    dplyr::select(which(not_all_NA(.)))
  
  if (!ncol(df)) {
    message("No SD columns available.")
    return (NULL)
  }
  
  dat_dir <- get_gl_dat_dir()
  
  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  col_select <- if (is.null(col_select)) 
    rlang::expr(Select)
  else 
    rlang::sym(col_select)

  col_order <- if (is.null(col_order)) 
    rlang::expr(Order)
  else 
    rlang::sym(col_order)

  if (col_select == rlang::expr(Sample_ID)) 
    stop(err_msg1, call. = FALSE)
  
  if (col_order == rlang::expr(Sample_ID)) 
    stop(err_msg1, call. = FALSE)
  
  label_scheme_full <- load_ls_group(dat_dir, label_scheme_full)
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  
  if (is.null(label_scheme_full[[col_select]])) 
    stop("Column \'", rlang::as_string(col_select), "\' does not exist.", 
         call. = FALSE)
  else if (sum(!is.na(label_scheme_full[[col_select]])) == 0) 
    stop("No samples were selected under column \'", rlang::as_string(col_select), "\'.",
         call. = FALSE)
  
  if (is.null(label_scheme_full[[col_order]])) {
    warning("Column \'", rlang::as_string(col_order), "\' does not exist.
			Samples will be arranged by the alphebatic order.", 
            call. = FALSE)
  } 
  else if (sum(!is.na(label_scheme_full[[col_order]])) == 0) {
    # warning("No orders under column \'", rlang::as_string(col_order), "\'.", call. = FALSE)
  }
  
  if (is_psm) {
    label_scheme_sub <- local({
      set_idx <- as.integer(gsub("^.*TMTset(\\d+)_.*", "\\1", filepath))
      injn_idx <- as.integer(gsub("^.*TMTset\\d+_LCMSinj(\\d+)_sd\\.png$", "\\1", filepath))

      label_scheme_full %>% 
        dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx) 
    })
  } 
  else {
    label_scheme_sub <- label_scheme
  }
  
  label_scheme_sub <- label_scheme_sub %>% 
    dplyr::mutate(new_id = paste0(TMT_Channel, " (", Sample_ID, ")")) %>% 
    dplyr::mutate(new_id = gsub("TMT-", "", new_id)) %>% 
    dplyr::select(Sample_ID, TMT_Set, new_id, !!col_select, !!col_order) %>%
    dplyr::filter(!is.na(!!col_select))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  
  if (!TMT_plex) {
    label_scheme_sub$new_id <- label_scheme_sub$Sample_ID
  }

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  if (rlang::is_missing(width)) width <- 8
  if (rlang::is_missing(height)) height <- 8
  
  ymax <- eval(dots$ymax, envir = rlang::caller_env())
  ybreaks <- eval(dots$ybreaks, envir = rlang::caller_env())
  
  flip_coord <- eval(dots$flip_coord, envir = rlang::caller_env())
  if (is.null(flip_coord)) flip_coord <- FALSE
  
  df <- df %>% 
    dplyr::filter(!duplicated(.[[id]]))
  
  if (("pep_scan_num" %in% names(df)) && (!all(is.na(df$pep_scan_num)))) {
    df <- df %>% 
      tidyr::unite(uniq_by., raw_file, pep_scan_num, sep = "@") %>% 
      dplyr::group_by(uniq_by.) %>% 
      dplyr::filter(row_number() == 1L) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-uniq_by.)
  }

  pat0 <- "[0-9]{3}[NC]{0,1}"
  
  df_sd <- df %>% 
    dplyr::select(id, grep(paste0("^sd_", type, pat0), names(.)))
  
  # all-NA first removed for finding all-NaN columns
  pat_log2_R <- paste0("^.*log2_R", pat0)

  df_sd <- df_sd %>% 
    dplyr::filter(rowSums(!is.na(.[grep(pat_log2_R, names(.))])) > 0L)
  
  if (adjSD) {
    SD <- df %>%
      dplyr::select(grep("^log2_R[0-9]{3}[NC]{0,1}|^I[0-9]{3}[NC]{0,1}", names(.))) %>%
      dblTrim(range_log2r = c(0, 100), range_int = c(0, 100), 
              type_r = "log2_R", type_int = "I")
    
    df_sd[, grep("^.*log2_R", names(df_sd))] <- 
      df_sd[, grep("^.*log2_R", names(df_sd)), drop = FALSE] %>% 
      sweep(., 2, sqrt(SD), "/")
    
    df_z <- df_sd %>% dplyr::select(grep(pat_log2_R, names(.)))
    nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
    df_z[, nan_cols] <- 0
    df_sd[, grep("^.*_log2_R[0-9]{3}", names(df_sd))] <- df_z
    
    rm(list = c("df_z", {nan_cols}, "SD"))
  }

  if (TMT_plex) {
    df_sd <- df_sd %>% 
      `names<-`(gsub("^.*log2_R", "", names(.))) 
  } 
  else {
    df_sd <- df_sd %>% 
      `names<-`(gsub("sd_log2_R000 \\((.*)\\)$", "\\1", names(.))) 
  }
  
  if (!is_psm) {
    df_sd <- df_sd %>% 
      dplyr::select(id, which(names(.) %in% label_scheme_sub$new_id))
  }

  Levels <- names(df_sd) %>% .[! . %in% id]

  if (length(Levels)) {
    df_sd <- df_sd %>%
      tidyr::gather(key = !!rlang::sym(id), value = "SD") %>%
      dplyr::rename(Channel := !!rlang::sym(id)) %>% 
      dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
      dplyr::filter(!is.na(SD))
    
    mean_sd <- df_sd %>% 
      dplyr::group_by(Channel) %>% 
      dplyr::summarise(SD = mean(SD, na.rm = TRUE)) %>% 
      dplyr::mutate(SD = round(SD, digits = 3L)) %T>% 
      readr::write_tsv(paste0(fn_prefix, ".txt"))
    
    p <- ggplot() +
      geom_violin(df_sd, mapping = aes(x = Channel, y = SD, fill = Channel), 
                  alpha = .8, size = .25, linetype = 0, 
                  draw_quantiles = c(.95, .99)) +
      geom_boxplot(df_sd, mapping = aes(x = Channel, y = SD), 
                   width = 0.1, lwd = .2, fill = "white") +
      stat_summary(df_sd, mapping = aes(x = Channel, y = SD), 
                   fun = "mean", geom = "point",
                   shape=23, size=2, fill="white", alpha=.5) +
      labs(title = expression(""), 
           x = expression("Channel"), 
           y = expression("SD ("*log[2]*"FC)")) + 
      geom_text(data = mean_sd, aes(x = Channel, label = SD, y = SD + 0.2), 
                size = 5, colour = "red", alpha = .5)
    
    if (!is.null(ymax)) {
      if (is.null(ybreaks)) {
        ybreaks <- ifelse(ymax > 1, 0.5, ifelse(ymax > 0.5, 0.2, 0.1))
      }

      p <- p + 
        scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, ybreaks))
    }

    if (is.null(theme) || purrr::is_function(theme)) 
      theme <- theme_psm_violin
    
    p <- p + theme

    if (flip_coord) {
      p <- p + coord_flip()
      width_temp <- width
      width <- height
      height <- width_temp
      rm(list = "width_temp")
    }
    
    ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot", "width", "height", 
                                           "limitsize"))
    
    my_call <- rlang::expr(ggplot2::ggsave(filename = !!filepath, plot = !!p, 
                                           width = !!width, height = !!height, 
                                           limitsize = FALSE, 
                                           !!!ggsave_dots))
    
    ok <- tryCatch(eval(my_call, rlang::caller_env()), 
                   error = function (e) NA_character_)
    
    if (is.na(ok)) {
      my_call$width <- 40
      suppressWarnings(try(eval(my_call, rlang::caller_env())))
    }
  }
}


#' Violin plots of reporter-ion intensity per TMT_Set and LCMS_injection
#' 
#' @inheritParams info_anal
#' @inheritParams sd_violin
rptr_violin <- function(df, filepath, width, height) 
{
  df <- df %>% 
    dplyr::select(which(not_all_NA(.)))
  
  if (!ncol(df)) {
    message("No intensity columns available.")
    return (NULL)
  }
  
  df_int <- df %>% 
    `names<-`(gsub("^N_I|^I", "", names(.))) 
  
  Levels <- names(df_int)
  
  df_int <- df_int %>%
    tidyr::gather(key = "Channel", value = "Intensity") %>%
    dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
    dplyr::filter(!is.na(Intensity))
  
  mean_int <- df_int %>% 
    dplyr::group_by(Channel) %>% 
    dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
    dplyr::mutate(Intensity = round(Intensity, digits = 1))
  
  p <- ggplot() +
    geom_violin(df_int, mapping = aes(x = Channel, y = log10(Intensity), fill = Channel), 
                alpha = .8, size = .25, linetype = 0) +
    geom_boxplot(df_int, mapping = aes(x = Channel, y = log10(Intensity)), 
                 width = 0.2, lwd = .2, fill = "white") +
    stat_summary(df_int, mapping = aes(x = Channel, y = log10(Intensity)), 
                 fun = "mean", geom = "point",
                 shape = 23, size = 2, fill = "white", alpha = .5) +
    labs(title = expression("Ions"), 
         x = expression("Channel"), 
         y = expression("Intensity ("*log[10]*")")) + 
    geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), 
              size = 5, colour = "red", alpha = .5) + 
    theme_psm_violin
  
  try(ggplot2::ggsave(filepath, p, width = width, height = height, units = "in"))
}


#' Geometric mean
#' 
#' @param x A data frame.
#' @param ... Additional arguments for \code{mean}.
my_geomean <- function (x, ...) 
{
  x <- mean(log10(x), ...)
  10^x
}


#' Counts phospho.
#'
#' Not currently used.
#'
#' @param filename The file name of peptide table.
#' @param collapses The modified residues to be collapsed with unmodified
#'   counterparts.
#' @param pat_mod Pattern of modification
#' @inheritParams normPSM
count_phosphopeps <- function(filename = "Peptide.txt", collapses = c("mn"), 
                              rm_allna = FALSE, pat_mod = "sty") 
{
  dat_dir <- get_gl_dat_dir()
  
  df <- read.csv(file.path(dat_dir, "Peptide", filename), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") 
  
  not_single_sample <- (sum(grepl("^log2_R[0-9]{3}", names(df))) > 1)
  
  if (not_single_sample && rm_allna) {
    df <- df %>% 
      dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
  }
  
  id <- match_call_arg(normPSM, group_psm_by)
  use_lowercase_aa <- match_call_arg(normPSM, use_lowercase_aa)
  
  if (use_lowercase_aa) {
    pat_mod2 <- paste0("[", pat_mod, "]")
    df_phos <- df %>% dplyr::filter(grepl(pat_mod2, .[[id]]))
    n_phos_peps <- nrow(df_phos)
    
    pat <- paste(strsplit(collapses, "")[[1]], collapse = "|") %>% 
      paste0("(", ., ")")
    
    df_phos <- df_phos %>% 
      dplyr::mutate(pep_seq_mod2 = gsub(pat, paste0("@", "\\1"), !!rlang::sym(id))) %>% 
      dplyr::mutate_at(vars("pep_seq_mod2"), ~ purrr::map_chr(.x, my_upper, "@")) %>% 
      dplyr::arrange(-pep_locprob) %>% 
      dplyr::filter(!duplicated(pep_seq_mod2)) %>% 
      dplyr::mutate(phos_class = dplyr::case_when(
        pep_locprob >= .75 ~ 1L,
        pep_locprob >= .50 ~ 2L,
        pep_locprob >= .25 ~ 3L,
        pep_locprob >= 0 ~ NA_integer_,
      ))
  }
  else {
    hdr_files <- list.files(file.path(dat_dir, "PSM/cache"), "_header\\.txt$", 
                            full.names = TRUE)
    
    if (!length(hdr_files)) {
      message("Mascot header file not found")
      return(NULL)
    }
    
    hdr_files <- hdr_files[[1]]
    hdr <- readLines(hdr_files)
    
    lvs <- local({
      l1 <- grep("Variable modifications,-----", hdr)[[1]]
      l2 <- grep("Search Parameters,------", hdr)[[1]]
      lvs <- hdr[(l1+1):(l2-1)]
      lvs <- lvs[lvs != "\"\""]
      lvs <- lvs[!grepl("-term", lvs)]
    })
    
    pat_mod <- lapply(strsplit(toupper(pat_mod), "")[[1]], function (x) {
      l <- lvs[grepl(paste0(" \\(.*", x), lvs)]
      l <- gsub("\"", "", l)
      gsub("([^\\,])\\,.*", "\\1", l)
    }) %>% 
      unlist() %>% 
      unique() %>% 
      paste(collapse = "")

    pat_mod2 <- paste0("[", pat_mod, "]")
    df_phos <- df %>% dplyr::filter(grepl(pat_mod2, .[[id]]))
    n_phos_peps <- nrow(df_phos)
    
    ## convert collapses: "mn" -> "24"
    collapses <- lapply(strsplit(toupper(collapses), "")[[1]], function (x) {
      l <- lvs[grepl(paste0(" \\(.*", x), lvs)]
      l <- gsub("\"", "", l)
      gsub("([^\\,])\\,.*", "\\1", l)
    }) %>% 
      unlist() %>% 
      unique() %>% 
      paste(collapse = "")
    
    pat <- paste(strsplit(collapses, "")[[1]], collapse = "|") %>% 
      paste0("(", ., ")")

    df_phos <- df_phos %>% 
      dplyr::mutate(pep_seq_mod2 = gsub(pat, "0", !!rlang::sym(id))) %>% 
      dplyr::arrange(-pep_locprob) %>% 
      dplyr::filter(!duplicated(pep_seq_mod2)) %>% 
      dplyr::mutate(phos_class = dplyr::case_when(
        pep_locprob >= .75 ~ 1L,
        pep_locprob >= .50 ~ 2L,
        pep_locprob >= .25 ~ 3L,
        pep_locprob >= 0 ~ NA_integer_,
      ))
  }

  n_phos_sites <- sum(stringr::str_count(df_phos[[id]], pat_mod2))
  
  n_classes_12 <- df_phos %>% 
    dplyr::filter(phos_class <= 2L) %$% 
    stringr::str_count(.[[id]], pat_mod2) %>% 
    sum()
  
  ans <- data.frame(n_peps = n_phos_peps, 
                    n_sites = n_phos_sites, 
                    n_12 = n_classes_12)
  out_name <- paste0("phospep_", gsub("\\.[^.]*$", "", filename), ".tsv")
  write.table(ans, file.path(dat_dir, "Peptide/cache", out_name), 
              sep = "\t", row.names = FALSE)

  ans
}


#' Peptide mis-cleavage counts
#' 
#' NOt currently used.
#' 
count_pepmiss <- function() 
{
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "PSM/cache"), recursive = TRUE, showWarnings = FALSE)
  
  rmPSMHeaders()
  
  filelist <- list.files(path = file.path(dat_dir, "PSM/cache"),
                        pattern = "^F[0-9]{6}_hdr_rm.csv$")
  
  if (!length(filelist)) 
    stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
  df <- purrr::map(filelist, ~ {
    data <- read.delim(file.path(dat_dir, "PSM/cache", .x), sep = ',', 
                       check.names = FALSE, header = TRUE, 
                       stringsAsFactors = FALSE, quote = "\"", 
                       fill = TRUE , skip = 0)
    
    data$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    data <- data %>% 
      dplyr::filter(!duplicated(pep_seq))
    
    tot <- nrow(data)
    mis <- data %>% dplyr::filter(pep_miss > 0) %>% nrow()
    
    tibble::tibble(total = tot, miscleavage = mis, percent = miscleavage/tot)
  }) %>% 
    do.call(rbind, .)
  
  write.csv(df, file.path(dat_dir, "PSM/cache/miscleavage_nums.csv"), 
            row.names = FALSE)
  
  invisible(df)
}


#' Helper of row filtration
#'
#' \code{contain_str}: contain a literal string; "PEPTIDES" contain_str "TIDE".
#' 
#' @param match A character string containing the pattern for matching.
#' @param vars A character string of the name of a variable. The default is
#'   FALSE.
#' @param ignore.case Logical; if TRUE, ignores case when matching.
#' @examples
#' \donttest{
#' pepHist(
#'   col_select = BI,
#'   scale_log2r = TRUE,
#'   filter_peps = exprs(contain_chars_in("sty", pep_seq_mod)),
#'   scale_y = FALSE,
#'   ncol = 4,
#'   filename = "BI_pSTY_scaley_no.png",
#' )
#' }
#' @export
contain_str <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  grepl(match, vars, fixed = TRUE, ignore.case)
}

#' Helper of row filtration
#'
#' \code{contain_chars_in}: contain some of the characters in a literal string;
#' "PEPTIDES" contain_chars_in "XP".
#' 
#' @rdname contain_str
#' @export
contain_chars_in <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Helper of row filtration
#'
#' \code{not_contain_str}" not contain a literal string; "PEPTIDES"
#' not_contain_str "TED".
#' 
#' @rdname contain_str
#' @export
not_contain_str <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  !grepl(match, vars, fixed = TRUE, ignore.case)
}

#' Helper of row filtration
#'
#' \code{not_contain_chars_in}: not contain any of the characters in a literal
#' string; "PEPTIDES" not_contain_chars_in  "CAB".
#' 
#' @rdname contain_str
#' @export
not_contain_chars_in <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  !grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}

#' Helper of row filtration
#'
#' \code{start_with_str}: start with a literal string. "PEPTIDES" start_with_str
#' "PEP".
#' 
#' @rdname contain_str
#' @export
start_with_str <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("^", match), vars, fixed = FALSE, ignore.case)
}

#' Helper of row filtration
#'
#' \code{end_with_str}: end with a literal string. "PEPTIDES" end_with_str
#' "TIDES".
#' 
#' @rdname contain_str
#' @export
end_with_str <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  grepl(paste0(match, "$"), vars, fixed = FALSE, ignore.case)
}

#' Helper of row filtration
#'
#' \code{start_with_chars_in}: start with one of the characters in a literal
#' string. "PEPTIDES" start_with_chars_in "XP".
#' 
#' @rdname contain_str
#' @export
start_with_chars_in <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  grepl(paste0("^[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Helper of row filtration
#'
#' \code{ends_with_chars_in}: end with one of the characters in a literal
#' string. "PEPTIDES" ends_with_chars_in "XS".
#' 
#' @rdname contain_str
#' @export
ends_with_chars_in <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[", match, "]$"), vars, fixed = FALSE, ignore.case)
}


#' Helper of row filtration
#'
#' \code{rows_are_all}: rows are all
#' @rdname contain_str
rows_are_all <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  !grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}


#' Helper of row filtration
#'
#' \code{rows_are_all}: rows are all
#' @rdname contain_str
rows_are_not_all <- function (match, vars, ignore.case = FALSE) 
{
  stopifnot(rlang::is_string(match), nchar(match) > 0L)
  grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}


#' Concatenates formula(s) to varargs of dots
#' 
#' @param fmls A character vector of formula(s)
#' @param dots A character vector of formula(s) in \code{dots}
#' @param fml_nms A character vector containing the names of \code{fmls}.
#' @inheritParams info_anal
concat_fml_dots <- function(fmls = NULL, fml_nms = NULL, 
                            dots = NULL, anal_type = "zzz") 
{
  dat_dir <- get_gl_dat_dir()
  
  if ((!purrr::is_empty(fmls)) && (anal_type == "GSEA")) 
    return(c(dots, fmls))
  
  if (!length(fmls)) {
    fml_file <-  file.path(dat_dir, "Calls/pepSig_formulas.rda")
    
    if (file.exists(fml_file)) {
      load(file = fml_file)
      
      if (!is.null(fml_nms)) {
        stopifnot(all(fml_nms %in% names(pepSig_formulas)))
        
        pepSig_formulas <- pepSig_formulas %>% 
          .[names(.) %in% fml_nms]
      }
      
      dots <- c(dots, pepSig_formulas)
    } else {
      stop("Run both `pepSig()` and `prnSig()` first.", 
           call. = FALSE)
    }
  } 
  else {
    match_fmls(fmls)
    dots <- c(dots, fmls)
  }
  
  dots
}


#' Rolls up genes
#' 
#' @param df A data frame
#' @param cols Column indexes
gn_rollup <- function (df, cols) 
{
  if (! "gene" %in% names(df)) 
    return(df)
  
  dfa <- df %>% 
    dplyr::select(gene, cols) %>% 
    dplyr::filter(!is.na(gene)) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(~ median(.x, na.rm = TRUE))

  dfb <- df %>% 
    dplyr::select(-cols) %>% 
    dplyr::select(-which(names(.) %in% c("prot_cover"))) %>% 
    dplyr::filter(!is.na(gene)) %>% 
    dplyr::filter(!duplicated(gene))
  
  if ("prot_cover" %in% names(df)) {
    dfc <- df %>% 
      dplyr::select(gene, prot_cover) %>% 
      dplyr::filter(!is.na(gene), !is.na(prot_cover)) %>% 
      dplyr::group_by(gene) 
    
    dfc <- suppressWarnings(
      dfc %>% dplyr::summarise_all(~ max(.x, na.rm = TRUE))
    )
  } 
  else {
    dfc <- df %>% 
      dplyr::select(gene) %>% 
      dplyr::mutate(prot_cover = NA)
  }
  
  df <- list(dfc, dfb, dfa) %>% 
    purrr::reduce(right_join, by = "gene") %>% 
    dplyr::filter(!is.na(gene), !duplicated(gene))
}


#' Compares dot-dot-dot between prior and current
#' 
#' @param call_nm The name of a function call.
#' @param curr_dots The values of the current dots.
#' @param pattern The pattern for comparison.
identical_dots <- function(call_nm, curr_dots, pattern) 
{
  dat_dir <- get_gl_dat_dir()
  
  file <- file.path(dat_dir, "Calls", paste0(call_nm, ".rda"))
  
  if (!file.exists(file)) 
    return(FALSE)
  
  load(file = file)
  identical(call_pars %>% .[grepl(pattern, names(.))], curr_dots)
}


#' Complete cases among sample IDs in label_scheme_sub, not label_scheme
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
my_complete_cases <- function (df, scale_log2r, label_scheme_sub) 
{
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  NorZ_ratios <- find_NorZ(scale_log2r)
  
  # in case that reference(s) are included in label_scheme_sub, 
  # the values are zero for single ref or non-NA/NaN for multi ref 
  # and will not affect the complete.cases assessment
  rows <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% label_scheme_sub$Sample_ID)) %>% 
    complete.cases(.)
  
  if (sum(rows) == 0) 
    stop("None of the cases are complete.", call. = FALSE)
  
  df[rows, ]
}


#' The union of named list 
#' 
#' names will be kept after the unification
#' 
#' @param x A list of values.
#' @param y Another list of values.
my_union <- function (x, y) 
{
  x %>% 
    .[! names(.) %in% names(y)] %>% 
    c(y)
}


#' Finds the base names of files 
#' 
#' (not currently used)
#' 
#' @param filenames A character vector of filenames
find_fn_bases <- function (filenames) 
{
  gsub("\\.[^.]*$", "", filenames) # %>% .[1]
}


#' Finds the extensions of files
#' 
#' @param filename A character string of filename
#' @param type the type of filename extension
find_fn_exts <- function (filename, type = "text") 
{
  purrr::map_chr(filename, ~ {
    if (!grepl("\\.", .x)) {
      fn_ext <- switch(
        text = "txt",
        graphics = "png",
        stop("Invalid extension in file name(s).", Call. = FALSE)
      )
    } 
    else {
      fn_ext <- gsub("^.*\\.([^.]*)$", "\\1", .x) # %>% .[1]
    }
  })
}


#' Checks duplicate argments in 'dots'
#' 
#' @param blacklist A character vector of variable names.
#' @param ... A list of arguments for checking.
check_dots <- function (blacklist = NULL, ...) 
{
  dots <- rlang::enexprs(...)
  dups <- purrr::map_lgl(names(dots), ~ .x %in% blacklist)
  nms <- names(dots[dups])
  
  if (!rlang::is_empty(nms)) 
    stop("Do not use argument(s): ", purrr::reduce(nms, paste, sep = ", "), 
         call. = FALSE)
}


#' Checks depreciated arguments
#' 
#' @param blacklist A list of paired character vectors. In each pair, the first
#'   element is the old argument name if incurred and the second is the new name
#'   for reminding.
#' @param ... A list of arguments for checking. The depreciated argument(s) are
#'   no longer in the formalArgs of a function. If present, they will be in
#'   \code{...}.
check_depreciated_args <- function (blacklist = NULL, ...) 
{
  dots <- rlang::enexprs(...)
  old_args <- purrr::map_chr(blacklist, `[[`, 1)
  new_args <- purrr::map_chr(blacklist, `[[`, 2)
  
  depreciated <- purrr::map_lgl(names(dots), ~ .x %in% old_args)
  nms <- names(dots[depreciated])
  
  ind <- which(old_args %in% nms)
  old_args <- old_args %>% .[ind]
  new_args <- new_args %>% .[ind]
  
  if (!rlang::is_empty(nms)) {
    message("Depreciated argument(s): ", 
            purrr::reduce(old_args, paste, sep = ", "))
    stop("Use replacement argument(s): ", 
         purrr::reduce(new_args, paste, sep = ", "), 
         call. = FALSE)
  }
}


#' Forces 'complete_cases = TRUE' at 'impute_na = FALSE'
#' 
#' @inheritParams prnHM
to_complete_cases <- function (complete_cases = FALSE, impute_na = FALSE) 
{
  warn_msg1 <- "Coerce `complete_cases = TRUE` at `impute_na = FALSE`."
  
  if (!(impute_na || complete_cases)) {
    complete_cases <- TRUE
    rlang::warn(warn_msg1)
  }
  
  invisible(complete_cases)
}


#' Checks \code{gset_nms}
#'
#' It is OK that \code{gset_nms} not in "go_sets", "c2_msig", "kegg_sets",
#' "kinsub" as custom gene sets may be applied.
#'
#' @inheritParams prnGSPA
check_gset_nms <- function (gset_nms) 
{
  customs <- gset_nms %>% 
    .[! . %in% c("go_sets", "c2_msig", "kegg_sets", "kinsub")]
  
  if (!purrr::is_empty(customs)) {
    message("Apply custom database(s): ", 
            purrr::reduce(customs, paste, sep = ", ", .init = ""))
  }

  soft_dep <- gset_nms %>% .[. == "kegg_sets"]
  
  if (!rlang::is_empty(soft_dep)) {
    warning("Data set(s) softely depreciated: ", 
            purrr::reduce(soft_dep, paste, sep = ", ", .init = ""), 
            call. = FALSE)
  }

  invisible(gset_nms)
}


#' Checks the suitability of existing param files for MGKernel
#' 
#' @param filepath A file path to \code{"MGKernel_params_N.txt"}
ok_existing_params <- function (filepath) 
{
  dat_dir <- get_gl_dat_dir()
  # load(file = file.path(dat_dir, "label_scheme.rda"))
  label_scheme <- load_ls_group(dat_dir, label_scheme)
  
  if (!file.exists(filepath)) 
    return(TRUE)
  
  params <- readr::read_tsv(file.path(filepath), 
                            col_types = cols(Component = col_double()))

  missing_samples <- local({
    par_samples <- unique(params$Sample_ID) %>% .[!is.na(.)]
    expt_samples <- unique(label_scheme$Sample_ID)
    setdiff(expt_samples, par_samples)
  })
  
  if (!length(missing_samples)) 
    return(TRUE)
  
  try(unlink(file.path(dat_dir, "Peptide/Histogram/MGKernel_params_[NZ].txt")))
  try(unlink(file.path(dat_dir, "Protein/Histogram/MGKernel_params_[NZ].txt")))

  warning(
    "\nEntries(s) in `expt_smry.xlsx` missing from `MGKernel_params` files: \n  ", 
    purrr::reduce(missing_samples, paste, sep = ", "), 
    "\nThis is probably due to the modification of values under the `Sample_ID` columns.", 
    "\n=========================================================", 
    "\n=== The previous `MGKernel_params` files are deleted. ===", 
    "\n=========================================================\n", 
    call. = FALSE)
  
  invisible(return)  
} 


#' Finds the environment of a function
#' 
#' @param name The name of a function.
#' @param env The environment.
#' @export
env_where <- function(name, env = rlang::caller_env()) 
{
  if (identical(env, rlang::empty_env())) 
    stop("Can't find ", name, call. = FALSE)
  else if (rlang::env_has(env, name)) 
    env
  else 
    env_where(name, rlang::env_parent(env))
}


#' Generates cut points.
#' 
#' @param x A sequence of numeric values with probable NA value.
#' @inheritParams prnHist
set_cutpoints <- function (cut_points = NULL, x = NULL) 
{
  if (is.null(cut_points) && is.null(x)) 
    stop("`cut_points` and `x` can not be both NULL.", 
         call. = FALSE)
    
  if (!is.null(x)) {
    oks <- x %>% .[!is.na(.)] %>% is.numeric() %>% all()
    
    if (!oks) 
      stop("All numeric `x` are required for setting cut points.", 
           call. = FALSE)
  }
  
  if (!is.null(cut_points)) {
    oks <- cut_points %>% is.numeric() %>% all()
    
    if (!oks) 
      stop("All numeric values are required for `cut_points`.", 
           call. = FALSE)
  }

  if (is.null(x) && !is.null(cut_points)) {
    cut_points <- cut_points %>% 
      c(-Inf, ., Inf) %>% 
      .[!duplicated(.)] %>% 
      .[order(.)]
    
    return(cut_points)
  }
  
  if (is.null(cut_points) && !is.null(x)) {
    cut_points <- x %>% 
      quantile() %>% 
      c(-Inf, ., Inf) %>% 
      round(digits = 1)
    
    return(cut_points)
  }
  
  cut_points %>% 
    c(-Inf, ., Inf) %>% 
    .[!duplicated(.)] %>% 
    .[order(.)]
}


#' A wrapper of \code{set_cutpoint}.
#' 
#' @param df A data frame.
#' @inheritParams prnHist
set_cutpoints2 <- function(cut_points, df) 
{
  if (is.null(cut_points) || is.infinite(cut_points)) 
    return(c(mean_lint = -Inf, mean_lint = Inf))
  
  if (any(is.na(cut_points))) {
    cut_points <- cut_points[1]
    
    if (is.null(names(cut_points))) 
      cut_nm <- "mean_lint"
    else 
      cut_nm <- names(cut_points)
    
    if (! cut_nm %in% names(df)) {
      warning("Column `", cut_nm, "` not found in `df`;", 
              "instead use column `mean_lint`.", 
              call. = FALSE)
      cut_nm <- "mean_lint"
    }
    
    if (!is.numeric(df[[cut_nm]])) 
      stop("Column `", cut_nm, "` is not numeric.", call. = FALSE)

    cut_points <- quantile(df[[cut_nm]], na.rm = TRUE) 
    seqs <- seq_along(cut_points)
    names(cut_points)[seqs] <- cut_nm
  } 
  else {
    cut_nm <- gsub("1$", "", names(cut_points[1]))
    
    if (!length(cut_nm)) 
      cut_nm <- "mean_lint"
    
    if (!cut_nm %in% names(df)) {
      warning("Column `", cut_nm, "` not found in `df`;", 
              "instead use column `mean_lint`.", 
              call. = FALSE)
      cut_nm <- "mean_lint"
    }
    
    if (!is.numeric(df[[cut_nm]])) 
      stop("Column `", cut_nm, "` is not numeric.", call. = FALSE)
    
    cut_points <- set_cutpoints(cut_points, df[[cut_nm]])
    seqs <- seq_along(cut_points)
    names(cut_points)[seqs] <- cut_nm
  }

  invisible(cut_points)
}


#' Loads cRAPs
#' 
#' Load the cRAP proteins by the types of accessions
#' @param acc_types The types of protein accessions
load_craps <- function(acc_types) 
{
  data(package = "proteoQ", prn_annot_crap, envir = environment())
  
  if (identical(acc_types, "other_acc")) {
    craps <- prn_annot_crap %>% 
      dplyr::select("fasta_name") %>% 
      dplyr::filter(!is.na(fasta_name), !duplicated(fasta_name)) %>% 
      unlist()
  } 
  else {
    craps <- prn_annot_crap %>% 
      dplyr::select(which(names(.) %in% acc_types)) %>% 
      tidyr::gather("acc_type", "prot_acc") %>% 
      dplyr::filter(!is.na(prot_acc), !duplicated(prot_acc)) %>% 
      dplyr::select(prot_acc) %>% 
      unlist()
  }
}


#' Finds the columns of reporter-ion intensity.
#' 
#' @inheritParams TMT_levels
find_int_cols <- function (TMT_plex = 10L) 
{
  oks <- c(0L, 6L, 10L, 11L, 16L, 18L)
  
  if (! TMT_plex %in% oks) 
    warning("`TMT_plex` is not one of ", paste(oks, collapse = ", "), ".", 
            call. = FALSE)
  
  col_int <- if (TMT_plex == 18L) {
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131N", "I131C", 
      "I132N", "I132C", "I133N", "I133C", "I134N", 
      "I134C", "I135N")
  } 
  else if (TMT_plex == 16L) {
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131N", "I131C", 
      "I132N", "I132C", "I133N", "I133C", "I134N")
  } 
  else if (TMT_plex == 11L) {
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131N", "I131C")
  } 
  else if (TMT_plex == 10L) {
    c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
      "I130N", "I130C", "I131")
  } 
  else if(TMT_plex == 6L) {
    c("I126", "I127", "I128", "I129", "I130", "I131")
  } 
  else {
    NULL
  }
}


#' Finds the columns of reporter-ion ratios.
#' 
#' @inheritParams TMT_levels
find_ratio_cols <- function (TMT_plex = 10L) 
{
  col_ratio <- if (TMT_plex == 18L) {
    c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
      "R130N", "R130C", "R131N", "R131C", 
      "R132N", "R132C", "R133N", "R133C", "R134N", 
      "R134C", "R135N")
  } 
  else if (TMT_plex == 16L) {
    c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
      "R130N", "R130C", "R131N", "R131C", 
      "R132N", "R132C", "R133N", "R133C", "R134N")
  } 
  else if (TMT_plex == 11L) {
    c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
      "R130N", "R130C", "R131N", "R131C")
  } 
  else if (TMT_plex == 10L) {
    c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
      "R130N", "R130C", "R131")
  } 
  else if(TMT_plex == 6L) {
    c("R127", "R128", "R129", "R130", "R131")
  } 
  else {
    NULL
  }
}


#' Processes MaxQuant mqpar.xml
#' 
#' @inheritParams load_expts
#' @param mqpar The name of .xml file. The default is "mqpar.xml".
#' @param filename The name of metadata file.
#' @param out The name of output .xml file.
make_mq_meta <- function (dat_dir = proteoQ:::get_gl_dat_dir(), mqpar = "mqpar.xml", 
                          filename = "mq_meta.txt", out = "new_mqpar.xml") 
{
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("\n====================================================================", 
         "\nNeed package \"xml2\" for this function to work.",
         "\n====================================================================",
         call. = FALSE)
  }
  
  lookup <- readr::read_tsv(file.path(dat_dir, filename)) %>% 
    dplyr::rename(RAW_File = Name, Sample_ID = Experiment) %>% 
    dplyr::filter(rowSums(!is.na(.)) > 0)

  n_smpls <- nrow(lookup)
  
  parent <- xml2::read_xml(file.path(dat_dir, mqpar))
  children <- xml2::xml_children(parent)
  contents <- parent %>% xml2::xml_contents() 
  
  filePaths <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "filePaths")]])
  len <- length(filePaths)
  if (len > n_smpls) {
    xml2::xml_text(filePaths[(n_smpls + 1):len]) <- ""
    filePaths <- filePaths[1:n_smpls]
  }
  pre <- filePaths %>% xml2::xml_text() %>% gsub("(.*\\\\|/).*", "\\1", .)
  xml2::xml_text(filePaths) <- paste0(pre, lookup$RAW_File)
  
  experiments <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "experiments")]])
  len <- length(experiments)
  if (len > n_smpls) {
    xml2::xml_text(experiments[(n_smpls + 1):len]) <- ""
    experiments <- experiments[1:n_smpls]
  }
  xml2::xml_text(experiments) <- lookup$Sample_ID
  
  fractions <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "fractions")]])
  len <- length(fractions)
  if (len > n_smpls) {
    xml2::xml_text(fractions[(n_smpls + 1):len]) <- ""
    fractions <- fractions[1:n_smpls]
  }
  
  ptms <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "ptms")]])
  len <- length(ptms)
  if (len > n_smpls) {
    xml2::xml_text(ptms[(n_smpls + 1):len]) <- ""
    ptms <- ptms[1:n_smpls]
  }
  
  paramGroupIndices <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "paramGroupIndices")]])
  referenceChannel <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "referenceChannel")]])
  
  xml2::write_xml(parent, file.path(dat_dir, out), options = "format")
}


#' Processes MaxQuant mqpar.xml
#' 
#' @inheritParams load_expts
#' @param mqpar The name of .xml file. The default is "mqpar.xml".
#' @param filename The name of metadata file.
#' @param out The name of output .xml file.
#' @examples
#' \dontrun{
#' dat_dir <- "Y:/qzhang/MaxQuant"
#' pre_mq_meta()
#' make_mq_meta2()
#' }
make_mq_meta2 <- function (dat_dir = proteoQ:::get_gl_dat_dir(), 
                           mqpar = "mqpar.xml", 
                           filename = "mq_meta.txt", out = "new_mqpar.xml") 
{
  fn_suffix <- gsub("^.*\\.([^\\.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^\\.]*$", "", filename)
  
  lookup <- if (fn_suffix %in% c("xls", "xlsx")) {
    readxl::read_excel(file.path(dat_dir, filename))
  } else if (fn_suffix == "csv") {
    readr::read_csv(file.path(dat_dir, filename))
  } else if (fn_suffix == "txt") {
    readr::read_tsv(file.path(dat_dir, filename))
  } else {
    stop(filename, " needs to be '.xls' or '.xlsx'.")
  }
  
  lookup <- lookup %>% 
    dplyr::filter(rowSums(!is.na(.)) > 0)

  if ((!"PTM" %in% names(lookup)) || all(is.na(lookup$PTM))) {
    lookup$PTM <- "False"
  }
  
  if ((!"Fraction" %in% names(lookup)) || all(is.na(lookup$Fraction))) {
    lookup$Fraction <- 1
  }
  
  if (!"Group" %in% names(lookup) || all(is.na(lookup$Group))) {
    lookup$Group <- 0
  }
  
  if ((!"Reference" %in% names(lookup)) || all(is.na(lookup$Reference))) {
    lookup$Reference <- ""
  }
  
  if (!grepl("\\.d$", lookup$Name[[1]])) {
    lookup$Name <- paste0(lookup$Name, ".d")
  }

  mqxml <- readLines(file.path(dat_dir, mqpar)) %>% 
    update_nodes("filePaths", "string", lookup$Name) %>% 
    update_nodes("experiments", "string", lookup$Experiment) %>% 
    update_nodes("fractions", "short", lookup$Fraction) %>% 
    update_nodes("ptms", "boolean", lookup$PTM) %>% 
    update_nodes("paramGroupIndices", "int", lookup$Group) %>% 
    update_nodes("referenceChannel", "string", lookup$Reference) %T>% 
    writeLines(file.path(dat_dir, out))
}


#' Helper of \link{make_mq_meta2}
#' 
#' @param mqxml The XML from mqpar.xml.
#' @param field The field to be updated.
#' @param type The type of data.
#' @param vals The Values of data.
update_nodes <- function(mqxml, field = "filePaths", type = "string", vals) 
{
  row_1 <- grep(paste0("<", field, ">"), mqxml)
  row_2 <- grep(paste0("</", field, ">"), mqxml)
  
  lines <- paste0(paste0("\t\t\t<", type, ">"), 
                  vals, 
                  paste0("</", type, ">"))
  
  if (field == "filePaths") {
    paths <- mqxml[(row_1+1):(row_2-1)] %>% 
      gsub("^.*<string>", "", .) %>% 
      gsub("</string>$", "", .) %>% 
      gsub("(.*\\\\).*", "\\1", .) %>% 
      gsub("(.*/).*", "\\1", .)
  } else {
    paths <- ""
  }
  
  lines <- paste0(paste0("\t\t\t<", type, ">"), 
                  paths, 
                  vals, 
                  paste0("</", type, ">"))
  
  bf <- mqxml[1:row_1]
  af <- mqxml[row_2:length(mqxml)]
  mqxml <- bf %>% append(lines) %>% append(af) 
}


#' Updates MaxQuant experimentalDesignTemplate.txt
#' 
#' @param dat_dir The working directory.
#' @param expt_smry The filename of proteoQ metadata.
#' @param filename The filename of MaxQuant metadata.
pre_mq_meta <- function (dat_dir = proteoQ:::get_gl_dat_dir(), 
                         expt_smry = "expt_smry.xlsx", 
                         filename = "combined/experimentalDesignTemplate.txt")
{
  meta_proteoq <- file.path(dat_dir, expt_smry) %>% 
    readxl::read_excel() %>% 
    dplyr::select(c("Sample_ID", "RAW_File")) %>% 
    dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File))
  
  meta_mq <- file.path(dat_dir, filename) %>% 
    readr::read_tsv() %>% 
    dplyr::left_join(meta_proteoq, by = c("Name" = "RAW_File")) %>% 
    dplyr::mutate(Experiment = Sample_ID) %>% 
    dplyr::select(-Sample_ID) %T>%
    readr::write_tsv(file.path(dat_dir, "mq_meta.txt"))
}


#' Adds numeric indices for peptides or proteins
#' 
#' @param df A data frame.
#' @param col_in The column key to be indexed.
#' @param col_out The name of the output column.
add_entry_ids <- function (df, col_in = "pep_seq", col_out = "pep_index") 
{
  if (col_out %in% names(df)) {
    return(reloc_col_after_last(df, col_out))
  }
  
  ids <- unique(df[[col_in]])
  
  seqs <- seq_along(ids) %>% `names<-`(ids)
  seqs <- tibble::tibble(!!col_in := names(seqs), !!col_out := seqs) 
  
  df %>% dplyr::left_join(seqs, by = col_in)
}


#' Checks conflicts in function arguments.
#' 
#' Conflicts between wrapper and function to be wrapped. Note that \code{f}
#' contains function name but not package name.
#' 
#' if `method` in the formalArgs and "m" in the "..." of `foo`, 
#' foo(m = 2) without `method` interpreted as `method`.
#' Call foo(method = kmeans, m =2) to avoid such ambiguity. 
#' 
#' @param w Expression; the wrapper function.
#' @param f Expression; the function to be wrapped.
#' @param excludes Character string; arguments from \code{f} to be excluded for
#'   checking.
check_formalArgs <- function(w = anal_prnTrend, f = cmeans, excludes = NULL) 
{
  w <- rlang::as_string(rlang::enexpr(w))
  f <- rlang::as_string(rlang::enexpr(f))
  
  args_w <- formalArgs(w) %>% .[! . == "..."]
  args_f <- formalArgs(f) %>% .[! . == "..."]

  if (!is.null(excludes)) {
    args_f <- args_f %>% .[! . %in% excludes]
  }
  
  if (length(args_f)) {
    purrr::walk(args_w, ~ {
      ambi <- stringr::str_detect(.x, paste0("^", args_f)) %>% 
        which()
      
      if (length(ambi)) {
        warning("Ambiguous arguments between \"", 
                purrr::reduce(args_f[ambi], paste, sep = ", "), 
                "\" from `", f, 
                "` and \"", .x, 
                "\" from `", w, 
                "`.\n", 
                "### Include explicitly \"", .x, "\"", 
                " when calling `", w, "`. ###\n",  
                call. = FALSE)
      }
    })
  }
}


#' Checks conflicts in function arguments.
#' 
#' Conflicts between wrapper and function to be wrapped. Note that \code{f}
#' contains both function name and package name.
#' 
#' @inheritParams check_formalArgs
check_formalArgs2 <- function(w = anal_prnTrend, f = NMF::nmf, excludes = NULL) 
{
  # for error messages
  nm_f <- rlang::enexpr(f)
  if (any(grepl("::", nm_f))) {
    nm_f <- nm_f %>% as.character() 
    stopifnot(length(nm_f) ==3)
    nm_f <- nm_f %>% `[`(3)
  }
  
  w <- rlang::as_string(rlang::enexpr(w))

  args_w <- formalArgs(w) %>% .[! . == "..."]
  args_f <- formalArgs(f) %>% .[! . == "..."]
    
  if (!is.null(excludes)) {
    args_f <- args_f %>% .[! . %in% excludes]
  }
  
  purrr::walk(args_w, ~ {
    ambi <- stringr::str_detect(.x, paste0("^", args_f)) %>% 
      which()
    
    if (!purrr::is_empty(ambi)) {
      warning("Ambiguous arguments between \"", 
              purrr::reduce(args_f[ambi], paste, sep = ", "), 
              "\" from `", nm_f, 
              "` and \"", .x, 
              "\" from `", w, 
              "`.\n", 
              "### Include explicitly \"", .x, "\"", 
              " when calling `", w, "`. ###\n",  
              call. = FALSE)
    }
  })
}


#' Sets dot-dot-dot for ggsave.
#' 
#' @param dots Arguments passed on.
#' @param dummies Arguments automated for ggsave.
set_ggsave_dots <- function(dots, dummies = c("filename", "plot", "device", "path")) 
{
  if (purrr::is_empty(dots)) return(dots)
  
  args <- formalArgs(ggplot2::ggsave) %>% .[! . == "..."]
  
  dots %>% 
    .[names(.) %in% args] %>% 
    .[! names(.) %in% dummies]
}


#' Finds the type of search engine
#' 
#' @inheritParams load_expts
find_search_engine <- function(dat_dir = NULL) 
{
  if (is.null(dat_dir)) 
    dat_dir <- get_gl_dat_dir()
  
  pat_mascot <- "^F[0-9]+.*\\.csv$"
  pat_mq <- "^msms.*\\.txt$"
  pat_mf <- "^psm.*\\.tsv$"
  pat_sm <- "^PSMexport.*\\.ssv$"
  pat_pq <- "^psm[QC]{1}.*\\.txt$"
  # pat_msgf <- "(^peptide)|(^protein)|(^psm)\\.tsv$"
  pat_msgf <- "^psmMSGF.*\\.txt$$"
  
  mascot <-list.files(path = file.path(dat_dir), pattern = pat_mascot) %>% 
    length() %>% `>`(0L)
  mq <- list.files(path = file.path(dat_dir), pattern = pat_mq) %>% 
    length() %>% `>`(0L)
  mf <- list.files(path = file.path(dat_dir), pattern = pat_mf) %>% 
    length() %>% `>`(0L)
  sm <- list.files(path = file.path(dat_dir), pattern = pat_sm) %>% 
    length() %>% `>`(0L)
  pq <- list.files(path = file.path(dat_dir), pattern = pat_pq) %>% 
    length() %>% `>`(0L)
  
  msgf <- list.files(path = file.path(dat_dir), pattern = pat_msgf) %>% 
    length() %>% `>`(0L)
  
  if (FALSE) {
    msgf <- local({
      files <- list.files(path = file.path(dat_dir), pattern = "\\.tsv$")
      files <- files[!grepl("(^peptide)|(^protein)|(^psm)", files)]
      length(files) > 0L
    })
  }

  engines <- c(mascot = mascot, mq = mq, mf = mf, sm = sm, pq = pq, msgf = msgf)
  oks <- names(engines[engines])

  if (!length(oks)) {
    stop("No PSM files found.\n", 
         "File name(s) need to be in the format of \n",
         paste(c(paste("  Mascot", pat_mascot, sep = ": "), 
                 paste("  MaxQuant", pat_mq, sep = ": "), 
                 paste("  MSFragger", pat_mf, sep = ": "), 
                 paste("  proteoQ", pat_pq, sep = ": "), 
                 paste("  MSGF", pat_msgf, sep = ": "), 
                 paste("  Spectrum Mill", pat_sm, sep = ": ")), 
               collapse = "\n"))
  }
  
  if (length(oks) > 1L)
    stop("Multiple data types found: ", paste(oks, collapse = ", "))
  
  invisible(oks)
}


#' Checks missing values in a ggplot2 object.
#' 
#' @param p A ggplot2 object.
check_ggplot_aes <- function(p) 
{
  q <- p %>% ggplot_build() %>% 
    `[[`(1) %>% 
    `[[`(1) %>% 
    dplyr::select(which(not_all_NA(.))) 
  
  cols <- q %>% 
    purrr::map_lgl(anyNA)
  
  if (any(cols)) {
    warning("\nMismatches in the lengths of vectors ", 
            "between selected samples and the aesthetics of `", 
            purrr::reduce(names(which(cols)), paste, sep = ", "), "`.\n", 
            "Fix the missing values in the corresponding columns in `expt_smry.xlsx`.\n", 
            "Otherwise, may run into errors.\n", 
            call. = FALSE)
  }
}


#' Standardizes range values.
#' 
#' @param from_to A numeric vector of length 2.
prep_range <- function(from_to = c(0, 1)) 
{
  stopifnot(length(from_to) == 2)
  
  stopifnot(from_to[2] > from_to[1], 
            from_to[1] >= 0,
            from_to[2] <= 100) 
  
  if (from_to[2] <= 1) from_to <- from_to * 100
  
  from_to
}


#' Finds the deliminator of a file.
#' 
#' @param filename A file name.
find_delim <- function (filename) 
{
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  
  if (fn_suffix %in% c("txt", "tsv")) {
    delim <- "\t"
  } else if (fn_suffix == "csv") {
    delim <- ","
  } else if (fn_suffix == "ssv") {
    delim <- ";"
  } else {
    stop("File extension needs to be one of 'txt', 'tsv', 'csv' or 'ssv'.", 
         call. = FALSE)
  }
  
  invisible(delim)
}



#' Flattens lists recursively
#' 
#' @param x Lists
recur_flatten <- function (x) 
{
  if (!inherits(x, "list")) {
    list(x)
  } else {
    unlist(c(purrr::map(x, recur_flatten)), recursive = FALSE)
  }
}


#' Subsets elements with attributes kept.
#'
#' @param data The data vector.
#' @param nm The names of elements in the vector for removals.
subset_keepattr <- function (data, nm = c("M", "N")) 
{
  attrs <- attributes(data)
  
  for (i in seq_along(nm)) {
    data <- data %>% .[! names(.) == nm[i]]
  }
  
  attrs$names <- names(data)
  attributes(data) <- attrs
  
  data
}


#' Lists to a data frame
#' 
#' A simple replacement of plyr::ldply. 
#' 
#' @param x Lists of vectors.
list_to_dataframe <- function (x) 
{
  nms <- names(x)
  lens <- unlist(lapply(x, length))
  mlen <- max(lens)
  
  for (i in seq_along(x)) {
    d <- mlen - lens[i]
    if (d) x[[i]] <- c(x[[i]], rep(NA, d))
  }
  
  ans <- do.call(rbind, x)
  
  if (is.null(nms)) {
    # rownames(ans) <- seq_len(nrow(ans))
    cns <- seq_len(mlen)
  } else {
    # rownames(ans) <- nms
    ans <- cbind(nms, ans)
    cns <- c(".id", seq_len(mlen))
  }

  ans <- data.frame(ans)
  colnames(ans) <- cns
  rownames(ans) <- seq_len(nrow(ans))

  invisible(ans)
}


#' Finds the list of psmQ files.
#' 
#' @inheritParams load_expts
find_psmQ_files <- function (dat_dir) 
{
  filelistQ <- list.files(path = file.path(dat_dir), pattern = "^psmQ.*\\.txt$")
  filelistC <- list.files(path = file.path(dat_dir), pattern = "^psmC.*\\.txt$")
  
  if (length(filelistQ) && length(filelistC)) {
    message("Both `psmQ[...].txt` and `psmC[...].txt` are present; ", 
            "\"psmQ[...].txt\" will be used in baseline analysis.\n", 
            "Need \"psmC[...].txt\" for optional matches between runs in LFQ.")
    
    filelist <- filelistQ
  } 
  else if (length(filelistQ)) {
    filelist <- filelistQ
  } 
  else if (length(filelistC)) {
    filelist <- filelistC
    
    stop(paste("No PSM files `psmQ[...].txt` under", file.path(dat_dir), ".",
               "\nMake sure that the names of PSM files start with `psmQ`."), 
         call. = FALSE)
  } 
  else {
    stop(paste("No PSM files `psmQ[...].txt` under", file.path(dat_dir), ".",
               "\nMake sure that the names of PSM files start with `psmQ`."), 
         call. = FALSE)
  }
  
  filelist
}


#' Gets column types.
#' 
#' @import readr
get_col_types <- function () 
{
  # just a remainder of mzion columns
  col_types_pq <- cols(
    prot_acc = col_character(), 
    prot_issig = col_logical(), 
    prot_isess = col_logical(),
    prot_tier = col_integer(), 
    prot_hit_num = col_integer(), 
    prot_family_member = col_integer(), 
    prot_es = col_number(), 
    prot_es_co = col_number(), 
    pep_seq = col_character(), 
    pep_n_ms2 = col_integer(), 
    pep_scan_title = col_character(), 
    pep_exp_mz = col_number(),
    pep_exp_mr = col_number(), 
    pep_exp_z = col_character(), 
    pep_calc_mr = col_number(), 
    pep_delta = col_number(),
    pep_tot_int = col_number(), 
    pep_ret_range = col_number(), 
    pep_scan_num = col_character(), # timsTOF
    pep_mod_group = col_integer(), 
    pep_frame = col_integer(), 
    pep_fmod = col_character(),
    pep_vmod = col_character(),
    pep_isdecoy = col_logical(),
    pep_ivmod = col_character(),
    pep_len = col_integer(), 
    
    pep_ms2_moverzs = col_character(),
    pep_ms2_ints = col_character(), 
    pep_ms2_theos = col_character(),
    pep_ms2_theos2 = col_character(),
    pep_ms2_exptints = col_character(),
    pep_ms2_exptints2 = col_character(),
    
    pep_n_matches = col_integer(), 
    pep_n_matches2 = col_integer(),
    pep_ms2_deltas = col_character(),
    pep_ms2_ideltas = col_character(),
    pep_ms2_deltas2 = col_character(), 
    pep_ms2_ideltas2 = col_character(),
    pep_ms2_deltas_mean = col_double(),
    pep_ms2_deltas_sd = col_double(),
    
    pep_issig = col_logical(),
    pep_score = col_double(),
    pep_expect = col_double(),
    pep_rank = col_integer(), 
    pep_locprob = col_double(),
    pep_locdiff = col_double(),
    pep_rank_nl = col_integer(), 
    pep_literal_unique = col_logical(),
    pep_razor_unique = col_logical(),
    raw_file = col_character(), )
  
  col_types <- cols(
    prot_hit_num = col_integer(),
    prot_family_member = col_integer(),
    prot_acc = col_character(),
    prot_desc = col_character(),
    prot_score = col_double(),
    prot_mass = col_double(),
    prot_matches = col_integer(),
    prot_matches_sig = col_integer(),
    prot_sequences = col_integer(),
    prot_sequences_sig = col_integer(),
    prot_n_psm = col_integer(),
    prot_n_pep = col_integer(),
    prot_len = col_integer(),
    prot_pi = col_double(),
    prot_tax_str = col_character(),
    prot_tax_id = col_integer(),
    prot_seq = col_character(),
    prot_empai = col_double(),
    prot_icover = col_double(),
    prot_cover = col_double(),
    
    prot_index = col_integer(),
    prot_issig = col_logical(),
    prot_isess = col_logical(),
    prot_tier = col_integer(),
    prot_es = col_double(),
    prot_es_co = col_double(),
    pep_query = col_integer(),
    pep_rank = col_integer(),
    pep_n_psm = col_integer(),
    pep_isbold = col_logical(),
    pep_isunique = col_logical(),
    pep_literal_unique = col_logical(),
    pep_razor_unique = col_logical(),
    pep_tot_int = col_double(),
    pep_unique_int = col_double(),
    pep_razor_int = col_double(),
    pep_exp_mz = col_double(),
    pep_exp_mr = col_double(),
    pep_exp_z = col_character(),
    pep_calc_mr = col_double(),
    pep_delta = col_double(),
    pep_score = col_double(),
    pep_homol = col_double(),
    pep_ident = col_double(),
    pep_expect = col_double(),
    
    pep_res_before = col_character(),
    pep_seq = col_character(),
    pep_seq_mod = col_character(),
    pep_res_after = col_character(),
    pep_start = col_integer(),
    pep_end = col_integer(),
    pep_len = col_integer(),
    pep_miss = col_integer(),
    pep_frame = col_integer(),
    pep_var_mod = col_character(),
    pep_var_mod_pos = col_character(),
    pep_summed_mod_pos = col_character(),
    pep_local_mod_pos = col_character(),
    pep_num_match = col_integer(), # Mascot
    pep_scan_title = col_character(),
    pep_index = col_integer(),
    pep_scan_range = col_character(), # timsTOF
    
    pep_n_exp_z = col_integer(), 
    pep_ret_sd = col_double(),
    pep_n_nl = col_integer(),
    
    pep_ret_range = col_double(),
    pep_ms2_sumint = col_double(),
    pep_n_ions = col_integer(),
    pep_locprob = col_double(),
    pep_locdiff = col_double(),
    pep_phospho_locprob = col_double(),
    pep_phospho_locdiff = col_double(),
    pep_ions_first = col_character(),
    pep_ions_second = col_character(),
    pep_ions_third = col_character(),
    pep_n_ms2 = col_integer(),
    pep_scan_num = col_character(),
    pep_mod_group = col_integer(),
    pep_fmod = col_character(),
    pep_vmod = col_character(),
    pep_isdecoy = col_logical(),
    pep_ivmod = col_character(),
    pep_issig = col_logical(),
    pep_rank_nl = col_integer(),
    
    gene = col_character(),
    fasta_name = col_character(),
    uniprot_acc = col_character(),
    uniprot_id = col_character(),
    refseq_acc = col_character(),
    other_acc = col_character(),
    entrez = col_integer(),
    species = col_character(),
    acc_type = col_character(),
    shared_prot_accs = col_character(),
    shared_genes = col_character(),
    kin_attr = col_logical(),
    kin_class = col_character(),
    kin_order = col_integer(),
    pep_istryptic = col_logical(),
    pep_semitryptic = col_character(),
    
    dat_file = col_character(),
    raw_file = col_character(),
    mean_lint = col_double(),
    count_nna = col_integer(),
    TMT_Set = col_integer(),
  )
  
  nms <- names(col_types$cols)
  
  stopifnot(length(nms) == length(unique(nms)))
  
  col_types
}


#' Wrapper of writing Excel
#' 
#' @param df A data frame.
#' @param sheetName An Excel sheet name.
#' @param filename An output file name.
#' @inheritParams load_expts
write_excel_wb <- function (df, sheetName = "Setup", dat_dir = NULL, filename) 
{
  if (is.null(dat_dir)) 
    dat_dir <- get_gl_dat_dir()
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = sheetName)
  openxlsx::writeData(wb, sheet = sheetName, df)
  openxlsx::saveWorkbook(wb, file.path(dat_dir, filename), overwrite = TRUE)
}


#' Loads label_scheme or label_scheme_group.
#'
#' Prefer label_scheme_group over label_scheme.
#'
#' @param filename A file name of "label_scheme" or "label_scheme_full".
#' @param prefer_group Logical; if TRUE, prefer label_scheme_group over
#'   label_scheme etc.
#' @inheritParams load_expts
load_ls_group <- function (dat_dir, filename = "label_scheme", prefer_group = TRUE) 
{
  filename <- rlang::as_string(rlang::enexpr(filename))
  
  if (grepl("\\.rda$", filename)) {
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
  }
  else {
    fn_prefix <- filename
    filename <- paste0(filename, ".rda")
  }
  
  file <- file.path(dat_dir, filename)
  file2 <- file.path(dat_dir, paste0(fn_prefix, "_group.rda"))
  fn_prefix2 <- paste0(fn_prefix, "_group")
  
  if (file.exists(file2) && prefer_group) {
    fi <- file2
    nm <- fn_prefix2
  }
  else if (file.exists(file)) {
    fi <- file
    nm <- fn_prefix
  }
  else {
    stop("`label_schem` not found under ", dat_dir)
  }
  
  load(fi)
  
  get(nm, envir = environment(), inherits = FALSE)
}


#' Parses \code{col_select}.
#' 
#' @inheritParams standPep
#' @param label_scheme Experiment summary
parse_col_select <- function (col_select, label_scheme) 
{
  sids <- label_scheme[["Sample_ID"]]
  
  if (is.null(sids))
    stop("Column `Sample_ID` not found in metadata.", call. = FALSE)
  
  if (col_select == "Sample_ID" && all(is.na(sids)))
    stop("All values under  column `Sample_ID` are NA.", call. = FALSE)
  
  if (col_select == "Sample_ID")
    return(col_select)
  
  if (is.null(label_scheme[[col_select]])) {
    col_select <- "Sample_ID"
    
    warning("Column `", col_select, "` not existed. ", 
            "Used column `Sample_ID` instead.", call. = FALSE)
  } 
  else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    col_select <- "Sample_ID"
    
    warning("No samples under column '", col_select, "`. ", 
            "Used column `Sample_ID` instead.", call. = FALSE)
  }
  
  col_select
}


#' Parses file name.
#' 
#' @param filename A file name.
#' @param dat_dir A working directory.
#' @param must_exists Logical; if TRUE, the file must be present.
parse_filename <- function (filename, dat_dir, must_exists = FALSE) 
{
  file <- file.path(dat_dir, filename)
  
  if (must_exists && !file.exists(file))
    stop(filename, " not found under '", dat_dir, "'.", call. = FALSE)
  
  if (!grepl("\\.", filename))
    stop("File name has no extension: ", filename)
  
  fn_suffix <- gsub("^.*\\.([^\\.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^\\.]*$", "", filename)
  
  if (nchar(fn_prefix) <= 0L)
    stop("File name has no base: ", filename)
  
  list(fn_prefix = fn_prefix, fn_suffix = fn_suffix)
}


#' Subsets Bruker's MGF data
#' 
#' Applied after the fixes of TITLE lines.
#'
#' @param file A file name
#' @param begin_offset The number of lines before a BEGIN line.
#' @param charge_offset The number lines after a BEGIN line to a following
#'   CHARGE line.
#' @param topn_ms2ions Top-n MS2 ions to be retained.
subsetBrukerMGF <- function (file, begin_offset = 5L, charge_offset = 5L, 
                             topn_ms2ions = Inf) 
{
  message("Processing: ", file)
  
  lines  <- readLines(file)
  begins <- .Internal(which(stringi::stri_startswith_fixed(lines, "BEGIN IONS")))
  ends   <- .Internal(which(stringi::stri_endswith_fixed(lines, "END IONS")))
  hdrs   <- 1:(begins[1]-begin_offset-1L)
  
  z_lns <- lines[begins+charge_offset]
  oks   <- grepl("^CHARGE", z_lns) & (z_lns != "CHARGE=1+")
  b_oks <- begins[oks] - begin_offset
  e_oks <- ends[oks]
  
  ranges <- mapply(function (x, y) x:y, b_oks, e_oks, SIMPLIFY = TRUE)
  ranges <- do.call(`c`, ranges)
  ranges <- c(hdrs, ranges)
  lines <- lines[ranges]
  rm(list = c("begins", "ends", "b_oks", "e_oks", "oks", "ranges", "z_lns"))

  if (is.infinite(topn_ms2ions)) {
    writeLines(lines, file)
    return(NULL)
  }
  
  pat_mgf <- mzion:::find_mgf_type(file)
  
  type_mgf <- pat_mgf$type
  n_bf_begin <- pat_mgf$n_bf_begin
  n_spacer <- pat_mgf$n_spacer
  n_hdr <- pat_mgf$n_hdr
  n_to_pepmass <- pat_mgf$n_to_pepmass
  n_to_title <- pat_mgf$n_to_title
  n_to_scan <- pat_mgf$n_to_scan
  n_to_rt <- pat_mgf$n_to_rt
  n_to_charge <- pat_mgf$n_to_charge
  sep_ms2s <- pat_mgf$sep_ms2s
  nfields_ms2s <- pat_mgf$nfields_ms2s
  sep_pepmass <- pat_mgf$sep_pepmass
  nfields_pepmass <- pat_mgf$nfields_pepmass
  raw_file <- pat_mgf$raw_file
  
  begins <- .Internal(which(stringi::stri_startswith_fixed(lines, "BEGIN IONS")))
  ends <- .Internal(which(stringi::stri_endswith_fixed(lines, "END IONS")))
  
  ## MS2
  # (-1L: one line above "END IONS")
  ms2s <- mapply(function (x, y) lines[(x + n_hdr) : (y - 1L)], 
                 begins, ends, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  ms2s <- lapply(ms2s, stringi::stri_split_fixed, pattern = sep_ms2s, 
                 n = nfields_ms2s, simplify = TRUE)
  
  ms2_moverzs <- lapply(ms2s, function (x) as.numeric(x[, 1]))
  # not as.integer; intensity may be > .Machine$integer.max
  ms2_ints <- lapply(ms2s, function (x) as.numeric(x[, 2]))
  rm(list = c("ms2s"))
  
  mz_n_int <- mzion:::sub_mgftopn(ms2_moverzs = ms2_moverzs, 
                                  ms2_ints = ms2_ints, 
                                  topn_ms2ions = topn_ms2ions, 
                                  min_ms2mass = 0L, 
                                  max_ms2mass = 5000L)
  
  ms2_moverzs <- mz_n_int[["ms2_moverzs"]]
  ms2_ints <- mz_n_int[["ms2_ints"]]
  
  ms2_out <- mapply(function (x, y) paste0(x, "\t", y, "\t"), 
                    ms2_moverzs, ms2_ints, SIMPLIFY = FALSE)
  
  aboves <- begins - (n_bf_begin + 1L)
  belows <- begins + (n_hdr - 1L)
  ms1_out <- mapply(function (x, y) lines[x:y], aboves, belows, SIMPLIFY = FALSE)
  
  ans <- mapply(function (m1, m2) c(m1, m2, "END IONS"), ms1_out, ms2_out)
  ans <- unlist(ans)
  writeLines(c(lines[hdrs], ans), file)
}


#' Parallel subsetBrukerMGF
#' 
#' @param filepath A file path to MGF.
#' @param n_cores The number of CPU cores.
#' @inheritParams subsetBrukerMGF
#' @export
msubsetBrukerMGF <- function (filepath, begin_offset = 5L, charge_offset = 5L, 
                              topn_ms2ions = Inf, n_cores = 1L) 
{
  files <- list.files(filepath, pattern = "\\.mgf$", full.names = TRUE, 
                      recursive = TRUE)
  
  len <- length(files)
  
  if (!len)
    stop("No MGF files found.")
  
  if (n_cores > 1L)
    n_cores <- min(parallel::detectCores(), n_cores, len)
  
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterApply(cl, files, subsetBrukerMGF, 
                           begin_offset = begin_offset, 
                           charge_offset = charge_offset, 
                           topn_ms2ions = topn_ms2ions)
    parallel::stopCluster(cl)
  }
  else
    lapply(files, subsetBrukerMGF, begin_offset = begin_offset, 
           charge_offset = charge_offset, topn_ms2ions = topn_ms2ions)
}


#' Reprocesses of Bruker MGF files.
#' 
#' To be deleted
#' 
#' @param file A file name with prepending path.
#' @examples
#' \donttest{
#' filepath = "F:/timsTOF/mgf"
#' mprocBrukerMGF(filepath)
#' }
procBrukerMGF_v1 <- function (file) 
{
  message("Processing: ", file)
  
  lines <- readLines(file)
  
  fi <- local({
    hdr <- lines[1:50L]
    ln <- hdr[grepl("###.*\\.d$", hdr)]
    nm <- gsub("###\t(.*\\.d$)", "\\1", ln)
    paste0("File:", nm)
  })
  
  rows <- which(stringi::stri_startswith_fixed(lines, "TITLE="))
  
  tits <- lapply(rows, function (i) 
    gsub("^([^,]*?), (.*)", paste("\\1", fi, "\\2", sep = ", "), lines[i]))
  
  lines[rows] <- tits
  
  writeLines(unlist(lines), file)
}


#' Reprocesses of Bruker MGF files.
#' 
#' @param file A file name with prepending path.
#' 
#' @examples
#' \donttest{
#' filepath = "F:/timsTOF/mgf"
#' mprocBrukerMGF(filepath)
#' }
procBrukerMGF <- function (file) 
{
  message("Processing: ", file)
  
  lines <- readLines(file)
  
  fi <- local({
    hdr <- lines[1:50L]
    ln <- hdr[grepl("###.*\\.d$", hdr)]
    nm <- gsub("###\t(.*\\.d$)", "\\1", ln)
    paste0("~File:", nm, "~")
  })
  
  rows <- which(stringi::stri_startswith_fixed(lines, "TITLE="))
  
  tits <- lapply(rows, function (i) 
    gsub("^([^,]*?),[ ]{0,1}(.*)", paste("\\1", fi, "\\2", sep = ", "), lines[i]))
  
  lines[rows] <- tits
  
  writeLines(unlist(lines), file)
}


#' Batch-reprocessing of Bruker MGF files.
#' 
#' @param filepath A file path to MGF.
#' @param n_cores The number of CPU cores.
#' @export
mprocBrukerMGF <- function (filepath, n_cores = 1L) 
{
  files <- list.files(filepath, pattern = "\\.mgf$", full.names = TRUE, 
                      recursive = TRUE)
  
  len <- length(files)
  
  if (!len)
    stop("No MGF files found.")
  
  if (n_cores > 1L)
    n_cores <- min(parallel::detectCores(), n_cores, len)
  
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterApply(cl, files, procBrukerMGF)
    parallel::stopCluster(cl)
  }
  else
    lapply(files, procBrukerMGF)
}


#' Fix the \code{File} field in \code{Title} lines.
#' 
#' The first tilde in the \code{Title} lines needs to before \code{File}.
#' 
#' @param file A file name with prepending path.
fixBrukerMGF <- function (file)
{
  lines <- readLines(file)
  rows <- which(stringi::stri_startswith_fixed(lines, "TITLE="))
  tits <- lapply(rows, function (i) gsub("File:~", "~File:", lines[i]))
  lines[rows] <- tits
  
  writeLines(unlist(lines), file)
}


#' Batch-fixing of Bruker MGF files.
#' 
#' The first tilde in the \code{Title} lines needs to before \code{File}.
#' 
#' @param filepath A file path to MGF.
#' @param n_cores The number of CPU cores.
#' @export
mfixBrukerMGF <- function (filepath, n_cores = 1L) 
{
  files <- list.files(filepath, pattern = "\\.mgf$", full.names = TRUE, 
                      recursive = TRUE)
  
  len <- length(files)
  
  if (!len)
    stop("No MGF files found.")
  
  if (n_cores > 1L)
    n_cores <- min(parallel::detectCores(), n_cores, len)
  
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterApply(cl, files, fixBrukerMGF)
    parallel::stopCluster(cl)
  }
  else
    lapply(files, fixBrukerMGF)
}


#' Checks lengths of aesthetics.
#'
#' The length between metadata and manual aesthetics.
#'
#' @param label_scheme A metadata.
#' @param x A column name in metadata.
#' @param vals Values of manual aesthetics.
#' @param aes The name of an aesthetics.
check_aes_length <- function (label_scheme = NULL, x = "Size", 
                              aes = "size_manual", vals = c(3, 5)) 
{
  len_1 <- length(unique(label_scheme[[x]]))
  len_2 <- length(vals)
  
  if (len_1 != len_2)
    stop("Unequal lengths: metadata ", x, " = ", len_1, ", ", aes, " = ", len_2)
}


#' Adds column \code{raw_file}.
#' 
#' @param df A PSM table.
add_col_rawfile <- function (df)
{
  if ("raw_file" %in% names(df)) {
    df <- df %>% 
      dplyr::rename(RAW_File = raw_file) %>% 
      dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File)) %>% 
      dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File))
    
    df <- df %>%
      dplyr::mutate(pep_scan_title = gsub("\\\\", "/", pep_scan_title))
  } 
  else {
    df <- df %>%
      dplyr::mutate(pep_scan_title = gsub("\\\\", "/", pep_scan_title)) %>% 
      dplyr::mutate(RAW_File = pep_scan_title) %>% 
      dplyr::mutate(RAW_File = gsub("^.*/([^/]*)\\.raw[\\\"]{0,1}; .*", "\\1", 
                                    RAW_File)) %>% 
      dplyr::mutate(RAW_File = gsub("^.*/([^/]*)\\.d[\\\"]{0,1}; .*", "\\1", 
                                    RAW_File))
  }
  
  # (with some refseq_acc)
  df <- df %>%
    dplyr::mutate(prot_acc = gsub("\\.[0-9]*$", "", prot_acc)) 
}


#' Loads psmC.txt
#' 
#' @param file The file name of \code{psmC[...].txt}.
load_psmC <- function(file = NULL, ...) 
{
  dat_dir <- get_gl_dat_dir()
  base_name <- gsub("\\.txt$", "", file)
  
  dfC <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, file), 
                    col_types = cols(
                      prot_acc = col_character(), 
                      prot_issig = col_logical(), 
                      # prot_isess = col_logical(),
                      # prot_tier = col_integer(), 
                      # prot_hit_num = col_integer(), 
                      # prot_family_member = col_integer(), 
                      prot_es = col_number(), 
                      prot_es_co = col_number(), 
                      pep_seq = col_character(), 
                      pep_n_ms2 = col_integer(), 
                      pep_scan_title = col_character(), 
                      pep_exp_mz = col_number(),
                      pep_exp_mr = col_number(), 
                      pep_exp_z = col_character(), 
                      pep_calc_mr = col_number(), 
                      pep_delta = col_number(),
                      pep_tot_int = col_number(), 
                      pep_ret_range = col_number(), 
                      pep_scan_num = col_character(), # timsTOF
                      pep_mod_group = col_integer(), 
                      pep_frame = col_integer(), 
                      pep_fmod = col_character(),
                      pep_vmod = col_character(),
                      pep_isdecoy = col_logical(),
                      pep_ivmod = col_character(),
                      pep_len = col_integer(), 
                      pep_issig = col_logical(),
                      pep_score = col_double(),
                      pep_expect = col_double(),
                      pep_rank = col_integer(), 
                      pep_locprob = col_double(),
                      pep_locdiff = col_double(),
                      # pep_rank_nl = col_integer(), 
                      # pep_literal_unique = col_logical(),
                      # pep_razor_unique = col_logical(),
                      raw_file = col_character(), ), 
                    show_col_types = FALSE)
  )
  
  ans <- dfC %>% 
    add_col_rawfile() %>% 
    dplyr::select(c("pep_seq", "pep_tot_int", "pep_ret_range", "pep_scan_num", 
                    # "pep_scan_title", 
                    "pep_n_ms2", "pep_exp_mz", "pep_exp_mr", "pep_delta", 
                    # "pep_calc_mr", "pep_mod_group", "pep_isdecoy", "pep_len", 
                    "pep_issig", "pep_score", "pep_locprob", "pep_locdiff", 
                    "pep_rank", "pep_expect", 
                    "RAW_File", "pep_ivmod", "pep_fmod", "pep_vmod", 
                    "pep_exp_z", )) %>% 
    dplyr::mutate(I000 = pep_tot_int) %>% 
    dplyr::mutate(R000 = I000/I000, 
                  R000 = ifelse(is.infinite(R000), NA_real_, R000)) %>% 
    dplyr::mutate(dat_file = base_name)
}


