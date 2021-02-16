#'Prepares data for analysis
#'
#'\code{prepDM} prepares a minimal adequate data frame for subsequent analysis.
#'
#'@inheritParams  prnEucDist
#'@inheritParams  info_anal
#'@param sub_grp Numeric.  A list of sample IDs that will be used in subsequent
#'  analysis.
#'@param type The type of data, for example ratio or intensity.
#'@return A data frame tailored for subsequent analysis.
#'
#' @examples
#' \donttest{tempData <- prepDM(df, entrez, scale_log2r, label_scheme_sub$Sample_ID)}
#'@import dplyr
#'@importFrom magrittr %>% %T>% %$% %<>% 
prepDM <- function(df, id, scale_log2r, sub_grp, type = "ratio", anal_type) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  stopifnot(nrow(df) > 0)  
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  if (anal_type %in% c("ESGAGE", "GSVA")) {
    id <- "entrez"
  }
  
  if ((anal_type %in% c("GSEA")) && (id != "gene")) {
    stop("Primary ID is not `gene`.", call. = FALSE)
  }
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  pattern <- 
    "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("

  df <- local({
    df <- df %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(!duplicated(!!rlang::sym(id)),
                    !is.na(!!rlang::sym(id)),) 
    
    if (nrow(df) == 0) {
      stop("All values are NA under the column `", id, "` in the input data.", 
           call. = FALSE)
    }
    
    df <- df %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(rowSums(!is.na(.[, grep(NorZ_ratios, names(.)), 
                                     drop = FALSE])) > 0) %>%
      reorderCols(endColIndex = grep(pattern, names(.)), col_to_rn = id)
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
    
    if (!purrr::is_empty(trivials)) {
      warning("Samples with all NA entries being dropped: ", 
              purrr::reduce(trivials, paste, sep = ", "), 
              call. = FALSE)
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
  
  return(list(log2R = dfR, Intensity = dfI, IR = rbind(tempR, tempI)))
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
`names_pos<-` <- function(x, pos, value) {
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
reorder_files <- function(filelist) {
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
reorderCols <- function(df, endColIndex, col_to_rn) {
	if (length(endColIndex) == 0) {
	  endColIndex <- grep("I[0-9]{3}|log2_R[0-9]{3}", names(df))
	}
	
	df <- cbind(df[, -endColIndex], df[, endColIndex])

	if(sum(duplicated(df[[col_to_rn]])) > 0) {
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
reorderCols2 <- function(df = NULL, pattern = NULL) {
  if (is.null(df)) stop("`df` cannot be `NULL`.", call. = FALSE)
  
  if (is.null(pattern)) {
    pattern <- 
      "I[0-9]{3}\\(|log2_R[0-9]{3}\\(|pVal\\s+\\(|adjP\\s+\\(|log2Ratio\\s+\\(|\\.FC\\s+\\("
  }

  endColIndex <- grep(pattern, names(df))
  
  if (length(endColIndex) > 0) 
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
#' @export
ins_cols_after <- function(df = NULL, idx_bf = ncol(df), idx_ins = NULL) {
  if (is.null(df)) stop("`df` cannot be `NULL`.", call. = FALSE)
  if (is.null(idx_ins)) return(df)
  if (idx_bf >= ncol(df)) return(df)
  
  col_nms <- names(df)[idx_ins]
  
  dplyr::bind_cols(
    df[, 1:idx_bf, drop = FALSE], 
    df[, idx_ins, drop = FALSE], 
    df[, (idx_bf + 1):ncol(df), drop = FALSE] %>% dplyr::select(-col_nms),
  )
}


#' Sort by the indexes of TMT_Set and LCMS_Injection
#' 
#' "8.3"    "8.5"   "11.1"   "11.3"   "99.2"  "99.10" 
#' Not "8.3"    "8.5"   "11.1"   "11.3"   "99.10"  "99.2" 
#' Nor "11.1"   "11.3"  "8.3"    "8.5"   "99.10"  "99.2" 
#' @param nms A list of names of "TMT_Set" and "LCMS_Injection" separated by dot.
sort_tmt_lcms <- function (nms) {
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


#' Replace zero intensity with NA
#'
#' \code{na_zeroIntensity} replaces zero intensity with NA to avoid -Inf in
#' log10 transformation.
#'
#' @param df A list of file names.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>% %T>% %$% %<>% 
na_zeroIntensity <- function(df) {
	ind <- grep("I[0-9]{3}|^LFQ\\sintensity\\s|^Intensity\\s", names(df))

	temp <- df[, ind]
	temp[temp == 0] <- NA
	df[, ind] <- temp

	return(df)
}

#' Find all zero rows
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
aggrNums <- function(f) {
	function (df, id, ...) {
		id <- rlang::as_string(rlang::enexpr(id))
		dots <- rlang::enexprs(...)

		df %>%
			dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise_all(~ f(.x, !!!dots))
	}
}


#' Sum top_n
#'
#' \code{aggrTopn} summarizes \code{log2FC} and \code{intensity} by the
#' descriptive statistics of \code{c("mean", "median", "sum")}. Note the
#' difference to \link{TMT_top_n}, which uses mean statistics.
#'
#' @param f A function for data summarization.
#' @examples \donttest{df_num <- aggrTopn(sum)(df, prot_acc, 3, na.rm = TRUE)}
aggrTopn <- function(f) {
  function (df, id, n, ...) {
    id <- rlang::as_string(rlang::enexpr(id))
    dots <- rlang::enexprs(...)
    
    # OK for LFQ: rowSums against single column
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
tmt_wtmean <- function (x, id, na.rm = TRUE, ...) {
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
TMT_top_n <- function (x, id, ...) {
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
not_all_zero <- function (df) {
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
not_all_NA <- function (df) {
  stopifnot(is.data.frame(df))
  colSums(!is.na(df), na.rm = TRUE) > 0
}


#' Finds non-all-NaN column(s)
#'
#' @param x A vector.
#' @param ... Additional arguments for \code{sum}.
#' @examples \donttest{purrr::map_lgl(mtcars, not_all_nan, na.rm = TRUE)}
not_all_nan <- function(x, ...) {
  stopifnot(is.vector(x))
  sum(is.nan(x), ...) != length(x)
}


#' Finds all-NaN column(s)
#' 
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param ... The same in \code{sum}.
#' @examples \donttest{purrr::map_lgl(mtcars, is_all_nan, na.rm = TRUE)}
is_all_nan <- function(x, ...) {
  stopifnot(is.vector(x))
  sum(is.nan(x), ...) == length(x)
}


#' Finds a trivial column
#' 
#' @param df A data frame of \code{log2FC} and \code{intensity}.
#' @param col The index of column.
#' @examples \donttest{is_trivial_col(mtcars, 1)}
is_trivial_col <- function (df, col) {
  stopifnot(length(col) == 1)
  
  x <- df[, col]
  
  x[x == 0] <- NA
  x[is.nan(x)] <- NA
  
  all(is.na(x))
}


#' Find a trivial vector
#' 
#' @param x A numeric vector
is_trivial_dbl <- function (x) {
  stopifnot(is.vector(x))
  
  if (!is.numeric(x)) return(FALSE)
  
  x[x == 0] <- NA
  x[is.nan(x)] <- NA
  
  all(is.na(x))
}


#' Replace a trivial vector with NA values.
#' 
#' @param x A numeric vector
replace_trivial_with_na <- function (x) {
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
colAnnot <- function (annot_cols = NULL, sample_ids = NULL, annot_colnames = NULL) {
	if (is.null(annot_cols)) return(NA)

  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme.rda"))
	exists <- annot_cols %in% names(label_scheme)

	if (sum(!exists) > 0) {
		warning(paste0("Column '", annot_cols[!exists], "'",
		               " not found in \'label_scheme\' and will be skipped."), 
		        call. = FALSE)
		annot_cols <- annot_cols[exists]
	}

	if (length(annot_cols) == 0) return(NA)

	x <- label_scheme %>%
		dplyr::filter(Sample_ID %in% sample_ids) %>%
		dplyr::select(annot_cols, Sample_ID) 
	
	idx <- !not_all_NA(x)
	if (sum(idx) > 0) stop("No aesthetics defined under column(s) `", 
	                       purrr::reduce(names(x)[idx], paste, sep = ", "), "`.", 
	                       call. = FALSE)

	x <- x %>%
		dplyr::filter(!grepl("^Empty\\.", .[["Sample_ID"]]),
		              !is.na(.[["Sample_ID"]])) %>%
		data.frame(check.names = FALSE) %>%
		`rownames<-`(.[["Sample_ID"]])

	if (any(duplicated(x[["Sample_ID"]]))) {
	  stop("Duplicated sample IDs found.\n", 
	       call. = FALSE)
	}

	if (!"Sample_ID" %in% annot_cols) x <- x %>% dplyr::select(-Sample_ID)

	if (ncol(x) == 0) return(NA)

	if ("TMT_Set" %in% names(x)) {
		x <- x %>%
			tibble::rownames_to_column() %>%
			dplyr::mutate(TMT_Set = as.factor(TMT_Set)) %>%
		  `rownames<-`(NULL) %>% 
			tibble::column_to_rownames(var = "rowname")
	}

	return(x)
}


#' Customizes the colors in column annotation
#' 
#' @param annotation_col The same as in \link[pheatmap]{pheatmap}.
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
setHMColor <- function (annotation_col) {
	ncol <- ncol(annotation_col)

	if(is.null(ncol)) return (NA)
	if(ncol == 0) return (NA)

	suppressWarnings(
	  annotation_col <- annotation_col %>%
	    dplyr::mutate_at(vars(one_of("TMT_Set")), ~ factor(TMT_Set)) %>%
	    dplyr::mutate_if(is.character, as.factor)  
	)

	palette <- lapply(names(annotation_col), function(x) {
		n <- nlevels(annotation_col[[x]])

		palette <- if(n <= 9 & n >= 3) {
		  brewer.pal(n, name = "Set2")
		} else if(n > 9) {
			colorRampPalette(brewer.pal(n = 7, "Set1"))(n)
		} else if(n == 2) {
		  c("#66C2A5", "#FC8D62")
		} else if(n == 1) {
			# c("#E41A1C")
		  c("#66C2A5")
		} else if(n == 0) {
			colorRampPalette(brewer.pal(n = 9, "YlOrBr"))(100)
		}

		names(palette) <- levels(annotation_col[[x]])

		return(palette)
	})

	names(palette) <- names(annotation_col)

	return(palette)
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
setHMlims <- function (x, xmin, xmax) {
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
ratio_toCtrl <- function(df, id, label_scheme_sub, nm_ctrl) {
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
imputeNA <- function (id, overwrite = FALSE, ...) {
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

		return(x)
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

		df[, grep("N_log2_R", names(df))] <- handleNA(df[, grep("N_log2_R", names(df))], ...)

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


#'Impute NA values for peptide data
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
pepImp <- function (...) {
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
prnImp <- function (...) {
  check_dots(c("id"), ...)
  id <- match_call_arg(normPSM, group_pep_by)
  imputeNA(id = !!id, ...)
}


#' Species lookup
#' 
#' @inheritParams load_dbs
sp_lookup <- function(species) {
  switch(species, 
         human = "hs",
         mouse = "mm",
         rat = "rn",
         unknown = "unknown",
         species
  )    
}


#' Taxonomy lookup
#' @inheritParams load_dbs
taxid_lookup <- function(species) {
  switch (species,
          human = 9606, 
          mouse = 10090,
          rat = 10116, 
          unknown = 999999, 
          999999
  )
}


#' Reversed taxonomy lookup
#' @inheritParams load_dbs
taxid_lookup_rev <- function(species) {
  switch (species,
          "9606" = "human", 
          "10090" = "mouse",
          "10116" = "rat", 
          "unknown"
  )
}


#' Species lookup UpperLower (Ul)
#' 
#' @inheritParams load_dbs
sp_lookup_Ul <- function(species) {
  switch(species, 
         human = "Hs",
         mouse = "Mm",
         rat = "Rn",
         unknown = "Unknown",
         "unknown"
  )    
}


#' Add gene IDs for RefSeq
#' 
#' @inheritParams add_entrez
add_refseq_gene <- function (acc_lookup, acc_type) {
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
                            Call. = FALSE)
  )
  
  if (all(species %in% c("human", "mouse", "rat"))) {
    filelist <- paste0(abbr_acc, "entrez_", abbr_sp)
    data(package = "proteoQ", list = filelist)
    
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


#' Add Entrez IDs
#' 
#' @param acc_lookup A data frame of protein accession lookups
#' @param acc_type The type of protein accessions
#' @inheritParams parse_fasta
add_entrez <- function (acc_lookup, acc_type, warns = TRUE) {
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
  
  if ((!all(species %in% def_sps))) {
    if (warns) {
      warning("No default `entrez` lookups for species other than: ", 
              purrr::reduce(def_sps, paste, sep = ", "), ".\n",
              "To annotate, provide the file path and name(s) via argument `entrez`.\n", 
              "See also `?Uni2Entrez` and `?Ref2Entrez` for custom `entrez` lookups.", 
              call. = FALSE)
      
      warning("At the `entrez = NULL` default, ignore non-default species: \n", 
              purrr::reduce(species %>% .[! . %in% def_sps], paste, sep = ", "), 
              call. = FALSE)
    }

    species <- species %>% .[. %in% def_sps]
  }
  
  if (purrr::is_empty(species)) {
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
                       stop("Unknown `accession type`.", Call. = FALSE)
  )
  
  # multiple uniprot_acc(s) can share the same entrez id
  filelist <- paste0(abbr_acc, "entrez_", abbr_sp)
  
  # four columns: 
  # 1. uniprot_acc/refseq_acc/other_acc  2. gene  3. entrez  4. species
  data(package = "proteoQ", list = filelist)
  
  db <- purrr::map(filelist, ~ get(.x)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(!duplicated(.[[entrez_key]])) %>% 
    dplyr::mutate(entrez = as.numeric(entrez),
                  !!entrez_key := as.character(!!rlang::sym(entrez_key))) %>% 
    dplyr::select(entrez_key, "entrez")
  
  acc_lookup <- dplyr::left_join(acc_lookup, db, by = entrez_key)
}


#' Add custom Entrez IDs and optional overwriting species
#' 
#' @inheritParams  add_entrez
#' @inheritParams normPSM
#' @inheritParams annotKin
add_custom_entrez <- function(acc_lookup, entrez, acc_type) {
  if (!all(file.exists(entrez))) {
    stop("Wrong `entrez` file path(s) or name(s).", 
         "\nMake sure that the file type is `.rds`, not `.rda`", 
         call. = FALSE)
  }
  
  entrez_key <- switch(acc_type, 
                       uniprot_acc = "uniprot_acc", 
                       uniprot_id = "uniprot_acc", 
                       refseq_acc = "refseq_acc", 
                       other_acc = "other_acc", 
                       stop("Unknown `accession type`.", 
                            Call. = FALSE)
  )
  
  cols_ess <- c(entrez_key, "entrez")
  cols_out <- c(entrez_key, "entrez", "species")
  
  db <- purrr::map(entrez, ~ {
    db <- readRDS(.x)
    
    if (!all(cols_ess %in% names(db))) {
      stop("Need columns '", 
           purrr::reduce(cols_ess, paste, sep = ", "), 
           "' in ", .x, ".", 
           call. = FALSE)
    }
    
    if (length(db) == 0) {
      stop("Empty file: ", .x, call. = FALSE)
    }
    
    db %>% 
      dplyr::select(which(names(.) %in% cols_out)) %>% 
      { if ("species" %in% names(.)) . else .[["species"]] <- NA_character_ } %>% 
      dplyr::mutate(entrez = as.numeric(entrez)) %>% 
      dplyr::select(cols_out)
  }) %>% 
    dplyr::bind_rows() 
  
  # higher priority for custom "species" in `db`
  if (all(is.na(db$species))) {
    db[["species"]] <- NULL
  } else {
    acc_lookup[["species"]] <- NULL
  }
  
  acc_lookup %>% dplyr::left_join(db, by = entrez_key)
} 


#' Determine the protein accession type from PSM tables
#'
#' Find the protein accession from a non-cRAP entry and parse it.
#' 
#'@param df An input data frame.
parse_acc <- function(df) {
  stopifnot("prot_acc" %in% names(df))
  
  prot_accs <- df %>%
    dplyr::select(prot_acc) %>%
    unlist %>% 
    unique()
  
  pat_uni_acc <- 
    "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  pat_uni_id <- 
    "^[[:alnum:]]+_[A-Z]{1,10}$"
  pat_ref_acc <- 
    "^[XNY]{1}[MRP]{1}_"
  
  acc_types <- purrr::map_chr(prot_accs, ~ {
    if (grepl(pat_uni_id, .x)) {
      acc_type <- "uniprot_id"
    } else if (grepl(pat_ref_acc, .x)) {
      acc_type <- "refseq_acc"
    } else if (grepl(pat_uni_acc, .x)) {
      acc_type <- "uniprot_acc"
    } else {
      acc_type <- "other_acc"
    }
  }, prot_accs)
  
  df %>% 
    dplyr::left_join(data.frame(prot_acc = prot_accs, 
                                acc_type = acc_types), 
                     by = "prot_acc")
}


#' Parse FASTA for accession lookups
#'
#' NA genes are replaced with prot_acc
#' 
#' @param df An input data frame.
#' @param warns Logical; if TRUE, show warning message(s).
#' @inheritParams add_entrez
#' @inheritParams normPSM
#' @seealso \code{\link{read_fasta}} for the definition of fasta_name(s).
#' @return A lookup table, \code{acc_lookup}.
parse_fasta <- function (df, fasta, entrez, warns = TRUE) {
  
  na_genes_by_acc <- function(acc_lookup) {
    if (nrow(acc_lookup) > 0 && 
        "prot_acc" %in% names(acc_lookup) && 
        "gene" %in% names(acc_lookup)) {
      na_gene <- (is.na(acc_lookup$gene)) | (stringr::str_length(acc_lookup$gene) == 0)
      acc_lookup$gene <- as.character(acc_lookup$gene)
      acc_lookup$gene[na_gene] <- as.character(acc_lookup$prot_acc[na_gene])
    }

    return(acc_lookup)
  }
  
  na_species_by_org <- function(acc_lookup) {
    if (nrow(acc_lookup) > 0 && 
        "organism" %in% names(acc_lookup) && 
        "species" %in% names(acc_lookup)) {
      na_species <- 
        (is.na(acc_lookup$species)) | (stringr::str_length(acc_lookup$species) == 0)
      acc_lookup$species <- as.character(acc_lookup$species)
      acc_lookup$species[na_species] <- as.character(acc_lookup$organism[na_species])
    }

    return(acc_lookup)
  }
  

  my_lookup <- c(
    "Homo sapiens" = "human",
    "Mus musculus" = "mouse",
    "Rattus norvegicus" = "rat"
  )

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
  #   read_fasta(pattern = ">([^ ]+?) .*", ...) 
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
  df <- df %>% 
    parse_acc() %>% 
    dplyr::filter(nchar(prot_acc) > 0)

  # load fasta_db
  fasta_db <- local({
    if (is.null(fasta)) {
      stop("FASTA file(s) are required.", call. = FALSE)
    }
    
    if (!all(file.exists(fasta))) {
      stop("Missing FASTA file(s): \n", 
           purrr::reduce(fasta %>% .[!file.exists(.)], paste, sep = "\n"), 
           call. = FALSE)
    }
    
    purrr::map(fasta, ~ read_fasta(.x)) %>% 
      do.call(`c`, .) %>% 
      `names<-`(gsub(">", "", names(.))) %>% 
      .[!duplicated(names(.))]
  })

  # accession lookups by acc_type's
  df_splits <- df %>% split(., .$acc_type, drop = TRUE)
  
  acc_lookup <- purrr::map(df_splits, ~ {
    acc_type <- .x[["acc_type"]] %>% .[1]
    
    if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
      acc_lookup <- names(fasta_db) %>% 
        gsub("^..\\|", "", .) %>% 
        stringr::str_split("\\|", simplify = TRUE) %>% 
        data.frame() %>% 
        `colnames<-`(paste("X", seq_len(ncol(.)), sep = ".")) %>% 
        `names_pos<-`(1:2, c("uniprot_acc", "uniprot_id")) %>% 
        dplyr::select(c("uniprot_acc", "uniprot_id")) %>% 
        dplyr::bind_cols(data.frame(fasta_name = names(fasta_db)), .) %>% 
        dplyr::filter(grepl(pat_uni_acc, uniprot_acc), 
                      grepl(pat_uni_id, uniprot_id)) %>%         
        dplyr::filter(.[[acc_type]] %in% unique(.x$prot_acc)) %>% 
        dplyr::filter(!duplicated(.[[acc_type]])) %>% 
        dplyr::mutate(refseq_acc = NA_character_, 
                      other_acc = NA_character_, 
                      acc_type = acc_type, 
                      prot_acc = !!rlang::sym(acc_type))
      
      local({
        n_fasta <- length(fasta_db)
        n_lookup <- nrow(acc_lookup)
        perc <- n_lookup/n_fasta 

        if (perc < .1 && warns) {
          warning("The portion of uniprot accessions being annotated with \n", 
                  purrr::reduce(fasta, paste, sep = "\n"), " is ", 
                  format(round(perc, 3), nsmall = 3), 
                  "\nCheck the choice of fasta database(s) for probable mismatches.", 
                  call. = FALSE)
        }
      })
    } else if (acc_type == "refseq_acc") {
      acc_lookup <- tibble::tibble(fasta_name = names(fasta_db), 
                                   refseq_acc = names(fasta_db)) %>% 
        dplyr::filter(grepl(pat_ref_acc, refseq_acc)) %>% 
        dplyr::filter(.[[acc_type]] %in% unique(.x$prot_acc)) %>% 
        dplyr::filter(!duplicated(.[[acc_type]])) %>% 
        dplyr::mutate(uniprot_acc = NA_character_, 
                      uniprot_id = NA_character_, 
                      other_acc = NA_character_, 
                      acc_type = acc_type, 
                      prot_acc = !!rlang::sym(acc_type)) %>% 
        dplyr::select(c("fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                        "other_acc", "acc_type", "prot_acc"))
      
      local({
        n_fasta <- length(fasta_db)
        n_lookup <- nrow(acc_lookup)
        perc <- n_lookup/n_fasta 
        
        if (perc < .1 && warns) {
          warning("The portion of refseq accessions being annotated with \n", 
                  purrr::reduce(fasta, paste, sep = "\n"), " is ", 
                  format(round(perc, 3), nsmall = 3), 
                  "\nInsepct the specification of fasta databases for probable mismatches", 
                  call. = FALSE)
        }
      })
    } else if (acc_type == "other_acc") {
      acc_lookup <- tibble::tibble(fasta_name = names(fasta_db), 
                                   other_acc = names(fasta_db)) %>% 
        dplyr::filter(other_acc %in% unique(.x$prot_acc)) %>% 
        dplyr::filter(!duplicated(other_acc)) %>% 
        dplyr::mutate(uniprot_acc = NA_character_, 
                      uniprot_id = NA_character_, 
                      refseq_acc = NA_character_, 
                      acc_type = acc_type, 
                      prot_acc = other_acc) %>% 
        dplyr::select(c("fasta_name", "uniprot_acc", "uniprot_id", "refseq_acc", 
                        "other_acc", "acc_type", "prot_acc")) 
    }

    invisible(acc_lookup)
  }) %>% 
    dplyr::bind_rows()
  
  # --- subset and annotate fasta_db with entries in acc_lookup ---
  fasta_db <- local({
    fasta_db <- fasta_db %>% 
      .[names(.) %in% unique(acc_lookup$fasta_name)] %>% 
      .[! duplicated(names(.))]
    
    if (length(fasta_db) == 0) {
      stop("No fasta entries match protein accessions; probably wrong fasta file.", 
           call. = FALSE)
    } else {
      # write_fasta(fasta_db, file.path(dat_dir, "my_project.fasta"))
    }
    
    acc_types <- unique(df$acc_type)
    
    purrr::map(acc_types, ~ {
      if (.x == "uniprot_acc") {
        fasta_names_sub <- acc_lookup %>% 
          dplyr::filter(!is.na(!!rlang::sym(.x))) %>% 
          .[["fasta_name"]] %>% 
          unique()
        
        fasta_db_sub <- fasta_db %>% .[names(.) %in% fasta_names_sub]
        
        # add attributes
        prot_acc <- gsub("^..\\|(.*)\\|.*$", "\\1", names(fasta_db_sub))
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_acc, ~ {
          attr(.x, "prot_acc") <- .y
          return(.x)
        })
        
        prot_desc <- fasta_db_sub %>%
          purrr::map(attr, "header") %>% 
          gsub("^.*\\|.*\\s+?", "", .)
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_desc, ~ {
          attr(.x, "prot_desc") <- .y
          return(.x)
        })
        
      } else if (.x == "uniprot_id") {
        fasta_names_sub <- acc_lookup %>% 
          dplyr::filter(!is.na(!!rlang::sym(.x))) %>% 
          .[["fasta_name"]] %>% 
          unique()
        
        fasta_db_sub <- fasta_db %>% .[names(.) %in% fasta_names_sub]

        prot_acc <- gsub("^.*\\|(.*)$", "\\1", names(fasta_db_sub))
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_acc, ~ {
          attr(.x, "prot_acc") <- .y
          return(.x)
        })
        
        prot_desc <- fasta_db_sub %>%
          purrr::map(attr, "header") %>% 
          gsub("^.*\\|.*\\s+?", "", .)
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_desc, ~ {
          attr(.x, "prot_desc") <- .y
          return(.x)
        })
        
      } else if (.x == "refseq_acc") {
        fasta_names_sub <- acc_lookup %>% 
          dplyr::filter(!is.na(!!rlang::sym(.x))) %>% 
          .[["fasta_name"]] %>% 
          unique()
        
        fasta_db_sub <- fasta_db %>% .[names(.) %in% fasta_names_sub]
        
        prot_acc <- gsub("\\.[0-9]*$", "", names(fasta_db_sub))
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_acc, ~ {
          attr(.x, "prot_acc") <- .y
          return(.x)
        })
        
        prot_desc <- fasta_db_sub %>%
          purrr::map(attr, "header") %>% 
          gsub("^.*_.*\\s+?", "", .)
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_desc, ~ {
          attr(.x, "prot_desc") <- .y
          return(.x)
        })
        
      } else {
        fasta_names_sub <- acc_lookup %>% 
          dplyr::filter(!is.na(other_acc)) %>% 
          .[["fasta_name"]] %>% 
          unique()
        
        fasta_db_sub <- fasta_db %>% .[names(.) %in% fasta_names_sub]

        prot_acc <- names(fasta_db_sub)
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_acc, ~ {
          attr(.x, "prot_acc") <- .y
          return(.x)
        })
        
        prot_desc <- fasta_db_sub %>% 
          purrr::map(attr, "header") 
        fasta_db_sub <- purrr::map2(fasta_db_sub, prot_desc, ~ {
          attr(.x, "prot_desc") <- .y
          return(.x)
        })
      }
      
      invisible(fasta_db_sub)
    }) %>% 
      do.call(`c`, .)
  }) 
  fasta_db <- fasta_db %>% .[!duplicated(names(.))]
  save(fasta_db, file = file.path(dat_dir, "fasta_db.rda"))

  # -- add columns prot_desc, prot_mass, prot_len ---
  acc_lookup <- suppressWarnings(
    fasta_smry <- dplyr::bind_cols(
      prot_acc = purrr::map_chr(fasta_db, attr, "prot_acc"), 
      prot_desc = purrr::map_chr(fasta_db, attr, "prot_desc"), 
      prot_mass = purrr::map_dbl(fasta_db, ~ calc_avg_pep(.x)), 
      prot_len = nchar(fasta_db)
    ) %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::mutate(prot_mass = round(prot_mass, digits = 0))
  ) %>% 
    dplyr::left_join(acc_lookup, ., by = "prot_acc")
  
  # -- add columns gene, organism, species, entrez ---
  df_splits <- df %>% split(., .$acc_type, drop = TRUE)
  
  acc_lookup <- purrr::map(df_splits, ~ {
    acc_type <- .x[["acc_type"]] %>% .[1]
    
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
        
        # columns gene and organism
        acc_lookup <- genes %>% 
          dplyr::filter(grepl("OS=", .$prot_desc)) %>% 
          dplyr::mutate(organism = gsub("^.*OS=(.*?)=.*$", "\\1", prot_desc)) %>% 
          dplyr::mutate(organism = gsub("\\s\\S*$", "", organism)) %>% 
          dplyr::bind_rows(., na_org) %>% 
          dplyr::select(-prot_desc) %>% 
          dplyr::right_join(acc_lookup, by = "prot_acc")
        
        # column species
        acc_lookup <- acc_lookup %>% 
          dplyr::mutate(species = my_lookup[.$organism]) %>% 
          na_species_by_org() 
      })
      
      # column entrez
      if (is.null(entrez)) {
        acc_lookup <- add_entrez(acc_lookup, acc_type, warns)
      } else {
        acc_lookup <- add_custom_entrez(acc_lookup, entrez, acc_type)
      }
      
      # added on 12-01-2020
      acc_lookup <- acc_lookup %>% dplyr::filter(!is.na(.[["uniprot_acc"]]))
    } else if (acc_type == "refseq_acc") {
      # columns organism and species
      acc_lookup <- acc_lookup %>% 
        dplyr::filter(grepl("^[XNY]{1}[MRP]{1}_", prot_acc)) %>% 
        dplyr::mutate(organism = gsub("^.*\\[(.*)\\]\\.*.*", "\\1", prot_desc)) %>% 
        dplyr::mutate(species = my_lookup[.$organism]) %>% 
        na_species_by_org()
      
      # column entrez
      if (is.null(entrez)) {
        acc_lookup <- add_entrez(acc_lookup, acc_type, warns)
      } else {
        acc_lookup <- add_custom_entrez(acc_lookup, entrez, acc_type)
      }
      
      # separately column gene
      acc_lookup <- add_refseq_gene(acc_lookup, acc_type)
      
      # added on 12-01-2020
      acc_lookup <- acc_lookup %>% dplyr::filter(!is.na(.[["refseq_acc"]]))
    } else {
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
  
  return(acc_lookup)
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
annotPrn <- function (df, fasta, entrez) {
  na_fasta_name_by_prot_acc <- function(df) {
    if (nrow(df) > 0 && "fasta_name" %in% names(df)) {
      na_fasta_name <- 
        (is.na(df$fasta_name)) | (stringr::str_length(df$fasta_name) == 0)
      df$fasta_name[na_fasta_name] <- df$prot_acc[na_fasta_name]
    }
    
    return(df)
  }
  
  na_acc_type_to_other <- function(df) {
    if (nrow(df) > 0 && "acc_type" %in% names(df)) {
      na_acc_type <- 
        (is.na(df$acc_type)) | (stringr::str_length(df$acc_type) == 0)
      df$acc_type[na_acc_type] <- "other_acc"
    }

    return(df)
  }
  
  acc_lookup <- parse_fasta(df, fasta, entrez)
	
	acc_lookup <- dplyr::bind_cols(
	  acc_lookup %>% 
	    dplyr::select(prot_acc), 
	  acc_lookup %>% 
	    dplyr::select(-which(names(.) %in% names(df))) %>% 
	    dplyr::select(-"organism"), 
	) 
	
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
	  reloc_col_after("prot_mass", "prot_acc")

	return(df)
}


#' Adds kinase annotation
#'
#' @inheritParams info_anal
#' @param acc_type Character string; the type of protein accessions in one of
#'   c("refseq_acc", "uniprot_acc", "uniprot_id")
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
annotKin <- function (df, acc_type) {
	if (is.null(df)) return(df)
  
  stopifnot ("prot_acc" %in% names(df))
	
  data(package = "proteoQ", kinase_lookup)
  
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

	return(df)
}


#' Saves the arguments in a function call
#'
#' @param call_pars Language.
#' @param fn The name of function being saved.
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
save_call <- function(call_pars, fn) {
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
match_fmls <- function(formulas) {
  dat_dir <- get_gl_dat_dir()
  fml_file <-  file.path(dat_dir, "Calls/pepSig_formulas.rda")
  
  if (file.exists(fml_file)) {
    load(file = fml_file)
  } else {
    stop("Run `pepSig()` first.")
  }
  
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
match_call_arg <- function (call_rda = "foo", arg = "scale_log2r") {
  dat_dir <- get_gl_dat_dir()
  
  call_rda <- rlang::as_string(rlang::enexpr(call_rda))
  arg <- rlang::as_string(rlang::enexpr(arg))
  
  rda <- paste0(call_rda, ".rda")
  file <- file.path(dat_dir, "Calls", rda)
  if (!file.exists(file)) stop(rda, " not found.")
  
  load(file = file)
  if (!arg %in% names(call_pars)) 
    stop(arg, " not found in the latest call to ", call_rda, 
         call. = FALSE)
  
  call_pars[[arg]]
}


#' Matches the name of GSPA result file (not currently used)
#'
#' @param anal_type Always \code{GSPA}; maybe different value for future uses.
#' @param subdir Character string of sub directory
#' @inheritParams prnHM
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
match_gspa_filename <- function (anal_type = "GSPA", subdir = NULL, 
                                 scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(!is.null(subdir))

  dat_dir <- get_gl_dat_dir()
  
  if (anal_type == "GSPA") {
    file <- file.path(dat_dir, "Calls/anal_prnGSPA.rda")
    if (!file.exists(file)) stop("Run `prnGSPA` first.", call. = FALSE)
    load(file = file)
  }
  
  filename <- call_pars$filename
  if (rlang::is_empty(filename)) {
    filename <- list.files(path = file.path(dat_dir, "Protein/GSPA", subdir), 
                           pattern = "^Protein_GSPA_.*\\.txt$", 
                           full.names = FALSE)
    
    if (rlang::is_empty(filename)) {
      stop("No result file found under `", subdir, "`", call. = FALSE)
    }
    
    if (scale_log2r) {
      filename <- filename %>% .[grepl("^Protein_GSPA_Z", .)]
    } else {
      filename <- filename %>% .[grepl("^Protein_GSPA_N", .)]
    }

    if (impute_na) {
      filename <- filename %>% .[grepl("Protein_GSPA_[NZ]_impNA\\.txt", .)]
    } else {
      filename <- filename %>% .[grepl("Protein_GSPA_[NZ]\\.txt", .)]
    }

    if (length(filename) > 1) {
      stop("More than one result file found under `", subdir, "`", call. = FALSE)
    }
    
    if (rlang::is_empty(filename)) 
      stop("No input files correspond to impute_na = ", 
           impute_na, ", scale_log2r = ", scale_log2r, 
           call. = FALSE)
  }

  return(filename)
}


#' Matches gset_nms to prnGSPA (not currently used)
#' 
#' @inheritParams prnGSPA
match_gset_nms <- function (gset_nms = NULL) {
  dat_dir <- get_gl_dat_dir()
  file <- file.path(dat_dir, "Calls/anal_prnGSPA.rda")
  
  if (is.null(gset_nms)) {
    if (file.exists(file)) {
      gset_nms <- match_call_arg(anal_prnGSPA, gset_nms)
    } else {
      stop("`gset_nms` is NULL and ", file, " not found.", 
           call. = FALSE)
    }
  } else {
    if (file.exists(file)) {
      gset_nms <- gset_nms %>% 
        .[. %in% match_call_arg(anal_prnGSPA, gset_nms)]
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
  
  if (purrr::is_empty(gset_nms)) {
    stop ("The `gset_nms` is EMPTY after parameter matching.", 
          "\n\tCheck the values of `gset_nms` in the latest call to `prnGSPA(...)`.", 
          call. = FALSE)
  }
  
  return(gset_nms)
}


#' Match scale_log2r
#'
#' \code{match_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in caller environment.
#' 
#' @inheritParams prnHist
match_scale_log2r <- function(scale_log2r) {
  stopifnot(rlang::is_logical(scale_log2r))
  
  global_var <-tryCatch(global_var <-get("scale_log2r", envir = .GlobalEnv),
                        error = function(e) "e")
  if (global_var != "e" & is.logical(global_var)) scale_log2r <- global_var
  
  return(scale_log2r)
}


#' Match to a global logical variable
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
match_logi_gv <- function(var, val) {
  oval <- val
  gvar <- tryCatch(gvar <- get(var, envir = .GlobalEnv), 
                   error = function(e) "e")
  
  if (gvar != "e") {
    stopifnot(rlang::is_logical(gvar))
    if (gvar != oval) {
      warning("Coerce ", var, " to ", gvar, 
              " after matching to the global setting.", 
              call. = FALSE)
    }
    return(gvar)
  } else {
    return(val)
  }
}


#' Match scale_log2r 
#'
#' \code{match_prnSig_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in the most recent prnSig at a given impute_na status
#' 
#' @inheritParams prnHM
match_prnSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  if (impute_na) {
    mscale_log2r <- match_call_arg(prnSig_impTRUE, scale_log2r)
  } else {
    mscale_log2r <- match_call_arg(prnSig_impFALSE, scale_log2r)
  }
  
  if (scale_log2r != mscale_log2r) {
    warning("scale_log2r = ", mscale_log2r, 
            " after matching to prnSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }
  
  scale_log2r <- mscale_log2r
}


#' Match scale_log2r 
#'
#' \code{match_pepSig_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in the most recent pepSig at a given impute_na status.
#' @inheritParams prnHM
match_pepSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  if (impute_na) {
    mscale_log2r <- match_call_arg(pepSig_impTRUE, scale_log2r)
  } else {
    mscale_log2r <- match_call_arg(pepSig_impFALSE, scale_log2r)
  }
  
  if (scale_log2r != mscale_log2r) {
    warning("scale_log2r = ", mscale_log2r, 
            " after matching to pepSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }
  
  scale_log2r <- mscale_log2r
}


#' Replaces NA genes (not currently used)
#' 
#' @inheritParams info_anal
#' @inheritParams annotKin
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
replace_na_genes <- function(df, acc_type) {
	acc_type <- tolower(acc_type)

	if (acc_type == "refseq_acc") {
		na_gene <- (is.na(df[, c("gene")])) | (str_length(df$gene) == 0)
		df$gene <- as.character(df$gene)
		df$gene[na_gene] <- as.character(df$prot_acc[na_gene])
	} else if (acc_type == "uniprot_id") {
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
		rm(temp)

		df <- df %>%
			dplyr::mutate(gene = gsub(".*\\|", "", gene))
	}

	return(df)
}


#' Replaces NA genes
#' 
#' @param df A data frame.
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
na_genes_by_acc <- function(df) {
  stopifnot("prot_acc" %in% names(df), "gene" %in% names(df))
  
  na_genes <- (is.na(df$gene)) | (stringr::str_length(df$gene) == 0)
  df$gene <- as.character(df$gene)
  df$gene[na_genes] <- as.character(df$prot_acc[na_genes])
  
  invisible(df)
}


#' Find peptide start and end positions
#'
#' \code{find_pep_pos} finds the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}.
#' 
#' @param fasta_name The fasta name
#' @param pep_seq Peptide sequence
#' @param fasta The database of fasta
#' @import dplyr purrr stringr tidyr
#' @importFrom magrittr %>% %T>% %$% %<>% 
find_pep_pos <- function (fasta_name, pep_seq, fasta) {
  this_fasta <- fasta %>% .[names(.) == fasta_name]
  pep_seq <- as.character(pep_seq)
  
  if (!rlang::is_empty(this_fasta)) {
    pep_pos <- stringr::str_locate(this_fasta, pattern = pep_seq)
    
    pos_bf <- pep_pos[1] - 1
    pos_af <- pep_pos[2] + 1
    
    pep_res_before <- stringr::str_sub(this_fasta, pos_bf, pos_bf)
    pep_res_after <- stringr::str_sub(this_fasta, pos_af, pos_af)

    # Mascot specialty of "X" residues (see also aa_residues.rda)
    # prot_acc: "XP_003960355", original "QERFCQXK" becomes "QERFCQVK"
    if (any(is.na(c(pep_res_before, pep_res_after)))) {
      pep_pos <- cbind(pep_seq, pep_res_before = NA, start = NA, end = NA, 
                       pep_res_after = NA, fasta_name = fasta_name, is_tryptic = NA)
      
      return(pep_pos)
    }
    
    if (any(nchar(pep_res_before) == 0)) pep_res_before <- "-"
    if (any(nchar(pep_res_after) == 0)) pep_res_after <- "-"
    
    # ADVSLPSMQGDLK|NP_612429: not "E.ADVSLPSMQGDLK.T" but "K.ADVSLPSMQGDLK.T"
    if (any(pep_res_before %in% c("K", "R", "-"))) { # the first match is tryptic
      pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, 
                       pep_res_after, fasta_name, is_tryptic = TRUE)
    } else if (pep_res_before == "M" && pep_pos[1] == 2) { # the first match also tryptic
      pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, 
                       pep_res_after, fasta_name, is_tryptic = TRUE)
    } else { # the first match is non-tryptic
      pep_seq_new <- paste0(c("K", "R"), pep_seq)
      pep_pos_new_all <- purrr::map(pep_seq_new, ~ stringr::str_locate(this_fasta, .x))
      ok_pos <- purrr::map_lgl(pep_pos_new_all, ~ !is.na(.x[[1]]))
      
      if (sum(ok_pos) > 0) { # tryptic match existed
        pep_pos_new <- pep_pos_new_all[[which(ok_pos)[1]]]
        
        pos_bf_new <- pep_pos_new[1]
        pos_af_new <- pep_pos_new[2] + 1
        
        pep_res_before_new <- stringr::str_sub(this_fasta, pos_bf_new, pos_bf_new)
        pep_res_after_new <- stringr::str_sub(this_fasta, pos_af_new, pos_af_new)
        
        pep_pos_new[1] <- pep_pos_new[1] + 1
        
        pep_pos <- cbind(pep_seq, pep_res_before_new, pep_pos_new, 
                         pep_res_after_new, fasta_name, is_tryptic = TRUE)        
      } else { # no tryptic matches
        pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, 
                         pep_res_after, fasta_name, is_tryptic = FALSE)
      }
    }
  } else { # no fasta matches
    pep_pos <- cbind(pep_seq, pep_res_before = NA_character_, 
                     start = NA_integer_, end = NA_integer_, 
                     pep_res_after = NA_character_, fasta_name = fasta_name, 
                     is_tryptic = FALSE)
  }
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
annotPeppos <- function (df){
  stopifnot(all(c("fasta_name", "pep_seq") %in% names(df)))
  
  load(file = file.path(dat_dir, "fasta_db.rda"))

  # ok cases that same `pep_seq` but different `prot_acc`
  # (K)	MENGQSTAAK	(L) NP_510965
  # (-)	MENGQSTAAK	(L) NP_001129505  
  
  df <- df %>% 
    dplyr::mutate(pep_prn = paste(pep_seq, fasta_name, sep = "@"))
  
  df_pep_prn <- df %>% 
    dplyr::filter(!duplicated(pep_prn)) %>% 
    dplyr::select(c("pep_seq", "fasta_name")) 
  
  fasta_db_sub <- fasta_db %>% .[names(.) %in% unique(df_pep_prn$fasta_name)]

  pep_pos_all <- suppressWarnings(
    purrr::map2(as.list(df_pep_prn$fasta_name), 
                as.list(df_pep_prn$pep_seq), 
                find_pep_pos, 
                fasta_db_sub)
  ) %>% 
    do.call(rbind, .) %>% 
    `colnames<-`(c("pep_seq", "pep_res_before", "pep_start", 
                   "pep_end", "pep_res_after", 
                   "fasta_name", "is_tryptic")) %>% 
    data.frame(check.names = FALSE) %>% 
    tidyr::unite(pep_prn, pep_seq, fasta_name, sep = "@", remove = TRUE) %>% 
    dplyr::mutate(pep_start = as.numeric(pep_start), pep_end = as.numeric(pep_end))
  
  df$pep_res_before <- NULL
  df$pep_res_after <- NULL
  df$pep_start <- NULL
  df$pep_end <- NULL
  
  if ("pep_res_before" %in% names(df)) pep_pos_all$pep_res_before <- NULL
  if ("pep_res_after" %in% names(df)) pep_pos_all$pep_res_after <- NULL
  if ("pep_start" %in% names(df)) pep_pos_all$pep_start <- NULL
  if ("pep_end" %in% names(df)) pep_pos_all$pep_end <- NULL
  
  df <- df %>% 
    dplyr::left_join(pep_pos_all, by = "pep_prn") %>% 
    dplyr::select(-pep_prn) %>% 
    reloc_col_before("pep_end", "pep_res_before") %>% 
    reloc_col_before("pep_start", "pep_end")
}


#' Subset fasta by accession type
#' 
#' @inheritParams info_anal
#' @inheritParams normPSM
#' @inheritParams annotKin
#' @import dplyr purrr  
#' @importFrom magrittr %>% %T>% %$% %<>% 
subset_fasta <- function (df, fasta, acc_type) {
  stopifnot("prot_acc" %in% names(df))
  
  if (! acc_type %in% c("refseq_acc", "uniprot_id", "uniprot_acc")) {
    stop("The type of protein accesion needs to one of ", 
         "\'uniprot_id\', \'uniprot_acc\' or \'refseq_acc\'.",
         call. = FALSE)
  }
  
  fasta <- purrr::map(fasta, ~ read_fasta(.x)) %>% do.call(`c`, .)
  
  if (acc_type == "uniprot_id") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|.*\\|(.*)$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } else if (acc_type == "uniprot_acc") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } else if (acc_type == "refseq_acc") {
    fasta <- fasta %>% 
      .[names(.) %in% unique(df$prot_acc)]
  }    
}


#' Find peptide abundance indexes.
#' 
#' @param max_len The maximum length of peptide sequence for consideration.
#' @inheritParams info_anal
#' @inheritParams find_shared_prots
add_prot_icover <- function (df, id = "gene", pep_id = "pep_seq", 
                             max_len = 50) {
  dat_dir <- get_gl_dat_dir()
  
  ok <- tryCatch(load(file = file.path(dat_dir, "fasta_db.rda")),
                 error = function(e) "e")
  if (ok != "fasta_db") {
    stop("`fasta_db.rda` not found under ", dat_dir, ".", 
         call. = FALSE)
  }
  
  fasta_db <- fasta_db %>% .[!duplicated(names(.))]
  
  id <- rlang::as_string(rlang::enexpr(id))
  if (id == "gene") {
    gn_rollup <- TRUE
    id <- "prot_acc"
  } else {
    gn_rollup <- FALSE
  }
  
  min_len <- min(nchar(df[[pep_id]])) - 1
  max_len2 <- max_len - 1
  
  zero_npeps <- df %>%
    dplyr::filter(pep_len <= max_len, pep_miss == 0) %>% 
    dplyr::select(c(pep_id, "fasta_name")) %>%
    dplyr::filter(!duplicated(!!rlang::sym(pep_id))) %>% 
    dplyr::group_by(fasta_name) %>%
    dplyr::summarise(prot_n_pep0 = n())
  
  max_npeps <- fasta_db %>% 
    purrr::map_int(~ {
      x <- stringr::str_split(.x, "[KR]")
      
      x %>% purrr::map_int(~ {
        len <- stringr::str_length(.x)
        sum(len >= min_len & len <= max_len2)
      })
    }) %>% 
    tibble::tibble(
      fasta_name = names(.), 
      prot_n_pepi = .,
    ) %>% 
    dplyr::left_join(zero_npeps, by = "fasta_name") %>% 
    dplyr::mutate(prot_icover = prot_n_pep0/prot_n_pepi, 
                  prot_icover = round(prot_icover, digits = 3)) %>% 
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
  } else {
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
calc_cover <- function(df, id) {
  dat_dir <- get_gl_dat_dir()
  
  stopifnot(all(c("prot_acc", "gene", "pep_start", "pep_end") %in% names(df)))
  
  df$prot_cover <- NULL
  
  load(file = file.path(dat_dir, "fasta_db.rda"))

  if (all(is.factor(df$pep_start))) {
    df$pep_start <- df$pep_start %>% as.character() %>% as.numeric()
  }
    
  if (all(is.factor(df$pep_end))) {
    df$pep_end <- df$pep_end %>% as.character() %>% as.numeric()
  }

  id <- rlang::as_string(rlang::enexpr(id))
  if (id == "gene") {
    gn_rollup <- TRUE
    id <- "prot_acc"
  } else {
    gn_rollup <- FALSE
  }
  
  load(file = file.path(dat_dir, "acc_lookup.rda"))
  
  if (length(fasta_db) == 0) {
    stop("No fasta entries match protein accessions.", 
         "Check the correctness of fasta file(s).", 
         call. = FALSE)
  }
  
  if (length(fasta_db) <= 200) {
    warning("Less than 200 entries in fasta matched by protein accessions.\n", 
            "Make sure the fasta file is correct.", 
            call. = FALSE)
  }
  
  df_sels <- df %>%
    dplyr::select(prot_acc, pep_start, pep_end) %>%
    dplyr::mutate(index = row_number()) %>% 
    dplyr::left_join(acc_lookup, by = "prot_acc") %>%
    dplyr::filter(!is.na(prot_len), !duplicated(index)) %>% 
    dplyr::select(-index)
  
  if (nrow(df_sels) == 0) {
    stop("Probably incorrect accession types in the fasta file(s).", 
         call. = FALSE)
  }
  
  df_sels <- df_sels %>%
    dplyr::filter(pep_start <= prot_len) %>%
    dplyr::filter(pep_end <= prot_len) %>%
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
    df_sels <- df %>% 
      dplyr::select(prot_acc, gene) %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_sels, by = "prot_acc") %>% 
      dplyr::select(-prot_acc) %>% 
      dplyr::filter(!is.na(gene)) %>% 
      dplyr::group_by(gene) 
    
    df_sels <- suppressWarnings(
      df_sels %>% dplyr::summarise_all(~ max(.x, na.rm = TRUE))
    )
    
    df_sels <- df %>% 
      dplyr::select(prot_acc, gene) %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_sels, by = "gene") %>% 
      dplyr::select(-gene) 
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
                    prot_icover = round(prot_icover, digits = 3))
  }
  
  df <- df %>% 
    dplyr::mutate(prot_cover = round(prot_cover, digits = 3)) 

  invisible(df)
}


#' Converts log2FC to linear fold changes
#' 
#' @inheritParams info_anal
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
to_linfc <- function(df) {
		nms <- rownames(df)

		df %>%
			purrr::map(~ {ifelse(.x > 0, 2^.x, -1/(2^.x))}) %>%
			data.frame(check.names = FALSE) %>%
			`rownames<-`(nms)
}


#' Remove single-value columns
#' 
#' @inheritParams setHMlims
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
rm_sglval_cols <- function (x) {
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


#' Combine data with metadata
#' 
#' @param data A data frame
#' @param metadata Another data frame 
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
cmbn_meta <- function(data, metadata) {
  data <- data %>% 
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::left_join(metadata, by = "Sample_ID") %>%
    dplyr::mutate_at(vars(one_of("Color", "Fill", "Shape", "Size", "Alpha")), 
                     ~ as.factor(.)) 
  
  col_nms <- data %>% 
    .[!not_all_NA(.)] %>% 
    names()
  
  if (!purrr::is_empty(col_nms)) {
    warning("columns with all-NA values for aesthetics excluded: \n", 
            purrr::reduce(col_nms, paste, sep = ", "), ".\n", 
            call. = FALSE)
  }
  
  data %>%
    dplyr::select(which(not_all_NA(.))) %>% 
    rm_sglval_cols()
}


#' Check file names for ggsave
#' @param filename Character string; An output file name.
gg_imgname <- function(filename) {
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


#' Check file names for ggsave
#' 
#' @inheritParams info_anal
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>% 
rm_pval_whitespace <- function(df) {
  df <- df %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.character) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ gsub("\\s*", "", .x) ) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ suppressWarnings(as.numeric(.x)))
}


#' Filter rows
#'
#' @param df a data frame. 
#' @param ... Arguments for \link[dplyr]{filter}
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
filters_in_call <- function (df, ...) {
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
    
    if (sum(row_vals) == 0) {
      stop("Zero row of data available after `filter_` vararg(s).\n", 
           "Examine both (a) the values under the selected column(s) and ", 
           "(b) the supplied logical condition(s).\n", 
           "For example, values are not all NAs or FALSE under the column(s).", 
           call. = FALSE)
    }
    
    df <- df[row_vals, , drop = FALSE]
  }
  
  return(df)
}


#' Arrange rows
#'
#' @param .df a data frame. 
#' @param .na.last The same as \link[base]{order}
#' @param ... Arguments for \link[dplyr]{arrange}
#'
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
arrangers_in_call <- function(.df, ..., .na.last = TRUE) {
  dots <- rlang::enexprs(...)
  nms <- names(dots)
  
  for (i in seq_along(dots)) {
    row_orders <- dots[[nms[i]]] %>% 
      rlang::eval_bare()
    
    if (!rlang::is_list(row_orders)) row_orders <- list(row_orders)
    
    row_orders <- purrr::map(row_orders, ~ {
      is_char <- is.character(.x)
      if (is_char) .x <- rlang::sym(.x)
      
      return(.x)
    })
    
    order_call <- rlang::expr(order(!!!row_orders, na.last = !!.na.last))
    ord <- rlang::eval_tidy(order_call, .df)
    stopifnot(length(ord) == nrow(.df))

    .df <- .df[ord, , drop = FALSE]
  }
  
  return(.df)
}


#' Calculate PSM SDs
#' 
#' The standard deviations are based on samples under each TMT set and LCMS.
#' 
#' @inheritParams info_anal
#' @inheritParams standPep
#' @inheritParams channelInfo
#' @inheritParams calcPeptide
calc_sd_fcts_psm <- function (df, range_log2r = c(5, 95), range_int = c(5, 95), 
                              set_idx, injn_idx) {
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


#' Calculate CV per TMT_Set and LCMS_injection
#' 
#' @inheritParams info_anal
#' @param type Character string; the type of data.
calcSD_Splex <- function (df, id, type = "log2_R") {
  if (type == "log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^log2_R[0-9]{3}", names(.)))
  } else if (type == "N_log2_R") {
    df <- df %>% 
      dplyr::select(!!rlang::sym(id), grep("^N_log2_R[0-9]{3}", names(.)))
  } else if (type == "Z_log2_R") {
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
                      theme = NULL, ...) {

  df <- df %>% 
    dplyr::select(which(not_all_NA(.)))
  
  if (ncol(df) == 0) {
    message("No SD columns available.")
    return (NULL)
  }
  
  dat_dir <- get_gl_dat_dir()
  
  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  col_select <- ifelse(is.null(col_select), rlang::expr(Select), rlang::sym(col_select))
  col_order <- ifelse(is.null(col_order), rlang::expr(Order), rlang::sym(col_order))
  
  if (col_select == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
  if (col_order == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  if (is.null(label_scheme_full[[col_select]])) {
    stop("Column \'", rlang::as_string(col_select), "\' does not exist.", 
         call. = FALSE)
  } else if (sum(!is.na(label_scheme_full[[col_select]])) == 0) {
    stop("No samples were selected under column \'", rlang::as_string(col_select), "\'.",
         call. = FALSE)
  }
  
  if (is.null(label_scheme_full[[col_order]])) {
    warning("Column \'", rlang::as_string(col_order), "\' does not exist.
			Samples will be arranged by the alphebatic order.", 
            call. = FALSE)
  } else if (sum(!is.na(label_scheme_full[[col_order]])) == 0) {
    # warning("No orders under column \'", rlang::as_string(col_order), "\'.", call. = FALSE)
  }
  
  if (is_psm) {
    label_scheme_sub <- local({
      set_idx <- filepath %>% 
        gsub("^.*TMTset(\\d+)_.*", "\\1", .) %>% 
        as.integer()
      
      injn_idx <- filepath %>% 
        gsub("^.*TMTset\\d+_LCMSinj(\\d+)_sd\\.png$", "\\1", .) %>% 
        as.integer()
    
      label_scheme_full %>% 
        dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx) 
    })
  } else {
    label_scheme_sub <- label_scheme
  }
  
  label_scheme_sub <- label_scheme_sub %>% 
    dplyr::mutate(new_id = paste0(TMT_Channel, " (", Sample_ID, ")")) %>% 
    dplyr::mutate(new_id = gsub("TMT-", "", new_id)) %>% 
    dplyr::select(Sample_ID, TMT_Set, new_id, !!col_select, !!col_order) %>%
    dplyr::filter(!is.na(!!col_select))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  if (TMT_plex == 0) {
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
  
  df <- df %>% dplyr::filter(!duplicated(.[[id]]))

  if (type == "log2_R") {
    df_sd <- df %>% dplyr::select(id, grep("^sd_log2_R[0-9]{3}[NC]*", names(.)))
  } else if (type == "N_l og2_R") {
    df_sd <- df %>% dplyr::select(id, grep("^sd_N_log2_R[0-9]{3}[NC]*", names(.)))
  } else if (type == "Z_log2_R") {
    df_sd <- df %>% dplyr::select(id, grep("^sd_Z_log2_R[0-9]{3}[NC]*", names(.)))
  }
  
  # all-NA first removed for finding all-NaN columns
  df_sd <- df_sd %>% 
    dplyr::filter(rowSums(!is.na(.[grep("^.*log2_R[0-9]{3}", names(.))])) > 0)
  
  if (adjSD) {
    SD <- df %>%
      dplyr::select(grep("^log2_R[0-9]{3}|^I[0-9]{3}", names(.))) %>%
      dblTrim(., range_log2r = c(0, 100), range_int = c(0, 100), 
              type_r = "log2_R", type_int = "I")
    
    df_sd[, grep("^.*log2_R", names(df_sd))] <- 
      df_sd[, grep("^.*log2_R", names(df_sd)), drop = FALSE] %>% 
      sweep(., 2, sqrt(SD), "/")
    
    df_z <- df_sd %>% dplyr::select(grep("^.*log2_R[0-9]{3}", names(.)))
    nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
    df_z[, nan_cols] <- 0
    df_sd[, grep("^.*_log2_R[0-9]{3}", names(df_sd))] <- df_z
    
    rm(df_z, nan_cols, SD)
  }

  if (TMT_plex > 0) {
    df_sd <- df_sd %>% `names<-`(gsub("^.*log2_R", "", names(.))) 
  } else {
    df_sd <- df_sd %>% 
      `names<-`(gsub("sd_log2_R000 \\((.*)\\)$", "\\1", names(.))) 
  }
  
  if (!is_psm) {
    df_sd <- df_sd %>% 
      dplyr::select(id, which(names(.) %in% label_scheme_sub$new_id))
  }

  Levels <- names(df_sd) %>% .[! . %in% id]

  if (!purrr::is_empty(Levels)) {
    df_sd <- df_sd %>%
      tidyr::gather(key = !!rlang::sym(id), value = "SD") %>%
      dplyr::rename(Channel := !!rlang::sym(id)) %>% 
      dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
      dplyr::filter(!is.na(SD))
    
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
           y = expression("SD ("*log[2]*"FC)")) 
    
    if (!is.null(ymax)) {
      if (is.null(ybreaks)) {
        ybreaks <- ifelse(ymax > 1, 0.5, ifelse(ymax > 0.5, 0.2, 0.1))
      }
      p <- p + scale_y_continuous(limits = c(0, ymax), 
                                  breaks = seq(0, ymax, ybreaks))
    }

    if (is.null(theme) || purrr::is_function(theme)) {
      theme <- theme_psm_violin
    }
    p <- p + theme

    if (flip_coord) {
      p <- p + coord_flip()
      width_temp <- width
      width <- height
      height <- width_temp
      rm(width_temp)
    }
    
    ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot", "width", "height", 
                                           "limitsize"))
    
    my_call <- rlang::expr(ggplot2::ggsave(filename = !!filepath, plot = !!p, 
                                           width = !!width, height = !!height, 
                                           limitsize = FALSE, 
                                           !!!ggsave_dots))

    suppressWarnings(try(eval(my_call, rlang::caller_env())))
  }
}


#' Violin plots of reporter-ion intensity per TMT_Set and LCMS_injection
#' 
#' @inheritParams info_anal
#' @inheritParams sd_violin
rptr_violin <- function(df, filepath, width, height) {
  df <- df %>% 
    dplyr::select(which(not_all_NA(.)))
  
  if (ncol(df) == 0) {
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


#' geometric mean
#' 
#' @param x A data frame.
#' @param ... The same in \code{mean}.
my_geomean <- function (x, ...) {
  x <- log10(x) %>% mean(...)
  10^x
}


#' phospho counts
#' 
count_phosphopeps <- function() {
  dat_dir <- get_gl_dat_dir()
  
  df <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
  
  id <- match_call_arg(normPSM, group_psm_by)
  
  df_phos <- df %>% dplyr::filter(grepl("[sty]", .[[id]]))
  
  n_phos_peps <- nrow(df_phos)
  n_phos_sites <- stringr::str_count(df_phos[[id]], "[sty]") %>% sum()
  
  write.csv(
    data.frame(n_peps = n_phos_peps, n_sites = n_phos_sites), 
    file.path(dat_dir, "Peptide/cache", "phos_pep_nums.csv"), 
    row.names = FALSE
  )
}


#' peptide mis-cleavage counts
#' 
count_pepmiss <- function() {
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "PSM/cache"), recursive = TRUE, showWarnings = FALSE)
  
  rmPSMHeaders()
  
  filelist <- list.files(path = file.path(dat_dir, "PSM/cache"),
                        pattern = "^F[0-9]{6}_hdr_rm.csv$")
  
  if (length(filelist) == 0) stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
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
  }) %>% do.call(rbind, .)
  
  write.csv(df, file.path(dat_dir, "PSM/cache/miscleavage_nums.csv"), 
            row.names = FALSE)
}


#' Row filtration helpers
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
contain_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(match, vars, fixed = TRUE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{contain_chars_in}: contain some of the characters in a literal string;
#' "PEPTIDES" contain_chars_in "XP".
#' 
#' @rdname contain_str
#' @export
contain_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{not_contain_str}" not contain a literal string; "PEPTIDES"
#' not_contain_str "TED".
#' 
#' @rdname contain_str
#' @export
not_contain_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(match, vars, fixed = TRUE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{not_contain_chars_in}: not contain any of the characters in a literal
#' string; "PEPTIDES" not_contain_chars_in  "CAB".
#' 
#' @rdname contain_str
#' @export
not_contain_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}

#' Row filtration helpers
#'
#' \code{start_with_str}: start with a literal string. "PEPTIDES" start_with_str
#' "PEP".
#' 
#' @rdname contain_str
#' @export
start_with_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("^", match), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{end_with_str}: end with a literal string. "PEPTIDES" end_with_str
#' "TIDES".
#' 
#' @rdname contain_str
#' @export
end_with_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0(match, "$"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{start_with_chars_in}: start with one of the characters in a literal
#' string. "PEPTIDES" start_with_chars_in "XP".
#' 
#' @rdname contain_str
#' @export
start_with_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("^[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{ends_with_chars_in}: end with one of the characters in a literal
#' string. "PEPTIDES" ends_with_chars_in "XS".
#' 
#' @rdname contain_str
#' @export
ends_with_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[", match, "]$"), vars, fixed = FALSE, ignore.case)
}


#' Row filtration helpers
#'
#' \code{rows_are_all}: rows are all
#' @rdname contain_str
rows_are_all <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}


#' Row filtration helpers
#'
#' \code{rows_are_all}: rows are all
#' @rdname contain_str
rows_are_not_all <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}


#' Concatenate formula(s) to varargs of dots
#' 
#' @param fmls A character vector of formula(s)
#' @param dots A character vector of formula(s) in \code{dots}
#' @param fml_nms A character vector containing the names of \code{fmls}.
#' @inheritParams info_anal
concat_fml_dots <- function(fmls = NULL, fml_nms = NULL, 
                            dots = NULL, anal_type = "zzz") {
  dat_dir <- get_gl_dat_dir()
  
  if ((!purrr::is_empty(fmls)) && (anal_type == "GSEA")) return(c(dots, fmls))
  
  if (purrr::is_empty(fmls)) {
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
  } else {
    match_fmls(fmls)
    dots <- c(dots, fmls)
  }
  
  return(dots)
}


#' Roll up genes
#' 
#' @param df A data frame
#' @param cols Column indexes
gn_rollup <- function (df, cols) {
  if (! "gene" %in% names(df)) return(df)
  
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
  } else {
    dfc <- df %>% 
      dplyr::select(gene) %>% 
      dplyr::mutate(prot_cover = NA)
  }
  
  df <- list(dfc, dfb, dfa) %>% 
    purrr::reduce(right_join, by = "gene") %>% 
    dplyr::filter(!is.na(gene), !duplicated(gene))
}


#' Compare dot-dot-dot between prior and current
#' 
#' @param call_nm The name of a function call.
#' @param curr_dots The values of the current dots.
#' @param pattern The pattern for comparison.
identical_dots <- function(call_nm, curr_dots, pattern) {
  dat_dir <- get_gl_dat_dir()
  
  file <- file.path(dat_dir, "Calls", paste0(call_nm, ".rda"))
  if (!file.exists(file)) return(FALSE)
  
  load(file = file)
  identical(call_pars %>% .[grepl(pattern, names(.))], curr_dots)
}


#' Complete cases among sample IDs in label_scheme_sub, not label_scheme
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
my_complete_cases <- function (df, scale_log2r, label_scheme_sub) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  # in case that reference(s) are included in label_scheme_sub, 
  # the values are zero for single ref or non-NA/NaN for multi ref 
  # and will not affect the complete.cases assessment
  rows <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% label_scheme_sub$Sample_ID)) %>% 
    complete.cases(.)
  
  if (sum(rows) == 0) {
    stop("None of the cases are complete.", call. = FALSE)
  }
  
  df <- df[rows, ]
}


#' my union of named list 
#' 
#' names will be kept after the unification
#' 
#' @param x A list of values.
#' @param y Another list of values.
my_union <- function (x, y) {
  x %>% 
    .[! names(.) %in% names(y)] %>% 
    c(y)
}


#' find the base names of files (not currently used)
#' 
#' @param filenames A character vector of filenames
find_fn_bases <- function (filenames) {
  gsub("\\.[^.]*$", "", filenames) # %>% .[1] 
}


#' find the extensions of files
#' 
#' @param filename A character string of filename
#' @param type the type of filename extension
find_fn_exts <- function (filename, type = "text") {
  purrr::map_chr(filename, ~ {
    if (!grepl("\\.", .x)) {
      fn_ext <- switch(
        text = "txt",
        graphics = "png",
        stop("Invalid extension in file name(s).", Call. = FALSE)
      )
    } else {
      fn_ext <- gsub("^.*\\.([^.]*)$", "\\1", .x) # %>% .[1]
    }
  })
}


#' check duplicate argments in 'dots'
#' 
#' @param blacklist A character vector of variable names.
#' @param ... A list of arguments for checking.
check_dots <- function (blacklist = NULL, ...) {
  dots <- rlang::enexprs(...)
  dups <- purrr::map_lgl(names(dots), ~ .x %in% blacklist)
  nms <- names(dots[dups])
  
  if (!rlang::is_empty(nms)) {
    stop("Do not use argument(s): ", purrr::reduce(nms, paste, sep = ", "), 
         call. = FALSE)
  }
}


#' check depreciated arguments
#' 
#' @param blacklist A list of paired character vectors. In each pair, the first
#'   element is the old argument name if incurred and the second is the new name
#'   for reminding.
#' @param ... A list of arguments for checking. The depreciated argument(s) are
#'   no longer in the formalArgs of a function. If present, they will be in
#'   \code{...}.
check_depreciated_args <- function (blacklist = NULL, ...) {
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


#' force 'complete_cases = TRUE' at 'impute_na = FALSE'
#' 
#' @inheritParams prnHM
to_complete_cases <- function (complete_cases = FALSE, impute_na = FALSE) {
  warn_msg1 <- "Coerce `complete_cases = TRUE` at `impute_na = FALSE`."
  
  if (!(impute_na || complete_cases)) {
    complete_cases <- TRUE
    rlang::warn(warn_msg1)
  }
  
  return(complete_cases)
}


#' check 'gset_nms'
#'
#' It is OK that `gset_nms` not in `c("go_sets", "c2_msig", "kegg_sets", "kinsub)` as
#' custom gene sets may be applied.
#' @inheritParams prnGSPA
check_gset_nms <- function (gset_nms) {
  customs <- gset_nms %>% .[! . %in% c("go_sets", "c2_msig", "kegg_sets", "kinsub")]
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


#' Check the suitability of existing param files for MGKernel
#' 
#' @param filepath A file path to \code{"MGKernel_params_N.txt"}
ok_existing_params <- function (filepath) {
  dat_dir <- get_gl_dat_dir()
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  if (!file.exists(filepath)) return(TRUE)
  
  params <- readr::read_tsv(file.path(filepath), 
                            col_types = cols(Component = col_double()))

  missing_samples <- local({
    par_samples <- params$Sample_ID %>% unique() %>% .[!is.na(.)]
    expt_samples <- label_scheme$Sample_ID %>% unique()
    setdiff(expt_samples, par_samples)
  })
  
  if (purrr::is_empty(missing_samples)) return(TRUE)
  
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


#' Find the environment of a function
#' 
#' This is from Hadley's Advanced R.
#' @param name The name of a function.
#' @param env The environment.
env_where <- function(name, env = rlang::caller_env()) {
  if (identical(env, rlang::empty_env())) {
    stop("Can't find ", name, call. = FALSE)
  } else if (rlang::env_has(env, name)) {
    env
  } else {
    env_where(name, rlang::env_parent(env))
  }
}


#' Generate cut points.
#' 
#' @param x A sequence of numeric values with probable NA value.
#' @inheritParams prnHist
set_cutpoints <- function (cut_points = NULL, x = NULL) {
  if (is.null(cut_points) && is.null(x)) {
    stop("`cut_points` and `x` can not be both NULL.", 
         call. = FALSE)
  }
    
  if (!is.null(x)) {
    oks <- x %>% .[!is.na(.)] %>% is.numeric() %>% all()
    
    if (!oks) {
      stop("All numeric `x` are required for setting cut points.", 
           call. = FALSE)
    }
  }
  
  if (!is.null(cut_points)) {
    oks <- cut_points %>% is.numeric() %>% all()
    
    if (!oks) {
      stop("All numeric values are required for `cut_points`.", 
           call. = FALSE)
    }
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
set_cutpoints2 <- function(cut_points, df) {
  if (is.null(cut_points) || is.infinite(cut_points)) {
    return(c(mean_lint = -Inf, mean_lint = Inf))
  }
  
  if (any(is.na(cut_points))) {
    cut_points <- cut_points[1]
    
    if (is.null(names(cut_points))) {
      cut_nm <- "mean_lint"
    } else {
      cut_nm <- names(cut_points)
    }
    
    if (! cut_nm %in% names(df)) {
      warning("Column `", cut_nm, "` not found in `df`;", 
              "instead use column `mean_lint`.", 
              call. = FALSE)
      cut_nm <- "mean_lint"
    }
    
    if (! is.numeric(df[[cut_nm]])) {
      stop("Column `", cut_nm, "` is not numeric.", 
           call. = FALSE)
    }
    
    cut_points <- quantile(df[[cut_nm]], na.rm = TRUE) 
    seqs <- seq_along(cut_points)
    names(cut_points)[seqs] <- cut_nm
  } else {
    cut_nm <- cut_points %>% 
      `[`(1) %>% 
      names() %>% 
      gsub("1$", "", .)
    
    if (purrr::is_empty(cut_nm)) cut_nm <- "mean_lint"
    
    if (! cut_nm %in% names(df)) {
      warning("Column `", cut_nm, "` not found in `df`;", 
              "instead use column `mean_lint`.", 
              call. = FALSE)
      cut_nm <- "mean_lint"
    }
    
    if (! is.numeric(df[[cut_nm]])) {
      stop("Column `", cut_nm, "` is not numeric.", call. = FALSE)
    }
    
    cut_points <- set_cutpoints(cut_points, df[[cut_nm]])
    seqs <- seq_along(cut_points)
    names(cut_points)[seqs] <- cut_nm
  }

  invisible(cut_points)
}


#' Load cRAPs
#' 
#' Load the cRAP proteins by the types of accessions
#' @param acc_types The types of protein accessions
load_craps <- function(acc_types) {
  data(package = "proteoQ", prn_annot_crap)
  
  if (identical(acc_types, "other_acc")) {
    craps <- prn_annot_crap %>% 
      dplyr::select("fasta_name") %>% 
      dplyr::filter(!is.na(fasta_name), !duplicated(fasta_name)) %>% 
      unlist()
  } else {
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
find_int_cols <- function (TMT_plex) {
  if (TMT_plex == 16) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C", 
                 "I132N", "I132C", "I133N", "I133C", "I134N")
  } else if (TMT_plex == 11) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_int <- NULL
  }
}


#' Finds the columns of reporter-ion ratios.
#' 
#' @inheritParams TMT_levels
find_ratio_cols <- function (TMT_plex) {
  if (TMT_plex == 16) {
    col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                   "R130N", "R130C", "R131N", "R131C", 
                   "R132N", "R132C", "R133N", "R133C", "R134N")
  } else if (TMT_plex == 11) {
    col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                   "R130N", "R130C", "R131N", "R131C")
  } else if (TMT_plex == 10) {
    col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                   "R130N", "R130C", "R131")
  } else if(TMT_plex == 6) {
    col_ratio <- c("R127", "R128", "R129", "R130", "R131")
  } else {
    col_ratio <- NULL
  }
}


#' Process MaxQuant mqpar.xml
#' 
#' @inheritParams load_expts
#' @param mqpar The name of .xml file. The default is "mqpar.xml".
#' @param filename The name of metadata file.
#' @param out The name of output .xml file.
make_mq_meta <- function (dat_dir = proteoQ:::get_gl_dat_dir(), mqpar = "mqpar.xml", 
                          filename = "mq_meta.txt", out = "new_mqpar.xml") {
  
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


#' Process MaxQuant mqpar.xml
#' 
#' @inheritParams load_expts
#' @param mqpar The name of .xml file. The default is "mqpar.xml".
#' @param filename The name of metadata file.
#' @param out The name of output .xml file.
#' @examples
#' \dontrun{
#' dat_dir <- "Y:/qzhang/MaxQuant"
#' make_mq_meta2()
#' }
make_mq_meta2 <- function (dat_dir = proteoQ:::get_gl_dat_dir(), mqpar = "mqpar.xml", 
                           filename = "mq_meta.txt", out = "new_mqpar.xml") {
  
  update_nodes <- function(mqxml, field = "filePaths", type = "string", 
                           vals = lookup$Name) {
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
  
  
  fn_suffix <- gsub("^.*\\.([^\\.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^\\.]*$", "", filename)
  
  if (fn_suffix %in% c("xls", "xlsx")) {
    lookup <- readxl::read_excel(file.path(dat_dir, filename))
  } else if (fn_suffix == "csv") {
    lookup <- readr::read_csv(file.path(dat_dir, filename))
  } else if (fn_suffix == "txt") {
    lookup <- readr::read_tsv(file.path(dat_dir, filename))
  } else {
    stop(filename, " needs to be '.xls' or '.xlsx'.")
  }
  
  lookup <- lookup %>% dplyr::filter(rowSums(!is.na(.)) > 0)

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
  
  
  mqxml <- readLines(file.path(dat_dir, mqpar)) %>% 
    update_nodes("filePaths", "string", lookup$Name) %>% 
    update_nodes("experiments", "string", lookup$Experiment) %>% 
    update_nodes("fractions", "short", lookup$Fraction) %>% 
    update_nodes("ptms", "boolean", lookup$PTM) %>% 
    update_nodes("paramGroupIndices", "int", lookup$Group) %>% 
    update_nodes("referenceChannel", "string", lookup$Reference) %T>% 
    writeLines(file.path(dat_dir, out))
}


#' Add numeric indices for peptides or proteins
#' 
#' @param df A data frame.
#' @param col_in The column key to be indexed.
#' @param col_out The name of the output column.
add_entry_ids <- function (df, col_in = "pep_seq", col_out = "pep_index") {
  if (col_out %in% names(df)) {
    return(reloc_col_after_last(df, col_out))
  }
  
  ids <- unique(df[[col_in]])
  
  seqs <- seq_along(ids) %>% `names<-`(ids)
  seqs <- tibble::tibble(!!col_in := names(seqs), !!col_out := seqs) 
  
  df %>% dplyr::left_join(seqs, by = col_in)
}


#' Check conflicts in function arguments.
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
check_formalArgs <- function(w = anal_prnTrend, f = cmeans, excludes = NULL) {
  w <- rlang::as_string(rlang::enexpr(w))
  f <- rlang::as_string(rlang::enexpr(f))
  
  args_w <- formalArgs(w) %>% .[! . == "..."]
  args_f <- formalArgs(f) %>% .[! . == "..."]

  if (!is.null(excludes)) {
    args_f <- args_f %>% .[! . %in% excludes]
  }
  
  if (!purrr::is_empty(args_f)) {
    purrr::walk(args_w, ~ {
      ambi <- stringr::str_detect(.x, paste0("^", args_f)) %>% 
        which()
      
      if (!purrr::is_empty(ambi)) {
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


#' Check conflicts in function arguments.
#' 
#' Conflicts between wrapper and function to be wrapped. Note that \code{f}
#' contains both function name and package name.
#' 
#' @inheritParams check_formalArgs
check_formalArgs2 <- function(w = anal_prnTrend, f = NMF::nmf, excludes = NULL) {
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


#' Set dot-dot-dot for ggsave.
#' 
#' @param dots Arguments passed on.
#' @param dummies Arguments automated for ggsave.
set_ggsave_dots <- function(dots, dummies = c("filename", "plot", "device", "path")) {
  if (purrr::is_empty(dots)) return(dots)
  
  args <- formalArgs(ggplot2::ggsave) %>% .[! . == "..."]
  
  dots %>% 
    .[names(.) %in% args] %>% 
    .[! names(.) %in% dummies]
}


#' Find the type of search engine
#' 
#' @inheritParams load_expts
find_search_engine <- function(dat_dir = NULL) {
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  
  pat_mascot <- "^F[0-9]+.*\\.csv$"
  pat_mq <- "^msms.*\\.txt$"
  pat_mf <- "^psm.*\\.tsv$"
  pat_sm <- "^PSMexport.*\\.ssv$"
  
  mascot <-list.files(path = file.path(dat_dir), pattern = pat_mascot) %>% 
    purrr::is_empty() %>% `!`
  mq <- list.files(path = file.path(dat_dir), pattern = pat_mq) %>% 
    purrr::is_empty() %>% `!`
  mf <- list.files(path = file.path(dat_dir), pattern = pat_mf) %>% 
    purrr::is_empty() %>% `!`
  sm <- list.files(path = file.path(dat_dir), pattern = pat_sm) %>% 
    purrr::is_empty() %>% `!`

  engines <- c(mascot = mascot, mq = mq, mf = mf, sm = sm)
  engine <- engines %>% .[.] %>% names()
  
  if (purrr::is_empty(engine)) {
    stop("No PSM files found.\n", 
         "File name(s) need to be in the format of \n",
         purrr::reduce(c(paste("Mascot", pat_mascot, sep = ": "), 
                         paste("MaxQuant", pat_mq, sep = ": "), 
                         paste("MSFragger", pat_mf, sep = ": "), 
                         paste("Spectrum Mill", pat_sm, sep = ": ")), 
                       paste, sep = " || "), 
         call. = FALSE)
  }
  
  if (length(engine) > 1) {
    stop("More than one data type found: ", 
         purrr::reduce(engine, paste, sep = ", "), ".", 
         call. = FALSE)
  }
  
  invisible(engine)
}


#' Check missing values in a ggplot2 object.
#' 
#' @param p A ggplot2 object.
check_ggplot_aes <- function(p) {
  q <- p %>% ggplot_build() %>% 
    `[[`(1) %>% 
    `[[`(1) %>% 
    dplyr::select(which(not_all_NA(.))) 
  
  cols <- q %>% 
    purrr::map_lgl(anyNA)
  
  if (any(cols)) {
    warning("NA values detected in the aesthetics of `", 
            purrr::reduce(names(which(cols)), paste, sep = ", "), "`.\n", 
            "Fix the missing values in the corresponding columns in `expt_smry.xlsx`.\n", 
            "Otherwise, may run into errors.\n", 
            call. = FALSE)
  }
}


#' Standardize range values.
#' 
#' @param from_to A numeric vector of length 2.
prep_range <- function(from_to = c(0, 1)) {
  stopifnot(length(from_to) == 2)
  
  stopifnot(from_to[2] > from_to[1], 
            from_to[1] >= 0,
            from_to[2] <= 100) 
  
  if (from_to[2] <= 1) from_to <- from_to * 100
  
  from_to
}


#' Find the deliminator of a file.
#' 
#' @param filename A file name.
find_delim <- function (filename) {
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
