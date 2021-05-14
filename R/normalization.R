#' Find sample IDS for used in the current fitting.
#' 
#' @param nm_a A name list from data table.
#' @param nm_b A name list from metadata.
find_fit_nms <- function(nm_a, nm_b) {
  ind <- purrr::map(nm_b, ~ grepl(paste0(" (", .x, ")"), nm_a, fixed = TRUE)) %>% 
    purrr::reduce(`|`) # .init = FALSE
  nm_a <- nm_a[ind]
}

#' SD for all non-trival samples.
#' 
#' @param df A data frame.
#' @param label_scheme Experiment summary
#' @inheritParams standPep
calc_sd_fcts <- function (df, range_log2r, range_int, label_scheme) {
  label_scheme_sd <- label_scheme %>%
    dplyr::filter(!Reference, !grepl("^Empty\\.", Sample_ID))	%>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = (.$Sample_ID)))
  
  SD <- df %>%
    dplyr::select(grep("^N_log2_R|^N_I", names(.))) %>%
    dblTrim(., range_log2r, range_int) %>%
    `names<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.)))
  
  cf_SD <- SD/mean(SD %>% .[names(.) %in% label_scheme_sd$Sample_ID], na.rm = TRUE)
  cf_SD <- cbind.data.frame(fct = cf_SD, SD) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
    dplyr::arrange(Sample_ID)
}


#' Update df after normalization.
#' 
#' @param df A data frame of peptide or protein table.
#' @param label_scheme_fit Experiment summary for samples being selected for
#'   fitting.
#' @param cf_x_fit A data frame containing the \code{x} positions for each
#'   samples indicated in label_scheme_fit.
#' @param sd_coefs_fit The standard deviations for each samples indicated in
#'   label_scheme_fit.
update_df <- function (df, label_scheme_fit, cf_x_fit, sd_coefs_fit) {
  nm_log2r_n <- names(df) %>% 
    .[grepl("^N_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
    find_fit_nms(label_scheme_fit$Sample_ID)
  
  nm_int_n <- names(df) %>% 
    .[grepl("^N_I[0-9]{3}[NC]*\\s+\\(", .)] %>% 
    find_fit_nms(label_scheme_fit$Sample_ID)
  
  nm_log2r_z <- names(df) %>% 
    .[grepl("^Z_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
    find_fit_nms(label_scheme_fit$Sample_ID)  
  
  df_z <- mapply(normSD, df[, nm_log2r_n, drop = FALSE], 
                 center = cf_x_fit$x, SD = sd_coefs_fit$fct, SIMPLIFY = FALSE) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
    `rownames<-`(rownames(df))    
  
  nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
  df_z[, nan_cols] <- 0
  rm(nan_cols)
  
  if (purrr::is_empty(nm_log2r_z)) {
    df <- cbind(df, df_z)
    
    nm_log2r_z <- names(df) %>% 
      .[grepl("^Z_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
      find_fit_nms(label_scheme_fit$Sample_ID)
  } else {
    df[, nm_log2r_z] <- df_z
  }
  
  # aligned log2FC and intensity after the calculation of "df_z"
  df[, nm_log2r_n] <- sweep(df[, nm_log2r_n, drop = FALSE], 2, cf_x_fit$x, "-")
  df[, nm_int_n] <- sweep(df[, nm_int_n, drop = FALSE], 2, 2^cf_x_fit$x, "/")    
  
  rlang::env_bind(caller_env(), nm_log2r_z = nm_log2r_z)
  
  return(df)
}


#'Data normalization
#'
#'\code{normMulGau} normalizes \code{log2FC} under the assumption of multi
#'Gaussian kernels.
#'
#'@param df An input data frame
#'@inheritParams prnHist
#'@inheritParams standPep
#'@return A data frame.
#'
#'@import dplyr purrr  
#'@importFrom magrittr %>% %T>% %$% %<>% 
normMulGau <- function(df, method_align, n_comp, seed = NULL, 
                       range_log2r, range_int, filepath, 
                       col_select = NULL, cut_points = Inf, ...) {

  # check n_comp
  find_n_comp <- function (n_comp, method_align) {
    if (is.null(n_comp)) {
      if (method_align == "MGKernel") {
        n_comp <- ifelse(nrow(df) > 3000, 3L, 2L)
        n_comp <- n_comp %>% as.integer()
      } else if (method_align == "MC") {
        n_comp <- 1L
      } else {
        n_comp <- 1L
      }
    } else {
      stopifnot(n_comp >= 1)
      if (method_align == "MGKernel") stopifnot(n_comp >= 2)
    }
    
    return(n_comp)
  }
  
  
  # find the x at max density
  my_which_max <- function (params, label_scheme) {
    fit <- params %>%
      split(., .$Sample_ID) %>%
      lapply(sumdnorm, xmin = -2, xmax = 2, by = 2/400) %>%
      do.call(rbind, .) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
      dplyr::arrange(Sample_ID)
    
    cf_x <- suppressWarnings(
      fit %>%
        dplyr::group_by(Sample_ID) %>% 
        dplyr::mutate(Max = max(Sum, na.rm = TRUE)) %>% 
        dplyr::mutate(Max = ifelse(is.infinite(Max), NA, Max)) %>% 
        dplyr::filter(Sum == Max) %>%
        dplyr::mutate(x = mean(x, na.rm = TRUE)) %>% # tie-breaking
        dplyr::filter(!duplicated(x)) %>% 
        dplyr::select(-Max)
    )
    
    cf_empty <- fit %>%
      dplyr::filter(! Sample_ID %in% cf_x$Sample_ID)
    
    if (nrow(cf_empty) > 0) {
      cf_empty <- cf_empty %>%
        dplyr::group_by(Sample_ID) %>%
        dplyr::filter(!duplicated(Sum)) %>%
        dplyr::mutate(x = 0)
    }
    
    # add back the empty samples
    cf_x <- dplyr::bind_rows(cf_x, cf_empty) %>% 
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
      dplyr::arrange(Sample_ID)
  }
  

  # compare to prior n_comp value
  ok_file_ncomp <- function(filepath, filename, n_comp) {
    if (file.exists(file.path(filepath, filename))) {
      params <- read.table(file.path(filepath, filename), 
                           check.names = FALSE, header = TRUE, 
                           comment.char = "#")
      
      n_comp == dplyr::n_distinct(params$Component)
    } else {
      return(FALSE)
    }
  }
  
  
  # calculates `x` positions by spline points (only for method_align = MC)
  spline_coefs <- function (df, label_scheme, label_scheme_fit, ...) {
    dots <- rlang::enexprs(...)
    
    slice_dots <- dots %>% 
      .[purrr::map_lgl(., is.language)] %>% 
      .[grepl("^slice_", names(.))]
    
    # initialization: NA for Empty samples; 0 for the remaining
    x_vals <- df %>%
      dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\(", names(.))) %>% 
      `colnames<-`(gsub("^N_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>%
      dplyr::summarise_all(~ median(.x, na.rm = TRUE)) %>%
      unlist() %>%
      data.frame(x = .) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::filter(.data$Sample_ID %in% label_scheme$Sample_ID) %>% 
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
      dplyr::arrange(Sample_ID) %>% 
      dplyr::mutate(x = ifelse(is.na(x), NA, 0))
    
    x_vals_fit <- df %>%
      filters_in_call(!!!slice_dots) %>% 
      dplyr::select(grep("^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\(", names(.))) %>% 
      `colnames<-`(gsub("^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\((.*)\\)$", "\\1", names(.))) %>%
      dplyr::select(which(names(.) %in% label_scheme_fit$Sample_ID)) %>% 
      dplyr::summarise_all(~ median(.x, na.rm = TRUE)) %>%
      unlist() %>%
      data.frame(x = .) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
      dplyr::arrange(Sample_ID)
    
    rows <- x_vals$Sample_ID %in% x_vals_fit$Sample_ID
    x_vals[rows, ] <- x_vals_fit
    
    invisible(x_vals)
  }
  
  
  # add mean-deviation info for `expt_smry::Select`ed samples
  add_mean_dev <- function (df, label_scheme_fit) {
    if (grepl("Peptide\\\\", filepath) || grepl("Peptide/", filepath)) {
      prefix <- "pep_mean_"
    } else if (grepl("Protein\\\\", filepath) || grepl("Protein/", filepath)) {
      prefix <- "prot_mean_"
    } else {
      return(df)
    }
    
    nm_log2r_raw <- names(df) %>% 
      .[grepl("^log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
      find_fit_nms(label_scheme_fit$Sample_ID)
    
    nm_log2r_n <- names(df) %>% 
      .[grepl("^N_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
      find_fit_nms(label_scheme_fit$Sample_ID)
    
    nm_log2r_z <- names(df) %>% 
      .[grepl("^Z_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
      find_fit_nms(label_scheme_fit$Sample_ID)
    
    df <- df %>% 
      dplyr::mutate(!!paste0(prefix, "raw") := rowMeans(.[, names(.) %in% nm_log2r_raw], 
                                                        na.rm = TRUE), 
                    !!paste0(prefix, "n") := rowMeans(.[, names(.) %in% nm_log2r_n], 
                                                      na.rm = TRUE), 
                    !!paste0(prefix, "z") := rowMeans(.[, names(.) %in% nm_log2r_z], 
                                                      na.rm = TRUE)) 
    
    df[, grepl(paste0("^", prefix), names(df))] <- 
      df[, grepl(paste0("^", prefix), names(df))] %>%
      dplyr::mutate_if(is.logical, as.numeric) %>%
      round(., digits = 3)

    if (prefix == "pep_mean_") {
      df <- dplyr::bind_cols(
        df %>% dplyr::select(grep("^prot_", names(.))), 
        df %>% dplyr::select(grep("^pep_", names(.))), 
        df %>% dplyr::select(-grep("^pep_|^prot_", names(.)))
      )
    } else {
      df <- local({
        new_nms <- paste0(prefix, c("raw", "n", "z"))
        
        nm_last <- names(df) %>% 
          .[! . %in% new_nms] %>% 
          .[grepl(paste0("^", prefix %>% gsub("mean_", "", .)), .)] %>% 
          .[length(.)]
        
        idx <- which(names(df) == nm_last)
        
        dplyr::bind_cols(
          df %>% dplyr::select(1:idx), 
          df %>% dplyr::select(new_nms), 
          df %>% dplyr::select((idx+1):ncol(.)) %>% dplyr::select(-new_nms)
        )    
      })
    }
  }
  
  
  plot_foo <- function () {
    # Pass arguments by row
    args <- params %>% 
      dplyr::select(c("lambda", "mean", "sd")) %>% pmap(list)
    
    # Define the function, i.e., normal density function with the weight of lambda
    wt_dnorm <- function(x, lambda, mean, sd) {
      lambda * length(x) * dnorm(x, mean, sd)
    }
    
    # The wrapper of stat_function(); easier to pass "x"
    stat_dnorm <- function(x, args) {
      stat_function(aes(x), fun = wt_dnorm, n = 100, args = args, size = .2)
    }
    
    # Pass the list of "args" to the wrapper function "stat_dnorm()" for plotting
    # ggplot() + lapply(args, stat_dnorm, x = seq(-2, 2, 0.1))
  }
  
  
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  dat_dir <- get_gl_dat_dir()

	dots <- rlang::enexprs(...)
	
	slice_dots <- dots %>% 
	  .[purrr::map_lgl(., is.language)] %>% 
	  .[grepl("^slice_", names(.))]
	
	nonslice_dots <- dots %>% 
	  .[! . %in% slice_dots]
	
	n_comp <- find_n_comp(n_comp, method_align)

	if (!purrr::is_empty(nonslice_dots)) {
	  data.frame(nonslice_dots) %>%
			dplyr::bind_cols(n_comp = n_comp) %>%
			write.table(file = file.path(filepath, "normalmixEM_pars.txt"),
			            sep = "\t", col.names = TRUE, row.names = FALSE)
	}

	# if different `n_comp` between two `method_align = MGKernel`, 
	#   force `col_select` to all samples
	# if `n_comp` is given but with `method_align = MC`, 
	#   ignore difference in n_comp
	
	ok_N_ncomp <- ok_file_ncomp(filepath, "MGKernel_params_N.txt", n_comp)
	ok_Z_ncomp <- ok_file_ncomp(filepath, "MGKernel_params_Z.txt", n_comp)
	
	if (method_align == "MGKernel") {
  	if ((!ok_N_ncomp) && (col_select != rlang::expr(sample_ID))) {
  	  col_select <- rlang::expr(Sample_ID)
  	}
  
  	if ((!ok_Z_ncomp) && (col_select != rlang::expr(sample_ID))) {
  	  col_select <- rlang::expr(Sample_ID)
  	}
	}

	load(file = file.path(dat_dir, "label_scheme.rda"))
	label_scheme_fit <- label_scheme %>% .[!is.na(.[[col_select]]), ]
	
	## `nm_log2r_z` the same as `nm_log2r_n` 
	#   if is the first `normMulGau` after `mergePep` or `mergePrn`
	
	nm_log2r_n <- names(df) %>% 
	  .[grepl("^N_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
	  find_fit_nms(label_scheme_fit$Sample_ID)
	
	nm_int_n <- names(df) %>% 
	  .[grepl("^N_I[0-9]{3}[NC]*\\s+\\(", .)] %>% 
	  find_fit_nms(label_scheme_fit$Sample_ID)

	if (method_align == "MGKernel") {
	  message("method_align = ", method_align)
	  message("n_comp = ", n_comp)
	  message("col_select = ", col_select)

	  params_sub <- df %>% 
	    filters_in_call(!!!slice_dots) %>% 
	    dplyr::select(nm_log2r_n) %>% 
	    `names<-`(gsub("^N_log2_R[0-9]{3}[NC]{0,1}\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
	    fitKernelDensity(n_comp = n_comp, seed = seed, !!!nonslice_dots) %>% 
	    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
	    dplyr::arrange(Sample_ID, Component)

    if (!ok_N_ncomp) {
      # earlierly forced `col_select = Sample_ID` if detected different `n_comp` 
      # so if not `ok`, `params_sub` must be for all samples
      params <- params_sub
    } else {
      params <- read.table(file.path(filepath, "MGKernel_params_N.txt"), 
                           check.names = FALSE, header = TRUE, comment.char = "#") %>% 
	      dplyr::select(names(params_sub)) %>% 
	      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
	      dplyr::arrange(Sample_ID, Component)
      
	    rows <- params$Sample_ID %in% params_sub$Sample_ID
	    params[rows, ] <- params_sub
    }
    
	  # profile widths based on all sample columns and data rows
	  sd_coefs <- calc_sd_fcts(df, range_log2r, range_int, label_scheme)

		cf_x <- my_which_max(params, label_scheme) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) 

		list(params, cf_x, sd_coefs) %>%
			purrr::reduce(dplyr::left_join, by = "Sample_ID") %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) %>%
			write.table(., file = file.path(filepath, "MGKernel_params_N.txt"), sep = "\t",
			            col.names = TRUE, row.names = FALSE)
		
		sd_coefs_fit <- sd_coefs %>% 
		  dplyr::filter(Sample_ID %in% label_scheme_fit$Sample_ID)
		
		cf_x_fit <- cf_x %>% 
		  dplyr::filter(Sample_ID %in% label_scheme_fit$Sample_ID)
		
		df <- update_df(df, label_scheme_fit, cf_x_fit, sd_coefs_fit)
		
		## non-empty `Z_log2_R` guaranteed after `update_df`
		
		# separate fits of Z_log2_R for updating curve parameters only
		if (!ok_Z_ncomp) {
		  params_z <- df %>% 
		    filters_in_call(!!!slice_dots) %>% 
		    dplyr::select(nm_log2r_z) %>% 
		    `names<-`(gsub("^Z_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
		    fitKernelDensity(n_comp = n_comp, seed = seed, !!!nonslice_dots) %>% 
		    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
		    dplyr::arrange(Sample_ID, Component) %>% 
		    dplyr::mutate(x = 0)
		} else {
		  params_z_sub <- df %>% 
		    filters_in_call(!!!slice_dots) %>% 
		    dplyr::select(nm_log2r_z) %>% 
		    `names<-`(gsub("^Z_log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
		    fitKernelDensity(n_comp = n_comp, seed, !!!nonslice_dots) %>% 
		    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
		    dplyr::arrange(Sample_ID, Component)
		  
		  params_z <- read.table(file.path(filepath, "MGKernel_params_Z.txt"), 
		                         check.names = FALSE, header = TRUE, comment.char = "#") %>% 
		    dplyr::select(names(params_z_sub)) %>% 
		    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
		    dplyr::arrange(Sample_ID, Component)
		  
		  rows_z <- params_z$Sample_ID %in% params_z_sub$Sample_ID
		  params_z[rows_z, ] <- params_z_sub
		  
		  params_z$x <- 0
		}	
		
		write.table(params_z, file = file.path(filepath, "MGKernel_params_Z.txt"),
		            sep = "\t", col.names = TRUE, row.names = FALSE)

	} else if (method_align == "MC") {
	  message("method_align = ", method_align)
	  message("n_comp = NULL")
	  message("col_select = ", col_select)
	  
	  # profile widths based on all sample columns and data rows
	  sd_coefs <- df %>% calc_sd_fcts(range_log2r, range_int, label_scheme)
	  
	  # cut_points = c(mean_lint = seq(4, 7, .5))
	  # cut_points = c(prot_icover = seq(.25, .75, .25))
	  # cut_points = c(prot_icover = Inf)
	  # cut_points = c(prot_icover = NULL)
	  # cut_points = c(prot_icover = NA)
	  # cut_points = Inf
	  # cut_points = NULL
	  # cut_points = NA
	  
	  cut_points <- set_cutpoints2(cut_points, df)

	  if (all(is.infinite(cut_points))) {
	    df <- list(df)
	  } else {
	    df <- df %>% 
	      dplyr::mutate(col_cut = !!rlang::sym(names(cut_points)[1])) %>% 
	      dplyr::mutate_at(.vars = "col_cut", cut, 
	                       breaks = cut_points, 
	                       labels = cut_points[1:(length(cut_points)-1)]) %>%
	      split(., .$col_cut, drop = TRUE) %>% 
	      purrr::map(~ .x %>% dplyr::select(-col_cut))
	  }

	  x_vals <- df %>% 
	    purrr::map(spline_coefs, label_scheme, label_scheme_fit, !!!slice_dots)
	  
	  df <- purrr::map2(df, x_vals, ~ update_df(.x, label_scheme, .y, sd_coefs)) %>% 
	    dplyr::bind_rows()
	}
  
	# (1) mean deviation based sample IDs in `label_scheme` not `label_scheme_fit`
	#  as sample IDs in `label_scheme_fit` may be only for mixed-bed normalization
	# (2) the `mean` also get updated after normalization against sample subsets 
	label_scheme_not_trivial <- label_scheme %>% 
	  dplyr::filter(!Reference, !grepl("^Empty\\.[0-9]+", Sample_ID))
	
	df <- df %>% add_mean_dev(label_scheme_not_trivial)
	
	return(df)
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
dblTrim <- function(df, range_log2r, range_int, type_r = "N_log2_R", type_int = "N_I") {
	df_trim <- df
	
	type_r <- paste0("^", type_r, "[0-9]{3}")
	type_int <- paste0("^", type_int, "[0-9]{3}")
	
	# trim by log2-ratios
	col_r <- grepl(type_r, names(df_trim))
	df_trim[, col_r] <- lapply(df_trim[, col_r], function (x) {
			q_ratio <- quantile(x, probs = range_log2r/100, na.rm = TRUE)
			x[x < q_ratio[1] | x > q_ratio[2]] <- NA
			return(x)
		}
	)

	# trim by intensity
	col_int <- grepl(type_int, names(df_trim))
	df_trim[, col_int] <- lapply(df_trim[, col_int], function (x) {
			q_intensity <- quantile(x, probs = range_int/100, na.rm = TRUE)
			x[x < q_intensity[1] | x > q_intensity[2]] <- NA
			return(x)
		}
	)

	# doubly trim
	df_trim[!is.na(df_trim)] <- 1

	df_trim <- mapply(`*`, df_trim[, grepl(type_r, names(df_trim))],
	                  df_trim[, grepl(type_int, names(df_trim))], SIMPLIFY = FALSE) %>%
		data.frame(check.names = FALSE)

	df_trim[] <- mapply(`*`, df[, grepl(type_r, names(df))] , df_trim, SIMPLIFY = FALSE)

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
sumdnorm <- function (x, xmin = -4, xmax = 4, by = xmax/200) {
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
normSD <- function (x, center = 0, SD = 1) {
	if (sum(is.na(x)) == length(x)) {
		x
	} else if ((sum(is.na(x)) + sum(x == 0, na.rm = TRUE)) == length(x)) {
		x[1:length(x)] <- NaN
	} else if (all(x == 0)) {
		x
	} else {
		x <- (x - center) / SD
	}
	return (x)
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
fitKernelDensity <- function (df, n_comp = 3, seed = NULL, ...) {

	nmix_params <- function (x, n_comp = 3, seed = seed, ...) {
		dots <- rlang::enexprs(...)
		x <- rlang::enexpr(x)
		
		if (!is.null(dots$k)) {
			cat(paste("k =", dots$k, "replaced by", paste("n_comp =", n_comp, "\n")))
			dots$k <- NULL
		}

		if (sum(is.na(x)) == length(x) |
					(sum(is.na(x)) + sum(x == 0, na.rm = TRUE)) == length(x) |
					all(x == 0)) {
			df_par <- data.frame(Component = c(1:n_comp), lambda = rep(NA, n_comp),
									mean = rep(NA, n_comp), sd = rep(NA, n_comp))
		} else {
			x <- x[!is.na(x)]

			stopifnot(n_comp > 1)

			mixEM_call <- rlang::quo(mixtools::normalmixEM(!!x, k = !!n_comp, !!!dots))
			if (!is.null(seed)) set.seed(seed) else set.seed(sample(.Random.seed, 1))
			
			quietly_out <- purrr::quietly(rlang::eval_tidy)(mixEM_call, caller_env())
			x_k2 <- quietly_out$result
			df_par <- data.frame(Component = 1:n_comp, 
			                     lambda = x_k2$lambda, 
			                     mean = x_k2$mu, 
			                     sd = x_k2$sigma)
		}

		return(df_par)
	}

	dots <- rlang::enexprs(...)

	min_n <- local({
	  ok_nan_cols <- df %>% 
      purrr::map_lgl(not_all_nan, na.rm = TRUE)
	  
	  min_n <- df %>% 
      .[, ok_nan_cols, drop = FALSE] %>% 
      .[, not_all_NA(.), drop = FALSE] %>% 
      purrr::map_dbl(., ~ (!is.na(.x)) %>% sum()) %>% 
      min()
	})

	if (min_n < 50) {
	  stop("Too few data points for fitting with multiple Gaussian functions.", 
	       call. = FALSE)
	}

	lapply(df, nmix_params, n_comp, seed = seed, !!!dots) %>%
		do.call(rbind, .) %>%
		dplyr::mutate(Sample_ID = rownames(.)) %>%
		dplyr::mutate(Sample_ID = gsub("(.*)\\.\\d+$", "\\1", Sample_ID)) %>%
		dplyr::mutate(Channel = rep(1:(nrow(.)/n_comp), each = n_comp)) %>%
		dplyr::arrange(Channel, -lambda) %>%
		dplyr::mutate(Component = rep(1:n_comp, nrow(.)/n_comp)) %>%
		dplyr::mutate(Height = .$lambda * dnorm(.$mean, mean = .$mean, sd = .$sd))
}

