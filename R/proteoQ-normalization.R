#'Data normalization
#'
#'\code{normMulGau} normalizes \code{log2-ratios} under the assumption of multi
#'Gaussian kernels.
#'
#' @param df An input data frame
#'@inheritParams mixtools::normalmixEM
#'@inheritParams normPep
#'@return A data frame.
#'
#'@import dplyr purrr rlang mixtools
#'@importFrom magrittr %>%
normMulGau <- function(df, method_align, n_comp, seed = NULL, range_log2r, range_int, filepath, 
                       col_refit = NULL, ...) {

  my_which_max <- function (params) {
    fit <- params %>%
      split(., .$Sample_ID) %>%
      lapply(sumdnorm, xmin = -2, xmax = 2, by = 2/400) %>%
      do.call(rbind, .) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
      dplyr::arrange(Sample_ID)
    
    ## Calibration coefficients for centering ratio profiles
    cf_x <- fit %>%
      dplyr::group_by(Sample_ID) %>%
      dplyr::filter(Sum == max(Sum, na.rm = TRUE)) %>%
      dplyr::mutate(x = mean(x, na.rm = TRUE)) %>% # tie-breaking
      dplyr::filter(!duplicated(x))
    
    cf_empty <- fit %>%
      dplyr::filter(! Sample_ID %in% cf_x$Sample_ID)
    
    if(nrow(cf_empty) > 0) {
      cf_empty <- cf_empty %>%
        dplyr::group_by(Sample_ID) %>%
        dplyr::filter(!duplicated(Sum)) %>%
        dplyr::mutate(x = 0)
    }
    
    cf_x <- dplyr::bind_rows(cf_x, cf_empty) %>% # add back the empty samples
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% # avoid mismatch
      dplyr::arrange(Sample_ID)
  }
 
  calc_sd_fcts <- function (df, range_log2r, range_int, label_scheme) {
    label_scheme_sd <- label_scheme %>%
      dplyr::filter(!Reference, !grepl("^Empty\\.", Sample_ID))	%>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = (.$Sample_ID)))
    
    SD <- df %>%
      dplyr::select(grep("^N_log2_R|^N_I", names(.))) %>%
      dblTrim(., range_log2r, range_int) %>%
      `names<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.)))
    
    cf_SD <- SD/mean(SD %>% .[names(.) %in% label_scheme_sd$Sample_ID], na.rm = TRUE)
    cf_SD <- cbind.data.frame(fct = cf_SD, SD) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
      dplyr::arrange(Sample_ID)
  }
  
  find_fit_nms <- function(nm_a, nm_b) {
    ind <- purrr::map(nm_b, ~ grepl(.x, nm_a)) %>% 
      purrr::reduce(`|`)
    nm_a <- nm_a[ind]
  }
  
	
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)

	dots <- rlang::enexprs(...)

	if(!purrr::is_empty(dots)) {
		data.frame(dots) %>%
			bind_cols(n_comp = n_comp) %>%
			write.table(., file = file.path(filepath, "normalmixEM_pars.txt"),
			            sep = "\t", col.names = TRUE, row.names = FALSE)
	}

	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	  
	label_scheme_fit <- label_scheme %>% 
	  .[!is.na(.[[col_refit]]), ]
	
	nm_log2r_n <- names(df) %>% 
	  .[grepl("^N_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
	  find_fit_nms(label_scheme_fit$Sample_ID)

	nm_int_n <- names(df) %>% 
	  .[grepl("^N_I[0-9]{3}[NC]*\\s+\\(", .)] %>% 
	  find_fit_nms(label_scheme_fit$Sample_ID)
	
	if (method_align == "MGKernel") {
		print(paste("Number of Gaussian components =", n_comp))
	  
	  # N_log2_R... in the inputs are median centered at the first pass
	  params_sub <- df[, nm_log2r_n, drop = FALSE] %>% 
	    `names<-`(gsub("^N_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>% 
	    fitKernelDensity(n_comp = n_comp, seed = seed, !!!dots) %>% 
	    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
	    dplyr::arrange(Sample_ID, Component)

    if (!file.exists(file.path(filepath, "MGKernel_params_N.txt"))) {
      # warning("First-pass normalization for the complete data set.")
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

	  # standard deviation
	  cf_SD <- calc_sd_fcts(df, range_log2r, range_int, label_scheme)

	  # x values at max density
		cf_x <- my_which_max(params) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) 

		list(params, cf_x, cf_SD) %>%
			purrr::reduce(left_join, by = "Sample_ID") %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID) %>%
			write.table(., file = file.path(filepath, "MGKernel_params_N.txt"), sep = "\t",
			            col.names = TRUE, row.names = FALSE)
		

		# ========================================================================================
		# not run
		#
		# Pass arguments by row
		args <- params %>%
						dplyr::select(c("lambda", "mean", "sd")) %>% pmap(list)

		# Define the function, i.e., normal density function with the weight of lambda
		wt_dnorm <- function(x, lambda, mean, sd) lambda * length(x) * dnorm(x, mean, sd)

		# The wrapper of stat_function(); easier to pass "x"
		stat_dnorm <- function(x, args) stat_function(aes(x), fun = wt_dnorm, n = 100, args = args, size = .2)

		# Pass the list of "args" to the wrapper function "stat_dnorm()" for plotting
		# ggplot() + lapply(args, stat_dnorm, x = seq(-2, 2, 0.1))
		#
		# ========================================================================================
		
		cf_SD_fit <- cf_SD %>% 
		  dplyr::filter(Sample_ID %in% label_scheme_fit[[col_refit]])
		
		cf_x_fit <- cf_x %>% 
		  dplyr::filter(Sample_ID %in% label_scheme_fit[[col_refit]])

		# data with MGKernel alignment and scaling normalization
		df_z <- mapply(normSD, df[, nm_log2r_n, drop = FALSE], 
		               center = cf_x_fit$x, SD = cf_SD_fit$fct, SIMPLIFY = FALSE) %>%
		  data.frame(check.names = FALSE) %>%
		  `colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
		  `rownames<-`(rownames(df))
		
		nm_log2r_z <- names(df) %>% 
		  .[grepl("^Z_log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
		  find_fit_nms(label_scheme_fit$Sample_ID)
		
		if (rlang::is_empty(nm_log2r_z)) {
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

		# separate fits of Z_log2_r for curve parameters only...
		if (!file.exists(file.path(filepath, "MGKernel_params_Z.txt"))) {
		  warning("First-pass normalization for the complete data set.")
		  
		  params_z <- df[, nm_log2r_z] %>%
		    `names<-`(gsub("^Z_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>% 
		    fitKernelDensity(n_comp = n_comp, seed, !!!dots) %>% 
		    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>% 
		    dplyr::arrange(Sample_ID, Component) %>% 
		    dplyr::mutate(x = 0)
		} else {
		  params_z_sub <- df[, nm_log2r_z, drop = FALSE] %>% 
		    `names<-`(gsub("^Z_log2_R[0-9]{3}.*\\((.*)\\)$", "\\1", names(.))) %>% 
		    fitKernelDensity(n_comp = n_comp, seed, !!!dots) %>% 
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
		cf_x <- df %>%
			dplyr::select(matches("^N_log2_R[0-9]{3}")) %>%
			`colnames<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.))) %>%
			dplyr::summarise_all(funs(median(., na.rm = TRUE))) %>%
			unlist() %>%
			data.frame(x = .) %>%
			tibble::rownames_to_column("Sample_ID") %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID)

		df <- mapply(normSD, df[,grepl("^N_log2_R[0-9]{3}", names(df))],
		             center = cf_x$x, SD = cf_SD$fct, SIMPLIFY = FALSE) %>% # scale to the same SD
				data.frame(check.names = FALSE) %>%
				`colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
				cbind(df, .)

		df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
		  sweep(df[, grepl("^N_log2_R[0-9]{3}", names(df))], 2, cf_x$x, "-")
		df[, grepl("^N_I[0-9]{3}", names(df))] <-
		  sweep(df[, grepl("^N_I[0-9]{3}", names(df))], 2, 2^cf_x$x, "/")

	} else { # method_align <- c("YWHAE", "CBR4")
		if(!all(method_align %in% df[["gene"]])) {
			stop(paste(setdiff(method_align, df[["gene"]]),
			           "in 'method_align' not found in gene names."), call. = TRUE)
		}

		cf_x <- df %>%
			dplyr::filter(.[["gene"]] %in% method_align) %>%
			dplyr::select(matches("^N_log2_R[0-9]{3}")) %>%
			`colnames<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.))) %>%
			dplyr::summarise_all(funs(median(., na.rm = TRUE))) %>%
			unlist() %>%
			data.frame(x = .) %>%
			tibble::rownames_to_column("Sample_ID") %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = label_scheme$Sample_ID)) %>%
			dplyr::arrange(Sample_ID)

		if (sum(!is.na(cf_x$x)) == 0) {
			stop(paste("Protein ID/Alignment method", method_align, "not found.") , call. = TRUE)
		} else {
			df <- mapply(normSD, df[,grepl("^N_log2_R[0-9]{3}", names(df))],
			             center = cf_x$x, SD = cf_SD$fct, SIMPLIFY = FALSE) %>% # scale to the same SD
					data.frame(check.names = FALSE) %>%
					`colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
					cbind(df, .)

			df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
			  sweep(df[, grepl("^N_log2_R[0-9]{3}", names(df))], 2, cf_x$x, "-")
			df[, grepl("^N_I[0-9]{3}", names(df))] <-
			  sweep(df[, grepl("^N_I[0-9]{3}", names(df))], 2, 2^cf_x$x, "/")
		}

	}

	return(df)
}


#' Data normalization
#'
#' \code{dblTrim} doubly trims the \code{log2FC} and reporter-ion intensity by
#' the given ranges.
#'
#' @inheritParams normPep
#' @return A data frame.
#'
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
dblTrim <- function(df, range_log2r, range_int) {
	df_trim <- df

	# trim by log2-ratios
	col_r <- grepl("^N_log2_R", names(df_trim))
	df_trim[, col_r] <- lapply(df_trim[, col_r], function (x) {
			q_ratio <- quantile(x, probs = range_log2r/100, na.rm = TRUE)
			x[x < q_ratio[1] | x > q_ratio[2]] <- NA
			return(x)
		}
	)

	# trim by intensity
	col_int <- grepl("^N_I", names(df_trim))
	df_trim[, col_int] <- lapply(df_trim[, col_int], function (x) {
			q_intensity <- quantile(x, probs = range_int/100, na.rm = TRUE)
			x[x < q_intensity[1] | x > q_intensity[2]] <- NA
			return(x)
		}
	)

	# doubly trim
	df_trim[!is.na(df_trim)] <- 1  # boolean matrix

	df_trim <- mapply(`*`, df_trim[, grepl("^N_log2_R", names(df_trim))],
	                  df_trim[, grepl("^N_I", names(df_trim))], SIMPLIFY = FALSE) %>%
		data.frame(check.names = FALSE) # doubly trimmed boolean matrix

	df_trim[] <- mapply(`*`, df[, grepl("^N_log2_R", names(df))] , df_trim, SIMPLIFY = FALSE)

	sapply(df_trim, sd, na.rm=T)
}


#' Data normalization
#'
#' \code{sumdnorm} calculates summed density from \code{normMulGau}.
#'
#' @return A data frame.
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
sumdnorm <- function (x, xmin = -4, xmax = 4, by = xmax/200) {
	wt_dnorm <- function (x, lambda, mean, sd) lambda * dnorm(x, mean = mean, sd = sd)

	args <- purrr::pmap(x[, names(x) %in% c("lambda", "mean", "sd")], list) %>%
		`names<-`(x$Sample_ID)

	# nm_comps <- paste("G", seq_len(length(args)), sep = ".")
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
#' \code{normSD} normalizes the SD of \code{log2-ratios}.
#'
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
#' @export
normSD <- function (x, center = 0, SD = 1) {
	if (sum(is.na(x)) == length(x)) {
		x
	} else if ((sum(is.na(x))+sum(x==0, na.rm=T)) == length(x)) {
		x[1:length(x)] <- NaN
	} else if (all(x == 0)) {
		x
	} else {
		x <- (x-center)/SD
	}
	return (x)
}


#' Data normalization
#'
#' \code{fitKernelDensity} calculates the fitted density of \code{log2-ratios}.
#'
#' @return A data frame.
#'
#' @import dplyr purrr rlang mixtools
#' @importFrom magrittr %>%
fitKernelDensity <- function (df, n_comp = 3, seed = NULL, ...) {

	nmix_params <- function (x, n_comp = 3, seed = seed, ...) {
		dots <- rlang::enexprs(...)
		x <- rlang::enexpr(x)
		
		if(!is.null(dots$k)) {
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
			x_k2 <- rlang::eval_tidy(mixEM_call, caller_env())
			df_par <- data.frame(Component = 1:n_comp, lambda = x_k2$lambda, mean = x_k2$mu, sd = x_k2$sigma)
		}

		return(df_par)
	}

	dots <- rlang::enexprs(...)

	lapply(df, nmix_params, n_comp, seed = seed, !!!dots) %>%
		do.call(rbind, .) %>%
		dplyr::mutate(Sample_ID = rownames(.)) %>%
		dplyr::mutate(Sample_ID = gsub("(.*)\\.\\d+$", "\\1", Sample_ID)) %>%
		dplyr::mutate(Channel = rep(1:(nrow(.)/n_comp), each = n_comp)) %>%
		dplyr::arrange(Channel, -lambda) %>%
		dplyr::mutate(Component = rep(1:n_comp, nrow(.)/n_comp)) %>%
		dplyr::mutate(Height = .$lambda * dnorm(.$mean, mean = .$mean, sd = .$sd))
}

