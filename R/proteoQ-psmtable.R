theme_psm_violin <- theme_bw() +
	theme(
	axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 24),
	axis.ticks.x  = element_blank(),
	axis.text.y  = element_text(angle = 0, vjust = 0.5, size = 24),
	axis.title.x = element_text(colour = "black", size = 24),
	axis.title.y = element_text(colour = "black", size = 24),
	plot.title = element_text(colour = "black", size  =24, hjust = .5, vjust = .5),

	strip.text.x = element_text(size = 18, colour = "black", angle = 0),
	strip.text.y = element_text(size = 18, colour = "black", angle = 90),

	panel.grid.major.x = element_blank(),
	panel.grid.minor.x = element_blank(),
	panel.grid.major.y = element_blank(),
	panel.grid.minor.y = element_blank(),

	legend.key = element_rect(colour = NA, fill = 'transparent'),
	legend.background = element_rect(colour = NA,  fill = "transparent"),
	legend.position = "none",
	legend.title = element_blank(),
	legend.text = element_text(colour = "black", size = 18),
	legend.text.align = 0,
	legend.box = NULL
)


#' Removes PSM headers
#'
#' \code{rmPSMHeaders} removes the header of PSM from
#' \code{\href{https://http://www.matrixscience.com/}{Mascot}} outputs. It also
#' removes the spacer columns in the fields of ratio and intensity values.
#'
#' @return Header file(s) and PSM table(s) without the header.
#'
#' @import dplyr
#' @importFrom purrr walk
#' @importFrom magrittr %>%
#' @export
rmPSMHeaders <- function () {
	dir.create(file.path(dat_dir, "PSM\\Violin"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\cache"), recursive = TRUE, showWarnings = FALSE)

	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	on.exit(message("Remove PSM headers --- Completed."), add = TRUE)

	options(max.print = 5000000)

	filelist = list.files(path = file.path(dat_dir), pattern = "^F[0-9]{6}\\.csv$")

	if(purrr::is_empty(filelist))
	  stop("Can't find PSM files(s) with `.csv` extension under ", dat_dir, call. = FALSE)

	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	TMT_plex <- TMT_plex(label_scheme)

	batchPSMheader <- function(filelist, TMT_plex) {
		temp <- readLines(file.path(dat_dir, filelist))

		first_nonheader_row <- grep("pep_seq", temp)
		temp_header <- temp[1 : (first_nonheader_row - 1)]
		temp_header <- gsub("\"", "", temp_header, fixed = TRUE)

		output_prefix <- gsub(".csv$", "", filelist)
		write.table(temp_header, file.path(dat_dir, "PSM\\cache",
		            paste0(output_prefix, "_header", ".txt")),
		            sep = "\t", col.names = FALSE, row.names = FALSE)
		rm(temp_header)

		# the remaining without header
		temp <- temp[first_nonheader_row : length(temp)]
		temp <- gsub("\"---\"", -1, temp,fixed = TRUE)
		temp <- gsub("\"###\"", -1, temp,fixed = TRUE)

		if(TMT_plex == 11) {
			temp[1] <- paste0(temp[1], paste(rep(",", 42), collapse = ''))
		} else if (TMT_plex == 10) {
			temp[1] <- paste0(temp[1], paste(rep(",", 38), collapse = ''))
		} else if(TMT_plex == 6) {
			temp[1] <- paste0(temp[1], paste(rep(",", 22), collapse = ''))
		}

		writeLines(temp, file.path(dat_dir, "PSM\\cache",
		                           paste0(output_prefix, "_hdr_rm.csv")))
		rm(temp)
	}

	purrr::walk(filelist, batchPSMheader, TMT_plex)
}


#' Extract RAW file names
#'
#' \code{extractRAW} extracts the RAW file names from the PSM result files under
#' the current working directory.
#'
#' @examples
#' extractRAW()
#'
#' @import dplyr tidyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
extractRAW <- function(rptr_intco = 1000, rm_craps = FALSE, plot_violins = TRUE) {
  dir.create(file.path(dat_dir, "PSM\\temp"), recursive = TRUE, showWarnings = FALSE)
  filelist = list.files(path = file.path(dat_dir), pattern = "^F[0-9]{6}\\.csv$")

  if(purrr::is_empty(filelist))
    stop("Can't find PSM files(s) with `.csv` extension under ", dat_dir, call. = FALSE)
  
  output_prefix <- gsub(".csv$", "", filelist)
  
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  TMT_plex <- TMT_plex(label_scheme)
  
  batchPSMheader_2 <- function(filelist, TMT_plex) {
    df <- readLines(file.path(dat_dir, filelist))
    
    first_nonheader_row <- grep("pep_seq", df)
    
    header <- df[1 : (first_nonheader_row - 1)]
    header <- gsub("\"", "", header, fixed = TRUE)
    
    output_prefix <- gsub(".csv$", "", filelist)
    write.table(header, file.path(dat_dir, "PSM\\temp",
                                  paste0(output_prefix, "_header", ".txt")), 
                sep = "\t", col.names = FALSE, row.names = FALSE)

    df <- df[first_nonheader_row : length(df)]
    df <- gsub("\"---\"", -1, df, fixed = TRUE)
    df <- gsub("\"###\"", -1, df, fixed = TRUE)
    
    if(TMT_plex == 11) {
      df[1] <- paste0(df[1], paste(rep(",", 42), collapse = ''))
    } else if (TMT_plex == 10) {
      df[1] <- paste0(df[1], paste(rep(",", 38), collapse = ''))
    } else if(TMT_plex == 6) {
      df[1] <- paste0(df[1], paste(rep(",", 22), collapse = ''))
    }
    
    writeLines(df, file.path(dat_dir, "PSM\\temp", paste0(output_prefix, "_hdr_rm.csv")))
  }
  
  purrr::walk(filelist, batchPSMheader_2, TMT_plex)

  df <- do.call(rbind,
                lapply(
                  list.files(path = file.path(dat_dir, "PSM\\temp"),
                             pattern = paste0("F[0-9]{6}", "_hdr_rm.csv"),
                             full.names = TRUE), read.delim, sep = ',',
                  check.names = FALSE, header = TRUE,
                  stringsAsFactors = FALSE, quote = "\"",
                  fill = TRUE , skip = 0
                )
  )
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  TMT_plex <- TMT_plex(label_scheme_full)
  
  # remove the spacer columns in Mascot outputs
  r_start <- which(names(df) == "pep_scan_title") + 1
  int_end <- ncol(df)
  if(int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
  
  if (TMT_plex == 11) {
    col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                   "R130N", "R130C", "R131N", "R131C")
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                   "R130N", "R130C", "R131")
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_ratio <- c("R127", "R128", "R129", "R130", "R131")
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_ratio <- NULL
    col_int <- NULL
  }
  
  df <- df[, "pep_scan_title", drop = FALSE] %>% 
    dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title))
  
  missing_msfs <- unique(df$RAW_File)
  
  if(!purrr::is_empty(missing_msfs)) cat(paste0(missing_msfs, "\n"))
  
  unlink(file.path(dat_dir, "PSM\\temp"), recursive = TRUE, force = TRUE)
}


#' Splits PSM tables
#'
#' \code{splitPSM} splits the PSM outputs after \code{rmPSMHeaders()}. It
#' separates PSM data by TMT experiment and LC/MS injection.
#'
#' @param rptr_intco Numeric; the threshold of reporter ion intensity.
#' @param rm_craps Logical; if TRUE, removes
#'   \code{\href{https://www.thegpm.org/crap/}{cRAP}} proteins.
#' @param plot_violins Logical; if TRUE, prepares the violin plots of
#'   reporter-ion intensities.
#' @examples
#' splitPSM(rptr_intco = 1000, rm_craps = TRUE, plot_violins = TRUE)
#'
#' @import dplyr tidyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
splitPSM <- function(rptr_intco = 1000, rm_craps = FALSE, plot_violins = TRUE, ...) {

	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."),
	        add = TRUE)

	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	load(file = file.path(dat_dir, "fraction_scheme.Rdata"))

	TMT_plex <- TMT_plex(label_scheme_full)

  filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
                        pattern = "^F[0-9]{6}\\_hdr_rm.csv$")

	if (length(filelist) == 0)
	  stop(paste("No PSM files were found under", file.path(dat_dir, "PSM")))

  df <- do.call(rbind,
		lapply(
			list.files(path = file.path(dat_dir, "PSM\\cache"),
			pattern = paste0("F[0-9]{6}", "_hdr_rm.csv"),
			full.names = TRUE), read.delim, sep = ',',
			check.names = FALSE, header = TRUE,
			stringsAsFactors = FALSE, quote = "\"",
			fill = TRUE , skip = 0
		)
	)
  
  df <- filters_in_call(df, ...)

  # remove the spacer columns in Mascot outputs
  r_start <- which(names(df) == "pep_scan_title") + 1
  int_end <- ncol(df)
	if(int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]

  if (TMT_plex == 11) {
		col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
		               "R130N", "R130C", "R131N", "R131C")
		col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
		             "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
		col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
		               "R130N", "R130C", "R131")
		col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
		             "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
		col_ratio <- c("R127", "R128", "R129", "R130", "R131")
		col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
		col_ratio <- NULL
		col_int <- NULL
	}

	if(TMT_plex > 0) {
		colnames(df)[r_start:(r_start+TMT_plex-2)] <- col_ratio
		colnames(df)[(r_start+TMT_plex-1):(r_start+TMT_plex+TMT_plex-2)] <- col_int
		rm(r_start, int_end, col_ratio, col_int)
	}

  if(rm_craps) df <- df %>% dplyr::filter(!grepl("\\|.*\\|$", prot_acc))

  df <- df %>% dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc))

  prn_acc <- parse_acc(df)
  label_scheme_full$Accession_Type <- prn_acc
  label_scheme$Accession_Type <- prn_acc
  save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.Rdata"))
  save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))

  if(length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ifelse(. == -1, NA, .)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ifelse(. <= rptr_intco, NA, .)) %>%
      dplyr::filter(pep_expect <=  0.10) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>%
      dplyr::mutate(RAW_File = gsub("^.*\\s{1}File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>%
      # a special case of redundant entries from Mascot
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::filter(pep_expect <=  0.10) %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>% 
			dplyr::mutate(RAW_File = gsub("^.*\\s{1}File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` under the experimental summary file."),
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")

	for (i in seq_along(df_split)) {
		df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)

		out_fn <- fn_lookup %>%
			dplyr::filter(TMT_inj == names(df_split)[i]) %>%
			dplyr::select(filename) %>%
			unique() %>%
			unlist() %>%
			paste0(., ".csv")

		df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)

		write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)

		if (plot_violins & TMT_plex > 0) {
			dir.create(file.path(dat_dir, "PSM\\Violin\\bf_olm"), recursive = TRUE, showWarnings = FALSE)

			TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()
			Levels <- TMT_levels %>% gsub("^TMT-", "I", .)
			df_int <- df_split[[i]][, grepl("^I[0-9]{3}", names(df_split[[i]]))] %>%
				tidyr::gather(key = "Channel", value = "Intensity") %>%
				dplyr::mutate(Channel = factor(Channel, levels = Levels))
			rm(Levels)
			
			mean_int <- df_int %>% 
			  dplyr::group_by(Channel) %>% 
			  dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
			  dplyr::mutate(Intensity = round(Intensity, digit = 1))
			
			Width <- 8
			Height <- 8

			filename <- gsub("\\.csv", "\\.png", out_fn)
			p <- ggplot() +
			  geom_violin(df_int, mapping = aes(x = Channel, y = log10(Intensity), fill = Channel), size = .25) +
			  geom_boxplot(df_int, mapping = aes(x = Channel, y = log10(Intensity)), width = 0.2, lwd = .2, fill = "white") +
			  stat_summary(df_int, mapping = aes(x = Channel, y = log10(Intensity)), fun.y = "mean", geom = "point",
			               shape = 23, size = 2, fill = "white", alpha = .5) +
			  labs(title = expression("Reporter ions"), x = expression("Channel"), y = expression("Intensity ("*log[10]*")")) +
			  geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
			  theme_psm_violin

			ggsave(file.path(dat_dir, "PSM\\Violin\\bf_olm", filename), p, width = Width, height = Height, units = "in")
		}

	}

}


#' Locates the positions of outliers
#' 
#' @param df A data frame containing the PSM table from Mascot.
#' @param range_colRatios The range of columns.
#' @return A data frame.
#' @examples
#' locate_outliers(df, 2:3)
locate_outliers <- function (df, range_colRatios) {
	for(col_index in range_colRatios) {
		counts <- colSums(!is.na(df[col_index]))
		if(counts > 25) df[, col_index] <- Rosner_outliers(df[, col_index])
		else if(counts > 2) df[, col_index] <- Dixon_outliers(df[, col_index])
	}

	return(df)
}


#' Outlier removals with Rosner's method
Rosner_outliers <- function(x) {

	if (length(unique(x)) < 5) return(x)
	# up to 9-number outliers; may get warnings with NA being an outlier
	gofOutlier_obj <- rosnerTest(as.numeric(x), 9)

	if(gofOutlier_obj$n.outliers > 0) {
		Index <- with(gofOutlier_obj$all.stat, Obs.Num[Outlier == TRUE])
		x[Index] <- NA
	}

	return(x)
}


#' Outlier removals with Dixon's method
Dixon_outliers <- function(x) {
	# x = c(0.0000000, 0.0000000, 1.0271542, 0.0000000, 0.2080097)
	# x = c(0.0000000, 0.0000000, NA, 0.0000000, 0.2080097)
	# x = c(0.0000000, 0.0000000, 0.0000000, 0.2080097)
	# x = c(NA, NA, NA, 0.2080097)

	newx <- x[!is.na(x)]
	len_newx <- length(newx)
	uni_newx <- length(unique(newx))

	if (len_newx > 2 & uni_newx > 1) {
		gofOutlier_obj <- dixon.test(as.numeric(x), type = 0)

		while(gofOutlier_obj$p.value < 0.05) {
			if (grepl("^highest", gofOutlier_obj$alternative)) x[which.max(x)] <- NA else x[which.min(x)] <- NA

			newx <- x[!is.na(x)]
			len_newx <- length(newx)
			uni_newx <- length(unique(newx))

			if (len_newx > 2 & uni_newx > 1)
			  gofOutlier_obj <- dixon.test(as.numeric(x), type = 0) else gofOutlier_obj$p.value <- 1
		}
	}

	return (x)
}


#' Outlier removals with Grubbs's method
Grubbs_outliers <- function(x, type = 10) {
	newx <- x[!is.na(x)]
	len_newx <- length(newx)
	uni_newx <- length(unique(newx))

	if (len_newx > 2 & uni_newx > 1) {
		gofOutlier_obj <- grubbs.test(as.numeric(x, type))

		while(gofOutlier_obj$p.value < 0.05) {
			if (grepl("^highest", gofOutlier_obj$alternative)) x[which.max(x)] <- NA else x[which.min(x)] <- NA

			newx <- x[!is.na(x)]
			len_newx <- length(newx)
			uni_newx <- length(unique(newx))

			if (len_newx > 2 & uni_newx > 1)
			  gofOutlier_obj <- grubbs.test(as.numeric(x), type = type) else gofOutlier_obj$p.value <- 1
		}
	}

	return (x)
}


#' Cleans Up PSM results
#'
#' \code{cleanupPSM} cleans up PSM tables with the option of outlier removals.
#' The outlier removals will be assessed at the basis of per peptide per TMT
#' channel.
#'
#' Dixon's method will be used when 2 < n <= 25; Rosner's method will be used
#' when n > 25
#'
#' @param rm_outliers Logical; if TRUE, performs outlier removals.
#' @import dplyr tidyr
#' @importFrom stringr str_split
#' @importFrom outliers dixon.test
#' @importFrom EnvStats rosnerTest
#' @export
cleanupPSM <- function(rm_outliers = FALSE) {

	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 5000000)

	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	TMT_plex <- TMT_plex(label_scheme_full)

	filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
	                      pattern = "^TMT.*LCMS.*\\.csv$")

	for (i in seq_along(filelist)) {
		df <- read.csv(file.path(dat_dir, "PSM\\cache", filelist[i]), check.names = FALSE,
		               header = TRUE, comment.char = "#")

		if(TMT_plex == 0) {
			# lable-free data
		  # re-save ".csv" as ".txt"
			fn <- paste0(gsub(".csv", "", filelist[i]), "_Clean.txt")
			write.table(df, file.path(dat_dir, "PSM\\cache", fn), sep = "\t", col.names = TRUE,
			            row.names = FALSE)
			cat(filelist[i], "processed\n")

			next
		}

		# remove all "-1" ratio columns
		N <- sum(grepl("^R[0-9]{3}", names(df)))

		df <- df %>%
			dplyr::mutate(n = rowSums(.[, grep("^R[0-9]{3}", names(.))] == -1)) %>%
			dplyr::filter(n != N) %>%
			dplyr::select(-n)

		channelInfo <- channelInfo(label_scheme, set_idx =
		                  as.integer(gsub("TMTset(\\d+)_.*", "\\1", filelist[i])))


		# add a column of "R126"
		pos_af <- min(grep("^R1[0-9]{2}", names(df)))

		df$R126 <- 1
		df <- cbind.data.frame(df[, 1:(pos_af-1)], R126 = df$R126, df[, (pos_af):(ncol(df)-1)]) %>%
				dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^I[0-9]{3}", names(.) )])) > 0) %>%
				dplyr::mutate_at(.vars = which(names(.)=="I126")-1+channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
				dplyr::mutate_at(.vars = which(names(.)=="R126")-1+channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) # note that "> 1" not "0"

		if (rm_outliers) {
			dfw_split <- df %>%
				dplyr::select(grep("^I[0-9]{3}", names(.))) %>%
				dplyr::mutate(RM = rowMeans(.[, grep("^I[0-9]{3}", names(.))[channelInfo$labeledChannels]],
				                            na.rm = TRUE)) %>%
				dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ log2(.x/RM)) %>%
				dplyr::select(-c("RM")) %>%
				`colnames<-`(gsub("I", "X", names(.))) %>%
				dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA)) %>%
				dplyr::bind_cols(df[, c("psm_index", "pep_seq")], .) %>%
				split(., .$pep_seq, drop = TRUE)

			range_colRatios <- grep("^X[0-9]{3}", names(dfw_split[[1]]))

			dfw_split <- do.call("rbind", lapply(dfw_split, locate_outliers, range_colRatios)) %>%
					dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA)) %>%
					tidyr::unite(pep_seq_i, pep_seq, psm_index, sep = ".") %>%
					dplyr::mutate_at(.vars = grep("^X[0-9]{3}", names(.)), ~ replace(.x, !is.na(.x), 1))

			df <- df %>%
					tidyr::unite(pep_seq_i, pep_seq, psm_index, sep = ".") %>%
					dplyr::left_join(., dfw_split, by = "pep_seq_i") %>%
					tidyr::separate(pep_seq_i, into = c("pep_seq", "psm_index"), remove = TRUE) %>%
					dplyr::select(-c("psm_index"))

			rm(dfw_split, range_colRatios)

			df[, grepl("^I[0-9]{3}", names(df))] <-
			  purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
			              as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
			  dplyr::bind_rows()

			df[, grepl("^R[0-9]{3}", names(df))] <-
			  purrr::map2(as.list(df[, grepl("^R[0-9]{3}", names(df))]),
			              as.list(df[, grepl("^X[0-9]{3}", names(df))]), `*`) %>%
			  dplyr::bind_rows()

			df <- cbind.data.frame(raw_file = df[, c("raw_file")],
					df[, !grepl("^R[0-9]{3}|^I[0-9]{3}|^X[0-9]{3}|^raw_file$", names(df))],
					df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df))]) %>%
					dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1) %>% # "> 1" as "R126 == 1"
					dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 + channelInfo$emptyChannels, ~ replace(.x, , 0))
		} else {
			df <- cbind.data.frame(raw_file = df[, c("raw_file")],
				df[, !grepl("^R[0-9]{3}|^I[0-9]{3}|^psm_index$|^raw_file$", names(df))],
				df[, grepl("^R[0-9]{3}|^I[0-9]{3}", names(df))]) %>%
			  dplyr::mutate_at(.vars = which(names(.)=="I126") - 1 +
				                   channelInfo$emptyChannels, ~ replace(.x, , 0)) %>%
				dplyr::mutate_at(.vars = which(names(.)=="R126") - 1 +
				                   channelInfo$emptyChannels, ~ replace(.x, , NA)) %>%
				dplyr::filter(rowSums(!is.na(.[, grep("^R[0-9]{3}", names(.) )])) > 1)

		}

		fn <- paste0(gsub(".csv", "", filelist[i]), "_Clean.txt")
		write.table(df, file.path(dat_dir, "PSM\\cache", fn), sep = "\t", col.names = TRUE,
		            row.names = FALSE)
		cat(filelist[i], "processed\n")
	}

}


#'Median-centering normalization of PSM data
#'
#'\code{mcPSM} adds fields of \code{log2_R, N_log2_R and N_I} to PSM tables.
#'
#'@import dplyr tidyr purrr
#'@importFrom magrittr %>%
mcPSM <- function(df, set_idx) {
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  label_scheme_sub <- label_scheme[label_scheme$TMT_Set == set_idx & 
                                     label_scheme$LCMS_Injection == 1, ]
  
  channelInfo <- channelInfo(label_scheme_sub, set_idx)
  
  dfw <- df[rowSums(!is.na(df[, grepl("^R[0-9]{3}", names(df)), drop = FALSE])) > 1, ] %>%
    dplyr::arrange(pep_seq, prot_acc) %>%
    dplyr::filter(pep_expect <=  0.10) %>%
    dplyr::mutate_at(.vars = which(names(.) == "I126") - 1 +
                       channelInfo$emptyChannels, ~ replace(.x, , NaN))
  
  col_sample <- grep("^I[0-9]{3}", names(dfw))
  
  if(length(channelInfo$refChannels) > 0) {
    ref_index <- channelInfo$refChannels
  } else {
    ref_index <- channelInfo$labeledChannels
  }
  
  dfw <- sweep(dfw[, col_sample], 1,
               rowMeans(dfw[, col_sample[ref_index], drop = FALSE], na.rm = TRUE), "/") %>%
    log2(.) %>%
    `colnames<-`(gsub("I", "log2_R", names(.)))	%>%
    cbind(dfw, .) %>%
    dplyr::mutate_at(.vars = grep("[I|R][0-9]{3}", names(.)), ~ replace(.x, is.infinite(.), NA))
  
  col_log2Ratio <- grepl("^log2_R[0-9]{3}", names(dfw))
  cf <- apply(dfw[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  
  dfw <- sweep(dfw[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- sweep(dfw[, grepl("^I[0-9]{3}", names(dfw)), drop = FALSE], 2, 2^cf, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(dfw, .)
  
  dfw <- dfw %>%
    reorderCols(endColIndex = grep("[RI][0-9]{3}", names(dfw)), col_to_rn = "pep_seq_mod") %>%
    na_zeroIntensity()
}


#' Annotates PSM results
#'
#' \code{annotPSM} adds fields of annotation to PSM tables.
#'
#'@inheritParams load_expts
#' @param rm_krts Logical; if TRUE, removes keratin entries from the output.
#' @param plot_violins Logical; if TRUE, prepares the violin plots of
#'   reporter-ion intensities.
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
annotPSM <- function(expt_smry = "expt_smry.xlsx", rm_krts = FALSE, plot_violins = TRUE) {
	dir.create(file.path(dat_dir, "PSM\\Individual_mods"), recursive = TRUE, showWarnings = FALSE)

	old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	options(max.print = 5000000)

	hd_fn <- list.files(path = file.path(dat_dir, "PSM\\cache"),
	                    pattern = "^F\\d+_header.txt$")
	assign("df_header", readLines(file.path(dat_dir, "PSM\\cache", hd_fn[1])))

	load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	n_TMT_sets <- n_TMT_sets(label_scheme_full)
	TMT_plex <- TMT_plex(label_scheme_full)

	filelist <- list.files(
								path = file.path(dat_dir, "PSM\\cache"),
	              pattern = "^TMT.*LCMS.*_Clean.txt$"
							) %>%
		reorder_files(., n_TMT_sets)

	for(set_idx in seq_len(n_TMT_sets)){
	  sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]

		out_fn <- data.frame(Filename =
		                       do.call('rbind', strsplit(as.character(sublist),
		                              '.txt', fixed = TRUE))) %>%
			dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))

	  channelInfo <- channelInfo(label_scheme, set_idx)

	  # injections under the same TMT experiment
	  for (injn_idx in seq_along(sublist)) {
			df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
			               check.names = FALSE, header = TRUE, sep = "\t",
			               comment.char = "#") %>%
				dplyr::mutate(pep_seq_mod = as.character(pep_seq)) %>%
				dplyr::select(which(names(.) == "pep_seq_mod"),
				              which(names(.) != "pep_seq_mod"))

			acc_type <- find_acctype()
			species <- find_df_species(df, acc_type)

		  label_scheme_full$Species <- species
		  label_scheme$Species <- species
		  save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.Rdata"))
		  save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))
		  write.table(label_scheme[1, c("Accession_Type", "Species")],
		              file.path(dat_dir, "acctype_sp.txt"), sep = "\t",
		              col.names = TRUE, row.names = FALSE)

		  load_dbs()

			if(rm_krts) {
				krts <- dbs$prn_annot %>%
					dplyr::filter(grepl("^krt[0-9]+", .$gene, ignore.case = TRUE)) %>%
					dplyr::filter(!is.na(.[[acc_type]])) %>%
					dplyr::select(acc_type) %>%
					unlist()

				df <- df %>%
					dplyr::filter(! .$prot_acc %in% krts)
			}

			# combine peptide sequences and modifications
			fixed_mods <- df_header[((grep("Fixed modifications", df_header))[1] + 3) :
															((grep("Variable modifications", df_header))[1] - 2)] %>%
				gsub("\"", "", ., fixed = TRUE) %>%
				data.frame() %>%
				tidyr::separate(".", sep = ",", c("Mascot_abbr", "Descption", "Delta_mass"))

			var_mods <- df_header[((grep("Variable modifications", df_header))[1] + 3) :
														((grep("Search Parameters", df_header))[1] - 2)] %>%
				gsub("\"", "", ., fixed = TRUE) %>%
				data.frame() %>%
				tidyr::separate(".", sep = ",", extra = "drop", c("Mascot_abbr", "Descption", "Delta_mass")) %>%
				dplyr::mutate(Abbr = Descption) %>%
				dplyr::mutate(Abbr = try(gsub("Acetyl (Protein N-term)", "_", .$Abbr, fixed = TRUE))) %>%
				dplyr::mutate(Abbr = try(gsub("Deamidated (N)", "n", .$Abbr, fixed = TRUE))) %>%
				dplyr::mutate(Abbr = try(gsub("Gln->pyro-Glu (N-term Q)", "q", .$Abbr, fixed = TRUE))) %>%
				dplyr::mutate(Abbr = try(gsub("Oxidation (M)", "m", .$Abbr, fixed = TRUE))) %>%
				dplyr::mutate(Abbr = try(gsub("Pyro-carbamidomethyl (N-term C)", "c", .$Abbr, fixed = TRUE))) %>%
				dplyr::mutate(Abbr = try(gsub("Carbamidomethyl (C)", "C", .$Abbr, fixed = TRUE))) %>%
				dplyr::mutate(Filename = Descption) %>%
				dplyr::mutate(Filename = try(gsub("Acetyl (Protein N-term)", "N_term_Ac", .$Filename, fixed = TRUE))) %>%
				dplyr::mutate(Filename = try(gsub("Deamidated (N)", "Deamidated_N", .$Filename, fixed = TRUE))) %>%
				dplyr::mutate(Filename = try(gsub("Gln->pyro-Glu (N-term Q)", "N_term_Q", .$Filename, fixed = TRUE))) %>%
				dplyr::mutate(Filename = try(gsub("Oxidation (M)", "Oxidation_M", .$Filename, fixed = TRUE))) %>%
				dplyr::mutate(Filename = try(gsub("Pyro-carbamidomethyl (N-term C)", "N_term_C", .$Filename, fixed = TRUE)))

			# phosphorylation mods
			phospho_mods <- var_mods %>% dplyr::filter(grepl("^Phospho", Descption))
			if(nrow(phospho_mods) > 0) var_mods <- var_mods %>%
			  dplyr::filter(!grepl("^Phospho", Descption))

			# (1) Non-acetylated N-term variable modifications
			NT_var_mods <- var_mods %>%
			  dplyr::filter(grepl("N-term", Descption, fixed = TRUE)) %>%
				dplyr::filter(!grepl("Acetyl (Protein N-term)", Descption, fixed = TRUE))

			for (mod in NT_var_mods$Mascot_abbr) {
				df_sub <- df[grepl(mod, df$pep_var_mod_pos), ]

				if (nrow(df_sub) > 0) {
					locales <- 1 # always "1" since a N-term mod

					lowers <- substr(df_sub$pep_seq_mod, locales, locales) %>% tolower()
					substr(df_sub$pep_seq_mod, locales, locales) <- lowers

					fn <- var_mods[var_mods$Mascot_abbr == mod, "Filename"]
					try(write.table(df_sub, file.path(dat_dir, "PSM\\Individual_mods", paste0(fn, ".txt")), 
					                sep = "\t", col.names = TRUE, row.names = FALSE))

					df <- rbind(df[!grepl(mod, df$pep_var_mod_pos), ], df_sub)
				}
			}

			rm (NT_var_mods, mod, df_sub, locales, lowers, fn)

			# (2) non-N-term variable modifications
			nNT_var_mods <- var_mods[!grepl("N-term", var_mods$Descption), ]

			for (mod in nNT_var_mods$Mascot_abbr) {
				df_sub <- df[grepl(mod, df$pep_var_mod_pos), ]

				if (nrow(df_sub) > 0) {
					# need to find all matches
				  pos_matrix  <- gregexpr(mod, df_sub$pep_var_mod_pos) %>%
						plyr::ldply(., rbind) %>%
				    # "-2" to account for the characters, "0.", "1."..., in 'pep_var_mod_pos'
				    purrr::map(function(x) {x - 2}) %>%
						data.frame(check.names = FALSE)

					for (k in 1:ncol(pos_matrix)) {
						rows <- !is.na(pos_matrix[, k])
						locales <- pos_matrix[rows, k]

						lowers <- substr(df_sub$pep_seq_mod[rows], locales, locales) %>% tolower()
						substr(df_sub$pep_seq_mod[rows], locales, locales) <- lowers
					}
				  
				  fn <- var_mods[var_mods$Mascot_abbr == mod, "Filename"]
				  try(write.table(df_sub, file.path(dat_dir, "PSM\\Individual_mods", paste0(fn, ".txt")), 
				                  sep = "\t", col.names = TRUE, row.names = FALSE))

					df <- rbind(df[!grepl(mod, df$pep_var_mod_pos), ], df_sub)
				}
			}

			rm (nNT_var_mods, mod, df_sub, pos_matrix, k, rows, locales, lowers, fn)

			# (3) phosphopeptides
			if(nrow(phospho_mods) > 0) {
				for (mod in phospho_mods$Mascot_abbr) {
					df_sub <- df[grepl(mod, df$pep_var_mod_pos), ]

					if (nrow(df_sub) > 0) {
						# need to find all matches
					  pos_matrix  <- gregexpr(mod, df_sub$pep_var_mod_pos) %>%
							plyr::ldply(., rbind) %>%
							purrr::map(function(x) {x - 2}) %>%
							data.frame(check.names = FALSE)

						for (k in 1:ncol(pos_matrix)) {
							rows <- !is.na(pos_matrix[, k])
							locales <- pos_matrix[rows, k]

							lowers <- substr(df_sub$pep_seq_mod[rows], locales, locales) %>%
							  tolower()
							substr(df_sub$pep_seq_mod[rows], locales, locales) <- lowers

							fn <- phospho_mods[phospho_mods$Mascot_abbr == mod, "Filename"]
							try(
							  write.table(df_sub, file.path(dat_dir, "PSM\\Individual_mods",paste0(fn, ".txt")), 
							              sep = "\t", col.names = TRUE, row.names = FALSE)							  
							)
						}

						df <- rbind(df[!grepl(mod, df$pep_var_mod_pos), ], df_sub)
					}
				}

				rm (phospho_mods, mod, df_sub, pos_matrix, k, rows, locales, lowers, fn)
			}

			# (4) replace numeric symbols with mod names for non-ace N-term mods: for example 2 <- n
			# nace_var_mods <- var_mods[!grepl("Acetyl (Protein N-term)", var_mods$Descption, fixed = TRUE), ]
			# for (mod in nace_var_mods$Mascot_abbr)
			#   df$pep_seq_mod <- gsub(mod, nace_var_mods[nace_var_mods$Mascot_abbr == mod, ]$Abbr, df$pep_seq_mod)

			# (5) Add "_" for "Acetyl (Protein N-term)"
			nt_ace_var_mods <- var_mods[grepl("Acetyl (Protein N-term)",
			                                  var_mods$Descption, fixed = TRUE), ]

			for (mod in nt_ace_var_mods$Mascot_abbr) {
				df_sub <- df[grepl(mod, df$pep_var_mod_pos), ]

				if (nrow(df_sub) > 0) {
					# df_sub$pep_seq_mod <- paste0(paste0(substring(df_sub$pep_seq_mod, 1, 1), nt_ace_var_mods$Abbr),
					#                              substring(df_sub$pep_seq_mod, 2)) # insert "_"
				  df_sub$pep_seq_mod <- paste0("_", df_sub$pep_seq_mod)
					fn <- nt_ace_var_mods[nt_ace_var_mods$Mascot_abbr == mod, "Filename"]
					try(write.table(df_sub, file.path(dat_dir, "PSM\\Individual_mods", paste0(fn, ".txt")), 
					                sep = "\t", col.names = TRUE, row.names = FALSE))

					df <- rbind(df[!grepl(mod, df$pep_var_mod_pos), ], df_sub)
				}
			}

			rm (nt_ace_var_mods, mod, df_sub, fn)

			# add "pep_res_before" and "pep_res_after"
			df <- df %>%
				dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = ".")) %>%
				dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = "."))

			if(TMT_plex > 0) df <- mcPSM(df, set_idx)

			# clean up cRAP accessions
			df$prot_acc <- df$prot_acc %>%
				gsub("\\|$", "", .) %>%
				gsub(".*\\|", "", .)
			# df$prot_acc <- gsub("\\.[0-9]+|.*\\|", "", df$prot_acc)

			write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
			            sep = "\t", col.names = TRUE, row.names = FALSE)

			if (plot_violins & TMT_plex > 0) {
				TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()
				Levels <- TMT_levels %>% gsub("^TMT-", "I", .)
				df_int <- df[, grepl("^I[0-9]{3}", names(df))] %>%
					tidyr::gather(key = "Channel", value = "Intensity") %>%
					dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
				  dplyr::filter(!is.na(Intensity))
				rm(Levels)
				
				mean_int <- df_int %>% 
				  dplyr::group_by(Channel) %>% 
				  dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
				  dplyr::mutate(Intensity = round(Intensity, digit = 1))

				Width <- 8
				Height <- 8
				filename <- paste0("Post_outlier_removals_set_", set_idx, "_inj_", injn_idx)

				p <- ggplot() +
					geom_violin(df_int, mapping=aes(x=Channel, y=log10(Intensity), fill=Channel), size=.25) +
					geom_boxplot(df_int, mapping=aes(x=Channel, y=log10(Intensity)), width=0.2, lwd=.2, fill="white") +
					stat_summary(df_int, mapping=aes(x=Channel, y=log10(Intensity)), fun.y="mean", geom="point",
					             shape=23, size=2, fill="white", alpha=.5) +
					labs(title=expression("Reporter ions"), x=expression("Channel"), y=expression("Intensity ("*log[10]*")")) + 
				  geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
					theme_psm_violin

				ggsave(file.path(dat_dir, "PSM\\Violin", paste0(filename, ".png")), p, width = Width, height = Height, units = "in")
			}

	  }

	}
}


#'Splits PSM tables
#'
#'\code{splitPSM} splits the PSM outputs after \code{rmPSMHeaders()}. It
#'separates PSM data by TMT experiment and LC/MS injection.
#'
#'@param fasta Character string(s) to the name(s) of fasta file(s) with
#'  prepended directory path. The \code{fasta} database(s) need to match those
#'  used in MS/MS ion search.
#'@param pep_unique_by A character string for use in the filtration of peptide
#'  data from \code{MaxQuant}. The choice is in one of \code{c("group",
#'  "protein", "none")} At the \code{group} default, only peptides unique by
#'  protein groups will be kept. At \code{protein}, only peptides unique by
#'  protein entries will be kept. At \code{none}, all peptides will be kept.
#'@param corrected_int Logical. If TRUE, values under columns "Reporter
#'  intensity corrected" in MaxQuant PSM results (\code{msms.txt}) will be used.
#'  Otherwise, "Reporter intensity" values without corrections will be used.
#'@param rm_reverses Logical; if TRUE, removes \code{Reverse} entries from
#'  MaxQuant peptide results.
#'@inheritParams splitPSM
#' @examples
#' splitPSM(rptr_intco = 1000, rm_craps = TRUE, plot_violins = TRUE)
#'
#'@import dplyr tidyr
#'@importFrom stringr str_split
#'@importFrom magrittr %>%
#'@export
splitPSM_mq <- function(fasta = NULL, pep_unique_by = "group", corrected_int = FALSE, rptr_intco = 1000, 
                        rm_craps = FALSE, rm_reverses = TRUE, plot_violins = TRUE, ...) {
  
  add_mascot_headers <- function (df) {
    df %>% 
      dplyr::mutate(pep_scan_title = RAW_File) %>% 
      dplyr::mutate(prot_hit_num = NA) %>% 
      dplyr::mutate(prot_family_member = NA) %>% 
      dplyr::mutate(prot_score = NA) %>% 
      dplyr::mutate(prot_matches = NA) %>% 
      dplyr::mutate(prot_matches_sig = NA) %>% 
      dplyr::mutate(prot_sequences = NA) %>% 
      dplyr::mutate(prot_sequences_sig = NA) %>% 
      dplyr::mutate(pep_query = NA) %>% 
      dplyr::mutate(pep_rank = NA) %>% 
      dplyr::mutate(pep_isbold = NA) %>% 
      dplyr::mutate(pep_exp_mr = NA) %>% 
      dplyr::mutate(pep_exp_z = NA) %>% 
      dplyr::mutate(pep_calc_mr = NA) %>% 
      dplyr::mutate(pep_exp_mz = NA) %>% 
      dplyr::mutate(pep_delta = NA) %>% 
      dplyr::mutate(pep_miss = NA) %>% 
      dplyr::mutate(pep_var_mod = NA) %>% 
      dplyr::mutate(pep_var_mod_pos = NA) %>% 
      dplyr::mutate(pep_summed_mod_pos = NA) %>% 
      dplyr::mutate(pep_local_mod_pos = NA)
  }
  
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  
  on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  load(file = file.path(dat_dir, "fraction_scheme.Rdata"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)

  filelist <- list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$")
  
  if (rlang::is_empty(filelist)) 
    stop(paste("No PSM files were found under", file.path(dat_dir), 
               "\nMake sure the names of PSM files start with `msms`."), call. = FALSE)
  
  df <- purrr::map(file.path(dat_dir, filelist), read.csv, 
                   check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::bind_rows() %>% 
    dplyr::rename(ID = id)
  
  df <- filters_in_call(df, ...)

  if (corrected_int) {
    df <- df %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}\\d+$", names(.)))
  } else {
    df <- df %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}corrected\\s{1}\\d+$", names(.)))
  }
  
  if (rm_craps) {
    df <- df %>% dplyr::filter(!grepl("^CON_", .[["Proteins"]]))
  }
  
  if (rm_reverses) {
    df <- df %>% dplyr::filter(.$Reverse != "+")
  }
  
  df <- df %>% 
    `names_pos<-`(grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+$", names(.)), 
                  gsub("TMT-", "I", as.character(TMT_levels)))
  
  if (TMT_plex == 11) {
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
  
  if (TMT_plex > 0) {
    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) %>% 
      dplyr::rename(
        pep_seq = Sequence, 
        prot_acc = Proteins, 
        pep_miss = `Missed cleavages`,
        pep_score = `Score`,
        pep_expect = PEP, 
        pep_var_mod = Modifications, 
        RAW_File = `Raw file`
      ) 
    
    if (pep_unique_by == "group") {
      df <- df %>% 
        dplyr::mutate(pep_isunique = ifelse(grepl(";", `Protein group IDs`), 0, 1))
    } else if (pep_unique_by == "protein") {
      df <- df %>% 
        dplyr::mutate(pep_isunique = ifelse(grepl(";", prot_acc), 0, 1))
    }
    
    df <- df %>% 
      dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc)) %>% 
      add_mascot_headers()
  }
  
  prn_acc <- parse_acc(df)
  label_scheme_full$Accession_Type <- prn_acc
  label_scheme$Accession_Type <- prn_acc
  save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.Rdata"))
  save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))
  
  df <- annotPrndesc(df, fasta)
  
  df <- df %>% 
    dplyr::filter(!duplicated(pep_seq)) %>% 
    dplyr::select(c("pep_seq", "prot_acc")) %>% 
    annotPeppos_mq(fasta) %>% 
    dplyr::select(-prot_acc) %>% 
    dplyr::right_join(., df, by = "pep_seq")
  

  df <- cbind.data.frame(
    df[, grepl("^[a-z]", names(df))], 
    df[, grepl("^[A-Z]", names(df)) & !grepl("^[IR][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[R][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[I][0-9]{3}[NC]*", names(df))]
  )

  if(length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      dplyr::filter(pep_expect <=  0.10) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>%
      # a special case of redundant entries from Mascot
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::filter(pep_expect <=  0.10) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` under the experimental summary file."),
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")
  
  
  for (i in seq_along(df_split)) {
    df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)
    
    out_fn <- fn_lookup %>%
      dplyr::filter(TMT_inj == names(df_split)[i]) %>%
      dplyr::select(filename) %>%
      unique() %>%
      unlist() %>%
      paste0(., ".csv")
    
    df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
    
    write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
    
    if (plot_violins & TMT_plex > 0) {
      dir.create(file.path(dat_dir, "PSM\\Violin\\bf_olm"), recursive = TRUE, showWarnings = FALSE)
      
      TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()
      Levels <- TMT_levels %>% gsub("^TMT-", "I", .)
      df_int <- df_split[[i]][, grepl("^I[0-9]{3}", names(df_split[[i]]))] %>%
        tidyr::gather(key = "Channel", value = "Intensity") %>%
        dplyr::mutate(Channel = factor(Channel, levels = Levels))
      rm(Levels)
      
      mean_int <- df_int %>% 
        dplyr::group_by(Channel) %>% 
        dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
        dplyr::mutate(Intensity = round(Intensity, digit = 1))
      
      Width <- 8
      Height <- 8
      
      filename <- gsub("\\.csv", "\\.png", out_fn)
      p <- ggplot() +
        geom_violin(df_int, mapping = aes(x = Channel, y = log10(Intensity), fill = Channel), size = .25) +
        geom_boxplot(df_int, mapping = aes(x = Channel, y = log10(Intensity)), width = 0.2, lwd = .2, fill = "white") +
        stat_summary(df_int, mapping = aes(x = Channel, y = log10(Intensity)), fun.y = "mean", geom = "point",
                     shape = 23, size = 2, fill = "white", alpha = .5) +
        labs(title = expression("Reporter ions"), x = expression("Channel"), y = expression("Intensity ("*log[10]*")")) +
        geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
        theme_psm_violin
      
      ggsave(file.path(dat_dir, "PSM\\Violin\\bf_olm", filename), p, width = Width, height = Height, units = "in")
    }
  }
  
}


#' Annotates PSM results
#'
#' \code{annotPSM} adds fields of annotation to PSM tables.
#'
#'@inheritParams load_expts
#'@inheritParams annotPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
annotPSM_mq <- function(expt_smry = "expt_smry.xlsx", rm_krts = FALSE, plot_violins = TRUE) {

  dir.create(file.path(dat_dir, "PSM\\Individual_mods"), recursive = TRUE, showWarnings = FALSE)
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  
  options(max.print = 5000000)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(., n_TMT_sets)
  
  for(set_idx in seq_len(n_TMT_sets)){
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename =
                           do.call('rbind', strsplit(as.character(sublist),
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    # LCMS injections under the same TMT experiment
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") %>%
        dplyr::rename(pep_seq_mod = `Modified sequence`) %>%
        dplyr::select(which(names(.) == "pep_seq_mod"),
                      which(names(.) != "pep_seq_mod"))
      
      acc_type <- find_acctype()
      species <- find_df_species(df, acc_type)
      
      label_scheme_full$Species <- species
      label_scheme$Species <- species
      save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.Rdata"))
      save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))
      write.table(label_scheme[1, c("Accession_Type", "Species")],
                  file.path(dat_dir, "acctype_sp.txt"), sep = "\t",
                  col.names = TRUE, row.names = FALSE)
      
      load_dbs()
      
      if (rm_krts) df <- rm_krts(df, acc_type)

      df <- df %>% 
        dplyr::mutate(pep_seq_mod = gsub("^_(.*)_$", "\\1", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("^\\(ac\\)", "_", pep_seq_mod)) %>% 
        # dplyr::mutate(pep_seq_mod = gsub("^\\(gl\\)", "q", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("A\\([^\\)]+\\)", "a", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("C\\([^\\)]+\\)", "c", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("D\\([^\\)]+\\)", "d", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("E\\([^\\)]+\\)", "e", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("F\\([^\\)]+\\)", "f", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("G\\([^\\)]+\\)", "g", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("H\\([^\\)]+\\)", "h", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("I\\([^\\)]+\\)", "i", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("K\\([^\\)]+\\)", "k", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("L\\([^\\)]+\\)", "l", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("M\\([^\\)]+\\)", "m", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("N\\([^\\)]+\\)", "n", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("P\\([^\\)]+\\)", "p", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("Q\\([^\\)]+\\)", "q", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("R\\([^\\)]+\\)", "r", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("s\\([^\\)]+\\)", "s", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("T\\([^\\)]+\\)", "t", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("V\\([^\\)]+\\)", "v", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("W\\([^\\)]+\\)", "w", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("Y\\([^\\)]+\\)", "y", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = ".")) %>% 
        dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = "."))
      
      if(TMT_plex > 0) df <- mcPSM(df, set_idx)
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_violins & TMT_plex > 0) {
        TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()
        Levels <- TMT_levels %>% gsub("^TMT-", "I", .)
        df_int <- df[, grepl("^I[0-9]{3}", names(df))] %>%
          tidyr::gather(key = "Channel", value = "Intensity") %>%
          dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
          dplyr::filter(!is.na(Intensity))
        rm(Levels)
        
        mean_int <- df_int %>% 
          dplyr::group_by(Channel) %>% 
          dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
          dplyr::mutate(Intensity = round(Intensity, digit = 1))
        
        Width <- 8
        Height <- 8
        filename <- paste0("Post_outlier_removals_set_", set_idx, "_inj_", injn_idx)
        
        p <- ggplot() +
          geom_violin(df_int, mapping=aes(x=Channel, y=log10(Intensity), fill=Channel), size=.25) +
          geom_boxplot(df_int, mapping=aes(x=Channel, y=log10(Intensity)), width=0.2, lwd=.2, fill="white") +
          stat_summary(df_int, mapping=aes(x=Channel, y=log10(Intensity)), fun.y="mean", geom="point",
                       shape=23, size=2, fill="white", alpha=.5) +
          labs(title=expression("Reporter ions"), x=expression("Channel"), y=expression("Intensity ("*log[10]*")")) + 
          geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
          theme_psm_violin
        
        ggsave(file.path(dat_dir, "PSM\\Violin", paste0(filename, ".png")), p, width = Width, height = Height, units = "in")
      }
      
    }
    
  }
}


#'Reports PSM results
#'
#'\code{normPSM} reports
#'\code{\href{https://www.ebi.ac.uk/pride/help/archive/search/tables}{PSM}}
#'results from \code{\href{https://en.wikipedia.org/wiki/Tandem_mass_tag}{TMT}}
#'experiments.
#'
#'The \code{PSM} files should be present under the working directory,
#'\code{"dat_dir"}, specified by end users.
#'
#'In each primary output file, "\code{...PSM_N.txt}", values under columns
#'\code{log2_R...} are relative to the \code{reference(s)} within each multiplex
#'TMT set or to the row means if no \code{reference(s)} are present. Values
#'under columns \code{N_log2_R...} are \code{log2_R...} with median-centering
#'alignment. Values under columns \code{N_I...} are normalized
#'\code{reporter-ion intensity}. Character strings under \code{pep_seq_mod}
#'denote peptide sequences with applicable variable modifications.
#'
#'@section \code{Mascot}: End users will export \code{PSM} data from
#'  \code{\href{https://http://www.matrixscience.com/}{Mascot}} at a \code{.csv}
#'  format. The header information should be included during the \code{.csv}
#'  export. The file name(s) should be defaulted by
#'  \code{\href{https://http://www.matrixscience.com/}{Mascot}}: starting with
#'  the letter \code{'F'}, followed by a six-digit number without space and
#'  ended with a \code{'.csv'} extension \code{(e.g., F004453.csv)}.
#'
#'  See \code{\link{normPrn}} for the description of column keys in the output.
#'
#'@section \code{MaxQuant}: End users will copy over \code{msms.txt}
#'  file(s) from \code{\href{https://www.maxquant.org/}{MaxQuant}} to the
#'  \code{dat_dir} directory. In the case of multiple \code{msms.txt} files for
#'  processing, the file names need to be compiled in that they all start with
#'  \code{'msms'} and end with a \code{'.txt'} extension.
#'
#'@param ... Additional parameters for row filtration of data. 
#'@inheritParams load_expts
#'@inheritParams splitPSM
#'@inheritParams splitPSM_mq
#'@inheritParams cleanupPSM
#'@inheritParams annotPSM
#'@seealso \code{\link{normPep}} for peptides and \code{\link{normPrn}} for
#'  proteins.
#'@return Outputs under \code{dat_dir\\PSM}. Primary results are in
#'  \code{TMTset1_LCMSinj1_PSM_N.txt, TMTset2_LCMSinj1_PSM_N.txt, ...} The
#'  indeces of TMT experiment and LC/MS injection are indicated in the file
#'  names.
#'
#' @examples
#' \dontrun{
#' # Mascot example
#' # set up a working directory
#' dat_dir <- c("C:\\The\\First\\Example")
#'
#' # copy fasta
#' library(proteoQDA)
#' copy_refseq_hs("~\\proteoQ\\dbs\\refseq")
#' copy_refseq_mm("~\\proteoQ\\dbs\\refseq")
#'
#' # copy Mascot PSM data
#' cptac_csv_1(dat_dir)
#'
#' # copy "expt_smry.xlsx" and "frac_smry.xlsx"
#' cptac_expt_1(dat_dir)
#' cptac_frac_1(dat_dir)
#'
#' # load experiments
#' library(proteoQ)
#' load_expts()
#'
#' # process Mascot PSMs
#' normPSM(
#'   rptr_intco = 3000,
#'   rm_craps = FALSE,
#'   rm_krts = FALSE,
#'   rm_outliers = FALSE,
#'   plot_violins = TRUE, 
#'   
#'   # `dot-dot-dot` filtration of Mascot PSMs
#'   pep_rank = expr(pep_rank == 1), 
#'   pep_miss = expr(pep_miss <= 2), 
#'   pep_exp_z = expr(pep_exp_z == 2 | pep_exp_z == 3), 
#' )
#'
#'
#' # MaxQuant example
#' dat_dir <- c("C:\\The\\MQ\\PSM_Example")
#'
#' # copy MaxQuant PSM data
#' library(proteoQDB)
#' cptac_mq_psm_1(dat_dir)
#'
#' # copy "expt_smry.xlsx" and "frac_smry.xlsx"
#' cptac_expt_1(dat_dir)
#' cptac_frac_1(dat_dir)
#'
#' # load experiments
#' library(proteoQ)
#' load_expts()
#'
#' # process MaxQuant PSMs
#' normPSM(
#'   fasta = c("~\\proteoQ\\dbs\\refseq\\refseq_hs_2013_07.fasta",
#'             "~\\proteoQ\\dbs\\refseq\\refseq_mm_2013_07.fasta"),
#'   rptr_intco = 3000, 
#'   
#'   # `dot-dot-dot` filtration of MaxQuant PSMs
#'   Charge = expr(Charge == "2" | Charge == "3"), 
#'   `Oxidation (M)` = expr(`Oxidation (M)` == 1),
#' )
#' }
#'
#'@import rlang dplyr purrr ggplot2 RColorBrewer
#'@importFrom stringr str_split
#'@importFrom magrittr %>%
#'@export
normPSM <- function(dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx", 
                    fasta = NULL, pep_unique_by = "group", corrected_int = FALSE, rm_reverses = TRUE, 
                    rptr_intco = 1000, rm_craps = FALSE, rm_krts = FALSE,
                    rm_outliers = FALSE, plot_violins = TRUE, ...) {
  
  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to `dat_dir` first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
  }
  
  dir.create(file.path(dat_dir, "PSM\\Violin"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM\\cache"), recursive = TRUE, showWarnings = FALSE)

  if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^F[0-9]{6}\\.csv$"))) {
    type <- "mascot"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$"))) {
    type <- "mq"
  } else {
    stop("Unknow data type or missing data files.", call. = FALSE)
	}

  pep_unique_by <- rlang::as_string(rlang::enexpr(pep_unique_by))
  
  mget(names(formals()), rlang::current_env()) %>% save_call("normPSM")
  
  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))
  reload_expts()
  
  if (type == "mascot") {
    rmPSMHeaders()
    splitPSM(rptr_intco, rm_craps, plot_violins, ...)
    cleanupPSM(rm_outliers)
    annotPSM(expt_smry, rm_krts, plot_violins)
  } else if (type == "mq") {
    splitPSM_mq(fasta, pep_unique_by, corrected_int, 
                rptr_intco, rm_craps, rm_reverses, plot_violins, ...)
    cleanupPSM(rm_outliers)
    annotPSM_mq(expt_smry, rm_krts, plot_violins)
  }
}


