#'Prepares data for analysis
#'
#'\code{prepDM} prepares a minimal adequate data frame for subsequent analysis.
#'
#'@inheritParams  proteoEucDist
#'@inheritParams  info_anal
#'@param sub_grp Numeric.  A list of sample IDs that will be used in subsequent
#'  analysis.
#'@return A data frame tailored for subsequent analysis.
#'
#' @examples
#' tempData <- prepDM(df, entrez, scale_log2r, sub_grp = label_scheme_sub$Sample_ID)
#'
#' \dontrun{
#' }
#'@import dplyr
#'@importFrom magrittr %>%
prepDM <- function(df, id, scale_log2r, sub_grp, type = "ratio", anal_type) {
  stopifnot(nrow(df) > 0)
  
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  id <- rlang::as_string(rlang::enexpr(id))
  if(anal_type %in% c("ESGAGE", "GSVA")) id <- "entrez"
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  # data filtration dominated by log2R, not Intensity
  df <- df %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, grep(NorZ_ratios, names(.))])) > 0) %>%
    reorderCols(endColIndex = grep("I[0-9]{3}|log2_R[0-9]{3}", names(.)), col_to_rn = id)
  
  Levels <- sub_grp %>%
    as.character(.) %>%
    .[!grepl("^Empty\\.[0-9]+", .)]
  
  dfR <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% sub_grp)) %>%
    dplyr::select(which(not_all_zero(.))) %>% # reference will drop with single reference
    dplyr::select(Levels[Levels %in% names(.)]) # ensure the same order
  
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
#' \code{names_pos<-} rename the columns at the indeces of \code{pos}.
#'
#' @param x A data frame.
#' @param pos Numeric.  The index of coloumns for name change.
#' @param value Characters.  The new column names.
#' @return The data frame with new names.
#'
#' @import dplyr
#' @importFrom magrittr %>%
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
#' @param n_TMT_sets the number of multiplex TMT experiments.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
reorder_files <- function(filelist, n_TMT_sets) {
  newlist <- NULL
  for (i in seq_len(n_TMT_sets))
    newlist <- c(newlist, filelist[grep(paste0("set*.",i), filelist, ignore.case = TRUE)])
  return(newlist)
}


#' Re-order columns in a data frame
#'
#' \code{reorderCols} re-orders columns in a data frame.
#'
#' @param df A data frame.
#' @param endColIndex the indeces of columns to be moved to the end of
#'   \code{df}.
#' @param col_to_rn the column identifier where the values under that column
#'   will be used as row names.
#' @import dplyr rlang
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
reorderCols <- function (df, endColIndex, col_to_rn) {
	if (length(endColIndex) == 0) endColIndex <- grep("I[0-9]{3}|log2_R[0-9]{3}", names(df))
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


#' Replace zero intensity with NA
#'
#' \code{na_zeroIntensity} replaces zero intensity with NA to avoid -Inf in
#' log10 transformation.
#'
#' @param df A list of file names.
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
na_zeroIntensity <- function (df) {
	ind <- grep("I1[0-9]{2}", names(df))

	temp <- df[, ind]
	temp[temp == 0] <- NA
	df[, ind] <- temp

	return(df)
}

#' Summarises numeric values
#'
#' \code{aggrNums} summarises \code{log2FC} and \code{intensity} by the
#' descriptive statistics of \code{c("mean", "median", "weighted.mean",
#' "top.3")}
#'
#' @param f A function for data summarisation.
#' @examples
#' df_num <- aggrNums(median)(df, prot_acc, na.rm = TRUE)
#'
#' @import dplyr rlang
#' @importFrom magrittr %>%
aggrNums <- function(f) {
	function (df, id, ...) {
		id <- rlang::as_string(rlang::enexpr(id))
		dots <- rlang::enexprs(...)

		df %>%
			dplyr::select(id, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise_all(~ f(., !!!dots))
	}
}


#' Calculates weighted mean
#'
#' \code{TMT_wt_mean} calculates the weighted mean of \code{log2FC} and
#' \code{intensity}.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param id The variable to summarise \code{log2FC}.
#' @import dplyr rlang
#' @importFrom stringr str_length
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
TMT_wt_mean <- function (x, id, ...) {
	id <- rlang::as_string(rlang::enexpr(id))
	dots <- rlang::enexprs(...)

	load(file = file.path(dat_dir, "label_scheme.Rdata"))
	TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()

	x_R <- x %>%
			dplyr::select(id, grep("^log2_R[0-9]{3}", names(.))) %>%
			tidyr::gather(key = variable, value = log2_R, -id)

	x_I <- x %>%
			dplyr::select(id, grep("^I[0-9]{3}", names(.))) %>%
			tidyr::gather(key = variable, value = Intensity, -id)

	x_wt <- cbind.data.frame(x_I, log2_R = x_R[, c("log2_R")]) %>%
			dplyr::mutate(variable = gsub("^I", "log2_R", variable)) %>%
			dplyr::mutate(TMT = gsub("^log2_R(.*)\\s+.*", "TMT-\\1", variable)) %>%
			dplyr::mutate(TMT = factor(TMT, levels = TMT_levels)) %>%
			dplyr::mutate(variable = factor(variable, levels = unique(variable))) %>%
			dplyr::mutate(n = stringr::str_length(.[[id]])) %>%
			dplyr::filter(n != 0) %>%
			dplyr::group_by(!!rlang::sym(id), variable) %>%
			dplyr::summarise(log2_R = weighted.mean(log2_R, log10(Intensity), !!!dots)) %>%
			tidyr::spread(variable, log2_R)

	x <- x %>%
			dplyr::select(id, grep("(N_log2_R)|(Z_log2_R)|(I)[0-9]{3}", names(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			dplyr::summarise_all(funs(mean(., !!!dots))) %>%
			dplyr::left_join(x_wt, by = id)

	rm(x_R, x_I, x_wt)

	cbind.data.frame(
		x[, names(x) == id],
		x[, grepl("^I[0-9]{3}", names(x))],
		x[, grepl("^N_I[0-9]{3}", names(x))],
		x[, grepl("^log2_R[0-9]{3}", names(x))],
		x[, grepl("^N_log2_R[0-9]{3}", names(x))],
		x[, grepl("^Z_log2_R[0-9]{3}", names(x))])
}


#' Calculate top.3
#'
#' \code{TMT_wt_mean} calculates the weighted mean of \code{log2FC} and
#' \code{intensity}.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @param id The variable to summarise \code{log2FC}.
#' @examples
#' df_num <- TMT_top_n(df, prot_acc, na.rm = TRUE)
#'
#' @import dplyr rlang
#' @importFrom magrittr %>%
TMT_top_n <- function (x, id, ...) {

	id_var <- rlang::enexpr(id)
	id_nm <- rlang::expr_name(id_var)

	dots <- rlang::enexprs(...)
	args <- lapply(dots, expr_name)

	x %>%
		dplyr::select(id_nm, grep("log2_R[0-9]{3}|I[0-9]{3}", names(.))) %>%
		dplyr::mutate(sum_Intensity =rowSums(.[grep("^I1[0-9]{2}", names(.))], na.rm = TRUE)) %>%
		dplyr::group_by(!!id_var) %>%
		dplyr::top_n(n = 3, wt = sum_Intensity) %>%
		dplyr::select(-sum_Intensity) %>%
		dplyr::summarise_all(funs(median(., !!!args)))
}


#' Finds all-zero column(s)
#'
#' \code{not_all_zero} identifies the column indeces with all NA values.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
not_all_zero <- function (x) (colSums(x != 0, na.rm = TRUE) > 0)


#' Finds all-NA column(s)
#'
#' \code{not_all_NA} identifies the column indeces with all NA values.
#'
#' @param x A data frame of \code{log2FC} and \code{intensity}.
#' @import dplyr rlang
#' @importFrom magrittr %>%
not_all_NA <- function (x) (colSums(!is.na(x), na.rm = TRUE) > 0)


#' Sets up the column annotation in heat maps
#'
#' @import dplyr rlang
#' @importFrom magrittr %>%
colAnnot <- function (annot_cols = NULL, sample_ids, annot_colnames = NULL) {
	if(is.null(annot_cols)) return(NA)

  load(file = file.path(dat_dir, "label_scheme.Rdata"))
	exists <- annot_cols %in% names(label_scheme)

	if(sum(!exists) > 0) {
		warning(paste0("Column '", annot_cols[!exists], "'",
		               " not found in \'label_scheme\' and will be skipped."))
		annot_cols <- annot_cols[exists]
	}

	if(length(annot_cols) == 0) return(NA)

	x <- label_scheme %>%
		dplyr::filter(Sample_ID %in% sample_ids) %>%
		dplyr::select(annot_cols, Sample_ID) %>%
		dplyr::select(which(not_all_NA(.))) %>%
		dplyr::filter(!grepl("^Empty\\.", .[["Sample_ID"]]),
		              !is.na(.[["Sample_ID"]])) %>%
		data.frame(check.names = FALSE) %>%
		`rownames<-`(.[["Sample_ID"]])

	if(any(duplicated(x[["Sample_ID"]]))) stop("Duplicated sample IDs found\n")

	if(!"Sample_ID" %in% annot_cols) x <- x %>% dplyr::select(-Sample_ID)

	if(ncol(x) == 0) return(NA)

	if("TMT_Set" %in% names(x)) {
		x <- x %>%
			tibble::rownames_to_column() %>%
			mutate(TMT_Set, as.factor(TMT_Set)) %>%
			tibble::column_to_rownames(var = "rowname")
	}

	return(x)
}


#' Customizes the colors in column annotation
#'
#' @import dplyr rlang
#' @importFrom magrittr %>%
setHMColor <- function (annotation_col) {
	ncol <- ncol(annotation_col)

	if(is.null(ncol)) return (NA)
	if(ncol == 0) return (NA)

	annotation_col <- annotation_col %>%
		dplyr::mutate_at(vars(one_of("TMT_Set")), ~ factor(TMT_Set)) %>%
		dplyr::mutate_if(is.character, as.factor)

	palette <- lapply(names(annotation_col), function(x) {
		n <- nlevels(annotation_col[[x]])

		palette <- if(n <= 9 & n >= 3) {
			# brewer.pal(n, name = "Set1")
		  brewer.pal(n, name = "Set2")
		} else if(n > 9) {
			colorRampPalette(brewer.pal(n = 7, "Set1"))(n)
		} else if(n == 2) {
			# c("#E41A1C", "#377EB8")
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
#' @examples
#' setHMlims(df, -1, 1)
#'
#' @import dplyr rlang
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
setHMlims <- function (x, xmin, xmax) {
	x[x < xmin] <- xmin
	x[x > xmax] <- xmax

	x
}


#' Calculate the log2-ratio to the control group of samples
#'
#' @examples
#' ratio_toCtrl(df, "gene", label_scheme_sub, Heatmap_Group)
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
#'
#'@inheritParams proteoHist
#'@param ... Parameters for \code{\link[mice]{mice}}
#'@return \code{Peptide_impNA.txt} for peptide data and \code{Protein_impNA.txt}
#'  for protein data.
#'@import dplyr purrr rlang mice
#'@importFrom magrittr %>%
#'@export
imputeNA <- function (id, ...) {

	my_mice <- function (data, ...) {
		data <- rlang::enexpr(data)
		dots <- rlang::enexprs(...)

		mice_call <- rlang::expr(mice(data = !!data, !!!dots))
		rlang::expr_print(mice_call)
		rlang::eval_bare(mice_call, env = caller_env())
	}

	handleNA <- function (x, ...) {
		ind <- purrr::map(x, is.numeric) %>% unlist()

		if(sum(ind) < ncol(x)) message("Not all of the columns are numeric.
		                               Only the numeric columns are taken for NA imputation.")

		# handle special column names
		nm_orig <- names(x[, ind])

		x[, ind] <- x[, ind] %>%
			data.frame(check.names = TRUE) %>%
			my_mice(...) %>%
			mice::complete(1) %>%
			`colnames<-`(nm_orig)

		return(x)
	}

	id <- rlang::as_string(rlang::enexpr(id))

	cat(paste0("id = \"", id, "\"", " by the current call", "\n"))
	id <- match_identifier(id)
	cat(paste0("id = \"", id, "\"", " after parameter matching to normPep() or normPrn()", "\n"))

	if (id %in% c("pep_seq", "pep_seq_mod")) {
		src_path <- file.path(dat_dir, "Peptide", "Peptide.txt")
		filename <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
	} else if (id %in% c("prot_acc", "gene")) {
		src_path <- file.path(dat_dir, "Protein", "Protein.txt")
		filename <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
	}

	if(!file.exists(filename)) {
		df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
		                        comment.char = "#"), error = function(e) NA)

		if(!is.null(dim(df))) {
			message(paste("File loaded:", gsub("\\\\", "/", src_path)))
		} else {
			stop(paste("No such file or directory:", gsub("\\\\", "/", src_path)))
		}

		df[, grep("N_log2_R", names(df))]  <- df[, grep("N_log2_R", names(df))]  %>% handleNA(...)

		fn_params <- file.path(dat_dir, "Protein\\Histogram", "MGKernel_params_N.txt")
		if(file.exists(fn_params)) {
			# back calculate Z_log2_R from N_log2_R
			cf_SD <-
				read.csv(fn_params,
								check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>%
				dplyr::filter(!duplicated(Sample_ID)) %>%
				dplyr::select(Sample_ID, fct)

			df[, grep("Z_log2_R", names(df))] <-
				mapply(normSD,
				       df[, grepl("^N_log2_R[0-9]{3}", names(df))],
				       center = 0, SD = cf_SD$fct, SIMPLIFY = FALSE) %>%
				data.frame(check.names = FALSE) %>%
				`colnames<-`(gsub("N_log2", "Z_log2", names(.)))
		} else {
			df[, grep("Z_log2_R", names(df))]  <- df[, grep("Z_log2_R", names(df))] %>% handleNA(...)
		}

		if(any(duplicated(df[[id]]))) {
			if(id == "pep_seq") {
				cat("\`pep_seq\` is not uqique for rownames; use \`pep_seq_mod\` instead.\n")
				rownames(df) <- df[["pep_seq_mod"]]
			}
			if(id == "gene") {
				cat("\`gene\` is not uqique for rownames; use \`prot_acc\` instead.\n")
				rownames(df) <- df[["prot_acc"]]
			}
		} else {
			rownames(df) <- df[[id]]
		}

		write.table(df, filename, sep = "\t", col.names = TRUE, row.names = FALSE)
	} else {
		cat("NA imputation has already been performed!\n")
	}
}


#'Impute NA values for peptide data
#'
#'\code{pepImp} is a wrapper of \code{\link{imputeNA}} for peptide data
#'
#'@rdname imputeNA
#'
#' @examples
#' \dontrun{
#' pepImp(
#'   m = 3,
#'   maxit = 3,
#' )
#' }
#'
#'@export
pepImp <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  imputeNA(id = pep_seq, ...)
}


#'Impute NA values for protein data
#'
#'\code{prnImp} is a wrapper of \code{\link{imputeNA}} for protein data
#'
#'@rdname imputeNA
#'
#' @examples
#' \dontrun{
#' prnImp(
#'   m = 5,
#'   maxit = 5,
#' )
#' }
#'
#'@export
prnImp <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  imputeNA(id = gene, ...)
}


#' Adds protein annotation
#'
#' \code{annotPrn} cross-referencing proteins among \code{uniprot_acc},
#' \code{uniprot_id}, \code{refseq} and \code{entrez}.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>% %$%
annotPrn <- function (df, acc_type) {
	stopifnot ("prot_acc" %in% names(df))
  
	acc_type <- find_acctype() %>% tolower()

	if(acc_type == "refseq_acc") {
		key <- "refseq_acc"
	} else if(acc_type == "uniprot_id") {
		key <- "uniprot_id"
	} else if (acc_type == "uniprot_acc") {
	  key <- "uniprot_acc"
	} else {
		stop("Unrecognized protein accesion type; need to one of \'uniprot_id\', \'uniprot_acc\' or \'refseq_acc\'",
		     call. = FALSE)
	}

	lookup <- dbs$prn_annot %>%
		dplyr::select(-status, -organism, -length) %>%
		dplyr::filter(!duplicated(.[[key]])) # %>% 
	  # dplyr::mutate(!!key := as.character(.[[key]]))

	if(any(names(df) == "prot_desc")) lookup <- lookup %>% dplyr::select(-prot_desc)

	df <- df %>% 
	  # dplyr::mutate(prot_acc = as.character(prot_acc)) %>% 
		dplyr::left_join(lookup, by = c("prot_acc" = key)) 
	
	ind <- is.na(df$gene) | str_length(df$gen) == 0
	if (acc_type %in% c("uniprot_id", "uniprot_acc") & sum(ind) > 0) {
	  substi_gns <- df[ind, "prot_desc"] %>% 
	    as.character(.) %>% 
	    strsplit(., ' GN=', fixed = TRUE) %>% 
	    plyr::ldply(., rbind) 
	  
	  if (ncol(substi_gns) > 1) {
	    substi_gns <- substi_gns %$% gsub("\\s+.*", "", .[, 2])
	    df[ind, "gene"] <- substi_gns
	  }

		rm(substi_gns)
	}
	
	# annotate NA genes with prot_acc
	ind <- is.na(df$gene) | str_length(df$gen) == 0
	if(sum(ind) > 0) df[ind, "gene"] <- df[ind, "prot_acc"]

	return(df)
}


#' Adds kinase annotation
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
annotKin <- function (df, acc_type) {
	stopifnot ("prot_acc" %in% names(df))
	
	# acc_type <- tolower(acc_type)
	acc_type <- find_acctype() %>% tolower()

	lookup <- kinase_lookup %>%
		dplyr::select(refseq_acc, gene, kin_attr, kin_class, kin_order) %>%
		dplyr::rename(gene_lookup = gene) %>%
		dplyr::filter(!duplicated(.[["refseq_acc"]]))

	if (acc_type == "refseq_acc") {
		key <- "prot_acc"
		df <- df %>% dplyr::left_join(lookup, by = c("prot_acc" = "refseq_acc"))
	} else if (acc_type %in% c("uniprot_id", "uniprot_acc")) {
		key <- "refseq_acc"
		if (!key %in% names(df)) stop("Column key, 'refseq_acc', is missing for kinase annotation.", 
		                              call. = FALSE)
		df <- df %>% dplyr::left_join(lookup, by = "refseq_acc")
	}

	# higher priority for genes from the lookup table
	kinases <- df[[key]] %in% lookup[["refseq_acc"]]
	if(sum(kinases) > 0) df$gene[kinases] <- as.character(df$gene_lookup[kinases])

	df <- df %>%
		dplyr::mutate(gene = as.factor(gene)) %>%
		dplyr::select(-gene_lookup)

	df$kin_attr[is.na(df$kin_attr)] <- FALSE

	return(df)
}


#' Saves the arguments in a function call
#'
#' @param call_pars Language.
#' @param fn The name of function being called.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
save_call <- function(call_pars, fn) {
	dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)
	save(fn, file = file.path(dat_dir, "Calls", paste0(fn, "_call.Rdata")))

	call_pars[names(call_pars) == "..."] <- NULL

	purrr::map(call_pars , as.character) %>%
		do.call('rbind', .) %>%
		data.frame() %>%
		`colnames<-`(paste("value", 1:ncol(.), sep = ".")) %>%
		tibble::rownames_to_column("var") %>%
		write.table(., file.path(dat_dir, "Calls", paste0(fn, ".txt")), sep = '\t',
		            col.names = TRUE, row.names = FALSE)
}


#' Matches the current id to the id in normPep or normPrn
#'
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_identifier <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene")) {
	if (id %in% c("prot_acc", "gene")) {
		fn_pars <- "normPrn.txt"
	} else if (id %in% c("pep_seq", "pep_seq_mod")) {
		fn_pars <- "normPep.txt"
	}

	call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
	                               header = TRUE, sep = "\t", comment.char = "#"),
	                      error = function(e) NA)

	if(!is.null(dim(call_pars))) {
		id <- call_pars %>%
			dplyr::filter(var == "id") %>%
			dplyr::select("value.1") %>%
			unlist() %>%
			as.character()
	}

	return(id)
}


#' Matches the file name containing the summary of TMT experiment
#'
#' The default file name is \code{expt_smry.xlsx}. The \code{match_expt} matches
#' the file name provided by users.
#'
#' @param fn_pars.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_expt <- function (fn_pars = "load_expts.txt") {
	call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
	                               header = TRUE, sep = "\t", comment.char = "#"),
	                      error = function(e) NA)

	if(!is.null(dim(call_pars))) {
		expt_smry <- call_pars %>%
			dplyr::filter(var == "expt_smry") %>%
			dplyr::select("value.1") %>%
			unlist() %>%
			as.character()
	} else {
		expt_smry <- "expt_smry.xlsx"
	}

	return(expt_smry)
}


#' Matches the file name containing the information of analyte fractionation
#'
#' The default file name is \code{frac_smry.xlsx}. The \code{match_frac} matches
#' the file name provided by users.
#'
#' @param fn_pars.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_frac <- function (fn_pars = "load_expts.txt") {
	call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
	                               header = TRUE, sep = "\t", comment.char = "#"),
	                      error = function(e) NA)

	if(!is.null(dim(call_pars))) {
		frac_smry <- call_pars %>%
			dplyr::filter(var == "frac_smry") %>%
			dplyr::select("value.1") %>%
			unlist() %>%
			as.character()
	} else {
		frac_smry <- "frac_smry.xlsx"
	}

	return(frac_smry)
}


#' Reload the "expt_smry.xlsx" and "frac_smry.xlsx"
#'
#' @import rlang
#' @importFrom magrittr %>%
#' @importFrom fs file_info
reload_expts <- function() {
  expt_smry <- match_expt()
  frac_smry <- match_frac()
  
  fi_xlsx <- fs::file_info(file.path(dat_dir, expt_smry))$change_time
  
  if(is.na(fi_xlsx)) stop("Time stamp of `expt_smry.xlsx` not available.")
  
  fi_rda <- fs::file_info(file.path(dat_dir, "label_scheme.Rdata"))$change_time
  if(fi_xlsx > fi_rda) {
    load_expts(dat_dir = dat_dir, expt_smry = !!expt_smry, frac_smry = !!frac_smry)
  }
}


#' Replaces NA genes
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
replace_na_genes <- function(df, acc_type) {
	acc_type <- tolower(acc_type)

	if(acc_type == "refseq_acc") {
		na_gene <- (is.na(df[, c("gene")])) | (str_length(df$gene) == 0)
		df$gene <- as.character(df$gene)
		df$gene[na_gene] <- as.character(df$prot_acc[na_gene])
	} else if(acc_type == "uniprot_id") {
		temp <- data.frame(do.call('rbind', strsplit(as.character(df$prot_desc), 'GN=', fixed = TRUE))) %>%
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


#' Adds protein description
#'
#' \code{annotPrndesc} annotates protein descriptions based on the \code{fasta}
#'
#' @import plyr dplyr purrr rlang stringr seqinr
#' @importFrom magrittr %>% %$%
annotPrndesc <- function (df, fasta){
  stopifnot("prot_acc" %in% names(df))
  
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  acc_type <- find_acctype() %>% tolower()
  
  if (acc_type == "refseq_acc") {
    key <- "refseq_acc"
  } else if (acc_type %in% "uniprot_id") {
    key <- "uniprot_id"
  } else if (acc_type %in% "uniprot_acc") {
    key <- "uniprot_acc"
  } else {
    warning("Unkown accession type.")
  }
  
  df <- df %>% 
    dplyr::mutate(prot_acc = gsub("^.*\\|(.*)\\|.*$", "\\1", prot_acc))
  
  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- purrr::map(fasta, ~ {
        seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
      }) %>% do.call(`c`, .) %>% 
        `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
        .[names(.) %in% unique(df$prot_acc)]
      
      if (length(fasta) == 0) {
        stop("No fasta entries match protein accessions; probably wrong fasta file.", 
             call. = FALSE)
      }
      
      lookup <- dplyr::bind_cols(
        # prot_acc = seqinr::getName(fasta), 
        prot_acc = names(fasta), 
        prot_desc = seqinr::getAnnot(fasta) %>% 
          purrr::map(., `[[`, 1) %>% 
          unlist(), 
        prot_mass = purrr::map_dbl(fasta, ~ {seqinr::getSequence(.x) %>% seqinr::pmw()})
      ) %>% 
        dplyr::filter(.$prot_acc %in% unique(df$prot_acc))
    } else {
      stop("Wrong FASTA file path or name.", call. = FALSE)
    }
  } else {
    stop("FASTA file not provided.")
  }
  
  rm(fasta)
  
  df %>% 
    dplyr::left_join(lookup, by = "prot_acc") %>% 
    dplyr::mutate(prot_mass = round(prot_mass, digits = 1))
}


#' Add peptide start and end positions
#'
#' \code{annotPeppos} annotates the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}.
#'
#' @import dplyr purrr rlang stringr seqinr
#' @importFrom magrittr %>% %$%
annotPeppos <- function (df, fasta){
  stopifnot("prot_acc" %in% names(df))
  stopifnot("pep_seq" %in% names(df))
  
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  acc_type <- find_acctype() %>% tolower()
  
  if (acc_type == "refseq_acc") {
    key <- "refseq_acc"
  } else if (acc_type %in% "uniprot_id") {
    key <- "uniprot_id"
  } else if (acc_type %in% "uniprot_acc") {
    key <- "uniprot_acc"
  } else {
    warning("Unkown accession type.")
  }
  
  df <- df %>% 
    dplyr::mutate(prot_acc = gsub("^.*\\|(.*)\\|.*$", "\\1", prot_acc))
  
  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- purrr::map(fasta, ~ {
        seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
      }) %>% do.call(`c`, .) %>% 
        `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
        .[names(.) %in% unique(df$prot_acc)]
      
      if (length(fasta) == 0) {
        stop("No fasta entries match protein accessions; probably wrong fasta file.", 
             call. = FALSE)
      }
      
      pep_pos_all <- purrr::map2(as.list(df$prot_acc), as.list(df$pep_seq), ~ {
        fasta_sub <- fasta %>% .[names(.) == .x]
        pep_seq <- as.character(.y)
        
        if (!rlang::is_empty(fasta_sub)) {
          pep_pos <- str_locate(fasta_sub, pattern = pep_seq)
          pep_pos <- cbind(pep_seq, pep_pos)
        } else {
          pep_pos <- cbind(pep_seq, start = NA, end = NA)
        }
      }) %>% 
        do.call(rbind, .) %>% 
        `colnames<-`(c("pep_seq", "pep_start", "pep_end")) %>% 
        data.frame(check.names = FALSE)
    } else {
      stop("Wrong FASTA file path or name.", call. = FALSE)
    }
  } else {
    stop("FASTA file not provided.")
  }
  
  rm(fasta)
  
  df %>% dplyr::left_join(pep_pos_all, by = "pep_seq")
}


#' Subset fasta by accession type
#'
#' @import plyr dplyr purrr rlang seqinr
#' @importFrom magrittr %>%
subset_fasta <- function (df, fasta, acc_type) {
  stopifnot("prot_acc" %in% names(df))
  
  if (acc_type == "refseq_acc") {
    key <- "refseq_acc"
  } else if (acc_type == "uniprot_id") {
    key <- "uniprot_id"
  } else if (acc_type == "uniprot_acc") {
    key <- "uniprot_acc"
  } else {
    warning("Unkown accession type.")
  }
  
  fasta <- purrr::map(fasta, ~ {
    seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
  }) %>% do.call(`c`, .)
  
  if (key == "uniprot_id") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|.*\\|(.*)$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } else if (key == "uniprot_acc") {
    fasta <- fasta %>% 
      `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
      .[names(.) %in% unique(df$prot_acc)]
  } else if (key == "refseq_acc") {
    fasta <- fasta %>% 
      .[names(.) %in% unique(df$prot_acc)]
  }    
}


#' Calculates protein percent coverage
#'
#' @import plyr dplyr purrr rlang seqinr
#' @importFrom magrittr %>%
calc_cover <- function(df, id, fasta = NULL) {
  id <- rlang::as_string(rlang::enexpr(id))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  acc_type <- find_acctype() %>% tolower()
  
  if (acc_type == "refseq_acc") {
    key <- "refseq_acc"
  } else if (acc_type %in% "uniprot_id") {
    key <- "uniprot_id"
  } else if (acc_type %in% "uniprot_acc") {
    key <- "uniprot_acc"
  } else {
    warning("Unkown accession type.")
  }

  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- subset_fasta(df, fasta, acc_type)
      
      if (length(fasta) == 0) {
        stop("No fasta entries matched the type of protein accession. Check the correctness of fasta file.", 
             call. = FALSE)
      }
      
      if (length(fasta) <= 200) {
        warning("Less than 200 entries in fasta match by protein accession. 
                Make sure the fasta file is correct.")
      }
      
      lookup <- data.frame(prot_acc = names(fasta), length = getLength(fasta)) %>%
        dplyr::rename(!!key := prot_acc) %>%
        dplyr::mutate(!!key := gsub(".*\\|", "", !!rlang::sym(key)))
    } else {
      cat(fasta, "\n")
      stop("Not all fasta files were found.", call. = FALSE)
    }
  } else {
    warning("Use pre-computed database to calculate protein coverages.")
    lookup <- dbs$prn_annot %>%
      dplyr::select(key, length) %>%
      dplyr::filter(!is.na(.[[key]]), !duplicated(.[[key]]))
  }
  
  df <- df %>%
    dplyr::select(prot_acc, pep_start, pep_end) %>%
    dplyr::left_join(lookup, by = c("prot_acc" = key)) %>%
    dplyr::filter(!is.na(length))
  
  if (nrow(df) == 0) stop("Probably incorrect accession types in the fasta file.", call. = FALSE)
  
  df <- df %>%
    dplyr::filter(.[["pep_start"]] <= .[["length"]]) %>%
    dplyr::filter(.[["pep_end"]] <= .[["length"]]) %>%
    split(.[["prot_acc"]], drop = TRUE) %>%
    purrr::map(function (x) {
      len <- x[1, "length"]
      aa_map <- rep(NA, len)
      for (i in 1:nrow(x)) aa_map[x[i, ]$pep_start:x[i, ]$pep_end] <- TRUE
      sum(aa_map, na.rm = TRUE)/len
    } ) %>%
    do.call("rbind", .) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`("prot_cover") %>%
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::mutate(prot_cover = ifelse(prot_cover > 1, 1, prot_cover)) %>%
    annotPrn(acc_type) %>%
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::select(!!rlang::sym(id), prot_cover) %>%
    dplyr::summarise_all(~max(., na.rm = TRUE)) %>%
    dplyr::mutate(prot_cover = round(prot_cover * 100, digits = 1)) %>%
    dplyr::mutate(prot_cover = paste0(prot_cover, "%"))
  
  return(df)
}


#' Matches formulas to those in calls to pepSig or prnSig
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_fmls <- function(formulas) {
  fml_file <-  file.path(dat_dir, "Calls\\prnSig_formulas.Rdata")

  if(file.exists(fml_file)) {
    load(file = fml_file)
  } else {
    stop("Run `prnSig()` first.")
  }
  
	fml_chr <- formulas %>%
		as.character() %>%
		gsub("\\s+", "", .)

	prnSig_chr <- prnSig_formulas %>%
		purrr::map(~ .[is_call(.)]) %>%
		as.character() %>%
		gsub("\\s+", "", .)

	ok <- purrr::map_lgl(fml_chr, ~ . %in% prnSig_chr)

	if(!all(ok))
		stop("Formula match failed: ", formulas[[which(!ok)]],
		     " not found in the latest call to 'prnSig(...)'.")
}


#' Converts log2FC to linear fold changes
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
to_linfc <- function(df) {
		nms <- rownames(df)

		df %>%
			purrr::map(~ {ifelse(.x > 0, 2^.x, -1/(2^.x))}) %>%
			data.frame(check.names = FALSE) %>%
			`rownames<-`(nms)
}


#' Extract RAW MS file names
#'
#' Extract a list of \code{RAW} file names that can be passed to \code{frac_smry.xlsx}
#'
#' @examples
#' \dontrun{
#' # Supposed that RAW MS files are stored under "C:\my_raw"
#' extract_raws("C:\\my_raw")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
#' @importFrom tools md5sum
extract_raws <- function(raw_dir) {
  fns <- names(tools::md5sum(dir(raw_dir, pattern = "\\.raw$", full.names = FALSE)))
  data.frame(Fraction = seq_along(fns), RAW_File = fns)
}


#' Remove single-value columns
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
rm_sglval_cols <- function (x) {
  sgl_val <- x %>% 
    summarise_all(funs(n_distinct(.))) %>% 
    purrr::map(~ .x == 1) %>% 
    purrr::flatten_lgl()
  
  x[, !sgl_val, drop = FALSE]
}


#' Combine data with metadata
#'
#' @import dplyr purrr
#' @importFrom magrittr %>%
cmbn_meta <- function(data, metadata) {
  data %>% 
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::left_join(metadata) %>%
    dplyr::mutate_at(vars(one_of("Color", "Fill", "Shape", "Size", "Alpha")), ~ as.factor(.)) %>%
    dplyr::select(which(not_all_NA(.))) %>% 
    rm_sglval_cols()
}


#' Check file names for ggsave()
#'
gg_imgname <- function(filename) {
  fn_prx <- gsub("\\..*$", "", filename)
  fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
  
  exts <- c("png", "eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg") 
  
  if(! fn_suffix %in% exts) {
    warning(paste0("Unrecognized file extenstion: '", fn_suffix, "'. Image will be saved as a 'png'.\n"))
    fn_suffix <- "png"
  }
  
  paste0(fn_prx, ".", fn_suffix)
}


#' Check file names for ggsave()
#' 
#' @import dplyr purrr
#' @importFrom magrittr %>%
rm_pval_whitespace <- function(df) {
  df <- df %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.character) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ gsub("\\s*", "", .x) ) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.numeric)
}


#' Match the database of gene sets
#' 
#' @import dplyr purrr
#' @importFrom magrittr %>%
match_gsets <- function(gset_nm = "go_sets") {
  stopifnot(all(gset_nm %in% c("go_sets", "kegg_sets", "c2_msig")))
  
  gsets <- dbs %>% .[names(.) %in% gset_nm]
  stopifnot(length(gsets) > 0)
  
  is_null <- purrr::map_lgl(gsets, ~ is.null(.x))
  gset_nm <- gset_nm[!is_null]
  
  
  purrr::walk2(is_null, names(is_null), 
               ~ if(.x) warning("Gene set: `", .y, "` not found", call. = FALSE))
  
  gsets[is_null] <- NULL
  gsets <- gsets %>% purrr::reduce(`c`)
}


#' Matches the current gene set names to those used in prnGSPA()
#'
#' @param id. 
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_gset_nm <- function (id = c("go_sets", "kegg_sets", "c2_msig")) {
  fn_pars <- "prnGSPA.txt"

  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if(!is.null(dim(call_pars))) {
    id <- call_pars %>%
      dplyr::filter(var == "gset_nm") %>% 
      dplyr::select(-1) %>% 
      unlist() %>%
      as.character()
  }
  
  return(id)
}



