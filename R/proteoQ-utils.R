#'Prepares data for analysis
#'
#'\code{prepDM} prepares a minimal adequate data frame for subsequent analysis.
#'
#'@inheritParams  prnEucDist
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
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  id <- rlang::as_string(rlang::enexpr(id))
  if (anal_type %in% c("ESGAGE", "GSVA")) id <- "entrez"
  if ((anal_type %in% c("GSEA")) & (id != "gene")) stop("Primary ID is not `gene`.")
  
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

	load(file = file.path(dat_dir, "label_scheme.rda"))
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


#' Finds all-NaN column(s)
not_all_nan <- function(x, ...) {
  sum(is.nan(x), ...) != length(x)
}


#' Finds all-NaN column(s)
is_all_nan <- function(x, ...) {
  sum(is.nan(x), ...) == length(x)
}


#' Sets up the column annotation in heat maps
#'
#' @import dplyr rlang
#' @importFrom magrittr %>%
colAnnot <- function (annot_cols = NULL, sample_ids = NULL, annot_colnames = NULL) {
	if (is.null(annot_cols)) return(NA)

  load(file = file.path(dat_dir, "label_scheme.rda"))
	exists <- annot_cols %in% names(label_scheme)

	if (sum(!exists) > 0) {
		warning(paste0("Column '", annot_cols[!exists], "'",
		               " not found in \'label_scheme\' and will be skipped."))
		annot_cols <- annot_cols[exists]
	}

	if (length(annot_cols) == 0) return(NA)

	x <- label_scheme %>%
		dplyr::filter(Sample_ID %in% sample_ids) %>%
		dplyr::select(annot_cols, Sample_ID) 
	
	idx <- !not_all_NA(x)
	if (sum(idx) > 0) stop("No aesthetics defined under column `", names(x)[idx], "`")

	x <- x %>%
		dplyr::filter(!grepl("^Empty\\.", .[["Sample_ID"]]),
		              !is.na(.[["Sample_ID"]])) %>%
		data.frame(check.names = FALSE) %>%
		`rownames<-`(.[["Sample_ID"]])

	if (any(duplicated(x[["Sample_ID"]]))) stop("Duplicated sample IDs found\n")

	if (!"Sample_ID" %in% annot_cols) x <- x %>% dplyr::select(-Sample_ID)

	if (ncol(x) == 0) return(NA)

	if ("TMT_Set" %in% names(x)) {
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
#'Users should avoid call the method directly, but instead use the following
#'wrappers.
#'
#'@param overwrite Logical. If true, overwrite the previous results. The default
#'  is FALSE.
#'@inheritParams prnHist
#'@param ... Parameters for \code{\link[mice]{mice}}
#'@return \code{Peptide_impNA.txt} for peptide data and \code{Protein_impNA.txt}
#'  for protein data.
#'
#'@import dplyr purrr rlang mice
#'@importFrom magrittr %>%
#'@export
imputeNA <- function (id, overwrite = FALSE, ...) {
	my_mice <- function (data, ...) {
		data <- rlang::enexpr(data)
		dots <- rlang::enexprs(...)

		mice_call <- rlang::expr(mice(data = !!data, !!!dots))
		rlang::eval_bare(mice_call, env = caller_env())
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

	id <- rlang::as_string(rlang::enexpr(id))

	if (id %in% c("pep_seq", "pep_seq_mod")) {
		src_path <- file.path(dat_dir, "Peptide", "Peptide.txt")
		filename <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
	} else if (id %in% c("prot_acc", "gene")) {
		src_path <- file.path(dat_dir, "Protein", "Protein.txt")
		filename <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
	}

	if ((!file.exists(filename)) || overwrite) {
		df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE, sep = "\t",
		                        comment.char = "#"), error = function(e) NA)

		if (!is.null(dim(df))) {
			message(paste("File loaded:", gsub("\\\\", "/", src_path)))
		} else {
			stop(paste("File or directory not found:", gsub("\\\\", "/", src_path)))
		}

		df[, grep("N_log2_R", names(df))]  <- df[, grep("N_log2_R", names(df))]  %>% handleNA(...)

		fn_params <- file.path(dat_dir, "Protein\\Histogram", "MGKernel_params_N.txt")
		if (file.exists(fn_params)) {
			cf_SD <- 
				read.csv(fn_params, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>%
				dplyr::filter(!duplicated(Sample_ID)) %>%
				dplyr::select(Sample_ID, fct)

			df[, grep("Z_log2_R", names(df))] <-
				mapply(normSD,
				       df[, grepl("^N_log2_R[0-9]{3}", names(df))], center = 0, SD = cf_SD$fct, SIMPLIFY = FALSE) %>%
				data.frame(check.names = FALSE) %>%
				`colnames<-`(gsub("N_log2", "Z_log2", names(.)))
		} else {
			df[, grep("Z_log2_R", names(df))]  <- df[, grep("Z_log2_R", names(df))] %>% handleNA(...)
		}

		if (any(duplicated(df[[id]]))) {
			if (id == "pep_seq") {
				warning("\`pep_seq\` is not uqique for rownames; use \`pep_seq_mod\` instead.\n", call. = FALSE)
				rownames(df) <- df[["pep_seq_mod"]]
			}
			if (id == "gene") {
				warning("\`gene\` is not uqique for rownames; use \`prot_acc\` instead.\n", call. = FALSE)
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
#' # ===================================
#' # Peptide NA imputation
#' # ===================================
#' pepImp(
#'   m = 3,
#'   maxit = 3,
#' )
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
#' # ===================================
#' # Protein NA imputation
#' # ===================================
#' prnImp(
#'   m = 5,
#'   maxit = 5,
#' )
#'
#'@export
prnImp <- function (...) {
  check_dots(c("id"), ...)
  id <- match_call_arg(normPSM, group_pep_by)
  imputeNA(id = !!id, ...)
}


#'Species lookup
sp_lookup <- function(species) {
  switch(species, 
         human = "hs",
         mouse = "mm",
         rat = "rn",
         fly = "dm", 
         bovine = "bt",
         dog = "cf", 
         crap = "crap", 
         stop("Unknown `species`.", Call. = FALSE)
  )    
}


#'Toxonomy lookup
taxid_lookup <- function(species) {
  switch (species,
          human = 9606, 
          mouse = 10090,
          rat = 10116, 
          fly = 7227, 
          bovine = 9913,
          dog = 9612, 
          crap = 000000, 
          stop("Unknown `species`.", Call. = FALSE)
  )
}


add_entrez <- function (acc_lookup) {
  # 1. find all speices and acc_types
  sp_map <- c(
    human = "hs",
    mouse = "mm",
    rat = "rn",
    fly = "dm", 
    bovine = "bt",
    dog = "cf", 
    crap = "crap"
  )
  
  acc_map <- c(
    uniprot_acc = "uniprot_",
    uniprot_id = "uniprot_", 
    refseq_acc = "refseq_"
  )
  
  acc_type <- unique(acc_lookup$acc_type)
  species <- unique(acc_lookup$species)
  abbr_sp <- sp_map[species]
  abbr_acc <- acc_map[acc_type]
  
  entrez_key <- switch(acc_type, 
    uniprot_acc = "uniprot_acc", 
    uniprot_id = "uniprot_acc", 
    refseq_acc = "refseq_acc", 
    stop("Unknown `accession type`.", Call. = FALSE)
  )

  stopifnot(length(acc_type) == 1)

  # multiple uniprot_acc(s) can share the same entrez id
  filelist <- paste0(abbr_acc, "entrez_", abbr_sp)
  data(package = "proteoQ", list = filelist)
  entrez <- purrr::map(filelist, ~ get(.x)) %>% bind_rows() 
  
  if (acc_type == "uniprot_id") {
    acc_lookup <- left_join(acc_lookup, entrez, by = c("uniprot_acc" = entrez_key))
  } else {
    acc_lookup <- left_join(acc_lookup, entrez, by = c("prot_acc" = entrez_key))
  }
  
  return(acc_lookup)
}


#' Determine the protein accession type from PSM tables
#'
#' Find the protein accession from a non-cRAP entry and parse it.
#'
#' \code{parse_acc} parse the protein accession.
parse_acc <- function(df) {
  stopifnot("prot_acc" %in% names(df))
  
  prn_acc <- df %>%
    dplyr::filter(!grepl("^REV__", prot_acc)) %>% 
    dplyr::filter(!grepl("^CON__", prot_acc)) %>% 
    dplyr::select(prot_acc) %>%
    unlist %>%
    .[1]

  if (grepl("_[A-Z]+", prn_acc)) {
    acc_type <- "uniprot_id"
  } else if (grepl("^NP_[0-9]+", prn_acc)) {
    acc_type <- "refseq_acc"
  } else if (grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", prn_acc)) {
    acc_type <- "uniprot_acc"
  } else {
    # stop("Unknown protein accession", call. = FALSE)
    acc_type <- "unknown"
  }
  
  return(acc_type)
}


parse_uniprot_fasta <- function (df, fasta) {
  
  na_genes_by_acc <- function(acc_lookup, acc_type) {
    acc_type <- tolower(acc_type)
    
    na_gene <- (is.na(acc_lookup$gene)) | (str_length(acc_lookup$gene) == 0)
    acc_lookup$gene <- as.character(acc_lookup$gene)
    acc_lookup$gene[na_gene] <- as.character(acc_lookup$prot_acc[na_gene])
    
    return(acc_lookup)
  }
  
  my_lookup <- c(
    "Homo sapiens" = "human",
    "Mus musculus" = "mouse",
    "Rattus norvegicus" = "rat",
    "Drosophila melanogaster" = "fly",
    "Bos taurus" = "bovine",
    "Canis lupus familiaris" = "dog",
    "Canis lupus" = "dog",
    "crap" = "crap"
  )
  
  
  stopifnot ("prot_acc" %in% names(df))

  acc_type <- parse_acc(df)
  
  if (acc_type == "unknown") {
    warning("Unknown accession; use hypothetical `refseq_acc`.")
    acc_type <- "refseq_acc"
  } else if (! acc_type %in% c("refseq_acc", "uniprot_id", "uniprot_acc")) {
    stop("The type of protein accesion needs to one of \'uniprot_id\', \'uniprot_acc\' or \'refseq_acc\'",
         call. = FALSE)
  }
  

  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- purrr::map(fasta, ~ {
        seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
      }) %>% do.call(`c`, .) 
      
      if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
        acc_lookup <- names(fasta) %>% 
          gsub("^..\\|", "", .) %>% 
          str_split("\\|", simplify = TRUE) %>% 
          data.frame() %>% 
          `colnames<-`(c("uniprot_acc", "uniprot_id")) %>% 
          bind_cols(data.frame(fasta_name = names(fasta)), .) %>% 
          dplyr::filter(.[[acc_type]] %in% unique(df$prot_acc)) %>% 
          dplyr::filter(!duplicated(.[[acc_type]]))
      } else if (acc_type == "refseq_acc") {
        acc_lookup <- tibble(fasta_name = names(fasta), refseq_acc = fasta_name) %>% 
          dplyr::filter(.[[acc_type]] %in% unique(df$prot_acc)) %>% 
          dplyr::filter(!duplicated(.[[acc_type]]))
      }
      
      fasta <- fasta %>% 
        .[names(.) %in% unique(acc_lookup$fasta_name)] %>% 
        .[! duplicated(names(.))]

      if (length(fasta) == 0) {
        stop("No fasta entries match protein accessions; probably wrong fasta file.", 
             call. = FALSE)
      } else {
        write.fasta(sequences = fasta, names = seqinr::getAnnot(fasta), nbchar = 80, 
                    file.out = file.path(dat_dir, "my_project.fasta"))
      }
      
      if (acc_type == "uniprot_acc") {
        names(fasta) <- gsub("^..\\|(.*)\\|.*$", "\\1", names(fasta))
      } else if (acc_type == "uniprot_id") {
        names(fasta) <- gsub("^.*\\|(.*)$", "\\1", names(fasta))
      } else if (acc_type == "refseq_acc") {
        names(fasta) <- gsub("\\.[^\\.]*$", "", names(fasta))
      }
      
      fasta_smry <- dplyr::bind_cols(
        prot_acc = names(fasta), 
        prot_desc = seqinr::getAnnot(fasta) %>% 
          purrr::map(., `[[`, 1) %>% 
          unlist(), 
        prot_mass = purrr::map_dbl(fasta, ~ {seqinr::getSequence(.x) %>% seqinr::pmw()}), 
        prot_len = getLength(fasta)
      ) %>% 
        dplyr::filter(!duplicated(.$prot_acc)) %>% 
        dplyr::mutate(acc_type = acc_type) %>% 
        dplyr::mutate(prot_mass = round(prot_mass, digits = 0))
      
      if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
        fasta_smry <- fasta_smry %>% 
          dplyr::mutate(prot_desc = gsub("^.*\\|.*\\s+?", "", prot_desc))
      } else if (acc_type == "refseq_acc") {
        fasta_smry <- fasta_smry %>% 
          dplyr::mutate(prot_desc = gsub("^.*_.*\\s+?", "", prot_desc))
      }
      
      acc_lookup <- fasta_smry %>% 
        dplyr::left_join(acc_lookup, by = c("prot_acc" = acc_type))
      rm(fasta_smry)
      
      if (acc_type %in% c("uniprot_acc", "uniprot_id")) {
        # add gene name
        genes <- acc_lookup %>% dplyr::select(prot_acc, prot_desc)
        
        na_genes <- genes %>% 
          dplyr::filter(!grepl("GN=", .$prot_desc)) %>% 
          dplyr::mutate(gene = NA)
        
        genes <- genes %>% 
          dplyr::filter(grepl("GN=", .$prot_desc)) %>% 
          dplyr::mutate(gene = gsub("^.*GN=(\\S+)\\s*.*", "\\1", prot_desc)) %>% 
          bind_rows(., na_genes)
        
        # add organism
        na_org <- genes %>% 
          dplyr::filter(!grepl("OS=", .$prot_desc)) %>% 
          dplyr::mutate(organism = NA)
        
        acc_lookup <- genes %>% 
          dplyr::filter(grepl("OS=", .$prot_desc)) %>% 
          dplyr::mutate(organism = gsub("^.*OS=(.*?)=.*$", "\\1", prot_desc)) %>% 
          dplyr::mutate(organism = gsub("\\s\\S*$", "", organism)) %>% 
          bind_rows(., na_org) %>% 
          dplyr::select(-prot_desc) %>% 
          dplyr::right_join(acc_lookup, by = "prot_acc")
        
        acc_lookup <- acc_lookup %>% 
          dplyr::mutate(species = my_lookup[.$organism]) %>% 
          add_entrez(.) %>% 
          na_genes_by_acc(acc_type)
      } else if (acc_type == "refseq_acc") {
        refseq_fns <- c("refseq_hs_crossref", "refseq_mm_crossref")
        data(package = "proteoQ", list = refseq_fns)
        
        refseq_crossref <- purrr::map(refseq_fns, ~ {
          get(.x)
        }) %>% do.call(`rbind`, .) %>% 
          dplyr::filter(!duplicated(.[[acc_type]]))
        
        stopifnot(acc_type %in% names(refseq_crossref))
        
        acc_lookup <- acc_lookup %>% 
          left_join(refseq_crossref, by = c("prot_acc" = acc_type)) %>% 
          dplyr::mutate(species = my_lookup[.$organism]) %>% 
          na_genes_by_acc(acc_type)
      }
      
    } else {
      stop("Wrong FASTA file path(s) or name(s).", call. = FALSE)
    }
  } else {
    stop("FASTA file not provided.")
  }
  
  save(acc_lookup, file = file.path(dat_dir, "acc_lookup.rda"))
  
  return(acc_lookup)
}


#' Adds protein annotation
#'
#' \code{annotPrn} cross-referencing proteins among \code{uniprot_acc},
#' \code{uniprot_id}, \code{refseq} and \code{entrez}.
#'
#' @import plyr dplyr purrr rlang seqinr stringr 
#' @importFrom magrittr %>% %$% %T>% 
annotPrn <- function (df, fasta) {
	acc_lookup <- parse_uniprot_fasta(df, fasta)
	
	acc_lookup <- bind_cols(
	  acc_lookup %>% 
	    dplyr::select(prot_acc), 
	  acc_lookup %>% 
	    dplyr::select(-which(names(.) %in% names(df))) %>% 
	    dplyr::select(-"organism", -"fasta_name"), 
	) 
	
	# (1) multiple uniprot_acc(s) can share the same entrez id
	# prot_acc    gene      acc_type   uniprot_acc species entrez
	# A1AT1_MOUSE Serpina1a uniprot_id P07758      mouse    20703
	# A1AT4_MOUSE Serpina1d uniprot_id Q00897      mouse    20703
		
	# (2) samee uniprot_acc can have different entrez ids
	# From	To
	# P02088	100503605
	# P02088	101488143
	# P02088	15129	
	
	# (3) ok that some uniprot_accs(s) have no corresponding entrez id

	df <- df %>% 
	  dplyr::mutate(psm_index = row_number()) %>% 
	  dplyr::left_join(acc_lookup, by = "prot_acc") %>% 
	  dplyr::filter(!duplicated(psm_index)) %>% 
	  dplyr::select(-psm_index)

	return(df)
}


#' Adds kinase annotation
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
annotKin <- function (df, acc_type) {
	stopifnot ("prot_acc" %in% names(df))
	
  data(package = "proteoQ", kinase_lookup)
  stopifnot(acc_type %in% names(kinase_lookup))

  lookup <- kinase_lookup %>% 
    dplyr::select(acc_type, kin_attr, kin_class, kin_order) %>%
    dplyr::filter(!duplicated(.)) %>% 
    dplyr::filter(!is.na(.[[acc_type]])) 
  
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
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
save_call <- function(call_pars, fn) {
	dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)
	call_pars[names(call_pars) == "..."] <- NULL
	save(call_pars, file = file.path(dat_dir, "Calls", paste0(fn, ".rda")))
}


#' Matches formulas to those in calls to pepSig or prnSig
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_fmls <- function(formulas) {
  fml_file <-  file.path(dat_dir, "Calls\\pepSig_formulas.rda")
  
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
         " not found in the latest call to 'pepSig(...)'.")
}


#' Matches the arg to anal_prnGSPA
#'
#' @param arg Argument to be matched.
#' @param call_rda the name of a rda
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
match_call_arg <- function (call_rda = "foo", arg = "scale_log2r") {
  call_rda <- rlang::as_string(rlang::enexpr(call_rda))
  arg <- rlang::as_string(rlang::enexpr(arg))
  
  rda <- paste0(call_rda, ".rda")
  file <- file.path(dat_dir, "Calls", rda)
  if (!file.exists(file)) stop(rda, " not found.")
  
  load(file = file)
  if (is.null(call_pars[[arg]])) 
    stop(arg, " not found in the latest call to ", call_rda, call. = FALSE)
  
  call_pars[[arg]]
}


#' Matches the name of GSPA result file
#'
#' @param anal_type Always \code{GSPA}; maybe different value for future uses.
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
match_gspa_filename <- function (anal_type = "GSPA", subdir = NULL, scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(!is.null(subdir))

  if (anal_type == "GSPA") {
    file <- file.path(dat_dir, "Calls\\anal_prnGSPA.rda")
    if (!file.exists(file)) stop("Run `prnGSPA` first.", call. = FALSE)
    load(file = file)
  }
  
  filename <- call_pars$filename
  if (rlang::is_empty(filename)) {
    filename <- list.files(path = file.path(dat_dir, "Protein\\GSPA", subdir), 
                           pattern = "^Protein_GSPA_.*\\.txt$", 
                           full.names = FALSE)
    
    if (rlang::is_empty(filename)) stop("No result file found under `", subdir, "`", call. = FALSE)
    
    if (scale_log2r) filename <- filename %>% .[grepl("^Protein_GSPA_Z", .)] else 
      filename <- filename %>% .[grepl("^Protein_GSPA_N", .)]
    
    if (impute_na) filename <- filename %>% .[grepl("Protein_GSPA_[NZ]_impNA\\.txt", .)] else 
      filename <- filename %>% .[grepl("Protein_GSPA_[NZ]\\.txt", .)]

    if (length(filename) > 1) stop("More than one result file found under `", subdir, "`", call. = FALSE)
    if (rlang::is_empty(filename)) 
      stop("No input files correspond to impute_na = ", impute_na, ", scale_log2r = ", scale_log2r, call. = FALSE)
  }

  return(filename)
}


#' Match scale_log2r
#'
#' \code{match_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in caller environment.
match_scale_log2r <- function(scale_log2r) {
  stopifnot(rlang::is_logical(scale_log2r))
  
  global_var <-tryCatch(global_var <-get("scale_log2r", envir = .GlobalEnv),
                        error = function(e) "e")
  if (global_var != "e" & is.logical(global_var)) scale_log2r <- global_var
  
  return(scale_log2r)
}


#' Match to a global logical variable
#'
#' @examples
#' scale_log2r <- TRUE
#' foo <- function(scale_log2r = FALSE) {
#'   match_logi_gv("scale_log2r", scale_log2r)
#' }
#' foo()
match_logi_gv <- function(var, val) {
  oval <- val
  gvar <-tryCatch(gvar <- get(var, envir = .GlobalEnv), error = function(e) "e")
  
  if (gvar != "e") {
    stopifnot(rlang::is_logical(gvar))
    if (gvar != oval) {
      warning("Coerce ", var, " to ", gvar, " after matching to the global setting.", 
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
match_prnSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  if (impute_na) {
    mscale_log2r <- match_call_arg(prnSig_impTRUE, scale_log2r)
  } else {
    mscale_log2r <- match_call_arg(prnSig_impFALSE, scale_log2r)
  }
  
  if (scale_log2r != mscale_log2r) {
    warning("scale_log2r = ", mscale_log2r, " after matching to prnSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }
  
  scale_log2r <- mscale_log2r
}


#' Match scale_log2r 
#'
#' \code{match_pepSig_scale_log2r} matches the value of \code{scale_log2r} to the value
#' in the most recent pepSig at a given impute_na status
match_pepSig_scale_log2r <- function(scale_log2r = TRUE, impute_na = FALSE) {
  stopifnot(rlang::is_logical(scale_log2r), rlang::is_logical(impute_na))
  
  if (impute_na) {
    mscale_log2r <- match_call_arg(pepSig_impTRUE, scale_log2r)
  } else {
    mscale_log2r <- match_call_arg(pepSig_impFALSE, scale_log2r)
  }
  
  if (scale_log2r != mscale_log2r) {
    warning("scale_log2r = ", mscale_log2r, " after matching to pepSig(impute_na = ", impute_na, ", ...)", 
            call. = FALSE)
  }
  
  scale_log2r <- mscale_log2r
}


#' Replaces NA genes
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
replace_na_genes <- function(df, acc_type) {
	acc_type <- tolower(acc_type)

	if (acc_type == "refseq_acc") {
		na_gene <- (is.na(df[, c("gene")])) | (str_length(df$gene) == 0)
		df$gene <- as.character(df$gene)
		df$gene[na_gene] <- as.character(df$prot_acc[na_gene])
	} else if (acc_type == "uniprot_id") {
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


#' Find peptide start and end positions
#'
#' \code{find_pep_pos} finds the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}.
#'
#' @import dplyr purrr rlang stringr seqinr tidyr
#' @importFrom magrittr %>% %$%
find_pep_pos <- function (prot_acc, pep_seq, fasta) {
  fasta_sub <- fasta %>% .[names(.) == prot_acc]
  pep_seq <- as.character(pep_seq)
  
  if (!rlang::is_empty(fasta_sub)) {
    pep_pos <- stringr::str_locate(fasta_sub, pattern = pep_seq)
    
    pos_bf <- pep_pos[1] - 1
    pos_af <- pep_pos[2] + 1
    
    pep_res_before <- stringr::str_sub(fasta_sub, pos_bf, pos_bf)
    pep_res_after <- stringr::str_sub(fasta_sub, pos_af, pos_af)
    
    # Mascot can alter the original sequence in fasta
    # prot_acc: "XP_003960355", original "QERFCQXK" becomes "QERFCQVK"
    if (any(is.na(c(pep_res_before, pep_res_after)))) {
      pep_pos <- cbind(pep_seq, pep_res_before = NA, start = NA, end = NA, 
                       pep_res_after = NA, prot_acc = prot_acc, is_tryptic = NA)
      
      return(pep_pos)
    }
    
    if (nchar(pep_res_before) == 0) pep_res_before <- "-"
    if (nchar(pep_res_after) == 0) pep_res_after <- "-"
    
    # ADVSLPSMQGDLK|NP_612429: not "E.ADVSLPSMQGDLK.T" but "K.ADVSLPSMQGDLK.T"
    if (pep_res_before %in% c("K", "R", "-")) { # the first match is tryptic
      pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, pep_res_after, prot_acc, is_tryptic = TRUE)
    } else if (pep_res_before == "M" & pep_pos[1] == 2) { # the first match is also tryptic
      pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, pep_res_after, prot_acc, is_tryptic = TRUE)
    } else { # the first match is non-tryptic
      pep_seq_new <- paste0(c("K", "R"), pep_seq)
      pep_pos_new_all <- purrr::map(pep_seq_new, ~ str_locate(fasta_sub, .x))
      ok_pos <- purrr::map_lgl(pep_pos_new_all, ~ !is.na(.x[[1]]))
      
      if (sum(ok_pos) > 0) { # tryptic match existed
        pep_pos_new <- pep_pos_new_all[[which(ok_pos)[1]]]
        
        pos_bf_new <- pep_pos_new[1]
        pos_af_new <- pep_pos_new[2] + 1
        
        pep_res_before_new <- stringr::str_sub(fasta_sub, pos_bf_new, pos_bf_new)
        pep_res_after_new <- stringr::str_sub(fasta_sub, pos_af_new, pos_af_new)
        
        pep_pos_new[1] <- pep_pos_new[1] + 1
        
        pep_pos <- cbind(pep_seq, pep_res_before_new, pep_pos_new, pep_res_after_new, prot_acc, is_tryptic = TRUE)        
      } else { # no tryptic matches
        pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, pep_res_after, prot_acc, is_tryptic = FALSE)
      }
    }
  } else { # no fasta matches
    pep_pos <- cbind(pep_seq, pep_res_before = NA, start = NA, end = NA, 
                     pep_res_after = NA, prot_acc = prot_acc, is_tryptic = FALSE)
  }
}


#' Annotation of peptide positions and adjacent amino acid residues
#'
#' \code{annotPeppos} annotates the start and the end positions of peptides in
#' ascribed proteins description based on the \code{fasta}. It also annotes the
#' preceding and the following AA residues.
#'
#' @import dplyr purrr rlang stringr seqinr tidyr
#' @importFrom magrittr %>% %$%
annotPeppos <- function (df, fasta){
  stopifnot(all(c("prot_acc", "pep_seq") %in% names(df)))
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)
  
  load(file = file.path(dat_dir, "label_scheme.rda"))

  # ok cases that same `pep_seq` but different `prot_acc`
  # (K)	MENGQSTAAK	(L) NP_510965
  # (-)	MENGQSTAAK	(L) NP_001129505  
  
  if (! "pep_seq_bare" %in% names(df)) {
    df <- df %>% 
      dplyr::mutate(pep_seq_bare = gsub("^.*\\.([^\\.]+)\\..*", "\\1", pep_seq))
  }

  df <- df %>% 
    dplyr::mutate(prot_acc = gsub("^.*\\|(.*)\\|.*$", "\\1", prot_acc)) %>% 
    dplyr::mutate(pep_prn = paste(pep_seq_bare, prot_acc, sep = "|"))

  df_pep_prn <- df %>% 
    dplyr::filter(!duplicated(pep_prn)) %>% 
    dplyr::select(c("pep_seq_bare", "prot_acc")) 
  
  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- purrr::map(fasta, ~ {
        seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
      }) %>% do.call(`c`, .) %>% 
        `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
        .[names(.) %in% unique(df_pep_prn$prot_acc)]
      
      if (length(fasta) == 0) {
        stop("No fasta entries matched protein accessions; probably wrong fasta file(s).", 
             call. = FALSE)
      }
      
      pep_pos_all <- purrr::map2(as.list(df_pep_prn$prot_acc), as.list(df_pep_prn$pep_seq_bare), 
                                 find_pep_pos, fasta) %>% 
        do.call(rbind, .) %>% 
        `colnames<-`(c("pep_seq_bare", "pep_res_before", "pep_start", "pep_end", "pep_res_after", 
                       "prot_acc", "is_tryptic")) %>% 
        data.frame(check.names = FALSE) %>% 
        tidyr::unite(pep_prn, pep_seq_bare, prot_acc, sep = "|", remove = TRUE)
    } else {
      stop("Wrong FASTA file path or name(s).", call. = FALSE)
    }
  } else {
    stop("FASTA file(s) not provided.")
  }
  
  rm(fasta)
  
  if ("pep_res_before" %in% names(df)) pep_pos_all$pep_res_before <- NULL
  if ("pep_res_after" %in% names(df)) pep_pos_all$pep_res_after <- NULL
  if ("pep_start" %in% names(df)) pep_pos_all$pep_start <- NULL
  if ("pep_end" %in% names(df)) pep_pos_all$pep_end <- NULL
  
  df$pep_seq_bare <- NULL

  df <- df %>% 
    dplyr::left_join(pep_pos_all, by = "pep_prn") %>% 
    dplyr::select(-pep_prn)
}


#' Subset fasta by accession type
#'
#' @import plyr dplyr purrr rlang seqinr
#' @importFrom magrittr %>%
subset_fasta <- function (df, fasta, acc_type) {
  stopifnot("prot_acc" %in% names(df))
  
  if (! acc_type %in% c("refseq_acc", "uniprot_id", "uniprot_acc")) {
    stop("The type of protein accesion needs to one of \'uniprot_id\', \'uniprot_acc\' or \'refseq_acc\'",
         call. = FALSE)
  }
  
  fasta <- purrr::map(fasta, ~ {
    seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
  }) %>% do.call(`c`, .)
  
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


#' Calculates protein percent coverage
#'
#' @import plyr dplyr purrr rlang seqinr
#' @importFrom magrittr %>%
calc_cover <- function(df, id, fasta = NULL) {
  stopifnot(all(c("prot_acc", "gene", "pep_start", "pep_end") %in% names(df)))

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
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "acc_lookup.rda"))
  
  if (length(fasta) == 0) {
    stop("No fasta entries matched the type of protein accession. Check the correctness of fasta file(s).", 
         call. = FALSE)
  }
  
  if (length(fasta) <= 200) {
    warning("Less than 200 entries in fasta matched by protein accession. 
            Make sure the fasta file is correct.")
  }
  
  df_sels <- df %>%
    dplyr::select(prot_acc, pep_start, pep_end) %>%
    dplyr::mutate(index = row_number()) %>% 
    dplyr::left_join(acc_lookup, by = "prot_acc") %>%
    dplyr::filter(!is.na(prot_len), !duplicated(index)) %>% 
    dplyr::select(-index)
  
  if (nrow(df_sels) == 0) stop("Probably incorrect accession types in the fasta file(s).", call. = FALSE)
  
  df_sels <- df_sels %>%
    dplyr::filter(pep_start <= prot_len) %>%
    dplyr::filter(pep_end <= prot_len) %>%
    split(.[["prot_acc"]], drop = TRUE) %>%
    purrr::map(function (.x) {
      len <- .x[1, "prot_len"]
      aa_map <- rep(NA, len)
      for (i in 1:nrow(.x)) aa_map[.x[i, ]$pep_start : .x[i, ]$pep_end] <- TRUE
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
      dplyr::group_by(gene) %>% 
      dplyr::summarise_all(~ max(.x, na.rm = TRUE))
    
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
    dplyr::select(-index) %>% 
    dplyr::mutate(prot_cover = round(prot_cover * 100, digits = 1)) %>%
    dplyr::mutate(prot_cover = paste0(prot_cover, "%"))

  return(df)
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
  suppressWarnings(
    data %>% 
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::left_join(metadata, by = "Sample_ID") %>%
      dplyr::mutate_at(vars(one_of("Color", "Fill", "Shape", "Size", "Alpha")), ~ as.factor(.)) %>%
      dplyr::select(which(not_all_NA(.))) %>% 
      rm_sglval_cols()    
  )
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
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ suppressWarnings(as.numeric(.x)))
}


#' Filter rows
#'
#' @param df; data frame. 
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
filters_in_call <- function (df, ...) {
  dots <- rlang::enexprs(...)
  nms <- names(dots)
  
  for (i in seq_along(dots)) {
    row_exprs <- dots[[nms[i]]] %>% 
      rlang::eval_bare()
    
    if (!rlang::is_list(row_exprs)) row_exprs <- list(row_exprs)
    
    row_vals <- row_exprs %>% 
      purrr::map(eval_tidy, df) %>% 
      purrr::reduce(`&`, .init = 1)
    
    row_vals <- row_vals %>% ifelse(is.na(.), FALSE, .)

    stopifnot(is.logical(row_vals))
    
    df <- df[row_vals, , drop = FALSE]
  }
  
  return(df)
}


#' Arrange rows
#'
#' @param df; data frame. 
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
arrangers_in_call <- function(.df, ..., .na.last = TRUE) {
  dots <- rlang::enexprs(...)
  nms <- names(dots)
  
  for (i in seq_along(dots)) {
    row_orders <- dots[[nms[i]]] %>% 
      rlang::eval_bare()
    
    if (!rlang::is_list(row_orders)) row_orders <- list(row_orders)
    
    order_call <- rlang::expr(order(!!!row_orders, na.last = !!.na.last))
    ord <- rlang::eval_tidy(order_call, .df)
    stopifnot(length(ord) == nrow(.df))

    .df <- .df[ord, , drop = FALSE]
  }
  
  return(.df)
}


#' Calculate PSM SDs
calc_sd_fcts_psm <- function (df, range_log2r, range_int, set_idx, injn_idx) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  label_scheme <- label_scheme %>% 
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
  
  run_scripts <- FALSE
  if (run_scripts) {
    df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::group_by(!!rlang::sym(id)) %>%
      dplyr::summarise_at(vars(starts_with(type)), ~ sd(.x, na.rm = TRUE))     
  }

  df %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id))) %>% 
    dplyr::arrange(!!rlang::sym(id)) %>% 
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarise_at(vars(starts_with(type)), ~ sd(.x, na.rm = TRUE)) # %>% 
    # dplyr::mutate(!!id := as.factor(!!rlang::sym(id)))
}



#' Violin plots of CV per TMT_Set and LCMS_injection
#' #'@import dplyr rlang ggplot2
sd_violin <- function(df = NULL, id = NULL, filepath = NULL, width = NULL, height = NULL, 
                      type = "log2_R", adjSD = FALSE, 
                      is_psm = FALSE, col_select = NULL, col_order = NULL, theme = NULL, ...) {

  err_msg1 <- paste0("\'Sample_ID\' is reserved. Choose a different column key.")
  
  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  
  col_select <- ifelse(is.null(col_select), rlang::expr(Select), rlang::sym(col_select))
  col_order <- ifelse(is.null(col_order), rlang::expr(Order), rlang::sym(col_order))
  
  if (col_select == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
  if (col_order == rlang::expr(Sample_ID)) stop(err_msg1, call. = FALSE)
  
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  if (is.null(label_scheme[[col_select]])) {
    stop("Column \'", rlang::as_string(col_select), "\' does not exist.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    stop("No samples were selected under column \'", rlang::as_string(col_select), "\'.",
         call. = FALSE)
  }
  
  if (is.null(label_scheme[[col_order]])) {
    warning("Column \'", rlang::as_string(col_order), "\' does not exist.
			Samples will be arranged by the alphebatic order.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_order]])) == 0) {
    warning("No samples were specified under column \'", rlang::as_string(col_order), "\'.",
            call. = FALSE)
  }
  
  label_scheme_sub <- label_scheme %>% 
    dplyr::mutate(new_id = paste0(TMT_Channel, " (", Sample_ID, ")")) %>% 
    dplyr::mutate(new_id = gsub("TMT-", "", new_id)) %>% 
    dplyr::select(Sample_ID, TMT_Set, new_id, !!col_select, !!col_order) %>%
    dplyr::filter(!is.na(!!col_select))

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  if (rlang::is_missing(width)) width <- 8
  if (rlang::is_missing(height)) height <- 8
  
  ymax <- eval(dots$ymax, env = caller_env())
  ybreaks <- eval(dots$ybreaks, env = caller_env())
  
  flip_coord <- eval(dots$flip_coord, env = caller_env())
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
      dblTrim(., range_log2r = c(0, 100), range_int = c(0, 100), type_r = "log2_R", type_int = "I")
    
    df_sd[, grep("^.*log2_R", names(df_sd))] <- df_sd[, grep("^.*log2_R", names(df_sd)), drop = FALSE] %>% 
      sweep(., 2, sqrt(SD), "/")
    
    df_z <- df_sd %>% dplyr::select(grep("^.*log2_R[0-9]{3}", names(.)))
    nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
    df_z[, nan_cols] <- 0
    df_sd[, grep("^.*_log2_R[0-9]{3}", names(df_sd))] <- df_z
    
    rm(df_z, nan_cols, SD)
  }

  df_sd <- df_sd %>% 
    `names<-`(gsub("^.*log2_R", "", names(.))) 
  
  if (!is_psm) {
    df_sd <- df_sd %>% dplyr::select(id, which(names(.) %in% label_scheme_sub$new_id))
  }

  Levels <- names(df_sd) %>% .[! . %in% id]

  if (!purrr::is_empty(Levels)) {
    df_sd <- df_sd %>%
      tidyr::gather(key = !!rlang::sym(id), value = "SD") %>%
      dplyr::rename(Channel := !!rlang::sym(id)) %>% 
      dplyr::ungroup(Channel) %>% 
      dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
      dplyr::filter(!is.na(SD))
    
    p <- ggplot() +
      geom_violin(df_sd, mapping = aes(x = Channel, y = SD, fill = Channel), size = .25, draw_quantiles = c(.95, .99)) +
      geom_boxplot(df_sd, mapping = aes(x = Channel, y = SD), width = 0.1, lwd = .2, fill = "white") +
      stat_summary(df_sd, mapping = aes(x = Channel, y = SD), fun.y = "mean", geom = "point",
                   shape=23, size=2, fill="white", alpha=.5) +
      labs(title = expression(""), x = expression("Channel"), y = expression("SD ("*log[2]*"FC)")) 
    
    if (!is.null(ymax)) {
      if (is.null(ybreaks)) {
        ybreaks <- ifelse(ymax > 1, 0.5, ifelse(ymax > 0.5, 0.2, 0.1))
      }
      p <- p + scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, ybreaks))
    }

    if (is.null(theme)) theme <- theme_psm_violin
    p <- p + theme

    if (flip_coord) {
      p <- p + coord_flip()
      width_temp <- width
      width <- height
      height <- width_temp
      rm(width_temp)
    }
    
    dots <- dots %>% .[! names(.) %in% c("width", "height", "in", "limitsize")]
    my_call <- rlang::expr(ggplot2::ggsave(filename = !!filepath, plot = !!p, 
                                           width = !!width, height = !!height, 
                                  units = "in", limitsize = FALSE, !!!dots))
    try(eval(my_call, caller_env()))
  }
}


#' Violin plots of reporter-ion intensity per TMT_Set and LCMS_injection
rptr_violin <- function(df, filepath, width, height) {
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
    dplyr::mutate(Intensity = round(Intensity, digit = 1))
  
  p <- ggplot() +
    geom_violin(df_int, mapping = aes(x = Channel, y = log10(Intensity), fill = Channel), size = .25) +
    geom_boxplot(df_int, mapping = aes(x = Channel, y = log10(Intensity)), width = 0.2, lwd = .2, fill = "white") +
    stat_summary(df_int, mapping = aes(x = Channel, y = log10(Intensity)), fun.y = "mean", geom = "point",
                 shape = 23, size = 2, fill = "white", alpha = .5) +
    labs(title = expression("Reporter ions"), x = expression("Channel"), y = expression("Intensity ("*log[10]*")")) + 
    geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
    theme_psm_violin
  
  try(ggplot2::ggsave(filepath, p, width = width, height = height, units = "in"))
}


#' geometric mean
my_geomean <- function (x, ...) {
  x <- log10(x) %>% mean(...)
  10^x
}


#' phospho counts
count_phosphopeps <- function() {
  df <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"), check.names = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(rowSums(!is.na( .[grep("^log2_R[0-9]{3}", names(.))] )) > 0)
  
  id <- match_call_arg(normPSM, group_psm_by)
  
  df_phos <- df %>% dplyr::filter(grepl("[sty]", .[[id]]))
  
  n_phos_peps <- nrow(df_phos)
  n_phos_sites <- stringr::str_count(df_phos[[id]], "[sty]") %>% sum()
  
  write.csv(
    data.frame(n_peps = n_phos_peps, n_sites = n_phos_sites), 
    file.path(dat_dir, "Peptide\\cache", "phos_pep_nums.csv"), 
    row.names = FALSE
  )
}


#' peptide miscleavage counts
count_pepmiss <- function() {
  dir.create(file.path(dat_dir, "PSM\\cache"), recursive = TRUE, showWarnings = FALSE)
  
  rmPSMHeaders()
  
  filelist = list.files(path = file.path(dat_dir, "PSM\\cache"),
                        pattern = "^F[0-9]{6}\\_hdr_rm.csv$")
  
  if (length(filelist) == 0) stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
  df <- purrr::map(filelist, ~ {
    data <- read.delim(file.path(dat_dir, "PSM\\cache", .x), sep = ',', check.names = FALSE, 
                       header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    
    data$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    data <- data %>% 
      dplyr::filter(!duplicated(pep_seq))
    
    tot <- nrow(data)
    mis <- data %>% dplyr::filter(pep_miss > 0) %>% nrow()
    
    tibble::tibble(total = tot, miscleavage = mis, percent = miscleavage/tot)
  }) %>% do.call(rbind, .)
  
  write.csv(df, file.path(dat_dir, "PSM\\cache\\miscleavage_nums.csv"), row.names = FALSE)
}


#' Row filtration helpers
#'
#' \code{contain_str}: contain a literal string; "PEPTIDES" contain_str "TIDE".
#' @param match A character string containing the pattern for matching.
#' @param vars A character string of the name of a variable. The default is
#'   FALSE.
#' @param ignore.case Logical; if TRUE, ignores case when matching.
#' @examples
#' \dontrun{
#' pepHist(
#'   col_select = BI,
#'   scale_log2r = TRUE,
#'   filter_peps = exprs(contain_chars_in("sty", pep_seq_mod)),
#'   scale_y = FALSE,
#'   ncol = 4,
#'   filename = "BI_pSTY_scaley_no.png",
#' )
#'
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
#' @rdname contain_str
#' @inheritParams purgePSM
#' @export
contain_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{not_contain_str}" not contain a literal string; "PEPTIDES"
#' not_contain_str "TED".
#' @rdname contain_str
#' @inheritParams purgePSM
#' @export
not_contain_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(match, vars, fixed = TRUE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{not_contain_chars_in}: not contain any of the characters in a literal
#' string; "PEPTIDES" not_contain_chars_in  "CAB".
#' @rdname contain_str
#' @inheritParams purgePSM
#' @export
not_contain_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(paste0("[", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}

#' Row filtration helpers
#'
#' \code{start_with_str}: start with a literal string. "PEPTIDES" start_with_str
#' "PEP".
#' @rdname contain_str
#' @inheritParams purgePSM
#' @export
start_with_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("^", match), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{end_with_str}: end with a literal string. "PEPTIDES" end_with_str
#' "TIDES".
#' @rdname contain_str
#' @inheritParams purgePSM
#' @export
end_with_str <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0(match, "$"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{start_with_chars_in}: start with one of the characters in a literal
#' string. "PEPTIDES" start_with_chars_in "XP".
#' @rdname contain_str
#' @inheritParams purgePSM
#' @export
start_with_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("^[", match, "]"), vars, fixed = FALSE, ignore.case)
}

#' Row filtration helpers
#'
#' \code{ends_with_chars_in}: end with one of the characters in a literal
#' string. "PEPTIDES" ends_with_chars_in "XS".
#' @rdname contain_str
#' @inheritParams purgePSM
#' @export
ends_with_chars_in <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[", match, "]$"), vars, fixed = FALSE, ignore.case)
}

rows_are_all <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  !grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}

rows_are_not_all <- function (match, vars, ignore.case = FALSE) {
  stopifnot(is_string(match), nchar(match) > 0)
  grepl(paste0("[^", match, "]"), vars, fixed = FALSE, ignore.case = FALSE)
}


#' Concatenate formula(e) to varargs of dots
concat_fml_dots <- function(fmls = NULL, fml_nms = NULL, dots = NULL, anal_type = "zzz") {
  if ((!is_empty(fmls)) & (anal_type == "GSEA")) return(c(dots, fmls))
  
  if (purrr::is_empty(fmls)) {
    fml_file <-  file.path(dat_dir, "Calls\\pepSig_formulas.rda")
    if (file.exists(fml_file)) {
      load(file = fml_file)
      
      if (!is.null(fml_nms)) {
        stopifnot(all(fml_nms %in% names(pepSig_formulas)))
        # stopifnot(all(names(pepSig_formulas) %in% fml_nms))
        # pepSig_formulas <- pepSig_formulas %>% .[map_dbl(names(.), ~ which(.x == fml_nms))]
        pepSig_formulas <- pepSig_formulas %>% .[names(.) %in% fml_nms]
      }
      
      dots <- c(dots, pepSig_formulas)
    } else {
      stop("Run both `pepSig()` and `prnSig()` first.")
    }
  } else {
    match_fmls(fmls)
    dots <- c(dots, fmls)
  }
  
  return(dots)
}


#' Roll up genes
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
      dplyr::group_by(gene) %>% 
      dplyr::mutate(prot_cover = as.numeric(sub("%", "", prot_cover))) %>% 
      dplyr::summarise_all(~ max(.x, na.rm = TRUE)) %>% 
      dplyr::mutate(prot_cover = paste0(prot_cover, "%"))
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
identical_dots <- function(call_nm, curr_dots, pattern) {
  file <- file.path(dat_dir, "Calls", paste0(call_nm, ".rda"))
  if (!file.exists(file)) return(FALSE)
  
  load(file = file)
  identical(call_pars %>% .[grepl(pattern, names(.))], curr_dots)
}


#' Complete cases among sample IDs in label_scheme_sub, not label_scheme
my_complete_cases <- function (df, scale_log2r, label_scheme_sub) {
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  rows <- df %>%
    dplyr::select(grep(NorZ_ratios, names(.))) %>%
    `colnames<-`(label_scheme$Sample_ID) %>%
    dplyr::select(which(names(.) %in% label_scheme_sub$Sample_ID)) %>% 
    complete.cases(.)
  
  if (sum(rows) == 0) stop("None of the cases are complete.", call. = FALSE)
  
  df <- df[rows, ]
}


#' my union of named list 
#' names will be kept after the unification
my_union <- function (x, y) {
  x %>% 
    .[! names(.) %in% names(y)] %>% 
    c(y)
}


#' find the base names of files
find_fn_bases <- function (filenames) {
  gsub("\\.[^.]*$", "", filenames) # %>% .[1] 
}


#' find the extensions of files
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
check_dots <- function (blacklist = NULL, ...) {
  dots <- rlang::enexprs(...)
  dups <- purrr::map_lgl(names(dots), ~ .x %in% blacklist)
  nms <- names(dots[dups])
  
  if (!rlang::is_empty(nms)) {
    stop("Do not use argument(s): ", reduce(nms, paste, sep = ", "), call. = FALSE)
  }
}


#' force 'complete_cases = TRUE' at 'impute_na = FALSE'
to_complete_cases <- function (complete_cases = FALSE, impute_na = FALSE) {
  warn_msg1 <- "Coerce `complete_cases = TRUE` at `impute_na = FALSE`."
  if (!(impute_na || complete_cases)) {
    complete_cases <- TRUE
    rlang::warn(warn_msg1)
  }
  return(complete_cases)
}
