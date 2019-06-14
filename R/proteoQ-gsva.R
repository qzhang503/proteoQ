#'MDS Plots
#'
#'\code{proteoMDS} visualizes the results from multidimensional scaling (MDS)
#'and principal component analysis (PCA).
#'
#'
#'@import dplyr rlang ggplot2 GSVA
#'@importFrom magrittr %>%
#'@export
prnGSVA <- function (id = gene,
										scale_log2r = TRUE, df = NULL, filepath = NULL, filename = NULL,
										impute_na = TRUE, complete_cases = FALSE, method = "limma",
										gset_nm = "go_sets", var_cutoff = .5, pval_cutoff = 1E-4,
										logFC_cutoff = log2(1.1), mx.diff = TRUE, ...) {

  scale_log2r <- match_logi_gv("scale_log2r")

  id <- rlang::enexpr(id)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	dots <- rlang::enexprs(...)

	if(purrr::is_empty(dots)) {
		load(file = file.path(dat_dir, "Calls\\prnSig_formulas.Rdata"))
		dots <- prnSig_formulas
	} else {
		match_fmls(dots)
	}

	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
	info_anal(df = !!df, id = !!id, scale_log2r = scale_log2r,
					filepath = !!filepath, filename = !!filename,
					impute_na = impute_na,
					anal_type = "GSVA")(complete_cases, method, gset_nm, var_cutoff, pval_cutoff,
					                    logFC_cutoff, mx.diff, !!!dots)
}


#' Perform GSVA tests
#'
#' logFC_cutoff subsets data for adjusted pvals
#'
#' @import limma stringr purrr tidyr dplyr rlang gage
#' @importFrom magrittr %>% %$%
#' @importFrom outliers grubbs.test
#' @importFrom broom.mixed tidy
gsvaTest <- function(df, id = "entrez", label_scheme_sub, filepath, filename, complete_cases = FALSE,
										method = "limma", gset_nm = "go_sets", var_cutoff = .5, pval_cutoff = 1E-4,
										logFC_cutoff = log2(1.1), mx.diff = TRUE, ...) {

		id <- rlang::as_string(rlang::enexpr(id))
		dots = rlang::enexprs(...)

		if(purrr::is_empty(dots))
			stop("Please supply formula(s) of contrasts for significance assessment of enchriment scores.",
			     call. = TRUE)

		stopifnot(all(gset_nm %in% c("go_sets", "kegg_sets", "c2_msig")))

		gsets <- dbs %>% .[names(.) %in% gset_nm]
		stopifnot(length(gsets) > 0)

		res_es <- df %>%
			filterData(var_cutoff = var_cutoff) %>%
			as.matrix()

		res_es <- purrr::map(
		  gsets, ~ GSVA::gsva(res_es, .x, min.sz = 10, mx.diff = mx.diff, verbose = FALSE, parallel.sz = 0))

		fn_prx <- gsub("\\.png$", "", filename) %>% paste0(., "_", gset_nm)

		purrr::walk2(res_es, fn_prx, ~ {
			.x %>%
				data.frame(check.names = FALSE) %>%
				tibble::rownames_to_column("term") %>%
				write.table(file.path(filepath, paste0(.y, "_ES.txt")), sep = "\t", col.names = TRUE,
				            row.names = FALSE)
		})

		for(i in seq_along(res_es)) {
			purrr::map(dots,
					~ model_onechannel(data.frame(res_es[[i]], check.names = FALSE), id = !!id,
					.x, label_scheme_sub, complete_cases, method = "limma",
					var_cutoff, pval_cutoff, logFC_cutoff)) %>%
				do.call("cbind", .) %>%
				tibble::rownames_to_column("term") %>%
				write.table(file.path(filepath, paste0(fn_prx[i], "_pVals.txt")), sep = "\t", col.names = TRUE,
				            row.names = FALSE)
		}

		files <- list.files(path = filepath, pattern = "Protein_GSVA.+_pVals.txt$", , full.names = TRUE)

		res <- do.call(rbind,
			lapply(
				files,
				read.csv, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#"
			)
		)

		kept_rows <- res %>%
			tibble::column_to_rownames("term") %>%
			dplyr::select(grep("pVal", names(.))) %>%
			dplyr::mutate(Kept = rowSums(!is.na (.)) > 0) %>%
			dplyr::select(Kept) %>%
			unlist

		pvals <- res[kept_rows, ] %>%
			`rownames<-`(.$term) %>%
			dplyr::select(grep("pVal", names(.))) %>%
			my_padj(pval_cutoff) %>%
			tibble::rownames_to_column("term")

		res <- res %>%
			tibble::column_to_rownames("term") %>%
			dplyr::select(-which(names(.) %in% names(pvals))) %>%
			tibble::rownames_to_column("term") %>%
			dplyr::right_join(pvals, by = "term")

		res %>%
			write.table(file.path(filepath, paste0(gsub("\\.png$", "", filename), "_pVals.txt")), sep = "\t",
			            col.names = TRUE, row.names = FALSE)

		purrr::walk(file.path(filepath, paste0(fn_prx, "_pVals.txt")), ~ file.remove(.x), full.names = TRUE)
}
