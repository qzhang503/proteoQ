#'MDS Plots
#'
#'\code{proteoMDS} visualizes the results from multidimensional scaling (MDS)
#'and principal component analysis (PCA).
#'
#'An Euclidean distance matrix of \code{log2-ratios} is returned by
#'\code{\link[base]{dist}}, followed by a metric (\code{\link[stats]{cmdscale}})
#'or non-metric (\code{\link[MASS]{isoMDS}}) MDS. The default is metric MDS with
#'the input dissimilarities being euclidean distances.
#'
#'\code{log2-ratios} are used in PCA (\code{\link[stats]{prcomp}}).
#'
#'@param id Character string to indicate the type of data (a column ID). Peptide
#'  data will be used at column \code{id = "pep_seq"} or \code{"pep_seq_mod"},
#'  and protein data at \code{id = "prot_acc"} or \code{"gene"}.
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param adjEucDist Logical; if TRUE, adjusts the inter-plex Euclidean distance
#'  by \eqn{1/sqrt(2)}. The option \code{adjEucDist = TRUE} may be suitable when
#'  \code{reference samples} in each TMT plex undergo approximately the same
#'  sample handling process as the rest of the samples. For instance,
#'  \code{reference samples} were split at the levels of protein lysates.
#'  Typically, \code{adjEucDist = FALSE} if \code{reference samples} were split
#'  near the end of a sample handling process, for instance, at the stages
#'  immediately before or after TMT labeling.
#'@param classical Logical; performs metric MDS at TRUE and non-metric MDS at
#'  FALSE.
#'@param show_ids Logical; if TRUE, shows the sample IDs in MDS and PCA plots.
#'@param annot_cols Column names available in \code{expt_smry.csv}. Values
#'  under the selected columns will be used to color-code sample IDs on the top
#'  of \code{EucDist} plots.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the name of the current function \code{call} and the value of
#'  \code{id}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... More parameters for plotting:
#'
#'  \code{mds_color}, character string to a column ID in \code{Label
#'  scheme.csv}; values under which will be used to color-code data points in
#'  \code{MDS/PCA} plots.
#'
#'  \code{mds_shape}, character string to a column ID in \code{Label
#'  scheme.csv}; values under which will be used to shape-code data points in
#'  \code{MDS/PCA} plots.
#'
#'  \code{mds_alpha}, character string to a column ID in \code{Label
#'  scheme.csv}; values under which will be used to control the transparency of
#'  data points in \code{MDS/PCA} plots.
#'@return The plots of MDS, PCA and Euclidean distance.
#'
#'@import dplyr rlang ggplot2 GSVA GSVAdata
#'@importFrom magrittr %>%
#'@export
prnGSVA <- function (id = gene, 
										scale_log2r = FALSE, df = NULL, filepath = NULL, filename = NULL, 
										impute_na = TRUE, complete_cases = FALSE, method = "limma", 
										gset_nm = "go_sets", var_cutoff = .5, pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), mx.diff = TRUE, ...) {
	
	id <- rlang::enexpr(id)
	
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
	info_anal(df = df, id = !!id, scale_log2r = scale_log2r, 
					filepath = filepath, filename = filename, 
					impute_na = impute_na, 
					anal_type = "GSVA")(complete_cases, method, gset_nm, var_cutoff, pval_cutoff, logFC_cutoff, mx.diff, !!!dots)
}


#' Perform GSVA tests
#' 
#' @import limma stringr purrr tidyr dplyr rlang gage
#' @importFrom magrittr %>% %$%
#' @importFrom outliers grubbs.test 
#' @importFrom broom.mixed tidy
gsvaTest <- function(df, id = "entrez", label_scheme_sub, filepath, filename, complete_cases = FALSE, 
										method = "limma", gset_nm = "go_sets", 
										var_cutoff = .5, pval_cutoff = 1E-4, logFC_cutoff = log2(1.1), # a subset for adjusted pVals
										mx.diff = TRUE, ...) {
		
		id <- rlang::as_string(rlang::enexpr(id))
		dots = rlang::enexprs(...)
		
		if(purrr::is_empty(dots)) 
			stop("Please supply formula(s) of contrasts for significance assessment of enchriment scores.", call. = TRUE)
		
		stopifnot(all(gset_nm %in% c("go_sets", "kegg_sets", "c2_msig")))
		
		gsets <- dbs %>% .[names(.) %in% gset_nm]
		stopifnot(length(gsets) > 0)

		res_es <- df %>% 
			filterData(var_cutoff = var_cutoff) %>% 
			as.matrix()
		
		res_es <- purrr::map(gsets, ~ GSVA::gsva(res_es, .x, min.sz = 10, mx.diff = mx.diff, verbose = FALSE, parallel.sz = 0)) # %>% do.call("rbind", .) 
		
		fn_prx <- gsub("\\.png$", "", filename) %>% paste0(., "_", gset_nm)
		
		purrr::walk2(res_es, fn_prx, ~ {
			.x %>% 
				data.frame(check.names = FALSE) %>% 
				tibble::rownames_to_column("term") %>% 
				write.table(file.path(filepath, paste0(.y, "_ES.txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
		})

		for(i in seq_along(res_es)) {
			purrr::map(dots, 
					~ model_onechannel(data.frame(res_es[[i]], check.names = FALSE), id = !!id, 
					.x, label_scheme_sub, complete_cases, method = "limma", 
					var_cutoff, pval_cutoff, logFC_cutoff)) %>% 
				do.call("cbind", .) %>% 
				tibble::rownames_to_column("term") %>% 
				write.table(file.path(filepath, paste0(fn_prx[i], "_pVals.txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
		}

		# --------
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
			write.table(file.path(filepath, paste0(gsub("\\.png$", "", filename), "_pVals.txt")), sep = "\t", col.names = TRUE, row.names = FALSE)

		purrr::walk(file.path(filepath, paste0(fn_prx, "_pVals.txt")), ~ file.remove(.x), full.names = TRUE)
}
