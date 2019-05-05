#'GAGE analysis
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
prnGAGE <- function (id = gene,
										scale_log2r = TRUE, df = NULL, filepath = NULL, filename = NULL,
										impute_na = TRUE, complete_cases = FALSE,
										method = "limma", gset_nm = "go_sets", var_cutoff = 1E-3,
										pval_cutoff = .05, logFC_cutoff = log2(1.1), ...) {

  # scale_log2r <- match_logi_gv(scale_log2r)

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
					anal_type = "ESGAGE")(complete_cases, method, gset_nm, var_cutoff,
					                      pval_cutoff, logFC_cutoff, !!!dots)
}


#' Perform significance tests
#'
#' @import limma stringr purrr tidyr dplyr rlang gage
#' @importFrom magrittr extract %>% %$%
#' @importFrom outliers grubbs.test
#' @importFrom broom.mixed tidy
gageTest <- function(df, id, label_scheme_sub, filepath, filename, complete_cases,
										method = "limma", gset_nm = "go_sets", var_cutoff = 1E-3,
										pval_cutoff = .05, logFC_cutoff = log2(1.1), ...) {

  # load(file = file.path(dat_dir, "label_scheme.Rdata"))

  id <- rlang::as_string(rlang::enexpr(id))

	dots = rlang::enexprs(...)

	add_feature <- function(x, fea_nm, col_nm) {
		x <- x %>% data.frame(check.names = FALSE)
		x[[col_nm]] <- fea_nm

		return(x)
	}

	rm_outliers <- function(x) {
		# ind <- colnames(x) %in% as.character(label_scheme$Sample_ID)
		ind <- colnames(x) %in% as.character(label_scheme_sub$Sample_ID)

		x$p_geomean <- x[, ind] %>%
			log10(.) %>%
			data.frame(check.names = FALSE) %>%
			purrr::pmap(list) %>%
			purrr::map(., ~ Grubbs_outliers(.)) %>%
			dplyr::bind_rows() %>%
			rowMeans(., na.rm = TRUE) %>%
			10^.

		return(x)
	}

	cut_pvals <- function(x, pval_cutoff) {
		topN <- 50
		topN <- pmin(topN, nrow(x))

		x %>%
			dplyr::arrange(p.val) %>%
			dplyr::filter(row_number(p.val) <= topN, p.val < pval_cutoff) %>%
			dplyr::mutate(q_geomean = p.adjust(.$p_geomean, "BH")) %>%
			# dplyr::select(-which(names(.) %in% label_scheme$Sample_ID)) %>%
		  dplyr::select(-which(names(.) %in% label_scheme_sub$Sample_ID)) %>%
			dplyr::mutate(
				p.val = format(p.val, scientific = TRUE, digits = 2),
				q.val = format(q.val, scientific = TRUE, digits = 2),
				p_geomean = format(p_geomean, scientific = TRUE, digits = 2),
				q_geomean = format(q_geomean, scientific = TRUE, digits = 2)
			)
	}

	essetGrp <- function (res, samp, ref, compare, fn, x, gsets, type, direction, key_col) {
		dir.create(file.path(filepath, type, "Esstentials"), recursive = TRUE, showWarnings = FALSE)

		if (type == "GO") {
			if (direction == "greater") {
				op_type <-  "go.up"
				test4up <- TRUE
			} else if (direction == "less") {
				op_type <- "go.dn"
				test4up <- FALSE
			}
		} else if (type == "KEGG") {
			if (direction == "greater") {
				op_type <-  "kegg.up"
				test4up <- TRUE
			} else if (direction == "less") {
				op_type <- "kegg.dn"
				test4up <- FALSE
			}
		}

		res_sig <- res[[direction]] %>%
			data.frame(check.names = FALSE) %>%
			tibble::rownames_to_column() %>%
			# dplyr::filter(q.val < 0.05, rowSums(is.na(.)) == 0) %>%
			dplyr::filter(rowSums(is.na(.)) == 0) %>%
			tibble::column_to_rownames(var = "rowname")

		n_row <- nrow(res_sig)
		if (is.null(n_row)) n_row <- 0

		if (n_row > 1) {
			res_esg <- esset.grp(res_sig, x, gsets = gsets, ref = ref, samp = samp,
					test4up = test4up, output = TRUE, cutoff = 0.05, pc = 10e-10,
					outname = file.path(filepath, type, "Esstentials", paste(fn, op_type)), compare = compare)

			nms <- as.character(label_scheme$Sample_ID)

			outputs <- read.csv(file.path(filepath, type, "Esstentials", paste(fn, paste0(op_type, ".esgp.txt"))), check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")

			outputs$p_val <- outputs %>%
				dplyr::select(which(names(.) %in% nms)) %>%
				log10(.) %>%
				purrr::pmap(list) %>%
				purrr::map(., ~ Grubbs_outliers(.)) %>%
				dplyr::bind_rows() %>%
				rowMeans(., na.rm = TRUE) %>%
				10^.

			outputs <- outputs %>%
				dplyr::arrange(p.val) %>%
				dplyr::mutate(q_val = p.adjust(.$p_val, "BH")) %>%
				dplyr::select(-which(names(.) %in% nms)) %>%
				dplyr::mutate(contrast = fn) %>%
				dplyr::mutate(key_col = key_col) %>%
				dplyr::mutate(p.val = format(p.val, scientific = TRUE, digits = 2),
											q.val = format(q.val, scientific = TRUE, digits = 2),
											p_val = format(p_val, scientific = TRUE, digits = 2),
											q_val = format(q_val, scientific = TRUE, digits = 2))

			write.csv(outputs, file.path(filepath, type, "Esstentials", paste(fn, paste0(op_type, ".esgp.csv"))), row.names = FALSE)
			file.remove(file.path(filepath, type, "Esstentials", paste(fn, paste0(op_type, ".esgp.txt"))))

			# an exception in "esset.grp" function
			if(length(res_esg$essentialSets) == 1) {
				outputs$essentialSets[1] <- res_esg$essentialSets
				write.csv(outputs, file.path(filepath, type, "Esstentials", paste(fn, paste0(op_type, ".esgp.csv"))), row.names = FALSE) # update the output file
			}

		} else if (n_row == 1) { # single row not suitable for essentialSets analysis
			outputs <- cbind(essentialSets = rownames(res_sig), data.frame(res_sig, check.names = FALSE), setGroup = rownames(res_sig)) %>%
				dplyr::select(-dplyr::one_of("exp1")) %>% # compare == "as.group"
				dplyr::mutate(contrast = fn) %>%
				dplyr::mutate(key_col = key_col)

			write.csv(outputs, file.path(filepath, type, "Esstentials", paste(fn, paste0(op_type, ".esgp.csv"))), row.names = FALSE)
		} else { # zero row not suitable for essentialSets analysis
			res_esg <- NA
			outputs <- cbind(data.frame(essentialSets = NA, p.geomean = NA, stat.mean = NA, p.val = NA, q.val = NA, set.size = NA),
											data.frame(t(data.frame(colnames(x))) %>% `colnames<-`(.[1, ]), check.names = FALSE),
											data.frame(setGroup = NA), contrast = NA, key_col = NA) %>%
								dplyr::select(which(!names(.) %in% names(x)))
			outputs <- outputs[0, ]

			write.csv(outputs, file.path(filepath, type, "Esstentials", paste(fn, paste0(op_type, ".esgp.csv"))), row.names = FALSE)
		}

	}

	my_gage <- function (df, id, formula, label_scheme_sub, complete_cases) {
		# formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]  + (1|TMT_Set) + (1|Duplicate)
		# formula = ~ Term["Ner-V", "Ner_PLUS_PD-PD", "(Ner_PLUS_PD-PD)-(Ner-V)"]
		# formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|TMT_Set)
		# formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|Duplicate)
		# formula = log2Ratio ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]
		# formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|TMT_Set)
		# formula = ~ Term["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|Duplicate)
		# formula = ~ Term["Ner-PD", "V-PD", "(Ner_PLUS_PD-PD)"] + (1|Duplicate)
		# formula = log2Ratio ~ Term
		# formula = ~ Term[~V] # no interaction terms
		# formula = ~ Term

		fml_ops <- prepFml(formula, label_scheme_sub)

		contr_mat <- fml_ops$contr_mat
		key_col <- fml_ops$key_col
		random_vars <- fml_ops$random_vars
		label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub

		# housekeeping proteins with all-zero log2Ratios will be removed by the "Variance" criteria
		df_sub <- df %>%
			filterData(var_cutoff = var_cutoff) %>%
			dplyr::select(as.character(label_scheme_sub_sub$Sample_ID))

		if(complete_cases) df_sub <- df_sub[complete.cases(df_sub), ]

		fml_rhs <- gsub("\\[.*\\]", "", formula) %>% .[length(.)]
		new_formula <- as.formula(paste("log2Ratio", "~", fml_rhs))
		contr_mat_lm <- list(Cdn = contr_mat) %>% `names<-`(key_col)

		samps <- purrr::map(data.frame(contr_mat, check.names = FALSE), ~ {
			names(.x) <- rownames(contr_mat)
			label_scheme_sub_sub[[key_col]] %in% names(.x[.x == 1]) %>% which()
		})

		refs <- purrr::map(data.frame(contr_mat, check.names = FALSE), ~ {
			names(.x) <- rownames(contr_mat)
			label_scheme_sub_sub[[key_col]] %in% names(.x[.x == -1]) %>% which()
		})

		compares <- purrr::map2(samps, refs, ~ length(.x) == length(.y)) %>%
			ifelse("paired", "as.group") %>% as.list()

		# force "compares" to "unpaired"
		compares <- purrr::map(compares, function(x) x <- "unpaired")

		if(!purrr::is_empty(random_vars))
			cat("Random effects modeling not yet available in gene set analysis and will be omitted.")

		gsets <- dbs %>% .[names(.) %in% gset_nm]
		stopifnot(length(gsets) > 0)

		res <- purrr::map(gsets, function(gset) {
			purrr::pmap(list(samps, refs, compares),
					function(.x, .y, .z) {
						gage::gage(df_sub, gsets = gset, samp = .x, ref = .y, compare = .z) %>%
						magrittr::extract(c("greater", "less")) %>%
						purrr::map(~ data.frame(.x, check.names = FALSE)) %>%
						purrr::map(~ tibble::rownames_to_column(.x, "term")) %>%
						purrr::map2(., names(.), add_feature, "direction") %>%
						do.call(rbind, .)
					}
				) %>%
				purrr::map2(., names(.), add_feature, "contrast") %>%
				purrr::map(~ rm_outliers(.x)) %>%
				purrr::map(~ cut_pvals(.x, pval_cutoff)) %>%
				do.call(rbind, .)
			}
		) %>%
		do.call(rbind, .) %>%
		# dplyr::filter(complete.cases(.)) %>%
		write.table(file.path(filepath, paste0(gsub("\\.png$", "", filename), "_pVals.txt")),
			sep = "\t", col.names = TRUE, row.names = FALSE)

		# purrr::pwalk(list(res = res_go, samp = samps, ref = refs, compare = compares, fn = names(res_go)), essetGrp, x = df_sub, gsets = dbs$go_sets, type = "GO", direction = "greater", key_col = key_col)
		# purrr::pwalk(list(res = res_go, samp = samps, ref = refs, compare = compares, fn = names(res_go)), essetGrp, x = df_sub, gsets = dbs$go_sets, type = "GO", direction = "less", key_col = key_col)
		# purrr::pwalk(list(res = res_kegg, samp = samps, ref = refs, compare = compares, fn = names(res_kegg)), essetGrp, x = df_sub, gsets = dbs$kegg_sets, type = "KEGG", direction = "greater", key_col = key_col)
		# purrr::pwalk(list(res = res_kegg, samp = samps, ref = refs, compare = compares, fn = names(res_kegg)), essetGrp, x = df_sub, gsets = dbs$kegg_sets, type = "KEGG", direction = "less", key_col = key_col)

	}

	if(purrr::is_empty(dots)) stop("Need to supply formula(s) of contrasts for enchriment analysis.", call. = TRUE)

	df_op <- purrr::map(dots, ~ my_gage(df, id, ., label_scheme_sub, complete_cases)) %>%
		do.call("cbind", .)
}
