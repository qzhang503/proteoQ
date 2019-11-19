#'Significance tests
#'
#'\code{proteoSigtest} performs significance tests aganist peptide or protein
#'\code{log2FC}. Users should avoid call the method directly, but instead use
#'the following wrappers.
#'
#'In general, special characters of \code{+} or \code{-} should be avoided from
#'contrast terms. Occasionally, such as in biological studies, it may be
#'convenient to use \code{A+B} to denote a condition of combined treatment of
#'\code{A} and \code{B} . In the case, one can put the term(s) containing
#'\code{+} or \code{-} into a pair of pointy brackets. The syntax in the
#'following hypothetical example will compare the effects of \code{A}, \code{B},
#'\code{A+B} and the average of \code{A} and \code{B} to control \code{C}:
#'
#'\code{prnSig(fml = ~ Term["A - C", "B - C", "<A + B> - C", "(A + B)/2 - C"])}
#'
#'Note that \code{<A + B>} stands for one sample and \code{(A + B)} has two
#'samples in it.
#'
#'@inheritParams  proteoHist
#'@inheritParams  proteoHM
#'@param filename A file name to output results. The default is
#'  \code{Peptide_pVals.txt} for peptides and \code{Protein_pVals} for proteins.
#'@param method Character string; the method of linear modeling. The default is
#'  \code{limma}. At \code{method = lm}, the \code{lm()} in base R will be used
#'  for models without random effects and the \code{\link[lmerTest]{lmer}} will
#'  be used for models with random effects.
#'@param var_cutoff Numeric; the cut-off in the variances of \code{log2FC}.
#'  Entries with variances smaller than the threshold will be removed from
#'  linear modeling.
#'@param pval_cutoff Numeric; the cut-off in significance \code{pVal}. Entries
#'  with \code{pVals} smaller than the threshold will be removed from multiple
#'  test corrections.
#'@param logFC_cutoff Numeric; the cut-off in \code{log2FC}. Entries with
#'  \code{log2FC} smaller than the threshold will be removed from multiple test
#'  corrections.
#'@param ... User-defined formulae for linear modeling. The syntax starts with a
#'  tilde, followed by the name of an available column key in
#'  \code{expt_smry.xlsx} and square brackets. The contrast groups are then
#'  quoted with multiple contrast groups separated by commas. Additive random
#'  effects are indicated by \code{+ (1|col_key_1) + (1|col_key_2)}... \cr \cr
#'  \code{filter_}: Logical expression(s) for the row filtration of data; also
#'  see \code{\link{normPSM}}.
#'@return The primary output is
#'  \code{~\\dat_dir\\Peptide\\Model\\Peptide_pVals.txt} for peptide data or
#'  \code{~\\dat_dir\\Protein\\Model\\Protein_pVals.txt} for protein data.
#'
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/pepseqmod_min.R
#'@example inst/extdata/examples/normPep_min.R
#'@example inst/extdata/examples/normPrn_min.R
#'@example inst/extdata/examples/imputeNA_examples.R
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoSigtest <- function (df = NULL, id = gene, scale_log2r = TRUE, filepath = NULL, filename = NULL,
											impute_na = TRUE, complete_cases = FALSE, method = "limma",
											var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), ...) {

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

	id <- rlang::enexpr(id)
	df <- rlang::enexpr(df)
	method <- rlang::enexpr(method)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()

	# force impute_na to TRUE for "lm" and "lmer"
	if (!impute_na & method != "limma") impute_na <- TRUE

	# Sample selection criteria:
	#   !is_reference under "Reference" ->
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]
	info_anal(df = !!df, id = !!id, scale_log2r = scale_log2r,
					filepath = !!filepath, filename = !!filename, impute_na = impute_na,
					anal_type = "Model")(complete_cases = complete_cases, method = !!method, 
					                     var_cutoff, pval_cutoff, logFC_cutoff, ...)
}


#' Row Variance
#'
#' @importFrom magrittr %>%
rowVars <- function (x, na.rm = TRUE) {
		sqr <- function(x) x * x
		n <- rowSums(!is.na(x))
		n[n <= 1] <- NA
		return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}


#' Filter rows by variance quantiles
#'
#' @import dplyr rlang
#' @importFrom magrittr %>%
filterData <- function (df, var_cutoff = 1E-3) {
	df <- df %>%
		dplyr::mutate(Variance = rowVars(.)) %>%
		dplyr::mutate(rowname = rownames(df))

	Quantile <- quantile(df$Variance, probs = var_cutoff, na.rm = TRUE)

	df %>%
		dplyr::filter(Variance >= pmax(Quantile, 1E-3)) %>%
		dplyr::select(-Variance) %>%
		tibble::column_to_rownames()
}


#' Prepare formulas
#'
#' @importFrom magrittr %>%
prepFml <- function(formula, label_scheme_sub, ...) {

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
  
  dots <- rlang::enexprs(...)

	fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
	len <- length(fml)

	key_col <- fml[len] %>% gsub("(.*)\\[\\s*\\~*.*\\].*", "\\1", .)

	label_scheme_sub <- label_scheme_sub %>% dplyr::filter(!is.na(!!sym(key_col)))

	if (grepl("\\[\\s*\\~\\s*", fml[len])) { # formula = ~ Term[ ~ V]
		base <- fml[len] %>% gsub(".*\\[\\s*\\~\\s*(.*)\\]", "\\1", .)
	} else {
		base <- NULL
	}

	if (!is.null(base)) { # formula = ~ Term[~V]
		new_levels <- label_scheme_sub[[key_col]] %>% levels()
		new_levels <- c(new_levels[new_levels == base], new_levels[new_levels != base])
		label_scheme_sub <- label_scheme_sub %>%
			dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col), levels = new_levels))
	} else if (!grepl("\\[", fml[len])) { # formula = log2Ratio ~ Term
		new_levels <- label_scheme_sub[[key_col]] %>% levels() # leveled by the alphebatic order
	} else { # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]
		new_levels <- NULL
	}

	if (!is.null(new_levels)) {
		contrs <- paste(new_levels[-1], new_levels[1], sep = "-")
		elements <- new_levels
	} else {
	  contrs <- fml[len] %>%
	    gsub(".*\\[(.*)\\].*", "\\1", .) %>%
	    gsub("\\\"", "", .) %>%
	    str_split(",\\s*", simplify = TRUE) %>%
	    as.character()

	  new_contrs <- fml[len] %>%
	    gsub("^.*\\[(.*)\\].*", "\\1", .) %>% # may have random terms at the end
	    gsub("\\\"", "", .) %>% 
	    str_split(",\\s*", simplify = TRUE) %>% 
	    gsub("\\s+", "", .) %>% 
	    gsub("<([^>]*?)\\+([^>]*?)>", "<\\1.plus.\\2>", .) %>% 
	    gsub("<([^>]*?)\\-([^>]*?)>", "<\\1.minus.\\2>", .) %>% 
	    gsub("[ <>]+", "", .)

	  new_elements <- new_contrs %>%
	    gsub("/[0-9]", "", .) %>% # (A+B+C)/3-D
	    gsub("[\\(\\)]", "", .) %>%
	    str_split("[\\+\\-]\\s*", simplify = TRUE) %>%
	    as.character() %>%
	    unique() %>% 
	    .[. != ""]
	  
	  elements <- new_elements %>% 
	    gsub(".plus.", "+", ., fixed = TRUE) %>% 
	    gsub(".minus.", "-", ., fixed = TRUE)
	  
	  cat("\ncontrs: ", contrs %>% as.character, "\n")
	  cat("new_contrs: ", new_contrs %>% as.character, "\n")
	  cat("elements: ", elements %>% as.character, "\n")
	  cat("new_elements: ", new_elements %>% as.character, "\n\n")		
	}

	label_scheme_sub_sub <- label_scheme_sub %>%
		dplyr::filter(!!sym(key_col) %in% elements) %>%
		dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col)))
	
	if (nrow(label_scheme_sub_sub) == 0) {
	  stop("No samples were found for formula ", formula, 
	       "\nCheck the terms under column ", key_col, call. = FALSE)
	}

	design <- model.matrix(~0+label_scheme_sub_sub[[key_col]]) %>%
		`colnames<-`(levels(label_scheme_sub_sub[[key_col]]))

	new_design_nms <- colnames(design) %>% 
	  gsub("+", ".plus.", ., fixed = TRUE) %>% 
	  gsub("-", ".minus.", ., fixed = TRUE)
	
	new_design <- design %>% 
	  `colnames<-`(new_design_nms)
	
	contr_mat <- makeContrasts(contrasts = new_contrs, levels = data.frame(new_design)) %>% 
	  `colnames<-`(contrs) %>% 
	  `rownames<-`(colnames(design))
	
	rm(new_design_nms, new_design)
	
	random_vars <- fml[len] %>%
		gsub("\\[.*\\]+?", "", .) %>%
		paste("~", .) %>%
		as.formula() %>%
		terms.formula(.) %>%
		attr(., "term.labels") %>%
		.[grepl("\\|", .)] %>%
		gsub("1\\s*\\|\\s*(.*)", "\\1", .)

	cat("random_vars: ", random_vars %>% as.character, "\n\n")
	
	return(list(design = design, contr_mat = contr_mat, key_col = key_col, random_vars = random_vars,
							label_scheme_sub_sub = label_scheme_sub_sub))
}


#' Adjusted pVal
#'
#' @importFrom magrittr %>%
my_padj <- function(df_pval, pval_cutoff) {
	df_pval %>%
		purrr::map(~ .x <= pval_cutoff)	%>%
		purrr::map(~ ifelse(!.x, NA, .x)) %>%
		purrr::map2(as.list(df_pval), `*`) %>%
		purrr::map(~ p.adjust(.x, "BH")) %>%
		dplyr::bind_cols() %>%
		`names<-`(gsub("pVal", "adjP", colnames(.))) %>%
		dplyr::mutate(rowname = rownames(df_pval)) %>%
		bind_cols(df_pval, .) %>%
		mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		tibble::column_to_rownames()
}


#' Model summary
#'
#' @importFrom magrittr %>% %$%
lm_summary <- function(pvals, log2rs, pval_cutoff, logFC_cutoff) {
	nms <- rownames(pvals)

	pass_pvals <- pvals %>% purrr::map(~ .x <= pval_cutoff)
	pass_fcs <- log2rs %>% purrr::map(~ abs(.x) >= log2(logFC_cutoff))
	pass_both <- purrr::map2(pass_pvals, pass_fcs, `&`) %>% 
	  purrr::map(~ ifelse(!.x, NA, .x))

	res_padj <- pvals %>%
		purrr::map2(pass_both, `*`) %>%
		purrr::map(~ p.adjust(.x, "BH")) %>%
		data.frame(check.names = FALSE) %>%
		`names<-`(gsub("pVal", "adjP", colnames(.))) %>%
		`rownames<-`(nms) %>%
		tibble::rownames_to_column() %>%
		dplyr::bind_cols(pvals, .) %>%
		mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>%
		tibble::column_to_rownames()

	log2rs <- log2rs %>%
		to_linfc() %>%
		`colnames<-`(gsub("log2Ratio", "FC", names(.))) %>%
		dplyr::bind_cols(log2rs, .) %>%
		dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2) %>%
		`rownames<-`(nms)

	cbind.data.frame(res_padj, log2rs)
}


#' Factorial model formula for interaction terms. All factors under one channel
#' in `label_scheme_sub`
#'
#' @importFrom MASS ginv
model_onechannel <- function (df, id, formula, label_scheme_sub, complete_cases, method, var_cutoff, 
                              pval_cutoff, logFC_cutoff, ...) {

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

	id <- rlang::as_string(rlang::enexpr(id))

	fml_ops <- prepFml(formula, label_scheme_sub, ...)
		contr_mat <- fml_ops$contr_mat
		design <- fml_ops$design
		key_col <- fml_ops$key_col
		random_vars <- fml_ops$random_vars
		label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub

	# keep the name list as some rows may drop in procedures such as filtration
	df_nms <- df %>%
		tibble::rownames_to_column(id) %>%
		dplyr::select(id)

	if (complete_cases) df <- df[complete.cases(df), ]

	# remove small-variance entries such as normalizers
	df <- df %>% filterData(var_cutoff)

	# ensure the same order in Sample_ID
	df <- df %>% 
	  dplyr::select(as.character(label_scheme_sub_sub$Sample_ID))

	if (length(random_vars) > 0) {
		# only use the first random variable
	  design_random <- label_scheme_sub_sub[[random_vars[1]]] 
		corfit <- duplicateCorrelation(df, design = design, block = design_random)
		fit <- df %>%
			lmFit(design = design, block = design_random, correlation = corfit$consensus) %>%
			contrasts.fit(contr_mat) %>%
			eBayes()
	} else {
		fit <- df %>%
			lmFit(design = design) %>%
			contrasts.fit(contr_mat) %>%
			eBayes()
	}

	print(design)
	print(contr_mat)
	
	# limma
	log2rs <- fit$coefficients %>%
		data.frame(check.names = FALSE) %>%
		`names<-`(paste0("log2Ratio (", names(.), ")"))

	pvals <- fit$p.value %>%
		data.frame(check.names = FALSE) %>%
		`names<-`(paste0("pVal (", names(.), ")"))

	res_lm <- lm_summary(pvals, log2rs, pval_cutoff, logFC_cutoff)

	if (method %in% c("lmer", "lme", "lm")) {
		fml_rhs <- gsub("\\[.*\\]", "", formula) %>% .[length(.)]
		new_formula <- as.formula(paste("log2Ratio", "~", fml_rhs))

		contr_mat_lm <- t(contr_mat) %>%
			MASS::ginv() %>%
			`colnames<-`(colnames(contr_mat)) %>%
			`rownames<-`(rownames(contr_mat))
		names(dimnames(contr_mat_lm)) <- c("Levels", "Contrasts")
		contr_mat_lm <- list(Cdn = contr_mat_lm) %>% `names<-`(key_col)

		smpl_levels <- names(df)
		contr_levels <- attributes(contr_mat_lm[[key_col]])$dimnames$Contrasts
		
		df_lm <- df %>%
			tibble::rownames_to_column(id) %>%
			tidyr::gather(-id, key = Sample_ID, value = log2Ratio) %>%
			dplyr::mutate(Sample_ID = factor(Sample_ID, levels = smpl_levels)) %>%
			dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", key_col, random_vars)], by = "Sample_ID") %>%
			dplyr::select(which(not_all_NA(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			tidyr::nest()

		if (!purrr::is_empty(random_vars)) {
			res_lm <- df_lm %>%
				dplyr::mutate(
				  model = purrr::map(
				    data, ~ lmerTest::lmer(data = .x, formula = new_formula, contrasts = contr_mat_lm))) %>%
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>%
			  tidyr::unnest(glance, keep_empty = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>%
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>%
				dplyr::mutate(term = gsub(key_col, "", term)) %>% 
			  dplyr::select(-data, -model) %>% 
				dplyr::mutate(term = factor(term, levels = contr_levels)) %>%
				tidyr::spread(term , p.value) %>%
				tibble::column_to_rownames(id) %>%
				`names<-`(paste0("pVal (", names(.), ")")) %>%
				lm_summary(log2rs, pval_cutoff, logFC_cutoff)
		} else {
			res_lm <- df_lm %>%
				dplyr::mutate(model = purrr::map(data, ~ lm(data = .x, formula = new_formula,
				                                            contrasts = contr_mat_lm))) %>%
				dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>%
			  tidyr::unnest(glance, keep_empty = TRUE) %>%	
				dplyr::filter(!grepl("Intercept", term)) %>%
				dplyr::select(-c("std.error", "estimate", "statistic")) %>% 
				dplyr::mutate(term = gsub(key_col, "", term)) %>% 
			  dplyr::select(-data, -model) %>% 
				dplyr::mutate(term = factor(term, levels = contr_levels))	%>% 
				tidyr::spread(term , p.value)	%>% 
				tibble::column_to_rownames(id) %>%
				`names<-`(paste0("pVal (", names(.), ")"))  %>% 
			  lm_summary(log2rs, pval_cutoff, logFC_cutoff)
		}
	}

	df_op <- res_lm %>% 
	  tibble::rownames_to_column(id) %>%
		dplyr::right_join(df_nms, by = id) %>%
		tibble::column_to_rownames(var = id)
}


#' Perform significance tests
#'
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$%
#' @importFrom broom.mixed tidy
sigTest <- function(df, id, label_scheme_sub, filepath, filename, complete_cases,
										method, var_cutoff, pval_cutoff, logFC_cutoff, ...) {

	id <- rlang::as_string(rlang::enexpr(id))
	method <- rlang::as_string(rlang::enexpr(method))

	dots <- rlang::enexprs(...)
	
	non_fml_dots <- dots[!map_lgl(dots, is_formula)]
	dots <- dots[map_lgl(dots, is_formula)]

	if (id %in% c("prot_acc", "gene")) {
		prnSig_formulas <- dots
		save(prnSig_formulas, file = file.path(dat_dir, "Calls", "prnSig_formulas.Rdata"))
		rm(prnSig_formulas)
	} else if (id %in% c("pep_seq", "pep_seq_mod")) {
		pepSig_formulas <- dots
		save(pepSig_formulas, file = file.path(dat_dir, "Calls", "pepSig_formulas.Rdata"))
		rm(pepSig_formulas)
	}

	df_op <- purrr::map(dots, 
	                    ~ model_onechannel(df, !!id, .x, label_scheme_sub, complete_cases, method, 
	                                       var_cutoff, pval_cutoff, logFC_cutoff, !!!non_fml_dots)) %>%
					do.call("cbind", .)
}


#'Significance tests of peptide \code{log2FC}
#'
#'\code{pepSig} is a wrapper of \code{\link{proteoSigtest}} for the
#'significance tests of peptide data
#'
#'@rdname proteoSigtest
#'
#'@import purrr
#'@export
pepSig <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\Model\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoSigtest)(id = pep_seq, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\Model\\log","pepSig_log.csv"), append = TRUE)
}


#'Significance tests of protein \code{log2FC}
#'
#'
#'\code{prnSig} is a wrapper of \code{\link{proteoSigtest}} for the significance
#'tests of protein data
#'
#'
#'@rdname proteoSigtest
#'
#'@import purrr
#' @examples
#' # ===================================
#' # Significance test
#' # ===================================
#' scale_log2r <- TRUE
#' 
#' # Mascot
#' # peptide significance tests
#' pepSig(
#'   impute_na = FALSE, 
#'   W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
#'   W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
#'   W16_vs_W2 = ~ Term_3["W16-W2"], 
#' )
#' 
#' prnSig(
#'   impute_na = FALSE,
#'   W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", "W2.JHU.TMT2-W2.JHU.TMT1", "W2.PNNL.TMT2-W2.PNNL.TMT1"], # batches
#'   W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # locations
#' )
#'
#' prnSig(
#'   impute_na = FALSE,
#'   W2_loc_2 = ~ Term["(W2.BI.TMT2+W2.BI.TMT1)/2 - (W2.JHU.TMT2+W2.JHU.TMT1)/2"], # locations with batch average
#' )
#'
#' prnSig(
#'   impute_na = TRUE,
#'   W2_vs_W16_fix = ~ Term_3["W16-W2"], # fixed effect
#'   W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set), # one fixed and one random effect
#' )
#' 
#' prnSig(
#'   impute_na = TRUE,
#'   method = lm,
#'   W2_vs_W16_fix = ~ Term_3["W16-W2"], # one fixed effect
#'   W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set), # one fixed and one random effect
#'   W2_vs_W16_mix_2 = ~ Term_3["W16-W2"] + (1|TMT_Set) + (1|Color), # one fixed and two random effects
#' )
#' 
#' 
#' # MaxQuant
#' # peptide significance tests
#' pepSig(
#'   impute_na = FALSE, 
#'   W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
#'   W16_vs_W2_course = ~ Term_2["W16-W2"], 
#' )
#' 
#' # protein significance tests
#' prnSig(
#'   impute_na = FALSE, 
#'   W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
#'   W16_vs_W2_course = ~ Term_2["W16-W2"], 
#' )
#' 
#'@import purrr
#'@export
prnSig <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Model\\log"), recursive = TRUE, showWarnings = FALSE)

  quietly_log <- purrr::quietly(proteoSigtest)(id = gene, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Model\\log","prnSig_log.csv"), append = TRUE)
}

