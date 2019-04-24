#' Significance tests
#'
#' \code{proteoSigtest} produces heat maps.  It is a wrapper round the call \code{info_anal(..., anal_type = "Heatmap")}.
#'
#' reads the data from either "\code{~\\Direcotry\\Peptide\\Peptide All.txt}" at \code{id = pep_seq_mod}, 
#'
#' @param id The name of a unique id (see \code{\link[proteoQ]{MDS}}).
#' @param scale_log2r Logical; if TRUE, rescales \code{log2-ratios} to the same scale of standard deviation for all samples.
#' @return Images stored under the file folders that are associated to \code{id}, \code{anal_type}.
#'
#' @examples
#' MA(
#' 	id = gene, 
#' 	scale_log2r = scale_log2r, 
#' )
#' 
#' # or use the form of functional factories 
#' 	my_prnMA <- info_anal(df = NULL, id = prot_acc, scale_log2r, filepath = NULL, filename = NULL, anal_type = "MA")
#' 	my_prnMA()
#' 	
#' 	my_pepMA <- info_anal(df = NULL, id = pep_seq_mod, scale_log2r, filepath = NULL, filename = NULL, anal_type = "MA")
#' 	my_pepMA()
#' 	
#' \dontrun{
#' }
#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
#' @export
proteoSigtest <- function (df = NULL, id = gene, scale_log2r = FALSE, filepath = NULL, filename = NULL, 
											impute_na = TRUE, complete_cases = FALSE, method = "limma", 
											var_cutoff = 1E-3, pval_cutoff = 1, logFC_cutoff = log2(1), ...) {
	

	id <- rlang::enexpr(id)
	
	# force impute_na to TRUE for "lm" and "lmer"
	if (!impute_na & method != "limma") impute_na <- TRUE
	
	# Sample selection criteria:
	#   !is_reference under "Reference"
	#   !is_empty & !is.na under the column specified by a formula e.g. ~Term["KO-WT"]	
	info_anal(df = df, id = !!id, scale_log2r = scale_log2r, 
					filepath = filepath, filename = filename, 
					impute_na = impute_na, 
					anal_type = "Model")(complete_cases, method, var_cutoff, pval_cutoff, logFC_cutoff, ...)
}


#' Row Variance
#' 
#' @importFrom magrittr %>%
rowVars <- function (x, na.rm = TRUE) { # row variance
		sqr <- function(x) x * x
		n <- rowSums(!is.na(x))
		n[n <= 1] = NA
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
prepFml <- function(formula, label_scheme_sub) {

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

		elements <- fml[len] %>% 
			gsub(".*\\[(.*)\\].*", "\\1", .) %>% 
			gsub("/[0-9]", "", .) %>% 
			gsub("\\\"", "", .) %>% 
			gsub("[\\(\\)]", "", .) %>% 
			str_split("[,\\+\\-]\\s*", simplify = TRUE) %>% 
			as.character() %>% 
			unique()
	}

	label_scheme_sub_sub <- label_scheme_sub %>% 
		dplyr::filter(!!sym(key_col) %in% elements) %>% 
		dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col)))

	design <- model.matrix(~0+label_scheme_sub_sub[[key_col]]) %>% 
		`colnames<-`(levels(label_scheme_sub_sub[[key_col]]))

	contr_mat <- makeContrasts(contrasts = contrs, levels = data.frame(design))

	random_vars <- fml[len] %>% 
		gsub("\\[.*\\]+?", "", .) %>% 
		paste("~", .) %>% 
		as.formula() %>% 
		terms.formula(.) %>% 
		attr(., "term.labels") %>% 
		.[grepl("\\|", .)] %>% 
		gsub("1\\s*\\|\\s*(.*)", "\\1", .)		

	return(list(design = design, contr_mat = contr_mat, key_col = key_col, random_vars = random_vars, 
							label_scheme_sub_sub = label_scheme_sub_sub))
}


#' Adjusted pVal
#' 
#' @importFrom magrittr %>% 
my_padj <- function(df_pval, pval_cutoff) {
	df_pval %>% 
		purrr::map(~ . <= pval_cutoff)	%>% 
		purrr::map(~ ifelse(!., NA, .)) %>% 
		purrr::map2(as.list(df_pval), `*`) %>% 
		purrr::map(~ p.adjust(., "BH")) %>% 
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
	
	pass_pvals <- pvals %>% purrr::map(~ . <= pval_cutoff)
	pass_fcs <- log2rs %>% purrr::map(~ abs(.) >= log2(logFC_cutoff))
	pass_both <- purrr::map2(pass_pvals, pass_fcs, `&`) %>% purrr::map(~ ifelse(!., NA, .))

	res_padj <- pvals %>% 
		purrr::map2(pass_both, `*`) %>% 
		purrr::map(~ p.adjust(., "BH")) %>% 
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

	
#' factorial model formula for interaction terms
#' all factors under one channel in label_scheme_sub
#' comparisons defined by contrasts
#' var_cutoff to remove low-variance entries
#' pval_cutoff for adjP
#' logFC_cutoff
#' 
#' #' @importFrom MASS ginv
model_onechannel <- function (df, id, formula, label_scheme_sub, complete_cases, method, var_cutoff, pval_cutoff, logFC_cutoff) {

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
	
	fml_ops <- prepFml(formula, label_scheme_sub)
		contr_mat <- fml_ops$contr_mat
		design <- fml_ops$design
		key_col <- fml_ops$key_col
		random_vars <- fml_ops$random_vars
		label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub

	# keep the name list as some rows may drop in procedures such as filtration
	df_nms <- df %>% 
		tibble::rownames_to_column(id) %>% 
		dplyr::select(id)
	
	if(complete_cases) df <- df[complete.cases(df), ]
	
	# remove small-variance entries such as normalizers
	df <- df %>% filterData(var_cutoff)

	# ensure the same order in Sample_ID
	df <- df %>% dplyr::select(as.character(label_scheme_sub_sub$Sample_ID)) 

	if(length(random_vars) > 0) {
		design_random <- label_scheme_sub_sub[[random_vars[1]]] # choose only the first random variable
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
	
	# limma
	log2rs <- fit$coefficients %>% 
		data.frame(check.names = FALSE) %>% 
		`names<-`(paste0("log2Ratio (", names(.), ")"))
		
	pvals <- fit$p.value %>% 
		data.frame(check.names = FALSE) %>% 
		`names<-`(paste0("pVal (", names(.), ")"))

	res_lm <- lm_summary(pvals, log2rs, pval_cutoff, logFC_cutoff)

	# lm or lmer
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

		if(!purrr::is_empty(random_vars)) {
			res_lm <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, ~ lmerTest::lmer(data = .x, formula = new_formula, contrasts = contr_mat_lm))) %>% 
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
				dplyr::mutate(term = gsub(key_col, "", term)) %>% 
				dplyr::mutate(term = factor(term, levels = contr_levels)) %>% 
				tidyr::spread(term , p.value) %>% 
				tibble::column_to_rownames(id) %>% 
				`names<-`(paste0("pVal (", names(.), ")")) %>% 
				lm_summary(log2rs, pval_cutoff, logFC_cutoff)
		} else {
			res_lm <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, ~ lm(data = .x, formula = new_formula, contrasts = contr_mat_lm))) %>% 
				dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term)) %>% 
				dplyr::select(-c("std.error", "estimate", "statistic")) %>% 
				dplyr::mutate(term = gsub(key_col, "", term)) %>% 
				dplyr::mutate(term = factor(term, levels = contr_levels)) %>% 
				tidyr::spread(term , p.value) %>% 
				tibble::column_to_rownames(id) %>% 
				`names<-`(paste0("pVal (", names(.), ")")) %>% 
				lm_summary(log2rs, pval_cutoff, logFC_cutoff)
		}	
	}
	
	df_op <- res_lm %>% tibble::rownames_to_column(id) %>% 
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
	
	dots = rlang::enexprs(...)

	if(id %in% c("prot_acc", "gene")) {
		prnSig_formulas <- dots
		save(prnSig_formulas, file = file.path(dat_dir, "Calls", "prnSig_formulas.Rdata"))
		rm(prnSig_formulas)
	} else if(id %in% c("pep_seq", "pep_seq_mod")) {
		pepSig_formulas <- dots
		save(pepSig_formulas, file = file.path(dat_dir, "Calls", "pepSig_formulas.Rdata"))
		rm(pepSig_formulas)
	}
	# load(file = file.path(dat_dir, "Calls\\prnSig_formulas.Rdata"))
	# prnSig_formulas <- prnSig_formulas %>% purrr::map(~ .[is_call(.)])

	df_op <- purrr::map(dots, ~ model_onechannel(df, id, ., label_scheme_sub, complete_cases,
							 method, var_cutoff, pval_cutoff, logFC_cutoff)) %>% 
					do.call("cbind", .) 
}


#'Significance Tests of Peptide \code{log2-ratios}
#'@seealso \code{\link{proteoSigtest}} for parameters
#'@export
pepSig <- function (...) {
	proteoSigtest(id = "pep_seq", ...)
}


#'Significance Tests of Protein \code{log2-ratios}
#'@seealso \code{\link{proteoSigtest}} for parameters
#'@export
prnSig <- function (...) {
	proteoSigtest(id = "gene", ...)
}
