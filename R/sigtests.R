#' Row Variance
#' 
#' @param x A data frame.
#' @param na.rm The same as in \code{mean}.
#' @importFrom magrittr %>% %T>% %$% %<>% 
rowVars <- function (x, na.rm = TRUE) 
{
  sqr <- function(x) x * x
  n <- rowSums(!is.na(x))
  n[n <= 1] <- NA
  
  rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1)
}


#' Filters rows by variance quantiles
#' 
#' @inheritParams info_anal
#' @inheritParams gn_rollup
#' @inheritParams prnSig
#' @import dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
filterData <- function (df, cols = NULL, var_cutoff = 1E-3) 
{
  if (is.null(cols)) cols <- 1:ncol(df)
  
  if (length(cols) > 1L) {
    df <- df %>% 
      dplyr::select(cols) %>% 
      dplyr::mutate(Variance = rowVars(.), 
                    rowname = rownames(.)) 
    
    Quantile <- quantile(df$Variance, probs = var_cutoff, na.rm = TRUE)
    
    df <- df %>%
      dplyr::filter(Variance >= pmax(Quantile, 1E-3)) %>%
      tibble::remove_rownames() %>% 
      tibble::column_to_rownames() %>% 
      dplyr::select(-Variance)
  }
  
  invisible(df)
}


#' Prepare formulas
#' 
#' @inheritParams gspaTest
#' @inheritParams model_onechannel
#' @importFrom magrittr %>% %T>% %$% %<>% 
prepFml <- function(formula, label_scheme_sub, ...) 
{
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

	label_scheme_sub <- label_scheme_sub %>% 
	  dplyr::filter(!is.na(!!sym(key_col)))

	# formula = ~ Term[ ~ V]
	base <- if (grepl("\\[\\s*\\~\\s*", fml[len])) 
		fml[len] %>% gsub(".*\\[\\s*\\~\\s*(.*)\\]", "\\1", .)
	else 
		NULL

	if (!is.null(base)) { # formula = ~ Term[~V]
		new_levels <- label_scheme_sub[[key_col]] %>% levels()
		new_levels <- c(new_levels[new_levels == base], 
		                new_levels[new_levels != base])
		label_scheme_sub <- label_scheme_sub %>%
			dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col), levels = new_levels))
	} 
	else if (!grepl("\\[", fml[len])) { # formula = log2Ratio ~ Term
		new_levels <- label_scheme_sub[[key_col]] %>% levels() # leveled by the alphabetic order
	}
	else { # formula = ~ Term["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]
		new_levels <- NULL
	}

	if (!is.null(new_levels)) {
		contrs <- paste(new_levels[-1], new_levels[1], sep = "-")
		elements <- new_levels
	} 
	else {
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
	  
	  message("\ncontrs: ", contrs %>% as.character, "\n")
	  message("new_contrs: ", new_contrs %>% as.character, "\n")
	  message("elements: ", paste(elements, collapse = ", "), "\n")
	  message("new_elements: ", paste(new_elements, collapse = ", "), "\n\n")
	}

	label_scheme_sub_sub <- label_scheme_sub %>%
		dplyr::filter(!!sym(key_col) %in% elements) %>%
		dplyr::mutate(!!sym(key_col) := factor(!!sym(key_col)))
	
	if (nrow(label_scheme_sub_sub) == 0) {
	  stop("No samples were found for formula ", formula, 
	       "\nCheck the terms under column ", key_col, 
	       call. = FALSE)
	}

	design <- model.matrix(~0+label_scheme_sub_sub[[key_col]]) %>%
		`colnames<-`(levels(label_scheme_sub_sub[[key_col]]))

	new_design_nms <- colnames(design) %>% 
	  gsub("+", ".plus.", ., fixed = TRUE) %>% 
	  gsub("-", ".minus.", ., fixed = TRUE)
	
	new_design <- design %>% 
	  `colnames<-`(new_design_nms)
	
	contr_mat <- makeContrasts(contrasts = new_contrs, 
	                           levels = data.frame(new_design)) %>% 
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

	message("random_vars: ", random_vars %>% as.character, "\n\n")
	
	invisible(list(design = design, 
	               contr_mat = contr_mat, 
	               key_col = key_col, 
	               random_vars = random_vars,
	               label_scheme_sub_sub = label_scheme_sub_sub))
}


#' Adjusted pVal
#' 
#' @param df_pval A data frame containing pVals
#' @inheritParams prnSig 
#' @importFrom magrittr %>% %T>% %$% %<>% 
my_padj <- function(df_pval, pval_cutoff) 
{
	df_pval %>%
		purrr::map(~ .x <= pval_cutoff)	%>%
		purrr::map(~ ifelse(!.x, NA, .x)) %>%
		purrr::map2(as.list(df_pval), `*`) %>%
		purrr::map(~ p.adjust(.x, "BH")) %>%
		dplyr::bind_cols() %>%
		`names<-`(gsub("pVal", "adjP", colnames(.))) %>%
		dplyr::mutate(rowname = rownames(df_pval)) %>%
    dplyr::bind_cols(df_pval, .) %>%
    dplyr::mutate_at(.vars = grep("pVal\\s+", names(.)), format, 
                     scientific = TRUE, digits = 2) %>%
    dplyr::mutate_at(.vars = grep("adjP\\s+", names(.)), format, 
                     scientific = TRUE, digits = 2) %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames()
}


#' Model summary
#' 
#' @param pvals A data frame of pVal
#' @param log2rs A data frame of log2FC
#' @inheritParams prnSig
#' @importFrom magrittr %>% %T>% %$% %<>% 
lm_summary <- function(pvals, log2rs, pval_cutoff, logFC_cutoff, 
                       padj_method = "BH") 
{
	nms <- rownames(pvals)

	pass_pvals <- pvals %>% purrr::map(~ .x <= pval_cutoff)
	pass_fcs <- log2rs %>% purrr::map(~ abs(.x) >= log2(logFC_cutoff))
	pass_both <- purrr::map2(pass_pvals, pass_fcs, `&`) %>% 
	  purrr::map(~ ifelse(!.x, NA, .x))

	res_padj <- pvals %>%
		purrr::map2(pass_both, `*`) %>%
		purrr::map(~ p.adjust(.x, padj_method)) %>%
		data.frame(check.names = FALSE) %>%
		`names<-`(gsub("pVal", "adjP", colnames(.))) %>%
		`rownames<-`(nms) %>%
	  tibble::rownames_to_column() %>% 
		dplyr::bind_cols(pvals, .) %>%
	  dplyr::mutate_at(.vars = grep("pVal\\s+", names(.)), format, 
	                   scientific = TRUE, digits = 2) %>%
	  dplyr::mutate_at(.vars = grep("adjP\\s+", names(.)), format, 
	                   scientific = TRUE, digits = 2) %>% 
	  tibble::remove_rownames() %>% 
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
#' @param formula Language; the formula in linear modeling.
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnSig
#' @importFrom MASS ginv
model_onechannel <- function (df, id, formula, label_scheme_sub, complete_cases, 
                              method, padj_method, var_cutoff, 
                              pval_cutoff, logFC_cutoff, ...) 
{
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
  
  options(warn = 1L)
  
  dots <- rlang::enexprs(...)
  # lmFit_dots <- dots %>% .[. %in% c("method")]
  # eBayes_dots <- dots %>% .[. %in% c("proportion")]

	id <- rlang::as_string(rlang::enexpr(id))

	fml_ops <- prepFml(formula, label_scheme_sub, ...)
	contr_mat <- fml_ops$contr_mat
	design <- fml_ops$design
	key_col <- fml_ops$key_col
	random_vars <- fml_ops$random_vars
	label_scheme_sub_sub <- fml_ops$label_scheme_sub_sub
	
	local({
	  ss <- (colSums(design) == 1)
	  bads <- ss[ss]

	  if (length(bads)) 
	    warning("Single sample condition: ", paste0(names(bads), collapse = ", "), 
	            " under `", formula, "`.", 
	            call. = FALSE)
	})

	# keep the name list as rows may drop in filtration
	df_nms <- df %>%
		tibble::rownames_to_column(id) %>%
		dplyr::select(id)

	if (complete_cases) df <- df[complete.cases(df), ]
	
	df <- df %>% 
	  filterData(var_cutoff = var_cutoff) %>% 
	  dplyr::select(as.character(label_scheme_sub_sub$Sample_ID))
	
	local({
	  ncol <- ncol(df)
	  
	  if (ncol < 4L) 
	    warning(formula, ": the total number of samples is ", ncol, ".\n", 
	            "May need more samples for statistical tests.", 
	            call. = FALSE)
	})

	if (length(random_vars)) {
	  if (length(random_vars) > 1) {
	    warning("Uses only the first random variable: ", random_vars[1], 
	            call. = FALSE)
	  }
	  
	  design_random <- label_scheme_sub_sub[[random_vars[1]]] 
		corfit <- duplicateCorrelation(df, design = design, block = design_random)
		
		fit <- suppressWarnings(
		  df %>%
		    lmFit(design = design, 
		          block = design_random, 
		          correlation = corfit$consensus) %>%
		    contrasts.fit(contr_mat) %>%
		    eBayes()
		)
	} else {
	  fit <- suppressWarnings(
	    df %>%
	      lmFit(design = design) %>%
	      contrasts.fit(contr_mat) %>%
	      eBayes()
	  )
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

	res_lm <- lm_summary(pvals, log2rs, pval_cutoff, logFC_cutoff, padj_method)

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
			dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", key_col, random_vars)], 
			                 by = "Sample_ID") %>%
			dplyr::select(which(not_all_NA(.))) %>%
			dplyr::group_by(!!rlang::sym(id)) %>%
			tidyr::nest()

		if (length(random_vars)) {
		  if (!requireNamespace("broom.mixed", quietly = TRUE)) {
		    stop("\n====================================================================", 
		         "\nNeed package \"broom.mixed\" for this function to work.",
		         "\n====================================================================",
		         call. = FALSE)
		  }
		  
		  if (!requireNamespace("lmerTest", quietly = TRUE)) {
		    stop("\n============================================================", 
		         "\nNeed package \"lmerTest\" for this function to work.",
		         "\n============================================================",
		         call. = FALSE)
		  }
		  
		  res_lm <- df_lm %>%
				dplyr::mutate(
				  model = purrr::map(data, 
				    ~ lmerTest::lmer(data = .x, formula = new_formula, 
				                     contrasts = contr_mat_lm))) %>%
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>%
			  tidyr::unnest(glance, keep_empty = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>%
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>%
				dplyr::mutate(term = gsub(key_col, "", term)) %>% 
			  dplyr::select(-data, -model) %>% 
				dplyr::mutate(term = factor(term, levels = contr_levels)) %>%
				tidyr::spread(term , p.value) %>%
		    `rownames<-`(NULL) %>% 
				tibble::column_to_rownames(id) %>%
				`names<-`(paste0("pVal (", names(.), ")")) %>%
				lm_summary(log2rs, pval_cutoff, logFC_cutoff, padj_method)
		} else {
		  if (!requireNamespace("broom", quietly = TRUE)) {
		    stop("\n============================================================", 
		         "\nNeed package \"broom\" for this function to work.",
		         "\n============================================================",
		         call. = FALSE)
		  }
		  
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
		    `rownames<-`(NULL) %>% 
				tibble::column_to_rownames(id) %>%
				`names<-`(paste0("pVal (", names(.), ")"))  %>% 
			  lm_summary(log2rs, pval_cutoff, logFC_cutoff, padj_method)
		}
	}

  df_op <- res_lm %>% 
    tibble::rownames_to_column(id) %>%
    dplyr::right_join(df_nms, by = id) %>%
    `rownames<-`(NULL) %>% 
    tibble::column_to_rownames(var = id)	
}


#' Perform significance tests
#' 
#' @param data_type The type of data being either \code{Peptide} or \code{Protein}.
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnSig
#' @import limma stringr purrr tidyr dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
sigTest <- function(df, id, label_scheme_sub, 
                    scale_log2r, complete_cases, impute_na, rm_allna, 
                    filepath, filename, 
										method, padj_method, var_cutoff, pval_cutoff, logFC_cutoff, 
										data_type, anal_type, ...) 
{
  dat_dir <- get_gl_dat_dir()
  
  stopifnot(vapply(c(var_cutoff, pval_cutoff, logFC_cutoff), is.numeric, 
                   logical(1L)))
  
	id <- rlang::as_string(rlang::enexpr(id))
	method <- rlang::as_string(rlang::enexpr(method))

	dots <- rlang::enexprs(...)
	
	filter_dots <- dots %>% 
	  .[purrr::map_lgl(., is.language)] %>% 
	  .[grepl("^filter_", names(.))]
	
	arrange_dots <- dots %>% 
	  .[purrr::map_lgl(., is.language)] %>% 
	  .[grepl("^arrange_", names(.))]
	
	dots <- dots %>% 
	  .[! . %in% c(filter_dots, arrange_dots)]
	
	non_fml_dots <- dots[!purrr::map_lgl(dots, is_formula)]
	dots <- dots[purrr::map_lgl(dots, is_formula)]
	
	if (id %in% c("pep_seq", "pep_seq_mod")) {
	  pepSig_formulas <- dots
	  save(pepSig_formulas, file = file.path(dat_dir, "Calls/pepSig_formulas.rda"))
	  rm(pepSig_formulas)
	} 
	else if (id %in% c("prot_acc", "gene")) {
	  if (rlang::is_empty(dots)) {
	    prnSig_formulas <- dots <- concat_fml_dots()
	  } 
	  else {
	    prnSig_formulas <- dots
	  }
	  save(prnSig_formulas, file = file.path(dat_dir, "Calls", "prnSig_formulas.rda"))
	  rm(prnSig_formulas)
	}	

	fn_prefix2 <- ifelse(impute_na, "_impNA_pVals.txt", "_pVals.txt")
	
	df_op <- local({
	  dfw <- df %>% 
	    filters_in_call(!!!filter_dots) %>% 
	    arrangers_in_call(!!!arrange_dots) %>% 
	    prepDM(id = !!id, 
	           scale_log2r = scale_log2r, 
	           sub_grp = label_scheme_sub$Sample_ID, 
	           anal_type = anal_type, 
	           rm_allna = rm_allna) %>% 
	    .$log2R
	  
	  # `complete_cases` depends on lm contrasts
	  purrr::map(dots, ~ model_onechannel(dfw, !!id, .x, label_scheme_sub, 
	                                             complete_cases, method, 
	                                             padj_method, var_cutoff, 
	                                             pval_cutoff, logFC_cutoff, 
	                                             !!!non_fml_dots)) %>% 
	    do.call("cbind", .)
	})
	
	# record the `scale_log2r` status; otherwise, need to indicate it in a way
	# for example, `_N` or `_Z` in file names
	local({
  	dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)	  
	  
	  if (data_type == "Peptide") {
	    type <- "pep"
	  } else if (data_type == "Protein") {
	    type <- "prn"
	  } else {
	    stop("`data_type` needs to be either `Peptide` or `Protein`.", 
	         call. = FALSE)
	  }
	  
	  file <- paste0(type, "Sig_imp", ifelse(impute_na, "TRUE", "FALSE"), ".rda")
	  
	  call_pars <- c(scale_log2r = scale_log2r, 
	                 complete_cases = complete_cases, 
	                 impute_na = impute_na) %>% 
	    as.list()
	  
	  save(call_pars, file = file.path(dat_dir, "Calls", file))
	})
	
	suppressWarnings(
  	df_op <- df_op %>%
  	  tibble::rownames_to_column(id) %>% 
  	  dplyr::mutate(!!id := forcats::fct_explicit_na(!!rlang::sym(id))) %>% 
  	  dplyr::right_join(df, ., by = id) %T>% 
  	  write.table(file.path(filepath, paste0(data_type, fn_prefix2)), sep = "\t",
  	              col.names = TRUE, row.names = FALSE)	  
	)

	wb <- createWorkbook("proteoQ")
	addWorksheet(wb, sheetName = "Results")
	openxlsx::writeData(wb, sheet = "Results", df_op)
	saveWorkbook(wb, file = file.path(filepath, paste0(data_type, "_pVals.xlsx")), 
	             overwrite = TRUE) 

	invisible(df_op)
}


#'Significance tests of peptide/protein log2FC
#'
#'\code{pepSig} performs significance tests against peptide log2FC. 
#'
#'@rdname prnSig
#'
#'@import purrr
#'@export
pepSig <- function (scale_log2r = TRUE, impute_na = FALSE, complete_cases = FALSE, 
                    rm_allna = FALSE, method = c("limma", "lm"), padj_method = "BH", 
                    var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), 
                    df = NULL, filepath = NULL, filename = NULL, ...) 
{
  on.exit({
    mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
      c(rlang::enexprs(...)) %>% 
      save_call("pepSig")
  }, add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  method <- rlang::enexpr(method)
  if (method == rlang::expr(c("limma", "lm"))) {
    method <- "limma"
  } else {
    method <- rlang::as_string(method)
    stopifnot(method %in% c("limma", "lm"), length(method) == 1L)
  }
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  padj_method <- rlang::as_string(rlang::enexpr(padj_method))

  reload_expts()
  
  if ((!impute_na) && (method != "limma")) {
    impute_na <- TRUE
    warning("Coerce `impute_na = ", impute_na, "` at method = ", method, 
            call. = FALSE)
  }
  
  stopifnot(vapply(c(scale_log2r, impute_na, complete_cases, rm_allna), 
                   rlang::is_logical, logical(1)))
  
  stopifnot(vapply(c(var_cutoff, pval_cutoff, logFC_cutoff), 
                   is.numeric, logical(1L)))
  
  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            filepath = !!filepath, 
            filename = !!filename, 
            anal_type = "Model")(method = method, 
                                 padj_method = padj_method, 
                                 var_cutoff = var_cutoff, 
                                 pval_cutoff = pval_cutoff, 
                                 logFC_cutoff = logFC_cutoff, 
                                 rm_allna = rm_allna, 
                                 ...)
}


#'Significance tests of peptide/protein log2FC.
#'
#'\code{prnSig} performs significance tests against protein log2FC.
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
#'@inheritParams prnHist
#'@inheritParams prnHM
#'@inheritParams normPSM
#' @param filename A file name to output results. The default is
#'   \code{Peptide_pVals.txt} for peptides and \code{Protein_pVals} for
#'   proteins.
#' @param method Character string; the method of linear modeling. The default is
#'   \code{limma}. At \code{method = lm}, the \code{lm()} in base R will be used
#'   for models without random effects and the \code{\link[lmerTest]{lmer}} will
#'   be used for models with random effects.
#' @param padj_method Character string; the method of multiple-test corrections
#'   for uses with \link[stats]{p.adjust}. The default is "BH". See
#'   ?p.adjust.methods for additional choices.
#' @param var_cutoff Numeric; the cut-off in the variances of \code{log2FC}.
#'   Entries with variances smaller than the threshold will be excluded from
#'   linear modeling. The default is 1E-3.
#' @param pval_cutoff Numeric; the cut-off in significance \code{pVal}. Entries
#'   with \code{pVals} smaller than the threshold will be excluded from multiple
#'   test corrections. The default is at \code{1} to include all entries.
#' @param logFC_cutoff Numeric; the cut-off in \code{log2FC}. Entries with
#'   absolute \code{log2FC} smaller than the threshold will be removed from
#'   multiple test corrections. The default is at \code{log2(1)} to include all
#'   entries.
#' @param ... User-defined formulas for linear modeling. The syntax starts with
#'   a tilde, followed by the name of an available column key in
#'   \code{expt_smry.xlsx} and square brackets. The contrast groups are then
#'   quoted with one to multiple contrast groups separated by commas. The
#'   default column key is \code{Term} in \code{expt_smry.xlsx}: \cr \code{~
#'   Term["A - C", "B - C"]}. \cr \cr Additive random effects are indicated by
#'   \code{+ (1|col_key_1) + (1|col_key_2)}... Currently only a syntax of single
#'   contrast are supported for uses with random effects: \cr \code{~ Term["A -
#'   C"] + (1|col_key_1) + (1|col_key_2)} \cr \cr \code{filter_}: Logical
#'   expression(s) for the row filtration against data in a primary file of
#'   \code{Peptide[_impNA].txt} or \code{Protein[_impNA].txt}. See also
#'   \code{\link{normPSM}} for the format of \code{filter_} statements.
#' @return The primary output is \code{.../Peptide/Model/Peptide_pVals.txt} for
#'   peptide data or \code{.../Protein/Model/Protein_pVals.txt} for protein
#'   data. At \code{impute_na = TRUE}, the corresponding outputs are
#'   \code{Peptide_impNA_pvals.txt} or \code{Protein_impNA_pvals.txt}.
#'
#' @example inst/extdata/examples/prnSig_.R
#' @seealso 
#'  \emph{Metadata} \cr 
#'  \code{\link{load_expts}} for metadata preparation and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr 
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
#'  \code{\link{mergePep}} for extended examples in peptide data merging \cr 
#'  \code{\link{standPep}} for extended examples in peptide data normalization \cr
#'  \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
#'  \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
#'  
#'  \emph{Variable arguments of `filter_...`} \cr 
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \emph{Missing values} \cr 
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  
#'  \emph{Informatics} \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA visualization \cr 
#'  \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#' @import dplyr ggplot2
#' @importFrom magrittr %>% %T>% %$% %<>% 
#'
#' @export
prnSig <- function (scale_log2r = TRUE, impute_na = FALSE, complete_cases = FALSE, 
                    rm_allna = FALSE, method = c("limma", "lm"), padj_method = "BH", 
                    var_cutoff = 1E-3, pval_cutoff = 1.00, logFC_cutoff = log2(1), 
                    df = NULL, filepath = NULL, filename = NULL, ...) 
{
  on.exit({
    load(file.path(get_gl_dat_dir(), "Calls/prnSig_formulas.rda"))
    dots <- my_union(rlang::enexprs(...), prnSig_formulas)
    mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
      c(dots) %>% 
      save_call("prnSig")
  }, add = TRUE)
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  method <- rlang::enexpr(method)
  if (method == rlang::expr(c("limma", "lm"))) {
    method <- "limma"
  } else {
    method <- rlang::as_string(method)
    stopifnot(method %in% c("limma", "lm"), length(method) == 1L)
  }
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  if ((!impute_na) && (method != "limma")) {
    impute_na <- TRUE
    warning("Coerce `impute_na = ", impute_na, "` at method = ", method, 
            call. = FALSE)
  }
  
  stopifnot(vapply(c(scale_log2r, impute_na, complete_cases, rm_allna), 
                   rlang::is_logical, logical(1L)))
  
  stopifnot(vapply(c(var_cutoff, pval_cutoff, logFC_cutoff), 
                   is.numeric, logical(1L)))
  
  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            filepath = !!filepath, 
            filename = !!filename, 
            anal_type = "Model")(method = method, 
                                 padj_method = padj_method, 
                                 var_cutoff = var_cutoff, 
                                 pval_cutoff = pval_cutoff, 
                                 logFC_cutoff = logFC_cutoff, 
                                 rm_allna = rm_allna, 
                                 ...)
}


