#' Perform significance tests
#' 
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$% 
#' @importFrom broom.mixed tidy
sigTest_v1 <- function(df, id, label_scheme_sub, filepath, filename, complete_cases, method, ...) {

	id <- rlang::ensym(id)
	Identifier <- rlang::as_string(id)
	
	dots = rlang::enexprs(...)
	
	restore_names <- function(x) {
		x %>% 
			gsub("_PLUS_", "+", .) %>% 
			gsub("_MINUS_", "-", .) %>% 
			gsub("_MULTIPLY_", "*", .) %>% 
			gsub("_DIVIDE_", "/", .)
	}
	
	change_names <- function(x) {
		x %>% 
			gsub("(.*)\\+(.*)", "\\1_PLUS_\\2", .) %>% 
			gsub("(.*)\\-(.*)", "\\1_MINUS_\\2", .) %>% 
			gsub("(.*)\\*(.*)", "\\1_MULTIPLY_\\2", .) %>% 
			gsub("(.*)\\/(.*)", "\\1_DIVIDE_\\2", .)
	}
		
	change_inquotes <- function (x) {
		x %>% 
			gsub("\\\"(.*)\\+(.*)\\\"+?", "\\\"\\1_PLUS_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\-(.*)\\\"+?", "\\\"\\1_MINUS_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\*(.*)\\\"+?", "\\\"\\1_MULTIPLY_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\/(.*)\\\"+?", "\\\"\\1_DIVIDE_\\2\\\"", .) %>% 
			gsub("\\\"", "", .)
	}
	
	parse_formula <- function (formula) {

		fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
		fml[1] <- fml[1] %>% change_names()
		fml[2] <- fml[2] %>% change_inquotes()
		
		new_formula <- as.formula(paste(fml[1], fml[2], sep = "~"))
		new_terms = attr(terms.formula(new_formula), "term.labels")
		new_terms_paren <- gsub("(1\\s+\\|\\s+.*)", "\\(\\1\\)", new_terms) # with parenthesis for random effects
		new_formula <- as.formula(paste(fml[1], sprintf(" ~ %s", paste(new_terms_paren, collapse = "+"))))

		# extract elements from the new_formula
		resp_var <- fml[1]		
		
		fixed_terms <- new_terms_paren %>% .[!grepl("\\|", .)] # variables and their interactions
		fixed_vars <- fixed_terms %>% .[!grepl("\\:", .)] # variables
		interac_vars <- fixed_terms %>% .[grepl("\\:", .)]
		new_formula_fixed <- as.formula(paste(resp_var, sprintf(" ~ %s", paste(fixed_terms, collapse = "+")))) # formula

		random_terms <- new_terms_paren %>% .[grepl("\\|", .)] # "1|variable"
		random_vars <- random_terms %>% gsub("\\(1\\s+\\|\\s+(.*)\\)$", "\\1", .) # variables
		random_terms_cmbn <- random_terms %>% paste(collapse = "+")
		
		if(purrr::is_empty(random_vars)) {
			lmer_formula <- as.formula(paste("log2Ratio", "~", "Model_Group"))
		} else {
			lmer_formula <- as.formula(paste("log2Ratio", "~", "Model_Group", "+", random_terms_cmbn))
		}

		label_scheme_sub_sub <- label_scheme_sub %>% 
			dplyr::mutate(Model_Group = change_names(Model_Group)) %>% 
			dplyr::filter(Model_Group %in% c(resp_var, fixed_vars))
			
		for (i in seq_along(label_scheme_sub_sub)) 
			assign(names(label_scheme_sub_sub[i]), as.vector(label_scheme_sub_sub[, i]))
		Model_Group <- Model_Group %>% as.character() %>% factor(levels = c(resp_var, fixed_vars))
		rm(i)
		
		design <- model.matrix(~0+Model_Group) %>% 
			data.frame(check.names = FALSE) %>% 
			`colnames<-`(gsub("Model_Group", "", colnames(.)))

		# ========================================
		# design matrix for lmFit
		# interaction terms will be removed for now
		if(purrr::is_empty(random_vars)) {
			limma_formula <- new_formula_fixed
			design_limma <- model.matrix(limma_formula, design)
		} else {
			limma_formula <- as.formula(paste(resp_var, sprintf(" ~ %s + %s", paste(fixed_terms, collapse = "+"), "Model_Dupcorr")))
			design_limma <- model.matrix(limma_formula, data.frame(design, Model_Dupcorr))
		}

		design_limma <- design_limma[, !duplicated(t(design_limma))] 
		ind <- design_limma %>% colSums() %>% `>`(0)
		design_limma <- design_limma[, ind]
		# ========================================

		# change names back
		label_scheme_sub_sub <- label_scheme_sub_sub %>% 
			dplyr::mutate(Model_Group = restore_names(Model_Group))
			
		colnames(design_limma) <- restore_names(colnames(design_limma))
		fixed_vars <- restore_names(fixed_vars)
		
		return(list(new_formula = new_formula, new_formula_fixed = new_formula_fixed, design_limma = design_limma, 
						lmer_formula = lmer_formula, resp_var = resp_var, fixed_vars = fixed_vars, 
						random_terms_cmbn = random_terms_cmbn, random_vars = random_vars, interac_vars = interac_vars, 
						label_scheme_sub_sub = label_scheme_sub_sub))
	}
	
	parse_formula_long <- function (formula) {

		fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
		
		if(fml[1] != "log2Ratio") {
			fml[1] <- fml[1] %>% change_names()
			fml[2] <- fml[2] %>% change_inquotes()
		}

		new_formula <- as.formula(paste(fml[1], fml[2], sep = "~"))
		new_terms = attr(terms.formula(new_formula), "term.labels")
		new_terms_paren <- gsub("(1\\s+\\|\\s+.*)", "\\(\\1\\)", new_terms) # with parenthesis for random effects
		new_formula <- as.formula(paste(fml[1], sprintf(" ~ %s", paste(new_terms_paren, collapse = "+"))))

		resp_var <- fml[1]		
		
		fixed_terms <- new_terms_paren %>% .[!grepl("\\|", .)] # variables and their interactions
		fixed_vars <- fixed_terms %>% .[!grepl("\\:", .)] # variables
		interac_vars <- fixed_terms %>% .[grepl("\\:", .)]
		new_formula_fixed <- as.formula(paste(resp_var, sprintf(" ~ %s", paste(fixed_terms, collapse = "+")))) # formula

		random_terms <- new_terms_paren %>% .[grepl("\\|", .)] # "1|variable"
		random_vars <- random_terms %>% gsub("\\(1\\s+\\|\\s+(.*)\\)$", "\\1", .) # variables
		random_terms_cmbn <- random_terms %>% paste(collapse = "+")
		
		if(purrr::is_empty(random_vars)) {
			new_limma_formula <- new_formula_fixed
		} else {
			new_limma_formula <-  as.formula(paste(resp_var, sprintf(" ~ %s", paste(paste(fixed_terms, collapse = "+"), paste(random_vars, collapse = "+"), sep = "+"))))
		}

		return(list(new_formula = new_formula, new_formula_fixed = new_formula_fixed, new_limma_formula = new_limma_formula, 
						resp_var = resp_var, fixed_vars = fixed_vars, 
						random_terms_cmbn = random_terms_cmbn, random_vars = random_vars, interac_vars = interac_vars))
	}
		
	rowVars <- function (x, na.rm = TRUE) { # row variance
			sqr <- function(x) x * x
			n <- rowSums(!is.na(x))
			n[n <= 1] = NA
			return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
	}
	
	my_lm <- function (f) {
		force(f)
		function (formula, data) try(f(formula, data))
	}
	
	longP <- function (df, formula, method, label_scheme_sub) {
		
		df_lm <- df %>% 
			dplyr::mutate(Variance = rowVars(.)) %>% 
			dplyr::mutate(!!id := rownames(df)) %>% 
			dplyr::filter(Variance >= 1E-3) %>% 
			tibble::column_to_rownames(Identifier) %>% 
			dplyr::select(which(names(.) %in% label_scheme_sub$Sample_ID)) %>% 
			tibble::rownames_to_column(Identifier) %>% 
			tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
			dplyr::left_join(label_scheme_sub, by = "Sample_ID") %>% 
			dplyr::select(which(not_all_NA(.))) %>% 
			dplyr::select(-Model_Sample) %>% 
			`names<-`(gsub("^Model_", "", names(.))) %>% 
			dplyr::group_by(!!id) %>% 
			tidyr::nest()

		if(method == "limma") {
			
			clean_colnames <- function (x, fixed_vars) {
				# x = colnames(design_limma)
				x_fixed <- x %>% .[!grepl(":", .)]

				for (i in seq_along(fixed_vars)) {
					x_fixed <- gsub(eval(expr(paste0("(", !!fixed_vars[[i]], ")", ".*"))), "\\1", x_fixed)
				}
				
			}

			formula_op <- parse_formula_long(formula)
			
			new_formula <- formula_op$new_formula
			new_terms_paren <- formula_op$random_terms_cmbn
			new_limma_formula <- formula_op$new_limma_formula
			random_vars <- formula_op$random_vars
			fixed_vars <- formula_op$fixed_vars
			
			design_limma <- model.matrix(new_limma_formula, df_lm$data[[1]])
			design_fixed_limma <- design_limma[, !colnames(design_limma) %in% random_vars, drop = FALSE]

			if((ncol(design_limma) - ncol(design_fixed_limma)) > 0) {
				design_random_limma <- design_limma[, colnames(design_limma) %in% random_vars, drop = FALSE] %>% 
					.[, 1, drop = FALSE] # choose the first one being the random effects

				corfit <- duplicateCorrelation(df, design = design_fixed_limma, block = design_random_limma)
				efit <- df %>% lmFit(design = design_fixed_limma, block = design_random_limma, correlation = corfit$consensus) %>% eBayes()
			} else {
				efit <- df %>% lmFit(design = design_fixed_limma) %>% eBayes()
			}

			df_pval <- efit$p.value %>% 
				data.frame(check.names = FALSE) %>% 
				dplyr::select(-grep("Intercept", colnames(.))) %>% 
				# dplyr::select( which(!names(.) %in% random_vars)) %>% 
				tibble::rownames_to_column(Identifier)

		} else if(method == "lm") {
			df_pval <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, my_lm(lm), formula = new_formula)) %>% 
				dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term)) %>% 
				dplyr::select(-c("estimate", "std.error", "statistic")) %>% 
				tidyr::spread(term , p.value)
		} else if(method %in% c("lme", "lmer")) {
			df_pval <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, my_lm(lmerTest::lmer), formula = new_formula)) %>% 
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
				tidyr::spread(term , p.value) 
		}
		
		# summarise results
		df_pval <- df_pval %>% 
			tibble::column_to_rownames(Identifier) %>% 
			`names<-`(paste0("pVal (", names(.), ")"))

		df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
			`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
			bind_cols()
			
		df_pval <- df_pval %>% 
			tibble::rownames_to_column(Identifier) %>% 
			bind_cols(df_padj)
			
		if (nrow(df_pval) == 1) df_pval$rowname <- rownames(df)[1]  # an exception in lmFit with only one entry for testing
		rm(df_padj)

		df_pval <- df_pval %>% 
			# mutate_at(.vars = grep("^log2Ratio|^FC \\(", names(.)), round, 2) %>% 
			mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			tibble::column_to_rownames(var = Identifier)
	}
	
	wideP <- function (df, formula, method, label_scheme_sub) {
		formula_op <- parse_formula(formula)
		
		label_scheme_sub_sub <- formula_op$label_scheme_sub_sub
		design_limma <- formula_op$design_limma
		resp_var <- formula_op$resp_var
		fixed_vars <- formula_op$fixed_vars
		interac_vars <- formula_op$interac_vars
		random_vars <- formula_op$random_vars
		
		new_formula <- formula_op$new_formula 
		new_formula_fixed <- formula_op$new_formula_fixed
		lmer_formula <- formula_op$lmer_formula

		df_sub <- df[, names(df) %in% label_scheme_sub_sub$Sample_ID]
		df_sub <- df_sub[, as.character(label_scheme_sub_sub$Sample_ID)]
		
		if(method == "limma") {
			
			# (1) without interaction terms
			design_fixed_limma <- design_limma[, !grepl("TMT_Set|Model_Dupcorr", colnames(design_limma))]

			if((ncol(design_limma) - ncol(design_fixed_limma)) > 0) {
				design_random_limma <- design_limma[, !colnames(design_limma) %in% c("(Intercept)", fixed_vars), drop = FALSE] %>% 
					.[, 1, drop = FALSE] # choose the first one being the random effects
					
				corfit <- duplicateCorrelation(df_sub, design = design_fixed_limma, block = design_random_limma)
				efit <- df_sub %>% lmFit(design = design_fixed_limma, block = design_random_limma, correlation = corfit$consensus) %>% eBayes()
			} else {
				efit <- df_sub %>% lmFit(design = design_fixed_limma) %>% eBayes()
			}
			
			rm(design_fixed_limma)
			
			# interaction terms
			if (!purrr::is_empty(interac_vars)) {
				my_model_matrix <- function (x) {
					out <- model.matrix(new_formula_fixed, x)
					attr(out, "assign") <- NULL
					
					out <- cbind(x[[resp_var]], out) 
					colnames(out)[1] <- resp_var
					out[, colnames(out) != "(Intercept)"]
				}
				
				replace_inf <- function (x) {
					replace(x, is.infinite(x), NA)
				}

				df_limma <- df_sub %>% 
					tibble::rownames_to_column(Identifier) %>% 
					tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
					dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", "Model_Group")], by = "Sample_ID") %>% 
					dplyr::mutate(Model_Group = factor(Model_Group, unique(Model_Group), levels = c(resp_var, fixed_vars))) %>% 
					dplyr::select(-Sample_ID) %>% 
					dplyr::group_by(!!id, Model_Group) %>% 
					dplyr::mutate(N = row_number()) %>% 
					dplyr::arrange(!!id, Model_Group, N) %>% # ensure the same order when spreading
					dplyr::ungroup(Model_Group) %>% 
					dplyr::mutate(Model_Group = as.character(Model_Group)) %>% 
					tidyr::complete(N, nesting(!!id, Model_Group)) %>% 
					tidyr::unite(key, !!id, N, remove = FALSE) %>% 
					tidyr::spread(Model_Group, log2Ratio) %>% 
					dplyr::select(-key, -N) %>% 
					`colnames<-`(change_names(colnames(.))) %>% 
					dplyr::mutate_all(funs(ifelse(is.na(.), -Inf, .))) %>% 
					dplyr::group_by(!!id) %>% 
					tidyr::nest() %>% 
					dplyr::mutate(data = purrr::map(data, my_model_matrix)) %>% 
					dplyr::mutate(data = purrr::map(data, replace_inf))

				nm <- rep(colnames(df_limma$data[[1]]), each = nrow(df_limma$data[[1]]))

				df_limma <- df_limma %>% 
					dplyr::mutate(flat_data = purrr::map(data, as.vector)) %>% 
					dplyr::mutate(flat_data = purrr::map(flat_data, set_names, nm)) 
					
				df_limma <- do.call("rbind", df_limma$flat_data) %>% 
					`rownames<-`(df_limma[[Identifier]]) %>% 
					`colnames<-`(restore_names(colnames(.)))
					
				nm_fct <- nm %>% 
					restore_names() %>% 
					factor(levels = c(resp_var, fixed_vars, interac_vars))

				design_interac <- model.matrix(~nm_fct) %>% 
					data.frame(check.names = FALSE) %>% 
					`colnames<-`(gsub("^nm_fct", "", colnames(.)))
					
				efit_interac <- df_limma %>% 
					lmFit(design = design_interac) %>% 
					eBayes()
					
				df_pval_interac <- efit_interac$p.value %>% 
					data.frame(check.names = FALSE) %>% 
					dplyr::select(-grep("Intercept", colnames(.))) %>% 
					dplyr::select(-c(fixed_vars)) %>% 
					tibble::rownames_to_column(Identifier)

			} else {
				efit_interac <- NULL
				df_pval_interac <- NULL
			}

			df_pval <- efit$p.value %>% 
				data.frame(check.names = FALSE) %>% 
				dplyr::select(-grep("Intercept", colnames(.))) %>% 
				dplyr::select(-grep("\\(TMT_Set\\)", names(.))) %>% 
				tibble::rownames_to_column(Identifier)
				
			if(!is.null(df_pval_interac)) df_pval <- df_pval %>% dplyr::left_join(df_pval_interac, by = Identifier)

		} else if(method == "lm") {
			my_lm <- function (x) {
				lm(new_formula_fixed, data = x)
			}
			
			df_lm <- df_sub %>% 
				tibble::rownames_to_column(Identifier) %>% 
				tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
				dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", "Model_Group")], by = "Sample_ID") %>% 
				# dplyr::mutate(Model_Group = factor(Model_Group, unique(Model_Group))) %>% 
				dplyr::mutate(Model_Group = factor(Model_Group, unique(Model_Group), levels = c(resp_var, fixed_vars))) %>% 
				dplyr::select(-Sample_ID) %>% 
				dplyr::group_by(!!id, Model_Group) %>% 
				dplyr::mutate(N = row_number()) %>% 
				dplyr::arrange(!!id, Model_Group, N) %>% # ensure the same order when spreading
				dplyr::ungroup(Model_Group) %>% 
				dplyr::mutate(Model_Group = as.character(Model_Group)) %>% 
				tidyr::complete(N, nesting(!!id, Model_Group)) %>% 
				tidyr::unite(key, !!id, N, remove = FALSE) %>% 
				tidyr::spread(Model_Group, log2Ratio) %>% 
				dplyr::select(-key, -N)

			df_lm <- df_lm %>% 
				`colnames<-`(change_names(colnames(.))) %>% 
				dplyr::group_by(!!id) %>% 
				tidyr::nest() %>% 
				dplyr::mutate(model = purrr::map(data, my_lm)) 

			df_pval <- df_lm %>% 
				dplyr::mutate(glance = map(model, broom::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term)) %>% 
				dplyr::select(-c("estimate", "std.error", "statistic")) %>% 
				tidyr::spread(term , p.value) %>% 
				`colnames<-`(restore_names(colnames(.)))
			
		} else if(method %in% c("lme", "lmer")) {
			my_lmer <- function (x) {
				lmerTest::lmer(lmer_formula, data = x)
			}

			df_lmer <- df_sub %>% 
				tibble::rownames_to_column(Identifier) %>% 
				tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
				dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", "Model_Group", random_vars)], by = "Sample_ID") %>% # assume random_vars exist...
				dplyr::mutate(Model_Group = factor(Model_Group, unique(Model_Group), levels = c(resp_var, fixed_vars))) %>% 
				dplyr::select(-Sample_ID) %>% 
				dplyr::group_by(!!id) %>% 
				tidyr::nest() %>% 
				dplyr::mutate(model = purrr::map(data, my_lmer))

			# library(broom.mixed)
			df_pval <- df_lmer %>% 
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
				tidyr::spread(term , p.value) %>% 
				`colnames<-`(gsub("^Model_Group", "", names(.)))
		}
				
		
		df_pval <- df_pval %>% 
			tibble::column_to_rownames(Identifier) %>% 
			`names<-`(paste0("pVal (", names(.), "/", resp_var, ")"))
			
			df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
				`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
				bind_cols()
				
			df_pval <- df_pval %>% 
				tibble::rownames_to_column(Identifier) %>% 
				bind_cols(df_padj)
				
			if (nrow(df_pval) == 1) df_pval$rowname <- rownames(df)[1]  # an exception in lmFit with only one entry for testing
			rm(df_padj)

			# may consider to calculate log2Ratio and FC using the original data before the imputation of NA
			resBase <- df_sub %>% 
				dplyr::select(as.character(label_scheme_sub_sub[label_scheme_sub_sub$Model_Group == resp_var, "Sample_ID"])) %>% 
				rowMeans(., na.rm = TRUE) # the mean log2Ratio of base samples

			df_op <- purrr::map(fixed_vars, function (x) {
				df_sub %>% 
					dplyr::select(as.character(label_scheme_sub_sub[label_scheme_sub_sub$Model_Group == x, "Sample_ID"])) %>% 
					rowMeans(., na.rm = TRUE) %>% 
					`-`(resBase)		
				}) %>% 
				data.frame() %>% 
				`names<-`( paste0("log2Ratio (", fixed_vars, "/", resp_var, ")")) 

			df_op <- cbind.data.frame(df_op, 
					df_op %>% 
						sapply(function(x) {ifelse(x > 0, 2^x, -1/(2^x))}) %>% 
						data.frame(check.names = FALSE) %>% 
						`colnames<-`(gsub("log2Ratio", "FC", names(.)))
				)
				
			df_op <- df_op %>% 
				tibble::rownames_to_column(Identifier) %>% 
				left_join(df_pval, by = Identifier)
			
			rm(df_pval)

			df_op <- df_op %>% 
				mutate_at(.vars = grep("^log2Ratio|^FC \\(", names(.)), round, 2) %>% 
				mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				tibble::column_to_rownames(var = Identifier)
	}


	
	
			not_run <- TRUE
			if (!not_run) {
				# https://stackoverflow.com/questions/22649536/model-matrix-with-all-pairwise-interactions-between-columns
				set.seed(1)
				dat <- data.frame(Y = rnorm(10), x = rnorm(10), y = rnorm(10), z = rnorm(10)+3)
				form <- Y ~ (x + y + z)^2
				design <- model.matrix(form, data = dat)
				efit <- lmFit(melt(dat)$value) %>% eBayes()
				
				# df = melt(design)
				# mod = lm(Var2 ~ value, df)

				#====================================================================================================
				# https://support.bioconductor.org/p/46153/
				set.seed(123)
				treat1 <- factor(rep(0:1, each = 5))
				treat2 <- sample(0:1, 10, TRUE)
				design <- model.matrix(~treat1*treat2)
				data <- c(sample(25:35, 5, TRUE), sample(35:55, 5, TRUE))
				efit <- lmFit(data, design) %>% eBayes()
				#====================================================================================================
				
				# makeContrasts
				modelMatrix(targets, ref="Ref")
				#  design <- modelMatrix(targets, ref="Ref")			
				# ------------------------------------------------------------
			}
			
			not_run <- TRUE
			if(!not_run) {
				fixed_expns <- fml[2] %>% gsub("\\+\\(.*", "\\1", .)

				if(grepl("\\\"", fixed_expns)) {
					# str_extract_all(fml[2], "\\\".*\\\"", simplify = TRUE)
					# ------------------------------------------------
					fixed_no_quo <- fixed_expns %>% str_split("\\+", simplify = TRUE) %>% .[!grepl("\\\"", .)]
					fixed_quo <- fixed_expns %>% str_split("\\\"", simplify = TRUE) %>% .[nchar(.) != 0] %>% .[!grepl("[\\+\\*]$", .)]
					fixed_expns <- c(fixed_no_quo, fixed_quo)
					# ------------------------------------------------
					
					random_expns <- fml[2] %>% 
						str_split("\\\"", simplify = TRUE) %>% 
						.[nchar(.) != 0] %>% 
						.[grepl("1\\|", .)] %>% 
						gsub("\\+?\\(1", "\\(1", .) %>% 
						str_split("\\(1\\|", simplify = TRUE) %>% 
						gsub("\\)$", "", .) %>% 
						.[nchar(.) != 0]
				} else {
					fixed_expns <- fixed_expns %>% str_split("\\+", simplify = TRUE)
					
					random_expns <- fml[2] %>% 
						str_split("\\\"", simplify = TRUE) %>% 
						str_split("[\\+]|[\\*]", simplify = TRUE) %>% 
						.[grepl("1\\|", .)] %>% 
						gsub("^\\+?", "", .) %>% 
						str_split("\\(1\\|", simplify = TRUE) %>% 
						gsub("\\)\\+?", "", .) %>% 
						.[nchar(.) != 0]
				}
				
			
			}
			

	
	df_op <- lapply(dots, 
		function(formula) 
		{
			formula_op <- parse_formula_long(formula)

			if(formula_op$resp_var == "log2Ratio") {
				df_op <- longP(df, formula, method, label_scheme_sub)
			} else {
				df_op <- wideP(df, formula, method, label_scheme_sub)
			}
		}
	) %>% do.call("cbind", .)
	
}



#' Perform significance tests
#' 
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$% 
#' @importFrom broom.mixed tidy
sigTest_master_col <- function(df, id, master_col, label_scheme_sub, filepath, filename, complete_cases, method, ...) {

	id <- rlang::ensym(id)
	Identifier <- rlang::as_string(id)
	
	master_col <- rlang::ensym(master_col) 
	master_col_nm <- master_col%>% rlang::as_string()

	dots = rlang::enexprs(...)
	
	restore_names <- function(x) {
		x %>% 
			gsub("_PLUS_", "+", .) %>% 
			gsub("_MINUS_", "-", .) %>% 
			gsub("_MULTIPLY_", "*", .) %>% 
			gsub("_DIVIDE_", "/", .)
	}
	
	change_names <- function(x) {
		x %>% 
			gsub("(.*)\\+(.*)", "\\1_PLUS_\\2", .) %>% 
			gsub("(.*)\\-(.*)", "\\1_MINUS_\\2", .) %>% 
			gsub("(.*)\\*(.*)", "\\1_MULTIPLY_\\2", .) %>% 
			gsub("(.*)\\/(.*)", "\\1_DIVIDE_\\2", .)
	}
		
	change_inquotes <- function (x) {
		x %>% 
			gsub("\\\"(.*)\\+(.*)\\\"+?", "\\\"\\1_PLUS_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\-(.*)\\\"+?", "\\\"\\1_MINUS_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\*(.*)\\\"+?", "\\\"\\1_MULTIPLY_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\/(.*)\\\"+?", "\\\"\\1_DIVIDE_\\2\\\"", .) %>% 
			gsub("\\\"", "", .)
	}
	
	parse_formula_wide <- function (formula, master_col_nm, label_scheme_sub) {

		fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
		fml[1] <- fml[1] %>% change_names()
		fml[2] <- fml[2] %>% change_inquotes()
		
		new_formula <- as.formula(paste(fml[1], fml[2], sep = "~"))
		new_terms = attr(terms.formula(new_formula), "term.labels")
		new_terms_paren <- gsub("(1\\s+\\|\\s+.*)", "\\(\\1\\)", new_terms) # with parenthesis for random effects
		new_formula <- as.formula(paste(fml[1], sprintf(" ~ %s", paste(new_terms_paren, collapse = "+"))))

		# extract elements from the new_formula
		resp_var <- fml[1]		
		
		fixed_terms <- new_terms_paren %>% .[!grepl("\\|", .)] # variables and their interactions
		fixed_vars <- fixed_terms %>% .[!grepl("\\:", .)] # variables
		interact_vars <- fixed_terms %>% .[grepl("\\:", .)]
		new_formula_fixed <- as.formula(paste(resp_var, sprintf(" ~ %s", paste(fixed_terms, collapse = "+")))) # formula

		random_terms <- new_terms_paren %>% .[grepl("\\|", .)] # "1|variable"
		random_vars <- random_terms %>% gsub("\\(1\\s+\\|\\s+(.*)\\)$", "\\1", .) # variables
		random_terms_cmbn <- random_terms %>% paste(collapse = "+")

		# -----------------------------------------------------------
		# design matrix for lmFit
		# interaction terms will be removed for now
		if(purrr::is_empty(random_vars)) {
			lmer_formula <- as.formula(paste("log2Ratio", "~", master_col_nm))
		} else {
			lmer_formula <- as.formula(paste("log2Ratio", "~", master_col_nm, "+", random_terms_cmbn))
		}
		
		label_scheme_sub_sub <- label_scheme_sub %>% 
			dplyr::mutate(!!master_col_nm := change_names(.[, master_col_nm])) %>% 
			# mutate_at(vars(master_col_nm), as.character) %>% 
			dplyr::filter(.[, master_col_nm] %in% c(resp_var, fixed_vars)) 
		
		design <- model.matrix(~0+label_scheme_sub_sub[[master_col_nm]]) %>% 
			data.frame(check.names = FALSE) %>% 
			`colnames<-`(gsub("label_scheme_sub_sub[[master_col_nm]]", "", colnames(.), fixed = TRUE))

		if(purrr::is_empty(random_vars)) {
			limma_formula <- new_formula_fixed
			design_limma <- model.matrix(limma_formula, design)
		} else {
			limma_formula <- as.formula(paste(resp_var, sprintf(" ~ %s + %s", paste(fixed_terms, collapse = "+"), paste(random_vars, collapse = "+"))))
			design_limma <- model.matrix(limma_formula, data.frame(design, label_scheme_sub_sub[, random_vars, drop = FALSE]))
		}
		
		design_limma <- design_limma[, !duplicated(t(design_limma))] 
		ind <- design_limma %>% colSums() %>% `>`(0)
		design_limma <- design_limma[, ind]
		# -----------------------------------------------------------

		# change names back
		label_scheme_sub_sub <- label_scheme_sub_sub %>% 
			dplyr::mutate(!!master_col_nm := restore_names(.[, master_col_nm]))
			
		colnames(design_limma) <- restore_names(colnames(design_limma))
		fixed_vars <- restore_names(fixed_vars)

		# for (i in seq_along(label_scheme_sub_sub)) assign(names(label_scheme_sub_sub[i]), as.vector(label_scheme_sub_sub[, i]))
		# Model_Group <- Model_Group %>% as.character() %>% factor(levels = c(resp_var, fixed_vars))
		# rm(i)

		return(list(new_formula = new_formula, new_formula_fixed = new_formula_fixed, design_limma = design_limma, 
						lmer_formula = lmer_formula, resp_var = resp_var, fixed_vars = fixed_vars, 
						random_terms_cmbn = random_terms_cmbn, random_vars = random_vars, interact_vars = interact_vars, 
						label_scheme_sub_sub = label_scheme_sub_sub))
	}
	
	parse_formula_long <- function (formula, label_scheme_sub) {

		fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
		
		if(fml[1] != "log2Ratio") {
			fml[1] <- fml[1] %>% change_names()
			fml[2] <- fml[2] %>% change_inquotes()
		}

		new_formula <- as.formula(paste(fml[1], fml[2], sep = "~"))
		new_terms = attr(terms.formula(new_formula), "term.labels")

		expr_vars <- new_terms %>% .[!grepl("\\:|\\|", .)] %>% gsub("(.*)\\[.*\\]", "\\1", .)

		# ----------------------------------------------------
		# label_scheme to be included for significance tests
		params <- new_terms %>% 
			.[!grepl("\\:|\\|", .)] %>% 
			gsub("\\]$|\\s+", "", .) %>% 
			stringr::str_split("[\\[|\\]|[,\\s+]]", simplify = TRUE) %>% 
			data.frame() %>%
			dplyr::mutate_if(is.factor, as.character) %>% 
			dplyr::mutate_all(funs(ifelse(nchar(.) == 0, NA, .))) %>% 
			unlist() %>% 
			.[!is.na(.)] %>% 
			.[! . %in% expr_vars]
		
		ind <- label_scheme_sub %>% 
			dplyr::select(which(names(.) %in% expr_vars)) %>% 
			dplyr::mutate_all(funs(ifelse(. %in% params, TRUE, NA))) %>% 
			purrr::pmap(list) %>% 
			purrr::map(anyNA) %>% 
			unlist() %>% 
			!.
		
		if(sum(ind) == 0) label_scheme_sub_sub <- label_scheme_sub else label_scheme_sub_sub <- label_scheme_sub[ind, ]
		# ---------------------------------------------

		fml[2] <- fml[2] %>% gsub("\\[.*\\]+?", "", .)
		
		new_formula <- as.formula(paste(fml[1], fml[2], sep = "~"))
		new_terms = attr(terms.formula(new_formula), "term.labels")
		new_terms_paren <- gsub("(1\\s+\\|\\s+.*)", "\\(\\1\\)", new_terms) # with parenthesis for random effects
		new_formula <- as.formula(paste(fml[1], sprintf(" ~ %s", paste(new_terms_paren, collapse = "+"))))

		resp_var <- fml[1]		
		
		fixed_terms <- new_terms_paren %>% .[!grepl("\\|", .)] # variables and their interactions
		fixed_vars <- fixed_terms %>% .[!grepl("\\:", .)] # variables
		interac_vars <- fixed_terms %>% .[grepl("\\:", .)]
		new_formula_fixed <- as.formula(paste(resp_var, sprintf(" ~ %s", paste(fixed_terms, collapse = "+")))) # formula

		random_terms <- new_terms_paren %>% .[grepl("\\|", .)] # "1|variable"
		random_vars <- random_terms %>% gsub("\\(1\\s+\\|\\s+(.*)\\)$", "\\1", .) # variables
		random_terms_cmbn <- random_terms %>% paste(collapse = "+")
		
		if(purrr::is_empty(random_vars)) {
			new_limma_formula <- new_formula_fixed
		} else {
			new_limma_formula <-  as.formula(paste(resp_var, sprintf(" ~ %s", paste(paste(fixed_terms, collapse = "+"), paste(random_vars, collapse = "+"), sep = "+"))))
		}

		return(list(new_formula = new_formula, new_formula_fixed = new_formula_fixed, new_limma_formula = new_limma_formula, 
						resp_var = resp_var, fixed_vars = fixed_vars, 
						random_terms_cmbn = random_terms_cmbn, random_vars = random_vars, interac_vars = interac_vars, 
						label_scheme_sub_sub = label_scheme_sub_sub))
	}
		
	rowVars <- function (x, na.rm = TRUE) { # row variance
			sqr <- function(x) x * x
			n <- rowSums(!is.na(x))
			n[n <= 1] = NA
			return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
	}
	
	my_lm <- function (f) {
		force(f)
		function (formula, data) try(f(formula, data))
	}
	
	longP <- function (df, formula, method, model_info) {
		df_lm <- df %>% 
			dplyr::mutate(Variance = rowVars(.)) %>% 
			dplyr::mutate(!!id := rownames(df)) %>% 
			dplyr::filter(Variance >= 1E-3) %>% 
			tibble::column_to_rownames(Identifier) %>% 
			dplyr::select(which(names(.) %in% model_info$Sample_ID)) %>% 
			tibble::rownames_to_column(Identifier) %>% 
			tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
			dplyr::left_join(model_info, by = "Sample_ID") %>% 
			dplyr::select(which(not_all_NA(.))) %>% 
			# dplyr::select(-Model_Sample) %>% 
			# `names<-`(gsub("^Model_", "", names(.))) %>% 
			dplyr::group_by(!!id) %>% 
			tidyr::nest()

		formula_op <- parse_formula_long(formula, model_info)
		new_formula <- formula_op$new_formula
		new_terms_paren <- formula_op$random_terms_cmbn
		new_limma_formula <- formula_op$new_limma_formula
		random_vars <- formula_op$random_vars
		fixed_vars <- formula_op$fixed_vars
		
		if(method == "limma") {
			design_limma <- model.matrix(new_limma_formula, df_lm$data[[1]])
			design_fixed_limma <- design_limma[, !colnames(design_limma) %in% random_vars, drop = FALSE]

			if((ncol(design_limma) > ncol(design_fixed_limma))) {
				design_random_limma <- design_limma[, colnames(design_limma) %in% random_vars, drop = FALSE] %>% 
					.[, 1, drop = FALSE] # choose the first one being the random effects

				corfit <- duplicateCorrelation(df, design = design_fixed_limma, block = design_random_limma)
				efit <- df %>% lmFit(design = design_fixed_limma, block = design_random_limma, correlation = corfit$consensus) %>% eBayes()
			} else {
				efit <- df %>% lmFit(design = design_fixed_limma) %>% eBayes()
			}

			df_pval <- efit$p.value %>% 
				data.frame(check.names = FALSE) %>% 
				dplyr::select(-grep("Intercept", colnames(.))) %>% 
				# dplyr::select( which(!names(.) %in% random_vars)) %>% 
				tibble::rownames_to_column(Identifier)

		} else if(method == "lm") {
			df_pval <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, my_lm(lm), formula = new_formula)) %>% 
				dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term)) %>% 
				dplyr::select(-c("estimate", "std.error", "statistic")) %>% 
				tidyr::spread(term , p.value)
		} else if(method %in% c("lme", "lmer")) {
			df_pval <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, my_lm(lmerTest::lmer), formula = new_formula)) %>% 
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
				tidyr::spread(term , p.value) 
		}
		
		# summarise results
		df_pval <- df_pval %>% 
			tibble::column_to_rownames(Identifier) %>% 
			`names<-`(paste0("pVal (", names(.), ")"))

		df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
			`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
			bind_cols()
			
		df_pval <- df_pval %>% 
			tibble::rownames_to_column(Identifier) %>% 
			bind_cols(df_padj)
			
		if (nrow(df_pval) == 1) df_pval$rowname <- rownames(df)[1]  # an exception in lmFit with only one entry for testing
		rm(df_padj)

		df_pval <- df_pval %>% 
			# mutate_at(.vars = grep("^log2Ratio|^FC \\(", names(.)), round, 2) %>% 
			mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			tibble::column_to_rownames(var = Identifier)
	}
	
	wideP <- function (df, formula, method, master_col, model_info) {
		master_col <- rlang::ensym(master_col) 
		master_col_nm <- master_col %>% rlang::as_string()

		formula_op <- parse_formula_wide(formula, master_col_nm, model_info)

		# model_info <- formula_op$model_info
		design_limma <- formula_op$design_limma
		resp_var <- formula_op$resp_var
		fixed_vars <- formula_op$fixed_vars
		interact_vars <- formula_op$interact_vars
		random_vars <- formula_op$random_vars
		
		new_formula <- formula_op$new_formula 
		new_formula_fixed <- formula_op$new_formula_fixed
		lmer_formula <- formula_op$lmer_formula

		if(method == "limma") {
			
			# (1) without interaction terms
			design_fixed_limma <- design_limma[, colnames(design_limma) %in% c("(Intercept)", fixed_vars)]

			if((ncol(design_limma) - ncol(design_fixed_limma)) > 0) {
				design_random_limma <- design_limma[, !colnames(design_limma) %in% c("(Intercept)", fixed_vars), drop = FALSE] %>% 
					.[, 1, drop = FALSE] # choose the first one being the random effects
					
				corfit <- duplicateCorrelation(df, design = design_fixed_limma, block = design_random_limma)
				efit <- df %>% lmFit(design = design_fixed_limma, block = design_random_limma, correlation = corfit$consensus) %>% eBayes()
			} else {
				efit <- df %>% lmFit(design = design_fixed_limma) %>% eBayes()
			}

			rm(design_fixed_limma)
			
			# interaction terms
			if (!purrr::is_empty(interact_vars)) {
				my_model_matrix <- function (x) {
					out <- model.matrix(new_formula_fixed, x)
					attr(out, "assign") <- NULL
					
					out <- cbind(x[[resp_var]], out) 
					colnames(out)[1] <- resp_var
					out[, colnames(out) != "(Intercept)"]
				}
				
				replace_inf <- function (x) {
					replace(x, is.infinite(x), NA)
				}

				df_limma <- df %>% 
					tibble::rownames_to_column(Identifier) %>% 
					tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
					dplyr::left_join(model_info[, c("Sample_ID", master_col_nm)], by = "Sample_ID") %>% 
					dplyr::mutate(!!master_col_nm := factor(.[, master_col_nm], unique(.[, master_col_nm]), levels = c(resp_var, fixed_vars))) %>% 
					dplyr::select(-Sample_ID) %>% 
					dplyr::group_by(!!id, !!master_col) %>% 
					dplyr::mutate(N = row_number()) %>% 
					dplyr::arrange(!!id, !!master_col, N) %>% # ensure the same order when spreading
					dplyr::ungroup(!!master_col) %>% 
					dplyr::mutate(!!master_col := as.character(!!sym(master_col))) %>% 
					tidyr::complete(N, nesting(!!id, !!master_col)) %>% 
					tidyr::unite(key, !!id, N, remove = FALSE) %>% 
					tidyr::spread(!!master_col, log2Ratio) %>% 
					dplyr::select(-key, -N) %>% 
					`colnames<-`(change_names(colnames(.))) %>% 
					dplyr::mutate_all(funs(ifelse(is.na(.), -Inf, .))) %>% 
					dplyr::group_by(!!id) %>% 
					tidyr::nest() %>% 
					dplyr::mutate(data = purrr::map(data, my_model_matrix)) %>% 
					dplyr::mutate(data = purrr::map(data, replace_inf))

				nm <- rep(colnames(df_limma$data[[1]]), each = nrow(df_limma$data[[1]]))

				df_limma <- df_limma %>% 
					dplyr::mutate(flat_data = purrr::map(data, as.vector)) %>% 
					dplyr::mutate(flat_data = purrr::map(flat_data, set_names, nm)) 
					
				df_limma <- do.call("rbind", df_limma$flat_data) %>% 
					`rownames<-`(df_limma[[Identifier]]) %>% 
					`colnames<-`(restore_names(colnames(.)))
					
				nm_fct <- nm %>% 
					restore_names() %>% 
					factor(levels = c(resp_var, fixed_vars, interact_vars))

				design_interac <- model.matrix(~nm_fct) %>% 
					data.frame(check.names = FALSE) %>% 
					`colnames<-`(gsub("^nm_fct", "", colnames(.)))
					
				efit_interac <- df_limma %>% 
					lmFit(design = design_interac) %>% 
					eBayes()
					
				df_pval_interac <- efit_interac$p.value %>% 
					data.frame(check.names = FALSE) %>% 
					dplyr::select(-grep("Intercept", colnames(.))) %>% 
					dplyr::select(-c(fixed_vars)) %>% 
					tibble::rownames_to_column(Identifier)

			} else {
				efit_interac <- NULL
				df_pval_interac <- NULL
			}

			df_pval <- efit$p.value %>% 
				data.frame(check.names = FALSE) %>% 
				dplyr::select(-grep("Intercept", colnames(.))) %>% 
				# dplyr::select(-grep("\\(TMT_Set\\)", names(.))) %>% 
				tibble::rownames_to_column(Identifier)
				
			if(!is.null(df_pval_interac)) df_pval <- df_pval %>% dplyr::left_join(df_pval_interac, by = Identifier)

		} else if(method == "lm") {
			
			# ---------------------------------
			# random effects not yet handled
			# ---------------------------------
			
			my_lm <- function (x) {
				lm(new_formula_fixed, data = x)
			}
			
			df_lm <- df %>% 
				tibble::rownames_to_column(Identifier) %>% 
				tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
				dplyr::left_join(model_info[, c("Sample_ID", master_col_nm)], by = "Sample_ID") %>% 
				dplyr::mutate(!!master_col_nm := factor(.[, master_col_nm], unique(.[, master_col_nm]), levels = c(resp_var, fixed_vars))) %>% 
				dplyr::select(-Sample_ID) %>% 
				dplyr::group_by(!!id, !!master_col) %>% 
				dplyr::mutate(N = row_number()) %>% 
				dplyr::arrange(!!id, !!master_col, N) %>% # ensure the same order when spreading
				dplyr::ungroup(!!master_col) %>% 
				dplyr::mutate(!!master_col := as.character(!!master_col)) %>% 
				tidyr::complete(N, nesting(!!id, !!master_col)) %>% 
				tidyr::unite(key, !!id, N, remove = FALSE) %>% 
				tidyr::spread(!!master_col, log2Ratio) %>% 
				dplyr::select(-key, -N)

			df_lm <- df_lm %>% 
				`colnames<-`(change_names(colnames(.))) %>% 
				dplyr::group_by(!!id) %>% 
				tidyr::nest() %>% 
				dplyr::mutate(model = purrr::map(data, my_lm)) 

			df_pval <- df_lm %>% 
				dplyr::mutate(glance = map(model, broom::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term)) %>% 
				dplyr::select(-c("estimate", "std.error", "statistic")) %>% 
				tidyr::spread(term , p.value) %>% 
				`colnames<-`(restore_names(colnames(.)))
			
		} else if(method %in% c("lme", "lmer")) {
			my_lmer <- function (x) {
				lmerTest::lmer(lmer_formula, data = x)
			}

			df_lmer <- df %>% 
				tibble::rownames_to_column(Identifier) %>% 
				tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
				dplyr::left_join(model_info[, c("Sample_ID", master_col_nm, random_vars)], by = "Sample_ID") %>% # assume random_vars exist...
				dplyr::mutate(!!master_col_nm := factor(.[, master_col_nm], unique(.[, master_col_nm]), levels = c(resp_var, fixed_vars))) %>% 
				dplyr::select(-Sample_ID) %>% 
				dplyr::group_by(!!id) %>% 
				tidyr::nest() %>% 
				dplyr::mutate(model = purrr::map(data, my_lmer))

			# library(broom.mixed)
			df_pval <- df_lmer %>% 
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
				tidyr::spread(term , p.value) %>% 
				`colnames<-`(gsub(master_col_nm, "", names(.)))
		}
				
		
		df_pval <- df_pval %>% 
			tibble::column_to_rownames(Identifier) %>% 
			`names<-`(paste0("pVal (", names(.), "/", resp_var, ")"))
			
			df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
				`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
				bind_cols()
				
			df_pval <- df_pval %>% 
				tibble::rownames_to_column(Identifier) %>% 
				bind_cols(df_padj)
				
			if (nrow(df_pval) == 1) df_pval$rowname <- rownames(df)[1]  # an exception in lmFit with only one entry for testing
			rm(df_padj)

			# may consider to calculate log2Ratio and FC using the original data before the imputation of NA
			resBase <- df %>% 
				dplyr::select(as.character(model_info[model_info[[master_col_nm]] == resp_var, "Sample_ID"])) %>% 
				rowMeans(., na.rm = TRUE) # the mean log2Ratio of base samples

			df_op <- purrr::map(fixed_vars, function (x) {
				df %>% 
					dplyr::select(as.character(model_info[model_info[[master_col_nm]] == x, "Sample_ID"])) %>% 
					rowMeans(., na.rm = TRUE) %>% 
					`-`(resBase)		
				}) %>% 
				data.frame() %>% 
				`names<-`( paste0("log2Ratio (", fixed_vars, "/", resp_var, ")")) 

			df_op <- cbind.data.frame(df_op, 
					df_op %>% 
						sapply(function(x) {ifelse(x > 0, 2^x, -1/(2^x))}) %>% 
						data.frame(check.names = FALSE) %>% 
						`colnames<-`(gsub("log2Ratio", "FC", names(.)))
				)
				
			df_op <- df_op %>% 
				tibble::rownames_to_column(Identifier) %>% 
				left_join(df_pval, by = Identifier)
			
			rm(df_pval)

			df_op <- df_op %>% 
				mutate_at(.vars = grep("^log2Ratio|^FC \\(", names(.)), round, 2) %>% 
				mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				tibble::column_to_rownames(var = Identifier)
	}


	df_op <- lapply(dots, 
			function(formula) 
			{
				if (master_col_nm == "log2Ratio") {
					formula_op <- parse_formula_long(formula, label_scheme_sub)
					Sample_ID <- formula_op$label_scheme_sub_sub$Sample_ID
					model_info <- formula_op$label_scheme_sub_sub
					
					df_op <- longP(df = df[, names(df) %in% Sample_ID], formula, method, model_info)
				} else {
					formula_op <- parse_formula_wide(formula, master_col_nm, label_scheme_sub)
					Sample_ID <- formula_op$label_scheme_sub_sub$Sample_ID
					model_info <- formula_op$label_scheme_sub_sub

					df_op <- wideP(df = df[, names(df) %in% Sample_ID], formula, method, !!master_col, model_info)
					
				}
			}
		) %>% 
		do.call("cbind", .)

}


#' Significance tests
#'
#' \code{proteoSigtest} produces heat maps.  It is a wrapper round the call \code{info_anal(..., anal_type = "Heatmap")}.
#'
#' reads the data from either "\code{~\\Direcotry\\Peptide\\Peptide All.txt}" at \code{id = pep_seq_mod}, 
#'
#' @param id The name of a unique identifier (see \code{\link[proteoQ]{MDS}}).
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
proteoSigtest_master_col <- function (df = NULL, id = gene, scale_log2r = FALSE, filepath = NULL, filename = NULL, 
											impute_na = TRUE, complete_cases = FALSE, method = "limma", master_col = "log2Ratio", ...) {
	id <- enexpr(id)
	
	if (!impute_na) impute_na <- TRUE
	
	# formulas <- enexprs(...)

	info_anal(df, id = !!id, scale_log2r = scale_log2r, filepath = filepath, filename = filename, 
						impute_na = TRUE, 
						anal_type = "Model")(complete_cases, method, master_col, ...)
					
}



#' Perform significance tests
#' 
#' @import limma stringr purrr tidyr dplyr rlang
#' @importFrom magrittr %>% %$% 
#' @importFrom broom.mixed tidy
sigTest_v2 <- function(df, id, master_col, label_scheme_sub, filepath, filename, complete_cases, method, ...) {

	id <- rlang::ensym(id)
	Identifier <- rlang::as_string(id)
	
	# master_col <- rlang::ensym(master_col) 
	# master_col_nm <- master_col%>% rlang::as_string()

	dots = rlang::enexprs(...)
	
	
	
	restore_names <- function(x) {
		x %>% 
			gsub("_PLUS_", "+", .) %>% 
			gsub("_MINUS_", "-", .) %>% 
			gsub("_MULTIPLY_", "*", .) %>% 
			gsub("_DIVIDE_", "/", .)
	}
	
	change_names <- function(x) {
		x %>% 
			gsub("(.*)\\+(.*)", "\\1_PLUS_\\2", .) %>% 
			gsub("(.*)\\-(.*)", "\\1_MINUS_\\2", .) %>% 
			gsub("(.*)\\*(.*)", "\\1_MULTIPLY_\\2", .) %>% 
			gsub("(.*)\\/(.*)", "\\1_DIVIDE_\\2", .)
	}
		
	change_inquotes <- function (x) {
		x %>% 
			gsub("\\\"(.*)\\+(.*)\\\"+?", "\\\"\\1_PLUS_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\-(.*)\\\"+?", "\\\"\\1_MINUS_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\*(.*)\\\"+?", "\\\"\\1_MULTIPLY_\\2\\\"", .) %>% 
			gsub("\\\"(.*)\\/(.*)\\\"+?", "\\\"\\1_DIVIDE_\\2\\\"", .) %>% 
			gsub("\\\"", "", .)
	}
	
	parse_formula_wide <- function (formula, master_col_nm, label_scheme_sub) {

		fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
		fml[1] <- fml[1] %>% change_names()
		fml[2] <- fml[2] %>% change_inquotes()
		
		new_formula <- as.formula(paste(fml[1], fml[2], sep = "~"))
		new_terms = attr(terms.formula(new_formula), "term.labels")
		new_terms_paren <- gsub("(1\\s+\\|\\s+.*)", "\\(\\1\\)", new_terms) # with parenthesis for random effects
		new_formula <- as.formula(paste(fml[1], sprintf(" ~ %s", paste(new_terms_paren, collapse = "+"))))

		# extract elements from the new_formula
		resp_var <- fml[1]		
		
		fixed_terms <- new_terms_paren %>% .[!grepl("\\|", .)] # variables and their interactions
		fixed_vars <- fixed_terms %>% .[!grepl("\\:", .)] # variables
		interact_vars <- fixed_terms %>% .[grepl("\\:", .)]
		new_formula_fixed <- as.formula(paste(resp_var, sprintf(" ~ %s", paste(fixed_terms, collapse = "+")))) # formula

		random_terms <- new_terms_paren %>% .[grepl("\\|", .)] # "1|variable"
		random_vars <- random_terms %>% gsub("\\(1\\s+\\|\\s+(.*)\\)$", "\\1", .) # variables
		random_terms_cmbn <- random_terms %>% paste(collapse = "+")

		# -----------------------------------------------------------
		# design matrix for lmFit
		# interaction terms will be removed for now
		if(purrr::is_empty(random_vars)) {
			lmer_formula <- as.formula(paste("log2Ratio", "~", master_col_nm))
		} else {
			lmer_formula <- as.formula(paste("log2Ratio", "~", master_col_nm, "+", random_terms_cmbn))
		}
		
		label_scheme_sub_sub <- label_scheme_sub %>% 
			dplyr::mutate(!!master_col_nm := change_names(.[, master_col_nm])) %>% 
			# mutate_at(vars(master_col_nm), as.character) %>% 
			dplyr::filter(.[, master_col_nm] %in% c(resp_var, fixed_vars)) 
		
		design <- model.matrix(~0+label_scheme_sub_sub[[master_col_nm]]) %>% 
			data.frame(check.names = FALSE) %>% 
			`colnames<-`(gsub("label_scheme_sub_sub[[master_col_nm]]", "", colnames(.), fixed = TRUE))

		if(purrr::is_empty(random_vars)) {
			limma_formula <- new_formula_fixed
			design_limma <- model.matrix(limma_formula, design)
		} else {
			limma_formula <- as.formula(paste(resp_var, sprintf(" ~ %s + %s", paste(fixed_terms, collapse = "+"), paste(random_vars, collapse = "+"))))
			design_limma <- model.matrix(limma_formula, data.frame(design, label_scheme_sub_sub[, random_vars, drop = FALSE]))
		}
		
		design_limma <- design_limma[, !duplicated(t(design_limma))] 
		ind <- design_limma %>% colSums() %>% `>`(0)
		design_limma <- design_limma[, ind]
		# -----------------------------------------------------------

		# change names back
		label_scheme_sub_sub <- label_scheme_sub_sub %>% 
			dplyr::mutate(!!master_col_nm := restore_names(.[, master_col_nm]))
			
		colnames(design_limma) <- restore_names(colnames(design_limma))
		fixed_vars <- restore_names(fixed_vars)

		# for (i in seq_along(label_scheme_sub_sub)) assign(names(label_scheme_sub_sub[i]), as.vector(label_scheme_sub_sub[, i]))
		# Model_Group <- Model_Group %>% as.character() %>% factor(levels = c(resp_var, fixed_vars))
		# rm(i)

		return(list(new_formula = new_formula, new_formula_fixed = new_formula_fixed, design_limma = design_limma, 
						lmer_formula = lmer_formula, resp_var = resp_var, fixed_vars = fixed_vars, 
						random_terms_cmbn = random_terms_cmbn, random_vars = random_vars, interact_vars = interact_vars, 
						label_scheme_sub_sub = label_scheme_sub_sub))
	}
	
	wideP <- function (df, formula, method, master_col, model_info) {
		master_col <- rlang::ensym(master_col) 
		master_col_nm <- master_col %>% rlang::as_string()

		# formula = V ~ Ner * PD + "Ner+PD" + (1|TMT_Set)
		
		formula_op <- parse_formula_wide(formula, master_col_nm, model_info)

		# model_info <- formula_op$model_info
		design_limma <- formula_op$design_limma
		resp_var <- formula_op$resp_var
		fixed_vars <- formula_op$fixed_vars
		interact_vars <- formula_op$interact_vars
		random_vars <- formula_op$random_vars
		
		new_formula <- formula_op$new_formula 
		new_formula_fixed <- formula_op$new_formula_fixed
		lmer_formula <- formula_op$lmer_formula

		if(method == "limma") {
			
			# (1) without interaction terms
			design_fixed_limma <- design_limma[, colnames(design_limma) %in% c("(Intercept)", fixed_vars)]

			if((ncol(design_limma) - ncol(design_fixed_limma)) > 0) {
				design_random_limma <- design_limma[, !colnames(design_limma) %in% c("(Intercept)", fixed_vars), drop = FALSE] %>% 
					.[, 1, drop = FALSE] # choose the first one being the random effects
					
				corfit <- duplicateCorrelation(df, design = design_fixed_limma, block = design_random_limma)
				efit <- df %>% lmFit(design = design_fixed_limma, block = design_random_limma, correlation = corfit$consensus) %>% eBayes()
			} else {
				efit <- df %>% lmFit(design = design_fixed_limma) %>% eBayes()
			}

			rm(design_fixed_limma)
			
			# interaction terms
			if (!purrr::is_empty(interact_vars)) {
				my_model_matrix <- function (x) {
					out <- model.matrix(new_formula_fixed, x)
					attr(out, "assign") <- NULL
					
					out <- cbind(x[[resp_var]], out) 
					colnames(out)[1] <- resp_var
					out[, colnames(out) != "(Intercept)"]
				}
				
				replace_inf <- function (x) {
					replace(x, is.infinite(x), NA)
				}

				df_limma <- df %>% 
					tibble::rownames_to_column(Identifier) %>% 
					tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
					dplyr::left_join(model_info[, c("Sample_ID", master_col_nm)], by = "Sample_ID") %>% 
					dplyr::mutate(!!master_col_nm := factor(.[, master_col_nm], unique(.[, master_col_nm]), levels = c(resp_var, fixed_vars))) %>% 
					dplyr::select(-Sample_ID) %>% 
					dplyr::group_by(!!id, !!master_col) %>% 
					dplyr::mutate(N = row_number()) %>% 
					dplyr::arrange(!!id, !!master_col, N) %>% # ensure the same order when spreading
					dplyr::ungroup(!!master_col) %>% 
					dplyr::mutate(!!master_col := as.character(!!sym(master_col))) %>% 
					tidyr::complete(N, nesting(!!id, !!master_col)) %>% 
					tidyr::unite(key, !!id, N, remove = FALSE) %>% 
					tidyr::spread(!!master_col, log2Ratio) %>% 
					dplyr::select(-key, -N) %>% 
					`colnames<-`(change_names(colnames(.))) %>% 
					dplyr::mutate_all(funs(ifelse(is.na(.), -Inf, .))) %>% 
					dplyr::group_by(!!id) %>% 
					tidyr::nest() %>% 
					dplyr::mutate(data = purrr::map(data, my_model_matrix)) %>% 
					dplyr::mutate(data = purrr::map(data, replace_inf))

				nm <- rep(colnames(df_limma$data[[1]]), each = nrow(df_limma$data[[1]]))

				df_limma <- df_limma %>% 
					dplyr::mutate(flat_data = purrr::map(data, as.vector)) %>% 
					dplyr::mutate(flat_data = purrr::map(flat_data, set_names, nm)) 
					
				df_limma <- do.call("rbind", df_limma$flat_data) %>% 
					`rownames<-`(df_limma[[Identifier]]) %>% 
					`colnames<-`(restore_names(colnames(.)))
					
				nm_fct <- nm %>% 
					restore_names() %>% 
					factor(levels = c(resp_var, fixed_vars, interact_vars))

				design_interac <- model.matrix(~nm_fct) %>% 
					data.frame(check.names = FALSE) %>% 
					`colnames<-`(gsub("^nm_fct", "", colnames(.)))
					
				efit_interac <- df_limma %>% 
					lmFit(design = design_interac) %>% 
					eBayes()
					
				df_pval_interac <- efit_interac$p.value %>% 
					data.frame(check.names = FALSE) %>% 
					dplyr::select(-grep("Intercept", colnames(.))) %>% 
					dplyr::select(-c(fixed_vars)) %>% 
					tibble::rownames_to_column(Identifier)

			} else {
				efit_interac <- NULL
				df_pval_interac <- NULL
			}

			df_pval <- efit$p.value %>% 
				data.frame(check.names = FALSE) %>% 
				dplyr::select(-grep("Intercept", colnames(.))) %>% 
				# dplyr::select(-grep("\\(TMT_Set\\)", names(.))) %>% 
				tibble::rownames_to_column(Identifier)
				
			if(!is.null(df_pval_interac)) df_pval <- df_pval %>% dplyr::left_join(df_pval_interac, by = Identifier)

		} else if(method == "lm") {
			
			# ---------------------------------
			# random effects not yet handled
			# ---------------------------------
			
			my_lm <- function (x) {
				lm(new_formula_fixed, data = x)
			}
			
			df_lm <- df %>% 
				tibble::rownames_to_column(Identifier) %>% 
				tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
				dplyr::left_join(model_info[, c("Sample_ID", master_col_nm)], by = "Sample_ID") %>% 
				dplyr::mutate(!!master_col_nm := factor(.[, master_col_nm], unique(.[, master_col_nm]), levels = c(resp_var, fixed_vars))) %>% 
				dplyr::select(-Sample_ID) %>% 
				dplyr::group_by(!!id, !!master_col) %>% 
				dplyr::mutate(N = row_number()) %>% 
				dplyr::arrange(!!id, !!master_col, N) %>% # ensure the same order when spreading
				dplyr::ungroup(!!master_col) %>% 
				dplyr::mutate(!!master_col := as.character(!!master_col)) %>% 
				tidyr::complete(N, nesting(!!id, !!master_col)) %>% 
				tidyr::unite(key, !!id, N, remove = FALSE) %>% 
				tidyr::spread(!!master_col, log2Ratio) %>% 
				dplyr::select(-key, -N)

			df_lm <- df_lm %>% 
				`colnames<-`(change_names(colnames(.))) %>% 
				dplyr::group_by(!!id) %>% 
				tidyr::nest() %>% 
				dplyr::mutate(model = purrr::map(data, my_lm)) 

			df_pval <- df_lm %>% 
				dplyr::mutate(glance = map(model, broom::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term)) %>% 
				dplyr::select(-c("estimate", "std.error", "statistic")) %>% 
				tidyr::spread(term , p.value) %>% 
				`colnames<-`(restore_names(colnames(.)))
			
		} else if(method %in% c("lme", "lmer")) {
			my_lmer <- function (x) {
				lmerTest::lmer(lmer_formula, data = x)
			}

			df_lmer <- df %>% 
				tibble::rownames_to_column(Identifier) %>% 
				tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
				dplyr::left_join(model_info[, c("Sample_ID", master_col_nm, random_vars)], by = "Sample_ID") %>% # assume random_vars exist...
				dplyr::mutate(!!master_col_nm := factor(.[, master_col_nm], unique(.[, master_col_nm]), levels = c(resp_var, fixed_vars))) %>% 
				dplyr::select(-Sample_ID) %>% 
				dplyr::group_by(!!id) %>% 
				tidyr::nest() %>% 
				dplyr::mutate(model = purrr::map(data, my_lmer))

			# library(broom.mixed)
			df_pval <- df_lmer %>% 
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
				tidyr::spread(term , p.value) %>% 
				`colnames<-`(gsub(master_col_nm, "", names(.)))
		}
				
		
		df_pval <- df_pval %>% 
			tibble::column_to_rownames(Identifier) %>% 
			`names<-`(paste0("pVal (", names(.), "/", resp_var, ")"))
			
			df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
				`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
				bind_cols()
				
			df_pval <- df_pval %>% 
				tibble::rownames_to_column(Identifier) %>% 
				bind_cols(df_padj)
				
			if (nrow(df_pval) == 1) df_pval$rowname <- rownames(df)[1]  # an exception in lmFit with only one entry for testing
			rm(df_padj)

			# may consider to calculate log2Ratio and FC using the original data before the imputation of NA
			resBase <- df %>% 
				dplyr::select(as.character(model_info[model_info[[master_col_nm]] == resp_var, "Sample_ID"])) %>% 
				rowMeans(., na.rm = TRUE) # the mean log2Ratio of base samples

			df_op <- purrr::map(fixed_vars, function (x) {
				df %>% 
					dplyr::select(as.character(model_info[model_info[[master_col_nm]] == x, "Sample_ID"])) %>% 
					rowMeans(., na.rm = TRUE) %>% 
					`-`(resBase)		
				}) %>% 
				data.frame() %>% 
				`names<-`( paste0("log2Ratio (", fixed_vars, "/", resp_var, ")")) 

			df_op <- cbind.data.frame(df_op, 
					df_op %>% 
						sapply(function(x) {ifelse(x > 0, 2^x, -1/(2^x))}) %>% 
						data.frame(check.names = FALSE) %>% 
						`colnames<-`(gsub("log2Ratio", "FC", names(.)))
				)
				
			df_op <- df_op %>% 
				tibble::rownames_to_column(Identifier) %>% 
				left_join(df_pval, by = Identifier)
			
			rm(df_pval)

			df_op <- df_op %>% 
				mutate_at(.vars = grep("^log2Ratio|^FC \\(", names(.)), round, 2) %>% 
				mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				tibble::column_to_rownames(var = Identifier)
	}

	
	
	
	

	my_lm <- function (f) {
		force(f)
		function (formula, data) try(f(formula, data))
	}
	
	# parse formula with design information 
	# factors under different channels in label_scheme_sub
	# model.matrix to intercept... need revisit
	# consider take contrasts in function input to reduce complexity
	parse_mc_formula <- function (formula, label_scheme_sub) {
		
		# formula = ~ Group * Treatment + (1 | TMT_Set)
		# formula = ~ Type[WT, DB] * Treatment[CTRL, TR] * Time[T1, T2, T3] - Type:Treatment:Time + (1|TMT_Set)
		# formula = log2Ratio ~ Type[WT, DB] * Treatment[CTRL, TR] * Time[T1, T2, T3] - Type:Treatment:Time + (1|TMT_Set)
		
		fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
		len <- length(fml)
		
		# ---------------------------------------------
		interim_formula <- as.formula(paste("log2Ratio", fml[len], sep = "~"))
		interim_terms = attr(terms.formula(interim_formula), "term.labels")
		
		expl_vars <- interim_terms %>% .[!grepl("\\:|\\|", .)] %>% gsub("(.*)\\[.*\\]", "\\1", .)

		# label_scheme to be included for significance tests
		params <- interim_terms %>% 
			.[!grepl("\\:|\\|", .)] %>% 
			gsub("\\]$|\\s+", "", .) %>% 
			stringr::str_split("[\\[|\\]|[,\\s+]]", simplify = TRUE) %>% 
			data.frame() %>%
			dplyr::mutate_if(is.factor, as.character) %>% 
			dplyr::mutate_all(funs(ifelse(nchar(.) == 0, NA, .))) %>% 
			unlist() %>% 
			.[!is.na(.)] %>% 
			.[! . %in% expl_vars]
		# purrr::is_empty(params)	
		
		
		ind <- label_scheme_sub %>% 
			dplyr::select(which(names(.) %in% expl_vars)) %>% 
			dplyr::mutate_all(funs(ifelse(. %in% params, TRUE, NA))) %>% 
			purrr::pmap(list) %>% 
			purrr::map(anyNA) %>% 
			unlist() %>% 
			!.
		
		if(sum(ind) == 0) label_scheme_sub_sub <- label_scheme_sub else label_scheme_sub_sub <- label_scheme_sub[ind, ]
		
		rm(interim_formula, interim_terms)
		rm(params, ind)
		# ---------------------------------------------

		fml[len] <- fml[len] %>% gsub("\\[.*\\]+?", "", .)
		new_formula <- as.formula(paste("log2Ratio", fml[len], sep = "~"))
		new_terms = attr(terms.formula(new_formula), "term.labels")
		new_terms_paren <- gsub("(1\\s+\\|\\s+.*)", "\\(\\1\\)", new_terms) # with parenthesis for random effects
		new_formula <- as.formula(paste("log2Ratio", sprintf(" ~ %s", paste(new_terms_paren, collapse = "+"))))

		fixed_terms <- new_terms_paren %>% .[!grepl("\\|", .)] # variables and their interactions
		fixed_vars <- fixed_terms %>% .[!grepl("\\:", .)] # variables
		interact_vars <- fixed_terms %>% .[grepl("\\:", .)]
		new_formula_fixed <- as.formula(paste("log2Ratio", sprintf(" ~ %s", paste(fixed_terms, collapse = "+")))) # formula

		random_terms <- new_terms_paren %>% .[grepl("\\|", .)] # "1|variable"
		random_vars <- random_terms %>% gsub("\\(1\\s+\\|\\s+(.*)\\)$", "\\1", .) # variables
		random_terms_cmbn <- random_terms %>% paste(collapse = "+")
		
		if(purrr::is_empty(random_vars)) {
			new_limma_formula <- new_formula_fixed
		} else {
			new_limma_formula <-  as.formula(paste("log2Ratio", sprintf(" ~ %s", paste(paste(fixed_terms, collapse = "+"), paste(random_vars, collapse = "+"), sep = "+"))))
		}

		return(list(new_formula = new_formula, new_formula_fixed = new_formula_fixed, new_limma_formula = new_limma_formula, 
						resp_var = "log2Ratio", fixed_vars = fixed_vars, 
						random_terms_cmbn = random_terms_cmbn, random_vars = random_vars, interact_vars = interact_vars, 
						label_scheme_sub_sub = label_scheme_sub_sub))
	}
	
	longP <- function (df, formula, method, model_info) {
		# now df contains the data for all samples...

		df_lm <- df %>% 
			dplyr::mutate(Variance = rowVars(.)) %>% 
			dplyr::mutate(!!id := rownames(df)) %>% 
			dplyr::filter(Variance >= 1E-3) %>% 
			tibble::column_to_rownames(Identifier) %>% 
			dplyr::select(which(names(.) %in% model_info$Sample_ID)) %>% 
			tibble::rownames_to_column(Identifier) %>% 
			tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
			dplyr::left_join(model_info, by = "Sample_ID") %>% 
			dplyr::select(which(not_all_NA(.))) %>% 
			# dplyr::select(-Model_Sample) %>% 
			# `names<-`(gsub("^Model_", "", names(.))) %>% 
			dplyr::group_by(!!id) %>% 
			tidyr::nest()

		if(method == "limma") {
			formula_op <- parse_mc_formula(formula, model_info)
			
			new_formula <- formula_op$new_formula
			new_terms_paren <- formula_op$random_terms_cmbn
			new_limma_formula <- formula_op$new_limma_formula
			random_vars <- formula_op$random_vars
			fixed_vars <- formula_op$fixed_vars

			design_limma <- model.matrix(new_limma_formula, df_lm$data[[1]])
			design_fixed_limma <- design_limma[, !colnames(design_limma) %in% random_vars, drop = FALSE]

			if((ncol(design_limma) - ncol(design_fixed_limma)) > 0) {
				design_random_limma <- design_limma[, colnames(design_limma) %in% random_vars, drop = FALSE] %>% 
					.[, 1, drop = FALSE] # choose the first one being the random effects

				corfit <- duplicateCorrelation(df, design = design_fixed_limma, block = design_random_limma)
				efit <- df %>% lmFit(design = design_fixed_limma, block = design_random_limma, correlation = corfit$consensus) %>% eBayes()
			} else {
				efit <- df %>% lmFit(design = design_fixed_limma) %>% eBayes()
			}

			df_pval <- efit$p.value %>% 
				data.frame(check.names = FALSE) %>% 
				dplyr::select(-grep("Intercept", colnames(.))) %>% 
				# dplyr::select( which(!names(.) %in% random_vars)) %>% 
				tibble::rownames_to_column(Identifier)

		} else if(method == "lm") {
			df_pval <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, my_lm(lm), formula = new_formula)) %>% 
				dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term)) %>% 
				dplyr::select(-c("estimate", "std.error", "statistic")) %>% 
				tidyr::spread(term , p.value)
		} else if(method %in% c("lme", "lmer")) {
			df_pval <- df_lm %>% 
				dplyr::mutate(model = purrr::map(data, my_lm(lmerTest::lmer), formula = new_formula)) %>% 
				dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
				tidyr::unnest(glance, .drop = TRUE) %>% 
				dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
				dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
				tidyr::spread(term , p.value) 
		}
		
		# summarise results
		df_pval <- df_pval %>% 
			tibble::column_to_rownames(Identifier) %>% 
			`names<-`(paste0("pVal (", names(.), ")"))

		df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
			`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
			bind_cols()
			
		df_pval <- df_pval %>% 
			tibble::rownames_to_column(Identifier) %>% 
			bind_cols(df_padj)
			
		if (nrow(df_pval) == 1) df_pval$rowname <- rownames(df)[1]  # an exception in lmFit with only one entry for testing
		rm(df_padj)

		df_pval <- df_pval %>% 
			# mutate_at(.vars = grep("^log2Ratio|^FC \\(", names(.)), round, 2) %>% 
			mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			tibble::column_to_rownames(var = Identifier)
	}
	

	

	
	
	
	
	
	
	
	
	
	rowVars <- function (x, na.rm = TRUE) { # row variance
			sqr <- function(x) x * x
			n <- rowSums(!is.na(x))
			n[n <= 1] = NA
			return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
	}
		
	limma_summary <- function(fit) {
		# summarise p-values
		df_pval <- fit$p.value %>% 
			data.frame(check.names = FALSE) %>% 
			`names<-`(paste0("pVal (", names(.), ")"))

		df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
			`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
			bind_cols()
			
		df_pval <- df_pval %>% 
			tibble::rownames_to_column() %>% 
			bind_cols(df_padj)
			
		if (nrow(df_pval) == 1) df_pval$rowname <- rownames(df)[1]  # an exception in lmFit with only one entry for testing
		rm(df_padj)

		df_pval <- df_pval %>% 
			mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
			tibble::column_to_rownames()
			
		# summarise log2Ratios
		df_FC <- fit$coefficients %>% 
			data.frame(check.names = FALSE) %>% 
			`names<-`(paste0("log2Ratio (", names(.), ")"))

		df_FC <- cbind.data.frame(df_FC, 
				df_FC %>% 
					sapply(function(x) {ifelse(x > 0, 2^x, -1/(2^x))}) %>% 
					data.frame(check.names = FALSE) %>% 
					`colnames<-`(gsub("log2Ratio", "FC", names(.)))
			) %>% 
			dplyr::mutate_at(.vars = grep("^log2Ratio|^FC\\s*\\(", names(.)), round, 2) 
		
		cbind.data.frame(df_FC, df_pval)
	}
	
	# factorial model formula for interaction terms
	# all factors under one channel in label_scheme_sub
	# comparisons defined by contrasts
	# no interaction terms appear in the formula
	limma_onechannel <- function (df, formula, label_scheme_sub, method) {

		# formula = log2Ratio ~ Group["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]  + (1|TMT_Set) + (1|Duplicate)
		# formula = ~ Group["Ner-V", "Ner_PLUS_PD-PD", "(Ner_PLUS_PD-PD)-(Ner-V)"]
		# formula = ~ Group["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|TMT_Set) 
		# formula = ~ Group["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]  + (1|Duplicate) 
		# formula = log2Ratio ~ Group["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "PD-V"]
		# formula = ~ Group["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|TMT_Set)
		# formula = ~ Group["Ner-V", "PD-V", "(Ner_PLUS_PD-V)"] + (1|Duplicate) 
		# formula = ~ Group["Ner-PD", "V-PD", "(Ner_PLUS_PD-PD)"] + (1|Duplicate) 
		# formula = ~ Group["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"] + (1|TMT_Set)
		# formula = log2Ratio ~ Group
		# formula = ~ Group[~V] # no interaction terms
		
		# master_col <- rlang::ensym(master_col) 
		# master_col_nm <- master_col %>% rlang::as_string()

		fml <- as.character(formula) %>% gsub("\\s+", "", .) %>% .[. != "~"]
		len <- length(fml)
		
		master_col_nm <- fml[len] %>% gsub("(.*)\\[\\s*\\~*.*\\].*", "\\1", .)
		master_col <- sym(master_col_nm)

		if (grepl("\\[\\s*\\~\\s*", fml[len])) { # formula = ~ Group[ ~ V]
			base <- fml[len] %>% gsub(".*\\[\\s*\\~\\s*(.*)\\]", "\\1", .)
		} else {
			base <- NULL
		}

		if (!is.null(base)) { # formula = ~ Group[~V]
			new_levels <- label_scheme_sub[[master_col_nm]] %>% levels()
			new_levels <- c(new_levels[new_levels == base], new_levels[new_levels != base])
			label_scheme_sub <- label_scheme_sub %>% 
				dplyr::mutate(!!master_col := factor(!!master_col, levels = new_levels))
		} else if (!grepl("\\[", fml[len])) { # formula = log2Ratio ~ Group
			new_levels <- label_scheme_sub[[master_col_nm]] %>% levels() # leveled by the alphebatic order
		} else {
			new_levels <- NULL # formula = ~ Group["(Ner+Ner_PLUS_PD)/2-V", "Ner_PLUS_PD-V", "Ner-V"]
		}
		
		if (!is.null(new_levels)) {
			contrs <- paste(new_levels[-1], new_levels[1], sep = "-")
			elements <- new_levels
		} else {
			contrs <- fml[len] %>% gsub(".*\\[(.*)\\].*", "\\1", .) %>% 
				gsub("\\\"", "", .) %>% str_split(",\\s*", simplify = TRUE) %>% 
				as.character()

			elements <- fml[len] %>% gsub(".*\\[(.*)\\].*", "\\1", .) %>% 
				gsub("/[0-9]", "", .) %>% 
				gsub("\\\"", "", .) %>% 
				gsub("[\\(\\)]", "", .) %>% 
				str_split("[,\\+\\-]\\s*", simplify = TRUE) %>% 
				as.character() %>% 
				unique()
		}

		label_scheme_sub_sub <- label_scheme_sub %>% 
			dplyr::filter(!!master_col %in% elements) %>% 
			dplyr::mutate(!!master_col := factor(!!master_col))

		design <- model.matrix(~0+label_scheme_sub_sub[[master_col_nm]]) %>% 
			`colnames<-`(levels(label_scheme_sub_sub[[master_col_nm]]))

		cont_matrix <- makeContrasts(contrasts = contrs, levels = data.frame(design))

		random_vars <- fml[len] %>% 
			gsub("\\[.*\\]+?", "", .) %>% 
			paste("~", .) %>% 
			as.formula() %>% 
			terms.formula(.) %>% 
			attr(., "term.labels") %>% 
			.[grepl("\\|", .)] %>% 
			gsub("1\\s*\\|\\s*(.*)", "\\1", .)		
		
		# ------------------------------------------------------------------------
		# limma
		df_sub <- df[, names(df) %in% label_scheme_sub_sub$Sample_ID]
		df_sub <- df_sub[, as.character(label_scheme_sub_sub$Sample_ID)]

		if(length(random_vars) > 0) {
			design_random <- label_scheme_sub_sub[[random_vars[1]]] # choose the first random variable for random effects
			corfit <- duplicateCorrelation(df_sub, design = design, block = design_random) 
			fit <- df_sub %>% lmFit(design = design, block = design_random, correlation = corfit$consensus) %>% contrasts.fit(cont_matrix) %>% eBayes()
		} else {
			fit <- df_sub %>% lmFit(design = design) %>% contrasts.fit(cont_matrix) %>% eBayes()
		}		
		
		df_op <- df_limma <- limma_summary(fit)
		
		# ------------------------------------------------------------------------
		
		
		# ------------------------------------------------------------------------
		if (method %in% c("lmer", "lme", "lm")) {
			df_lm <- df %>% 
				dplyr::mutate(Variance = rowVars(.)) %>% 
				dplyr::mutate(!!id := rownames(df)) %>% 
				dplyr::filter(Variance >= 1E-3) %>% 
				tibble::column_to_rownames(Identifier) %>% 
				dplyr::select(which(names(.) %in% label_scheme_sub_sub$Sample_ID)) %>% 
				tibble::rownames_to_column(Identifier) %>% 
				tidyr::gather(-Identifier, key = Sample_ID, value = log2Ratio) %>% 
				dplyr::left_join(label_scheme_sub_sub[, c("Sample_ID", master_col_nm, random_vars)], by = "Sample_ID") %>% 
				dplyr::select(which(not_all_NA(.))) %>% 
				# dplyr::select(-Model_Sample) %>% 
				# `names<-`(gsub("^Model_", "", names(.))) %>% 
				dplyr::group_by(!!id) %>% 
				tidyr::nest()
				
			fml_rhs <- gsub("\\[.*\\]", "", formula) %>% .[length(.)]
			new_formula <- as.formula(paste("log2Ratio", "~", fml_rhs))
			cont_matrix_lm <- list(Cdn = cont_matrix) %>% `names<-`(master_col_nm)
			
			if(!purrr::is_empty(random_vars)) {
				df_pval <- df_lm %>% 
					dplyr::mutate(model = purrr::map(data, function (x) lmerTest::lmer(data = x, formula = new_formula, contrasts = cont_matrix_lm))) %>% 
					# dplyr::mutate(model = purrr::map(.$data, function (x) lmerTest::lmer(data = x, formula = new_formula, contrasts = cont_matrix_lm))) %>% 
					dplyr::mutate(glance = purrr::map(model, broom.mixed::tidy)) %>% 
					tidyr::unnest(glance, .drop = TRUE) %>% 
					dplyr::filter(!grepl("Intercept", term), effect != "ran_pars") %>% 
					dplyr::select(-c("group", "effect", "estimate", "std.error", "statistic", "df")) %>% 
					dplyr::mutate(term = gsub(master_col_nm, "", term)) %>% 
					tidyr::spread(term , p.value) %>% 
					tibble::column_to_rownames(Identifier) %>% 
					`names<-`(paste0("pVal (", names(.), ")"))

				# -------------
				# temp = df_lm$data[[825]] %>% lmerTest::lmer(data = ., formula = new_formula, contrasts = cont_matrix_lm)
				# -------------
			} else {
				df_pval <- df_lm %>% 
					dplyr::mutate(model = purrr::map(data, function (x) lm(data = x, formula = new_formula, contrasts = cont_matrix_lm))) %>% 
					# dplyr::mutate(model = purrr::map(.$data, function (x) lm(data = x, formula = new_formula, contrasts = cont_matrix_lm))) %>% 
					dplyr::mutate(glance = purrr::map(model, broom::tidy)) %>% 
					tidyr::unnest(glance, .drop = TRUE) %>% 
					dplyr::filter(!grepl("Intercept", term)) %>% 
					dplyr::select(-c("std.error", "estimate", "statistic")) %>% 
					dplyr::mutate(term = gsub(master_col_nm, "", term)) %>% 
					tidyr::spread(term , p.value) %>% 
					tibble::column_to_rownames(Identifier) %>% 
					`names<-`(paste0("pVal (", names(.), ")"))
				
				# -------------
				# plot(-log10(df_pval[,1]) ~ -log10(fit$p.value[, 1]))
				# plot(-log10(df_pval[,2]) ~ -log10(fit$p.value[, 2]))
				# plot(-log10(df_pval[,3]) ~ -log10(fit$p.value[, 2]))
				# Index = 825
				# Index = 3696
				# temp = df_lm$data[[Index]] %>% lm(data = ., formula = new_formula, contrasts = cont_matrix_lm)
				# temp = lm(data =  df_lm$data[[Index]], formula = new_formula, contrasts = cont_matrix_lm)
				# temp = lm(data =  df_lm$data[[Index]] %>% dplyr::filter(Group %in% c("PD", "V")), formula = new_formula)
				# fit$coefficient %>% data.frame(check.names = FALSE) %>% filter(rownames(.) == "CDK1")
				# fit$p.value %>% data.frame(check.names = FALSE) %>% filter(rownames(.) == "CDK1")
				# fit$p.value %>% data.frame(check.names = FALSE) %>% filter(rownames(.) == "PCNA")
				# -------------

			}
			
			df_padj <- purrr::map(df_pval, p.adjust, "BH") %>% 
				`names<-`(gsub("pVal", "adjP", colnames(df_pval))) %>% 
				bind_cols()
				
			df_pval <- df_pval %>% 
				tibble::rownames_to_column(Identifier) %>% 
				bind_cols(df_padj)
				
			df_pval <- df_pval %>% 
				mutate_at(.vars = grep("pVal\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				mutate_at(.vars = grep("adjP\\s+", names(.)), format, scientific = TRUE, digits = 2) %>% 
				tibble::column_to_rownames(var = Identifier)
		}
		
		df_op <- dplyr::bind_cols(df_limma %>% dplyr::select(-grep("pVal|adjP", names(.))), df_pval)
		# df_op2 <- dplyr::bind_cols(df_limma %>% dplyr::select(-grep("pVal|adjP", names(.))), df_pval)
		# ------------------------------------------------------------------------
	
	return(df_op)
	}
	

	df_op <- lapply(dots, 
			function(formula) 
			{
				df_op <- limma_onechannel(df, formula, label_scheme_sub, method)
				
			}
		) %>% 
		do.call("cbind", .)
		
		
		
		
		
		
		
		
		
		
		
		

}


