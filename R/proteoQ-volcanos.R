#' Significance tests
#'
#' \code{proteoSigtest} produces heat maps.  It is a wrapper round the call
#' \code{info_anal(..., anal_type = "Heatmap")}.
#'
#' reads the data from either "\code{~\\Direcotry\\Peptide\\Peptide All.txt}" at
#' \code{id = pep_seq_mod},
#'
#' @param id The name of a unique id (see \code{\link[proteoQ]{MDS}}).
#' @param scale_log2r Logical; if TRUE, rescales \code{log2-ratios} to the same
#'   scale of standard deviation for all samples.
#' @return Images stored under the file folders that are associated to
#'   \code{id}, \code{anal_type}.
#'
#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
#' @export
proteoVolcano <- function (id = "gene", anal_type = "Volcano", df = NULL, scale_log2r = TRUE,
													filepath = NULL, filename = NULL, impute_na = TRUE, adjP = FALSE,
													show_labels = TRUE, use_gagep = TRUE, pval_cutoff = 5E-2, show_sig = NULL, ...) {
	
  options(scipen=999)
  
  old_opt <- options(max.print = 99999)
	on.exit(options(old_opt), add = TRUE)

	old_dir <- getwd()
	on.exit(setwd(old_dir), add = TRUE)

	err_msg_1 <- "Unrecognized 'id'; needs to be \"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\" or \"term\""
	err_msg_2 <- "Unrecognized 'anal_type'; needs to be \"Volcano\" or \"GSVA\""
	err_msg_3 <- "Volcano plots of peptides not available for GSVA."
	err_msg_4 <- "GSVA results not found. Please perform prnGSVA() first."

	scale_log2r <- match_logi_gv("scale_log2r")

	id <- rlang::as_string(rlang::enexpr(id))
	cat(paste0("id = \"", id, "\"", " by the current call\n"))
	id <- match_identifier(id)
	
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)	

	if(id %in% c("prot_acc", "gene")) {
	  cat(paste0("id = \"", id, "\"", " after parameter matching to normPrn()\n"))
	} else if(id %in% c("pep_seq", "pep_seq_mod")) {
	  cat(paste0("id = \"", id, "\"", " after parameter matching to normPep()\n"))
	} else if(id == "term") {
		cat("Enrichment terms being \'id\'\n")
	} else {
		stop(err_msg_1, call. = TRUE)
	}

	if(!anal_type %in% c("Volcano", "GSVA", "mapGAGE", "mapGSVA")) stop(err_msg_2, call. = TRUE)

	if (id %in% c("prot_acc", "gene", "term")) {
		data_type <- "Protein"
	} else if (id %in% c("pep_seq", "pep_seq_mod")) {
		data_type <- "Peptide"
	}

	if(is.null(filepath)) {
		filepath = file.path(dat_dir, data_type, anal_type)
		dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
	}

	if(is.null(filename)) {
		fn_prx <- paste(data_type, anal_type, sep = "_")
		fn_suffix <- "png"
	} else {
		fn_prx <- gsub("\\..*$", "", filename)
		fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
	}
	filename <- paste0(fn_prx, ".", fn_suffix)

	if(is.null(df)) {
		if(anal_type %in% c("Volcano", "ESGAGE", "mapGAGE", "mapGSVA")) {
			if(id %in% c("pep_seq", "pep_seq_mod")) {
				fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
				fn_imp <- file.path(dat_dir, "Peptide", "Peptide_impNA.txt")
				fn_raw <- file.path(dat_dir, "Peptide", "Peptide.txt")

				if(file.exists(fn_p)) {
					src_path <- fn_p
				} else {
					src_path <- ifelse(impute_na, fn_imp, fn_raw)
				}
			} else if(id %in% c("prot_acc", "gene")) {
				fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
				fn_imp <- file.path(dat_dir, "Protein", "Protein_impNA.txt")
				fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")

				if(file.exists(fn_p)) {
					src_path <- fn_p
				} else {
					src_path <- ifelse(impute_na, fn_imp, fn_raw)
				}
			}

		} else if(anal_type %in% c("GSVA")) {
			if(id %in% c("pep_seq", "pep_seq_mod")) {
				stop(err_msg_3, call. = TRUE)
			} else if(id %in% c("prot_acc", "gene", "term")) {
				fn_p <- file.path(dat_dir, "Protein\\GSVA", "Protein_GSVA_pVals.txt")
				if(file.exists(fn_p)) src_path <- fn_p else stop(err_msg_4)
			}
		}

		df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE,
		                        sep = "\t", comment.char = "#"), error = function(e) NA)
		if(!is.null(dim(df))) {
			message(paste("File loaded:", gsub("\\\\", "/", src_path)))
		} else {
			stop(paste("No such file or directory:", gsub("\\\\", "/", src_path)))
		}

		df <- df %>%
			`rownames<-`(.[, id]) %>%
			dplyr::select(-grep("I[0-9]{3}|log2_R[0-9]{3}|^FC", names(.)))
	} else {
	  if (id %in% c("pep_seq", "pep_seq_mod")) {
	    fn_raw <- file.path(dat_dir, "Peptide\\Model", df)
	  } else if (id %in% c("prot_acc", "gene")) {
	    fn_raw <- file.path(dat_dir, "Protein\\Model", df)
	  }
	  
	  df <- tryCatch(read.csv(fn_raw, check.names = FALSE, header = TRUE, sep = "\t",
	                          comment.char = "#"), error = function(e) NA)
	  
	  if(!is.null(dim(df))) {
	    message(paste("File loaded:", gsub("\\\\", "/", fn_raw)))
	  } else {
	    stop(paste("Non-existed file or directory:", gsub("\\\\", "/", fn_raw)))
	  }
	}
	
	# remove white space before NA in scientific notation
	df <- df %>% 
	  dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.character) %>% 
	  dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ gsub("\\s*", "", .x) ) %>% 
	  dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.numeric)

	load(file = file.path(dat_dir, "label_scheme.Rdata"))

	# check the exisitence of "kinase" column
	if(!any(names(df) == "kin_attr")) {
		cat("Columns of kinase annoation not found.\n")
		cat("Set 'annot_kinases = FALSE'.\n")
		annot_kinases <- FALSE
	}	else {
		cat("Columns of kinase annoation found.\n")
		cat("Set 'annot_kinases = TRUE'.\n")
		annot_kinases <- TRUE
	}

	if(!adjP) {
		df <- df %>%
			dplyr::select(-contains("adjP"))
	} else {
		df <- df %>%
			dplyr::select(-contains("pVal")) %>%
			`names<-`(gsub("adjP", "pVal", names(.)))
	}

	plotVolcano(df, !!id, filepath, filename, adjP, show_labels, anal_type,
	            use_gagep, pval_cutoff, show_sig, ...)

	if(annot_kinases & anal_type == "Volcano") {
		dir.create(file.path(filepath, "Kinases"), recursive = TRUE, showWarnings = FALSE)
		df_kinase <- df %>% dplyr::filter(kin_attr)

		plotVolcano(df = df_kinase,
			id = !!id, filepath = file.path(filepath, "Kinases"),
			filename = paste0(fn_prx, "_Kinases.", fn_suffix),
			adjP, show_labels = show_labels, anal_type = anal_type, use_gagep, show_sig, ...)
	}

}


#' Perform significance tests
#'
#' @import limma stringr purrr dplyr rlang grid gridExtra gtable
#' @importFrom magrittr %>%
plotVolcano <- function(df = NULL, id = "gene", filepath = NULL, filename = NULL,
												adjP = FALSE, show_labels = TRUE, anal_type, use_gagep, pval_cutoff, show_sig, ...) {

	fml_volcano <- function (df, formula, ...) {
		ind <- purrr::map(formulas, ~ grepl(., names(df))) %>%
			dplyr::bind_cols() %>%
			rowSums() %>%
			`>`(0)

		df <- df %>%
			dplyr::select(grep(formula, names(.), fixed = TRUE)) %>%
			`colnames<-`(gsub(paste0(formula, "."), "", names(.))) %>%
			dplyr::bind_cols(df[, !ind, drop = FALSE], .)

		plotVolcano_sub(df = df, id = !!id, filepath = file.path(filepath, formula),
		                filename = filename, adjP = adjP, show_labels = show_labels,
		                anal_type = anal_type)(use_gagep, pval_cutoff, show_sig, ...)
	}

	id <- rlang::as_string(rlang::enexpr(id))
	id <- match_identifier(id)

	dots <- rlang::enexprs(...)

	formulas <- names(df) %>% .[grepl("pVal", .)] %>% gsub("(.*)\\.pVal.*", "\\1", .) %>% unique()
	sub_dirs <- file.path(filepath, formulas)
	purrr::walk(sub_dirs, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

	if (length(formulas) > 0) purrr::walk(formulas, fml_volcano, df = df, !!!dots)
}


#' Perform significance tests
#'
#' @import limma stringr purrr dplyr rlang grid gridExtra gtable
#' @importFrom magrittr %>%
plotVolcano_sub <- function(df = NULL, id = "gene", filepath = NULL, filename = NULL,
                            adjP = FALSE, show_labels = TRUE, anal_type) {

  id <- rlang::as_string(rlang::enexpr(id))

	volcano_theme <- theme_bw() +
	  theme(
		axis.text.x  = element_text(angle=0, vjust=0.5, size=24),
		axis.ticks.x  = element_blank(),
		axis.text.y  = element_text(angle=0, vjust=0.5, size=24),
		axis.title.x = element_text(colour="black", size=24),
		axis.title.y = element_text(colour="black", size=24),
		plot.title = element_text(face="bold", colour="black", size=20, hjust=.5, vjust=.5),

		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),

		strip.text.x = element_text(size = 16, colour = "black", angle = 0),
		strip.text.y = element_text(size = 16, colour = "black", angle = 90),

		legend.key = element_rect(colour = NA, fill = 'transparent'),
		legend.background = element_rect(colour = NA,  fill = "transparent"),
		legend.position = "none",
		legend.title = element_text(colour="black", size=18),
		legend.text = element_text(colour="black", size=18),
		legend.text.align = 0,
		legend.box = NULL
	 )

	contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
	  gsub("^log2Ratio\\s+\\(|\\)$", "", .)

	if(anal_type %in% c("Volcano", "GSVA")) { # overall volcano
		function(...) {
			on.exit(message("Volcano plots --- Completed."), add = TRUE)

			# "id" only being used for drawing Venn
			fullVolcano(df = df, id = !!id, contrast_groups = contrast_groups,
				volcano_theme = volcano_theme, filepath = filepath, filename = filename,
				adjP = adjP, show_labels = show_labels, ...)
		}
	} else if(anal_type == "mapGSVA") {
		function(use_gagep = FALSE, pval_cutoff = 1E-6, show_sig = NULL, ...) {
			gsea_res <- prep_gsva(contrast_groups, pval_cutoff)
			gsets <- c(dbs$go_sets, dbs$kegg_sets, dbs$c2_msig)

			gsVolcano(df, contrast_groups, gsea_res, gsea_key = "term", gsets, volcano_theme,
				filepath, filename, adjP, show_labels, show_sig, pval_cutoff, ...)
		}
	} else if(anal_type == "mapGAGE") {
		function(use_gagep = FALSE, pval_cutoff = 1E-6, show_sig = NULL, ...) {
			gsea_res <- prep_gage(contrast_groups, key = "term", pval_cutoff, use_gagep)
			gsets <- c(dbs$go_sets, dbs$kegg_sets, dbs$c2_msig)

			gsVolcano(df, contrast_groups, gsea_res, gsea_key = "term", gsets, volcano_theme,
				filepath, filename, adjP, show_labels, show_sig, pval_cutoff, ...)
		}
	}

}


#' Plot all species
#'
#' @import dplyr purrr rlang ggplot2
#' @importFrom magrittr %>%
fullVolcano <- function(df, id = "gene", contrast_groups, volcano_theme = volcano_theme,
                        filepath = NULL, filename = NULL, adjP = FALSE, show_labels = TRUE, ...) {

  # options(scipen = 999)
  
  id <- rlang::as_string(rlang::enexpr(id))

	dots <- rlang::enexprs(...)
	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)

	x_label <- expression("Ratio ("*log[2]*")")

	stopifnot(length(contrast_groups) > 0)

	dfw <- do.call(rbind, purrr::map(contrast_groups, ~ {
			df[, grepl(paste0(" (", .x, ")"), names(df), fixed = TRUE)] %>%
				`colnames<-`(gsub("\\s+\\(.*\\)$", "", names(.))) %>%
				mutate(Contrast = .x) %>%
				bind_cols(df[, !grepl("^pVal\\s+|^adjP\\s+|^log2Ratio\\s+", names(df)), drop = FALSE], .)
		} )) %>% 
	  dplyr::mutate(
	    Contrast = factor(Contrast, levels = contrast_groups),
	    # pVal = as.numeric(as.character(.$pVal)), 
	    valence = ifelse(.$log2Ratio > 0, "pos", "neg")
	  ) %>%
	  dplyr::filter(!is.na(pVal))
	  
	dfw_sub <- dfw %>%
		dplyr::filter(pVal < yco & abs(log2Ratio) > log2(xco)) %>%
		dplyr::arrange(Contrast, pVal) %>%
		dplyr::group_by(Contrast) %>%
		dplyr::mutate(Index = row_number()) %>%
		data.frame (check.names = FALSE)

	dfw_sub_top20 <- dfw_sub %>%
		dplyr::group_by(Contrast) %>%
		dplyr::top_n(n = -20, wt = pVal) %>%
		data.frame (check.names = FALSE)

	# require(gridExtra)
	# tt <- gridExtra::ttheme_minimal(core = list(fg_params=list(cex = .7)), colhead = list(fg_params=list(cex = .7), parse=TRUE), rowhead = list(fg_params=list(cex = .7)))
	# nrow <- max(min(nrow(dfw_sub),20), 1)
	# tbl <- tableGrob(dfw_sub_top20[,c("Contrast", "Index","gene")], rows=NULL, col=NULL, theme=tt)
	# tbl$heights <- tbl$heights*.6

	# data table for labels
	dt <- purrr::map(contrast_groups, ~ {
			to_csv_(dfw_sub_top20 %>%
									dplyr::filter(Contrast == .x) %>%
									dplyr::select(c("Index", id))) %>%
			{if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .} # the exception of zero entry
		})  %>%
		do.call(rbind, .) %>%
		data.frame(Contrast = contrast_groups, id = ., stringsAsFactors = FALSE) %>%
		dplyr::rename(!!rlang::sym(id) := id) %>%
		dplyr::mutate(Contrast = factor(Contrast,  levels = contrast_groups))

	fn_prx <- gsub("\\..*$", "", filename)

	myPalette <- c("#377EB8", "#E41A1C")
	xmax <- ceiling(pmax(abs(min(dfw$log2Ratio, na.rm = TRUE)), max(dfw$log2Ratio, na.rm = TRUE)))
	ymax <- ceiling(max(-log10(dfw$pVal), na.rm = TRUE)) * 1.1
	dfw_greater <- dfw_sub[dfw_sub$valence == "pos", ]
	dfw_less <- dfw_sub[dfw_sub$valence == "neg", ]

	if(is.null(dots$nrow)) {
		nrow <- ifelse(length(unique(dfw$Contrast)) > 3, 2, 1)
	} else {
		nrow <- dots$nrow
	}

	p <-ggplot() +
		geom_point(data = dfw, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, colour = "gray", shape = 20, alpha = .5) +
		geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[2], shape = 20, alpha = .8) +
		geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[1], shape = 20, alpha = .8) +
		geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = log2(xco), linetype = "longdash", size = .5) +
		scale_x_continuous(limits = c(-xmax, xmax)) +
		scale_y_continuous(limits = c(0, ymax)) +
		labs(title = "", x = x_label, y = expression("P-value ("*-log[10]*")")) +
		# scale_colour_manual(name="Experimental\nCondition", values=(myPalette), guide = FALSE) +
		# scale_color_gradientn(colours = rainbow(5)) +
		# annotation_custom(tbl, xmin=-xmax*.95, xmax=-xmax*.75, ymin=-Inf, ymax=Inf) +
		volcano_theme

	# gene
		if (show_labels) {
			p <- p + geom_text(data = dfw_sub_top20, mapping = aes(x = log2Ratio, y = -log10(pVal), label = Index, color = Index),
					size = 3, alpha = .5, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05)
			p <- p + facet_wrap(~Contrast, nrow = nrow, labeller = label_value)
			p <- p + geom_table(data = dt, aes(table = !!rlang::sym(id)), x = -xmax*.85, y = ymax/2)
		} else {
			p <- p + facet_wrap(~Contrast, nrow = nrow, labeller = label_value)
		}

		if(is.null(dots$width)) {
			width <- ifelse(nrow > 1, 6*length(unique(dfw$Contrast))/nrow + 1, 6*length(unique(dfw$Contrast))/nrow)
		} else {
			width <- dots$width
		}

		if(is.null(dots$height)) {
			height <- 6*nrow
		} else {
			height <- dots$height
		}

		ggsave(file.path(filepath, filename), p, width = width, height = height, units = "in")


		# Venn diagram
		summ_venn <- function(df, id, contrast_groups) {
			stopifnot(length(contrast_groups) <= 5)

			universe <- sort(unique(df[[id]]))

			Counts <- matrix(0, nrow = length(universe), ncol = length(contrast_groups)) %>%
				`rownames<-`(universe) %>%
				`colnames<-`(contrast_groups)

			for (Group in contrast_groups) {
				ids <- df %>% dplyr::filter(Contrast == Group) %>% dplyr::select(id) %>% unlist()
				Counts[, Group] <- rownames(Counts) %in% ids
			}

			return(Counts)
		}

		plot_venn <- function(Counts, filepath, direction) {
			myPalette <- brewer.pal(n = 9, name = "Set1")
			Width <- 3
			Height <- 3

			fn_prx <- paste0("Venn_", direction)

			png(file.path(filepath, paste0(fn_prx, ".png")), width = Width, height = Height, units = "in", res = 300)
				limma::vennDiagram(limma::vennCounts(Counts), circle.col = myPalette, cex = .5)
			dev.off()
			write.csv(Counts, file.path(filepath, paste0(fn_prx, ".csv")), row.names = TRUE)
		}


		summ_venn(dfw_greater, id, contrast_groups) %>% plot_venn(filepath, "greater")
		summ_venn(dfw_less, id, contrast_groups) %>% plot_venn(filepath, "less")
}


#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
prep_gsva <- function(contrast_groups, pval_cutoff = 1E-2) {
	do.call(rbind,
		lapply(
			list.files(path = file.path(dat_dir, "Protein", "GSVA"), pattern = "_pVals.txt$", , full.names = TRUE),
			readr::read_tsv
		)
	) %>%
	dplyr::select(-grep("\\.FC\\s+\\(", names(.))) %>%
	dplyr::select(-grep("\\.log2Ratio\\s+\\(", names(.))) %>%
	dplyr::select(-grep("\\.adjP\\s+\\(", names(.))) %>%
	`names<-`(gsub(".*\\.pVal\\s+\\((.*)\\)$", "\\1", names(.))) %>%
	tidyr::gather(key = contrast, value = p_val, -term) %>%
	dplyr::filter(p_val < pval_cutoff) %>%
	dplyr::mutate(q_val = p.adjust(.$p_val, "BH")) %>%
	dplyr::mutate(term = factor(term), contrast = factor(contrast, levels = contrast_groups)) %>%
	dplyr::arrange(contrast) %>%
	tidyr::complete(term, contrast)
}


#' @importFrom magrittr %>%
prep_gage <- function(contrast_groups, key = "term", pval_cutoff = 1E-2, use_gagep = TRUE) {
	gsea_res <- do.call(rbind,
		lapply(
			list.files(path = file.path(dat_dir, "Protein", "ESGAGE"), pattern = "_pVals.txt$", , full.names = TRUE),
			readr::read_tsv
		)
	) %>%
	dplyr::select(-stat.mean, -set.size, -direction, -p.geomean, -q_geomean, -q.val) %>%
	dplyr::filter(p.val < pval_cutoff)

	if(use_gagep) {
		gsea_res <- gsea_res %>%
			dplyr::select(-p_geomean)
	} else {
		gsea_res <- gsea_res %>%
			dplyr::select(-p.val) %>%
			dplyr::rename(p.val = p_geomean)
	}

	gsea_res <- gsea_res %>%
		dplyr::rename(p_val = p.val) %>%
		dplyr::mutate(q_val = p.adjust(.$p_val, "BH")) %>%
		dplyr::mutate(term = factor(term), contrast = factor(contrast, levels = contrast_groups)) %>%
		dplyr::arrange(contrast) %>%
		tidyr::complete(term, contrast)
}


#' Volcano Plots of Protein \code{log2-ratios} under Given Gene Sets
#'
#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
gsVolcano <- function(df, contrast_groups, gsea_res, gsea_key = "term", gsets, volcano_theme,
												filepath = NULL, filename = NULL, adjP = FALSE, show_labels = TRUE,
												show_sig = NULL, pval_cutoff = 1E-6, ...) {

	dots <- rlang::enexprs(...)
	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)
	x_label <- expression("Ratio ("*log[2]*")")
	y_label <- ifelse(adjP, expression("pVal ("*-log[10]*")"), expression("adjP ("*-log[10]*")"))

	terms <- gsea_res %>%
		dplyr::filter(p_val <= pval_cutoff) %>%
		dplyr::select(gsea_key) %>%
		unique %>%
		unlist

	# stopifnot(length(terms) > 0)

	# continue with the next file at no terms
	if (length(terms) > 0) {
  	dfw <- do.call(rbind,
  		purrr::map(contrast_groups, ~ {
  			df[, grepl(paste0(" (", .x, ")"), names(df), fixed = TRUE)] %>%
  				`colnames<-`(gsub("\\s+\\(.*\\)$", "", names(.))) %>%
  				mutate(Contrast = .x) %>%
  				bind_cols(df[, !grepl("^pVal\\s+|^adjP\\s+|^log2Ratio\\s+", names(df))], .)
  		} )) %>%
  		dplyr::mutate(Contrast = factor(Contrast, levels = contrast_groups),
  			pVal = as.numeric(pVal),
  			valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>%
  		dplyr::filter(!is.na(pVal))

  	lapply(terms, function(gt) {
  		# gt = terms[6]
  		gsets_sub <- gsets %>% .[names(.) == gt]
  		stopifnot(length(gsets_sub) > 0)

  		fn <- gsub(":", "~", gsub("/", "or", names(gsets_sub)[[1]]), fixed = TRUE)

  		res_sub <- gsea_res[as.character(gsea_res$term) == gt, ] %>% data.frame(check.names = FALSE)

  		dfw_sub <- dfw[as.character(dfw$entrez) %in% gsets_sub[[1]], ]
  		stopifnot(nrow(dfw_sub) > 0)

  		xmax <- ceiling(pmax(abs(min(dfw_sub$log2Ratio)), max(dfw_sub$log2Ratio)))
  		ymax <- ceiling(max(-log10(dfw_sub$pVal))) * 1.1

  		# ensure the same levels between "Levels" and "newLevels"
  		Levels <- levels(dfw_sub$Contrast)

  		dfw_sub <- dfw_sub %>%
  			dplyr::arrange(Contrast, pVal) %>%
  			dplyr::group_by(Contrast) %>%
  			dplyr::mutate(Index = row_number()) %>%
  			data.frame(check.names = FALSE) %>%
  			dplyr::left_join(., res_sub, by = c("Contrast" = "contrast")) %>%
  			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
  			dplyr::arrange(Contrast) %>%
  			dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>%
  			dplyr::mutate(p_val = as.numeric(p_val)) %>%
  			dplyr::mutate(q_val = format(q_val, scientific = TRUE, digits = 2)) %>%
  			dplyr::mutate(q_val = as.numeric(q_val)) %>%
  			# mutate(sig_level = ifelse(.$q.val > 0.05, "n.s.", ifelse(.$q.val > 0.005, "*", "**"))) %>%
  			# mutate(newContrast = paste0(Contrast, " (", sig_level, ")"))
  			dplyr::mutate(newContrast = Contrast)

  		if(!is.null(show_sig)) {
  			if(show_sig == "pVal") {
  				dfw_sub <- dfw_sub %>%
  					dplyr::mutate(newContrast = paste0(Contrast, " (p = ", p_val, ")"))
  			} else if(show_sig == "qVal") {
  				dfw_sub <- dfw_sub %>%
  					dplyr::mutate(newContrast = paste0(Contrast, " (q = ", q_val, ")"))
  			}
  		}

  		newLevels <- unique(dfw_sub$newContrast)

  		dfw_sub <- dfw_sub %>%
  			dplyr::mutate(newContrast = factor(newContrast, levels = newLevels)) %>%
  			dplyr::arrange(newContrast)

  		dfw_sub_top20 <- dfw_sub %>%
  			dplyr::group_by(newContrast) %>%
  			dplyr::top_n(n = -20, wt = pVal)

  		dt <- purrr::map(newLevels, ~ {
  				dfw_sub_top20 %>%
  					dplyr::filter(newContrast == .x) %>%
  					data.frame(check.names = FALSE) %>%
  					dplyr::select(c("Index", "gene")) %>%
  					to_csv_() %>%
  					{if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .}
  			}) %>%
  			do.call(rbind, .) %>%
  			data.frame(newContrast = newLevels, Contrast = Levels, Gene = ., stringsAsFactors = FALSE) %>%
  			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
  			dplyr::mutate(newContrast = factor(newContrast, levels = newLevels))

  		dfw_greater <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio > log2(xco))
  		dfw_less <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio < -log2(xco))

  		dt_pos <- ifelse(nrow(dfw_greater) > nrow(dfw_less), -xmax*.85, xmax*.6) # table position
  		myPalette <- c("#377EB8", "#E41A1C")

  		p <- ggplot() +
  			geom_point(data = dfw_sub, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, colour = "gray", shape = 20, alpha = .5) +
  			geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[2], shape = 20, alpha = .8) +
  			geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[1], shape = 20, alpha = .8) +
  			# geom_text(data = dfw_sub_top20, mapping = aes(x = log2Ratio, y = -log10(pVal), label = Index, color = Index), size = 2, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05) +
  			geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = log2(xco), linetype = "longdash", size =.5) +
  			labs(title = names(gsets_sub), x = x_label, y = y_label) +
  			scale_x_continuous(limits = c(-xmax, xmax)) +
  			scale_y_continuous(limits = c(0, ymax)) +
  			volcano_theme
  		p <- p + facet_wrap(~ newContrast, nrow = 1, labeller = label_value)

  		if(show_labels)
  			p <- p + geom_text(data = dfw_sub_top20, mapping = aes(x = log2Ratio, y = -log10(pVal),
  					label = Index, color = Index), size = 2, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05) +
  				geom_table(data = dt, aes(table = Gene), x = dt_pos, y = ymax/2)

  		if(nchar(fn) > 50) fn <- paste0(str_sub(fn, 1, 50), "...") # to avoid long fns for pdf()
  		Width <- (5*length(unique(dfw_sub$newContrast))+1)
  		Height <- 6
  		ggsave(file.path(filepath, paste0(fn, ".png")), p, width = Width, height = Height, dpi = 300, units = "in")
  		write.csv(dfw_sub, file = file.path(filepath, paste0(fn, ".csv")), row.names = FALSE)
  	})
	}
}




#'Volcano Plots of Protein \code{log2-ratios}
#'
#'@seealso \code{\link{proteoVolcano}} for parameters
#'@export
pepVol <- function (...) {
	proteoVolcano(id = "pep_seq", anal_type = "Volcano", ...)
}


#'Volcano Plots of Protein \code{log2-ratios}
#'@seealso \code{\link{proteoVolcano}} for parameters
#'@export
prnVol <- function (...) {
		proteoVolcano(id = "gene", anal_type = "Volcano", ...)
}


#'Volcano Plots of GSVA Terms
#'
#'@seealso \code{\link{proteoVolcano}} for parameters
#'@export
gsvaVol <- function (...) {
	proteoVolcano(id = "term", anal_type = "GSVA", ...)
}


#'Volcano Plots of Protein \code{log2-ratios} under Given Gene Sets
#'
#'@seealso \code{\link{proteoVolcano}} for parameters
#'@export
gsvaMap <- function (...) {
	proteoVolcano(id = "gene", anal_type = "mapGSVA", ...)
}


#'Volcano Plots of Protein \code{log2-ratios} under Given Gene Sets
#'
#'@seealso \code{\link{proteoVolcano}} for parameters
#'@export
gageMap <- function (...) {
	proteoVolcano(id = "gene", anal_type = "mapGAGE", ...)
}



GeomTable <- ggproto(
	"GeomTable",
	Geom,
	required_aes = c("x", "y",  "table"),
	default_aes = aes(
	widthx = 10,
	widthy = 10,
	rownames = NA
	),
	draw_key = draw_key_blank,

	draw_panel = function(data, panel_scales, coord) {
	if (nrow(data) != 1) {
		stop(
		sprintf(
			"only one table per panel allowed, got %s (%s)",
			nrow(data),
			as.character(data)
		),
		call. = FALSE
		)
	}
	wx = data$widthx / 2
	wy = data$widthy / 2

	corners <-
		data.frame(x = c(data$x - wx, data$x + wx),
				 y = c(data$y - wy, data$y + wy))
	d <- coord$transform(corners, panel_scales)

	table = read.csv(text = data$table, header = TRUE)
	if (!is.na(data$rownames)) {
		rownames(table) <-
		unlist(strsplit(data$rownames, "|", fixed = TRUE))
	}

	x_rng <- range(d$x, na.rm = TRUE)
	y_rng <- range(d$y, na.rm = TRUE)

	vp <-
		viewport(
		x = mean(x_rng),
		y = mean(y_rng),
		width = diff(x_rng),
		height = diff(y_rng),
		just = c("center", "center")
		)

	grob <- tableGrob(table, rows = NULL, col = NULL,
	                  theme = ttheme_minimal(core = list(fg_params=list(cex = .7)),
	                                         colhead = list(fg_params=list(cex = .7), parse=TRUE),
	                                         rowhead = list(fg_params=list(cex = .7))))
	grob$heights <- grob$heights*.6

	## add a line across the header
	# grob <- gtable_add_grob(
	#   grob,
	#   grobs = segmentsGrob(y1 = unit(0, "npc"), gp = gpar(lwd = 2.0)),
	#   t = 1,
	#   b = 1,
	#   l = 1,
	#   r = ncol(d) + 1
	# )
	editGrob(grob, vp = vp, name = paste(grob$name, facet_id()))
	}
)


facet_id <- local({
	i <- 1
	function() {
  	i <<- i + 1
  	i
	}
})


geom_table <- function(mapping = NULL, data = NULL, stat = "identity",
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, ...) {
  layer(geom = GeomTable, mapping = mapping, data = data, stat = stat, position = position,
        show.legend = show.legend, inherit.aes = inherit.aes, params = list(na.rm = na.rm, ...)
  )
}


to_csv_ <- function(x) {
	paste(capture.output(write.csv(x, stdout(), row.names = F)), collapse = "\n")
}














#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
vc_ESGAGE_old <- function(df, contrast_groups, method = "ESGAGE", type = "GO", key = "essentialSets", pval_cutoff = 1E-2, volcano_theme = volcano_theme,
												filepath = NULL, filename = NULL, adjP = FALSE, show_labels = TRUE, use_gagep = FALSE,
												show_sig = NULL, ...) {

	read_gage_old <- function (method = "ESGAGE", type = "GO", esg = FALSE) {
		if(esg) {
			path = file.path(dat_dir, "Protein", method, type, "Esstentials")
			pattern = "esgp\\.csv$"
		} else {
			path = file.path(dat_dir, "Protein", method, type)
			pattern = "(greater)|(less)\\.csv$"
		}

		files <- list.files(path = path, pattern = pattern, full.names = TRUE)
		if(purrr::is_empty(files)) {
			stop("No GAGE results available for volcano plot visualizations.", call. = TRUE)
		} else {
			res <- do.call(rbind,
				lapply(list.files(path = path, pattern = pattern, full.names = TRUE),
				read.csv, check.names = FALSE, header = TRUE, comment.char = "#"))
		}
	}


	prep_gage_old <- function(contrast_groups, type = "GO", key = "essentialSets", pval_cutoff = 1E-2, use_gagep = FALSE) {
		esg_sets <- read_gage_old(method = "ESGAGE", type = type, esg = TRUE) %>%
			dplyr::filter(q.val < pval_cutoff) %>%
			dplyr::mutate(contrast = factor(contrast, levels = contrast_groups)) %>%
			dplyr::arrange(contrast)

		if(nrow(esg_sets) == 0) abort("No gene sets pass significance threshold.")

		terms <- unique(esg_sets[[key]])

		all_sets <- read_gage_old(method = "ESGAGE", type = type, esg = FALSE) %>%
			dplyr::filter(.$term %in% terms) %>%
			dplyr::mutate(term = factor(term) ) %>%
			tidyr::complete(term, contrast) %>%
			dplyr::group_by(term, contrast) %>%
			dplyr::arrange(p.val)

		if(use_gagep) {
			all_sets <- all_sets %>%
				dplyr::summarise(q_val = min(q.val), p_val = min(p.val))
		} else {
			all_sets <- all_sets %>%
				dplyr::summarise(q_val = min(q_val), p_val = min(p_val))
		}

		all_sets <- all_sets %>%
			dplyr::arrange(p_val) %>%
			data.frame(check.names = FALSE)

		gsets <- if(type == "GO") dbs$go_sets else if(type == "KEGG") dbs$kegg_sets

		return(list(esg_sets = esg_sets, all_sets = all_sets, gsets = gsets, terms = terms))
	}


	dots <- rlang::enexprs(...)
	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)

	x_label <- expression("Ratio ("*log[2]*")")

	res <- prep_gage_old(contrast_groups, type, key, pval_cutoff, use_gagep)
	res_esg <- res$esg_sets
	res_all <- res$all_sets
	gsets <- res$gsets
	terms <- res$terms
	rm(res)


	stopifnot(nrow(res_esg) > 0)

	dfw <- do.call(rbind,
		purrr::map(contrast_groups, ~ {
			df[, grepl(paste0(" (", .x, ")"), names(df), fixed = TRUE)] %>%
				`colnames<-`(gsub("\\s+\\(.*\\)$", "", names(.))) %>%
				mutate(Contrast = .x) %>%
				bind_cols(df[, !grepl("^pVal\\s+|^adjP\\s+|^log2Ratio\\s+", names(df))], .)
		} )) %>%
		dplyr::mutate(Contrast = factor(Contrast, levels = contrast_groups),
			pVal = as.numeric(pVal),
			valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>%
		dplyr::filter(!is.na(pVal))

	lapply(terms, function(gt) {
		# gt = terms[4]
		gsets_sub <- gsets %>% .[names(.) == gt]
		fn <- gsub(":", "~", gsub("/", "or", names(gsets_sub)[[1]]), fixed = TRUE)

		res_sub <- res_all[as.character(res_all$term) == gt, ] %>% data.frame(check.names = FALSE)

		stopifnot(length(gsets_sub) > 0)

		dfw_sub <- dfw[as.character(dfw$entrez) %in% gsets_sub[[1]], ]

		stopifnot(nrow(dfw_sub) > 0)

		xmax <- ceiling(pmax(abs(min(dfw_sub$log2Ratio)), max(dfw_sub$log2Ratio)))
		ymax <- ceiling(max(-log10(dfw_sub$pVal))) * 1.1

		# ensure the same levels between "Levels" and "newLevels"
		Levels <- levels(dfw_sub$Contrast)

		dfw_sub <- dfw_sub %>%
			dplyr::arrange(Contrast, pVal) %>%
			dplyr::group_by(Contrast) %>%
			dplyr::mutate(Index = row_number()) %>%
			data.frame(check.names = FALSE) %>%
			dplyr::left_join(., res_sub, by = c("Contrast" = "contrast")) %>%
			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
			dplyr::arrange(Contrast) %>%
			dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>%
			dplyr::mutate(p_val = as.numeric(p_val)) %>%
			dplyr::mutate(q_val = format(q_val, scientific = TRUE, digits = 2)) %>%
			dplyr::mutate(q_val = as.numeric(q_val)) %>%
			# mutate(sig_level = ifelse(.$q.val > 0.05, "n.s.", ifelse(.$q.val > 0.005, "*", "**"))) %>%
			# mutate(newContrast = paste0(Contrast, " (", sig_level, ")"))
			dplyr::mutate(newContrast = Contrast)

		# ----------------------
		if(!is.null(show_sig)) {
			if(show_sig == "pVal") {
				dfw_sub <- dfw_sub %>%
					dplyr::mutate(newContrast = paste0(Contrast, " (p = ", p_val, ")"))
			} else if(show_sig == "qVal") {
				dfw_sub <- dfw_sub %>%
					dplyr::mutate(newContrast = paste0(Contrast, " (q = ", q_val, ")"))
			}
		}

		newLevels <- unique(dfw_sub$newContrast)

		dfw_sub <- dfw_sub %>%
			dplyr::mutate(newContrast = factor(newContrast, levels = newLevels)) %>%
			dplyr::arrange(newContrast)

		dfw_sub_top20 <- dfw_sub %>%
			dplyr::group_by(newContrast) %>%
			dplyr::top_n(n = -20, wt = pVal)

		dt <- purrr::map(newLevels, ~ {
				dfw_sub_top20 %>%
					dplyr::filter(newContrast == .x) %>%
					data.frame(check.names = FALSE) %>%
					dplyr::select(c("Index", "gene")) %>%
					to_csv_() %>%
					{if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .}
			}) %>%
			do.call(rbind, .) %>%
			data.frame(newContrast = newLevels, Contrast = Levels, Gene = ., stringsAsFactors = FALSE) %>%
			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
			dplyr::mutate(newContrast = factor(newContrast, levels = newLevels))

		dfw_greater <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio > log2(xco))
		dfw_less <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio < -log2(xco))

		dt_pos <- ifelse(nrow(dfw_greater) > nrow(dfw_less), -xmax*.85, xmax*.6) # table position
		myPalette <- c("#377EB8", "#E41A1C")

		p <- ggplot() +
			geom_point(data = dfw_sub, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, colour = "gray", shape = 20, alpha = .5) +
			geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[2], shape = 20, alpha = .8) +
			geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[1], shape = 20, alpha = .8) +
			# geom_text(data = rbind(dfw_greater, dfw_less), mapping = aes(x = log2Ratio, y = -log10(pVal), label = Index, color = Index), size = 2, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05) +
			geom_text(data = dfw_sub_top20, mapping = aes(x = log2Ratio, y = -log10(pVal), label = Index, color = Index), size = 2, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05) +
			geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
			geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
			geom_vline(xintercept = log2(xco), linetype = "longdash", size =.5) +
			labs(title = names(gsets_sub), x = x_label, y = expression("P-value ("*-log[10]*")")) +
			scale_x_continuous(limits = c(-xmax, xmax)) +
			scale_y_continuous(limits = c(0, ymax)) +
			# scale_color_gradientn(colours = rainbow(5)) +
			# annotation_custom(tbl, xmin = -xmax*.95, xmax = -xmax*.75, ymin = -Inf, ymax = Inf) +
			volcano_theme
		p <- p + facet_wrap(~ newContrast, nrow = 1, labeller = label_value)
		p <- p + geom_table(data = dt, aes(table = Gene), x = dt_pos, y = ymax/2)

		if(nchar(fn) > 50) fn <- paste0(str_sub(fn, 1, 50), "...") # to avoid long fns for pdf()
		Width <- (5*length(unique(dfw_sub$newContrast))+1)
		Height <- 6
		ggsave(file.path(filepath, paste0(fn, ".png")), p, width = Width, height = Height, dpi = 300, units = "in")
		write.csv(dfw_sub, file = file.path(filepath, paste0(fn, ".csv")), row.names = FALSE)
	} )

}
