#'Volcano plot visualization
#'
#'\code{proteoVolcano} visualizes the volcano plots of peptide or protein data,
#'or protein subgroups under the same gene sets. Users should avoid call the
#'method directly, but instead use the following wrappers.
#'
#'By default, the value of \code{gset_nms} in \code{gspaMap(...)} will match to
#'the one in the latest call to \code{\link{prnGSPA}}.
#'
#'@inheritParams  proteoEucDist
#'@inheritParams  proteoHM
#'@inheritParams  info_anal
#'@inheritParams  prnGSPA
#'@param adjP Logical; if TRUE, use Benjamini-Hochberg pVals.
#'@param show_labels Logical; if TRUE, shows the labels of top twenty entries.
#'@param show_sig Character string indicating the type of significance values to
#'  be shown on the volcano plots of gene sets. The default is \code{"none"}.
#'  Additional choices are from \code{c("pVal", "qVal")} where \code{pVal} or
#'  \code{qVal} will be shown, respectively, in the facet grid of the plots. The
#'  argument is not used in \code{prnVol} and \code{pepVol}.
#'@param pval_cutoff Numeric value or vector. \code{Gene sets} with enrichment
#'  \code{pVals} less significant than the threshold will be excluded for
#'  volcano plot visualization. The argument is not used in \code{prnVol} and
#'  \code{pepVol}.
#'@param logFC_cutoff Numeric value or vector. \code{Gene sets} with absolute
#'  enrichment \code{log2FC} less than the threshold will be excluded for
#'  volcano plot visualization. The cut-off is in a logarithmic base of 2, not
#'  in a linear scale. The argument is not used in \code{prnVol} and
#'  \code{pepVol}.
#'@param gset_nms Character string or vector containing the name(s) of gene
#'  sets. The argument is not used in \code{prnVol} and \code{pepVol}.
#'@inheritParams proteoSigtest
#'@inheritParams proteoGSPA
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}. \cr \cr Additional parameters for
#'  plotting: \cr \code{xco}, the cut-off lines of fold changes at position
#'  \code{x}; the default is at \eqn{-1.2} and \eqn{+1.2}. \cr \code{yco}, the
#'  cut-off line of \code{pVal} at position \code{y}; the default is \eqn{0.05}.
#'  \cr \code{width}, the width of plot; \cr \code{height}, the height of plot.
#'  \cr \code{nrow}, the number of rows in a plot.
#'
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'
#'@example inst/extdata/examples/fasta_psm.R
#'@example inst/extdata/examples/pepseqmod_min.R
#'@example inst/extdata/examples/normPep_min.R
#'@example inst/extdata/examples/normPrn_min.R
#'@example inst/extdata/examples/imputeNA_examples.R
#'@example inst/extdata/examples/sigtest_min.R
#'@seealso \code{\link{prnGSPA}} for enrichment analysis against gene sets,
#'  \code{\link{gspaMap}} for the visualization of gene sets under volcano
#'  plots.
#'@export
proteoVolcano <- function (id = "gene", anal_type = "Volcano", df = NULL, scale_log2r = TRUE,
                           filepath = NULL, filename = NULL, fml_nms = NULL, 
                           impute_na = FALSE, adjP = FALSE, show_labels = TRUE, 
                           pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), 
                           show_sig = "none", gset_nms = c("go_sets", "kegg_sets"), 
                           ...) {
  
  old_opt <- options(
    scipen = 0, 
    warn = 0, 
    max.print = 99999
  )
  
  options(
    scipen = 999,
    warn = 1,
    max.print = 2000000
  )
  
	on.exit(options(old_opt), add = TRUE)

	stopifnot(rlang::is_logical(scale_log2r))
	stopifnot(rlang::is_logical(impute_na))
	stopifnot(rlang::is_logical(adjP))
	stopifnot(rlang::is_logical(show_labels))
	stopifnot(rlang::is_double(pval_cutoff))
	stopifnot(rlang::is_double(logFC_cutoff))

	err_msg_1 <- "Unrecognized 'id'; needs to be one of \"pep_seq\", \"pep_seq_mod\", \"prot_acc\", \"gene\" or \"term\""
	err_msg_2 <- "Unrecognized 'anal_type'; needs to be one of \"Volcano\" or \"GSPA\""
	err_msg_3 <- "Volcano plots of peptides not available for GSPA."
	err_msg_4 <- "GSPA results not found. Perform prnGSPA() first."
	err_msg_5 <- "Perform pepSig() first."
	err_msg_6 <- "Perform both pepImp() and pepSig(impute_na = TRUE) first."
	err_msg_7 <- "Perform prnSig() first."
	err_msg_8 <- "Perform both prnImp() and prnSig(impute_na = TRUE) first."

	scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

	id <- rlang::as_string(rlang::enexpr(id))
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)	
	show_sig <- rlang::as_string(rlang::enexpr(show_sig))
	
	gset_nms <- gset_nms %>% .[. %in% match_gset_nms(gset_nms)]
	if (is.null(gset_nms)) stop ("Unknown gene sets.")

	if (! id %in% c("prot_acc", "gene", "pep_seq", "pep_seq_mod", "term")) stop(err_msg_1, call. = FALSE)
	if (! anal_type %in% c("Volcano", "GSPA", "mapGSPA")) stop(err_msg_2, call. = FALSE)

	if (id %in% c("prot_acc", "gene", "term")) {
		data_type <- "Protein"
	} else if (id %in% c("pep_seq", "pep_seq_mod")) {
		data_type <- "Peptide"
	}

	if (is.null(filepath)) {
		filepath = file.path(dat_dir, data_type, anal_type)
		dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
	}

	if (is.null(filename)) {
		fn_prefix <- paste(data_type, anal_type, sep = "_")
		fn_suffix <- "png"
	} else {
		fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
		fn_prefix <- gsub("\\.[^.]*$", "", filename)
	}
	filename <- paste0(fn_prefix, ".", fn_suffix)

	if (is.null(df)) {
		if (anal_type %in% c("Volcano", "mapGSPA")) {
			if (id %in% c("pep_seq", "pep_seq_mod")) {
			  fn_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_pVals.txt")
			  fn_imp_p <- file.path(dat_dir, "Peptide\\Model", "Peptide_impNA_pvals.txt")
			  src_path <- ifelse(impute_na, fn_imp_p, fn_p)
			  
			  if (!file.exists(src_path)) {
			    if (!impute_na) stop(err_msg_5, call. = FALSE) else stop(err_msg_6, call. = FALSE) 
			  }
			} else if (id %in% c("prot_acc", "gene")) {
				fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
				fn_imp_p <- file.path(dat_dir, "Protein\\Model", "Protein_impNA_pvals.txt")
				src_path <- ifelse(impute_na, fn_imp_p, fn_p)
				
				if (!file.exists(src_path)) {
				  if (!impute_na) stop(err_msg_7, call. = FALSE) else stop(err_msg_8, call. = FALSE) 
				}
			}
		} else if (anal_type %in% c("GSPA")) {
		  if (id %in% c("pep_seq", "pep_seq_mod")) {
		    stop(err_msg_3, call. = FALSE)
		  } else if (id %in% c("prot_acc", "gene", "term")) {
		    fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
		    fn_imp_p <- file.path(dat_dir, "Protein\\Model", "Protein_impNA_pvals.txt")
		    src_path <- ifelse(impute_na, fn_imp_p, fn_p)

		    if (!file.exists(src_path)) {
		      if (!impute_na) stop(err_msg_7, call. = FALSE) else stop(err_msg_8, call. = FALSE) 
		    }
		  }
		}

		df <- tryCatch(read.csv(src_path, check.names = FALSE, header = TRUE,
		                        sep = "\t", comment.char = "#"), error = function(e) NA)
		
		if (!is.null(dim(df))) {
			message(paste("File loaded:", gsub("\\\\", "/", src_path)))
		} else {
			stop(paste("Non-existed file or directory:", gsub("\\\\", "/", src_path)))
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
	  
	  if (!is.null(dim(df))) {
	    message(paste("File loaded:", gsub("\\\\", "/", fn_raw)))
	  } else {
	    stop(paste("Non-existed file or directory:", gsub("\\\\", "/", fn_raw)))
	  }
	}
	
	df <- df %>% rm_pval_whitespace()

	species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
	if (!is_empty(species)) load_dbs(gset_nms = gset_nms, species = species)
	
	dots <- rlang::enexprs(...)
	fmls <- dots %>% .[grepl("^\\s*~", .)]
	dots <- dots[!names(dots) %in% names(fmls)]
	dots <- concat_fml_dots(fmls, fml_nms, dots)
	
	filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
	arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
	select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
	dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]

	df <- df %>% 
	  filters_in_call(!!!filter_dots) %>% 
	  arrangers_in_call(!!!arrange_dots)

	load(file = file.path(dat_dir, "label_scheme.rda"))

	if (!adjP) {
		df <- df %>%
			dplyr::select(-contains("adjP"))
	} else {
		df <- df %>%
			dplyr::select(-contains("pVal")) %>%
			`names<-`(gsub("adjP", "pVal", names(.)))
	}

	plotVolcano(df, !!id, filepath, filename, adjP, show_labels, anal_type,
							pval_cutoff, logFC_cutoff, show_sig, gset_nms, !!!dots)
}


#' Plot volcanos
#'
#' @import limma stringr purrr dplyr rlang grid gridExtra gtable
#' @importFrom magrittr %>%
plotVolcano <- function(df = NULL, id = "gene", filepath = NULL, filename = NULL,
                        adjP = FALSE, show_labels = TRUE, anal_type = "Volcano", 
												pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), show_sig = "none", 
												gset_nms = c("go_sets", "kegg_sets"), ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  fml_nms <- names(df) %>% .[grepl("pVal", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique()  %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls %>% .[names(.) %in% fml_nms]
  fml_nms <- fml_nms %>% .[map_dbl(., ~ which(.x == names(fmls)))]
  
  if (is_empty(fml_nms)) stop("No formula matched; compare the formula name(s) with those in `prnSig(...)`")
  
  sub_dirs <- file.path(filepath, fml_nms)
  purrr::walk(sub_dirs, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))
  
  col_ind <- purrr::map(fml_nms, ~ grepl(.x, names(df))) %>%
    dplyr::bind_cols() %>%
    rowSums() %>%
    `>`(0)
  
  if (length(fml_nms) > 0) purrr::pwalk(list(fml_nms, pval_cutoff, logFC_cutoff), fml_volcano, 
                                        df = df, col_ind = col_ind, id = !!id, 
                                        filepath, filename, adjP, show_labels, anal_type, 
                                        show_sig, gset_nms, !!!dots)
}


#' Formula specific volcano plots
#'
#' @import purrr dplyr rlang
#' @importFrom magrittr %>%
fml_volcano <- function (fml_nm, pval_cutoff, logFC_cutoff, df, col_ind, id, 
                         filepath, filename, adjP, show_labels, anal_type, 
												 show_sig, gset_nms, ...) {

  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- df %>%
    dplyr::select(grep(fml_nm, names(.), fixed = TRUE)) %>%
    `colnames<-`(gsub(paste0(fml_nm, "."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .)
  
  plotVolcano_sub(df = df, id = !!id, filepath = file.path(filepath, fml_nm), 
                  filename = filename, adjP = adjP, show_labels = show_labels,
                  anal_type = anal_type, gset_nms)(pval_cutoff, logFC_cutoff, show_sig, subdir = fml_nm, ...)
}


#' Plot volcanos
#'
#' @import limma stringr purrr dplyr rlang grid gridExtra gtable
#' @importFrom magrittr %>%
plotVolcano_sub <- function(df = NULL, id = "gene", filepath = NULL, filename = NULL,
                            adjP = FALSE, show_labels = TRUE, anal_type = "Volcano", gset_nms = "go_sets", ...) {

  id <- rlang::as_string(rlang::enexpr(id))

	volcano_theme <- theme_bw() +
	  theme(
		axis.text.x = element_text(angle = 0, vjust = 0.5, size = 24),
		axis.ticks.x = element_blank(),
		axis.text.y = element_text(angle = 0, vjust = 0.5, size = 24),
		axis.title.x = element_text(colour = "black", size = 24),
		axis.title.y = element_text(colour="black", size = 24),
		plot.title = element_text(face = "bold", colour = "black", size = 14, hjust = .5, vjust = .5),

		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),

		strip.text.x = element_text(size = 16, colour = "black", angle = 0),
		strip.text.y = element_text(size = 16, colour = "black", angle = 90),

		legend.key = element_rect(colour = NA, fill = 'transparent'),
		legend.background = element_rect(colour = NA,  fill = "transparent"),
		legend.position = "none",
		legend.title = element_text(colour="black", size = 18),
		legend.text = element_text(colour="black", size = 18),
		legend.text.align = 0,
		legend.box = NULL
	 )

	contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
	  gsub("^log2Ratio\\s+\\(|\\)$", "", .)

	if (anal_type %in% c("Volcano")) {
		function(pval_cutoff = NULL, logFC_cutoff = NULL, show_sig = NULL, subdir = NULL, ...) {
			rm(pval_cutoff, logFC_cutoff, show_sig, subdir)

		  # "id" only being used for drawing Venn
			fullVolcano(df = df, id = !!id, contrast_groups = contrast_groups,
				volcano_theme = volcano_theme, filepath = filepath, filename = filename,
				adjP = adjP, show_labels = show_labels, ...)
		}
	} else if (anal_type == "mapGSPA")
	  function(pval_cutoff = 5E-2, logFC_cutoff = log2(1.2), show_sig = "none", subdir = NULL, ...) {
	    stopifnot(!is.null(subdir))
	    
	    filename <- match_gspa_filename("GSPA", subdir)
	    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% unique()
	    fn_prefix <- gsub("\\.[^.]*$", "", filename) %>% unique()
	    
	    stopifnot(length(fn_suffix) == 1)
	    
	    gsea_res <- do.call(rbind,
	                        purrr::map(
	                          list.files(path = file.path(dat_dir, "Protein\\GSPA", subdir), 
	                                     pattern = paste0("^", fn_prefix, ".", fn_suffix, "$"), 
	                                     full.names = TRUE),
	                          readr::read_csv
	                        )) %>% 
	      dplyr::filter(.$term %in% names(gsets))

	    gsVolcano(df, contrast_groups, gsea_res, gsea_key = "term", gsets, volcano_theme,
	              filepath, filename, adjP, show_labels, show_sig, pval_cutoff, logFC_cutoff, ...)
	  }
}


#' Volcano plots for all proteins or peptides in a data set
#'
#' @import dplyr purrr rlang ggplot2
#' @importFrom magrittr %>%
#' @importFrom limma vennDiagram
fullVolcano <- function(df = NULL, id = "gene", contrast_groups = NULL, volcano_theme = volcano_theme,
                        filepath = NULL, filename = NULL, adjP = FALSE, show_labels = TRUE, ...) {

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

	# tt <- gridExtra::ttheme_minimal(core = list(fg_params=list(cex = .7)), 
	#                                 colhead = list(fg_params=list(cex = .7), parse=TRUE), 
	#                                 rowhead = list(fg_params=list(cex = .7)))
	# nrow <- max(min(nrow(dfw_sub),20), 1)
	# tbl <- tableGrob(dfw_sub_top20[,c("Contrast", "Index","gene")], rows=NULL, col=NULL, theme=tt)
	# tbl$heights <- tbl$heights*.6

	# data table for labels
	dt <- purrr::map(contrast_groups, ~ {
	  to_csv_(dfw_sub_top20 %>%
            dplyr::filter(Contrast == .x) %>%
            dplyr::select(c("Index", id))) %>%
            {if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .} # a zero-entry exception
	  })  %>%
	  do.call(rbind, .) %>%
	  data.frame(Contrast = contrast_groups, id = ., stringsAsFactors = FALSE) %>%
	  dplyr::rename(!!rlang::sym(id) := id) %>%
	  dplyr::mutate(Contrast = factor(Contrast,  levels = contrast_groups))

	fn_prefix <- gsub("\\.[^.]*$", "", filename)

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
		geom_point(data = dfw, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, colour = "gray", shape = 20, alpha = .5) +
		geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, color = myPalette[2], shape = 20, alpha = .8) +
		geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, color = myPalette[1], shape = 20, alpha = .8) +
		geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = log2(xco), linetype = "longdash", size = .5) +
		scale_x_continuous(limits = c(-xmax, xmax)) +
		scale_y_continuous(limits = c(0, ymax)) +
		labs(title = "", x = x_label, y = expression("P-value ("*-log[10]*")")) +
		volcano_theme

	if (show_labels) {
		p <- p + geom_text(data = dfw_sub_top20, 
		                   mapping = aes(x = log2Ratio, y = -log10(pVal), label = Index, color = Index),
		                   size = 3, alpha = .5, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05)
		p <- p + facet_wrap(~ Contrast, nrow = nrow, labeller = label_value)
		p <- p + geom_table(data = dt, aes(table = !!rlang::sym(id)), x = -xmax*.85, y = ymax/2)
	} else {
		p <- p + facet_wrap(~ Contrast, nrow = nrow, labeller = label_value)
	}

	if(is.null(dots$width)) {
		width <- ifelse(nrow > 1, 6*length(unique(dfw$Contrast))/nrow + 1, 
		                6*length(unique(dfw$Contrast)) / nrow)
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

		fn_prefix <- paste0("Venn_", direction)

		png(file.path(filepath, paste0(fn_prefix, ".png")), width = Width, height = Height, units = "in", res = 300)
			limma::vennDiagram(limma::vennCounts(Counts), circle.col = myPalette, cex = .5)
		dev.off()
		write.csv(Counts, file.path(filepath, paste0(fn_prefix, ".csv")), row.names = TRUE)
	}

	summ_venn(dfw_greater, id, contrast_groups) %>% plot_venn(filepath, "greater")
	summ_venn(dfw_less, id, contrast_groups) %>% plot_venn(filepath, "less")
}


#' Volcano plots of protein \code{log2FC} under given gene sets
#'
#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
gsVolcano <- function(df = NULL, contrast_groups = NULL, gsea_res = NULL, gsea_key = "term", gsets = NULL, 
                      volcano_theme = NULL, filepath = NULL, filename = NULL, adjP = FALSE, show_labels = TRUE,
                      show_sig = "none", pval_cutoff = 1E-6, logFC_cutoff = log2(1.2), ...) {

	dots <- rlang::enexprs(...)
	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)
	x_label <- expression("Ratio ("*log[2]*")")
	y_label <- ifelse(adjP, expression("pVal ("*-log[10]*")"), expression("adjP ("*-log[10]*")"))
	
	if (nrow(gsea_res) == 0) {
	  stop("No GSEA terms availbale for volcano plots. Consider relax `pval_cutoff` and/or `logFC_cutoff`.")
	}

	N <- pmin(dplyr::n_distinct(gsea_res[[gsea_key]]), 100)
	terms <- gsea_res %>%
	  dplyr::arrange(p_val) %>% 
	  dplyr::slice(1:(N*length(contrast_groups))) %>% 
	  dplyr::filter(p_val <= pval_cutoff, abs(log2fc) >= logFC_cutoff) %>%
	  dplyr::select(gsea_key) %>%
	  unique %>%
	  unlist
	
	gsea_res <- gsea_res %>% 
	  dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>% 
	  dplyr::mutate(q_val = format(p_val, scientific = TRUE, digits = 2)) %>% 
	  dplyr::mutate(log2fc = round(log2fc, digits = 2))

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
  		gsets_sub <- gsets %>% .[names(.) == gt]
  		stopifnot(length(gsets_sub) > 0)

  		fn <- gsub(":", "~", gsub("/", "or", names(gsets_sub)[[1]]), fixed = TRUE)

  		res_sub <- gsea_res[as.character(gsea_res$term) == gt, ] %>% 
  		  data.frame(check.names = FALSE)

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
  			# dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>%
  			# dplyr::mutate(p_val = as.numeric(p_val)) %>%
  			# dplyr::mutate(q_val = format(q_val, scientific = TRUE, digits = 2)) %>%
  			# dplyr::mutate(q_val = as.numeric(q_val)) %>%
  			# mutate(sig_level = ifelse(.$q.val > 0.05, "n.s.", ifelse(.$q.val > 0.005, "*", "**"))) %>%
  			# mutate(newContrast = paste0(Contrast, " (", sig_level, ")"))
  			dplyr::mutate(newContrast = Contrast)
				
  		if (show_sig != "none") {
  			if (grepl("^p", show_sig)) {
  				dfw_sub <- dfw_sub %>%
  					dplyr::mutate(newContrast = paste0(Contrast, " (p = ", p_val, ")"))
  			} else if (grepl("^q", show_sig)) {
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
  		
  		if(is.null(dots$nrow)) {
  		  nrow <- ifelse(length(unique(dfw_sub$Contrast)) > 3, 2, 1)
  		} else {
  		  nrow <- dots$nrow
  		}

  		p <- ggplot() +
  			geom_point(data = dfw_sub, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, colour = "gray", shape = 20, alpha = .5) +
  			geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[2], shape = 20, alpha = .8) +
  			geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), size = 3, color = myPalette[1], shape = 20, alpha = .8) +
  			geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = log2(xco), linetype = "longdash", size =.5) +
  			labs(title = names(gsets_sub), x = x_label, y = y_label) +
  			scale_x_continuous(limits = c(-xmax, xmax)) +
  			scale_y_continuous(limits = c(0, ymax)) +
  			volcano_theme
  		p <- p + facet_wrap(~ newContrast, nrow = nrow, labeller = label_value)

  		if(show_labels)
  			p <- p + geom_text(data = dfw_sub_top20, mapping = aes(x = log2Ratio, y = -log10(pVal),
  					label = Index, color = Index), size = 2, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05) +
  				geom_table(data = dt, aes(table = Gene), x = dt_pos, y = ymax/2)

  		if (nchar(fn) > 50) fn <- paste0(str_sub(fn, 1, 50), "...") # avoid long fns for pdf()
  		
  		if (is.null(dots$width)) {
  		  width <- ifelse(nrow > 1, 6*length(unique(dfw_sub$newContrast))/nrow + 1, 
  		                  6*length(unique(dfw$Contrast)) / nrow)
  		} else {
  		  width <- dots$width
  		}
  		
  		if (is.null(dots$height)) {
  		  height <- 6*nrow
  		} else {
  		  height <- dots$height
  		}
  		
  		# Width <- (5*length(unique(dfw_sub$newContrast))+1)
  		# Height <- 6
  		ggsave(file.path(filepath, paste0(fn, ".png")), p, width = width, height = height, dpi = 300, units = "in")
  		write.csv(dfw_sub, file = file.path(filepath, paste0(fn, ".csv")), row.names = FALSE)
  	})
	}
}



#'Volcano plots of peptide \code{log2FC}
#'
#'\code{pepVol} is a wrapper of \code{\link{proteoVolcano}} for peptide data
#'
#'@rdname proteoVolcano
#'
#'@import purrr
#'@export
pepVol <- function (...) {
  err_msg <- "Don't call the function with arguments `id` and/or `anal_type`.\n"
  warn_msg <- "Parameters \"pval_cutoff\" and \"show_sig\" are not used\n."
  
  dots <- rlang::enexprs(...)
  if (any(names(dots) %in% c("id", "anal_type"))) stop(err_msg)
  if (any(names(dots) %in% c("pval_cutoff", "show_sig"))) warning(warn_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\Volcano\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_pepid()
  
  quietly_log <- purrr::quietly(proteoVolcano)(id = !!id, anal_type = "Volcano", ...)
  purrr::walk(quietly_log[-1], write, file.path(dat_dir, "Peptide\\Volcano\\log\\pepnVol_log.csv"), 
              append = TRUE)
}


#'Volcano plots of protein \code{log2FC}
#'
#'\code{prnVol} is a wrapper of \code{\link{proteoVolcano}} for protein data
#'
#'@rdname proteoVolcano
#'
#' @examples
#' # ===================================
#' # Volcano plot
#' # ===================================
#' scale_log2r <- TRUE
#' 
#' # all peptides
#' pepVol()
#' 
#' # all proteins
#' prnVol(
#'   show_labels = TRUE,
#'   xco = 1.2,
#'   yco = 0.01,
#' )
#'
#' # kinases and prot_n_pep >= 2
#' prnVol(
#'   show_labels = TRUE,
#'   xco = 1.2,
#'   yco = 0.01,
#'   filter_by_kin = exprs(kin_attr, prot_n_pep >= 2),
#'   filename = "prnvol_kin_npep2.png"
#' )
#'
#'
#' # protein subgroups by gene sets
#' # filtered by proteins with two or more identifying peptides for visualization
#' gspaMap(
#'   pval_cutoff = 5E-3,
#'   logFC_cutoff = log2(1.2),
#'   gset_nms = c("go_sets"),
#'   show_sig = p,
#'   show_labels = TRUE,
#'   yco = 0.01,
#'   filter_by_npep = exprs(prot_n_pep >= 2),
#'   # `filename`(s) will be automated, i.e., by gene-set names
#' )
#'
#' # customized thresholds for the corresponding formulae in `prnSig()`
#' gspaMap(
#'   fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
#'   pval_cutoff = c(5E-2, 5E-2, 1E-10),
#'   logFC_cutoff = log2(1.2),
#'
#'   show_sig = pVal,
#'   show_labels = TRUE,
#'   yco = 0.05,
#'   filter_by_npep = exprs(prot_n_pep >= 2),
#' )
#' 
#'@import purrr
#'@export
prnVol <- function (...) {
  err_msg <- "Don't call the function with arguments `id` and/or `anal_type`.\n"
	warn_msg <- "Parameters \"pval_cutoff\" and \"show_sig\" are not used.\n"
	
  dots <- rlang::enexprs(...)
  if (any(names(dots) %in% c("id", "anal_type"))) stop(err_msg)
  if (any(names(dots) %in% c("pval_cutoff", "show_sig"))) warning(warn_msg)

  dir.create(file.path(dat_dir, "Protein\\Volcano\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_protid()

  quietly_log <- purrr::quietly(proteoVolcano)(id = !!id, anal_type = "Volcano", ...)
  purrr::walk(quietly_log[-1], write, file.path(dat_dir, "Protein\\Volcano\\log\\prnVol_log.csv"), 
              append = TRUE)
}


#'Volcano plots of protein \code{log2FC} under gene sets
#'
#'\code{gspaMap} is a wrapper of \code{\link{proteoVolcano}} for
#'mapping of protein data by gene sets.
#'
#'@rdname proteoVolcano
#'
#'@export
gspaMap <- function (...) {
  err_msg <- "Don't call the function with arguments `id` and/or `anal_type`.\n"
  if (any(names(rlang::enexprs(...)) %in% c("id", "anal_type"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\mapGSPA\\log"), recursive = TRUE, showWarnings = FALSE)
  
  id <- match_normPSM_protid()

  quietly_log <- purrr::quietly(proteoVolcano)(id = !!id, anal_type = "mapGSPA", ...)
  purrr::walk(quietly_log[-1], write, file.path(dat_dir, "Protein\\mapGSPA\\log\\mapGSPA_log.csv"), 
              append = TRUE)
}


GeomTable <- ggproto(
	"GeomTable",
	Geom,
	required_aes = c("x", "y",  "table"),
	default_aes = aes(widthx = 10, widthy = 10, rownames = NA),
	draw_key = draw_key_blank,

	draw_panel = function(data, panel_scales, coord) {
  	if (nrow(data) != 1) {
  		stop(sprintf("only one table per panel allowed, got %s (%s)", 
  		             nrow(data), 
  		             as.character(data)), 
  		     call. = FALSE)
  	}
	  
	  wx = data$widthx / 2
	  wy = data$widthy / 2

  	corners <- data.frame(x = c(data$x - wx, data$x + wx), y = c(data$y - wy, data$y + wy))
  	d <- coord$transform(corners, panel_scales)
  
  	table = read.csv(text = data$table, header = TRUE)
  	if (!is.na(data$rownames)) {
  		rownames(table) <-
  		unlist(strsplit(data$rownames, "|", fixed = TRUE))
  	}
  
  	x_rng <- range(d$x, na.rm = TRUE)
  	y_rng <- range(d$y, na.rm = TRUE)

  	vp <- viewport(x = mean(x_rng), y = mean(y_rng), width = diff(x_rng), height = diff(y_rng), 
  	               just = c("center", "center"))
  	  
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



