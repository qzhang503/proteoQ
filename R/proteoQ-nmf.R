#' Plots results from NMF analysis
#'
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMF <- function(df, id, r, nrun, col_group, label_scheme_sub, filepath, filename,
                    complete_cases, xmin = -1, xmax = 1, x_margin = .1, annot_cols = NULL,
										width_consensus = 6, height_consensus = 6, width_coefmap = 6,
										height_coefmap = 6, ...) {

	stopifnot(nrow(label_scheme_sub) > 0)

  sample_ids <- label_scheme_sub$Sample_ID

	col_group <- rlang::enexpr(col_group)
	id <- rlang::as_string(rlang::enexpr(id))
	dots <- rlang::enexprs(...)

	n_color <- 500
	color_breaks <- c(seq(xmin, -x_margin, length = n_color/2)[1:(n_color/2-1)],
	                  seq(-x_margin, x_margin, length = 3),
	                  seq(x_margin, xmax, length = n_color/2)[2:(n_color/2)])

	if (is.null(dots$color))
	  mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
	else
	  mypalette <- eval(dots$color, env = caller_env())

	if (is.null(annot_cols))
	  annotation_col <- NA
	else
	  annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)

	if (is.null(dots$annotation_colors)) {
	  annotation_colors <- setHMColor(annotation_col)
	} else if (is.na(dots$annotation_colors)) {
	  annotation_colors <- NA
	} else {
	  annotation_colors <- eval(dots$annotation_colors, env = caller_env())
	}

	nm_idx <- names(dots) %in% c("mat", "filename", "annotation_col", "color",
	                             "annotation_colors", "breaks")
	dots[nm_idx] <- NULL
	dots$width <- NULL
	dots$height <- NULL

	x_label <- expression("Ratio ("*log[2]*")")

	fn_prx <- gsub("\\..*$", "", filename)
	fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)

	exprs_data <- data.matrix(2^df)

	pData <- label_scheme_sub %>%
		dplyr::filter(!is.na(!!col_group)) %>%
		dplyr::select(Sample_ID, !!col_group) %>%
		dplyr::rename(Group := !!col_group) %>%
		tibble::column_to_rownames("Sample_ID")

	metadata <- data.frame(labelDescription = c("Case/control status"), row.names = c("Group"))
	phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)
	experimentData <- new("MIAME", name = "Pierre Fermat", lab = "", contact = "",
	                      title = "", abstract = "", url = "", other = list(notes = ""))
	exampleSet <- ExpressionSet(assayData = exprs_data, phenoData = phenoData,
	                            experimentData = experimentData)

	# featureNames(exampleSet)[1:5]
	# phenoData(exampleSet)
	# mat <- Biobase::exprs(exampleSet)

	if(!file.exists(file.path(filepath, paste0(fn_prx, ".Rdata")))) {
  	res <- nmf(exampleSet, r, nrun = nrun)
  	save(res, file = file.path(filepath, paste0(fn_prx, ".Rdata")))
	} else {
	  load(file = file.path(filepath, paste0(fn_prx, ".Rdata")))
	}

	x_coef <- coef(res) %>% as.matrix
	x_consensus <- res@consensus %>% as.matrix

	pdf(file.path(filepath, paste0(fn_prx, "_consensus", "-r_", r, ".pdf")), width = width_consensus,
	    height = height_consensus)
	  consensusmap(res, annCol = exampleSet,
	               annColor=list(Type = 'Spectral', basis = 'Set3',consensus = 'YlOrRd:50'),
	               tracks = c("basis:"), main = '', sub = '')
	dev.off()

	pdf(file.path(filepath, paste0(fn_prx, "_coef", "-r_", r, ".pdf")), width = width_coefmap,
	    height = height_coefmap)
		coefmap(res, annColor = list(Type = 'Spectral', basis = 'Set3', consensus = 'YlOrRd:50'),
		        tracks = c("basis:"))
	dev.off()

	# metagene-specific features
	V_hat <- fitted(res)
	s <- extractFeatures(res)
	# s <- featureScore(res)
	# w <- basis(res) # get matrix W
	# h <- coef(res) # get matrix H

	for (i in 1:r) {
		df_sub <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ]

		nrow <- nrow(df_sub)
		if (nrow > 0) {
			fn_sub <- paste0(fn_prx, "_metagene_", r, "-", i, ".", fn_suffix)
			width <- ncol(df_sub) * 2 + 2

			if (nrow > 300) {
				height <- width * 1.5
				dots$show_rownames = FALSE
			} else {
				cellheight <- 5
				height <- cellheight * nrow + 8
				fontsize_row <- 5

				dots$cellheight <- cellheight
				dots$fontsize_row <- fontsize_row
			}

			p <- my_pheatmap(
				mat = df_sub,
				filename = file.path(filepath, fn_sub),
				annotation_col = annotation_col,
				color = mypalette,
				annotation_colors = annotation_colors,
				breaks = color_breaks,
				!!!dots
			)

			df_op <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ] %>%
				tibble::rownames_to_column(id)

			write.csv(df_op, file.path(filepath, paste0(fn_prx, "__metagene_", r, "-", i, ".csv")),
			          row.names = FALSE)
		}
	}
}


#'NMF analysis
#'
#'Analyzes and visualizes the NMF clustering of peptide or protein
#'\code{log2-ratios}
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}
#'
#'@param id Character string to indicate the type of data. Peptide data will be
#'  used at \code{id = pep_seq} or \code{pep_seq_mod}, and protein data at
#'  \code{id = prot_acc} or \code{gene}.
#'@param  col_select Character string to a column key in \code{expt_smry.xlsx}.
#'  The default key is \code{Select}. Samples corresponding to non-empty entries
#'  under \code{col_select} will be included in the indicated analysis.
#'@param impute_na Logical; if TRUE, imputes missing values.
#'@param complete_cases Logical; if TRUE, only cases that are complete with no
#'  missing values will be used for visualization.
#'@param r Numeric; the factorization rank (\code{\link[NMF]{nmf}}).
#'@param nrun Numeric; the number of runs to perform (\code{\link[NMF]{nmf}}).
#'@param scale_log2r Logical; if TRUE, adjusts \code{log2-ratios} to the same
#'  scale of standard deviation for all samples.
#'@param df The filename of input data. By default, it will be determined by the
#'  value of \code{id}.
#'@param filepath The filepath to output results. By default, it will be
#'  determined by the names of the current functional \code{call}.
#'@param filename A representative filename to output images. By default, it
#'  will be determined by the names of the current \code{call}.
#'@param ... Additional parameters for plotting: \cr \code{xmin}, the minimum x; \cr
#'  \code{xmax}, the maximum x; \cr \code{x_margin}, the margin in heat
#'  scales;\cr \code{annot_cols}, the column keys in \code{expt_smry.xlsx} for
#'  use in the color coding of samples;\cr \code{width_}, the width of plot; \cr
#'  \code{height_}, the height of plot; \cr additional arguments inherited from
#'  \code{\link[pheatmap]{pheatmap}}.
#'@return NMF plots.
#'@import NMF dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoNMF <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), col_select = NULL,
											scale_log2r = FALSE, impute_na = TRUE, complete_cases = FALSE,
											df = NULL, filepath = NULL, filename = NULL, r = 4, nrun = 200,
											xmin = -1, xmax = 1, x_margin = 0.1, annot_cols = NULL,
											width_consensus = 6, height_consensus = 6,
											width_coefmap = 6, height_coefmap = 6, ...) {

	id <- rlang::enexpr(id)
	if(length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))

	col_select <- rlang::enexpr(col_select)

	if(!impute_na) complete_cases <- TRUE

	info_anal(id = !!id, col_select = !!col_select, scale_log2r = scale_log2r, impute_na = impute_na,
					df = df, filepath = filepath, filename = filename,
					anal_type = "NMF")(r = r, nrun = nrun, complete_cases = complete_cases,
					xmin = xmin, xmax = xmax, x_margin = x_margin, annot_cols = annot_cols,
					width_consensus = width_consensus, height_consensus = height_consensus,
					width_coefmap = width_coefmap, height_coefmap = height_coefmap, ...)
}


#'Peptide NMF Classification
#'
#'NMF Classification of peptide \code{log2-ratios}
#'
#' @examples
#' pepNMF(
#'   scale_log2r = TRUE,
#'   xmin = -1,
#'   xmax = 1,
#'   x_margin = 0.1,
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   r = 6,
#'   nrun = 200
#' )
#'
#'@seealso \code{\link{proteoNMF}} for parameters
#'@export
pepNMF <- function (...) {
	proteoNMF(id = pep_seq, ...)
}

#'Protein NMF Classification
#'
#'NMF Classification of protein \code{log2-ratios}
#'
#' @examples
#' prnNMF(
#'   scale_log2r = TRUE,
#'   xmin = -1,
#'   xmax = 1,
#'   x_margin = 0.1,
#'   annot_cols = c("Peptide_Yield", "TMT_Set", "Group"),
#'   r = 6,
#'   nrun = 200
#' )
#'
#'
#'@seealso \code{\link{proteoNMF}} for parameters
#'@export
prnNMF <- function (...) {
	proteoNMF(id = gene, ...)
}

