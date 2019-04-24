#' Plots results from NMF analysis
#' 
#' @import NMF dplyr rlang Biobase
#' @importFrom magrittr %>%
plotNMF <- function(df, id, r, nrun, col_group, label_scheme_sub, filepath, filename, complete_cases, xmin = -1, xmax = 1, x_margin = .1, 
										annot_cols = NULL, 
										width_consensus = 6, height_consensus = 6, width_coefmap = 6, height_coefmap = 6, ...) {

	stopifnot(nrow(label_scheme_sub) > 0)

	# col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	# col_color <- rlang::enexpr(col_color) # %>% rlang::as_string()
	# col_shape <- rlang::enexpr(col_shape) # %>% rlang::as_string()
	# col_size <- rlang::enexpr(col_size) # %>% rlang::as_string()
	# col_alpha <- rlang::enexpr(col_alpha) # %>% rlang::as_string()

	id <- rlang::as_string(rlang::enexpr(id))
	
	dots <- rlang::enexprs(...)
	
	fn_prx <- gsub("\\..*$", "", filename)
	fn_suffix <- gsub(".*\\.(.*)$", "\\1", filename)
	
	# (3) anti-log2 transformation to generate non-negative ratio values for NMF analysis
	exprs_data <- data.matrix(2^df)
	
	# (4) create an expression object
	# -------------------
	pData <- label_scheme_sub %>%
		dplyr::filter(!is.na(!!col_group)) %>% 
		dplyr::select(Sample_ID, !!col_group) %>% 
		dplyr::rename(Group := !!col_group) %>% 
		tibble::column_to_rownames("Sample_ID") 
	
	sample_ids <- label_scheme_sub$Sample_ID

	metadata <- data.frame(labelDescription = c("Case/control status"), row.names = c("Group"))
	phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)
	experimentData <- new("MIAME", name = "Pierre Fermat", lab = "", contact = "", title = "", abstract = "", url = "", other = list(notes = ""))
	exampleSet <- ExpressionSet(assayData = exprs_data, phenoData = phenoData, experimentData = experimentData)
	# featureNames(exampleSet)[1:5]
	# phenoData(exampleSet)
	# mat <- Biobase::exprs(exampleSet)
	# -------------------

	# (6) NMF fitting at a given r
	res <- nmf(exampleSet, r, nrun = nrun)
	save(res, file = file.path(filepath, paste0(fn_prx, ".Rdata")))
	
	pdf(file.path(filepath, paste0(fn_prx, "_consensus", "-r_", r, ".pdf")), width = width_consensus, height = height_consensus)
		consensusmap(res, annCol = exampleSet, annColor=list(Type = 'Spectral', basis = 'Set3', consensus = 'YlOrRd:50'), tracks = c("basis:"), main = '', sub = '')
	dev.off()
	
	pdf(file.path(filepath, paste0(fn_prx, "_coef", "-r_", r, ".pdf")), width = width_coefmap, height = height_coefmap)
		coefmap(res, annColor = list(Type = 'Spectral', basis = 'Set3', consensus = 'YlOrRd:50'), tracks = c("basis:")) # mixture coefficients
	dev.off()

	# extracting metagene-specific features
	V_hat <- fitted(res) # estimated target matrix
	# s <- featureScore(res)
	s <- extractFeatures(res)	
	# write.csv(V_hat, file.path(dat_dir, "Peptide\\NMF", paste0(filename, ".csv")), row.names=TRUE) # save results
	# w <- basis(res) # get matrix W
	# h <- coef(res) # get matrix H

	# set up color palette
	n_color <- 500	
	color_breaks <- c(seq(xmin, -x_margin, length = n_color/2)[1:(n_color/2-1)], 
										seq(-x_margin, x_margin, length = 3), 
										seq(x_margin, xmax, length = n_color/2)[2:(n_color/2)])
	# mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
	if (is.null(dots$color)) mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color) else mypalette <- eval(dots$color, env = caller_env())
	
	# set up annotation columns
	if (is.null(annot_cols)) annotation_col <- NA else annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids) 

	# set up the color palettes of the annotation columns
	if (is.null(dots$annotation_colors)) {
		annotation_colors <- setHMColor(annotation_col)
	} else if (is.na(dots$annotation_colors)) {
		annotation_colors <- NA
	} else {
		annotation_colors <- eval(dots$annotation_colors, env = caller_env())
	}
	
	# parameter(s) in "dots" disabled 
	nm_idx <- names(dots) %in% c("mat", "filename", "annotation_col", "color", "annotation_colors", "breaks")
	dots[nm_idx] <- NULL
	dots$width <- NULL
	dots$height <- NULL

	# set up additional labels
	# sample_ids <- label_scheme_sub$Sample_ID
	x_label <- expression("Ratio ("*log[2]*")")	

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
				# width = width, # to allow settings in cellheight et al., do not specify width and height
				# height = width, 
				!!!dots
			)

			df_op <- df[rownames(df) %in% rownames(V_hat[s[[i]], ]), ] %>% 
				tibble::rownames_to_column(id)

			write.csv(df_op, file.path(filepath, paste0(fn_prx, "__metagene_", r, "-", i, ".csv")), row.names = FALSE)
		}
	}
	
	# return(res)
}


#' NMF
#'
#' \code{proteoNMF} produces heat maps.
#'
#' reads the data from either "\code{~\\Direcotry\\Peptide\\Peptide All.txt}" at
#' \code{id = pep_seq_mod}, or "\code{~\\Direcotry\\Protein\\Protein All by
#' Accession.txt}" at \code{id = prot_acc} or
#' "\code{~\\Direcotry\\Protein\\Protein All by gene.txt}" at \code{id = gene}.
#'
#' @param id The name of a unique identifier (see \code{\link[proteoQ]{MDS}}).
#' @param scale_log2r Logical; if TRUE, rescales \code{log2-ratios} to the same
#'   scale of standard deviation for all samples.
#' @param annot_kinases Logical; if TRUE, annotates proteins being kinases or
#'   not.
#' @return Images stored under the file folders that are associated to
#'   \code{id}, \code{anal_type} and \code{annot_kinases}.
#'
#' @examples
#' MA(
#' 	id = gene,
#' 	scale_log2r = scale_log2r,
#' 	annot_kinases = annot_kinases,
#' )
#'
#' # or use the form of functional factories
#' 	my_prnMA <- info_anal(df = NULL, id = prot_acc, scale_log2r, filepath = NULL, filename = NULL, annot_kinases, anal_type = "MA")
#' 	my_prnMA()
#'
#' 	my_pepMA <- info_anal(df = NULL, id = pep_seq_mod, scale_log2r, filepath = NULL, filename = NULL, annot_kinases, anal_type = "MA")
#' 	my_pepMA()
#'
#' \dontrun{
#' }
#' @import NMF dplyr rlang ggplot2
#' @importFrom magrittr %>%
#' @export
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


#'NMF Analysis of Peptide \code{log2-ratios}
#'@seealso \code{\link{proteoNMF}} for parameters
#'@export
pepNMF <- function (...) {
	proteoNMF(id = pep_seq, ...)
}

#'NMF Analysis of Protein \code{log2-ratios}
#'@seealso \code{\link{proteoNMF}} for parameters
#'@export
prnNMF <- function (...) {
	proteoNMF(id = gene, ...)
}

