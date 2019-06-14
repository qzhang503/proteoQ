
#' Make MA plots
#'
#' @import dplyr rlang ggplot2 RColorBrewer
#' @importFrom boot inv.logit
#' @importFrom magrittr %>%
plotMA <- function (df, col_group, label_scheme_sub, filepath, filename, ...) {
	dots <- rlang::enexprs(...)

	col_group <- rlang::enexpr(col_group)

	fn_prx <- gsub("\\..*$", "", filename)

	xmin <- 10^3.5
	xmax <- 10^6.05
	ymin <- -1
	ymax <- 1
	y_breaks <- 1
	x_label <- expression("Ratio ("*log[2]*")")
	logit_x_label <- expression("Ratio (inversed logit "*log[2]*")")

	my_theme <- theme_bw() + theme(
		axis.text.x  = element_text(angle=0, vjust=0.5, size=20),
		# axis.ticks.x  = element_blank(), # x-axis ticks
		axis.text.y  = element_text(angle=0, vjust=0.5, size=20),
		axis.title.x = element_text(colour="black", size=24),
		axis.title.y = element_text(colour="black", size=24),
		plot.title = element_text(face="bold", colour="black", size=24, hjust=0.5, vjust=0.5),

		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),

		strip.text.x = element_text(size = 12, colour = "black", angle = 0),
		strip.text.y = element_text(size = 12, colour = "black", angle = 90),   # Facet title size
		# strip.background = element_rect(colour="red", fill="#fdd49e"),

		legend.key = element_rect(colour = NA, fill = 'transparent'),
		legend.background = element_rect(colour = NA,  fill = "transparent"),
		# legend.position = c(0.2, 0.3),
		# legend.position = "none",  # to remove all legends
		# legend.title = element_blank(),
		legend.title = element_text(colour="black", size=16),
		legend.text = element_text(colour="black", size=14),
		legend.text.align = 0,
		legend.box = NULL
	)

	my_ma_theme_logit <- theme_bw() + theme(
		axis.text.x  = element_text(angle=0, vjust=0.5, size=14),   # rotate x-axis label and define label size
		# axis.ticks.x  = element_blank(), # x-axis ticks
		axis.text.y  = element_text(angle=0, vjust=0.5, size=14),   # rotate y-axis label and define label size
		axis.title.x = element_text(colour="black", size=16),  # x-axis title size
		axis.title.y = element_text(colour="black", size=16),  # y-axis title size
		# plot.title = element_text(face="bold", colour="black", size=24, hjust=.5, vjust=.5),   # main title size
		plot.title = element_blank(),
		plot.margin = unit(c(0, 0.2, 0, 0), "cm"), # top, right, bottom, left

		strip.text.x = element_text(size = 14, colour = "black", angle = 0),
		strip.text.y = element_text(size = 14, colour = "black", angle = 90),   # facet title size

		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.y = element_blank(),
		panel.grid.minor.y = element_blank(),

		legend.key = element_rect(colour = NA, fill = 'transparent'),
		legend.background = element_rect(colour = NA,  fill = "transparent"),
		# legend.position = c(0.15, 0.78),
		legend.position = "none",  # to remove all legends
		legend.title = element_blank(),
		# legend.title = element_text(colour="black", size=18),
		legend.text = element_text(colour="black", size=18),
		legend.text.align = 0,
		legend.box = NULL
	 )

	# Go through sample groups
	label_scheme_sub <- label_scheme_sub %>%
		dplyr::mutate(MA_Group := !!col_group)

	if(nlevels(label_scheme_sub$MA_Group) > 0) {
		purrr::walk(levels(label_scheme_sub$MA_Group), function(group) {
			label_scheme_sub_sub <- label_scheme_sub[label_scheme_sub$MA_Group == group, ]

			dfI_sub <- df$Intensity %>%
				tibble::rownames_to_column() %>%
				dplyr::select(rowname, as.character(label_scheme_sub_sub$Sample_ID)) %>%
				tibble::column_to_rownames(var = "rowname")

			dfR_sub <- df$log2R %>%
				tibble::rownames_to_column() %>%
				dplyr::select(rowname, as.character(label_scheme_sub_sub$Sample_ID)) %>%
				tibble::column_to_rownames(var = "rowname")

			sample_pairs <- combn(names(dfR_sub), 2)
			lcdn <- ncol(sample_pairs)

			### Obtain experimental quantile
			## (1-2) Bin intensity values with finer grids for unparametric confidence interval calculation
			# https://stackoverflow.com/questions/23901907/create-a-log-sequence-across-multiple-orders-of-magnitude?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
			# Seq_int2 <- c(-Inf, 10^log10(c(seq(2,10, 1) %o% 10^(log10(xmin):log10(xmax)))), Inf)
			# Seq_int2 <- c(-Inf, c(seq(2,10, 10/5) %o% 10^(3:6)), Inf)
			Seq_int <- c(-Inf, c(seq(2, 10, 10/10) %o% 10^(3:4)), c(seq(2, 10, 10/5) %o% 10^(5:6)), Inf)

			dfw_ma_sub <- do.call(rbind, lapply (seq_len(lcdn), function (i) {
				tempR <- dfR_sub %>%
					dplyr::select(sample_pairs[, i]) %>%
					dplyr::mutate(log2Ratio = .[, 2] - .[, 1])

				tempI <- dfI_sub %>%
					dplyr::select(sample_pairs[, i]) %>%
					dplyr::mutate(log10Int = log10(rowMeans(.)))

				cbind.data.frame(tempR[, c("log2Ratio"), drop = FALSE], tempI[, c("log10Int"), drop = FALSE]) %>%
					dplyr::mutate(Type = paste(sample_pairs[2, i], "vs", sample_pairs[1, i]))

				# temp_log2R <- dfw_log2R %>%
				# 	dplyr::select(sample_pairs[, i]) %>%
				# 	dplyr::mutate(log2Ratio = .[, 2] - .[, 1])

				# temp_int <- dfw_int %>%
				# 	dplyr::select(sample_pairs[, i]) %>%
				# 	dplyr::mutate(log10Int = log10(rowMeans(.)))

				# cbind.data.frame(temp_log2R[, c("log2Ratio"), drop = FALSE], temp_int[, c("log10Int"), drop = FALSE]) %>%
				# 	dplyr::mutate(Type = paste(sample_pairs[2, i], "vs", sample_pairs[1, i]))
				}
			)) %>%
			dplyr::mutate(Intensity = 10^log10Int, logit_log2Ratio = boot::inv.logit(log2Ratio), Type = factor(Type)) %>%
			dplyr::filter(!is.na(log10Int)) %>%
			dplyr::mutate(Int_index = cut(Intensity, Seq_int))

			dfw_ma_sub <- dfw_ma_sub %>%
				dplyr::group_by(Int_index) %>%
				dplyr::summarise_at(vars(c("Intensity")), ~ mean(., na.rm = TRUE)) %>%
				dplyr::rename(Mean_Int = Intensity) %>%
				dplyr::right_join(., dfw_ma_sub, by = "Int_index") %>%
				dplyr::select(-c("Int_index")) %>%
				dplyr::mutate(key = paste(Mean_Int, Type, sep = "_"))

			## (1-3) calculate the 2.5 -> 97.5 percentile of log2Ratio for each bin (to be used as y values for connecting lines)
			# https://stackoverflow.com/questions/38177908/split-a-data-frame-column-containing-a-list-into-multiple-columns-using-dplyr-o
			q_logit <- dfw_ma_sub %>%
				## (1) much easier use plyr::ddply
				# ddply(.(Mean_Int, Type), function(x) quantile(x$logit_log2Ratio, probs = c(2.5, 97.5)/100, na.rm = TRUE)) %>% # if use
				## (2) an alternative in using purrr::map
				# split(list(.$Mean_Int, .$Type)) %>% purrr::map(~ quantile(.$logit_log2Ratio, probs = c(2.5, 97.5)/100, na.rm = TRUE, data = .)) %>% purrr::map(as.list) %>% dplyr::bind_rows(.id = "key") %>%
				dplyr::group_by(Mean_Int, Type) %>%
				do(data.frame(as.list(quantile(.$logit_log2Ratio, probs = c(2.5, 97.5)/100, na.rm = TRUE)), check.names = FALSE)) %>%
				dplyr::rename(lwr_exp_logit_r = "2.5%", upr_exp_logit_r = "97.5%") %>%
				tidyr::unite(key, Mean_Int, Type)

			## very descent to use data.table
			# library(data.table)
			# q_logit <- as.data.table(dfw_ma_sub)[, as.list(quantile(logit_log2Ratio, probs = c(2.5, 97.5)/100, na.rm = TRUE)), by = list(Mean_Int, Type)] %>%
			# 	dplyr::rename(lwr_exp_logit_r = "2.5%", upr_exp_logit_r = "97.5%") %>%
			# 	tidyr::unite(key, Mean_Int, Type)

			dfw_ma_sub <- dfw_ma_sub %>%
				dplyr::left_join(q_logit, by = "key") %>%
				dplyr::select(-c("key")) %>%
				na.omit(.)

			# (3) Plot results
			# Initialization
			dil_factor_min <- 1 # the lower limit dilution factor for Intensity values
			dil_factor_max <- 500000 # the upper limit dilution factor for Intensity values
			dil_factor <- dil_factor_min # the intial value of dilution factor
			fdr <- 1 # the initial FDR
			n_row <- nrow(dfw_ma_sub)

			N <- dfw_ma_sub$Intensity/dil_factor # Scaled intensity for use as an analogous number of coin tosses
			dfw_ma_sub <- within(dfw_ma_sub, {logit_log2Ratio.se <- sqrt(logit_log2Ratio*(1-logit_log2Ratio)/N)}) # standard error at each N
			p_fem <- weighted.mean(dfw_ma_sub$logit_log2Ratio, 1/dfw_ma_sub$logit_log2Ratio.se^2, na.rm=TRUE) # FEM: a fixed, overall probability from meta-data
			p_fem <- 0.5

			# Determine the value of dil_factor to achieve 5% FDR
			tolerance <- 1E-3 # the threshold for FDR
			count <- 0
			while(abs(fdr - .05) > tolerance & count < 100) {
				count <- count + 1
				# Compare experimental logit_log2Ratio and predicted logit_log2Ratio
				dil_factor <- floor((dil_factor_min + dil_factor_max)/2) # new dilution factor
				pred_ll95 <- (p_fem - 1.96 * sqrt((p_fem*(1-p_fem)) / (dfw_ma_sub$Intensity/dil_factor)))
				pred_ul95 <- (p_fem + 1.96 * sqrt((p_fem*(1-p_fem)) / (dfw_ma_sub$Intensity/dil_factor)))
				dfw_ma_sub$is_outlier <- (dfw_ma_sub$logit_log2Ratio < pred_ll95 | dfw_ma_sub$logit_log2Ratio > pred_ul95) # experimental logit_ratios that are outside the range of 95% predicted CI
				N_outliers <- sum(dfw_ma_sub$is_outlier)

				# Update fdr, dil_factor_min and/or dil_factor_max
				fdr <- N_outliers/n_row
				if (fdr < .05) { # dil_factor too large
					dil_factor_max <- dil_factor
				} else { # dil_factor too small
					dil_factor_min <- dil_factor
				}
			}
			rm(count)

			N <- dfw_ma_sub$Intensity/dil_factor # Update the scaled intensity with the dilution factor at 5% FDR
			dfw_ma_sub <- within(dfw_ma_sub, {logit_log2Ratio.se <- sqrt(logit_log2Ratio*(1-logit_log2Ratio)/N)}) # update the standard error at each N
			p_fem <- weighted.mean(dfw_ma_sub$logit_log2Ratio, 1/dfw_ma_sub$logit_log2Ratio.se^2, na.rm=TRUE) # update the FEM probability
			p_fem <- 0.5

			## lower and upper limits for 95% and 99.9% CI, based on FEM estimator
			# number_seq <- seq(0.0001, max(N), 1)
			number_seq <- c(seq(2,10, 10/100) %o% 10^(3:6)) # the series of x values for funnel plot
			number_ll95 <- p_fem - 1.96 * sqrt((p_fem*(1-p_fem)) / (number_seq/dil_factor))
			number_ul95 <- p_fem + 1.96 * sqrt((p_fem*(1-p_fem)) / (number_seq/dil_factor))
			number_ll999 <- p_fem - 3.29 * sqrt((p_fem*(1-p_fem)) / (number_seq/dil_factor))
			number_ul999 <- p_fem + 3.29 * sqrt((p_fem*(1-p_fem)) / (number_seq/dil_factor))
			dfCI <- data.frame(number_ll95, number_ul95, number_ll999, number_ul999, number_seq, p_fem)
			dfCI <- dfCI[dfCI$number_ll95 >= 0 & dfCI$number_ul95 <= 1,]

			p_ma <- ggplot(data = dfw_ma_sub, aes(x = Intensity, y = logit_log2Ratio)) +
				geom_point(color = "darkgray", alpha = .9, size = .01) +
				# geom_line(aes(x = Mean_Int, y = upr_exp_logit_r), color="#252525") +
				# geom_line(aes(x = Mean_Int, y = lwr_exp_logit_r), color="#252525") +
				# stat_density2d(aes(fill = ..level.., alpha = ..level..), geom = "polygon") +
				## stat_density2d(aes(fill = ..density.., alpha = ..density..), geom = "raster", contour=FALSE) +
				scale_fill_distiller(limits = c(0, 1), palette = "Blues", direction = -1, na.value = brewer.pal(n = 9, name = "Blues")[1]) +
				# scale_fill_continuous(terrain.colors(10)) +
				geom_hline(yintercept = p_fem, linetype = "longdash", color = "darkgray", size = .5) +
				# stat_smooth(data=df.ma.trim.sub, se=FALSE, linetype="solid", color="#fb6a4a", size=1) +
				# stat_smooth(data=dfw_ma_sub, se=FALSE, linetype="solid", color="#fb6a4a", size=1) +
				# annotate("text", x = xmin+0.2, y = ymax-.4, hjust = 0, label = Sample, size=6, color="red") +
				labs(title = "", x = expression("Intensity"), y = logit_x_label) +
				# geom_line(aes(x = Intensity, y = lwr_r_logit), color="red", linetype=2) +
				# geom_line(aes(x = Intensity, y = upr_r_logit), color="red", linetype=2) +
				# geom_line(data=dfw_ma_sub, aes(x = Intensity, y = fit_upr_loess), color="green", linetype=2) +
				# geom_line(data=dfw_ma_sub, aes(x = Intensity, y = fit_lwr_loess), color="green", linetype=2) +
				geom_line(aes(x = number_seq, y = number_ll95), color = "#fb6a4a", data = dfCI) +
				geom_line(aes(x = number_seq, y = number_ul95), color = "#fb6a4a", data = dfCI) +
				geom_line(aes(x = number_seq, y = number_ll999), color = "#fb6a4a", linetype = "dashed", data = dfCI) +
				geom_line(aes(x = number_seq, y = number_ul999), color = "#fb6a4a", linetype = "dashed", data = dfCI) +
				# scale_x_continuous(limits=c(xmin, xmax), breaks=c(5 %o% 10^(3:6)), labels = scales::trans_format("log10", scales::math_format(10^.x)), expand=c(0,0)) +
				# scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
				scale_x_continuous(breaks = c(1e+04, 1e+05, 10^5.5, 1e+06), labels = scales::trans_format("log10", scales::math_format(10^.x)), expand = c(0,0)) +
				coord_cartesian(xlim = c(xmin, xmax)) +
				scale_y_continuous(limits = c(0, 1), expand=c(0,0)) +
				facet_wrap(~ Type, ncol = 3, labeller = label_value) +
				my_ma_theme_logit

			Height <- ifelse(lcdn > 2, 4*(lcdn-1)/2, 3)
			ggsave(file.path(filepath, paste0(fn_prx, "~", gsub(":", "", group), "_logit.png")), p_ma, width = 3*lcdn*1.5, height = Height*1.5, units = "in")

		} )

	}

}




#' MA Plots
#'
#' \code{proteoMA} produces MA plots.
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
#' \dontrun{
#' }
#' @import dplyr rlang ggplot2
#' @importFrom magrittr %>%
#' @export
proteoMA <- function (id = gene, col_select = NULL, col_group = NULL, scale_log2r = TRUE,
                      df = NULL, filepath = NULL, filename = NULL, ...) {

  scale_log2r <- match_logi_gv("scale_log2r")

  id <- rlang::enexpr(id)
	if(length(id) != 1) id <- rlang::expr(gene)
	stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))

	col_select <- rlang::enexpr(col_select)
	col_group <- rlang::enexpr(col_group)
	df <- rlang::enexpr(df)
	filepath <- rlang::enexpr(filepath)
	filename <- rlang::enexpr(filename)
	
	reload_expts()

	info_anal(id = !!id, col_select = !!col_select, col_group = !!col_group,
	          scale_log2r = scale_log2r, impute_na = FALSE, df = !!df,
	          filepath = !!filepath, filename = !!filename, anal_type = "MA")(...)
}


#'MA plots of Peptide data
#'@seealso \code{\link{proteoMA}} for parameters
#'@export
pepMA <- function (...) {
	proteoMA(id = pep_seq, ...)
}

#'MA plots of Protein data
#'@seealso \code{\link{proteoMA}} for parameters
#'@export
prnMA <- function (...) {
	proteoMA(id = gene, ...)
}
