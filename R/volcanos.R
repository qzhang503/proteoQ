#' Volcano plots
#' 
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @import limma stringr purrr dplyr   
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotVolcano <- function(df = NULL, df2 = NULL, id = "gene", 
                        adjP = FALSE, topn_labels = 20, 
                        anal_type = "Volcano", 
                        gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), topn_gsets = Inf, 
                        show_sig = "none", fml_nms = NULL, gset_nms = "go_sets", 
                        scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                        filepath = NULL, filename = NULL, theme = NULL, ...) {
  
  stopifnot(vapply(c(scale_log2r, complete_cases, impute_na, adjP), 
                   rlang::is_logical, logical(1)))
  
  stopifnot(vapply(c(gspval_cutoff, gslogFC_cutoff, topn_gsets), is.numeric, 
                   logical(1)))
  
  stopifnot(is.numeric(topn_labels), topn_labels >= 0)
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- df %>%
    `rownames<-`(.[, id]) %>%
    dplyr::select(-grep("I[0-9]{3}|log2_R[0-9]{3}|^.*\\.FC \\(.*\\)$", names(.)))
  
  dots <- rlang::enexprs(...)
  dots <- concat_fml_dots(
    fmls = dots %>% .[grepl("^\\s*~", .)], 
    fml_nms = fml_nms, 
    dots = dots %>% .[!grepl("^\\s*~", .)], 
  )
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  arrange_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange_", names(.))]
  
  select_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^select_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  rm(filter_dots, arrange_dots, select_dots)
  
  if (!adjP) {
    df <- df %>%
      dplyr::select(-contains("adjP"))
  } else {
    df <- df %>%
      dplyr::select(-contains("pVal")) %>%
      `names<-`(gsub("adjP", "pVal", names(.)))
  }
  
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots %>% .[! names(.) %in% names(fmls)]
  
  fml_nms <- names(df) %>% 
    .[grepl("pVal", .)] %>% 
    gsub("(.*)\\.pVal.*", "\\1", .) %>% 
    unique()  %>% 
    .[. %in% names(fmls)]
  
  fmls <- fmls %>% .[names(.) %in% fml_nms]
  fml_nms <- fml_nms %>% .[map_dbl(., ~ which(.x == names(fmls)))]
  
  if (rlang::is_empty(fml_nms)) {
    stop("No formula (names) matched to those in `pepSig(...)` or `prnSig(...)`", 
         call. = FALSE)
  }
  
  purrr::walk(file.path(filepath, fml_nms), 
              ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))
  
  col_ind <- purrr::map(fml_nms, ~ grepl(.x, names(df))) %>%
    `names<-`(paste0("nm_", seq_along(.))) %>% 
    dplyr::bind_cols() %>%
    rowSums() %>%
    magrittr::is_greater_than(0)
  
  species <- df$species %>% unique() %>% .[!is.na(.)] %>% as.character()
  if (!(purrr::is_empty(species) || is.null(gset_nms))) {
    load_dbs(gset_nms = gset_nms, species = species)
  }
  
  proteoq_volcano_theme <- theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.5, size = 24),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 0, vjust = 0.5, size = 24),
      axis.title.x = element_text(colour = "black", size = 24),
      axis.title.y = element_text(colour="black", size = 24),
      plot.title = element_text(face = "bold", colour = "black", 
                                size = 14, hjust = .5, vjust = .5),
      
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
  
  if (is.null(theme)) theme <- proteoq_volcano_theme

  if (length(fml_nms) > 0) purrr::pwalk(list(fml_nms, 
                                             gspval_cutoff, 
                                             gslogFC_cutoff, 
                                             topn_gsets), 
                                        byfml_volcano, 
                                        df = df, 
                                        df2 = df2, 
                                        col_ind = col_ind, 
                                        id = !!id, 
                                        filepath = filepath, 
                                        filename = filename, 
                                        adjP = adjP, 
                                        topn_labels = topn_labels, 
                                        anal_type = anal_type, 
                                        show_sig = show_sig, 
                                        gset_nms = gset_nms, 
                                        scale_log2r = scale_log2r, 
                                        complete_cases = complete_cases, 
                                        impute_na = impute_na, 
                                        theme = theme, 
                                        !!!dots)
}


#' Formula specific volcano plots
#' 
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import purrr dplyr 
#' @importFrom magrittr %>% %T>% %$% %<>% 
byfml_volcano <- function (fml_nm, gspval_cutoff, gslogFC_cutoff, topn_gsets, df, df2, 
                           col_ind, id, 
                           filepath, filename, adjP, topn_labels, 
                           anal_type, 
                           show_sig, gset_nms, scale_log2r, complete_cases, impute_na, 
                           theme = NULL, ...) {

  pval_complete_cases <- function (df) {
    rows <- df %>% 
      dplyr::select(grep("^pVal\\s+", names(.))) %>% 
      complete.cases()
    
    if (sum(rows) == 0) {
      stop("None of the cases are complete.", call. = FALSE)
    }
    
    df <- df[rows, ]
  }  
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  df <- df %>%
    dplyr::select(grep(paste0("^", fml_nm, "\\."), names(.))) %>%
    `colnames<-`(gsub(paste0("^", fml_nm, "\\."), "", names(.))) %>%
    dplyr::bind_cols(df[, !col_ind, drop = FALSE], .) 
  
  if (complete_cases) df <- df %>% pval_complete_cases()

  byfile_plotVolcano(df = df, df2 = df2, id = !!id, fml_nm = fml_nm, 
                     filepath = filepath, filename = filename, 
                     adjP = adjP, topn_labels = topn_labels, 
                     anal_type = anal_type, gset_nms = gset_nms, 
                     scale_log2r = scale_log2r, 
                     impute_na = impute_na)(gspval_cutoff = gspval_cutoff, 
                                            gslogFC_cutoff = gslogFC_cutoff, 
                                            topn_gsets = topn_gsets, 
                                            show_sig = show_sig, 
                                            theme = theme, ...)
}


#' Plot volcanos
#' 
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import limma stringr purrr dplyr   
#' @importFrom magrittr %>% %T>% %$% %<>% 
byfile_plotVolcano <- function(df = NULL, df2 = NULL, id = "gene", fml_nm = NULL, 
                               filepath = NULL, filename = NULL, 
                               adjP = FALSE, topn_labels = 20, 
                               anal_type = "Volcano", gset_nms = "go_sets", 
                               scale_log2r, impute_na, theme = NULL, ...) {

  id <- rlang::as_string(rlang::enexpr(id))

	contrast_groups <- names(df[grep("^log2Ratio\\s+\\(", names(df))]) %>%
	  gsub("^log2Ratio\\s+\\(|\\)$", "", .)

	if (anal_type %in% c("Volcano")) {
		function(gspval_cutoff = 1, gslogFC_cutoff = 0, topn_gsets = 0, 
		         show_sig = "none", theme = theme, ...) {
			rm(gspval_cutoff, gslogFC_cutoff, show_sig)
      
			fullVolcano(df = df, 
			            id = !!id, 
			            contrast_groups = contrast_groups,
			            theme = theme, 
			            fml_nm = fml_nm, 
			            filepath = filepath, 
			            filename = filename, 
			            adjP = adjP, 
			            topn_labels = topn_labels, 
			            ...)
		}
	} else if (anal_type == "mapGSPA")
	  function(gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), topn_gsets = Inf, 
	           show_sig = "none", theme = theme, ...) {
	    stopifnot(!is.null(fml_nm))
	    
	    filepath_fml <- file.path(filepath, fml_nm)
	    
	    in_names <- list.files(path = filepath_fml, pattern = "_GSPA_[NZ]{1}.*\\.txt$")
	    if (purrr::is_empty(in_names)) {
	      warning("No inputs under ", filepath_fml, call. = FALSE)
	      return(NULL)
	    }
	    
	    if (is.null(df2)) {
  	    in_names <- in_names %>% 
  	      .[!grepl("_essmap|_essmeta|_resgreedy", .)] %>% 
  	      {if (impute_na) .[grepl("_impNA", .)] else .[!grepl("_impNA", .)]} %>% 
  	      {if (scale_log2r) .[grepl("_GSPA_Z", .)] else .[grepl("_GSPA_N", .)]}

  	    if (rlang::is_empty(in_names)) {
  	      stop("No inputs correspond to impute_na = ", impute_na, 
  	           ", scale_log2r = ", scale_log2r, 
  	           " at fml_nms = ", fml_nm, 
  	           call. = FALSE)
  	    }
  	    
  	    df2 <- in_names
	    } else {
	      local({
	        if (grepl("_essmap|_essmeta|_resgreedy", df2)) {
	          stop("Do not use `_essmap`, `_essmeta` or `_resgreedy` for `df2`", 
	               call. = FALSE)
	        }
	        
	        non_exists <- df2 %>% .[! . %in% in_names]
	        if (!purrr::is_empty(non_exists)) {
	          stop("Missing file(s): ", purrr::reduce(non_exists, paste, sep = ", "), 
	               call. = FALSE)
	        }
	      })
	      
	      if (purrr::is_empty(df2)) {
	        stop("File(s) not found under ", filepath_fml, 
	             call. = FALSE)
	      }
	    }

	    # plot data ---------------------------
	    purrr::walk(df2, gsVolcano, 
	                df = df, 
	                contrast_groups = contrast_groups, 
	                gsea_key = "term", 
	                gsets = gsets, 
	                theme = theme, 
	                fml_nm = fml_nm, 
	                filepath = filepath, 
	                filename = filename, 
	                adjP = adjP, 
	                topn_labels = topn_labels, 
	                show_sig = show_sig, 
	                gspval_cutoff = gspval_cutoff, 
	                gslogFC_cutoff = gslogFC_cutoff, 
	                topn_gsets = topn_gsets, 
	                ...)
	  }
}


#' Volcano plots for all proteins or peptides in a data set
#' 
#' @param contrast_groups The contrast groups defined under a formula at \code{fml_nm}.
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @import dplyr purrr ggplot2
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom limma vennDiagram
fullVolcano <- function(df = NULL, id = "gene", contrast_groups = NULL, theme = NULL,
                        fml_nm = NULL, filepath = NULL, filename = NULL, adjP = FALSE, 
                        topn_labels = 20, ...) {

  id <- rlang::as_string(rlang::enexpr(id))

	dots <- rlang::enexprs(...)
	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)
	
	title <- dots[["title"]]
  
	if (is.null(dots$x_label)) {
	  x_label <- expression("Ratio ("*log[2]*")")
	} else {
	  x_label <- dots$x_label
	}

	if (is.null(dots$y_label)) {
	  y_label <- ifelse(adjP, 
	                    expression("adjP ("*-log[10]*")"), 
	                    expression("pVal ("*-log[10]*")"))
	} else {
	  y_label <- dots$y_label
	}
	
	stopifnot(length(contrast_groups) > 0)

	dfw <- do.call(rbind, purrr::map(contrast_groups, ~ {
			df[, grepl(paste0(" (", .x, ")"), names(df), fixed = TRUE)] %>%
				`colnames<-`(gsub("\\s+\\(.*\\)$", "", names(.))) %>%
				dplyr::mutate(Contrast = .x) %>%
				dplyr::bind_cols(df[, !grepl("^pVal\\s+|^adjP\\s+|^log2Ratio\\s+", names(df)), 
				                    drop = FALSE], .)
		} )) %>% 
	  # dplyr::filter(rowSums(!is.na(.[, grepl("\\.log2Ratio \\(", names(.))])) > 0) %>% 
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
		dplyr::top_n(n = -topn_labels, wt = pVal) %>%
		data.frame (check.names = FALSE)

	# data table for labels
	dt <- purrr::map(contrast_groups, ~ {
	  to_csv_(dfw_sub_top20 %>%
            dplyr::filter(Contrast == .x) %>%
            dplyr::select(c("Index", id))) %>%
            {if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .} # zero-entry exception
	  })  %>%
	  do.call(rbind, .) %>%
	  data.frame(Contrast = contrast_groups, id = ., stringsAsFactors = FALSE) %>%
	  dplyr::rename(!!rlang::sym(id) := id) %>%
	  dplyr::mutate(Contrast = factor(Contrast,  levels = contrast_groups))

	fn_prefix <- gsub("\\.[^.]*$", "", filename)

	myPalette <- c("#377EB8", "#E41A1C")
	dfw_greater <- dfw_sub[dfw_sub$valence == "pos", ]
	dfw_less <- dfw_sub[dfw_sub$valence == "neg", ]

	if (is.null(dots$nrow)) {
		nrow <- ifelse(length(unique(dfw$Contrast)) > 3, 2, 1)
	} else {
		nrow <- dots$nrow
	}
	
	if (is.null(dots$xmax)) {
	  xmax <- ceiling(pmax(abs(min(dfw$log2Ratio, na.rm = TRUE)), 
	                       max(dfw$log2Ratio, na.rm = TRUE)))
	} else {
	  xmax <- eval(dots$xmax)
	  stopifnot(xmax > 0)
	}
	
	if (is.null(dots$xmin)) {
	  xmin <- -xmax
	} else {
	  xmin <- eval(dots$xmin)
	  stopifnot(xmin < 0)
	}

	if (is.null(dots$ymax)) {
	  ymax <- ceiling(max(-log10(dfw$pVal), na.rm = TRUE)) * 1.1
	} else {
	  ymax <- eval(dots$ymax)
	  stopifnot(ymax > 0)
	}	

	if (is.null(dots$ymin)) {
	  ymin <- 0
	} else {
	  ymin <- eval(dots$ymin)
	  stopifnot(ymin < ymax)
	}	
	
	p <- ggplot() +
		geom_point(data = dfw, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, colour = "#252525", shape = 20, alpha = .5) +
		geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, color = myPalette[2], shape = 20, alpha = .8) +
		geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
		           size = 3, color = myPalette[1], shape = 20, alpha = .8) +
		geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
		geom_vline(xintercept = log2(xco), linetype = "longdash", size = .5) +
		scale_x_continuous(limits = c(xmin, xmax)) +
		scale_y_continuous(limits = c(ymin, ymax)) +
		labs(title = title, x = x_label, y = y_label) +
		theme

	if (nrow(dfw_sub_top20) > 0) {
		p <- p + geom_text(data = dfw_sub_top20, 
		                   mapping = aes(x = log2Ratio, 
		                                 y = -log10(pVal), 
		                                 label = Index, 
		                                 color = Index),
		                   size = 3, 
		                   alpha = .5, 
		                   hjust = 0, 
		                   nudge_x = 0.05, 
		                   vjust = 0, 
		                   nudge_y = 0.05, 
		                   na.rm = TRUE)
		p <- p + facet_wrap(~ Contrast, nrow = nrow, labeller = label_value)
		p <- p + geom_table(data = dt, aes(table = !!rlang::sym(id)), 
		                    x = -xmax*.85, y = ymax/2)
	} else {
		p <- p + facet_wrap(~ Contrast, nrow = nrow, labeller = label_value)
	}

	if(is.null(dots$width)) {
		width <- ifelse(nrow > 1, 6*length(unique(dfw$Contrast))/nrow + 1, 
		                6*length(unique(dfw$Contrast)) / nrow)
	} else {
		width <- eval(dots$width)
	}

	if(is.null(dots$height)) {
		height <- 6*nrow
	} else {
		height <- eval(dots$height)
	}

	ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
	
	rlang::eval_tidy(rlang::quo(ggsave(filename = file.path(filepath, fml_nm, filename),
	                                   plot = p, 
	                                   width = width, 
	                                   height = height, 
	                                   !!!ggsave_dots)))

	# Venn
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

		fn_prefix <- paste0(direction, "_venn")

		png(file.path(filepath, fml_nm, paste0(fn_prefix, ".png")), 
		    width = Width, height = Height, units = "in", res = 300)
			limma::vennDiagram(limma::vennCounts(Counts), circle.col = myPalette, cex = .5)
		dev.off()
		write.table(Counts, file.path(filepath, fml_nm, paste0(fn_prefix, ".txt")), sep = "\t", 
		            col.names = TRUE, row.names = TRUE, quote = FALSE)	
	}

	summ_venn(dfw_greater, id, contrast_groups) %>% 
	  plot_venn(filepath, paste0(fn_prefix, "_greater"))
	
	summ_venn(dfw_less, id, contrast_groups) %>% 
	  plot_venn(filepath, paste0(fn_prefix, "_less"))

	saveRDS(
	  list(
	    data = dfw, 
	    greater = dfw_greater, 
	    less = dfw_less, 
	    topns = dfw_sub_top20, 
	    topn_labels = dt, 
	    palette = myPalette, 
	    xco = xco, 
	    yco = yco,
	    xmin = xmin, 
	    xmax = xmax,
	    ymin = ymin,
	    ymax = ymax, 
	    title = title, 
	    x_label = x_label,
	    y_label = y_label,
	    theme = theme
	  ), 
	  file.path(filepath, fml_nm, paste0(fn_prefix, ".rds"))
	)
	
}


#' Volcano plots of protein \code{log2FC} under given gene sets
#'
#' @param gsea_key Character string; the column key indicating the terms of gene sets.
#' @param gsets The gene sets.
#' @inheritParams info_anal
#' @inheritParams prnVol
#' @inheritParams gspaMap
#' @inheritParams fml_gspa
#' @inheritParams fullVolcano
#' @import dplyr ggplot2
#' @importFrom magrittr %>% %T>% %$% %<>% 
gsVolcano <- function(df2 = NULL, df = NULL, contrast_groups = NULL, 
                      gsea_key = "term", gsets = NULL, 
                      theme = NULL, 
                      fml_nm = NULL, 
                      filepath = NULL, filename = NULL, adjP = FALSE, 
                      topn_labels = 20, show_sig = "none", 
                      gspval_cutoff = 1E-6, gslogFC_cutoff = log2(1.2), 
                      topn_gsets = Inf, ...) {
  
  dat_dir <- get_gl_dat_dir()
  par_filepath <- file.path(dat_dir, "Calls", 
                            gsub("\\.txt$", "@", df2) %>% paste0(fml_nm, ".rda"))
  
  if (file.exists(par_filepath)) load(par_filepath) else call_pars <- NULL
  
  df <- local({
    par_filter_dots <- call_pars %>% 
      .[purrr::map_lgl(., is.language)] %>% 
      .[grepl("^filter_", names(.))]
    
    if (!purrr::is_empty(par_filter_dots)) {
      message(
        "\nApply the following vararg(s) from matching `prnGSPA` to: ", fml_nm, ".", 
        paste(par_filter_dots, "\n")
      )
      df %>% filters_in_call(!!!par_filter_dots)
    } else {
      message("\nNo `filter_` varargs with formula `", fml_nm, "`.")
      invisible(df)
    }
  })
  
  local({
    use_adjP <- call_pars$use_adjP
    
    if (use_adjP != adjP) {
      warning("\n", 
              "Current analysis: `adjP = ", adjP, "`.\n", 
              "Data `", df2, "@", fml_nm, 
              "`: based on `prnGSPA(use_adjP = ", use_adjP, ")`\n", 
              "May consider `gspaMap(adjP = ", use_adjP, ")` for ", 
              df2, "@", fml_nm, ".\n", 
              "  See also `?gspaMap` for analyses against a specific ", 
              "file and and formula.\n\n", 
              call. = FALSE)
    }
  })
  
  pval_cutoff <- local({
    par_pval_cutoff <- call_pars$pval_cutoff
    if (!purrr::is_empty(par_pval_cutoff)) par_pval_cutoff else 0.05
  })
  
  logFC_cutoff <- local({
    par_logFC_cutoff <- call_pars$logFC_cutoff
    if (!purrr::is_empty(par_logFC_cutoff)) par_logFC_cutoff else log2(1.2)
  })
  
  custom_prefix <- gsub("(.*_{0,1})Protein_GSPA.*", "\\1", df2)
  dir.create(path = file.path(filepath, fml_nm, custom_prefix), 
             recursive = TRUE, showWarnings = FALSE)

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) 
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  filename <- paste0(custom_prefix, fn_prefix, ".", fn_suffix)
  
  dots <- rlang::enexprs(...)
  
	xco <- ifelse(is.null(dots$xco), 1.2, dots$xco)
	yco <- ifelse(is.null(dots$yco), .05, dots$yco)

	if (is.null(dots$x_label)) {
	  x_label <- expression("Ratio ("*log[2]*")")
	} else {
	  x_label <- dots$x_label
	}
	
	if (is.null(dots$y_label)) {
	  y_label <- ifelse(adjP, 
	                    expression("adjP ("*-log[10]*")"),
	                    expression("pVal ("*-log[10]*")"))
	} else {
	  y_label <- dots$y_label
	}
	
	gsea_res <- tryCatch(
	  readr::read_tsv(file.path(filepath, fml_nm, df2), 
	                  col_types = cols(term = col_character(),
	                                   is_essential = col_logical(),
	                                   size = col_double(),
	                                   ess_size = col_double(),
	                                   contrast = col_character(),
	                                   p_val = col_double(),
	                                   q_val = col_double(),
	                                   log2fc = col_double())), 
	  error = function(e) NA
	)
	
	message("Secondary file loaded: ", file.path(filepath, fml_nm, df2))

	filter2_dots <- dots %>% 
	  .[purrr::map_lgl(., is.language)] %>% 
	  .[grepl("^filter2_", names(.))]
	
	arrange2_dots <- dots %>% 
	  .[purrr::map_lgl(., is.language)] %>% 
	  .[grepl("^arrange2_", names(.))]
	
	select2_dots <- dots %>% 
	  .[purrr::map_lgl(., is.language)] %>% 
	  .[grepl("^select2_", names(.))]
	
	dots <- dots %>% 
	  .[! . %in% c(filter2_dots, arrange2_dots, select2_dots)]
	
	gsea_res <- gsea_res %>% 
	  dplyr::arrange(p_val, abs(log2fc)) %>% 
	  filters_in_call(!!!filter2_dots) %>% 
	  arrangers_in_call(!!!arrange2_dots)
	rm(filter2_dots, arrange2_dots, select2_dots)

	if (nrow(gsea_res) == 0) {
	  stop("No GSPA terms available after data filtration.", 
	       call. = FALSE)
	}

	topn_gsets <- pmin(dplyr::n_distinct(gsea_res[[gsea_key]]), topn_gsets)
	terms <- gsea_res %>%
	  dplyr::arrange(p_val) %>% 
	  dplyr::slice(1:(topn_gsets*length(contrast_groups))) %>% 
	  dplyr::filter(p_val <= gspval_cutoff, 
	                abs(log2fc) >= gslogFC_cutoff) %>%
	  dplyr::select(gsea_key) %>%
	  unique() %>%
	  unlist() %>% 
	  .[. %in% names(gsets)]
	
	gsea_res <- gsea_res %>% 
	  dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>% 
	  dplyr::mutate(q_val = format(p_val, scientific = TRUE, digits = 2)) %>% 
	  dplyr::mutate(log2fc = round(log2fc, digits = 2))

	if (length(terms) > 0) {
  	dfw <- do.call(rbind,
  		purrr::map(contrast_groups, ~ {
  			df[, grepl(paste0(" (", .x, ")"), names(df), fixed = TRUE)] %>%
  				`colnames<-`(gsub("\\s+\\(.*\\)$", "", names(.))) %>%
  				dplyr::mutate(Contrast = .x) %>%
  		    dplyr::bind_cols(df[, !grepl("^pVal\\s+|^adjP\\s+|^log2Ratio\\s+", names(df))], .)
  		} )) %>%
  		dplyr::mutate(Contrast = factor(Contrast, levels = contrast_groups),
  			pVal = as.numeric(pVal),
  			valence = ifelse(.$log2Ratio > 0, "pos", "neg")) %>%
  		dplyr::filter(!is.na(pVal))

  	lapply(terms, function(gt) {
  	  # some results may be based on gene sets from an older database, 
  	  # which become missing terms in the current
  	  
  	  gsets_sub <- gsets %>% .[names(.) == gt]
  		if (length(gsets_sub) == 0) return(NULL)

  		fn <- gsub(":", "~", gsub("/", "or", names(gsets_sub)[[1]]), fixed = TRUE)

  		res_sub <- gsea_res[as.character(gsea_res$term) == gt, ] %>% 
  		  data.frame(check.names = FALSE)

  		dfw_sub <- dfw[as.character(dfw$entrez) %in% gsets_sub[[1]], ]

  		if (nrow(dfw_sub) == 0) return(NULL) 

  		if (is.null(dots$xmax)) {
  		  xmax <- ceiling(pmax(abs(min(dfw_sub$log2Ratio, na.rm = TRUE)), 
  		                       max(dfw_sub$log2Ratio, na.rm = TRUE)))
  		} else {
  		  xmax <- eval(dots$xmax)
  		  stopifnot(xmax > 0)
  		}
  		
  		if (is.null(dots$xmin)) {
  		  xmin <- -xmax
  		} else {
  		  xmin <- eval(dots$xmin)
  		  stopifnot(xmin < 0)
  		}
  		
  		if (is.null(dots$ymax)) {
  		  ymax <- ceiling(max(-log10(dfw_sub$pVal))) * 1.1
  		} else {
  		  ymax <- eval(dots$ymax)
  		  stopifnot(ymax > 0)
  		}
  		
  		if (is.null(dots$ymin)) {
  		  ymin <- 0
  		} else {
  		  ymin <- eval(dots$ymin)
  		  stopifnot(ymin < ymax)
  		}	

  		# ensure the same levels between "Levels" and "newLevels"
  		Levels <- levels(dfw_sub$Contrast)

  		dfw_sub <- dfw_sub %>%
  			dplyr::arrange(Contrast, pVal) %>%
  			dplyr::group_by(Contrast) %>%
  			dplyr::mutate(Index = row_number()) %>%
  			data.frame(check.names = FALSE) %>% 
  		  dplyr::mutate(Contrast = as.character(Contrast)) %>% 
  			dplyr::left_join(., res_sub, by = c("Contrast" = "contrast")) %>%
  			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
  			dplyr::arrange(Contrast) %>%
  			# dplyr::mutate(p_val = format(p_val, scientific = TRUE, digits = 2)) %>%
  			# dplyr::mutate(p_val = as.numeric(p_val)) %>%
  			# dplyr::mutate(q_val = format(q_val, scientific = TRUE, digits = 2)) %>%
  			# dplyr::mutate(q_val = as.numeric(q_val)) %>%
  			# dplyr::mutate(sig_level = ifelse(.$q.val > 0.05, "n.s.", ifelse(.$q.val > 0.005, "*", "**"))) %>%
  			# dplyr::mutate(newContrast = paste0(Contrast, " (", sig_level, ")"))
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
  			dplyr::top_n(n = -topn_labels, wt = pVal)

  		dt <- purrr::map(newLevels, ~ {
  				dfw_sub_top20 %>%
  					dplyr::filter(newContrast == .x) %>%
  					data.frame(check.names = FALSE) %>%
  					dplyr::select(c("Index", "gene")) %>%
  					to_csv_() %>%
  					{if(!grepl("\n", .)) . <- paste0(.,"\n1,\"NA\"") else .}
  			}) %>%
  			do.call(rbind, .) %>%
  			data.frame(newContrast = newLevels, 
  			           Contrast = Levels, 
  			           Gene = ., 
  			           stringsAsFactors = FALSE) %>%
  			dplyr::mutate(Contrast = factor(Contrast, levels = Levels)) %>%
  			dplyr::mutate(newContrast = factor(newContrast, levels = newLevels))

  		dfw_greater <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio > log2(xco))
  		dfw_less <- dfw_sub %>% dplyr::filter(pVal < yco & log2Ratio < -log2(xco))

  		dt_pos <- ifelse(nrow(dfw_greater) > nrow(dfw_less), -xmax*.85, xmax*.6)
  		myPalette <- c("#377EB8", "#E41A1C")
  		
  		if(is.null(dots$nrow)) {
  		  nrow <- ifelse(length(unique(dfw_sub$Contrast)) > 3, 2, 1)
  		} else {
  		  nrow <- dots$nrow
  		}

  		p <- ggplot() +
  			geom_point(data = dfw_sub, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
  			           size = 3, colour = "gray", shape = 20, alpha = .5) +
  			geom_point(data = dfw_greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
  			           size = 3, color = myPalette[2], shape = 20, alpha = .8) +
  			geom_point(data = dfw_less, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
  			           size = 3, color = myPalette[1], shape = 20, alpha = .8) +
  			geom_hline(yintercept = -log10(yco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = -log2(xco), linetype = "longdash", size = .5) +
  			geom_vline(xintercept = log2(xco), linetype = "longdash", size =.5) +
  		  geom_hline(yintercept = -log10(pval_cutoff), 
  		             linetype = "longdash", size = .5, color = "#fc9272") +
  		  geom_vline(xintercept = -logFC_cutoff, 
  		             linetype = "longdash", size = .5, color = "#fc9272") +
  		  geom_vline(xintercept = logFC_cutoff, 
  		             linetype = "longdash", size =.5, color = "#fc9272") +
  			labs(title = names(gsets_sub), x = x_label, y = y_label) +
  			scale_x_continuous(limits = c(xmin, xmax)) +
  			scale_y_continuous(limits = c(ymin, ymax)) +
  			theme
  		p <- p + facet_wrap(~ newContrast, nrow = nrow, labeller = label_value)

  		if (nrow(dfw_sub_top20) > 0) {
  		  p <- p + geom_text(data = dfw_sub_top20, 
  		                     mapping = aes(x = log2Ratio, y = -log10(pVal),
  		                                   label = Index, color = Index), 
  		                     size = 2, hjust = 0, nudge_x = 0.05, vjust = 0, nudge_y = 0.05) +
  		    geom_table(data = dt, aes(table = Gene), x = dt_pos, y = ymax/2)
  		}

  		if (nchar(fn) > 50) fn <- paste0(str_sub(fn, 1, 50), "...") 
  		
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
  		
  		ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
  		
  		rlang::eval_tidy(
  		  rlang::quo(
  		    ggsave(filename = file.path(filepath, fml_nm, custom_prefix, 
  		                                paste0(fn, ".", fn_suffix)),
  		           plot = p, 
  		           width = width, 
  		           height = height, 
  		           !!!ggsave_dots)
  		  )
  		)
  		
  		write.table(dfw_sub, 
  		            file = file.path(filepath, fml_nm, custom_prefix, paste0(fn, ".txt")), 
  		            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)	
  		            
  	})
	}
}



#'Volcano plots
#'
#'\code{pepVol} visualizes the volcano plots of peptide data.
#'
#'@rdname prnVol
#'
#'@import purrr
#'@export
pepVol <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjP = FALSE, topn_labels = 20, 
                    df = NULL, filepath = NULL, filename = NULL, 
                    fml_nms = NULL, theme = NULL, ...) {
  
  check_dots(c("id", "anal_type", "df2"), ...)
  check_depreciated_args(list(c("show_labels", "topn_labels")), ...)
  
  id <- match_call_arg(normPSM, group_psm_by)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod"), 
            length(id) == 1)
  
  scale_log2r <- match_pepSig_scale_log2r(scale_log2r = scale_log2r, 
                                          impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)	
  
  reload_expts()

  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            filepath = !!filepath, 
            filename = !!filename, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            anal_type = "Volcano")(fml_nms = fml_nms, 
                                   adjP = adjP, 
                                   topn_labels = topn_labels, 
                                   theme = theme, 
                                   ...)
}


#'Volcano plots
#'
#'\code{prnVol} visualizes the volcano plots of protein data.
#'
#'@inheritParams prnGSPA
#'@inheritParams prnHist
#'@param adjP Logical; if TRUE, use Benjamini-Hochberg pVals in volcano plots.
#'  The default is FALSE.
#'@param topn_labels A non-negative integer; the top-n species for labeling in a
#'  plot. At \code{topn_labels = 0}, no labels of proteins/peptides will be
#'  shown. The default is to label the top-20 species with the lowest p-values.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  Additional parameters for plotting: \cr \code{xco}, the cut-off lines of
#'  fold changes at position \code{x}; the default is at \eqn{-1.2} and
#'  \eqn{+1.2}. \cr \code{yco}, the cut-off line of \code{pVal} at position
#'  \code{y}; the default is \eqn{0.05}. \cr \code{width}, the width of plot;
#'  \cr \code{height}, the height of plot. \cr \code{nrow}, the number of rows
#'  in a plot. \cr \code{xmin}, the minimum \code{x}. \cr \code{xmax}, the
#'  maximum \code{x}. \cr \code{ymin}, the minimum \code{y}. \cr \code{ymax},
#'  the maximum \code{y}. \cr \code{x_label}, the label on \code{x}. \cr
#'  \code{y_label}, the label on \code{y}.
#'
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'
#'@example inst/extdata/examples/prnVol_.R
#'
#'@seealso 
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
#'@export
prnVol <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                    adjP = FALSE, topn_labels = 20, 
                    df = NULL, filepath = NULL, filename = NULL, 
                    fml_nms = NULL, theme = NULL, ...) {
  
  check_dots(c("id", "anal_type", "df2"), ...)
  check_depreciated_args(list(c("show_labels", "topn_labels")), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1)
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, 
                                          impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)	
  
  reload_expts()

  info_anal(df = !!df, 
            df2 = NULL, 
            id = !!id, 
            filepath = !!filepath, 
            filename = !!filename, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            anal_type = "Volcano")(fml_nms = fml_nms, 
                                   adjP = adjP, 
                                   topn_labels = topn_labels, 
                                   theme = theme, 
                                   ...)
}


#'Volcano plots of protein \code{log2FC} under gene sets
#'
#'\code{gspaMap} visualizes the volcano plots of protein subgroups under the
#'same gene sets.
#'
#'@inheritParams prnGSPA
#'@inheritParams prnVol
#'@inheritParams prnHist
#'@inheritParams plot_prnTrend
#'@inheritParams anal_prnTrend
#'@param filename Use system default for each gene set.
#'@param show_sig Character string indicating the type of significance values to
#'  be shown with \code{\link{gspaMap}}. The default is \code{"none"}.
#'  Additional choices are from \code{c("pVal", "qVal")} where \code{pVal} or
#'  \code{qVal} will be shown, respectively, in the facet grid of the plots.
#'@param gspval_cutoff Numeric value or vector for uses with
#'  \code{\link{gspaMap}}. \code{Gene sets} with enrichment \code{pVals} less
#'  significant than the threshold will be excluded from volcano plot
#'  visualization. The default significance is 0.05 for all formulas matched to
#'  or specified in argument \code{fml_nms}. Formula-specific threshold is
#'  allowed by supplying a vector of cut-off values.
#'@param gslogFC_cutoff Numeric value or vector for uses with
#'  \code{\link{gspaMap}}. \code{Gene sets} with absolute enrichment
#'  \code{log2FC} less than the threshold will be excluded from volcano plot
#'  visualization. The default magnitude is \code{log2(1.2) } for all formulas
#'  matched to or specified in argument \code{fml_nms}. Formula-specific
#'  threshold is allowed by supplying a vector of absolute values in
#'  \code{log2FC}.
#'@param topn_gsets Numeric value or vector; top entries in gene sets ordered by
#'  increasing \code{pVal} for visualization. The default is to use all
#'  available entries.
#'
#'  Note that it is users' responsibility to ensure that the custom gene sets
#'  contain terms that can be found from the one or multiple preceding analyses
#'  of \code{\link{prnGSPA}}. For simplicity, it is generally applicable to
#'  include \emph{all} of the data bases that have been applied to
#'  \code{\link{prnGSPA}} and in that way no terms will be missed for
#'  visualization. See also \code{\link{prnGSPA}} for examples of custom data
#'  bases.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements and the
#'  association between \code{filter_} and \code{df}. \cr \cr \code{filter2_}:
#'  Variable argument statements for the row filtration against data in
#'  secondary file(s) linked to \code{df2}. See also \code{\link{prnGSPAHM}} for
#'  the format of \code{filter2_}, \code{normPSM} for the association between
#'  \code{filter_} and \code{df}. \cr \cr Additional parameters for plotting:
#'  \cr \code{xco}, the cut-off lines of fold changes at position \code{x}; the
#'  default is at \eqn{-1.2} and \eqn{+1.2}. \cr \code{yco}, the cut-off line of
#'  \code{pVal} at position \code{y}; the default is \eqn{0.05}. \cr
#'  \code{width}, the width of plot; \cr \code{height}, the height of plot. \cr
#'  \code{nrow}, the number of rows in a plot.
#'
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'
#'@example inst/extdata/examples/prnVol_.R
#'
#'@seealso 
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
#'@export
gspaMap <- function (gset_nms = c("go_sets", "c2_msig", "kinsub"), 
                     scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                     df = NULL, df2 = NULL, 
                     filepath = NULL, filename = NULL, fml_nms = NULL, 
                     adjP = FALSE, topn_labels = 20, 
                     show_sig = "none", 
                     gspval_cutoff = 5E-2, gslogFC_cutoff = log2(1.2), topn_gsets = Inf, 
                     theme = NULL, ...) {
  check_dots(c("id", "anal_type"), ...)
  check_depreciated_args(list(c("pval_cutoff", "gspval_cutoff"), 
                              c("logFC_cutoff", "gslogFC_cutoff"), 
                              c("show_labels", "topn_labels")), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1)
  
  scale_log2r <- match_prnSig_scale_log2r(scale_log2r = scale_log2r, 
                                          impute_na = impute_na)
  
  df <- rlang::enexpr(df)
  df2 <- rlang::enexpr(df2)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  show_sig <- rlang::as_string(rlang::enexpr(show_sig))
  stopifnot(show_sig %in% c("none", "pVal", "qVal"), 
            length(show_sig) == 1)

  check_gset_nms(gset_nms)

  reload_expts()

  info_anal(df = !!df, 
            df2 = !!df2, 
            id = !!id, 
            filepath = !!filepath, 
            filename = !!filename, 
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na, 
            anal_type = "mapGSPA")(fml_nms = fml_nms, 
                                   adjP = adjP, 
                                   topn_labels = topn_labels, 
                                   gspval_cutoff = gspval_cutoff, 
                                   gslogFC_cutoff = gslogFC_cutoff, 
                                   topn_gsets = topn_gsets, 
                                   show_sig = show_sig,
                                   gset_nms = gset_nms, 
                                   theme = theme, 
                                   ...)
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

  	vp <- grid::viewport(x = mean(x_rng), y = mean(y_rng), 
  	                     width = diff(x_rng), height = diff(y_rng), 
  	               just = c("center", "center"))
  	  
  	grob <- gridExtra::tableGrob(table, rows = NULL, cols = NULL,
  	                  theme = 
  	                    gridExtra::ttheme_minimal(core = list(fg_params=list(cex = .7)),
  	                                         colhead = list(fg_params=list(cex = .7), 
  	                                                        parse=TRUE),
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
  	grid::editGrob(grob, vp = vp, name = paste(grob$name, facet_id()))
	}
)


facet_id <- local({
	i <- 1
	function() {
  	i <<- i + 1
  	i
	}
})

#' Print table in ggplot2 images
#' 
#' @param mapping The same as ggplot2.
#' @param data same as ggplot2.
#' @param stat The same as ggplot2.
#' @param position same as ggplot2.
#' @param na.rm The same as ggplot2.
#' @param show.legend same as ggplot2.
#' @param inherit.aes The same as ggplot2.
#' @param ... same as ggplot2.
#' @export
geom_table <- function(mapping = NULL, data = NULL, stat = "identity",
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, ...) {
  layer(geom = GeomTable, 
        mapping = mapping, 
        data = data, 
        stat = stat, 
        position = position,
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(na.rm = na.rm, ...)
  )
}


#' Convert a data frame column to a csv vector
#' 
#' @param x A data frame column.
to_csv_ <- function(x) {
	paste(capture.output(write.csv(x, stdout(), row.names = F)), collapse = "\n")
}

