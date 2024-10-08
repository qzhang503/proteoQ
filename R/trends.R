#' Check the dot-dot-dot of \link{anal_prnTrend}.
#' 
#' Stop if \code{...} contains disabled argument(s) for a given \code{choice}.
#' 
#' @inheritParams anal_prnTrend
#' @param ... A character vector of argument(s) in \code{dots}.
checkdots_analTrend <- function(choice, ...) 
{
  ## error
  # anal_prnTrend(choice = cmeans, centers = 3)
  
  ## `k` later excluded as not in `formalArgs(cmeans)`
  # anal_prnTrend(choice = cmeans, k = 3)
  
  dots <- rlang::enexprs(...)
  
  msg1 <- "determined automatically."
  msg2 <- "replaced with `n_clust`."
  
  if (choice %in% c("cmeans", "kmeans")) {
    purrr::walk2(c("x", "centers"), c(msg1, msg2), ~ {
      if (.x %in% names(dots)) stop("`", .x, "` ", .y, call. = FALSE)
    })
  } 
  else if (choice %in% c("clara", "pam", "fanny")) {
    purrr::walk2(c("x", "k"), c(msg1, msg2), ~ {
      if (.x %in% names(dots)) stop("`", .x, "` ", .y, call. = FALSE)
    })
  }
}


#' Trend analysis
#' 
#' @inheritParams anal_prnTrend
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr purrr  
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom cluster clusGap
#' @importFrom magrittr %>% %T>% %$% %<>% 
analTrend <- function (df, id, col_group, col_order, label_scheme_sub, 
                       choice, n_clust,
                       scale_log2r, complete_cases, impute_na, 
                       filepath, filename, anal_type, ...) 
{
  if (!nrow(label_scheme_sub))
    stop("Empty metadata.")
  
  complete_cases <- to_complete_cases(complete_cases = complete_cases, 
                                      impute_na = impute_na)
  if (complete_cases)
    df <- my_complete_cases(df, scale_log2r, label_scheme_sub)
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)
  
  filter_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^filter_", names(.))]
  
  arrange_dots <- dots %>% 
    .[purrr::map_lgl(., is.language)] %>% 
    .[grepl("^arrange_", names(.))]
  
  dots <- dots %>% 
    .[! . %in% c(filter_dots, arrange_dots)]
  
  df_num <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    prepDM(id = !!id, scale_log2r = scale_log2r, 
           sub_grp = label_scheme_sub$Sample_ID, 
           anal_type = anal_type, 
           rm_allna = TRUE) %>% 
    .$log2R
  
  label_scheme_sub <- label_scheme_sub %>% 
    dplyr::filter(Sample_ID %in% colnames(df_num))
  
  cols_lsnum <- unlist(lapply(label_scheme_sub, is.numeric))
  cols_lsnum <- names(cols_lsnum[cols_lsnum])
  
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  
  lobal({
    grps <- label_scheme_sub[[as.character(col_group)]]
    ords <- label_scheme_sub[[as.character(col_order)]]
    grp_nas <- is.na(grps)
    ord_nas <- is.na(ords)
    
    if (!identical(grp_nas, ord_nas)) {
      nas <- tibble::tibble(
        !!col_group := grps, !!col_order := ords,
      )
      warning("Mismatches in metadata columns: \n", )
      print(nas)
    }
    
    grps <- grps[!grp_nas]
    ords <- ords[!ord_nas]
    
    oks <- tibble::tibble(
      !!col_group := grps, !!col_order := ords,
    )
    
    message("Summary of data groups and orders: \n")
    print(oks)
  })

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  df_mean <- t(df_num) %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::left_join(label_scheme_sub, by = "Sample_ID") %>%
    dplyr::filter(!is.na(!!col_group)) %>%
    dplyr::mutate(Group := !!col_group) %>% 
    dplyr::group_by(Group)
  
  df_mean <- suppressWarnings(
    df_mean %>%
      dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
  )
  
  if (rlang::as_string(col_order) %in% names(df_mean)) {
    df_mean <- df_mean %>%
      dplyr::arrange(!!col_order) %>%
      dplyr::select(-col_order) %>%
      dplyr::select(-which(names(.) %in% cols_lsnum))
  }
  
  df_mean <- df_mean %>%
    dplyr::ungroup() %>% 
    dplyr::select(-Group) %>% 
    t() %>% 
    data.frame(check.names = FALSE) %>%
    `colnames<-`(df_mean$Group) %>%
    tibble::rownames_to_column(id) %>%
    dplyr::filter(complete.cases(.[, !grepl(id, names(.))])) %>%
    `rownames<-`(NULL) %>% 
    tibble::column_to_rownames(id)
  
  if (is.null(n_clust)) {
    n_clust <- local({
      gap_stat <- cluster::clusGap(df_mean, kmeans, 10, B = 100)
      cluster::maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])
    })
    
    message("Set `n_clust` to ", n_clust, ".")
    
    fn_prefix <- paste0(fn_prefix, n_clust)
  } 
  else {
    stopifnot(all(n_clust >= 2L), all(n_clust %% 1 == 0L))
  }
  
  dots <- local({
    if (choice == "cmeans") {
      # (1) overwriting args: to this <- from `cmeans`
      if (!is.null(dots$centers)) {
        n_clust <- centers
        dots$centers <- NULL
        warning("Replace the value of `n_clust` with `centers`.")
      }
      
      # (2) excluded formalArgs: `x` and `centers`; 
      #     mistaken `k` intended for `cluster::clara` excluded as well; 
      #     no mismatch between "m" and "method" with `%in%`
      dots <- dots %>% 
        .[names(.) %in% c("iter.max", 
                          "verbose",
                          "dist",
                          "method", 
                          "m", 
                          "rate.par", 
                          "weights", 
                          "control")]
      
      # (3) conversion: expr to character string
      purrr::map(c("dist", "method"), ~ {
        if (.x %in% names(dots)) dots[[.x]] <<- rlang::as_string(dots[[.x]])
      })
      
      # (4) overwriting from this -> to `cmeans` defaults
      if (is.null(dots$m)) {
        dots$m <- local({
          N <- nrow(df_mean)
          D <- ncol(df_mean)
          m <- 1 + (1418/N + 22.05) * D^(-2) + (12.33/N + 0.243) *
            D^(-0.0406 * log(N) - 0.1134)
          
          message("Set `m` to ", round(m, digits = 2), ".")
          
          invisible(m)
        })
      }
    } 
    else if (choice == "kmeans") {
      if (!is.null(dots$centers)) {
        n_clust <- centers
        dots$centers <- NULL
        warning("Replace the value of `n_clust` with `centers`.")
      }
      
      dots <- dots %>% 
        .[names(.) %in% c("iter.max", 
                          "nstart",
                          "algorithm",
                          "trace")]
      
      purrr::map(c("algorithm"), ~ {
        if (.x %in% names(dots)) dots[[.x]] <<- rlang::as_string(dots[[.x]])
      })
      
      if (is.null(dots$nstart)) {
        dots$nstart <- 20
        message("Set `n_start` to 20.")
      }
    } 
    else if (choice == "clara") {
      if (!is.null(dots[["k"]])) {
        n_clust <- k
        dots[["k"]] <- NULL
        warning("Replace the value of `n_clust` with `k`.")
      }
      
      dots <- dots %>% 
        .[names(.) %in% c("metric", 
                          "stand",
                          "samples",
                          "sampsize", 
                          "trace", 
                          "medoids.x",
                          "keep.data",
                          "rngR", 
                          "pamLike", 
                          "correct.d")]
      
      purrr::map(c("metric"), ~ {
        if (.x %in% names(dots)) dots[[.x]] <<- rlang::as_string(dots[[.x]])
      })
      
      if (is.null(dots$samples)) {
        dots$samples <- 50
        message("Set `samples` to 50.")
      }
    } 
    else if (choice == "pam") {
      if (!is.null(dots[["k"]])) {
        n_clust <- k
        dots[["k"]] <- NULL
        warning("Replace the value of `n_clust` with `k`.")
      }
      
      dots <- dots %>% 
        .[names(.) %in% c("diss",
                          "metric", 
                          "medoids",
                          "stand",
                          "cluster.only",
                          "do.swap",
                          "keep.diss",
                          "keep.data",
                          "pamonce", 
                          "trace.lev")]
      
      purrr::map(c("metric"), ~ {
        if (.x %in% names(dots)) dots[[.x]] <<- rlang::as_string(dots[[.x]])
      })
    } 
    else if (choice == "fanny") {
      if (!is.null(dots[["k"]])) {
        n_clust <- k
        dots[["k"]] <- NULL
        warning("Replace the value of `n_clust` with `k`.")
      }
      
      dots <- dots %>% 
        .[names(.) %in% c("diss",
                          "memb.exp",
                          "metric", 
                          "stand",
                          "iniMem.p",
                          "cluster.only",
                          "keep.diss",
                          "keep.data",
                          "maxit",
                          "tol",
                          "trace.lev")]
      
      purrr::map(c("metric"), ~ {
        if (.x %in% names(dots)) dots[[.x]] <<- rlang::as_string(dots[[.x]])
      })
      
      if (is.null(dots$maxit)) {
        dots$maxit <- 50
        message("Set `maxit` to 50.")
      }
    } 
    else {
      stop("Unknown `choice = ", choice, "`.", call. = FALSE)
    }
    
    invisible(dots)
  })
  
  purrr::walk(fn_prefix, ~ {
    n_clust <- gsub("^.*_nclust(\\d+).*", "\\1", .x) %>% 
      as.numeric()
    
    filename <- paste0(.x, ".txt")
    
    args <- if (choice %in% c("cmeans", "kmeans"))
      c(list(x = as.matrix(df_mean), centers = n_clust), dots)
    else if (choice %in% c("clara", "pam", "fanny"))
      c(list(x = as.matrix(df_mean), k = n_clust), dots)
    else
      stop("Unknown `choice` ", choice, ".")
    
    if (choice == "cmeans") {
      # Warning of cmeans: partial argument match of 'length' to 'length.out'
      cl <- suppressWarnings(do.call(e1071::cmeans, args))
    } 
    else if (choice == "kmeans") {
      cl <- do.call(stats::kmeans, args)
    } 
    else if (choice == "clara") {
      cl <- do.call(cluster::clara, args)
      names(cl)[which(names(cl) == "clustering")] <- "cluster"
    } 
    else if (choice == "pam") {
      cl <- do.call(cluster::pam, args)
      names(cl)[which(names(cl) == "clustering")] <- "cluster"
    } 
    else if (choice == "fanny") {
      cl <- do.call(cluster::fanny, args)
      names(cl)[which(names(cl) == "clustering")] <- "cluster"
    } 
    else {
      stop("Unknown `choice` ", choice, ".", call. = FALSE)
    }
    
    res_cl <- data.frame(cluster = cl$cluster) |>
      `rownames<-`(rownames(df_mean)) |>
      tibble::rownames_to_column() |>
      dplyr::rename(!!id := rowname)

    Levels <- names(df_mean)
    
    df_mean %>% 
      tibble::rownames_to_column(id) %>% 
      dplyr::left_join(res_cl, by = id) %>% 
      tidyr::gather(key = variable, value = value, -id, -cluster) %>%
      dplyr::mutate(variable = factor(variable, levels = Levels)) %>%
      dplyr::arrange(variable) %>% 
      dplyr::rename(group = variable, log2FC = value) %>% 
      dplyr::left_join(df[, c(id, "species")], by = id) %>% 
      dplyr::mutate(col_group = rlang::as_string(col_group), 
                    col_order = rlang::as_string(col_order)) %T>% 
      write.table(file.path(filepath, filename), sep = "\t", 
                  col.names = TRUE, row.names = FALSE, quote = FALSE)	
  })
}


#' Plots trends
#'
#' @inheritParams plot_prnTrend
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @import dplyr purrr ggplot2 RColorBrewer
#' @importFrom tidyr gather
#' @importFrom e1071 cmeans
#' @importFrom magrittr %>% %T>% %$% %<>% 
plotTrend <- function(id, col_group, col_order, label_scheme_sub, n_clust, 
                      scale_log2r, complete_cases, impute_na, 
                      df2 = NULL, filepath, filename, theme, ...) 
{
  if (!nrow(label_scheme_sub))
    stop("Empty metadata.")
  
  id <- rlang::as_string(rlang::enexpr(id))
  dots <- rlang::enexprs(...)

  # find input df2 ---------------------------
  ins <- list.files(path = filepath, pattern = "Trend_[ONZ]_.*nclust\\d+.*\\.txt$")
  
  if (!length(ins)) 
    stop("No inputs under ", filepath, call. = FALSE)
  
  if (is.null(df2)) {
    ins <- if (impute_na) ins[grepl("_impNA", ins)] else ins[!grepl("_impNA", ins)]

    ins <- if (is.na(scale_log2r))
      ins[grepl("_Trend_O_", ins)]
    else if (scale_log2r)
      ins[grepl("_Trend_Z_", ins)]
    else
      ins[grepl("_Trend_N_", ins)]

    if (!length(ins)) 
      stop("No inputs correspond to impute_na = ", impute_na, 
           ", scale_log2r = ", scale_log2r)
    
  	if (is.null(n_clust)) {
  	  df2 <- ins
  	} 
    else {
  	  stopifnot(all(n_clust >= 2L), all(n_clust %% 1 == 0L))
  	  
  	  df2 <- local({
    	  possibles <- ins %>% 
    	    gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", .) %>% 
    	    as.numeric() %>% 
    	    `names<-`(ins)
    	  
    	  n_clust2 <- n_clust %>% .[. %in% possibles]
    	  
    	  df2 <- possibles %>% 
    	    .[. %in% n_clust2] %>% 
    	    names(.)	    
  	  })
  	  
  	  if (!length(df2)) 
  	    stop("No input files correspond to impute_na = ", impute_na, 
  	         ", scale_log2r = ", scale_log2r, 
  	         " at n_clust = ", paste0(n_clust, collapse = ","))
  	}    
  } 
  else {
    local({
      non_exists <- df2 %>% .[! . %in% ins]
      
      if (length(non_exists))
        stop("Missing trend file(s): ", paste(non_exists, collapse = ", "))
    })
    
    if (!length(df2))
      stop("File(s) not found under ", filepath)
  }

  # prepare output filename ---------------------------	
  custom_prefix <- if (id %in% c("pep_seq", "pep_seq_mod"))
    purrr::map_chr(df2, ~ gsub("(.*_{0,1})Peptide_Trend.*", "\\1", .x))
  else if (id %in% c("prot_acc", "gene"))
    purrr::map_chr(df2, ~ gsub("(.*_{0,1})Protein_Trend.*", "\\1", .x))
  else
    stop("Unknown id = ", id)

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename) %>% .[1]
  fn_prefix <- gsub("\\.[^.]*$", "", filename)

  # plot data ---------------------------
  proteoq_trend_theme <- theme_bw() + theme(
    axis.text.x  = element_text(angle=60, vjust=0.5, size=24),
    axis.ticks.x  = element_blank(), 
    axis.text.y  = element_text(angle=0, vjust=0.5, size=24),
    axis.title.x = element_text(colour="black", size=24),
    axis.title.y = element_text(colour="black", size=24),
    plot.title = element_text(face="bold", colour="black",
                              size=20, hjust=.5, vjust=.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = '#0868ac', colour = 'red'),
    
    strip.text.x = element_text(size = 24, colour = "black", angle = 0),
    strip.text.y = element_text(size = 24, colour = "black", angle = 90),
    
    plot.margin = unit(c(5.5, 55, 5.5, 5.5), "points"), 
    
    legend.key = element_rect(colour = NA, fill = 'transparent'),
    legend.background = element_rect(colour = NA,  fill = "transparent"),
    legend.position = "none",
    legend.title = element_text(colour="black", size=18),
    legend.text = element_text(colour="black", size=18),
    legend.text.align = 0,
    legend.box = NULL
  )
  
  if (is.null(theme)) 
    theme <- proteoq_trend_theme

  purrr::walk2(df2, custom_prefix, ~ {
    n <- gsub(".*_nclust(\\d+)[^\\d]*\\.txt$", "\\1", .x) %>% 
      as.numeric()
    
    out_nm <- paste0(.y, fn_prefix, "_nclust", n, ".", fn_suffix)
    src_path <- file.path(filepath, .x)

    df <- tryCatch(
      readr::read_tsv(src_path, col_types = cols(group = col_factor())), 
      error = function(e) NA)

    if (!is.null(dim(df)))
      message(paste("File loaded:", src_path))
    else
      stop(paste("Non-exist file or directory:", src_path))

    col_group <- df[["col_group"]][1]
    col_order <- df[["col_order"]][1]
    
    Levels <- label_scheme_sub %>% 
      dplyr::arrange(!!rlang::sym(col_order)) %>% 
      dplyr::select(!!rlang::sym(col_group)) %>% 
      unique() %>% 
      unlist()
    
    local({
      levs_df <- levels(df$group)
      mis_levs <- levs_df %>% .[! . %in% Levels]

      if (length(mis_levs)) {
        if (length(mis_levs) > 12) 
          mis_levs <- c(mis_levs[1:12], "...")
        
        stop("\n--- Mismatches in data levels ---\n\n", 
                "Levels in `", .x, "`:\n",
             paste(mis_levs, collapse = ", "), 
             "\n\n", 
             "Levels by `col_group = ", rlang::as_string(col_group), "`:\n", 
             paste(Levels, collapse = ", "), "\n\n", 
             "??? Check for consistency in the setting of `anal_prnTrend(col_group = ...)` ", 
             "and `plot_prnTrend(col_group = ...)` for file `", 
             .x, "`.", 
             call. = FALSE)
      }
    })
    
    if (length(dots)) {
      if (any(grepl("^filter_", names(dots))))
        stop("Primary `filter_` depreciated; use secondary `filter2_`.")
    }

    filter2_dots <- dots %>% 
      .[purrr::map_lgl(., is.language)] %>% 
      .[grepl("^filter2_", names(.))]
    
    arrange2_dots <- dots %>% 
      .[purrr::map_lgl(., is.language)] %>% 
      .[grepl("^arrange2_", names(.))]
    
    dots <- dots %>% 
      .[! . %in% c(filter2_dots, arrange2_dots)]
    
    df <- df %>% 
      dplyr::filter(group %in% Levels) %>% 
      filters_in_call(!!!filter2_dots) %>% 
      arrangers_in_call(!!!arrange2_dots) %>% 
      dplyr::mutate(group = factor(group, levels = Levels))
    
    rm(list = c("filter2_dots", "arrange2_dots"))

    if (complete_cases) 
      df <- df[complete.cases(df), ]
    
    ymin <- eval(dots$ymin, envir = rlang::caller_env())
    ymax <- eval(dots$ymax, envir = rlang::caller_env())
    ybreaks <- eval(dots$ybreaks, envir = rlang::caller_env())
    ncol <- dots$ncol
    nrow <- dots$nrow
    width <- dots$width
    height <- dots$height
    color <- dots$color
    alpha <- dots$alpha
    
    if (is.null(ymin)) ymin <- -2
    if (is.null(ymax)) ymax <- 2
    if (is.null(ybreaks)) ybreaks <- 1
    if (is.null(ncol)) ncol <- 1
    if (is.null(nrow)) nrow <- 2
    if (is.null(color)) color <- "#f0f0f0"
    if (is.null(alpha)) alpha <- .25
    
    dots$ymin <- NULL
    dots$ymax <- NULL
    dots$ybreaks <- NULL
    dots$ncol <- NULL
    dots$nrow <- NULL
    dots$color <- NULL
    dots$alpha <- NULL
    
    x_label <- expression("Ratio ("*log[2]*")")
    
    n_clust <- length(unique(df$cluster))
    if (is.null(width)) width <- n_clust * 8 / nrow + 2
    if (is.null(height)) height <- 8 * nrow

    dots$width <- NULL
    dots$height <- NULL

    p <- ggplot(data = df,
                mapping = aes(x = group, y = log2FC, group = !!rlang::sym(id))) +
      geom_line(colour = color, alpha = alpha) + 
      scale_y_continuous(limits = c(ymin, ymax), breaks = c(ymin, 0, ymax)) +
      labs(title = "", x = "", y = x_label) +
      theme
    
    p <- p + facet_wrap(~ cluster, nrow = nrow, labeller = label_value)
    
    ggsave_dots <- set_ggsave_dots(dots, c("filename", "plot", "width", "height"))
    
    rlang::quo(ggsave(filename = file.path(filepath, gg_imgname(out_nm)),
                      plot = p, 
                      width = width, 
                      height = height, 
                      !!!ggsave_dots)) %>% 
      rlang::eval_tidy()
  }, complete_cases = complete_cases)
}


#'Trend analysis of protein data
#'
#'\code{anal_prnTrend} performs unsupervised clustering of protein
#'\code{log2FC}.
#'
#'The option of \code{complete_cases} will be forced to \code{TRUE} at
#'\code{impute_na = FALSE}
#'
#'@inheritParams prnCorr_logFC
#'@inheritParams anal_prnNMF
#'@inheritParams prnHM
#'@param choice Character string; the clustering method in \code{c("cmeans",
#'  "clara", "kmeans", "pam", "fanny")}. The default is "cmeans".
#'@param n_clust Numeric vector; the number(s) of clusters that data will be
#'  divided into. At the NULL default, it will be determined by the gap method
#'  in \code{\link[cluster]{clusGap}}. The \code{n_clust} overwrites the
#'  argument \code{centers} in \code{\link[e1071]{cmeans}}.
#'@param impute_na Logical; if TRUE, data with the imputation of missing values
#'  will be used. The default is FALSE.
#'@param filepath Use system default.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements. \cr \cr Additional arguments for
#'  \link[e1071]{cmeans}, \link[stats]{kmeans}, \link[cluster]{clara},
#'  \link[cluster]{pam}. Note that \code{centers} in \link[e1071]{cmeans} or
#'  \link[stats]{kmeans} is replaced with \code{n_clust}. The same applies to
#'  \code{k} in \link[cluster]{clara} or \link[cluster]{pam}. \cr With
#'  \link[e1071]{cmeans}, \code{m} is according to Schwaemmle and Jensen if not
#'  provided; \cr \code{x} is disabled with input data being determined
#'  automatically.
#'@return Classified \code{log2FC}.
#'@import dplyr ggplot2
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'
#'@example inst/extdata/examples/prnTrend_.R
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
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend
#'  analysis and visualization \cr \code{\link{cluego}} for the visualization of
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} via
#'  \code{Cytoscape/ClueGO} \cr \code{\link{anal_pepNMF}},
#'  \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}},
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}},
#'  \code{\link{plot_prnNMFCoef}} and \code{\link{plot_metaNMF}} for NMF
#'  analysis and visualization \cr
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
anal_prnTrend <- function (col_select = NULL, col_group = NULL, col_order = NULL, 
                           choice = c("cmeans", "clara", "kmeans", "pam", "fanny"), 
                           n_clust = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, 
                           impute_na = FALSE, 
                           df = NULL, filepath = NULL, filename = NULL, ...) 
{
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call(paste0("anal", "_pepTrend"))
    } 
    else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(rlang::enexprs(...)) %>% 
        save_call(paste0("anal", "_prnTrend"))
    }
    , add = TRUE
  )
  
  old_opts <- options()
  options(warn = 1, warnPartialMatchArgs = TRUE)
  on.exit(options(old_opts), add = TRUE)
  
  check_dots(c("id", "df2", "anal_type"), ...)
  
  purrr::walk(formals(anal_prnTrend)[["choice"]] %>% eval(), 
              ~ check_formalArgs(anal_prnTrend, !!.x))
  
  choice <- rlang::enexpr(choice)
  choice <- if (length(choice) > 1L) "cmeans" else rlang::as_string(choice)
  checkdots_analTrend(choice, ...)
  
  ##############################################################
  # unlike `prnHM` that only wraps `dist` and `pheatmap`;
  # difficult to make every dummy variables in the instance
  ##############################################################
  
  dir.create(file.path(get_gl_dat_dir(), "Protein/Trend/log"), 
             recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)

  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)

  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_group = !!col_group, 
            col_order = !!col_order,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "Trend")(choice = choice, 
                                 n_clust = n_clust, 
                                 ...)
}


#'Visualization of trend results
#'
#'\code{plot_prnTrend} plots the trends of protein expressions from
#'\code{\link{anal_prnTrend}}.
#'
#'The function reads \code{Protein_Trend_[...].txt} files under the
#'\code{.../Protein/Trend} directory.
#'
#'@section \code{Protein_Trend_[...].txt}:
#'
#'  \tabular{ll}{ \strong{Key}   \tab \strong{Descrption}\cr id \tab a gene name
#'  or an acession number for protein data \cr cluster \tab a cluster ID
#'  assigned to an \code{id} \cr group \tab a name of the sample group for a
#'  \code{id} \cr log2FC \tab the mean \code{log2FC} of an \code{id} under a
#'  \code{group} at a given \code{cluster} \cr }
#'
#'@inheritParams anal_prnNMF
#'@inheritParams prnCorr_logFC
#'@inheritParams prnHist
#'@param df2 Character vector or string; the name(s) of secondary data file(s).
#'  An informatic task, i.e. \code{anal_prnTrend(...)} against a primary
#'  \code{df} generates secondary files such as
#'  \code{Protein_Trend_Z_nclust6.txt} etc. See also \code{\link{prnHist}} for
#'  the description of a primary \code{df}; \code{\link{normPSM}} for the lists
#'  of \code{df} and \code{df2}.
#'@param scale_log2r Logical; at the TRUE default, input files with
#'  \code{_Z[...].txt} in name will be used. Otherwise, files with
#'  \code{_N[...].txt} in name will be taken. An error will be thrown if no
#'  files are matched under given conditions.
#'@param impute_na Logical; at the TRUE default, input files with
#'  \code{_impNA[...].txt} in name will be loaded. Otherwise, files without
#'  \code{_impNA} in name will be taken. An error will be thrown if no files are
#'  matched under given conditions.
#'@param n_clust Numeric vector; the cluster ID(s) corresponding to
#'  \code{\link{anal_prnTrend}} for visualization. At the NULL default, all
#'  available cluster IDs will be used.
#'@param ...  \code{filter2_}: Variable argument statements for the row
#'  filtration against data in secondary file(s) of
#'  \code{[...]Protein_Trend_[...].txt}. See also \code{\link{prnGSPAHM}} for
#'  the format of \code{filter2_} statements. \cr \cr \code{arrange2_}: Variable
#'  argument statements for the row ordering against data in secondary file(s)
#'  of \code{[...]Protein_Trend_[...].txt}. See also \code{\link{prnGSPAHM}} for
#'  the format of \code{arrange2_} statements. \cr \cr Additional parameters for
#'  use in \code{plot_} functions: \cr \code{ymin}, the minimum y at \code{log2}
#'  scale; \cr \code{ymax}, the maximum y at \code{log2} scale; \cr
#'  \code{ybreaks}, the breaks in y-axis at \code{log2} scale; \cr \code{nrow},
#'  the number of rows; \cr \code{width}, the width of plot; \cr \code{height},
#'  the height of plot; \cr \code{color}, the color of lines; \cr \code{alpha},
#'  the transparency of lines.
#'@import purrr
#'
#'@example inst/extdata/examples/prnTrend_.R
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
plot_prnTrend <- function (col_select = NULL, col_order = NULL, n_clust = NULL, 
                           scale_log2r = TRUE, complete_cases = FALSE, 
                           impute_na = FALSE, 
                           df2 = NULL, filename = NULL, theme = NULL, ...) 
{
  check_dots(c("id", "anal_type", "df", "col_group", "filepath"), ...)

  dir.create(file.path(get_gl_dat_dir(), "Protein/Trend/log"), 
             recursive = TRUE, showWarnings = FALSE)

  id <- match_call_arg(normPSM, group_pep_by)
  
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), 
            length(id) == 1L)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  col_select <- rlang::enexpr(col_select)
  col_order <- rlang::enexpr(col_order)
  filename <- rlang::enexpr(filename)
  df2 <- rlang::enexpr(df2)

  reload_expts()
  
  info_anal(id = !!id, 
            col_select = !!col_select, 
            col_group = NULL, 
            col_order = !!col_order,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = NULL, 
            df2 = !!df2, 
            filepath = NULL, 
            filename = !!filename,
            anal_type = "Trend_line")(n_clust = n_clust, 
                                      theme = theme, ...)
}

