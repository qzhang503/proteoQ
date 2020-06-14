#' TMT labeling efficiency
#' 
#' @import dplyr rlang ggplot2
#' @inheritParams load_expts
#' @inheritParams splitPSM
#' @inheritParams cleanupPSM
#' @inheritParams annotPSM
#' @inheritParams normPSM
#' @examples 
#' \donttest{
#' res <- labEffPSM(
#'   fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_mm_2014_07.fasta"),
#' )
#' }
#' @export
labEffPSM <- function(group_psm_by = c("pep_seq", "pep_seq_mod"), group_pep_by = c("prot_acc", "gene"), 
                      dat_dir = NULL, expt_smry = "expt_smry.xlsx", frac_smry = "frac_smry.xlsx", 
                      fasta = NULL, entrez = NULL, 
                      pep_unique_by = "group", corrected_int = TRUE, rm_reverses = TRUE, 
                      rptr_intco = 1000, rm_craps = FALSE, rm_krts = FALSE, rm_outliers = FALSE, 
                      annot_kinases = FALSE, plot_rptr_int = TRUE, plot_log2FC_cv = TRUE, 
                      use_lowercase_aa = TRUE, ...) {
  
  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) 1)
    if (dat_dir == 1) 
      stop("Variable `dat_dir` not found; assign the working directory to the variable first.", call. = FALSE)
  } else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
    message("Variable `dat_dir` added to the Global Environment.")
  }
  
  if (is.null(fasta)) stop("Path(s) to fasta file(s) not found.", call. = FALSE)
  
  group_psm_by <- rlang::enexpr(group_psm_by)
  if (group_psm_by == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
    group_psm_by <- "pep_seq"
  } else {
    group_psm_by <- rlang::as_string(group_psm_by)
    stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  }
  
  group_pep_by <- rlang::enexpr(group_pep_by)
  if (group_pep_by == rlang::expr(c("prot_acc", "gene"))) {
    group_pep_by <- "prot_acc"
  } else {
    group_pep_by <- rlang::as_string(group_pep_by)
    stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  }
  
  dir.create(file.path(dat_dir, "PSM/cache"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/rprt_int/raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/rprt_int/mc"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/log2FC_cv/raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/log2FC_cv/purged"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/individual_mods"), recursive = TRUE, showWarnings = FALSE)
  
  if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^F[0-9]+\\.csv$"))) {
    type <- "mascot"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$"))) {
    type <- "mq"
  } else if (!purrr::is_empty(list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$"))) {
    type <- "sm"
  } else {
    stop("Unknow data type or missing data files.", call. = FALSE)
  }
  
  pep_unique_by <- rlang::as_string(rlang::enexpr(pep_unique_by))
  
  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))
  reload_expts()
  
  rmPSMHeaders()
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  load(file = file.path(dat_dir, "fraction_scheme.rda"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist = list.files(path = file.path(dat_dir, "PSM/cache"),
                        pattern = "^F[0-9]{6}_hdr_rm.csv$")
  
  if (length(filelist) == 0) stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
  if (length(filelist) != 4) stop("Nee all four PSM files using search parameters of:\n", 
                                  "1. TMT at both N-terminal and K\n", 
                                  "2. TMT at N-terminal but not K\n", 
                                  "3. TMT at K but not N-terminal\n",
                                  "4. TMT on neither N-terminal or K",
                                  call. = FALSE)
  
  df <- purrr::map(filelist, ~ {
    df <- read.delim(file.path(dat_dir, "PSM/cache", .x), sep = ',', check.names = FALSE, 
                     header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    df$dat_file <- gsub("_hdr_rm\\.csv", "", .x)
    
    r_start <- which(names(df) == "pep_scan_title") + 1
    int_end <- ncol(df) - 1
    if (int_end > r_start) df <- df[, -c(seq(r_start, int_end, 2))]
    
    if (TMT_plex == 16) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131N", "R131C", 
                     "R132N", "R132C", "R133N", "R133C", "R134N")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131N", "I131C", 
                   "I132N", "I132C", "I133N", "I133C", "I134N")
    } else if (TMT_plex == 11) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131N", "R131C")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131N", "I131C")
    } else if (TMT_plex == 10) {
      col_ratio <- c("R127N", "R127C", "R128N", "R128C", "R129N", "R129C",
                     "R130N", "R130C", "R131")
      col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                   "I130N", "I130C", "I131")
    } else if(TMT_plex == 6) {
      col_ratio <- c("R127", "R128", "R129", "R130", "R131")
      col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
    } else {
      col_ratio <- NULL
      col_int <- NULL
    }

    if (TMT_plex > 0) {
      colnames(df)[r_start:(r_start+TMT_plex-2)] <- col_ratio
      colnames(df)[(r_start+TMT_plex-1):(r_start+TMT_plex+TMT_plex-2)] <- col_int
      rm(r_start, int_end, col_ratio, col_int)
    }
    
    dat_id <- df$dat_file %>% unique()
    dat_file <- file.path(dat_dir, "PSM/cache", paste0(dat_id, "_header.txt"))
    stopifnot(length(dat_id)== 1, file.exists(dat_file))
    
    df_header <- readLines(dat_file)
    
    tmt_line <- df_header %>% 
      .[grep("Fixed modifications", .)] %>% 
      .[grep("TMT[0-9]{1}plex", .)]
    
    if (rlang::is_empty(tmt_line)) {
      tmt_type <- "neither"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(K\\)", tmt_line) & grepl("TMT[0-9]{1}plex\\s\\(N-term\\)", tmt_line)) {
      tmt_type <- "both"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(K\\)", tmt_line)) {
      tmt_type <- "k_only"
    } else if (grepl("TMT[0-9]{1}plex\\s\\(N-term\\)", tmt_line)) {
      tmt_type <- "nt_only"
    } 
    
    if (tmt_type == "nt_only") {
      df <- df %>% 
        dplyr::filter(grepl("K", .$pep_seq))
    }
    
    df <- df %>% dplyr::mutate(tmt_type = tmt_type)
    
    return(df)
  }) %>% 
    do.call(rbind, .)
  
  # convenience craps removals where their uniprot afasta names ended with "|"
  if (rm_craps) df <- df %>% dplyr::filter(!grepl("\\|.*\\|$", prot_acc))
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% filter_dots]
  
  # note pep_seq: from such as MENGQSTAAK to K.MENGQSTAAK.L
  df <- df %>% 
    dplyr::mutate(pep_len = str_length(pep_seq)) %>% 
    split(., .$dat_file, drop = TRUE) %>% 
    purrr::map(add_mascot_pepseqmod, use_lowercase_aa) %>% 
    bind_rows()
  
  df <- df %>% 
    dplyr::mutate(prot_acc_orig = prot_acc) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc)) %>% 
    annotPrn(fasta, entrez) %>% 
    dplyr::mutate(prot_acc = prot_acc_orig) %>% 
    dplyr::select(-prot_acc_orig)
  
  prot_matches_sig <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_matches_sig_new = n())
  
  prot_sequences_sig <- df %>%
    dplyr::select(!!rlang::sym(group_psm_by), !!rlang::sym(group_pep_by)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(group_psm_by))) %>% 
    dplyr::group_by(!!rlang::sym(group_pep_by)) %>%
    dplyr::summarise(prot_sequences_sig_new = n())
  
  df <- list(df, prot_matches_sig, prot_sequences_sig) %>% 
    purrr::reduce(left_join, by = group_pep_by) %>% 
    dplyr::mutate(prot_matches_sig = prot_matches_sig_new, 
                  prot_sequences_sig = prot_sequences_sig_new) %>%
    dplyr::select(-prot_matches_sig_new, -prot_sequences_sig_new)
  
  rm(prot_matches_sig, prot_sequences_sig)
  
  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc))
  
  # re-apply craps after annotation
  # 'acc_type' will be NA for entries not found in fasta
  acc_type <- unique(df$acc_type) %>% .[!is.na(.)]
  
  stopifnot(length(acc_type) == 1)
  
  if (rm_craps) {
    data(package = "proteoQ", prn_annot_crap)
    
    craps <- prn_annot_crap %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::select(acc_type) %>% 
      unlist()
    
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) df <- annotKin(df, acc_type)
  
  # `pep_start`, `pep_end` and `gene` will be used for protein percent coverage calculation
  if (!all(c("pep_start", "pep_end", "gene") %in% names(df))) df <- df %>% annotPeppos(fasta)
  
  if (!("prot_cover" %in% names(df) & length(filelist) == 1)) {
    df$prot_cover <- NULL
    
    df <- df %>% 
      calc_cover(id = !!rlang::sym(group_pep_by), 
                 fasta = seqinr::read.fasta(file.path(dat_dir, "my_project.fasta"), 
                                            seqtype = "AA", as.string = TRUE, set.attributes = TRUE))
  } 
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^pep_", names(.))), 
    df %>% dplyr::select(-grep("^pep_", names(.))), 
  )
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(grep("^prot_", names(.))), 
    df %>% dplyr::select(-grep("^prot_", names(.))), 
  )
  
  if (length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      # dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      # dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>%
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) # %>%
    # dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      dplyr::mutate(RAW_File = gsub('^(.*)\\\\(.*)\\.raw.*', '\\2', .$pep_scan_title)) %>% 
      dplyr::mutate(RAW_File = gsub("^.*File:~(.*)\\.d~?.*", '\\1', .$RAW_File)) %>% # Bruker
      dplyr::mutate(prot_acc = gsub("\\d::", "", .$prot_acc)) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  df_split <- df_split %>% split(., .$RAW_File, drop = TRUE)
  
  dat_files <- df %>% 
    dplyr::select(dat_file, tmt_type) %>% 
    dplyr::filter(!duplicated(dat_file))
  
  lab_eff <- purrr::map2_dbl(df_split, names(df_split), ~ {
    res <- .x %>% 
      dplyr::group_by(dat_file) %>% 
      dplyr::summarize(N = n()) %>% 
      dplyr::right_join(dat_files, by = "dat_file") %>% 
      dplyr::arrange(-N, dat_file) %T>% 
      write.csv(file.path(dat_dir, "PSM/cache", paste0(.y, ".csv")), row.names = FALSE)
    
    both <- res %>% 
      dplyr::filter(tmt_type == "both") %>% 
      dplyr::mutate(N = 2 *N) %>% 
      .[["N"]]
    
    k_only <- res %>% 
      dplyr::filter(tmt_type == "k_only") %>% 
      .[["N"]]
    
    nt_only <- res %>% 
      dplyr::filter(tmt_type == "nt_only") %>% 
      .[["N"]]
    
    neither <- res %>% 
      dplyr::filter(tmt_type == "neither") %>% 
      dplyr::mutate(N = 2 *N) %>% 
      .[["N"]]
    
    1 - (k_only + nt_only + neither) / (both + k_only + nt_only)
  })
}


#' Makes heat maps
#' 
#' @param df_meta A file name of meta data.
#' @inheritParams prnHM
#' @inheritParams info_anal
#' @inheritParams gspa_colAnnot
#' 
#' @examples
#' \donttest{
#' proteo_hm(
#'   df = Protein_delta.txt, 
#'   id = gene, 
#'   df_meta = hm_meta.xlsx, 
#'   filepath = file.path(dat_dir, "Protein/Heatmap"), 
#'   filename = "kin_delta.png",
#'   complete_cases = FALSE, 
#'   annot_cols = NULL, 
#'   annot_colnames = NULL, 
#'   annot_rows = c("kin_class"), 
#'   cluster_rows = FALSE, 
#'   xmin = -1, 
#'   xmax = 1, 
#'   xmargin = .1, 
#'   width = 5, 
#'   height = 12,
#'   arrange2_by = exprs(kin_class, gene), 
#' )
#' }
#' 
#' @import stringr dplyr magrittr readr readxl rlang ggplot2 RColorBrewer pheatmap
#' @export
proteo_hm <- function(df = NULL, id = NULL, df_meta = NULL, sample_ids = NULL, 
                      filepath = NULL, filename = NULL, complete_cases = FALSE, 
                      annot_cols = NULL, annot_colnames = NULL, annot_rows = NULL, 
                      xmin = -1, xmax = 1, xmargin = .1, ...) {
  
  dir.create(file.path(filepath), recursive = TRUE, showWarnings = FALSE)
  
  id <- rlang::enexpr(id)
  df <- rlang::enexpr(df)
  df_meta <- rlang::enexpr(df_meta)

  if (is.null(df)) stop("Data file `df` cannot be NULL.", call. = FALSE)
  if (is.null(df_meta)) stop("Metadata file `df_meta` cannot be NULL.", call. = FALSE)
  if (is.null(id)) stop("Column key `id` cannot be NULL.", call. = FALSE)
  if (is.null(filepath)) stop("`filepath` cannot be NULL.", call. = FALSE)
  if (is.null(filename)) stop("`filename` cannot be NULL.", call. = FALSE)

  id <- rlang::as_string(id)
  df <- rlang::as_string(df)
  df_meta <- rlang::as_string(df_meta)
  
  message("Use the default sheet name `Sheet1` for metadata Excel.")
  
  df_path <- file.path(filepath, df)
  if (file.exists(df_path)) {
    df <- readr::read_tsv(df_path)
  } else {
    stop("File not found: ", df_path, call. = FALSE)
  }
  
  df_meta_path <- file.path(filepath, df_meta)
  if (file.exists(df_meta_path)) {
    df_meta <- readxl::read_excel(df_meta_path) %>% dplyr::filter(rowSums(!is.na(.)) > 0)
    sample_ids <- df_meta$Sample_ID
  } else {
    if (is.null(sample_ids)) {
      stop("Provide sample IDs under either column `Sample_ID` in the Excel indicated by `df_meta` 
           or in the vector of `sample_ids`.", call. = FALSE)
    }
  }

  dots <- rlang::enexprs(...)
  filter2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter2_", names(.))]
  arrange2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange2_", names(.))]
  select2_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select2_", names(.))]
  dots <- dots %>% .[! . %in% c(filter2_dots, arrange2_dots, select2_dots)]
  
  # needed defaults before calling `pheatmap`
  if (is.null(dots$cluster_rows)) {
    cluster_rows <- TRUE
  } else {
    cluster_rows <- dots$cluster_rows
  }
  
  if (is.null(dots$cluster_cols)) {
    cluster_cols <- TRUE
  } else {
    cluster_cols <- dots$cluster_cols
  }
  
  if (is.null(dots$clustering_distance_rows)) {
    clustering_distance_rows <- "euclidean"
  } else {
    clustering_distance_rows <- dots$clustering_distance_rows
  }
  
  if (is.null(dots$clustering_distance_cols)) {
    clustering_distance_cols <- "euclidean"
  } else {
    clustering_distance_cols <- dots$clustering_distance_cols
  }
  
  n_color <- 500
  if (is.null(dots$breaks)) {
    color_breaks <- c(seq(from = xmin, -xmargin, length = n_color/2)[1:(n_color/2-1)],
                      seq(-xmargin, xmargin, length = 3),
                      seq(xmargin, xmax, length = n_color/2)[2:(n_color/2)])
  } else if (is.na(dots$breaks)) {
    color_breaks <- NA
  } else {
    color_breaks <- eval(dots$breaks, env = caller_env())
  }
  
  if (is.null(dots$color)) {
    mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
  } else if (is.na(dots$color)) {
    mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  } else {
    mypalette <- eval(dots$color, env = caller_env())
  }
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  
  x_label <- expression("Ratio ("*log[2]*")")

  df <- df %>%
    dplyr::mutate_at(vars(which(names(.) %in% sample_ids)), as.numeric) %>%
    dplyr::mutate_at(vars(which(names(.) %in% sample_ids)), ~ setHMlims(.x, xmin, xmax)) %>%
    dplyr::filter(!duplicated(!!rlang::sym(id)),
                  !is.na(!!rlang::sym(id)),
                  rowSums(!is.na(.[, which(names(.) %in% sample_ids)])) > 0) 
  
  df <- df %>% 
    filters_in_call(!!!filter2_dots) %>% 
    arrangers_in_call(!!!arrange2_dots)
  
  if (nrow(df) == 0) stop("No rows available after data filtratin.", call. = FALSE)
  
  df <- df %>%
    dplyr::filter(!is.na(.[[id]])) %>% 
    `rownames<-`(.[[id]])
  
  # generate `annotation_col` from keys in `annot_col`, et al.
  if (is.null(annot_cols)) {
    annotation_col <- NA
  } else {
    annotation_col <- colAnnot(annot_cols = annot_cols, sample_ids = sample_ids)
  }
  
  if (!is.null(annot_colnames) & length(annot_colnames) == length(annot_cols)) {
    colnames(annotation_col) <- annot_colnames
  }
  
  if (is.null(annot_rows)) {
    annotation_row <- NA
  } else {
    annotation_row <- df %>% dplyr::select(annot_rows) %>% data.frame(check.names = FALSE)
  }
  
  # column annotations
  if (is.null(dots$annotation_colors)) {
    annotation_colors <- setHMColor(annotation_col)
  } else if (is.na(dots$annotation_colors)) {
    annotation_colors <- NA
  } else {
    annotation_colors <- eval(dots$annotation_colors, env = caller_env())
  }
  
  if (complete_cases) {
    df_hm <- df %>%
      dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
  } else {
    df_hm <- df
  }
  
  df_hm <- df_hm %>%
    `rownames<-`(.[[id]])	%>%
    dplyr::select(which(names(.) %in% sample_ids))
  
  if (cluster_rows) {
    d <- dist(df_hm, method = clustering_distance_rows)
    d[is.na(d)] <- .5 * max(d, na.rm = TRUE)
    h <- hclust(d)
    dots$cluster_rows <- h
    rm(d, h)
  } else {
    dots$cluster_rows <- FALSE
  }
  
  if (cluster_cols) {
    d_cols <- dist(t(df_hm), method = clustering_distance_cols)
    d_cols[is.na(d_cols)] <- .5 * max(d_cols, na.rm = TRUE)
    h_cols <- hclust(d_cols)
    dots$cluster_cols <- h_cols
    # rm(d_cols, h_cols) # h_cols also for subtrees
  } else {
    dots$cluster_cols <- FALSE
  }
  
  filename <- gg_imgname(filename)
  
  # form `annotation_col` and `annotation_row` from `annot_col` and `annot_row` 
  dots <- dots %>% 
    .[! names(.) %in% c("mat", "filename", "annotation_col", "annotation_row", 
                        "clustering_distance_rows", "clustering_distance_cols", 
                        "color", "annotation_colors", "breaks")]
  
  p <- my_pheatmap(
    mat = df_hm,
    filename = file.path(filepath, filename),
    annotation_col = annotation_col,
    annotation_row = annotation_row, 
    color = mypalette,
    annotation_colors = annotation_colors,
    breaks = color_breaks,
    !!!dots
  )
}

