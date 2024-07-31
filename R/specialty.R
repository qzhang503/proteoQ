#' TMT labeling efficiency
#' 
#' @import dplyr ggplot2
#' @inheritParams load_expts
#' @inheritParams splitPSM_ma
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
labEffPSM <- function(group_psm_by = c("pep_seq", "pep_seq_mod"), 
                      group_pep_by = c("prot_acc", "gene"), 
                      dat_dir = NULL, expt_smry = "expt_smry.xlsx", 
                      frac_smry = "frac_smry.xlsx", 
                      fasta = NULL, entrez = NULL, 
                      pep_unique_by = "group", corrected_int = TRUE, 
                      rm_reverses = TRUE, 
                      rptr_intco = 1000, rm_craps = FALSE, rm_krts = FALSE, 
                      rm_outliers = FALSE, 
                      annot_kinases = FALSE, plot_rptr_int = TRUE, 
                      plot_log2FC_cv = TRUE, use_lowercase_aa = TRUE, ...) 
{
  if (is.null(dat_dir)) {
    dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv), error = function(e) NA)
    
    if (is.na(dat_dir)) 
      stop("Variable `dat_dir` not found; ", 
           "assign the working directory to the variable first.")
  } 
  else {
    assign("dat_dir", dat_dir, envir = .GlobalEnv)
    message("Variable `dat_dir` added to the Global Environment.")
  }
  
  if (is.null(fasta)) 
    stop("Path(s) to fasta file(s) not found.")
  
  group_psm_by <- rlang::enexpr(group_psm_by)
  
  if (group_psm_by == rlang::expr(c("pep_seq", "pep_seq_mod"))) {
    group_psm_by <- "pep_seq"
  } 
  else {
    group_psm_by <- rlang::as_string(group_psm_by)
    stopifnot(group_psm_by %in% c("pep_seq", "pep_seq_mod"))
  }
  
  group_pep_by <- rlang::enexpr(group_pep_by)
  
  if (group_pep_by == rlang::expr(c("prot_acc", "gene"))) {
    group_pep_by <- "prot_acc"
  } 
  else {
    group_pep_by <- rlang::as_string(group_pep_by)
    stopifnot(group_pep_by %in% c("prot_acc", "gene"))
  }
  
  dir.create(file.path(dat_dir, "PSM/cache"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/rprt_int/raw"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/rprt_int/mc"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/log2FC_cv/raw"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/log2FC_cv/purged"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "PSM/individual_mods"), 
             recursive = TRUE, showWarnings = FALSE)
  
  if (length(list.files(path = file.path(dat_dir), pattern = "^F[0-9]+\\.csv$"))) {
    type <- "mascot"
  } 
  else if (length(list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$"))) {
    type <- "mq"
  } 
  else if (length(list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$"))) {
    type <- "sm"
  } 
  else {
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
  filepath <- file.path(dat_dir, "PSM/cache")
  filelist = list.files(path = filepath, pattern = "^F[0-9]{6}_hdr_rm.csv$")
  
  if (!length(filelist)) 
    stop(paste("No PSM files under", file.path(dat_dir, "PSM")))
  
  lapply(filelist, procTMT0) %>% unlist()
}


#' Helper in reading TMTzero
#' 
#' @param file A file name.
procTMT0 <- function (file)
{
  read.delim(file.path(dat_dir, "PSM/cache", file), sep = ',', 
                   check.names = FALSE, header = TRUE, stringsAsFactors = FALSE, 
                   quote = "\"",fill = TRUE , skip = 0) %>% 
    add_mascot_raw() %>% 
    dplyr::rename(raw_file = RAW_File) %>% 
    split(.$raw_file) %>% 
    lapply(calcTMTLabs, file)
}


#' Calculates labeling efficiency
#' 
#' @param df A PSM table
#' @param file A file name.
calcTMTLabs <- function (df, file)
{
  df <- df %>% 
    filter(!is.na(prot_family_member), 
           pep_rank == 1L) %>% 
    mutate(pep_scan_num = gsub(".* scan\\=(.*)~$", "\\1", pep_scan_title)) %>% 
    tidyr::unite(id, raw_file, pep_scan_num, sep = "@") %>% 
    group_by(id) %>% 
    filter(row_number() == 1L) %>% 
    ungroup()
  
  df <- df %>% 
    filter(!grepl("Acetyl (Protein N-term)", pep_var_mod, fixed = TRUE))
  
  df$dat_file <- gsub("_hdr_rm\\.csv", "", file)
  dat_id <- unique(df$dat_file)
  dat_file <- file.path(dat_dir, "PSM/cache", paste0(dat_id, "_header.txt"))
  stopifnot(length(dat_id)== 1, file.exists(dat_file))
  
  var_mods <- readLines(dat_file) %>% find_mascot_vmods()
  desc <- var_mods$Description
  idx_k <- var_mods$Mascot_abbr[grep("TMT[^ ]* \\(K\\)", desc)]
  idx_nt <- var_mods$Mascot_abbr[grep("TMT[^ ]* \\(N-term\\)", desc)]
  # other_nt <- var_mods$Mascot_abbr[grep("Acetyl (Protein N-term)", desc, fixed = TRUE)]
  rm(list = c("desc", "var_mods"))
  
  if (!length(idx_k))
    stop("TMT (K) need to be a variable modification.")
  
  if (!length(idx_nt))
    stop("TMT (N-term) need to be a variable modification.")
  
  dfNT1 <- df %>% filter(grepl(idx_nt, pep_var_mod_pos)) # 127539, +NT_lab
  dfNT1rK1 <- dfNT1 %>% filter(grepl("K", pep_seq)) # 78237, +NT_lab, with residue K
  dfNT1rK1_0 <- dfNT1rK1 %>% filter(!grepl(idx_k, pep_var_mod_pos)) # 23, +NT_lab, -K_lab
  n_k0 <- nrow(dfNT1rK1_0) # 23
  
  # No N-term; 
  # may need to remove Acetyl (Protein N-term)
  dfNT0 <- df %>% filter(!grepl(idx_nt, pep_var_mod_pos)) # 2780
  dfNT0_rK0 <- dfNT0 %>% filter(!grepl("K", pep_seq)) # 871, -NT_lab, no residue K -> done
  dfNT0_rK1 <- dfNT0 %>% filter(grepl("K", pep_seq)) # 1909, -NT_lab, with residue K
  dfNT0_rK1_0 <- dfNT0_rK1 %>% filter(!grepl(idx_k, pep_var_mod_pos)) # 54,   -NT_lab, -K_lab
  dfNT0_rK1_1 <- dfNT0_rK1 %>% filter( grepl(idx_k, pep_var_mod_pos)) # 1855, -NT_lab, +K_lab
  n_nt0 <- nrow(dfNT0_rK0) + nrow(dfNT0_rK1_0) * 2 + nrow(dfNT0_rK1_1) # 871 + 2 * 54 + 1855 -> 2834
  
  1 - (n_k0 + n_nt0)/nrow(df)
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
#' @import stringr dplyr readr readxl ggplot2 RColorBrewer pheatmap
#' @importFrom magrittr %>% %T>% %$% %<>% 
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
    
    color_breaks <- color_breaks %>% 
      .[. >= xmin] %>% 
      unique()
  } else if (is.na(dots$breaks)) {
    color_breaks <- NA
  } else {
    color_breaks <- eval(dots$breaks, envir = rlang::caller_env())
  }
  
  if (is.null(dots$color)) {
    mypalette <- colorRampPalette(c("blue", "white", "red"))(n_color)
  } else if (is.na(dots$color)) {
    mypalette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  } else {
    mypalette <- eval(dots$color, envir = rlang::caller_env())
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
    annotation_colors <- eval(dots$annotation_colors, envir = rlang::caller_env())
  }
  
  if (complete_cases) {
    df_hm <- df %>%
      dplyr::filter(complete.cases(.[, names(.) %in% sample_ids]))
  } else {
    df_hm <- df
  }
  
  df_hm <- data.frame(df_hm, check.names = FALSE) %>%  
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

