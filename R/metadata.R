#' Processes the metadata of TMT experiments
#'
#' @section \code{Empty samples}: Under the same TMT_Set and TMT_Channel, the
#'   Sample_ID can be non-trivial under one LCMS_Injection but Empty under
#'   another. During the \link{PSM2Pep} roll-ups, the values under Empty samples
#'   will be all NA and thus has no contributions whereas the values of the
#'   non-trivial samples will be used normally. 
#'
#' @inheritParams load_expts
#' @import dplyr tidyr purrr openxlsx
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @importFrom readxl read_excel
prep_label_scheme <- function(dat_dir = NULL, expt_smry = "expt_smry.xlsx", 
                              frac_smry = "frac_smry.xlsx") 
{
  if (is.null(dat_dir)) 
    dat_dir <- get_gl_dat_dir()
  
  ext <- parse_filename(expt_smry, dat_dir)$fn_suffix

  # (do not remove empty columns)
  label_scheme_full <- read_metadata(expt_smry, dat_dir, ext, "expt") %>% 
    check_unnamed_cols(expt_smry) %>% 
    check_empty_rows() %>% 
    check_required_cols(expt_smry) %>% 
    check_optional_cols() %>% 
    check_tmt126_row() %>% 
    check_tmt_nc()
  
  # BAD: only with MQ timsTOF missing PSM intensities
  label_scheme_full <- local({
    file <- file.path(dat_dir, "peptides.txt")
    
    if (file.exists(file)) {
      df <- readr::read_tsv(file)
      
      sample_ids <- names(df) %>% 
        .[grepl("^Experiment ", .)] %>% 
        gsub("^Experiment ", "", .)
      
      meta_ids <- label_scheme_full$Sample_ID
      missings <- setdiff(sample_ids, meta_ids)
      
      if (length(missings))
        stop("Sample IDs in MaxQuant 'peptides.txt' not found in '", expt_smry, 
                "':\n", paste(missings, collapse = "\n"))
      
      label_scheme_full <- label_scheme_full %>% 
        dplyr::mutate(Sample_ID = factor(Sample_ID, levels = sample_ids)) %>% 
        dplyr::arrange(Sample_ID) %>% 
        dplyr::mutate(Sample_ID = as.character(Sample_ID))
    }
    else
      label_scheme_full
  })

  TMT_plex <- TMT_plex(label_scheme_full) # after check_tmt_nc
  label_scheme_full <- check_tmt_plex(label_scheme_full, TMT_plex)
  TMT_levels <- TMT_levels(TMT_plex) # after check_tmt_plex
  is_lfq <- is.null(TMT_levels)
  is_tmt <- !is.null(TMT_levels)
  
  label_scheme_full <- label_scheme_full %>% 
    dplyr::mutate(TMT_Channel = check_channel_prefix(TMT_Channel, is_tmt)) %>% 
    dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File, ignore.case = TRUE)) %>% 
    dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File, ignore.case = TRUE)) %>% 
    dplyr::mutate(Reference = check_tmt_ref(Reference)) %>% 
    dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ as.numeric(.x)) %>%
    dplyr::mutate_at(vars(one_of("Peptide_Yield")), ~ round(.x, digits = 2L)) %>%
    tidyr::fill(one_of("TMT_Set", "LCMS_Injection", "RAW_File")) %>%
    tidyr::complete(TMT_Set, LCMS_Injection, TMT_Channel) 
  
  label_scheme_full <- label_scheme_full %>% 
    rm_fully_empty_tmt_sets(frac_smry, TMT_plex, dat_dir, ext) 
  
  label_scheme_full <- label_scheme_full %>% 
    dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
    dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel) %>% 
    check_metadata_integers(is_tmt, expt_smry)

  if (is_tmt) {
    label_scheme_full <- label_scheme_full %>% 
      syn_empty_channels() %>% 
      recheck_tmt_channels(TMT_plex) %>% 
      check_tmt_rawfiles() %>% 
      check_tmt_sampleids(TMT_plex) %>% 
      check_tmt_emptyids() %>% 
      check_tmt_chans_vs_levs(TMT_levels, expt_smry) %>% 
      check_tmt_mixplexes(TMT_levels, TMT_plex)
  } 
  else {
    # Empty.xxx removed during `rm_fully_empty_tmt_sets`
    label_scheme_full <- label_scheme_full %>% 
      check_dups_at_lcms_and_sid() %>% 
      check_lfq_exptraws()
  }
  
  ###
  # may check again... no NA for TMT_Set, LCMS_Injection...
  ###
  
  write_metadata(label_scheme_full, dat_dir, expt_smry, metatype = "expt")
  save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.rda"))
  simple_label_scheme(dat_dir, label_scheme_full)
}


#' Loads the information of analyte prefractionation
#'
#' @inheritParams load_expts
#' @inheritParams prnHist
#' @import dplyr purrr tidyr openxlsx
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom readxl read_excel
prep_fraction_scheme <- function(dat_dir = NULL, expt_smry = "expt_smry.xlsx", 
                                 frac_smry = "frac_smry.xlsx") 
{
  old_opts <- options()
  on.exit(options(old_opts), add = TRUE)
  options(warning.length = 5000L)
  
  tbl_mascot <- c(
    "\nMascot (similarly for MaxQuant and Spectrum Mill):\n",
    "-----------------------------------------------------------------------\n", 
    "     TMT_Set    | LCMS_Injection  |     RAW_File     |     PSM_File    \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "        1       |        1        | TMT1_Inj1_F1.raw |    F000001.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "                |                 | TMT1_Inj1_F2.raw |    F000001.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "       ...      |       ...       |        ...       |        ...      \n",
    "----------------|-----------------|------------------|-----------------\n",
    "        2       |        1        | TMT2_Inj1_F1.raw |    F000002.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n",
    "                |                 | TMT2_Inj1_F2.raw |    F000002.csv  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "       ...      |       ...       |        ...       |        ...      \n",
    "-----------------------------------------------------------------------\n"
  )
  
  tbl_mq <- c(
    "\nExamplary MaxQuant:\n",
    "-----------------------------------------------------------------------\n", 
    "     TMT_Set    | LCMS_Injection  |     RAW_File     |     PSM_File    \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "        1       |        1        | dup_msfile_1.raw |   msms_xxx.txt  \n", 
    "----------------|-----------------|------------------|-----------------\n",
    "        2       |        1        | dup_msfile_1.raw |   msms_yyy.txt  \n", 
    "----------------|-----------------|------------------|-----------------\n", 
    "      ...       |       ...       |        ...       |        ...      \n",
    "-----------------------------------------------------------------------\n"
  )
  
  tbl_lfq <- c(
    "\nExamplary LFQ:\n",
    "-------------------------------------------------------------------------------\n", 
    "    Sample_ID   |    TMT_Set    | LCMS_Injection |    RAW_File    |  Fraction  \n", 
    "----------------|---------------|----------------|----------------|------------\n", 
    "       s1       |  (automated)  |       1        |  msfile_1.raw  |      1     \n", 
    "----------------|---------------|----------------|----------------|------------\n", 
    "       s1       |               |       1        |  msfile_2.raw  |      2     \n", 
    "----------------|---------------|----------------|----------------|------------\n",
    "       ...      |               |      ...       |       ...      |     ...    \n",
    "-------------------------------------------------------------------------------\n"
  )
  
  
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  load(file.path(dat_dir, "label_scheme_full.rda"))
  TMT_plex <- TMT_plex2()
  TMT_levels <- TMT_levels(TMT_plex)
  is_lfq <- is.null(TMT_levels)
  is_tmt <- !is.null(TMT_levels)
  
  ext <- parse_filename(frac_smry, dat_dir)$fn_suffix

  if (file.exists(file.path(dat_dir, frac_smry))) {
    fraction_scheme <- update_frac_smry(frac_smry = frac_smry, 
                                        expt_smry = expt_smry, 
                                        dat_dir = dat_dir, 
                                        label_scheme_full = label_scheme_full, 
                                        ext = ext, 
                                        tbl_mascot = tbl_mascot, 
                                        tbl_lfq = tbl_lfq, 
                                        is_tmt = is_tmt)
  } 
  else {
    fraction_scheme <- make_frac_smry(frac_smry = frac_smry, 
                                      expt_smry = expt_smry, 
                                      dat_dir = dat_dir, 
                                      label_scheme_full = label_scheme_full, 
                                      tbl_mascot = tbl_mascot, 
                                      is_tmt = is_tmt)
  }
  
  # LFQ: Sample_ID in .xlsx, but not .rda
  # TMT: no Sample_ID in either
  fraction_scheme <- 
    suppressWarnings(dplyr::select(fraction_scheme, -one_of("Sample_ID")))
  
  save(fraction_scheme, file = file.path(dat_dir, "fraction_scheme.rda"))
  
  invisible(fraction_scheme)
}


#' Updates pre-existed frac_smry.xlsx.
#' 
#' @param label_scheme_full The metadata of expt_smry.
#' @param is_tmt Logical; is TMT experiments or not.
#' @param tbl_mascot A table of message for Mascot workflow.
#' @param tbl_lfq A table of message for LFQ workflow.
#' @inheritParams read_metadata
#' @inheritParams load_expts
update_frac_smry <- function (frac_smry = "frac_smry.xlsx", expt_smry = "expt_smry", 
                              dat_dir = NULL, label_scheme_full = NULL, ext = "xlsx", 
                              tbl_mascot = NULL, tbl_lfq = NULL, is_tmt = TRUE) 
{
  fraction_scheme <- read_metadata(frac_smry, dat_dir, ext, "frac")
  
  if (all(is.na(fraction_scheme$LCMS_Injection)))
    fraction_scheme$LCMS_Injection <- 1L
  
  fraction_scheme <- fraction_scheme %>% 
    dplyr::filter(!is.na(RAW_File)) %>% 
    dplyr::mutate(RAW_File = gsub("\\.raw$", "", RAW_File, ignore.case = TRUE)) %>% 
    dplyr::mutate(RAW_File = gsub("\\.d$", "", RAW_File, ignore.case = TRUE))
  
  if (any(duplicated(fraction_scheme$RAW_File)) && 
      is.null(fraction_scheme[["PSM_File"]])) {
    stop("\nDuplicated `RAW_File` names in `", frac_smry, "`:\n", 
         "This can occur with one RAW file corresponding to ", 
         "multiple PSM files \n", 
         "(e.g. searching the same RAW file with different parameter sets).\n", 
         "To distinguish, add PSM file names to column `PSM_File`:\n", 
         tbl_mascot)
  }	
  
  if (!is.null(fraction_scheme[["PSM_File"]]))
    fraction_scheme <- check_frac_multipsms(fraction_scheme, frac_smry, dat_dir)
  
  fraction_scheme <- check_exptfrac_raws(fraction_scheme, label_scheme_full, 
                                         frac_smry, dat_dir)
  
  if (is_tmt && !rlang::is_integerish(fraction_scheme$TMT_Set))
    stop("Values under ", frac_smry, "::TMT_Set need to be integers.")
  
  if (is_tmt) {
    stopifnot("TMT_Set" %in% names(fraction_scheme))
    
    fraction_scheme <- fraction_scheme %>% 
      tidyr::fill(TMT_Set, LCMS_Injection) %>% 
      dplyr::group_by(TMT_Set, LCMS_Injection) %>% 
      dplyr::mutate(Fraction = row_number()) %>% 
      dplyr::arrange(TMT_Set, LCMS_Injection, Fraction)
  } 
  else {
    # Although one-to-one between `Sample_ID` and `TMT_Set` in LFQ, 
    # the use of `TMT_Set` may be counter-intuitive and 
    # users may leave `TMT_Set` column blank. Therefore, column `Sample_ID` enforced.
    if (! "Sample_ID" %in% names(fraction_scheme)) 
      stop("Need column `Sample_ID` in ", frac_smry, " for LFQ.\n", tbl_lfq)
    
    fraction_scheme <- fraction_scheme %>% 
      tidyr::fill(Sample_ID, LCMS_Injection) %>% 
      dplyr::group_by(Sample_ID, LCMS_Injection) %>% 
      dplyr::mutate(Fraction = row_number()) %>% 
      dplyr::arrange(Sample_ID, LCMS_Injection, Fraction) %>% 
      dplyr::ungroup()
    
    # updates TMT_Set
    ls_sub <- label_scheme_full %>% 
      dplyr::select(Sample_ID, TMT_Set) %>% 
      dplyr::filter(!duplicated(Sample_ID))
    
    fraction_scheme <- fraction_scheme %>% 
      # check_exptfrac_raws(label_scheme_full, frac_smry, dat_dir) %>% 
      dplyr::select(-TMT_Set) %>% 
      dplyr::left_join(ls_sub, by = "Sample_ID") %>% 
      dplyr::select(names(fraction_scheme))
  }
  
  if (!rlang::is_integerish(fraction_scheme$LCMS_Injection)) {
    stop("Values under `frac_smry.xlsx::LCMS_Injection` need to be integers.", 
         call. = FALSE)
  }
  
  if (!rlang::is_integerish(fraction_scheme$Fraction)) {
    stop("Values under `frac_smry.xlsx::Fraction` need to be integers.", 
         call. = FALSE)
  }
  
  write_metadata(fraction_scheme, dat_dir, frac_smry, "frac")
  
  invisible(fraction_scheme)
}


#' Makes new frac_smry.xlsx.
#'
#' warning: data in a auto-generated \code{frac_smry.xlsx} currently may not fix
#' itself if they were based on wrong information from \code{expt_smry.xlsx}.
#' May need more thorough checking of \code{frac_smry.xlsx} upon any update or
#' correction of \code{expt_smry.xlsx}
#'
#' @inheritParams update_frac_smry
make_frac_smry <- function (frac_smry = "frac_smry.xlsx", expt_smry = "expt_smry", 
                            dat_dir = NULL, label_scheme_full = NULL, 
                            tbl_mascot = NULL, is_tmt = TRUE) 
{
  # in case forget to enter RAW_File names 
  # (if frac_smry can be automated, RAW_File must be present in expt_smry.xlsx)
  if (anyNA(label_scheme_full$RAW_File)) {
    stop("Some NA values in RAW_File detected.\n", 
         "Need RAW_File in ", expt_smry, " to automate ", frac_smry, ".")
  }
  
  cols <- c("TMT_Set", "LCMS_Injection", "RAW_File")
  if (!is_tmt) cols <- c("Sample_ID", cols)
  
  fraction_scheme <- label_scheme_full %>%
    tidyr::unite(tmt_lcms, c("TMT_Set", "LCMS_Injection"), remove = FALSE) %>% 
    dplyr::filter(!duplicated(tmt_lcms)) %>%
    dplyr::select(cols) %>%
    dplyr::group_by(TMT_Set, LCMS_Injection) %>%
    dplyr::mutate(Fraction = row_number())
  
  # the same RAW file can go into different searches
  # e.g. the same RAW but different TMT_Set
  if (any(duplicated(fraction_scheme$RAW_File)) && 
      is.null(fraction_scheme[["PSM_File"]])) {
    stop("\nDuplicated `RAW_File` names detected during the auto-generation of `", 
         frac_smry, "`.\n",
         "This could occur with one RAW file corresponding to multiple PSM files \n", 
         "(e.g. searching the same RAW file with different parameter sets).\n", 
         
         "To distinguish, add PSM file names to column `PSM_File`:\n", 
         tbl_mascot)
  }	
  
  write_metadata(fraction_scheme, dat_dir, frac_smry, "frac")
  
  invisible(fraction_scheme)
}


#' Checks the PSM_File in frac_smry.xlsx
#' 
#' Also made changes to the PSM_File. 
#' 
#' @param df A data frame of metadata.
#' @param filename A file name of fraction_scheme.
#' @inheritParams load_expts
check_frac_multipsms <- function (df = NULL, filename = "frac_smry.xlsx", 
                                  dat_dir = NULL) 
{
  df <- df %>% tidyr::fill(PSM_File) 
  
  psm_files <- unique(df$PSM_File) 
  
  if (!is.null(psm_files)) {
    exts <- c("txt", "csv", "ssv", "tsv")
    
    bads <- purrr::map(psm_files, 
                       ~ list.files(dat_dir, pattern = paste0(.x, ".", exts))) %>% 
      purrr::map_lgl(purrr::is_empty)

    if (any(bads)) {
      missings <- psm_files[bads]
      
      stop("PSM file (names) in `", filename, "` not found under ", dat_dir, ":\n", 
           paste(missings, collapse = "\n"))
    }
  }
  
  df <- df %>% 
    dplyr::mutate(PSM_File = gsub("\\.csv$|\\.txt$|\\.ssv$|\\.tsv$", "", PSM_File))
  
  raw_psm <- df %>% 
    tidyr::unite(RAW_File., RAW_File, PSM_File, sep = "@") 
  
  if (any(duplicated(raw_psm$RAW_File.))) {
    stop("The combination of `RAW_File` and `PSM_File` is not unique in `", 
         filename, "`.")
  }
  
  invisible(df)
}


#' Checks the consistency in RAW_File between expt_smry.xlsx and frac_smry.xlsx.
#'
#' @param label_scheme_full The metadata of expt_smry.
#' @param fraction_scheme The metadata of frac_smry.
#' @param filename A file name of fraction_scheme.
#' @inheritParams load_expts
check_exptfrac_raws <- function (fraction_scheme = NULL, label_scheme_full = NULL, 
                                 filename = "frac_smry.xlsx", dat_dir = NULL) 
{
  expt_raws <- unique(label_scheme_full$RAW_File) %>% 
    gsub("\\.raw$", "", .) %>% 
    gsub("\\.d$", "", .)
  
  frac_raws <- unique(fraction_scheme$RAW_File) %>% 
    gsub("\\.raw$", "", .) %>% 
    gsub("\\.d$", "", .)
  
  # no prefractionation
  if (!all(is.na(expt_raws))) {
    not_oks <- frac_raws %>% .[! . %in% expt_raws]
    
    if (length(not_oks)) 
      stop("Remove superfluous file(s) in `", filename, 
           "` not in `expt_smry.xlsx`:\n", paste(not_oks, collapse = ", "))
    
    not_oks_2 <- expt_raws %>% .[! . %in% frac_raws]
    
    if (length(not_oks_2)) 
      stop("File(s) in `expt_smry` not in `", filename, "`:\n", 
           paste(not_oks_2, collapse = ", "))
  }
  
  invisible(fraction_scheme)
}


#' Reads the metadata of expt_smry or frac_smry.
#' 
#' @param filename The filename of metadata.
#' @param ext The extension of a file name.
#' @param metatype The type of metadata in one of "expt" or "frac".
#' @inheritParams load_expts
read_metadata <- function (filename = "expt_smry.xlsx", dat_dir = NULL, 
                           ext = "xlsx", metatype = "expt") 
{
  stopifnot(length(filename) == 1L, 
            length(ext) == 1L, 
            length(metatype) == 1L)
  
  file <- file.path(dat_dir, filename)
  
  if (!file.exists(file))
    stop("File not found: ", file)
  
  if (ext %in% c("xls", "xlsx")) {
    if (metatype == "expt") {
      df <- tryCatch(readxl::read_excel(file, sheet = "Setup"), 
                     error = function(e) NULL)
      
      if (is.null(df))
        stop("File sheet `Setup` not found in Excel.")
    }
    else if (metatype == "frac") {
      df <- tryCatch(readxl::read_excel(file, sheet = "Fractions"), 
                     error = function(e) NULL)
      
      if (is.null(df))
        stop("File sheet `Fractions` not found in Excel.")
    }
    else {
      stop("`metatype` is not one of \"expt\" or \"frac\".")
    }
  } 
  else if (ext == "csv") {
    df <- read.csv(file, check.names = TRUE, header = TRUE, 
                   comment.char = "#", na.strings = c("", "NA")) 
  } 
  else {
    stop("File extension is not one of '.xls', '.xlsx' or '.csv': ", filename, 
         call. = FALSE)
  }
  
  df %>% 
    dplyr::filter(rowSums(!is.na(.)) > 0)
}


#' Writes the metadata of label_scheme or fraction_scheme.
#' 
#' @param df A data frame of metadata.
#' @inheritParams load_expts
#' @inheritParams read_metadata
write_metadata <- function (df, dat_dir = NULL, filename = "expt_smry.xlsx", 
                            metatype = "expt") 
{
  stopifnot(length(filename) == 1L, 
            length(metatype) == 1L, 
            is.data.frame(df))

  if (is.null(dat_dir)) 
    dat_dir <- get_gl_dat_dir()
  
  file <- file.path(dat_dir, filename)
  ext <- parse_filename(filename, dat_dir)$fn_suffix
  
  if (ext %in% c("xls", "xlsx")) {
    wb <- openxlsx::createWorkbook()
    
    sheet <- if (metatype == "expt") 
      "Setup"
    else if (metatype == "frac") 
      "Fractions"
    else 
      stop("`metatype` is not one of \"expt\" or \"frac\".")

    openxlsx::addWorksheet(wb, sheetName = sheet)
    openxlsx::writeData(wb, sheet = sheet, df)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  } 
  else if (ext == "csv") {
    write.csv(df, file, row.names = FALSE)
  } 
  else {
    stop("File extension is not one of '.xls', '.xlsx' or '.csv': ", filename, 
         call. = FALSE)
  }
  
  invisible(df)
}


#' Checks the columns without names in metadata.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param filename A name of metadata file (e.g., expt_smry.xlsx).
check_unnamed_cols <- function (df, filename = NULL) 
{
  nms <- names(df)
  nms <- nms[grepl("^\\.{3}", nms)]
  
  if (length(nms)) 
    warning("Missing column name in `", filename, "`: ", 
            paste(nms, collapse = ", "), 
            call. = FALSE)
  
  invisible(df)
}


#' Checks the presence of required columns in metadata.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param filename A name of metadata file (e.g., expt_smry.xlsx).
check_required_cols <- function (df, filename = NULL) 
{
  must_have <- c("TMT_Channel", "TMT_Set", "LCMS_Injection", "RAW_File",
                 "Sample_ID", "Reference")
  
  missing_cols <- must_have %>% .[! . %in% names(df)]
  
  if (length(missing_cols)) {
    lapply(missing_cols, function (x) {
      message(paste0("\'", x, "\' must be present in \'", filename, "\'\n"))
    })
    
    stop("Not all required columns are present in \'", filename, "\'", 
         call. = FALSE)
  }
  
  invisible(df)
}


#' Checks the default columns in metadata.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_optional_cols <- function (df) 
{
  default_names <- c("Select", "Group", "Order", "Fill",  "Color", 
                     "Shape", "Size", "Alpha", "Peptide_Yield")
  
  nms <- names(df)
  
  for (x in default_names) {
    if (! x %in% nms) {
      message("Column \'", x, "\' added to \'", filename, "\'")
      
      df[[x]] <- if (x %in% c("Order"))
        NA_integer_
      else if (x %in% c("Size", "Alpha", "Peptide_Yield"))
        NA_real_
      else if (x %in% c("Select", "Group", "Fill",  "Color", "Shape"))
        NA_character_
      else
        stop("`", x, "` is not one of the default names.", call. = FALSE)
    }
  }
  
  invisible(df)
}


#' Removes rows of all NA values
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_empty_rows <- function (df) 
{
  rows <- (rowSums(is.na(df)) < ncol(df))
  df <- df[rows, ]
  
  nms <- colnames(df)
  tier_one <- c("Sample_ID", "TMT_Channel", "TMT_Set", "LCMS_Injection", 
                "RAW_File", "Reference", "Select")
  tier_oth <- nms[!nms %in% tier_one]
  tier_one <- tier_one[tier_one %in% nms]
  
  df1 <- df[, tier_one]
  rows1 <- (rowSums(is.na(df1)) < ncol(df1))
  df1 <- df1[rows1, ]
  
  if (nrow(df1) < nrow(df)) {
    bads <- which(!rows1) %>% 
      purrr::reduce(paste, sep = ", ") %>% 
      warning("Check the \"expt_smry\" for unbalanced row(s):\n", ., ".")
  }

  df
}


#' Checks the rows TMT-126 in metadata.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_tmt126_row <- function (df) 
{
  if (all(is.na(df$TMT_Channel)))
    return(df)
  
  tmt126 <- df %>% 
    dplyr::filter(TMT_Channel == "TMT-126" | TMT_Channel == "126") %>% 
    dplyr::filter(is.na(TMT_Set) | is.na(LCMS_Injection))
  
  if (nrow(tmt126)) {
    stop("The indexes of `TMT_Set` and/or `LCMS_Injection` corresponding to \n", 
         "  the `TMT-126` rows in `expt_smry.xlsx` cannot be empty.")
  }
  
  invisible(df)
}


#' Standardizes TMT channel name.
#' 
#' For trailing N or C.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_tmt_nc <- function (df) 
{
  tmt_channels <- df$TMT_Channel
  
  if (all(is.na(tmt_channels)))
    return(df)
  
  tmt_pairs <- tibble::tibble(
    tmt127 = c("127", "127N"), 
    tmt128 = c("128", "128N"), 
    tmt129 = c("129", "129N"),
    tmt130 = c("130", "130N"), 
    tmt131 = c("131", "131N"),
  )
  
  unique_tmt <- gsub("^TMT-", "", unique(tmt_channels))
  
  for (p in tmt_pairs) {
    p <- unlist(p)
    
    if (p[1] %in% unique_tmt && p[2] %in% unique_tmt) {
      df <- df %>% 
        dplyr::mutate(TMT_Channel = gsub(paste0(p[1], "$"), p[2], TMT_Channel))
    }
  }
  
  df <- df %>% 
    dplyr::mutate(TMT_Channel = gsub("134$", "134N", TMT_Channel)) %>% 
    dplyr::mutate(TMT_Channel = gsub("126N$", "126", TMT_Channel))
  
  invisible(df)
}


#' Standardizes TMT-plexes.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param TMT_plex The multiplexity of TMT.
check_tmt_plex <- function (df, TMT_plex) 
{
  if (!TMT_plex %in% c(0L, 6L, 10L, 11L, 16L, 18L)) {
    stop("Apparent `TMT_plex = ", TMT_plex, 
         "` is not one of c(0, 6, 10, 11, 16, 18).\n", 
         "For non-TMT experiments (LFQ or simple qualitative processsing), \n", 
         "  leave the fields under `TMT_Channel` blank.\n", 
         "For other TMT-plexes, ", 
         "pad the `TMT_Channel` upstream to an available multiplicity: \n", 
         "  e.g. 1 -> 6, 8 -> 10 (regular TMT) or 9 -> 18, 15 -> 18 (TMTpro).", 
         call. = FALSE)
  }
  
  if (TMT_plex == 10L) {
    df <- df %>% 
      dplyr::mutate(TMT_Channel = gsub("126N$", "126", TMT_Channel)) %>% 
      dplyr::mutate(TMT_Channel = gsub("131N$", "131", TMT_Channel))
  } 
  else if (TMT_plex == 6L) {
    df <- df %>% 
      dplyr::mutate(TMT_Channel = gsub("(12[6-9]{1})N$", "\\1", TMT_Channel)) %>% 
      dplyr::mutate(TMT_Channel = gsub("(13[0-1]{1})N$", "\\1", TMT_Channel)) 
  }
  
  invisible(df)
}


#' Adds the \code{TMT-} prefix.
#' 
#' @param x Values under \code{TMT_Channel} in metadata.
#' @param is_tmt Logical; is TMT experiments or not.
check_channel_prefix <- function (x, is_tmt = TRUE) 
{
  if (is_tmt) {
    x <- as.character(x)
    pos <- !grepl("^TMT", x)
    x[pos] <- paste0("TMT-", x[pos])
  }

  invisible(x)
}


#' Checks the \code{Reference} channels in metadata
#'
#' @param x A vector.
check_tmt_ref <- function (x) 
{
  # (!is.na(x)) & x & (x != 0)
  (!is.na(x)) & (x != 0)
}


#' Removes fully empty TMT sets.
#'
#' If all \code{Sample_ID}'s are \code{Empty} within a \code{TMT_Set}. The whole
#' set of experiment, namely all \code{TMT_Channel}'s, will be removed.
#' Otherwise, \code{Empty} sample IDs will still be kept (in that the number of
#' Sample_IDs within a \code{TMT_Set} will be always equal to the multiplexity
#' of a TMT experiment).
#'
#' The utility is also applicable to LFQ with a special handling at TMT_plex = 0
#' (not 1).
#'
#' Sample IDs started with \code{Empty} are at first replaced with NA.
#'
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param TMT_plex The multiplexity of TMT.
#' @inheritParams load_expts
#' @inheritParams read_metadata
rm_fully_empty_tmt_sets <- function (df, frac_smry = NULL, TMT_plex = 0L, 
                                     dat_dir = NULL, ext = NULL) 
{
  df <- df %>% 
    replace_empty_with_na() %>% 
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = FALSE) %>% 
    dplyr::group_by(TMT_inj) %>% 
    dplyr::mutate(count_na. = sum(is.na(Sample_ID))) 
  
  n_bf <- nrow(df)
  
  df <- if (TMT_plex) {
    df %>% dplyr::filter(count_na. != TMT_plex)
  }
  else {
    df %>% dplyr::filter(count_na. == 0L)
  }
  
  n_af <- nrow(df)
  
  if (!n_af)
    stop("Zero zero in metadata; probably forget to set \"TMT_Set\" for LFQ.")
  
  if (n_af < n_bf) {
    file <- file.path(dat_dir, frac_smry)
    
    if (file.exists(file)) {
      fraction_scheme <- read_metadata(frac_smry, dat_dir, ext, "frac") %>% 
        tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = FALSE)
      
      rows <- fraction_scheme$TMT_inj %in% df$TMT_inj
      
      if (!all(rows)) {
        fraction_scheme <- fraction_scheme[rows, ] %>% 
          dplyr::select(-TMT_inj) %T>% 
          write_metadata(dat_dir, frac_smry, "frac")
      }
    }
  }
  
  df <- df %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-TMT_inj, -count_na.)

  invisible(df)
}


#' Replaces \code{Empty} \code{Sample_ID}'s with NA.
#' 
#' @param df A data frame of metadata.
replace_empty_with_na <- function (df) 
{
  df %>% dplyr::mutate(Sample_ID = ifelse(grepl("^Empty\\.[0-9]+", Sample_ID), 
                                          NA_character_, Sample_ID))
}


#' Replaces NA sample IDs with Empty.1 etc.
#' 
#' @param x A vector of Sample_ID.
#' @param prefix A character string prefixing \code{Empty} samples.
replace_na_with_empty <- function(x, prefix = "Empty") 
{
  nas <- is.na(x)
  seqs <- seq_len(sum(nas))
  nms <- paste0(prefix, ".", seqs)
  
  replace(as.character(x), nas, nms)
}


#' Checks the integerish of columns.
#' 
#' Used integerish (to handle NA values).
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param filename A name of metadata file (e.g., expt_smry.xlsx).
#' @param is_tmt Logical; is TMT experiments or not.
check_metadata_integers <- function (df, is_tmt, filename = "expt_smry.xlsx") 
{
  if (is_tmt && !rlang::is_integerish(df$TMT_Set))
    stop("Values under ", filename, "::TMT_Set need to be integers.")
  
  if (is_tmt) {
    if (!rlang::is_integerish(df$LCMS_Injection))
      stop("Values under ", filename, "::LCMS_Injection need to be integers.")
  }
  else {
    if (all(is.na(df$LCMS_Injection)))
      df$LCMS_Injection <- 1L
  }
  
  invisible(df)
}


#' Synchronizes Sample_ID(s) for Empty samples.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
syn_empty_channels <- function (df) 
{
  df <- replace_empty_with_na(df)
  rows <- is.na(df$Sample_ID)
  df_na <- df[rows, ]
  df_ok <- df[!rows, ]
  rm(list = c("rows"))

  # the same "TMT_Set + TMT_Channel" but different LCMS_Injection(s)
  # -> keeps temporarily the first LCMS_Injection
  #    as unique by "TMT_Set + TMT_Channel"
  emp <- df_na %>%
    dplyr::select(TMT_Set, LCMS_Injection, TMT_Channel, Sample_ID) %>% 
    tidyr::unite(Chan_Set, TMT_Channel, TMT_Set, remove = FALSE)
  
  rows <- duplicated(emp$Chan_Set)
  emp_dup <- emp[rows, ]
  emp_uni <- emp[!rows, ]
  rm(list = c("rows"))

  # TMT-128C_3_2: [TMT_Channel]_[TMT_Set]_[LCMS_Injection]
  # TMT-128N_4_1, TMT-128N_4_2, TMT-128N_4_3: 
  # if all LCMS_Injection's under the same [TMT_Channel]_[TMT_Set] are Empty 
  #   -> only keep the first (TMT-128N_4_2, _3 temporarily left out)
  emp_uni <- emp_uni %>% 
    dplyr::mutate(Sample_ID_new = replace_na_with_empty(Sample_ID, "Empty")) %>% 
    dplyr::select(Chan_Set, TMT_Set, LCMS_Injection, TMT_Channel, Sample_ID_new) %>% 
    tidyr::unite(Chan_Set_LCMS, TMT_Channel, TMT_Set, LCMS_Injection, remove = TRUE)

  # bring back TMT-128N_4_2 etc.
  emp_dup <- emp_dup %>% 
    dplyr::rename(Sample_ID_new = Sample_ID) %>% 
    tidyr::unite(Chan_Set_LCMS, TMT_Channel, TMT_Set, LCMS_Injection, remove = TRUE)
  
  stopifnot(identical(names(emp_uni), names(emp_dup)))
  
  emp <- dplyr::bind_rows(emp_uni, emp_dup) %>% 
    dplyr::group_by(Chan_Set) %>% 
    tidyr::fill(one_of("Sample_ID_new")) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(Sample_ID_new, Chan_Set_LCMS) %>% 
    dplyr::select(-Chan_Set)

  # updates Sample_ID's for Empty's
  df_na <- df_na %>% 
    tidyr::unite(Chan_Set_LCMS, TMT_Channel, TMT_Set, LCMS_Injection, remove = FALSE) %>% 
    dplyr::left_join(emp, by = "Chan_Set_LCMS") %>% 
    dplyr::mutate(Sample_ID = Sample_ID_new) %>% 
    dplyr::select(-c("Sample_ID_new", "Chan_Set_LCMS"))

  # outputs
  stopifnot(identical(names(df_ok), names(df_na)))
  
  ans <- dplyr::bind_rows(df_ok, df_na) %>% 
    dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel) %>% 
    triv_optcols_at_empty()
  
  dplyr::bind_cols(
    ans %>% dplyr::select(Sample_ID), 
    ans %>% dplyr::select(-Sample_ID))
}


#' Trivializes optional columns for Empty samples.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
triv_optcols_at_empty <- function (df) 
{
  opt_cols <- names(df) %>% 
    .[! . %in% c("TMT_Channel", "TMT_Set", "LCMS_Injection", 
                 "RAW_File", "Sample_ID", "Reference")]
  
  emps <- df %>% 
    dplyr::filter(grepl("^Empty\\.[0-9]+", Sample_ID))
  
  sids <- emps$Sample_ID
  
  if (length(sids)) {
    emps[, opt_cols] <- NA
    emps[, "Reference"] <- FALSE  	  
  }
  
  df0 <- df %>% 
    dplyr::filter(! Sample_ID %in% sids)
  
  stopifnot(identical(names(df0), names(emps)))
  
  df <- df0 %>% 
    dplyr::bind_rows(emps) %>% 
    dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)
}


#' Rechecks the completeness of TMT_Channel.
#'
#' Rows may be deleted later. Checks the completeness of TMT_Channel (RAW_File
#' is either all NA or all filled)
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param TMT_plex Numeric; the multiplexity of TMT.
recheck_tmt_channels <- function (df, TMT_plex) 
{
  check_tmt <- df %>%
    tidyr::complete(TMT_Set, LCMS_Injection, TMT_Channel) %>%
    dplyr::filter(is.na(RAW_File)) %>%
    dplyr::group_by(TMT_Set, LCMS_Injection) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(n != TMT_plex)
  
  if (nrow(check_tmt)) {
    message("=======================================")
    message("=== Mismatched TMT channel & plexes ===")
    check_tmt %>%
      dplyr::select(TMT_Set, LCMS_Injection, TMT_Channel) %>%
      data.frame(check.names = FALSE) %>% 
      print()
    message("=======================================")
    
    stop("Some parsed `TMT_plex` are different to `", TMT_plex, "`.", 
         "\n(E.g., Use `TMT-131`, instead of `TMT-131N`, ", 
         "for 10-plex experiment(s).)", 
         call. = FALSE)
  }
  
  invisible(df)
}


#' Checks the uniqueness of RAW_File per TMT_Set and LCMS_Injection.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_tmt_rawfiles <- function (df) 
{
  check_fn <- df %>%
    dplyr::group_by(TMT_Set, LCMS_Injection) %>%
    dplyr::summarise(count. = n_distinct(RAW_File)) %>%
    dplyr::filter(count. > 1L) %>%
    dplyr::select(-count.)
  
  if (nrow(check_fn)) {
    print(check_fn)
    
    stop("More than one unique RAW filename in the above combination of ", 
         "TMT sets and LCMS injections.")
  }
  
  invisible(df)
}


#' Checks the counts of Sample_ID under each TMT_Set and LCMS_Injection.
#'  
#' Should equal to the TMT_plex.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param TMT_plex Numeric; the multiplexity of TMT.
check_tmt_sampleids <- function (df, TMT_plex) 
{
  check_smpls <- df %>%
    dplyr::group_by(TMT_Set, LCMS_Injection) %>%
    dplyr::summarise(count. = n_distinct(Sample_ID)) %>%
    dplyr::filter(count. != TMT_plex)
  
  if (nrow(check_smpls) && TMT_plex) {
    print(check_smpls)
    
    stop("Need ", TMT_plex, " unique sample IDs ", 
         "in the above combination of TMT sets and LCMS injections.")
  }
  
  invisible(df)
}


#' Checks the uniqueness of Sample_ID for Empty's.
#'
#' Under the same TMT_Set but different LCMS_Injection, it is allowed that one
#' is \code{sample_1} and another is \code{Empty.xxx}.
#'
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_tmt_emptyids <- function (df) 
{
  dups <- df %>% 
    dplyr::filter(!grepl("^Empty\\.[0-9]+$", Sample_ID)) %>% 
    dplyr::group_by(TMT_Set, TMT_Channel) %>% 
    dplyr::summarise(count. = n_distinct(Sample_ID)) %>% 
    dplyr::filter(count. > 1) 
  
  if (nrow(dups) >= 2L) {
    stop("Fix the above duplication(s) in sample IDs.\n", 
         print(data.frame(dups)))
  }
  
  invisible(df)
}


#' Checks the consistency between TMT_Channel's and TMT_level's.
#'
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param TMT_levels The levels of TMT channels.
#' @param filename A name of metadata file (e.g., expt_smry.xlsx).
check_tmt_chans_vs_levs <- function (df, TMT_levels, filename = "expt_smry.xlsx") 
{
  if (!all(unique(df$TMT_Channel) %in% TMT_levels)) {
    stop("Not all TMT_Channel in ", filename, " are among \n", 
         paste(TMT_levels, collapse = ", "))
  }
  
  invisible(df)
}


#' Checks TMT_Channel after \link{check_channel_prefix} where "126" becomes
#' "TMT-126" etc.
#'
#' For example, 16-plex TMT may become 10-plex after Sample_ID removals. whereas
#' TMT_Channel may include 134 etc. The utility checks the compatibilty of
#' TMT_Channel's at mix-plexes.
#'
#' @param df A data frame of metadata (e.g., label_scheme_full).
#' @param TMT_levels The levels of TMT channels.
check_tmt_mixplexes <- function (df, TMT_levels, TMT_plex) 
{
  ls_channels <- unique(df$TMT_Channel)
  wrongs <- ls_channels %>% .[! . %in% TMT_levels]
  
  if (length(wrongs)) {
    stop("Channels not belong to TMT-", TMT_plex, ":\n", 
         paste(paste(wrongs, collapse = ", ")))
  }
  
  invisible(df)
}


#' Checks duplicated Sample_IDs in different LCMS_Injection's for LFQ.
#'
#' Supposed under a "LCMS_Injection + Sample_ID", there should be only one row.
#' If there is more than one row, likely the RAW_File is duplicated.
#'
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_dups_at_lcms_and_sid <- function (df) 
{
  dups <- df %>% 
    tidyr::unite(key, LCMS_Injection, Sample_ID, remove = FALSE) %>% 
    dplyr::filter(duplicated(key))
  
  if (nrow(dups)) {
    print(dups[, c("LCMS_Injection", "Sample_ID")])
    
    stop("Probably non-unique RAW_File in the above combination of ", 
         "LCMS_Injection and Sample_ID.", 
         call. = FALSE)
  }
  
  invisible(df)
}


#' Checks duplicated RAW_Files in expt_smry.xlsx for LFQ.
#' 
#' Also makes changes to the input \code{df}.
#' 
#' @param df A data frame of metadata (e.g., label_scheme_full).
check_lfq_exptraws <- function (df) 
{
  df <- df %>% 
    dplyr::mutate(Sample_ID = replace_na_with_empty(Sample_ID, "Empty")) %>% 
    dplyr::arrange(LCMS_Injection) %>% 
    dplyr::group_by(LCMS_Injection) %>% 
    dplyr::mutate(TMT_Set = row_number()) %>% 
    dplyr::ungroup()
  
  if (!all(is.na(df$RAW_File))) {
    # if any RAW_File is provided -> safe to remove NA RAW_File(s)
    df <- df %>% 
      dplyr::filter(!is.na(RAW_File)) 
    
    dups <- df %>% 
      dplyr::filter(duplicated(RAW_File)) %>% 
      dplyr::select(RAW_File) %>% 
      unlist()

    if (length(dups)) {
      stop("Duplicaed entries under `RAW_File`:\n", 
           paste(dups, sep = ", "))
    }
  }
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(Sample_ID), 
    df %>% dplyr::select(-Sample_ID))
  
  invisible(df)
}



#'Loads species-specific Databases
#'
#'A function loads a set of precompiled gene sets of 
#'\href{http://current.geneontology.org/products/pages/downloads.html}{GO}
#'and
#'\href{http://software.broadinstitute.org/gsea/msigdb}{molecular signatures}.
#'@seealso \code{\link{load_expts}} for supported species.
#'
#' @examples
#' \donttest{load_dbs("go_sets", "human")}
#'
#'@param species Character string; the name of a species. 
#'@inheritParams prnGSPA
#'@import dplyr
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@export
load_dbs <- function (gset_nms = NULL, species = NULL) 
{
  if (is.null(gset_nms)) 
    stop("`gset_nms` cannot be NULL.", call. = FALSE)
  if (is.null(species)) 
    stop("`species` cannot be NULL.", call. = FALSE)
  
  defaults <- c("go_sets", "kegg_sets", "c2_msig", "kinsub")
  sys_defs <- gset_nms %>% .[. %in% defaults]
  not_sys_defs <- gset_nms %>% .[! . %in% defaults]

  if (!purrr::is_empty(sys_defs)) {
    abbr_sp <- purrr::map_chr(species, sp_lookup)
    
    if (purrr::is_empty(abbr_sp)) {
      warning("Unknown species.", call. = FALSE)
    }
    
    filelist <- purrr::map(abbr_sp, ~ paste0(sys_defs, "_", .x)) %>% 
      unlist()

    suppressWarnings(data(package = "proteoQ", list = filelist, 
                          envir = environment()))
    gsets <- purrr::map(filelist, ~ try(get(.x), silent = TRUE)) %>% 
      do.call(`c`, .)
    # suppressWarnings(rm(list = filelist, envir = .GlobalEnv))
    
    if (length(gsets) > 0) names(gsets) <- gsub("/", "-", names(gsets))      
  } else {
    gsets <- NULL
  }
  
  if (!purrr::is_empty(not_sys_defs)) {
    if (!all(grepl("\\.rds$", not_sys_defs))) {
      stop("Custom gene set files indicated by `gset_nms` must ", 
           "end with extension `.rds`.", 
           call. = FALSE)
    }

    not_oks <- not_sys_defs %>% .[!file.exists(not_sys_defs)]
    if (!purrr::is_empty(not_oks)) {
      stop("File not found: \n", purrr::reduce(not_oks, paste, sep = ", \n"), 
           call. = FALSE)
    }
    
    gsets2 <- purrr::map(not_sys_defs, readRDS) %>% do.call(`c`, .)
    
    if (length(gsets2))  {
      names(gsets2) <- gsub("/", "-", names(gsets2))
    } else {
      stop("Empty data file in: \n", purrr::reduce(not_sys_defs, paste, sep = ", \n"), 
           call. = FALSE)
    }
  } else {
    gsets2 <- NULL
  }
  
  gsets <- c(gsets, gsets2) %>% .[!duplicated(.)]
  
  if (!length(gsets)) {
    warning("Zero entries in `gsets`.", call. = FALSE)
  }

  assign("gsets", gsets, envir = .GlobalEnv)
} 


#'Loads TMT or LFQ experiments
#'
#'\code{load_expts} processes \code{.xlsx} or \code{.csv} files containing the
#'metadata of TMT or LFQ experiments. For simplicity, \code{.xlsx} will be
#'assumed in the document.
#'
#'@section \code{expt_smry.xlsx}: The \code{expt_smry.xlsx} files should be
#'  placed immediately under the file folder defined by \code{dat_dir}. The tab
#'  containing the metadata of TMT or LFQ experiments should be named
#'  \code{Setup}. The \code{Excel} spread sheet therein is comprised of three
#'  tiers of fields: (1) essential, (2) optional default and (3) optional open.
#'  The \code{essential} columns contain the mandatory information of the
#'  experiments. The \code{optional default} columns serve as the fields for
#'  default lookups in sample selection, grouping, ordering, aesthetics, etc.
#'  The \code{optional open} fields allow users to define their own analysis,
#'  aesthetics, etc.
#'
#'  \tabular{ll}{ \strong{Essential column}   \tab \strong{Descrption}\cr
#'  Sample_ID \tab Unique sample IDs \cr TMT_Channel \tab TMT channel names:
#'  \code{126}, \code{127N}, \code{127C} etc. (left void for LFQ) \cr TMT_Set
#'  \tab TMT experiment indexes 1, 2, 3, ... (auto-filled for LFQ) \cr
#'  LCMS_Injection   \tab LC/MS injection indexes 1, 2, 3, ... under a
#'  \code{TMT_Set} \cr RAW_File \tab MS data file names originated by \code{MS}
#'  software(s) \cr Reference \tab Labels indicating reference samples in TMT or
#'  LFQ experiments \cr }
#'
#'  \code{Sample_ID}: values should be unique for entries at a unique
#'  combination of \code{TMT_Channel} and \code{TMT_Set}, or voided for unused
#'  entries. Samples with the same indexes of \code{TMT_Channel} and
#'  \code{TMT_Set} but different indexes of \code{LCMS_Injection} should have
#'  the same value in \code{Sample_ID}. No white space or special characters are
#'  allowed. See also posts for
#'  \href{https://proteoq.netlify.app/post/sample-exlusion-from-metadata}{sample
#'  exclusion}.
#'
#'  \code{RAW_File}: (a) for analysis with off-line fractionation of peptides
#'  before LC/MS, values under the \code{RAW_File} column should be left void.
#'  Instead, the correspondence between the fraction numbers and \code{RAW_File}
#'  names should be specified in a separate file, for example,
#'  \code{frac_smry.xlsx}. (2) For analysis without off-line fractionation, it
#'  is recommended as well to leave the field under the \code{RAW_File} column
#'  blank and instead enter the MS file names in \code{frac_smry.xlsx}.
#'
#'  The set of \code{RAW_File} names in metadata needs to be identifiable in PSM
#'  data. Impalpable mismatches might occur when \code{OS} file names were
#'  altered by MS users and thus different to those recorded internally in MS
#'  data for parsing by search engine(s). In the case, machine-generated MS file
#'  names should be used. In addition, MS files may occasionally have no
#'  contributions to PSM findings. In the case, users will be prompted to remove
#'  these MS file names.
#'
#'  Utilities \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} may
#'  aid matching MS file names between metadata and PSM data. Utility
#'  \code{\link{extract_raws}} extracts the names of MS files under a file
#'  folder. Utility \code{\link{extract_psm_raws}} extracts the names of MS
#'  files that are available in PSM data.
#'
#'  \code{Reference}: reference entry(entries) are indicated with non-void
#'  string(s).
#'
#'  \tabular{ll}{ \strong{Optional default column}   \tab \strong{Descrption}\cr
#'  Select \tab Samples to be selected for indicated analysis \cr Group \tab
#'  Aesthetic labels annotating the prior knowledge of sample groups, e.g.,
#'  Ctrl_T1, Ctrl_T2, Disease_T1, Disease_T2, ...\cr Order \tab Numeric labels
#'  specifying the order of sample \code{groups} \cr Fill \tab Aesthetic labels
#'  for sample annotation by filled color\cr Color \tab Aesthetic labels for
#'  sample annotation by edge color\cr Shape \tab Aesthetic labels for sample
#'  annotation by shape\cr Size \tab Aesthetic labels for sample annotation by
#'  size \cr Alpha \tab Aesthetic labels for sample annotation by transparency
#'  \cr \cr}
#'
#'  \tabular{ll}{ \strong{Exemplary optional open column}   \tab
#'  \strong{Descrption}\cr Term \tab Categorical terms for statistical modeling.
#'  \cr Peptide_Yield \tab Yields of peptides in sample handling \cr}
#'
#'
#'@section \code{frac_smry.xlsx}: \tabular{ll}{ \strong{Column}   \tab
#'  \strong{Descrption}\cr Sample_ID \tab Unique sample IDs (only required with
#'  LFQ) \cr TMT_Set \tab TMT experiment indexes (auto-filled for LFQ) \cr
#'  LCMS_Injection \tab LC/MS injection indexes \cr Fraction \tab Fraction
#'  indexes under a \code{TMT_Set} \cr RAW_File \tab MS data file names \cr
#'  PSM_File \tab Names of PSM files. Required only when one \code{RAW_File} can
#'  be linked to multiple PSM files (e.g. F012345.csv and F012346.csv both from
#'  ms_1.raw). }
#'
#'@family normalization functions
#'@seealso \emph{Data normalization} \cr \code{\link{normPSM}} for extended
#'  examples in PSM data normalization \cr \code{\link{PSM2Pep}} for extended
#'  examples in PSM to peptide summarization \cr \code{\link{mergePep}} for
#'  extended examples in peptide data merging \cr \code{\link{standPep}} for
#'  extended examples in peptide data normalization \cr \code{\link{Pep2Prn}}
#'  for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization.
#'  \cr \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples
#'  in data purging \cr \code{\link{pepHist}} and \code{\link{prnHist}} for
#'  extended examples in histogram visualization. \cr \code{\link{extract_raws}}
#'  and \code{\link{extract_psm_raws}} for extracting MS file names \cr
#'
#'@family data row filtration
#'@seealso \emph{User-friendly utilities for variable arguments of `filter_...`}
#'  \cr \code{\link{contain_str}}, \code{\link{contain_chars_in}},
#'  \code{\link{not_contain_str}}, \code{\link{not_contain_chars_in}},
#'  \code{\link{start_with_str}}, \code{\link{end_with_str}},
#'  \code{\link{start_with_chars_in}} and \code{\link{ends_with_chars_in}} for
#'  data subsetting by character strings \cr
#'
#'@family missing value imputation
#'@seealso \emph{Missing values} \cr \code{\link{pepImp}} and
#'  \code{\link{prnImp}} for missing value imputation \cr
#'
#'@family basic informatics
#'@seealso \emph{Informatics} \cr \code{\link{pepSig}} and \code{\link{prnSig}}
#'  for significance tests \cr \code{\link{pepVol}} and \code{\link{prnVol}} for
#'  volcano plot visualization \cr \code{\link{prnGSPA}} for gene set enrichment
#'  analysis by protein significance pVals \cr \code{\link{gspaMap}} for mapping
#'  GSPA to volcano plot visualization \cr \code{\link{prnGSPAHM}} for heat map
#'  and network visualization of GSPA results \cr \code{\link{prnGSVA}} for gene
#'  set variance analysis \cr \code{\link{prnGSEA}} for data preparation for
#'  online GSEA. \cr \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS
#'  visualization \cr \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA
#'  visualization \cr \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA
#'  visualization \cr \code{\link{pepHM}} and \code{\link{prnHM}} for heat map
#'  visualization \cr \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}},
#'  \code{\link{pepCorr_logInt}} and \code{\link{prnCorr_logInt}}  for
#'  correlation plots \cr \code{\link{anal_prnTrend}} and
#'  \code{\link{plot_prnTrend}} for trend analysis and visualization \cr
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}},
#'  \code{\link{plot_pepNMFCon}}, \code{\link{plot_prnNMFCon}},
#'  \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr
#'
#'  \emph{Custom databases} \cr \code{\link{Uni2Entrez}} for lookups between
#'  UniProt accessions and Entrez IDs \cr \code{\link{Ref2Entrez}} for lookups
#'  among RefSeq accessions, gene names and Entrez IDs \cr \code{\link{prepGO}}
#'  for
#'  \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'   ontology}} \cr \code{\link{prepMSig}} for
#'  \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular
#'   signatures} \cr \code{\link{prepString}} and \code{\link{anal_prnString}}
#'  for STRING-DB \cr
#'
#'  \emph{Workflow scripts} \cr # TMT \cr system.file("extdata",
#'  "workflow_tmt_base.R", package = "proteoQ") \cr system.file("extdata",
#'  "workflow_tmt_ext.R", package = "proteoQ") \cr
#'
#'  # LFQ \cr system.file("extdata", "workflow_lfq_base.R", package = "proteoQ")
#'  \cr
#'
#'  \emph{Metadata files} \cr # TMT, no fractionation --- OK without
#'  `frac_smry.xlsx` \cr # (a. no references) \cr system.file("extdata",
#'  "expt_smry_no_prefrac.xlsx", package = "proteoQDA") \cr # (b. W2 and W16
#'  references) \cr system.file("extdata",
#'  "expt_smry_no_prefrac_ref_w2_w16.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation \cr # (a. no references) \cr system.file("extdata",
#'  "expt_smry_gtmt.xlsx", package = "proteoQDA") \cr system.file("extdata",
#'  "frac_smry_gtmt.xlsx", package = "proteoQDA") \cr
#'
#'  # (b. W2 references) \cr system.file("extdata", "expt_smry_ref_w2.xlsx",
#'  package = "proteoQDA") \cr system.file("extdata", "frac_smry_gtmt.xlsx",
#'  package = "proteoQDA") \cr
#'
#'  # (c. W2 and W16 references) \cr system.file("extdata",
#'  "expt_smry_ref_w2_w16.xlsx", package = "proteoQDA") \cr
#'  system.file("extdata", "frac_smry_gtmt.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation (global + phospho) \cr system.file("extdata",
#'  "expt_smry_tmt_cmbn.xlsx", package = "proteoQDA") \cr system.file("extdata",
#'  "frac_smry_tmt_cmbn.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation, one MS to multiple PSM files \cr
#'  system.file("extdata", "expt_smry_psmfiles.xlsx", package = "proteoQDA") \cr
#'  system.file("extdata", "frac_smry_psmfiles.xlsx", package = "proteoQDA") \cr
#'
#'  # TMT, prefractionation, mixed-plexes \cr # (column PSM_File needed; as with
#'  this example, \cr #  mixed-plexes results are actually from the same MS
#'  files \cr #  but searched separately at 6- and 10-plex settings!) \cr
#'  system.file("extdata", "expt_smry_mixplexes.xlsx", package = "proteoQDA")
#'  \cr system.file("extdata", "frac_smry_mixplexes.xlsx", package =
#'  "proteoQDA") \cr
#'
#'  # LFQ, prefractionation \cr system.file("extdata", "expt_smry_plfq.xlsx",
#'  package = "proteoQDA") \cr system.file("extdata", "frac_smry_plfq.xlsx",
#'  package = "proteoQDA") \cr
#'
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@family Monoisotopic masses of peptides
#'@seealso \emph{MS1 peptide masses} \cr \code{\link[mzion]{calc_pepmasses}}
#'  for mono-isotopic masses of peptides from fasta databases \cr
#'  \code{\link[mzion]{calc_monopeptide}} for mono-isotopic masses of peptides
#'  from individual sequences \cr \code{\link[mzion]{parse_unimod}} for
#'  parsing \href{https://www.unimod.org/}{Unimod} fixed modifications, variable
#'  modifications and neutral losses. \cr \code{\link[mzion]{find_unimod}} for
#'  finding a Unimod
#'
#'
#'
#'@param dat_dir A character string to the working directory. The default is to
#'  match the value under the global environment.
#'@param expt_smry A character string to a \code{.xlsx} file containing the
#'  metadata of TMT or LFQ experiments. The default is \code{expt_smry.xlsx}.
#'@param frac_smry A character string to a \code{.xlsx} file containing peptide
#'  fractionation summary. The default is \code{frac_smry.xlsx}.
#'
#'@example inst/extdata/examples/load_expts_.R
#'
#'@import dplyr fs
#'@rawNamespace import(rlang, except = c(list_along, invoke, flatten_raw,
#'  modify, as_function, flatten_dbl, flatten_lgl, flatten_int, flatten_chr,
#'  splice, flatten, prepend, "%@%", ":="))
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@export
load_expts <- function (dat_dir = NULL, expt_smry = "expt_smry.xlsx", 
                        frac_smry = "frac_smry.xlsx") 
{
  on.exit(mget(names(formals()), envir = environment(), inherits = FALSE) %>% 
            save_call(fun))
  
  this_call <- match.call()
  fun <- as.character(this_call[[1]])

  expt_smry <- rlang::as_string(rlang::enexpr(expt_smry))
  frac_smry <- rlang::as_string(rlang::enexpr(frac_smry))

  dat_dir <- set_dat_dir(dat_dir)

  prep_label_scheme(dat_dir, expt_smry, frac_smry)
  prep_fraction_scheme(dat_dir, expt_smry, frac_smry)
  
  session_info <- sessionInfo()
  dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)
  try(save(session_info, file = file.path(dat_dir, "Calls", "proteoQ.rda")))
}


#' Reloads the "expt_smry.xlsx" and "frac_smry.xlsx"
#'
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @importFrom fs file_info
reload_expts <- function() 
{
  dat_dir <- get_gl_dat_dir()
  
  expt_smry <- match_call_arg(load_expts, expt_smry)
  frac_smry <- match_call_arg(load_expts, frac_smry)
  
  fi_xlsx <- fs::file_info(file.path(dat_dir, expt_smry))$change_time
  
  if (is.na(fi_xlsx)) 
    stop("Time stamp of ", expt_smry, " not available.")
  
  fi_rda <- 
    fs::file_info(file.path(dat_dir, "label_scheme_full.rda"))$change_time
  
  if (fi_xlsx > fi_rda) 
    load_expts(dat_dir = dat_dir, 
               expt_smry = !!expt_smry, 
               frac_smry = !!frac_smry)
}


#' Extracts the channel information in TMT experiments
#'
#' A function returns the indexes of TMT channels that are associated to
#' reference(s), sample(s) and probable unused void(s).
#'
#' @param set_idx Numeric.  The index of a multiplex TMT experiment in metadata
#'   files such as \code{label_scheme.xlsx} and \code{frac_scheme.xlsx}.
#' @param injn_idx Numeric. The index of \code{LCMS_Inj} in metadata files such
#'   as \code{label_scheme.xlsx} and \code{frac_scheme.xlsx}.
#' @inheritParams load_expts
#' @return Three lists of indexes: \code{refChannels}, reference channels(s);
#'   \code{emptyChannels}, empty channel(s) that were not used for sample
#'   labeling; \code{labeledChannels}, non-empty channels including both
#'   reference(s) and sample(s).
#' 
#' @importFrom dplyr select filter
channelInfo <- function (dat_dir = NULL, set_idx = NULL, injn_idx = 1) 
{
  stopifnot(length(set_idx) == 1L)
  
  if (is.null(set_idx)) 
    stop("Need to specify `set_idx`.", call. = FALSE)
  if (is.null(dat_dir)) 
    stop("Need to specify `dat_dir`.", call. = FALSE)
  
  load(file = file.path(dat_dir, "label_scheme_full.rda"))
  
  label_scheme_sub <- label_scheme_full %>%
    dplyr::filter(TMT_Set == set_idx, LCMS_Injection == injn_idx)
  
  if (nrow(label_scheme_sub) == 0) 
    stop("No samples at TMT_Set ", set_idx, " and LCMS_Injection ", injn_idx, ".", 
         call. = FALSE)
  
  ref <- label_scheme_sub$Reference
  
  empty_channels <- is.na(label_scheme_sub$Sample_ID) |
    grepl("^Empty\\.[0-9]+$", label_scheme_sub$Sample_ID, ignore.case = TRUE)
  
  labeled_channels <- !empty_channels
  
  out <- list(
    refChannels = ref,
    emptyChannels = empty_channels,
    labeledChannels = labeled_channels
  )
  
  lapply(out, which)
}


#' Finds the number of multiplex TMT experiments
#'
#' @param label_scheme_full The label_scheme with the probable inclusion of
#'   different LCMS_inj under the same TMT_Set.
n_TMT_sets <- function (label_scheme_full) {
  length(unique(label_scheme_full$TMT_Set))
}


#' Finds the maximum number of LCMS injections under the same TMT_Set.
#'
#' @inheritParams n_TMT_sets
n_LCMS <- function (label_scheme_full) 
{
  label_scheme_full %>% 
    dplyr::select(TMT_Set, LCMS_Injection) %>% 
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = FALSE) %>% 
    dplyr::filter(!duplicated(TMT_inj)) %>% 
    dplyr::group_by(TMT_Set) %>% 
    dplyr::mutate(n_LCMS = n()) %>% 
    dplyr::select(TMT_Set, n_LCMS) %>% 
    dplyr::filter(!duplicated(TMT_Set)) %>% 
    dplyr::ungroup() 
}


#' Finds the multiplexity of TMT labels
#'
#' \code{TMT_plex} returns the multiplexity of TMT labels.
#' @inheritParams n_TMT_sets
TMT_plex <- function (label_scheme_full) {
  nlevels(as.factor(label_scheme_full$TMT_Channel))
}


#' Finds the multiplexity of TMT labels from previously made
#' label_scheme_full.rda
#'
#' \code{TMT_plex} returns the multiplexity of TMT labels.
TMT_plex2 <- function () 
{
  load(file.path(get_gl_dat_dir(), "label_scheme_full.rda"))
  TMT_plex(label_scheme_full)
}


#' Finds the factor levels of TMT labels
#'
#' \code{TMT_levels} returns the factor levels of TMT labels.
#' @param TMT_plex Numeric; the multiplexity of TMT, i.e., 10, 11 etc.
TMT_levels <- function (TMT_plex) 
{
  TMT_levels <- if (TMT_plex == 18L) 
    c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
      "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
      "TMT-131N", "TMT-131C", "TMT-132N", "TMT-132C", 
      "TMT-133N", "TMT-133C", "TMT-134N", "TMT-134C", "TMT-135N")
  else if (TMT_plex == 16L) 
    c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
      "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
      "TMT-131N", "TMT-131C", "TMT-132N", "TMT-132C", 
      "TMT-133N", "TMT-133C", "TMT-134N")
  else if (TMT_plex == 11L) 
    c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
      "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", 
      "TMT-131N", "TMT-131C")
  else if (TMT_plex == 10L) 
    c("TMT-126", "TMT-127N", "TMT-127C", "TMT-128N", "TMT-128C", 
      "TMT-129N", "TMT-129C", "TMT-130N", "TMT-130C", "TMT-131")
  else if (TMT_plex == 6L) 
    c("TMT-126", "TMT-127", "TMT-128", "TMT-129", "TMT-130", "TMT-131")
  else if (TMT_plex == 1L) 
    c("TMT-126")
  else if (TMT_plex == 0L) 
    NULL
}


#' Simplifies label schemes from \code{label_scheme_full}
#'
#' Removes duplicated sample entries under different LC/MS injections.
#' 
#' @inheritParams load_expts
#' @inheritParams n_TMT_sets
simple_label_scheme <- function (dat_dir, label_scheme_full) 
{
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  # OK Sample_ID of "sample_x" and "Empty.1" at the same TMT_Set and TMT_Channel,
  # but different LCMS_Injection
  label_scheme <- label_scheme_full %>%
    dplyr::filter(!duplicated(Sample_ID), !is.na(Sample_ID)) %>%
    dplyr::mutate(TMT_Channel = factor(TMT_Channel, levels = TMT_levels)) %>%
    dplyr::arrange(TMT_Set, LCMS_Injection, TMT_Channel)
  
  multi_dips <- label_scheme %>% 
    tidyr::unite(key, TMT_Channel, TMT_Set, remove = FALSE) %>% 
    dplyr::group_by(key) %>% 
    dplyr::mutate(n = n()) %>% 
    dplyr::filter(n > 1L) %>% 
    dplyr::ungroup()
  
  label_scheme <- label_scheme %>% 
    dplyr::filter(! .data$Sample_ID %in% multi_dips$Sample_ID)
  
  label_scheme <- multi_dips %>% 
    dplyr::filter(!grepl("^Empty\\.[0-9]+$", Sample_ID)) %>% 
    dplyr::select(names(label_scheme)) %>% 
    dplyr::bind_rows(label_scheme) %>% 
    dplyr::mutate(LCMS_Injection = 1L) %>% 
    dplyr::arrange(TMT_Set, TMT_Channel)
  
  if (nrow(label_scheme) < (TMT_plex * n_TMT_sets(label_scheme))) 
    stop("Duplicated sample ID(s) in `expt_smry.xlsx`", call. = FALSE)
  
  save(label_scheme, file = file.path(dat_dir, "label_scheme.rda"))
  
  invisible(label_scheme)
}


#' Finds the TMT-plex from Mascot PSM exports
#' 
#' @param df A PSM export from Mascot
#' @param pep_seq_rows The row(s) contain character string "pep_seq".
find_mascot_tmtplex <- function(df, pep_seq_rows = NULL) 
{
  if (is.null(pep_seq_rows)) 
    pep_seq_rows <- grep("pep_seq", df)
  
  if (!length(pep_seq_rows)) 
    stop("The row of `pep_seq` not found.", call. = FALSE)
  
  first_line <- df[pep_seq_rows[1] + 1]
  
  mascot_tmtplex <- if (grepl("\"135C\"", first_line, fixed = TRUE) || 
      grepl("\"135\"", first_line, fixed = TRUE)) 
    18L 
  else if (grepl("\"134N\"", first_line, fixed = TRUE) || 
      grepl("\"134\"", first_line, fixed = TRUE)) 
    16L
  else if (grepl("\"131C\"", first_line, fixed = TRUE)) 
    11L
  else if (grepl("\"131\"", first_line, fixed = TRUE) && 
             grepl("\"130C\"", first_line, fixed = TRUE)) 
    10L
  else if (grepl("\"131\"", first_line, fixed = TRUE)) 
    6L
  else 
    0L
}

