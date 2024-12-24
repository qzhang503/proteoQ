#' Finds mod_indexes.
#' 
#' @inheritParams load_expts
#' @return A named vector with modifications in names and indexes in values.
find_mod_indexesQ <- function (dat_dir) 
{
  rda <- file.path(dat_dir, "Calls", "ms2match.rda")
  
  # single psmQ.txt with mzion and proteoQ under the same folder
  if (file.exists(rda)) {
    load(rda)
    mod_indexes <- call_pars$mod_indexes
    
    if (is.null((mod_indexes)))
      warning("Corrupted: ", rda)
    else
      return(mod_indexes)
  }
  
  message("Deduce `mod_indexes` from psmQ[...].txt.")
  
  files <- suppressMessages(find_psmQ_files(dat_dir))
  mods <- lapply(files, deduce_mod_indexes, dat_dir)
  
  if (!identical(mods, unique(mods))) {
    stop("Not all sets of deduced `mod_indexes` are identical.\n", 
         "Probably different sets of variable modifications in mzion searches.")
  }
  
  mods[[1]]
}


#' Deduces mod_indexes.
#' 
#' Not applicable to fixed modifications or terminal modifications.
#' 
#' @param file A psmQ file.
#' @inheritParams load_expts
deduce_mod_indexes <- function (file = "psmQ.txt", dat_dir) 
{
  df <- suppressWarnings(
    readr::read_tsv(file.path(dat_dir, file), 
                    col_types = get_col_types(), 
                    show_col_types = FALSE))
  
  dfs <- df %>% 
    dplyr::select(c("pep_vmod", "pep_ivmod")) %>% 
    dplyr::mutate(pep_ivmod = gsub(" .*", "", pep_ivmod)) %>% 
    dplyr::filter(!grepl(",", pep_vmod), 
                  !is.na(pep_vmod)) %>% 
    unique() %>% 
    split(.$pep_vmod)
  
  are_terms <- grepl("[NC]-term", names(dfs))
  dfs <- dfs[!are_terms]
  
  if (!length(dfs))
    return(NULL)
  
  ids <- sort(unlist(lapply(dfs, hdeduce_mod_indexes)))
  
  # fn_prefix <- gsub("\\.[^.]*$", "", file)
  # saveRDS(ids, file.path(dat_dir, paste0("mod_indexes_", fn_prefix, ".rds")))
  
  invisible(ids)
}

#' Helper of \link{deduce_mod_indexes}.
#' 
#' @param df A subset of psmQ table with single modification.
hdeduce_mod_indexes <- function (df) 
{
  if (nrow(df)) {
    digits <- unlist(strsplit(df$pep_ivmod[1], ""))
    id <- unique(digits[digits != "0"])
    
    # should not occur
    if (length(id) != 1L) 
      stop("No matched `mod_index`.")
  }
  # should not occur
  else {
    stop("No matched `mod_index`.")
  }
  
  invisible(as.integer(id))
}


#' Find the path of MS1 data
#' 
#' Depreciated.
#' 
#' @param dat_dir A working directory.
#' @param pat A pattern for recognizing full MS1 files.
#' @param n_files The number of RAW files.
find_path_ms1 <- function (dat_dir = NULL, pat = "ms1full", n_files = 0L)
{
  if (is.null(dat_dir)) {
    dat_dir <- get_gl_dat_dir()
  }
  
  ms1files <- list.files(dat_dir, pattern = "^ms1full_.*\\.rds$")
  n_ms1fis <- length(ms1files)
  
  if (n_ms1fis) {
    path_ms1 <- dat_dir
  }
  else {
    path_ms1 <- file.path(dat_dir, "ms1data")
    ms1files <- list.files(path_ms1, pattern = "^ms1full_.*\\.rds$")
    n_ms1fis <- length(ms1files)
    
    if (!n_ms1fis) {
      # stop: otherwise df is unique by pep_seq_modz and need to be collapsed
      # to pep_seq_mod
      warning("MS1 peak lists of `^ms1full_[...].rds` not found for LFQ.\n", 
              "Please copy the entire Mzion folder of `ms1data` ", 
              "to your working directory.")
      return(NULL)
    }
  }
  
  if (n_ms1fis != n_files) {
    warning("The number of `^ms1full_[...].rds` files is different ", 
            "to the number of sample IDs. ", "LFQ disabled.")
    path_ms1 <- NULL
  }
  
  path_ms1
}



