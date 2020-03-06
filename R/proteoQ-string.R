#'Downloads STRING databases
#'
#'@param species Character string; the species. The currently supported species
#'  include \code{human, mouse, rat, fly, cow, dog}. The default is
#'  \code{human}.
#'@param overwrite Logical; if TRUE, overwrite the downloaded database(s). The
#'  default is FALSE.
#'@inheritParams anal_prnString
#'@import rlang dplyr magrittr purrr fs downloader
#'@seealso \code{\link{anal_prnString}} for protein-protein interaction
#'  networks.
#'@export
dl_stringdbs <- function(species = "human", db_path = "~\\proteoQ\\dbs\\string", overwrite = FALSE) {
  species <- rlang::as_string(rlang::enexpr(species))

  if (!fs::dir_exists(db_path)) {
    new_db_path <- fs::path_expand_r(db_path)
    new_db_path2 <- fs::path_expand(db_path)
    
    if (fs::dir_exists(new_db_path)) {
      db_path <- new_db_path
    } else if (fs::dir_exists(new_db_path2)) {
      db_path <- new_db_path2
    } else {
      stop(db_path, " not existed.", call. = FALSE)
    }
    
    rm(new_db_path, new_db_path2)
  }

  urls <- switch(species, 
                 human = c(
                   "9606.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz",
                   "9606.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz", 
                   "9606.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz"
                 ), 
                 mouse = c(
                   "10090.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/10090.protein.links.full.v11.0.txt.gz",
                   "10090.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz",
                   "10090.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/10090.protein.info.v11.0.txt.gz"
                 ), 
                 rat = c(
                   "10116.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/10116.protein.links.full.v11.0.txt.gz",
                   "10116.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/10116.protein.aliases.v11.0.txt.gz", 
                   "10116.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/10116.protein.info.v11.0.txt.gz"
                 ), 
                 fly = c(
                   "7227.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/7227.protein.links.full.v11.0.txt.gz",
                   "7227.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/7227.protein.aliases.v11.0.txt.gz", 
                   "7227.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/7227.protein.info.v11.0.txt.gz"
                 ), 
                 cow = c(
                   "9913.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/9913.protein.links.full.v11.0.txt.gz", 
                   "9913.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/9913.protein.aliases.v11.0.txt.gz",
                   "9913.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/9913.protein.info.v11.0.txt.gz"
                 ), 
                 dog = c(
                   "9612.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/9612.protein.links.full.v11.0.txt.gz", 
                   "9612.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/9612.protein.aliases.v11.0.txt.gz", 
                   "9612.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/9612.protein.info.v11.0.txt.gz"
                 ), 
                 stop("Unknown `species`.", Call. = FALSE)
  )
  
  abbr_species <- sp_lookup(species) 
  taxid <- taxid_lookup(species)
  
  db_path2 <- file.path(db_path, abbr_species)
  dir.create(file.path(db_path2, "zip"), recursive = TRUE, showWarnings = FALSE)

  for(i in seq_along(urls)) {
    url <- urls[[i]]
    fn_zip <- names(urls)[[i]]
    fn_tsv <- gsub("\\.gz$", "", fn_zip)
    filepath <- file.path(db_path2, fn_zip)
    
    if ((!file.exists(filepath)) | overwrite) {
      downloader::download(url, filepath, mode = "wb")
      con <- gzfile(path.expand(filepath))
      
      if (grepl("protein\\.links", filepath)) {
        read.csv(con, sep = " ", check.names = FALSE, header = TRUE, comment.char = "#") %>% 
          write.table(file.path(db_path2, fn_tsv), sep = "\t", col.names = TRUE, row.names = FALSE)
      } else if (grepl("protein\\.info", filepath)) {
        temp <- readLines(con)
        for (idx in seq_along(temp)) {
          temp[idx] <- gsub("^(.*)\t[^\t].*$", "\\1", temp[idx])
        }
        temp[1] <- "protein_external_id\tpreferred_name\tprotein_size"
        writeLines(temp, file.path(db_path2, fn_tsv))
        
        rm(temp, idx)
      } else if (grepl("protein\\.alias", filepath)) {
        temp <- read.csv(con, sep = "\t", check.names = FALSE, header = TRUE, comment.char = "#")
        
        col_nms <- c("string_protein_id", "alias", "source")
        first_row <- names(temp) %>% 
          data.frame() %>% 
          t() %>% 
          `colnames<-`(col_nms)
        
        temp %>% 
          `colnames<-`(col_nms) %>% 
          dplyr::mutate_all(as.character) %>% 
          rbind(first_row, .) %>% 
          write.table(file.path(db_path2, fn_tsv), sep = "\t", col.names = TRUE, row.names = FALSE)
        
        rm(temp)
      }
      
      # close(con)
    }
  }

}


#'Annotates protein STRING ids by species
#'
#'@inheritParams info_anal
#'@inheritParams dl_stringdbs
#'@inheritParams anal_prnString
#'@import rlang dplyr purrr magrittr
annot_stringdb <- function(species, df, db_path, id, score_cutoff, 
                           filepath = NULL, filename = NULL, ...) {
  abbr_species <- sp_lookup(species) 
  taxid <- taxid_lookup(species)
  dl_stringdbs(!!species)
  
  db_path2 <- file.path(db_path, abbr_species)
  filelist_info <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.info", ".*.txt$"))
  filelist_link <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.links", ".*.txt$"))
  filelist_alias <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.aliases", ".*.txt$"))
  
  df <- df %>% dplyr::mutate(!!id := as.character(!!rlang::sym(id)))

  prn_info <- read.csv(file.path(db_path2, filelist_info), sep = "\t", check.names = FALSE, 
                       header = TRUE, comment.char = "#") %>% 
    dplyr::select(protein_external_id, preferred_name) %>% 
    dplyr::rename(!!id := preferred_name) %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id)))
  
  prn_links <- read.csv(file.path(db_path2, filelist_link), sep = "\t", 
                        check.names = FALSE, header = TRUE, comment.char = "#") %>% 
    dplyr::mutate(protein1 = as.character(protein1), protein2 = as.character(protein2))
  
  prn_alias <- read.csv(file.path(db_path2, filelist_alias), sep = "\t", 
                        check.names = FALSE, header = TRUE, comment.char = "#")
  
  prn_alias_sub <- prn_alias %>% 
    dplyr::filter(.$alias %in% df[[id]]) %>% 
    dplyr::filter(!duplicated(alias)) %>% 
    dplyr::select(-source) %>% 
    dplyr::rename(!!id := alias) %>% 
    dplyr::rename(protein_external_id = string_protein_id) %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id)))

  if (id == "gene") {
    string_map <- df %>%
      dplyr::select(id) %>% 
      dplyr::left_join(prn_info, by = id) %>% 
      dplyr::filter(!is.na(protein_external_id))
  } else {
    string_map <- df %>%
      dplyr::select(id) %>% 
      dplyr::left_join(prn_alias_sub, by = id) %>% 
      dplyr::filter(!is.na(protein_external_id))
  }
  string_map <- string_map %>% 
    dplyr::mutate(protein_external_id = as.character(protein_external_id))
  
  prn_links_sub <- prn_links %>% 
    dplyr::filter(protein1 %in% string_map$protein_external_id) %>% 
    dplyr::left_join(string_map, by = c("protein1" = "protein_external_id")) %>% 
    dplyr::rename(node1 = !!rlang::sym(id)) %>% 
    dplyr::left_join(string_map, by = c("protein2" = "protein_external_id")) %>% 
    dplyr::rename(node2 = !!rlang::sym(id))
  
  first_four <- c("node1", "node2", "protein1", "protein2")
  ppi <- dplyr::bind_cols(
    prn_links_sub[, first_four], 
    prn_links_sub[, !names(prn_links_sub) %in% first_four]
  ) %>% 
    dplyr::filter(!is.na(node1), !is.na(node2)) %>% 
    `names_pos<-`(1:4, c("#node1", "node2", "node1_external_id", "node2_external_id")) %>% 
    dplyr::filter(combined_score > score_cutoff)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  def_pre <- "^Protein_String_[NZ]"
  if (grepl(def_pre, fn_prefix)) {
    outnm_ppi <- paste0(fn_prefix, "_", species, "_ppi.tsv")
    outnm_expr <- paste0(fn_prefix, "_", species, "_expr.tsv")
  } else if (grepl(species, fn_prefix)) {
    outnm_ppi <- paste0(fn_prefix, "_ppi.tsv")
    outnm_expr <- paste0(fn_prefix, "_expr.tsv")
  } else {
    outnm_ppi <- paste0(fn_prefix, "_", species, "_ppi.tsv")
    outnm_expr <- paste0(fn_prefix, "_", species, "_expr.tsv")
  }
  
  write.table(ppi, file.path(dat_dir, "Protein\\String", outnm_ppi), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # expression data
  gns <- c(ppi[["#node1"]], ppi[["node2"]]) %>% unique()
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(id), 
    df %>% dplyr::select(grep("pVal|adjP|log2Ratio", names(.))),
    df %>% dplyr::select(-grep("pVal|adjP|log2Ratio", names(.)), -id)
  ) %>% 
    dplyr::filter(!!rlang::sym(id) %in% gns)

  suppressWarnings(df[is.na(df)] <- "") # Cytoscape compatibility

  write.table(df, file.path(dat_dir, "Protein\\String", outnm_expr), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
}


#' String analysis
#'
#' The input contains pVal fields
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnCorr_logFC
#' @inheritParams anal_prnString
#' @import dplyr purrr rlang fs
#' @importFrom magrittr %>%
stringTest <- function(df = NULL, id = gene, col_group = Group, col_order = Order, 
                       label_scheme_sub = NULL, 
                       db_path = "~\\proteoQ\\dbs\\string", score_cutoff = .7, 
                       scale_log2r = TRUE, complete_cases = FALSE, 
                       filepath = NULL, filename = NULL, ...) {
  
  # `scale_log2r` not used; both `_N` and `_Z` columns will be kept
  
  stopifnot(nrow(label_scheme_sub) > 0)
  stopifnot(rlang::is_double(score_cutoff))
  
  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots)
  
  if (score_cutoff <= 1) score_cutoff <- score_cutoff * 1000
  
  dir.create(file.path(dat_dir, "Protein\\String\\cache"), recursive = TRUE, showWarnings = FALSE)
  
  if (!fs::dir_exists(db_path)) {
    new_db_path <- fs::path_expand_r(db_path)
    new_db_path2 <- fs::path_expand(db_path)
    
    if (fs::dir_exists(new_db_path)) {
      db_path <- new_db_path
    } else if (fs::dir_exists(new_db_path2)) {
      db_path <- new_db_path2
    } else {
      stop(db_path, " not existed.", call. = FALSE)
    }
    
    rm(new_db_path, new_db_path2)
  }
  
  id <- rlang::as_string(rlang::enexpr(id))
  
  col_group <- rlang::enexpr(col_group)
  col_order <- rlang::enexpr(col_order)
  
  stopifnot(id %in% names(df))
  stopifnot("species" %in% names(df))
  
  species <- unique(df$species) %>% 
    .[!is.na(.)] %>% 
    as.character()
  
  stopifnot(length(species) >= 1)
  
  run_scripts <- FALSE
  if (run_scripts) {
    if (scale_log2r) {
      df <- df %>% 
        dplyr::select(-grep("N_log2_R", names(.)))
    } else {
      df <- df %>% 
        dplyr::select(-grep("Z_log2_R", names(.)))
    }    
  }

  df <- df %>% rm_pval_whitespace()

  purrr::walk(species, annot_stringdb, df, db_path, id, score_cutoff, filepath, filename)
}


#'STRING outputs of protein-protein interactions
#'
#'\code{anal_prnString} prepares the data of both
#'\href{https://string-db.org/}{STRING} protein-protein interactions
#'(ppi) and companion protein expressions.
#'
#'Convenience features, such as data row filtration via \code{filter_} varargs,
#'are available. The ppi file, \code{..._ppi.tsv}, and the expression file,
#'\code{..._expr.tsv}, are also compatible with a third-party
#'\href{https://cytoscape.org/}{Cytoscape}.
#'
#'@inheritParams anal_prnTrend
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~\\proteoQ\\dbs\\string"}.
#'@param score_cutoff Numeric; the threshold in the \code{combined_score} of
#'  protein-protein interaction. The default is 0.7.
#'@param scale_log2r Not currently used. Values before and after scaling will be
#'  both reported.
#'@param filename Use system default. Otherwise, the basename will be prepended
#'  to \code{_[species]_ppi.tsv} for network data and \code{_[species]_expr.tsv}
#'  for expression data.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration
#'  against data in a primary file of \code{\\Model\\Protein[_impNA]_pVals.txt}.
#'  See also \code{\link{normPSM}} for the format of \code{filter_} statements.
#'  \cr \cr \code{arrange_}: Variable argument statements for the row ordering
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{prnHM}} for the format of \code{arrange_} statements.
#'@seealso \code{\link{dl_stringdbs}} for database downloads. \cr
#'  \code{\link{prnSig}} for significance tests \cr
#'@example inst/extdata/examples/getStringDB_.R
#'
#'@import dplyr purrr rlang fs
#'@export
anal_prnString <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                            db_path = "~\\proteoQ\\dbs\\string", score_cutoff = .7, 
                            df = NULL, filepath = NULL, filename = NULL, ...) {
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_pepString"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), current_env()) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_prnString"))
    }
    , add = TRUE
  )
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, col_select = NULL, col_group = NULL, col_order = NULL,
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na,
            df = !!df, df2 = NULL, filepath = !!filepath, filename = !!filename,
            anal_type = "String")(db_path = db_path, score_cutoff = score_cutoff, adjP = adjP, 
                                  ...)
}




