#'Downloads STRING databases
#'
#'@param species Character string; the species. The currently supported species
#'  include \code{human, mouse, rat, fly, bovine, dog}. The default is
#'  \code{human}.
#'@param db_path Character string; the local path for database(s). The default is \code{"~\\proteoQ\\dbs\\string"}.
#'@param overwrite Logical; if TRUE, overwrite the databse(s). The default is FALSE.
#'@import rlang dplyr magrittr purrr fs downloader
#'@seealso \code{\link{getStringDB}} for protein-protein interaction networks.
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
                 bovine = c(
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
#'@inheritParams getStringDB
#'@import dplyr purrr
annot_stringdb <- function(species, df, db_path, id, score_cutoff) {
  abbr_species <- sp_lookup(species) 
  taxid <- taxid_lookup(species)
  dl_stringdbs(!!species)
  
  db_path2 <- file.path(db_path, abbr_species)
  filelist_info <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.info", ".*.txt$"))
  filelist_link <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.links", ".*.txt$"))
  filelist_alias <- list.files(path = db_path2, pattern = paste0(taxid, ".protein.aliases", ".*.txt$"))
  
  prn_info <- read.csv(file.path(db_path2, filelist_info), sep = "\t", check.names = FALSE, 
                       header = TRUE, comment.char = "#") %>% 
    dplyr::select(protein_external_id, preferred_name) %>% 
    dplyr::rename(!!id := preferred_name)
  
  prn_links <- read.csv(file.path(db_path2, filelist_link), sep = "\t", 
                        check.names = FALSE, header = TRUE, comment.char = "#")
  
  prn_alias <- read.csv(file.path(db_path2, filelist_alias), sep = "\t", 
                        check.names = FALSE, header = TRUE, comment.char = "#")
  
  prn_alias_sub <- prn_alias %>% 
    dplyr::filter(.$alias %in% df[[id]]) %>% 
    dplyr::filter(!duplicated(alias)) %>% 
    dplyr::select(-source) %>% 
    dplyr::rename(!!id := alias) %>% 
    dplyr::rename(protein_external_id = string_protein_id)
  
  if (id == "gene") {
    string_map <- df %>%
      dplyr::select(id) %>% 
      dplyr::left_join(prn_info) %>% 
      dplyr::filter(!is.na(protein_external_id))
  } else {
    string_map <- df %>%
      dplyr::select(id) %>% 
      dplyr::left_join(prn_alias_sub) %>% 
      dplyr::filter(!is.na(protein_external_id))
  }
  
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
  
  write.table(ppi, file.path(dat_dir, "Protein\\String", paste0("Protein_STRING_ppi_", species, "_.tsv")), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # expression data file
  gns <- c(ppi[["#node1"]], ppi[["node2"]]) %>% unique()
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(id), 
    df %>% dplyr::select(grep("pVal|log2Ratio", names(.)))
  ) %>% 
    dplyr::filter(!!rlang::sym(id) %in% gns)

  df[is.na(df)] <- "" # Cytoscape compatibility
  
  write.table(df, file.path(dat_dir, "Protein\\String", paste0("Protein_STRING_expr_", species, "_.tsv")), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
}


#'STRING outputs
#'
#'\code{getStringDB} prepares locally the
#'\code{\href{https://string-db.org/}{STRING}} results of protein-protein
#'interaction and protein expressions. The interaction file,
#'\code{Protein_STRING_ppi.tsv}, and the expression file,
#'\code{Protein_STRING_expr.tsv}, are compatible with
#'\code{\href{https://cytoscape.org/}{Cytoscape}}.
#'
#'@inheritParams dl_stringdbs
#'@param score_cutoff Numeric; the threshold in the \code{combined_score} of
#'  protein-protein interaction. The default is 0.7. 
#'@param adjP Logical. At the FALSE default, unadjusted pVals from
#'  \code{\link{prnSig}} will be used in PPI reports; otherwise, the
#'  Benjamini-Hochberg pVals will be used in the reports.
#'@param ... \code{filter_}: Logical expression(s) for the row filtration of
#'  data; also see \code{\link{normPSM}}.
#'@seealso \code{\link{dl_stringdbs}} for database downloads. \cr 
#' \code{\link{prnSig}} for significance tests \cr 
#'@example inst/extdata/examples/getStringDB_.R
#'
#'@inheritParams proteoVolcano
#'@import dplyr purrr fs
#'@export
getStringDB <- function(db_path = "~\\proteoQ\\dbs\\string", score_cutoff = .7, adjP = FALSE, ...) {
  stopifnot(rlang::is_logical(adjP))
  stopifnot(rlang::is_double(score_cutoff))

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
  
  fn_p <- file.path(dat_dir, "Protein\\Model", "Protein_pVals.txt")
  fn_raw <- file.path(dat_dir, "Protein", "Protein.txt")
  
  if (file.exists(fn_p)) {
    src_path <- fn_p
  } else if (file.exists(fn_raw)) {
    src_path <- fn_raw
  } else {
    stop("Protein results not found.", call. = FALSE)
  }
  
  df <- read.csv(src_path, check.names = FALSE, stringsAsFactors = FALSE, 
                 header = TRUE, sep = "\t", comment.char = "#")
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(id %in% names(df))
  stopifnot("species" %in% names(df))
  
  species <- unique(df$species) %>% 
    .[!is.na(.)] %>% 
    as.character()
  
  stopifnot(length(species) >= 1)

  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  select_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^select_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots, select_dots)]
  
  # cat("Available column keys for data filtration: \n")
  # cat(paste0(names(df), "\n"))

  df <- df %>% 
    dplyr::filter(!is.na(id), nchar(id) > 0) %>% 
    filters_in_call(!!!filter_dots)

  if(!adjP) {
    df <- df %>%
      dplyr::select(-contains("adjP"))
  } else {
    df <- df %>%
      dplyr::select(-contains("pVal")) %>%
      `names<-`(gsub("adjP", "pVal", names(.)))
  }
  
  df <- df %>% 
    dplyr::select(id, grep("log2Ratio|pVal", names(.)))
  
  # remove white space before NA in scientific notation
  df <- df %>% rm_pval_whitespace()

  # check db_path exists?

  purrr::walk(species, annot_stringdb, df, db_path, id, score_cutoff)
}

