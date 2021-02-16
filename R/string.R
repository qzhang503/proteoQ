#'Annotates protein STRING ids
#'
#'@inheritParams info_anal
#'@inheritParams anal_prnString
#'@import dplyr purrr 
#'@importFrom magrittr %>% %T>% %$% %<>%
annot_stringdb <- function(df, db_nms, id, score_cutoff, filepath = NULL, filename = NULL, ...) {
  dat_dir <- get_gl_dat_dir()
  
  df <- df %>% dplyr::mutate(!!id := as.character(!!rlang::sym(id)))
  dbs <- load_stringdbs(db_nms)

  prn_info <- dbs$info %>% 
    dplyr::select(protein_external_id, preferred_name) %>% 
    dplyr::rename(!!id := preferred_name) %>% 
    dplyr::mutate(!!id := as.character(!!rlang::sym(id)))

  prn_aliases <- dbs$aliases %>% 
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
      dplyr::left_join(prn_aliases, by = id) %>% 
      dplyr::filter(!is.na(protein_external_id))
  }
  string_map <- string_map %>% 
    dplyr::mutate(protein_external_id = as.character(protein_external_id))
  
  # --- ppi data ---
  prn_links <- dbs$links %>% 
    dplyr::mutate(protein1 = as.character(protein1), protein2 = as.character(protein2)) %>% 
    dplyr::filter(protein1 %in% string_map$protein_external_id) %>% 
    dplyr::left_join(string_map, by = c("protein1" = "protein_external_id")) %>% 
    dplyr::rename(node1 = !!rlang::sym(id)) %>% 
    dplyr::left_join(string_map, by = c("protein2" = "protein_external_id")) %>% 
    dplyr::rename(node2 = !!rlang::sym(id))
  
  first_four <- c("node1", "node2", "protein1", "protein2")
  ppi <- dplyr::bind_cols(
    prn_links[, first_four], 
    prn_links[, !names(prn_links) %in% first_four]
  ) %>% 
    dplyr::filter(!is.na(node1), !is.na(node2)) %>% 
    `names_pos<-`(1:4, c("#node1", "node2", "node1_external_id", "node2_external_id")) %>% 
    dplyr::filter(combined_score > score_cutoff)
  
  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
  fn_prefix <- gsub("\\.[^.]*$", "", filename)
  outnm_ppi <- paste0(fn_prefix, "_ppi.tsv")
  outnm_expr <- paste0(fn_prefix, "_expr.tsv")

  write.table(ppi, file.path(dat_dir, "Protein/String", outnm_ppi), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # --- expression data ---
  gns <- c(ppi[["#node1"]], ppi[["node2"]]) %>% unique()
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(id), 
    df %>% dplyr::select(grep("pVal|adjP|log2Ratio", names(.))),
    df %>% dplyr::select(-grep("pVal|adjP|log2Ratio", names(.)), -id)
  ) %>% 
    dplyr::filter(!!rlang::sym(id) %in% gns)

  suppressWarnings(df[is.na(df)] <- "") # Cytoscape compatibility

  write.table(df, file.path(dat_dir, "Protein/String", outnm_expr), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
}


#' String analysis
#'
#' The input \code{df} contains pVal fields.
#' 
#' The argument \code{scale_log2r} is not used in that both `_N` and `_Z`
#' columns from primary \code{df} will be kept. The argument \code{species} is
#' used for the generation of separate outputs by \code{species}.
#' 
#' @inheritParams info_anal
#' @inheritParams gspaTest
#' @inheritParams prnCorr_logFC
#' @inheritParams anal_prnString
#' @import dplyr purrr fs
#' @importFrom magrittr %>% %T>% %$% %<>% 
stringTest <- function(df = NULL, id = gene, 
                       label_scheme_sub = NULL, 
                       db_nms = NULL, score_cutoff = .7, 
                       scale_log2r = TRUE, complete_cases = FALSE, 
                       filepath = NULL, filename = NULL, ...) {
  dat_dir <- get_gl_dat_dir()
  id <- rlang::as_string(rlang::enexpr(id))
  
  stopifnot(id %in% names(df), 
            nrow(label_scheme_sub) > 0, 
            vapply(c(score_cutoff), rlang::is_double, logical(1)))

  if (complete_cases) df <- df %>% my_complete_cases(scale_log2r, label_scheme_sub)
  
  dots <- rlang::enexprs(...)
  filter_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  arrange_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^arrange_", names(.))]
  dots <- dots %>% .[! . %in% c(filter_dots, arrange_dots)]

  df <- df %>% 
    filters_in_call(!!!filter_dots) %>% 
    arrangers_in_call(!!!arrange_dots) %>% 
    rm_pval_whitespace()
  
  if (score_cutoff <= 1) score_cutoff <- score_cutoff * 1000
  
  dir.create(file.path(dat_dir, "Protein/String/cache"), recursive = TRUE, showWarnings = FALSE)

  annot_stringdb(df, db_nms, id, score_cutoff, filepath, filename)
}


#'STRING outputs of protein-protein interactions
#'
#'\code{anal_prnString} prepares the data of both
#'\href{https://string-db.org/}{STRING} protein-protein interactions
#'(ppi) and companion protein expressions.
#'
#'The ppi file, \code{[...]_ppi.tsv}, and the expression file,
#'\code{[...]_expr.tsv}, are also compatible with
#'\href{https://cytoscape.org/}{Cytoscape}.
#'
#'@inheritParams anal_prnTrend
#'@param db_nms  Character string(s) to the name(s) of STRING database(s) with
#'  prepended directory path. The \code{STRING} database(s) need to match those
#'  generated from \code{\link{prepString}}. There is no default and users need
#'  to provide the correct file path(s) and name(s).
#'@param score_cutoff Numeric; the threshold in the \code{combined_score} of
#'  protein-protein interaction. The default is 0.7.
#'@param scale_log2r Not currently used. Values before and after scaling will be
#'  both reported.
#'@param filename Use system default. Otherwise, the user-provided basename will
#'  be prepended with \code{_ppi.tsv} for network data and \code{_expr.tsv} for
#'  expression data.
#'@param ... \code{filter_}: Variable argument statements for the row filtration
#'  against data in a primary file linked to \code{df}. See also
#'  \code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
#'  \code{arrange_}: Variable argument statements for the row ordering against
#'  data in a primary file linked to \code{df}. See also \code{\link{prnHM}} for
#'  the format of \code{arrange_} statements.
#'@example inst/extdata/examples/prnString_.R
#'@seealso \code{\link{prepString}} for database downloads and preparation. \cr
#'@import dplyr purrr fs
#'@export
anal_prnString <- function (scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                            db_nms = NULL, score_cutoff = .7, 
                            df = NULL, filepath = NULL, filename = NULL, ...) {
  on.exit(
    if (id %in% c("pep_seq", "pep_seq_mod")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_pepString"))
    } else if (id %in% c("prot_acc", "gene")) {
      mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
        c(enexprs(...)) %>% 
        save_call(paste0("anal", "_prnString"))
    }, 
    add = TRUE
  )
  
  check_dots(c("id", "anal_type", "df2"), ...)
  
  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"), length(id) == 1)
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)

  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, 
            col_select = NULL, 
            col_group = NULL, 
            col_order = NULL,
            scale_log2r = scale_log2r, 
            complete_cases = complete_cases, 
            impute_na = impute_na,
            df = !!df, 
            df2 = NULL, 
            filepath = !!filepath, 
            filename = !!filename,
            anal_type = "String")(db_nms = db_nms, 
                                  score_cutoff = score_cutoff, 
                                  ...)
}





#'Downloads and prepares STRING databases
#'
#'\code{prepString} downloads and prepares the data sets of
#'\href{https://string-db.org/}{STRING} interactions, accession lookups and gene
#'name lookups.
#'
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~/proteoQ/dbs/string"}.
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of STRING databases. The species available
#'  for the convenience feature is in one of \code{c("human", "mouse", "rat")}
#'  with \code{"human"} being the default. The argument is not required for
#'  other species; instead, users will provide values under arguments
#'  \code{links_url}, \code{aliases_url} and \code{info_url}.
#'@param links_url A URL to \href{http://string-db.org/}{STRING} interaction
#'  data. A valid web address is required for species other than
#'  \code{c("human", "mouse", "rat")}. At the NULL default and the species in
#'  one of \code{c("human", "mouse", "rat")}, the link will be determined
#'  automatically; note that users can overwrite the default interaction data by
#'  providing their own URL.
#'@param aliases_url A URL to \href{http://string-db.org/}{STRING} aliases
#'  between \code{string_protein_ids} and \code{protein accessions}. A valid web
#'  address is required for species other than \code{c("human", "mouse",
#'  "rat")}. At the NULL default and the species in one of \code{c("human",
#'  "mouse", "rat")}, the link will be determined automatically; note that users
#'  can overwrite the default alias data by providing their own URL.
#'@param info_url A URL to \href{http://string-db.org/}{STRING} display names
#'  between \code{string_protein_ids} and \code{gene names}. A valid web address
#'  is required for species other than \code{c("human", "mouse", "rat")}. At the
#'  NULL default and the species in one of \code{c("human", "mouse", "rat")},
#'  the link will be determined automatically; note that users can overwrite the
#'  default \code{info} data by providing their own URL.
#'@param overwrite Logical; if TRUE, overwrite the downloaded database(s). The
#'  default is FALSE.
#'@inheritParams anal_prnString
#'@examples
#'prepString(
#'  # species = fly, 
#'  links_url = "https://stringdb-static.org/download/protein.links.full.v11.0/7227.protein.links.full.v11.0.txt.gz",
#'  aliases_url = "https://stringdb-static.org/download/protein.aliases.v11.0/7227.protein.aliases.v11.0.txt.gz",
#'  info_url = "https://stringdb-static.org/download/protein.info.v11.0/7227.protein.info.v11.0.txt.gz", 
#'  filename = string_dm.rds,
#')
#'@import dplyr purrr fs
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@seealso \code{\link{anal_prnString}} for protein-protein interaction
#'  networks.
#'@export
prepString <- function(species = "human", # abbr_species = NULL, 
                       links_url = NULL, aliases_url = NULL, info_url = NULL, 
                       db_path = "~/proteoQ/dbs/string", filename = NULL, overwrite = FALSE) {
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  if (!requireNamespace("downloader", quietly = TRUE)) {
    stop("\n====================================================================", 
         "\nNeed package \"downloader\" for this function to work.",
         "\n====================================================================",
         call. = FALSE)
  }
  
  species <- rlang::as_string(rlang::enexpr(species))
  db_path <- create_db_path(db_path)

  if (is.null(links_url)) {
    links_url <- switch(species, 
                        human = c("9606.protein.links.full.v11.0.txt.gz" = 
                                    "https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz"), 
                        mouse = c("10090.protein.links.full.v11.0.txt.gz" = 
                                    "https://stringdb-static.org/download/protein.links.full.v11.0/10090.protein.links.full.v11.0.txt.gz"), 
                        rat = c("10116.protein.links.full.v11.0.txt.gz" = 
                                  "https://stringdb-static.org/download/protein.links.full.v11.0/10116.protein.links.full.v11.0.txt.gz"), 
                        stop("`species` need to be one of `human`, `mouse` or `rat` for an auto lookup of `links` files.", call. = FALSE))
    fn_links <- names(links_url)
    taxid_links <- fn_links %>% gsub("^([0-9]+)\\..*", "\\1", .)
  } else {
    fn_links <- links_url %>% gsub("^.*/(.*)$", "\\1", .)
    taxid_links <- fn_links %>% gsub("^([0-9]+)\\..*", "\\1", .)
    
    if (taxid_links != "9606" && species == "human") {
      species <- local({
        data(uniprot_species, package = "proteoQ")
        
        species_string <- uniprot_species %>% 
          dplyr::filter(.data$taxid == .env$taxid_links) %>% 
          dplyr::select(organism) %>% 
          unlist()
        
        if (is.null(species_string)) {
          species_string <- "unknown"
        } else if (species_string != species) {
          message("The species is `", species_string, "`.")
        }
        
        rm(uniprot_species, envir = .GlobalEnv)
        species <- species_string  
      })
    }
  }
  
  if (is.null(aliases_url)) {
    aliases_url <- switch(species, 
                          human = c("9606.protein.aliases.v11.0.txt.gz" = 
                                      "https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz"), 
                          mouse = c("10090.protein.aliases.v11.0.txt.gz" = 
                                      "https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz"), 
                          rat = c("10116.protein.aliases.v11.0.txt.gz" = 
                                    "https://stringdb-static.org/download/protein.aliases.v11.0/10116.protein.aliases.v11.0.txt.gz"), 
                          stop("`species` need to be one of `human`, `mouse` or `rat` for an auto lookup of `alias` files.", call. = FALSE))
    fn_aliases <- names(aliases_url)
  } else {
    fn_aliases <- aliases_url %>% gsub("^.*/(.*)$", "\\1", .)
  }
  taxid_aliases <- fn_aliases %>% gsub("^([0-9]+)\\..*", "\\1", .)
  
  if (is.null(info_url)) {
    info_url <- switch(species, 
                       human = c("9606.protein.info.v11.0.txt.gz" = 
                                   "https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz"), 
                       mouse = c("10090.protein.info.v11.0.txt.gz" = 
                                   "https://stringdb-static.org/download/protein.info.v11.0/10090.protein.info.v11.0.txt.gz"), 
                       rat = c("10116.protein.info.v11.0.txt.gz" = 
                                 "https://stringdb-static.org/download/protein.info.v11.0/10116.protein.info.v11.0.txt.gz"), 
                       stop("`species` need to be one of `human`, `mouse` or `rat` ", 
                            "for an auto lookup of `info` files.", 
                            call. = FALSE))
    fn_info <- names(info_url)
  } else {
    fn_info <- info_url %>% gsub("^.*/(.*)$", "\\1", .)
  }
  taxid_info <- fn_info %>% gsub("^([0-9]+)\\..*", "\\1", .)
  
  if (length(unique(c(taxid_links, taxid_aliases, taxid_info))) > 1) 
    stop("Different `taxid` detected among `links_url`, `aliases_url` and `info_url`", 
         call. = FALSE)
  
  if (length(unique(c(links_url, aliases_url, info_url))) != 3) 
    stop("Duplicated `urls` detected.", 
         call. = FALSE)
  
  if (!grepl("links", links_url)) 
    warning("The `", links_url, "` does not contain character string `links`.", 
            call. = FALSE)
  
  if (!grepl("aliases", aliases_url)) 
    warning("The `", aliases_url, "` does not contain character string `aliases`.", 
            call. = FALSE)
  
  if (!grepl("info", info_url)) 
    warning("The `", info_url, "` does not contain character string `info`.", 
            call. = FALSE)

  df_links <- local({
    filepath <- file.path(db_path, fn_links)
    if ((!file.exists(filepath)) || overwrite) {
      downloader::download(links_url, filepath, mode = "wb")
    }
    
    con <- gzfile(path.expand(filepath))
    read.csv(con, sep = " ", check.names = FALSE, header = TRUE, comment.char = "#")
  })
  
  df_aliases <- local({
    filepath <- file.path(db_path, fn_aliases)
    if ((!file.exists(filepath)) || overwrite) {
      downloader::download(aliases_url, filepath, mode = "wb")
    }
    
    con <- gzfile(path.expand(filepath))
    temp <- read.csv(con, sep = "\t", check.names = FALSE, header = TRUE, comment.char = "#")
    
    col_nms <- c("string_protein_id", "alias", "source")
    first_row <- names(temp) %>% 
      data.frame() %>% 
      t() %>% 
      `colnames<-`(col_nms)
    
    temp %>% 
      `colnames<-`(col_nms) %>% 
      dplyr::mutate_all(as.character) %>% 
      rbind(first_row, .) 
  })
  
  df_info <- local({
    filepath <- file.path(db_path, fn_info)
    filepath2 <- file.path(db_path, gsub("\\.gz$", "", fn_info))
    
    if ((!file.exists(filepath)) || overwrite) {
      downloader::download(info_url, filepath, mode = "wb")
    }
    
    con <- gzfile(path.expand(filepath))
    temp <- readLines(con)
    for (idx in seq_along(temp)) {
      temp[idx] <- gsub("^(.*)\t[^\t].*$", "\\1", temp[idx])
    }
    temp[1] <- "protein_external_id\tpreferred_name\tprotein_size"
    writeLines(temp, filepath2)
    try(close(con))
    
    temp <- read.csv(filepath2, sep = "\t", 
             check.names = FALSE, header = TRUE, comment.char = "#")
    
    unlink(filepath2)
    return(temp)    
  })
  
  filename <- set_db_outname(!!rlang::enexpr(filename), species, "string")
  
  saveRDS(list(links = df_links, aliases = df_aliases, info = df_info), 
          file.path(db_path, filename))
  
  invisible(file.path(db_path, filename))
}



#'Loads species-specific databases of STRING
#'
#'A function loads a set of precompiled gene sets of 
#'\href{http://current.geneontology.org/products/pages/downloads.html}{GO}
#'and
#'\href{http://string-db.org/}{String}.
#'@seealso \code{\link{load_expts}} for supported species.
#'
#' @examples
#' \donttest{
#' prepString(human)
#' prepString(mouse)
#' 
#' load_stringdbs(
#'   db_nms = c("~/proteoQ/dbs/string/string_hs.rds",
#'              "~/proteoQ/dbs/string/string_mm.rds")
#' )
#' }
#'
#'@inheritParams anal_prnString
#'@import dplyr purrr fs
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@seealso \code{\link{load_dbs}} for loading databases of \code{GO} and
#'  \code{MSig}.
load_stringdbs <- function (db_nms = NULL) {
  if (is.null(db_nms)) stop("`db_nms` cannot be NULL.", call. = FALSE)
  
  if (!all(grepl("\\.rds$", db_nms))) {
    stop("Custom gene set files indicated by `db_nms` must end with the `.rds` extension.", 
         call. = FALSE)
  }
  
  local({
    not_oks <- db_nms %>% .[!file.exists(db_nms)]
    if (!purrr::is_empty(not_oks)) {
      stop("File(s) not found: \n", purrr::reduce(not_oks, paste, sep = ", \n"), 
           call. = FALSE)
    }    
  })

  dbs <- purrr::map(db_nms, readRDS)

  local({
    lens <- purrr::map(dbs,length)
    not_oks <- which(lens != 3)
    
    if (!purrr::is_empty(not_oks)) {
      stop("File(s) not containing all three pieces of `links`, `alias` and `info`: \n", 
           purrr::reduce(db_nms[not_oks], paste, sep = ", \n"), 
           call. = FALSE)
    }
  })

  links <- suppressWarnings(purrr::map(dbs, `[[`, 1) %>% dplyr::bind_rows())
  aliases <- suppressWarnings(purrr::map(dbs, `[[`, 2) %>% dplyr::bind_rows())
  info <- suppressWarnings(purrr::map(dbs, `[[`, 3) %>% dplyr::bind_rows())  
  
  invisible(list(links = links, aliases = aliases, info = info))
} 



