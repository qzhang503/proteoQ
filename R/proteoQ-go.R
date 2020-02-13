#' Helper to save `obo` without header
proc_obo <- function(db_path, fn_obo, type = c("biological_process", "cellular_component", "molecular_function")) {
  filepath <- file.path(db_path, "cache", fn_obo)
  if (!file.exists(filepath)) stop("File not found ", filepath, ".", call. = FALSE)

  suppressWarnings(df <- readLines(filepath))
  first_row <- grep("\\[Term\\]", df)[1]
  last_row <- grep("\\[Typedef\\]", df)[1] - 2
  df <- df[first_row:last_row]
  
  go_ids <- df %>% .[grepl("^id:\\s{1}", .)] %>% gsub("^id:\\s{1}", "", .)
  go_nms <- df %>% .[grepl("^name:\\s{1}", .)] %>% gsub("^name:\\s{1}", "", .)
  go_type <- df %>% .[grepl("^namespace:\\s{1}", .)] %>% gsub("^namespace:\\s{1}", "", .)
  
  df <- tibble::tibble(go_id = go_ids, go_name = go_nms, go_space = go_type) %>% 
    dplyr::filter(go_space %in% type) %>% 
    dplyr::select(-go_space) 
}


#' Helper to save `gaf` without header
proc_gaf <- function(db_path, fn_gaf) {
  filepath <- file.path(db_path, "cache", fn_gaf)
  if (!file.exists(filepath)) stop("File not found ", filepath, ".", call. = FALSE)

  con <- gzfile(path.expand(filepath))
  
  suppressWarnings(df <- readLines(con))
  first_row <- grep("UniProtKB", df)[1]
  last_row <- length(df)
  df <- df[first_row:last_row]
  
  out_nm <- gsub("\\.gz$", "_hdr_rm.txt", fn_gaf)
  writeLines(df, file.path(db_path, out_nm))
  close(con)
  
  df <- readr::read_tsv(file.path(db_path, out_nm), col_names = FALSE, 
                        col_types = cols(
                          X1 = col_character(),
                          X2 = col_character(),
                          X3 = col_character(),
                          X4 = col_character(),
                          X5 = col_character(),
                          X6 = col_character(),
                          X7 = col_character(),
                          X8 = col_character(),
                          X9 = col_character(),
                          X10 = col_character(),
                          X11 = col_character(),
                          X12 = col_character(),
                          X13 = col_character(),
                          X14 = col_character(),
                          X15 = col_character(),
                          X16 = col_character(),
                          X17 = col_character()
                        )) 
  
  unlink(file.path(db_path, out_nm))
  
  df <- df %>% 
    dplyr::select(X2:X3, X5) %>% 
    `colnames<-`(c("uniprot_acc", "gene", "go_id"))
}


#' Helper to map `SYMBOL` to `ENTREZID` 
#' @import dplyr purrr AnnotationDbi
annot_from_to <- function(abbr_species = "Hs", keys = NULL, from = "SYMBOL", to = "ENTREZID") {
  if (all(is.null(keys))) stop("Argument `keys` cannot be NULL", call. = FALSE)

  pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
  
  if (!requireNamespace(pkg_nm, quietly = TRUE)) {
    stop("Run `BiocManager::install(\"", pkg_nm, "\")` first.", call. = FALSE)
  }

  x <- tryCatch(
    get(pkg_nm), 
    error = function(e) 1
  )
  
  if (!is.object(x)) {
    if (x == 1) stop("Did you forget to run `library(", pkg_nm, ")`?", call. = FALSE)
  }

  accessions <- AnnotationDbi::select(
    get(pkg_nm), 
    keys = keys,
    columns = to,
    keytype = from) 
}


#' Helper to look up GO species 
sp_lookup_go <- function(species) {
  switch(species, 
         human = "human",
         mgi = "mouse",
         rgd = "rat",
         fb = "fly", 
         cow = "cow",
         dog = "dog", 
         crap = "crap")
}


#'Download and prepare gene ontology
#'
#'\code{prepGO} prepares data bases of
#'\code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'ontology}} (GO).
#'
#'@import rlang dplyr magrittr purrr fs readr downloader org.Hs.eg.db
#'  org.Mm.eg.db org.Rn.eg.db org.Dm.eg.db org.Bt.eg.db org.Cf.eg.db
#'@inheritParams dl_stringdbs
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of GO. The species available for the
#'  convenience feature is in one of \code{c("human", "mouse", "rat")} with
#'  \code{"human"} being the default. The argument is not required for other
#'  species. Instead, users will provide values under arguments
#'  \code{abbr_species}, \code{gaf_url} and \code{obo_url}.
#'@param abbr_species Two-letter character string; the abbreviated name of
#'  species for uses with
#'  \code{\href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Xx.eg.db}}.
#'   If the species is in one of \code{c("human", "mouse", "rat")}, the value of
#'  \code{abbr_species} will be determined automatically. Otherwise, for
#'  example, users need to provide \code{abbr_species = Ce} for fetching the
#'  \code{org.Ce.eg.db} package in the name space of proteoQ.
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~\\proteoQ\\dbs\\go"}.
#'@param gaf_url A URL to
#'  \href{http://current.geneontology.org/products/pages/downloads.html}{GO
#'  Annotation File} (GAF). A valid web address is required for species other
#'  than \code{c("human", "mouse", "rat")}. At the NULL default and the species
#'  in one of \code{c("human", "mouse", "rat")}, the link will be determined
#'  automatically; note that users can overwrite the default GAF by providing
#'  their own URL.
#'@param obo_url A URL link to GO terms in an OBO format. At the NULL default,
#'  the web address will be determined automatically. Users can overwrite the
#'  default OBO by providing their own URL.
#'@param filename Character string; An output file name. At the NULL default,
#'  the name will be determined automatically at a given \code{species}; i.e.,
#'  \code{go_hs.rds} for human data. The file is saved as a \code{.rds} object
#'  for uses with \code{\link{prnGSPA}}.
#'@param type Character vector. The name space in gene ontology to be included.
#'  The default is to include all in \code{c("biological_process",
#'  "cellular_component", "molecular_function")}. In the example of \code{type =
#'  c("biological_process", "cellular_component")}, terms of
#'  \code{molecular_function} will be excluded.
#' @examples
#' library(proteoQ)
#'
#' # `human` and `mouse` with a default OBO; 
#' # outputs under `db_path`
#' prepGO(human)
#' prepGO(mouse)
#'
#' # `mouse` with a slim OBO
#' prepGO(
#'   species = mouse,
#'   obo_url = "http://current.geneontology.org/ontology/subsets/goslim_mouse.obo",
#'   filename = mm_slim.rds,
#' )
#'
#' # `worm` not available for default GO preparation
#' if (!requireNamespace("BiocManager", quietly = TRUE))
#'     install.packages("BiocManager")
#' BiocManager::install("org.Ce.eg.db")
#'
#' library(org.Ce.eg.db)
#'
#' prepGO(
#'   # species = worm,
#'   abbr_species = Ce,
#'   gaf_url = "http://current.geneontology.org/annotations/wb.gaf.gz",
#'   obo_url = "http://purl.obolibrary.org/obo/go/go-basic.obo",
#' )
#'
#' \dontrun{
#' gsets <- readRDS(file.path("~\\proteoQ\\dbs\\go", "mm_slim.rds"))
#' }
#'@export
prepGO <- function(species = "human", abbr_species = NULL, gaf_url = NULL, obo_url = NULL, 
                   db_path = "~\\proteoQ\\dbs\\go", 
                   type = c("biological_process", "cellular_component", "molecular_function"), 
                   filename = NULL, overwrite = FALSE) {
  
  old_opt <- options(warn = 0)
  options(warn = 1)
  on.exit(options(old_opt), add = TRUE)
  
  species <- rlang::as_string(rlang::enexpr(species))
  
  db_path <- local({
    if (!fs::dir_exists(db_path)) {
      new_db_path <- fs::path_expand_r(db_path)
      new_db_path2 <- fs::path_expand(db_path)
      
      if (fs::dir_exists(new_db_path)) {
        db_path <- new_db_path
      } else if (fs::dir_exists(new_db_path2)) {
        db_path <- new_db_path2
      } else {
        dir.create(file.path(db_path, "cache"), recursive = TRUE, showWarnings = FALSE)
        db_path <- new_db_path
      }
    }
  })
  
  if (is.null(gaf_url)) {
    gaf_url <- switch(species, 
                      human = c("goa_human.gaf.gz" = "http://current.geneontology.org/annotations/goa_human.gaf.gz"), 
                      mouse = c("goa_mouse.gaf.gz" = "http://current.geneontology.org/annotations/mgi.gaf.gz"), 
                      rat = c("goa_rat.gaf.gz" = "http://current.geneontology.org/annotations/rgd.gaf.gz"), 
                      stop("`species` need to be one of `human`, `mouse` or `rat` for an auto lookup of GAF files.", 
                           call. = FALSE))
    fn_gaf <- names(gaf_url)
  } else {
    fn_gaf <- gaf_url %>% gsub("^.*/(.*)$", "\\1", .)
    
    species <- local({
      species_go <- fn_gaf %>% 
        gsub("([^\\.]*)[\\.].*", "\\1", .) %>% 
        gsub("^goa_", "", .) %>% 
        sp_lookup_go()
      
      if (is.null(species_go)) {
        species_go <- "unknown"
      } else if (species_go != species) {
        message("The species is `", species_go, "`.")
      }
      
      species <- species_go      
    })
  }
  
  if (is.null(obo_url)) {
    obo_url <- c("go-basic.obo" = "http://purl.obolibrary.org/obo/go/go-basic.obo")
    fn_obo <- names(obo_url)
  } else {
    fn_obo <- obo_url %>% gsub("^.*/(.*)$", "\\1", .)
  }

  if ((!file.exists(file.path(db_path, "cache", fn_gaf))) | overwrite)  {
    downloader::download(gaf_url, file.path(db_path, "cache", fn_gaf), mode = "wb")
  }
  
  if ((!file.exists(file.path(db_path, "cache", fn_obo))) | overwrite) {
    downloader::download(obo_url, file.path(db_path, "cache", fn_obo), mode = "wb")
  }
  
  abbr_species <- rlang::enexpr(abbr_species)
  if (is.null(abbr_species)) {
    abbr_species <- sp_lookup_Ul(species)
    abbr_species_lwr <- sp_lookup(species) 
  } else {
    abbr_species <- rlang::as_string(abbr_species)
    abbr_species_lwr <- tolower(abbr_species)
  }

  if (abbr_species_lwr == "unknown") {
    stop("Provide the two-letter abbreviation of `abbr_species = Xx`; 
         i.e., `abbr_species = Ce` for uses with `org.Ce.eg.db` annotation.", 
         call. = FALSE)
  }
  
  filename <- rlang::enexpr(filename)
  if (is.null(filename)) {
    fn_prefix <- paste0("go_", abbr_species_lwr)
    filename <- paste0(fn_prefix, ".rds")
  } else {
    filename <- rlang::as_string(filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    if (fn_prefix == fn_suffix) stop("No '.' in the file name.", call. = FALSE)
    if (fn_suffix != "rds") warning("File extension must be `.rds`.", call. = FALSE)
    filename <- paste0(fn_prefix, ".rds")
  }
  
  df <- local({
    df_gaf <- proc_gaf(db_path, fn_gaf)
    df_obo <- proc_obo(db_path, fn_obo, type)
    
    dplyr::left_join(df_gaf, df_obo, by = "go_id") %>% 
      dplyr::mutate(go_name = paste(go_id, go_name)) %>% 
      dplyr::select(go_name, gene)
  })
  
  accessions <- annot_from_to(abbr_species = abbr_species, 
                              keys = unique(df$gene), from = "SYMBOL", to = "ENTREZID") %>% 
    dplyr::rename(gene = SYMBOL, entrez = ENTREZID) %>% 
    dplyr::filter(!is.na(entrez), !is.na(gene))  %>% 
    dplyr::filter(!duplicated(gene))
  
  gsets <- accessions %>% 
    dplyr::left_join(df, by = "gene") %>% 
    dplyr::select(-gene) %>% 
    split(., .$go_name, drop = TRUE) %>% 
    purrr::map(`[[`, 1) %>% 
    `names<-`(paste(abbr_species_lwr, names(.), sep = "_"))
  
  saveRDS(gsets, file.path(db_path, filename))
}














#' Map uniprot or refseq to entrez (not currently used)
#'
#' @import dplyr purrr tidyr plyr reshape2 org.Hs.eg.db org.Mm.eg.db
#'   org.Rn.eg.db
#' @examples
#' \dontrun{
#' # available species: c("human", "mouse", "rat")
#' res <- map_to_entrez(species = "human", from = "egUNIPROT")
#' dir.create(file.path("~\\proteoQ\\dbs\\entrez\\to_unirpot"), recursive = TRUE, showWarnings = FALSE)
#' write.table(res, file.path("~\\proteoQ\\dbs\\entrez\\to_unirpot", paste0("uniprot_entrez_", species, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
#'
#' res <- map_entrez(species, from = "egREFSEQ")
#' dir.create(file.path("~\\proteoQ\\dbs\\entrez\\to_refseq"), recursive = TRUE, showWarnings = FALSE)
#' write.table(res, file.path("~\\proteoQ\\dbs\\entrez\\to_refseq", paste0("refseq_entrez_", species, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
#' }
map_to_entrez <- function(species, from) {
  taxid <- taxid_lookup(species)
  abbr_species <- sp_lookup_Ul(species) 
  lwr_species <- tolower(abbr_species)
  
  pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")

  if (!requireNamespace(pkg_nm, quietly = TRUE)) {
    stop("Run `BiocManager::install(\"", pkg_nm, "\")` first.", call. = FALSE)
  }

  x <- get(paste("org", abbr_species, from, sep = "."))
  mapped_genes <- mappedkeys(x) 
  
  accessions <- as.list(x[mapped_genes]) %>% 
    plyr::ldply(., rbind) %>% 
    `names_pos<-`(., 1, c("entrez")) %>% 
    `names_pos<-`(., 2:ncol(.), paste(from, 1:(length(.)-1), sep = ".")) %>% 
    mutate_at(.vars = grep("^eg", names(.)), ~ as.character(.x)) %>% 
    reshape2::melt(id = "entrez") %>% 
    dplyr::filter(!is.na(entrez), !is.na(value)) %>% 
    dplyr::select(-c("variable"))

  if (from == "egUNIPROT") {
    accessions <- accessions %>% dplyr::rename(uniprot_acc = value)
    # path <- file.path(db_dir, "temp", "to_uniprot")
    # out_nm <- paste0("uniprot_entrez_", lwr_species, ".txt")
  } else if (from == "egREFSEQ") {
    accessions <- accessions %>% dplyr::rename(refseq_acc = value)
    # path <- file.path(db_dir, "temp", "to_refseq")
    # out_nm <- paste0("refseq_entrez_", lwr_species, ".txt")
  }

  # dir.create(file.path(path), recursive = TRUE, showWarnings = FALSE)
  write.table(accessions, file.path(db_dir, "temp", paste0(lwr_species, ".txt")), 
              sep = "\t", col.names = TRUE, row.names = FALSE)
  
  invisible(accessions)
}



