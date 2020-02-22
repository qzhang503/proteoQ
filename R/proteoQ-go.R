#' Helper to create `db_path`
create_db_path <- function (db_path) {
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
}


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

  accessions <- purrr::quietly(AnnotationDbi::select)(
    get(pkg_nm), 
    keys = keys,
    columns = to,
    keytype = from
  )$result
}


#' Helper to get a complete list of `ENTREZID` from `egUNIPROT`
get_full_entrez <- function(species, from = "egUNIPROT") {
  abbr_species <- sp_lookup_Ul(species)
  
  pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
  pkg_nm_from <- paste("org", abbr_species, from, sep = ".")
  
  if (!requireNamespace(pkg_nm, quietly = TRUE)) {
    stop("Run `BiocManager::install(\"", pkg_nm, "\")` first.", call. = FALSE)
  }
  
  pkg_nm_from %>% get() %>% mappedkeys()
}


#' Helper to find human orthologs
#' #' @examples
#' res <- find_human_orthologs(mouse)
find_human_orthologs <- function(species, ortho_mart) {
  if (species == "human") stop("Ortholog `species` needs to be different to `human`.", call. = FALSE)
  
  out_nm <- paste0("ortho_hs", sp_lookup(species))

  data(package = "proteoQ", mart_hs)
  
  if (species == "mouse") {
    data(package = "proteoQ", mart_mm)
    martL <- mart_mm
  } else if (species == "rat") {
    data(package = "proteoQ", mart_rn)
    martL <- mart_rn
  } else {
    martL <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ortho_mart)
  }
  
  res <- biomaRt::getLDS(
    attributes = c("entrezgene_id"),
    filters = "entrezgene_id", 
    values = get_full_entrez(species = "human", from = "egUNIPROT"), 
    mart = mart_hs,
    attributesL = c("entrezgene_id"), 
    martL = martL
  ) %>% 
    `colnames<-`(c("human", species)) 
  
  assign(out_nm, res)
  save(list = out_nm, file = file.path(dat_dir, paste0(out_nm, ".rda")))

  invisible(res)
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


#' Helper to find the two-letter character string of abbreviated species
find_abbr_species <- function(species = "human", abbr_species = NULL) {
  species <- rlang::as_string(rlang::enexpr(species))
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
  
  return(abbr_species)
}


#' Helper to set the output file name for a data base
set_db_outname <- function(filename = NULL, species = "human", siganature) {
  filename <- rlang::enexpr(filename)
  
  if (is.null(filename)) {
    filename <- paste0(siganature, "_", tolower(sp_lookup_Ul(species)), ".rds")
  } else {
    filename <- rlang::as_string(filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    if (fn_prefix == fn_suffix) stop("No '.' in the file name.", call. = FALSE)
    if (fn_suffix != "rds") stop("File extension must be `.rds`.", call. = FALSE)
    filename <- paste0(fn_prefix, ".rds")
  }
}


#' Helper to download `gmt` 
dl_msig <- function(msig_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.all.v7.0.entrez.gmt", 
                    db_path = "~\\proteoQ\\dbs\\msig", overwrite = FALSE) {
  db_path <- create_db_path(db_path)
  
  if (!grepl("\\.entrez\\.gmt$", msig_url)) {
    stop("Use a link to a `.entrez.gmt` file; not `.symbols.gmt`.", call. = FALSE)
  }
  
  fn_msig <- msig_url %>% gsub("^.*/(.*)$", "\\1", .)
  
  if ((!file.exists(file.path(db_path, "cache", fn_msig))) | overwrite)  {
    downloader::download(msig_url, file.path(db_path, "cache", fn_msig), mode = "wb")
  }
  
  return(fn_msig)
}


#' Helper to save `gmt` 
proc_gmt <- function(species, abbr_species, ortho_mart, fn_gmt, db_path, filename) {
  filepath <- file.path(db_path, "cache", fn_gmt)
  if (!file.exists(filepath)) stop("File not found ", filepath, ".", call. = FALSE)
  
  df <- suppressWarnings(readr::read_tsv(filepath, col_names = FALSE)) %>% 
    `names_pos<-`(1, "term")
  
  df <- local({
    cols_entrez <- purrr::map_lgl(df, is.numeric) %>% which()
    
    df[, c(1, cols_entrez)] %>% 
      tidyr::gather("col", "entrez", -term) %>% 
      dplyr::select(-col)  %>% 
      dplyr::filter(!is.na(entrez)) 
  })
  
  if (species != "human") {
    df <- local({
      orthos <- find_human_orthologs(species, ortho_mart) %>% 
        `names<-`(c("from", "to"))
      
      df <- df %>% 
        dplyr::left_join(orthos, by = c("entrez" = "from")) %>% 
        dplyr::filter(!is.na(to))  %>% 
        dplyr::select(-entrez) %>% 
        dplyr::rename(entrez = to)      
    })
  }
  
  gsets <- df %>% 
    split(., .$term, drop = TRUE) %>% 
    purrr::map(`[[`, 2) %>% 
    `names<-`(paste(tolower(abbr_species), names(.), sep = "_"))
  
  saveRDS(gsets, file.path(db_path, filename))
}



#'Download and prepare gene ontology
#'
#'\code{prepGO} downloads and prepares data bases of
#'\code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'ontology}} (GO) for enrichment analysis by gene sets.
#'
#'@import rlang dplyr magrittr purrr fs readr downloader org.Hs.eg.db
#'  org.Mm.eg.db org.Rn.eg.db 
#'@inheritParams dl_stringdbs
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of GO. The species available for the
#'  convenience feature is in one of \code{c("human", "mouse", "rat")} with
#'  \code{"human"} being the default. The argument is not required for other
#'  species; instead, users will provide values under arguments
#'  \code{abbr_species}, \code{gaf_url} and \code{obo_url}.
#'@param abbr_species Two-letter character string; the abbreviated name of
#'  species used with
#'  \code{\href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Xx.eg.db}}.
#'   The value of \code{abbr_species} will be determined automatically if the
#'  species is in one of \code{c("human", "mouse", "rat")}. Otherwise, for
#'  example, users need to provide \code{abbr_species = Ce} for fetching the
#'  \code{org.Ce.eg.db} package in the name space of proteoQ. The argument is
#'  further applied to differentiate the same biological terms of
#'  \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'   ontology}},
#'  \code{\href{https://www.gsea-msigdb.org/gsea/index.jsp}{Molecular
#'  Signatures}} etc. under different species.
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

  db_path <- create_db_path(db_path)

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
  
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, "go")

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
    `names<-`(paste(tolower(abbr_species), names(.), sep = "_"))
  
  saveRDS(gsets, file.path(db_path, filename))
}


#'Preparation of Molecular Signatures
#'
#'\code{prepMSig} downloads and prepares data bases of
#'\code{\href{https://www.gsea-msigdb.org/gsea/index.jsp}{Molecular Signatures}}
#'(MSig) for enrichment analysis by gene sets.
#'
#'@import rlang dplyr magrittr purrr fs readr downloader biomaRt org.Hs.eg.db
#'@inheritParams dl_stringdbs
#'@inheritParams prepGO
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of \code{MSig} data bases. The species
#'  available for the convenience feature is in one of \code{c("human", "mouse",
#'  "rat")} with \code{"human"} being the default. The argument is not required
#'  for other species; instead, users will provide values under arguments
#'  \code{ortho_mart} for the lookup of orthologs to human.
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~\\proteoQ\\dbs\\msig"}.
#'@param msig_url A URL to
#'  \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{MSig
#'  }. At the \code{NULL} default, a \code{c2.all.v[...].entrez.gmt} data will
#'  be used for all species. A valid web address is required for a custom data
#'  base, and for simplicity, only files with entrez IDs will be handled.
#'@param ortho_mart Character string; a dataset name from
#'  \code{\link[biomaRt]{useMart}} and/or \code{\link[biomaRt]{listDatasets}}
#'  for the lookup of orthologs to \code{human} genes. For species in
#'  \code{c("human", "mouse", "rat")}, the value will be determined
#'  automatically unless otherwise specified.
#'@param filename Character string; An output file name. At the \code{NULL}
#'  default, the name will be determined automatically at a given
#'  \code{species}; i.e., \code{msig_hs.rds} for \code{human} data. The file is
#'  saved as a \code{.rds} object for uses with \code{\link{prnGSPA}}.
#' @examples
#' library(proteoQ)
#'
#' ## the default `MSig` is `c2.all`
#' # `human`; outputs under `db_path`
#' prepMSig()
#' head(readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_hs.rds")))
#'
#' # `mouse`
#' prepMSig(species = mouse, filename = msig_mm.rds)
#' head(readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_mm.rds")))
#'
#' # `rat`
#' prepMSig(species = rat, filename = msig_rn.rds)
#' head(readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_rn.rds")))
#'
#' # `dog`; need `ortho_mart` for species other than `human`, `mouse` and `rat`
#' # (try `?biomaRt::useMart` for a list of marts)
#' prepMSig(
#'   # species = dog,
#'   abbr_species = Cf, 
#'   ortho_mart = cfamiliaris_gene_ensembl,
#'   filename = msig_cf.rds,
#' )
#'
#' msig_cf <- readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_cf.rds"))
#'
#' # also `dog`
#' prepMSig(
#'   species = my_dog,
#'   abbr_species = Cf, 
#'   ortho_mart = cfamiliaris_gene_ensembl,
#'   filename = msig_cf2.rds,
#' )
#'
#' msig_cf2 <- readRDS(file.path("~\\proteoQ\\dbs\\msig\\msig_cf2.rds"))
#'
#'identical(msig_cf, msig_cf2)
#'
#' ## use an `MSig`other than the default of `c2.all`
#' prepMSig(
#'   msig_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.cgp.v7.0.entrez.gmt",
#'   species = human,
#'   filename = c2_cgp_hs.rds,
#' )
#'
#' prepMSig(
#'   msig_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.cgp.v7.0.entrez.gmt",
#'   species = dog,
#'   ortho_mart = cfamiliaris_gene_ensembl,
#'   filename = c2_cgp_cf.rds,
#' )
#'
#'@export
prepMSig <- function(msig_url = NULL, 
                     species = "human", 
                     abbr_species = NULL, 
                     ortho_mart = switch(species, 
                                         mouse = "mmusculus_gene_ensembl", 
                                         rat = "rnorvegicus_gene_ensembl", 
                                         human = "zzz",
                                         "unknown"), 
                     db_path = "~\\proteoQ\\dbs\\msig", filename = NULL, overwrite = FALSE) {

  old_opt <- options(warn = 0)
  options(warn = 1)
  on.exit(options(old_opt), add = TRUE)
  
  if (is.null(msig_url)) {
    msig_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.all.v7.0.entrez.gmt"
  } else {
    msig_url <- rlang::as_string(rlang::enexpr(msig_url))
  }
  
  fn_gmt <- dl_msig(msig_url, db_path, overwrite)
  
  species <- rlang::as_string(rlang::enexpr(species))
  
  ortho_mart <- local({
    ok <- tryCatch(force(ortho_mart), error = function(e) 1)
    
    if (ok == 1) {
      ortho_mart <- rlang::as_string(rlang::enexpr(ortho_mart))
    } else {
      ortho_mart <- ok
    }
  })

  if (ortho_mart == "unknown") stop("Specify the value of `ortho_mart` for species `", species, "`.", 
                                    call. = FALSE)
  
  if (species == "human" && ortho_mart != "zzz") species <- "unknown"
  
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, "msig")
  proc_gmt(species, abbr_species, ortho_mart, fn_gmt, db_path, filename)
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
  } else if (from == "egREFSEQ") {
    accessions <- accessions %>% dplyr::rename(refseq_acc = value)
  }

  write.table(accessions, file.path(db_dir, "temp", paste0(lwr_species, ".txt")), 
              sep = "\t", col.names = TRUE, row.names = FALSE)
  
  invisible(accessions)
}



