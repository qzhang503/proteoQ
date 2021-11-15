#' Helper to create `db_path`
#' 
#' @inheritParams prepGO
create_db_path <- function (db_path) 
{
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
  
  invisible(db_path)
}


#' Helper to save `obo` without header
#' 
#' @param fn_obo filename according to \code{obo_url}
#' @inheritParams prepGO
proc_obo <- function(db_path, fn_obo, 
                     type = c("biological_process", "cellular_component", "molecular_function")) 
{
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
#' 
#' @param fn_gaf filename according to \code{gaf_url}
#' @inheritParams prepGO
proc_gaf <- function(db_path, fn_gaf) 
{
  filepath <- file.path(db_path, "cache", fn_gaf)
  if (!file.exists(filepath)) stop("File not found ", filepath, ".", call. = FALSE)

  con <- gzfile(path.expand(filepath))
  
  suppressWarnings(df <- readLines(con))
  first_row <- grep("UniProtKB", df)[1]
  last_row <- length(df)
  df <- df[first_row:last_row]
  
  out_nm <- gsub("\\.gz$", "_hdr_rm.txt", fn_gaf)
  writeLines(df, file.path(db_path, out_nm))
  try(close(con))
  
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
#' 
#' @param keys Identifier such as gene names.
#' @param from the type of \code{keys}
#' @param to the type of target IDs
#' @inheritParams prepGO
#' @import dplyr purrr 
annot_from_to <- function(abbr_species = "Hs", keys = NULL, from = "SYMBOL", to = "ENTREZID") 
{
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
#' 
#' @inheritParams prepGO
#' @inheritParams annot_from_to
get_full_entrez <- function(species, from = "egUNIPROT") 
{
  abbr_species <- sp_lookup_Ul(species)
  
  pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
  pkg_nm_from <- paste("org", abbr_species, from, sep = ".")
  
  if (!requireNamespace(pkg_nm, quietly = TRUE)) {
    stop("Run `BiocManager::install(\"", pkg_nm, "\")` first.", call. = FALSE)
  }
  
  pkg_nm_from %>% get() %>% AnnotationDbi::mappedkeys()
}


#' Helper to find human orthologs
#' 
#' @inheritParams prepMSig
#' @examples \donttest{res <- find_human_orthologs(mouse)}
find_human_orthologs <- function(species, ortho_mart) 
{
  if (species == "human") {
    stop("Ortholog `species` needs to be different to `human`.", call. = FALSE)
  }
  
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
  
  invisible(res)
}


#' Helper to look up GO species 
#' 
#' @inheritParams prepMSig
sp_lookup_go <- function(species) 
{
  switch(species, 
         human = "human",
         mgi = "mouse",
         rgd = "rat")
}


#' Helper to find the two-letter character string of abbreviated species
#' 
#' @import stringr
#' @inheritParams prepMSig
find_abbr_species <- function(species = "human", abbr_species = NULL) 
{
  species <- rlang::as_string(rlang::enexpr(species))
  abbr_species <- rlang::enexpr(abbr_species)
  
  if (is.null(abbr_species) || species %in% c("human", "mouse", "rat")) {
    abbr_species <- switch(species, 
           human = "Hs",
           mouse = "Mm",
           rat = "Rn",
           stop("`species` not in one of `human`, `mouse` or `rat`.", 
                "\nThus users need to provide a two-letter abbreviation of `abbr_species`, ", 
                "\ni.e., `abbr_species = Ce` for later uses with `org.Ce.eg.db` annotation.", 
                call. = FALSE)
    )
  } else {
    abbr_species <- rlang::as_string(abbr_species)
    
    if (stringr::str_length(abbr_species) != 2) {
      warning("The number of characters is typically `2` for `abbr_species`.", call. = FALSE)
    }
    
    if (abbr_species != stringr::str_to_title(abbr_species)) {
      warning("An `abbr_species` is typically in Title case, i.e., `Xx`.", call. = FALSE)
    }    
  }

  return(abbr_species)
}


#' Helper to set the output file name for a data base
#' 
#' @param signature A character string, i.e. "go" for uses in an output filename.
#' @inheritParams prepMSig
set_db_outname <- function(filename = NULL, species = "human", signature) 
{
  filename <- rlang::enexpr(filename)
  
  if (is.null(filename)) {
    abbr_species_lwr <- switch(species, 
                               human = "hs", 
                               mouse = "mm", 
                               rat = "rn", 
                               stop("`species` not in one of `human`, `mouse` or `rat`.", 
                                    "\nThus users need to provide a `filename`.", 
                                    call. = FALSE))
    
    filename <- paste0(signature, "_", abbr_species_lwr, ".rds")
  } else {
    filename <- rlang::as_string(filename)
    fn_prefix <- gsub("\\.[^.]*$", "", filename)
    fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filename)
    if (fn_prefix == fn_suffix) stop("No '.' to separate a basename and an extension.", call. = FALSE)
    if (fn_suffix != "rds") stop("File extension must be `.rds`.", call. = FALSE)
    filename <- paste0(fn_prefix, ".rds")
  }
}


#' Helper to download `gmt`
#' 
#' @inheritParams prepMSig
dl_msig <- function(msig_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.all.v7.0.entrez.gmt", 
                    db_path = "~/proteoQ/dbs/msig", overwrite = FALSE) 
{
  if (!grepl("\\.entrez\\.gmt$", msig_url)) {
    stop("Use a link to a `.entrez.gmt` file; not `.symbols.gmt`.", call. = FALSE)
  }
  
  if (!requireNamespace("downloader", quietly = TRUE)) {
    stop("\n====================================================================", 
         "\nNeed package \"downloader\" for this function to work.",
         "\n====================================================================",
         call. = FALSE)
  }

  fn_msig <- msig_url %>% gsub("^.*/(.*)$", "\\1", .)
  
  if ((!file.exists(file.path(db_path, "cache", fn_msig))) | overwrite)  {
    downloader::download(msig_url, file.path(db_path, "cache", fn_msig), mode = "wb")
  }
  
  return(fn_msig)
}


#' Helper to save `gmt` 
#' 
#' @param fn_gmt filename of downloaded gmt results.
#' @inheritParams prepMSig
proc_gmt <- function(species, abbr_species, ortho_mart, fn_gmt, db_path, filename) 
{
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
  
  invisible(file.path(db_path, filename))
}


#'Download and prepare gene ontology
#'
#'\code{prepGO} downloads and prepares data bases of
#'\href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'ontology} (GO) for enrichment analysis by gene sets.
#'
#'@import dplyr purrr fs readr org.Hs.eg.db org.Mm.eg.db
#'  org.Rn.eg.db
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@param overwrite Logical; if TRUE, overwrite the downloaded database(s). The
#'  default is FALSE.
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of GO. The species available for the
#'  convenience feature is in one of \code{c("human", "mouse", "rat")} with
#'  \code{"human"} being the default. The argument is not required for other
#'  species; instead, users will provide values under arguments
#'  \code{abbr_species}, \code{gaf_url} and \code{obo_url}.
#'@param abbr_species Two-letter character string; the abbreviated name of
#'  species used with
#'  \href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Xx.eg.db}.
#'  The value of \code{abbr_species} will be determined automatically if the
#'  species is in one of \code{c("human", "mouse", "rat")}. Otherwise, for
#'  example, users need to provide \code{abbr_species = Ce} for fetching the
#'  \code{org.Ce.eg.db} package in the name space of proteoQ. 
#'  
#'  For analysis against
#'  \href{http://current.geneontology.org/products/pages/downloads.html}{gene
#'  ontology} and \href{https://www.gsea-msigdb.org/gsea/index.jsp}{Molecular
#'  Signatures}, the argument is further applied to differentiate the same
#'  biological terms under different species; e.g., \code{GO~0072686 mitotic
#'  spindle} becomes \code{hs_GO~0072686 mitotic spindle} for human and
#'  \code{mm_GO~0072686 mitotic spindle} for mouse.
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~/proteoQ/dbs/go"}.
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
#' \donttest{
#' library(proteoQ)
#'
#' # `human` and `mouse` with a default OBO;
#' # outputs under `db_path`
#' prepGO(human)
#' prepGO(mouse)
#' 
#' # head(readRDS(file.path("~/proteoQ/dbs/go/go_hs.rds")))
#' # head(readRDS(file.path("~/proteoQ/dbs/go/go_mm.rds")))
#' 
#' # enrichment analysis with custom `GO`
#' prnGSPA(
#'   gset_nms = c("~/proteoQ/dbs/go/go_hs.rds",
#'                "~/proteoQ/dbs/go/go_mm.rds"),
#' )
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
#'   filename = go_ce.rds,
#' )
#' }
#' 
#'@export
prepGO <- function(species = "human", abbr_species = NULL, gaf_url = NULL, obo_url = NULL, 
                   db_path = "~/proteoQ/dbs/go", 
                   type = c("biological_process", "cellular_component", "molecular_function"), 
                   filename = NULL, overwrite = FALSE) 
{
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  if (!requireNamespace("downloader", quietly = TRUE)) {
    stop("\n====================================================================", 
         "\nNeed package \"downloader\" for this function to work.",
         "\n====================================================================",
         call. = FALSE)
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("\n====================================================================", 
         "\nNeed package \"AnnotationDbi\" for this function to work.\n",
         "\nif (!requireNamespace(\"BiocManager\", quietly = TRUE)) ", 
         "\n\tinstall.packages(\"BiocManager\")",
         "\nBiocManager::install(\"AnnotationDbi\")", 
         "\n====================================================================",
         call. = FALSE)
  }
  

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

  if ((!file.exists(file.path(db_path, "cache", fn_gaf))) || overwrite)  {
    downloader::download(gaf_url, file.path(db_path, "cache", fn_gaf), mode = "wb")
    # download.file(gaf_url, file.path(db_path, "cache", fn_gaf), mode = "w")
  }
  
  if ((!file.exists(file.path(db_path, "cache", fn_obo))) | overwrite) {
    downloader::download(obo_url, file.path(db_path, "cache", fn_obo), mode = "wb")
    # download.file(obo_url, file.path(db_path, "cache", fn_obo), mode = "w")
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
  
  invisible(file.path(db_path, filename))
}


#'Preparation of Molecular Signatures
#'
#'\code{prepMSig} downloads and prepares data bases of
#'\href{https://www.gsea-msigdb.org/gsea/index.jsp}{Molecular Signatures}
#'(MSig) for enrichment analysis by gene sets.
#'
#'@import dplyr purrr fs readr org.Hs.eg.db
#'@importFrom magrittr %>% %T>% %$% %<>% 
#'@inheritParams prepGO
#'@param species Character string; the name of a species for the
#'  \emph{conveninent} preparation of \code{MSig} data bases. The species
#'  available for the convenience feature is in one of \code{c("human", "mouse",
#'  "rat")} with \code{"human"} being the default. The argument is not required
#'  for other species; instead, users will provide values under arguments
#'  \code{ortho_mart} for the lookup of orthologs to human.
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~/proteoQ/dbs/msig"}.
#'@param msig_url A URL to
#'  \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{MSig
#'  }. At the \code{NULL} default, a \code{c2.all.v[...].entrez.gmt} data will
#'  be used for all species. A valid web address is required for a custom data
#'  base. For simplicity, only files with entrez IDs will be handled; 
#'  files of \code{c2.all.v[...].symbols.gmt} will not be parsed. 
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
#' \donttest{
#' library(proteoQ)
#'
#' ## the default `MSig` is `c2.all`
#' # `human`; outputs under `db_path`
#' prepMSig()
#' head(readRDS(file.path("~/proteoQ/dbs/msig/msig_hs.rds")))
#' 
#' prnGSPA(
#'   gset_nms = file.path("~/proteoQ/dbs/msig/msig_hs.rds"), 
#' )
#'
#' # `mouse`
#' prepMSig(species = mouse, filename = msig_mm.rds)
#' head(readRDS(file.path("~/proteoQ/dbs/msig/msig_mm.rds")))
#'
#' # `rat`
#' prepMSig(species = rat, filename = msig_rn.rds)
#' head(readRDS(file.path("~/proteoQ/dbs/msig/msig_rn.rds")))
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
#' # also `dog`
#' prepMSig(
#'   species = my_dog,
#'   abbr_species = Cf, 
#'   ortho_mart = cfamiliaris_gene_ensembl,
#'   filename = msig_cf2.rds,
#' )
#'
#' msig_cf <- readRDS(file.path("~/proteoQ/dbs/msig/msig_cf.rds"))
#' msig_cf2 <- readRDS(file.path("~/proteoQ/dbs/msig/msig_cf2.rds"))
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
#' }
#' 
#' \dontrun{
#' # enrichment analysis with custom `MSig`
#' prnGSPA(
#'   gset_nms = c("~/proteoQ/dbs/msig/msig_hs.rds",
#'                "~/proteoQ/dbs/msig/msig_mm.rds"),
#' )
#' }
#'
#'@export
prepMSig <- function(species = "human", msig_url = NULL, abbr_species = NULL, 
                     ortho_mart = switch(species, 
                                         mouse = "mmusculus_gene_ensembl", 
                                         rat = "rnorvegicus_gene_ensembl", 
                                         human = "to_itself",
                                         "unknown"), 
                     db_path = "~/proteoQ/dbs/msig", filename = NULL, overwrite = FALSE) 
{
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("\n====================================================================", 
         "\nNeed package \"biomaRt\" for this function to work.\n",
         "\nif (!requireNamespace(\"BiocManager\", quietly = TRUE)) ", 
         "\n\tinstall.packages(\"BiocManager\")",
         "\nBiocManager::install(\"biomaRt\")", 
         "\n====================================================================",
         call. = FALSE)
  }

  if (is.null(msig_url)) {
    msig_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/c2.all.v7.0.entrez.gmt"
  } else {
    msig_url <- rlang::as_string(rlang::enexpr(msig_url))
  }
  
  db_path <- create_db_path(db_path)
  
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
  
  if (species == "human" && ortho_mart != "to_itself") species <- "unknown"

  if (ortho_mart == "unknown") 
    stop(
      "Specify the value of `ortho_mart` for species other than \"human\", \"mouse\", and \"rat\".", 
      call. = FALSE
    )

  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, "msig")
  out_path <- proc_gmt(species, abbr_species, ortho_mart, fn_gmt, db_path, filename)
  invisible(out_path)
}


#' Map UniProt or Refseq accessions to Entrez IDs
#'
#' @param from Character string; the type of accession keys in c("UNIPROT",
#'   "REFSEQ", "ACCNUM").
#' @inheritParams prepMSig
#' @inheritParams annot_from_to
#' @import dplyr purrr tidyr reshape2 org.Hs.eg.db org.Mm.eg.db
#'   org.Rn.eg.db
#' @importFrom plyr ldply
#' @export
map_to_entrez <- function(species = "human", abbr_species = NULL, from = "UNIPROT", 
                          filename = NULL, db_path = "~/proteoQ/dbs/entrez", 
                          overwrite = FALSE) 
{
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  db_path <- create_db_path(db_path)
  species <- rlang::as_string(rlang::enexpr(species))
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  from <- rlang::as_string(rlang::enexpr(from))
  filename <- set_db_outname(!!rlang::enexpr(filename), 
                             species, paste(tolower(from), 
                                            "entrez", 
                                            sep = "_" ))
  
  if ((!file.exists(file.path(db_path, filename))) || overwrite)  {
    pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
    
    if (!requireNamespace(pkg_nm, quietly = TRUE)) {
      stop("Run `BiocManager::install(\"", pkg_nm, "\")` first, 
           then `library(", pkg_nm, ")`", call. = FALSE)
    }
    
    new_from <- paste0("eg", from)
    x <- tryCatch(
      get(paste("org", abbr_species, new_from, sep = ".")),
      error = function(e) 1
    )
    
    if (!is.object(x)) {
      if (x == 1) stop("Did you forget to run `library(", pkg_nm, ")`?", 
                       call. = FALSE)
    }
    
    entrez_ids <- AnnotationDbi::mappedkeys(x) 
    
    accessions <- as.list(x[entrez_ids]) %>% 
      plyr::ldply(., rbind) %>% 
      `names_pos<-`(., 1, c("entrez")) %>% 
      `names_pos<-`(., 2:ncol(.), paste(new_from, 1:(length(.)-1), sep = ".")) %>% 
      dplyr::mutate_at(.vars = grep("^eg", names(.)), ~ as.character(.x)) %>% 
      tidyr::gather("variable", "value", -entrez) %>% 
      dplyr::filter(!is.na(entrez), !is.na(value)) %>% 
      dplyr::select(-c("variable")) %>% 
      dplyr::mutate(species = species)
    
    if (from == "UNIPROT") {
      accessions <- accessions %>% 
        dplyr::rename(uniprot_acc = value)
      
      accessions <- annot_from_to(abbr_species, unique(accessions$uniprot_acc), 
                             "UNIPROT", "SYMBOL") %>% 
        dplyr::rename(uniprot_acc = UNIPROT, gene = SYMBOL) %>% 
        dplyr::filter(!is.na(uniprot_acc), !is.na(gene))  %>% 
        dplyr::left_join(accessions, by = "uniprot_acc") %>% 
        dplyr::mutate(entrez = as.numeric(entrez)) %>% 
        dplyr::select(c("uniprot_acc", "gene", "entrez", "species"))
    } else if (from == "REFSEQ") {
      accessions <- accessions %>% dplyr::rename(refseq_acc = value)
      
      accessions <- annot_from_to(abbr_species, unique(accessions$refseq_acc), 
                                  "REFSEQ", "SYMBOL") %>% 
        dplyr::rename(refseq_acc = REFSEQ, gene = SYMBOL) %>% 
        dplyr::filter(!is.na(refseq_acc), !is.na(gene))  %>% 
        dplyr::left_join(accessions, by = "refseq_acc") %>% 
        dplyr::mutate(entrez = as.numeric(entrez)) %>% 
        dplyr::select(c("refseq_acc", "gene", "entrez", "species"))
    } else {
      stop("Variable `from` needs to be either `UNIPROT` or `REFSEQ`.", call. = FALSE)
    }
    
    saveRDS(accessions, file.path(db_path, filename))

    invisible(accessions)
  } else {
    invisible(NULL)
  }
}


#'Map UniProt accessions to Entrez IDs
#'
#'\code{Uni2Entrez} prepares lookup tables between UniProt accessions and
#'Entrez IDs for uses with \link{normPSM} and downstream gene-set analysis such
#'as \link{prnGSPA}. The utility is optional for \code{human}, \code{mouse} and
#'\code{rat} data. It is \strong{required} for other species with \link{prnGSPA}
#'in users' workflows. It can also be used to update and overrule the lookups
#'for \code{human}, \code{mouse} and \code{rat} that are defaulted by
#'\code{proteoQ}.
#'
#'@param species Character string; the name of a species. 
#'@param db_path Character string; the local path for database(s). The default
#'  is \code{"~/proteoQ/dbs/entrez"}.
#'@inheritParams prepMSig
#'@import dplyr purrr tidyr reshape2 org.Hs.eg.db org.Mm.eg.db org.Rn.eg.db
#'@example inst/extdata/examples/prepEntrez_.R
#'@export
Uni2Entrez <- function(species = "human", abbr_species = NULL, filename = NULL, 
                       db_path = "~/proteoQ/dbs/entrez", overwrite = FALSE) 
{
  map_to_entrez(!!rlang::enexpr(species), 
                !!rlang::enexpr(abbr_species), 
                "UNIPROT", 
                !!rlang::enexpr(filename), 
                db_path, 
                overwrite)
}


#'Map RefSeq accessions to Entrez IDs and gene names
#'
#'\code{Ref2Entrez} prepares lookup tables between RefSeq accessions and
#'Entrez IDs and gene names for uses with \link{normPSM} and downstream gene-set
#'analysis such as \link{prnGSPA}. The utility is optional for \code{human} and
#'\code{mouse} data. It is \strong{required} for other species with
#'\link{prnGSPA} in users' workflows. It can also be used to update and overrule
#'the lookups for \code{human} and \code{mouse} that are defaulted by
#'\code{proteoQ}.
#'
#'@rdname Uni2Entrez
#'@import dplyr purrr tidyr reshape2 org.Hs.eg.db org.Mm.eg.db org.Rn.eg.db
#'@example inst/extdata/examples/prepEntrez_.R
#'@export
Ref2Entrez <- function(species = "human", abbr_species = NULL, filename = NULL, 
                       db_path = "~/proteoQ/dbs/entrez", overwrite = FALSE) 
{
  map_to_entrez(!!rlang::enexpr(species), 
                !!rlang::enexpr(abbr_species), 
                "REFSEQ", 
                !!rlang::enexpr(filename), 
                db_path, 
                overwrite)
}


#' Map uniprot or refseq to entrez (not currently used)
#'
#' @param os_name An organism name by UniProt.
#' @inheritParams Uni2Entrez
#' @inheritParams annot_from_to
#' @import dplyr purrr tidyr reshape2 org.Hs.eg.db org.Mm.eg.db org.Rn.eg.db
#' @importFrom plyr ldply
map_to_entrez_os_name <- function(species = "human", abbr_species = NULL, 
                                  os_name = "Homo sapiens", 
                                  from = "UNIPROT", 
                                  filename = NULL, db_path = "~/proteoQ/dbs/entrez", 
                                  overwrite = FALSE) 
{
  # the value of species will be used under column `species` in Protein.txt etc
  # add column species in the rds output
  
  db_path <- create_db_path(db_path)
  species <- rlang::as_string(rlang::enexpr(species))
  os_name <- rlang::as_string(rlang::enexpr(os_name))
  
  create_os_lookup(species, os_name, overwrite)
  
  abbr_species <- find_abbr_species(!!species, !!rlang::enexpr(abbr_species))
  filename <- set_db_outname(!!rlang::enexpr(filename), species, paste(tolower(from), "entrez", sep = "_" ))
  
  if ((!file.exists(file.path(db_path, filename))) || overwrite)  {
    pkg_nm <- paste("org", abbr_species, "eg.db", sep = ".")
    if (!requireNamespace(pkg_nm, quietly = TRUE)) {
      stop("Run `BiocManager::install(\"", pkg_nm, "\")` first, 
           then `library(", pkg_nm, ")`", call. = FALSE)
    }
    
    new_from <- paste0("eg", from)
    x <- tryCatch(
      get(paste("org", abbr_species, new_from, sep = ".")),
      error = function(e) 1
    )
    
    if (!is.object(x)) {
      if (x == 1) stop("Did you forget to run `library(", pkg_nm, ")`?", call. = FALSE)
    }
    
    entrez_ids <- AnnotationDbi::mappedkeys(x) 
    
    accessions <- as.list(x[entrez_ids]) %>% 
      plyr::ldply(., rbind) %>% 
      `names_pos<-`(., 1, c("entrez")) %>% 
      `names_pos<-`(., 2:ncol(.), paste(new_from, 1:(length(.)-1), sep = ".")) %>% 
      dplyr::mutate_at(.vars = grep("^eg", names(.)), ~ as.character(.x)) %>% 
      tidyr::gather("variable", "value", -entrez) %>% 
      dplyr::filter(!is.na(entrez), !is.na(value)) %>% 
      dplyr::select(-c("variable")) # %>% 
    # dplyr::mutate(species = species)
    
    if (from == "UNIPROT") {
      accessions <- accessions %>% dplyr::rename(uniprot_acc = value)
    } else if (from == "REFSEQ") {
      accessions <- accessions %>% dplyr::rename(refseq_acc = value)
    } else {
      stop("Variable `from` needs to be either `UNIPROT` or `REFSEQ`.", call. = FALSE)
    }
    
    saveRDS(accessions, file.path(db_path, filename))
    
    invisible(accessions)
  } else {
    invisible(NULL)
  }
}



#' create a lookup table and remove duplicated entries
#' 
#' @inheritParams map_to_entrez_os_name
create_os_lookup <- function(species, os_name, overwrite = FALSE) 
{
  my_lookup <- c(
    "Homo sapiens" = "human",
    "Mus musculus" = "mouse",
    "Rattus norvegicus" = "rat"
  )

  # check if conflict between species and os_name
  convert_default_os <- function (species, os_name) {
    if (species %in% my_lookup) {
      lookup_nm <- my_lookup %>% .[. %in% species] %>% names()
      if (!purrr::is_empty(lookup_nm) && lookup_nm != os_name) {
        os_name <- lookup_nm
      }
    }
    
    return(os_name)
  }
  
  dat_dir <- get_gl_dat_dir()
  
  os_name <- convert_default_os(species, os_name)
  curr_lookup <- uniprot_entrez_lookup <- setNames(species, os_name)

  if (os_name == "Homo sapiens" && species != "human") 
    stop("`os_name = Homo sapiens` is reversed for `species = human`.", call. = FALSE)
  
  if (os_name == "Mus musculus" && species != "mouse") 
    stop("`os_name = Mus musculus` is reversed for `species = mouse`.", call. = FALSE)
  
  if (os_name == "Rattus norvegicus" && species != "rat") 
    stop("`os_name = Rattus norvegicus` is reversed for `species = rat`.", call. = FALSE)  
  
  # check if the same os_name but diferent species
  file <- file.path(dat_dir, "uniprot_entrez_lookup.rda")
  
  if (file.exists(file)) {
    load(file)
    
    if (overwrite) {
      uniprot_entrez_lookup <- uniprot_entrez_lookup %>% 
        .[! names(.) == os_name] %>% 
        .[! . == species]
    } else {
      old_sp <- uniprot_entrez_lookup %>% .[names(.) == os_name]
      curr_sp <- curr_lookup

      if (!purrr::is_empty(old_sp) && (curr_sp != old_sp)) {
        stop("`", names(curr_lookup), "` was previously linked to `species = ", old_sp, "`.", 
             "\nTo overwrite the value of `os_name`, set `overwrite = TRUE`.", 
             call. = FALSE)
      }
      
      old_nm <- uniprot_entrez_lookup %>% .[. == species] %>% names()
      curr_nm <- os_name

      if (!purrr::is_empty(old_nm) && (curr_nm != old_nm)  && !overwrite) {
        stop("`", species, "` was previously linked to `os_name = ", old_nm, "`.", 
             "\nTo overwrite, set `overwrite = TRUE`.", 
             call. = FALSE)
      }
    }

    uniprot_entrez_lookup <- c(curr_lookup, uniprot_entrez_lookup) %>% 
      .[!duplicated(.)]
  }
  
  save(uniprot_entrez_lookup, file = file)
}


