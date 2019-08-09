#'STRING outputs
#'
#'\code{getStringDB} reports the \code{\href{https://string-db.org/}{STRING}}
#'results of protein-protein interaction and protein expressions. The
#'interaction file, \code{Protein_STRING_ppi.tsv}, and the expression file,
#'\code{Protein_STRING_expr.tsv}, can be used by
#'\code{\href{https://cytoscape.org/}{Cytoscape}} for further processing and
#'analysis.
#'
#'@param score_cutoff Numeric; the threshold of the \code{combined_score} of
#'  protein-protein interaction.
#'@param nseq_cutoff Positive integer. Proteins with the number of identifying
#'  peptides smaller than the cut-off value will be removed from further
#'  analysis.
#' @examples
#' getStringDB(
#'   score_cutoff = .9,
#'   nseq_cutoff = 2,
#'   adjP = FALSE
#' )
#'
#'@inheritParams proteoVolcano
#'@inheritParams proteoPurge
#'@import dplyr purrr STRINGdb
#'@export
getStringDB <- function(score_cutoff = .7, nseq_cutoff = 2, adjP = FALSE) {
  stopifnot(rlang::is_logical(adjP))
  stopifnot(rlang::is_double(score_cutoff))
  stopifnot(rlang::is_double(nseq_cutoff))
  
  id <- match_identifier("gene")
  
  taxid <- switch(find_species(), 
                  human = 9606, 
                  mouse = 10090, 
                  rat = 10116, 
                  stop("Unknown `species`.", Call. = FALSE)
  )
  
  if (score_cutoff <= 1) score_cutoff <- score_cutoff * 100
  
  dir.create(file.path(dat_dir, "Protein\\String\\cache"), recursive = TRUE, showWarnings = FALSE)

  df <- read.csv(file.path(dat_dir, "Protein\\Model\\Protein_pVals.txt"), check.names = FALSE, 
                 stringsAsFactors = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::filter(!is.na(id), nchar(id) > 0, n_pep >= nseq_cutoff) 
  
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
  df <- df %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.character) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), ~ gsub("\\s*", "", .x) ) %>% 
    dplyr::mutate_at(vars(grep("pVal|adjP", names(.))), as.numeric)

  string_db <- STRINGdb$new(version = "10", species = taxid, score_threshold = 0, 
                            input_directory = file.path(dat_dir, "Protein\\String\\cache") %>% 
                              gsub("\\/", "\\\\", .))

  example1_mapped <- string_db$map(df, id, removeUnmappedRows = TRUE)
  
  nm_string <- file.path(dat_dir, "Protein\\String\\cache\\Protein_STRING_interactions.tsv")
  if (file.exists(nm_string)) {
    ppi <- read.csv(nm_string, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")
  } else {
    ppi <- string_db$get_interactions(example1_mapped$STRING_id)
    write.table(ppi, nm_string, quote = FALSE, sep = "\t", row.names = FALSE)
  }

  ppi <- ppi[ppi$combined_score > score_cutoff, ]
  
  #node1
  ppi <- example1_mapped[, c(id, "STRING_id")] %>% 
    `names_pos<-`(2, "from") %>% 
    dplyr::full_join(ppi, by = "from")
  
  ppi <- bind_cols(
    ppi %>% dplyr::select(id), 
    ppi %>% dplyr::select(-id)
  ) %>% 
    `names_pos<-`(1, "node1")
  
  #node2
  ppi <- example1_mapped[, c(id, "STRING_id")] %>% 
    `names_pos<-`(2, "to") %>% 
    dplyr::full_join(ppi, by = "to")
  
  ppi <- bind_cols(
    ppi %>% dplyr::select(id), 
    ppi %>% dplyr::select(-id)
  ) %>% 
    `names_pos<-`(1, "node2")
  
  # re-arrange
  ppi <- dplyr::bind_cols(
    ppi %>% dplyr::select(c("node1", "node2", "from", "to")), 
    ppi %>% dplyr::select(-c("node1", "node2", "from", "to"))
  ) %>% 
    `names_pos<-`(c(3,4), c("node1_external_id", "node2_external_id")) %>% 
    dplyr::mutate_at(vars(5:ncol(.)), ~ .x/1000) %>% 
    `names_pos<-`(1, "#node1")

  write.table(ppi, file.path(dat_dir, "Protein\\String\\Protein_STRING_ppi.tsv"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # expression data file
  df <- dplyr::bind_cols(
    df %>% dplyr::select(id), 
    df %>% dplyr::select(grep("pVal|_log2_R", names(.)))
  )
  
  df[is.na(df)] <- "" # Cytoscape compatibility
  write.table(df, file.path(dat_dir, "Protein\\String", "Protein_STRING_expr.tsv"), 
              quote = FALSE, sep = "\t", row.names = FALSE)

}
