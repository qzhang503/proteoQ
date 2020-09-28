#' Generate simulated UniProt results from RefSeq PSM data 
#' 
#' @inheritParams load_expts
#' @param type The type of PSM data in \code{c("Masoct", "MaxQuant")}.
#' @import dplyr tidyr stringr
#' @importFrom magrittr %>% %T>% %$% %<>% 
#' @export
simulUniprotPSM <- function(type = "Mascot", dat_dir = NULL) {
  type <- rlang::as_string(rlang::enexpr(type)) %>% 
    tolower()
  
  if (is.null(dat_dir)) dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv),
                                            error = function(e) 1)
  
  if (dat_dir == 1) stop("Set up the working directory first.", call. = FALSE)
  
  switch(type,
         mascot = mascot_refseq2uniprot(dat_dir), 
         maxquant = maxquant_refseq2uniprot(dat_dir), 
         stop("Unknown `type = ", type, "`", call. = FALSE))
}


#' Helper function 
#' 
#' @inheritParams load_expts
mascot_refseq2uniprot <- function(dat_dir) {
  dir.create(file.path(dat_dir, "Mascot_refseq"))
  local({
    filelist <- list.files(path = file.path(dat_dir), pattern = "^F[0-9]+.csv$")
    file.copy(file.path(dat_dir, filelist), file.path(dat_dir, "Mascot_refseq", filelist))
  })
  
  suppressMessages(rmPSMHeaders())
  
  filelist <- list.files(path = file.path(dat_dir, "PSM/cache"), pattern = "^F[0-9]{6}_hdr_rm.csv$")
  filelist_hdr <- list.files(path = file.path(dat_dir, "PSM/cache"), pattern = "^F[0-9]{6}_header.txt$")
  
  purrr::walk2(filelist, filelist_hdr, ~ {
    df <- read.delim(file.path(dat_dir, "PSM/cache", .x), sep = ',', check.names = FALSE, 
                     header = TRUE, stringsAsFactors = FALSE, quote = "\"",fill = TRUE , skip = 0)
    df$psm_index <- seq_along(1:nrow(df))
    
    df_acc <- data.frame(prot_acc_orig = df$prot_acc, 
                         database = gsub("(.*)::.*", "\\1", df$prot_acc), 
                         prot_acc = gsub("[1-9]{1}::", "", df$prot_acc), 
                         prot_desc = df$prot_desc, 
                         psm_index = df$psm_index) %>% 
      dplyr::mutate(prot_acc = as.character(prot_acc))
      
    
    df_acc_hs <- df_acc %>% dplyr::filter(grepl("\\[Homo sapiens\\]", prot_desc))
    df_acc_mm <- df_acc %>% dplyr::filter(grepl("\\[Mus musculus\\]", prot_desc))
    df_others <- df_acc %>% dplyr::filter(!grepl("\\[Homo sapiens\\]|\\[Mus musculus\\]", prot_desc))
    
    uniprot_acc_hs <- local({
      accessions <- annot_from_to(abbr_species = "Hs", 
                                  keys = as.character(df_acc_hs$prot_acc), 
                                  from = "REFSEQ", 
                                  to = "UNIPROT") %>% 
        dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
        dplyr::filter(!duplicated(prot_acc))
      
      dplyr::left_join(df_acc_hs, accessions, by = "prot_acc") %>% 
        dplyr::filter(!duplicated(psm_index)) %>% 
        dplyr::mutate(uniprot_acc = paste(database, uniprot_acc, sep = "::")) %>% 
        dplyr::select(uniprot_acc, psm_index) %>% 
        dplyr::rename(prot_acc = uniprot_acc)
    })
    
    uniprot_acc_mm <- local({
      accessions <- annot_from_to(abbr_species = "Mm", 
                                  keys = as.character(df_acc_mm$prot_acc), 
                                  from = "REFSEQ", 
                                  to = "UNIPROT") %>% 
        dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
        dplyr::filter(!duplicated(prot_acc))
      
      dplyr::left_join(df_acc_mm, accessions, by = "prot_acc") %>% 
        dplyr::filter(!duplicated(psm_index)) %>% 
        dplyr::mutate(uniprot_acc = paste(database, uniprot_acc, sep = "::")) %>% 
        dplyr::select(uniprot_acc, psm_index) %>% 
        dplyr::rename(prot_acc = uniprot_acc)
    })
    
    uniprot_acc_others <- df_others %>% 
      dplyr::select(prot_acc_orig, psm_index) %>% 
      dplyr::rename(prot_acc = prot_acc_orig) %>% 
      dplyr::mutate(prot_acc = as.character(prot_acc))
    
    uniprot_acc <- dplyr::bind_rows(
      uniprot_acc_hs,
      uniprot_acc_mm,
      uniprot_acc_others
    ) %>% 
      dplyr::arrange(psm_index) %>% 
      dplyr::select(prot_acc)
    
    df[["prot_acc"]] <- unlist(uniprot_acc)
    df[["psm_index"]] <- NULL
    
    dir.create(file.path(dat_dir, "PSM/cache/temp"))
    write.table(df, file = file.path(dat_dir, "PSM/cache/temp", .x), sep = ",")
    df <- readLines(file.path(dat_dir, "PSM/cache/temp", .x))
    hdr <- readLines(file.path(dat_dir, "PSM/cache", .y))
    
    write(append(hdr, df), file = file.path(dat_dir, gsub("_hdr_rm", "", .x)))
  })
  
  unlink(file.path(dat_dir, "PSM/cache/temp"), recursive = TRUE)
}


#' Helper function 
#' 
#' @inheritParams load_expts
maxquant_refseq2uniprot <- function(dat_dir) {
  dir.create(file.path(dat_dir, "Maxquant_refseq"))
  filelist <- list.files(path = file.path(dat_dir), pattern = "^msms.*\\.txt$")
  file.copy(file.path(dat_dir, filelist), file.path(dat_dir, "Maxquant_refseq", filelist))
  
  purrr::walk(filelist, ~ {
    df <- read.csv(file.path(dat_dir, .x), sep = "\t", check.names = FALSE, header = TRUE, comment.char = "#")
    df$psm_index <- seq_along(1:nrow(df))
    
    suppressWarnings(
      df_acc <- df$Proteins %>% 
        stringr::str_split(";", simplify = TRUE) %>% 
        `colnames<-`(paste0("prot_acc", 1:ncol(.))) %>% 
        data.frame() %>% 
        dplyr::mutate(psm_index = df$psm_index) %>% 
        tidyr::gather("Index", "prot_acc", -psm_index) %>% 
        dplyr::filter(nchar(prot_acc) != 0) %>% 
        dplyr::select(-Index) %>% 
        dplyr::arrange(psm_index)      
    )

    accessions_hs <- annot_from_to(abbr_species = "Hs", 
                                   keys = as.character(df_acc$prot_acc), 
                                   from = "REFSEQ", 
                                   to = "UNIPROT") %>% 
      dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_acc, ., by = "prot_acc") %>% 
      dplyr::filter(!is.na(uniprot_acc)) %>% 
      dplyr::select(-prot_acc) 
    
    accessions_mm <- annot_from_to(abbr_species = "Mm", 
                                   keys = as.character(df_acc$prot_acc), 
                                   from = "REFSEQ", 
                                   to = "UNIPROT") %>% 
      dplyr::rename(prot_acc = "REFSEQ", uniprot_acc = "UNIPROT") %>% 
      dplyr::filter(!duplicated(prot_acc)) %>% 
      dplyr::left_join(df_acc, ., by = "prot_acc") %>% 
      dplyr::filter(!is.na(uniprot_acc)) %>% 
      dplyr::select(-prot_acc) 
    
    accessions <- dplyr::bind_rows(accessions_hs, accessions_mm) %>% 
      dplyr::arrange(psm_index) %>% 
      tidyr::unite(id, c("uniprot_acc", "psm_index"), remove = FALSE) %>% 
      dplyr::filter(!duplicated(id)) %>% 
      dplyr::select(-id) %>% 
      dplyr::group_by(psm_index) %>% 
      dplyr::summarise(uniprot_acc = paste(uniprot_acc, collapse = ";"))
    
    df_others <- df %>% 
      dplyr::filter(!psm_index %in% unique(accessions$psm_index)) %>% 
      dplyr::select(-psm_index) %>% 
      dplyr::mutate(Proteins = as.character(Proteins))
    
    df_hsmm <- df %>% 
      dplyr::filter(psm_index %in% unique(accessions$psm_index)) %>% 
      dplyr::left_join(accessions) %>% 
      dplyr::mutate(Proteins = uniprot_acc) %>% 
      dplyr::select(-uniprot_acc, -psm_index)
    
    write.table(dplyr::bind_rows(df_hsmm, df_others), file = file.path(dat_dir, .x), sep = "\t", row.names = FALSE)
  })
}


