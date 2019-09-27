#' Adds protein description
#'
#' \code{annotPrndesc} annotates protein descriptions based on the \code{fasta}
#'
#' @import plyr dplyr purrr rlang stringr seqinr
#' @importFrom magrittr %>% %$%
annotPrndesc <- function (df, fasta){
  stopifnot("prot_acc" %in% names(df))
  
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  acc_type <- find_acctype() %>% tolower()
  
  if (acc_type == "refseq_acc") {
    key <- "refseq_acc"
  } else if (acc_type %in% "uniprot_id") {
    key <- "uniprot_id"
  } else if (acc_type %in% "uniprot_acc") {
    key <- "uniprot_acc"
  } else {
    warning("Unkown accession type.")
  }
  
  df <- df %>% 
    dplyr::mutate(prot_acc = gsub("^.*\\|(.*)\\|.*$", "\\1", prot_acc))
  
  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- purrr::map(fasta, ~ {
        seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
      }) %>% do.call(`c`, .) %>% 
        `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
        .[names(.) %in% unique(df$prot_acc)]
      
      if (length(fasta) == 0) {
        stop("No fasta entries match protein accessions; probably wrong fasta file.", 
             call. = FALSE)
      }
      
      lookup <- dplyr::bind_cols(
        # prot_acc = seqinr::getName(fasta), 
        prot_acc = names(fasta), 
        prot_desc = seqinr::getAnnot(fasta) %>% 
          purrr::map(., `[[`, 1) %>% 
          unlist(), 
        prot_mass = purrr::map_dbl(fasta, ~ {seqinr::getSequence(.x) %>% seqinr::pmw()})
      ) %>% 
        dplyr::filter(.$prot_acc %in% unique(df$prot_acc))
    } else {
      stop("Wrong FASTA file path or name.", call. = FALSE)
    }
  } else {
    stop("FASTA file not provided.")
  }
  
  rm(fasta)
  
  df %>% 
    dplyr::left_join(lookup, by = "prot_acc") %>% 
    dplyr::mutate(prot_mass = round(prot_mass, digits = 1))
}


#' Add start and end positions and preceding and succeeding residues of peptides
#'
#' \code{annotPeppos_mq} annotates start and end positions and the preceding and
#' succeeding residues of peptides.
#'
#' @import dplyr purrr rlang stringr seqinr
#' @importFrom magrittr %>% %$%
annotPeppos_mq <- function (df, fasta){
  stopifnot(all(c("prot_acc", "pep_seq") %in% names(df)))
  
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)
  
  if (acc_type == "refseq_acc") {
    key <- "refseq_acc"
  } else if (acc_type %in% "uniprot_id") {
    key <- "uniprot_id"
  } else if (acc_type %in% "uniprot_acc") {
    key <- "uniprot_acc"
  } else {
    warning("Unkown accession type.")
  }
  
  df <- df %>% 
    dplyr::mutate(prot_acc = gsub("^.*\\|(.*)\\|.*$", "\\1", prot_acc))
  
  if (!is.null(fasta)) {
    if (all(file.exists(fasta))) {
      fasta <- purrr::map(fasta, ~ {
        seqinr::read.fasta(.x, seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
      }) %>% do.call(`c`, .) %>% 
        `names<-`(gsub("^.*\\|(.*)\\|.*$", "\\1", names(.))) %>% 
        .[names(.) %in% unique(df$prot_acc)]
      
      if (length(fasta) == 0) {
        stop("No fasta entries match protein accessions; probably wrong fasta files.", 
             call. = FALSE)
      }
      
      pep_pos_all <- purrr::map2(as.list(df$prot_acc), as.list(df$pep_seq), ~ {
        fasta_sub <- fasta %>% .[names(.) == .x]
        pep_seq <- as.character(.y)
        
        if (!rlang::is_empty(fasta_sub)) {
          pep_pos <- str_locate(fasta_sub, pattern = pep_seq)
          pos_bf <- pep_pos[1] - 1
          pos_af <- pep_pos[2] + 1
          
          pep_res_before <- str_sub(fasta_sub, pos_bf, pos_bf)
          pep_res_after <- str_sub(fasta_sub, pos_af, pos_af)
          
          if (nchar(pep_res_before) == 0) pep_res_before <- "-"
          if (nchar(pep_res_after) == 0) pep_res_after <- "-"
          
          pep_pos <- cbind(pep_seq, pep_res_before, pep_pos, pep_res_after)
        } else {
          pep_pos <- cbind(pep_seq, pep_res_before = NA, start = NA, end = NA, pep_res_after = NA)
        }
      }) %>% 
        do.call(rbind, .) %>% 
        `colnames<-`(c("pep_seq", "pep_res_before", "pep_start", "pep_end", "pep_res_after")) %>% 
        data.frame(check.names = FALSE)
    } else {
      stop("Wrong FASTA file path or names.", call. = FALSE)
    }
  } else {
    stop("FASTA file not provided.")
  }
  
  rm(fasta)
  
  df %>% dplyr::left_join(pep_pos_all, by = "pep_seq")
}


#' Keratin removals based on dbs
#'
#' @param id. 
#'
#' @import dplyr purrr rlang
#' @importFrom magrittr %>%
rm_krts <- function (df, acc_type) {
  stopifnot("prot_acc" %in% names(df))
  
  krts <- dbs$prn_annot %>%
    dplyr::filter(grepl("^krt[0-9]+", .$gene, ignore.case = TRUE)) %>%
    dplyr::filter(!is.na(.[[acc_type]])) %>%
    dplyr::select(acc_type) %>%
    unlist()
  
  df %>%
    dplyr::filter(! .$prot_acc %in% krts)  
}




#' Peptide reports for individual TMT experiments
#'
#' \code{normPepSplex_mqpep} processes peptide data from MaxQuant peptide
#' table(s). It splits the original tables by TMT set numbers. Different LCMS
#' injections under the same TMT experiment will be handled as different
#' fractions. Therefore, the injection index is always "1".
#'
#' @inheritParams normPep
#' @return Results in \code{.txt} files for each of TMT experiments and LC/MS
#'   injections.
#'
#' @import stringr dplyr purrr rlang  magrittr
normPepSplex_mqpep <- function (id = "pep_seq_mod", fasta = fasta, 
                                pep_unique_by = "group", corrected_mq_int = TRUE, 
                                rm_mq_krts = TRUE, rm_mq_craps = TRUE, rm_mq_reverses = TRUE) {
  
  id <- rlang::as_string(rlang::enexpr(id))
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  TMT_plex <- TMT_plex(label_scheme)
  TMT_levels <- TMT_levels(TMT_plex)
  
  if (id == "pep_seq_mod") {
    filelist <- list.files(path = file.path(dat_dir), pattern = "^modificationSpecificPeptides.*\\.txt$")
  } else if (id == "pep_seq") {
    filelist <- list.files(path = file.path(dat_dir), pattern = "^peptides.*\\.txt$")
  }
  
  if (rlang::is_empty(filelist)) {
    stop("No MaxQuant peptide results were found. 
         Make sure the file names start with `modificationSpecificPeptides` at `id = pep_seq_mod` or 
         `peptides` at `id = pep_seq`.", call. = FALSE)
  }
  
  purrr::walk(filelist, ~ {
    df <- read.csv(file.path(dat_dir, .x), check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") 
    
    df_psm <- df %>% 
      dplyr::select(Sequence, grep("^Experiment\\s{1}", names(.)))
    
    df <- df %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}count\\s{1}", names(.))) %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}\\d+$", names(.))) %>% 
      dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}corrected\\s{1}\\d+$", names(.))) %>% 
      dplyr::select(-grep("^Intensity\\s{1}", names(.))) %>% 
      dplyr::select(-grep("^Experiment\\s{1}", names(.))) %>% 
      dplyr::select(-grep("^Fraction\\s{1}", names(.)))
    
    if (id == "pep_seq_mod") {
      if (all(names(df) != "Modifications")) {
        stop("Column `Modifications` not found; use modification-specific peptide table(s) or set `id = pep_seq`.", 
             call. = FALSE)
      }
    } else if (id == "pep_seq") {
      if (any(names(df) == "Modifications")) {
        stop("Use peptide table(s) without column `Modifications` or set `id = pep_seq_mod`.", 
             call. = FALSE)
      }
    }
    
    if (corrected_mq_int) {
      df <- df %>% 
        dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}\\d+\\s{1}", names(.)))
    } else {
      df <- df %>% 
        dplyr::select(-grep("^Reporter\\s{1}intensity\\s{1}corrected\\s{1}\\d+[^\\d]", names(.)))
    }
    
    if (pep_unique_by == "group") {
      df <- df %>% 
        dplyr::mutate(pep_isunique = ifelse(`Unique (Groups)` == "yes", 1, 0))
    } else if (pep_unique_by == "protein") {
      df <- df %>% 
        dplyr::mutate(pep_isunique = ifelse(`Unique (Proteins)` == "yes", 1, 0))
    } else {
      df %>% 
        dplyr::mutate(pep_isunique == 1)
    }
    
    mq_grps_in_df <- names(df) %>% 
      .[grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+", .)] %>% 
      gsub("Reporter\\s{1}intensity\\s{1}.*[0-9]+\\s{1}(.*)", "\\1", .) %>% 
      unique()
    
    for (set_idx in unique(label_scheme$TMT_Set)) {
      mq_grp_allowed <- label_scheme %>% 
        dplyr::filter(TMT_Set == set_idx) %>% 
        dplyr::select(MQ_Experiment) %>% 
        unlist() %>% 
        unique()
      
      mq_grps_matched <- mq_grps_in_df %>% 
        .[. %in% mq_grp_allowed]
      
      if (purrr::is_empty(mq_grps_matched)) next
      
      df_int <- df[, grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+", names(df)), drop = FALSE] %>% 
        dplyr::select(grep(paste0(mq_grps_matched, "$"), names(.))) %>% 
        `names<-`(gsub("TMT-", "I", as.character(TMT_levels)))
      
      stopifnot(ncol(df_int) == TMT_plex)
      
      df_int <- df_int %>% 
        logfcPep(label_scheme, set_idx)
      
      df_psm_sub <- df_psm %>% 
        dplyr::select(grep(paste0(mq_grps_matched, "$"), names(.))) %>% 
        `colnames<-`("n_psm")
      
      df_sub <- dplyr::bind_cols(
        df[, !grepl("Reporter\\s{1}intensity\\s{1}.*[0-9]+", names(df)), drop = FALSE], 
        df_psm_sub, 
        df_int
      ) %>% 
        dplyr::mutate(TMT_Set = set_idx)
      
      if (id == "pep_seq_mod") {
        df_sub <- df_sub %>% 
          dplyr::rename(pep_seq = Sequence, 
                        prot_acc = Proteins, 
                        pep_score = Score, 
                        pep_expect = PEP) %>% 
          dplyr::mutate(prot_acc = gsub("\\;.*", "", prot_acc)) %>% 
          dplyr::mutate(pep_seq_mod = paste0(pep_seq, " [", Modifications, "]"))
      } else if (id == "pep_seq") {
        df_sub <- df_sub %>% 
          dplyr::rename(pep_seq = Sequence, 
                        prot_acc = `Leading razor protein`, 
                        pep_score = Score, 
                        pep_expect = PEP) %>% 
          dplyr::mutate(pep_seq_mod = pep_seq)
      }
      
      df_sub <- df_sub %>% 
        dplyr::filter(.$pep_isunique == 1) %>% 
        dplyr::select(-pep_isunique)
      
      df_sub <- annotPrndesc(df_sub, fasta)
      
      if(rm_mq_krts) {
        df_sub <- df_sub %>% 
          dplyr::mutate(is_krt = grepl("^.*\\s+krt[0-9]+", prot_desc, ignore.case = TRUE)) %>% 
          dplyr::filter(!is_krt) %>% 
          dplyr::select(-is_krt)
      }
      
      df_sub <- annotPeppos(df_sub, fasta)
      
      if (rm_mq_craps) {
        df_sub <- df_sub %>% 
          dplyr::filter(.$`Potential contaminant` != "+")
      }
      
      if (rm_mq_reverses) {
        df_sub <- df_sub %>% 
          dplyr::filter(.$Reverse != "+")
      }
      
      df_sub <- df_sub %>% 
        dplyr::select(pep_seq, pep_seq_mod, n_psm, prot_acc, prot_desc, 
                      prot_mass, pep_start, pep_end, pep_score, pep_expect) %>% 
        dplyr::bind_cols(df_sub %>% dplyr::select(grep("^log2_R[0-9]{3}", names(df_sub)))) %>% 
        dplyr::bind_cols(df_sub %>% dplyr::select(grep("^I[0-9]{3}", names(df_sub)))) %>% 
        dplyr::bind_cols(df_sub %>% dplyr::select(grep("^N_log2_R[0-9]{3}", names(df_sub)))) %>% 
        dplyr::bind_cols(df_sub %>% dplyr::select(grep("^N_I[0-9]{3}", names(df_sub)))) %>% 
        dplyr::bind_cols(df_sub %>% dplyr::select(TMT_Set)) %>% 
        dplyr::filter(!grepl("^REV_", .$prot_acc)) %>% 
        dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0)
      
      if(id == "pep_seq_mod") {
        df_sub <- df_sub %>% dplyr::select(-pep_seq)
      } else if (id == "pep_seq") {
        df_sub <- df_sub %>% dplyr::select(-pep_seq_mod)
      } else {
        stop("Unknown `id`.", call. = FALSE)
      }
      
      injn_idx <- 1
      write.table(df_sub, 
                  file.path(dat_dir, "Peptide", paste0("TMTset", set_idx, "_LCMSinj", injn_idx, "_Peptide_N.txt")), 
                  sep = "\t", col.names = TRUE, row.names = FALSE)
    }
  })
  
}


