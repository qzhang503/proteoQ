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




#' Violin plots of SDs
#'
#' \code{sd_violin_full} visualizes the SD distribution of SD
#'
#' @import dplyr purrr rlang ggplot2
#' @importFrom magrittr %>%
sd_violin_full <- function(df_sd, id, label_scheme, filepath, filename) {
  id <- rlang::as_string(rlang::enexpr(id))
  
  Levels <- names(df_sd) %>% 
    .[grepl("^log2_R[0-9]{3}[NC]*\\s+\\(", .)] %>% 
    gsub("^log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", .)  
  
  df_sd <- df_sd %>%
    `names<-`(gsub("^log2_R[0-9]{3}[NC]*\\s+\\((.*)\\)$", "\\1", names(.))) %>% 
    tidyr::gather(key = !!rlang::sym(id), value = "SD") %>%
    dplyr::rename(Channel := !!rlang::sym(id)) %>% 
    dplyr::ungroup(Channel) %>% 
    dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
    dplyr::filter(!is.na(SD))
  
  p <- ggplot() +
    geom_violin(df_sd, mapping = aes(x = Channel, y = SD, fill = Channel), size = .25) +
    geom_boxplot(df_sd, mapping = aes(x = Channel, y = SD), width = 0.1, lwd = .2, fill = "white") +
    stat_summary(df_sd, mapping = aes(x = Channel, y = SD), fun.y = "mean", geom = "point",
                 shape=23, size=2, fill="white", alpha=.5) +
    labs(title = expression(""), x = expression("Channel"), y = expression("SD ("*log[2]*"FC)")) +
    scale_y_continuous(limits = c(0, .6), breaks = seq(0, .6, .2)) +
    theme_psm_violin
  
  try(ggsave(file.path(filepath, filename), p, width = 7* n_TMT_sets(label_scheme), height = 7, units = "in"))
}


#'Rmove SD outliers
#'
#'\code{purgeData} removes entries with SD > cv_cutoff
#'
#'@import dplyr purrr rlang
#'@importFrom magrittr %>%
#'@importFrom magrittr %T>%
purgeData <- function(df, id, cv_cutoff = NULL, nseq_cutoff = 1, ...) {
  id <- rlang::as_string(rlang::enexpr(id))
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  
  stopifnot(nseq_cutoff > 0 & nseq_cutoff%%1 == 0)
  
  df <- df %>% dplyr::arrange(!!rlang::sym(id))
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    filelist <- list.files(path = file.path(dat_dir, "PSM"), pattern = "*_PSM_N\\.txt$") %>%
      reorder_files(n_TMT_sets(label_scheme_full))
    
    if (is_empty(filelist)) {
      stop("PSM files not found; may be you are starting with MaxQuant peptide outputs.", call. = FALSE)
    }
    
    df_sd <- purrr::map(as.list(filelist), ~ {
      df <- read.csv(file.path(dat_dir, "PSM", .x), check.names = FALSE, header = TRUE,
                     sep = "\t", comment.char = "#") %>% 
        dplyr::mutate(TMT_Set = as.integer(gsub("^TMTset(\\d+)_LCMSinj.*", "\\1", .x)))
      
      return(df)
    }) %>% do.call("rbind", .)
    
    df_sd <- df_sd %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::select(!!rlang::sym(id), TMT_Set, grep("^log2_R[0-9]{3}", names(.))) %>%
      dplyr::group_by(!!rlang::sym(id), TMT_Set) %>%
      dplyr::summarise_at(vars(starts_with("log2_R")), ~ sd(.x, na.rm = TRUE)) %>% 
      dplyr::arrange(TMT_Set) %>% 
      tidyr::gather(grep("log2_R[0-9]{3}", names(.)), key = ID, value = value) %>%
      tidyr::unite(ID, ID, TMT_Set)
    
    Levels <- unique(df_sd$ID)
    df_sd <- df_sd %>%
      dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
      tidyr::spread(ID, value) %>% 
      dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3))
    
    for (set_idx in seq_len(n_TMT_sets(label_scheme))) {
      df_sd <- newColnames(set_idx, df_sd, label_scheme)
    }
    
    df_sd <- df_sd %>% 
      dplyr::select(!!rlang::sym(id), grep("log2_R[0-9]{3}[NC]*", names(.))) %>% 
      dplyr::arrange(!!rlang::sym(id)) %T>%
      write.csv(file.path(dat_dir, "Peptide\\cache", "pep_sd.csv"), row.names = FALSE)
  } else if (id %in% c("prot_acc", "gene")) {
    df_sd <- read.csv(file.path(dat_dir, "Peptide", "Peptide.txt"), check.names = FALSE, 
                      header = TRUE, sep = "\t", comment.char = "#") %>% 
      filter(rowSums(!is.na(.[grep("^log2_R[0-9]{3}", names(.))])) > 0) %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::select(!!rlang::sym(id), grep("^log2_R[0-9]{3}", names(.))) %>%
      dplyr::group_by(!!rlang::sym(id)) %>%
      dplyr::summarise_at(vars(starts_with("log2_R")), ~ sd(.x, na.rm = TRUE)) %T>%
      write.csv(file.path(dat_dir, "Protein\\cache", "prn_sd.csv"), row.names = FALSE) %>% 
      tidyr::gather(grep("log2_R[0-9]{3}", names(.)), key = ID, value = value)
    
    Levels <- unique(df_sd$ID)
    df_sd <- df_sd %>%
      dplyr::mutate(ID = factor(ID, levels = Levels)) %>%
      tidyr::spread(ID, value) %>% 
      dplyr::mutate_at(vars(grep("^log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %>% 
      dplyr::arrange(!!rlang::sym(id))
  }
  
  # purging
  if (!is.null(cv_cutoff)) {
    stopifnot(is.numeric(cv_cutoff))
    
    df_sd <- df_sd %>% 
      dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, .x > cv_cutoff, NA)) 
    
    df_sd_lgl <- df_sd%>% 
      dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ replace(.x, !is.na(.x), 1)) %>% 
      dplyr::filter(!!rlang::sym(id) %in% df[[id]]) %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      tibble::column_to_rownames(id) %>% 
      `names<-`(paste0("sd_", names(.))) %>% 
      tibble::rownames_to_column(id)
    
    df <- df %>% 
      dplyr::arrange(!!rlang::sym(id)) %>% 
      dplyr::left_join(df_sd_lgl, by = id)
    
    df[, grepl("^log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^Z_log2_R[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^Z_log2_R[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df[, grepl("^N_I[0-9]{3}", names(df))] <-
      purrr::map2(as.list(df[, grepl("^N_I[0-9]{3}", names(df))]),
                  as.list(df[, grepl("^sd_log2_R[0-9]{3}", names(df))]), `*`) %>%
      dplyr::bind_rows()
    
    df <- df %>% 
      dplyr::select(-grep("^sd_log2_R", names(.)))
  }
  
  if (id %in% c("pep_seq", "pep_seq_mod")) {
    filepath <- file.path(dat_dir, "Peptide\\Purge")
    filename <- "Peptide_SD.png"
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
    df <- df %>% dplyr::filter(n_psm >= nseq_cutoff)
  } else if (id %in% c("prot_acc", "gene")) {
    filepath <- file.path(dat_dir, "Protein\\Purge")
    filename <- "Protein_SD.png"
    dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
    df <- df %>% dplyr::filter(n_pep >= nseq_cutoff)
  }
  
  df <- df %>% 
    dplyr::filter(rowSums(!is.na(.[, grep("^log2_R[0-9]{3}", names(.) )])) > 0)
  
  # violin plots
  df_sd <- df_sd %>% 
    dplyr::filter(!!rlang::sym(id) %in% df[[id]])
  
  sd_violin_full(df_sd, !!id, label_scheme, filepath, filename)
  
  return(df)
}


#'Purge data
#'
#'\code{proteoPurge} provides additional cleanup of protein or peptide data.
#'
#'The function matches the current \code{id} to those in the latest \code{call}
#'to \code{\link{normPep}} or \code{\link{normPrn}}.  For example, if
#'\code{pep_seq} was used in \code{\link{normPep}}, the current \code{id =
#'pep_seq_mod} will be matched to \code{id = pep_seq}.
#'
#'@inheritParams proteoHist
#'@param cv_cutoff Numeric; the cut-off of CV. The CVs are from the ascribing
#'  PSMs for peptide data or ascribing peptides for protein data.
#'@param nseq_cutoff Positive integer. When calling from \code{purPep}, peptide
#'  entries in \code{Peptide.txt} with the number of identifying PSMs smaller
#'  than \code{nseq_cutoff} will be replaced with NA. When calling from
#'  \code{purPrn}, protein entries in \code{Protein.txt} with the number of
#'  identifying peptides smaller than \code{nseq_cutoff} will be replaced with
#'  NA.
#'@param ... Additional parameters for plotting yet to be defined.
#'@import dplyr rlang ggplot2
#'@importFrom magrittr %>%
#'@export
proteoPurge <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene"), 
                         cv_cutoff = NULL, nseq_cutoff = 1, 
                         df = NULL, filepath = NULL, filename = NULL, ...) {
  
  id <- rlang::enexpr(id)
  if(length(id) != 1) id <- rlang::expr(gene)
  stopifnot(rlang::as_string(id) %in% c("pep_seq", "pep_seq_mod", "prot_acc", "gene"))
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)
  
  reload_expts()
  
  info_anal(id = !!id, 
            df = !!df, filepath = !!filepath, filename = !!filename,
            anal_type = "Purge")(cv_cutoff = cv_cutoff, nseq_cutoff = nseq_cutoff, ...)
}


#'Purge peptide data
#'
#'\code{purPep} is a wrapper of \code{\link{proteoPurge}} for peptide data
#'
#'@rdname proteoPurge
#'
#' @examples
#' \dontrun{
#' purPep(
#'   cv_cutoff = 0.4,
#'   nseq_cutoff = 1,
#' )
#' }
#'
#'@import purrr
#'@export
purPep <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Peptide\\Purge\\log"), recursive = TRUE, showWarnings = FALSE)
  quietly_log <- purrr::quietly(proteoPurge)(id = pep_seq, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Peptide\\Purge\\log","pepPurge_log.csv"), append = TRUE)
}


#'Purge protein data
#'
#'\code{purPrn} is a wrapper of \code{\link{proteoPurge}} for protein data
#'
#'@rdname proteoPurge
#'
#' @examples
#' \dontrun{
#' purPrn(
#'   cv_cutoff = 0.4,
#'   nseq_cutoff = 1,
#' )
#' }
#'
#'@import purrr
#'@export
purPrn <- function (...) {
  err_msg <- "Don't call the function with argument `id`.\n"
  if(any(names(rlang::enexprs(...)) %in% c("id"))) stop(err_msg)
  
  dir.create(file.path(dat_dir, "Protein\\Purge\\log"), recursive = TRUE, showWarnings = FALSE)
  
  quietly_log <- purrr::quietly(proteoPurge)(id = gene, ...)
  purrr::walk(quietly_log, write, 
              file.path(dat_dir, "Protein\\Purge\\log","prnPurge_log.csv"), append = TRUE)
}








