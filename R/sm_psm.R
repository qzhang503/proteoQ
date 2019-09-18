#' Splits PSM tables
#'
#' \code{splitPSM} splits the PSM outputs after \code{rmPSMHeaders()}. It
#' separates PSM data by TMT experiment and LC/MS injection.
#'
#' @inheritParams splitPSM
#' @inheritParams splitPSM_mq
#' @examples
#' splitPSM_sm(rptr_intco = 1000, rm_craps = TRUE, plot_rptr_int = TRUE)
#'
#' @import dplyr tidyr readr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
splitPSM_sm <- function(fasta = NULL, rm_craps = FALSE, rm_krts = FALSE, rptr_intco = 1000, 
                        annot_kinases = FALSE, plot_rptr_int = TRUE, ...) {
                        
  add_mascot_headers <- function (df) {
    df %>% 
      dplyr::mutate(pep_scan_title = RAW_File) %>% 
      dplyr::mutate(prot_hit_num = NA) %>% 
      dplyr::mutate(prot_family_member = NA) %>% 
      dplyr::mutate(prot_score = NA) %>% 
      dplyr::mutate(prot_matches = NA) %>% 
      dplyr::mutate(prot_matches_sig = NA) %>% 
      dplyr::mutate(prot_sequences = NA) %>% 
      dplyr::mutate(prot_sequences_sig = NA) %>% 
      dplyr::mutate(pep_query = NA) %>% 
      dplyr::mutate(pep_rank = NA) %>% 
      dplyr::mutate(pep_isbold = NA) %>% 
      dplyr::mutate(pep_exp_mr = NA) %>% 
      dplyr::mutate(pep_exp_z = NA) %>% 
      dplyr::mutate(pep_calc_mr = NA) %>% 
      dplyr::mutate(pep_exp_mz = NA) %>% 
      dplyr::mutate(pep_delta = NA) %>% 
      dplyr::mutate(pep_miss = NA) %>% 
      dplyr::mutate(pep_var_mod = NA) %>% 
      dplyr::mutate(pep_var_mod_pos = NA) %>% 
      dplyr::mutate(pep_summed_mod_pos = NA) %>% 
      dplyr::mutate(pep_local_mod_pos = NA)
  }
  
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  
  on.exit(message("Split PSM by sample IDs and LCMS injections --- Completed."), add = TRUE)
  
  if (is.null(fasta)) stop("FASTA file not provided.", call. = FALSE)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  load(file = file.path(dat_dir, "fraction_scheme.Rdata"))
  
  TMT_plex <- TMT_plex(label_scheme_full)
  TMT_levels <- TMT_levels(TMT_plex)
  
  filelist <- list.files(path = file.path(dat_dir), pattern = "^PSMexport.*\\.ssv$")
  
  if (rlang::is_empty(filelist)) 
    stop(paste("No PSM files were found under", file.path(dat_dir), 
               "\nMake sure the names of PSM files start with `PSMexport`."), call. = FALSE)
  
  df <- purrr::map(file.path(dat_dir, filelist), readr::read_delim, delim = ";") %>% 
    dplyr::bind_rows() 

  dots <- rlang::enexprs(...)
  lang_dots <- dots %>% .[purrr::map_lgl(., is.language)] %>% .[grepl("^filter_", names(.))]
  dots <- dots %>% .[! . %in% lang_dots]
  
  cat("Available column keys for data filtration: \n")
  cat(paste0(names(df), "\n"))
  
  df <- df %>% 
    filters_in_call(!!!lang_dots) %>% 
    dplyr::rename(prot_acc = accession_number) %>% 
    dplyr::mutate(prot_acc = gsub("[1-9]{1}::", "", prot_acc)) %>% 
    annotPrn(fasta)  

  if (TMT_plex == 11) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131N", "I131C")
  } else if (TMT_plex == 10) {
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C",
                 "I130N", "I130C", "I131")
  } else if(TMT_plex == 6) {
    col_int <- c("I126", "I127", "I128", "I129", "I130", "I131")
  } else {
    col_int <- NULL
  }  
  
  if (TMT_plex > 0) {
    df <- df %>% 
      `names<-`(gsub("^TMT_([0-9]{3}[NC]?)", "I\\1", names(.))) %>% 
      dplyr::select(-grep("^I[0-9]{3}[NC]?_[0-9]{3}[NC]?$", names(.))) %>% 
      as.data.frame()

    df <- sweep(df[, col_int, drop = FALSE], 1, df[, "I126"], "/") %>% 
      `colnames<-`(gsub("I", "R", names(.))) %>% 
      dplyr::select(-R126) %>% 
      dplyr::mutate_at(vars(grep("^R[0-9]{3}", names(.))), ~ replace(.x, is.infinite(.x), NA)) %>% 
      dplyr::bind_cols(df, .) 
    
    # if `sequence` contains lower-case letters...
    # parse out `variableSites`
    df <- df %>% 
      dplyr::rename(
        pep_seq = sequence, 
        pep_exp_z = parent_charge, 
        pep_score = `score`,
        pep_var_mod = modifications, 
        RAW_File = `filename`
      ) %>% 
      dplyr::mutate(RAW_File = gsub("\\.[0-9]+\\.[0-9]+\\.[0-9]+$", "", RAW_File)) %>% 
      dplyr::mutate(pep_seq = toupper(pep_seq)) %>% # for now
      dplyr::mutate(pep_miss = ifelse(.$next_aa == "(-)", str_count(pep_seq, "[KR]"), str_count(pep_seq, "[KR]") - 1))      
  }
  
  acc_type <- df$acc_type %>% unique() %>% .[!is.na(.)] %>% as.character()
  stopifnot(length(acc_type) == 1)

  if (rm_craps) {
    data(package = "proteoQ", prn_annot_crap)
    
    craps <- prn_annot_crap %>% 
      dplyr::filter(!duplicated(.[[acc_type]])) %>% 
      dplyr::select(acc_type) %>% 
      unlist()
    
    df <- df %>% dplyr::filter(! prot_acc %in% craps)
  }
  
  if (rm_krts) {
    df <- df %>% dplyr::filter(!grepl("^krt[0-9]+", gene, ignore.case = TRUE))
  }
  
  if (annot_kinases) df <- annotKin(df, acc_type)
  if (!any(c("pep_start", "pep_end") %in% names(df))) df <- df %>% annotPeppos(fasta)

  # first to annotate pep_seq to n.sequence.c to avoid the following ambiguity
  # need to know if is n-termal or not
  # M._sequence.c; -._sequence.c; n.sequence.c; -.sequence.c

  # temp <- df %>% filter(pep_seq == "MENGQSTAAK")
  
  browser()

  df <- df %>% 
    # annotPrndesc(fasta) %>% 
    # dplyr::select(-entry_name) %>% 
    dplyr::mutate(pep_seq_mod = pep_seq) %>% 
    dplyr::mutate(pep_seq = toupper(pep_seq))
  
  df <- df %>% 
    dplyr::filter(!duplicated(pep_seq)) %>% 
    dplyr::select(c("pep_seq", "prot_acc")) %>% 
    annotPeppos_mq(fasta) %>% 
    dplyr::select(-prot_acc) %>% 
    dplyr::right_join(., df, by = "pep_seq") %>% 
    dplyr::select(-previous_aa, -next_aa)
  
  
  
  
  
  
  df <- cbind.data.frame(
    df[, grepl("^[a-z]", names(df))], 
    df[, grepl("^[A-Z]", names(df)) & !grepl("^[IR][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[R][0-9]{3}[NC]*", names(df))], 
    df[, grepl("^[I][0-9]{3}[NC]*", names(df))]
  )
  
  if(length(grep("^R[0-9]{3}", names(df))) > 0) {
    df_split <- df %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}|^R[0-9]{3}", names(.)), as.numeric) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x == -1, NA, .x)) %>%
      dplyr::mutate_at(.vars = grep("^I[0-9]{3}", names(.)), ~ ifelse(.x <= rptr_intco, NA, .x)) %>%
      # dplyr::filter(pep_expect <=  0.10) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^R[0-9]{3}", names(.))])) > 0) %>%
      dplyr::filter(rowSums(!is.na(.[grep("^I[0-9]{3}", names(.))])) > 0) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc) %>%
      # a special case of redundant entries from Mascot
      dplyr::filter(!duplicated(.[grep("^pep_seq$|I[0-9]{3}", names(.))]))
  } else {
    df_split <- df %>%
      # dplyr::filter(pep_expect <=  0.10) %>%
      dplyr::arrange(RAW_File, pep_seq, prot_acc)
  }
  
  tmtinj_raw_map <- check_raws(df_split)
  
  df_split <- df_split %>%
    dplyr::left_join(tmtinj_raw_map, id = "RAW_File") %>%
    dplyr::group_by(TMT_inj) %>%
    dplyr::mutate(psm_index = row_number()) %>%
    data.frame(check.names = FALSE) %>%
    split(., .$TMT_inj, drop = TRUE)
  
  missing_tmtinj <- setdiff(names(df_split), unique(tmtinj_raw_map$TMT_inj))
  if (!purrr::is_empty(missing_tmtinj)) {
    cat("The following TMT sets and LC/MS injections do not have corresponindg PSM files:\n")
    cat(paste0("\tTMT.LCMS: ", missing_tmtinj, "\n"))
    
    stop(paste("Remove mismatched `TMT_Set` and/or `LC/MS` under the experimental summary file."),
         call. = FALSE)
  }
  
  fn_lookup <- label_scheme_full %>%
    dplyr::select(TMT_Set, LCMS_Injection, RAW_File) %>%
    dplyr::mutate(filename = paste(paste0("TMTset", .$TMT_Set),
                                   paste0("LCMSinj", .$LCMS_Injection), sep = "_")) %>%
    dplyr::filter(!duplicated(filename)) %>%
    tidyr::unite(TMT_inj, TMT_Set, LCMS_Injection, sep = ".", remove = TRUE) %>% 
    dplyr::select(-RAW_File) %>%
    dplyr::left_join(tmtinj_raw_map, by = "TMT_inj")
  
  
  for (i in seq_along(df_split)) {
    df_split[[i]] <- df_split[[i]] %>% dplyr::select(-TMT_inj)
    
    out_fn <- fn_lookup %>%
      dplyr::filter(TMT_inj == names(df_split)[i]) %>%
      dplyr::select(filename) %>%
      unique() %>%
      unlist() %>%
      paste0(., ".csv")
    
    df_split[[i]] <- df_split[[i]] %>% dplyr::rename(raw_file = RAW_File)
    
    write.csv(df_split[[i]], file.path(dat_dir, "PSM\\cache", out_fn), row.names = FALSE)
    
    if (plot_rptr_int & TMT_plex > 0) {
      TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()
      Levels <- TMT_levels %>% gsub("^TMT-", "I", .)
      df_int <- df_split[[i]][, grepl("^I[0-9]{3}", names(df_split[[i]]))] %>%
        tidyr::gather(key = "Channel", value = "Intensity") %>%
        dplyr::mutate(Channel = factor(Channel, levels = Levels))
      rm(Levels)
      
      mean_int <- df_int %>% 
        dplyr::group_by(Channel) %>% 
        dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
        dplyr::mutate(Intensity = round(Intensity, digit = 1))
      
      Width <- 8
      Height <- 8
      
      filename <- gsub("\\.csv", "\\.png", out_fn)
      p <- ggplot() +
        geom_violin(df_int, mapping = aes(x = Channel, y = log10(Intensity), fill = Channel), size = .25) +
        geom_boxplot(df_int, mapping = aes(x = Channel, y = log10(Intensity)), width = 0.2, lwd = .2, fill = "white") +
        stat_summary(df_int, mapping = aes(x = Channel, y = log10(Intensity)), fun.y = "mean", geom = "point",
                     shape = 23, size = 2, fill = "white", alpha = .5) +
        labs(title = expression("Reporter ions"), x = expression("Channel"), y = expression("Intensity ("*log[10]*")")) +
        geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
        theme_psm_violin
      
      ggsave(file.path(dat_dir, "PSM\\Violin\\bfr_mc", filename), p, width = Width, height = Height, units = "in")
    }
  }
  
}







#' Annotates PSM results
#'
#' \code{annotPSM_sm} adds fields of annotation to PSM tables.
#'
#'@inheritParams load_expts
#'@inheritParams annotPSM
#' @import dplyr tidyr purrr ggplot2 RColorBrewer
#' @importFrom stringr str_split
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
annotPSM_sm <- function(fasta = NULL, expt_smry = "expt_smry.xlsx", rm_krts = FALSE, plot_rptr_int = TRUE) {
  
  old_opt <- options(max.print = 99999)
  on.exit(options(old_opt), add = TRUE)
  
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  
  options(max.print = 5000000)
  
  load(file = file.path(dat_dir, "label_scheme_full.Rdata"))
  load(file = file.path(dat_dir, "label_scheme.Rdata"))
  n_TMT_sets <- n_TMT_sets(label_scheme_full)
  TMT_plex <- TMT_plex(label_scheme_full)
  
  filelist <- list.files(
    path = file.path(dat_dir, "PSM\\cache"),
    pattern = "^TMT.*LCMS.*_Clean.txt$"
  ) %>%
    reorder_files(., n_TMT_sets)
  
  for(set_idx in seq_len(n_TMT_sets)){
    sublist <- filelist[grep(paste0("set*.", set_idx), filelist, ignore.case = TRUE)]
    
    out_fn <- data.frame(Filename =
                           do.call('rbind', strsplit(as.character(sublist),
                                                     '.txt', fixed = TRUE))) %>%
      dplyr::mutate(Filename = gsub("_Clean", "_PSM_N", Filename))
    
    channelInfo <- channelInfo(label_scheme, set_idx)
    
    # LCMS injections under the same TMT experiment
    for (injn_idx in seq_along(sublist)) {
      df <- read.csv(file.path(dat_dir, "PSM\\cache", sublist[injn_idx]),
                     check.names = FALSE, header = TRUE, sep = "\t",
                     comment.char = "#") %>%
        # dplyr::rename(pep_seq_mod = `Modified sequence`) %>%
        dplyr::select(which(names(.) == "pep_seq_mod"),
                      which(names(.) != "pep_seq_mod"))
      
      acc_type <- find_acctype()
      species <- find_df_species(df, acc_type)
      
      label_scheme_full$Species <- species
      label_scheme$Species <- species
      save(label_scheme_full, file = file.path(dat_dir, "label_scheme_full.Rdata"))
      save(label_scheme, file = file.path(dat_dir, "label_scheme.Rdata"))
      write.table(label_scheme[1, c("Accession_Type", "Species")],
                  file.path(dat_dir, "acctype_sp.txt"), sep = "\t",
                  col.names = TRUE, row.names = FALSE)
      
      load_dbs()
      
      if (rm_krts) df <- rm_krts(df, acc_type)
      
      browser()
      browser()
      browser()
      browser()
      browser()
      browser()
      browser()
      browser()
      browser()
      browser()
      
      # nterm == "Acetyl"
      # (5) Add "_" for "Acetyl (Protein N-term)"
      mod <- "Acetyl"
      df_sub <- df[grepl(mod, df$nterm), ]
      
      if (nrow(df_sub) > 0) {
        df_sub$pep_seq_mod <- paste0("_", df_sub$pep_seq_mod)
        df <- rbind(df[!grepl(mod, df$nterm), ], df_sub)
      }
      
      nt_ace_var_mods <- var_mods[grepl("Acetyl (Protein N-term)",
                                        var_mods$Descption, fixed = TRUE), ]
      
      for (mod in nt_ace_var_mods$Mascot_abbr) {
        df_sub <- df[grepl(mod, df$pep_var_mod_pos), ]
        
        if (nrow(df_sub) > 0) {
          # df_sub$pep_seq_mod <- paste0(paste0(substring(df_sub$pep_seq_mod, 1, 1), nt_ace_var_mods$Abbr),
          #                              substring(df_sub$pep_seq_mod, 2)) # insert "_"
          df_sub$pep_seq_mod <- paste0("_", df_sub$pep_seq_mod)
          fn <- nt_ace_var_mods[nt_ace_var_mods$Mascot_abbr == mod, "Filename"]
          try(write.table(df_sub, file.path(dat_dir, "PSM\\individual_mods", paste0(fn, ".txt")), 
                          sep = "\t", col.names = TRUE, row.names = FALSE))
          
          df <- rbind(df[!grepl(mod, df$pep_var_mod_pos), ], df_sub)
        }
      }
      
      rm (nt_ace_var_mods, mod, df_sub, fn)
      
      
      
      df <- df %>% 
        dplyr::mutate(pep_seq_mod = gsub("^_(.*)_$", "\\1", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("^\\(ac\\)", "_", pep_seq_mod)) %>% 
        # dplyr::mutate(pep_seq_mod = gsub("^\\(gl\\)", "q", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("A\\([^\\)]+\\)", "a", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("C\\([^\\)]+\\)", "c", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("D\\([^\\)]+\\)", "d", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("E\\([^\\)]+\\)", "e", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("F\\([^\\)]+\\)", "f", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("G\\([^\\)]+\\)", "g", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("H\\([^\\)]+\\)", "h", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("I\\([^\\)]+\\)", "i", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("K\\([^\\)]+\\)", "k", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("L\\([^\\)]+\\)", "l", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("M\\([^\\)]+\\)", "m", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("N\\([^\\)]+\\)", "n", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("P\\([^\\)]+\\)", "p", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("Q\\([^\\)]+\\)", "q", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("R\\([^\\)]+\\)", "r", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("s\\([^\\)]+\\)", "s", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("T\\([^\\)]+\\)", "t", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("V\\([^\\)]+\\)", "v", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("W\\([^\\)]+\\)", "w", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = gsub("Y\\([^\\)]+\\)", "y", pep_seq_mod)) %>% 
        dplyr::mutate(pep_seq_mod = paste(pep_res_before, pep_seq_mod, pep_res_after, sep = ".")) %>% 
        dplyr::mutate(pep_seq = paste(pep_res_before, pep_seq, pep_res_after, sep = "."))
      
      if(TMT_plex > 0) df <- mcPSM(df, set_idx)
      
      write.table(df, file.path(dat_dir, "PSM", paste0(out_fn[injn_idx, 1], ".txt")),
                  sep = "\t", col.names = TRUE, row.names = FALSE)
      
      if (plot_rptr_int & TMT_plex > 0) {
        TMT_levels <- label_scheme %>% TMT_plex() %>% TMT_levels()
        Levels <- TMT_levels %>% gsub("^TMT-", "I", .)
        df_int <- df[, grepl("^I[0-9]{3}", names(df))] %>%
          tidyr::gather(key = "Channel", value = "Intensity") %>%
          dplyr::mutate(Channel = factor(Channel, levels = Levels)) %>% 
          dplyr::filter(!is.na(Intensity))
        rm(Levels)
        
        mean_int <- df_int %>% 
          dplyr::group_by(Channel) %>% 
          dplyr::summarise(Intensity = mean(log10(Intensity), na.rm = TRUE)) %>% 
          dplyr::mutate(Intensity = round(Intensity, digit = 1))
        
        Width <- 8
        Height <- 8
        filename <- paste0("Post_outlier_removals_set_", set_idx, "_inj_", injn_idx)
        
        p <- ggplot() +
          geom_violin(df_int, mapping=aes(x=Channel, y=log10(Intensity), fill=Channel), size=.25) +
          geom_boxplot(df_int, mapping=aes(x=Channel, y=log10(Intensity)), width=0.2, lwd=.2, fill="white") +
          stat_summary(df_int, mapping=aes(x=Channel, y=log10(Intensity)), fun.y="mean", geom="point",
                       shape=23, size=2, fill="white", alpha=.5) +
          labs(title=expression("Reporter ions"), x=expression("Channel"), y=expression("Intensity ("*log[10]*")")) + 
          geom_text(data = mean_int, aes(x = Channel, label = Intensity, y = Intensity + 0.2), size = 5, colour = "red", alpha = .5) + 
          theme_psm_violin
        
        ggsave(file.path(dat_dir, "PSM\\Violin", paste0(filename, ".png")), p, width = Width, height = Height, units = "in")
      }
      
    }
    
  }
}

