#' Finds the indexes of top-n entries without re-ordering.
#'
#' @param x A numeric vector.
#' @param n The number of top entries to keep.
#' @return The indexes of the top-n entries.
which_topx <- function(x, n = 50L, ...) {
  len <- length(x)
  p <- len - n
  
  if (p  <= 0L) return(seq_along(x))
  
  xp <- sort(x, partial = p, ...)[p]
  
  which(x > xp)
}


#' Finds the top-n entries without re-ordering.
#' 
#' @inheritParams which_topx
#' @return The top-n entries.
topx <- function(x, n = 50L, ...) {
  len <- length(x)
  p <- len - n
  
  if (p  <= 0L) return(x)
  
  xp <- sort(x, partial = p, ...)[p]
  
  x[x > xp]
}


#' Finds the numeric difference in ppm.
#' 
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The difference between \eqn{x} and \eqn{y} in ppm.
find_ppm_error <- function (x = 1000, y = 1000.01) {
  (y - x)/y * 1E6
}


#' Finds the error range of a number.
#'
#' Assumes \eqn{x} is positive without checking.
#' 
#' @param x A numeric value.
#' @param ppm Numeric; the ppm allowed from \code{x}.
#' @return The lower and the upper bound to \eqn{x} by \eqn{ppm}.
find_mass_error_range <- function (x = 500, ppm = 20) {
  d <- x * ppm/1E6
  c(x-d, x+d)
}


#' Splits data by groups then into chunks.
#' 
#' @inheritParams chunksplit
#' @param f A factor; see also base \code{split}.
chunk_groupsplit <- function (data, f, n_chunks) {
  
  if (n_chunks <= 1L) return(data)

  data <- split(data, f)
  len <- length(data)

  labs <- levels(cut(1:len, n_chunks))
  
  x <- cbind(
    lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
    upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))
  
  grps <- findInterval(1:len, x[, 1])
  
  data <- split(data, grps)
  
  map(data, ~ do.call(rbind, .x))
}


#' Splits data into chunks by length.
#' 
#' @param data Input data.
#' @param n_chunks The number of chunks.
#' @param type The type of data for splitting.
#' @export
chunksplit <- function (data, n_chunks = 5L, type = "list") {
  
  stopifnot(type %in% c("list", "row"))
  
  if (n_chunks <= 1L) return(data)
  
  if (type == "list") {
    len <- length(data)
  } else if (type == "row") {
    len <- nrow(data)
  } else {
    stop("Unknown type.", call. = TRUE)
  }
  
  if (len == 0L) return(data)
  
  labs <- levels(cut(1:len, n_chunks))
  
  x <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
             upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))
  
  # x[x[, 1] < 0, 1] <- 1
  # x[x[, 2] > len, 2] <- len
  
  grps <- findInterval(1:len, x[, 1])
  split(data, grps)
}


#' Splits data into chunks with approximately equal sizes.
#'
#' @param nx Positive integer; an arbitrarily large number for data to be split
#'   into for estimating the cumulative sizes.
#' @inheritParams chunksplit
#' @export
chunksplitLB <- function (data, n_chunks = 5L, nx = 100L, type = "list") {
  stopifnot(type %in% c("list", "row"))
  
  if (n_chunks <= 1) return(data)
  
  if (type == "list") {
    len <- length(data)
  } else if (type == "row") {
    len <- nrow(data)
  } else {
    stop("Unknown type.", call. = TRUE)
  }
  
  if (len == 0L) return(data)
  
  # The finer groups by 'nx'
  grps_nx <- local({
    labsx <- levels(cut(1:len, nx))
    
    xx <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1", 
                                              labsx))),
                upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", 
                                                labsx))))
    
    findInterval(1:len, xx[, 1])
  })
  
  # The equated size for a chunk
  size_chunk <- local({
    size_nx <- data %>% 
      split(., grps_nx) %>% 
      purrr::map(object.size) %>% 
      cumsum() 
    
    size_nx[length(size_nx)]/n_chunks
  })
  
  #  Intervals
  grps <- local({
    size_data <- data %>% 
      purrr::map(object.size) %>% 
      cumsum()
    
    # the position indexes
    ps <- purrr::map_dbl(1:(n_chunks-1), ~ {
      which(size_data < size_chunk * .x) %>% `[`(length(.))
    })
    
    grps <- findInterval(1:len, ps)
  })
  
  split(data, grps)
}


#' Searches MS ions.
#'
#' @inheritParams calc_ms2ions
#' @inheritParams calc_pepmasses
#' @inheritParams search_mgf_frames
#' @inheritParams readMGF
#' @inheritParams normPSM
#' @inheritParams calc_pepfdr
#' @inheritParams calc_tmtint
#' @inheritParams load_fasta2
#' @param mgf_path The file path to a list of MGF files. There is no default and
#'   the experimenters need to supply the files.
#' @param fasta Character string(s) to the name(s) of fasta file(s) with
#'   prepended directory path. There is no default and the experimenters need to
#'   supply the files.
#' @param acc_pattern Regular expression(s) describing the patterns to separate
#'   the header lines of fasta entries. At the \code{NULL} default, the pattern
#'   will be automated when \code{acc_type} are among c("uniprot_acc",
#'   "uniprot_id", "refseq_acc", "other"). See also \link{load_fasta2} for
#'   custom examples.
#' @param target_fdr Numeric; the levels of false-discovery rate (FDR) at the
#'   levels of PSMs or peptides. The same level applies further to protein FDR
#'   and indicated under the column \code{prot_issig} in outputs. This results
#'   in FDR controls of either PSM <-> protein or peptide <-> protein at a given
#'   \code{target_fdr}.
#' @param fdr_type Character string; the type of FDR controlling. The value is
#'   in one of c("psm", "peptide"). Note that protein FDR is in conjunction with
#'   "psm" or "peptide". Separate protein FDR may be available in the future.
#' @param combine_tier_three Logical; if TRUE, combines tier-3 proteins in the
#'   output of \code{psmQ.txt}. If FALSE, saves tier-3 proteins separately in
#'   \code{psmT3.txt}.
#' @seealso \link{load_fasta2} for setting the values of \code{acc_type} and
#'   \code{acc_pattern}. \link{parse_unimod} for the grammar of Unimod.
#' @export
matchMS <- function (out_path = "~/proteoQ/outs", 
                     mgf_path = file.path(out_path, "mgf"), 
                     fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
                               "~/proteoQ/dbs/fasta/crap/crap.fasta"), 
                     acc_type = c("uniprot_acc", "other"),
                     acc_pattern = NULL, 
                     fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)"), 
                     varmods = c("TMT6plex (N-term)", "Acetyl (Protein N-term)", 
                                 "Oxidation (M)", "Deamidated (N)", 
                                 "Gln->pyro-Glu (N-term = Q)"), 
                     include_insource_nl = FALSE, 
                     enzyme = c("trypsin"), 
                     maxn_fasta_seqs = 200000L,
                     maxn_vmods_setscombi = 64L, 
                     maxn_vmods_per_pep = 5L,
                     maxn_sites_per_vmod = 3L, 
                     maxn_vmods_sitescombi_per_pep = 64L, 
                     min_len = 7L, max_len = 100L, max_miss = 2L, 
                     type_ms2ions = "by", 
                     topn_ms2ions = 100L, 
                     minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                     ppm_reporters = 10L, 
                     quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"), 
                     target_fdr = 0.01, 
                     fdr_type = c("psm", "peptide", "protein"), 
                     combine_tier_three = FALSE, 
                     digits = 4L) {
  
  options(digits = 9L)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), 
                   fun = as.character(match.call()[[1]]))
      }
    }, 
    add = TRUE
  )
  
  out_path <- create_dir(out_path)
  
  local({
    filelist <- list.files(path = file.path(mgf_path), pattern = "\\.mgf$")
    
    if (purrr::is_empty(filelist)) {
      stop("No `.mgf` files under ", mgf_path, 
           call. = FALSE)
    }
  })
  
  ## Indexes of modifications
  mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
    as.hexmode() %>% 
    `names<-`(c(fixedmods, varmods))
  
  ## Theoretical MS1 masses
  res <- calc_pepmasses2(
    fasta = fasta, 
    acc_type = acc_type,
    acc_pattern = acc_pattern, 
    fixedmods = fixedmods, 
    varmods = varmods, 
    include_insource_nl = include_insource_nl, 
    index_mods = FALSE, 
    enzyme = enzyme, 
    maxn_fasta_seqs = maxn_fasta_seqs, 
    maxn_vmods_setscombi = maxn_vmods_setscombi, 
    maxn_vmods_per_pep = maxn_vmods_per_pep, 
    maxn_sites_per_vmod = maxn_sites_per_vmod, 
    min_len = min_len, 
    max_len = max_len, 
    max_miss = max_miss, 
    out_path = out_path, 
    digits = digits, 
    parallel = TRUE
  )
  
  ## Hash tables
  
  
  
  
  ## AA masses
  aa_masses_all <- res$fwd %>% 
    purrr::map(~ {
      attr(.x, "data") <- NULL
      .x
    })
  
  ## Mass range
  min_mass <- 500L
  
  max_mass <- res$fwd %>% 
    purrr::map(attr, "data") %>% 
    unlist(use.names = FALSE) %>% 
    max(na.rm = TRUE)

  ## Bin theoretical peptides
  bin_ms1masses(res, min_mass, max_mass, ppm_ms1, out_path)
  rm(list = c("res"))
  
  ## MGFs
  mgf_frames <- load_mgfs(mgf_path = mgf_path, 
                          min_mass = min_mass, 
                          min_ms2mass = 110L, 
                          topn_ms2ions = topn_ms2ions, 
                          ppm_ms1 = ppm_ms1, 
                          ppm_ms2 = ppm_ms2, 
                          index_ms2 = FALSE)
  
  gc()
  
  ## MSMS matches
  out <- ms2match(mgf_path = mgf_path, 
                  aa_masses_all = aa_masses_all, 
                  out_path = out_path, 
                  mod_indexes = mod_indexes, 
                  type_ms2ions = type_ms2ions, 
                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                  maxn_vmods_sitescombi_per_pep = 
                    maxn_vmods_sitescombi_per_pep, 
                  minn_ms2 = minn_ms2, 
                  ppm_ms1 = ppm_ms1, 
                  ppm_ms2 = ppm_ms2, 
                  min_ms2mass = 110L, 
                  
                  # dummy for argument matching
                  fasta = fasta, 
                  acc_type = acc_type, 
                  acc_pattern = acc_pattern,
                  topn_ms2ions = topn_ms2ions,
                  fixedmods = fixedmods, 
                  varmods = varmods, 
                  include_insource_nl = include_insource_nl, 
                  enzyme = enzyme,
                  maxn_fasta_seqs = maxn_fasta_seqs,
                  maxn_vmods_setscombi = maxn_vmods_setscombi, 
                  min_len = min_len, 
                  max_len = max_len, 
                  max_miss = max_miss, 
                  # target_fdr = target_fdr, 
                  # fdr_type = fdr_type, 
                  # quant = !!enexpr(quant), 
                  # ppm_reporters = ppm_reporters, 
                  
                  digits = digits)
  
  # Peptide scores
  fdr_type <- rlang::enexpr(fdr_type)
  oks <- eval(formals()[["fdr_type"]])
  
  if (length(fdr_type) > 1L) {
    fdr_type <- oks[[1]]
  } else {
    fdr_type <- rlang::as_string(fdr_type)
  }
  
  stopifnot(fdr_type %in% oks)
  rm(list = c("oks"))
  
  # TMT intensities
  quant <- rlang::enexpr(quant)
  oks <- eval(formals()[["quant"]])
  
  if (length(quant) > 1L) {
    quant <- oks[[1]]
  } else {
    quant <- rlang::as_string(quant)
  }
  
  stopifnot(quant %in% oks, length(quant) == 1L)
  rm(list = c("oks"))
  
  reporters <- out %>% 
    imap( ~ {
      .x %>% 
        calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
        tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                     remove = TRUE) %>% 
        dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) 
    }) %>% 
    dplyr::bind_rows()

  # ---
  out <- calc_pepscores(topn_ms2ions = topn_ms2ions, 
                        type_ms2ions = type_ms2ions, 
                        target_fdr = target_fdr, 
                        fdr_type = fdr_type, 
                        min_len = min_len, 
                        max_len = max_len, 
                        penalize_sions = TRUE, 
                        ppm_ms2 = ppm_ms2, 
                        out_path = out_path, 
                        digits = digits)
  
  gc()

  # Protein accessions, score cut-offs
  out <- out %>% add_prot_acc() 
  out <- out %>% calc_protfdr(target_fdr)

  gc()
  
  # Clean-ups
  data.frame(Abbr = as.character(mod_indexes), Desc = names(mod_indexes)) %>% 
    readr::write_tsv(file.path(out_path, "mod_indexes.txt"))
  
  # Brings back reporter intensities
  if (quant != "none") {
    out <- out %>% 
      tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                   remove = FALSE) %>% 
      dplyr::left_join(reporters, by = "uniq_id") %>% 
      dplyr::select(-uniq_id)
  }
  
  out <- out %>% 
    dplyr::mutate(pep_ms1_delta = ms1_mass - theo_ms1) %>% 
    dplyr::rename(pep_scan_title = scan_title, 
                  pep_exp_mz = ms1_moverz, 
                  pep_exp_mr = ms1_mass, 
                  pep_exp_z = ms1_charge, 
                  pep_calc_mr = theo_ms1, 
                  pep_delta = pep_ms1_delta, 
                  
                  pep_tot_int = ms1_int, 
                  pep_ret_time = ret_time, 
                  pep_scan_num = scan_num, 
                  pep_ms2_n = ms2_n, 
                  pep_frame = frame)
  
  out <- dplyr::bind_cols(
    out %>% .[grepl("^prot_", names(.))], 
    out %>% .[grepl("^pep_", names(.))], 
    out %>% .[grepl("^psm_", names(.))], 
    out %>% .[!grepl("^prot_|^pep_|^psm_", names(.))], 
  ) %>% 
    reloc_col_after("pep_exp_z", "pep_exp_mr") %>% 
    reloc_col_after("pep_calc_mr", "pep_exp_z") %>% 
    reloc_col_after("pep_delta", "pep_calc_mr") 
  
  # ---
  out <- out %T>% 
    readr::write_tsv(file.path(out_path, "psmC.txt")) %>% 
    dplyr::filter(pep_issig, !pep_isdecoy, pep_rank <= 3L, 
                  !grepl("^-", prot_acc))
  
  # set aside one-hit wonders
  out0 <- out %>% 
    dplyr::filter(!prot_issig, prot_n_pep == 1L) %>% 
    dplyr::mutate(prot_tier = 3L)
  
  out <- dplyr::bind_rows(
    out %>% dplyr::filter(prot_issig), 
    out %>% filter(!prot_issig, prot_n_pep >= 2L)
  ) %>% 
    dplyr::mutate(prot_tier = ifelse(prot_issig, 1L, 2L))
  
  gc()
  
  # Protein groups
  out <- out %>% grp_prots(file.path(out_path, "temp1")) 
  out0 <- out0 %>% grp_prots(file.path(out_path, "temp0")) 
  
  out <- dplyr::bind_cols(
    out %>% .[grepl("^prot_", names(.))], 
    out %>% .[grepl("^pep_", names(.))], 
    out %>% .[grepl("^psm_", names(.))], 
    out %>% .[!grepl("^prot_|^pep_|^psm_", names(.))], 
  ) %>% 
    reloc_col_after("prot_es", "prot_family_member") %>% 
    reloc_col_after("prot_es_co", "prot_es") %>% 
    reloc_col_after("prot_tier", "prot_isess")
  
  out <- local({
    message("Tier-3 proteins are low-confidence identifications.")
    
    max <- max(out$prot_hit_num, na.rm = TRUE)
    
    if (combine_tier_three) {
      out <- out %>% 
        dplyr::bind_rows(out0) %>% 
        dplyr::arrange(prot_acc, pep_seq) %T>% 
        readr::write_tsv(file.path(out_path, "psmQ.txt"))
    } else {
      message("Tier-3 proteins saved separately in `psmT3.txt`.")
      
      out0 <- out0[names(out)] %>% 
        dplyr::mutate(prot_hit_num = prot_hit_num + max)  %T>% 
        readr::write_tsv(file.path(out_path, "psmT3.txt"))
      
      out <- out %>% 
        dplyr::arrange(prot_acc, pep_seq) %T>% 
        readr::write_tsv(file.path(out_path, "psmQ.txt"))
    }
    
    invisible(out)
  })
  
  message("Search completed.")
  
  # ---
  rm(list = c(".path_cache", ".path_fasta", ".time_stamp"), envir = .GlobalEnv)
  
  .savecall <- TRUE
  
  invisible(out)
}



#' Subsets the frames of theoretical peptides.
#' 
#' @inheritParams search_mgf_frames_d
subset_theoframes <- function (mgf_frames, theopeps) {
  
  if (purrr::is_empty(mgf_frames) || purrr::is_empty(theopeps)) {
    return(NULL)
  }
  
  frames <- as.integer(names(mgf_frames))
  breaks <- which(diff(frames) != 1) + 1
  grps <- findInterval(frames, frames[breaks])
  frames <- split(frames, grps)
  
  frames <- frames %>% 
    purrr::map(~ c(.x[1]-1, .x, .x[length(.x)]+1)) %>% 
    unlist(use.names = FALSE) %>% 
    .[!duplicated(.)]
  
  theopeps[as.character(frames)]
}


#' Matches theoretical peptides (parallel by mgf chunks).
#'
#' All files under `out_path` are removed if incur \code{calc_pepmasses} in the
#' upstream.
#'
#' @param aa_masses_all A list of amino acid lookups for all the combination of
#'   fixed and variable modifications.
#' @param n_cores Integer; the number of computer cores.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @inheritParams matchMS
#' @inheritParams readMGF
#' @inheritParams search_mgf_frames_d
#' @import parallel
pmatch_bymgfs <- function (mgf_path, aa_masses_all, n_cores, out_path, 
                           mod_indexes, type_ms2ions, maxn_vmods_per_pep, 
                           maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                           minn_ms2, ppm_ms1, ppm_ms2, 
                           
                           # dummy
                           fasta, acc_type, acc_pattern,
                           topn_ms2ions, fixedmods, varmods, 
                           include_insource_nl, enzyme, 
                           maxn_fasta_seqs, maxn_vmods_setscombi, 
                           min_len, max_len, max_miss, 
                           # topn_ms2ions, target_fdr, fdr_type, quant, ppm_reporters, 
                           
                           digits) {

  on.exit(
    if (exists(".savecall", envir = current_env())) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = "pmatch_bymgfs")
      }
    }, 
    add = TRUE
  )
  
  # Check cached 
  fun <- as.character(match.call()[[1]])
  args_except <- c("n_cores")
  
  cache_pars <- find_callarg_vals(time = NULL, 
                                  path = file.path(out_path, "Calls"), 
                                  fun = paste0(fun, ".rda"), 
                                  arg = names(formals(fun))) %>% 
    .[! . %in% args_except] %>% 
    .[sort(names(.))]
  
  call_pars <- mget(names(formals()), envir = rlang::current_env(), 
                    inherits = FALSE) %>% 
    .[! . %in% args_except] %>% 
    .[sort(names(.))]

  if (identical(cache_pars, call_pars)) {
    files <- c(
      file.path(out_path, "ion_matches.rds"), 
      file.path(out_path, "ion_matches_rev.rds")
    )
    
    if (all(file.exists(files))) {
      message("Loading ion matches from cache.")
      
      out <- readRDS(file.path(out_path, "ion_matches.rds"))
      
      .savecall <- FALSE
      
      return(out)
    }
  }
  
  rm(list = c("fun", "args_except", "cache_pars", "call_pars"))
  
  delete_files(out_path, ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                                     "\\.xls$", "\\.csv$", "\\.txt$", 
                                     "^mgf$", "^mgfs$"))

  # Targets 
  out <- vector("list", length(aa_masses_all))
  
  out <- map2(seq_along(aa_masses_all), aa_masses_all, 
              pmatch_bymgfs_i, 
              mgf_path = mgf_path, 
              n_cores = n_cores, 
              out_path = out_path, 
              mod_indexes = mod_indexes, 
              type_ms2ions = type_ms2ions, 
              maxn_vmods_per_pep = maxn_vmods_per_pep, 
              maxn_sites_per_vmod = maxn_sites_per_vmod, 
              maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
              minn_ms2 = minn_ms2, 
              ppm_ms1 = ppm_ms1, 
              ppm_ms2 = ppm_ms2, 
              digits = digits) %>% 
    `names<-`(seq_along(.)) %T>% 
    saveRDS(file.path(out_path, "ion_matches.rds")) 

  # Decoys
  maxs <- map_dbl(out, object.size) %>% which_topx(1)
  maxs2 <- maxs %>% paste0("rev_", .)

  rev  <- map2(maxs2, aa_masses_all[maxs], 
               pmatch_bymgfs_i, 
               mgf_path = mgf_path, 
               n_cores = n_cores, 
               out_path = out_path, 
               mod_indexes = mod_indexes, 
               type_ms2ions = type_ms2ions, 
               maxn_vmods_per_pep = maxn_vmods_per_pep, 
               maxn_sites_per_vmod = maxn_sites_per_vmod, 
               maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
               minn_ms2 = minn_ms2, 
               ppm_ms1 = ppm_ms1, 
               ppm_ms2 = ppm_ms2, 
               digits = digits) %>% 
    `names<-`(maxs2) %T>% 
    saveRDS(file.path(out_path, "ion_matches_rev.rds")) 
  
  .savecall <- TRUE

  invisible(out)
}


#' Helper of \link{pmatch_bymgfs}.
#'
#' @param i Integer; the index for a set of corresponding aa_masses and
#'   theoretical peptides.
#' @param is_target Logical; indicator of target or decoy peptides for
#'   searching.
#' @inheritParams pmatch_bymgfs
pmatch_bymgfs_i <- function (i, aa_masses, mgf_path, n_cores, out_path, 
                             mod_indexes, type_ms2ions, maxn_vmods_per_pep, 
                             maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                             minn_ms2, ppm_ms1, ppm_ms2, 
                             digits) {
  
  # loads freshly mgfs (as will be modified)
  mgf_frames <- readRDS(file.path(mgf_path, "mgf_queries.rds")) %>% 
    dplyr::group_by(frame) %>% 
    dplyr::group_split() %>% 
    setNames(purrr::map_dbl(., ~ .x$frame[1]))
  
  # x = mgf_frames %>% .[names() == 80964]
  # which(names(x)) == 56090
  
  mgf_frames <- local({
    labs <- levels(cut(1:length(mgf_frames), n_cores^2))
    
    x <- cbind(
      lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
      upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))
    
    grps <- findInterval(1:length(mgf_frames), x[, 1])
    
    split(mgf_frames, grps)
  })
  
  # parses aa_masses
  nm_fmods <- attr(aa_masses, "fmods", exact = TRUE)
  nm_vmods <- attr(aa_masses, "vmods", exact = TRUE)
  msg_end <- if (grepl("^rev_", i)) " (decoy)." else "." 

  message("Matching against: ", 
          paste0(nm_fmods, 
                 nm_vmods %>% { if (nchar(.) > 0) paste0(" + ", .) else . }, 
                 msg_end))

  # reads theoretical peptide data
  .path_fasta <- get(".path_fasta", envir = .GlobalEnv)
  .time_stamp <- get(".time_stamp", envir = .GlobalEnv)
  
  theopeps <- readRDS(file.path(.path_fasta, "pepmasses", .time_stamp, 
                                paste0("binned_theopeps_", i, ".rds")))
  
  if (ppm_ms1 < 1000) {
    theopeps <- theopeps %>% 
      purrr::map(~ {
        rows <- !duplicated(.x$pep_seq)
        .x[rows, c("pep_seq", "mass")]
      })
  }
  
  # (1) for a given aa_masses_all[[i]], some mgf_frames[[i]] 
  #     may not be found in theopeps[[i]] 
  mgf_frames <- purrr::map(mgf_frames, ~ {
    x <- .x
    
    oks <- names(x) %in% names(theopeps)
    x <- x[oks]
    
    empties <- purrr::map_lgl(x, purrr::is_empty)
    x[!empties]
  })
  
  # (2) splits `theopeps` in accordance to `mgf_frames` with 
  #     preceding and following frames: (o)|range of mgf_frames[[1]]|(o)
  theopeps <- local({
    frames <- purrr::map(mgf_frames, ~ as.integer(names(.x)))
    
    mins <- purrr::map_dbl(frames, ~ {
      if (length(.x) == 0) x <- 0 else x <- min(.x, na.rm = TRUE)
    })
    
    maxs <- purrr::map_dbl(frames, ~ {
      if (length(.x) == 0) x <- 0 else x <- max(.x, na.rm = TRUE)
    })
    
    nms <- as.integer(names(theopeps))
    
    theopeps <- purrr::map2(mins, maxs, ~ {
      theopeps[which(nms >= (.x - 1) & nms <= (.y + 1))]
    })
    
    # --- may cause uneven length between `mgf_frames` and `theopeps`
    # empties <- map_lgl(theopeps, is_empty)
    # theopeps[!empties]
  })
  
  # (3) removes unused frames of `theopeps`
  theopeps <- purrr::map2(mgf_frames, theopeps, subset_theoframes)
  
  # (4) removes empties (zero overlap between mgf_frames and theopeps)
  oks <- purrr::map_lgl(mgf_frames, ~ !purrr::is_empty(.x)) | 
    purrr::map_lgl(theopeps, ~ !purrr::is_empty(.x))
  
  mgf_frames <- mgf_frames[oks]
  theopeps <- theopeps[oks]
  
  rm(list = "oks")
  
  if (is_empty(mgf_frames) || is_empty(theopeps)) {
    return(NULL)
  }
  
  cl <- makeCluster(getOption("cl.cores", n_cores))
  
  clusterExport(cl, list("%>%"), 
                envir = environment(magrittr::`%>%`))
  clusterExport(cl, list("search_mgf_frames_d"), 
                envir = environment(proteoQ:::search_mgf_frames_d))
  clusterExport(cl, list("search_mgf_frames"), 
                envir = environment(proteoQ:::search_mgf_frames))
  clusterExport(cl, list("search_mgf"), 
                envir = environment(proteoQ:::search_mgf))
  clusterExport(cl, list("find_ppm_outer_bypep"), 
                envir = environment(proteoQ:::find_ppm_outer_bypep))
  clusterExport(cl, list("find_ppm_outer_bycombi"), 
                envir = environment(proteoQ:::find_ppm_outer_bycombi))
  
  out <- clusterMap(cl, search_mgf_frames_d, 
                    mgf_frames, theopeps, 
                    MoreArgs = list(aa_masses = aa_masses, 
                                    mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = 
                                      maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = 
                                      maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    minn_ms2 = minn_ms2, 
                                    ppm_ms1 = ppm_ms1, 
                                    ppm_ms2 = ppm_ms2, 
                                    digits = digits), 
                    .scheduling = "dynamic") %>% 
    dplyr::bind_rows() %>% # across nodes
    dplyr::mutate(pep_fmod = nm_fmods, 
                  pep_vmod = nm_vmods, 
                  pep_mod_group = as.character(i)) %>% 
    { if (msg_end == ".") dplyr::mutate(., pep_isdecoy = FALSE) else 
      dplyr::mutate(., pep_isdecoy = TRUE) }
  
  out <- out %>%
    dplyr::mutate(raw_file = gsub("\\\\", "/", scan_title)) %>% 
    dplyr::mutate(raw_file = gsub("^.*/(.*)\\.(raw|RAW)[\\\"]{0,1}; .*", "\\1", 
                                  raw_file)) %>% 
    dplyr::mutate(raw_file = gsub("^.*/(.*)\\.d[\\\"]{0,1}; .*", "\\1", 
                                  raw_file)) %>% 
    reloc_col_after("raw_file", "scan_num") %>% 
    reloc_col_after("pep_mod_group", "raw_file")

  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, file.path(out_path, "temp", paste0("ion_matches_", i, ".rds")))

  stopCluster(cl)
  
  gc()
  
  invisible(out)
}



#' Helper of \link{search_mgf_frames}
#'
#' Searches MGFs in a frame at a given combination of fixed and variable
#' modifications.
#'
#' @param theopeps Binned theoretical peptides at a given combination of fixed
#'   and variable.
#' @param aa_masses Amino-acid lookup at a given combination of fixed and
#'   variable.
#' @inheritParams mcalc_monopep
#' @inheritParams search_mgf_frames
#' @export
search_mgf_frames_d <- function (mgf_frames, theopeps, aa_masses, 
                                 mod_indexes, type_ms2ions = "by", 
                                 maxn_vmods_per_pep = 5, 
                                 maxn_sites_per_vmod = 3, 
                                 maxn_vmods_sitescombi_per_pep = 32, 
                                 minn_ms2 = 7, 
                                 ppm_ms1 = 20, ppm_ms2 = 25, digits = 5) {
  
  # `res[[i]]` contains results for multiple mgfs within a frame
  # (the number of entries equals to the number of mgf frames)
  res <- search_mgf_frames(mgf_frames = mgf_frames, 
                           theopeps = theopeps, 
                           aa_masses = aa_masses, 
                           mod_indexes = mod_indexes, 
                           type_ms2ions = type_ms2ions, 
                           maxn_vmods_per_pep = maxn_vmods_per_pep, 
                           maxn_sites_per_vmod = maxn_sites_per_vmod, 
                           maxn_vmods_sitescombi_per_pep = 
                             maxn_vmods_sitescombi_per_pep, 
                           minn_ms2 = minn_ms2, 
                           ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
                           digits = digits) 

  # flatten mgfs within each frame
  # (the number of entries equals to the number of mgfs)
  res <- res %>% unlist(recursive = FALSE)
  
  empties <- purrr::map_lgl(res, purrr::is_empty)
  
  # !!!dplyr::bind_rows() temporarily not working for list-columns!!!
  # (contains attr(,"theo_ms1") for theoretical MS1 mass)
  res <- do.call(rbind, mgf_frames) %>% 
    dplyr::mutate(matches = res) 
  
  # res <- mgf_frames %>% 
  #   dplyr::bind_rows() %>% 
  #   dplyr::mutate(matches = res) 
  
  res <- res[!empties, ]
  
  rm(list = "mgf_frames", "theopeps")
  
  invisible(res)
}


#' Searches MGFs in a frame.
#'
#' It reads and searches one frame of MGFs against \code{theopeps}. The frame
#' number links experimental and theoretical spectra by MS1 Masses.
#'
#' @param theopeps Binned theoretical peptides corresponding to an i-th
#'   \code{aa_masses}.
#' @param mgf_frames MGFs in frames. Each frame contains one to multiple MGFs
#'   whose MS1 masses are in the same interval.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @inheritParams mcalc_monopep
#' @inheritParams calc_ms2ions
#' @return Matches to each MGF as a list elements. The length of the output is
#'   equal to the number of MGFs in the given frame.
#' @export
search_mgf_frames <- function (mgf_frames, theopeps, aa_masses, mod_indexes, 
                               type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5, 
                               maxn_sites_per_vmod = 3, 
                               maxn_vmods_sitescombi_per_pep = 32, 
                               minn_ms2 = 7, 
                               ppm_ms1 = 20, ppm_ms2 = 25, digits = 5) {
                                
  len <- length(mgf_frames)
  out <- vector("list", len) 
  
  ## --- initiation ---
  mgfs_cr <- mgf_frames[[1]]
  frame <- mgfs_cr[["frame"]][[1]]
  
  theos_bf_ms1 <- theopeps[[as.character(frame-1)]]
  theopeps_bf_ms1 <- theos_bf_ms1$pep_seq
  theomasses_bf_ms1 <- theos_bf_ms1$mass

  theos_cr_ms1 <- theopeps[[as.character(frame)]]
  theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
  theomasses_cr_ms1 <- theos_cr_ms1$mass

  theos_bf_ms2 <- purrr::map2(theopeps_bf_ms1, 
                              theomasses_bf_ms1, 
                              # calc_ms2ions, 
                              calc_ms2ions2, 
                              aa_masses = aa_masses, 
                              mod_indexes = mod_indexes, 
                              type_ms2ions = type_ms2ions, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              maxn_vmods_sitescombi_per_pep = 
                                maxn_vmods_sitescombi_per_pep, 
                              digits = digits) %>% 
    `names<-`(theopeps_bf_ms1)
  
  theos_cr_ms2 <- purrr::map2(theopeps_cr_ms1, 
                              theomasses_cr_ms1, 
                              # calc_ms2ions, 
                              calc_ms2ions2, 
                              aa_masses = aa_masses, 
                              mod_indexes = mod_indexes, 
                              type_ms2ions = type_ms2ions, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              maxn_vmods_sitescombi_per_pep = 
                                maxn_vmods_sitescombi_per_pep, 
                              digits = digits) %>% 
    `names<-`(theopeps_cr_ms1)
  
  ## --- iteration ---
  for (i in seq_len(len)) {
    exptmasses_ms1 <- mgfs_cr[["ms1_mass"]]
    exptmoverzs_ms2 <- mgfs_cr[["ms2_moverz"]]
    
    theos_af_ms1 <- theopeps[[as.character(frame+1)]]
    theopeps_af_ms1 <- theos_af_ms1$pep_seq
    theomasses_af_ms1 <- theos_af_ms1$mass

    theos_af_ms2 <- purrr::map2(theopeps_af_ms1, 
                                theomasses_af_ms1, 
                                # calc_ms2ions, 
                                calc_ms2ions2, 
                                aa_masses = aa_masses, 
                                mod_indexes = mod_indexes, 
                                type_ms2ions = type_ms2ions, 
                                maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                maxn_vmods_sitescombi_per_pep = 
                                  maxn_vmods_sitescombi_per_pep, 
                                digits = digits) %>% 
      `names<-`(theopeps_af_ms1)
    
    # each `out` for the results of multiple mgfs in one frame
    
    # Browse[4]> exptmasses_ms1
    # [[1]]
    # [1] 748.426367
    
    # [[2]]
    # [1] 748.427407
    
    # Browse[4]> out[[i]]
    # [[1]]
    # named list()
    
    # [[2]]
    # named list()
    
    out[[i]] <- purrr::map2(exptmasses_ms1, exptmoverzs_ms2, 
                            search_mgf, 
                            theomasses_bf_ms1, 
                            theomasses_cr_ms1, 
                            theomasses_af_ms1, 
                            theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                            minn_ms2, ppm_ms1, ppm_ms2) 
    
    # advance to the next frame
    if (i == len) {
      break
    }
    
    mgfs_cr <- mgf_frames[[i+1]]
    new_frame <- mgfs_cr[["frame"]][[1]]
    
    if (isTRUE(new_frame == (frame+1))) {
      theos_bf_ms1 <- theos_cr_ms1
      # theopeps_bf_ms1 <- theopeps_cr_ms1 
      theomasses_bf_ms1 <- theomasses_cr_ms1
      theos_bf_ms2 <- theos_cr_ms2

      theos_cr_ms1 <- theos_af_ms1
      # theopeps_cr_ms1 <- theopeps_af_ms1 
      theomasses_cr_ms1 <- theomasses_af_ms1
      theos_cr_ms2 <- theos_af_ms2
    } else if (isTRUE(new_frame == (frame+2))) {
      theos_bf_ms1 <- theos_af_ms1
      # theopeps_bf_ms1 <- theopeps_af_ms1 
      theomasses_bf_ms1 <- theomasses_af_ms1
      theos_bf_ms2 <- theos_af_ms2

      theos_cr_ms1 <- theopeps[[as.character(new_frame)]]
      theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
      theomasses_cr_ms1 <- theos_cr_ms1$mass

      theos_cr_ms2 <- purrr::map2(theopeps_cr_ms1, 
                                  theomasses_cr_ms1, 
                                  # calc_ms2ions, 
                                  calc_ms2ions2, 
                                  aa_masses = aa_masses, 
                                  mod_indexes = mod_indexes, 
                                  type_ms2ions = type_ms2ions, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  maxn_vmods_sitescombi_per_pep = 
                                    maxn_vmods_sitescombi_per_pep, 
                                  digits = digits) %>% 
        `names<-`(theopeps_cr_ms1)
    } else {
      theos_bf_ms1 <- theopeps[[as.character(new_frame-1)]]
      theopeps_bf_ms1 <- theos_bf_ms1$pep_seq
      theomasses_bf_ms1 <- theos_bf_ms1$mass

      theos_cr_ms1 <- theopeps[[as.character(new_frame)]]
      theopeps_cr_ms1 <- theos_cr_ms1$pep_seq
      theomasses_cr_ms1 <- theos_cr_ms1$mass

      theos_bf_ms2 <- purrr::map2(theopeps_bf_ms1, 
                                  theomasses_bf_ms1, 
                                  # calc_ms2ions, 
                                  calc_ms2ions2, 
                                  aa_masses = aa_masses, 
                                  mod_indexes = mod_indexes, 
                                  type_ms2ions = type_ms2ions, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  maxn_vmods_sitescombi_per_pep = 
                                    maxn_vmods_sitescombi_per_pep, 
                                  digits = digits) %>% 
        `names<-`(theopeps_bf_ms1)
      
      theos_cr_ms2 <- purrr::map2(theopeps_cr_ms1, 
                                  theomasses_cr_ms1, 
                                  # calc_ms2ions, 
                                  calc_ms2ions2, 
                                  aa_masses = aa_masses, 
                                  mod_indexes = mod_indexes, 
                                  type_ms2ions = type_ms2ions, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  maxn_vmods_sitescombi_per_pep = 
                                    maxn_vmods_sitescombi_per_pep, 
                                  digits = digits) %>% 
        `names<-`(theopeps_cr_ms1)
    }
    
    frame <- new_frame
  }
  
  rm(list = c("mgf_frames", "theopeps", "theos_bf_ms1", "theos_cr_ms1", 
              "theos_af_ms1", "theomasses_bf_ms1", "theomasses_cr_ms1", 
              "theomasses_af_ms1", "theopeps_bf_ms1", "theopeps_cr_ms1", 
              "theopeps_af_ms1", "theos_bf_ms2", "theos_cr_ms2", 
              "theos_af_ms2", "exptmasses_ms1", "exptmoverzs_ms2", 
              "mgfs_cr", "new_frame", "frame"))
  
  invisible(out)
}


#' Searches a single MGF.
#'
#' @param expt_mass_ms1 Numeric; the experimental MS1 mass.
#' @param expt_moverz_ms2 A numeric list; the experimental MS2 m/z's.
#' @param theomasses_bf_ms1 Numeric vector; the theoretical MS1 masses at the
#'   preceding \code{-1} frame.
#' @param theomasses_cr_ms1 Numeric vector; the theoretical MS1 masses at the
#'   current frame.
#' @param theomasses_af_ms1 Numeric vector; the theoretical MS1 masses at the
#'   following \code{+1} frame.
#' @param theos_bf_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   preceding \code{-1} frame.
#' @param theos_cr_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   current frame.
#' @param theos_af_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   following \code{+1} frame.
#' @inheritParams search_mgf_frames
#' @import dplyr
#' @import purrr
#' @return Lists of tibbles.
#' @export 
search_mgf <- function (expt_mass_ms1, expt_moverz_ms2, 
                        theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
                        theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                        minn_ms2 = 7, ppm_ms1 = 20, ppm_ms2 = 25) {

  # --- subsets from the `before` and the `after` by MS1 mass tolerance 
  mass_ranges <- find_mass_error_range(expt_mass_ms1, ppm_ms1)
  bf_allowed <- which(theomasses_bf_ms1 >= mass_ranges[1])
  af_allowed <- which(theomasses_af_ms1 <= mass_ranges[2])
  
  # not used but kept for tidiness
  theomasses_bf_ms1 <- theomasses_bf_ms1[bf_allowed]
  theomasses_af_ms1 <- theomasses_af_ms1[af_allowed]

  theos_bf_ms2 <- theos_bf_ms2[bf_allowed]
  theos_af_ms2 <- theos_af_ms2[af_allowed]

  # --- find MS2 matches ---
  if (is_empty(theos_bf_ms2)) {
    x_bf <- theos_bf_ms2
  } else {
    x_bf <- purrr::map(theos_bf_ms2, find_ppm_outer_bypep, 
                       expt_moverz_ms2, ppm_ms2)
  }
  
  if (is_empty(theos_cr_ms2)) {
    x_cr <- theos_cr_ms2
  } else {
    x_cr <- purrr::map(theos_cr_ms2, find_ppm_outer_bypep, 
                       expt_moverz_ms2, ppm_ms2)
  }
  
  if (is_empty(theos_af_ms2)) {
    x_af <- theos_af_ms2
  } else {
    x_af <- purrr::map(theos_af_ms2, find_ppm_outer_bypep, 
                       expt_moverz_ms2, ppm_ms2)
  }
  
  x <- c(x_bf, x_cr, x_af)
  
  # cleans up
  rows <- map(x, ~ {
    this <- .x
    map_lgl (this, ~ sum(!is.na(.x[["expt"]])) >= minn_ms2)
  })
  x <- map2(x, rows, ~ .x[.y])
  
  empties <- map_lgl(x, is_empty)
  x <- x[!empties]
  
  # ---
  theomasses_ms1 <- c(theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1)
  theomasses_ms1 <- theomasses_ms1[!empties]
  x <- map2(x, theomasses_ms1, ~ {
    attr(.x, "theo_ms1") <- .y
    .x
  })
  
  # ---
  # length(x) == N(theos_peps) within the ppm window
  # 
  # ATIPIFFDMMLCEYQR
  # (1) ATIPIFFDMMLCEYQR$`0000000050000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #     1  173.  173.
  #     2  175.  175.
  # (2) ATIPIFFDMMLCEYQR$`0000000005000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #     1  173.  173.
  #     2  175.  175.
  # $KADEQMESMTYSTER
  # ...
  
  ## No evidence of M
  # 
  # $ATIPIFFDMMLCEYQR
  # $ATIPIFFDMMLCEYQR$`0000000050000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #   1  173.  173.
  #   2  175.  175.
  #   3  643.  643.
  #   4  790.  790.
  #   5  868.  868.
  #   6 1297. 1297.
  # 
  # $ATIPIFFDMMLCEYQR$`0000000005000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #   1  173.  173.
  #   2  175.  175.
  #   3  643.  643.
  #   4  790.  790.
  #   5  868.  868.
  #   6 1297. 1297.
  
  invisible(x)
}


#' Helper: finds the the outer products for vectors of MS2 ions.
#'
#' The same theoretical peptides at different position permutations.
#'
#' @param expts Numeric vector; one series of experimental MS2s.
#' @param theos Numeric vector; one to multiple series of theoretical MS2s.
#' @importFrom purrr map
#' @inheritParams search_mgf_frames
find_ppm_outer_bypep <- function (theos, expts, ppm_ms2) {
  
  if (!is.list(theos)) {
    theos <- list(theos)
  }
  
  map(theos, find_ppm_outer_bycombi, expts, ppm_ms2)
}


#' Finds the the outer products for a vector of MS2 ions at a given ion series.
#'
#' A theoretical peptide at a given position permutation.
#'
#' @param expts Numeric vector; one series experimental MS2s.
#' @param theos Numeric vector; one series of theoretical MS2s.
#' @importFrom dplyr bind_cols
#' @inheritParams find_ppm_outer_bypep
#' @examples 
#' \donttest{
#' expts <- c(101.0714, 102.0554, 110.0717, 115.0505, 126.1279, 127.0504, 
#'            127.1249, 127.1312, 128.1283, 128.1346, 129.0660, 129.1316, 
#'            129.1379, 130.1350, 130.1412, 131.1383, 133.0431, 133.0609, 
#'            134.0448, 136.0757, 145.0608, 158.0924, 173.1496, 175.1190, 
#'            176.1594, 191.0663, 193.0971, 198.0873, 201.0869, 201.1233, 
#'            202.0821, 230.1702, 248.0876, 248.1806, 257.1527, 298.6711, 
#'            312.6687, 335.1190, 361.1768, 361.6685, 367.2292, 369.6903, 
#'            370.1821, 370.6830, 376.2756, 377.2790, 384.6952, 412.7264, 
#'            420.2143, 423.2597, 475.2341, 475.3442, 476.3469, 484.2382, 
#'            484.7389, 487.7817, 488.2842, 496.2717, 523.3182, 537.3159, 
#'            537.8178, 572.8343, 573.3373, 623.3586, 624.3298, 636.8528, 
#'            637.3558, 637.8575, 638.4089, 645.3453, 645.8500, 646.3501, 
#'            687.8795, 688.3799, 692.8696, 695.4283, 701.3759, 701.8768, 
#'            702.3782, 702.8797, 737.3956, 737.8962, 738.3671, 739.3572, 
#'            740.3575, 745.9085, 746.4121, 752.4502, 801.4272, 809.4722, 
#'            810.4739, 852.4418, 967.4690, 968.4714, 969.5046, 1084.5293, 
#'            1132.5676, 1133.5686, 1134.5735, 1135.5786)
#'            
#' theos <- c(116.0342, 203.0662, 290.0983, 418.1569, 475.1783, 572.2311, 
#'            732.2617, 861.3043, 958.3571, 1071.4412, 1168.4939, 1225.5154, 
#'            1322.5681, 1435.6522, 1536.6999, 1664.7585, 1761.8112, 1917.9123, 
#'            175.1190, 272.1717, 400.2303, 501.2780, 614.3620, 711.4148, 
#'            768.4363, 865.4890, 978.5731, 1075.6259, 1204.6684, 1364.6991, 
#'            1461.7519, 1518.7733, 1646.8319, 1733.8639, 1820.8960, 1935.9229)
#'            
#' names(theos) <- c("D", "S", "S", "Q", "G", "P", "C", "E", "P", 
#'                   "L", "P", "G", "P", "L", "T", "Q", "P", "R", 
#'                   "R", "P", "Q", "T", "L", "P", "G", "P", "L", 
#'                   "P", "E", "C", "P", "G", "Q", "S", "S", "D")
#' 
#' find_ppm_outer_bycombi(theos, expts)
#' }
#' 
#' 
find_ppm_outer_bycombi <- function (theos, expts, ppm_ms2 = 25) {
  
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_ms2, arr.ind = TRUE)
  
  e1 <- expts[row_cols[, 2]] %>% 
    `names<-`(row_cols[, 1])
  
  len <- length(theos)
  
  es <- rep(NA, len) %>% 
    `names<-`(seq_len(len))
  es[names(e1)] <- e1
  
  # the first half are b-ions and the second half are y-ions
  list(theo = theos, expt = es)
}



