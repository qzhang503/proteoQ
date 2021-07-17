#' Generates and Calculates the masses of tryptic peptides from a fasta
#' database.
#'
#' @param index_mods Logical; if TRUE, converts the names of modifications to
#'   indexes. Not currently used.
#' @param enzyme A character string; the proteolytic specificity of the assumed
#'   enzyme will be used to generate peptide sequences from proteins. The enzyme
#'   is currently \code{trypsin}.
#' @param maxn_fasta_seqs Integer; the maximum number of protein sequences in
#'   fasta files.
#' @param min_len Integer; the minimum length of peptides. Shorter peptides will
#'   be excluded.
#' @param max_len Integer; the maximum length of peptides. Longer peptides will
#'   be excluded.
#' @param max_miss The maximum number of mis-cleavages per peptide sequence.
#' @param digits Integer; the number of decimal places to be used.
#' @inheritParams splitPSM
#' @inheritParams binTheoPeps
#' @inheritParams calc_aamasses
#' @inheritParams mcalc_monopep
#' @inheritParams calc_monopep
#' @inheritParams load_fasta2
#' @import parallel
#' @importFrom rlang current_env
#' @examples
#' \donttest{
#' res <- calc_pepmasses()
#'
#' library(purrr)
#' library(magrittr)
#'
#' res_attrs <- map(res, attributes)
#' map(res_attrs, names)
#' res_attrs %>% map(`[[`, "vmods")
#' res_mods <- map(res_attrs, `[`,
#'                 c("fmods", "fmods_ps", "fmods_neuloss",
#'                   "vmods", "vmods_ps", "vmods_neuloss"))
#'
#' res_data <- map(res_attrs, `[[`, "data")
#' peps_combi_1 <- res_data[[1]]
#'
#' # base: fixedmods without neulosses
#' unlist(res_data[[1]], use.names = FALSE) %>% length()
#'
#' # fixedmods, fixedmods + fixedmods_neulosses, varmods, varmods_neulosses
#' unlist(res_data, use.names = FALSE) %>% length()
#'
#' }
#'
#' @export
calc_pepmasses2 <- function (
  fasta = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
  acc_type = "uniprot_acc", 
  acc_pattern = NULL, 
  fixedmods = c("TMT6plex (K)", 
                "Carbamidomethyl (C)"), 
  varmods = c("TMT6plex (N-term)", 
              "Acetyl (Protein N-term)", 
              "Oxidation (M)", 
              "Deamidated (N)", 
              "Gln->pyro-Glu (N-term = Q)"), 
  include_insource_nl = FALSE, 
  index_mods = FALSE, 
  enzyme = c("trypsin"), 
  maxn_fasta_seqs = 50000L,
  maxn_vmods_setscombi = 64L, 
  maxn_vmods_per_pep = 5L,
  maxn_sites_per_vmod = 3L, 
  min_len = 7L, max_len = 100L, max_miss = 2L, 
  out_path = "~/proteoQ/outs", 
  digits = 4L, 
  parallel = TRUE) {
  
  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = current_env())) {
      if (.savecall) {
        save_call2(path = file.path(.path_cache), fun = "calc_pepmasses", 
                   time = .time_stamp)
      }
    } else {
      # warning("terminated abnormally.", call. = TRUE)
    }, 
    add = TRUE
  )
  
  # ---
  stopifnot(vapply(c(min_len, max_len, max_miss), is.numeric, 
                   logical(1L)))
  
  stopifnot(min_len >= 0L, max_len >= min_len, max_miss <= 100L)
  
  # ---
  .path_cache <- create_dir("~/proteoQ/.MSearch/Cache/Calls")
  .time_stamp <- match_calltime(path = .path_cache, 
                                fun = "calc_pepmasses", 
                                nms = c("fasta", "acc_type", "acc_pattern", 
                                        "fixedmods", "varmods", 
                                        "include_insource_nl", "enzyme", 
                                        "maxn_fasta_seqs", "maxn_vmods_setscombi", 
                                        "maxn_vmods_per_pep", "maxn_sites_per_vmod", 
                                        "min_len", "max_len", "max_miss"))
  
  .path_fasta <- fasta %>% 
    gsub("\\\\", "/", .) %>% 
    gsub("(.*/)[^/]*", "\\1", .) %>% 
    map_chr(find_dir) %>% 
    `[[`(1)
  
  # ---
  if (!is_empty(.time_stamp)) {
    message("Loading peptide masses from cache.")
    
    out <- readRDS(file.path(.path_fasta, "pepmasses", .time_stamp, 
                             "pepmasses.rds"))
    
    rev_out <- readRDS(file.path(.path_fasta, "pepmasses", .time_stamp, 
                                 "pepmasses_rev.rds"))

    .savecall <- FALSE
  } else {
    delete_files(out_path, ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                                       "\\.xls$", "\\.csv$", "\\.txt$", 
                                       "^mgf$", "^mgfs$"))
    
    .time_stamp <- format(Sys.time(), ".%Y-%m-%d_%H%M%S")

    # --- 
    aa_masses <- local({
      message("Computing the combinations of fixed and variable modifications.")
      
      if (index_mods) {
        mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
          as.hexmode() %>% 
          `names<-`(c(fixedmods, varmods))
      } else {
        mod_indexes <- NULL
      }
      
      aa_masses <- calc_aamasses(fixedmods = fixedmods, 
                                 varmods = varmods, 
                                 maxn_vmods_setscombi = maxn_vmods_setscombi, 
                                 mod_indexes = mod_indexes, 
                                 add_varmasses = FALSE, 
                                 add_nlmasses = FALSE)
    })
    
    gc()
    
    # --- forward and reversed sequences  ---
    # (not yet concatenation by the number of mis-cleavages)
    seqs_0 <- split_fastaseqs(fasta = fasta, 
                              acc_type = acc_type, 
                              acc_pattern = acc_pattern, 
                              maxn_fasta_seqs = maxn_fasta_seqs, 
                              max_miss = max_miss)

    # --- Precursor masses (bared) --- 
    # (H and OH for terminals; no varmods, no NLs)
    message("Calculating peptide masses (fixed modifications) ...")
    
    fwd_peps <- ms1masses_bare(seqs = seqs_0, 
                               aa_masses = aa_masses[[1]], 
                               max_miss = max_miss, 
                               min_len = min_len, 
                               max_len = max_len, 
                               maxn_vmods_per_pep = maxn_vmods_per_pep, 
                               maxn_sites_per_vmod = maxn_sites_per_vmod, 
                               parallel = parallel, 
                               digits = digits) 
    
    message("\tCompleted target peptide masses: ", 
            paste(attributes(aa_masses[[1]])$fmods, 
                  attributes(aa_masses[[1]])$vmods, 
                  collapse = ", "))
    
    # ---
    message("Distributing target peptides distributions ", 
            "by fixed and variable modifications.")
    
    fwd_peps <-  distri_peps(aa_masses = aa_masses, 
                             peps = fwd_peps, 
                             max_miss = max_miss) # %T>% 
    # saveRDS(file.path(out_path, "fwds0_masses.rds"))
    
    message("\tCompleted peptides distributions.")

    rm(list = c("seqs_0"))
    gc()

    # --- initialization ---
    len <- length(aa_masses)
    types <- purrr::map_chr(aa_masses, attr, "type", exact = TRUE)
    
    out <- vector("list", len)
    out[[1]] <- fwd_peps[[1]]

    # --- add terminal masses delta ---
    message("Adding terminal masses (fixed and variable modifications) ...")
    
    inds <- grep("tmod+", types, fixed = TRUE)
    
    if (length(inds) > 0L) {
      for (i in inds) {
        out[[i]] <- add_term_mass(aa_masses[[i]], fwd_peps[[i]])
      }
    }
    
    # --- combinatorial ---
    message("Calculating peptide masses (variable modifications) ...")
    
    n_cores <- detectCores()
    cl <- makeCluster(getOption("cl.cores", n_cores))
    
    fmods_ps <- purrr::map(aa_masses, attr, "fmods_ps", exact = TRUE)
    vmods_ps <- purrr::map(aa_masses, attr, "vmods_ps", exact = TRUE)
    fmods_nl <- purrr::map(aa_masses, attr, "fmods_nl", exact = TRUE)
    vmods_nl <- purrr::map(aa_masses, attr, "vmods_nl", exact = TRUE)
    amods <- purrr::map(aa_masses, attr, "amods", exact = TRUE)
    tmod <- purrr::map(aa_masses, attr, "tmod", exact = TRUE)
    
    # `amods-` and `fnl+` (must be vnl- since amods-)
    # 
    # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
    
    if (include_insource_nl) {
      inds <- which(types %in% c("amods- tmod- vnl- fnl+", 
                                 "amods- tmod+ vnl- fnl+"))
      
      if (length(inds) > 0L) {
        for (i in inds) {
          fmods_ps_i <- fmods_ps[[i]]
          vmods_ps_i <- vmods_ps[[i]]
          fmods_nl_i <- fmods_nl[[i]]
          vmods_nl_i <- vmods_nl[[i]]
          amods_i <- amods[[i]]
          tmod_i <- tmod[[i]]
          
          aa_masses_i <- aa_masses[[i]]
          ntmod_i <- attr(aa_masses_i, "ntmod", exact = TRUE)
          ctmod_i <- attr(aa_masses_i, "ctmod", exact = TRUE)
          
          fwd_peps_i <- fwd_peps[[i]]

          fnl_combi_i <- expand.grid(fmods_nl_i)
          
          clusterExport(cl, list("ms1_a0_fnl1_byprot"), 
                        envir = environment(proteoQ:::ms1_a0_fnl1_byprot))
          clusterExport(cl, list("ms1_a0_fnl1_bypep"), 
                        envir = environment(proteoQ:::ms1_a0_fnl1_bypep))
          clusterExport(cl, list("delta_ms1_a0_fnl1"), 
                        envir = environment(proteoQ:::delta_ms1_a0_fnl1))
          
          out[[i]] <- parLapply(cl, fwd_peps_i, ms1_a0_fnl1_byprot, 
                                fnl_combi_i, aa_masses_i, digits = digits)
        }
      }
    }

    # `amods+`; (9-14) nested under (7-8)
    # 
    # (7-8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
    #   (9-10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
    #   (11-12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
    #   (13-14) "amods+ tmod- vnl+ fnl+", "amods+ tmod+ vnl+ fnl+"
    
    inds <- which(types %in% c("amods+ tmod- vnl- fnl-", 
                               "amods+ tmod+ vnl- fnl-", 
                               "amods+ tmod- vnl+ fnl-", 
                               "amods+ tmod+ vnl+ fnl-", 
                               "amods+ tmod- vnl- fnl+", 
                               "amods+ tmod+ vnl- fnl+", 
                               "amods+ tmod- vnl+ fnl+", 
                               "amods+ tmod+ vnl+ fnl+"))
    
    if (length(inds) > 0L) {
      clusterExport(cl, list("ms1_a1_vnl0_fnl0_byprot"), 
                    envir = environment(proteoQ:::ms1_a1_vnl0_fnl0_byprot))
      clusterExport(cl, list("ms1_a1_vnl0_fnl0_bypep"), 
                    envir = environment(proteoQ:::ms1_a1_vnl0_fnl0_bypep))
      
      for (i in inds) {
        amods_i <- amods[[i]]
        aa_masses_i <- aa_masses[[i]]
        
        fwd_peps_i <- fwd_peps[[i]]

        vmods_nl_i = vmods_nl[[i]]
        fmods_nl_i = fmods_nl[[i]]

        out[[i]] <- parLapply(cl, fwd_peps_i, ms1_a1_vnl0_fnl0_byprot, 
                              amods_i, aa_masses_i, 
                              vmods_nl_i, fmods_nl_i, include_insource_nl, 
                              maxn_vmods_per_pep, 
                              maxn_sites_per_vmod, 
                              digits = digits)

        message("\tCompleted target peptide masses: ", 
                paste(attributes(aa_masses[[i]])$fmods, 
                      attributes(aa_masses[[i]])$vmods, 
                      collapse = ", "))
      }
    }
    
    stopCluster(cl)
    gc()
    
    # --- Reversed ---
    rev_out <- purrr::map(out, rev_pepseqs)
    
    # --- outputs ---
    out <- purrr::map2(aa_masses, out, ~ {
      attr(.x, "data") <- .y
      .x
    })
    
    rev_out <- purrr::map2(aa_masses, rev_out, ~ {
      attr(.x, "data") <- .y
      .x
    })
    
    out_path <- create_dir(file.path(.path_fasta, "pepmasses", .time_stamp))
    saveRDS(out, file.path(out_path, "pepmasses.rds"))
    saveRDS(rev_out, file.path(out_path, "pepmasses_rev.rds"))
    
    # ---
    .savecall <- TRUE
  }
  
  assign(".path_cache", .path_cache, envir = .GlobalEnv)
  assign(".time_stamp", .time_stamp, envir = .GlobalEnv)
  assign(".path_fasta", .path_fasta, envir = .GlobalEnv)
  
  invisible(list(fwd = out, rev = rev_out))
}




#' Helper in calculating peptide masses.
#' 
#' (2) "amods- tmod+ vnl- fnl-".
#' 
#' @inheritParams distri_peps
add_term_mass <- function (aa_masses, peps) {
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  
  # No: is_empty(ntmod) && is_empty(ctmod)
  
  if (!(is_empty(ntmod) || is_empty(ctmod))) {
    delta <- aa_masses[names(ntmod)] + aa_masses[names(ctmod)]
  } else if (!is_empty(ntmod)) {
    delta <- aa_masses[names(ntmod)]
  } else if (!is_empty(ctmod)) {
    delta <- aa_masses[names(ctmod)]
  }
  
  purrr::map(peps, `+`, delta)
}


#' Reverses peptide sequences.
#' 
#' @param xs Lists of proteins with named peptide masses under each entry.
rev_pepseqs <- function (xs) {
  
  names(xs) <- paste0("-", names(xs))
  
  rev_seqs <- purrr::map(xs, ~ {
    seqs <- names(.x)
    a <- stringi::stri_sub(seqs, 1, 1)
    b <- stringi::stri_sub(seqs, -1, -1)
    
    lens <- stri_length(seqs)
    
    revs <- stringi::stri_reverse(seqs)
    substring(revs, 1) <- a
    substring(revs, lens) <- b
    
    revs
  })
  
  out <- purrr::map2(xs, rev_seqs, ~ {
    names(.x) <- .y
    .x
  })
}


#' Make peptide sequences from FASTA databases.
#' 
#' Not yet concatenating peptides by the number of mis-cleavages.
#' 
#' @param fasta_db Fasta database(s).
#' @inheritParams calc_pepmasses
make_fastapeps0 <- function (fasta_db, max_miss = 2L) {
  
  inds_m <- grep("^M", fasta_db)
  
  fasta_db <- fasta_db %>% 
    purrr::map(~ gsub("([KR]{1})", paste0("\\1", "@"), .x) %>% 
                 paste0("-", ., "-")) 
  
  fasta_dbm <- fasta_db[inds_m] %>% 
    purrr::map(~ gsub("^-M", "-", .x))
  
  # --- with protein N-term (initiator) methionine ---
  peps <- fasta_db %>% 
    purrr::map(~ .x %>% stringr::str_split("@", simplify = TRUE))
  
  # --- without protein N-term (initiator) methionine ---
  peps_m <- fasta_dbm %>% 
    purrr::map(~ .x %>% 
                 stringr::str_split("@", simplify = TRUE) %>% 
                 keep_n_misses(max_miss))
  
  # Note: NA sequences -> NULL during mass calculations
  peps[inds_m] <- purrr::map2(peps_m, peps[inds_m], list)
  peps[-inds_m] <- purrr::map(peps[-inds_m], ~ list(NA, .x))
  
  rm(list = c("inds_m", "fasta_db", "fasta_dbm", "peps_m"))
  
  invisible(peps)
}


#' Helper of \link{calc_pepmasses}.
#'
#' Prior to the calculation of peptide masses; for the base with fixed
#' modification only.
#' 
#' @inheritParams calc_pepmasses
#' @inheritParams mcalc_monopep
#' @importFrom stringi stri_reverse
#' @return Two named list of "fwds" and "revs". List "fwds" contains peptide
#'   sequences split from forward fasta and "revs" from reversed fasta.
split_fastaseqs <- function (fasta, acc_type, acc_pattern, maxn_fasta_seqs, 
                             max_miss = 2L) {
  
  # ---
  message("Loading fasta databases.")
  
  fasta_db <- load_fasta2(fasta, acc_type, acc_pattern) 
  
  if (length(fasta_db) > maxn_fasta_seqs) {
    stop("More than `", maxn_fasta_seqs, "` sequences in fasta files.\n",
         "  May consider a higher `maxn_fasta_seqs`.", 
         call. = FALSE)
  }
  
  # ---
  # message("Reversing protein sequences.")
  
  n_cores <- detectCores()
  cl <- makeCluster(getOption("cl.cores", n_cores))
  
  clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  clusterExport(cl, list("keep_n_misses"), envir = environment(proteoQ:::keep_n_misses))
  
  fasta_db <- chunksplit(fasta_db, n_cores)
  
  # rev_fasta_db <- clusterApply(cl, fasta_db, stringi::stri_reverse) %>% 
  #   purrr::map2(., fasta_db, ~ {
  #     names(.x) <- paste0("-", names(.y))
  #     .x
  #   })

  # ---
  message("Splitting fasta sequences.")
  
  peps <- clusterApply(cl, fasta_db, make_fastapeps0, max_miss) %>% 
    purrr::flatten()
  
  # rev_peps <- clusterApply(cl, rev_fasta_db, make_fastapeps0, max_miss) %>% 
  #   purrr::flatten()
  
  # rm(list = c("fasta_db", "rev_fasta_db"))
  rm(list = c("fasta_db"))
  stopCluster(cl)
  
  # invisible(list(fwds = peps, revs = rev_peps))
  invisible(peps)
}


#' Distributes peptides by variable modifications.
#'
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @param peps Lists of peptides, either "fwds" or "revs", from
#'   \link{split_fastaseqs}.
#' @inheritParams calc_pepmasses
distri_peps <- function (aa_masses, peps, max_miss) {
  nms <- purrr::map(peps, names)
  
  out <- aa_masses %>% 
    purrr::map(subpeps_by_vmods, nms) 
  
  out <- out %>% 
    purrr::map(~ {
      xs <- .x
      idxes <- map2(xs, nms, fastmatch::fmatch)
      
      map2(peps, idxes, ~ .x[.y])
    })
  
  n2 <- (max_miss + 1L)
  n1 <- n2 * 2L
  
  out <- out %>% 
    purrr::map(~ {
      xs <- .x
      
      # ZN207_HUMAN: MGRKKKK (no N-term pep_seq at 2 misses and min_len >= 7L)
      len <- purrr::map_int(xs, length)
      xs <- xs[len > 0L]
      
      xs %>% 
        purrr::map(rm_char_in_nfirst2, char = "^-", n = n1) %>% 
        purrr::map(rm_char_in_nlast2, char = "-$", n = n2)
        # purrr::map(rm_char_in_nlast2, char = "-$", n = (max_miss + 1) * 2L)
    })
}


#' Remove a starting character from the first \code{n} entries.
#' 
#' @param x A list of character strings.
#' @param char A starting character to be removed.
#' @param n The number of beginning entries to be considered.
rm_char_in_nfirst2 <- function (x, char = "^-", n = (max_miss + 1L) * 2L) {
  nms <- names(x)
  
  len <- length(nms)
  n <- min(len, n)
  
  seqs <- seq_len(n)
  
  nms[seqs] <- gsub(char, "", nms[seqs])
  names(x) <- nms
  
  x
}


#' Remove a trailing character from the last \code{n} entries.
#' 
#' @param char A trailing character to be removed.
#' @inheritParams rm_char_in_nfirst
rm_char_in_nlast2 <- function (x, char = "-$", n = (max_miss + 1L) * 2L) {
  nms <- names(x)
  
  len <- length(nms)
  n <- min(len, n)
  
  seqs <- (len-n+1):len
  
  nms[seqs] <- gsub(char, "", nms[seqs])
  names(x) <- nms

  x
}


#' Calculates mono-isotopic masses of peptide sequences.
#'
#' @param seqs Sequences of peptides from FASTAs by protein accessions. Each
#'   list contains two lists of sequences: (1) without and (2) with N-terminal
#'   methionine.
#' @inheritParams calc_pepmasses
#' @inheritParams mcalc_monopep
ms1masses_bare <- function (seqs, aa_masses, max_miss = 2L, 
                            min_len = 7L, max_len = 100L, 
                            maxn_vmods_per_pep = 5L, 
                            maxn_sites_per_vmod = 3L, parallel = TRUE, 
                            digits = 4L) {
  
  # (1) before rolling sum (not yet terminal H2O)
  # (1.1) without N-term methionine
  data_1 <- seqs %>% 
    purrr::map(`[[`, 1) %>% 
    ms1masses_noterm(aa_masses = aa_masses, 
                     maxn_vmods_per_pep = maxn_vmods_per_pep, 
                     maxn_sites_per_vmod = maxn_sites_per_vmod, 
                     parallel = parallel, 
                     digits = digits) %>% 
    attr("data") 
  
  # (1.2) with N-term methionine
  data_2 <- seqs %>% 
    purrr::map(`[[`, 2) %>% 
    ms1masses_noterm(aa_masses = aa_masses, 
                     maxn_vmods_per_pep = maxn_vmods_per_pep, 
                     maxn_sites_per_vmod = maxn_sites_per_vmod, 
                     parallel = parallel, 
                     digits = digits) %>% 
    attr("data") 
  
  # (2) rolling sum (not yet terminal H2O)
  n_cores <- detectCores()
  cl <- makeCluster(getOption("cl.cores", n_cores))
  
  clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  clusterExport(cl, list("roll_sum"), envir = environment(proteoQ:::roll_sum))
  
  ms_1 <- parLapply(cl, data_1, roll_sum, max_miss, include_cts = FALSE)
  ms_2 <- parLapply(cl, data_2, roll_sum, max_miss, include_cts = TRUE)
  
  stopCluster(cl)
  
  # (3) putting together (+ terminal H2O)
  ms <- purrr::map2(ms_1, ms_2, `c`)
  
  if (min_len > 1L && !is.infinite(max_len)) {
    ms <- ms %>% 
      purrr::map(
        ~ .x %>% .[str_exclude_count(names(.)) >= min_len & 
                     str_exclude_count(names(.)) <= max_len]) 
  }
  
  # "HD101_HUMAN" etc. has no tryptic peptides
  lens <- purrr::map_int(ms, length)
  
  # adding H2O mass
  ms <- ms %>% 
    .[lens > 0L] %>% 
    purrr::map(~ .x[!duplicated(names(.x))]) %>% 
    purrr::map(~ .x + 18.010565)
  
  invisible(ms)
}



#' Helper of \link{ms1masses_bare}.
#' 
#' For either forward or reversed sequences.
#' 
#' @inheritParams mcalc_monopep
ms1masses_noterm <- function (aa_seqs, aa_masses, 
                              maxn_vmods_per_pep = 5L, 
                              maxn_sites_per_vmod = 3L, 
                              parallel = TRUE, 
                              digits = 4L) {
  
  options(digits = 9)
  
  n_cores <- detectCores()
  
  aa_seqs <- suppressWarnings(split(aa_seqs, seq_len(n_cores)))
  
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(cl, list("%>%"), 
                          envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("calcms1mass_noterm"), 
                          envir = environment(proteoQ:::calcms1mass_noterm))
  parallel::clusterExport(cl, list("calcms1mass_noterm_byprot"), 
                          envir = environment(proteoQ:::calcms1mass_noterm_byprot))
  parallel::clusterExport(cl, list("calcms1mass_noterm_bypep"), 
                          envir = environment(proteoQ:::calcms1mass_noterm_bypep))
  
  out <- parallel::clusterApply(cl, aa_seqs, calcms1mass_noterm, 
                                aa_masses = aa_masses, 
                                maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                digits = digits) %>% 
    purrr::flatten()
  
  parallel::stopCluster(cl)
  
  attr(aa_masses, "data") <- out
  
  # message("\tCompleted peptide masses: ", 
  #         paste(attributes(aa_masses)$fmods, 
  #               attributes(aa_masses)$vmods, 
  #               collapse = ", "))
  
  rm(list = c("aa_seqs", "out"))
  gc()
  
  invisible(aa_masses)
}


#' Helper function for parallel calculations peptide masses by proteins.
#' 
#' For each split of multiple proteins; no terminal masses.
#' 
#' @param prot_peps Lists of peptides under a proteins.
#' @inheritParams calc_monopep
calcms1mass_noterm <- function (aa_seqs, aa_masses, 
                                maxn_vmods_per_pep = 5L, maxn_sites_per_vmod = 3L, 
                                digits = 4L) {
  
  purrr::map(aa_seqs, ~ {
    prot_peps <- .x
    
    calcms1mass_noterm_byprot(prot_peps = prot_peps, 
                              aa_masses = aa_masses, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              digits = digits)
  })
}


#' Helper of \link{calcms1mass_noterm}.
#' 
#' For single protein.
#' 
#' @param prot_peps Lists of peptides under a proteins.
#' @inheritParams calc_monopep
calcms1mass_noterm_byprot <- function (prot_peps, aa_masses, 
                                       maxn_vmods_per_pep = 5L, 
                                       maxn_sites_per_vmod = 3L, 
                                       digits = 4L) {
  
  prot_peps %>% 
    purrr::map(calcms1mass_noterm_bypep, 
               aa_masses = aa_masses, 
               maxn_vmods_per_pep = maxn_vmods_per_pep, 
               maxn_sites_per_vmod = maxn_sites_per_vmod, 
               digits = digits) %>% # by peptides
    unlist(use.names = TRUE)
}


#' Helper of \link{calcms1mass_noterm_byprot}.
#' 
#' For single protein.
#' 
#' @inheritParams calc_monopep
#' @importFrom stringr str_split
calcms1mass_noterm_bypep <- function (aa_seq, aa_masses, 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      digits = 4L) {
  
  if (is.na(aa_seq)) return(NULL)
  
  aa_seq %>% 
    str_split("", simplify = TRUE) %>% 
    aa_masses[.] %>% 
    sum() %>% 
    setNames(aa_seq) %>% 
    round(digits = digits)
}












#' Converts protein accessions and peptide sequences to integers.
#'
#' @param pep_seqs Lists of peptide sequences by protein accessions. Peptide
#'   sequences are names and monoisotopic masses in values.
#' @param out_path An output path
#' @param type An indicator of either forward or reversed sequences.
#' @examples
#' \donttest{
#' pep_seqs <- list(ALG2_HUMAN = c("-MAEEQGR" = 819.35450, 
#'                                 "-MAEEQGRER" = 1122.50877 ))
#' 
#' out <- conv_to_pepids(pep_seqs, out_path, "rev")
#' }
conv_to_pepids <- function (data, out_path, type = "fwd") {
  
  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
  
  # pep_seqs <- purrr::map(data, attr, "data")
  
  # two-stage matches
  #   prot_acc -> pep_seq(s); find position indexes
  
  lens <- purrr::map_int(pep_seqs, length)
  
  prots <- purrr::imap(lens, ~ rep(.y, .x)) %>% 
    unlist(use.names = FALSE)
  
  peps <- pep_seqs %>% 
    purrr::map(names) %>% 
    unlist(use.names = FALSE)
  
  df_protids <- local({
    prs <- names(pep_seqs)
    data.frame(prot_acc = prs, prot_index = seq_along(prs))
  })
  
  df_pepids <- local({
    pes <- unique(peps)
    data.frame(pep_seq = pes, pep_index = seq_along(pes))
  })
  
  df <- data.frame(prot_acc = prots, pep_seq = peps) %>% 
    dplyr::left_join(df_protids, by = "prot_acc") %>% 
    dplyr::left_join(df_pepids, by = "pep_seq") %>% 
    dplyr::group_by(pep_seq) %>% 
    dplyr::mutate(n = row_number()) %T>% 
    saveRDS(file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  stopifnot(identical(unique(df$prot_acc), names(pep_seqs)))
  
  dfs <- df %>% split(.$prot_index)
  # saveRDS(dfs, file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  names(pep_seqs) <- names(dfs)
  
  pep_seqs <- purrr::map2(pep_seqs, dfs, ~ {
    names(.x) <- .y[["pep_index"]]
    .x
  })
}








conv_to_pepids_v1 <- function (pep_seqs, out_path, type = "fwd") {
  
  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
  
  # pep_seqs <- purrr::map(pep_seqs, attr, "data")
  
  lens <- purrr::map_int(pep_seqs, length)
  
  prots <- purrr::imap(lens, ~ rep(.y, .x)) %>% 
    unlist(use.names = FALSE)
  
  peps <- pep_seqs %>% 
    purrr::map(names) %>% 
    unlist(use.names = FALSE)
  
  df_protids <- local({
    prs <- names(pep_seqs)
    data.frame(prot_acc = prs, prot_index = seq_along(prs))
  })
  
  df_pepids <- local({
    pes <- unique(peps)
    data.frame(pep_seq = pes, pep_index = seq_along(pes))
  })
  
  df <- data.frame(prot_acc = prots, pep_seq = peps) %>% 
    dplyr::left_join(df_protids, by = "prot_acc") %>% 
    dplyr::left_join(df_pepids, by = "pep_seq") %>% 
    dplyr::group_by(pep_seq) %>% 
    dplyr::mutate(n = row_number()) %T>% 
    saveRDS(file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  stopifnot(identical(unique(df$prot_acc), names(pep_seqs)))
  
  dfs <- df %>% split(.$prot_index)
  # saveRDS(dfs, file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  names(pep_seqs) <- names(dfs)
  
  pep_seqs <- purrr::map2(pep_seqs, dfs, ~ {
    names(.x) <- .y[["pep_index"]]
    .x
  })
}



