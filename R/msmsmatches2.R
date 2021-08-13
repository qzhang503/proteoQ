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
#' @inheritParams hms2_base
#' @import parallel
ms2match <- function (mgf_path, aa_masses_all, out_path, 
                      mod_indexes, type_ms2ions, maxn_vmods_per_pep, 
                      maxn_sites_per_vmod, maxn_vmods_sitescombi_per_pep, 
                      minn_ms2, ppm_ms1, ppm_ms2, min_ms2mass, 
                      quant, ppm_reporters, 

                      # dummies
                      fasta, acc_type, acc_pattern,
                      topn_ms2ions, fixedmods, varmods, 
                      include_insource_nl, enzyme, 
                      maxn_fasta_seqs, maxn_vmods_setscombi, 
                      min_len, max_len, max_miss, 

                      digits) {
  
  options(digits = 9L)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        save_call2(path = file.path(out_path, "Calls"), fun = "ms2match")
      }
    }, 
    add = TRUE
  )
  
  # Check cached 
  fun <- as.character(match.call()[[1]])
  args_except <- NULL
  
  cache_pars <- find_callarg_vals(time = NULL, 
                                  path = file.path(out_path, "Calls"), 
                                  fun = paste0(fun, ".rda"), 
                                  args = names(formals(fun))) %>% 
    .[! . %in% args_except] %>% 
    .[sort(names(.))]
  
  call_pars <- mget(names(formals()), envir = rlang::current_env(), 
                    inherits = FALSE) %>% 
    .[! . %in% args_except] %>% 
    .[sort(names(.))]
  
  if (identical(cache_pars, call_pars)) {
    len <- length(aa_masses_all)
    
    fions <- list.files(path = file.path(out_path, "temp"), 
                        pattern = "ion_matches_[0-9]+\\.rds$")
    
    ok_ions <- (length(fions) == len)
    
    if (grepl("^tmt[0-9]+$", quant)) {
      ftmt <- list.files(path = file.path(out_path, "temp"), 
                          pattern = "reporters_[0-9]+\\.rds$")
      
      ok_tmt <- (length(ftmt) == len)
    } else {
      ok_tmt <- TRUE
    }
    
    if (ok_ions && ok_tmt) {
      message("Found cached ion matches.")
      .savecall <- FALSE
      
      return(NULL)
    }
  }
  
  rm(list = c("fun", "args_except", "cache_pars", "call_pars"))
  
  delete_files(out_path, ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                                     "\\.xls$", "\\.csv$", "\\.txt$", 
                                     "^mgf$", "^mgfs$"))

  ## Targets 
  obj_sizes <- numeric(length(aa_masses_all))
  types <- purrr::map_chr(aa_masses_all, attr, "type", exact = TRUE)
  
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  inds <- which(types %in% c("amods- tmod- vnl- fnl-", 
                             "amods- tmod+ vnl- fnl-"))
  
  if (length(inds) > 0L) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      
      # need ntmod and ctmod for `amod+_...` for 
      #   excluding additive terminal mod and anywhere mod
      
      # aa_masses["N-term"] not necessary H; 
      #   e.g., `Acetyl + hydrogen` in case of FIXED Protein N-term 
      #   or fixed N-term `TMT + hydrogen`
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod) == 0L) {
        ntmass <- aa_masses["N-term"] - 0.000549 # - electron
      } else {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647 # + proton
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod) == 0L) {
        ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
      } else {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      }
      
      # (`map` against groups of frames)
      out <- ms2match_base(
        i = i, 
        aa_masses = aa_masses, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
      gc()
    }
  }
  
  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  #        (mutual exclusive btw. (1, 2) and (5, 6)
  #         "ANY" fmod has neuloss -> 5, 6;
  #         "ALL" fmods have no neuloss -> 1, 2)
  
  inds <- which(types %in% c("amods- tmod- vnl- fnl+", 
                             "amods- tmod+ vnl- fnl+"))
  
  if (length(inds) > 0L) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod) == 0L) {
        ntmass <- aa_masses["N-term"] - 0.000549
      } else {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod) == 0L) {
        ctmass <- aa_masses["C-term"] + 2.01510147
      } else {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      }
      
      fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
      
      out <- ms2match_a0_vnl0_fnl1(
        i = i, 
        aa_masses = aa_masses, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        fmods_nl = fmods_nl, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
      gc()
    }
  }
  
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  #        (ALL amods are vnl-)
  
  inds <- which(types %in% c("amods+ tmod- vnl- fnl-", 
                             "amods+ tmod+ vnl- fnl-"))
  
  if (length(inds) > 0L) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod) == 0L) {
        ntmass <- aa_masses["N-term"] - 0.000549
      } else {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod) == 0L) {
        ctmass <- aa_masses["C-term"] + 2.01510147
      } else {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      }
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      
      out <- ms2match_a1_vnl0_fnl0(
        i = i, 
        aa_masses = aa_masses, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
      gc()
    }
  }
  
  # (9, 10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
  #         (ANY amod is vnl+)
  
  inds <- which(types %in% c("amods+ tmod- vnl+ fnl-", 
                             "amods+ tmod+ vnl+ fnl-"))
  
  if (length(inds) > 0L) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod) == 0L) {
        ntmass <- aa_masses["N-term"] - 0.000549
      } else {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod) == 0L) {
        ctmass <- aa_masses["C-term"] + 2.01510147
      } else {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      }
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
      
      out <- ms2match_a1_vnl1_fnl0(
        i = i, 
        aa_masses = aa_masses, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        vmods_nl = vmods_nl, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
      gc()
    }
  }
  
  # (11, 12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
  #          (mutual exclusive btw. (11, 12) and (7, 8);
  #           logicial ANY versus ALL)
  
  inds <- which(types %in% c("amods+ tmod- vnl- fnl+", 
                             "amods+ tmod+ vnl- fnl+"))
  
  if (length(inds) > 0L) {
    for (i in inds) {
      aa_masses <- aa_masses_all[[i]]
      
      ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
      if (length(ntmod) == 0L) {
        ntmass <- aa_masses["N-term"] - 0.000549
      } else {
        ntmass <- aa_masses[names(ntmod)] + 1.00727647
      }
      
      ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
      if (length(ctmod) == 0L) {
        ctmass <- aa_masses["C-term"] + 2.01510147
      } else {
        ctmass <- aa_masses[names(ctmod)] + 2.01510147
      }
      
      amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
      fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
      
      out <- ms2match_a1_vnl0_fnl1(
        i = i, 
        aa_masses = aa_masses, 
        ntmod = ntmod, 
        ctmod = ctmod, 
        ntmass = ntmass, 
        ctmass = ctmass, 
        amods = amods, 
        fmods_nl = fmods_nl, 
        mod_indexes = mod_indexes, 
        mgf_path = mgf_path, 
        out_path = out_path, 
        type_ms2ions = type_ms2ions, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        minn_ms2 = minn_ms2, 
        ppm_ms1 = ppm_ms1, 
        ppm_ms2 = ppm_ms2, 
        min_ms2mass = min_ms2mass, 
        digits = digits)
      
      obj_sizes[i] <- object.size(out)
      
      if (grepl("^tmt[0-9]+$", quant)) {
        out <- out %>% 
          calc_tmtint(quant = quant, ppm_reporters = ppm_reporters) %>% 
          tidyr::unite(uniq_id, raw_file, pep_mod_group, scan_num, sep = ".", 
                       remove = TRUE) %>% 
          dplyr::select(uniq_id, grep("^I[0-9]{3}[NC]{0,1}$", names(.))) %T>% 
          saveRDS(file.path(out_path, "temp", paste0("reporters_", i, ".rds")))
      }
      
      rm(list = c("out"))
      gc()
    }
  }
  
  ## Decoys
  maxs <- which.max(obj_sizes)
  maxs2 <- paste0("rev_", maxs)

  aa_masses <- aa_masses_all[[maxs]]
  
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  if (length(ntmod) == 0L) {
    ntmass <- aa_masses["N-term"] - 0.000549
  } else {
    ntmass <- aa_masses[names(ntmod)] + 1.00727647
  }
  
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  if (length(ctmod) == 0L) {
    ctmass <- aa_masses["C-term"] + 2.01510147
  } else {
    ctmass <- aa_masses[names(ctmod)] + 2.01510147
  }
  
  amods <- attr(aa_masses, "amods", exact = TRUE) # variable anywhere
  
  if (length(amods) == 0L) { # (1, 2)
    rev <- ms2match_base(
      i = maxs2, 
      aa_masses = aa_masses, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      mod_indexes = mod_indexes, 
      mgf_path = mgf_path, 
      out_path = out_path, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      minn_ms2 = minn_ms2, 
      ppm_ms1 = ppm_ms1, 
      ppm_ms2 = ppm_ms2, 
      min_ms2mass = min_ms2mass, 
      digits = digits)
  } else { # (7, 8)
    rev <- ms2match_a1_vnl0_fnl0(
      i = maxs2, 
      aa_masses = aa_masses, 
      ntmod = ntmod, 
      ctmod = ctmod, 
      ntmass = ntmass, 
      ctmass = ctmass, 
      amods = amods, 
      mod_indexes = mod_indexes, 
      mgf_path = mgf_path, 
      out_path = out_path, 
      type_ms2ions = type_ms2ions, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      minn_ms2 = minn_ms2, 
      ppm_ms1 = ppm_ms1, 
      ppm_ms2 = ppm_ms2, 
      min_ms2mass = min_ms2mass, 
      digits = digits)
  }
  
  saveRDS(rev, file.path(out_path, "temp", paste0("ion_matches_", maxs2, ".rds"))) 
  
  rm(list = c("rev"))
  gc()
  
  .savecall <- TRUE
  
  invisible(NULL)
}


