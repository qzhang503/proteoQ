#' Matching MS2 ions.
#' 
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#' 
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' @rdname ms2match_base
ms2match_a0_vnl0_fnl1 <- function (i, aa_masses, ntmass, ctmass, fmods_nl, 
                                   mod_indexes, n_cores, cl, mgf_path, out_path, 
                                   type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                                   maxn_sites_per_vmod = 3L, 
                                   maxn_vmods_sitescombi_per_pep = 32L, 
                                   minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                                   min_ms2mass = 110L, digits = 4L) {
  
  tempdata <- purge_search_space(i, aa_masses, mgf_path, n_cores, ppm_ms1)
  mgf_frames <- tempdata$mgf_frames
  theopeps <- tempdata$theopeps
  rm(list = c("tempdata"))
  
  if (is_empty(mgf_frames) || is_empty(theopeps)) return(NULL)
  
  # ---
  out <- hms2_a0_vnl0_fnl1(mgf_frames[[1]], theopeps[[1]], 
                           aa_masses = aa_masses, 
                           ntmass = ntmass, 
                           ctmass = ctmass, 
                           fmods_nl = fmods_nl, 
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
                           min_ms2mass = min_ms2mass, 
                           digits = digits)

  # ---
  out <- clusterMap(cl, hms2_a0_vnl0_fnl1, 
                    mgf_frames, theopeps, 
                    MoreArgs = list(aa_masses = aa_masses, 
                                    ntmass = ntmass, 
                                    ctmass = ctmass, 
                                    fmods_nl = fmods_nl, 
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
                                    min_ms2mass = min_ms2mass, 
                                    digits = digits), 
                    .scheduling = "dynamic") %>% 
    dplyr::bind_rows() %>% # across nodes
    post_ms2match(i, aa_masses, out_path)
  
  rm(list = c("mgf_frames", "theopeps"))
  gc()
  
  invisible(out)
}


#' Searches MGF frames.
#'
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#'
#' @inheritParams ms2match_base
#' @rdname hms2_base
hms2_a0_vnl0_fnl1 <- function (mgf_frames, theopeps, aa_masses, ntmass, ctmass, 
                               fmods_nl, mod_indexes, type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L, 
                               maxn_vmods_sitescombi_per_pep = 32L, 
                               minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                               min_ms2mass = 110L, digits = 4L) {
  
  # `res[[i]]` contains results for multiple mgfs within a frame
  # (the number of entries equals to the number of mgf frames)
  res <- frames_adv_a0_vnl0_fnl1(mgf_frames = mgf_frames, 
                                 theopeps = theopeps, 
                                 aa_masses = aa_masses, 
                                 ntmass = ntmass, 
                                 ctmass = ctmass, 
                                 fmods_nl = fmods_nl, 
                                 mod_indexes = mod_indexes, 
                                 type_ms2ions = type_ms2ions, 
                                 maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                 maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                 maxn_vmods_sitescombi_per_pep = 
                                   maxn_vmods_sitescombi_per_pep, 
                                 minn_ms2 = minn_ms2, 
                                 ppm_ms1 = ppm_ms1, 
                                 ppm_ms2 = ppm_ms2, 
                                 min_ms2mass = min_ms2mass, 
                                 digits = digits) %>% 
    post_frame_adv(mgf_frames)
  
  rm(list = "mgf_frames", "theopeps")
  
  invisible(res)
}


#' Frames advancement.
#'
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#'
#' @rdname frames_adv_base
frames_adv_a0_vnl0_fnl1 <- function (mgf_frames, theopeps, aa_masses, 
                                     ntmass, ctmass, fmods_nl, mod_indexes, 
                                     type_ms2ions = "by", 
                                     maxn_vmods_per_pep = 5L, 
                                     maxn_sites_per_vmod = 3L, 
                                     maxn_vmods_sitescombi_per_pep = 32L, 
                                     minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                                     min_ms2mass = 110L, digits = 4L) {
  
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
                              gen_ms2ions_a0_vnl0_fnl1, 
                              aa_masses = aa_masses, 
                              ntmass = ntmass, 
                              ctmass = ctmass, 
                              fmods_nl = fmods_nl, 
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
                              gen_ms2ions_a0_vnl0_fnl1, 
                              aa_masses = aa_masses, 
                              ntmass = ntmass, 
                              ctmass = ctmass, 
                              fmods_nl = fmods_nl, 
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
                                gen_ms2ions_a0_vnl0_fnl1, 
                                aa_masses = aa_masses, 
                                ntmass = ntmass, 
                                ctmass = ctmass, 
                                fmods_nl = fmods_nl, 
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
                            search_mgf2, 
                            theomasses_bf_ms1, 
                            theomasses_cr_ms1, 
                            theomasses_af_ms1, 
                            theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                            minn_ms2, ppm_ms1, ppm_ms2, min_ms2mass) 

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
                                  gen_ms2ions_a0_vnl0_fnl1, 
                                  aa_masses = aa_masses, 
                                  ntmass = ntmass, 
                                  ctmass = ctmass, 
                                  fmods_nl = fmods_nl, 
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
                                  gen_ms2ions_a0_vnl0_fnl1, 
                                  aa_masses = aa_masses, 
                                  ntmass = ntmass, 
                                  ctmass = ctmass, 
                                  fmods_nl = fmods_nl, 
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
                                  gen_ms2ions_a0_vnl0_fnl1, 
                                  aa_masses = aa_masses, 
                                  ntmass = ntmass, 
                                  ctmass = ctmass, 
                                  fmods_nl = fmods_nl, 
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


#' Calculates the masses of MS2 ion series.
#'
#' (5) "amods- tmod- vnl- fnl+", (6) "amods- tmod+ vnl- fnl+"
#' 
#' @rdname gen_ms2ions_base
#' 
#' @examples 
#' \donttest{
#' # (5) "amods- tmod+ vnl- fnl+"
#' fixedmods <- c("TMT6plex (N-term)", "Oxidation (M)", "dHex (S)")
#' varmods <- c("Acetyl (Protein N-term)")
#' 
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'   
#' aa_masses_all <- calc_aamasses(fixedmods, varmods, add_varmasses = FALSE, 
#'                                add_nlmasses = FALSE)
#' 
#' aa_masses <- aa_masses_all[[3]]
#' 
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#' 
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#' 
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#' 
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#' 
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#' 
#' # variable `TMT6plex (N-term)` + `fixed Oxidation (M)`
#' # (additive varmod on top of fixedmod allowed)
#' 
#' out <- gen_ms2ions_a0_vnl0_fnl1(aa_seq, NULL, aa_masses, ntmass, ctmass, 
#'                                 fmods_nl, mod_indexes)
#' 
#' }
#' 
gen_ms2ions_a0_vnl0_fnl1 <- function (aa_seq, ms1_mass = NULL, aa_masses, 
                                      ntmass = NULL, ctmass = NULL, fmods_nl = NULL, 
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_vmods_sitescombi_per_pep = 32L, 
                                      digits = 4L) {
  
  # (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-" 
  # (no pep_seq dispatching by Anywhere fmod residues -> possible no matched sites)
  
  sites <- names(fmods_nl)
  pattern <- paste(sites, collapse = "|")
  
  if (!grepl(pattern, aa_seq)) {
    out <- gen_ms2ions_base(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                            aa_masses = aa_masses, 
                            ntmass = ntmass, ctmass = ctmass, 
                            mod_indexes = mod_indexes, 
                            type_ms2ions = type_ms2ions, 
                            maxn_vmods_per_pep = maxn_vmods_per_pep, 
                            maxn_sites_per_vmod = maxn_sites_per_vmod, 
                            maxn_vmods_sitescombi_per_pep = 
                              maxn_vmods_sitescombi_per_pep, 
                            digits = digits)
    return(out)
  }

  # (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+" 
  aas <- aa_seq %>% str_split("", simplify = TRUE)
  
  idxes <- which(aas %in% names(fmods_nl))
  
  # At varmods "Oxidation (M)", pep_seq(s) must contain "M" 
  #   (with an additional entry of "Oxidation (M)" in aa_masses)
  # 
  # At fixedmods "Oxidation (M)", pep_seq(s) may not contain "M"; 
  #   (as `distri_peps` does not filter pep_seq by fixedmods)
  
  len_i <- length(idxes)

  if (len_i > maxn_vmods_per_pep) {
    len_i <- maxn_vmods_per_pep
    idxes <- idxes[1:len_i]
  }
  
  # to check maxn_sites_per_vmod (R table slow)?
  
  # ---
  fmods_combi <- aas[idxes]
  names(fmods_combi) <- idxes
  fnl_combi <- expand.grid(fmods_nl[fmods_combi])
  
  len <- nrow(fnl_combi)
  
  if (len > maxn_vmods_sitescombi_per_pep) {
    len <- maxn_vmods_sitescombi_per_pep
    fnl_combi <- fnl_combi[1:len, ]
  }
  
  out <- vector("list", len)
  
  aas2 <- aa_masses[aas]

  # do not change to `2:len`
  for (i in 1:len) {
    aas2_i <- aas2
    delta <- unlist(fnl_combi[i, ], use.names = FALSE)
    aas2_i[idxes] <- aas2_i[idxes] - delta
    out[[i]] <- ms2ions_by_type(aas2_i, ntmass, ctmass, type_ms2ions, digits)
  }
  
  nm <- rep(0, length(aas)) %>% paste0(collapse = "")
  names(out) <- paste0(nm, " [", seq_len(nrow(fnl_combi)), "]")
  
  invisible(out)
}



