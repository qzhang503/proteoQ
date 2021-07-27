#' Matching MS2 ions.
#' 
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#' 
#' @param i Integer; the index for a set of corresponding aa_masses and
#'   theoretical peptides.
#' @param aa_masses Amino-acid lookup at a given combination of fixed and
#'   variable.
#' @param ntmass The mass of N-terminal.
#' @param ctmass The mass of C-terminal.
#' @param n_cores The number of CPU cores.
#' @param cl The cluster nodes.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @inheritParams matchMS
#' @inheritParams ms2match
ms2match_base <- function (i, aa_masses, ntmass, ctmass, mod_indexes, 
                           n_cores, cl, mgf_path, out_path, 
                           type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                           maxn_sites_per_vmod = 3L, 
                           maxn_vmods_sitescombi_per_pep = 32L, 
                           minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                           digits = 4L) {
  
  tempdata <- purge_search_space(i, aa_masses, mgf_path, n_cores, ppm_ms1)
  mgf_frames <- tempdata$mgf_frames
  theopeps <- tempdata$theopeps
  rm(list = c("tempdata"))
  
  if (is_empty(mgf_frames) || is_empty(theopeps)) return(NULL)
  
  out <- clusterMap(cl, hms2_base, 
                    mgf_frames, theopeps, 
                    MoreArgs = list(aa_masses = aa_masses, 
                                    ntmass = ntmass, 
                                    ctmass = ctmass, 
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
    post_ms2match(i, aa_masses, out_path)
  
  rm(list = c("mgf_frames", "theopeps"))
  gc()
  
  invisible(out)
}


#' Searches MGF frames.
#'
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#'
#' @param mgf_frames A group of mgf frames. Each frame contains one to multiple
#'   mgf queries with the same frame number.
#' @param theopeps Binned theoretical peptides corresponding to an i-th
#'   \code{aa_masses}.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams ms2match_base
hms2_base <- function (mgf_frames, theopeps, aa_masses, ntmass, ctmass, 
                       mod_indexes, type_ms2ions = "by", 
                       maxn_vmods_per_pep = 5L, 
                       maxn_sites_per_vmod = 3L, 
                       maxn_vmods_sitescombi_per_pep = 32L, 
                       minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                       digits = 4L) {
  
  # `res[[i]]` contains results for multiple mgfs within a frame
  # (the number of entries equals to the number of mgf frames)
  
  res <- frames_adv_base(mgf_frames = mgf_frames, 
                         theopeps = theopeps, 
                         aa_masses = aa_masses, 
                         ntmass = ntmass, 
                         ctmass = ctmass, 
                         mod_indexes = mod_indexes, 
                         type_ms2ions = type_ms2ions, 
                         maxn_vmods_per_pep = maxn_vmods_per_pep, 
                         maxn_sites_per_vmod = maxn_sites_per_vmod, 
                         maxn_vmods_sitescombi_per_pep = 
                           maxn_vmods_sitescombi_per_pep, 
                         minn_ms2 = minn_ms2, 
                         ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
                         digits = digits) %>% 
    post_frame_adv(mgf_frames)

  rm(list = "mgf_frames", "theopeps")
  
  invisible(res)
}


#' Frames advancement.
#' 
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#' 
#' @param mgf_frames MGFs in frames. Each frame contains one to multiple MGFs
#'   whose MS1 masses are in the same interval.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams ms2match_base
#' @inheritParams hms2_base
#' @return Matches to each MGF as a list elements. The length of the output is
#'   equal to the number of MGFs in the given frame.
frames_adv_base <- function (mgf_frames, theopeps, aa_masses, 
                             ntmass, ctmass, mod_indexes, 
                             type_ms2ions = "by", 
                             maxn_vmods_per_pep = 5L, 
                             maxn_sites_per_vmod = 3L, 
                             maxn_vmods_sitescombi_per_pep = 32L, 
                             minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                             digits = 4L) {
  
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
                              gen_ms2ions_base, 
                              aa_masses = aa_masses, 
                              ntmass = ntmass, 
                              ctmass = ctmass, 
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
                              gen_ms2ions_base, 
                              aa_masses = aa_masses, 
                              ntmass = ntmass, 
                              ctmass = ctmass, 
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
                                gen_ms2ions_base, 
                                aa_masses = aa_masses, 
                                ntmass = ntmass, 
                                ctmass = ctmass, 
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
                                  gen_ms2ions_base, 
                                  aa_masses = aa_masses, 
                                  ntmass = ntmass, 
                                  ctmass = ctmass, 
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
                                  gen_ms2ions_base, 
                                  aa_masses = aa_masses, 
                                  ntmass = ntmass, 
                                  ctmass = ctmass, 
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
                                  gen_ms2ions_base, 
                                  aa_masses = aa_masses, 
                                  ntmass = ntmass, 
                                  ctmass = ctmass, 
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
#' (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-"
#' 
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @inheritParams ms2match_base
#' 
#' @seealso \link{bions_base}, \link{yions_base}.
#' 
#' @examples
#' \donttest{
#' # (2) "amods- tmod+ vnl- fnl-"
#' fixedmods <- c("TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("TMT6plex (N-term)", "Acetyl (Protein N-term)", "Oxidation (M)",
#'              "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses = aa_masses_all[[2]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] + 1.00727647 # + proton
#' }
#'
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#'
#' out <- gen_ms2ions_base(aa_seq, ms1_mass, aa_masses,
#'                         ntmass, ctmass, mod_indexes)
#'                         
#' # (1) "amods- tmod- vnl- fnl-"
#' fixedmods <- c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods <- c("Oxidation (M)", "Deamidated (N)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses <- aa_masses_all[[1]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] + 1.00727647
#' }
#'
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#'
#' out <- gen_ms2ions_base(aa_seq, ms1_mass, aa_masses,
#'                         ntmass, ctmass, mod_indexes)
#' 
#' }
gen_ms2ions_base <- function (aa_seq, ms1_mass = NULL, aa_masses, 
                              ntmass = NULL, ctmass = NULL, mod_indexes = NULL, 
                              type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                              maxn_sites_per_vmod = 3L, 
                              maxn_vmods_sitescombi_per_pep = 32L, 
                              digits = 4L) {
  
  aas <- aa_seq %>% str_split("", simplify = TRUE)
  aas2 <- aa_masses[aas]
  
  out <- ms2ions_by_type(aas2, ntmass, ctmass, type_ms2ions, digits)
  
  nm <- rep(0, length(aas)) %>% paste0(collapse = "")
  out <- list(out) %>% `names<-`(nm)
}



