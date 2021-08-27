#' Matching MS2 ions.
#' 
#' (7) "amods+ tmod- vnl- fnl-", (8) "amods+ tmod+ vnl- fnl-"
#' 
#' @param amods The attribute of \code{amods} from an \code{aa_masses}.
#' @rdname ms2match_base
#' @import purrr
#' @import parallel
#' @import dplyr
ms2match_a1_vnl0_fnl0 <- function (i, aa_masses, ntmod = NULL, ctmod = NULL, 
                                   ntmass = NULL, ctmass = NULL, amods, 
                                   mod_indexes, mgf_path, out_path, 
                                   type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                                   maxn_sites_per_vmod = 3L, 
                                   maxn_vmods_sitescombi_per_pep = 32L, 
                                   minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                                   min_ms2mass = 110L, digits = 4L) {
  
  n_cores <- detect_cores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  parallel::clusterExport(cl, list("%fin%"), envir = environment(fastmatch::`%fin%`))
  parallel::clusterExport(cl, list("fmatch"), envir = environment(fastmatch::fmatch))
  parallel::clusterExport(cl, list("post_frame_adv"), envir = environment(proteoQ:::post_frame_adv))
  
  tempdata <- purge_search_space(i, aa_masses, mgf_path, n_cores, ppm_ms1)
  mgf_frames <- tempdata$mgf_frames
  theopeps <- tempdata$theopeps
  rm(list = c("tempdata"))
  
  if (length(mgf_frames) == 0L || length(theopeps) == 0L) return(NULL)
  
  out <- parallel::clusterMap(
    cl, hms2_a1_vnl0_fnl0, 
    mgf_frames, theopeps, 
    MoreArgs = list(aa_masses = aa_masses, 
                    ntmod = ntmod, 
                    ctmod = ctmod, 
                    ntmass = ntmass, 
                    ctmass = ctmass, 
                    amods = amods, 
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
  
  parallel::stopCluster(cl)
  
  rm(list = c("mgf_frames", "theopeps"))
  gc()
  
  invisible(out)
}


#' Searches MGF frames.
#'
#' (7) "amods+ tmod- vnl- fnl-", (8) "amods+ tmod+ vnl- fnl-"
#' 
#' @inheritParams ms2match_a1_vnl0_fnl0
#' @rdname hms2_base
hms2_a1_vnl0_fnl0 <- function (mgf_frames, theopeps, aa_masses, 
                               ntmod = NULL, ctmod = NULL, 
                               ntmass, ctmass, 
                               amods, mod_indexes, type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L, 
                               maxn_vmods_sitescombi_per_pep = 32L, 
                               minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                               min_ms2mass = 110L, digits = 4L) {
  
  # `res[[i]]` contains results for multiple mgfs within a frame
  # (the number of entries equals to the number of mgf frames)
  
  res <- frames_adv_a1_vnl0_fnl0(mgf_frames = mgf_frames, 
                                 theopeps = theopeps, 
                                 aa_masses = aa_masses, 
                                 ntmod = ntmod, 
                                 ctmod = ctmod, 
                                 ntmass = ntmass, 
                                 ctmass = ctmass, 
                                 amods = amods, 
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
#' (7) "amods+ tmod- vnl- fnl-", (8) "amods+ tmod+ vnl- fnl-"
#'
#' @rdname frames_adv_base
#' @import purrr
frames_adv_a1_vnl0_fnl0 <- function (mgf_frames, theopeps, aa_masses, 
                                     ntmod = NULL, ctmod = NULL, 
                                     ntmass, ctmass, amods, mod_indexes, 
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
  
  theos_bf_ms2 <- map2(theopeps_bf_ms1, 
                       theomasses_bf_ms1, 
                       gen_ms2ions_a1_vnl0_fnl0, 
                       aa_masses = aa_masses, 
                       ntmod = ntmod, 
                       ctmod = ctmod, 
                       ntmass = ntmass, 
                       ctmass = ctmass, 
                       amods = amods, 
                       mod_indexes = mod_indexes, 
                       type_ms2ions = type_ms2ions, 
                       maxn_vmods_per_pep = maxn_vmods_per_pep, 
                       maxn_sites_per_vmod = maxn_sites_per_vmod, 
                       maxn_vmods_sitescombi_per_pep = 
                         maxn_vmods_sitescombi_per_pep, 
                       digits = digits) %>% 
    `names<-`(theopeps_bf_ms1)
  
  theos_cr_ms2 <- map2(theopeps_cr_ms1, 
                       theomasses_cr_ms1, 
                       gen_ms2ions_a1_vnl0_fnl0, 
                       aa_masses = aa_masses, 
                       ntmod = ntmod, 
                       ctmod = ctmod, 
                       ntmass = ntmass, 
                       ctmass = ctmass, 
                       amods = amods, 
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
    
    theos_af_ms2 <- map2(theopeps_af_ms1, 
                         theomasses_af_ms1, 
                         gen_ms2ions_a1_vnl0_fnl0, 
                         aa_masses = aa_masses, 
                         ntmod = ntmod, 
                         ctmod = ctmod, 
                         ntmass = ntmass, 
                         ctmass = ctmass, 
                         amods = amods, 
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
    
    ## `theos_xx_ms2` may contain empty entries and handled in `find_ms2_bypep`: 
    #   e.g. the matched one is after `maxn_vmods_sitescombi_per_pep` 
    #   and never get matched.
    
    out[[i]] <- map2(exptmasses_ms1, exptmoverzs_ms2, 
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
      
      theos_cr_ms2 <- map2(theopeps_cr_ms1, 
                           theomasses_cr_ms1, 
                           gen_ms2ions_a1_vnl0_fnl0, 
                           aa_masses = aa_masses, 
                           ntmod = ntmod, 
                           ctmod = ctmod, 
                           ntmass = ntmass, 
                           ctmass = ctmass, 
                           amods = amods, 
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
      
      theos_bf_ms2 <- map2(theopeps_bf_ms1, 
                           theomasses_bf_ms1, 
                           gen_ms2ions_a1_vnl0_fnl0, 
                           aa_masses = aa_masses, 
                           ntmod = ntmod, 
                           ctmod = ctmod, 
                           ntmass = ntmass, 
                           ctmass = ctmass, 
                           amods = amods, 
                           mod_indexes = mod_indexes, 
                           type_ms2ions = type_ms2ions, 
                           maxn_vmods_per_pep = maxn_vmods_per_pep, 
                           maxn_sites_per_vmod = maxn_sites_per_vmod, 
                           maxn_vmods_sitescombi_per_pep = 
                             maxn_vmods_sitescombi_per_pep, 
                           digits = digits) %>% 
        `names<-`(theopeps_bf_ms1)
      
      theos_cr_ms2 <- map2(theopeps_cr_ms1, 
                           theomasses_cr_ms1, 
                           gen_ms2ions_a1_vnl0_fnl0, 
                           aa_masses = aa_masses, 
                           ntmod = ntmod, 
                           ctmod = ctmod, 
                           ntmass = ntmass, 
                           ctmass = ctmass, 
                           amods = amods, 
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
#' (7) "amods+ tmod- vnl- fnl-", (8) "amods+ tmod+ vnl- fnl-"
#' 
#' @rdname gen_ms2ions_base
#' 
#' 
#' @param amods The attribute \code{amods} from a \code{aa_masses}.
#' @param ntmod The attribute \code{ntmod} from a \code{aa_masses} (for MS1
#'   calculations).
#' @param ctmod The attribute \code{ctmod} from a \code{aa_masses} (for MS1
#'   calculations).
#'  @rdname gen_ms2ions_base
#'
#' @examples
#' \donttest{
#' # (8) "amods+ tmod+ vnl- fnl-"
#' fixedmods <- c("TMT6plex (K)")
#' varmods <- c("Deamidated (N)", "Carbamidomethyl (S)",
#'              "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses <- aa_masses_all[[8]]
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
#' amods <- attr(aa_masses, "amods", exact = TRUE)
#'
#' aa_seq <- "MHQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide("MHQGVMNVGMGQKMNS",
#'                                fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[8]][2] # 2077.9256
#'
#' out <- gen_ms2ions_a1_vnl0_fnl0(aa_seq, ms1_mass, aa_masses,
#'                                 ntmod, ctmod,
#'                                 ntmass, ctmass, amods,
#'                                 mod_indexes)
#'
#' }
gen_ms2ions_a1_vnl0_fnl0 <- function (aa_seq, ms1_mass = NULL, aa_masses, 
                                      ntmod = NULL, ctmod = NULL, 
                                      ntmass = NULL, ctmass = NULL, amods = NULL, 
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_vmods_sitescombi_per_pep = 32L, 
                                      digits = 4L) {

  aas <- aa_seq %>% stringr::str_split("", simplify = TRUE)
  aas2 <- aa_masses[aas]

  vmods_combi <- combi_mvmods2(amods = amods, 
                               aas = aas, 
                               aa_masses = aa_masses, 
                               maxn_vmods_per_pep = maxn_vmods_per_pep, 
                               maxn_sites_per_vmod = maxn_sites_per_vmod, 
                               maxn_vmods_sitescombi_per_pep = 
                                 maxn_vmods_sitescombi_per_pep, 
                               digits = digits) %>% 
    find_intercombi_p2(maxn_vmods_sitescombi_per_pep)
  
  # filtered by MS1 masses
  if (!(purrr::is_empty(vmods_combi) || is.null(ms1_mass))) {
    idxes <- purrr::map_lgl(vmods_combi, check_ms1_mass_vmods2, aas2, aa_masses, 
                            ntmod, ctmod, ms1_mass)
    
    vmods_combi <- vmods_combi[idxes]
    rm(list = c("idxes"))
  }

  # ---
  out <- purrr::map(vmods_combi, 
                    calc_ms2ions_a1_vnl0_fnl0, 
                    aas2, aa_masses, ntmass, ctmass, 
                    type_ms2ions, digits = digits)
  
  out <- add_hexcodes(out, vmods_combi, length(aas), mod_indexes)

  invisible(out)
}


#' Calculates 
#' 
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' 
calc_ms2ions_a1_vnl0_fnl0 <- function (vmods_combi, aas2, aa_masses, 
                                       ntmass, ctmass, type_ms2ions, digits) {

  # find delta
  delta <- aa_masses[vmods_combi]
  idxes <- as.numeric(names(vmods_combi))
  aas2[idxes] <- aas2[idxes] + delta
  
  ms2ions_by_type(aas2, ntmass, ctmass, type_ms2ions, digits)
}


#' Check the MS1 mass for proceeding with MS2 matches.
#'
#' Maybe missed if the "matched" occurred after `maxn_vmods_sitescombi_per_pep`.
#'
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @param tol The tolerance in mass.
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @importFrom purrr is_empty
check_ms1_mass_vmods2 <- function (vmods_combi, aas2, aa_masses, ntmod, ctmod, 
                                   ms1_mass, tol = 1e-3) {

  # bare + 18.010565 + terminals + anywhere
  
  bare <- aas2 %>% 
    sum() %>% 
    `+`(18.010565) 
  
  # No need of is_empty(ntmod) && is_empty(ctmod)
  
  if (purrr::is_empty(ntmod) && purrr::is_empty(ctmod)) {
    delta <- 0
  } else if (!(purrr::is_empty(ntmod) || purrr::is_empty(ctmod))) {
    delta <- aa_masses[names(ntmod)] + aa_masses[names(ctmod)]
  } else if (!purrr::is_empty(ntmod)) {
    delta <- aa_masses[names(ntmod)]
  } else if (!purrr::is_empty(ctmod)) {
    delta <- aa_masses[names(ctmod)]
  }
  
  vmass <- sum(aa_masses[vmods_combi])
  
  ok_mass <- bare + delta + vmass
  
  abs(ms1_mass - ok_mass) <= tol
}


#' Finds the combinations of variable modifications (multiple sites)
#'
#' For all the \code{Anywhere} modifications specified in \code{amods}.
#'
#' @param amods \code{Anywhere} variable modifications.
#' @param famods \code{Anywhere} fixed modifications.
#' @param fntmod The parsed \emph{fixed} \code{N-term} modifications from
#'   \code{aa_masses}.
#' @param fctmod The parsed \emph{fixed} \code{C-term} modifications from
#'   \code{aa_masses}.
#' @inheritParams calc_monopep
#' @inheritParams mcalc_monopep
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @return Lists by residues in \code{amods}.
combi_mvmods2 <- function (amods, 
                           aas, 
                           aa_masses, 
                           maxn_vmods_per_pep = 5L, 
                           maxn_sites_per_vmod = 3L, 
                           maxn_vmods_sitescombi_per_pep = 32L, 
                           digits = 4L) {
  
  ## Split by residues
  # Oxidation (M) Carbamidomethyl (M) 
  # "M"                 "M" 
  # 
  # $S
  # dHex (S) 
  # "S" 
  
  residue_mods <- unlist(amods, use.names = FALSE) %>% 
    `names<-`(names(amods)) %>% 
    split(., .)

  purrr::map(residue_mods, ~ combi_vmods2(
    aas, .x, 
    aa_masses, 
    maxn_vmods_per_pep, 
    maxn_sites_per_vmod, 
    maxn_vmods_sitescombi_per_pep, 
    digits
  ))
}


#' The combinations of variable modifications (single site)
#'
#' @param residue_mods A residue with \code{Anywhere} modification(s).
#' @inheritParams combi_mvmods
#' @import purrr
#' @importFrom stringr str_locate_all
combi_vmods2 <- function (aas, 
                          residue_mods, 
                          aa_masses, 
                          maxn_vmods_per_pep = 5L, 
                          maxn_sites_per_vmod = 3L, 
                          maxn_vmods_sitescombi_per_pep = 32L, 
                          digits = 4L) {
  
  ##################################################################
  # values: n (modifications)
  # names: p (positions)
  # 
  # n = LETTERS[1:2]; p = c(1, 3, 16)
  # n = c("Carbamidomethyl (M)",  "Oxidation (M)"); p = c(1, 3, 16)
  # 2*3, 4*3, 8*1
  # l = length(p)
  # n^1 * combn(p, 1) + n^2 * combn(p, 2) + ... + n^l * combn(p, l)
  # 
  ##################################################################
  
  ##################################################################
  # !!! Danger !!!
  # combn(3, 1) is combn(1:3, 1) not combn("3", 1)
  ##################################################################
  
  residue <- residue_mods[[1]]
  
  n <- names(residue_mods)
  p <- which(aas == residue)
  
  # (1) btw Anywhere "M" and "Acetyl N-term" where "M" on the "N-term"
  # MFGMFNVSMR cannot have three `Oxidation (M)` and `Acetyl (N-term)`
  # (2) the same for fixed terminal mod: `TMT6plex (N-term)` 
  
  # p <- check_tmod_p(aas, residue, p, ntmod, ctmod)
  # p <- check_tmod_p(aas, residue, p, fntmod, fctmod)
  
  len_n <- length(n)
  len_p <- length(p)
  
  if (len_n > len_p) {
    return(NULL)
  }
  
  if (len_p == 1L) {
    names(n) <- p
    return(n)
  }
  
  # --- combinations ---
  if (len_n == 1L) { # "Oxidation (M)"
    len <- length(p)
    out <- vector("list", len)
    
    for (m in seq_along(p)) {
      ns <- rep(n, m)
      ps <- combn(p, m)
      
      ncol <- ncol(ps)
      ns <- rep(list(ns), ncol)
      
      for (i in 1:ncol) {
        names(ns[[i]]) <- ps[, i]
      }
      
      out[[m]] <- ns
    }
  } else { # "Oxidation (M)" and "Carbamidomethyl (M)"
    ns <- purrr::map(seq_along(p), ~ {
      expand.grid(rep(list(n), length(p[1:.x])), KEEP.OUT.ATTRS = FALSE, 
                  stringsAsFactors = FALSE)
    }, n, p) 
    
    ps <- purrr::map(seq_along(p), ~ {
      combn(as.character(p), .x)
    }, p)
    
    # out <- mvmcombC(ns, ps)
    
    out <- purrr::map2(ns, ps, ~ {
      lns <- nrow(.x)
      lps <- ncol(.y)
      np <- vector("list", lns * lps)
      
      k <- 1
      for (i in seq_len(lns)) {
        for (j in seq_len(lps)) {
          x <- .x[i, ] 
          names(x) <- .y[, j]
          np[[k]] <- x
          
          k <- k + 1
        }
      }
      
      # otherwise is data.frame
      # use names (positions) -> TRUE
      np <- purrr::map(np, unlist, use.names = TRUE)
    }) 
  }
  
  # ---
  out <- out %>% 
    purrr::flatten() 
  
  rows <- purrr::map_lgl(out, ~ length(.x) > maxn_vmods_per_pep)
  out <- out[!rows]
  
  maxn_vmod <- out %>% 
    # map(unlist) %>% 
    purrr::map(table) %>% 
    purrr::map(max)
  rows <- maxn_vmod > maxn_sites_per_vmod
  out <- out[!rows]
  
  if (!is_empty(out)) {
    out <- out %>% 
      `[`(1:min(length(.), maxn_vmods_sitescombi_per_pep))
  }
  
  invisible(out)
}


#' Finds the combinations of positions and sites across residues.
#' 
#' For multiple residues (each residue one to multiple modifications).
#' 
#' @param intra_combis The results from \link{unique_mvmods}.
#' @inheritParams calc_ms2ions
#' @importFrom purrr is_empty map map_lgl flatten 
#' @seealso \link{find_intercombi}
find_intercombi_p2 <- function (intra_combis, maxn_vmods_sitescombi_per_pep) {
  
  if (is_empty(intra_combis) || any(map_lgl(intra_combis, is_empty))) { # scalar or list
    v_out <- list() 
  } else if (length(intra_combis) == 1L) { # M, one to multiple positions; Oxidation and/or Carbamidomethyl
    if (length(intra_combis[[1]]) == 1L) { # 2: "Oxidation (M)"
      v_out <- intra_combis
    } else { # 2: "Oxidation (M)"; 3: "Oxidation (M)"; 2: "Oxidation (M)", 3: "Oxidation (M)"; ... Carbamidomethyl
      v_out <- flatten(intra_combis)
    }
  } else { # M, N
    p_combis <- map(intra_combis, ~ {
      x <- .x
      
      if (length(x) > 1L) {
        map(x, names)
      } else {
        names(x)
      }
    }) 
    
    v_combis <- map(intra_combis, ~ {
      x <- .x
      map(x, reduce, `c`)
    })
    
    # max_combi <- map(p_combis, ~ min(length(.x), maxn_vmods_sitescombi_per_pep))
    # p_combis <- map2(p_combis, max_combi, ~ .x[seq_len(.y)])
    # v_combis <- map2(v_combis, max_combi, ~ .x[seq_len(.y)])

    p_combis <- expand.grid(p_combis, KEEP.OUT.ATTRS = FALSE, 
                            stringsAsFactors = FALSE)
    v_combis <- expand.grid(v_combis, KEEP.OUT.ATTRS = FALSE, 
                            stringsAsFactors = FALSE)
    
    nrow <- min(nrow(v_combis), maxn_vmods_sitescombi_per_pep)
    
    p_out <- v_out <- vector("list", nrow)
    
    for (i in seq_len(nrow)) {
      p_out[[i]] <- unlist(p_combis[i, ], use.names = FALSE)
      v_out[[i]] <- unlist(v_combis[i, ], use.names = FALSE)
      
      names(v_out[[i]]) <- p_out[[i]]
      v_out[[i]] <- v_out[[i]][order(as.numeric(names(v_out[[i]])))]
    }
    
    # use names (positions) -> TRUE
    # v_out <- map(v_out, unlist, use.names = TRUE)
  }
  
  ## use names (positions)
  # v_out <- map(v_out, unlist, use.names = TRUE)
  
  invisible(v_out)
}



