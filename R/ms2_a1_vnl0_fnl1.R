#' Matching MS2 ions.
#' 
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
#' 
#' @param amods The attribute of \code{amods} from an \code{aa_masses}.
#' @param fmods_nl The attribute of \code{fmods_nl} from an \code{aa_masses}.
#' 
#' @rdname ms2match_base
#' @import purrr
#' @import parallel
#' @import dplyr
ms2match_a1_vnl0_fnl1 <- function (i, aa_masses, ntmod = NULL, ctmod = NULL, 
                                   ntmass = NULL, ctmass = NULL, amods, fmods_nl, 
                                   mod_indexes, mgf_path, out_path, 
                                   type_ms2ions = "by", maxn_vmods_per_pep = 5L, 
                                   maxn_sites_per_vmod = 3L, 
                                   maxn_vmods_sitescombi_per_pep = 32L, 
                                   minn_ms2 = 6L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                                   min_ms2mass = 110L, digits = 4L) {
  
  n_cores <- parallel::detectCores()
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
    cl, hms2_a1_vnl0_fnl1, 
    mgf_frames, theopeps, 
    MoreArgs = list(aa_masses = aa_masses, 
                    ntmod = ntmod, 
                    ctmod = ctmod, 
                    ntmass = ntmass, 
                    ctmass = ctmass, 
                    amods = amods, 
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
  
  parallel::stopCluster(cl)
  
  rm(list = c("mgf_frames", "theopeps"))
  gc()
  
  invisible(out)
}


#' Searches MGF frames.
#'
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
#' 
#' @inheritParams ms2match_a1_vnl0_fnl0
#' @rdname hms2_base
hms2_a1_vnl0_fnl1 <- function (mgf_frames, theopeps, aa_masses, 
                               ntmod = NULL, ctmod = NULL, 
                               ntmass, ctmass, 
                               amods, fmods_nl, mod_indexes, type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5L, 
                               maxn_sites_per_vmod = 3L, 
                               maxn_vmods_sitescombi_per_pep = 32L, 
                               minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                               min_ms2mass = 110L, digits = 4L) {
  
  res <- frames_adv_a1_vnl0_fnl1(mgf_frames = mgf_frames, 
                                 theopeps = theopeps, 
                                 aa_masses = aa_masses, 
                                 ntmod = ntmod, 
                                 ctmod = ctmod, 
                                 ntmass = ntmass, 
                                 ctmass = ctmass, 
                                 amods = amods, 
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
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
#'
#' @rdname frames_adv_base
#' @import purrr
frames_adv_a1_vnl0_fnl1 <- function (mgf_frames, theopeps, aa_masses, 
                                     ntmod = NULL, ctmod = NULL, 
                                     ntmass, ctmass, amods, fmods_nl, mod_indexes, 
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
                       gen_ms2ions_a1_vnl0_fnl1, 
                       aa_masses = aa_masses, 
                       ntmod = ntmod, 
                       ctmod = ctmod, 
                       ntmass = ntmass, 
                       ctmass = ctmass, 
                       amods = amods, 
                       fmods_nl = fmods_nl, 
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
                       gen_ms2ions_a1_vnl0_fnl1, 
                       aa_masses = aa_masses, 
                       ntmod = ntmod, 
                       ctmod = ctmod, 
                       ntmass = ntmass, 
                       ctmass = ctmass, 
                       amods = amods, 
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
    
    theos_af_ms2 <- map2(theopeps_af_ms1, 
                         theomasses_af_ms1, 
                         gen_ms2ions_a1_vnl0_fnl1, 
                         aa_masses = aa_masses, 
                         ntmod = ntmod, 
                         ctmod = ctmod, 
                         ntmass = ntmass, 
                         ctmass = ctmass, 
                         amods = amods, 
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
                           gen_ms2ions_a1_vnl0_fnl1, 
                           aa_masses = aa_masses, 
                           ntmod = ntmod, 
                           ctmod = ctmod, 
                           ntmass = ntmass, 
                           ctmass = ctmass, 
                           amods = amods, 
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
      
      theos_bf_ms2 <- map2(theopeps_bf_ms1, 
                           theomasses_bf_ms1, 
                           gen_ms2ions_a1_vnl0_fnl1, 
                           aa_masses = aa_masses, 
                           ntmod = ntmod, 
                           ctmod = ctmod, 
                           ntmass = ntmass, 
                           ctmass = ctmass, 
                           amods = amods, 
                           fmods_nl = fmods_nl, 
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
                           gen_ms2ions_a1_vnl0_fnl1, 
                           aa_masses = aa_masses, 
                           ntmod = ntmod, 
                           ctmod = ctmod, 
                           ntmass = ntmass, 
                           ctmass = ctmass, 
                           amods = amods, 
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
#' (11) "amods+ tmod- vnl- fnl+", (12) "amods+ tmod+ vnl- fnl+"
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
#' # (11) "amods+ tmod+ vnl+ fnl-"
#' fixedmods <- c("TMT6plex (K)", "Oxidation (M)", "dHex (S)")
#' varmods <- c("Deamidated (N)", "Acetyl (Protein N-term)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods,
#'                                add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses <- aa_masses_all[[2]]
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
#' fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
#'
#' aa_seq <- "HQGVMNVGMGQKMNS"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' 
#' ms1_mass <- ms1_masses$mass[[2]][2] # 2041.8958
#'
#' out <- gen_ms2ions_a1_vnl0_fnl1(aa_seq, ms1_mass, aa_masses,
#'                                 ntmod, ctmod,
#'                                 ntmass, ctmass, amods, fmods_nl, 
#'                                 mod_indexes)
#' 
#' # No "M", no "S"
#' aa_seq <- "HQGVNVGGQKN"
#' ms1_masses <- calc_monopeptide(aa_seq, fixedmods, varmods)
#' ms1_mass <- ms1_masses$mass[[2]][2] # 1367.6996
#' 
#' }
gen_ms2ions_a1_vnl0_fnl1 <- function (aa_seq, ms1_mass = NULL, aa_masses, 
                                      ntmod = NULL, ctmod = NULL, 
                                      ntmass = NULL, ctmass = NULL, 
                                      amods = NULL, fmods_nl, 
                                      mod_indexes = NULL, type_ms2ions = "by", 
                                      maxn_vmods_per_pep = 5L, 
                                      maxn_sites_per_vmod = 3L, 
                                      maxn_vmods_sitescombi_per_pep = 32L, 
                                      digits = 4L) {
  
  # (7, 8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
  # (no pep_seq dispatching by fmod residues -> possible no matched sites)
  sites <- names(fmods_nl)
  pattern <- paste(sites, collapse = "|")

  if (!grepl(pattern, aa_seq)) {
    out <- gen_ms2ions_a1_vnl0_fnl0(aa_seq = aa_seq, ms1_mass = ms1_mass, 
                                    aa_masses = aa_masses, 
                                    ntmod = ntmod, ctmod = ctmod, 
                                    ntmass = ntmass, ctmass = ctmass, 
                                    amods = amods, mod_indexes = mod_indexes, 
                                    type_ms2ions = type_ms2ions, 
                                    maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                    maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                    maxn_vmods_sitescombi_per_pep = 
                                      maxn_vmods_sitescombi_per_pep, 
                                    digits = digits)
    return(out)
  }
  
  # (11, 12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
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

  # NLs of fixedmods
  fnl_idxes <- which(aas %in% names(fmods_nl))
  
  # if (length(fnl_idxes) == 0L) return (NULL)
  
  fmods_combi <- aas[fnl_idxes]
  names(fmods_combi) <- fnl_idxes
  fnl_combi <- expand.grid(fmods_nl[fmods_combi])
  
  # go through each vmods_combi
  out <- purrr::map(vmods_combi, 
                    calc_ms2ions_a1_vnl0_fnl1, 
                    fnl_combi, fnl_idxes, aas2, aa_masses, ntmass, ctmass, 
                    type_ms2ions, digits = digits)
  
  out <- map2(out, vmods_combi, add_hexcodes_fnl2, length(aas), mod_indexes) %>% 
    flatten()
  
  invisible(out)
}


#' Calculates 
#' 
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' 
calc_ms2ions_a1_vnl0_fnl1 <- function (vmods_combi, fnl_combi, fnl_idxes, 
                                       aas2, aa_masses, 
                                       ntmass, ctmass, type_ms2ions, digits) {

  # updates amod masses
  delta_amod <- aa_masses[vmods_combi]
  amod_idxes <- as.numeric(names(vmods_combi))
  aas2[amod_idxes] <- aas2[amod_idxes] + delta_amod
  
  # updates fnl masses
  len <- nrow(fnl_combi)
  
  out <- vector("list", len)

  for (i in 1:len) {
    aas2_i <- aas2
    delta_nl <- unlist(fnl_combi[i, ], use.names = FALSE)
    aas2_i[fnl_idxes] <- aas2_i[fnl_idxes] - delta_nl
    out[[i]] <- ms2ions_by_type(aas2_i, ntmass, ctmass, type_ms2ions, digits)
  }

  invisible(out)
}


#' Adds hex codes (with variable NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @inheritParams add_hexcodes
add_hexcodes_fnl2 <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) {

  hex_mods = rep("0", len)
  hex_mods[as.numeric(names(vmods_combi))] <- mod_indexes[unlist(vmods_combi)]
  hex_mods <- paste0(hex_mods, collapse = "")
  
  # `(` for `vnl` and `[` for fnl
  names(ms2ions) <- paste0(hex_mods, " [", seq_along(ms2ions), "]")
  
  invisible(ms2ions)
}


