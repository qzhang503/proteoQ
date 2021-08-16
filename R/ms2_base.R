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
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz". The default is "by" for b-
#'   and y-ions.
#' @inheritParams matchMS
#' @inheritParams ms2match
#' @import purrr
#' @import parallel
ms2match_base <- function (i, aa_masses, ntmass, ctmass, mod_indexes, 
                           mgf_path, out_path, 
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
  gc()
  
  if (length(mgf_frames) == 0L || length(theopeps) == 0L) return(NULL)
  
  out <- parallel::clusterMap(
    cl, hms2_base, 
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
                       min_ms2mass = 110L, digits = 4L) {

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
                         min_ms2mass = min_ms2mass, 
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


#' Matches by indexes
#' 
#' @examples 
#' \donttest{
#' expt_ms2 <- 
#'   c(1628,3179,7677,9129,13950,14640,18571,19201,19205,19830,19833,20454,
#'     20457,21030,21073,21077,21687,24644,25232,37146,42042,43910,43920,44811,
#'     44824,45298,47494,55080,55901,56677,59014,66693,72396,72402,72720,73043,
#'     82411,91067,91838,93101,95572,98301,98665,100270,102081,102305,102744,106013,
#'     107998,108102,113713,113898,115045,115140,117669,119131,120730,123859,124029,124200,
#'     126199,126208,126610,126693,126775,128157,129447,129603,132396,135402,135475,138158,
#'     140397,141566,141634,141702,142183,142580,144189,147799,147926,148678,148860,149911,
#'     149973,153047,155607,158520,158631,162612,162717,163346,169537,170401,171249,171344,
#'     178012,178620,181980,188455)
#' 
#' theo_ms2 <- 
#'   c(-26231,62754,105787,129278,151731,161552,174924,184489,196534,204867,212917,219771,
#'     236270,106013,129447,148679,163242,178619,187776,197630,203310,212976,219825,227451,
#'     234026,237018)
#' cr <- which(expt_ms2 %fin% theo_ms2)
#' pr <- which((expt_ms2-1) %fin% theo_ms2)
#' af <- which((expt_ms2+1) %fin% theo_ms2)
#' c(cr, pr, af)
#' 
#' }
search_mgf2 <- function (expt_mass_ms1, expt_moverz_ms2, 
                         theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
                         theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                         minn_ms2 = 7L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                         min_ms2mass = 110L) {
  
  # --- subsets from the `before` and the `after` by MS1 mass tolerance 
  ms1_range <- find_mass_error_range(expt_mass_ms1, ppm_ms1) # 3.1 us
  bf_allowed <- which(theomasses_bf_ms1 >= ms1_range[1]) # 2 us
  af_allowed <- which(theomasses_af_ms1 <= ms1_range[2]) # 1.9 us
  
  # not used but kept for tidiness
  theomasses_bf_ms1 <- theomasses_bf_ms1[bf_allowed]
  theomasses_af_ms1 <- theomasses_af_ms1[af_allowed]
  
  theos_bf_ms2 <- theos_bf_ms2[bf_allowed]
  theos_af_ms2 <- theos_af_ms2[af_allowed]
  
  # --- find MS2 matches ---
  if (is_empty(theos_bf_ms2)) {
    x_bf <- theos_bf_ms2
  } else {
    x_bf <- purrr::map(theos_bf_ms2, find_ms2_bypep, expt_moverz_ms2, ppm_ms2, 
                       min_ms2mass)
  }
  
  if (is_empty(theos_cr_ms2)) {
    x_cr <- theos_cr_ms2
  } else {
    x_cr <- purrr::map(theos_cr_ms2, find_ms2_bypep, expt_moverz_ms2, ppm_ms2, 
                       min_ms2mass)
  }
  
  if (is_empty(theos_af_ms2)) {
    x_af <- theos_af_ms2
  } else {
    x_af <- purrr::map(theos_af_ms2, find_ms2_bypep, expt_moverz_ms2, ppm_ms2, 
                       min_ms2mass)
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
  theomasses_ms1 <- c(theomasses_bf_ms1, 
                      theomasses_cr_ms1, 
                      theomasses_af_ms1)
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
#' Run time about 41 us.
#'
#' @param expts Numeric vector; one series of experimental MS2s.
#' @param theos Numeric vector; one to multiple series of theoretical MS2s.
#' @inheritParams search_mgf_frames
#' @importFrom purrr map
#' @importFrom fastmatch fmatch %fin% 
#' @examples 
#' \donttest{
#' # experimental `322.18704` fit to both b- and y- ions
#' expts <- c(101.07140,102.05540,107.04956,110.07165,111.07500,
#'            112.08729,113.07130,115.08693,116.07092,120.08105,
#'            121.08428,126.12794,127.12494,127.13121,128.12825,
#'            128.13452,129.13164,129.13789,130.06493,130.09772,
#'            130.13495,130.14124,131.13831,132.14160,134.07635,
#'            136.07573,157.10826,158.09238,159.09117,170.12067,
#'            173.12846,173.14980,175.11896,176.12245,176.15956,
#'            184.11806,186.15303,188.15977,190.09743,193.63914,
#'            207.11292,210.12289,229.16670,230.17030,231.17410,
#'            235.10779,240.14383,248.18073,262.15305,265.67874,
#'            273.21210,285.13416,301.20700,305.16055,312.17740,
#'            314.69138,322.18704,365.23856,369.24496,371.70316,
#'            374.18283,376.27573,392.19308,393.23337,394.23743,
#'            399.20920,400.27567,400.72501,401.22650,401.27139,
#'            401.58405,401.72778,409.21918,410.22241,423.24564,
#'            433.29709,452.27066,462.25473,480.26578,481.26828,
#'            498.27670,530.35126,572.28278,573.28625,599.33899,
#'            600.34039,609.32202,626.29218,627.29303,627.33948,
#'            628.33417,629.30463,630.30219,643.31946,644.31763,
#'            644.35913,645.32520,646.32880,647.32825,648.32892)
#' 
#' theos <- c(114.0913,251.1503,322.1874,423.2350,480.2565,
#'            627.3249,783.4260,175.1190,322.1874,379.2088,
#'            480.2565,551.2936,688.3525,801.4366)
#' 
#' nms <- unlist(stringr::str_split(pep, ""))
#' 
#' names(theos) <- c(nms, rev(nms))
#' theos <- list(`0000000` = theos)
#' 
#' x <- find_ms2_bypep(theos, expts)
#' }
find_ms2_bypep <- function (theos, expts, ppm_ms2 = 25L, min_ms2mass = 110L) {
  
  # `theos` may be empty: 
  #   e.g. the matched one is after `maxn_vmods_sitescombi_per_pep` 
  #   and never get matched.
  
  len <- length(theos)
  
  if (len == 0L) return(list(theo = NULL, expt = NULL))
  
  if (!is.list(theos)) theos <- list(theos)
  
  ex <- find_ms1_interval(expts, from = min_ms2mass, ppm = ppm_ms2)
  
  # 100 us slower compared to direct for loops 
  # (1.2 ms vs 1.1 ms)
  # out <- map(theos, find_ms2_bycombi, expts, ex, ppm_ms2, min_ms2mass)
  
  # ---
  out <- vector("list", len)
  
  for (i in 1:len) {
    theos_i <- theos[[i]]

    th_i <- find_ms1_interval(theos_i, from = min_ms2mass, ppm = ppm_ms2)
    ps <- fuzzy_match_one(th_i, ex) # 5 us
    
    if (sum(ps) == 0L) {
      es <- rep(NA, length(theos_i))
      out[[i]] <- list(theo = theos_i, expt = es)
      # out[[i]] <- list(theo = NA, expt = NA)
    } else {
      # first filled by theoretical values
      es <- theos_i # .1 us
      es[!ps] <- NA # .9 us
      
      # separate b and y matches (to handled double-dipping btw b and y)
      lth <- length(ps)
      mid <- lth/2L
      
      bps <- fuzzy_match_one(ex, th_i[1:mid])
      yps <- fuzzy_match_one(ex, th_i[(mid+1L):lth])
      
      es[ps] <- c(expts[bps], expts[yps])
      out[[i]] <- list(theo = theos_i, expt = es)
    }
  }
  
  names(out) <- names(theos)
  
  out
}


#' Fuzzy matches with a +/-1 window.
#' 
#' @param x A vector to be matched.
#' @param y A vector to be matched against.
#' @importFrom fastmatch fmatch %fin% 
fuzzy_match_one <- function (x, y) {
  
  mi <- x %fin% y # 1.1 us
  bf <- (x - 1L) %fin% y # 1.7 us
  af <- (x + 1L) %fin% y # 1.6 us
  
  ps <- mi | bf | af # .5 us
}


#' Helper of \link{find_ms2_bypep}
#' 
#' @inheritParams find_ms2_bypep
#' @param ex Indexed \code{expts}.
find_ms2_bycombi <- function (theos, expts, ex, ppm_ms2, min_ms2mass) {
  
  th <- find_ms1_interval(theos, from = min_ms2mass, ppm = ppm_ms2)
  ps <- fuzzy_match_one(th, ex)
  
  if (sum(ps) == 0L) {
    es <- rep(NA, length(theos))
    out <- list(theo = theos, expt = es)
  } else {
    es <- theos
    es[!ps] <- NA
    
    lth <- length(ps)
    mid <- lth/2L
    
    bps <- fuzzy_match_one(ex, th[1:mid])
    yps <- fuzzy_match_one(ex, th[(mid+1L):lth])
    
    es[ps] <- c(expts[bps], expts[yps])
    out <- list(theo = theos, expt = es)
  }

  out
}

