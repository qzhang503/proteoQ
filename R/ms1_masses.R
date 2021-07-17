#' Helper of \link{calcms1mass_noterm}.
#'
#' By individual proteins.
#'
#' @param fnl_combi A data.frame of combinations of neutral losses for fixed
#'   modifications. Each row corresponds to a set of neutral loss. The first row
#'   corresponds to the combination without NLs (all zeros).
#' @param prot_peps Lists of named peptide masses under a protein.
#' @inheritParams calc_monopep
ms1_a0_fnl1_byprot <- function (prot_peps, fnl_combi, aa_masses, digits = 4L) {
  
  out <- purrr::map2(prot_peps, names(prot_peps), ms1_a0_fnl1_bypep, 
                     fnl_combi = fnl_combi, 
                     aa_masses = aa_masses, 
                     digits = digits) 
  
  nms <- purrr::imap(out, ~ rep(.y, length(.x))) %>% 
    unlist(recursive = FALSE, use.names = FALSE)
  
  out <- out %>% 
    unlist(recursive = FALSE, use.names = FALSE) %>% 
    `names<-`(nms)
}


#' Helper of \link{calcms1mass_noterm_byprot}.
#' 
#' By individual peptides.
#' 
#' @param mass The mass of a peptide.
#' @inheritParams ms1_a0_fnl1_byprot
#' @inheritParams calc_monopep
#' @importFrom stringr str_split
ms1_a0_fnl1_bypep <- function (mass, aa_seq, fnl_combi, aa_masses, digits = 4L) {

  aas <- str_split::str_split(aa_seq, "", simplify = TRUE)

  delta <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)
  
  round(mass - delta, digits = digits)
}


#' Helper of peptide-mass calculation..
#'
#' (5) "amods- tmod+ vnl- fnl+"; (6) "amods- tmod- vnl- fnl+".
#'
#' The calculation goes through the rows in \code{fnl_combi}. 
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams ms1_a0_fnl1_byprot
#' @return A numeric vector
delta_ms1_a0_fnl1 <- function (fnl_combi, aas, aa_masses) {
  
  # not an option but always: include_insource_nl = TRUE
  # bring back the option later...
  
  nms <- colnames(fnl_combi)
  oks <- aas[aas %in% nms]
  
  if (is_empty(oks)) return (0L)
  
  len <- nrow(fnl_combi)
  out <- vector("numeric", len)
  out[[1]] <- 0L

  for (i in 2:len) {
    aa_masses[nms] <- unlist(fnl_combi[i, ])
    oks <- aas[aas %in% nms]
    out[[i]] <- sum(aa_masses[oks])
  }
  
  out <- unique(out)
}


#' Helper by individual proteins.
#'
#' (7) "amods+ tmod- vnl- fnl-"; (8) "amods+ tmod-+ vnl- fnl-".
#'
#' @param amods \code{Anywhere} variable modifications.
#' @inheritParams ms1_a0_fnl1_byprot
#' @inheritParams calc_monopep
#' 
#' @examples 
#' \donttest{
#' m0 <- calc_monopeptide("HQGVMCNVGMGQKMNSC", NULL, NULL)
#' stopifnot(unlist(m0$mass, use.names = FALSE) - 1822.7405 < 1e-4)
#' 
#' m1 <- calc_monopeptide("HQGVMCNVGMGQKMNSC", "TMT6plex (N-term)", NULL)
#' stopifnot(unlist(m1$mass, use.names = FALSE) - 2051.9035 < 1e-4)
#' 
#' # (7) "amods+ tmod- vnl- fnl-"
#' aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)"), 
#'                                varmods = c("Deamidated (N)", 
#'                                            "Carbamidomethyl (C)"), 
#'                                add_varmasses = FALSE, 
#'                                add_nlmasses = FALSE)
#' 
#' pep <- c("HQGVMCNVGMGQKMNSC" = 2051.90346)
#' amods <- list(`Deamidated (N)` = c("Anywhere" = "N"), 
#'               `Carbamidomethyl (C)` = c("Anywhere" = "C"))
#' 
#' x <- ms1_a1_vnl0_fnl0_byprot(pep, amods, aa_masses_all[[4]])
#' 
#' stopifnot(x[[1]] - 2109.9089 < 1e-4, 
#'           x[[2]] - 2166.9304 < 1e-4, 
#'           x[[3]] - 2110.8930 < 1e-4, 
#'           x[[4]] - 2167.9144 < 1e-4)
#'           
#' 
#' # (8) "amods+ tmod+ vnl- fnl-"
#' aa_masses_all <- calc_aamasses(fixedmods = "TMT6plex (K)", 
#'                                varmods = c("Deamidated (N)", 
#'                                            "Carbamidomethyl (S)", 
#'                                            "Acetyl (Protein N-term)"), 
#'                                add_varmasses = FALSE, 
#'                                add_nlmasses = FALSE)
#' 
#' pep <- c("HQGVMNVGMGQKSMNS" = 1932.9171) # + TMT6plex (K)
#' amods <- list(`Deamidated (N)` = c("Anywhere" = "N"), 
#'               `Carbamidomethyl (S)` = c("Anywhere" = "S"))
#'               
#' x <- ms1_a1_vnl0_fnl0_byprot(pep, amods, aa_masses_all[[8]])
#' 
#' stopifnot(x[[1]] - 1990.9226 < 1e-4, 
#'           x[[2]] - 1991.9066 < 1e-4, 
#'           x[[3]] - 2047.9440 < 1e-4, 
#'           x[[4]] - 2048.9281 < 1e-4)
#' 
#' }
ms1_a1_vnl0_fnl0_byprot <- function (prot_peps, amods, aa_masses, 
                                     vmods_nl = NULL, fmods_nl = NULL, 
                                     include_insource_nl = FALSE, 
                                     maxn_vmods_per_pep = 5L,
                                     maxn_sites_per_vmod = 3L, 
                                     digits = 4L) {

  out <- purrr::map2(prot_peps, names(prot_peps), ms1_a1_vnl0_fnl0_bypep, 
                     amods = amods, aa_masses = aa_masses, 
                     vmods_nl = vmods_nl, fmods_nl = fmods_nl, 
                     include_insource_nl = include_insource_nl, 
                     maxn_vmods_per_pep = maxn_vmods_per_pep,
                     maxn_sites_per_vmod = maxn_sites_per_vmod, 
                     digits = digits) 

  nms <- purrr::imap(out, ~ rep(.y, length(.x))) %>% 
    unlist(recursive = FALSE, use.names = FALSE)
  
  out <- out %>% 
    unlist(recursive = FALSE, use.names = FALSE) %>% 
    `names<-`(nms)
}


#' Helper by individual peptides.
#' 
#' (7) "amods+ tmod- vnl- fnl-"; (8) "amods+ tmod-+ vnl- fnl-".
#' 
#' @param mass The mass of a peptide.
#' @inheritParams ms1_a1_vnl0_fnl0_byprot
#' @inheritParams calc_monopep
#' @importFrom stringr str_split
ms1_a1_vnl0_fnl0_bypep <- function (mass, aa_seq, amods, aa_masses, 
                                    vmods_nl = NULL, fmods_nl = NULL, 
                                    include_insource_nl = FALSE, 
                                    maxn_vmods_per_pep = 5L,
                                    maxn_sites_per_vmod = 3L, 
                                    digits = 4L) {
  
  aas <- stringr::str_split(aa_seq, "", simplify = TRUE)
  
  vmods_combi <- unique_mvmods(amods = amods, ntmod = NULL, ctmod = NULL, 
                               aa_masses = aa_masses, aas = aas, 
                               maxn_vmods_per_pep = maxn_vmods_per_pep, 
                               maxn_sites_per_vmod = maxn_sites_per_vmod, 
                               digits = digits) %>% 
    find_intercombi()
  
  deltas <- purrr::map(vmods_combi, ~ sum(aa_masses[.x]))
  
  out <- purrr::map_dbl(deltas, ~ round(mass + .x, digits = digits))
  
  if (include_insource_nl) {
    if (length(vmods_nl) > 0L) {
      vnl_combi <- purrr::map(vmods_combi, ~ expand.grid(vmods_nl[.x]))
      deltas_vnls <- purrr::map(vnl_combi, ~ unique(rowSums(.x)))
      
      out <- purrr::map2(out, deltas_vnls, `-`)
    }
    
    if (length(fmods_nl) > 0L) {
      fnl_combi <- expand.grid(fmods_nl)
      deltas_fnls <- delta_ms1_a0_fnl1(fnl_combi, aas, aa_masses)
      
      out <- purrr::map2(out, deltas_fnls, `-`)
    }
  }

  invisible(out)
}




