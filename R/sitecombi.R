#' Helper: switches among ion types for calculating MS2 masses.
#'
#' @param aas2 A sequence of amino-acid residues with \emph{masses}. Residues
#'   are in names and masses in values (note that argument \code{aas}
#'   corresponds to residues without masses).
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_by_iontype <- function (aas2, aa_masses, ntmod, ctmod, type_ms2ions, 
                                digits) {
  switch(type_ms2ions, 
         by = calc_byions(ntmod, ctmod, aas2, aa_masses, digits), 
         cz = calc_czions(ntmod, ctmod, aas2, aa_masses, digits), 
         ax = calc_axions(ntmod, ctmod, aas2, aa_masses, digits), 
         # by2 = calc_by2ions(ntmod, ctmod, aas2, aa_masses, digits), 
         # bystar = calc_bystarions(ntmod, ctmod, aas2, aa_masses, digits), 
         # bystar2 = calc_bystar2ions(ntmod, ctmod, aas2, aa_masses, digits), 
         # by0 = calc_by0ions(ntmod, ctmod, aas2, aa_masses, digits), 
         # by02 = calc_by02ions(ntmod, ctmod, aas2, aa_masses, digits), 
         # cz2 = calc_cz2ions(ntmod, ctmod, aas2, aa_masses, digits), 
         # ax2 = calc_ax2ions(ntmod, ctmod, aas2, aa_masses, digits), 
         stop("Unknown type.", call. = FALSE))
}


#' Check the MS1 mass for proceeding with MS2 matches.
#' 
#' @param tol The tolerance in mass.
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
check_ms1_mass_vmods <- function (vmods_combi, aas, aa_masses, ntmod, ctmod, ms1_mass, 
                            tol = 2e-5) {
  aas[as.numeric(names(vmods_combi))] <- vmods_combi
  ok_mass <- calcpepmass_nl0(aas, aa_masses, ntmod, ctmod)
  
  abs(ms1_mass - ok_mass) <= tol
}


#' Helper: masses of a series of MS2 ions.
#'
#' Updates \code{aas} (since `amods+`) and calculates masses for moduels (7)
#' "amods+ tmod- vnl- fnl-" and (8) "amods+ tmod+ vnl- fnl-".
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
mcalcms2_a1_vnl0 <- function (vmods_combi, ntmod, ctmod, aas, aa_masses, 
                              type_ms2ions, digits = 5) {
  
  if (!is.null(vmods_combi)) {
    aas[as.numeric(names(vmods_combi))] <- vmods_combi
  }
  
  calcms2_by_iontype(aa_masses[aas], aa_masses, ntmod, ctmod, type_ms2ions, 
                     digits)
}


#' Helper: masses of a series of MS2 ions.
#'
#' Updates `aas` and calculates masses for moduels (9) "amods+ tmod- vnl+ fnl-"
#' and (10) "amods+ tmod+ vnl+ fnl-".
#'
#' It first updates the names of residues in \code{aas} (since `amods+`), then
#' runs column-wisely through the table of \code{nl_combi} (variable) for
#' corresponding masses of residues.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
mcalcms2_a1_vnl1 <- function (vmods_combi, nl_combi, ntmod, ctmod, 
                              aas, aa_masses, type_ms2ions, digits) {
  aas2 <- aa_masses[aas]
  
  if (!is.null(vmods_combi)) {
    aas[as.numeric(names(vmods_combi))] <- vmods_combi
  }
  
  len <- ncol(nl_combi)
  out <- vector("list", len)
  
  nms <- rownames(nl_combi)
  idxes <- which(aas %in% nms)
  
  for (i in 1:len) {
    aas2_af <- aas2
    aas2_af[idxes] <- nl_combi[, i]
    names(aas2_af)[idxes] <- nms
    
    out[[i]] <- calcms2_by_iontype(aas2_af, aa_masses, ntmod, ctmod, 
                                   type_ms2ions, digits)
  }
  
  out
} 


#' Helper: masses of a series of MS2 ions.
#'
#' Updates `aas` and calculates masses for moduels (11) "amods+ tmod- vnl- fnl+"
#' and (12) "amods+ tmod+ vnl- fnl+".
#'
#' It first updates the names of residues in \code{aas} (since `amods+`), then
#' runs column-wisely through the table of \code{nl_combi} (fixed) for
#' corresponding masses of residues.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
mcalcms2_a1_fnl1 <- function (vmods_combi, nl_combi, ntmod = NULL, ctmod = NULL, 
                              aas, aa_masses, type_ms2ions, digits) {
  
  if (!is.null(vmods_combi)) {
    aas[as.numeric(names(vmods_combi))] <- unlist(vmods_combi)
  }

  ms2ions_cols_fn11(aas, nl_combi, aa_masses, ntmod, ctmod, 
                    type_ms2ions, digits) 
}


#' Helper in calculating MS2 masses with "fnl+".
#'
#' For \code{fnl+} only (modules 5, 6, 11, 12).
#'
#' The calculation runs column-wisely through the masses in \code{nl_combi}
#' (fixed). Note that the names of residues (\code{aas2}) have already been
#' updated before the utility was called.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @return A numeric vector
ms2ions_cols_fn11 <- function (aas, nl_combi, aa_masses, ntmod, ctmod, 
                               type_ms2ions, digits) {
  len <- ncol(nl_combi)
  out <- vector("list", len)
  nms <- rownames(nl_combi)
  
  for (i in 1:len) {
    aa_masses[nms] <- nl_combi[, i]
    aas2 <- aa_masses[aas]
    out[[i]] <- calcms2_by_iontype(aas2, aa_masses, ntmod, ctmod, 
                                   type_ms2ions, digits)
  }
  
  invisible(out)
}


#' Helper in calculating MS2 masses.
#' 
#' (1) "amods- tmod- vnl- fnl-".
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_0 <- function (aas, ms1_mass = NULL, aa_masses, type_ms2ions, digits) {
  
  out <- calcms2_by_iontype(aa_masses[aas], aa_masses, 
                            ntmod = NULL, ctmod = NULL, type_ms2ions, digits)
  
  nm <- rep(0, length(aas)) %>% paste0(collapse = "")
  out <- list(out) %>% `names<-`(nm)
}


#' Helper in calculating MS2 masses.
#' 
#' (2) "amods- tmod+ vnl- fnl-".
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a0_t1_nl0 <- function (aas, ms1_mass = NULL, aa_masses, ntmod, ctmod, 
                               type_ms2ions, digits) {
  
  out <- calcms2_by_iontype(aa_masses[aas], aa_masses, 
                            ntmod, ctmod, type_ms2ions, digits)
  
  nm <- rep(0, length(aas)) %>% paste0(collapse = "")
  out <- list(out) %>% `names<-`(nm)
}


#' Helper in calculating MS2 masses.
#' 
#' (3) "amods- tmod+ vnl+ fnl-".
calcms2_a0_t1_vnl1 <- function () {
  message("Combination not possible: `amods- tmod+ vnl+`.")
}


#' Helper in calculating MS2 masses.
#' 
#' (4) "amods- tmod- vnl+ fnl-".
calcms2_a0_t0_vnl1 <- function () {
  message("Combination not possible: `amods- tmod- vnl+`.")
}


#' Helper in calculating MS2 masses.
#'
#' (5) "amods- tmod+ vnl- fnl+".
#'
#' Note that no need to update the \code{aas} (since `amods-`). The utility runs
#' column-wisely through the table of \code{nl_combi} (fixed) for corresponding
#' masses of residues.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a0_t1_fnl1 <- function (nl_combi, aas, aa_masses, ntmod, ctmod, 
                                type_ms2ions, digits) {

  out <- ms2ions_cols_fn11(aas, nl_combi, aa_masses, 
                           ntmod, ctmod, type_ms2ions, digits) 

  nm <- rep(0, length(aas)) %>% paste0(collapse = "")
  names(out) <- paste0(nm, " [", seq_len(ncol(nl_combi)), "]")

  invisible(out)
}


#' Helper in calculating MS2 masses.
#' 
#' (6) "amods- tmod- vnl- fnl+".
#'
#' Note that no need to update the \code{aas} (since `amods-`). The utility runs
#' column-wisely through the table of \code{nl_combi} (fixed) for corresponding
#' masses of residues.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a0_t0_fnl1 <- function (nl_combi, aas, aa_masses, ntmod, ctmod, 
                                type_ms2ions, digits) {
  
  out <- ms2ions_cols_fn11(aas, nl_combi, aa_masses, 
                           ntmod = NULL, ctmod = NULL, type_ms2ions, digits)
  
  nm <- rep(0, length(aas)) %>% paste0(collapse = "")
  names(out) <- paste0(nm, " [", seq_len(ncol(nl_combi)), "]")
  
  invisible(out)
}


#' Helper in calculating MS2 masses.
#' 
#' (7) "amods+ tmod- vnl- fnl-".
#'
#' Updates \code{aas} (since `amods+`) and calculates masses.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a1_t0_nl0 <- function (vmods_combi, nl_combi = NULL, 
                               ntmod = NULL, ctmod = NULL, aas, 
                               aa_masses, type_ms2ions, mod_indexes = NULL, 
                               digits) {

  out <- map(vmods_combi, mcalcms2_a1_vnl0, 
             ntmod = NULL, ctmod = NULL, aas, aa_masses, 
             type_ms2ions, digits = 5)

  out <- add_hexcodes(out, vmods_combi, length(aas), mod_indexes)
}


#' Helper in calculating MS2 ion masses.
#' 
#' (8) "amods+ tmod+ vnl- fnl-".
#' 
#' Updates \code{aas} (since `amods+`) and calculates masses.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a1_t1_nl0 <- function (vmods_combi, nl_combi, ntmod, ctmod, 
                               aas, aa_masses, type_ms2ions, mod_indexes = NULL, 
                               digits) {
  
  # Conflicts btw. amods+ and tmod+
  # NLTLALEALVQLR: the only N on the N-term and cannot be `Deamidated (N)`
  
  # still an empty list() at the end but force an early return here
  if (is_empty(vmods_combi)) return(NULL)
  
  out <- map(vmods_combi, mcalcms2_a1_vnl0, 
             ntmod, ctmod, aas, aa_masses, 
             type_ms2ions, digits = 5)

  out <- add_hexcodes(out, vmods_combi, length(aas), mod_indexes)
}


#' Helper in calculating MS2 ion masses.
#'
#' (9) "amods+ tmod- vnl+ fnl-".
#'
#' Updates \code{aas} (since `amods+`) and runs column-wisely through the table
#' of \code{nl_combi} (variable) for corresponding masses of residues.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a1_t0_nl1 <- function (vmods_combi, nl_combi, ntmod = NULL, 
                               ctmod = NULL, aas, aa_masses, 
                               type_ms2ions, mod_indexes = NULL, digits) {

  out <- map2(vmods_combi, nl_combi, mcalcms2_a1_vnl1, 
              ntmod = NULL, ctmod = NULL, aas, aa_masses, 
              type_ms2ions, digits)

  out <- map2(out, vmods_combi, add_hexcodes_vnl, length(aas), mod_indexes)

  out <- out %>% flatten()
}


#' Helper in calculating MS2 ion masses.
#' 
#' (10) "amods+ tmod+ vnl+ fnl-".
#'
#' Updates \code{aas} (since `amods+`) and runs column-wisely through the table
#' of \code{nl_combi} (variable) for corresponding masses of residues.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a1_t1_nl1 <- function (vmods_combi, nl_combi, ntmod, 
                               ctmod, aas, aa_masses, 
                               type_ms2ions, mod_indexes = NULL, digits) {
  
  # Conflicts between amods+ and tmod+
  # NLTLALEALVQLR: the only N on the N-term and cannot be `Deamidated (N)`

  # still an empty list() at the end if without the early return,
  # but already handled at `calc_ms2ions`
  # if (is_empty(vmods_combi)) return(NULL)
  
  out <- map2(vmods_combi, nl_combi, mcalcms2_a1_vnl1, 
              ntmod, ctmod, aas, aa_masses, 
              type_ms2ions, digits)

  out <- map2(out, vmods_combi, add_hexcodes_vnl, length(aas), mod_indexes)
  
  out <- out %>% flatten()
}


#' Helper in calculating MS2 ion masses.
#' 
#' (11) "amods+ tmod- vnl- fnl+".
#'
#' Updates \code{aas} (since `amods+`) and runs column-wisely through the table
#' of \code{nl_combi} (fixed) for corresponding masses of residues.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a1_t0_fnl1 <- function (vmods_combi, nl_combi, ntmod = NULL, 
                                ctmod = NULL, aas, aa_masses, 
                                type_ms2ions, mod_indexes = NULL, digits) {
  
  out <- map(vmods_combi, mcalcms2_a1_fnl1, 
             nl_combi, ntmod = NULL, ctmod = NULL, 
             aas, aa_masses, type_ms2ions, digits)

  out <- map2(out, vmods_combi, add_hexcodes_fnl, length(aas), mod_indexes)
  
  out <- flatten(out)
}


#' Helper in calculating MS2 ion masses.
#' 
#' (12) "amods+ tmod+ vnl- fnl+".
#'
#' Updates \code{aas} (since `amods+`) and runs column-wisely through the table
#' of \code{nl_combi} (fixed) for corresponding masses of residues.
#'
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calcms2_a1_t1_fnl1 <- function (vmods_combi, nl_combi, ntmod, 
                                ctmod, aas, aa_masses, 
                                type_ms2ions, mod_indexes = NULL, digits) {
  
  # Conflicts between amods+ and tmod+
  # NLTLALEALVQLR: the only N on the N-term and cannot be `Deamidated (N)`
  
  # still an empty list() at the end but force an early return here
  if (is_empty(vmods_combi)) return(NULL)
  
  out <- map(vmods_combi, mcalcms2_a1_fnl1, 
             nl_combi, ntmod, ctmod, 
             aas, aa_masses, type_ms2ions, digits)

  out <- map2(out, vmods_combi, add_hexcodes_fnl, length(aas), mod_indexes)
  
  out <- flatten(out)
}


#' Calculates the masses of MS2 ion series.
#'
#' For a given type of fragmentation. Minimal error handling for speeds.
#' 
#' @param ms1_mass The mass of a theoretical MS1 (for subsetting).
#' @param maxn_vmods_sitescombi_per_pep Integer; the maximum number of
#'   combinatorial variable modifications per peptide sequence.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz".
#' @inheritParams calc_monopep
#' @inheritParams calc_aamasses
#' @import purrr
#' 
#' @examples
#' \donttest{
#' ## No variable modifications
#' # (1)
#' library(magrittr)
#' 
#' fixedmods = NULL
#' varmods = NULL
#' 
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
#'   as.hexmode() %>% 
#'   `names<-`(c(fixedmods, varmods))
#' aa_masses_all <- calc_aamasses(fixedmods, varmods)
#' 
#' x <- calc_ms2ions("MAKEMASSPECFUN", NULL, aa_masses_all[[1]], mod_indexes)
#' 
#' }
calc_ms2ions <- function (aa_seq, ms1_mass = NULL, aa_masses, mod_indexes = NULL, 
                          type_ms2ions = "by", maxn_vmods_per_pep = 5, 
                          maxn_sites_per_vmod = 3, 
                          maxn_vmods_sitescombi_per_pep = 32, digits = 5) {
  
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # carbon <- 12.0
  # oxygen <- 15.99491462
  # nitrogen <- 14.003074
  # h3o_p <- 19.0178415
  # electron <- 0.000549
  
  # moverz <- (mass + h2o + z*proton)/z
  # moverz*z 
  
  if (is.na(aa_seq) || is.null(aa_seq)) return(NULL)
  
  aas <- aa_seq %>% str_split("", simplify = TRUE)
  type <- attr(aa_masses, "type", exact = TRUE)
  
  if (type == "amods- tmod- vnl- fnl-") {
    return(calcms2_0(aas, ms1_mass = NULL, aa_masses, type_ms2ions, digits))
  }
  
  fntmod <- attr(aa_masses, "fntmod", exact = TRUE)
  fctmod <- attr(aa_masses, "fctmod", exact = TRUE)
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  
  if (type == "amods- tmod+ vnl- fnl-") {
    return(calcms2_a0_t1_nl0(aas, ms1_mass = NULL, aa_masses, ntmod, ctmod, 
                             type_ms2ions, digits))
  } 
  
  # --- combinatorial ---
  fmods_ps <- attr(aa_masses, "fmods_ps", exact = TRUE)
  vmods_ps <- attr(aa_masses, "vmods_ps", exact = TRUE)
  fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
  vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
  famods <- attr(aa_masses, "famods", exact = TRUE)
  ftmod <- attr(aa_masses, "ftmod", exact = TRUE)
  amods <- attr(aa_masses, "amods", exact = TRUE)
  tmod <- attr(aa_masses, "tmod", exact = TRUE)
  
  # need ftmod in aa_masses
  # fixed TMT6plex (N-term) can still conflict with Oxidation (M) on N-term
  
  # at least one of `vmods_combi`, `fnl_combi`, `vnl_combi`
  vmods_combi <- combi_mvmods(amods = amods, ntmod = ntmod, ctmod = ctmod, 
                              famods = famods, fntmod = fntmod, fctmod = fctmod, 
                              aas = aas, 
                              aa_masses = aa_masses, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              maxn_vmods_sitescombi_per_pep = 
                                maxn_vmods_sitescombi_per_pep, 
                              digits = digits) %>% 
    find_intercombi_p(maxn_vmods_sitescombi_per_pep) 
  
  # excludes the combination of variable mods with incorrect MS1 masses
  if (!(is_empty(vmods_combi) || is.null(ms1_mass))) {
    idxes <- map_lgl(vmods_combi, check_ms1_mass_vmods, aas, aa_masses, 
                     ntmod, ctmod, ms1_mass)
    
    vmods_combi <- vmods_combi[idxes]
    rm(idxes)
  }

  # Conflicts between `amod+` and `tmods+` (only for "amods+ tmod+ vnl+")
  # "MAEEQGR" - vmods_combi can be NULL: 
  #   cannot be simultaneous `Oxidation (M)` and `Acetyl (Protein N-term)`
  #   (the following only handles early for (10) "amods+ tmod+ vnl+"; 
  #    additional handling of "(8) amods+ tmod+ vnl-" and 
  #    (12) "amods+ tmod+ vnl-" by corresponding functions)
  
  if (is_empty(vmods_combi) && !is_empty(vmods_nl)) return(NULL)
  
  if (is_empty(vmods_nl)) {
    vnl_combi <- NULL
  } else {
    vnl_combi <- map(vmods_combi, ~ expand.grid(vmods_nl[.x]) %>% t())
  }
  
  # not yet subset vnl_combi by ms1_mass...
  # no need if `include_insource_nl = FALSE`

  if (is_empty(fmods_nl)) {
    fnl_combi <- NULL
  } else {
    # no full combinatorial NL for fixed mods
    # (e.g., `Oxidation (M)` either 147 or 83 for every M)
    fnl_combi <- expand.grid(fmods_nl) %>% t()
  }
  
  # not yet subset fnl_combi by ms1_mass...
  # no need if `include_insource_nl = FALSE`

  out <- switch(type, 
                # "amods- tmod- vnl- fnl-" = calcms2_0(aas, ms1_mass = NULL, aa_masses, type_ms2ions, digits), 
                # "amods- tmod+ vnl- fnl-" = calcms2_a0_t1_nl0(aas, ms1_mass = NULL, aa_masses, ntmod, ctmod, type_ms2ions, digits), 
                # "amods- tmod+ vnl+ fnl-" = calcms2_a0_t1_vnl1(), 
                # "amods- tmod- vnl+ fnl-" = calcms2_a0_t0_vnl1(),
                
                "amods- tmod+ vnl- fnl+" = 
                  calcms2_a0_t1_fnl1(fnl_combi, aas, aa_masses, ntmod, ctmod, type_ms2ions, digits), 
                "amods- tmod- vnl- fnl+" = 
                  calcms2_a0_t0_fnl1(fnl_combi, aas, aa_masses, NULL, NULL, type_ms2ions, digits), 
                
                "amods+ tmod- vnl- fnl-" = 
                  calcms2_a1_t0_nl0(vmods_combi, NULL, NULL, NULL, aas, aa_masses, type_ms2ions, mod_indexes, digits), 
                "amods+ tmod+ vnl- fnl-" = 
                  calcms2_a1_t1_nl0(vmods_combi, NULL, ntmod, ctmod, aas, aa_masses, type_ms2ions, mod_indexes, digits), 
                
                "amods+ tmod- vnl+ fnl-" = 
                  calcms2_a1_t0_nl1(vmods_combi, vnl_combi, NULL, NULL, aas, aa_masses, type_ms2ions, mod_indexes, digits), 
                "amods+ tmod+ vnl+ fnl-" = 
                  calcms2_a1_t1_nl1(vmods_combi, vnl_combi, ntmod, ctmod, aas, aa_masses, type_ms2ions, mod_indexes, digits), 
                
                "amods+ tmod- vnl- fnl+" = 
                  calcms2_a1_t0_fnl1(vmods_combi, fnl_combi, NULL, NULL, aas, aa_masses, type_ms2ions, mod_indexes, digits), 
                "amods+ tmod+ vnl- fnl+" = 
                  calcms2_a1_t1_fnl1(vmods_combi, fnl_combi, ntmod, ctmod, aas, aa_masses, type_ms2ions, mod_indexes, digits), 
                
                # "amods- tmod- vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
                # "amods+ tmod- vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
                # "amods- tmod+ vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
                # "amods+ tmod+ vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
                # message("Ignores `vnl+ fnl+`.")
                NULL)
}


#' Adds hex codes (without NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @param ms2ions A series of MS2 ions with masses.
#' @param len The number of amino acid residues for the sequence indicated in
#'   \code{ms2ions}.
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_aamasses
add_hexcodes <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) {
  hex_mods = rep("0", len)
  rows <- map_lgl(ms2ions, ~ !is.null(.x))
  
  vmods_combi <- vmods_combi[rows]
  ms2ions <- ms2ions[rows]
  
  hex_mods <- map(vmods_combi, ~ {
    hex_mods[as.numeric(names(.x))] <- mod_indexes[unlist(.x)]
    hex_mods <- paste0(hex_mods, collapse = "")
    hex_mods
  }, hex_mods)
  
  names(ms2ions) <- hex_mods
  
  invisible(ms2ions)
}


#' Adds hex codes (with fixed NLs).
#'
#' To indicate the variable modifications of an amino acid sequence. The only
#' difference to \link{add_hexcodes_vnl} is square brackets (for fixed NLs) or
#' round parenthesis (for variable NLs).
#'
#' @inheritParams add_hexcodes
add_hexcodes_fnl <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) {
  hex_mods = rep("0", len)
  hex_mods[as.numeric(names(vmods_combi))] <- mod_indexes[unlist(vmods_combi)]
  hex_mods <- paste0(hex_mods, collapse = "")
  
  names(ms2ions) <- paste0(hex_mods, " [", seq_along(ms2ions), "]")
  
  invisible(ms2ions)
}


#' Adds hex codes (with variable NLs).
#' 
#' To indicate the variable modifications of an amino acid sequence.
#' 
#' @inheritParams add_hexcodes
add_hexcodes_vnl <- function (ms2ions, vmods_combi, len, mod_indexes = NULL) {
  hex_mods = rep("0", len)
  hex_mods[as.numeric(names(vmods_combi))] <- mod_indexes[unlist(vmods_combi)]
  hex_mods <- paste0(hex_mods, collapse = "")
  
  names(ms2ions) <- paste0(hex_mods, " (", seq_along(ms2ions), ")")
  
  invisible(ms2ions)
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
combi_mvmods <- function (amods, ntmod, ctmod, 
                          famods, fntmod, fctmod, 
                          aas, 
                          aa_masses, 
                          maxn_vmods_per_pep = 5, 
                          maxn_sites_per_vmod = 3, 
                          maxn_vmods_sitescombi_per_pep = 32, 
                          digits = 5) {
  # (6) "amods- tmod- vnl- fnl+"
  if (is_empty(amods)) return(NULL)
  
  residue_mods <- unlist(amods, use.names = FALSE) %>% 
    `names<-`(names(amods)) %>% 
    split(., .)
  
  # Oxidation (M) Carbamidomethyl (M) 
  # "M"                 "M" 
  # 
  # $S
  # dHex (S) 
  # "S" 
  
  map(residue_mods, ~ combi_vmods(
    aas, .x, 
    ntmod, ctmod, 
    fntmod, fctmod, 
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
combi_vmods <- function (aas, 
                         residue_mods, 
                         ntmod, ctmod, 
                         fntmod, fctmod, 
                         aa_masses, 
                         maxn_vmods_per_pep = 5, 
                         maxn_sites_per_vmod = 3, 
                         maxn_vmods_sitescombi_per_pep = 32, 
                         digits = 5) {
  
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

  p <- check_tmod_p(aas, residue, p, ntmod, ctmod)
  p <- check_tmod_p(aas, residue, p, fntmod, fctmod)
  
  if (length(p) == 1) {
    names(n) <- p
    return(n)
  }

  # --- combinations ---
  if (length(n) == 1) { # "Oxidation (M)"
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
    ns <- map(seq_along(p), ~ {
      expand.grid(rep(list(n), length(p[1:.x])), KEEP.OUT.ATTRS = FALSE, 
                  stringsAsFactors = FALSE) %>% t()
    }, n, p) 
    
    ps <- map(seq_along(p), ~ {
      combn(as.character(p), .x)
    }, p)
    
    out <- mvmcombC(ns, ps)
  }
  
  # ---
  out <- out %>% 
    flatten() 
  
  rows <- map_lgl(out, ~ length(.x) > maxn_vmods_per_pep)
  out <- out[!rows]
  
  maxn_vmod <- out %>% 
    map(unlist) %>% 
    map(table) %>% 
    map(max)
  rows <- maxn_vmod > maxn_sites_per_vmod
  out <- out[!rows]
  
  if (!is_empty(out)) {
    out <- out %>% 
      `[`(1:min(length(.), maxn_vmods_sitescombi_per_pep))
  }
  
  invisible(out)
}


#' Checks the compatibility between terminal and anywhere modifications.
#' 
#' For both fixed and variable terminal modifications
#' 
#' @param residue The residue (e.g. "M").
#' @param p The position indexes of \code{residue}.
#' @param ntmod The parsed fixed or variable \code{N-term} modifications from
#'   \code{aa_masses}.
#' @param ctmod The parsed fixed or variable \code{C-term} modifications from
#'   \code{aa_masses}.
#' @inheritParams combi_vmods
#'   
#' @return Updated position \code{p}.
check_tmod_p <- function (aas, residue, p, ntmod, ctmod) {
  if (!(is_empty(ntmod) || is_empty(ctmod))) {
    len_aas <- length(aas)
    
    if (aas[1] == residue && aas[len_aas] == residue) {
      p <- p[-length(p)]
      p <- p[-1]
    } else if (aas[1] == residue) {
      p <- p[-1]
    } else if (aas[len_aas] == residue) {
      p <- p[-length(p)]
    }
  } else if (!is_empty(ntmod)) {
    if (aas[1] == residue) {
      p <- p[-1]
    }
  } else if (!is_empty(ctmod)) {
    if (aas[length(aas)] == residue) {
      p <- p[-length(p)]
    }
  }
  
  invisible(p)
}


#' Finds the combinations of positions and sites across residues.
#' 
#' For multiple residues (each residue one to multiple modifications).
#' 
#' @param intra_combis The results from \link{unique_mvmods}.
#' @inheritParams calc_ms2ions
#' @seealso \link{find_intercombi}
find_intercombi_p <- function (intra_combis, maxn_vmods_sitescombi_per_pep) {
  # i.e. conflicting M, N-term where the only M vmod on N-term

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
    
    # p_combis <- p_combis %>% reduce(expand.grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    # v_combis <- v_combis %>% reduce(expand.grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    
    p_combis <- expand.grid(p_combis, KEEP.OUT.ATTRS = FALSE, 
                            stringsAsFactors = FALSE)
    v_combis <- expand.grid(v_combis, KEEP.OUT.ATTRS = FALSE, 
                            stringsAsFactors = FALSE)
    
    nrow <- min(nrow(v_combis), maxn_vmods_sitescombi_per_pep)
    
    p_out <- v_out <- vector("list", nrow)
    
    for (i in seq_len(nrow)) {
      # p_out[[i]] <- reduce(p_combis[i, ], `c`) %>% unlist()
      # v_out[[i]] <- reduce(v_combis[i, ], `c`) %>% unlist()
      
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



#' Masses of singly-charged b- and y-ions.
#'
#' @inheritParams calc_bions
#' @inheritParams calc_yions
#' @seealso \link{add_complement_ions}
calc_byions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  bs <- calc_bions(ntmod, aas, aa_masses, digits)
  ys <- calc_yions(ctmod, aas, aa_masses, digits)
  c(bs, ys)
}


#' Masses of singly-charged b0- and y0-ions.
#' 
#' @inheritParams calc_bions
#' @seealso \link{add_complement_ions}
calc_by0ions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  bs <- calc_b0ions(ntmod, aas, aa_masses, digits)
  ys <- calc_y0ions(ctmod, aas, aa_masses, digits)
  c(bs, ys)
}


#' Masses of singly-charged b*- and y*-ions.
#' 
#' @inheritParams calc_bions
#' @seealso \link{add_complement_ions}
calc_bystarions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  bs <- calc_bstarions(ntmod, aas, aa_masses, digits)
  ys <- calc_ystarions(ctmod, aas, aa_masses, digits)
  c(bs, ys)
}



#' Masses of doubly-charged b- and y-ions.
#' 
#' @inheritParams calc_bions
#' @seealso \link{add_complement_ions}
calc_by2ions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  bs <- calc_b2ions(ntmod, aas, aa_masses, digits)
  ys <- calc_y2ions(ctmod, aas, aa_masses, digits)
  c(bs, ys)
}


#' Masses of doubly-charged b0- and y0-ions.
#' 
#' @inheritParams calc_bions
#' @seealso \link{add_complement_ions}
calc_by02ions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  bs <- calc_b02ions(ntmod, aas, aa_masses, digits)
  ys <- calc_y02ions(ctmod, aas, aa_masses, digits)
  c(bs, ys)
}


#' Masses of doubly-charged b*- and y*-ions.
#' 
#' @inheritParams calc_bions
#' @seealso \link{add_complement_ions}
calc_bystar2ions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  bs <- calc_bstar2ions(ntmod, aas, aa_masses, digits)
  ys <- calc_ystar2ions(ctmod, aas, aa_masses, digits)
  c(bs, ys)
}


#' Masses of singly-charged b-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_bions <- function (ntmod, aas, aa_masses, digits = 5) {
  # aa_masses["N-term"] is neutral, e.g. H (1.007825 ), not H+ (1.00727647)
  electron <- 0.000549

  if (is_empty(ntmod)) {
    ions <- c(aa_masses["N-term"] - electron, aas)
  } else {
    ions <- c(aa_masses[names(ntmod)] - electron, aas) 
  }
  
  ions %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}



#' Masses of doubly-charged b-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_b2ions <- function (ntmod, aas, aa_masses, digits = 5) {
  ions <- calc_bions(ntmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Masses of singly-charged b*-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_bstarions <- function (ntmod, aas, aa_masses, digits = 5) {
  
  # electron <- 0.000549
  # -NH3:17.026549
  # -17.027098 = -17.026549 - 0.000549
  
  if (is_empty(ntmod)) {
    ions <- c(aa_masses["N-term"] - 17.027098, aas)
  } else {
    ions <- c(aa_masses[names(ntmod)] - 17.027098, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of doubly-charged b*-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_bstar2ions <- function (ntmod, aas, aa_masses, digits = 5) {
  ions <- calc_bstarions(ntmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Masses of singly-charged b^0-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_b0ions <- function (ntmod, aas, aa_masses, digits = 5) {
  
  # -H2O 18.010565
  # electron <- 0.000549
  # -18.011114 = -18.010565 -0.000549 
  
  if (is_empty(ntmod)) {
    ions <- c(aa_masses["N-term"] - 18.011114, aas)
  } else {
    ions <- c(aa_masses[names(ntmod)] - 18.011114, aas)
  } 
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of doubly-charged b0-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calc_b02ions <- function (ntmod, aas, aa_masses, digits = 5) {
  ions <- calc_b0ions(ntmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}



#' Masses of singly-charged y-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_yions <- function (ctmod, aas, aa_masses, digits = 5) {
  
  # (1) OH (C-term), + H (N-term) + H+
  # (2) Other C-term (other than OH) + H + H+: X + 1.007825 + 1.00727647

  aas <- rev(aas)
  
  if (is_empty(ctmod)) { 
    ions <- c(aa_masses["C-term"] + 2.01510147, aas)
  } else { 
    ions <- c(aa_masses[names(ctmod)] + 2.01510147, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Masses of doubly-charged y-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_y2ions <- function (ctmod, aas, aa_masses, digits = 5) {
  ions <- calc_yions(ctmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Masses of singly-charged y*-ions 
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_ystarions <- function (ctmod, aas, aa_masses, digits = 5) {
  # (1) OH (C-term), + H (N-term) + H+
  # (2) Other C-term + H + H+: X + 1.007825 + 1.00727647
  
  # -NH3 -17.026549
  # (1) 19.0178415 - 17.026549
  # (2) 2.01510147 - 17.026549
  
  aas <- rev(aas)
  
  if (is_empty(ctmod)) { 
    ions <- c(1.9912925, aas)
  } else { 
    ions <- c(aa_masses[names(ctmod)] - 15.0114475, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of doubly-charged y*-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
calc_ystar2ions <- function (ctmod, aas, aa_masses, digits = 5) {
  ions <- calc_ystarions(ctmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Masses of singly-charged y0-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_y0ions <- function (ctmod, aas, aa_masses, digits = 5) {
  # (1) OH (C-term), + H (N-term) + H+
  # (2) Other C-term + H + H+: X + 1.007825 + 1.00727647
  
  # -H2O - 18.010565
  # (1) 19.0178415 - 18.010565
  # (2) 2.01510147 - 18.010565
  
  aas <- rev(aas)
  
  if (is_empty(ctmod)) { 
    ions <- c(1.0072765, aas)
  } else { 
    ions <- c(aa_masses[names(ctmod)] - 15.9954635, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of doubly-charged y0-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_y02ions <- function (ctmod, aas, aa_masses, digits = 5) {
  ions <- calc_y0ions(ctmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}





#' Masses of singly-charged c- and z-ions.
#' 
#' @inheritParams calc_bions
#' @inheritParams calc_yions
#' @seealso \link{add_complement_ions}
calc_czions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  cs <- calc_cions(ntmod, aas, aa_masses, digits)
  zs <- calc_zions(ctmod, aas, aa_masses, digits)
  c(cs, zs)
}


#' Masses of doubly-charged c- and z-ions.
#' 
#' @inheritParams calc_bions
#' @seealso \link{add_complement_ions}
calc_cz2ions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  cs <- calc_c2ions(ntmod, aas, aa_masses, digits)
  zs <- calc_z2ions(ctmod, aas, aa_masses, digits)
  c(cs, zs)
}


#' Masses of singly-charged c-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_cions <- function (ntmod, aas, aa_masses, digits = 5) {
  # aa_masses["N-term"] is H (1.007825 ), not H+ (1.00727647)
  # H+ + NH3: 18.0338255 == 1.00727647 + 17.026549
  # NH3 - e: 17.026 == 17.026549 - 0.000549

  if (is_empty(ntmod)) {
    ions <- c(aa_masses["N-term"] + 17.026, aas)
  } else {
    ions <- c(aa_masses[names(ntmod)] + 17.026, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Masses of doubly-charged c-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_c2ions <- function (ntmod, aas, aa_masses, digits = 5) {
  ions <- calc_cions(ntmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Masses of singly-charged z-ions 
#' 
#' Singly charged.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_zions <- function (ctmod, aas, aa_masses, digits = 5) {
  # (1) [H3O+ - NH3], 1.9912925 == 19.0178415 - 17.026549
  # (2) Other C-term (other than OH) + H + H+ - NH3
  # 15.0114475 == 1.007825 + 1.00727647 - 17.026549
  
  aas <- rev(aas)
  
  if (is_empty(ctmod)) {
    ions <- c(aa_masses["C-term"] - 15.0114475, aas)
  } else {
    ions <- c(aa_masses[names(ctmod)] - 15.0114475, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Masses of doubly-charged z-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_z2ions <- function (ctmod, aas, aa_masses, digits = 5) {
  ions <- calc_zions(ctmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}





#' Masses of singly-charged a- and x-ions.
#' 
#' @inheritParams calc_bions
#' @inheritParams calc_yions
#' @seealso \link{add_complement_ions}
calc_axions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  as <- calc_aions(ntmod, aas, aa_masses, digits)
  xs <- calc_xions(ctmod, aas, aa_masses, digits)
  c(as, xs)
}


#' Masses of doubly-charged a- and x-ions.
#' 
#' @inheritParams calc_bions
#' @seealso \link{add_complement_ions}
calc_ax2ions <- function (ntmod, ctmod, aas, aa_masses, digits = 5) {
  as <- calc_a2ions(ntmod, aas, aa_masses, digits)
  xs <- calc_x2ions(ctmod, aas, aa_masses, digits)
  c(as, xs)
}


#' Masses of singly-charged a-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_aions <- function (ntmod, aas, aa_masses, digits = 5) {
  # aa_masses["N-term"] is H (1.007825), not H+ (1.00727647)
  # H+ - CO:  26.9876381 == 1.00727647 - 27.9949146
  # CO - e: 27.9943656 == 27.9949146 - 0.000549
  # e + CO = 27.9954636

  if (is_empty(ntmod)) {
    ions <- c(aa_masses["N-term"] - 27.9954636, aas)
  } else {
    ions <- c(aa_masses[names(ntmod)] - 27.9943656, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Masses of doubly-charged a-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_a2ions <- function (ntmod, aas, aa_masses, digits = 5) {
  ions <- calc_aions(ntmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Masses of singly-charged a*-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_astarions <- function (ntmod, aas, aa_masses, digits = 5) {
  # H+ - CO - NH3:  -44.0141871 == 1.00727647 - 27.9949146 - 17.026549
  # CO -NH3 - e: -45.0220126 == -27.9949146 - 17.026549 - 0.000549

  if (is_empty(ntmod)) {
    ions <- c(-44.0141871, aas)
  } else {
    ions <- c(aa_masses[names(ntmod)] - 45.0220126, aas)
  } 
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of doubly-charged a*-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_astar2ions <- function (ntmod, aas, aa_masses, digits = 5) {
  ions <- calc_astarions(ntmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Masses of singly-charged a0-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_a0ions <- function (ntmod, aas, aa_masses, digits = 5) {
  # H+ - CO - H2O:  -44.9982031 == 1.00727647 - 27.9949146 - 18.010565
  # CO -H2O - e: -46.0060286 == -27.9949146 - 18.010565 - 0.000549

  if (is_empty(ntmod)) {
    ions <- c(-44.9982031, aas)
  } else {
    ions <- c(aa_masses[names(ntmod)] - 46.0060286, aas)
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of doubly-charged a0-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_a02ions <- function (ntmod, aas, aa_masses, digits = 5) {
  ions <- calc_a0ions(ntmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}




#' Masses of singly-charged x-ions.
#' 
#' Singly charged.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @import purrr
#' @seealso \link{add_complement_ions}
calc_xions <- function (ctmod, aas, aa_masses, digits = 5) {
  # (1) [H3O+ + CO - H2], 44.9971061 == 19.0178415 + 27.9949146 - 2 * 1.007825
  # (2) Other C-term (other than OH) + H + H+ + CO - 2H
  # 27.9943661 == 1.007825 + 1.00727647 + 27.9949146 - 2 * 1.007825
  #            == CO - e
  
  aas <- rev(aas)
  
  if (is_empty(ctmod)) {
    ions <- c(aa_masses["C-term"] + 27.9943661, aas)
  } else {
    ions <- c(aa_masses[names(ctmod)] + 27.9943661, aas)
  } 
  
  ions %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Masses of doubly-charged x-ions.
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
#' @seealso \link{add_complement_ions}
calc_x2ions <- function (ctmod, aas, aa_masses, digits = 5) {
  ions <- calc_xions(ctmod, aas, aa_masses, digits)
  (ions + 1.00727647)/2
}


#' Adds the complimentary series of ions.
#' 
#' @param ms1 The MS1 mass of precursor (neutral from \link{calc_pepmasses}).
#' @param ms2 The series of MS2 product ions in one of a-, b- or c-ions.
#' @importFrom purrr map map2 map_dbl
add_complement_ions <- function (ms1, ms2) {
  # + 2 protons
  # (the first to offset the H+ in `ms2`)
  # (the second for carrying charge)
  
  ms2c <- map_dbl(ms2, ~ ms1 + 2.01455294 - .x)
  
  len <- length(ms2c)
  names(ms2c)[1:(len-1)] <- names(ms2c)[2:len]
  
  c(ms2[-1], ms2c[-len])
}




#' Calculates the mono-isotopic mass of a MS2 ions.
#'
#' For direct uses from an R console (with trade-offs in speed).
#'
#' @inheritParams calc_ms2ions
#' @inheritParams calc_monopeptide
#' @examples
#' \donttest{
#' ## No variable modifications
#' # (1)
#' x <- calc_ms2ionseries("MAKEMASSPECFUN", 
#'                        fixedmods = NULL, 
#'                        varmods = NULL)
#' 
#' x$mass
#' 
#' # (2) no combinatorial NL for fixed modifications
#' x <- calc_ms2ionseries("MAKEMASSPECFUN", 
#'                        fixedmods = "Oxidation (M)", 
#'                        varmods = NULL)
#'                        
#' x$mass
#' 
#' ## With variable modifications
#' # (3) combinatorial sites and NL available
#' x <- calc_ms2ionseries("MAKEMASSPECFUN", 
#'                        fixedmods = NULL, 
#'                        varmods = "Oxidation (M)")
#' 
#' x$mass
#' # x$vmods_ps
#' 
#' # (4)
#' x <- calc_ms2ionseries("MAKEMASSPECFUN",
#'                        c("TMT6plex (N-term)", 
#'                          "TMT6plex (K)", 
#'                          "Carbamidomethyl (C)"),
#'                        c("Acetyl (N-term)", 
#'                          "Gln->pyro-Glu (N-term = Q)", 
#'                          "Oxidation (M)"))
#'                       
#' x$mass
#' 
#' # The N-term M realizes with acetylation
#' x$vmods_ps[[5]]
#' 
#' # (5) Neutral losses for occurrences of both fixed 
#' #     and variable modifications ignored
#' x <- calc_ms2ionseries("MAKEMASSPECFUN",
#'                        c("TMT6plex (N-term)", 
#'                          "Oxidation (M)", 
#'                          "Deamidated (N)"), 
#'                        c("dHex (S)"))
#'                       
#' x$mass[[2]]
#' 
#' # Change from fixed to variable for full combinatorials
#' x <- calc_ms2ionseries("MAKEMASSPECFUN",
#'                        c("TMT6plex (N-term)", 
#'                          "Deamidated (N)"), 
#'                        c("Oxidation (M)", 
#'                          "dHex (S)"))
#'                       
#' x$mass[[4]]
#' x$vmods_ps[[4]]
#' }
#' @export
calc_ms2ionseries <- function (aa_seq, fixedmods, varmods, 
                               type_ms2ions = "by", ms1_mass = NULL, 
                               maxn_vmods_setscombi = 64,
                               maxn_vmods_per_pep = Inf, 
                               maxn_sites_per_vmod = Inf, 
                               maxn_vmods_sitescombi_per_pep = 32, 
                               digits = 5) {
  
  options(digits = 9)
  
  aa_masses_all <- calc_aamasses(fixedmods, varmods, maxn_vmods_setscombi)
  
  peps <- purrr::map(aa_masses_all, subpeps_by_vmods, aa_seq) %>% 
    purrr::flatten()
  
  oks <- purrr::map_lgl(peps, ~ !purrr::is_empty(.x))
  peps <- peps[oks]
  aa_masses_all <- aa_masses_all[oks]
  
  mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
    as.hexmode() %>% 
    `names<-`(c(fixedmods, varmods))
  
  ms <- purrr::map2(peps, aa_masses_all, ~ {
    pri <- calc_ms2ions(.x, ms1_mass, .y, mod_indexes, type_ms2ions, 
                        maxn_vmods_per_pep, maxn_sites_per_vmod, 
                        maxn_vmods_sitescombi_per_pep, digits)
    
    sec <- map(pri, add_seions, type_ms2ions, digits)
    
    list(pri = pri, sec = sec)
  })
  
  attrs <- purrr::map(aa_masses_all, attributes)
  vmods_ps <- map(attrs, `[[`, "vmods_ps")
  
  list(mass = map(ms, `[[`, "pri"), 
       sec_mass = map(ms, `[[`, "sec"), 
       vmods_ps = vmods_ps)
}



