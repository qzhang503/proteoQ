#' Calculates the masses of MS2 ion series.
#'
#' For a given type of fragmentation. Minimal error handling for speeds.
#'
#' @param maxn_vmods_sitescombi_per_pep Integer; the maximum number of
#'   combinatorial variable modifications per peptide sequence.
#' @param type Character; the type of fragment ions.
#' @param z Integer; the charge state of an MS1 ion. Currently only for postive
#'   charges.
#' @inheritParams calc_monopep
#' @import purrr
calc_ms2ionseries <- function (aa_seq = NA, 
                               aa_masses = NA, 
                               type = "b", 
                               maxn_vmods_per_pep = 5, 
                               maxn_sites_per_vmod = 3, 
                               maxn_vmods_sitescombi_per_pep = 32, 
                               z = 1, 
                               digits = 5) {
  
  options(digits = 9)
  
  vmods_ps <- aa_masses %>% 
    attributes() %>% 
    `[[`("vmods_ps")
  
  # multiple mods to [NC]-term already excluded from aa_masses
  amods <- local({
    sites <- vmods_ps %>% 
      map(~ .x[grepl("Anywhere", names(.x))]) 
    
    empties <- sites %>% map_lgl(is_empty)
    
    sites <- sites[!empties] 
    
    sites
  })
  
  tmod <- vmods_ps %>% .[! . %in% amods]
  if (is_empty(tmod)) {
    tmod <- NULL
  } else if (tmod == "") {
    tmod <- NULL
  }
  
  ntmod <- tmod %>% .[. == "N-term"]
  ctmod <- tmod %>% .[. == "C-term"]
  
  calc_ms2ions(aa_seq = aa_seq, 
               type = type, 
               aa_masses = aa_masses, 
               vmods_ps = vmods_ps, 
               amods = amods, 
               tmod = tmod, 
               ntmod = ntmod, 
               ctmod = ctmod, 
               maxn_vmods_per_pep = maxn_vmods_per_pep, 
               maxn_sites_per_vmod = maxn_sites_per_vmod, 
               z = z, 
               digits = digits)
}



#' Helper With position combinations
#'
#' @inheritParams calc_monopep
#' @inheritParams calc_ms2ionseries
#' @import purrr
#' @importFrom stringr str_split
calc_ms2ions <- function (aa_seq, type, aa_masses, vmods_ps, 
                          amods, tmod, ntmod, ctmod, 
                          maxn_vmods_per_pep = 3, 
                          maxn_sites_per_vmod = 5, 
                          maxn_vmods_sitescombi_per_pep = 32, 
                          z = 1, 
                          digits = 5) {
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # carbon <- 12.0
  # oxygen <- 15.99491462
  # nitrogen <- 14.003074
  # h3o <- 19.0178415
  # electron <- 0.000549
  
  # moverz <- (mass + h2o + z*proton)/z
  
  options(digits = 9)
  
  if (is.na(aa_seq)) {
    return(NULL)
  }
  
  aas <- aa_seq %>% 
    str_split("", simplify = TRUE)
  
  if (is_empty(amods) && is_empty(tmod)) {
    return(
      calc_mvmods_ms2masses(vmods = NULL, ntmod = NULL, ctmod = NULL, 
                            aa_seq = aa_seq, aas = aas, aa_masses = aa_masses, 
                            type = type, z = z, digits = digits)
    )
  } else if (is_empty(amods)) { 
    # tmod only
    return(
      calc_mvmods_ms2masses(vmods = NULL, ntmod = ntmod, ctmod = ctmod, 
                            aa_seq = aa_seq, aas = aas, aa_masses = aa_masses, 
                            type = type, z = z, digits = digits)
    )
  }
  
  # amods and/or tmods
  vmods_combi <- local({
    intra_combis <- combi_mvmods(amods = amods, ntmod = ntmod, ctmod = ctmod, 
                                 aas = aas, 
                                 aa_masses = aa_masses, 
                                 maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                 maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                 maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
                                 digits = digits)
    
    # i.e. conflicting M, N-term where an M vmod on N-term
    if (any(map_lgl(intra_combis, is_empty))) { 
      v_out <- list()
    } else if (length(intra_combis) == 1L) { # M, one to multiple positions
      if (length(intra_combis[[1]]) == 1L) { # 2: "Oxidation (M)"
        v_out <- intra_combis
      } else { # 2: "Oxidation (M)"; 3: "Oxidation (M)"; 2: "Oxidation (M)", 3: "Oxidation (M)"
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
      
      p_combis <- p_combis %>% 
        reduce(expand.grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      
      v_combis <- v_combis %>% 
        reduce(expand.grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      
      nrow <- min(nrow(v_combis), maxn_vmods_sitescombi_per_pep)
      
      p_out <- v_out <- vector("list", nrow)
      
      for (i in seq_len(nrow)) {
        p_out[[i]] <- reduce(p_combis[i, ], `c`) %>% unlist()
        v_out[[i]] <- reduce(v_combis[i, ], `c`) %>% unlist()
        names(v_out[[i]]) <- p_out[[i]]
        v_out[[i]] <- v_out[[i]][order(as.numeric(names(v_out[[i]])))]
      }
      
      # v_out <- map2(p_combis[[1]], v_combis[[1]], ~ {
      #   names(.y) <- .x
      #   .y
      # })
    }

    v_out
  })
  
  vmods_combi %>% 
    map(calc_mvmods_ms2masses, ntmod, ctmod, aa_seq, aas, aa_masses, type, z, digits)
}


#' Finds the combinations of variable modifications (multiple sites)
#'
#' For all the \code{Anywhere} modifications specified in \code{amods}.
#'
#' @param amods \code{Anywhere} modifications.
#' @inheritParams calc_monopep
#' @inheritParams mcalc_monopep
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
#' @import purrr
#' @return Lists by residues in \code{amods}.
combi_mvmods <- function (amods, ntmod, ctmod, 
                          aas, 
                          aa_masses, 
                          maxn_vmods_per_pep = 5, 
                          maxn_sites_per_vmod = 3, 
                          maxn_vmods_sitescombi_per_pep = 32, 
                          digits) {
  
  residue_mods <- unlist(amods, use.names = FALSE) %>% 
    `names<-`(names(amods)) %>% 
    split(., .)
  
  map(residue_mods, ~ combi_vmods(
    aas, .x, ntmod, ctmod, 
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
                         ntmod, 
                         ctmod, 
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
  
  # i.e., btw Anywhere "M" and "Acetyl N-term" where "M" on the "N-term"
  # MFGMFNVSMR cannot have three `Oxidation (M)` and `Acetyl (N-term)`
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
  
  # test: combi only with more than one site
  # "QMNLAAENR"
  if (length(p) == 1) {
    names(n) <- p
    return(n)
  }

  # seq_len?
  ns <- map(seq_along(p), ~ {
    expand.grid(rep(list(n), length(p[1:.x])), KEEP.OUT.ATTRS = FALSE, 
                stringsAsFactors = FALSE)
  }, n, p) 
  
  ps <- map(seq_along(p), ~ {
    combn(p, .x)
  }, p)
  
  out <- map2(ns, ps, ~ {
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
    
    np
  }) %>% 
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


#' Masses of a series of MS2 ions.
#' 
#' Singly charged for now. 
#' 
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
calc_mvmods_ms2masses <- function (vmods, ntmod, ctmod, 
                                   aa_seq, aas, aa_masses, type, z = 1, digits) {
  if (!is.null(vmods)) {
    aas[as.numeric(names(vmods))] <- unlist(vmods)
  }

  switch(type, 
         b = calc_bions(ntmod, ctmod, aas, aa_masses, z, digits), 
         y = calc_yions(ntmod, ctmod, aas, aa_masses, z, digits), 
         c = calc_cions(ntmod, ctmod, aas, aa_masses, z, digits), 
         z = calc_zions(ntmod, ctmod, aas, aa_masses, z, digits), 
         a = calc_aions(ntmod, ctmod, aas, aa_masses, z, digits), 
         x = calc_xions(ntmod, ctmod, aas, aa_masses, z, digits), 
         stop("Unknown type.", call. = FALSE))

}


#' Masses of b-ions 
#' 
#' Singly charged.
#' 
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
#' @import purrr
calc_bions <- function (ntmod, ctmod, aas, aa_masses, z, digits) {
  if (is_empty(ntmod) && is_empty(ctmod)) {
    # Anywhere mods, no [NC]-term mods
    ions <- c(aa_masses["N-term"], aa_masses[aas])
  } else if (!(is_empty(ntmod) || is_empty(ctmod))) {
    # Anywhere, both N- and C-term mods
    ions <- c(aa_masses[names(ntmod)], aa_masses[aas])
  } else if (!is_empty(ntmod)) {
    # Anywhere, N-term
    ions <- c(aa_masses[names(ntmod)], aa_masses[aas])
  } else if (!is_empty(ctmod)) {
    # Anywhere C-term
    ions <- c(aa_masses["N-term"], aa_masses[aas])
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of a-ions 
#' 
#' Singly charged.
#' 
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
#' @import purrr
calc_aions <- function (ntmod, ctmod, aas, aa_masses, z, digits) {
  if (is_empty(ntmod) && is_empty(ctmod)) {
    # Anywhere mods, no [NC]-term mods
    ions <- c(aa_masses["N-term"] - 27.9949146, aa_masses[aas])
  } else if (!(is_empty(ntmod) || is_empty(ctmod))) {
    # Anywhere, both N- and C-term mods
    ions <- c(aa_masses[names(ntmod)] - 27.9949146, aa_masses[aas])
  } else if (!is_empty(ntmod)) {
    # Anywhere, N-term
    ions <- c(aa_masses[names(ntmod)] - 27.9949146, aa_masses[aas])
  } else if (!is_empty(ctmod)) {
    # Anywhere C-term
    ions <- c(aa_masses["N-term"] - 27.9949146, aa_masses[aas])
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of c-ions 
#' 
#' Singly charged.
#' 
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
#' @import purrr
calc_cions <- function (ntmod, ctmod, aas, aa_masses, z, digits) {
  # + NH3:17.026549
  
  if (is_empty(ntmod) && is_empty(ctmod)) {
    # Anywhere mods, no [NC]-term mods
    ions <- c(aa_masses["N-term"] + 17.026549, aa_masses[aas])
  } else if (!(is_empty(ntmod) || is_empty(ctmod))) {
    # Anywhere, both N- and C-term mods
    ions <- c(aa_masses[names(ntmod)] + 17.026549, aa_masses[aas])
  } else if (!is_empty(ntmod)) {
    # Anywhere, N-term
    ions <- c(aa_masses[names(ntmod)] + 17.026549, aa_masses[aas])
  } else if (!is_empty(ctmod)) {
    # Anywhere C-term
    ions <- c(aa_masses["N-term"] + 17.026549, aa_masses[aas])
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of y-ions 
#' 
#' Singly charged.
#' 
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
#' @import purrr
calc_yions <- function (ntmod, ctmod, aas, aa_masses, z, digits) {
  # (1) OH (C-term), + H (N-term) + H+
  # (2) Other C-term + H + H+: X + 1.007825 + 1.00727647
  
  if (is_empty(ntmod) && is_empty(ctmod)) { # anywhere mods, no term mods
    ions <- c(19.0178415, aa_masses[rev(aas)])
  } else if (!(is_empty(ntmod) || is_empty(ctmod))) { # anywhere, both N and C-term mods
    ions <- c(aa_masses[names(ctmod)] + 2.01510147, aa_masses[rev(aas)])
  } else if (!is_empty(ntmod)) { # anywhere, N-term
    ions <- c(19.0178415, aa_masses[rev(aas)])
  } else if (!is_empty(ctmod)) { # anywhere C-term
    ions <- c(aa_masses[names(ctmod)] + 2.01510147, aa_masses[rev(aas)])
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of x-ions 
#' 
#' Singly charged.
#' 
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
#' @import purrr
calc_xions <- function (ntmod, ctmod, aas, aa_masses, z, digits) {
  # (1) [H3O+ + CO - H2], 44.9971061: 19.0178415 + 27.9949146 - 2 * 1.007825
  # (2) Other C-term + H + H+ + CO - 2H
  # X + 1.007825 + 1.00727647 + 27.9949146 - 2 * 1.007825

  if (is_empty(ntmod) && is_empty(ctmod)) { # anywhere mods, no term mods
    ions <- c(44.9971061, aa_masses[rev(aas)])
  } else if (!(is_empty(ntmod) || is_empty(ctmod))) { # anywhere, both N and C-term mods
    ions <- c(aa_masses[names(ctmod)] + 27.9943661, aa_masses[rev(aas)])
  } else if (!is_empty(ntmod)) { # anywhere, N-term
    ions <- c(44.9971061, aa_masses[rev(aas)])
  } else if (!is_empty(ctmod)) { # anywhere, C-term
    ions <- c(aa_masses[names(ctmod)] + 27.9943661, aa_masses[rev(aas)])
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}


#' Masses of z-ions 
#' 
#' Singly charged.
#' 
#' @inheritParams calc_mvmods_masses
#' @inheritParams calc_ms2ionseries
#' @import purrr
calc_zions <- function (ntmod, ctmod, aas, aa_masses, z, digits) {
  # (1) [H3O+ - NH3], 1.9912925: 19.0178415 - 17.026549
  # (2) Other C-term + H + H+ - NH3
  # X + 1.007825 + 1.00727647 - 17.026549

  if (is_empty(ntmod) && is_empty(ctmod)) { # anywhere mods, no term mods
    ions <- c(1.9912925, aa_masses[rev(aas)])
  } else if (!(is_empty(ntmod) || is_empty(ctmod))) { # anywhere, both N and C-term mods
    ions <- c(aa_masses[names(ctmod)] - 15.0114475, aa_masses[rev(aas)])
  } else if (!is_empty(ntmod)) { # anywhere, N-term
    ions <- c(1.9912925, aa_masses[rev(aas)])
  } else if (!is_empty(ctmod)) { # anywhere C-term
    ions <- c(aa_masses[names(ctmod)] - 15.0114475, aa_masses[rev(aas)])
  }
  
  ions %>% 
    cumsum() %>% 
    round(digits = digits)
}

