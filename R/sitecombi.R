#' With position combinations
#'
#' For future usage.
#' @param maxn_vmods_sitescombi_per_pep Integer; the maximum number of
#'   combinatorial variable modifications per peptide sequence.
#' @inheritParams calc_monopep
calc_monopep2 <- function (aa_seq, aa_masses, vmods_ps, amods, tmod, ntmod, ctmod, 
                           mod_indexes = NULL, maxn_vmods_per_resid = 3, 
                           maxn_vmods_sitescombi_per_pep = 64, digits = 5) {
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # oxygen <- 15.99491462
  
  # moverz <- (mass + h2o + z*proton)/z
  
  options(digits = 9)
  
  if (is.na(aa_seq)) {
    return(NULL)
  }
  
  # no "-" in `aa_seq`
  aas <- aa_seq %>% 
    stringr::str_split("", simplify = TRUE)
  
  if (purrr::is_empty(amods) && purrr::is_empty(tmod)) {
    return(
      aas %>% 
        aa_masses[.] %>% 
        sum() %>% 
        `+`(aa_masses["N-term"]) %>% 
        `+`(aa_masses["C-term"]) %>% 
        setNames(aa_seq) %>% 
        round(digits = digits)
    )
  } else if (purrr::is_empty(amods)) { # tmod only
    return(
      calc_mvmods_masses(vmods = NULL, amods, ntmod, ctmod, aa_seq, aas, aa_masses, digits)
    )
  }
  
  # ---
  # amods and/or tmods
  vmods_combi <- local({
    intra_combis <- combi_mvmods(amods, aa_seq, aa_masses, maxn_vmods_per_resid, digits)
    
    # expand outs
    p_combis <- purrr::map(intra_combis, ~ {
      x <- .x
      purrr::map(x, names)
    }) 
    
    v_combis <- purrr::map(intra_combis, ~ {
      x <- .x
      purrr::map(x, purrr::reduce, `c`)
    })
    
    if (length(p_combis) > 1L) {
      max_combi <- purrr::map(p_combis, ~ pmin(length(.x), maxn_vmods_sitescombi_per_pep))
      p_combis <- purrr::map2(p_combis, max_combi, ~ .x[seq_len(.y)])
      v_combis <- purrr::map2(v_combis, max_combi, ~ .x[seq_len(.y)])
      
      p_combis <- p_combis %>% 
        purrr::reduce(expand.grid, KEEP.OUT.ATTRS = FALSE)
      
      v_combis <- v_combis %>% 
        purrr::reduce(expand.grid, KEEP.OUT.ATTRS = FALSE)
      
      nrow <- nrow(v_combis)
      p_out <- v_out <- vector("list", nrow)
      
      for (i in 1:nrow) {
        p_out[[i]] <- reduce(p_combis[i, ], `c`) %>% unlist
        v_out[[i]] <- reduce(v_combis[i, ], `c`) %>% unlist
        names(v_out[[i]]) <- p_out[[i]]
        v_out[[i]] <- v_out[[i]][order(as.numeric(names(v_out[[i]])))]
      }
    } else {
      v_out <- purrr::map2(p_combis[[1]], v_combis[[1]], ~ {
        names(.y) <- .x
        .y
      })
    }
    
    v_out
  })
  
  # i.e., btw Anywhere "M" and "Acetyl N-term" where "M" on the "N-term"
  vmods_combi <- local({
    if (!purrr::is_empty(ntmod)) {
      # the "start" may be "2", instead of "1"
      vmods_combi <- purrr::map(vmods_combi, ~ {
        if (names(.x[1]) == "1") .x <- NULL
        .x
      })
      
      rows <- purrr::map_lgl(vmods_combi, is.null)
      
      vmods_combi <- vmods_combi[!rows]
    }
    
    vmods_combi
  })
  
  vmods_combi <- local({
    if (!purrr::is_empty(ctmod)) {
      len <- length(aas)
      # if (aas[len] == "-") { len <- len - 1 }
      len <- as.character(len)
      
      vmods_combi <- purrr::map(vmods_combi, ~ {
        if (names(.x[length(.x)]) == len) {
          .x <- NULL
        }
        
        .x
      })
      
      rows <- purrr::map_lgl(vmods_combi, is.null)
      
      vmods_combi <- vmods_combi[!rows]
    }
    
    vmods_combi
  })
  
  if (length(vmods_combi) > maxn_vmods_sitescombi_per_pep) {
    vmods_combi <- vmods_combi[seq_len(maxn_vmods_sitescombi_per_pep)]
  }
  
  # the same peptide at different combinations are shown as different entries
  # may be remove the duplication
  vmods_combi %>% 
    purrr::map(calc_mvmods_masses2, ntmod, ctmod, aa_seq, aas, 
               aa_masses, digits)
}


#' The combinations of variable modifications (multiple sites)
#'
#' For all the \code{Anywhere} modifications specified in \code{amods}.
#'
#' @param amods \code{Anywhere} modifications.
#' @inheritParams calc_monopep
#' @inheritParams mcalc_monopep
#' @return Lists by residues in \code{amods}.
combi_mvmods <- function (amods, aa_seq, aa_masses, maxn_vmods_per_resid, digits) {
  residue_mods <- do.call(rbind, amods) %>% 
    data.frame() %>% 
    split(., .[["Anywhere"]])
  
  purrr::map(residue_mods, 
             ~ combi_vmods(aa_seq, .x, aa_masses, maxn_vmods_per_resid, digits))
}


#' The combinations of variable modifications (single site)
#'
#' @param residue_mods A residue with \code{Anywhere} modification(s).
#' @inheritParams combi_mvmods
combi_vmods <- function (aa_seq, residue_mods, aa_masses, maxn_vmods_per_resid, digits) {
  # n = LETTERS[1:2]; p = c(1, 3, "e")
  # 2*3, 4*3, 8*1
  # l = length(p)
  # n^1 * combn(p, 1) + n^2 * combn(p, 2) + ... + n^l * combn(p, l) 
  
  residue <- residue_mods[1, 1]
  
  n <- rownames(residue_mods)
  p <- stringr::str_locate_all(aa_seq, pattern = residue) %>% 
    `[[`(1) 
  p <- p[, 1] %>% as.character()
  
  ns <- purrr::map(seq_along(p), ~ {
    expand.grid(rep(list(n), length(p[1:.x])), KEEP.OUT.ATTRS = FALSE)
  }, n, p)
  
  ps <- purrr::map(seq_along(p), ~ {
    combn(p, .x)
  }, p)
  
  # values: n (modifications)
  # names: p (positions)
  out <- purrr::map2(ns, ps, ~ {
    lns <- nrow(.x)
    lps <- ncol(.y)
    np <- vector("list", lns * lps)
    
    k <- 1
    for (i in seq_len(lns)) {
      for (j in seq_len(lps)) {
        x <- .x[i, ] 
        x <- sapply(x, as.character)
        names(x) <- .y[, j]
        np[[k]] <- x
        
        k <- k + 1
      }
    }
    
    np
  }) %>% 
    purrr::flatten() 
  
  rows <- purrr::map_lgl(out, ~ length(.x) > maxn_vmods_per_resid)
  out[!rows]
}


#' With position combinations
#' 
#' For future usage.
#' @inheritParams calc_mvmods_masses
calc_mvmods_masses2 <- function (vmods, ntmod, ctmod, aa_seq, aas, aa_masses, digits) {
  # unimods <- vmods %>% purrr::map(parse_unimod) 
  # sites <- unimods %>% purrr::map_chr(`[[`, "site")
  
  idxes <- names(vmods) %>% as.numeric()
  aas[idxes] <- vmods
  
  out <- aas %>% 
    aa_masses[.] %>% 
    sum() 
  
  if (is_empty(ntmod) && is_empty(ctmod)) {
    out <- out %>% 
      `+`(aa_masses["N-term"]) %>% 
      `+`(aa_masses["C-term"]) 
  } else if (!(is_empty(ntmod) || is_empty(ctmod))) {
    out <- out %>% 
      `+`(aa_masses[names(ntmod)]) %>% 
      `+`(aa_masses[names(ctmod)])
  } else if (!is_empty(ntmod)) {
    out <- out %>% 
      `+`(aa_masses[names(ntmod)]) %>% 
      `+`(aa_masses["C-term"]) 
  } else if (!is_empty(ctmod)) {
    out <- out %>% 
      `+`(aa_masses["N-term"]) %>% 
      `+`(aa_masses[names(ctmod)])
  }
  
  out %>% 
    setNames(aa_seq) %>% 
    round(digits = digits)
}

