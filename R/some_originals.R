#' Calculates the masses of MS2 ion series.
#'
#' For a given type of fragmentation. Minimal error handling for speeds.
#' 
#' @param mass The mass of a theoretical MS1.
#' @param maxn_vmods_sitescombi_per_pep Integer; the maximum number of
#'   combinatorial variable modifications per peptide sequence.
#' @param type_ms2ions Character; the type of
#'   \href{http://www.matrixscience.com/help/fragmentation_help.html}{ MS2
#'   ions}. Values are in one of "by", "ax" and "cz".
#' @param z Integer; the charge state of an MS1 ion. Currently only for positive
#'   charges (not currently used).
#' @inheritParams calc_monopep
#' @import purrr
calc_ms2ionseries_orig <- function (aa_seq = NA, 
                               mass = NA, 
                               aa_masses = NA, 
                               mod_indexes = NULL, 
                               type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5, 
                               maxn_sites_per_vmod = 3, 
                               maxn_vmods_sitescombi_per_pep = 32, 
                               digits = 5) {
  
  options(digits = 9)
  
  run_scripts <- FALSE
  if (run_scripts) {
    vmods_ps <- aa_masses %>% attributes() %>% `[[`("vmods_ps")
    
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
  }
  
  vmods_ps <- attr(aa_masses, "vmods_ps", exact = TRUE)
  amods <- attr(aa_masses, "amods", exact = TRUE)
  tmod <- attr(aa_masses, "tmod", exact = TRUE)
  
  # multiple mods to [NC]-term already excluded from aa_masses
  ntmod <- tmod %>% .[. == "N-term"]
  ctmod <- tmod %>% .[. == "C-term"]
  
  calc_ms2ions(aa_seq = aa_seq, 
               mass = mass, 
               aa_masses = aa_masses, 
               mod_indexes = mod_indexes, 
               vmods_ps = vmods_ps, 
               amods = amods, 
               tmod = tmod, 
               ntmod = ntmod, 
               ctmod = ctmod, 
               maxn_vmods_per_pep = maxn_vmods_per_pep, 
               maxn_sites_per_vmod = maxn_sites_per_vmod, 
               type_ms2ions = type_ms2ions, 
               digits = digits)
}


#' Helper With position combinations
#'
#' @inheritParams calc_monopep
#' @inheritParams calc_ms2ionseries
#' @import purrr
#' @importFrom stringr str_split
calc_ms2ions_orig <- function (aa_seq, mass, aa_masses, mod_indexes, 
                          type_ms2ions = "by", 
                          vmods_ps, amods, tmod, ntmod, ctmod, 
                          maxn_vmods_per_pep = 3, 
                          maxn_sites_per_vmod = 5, 
                          maxn_vmods_sitescombi_per_pep = 32, 
                          digits = 5) {
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
  
  options(digits = 9)
  
  if (is.na(aa_seq)) {
    return(NULL)
  }
  
  aas <- aa_seq %>% 
    str_split("", simplify = TRUE)
  
  if (is_empty(amods) && is_empty(tmod)) {
    out <- calc_mvmods_ms2masses(vmods = NULL, ntmod = NULL, ctmod = NULL, 
                                 aas = aas, mass = mass, aa_masses = aa_masses, 
                                 type_ms2ions = type_ms2ions, digits = digits) 
    
    nm <- rep(0, length(aas)) %>% paste0(collapse = "")
    out <- list(out) %>% `names<-`(nm)
    
    return(out)
  } else if (is_empty(amods)) { 
    # tmod only
    out <- calc_mvmods_ms2masses(vmods = NULL, ntmod = ntmod, ctmod = ctmod, 
                                 aas = aas, mass = mass, aa_masses = aa_masses, 
                                 type_ms2ions = type_ms2ions, digits = digits) 
    
    nm <- rep(0, length(aas)) %>% paste0(collapse = "")
    out <- list(out) %>% `names<-`(nm)
    
    return(out)
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
    }
    
    v_out
  })
  
  out <- vmods_combi %>% 
    map(calc_mvmods_ms2masses, ntmod, ctmod, aas, mass = mass, aa_masses, 
        type_ms2ions, digits)
  
  # assign indexes to mods
  # (no need for terminal mods: non additive and unambiguously from aa_masses)
  hex_mods = rep("0", length(aas))
  rows <- map_lgl(out, ~ !is.null(.x))
  
  hex_mods <- map(vmods_combi[rows], ~ {
    hex_mods[as.numeric(names(.x))] <- mod_indexes[unlist(.x)]
    hex_mods <- paste0(hex_mods, collapse = "")
    hex_mods
  }, hex_mods)
  
  names(out)[rows] <- hex_mods
  
  invisible(out[rows])
}




