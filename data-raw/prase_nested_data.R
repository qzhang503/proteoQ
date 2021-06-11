#' Finds an entry from the results of ion matches with multiple peptides and
#' nested varmod positions.
foo_find_peps <- function () {
  # res <- readRDS(file.path("~/proteoQ/outs", "ion_matches_8.rds"))
  res <- readRDS(file.path("~/proteoQ/outs", "ion_matches.rds"))[[7]]

  # find the nested case
  mts = res$matches
  lens <- map(mts, length)
  rows <- which(lens > 1)
  mts <- mts[rows]
  
  x = unlist(mts, recursive = F)
  lens <- map(x, length)
  rows <- which(lens > 1)
  pep <- names(x[rows[1]])
  # TNLAMMR
  
  pep <- names(x[rows[2]])
  # IILQMMR
  
  # ---
  mts = res$matches
  
  rows <- map_lgl(mts, ~ {
    any(names(.x) == pep)
  })
  
  # 63
  mts <- mts[rows]
  
  # which(names(mts) == x[rows[1]])
}


foo_find_ms2ions <- function () {
  i <- 8
  aa_masses <- aa_masses_all[[i]]
  
  theopeps <- readRDS(file.path(out_path, "pepmasses/", 
                                paste0("binned_theopeps_", i, ".rds")))
  
  theopeps <- theopeps %>% 
    purrr::map(~ {
      rows <- !duplicated(.x$pep_seq)
      .x[rows, c("pep_seq", "mass")]
    })
  
  # off by one frame
  theopeps <- theopeps[names(theopeps) %in% c(36934:36936)]
  # theopeps <- theopeps[names(theopeps) == 36935]
  
  # in `search_mgf_frames`
  # TNLAMMR
  theos_cr_ms2[15] 
  # $TNLAMMR --- Level 1
  # $TNLAMMR$`0000500` --- Level 2
  # $TNLAMMR$`0000500`$bs --- Level 3 (List of 12: bs, ys etc.)
  # ...
  # $TNLAMMR$`0000050` --- Level 2
  # $TNLAMMR$`0000050`$bs --- Level 3 (List of 12: bs, ys etc.)
  # ...
}


#' Finds peptide DDDIAALVVDNGSGMCK
find_DDDIAALVVDNGSGMCK <- function () {
  theopeps <- readRDS(file.path(out_path, "pepmasses/", 
                                paste0("binned_theopeps_", i, ".rds")))
  
  rows <- map_lgl(theopeps, ~ {
    x <- which(.x$pep_seq == "DDDIAALVVDNGSGMCK")
    !is_empty(x)
  }) %>% 
    which()

  # 68951 --- frame name
  # 12015 --- entry index
  
  # theopeps <- theopeps[names(theopeps) == names(rows)]
  theopeps <- theopeps[names(theopeps) == 68951]
  
  # ...
  
  x <- search_mgf_frames_d(mgf_frames[[1]], theopeps[[1]], aa_masses = aa_masses, 
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
                           digits = digits)
}


foo_byions <- function () {
  # ...
  assign("aa_masses_all", aa_masses_all, envir = .GlobalEnv)
  
  aa_masses <- aa_masses_all[[8]]
  
  aas <- str_split("HQGVMVGMGQK", "", simplify = T)
  digits <- 5
  
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
  
  maxn_vmods_setscombi = 64
  maxn_vmods_per_pep = 5
  maxn_sites_per_vmod = 3
  maxn_vmods_sitescombi_per_pep = 32
  
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
  
  # ---
  mass <- calc_monopeptide("HQGVMVGMGQK", aa_masses)
  mass <- mass[[1]]
  mass <- unname(mass)
  type_ms2ions <- "by"
  
  out <- vmods_combi %>% 
    map(calc_mvmods_ms2masses, ntmod, ctmod, aas, mass = mass, aa_masses, 
        type_ms2ions, digits)
  
  # ---
  vmods = vmods_combi[[1]]
    aas <- str_split("HQGVMVGMGQK", "", simplify = T)
  aas[as.numeric(names(vmods))] <- unlist(vmods)
  
  bs <- calc_bions(ntmod, ctmod, aas, aa_masses, digits = 5)
  ys <- calc_yions(ntmod, ctmod, aas, aa_masses, digits = 5)
  bs + rev(ys)
  
  add_complement_ions(mass, bs)
}


foo_neulosses <- function () {
  neulosses <- list(`Oxidation (M)` = c(0, 63.998285), 
                    `Carbamidomethyl (M)` = c(0, 105.024835))
  
  neulosses <- list(`TMT6plex (K)` = c(0, 50, 100), 
                    `Acetyl (Protein N-term)` = c(0, 80))
  
  fmods_nl <- list(`N-term` = c(0, 50, 100))
  vmods_nl <- list(`Oxidation (M)` = c(0, 63.998285), 
                   `Carbamidomethyl (M)` = c(0, 105.024835))
}

foo_combis <- function () {
  # HQGVMNMVGMGQK
  
  # Browse[2]> intra_combis
  # $M
  # $M[[1]]
  # [1] "Oxidation (M)"       "Carbamidomethyl (M)"
  # 
  # $M[[2]]
  # [1] "Oxidation (M)"       "Carbamidomethyl (M)" "Carbamidomethyl (M)"
  # 
  # $M[[3]]
  # [1] "Oxidation (M)"       "Carbamidomethyl (M)" "Oxidation (M)"      
  # 
  # 
  # $N
  # $N[[1]]
  # [1] "Deamidated (N)"
  
  
}