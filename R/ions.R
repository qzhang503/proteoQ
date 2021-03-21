#' Function factories for finding the amino-acid residue at a position.
#'
#' If the the value in a \code{vmods} is one of the LETTERS -> a site is
#' specified in the corresponding variable modification. Otherwise, a site is in
#' \code{c("C-term", "N-term")}.
#'
#' @param pos The Unimod position of a variable modification.
#' @return A function for finding the residue at the position specified by the
#'   argument \code{pos}. For each function, it takes a list of variable
#'   modifications specified by argument \code{vmods} as inputs.
#' @export
find_pos_site <- function (pos) {
  force(pos)
  
  stopifnot(pos %in% c("Protein N-term", "Protein C-term", 
                       "Any N-term", "Any C-term", 
                       "Anywhere"))
  
  pos <- paste0("^", pos)
  
  function (vmods) {
    idxes <- vmods %>% 
      purrr::map_lgl(~ {
        grepl(pos, names(.x)) & .x %in% LETTERS
      }) 
    
    vmods[idxes]
  }
}


#' Finds the sites of amino-acid residues at corresponding variable
#' modifications.
#'
#' Flowchart (1-nt): `Dimethyl (Protein N-term = P)`
#'
#' @param vmods A named list of variable modifications.
#' @seealso contain_protntany
#' @examples
#' # `Protein N-term = P`
#' sites <- list(`Dimethyl (Protein N-term = P)` = "P", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_protntsite(vmods)
#' find_protntsite(vmods)
#' @export
find_protntsite <- find_pos_site("Protein N-term")


#' Finds the sites of amino-acid residues at corresponding variable
#' modifications.
#' 
#' Flowchart (3-nt): `Gln->pyro Glu (N-term = Q)`
#' 
#' @rdname  find_protntsite
#' @examples 
#' # Gln->pyro Glu (N-term = Q)
#' sites <- list(`Gln->pyro Glu (N-term = Q)` = "Q", 
#'               `Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_anyntsite(vmods)
#' find_anyntsite(vmods)
#' @export
find_anyntsite <- find_pos_site("Any N-term")


#' Finds aa sites at given variable modifications.
#' 
#' Flowchart (5): `Oxidation (M)`
#' 
#' @rdname  find_protntsite
#' 
#' @examples 
#' # `Oxidation (M)` and `Deamidated (N)`
#' sites <- list(`Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_anysite(vmods)
#' find_anysite(vmods)
#' @export
find_anysite <- find_pos_site("Anywhere")


#' Finds aa sites at given variable modifications.
#' 
#' Flowchart f(1-ct): `Dehydrated (Protein C-term = N)`
#' 
#' @rdname  find_protntsite
#' 
#' @examples 
#' # `Dehydrated (Protein C-term = N)`
#' sites <- list(`Dehydrated (Protein C-term = N)` = "N", 
#'               `Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein C-term", "Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_protctsite(vmods)
#' find_protctsite(vmods)
#' @export
find_protctsite <- find_pos_site("Protein C-term")


#' Finds aa sites at given variable modifications.
#' 
#' Flowchart f(3-ct): `Oxidation (C-term = G)`
#' 
#' @rdname  find_protntsite
#' 
#' @examples 
#' # `Oxidation (C-term = G)`
#' sites <- list(`Oxidation (C-term = G)` = "G", 
#'               `Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any C-term", "Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_anyctsite(vmods)
#' find_anyctsite(vmods)
#' @export
find_anyctsite <- find_pos_site("Any C-term")



#' Function factories for checking the existence of an amino-acid residue at a
#' position.
#' 
#' If the the value in a \code{vmods} is one of the LETTERS -> a site is
#' specified in the corresponding variable modification. Otherwise, a site is in
#' \code{c("C-term", "N-term")}.
#'
#' @inheritParams find_pos_site
#' @return A function for checking the existence of a residue at the position
#'   specified by the argument \code{pos}. For each function, it takes a list of
#'   variable modifications specified by argument \code{vmods} as inputs.
contain_pos_site <- function (pos) {
  force(pos)
  
  stopifnot(pos %in% c("Protein N-term", "Protein C-term", 
                       "Any N-term", "Any C-term", 
                       "Anywhere"))
  
  pos <- paste0("^", pos)

  function (vmods) {
    if (length(vmods) == 0) return(FALSE)
    if (length(vmods) == 1 && vmods == "") return(FALSE)
    
    vmods %>% 
      purrr::map_lgl(~ {
        grepl(pos, names(.x)) & .x %in% LETTERS 
      }) %>% 
      any()
  }
}


#' Flowchart (1-nt): `Dimethyl (Protein N-term = P)`
#'
#' @rdname find_protntsite
contain_protntsite <- contain_pos_site("Protein N-term")


#' Flowchart (3-nt): `Gln->pyro Glu (N-term = Q)`
#' 
#' @rdname  find_protntsite
contain_anyntsite <- contain_pos_site("Any N-term")


#' Flowchart (5): `Oxidation (M)`
#' 
#' @rdname find_protntsite
contain_anysite <- contain_pos_site("Anywhere")


#' Flowchart (1-ct): 'Dehydrated (Protein C-term = N)'
#' 
#' @rdname find_protntsite
contain_protctsite <- contain_pos_site("Protein C-term")


#' Flowchart (3-ct): 'Oxidation (C-term = G)'
#' 
#' @rdname find_protntsite
contain_anyctsite <- contain_pos_site("Any C-term")



#' Function factories for checking the existence of an amino-acid residue at a
#' \code{terminal} position.
#'
#' A \code{pos} must be terminal in one of \code{c("Protein N-term", "Any
#' N-term", "Protein C-term", "Any C-term")}.
#'
#' @param pos The position. It must be terminal in \code{c("Protein N-term",
#'   "Any N-term", "Protein C-term", "Any C-term")}.
#' @return A function for checking the existence of a residue at the
#'   \code{terminal} position specified by the argument \code{pos}. For each
#'   function, it takes a list of variable modifications specified by argument
#'   \code{vmods} as inputs.
contain_termpos_any <- function (pos) {
  force(pos)

  stopifnot(pos %in% c("Protein N-term", "Any N-term", 
                       "Protein C-term", "Any C-term"))

  pos <- paste0("^", pos)

  function (vmods) {
    if (length(vmods) == 0) return(FALSE)
    if (length(vmods) == 1 && vmods == "") return(FALSE)
    
    vmods %>% 
      purrr::map_lgl(~ {
        grepl(pos, names(.x)) 
      }) %>% 
      any()
  }
}


#' Checks if variable modifications are on terminal positions.
#'
#' Flowchart (2-nt): `Acetyl (Protein N-term)`
#'
#' @inheritParams find_nmodtree
#'
#' @examples
#' # `Acetyl (Protein N-term)`
#' sites <- list(`Acetyl (Protein N-term)` = "N-term", `Oxidation (M)` = "M", `Deamidated (N)` = "N")
#' positions <- c("Protein N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_protntany(vmods)
#' @export
contain_protntany <- contain_termpos_any("Protein N-term")


#' Checks if variable modifications are on terminal positions.
#' 
#' Flowchart (4-nt): `Acetyl (N-term)`
#' 
#' @rdname contain_protntany
#' @examples
#' # `Acetyl (N-term)`
#' sites <- list(`Acetyl (N-term)` = "N-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any N-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_anyntany(vmods)
#' @export
contain_anyntany <- contain_termpos_any("Any N-term")


#' Checks if variable modifications are on terminal positions.
#' 
#' Flowchart (2-ct): `Amidated (Protein C-term)`
#' 
#' @rdname contain_protntany
#' @examples
#' # `Amidated (Protein C-term)`
#' sites <- list(`Amidated (Protein C-term)` = "C-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Protein C-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_protctany(vmods)
#' @export
contain_protctany <- contain_termpos_any("Protein C-term")


#' Checks if variable modifications are on terminal positions.
#' 
#' Flowchart (4-ct): `Amidated (C-term)`
#' 
#' @rdname contain_protntany
#' @examples
#' # `Amidated (C-term)`
#' sites <- list(`Amidated (C-term)` = "C-term", 
#'               `Oxidation (M)` = "M", 
#'               `Deamidated (N)` = "N")
#' positions <- c("Any C-term", "Anywhere", "Anywhere")
#' vmods <- purrr::map2(sites, positions, ~ setNames(.x, .y))
#' contain_anyctany(vmods)
#' @export
contain_anyctany <- contain_termpos_any("Any C-term")



#' Subsets proteins by variable modifications.
#'
#' Flowchart (1-nt): `Dimethyl (Protein N-term = P)`
#'
#' @param sites A list of sites. No need to specify for any terminal sites.
#' @inheritParams find_nmodtree
subset_protntsite <- function (peps, sites = "M") {
  peps %>% 
    purrr::map(~ {
      peps_i <- .x
      
      # maybe use `reduce(sites, f)`
      purrr::walk(sites, ~ {
        peps_i <<- peps_i %>% 
          .[grepl(paste0("^-", .x), .)]
      })
      
      peps_i
    })
  
  # peps %>% purrr::map(~ .x %>% .[grepl(paste0("^-", sites), .)])
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (2-nt): `Acetyl (Protein N-term)`
#' 
#' @rdname subset_protntsite
subset_protntany <- function (peps) {
  peps %>% 
    purrr::map(~ .x %>% .[grepl("^-", .)])
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (3-nt): `Gln->pyro Glu (N-term = Q)`
#' @rdname subset_protntsite
subset_anyntsite <- function (peps, sites = "Q") {
  peps %>% 
    purrr::map(~ {
      peps_i <- .x
      
      purrr::walk(sites, ~ {
        peps_i <<- peps_i %>% 
          .[grepl(paste0("^", .x), .)]
      })
      
      peps_i
    })

  # peps %>% purrr::map(~ .x %>% .[grepl(paste0("^", sites), .)])
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart f(4-nt): `Acetyl (N-term)`
#' @rdname subset_protntsite
subset_anyntany <- function (peps) {
  peps 
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (5): `Oxidation (M)`
#' @rdname subset_protntsite
subset_anysite <- function (peps, sites = "N") {
  peps %>% 
    purrr::map(~ {
      peps_i <- .x
      
      purrr::walk(sites, ~ {
        peps_i <<- peps_i %>% .[grepl(.x, .)]
      })
      
      peps_i
    })
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (1-ct): `Dehydrated (Protein C-term = N)`
#' @rdname subset_protntsite
subset_protctsite <- function (peps, sites = "V") {
  peps %>% 
    purrr::map(~ {
      peps_i <- .x
      
      purrr::walk(sites, ~ {
        peps_i <<- peps_i %>% 
          .[grepl(paste0(.x, "-$"), .)]
      })
      
      peps_i
    })
  
  # peps %>% purrr::map(~ .x %>% .[grepl(paste0(sites, "-$"), .)])
}

# f(2-ct): 
subset_protctany <- function (peps) {
  peps %>% 
    purrr::map(~ .x %>% .[grepl("-$", .)])
}

#' Subsets proteins by variable modifications.
#' 
#' Flowchart (3-ct): `Oxidation (C-term = G)`
#' @rdname subset_protntsite
subset_anyctsite <- function (peps, sites = "Q") {
  peps %>% 
    purrr::map(~ {
      peps_i <- .x
      
      purrr::walk(sites, ~ {
        peps_i <<- peps_i %>% 
          .[grepl(paste0(.x, "$"), .)]
      })
      
      peps_i
    })
  
  # peps %>% purrr::map(~ .x %>% .[grepl(paste0(sites, "$"), .)])
}


#' Subsets proteins by variable modifications.
#' 
#' Flowchart (4-ct): `Amidated (C-term)`
#' @rdname subset_protntsite
subset_anyctany <- function (peps) {
  peps 
}



#' Find and subset peptides.
#'
#' (1-nt) Protein N-term + site -> (2-nt) Any protein N-term -> (3-nt) Any
#' N-term + site -> (4-nt) Any N-term -> (5) Anywhere + site.
#'
#' @param vmods A named list of variable modifications.
#' @inheritParams concat_peps
find_nmodtree <- function (peps, vmods) {
  if (contain_protntsite(vmods)) { # level_1: Protein N-term + Site
    peps <- peps %>% 
      subset_protntsite(find_protntsite(vmods))
    if (contain_anysite(vmods)) {
      # (1) -|* .. |
      peps <- peps %>% 
        subset_anysite(find_anysite(vmods))
    } else {
      # (2) -|*    |
      peps <- peps
    }
  } else {
    if (contain_protntany(vmods)) { # level_2: Protein N-term 
      peps <- peps %>% 
        subset_protntany()
      if (contain_anysite(vmods)) {
        # (3) -|o .. |
        peps <- peps %>% 
          subset_anysite(find_anysite(vmods))
      } else {
        # (4) -|o    |
        peps <- peps
      }
    } else { 
      if (contain_anyntsite(vmods)) { # level_3: Any N-term + Site
        peps <- peps %>% 
          subset_anyntsite(find_anyntsite(vmods))
        if (contain_anysite(vmods)) {
          # (5) |* .. |
          peps <- peps %>% 
            subset_anysite(find_anysite(vmods))
        } else {
          # (6) |*    |
          peps <- peps
        }
      } else { 
        if (contain_anyntany(vmods)) { # level_4: Any N-term 
          peps <- peps 
          if (contain_anysite(vmods)) {
            # (7) |o .. |
            peps <- peps %>% 
              subset_anysite(find_anysite(vmods))
          } else {
            # (8) |o    |
            peps <- peps
          }
        } else { 
          if (contain_anysite(vmods)) { # level_5: Anywhere
            # (9) |  .. |
            peps <- peps %>% 
              subset_anysite(find_anysite(vmods))
          } else {
            # (10) |     |
            peps <- peps
          }
        }
      }
    }
  }
}


#' Find and subset peptides.
#'
#' (1-ct) Protein C-term + site -> (2-ct) Any protein C-term -> (3-ct) Any
#' C-term + site -> (4-ct) Any C-term.
#'
#' @inheritParams find_nmodtree
find_cmodtree <- function (peps, vmods) {
  if (contain_protctsite(vmods)) { # level_1: Protein C-term + Site
    # (1) -|* .. *|-, (2) -|*    *|-, (3) -|o .. *|-, (4) -|o    *|-, (5) |* .. *|-, 
    # (6) |*    *|-, (7) |o .. *|-, (8) |o    *|-, (9) |  .. *|-, (10) |     *|-
    peps <- peps %>% 
      subset_protctsite(find_protctsite(vmods))
  } else {
    if (contain_protctany(vmods)) { # level_2: Protein C-term 
      # (1) -|* .. o|-, (2) -|*    o|-, (3) -|o .. o|-, # (4) -|o    o|-, # (5) |* .. o|-, 
      # (6) |*    o|-, (7) |o .. o|-, (8) |o    o|-, (9) |  .. o|-, (10) |     o|-
      peps <- peps %>% 
        subset_protctany()
    } else {
      if (contain_anyctsite(vmods)) { # level_3: Any C-term + Site
        # (1) -|* .. *|, (2) -|*    *|, (3) -|o .. *|, # (4) -|o    *|, # (5) |* .. *|, 
        # (6) |*    *|, (7) |o .. *|, (8) |o    *|, (9) |  .. *|, (10) |     *|
        peps <- peps %>% 
          subset_anyctsite(find_anyctsite(vmods))
      } else {
        if (contain_anyctany(vmods)) { # level_4: Any C-term
          # (1) -|* .. o|, (2) -|*    o|, (3) -|o .. o|, # (4) -|o    o|, # (5) |* .. o|, 
          # (6) |*    o|, (7) |o .. o|, (8) |o    o|, (9) |  .. o|, (10) |     o|
          peps <- peps 
        } # else {
          # No use
        # }
      }
    }
  }
  
  invisible(peps)
}


#' Subset peptides by variable modifications.
#' 
#' From N-term to C-term. 
#' 
#' @inheritParams concat_peps
#' @inheritParams add_fixvar_masses
subpeps_by_vmods <- function(aa_masses, peps) {
  vmods <- attr(aa_masses, "vmods_ps") 
  
  if (is.list(vmods)) purrr::flatten_chr(vmods) else vmods
  
  peps <- peps %>% 
    find_nmodtree(vmods) %>% 
    find_cmodtree(vmods) 
}

