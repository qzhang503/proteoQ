#' Reads a file in fasta format
#'
#' Reads a file in fasta format by line.
#'
#' @param file A character string to the name of a protein fasta file.
#' @param pattern A regular expression describing the pattern to separate the
#'   header lines of fasta entries. The default is to separate a header and keep
#'   the character string before the first space where the so kept will be used
#'   as the name of an entry. The character ">" at the beginning of the header
#'   will also be removed.
#' @param comment_char Character: a character or an empty string. Use "" to turn
#'   off the interpretation of comment lines.
#' @examples
#' \donttest{
#' # assume the file and location of "uniprot_hs_2020_05.fasta"
#' fasta <- read_fasta("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' head(names(fasta))
#' 
#' # use the first fifty characters
#' fasta <- read_fasta("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
#'                     pattern = ">(.{50}).*")
#' head(names(fasta))
#' 
#' # uniprot_acc
#' fasta <- read_fasta("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
#'                     pattern = ">..\\|([^\\|]+)\\|.*")
#' head(names(fasta))
#' 
#' # use every in the header
#' fasta <- read_fasta("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
#'                     pattern = ">(.*)")
#' head(names(fasta))
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @seealso \code{\link{write_fasta}}
#' @export
read_fasta <- function (file, pattern = ">([^ ]+?) .*", comment_char = "") {
  lines <- readLines(file)
  
  if (nchar(comment_char) > 0) {
    lines <- lines %>% .[!grepl(paste0("^", comment_char), .)]
  }
  
  headers <- grep(">", lines)
  begins <- headers + 1
  ends <- (headers[-1] - 1) %>% `c`(length(lines))
  
  seqs <- purrr::map2(begins, ends, ~ lines[.x : .y] %>% purrr::reduce(paste0), lines)
  hdrs <- lines[headers]
  
  db <- purrr::map2(seqs, hdrs, ~ {
    attr(.x, "header") <- .y
    return(.x)
  })
  
  # maybe also first 80 characters...
  names(db) <- hdrs %>% gsub(pattern, "\\1", .)

  return(db)
}


#' Writes fasta
#'
#' Writes a fasta file.
#'
#' @param fasta_db A list of protein entries from \code{\link{read_fasta}}.
#' @inheritParams read_fasta
#' @examples
#' \donttest{
#' fasta_db <- read_fasta(file = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' write_fasta(fasta_db, "~/proteoQ/examples/my.fasta")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @export
write_fasta <- function (fasta_db, file) {
  filepath <- gsub("(^.*/).*$", "\\1", file)
  dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
  
  purrr::map(fasta_db, ~ paste(attr(.x, "header"), .x, sep = "\n")) %>% 
    unlist() %>% 
    writeLines(file)
}


#' Loads fasta
#' 
#' @inheritParams splitPSM
#' @examples
#' \donttest{
#' fasta_db <- load_fasta("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' }
load_fasta <- function (fasta = NULL) {
  if (is.null(fasta)) {
    stop("FASTA file(s) are required.", call. = FALSE)
  }
  
  if (!all(file.exists(fasta))) {
    stop("Missing FASTA file(s): \n", 
         purrr::reduce(fasta %>% .[!file.exists(.)], paste, sep = "\n"), 
         call. = FALSE)
  }
  
  purrr::map(fasta, ~ read_fasta(.x)) %>% 
    do.call(`c`, .) %>% 
    `names<-`(gsub(">", "", names(.))) %>% 
    .[!duplicated(names(.))]
}


#' Averaged peptide molecular weight (for proteins)
#'
#' Calculates the average molecular weight of a polypeptide (neutral species).
#' 
#' @inheritParams calc_monopep
#' @examples
#' \donttest{
#' calc_avgpep("AAIDWFDGKEFSGNPIK")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
calc_avgpep <- function (aa_seq, digits = 4) {
  options(digits = 9)
  
  # data(aa_residues)
  # aa_masses <- aa_residues$average_da %>% `names<-`(aa_residues$one_letter)
  
  aa_masses <- c(
    A = 71.0779, R = 156.1857, N = 114.1026, D = 115.0874, 
    C = 103.1429, E = 129.1140, Q = 128.1292, G = 57.0513, 
    H = 137.1393, I = 113.1576,L = 113.1576, K = 128.1723, 
    M = 131.1961, F = 147.1739, P = 97.1152, S = 87.0773, 
    T = 101.1039, W = 186.2099,Y = 163.1733, V = 99.1311, 
    "N-term" = 1.0079, "C-term" = 17.0073, 
    U = 150.0379, B = 114.5950, X = 111.0000, Z = 128.6216)
  
  aa_seq %>% 
    stringr::str_split("", simplify = TRUE) %>% 
    aa_masses[.] %>% 
    sum() %>% 
    `+`(18.010565) %>% 
    round(digits = digits) 
}


#' Calculates multiply the molecular weight of a polypeptides ([MH]+).
#'
#' The calculations iterate through protein accessions and peptide sequences.
#'
#' @param aa_seqs Character string; a vector of peptide sequences with
#'   one-letter representation of amino acids.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @param maxn_sites_per_vmod Integer; the maximum number of combinatorial
#'   variable modifications per site in a per peptide sequence.
#' @param maxn_vmods_per_pep The maximum number of variable modifications per
#'   peptide.
#' @inheritParams add_fixvar_masses
#' @inheritParams calc_pepmasses
#' @examples
#' \dontrun{
#' data(package = "proteoQ", aa_residues)
#' aa_masses <- aa_residues %>%
#'   dplyr::select(c("one_letter", "monoisotopic_da"))
#' aa_masses <- aa_masses$monoisotopic_da %>% `names<-`(aa_masses$one_letter)
#'
#' x <- mcalc_monopep(list(protein_a = c("AAIDWFDGKEFSGNPIK", "MAAIDWFDGKEFSGNPIK"),
#'                         protein_b = c("AAIDWFDGKEFSGNPIK")),
#'                    aa_masses)
#'
#' aa_masses_all <- calc_aamasses()
#'
#' # Error if peptide sequences not 'compatible' with the given 'aa_masses'
#' x <- mcalc_monopep(list(protein_a = c("AAIMDWFDMGKEFSGNPNIK", "MAAIDWFDGKEFSGNPNIK"),
#'                         protein_b = c("AAIDWFDGKEMFSGNPNIK")),
#'                    aa_masses_all[[11]])
#' }
mcalc_monopep <- function (aa_seqs, aa_masses, mod_indexes, 
                           maxn_vmods_per_pep = 5, 
                           maxn_sites_per_vmod = 3, 
                           digits = 5) {
  vmods_ps <- aa_masses %>% 
    attributes() %>% 
    `[[`("vmods_ps")
  
  # multiple mods to [NC]-term already excluded from aa_masses
  amods <- local({
    sites <- vmods_ps %>% 
      purrr::map(~ .x[grepl("Anywhere", names(.x))]) 
    
    empties <- sites %>% purrr::map_lgl(purrr::is_empty)
    
    sites <- sites[!empties] 
    # names(sites) <- mod_indexes[names(sites)]
    sites
  })
  
  tmod <- vmods_ps %>% .[! . %in% amods]
  if (purrr::is_empty(tmod)) {
    tmod <- NULL
  } else if (tmod == "") {
    tmod <- NULL
  }
  
  ntmod <- tmod %>% .[. == "N-term"]
  ctmod <- tmod %>% .[. == "C-term"]
  
  attr(aa_masses, "data") <- aa_seqs %>% 
    purrr::map(~ { # by prot_acc
      .x %>% 
        purrr::map(calc_monopep, 
                   aa_masses, vmods_ps, amods, tmod, ntmod, ctmod, 
                   mod_indexes, 
                   maxn_vmods_per_pep, 
                   maxn_sites_per_vmod, 
                   digits) %>% # by peptides
        unlist() 
    })
}


#' Helper: calculates the mono-isotopic mass of a peptide sequence.
#'
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param vmods_ps A named list of positions and sites with positions being names
#'   and sites being values.
#' @param amods Anywhere modifications subset from \code{vmods_ps}.
#' @param tmod A terminal modification subset from \code{vmod_ps}.
#' @param ntmod A N-term modification subset from \code{tmod}.
#' @param ctmod A C-term modification subset from \code{tmod}.
#' @inheritParams mcalc_monopep
calc_monopep <- function (aa_seq, aa_masses, vmods_ps, amods, tmod, ntmod, ctmod, 
                          mod_indexes = NULL, 
                          maxn_vmods_per_pep = 5, 
                          maxn_sites_per_vmod = 3, 
                          digits = 5) {

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
  } else if (purrr::is_empty(amods)) { 
    # tmod only
    return(
      calc_mvmods_masses(vmods = NULL, amods, ntmod, ctmod, 
                         aa_seq, aas, aa_masses, digits)
    )
  }

  # amods and/or tmods
  vmods_combi <- local({
    intra_combis <- unique_mvmods(amods = amods, aas = aas, aa_masses = aa_masses, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  digits = digits)
    
    if (length(intra_combis) > 1L) {
      inter_combi <- intra_combis %>% 
        purrr::reduce(expand.grid, KEEP.OUT.ATTRS = FALSE)
      
      nrow <- nrow(inter_combi)
      v_out <- vector("list", nrow)
      
      for (i in 1:nrow) {
        v_out[[i]] <- purrr::reduce(inter_combi[i, ], `c`) %>% 
          unlist()
      }
    } else {
      v_out <- purrr::flatten(intra_combis)
    }
    
    invisible(v_out)
  })
  
  # allows some 'conflicting' combinations for speed
  # i.e., btw Anywhere "M" and "Acetyl N-term" where "M" on the "N-term"
  # otherwise need to count the number of M's to decide the availability of Acetyl N-term

  vmods_combi %>% 
    purrr::map(calc_mvmods_masses, amods, ntmod, ctmod, aa_seq, aas, 
               aa_masses, digits)
}


#' Helper of the calculations of peptide masses.
#'
#' Currently, the redundancies are kept for the different combinations (i.e.
#' M[1]N[2], N[1]M[2]) of the same peptide. 
#' 
#' @param vmods A list of variable modifications. Site indexes in names;
#'   modifications in values.
#' @param ntmod The \code{N-term} modification.
#' @param ctmod The \code{C-term} modification.
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
#' @inheritParams calc_monopep
calc_mvmods_masses <- function (vmods, amods, ntmod, ctmod, 
                                aa_seq, aas, aa_masses, digits) {
  if (!purrr::is_empty(amods)) {
    nms <- purrr::imap(amods, ~ {
      names(.x) <- .y
      .x
    }) %>% purrr::flatten() %>% 
      unlist()
    
    n_per_site <- purrr::map(names(nms), ~ sum(vmods == .x)) %>% 
      `names<-`(names(nms)) %>% 
      unlist()
    
    for (i in seq_along(nms)) {
      x <- nms[i]
      y <- n_per_site[i]
      
      idxes <- which(aas == x)[1:y]
      aas[idxes] <- rep(names(y), y)
    }
    
    rm(i, x, y, nms, n_per_site)
  }

  out <- aas %>% 
    aa_masses[.] %>% 
    sum() 
  
  if (purrr::is_empty(ntmod) && purrr::is_empty(ctmod)) {
    out <- out %>% 
      `+`(aa_masses["N-term"]) %>% 
      `+`(aa_masses["C-term"]) 
  } else if (!(purrr::is_empty(ntmod) || purrr::is_empty(ctmod))) {
    out <- out %>% 
      `+`(aa_masses[names(ntmod)]) %>% 
      `+`(aa_masses[names(ctmod)])
  } else if (!purrr::is_empty(ntmod)) {
    out <- out %>% 
      `+`(aa_masses[names(ntmod)]) %>% 
      `+`(aa_masses["C-term"]) 
  } else if (!purrr::is_empty(ctmod)) {
    out <- out %>% 
      `+`(aa_masses["N-term"]) %>% 
      `+`(aa_masses[names(ctmod)])
  }
  
  out %>% 
    setNames(aa_seq) %>% 
    round(digits = digits)
}


#' The unique entries of variable modifications (multiple sites)
#'
#' For all the \code{Anywhere} modifications specified in \code{amods}.
#'
#' @param amods Anywhere modifications.
#' @inheritParams calc_monopep
#' @inheritParams mcalc_monopep
#' @inheritParams calc_mvmods_masses
#' @return Lists by residues in \code{amods}.
unique_mvmods <- function (amods, aas, aa_masses, 
                           maxn_vmods_per_pep = 5, 
                           maxn_sites_per_vmod = 3, 
                           digits = 5) {
  residue_mods <- do.call(rbind, amods) %>% 
    data.frame() %>% 
    split(., .[["Anywhere"]])
  
  purrr::map(residue_mods, 
             ~ vmods_elements(aas, .x, aa_masses, 
                              maxn_vmods_per_pep, 
                              maxn_sites_per_vmod, 
                              digits)) 
}


#' Find the sets of variable modifications.
#'
#' Excluding position differernces, i.e., \code{A, B} and \code{B, A} is the
#' same set.
#'
#' @param residue_mods Amino-acid residues with Unimod names. For example
#'   rownames of \code{Carbamidomethyl (M)} and \code{Oxidation (M)} and a
#'   column residues of \code{M, M}.
#' @inheritParams unique_mvmods
vmods_elements <- function (aas, 
                            residue_mods, aa_masses, 
                            maxn_vmods_per_pep = 5, 
                            maxn_sites_per_vmod = 3, 
                            digits = 5) {
  residue <- residue_mods[1, 1]
  
  ns <- rownames(residue_mods)
  ps <- which(aas == residue)
  
  len_n <- length(ns)
  len_p <- length(ps)

  if (len_p > len_n) {
    x <- purrr::map((length(ns) + 1):length(ps), 
                    ~ fnd_unique_sets(ps[seq_len(.x)], ns)) %>% 
      recur_flatten() %>% 
      `c`(list(ns), .)
  } else {
    x <- list(ns)
  }

  len_x <- purrr::map(x, length)
  rows <- purrr::map_lgl(x, ~ length(.x) > maxn_vmods_per_pep)
  x <- x[!rows]
  
  maxn_vmod <- x %>% 
    purrr::map(table) %>% 
    purrr::map(max)
  rows <- maxn_vmod > maxn_sites_per_vmod
  x <- x[!rows]
}


#' Finds the sets of \code{n} elements in \code{positions}.
#' 
#' At least one occupancy for each elements in \code{ns}.
#' 
#' @param ps A vector of positions.
#' @param ns The names to be filled into \code{p}.
#' @importFrom gtools combinations
fnd_unique_sets <- function (ps = c(1:5), ns = c("A", "B", "C")) {
  len_p <- length(ps)
  len_n <- length(ns)
  
  x <- gtools::combinations(len_n, len_p - len_n, ns, repeats = TRUE)
  
  purrr::map(1:nrow(x), ~ c(ns, x[.x, ]))
}


#' Helper to add modification masses to amino-acid residues.
#'
#' It adds the masses of fixed, variable, and neutral-loss modifications to
#' amino-acid residues.
#'
#' @param mods A list of modifications.
#' @param mod_type The type of modification in one of \code{c("fmods", "vmods")}
#'   where \code{fmods}: fixed modifications and \code{vmods}: variable
#'   modifications.
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @return Lists of of amino-acid residues with modified mono-isotopic masses
#'   being incorporated.
add_fixvar_masses <- function (mods, mod_type, aa_masses) {
  stopifnot(mod_type %in% c("fmods", "vmods"), 
            length(mod_type) == 1)
  
  if (!is.null(mods)) {
    all_mods <- purrr::reduce(mods, paste, sep = ", ")
  } else {
    all_mods <- ""
  }

  res <- mods %>% 
    purrr::map(find_unimod) %>% 
    `names<-`(mods)
  
  mod_masses <- res %>% purrr::map(`[[`, 1)
  positions_sites <- res %>% purrr::map(`[[`, 2)
  neulosses <- res %>% purrr::map(`[[`, 3)
  rm(res)
  
  # the same `site` with different fixedmods
  local({
    if (mod_type == "fmods" && length(positions_sites) > 1) {
      dups <- purrr::reduce(positions_sites, `c`) %>% 
        .[duplicated(.)]
      
      if (!purrr::is_empty(dups)) {
        dups_in_each <- positions_sites %>% 
          purrr::map(~ .x[.x == dups])
        
        warning("Conflicts in fixed modifications: \n", 
                purrr::reduce(names(dups_in_each), paste, sep = "\n"), "\n",
                "May consider change from fixed to variable modifications(s); \n",
                "or create a new Unimod for joint modifications.", 
                call. = FALSE)
      }
    }
  })

  if (mod_type == "fmods") {
    # Add mod_masses of fixed mods
    purrr::walk2(positions_sites, mod_masses, ~ {
      # 'N-term' not 'Q' for Gln-> pyro-Glu (N-term = Q)
      if (grepl("[NC]{1}-term", names(.x))) {
        site <- names(.x) %>% 
          gsub("(Protein|Any) ([NC]{1}-term)", "\\2", .)
      } else {
        site <- .x
      }
      
      m <- aa_masses[site]
      aa_masses[site] <<- m + .y
    })
  } else {
    #  Add mod_masses of variable mods (multiple lists)
    aas <- purrr::map2(positions_sites, mod_masses, ~ {
      if (grepl("[NC]{1}-term", names(.x))) {
        site <- names(.x) %>% 
          gsub("(Protein|Any) ([NC]{1}-term)", "\\2", .)
      } else {
        site <- .x
      }
      
      aa_masses[site] <- aa_masses[site] + .y
      
      aa_masses
    }, aa_masses)
    
    # Flatten the lists
    aa_masses <- local({
      sites <- purrr::map2_dbl(positions_sites, aas, ~ .y[.x])
      attrs <- attributes(aa_masses)
      aa_masses <- c(aa_masses, sites)
      attrs$names <- names(aa_masses)
      attributes(aa_masses) <- attrs
      
      aa_masses
    })
    
    rm(aas)
  }
  
  attr(aa_masses, mod_type) <- all_mods
  attr(aa_masses, paste0(mod_type, "_ps")) <- positions_sites
  attr(aa_masses, paste0(mod_type, "_mass")) <- mod_masses

  # add mod_masses of neutral losses
  no_nls <- neulosses %>% 
    purrr::map_lgl(~ all(.x == 0)) %>% 
    all()
  
  if (no_nls) {
    return(list(aa_masses))
  }
  
  ## Full nams for NLs
  #                          nl_1      nl_2       nl_3
  # Carbamidomethyl (M) 105.024835  0.000000 105.024835
  # Oxidation (M)         0.000000 63.998285  63.998285
  # Methyl (. = N)        0.000000  0.000000   0.000000
  # Deamidated (. = N)    0.000000  0.000000   0.000000
  
  nl_combi <- neulosses %>% 
    expand.grid() %>% 
    t() %>% 
    data.frame() %>% 
    dplyr::select(which(not_all_zero(.))) %>% 
    `colnames<-`(paste0("nl_", 1:ncol(.)))
  
  # Subtract NL masses
  # (NLs are realized and reflected in the current aa_masses)
  aa_masses_nl <- purrr::imap(nl_combi, ~ {
    sites <- purrr::map_chr(positions_sites, ~ {
      if (grepl("[NC]{1}-term", names(.x))) {
        site <- names(.x) %>% 
          gsub("(Protein|Any) ([NC]{1}-term)", "\\2", .)
      } else {
        site <- .x
      }
    }) 
    
    nms <- purrr::map2_chr(sites, names(sites), ~ {
      if (.y %in% names(aa_masses)) {
        .y
      } else {
        .x
      }
    })
    
    m <- aa_masses[nms]
    aa_masses[names(m)] <- m - .x
    
    attr(aa_masses, paste0(mod_type, "_neuloss")) <- .y
    attr(aa_masses, paste0(mod_type, "_nlmass")) <- .x
    
    aa_masses
  }, aa_masses, positions_sites) 
  
  aa_masses <- list(aa_masses)
  
  c(aa_masses, aa_masses_nl)
}


#' Calculates molecular weight of a polypeptide ([MH]+).
#'
#' @param fixedmods A character vector of fixed modifications. See also
#'   \link{parse_unimod} for grammars.
#' @param varmods A character vector of variable modifications.
#' @param maxn_vmods_setscombi Integer; the maximum number of combinatorial variable
#'   modifications and neutral losses.
#' @examples
#' \donttest{
#' library(purrr)
#'
#' x <- calc_aamasses()
#' x_att <- map(x, attributes)
#' names(x_att[[1]])
#' x_vmods <- map(x_att, `[`, c("vmods"))
#'
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (C)"), c("Acetyl (N-term)", 
#'                      "Gln->pyro-Glu (N-term = Q)", "Oxidation (M)"))
#'
#' x <- calc_aamasses(fixedmods = NULL)
#' x <- calc_aamasses(fixedmods = NULL, varmods = NULL)
#'
#' # Fixed mod, no NL
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                                  "Carbamidomethyl (. = C)"), varmods = NULL)
#'
#' # Fixed mod + NL
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                                  "Carbamidomethyl (. = M)"), varmods = NULL)
#'
#' # Fixed mod, no NL; var mod, no NL
#' x <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                                  "Carbamidomethyl (. = C)"), 
#'                    varmods = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)"))
#'
#' # Fixed mod + NL; var mod + NL
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)", 
#'                      "Carbamidomethyl (. = M)", 
#'                      "Deamidated (. = R)"), 
#'                    c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)",
#'                      "Hex(5)HexNAc(2) (N)"))
#' 
#' x <- calc_aamasses(c("TMT6plex (N-term)", "TMT6plex (K)", 
#'        "Carbamidomethyl (. = M)", "Deamidated (. = R)"), 
#'        c("Acetyl (N-term)", "Carbamyl (. = M)", 
#'        "Gln->pyro-Glu (N-term = Q)", "Hex(5)HexNAc(2) (N)"))
#' }
#' \dontrun{
#' # conflicts
#' x <- calc_aamasses(c("Carbamidomethyl (N-term)", "TMT2plex (N-term)"), NULL)
#' }
#' @export
calc_aamasses <- function (fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                                        "Carbamidomethyl (. = C)"), 
                           varmods = c("Acetyl (Protein N-term)", 
                                       "Oxidation (M)", 
                                       "Deamidated (N)", 
                                       "Gln->pyro-Glu (N-term = Q)"), 
                           maxn_vmods_setscombi = 64) {

  # title (position = site); 
  # . stands for (a) anywhere in position or (b) any residue in site or both
  # Acetyl (Protein N-term) <-> Acetyl (Protein N-term = .)
  # Acetyl (N-term) <-> Acetyl (N-term = .)
  # Carbamidomethyl (C) <-> Carbamidomethyl (. = C)
  # Carbamidomethyl <-> Carbamidomethyl(. = .)
  # Gln->pyro-Glu (N-term Q) <-> Gln->pyro-Glu (N-term = Q)
  # anywhere, can be skipped: title (site)
  # N-term, C-term, Protein N-term, Protein C-term
  
  ## The same site but different mods
  #  (a) Among fixedmods; Failure (only 'one_of', prompt to devise of a joint Unimod)
  #      (a1) TMT6plex (N-term)
  #      (a2) Biotin (N-term)
  #  (b) Between fixedmods and varmods; Relaxation (from fixedmods to varmods)
  #      (b1) Fixed TMT (N-term) -> Variable TMT (N-term) 
  #           Variable Acetyl (N-term)
  #      (b2) Fixed Oxidation (M) -> Variable Oxidation (M)
  #           variable Met->Ala (M)
  #  (c) Among varmods; OK if is 'one_of' in combination but not 'multiple'  
  #      (c1) Oxidation (M)
  #      (c2) Met->Ala (M)
  #      excludes (c1) + (c2)
  #      OK if different mods to the same residue at different indexes of a peptide
  #      (c3) dHex(1)Hex(1) (S): VS(3)SALSPSK
  #      (c4) Phospho (S): VSS(4)ALSPSK
  
  options(digits = 9)
  
  local({
    ## (0) Duplicated mods
    dup_mods <- intersect(fixedmods, varmods)
    
    if (!purrr::is_empty(dup_mods)) {
      stop("Modifications cannot be both 'fixed' and 'variable' at the same time: \n", 
           purrr::reduce(dup_mods, paste, sep = ", "), 
           call. = FALSE)
    }
  })
  
  new_mods <- local({
    # (a) Different fixedmods to the same site not allowed
    fmods_ps <- fixedmods %>% 
      purrr::map(find_unimod) %>% 
      purrr::map(`[[`, "position_site") %>% 
      `names<-`(fixedmods) %>% 
      purrr::flatten() 
    
    dup_fixedmods <- fmods_ps %>% 
      .[duplicated(.)]
    
    if (!purrr::is_empty(dup_fixedmods)) {
      stop("Multiple fixed modifications to the same site: \n", 
           "'", purrr::reduce(dup_fixedmods, paste, sep = ", "), "'",
           call. = FALSE)
    }
    
    # (b) Coercion from fixedmods to varmods
    vmods_ps <- varmods %>% 
      purrr::map(find_unimod) %>% 
      purrr::map(`[[`, "position_site") %>% 
      `names<-`(varmods) %>% 
      purrr::flatten() 
    
    dup_mods <- intersect(unlist(fmods_ps), unlist(vmods_ps)) 
    
    if (!purrr::is_empty(dup_mods)) {
      f_to_v <- local({
        idxes <- dup_mods %>% purrr::map(~ fmods_ps == .x)
        idxes %>% purrr::map_chr(~ fixedmods[.x]) 
      })
      
      varmods <- c(varmods, f_to_v)
      
      fixedmods <- local({
        idxes <- fixedmods %>% map_lgl(~ .x %in% f_to_v)
        fixedmods[!idxes]
      }) 
      
      warning("Coerce '", 
              purrr::reduce(f_to_v, paste, sep = ", "), "'", 
              " to variable modifications.", 
              call. = FALSE)
    }
    
    default_mods <- c("initiator methionine from protein N-terminus")
    if (any(grepl(default_mods, c(fixedmods, varmods)))) {
      warning("Modifications defaulted and no need to specify: `\n", 
              default_mods, "`.\n", 
              call. = FALSE)
      
    }
    
    invisible(list(fixedmods = fixedmods, varmods = varmods))
  })
  
  fixedmods <- new_mods$fixedmods
  varmods <- new_mods$varmods
  rm(new_mods)

  aa_masses <- c(
    A = 71.037114, R = 156.101111, N = 114.042927, D = 115.026943, 
    C = 103.009185, E = 129.042593, Q = 128.058578, G = 57.021464, 
    H = 137.058912, I = 113.084064, L = 113.084064, K = 128.094963,
    M = 131.040485, F = 147.068414, P = 97.052764, S = 87.032028,
    T = 101.047679, W = 186.079313, Y = 163.063329, V = 99.068414, 
    "N-term" = 1.007825, "C-term" = 17.002740, 
    U = 150.953633, B = 114.534940, X = 111.000000, Z = 128.550590, 
    "-" = 0)

  ## (1) add fixed mods + NL
  aa_masses_fi2 <- add_fixvar_masses(fixedmods, "fmods", aa_masses) 

  aa_masses_fi2 <- aa_masses_fi2 %>% 
    purrr::map(~ {
      if (is.null(attr(.x, "fmods"))) {
        attr(.x, "fmods") <- ""
      }
      
      if (is.null(attr(.x, "fmods_neuloss"))) {
        attr(.x, "fmods_neuloss") <- ""
      }
      
      if (is.null(attr(.x, "fmods_mass"))) {
        attr(.x, "fmods_mass") <- 0
      }
      
      if (is.null(attr(.x, "fmods_nlmass"))) {
        attr(.x, "fmods_nlmass") <- 0
      }
      
      .x
    })
  
  ## (2) add variable mods + NL
  varmods_comb <- local({
    # (c) Remove entries with multiple terminal mods
    varmods_comb <- seq_along(varmods) %>% 
      purrr::map(~ combn(varmods, .x, simplify = FALSE)) %>% 
      purrr::flatten()
    
    # Check if multiple terminal varmods by names
    #   Gln->pyro Glu (N-term = Q) and Acetyl (Protein N-term = N-term)
    #   have different sites but 'N-term' in names; and also need to be excluded
    vmods_ps <- varmods %>% 
      purrr::map(find_unimod) %>% 
      purrr::map(`[[`, "position_site") %>% 
      `names<-`(varmods) %>% 
      purrr::flatten() 
    
    vmods_ps_combi <- seq_along(vmods_ps) %>% 
      purrr::map(~ combn(vmods_ps, .x, simplify = FALSE)) %>% 
      purrr::flatten()
    
    dup_terms <- 
      vmods_ps_combi %>% purrr::map_lgl(~ sum(grepl("N-term", names(.x))) >= 2L) | 
      vmods_ps_combi %>% purrr::map_lgl(~ sum(grepl("C-term", names(.x))) >= 2L)

    # Not currently used
    dup_anywhere <- vmods_ps_combi %>% 
      purrr::map_lgl(~ {
        .x <- .x %>% 
          .[!grepl("[NC]{1}-term", names(.))]
        
        .x %>% 
          duplicated() %>% 
          any()
      })

    # Allow entries with different Anywhere mods to the same site
    #   (1) dHex(1)Hex(1) (S) and (2) Phospho (S)
    #   VS(1)S(2)ALSPSK
    
    varmods_comb <- varmods_comb  %>% 
      .[!(dup_terms)] %>% 
      purrr::map(unlist)
  })
  
  aa_masses_var2 <- purrr::map(varmods_comb, ~ {
    varmods_i <- .x
    
    aa_masses_fi2 %>% 
      purrr::map(~ add_fixvar_masses(varmods_i, "vmods", .x)) %>% 
      purrr::flatten()
  }, aa_masses_fi2) %>% 
    purrr::flatten()

  aa_masses_var2 <- aa_masses_var2 %>% 
    purrr::map(~ {
      if (is.null(attr(.x, "vmods"))) {
        attr(.x, "vmods") <- ""
      }
      
      if (is.null(attr(.x, "vmods_neuloss"))) {
        attr(.x, "vmods_neuloss") <- ""
      }
      
      if (is.null(attr(.x, "vmods_mass"))) {
        attr(.x, "vmods_mass") <- 0
      }
      
      if (is.null(attr(.x, "vmods_nlmass"))) {
        attr(.x, "vmods_nlmass") <- 0
      }
      
      .x
    })
  
  ## (3) complete 'vmods', 'vmods_neuloss' and 'vmods_ps' to fixedmods
  aa_masses_fi2 <- aa_masses_fi2 %>% 
    purrr::map(~ {
      if (is.null(attr(.x, "vmods"))) {
        attr(.x, "vmods") <- ""
      }
      
      if (is.null(attr(.x, "vmods_neuloss"))) {
        attr(.x, "vmods_neuloss") <- ""
      }
      
      if (is.null(attr(.x, "vmods_ps"))) {
        attr(.x, "vmods_ps") <- ""
      }
      
      if (is.null(attr(.x, "vmods_mass"))) {
        attr(.x, "vmods_mass") <- 0
      }
      
      if (is.null(attr(.x, "vmods_nlmass"))) {
        attr(.x, "vmods_nlmass") <- 0
      }
      
      .x
    })
  
  if (length(aa_masses_var2) >= maxn_vmods_setscombi) {
    warning("The number of ways of fixed and variable modifications is greater than ", 
            maxn_vmods_setscombi, ".\n", 
            "Some combinations may be dropped.", 
            call. = FALSE)
    aa_masses_var2 <- aa_masses_var2[1:maxn_vmods_setscombi]
  }
  
  c(aa_masses_fi2, aa_masses_var2)
}


#' Concatenates adjacent peptides in a list.
#'
#' Concatenates adjacent peptides, for example, to form peptides with
#' mis-cleavages.
#'
#' @param peps A list of peptide sequences with a one-letter representation of
#'   amino acid residues.
#' @param n The number of mis-cleavages for consideration.
#' @param include_cts Logical; the list, \code{peps}, includes the protein
#'   C-terminal sequence or not. At the default of TRUE, mis-cleaved peptides at
#'   the end of the protein C-terms will be added as they should. The arguments
#'   would be typically at FALSE, for example, when used for generating
#'   mis-cleaved peptides from the N-terminal of peptides with the removal of a
#'   starting residue \code{M}.
#'
#' @examples
#' \donttest{
#' peps <- c("MESNHK", "SGDGLSGTQK", "EAALR", "ALVQR", "TGYSLVQENGQR")
#' res <- concat_peps(peps, 2)
#' 
#' peps2 <- gsub("^M", "", peps[1:(2+1)])
#' res2 <- concat_peps(peps2, 2, FALSE)
#' 
#' c(res2, res)
#' 
#' # a short protein
#' concat_peps(c(TRDD1_HUMAN = "EI"))
#' }
concat_peps <- function (peps, n = 2, include_cts = TRUE) {
  len <- length(peps)
  
  if (n >= len) n <- len - 1

  res <- purrr::map(seq_len((len - n)), ~ {
    peps[.x:(.x + n)] %>% purrr::accumulate(paste0) 
  }) %>% unlist()
  
  if (include_cts && n >= 1) {
    res_cts <- local({
      cts <- peps[(len - n + 1):len]
      
      purrr::map(n:1, ~ {
        tail(cts, .x) %>% purrr::accumulate(paste0) 
      }) %>% unlist()
    })
  } else {
    res_cts <- NULL
  }

  # res_cts_a <- res_cts %>% .[!grepl("-$", .)]
  # res_cts_b <- res_cts %>% .[grepl("-$", .)]
  
  # c(res[-length(res)], res_cts_a, res[length(res)], res_cts_b)
  c(res, res_cts)
}


#' Parse the name of a Unimod
#' 
#' The general format: \code{parse_unimod("title (position = site)")}.
#' 
#' @param unimod The name of a \href{https://www.unimod.org/}{Unimod} modification. 
#' @examples
#' \donttest{
#' # "dot" for anywhere (either position or site)
#' x1 <- parse_unimod("Carbamidomethyl (. = C)")
#' x2 <- parse_unimod("Carbamidomethyl (Anywhere = C)")
#' x3 <- parse_unimod("Carbamidomethyl (C)")
#' 
#' identical(x1, x2); identical(x2, x3)
#' 
#' # Any residue on protein N-term 
#' x1 <- parse_unimod("Acetyl (Protein N-term = .)")
#' x2 <- parse_unimod("Acetyl (Protein N-term)")
#' 
#' identical(x1, x2)
#' 
#' # Any N-term residue
#' x1 <- parse_unimod("Acetyl (N-term)")
#' x2 <- parse_unimod("Acetyl (N-term = .)")
#' 
#' identical(x1, x2)
#' 
#' # N-term Q
#' x1 <- parse_unimod("Gln->pryo-Glu (N-term = Q)")
#' x2 <- parse_unimod("Gln->pryo-Glu (N-term Q)")
#' 
#' identical(x1, x2)
#' 
#' # ok with parenthesis in the 'title' 
#' x <- parse_unimod("Hex(5)HexNAc(2) (N)")
#' }
#' 
#' \dontrun{
#' # No modification of anywhere and anything 
#' # (to every position and site)
#' x <- parse_unimod("Carbamidomethyl (Anywhere = .)")
#' x <- parse_unimod("Carbamidomethyl")
#' x <- parse_unimod("Carbamidomethyl (. = .)")
#' 
#' # Prefer an "=" sign between 'N-term' and 'Q'
#' x <- parse_unimod("Gln->pyro-Glu (N-term Q)")
#' }
#' @export
parse_unimod <- function (unimod) {
  # unimod = "Carbamidomethyl (Protein N-term = C)" # --> pos_site = "Protein N-term = C"
  # unimod = "Carbamidomethyl (Any N-term = C)" # --> pos_site = "Any N-term = C"
  # unimod = "Carbamidomethyl (N-term = C)" # --> pos_site = "N-term = C"
  # unimod = "Carbamidomethyl (. = C)" # --> pos_site = ". = C"
  # unimod = "Carbamidomethyl (C)" # --> pos_site = "C"
  # unimod = "Carbamidomethyl ()" # --> pos_site = ""
  # unimod = "Carbamidomethyl" # --> pos_site = ""
  # unimod = "" # --> pos_site = ""
  
  ## any N-term residue
  # unimod = "Carbamidomethyl (Protein N-term = .)" 
  
  # unimod = "Hex(5)HexNAc(2) (N)"
  
  ## dual parentheses
  # unimod = "Carbamidomethyl ((. = C))" # --> pos_site = ". = C"
  
  if (grepl("([NC]{1}-term|Anywhere) [A-Z]{1}", unimod)) {
    unimod <- unimod %>% 
      gsub("^(.*[NC]{1}-term|.*Anywhere)\\s*([A-Z]{1})", "\\1 = \\2", .)
  }

  title <- unimod %>% 
    gsub("^([^ ]+?) .*", "\\1", .)

  pos_site <- unimod %>% 
    gsub("^[^ ]+", "", .) %>% 
    gsub("^[^\\(]+[\\(]*([^\\)]*)[\\)]*$", "\\1", .)
  
  if (grepl("=", pos_site)) {
    pos <- pos_site %>% 
      gsub("^([^=]+?)[=].*", "\\1", .) %>% 
      gsub("^[ ]*", "\\1", .) %>% 
      gsub(" *$", "", .)
    
    site <- pos_site %>% 
      gsub("^[^=]+?[=](.*)", "\\1", .) %>% 
      gsub("^[ ]*", "\\1", .) %>% 
      gsub(" *$", "", .)
  } else {
    pos <- "."
    site <- pos_site
  }
  
  if (site == "") {
    site = "."
  }
  
  if (site %in% c("Protein N-term", "Protein C-term", 
                  "Anywhere N-term", "Anywhere C-term", 
                  "N-term", "C-term")) {
    pos <- site
    site <- site %>% gsub("^(Protein|Anywhere) ", "", .)
  } 
  
  # standardize `position`
  pos <- pos %>% 
    gsub("^([NC]){1}-term", "Any \\1-term", .)
  
  if (pos %in% c(".", "")) {
    pos <- "Anywhere"
  }
  
  pos_allowed <- c("Anywhere", "Protein N-term", "Protein C-term", 
                   "Any N-term", "Any C-term")

  stopifnot(pos %in% pos_allowed)
  
  # standardize terminal sites 
  if (site == ".") {
    if (pos %in% c("Protein N-term", "Any N-term")) {
      site <- "N-term"
    } else if (pos %in% c("Protein C-term", "Any C-term")) {
      site <- "C-term"
    }
  }
  
  if (pos == "Anywhere" && site == ".") {
    stop("'position' or 'site' cannot be both 'Anywhere'.", 
         call. = FALSE)
  }

  invisible(list(title = title, position = pos, site = site))
}


#' Find a Unimod.
#'
#' Find the mono-isotopic mass, position, site and neutral losses of a
#' modification. 
#' 
#' In the field of \code{position_site}, \code{position} is the name and
#' \code{site} is the value.
#'
#' @inheritParams parse_unimod
#' @examples
#' \donttest{
#' x <- find_unimod("Carbamidomethyl (C)")
#' x <- find_unimod("Carbamidomethyl (M)")
#' x <- find_unimod("Acetyl (Protein N-term)")
#' x <- find_unimod("Gln->pyro-Glu (N-term = Q)")
#' x <- find_unimod("Hex(5)HexNAc(2) (N)")
#' }
#'
#' \dontrun{
#' # Prefer an "=" sign between 'N-term' and 'Q'
#' x <- find_unimod("Gln->pyro-Glu (N-term Q)")
#' }
#' @export
find_unimod <- function (unimod = "Carbamidomethyl (C)") {
  options(digits = 9)
  
  res <- parse_unimod(unimod)
  title <- res$title
  position <- res$position
  site <- res$site
  rm(res)

  parent <- system.file("extdata", "master.xml", package = "proteoQ") %>% 
    xml2::read_xml()
  
  # <umod:elements>
  # <umod:modifications>
  # <umod:amino_acids>
  # <umod:mod_bricks>
  
  children <- xml2::xml_children(parent)
  contents <- parent %>% xml2::xml_contents() 
  
  elements <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "elements")]])
  modifications <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "modifications")]])
  amino_acids <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "amino_acids")]])
  mod_bricks <- 
    xml2::xml_children(children[[which(xml2::xml_name(contents) == "mod_bricks")]])
  
  this_mod <- local({
    idx <- which(xml2::xml_attr(modifications, "title") == title)
    
    if (purrr::is_empty(idx)) {
      stop("Modification not found: '", title, "'.\n", 
           "For example, use 'Acetyl' instead of 'Acetylation'.", 
           call. = FALSE)
    }
    
    this_mod <- modifications[[idx]]
  })
  
  stopifnot(xml2::xml_attrs(this_mod) %>% .["title"] == title)
  
  # --- children of `this_mod` ---
  modch <- xml2::xml_children(this_mod)
  
  # --- find mass --- 
  monomass <- xml2::xml_attrs(modch) %>% 
    purrr::map(`[`, "mono_mass") %>%
    `[`(!is.na(.)) %>% 
    unlist() %>% 
    as.numeric()
  
  stopifnot(length(monomass) == 1)
  
  # --- find sites and positions ---
  positions_sites <- xml2::xml_attrs(modch) %>% 
    purrr::map(`[`, c("site", "position")) %>% 
    purrr::map(~ invisible(setNames(.x[1], .x[2]))) 
  
  sites <- xml2::xml_attrs(modch) %>% 
    map(`[`, c("site")) %>% 
    unlist() 
    
  positions <- xml2::xml_attrs(modch) %>% 
    map(`[`, c("position")) %>% 
    unlist() 

  # --- neutral loss ---
  idx_nl <- grep("NeutralLoss", modch)
  
  if (!purrr::is_empty(idx_nl)) {
    nls <- purrr::map(idx_nl, ~ {
      modnl <- xml2::xml_children(modch[.x]) 
      
      xml2::xml_attrs(modnl) %>% 
        purrr::map(`[`, "mono_mass") %>%
        `[`(!is.na(.)) %>% 
        as.numeric() %>% 
        setNames(paste("nl", 1:length(.), sep = "."))
    }) %>% 
      setNames(rep("nl", length(.)))
    
    positions_sites[idx_nl] <- 
      purrr::map2(idx_nl, nls, ~ {c(positions_sites[[.x]], .y)}) 
  } 
  
  positions_sites <- positions_sites[sites == site & positions == position] %>% 
    .[purrr::map_lgl(., ~ !is.null(.x))]
  
  if (purrr::is_empty(positions_sites)) {
    stop("'", unimod, "' not found.", call. = FALSE)
  }

  neulosses <- positions_sites[[1]][-1]
  if (purrr::is_empty(neulosses)) {
    neulosses <- 0
  } else {
    neulosses <- as.numeric(neulosses)
  }
  
  invisible(list(monomass = monomass, 
                 position_site = positions_sites[[1]][1], 
                 nl = neulosses))
}


#' Find mis-cleavages in a vector.
#'
#' A convenience utility may be used to extract the first \eqn{n+1} peptides
#' from 0 to n mis-cleavages. It also assumes that the data were already sorted
#' in a desirable way.
#'
#' @param x A vector of data.
#' @param n Integer. The number of mis-cleavages.
keep_n_misses <- function (x, n) {
  len <- length(x)
  
  stopifnot(n >= 0L, len >= 1L)
  
  x[1:pmin(n + 1, len)]
} 


#' Exclude mis-cleavages in a vector.
#'
#' @inheritParams keep_n_misses
#' @seealso keep_n_misses
exclude_n_misses <- function (x, n) {
  len <- length(x)
  
  stopifnot(n >= 0L, len >= 1L)

  x[-(1:pmin(n + 1, len))]
}


#' Excludes a character in string counting
#' 
#' @param x A character string
#' @param char A character to be excluded for counting.
#' @importFrom stringi stri_length stri_count_fixed
str_exclude_count <- function (x, char = "-") {
  stringi::stri_length(x) - stringi::stri_count_fixed(x, char)
}


#' Remove a starting character from the first \code{n} entries.
#' 
#' @param x A list of character strings.
#' @param char A starting character to be removed.
#' @param n The number of beginning entries to be considered.
rm_char_in_nfirst <- function (x, char = "^-", n = (max_miss + 1) * 2) {
  n <- pmin(length(x), n)
  x[seq_len(n)] <- gsub(char, "", x[seq_len(n)])
  
  x
}


#' Remove a trailing character from the last \code{n} entries.
#' 
#' @param char A trailing character to be removed.
#' @inheritParams rm_char_in_nfirst
rm_char_in_nlast <- function (x, char = "-$", n = (max_miss + 1) * 2) {
  len <- length(x)
  n <- pmin(len, n)
  x[(len-n+1):len] <- gsub(char, "", x[(len-n+1):len])
  
  x
}


#' Calculates the mono-isotopic mass of a peptide sequence.
#' 
#' For direct uses from an R console (with trade-offs in speed).
#'
#' @inheritParams calc_monopep
#'
#' @examples
#' \dontrun{
#' ## No modifications
#' data(package = "proteoQ", aa_residues)
#' aa_masses <- aa_residues %>%
#'   dplyr::select(c("one_letter", "monoisotopic_da"))
#' aa_masses <- aa_masses$monoisotopic_da %>% `names<-`(aa_masses$one_letter)
#' aa_masses["-"] <- 0
#'
#' x <- calc_monopeptide("AAIDWFDGKEFSGNPIK", aa_masses)
#' x <- calc_monopeptide("-AAAAAAAGDSDSWDADAFSVEDPVR", aa_masses)
#'
#' ## With possible modifications
#' library(purrr)
#'
#' aa_masses_all <- calc_aamasses()
#'
#' calc_monopeptide("GFGFVSNFER", aa_masses_all[[1]])
#' calc_monopeptide("GFGMFVSNFER", aa_masses_all[[1]])
#' calc_monopeptide("GFGMFVSNMFER", aa_masses_all[[1]])
#' calc_monopeptide("GFGMFVSNMFER", aa_masses_all[[11]])
#'
#' # Error if modifications not in 'aa_masses_all'
#' # (`subpeps_by_vmods()` subsets peptides by 
#' #  variable modifications in an 'aa_masses')
#' map(aa_masses_all, ~ calc_monopeptide("GFGFVTFSNSMAEVDAAMAAR", .x))
#'
#' calc_monopeptide("GFGFNNTFSSMAEVDAMMANAR", aa_masses_all[[11]])
#' map(aa_masses_all, ~ calc_monopeptide("GFGFNNTFSSMAEVDAMMANAR", .x))
#' }
#' @export
calc_monopeptide <- function (aa_seq, aa_masses, 
                              mod_indexes = NULL, maxn_vmods_per_pep = 5, 
                              maxn_sites_per_vmod = 3, digits = 5) {
  
  options(digits = 9)
  
  vmods_ps <- aa_masses %>% 
    attributes() %>% 
    `[[`("vmods_ps")
  
  # multiple mods to [NC]-term already excluded from aa_masses
  amods <- local({
    sites <- vmods_ps %>% 
      purrr::map(~ .x[grepl("Anywhere", names(.x))]) 
    
    empties <- sites %>% purrr::map_lgl(purrr::is_empty)
    
    sites <- sites[!empties] 
    # names(sites) <- mod_indexes[names(sites)]
    sites
  })
  
  tmod <- vmods_ps %>% .[! . %in% amods]
  if (purrr::is_empty(tmod)) {
    tmod <- NULL
  } else if (tmod == "") {
    tmod <- NULL
  }
  
  ntmod <- tmod %>% .[. == "N-term"]
  ctmod <- tmod %>% .[. == "C-term"]
  
  calc_monopep(aa_seq = aa_seq, 
               aa_masses = aa_masses, 
               vmods_ps = vmods_ps, 
               amods = amods, 
               tmod = tmod, 
               ntmod = ntmod, 
               ctmod = ctmod, 
               mod_indexes = mod_indexes, 
               maxn_vmods_per_pep = maxn_vmods_per_pep, 
               maxn_sites_per_vmod = maxn_sites_per_vmod, 
               digits = digits)
}


#' Generates and Calculates the masses of tryptic peptides from a fasta
#' database.
#'
#' @param enzyme A character string; the proteolytic specificity of the assumed
#'   enzyme will be used to generate peptide sequences from proteins. The enzyme
#'   is currently \code{trypsin}.
#' @param maxn_fasta_seqs Integer; the maximum number of protein sequences in
#'   fasta files.
#' @param min_len Integer; the minimum length of peptides. Shorter peptides will
#'   be excluded.
#' @param max_len Integer; the maximum length of peptides. Longer peptides will
#'   be excluded.
#' @param max_miss The maximum number of mis-cleavages per peptide sequence.
#' @param out_path The file name with prepending path for outputs.
#' @param digits Integer; the number of decimal places to be used.
#' @inheritParams splitPSM
#' @inheritParams calc_aamasses
#' @inheritParams mcalc_monopep
#' @inheritParams calc_monopep
#' @examples
#' \donttest{
#' res <- calc_pepmasses()
#'
#' library(purrr)
#' library(magrittr)
#'
#' res_attrs <- map(res, attributes)
#' map(res_attrs, names)
#' res_attrs %>% map(`[[`, "vmods")
#' res_mods <- map(res_attrs, `[`,
#'                 c("fmods", "fmods_ps", "fmods_neuloss",
#'                   "vmods", "vmods_ps", "vmods_neuloss"))
#'
#' res_data <- map(res_attrs, `[[`, "data")
#' peps_combi_1 <- res_data[[1]]
#'
#' # base: fixedmods without neulosses
#' map(res_data[[1]], length) %>% reduce(sum)
#'
#' # fixedmods, fixedmods + fixedmods_neulosses, varmods, varmods_neulosses
#' map(res_data, ~ map(.x, length) %>% reduce(sum)) %>% reduce(sum)
#'
#' }
#'
#' @export
calc_pepmasses <- function (fasta = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
                            fixedmods = c("TMT6plex (K)", 
                                          "Carbamidomethyl (C)"), 
                            varmods = c("TMT6plex (N-term)", 
                                        "Acetyl (Protein N-term)", 
                                        "Oxidation (M)", 
                                        "Deamidated (N)", 
                                        "Gln->pyro-Glu (N-term = Q)"), 
                            enzyme = c("trypsin"), 
                            maxn_fasta_seqs = 50000,
                            maxn_vmods_setscombi = 64, 
                            maxn_vmods_per_pep = 5,
                            maxn_sites_per_vmod = 3, 
                            min_len = 7, max_len = 100, max_miss = 2, 
                            digits = 5, 
                            out_path = "~/proteoQ/dbs/fasta/uniprot/pepmass/pepmasses.rds") {

  stopifnot(vapply(c(min_len, max_len, max_miss), is.numeric, 
                   logical(1)))
  
  stopifnot(min_len >= 0, max_len >= min_len, max_miss <= 100)
  
  mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
    as.hexmode() %>% 
    `names<-`(c(fixedmods, varmods))
  
  aa_masses <- calc_aamasses(fixedmods, varmods, maxn_vmods_setscombi)
  
  fasta_db <- fasta %>% load_fasta() 

  if (length(fasta_db) > maxn_fasta_seqs) {
    stop("More than `", maxn_fasta_seqs, "` sequences in fasta files.", 
         call. = FALSE)
  }
  
  inds_m <- grep("^M", fasta_db)
  
  fasta_db <- fasta_db %>% 
    purrr::map(~ gsub("([KR]{1})", paste0("\\1", "@"), .x) %>% 
                 paste0("-", ., "-")) 
  
  run_scripts <- FALSE
  if (run_scripts) {
    peptide <- paste0("-WASQVSENR@PVCK@AIIQGK@QFEGLVDTGADVSIIALNQWPK@NWPK@QK@", 
                      "AVTGLVGIGTASEVYQSTEILHCLGPDNQESTVQPMITSIPLNLWGR@", 
                      "DLLQQWGAEITMPAPLYSPTSQK@IMTK@MGYIPGK@GLGK@NEDGIK@", 
                      "IPFEAK@INQK@R@EGIGYPF-")

    gsub(paste0("^(", 
                rep("[^@]*?@{1}", max_miss+1) %>% purrr::reduce(paste0), ")", 
                "([^@]*?@{1}).*$"), 
         "\\1", peptide)
  }
  
  fasta_dbm <- fasta_db[inds_m] %>% 
    purrr::map(~ gsub("^-M", "-", .x))

  # --- Protein N-term initiator methionine kept ---
  peps <- fasta_db %>% 
    purrr::map(~ .x %>% stringr::str_split("@", simplify = TRUE)) %>% 
    purrr::map(concat_peps, max_miss) %>% 
    purrr::map(
      ~ .x %>% .[str_exclude_count(.) >= min_len & str_exclude_count(.) <= max_len]) 
  
  # --- Protein N-term initiator methionine removed ---
  peps_m <- fasta_dbm %>% 
    purrr::map(~ .x %>% 
                 stringr::str_split("@", simplify = TRUE) %>% 
                 keep_n_misses(max_miss)) %>% 
    purrr::map(~ concat_peps(.x, max_miss, include_cts = FALSE)) %>% 
    purrr::map(
      ~ .x %>% .[str_exclude_count(.) >= min_len & str_exclude_count(.) <= max_len]) 
  
  peps[inds_m] <- purrr::map2(peps_m, peps[inds_m], `c`)
  
  pep_masses <- aa_masses %>% 
    purrr::map(subpeps_by_vmods, peps) %>% 
    purrr::map(~ {
      xs <- .x
      
      xs %>% 
        purrr::map(rm_char_in_nfirst, char = "^-", n = (max_miss + 1) * 2) %>% 
        purrr::map(rm_char_in_nlast, char = "-$", n = (max_miss + 1) * 2)
    })
  
  pep_masses <- pep_masses %>% 
    purrr::map2(aa_masses, mcalc_monopep, mod_indexes, 
                maxn_vmods_per_pep, 
                maxn_sites_per_vmod, 
                digits)
  
  out <- purrr::map2(aa_masses, pep_masses, ~ {
    attr(.x, "data") <- .y
    return(.x)
  })

  dir.create(gsub("(^.*/).*$", "\\1", out_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, file.path(out_path))
  
  invisible(out)
}

