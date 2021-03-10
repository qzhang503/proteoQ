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


#' Calculates multiply the molecular weight of a polypeptides ([MH]+).
#'
#' The calculations iterate through protein accessions and peptide sequences.
#'
#' @param aa_seqs Character string; a vector of peptide sequences with
#'   one-letter representation of amino acids.
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
#' aa_masses_all <- ext_aa_masses()
#' x <- mcalc_monopep(list(protein_a = c("AAIMDWFDMGKEFSGNPIK", "MAAIDWFDGKEFSGNPIK"), 
#'                         protein_b = c("AAIDWFDGKEFSGNPIK")), 
#'                    aa_masses_all[[11]], aa_masses_all[[1]])
#' }
mcalc_monopep <- function (aa_seqs, aa_masses, aa_masses_ref, digits = 5) {
  attr(aa_masses, "data") <- aa_seqs %>% 
    purrr::map(~ { # by prot_acc
      .x %>% 
        purrr::map(calc_monopep, aa_masses, aa_masses_ref, digits) %>% # by peptides
        unlist() 
    })
}


#' Helper: calculates the mono-isotopic mass of a peptide sequence.
#'
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param aa_masses_ref A named list containing the reference (mono-isotopic)
#'   masses of amino acid residues.
#' @inheritParams mcalc_monopep
#'
#' @examples
#' \dontrun{
#' data(package = "proteoQ", aa_residues)
#' aa_masses <- aa_residues %>%
#'   dplyr::select(c("one_letter", "monoisotopic_da"))
#' aa_masses <- aa_masses$monoisotopic_da %>% `names<-`(aa_masses$one_letter)
#' aa_masses["-"] <- 0
#'
#' library(purrr)
#' 
#' x <- calc_monopep("AAIDWFDGKEFSGNPIK", aa_masses)
#' x <- calc_monopep("-AAAAAAAGDSDSWDADAFSVEDPVR", aa_masses)
#'
#'
#' aa_masses_all <- ext_aa_masses()
#' calc_monopep("GFGFVSFER", aa_masses_all[[1]])
#'
#' calc_monopep("GFGFVTFSSMAEVDAAMAAR", aa_masses_all[[11]], aa_masses_all[[1]])
#' map(aa_masses_all, ~ calc_monopep("GFGFVTFSSMAEVDAAMAAR", .x, aa_masses_all[[1]]))
#' 
#' calc_monopep("GFGFNNTFSSMAEVDAMMANAR", aa_masses_all[[11]], aa_masses_all[[1]])
#' map(aa_masses_all, ~ calc_monopep("GFGFNNTFSSMAEVDAMMANAR", .x, aa_masses_all[[1]]))
#' }
calc_monopep <- function (aa_seq, aa_masses, aa_masses_ref = aa_masses, digits = 5) {
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # oxygen <- 15.99491462
  
  # moverz <- (mass + h2o + z*proton)/z
  
  # data(package = "proteoQ", aa_residues)
  # aa_masses <- aa_residues %>% 
  #   dplyr::select(c("one_letter", "monoisotopic_da")) 
  # aa_masses <- aa_masses$monoisotopic_da %>% `names<-`(aa_masses$one_letter)
  
  options(digits=9)
  
  monomass <- aa_seq %>% 
    stringr::str_split("", simplify = TRUE) %>% 
    aa_masses[.] %>% 
    sum() %>% 
    `+`(aa_masses["N-term"]) %>% 
    `+`(aa_masses["C-term"]) %>% 
    setNames(aa_seq) %>% 
    round(digits = digits) 
  
  # ---
  vmods_ps <- aa_masses %>% 
    attributes() %>% 
    `[[`("vmods_ps")
  
  sites <- local({
    sites <- vmods_ps %>% 
      purrr::map(~ .x[grepl("Anywhere", names(.x))]) 
    
    empties <- sites %>% purrr::map_lgl(purrr::is_empty)
    
    sites[!empties]
  })
  
  if (purrr::is_empty(sites)) return(monomass)
  
  sites <- sites %>% purrr::flatten() %>% unlist()
  
  monomasses <- monomass
  for (i in seq_along(sites)) {
    site = sites[i]
    
    if (names(site) == "Anywhere") {
      mass_delta <- (aa_masses - aa_masses_ref) %>% 
        .[names(.) == site]
      
      monomasses <- sub_modmasses(aa_seq, monomasses, site, mass_delta, digits)
    } 
  }
  
  names(monomasses) <- rep(aa_seq, length(monomasses))
  
  monomasses
}


#' Calculates of the mono-isotopic masses for \code{Anywhere} variable
#' modifications.
#'
#' At all counts of variable modifications for a given \code{site}. For example,
#' calculations of masses at 1 to 3 oxidations of M. Only counts, not
#' permutations.
#'
#' @param monomasses A vector of mono-isotopic masses for peptides.
#' @param site An amino-acide site.
#' @param mass_delta The mass delta between current aa_masses table and the
#'   reference ass_masses table.
#' @inheritParams calc_monopep
sub_modmasses <- function (aa_seq, monomasses, site, mass_delta, digits = 5) {
  counts <- stringr::str_count(aa_seq, site)
  
  if (counts >= 2L) {
    out <- local({
      out <- vector("list", counts)
      out[[1]] <- monomasses
      
      for (i in 2:counts) {
        out[[i]] <- out[[i-1]] - mass_delta
      }
      
      out
    })
  } else {
    out <- monomasses
  }
  
  unlist(out, use.names = FALSE)
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
  options(digits=9)
  
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
  
  masses <- res %>% purrr::map(`[[`, 1)
  positions_sites <- res %>% purrr::map(`[[`, 2)
  neulosses <- res %>% purrr::map(`[[`, 3)
  rm(res)
  
  local({
    if (mod_type == "fmods" && length(positions_sites) > 1) {
      dups <- purrr::reduce(positions_sites, `c`) %>% 
        .[duplicated(.)]
      
      if (!purrr::is_empty(dups)) {
        dups_in_each <- positions_sites %>% 
          purrr::map(~ .x[.x == dups])
        
        stop("Conflicts in fixed modifications: \n", 
             purrr::reduce(names(dups_in_each), paste, sep = "\n"), "\n",
             "May consider change from fixed to variable modifications(s).", 
             call. = FALSE)
      }
    }
  })

  # add masses of (a) fixed or (b) variable mods
  purrr::walk2(positions_sites, masses, ~ {
    if (grepl("[NC]{1}-term", names(.x))) {
      site <- names(.x) %>% 
        gsub("(Protein|Any) ([NC]{1}-term)", "\\2", .)
    } else {
      site <- .x
    }
    
    m <- aa_masses %>% .[names(.) == site]
    aa_masses[names(aa_masses) == site] <<- m + .y
  })
  
  attr(aa_masses, mod_type) <- all_mods
  attr(aa_masses, paste0(mod_type, "_ps")) <- positions_sites

  # add masses of neutral losses
  no_nls <- neulosses %>% 
    purrr::map_lgl(~ all(.x == 0)) %>% 
    all()
  
  if (no_nls) {
    return(list(aa_masses))
  }
  
  nl_combi <- neulosses %>% 
    expand.grid() %>% 
    `colnames<-`(positions_sites) %>% 
    t() %>% 
    data.frame() %>% 
    dplyr::select(which(not_all_zero(.))) %>% 
    `colnames<-`(paste0("nl_", 1:ncol(.)))

  aa_masses_nl <- purrr::imap(nl_combi, ~ {
    m <- aa_masses %>% 
      .[names(.) %in% unlist(positions_sites)]
    
    aa_masses[names(m)] <- m - .x
    
    attr(aa_masses, paste0(mod_type, "_neuloss")) <- .y
    
    aa_masses
  }, aa_masses) 

  aa_masses <- list(aa_masses)
  
  c(aa_masses, aa_masses_nl)
}


#' Calculates molecular weight of a polypeptide ([MH]+).
#'
#' @param max_ncomb_mods Integer; the maximum number of combinatorial
#'   variable modifications allowed per peptide sequence.
#' @param digits Integer; the number of decimal places to be used.
#' @examples
#' \donttest{
#' library(purrr)
#' 
#' x <- ext_aa_masses()
#' x_att <- map(aa_masses_all, attributes)
#' names(x_att[[1]])
#' x_vmods <- map(x_att, `[`, c("vmods"))
#' 
#' x <- ext_aa_masses(c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)"), c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)", "Oxidation (M)"))
#' 
#' x <- ext_aa_masses(fixedmods = NULL)
#' x <- ext_aa_masses(fixedmods = NULL, varmods = NULL)
#'
#' # Fixed mod, no NL
#' x <- ext_aa_masses(fixedmods = c("TMT6plex (N-term)", "Carbamidomethyl (. = C)"), varmods = NULL)
#'
#' # Fixed mod + NL
#' x <- ext_aa_masses(fixedmods = c("TMT6plex (N-term)", "Carbamidomethyl (. = M)"), varmods = NULL)
#'
#' # Fixed mod, no NL; var mod, no NL
#' x <- ext_aa_masses(fixedmods = c("TMT6plex (N-term)", "Carbamidomethyl (. = C)"), varmods = c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)"))
#' 
#' # Fixed mod + NL; var mod + NL
#' x <- ext_aa_masses(c("TMT6plex (N-term)", "Carbamidomethyl (. = M)",
#' "Deamidated (. = R)"), c("Acetyl (N-term)", "Gln->pyro-Glu (N-term = Q)",
#' "Hex(5)HexNAc(2) (N)"))
#' }
#' \dontrun{
#' # conflicts
#' x <- ext_aa_masses(c("Carbamidomethyl (N-term)", "TMT2plex (N-term)"), NULL)
#' }
ext_aa_masses <- function (fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                                         "Carbamidomethyl (. = C)"), 
                           varmods = c("Acetyl (Protein N-term)", 
                                       "Oxidation (M)", 
                                       "Deamidated (N)", 
                                       "Gln->pyro-Glu (N-term = Q)"), 
                           max_ncomb_mods = 64) {

  # title (position = site); 
  # . stands for (a) anywhere in position or (b) any residue in site or both
  # Acetyl (Protein N-term) <-> Acetyl (Protein N-term = .)
  # Acetyl (N-term) <-> Acetyl (N-term = .)
  # Carbamidomethyl (C) <-> Carbamidomethyl (. = C)
  # Carbamidomethyl <-> Carbamidomethyl(. = .)
  # Gln->pyro-Glu (N-term Q) <-> Gln->pyro-Glu (N-term = Q)
  # anywhere, can be skipped: title (site)
  # N-term, C-term, Protein N-term, Protein C-term
  
  options(digits=9)
  
  local({
    dup_mods <- intersect(fixedmods, varmods)
    
    if (!purrr::is_empty(dup_mods)) {
      stop("Modifications cannot be both 'fixed' and 'variable' at the same time: \n", 
           purrr::reduce(dup_mods, paste, sep = ", "), 
           call. = FALSE)
    }
  })
  
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
      
      .x
    })
  
  ## (2) add variable mods + NL
  varmods_comb <- seq_len(pmin(max_ncomb_mods, length(varmods))) %>% 
    purrr::map(~ combn(varmods, .x, simplify = FALSE)) %>% 
    purrr::flatten()
  
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
      
      .x
    })
  
  # map(aa_masses_var2, ~ which(unlist(.x) != aa_masses_fi2[[1]]))
  
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
find_unimod <- function (unimod = "Carbamidomethyl (C)") {
  options(digits=9)
  
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


#' Find the first \eqn{n + 1} elements in a vector.
#'
#' A convenience utility may be used to extract the first \code{n} peptides with
#' mis-cleavages, plus the one without mis-clevages (a total of \code{n + 1}).
#' Assumed the data were already sorted in a desirable way.
#'
#' @param x A vector of data.
#' @param n Integer. The number of mis-cleavages.
keep_one_to_n <- function (x, n) {
  stopifnot(n >= 0)
  x[1:pmin(n + 1, length(x))]
} 


#' Exclude the first \eqn{n + 1} elements in a vector.
#' 
#' @inheritParams keep_one_to_n
exclude_one_to_n <- function (x, n) {
  stopifnot(n >= 0)
  x[-(1:pmin(n+1, length(x)))]
}


#' Generates and Calculates the masses of tryptic peptides from a fasta
#' database.
#'
#' @param fixedmods A character vector of fixed modifications. See also
#'   \link{parse_unimod} for grammars.
#' @param varmods A character vector of variable modifications.
#' @param min_len The minimum length of peptides. Shorter peptides will be
#'   excluded.
#' @param max_miss The maximum length of peptides. Longer peptides will be
#'   excluded.
#' @param digits Integer; the number of decimal places to be used.
#' @param out_path The file path for outputs.
#' @param out_name The file name for outputs.
#' @inheritParams splitPSM
#' @inheritParams ext_aa_masses
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
#' res_mods <- map(res_attrs, `[`, c("fmods", "fmods_ps", "fmods_neuloss", "vmods", "vmods_ps", "vmods_neuloss"))
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
                            fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", 
                                          "Carbamidomethyl (C)"), 
                            varmods = c(c("Acetyl (N-term)", "Oxidation (M)", 
                                          "Deamidated (N)")), 
                            max_ncomb_mods = 64, min_len = 7, max_len = 100, max_miss = 2, 
                            digits = 5, 
                            out_path = "~/proteoQ/dbs/fasta/uniprot/pepmasses/pepmasses.rds") {

  stopifnot(vapply(c(min_len, max_len, max_miss), is.numeric, 
                   logical(1)))
  
  stopifnot(min_len >= 0, max_len >= min_len, max_miss <= 100)
  
  fasta_db <- fasta %>% 
    load_fasta() 

  # fasta_db <- fasta_db[c(760:762, 950:954, 5020:5099)]
  # fasta_db <- fasta_db[which(names(fasta_db) == "NP_061185")]
  # fasta_db <- fasta_db[which(names(fasta_db) == "NP_000691")]
  # 
  
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

  # --- N-term methionine kept ---
  peps <- fasta_db %>% 
    purrr::map(~ .x %>% stringr::str_split("@", simplify = TRUE)) %>% 
    purrr::map(concat_peps, max_miss) %>% 
    purrr::map(
      ~ .x %>% .[stringr::str_count(.) >= min_len & stringr::str_count(.) <= max_len]) 
  
  # --- N-term methionine removed ---
  peps_m <- fasta_dbm %>% 
    purrr::map(~ .x %>% 
                 stringr::str_split("@", simplify = TRUE) %>% 
                 keep_one_to_n(max_miss)) %>% 
    purrr::map(~ concat_peps(.x, max_miss, include_cts = FALSE)) %>% 
    purrr::map(
      ~ .x %>% .[stringr::str_count(.) >= min_len & stringr::str_count(.) <= max_len]) 
  
  peps[inds_m] <- purrr::map2(peps_m, peps[inds_m], `c`)
  
  aa_masses <- ext_aa_masses(fixedmods, varmods, max_ncomb_mods)


  pep_masses <- aa_masses %>% 
    purrr::map(subpeps_by_vmods, peps) %>% 
    purrr::map2(aa_masses, mcalc_monopep, aa_masses[[1]], digits)

  # Anywhere varmods
  # one to multiple mods

  out <- purrr::map2(aa_masses, pep_masses, ~ {
    attr(.x, "data") <- .y
    return(.x)
  })

  dir.create(gsub("(^.*/).*$", "\\1", out_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, file.path(out_path))
  
  invisible(out)
}


