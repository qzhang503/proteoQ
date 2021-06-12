#' Reads a file in fasta format
#'
#' Reads a file in fasta format by line.
#'
#' @param file A character string to the name of a protein fasta file.
#' @param acc_pattern A regular expression describing the pattern to separate
#'   the header lines of fasta entries. The default is to separate a header and
#'   keep the character string before the first space where the so kept will be
#'   used as the name of an entry. The character ">" at the beginning of the
#'   header will also be removed.
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
#'                     ">(.{50}).*")
#' head(names(fasta))
#'
#' # uniprot_acc
#' fasta <- read_fasta("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                     ">..\\|([^\\|]+)\\|.*")
#' head(names(fasta))
#'
#' # use every in the header
#' fasta <- read_fasta("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                     ">(.*)")
#' head(names(fasta))
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @seealso \code{\link{write_fasta}}
#' @export
read_fasta <- function (file, acc_pattern = ">([^ ]+?) .*", comment_char = "") {
  lines <- readLines(file)
  
  if (nchar(comment_char) > 0) {
    lines <- lines %>% .[!grepl(paste0("^", comment_char), .)]
  }
  
  headers <- grep(">", lines)
  begins <- headers + 1
  ends <- (headers[-1] - 1) %>% `c`(length(lines))
  
  seqs <- purrr::map2(begins, ends, ~ lines[.x : .y] %>% 
                        purrr::reduce(paste0), lines)
  hdrs <- lines[headers]
  
  db <- purrr::map2(seqs, hdrs, ~ {
    attr(.x, "header") <- .y
    return(.x)
  })
  
  names(db) <- hdrs %>% gsub(acc_pattern, "\\1", .)

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


#' Loads fasta (with parsing rule).
#'
#' The length of \code{acc_type} needs to match the length of \code{fasta};
#' otherwise, the first value will be used for all \code{fasta} files.
#'
#' @param acc_type Character string(s); the types of protein accessions in one
#'   of c("uniprot_acc", "uniprot_id", "refseq_acc", "other"). For custom names,
#'   the corresponding regular expressions need to be supplied via argument
#'   \code{acc_pattern}.
#' @param acc_pattern Regular expression(s) describing the patterns to separate
#'   the header lines of fasta entries. At the \code{NULL} default, the pattern
#'   will be automated when \code{acc_type} are among c("uniprot_acc",
#'   "uniprot_id", "refseq_acc", "other").
#' @inheritParams read_fasta
#' @inheritParams splitPSM
#' @examples
#' \donttest{
#' fasta_db <- load_fasta2(
#'               c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                 "~/proteoQ/dbs/fasta/crap/crap.fasta"),
#'               c("uniprot_acc", "other")
#' )
#'
#' # Need `acc_pattern` as "crap" is not one of the default acc_type
#' load_fasta2(
#'    c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'      "~/proteoQ/dbs/fasta/crap/crap.fasta"),
#'    c("uniprot_acc", "crap")
#' )
#'
#' # ok
#' fasta_db2 <- load_fasta2(
#'                c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                         "~/proteoQ/dbs/fasta/crap/crap.fasta"),
#'                c("uniprot_acc", "crap"),
#'                c("^>..\\|([^\\|]+)\\|[^\\|]+", "(.*)")
#' )
#'
#' fasta_db3 <- load_fasta2(
#'                c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta",
#'                         "~/proteoQ/dbs/fasta/crap/crap.fasta"),
#'                c("my_acc", "crap"),
#'                c("^>..\\|([^\\|]+)\\|[^\\|]+", "(.*)")
#' )
#'
#' stopifnot(identical(fasta_db, fasta_db2),
#'           identical(fasta_db, fasta_db3))
#' }
#' @export
load_fasta2 <- function (fasta = NULL, acc_type = NULL, acc_pattern = NULL) {
  if (is.null(fasta)) {
    stop("FASTA file(s) are required.", call. = FALSE)
  }
  
  if (!all(file.exists(fasta))) {
    stop("Missing FASTA file(s): \n", 
         paste(fasta %>% .[!file.exists(.)], collapse = "\n"), 
         call. = FALSE)
  }
  
  len_f <- length(fasta)
  len_a <- length(acc_type)
  len_p <- length(acc_pattern)
  
  if (len_f < len_a) {
    stop("More accession types than fasta files.", 
         call. = FALSE)
  }
  
  if (len_f < len_p) {
    stop("More acc_pattern types than fasta files.", 
         call. = FALSE)
  }
  
  if (len_a > 0L && len_a < len_f) {
    warning("More fasta files than accession types; ", 
            "the first accession type will be used for all fastas.", 
            call. = FALSE)
    acc_type <- rep(acc_type[1], len_f)
  }
  
  if (len_p > 0L && len_p < len_f) {
    warning("More fasta files than acc_pattern expressions; ", 
            "the first acc_pattern expression will be used for all fastas.", 
            call. = FALSE)
    acc_pattern <- rep(acc_pattern[1], len_f)
  }
  
  if (! (is.null(acc_type) || is.null(acc_pattern))) {
    acc_type <- acc_type
    acc_pattern <- acc_pattern
  } else if (is.null(acc_type) && is.null(acc_pattern)) {
    acc_type <- rep("other", len_f)
    acc_pattern <- rep("(.*)", len_f)
  } else if (!is.null(acc_type)) {
    acc_pattern <- map_chr(acc_type, find_acc_pattern)
  } else {
    acc_type <- map_chr(acc_pattern, find_acc_type)
  }
  
  stopifnot(length(acc_pattern) == len_f)

  purrr::map2(fasta, acc_pattern, ~ read_fasta(.x, .y)) %>% 
    do.call(`c`, .) %>% 
    `names<-`(gsub(">", "", names(.))) %>% 
    .[!duplicated(names(.))]
}


#' Helper for \link{load_fasta2}.
#' 
#' Not used for custom acc_type, i.e. acc_type = "my_acctype".
#' 
#' @inheritParams load_fasta2
find_acc_pattern <- function (acc_type) {
  stopifnot(length(acc_type) == 1L)
  stopifnot(acc_type %in% c("uniprot_acc", "uniprot_id", "refseq_acc", "other"))
  
  if (acc_type == "uniprot_acc") {
    acc_pattern <- "^>..\\|([^\\|]+)\\|[^\\|]+"
  } else if (acc_type == "uniprot_id") {
    acc_pattern <- "^>..\\|[^\\|]+\\|([^ ]+) .*"
  } else if (acc_type == "refseq_acc") {
    acc_pattern <- "^>([^ ]+?) .*"
  } else if (acc_type == "other") {
    acc_type <- "other"
    acc_pattern <- "(.*)"
  } else {
    stop("Unknown `acc_type`.", 
         call. = FALSE)
  }
  
  invisible(acc_pattern)
}


#' Helper for \link{load_fasta2}.
#' 
#' Not used for custom acc_pattern, i.e. acc_pattern = "...".
#' 
#' @inheritParams load_fasta2
find_acc_type <- function (acc_pattern) {
  stopifnot(length(acc_pattern) == 1L)
  
  pat_upacc <- "^>..\\|([^\\|]+)\\|[^ ]+?"
  pat_upid <- "^>..\\|[^\\|]+\\|([^ ]+?)"
  pat_rsacc <- "^>([^ ]+?) "
  pat_other <- "(.*)"
  
  stopifnot(acc_pattern %in% c("pat_upacc", "pat_upid", "pat_rsacc", "pat_other"))
  
  if (acc_pattern == pat_upacc) {
    acc_type <- "uniprot_acc"
  } else if (acc_pattern == pat_upid) {
    acc_type <- "uniprot_id"
  } else if (acc_pattern == pat_rsacc) {
    acc_type <- "refseq_acc"
  } else if (acc_pattern == pat_other) {
    acc_type <- "other"
  } else {
    stop("Unknown `acc_pattern`.", 
         call. = FALSE)
  }
  
  invisible(acc_type)
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
#' @param maxn_sites_per_vmod Integer; the maximum number of combinatorial
#'   variable modifications per site in a per peptide sequence.
#' @param maxn_vmods_per_pep The maximum number of variable modifications per
#'   peptide.
#' @inheritParams add_fixvar_masses
#' @inheritParams calc_pepmasses
#' @examples
#' \dontrun{
#' aa_masses_all <- calc_aamasses()
#'
#' # Error if peptide sequences not 'compatible' with the given 'aa_masses'
#' x <- mcalc_monopep(list(protein_a = c("AAIMDWFDMGKEFSGNPNIK", "MAAIDWFDGKEFSGNPNIK"),
#'                         protein_b = c("AAIDWFDGKEMFSGNPNIK")),
#'                    aa_masses_all[[11]])
#' }
mcalc_monopep <- function (aa_seqs, aa_masses, 
                           include_insource_nl = FALSE, 
                           maxn_vmods_per_pep = 5, 
                           maxn_sites_per_vmod = 3, 
                           parallel = TRUE, 
                           digits = 5) {
  options(digits = 9)
  
  if (parallel) {
    n_cores <- detectCores()
    
    aa_seqs <- suppressWarnings(split(aa_seqs, seq_len(n_cores)))
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    
    parallel::clusterExport(cl, list("%>%"), 
                            envir = environment(magrittr::`%>%`))
    parallel::clusterExport(cl, list("calc_monopep"), 
                            envir = environment(proteoQ:::calc_monopep))
    parallel::clusterExport(cl, list("calc_prot_pepmasses"), 
                            envir = environment(proteoQ:::calc_prot_pepmasses))
    parallel::clusterExport(cl, list("calc_prots_pepmasses"), 
                            envir = environment(proteoQ:::calc_prots_pepmasses))
    
    out <- parallel::clusterApply(cl, aa_seqs, calc_prots_pepmasses, 
                                  aa_masses = aa_masses, 
                                  include_insource_nl = include_insource_nl, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  digits = digits) %>% 
      purrr::flatten()
    
    parallel::stopCluster(cl)
  } else {
    out <- aa_seqs %>% 
      purrr::map(~ { # by prot_acc
        x <- .x %>% 
          purrr::map(calc_monopep, 
                     aa_masses = aa_masses, 
                     include_insource_nl = include_insource_nl,
                     maxn_vmods_per_pep = maxn_vmods_per_pep, 
                     maxn_sites_per_vmod = maxn_sites_per_vmod, 
                     digits = digits) %>% # by peptides
          unlist(use.names = TRUE)
      })
  }
  
  attr(aa_masses, "data") <- out
  
  message("\tCompleted peptide masses: ", 
          paste(attributes(aa_masses)$fmods, 
                attributes(aa_masses)$vmods, 
                collapse = ", "))
  
  rm(aa_seqs, out)
  gc()
  
  invisible(aa_masses)
}


#' Helper function for parallel calculations peptide masses by proteins.
#' 
#' For each split of multiple proteins.
#' 
#' @param prot_peps Lists of peptides under a proteins.
#' @inheritParams calc_monopep
calc_prots_pepmasses <- function (aa_seqs, aa_masses, include_insource_nl = FALSE, 
                                  maxn_vmods_per_pep, maxn_sites_per_vmod, digits) {
  purrr::map(aa_seqs, ~ {
    prot_peps <- .x
    
    calc_prot_pepmasses(prot_peps = prot_peps, 
                        aa_masses = aa_masses, 
                        include_insource_nl = include_insource_nl, 
                        maxn_vmods_per_pep = maxn_vmods_per_pep, 
                        maxn_sites_per_vmod = maxn_sites_per_vmod, 
                        digits = digits)
  })
}

#' Helper function for parallel calculations peptide masses by proteins.
#' 
#' For single protein.
#' 
#' @param prot_peps Lists of peptides under a proteins.
#' @inheritParams calc_monopep
calc_prot_pepmasses <- function (prot_peps, aa_masses, 
                                 include_insource_nl, 
                                 maxn_vmods_per_pep, 
                                 maxn_sites_per_vmod, 
                                 digits) {
  prot_peps %>% 
    purrr::map(calc_monopep, 
               aa_masses = aa_masses, 
               include_insource_nl = include_insource_nl, 
               maxn_vmods_per_pep = maxn_vmods_per_pep, 
               maxn_sites_per_vmod = maxn_sites_per_vmod, 
               digits = digits) %>% # by peptides
    unlist(use.names = TRUE)
}


#' Calculates the mono-isotopic mass of a peptide sequence.
#'
#' Typically coupled to \link{subpeps_by_vmods} for automatic dispatching of
#' peptide sequences by sets of fixed and variable modifications. For manual
#' calculations, uses \link{calc_monopeptide}.
#'
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param include_insource_nl Logical Logical; if TRUE, includes MS1 precursor
#'   masses with the losses of neutral species prior to MS2 fragmentation. The
#'   default is FALSE.
#'
#'   Note that there is no combination of neutral losses for fixed modifications
#'   at the precursor levels. Changes from fixed to variable modifications for
#'   complete combinations. 
#' @inheritParams mcalc_monopep
#' @import purrr
#' @importFrom stringr str_split
calc_monopep <- function (aa_seq, aa_masses, 
                          include_insource_nl = FALSE, 
                          maxn_vmods_per_pep = 5, 
                          maxn_sites_per_vmod = 3, 
                          digits = 5) {
  
  if (is.na(aa_seq)) return(NULL)
  
  aas <- aa_seq %>% str_split("", simplify = TRUE)
  type <- attr(aa_masses, "type", exact = TRUE)
  
  if (type == "amods- tmod- vnl- fnl-") {
    return(
      aas %>% 
        aa_masses[.] %>% 
        sum() %>% 
        `+`(aa_masses["N-term"]) %>% 
        `+`(aa_masses["C-term"]) %>% 
        setNames(aa_seq) %>% 
        round(digits = digits)
    )
  }
  
  ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
  ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
  
  if (type == "amods- tmod+ vnl- fnl-") {
    return(calcpep_a0_t1_nl0(aa_seq, aas, aa_masses, ntmod, ctmod, digits))
  } 
  
  # --- combinatorial ---
  fmods_ps <- attr(aa_masses, "fmods_ps", exact = TRUE)
  vmods_ps <- attr(aa_masses, "vmods_ps", exact = TRUE)
  fmods_nl <- attr(aa_masses, "fmods_nl", exact = TRUE)
  vmods_nl <- attr(aa_masses, "vmods_nl", exact = TRUE)
  amods <- attr(aa_masses, "amods", exact = TRUE)
  tmod <- attr(aa_masses, "tmod", exact = TRUE)
  
  # at least one of `vmods_combi`, `fnl_combi`, `vnl_combi`
  vmods_combi <- unique_mvmods(amods = amods, ntmod = ntmod, ctmod = ctmod, 
                               aas = aas, aa_masses = aa_masses, 
                               maxn_vmods_per_pep = maxn_vmods_per_pep, 
                               maxn_sites_per_vmod = maxn_sites_per_vmod, 
                               digits = digits) %>% 
    find_intercombi()
  
  # (1) if called manually (without automatic sequences dispatching): 
  #   may be zero-intersect between an `aa_seq` and a given `aa_masses`.
  # 
  # (2) "MAEEQGR" - vmods_combi is NULL even with dispatching: 
  #   cannot be simultaneous `Oxidation (M)` and `Acetyl (Protein N-term)`
  #   (the following only handles early for (10) "amods+ tmod+ vnl+"; 
  #    additional handling of "(8) amods+ tmod+ vnl-" and 
  #    (12) "amods+ tmod+ vnl-" by corresponding functions)
  if (is_empty(vmods_combi) && !is_empty(vmods_nl)) return(NULL)
  
  if (is_empty(vmods_nl) || !include_insource_nl) {
    vnl_combi <- NULL
  } else {
    vnl_combi <- map(vmods_combi, ~ expand.grid(vmods_nl[.x]) %>% t())
  }
  
  # No `fmods_combi` as only one [NC]-term -> no `map` in `expand.grid`.
  # Thus unlike length(vmods_combi) == length(vnl_combi), 
  #   length(fnl_combi) == 1 table (with n columns of NLs); 
  # nrow(fnl_combi) == length(fmods_nl) == number of sites with fixed mods
  # 
  # For Anywhere fmods -> make it vmods for complete neutral losses; otherwise 
  #       [,1]      [,2]
  #    M 147.0354 83.037115
  # no combination of 147 and 83 (uniformly 147 or 83)
  
  # The difference also affect the mapping of peptides in mass calculations: 
  # (a) vnl+: mcalcpepmass_a1_vnl1 -> calcpepmass_vnl1
  # (b) fnl+: mcalcpepmass_a1_fnl1 -> calcpepmass_fnl1
  
  if (is_empty(fmods_nl) || !include_insource_nl) {
    fnl_combi <- NULL
  } else {
    fnl_combi <- expand.grid(fmods_nl) %>% t()
  }
  
  out <- switch(type, 
            # "amods- tmod- vnl- fnl-" = calcpep_0(aa_seq, aas, aa_masses, digits), 
            # "amods- tmod+ vnl- fnl-" = calcpep_a0_t1_nl0(aa_seq, aas, aa_masses, ntmod, ctmod, digits), 
            # "amods- tmod+ vnl+ fnl-" = calcpep_a0_t1_vnl1(vnl_combi, aa_seq, aas, aa_masses, ntmod, ctmod, digits), 
            # "amods- tmod- vnl+ fnl-" = calcpep_a0_t0_vnl1(vmods_combi, vnl_combi, NULL, NULL, NULL, aa_seq, aas, aa_masses, digits),
            
            "amods- tmod+ vnl- fnl+" = 
              calcpep_a0_t1_fnl1(fnl_combi, aa_seq, aas, aa_masses, ntmod, ctmod, digits), 
            "amods- tmod- vnl- fnl+" = 
              calcpep_a0_t0_fnl1(fnl_combi, aa_seq, aas, aa_masses, NULL, NULL, digits), 
            
            "amods+ tmod- vnl- fnl-" = 
              calcpep_a1_t0_nl0(vmods_combi, NULL, amods, NULL, NULL, aa_seq, aas, aa_masses, digits), 
            "amods+ tmod+ vnl- fnl-" = 
              calcpep_a1_t1_nl0(vmods_combi, NULL, amods, ntmod, ctmod, aa_seq, aas, aa_masses, digits), 
            
            "amods+ tmod- vnl+ fnl-" = 
              calcpep_a1_t0_nl1(vmods_combi, vnl_combi, amods, NULL, NULL, aa_seq, aas, aa_masses, digits), 
            "amods+ tmod+ vnl+ fnl-" = 
              calcpep_a1_t1_nl1(vmods_combi, vnl_combi, amods, ntmod, ctmod, aa_seq, aas, aa_masses, digits), 
            
            "amods+ tmod- vnl- fnl+" = 
              calcpep_a1_t0_fnl1(vmods_combi, fnl_combi, amods, NULL, NULL, aa_seq, aas, aa_masses, digits), 
            "amods+ tmod+ vnl- fnl+" = 
              calcpep_a1_t1_fnl1(vmods_combi, fnl_combi, amods, ntmod, ctmod, aa_seq, aas, aa_masses, digits), 
            
            # "amods- tmod- vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
            # "amods+ tmod- vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
            # "amods- tmod+ vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
            # "amods+ tmod+ vnl+ fnl+" = foo(aa_seq, aa_masses, digits), 
            # message("Ignores `vnl+ fnl+`.")
            NULL)
}


#' Calculates the mono-isotopic mass of a peptide sequence.
#'
#' For direct uses from an R console (with trade-offs in speed).
#'
#' @inheritParams calc_monopep
#' @inheritParams calc_aamasses
#'
#' @examples
#' \donttest{
#' ## No variable modifications
#' # (1)
#' x <- calc_monopeptide("MAKEMASSPECFUN", 
#'                       fixedmods = NULL, 
#'                       varmods = NULL)
#' 
#' x$mass
#' 
#' # (2-a) 
#' x <- calc_monopeptide("MAKEMASSPECFUN", 
#'                       fixedmods = "Oxidation (M)", 
#'                       varmods = NULL)
#' 
#' x$mass
#' 
#' # (2-b) no combinatorial NL for fixed modifications
#' # (see 3-b) for full combinations
#' x <- calc_monopeptide("MAKEMASSPECFUN", 
#'                       fixedmods = "Oxidation (M)", 
#'                       varmods = NULL, 
#'                       include_insource_nl = TRUE)
#' 
#' x$mass
#' 
#' ## With variable modifications
#' # (3-a)
#' x <- calc_monopeptide("MAKEMASSPECFUN", 
#'                       fixedmods = NULL, 
#'                       varmods = "Oxidation (M)")
#' 
#' x$mass
#' # x$vmods_ps
#' 
#' # (3-b)
#' x <- calc_monopeptide("MAKEMASSPECFUN", 
#'                       fixedmods = NULL, 
#'                       varmods = "Oxidation (M)", 
#'                       include_insource_nl = TRUE)
#' 
#' x$mass
#' 
#' # (4-a)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       c("TMT6plex (N-term)", 
#'                         "TMT6plex (K)", 
#'                         "Carbamidomethyl (C)"),
#'                       c("Acetyl (N-term)", 
#'                         "Gln->pyro-Glu (N-term = Q)", 
#'                         "Oxidation (M)"))
#'                       
#' x$mass
#' 
#' # The N-term M realizes with acetylation
#' x$vmods_ps[[5]]
#' 
#' # The N-term M realizes with TMT
#' x$vmods_ps[[6]]
#' 
#' 
#' # (4-b)
#' x <- calc_monopeptide("MAKEMASSPECFUN",
#'                       c("TMT6plex (N-term)", 
#'                         "TMT6plex (K)", 
#'                         "Carbamidomethyl (C)"),
#'                       c("Acetyl (N-term)", 
#'                         "Gln->pyro-Glu (N-term = Q)", 
#'                         "Oxidation (M)"), 
#'                         include_insource_nl = TRUE)
#'                       
#' x$mass
#' }
#' @export
calc_monopeptide <- function (aa_seq, fixedmods, varmods, 
                              include_insource_nl = FALSE, 
                              maxn_vmods_setscombi = 64,
                              maxn_vmods_per_pep = Inf, 
                              maxn_sites_per_vmod = Inf, 
                              digits = 5) {
  options(digits = 9)
  
  aa_masses_all <- calc_aamasses(fixedmods, varmods, maxn_vmods_setscombi)
  
  peps <- purrr::map(aa_masses_all, subpeps_by_vmods, aa_seq) %>% 
    purrr::flatten()
  
  oks <- purrr::map_lgl(peps, ~ !purrr::is_empty(.x))
  peps <- peps[oks]
  aa_masses_all <- aa_masses_all[oks]

  ms <- purrr::map2(peps, aa_masses_all, ~ {
    calc_monopep(.x, .y, 
                 include_insource_nl = include_insource_nl, 
                 maxn_vmods_per_pep = maxn_vmods_per_pep, 
                 maxn_sites_per_vmod = maxn_sites_per_vmod, 
                 digits = digits)
  })
  
  attrs <- purrr::map(aa_masses_all, attributes)
  vmods_ps <- map(attrs, `[[`, "vmods_ps")
  
  list(mass = ms, vmods_ps = vmods_ps)
}


#' Helper in calculating the mass of an unmodified peptide.
#' 
#' (1) "amods- tmod- vnl- fnl-".
#' 
#' @inheritParams calc_monopep
#' @inheritParams calcpep_a1_t1_nl1
calcpep_0 <- function (aa_seq, aas, aa_masses, digits) {
  aas %>% 
    aa_masses[.] %>% 
    sum() %>% 
    `+`(aa_masses["N-term"]) %>% 
    `+`(aa_masses["C-term"]) %>% 
    setNames(aa_seq) %>% 
    round(digits = digits)
}


#' Helper in calculating peptide masses.
#' 
#' (2) "amods- tmod+ vnl- fnl-".
#' 
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a0_t1_nl0 <- function (aa_seq, aas, aa_masses, ntmod, ctmod, digits) {
  out <- calcpepmass_nl0(aas, aa_masses, ntmod, ctmod) 
  
  out %>% 
    setNames(aa_seq) %>% 
    round(digits = digits)
}


#' Helper in calculating peptide masses.
#' 
#' (3) "amods- tmod+ vnl+ fnl-".
#' 
#' Not existed, `amods-` -> `vnl+` is NULL.
#' 
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a0_t1_vnl1 <- function (nl_combi, aa_seq, aas, aa_masses, ntmod, ctmod, 
                                digits) {
  message("Combination not possible: `amods- tmod+ vnl+`.")
}

#' Helper in calculating peptide masses.
#'
#' (4) "amods- tmod- vnl+ fnl-". 
#' 
#' Not existed, `amods-` -> `vnl+` is NULL.
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a0_t0_vnl1 <- function (vmods_combi, nl_combi, amods, ntmod, ctmod, 
                                aa_seq, aas, aa_masses, digits) {
  message("Combination not possible: `amods- tmod- vnl+`.")
}


#' Helper in calculating peptide masses.
#' 
#' (5) "amods- tmod+ vnl- fnl+".
#' 
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a0_t1_fnl1 <- function (nl_combi, aa_seq, aas, aa_masses, ntmod, ctmod, 
                                digits) {
  # `amods-` -> no need to `update_aas_anywhere`
  
  # `nl_combi` is `fnl_combi`
  # go by columns in `nl_combi`
  
  if (is.null(nl_combi)) {
    out <- calcpep_a0_t1_nl0(aa_seq, aas, aa_masses, ntmod, ctmod, digits)
  } else {
    out <- calcpepmass_fnl1(aas, nl_combi, aa_masses, ntmod, ctmod) %>% 
      setNames(rep(aa_seq, length(.))) %>% 
      round(digits = digits)
  }
  
  invisible(out)
}


#' Helper in calculating peptide masses.
#' 
#' (6) "amods- tmod- vnl- fnl+".
#' 
#' Identical to (5) "amods- tmod+ vnl- fnl+" other than `ntmod` and `ctmod`.
#'
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a0_t0_fnl1 <- function (nl_combi, aa_seq, aas, aa_masses, 
                                ntmod = NULL, ctmod = NULL, digits) {
  
  if (is.null(nl_combi)) {
    out <- calcpep_0(aa_seq, aas, aa_masses, digits)
  } else {
    out <- calcpepmass_fnl1(aas, nl_combi, aa_masses, 
                            ntmod = NULL, ctmod = NULL) %>% 
      setNames(rep(aa_seq, length(.))) %>% 
      round(digits = digits)
  }
  
  invisible(out)
}


#' Helper in calculating peptide masses.
#' 
#' (7) "amods+ tmod- vnl- fnl-".
#'
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a1_t0_nl0 <- function (vmods_combi, nl_combi = NULL, amods, 
                               ntmod = NULL, ctmod = NULL, aa_seq, aas, 
                               aa_masses, digits) {
  
  # `amods+` -> `update_aas_anywhere` by `vmods_combi`
  
  aases <- vmods_combi %>% 
    map(update_aas_anywhere, aas, amods)
  
  out <- aases %>% 
    map(calcpepmass_nl0, aa_masses, ntmod, ctmod) %>% 
    unlist(use.names = FALSE) # explicit unlist instead of map_dbl
  
  out <- out %>% 
    unique() %>% 
    setNames(rep(aa_seq, length(.))) %>% 
    round(digits = digits)
}

#' Helper in calculating peptide masses.
#'
#' (8) "amods+ tmod+ vnl- fnl-".
#'
#' Identical to "(7) amods+ tmod- vnl- fnl-" other than `ntmod` and `ctmod`.
#'
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a1_t1_nl0 <- function (vmods_combi, nl_combi = NULL, amods, 
                               ntmod, ctmod, aa_seq, aas, aa_masses, 
                               digits) {
  
  # A case of NULL `out` (and cannot setNames)
  # conflicts between `amods+` and `tmod+`
  # NLTLALEALVQLR: the only N or the N-term
  # if (is_empty(vmods_combi)) return(NULL)
  
  aases <- vmods_combi %>% 
    map(update_aas_anywhere, aas, amods)
  
  out <- aases %>% 
    map(calcpepmass_nl0, aa_masses, ntmod, ctmod) %>% 
    unlist(use.names = FALSE)
  
  out <- out %>% unique() 

  # No conflict between `amods+` and `tmod+`
  if (!is_empty(out)) {
    out <- out %>% 
      setNames(rep(aa_seq, length(.))) %>% 
      round(digits = digits)
  }
  
  invisible(out)
}


#' Helper in calculating peptide masses.
#' 
#' (9) "amods+ tmod- vnl+ fnl-".
#'
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a1_t0_nl1 <- function (vmods_combi, nl_combi, amods, ntmod = NULL, 
                               ctmod = NULL, aa_seq, aas, aa_masses, 
                               digits) {
  if (is.null(nl_combi)) {
    out <- calcpep_a1_t0_nl0(vmods_combi, nl_combi = NULL, amods, 
                             ntmod = NULL, ctmod = NULL, aa_seq, aas, 
                             aa_masses, digits)
  } else {
    out <- map2(vmods_combi, nl_combi, mcalcpepmass_a1_vnl1, 
                amods, ntmod = NULL, ctmod = NULL, aas, aa_masses, digits) %>% 
      unlist(use.names = FALSE) # no map2_dbl
    
    out <- out %>% 
      unique() %>% 
      setNames(rep(aa_seq, length(.))) %>% 
      round(digits = digits)
  }
  
  invisible(out)
}


#' Helper in calculating peptide masses.
#' 
#' (10) "amods+ tmod+ vnl+ fnl-".
#'
#' @param vmods_combi Lists of variable modifications.
#' @param nl_combi Lists of combinations of neutral losses for corresponding
#'   \code{vmods_combi}. Each list contains a table where each column
#'   corresponds to a set of neutral loss. The first column corresponds to the
#'   combination without NLs.
#' @param amods The parsed \emph{variable} \code{Anywhere} modifications from
#'   \code{aa_masses}. The value can be \code{NULL}, e.g., at "amods- tmod- vnl-
#'   fnl+".
#' @param ntmod The parsed \emph{variable} \code{N-term} modifications from
#'   \code{aa_masses}.
#' @param ctmod The parsed \emph{variable} \code{C-term} modifications from
#'   \code{aa_masses}.
#' @param aas \code{aa_seq} split in a sequence of LETTERS.
#' @param aa_masses A named list containing the (mono-isotopic) masses of amino
#'   acid residues.
#' @inheritParams calc_monopep
#' @inheritParams calc_pepmasses
calcpep_a1_t1_nl1 <- function (vmods_combi, nl_combi, amods, 
                               ntmod, ctmod, aa_seq, aas, aa_masses, 
                               digits) {
  if (is.null(nl_combi)) {
    out <- calcpep_a1_t1_nl0(vmods_combi, nl_combi = NULL, amods, 
                             ntmod, ctmod, aa_seq, aas, 
                             aa_masses, digits)
  } else {
    out <- map2(vmods_combi, nl_combi, mcalcpepmass_a1_vnl1, 
                amods, ntmod, ctmod, aas, aa_masses, digits) %>% 
      unlist(use.names = FALSE) # no map2_dbl
    
    ## A case of NULL `out`
    # conflicts btw. amods+ and tmod+
    #   i.e. NLTLALEALVQLR: the only N on the N-term
    #   cannot be both `Deamidated (N)` and `Acetyl (N-term)`
    # if (is_empty(vmods_combi)) return(NULL)
    
    out <- out %>% unique() 
    
    if (!is_empty(out)) {
      out <- out %>% 
        setNames(rep(aa_seq, length(.))) %>% 
        round(digits = digits)
    }
  }

  invisible(out)
}


#' Helper in calculating peptide masses.
#' 
#' (11) "amods+ tmod- vnl- fnl+"
#' 
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a1_t0_fnl1 <- function (vmods_combi, nl_combi, amods, ntmod = NULL, 
                                ctmod = NULL, aa_seq, aas, aa_masses, 
                                digits) {
  
  # `nl_combi` is `fnl_combi` and is a length(1) table
  # (length(vmods_combi) >= 1; length(fnl_combi) == 1 table of n columns)
  
  if (is.null(nl_combi)) {
    out <- calcpep_a1_t0_nl0(vmods_combi, nl_combi = NULL, amods, 
                             ntmod = NULL, ctmod = NULL, aa_seq, aas, 
                             aa_masses, digits)
  } else {
    out <- map(vmods_combi, mcalcpepmass_a1_fnl1, nl_combi, 
               amods, ntmod = NULL, ctmod = NULL, aas, aa_masses, digits) %>% 
      unlist(use.names = FALSE)
    
    out <- out %>% 
      unique() %>% 
      setNames(rep(aa_seq, length(.))) %>% 
      round(digits = digits)
  }

  invisible(out)
}


#' Helper in calculating peptide masses.
#' 
#' (12) "amods+ tmod+ vnl- fnl+"
#' 
#' Identical to (11) only differed by `ntmod` and `ctmod`.
#'
#' @inheritParams calcpep_a1_t1_nl1
calcpep_a1_t1_fnl1 <- function (vmods_combi, nl_combi, amods, 
                                ntmod, ctmod, aa_seq, aas, aa_masses, 
                                digits) {
  
  if (is.null(nl_combi)) {
    out <- calcpep_a1_t1_nl0(vmods_combi, nl_combi = NULL, amods, 
                             ntmod, ctmod, aa_seq, aas, aa_masses, 
                             digits)
  } else {
    # Again can be NULL `out` at `amods+` and `tmod+`
    out <- map(vmods_combi, mcalcpepmass_a1_fnl1, nl_combi, 
               amods, ntmod, ctmod, aas, aa_masses, digits) %>% 
      unlist(use.names = FALSE)
    
    out <- out %>% unique() 
    
    if (!is_empty(out)) {
      out <- out %>% 
        setNames(rep(aa_seq, length(.))) %>% 
        round(digits = digits)
    }
  }

  invisible(out)
}


#' Helper: multiple mode of peptide-mass calculation.
#'
#' Specific to "amods+" and "vnl+" (and fnl-).
#'
#' \code{amods} is not empty; either \code{vnl} or \code{fnl} is not empty.
#' @inheritParams calcpep_a1_t1_nl1
mcalcpepmass_a1_vnl1 <- function (vmods_combi, nl_combi, amods, ntmod, ctmod, 
                                 aas, aa_masses, digits) {
  # conflicts btw. vmods+ and tmod+
  if (is_empty(vmods_combi)) return(NULL)
  
  # since `amods+`
  aas <- update_aas_anywhere(vmods_combi, aas, amods)
  
  # no `map` since by the columns of nl_combi
  out <- calcpepmass_vnl1(aas, nl_combi, aa_masses, ntmod, ctmod)
}


#' Helper: multiple mode of peptide-mass calculation.
#'
#' Specific to "amods+" and "fnl+" (and vnl-).
#'
#' \code{amods} is not empty; either \code{vnl} or \code{fnl} is not empty.
#' @inheritParams calcpep_a1_t1_nl1
mcalcpepmass_a1_fnl1 <- function (vmods_combi, nl_combi, amods, ntmod, ctmod, 
                                  aas, aa_masses, digits) {
  # conflicts btw. vmods+ and tmod+
  # if (is_empty(vmods_combi)) return(NULL)
  
  # since `amods+`
  aas <- update_aas_anywhere(vmods_combi, aas, amods)
  
  # no `map` since by the columns of `nl_combi`
  out <- calcpepmass_fnl1(aas, nl_combi, aa_masses, ntmod, ctmod)
}


#' Helper of peptide-mass calculation..
#' 
#' Without NLs: "fnl-" or "vnl-".
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @return A numeric scalar
calcpepmass_nl0 <- function (aas, aa_masses, ntmod, ctmod) {
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
  
  invisible(out)
}


#' Helper of peptide-mass calculation..
#'
#' With NLs: "vnl+".
#'
#' The calculation goes through the columns in \code{nl_combi}. 
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @return A numeric vector
calcpepmass_vnl1 <- function (aas, nl_combi, aa_masses, ntmod, ctmod) {
  col_sums <- colSums(nl_combi)
  upds <- aas %in% rownames(nl_combi)
  out0 <- sum(aa_masses[aas[!upds]]) 

  len <- ncol(nl_combi)
  out <- vector("numeric", len)

  for (i in 1:len) {
    out[[i]] <- out0 + col_sums[i]
    
    if (is_empty(ntmod) && is_empty(ctmod)) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses["N-term"]) %>% 
        `+`(aa_masses["C-term"]) 
    } else if (!(is_empty(ntmod) || is_empty(ctmod))) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses[names(ntmod)]) %>% 
        `+`(aa_masses[names(ctmod)])
    } else if (!is_empty(ntmod)) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses[names(ntmod)]) %>% 
        `+`(aa_masses["C-term"]) 
    } else if (!is_empty(ctmod)) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses["N-term"]) %>% 
        `+`(aa_masses[names(ctmod)])
    }
  }
  
  out <- unique(out)
}


#' Helper of peptide-mass calculation..
#'
#' With NLs: "fnl+".
#'
#' The calculation goes through the columns in \code{nl_combi}. 
#' 
#' @inheritParams calcpep_a1_t1_nl1
#' @return A numeric vector
calcpepmass_fnl1 <- function (aas, nl_combi, aa_masses, ntmod, ctmod) {
  len <- ncol(nl_combi)
  out <- vector("numeric", len)
  nms <- rownames(nl_combi)

  for (i in 1:len) {
    # no combination of 147 for some "M" and 83 for others;
    # uniformly either 147 or 83 
    
    aa_masses[nms] <- nl_combi[, i]
    out[[i]] <- aas %>% aa_masses[.] %>% sum() 
    
    if (is_empty(ntmod) && is_empty(ctmod)) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses["N-term"]) %>% 
        `+`(aa_masses["C-term"]) 
    } else if (!(is_empty(ntmod) || is_empty(ctmod))) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses[names(ntmod)]) %>% 
        `+`(aa_masses[names(ctmod)])
    } else if (!is_empty(ntmod)) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses[names(ntmod)]) %>% 
        `+`(aa_masses["C-term"]) 
    } else if (!is_empty(ctmod)) {
      out[[i]] <- out[[i]] %>% 
        `+`(aa_masses["N-term"]) %>% 
        `+`(aa_masses[names(ctmod)])
    }
  }
  
  out <- unique(out)
}


#' Updates the labels of amino-acid residues.
#'
#' The mappings are according to fixed and variable modifications. The
#' \code{aas} in a \emph{single} sequence and \code{vmods_combi} is also a
#' \emph{single} vector.
#'
#' The mappings are for the count combination of \code{Anywhere} sites. It does
#' not apply to terminal sites (as there is only one N- or C-term).
#'
#' @inheritParams calcpep_a1_t1_nl1
update_aas_anywhere <- function (vmods_combi, aas, amods) {
  nms <- imap(amods, ~ {
    names(.x) <- .y
    .x
  }) %>% flatten() %>% 
    unlist()
  
  n_per_site <- map(names(nms), ~ sum(vmods_combi == .x)) %>% 
    `names<-`(names(nms)) %>% 
    unlist()
  
  for (i in seq_along(nms)) {
    x <- nms[i]
    y <- n_per_site[i]
    
    idxes <- which(aas == x)[1:y]
    aas[idxes] <- rep(names(y), y)
  }
  
  invisible(aas)
}


#' The unique entries of variable modifications (multiple sites)
#'
#' For all the \code{Anywhere} modifications specified in \code{amods}.
#'
#' @param amods Anywhere modifications.
#' @inheritParams calc_monopep
#' @inheritParams calcpep_a1_t1_nl1
#' @import purrr
#' @return Lists by residues in \code{amods}.
unique_mvmods <- function (amods, ntmod, ctmod, aas, aa_masses, 
                           maxn_vmods_per_pep = 5, 
                           maxn_sites_per_vmod = 3, 
                           digits = 5) {
  # (6) "amods- tmod- vnl- fnl+"
  if (is_empty(amods)) return(NULL)
  
  residue_mods <- unlist(amods, use.names = FALSE) %>% 
    `names<-`(names(amods)) %>% 
    split(., .)
  
  map(residue_mods, 
      ~ vmods_elements(aas, .x, ntmod, ctmod, 
                       aa_masses, 
                       maxn_vmods_per_pep, 
                       maxn_sites_per_vmod, 
                       digits))
}


#' Find the sets of variable modifications.
#'
#' Excluding position differences, i.e., \code{A, B} and \code{B, A} is the
#' same set.
#'
#' @param residue_mods Amino-acid residues with Unimod names. For example
#'   rownames of \code{Carbamidomethyl (M)} and \code{Oxidation (M)} and a
#'   column residues of \code{M, M}.
#' @inheritParams unique_mvmods
#' @import purrr
vmods_elements <- function (aas, 
                            residue_mods, 
                            ntmod, 
                            ctmod, 
                            aa_masses, 
                            maxn_vmods_per_pep = 5, 
                            maxn_sites_per_vmod = 3, 
                            digits = 5) {
  
  residue <- residue_mods[[1]]
  
  ns <- names(residue_mods)
  len_n <- length(ns)
  
  # no need of ps[seq_len(.x)] as the exact positions not needed
  # ps <- which(aas == residue)
  # len_p <- length(ps)
  
  len_p <- sum(aas == residue)
  
  # i.e., btw Anywhere "M" and "Acetyl N-term" where "M" on the "N-term"
  # MFGMFNVSMR cannot have three `Oxidation (M)` and `Acetyl (N-term)`
  
  if (!(is_empty(ntmod) || is_empty(ctmod))) {
    len_aas <- length(aas)
    
    if (aas[1] == residue && aas[len_aas] == residue) {
      len_p <- len_p - 2
    } else if ((aas[1] == residue) || (aas[len_aas] == residue)) {
      len_p <- len_p - 1
    }
  } else if (!is_empty(ntmod)) {
    if (aas[1] == residue) {
      len_p <- len_p - 1
    }
  } else if (!is_empty(ctmod)) {
    if (aas[length(aas)] == residue) {
      len_p <- len_p - 1
    }
  }

  if (len_p <= 0) return(list())

  if (len_p > len_n) {
    x <- 
      # will confuse gtools::combinations...
      # map((len_n + 1):len_p, ~ find_unique_sets(ps[seq_len(.x)], ns)) %>% 
      
      map((len_n + 1):len_p, ~ find_unique_sets(seq_len(.x), ns)) %>% 
      recur_flatten() %>% 
      `c`(list(ns), .)
  } else {
    x <- list(ns)
  }

  rows <- map_lgl(x, ~ length(.x) > maxn_vmods_per_pep)
  x <- x[!rows]
  
  maxn_vmod <- x %>% 
    map(table) %>% 
    map(max)
  rows <- maxn_vmod > maxn_sites_per_vmod
  x <- x[!rows]
  
  invisible(x)
}


#' Finds the sets of \code{n} elements in \code{positions}.
#' 
#' At least one occupancy for each elements in \code{ns}.
#' 
#' @param ps A vector of positions.
#' @param ns The names to be filled into \code{p}.
#' @importFrom gtools combinations
find_unique_sets <- function (ps = c(1:5), ns = c("A", "B", "C")) {
  lp <- length(ps)
  ln <- length(ns)
  r <- lp - ln
  
  if (r == 0) return(list(ns))
  
  x <- combinations(ln, r, ns, repeats = TRUE)
  
  n_row <- nrow(x)
  out <- vector("list", n_row)
  
  for (i in 1:n_row) {
    out[[i]] <- c(ns, x[i, ])
  }
  
  out
}


#' Finds the combinations across residues.
#' 
#' For multiple residues (each residue one to multiple modifications).
#' 
#' @param intra_combis The results from \link{unique_mvmods}.
find_intercombi <- function (intra_combis) {
  
  if (is_empty(intra_combis)) { # scalar
    v_out <- list()
  } else if (any(map_lgl(intra_combis, is_empty))) { # list
    v_out <- list()
  } else if (length(intra_combis) > 1L) {
    inter_combi <- expand.grid(intra_combis, KEEP.OUT.ATTRS = FALSE)
    
    nrow <- nrow(inter_combi)
    v_out <- vector("list", nrow)
    
    for (i in 1:nrow) {
      v_out[[i]] <- unlist(inter_combi[i, ], use.names = FALSE)
    }
  } else {
    v_out <- flatten(intra_combis)
  }

  invisible(v_out)
}


#' Finds the site of an AA residue.
#' 
#' @param pos_site A named value. Position in name and site in value.
find_aa_site <- function (pos_site) {
  
  # e.g. 'N-term' in `aa_masses`, not 'Q', for Gln-> pyro-Glu (N-term = Q)
  
  if (grepl("[NC]{1}-term", names(pos_site))) {
    site <- names(pos_site) %>% 
      gsub("(Protein|Any) ([NC]{1}-term)", "\\2", .)
  } else {
    site <- pos_site
  }
}


#' Calculates molecular weight of a polypeptide ([MH]+).
#'
#' @param fixedmods A character vector of fixed modifications. See also
#'   \link{parse_unimod} for grammars.
#' @param varmods A character vector of variable modifications.
#' @param maxn_vmods_setscombi Integer; the maximum number of combinatorial variable
#'   modifications and neutral losses.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
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
#' 
#' # need separate S and T
#' x <- calc_aamasses(NULL, "Phospho (ST)")
#' }
#' @export
calc_aamasses <- function (fixedmods = c("TMT6plex (K)", 
                                         "Carbamidomethyl (. = C)"), 
                           varmods = c("TMT6plex (N-term)", 
                                       "Acetyl (Protein N-term)", 
                                       "Oxidation (M)", 
                                       "Deamidated (N)", 
                                       "Gln->pyro-Glu (N-term = Q)"), 
                           maxn_vmods_setscombi = 64, 
                           mod_indexes = NULL) {
  
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
      
      .x
    })
  
  if (length(aa_masses_var2) >= maxn_vmods_setscombi) {
    warning("The number of ways of fixed and variable modifications is greater than ", 
            maxn_vmods_setscombi, ".\n", 
            "Some combinations may be dropped.", 
            call. = FALSE)
    aa_masses_var2 <- aa_masses_var2[1:maxn_vmods_setscombi]
  }
  
  if (!is.null(mod_indexes)) {
    aa_masses_var2 <- aa_masses_var2 %>% 
      purrr::map(~ {
        nms <- mod_indexes %>% 
          .[names(.) %in% names(.x)]
        
        idxes <- which(names(.x) %in% names(nms))
        names(.x)[idxes] <- nms
        
        .x
      }, mod_indexes)
    
    aa_masses_fi2 <- aa_masses_fi2 %>% 
      purrr::map(~ {
        names(attr(.x, "fmods_ps")) <- mod_indexes[fixedmods]
        names(attr(.x, "fmods_mass")) <- mod_indexes[fixedmods]
        
        .x
      })
    
    aa_masses_var2 <- aa_masses_var2 %>% 
      purrr::map(~ {
        v_nms <- names(attr(.x, "vmods_ps"))
        names(attr(.x, "vmods_ps")) <- mod_indexes[v_nms]
        names(attr(.x, "vmods_mass")) <- mod_indexes[v_nms]
        
        f_nms <- names(attr(.x, "fmods_ps"))
        names(attr(.x, "fmods_ps")) <- mod_indexes[f_nms]
        names(attr(.x, "fmods_mass")) <- mod_indexes[f_nms]
        
        .x
      })
  }
  
  aa_masses_all <- c(aa_masses_fi2, aa_masses_var2)
  
  aa_masses_all <- map(aa_masses_all, parse_aamasses)
}


#' Parses \code{aa_masses}.
#' 
#' @inheritParams add_fixvar_masses
parse_aamasses <- function (aa_masses) {
  fmods_ps <- attr(aa_masses, "fmods_ps", exact = TRUE)
  vmods_ps <- attr(aa_masses, "vmods_ps", exact = TRUE)
  
  fmods_nl <- local({
    neulosses <- attr(aa_masses, "fmods_neuloss", exact = TRUE)

    if (all(neulosses == "")) return(character())
    
    # add `0` if absent
    no_zero <- map_lgl(neulosses, ~ !any(.x == 0)) %>% 
      which()
    
    if (!is_empty(no_zero)) {
      neulosses[[no_zero]] <- c(0, neulosses[[no_zero]])
    }
    
    ## entries with NL = 0 also kept 
    ## (beneficial when calling `find_intercombi`)
    # idx <- map_lgl(neulosses, ~ length(.x) > 1)
    # neulosses <- neulosses[idx]
    
    ## In fixedmods: `M` instead of `Oxidation (M)`
    # fmods_ps <- fmods_ps[idx]
    
    names(neulosses) <- fmods_ps
    
    resids <- as.list(aa_masses)
    
    purrr::imap(neulosses, ~ {
      resids[[.y]] - .x
    })
  })
  
  vmods_nl <- local({
    neulosses <- attr(aa_masses, "vmods_neuloss", exact = TRUE)
    
    if (all(neulosses == "")) return(character())
    
    # add `0` if absent
    no_zero <- map_lgl(neulosses, ~ !any(.x == 0)) %>% which()
    
    if (!is_empty(no_zero)) {
      neulosses[[no_zero]] <- c(0, neulosses[[no_zero]])
    }
    
    ## entries of "NL = 0" kept
    # idx <- map_lgl(neulosses, ~ length(.x) > 1)
    # neulosses <- neulosses[idx]
    
    resids <- as.list(aa_masses)
    
    purrr::imap(neulosses, ~ {
      resids[[.y]] - .x
    })
  })
  
  ## variable mods
  # multiple mods to [NC]-term already excluded from aa_masses
  amods <- local({
    sites <- vmods_ps %>% 
      purrr::map(~ .x[grepl("Anywhere", names(.x))]) 
    
    empties <- sites %>% purrr::map_lgl(purrr::is_empty)
    
    sites <- sites[!empties] 
    
    sites
  })
  
  tmod <- vmods_ps %>% .[! . %in% amods]
  if (purrr::is_empty(tmod)) {
    tmod <- NULL
  } else if (tmod == "") {
    tmod <- NULL
  }
  
  # variable N-term, C-term
  ntmod <- tmod %>% .[. == "N-term"]
  ctmod <- tmod %>% .[. == "C-term"]
  
  ## fixed mods
  famods <- local({
    sites <- fmods_ps %>% 
      purrr::map(~ .x[grepl("Anywhere", names(.x))]) 
    
    empties <- sites %>% purrr::map_lgl(purrr::is_empty)
    
    sites <- sites[!empties] 
    
    sites
  })
  
  ftmod <- fmods_ps %>% .[! . %in% famods]
  if (purrr::is_empty(ftmod)) {
    ftmod <- NULL
  } else if (ftmod == "") {
    ftmod <- NULL
  }
  
  # fixed N-term, C-term
  fntmod <- ftmod %>% .[. == "N-term"]
  fctmod <- ftmod %>% .[. == "C-term"]
  
  
  # "amods- tmod- vnl- fnl-"
  if (is_empty(fmods_nl)) {
    type <- "fnl-" 
    if (is_empty(vmods_nl)) {
      type <- paste("vnl-", type)
      if (is_empty(tmod)) {
        type <- paste("tmod-", type)
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 1
        } else {
          type <- paste("amods+", type) # 2
        }
      } else {
        type <- paste("tmod+", type)
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 3
        } else {
          type <- paste("amods+", type) # 4
        }
      }
    } else {
      type <- paste("vnl+", type)
      if (is_empty(tmod)) {
        type <- paste("tmod-", type)
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 5
        } else {
          type <- paste("amods+", type) # 6
        }
      } else {
        type <- paste("tmod+", type)
        
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 7
        } else {
          type <- paste("amods+", type) # 8
        }
      }
    }
  } else {
    type <- "fnl+"
    if (is_empty(vmods_nl)) {
      type <- paste("vnl-", type)
      if (is_empty(tmod)) {
        type <- paste("tmod-", type)
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 1
        } else {
          type <- paste("amods+", type) # 2
        }
      } else {
        type <- paste("tmod+", type)
        
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 3
        } else {
          type <- paste("amods+", type) # 4
        }
      }
    } else {
      type <- paste("vnl+", type)
      if (is_empty(tmod)) {
        type <- paste("tmod-", type)
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 5
        } else {
          type <- paste("amods+", type) # 6
        }
      } else {
        type <- paste("tmod+", type)
        
        if (is_empty(amods)) {
          type <- paste("amods-", type) # 7
        } else {
          type <- paste("amods+", type) # 8
        }
      }
    }
  }
  
  attr(aa_masses, "type") <- type
  
  attr(aa_masses, "fmods_nl") <- fmods_nl
  attr(aa_masses, "famods") <- famods
  attr(aa_masses, "ftmod") <- ftmod
  attr(aa_masses, "fntmod") <- fntmod
  attr(aa_masses, "fctmod") <- fctmod
  
  attr(aa_masses, "vmods_nl") <- vmods_nl
  attr(aa_masses, "amods") <- amods
  attr(aa_masses, "tmod") <- tmod
  attr(aa_masses, "ntmod") <- ntmod
  attr(aa_masses, "ctmod") <- ctmod
  
  invisible(aa_masses)
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
  
  if (len == 0) return(NULL)
  
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

  c(res, res_cts)
}


#' Concatenates adjacent peptides in a list (with mass).
#' 
#' @inheritParams concat_peps
#' @examples
#' \donttest{
#' peps <- c(a = 1, b = 2, c = 3)
#' res <- roll_sum(peps, 2)
#' }
roll_sum <- function (peps, n = 2, include_cts = TRUE) {
  len <- length(peps)
  
  if (len == 0) return(NULL)
  
  if (n >= len) n <- len - 1
  
  res <- purrr::map(seq_len((len - n)), ~ {
    ranges <- .x:(.x + n)
    
    nms <- names(peps)[ranges] %>% 
      purrr::accumulate(paste0) 
    
    vals <- cumsum(peps[ranges]) %>% 
      `names<-`(nms)
  }) %>% 
    unlist()
  
  if (include_cts && n >= 1) {
    res_cts <- local({
      cts <- peps[(len - n + 1):len]
      
      purrr::map(n:1, ~ {
        x <- tail(cts, .x)
        
        nms <- names(x) %>% 
          purrr::accumulate(paste0)
        
        vals <- cumsum(x) %>% 
          `names<-`(nms)
      }) %>% unlist()
    })
  } else {
    res_cts <- NULL
  }
  
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
#' 
#' # ok with spaces in the 'title' 
#' x <- parse_unimod("Met-loss (Protein N-term = M)")
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

  # (assumed) no space in `title`
  # title <- unimod %>% gsub("(.*)\\s\\([^\\(]*\\)$", "\\1", .)
  title <- unimod %>% gsub("^([^ ]+?) .*", "\\1", .)

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
           "For example, use 'Acetyl' (title) instead of 'Acetylation' (full_name).", 
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
  
  x[1:min(n + 1, len)]
} 


#' Exclude mis-cleavages in a vector.
#'
#' @inheritParams keep_n_misses
#' @seealso keep_n_misses
exclude_n_misses <- function (x, n) {
  len <- length(x)
  
  stopifnot(n >= 0L, len >= 1L)

  x[-(1:min(n + 1, len))]
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
  n <- min(length(x), n)
  x[seq_len(n)] <- gsub(char, "", x[seq_len(n)])
  
  x
}


#' Remove a trailing character from the last \code{n} entries.
#' 
#' @param char A trailing character to be removed.
#' @inheritParams rm_char_in_nfirst
rm_char_in_nlast <- function (x, char = "-$", n = (max_miss + 1) * 2) {
  len <- length(x)
  n <- min(len, n)
  x[(len-n+1):len] <- gsub(char, "", x[(len-n+1):len])
  
  x
}


#' Make peptide sequences from fastas.
#' 
#' @param fasta_db Fasta database.
#' @inheritParams calc_pepmasses
make_fastapeps <- function (fasta_db, max_miss = 2, min_len = 1, 
                            max_len = 100) {
  
  inds_m <- grep("^M", fasta_db)
  
  fasta_db <- fasta_db %>% 
    purrr::map(~ gsub("([KR]{1})", paste0("\\1", "@"), .x) %>% 
                 paste0("-", ., "-")) 
  
  fasta_dbm <- fasta_db[inds_m] %>% 
    purrr::map(~ gsub("^-M", "-", .x))
  
  # --- Protein N-term initiator methionine kept ---
  peps <- fasta_db %>% 
    purrr::map(~ .x %>% stringr::str_split("@", simplify = TRUE)) %>% 
    purrr::map(concat_peps, max_miss) 
  
  if (min_len > 1 && !is.infinite(max_len)) {
    peps <- peps %>% 
      purrr::map(
        ~ .x %>% .[str_exclude_count(.) >= min_len & str_exclude_count(.) <= max_len]) 
  }
  
  # --- Protein N-term initiator methionine removed ---
  peps_m <- fasta_dbm %>% 
    purrr::map(~ .x %>% 
                 stringr::str_split("@", simplify = TRUE) %>% 
                 keep_n_misses(max_miss)) %>% 
    purrr::map(~ concat_peps(.x, max_miss, include_cts = FALSE)) 
  
  if (min_len > 1 && !is.infinite(max_len)) {
    peps_m <- peps_m %>% 
      purrr::map(
        ~ .x %>% .[str_exclude_count(.) >= min_len & str_exclude_count(.) <= max_len]) 
  }
  
  peps[inds_m] <- purrr::map2(peps_m, peps[inds_m], `c`)
  
  rm(inds_m, fasta_db, fasta_dbm, peps_m)
  
  invisible(peps)
}


#' Generates and Calculates the masses of tryptic peptides from a fasta
#' database.
#'
#' @param index_mods Logical; if TRUE, converts the names of modifications to
#'   indexes. Not currently used.
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
#' @param digits Integer; the number of decimal places to be used.
#' @inheritParams splitPSM
#' @inheritParams binTheoPeps
#' @inheritParams calc_aamasses
#' @inheritParams mcalc_monopep
#' @inheritParams calc_monopep
#' @inheritParams load_fasta2
#' @import parallel
#' @importFrom rlang current_env
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
#' unlist(res_data[[1]], use.names = FALSE) %>% length()
#'
#' # fixedmods, fixedmods + fixedmods_neulosses, varmods, varmods_neulosses
#' unlist(res_data, use.names = FALSE) %>% length()
#'
#' }
#'
#' @export
calc_pepmasses <- function (
  fasta = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
  acc_type = "uniprot_acc", 
  acc_pattern = NULL, 
  fixedmods = c("TMT6plex (K)", 
                "Carbamidomethyl (C)"), 
  varmods = c("TMT6plex (N-term)", 
              "Acetyl (Protein N-term)", 
              "Oxidation (M)", 
              "Deamidated (N)", 
              "Gln->pyro-Glu (N-term = Q)"), 
  include_insource_nl = FALSE, 
  index_mods = FALSE, 
  enzyme = c("trypsin"), 
  maxn_fasta_seqs = 50000,
  maxn_vmods_setscombi = 64, 
  maxn_vmods_per_pep = 5,
  maxn_sites_per_vmod = 3, 
  min_len = 7, max_len = 100, max_miss = 2, 
  out_path = "~/proteoQ/outs", 
  digits = 5, 
  parallel = TRUE) {

  old_opts <- options()
  options(warn = 1)
  on.exit(options(old_opts), add = TRUE)
  
  on.exit(
    if (exists(".savecall", envir = current_env())) {
      if (.savecall) {
        save_call2(path = file.path(.path_cache), fun = "calc_pepmasses", 
                   time = .time_stamp)
      }
    } else {
      # warning("terminated abnormally.", call. = TRUE)
    }, 
    add = TRUE
  )
  
  stopifnot(vapply(c(min_len, max_len, max_miss), is.numeric, 
                   logical(1)))
  
  stopifnot(min_len >= 0, max_len >= min_len, max_miss <= 100)
  
  # ---
  .path_cache <- create_dir("~/proteoQ/.MSearch/Cache/Calls")
  .time_stamp <- match_calltime(path = .path_cache, 
                                fun = "calc_pepmasses", 
                                nms = c("fasta", "acc_type", "acc_pattern", 
                                        "fixedmods", "varmods", 
                                        "include_insource_nl", "enzyme", 
                                        "maxn_fasta_seqs", "maxn_vmods_setscombi", 
                                        "maxn_vmods_per_pep", "maxn_sites_per_vmod", 
                                        # "maxn_vmods_sitescombi_per_pep", 
                                        "min_len", "max_len", "max_miss"))
  
  .path_fasta <- fasta %>% 
    gsub("\\\\", "/", .) %>% 
    gsub("(.*/)[^/]*", "\\1", .) %>% 
    map_chr(find_dir) %>% 
    `[[`(1)
  
  # ---
  if (!is_empty(.time_stamp)) {
    message("Loading peptide masses from cache.")
    
    pep_masses <- readRDS(file.path(.path_fasta, "pepmasses", .time_stamp, 
                                    "pepmasses.rds"))

    .savecall <- FALSE
  } else {
    delete_files(out_path, ignores = c("\\.[Rr]$", "\\.(mgf|MGF)$", "\\.xlsx$", 
                                       "\\.xls$", "\\.csv$", "\\.txt$", 
                                       "^mgf$", "^mgfs$"))

    .time_stamp <- format(Sys.time(), ".%Y-%m-%d_%H%M%S")
    
    # ---
    aa_masses <- local({
      message("Computing the combinations of fixed and variable modifications.")
      
      if (index_mods) {
        mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
          as.hexmode() %>% 
          `names<-`(c(fixedmods, varmods))
      } else {
        mod_indexes <- NULL
      }
      
      aa_masses <- calc_aamasses(fixedmods = fixedmods, 
                                 varmods = varmods, 
                                 maxn_vmods_setscombi = maxn_vmods_setscombi, 
                                 mod_indexes = mod_indexes)
    })

    # ---
    message("Preparing peptide sequences.")
    
    pre_masses <- pre_pepmasses(fasta = fasta, 
                                acc_type = acc_type,
                                acc_pattern = acc_pattern, 
                                maxn_fasta_seqs = maxn_fasta_seqs, 
                                aa_masses = aa_masses, 
                                max_miss = max_miss, 
                                min_len = min_len, 
                                max_len = max_len)

    # ---
    message("Calculating peptide masses (target) ...")
    
    pep_masses <- purrr::map2(pre_masses$fwds, 
                              aa_masses, 
                              mcalc_monopep, 
                              include_insource_nl = include_insource_nl, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              parallel = parallel, 
                              digits = digits)
    
    message("Calculating peptide masses (decoy) ...")
    
    rev_masses <- purrr::map2(pre_masses$revs, 
                              aa_masses, 
                              mcalc_monopep, 
                              include_insource_nl = include_insource_nl, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              parallel = parallel, 
                              digits = digits)
    
    rm(pre_masses)
    
    # ---
    out_path <- create_dir(file.path(.path_fasta, "pepmasses", .time_stamp))
    saveRDS(pep_masses, file.path(out_path, "pepmasses.rds"))
    saveRDS(rev_masses, file.path(out_path, "pepmasses_rev.rds"))

    # ---
    .savecall <- TRUE
  }
  
  assign(".path_cache", .path_cache, envir = .GlobalEnv)
  assign(".time_stamp", .time_stamp, envir = .GlobalEnv)
  assign(".path_fasta", .path_fasta, envir = .GlobalEnv)

  invisible(pep_masses)
}


#' Helper of \link{calc_pepmasses}.
#' 
#' Prior to the calculation of peptide masses.
#' 
#' @inheritParams calc_pepmasses
#' @inheritParams mcalc_monopep
pre_pepmasses <- function (fasta, acc_type, acc_pattern, maxn_fasta_seqs, aa_masses, 
                           max_miss, min_len, max_len) {
  # ---
  message("Loading fasta.")
  
  fasta_db <- load_fasta2(fasta, acc_type, acc_pattern) 
  
  if (length(fasta_db) > maxn_fasta_seqs) {
    stop("More than `", maxn_fasta_seqs, "` sequences in fasta files.\n",
         "  May consider a higher `maxn_fasta_seqs`.", 
         call. = FALSE)
  }

  n_cores <- detectCores()
  cl <- makeCluster(getOption("cl.cores", n_cores))
  
  clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  clusterExport(cl, list("concat_peps"), envir = environment(proteoQ:::concat_peps))
  
  fasta_db <- chunksplit(fasta_db, n_cores)

  rev_fasta_db <- clusterApply(cl, fasta_db, stringi::stri_reverse) %>% 
    purrr::map2(., fasta_db, ~ {
      names(.x) <- paste0("-", names(.y))
      .x
    })
  
  # ---
  message("Generating peptide sequences (target and decoy).")
  
  peps <- clusterApply(cl, fasta_db,  make_fastapeps, 
                       max_miss, min_len, max_len) %>% 
    purrr::flatten()
  
  rev_peps <- clusterApply(cl, rev_fasta_db, make_fastapeps, 
                           max_miss, min_len, max_len) %>% 
    purrr::flatten()
  
  rm(fasta_db, rev_fasta_db)
  stopCluster(cl)

  # ---
  message("Distributing peptides by fixed and variable modifications.")
  
  fwds <- aa_masses %>% 
    purrr::map(subpeps_by_vmods, peps) %>% 
    purrr::map(~ {
      xs <- .x
      
      xs %>% 
        purrr::map(rm_char_in_nfirst, char = "^-", n = (max_miss + 1) * 2) %>% 
        purrr::map(rm_char_in_nlast, char = "-$", n = (max_miss + 1) * 2)
    })
  
  revs <- aa_masses %>% 
    purrr::map(subpeps_by_vmods, rev_peps) %>% 
    purrr::map(~ {
      xs <- .x
      
      xs %>% 
        purrr::map(rm_char_in_nfirst, char = "^-", n = (max_miss + 1) * 2) %>% 
        purrr::map(rm_char_in_nlast, char = "-$", n = (max_miss + 1) * 2)
    })
  
  rm(peps, rev_peps)
  
  invisible(list(fwds = fwds, revs = revs))
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
  
  if (!purrr::is_empty(mods)) {
    all_mods <- paste(mods, collapse = ", ")
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
      site <- find_aa_site(.x)
      m <- aa_masses[site]
      aa_masses[site] <<- m + .y
    })
  } else {
    #  Add mod_masses of variable mods (multiple lists)
    aas <- purrr::map2(positions_sites, mod_masses, ~ {
      site <- find_aa_site(.x)
      aa_masses[site] <- aa_masses[site] + .y
      aa_masses
    }, aa_masses)
    
    # Flatten the lists (with attributes being kept)
    aa_masses <- local({
      masses <- purrr::map2_dbl(positions_sites, aas, ~ .y[.x])
      attrs <- attributes(aa_masses)
      aa_masses <- c(aa_masses, masses)
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
  
  attr(aa_masses, paste0(mod_type, "_neuloss")) <- neulosses

  list(aa_masses)
}


