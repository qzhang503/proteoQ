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
#' \dontrun{
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
#' \dontrun{
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
#' \dontrun{
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


#' Calculates molecular weight of a polypeptide ([MH]+).
#'
#' @param aa_seq Character string; the sequence of a peptide with one-letter
#'   representation of amino acids.
#' @param digits Integer; the number of decimal places to be used. Not currently
#'   used.
#' @examples
#' \dontrun{
#' x <- calc_monopep("AAIDWFDGKEFSGNPIK")
#' }
calc_monopep <- function (aa_seq, digits = 5) {
  options(digits=9)
  
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # oxygen <- 15.99491462
  
  # moverz <- (mass + h2o + z*proton)/z
  

  
  run_scripts <- FALSE
  if (run_scripts) {
    data(package = "proteoQ", aa_residues)
    aa_masses <- aa_residues %>% 
      dplyr::select(c("one_letter", "monoisotopic_da")) 
    aa_masses <- aa_masses$monoisotopic_da %>% `names<-`(aa_masses$one_letter)
  }
  
  aa_masses <- c(
    A = 71.037114, R = 156.101111, N = 114.042927, D = 115.026943, 
    C = 103.009185, E = 129.042593, Q = 128.058578, G = 57.021464, 
    H = 137.058912, I = 113.084064, L = 113.084064, K = 128.094963,
    M = 131.040485, F = 147.068414, P = 97.052764, S = 87.032028,
    T = 101.047679, W = 186.079313, Y = 163.063329, V = 99.068414, 
    "N-term" = 1.007825, "C-term" = 17.002740, 
    U = 150.953633, B = 114.534940, X = 111.000000, Z = 128.550590)
  
  aa_masses <- c(aa_masses, "-" = 0)
  
  # ---
  run_scripts <- FALSE
  if (run_scripts) {
    x <- find_modmass("Carbamidomethyl")
  }

  # ---
  aa_seq %>% 
    stringr::str_split("", simplify = TRUE) %>% 
    aa_masses[.] %>% 
    sum() %>% 
    `+`(19.01784) %>% # H3O+
    setNames(aa_seq) %>% 
    round(digits = digits) 
}


#' Averaged peptide molecular weight (for proteins)
#'
#' Calculates the average molecular weight of a polypeptide (neutral species).
#' 
#' @inheritParams calc_monopep
#' @examples
#' \dontrun{
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


#' Concatenates two adjacent peptides in a list.
#'
#' @param peps A list of peptide sequences with a one-letter representation of
#'   amino acid residues.
#' @param n The number of miss cleavages.
#'
#' @examples
#' \dontrun{
#' peps <- c("MESNHK", "SGDGLSGTQK", "EAALR", "ALVQR", "TGYSLVQENGQR")
#' concat_peps(peps)
#' }
concat_peps <- function (peps, n = 1) {
  len <- length(peps)
  
  if (n >= len) n <- len - 1
  
  purrr::map(seq_len((len - n)), ~ {
    peps[.x:(.x + n)] %>% purrr::accumulate(paste0) 
  }) %>% unlist()
}


#' Find the mono-isotopic mass of a modification
#' 
#' @param mod_nm The name of a modification.
find_modmass <- function (mod_nm = "Carbamidomethyl") {
  options(digits=9)
  
  file <- system.file("extdata", "master.xml", package = "proteoQ")
  parent <- xml2::read_xml(file)
  
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
  
  mod <- local({
    idx <- which(xml2::xml_attr(modifications, "title") == mod_nm)
    
    if (purrr::is_empty(idx)) {
      stop("Modification not found: '", mod_nm, "'.\n", 
           "For example, use 'Acetyl' instead of 'Acetylation'.", 
           call. = FALSE)
    }
    
    mod <- modifications[[idx]]
  })
  
  stopifnot(xml2::xml_attrs(mod) %>% .["title"] == mod_nm)
  
  # --- children of this `mod` ---
  modch <- xml2::xml_children(mod)
  
  monomass <- xml2::xml_attrs(modch) %>% 
    purrr::map(`[`, "mono_mass") %>%
    `[`(!is.na(.)) %>% 
    unlist() %>% 
    as.numeric()
  
  stopifnot(length(monomass) == 1)
  
  # --- find sites and positions ---
  sites_positions <- xml2::xml_attrs(modch) %>% 
    map(`[`, c("site", "position")) %>% 
    purrr::map(~ invisible(setNames(.x[1], .x[2]))) 
  
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
    
    sites_positions[idx_nl] <- 
      purrr::map2(idx_nl, nls, ~ {c(sites_positions[.x], .y)}) 
    
    sites_positions <- sites_positions %>% 
      .[!is.na(.)]
  }

  invisible(list(mass = monomass, mods = sites_positions))
}


#' Find the first n elements in a vector
#'
#' @param x A vectot of data.
#' @param n Integer. The first n elements by index in \code{x} to be kept.
keep_one_to_n <- function (x, n) {
  x[1:pmin(n+1, length(x))]
} 


#' Calculates the masses of tryptic peptides from a fasta database.
#'
#' @param min_len The minimum length of peptides. Shorter peptides will be
#'   excluded.
#' @param max_miss The maximum length of peptides. Longer peptides will be
#'   excluded.
#' @param out_path The file path to outputs.
#' @param out_name The file name of outputs.
#' @inheritParams splitPSM
#' @examples
#' \dontrun{
#' res <- calc_pepmasses()
#' }
calc_pepmasses <- function (fasta = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
                            out_path = "~/proteoQ/dbs/fasta/uniprot/pepmasses/pepmasses.rds", 
                            min_len = 7, max_len = 100, max_miss = 2) {

  stopifnot(vapply(c(min_len, max_len, max_miss), is.numeric, 
                   logical(1)))
  
  stopifnot(min_len >= 0, max_len >= min_len, max_miss <= 100)
  
  fasta_db <- fasta %>% load_fasta() 
  inds_m <- grep("^M", fasta_db)

  fasta_db <- fasta_db %>% 
    purrr::map(~ gsub("([KR]{1})", paste0("\\1", "@"), .x) %>% paste0("-", ., "-")) 
  
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

  # --- methionine kept ---
  peps <- fasta_db %>% 
    purrr::map(~ .x %>% stringr::str_split("@", simplify = TRUE)) %>% 
    purrr::map(concat_peps, max_miss) %>% # a slow step
    purrr::map(~ .x %>% 
                 .[stringr::str_count(.) >= min_len & stringr::str_count(.) <= max_len]) 
  
  peps <- peps %>% 
    purrr::map(~ {
      .x %>% 
        purrr::map(calc_monopep) %>% 
        unlist() 
    })
  
  # --- methionine removed ---
  peps_m <- fasta_dbm %>% 
    purrr::map(~ .x %>% 
                 stringr::str_split("@", simplify = TRUE) %>% 
                 keep_one_to_n(max_miss)) %>% 
    purrr::map(~ concat_peps(.x, max_miss)) %>% 
    purrr::map(~ .x %>% 
                 .[stringr::str_count(.) >= min_len & stringr::str_count(.) <= max_len]) %>% 
    purrr::map(~ {
      .x %>% 
        purrr::map(calc_monopep) %>% 
        unlist() 
    }) 
  
  # --- combined ---
  peps[inds_m] <- purrr::map2(peps_m, peps[inds_m], `c`)

  # peptides: 
  # 1636749 at 1 misses
  # 2831047 at 2 misses
  # 4027633 at 3 misses
  # 5178789 at 4 misses
  
  dir.create(gsub("(^.*/).*$", "\\1", out_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(peps, file.path(out_path))

  invisible(peps)
}


#' Add the masses of fixed and/or variable modifications to tryptic peptides.
#' 
#' @param file A character string to the name of a peptide-mass file.
#' @examples
#' \dontrun{
#' pepmasses <- add_modmasses("~/proteoQ/dbs/fasta/uniprot/pepmasses/pepmasses.rds")
#' }
add_modmasses <- function (file = NULL, varmods = c("Acetyl (N-term)", "Oxidation (M)"), 
                           fixmods = c("Carbamidomethyl (C)")) {
  if (is.null(file)) {
    stop("A file of peptide masses is required.", call. = FALSE)
  }
  
  if (length(file) > 1) {
    stop("One 'file' only.", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Missing file: \n", 
         purrr::reduce(file %>% .[!file.exists(.)], paste, sep = "\n"), 
         call. = FALSE)
  }
  
  peps <- readRDS(file)
  
  # need to know if a peptide sequence contains the site
  
}







#' Calculates molecular weight of a polypeptide ([MH]+).
#'
#' @param aa_seq Character string; the sequence of a peptide with one-letter
#'   representation of amino acids.
#' @param digits Integer; the number of decimal places to be used. Not currently
#'   used.
#' @examples
#' \dontrun{
#' calc_monopep("AAIDWFDGKEFSGNPIK")
#' }
calc_monopep_orig <- function (aa_seq, digits = 4) {
  options(digits=9)
  
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # oxygen <- 15.99491462
  
  # moverz <- (mass + h2o + z*proton)/z
  
  run_scripts <- FALSE
  if (run_scripts) {
    data(package = "proteoQ", aa_residues)
    aa_masses <- aa_residues %>% 
      dplyr::select(c("one_letter", "monoisotopic_da")) 
    aa_masses <- aa_masses$monoisotopic_da %>% `names<-`(aa_masses$one_letter)
  }
  
  aa_masses <- c(
    A = 71.037114, R = 156.101111, N = 114.042927, D = 115.026943, 
    C = 103.009185, E = 129.042593, Q = 128.058578, G = 57.021464, 
    H = 137.058912, I = 113.084064, L = 113.084064, K = 128.094963,
    M = 131.040485, F = 147.068414, P = 97.052764, S = 87.032028,
    T = 101.047679, W = 186.079313, Y = 163.063329, V = 99.068414, 
    "N-term" = 1.007825, "C-term" = 17.002740, 
    U = 150.953633, B = 114.534940, X = 111.000000, Z = 128.550590)
  
  aa_masses <- c(aa_masses, "-" = 0)
  
  aa_seq %>% 
    stringr::str_split("", simplify = TRUE) %>% 
    aa_masses[.] %>% 
    sum() %>% 
    `+`(19.01784) %>% # H3O+
    setNames(aa_seq) %>% 
    round(digits = digits) 
}


