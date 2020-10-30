#' Reads fasta
#'
#' Reads a fasta file.
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
#' fasta <- read_fasta(file = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' head(names(fasta))
#' 
#' # use the first fifty characters
#' fasta <- read_fasta(file = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
#'                     pattern = ">(.{50}).*")
#' head(names(fasta))
#' 
#' # uniprot_acc
#' fasta <- read_fasta(file = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
#'                     pattern = ">..\\|([^\\|]+)\\|.*")
#' head(names(fasta))
#' 
#' # use every in the header
#' fasta <- read_fasta(file = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
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
#' @param fasta A list of protein entries from \code{\link{read_fasta}}.
#' @inheritParams read_fasta
#' @examples
#' \dontrun{
#' fasta <- read_fasta(file = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
#' write_fasta(fasta, "~/proteoQ/examples/my.fasta")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
#' @export
write_fasta <- function (fasta, file) {
  purrr::map(fasta, ~ paste(attr(.x, "header"), .x, sep = "\n")) %>% 
    unlist() %>% 
    writeLines(file)
}


#' Peptide molecular weight
#'
#' Calculates the mono-isotopic molecular weight of a polypeptide (neutral
#' species).
#' 
#' @param aa_seq Character string; the sequence of a peptide with one-letter
#'   representation of amino acids.
#' @param digits Integer; the number of decimal places to be used.
#' @examples
#' \dontrun{
#' calc_mono_pep("AAIDWFDGKEFSGNPIK")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
calc_mono_pep <- function (aa_seq, digits = 4) {
  # tmt6_mass <- 229.162932
  # tmtpro_mass <- 304.207146
  # h2o <- 18.010565
  # proton <- 1.00727647
  # hydrogen <- 1.007825
  # oxygen <- 15.99491462
  
  # data(aa_residues)
  # aa_masses <- aa_residues$monoisotopic_da %>% `names<-`(aa_residues$one_letter)

  aa_masses <- c(
    A = 71.037114, R = 156.101111, N = 114.042927, D = 115.026943, 
    C = 103.009185, E = 129.042593, Q = 128.058578, G = 57.021464, 
    H = 137.058912, I = 113.084064, L = 113.084064, K = 128.094963,
    M = 131.040485, F = 147.068414, P = 97.052764, S = 87.032028,
    T = 101.047679, W = 186.079313, Y = 163.063329, V = 99.068414, 
    "N-term" = 1.007825, "C-term" = 17.002740, 
    U = 150.953633, B = 114.534940, X = 111.000000, Z = 128.550590)
  
  aa_seq %>% 
    stringr::str_split("", simplify = TRUE) %>% 
    aa_masses[.] %>% 
    sum() %>% 
    `+`(18.010565) %>% 
    round(digits = digits) 
}


#' Peptide molecular weight
#'
#' Calculates the average molecular weight of a polypeptide (neutral species).
#' 
#' @inheritParams calc_mono_pep
#' @examples
#' \dontrun{
#' calc_avg_pep("AAIDWFDGKEFSGNPIK")
#' }
#'
#' @import dplyr purrr
#' @importFrom magrittr %>% %T>% %$% %<>%
calc_avg_pep <- function (aa_seq, digits = 4) {
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


