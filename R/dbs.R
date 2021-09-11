#' Reads a file in fasta format
#'
#' Reads a file in fasta format by line (used in both
#' \link[proteoQ]{load_fasta} and \link[proteoM]{load_fasta2}).
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
#' Writes a fasta file (Not yet used).
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
#' @param aa_seq Character string; a peptide sequences with one-letter
#'   representation of amino acids.
#' @param digits Integer; the number of decimal places to be used.
#' 
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


