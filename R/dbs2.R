#' Reverses peptide sequences.
#' 
#' @param xs Lists of proteins with named peptide masses under each entry.
rev_pepseqs <- function (xs) {
  
  names(xs) <- paste0("-", names(xs))
  
  rev_seqs <- purrr::map(xs, ~ {
    seqs <- names(.x)
    a <- stringi::stri_sub(seqs, 1, 1)
    b <- stringi::stri_sub(seqs, -1, -1)
    
    lens <- stri_length(seqs)
    
    revs <- stringi::stri_reverse(seqs)
    substring(revs, 1) <- a
    substring(revs, lens) <- b
    
    revs
  })
  
  out <- purrr::map2(xs, rev_seqs, ~ {
    names(.x) <- .y
    .x
  })
}


#' Remove a starting character from the first \code{n} entries.
#' 
#' @param x A list of character strings.
#' @param char A starting character to be removed.
#' @param n The number of beginning entries to be considered.
rm_char_in_nfirst2 <- function (x, char = "^-", n = (max_miss + 1L) * 2L) {
  nms <- names(x)
  
  len <- length(nms)
  n <- min(len, n)
  
  seqs <- seq_len(n)
  
  nms[seqs] <- gsub(char, "", nms[seqs])
  names(x) <- nms
  
  x
}


#' Remove a trailing character from the last \code{n} entries.
#' 
#' @param char A trailing character to be removed.
#' @inheritParams rm_char_in_nfirst
rm_char_in_nlast2 <- function (x, char = "-$", n = (max_miss + 1L) * 2L) {
  nms <- names(x)
  
  len <- length(nms)
  n <- min(len, n)
  
  seqs <- (len-n+1):len
  
  nms[seqs] <- gsub(char, "", nms[seqs])
  names(x) <- nms

  x
}


#' Converts protein accessions and peptide sequences to integers.
#'
#' @param pep_seqs Lists of peptide sequences by protein accessions. Peptide
#'   sequences are names and monoisotopic masses in values.
#' @param out_path An output path
#' @param type An indicator of either forward or reversed sequences.
#' @examples
#' \donttest{
#' pep_seqs <- list(ALG2_HUMAN = c("-MAEEQGR" = 819.35450, 
#'                                 "-MAEEQGRER" = 1122.50877 ))
#' 
#' out <- conv_to_pepids(pep_seqs, out_path, "rev")
#' }
conv_to_pepids <- function (data, out_path, type = "fwd") {
  
  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
  
  # pep_seqs <- purrr::map(data, attr, "data")
  
  # two-stage matches
  #   prot_acc -> pep_seq(s); find position indexes
  
  lens <- purrr::map_int(pep_seqs, length)
  
  prots <- purrr::imap(lens, ~ rep(.y, .x)) %>% 
    unlist(use.names = FALSE)
  
  peps <- pep_seqs %>% 
    purrr::map(names) %>% 
    unlist(use.names = FALSE)
  
  df_protids <- local({
    prs <- names(pep_seqs)
    data.frame(prot_acc = prs, prot_index = seq_along(prs))
  })
  
  df_pepids <- local({
    pes <- unique(peps)
    data.frame(pep_seq = pes, pep_index = seq_along(pes))
  })
  
  df <- data.frame(prot_acc = prots, pep_seq = peps) %>% 
    dplyr::left_join(df_protids, by = "prot_acc") %>% 
    dplyr::left_join(df_pepids, by = "pep_seq") %>% 
    dplyr::group_by(pep_seq) %>% 
    dplyr::mutate(n = row_number()) %T>% 
    saveRDS(file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  stopifnot(identical(unique(df$prot_acc), names(pep_seqs)))
  
  dfs <- df %>% split(.$prot_index)
  # saveRDS(dfs, file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  names(pep_seqs) <- names(dfs)
  
  pep_seqs <- purrr::map2(pep_seqs, dfs, ~ {
    names(.x) <- .y[["pep_index"]]
    .x
  })
}








conv_to_pepids_v1 <- function (pep_seqs, out_path, type = "fwd") {
  
  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
  
  # pep_seqs <- purrr::map(pep_seqs, attr, "data")
  
  lens <- purrr::map_int(pep_seqs, length)
  
  prots <- purrr::imap(lens, ~ rep(.y, .x)) %>% 
    unlist(use.names = FALSE)
  
  peps <- pep_seqs %>% 
    purrr::map(names) %>% 
    unlist(use.names = FALSE)
  
  df_protids <- local({
    prs <- names(pep_seqs)
    data.frame(prot_acc = prs, prot_index = seq_along(prs))
  })
  
  df_pepids <- local({
    pes <- unique(peps)
    data.frame(pep_seq = pes, pep_index = seq_along(pes))
  })
  
  df <- data.frame(prot_acc = prots, pep_seq = peps) %>% 
    dplyr::left_join(df_protids, by = "prot_acc") %>% 
    dplyr::left_join(df_pepids, by = "pep_seq") %>% 
    dplyr::group_by(pep_seq) %>% 
    dplyr::mutate(n = row_number()) %T>% 
    saveRDS(file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  stopifnot(identical(unique(df$prot_acc), names(pep_seqs)))
  
  dfs <- df %>% split(.$prot_index)
  # saveRDS(dfs, file.path(out_path, "temp", paste0("pep_indexes_", type, ".rds")))
  
  names(pep_seqs) <- names(dfs)
  
  pep_seqs <- purrr::map2(pep_seqs, dfs, ~ {
    names(.x) <- .y[["pep_index"]]
    .x
  })
}



