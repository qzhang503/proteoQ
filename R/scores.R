#' Calculates peptide scores by single MGF entries (scan numbers).
#' 
#' Each entry corresponds to a row.
#' 
#' @param entry A row of data from \link{pmatch_bymgfs}.
#' @inheritParams matchMS
scalc_pepscores <- function (entry, topn_ms2ions = 100) {
  mts <- entry$matches[[1]]
  
  # N - the total number of features (white and black balls)
  # k - the number of sampled features
  # m(s) - the numbers of theoretical features (white balls)
  # n(s) - the number of noise (black balls)
  
  N <- entry$ms2_n[[1]]
  k <- min(topn_ms2ions, N)
  
  ms <- stringi::stri_length(names(mts)) * 2
  ms[ms > N] <- N
  ns <- N - ms

  xs <- purrr::map(mts, ~ {
    pep <- .x 
    
    if (is.data.frame(pep)) {
      nrow(pep)
    } else {
      # with multiple mod positions
      purrr::map(pep, nrow)
    }
  })
  
  lens <- purrr::map(xs, length)
  ms <- purrr::map2(ms, lens, rep)
  ns <- purrr::map2(ns, lens, rep)

  scores <- mapply(mdhyper, x = xs, m = ms, n = ns, k  = k, SIMPLIFY = FALSE) %>% 
    purrr::map(~ -log10(.x) * 10)
  
  scores <- scores %>% 
    purrr::map(~ {
      .x[.x > 250] <- 250
      .x
    })
}


#' The helper of \link[stats]{dhyper}.
#' 
#' @inheritParams stats::dhyper
mdhyper <- function (x, m, n, k) {
  mapply(dhyper, x, m, n, k)
}


#' Flattens the outputs of ion matches.
#' 
#' @param data The data from ion matches.
#' @param outcol The output column name.
#' @examples
#' \donttest{
#' x <- flatten_pepouts(res$matches, "matches")
#' scores <- flatten_pepouts(scores, "pep_score")
#' }
flatten_pepouts <- function (data, outcol = "matches") {
  purrr::map(data, ~ {
    tib <- .x
    
    tib <- purrr::map(tib, ~ {
      tibble(pep_mod = names(.x), !!outcol := .x)
    }) 
    
    nm <- names(tib)
    
    tib <- purrr::map2(tib, nm, ~ {
      .x$pep_seq <- .y
      .x
    }) %>% 
      dplyr::bind_rows()
  }) 
}


#' Calculates the scores of peptides at an \code{aa_masses}.
#' 
#' @inheritParams calc_pepscores
calc_pepscores_i <- function (res, topn_ms2ions = 100, 
                              out_path = "~/proteoQ/outs") {

  scores <- res %>% 
    dplyr::select(c("matches", "ms2_n")) %>% 
    split(., seq_len(nrow(.))) %>% 
    purrr::map(scalc_pepscores, topn_ms2ions) %>% 
    flatten_pepouts("pep_score")

  # with list columns
  res <- res %>% 
    dplyr::mutate(pep_score = scores) %T>% 
    saveRDS(file.path(out_path, "scores.rds")) %>% 
    dplyr::select(-c("ms2_moverz", "ms2_int", "matches", "pep_score"))

  # without list columns
  scores <- scores %>% 
    purrr::map2(res$scan_num, ~ {
    dplyr::bind_cols(scan_num = .y, .x)
  }) %>% 
    dplyr::bind_rows() 
  
  if (nrow(scores) > 0) {
    scores <- scores %>% 
      dplyr::group_by(scan_num) %>% 
      dplyr::arrange(-pep_score) %>% 
      dplyr::mutate(pep_rank = row_number()) %>% 
      dplyr::mutate(pep_score = round(pep_score, digits = 1)) %>% 
      dplyr::ungroup() 
    
    res <- res %>% 
      purrr::map(unlist, use.names = FALSE) %>% 
      dplyr::bind_cols() %>% 
      dplyr::left_join(scores, by = "scan_num") 
  } else {
    scores <- tibble(pep_mod = character(), 
                     pep_score = numeric(), 
                     pep_seq = character(), 
                     pep_rank = integer())

    res <- res %>% dplyr::bind_cols(scores)
  }

  invisible(res)
}


#' Calculates the scores of peptides.
#' 
#' @param res The result from \link{pmatch_bymgfs}.
#' @inheritParams matchMS
calc_pepscores <- function (res = NULL, topn_ms2ions = 100, 
                            out_path = "~/proteoQ/outs") {
  
  if (is.null(res)) {
    res <- readRDS(file.path(out_path, "ion_matches.rds"))
  }

  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  out <- parallel::clusterMap(cl, calc_pepscores_i, 
                              res = res, 
                              MoreArgs = list(topn_ms2ions = topn_ms2ions, 
                                              out_path = out_path), 
                              .scheduling = "dynamic")
  
  parallel::stopCluster(cl)

  oks <- purrr::map_dbl(out, ~ nrow(.x) > 0)
  out <- out[oks] %>% dplyr::bind_rows()
}

