#' Extracts the lists of theoretical m-over-z's.
#'
#' @param data Nested list of table from ion matches. Each table contains
#'   columns \code{theo} and \code{expt}.
#' @import purrr
extract_theos <- function (data) {
  purrr::map(data, ~ {
    x <- .x %>% # by sequence
      purrr::map(`[[`, "theo") # by varmod positions
    
    x <- map(x, ~ {
      names(.x) <- NULL
      .x
    })
  }) 
}


#' The helper of \link[stats]{dhyper}.
#' 
#' @inheritParams stats::dhyper
#' @importFrom purrr map map2
mdhyper <- function (x, m, n, k) {
  k <- map2(m, n, ~ min(.x + .y, k)) 
  mapply(dhyper, x, m, n, k)
}


#' Flattens the outputs of ion matches.
#' 
#' @param outcol The output column name.
#' @inheritParams extract_theos
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


#' Matches against secondary ions
#' 
#' @param theos A list of theoretical values.
#' @param expts A list of experimental values
#' @inheritParams calc_ms2ionseries
#' @inheritParams matchMS
match_secions <- function (theos, expts, type_ms2ions = "by", ppm_ms2 = 25, 
                           digits = 5) {
  purrr::map(theos, ~ {
    theos_i <- .x # by peptide
    
    purrr::map(theos_i, ~ { # by varmods
      .x %>% 
        add_seions(type_ms2ions, digits) %>% 
        find_ppm_outer_bycombi(expts, ppm_ms2)
    }) 
  }) 
}


#' Adds secondary ions of b0, y0 etc.
#' 
#' @param ms2s A vector of theoretical MS2 m-over-z values. 
#' @inheritParams match_secions
add_seions <- function (ms2s, type_ms2ions = "by", digits = 5) {
  len <- length(ms2s)
  
  if (type_ms2ions == "by") {
    proton <- 1.00727647
    h2o <- 18.010565
    nh3 <- 17.026549
    
    bs <- ms2s[1:(len/2)]
    ys <- ms2s[(len/2+1):len]
    
    b2s <- (bs + proton)/2
    bstars <- bs - nh3
    bstar2s <- (bstars + proton)/2
    b0s <- bs - h2o
    b02s <- (b0s + proton)/2
    
    y2s <- (ys + proton)/2
    ystars <- ys - nh3
    ystar2s <- (ystars + proton)/2
    y0s <- ys - h2o
    y02s <- (y0s + proton)/2
    
    round(c(b2s, bstars, bstar2s, b0s, b02s, y2s, ystars, ystar2s, y0s, y02s), 
          digits = digits)
  } else if (type_ms2ions == "ax") {
    proton <- 1.00727647
    h2o <- 18.010565
    nh3 <- 17.026549
    
    as <- ms2s[1:(len/2)]
    xs <- ms2s[(len/2+1):len]
    
    a2s <- (as + proton)/2
    astars <- as - nh3
    astar2s <- (astars + proton)/2
    a0s <- as - h2o
    a02s <- (a0s + proton)/2
    
    x2s <- (xs + proton)/2
    
    round(c(a2s, astars, astar2s, a0s, a02s, x2s), digits = digits)
  } else if (type_ms2ions == "cz") {
    proton <- 1.00727647
    
    cs <- ms2s[1:(len/2)]
    zs <- ms2s[(len/2+1):len]
    
    c2s <- (cs + proton)/2
    z2s <- (zs + proton)/2
    
    round(c(c2s, z2s), digits = digits)
  }
}


#' Calculates peptide scores by single MGF entries (scan numbers).
#' 
#' Each entry corresponds to a row in \code{ion_matches.rds}.
#' 
#' @param entry A row of data from \link{pmatch_bymgfs}.
#' @inheritParams matchMS
#' @import purrr
scalc_pepscores <- function (entry, topn_ms2ions = 100, type_ms2ions = "by", 
                             ppm_ms2 = 25, digits) {
  
  ## only one experimental and thus `[[1]]`
  expts <- entry$ms2_moverz[[1]] 
  
  ## matches btw. theoretical ands experimentals
  
  # [[1]] --- `entry$matches` (always at level-one and can be unlisted)
  # [[1]]$AMMASIGR --- `entry$matches[[1]]` (1:i peptides)
  # [[1]]$AMMASIGR$`00500000` --- `entry$matches[[1]][[1]]` (1:j positions)
  # A tibble: 18 x 2
  # theo  expt
  # <dbl> <dbl>
  #   1  175.   175.
  #   2  230.   230.
  #   3  301.   301.
  # 
  # [[1]]$TNLAMMR
  # [[1]]$TNLAMMR$`0000500`
  # A tibble: 16 x 2
  # theo  expt
  # <dbl> <dbl>
  #   1  175.   175.
  #   2  230.   230.
  #   3  331.   331.
  # 
  # [[1]]$TNLAMMR$`0000050`
  # A tibble: 16 x 2
  # theo  expt
  # <dbl> <dbl>
  #   1  175.   175.
  #   2  230.   230.
  #   3  331.   331.
  
  # flattens by one level
  mts <- entry$matches[[1]]
  
  # lists of theoretical values (primary only)
  theos <- extract_theos(mts)
  
  # $AMMASIGR
  # $AMMASIGR$`00500000`
  # [1]  175.11895  230.17021  301.20732 ...
  
  # $TNLAMMR
  # $TNLAMMR$`0000500`
  # [1]  175.11895  230.17021  331.21789 ...
  
  # $TNLAMMR$`0000050`
  # [1]  175.11895  230.17021  331.21789 ...
  
  # matches additionally against secondary ions
  mts2 <- match_secions(theos, expts, type_ms2ions, ppm_ms2, digits)
  
  # the number of secondary matches
  ms2 <- map(mts2, ~ {
    mts2_i <- .x
    map(mts2_i, ~ sum(!is.na(.x$expt))) 
  })
  
  # the number of primary matches
  xs <- map(mts, ~ {
    mts_i <- .x
    map(mts_i, ~ sum(!is.na(.x$expt))) 
  })
  
  # N - the total number of features (white and black balls)
  # k - the number of sampled features
  # m(s) - the numbers of theoretical features (white balls)
  # n(s) - the number of noise (black balls)
  
  N <- entry$ms2_n[[1]]
  k <- min(topn_ms2ions, N)
  
  ms <- map(mts, ~ {
    mts_i <- .x
    map(mts_i, ~ {
      # x <- nrow(.x) - 2
      x <- nrow(.x) # terminal removed
      x[x > N] <- N
      x
    })
  })
  
  # subtracts the counts of secondary b0, y0 matches etc. from noise
  ns <- map2(ms, ms2, ~ {
    ms_i <- .x
    ms2_i <- .y
    
    map2(ms_i, ms2_i, ~ {
      x <- N - .x - .y
    })
  })
  
  scores <- mapply(mdhyper, x = xs, m = ms, n = ns, k  = k, SIMPLIFY = FALSE) %>% 
    purrr::map(~ -log10(.x) * 10)
  
  scores <- scores %>% 
    purrr::map(~ {
      .x[.x > 250] <- 250
      .x
    })
  
  invisible(list(sec_matches = mts2, scores = scores))
}


#' Calculates the scores of peptides at an \code{aa_masses}.
#' 
#' @param i Integer; the index of the i-th \code{aa_masses} table.
#' @inheritParams calc_pepscores
calc_pepscores_i <- function (res, i, topn_ms2ions = 100, type_ms2ions = "by", 
                              ppm_ms2 = 25, out_path = "~/proteoQ/outs", 
                              digits = 5) {

  # res <- res[62:64, ]
  
  # add scores and secondary matches
  sc_mt2s <- res %>% 
    split(., seq_len(nrow(.))) %>% 
    purrr::map(scalc_pepscores, topn_ms2ions, type_ms2ions, ppm_ms2, digits)

  sec_matches <- sc_mt2s %>% 
    purrr::map(`[[`, "sec_matches") %>% 
    flatten_pepouts("sec_matches")
  
  scores <- sc_mt2s %>% 
    purrr::map(`[[`, "scores") %>% 
    flatten_pepouts("pep_score")
  
  rm(sc_mt2s)

  # outputs with list columns
  res <- res %>% 
    dplyr::mutate(sec_matches = sec_matches, pep_score = scores) %T>% 
    saveRDS(file.path(out_path, paste0("score_", i, ".rds"))) %>% 
    dplyr::select(-c("ms2_moverz", "ms2_int", "matches", 
                     "sec_matches", "pep_score"))

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
#' @param res The result from \link{pmatch_bymgfs}. At the \code{NULL} default,
#'   the file \code{ion_matches.rds} will be used.
#' @inheritParams matchMS
#' @export
calc_pepscores <- function (res = NULL, topn_ms2ions = 100, type_ms2ions = "by", 
                            ppm_ms2 = 25, out_path = "~/proteoQ/outs", 
                            digits = 5) {

  if (is.null(res)) {
    res <- readRDS(file.path(out_path, "ion_matches.rds"))
  }
  
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  out <- parallel::clusterMap(cl, calc_pepscores_i, 
                              res = res, seq_along(res), 
                              MoreArgs = list(topn_ms2ions = topn_ms2ions, 
                                              type_ms2ions = type_ms2ions, 
                                              ppm_ms2 = ppm_ms2, 
                                              out_path = out_path, 
                                              digits = digits), 
                              .scheduling = "dynamic")
  
  parallel::stopCluster(cl)
  
  oks <- purrr::map_lgl(out, ~ nrow(.x) > 0)
  out <- out[oks] %>% dplyr::bind_rows()
  
  saveRDS(out, file.path(out_path, "scores_simple.rds"))
  
  invisible(out)
}



