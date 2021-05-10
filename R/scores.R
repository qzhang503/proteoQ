#' Extracts the lists of theoretical or experimental m-over-z's.
#'
#' @param data Nested list of table from ion matches. Each table contains
#'   columns \code{theo} and \code{expt}.
#' @param col The name of a column where data will be extracted.
#' @import purrr
extract_matches_col <- function (data, col = "theo") {
  purrr::map(data, ~ {
    x <- .x %>% # by sequence
      purrr::map(`[[`, col) # by varmod positions
    
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
#' @inheritParams extract_matches_col
#' @examples
#' \donttest{
#' x <- flatten_pepouts(res$matches, "matches")
#' scores <- flatten_pepouts(scores, "pep_score")
#' }
flatten_pepouts <- function (data, outcol = "matches") {
  purrr::map(data, ~ {
    tib <- .x
    
    tib <- purrr::map(tib, ~ {
      tibble(pep_ivmod = names(.x), !!outcol := .x)
    }) 
    
    nm <- names(tib)
    
    tib <- purrr::map2(tib, nm, ~ {
      .x$pep_seq <- .y
      .x
    }) %>% 
      dplyr::bind_rows()
  }) 
}


#' Combines peptide scores to table.
#' 
#' @inheritParams flatten_pepouts
combine_pepvecs <- function (data, outcol = "pep_score") {
  purrr::imap(data, ~ {
    tibble(pep_seq = .y, pep_ivmod = names(.x), !!outcol := .x)
  }) %>% 
    dplyr::bind_rows()
}


#' Matches against secondary ions
#' 
#' @param theos A list of theoretical values.
#' @param expts A list of experimental values
#' @inheritParams calc_ms2ions
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
  
  ## matches btw. theoreticals and experimentals
  
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
  theos <- extract_matches_col(mts, "theo")
  
  # matches additionally against secondary ions
  mts2 <- match_secions(theos, expts, type_ms2ions, ppm_ms2, digits)
  
  # the number of secondary matches
  xs2 <- map(mts2, ~ {
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
      x <- nrow(.x) 
      x[x > N] <- N
      x
    })
  })
  
  # subtracts the counts of secondary b0, y0 matches etc. from noise
  ns <- map2(ms, xs2, ~ {
    ms_i <- .x
    xs2_i <- .y
    map2(ms_i, xs2_i, ~ N - .x - .y)
  })
  
  scores <- mapply(mdhyper, x = xs, m = ms, n = ns, k  = k, SIMPLIFY = FALSE) %>% 
    purrr::map(~ -log10(.x) * 10)
  
  scores <- scores %>% 
    purrr::map(~ {
      .x[.x > 250] <- 250
      .x
    })
  
  scores <- scores %>%
    combine_pepvecs("pep_score") %>% 
    dplyr::mutate(scan_num = unlist(entry$scan_num))
  
  dplyr::bind_cols(
    scores, 
    tibble(pri_matches = purrr::flatten(mts)), 
    tibble(sec_matches = purrr::flatten(mts2)))
}


#' Calculates the scores of peptides at an \code{aa_masses}.
#' 
#' @inheritParams calc_pepscores
calc_pepscores_i <- function (res, topn_ms2ions = 100, type_ms2ions = "by", 
                              ppm_ms2 = 25, out_path = "~/proteoQ/outs", 
                              digits = 5) {

  if (nrow(res) == 0) {
    scores <- tibble::tibble(
      pep_seq = as.character(), 
      pep_ivmod = as.character(), 
      pep_score = as.numeric(), 
      scan_num = as.integer(),
      pri_matches = list(), 
      sec_matches = list())
  } else {
    scores <- res %>% 
      split(., seq_len(nrow(.))) %>% 
      purrr::map(scalc_pepscores, topn_ms2ions, type_ms2ions, ppm_ms2, digits) %>% 
      dplyr::bind_rows()
  }
  
  res <- res %>% 
    dplyr::select(-c("matches")) %>% 
    dplyr::mutate(scan_num = as.numeric(scan_num)) %>% 
    dplyr::mutate(scan_title = as.character(scan_title), 
                  ms1_moverz = as.numeric(ms1_moverz), 
                  ms1_mass = as.numeric(ms1_mass), 
                  ms1_int = as.numeric(ms1_int), 
                  ms1_charge = as.character(ms1_charge), 
                  ret_time = as.integer(ret_time), 
                  ms2_n = as.integer(ms2_n), 
                  pep_fmod = as.character(pep_fmod), 
                  pep_vmod = as.character(pep_vmod), 
                  pep_nl = as.character(pep_nl)) %>% 
    dplyr::left_join(scores, by = "scan_num")
  
  res <- res %>% 
    dplyr::group_by(scan_num) %>% 
    dplyr::arrange(-pep_score) %>% 
    dplyr::mutate(pep_rank = row_number()) %>% 
    dplyr::mutate(pep_score = round(pep_score, digits = 1)) %>% 
    dplyr::ungroup() 
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
                              res = res, 
                              MoreArgs = list(topn_ms2ions = topn_ms2ions, 
                                              type_ms2ions = type_ms2ions, 
                                              ppm_ms2 = ppm_ms2, 
                                              out_path = out_path, 
                                              digits = digits), 
                              .scheduling = "dynamic")
  
  parallel::stopCluster(cl)
  
  oks <- purrr::map_lgl(out, ~ nrow(.x) > 0)
  out <- out[oks] %>% dplyr::bind_rows()
  
  saveRDS(out, file.path(out_path, "scores.rds"))

  invisible(out)
}



