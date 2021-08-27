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
add_seions <- function (ms2s, type_ms2ions = "by", digits = 5L) {
  
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


#' Matches two lists.
#' 
#' Without making a data frame.
#' 
#' @param a The left vector.
#' @param b The right vector.
#' @examples
#' \donttest{
#' a <- c(3, 4, 1, 2, 5)
#' b <- 2
#' 
#' list_leftmatch(a, b)
#' }
list_leftmatch <- function (a, b) {
  
  ord <- order(a, decreasing = TRUE)
  a <- a[ord]
  
  oks <- a %in% b
  
  b2 <- rep(NA, length(a))
  b2[oks] <- b
  
  b2
}


#' Helper for score calculations
#'
#' By the positions of variable modifications.
#' 
#' @param df Two lists of \code{theo} and matched \code{expt} m-over-z.
#' @param nms The names (character strings indicating the names and position of
#'   variable modifications).
#' @inheritParams calc_probi
#' @import dplyr
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @examples
#' \donttest{
#' df <- list(theo = c(390.2009, 550.2315, 710.2622, 809.3306, 880.3677, 
#'                     995.3946, 1151.4957, 175.1190, 290.1459, 361.1830, 
#'                     460.2514, 620.2821, 780.3127, 940.3434), 
#'            expt = c(390.2008, 550.2323, 710.2624, 809.3301, 880.3662, 
#'                     995.3970, NA, 175.1191, 290.1458, 361.1832, 
#'                     460.2517, 620.2880, 780.3126, 940.3438))
#' 
#' expt_moverzs <- c(110.0717, 112.0509, 112.0873, 113.0713, 115.0869, 
#'                   116.0167, 116.0709, 126.1280, 127.1250, 127.1313, 
#'                   128.1284, 128.1346, 129.1317, 129.1380, 130.0978, 
#'                   130.1350, 130.1413, 131.1384, 133.0432, 136.0619, 
#'                   139.0504, 157.1083, 158.0563, 158.0925, 159.0763, 
#'                   173.1498, 175.1191, 176.1223, 176.1600, 178.0646, 
#'                   186.1531, 188.1597, 212.1031, 230.1144, 230.1702, 
#'                   232.1115, 248.1808, 255.1082, 261.6220, 264.6251, 
#'                   273.1193, 275.6198, 284.1326, 284.6348, 290.1458, 
#'                   310.1302, 321.0684, 329.1280, 344.1566, 359.6652, 
#'                   361.1832, 362.2062, 390.2008, 407.2268, 420.1369, 
#'                   459.2220, 460.2517, 481.0993, 491.1726, 522.2371, 
#'                   539.7526, 540.2550, 550.2323, 551.2338, 567.2598, 
#'                   576.7457, 584.7561, 585.2570, 585.7582, 586.2587, 
#'                   619.2524, 620.2880, 682.2661, 683.2708, 710.2624, 
#'                   711.2637, 718.3199, 780.3126, 781.3342, 782.3379, 
#'                   809.3301, 810.3351, 880.3662, 881.3693, 921.3688, 
#'                   922.3726, 923.3927, 924.3849, 940.3438, 941.3491, 
#'                   995.3970, 996.3967, 997.3690, 998.3657, 1011.3803, 
#'                   1012.3803, 1013.3842, 1014.3911, 1015.3893, 1016.3904)
#' 
#' expt_ints <- c(12810.80, 14142.40, 58754.70, 12451.00, 29055.70, 
#'                45291.00, 63865.00, 250674.00, 261949.00, 179089.00, 
#'                253049.00, 190448.00, 240766.00, 275813.00, 28354.30, 
#'                219360.00, 189991.00, 229268.00, 60450.10, 12415.10, 
#'                11351.30, 17766.50, 29119.50, 105925.00, 10832.10, 
#'                15792.60, 707208.00, 29073.60, 12632.30, 18499.40, 
#'                18826.60, 33715.50, 12418.70, 18046.80, 164112.00, 
#'                15920.80, 13090.70, 24475.10, 14995.90, 49102.20, 
#'                56960.70, 17143.40, 93462.60, 15536.50, 23416.20, 
#'                15584.90, 30465.80, 12715.50, 31551.40, 12031.20, 
#'                22367.40, 55729.80, 77277.60, 19092.70, 29280.50, 
#'                17574.20, 21419.20, 10681.80, 8850.03, 27071.20, 
#'                25317.90, 13518.70, 55593.70, 14856.20, 30853.00, 
#'                11551.30, 22966.50, 515863.00, 243235.00, 11709.10, 
#'                47070.90, 19336.10, 57334.40, 11747.80, 147943.00, 
#'                42007.30, 18287.80, 51689.70, 62069.30, 13403.00, 
#'                84499.10, 24180.00, 47260.80, 13985.00, 14132.90, 
#'                10097.10, 12578.40, 13326.10, 49003.80, 12951.10, 
#'                98873.10, 38704.30, 12971.00, 8924.51, 89986.00, 
#'                162963.00, 119860.00, 114223.00, 125885.00, 32305.60)
#' 
#' calc_probi_byvmods(df, nms = "0000000", expt_moverzs, expt_ints, N = 190)
#' 
#' # 
#' df2 <- df
#' df2$expt[8] <- NA
#' calc_probi_byvmods(df2, nms = "0000000", expt_moverzs, expt_ints, N = 190)
#' }
calc_probi_byvmods <- function (df, nms, expt_moverzs, expt_ints, 
                                N, type_ms2ions = "by", topn_ms2ions = 100L, 
                                penalize_sions = TRUE, ppm_ms2 = 25L, 
                                digits = 5L) {
  
  # N - the total number of features (white and black balls)
  # k - the number of sampled features
  # m - the numbers of theoretical features (white balls)
  # n - the number of noise (black balls)
  
  m <- length(df[["theo"]])
  
  # OK: (N < m) -> (n < 0L)
  # if (N < m) N <- m
  
  ## matches additionally against secondary ions
  
  df2 <- add_seions(df[["theo"]], type_ms2ions = type_ms2ions, digits = digits) %>% 
    find_ppm_outer_bycombi(expt_moverzs, ppm_ms2) # 132 us
  
  df2[["theo"]] <- round(df2[["theo"]], digits = digits) # 4.4 us
  
  # subtracts `m` and the counts of secondary b0, y0 matches etc. from noise
  # (OK n < 0L)
  
  n <- N - m - sum(!is.na(df2$expt))
  
  ## step 0: the original tidyverse approach
  
  # m2 <- nrow(df2)
  
  # y2 <- df2 %>% 
  #   left_join(expts, by = "expt") %>% 
  #   `[[`("int") %>% 
  #   split(rep(seq_len(m2/m), each = m)) %>% 
  #   Reduce(`%+%`, .) %>% 
  #   data.frame(idx = seq_len(m), int2 = .)
  
  # y <- left_join(expts, df %>% mutate(idx = row_number()), by = "expt") %>% 
  #   dplyr::left_join(y2, by = "idx") %>% 
  #   mutate(int = ifelse(is.na(int2), int, int + int2)) %>% 
  #   select(-c("int2", "idx")) %>% 
  #   arrange(-int) %>% 
  #   mutate(k = row_number(), x = k - cumsum(is.na(theo))) %>% 
  #   filter(!is.na(theo))
  
  ## step 1: compiles secondary intensities
  
  # (1.1)
  len <- length(df2$expt) # .6 us
  
  i_se <- match(df2$expt, expt_moverzs) # secondary ions (df2) in expts; 4 us
  i_s <- which(!is.na(i_se)) # indexes in df2; 1.7 us
  i_e <- i_se[i_s] # indexes of the matches in expts; .4 us
  
  df2$int <- rep(NA, len) # matched intensities; .7 us
  df2$int[i_s] <- expt_ints[i_e] # the corresponding intensities from expts; 1.4 us
  
  # (1.2) collapse b0, b*, b2 etc.
  int2 <- df2$int %>% 
    split(rep(seq_len(len/m), each = m)) %>% 
    Reduce(`%+%`, .)
  
  y2 <- list(idx = 1:m, int2 = int2)
  
  ## step 2 (join expts and "theo", "idx" from the primary df): 
  # (y <- left_join(expts, df %>% mutate(idx = row_number()), by = "expt"))
  # 
  #   expts: "expt" (expt_moverzs), "int" (expt_ints)
  #   df: "theo" (m/z), "expt" (m/z), "idx" (indexes)
  #   -> y: "theo", "expt", "int" "idx"
  
  # (2.1)
  df$idx <- 1:m
  
  i_ep <- match(expt_moverzs, df$expt) # expts in primary ions (df); 4 us
  i_e <- which(!is.na(i_ep)) # indexes in expts; 1.7 us
  i_p <- i_ep[i_e] # indexes of matches in primary df; .4 us
  
  # (2.2)
  nu <- rep(NA, topn_ms2ions)
  y <- list(expt = expt_moverzs, int = expt_ints)
  
  y$theo <- nu
  y$idx <- nu
  
  y$theo[i_e] <- df$theo[i_p] # .5 us
  y$idx[i_e] <- df$idx[i_p]
  
  ## step 3: add `int2` from `y2`
  #  (dplyr::left_join(y, y2, by = "idx"))
  
  i_yy2 <- match(y$idx, y2$idx)
  i_y <- which(!is.na(i_yy2))
  i_y2 <- i_yy2[i_y]
  
  y$int2 <- nu
  y$int2[i_y] <- y2$int2[i_y2]
  
  ## step 4: collapses `int2` to `int`
  #  (mutate(y, int = ifelse(is.na(int2), int, int + int2)))

  y[["int"]] <- y[["int"]] %+% y[["int2"]]
  y[["int2"]] <- NULL
  y[["idx"]] <- NULL
  
  ## step 5: arrange(-int)
  
  idx <- order(y[["int"]], decreasing = TRUE, method = "radix", na.last = TRUE) # 16.4 us
  y[["expt"]] <- y[["expt"]][idx]
  y[["int"]] <- y[["int"]][idx]
  y[["theo"]] <- y[["theo"]][idx]
  
  ## step 6: mutate(k = row_number(), x = k - cumsum(is.na(theo)))
  
  k <- 1:topn_ms2ions
  x <- k - cumsum(is.na(y[["theo"]]))
  
  y$k <- k
  y$x <- x
  
  ## step 7: filter(!is.na(theo))
  
  idx <- !is.na(y[["theo"]])
  
  y[["expt"]] <- y[["expt"]][idx]
  y[["int"]] <- y[["int"]][idx]
  y[["theo"]] <- y[["theo"]][idx]
  y[["k"]] <- y[["k"]][idx]
  y[["x"]] <- y[["x"]][idx]
  
  ## Probability
  # note: x <= k <= x + n
  
  x <- y$x
  k <- y$k
  
  # (to have sufficient counts of noise)
  # (also guaranteed n > 0L)
  
  n <- max(n, topn_ms2ions + k[length(k)])
  
  pr <- min(mapply(dhyper, x[-c(1:2)], m, n, k[-c(1:2)]), na.rm = TRUE)
  
  list(pep_ivmod = nms, 
       pep_prob = pr, 
       pri_matches = list(df), 
       sec_matches = list(df2))
}


#' Helper for score calculations
#'
#' By peptides.
#' 
#' @param nms The names (of peptides).
#' @inheritParams calc_probi
#' @import dplyr
#' @importFrom purrr map2
#' @importFrom tibble tibble
calc_probi_bypep <- function (mts, nms, expt_moverzs, expt_ints, 
                              N, type_ms2ions, topn_ms2ions, 
                              penalize_sions, ppm_ms2, digits) {
  
  ## for different positions: $TNLAMMR$`0000500`, $TNLAMMR$`0000050`
  #    the same `pep_seq`, `theo_ms1` for different mod positions
  #    different `pep_ivmod`, `pep_prob`, `pri_matches`, `sec_matches`
  
  res <- map2(mts, names(mts), calc_probi_byvmods, 
              expt_moverzs = expt_moverzs, 
              expt_ints = expt_ints, 
              N = N, 
              type_ms2ions = type_ms2ions, 
              topn_ms2ions = topn_ms2ions, 
              penalize_sions = penalize_sions, 
              ppm_ms2 = ppm_ms2, 
              digits = digits)
  
  theo_ms1 <- attr(mts, "theo_ms1")
  
  len <- length(res)
  out <- vector("list", len)
  
  for (i in 1:len) {
    res_i <- res[[i]]
    
    out[[i]] <- list(
      pep_seq = nms,
      theo_ms1 = theo_ms1, 
      pep_ivmod = res_i$pep_ivmod, 
      pep_prob = res_i$pep_prob, 
      pri_matches = res_i["pri_matches"], 
      sec_matches = res_i["sec_matches"]
    )
  }
  
  out
}


#' Helper for score calculations
#'
#' @param mts Nested data frame of \code{theo} and matched \code{expt} m-over-z.
#' @param expt_moverzs Nested list of match and unmatched experimental m-over-z.
#' @param expt_ints Nested list of match and unmatched experimental intensity.
#' @param N Numeric; the number of MS2 features in an MGF query.
#' @inheritParams matchMS
#' @inheritParams calc_pepscores
#' @import dplyr
#' @importFrom purrr map
calc_probi <- function (mts, expt_moverzs, expt_ints, 
                        N, type_ms2ions = "by", topn_ms2ions = 100L, 
                        penalize_sions = TRUE, ppm_ms2 = 25L, digits = 5L) {
  
  out <- map2(mts, names(mts), calc_probi_bypep, 
              expt_moverzs = expt_moverzs, 
              expt_ints = expt_ints, 
              N = N, 
              type_ms2ions = type_ms2ions, 
              topn_ms2ions = topn_ms2ions, 
              penalize_sions = penalize_sions, 
              ppm_ms2 = ppm_ms2, 
              digits = digits)
  
  out <- unlist(out, recursive = FALSE, use.names = FALSE)
}


#' Calculates peptide scores by single MGF entries (scan numbers).
#' 
#' Each entry corresponds to a row in \code{ion_matches.rds}.
#' 
#' @param entry A row of data from \link{pmatch_bymgfs}.
#' @inheritParams matchMS
#' @import purrr
scalc_pepprobs <- function (entry, topn_ms2ions = 100L, type_ms2ions = "by", 
                            penalize_sions = FALSE, ppm_ms2 = 25L, digits = 4L) {

  # only one experimental set of values and thus `[[1]]`
  expt_moverzs <- entry$ms2_moverz[[1]]
  expt_ints <- entry[["ms2_int"]][[1]]

  ## matches between theoreticals and experimentals
  
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
  
  # (flattens by one level as is a list-column)
  mts <- entry$matches[[1]]
  
  N <- entry$ms2_n[[1]]
  topn_ms2ions <- min(topn_ms2ions, N)
  
  out <- calc_probi(mts = mts, 
                    expt_moverzs = expt_moverzs, 
                    expt_ints = expt_ints, 
                    N = N, 
                    type_ms2ions = type_ms2ions, 
                    topn_ms2ions = topn_ms2ions, 
                    penalize_sions = penalize_sions, 
                    ppm_ms2 = ppm_ms2, 
                    digits = digits)

  uniq_id <- unlist(entry$uniq_id) # 2.3 us
  out <- map(out, ~ {.x$uniq_id <- uniq_id; .x}) # 23.5 us

  invisible(out)
}


#' Calculates the scores of peptides at an \code{aa_masses}.
#' 
#' @inheritParams calc_pepscores
calc_pepprobs_i <- function (res, topn_ms2ions = 100L, type_ms2ions = "by", 
                             penalize_sions = FALSE, ppm_ms2 = 25L, 
                             out_path = "~/proteoQ/outs", digits = 5L) {

  if (nrow(res) == 0L) {
    probs <- data.frame(
      pep_seq = as.character(), 
      pep_ivmod = as.character(), 
      pep_prob = as.numeric(), 
      pri_matches = list(), 
      sec_matches = list(), 
      scan_num = as.integer())
  } else {
    probs <- res %>% 
      split(., seq_len(nrow(.))) %>% 
      lapply(scalc_pepprobs, 
             topn_ms2ions = topn_ms2ions, 
             type_ms2ions = type_ms2ions, 
             penalize_sions = penalize_sions, 
             ppm_ms2 = ppm_ms2, 
             digits = digits) # 2.79 ms

    probs <- unlist(probs, recursive = FALSE, use.names = FALSE)
    probs <- dplyr::bind_rows(probs) # 276 us
  }
  
  invisible(probs)
}


#' Calculates the scores of peptides.
#'
#' @param res The result from \link{pmatch_bymgfs}. At the \code{NULL} default,
#'   the file \code{ion_matches.rds} will be used.
#' @param penalize_sions Logical; if TRUE, penalizes secondary ions of b0, y0
#'   etc. with lower weights in peptide scoring.
#' @inheritParams matchMS
#' @import parallel
#' @export
calc_pepscores <- function (topn_ms2ions = 100L, type_ms2ions = "by", 
                            target_fdr = 0.01, fdr_type = "psm", 
                            min_len = 7L, max_len = 100L, 
                            penalize_sions = FALSE, ppm_ms2 = 25L, 
                            out_path = "~/proteoQ/outs", digits = 5L) {

  dir.create(file.path(out_path, "temp"), recursive = TRUE, showWarnings = FALSE)
  message("Calculating peptide scores.")
  
  ## --- Target ---
  list_t <- list.files(path = file.path(out_path, "temp"), 
                       pattern = "^ion_matches_\\d+\\.rds$")
  
  if (length(list_t) == 0L) stop("Target matches not found.", call. = FALSE)
  
  nms_t <- gsub("^ion_matches_(\\d+)\\.rds$", "\\1", list_t)
  
  ord <- order(as.integer(nms_t))
  nms_t <- nms_t[ord]
  list_t <- list_t[ord]
  
  for (i in seq_along(list_t)) {
    calcpepsc(file = list_t[i], 
              topn_ms2ions = topn_ms2ions, 
              type_ms2ions = type_ms2ions, 
              penalize_sions = penalize_sions, 
              ppm_ms2 = ppm_ms2, 
              out_path = out_path, 
              digits = digits)
    
    gc()
  }
  
  ## --- Decoy ---
  list_d <- list.files(path = file.path(out_path, "temp"), 
                       pattern = "^ion_matches_rev_[0-9]+\\.rds$")
  
  if (length(list_d) == 0L) stop("Decoy matches not found.", call. = FALSE)
  
  if (length(list_d) > 1L) {
    warning("More than one decoy matches found: '", decoy, "'.\n", 
            "The first will be used.", call. = FALSE)
    list_d <- list_d[[1]]
  }
  
  nms_d <- gsub("^ion_matches_(rev_[0-9]+)\\.rds$", "\\1", list_d)
  
  rev <- calcpepsc(file = list_d, 
                   topn_ms2ions = topn_ms2ions, 
                   type_ms2ions = type_ms2ions, 
                   penalize_sions = penalize_sions, 
                   ppm_ms2 = ppm_ms2, 
                   out_path = out_path, 
                   digits = digits) %>% 
    list() %>% 
    `names<-`(nms_d)

  ## --- FDR (target + decoy) --- 
  # Reload targets
  out <- local({
    files <- list.files(path = file.path(out_path, "temp"), 
                        pattern = "^pepscores_\\d+\\.rds$")
    
    if (length(files) == 0L) stop("Target scores not found.", call. = FALSE)
    
    nms <- gsub("^pepscores_(\\d+)\\.rds$", "\\1", files)
    
    ord <- order(as.integer(nms))
    nms <- nms[ord]
    files <- files[ord]
    
    out <- purrr::map(files, ~ readRDS(file.path(out_path, "temp", .x))) %>% 
      `names<-`(nms)
    
    out <- c(out, rev)
  })
  
  rm(list = c("rev"))
  gc()
  
  ## Probability cut-offs
  prob_cos <- calc_pepfdr(out, 
                          nms = nms_d, 
                          target_fdr = target_fdr, 
                          fdr_type = fdr_type, 
                          min_len = min_len, 
                          max_len = max_len, 
                          out_path = out_path)
  
  # homolog co
  prob_cos <- local({
    idxes <- prob_cos %>% .[. <= target_fdr]
    
    if (length(idxes) > 0L) {
      fct_homol <- target_fdr/max(idxes, na.rm = TRUE)
    } else {
      fct_homol <- 1L
    }
    
    prob_cos <- prob_cos * fct_homol
  })
  
  ## Outputs
  prob_cos <- prob_cos%>% 
    data.frame(pep_len = as.numeric(names(.)), pep_prob_co = .) 

  fct_score <- 10
  
  oks <- purrr::map_lgl(out, ~ nrow(.x) > 0L)

  out <- out[oks] %>% 
    dplyr::bind_rows()
  
  # Adjusted p-values
  out <- out %>% 
    dplyr::left_join(prob_cos, by = "pep_len") %>% 
    dplyr::mutate(pep_issig = ifelse(pep_prob <= pep_prob_co, TRUE, FALSE), 
                  pep_adjp = p.adjust(pep_prob, "BH"))
  
  prob_cos <- purrr::map_dbl(prob_cos$pep_prob_co, ~ {
    row <- abs(log10(out$pep_prob/.x)) %>% which.min() 
    out[row, ]$pep_adjp
  }) %>% 
    dplyr::bind_cols(prob_cos, pep_adjp_co = .) %T>% 
    saveRDS(file.path(out_path, "temp", "pep_prob_cos.rds"))
  
  out <- out %>% 
    dplyr::left_join(prob_cos[, c("pep_len", "pep_adjp_co")], by = "pep_len") %>% 
    dplyr::mutate(pep_score = -log10(pep_adjp) * fct_score, 
                  pep_score = ifelse(pep_score > 250, 250, pep_score), 
                  pep_score_co = -log10(pep_adjp_co) * fct_score) %>% 
    dplyr::select(-c("pep_prob", "pep_adjp", "pep_prob_co", "pep_adjp_co"))
  
  invisible(out)
}


#' Helper of \link{calc_pepscores}.
#' 
#' @param file A file name of \code{ion_matches_}.
#' @inheritParams matchMS
calcpepsc <- function (file, topn_ms2ions = 100L, type_ms2ions = "by", 
                       penalize_sions = FALSE, ppm_ms2 = 25L, out_path = NULL, 
                       digits = 4L) {
  
  n_cores <- detect_cores()
  n_cores2 <- n_cores^2L
  
  df <- readRDS(file.path(out_path, "temp", file))
  df[["uniq_id"]] <- paste(df[["scan_num"]], df[["raw_file"]], sep = "@")
  
  esscols <- c("ms2_moverz", "ms2_int", "matches", "ms2_n", "uniq_id")
  df2 <- df[, -which(names(df) %in% esscols)]
  df <- df[, esscols, drop = FALSE]
  
  # otherwise, chunksplit return NULL
  #   -> res[[i]] <- NULL 
  #   -> length(res) shortened by 1
  
  if (!is.null(df)) {
    df <- suppressWarnings(chunksplit(df, n_cores2, "row"))
    gc()
  }
  
  if (length(df) >= n_cores2) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    parallel::clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
    clusterExport(cl, list("scalc_pepprobs"), envir = environment(proteoQ:::scalc_pepprobs))
    
    probs <- parallel::clusterApplyLB(cl, df, 
                                      calc_pepprobs_i, 
                                      topn_ms2ions = topn_ms2ions, 
                                      type_ms2ions = type_ms2ions, 
                                      penalize_sions = penalize_sions, 
                                      ppm_ms2 = ppm_ms2,
                                      out_path = out_path, 
                                      digits = digits)
    
    parallel::stopCluster(cl)
    gc()
    
    probs <- probs %>% dplyr::bind_rows()
  } else {
    if (is.data.frame(df)) df <- list(df)
    
    probs <- purrr::map(df, calc_pepprobs_i, 
                        topn_ms2ions = topn_ms2ions, 
                        type_ms2ions = type_ms2ions, 
                        penalize_sions = penalize_sions, 
                        ppm_ms2 = ppm_ms2,
                        out_path = out_path, 
                        digits = digits) %>% 
      dplyr::bind_rows()
  }
  
  
  gc()
  
  ## Reassemble `df`
  df <- dplyr::bind_rows(df)
  df <- df[, -which(names(df) == "matches"), drop = FALSE]
  
  df <- dplyr::bind_cols(df, df2)
  df <- quick_rightjoin(df, probs, "uniq_id")
  df <- df[, -which(names(df) == "uniq_id"), drop = FALSE]
  df <- post_pepscores(df)
  
  rm(list = c("df2", "probs"))
  gc()
  
  ## Outputs 
  # (can be decoy => .*, not \\d+)
  idx <- gsub("^ion_matches_(.*)\\.rds$", "\\1", file)
  
  cols_a <- c("raw_file", "pep_mod_group", "scan_num")
  cols_b <- c("ms2_moverz", "ms2_int", "pri_matches", "sec_matches")
  
  saveRDS(df[, c(cols_a, cols_b), drop = FALSE], 
          file.path(out_path, "temp", paste0("list_table_", idx, ".rds")))
  
  df <- df[, -which(names(df) %in% cols_b), drop = FALSE] %T>% 
    saveRDS(file.path(out_path, "temp", paste0("pepscores_", idx, ".rds")))
  
  invisible(df)
}


#' Cleanups post \link{calc_pepscores}.
#' 
#' @param df A results after pep_scores.
post_pepscores <- function (df) {
  df[["scan_num"]] <- as.numeric(df[["scan_num"]]) # 70 us
  df[["pep_len"]] <- stringr::str_length(df[["pep_seq"]]) # 6 us
  df[["scan_title"]] <- as.character(df[["scan_title"]]) # 69 us
  df[["ms1_moverz"]] <- as.numeric(df[["ms1_moverz"]]) # 4.5 us
  df[["ms1_mass"]] <- as.numeric(df[["ms1_mass"]]) # 5.4 us
  df[["ms1_int"]] <- as.numeric(df[["ms1_int"]]) # 5.5 us
  df[["ms1_charge"]] <- as.character(df[["ms1_charge"]]) # 5.9 us
  df[["ret_time"]] <- as.integer(df[["ret_time"]]) # 5.6 us
  df[["ms2_n"]] <- as.integer(df[["ms2_n"]]) # 5.3 us
  df[["pep_fmod"]] <- as.character(df[["pep_fmod"]]) # 4.6 us
  df[["pep_vmod"]] <- as.character(df[["pep_vmod"]]) # 4.8 us
  
  invisible(df)
}


#' Helper of \link{calc_pepfdr}.
#'
#' Calculates the probability cut-off for target-decoy pairs at a given peptide
#' length.
#' 
#' @param td A target-decoy pair.
#' @param len Numeric; the length of peptides.
#' @inheritParams matchMS
probco_bypeplen <- function (len, td, fdr_type, target_fdr, out_path) {
  
  td <- td %>% dplyr::filter(pep_len == len)
  
  if (fdr_type %in% c("peptide", "protein")) {
    td <- td %>% 
      dplyr::arrange(pep_seq, pep_prob) %>% 
      dplyr::group_by(pep_seq) %>% 
      dplyr::filter(row_number() == 1L) %>% 
      dplyr::ungroup()
  }
  
  td <- td %>% 
    dplyr::select(pep_prob, pep_isdecoy) %>% 
    dplyr::arrange(pep_prob) %>% 
    dplyr::mutate(total = row_number()) %>% 
    dplyr::mutate(decoy = cumsum(pep_isdecoy)) %>% 
    dplyr::mutate(fdr = decoy/total) %>% 
    dplyr::mutate(pep_score = -log10(pep_prob) * 10)
  
  # ---
  count <- nrow(td)
  
  if (count <= 200L) return(NA)
  
  # ---
  rows <- which(td$fdr <= target_fdr) 
  
  if (length(rows) > 0L) {
    row <- max(rows, na.rm = TRUE)
    score_co <- td[row, "pep_score"] %>% unlist(use.names = FALSE)
    
    # fittings
    # (the data range can affect the fitting)
    df <- data.frame(x = td[["pep_score"]], y = td[["fdr"]]) %>% 
      { if (target_fdr <= .05) dplyr::filter(., y <= .05) else . }

    fit <- suppressWarnings(
      tryCatch(
        nls(y ~ SSlogis(x, Asym, xmid, scal), data = df, 
            control = list(tol = 1e-03, warnOnly = TRUE), 
            algorithm = "port"), 
        error = function (e) NA)
    )
    
    # p <- suppressWarnings(
    #  ggplot(df, aes(x = x, y = y)) + 
    #     geom_point() + 
    #     stat_smooth(method = "nls", formula = y ~ SSlogis(x, Asym, xmid, scal), 
    #                 se = FALSE) + 
    #     labs(title = paste("pep_len = ", len))
    # )

    # try(suppressMessages(
    #     ggsave(filename = file.path(out_path, "temp", paste0("peplen_sco_", len, ".png")))
    #   ))
    
    if (!all(is.na(fit))) {
      newx <- min(df$x, na.rm = TRUE):max(df$x, na.rm = TRUE)
      newy <- predict(fit, data.frame(x = newx)) %>% `names<-`(newx)
      
      # NA if not existed
      score_co2 <- which(newy <= target_fdr)[1] %>% names() %>% as.numeric()
      score_co <- min(score_co, score_co2, na.rm = TRUE)
      
      rm(list = c("newx", "newy", "score_co2"))
    }
    
    rm(list = c("df", "fit"))
    
    prob_co <- 10^(-score_co/10)
    
  } else {
    prob_co <- NA
  }
  
  names(prob_co) <- count

  invisible(prob_co)
}


#' Calculates the cut-off score at a peptide FDR.
#'
#' Needs \code{min_len} and \code{max_len} since the target-decoy pair may not
#' cover all \code{pep_len} values.
#'
#' @param out The output from \link{calc_pepscores}.
#' @param nms The name(s) of \code{out} that correspond(s) to decoy results.
#' @param target_fdr Numeric; the levels of false-discovery rate (FDR).
#' @param fdr_type Character string; the type of FDR for controlling.
#' @inheritParams matchMS
calc_pepfdr <- function (out, nms = "rev_2", target_fdr = .01, fdr_type = "psm", 
                         min_len = 7L, max_len = 100L, out_path) {

  find_optlens <- function (all_lens, counts, min_count = 2000L) {
    idxes <- which(counts >= min_count)
    
    if (length(idxes) > 0L) {
      return(all_lens[idxes])
    } else {
      find_optlens(all_lens, min_count/2L)
    }
  }
  
  # ---
  if (!is.null(nms)) {
    nms_t <- gsub("^rev_", "", nms)
    
    # keeps the best hit for each `scan_num`
    # (separated bests for targets and decoys)
    
    out[[nms_t]] <- out[[nms_t]] %>% 
      dplyr::group_by(scan_num, raw_file) %>% 
      dplyr::arrange(pep_prob) %>% 
      dplyr::filter(row_number() == 1L) %>% 
      dplyr::ungroup()
    
    out[[nms]] <- out[[nms]] %>% 
      dplyr::group_by(scan_num, raw_file) %>% 
      dplyr::arrange(pep_prob) %>% 
      dplyr::filter(row_number() == 1L) %>% 
      dplyr::ungroup()
    
    td <- out[c(nms_t, nms)]
    
    # decoy sequences may be present in targets
    
    dpeps <- unique(td[[nms]][["pep_seq"]])
    tpeps <- unique(td[[nms_t]][["pep_seq"]])
    
    dpeps <- dpeps %>% 
      .[! . %in% tpeps]
    
    td[[nms]] <- td[[nms]] %>% 
      dplyr::filter(pep_seq %in% dpeps)
    
    rm(list = c("dpeps", "tpeps", "out"))
    gc()
    
    # ---
    #  keeps the best hit for each `scan_num`
    # (combined best for targets and decoys)
    
    td <- td %>% 
      dplyr::bind_rows() %>% 
      dplyr::group_by(scan_num, raw_file) %>% 
      dplyr::arrange(pep_prob) %>% 
      dplyr::filter(row_number() == 1L) %>% 
      dplyr::ungroup()
    
    all_lens <- sort(unique(td$pep_len))
    
    prob_cos <- all_lens %>% 
      lapply(probco_bypeplen, td, fdr_type, target_fdr, out_path) %>% 
      unlist()
    
    if (all(is.na(prob_cos))) {
      stop("Cannot calculate peptide FDR; contact the developer.")
    }
    
    counts <- as.numeric(names(prob_cos))
    names(counts) <- all_lens
    names(prob_cos) <- all_lens
    
    lens <- find_optlens(all_lens, counts, 2000L)
    prob_cos <- prob_cos %>% .[names(.) %in% lens]
    counts <- counts %>% .[names(.) %in% lens]
    
    # ---
    best_score_co <- prob_cos %>% 
      .[which_topx(., n = 1L)] %>% 
      log10() %>% `-`
    
    valley <- as.numeric(names(best_score_co))
    
    # ---
    df <- data.frame(x = as.numeric(names(prob_cos)), y = -log10(prob_cos))
    fit <- lm(y ~ splines::ns(x, 4), df)
    
    # ggplot(df, aes(x = x, y = y)) + geom_point() +
    #   stat_smooth(method = "lm", formula = y ~ splines::ns(x, 4), se = FALSE)
    
    newx <- min_len : max_len
    newy <- predict(fit, data.frame(x = newx)) %>% 
      `names<-`(newx)
    newy[which(names(newy) == valley):length(newy)] <- best_score_co
    
    prob_cos <- 10^-newy
    
  } else {
    seqs <- min_len : max_len
    prob_cos <- rep(.05, length(seqs))
    names(prob_cos) <- seqs
  }
  
  invisible(prob_cos)
}


#' Calculates the cut-offs of protein scores.
#'
#' @param out An output from upstream steps.
#' @inheritParams calc_pepfdr
calc_protfdr <- function (out, target_fdr = .01) {
  
  message("Calculating peptide-protein FDR.")
  
  # target-decoy pair
  nms_d <- unique(out$pep_mod_group) %>% 
    .[grepl("^rev_\\d+", .)]
  
  nms_t <- gsub("^rev_", "", nms_d)
  
  # out2 <- out %>% dplyr::filter(pep_mod_group %in% c(nms_t, nms_d), pep_issig, pep_rank == 1L)
  out2 <- out %>% dplyr::filter(pep_mod_group %in% c(nms_t, nms_d), pep_issig)

  # score cut-offs as a function of prot_n_pep
  max_n_pep <- max(out$prot_n_pep, na.rm = TRUE)
  all_n_peps <- unique(out$prot_n_pep)

  # protein enrichment score cut-offs
  score_co <- out2 %>% 
    split(.$prot_n_pep) %>% 
    purrr::map_dbl(calc_protfdr_i, target_fdr) %>% 
    fit_protfdr(max_n_pep) %>% 
    dplyr::filter(prot_n_pep %in% all_n_peps) %>% 
    dplyr::rename(prot_es_co = prot_score_co)
  
  # add protein enrichment score
  prot_es <- out %>% 
    dplyr::group_by(prot_acc, pep_seq) %>% 
    dplyr::arrange(-pep_score) %>% 
    dplyr::filter(row_number() == 1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(pep_issig) %>% 
    dplyr::mutate(pep_es = pep_score - pep_score_co) %>% 
    dplyr::group_by(prot_acc) %>% 
    dplyr::summarise(prot_es = max(pep_es, na.rm = TRUE))
  
  out <- out %>% 
    dplyr::left_join(prot_es, by = "prot_acc")

  out <- out %>% 
    dplyr::left_join(score_co, by = "prot_n_pep") %>% 
    dplyr::mutate(prot_issig = ifelse(prot_es >= prot_es_co, TRUE, FALSE)) %>% 
    dplyr::mutate(pep_score = round(pep_score, digits = 1), 
                  pep_score_co = round(pep_score_co, digits = 1), 
                  prot_es = round(prot_es, digits = 1), 
                  prot_es_co = round(prot_es_co, digits = 1))
  
  rm(list = c("out2", "prot_es"))
  gc()

  invisible(out)
}


#' Helper of \link{calc_protfdr}.
#' 
#' The \code{prot_score} are only for a pair of target-decoy and thus not
#' exported. The offset of \code{base_score} not applied either.
#' 
#' @param td A data frame with paired target-decoys.
#' @inheritParams calc_pepfdr
calc_protfdr_i <- function (td, target_fdr = .01) {

  td <- td %>% 
    dplyr::group_by(prot_acc, pep_seq) %>% 
    dplyr::arrange(-pep_score) %>% 
    dplyr::filter(row_number() == 1L) %>% 
    dplyr::ungroup()

  # no decoys
  if (sum(td$pep_isdecoy) == 0L) {
    return(0L)
  }
  
  # all decoys
  if (sum(!td$pep_isdecoy) == 0L) {
    if (nrow(td) <= 5L) {
      return(0L)
    } else {
      return(20L)
    }
  }
  
  # both targets and decoys
  if (nrow(td) <= 20L) {
    return(1L)
  }
  
  prot_scores <- td %>% 
    dplyr::mutate(prot_es = pep_score - pep_score_co) %>% 
    dplyr::group_by(prot_acc) %>% 
    dplyr::summarise(prot_es = max(prot_es, na.rm = TRUE)) 
  
  td <- td %>% 
    dplyr::left_join(prot_scores, by = "prot_acc") %>% 
    # dplyr::filter(!(pep_isdecoy & prot_es == 250L)) %>% 
    dplyr::arrange(-prot_es) %>% 
    dplyr::mutate(total = row_number()) %>% 
    dplyr::mutate(decoy = cumsum(pep_isdecoy)) %>% 
    dplyr::mutate(fdr = decoy/total)

  rm(list = "prot_scores")
  
  # row <- which(td$fdr <= target_fdr) %>% max(na.rm = TRUE)
  rows <- which(td$fdr <= target_fdr)
  if (length(rows) > 0L) {
    row <- max(rows, na.rm = TRUE)
  } else {
    row <- -Inf
  }

  if (row == -Inf) {
    score_co <- 0L
    score_co2 <- score_co
  } else {
    score_co <- td[row, ]$prot_es

    score_co2 <- local({
      df <- data.frame(x = td[["prot_es"]], y = td[["fdr"]])
      
      fit <- suppressWarnings(
        tryCatch(
          nls(y ~ SSlogis(x, Asym, xmid, scal), data = df, 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA)
      )
      
      # ggplot(df, aes(x = x, y = y)) + geom_point() + 
      #  stat_smooth(method = "nls", formula = y ~ SSlogis(x, Asym, xmid, scal), 
      #  se = FALSE)

      if (all(is.na(fit))) {
        score_co2 <- score_co
      } else {
        min_score <- min(df$x, na.rm = TRUE)
        max_score <- max(df$x, na.rm = TRUE)
        newx <- min_score : max_score
        newy <- predict(fit, data.frame(x = newx)) %>% `names<-`(newx)
        
        # NA if not existed
        score_co2 <- which(newy <= target_fdr)[1] %>% names() %>% as.numeric()
      }

      invisible(score_co2)
    })
  }
  
  invisible(min(score_co, score_co2, na.rm = TRUE))
}


#' Fits the cut-offs of protein scores.
#'
#' Assumed a sigmoidal function.
#'
#' @param vec Named numeric vector. The values are the score cut-offs from
#'   \link{calc_protfdr}. The names correspond to the number of peptides being
#'   identified.
#' @param max_n_pep Integer; the maximum value of \code{prot_n_pep} for
#'   prediction.
fit_protfdr <- function (vec, max_n_pep = 1000L) {
  
  if (length(vec) <= 10L) {
    return(data.frame(prot_n_pep = as.numeric(names(vec)), 
                      prot_score_co = vec))
  }

  rv <- rev(vec)
  df <- data.frame(x = as.numeric(names(rv)), y = rv)
  elbow <- min(df[which(df$y == min(df$y, na.rm = TRUE)), "x"], na.rm = TRUE)
  amp <- max(df$y, na.rm = TRUE) * .8
  sca <- 0.5
  
  f <- function (x, m = 0, s = 1, a = 1) { a - a / (1 + exp(-(x-m)/s)) }
  
  fit <- suppressWarnings(
    tryCatch(
      nls(y ~ f(x, m, s, a), data = df, 
          start = list(a = amp, m = elbow, s = sca), 
          control = list(tol = 1e-03, warnOnly = TRUE), 
          algorithm = "port"), 
      error = function (e) NA)
  )
  
  # should not occur
  if (all(is.na(fit))) {
    fits <- suppressWarnings(
      map(seq_len(elbow-1), ~ {
        tryCatch(
          nls(y ~ f(x, m, s, a), data = df, 
              start = list(a = amp, m = .x, s = sca), 
              control = list(tol = 1e-03, warnOnly = TRUE), 
              algorithm = "port"), 
          error = function (e) NA)
      })
    )
    
    if (all(is.na(fits))) {
      fit <- NA
    } else {
      fit <- fits %>% .[!is.na(.)] %>% 
        .[[length(.)]]
    }
  }
  
  # ---
  newx <- seq(1, max_n_pep, by = 1)
  
  out <- data.frame(
    prot_n_pep = newx, 
    prot_score_co = predict(fit, data.frame(x = newx))) # %>% 
    # dplyr::mutate(prot_score_co = ifelse(prot_n_pep >= 10L, 0L, prot_score_co))

  invisible(out)
}


#' Helper of \link{calc_protfdr}.
#' 
#' @param df A data subset at a given \code{prot_n_pep}.
calc_protscore_i <- function (df) {
  
  df %>% 
    dplyr::group_by(prot_acc, pep_seq) %>% 
    dplyr::arrange(-pep_score) %>% 
    dplyr::filter(row_number() == 1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(prot_acc) %>% 
    dplyr::summarise(prot_score = max(pep_score)) %>% 
    dplyr::right_join(df, by = "prot_acc") %>% 
    dplyr::mutate(prot_issig = ifelse(prot_score >= prot_score_co, TRUE, FALSE))
}


#' Calculates peptide scores by single MGF entries (scan numbers).
#'
#' Not currently used. Each entry corresponds to a row in
#' \code{ion_matches.rds}.
#'
#' @param entry A row of data from \link{pmatch_bymgfs}.
#' @inheritParams matchMS
#' @import purrr
scalc_pepscores_static <- function (entry, topn_ms2ions = 100, type_ms2ions = "by", 
                                    ppm_ms2 = 25, digits) {
  
  ## only one experimental set of values and thus `[[1]]`
  
  expts <- entry$ms2_moverz[[1]]

  ## matches between theoreticals and experimentals
  
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
  
  scores <- mapply(mdhyper, x = xs, m = ms, n = ns, k  = k, 
                   SIMPLIFY = FALSE) %>% 
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


