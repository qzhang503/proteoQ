#' Reporter-ion quantitation.
#' 
#' @param scores The results from \link{calc_pepscores}.
#' @param quant A quantitation method. The default is "none". 
#' @param out_path The output path.
#' @inheritParams matchMS
calc_tmtint <- function (scores, 
                         quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"), 
                         ppm_ms2 = 25, 
                         out_path = "~/proteoQ/outs") {
  
  # file <- file.path(out_path, "scores.rds")
  # if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  # scores <- readRDS(file)
  
  val <- rlang::enexpr(quant)
  f <- match.call()[[1]]
  val <- match_valexpr(f = !!f, arg = "quant", val = !!val)
  
  if (val == "none") {
    out <- scores
  } else {
    nms_tmt6 <- c("126", "127N", "128N", "129N", "130N", "131N")
    
    nms_tmt10 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", 
                   "130N", "130C", "131N")
    
    nms_tmt11 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", 
                   "130N", "130C", "131N", "131C")
    
    nms_tmt16 <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", 
                   "130N", "130C", "131N", "131C", "132N", "132C", 
                   "133N", "133C", "134N")
    
    tmts <- c(
      `126` = 126.127726, `127N` = 127.124761, `127C` = 127.131080, 
      `128N` = 128.128115, `128C` = 128.134435, `129N` = 129.131470, 
      `129C` = 129.137790, `130N` = 130.134825, `130C` = 130.141145, 
      `131N` = 131.138180, `131C` = 131.144499, `132N` = 132.141535, 
      `132C` = 132.147855, `133N` = 133.14489, `133C` = 133.15121, 
      `134N` = 134.148245)
    
    theos <- switch(val, 
                    tmt6 = tmts %>% .[names(.) %in% nms_tmt6], 
                    tmt10 = tmts %>% .[names(.) %in% nms_tmt10], 
                    tmt11 = tmts %>% .[names(.) %in% nms_tmt11], 
                    tmt16 = tmts %>% .[names(.) %in% nms_tmt16], 
                    stop("Unknown TMt type.", call. = FALSE))
    
    ul <- switch(val, 
                 tmt6 = c(126.1, 131.2), 
                 tmt10 = c(126.1, 131.2), 
                 tmt11 = c(126.1, 131.2), 
                 tmt16 = c(126.1, 134.2), 
                 stop("Unknown TMt type.", call. = FALSE))
    
    out <- map2(scores$ms2_moverz, scores$ms2_int, 
                find_reporter_ints, 
                theos = theos, 
                ul = ul, 
                ppm_ms2 = ppm_ms2, 
                len = length(theos), 
                nms = names(theos)) %>% 
      bind_rows() %>% 
      bind_cols(scores, .)
  }

   out <- out %T>% 
    saveRDS(file.path(out_path, "peptides.rds"))
}


#' Finds the intensities of reporter-ions.
#' 
#' @param ms2_moverzs Numeric vector; a series of experimental MS2 m-over-z's
#'   (in the region of reporter ions).
#' @param ms2_ints Numeric vector; a series of experimental MS2 intensities (in
#'   the region of reporter ions).
#' @param theos The theoretical m-over-z of reporter ions.
#' @param ul The upper and lower bound for reporter-ion m-over-z's.
#' @param len The length of reporter-ion plexes.
#' @param nms The names of reporter-ion channels.
#' @inheritParams matchMS
find_reporter_ints <- function (ms2_moverzs, ms2_ints, theos, ul, ppm_ms2, 
                                len, nms) {
  
  range <- findInterval(ul, ms2_moverzs)

  ms <-ms2_moverzs[range[1]:range[2]]
  is <-ms2_ints[range[1]:range[2]]
  
  idxes <- find_reporters_ppm(theos, ms, ppm_ms2, len, nms)
  
  rptr_ints <- is[idxes] %>% 
    `names<-`(names(idxes))
  
  if (length(rptr_ints) < len) {
    es <- rep(NA, len) %>% 
      `names<-`(nms)
    
    es[names(rptr_ints)] <- rptr_ints
  }
  
  es
}

#' Finds the indexes of reporter ions.
#'
#' @param expts Numeric vector; a series of experimental MS2s (in the region of
#'   reporter ions).
#' @inheritParams find_reporter_ints
#' @return A vector of indexes
find_reporters_ppm <- function (theos, expts, ppm_ms2 = 25, len, nms) {
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_ms2, arr.ind = TRUE)

  idxes <- row_cols[, 2]
}



