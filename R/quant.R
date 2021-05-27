#' Reporter-ion quantitation.
#' 
#' @param data An upstream result from \link{matchMS}.
#' @param quant A quantitation method. The default is "none". 
#' @param ppm_reporters The mass tolerance of MS2 reporter ions.
calc_tmtint <- function (data = NULL, 
                         quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"), 
                         ppm_reporters = 10) {
  
  val <- rlang::enexpr(quant)
  f <- match.call()[[1]]
  val <- match_valexpr(f = !!f, arg = "quant", val = !!val)
  
  if (val == "none") {
    out <- data
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
    
    out <- map2(data$ms2_moverz, data$ms2_int, 
                find_reporter_ints, 
                theos = theos, 
                ul = ul, 
                ppm_reporters = ppm_reporters, 
                len = length(theos), 
                nms = names(theos)) %>% 
      bind_rows() %>% 
      bind_cols(data, .)
  }

   invisible(out)
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
#' @examples 
#' \donttest{
#' ms2_moverzs <- c(112.0873, 126.1280, 127.1251, 127.1313, 128.1250, 
#'                  128.1284, 128.1347, 129.1317, 129.1380, 130.0654, 
#'                  130.1351, 130.1413, 131.1384)
#'                  
#' ms2_ints <- c(5113.79, 135569.00, 120048.00, 122599.00, 3397.98, 
#'               140551.00, 144712.00, 103166.00, 145452.00, 3851.82, 
#'               148218.00, 135393.00, 131215.00)
#'               
#' theos <- c(126.1277, 127.1248, 127.1311, 128.1281, 128.1344, 
#'            129.1315, 129.1378, 130.1348, 130.1411, 131.1382)
#' names(theos) <- c("126", "127N", "127C", "128N", "128C", 
#'                   "129N", "129C", "130N", "130C", "131N")
#' 
#' ppm_reporters <- 10
#' ul <- c(126.1, 131.2)
#' len <- 10
#' nms <- names(theos)
#' 
#' x <- find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 10, 
#'                         len , nms)
#'                         
#' x <- find_reporter_ints(ms2_moverzs, ms2_ints, theos, ul, ppm_reporters = 25, 
#'                         len , nms)
#' }
find_reporter_ints <- function (ms2_moverzs, ms2_ints, theos, ul, 
                                ppm_reporters = 10, len, nms) {

  range <- findInterval(ul, ms2_moverzs)

  ms <-ms2_moverzs[range[1]:range[2]]
  is <-ms2_ints[range[1]:range[2]]

  idxes <- find_reporters_ppm(theos, ms, ppm_reporters, len, nms)
  
  if (is_empty(idxes)) {
    return(rep(NA, len) %>% `names<-`(nms))
  }
  
  # 126      127N      127C      128N      128N      128C      
  # 135569.00 120048.00 122599.00   3397.98 140551.00 144712.00 
  
  if (anyDuplicated(names(idxes))) {
    idxes <- idxes %>% 
      split(., names(.)) %>% 
      imap_int(~ {
        if (length(.x) > 1L) {
          p <- which.min(abs(ms[.x] - theos[.y]))
          .x <- .x[p]
        } 
        
        .x
      }) %>% 
      .[names(theos)]
  }
  
  rptr_ints <- is[idxes] %>% 
    `names<-`(names(idxes))
  
  if (length(rptr_ints) < len) {
    es <- rep(NA, len) %>% 
      `names<-`(nms)
    
    es[names(rptr_ints)] <- rptr_ints
  } else {
    es <- rptr_ints
  }
  
  es
}

#' Finds the indexes of reporter ions.
#'
#' @param expts Numeric vector; a series of experimental MS2s (in the region of
#'   reporter ions).
#' @inheritParams find_reporter_ints
#' @return A vector of indexes
find_reporters_ppm <- function (theos, expts, ppm_reporters = 10, len, nms) {
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_reporters, arr.ind = TRUE)

  row_cols[, 2]
}


#' Adds prot_acc to a peptide table
#' 
#' @param out_path An output path.
#' @param df The results after scoring.
add_prot_acc <- function (df, out_path = "~/proteoQ/outs") {
  # theoretical
  bins <- list.files(path = file.path(.path_fasta, "pepmasses", .time_stamp), 
                     pattern = "binned_theopeps_\\d+\\.rds$", 
                     full.names = TRUE)
  
  theopeps <- purrr::map(bins, ~ {
    x <- readRDS(.x) %>% 
      dplyr::bind_rows() %>% 
      dplyr::select(c("prot_acc", "pep_seq"))
  }) %>% 
    dplyr::bind_rows()
  
  theopeps <- theopeps %>% 
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>% 
    dplyr::filter(!duplicated(pep_prot.)) %>% 
    dplyr::select(-pep_prot.) %>% 
    dplyr::select(c("prot_acc", "pep_seq"))
  
  # adds `prot_acc`
  # df <- readr::read_tsv(file.path(out_path, "peptides.txt"))
  
  pep_seqs <- unique(df$pep_seq)
  
  # decoys kept
  theopeps %>% 
    dplyr::filter(.data$pep_seq %in% .env$pep_seqs) %>% 
    dplyr::right_join(df, by = "pep_seq")
}


