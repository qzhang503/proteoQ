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
    message("Calculating reporter-ion intensities.")
    
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
  
  # adds `prot_acc` (with decoys being kept)
  theopeps %>% 
    dplyr::filter(pep_seq %in% unique(df$pep_seq)) %>% 
    dplyr::right_join(df, by = "pep_seq")
}


#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
groupProts <- function (df, out_path = "~/proteoQ/outs") {
  message("Grouping proteins by families.")
  
  # --- protein ~ peptide map ---
  mat <- local({
    prp <- df[, c("prot_acc", "pep_seq")] %>% 
      tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>% 
      dplyr::filter(!duplicated(pep_prot.)) %>% 
      dplyr::select(-pep_prot.)
    
    mat <- model.matrix(~ 0 + prot_acc, prp)
    colnames(mat) <- gsub("prot_acc", "", colnames(mat))
    rownames(mat) <- prp$pep_seq
    
    # collapse rows of the same pep_seq; may use `sum`
    # 
    #      pep_seq prot_acc
    # 1       A        X
    # 2       A        Y
    # 3       B        X
    # 4       C        Y
    # 
    #   X Y
    # 1 1 0
    # 2 0 1
    # 3 1 0
    # 4 0 1
    # 
    # A tibble: 3 x 3
    # pep_seq X     Y    
    # <chr>   <lgl> <lgl>
    # 1 A       TRUE  TRUE 
    # 2 B       TRUE  FALSE
    # 3 C       FALSE TRUE 
    
    mat <- mat %>% 
      data.frame(check.names = FALSE) %>% 
      dplyr::mutate_all(as.logical) %>% 
      dplyr::mutate(pep_seq = prp$pep_seq) %>% 
      dplyr::arrange(pep_seq) %>% 
      dplyr::group_by(pep_seq) %>% 
      dplyr::summarise_all(any) %>% 
      tibble::column_to_rownames("pep_seq")
    
    saveRDS(mat, file.path(out_path, "prot_pep_map.rds"))
    
    invisible(mat)
  })

  # --- set aside df0 ---
  sets <- purrr::quietly(RcppGreedySetCover::greedySetCover)(
    df[, c("prot_acc", "pep_seq")], FALSE)$result
  
  df <- df %>% dplyr::mutate(prot_isess = prot_acc %in% sets$prot_acc) 
  df0 <- df %>% filter(!prot_isess)
  df1 <- df %>% filter(prot_isess)

  mat_ess <- mat[colnames(mat) %in% df1$prot_acc]
  
  peps_uniq <- local({
    rsums <- rowSums(mat)
    rsums2 <- rowSums(mat_ess)
    
    peps <- data.frame(pep_seq = rownames(mat)) %>% 
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>% 
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L)) 
  })
  
  # --- protein ~ protein distance map
  pep_seqs <- rownames(mat_ess)
  mat_ess <- as.list(mat_ess)
  len <- length(mat_ess)
  
  if (len <= 200) {
    out <- vector("list", len)
    
    for (i in seq_len(len)) {
      out[[i]] <- map_dbl(mat_ess[i:len], ~ sum(.x & mat_ess[[i]]))
      out[[i]] <- c(out[seq_len(i-1)] %>% map_dbl(`[[`, i), out[[i]])
    }
  } else {
    out <- parDist(mat_ess)
  }
  
  out <- do.call(rbind, out) 
  rownames(out) <- colnames(out)
  
  stopifnot(identical(out, t(out)))
  
  saveRDS(out, file.path(out_path, "prot_dist.rds"))
  
  # --- finds protein groups
  out[out == 0L] <- 1000000
  out[out < 1000000] <- 0
  out <- out %>% as.dist(diag = TRUE, upper = TRUE)
  
  hc <- hclust(out, method = "single")
  
  grps <- data.frame(prot_hit_num = cutree(hc, h = 1)) %>% 
    tibble::rownames_to_column("prot_acc") %>% 
    dplyr::group_by(prot_hit_num) %>% 
    dplyr::mutate(prot_family_member = row_number()) %>% 
    dplyr::ungroup()

  # --- put together
  df0 <- df0 %>% 
    dplyr::mutate(prot_hit_num = NA, prot_family_member = NA)

  df1 <- df1 %>% 
    dplyr::left_join(grps, by = "prot_acc") %>% 
    dplyr::bind_rows(df0) %>% 
    dplyr::left_join(peps_uniq, by = "pep_seq")

  invisible(df1)
}


#' Helper of \link{parDist}.
#' 
#' @param cols The column indexes of \code{mat}.
#' @inheritParams parDist
par_dist <- function (cols, mat) {
  len_m <- length(mat)
  len_c <- length(cols)
  
  stopifnot(len_c >= 1L)
  
  out <- vector("list", len_c)
  
  for (i in 1:len_c) {
    col <- cols[i]
    out[[i]] <- purrr::map_dbl(mat[col:len_m], ~ sum(.x & mat[[col]]))
  }
  
  invisible(out)
}


#' Parallel distance calculations.
#' 
#' @param mat A bool matrix.
parDist <- function (mat) {
  n_cores <- detectCores()
  idxes <- chunksplit(seq_along(mat), 2 * n_cores)
  
  cl <- makeCluster(getOption("cl.cores", n_cores))
  clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
  out <- clusterApply(cl, idxes, par_dist, mat) 
  stopCluster(cl)
  
  out <- out %>% purrr::flatten()
  
  len <- length(mat)
  
  if (len > 1L) {
    for (i in 2:len) {
      out[[i]] <- c(out[1:(i-1)] %>% map_dbl(`[[`, i), out[[i]])
    }
  }

  invisible(out)
}

