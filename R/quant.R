#' Reporter-ion quantitation.
#'
#' @param data An upstream result from \link{matchMS}.
#' @param quant A quantitation method. The default is "none".
#' @param ppm_reporters The mass tolerance of MS2 reporter ions.
calc_tmtint <- function (data = NULL,
                         quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"),
                         ppm_reporters = 10) {

  # val <- rlang::enexpr(quant)
  # f <- match.call()[[1]]
  # val <- match_valexpr(f = !!f, arg = "quant", val = !!val)

  if (quant == "none") {
    out <- data
  } else {
    # message("Calculating reporter-ion intensities.")

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

    theos <- switch(quant,
                    tmt6 = tmts %>% .[names(.) %in% nms_tmt6],
                    tmt10 = tmts %>% .[names(.) %in% nms_tmt10],
                    tmt11 = tmts %>% .[names(.) %in% nms_tmt11],
                    tmt16 = tmts %>% .[names(.) %in% nms_tmt16],
                    stop("Unknown TMt type.", call. = FALSE))

    ul <- switch(quant,
                 tmt6 = c(126.1, 131.2),
                 tmt10 = c(126.1, 131.2),
                 tmt11 = c(126.1, 131.2),
                 tmt16 = c(126.1, 134.2),
                 stop("Unknown TMt type.", call. = FALSE))

    stopifnot(all(c("ms2_moverz", "ms2_int") %in% names(data)))

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

  names(out)[grep("^([0-9]{3}[NC]{0,1})", names(out))] <-
    find_int_cols(length(theos))

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
#'
#' # Two `129C`, no `127N` etc.
#' ms2_moverzs <- c(105.1503, 107.0428, 111.7716, 120.0811, 126.1281, 127.1312,
#'                  128.1282, 128.1349, 129.1317, 129.1365, 129.1382, 230.1694,
#'                  233.4857, 233.4964, 337.3533, 352.1844, 376.2764, 463.3083,
#'                  525.2150, 562.3732, 569.3899, 591.2545, 596.0308, 632.3300,
#'                  636.3959, 703.3637, 789.0423, 816.4487, 817.4516, 839.9531,
#'                  864.3056, 914.7645, 921.5302, 1479.9816)
#'
#' ms2_ints <- c(1201.79, 1319.32, 1603.45, 1595.34, 2148.66, 1785.74, 1254.24,
#'               1986.43, 10127.40, 1522.60, 1562.71, 2926.01, 1590.48, 1692.17,
#'               1347.88, 1412.64, 3050.10, 3231.10, 1355.21, 2424.18, 1783.26,
#'               1365.32, 1727.12, 2661.72, 1660.05, 5525.95, 1399.96, 4654.03,
#'               1990.57, 1758.72, 1655.09, 1460.68, 1641.39, 1721.33)
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
      .[nms]
  }

  # missing channels:
  # 126 <NA> 127C 128N 128C 129N 129C <NA> <NA> <NA>
  #  2   NA    3    4    5    6    8   NA   NA   NA

  if (anyNA(names(idxes))) {
    names(idxes) <- nms
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
  # Targets, theoretical
  bins <- list.files(path = file.path(.path_fasta, "pepmasses", .time_stamp),
                     pattern = "binned_theopeps_\\d+\\.rds$",
                     full.names = TRUE)

  theopeps <- purrr::map(bins, ~ {
    x <- readRDS(.x) %>%
      dplyr::bind_rows() %>%
      dplyr::select(c("prot_acc", "pep_seq"))
  }) %>%
    dplyr::bind_rows() %>%
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>%
    dplyr::filter(!duplicated(pep_prot.)) %>%
    dplyr::select(-pep_prot.) %>%
    dplyr::select(c("prot_acc", "pep_seq")) %>%
    dplyr::filter(pep_seq %in% unique(df$pep_seq))

  # Decoys, theoretical
  bins_rev <- list.files(path = file.path(.path_fasta, "pepmasses", .time_stamp),
                         pattern = "binned_theopeps_rev_\\d+\\.rds$",
                         full.names = TRUE)

  theopeps_rev <- purrr::map(bins_rev, ~ {
    x <- readRDS(.x) %>%
      dplyr::bind_rows() %>%
      dplyr::select(c("prot_acc", "pep_seq"))
  }) %>%
    dplyr::bind_rows() %>%
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>%
    dplyr::filter(!duplicated(pep_prot.)) %>%
    dplyr::select(-pep_prot.) %>%
    dplyr::select(c("prot_acc", "pep_seq")) %>%
    dplyr::filter(pep_seq %in% unique(df$pep_seq))

  # adds `prot_acc` (with decoys being kept)
  out <- bind_rows(theopeps, theopeps_rev) %>%
    dplyr::right_join(df, by = "pep_seq")

  # adds prot_n_psm, prot_n_pep for protein FDR
  prot_n_psm <- out %>%
    dplyr::select(prot_acc) %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::summarise(prot_n_psm = n())

  # keep duplicated pep_seq; otherwise, some NA prot_acc after joining
  prot_n_pep <- out %>%
    dplyr::select(pep_seq, prot_acc) %>%
    # dplyr::filter(!duplicated(pep_seq)) %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::summarise(prot_n_pep = n())

  out <- list(out, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = "prot_acc") %>%
    dplyr::arrange(-prot_n_pep, -prot_n_psm)
}




#' Helper of \link{groupProts}.
#'
#' @param out The data frame from upstream steps.
#' @param out_path The output path.
grp_prots <- function (out, out_path) {
  out <- out %>%
    dplyr::arrange(pep_seq)

  # essentials
  # (not to add filter `prot_issig`)
  rows <- (out$pep_issig & (!out$pep_isdecoy) & (!grepl("^-", out$prot_acc)))

  df <- out[rows, ]

  if (nrow(df) > 1L) {
    df <- df %>% groupProts(out_path)
  } else {
    df <- df %>%
      dplyr::mutate(prot_isess = TRUE,
                    prot_hit_num = 1L,
                    prot_family_member = 1L)
  }

  # non-essentials
  prot_accs <- df %>%
    dplyr::filter(!duplicated(prot_acc), prot_isess) %>%
    `[[`("prot_acc")

  df2 <- out[!rows, ] %>%
    dplyr::filter(!duplicated(prot_acc)) %>%
    dplyr::select(prot_acc) %>%
    dplyr::mutate(prot_isess = ifelse(prot_acc %in% prot_accs, TRUE, FALSE)) %>%
    dplyr::left_join(df[, c("prot_acc", "prot_hit_num", "prot_family_member")] %>%
                       dplyr::filter(!duplicated(prot_acc)),
                     by = "prot_acc") %>%
    dplyr::right_join(out[!rows, ], by = "prot_acc") %>%
    dplyr::mutate(pep_literal_unique = NA, pep_razor_unique = NA) %>%
    dplyr::select(names(df))

  rm(prot_accs)

  out <- dplyr::bind_rows(df, df2) %>%
    dplyr::select(-which(names(.) %in% c("prot_n_psm", "prot_n_pep")))
}


#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
groupProts <- function (df, out_path = "~/proteoQ/outs") {

  # pep_seq in df are all from target and significant;
  # yet target pep_seq can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR

  # --- (1) protein ~ peptide map ---
  mat <- local({
    message("Building protein-peptide maps.")

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
      dplyr::group_by(pep_seq)

    run_scripts <- FALSE
    if (run_scripts) {
      mat <- mat %>%
        dplyr::summarise_all(any) %>%
        tibble::column_to_rownames("pep_seq")
    }

    n_cores <- detectCores()
    cl <- multidplyr::new_cluster(n_cores)

    mat <- mat %>%
      multidplyr::partition(cl) %>%
      dplyr::summarise_all(any) %>%
      dplyr::collect() %>%
      tibble::column_to_rownames("pep_seq") %>%
      .[order(rownames(.)), ]

    rm(cl)
    gc()

    saveRDS(mat, file.path(out_path, "prot_pep_map.rds"))

    invisible(mat)
  })

  # --- (2) protein ~ protein groups by distance map ---
  grps <- cut_protgrps(mat, out_path)

  # --- (3) lists of protein ~ peptide maps (by the same hit number) ---
  mats <- split.default(mat, grps$prot_hit_num) %>%
    map(~ .x %>% .[rowSums(.) > 0L, , drop = FALSE]) %T>%
    saveRDS(file.path(out_path, "prot_pep_famimaps.rds"))

  # --- (4) set covers by groups ---
  message("Grouping proteins by families.")

  if (length(mats) <= 500L) {
    sets <- map(mats, find_ess_prots)
  } else {
    n_cores <- detectCores()
    cl <- makeCluster(getOption("cl.cores", n_cores))
    
    clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
    clusterExport(cl, list("find_ess_prots"), 
                  envir = environment(proteoQ:::find_ess_prots))
    clusterExport(cl, list("greedysetcover"), 
                  envir = environment(proteoQ:::greedysetcover))

    sets <- parLapply(cl, mats, find_ess_prots)

    stopCluster(cl)
    rm(n_cores, cl)
    gc()
  }

  sets <- sets %>%
    bind_rows() %T>%
    saveRDS(file.path(out_path, "prot_sets.rds")) %>%
    `[[`("prot_acc") %>%
    unique()

  # --- set aside df0 ---
  df <- df %>% dplyr::mutate(prot_isess = prot_acc %in% sets)
  df0 <- df %>% filter(!prot_isess)
  df1 <- df %>% filter(prot_isess)
  rm(df)

  mat_ess <- mat[colnames(mat) %in% df1$prot_acc]

  peps_uniq <- local({
    rsums <- rowSums(mat)
    rsums2 <- rowSums(mat_ess)

    peps <- data.frame(pep_seq = rownames(mat)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })

  # rm(mat)

  # --- put together
  df0 <- df0 %>%
    dplyr::mutate(prot_hit_num = NA, prot_family_member = NA)

  df1 <- df1 %>%
    dplyr::left_join(grps, by = "prot_acc") %>%
    dplyr::bind_rows(df0) %>%
    dplyr::left_join(peps_uniq, by = "pep_seq")

  invisible(df1)
}


#' Cuts proteins into groups.
#'
#' By the number of shared peptides.
#'
#' @param mat A logical matrix; peptides in rows and proteins in columns.
#' @param out_path A file pth to outputs.
cut_protgrps <- function (mat, out_path) {

  cns <- colnames(mat)

  mat <- as.list(mat)
  len <- length(mat)

  if (len <= 200L) {
    out <- vector("list", len)

    for (i in seq_len(len)) {
      out[[i]] <- map_dbl(mat[i:len], ~ sum(.x & mat[[i]]))
      out[[i]] <- c(out[seq_len(i-1)] %>% map_dbl(`[[`, i), out[[i]])
    }
  } else {
    out <- parDist(mat)
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
    dplyr::ungroup() %T>%
    saveRDS(file.path(out_path, "prot_grps.rds"))

  stopifnot(identical(grps$prot_acc, cns))

  invisible(grps)
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
    y <- mat[[col]]
    out[[i]] <- purrr::map_int(mat[col:len_m], ~ sum(.x & y))

    # out[[i]] <- purrr::map_int(col:len_m, ~ {
    #   x <- mat[[.x]]
    #   sum(x & y)
    # }) %>%
    #   `names<-`(names(mat[col:len_m]))
  }

  invisible(out)
}


#' Helper of \link{parDist}.
#'
#' @param cols The column indexes of \code{mat}.
#' @inheritParams parDist
#' @examples
#' \donttest{
#' x <- par_dist1(idxes[[i]], mat[idxes[[i]][1]:len])
#' }
par_dist1 <- function (cols, mat) {

  nms <- names(mat)
  len_m <- length(mat)
  len_c <- length(cols)

  stopifnot(len_c >= 1L, len_m >= len_c)

  out <- vector("list", len_c)

  for (i in 1:len_c) {
    y <- mat[[i]]

    js <- i:len_m
    outj <- integer(length(js))

    k <- 1
    for (j in js) {
      outj[k] <- sum(mat[[j]] & y)
      k <- k + 1
    }

    names(outj) <- nms
    out[[i]] <- outj
    nms <- nms[-1]

    if( i %% 100 == 0) gc()
  }

  invisible(out)
}


#' Helper of \link{parDist}.
#'
#' @param cols The column indexes of \code{mat}.
#' @inheritParams parDist
par_dist2 <- function (cols, mat) {

  out <- par_distC(cols, mat)

  nms <- names(mat)[cols[1]:length(mat)]

  out <- out %>%
    map(~ {
      names(.x) <- nms
      nms <<- nms[-1]

      .x
    })
}


#' Parallel distance calculations.
#'
#' @param mat A bool matrix.
#' @import parallel
parDist <- function (mat) {

  message("Calculating distance matrix.")
  
  gc()

  size <- object.size(mat)/1024^3
  mem <- memory.limit() *.45/1024
  n_cores <- floor(min(mem/size, detectCores()))

  if (n_cores <= 1L) {
    stop("Not enough memory for parallel distance calculation.", call. = FALSE)
  }

  idxes <- chunksplit(seq_along(mat), 2 * n_cores, "list")
  len <- length(mat)
  nms <- names(mat)

  cl <- makeCluster(getOption("cl.cores", n_cores))
  clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))


  out <- clusterApplyLB(cl, idxes, proteoCpp::par_distC, mat) %>%
    purrr::flatten()

  stopCluster(cl)
  rm(mat)
  gc()
  
  # out <- out %>%
  #   map(~ {
  #     names(.x) <- nms
  #     nms <<- nms[-1]
  # 
  #     .x
  #   })

  # out <- clusterApplyLB(cl, idxes, par_dist, mat) %>% purrr::flatten()

  if (len > 1L) {
    for (i in 2:len) {
      out[[i]] <- c(out[1:(i-1)] %>% map_dbl(`[[`, i), out[[i]])
    }
  }

  out <- out %>%
    map(~ {
      names(.x) <- nms
      .x
    })

  invisible(out)
}


#' Greedy set cover.
#'
#' @param df A two-column data frame.
greedysetcover <- function (df) {

  stopifnot(ncol(df) == 2L)

  len <- length(unique(df[[1]]))

  if (len == 1L) {
    # assume no duplicated entries (for speeds)
    # df <- df %>% filter(!duplicated(.[[2]]))
    return(df)
  }

  nms <- colnames(df)
  colnames(df) <- c("s", "a")

  df <- df %>%
    tidyr::unite(sa, s, a, sep = "@", remove = FALSE) %>%
    dplyr::filter(!duplicated(sa)) %>%
    dplyr::select(-sa)

  # ---
  cts <- df %>%
    group_by(s) %>%
    summarise(n = n()) %>%
    dplyr::arrange(-n)

  df <- left_join(cts, df, by = "s")

  sets <- NULL

  while(nrow(df) > 0L) {
    s <- df[[1, "s"]]
    sa <- df[df$s == s, c("s", "a")]

    sets <- rbind(sets, sa)

    df <- df %>%
      filter(! a %in% sa[["a"]]) %>%
      group_by(s) %>%
      mutate(n = n()) %>%
      dplyr::arrange(-n)
  }

  colnames(sets) <- nms

  invisible(sets)
}


#' Helper of \link{greedysetcover}.
#'
#' @param mat A logical matrix; protein accession in column names and peptide
#'   sequences in row names.
find_ess_prots <- function (mat) {
  mat %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("pep_seq") %>%
    gather("prot_acc", "presence", -pep_seq) %>%
    filter(presence) %>%
    select(c("prot_acc", "pep_seq")) %>%
    greedysetcover()
}





#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
groupProts_orig <- function (df, out_path = "~/proteoQ/outs") {

  # pep_seq in df are all from target and significant;
  # yet target pep_seq can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR

  # --- protein ~ peptide map ---
  mat <- local({
    message("Building protein-peptide maps.")

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
      dplyr::group_by(pep_seq)

    run_scripts <- FALSE
    if (run_scripts) {
      mat <- mat %>%
        dplyr::summarise_all(any) %>%
        tibble::column_to_rownames("pep_seq")
    }

    n_cores <- detectCores()
    cl <- multidplyr::new_cluster(n_cores)

    mat <- mat %>%
      multidplyr::partition(cl) %>%
      dplyr::summarise_all(any) %>%
      dplyr::collect() %>%
      tibble::column_to_rownames("pep_seq") %>%
      .[order(rownames(.)), ]

    rm(cl)
    gc()

    saveRDS(mat, file.path(out_path, "prot_pep_map.rds"))

    invisible(mat)
  })

  # --- set aside df0 ---
  message("Grouping proteins by families.")

  sets <- df %>%
    dplyr::select(c("prot_acc", "pep_seq")) %>%
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>%
    dplyr::filter(!duplicated(pep_prot.)) %>%
    dplyr::select(-pep_prot.) %>%
    greedysetcover() %T>%
    saveRDS(file.path(out_path, "prot_sets.rds")) %>%
    `[[`("prot_acc") %>%
    unique()

  df <- df %>% dplyr::mutate(prot_isess = prot_acc %in% sets)
  df0 <- df %>% filter(!prot_isess)
  df1 <- df %>% filter(prot_isess)
  rm(df)

  mat_ess <- mat[colnames(mat) %in% df1$prot_acc]

  peps_uniq <- local({
    rsums <- rowSums(mat)
    rsums2 <- rowSums(mat_ess)

    peps <- data.frame(pep_seq = rownames(mat)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })

  # --- protein ~ protein distance map
  run_scripts <- FALSE
  if (run_scripts) {
    pep_seqs <- rownames(mat_ess)
    mat_ess <- as.list(mat_ess)
    len <- length(mat_ess)

    gc()

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
  }

  grps <- cut_protgrps(mat_ess, out_path)

  # --- put together
  df0 <- df0 %>%
    dplyr::mutate(prot_hit_num = NA, prot_family_member = NA)

  df1 <- df1 %>%
    dplyr::left_join(grps, by = "prot_acc") %>%
    dplyr::bind_rows(df0) %>%
    dplyr::left_join(peps_uniq, by = "pep_seq")

  invisible(df1)
}

