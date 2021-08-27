#' Reporter-ion quantitation.
#'
#' @param data An upstream result from \link{matchMS}.
#' @param quant A quantitation method. The default is "none". Additional choices
#'   include \code{tmt6} etc. For other multiplicities of \code{tmt}, use the
#'   compatible higher plexes, for example, \code{tmt16} for \code{tmt12} etc.
#'   and \code{tmt10} for \code{tmt8} etc.
#' @param ppm_reporters The mass tolerance of MS2 reporter ions.
calc_tmtint <- function (data = NULL,
                         quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16"),
                         ppm_reporters = 10) {

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
  
  rm(list = c("theopeps", "theopeps_rev"))
  gc()

  # adds prot_n_psm, prot_n_pep for protein FDR
  x <- out %>% 
    dplyr::filter(pep_issig)
  
  prot_n_psm <- x %>%
    dplyr::select(prot_acc) %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::summarise(prot_n_psm = n())
  
  prot_n_pep <- x %>%
    dplyr::select(pep_seq, prot_acc) %>% 
    tidyr::unite(pep_prot, c("pep_seq", "prot_acc"), sep = ".", remove = FALSE) %>% 
    dplyr::filter(!duplicated(pep_prot)) %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::summarise(prot_n_pep = n())

  out <- list(out, prot_n_psm, prot_n_pep) %>%
    purrr::reduce(dplyr::left_join, by = "prot_acc") %>%
    dplyr::arrange(-prot_n_pep, -prot_n_psm)
  
  rm(list = c("x", "prot_n_psm", "prot_n_pep"))
  gc()
  
  invisible(out)
}


#' Helper of \link{groupProts}.
#'
#' Builds the logical map between peptide (in rows) and proteins (in columns).
#'
#' @param df The data frame from upstream steps. It must contains the two
#'   columns of \code{prot_acc} and \code{pep_seq}.
#' @param out_path An output path.
#' @examples 
#' \donttest{
#' df <- data.frame(prot_acc = character(2000), pep_seq = character(2000))
#' set.seed(100)
#' df$prot_acc <- sample(LETTERS[1:20], 2000, replace = TRUE)
#' df$pep_seq <- sample(letters[1:26], 20, replace = TRUE)
#' df <- df[!duplicated(df), ] 
#' 
#' out <- map_pepprot(df)
#' }
map_pepprot <- function (df, out_path = NULL) {
  
  message("Building protein-peptide maps.")
  
  df <- df[, c("prot_acc", "pep_seq")] %>%
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>%
    dplyr::filter(!duplicated(pep_prot.)) %>%
    dplyr::select(-pep_prot.)
  gc()
  
  # ---
  nrow <- nrow(df)
  max_rows <- 30720L
  n_cores <- detect_cores()

  if (nrow > max_rows) {
    prots <- data.frame(prot_acc = unique(df$prot_acc))
    
    n_chunks <- ceiling(nrow/max_rows)
    df <- chunk_groupsplit(df, df$prot_acc, n_chunks)
    out <- vector("list", n_chunks)
    gc()
    
    for (i in seq_along(df)) {
      x <- df[[i]]
      x <- suppressMessages(x %>% dplyr::right_join(prots))
      gc()
      
      mat <- model.matrix(~ 0 + prot_acc, x)
      colnames(mat) <- gsub("prot_acc", "", colnames(mat))
      gc()
      
      peps <- x$pep_seq
      rows <- !is.na(peps)
      peps <- peps[rows]
      mat <- mat[rows, ]
      rownames(mat) <- peps
      gc()
      
      mat <- mat == 1L
      rm(list = "x")
      gc()
      
      mat <- mat %>%
        data.frame(check.names = FALSE) %>%
        dplyr::mutate(pep_seq = peps) %>%
        dplyr::group_by(pep_seq)
      gc()
      
      cl <- multidplyr::new_cluster(n_cores)
      
      out[[i]] <- mat %>%
        multidplyr::partition(cl) %>%
        dplyr::summarise_all(any) %>%
        dplyr::collect() 
      gc()
      
      out[[i]] <- out[[i]] %>%
        tibble::column_to_rownames("pep_seq")
      rm(list = c("cl", "mat"))
      gc()
    }
    
    out <- out %>% 
      dplyr::bind_rows() %>% 
      .[order(rownames(.)), ]
    
  } else {
    mat <- model.matrix(~ 0 + prot_acc, df)
    colnames(mat) <- gsub("prot_acc", "", colnames(mat))
    rownames(mat) <- df$pep_seq
    mat <- mat == 1L

    mat <- mat %>%
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(pep_seq = df$pep_seq) %>%
      dplyr::group_by(pep_seq)
    
    n_cores <- detect_cores()
    cl <- multidplyr::new_cluster(n_cores)
    
    mat <- mat %>%
      multidplyr::partition(cl) %>%
      dplyr::summarise_all(any) %>%
      dplyr::collect() %>%
      tibble::column_to_rownames("pep_seq") %>%
      .[order(rownames(.)), ]
    
    rm(list = c("cl"))
    gc()
    
    out <- mat
  }
  
  
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

  if (!is.null(out_path)) {
    saveRDS(out, file.path(out_path, "prot_pep_map.rds"))
  }

  invisible(out)
}


#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
groupProts <- function (df, out_path = NULL) {

  # `pep_seq` in `df` are all from target and significant;
  # yet target `pep_seq` can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR

  # --- (1) protein ~ peptide map ---
  mat <- map_pepprot(df[, c("prot_acc", "pep_seq")], out_path)
  gc()

  # --- (2) protein ~ protein groups by distance map ---
  grps <- cut_protgrps(mat, out_path)

  # --- (3) lists of protein ~ peptide maps (by the same hit number) ---
  mats <- split.default(mat, grps$prot_hit_num) %>%
    map(~ .x %>% .[rowSums(.) > 0L, , drop = FALSE]) 
  
  if (!is.null(out_path)) {
    saveRDS(mats, file.path(out_path, "prot_pep_famimaps.rds"))
  }

  # --- (4) set covers by groups ---
  message("Grouping proteins by families.")

  if (length(mats) <= 500L) {
    sets <- map(mats, find_ess_prots)
  } else {
    n_cores <- detect_cores()
    cl <- makeCluster(getOption("cl.cores", n_cores))
    
    clusterExport(cl, list("%>%"), envir = environment(magrittr::`%>%`))
    clusterExport(cl, list("find_ess_prots"), 
                  envir = environment(proteoQ:::find_ess_prots))
    clusterExport(cl, list("greedysetcover"), 
                  envir = environment(proteoQ:::greedysetcover))

    sets <- parLapply(cl, mats, find_ess_prots)

    stopCluster(cl)
    rm(list = c("cl", "n_cores"))
    gc()
  }

  sets <- sets %>%
    bind_rows() 
  
  if (!is.null(out_path)) {
    saveRDS(sets, file.path(out_path, "prot_sets.rds")) 
  }
  
  sets <- sets %>% 
    `[[`("prot_acc") %>%
    unique()

  # --- set aside df0 ---
  df <- df %>% dplyr::mutate(prot_isess = prot_acc %in% sets)
  df0 <- df %>% filter(!prot_isess)
  df1 <- df %>% filter(prot_isess)
  rm(list = c("df"))

  mat_ess <- mat[colnames(mat) %in% df1$prot_acc]

  peps_uniq <- local({
    rsums <- rowSums(mat)
    rsums2 <- rowSums(mat_ess)

    peps <- data.frame(pep_seq = rownames(mat)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })

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
cut_protgrps <- function (mat, out_path = NULL) {

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

  if (!is.null(out_path)) {
    # saveRDS(out, file.path(out_path, "prot_dist.rds"))
  }

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
  
  if (!is.null(out_path)) {
    # saveRDS(grps, file.path(out_path, "prot_grps.rds"))
  }

  stopifnot(identical(grps$prot_acc, cns))

  invisible(grps)
}


#' Helper of \link{parDist}.
#' 
#' R version of parallel distance. Not currently used.
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
  }

  invisible(out)
}


#' Helper of \link{parDist}.
#'
#' R version of parallel distance. Uses a subset of mat for memory efficiency.
#' Not currently used.
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
#' For coupling to par_distC; not currently used.
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
  n_cores <- floor(min(mem/size, detect_cores()))
  
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
  rm(list = c("mat"))
  gc()
  
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

  ## assume: (1) two columns 
  # stopifnot(ncol(df) == 2L)

  len <- length(unique(df[[1]]))

  if (len == 1L) {
    return(df)
  }

  nms <- colnames(df)
  colnames(df) <- c("s", "a")

  # (2) no duplicated entries (for speeds)
  # df <- df %>%
  #   tidyr::unite(sa, s, a, sep = "@", remove = FALSE) %>%
  #   dplyr::filter(!duplicated(sa)) %>%
  #   dplyr::select(-sa)

  # ---
  cts <- df %>%
    group_by(s) %>%
    summarise(n = n()) 

  df <- left_join(cts, df, by = "s") %>%
    dplyr::arrange(-n)

  sets <- NULL

  while(nrow(df) > 0L) {
    s <- df[[1, "s"]]
    sa <- df[df$s == s, c("s", "a")]

    sets <- rbind(sets, sa)

    # may consider partial sorting
    df <- df %>%
      filter(! a %in% sa[["a"]]) %>%
      group_by(s) %>%
      mutate(n = n()) %>%
      dplyr::arrange(-n)
  }

  colnames(sets) <- nms

  invisible(sets)
}


#' Greedy set cover.
#' 
#' A matrix input. Output essential sets only (no elements).
#' 
#' @param mat A matrix of protein (cols)-peptide (rows) map.
#' @return A list of proteins. 
greedysetcover2 <- function (mat) {
  
  sets <- NULL
  
  while(!is.null(dim(mat))) {
    max <- which.max(Matrix::colSums(mat))
    
    if (max == 0L) break
    
    prot <- colnames(mat)[max]
    sets <- c(sets, prot)
    
    rows <- mat[, max]
    mat <- mat[!rows, -max]
  }
  
  invisible(sets)
}


#' Greedy set cover.
#'
#' Expands a two-column input to a matrix input. The output table contains both
#' essential sets and elements.
#'
#' @param mat A matrix of protein (cols)-peptide (rows) map.
#' @return A two-column data frame of prot_acc and pep_seq.
greedysetcoverM <- function (df) {
  
  nms <- colnames(df)
  colnames(df) <- c("s", "a")
  
  sa <- df[, c("s", "a")] %>%
    tidyr::unite(sa., s, a, sep = "@", remove = FALSE) %>%
    dplyr::filter(!duplicated(sa.)) %>%
    dplyr::select(-sa.)
  
  mat <- model.matrix(~ 0 + s, sa)
  colnames(mat) <- gsub("s", "", colnames(mat))
  rownames(mat) <- sa$a
  
  mat <- mat %>%
    data.frame(check.names = FALSE) %>%
    dplyr::mutate_all(as.logical) %>%
    dplyr::mutate(a = sa$a) %>%
    dplyr::arrange(a) %>%
    dplyr::group_by(a)

  # collapse rows
  n_cores <- detect_cores()
  cl <- multidplyr::new_cluster(n_cores)
  
  mat <- mat %>%
    multidplyr::partition(cl) %>%
    dplyr::summarise_all(any) %>%
    dplyr::collect() %>%
    tibble::column_to_rownames("a") %>%
    .[order(rownames(.)), ]
  
  rm(list = c("cl"))
  gc()
  
  # ---
  mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  
  s_out <- NULL
  a_out <- NULL
  
  while(nrow(mat)) {
    max <- which.max(Matrix::colSums(mat))
    
    if (max == 0L) break
    
    s <- names(max)
    rows <- which(mat[, max] == 1L)
    a <- names(rows)

    s_out <- c(s_out, rep(s, length(a)))
    a_out <- c(a_out, a)
    
    mat <- mat[-rows, -max, drop = FALSE]
  }
  
  gc()
  
  dplyr::bind_cols(!!nms[1] := s_out, !!nms[2] := a_out)
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


