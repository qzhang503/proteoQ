#' Helper of adding apex retention times for all PSM tables
#'
#' Inputs are at the levels of all TMTSet_LCMSInj.txt files.
#' 
#' @param filelist The file names of PSM tables, e.g.,
#'   TMTSet1_LCMSinj1_PSM_N.txt etc.
#' @param ms1full_files The file names of \code{ms1full_.rds}.
#' @param dat_dir The working directory.
#' @param path_ms1 The file path to \code{ms1full_.rds}.
#' @param max_n_apexes The maximum number of apexes for consideration.
#' @param data_type The type of MS data.
haddApexRTs_allsets <- function (filelist = NULL, ms1full_files = NULL, 
                                 dat_dir = NULL, path_ms1 = NULL, 
                                 max_n_apexes = 2L, data_type = "raw")
{
  if (data_type == "raw") {
    yco   <- 100
    min_y <- 2E6
    step  <- 1E-5
  }
  else if (data_type == "pasef") {
    yco   <- 10
    min_y <- 500
    step  <- 1.5E-5
  }
  else {
    yco   <- 100
    min_y <- 2E6
    step  <- 1E-5
  }
  
  n_cores <- min(parallel::detectCores(), length(filelist))
  
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
    ans <- parallel::clusterMap(
      cl, haddApexRTs_oneset, file_psm = filelist, files_ms1 = ms1full_files, 
      MoreArgs = list(dat_dir = dat_dir, path_ms1 = path_ms1, yco = yco, 
                      min_y = min_y, step = step, n_cores = n_cores), 
      SIMPLIFY = FALSE, USE.NAMES = TRUE)
    parallel::stopCluster(cl)
  }
  else {
    ans <- mapply(
      haddApexRTs_oneset, file_psm = filelist, files_ms1 = ms1full_files, 
      MoreArgs = list(dat_dir = dat_dir, path_ms1 = path_ms1, yco = yco, 
                      min_y = min_y, step = step, n_cores = 1L),
      SIMPLIFY = FALSE, USE.NAMES = TRUE)
  }
  
  ans
}


#' Helper of adding apex retention times for one PSM table
#' 
#' Inputs are at the levels of one TMTSet_LCMSInj.txt file.
#'
#' @param file_psm The file names of a PSM table, e.g.,
#'   TMTSet1_LCMSinj1_PSM_N.txt.
#' @param files_ms1 The file name of an \code{ms1full_.rds}.
#' @param dat_dir The working directory.
#' @param path_ms1 The file path to \code{ms1full_.rds}.
#' @param max_n_apexes The maximum number of apexes for consideration.
#' @param n_cores The number of parallel processes calling the current utility.
#' @param from A starting mass.
#' @param step A step size.
#' @param rt_size The width of each LC retention times in seconds.
#' @param yco The cut-off in intensities.
#' @param min_y The minimum peak area.
haddApexRTs_oneset <- function (file_psm = NULL, files_ms1 = NULL, 
                                dat_dir = NULL, path_ms1 = NULL, 
                                max_n_apexes = 2L, from = 115, step = 1E-5, 
                                rt_size = 180, yco = 100, min_y = 2E6, 
                                n_cores = 1L)
{
  df <- readr::read_tsv(
    file.path(dat_dir, "PSM", file_psm), col_types = get_col_types(), 
    show_col_types = FALSE) |> 
    suppressWarnings()
  df <- df[with(df, order(pep_scan_num)), ]
  df$id. <- seq_len(nrow(df))
  
  dfs  <- split(df[, c("pep_exp_mz", "pep_ret_range", "id.")], df$raw_file)
  raws <- gsub("^ms1full_(.*)\\.rds$", "\\1", files_ms1)
  raws <- gsub("\\.(d|raw)$", "", raws)
  files_ms1 <- files_ms1[match(names(dfs), raws)]
  
  n_pcs   <- parallel::detectCores() - 1L
  n_para  <- min(floor(n_pcs / n_cores), length(dfs))
  n_para2 <- min(n_pcs / (n_cores * n_para), 6L)
  
  if (n_para > 1L) { # multiple RAWs under a set
    cl <- parallel::makeCluster(getOption("cl.cores", n_para))
    ans <- parallel::clusterMap(
      cl, addApexRTs, df = dfs, file_ms1 = files_ms1, 
      MoreArgs = 
        list(from = from, rt_size = rt_size, step = step, path_ms1 = path_ms1, 
             yco = yco, min_y = min_y, n_para = n_para2), 
      SIMPLIFY = FALSE, USE.NAMES = TRUE)
    parallel::stopCluster(cl)
  }
  else {
    ans <- mapply(
      addApexRTs, df = dfs, file_ms1 = files_ms1, 
      MoreArgs = 
        list(from = from, rt_size = rt_size, step = step, path_ms1 = path_ms1, 
             yco = yco, min_y = min_y, n_para = n_para2), 
      SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  
  dfs <- mapply(function (x, y) dplyr::bind_cols(x, y), dfs, ans, 
                SIMPLIFY = FALSE, USE.NAMES = FALSE)
  dfs <- lapply(dfs, function (df) {
    df$pep_exp_mz <- df$pep_ret_range <- NULL
    df
  })
  
  out <- dplyr::left_join(df, dplyr::bind_rows(dfs), by = "id.")
  out$id. <- NULL
  rows <- which(out$apex_p == 0L)
  out$apex_p[rows] <- out$apex_n[rows] <- NA_integer_
  out$apex_t[rows] <- out$apex_y[rows] <- NA_real_
  
  out$pep_tot_int <- out$apex_y
  out$pep_apex_ret <- out$apex_t
  out$pep_apex_scan <- out$apex_p
  out$N_I000 <- out$I000 <- out$pep_tot_int
  out$apex_y <- out$apex_t <- out$apex_p <- NULL
  
  out
}


#' Add apex retention times for one RAW
#'
#' Inputs are at the levels of one RAW file.
#'
#' @param df A subset of PSM table under one RAW.
#' @param file_ms1 The file name of an \code{ms1full_.rds}.
#' @param path_ms1 The file path to \code{ms1full_.rds}.
#' @param from The starting point for calculating the bin indexes of masses.
#' @param step A step size in mass for binning.
#' @param max_n_apexes The maximum number of apexes for consideration.
#' @param yco The cut-off in intensities.
#' @param min_y The minimum peak area.
#' @param rt_size The width of each LC retention times in seconds.
#' @param rt_margin The bracketing margin before and after an LC retention time
#'   window.
#' @param n_para The number of parallel processes.
addApexRTs <- function (df, file_ms1, path_ms1, from = 115, step = 1E-5, 
                        rt_size = 180, rt_margin = 120, yco = 100, min_y = 2E6, 
                        n_para = 6L)
{
  # low redundancy: ELEEIRK ~ qIqELRK; TENNDHINLK ~ TENNnHINLK
  # df$uid <- with(df, paste0(raw_file, "@", pep_scan_num ))
  # bads <- df$uid[duplicated(df$uid)]
  # dfx <- df[df$uid %in% bads, ]
  
  ## (1) subset the  MS1 trace space by PSM results
  trs <- qs::qread(file.path(path_ms1, file_ms1))
  tss <- trs$ret_time
  xss <- trs$msx_moverzs
  yss <- trs$msx_ints
  sss <- trs$scan_num
  
  rts <- df$pep_ret_range
  mzs <- df$pep_exp_mz
  unv <- unique(mzion:::index_mz(unique(mzs), from = 115, d = step))
  
  xys <- mapply(function (xs, ys) {
    ixs <- mzion:::index_mz(xs, from = 115, d = step)
    # allow duplicated(ixs)
    oks <- ixs %fin% unv | (ixs + 1L) %fin% unv | (ixs - 1L) %fin% unv
    list(xs = xs[oks], ys = ys[oks])
  }, xss, yss, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  xss  <- lapply(xys, `[[`, "xs")
  yss  <- lapply(xys, `[[`, "ys")
  trs$msx_moverzs <- xss
  trs$msx_ints <- yss
  rm(list = c("unv", "xys"))
  
  ## (2) split xss and yss into chunks with bracketed scans
  if (n_para > 1L) {
    # gap = mzion:::estimate_rtgap(tss, d = 180)
    df1s <- mzion:::sep_df1_byRTs(
      df1 = trs, col_rt = "ret_time", rt_size = rt_size, rt_margin = rt_margin)
    gaps_bf  <- attr(df1s, "gaps_bf",  exact = TRUE)
    gaps_af  <- attr(df1s, "gaps_af",  exact = TRUE)
    min_rts  <- attr(df1s, "min_rts",  exact = TRUE)
    max_rts  <- attr(df1s, "max_rts",  exact = TRUE)
    end1s    <- attr(df1s, "end1s",    exact = TRUE)
    sta1s    <- attr(df1s, "sta1s",    exact = TRUE)
    n_chunks <- attr(df1s, "n_chunks", exact = TRUE)
    
    dfs <- vector("list", n_chunks)
    for (i in 1:n_chunks) {
      # some MS2-RTs are dangling between max_rts[[i-1]] and min_rts[[i]]
      if (i == 1L) {
        dfs[[i]] <- df |> 
          dplyr::filter(pep_ret_range <= max_rts[[i]])
      }
      else if (i == n_chunks) {
        dfs[[i]] <- df |> 
          dplyr::filter(pep_ret_range > max_rts[[i-1L]])
      }
      else {
        dfs[[i]] <- df |> 
          dplyr::filter(pep_ret_range >  max_rts[[i-1L]], 
                        pep_ret_range <= max_rts[[i]])
      }
    }
    
    cl <- parallel::makeCluster(getOption("cl.cores", n_para))
    ans <- parallel::clusterMap(
      cl, addApexes, 
      rts = lapply(dfs, `[[`, "pep_ret_range"), 
      mzs = lapply(dfs, `[[`, "pep_exp_mz"), 
      trs = df1s, gap_bf = gaps_bf, gap_af = gaps_af, 
      MoreArgs = list(step = step, yco = yco, min_y = min_y), 
      SIMPLIFY = FALSE, USE.NAMES = TRUE)
    parallel::stopCluster(cl)
    ans <- dplyr::bind_rows(ans)
  }
  else {
    ans <- addApexes(rts = rts, mzs = mzs, trs = trs, gap_bf = rt_size, 
                     gap_af = rt_size, step = step, yco = yco, min_y = min_y)
  }
  
  ans
}


#' Add retention time apex for a chunk of PSMs under a single RAW
#'
#' For parallel processing.
#'
#' @param rts A vector of retention times at MS2 events (pep_ret_range).
#' @param mzs A vector of m-over-z values for PSMs.
#' @param gap_bf A leading gap of maximum allowance in retention times between
#'   an MS2 scan and a peak apex.
#' @param gap_af A following gap of maximum allowance in retention times between
#'   an MS2 scan and a peak apex.
#' @param step A step size in mass for mass error tolerance.
#' @param yco The cut-off in intensities.
#' @param min_y The minimum peak area.
addApexes <- function(rts, mzs, trs, gap_bf, gap_af, step = 1E-5, yco = 100, 
                      min_y = 2E6)
{
  if (!(nr <- length(rts))) {
    return(NULL)
  }
  
  # stopifnot(nr == length(mzs))
  
  tss <- trs$ret_time
  xss <- trs$msx_moverzs
  yss <- trs$msx_ints
  sss <- trs$scan_num
  
  # a second chance for small percents of intensity interpreted (partial peaks)
  step2 <- 1.5 * step
  y1s   <- t1s <- vector("numeric", nr)
  ap1s  <- n1s <- vector("integer", nr)
  scans <- tvals <- xvals <- yvals <- ns <- ranges <- fwhms <- vector("list", nr)
  
  for (i in seq_len(nr)) {
    # i <- which(abs(mzs - 459.2624) < .0001)[[1]]
    xref <- mzs[[i]] # pep_exp_mz
    tref <- rts[[i]] # pep_ret_range
    rgi  <- .Internal(which(tss >= (tref - gap_bf) & tss <= (tref + gap_af)))
    if (!length(rgi)) { next }
    xsi  <- xss[rgi]
    ysi  <- yss[rgi]
    xys  <- extract_mbry(xs = xsi, ys = ysi, mbr_mz = xref, step = step)
    xsi  <- xys[["x"]]
    ysi  <- xys[["y"]]
    if (all(is.na(ysi))) { next }
    ssi  <- sss[rgi]
    tsi  <- tss[rgi]
    
    gates <- mzion:::find_lc_gates(
      xs = xsi, ys = ysi, ts = tsi, yco = yco, ytot_co = 2E5, n_dia_scans = 6L)
    if (is.null(gates)) { next }
    rows  <- gates[["apex"]]
    sci   <- ssi[rows]
    rti   <- tsi[rows]
    yints <- gates[["yints"]]
    rngs  <- gates[["ranges"]]
    fwhmi <- gates[["fwhm"]]
    # ns[[i]] <- gates[["ns"]]

    # oks <- !is.na(ysi)
    # csum <- cumsum(ysi[oks] * c(diff(tsi[oks]), .5))
    # ysum1[[i]] <- sum(yints)
    # csum[[length(csum)]] # area
    # 
    # ysum0[[i]] <- sum(ysi, na.rm = TRUE)
    # ysum1[[i]] <- sum(ysi[unlist(rngs, recursive = FALSE, use.names = FALSE)], na.rm = TRUE)
    
    if (FALSE && 
        sum(ysi[unlist(rngs)], na.rm = TRUE) / sum(ysi, na.rm = TRUE) < .02) {
      xys  <- extract_mbry(xs = xsi, ys = ysi, mbr_mz = xref, step = step2)
      xsi  <- xys[["x"]]
      ysi  <- xys[["y"]]
      ssi  <- sss[rgi]
      tsi  <- tss[rgi]
      gates <- mzion:::find_lc_gates(
        xs = xsi, ys = ysi, ts = tsi, yco = yco, ytot_co = 2E5, n_dia_scans = 6L)
      if (is.null(gates)) { next }
      
      rows  <- gates[["apex"]]
      sci   <- ssi[rows]
      rti   <- tsi[rows]
      yints <- gates[["yints"]]
      rngs  <- gates[["ranges"]]
      fwhmi <- gates[["fwhm"]]
    }

    anx <- mzion:::find_best_apex(
      xvs = xsi, yvs = ysi, ss = ssi, 
      yints = yints, aps = sci, rts = rti, rngs = rngs, fwhms = fwhmi, 
      xref = xref, ret_ms1 = tref, step_tr = step, min_y = min_y, 
      max_rt_delta = 180, min_n1 = 10L, min_n2 = 20L, min_n3 = 15L)
    
    if (is.null(anx)) {
      next
    }
    
    ap1s[[i]] <- anx[["ap"]]   # apex_scan
    y1s[[i]]  <- anx[["y"]]    # apex_int
    t1s[[i]]  <- anx[["rt"]]   # apex_rt
    n1s[[i]]  <- anx[["n"]]    # number of MS1 scans
    scans[[i]] <- anx[["aps"]] # all apex scans
    tvals[[i]] <- anx[["rts"]] # all apex retention times
    ns[[i]]    <- anx[["ns"]]  # numbers of MS1 scans
    xvals[[i]] <- anx[["xs"]]  # all X values
    yvals[[i]] <- anx[["ys"]]  # all Y values
    fwhms[[i]] <- anx[["fwhms"]]
    # ranges[[i]] <- anx[["rngs"]] # indexes along tsi
    # tsi[anx$rngs[[3]]] # retention-time range for MBR; 
  }
  
  ans <- tibble::tibble(
    apex_p     = ap1s,
    apex_t     = t1s,
    apex_y     = y1s,
    apex_n     = n1s, 
    apex_ps    = scans,
    apex_ts    = tvals,
    apex_xs    = xvals, 
    apex_ys    = yvals,
    # differ by <= 20% for MBR
    apex_fwhms = fwhms, 
    apex_ns    = ns)
  
  if (nrow(ans) != nr) {
    stop("Developer: check for row dropping.")
  }
  
  ans
}


