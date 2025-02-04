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
  
  out$pep_tot_int   <- out$apex_y
  out$pep_apex_ret  <- out$apex_t
  out$pep_apex_scan <- out$apex_p
  out$pep_apex_n    <- out$apex_n
  out$pep_apex_fwhm <- out$apex_fwhm
  out$N_I000 <- out$I000 <- out$pep_tot_int
  out$apex_y <- out$apex_t <- out$apex_p <- out$apex_n <- out$apex_fwhm <- NULL
  
  out <- out |>
    reloc_col_after("pep_apex_n", "pep_apex_scan") |>
    reloc_col("pep_apex_fwhm", "pep_apex_n")
  
  nms <- names(out)
  out <- dplyr::bind_cols(
    out[grepl("^prot_", nms), drop = FALSE],
    out[grepl("^pep_",  nms), drop = FALSE],
    out[grepl("^psm_",  nms), drop = FALSE],
    out[!grepl("^prot_|^pep_|^psm_", nms), drop = FALSE])
  
  cols <- grepl("[IR]{1}[0-9]{3}[NC]{0,1}", names(out))
  out  <- dplyr::bind_cols(
    out[, !cols, drop = FALSE], 
    out[,  cols, drop = FALSE])

  out
}


#' Helper of adding apex retention times for one RAW
#'
#' Inputs are at the levels of one RAW file.
#'
#' @param df A subset of PSM table under one RAW.
#' @param file_ms1 The file name of an \code{ms1full_.rds}.
#' @param path_ms1 The file path to \code{ms1full_.rds}.
#' @param from The starting point for calculating the bin indexes of masses.
#' @param step A step size in mass for binning.
#' @param yco The cut-off in intensities.
#' @param min_y The minimum peak area.
#' @param rt_size The width of each LC retention times in seconds.
#' @param rt_margin The bracketing margin before and after an LC retention time
#'   window.
#' @param n_para The number of parallel processes.
addApexRTs <- function (df, file_ms1, path_ms1, from = 115, step = 1E-5, 
                        rt_size = 180, rt_margin = 180, 
                        yco = 100, min_y = 2E6, n_para = 6L)
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
  unv <- unique(mzion::index_mz(unique(mzs), from = from, d = step))
  
  xys <- mapply(function (xs, ys) {
    ixs <- mzion::index_mz(xs, from = from, d = step)
    # allow duplicated(ixs)
    oks <- ixs %fin% unv | (ixs + 1L) %fin% unv | (ixs - 1L) %fin% unv
    list(xs = xs[oks], ys = ys[oks])
  }, xss, yss, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  xss  <- trs$msx_moverzs <- lapply(xys, `[[`, "xs")
  yss  <- trs$msx_ints    <- lapply(xys, `[[`, "ys")
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
#' @param trs MS1 trace data.
#' @param gap_bf A leading gap of maximum allowance in retention times between
#'   an MS2 scan and a peak apex.
#' @param gap_af A following gap of maximum allowance in retention times between
#'   an MS2 scan and a peak apex.
#' @param step A step size in mass for mass error tolerance.
#' @param yco The cut-off in intensities.
#' @param min_y The minimum peak area.
addApexes <- function(rts, mzs, trs, gap_bf = 250L, gap_af = 250L, 
                      step = 1E-5, yco = 100, min_y = 2E6)
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
  fw1s  <- y1s <- t1s <- vector("numeric", nr)
  ap1s  <- n1s <- vector("integer", nr)
  scans <- tvals <- xvals <- yvals <- ns <- ranges <- fwhms <- vector("list", nr)
  
  for (i in seq_len(nr)) {
    # i <- which(abs(mzs - 897.4539) < .0001)[[1]]
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

    ## for finding percent of MS1 area interpreted
    # oks <- !is.na(ysi)
    # csum <- cumsum(ysi[oks] * c(diff(tsi[oks]), .5))
    # ysum1[[i]] <- sum(yints)
    # csum[[length(csum)]] # area
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
    
    ap1s[[i]]  <- anx[["ap"]]   # apex_scan
    y1s[[i]]   <- anx[["y"]]    # apex_int
    t1s[[i]]   <- anx[["rt"]]   # apex_rt
    n1s[[i]]   <- anx[["n"]]    # number of MS1 scans
    fw1s[[i]]  <- anx[["fwhm"]]
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
    apex_fwhm  = fw1s, 
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


#' Spreads Mzion peptide numbers.
#'
#' No aggregation of data at the same TMT_Set but different LCMS_Injs.
#'
#' Spreads fields of numeric values: sd_log2_R, log2_R, log2_R, I, N_I by TMT
#' sets.
#'
#' Also works for LFQ as each sample corresponds to a TMT set.
#'
#' For single SILAC sample, the values of log2Ratios spreads into
#' \emph{MULTIPLE} columns of heavy, light etc. Despite, log2Ratios remains NA,
#' just like regular single-sample LFQ. The log2Ratios will be later calculated
#' with \link{calcLFQPepNums} that are based on intensity values.
#'
#' @param ms1files The names of MS1 data files.
#' @param path_ms1 The path to MS1 data files.
#' @param temp_dir The path to temporary files.
#' @param filelist The name of peptide files (TMTset1_LCMSinj1_Peptide_N.txt)
#'   with prepending path.
#' @param basenames The basenames corresponding to \code{filelist}.
#' @param set_idxes The TMT_Set indexes corresponding to \code{basenames}.
#' @param injn_idxes The LCMS_Injection indexes corresponding to
#'   \code{basenames}.
#' @param are_refs A logical vector corresponding to \code{basenames},
#'   indicating the reference status.
#' @param are_smpls A logical vector corresponding to \code{basenames},
#'   indicating the sample status.
#' @param dfs Lists of peptide tables in correspondence to \code{filelist}.
#' @param df_sps A look-up table between peptide sequences and species.
#' @param prot_spec_counts A data frame of protein spectrum counts.
#' @param rt_tol Error tolerance in retention times.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param max_mbr_fold Not used. The maximum absolute fold change in MBR.
#' @param dat_dir The working directory.
#' @param sp_centers_only Logical; for a side-effect to return only the values
#'   of species centers.
#' @param imp_refs Logical; impute missing references or not.
#' @param group_psm_by Group PSMs by.
#' @param group_pep_by Group peptides by.
#' @inheritParams normPSM
hpepLFQ <- function (filelist, basenames, set_idxes, injn_idxes, 
                     are_refs, are_smpls, 
                     dfs = NULL, df_sps = NULL, prot_spec_counts = NULL, 
                     dat_dir = NULL, path_ms1 = NULL, ms1files = NULL, 
                     temp_dir = NULL, rt_tol = 25, mbr_ret_tol = 25, 
                     sp_centers_only = FALSE, imp_refs = FALSE, 
                     group_psm_by = "pep_seq_modz", group_pep_by = "gene", 
                     lfq_mbr = FALSE, max_mbr_fold = 20)
{
  if ((len <- length(basenames)) < 2L) {
    warning("Requires at least two files for MBR.")
    return(NULL)
  }
  
  if (file.exists(fi_datatype <- file.path(path_ms1, "data_type.rds"))) {
    data_type <- qs::qread(fi_datatype)
  }
  else {
    warning("File not found: ", fi_datatype, ". Assume Thermo's data.")
    data_type <- "raw"
  }
  
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
  
  message("Starting LFQ.")
  
  ## (1) group-split of MS1 RAWs corresponding to basenames (peptide tables)
  ms1grps <- sep_mbrfiles(
    b_nms = basenames, dat_dir = dat_dir, 
    ms1files = ms1files, type = "calibms1full")
  
  ## (2) find species centers (based on the tentative apexes)
  list_cols <- 
    c("pep_apex_ps", "pep_apex_ts", "pep_apex_xs", "pep_apex_ys", 
      "pep_apex_fwhms", "pep_apex_ns")
  cols<- 
    c("pep_seq_modz", "pep_tot_int", "pep_apex_ret", "pep_apex_scan", 
      "pep_exp_mz", "pep_apex_fwhm", "pep_apex_n", list_cols)
  
  dfx <- lapply(dfs, function (df) {
    df <- df[, cols]
    # df <- df[with(df, !is.na(pep_tot_int)), ]
    # df <- df[!duplicated(df$pep_seq_modz), ]
    
    df[list_cols] <- lapply(list_cols, function (col) {
      vs <- stringr::str_split(df[[col]], ";")
      
      vs <- if (col == "pep_apex_ns") {
        lapply(vs, as.integer)
      }
      else {
        lapply(vs, as.numeric)
      }
    })
    
    df
  })
  
  mats <- makeLFQXYTS(
    pep_seq_modzs  = lapply(dfx, `[[`, "pep_seq_modz"), 
    pep_exp_mzs    = lapply(dfx, `[[`, "pep_exp_mz"), 
    pep_exp_ints   = lapply(dfx, `[[`, "pep_tot_int"), 
    pep_apex_rets  = lapply(dfx, `[[`, "pep_apex_ret"), 
    pep_apex_scans = lapply(dfx, `[[`, "pep_apex_scan"), 
    pep_apex_fwhm  = lapply(dfx, `[[`, "pep_apex_fwhm"), 
    pep_apex_n     = lapply(dfx, `[[`, "pep_apex_n"), 
    add_colnames   = FALSE)
  ymat <- mats[["y"]]
  xmat <- mats[["x"]]
  tmat <- mats[["t"]]
  smat <- mats[["s"]]
  fmat <- mats[["f"]]
  nmat <- mats[["n"]] # not to filter by n: many MS1s without MS2s -> large n
  
  mbr_peps  <- mats[["unv"]]
  mbr_mzs   <- colMeans(xmat, na.rm = TRUE)
  mbr_ys    <- colMeans(ymat[are_refs, , drop = FALSE], na.rm = TRUE)
  mbr_rets  <- colMeans(tmat, na.rm = TRUE)
  # mbr_fwhms <- colMeans(fmat, na.rm = TRUE)
  # mbr_ns    <- colMeans(nmat, na.rm = TRUE)
  
  mbr_rets <- 
    fillMBRnaRTs(mbr_rets = mbr_rets, mbr_peps = mbr_peps, dat_dir = dat_dir) 

  # should contain no NA matches since df_sps$pep_seq_modz is complete
  if (any(is.na(mts <- match(mbr_peps, df_sps$pep_seq_modz)))) {
    stop("Developer: check for droppings of peptide entries.")
  }
  mbr_sps    <- df_sps[["species"]][mts]
  sp_centers <- find_species_centers2(
    ymat = ymat, are_refs = are_refs, are_smpls = are_smpls, sps = mbr_sps)
  mbr_sps    <- match(mbr_sps, names(sp_centers)) # convert to integers
  
  if (file.exists(parfile <- file.path(temp_dir, "pars_rt.rds"))) {
    fwhm_co <- qs::qread(parfile)$fwhm_co
  }
  else {
    warning("Parameter file not found: ", parfile)
    fwhm_co <- .5
  }
  
  ## (3) optional: make spectrum-count matrix
  message("  Adding spectrum-count results.")
  cmat <- local({
    cmat <- matrix(rep_len(NA_integer_, prod(dim(ymat))), nrow = nrow(ymat))
    colnames(cmat) <- mbr_peps
    rownames(cmat) <- basenames
    dfsc <- prot_spec_counts |> dplyr::left_join(df_sps, by = group_pep_by)
    dfsc <- dfsc[, c(group_psm_by, "TMT_Set", "LCMS_Injection", "prot_n_specs")]
    dfsc <- split(dfsc, with(dfsc, paste0(TMT_Set, "_", LCMS_Injection)))
    
    names(dfsc) <- gsub(
      "(\\d+)\\_(\\d+)", 
      paste0("TMTset", "\\1", "_LCMSinj", "\\2", "_Peptide_N.txt"), 
      names(dfsc))
    dfsc <- dfsc[match(basenames, names(dfsc))]
    dfsc <- lapply(dfsc, `[`, c(group_psm_by, "prot_n_specs"))
    
    for (i in seq_along(dfsc)) {
      dfi <- dfsc[[i]]
      mts <- match(dfi[[group_psm_by]], mbr_peps)
      oks <- !is.na(mts)
      cmat[i, mts[oks]] <- dfi[["prot_n_specs"]][oks]
    }
    
    cmat
  })
  rm(list = c("mats", "mts"))
  
  ## (4) extract local MBR patterns
  message("  Extracting local apex patterns.")
  matss <- makeLFQXYTS(
    pep_seq_modzs  = lapply(dfx, `[[`, "pep_seq_modz"), 
    pep_exp_mzs    = lapply(dfx, `[[`, "pep_apex_xs"), 
    pep_exp_ints   = lapply(dfx, `[[`, "pep_apex_ys"), 
    pep_apex_rets  = lapply(dfx, `[[`, "pep_apex_ts"), 
    pep_apex_scans = lapply(dfx, `[[`, "pep_apex_ps"), 
    pep_apex_fwhm  = lapply(dfx, `[[`, "pep_apex_fwhms"), 
    pep_apex_n     = lapply(dfx, `[[`, "pep_apex_ns"), 
    add_colnames   = FALSE)
  ysmat <- matss[["y"]]
  xsmat <- matss[["x"]]
  tsmat <- matss[["t"]]
  ssmat <- matss[["s"]]
  fsmat <- matss[["f"]]
  nsmat <- matss[["n"]]
  rm(list = c("dfx", "matss"))
  
  ###
  # TODO: each column -> any > species-specific contour -> rematch
  ###
  
  ## (4) MBR
  if (lfq_mbr) {
    message("  Extracting MBR data.")
    
    if ((n_cores <- min(parallel::detectCores(), len)) > 2L) {
      cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
      ans_mbr <- parallel::clusterMap(
        cl, saddMBR, 
        base_name = basenames, 
        ms1files  = file.path(path_ms1, ms1grps), 
        row_id    = seq_along(basenames), 
        MoreArgs  = list(
          xmat = xsmat, ymat = ysmat, tmat = tsmat, smat = ssmat, 
          fmat = fsmat, nmat = nsmat, mbr_mzs = mbr_mzs, mbr_rets = mbr_rets, 
          dat_dir = dat_dir, mbr_ret_tol = mbr_ret_tol, min_y = min_y, 
          yco = yco, fwhm_co = fwhm_co))
      parallel::stopCluster(cl)
    }
    else {
      ans_mbr <- mapply(
        saddMBR, 
        base_name = basenames, 
        ms1files  = file.path(path_ms1, ms1grps), 
        row_id    = seq_along(basenames), 
        MoreArgs  = list(
          xmat = xsmat, ymat = ysmat, tmat = tsmat, smat = ssmat, 
          fmat = fsmat, nmat = nsmat, mbr_mzs = mbr_mzs, mbr_rets = mbr_rets, 
          dat_dir = dat_dir, mbr_ret_tol = mbr_ret_tol, min_y = min_y, 
          yco = yco, fwhm_co = fwhm_co), 
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }

    # xsmat is not changed
    ysmat <- lapply(ans_mbr, `[[`, "y")
    tsmat <- lapply(ans_mbr, `[[`, "t")
    ssmat <- lapply(ans_mbr, `[[`, "s")
    fsmat <- lapply(ans_mbr, `[[`, "f")
    nsmat <- lapply(ans_mbr, `[[`, "n")
    ysmat <- do.call(rbind, ysmat)
    tsmat <- do.call(rbind, tsmat)
    ssmat <- do.call(rbind, ssmat)
    fsmat <- do.call(rbind, fsmat)
    nsmat <- do.call(rbind, nsmat)
    rownames(ysmat) <- NULL
    rownames(tsmat) <- NULL
    rownames(ssmat) <- NULL
    rownames(fsmat) <- NULL
    rownames(nsmat) <- NULL
    rm(list = "ans_mbr")
  }

  ## (5) LFQ
  message("  Executing LFQ.")
  
  out <- pepLFQ(
    basenames = basenames, are_refs = are_refs, are_smpls = are_smpls, 
    xsmat = xsmat, ysmat = ysmat, tsmat = tsmat, ssmat = ssmat, cmat = cmat, 
    mbr_peps = mbr_peps, mbr_sps = mbr_sps, sp_centers = sp_centers, 
    rt_tol = rt_tol, mbr_ret_tol = mbr_ret_tol, step = 2E-3)
  yout <- do.call(cbind, out[["y"]])
  tout <- do.call(cbind, out[["t"]])
  sout <- do.call(cbind, out[["s"]])
  # rm(list = "out")
  
  if (ncol(yout) != ncol(tsmat)) {
    stop("Developer: check for data drops.")
  }
  
  ## (6) optional: handle missing reference values
  if (imp_refs && any(are_refs) && !identical(are_refs, are_smpls)) {
    message("  Handling missing reference values.")
    
    yvs <- colMeans(yout[are_refs, , drop = FALSE], na.rm = TRUE)
    
    if (length(i_narefs <- which(is.nan(yvs)))) {
      ynas <- colMeans(yout[are_smpls, i_narefs, drop = FALSE], na.rm = TRUE)
      vals <- sp_centers[mbr_sps[i_narefs]]
      vals[is.na(names(vals))] <- sp_centers[names(sp_centers) == ""]
      ynas <- ynas * (2^-vals)
      yout[which(are_refs)[[1]], i_narefs] <- ynas
    }
  }
  
  ## (7) outputs
  out <- vector("list", nrow(yout))
  
  for (i in seq_along(out)) {
    yda <- yout[i, ]
    tda <- tout[i, ]
    sda <- sout[i, ]
    # rows <- !is.na(yda)

    out[[i]] <- outi <- data.frame(
      pep_seq_modz  = mbr_peps, 
      pep_apex_ret  = tda, 
      pep_apex_scan = sda, 
      pep_exp_mz    = mbr_mzs, 
      pep_tot_int   = yda)
    
    readr::write_tsv(
      outi, file.path(dat_dir, "Peptide/cache", paste0("MBR_", basenames[[i]])))
  }
  
  ###
  # need to exclude: the same peak to different peptides???
  ###

  message("Completed LFQ.")
  
  attr(out, "sp_centers") <- sp_centers
  
  out
}


#' Add MBR peptide intensities
#' 
#' @param file A file name with prepending path to a peptide table.
#' @param pep_tbl A data frame corresponding to \code{file}.
#' @param dat_dir A working directory.
#' @inheritParams normPSM
addMBRpeps <- function (file, pep_tbl = NULL, dat_dir, group_psm_by = "pep_seq_mod")
{
  base_name <- basename(file)
  mbr_file  <- file.path(dat_dir, "Peptide/cache", paste0("MBR_", base_name))
  
  cols <- c("pep_seq_modz", "pep_seq_mod", "pep_seq", "pep_exp_mz", "pep_exp_z", 
            "pep_apex_ret", "pep_apex_scan", "pep_tot_int", "pep_score")
  
  if (is.null(pep_tbl)) {
    pep_tbl <- file |>
      readr::read_tsv(col_types = get_col_types(), show_col_types = FALSE) |> 
      suppressWarnings()
  }
  
  # is.na(pep_apex_ret): no MBR
  # is.na(pep_score): both MBR and no MBR
  
  df <- readr::read_tsv(mbr_file) |> 
    dplyr::left_join(pep_tbl[, c("pep_seq_modz", "pep_score")], 
                     by = "pep_seq_modz") |>
    tidyr::separate(pep_seq_modz, into = c("pep_seq_mod", "pep_exp_z"), 
                    sep = "@", remove = FALSE) |>
    dplyr::mutate(pep_seq = gsub("^[\\^_~]", "", pep_seq_mod), 
                  pep_seq = gsub("[\\^_~]$", "", pep_seq), 
                  pep_seq = toupper(pep_seq), 
                  pep_exp_z = as.integer(pep_exp_z))
  
  # for SD calculations
  rows <- is.na(df$pep_apex_ret)
  df   <- df[!rows, ]
  
  if (!all(cols %in% names(df))) {
    stop("Developer: missing columns.")
  }
  
  # (2) add `pep_ret_sd`
  df <- df |> 
    dplyr::group_by_at("pep_seq_mod") |> 
    dplyr::mutate(N = dplyr::n())
  
  rows <- df$N == 1L
  df0  <- df[rows, ]
  df   <- df[!rows, ]
  df0$N <- df$N <- NULL
  df0$pep_ret_sd <- 0
  
  sds <- calc_pep_retsd(df, group_psm_by = "pep_seq_mod")
  df  <- sds |>
    dplyr::left_join(df, by = "pep_seq_mod")
  df <- df[, names(df0)]
  
  # ASSEGTIPQVQR; 636.8285; partial peak caused by 6 ppm in ms1 tracing; need 10 ppm
  # re-exam large differences in MBR, increase tracing tolerance to 10 ppm
  
  if (FALSE) {
    out_cols <- c("pep_seq_mod", "pep_seq", "pep_tot_int", "pep_apex_ret", "pep_apex_scan")
    list(pep_tbl = pep_tbl, mbr0 = df0[, out_cols], mbr1 = df)
  }
  
  # (3) delay: pep_seq_modz -> pep_seq_mod
  if (FALSE) {
    mbr <- groupMZPepZ(df, group_psm_by = group_psm_by)
  }
  else {
    mbr <- df[, c("pep_seq_modz", "pep_seq_mod", "pep_seq", "pep_tot_int", 
                  "pep_apex_ret", "pep_apex_scan")]
  }
  
  # some pep_seq_mod in drt0[, names(mbr)] can be in mbr
  mbr <- dplyr::bind_rows(df0[, names(mbr)], mbr)
  out_name <- file.path(dat_dir, "Peptide/cache", paste0("MBRpeps_", base_name))
  readr::write_tsv(mbr, out_name)
  
  # (4) side effect to return the Peptide.txt
  pep_tbl
}


#' Fill NA retention times with the MS2 triggering retention times
#' 
#' To obtain a reference point for reassess apexes and MBR.
#' 
#' @param mbr_rets A vector of retention times.
#' @param mbr_peps A vector of \code{pep_seq_modz}.
#' @param dat_dir The working directory.
fillMBRnaRTs <- function (mbr_rets, mbr_peps, dat_dir) 
{
  if (!length(nas <- which(is.nan(mbr_rets)))) {
    return(mbr_rets)
  }
  
  peps_na  <- mbr_peps[nas]
  filelist <- list.files(
    path = file.path(dat_dir, "PSM"), pattern = "_PSM_N\\.txt$", full.names = TRUE)
  
  if (!(n_files <- length(filelist))) {
    stop("Files of \"_PSM_N.txt\" not found.")
  }
  
  df <- lapply(filelist, function (x) {
    df <- readr::read_tsv(
      x, col_types = get_col_types(), show_col_types = FALSE) |> 
      suppressWarnings()
    
    if (!"pep_seq_modz" %in% names(df)) {
      df[["pep_seq_modz"]] <- paste0(df[["pep_seq_mod"]], "@", df[["pep_exp_z"]])
    }
    
    df <- df[, c("pep_seq_modz", "pep_ret_range")]
    df <- df[with(df, df[["pep_seq_modz"]] %in% peps_na), ]
  }) |>
    dplyr::bind_rows() |>
    dplyr::group_by(pep_seq_modz) |>
    dplyr::summarise(rt = median(pep_ret_range))
  
  mbr_rets[nas] <- df[["rt"]][match(peps_na, df[["pep_seq_modz"]])]
  
  mbr_rets
}


#' LFQ of peptides
#'
#' @param basenames The base names of a peptide table files
#'   (\code{TMTset[i]_LCMSinj[j]_Peptide_N.txt}).
#' @param are_refs Logical; are references or not.
#' @param are_smpls Logical; are samples or not.
#' @param xsmat A matrix of m-over-z values.
#' @param ysmat A matrix of apex intensities.
#' @param tsmat A matrix of apex retention times.
#' @param ssmat A matrix of apex scan numbers.
#' @param cmat A matrix of spectrum counts.
#' @param mbr_peps A vector of \code{pep_seq_modz} sequences in the universe.
#' @param mbr_sps A vector of peptide species in the universe.
#' @param sp_centers The centers of log2FC for each species; names: species,
#'   values: log2FC.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param step The retention-time error in \code{ppm / 1e6}.
#' @param err_log2r Not yet used. Error tolerance in log2FC.
#' @param rt_tol Error tolerance in retention times.
pepLFQ <- function (basenames, are_refs, are_smpls, xsmat, ysmat, tsmat, ssmat, 
                    cmat, mbr_peps, mbr_sps, sp_centers, rt_tol = 25, # was 15
                    mbr_ret_tol = 25, step = 2E-3, err_log2r = .25)
  
{
  # check again if multiple fractions under a set of TMTset1_LCMSinj1_Peptide_N
  # need to compile retention times, moverzs and intensities across fractions...
  
  # if (is_ref) { sp_centers <- -sp_centers }
  cent0 <- sp_centers[names(sp_centers) == ""]
  cents <- sp_centers[mbr_sps]
  cents[is.na(cents)] <- cent0
  
  yout <- tout <- sout <- 
    rep_len(list(vector("numeric", n_row <- nrow(tsmat))), n_col <- ncol(tsmat))
  
  null_dbl <- rep_len(NA_real_, n_row)
  null_int <- rep_len(NA_integer_, n_row)

  for (i in 1:n_col) {
    ## (1) collapse bins of T, Y and S values
    ans <- collapseSTY(
      xs = tsmat[, i], ys = ysmat[, i], zs = ssmat[, i], lwr = 10, step = step)
    tmat <- ans[["x"]]
    ymat <- ans[["y"]]
    smat <- ans[["z"]]
    unv  <- ans[["u"]]
    spc  <- cents[[i]]
    
    nc <- ncol(tmat)
    
    if (!nc) {
      yout[[i]] <- null_dbl
      tout[[i]] <- null_dbl
      sout[[i]] <- null_int
      
      next
    }
    
    if (nc == 1L) {
      yout[[i]] <- ymat[, 1]
      tout[[i]] <- tmat[, 1]
      sout[[i]] <- smat[, 1]
      
      next
    }
    
    ## (2) remove retention outliers
    for (j in 1:nc) {
      tvals <- tmat[, j]
      tbads <- abs(tvals - median(tvals, na.rm = TRUE)) > rt_tol
      ibads <- .Internal(which(tbads))
      
      if (length(ibads)) {
        tmat[ibads, j] <- NA_real_
        ymat[ibads, j] <- NA_real_
        smat[ibads, j] <- NA_real_
      }
    }

    ## (3) calculate log2FCs
    ys_smpls <- ymat[are_smpls, , drop = FALSE]
    ys_refs  <- ymat[are_refs,  , drop = FALSE]
    ymeans   <- colMeans(ys_refs, na.rm = TRUE)
    
    if (all(is.nan(ymeans))) { # all are without reference values
      it1 <- .Internal(which.min(calc_mat_sds(tmat[are_smpls, , drop = FALSE])))
      
      if (length(it1)) {
        yout[[i]] <- ymat[, it1]
        tout[[i]] <- tmat[, it1]
        sout[[i]] <- smat[, it1]
        
        next
      }
      else { # e.g. all with single sample values
        # or use the one closest to mbr_rets..
        yout[[i]] <- ymat[, 1]
        tout[[i]] <- tmat[, 1]
        sout[[i]] <- smat[, 1]
        
        next
      }
    }
    
    rmat <- sweep(ys_smpls, 2, ymeans, "/")
    rmat <- log2(rmat)
    
    ## (4) find the best Y column
    ds <- abs(colMeans(rmat, na.rm = TRUE) - spc)
    p  <- .Internal(which.min(ds))
    
    if (length(p)) {
      yout[[i]] <- ymat[, p]
      tout[[i]] <- tmat[, p]
      sout[[i]] <- smat[, p]
    }
    else { # e.g., one column exclusive samples, the other exclusive references
      iy1 <- .Internal(which.min(calc_mat_sds(ymat)))
      if (length(iy1)) {
        yout[[i]] <- ymat[, iy1]
        tout[[i]] <- tmat[, iy1]
        sout[[i]] <- smat[, iy1]
      }
      else { # e.g. all are with references only, or all single values...
        # may be use the column closest in RT to MS2 scan
        iy2 <- .Internal(which.max(colSums(ymat, na.rm = TRUE)))
        if (length(iy2)) {
          yout[[i]] <- ymat[, iy2]
          tout[[i]] <- tmat[, iy2]
          sout[[i]] <- smat[, iy2]
        }
        else {
          yout[[i]] <- ymat[, 1]
          tout[[i]] <- tmat[, 1]
          sout[[i]] <- smat[, 1]
        }
      }
    }
  }
  
  out <- list(y = yout, t = tout, s = sout)
}


#' Make MS1 S, T and Y matrices across adjacent scans
#'
#' Near identical to \code{mzion:::collapse_mms1ints} with three fields instead
#' of two.
#'
#' @param xs Vectors of X (retention time) values in an ASCENDING order.
#' @param ys Vectors of Y (intensity) values corresponding to \code{xs}.
#' @param zs Vectors of Z (scan numbers) values corresponding to \code{xs}.
#' @param lwr The lower (retention time) limit.
#' @param step The bin size in converting numeric X values to integers.
#' @param coll Logical; collapse adjacent columns or not.
#' @param look_back Logical; look up the preceding MS bin or not.
#' @importFrom fastmatch %fin%
collapseSTY <- function (xs = NULL, ys = NULL, zs = NULL, lwr = 10, step = 2e-3, 
                         coll = TRUE, look_back = TRUE)
{
  ixs <- lapply(xs, mzion::index_mz, lwr, step)
  
  for (i in seq_along(xs)) {
    ix <- ixs[[i]]
    ps <- .Internal(which(duplicated(ix)))
    
    if (length(ps)) {
      ixs[[i]] <- ix[-ps]
      xs[[i]]  <- xs[[i]][-ps]
      ys[[i]]  <- ys[[i]][-ps]
      zs[[i]]  <- zs[[i]][-ps]
    }
  }
  
  ## (1) make matrices
  unv  <- .Internal(unlist(ixs, recursive = FALSE, use.names = FALSE))
  unv  <- sort(unique(unv))
  lenu <- length(unv)
  lenx <- length(xs)
  ups  <- lapply(ixs, function (x) unv %fin% x)
  
  # note one-to-one correspondence between ixs and xs
  xmat <- mzion::mapcoll_xyz(
    vals = xs, ups = ups, lenx = lenx, lenu = lenu, direct_out = TRUE)
  ymat <- mzion::mapcoll_xyz(
    vals = ys, ups = ups, lenx = lenx, lenu = lenu, direct_out = TRUE)
  zmat <- mzion::mapcoll_xyz(
    vals = zs, ups = ups, lenx = lenx, lenu = lenu, direct_out = TRUE)
  
  if (!coll) { return(list(x = xmat, y = ymat, z = zmat, u = unv)) }
  
  ## (2) collapses adjacent entries with +/-1 in bin indexes
  if (is.null(ps <- mzion::find_gates(unv))) { # all discrete values
    return(list(x = xmat, y = ymat, z = zmat, u = unv))
  }
  
  # for recording matrix columns to be dropped
  ps1 <- vector("integer", lenp <- length(ps))
  
  for (i in 1:lenp) {
    c12 <- ps[[i]]
    c1  <- c12[[1]]
    c2  <- c12[[2]]
    
    if (c2 == 0L) { # no following adjacent
      next
    }
    
    rows <- .Internal(which(!is.na(xmat[, c1])))
    xmat[rows, c2] <- xmat[rows, c1]
    ymat[rows, c2] <- ymat[rows, c1]
    zmat[rows, c2] <- zmat[rows, c1]
    # no need since columns traverse from left to right and c1 will be dropped
    xmat[rows, c1] <- NA_real_
    ymat[rows, c1] <- NA_real_
    zmat[rows, c1] <- NA_real_
    
    if (look_back && i < lenp && (af <- ps[[i+1]][[1]]) == (c2 + 1L)) {
      # with NA rows in c2
      if (length(rows2 <- .Internal(which(is.na(xmat[, c2]))))) {
        xmat[rows2, c2] <- xmat[rows2, af]
        ymat[rows2, c2] <- ymat[rows2, af]
        zmat[rows2, c2] <- zmat[rows2, af]
        # about ok to enable, but may be unnecessary
        xmat[rows2, af] <- NA_real_
        ymat[rows2, af] <- NA_real_
        zmat[rows2, af] <- NA_real_
      }
    }
    
    ps1[[i]] <- c1
  }
  
  # note `-ps1` ok in that at least one ps1 is not 0
  # identical(xmat[, -c(0, 2:3), drop = FALSE], xmat[, -c(2:3), drop = FALSE])
  # !identical(xmat[, 0, drop = FALSE], xmat)
  xmat <- xmat[, -ps1, drop = FALSE]
  ymat <- ymat[, -ps1, drop = FALSE]
  zmat <- zmat[, -ps1, drop = FALSE]
  unv  <- unv[-ps1]
  
  # possible additional all-NA columns by look_back
  if (look_back && 
      length(nas <- .Internal(which(colSums(is.na(xmat)) == lenx)))) {
    xmat <- xmat[, -nas, drop = FALSE]
    ymat <- ymat[, -nas, drop = FALSE]
    zmat <- zmat[, -nas, drop = FALSE]
    unv  <- unv[-nas]
  }
  
  list(x = xmat, y = ymat, z = zmat, u = unv)
}


#' Calculate SD for matrix data
#' 
#' @param mat A data matrix.
calc_mat_sds <- function (mat)
{
  if (!(n_col <- ncol(mat))) {
    return(NA_real_)
  }
  
  sds <- vector("numeric", )
  for (i in 1:ncol(mat)) {
    sds[[i]] <- sd(mat[, i], na.rm = TRUE)
  }
  
  sds
}


#' MBR of peptides
#'
#' @param base_name The base name of a peptide table file
#'   (\code{TMTset[i]_LCMSinj[j]_Peptide_N.txt}).
#' @param ms1files A list of file names of retention-time calibrated RAW MS1
#'   data corresponding to the peptide table at \code{base_name} (More than one
#'   file if with analyte pre-fractionation for a sample).
#' @param row_id An integer corresponding to \code{base_name}, as well as the
#'   row number in \code{xmat} etc.
#' @param xmat A matrix of X (m-over-z) values under the peptide table at the
#'   \code{base_name}. The length is the number of peptides in the universe. NA
#'   values correspond to peptides (\code{pep_seq_modz}) that are not present in
#'   the current peptide table but found in others. Each cell in the matrix
#'   contains a list of all possible m-over-z candidates.
#' @param ymat A matrix of apex intensities.
#' @param tmat A matrix of apex retention times.
#' @param smat A matrix of apex scan numbers.
#' @param fmat A matrix of FWHMs.
#' @param nmat A matrix of apex width in counts of MS2 scans.
#' @param mbr_mzs A vector of m-over-z values in the universe (sample average).
#' @param mbr_rets A vector of retention times in the universe (sample average).
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param step The mass error in \code{ppm / 1e6}.
#' @param gap A gap in retention-time window in seconds for finding gates in the
#'   Y values of LC. The gap value serves as a margin to trace out a whole peak
#'   with its apex within the \code{mbr_ret_tol}.
#' @param min_y The cut-off of intensity values in MBR.
#' @param yco The cut-off in Y-values.
#' @param fwhm_co The cut-off in FWHM values.
#' @param dat_dir The working data directory.
saddMBR <- function (base_name, ms1files, row_id = 1L, 
                     xmat = NULL, ymat = NULL, tmat = NULL, 
                     smat = NULL, fmat = NULL, nmat = NULL, 
                     mbr_mzs, mbr_rets, dat_dir, mbr_ret_tol = 25.0, 
                     gap = mbr_ret_tol + 90.0, min_y = 2e6, yco = 100, 
                     step = 1e-5, fwhm_co = .5)
{
  # check again if multiple fractions under a set of TMTset1_LCMSinj1_Peptide_N
  # need to compile retention times, moverzs and intensities across fractions...
  
  ms1full <- dplyr::bind_rows(lapply(ms1files, qs::qread))
  tss <- unname(ms1full$ret_time)
  xss <- ms1full$msx_moverzs
  yss <- ms1full$msx_ints
  sss <- ms1full$scan_num
  # oss <- ms1full$orig_scan
  rm(list = c("ms1full"))
  
  xmi <- xmat[row_id, ] # no changes
  ymi <- ymat[row_id, ]
  tmi <- tmat[row_id, ]
  smi <- smat[row_id, ]
  fmi <- fmat[row_id, ]
  nmi <- nmat[row_id, ]

  cols  <- which(nas <- is.na(tmat[row_id, ]))
  for (i in seq_along(cols)) {
    col     <- cols[[i]]
    mbr_mz  <- mbr_mzs[[col]]
    mbr_ret <- mbr_rets[[col]]
    rng     <- .Internal(which(tss >= (mbr_ret - gap) & tss <= (mbr_ret + gap)))
    xys     <- 
      extract_mbry(xs = xss[rng], ys = yss[rng], mbr_mz = mbr_mz, step = step)
    if (length(xys) == 1L && is.na(xys)) { next }
    yhats   <- xys[["y"]]
    if (all(is.na(yhats))) { next }
    xhats   <- xys[["x"]]
    tsx     <- tss[rng]
    
    gates <- mzion:::find_lc_gates(
      ys = yhats, ts = tsx, yco = yco, n_dia_scans = 6L, ytot_co = min_y)
    if (is.null(gates)) { next }
    
    apexs  <- gates[["apex"]]
    ymi[[col]] <- gates[["yints"]]
    tmi[[col]] <- tsx[apexs]
    smi[[col]] <- sss[rng][apexs]
    fmi[[col]] <- gates[["fwhm"]]
    nmi[[col]] <- gates[["ns"]]
  }
  
  list(x = xmi, y = ymi, t = tmi, s = smi, f = fmi, n = nmi)
}


