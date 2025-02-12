#' Helper of adding PSM apexes
#'
#' By each PSM file of \code{TMTSet1_LCMSinj1.txt}, \code{TMTSet2_LCMSinj1.txt}
#' etc.
#'
#' @param file A PSM file name, e.g., \code{TMTSet1_LCMSinj1_PSM_N.txt}.
#' @param dfs A list split of PSM data under the \code{file} at different
#'   \code{RAW_File}.
#' @param trace_files The file names of \code{ms1apexes_} corresponding to
#'   \code{dfs}.
#' @param fits The models of retention-time fittings corresponding to
#'   \code{dfs}.
#' @param dat_dir A working directory.
#' @param path_ms1 The file path to MS1 data from Mzion.
#' @param lfq_ret_tol Retention time tolerance for LFQ.
#' @param ok_mbr Logical; OK to perform MBR or not.
hadd_psm_apexes <- function (file = NULL, dfs = NULL, trace_files = NULL, 
                             fits = NULL, dat_dir = NULL, path_ms1 = NULL, 
                             lfq_ret_tol = 60L, ok_mbr = TRUE)
{
  out <- mapply(add_psm_apexes, df = dfs, trace_file = trace_files, 
                MoreArgs = list(dat_dir = dat_dir, path_ms1 = path_ms1), 
                SIMPLIFY = FALSE, USE.NAMES = TRUE)
}


#' Add PSM apexes by each \code{RAW_File}.
#'
#' @param df A fraction of PSM table, e.g. \code{TMTSet1_LCMSinj1_PSM_N.txt},
#'   under a single \code{RAW_File}.
#' @param trace_file A file name of MS1 traces (\code{ms1apexes_[...].rds}).
#' @param dat_dir A working directory.
#' @param path_ms1 The path to the data of MS1 traces.
#' @param cols Column keys
add_psm_apexes <- function (df = NULL, trace_file = NULL, dat_dir = NULL, 
                            path_ms1 = NULL, 
                            # no need of pep_exp_z since the same pep_seq_mod at 
                            # different pep_exp_mz must have different pep_exp_z?
                            # later to simply the keys in `cols`
                            cols = c("raw_file", "pep_tot_int", "pep_exp_mz", "pep_seq_mod", 
                                     "pep_ret_range", "pep_apex_ret", "pep_scan_num", 
                                     "pep_orig_scan", "pep_apex_scan", "pep_n_apexes")) 
{
  ## MS1 traces
  #  level-1: (chimeric) precursors
  #    level-2: apexes for each precursor
  trs <- qs::qread(file.path(path_ms1, trace_file))
  mts <- fastmatch::fmatch(df$pep_orig_scan, trs$orig_scan)
  # mts <- mts[!is.na(mts)] # shouldn't be any
  
  ans  <- df[, cols] |>
    dplyr::bind_cols(trs[mts, !names(trs) %in% c("orig_scan", "ms_level")])
  lens <- lengths(ans[["apex_ts"]])
  rowx <- lens > 1L # with chimeric precursors
  rows <- !rowx     # with single precursor
  ans1 <- ans[rows, ]
  ans1[["apex_ts"]] <- lapply(ans1[["apex_ts"]],  `[[`, 1L) # unlist
  ans1[["apex_xs"]] <- lapply(ans1[["apex_xs"]],  `[[`, 1L)
  ans1[["apex_ys"]] <- lapply(ans1[["apex_ys"]],  `[[`, 1L)
  ans1[["apex_ps"]] <- lapply(ans1[["apex_ps"]],  `[[`, 1L)
  
  ansn <- cleanRTChim(df = ans[rowx, ]) # by the "correctness" of `pep_exp_mz`
  
  if (FALSE) {
    lens <- lengths(ansn$apex_ps)
    rowx <- lens > 3L
    ansx <- ansn[rowx, ]
  }
  
  ans[rows, ] <- ans1
  ans[rowx, ] <- ansn
  
  nms  <- names(ans)
  colx <- nms[!nms %in% cols]
  df[, cols] <- ans[, cols] # identical?
  
  df <- dplyr::bind_cols(df, ans[, colx])
}


#' Find the MBR intensity
#' 
#' Not used.
#' 
#' @param ys A vector of intensity values.
#' @param ts A vector of time values.
#' @param ss A vector of scan numbers.
#' @param mbr_ret The reference MBR retention time.
#' @param mbr_y The reference MBR intensity.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param sp_cent A species center (at log2 scale).
#' @param max_mbr_fold The maximum absolute fold change in MBR.
#' @param n_dia_scans The maximum number of zero-intensity scans for gap filling
#'   in peak tracing. Not adjustable by users but synchronized with mzion.
#' @param min_y The cut-off of intensity values in MBR. Change to a smaller
#'   value with PASEF.
#' @param yco The cut-off of intensity values.
#' @param fwhm_co The cut-off in FWHM values.
#' @param zero A zero value for replacements.
find_mbr_int <- function (ys, ts, ss, mbr_ret, mbr_y, mbr_ret_tol = 25, 
                          sp_cent = 0.0, 
                          max_mbr_fold = 20L, n_dia_scans = 6L, min_y = 2e6, 
                          yco = 100, fwhm_co = .5, zero = NA_real_)
{
  nout  <- list(y = zero, t = zero, s = zero)
  
  if (is.na(mbr_y)) {
    return(nout)
  }
  
  gates <- mzion::find_lc_gates(
    ys = ys, ts = ts, yco = yco, n_dia_scans = n_dia_scans)
  
  if (is.null(gates)) {
    return(nout)
  }
  
  fwhms  <- gates[["fwhm"]]
  apexs  <- gates[["apex"]]
  ranges <- gates[["ranges"]]
  yints  <- gates[["yints"]]
  ns     <- gates[["ns"]]
  xstas  <- gates[["xstas"]]
  lenp   <- length(apexs)
  
  oksfw <- fwhms > fwhm_co
  if (!(noksfw <- length(oksfw))) {
    return(nout)
  }
  
  if (noksfw < lenp) {
    apexs  <- apexs[oksfw]
    ranges <- ranges[oksfw]
    yints  <- yints[oksfw]
    ns     <- ns[oksfw]
    xstas  <- xstas[oksfw]
    lenp   <- length(oksfw)
  }
  
  if (lenp == 0L) {
    return(nout)
  }
  
  # (3) remove one-hit-wonders and spikes
  oks1 <- .Internal(which(ns > 20L))
  oks2 <- .Internal(which(ns > 15L))
  
  if (noks1 <- length(oks1)) {
    if (noks1 < lenp) {
      apexs  <- apexs[oks1]
      ranges <- ranges[oks1]
      yints  <- yints[oks1]
      ns     <- ns[oks1]
      xstas  <- xstas[oks1]
      lenp   <- length(apexs)
    }
  }
  else if (noks2 <- length(oks2)) {
    if (noks2 < lenp) {
      apexs  <- apexs[oks2]
      ranges <- ranges[oks2]
      yints  <- yints[oks2]
      ns     <- ns[oks2]
      xstas  <- xstas[oks2]
      lenp   <- length(apexs)
    }
  }
  
  upr <- mbr_y * max_mbr_fold
  lwr <- mbr_y / max_mbr_fold
  
  if (lenp == 1L) {
    
    # if no log2rs comparable to sp_cent -> expand to local pattern look-ups
    
    tval <- ts[[apexs]]
    
    if (abs(tval - mbr_ret) <= mbr_ret_tol) {
      # may be removed...
      if (yints <= upr && yints >= lwr) {
        return(list(y = yints, t = tval, s = ss[[apexs]]))
      }
      else {
        # return(list(y = mbr_y, t = tval, s = ss[[apexs]]))
        return(nout)
      }
    }
    else {
      return(nout)
    }
  }
  
  tvals <- ts[apexs]
  scans <- ss[apexs]
  tdiff <- abs(tvals - mbr_ret)
  ydiff <- abs(log2(yints / mbr_y) - sp_cent)
  idxt  <- .Internal(which.min(tdiff))
  idxy  <- .Internal(which.min(ydiff))
  
  if (tdiff[idxy] <= mbr_ret_tol) {
    yi <- yints[idxy]
    
    if (yi <= upr && yi >= lwr) {
      return(list(y = yi, t = tvals[idxy], s = scans[idxy]))
    }
    else {
      return(nout)
    }
  }
  
  if (tdiff[idxt] <= mbr_ret_tol) {
    yi <- yints[idxt]
    
    if (yi <= upr && yi >= lwr) {
      return(list(y = yi, t = tvals[idxt], s = scans[idxt]))
    }
    else {
      return(nout)
    }
  }
  
  return(nout)
}


