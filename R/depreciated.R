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
