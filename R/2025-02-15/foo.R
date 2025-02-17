#' LFQ of peptides
#' 
#' Peptides are unique to species.
#'
#' @param basenames The base names of a peptide table files
#'   (\code{TMTset[i]_LCMSinj[j]_Peptide_N.txt}).
#' @param are_refs Logical; are references or not.
#' @param are_smpls Logical; are samples or not.
#' @param xsmat A matrix of m-over-z values.
#' @param ysmat A matrix of apex intensities.
#' @param tsmat A matrix of apex retention times.
#' @param ssmat A matrix of apex scan numbers.
#' @param ymat A matrix of (first-pass) apex Y values. The binning of local
#'   patterns by retention time intervals can have boundary effects by
#'   separating close retention times into bin indexes \eqn{ge +1} or \eqn{le
#'   -1}. If the log2FC is by ymat is closer to a species center than newly
#'   derived log2FC from local patterns, use the Y values from the first pass.
#' @param tmat A matrix of (first-pass) apex retention times.
#' @param smat A matrix of (first-pass) apex scan numbers.
#' @param cmat A matrix of spectrum counts.
#' @param mbr_peps A vector of \code{pep_seq_modz} sequences in the universe.
#' @param mbr_sps A vector of peptide species in the universe.
#' @param sp_centers The centers of log2FC for each species; names: species,
#'   values: log2FC.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param step The step size in binning retention times.
#' @param err_log2r Not yet used. Error tolerance in log2FC.
#' @param new_na_species A replace value for NA species.
#' @param rt_tol Error tolerance in retention times.
#' @param reference_spec_counts Logical; reference protein spectrum counts or
#'   not.
pepLFQ_orig <- function (basenames, are_refs, are_smpls, xsmat, ysmat, tsmat, ssmat, 
                    ymat, tmat, smat, cmat, mbr_peps, mbr_sps, sp_centers, 
                    new_na_species = ".other", rt_tol = 25, mbr_ret_tol = 25, 
                    reference_spec_counts = FALSE, step = 5E-3, err_log2r = .25)
  
{
  # check again if multiple fractions under a set of TMTset1_LCMSinj1_Peptide_N
  # need to compile retention times, moverzs and intensities across fractions...
  
  n_row <- nrow(tsmat)
  n_col <- ncol(tsmat)
  yout  <- tout <- sout <- rep_len(list(vector("numeric", n_row)), n_col)
  null_dbl <- rep_len(NA_real_, n_row)
  null_int <- rep_len(NA_integer_, n_row)
  
  # distances of reference log2FCs from the "first pass" to cents
  cents <- sp_centers[mbr_sps]
  dys1  <- calc_mat_log2s_to_refs(ymat, are_smpls, are_refs) - cents
  
  if (reference_spec_counts) {
    cs1 <- calc_mat_log2s_to_refs(cmat, are_smpls, are_refs)
  }
  else {
    cs1 <- rep_len(Inf, n_col)
  }
  rownames(cmat) <- colnames(cmat) <- NULL
  
  for (i in 1:n_col) {
    if (is.na(mbr_sps[[i]])) {
      yout[[i]] <- null_dbl
      tout[[i]] <- null_dbl
      sout[[i]] <- null_int
      
      next
    }
    
    ## (1) collapse bins of T, Y and S values
    # i <- which(mbr_peps == "LLDVVHPAAK@3")
    ans <- collapseSTY(
      xs = tsmat[, i], ys = ysmat[, i], zs = ssmat[, i], lwr = 10, step = step)
    anst <- ans[["x"]]
    ansy <- ans[["y"]]
    anss <- ans[["z"]]
    unv  <- ans[["u"]]
    
    if (!(nc <- ncol(anst))) {
      yout[[i]] <- null_dbl
      tout[[i]] <- null_dbl
      sout[[i]] <- null_int
      
      next
    }
    
    # border effect: may disable since values can traverse for every +2
    # may check +2 +2 +2 senario...
    if (nc > 1L && length(us <- which(diff(unv) == 2L))) {
      for (u in us) {
        rows <- which(is.na(ansy[, u2 <- u + 1L]) & !is.na(ansy[, u]))
        
        if (length(rows)) {
          ansy[rows, u2] <- ansy[rows, u]
          anst[rows, u2] <- anst[rows, u]
          anss[rows, u2] <- anss[rows, u]
          # ansy[rows, u] <- NA_real_
          # anst[rows, u] <- NA_real_
          # anss[rows, u] <- NA_real_
        }
      }
    }
    
    spc   <- cents[[i]]
    
    # for cross-referencing
    dy1i  <- dys1[[i]]
    dc1i  <- cs1[[i]] - spc
    ady1i <- abs(dy1i)
    adc1i <- abs(dc1i)
    
    if (nc == 1L) {
      p      <- 1L
      ysi    <- ansy[, p]
      ysmpls <- ysi[are_smpls]
      yrefs  <- ysi[are_refs]
      ybar0  <- mean(yrefs, na.rm = TRUE)
      di     <- mean(log2(ysmpls / ybar0) - spc, na.rm = TRUE)
      adi    <- abs(di)
      bst    <- find_best_lfqvec(d1 = di, d2 = dy1i, dc = dc1i)
      
      if (isTRUE(bst == 1L)) {
        yout[[i]] <- ysi
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
      }
      else if (isTRUE(bst == 2L)) {
        yout[[i]] <- ymat[, i]
        tout[[i]] <- tmat[, i]
        sout[[i]] <- smat[, i]
      }
      else if (isTRUE(bst == 3L) && !(is.nan(di) || is.na(di))) {
        ybar1 <- mean(ysmpls / 2^(dc1i + spc), na.rm = TRUE)
        ysi[are_refs] <- yrefs * ybar1 / ybar0
        yout[[i]] <- ysi
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
      }
      else {
        yout[[i]] <- ysi
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
      }
      
      next
    }
    
    ## (2) remove retention outliers
    for (j in 1:nc) {
      tvals <- anst[, j]
      # reduce the value to be consistent with step...
      tbads <- abs(tvals - median(tvals, na.rm = TRUE)) > rt_tol
      ibads <- .Internal(which(tbads))
      
      if (length(ibads)) {
        anst[ibads, j] <- NA_real_
        ansy[ibads, j] <- NA_real_
        anss[ibads, j] <- NA_real_
      }
    }
    
    ## (3) calculate log2FCs
    ys_smpls <- ansy[are_smpls, , drop = FALSE]
    ys_refs  <- ansy[are_refs,  , drop = FALSE]
    ymeans   <- colMeans(ys_refs, na.rm = TRUE)
    
    # all are without reference values
    if (all(is.nan(ymeans))) {
      p <- .Internal(which.min(calc_mat_sds(anst[are_smpls, , drop = FALSE])))
      
      if (length(p)) {
        ysi <- ansy[, p]
        yout[[i]] <- ysi
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
        
        if (FALSE) {
          ysi <- ansy[, p]
          # all NA?
          adi <- abs(calc_vec_log2s_to_refs(ysi, are_smpls, are_refs) - spc)
          bst <- .Internal(which.min(c(adi, ady1i)))
          
          if (isTRUE(bst == 1L)) {
            yout[[i]] <- ysi
            tout[[i]] <- anst[, p]
            sout[[i]] <- anss[, p]
          }
          else if (isTRUE(bst == 2L)) {
            yout[[i]] <- ymat[, i]
            tout[[i]] <- tmat[, i]
            sout[[i]] <- smat[, i]
          }
          else {
            yout[[i]] <- ysi
            tout[[i]] <- anst[, p]
            sout[[i]] <- anss[, p]
          }
        }
      }
      else {
        # e.g. all with single sample values (may be due to the boundary effect)
        if (is.na(ady1i)) {
          yout[[i]] <- ansy[, 1]
          tout[[i]] <- anst[, 1]
          sout[[i]] <- anss[, 1]
        }
        else {
          yout[[i]] <- ymat[, i]
          tout[[i]] <- tmat[, i]
          sout[[i]] <- smat[, i]
        }
      }
      
      next
    }
    
    rmat <- sweep(ys_smpls, 2, ymeans, "/")
    rmat <- log2(rmat)
    
    ## (4) find the best Y column
    ds  <- colMeans(rmat, na.rm = TRUE) - spc
    ads <- abs(ds)
    p   <- .Internal(which.min(ads))
    
    if (length(p)) {
      # if any NA in ansy[, p] && and log2FC outside the contour...
      di <- ds[[p]]
      bst <- find_best_lfqvec(d1 = di, d2 = dy1i, dc = dc1i)
      
      if (isTRUE(bst == 1L)) {
        yout[[i]] <- ansy[, p]
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
      }
      else if (isTRUE(bst == 2L)) {
        yout[[i]] <- ymat[, i]
        tout[[i]] <- tmat[, i]
        sout[[i]] <- smat[, i]
      }
      else if (isTRUE(bst == 3L) && !(is.nan(di) || is.na(di))) {
        ysi    <- ansy[, p]
        ysmpls <- ysi[are_smpls]
        yrefs  <- ysi[are_refs]
        ybar0  <- mean(yrefs, na.rm = TRUE)
        ybar1  <- ysmpls / 2^(dc1i + spc)
        ysi[are_refs] <- yrefs * ybar1 / ybar0
        yout[[i]] <- ysi
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
      }
      else {
        yout[[i]] <- ansy[, p]
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
      }
      
      next
    }
    
    # e.g., one column exclusive samples, the other exclusive references
    p <- .Internal(which.min(calc_mat_sds(ansy)))
    
    if (length(p)) {
      ysi <- ansy[, p]
      adi  <- abs(calc_vec_log2s_to_refs(ysi, are_smpls, are_refs) - spc)
      bst <- .Internal(which.min(c(adi, ady1i)))
      
      if (isTRUE(bst == 1L)) {
        yout[[i]] <- ysi
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
      }
      else {
        yout[[i]] <- ymat[, i]
        tout[[i]] <- tmat[, i]
        sout[[i]] <- smat[, i]
      }
      
      next
    }
    
    # e.g. all are with references only, or all single values...
    # may be use the column closest in RT to MS2 scan
    p <- .Internal(which.max(colSums(ansy, na.rm = TRUE)))
    
    if (length(p)) {
      ysi <- ansy[, p]
      yout[[i]] <- ysi
      tout[[i]] <- anst[, p]
      sout[[i]] <- anss[, p]
      
      next
      
      if (FALSE) {
        ysi <- ansy[, p]
        # all NA? May be just check ymat[, i] see if is due to border effect 
        adi <- abs(calc_vec_log2s_to_refs(ysi, are_smpls, are_refs) - spc)
        bst <- .Internal(which.min(c(adi, ady1i)))
        
        if (isTRUE(bst == 1L)) {
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        else if (isTRUE(bst == 2L)) {
          yout[[i]] <- ymat[, i]
          tout[[i]] <- tmat[, i]
          sout[[i]] <- smat[, i]
        }
        else {
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        
        next
      }
    }
    
    yout[[i]] <- ansy[, 1]
    tout[[i]] <- anst[, 1]
    sout[[i]] <- anss[, 1]
  }
  
  out <- list(y = yout, t = tout, s = sout)
}


#' LFQ of peptides
#' 
#' Peptides are unique to species.
#'
#' @param basenames The base names of a peptide table files
#'   (\code{TMTset[i]_LCMSinj[j]_Peptide_N.txt}).
#' @param are_refs Logical; are references or not.
#' @param are_smpls Logical; are samples or not.
#' @param xsmat A matrix of m-over-z values.
#' @param ysmat A matrix of apex intensities.
#' @param tsmat A matrix of apex retention times.
#' @param ssmat A matrix of apex scan numbers.
#' @param ymat A matrix of (first-pass) apex Y values. The binning of local
#'   patterns by retention time intervals can have boundary effects by
#'   separating close retention times into bin indexes \eqn{ge +1} or \eqn{le
#'   -1}. If the log2FC is by ymat is closer to a species center than newly
#'   derived log2FC from local patterns, use the Y values from the first pass.
#' @param tmat A matrix of (first-pass) apex retention times.
#' @param smat A matrix of (first-pass) apex scan numbers.
#' @param cmat A matrix of spectrum counts.
#' @param mbr_peps A vector of \code{pep_seq_modz} sequences in the universe.
#' @param mbr_sps A vector of peptide species in the universe.
#' @param sp_centers The centers of log2FC for each species; names: species,
#'   values: log2FC.
#' @param mbr_ret_tol The tolerance in MBR retention time in seconds.
#' @param step The step size in binning retention times.
#' @param err_log2r Not yet used. Error tolerance in log2FC.
#' @param new_na_species A replace value for NA species.
#' @param rt_tol Error tolerance in retention times.
#' @param reference_spec_counts Logical; reference protein spectrum counts or
#'   not.
pepLFQ_x1 <- function (basenames, are_refs, are_smpls, xsmat, ysmat, tsmat, ssmat, 
                    ymat, tmat, smat, cmat, mbr_peps, mbr_sps, sp_centers, 
                    new_na_species = ".other", rt_tol = 25, mbr_ret_tol = 25, 
                    reference_spec_counts = FALSE, step = 5E-3, err_log2r = .25)
  
{
  # check again if multiple fractions under a set of TMTset1_LCMSinj1_Peptide_N
  # need to compile retention times, moverzs and intensities across fractions...
  
  n_row <- nrow(tsmat)
  n_col <- ncol(tsmat)
  yout  <- tout <- sout <- rep_len(list(vector("numeric", n_row)), n_col)
  null_dbl <- rep_len(NA_real_, n_row)
  null_int <- rep_len(NA_integer_, n_row)
  
  # distances of reference log2FCs from the "first pass" to cents
  cents <- sp_centers[mbr_sps]
  dys1  <- calc_mat_log2s_to_refs(ymat, are_smpls, are_refs) - cents
  
  if (reference_spec_counts) {
    cs1 <- calc_mat_log2s_to_refs(cmat, are_smpls, are_refs)
  }
  else {
    cs1 <- rep_len(Inf, n_col)
  }
  rownames(cmat) <- colnames(cmat) <- NULL
  
  for (i in 1:n_col) {
    if (is.na(mbr_sps[[i]])) {
      yout[[i]] <- null_dbl
      tout[[i]] <- null_dbl
      sout[[i]] <- null_int
      
      next
    }
    
    ## (1) collapse bins of T, Y and S values
    # i <- which(mbr_peps == "LLDVVHPAAK@3")
    ans <- collapseSTY(
      xs = tsmat[, i], ys = ysmat[, i], zs = ssmat[, i], lwr = 10, step = step)
    anst <- ans[["x"]]
    ansy <- ans[["y"]]
    anss <- ans[["z"]]
    unv  <- ans[["u"]]
    
    if (!(nc <- ncol(anst))) {
      yout[[i]] <- null_dbl
      tout[[i]] <- null_dbl
      sout[[i]] <- null_int
      
      next
    }
    
    # border effect: may disable since values can traverse for every +2
    # may check +2 +2 +2 senario...
    if (nc > 1L && length(us <- which(diff(unv) == 2L))) {
      for (u in us) {
        rows <- which(is.na(ansy[, u2 <- u + 1L]) & !is.na(ansy[, u]))
        
        if (length(rows)) {
          ansy[rows, u2] <- ansy[rows, u]
          anst[rows, u2] <- anst[rows, u]
          anss[rows, u2] <- anss[rows, u]
          # ansy[rows, u] <- NA_real_
          # anst[rows, u] <- NA_real_
          # anss[rows, u] <- NA_real_
        }
      }
    }
    
    spc   <- cents[[i]]
    
    # for cross-referencing
    dy1i  <- dys1[[i]]
    dc1i  <- cs1[[i]] - spc
    ady1i <- abs(dy1i)
    adc1i <- abs(dc1i)
    
    if (nc == 1L) {
      yout[[i]] <- ansy[, 1L]
      tout[[i]] <- anst[, 1L]
      sout[[i]] <- anss[, 1L]
      
      next
      
      if (FALSE) {
        p      <- 1L
        ysi    <- ansy[, p]
        ysmpls <- ysi[are_smpls]
        yrefs  <- ysi[are_refs]
        ybar0  <- mean(yrefs, na.rm = TRUE)
        di     <- mean(log2(ysmpls / ybar0) - spc, na.rm = TRUE)
        adi    <- abs(di)
        bst    <- find_best_lfqvec(d1 = di, d2 = dy1i, dc = dc1i)
        
        if (isTRUE(bst == 1L)) {
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        else if (isTRUE(bst == 2L)) {
          yout[[i]] <- ymat[, i]
          tout[[i]] <- tmat[, i]
          sout[[i]] <- smat[, i]
        }
        else if (isTRUE(bst == 3L) && !(is.nan(di) || is.na(di))) {
          ybar1 <- mean(ysmpls / 2^(dc1i + spc), na.rm = TRUE)
          ysi[are_refs] <- yrefs * ybar1 / ybar0
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        else {
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        
        next
      }
    }
    
    ## (2) remove retention outliers
    for (j in 1:nc) {
      tvals <- anst[, j]
      # reduce the value to be consistent with step...
      tbads <- abs(tvals - median(tvals, na.rm = TRUE)) > rt_tol
      ibads <- .Internal(which(tbads))
      
      if (length(ibads)) {
        anst[ibads, j] <- NA_real_
        ansy[ibads, j] <- NA_real_
        anss[ibads, j] <- NA_real_
      }
    }
    
    ## (3) calculate log2FCs
    ys_smpls <- ansy[are_smpls, , drop = FALSE]
    ys_refs  <- ansy[are_refs,  , drop = FALSE]
    ymeans   <- colMeans(ys_refs, na.rm = TRUE)
    
    # all are without reference values
    if (all(is.nan(ymeans))) {
      p <- .Internal(which.min(calc_mat_sds(anst[are_smpls, , drop = FALSE])))
      
      if (length(p)) {
        ysi <- ansy[, p]
        yout[[i]] <- ysi
        tout[[i]] <- anst[, p]
        sout[[i]] <- anss[, p]
        
        if (FALSE) {
          ysi <- ansy[, p]
          # all NA?
          adi <- abs(calc_vec_log2s_to_refs(ysi, are_smpls, are_refs) - spc)
          bst <- .Internal(which.min(c(adi, ady1i)))
          
          if (isTRUE(bst == 1L)) {
            yout[[i]] <- ysi
            tout[[i]] <- anst[, p]
            sout[[i]] <- anss[, p]
          }
          else if (isTRUE(bst == 2L)) {
            yout[[i]] <- ymat[, i]
            tout[[i]] <- tmat[, i]
            sout[[i]] <- smat[, i]
          }
          else {
            yout[[i]] <- ysi
            tout[[i]] <- anst[, p]
            sout[[i]] <- anss[, p]
          }
        }
      }
      else {
        # e.g. all with single sample values (may be due to the boundary effect)
        yout[[i]] <- ansy[, 1L]
        tout[[i]] <- anst[, 1L]
        sout[[i]] <- anss[, 1L]
        
        if (FALSE) {
          # e.g. all with single sample values (may be due to the boundary effect)
          if (is.na(ady1i)) {
            yout[[i]] <- ansy[, 1]
            tout[[i]] <- anst[, 1]
            sout[[i]] <- anss[, 1]
          }
          else {
            yout[[i]] <- ymat[, i]
            tout[[i]] <- tmat[, i]
            sout[[i]] <- smat[, i]
          }
        }
      }
      
      next
    }
    
    rmat <- sweep(ys_smpls, 2, ymeans, "/")
    rmat <- log2(rmat)
    
    ## (4) find the best Y column
    ds  <- colMeans(rmat, na.rm = TRUE) - spc
    ads <- abs(ds)
    p   <- .Internal(which.min(ads))
    
    if (length(p)) {
      # if any NA in ansy[, p] && and log2FC outside the contour...
      di <- ds[[p]]
      yout[[i]] <- ansy[, p]
      tout[[i]] <- anst[, p]
      sout[[i]] <- anss[, p]
      
      next
      
      if (FALSE) {
        di <- ds[[p]]
        bst <- find_best_lfqvec(d1 = di, d2 = dy1i, dc = dc1i)
        
        if (isTRUE(bst == 1L)) {
          yout[[i]] <- ansy[, p]
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        else if (isTRUE(bst == 2L)) {
          yout[[i]] <- ymat[, i]
          tout[[i]] <- tmat[, i]
          sout[[i]] <- smat[, i]
        }
        else if (isTRUE(bst == 3L) && !(is.nan(di) || is.na(di))) {
          ysi    <- ansy[, p]
          ysmpls <- ysi[are_smpls]
          yrefs  <- ysi[are_refs]
          ybar0  <- mean(yrefs, na.rm = TRUE)
          ybar1  <- ysmpls / 2^(dc1i + spc)
          ysi[are_refs] <- yrefs * ybar1 / ybar0
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        else {
          yout[[i]] <- ansy[, p]
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        
        next
      }
    }
    
    # e.g., one column exclusive samples, the other exclusive references
    p <- .Internal(which.min(calc_mat_sds(ansy)))
    
    if (length(p)) {
      ysi <- ansy[, p]
      yout[[i]] <- ysi
      tout[[i]] <- anst[, p]
      sout[[i]] <- anss[, p]
      
      next
      
      if (FALSE) {
        ysi <- ansy[, p]
        adi  <- abs(calc_vec_log2s_to_refs(ysi, are_smpls, are_refs) - spc)
        bst <- .Internal(which.min(c(adi, ady1i)))
        
        if (isTRUE(bst == 1L)) {
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        else {
          yout[[i]] <- ymat[, i]
          tout[[i]] <- tmat[, i]
          sout[[i]] <- smat[, i]
        }
        
        next
      }
    }
    
    # e.g. all are with references only, or all single values...
    # may be use the column closest in RT to MS2 scan
    p <- .Internal(which.max(colSums(ansy, na.rm = TRUE)))
    
    if (length(p)) {
      ysi <- ansy[, p]
      yout[[i]] <- ysi
      tout[[i]] <- anst[, p]
      sout[[i]] <- anss[, p]
      
      next
      
      if (FALSE) {
        ysi <- ansy[, p]
        # all NA? May be just check ymat[, i] see if is due to border effect 
        adi <- abs(calc_vec_log2s_to_refs(ysi, are_smpls, are_refs) - spc)
        bst <- .Internal(which.min(c(adi, ady1i)))
        
        if (isTRUE(bst == 1L)) {
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        else if (isTRUE(bst == 2L)) {
          yout[[i]] <- ymat[, i]
          tout[[i]] <- tmat[, i]
          sout[[i]] <- smat[, i]
        }
        else {
          yout[[i]] <- ysi
          tout[[i]] <- anst[, p]
          sout[[i]] <- anss[, p]
        }
        
        next
      }
    }
    
    yout[[i]] <- ansy[, 1]
    tout[[i]] <- anst[, 1]
    sout[[i]] <- anss[, 1]
  }
  
  out <- list(y = yout, t = tout, s = sout)
}

