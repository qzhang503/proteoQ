#' Finds the numeric difference in ppm.
#' 
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The difference between \eqn{x} and \eqn{y} in ppm.
find_ppm_error <- function (x = 1000, y = 1000.01) {
  (y - x)/y * 1E6
}


#' Finds the error range of a number.
#'
#' Assumes \eqn{x} is position.
#' 
#' @param x A numeric value.
#' @param ppm Numeric; the ppm allowed from \code{x}.
#' @return The lower and the upper bound to \eqn{x} by \eqn{ppm}.
find_mass_error_range <- function (x = 500, ppm = 20) {
  d <- x * ppm/1E6
  range(x-d, x+d)
}


#' Find the cut-points of MS1 masses for binning.
#'
#' The cut-points will be used as the boundary for data binning. Note that the
#' upper bound is open.
#'
#' @param from Numeric; the starting MS1 mass.
#' @param to Numeric; the ending MS1 mass.
#' @param ppm Numeric; the ppm for data binning.
find_ms1_cutpoints <- function (from = 350, to = 1700, ppm = 20) {
  d <- ppm/1e6
  n <- ceiling(log(to/from)/log(1+d))
  
  x <- vector("numeric", n)
  x[1] <- from
  
  for (i in seq_len(n-1)) {
    x[i+1] <- x[i] * (1 + d)
  }
  
  x
}


#' Calculates the frame number for an experimental mass by intervals.
#' 
#' Needs correct \code{from}.
#' @param mass Numeric; an MS1 mass.
#' @inheritParams find_ms1_cutpoints
find_ms1_interval <- function (from = 350, mass = 1714.81876, ppm = 20) {
  d <- ppm/1e6
  ceiling(log(mass/from)/log(1+d))
}


#' Helper: separates theoretical peptides into mass groups.
#' 
#' @param peps A list of theoretical peptides with masses.
#' @inheritParams binTheoPeps
bin_theopeps <- function (peps, min_mass = 350, max_mass = 1700, ppm = 20) {
  ps <- find_ms1_cutpoints(min_mass, max_mass, ppm)
  
  purrr::imap(peps, ~ {
    prot_peps <- .x
    
    prot_peps %>% 
      findInterval(ps) %>% 
      bind_cols(pep_seq = names(prot_peps), mass = prot_peps, frame = .) %>% 
      dplyr::mutate(prot_acc = .y)
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(frame, pep_seq, prot_acc) %>% 
    split(., .$frame, drop = FALSE)
}


#' Separates theoretical peptides into mass groups.
#'
#' @param res Lists of data containing theoretical peptides and masses from
#'   \link{readRDS}.
#' @param min_mass Numeric; the minimum MS1 mass.
#' @param max_mass Numeric; the maximum MS1 mass.
#' @param ppm Numeric; the eror tolerance of MS1 mass in ppm.
#' @examples
#' \donttest{
#' res <- readRDS("~/proteoQ/dbs/fasta/uniprot/pepmass/uniprot_hs_2020_05_2miss.rds")
#' theopeps <- binTheoPeps(res)
#' }
#' @return Lists of theoretical peptides binned by MS1 masses. The lists
#'   correspond to the lists of \code{res}.
binTheoPeps <- function (res, min_mass = 516.24046, max_mass = 10000, 
                         ppm = 20) {
  data <- res %>% 
    purrr::map(attributes) %>% 
    purrr::map(`[[`, "data")

  message("Bin MS1 peptide masses (theoretical).")
  
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(cl, list("find_ms1_cutpoints"), 
                          envir = env_where("find_ms1_cutpoints"))
  
  out <- parallel::clusterMap(cl, bin_theopeps, data, 
                              MoreArgs = list(min_mass, max_mass, ppm))

  parallel::stopCluster(cl)
  
  invisible(out)
}


#' Add the MS file name to timsTOF mgf files.
#'
#' Uses base R readLines (to avoid memory failure with stringi).
#' 
#' @inheritParams readMGF
#' @import stringi
#' @examples
#' \donttest{
#' mgf_queries <- readMGF()
#' }
readMGF_br <- function (filepath = "~/proteoQ/mgf_test", min_mass = 516.2405, 
                     ret_range = c(0, Inf)) {
  # path <- gsub("(^.*/)[^/]*\\.mgf$", "\\1", filepath)
  # files <- gsub("^.*/([^/]*\\.mgf$)", "\\1", filepath)
  
  filelist <- list.files(path = file.path(filepath), 
                         pattern = "^.*\\.mgf$")
  
  if (purrr::is_empty(filelist)) {
    stop("No mgf files under ", filepath, 
         call. = FALSE)
  }
  
  lines <- purrr::map(filelist, ~ readLines(file.path(filepath, .x))) %>% 
    purrr::reduce(`c`)
  
  # MS2 ions
  begins <- grep("^SCANS", lines) + 1
  ends <- which(lines == "END IONS") - 1
  
  ms2s <- purrr::map2(begins, ends, ~ lines[.x : .y]) %>% 
    purrr::map(stringr::str_split_fixed, " ", n = 2)
  
  ms2_moverzs <- ms2s %>% 
    purrr::map(~ .x[, 1] %>% as.numeric()) 
  
  ms2_ints <- ms2s %>% 
    purrr::map(~ .x[, 2] %>% as.numeric()) 
  
  rm(begins, ends, ms2s)
  
  # MS1 ions
  ms1s <- lines %>% 
    .[grepl("PEPMASS", .)] %>% 
    gsub("PEPMASS=", "", .) %>% 
    purrr::map(stringr::str_split, " ", simplify = TRUE)
  
  ms1_moverzs <- ms1s %>% 
    purrr::map(~ .x[, 1] %>% as.numeric())
  
  ms1_ints <- ms1s %>% 
    purrr::map(~ .x[, 2] %>% as.numeric())
  
  rm(ms1s)
  
  # Others
  scan_titles <- lines %>% 
    .[grepl("TITLE", .)] %>% 
    gsub("TITLE=", "", .) %>% 
    as.list()
  
  scan_nums <- lines %>% 
    .[grepl("SCANS", .)] %>% 
    gsub("SCANS=", "", .) %>% 
    as.numeric() %>% 
    as.list()
  
  ret_times <- lines %>% 
    .[grepl("RTINSECONDS", .)] %>% 
    gsub("RTINSECONDS=", "", .) %>% 
    as.numeric() %>% 
    as.list()
  
  ms1_charges <- lines %>% 
    .[grepl("CHARGE", .)] %>% 
    gsub("CHARGE=", "", .) %>% 
    as.list()
  
  # --- MS1 neutral masses
  proton <- 1.00727647
  
  charges <- ms1_charges %>% 
    purrr::map(stringi::stri_reverse) %>% 
    purrr::map(as.numeric)
  
  ms1_masses <- purrr::map2(ms1_moverzs, charges, ~ {
    .x * .y - .y * proton # + proton
  })
  
  frames <- purrr::map(ms1_masses, find_ms1_interval, from = min_mass)
  
  # subset
  rows <- (ret_times >= ret_range[1] & ret_times <= ret_range[2])
  
  scan_titles <- scan_titles[rows]
  ms1_moverzs <- ms1_moverzs[rows]
  ms1_masses <- ms1_masses[rows]
  ms1_ints <- ms1_ints[rows]
  ms1_charges <- ms1_charges[rows]
  ret_times <- ret_times[rows]
  scan_nums <- scan_nums[rows]
  ms2_moverzs <- ms2_moverzs[rows]
  ms2_ints <- ms2_ints[rows]
  frames <- frames[rows]
  
  # Ordered outputs
  orders <- order(unlist(ms1_masses))
  
  out <- tibble(scan_title = scan_titles[orders], 
                ms1_moverz = ms1_moverzs[orders], 
                ms1_mass = ms1_masses[orders], 
                ms1_int = ms1_ints[orders], 
                ms1_charge = ms1_charges[orders],
                ret_time = ret_times[orders], 
                scan_num = scan_nums[orders],
                ms2_moverz = ms2_moverzs[orders], 
                ms2_int = ms2_ints[orders], 
                frame = frames[orders])
}



#' Add the MS file name to timsTOF mgf files.
#'
#' @param filepath The file path to a list of mgf files.
#' @param min_mass Numeric; the minimum mass of MS1 species. The value needs to
#'   match the one in  \link{binTheoPeps}.
#' @param ret_range The range of retention time in seconds.
#' @import stringi
#' @examples
#' \donttest{
#' mgf_queries <- readMGF()
#' }
readMGF <- function (filepath = "~/proteoQ/mgf_test", min_mass = 516.2405, 
                     ret_range = c(0, Inf)) {
  
  filelist <- list.files(path = file.path(filepath), 
                         pattern = "^.*\\.mgf$")
  
  if (purrr::is_empty(filelist)) {
    stop("No mgf files under ", filepath, 
         call. = FALSE)
  }
  
  # lines <- purrr::map(filelist, ~ stri_read_lines(file.path(filepath, .x))) %>% 
  #   purrr::reduce(`c`)

  lines <- purrr::map(filelist, ~ readLines(file.path(filepath, .x))) %>% 
    purrr::reduce(`c`)

  # MS2 ions
  begins <- which(stri_startswith_fixed(lines, "SCANS")) + 1
  ends <- which(stri_endswith_fixed(lines, "END IONS")) - 1
  
  ms2s <- purrr::map2(begins, ends, ~ lines[.x : .y]) %>% 
    purrr::map(stri_split_fixed, " ", n = 2, simplify = TRUE) 

  ms2_moverzs <- ms2s %>% 
    purrr::map(~ .x[, 1] %>% as.numeric()) 
  
  ms2_ints <- ms2s %>% 
    purrr::map(~ .x[, 2] %>% as.numeric()) 
  
  rm(begins, ends, ms2s)

  # MS1 ions
  ms1s <- lines %>% 
    .[grepl("^PEPMASS", .)] %>% 
    stri_replace_first_fixed("PEPMASS=", "") %>% 
    purrr::map(stri_split_fixed, " ", n = 2, simplify = TRUE)

  ms1_moverzs <- ms1s %>% 
    purrr::map(~ .x[, 1] %>% as.numeric())
  
  ms1_ints <- ms1s %>% 
    purrr::map(~ .x[, 2] %>% as.numeric())
  
  rm(ms1s)
  
  # Others
  scan_titles <- lines %>% 
    .[grepl("^TITLE", .)] %>% 
    stri_replace_first_fixed("TITLE=", "") %>% 
    as.list()
  
  scan_nums <- lines %>% 
    .[grepl("^SCANS", .)] %>% 
    stri_replace_first_fixed("SCANS=", "") %>% 
    as.numeric() %>% 
    as.list()
  
  ret_times <- lines %>% 
    .[grepl("^RTINSECONDS", .)] %>% 
    stri_replace_first_fixed("RTINSECONDS=", "") %>% 
    as.numeric() %>% 
    as.list()

  ms1_charges <- lines %>% 
    .[grepl("^CHARGE", .)] %>% 
    stri_replace_first_fixed("CHARGE=", "") %>% 
    as.list()

  # MS1 neutral masses
  proton <- 1.00727647
  
  charges <- ms1_charges %>% 
    purrr::map(stri_reverse) %>% 
    purrr::map(as.numeric)
  
  ms1_masses <- purrr::map2(ms1_moverzs, charges, ~ {
    .x * .y - .y * proton 
  })
  
  frames <- purrr::map(ms1_masses, find_ms1_interval, from = min_mass)
  
  # Subsetting
  rows <- (ret_times >= ret_range[1] & ret_times <= ret_range[2])
  
  scan_titles <- scan_titles[rows]
  ms1_moverzs <- ms1_moverzs[rows]
  ms1_masses <- ms1_masses[rows]
  ms1_ints <- ms1_ints[rows]
  ms1_charges <- ms1_charges[rows]
  ret_times <- ret_times[rows]
  scan_nums <- scan_nums[rows]
  ms2_moverzs <- ms2_moverzs[rows]
  ms2_ints <- ms2_ints[rows]
  frames <- frames[rows]
  
  # Ordering
  orders <- order(unlist(ms1_masses))

  out <- tibble(scan_title = scan_titles[orders], 
                ms1_moverz = ms1_moverzs[orders], 
                ms1_mass = ms1_masses[orders], 
                ms1_int = ms1_ints[orders], 
                ms1_charge = ms1_charges[orders],
                ret_time = ret_times[orders], 
                scan_num = scan_nums[orders],
                ms2_moverz = ms2_moverzs[orders], 
                ms2_int = ms2_ints[orders], 
                frame = frames[orders])
}


#' Searches MGFs in a frame.
#'
#' It reads and searches one frame of MGFs against \code{theopeps}. The frame
#' number links experimental and theoretical spectra by MS1 Masses.
#'
#' @param theopeps Binned theoretical peptides corresponding to an i-th
#'   \code{aa_masses}.
#' @param mgf_frames MGFs in frames. Each frame contains one to multiple MGFs
#'   whose MS1 masses are in the same interval.
#' @param minn_ms2 Integer; the minimum number of MS2 ions for consideration as
#'   a hit.
#' @param ppm_ms1 The mass tolerance of MS1 species.
#' @param ppm_ms2 The mass tolerance of MS2 species.
#' @inheritParams calc_ms2ionseries
search_mgf_frames <- function (theopeps, aa_masses, mgf_frames, 
                               type = "by", 
                               maxn_vmods_per_pep = 5, 
                               maxn_sites_per_vmod = 3, 
                               maxn_vmods_sitescombi_per_pep = 32, 
                               minn_ms2 = 3, 
                               ppm_ms1 = 20, ppm_ms2 = 25, digits = 5) {
  len <- length(mgf_frames)
  out <- vector("list", len) 
  
  ## --- initialization ---
  mgfs_cr <- mgf_frames[[1]]
  frame <- mgfs_cr[["frame"]][[1]]

  theos_bf_ms1 <- theopeps[[as.character(frame-1)]]
  theos_cr_ms1 <- theopeps[[as.character(frame)]]
  
  theomasses_bf_ms1 <- theos_bf_ms1$mass
  theomasses_cr_ms1 <- theos_cr_ms1$mass
  
  theos_bf_ms2 <- theos_bf_ms1$pep_seq %>% 
    map(calc_ms2ionseries, 
        aa_masses = aa_masses, 
        type = type, 
        maxn_vmods_per_pep = maxn_vmods_per_pep, 
        maxn_sites_per_vmod = maxn_sites_per_vmod, 
        maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
        digits = digits)

  theos_cr_ms2 <- theos_cr_ms1$pep_seq %>% 
    map(
      calc_ms2ionseries, 
      aa_masses = aa_masses, 
      type = type, 
      maxn_vmods_per_pep = maxn_vmods_per_pep, 
      maxn_sites_per_vmod = maxn_sites_per_vmod, 
      maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
      digits = digits)
  
  ## --- iterating ---
  for (i in seq_len(len)) {
    exptmasses_ms1 <- mgfs_cr[["ms1_mass"]]
    exptmoverzs_ms2 <- mgfs_cr[["ms2_moverz"]] # matches with higher charge states later
    
    theos_af_ms1 <- theopeps[[as.character(frame+1)]]
    theomasses_af_ms1 <- theos_af_ms1$mass
    
    theos_af_ms2 <- theos_af_ms1$pep_seq %>% 
      map(calc_ms2ionseries, 
          aa_masses = aa_masses, 
          type = type, 
          maxn_vmods_per_pep = maxn_vmods_per_pep, 
          maxn_sites_per_vmod = maxn_sites_per_vmod, 
          maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
          digits = digits)

    # each `out` for the results of multiple mgfs in one frame
    out[[i]] <- map2(exptmasses_ms1, exptmoverzs_ms2, search_mgf, 
                     theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
                     theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                     minn_ms2, ppm_ms1, ppm_ms2)
    
    # advance to the next frame
    if (i == len) {
      break
    }
    
    mgfs_cr <- mgf_frames[[i+1]]
    new_frame <- mgfs_cr[["frame"]][[1]]
    
    if (new_frame == (frame+1)) {
      theos_bf_ms1 <- theos_cr_ms1
      theos_cr_ms1 <- theos_af_ms1
      
      theomasses_bf_ms1 <- theomasses_cr_ms1
      theomasses_cr_ms1 <- theomasses_af_ms1
      
      theos_bf_ms2 <- theos_cr_ms2
      theos_cr_ms2 <- theos_af_ms2
    } else if (new_frame == (frame+2)) {
      theos_bf_ms1 <- theos_af_ms1
      theos_cr_ms1 <- theopeps[[as.character(new_frame)]]
      
      theomasses_bf_ms1 <- theomasses_af_ms1
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_bf_ms2 <- theos_af_ms2
      
      theos_cr_ms2 <- theos_cr_ms1$pep_seq %>% 
        map(calc_ms2ionseries, 
            aa_masses = aa_masses, 
            type = type, 
            maxn_vmods_per_pep = maxn_vmods_per_pep, 
            maxn_sites_per_vmod = maxn_sites_per_vmod, 
            maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
            digits = digits)
    } else {
      theos_bf_ms1 <- theopeps[[as.character(new_frame-1)]]
      theos_cr_ms1 <- theopeps[[as.character(new_frame)]]
      
      theomasses_bf_ms1 <- theos_bf_ms1$mass
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_bf_ms2 <- theos_bf_ms1$pep_seq %>% 
        map(calc_ms2ionseries, 
            aa_masses = aa_masses, 
            type = type, 
            maxn_vmods_per_pep = maxn_vmods_per_pep, 
            maxn_sites_per_vmod = maxn_sites_per_vmod, 
            maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
            digits = digits)
      
      theos_cr_ms2 <- theos_cr_ms1$pep_seq %>% 
        map(calc_ms2ionseries, 
            aa_masses = aa_masses, 
            type = type, 
            maxn_vmods_per_pep = maxn_vmods_per_pep, 
            maxn_sites_per_vmod = maxn_sites_per_vmod, 
            maxn_vmods_sitescombi_per_pep = maxn_vmods_sitescombi_per_pep, 
            digits = digits)
    }
    
    frame <- new_frame
  }
  
  out
}


#' Searches a single MGF.
#'
#' @param expt_mass_ms1 Numeric; the experimental MS1 mass.
#' @param expt_moverz_ms2 A numeric list; the experimental MS2 m/z's.
#' @param theomasses_bf_ms1 Numeric vector; the theoretical MS1 masses at the
#'   preceding \code{-1} frame.
#' @param theomasses_cr_ms1 Numeric vector; the theoretical MS1 masses at the
#'   current frame.
#' @param theomasses_af_ms1 Numeric vector; the theoretical MS1 masses at the
#'   following \code{+1} frame.
#' @param theos_bf_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   preceding \code{-1} frame.
#' @param theos_cr_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   current frame.
#' @param theos_af_ms2 Numeric vector; the theoretical MS2 m/z's at the
#'   following \code{+1} frame.
#' @inheritParams search_mgf_frames
#' @import dplyr
search_mgf <- function (expt_mass_ms1, expt_moverz_ms2, 
                        theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
                        theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                        minn_ms2 = 3, ppm_ms1 = 20, ppm_ms2 = 25) {

  # --- subset from the `before` and the `after` by MS1 mass tolerance 
  mass_ranges <- find_mass_error_range(expt_mass_ms1, ppm_ms1)
  bf_allowed <- which(theomasses_bf_ms1 >= mass_ranges[1])
  af_allowed <- which(theomasses_af_ms1 <= mass_ranges[2])

  # --- find MS2 matches ---
  x <- purrr::map(c(theos_bf_ms2[bf_allowed], 
                    theos_cr_ms2, 
                    theos_af_ms2[af_allowed]), ~ {
    d <- outer(expt_moverz_ms2, .x, "find_ppm_error")
    row_cols <- which(abs(d) <= ppm_ms2, arr.ind = TRUE)
    bind_cols(expt = expt_moverz_ms2[row_cols[, 1]], theo = .x[row_cols[, 2]])
  }) %>% 
    `names<-`(c(names(bf_allowed), 
                names(theomasses_cr_ms1), 
                names(af_allowed)))
  
  rows <- map(x, nrow) > minn_ms2
  x <- x[rows]
  
  rm(mass_ranges, bf_allowed, af_allowed, rows)
  
  invisible(x)
}







foo <- function (filepath = "~/proteoQ/mgf_test") {
  res <- readRDS("~/proteoQ/dbs/fasta/uniprot/pepmass/uniprot_hs_2020_05_2miss.rds")
  
  # [2] "TMT6plex (N-term)"
  # [3] "Acetyl (Protein N-term)"
  aa_masses_all <- res %>% 
    purrr::map(~ {
      attr(.x, "data") <- NULL
      .x
    })

  tempdata <- res %>% 
    purrr::map(attributes) %>% 
    purrr::map(`[[`, "data") %>% 
    unlist(use.names = FALSE)
  min_mass <- min(tempdata)
  max_mass <- max(tempdata)
  rm(tempdata)

  binned_theopeps <- binTheoPeps(res, min_mass = min_mass, max_mass = max_mass)
  # saveRDS(binned_theopeps, "~/proteoQ/mgf_test/binned_theopeps.rds")
  # binned_theopeps <- readRDS("~/proteoQ/mgf_test/binned_theopeps.rds")
  
  res <- res[3]
  aa_masses_all <- aa_masses_all[3]
  binned_theopeps <- binned_theopeps[3]
  rm(res)
  
  # --- read experimental data
  # mgf_queries <- readMGF(min_mass = min_mass)
  # saveRDS(mgf_queries, "~/proteoQ/mgf_test/mgf_queries.rds")
  mgf_queries <- readRDS("~/proteoQ/mgf_test/mgf_queries.rds")
  
  # ret_time 11050
  # mgf_queries <- mgf_queries[39667:49667, ]
  
  # split by frame
  mgf_frames <- mgf_queries %>% 
    group_by(frame) %>% 
    group_split()
  
  out <- purrr::map2(binned_theopeps, aa_masses_all, search_mgf_frames, mgf_frames) %>% 
    flatten() %>% 
    map2(mgf_frames, ., ~ {
      .x$out <- .y
      .x
    })

}