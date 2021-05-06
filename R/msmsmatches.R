#' Finds the indexes of top-n entries without re-ordering.
#'
#' @param x A numeric vector.
#' @param n The number of top entries to keep.
#' @return The indexes of the top-n entries.
which_topx <- function(x, n = 50, ...) {
  len <- length(x)
  p <- len - n
  
  if (p  <= 0) return(seq_along(x))
  
  xp <- sort(x, partial = p, ...)[p]
  
  which(x > xp)
}


#' Finds the top-n entries without re-ordering.
#' 
#' @inheritParams which_topx
#' @return The top-n entries.
topx <- function(x, n = 50, ...) {
  len <- length(x)
  p <- len - n
  
  if (p  <= 0) return(x)
  
  xp <- sort(x, partial = p, ...)[p]
  
  x[x > xp]
}


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
#' Assumes \eqn{x} is positive without checking.
#' 
#' @param x A numeric value.
#' @param ppm Numeric; the ppm allowed from \code{x}.
#' @return The lower and the upper bound to \eqn{x} by \eqn{ppm}.
find_mass_error_range <- function (x = 500, ppm = 20) {
  d <- x * ppm/1E6
  range(x-d, x+d)
}


#' Finds the cut-points of MS1 masses for binning.
#'
#' The cut-points will be used as the boundary in data binning. Note that the
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


#' Calculates the frame number for an experimental MS1 mass by intervals.
#' 
#' Needs correct \code{from}.
#' @param mass Numeric; a list of MS1 masses.
#' @inheritParams find_ms1_cutpoints
find_ms1_interval <- function (mass = 1714.81876, from = 350, ppm = 20) {
  d <- ppm/1e6
  ceiling(log(unlist(mass, recursive = FALSE, use.names = FALSE)/from)/log(1+d))
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
      dplyr::bind_cols(pep_seq = names(prot_peps), 
                       mass = prot_peps, frame = .) %>% 
      dplyr::mutate(prot_acc = .y)
  })
}


#' Separates theoretical peptides into mass groups.
#'
#' @param res Lists of data containing theoretical peptides and masses from
#'   \link{readRDS}.
#' @param min_mass Numeric; the minimum MS1 mass.
#' @param max_mass Numeric; the maximum MS1 mass.
#' @param ppm Numeric; the eror tolerance of MS1 mass in ppm.
#' @param out_path The output path.
#' @examples
#' \donttest{
#' res <- readRDS("~/proteoQ/dbs/fasta/uniprot/pepmass/uniprot_hs_2020_05_2miss.rds")
#' theopeps <- binTheoPeps(res)
#' }
#' @return Lists of theoretical peptides binned by MS1 masses. The lists
#'   correspond to the lists of \code{res}.
#' @import parallel
#' @export
binTheoPeps <- function (res, min_mass = 516.24046, max_mass = 10000, ppm = 20, 
                         out_path = file.path("~/proteoQ/outs", 
                                              "pepmasses/binned_theopeps.rds")) {
  res <- res %>% 
    purrr::map(attributes) %>% 
    purrr::map(`[[`, "data")
  
  message("Binning MS1 peptide masses (theoretical).")

  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  parallel::clusterExport(cl, list("find_ms1_cutpoints"), 
                          envir = environment(proteoQ:::find_ms1_cutpoints))

  parallel::clusterExport(cl, list("min_mass", "max_mass", "ppm"), 
                          envir = environment())
  
  out <- vector("list", length(res))
  
  for (i in seq_along(out)) {
    res[[i]] <- chunksplit(res[[i]], n_cores)
    
    out[[i]] <- parallel::clusterApplyLB(cl, res[[i]], bin_theopeps, 
                                         min_mass, max_mass, ppm) %>% 
      purrr::flatten()
    
    out[[i]] <- out[[i]] %>% 
      dplyr::bind_rows() %>% 
      dplyr::arrange(frame, pep_seq, prot_acc) %>% 
      split(., .$frame, drop = FALSE)
  }
  
  rm(res)
  
  local({
    out_dir <- gsub("(^.*/).*$", "\\1", out_path)
    out_nms <- gsub("^.*/(.*)\\.[^\\.].*$", "\\1", out_path) %>% 
      paste(seq_along(out), sep = "_") %>% 
      paste0(".rds")
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    purrr::walk2(out, out_nms, 
                 ~ saveRDS(.x, file.path(out_dir, .y)))
  })

  parallel::stopCluster(cl)

  invisible(out)
}


#' Finds directory.
#' 
#' @param dir A file path.
find_dir <- function (dir) {
  dir_1 <- fs::path_expand_r(dir)
  dir_2 <- fs::path_expand(dir)
  
  if (fs::dir_exists(dir_1)) {
    dir <- dir_1
  } else if (fs::dir_exists(dir_2)) {
    dir <- dir_2
  } else {
    dir <- NULL
  }
  
  dir
}


#' Helper in processing MGF entries in chunks.
#' 
#' @inheritParams proc_mgf_chunks
#' @import stringi
proc_mgfs <- function (lines, topn_ms2ions = 100, ret_range = c(0, Inf)) {
  
  # MS2 ions
  begins <- which(stri_startswith_fixed(lines, "BEGIN IONS"))
  ends <- which(stri_endswith_fixed(lines, "END IONS"))
  
  ms2s <- purrr::map2(begins, ends, ~ lines[(.x + 6) : (.y - 1)]) %>% 
    purrr::map(stri_split_fixed, " ", n = 2, simplify = TRUE) 
  
  ms2_moverzs <- ms2s %>% 
    purrr::map(~ .x[, 1] %>% as.numeric()) 
  
  ms2_ints <- ms2s %>% 
    purrr::map(~ .x[, 2] %>% as.numeric()) 
  
  lens <- purrr::map(ms2_moverzs, length)
  
  if (topn_ms2ions < Inf) {
    rows <- purrr::map(ms2_ints, which_topx, topn_ms2ions)
    ms2_ints <- purrr::map2(ms2_ints, rows, ~ .x[.y])
    ms2_moverzs <- purrr::map2(ms2_moverzs, rows, ~ .x[.y])
    rm(rows, ms2s)
  }
  
  # MS1 ions
  ms1s <- lines[begins+2] %>% 
    stri_replace_first_fixed("PEPMASS=", "") %>% 
    purrr::map(stri_split_fixed, " ", n = 2, simplify = TRUE)
  
  ms1_moverzs <- ms1s %>% 
    purrr::map(~ .x[, 1] %>% as.numeric())
  
  ms1_ints <- ms1s %>% 
    purrr::map(~ .x[, 2] %>% as.numeric())
  
  rm(ms1s)
  
  # Others
  scan_titles <- lines[begins+1] %>% 
    stri_replace_first_fixed("TITLE=", "") %>% 
    as.list()
  
  scan_nums <- lines[begins+5] %>% 
    stri_replace_first_fixed("SCANS=", "") %>% 
    as.numeric() %>% 
    as.list()
  
  ret_times  <- lines[begins+4] %>% 
    stri_replace_first_fixed("RTINSECONDS=", "") %>% 
    as.numeric() %>% 
    as.list()
  
  ms1_charges <- lines[begins+3] %>% 
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
  
  out <- tibble::tibble(scan_title = scan_titles, 
                        ms1_moverz = ms1_moverzs, 
                        ms1_mass = ms1_masses, 
                        ms1_int = ms1_ints, 
                        ms1_charge = ms1_charges,
                        ret_time = ret_times, 
                        scan_num = scan_nums,
                        ms2_moverz = ms2_moverzs, 
                        ms2_int = ms2_ints, 
                        ms2_n = lens)
}


#' Processes MGF entries in chunks.
#' 
#' @param lines MGF lines.
#' @inheritParams readMGF
#' @import stringi
proc_mgf_chunks <- function (lines, topn_ms2ions = 100, ret_range = c(0, Inf), 
                             filepath = file.path("~/proteoQ/mgfs/temp")) {

  basename <- gsub("\\.[^.]*$", "", filepath)
  
  begins <- which(stri_startswith_fixed(lines, "BEGIN IONS"))
  ends <- which(stri_endswith_fixed(lines, "END IONS"))
  
  af <- local({
    le <- ends[length(ends)]
    lb <- begins[length(begins)]
    
    if (lb > le) {
      af <- lines[(le+2):length(lines)]
    } else {
      af <- NULL
    }
    
    write(af, file.path(paste0(basename, "_af.mgf")))
    
    af
  })
  
  bf <- local({
    le <- ends[1]
    lb <- begins[1]
    
    if (lb > le) {
      bf <- lines[1:(le+1)]
    } else {
      bf <- NULL
    }
    
    write(bf, file.path(paste0(basename, "_bf.mgf")))
    
    bf
  })
  
  if (!is.null(af)) {
    lines <- lines[1:(begins[length(begins)]-1)]
  }
  
  if (!is.null(bf)) {
    lines <- lines[-c(1:(ends[1]+1))]
  }
  
  out <- proc_mgfs(lines, topn_ms2ions = topn_ms2ions, ret_range = ret_range)
}


#' Reads mgfs in chunks.
#' 
#' @inheritParams readMGF
#' @import stringi
read_mgf_chunks <- function (filepath = "~/proteoQ/mgfs", 
                             topn_ms2ions = 100, ret_range = c(0, Inf)) {
  
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")
  
  if (purrr::is_empty(filelist)) {
    stop("No mgf files under ", filepath, 
         call. = FALSE)
  }
  
  out <- purrr::map(filelist, ~ {
    message("Parsing '", .x, "'.")
    
    filepath <- file.path(filepath, .x)
    
    stri_read_lines(filepath) %>% 
      proc_mgf_chunks(topn_ms2ions = topn_ms2ions, 
                      ret_range = ret_range, 
                      filepath = filepath)
  }) %>% 
    dplyr::bind_rows()
  
  # adds back broken mgf entries
  afs <- local({
    afs <- list.files(path = file.path(filepath), pattern = "^.*\\_af.mgf$")
    
    idxes <- afs %>% 
      gsub("^chunk_(\\d+)_af\\.mgf", "\\1", .) %>% 
      as.integer() %>% 
      sort()
    
    paste0("chunk_", idxes, "_af.mgf") %>% 
      .[-length(.)]
  })
  
  bfs <- local({
    bfs <- list.files(path = file.path(filepath), pattern = "^.*\\_bf.mgf$")
    
    idxes <- bfs %>% 
      gsub("^chunk_(\\d+)_bf\\.mgf", "\\1", .) %>% 
      as.integer() %>% 
      sort()
    
    paste0("chunk_", idxes, "_bf.mgf") %>% 
      .[-1]
  })
  
  stopifnot(length(afs) == length(bfs))
  
  gaps <- purrr::map2(afs, bfs, ~ {
    af <- stri_read_lines(file.path(filepath, .x))
    bf <- stri_read_lines(file.path(filepath, .y))
    append(af, bf)
  }) %>% 
    unlist(use.names = FALSE) %T>% 
    write(file.path(filepath, "gaps.mgf"))
  
  local({
    nms <- list.files(path = file.path(filepath), pattern = "^.*\\_[ab]f.mgf$")
    
    if (!purrr::is_empty(nms)) {
      suppressMessages(file.remove(file.path(filepath, nms)))
    }
  })
  
  if (!is.null(gaps)) {
    out <- dplyr::bind_rows(
      out, 
      proc_mgfs(gaps, topn_ms2ions = topn_ms2ions, ret_range = ret_range)
    )
  }

  invisible(out)
}


#' Reads MGF files in chunks.
#'
#' @param filepath The file path to a list of mgf files.
#' @param min_mass Numeric; the minimum mass of MS1 species. The value needs to
#'   match the one in  \link{binTheoPeps}.
#' @param topn_ms2ions A non-negative integer; the top-n species for uses in
#'   MS2 ion searches.
#' @param ret_range The range of retention time in seconds.
#' @inheritParams matchMS
#' @import stringi
#' @examples
#' \donttest{
#' mgf_queries <- readMGF()
#' }
#' @export
readMGF <- function (filepath = "~/proteoQ/mgfs", min_mass = 516.2405, 
                     topn_ms2ions = 100, ret_range = c(0, Inf), ppm_ms1 = 20, 
                     out_path = file.path(filepath, "mgf_queries.rds")) {

  f <- function(x, pos) {
    nm <- file.path(filepath, "temp", paste0("chunk", "_", pos, ".mgf"))
    writeLines(x, nm)
  }
  
  
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")

  if (purrr::is_empty(filelist)) {
    stop("No '.mgf' files under ", filepath, 
         call. = FALSE)
  }
  
  # by mgf files
  out <- vector("list", length(filelist))
  
  for (i in seq_along(filelist)) {
    temp_dir <- local({
      temp_dir <- find_dir(file.path(filepath, "temp"))
      
      if (!is.null(temp_dir)) {
        fs::file_delete(temp_dir)
      }
      
      dir.create(file.path(filepath, "temp"), showWarnings = FALSE)
      temp_dir <- find_dir(file.path(filepath, "temp"))
    })

    message("Loading '", filelist[i], "'.")
    
    readr::read_lines_chunked(file.path(filepath, filelist[i]), 
                              SideEffectChunkCallback$new(f), 
                              chunk_size = 10000000)
    
    out[[i]] <- read_mgf_chunks(temp_dir, 
                                topn_ms2ions = topn_ms2ions, 
                                ret_range = ret_range)

    local({
      temp_dir2 <- file.path(filepath, gsub("\\.[^.]*$", "", filelist[i]))
      dir.create(temp_dir2, showWarnings = FALSE)
      temp_dir2 <- find_dir(temp_dir2)
      
      if (file_exists(temp_dir2)) {
        fs::file_delete(temp_dir2)
      }
      
      fs::file_move(temp_dir, temp_dir2)
    })
  }

  out <- out %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(ms1_mass) %>% 
    dplyr::mutate(frame = find_ms1_interval(ms1_mass, from = min_mass, 
                                            ppm = ppm_ms1)) %T>% 
    saveRDS(., out_path)

  invisible(out)
}


#' Splits data into chunks by length.
#' 
#' @param data Input data.
#' @param n_chunks The number of chunks.
#' @export
chunksplit <- function (data, n_chunks = 5) {
  if (n_chunks <= 1) return(data)
  
  len <- length(data)
  labs <- levels(cut(1:len, n_chunks))
  
  x <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
             upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))
  
  # x[x[, 1] < 0, 1] <- 1
  # x[x[, 2] > len, 2] <- len
  
  grps <- findInterval(1:len, x[, 1])
  split(data, grps)
}


#' Splits data into chunks with approximately equal sizes.
#'
#' @param nx Positive integer; an arbitrarily large number for data to be split
#'   into for estimating the cumulative sizes.
#' @inheritParams chunksplit
#' @export
chunksplitLB <- function (data, n_chunks = 5, nx = 100) {
  if (n_chunks <= 1) return(data)
  
  len <- length(data)
  
  # The finer groups by 'nx'
  grps_nx <- local({
    labsx <- levels(cut(1:len, nx))
    
    xx <- cbind(lower = floor(as.numeric( sub("\\((.+),.*", "\\1", 
                                              labsx))),
                upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", 
                                                labsx))))
    
    findInterval(1:len, xx[, 1])
  })
  
  # The equated size for a chunk
  size_chunk <- local({
    size_nx <- data %>% 
      split(., grps_nx) %>% 
      purrr::map(object.size) %>% 
      cumsum() 
    
    size_nx[length(size_nx)]/n_chunks
  })
  
  #  Intervals
  grps <- local({
    size_data <- data %>% 
      purrr::map(object.size) %>% 
      cumsum()
    
    # the position indexes
    ps <- purrr::map_dbl(1:(n_chunks-1), ~ {
      which(size_data < size_chunk * .x) %>% `[`(length(.))
    })
    
    grps <- findInterval(1:len, ps)
  })
  
  split(data, grps)
}


#' Searches MS ions.
#'
#' @inheritParams calc_ms2ionseries
#' @inheritParams calc_pepmasses
#' @inheritParams search_mgf_frames
#' @inheritParams readMGF
#' @inheritParams normPSM
#' @param mgf_path The file path to a list of MGF files.
#' @export
matchMS <- function (mgf_path = "~/proteoQ/mgfs", 
                     out_path = "~/proteoQ/outs", 
                     fasta = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta", 
                     fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)"), 
                     varmods = c("TMT6plex (N-term)", "Acetyl (Protein N-term)", 
                                 "Oxidation (M)", "Deamidated (N)", 
                                 "Gln->pyro-Glu (N-term = Q)"), 
                     enzyme = c("trypsin"), 
                     maxn_fasta_seqs = 50000,
                     maxn_vmods_setscombi = 64, 
                     maxn_vmods_per_pep = 5,
                     maxn_sites_per_vmod = 3, 
                     maxn_vmods_sitescombi_per_pep = 32, 
                     min_len = 7, max_len = 100, max_miss = 2, 
                     type_ms2ions = "by", 
                     topn_ms2ions = 100, 
                     minn_ms2 = 7, ppm_ms1 = 20, ppm_ms2 = 25, 
                     digits = 5) {
  
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  
  local({
    filelist <- list.files(path = file.path(mgf_path), pattern = "\\.mgf$")
    
    if (purrr::is_empty(filelist)) {
      stop("No `.mgf` files under ", mgf_path, 
           call. = FALSE)
    }
  })
  
  ## Indexes of modifications
  mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
    as.hexmode() %>% 
    `names<-`(c(fixedmods, varmods))
  
  ## Theoretical MS1 masses
  res <- calc_pepmasses(
    fasta = fasta, fixedmods = fixedmods, varmods = varmods, 
    index_mods = FALSE, enzyme = enzyme, 
    maxn_fasta_seqs = maxn_fasta_seqs, 
    maxn_vmods_setscombi = maxn_vmods_setscombi, 
    maxn_vmods_per_pep = maxn_vmods_per_pep, 
    maxn_sites_per_vmod = maxn_sites_per_vmod, 
    min_len = min_len, max_len = max_len, max_miss = max_miss, 
    digits = digits, 
    parallel = TRUE, 
    add_masses = TRUE, 
    out_path = file.path(out_path, "pepmasses/pepmasses.rds")
  )

  ## Mass range
  tempdata <- res %>% 
    purrr::map(attributes) %>% 
    purrr::map(`[[`, "data") %>% 
    unlist(use.names = FALSE)
  min_mass <- min(tempdata, na.rm = TRUE)
  max_mass <- max(tempdata, na.rm = TRUE)
  rm(tempdata)
  
  ## AA masses
  aa_masses_all <- res %>% 
    purrr::map(~ {
      attr(.x, "data") <- NULL
      .x
    })
  
  ## Bin theoretical peptides
  binTheoPeps(res = res, min_mass = min_mass, max_mass = max_mass, 
              ppm = ppm_ms1, 
              out_path = file.path(out_path, "pepmasses/binned_theopeps.rds"))
  
  rm(res)
  
  ## MGFs
  message("Splitting MGFs.")
  
  mgf_frames <- readMGF(filepath = mgf_path, min_mass = min_mass, 
                        topn_ms2ions = topn_ms2ions, ret_range = c(0, Inf), 
                        ppm_ms1 = ppm_ms1, 
                        out_path = file.path(mgf_path, "mgf_queries.rds")) %>% 
    dplyr::group_by(frame) %>% 
    dplyr::group_split() %>% 
    setNames(purrr::map_dbl(., ~ .x$frame[1]))

  n_cores <- parallel::detectCores()
  if (n_cores <= 1 || length(mgf_frames) < n_cores) n_cores <- 1

  ## Ion searches
  out <- pmatch_bymgfs(mgf_frames = mgf_frames, 
                       aa_masses_all = aa_masses_all, 
                       n_cores = n_cores, 
                       out_path = out_path, 
                       mod_indexes = mod_indexes, 
                       type_ms2ions = type_ms2ions, 
                       maxn_vmods_per_pep = maxn_vmods_per_pep, 
                       maxn_sites_per_vmod = maxn_sites_per_vmod, 
                       maxn_vmods_sitescombi_per_pep = 
                         maxn_vmods_sitescombi_per_pep, 
                       minn_ms2 = minn_ms2, 
                       ppm_ms1 = ppm_ms1, 
                       ppm_ms2 = ppm_ms2, 
                       digits = digits) %T>%
    saveRDS(file.path(out_path, "ion_matches.rds")) 
  
  ## Peptide scores
  message("Calculating peptide scores.")
  
  out <- out %>% 
    calc_pepscores(topn_ms2ions, type_ms2ions, ppm_ms2, out_path, digits)

  invisible(out)
}


#' Subsets the frames of theoretical peptides.
#' 
#' @inheritParams search_mgf_frames_d
subset_theoframes <- function (mgf_frames, theopeps) {
  if (purrr::is_empty(mgf_frames) || purrr::is_empty(theopeps)) {
    return(NULL)
  }
  
  frames <- as.integer(names(mgf_frames))
  breaks <- which(diff(frames) != 1) + 1
  grps <- findInterval(frames, frames[breaks])
  frames <- split(frames, grps)
  
  frames <- frames %>% 
    purrr::map(~ c(.x[1]-1, .x, .x[length(.x)]+1)) %>% 
    unlist(use.names = FALSE) %>% 
    .[!duplicated(.)]
  
  theopeps[as.character(frames)]
}


#' Matches theoretical peptides (parallel by mgf chunks).
#' 
#' @param aa_masses_all A list of amino acid lookups for all the combination of
#'   fixed and variable modifications.
#' @param n_cores Integer; the number of computer cores.
#' @param mod_indexes Integer; the indexes of fixed and/or variable
#'   modifications.
#' @inheritParams matchMS
#' @inheritParams search_mgf_frames_d
#' @import parallel
pmatch_bymgfs <- function (mgf_frames, aa_masses_all, n_cores, out_path, 
                           mod_indexes, type_ms2ions, maxn_vmods_per_pep, 
                           maxn_sites_per_vmod, 
                           maxn_vmods_sitescombi_per_pep, 
                           minn_ms2, ppm_ms1, ppm_ms2, 
                           digits) {

  mgf_frames <- local({
    labs <- levels(cut(1:length(mgf_frames), n_cores^2))
    
    x <- cbind(
      lower = floor(as.numeric( sub("\\((.+),.*", "\\1", labs))),
      upper = ceiling(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))))
    
    grps <- findInterval(1:length(mgf_frames), x[, 1])
    mgf_frames <- split(mgf_frames, grps)
  })
  
  out <- vector("list", length(aa_masses_all))
  
  for (i in seq_along(out)) {
    # i = 4
    aa_masses <- aa_masses_all[[i]]
    
    # nm_fmods <- attributes(aa_masses)$fmods
    # nm_vmods <- attributes(aa_masses)$vmods
    # nm_nls <- attributes(aa_masses)$vmods_neuloss
    nm_fmods <- attr(aa_masses, "fmods", exact = TRUE)
    nm_vmods <- attr(aa_masses, "vmods", exact = TRUE)
    nm_nls <- attr(aa_masses, "vmods_neuloss", exact = TRUE)
    
    message("Matching against: ", 
            paste0(nm_fmods, 
                   nm_vmods %>% { if (nchar(.) > 0) paste0(" + ", .) else . }, 
                   nm_nls %>% { if (nchar(.) > 0) paste0(" + ", .) else . }))

    theopeps <- readRDS(file.path(out_path, "pepmasses/", 
                                  paste0("binned_theopeps_", i, ".rds")))
    
    # ---
    # temporarily remove `prot_acc` and duplicated `pep_seq` 
    # (frame numbers remain available for quick rejoin later)
    # frame number may be off by +/-1 in three-frame searches
    if (ppm_ms1 < 1000) {
      theopeps <- theopeps %>% 
        purrr::map(~ {
          rows <- !duplicated(.x$pep_seq)
          .x[rows, c("pep_seq", "mass")]
        })
    }
    
    # (1) for a given aa_masses_all[[i]], some mgf_frames[[i]] 
    #   may not be found in theopeps[[i]] 
    mgf_frames <- purrr::map(mgf_frames, ~ {
      x <- .x
      
      oks <- names(x) %in% names(theopeps)
      x <- x[oks]
      
      empties <- purrr::map_lgl(x, purrr::is_empty)
      x[!empties]
    })
    
    # (2) splits `theopeps` in accordance to `mgf_frames` with 
    #   preceding and following frames: (o)|range of mgf_frames[[1]]|(o)
    theopeps <- local({
      frames <- purrr::map(mgf_frames, ~ as.integer(names(.x)))
      
      mins <- purrr::map_dbl(frames, ~ {
        if (length(.x) == 0) x <- 0 else x <- min(.x, na.rm = TRUE)
      })
      
      maxs <- purrr::map_dbl(frames, ~ {
        if (length(.x) == 0) x <- 0 else x <- max(.x, na.rm = TRUE)
      })
      
      nms <- as.integer(names(theopeps))
      
      theopeps <- purrr::map2(mins, maxs, ~ {
        theopeps[which(nms >= (.x - 1) & nms <= (.y + 1))]
      })
    })
    
    # (3) removes unused frames of `theopeps`
    theopeps <- purrr::map2(mgf_frames, theopeps, subset_theoframes)
    
    # (4) removes empties (zero overlap between mgf_frames and theopeps)
    oks <- purrr::map_lgl(mgf_frames, ~ !purrr::is_empty(.x)) | 
      purrr::map_lgl(theopeps, ~ !purrr::is_empty(.x))
    
    mgf_frames <- mgf_frames[oks]
    theopeps <- theopeps[oks]
    
    rm(oks)
    
    cl <- makeCluster(getOption("cl.cores", n_cores))
    
    clusterExport(cl, list("%>%"), 
                  envir = environment(magrittr::`%>%`))
    clusterExport(cl, list("search_mgf_frames_d"), 
                  envir = environment(proteoQ:::search_mgf_frames_d))
    clusterExport(cl, list("search_mgf_frames"), 
                  envir = environment(proteoQ:::search_mgf_frames))
    clusterExport(cl, list("search_mgf"), 
                  envir = environment(proteoQ:::search_mgf))
    clusterExport(cl, list("find_ppm_outer_bypep"), 
                  envir = environment(proteoQ:::find_ppm_outer_bypep))
    clusterExport(cl, list("find_ppm_outer_bycombi"), 
                  envir = environment(proteoQ:::find_ppm_outer_bycombi))

    out[[i]] <- clusterMap(cl, search_mgf_frames_d, 
                           mgf_frames, theopeps, 
                           MoreArgs = list(aa_masses = aa_masses, 
                                           mod_indexes = mod_indexes, 
                                           type_ms2ions = type_ms2ions, 
                                           maxn_vmods_per_pep = 
                                             maxn_vmods_per_pep, 
                                           maxn_sites_per_vmod = 
                                             maxn_sites_per_vmod, 
                                           maxn_vmods_sitescombi_per_pep = 
                                             maxn_vmods_sitescombi_per_pep, 
                                           minn_ms2 = minn_ms2, 
                                           ppm_ms1 = ppm_ms1, 
                                           ppm_ms2 = ppm_ms2, 
                                           digits = digits), 
                           .scheduling = "dynamic") %>% 
      dplyr::bind_rows() %>% # across nodes
      dplyr::mutate(pep_fmod = nm_fmods, 
                    pep_vmod = nm_vmods, 
                    pep_nl = nm_nls)

    # saveRDS(out[[i]], file.path(out_path, paste0("ion_matches_", i, ".rds")))
    
    stopCluster(cl)
    
    gc()
  }
  
  invisible(out)
}


#' Helper of \link{search_mgf_frames}
#'
#' Searches MGFs in a frame at a given combination of fixed and variable
#' modifications.
#'
#' @param theopeps Binned theoretical peptides at a given combination of fixed
#'   and variable.
#' @param aa_masses Amino-acid lookup at a given combination of fixed and
#'   variable.
#' @inheritParams mcalc_monopep
#' @inheritParams search_mgf_frames
#' @export
search_mgf_frames_d <- function (mgf_frames, theopeps, aa_masses, 
                                 mod_indexes, type_ms2ions = "by", 
                                 maxn_vmods_per_pep = 5, 
                                 maxn_sites_per_vmod = 3, 
                                 maxn_vmods_sitescombi_per_pep = 32, 
                                 minn_ms2 = 7, 
                                 ppm_ms1 = 20, ppm_ms2 = 25, digits = 5) {
  # `res[[i]]` contains results for multiple mgfs within a frame
  # (the number of entries equals to the number of mgf frames)
  res <- search_mgf_frames(mgf_frames = mgf_frames, 
                           theopeps = theopeps, 
                           aa_masses = aa_masses, 
                           mod_indexes = mod_indexes, 
                           type_ms2ions = type_ms2ions, 
                           maxn_vmods_per_pep = maxn_vmods_per_pep, 
                           maxn_sites_per_vmod = maxn_sites_per_vmod, 
                           maxn_vmods_sitescombi_per_pep = 
                             maxn_vmods_sitescombi_per_pep, 
                           minn_ms2 = minn_ms2, 
                           ppm_ms1 = ppm_ms1, ppm_ms2 = ppm_ms2, 
                           digits = digits) 

  # flatten mgfs within each frame
  # (the number of entries equals to the number of mgfs)
  res <- res %>% unlist(recursive = FALSE)
  
  empties <- purrr::map_lgl(res, purrr::is_empty)
  
  # !!!dplyr::bind_rows() temporarily not working for list-columns!!!
  res <- do.call(rbind, mgf_frames) %>% 
    dplyr::mutate(matches = res) 
  
  # res <- mgf_frames %>% 
  #   dplyr::bind_rows() %>% 
  #   dplyr::mutate(matches = res) 
  
  res <- res[!empties, ]
  
  rm(mgf_frames, theopeps)
  
  invisible(res)
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
#' @inheritParams mcalc_monopep
#' @inheritParams calc_ms2ionseries
#' @return Matches to each MGF as a list elements. The length of the output is
#'   equal to the number of MGFs in the given frame.
#' @export
search_mgf_frames <- function (mgf_frames, theopeps, aa_masses, mod_indexes, 
                               type_ms2ions = "by", 
                               maxn_vmods_per_pep = 5, 
                               maxn_sites_per_vmod = 3, 
                               maxn_vmods_sitescombi_per_pep = 32, 
                               minn_ms2 = 7, 
                               ppm_ms1 = 20, ppm_ms2 = 25, digits = 5) {
                                
  len <- length(mgf_frames)
  out <- vector("list", len) 
  
  ## --- initiation ---
  mgfs_cr <- mgf_frames[[1]]
  frame <- mgfs_cr[["frame"]][[1]]
  
  theos_bf_ms1 <- theopeps[[as.character(frame-1)]]
  theos_cr_ms1 <- theopeps[[as.character(frame)]]
  
  theomasses_bf_ms1 <- theos_bf_ms1$mass
  theomasses_cr_ms1 <- theos_cr_ms1$mass
  
  theos_bf_ms2 <- purrr::map2(theos_bf_ms1$pep_seq, theomasses_bf_ms1, 
                              calc_ms2ionseries, 
                              aa_masses = aa_masses, 
                              mod_indexes = mod_indexes, 
                              type_ms2ions = type_ms2ions, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              maxn_vmods_sitescombi_per_pep = 
                                maxn_vmods_sitescombi_per_pep, 
                              digits = digits) %>% 
    `names<-`(names(theomasses_bf_ms1))
  
  theos_cr_ms2 <- purrr::map2(theos_cr_ms1$pep_seq, theomasses_cr_ms1, 
                              calc_ms2ionseries, 
                              aa_masses = aa_masses, 
                              mod_indexes = mod_indexes, 
                              type_ms2ions = type_ms2ions, 
                              maxn_vmods_per_pep = maxn_vmods_per_pep, 
                              maxn_sites_per_vmod = maxn_sites_per_vmod, 
                              maxn_vmods_sitescombi_per_pep = 
                                maxn_vmods_sitescombi_per_pep, 
                              digits = digits) %>% 
    `names<-`(names(theomasses_cr_ms1))
  
  ## --- iteration ---
  for (i in seq_len(len)) {
    exptmasses_ms1 <- mgfs_cr[["ms1_mass"]]
    exptmoverzs_ms2 <- mgfs_cr[["ms2_moverz"]]
    
    theos_af_ms1 <- theopeps[[as.character(frame+1)]]
    theomasses_af_ms1 <- theos_af_ms1$mass
    
    theos_af_ms2 <- purrr::map2(theos_af_ms1$pep_seq, theomasses_af_ms1, 
                                calc_ms2ionseries, 
                                aa_masses = aa_masses, 
                                mod_indexes = mod_indexes, 
                                type_ms2ions = type_ms2ions, 
                                maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                maxn_vmods_sitescombi_per_pep = 
                                  maxn_vmods_sitescombi_per_pep, 
                                digits = digits) %>% 
      `names<-`(names(theomasses_af_ms1))
    
    # each `out` for the results of multiple mgfs in one frame
    
    # Browse[4]> exptmasses_ms1
    # [[1]]
    # [1] 748.426367
    
    # [[2]]
    # [1] 748.427407
    
    # Browse[4]> out[[i]]
    # [[1]]
    # named list()
    
    # [[2]]
    # named list()
    
    out[[i]] <- purrr::map2(exptmasses_ms1, exptmoverzs_ms2, 
                            search_mgf, 
                            theomasses_bf_ms1, 
                            theomasses_cr_ms1, 
                            theomasses_af_ms1, 
                            theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                            minn_ms2, ppm_ms1, ppm_ms2) 
    
    # advance to the next frame
    if (i == len) {
      break
    }
    
    mgfs_cr <- mgf_frames[[i+1]]
    new_frame <- mgfs_cr[["frame"]][[1]]
    
    if (isTRUE(new_frame == (frame+1))) {
      theos_bf_ms1 <- theos_cr_ms1
      theos_cr_ms1 <- theos_af_ms1
      
      theomasses_bf_ms1 <- theomasses_cr_ms1
      theomasses_cr_ms1 <- theomasses_af_ms1
      
      theos_bf_ms2 <- theos_cr_ms2
      theos_cr_ms2 <- theos_af_ms2
    } else if (isTRUE(new_frame == (frame+2))) {
      theos_bf_ms1 <- theos_af_ms1
      theos_cr_ms1 <- theopeps[[as.character(new_frame)]]
      
      theomasses_bf_ms1 <- theomasses_af_ms1
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_bf_ms2 <- theos_af_ms2
      
      theos_cr_ms2 <- purrr::map2(theos_cr_ms1$pep_seq, theomasses_cr_ms1, 
                                  calc_ms2ionseries, 
                                  aa_masses = aa_masses, 
                                  mod_indexes = mod_indexes, 
                                  type_ms2ions = type_ms2ions, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  maxn_vmods_sitescombi_per_pep = 
                                    maxn_vmods_sitescombi_per_pep, 
                                  digits = digits) %>% 
        `names<-`(names(theomasses_cr_ms1))
    } else {
      theos_bf_ms1 <- theopeps[[as.character(new_frame-1)]]
      theos_cr_ms1 <- theopeps[[as.character(new_frame)]]
      
      theomasses_bf_ms1 <- theos_bf_ms1$mass
      theomasses_cr_ms1 <- theos_cr_ms1$mass
      
      theos_bf_ms2 <- purrr::map2(theos_bf_ms1$pep_seq, theomasses_bf_ms1, 
                                  calc_ms2ionseries, 
                                  aa_masses = aa_masses, 
                                  mod_indexes = mod_indexes, 
                                  type_ms2ions = type_ms2ions, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  maxn_vmods_sitescombi_per_pep = 
                                    maxn_vmods_sitescombi_per_pep, 
                                  digits = digits) %>% 
        `names<-`(names(theomasses_bf_ms1))
      
      theos_cr_ms2 <- purrr::map2(theos_cr_ms1$pep_seq, theomasses_cr_ms1, 
                                  calc_ms2ionseries, 
                                  aa_masses = aa_masses, 
                                  mod_indexes = mod_indexes, 
                                  type_ms2ions = type_ms2ions, 
                                  maxn_vmods_per_pep = maxn_vmods_per_pep, 
                                  maxn_sites_per_vmod = maxn_sites_per_vmod, 
                                  maxn_vmods_sitescombi_per_pep = 
                                    maxn_vmods_sitescombi_per_pep, 
                                  digits = digits) %>% 
        `names<-`(names(theomasses_cr_ms1))
    }
    
    frame <- new_frame
  }
  
  rm(
    mgf_frames, theopeps, 
    theos_bf_ms1, theos_cr_ms1, theos_af_ms1, 
    theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
    theos_bf_ms2, theos_cr_ms2, theos_af_ms2,
    exptmasses_ms1, exptmoverzs_ms2, 
    mgfs_cr, new_frame, frame
  )
  
  invisible(out)
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
#' @import purrr
#' @return Lists of tibbles.
#' @export 
search_mgf <- function (expt_mass_ms1, expt_moverz_ms2, 
                        theomasses_bf_ms1, theomasses_cr_ms1, theomasses_af_ms1, 
                        theos_bf_ms2, theos_cr_ms2, theos_af_ms2, 
                        minn_ms2 = 7, ppm_ms1 = 20, ppm_ms2 = 25) {

  # --- subsets from the `before` and the `after` by MS1 mass tolerance 
  mass_ranges <- find_mass_error_range(expt_mass_ms1, ppm_ms1)
  bf_allowed <- which(theomasses_bf_ms1 >= mass_ranges[1])
  af_allowed <- which(theomasses_af_ms1 <= mass_ranges[2])
  
  # not used but kept for tidiness
  theomasses_bf_ms1 <- theomasses_bf_ms1[bf_allowed]
  theomasses_af_ms1 <- theomasses_af_ms1[af_allowed]

  theos_bf_ms2 <- theos_bf_ms2[bf_allowed]
  theos_af_ms2 <- theos_af_ms2[af_allowed]

  # --- find MS2 matches ---
  if (is_empty(theos_bf_ms2)) {
    x_bf <- theos_bf_ms2
  } else {
    x_bf <- purrr::map(theos_bf_ms2, find_ppm_outer_bypep, 
                       expt_moverz_ms2, ppm_ms2)
  }
  
  if (is_empty(theos_cr_ms2)) {
    x_cr <- theos_cr_ms2
  } else {
    x_cr <- purrr::map(theos_cr_ms2, find_ppm_outer_bypep, 
                       expt_moverz_ms2, ppm_ms2)
  }
  
  if (is_empty(theos_af_ms2)) {
    x_af <- theos_af_ms2
  } else {
    x_af <- purrr::map(theos_af_ms2, find_ppm_outer_bypep, 
                       expt_moverz_ms2, ppm_ms2)
  }
  
  x <- c(x_bf, x_cr, x_af)
  
  # cleans up
  rows <- map(x, ~ {
    this <- .x
    map_lgl (this, ~ sum(!is.na(.x[["expt"]])) >= minn_ms2)
  })
  x <- map2(x, rows, ~ .x[.y])
  
  empties <- map_lgl(x, is_empty)
  x <- x[!empties]
  
  # length(x) == N(theos_peps) within the ppm window
  # 
  # ATIPIFFDMMLCEYQR
  # (1) ATIPIFFDMMLCEYQR$`0000000050000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #     1  173.  173.
  #     2  175.  175.
  # (2) ATIPIFFDMMLCEYQR$`0000000005000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #     1  173.  173.
  #     2  175.  175.
  # $KADEQMESMTYSTER
  # ...
  
  ## No evidence of M
  # 
  # $ATIPIFFDMMLCEYQR
  # $ATIPIFFDMMLCEYQR$`0000000050000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #   1  173.  173.
  #   2  175.  175.
  #   3  643.  643.
  #   4  790.  790.
  #   5  868.  868.
  #   6 1297. 1297.
  # 
  # $ATIPIFFDMMLCEYQR$`0000000005000000`
  #   A tibble: 6 x 2
  #   theo  expt
  #   <dbl> <dbl>
  #   1  173.  173.
  #   2  175.  175.
  #   3  643.  643.
  #   4  790.  790.
  #   5  868.  868.
  #   6 1297. 1297.
  
  invisible(x)
}


#' Helper: finds the the outer products for vectors of MS2 ions.
#'
#' The same theoretical peptides at different position permutations.
#'
#' @param expts Numeric vector; one series of experimental MS2s.
#' @param theos Numeric vector; one to multiple series of theoretical MS2s.
#' @importFrom purrr map
#' @inheritParams search_mgf_frames
find_ppm_outer_bypep <- function (theos, expts, ppm_ms2) {
  if (!is.list(theos)) {
    theos <- list(theos)
  }
  
  map(theos, find_ppm_outer_bycombi, expts, ppm_ms2)
}


#' Finds the the outer products for a vector of MS2 ions at a given ion series.
#'
#' A theoretical peptide at a given position permutation.
#'
#' @param expts Numeric vector; one series experimental MS2s.
#' @param theos Numeric vector; one series of theoretical MS2s.
#' @importFrom dplyr bind_cols
#' @inheritParams find_ppm_outer_bypep
find_ppm_outer_bycombi <- function (theos, expts, ppm_ms2) {
  d <- outer(theos, expts, "find_ppm_error")
  row_cols <- which(abs(d) <= ppm_ms2, arr.ind = TRUE)
  
  e1 <- expts[row_cols[, 2]] %>% 
    `names<-`(row_cols[, 1])
  
  len <- length(theos)
  
  es <- rep(NA, len) %>% 
    `names<-`(seq_len(len))
  es[names(e1)] <- e1
  
  # the first half are b-ions and the second half are y-ions
  bind_cols(theo = theos, expt = es)
}


