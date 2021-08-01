#' Helper in loading MGFs.
#' 
#' @param min_mass A minimum mass of precursors.
#' @param min_ms2mass A minimum m/z of MS2 ions.
#' @param index_ms2 Logical; if TRUE, converts MS2 m/z values to indices.
#' @inheritParams matchMS
load_mgfs <- function (mgf_path, min_mass = 500L, min_ms2mass = 110L, 
                       topn_ms2ions = 100L, ppm_ms1 = 20L, ppm_ms2 = 25L, 
                       index_ms2 = FALSE) {
  
  rds <- file.path(mgf_path, "mgf_queries.rds")
  
  if (file.exists(rds)) {
    message("Loading cached MGFs: `", rds, "`.")
    mgf_frames <- readRDS(rds)
  } else {
    message("Processing raw MGFs.")
    mgf_frames <- readMGF(filepath = mgf_path, 
                          min_mass = min_mass, 
                          min_ms2mass = min_ms2mass, 
                          topn_ms2ions = topn_ms2ions, 
                          ret_range = c(0, Inf), 
                          ppm_ms1 = ppm_ms1, 
                          ppm_ms2 = ppm_ms2, 
                          index_ms2 = index_ms2, 
                          out_path = rds) 
  }
}


#' Reads MGF files in chunks.
#'
#' @param filepath The file path to a list of mgf files.
#' @param min_mass Numeric; the minimum mass of MS1 species. The value needs to
#'   match the one in  \link{binTheoPeps}.
#' @param topn_ms2ions A non-negative integer; the top-n species for uses in
#'   MS2 ion searches. The default is to use the top-100 ions in an MS2 event.
#' @param ret_range The range of retention time in seconds.
#' @inheritParams matchMS
#' @import stringi
#' @examples
#' \donttest{
#' mgf_queries <- readMGF()
#' }
#' @export
readMGF <- function (filepath = "~/proteoQ/mgfs", 
                     min_mass = 500L, min_ms2mass = 110L, topn_ms2ions = 100L, 
                     ret_range = c(0, Inf), ppm_ms1 = 20L, ppm_ms2 = 25L, 
                     index_ms2 = FALSE, 
                     out_path = file.path(filepath, "mgf_queries.rds")) {
  
  f <- function(x, pos) {
    nm <- file.path(filepath, "temp", paste0("chunk", "_", pos, ".mgf"))
    writeLines(x, nm)
  }
  
  # ---
  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")
  
  if (length(filelist) == 0L) {
    stop("No '.mgf' files under ", filepath, call. = FALSE)
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
                              chunk_size = 5000000L)
    
    out[[i]] <- read_mgf_chunks(filepath = temp_dir, 
                                topn_ms2ions = topn_ms2ions, 
                                ret_range = ret_range, 
                                ppm_ms2 = ppm_ms2, 
                                min_ms2mass = min_ms2mass, 
                                index_ms2 = index_ms2)
    
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
    dplyr::mutate(frame = find_ms1_interval(ms1_mass, 
                                            from = min_mass, 
                                            ppm = ppm_ms1)) %T>% 
    saveRDS(., out_path)
  
  invisible(out)
}


#' Reads mgfs in chunks.
#' 
#' @inheritParams readMGF
#' @import stringi
read_mgf_chunks <- function (filepath = "~/proteoQ/mgfs", 
                             topn_ms2ions = 100L, ret_range = c(0L, Inf), 
                             ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE) {

  filelist <- list.files(path = file.path(filepath), pattern = "^.*\\.mgf$")
  
  if (length(filelist) == 0L) {
    stop("No mgf files under ", filepath, call. = FALSE)
  }
  
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  out <- parallel::clusterApply(cl, file.path(filepath, filelist), 
                                proc_mgf_chunks_i, 
                                topn_ms2ions = topn_ms2ions, 
                                ret_range = ret_range, 
                                ppm_ms2 = ppm_ms2, 
                                min_ms2mass = min_ms2mass, 
                                index_ms2 = index_ms2) %>% 
    dplyr::bind_rows()
  
  parallel::stopCluster(cl)
  
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
      proc_mgfs(gaps, 
                topn_ms2ions = topn_ms2ions, 
                ret_range = ret_range, 
                index_ms2 = index_ms2, 
                ppm_ms2 = ppm_ms2, 
                min_ms2mass = min_ms2mass)
    )
  }
  
  invisible(out)
}


#' Helper of \link{proc_mgf_chunks}.
#' 
#' @param file A mgf chunk (chunk_1.mgf etc.) with prepending file path.
#' @inheritParams read_mgf_chunks
#' @import stringi
proc_mgf_chunks_i <- function (file, topn_ms2ions = 100L, ret_range = c(0L, Inf), 
                               ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE) {

  message("Parsing '", file, "'.")
  
  x <- stri_read_lines(file) %>% 
    proc_mgf_chunks(topn_ms2ions = topn_ms2ions, 
                    ret_range = ret_range, 
                    ppm_ms2 = ppm_ms2, 
                    min_ms2mass = min_ms2mass, 
                    index_ms2 = index_ms2, 
                    filepath = file)
}


#' Processes MGF entries in chunks.
#' 
#' @param lines MGF lines.
#' @inheritParams readMGF
#' @import stringi
proc_mgf_chunks <- function (lines, topn_ms2ions = 100L, ret_range = c(0L, Inf), 
                             ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE, 
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
  
  out <- proc_mgfs(lines, 
                   topn_ms2ions = topn_ms2ions, 
                   ret_range = ret_range, 
                   ppm_ms2 = ppm_ms2, 
                   min_ms2mass = min_ms2mass, 
                   index_ms2 = index_ms2)
}


#' Helper in processing MGF entries in chunks.
#' 
#' @inheritParams proc_mgf_chunks
#' @import stringi
proc_mgfs <- function (lines, topn_ms2ions = 100L, ret_range = c(0L, Inf), 
                       ppm_ms2 = 25L, min_ms2mass = 110L, index_ms2 = FALSE) {
  
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
    rm(list = c("rows", "ms2s"))
  }
  
  # MS1 ions
  ms1s <- lines[begins+2] %>% 
    stri_replace_first_fixed("PEPMASS=", "") %>% 
    purrr::map(stri_split_fixed, " ", n = 2, simplify = TRUE)
  
  ms1_moverzs <- ms1s %>% 
    purrr::map(~ .x[, 1] %>% as.numeric())
  
  ms1_ints <- ms1s %>% 
    purrr::map(~ .x[, 2] %>% as.numeric())
  
  rm(list = c("ms1s"))
  
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
  
  # Rounding
  # ms1_ints <- ms1_ints %>% map(as.integer)
  # ms2_ints <- ms2_ints %>% map(as.integer)
  
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
  
  if (index_ms2) {
    ms2_moverzs <- purrr::map(ms2_moverzs, 
                              find_ms1_interval, 
                              from = min_ms2mass, 
                              ppm = ppm_ms2)
  }
  
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
                       


#' Calculates the frame numbers for a list of experimental MS1 mass by
#' intervals.
#'
#' Needs correct \code{from}.
#' @param mass Numeric; a list of MS1 masses.
#' @inheritParams find_ms1_cutpoints
#' @examples 
#' \donttest{
#' find_ms1_interval(c(500, 800.1))
#' }
#' @return Frame numbers.
#' @seealso find_ms1_cutpoints
find_ms1_interval <- function (mass = 1714.81876, from = 350L, ppm = 20L) {
  d <- ppm/1e6
  ceiling(log(unlist(mass, recursive = FALSE, use.names = FALSE)/from)/log(1+d))
}




