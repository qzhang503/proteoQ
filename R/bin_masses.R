#' Helper in binning precursor masses
#' 
#' @param res The results from \link{calc_pepmasses2}.
#' @param min_mass A minimum mass of precursors.
#' @param max_mass A maximum mass of precursors.
#' @inheritParams matchMS
bin_ms1masses <- function (res, min_mass = 500L, max_mass = 10000L, 
                           ppm_ms1 = 20L, out_path) {
  
  .path_cache <- get(".path_cache", envir = .GlobalEnv)
  .path_fasta <- get(".path_fasta", envir = .GlobalEnv)
  .time_stamp <- get(".time_stamp", envir = .GlobalEnv)
  
  # Targets
  bins <- list.files(path = file.path(.path_fasta, "pepmasses", .time_stamp), 
                     pattern = "binned_theopeps_\\d+\\.rds$")
  
  if (length(bins) == 0L) {
    message("Binning MS1 masses (theoretical target).")
    
    binTheoPeps(res = res$fwd, min_mass = min_mass, max_mass = max_mass, 
                ppm = ppm_ms1, 
                out_path = file.path(.path_fasta, "pepmasses", .time_stamp, 
                                     "binned_theopeps.rds"))
  } else {
    message("Loading bins of MS1 masses from cache (theoretical target).")
  }
  
  # Decoys
  bins2 <- list.files(path = file.path(.path_fasta, "pepmasses", .time_stamp), 
                      pattern = "binned_theopeps_rev_\\d+\\.rds$")
  
  if (length(bins2) == 0L) {
    message("Binning MS1 masses (theoretical decoy).")
    
    binTheoPeps(res = res$rev, min_mass = min_mass, max_mass = max_mass, 
                ppm = ppm_ms1, 
                out_path = file.path(.path_fasta, "pepmasses", .time_stamp, 
                                     "binned_theopeps_rev.rds"))
  } else {
    message("Loading bins of MS1 masses from cache (theoretical decoy).")
  }
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
binTheoPeps <- function (res, min_mass = 500L, max_mass = 10000L, ppm = 20L, 
                         out_path = file.path("~/proteoQ/outs", 
                                              "pepmasses/binned_theopeps.rds")) {
  
  out_dir <- create_dir(gsub("(^.*/).*$", "\\1", out_path))
  
  out_nms <- gsub("^.*/(.*)\\.[^\\.].*$", "\\1", out_path) %>% 
    paste(seq_along(res), sep = "_") %>% 
    paste0(".rds")
  
  res <- res %>% 
    purrr::map(attributes) %>% 
    purrr::map(`[[`, "data")
  
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(getOption("cl.cores", n_cores))
  
  out <- purrr::map(res, ~ {
    .x <- chunksplit(.x, n_cores)
    
    parallel::clusterApplyLB(cl, .x, bin_theopeps, min_mass, max_mass, ppm) %>% 
      purrr::flatten()
  }) %>% 
    parallel::clusterMap(cl, cbind_theopepes, ., file.path(out_dir, out_nms))
  
  parallel::stopCluster(cl)
  
  rm(list = c("res"))
  gc()
  
  invisible(out)
}


#' Helper: separates theoretical peptides into mass groups.
#' 
#' @param peps A list of theoretical peptides with masses.
#' @inheritParams binTheoPeps
bin_theopeps <- function (peps, min_mass = 500L, max_mass = 10000L, ppm = 20L) {
  
  ps <- find_ms1_cutpoints(min_mass, max_mass, ppm)
  
  out <- purrr::imap(peps, ~ {
    prot_peps <- .x
    
    frames <- findInterval(prot_peps, ps)
    
    list(pep_seq = names(prot_peps), 
         mass = prot_peps, 
         frames = frames, 
         prot_acc = .y)
  })
}


#' Finds the cut-points of MS1 masses for binning.
#'
#' The cut-points will be used as the boundary in data binning. Note that the
#' upper bound is open.
#'
#' @param from Numeric; the starting MS1 mass.
#' @param to Numeric; the ending MS1 mass.
#' @param ppm Numeric; the ppm for data binning.
#' @return Cut points.
#' @seealso find_ms1_interval
find_ms1_cutpoints <- function (from = 500L, to = 10000L, ppm = 20L) {
  
  d <- ppm/1e6
  n <- ceiling(log(to/from)/log(1+d))
  
  x <- vector("numeric", n)
  x[1] <- from
  
  for (i in seq_len(n-1)) {
    x[i+1] <- x[i] * (1 + d)
  }
  
  x
}


#' Combines theoretical peptides after binning.
#' 
#' @param out A list of binned theoretical peptides.
#' @param out_nm The output file path and name.
cbind_theopepes <- function (out, out_nm) {
  
  prot_acc <- purrr::imap(out, ~ rep(.y, length(.x$pep_seq))) %>% 
    do.call(`c`, .) %>% 
    unname()
  
  pep_seq <- purrr::map(out, `[[`, "pep_seq") %>% 
    do.call(`c`, .) %>% 
    unname()
  
  mass <- purrr::map(out, `[[`, "mass") %>% 
    do.call(`c`, .) %>% 
    unname()
  
  frame <- purrr::map(out, `[[`, "frames") %>% 
    do.call(`c`, .) %>% 
    unname()
  
  out <- data.frame(pep_seq = pep_seq, 
                    mass = mass, 
                    frame = frame, 
                    prot_acc = prot_acc) %>% 
    dplyr::arrange(frame, pep_seq, prot_acc) %>% 
    split(., .$frame, drop = FALSE) %T>% 
    saveRDS(., out_nm)
  
  invisible(out)
}


