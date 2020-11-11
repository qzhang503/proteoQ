#' Set global \code{dat_dir}
#' 
#' @inheritParams load_expts
set_dat_dir <- function(dat_dir = NULL) {
  if (is.null(dat_dir)) dat_dir <- get_gl_dat_dir()
  
  if (!fs::dir_exists(dat_dir)) {
    new_dat_dir <- fs::path_expand_r(dat_dir)
    new_dat_dir2 <- fs::path_expand(dat_dir)
    
    if (fs::dir_exists(new_dat_dir)) {
      dat_dir <- new_dat_dir
    } else if (fs::dir_exists(new_dat_dir2)) {
      dat_dir <- new_dat_dir2
    } else {
      stop(dat_dir, " not existed.", call. = FALSE)
    }
  }
  
  nm <- names(formals())[1]
  
  assign(nm, dat_dir, envir = .GlobalEnv, inherits = FALSE)

  invisible(dat_dir)
}


#' Fetch global \code{dat_dir}
get_gl_dat_dir <- function() {
  dat_dir <- tryCatch(get("dat_dir", envir = .GlobalEnv, inherits = FALSE), 
                      error = function(e) 1)
  if (dat_dir == 1) {
    stop("Unknown working directory; ", 
         "run `load_expts(\"my/fabulous/working/directory\")` first.", 
         call. = FALSE)
  }
  
  invisible(dat_dir)
}

