#' Finds a file directory.
#' 
#' With the option of creating a direcotry.
#' 
#' @param path A file path.
#' @param create Logical; if TRUE create the path if not yet existed.
find_dir <- function (path, create = FALSE) {
  stopifnot(length(path) == 1L)
  
  # path <- gsub("\\\\", "/", path)
  p1 <- fs::path_expand_r(path)
  p2 <- fs::path_expand(path)
  
  if (fs::dir_exists(p1)) {
    path <- p1
  } else if (fs::dir_exists(p2)) {
    path <- p2
  } else {
    if (create) {
      dir.create(file.path(path), recursive = TRUE, showWarnings = FALSE)
      path <- p1
    } else {
      message(path, " not found.")
      path <- NULL
    }
  }
  
  path
}


#' Creates a file directory.
#' 
#' @inheritParams find_dir
create_dir <- function (path) {
  find_dir(path, create = TRUE)
}


#' Saves the arguments in a function call
#' 
#' @param path A (system) file path.
#' @param fun The name of function being saved.
#' @param time The time stamp.
#' @importFrom rlang caller_env
save_call2 <- function(path, fun, time) {
  
  stopifnot(length(path) == 1L, 
            length(fun) == 1L, 
            length(time) = 1L)
  
  p2 <- create_dir(file.path(path, fun))
  
  call_pars <- mget(names(formals(fun)), envir = caller_env(), inherits = FALSE)
  call_pars[names(call_pars) == "..."] <- NULL

  save(call_pars, file = file.path(p2, paste0(time, ".rda")))
}


#' Finds the setting of a arguments.
#' 
#' Not currently used.
#' 
#' @param arg Argument to be matched.
#' @inheritParams save_call2
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$%
find_callarg_val <- function (time = ".2021-05-21_211227", 
                              path = "~/proteoQ/.MSearch/Cache/Calls", 
                              fun = "calc_pepmasses", arg = "fasta") {

  stopifnot(length(arg) == 1L)
  
  find_callarg_vals(time, path, fun, arg) %>% 
    purrr::flatten()
}


#' Finds the values of a list of arguments.
#'
#' @param args Arguments to be matched.
#' @inheritParams save_call2
#' @import dplyr purrr 
#' @importFrom magrittr %>% %T>% %$%
#' @return The time stamp of a matched cache results.
find_callarg_vals <- function (time = ".2021-05-21_211227.rda", 
                               path = "~/proteoQ/.MSearch/Cache/Calls", 
                               fun = "calc_pepmasses", 
                               args = c("fasta", "max_miss")) {
  
  stopifnot(length(path) == 1L, 
            length(fun) == 1L, 
            length(time) = 1L)
  
  file <- file.path(path, fun, time)
  
  if (!file.exists(file)) stop(file, " not found.", call. = FALSE)
  
  load(file = file)
  
  nots <- which(! args %in% names(call_pars))
  
  if (!purrr::is_empty(nots)) {
    stop("Arguments '", paste(args[nots], collapse = ", "), 
         "' not found in the latest call to ", fun, call. = FALSE)
  }
  
  call_pars[args]
}


#' Finds the time stamp of a matched call from cached results.
#' 
#' @param excludes Arguments excluded from matches.
#' @inheritParams find_callarg_vals
#' @importFrom rlang caller_env
#' @return An empty object if no matches.
match_calltime <- function (path = "~/proteoQ/.MSearch/Cache/Calls", 
                            fun = "calc_pepmasses", 
                            excludes = c("parallel", "out_path")) {
  
  stopifnot(length(path) == 1L, length(fun) == 1L)

  args <- mget(names(formals(fun)), envir = caller_env(), inherits = FALSE) %>% 
    .[! names(.) %in% excludes]

  if (is_empty(args)) stop("`args` cannot be empty.", call. = FALSE)
  
  times <- list.files(path = file.path(path, fun), 
                      pattern = "\\.rda$", 
                      all.files = TRUE) 

  cached <- map(times, find_callarg_vals, path = path, fun = fun, 
                args = names(args))
  
  oks <- map_lgl(cached, identical, args)
  
  times[oks] %>% gsub("\\.rda$", "", .)
}
