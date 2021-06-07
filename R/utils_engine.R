#' Finds a file directory.
#' 
#' With the option of creating a direcotry.
#' 
#' @param path A file path.
#' @param create Logical; if TRUE create the path if not yet existed.
find_dir <- function (path, create = FALSE) {
  stopifnot(length(path) == 1L)
  
  path <- gsub("\\\\", "/", path)
  
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
save_call2 <- function(path, fun, time = NULL) {
  
  stopifnot(length(path) == 1L, length(fun) == 1L)

  call_pars <- mget(names(formals(fun)), envir = caller_env(), inherits = FALSE)
  call_pars[names(call_pars) == "..."] <- NULL
  
  if (is.null(time)) {
    p2 <- create_dir(path)
    save(call_pars, file = file.path(p2, paste0(fun, ".rda")))
  } else {
    stopifnot(length(time) == 1L)
    p2 <- create_dir(file.path(path, fun))
    save(call_pars, file = file.path(p2, paste0(time, ".rda")))
  }
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
find_callarg_vals <- function (time = NULL, path = NULL, fun = NULL, 
                               args = NULL) {

  stopifnot(length(path) == 1L, length(fun) == 1L)

  if (is.null(time)) {
    file <- file.path(path, fun)
  } else {
    stopifnot(length(time) == 1L)
    file <- file.path(path, fun, time)
  }
  
  if (!file.exists(file)) {
    # waring(file, " not found.", call. = FALSE)
    return(NULL)
  }
  
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
#' @param nms Names of arguments to be included in or excluded from matching.
#' @param type Logical; if TRUE, includes all arguments in \code{nms}. At FALSE,
#'   excludes all \code{nms}.
#' @inheritParams find_callarg_vals
#' @importFrom rlang caller_env
#' @return An empty object if no matches.
match_calltime <- function (path = "~/proteoQ/.MSearch/Cache/Calls", 
                            fun = "calc_pepmasses", 
                            nms = c("parallel", "out_path"), 
                            type = c(TRUE, FALSE)) {
  
  stopifnot(length(path) == 1L, length(fun) == 1L)
  
  if (length(type) > 1L) type <- TRUE

  # current
  args <- mget(names(formals(fun)), envir = caller_env(), inherits = FALSE) %>% 
    { if (type) .[names(.) %in% nms] else .[! names(.) %in% nms] }
  
  if (is_empty(args)) stop("Arguments for matching is empty.", call. = FALSE)
  
  args <- args %>% map(sort)
  
  times <- list.files(path = file.path(path, fun), 
                      pattern = "\\.rda$", 
                      all.files = TRUE) 

  # cached
  cached <- map(times, find_callarg_vals, path = path, fun = fun, 
                args = names(args))
  
  cached <- cached %>% map(~ map(.x, sort))
  
  # matched
  oks <- map_lgl(cached, identical, args)
  
  times[oks] %>% gsub("\\.rda$", "", .)
}


#' Matches the value of a function argument.
#' 
#' @param f The function where one of its argument value will be matched.
#' @param arg The argument where its value will be matched.
#' @param val The value to be matched.
#' @param default The default value of \code{arg}.
#' @examples 
#' \donttest{
#' foo <- function (quant = c("none", "tmt6", "tmt10", "tmt11", "tmt16")) {
#'  val <- rlang::enexpr(quant)
#'  val <- match_valexpr(f = !!match.call()[[1]], arg = quant, val = !!val)
#' }
#' 
#' default <- foo()
#' custom <- foo(quant = tmt6)
#' }
match_valexpr <- function (f = NULL, arg = NULL, val = NULL, default = NULL) {
  
  f <- rlang::enexpr(f)
  arg <- rlang::enexpr(arg)
  
  stopifnot(length(f) == 1L, length(arg) == 1L)
  
  f <- rlang::as_string(f)
  arg <- rlang::as_string(arg)
  ok_vals <- as.character(formals(f)[[arg]])[-1]
  
  # may be plural, not yet to string
  val <- rlang::enexpr(val)
  
  if (is_empty(val)) {
    stop("`val` cannot be empty.", call. = FALSE)
  }
  
  if (length(val) > 1) {
    if (is.null(default)) {
      if (is.call(val[1])) {
        val <- val[[2]]
      } else {
        val <- val[[1]]
      }
    } else {
      stopifnot(length(default) == 1L)
      val <- default
    }
  } else {
    val <- rlang::as_string(val)
  }
  
  if (! val %in% ok_vals) {
    stop("`", val, "` is not an option for `", f, "`.", call. = FALSE)
  }
  
  val
}


#' Deletes files under a directory.
#' 
#' The directory will be kept. 
#' @param path A file path.
delete_files <- function (path, ...) {
  dots <- rlang::enexprs(...)
  recursive <- dots[["recursive"]]
  
  if (is.null(recursive)) recursive <- TRUE
  
  stopifnot(is.logical(recursive))
  
  nms <- list.files(path = file.path(path), recursive = recursive, ...)
  
  if (!purrr::is_empty(nms)) {
    suppressMessages(file.remove(file.path(path, nms)))
  }
  
  invisible(NULL)
}
