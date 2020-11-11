#' A wrapper of \link[stats]{dist} with the handling of partial argument
#' matches.
#' 
my_dist <- function (...) {
  dots <- rlang::enexprs(...)
  
  dummies <- c("p")
  
  msgs <- c(
    "`p` in `dist()` not used."
  )
  
  stopifnot(length(dummies) == length(msgs))
  
  purrr::walk2(dummies, msgs, ~ {
    if (.x %in% names(dots)) {
      warning(.y, call. = FALSE)
      dots[[.x]] <<- NULL
    } 
  })
  
  rlang::expr(stats::dist(!!!dots)) %>% 
    rlang::eval_bare(env = caller_env())
}

