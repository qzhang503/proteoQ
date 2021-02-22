#' A wrapper of \link[stats]{dist} with the handling of partial argument
#' matches.
#' 
#' @param ... Arguments for \link[stats]{dist}
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


#' Cosine similarity.
#' 
#' @param M An input Matrix.
#' 
#' Vectors are in rows: mtcars[1, ].
#' 
#' @examples 
#' sim <- cos_sim(as.matrix(mtcars[1:3, ]))
#' 
#' # distances (against normalized vectors)
#' as.dist(1 - sim)
cos_sim <- function (M) {
  stopifnot(is.matrix(M))
  
  L <- sqrt(rowSums(M * M)) # vector lengths; no `mean` subtraction
  Mn <- M / L # normalized M 
  Mn %*% t(Mn) # dot products; vectors in the rows of Mn
}

