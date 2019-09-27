
.onLoad <- function(libname, pkgname) {
  data(
    # "go_sets_hs", 
		package = pkgname, envir = parent.env(environment())
	)
}
