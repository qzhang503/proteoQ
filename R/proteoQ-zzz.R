
.onLoad <- function(libname, pkgname) {
  data(
		package = pkgname, envir = parent.env(environment())
	)
}
