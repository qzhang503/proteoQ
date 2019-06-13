
.onLoad <- function(libname, pkgname) {
  data(
    "prn_annot_hs", 
    "go_sets_hs", 
    "kegg_sets_hs", 
		package = pkgname, envir = parent.env(environment()))
}
