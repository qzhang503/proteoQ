
.onLoad <- function(libname, pkgname) {
  data("prn_annot_hs", "prn_annot_mm", "prn_annot_rn", "prn_annot_dm", 
			"go_sets_hs", "go_sets_mm", "go_sets_rn", 
			"kegg_sets_hs", "kegg_sets_mm", "kegg_sets_rn", 
		package = pkgname, envir = parent.env(environment()))
}

