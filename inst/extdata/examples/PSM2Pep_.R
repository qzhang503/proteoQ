\donttest{
# ===================================
# PSMs to peptides
# ===================================

## !!!require the brief working example in `?load_expts`

# default median statistics
PSM2Pep()

# mean statistics
PSM2Pep(method_psm_pep = mean)  

## Not run: 
# cut-offs in precursor intensity
# (error if `pep_tot_int` are all NAs)
PSM2Pep(filter_ms1int = rlang::exprs(pep_tot_int >= 1E4)) 

# cut-offs in precursor intensity by percentiles
PSM2Pep(filter_ms1int = rlang::exprs(
  pep_tot_int >= quantile(pep_tot_int, probs = .05, na.rm = TRUE), 
  pep_tot_int <= quantile(pep_tot_int, probs = .99, na.rm = TRUE), 
))
## End(Not run)
}
