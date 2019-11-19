# ===================================
# Peptide
# ===================================
normPep(
  method_psm_pep = median,
  method_align = MGKernel,
  range_log2r = c(5, 95),
  range_int = c(5, 95),
  n_comp = 3,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  # filter_by_sp = exprs(species == "human"),
)


# ===================================
# Protein
# ===================================
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >=2, pep_isunique),
)

# renormalization against selected samples
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >=2, pep_isunique),
  
  # selected samples
  col_refit = Select_sub,	
)

# renormalization against selected samples
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >=2, pep_isunique),
  
  # selected rows 
  slice_at = exprs(prot_n_psm >= 10), 
)

# renormalization against selected samples using partial data
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >=2, pep_isunique),
  slice_at = exprs(prot_n_psm >= 10), 
  col_refit = Select_sub,
)








