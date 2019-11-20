# ===================================
# Protein normalization
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
  filter_by = exprs(prot_n_pep >= 2, pep_isunique == TRUE),
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
  filter_by = exprs(prot_n_pep >= 2, pep_isunique == TRUE),
  
  # selected samples
  col_refit = Select_sub,	
)

# renormalization based on selected rows
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >= 2, pep_isunique == TRUE),
  
  # selected rows 
  slice_at = exprs(prot_n_psm >= 10), 
)

# renormalization against selected sample columns and data rows
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >= 2, pep_isunique),
	
  slice_at = exprs(prot_n_psm >= 10), 
  col_refit = Select_sub,
)

# addtive renormalization against selected sample columns and data rows
# first
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >= 2, pep_isunique),
	
  slice_at = exprs(prot_n_psm >= 10), 
  col_refit = W2,
)

# then
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_by = exprs(prot_n_pep >= 2, pep_isunique),
	
  slice_at = exprs(prot_n_psm >= 5), 
  col_refit = W16,
)







