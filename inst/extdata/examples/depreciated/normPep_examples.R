# ===================================
# More in peptide normalization
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

# renormalization against selected samples
# see README at https://github.com/qzhang503/proteoQ for details about 
#   the selection of sample columns and data rows
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  
  # selected samples
  col_refit = Select_sub,
)

# renormalization based on selected rows
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  
  # selected rows 
  slice_at = exprs(prot_n_psm >= 10, pep_n_psm >= 3), 
)

# renormalization against selected sample columns and data rows
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
	
  col_refit = Select_sub,
  slice_at = exprs(prot_n_psm >= 10, pep_n_psm >= 3), 
)

# addtive renormalization against selected sample columns and data rows
# first
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
	
  col_refit = W2,
  slice_at = exprs(prot_n_psm >= 5, pep_n_psm >= 3), 
)

# then
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
	
  col_refit = W16,
  slice_at = exprs(prot_n_psm >= 8, pep_n_psm >= 4), 
)

# ===================================
# Mixed-bed (1): begin with MGKernel
# ===================================
# (1-1) aligned by `MGKernel`
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(20, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 883, 
  maxit = 200, 
  epsilon = 1e-05, 
  filter_peps_by = exprs(pep_len <= 50),
)

pepHist(scale_log2r = TRUE)
pepHist(scale_log2r = FALSE)

# (1-1) followed by median centering for samples specified by `col_refit`
normPep(
  range_log2r = c(5, 95),
  range_int = c(5, 95),  
  
  n_comp = 3,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  col_refit = Select_sub, 

  filter_peps_by = exprs(pep_len <= 50),
)

# `MGKernel_params_N.txt` available for side-effects
pepHist(scale_log2r = TRUE, filename = "sel_mc_z.png",)
pepHist(scale_log2r = FALSE, filename = "sel_mc_n.png",)

# (1-3) aligned by median-centering for all samples
normPep()

# `MGKernel_params_N.txt` available for side-effects
pepHist(scale_log2r = TRUE, filename = "mc_z.png",)
pepHist(scale_log2r = FALSE, filename = "mc_n.png",)

# (1-4) followed by `MGKernel` for samples specified by `col_refit`
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(20, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 883, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_refit = Select_sub, 
  filter_peps_by = exprs(pep_len <= 50),
)

# updated `MGKernel_params_N.txt` for selected samples
pepHist(scale_log2r = TRUE, filename = "sel_mG_z.png",)
pepHist(scale_log2r = FALSE, filename = "sel_mG_n.png",)

# (1-5) fresh start with a different `n_comp = 2`
# will overrule `col_refit` to all samples
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(20, 95), 
  range_int = c(5, 95), 
  n_comp = 2, 
  seed = 883, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_refit = Select_sub, 
  filter_peps_by = exprs(pep_len <= 50),
)

# updated `MGKernel_params_N.txt` for all samples
pepHist(scale_log2r = TRUE, filename = "all_mG2_z.png",)
pepHist(scale_log2r = FALSE, filename = "all_mG2_n.png",)


# ===================================
# Mixed-bed (2): begin with MC
# ===================================
unlink(file.path(dat_dir, "Peptide\\Histogram"), recursive  = TRUE)
unlink(file.path(dat_dir, "Peptide\\Peptide.txt"))

# (2-1) median centering
normPep()

# `MGKernel_params_N.txt` not yet available
pepHist(scale_log2r = TRUE, filename = "mc_z.png",)
pepHist(scale_log2r = FALSE, filename = "mc_n.png",)

# (2-2) this is the first `MGKernel`: 
#   `MGKernel_params_N.txt` not yet available
#   so will overrule `col_refit` to all samples
normPep(
  method_psm_pep = median, 
  method_align = MGKernel, 
  range_log2r = c(20, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 883, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_refit = Select_sub, 
  filter_peps_by = exprs(pep_len <= 50),
)

# `MGKernel_params_N.txt` available
pepHist(scale_log2r = TRUE, filename = "all_mG_z.png",)
pepHist(scale_log2r = FALSE, filename = "all_mG_n.png",)




# ===================================
# Mixed-bed (3): begin with MC
# ===================================
# following Mixed-bed (2) examples
# normalized against housekeeping proteins for selected samples
normPep(
  range_log2r = c(20, 95), 
  range_int = c(5, 95), 
  col_refit = Select_sub,
  slice_hskp = exprs(gene %in% c("ACTB", "GAPDH")),
)

# `MGKernel_params_N.txt` available for side effects
pepHist(scale_log2r = TRUE, filename = "sel_hskp_z.png",)
pepHist(scale_log2r = FALSE, filename = "sel_hskp_n.png",)
