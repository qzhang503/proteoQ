# ===================================
# Peptide normalization
# ===================================

## !!! run the minimal working example in `?load_expts` before the following extended examples

# ===================================
# (1) `MGKernel`
# ===================================
# (1.1) first-pass normalization: against all sample columns and all data rows  
# (default double trimming by log2FC and intensity percentiles also apply)
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
)

# visualization of histograms containing all data rows
# (for samplicity only for samples under `expt_smry.xlsx::BI_1` and with ratio scaling)
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = bi1.png)

# visualization of the human subset aided by `filter_by` 
# (a little off in relative to the whole data set)
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = bi1_human.png)

# visualization of the mouse subset 
# (a lot off in relative to the whole data set)
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "mouse"), filename = bi1_mouse.png)

# (1.2) additive normalization to (1.1): against `expt_smry.xlsx::W2` samples and based on their human subset
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_refit = W2,  
  slice_peps_by = exprs(species == "human"),
)

# visualization of human subset 
# (W2 samples are now aligned)
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = bi1_human_slicehu.png)

# (1.3) additive normalization to (1.2): against `expt_smry.xlsx::W16` samples and based on their human subset
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_refit = W16,  
  slice_peps_by = exprs(species == "human"),
)

# visualization of human subset 
# (W16 samples are now also aligned)
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = bi1_human_slicehu2.png)

# side effects: to recaptulate the misalignment between human- and all-data rows
# (curves based on the latest `standPep` at `method_align = MGKernel`)
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = bi1_recap.png)


# ===================================
# (2) Mixed-bed
# ===================================
library(proteoQ)
dat_dir <- "~\\proteoQ\\examples"
load_expts()

# start over for demonstration purpose
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

# for simplicity, only visualize samples under `expt_smry.xlsx::BI_1`
# (initial `Peptide.txt` from `mergePep()` always aligned by median centering (MC))
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (2.1) the first `MGKernel` alignment against all samples using all data rows
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 883, 
  maxit = 200, 
  epsilon = 1e-05, 
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG.png)

# (2.2) followed by MC against selected samples using all data rows
standPep(
  method_align = MC, 
  n_comp = 3,
  seed = 883,
  maxit = 200,
  epsilon = 1e-05,
  col_refit = Select_sub, 
)

# the net result is mixed-bed alignment of MGKernel and MC
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mix.png)


# ===================================
# (3) more MC to MGKernel examples
# ===================================
# start over
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (3.1) also MC alignment for all samples but using only selected data rows
# (at `method_align = MC`, `n_comp`, `seed`, `maxit`, `epsilon` will be ignored)
standPep(
  method_align = MC, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_peps_by = exprs(prot_n_psm >= 10),
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc_selrows.png)

# (3.2) change to `MGKernel` for all samples using selected data rows
# (first-pass `MGKernel` will apply to all samples)
standPep(
  method_align = MGKernel, 
  
  # will be forced to all samples since this is the first `MGKernel`
  col_refit = Select_sub, 
  
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_peps_by = exprs(prot_n_psm >= 10),
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG_selrows.png)

# (3.3) back to MC for selected samples using selected data rows
# mixed-bed again: MC for samples indicated by Select_sub and MGKernel for the rest
standPep(
  method_align = MC, 
  n_comp = 3, 
  seed = 883, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_refit = Select_sub, 
  slice_peps_by = exprs(prot_n_psm >= 10),
)

# with the first `MGKernel` in (3.2), density curves are available for 
# the side-effects comparing `MC` and `MGKernel`
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mix_selcols_selrows.png)


# ===================================
# (4) changed `n_comp`
# ===================================
# start over
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (4.1) change to `MGKernel` for all samples
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG3.png)

# (4.2) fresh start with a different `n_comp`
# will ignore `col_refit = Select_sub` and apply to all samples 
standPep(
  method_align = MGKernel, 
  col_refit = Select_sub, 
  n_comp = 2, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG2.png)


# ===================================
# (5) housekeeping normalizers
# ===================================
# start over
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

# (initial `Peptide.txt` from `mergePep()` aligned by MC)
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

    # (5.1) change to `MGKernel`
    standPep(
      method_align = MGKernel, 
      n_comp = 3, 
      seed = 400, 
      maxit = 200, 
      epsilon = 1e-05, 
    )
    
    pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG.png)

    # (5.2) renormalize against data from selected proteins for selected sample(s)
    # (error: too few entries for fitting with multiple Gaussians)
    standPep(
      method_align = MGKernel, 
      col_refit = W2, 
      n_comp = 3, 
      seed = 400, 
      maxit = 200, 
      epsilon = 1e-05, 
      slice_hskp = exprs(gene %in% c("ACTB", "GAPDH")),
    )

# (5.3) try again with `MC` for W2 samples
standPep(
  method_align = MC, 
  # col_refit = W2, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_hskp = exprs(gene %in% c("ACTB", "GAPDH")),
)

# suggest "ACTB" and "GAPDH" less abundant in W16 
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc_hskp.png)
# pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc_selcols_selrowshskp.png)

