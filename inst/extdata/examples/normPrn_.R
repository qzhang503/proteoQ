\donttest{
# ===================================
# Protein normalization
# ===================================

## !!!require the brief working example in `?load_expts`

# ===================================
# (1) `MGKernel`
# ===================================

# !!! Initial `Protein.txt` results from `Pep2Prn()` are in median centering.
# !!! The first `MGKernel` normalization will always be applied to all samples.
# !!! If changing `n_comp`, succeeding `MGKernel` normalization(s) will be 
#     applied to all samples (a fresh-start principle, see also section 4)

# fresh start of `Protein.txt` (for demonstration)
unlink(file.path(dat_dir, "Protein"), recursive = TRUE, force = TRUE)
Pep2Prn(use_unique_pep = TRUE)

# data in initial `Protein.txt` from `Pep2Prn()` are aligned by MC
# (for samplicity only plot against samples indicated under `expt_smry.xlsx::BI_1`)
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (1.1) the first `MGKernel` after `Pep2Prn()`
# (default double trimming by log2FC and intensity percentiles also apply)
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
)

# histograms for all data rows
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mG.png)

# human subset
# (a little off in relative to all data rows)
prnHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = mG_human.png)

# mouse subset 
# (a lot off in relative to all data rows)
prnHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "mouse"), filename = mG_mouse.png)

# (1.2) additive step to (1.1): against `expt_smry.xlsx::W2` samples and based on their human subset
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_select = W2,  
  slice_peps_by = exprs(species == "human"),
)

# human subset 
# (W2 samples are now aligned)
prnHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = mG_bi1_human_slicehuw2.png)

# (1.3) additive to (1.2): against `expt_smry.xlsx::W16` samples and based on their human subset
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  col_select = W16,  
  slice_peps_by = exprs(species == "human"),
)

# human subset 
# (W16 samples are now also aligned)
prnHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = mG_human_slicehuw16.png)

# side effects: to recaptulate the misalignment between human data and human + mouse data
# (this is because density curves are based on the latest `standPrn` at `method_align = MGKernel`)
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mG_recap.png)


# ===================================
# (2) Mixed-bed
# ===================================
# start over
unlink(file.path(dat_dir, "Protein"), recursive = TRUE, force = TRUE)
Pep2Prn(use_unique_pep = TRUE)

# data in initial `Protein.txt` from `Pep2Prn()` are aligned by MC
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (2.1) the first `MGKernel` alignment against all samples using all data rows
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 883, 
  maxit = 200, 
  epsilon = 1e-05, 
)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mG.png)

# (2.2) followed by MC against selected samples using all data rows
# (arguments `n_comp`, `seed`, `maxit`, `epsilon` have no effects at `method_align = MC`)
standPrn(
  method_align = MC,
  
  n_comp = 3,
  seed = 883,
  maxit = 200,
  epsilon = 1e-05,
  
  col_select = Select_sub, 
)

# the net result is mixed-bed alignment of MGKernel and MC
# (with samples indicated by `Select_sub` aligned by MC and the remaining by MGKernel)
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mix.png)


# ===================================
# (3) Mixed-bed against data subset
# ===================================
# start over
unlink(file.path(dat_dir, "Protein"), recursive = TRUE, force = TRUE)
Pep2Prn(use_unique_pep = TRUE)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (3.1) MC alignment for all samples, but only selected data rows used for normalization
standPrn(
  method_align = MC, 
  slice_peps_by = exprs(prot_n_psm >= 10),
)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mc_selrows.png)

# (3.2) first `MGKernel` for all samples using selected data rows
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05,
  
  # will be forced to all samples since this is the first `MGKernel`
  # col_select = Select_sub, 

  slice_peps_by = exprs(prot_n_psm >= 10),
)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mG_selrows.png)

# (3.3) back to MC for selected samples using selected data rows
# (mixed-bed again, but based on data subset by `prot_n_psm >= 10`)
standPrn(
  method_align = MC, 
  col_select = Select_sub, 
  slice_peps_by = exprs(prot_n_psm >= 10),
)

# side-effects comparing `MC` and `MGKernel`
# (density curves are from the preceding `MGKernel` in (3.2))
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mix_selcols_selrows.png)


# ===================================
# (4) Modified `n_comp`
# ===================================
# start over
unlink(file.path(dat_dir, "Protein"), recursive = TRUE, force = TRUE)
Pep2Prn(use_unique_pep = TRUE)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (4.1) first `MGKernel` at `n_comp = 3` for all samples
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mG3.png)

# (4.2) a fresh start since changing `n_comp`
# (e.g. `col_select = Select_sub` ignored; instead apply `MGKernel` to all samples)
standPrn(
  method_align = MGKernel, 
  n_comp = 2, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05,
  
  # ignored
  col_select = Select_sub, 
)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mG2.png)


# ===================================
# (5) housekeeping normalizers: 
#     suggest `method_align = MC`
# ===================================
# start over
unlink(file.path(dat_dir, "Protein"), recursive = TRUE, force = TRUE)
Pep2Prn(use_unique_pep = TRUE)

# initial `Protein.txt` from `Pep2Prn()` is aligned by MC
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (GAPDH more abundant in W16 according to MC alignment)
prnHist(scale_log2r = TRUE, col_select = BI_1, filter_ = exprs(gene == "GAPDH"), filename = mcGAPDH.png)

# heat map
# (outputs under `Protein/Heatmap` folder; for help, ?pepHM)
prnHM(
  col_select = BI_1, 
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = FALSE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  
  filter_by = exprs(gene %in% c("GAPDH")),
  filename = "mcGAPDH.png",
)

# renormalize against GAPDH
# (not `MGkernel` with few data points)
standPrn(
  method_align = MC, 
  slice_hskp = exprs(gene %in% c("GAPDH")),
)

# (now log2FC profiles aligned by GAPDH)
prnHist(scale_log2r = TRUE, col_select = BI_1, filename = wrong_hskp.png)

# not to keep the above example with no data under `JHU_TMT1` and `PNNL_TMT1`
# (surely need different normalizer(s))
unlink(file.path(dat_dir, "Protein"), recursive = TRUE, force = TRUE)


## Not run: 
# change to `MGKernel`
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
)

prnHist(scale_log2r = TRUE, col_select = BI_1, filename = mG.png)

# then renormalize against data from GAPDH
# (error: too few entries for fitting with multiple Gaussians)
standPrn(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
  
  col_select = W2, 
  slice_hskp = exprs(gene %in% c("GAPDH")),
)
## End(Not run)
}
