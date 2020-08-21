\donttest{
# ===================================
# Peptide normalization
# ===================================

## !!!require the brief working example in `?load_expts`

# ===================================
# (1) `MGKernel`
# ===================================

# !!! Initial `Peptide.txt` results from `mergePep()` are in median centering.
# !!! The first `MGKernel` normalization will always be applied to all samples.
# !!! If changing `n_comp`, succeeding `MGKernel` normalization(s) will be 
#     applied to all samples (a fresh-start principle, see also section 4)

# (1.1) the first `MGKernel` after `mergePep()`
# (default double trimming by log2FC and intensity percentiles also applied)
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
)

# histograms for all data rows
# (for samplicity only plot against samples indicated under `expt_smry.xlsx::BI_1`)
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = bi1.png)

# human subset
# (a little off in relative to all data rows)
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = bi1_human.png)

# mouse subset 
# (a lot off in relative to all data rows)
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "mouse"), filename = bi1_mouse.png)

# (1.2) additive step to (1.1): against `expt_smry.xlsx::W2` samples and based on their human subset
standPep(
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
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = bi1_human_slicehuw2.png)

# (1.3) additive to (1.2): against `expt_smry.xlsx::W16` samples and based on their human subset
standPep(
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
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_by = exprs(species == "human"), filename = bi1_human_slicehuw16.png)

# a side effect: to recapitulate the misalignment between human data and human + mouse data
# (this is because density curves are based on the latest `standPep` at `method_align = MGKernel`)
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = bi1_recap.png)


# ===================================
# (2) Mixed-bed
# ===================================
# start over
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

# data in initial `Peptide.txt` from `mergePep()` are aligned by MC
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
# (arguments `n_comp`, `seed`, `maxit`, `epsilon` have no effects at `method_align = MC`)
standPep(
  method_align = MC,
  
  n_comp = 3,
  seed = 883,
  maxit = 200,
  epsilon = 1e-05,
  
  col_select = Select_sub, 
)

# the net result is mixed-bed alignment of MGKernel and MC
# (with samples indicated by `Select_sub` aligned by MC and the remaining by MGKernel)
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mix.png)


# ===================================
# (3) Mixed-bed against data subset
# ===================================
# start over
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (3.1) MC alignment for all samples, but only selected data rows used for normalization
standPep(
  method_align = MC, 
  slice_peps_by = exprs(prot_n_psm >= 10),
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc_selrows.png)

# (3.2) the first `MGKernel` for all samples using selected data rows
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05,
  
  # will be forced to all samples since this is the first `MGKernel`
  # col_select = Select_sub, 

  slice_peps_by = exprs(prot_n_psm >= 10),
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG_selrows.png)

# (3.3) back to MC for selected samples using selected data rows
# (mixed-bed again, but based on data subset by `prot_n_psm >= 10`)
standPep(
  method_align = MC, 
  col_select = Select_sub, 
  slice_peps_by = exprs(prot_n_psm >= 10),
)

# A side effect in comparing `MC` and `MGKernel`
# (density curves are from the preceding `MGKernel` in (3.2))
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mix_selcols_selrows.png)


# ===================================
# (4) Modified `n_comp`
# ===================================
# start over
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (4.1) first `MGKernel` at `n_comp = 3` for all samples
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG3.png)

# (4.2) a fresh start since changing `n_comp`
# (e.g. `col_select = Select_sub` ignored; instead apply `MGKernel` to all samples)
standPep(
  method_align = MGKernel, 
  n_comp = 2, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05,
  
  # ignored
  col_select = Select_sub, 
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG2.png)


# ===================================
# (5) housekeeping normalizers: 
#     suggest `method_align = MC`
# ===================================
# start over
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)
PSM2Pep()
mergePep()

# initial `Peptide.txt` from `mergePep()` is aligned by MC
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mc.png)

# (GAPDH more abundant in W16 according to MC alignment)
pepHist(scale_log2r = TRUE, col_select = BI_1, filter_ = exprs(gene == "GAPDH"), filename = mcGAPDH.png)

# heat map for `BI_1` samples
# (outputs under `Peptide/Heatmap` folder; for help, ?pepHM)
pepHM(
  col_select = BI_1, 
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
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

# heat map for all samples
# (no GAPDH under `JHU_TMT1` and `PNNL_TMT1`; will be problematic in data alignment next)
pepHM(
  # col_select = BI_1, 
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  
  filter_by = exprs(gene %in% c("GAPDH")),
  filename = "mcGAPDH_all.png",
)

# first renormalize against GAPDH for all samples
# (not to use `MGkernel` as may not have enough data entries from housekeeper(s))
standPep(
  method_align = MC, 
  slice_hskp = exprs(gene %in% c("GAPDH")),
)

# (now log2FC profiles aligned by GAPDH for `BI_1` samples)
pepHist(scale_log2r = TRUE, col_select = BI_1, filename = wrong_hskp.png)

# (no histograms for `JHU_TMT1` and `PNNL_TMT1` samples as no GAPDH underneath)
pepHist(scale_log2r = TRUE, filename = wrong_hskp_all.png)

# not to keep the above example with no data under `JHU_TMT1` and `PNNL_TMT1`
# (surely need different normalizer(s))
unlink(file.path(dat_dir, "Peptide"), recursive = TRUE, force = TRUE)


## Not run: 
PSM2Pep()
mergePep()

# change to `MGKernel`
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 400, 
  maxit = 200, 
  epsilon = 1e-05, 
)

pepHist(scale_log2r = TRUE, col_select = BI_1, filename = mG.png)

# then renormalize against data from GAPDH
# (error: too few entries for fitting with multiple Gaussians)
standPep(
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

