# ===================================
# Volcano plots
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

## for all peptides or proteins
# all peptides
pepVol()

# all proteins
prnVol(
  show_labels = TRUE,
  xco = 1.2,
  yco = 0.01,
)

# kinases and prot_n_pep >= 2
prnVol(
  show_labels = TRUE,
  xco = 1.2,
  yco = 0.01,
  filter_prots_by_kin = exprs(kin_attr, prot_n_pep >= 2),
  filename = "kin_npep2.png"
)


## protein subgroups by gene sets
# GSPA
prnGSPA(
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c("go_sets", "kegg_sets"),
  impute_na = FALSE,
)

# filtered by proteins with two or more identifying peptides for visualization
gspaMap(
  pval_cutoff = 5E-3,
  logFC_cutoff = log2(1.2),
  gset_nms = c("go_sets"),
  show_sig = pVal,
  show_labels = TRUE,
  yco = 0.01,
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
  # `filename`(s) will be automated, i.e., by gene-set names
)

# customized thresholds for the corresponding formulae in `pepSig` or `prnSig()`
gspaMap(
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  pval_cutoff = c(5E-2, 5E-2, 1E-10),
  logFC_cutoff = log2(1.2),
  
  show_sig = pVal,
  show_labels = TRUE,
  yco = 0.05,
  filter_by_npep = exprs(prot_n_pep >= 2),
)

