# ===================================
# GSPA tests
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

## base
# prerequisites in significance tests (impute_na = FALSE)
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

# GSPA analysis
prnGSPA(
  impute_na = FALSE,
  pval_cutoff = 5E-2, # protein pVal threshold
  logFC_cutoff = log2(1.2), # protein log2FC threshold
  gspval_cutoff = 5E-2, # gene-set pVal threshold
  gslogFC_cutoff = log2(1.2), # gene-set log2FC threshold
  gset_nms = c("go_sets", "kegg_sets"), 
)

# GSPA visualization under volcano plots
gspaMap(
  impute_na = FALSE,
  show_labels = FALSE, 
  gspval_cutoff = 5E-2, 
  gslogFC_cutoff = log2(1.2), 
  show_sig = pVal, 
  xco = 1.2, # position of x-axis cut-off lines
  yco = 0.05, # position of a y-axis cut-off line 
)

## row filtration
prnGSPA(
  impute_na = FALSE,
  gset_nms = c("go_sets", "kegg_sets"),
  filter_prots_by_npep = exprs(prot_n_pep >= 3),
  filename = unifil.txt,
)

# by default, volcano-plot visualization for all GSPA result files, 
#   which are "Protein_GSPA_Z.txt" and "unifil_Protein_GSPA_Z.txt" 
#   up to this point
gspaMap(
  impute_na = FALSE,
  show_labels = FALSE, 
  show_sig = pVal, 
  xco = 1.2, 
  yco = 0.05, 
)

# or manually supply specific secondary file(s) with `df2`; 
# may consider the same primary `fliter_` varargs as those
#   in the corresponding `prnGSPA(...)`
gspaMap(
  df2 = "unifil_Protein_GSPA_Z.txt", 
  filter_prots_by_npep = exprs(prot_n_pep >= 3),
  impute_na = FALSE,
  show_labels = FALSE, 
  show_sig = pVal, 
  xco = 1.2, 
  yco = 0.05, 
)

# possible to visualize all proteins without primary `filter_` varargs, 
#   but it is users' responsibility to note that 
#   "unifil_Protein_GSPA_Z.txt" is from 
#   prnGSA(filter_prots_by_npep = exprs(prot_n_pep >= 3), ...)
gspaMap(
  df2 = "unifil_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  show_labels = FALSE, 
  show_sig = pVal, 
  xco = 1.2, 
  yco = 0.05, 
)

# may apply both primary `filter_` varargs against `df` 
#   and secondary `fliter2_` varargs aganist `df2`
gspaMap(
  df2 = "unifil_Protein_GSPA_Z.txt", 
  filter_prots_by_npep = exprs(prot_n_pep >= 3),
  filter2_ess_size = exprs(ess_size >= 1),   
  impute_na = FALSE,
  show_labels = FALSE, 
  show_sig = pVal, 
  xco = 1.2, 
  yco = 0.05, 
)

# can also visualize results for specific formula(e) and/or 
#   specific `df2`
gspaMap(
  fml_nms = "W16_vs_W2",
  df2 = "unifil_Protein_GSPA_Z.txt", 
  filter_prots_by_npep = exprs(prot_n_pep >= 3),
  filter2_ess_size = exprs(ess_size >= 1),   
  impute_na = FALSE,
  show_labels = FALSE, 
  show_sig = pVal, 
  xco = 1.2, 
  yco = 0.05, 
)

## customized thresholds
# (formula names defined a priori in `pepSig(...)`)
prnGSPA(
  impute_na = FALSE, 
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  pval_cutoff = c(5E-2, 5E-2, 1E-6), 
  logFC_cutoff = c(log2(1.1), log2(1.1), log2(1.2)), 
  gspval_cutoff = c(5E-2, 5E-2, 1E-4), 
  max_size = c(Inf, Inf, 120), 
  gset_nms = c("go_sets"), 
  filter_by_npep = exprs(prot_n_pep >= 3), 
  filename = diffil.txt,
)

gspaMap(
  df2 = "diffil_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  show_labels = FALSE, 
  show_sig = pVal, 
  xco = 1.2, 
  yco = 0.05, 
)

## NA imputation
# if not yet, run prerequisitive NA imputation
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# if not yet, run prerequisitive significance tests at `impute_na = TRUE`
pepSig(
  impute_na = TRUE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = TRUE)

prnGSPA(
  impute_na = TRUE,
  gset_nms = c("go_sets", "kegg_sets"),
  filter_prots_by_npep = exprs(prot_n_pep >= 3),
)

# ===================================
# Distance heat maps of gene sets
# (also interactive networks)
# ===================================
# a `term` is a subset of an `ess_term` if the distance is zero
prnGSPAHM(
  filter2_by = exprs(distance <= .6),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  annot_rows = "ess_size",
  filename = show_some_redundancy.png,
)

# human terms only
prnGSPAHM(
  filter2_num = exprs(distance <= .95),
  filter2_sp = exprs(start_with_str("hs", term)),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  filename = show_more_connectivity.png,
)

# custom color palette
prnGSPAHM(
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter2_by = exprs(distance <= .95),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = custom_colors.png,
)

# can also visualize results for specific `df2` and formula(e)
prnGSPAHM(
  df2 = "diffil_Protein_GSPA_Z_essmap.txt", 
  fml_nms = "W16_vs_W2", 
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter2_by = exprs(distance <= .95),
  filename = w16_vs_w2.png,
)
