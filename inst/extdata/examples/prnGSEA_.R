# ===================================
# GSEA
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

## prerequisites in significance tests
# (see also ?prnSig)
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"], # location effects
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

# all human proteins
prnGSEA(
  var_cutoff = 0, 
  pval_cutoff = 1, 
  logFC_cutoff = log2(1), 
  filter_by_sp = exprs(species == "human"), 
)

# prefiltration by variances, pVals and logFCs
prnGSEA(
  var_cutoff = 0.5,     
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  filter_by_sp = exprs(species == "human"), 
  filename = hu_prefil.txt,
)

# cases that are complete with no missing values
prnGSEA(
  var_cutoff = 0.5,     
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  complete_cases = TRUE, 
  filter_by_sp = exprs(species == "human"), 
  filename = cc.txt,
)

# customized thresholds for the corresponding formulae in `pepSig` or `prnSig()`
prnGSEA(
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  var_cutoff = c(0, 0.2, 0.5), 
  pval_cutoff = c(5E-2, 5E-2, 1E-5),
  logFC_cutoff = c(log2(1.1), log2(1.1), log2(1.2)),
  filter_by_sp = exprs(species == "human"), 
  filename = custom_fil.txt,
)

