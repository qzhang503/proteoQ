# ===================================
# GSPA tests
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


## gene-set enrichment tests
prnGSPA(
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c("go_sets", "kegg_sets"),
  impute_na = FALSE,
)

# only include proteins with two or more identifying peptides
prnGSPA(
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c("go_sets", "kegg_sets"),
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
  impute_na = FALSE,
)

# customized thresholds
# (formula names defined a priori in `pepSig(...)`)
prnGSPA(
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  pval_cutoff = c(5E-2, 5E-2, 1E-10), 
  logFC_cutoff = c(log2(1.05), log2(1.1), log2(1.2)), 
  gspval_cutoff = c(5E-2, 5E-2, 1E-5), 
  max_size = c(Inf, Inf, 120), 
  gset_nms = c("go_sets", "kegg_sets"), 
  filter_by_npep = exprs(prot_n_pep >= 2), 
  impute_na = FALSE, 
)

