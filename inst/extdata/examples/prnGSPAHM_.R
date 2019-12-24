# ===================================
# Heat maps of GSPA
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

## prerequisites in significance and enrichment tests
# (see also ?prnSig, ?prnGSPA)
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

prnGSPA(
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c("go_sets", "kegg_sets"),
  impute_na = FALSE,
)

# ===================================
# Distance heat maps of gene sets
# (also interactive networks)
# ===================================
# a `term` is a subset of an `ess_term` if the distance is zero
prnGSPAHM(
  filter_by = exprs(distance <= .6),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  annot_rows = "ess_size",
  filename = show_some_redundancy.png,
)

# human terms only
prnGSPAHM(
  filter_num = exprs(distance <= .95),
  filter_sp = exprs(start_with_str("hs", term)),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  filename = show_more_connectivity.png,
)

# custom color palette
prnGSPAHM(
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter_by = exprs(distance <= .95),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = "custom_colors.png",
)


