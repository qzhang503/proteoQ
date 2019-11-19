# ===================================
# Significance test
# ===================================
scale_log2r <- TRUE

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  # Mascot
  # peptide significance tests
  pepSig(
    impute_na = FALSE, 
    W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
    W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
    W16_vs_W2 = ~ Term_3["W16-W2"], 
  )
  
  # protein significance tests
  prnSig(
    impute_na = FALSE, 
    W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
    W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
    W16_vs_W2 = ~ Term_3["W16-W2"], # between two WHIMs
  )
  
  	
  # MaxQuant
  # peptide significance tests
  pepSig(
    impute_na = FALSE, 
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
    W16_vs_W2_course = ~ Term_2["W16-W2"], 
  )
  
  # protein significance tests
  prnSig(
    impute_na = FALSE, 
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
    W16_vs_W2_course = ~ Term_2["W16-W2"], 
  )
}
## END of RUN `Mascot or Maxquant but not both`

# gene-set enrichment tests
prnGSPA(
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c("go_sets", "kegg_sets"),
  impute_na = FALSE,
)



