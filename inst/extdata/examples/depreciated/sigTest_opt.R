# ===================================
# Optional significance tests
# for data filtration by pVals
# ===================================

## !!!require the reduced working example in `?load_expts`

pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
  W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

