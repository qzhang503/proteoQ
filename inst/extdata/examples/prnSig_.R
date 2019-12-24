# ===================================
# Significance tests
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

## peptides (`Term` etc. are user-defined column keys in expt_smry.xlsx)
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                  "W2.JHU.TMT2-W2.JHU.TMT1", 
                  "W2.PNNL.TMT2-W2.PNNL.TMT1"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

pepVol()

## proteins (formulae matched to `pepSig` by default)
prnSig(impute_na = FALSE)
prnVol()

# note the incongruity in peptide and protein fold changes
# (no measures for peptides but for proteins)

#  sequence |   ref    | sample_1 | sample_2 |  log2FC   
# -------------------------------------------------------
# prnX_pep1 |    0     |   1.15   |    NA    |    NA
# prnX_pep2 |    0     |    NA    |   0.05   |    NA


#  protein  |   ref    | sample_1 | sample_2 |  log2FC   
# -------------------------------------------------------
#   prnX    |    0     |   1.15   |    0.05  |   1.10

## averaged batch effect
# (suggest run both `pepSig` and `prnSig` for consistency)
pepSig(
  impute_na = FALSE,
  W2_loc_2 = ~ Term["(W2.BI.TMT2+W2.BI.TMT1)/2 - (W2.JHU.TMT2+W2.JHU.TMT1)/2"], 
)

prnSig(impute_na = FALSE)

pepVol()
prnVol()

## random effects
# NA imputation (suggested for models with random effects)
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# single
pepSig(
  impute_na = TRUE,
  W2_vs_W16_fix = ~ Term_3["W16-W2"],
  W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set),
)

prnSig(impute_na = TRUE)

pepVol()
prnVol()

# one to multiple (method `lm` for multiple random)
pepSig(
  impute_na = TRUE,
  method = lm,
  W2_vs_W16_fix = ~ Term_3["W16-W2"],
  W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set),
  W2_vs_W16_mix_2 = ~ Term_3["W16-W2"] + (1|TMT_Set) + (1|Color),
)

prnSig(
  impute_na = TRUE,
  method = lm,
)

pepVol()
prnVol()

